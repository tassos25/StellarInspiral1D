! ***********************************************************************
!
!   This file is part of a mesa extension.
!   Authors of this file: Tassos Fragos, Jeff J. Andrews, Matthias U. Kruckow
!
! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module CE_orbit

      use star_def
      use const_def
      use math_lib

!      use binary_evolve, only: eval_rlobe

      implicit none

      contains

! ***********************************************************************
      subroutine CE_orbit_adjust(id, ierr)
         use const_def, only: standard_cgrav, Msun, Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_mass
         real(dp) :: J_tmp, J_init, J_final
         real(dp) :: E_init, E_loss, E_final, E_tmp
         real(dp) :: M_inner, R_inner, M_outer, R_outer, M_final, R_final
         real(dp) :: M_slope, R_slope, M_int, R_int
         real(dp) :: top, bottom, k_final
         real(dp) :: orbital_ang_mom_lost
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp) :: v_rel, v_rel_div_csound, M_encl, rho_at_companion, scale_height_at_companion

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_energy_rate = s% xtra(1)
         CE_companion_position = s% xtra(2)
         CE_companion_mass = s% xtra(4)

         ! If the star is in the initial relaxation phase, skip orbit calculations
         if (s% doing_relax) return
         ! If companion is outside star, skip orbit calculations
         if (CE_companion_position*Rsun > s% r(1)) return

         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)
         M_encl = s% xtra(15)
         v_rel = s% xtra(16)
         v_rel_div_csound = s% xtra(17)
         rho_at_companion = s% xtra(18)
         scale_height_at_companion = s% xtra(19)

         ! Calculate the angular momentum
         J_tmp = (CE_companion_mass * Msun)**2 * M_encl**2 / (CE_companion_mass * Msun + M_encl)
         J_init = sqrt(standard_cgrav * J_tmp * CE_companion_position * Rsun)
         ! Calculate the energies
         E_init = -standard_cgrav * CE_companion_mass * Msun * M_encl / (2.0 * CE_companion_position * Rsun)
         E_loss = CE_energy_rate * s% dt
         E_final = E_init - E_loss

         ! Move from outside of star in to find cell containing companion
         E_tmp = 0d0
         k = 1
         do while (E_tmp > E_final)
            k = k + 1
            M_inner = s% m(k)
            R_inner = s% r(k)
            E_tmp = -standard_cgrav * CE_companion_mass * Msun * M_inner / (2.0 * R_inner)
         end do

         ! If companion is outside star, set k to 2
         if (k < 2) k = 2

         ! save end points of cell containing companion
         M_inner = s% m(k)
         R_inner = s% r(k)
         M_outer = s% m(k-1)
         R_outer = s% r(k-1)

         ! We could choose to interpolate for R using M as the independent variable. Instead, here we
         ! linearly interpolate across cell (using k as the independent variable)
         ! M_slope = (m(k) - m(k-1))/((k) - (k-1)) = m(k) - m(k-1) = M_inner - M_outer
         ! R_slope = (r(k) - r(k-1))/((k) - (k-1)) = r(k) - r(k-1) = R_inner - R_outer
         M_slope = M_inner - M_outer
         R_slope = R_inner - R_outer
         ! The final values will be:
         ! M_final = m(k) - (k - k_final)/((k) - (k-1)) * (m(k) - m(k-1)) = M_inner - (k - k_final) * M_slope
         ! R_final = r(k) - (k - k_final)/((k) - (k-1)) * (r(k) - r(k-1)) = R_inner - (k - k_final) * R_slope

         ! Given the final energy, E_final, determine the resulting k that solves the equation
         ! 0 = 2*E_final*R_final + G*M_2*M_final = 2*E_final*(R_inner - (k - k_final)*R_slope) + G*M_2*(M_inner - (k - k_final)*M_slope)
         ! k_final = -(2*E_final*(R_inner - k*R_slope) + G*M_2*(M_inner - k*M_slope))/(2*E_final*R_slope + G*M_2*M_slope)
         !         = -(2*E_final*R_int + G*M_2*M_int)/(2*E_final*R_slope + G*M_2*M_slope)
         ! where M_int = M_inner - k*M_slope and R_int = R_inner - k*R_slope
         M_int = M_inner - real(k)*M_slope
         R_int = R_inner - real(k)*R_slope
         top = 2.0 * E_final * R_int + standard_cgrav * CE_companion_mass * Msun * M_int
         bottom = 2.0 * E_final * R_slope + standard_cgrav * CE_companion_mass * Msun * M_slope
         k_final = -top / bottom

         ! Now use the interpolations and the derived k_final, determine the resulting enclosed mass and separation
         ! M_final = M_inner - (k - k_final)*M_slope = M_int + k_final*M_slope
         ! R_final = R_inner - (k - k_final)*R_slope = R_int + k_final*R_slope
         M_final = M_int + M_slope * k_final
         R_final = R_int + R_slope * k_final

         ! save enclosed mass and separation in xtra(9) and xtra(2), respectively
         ! so that we output it in the history data
         s% xtra(2) = R_final/Rsun
         s% xtra(9) = M_final/Msun

         ! Calculate the angular momentum lost to the star's envelope
         J_tmp = (CE_companion_mass * Msun)**2 * M_final**2 / (CE_companion_mass * Msun + M_final)
         J_final = sqrt(standard_cgrav * J_tmp * R_final)

         ! The angular momentum that is lost from the orbit of the companion
         ! is added to the envelope of the donor.
         orbital_ang_mom_lost = J_final - J_init
         !We save in s% xtra(6) the total torque that will be applied to the Envelope
         s% xtra(6) = max(0.0d0, -orbital_ang_mom_lost/s% dt) 

         ! Keep track of orbital energy and angular momentum
         s% xtra(8) = E_final
         s% xtra(10) = J_final

         ! For diagnostics

         write(*,*) "Final k:                  ", k_final, &
                     " Step Time:             ", s% dt
         write(*,*) "Previous Enclosed Mass:   ", M_encl/Msun, &
                     " Final Enclosed Mass:   ", M_final/Msun
         write(*,*) "Previous Separation:      ", CE_companion_position, &
                     " Final Separation:      ", R_final/Rsun, &
                     " Stellar Radius:", s% r(1)/Rsun
         write(*,*) "Inner Accretion Radius:   ", R_acc_low/Rsun, &
                     " Outer Accretion Radius:", R_acc_high/Rsun
         write(*,*) "Relative Velocity:        ", v_rel, &
                     " Mach Number:           ", v_rel_div_csound
         write(*,*) "Previous Orbital Energy:  ", E_init, &
                     " Final Orbital Energy:  ", E_final, &
                     " Dissipated Energy Rate:", s% xtra(1)
         write(*,*) "Total Stellar Energy:     ", s% total_energy, "                   ", &
                     " Dissipated Rotational Energy:", s% xtra(6) * s% omega(k), &
                     s% xtra(6)*2.*3.14/AtoP(M_encl, CE_companion_mass*Msun, CE_companion_position*Rsun)
         write(*,*) "Previous Angular momentum:", J_init, &
                     " Final Angular momentum:", J_final, &
                     " Dissipated Angular Momentum Rate:", s% xtra(6)

         ! After adjusting the orbit, let's call the check_merger routine
         call check_merger(id, ierr)
         if (s% lxtra(1)) write(*,*) "MERGER!!!"
      end subroutine CE_orbit_adjust


! ***********************************************************************
      subroutine check_merger(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_companion_position, CE_companion_radius, CE_companion_mass, rl_core, rl_companion

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_companion_position = s% xtra(2)
         CE_companion_radius = s% xtra(3)
         CE_companion_mass = s% xtra(4)

         rl_core = eval_rlobe(s% he_core_mass,CE_companion_mass,CE_companion_position)*CE_companion_position
         rl_companion = eval_rlobe(CE_companion_mass,s% he_core_mass,CE_companion_position)*CE_companion_position

         ! Merger condition
         ! Merge if either the companion or the core fills its rochelobe
         if ((rl_companion .le. CE_companion_radius) .or. (rl_core .le. s% he_core_radius)) then
            s% lxtra(1) = .true.
            s% xtra(2) = 0.0
         end if
      end subroutine check_merger


! ***********************************************************************
      subroutine calc_quantities_at_comp_position(id, ierr)

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, k
         real(dp) :: CE_companion_position, CE_companion_mass, omega_at_companion, csound_at_companion, P, M2
         real(dp) :: v_rel, v_rel_div_csound, v_rel_at_companion
         real(dp) :: R_rel, R_acc, R_acc_low, R_acc_high
         real(dp) :: M_encl, rho_at_companion, scale_height_at_companion
         real(dp) :: C_1, C_2, M_env, t_kh_env, log_eta
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_companion_position = s% xtra(2)
         CE_companion_mass = s% xtra(4)
         ! Calculate quantities at the position of the companion

         ! If companion is outside star the set default values and return
         if (CE_companion_position*Rsun > s% r(1)) then
            !saving these values to xtra variable so that they are used in different CE_inject cases,
            ! in the torque calculations, and saved in the history file
            s% xtra(12) = 0.0d0 !R_acc
            s% xtra(13) = 0.0d0 !R_acc_low
            s% xtra(14) = 0.0d0 !R_acc_high
            s% xtra(15) = s% m(1) !M_encl
            s% xtra(16) = 0.0d0 !v_rel
            s% xtra(17) = 0.0d0 !v_rel_div_csound
            s% xtra(18) = 0.0d0 !rho_at_companion
            s% xtra(19) = 0.0d0 !scale_height_at_companion
            return
         endif

         ! Calculate the wind enhancement term due to giant star pulsations
         ! Formulation from Yoon & Cantiello (2010)
         C_1 = 9.219e-6
         C_2 = 0.0392844
         M_env = (s% mstar - s% he_core_mass * Msun)
         t_kh_env = standard_cgrav * s% mstar * M_env / (s% L(1) * s% r(1)) / secyer
         log_eta = C_1 * (s% r(1)/Rsun)**2 / (s% mstar / Msun) * t_kh_env**(-0.315) - C_2
         s% xtra(21) = 10.0**log_eta

         ! Keplerian velocity calculation depends on mass contained within a radius
         ! Include all the enclosed cells
         ! Add to it the enclosed mass of the current cell
         k = 1
         do while (s% r(k) > CE_companion_position * Rsun)
            k = k + 1
         end do

         M_encl = s% m(k)
         M_encl = M_encl + s% dm(k-1) * (CE_companion_position*Rsun - s% r(k)) / (s% r(k-1) - s% r(k))

         M2 = CE_companion_mass * Msun

         ! Determine orbital period in seconds
         P = AtoP(M_encl, M2, CE_companion_position*Rsun)

         ! Determine Keplerian velocity. Then subtract the local rotation velocity
         if (s%rotation_flag) then
            omega_at_companion = s% omega(k) + (CE_companion_position*Rsun - s% r(k)) * &
                 (s% omega(k-1)-s% omega(k)) / (s% r(k-1) - s% r(k))
            v_rel_at_companion = 2.0 * pi * CE_companion_position*Rsun / P - &
                 omega_at_companion * CE_companion_position*Rsun ! local rotation velocity = omega * r
         else
            v_rel_at_companion = 2.0 * pi * CE_companion_position*Rsun / P
         endif

         rho_at_companion = s% rho(k) + (CE_companion_position*Rsun - s% r(k)) * &
              (s% rho(k-1)-s% rho(k)) / (s% r(k-1) - s% r(k))

         ! Determine Mach number
         csound_at_companion = s% csound(k) + (CE_companion_position*Rsun - s% r(k)) * &
              (s% csound(k-1)-s% csound(k)) / (s% r(k-1) - s% r(k))
         v_rel_div_csound = v_rel_at_companion / csound_at_companion

         ! Determine scale height of envelope at companion's radius
         scale_height_at_companion =  s% scale_height(k) + (CE_companion_position*Rsun - s% r(k)) * &
              (s% scale_height(k-1)-s% scale_height(k)) / (s% r(k-1) - s% r(k))

         ! Determine accretion radius

         ! Accretion radius for a constant density medium
         R_acc = 2.0 * standard_cgrav * M2 / &
             ((v_rel_at_companion*v_rel_at_companion) + csound_at_companion*csound_at_companion)

         ! To be done appropriately, inner R_acc needs to be calculated separately from the outer R_acc
         ! Find lower R_acc, starting at the star's position
         j = k
         do while (j < s% nz)

           R_rel = abs(s% r(j) - CE_companion_position*Rsun)
           v_rel = 2.0 * pi * CE_companion_position*Rsun / P - s% omega(j) * s% r(j)
           if ((R_rel) > 2.0 * standard_cgrav * M2 / ((v_rel*v_rel) + s% csound(j) * s% csound(j))) exit
           j = j + 1
         end do
         R_acc_low = R_rel

         j = k
         do while (j > 1)
           R_rel = abs(s% r(j) - CE_companion_position*Rsun)
           v_rel = 2.0 * pi * CE_companion_position*Rsun / P - s% omega(j) * s% r(j)
           if ((R_rel) > 2.0 * standard_cgrav * M2 / ((v_rel*v_rel) + s% csound(j) * s% csound(j))) exit
           j = j - 1
         end do
         R_acc_high = R_rel

         !saving these values to xtra variable so that tey are used in different CE_inject cases,
         ! in the torque calculations, and saved in the history file
         s% xtra(12) = R_acc
         s% xtra(13) = R_acc_low
         s% xtra(14) = R_acc_high
         s% xtra(15) = M_encl
         s% xtra(16) = v_rel_at_companion
         s% xtra(17) = v_rel_div_csound
         s% xtra(18) = rho_at_companion
         s% xtra(19) = scale_height_at_companion
      end subroutine calc_quantities_at_comp_position


! ***********************************************************************
      real (dp) function AtoP(M1, M2, A)
         real(dp), intent(in) :: M1 ! in g
         real(dp), intent(in) :: M2 ! in g
         real(dp), intent(in) :: A ! in cm

         ! Kepler's 3rd Law - return orbital period in seconds
         AtoP = 2.0*pi * sqrt(A*A*A / (standard_cgrav * (M1+M2)))
      end function AtoP


! ***********************************************************************
      real(dp) function TukeyWindow(x,a)
         use const_def, only: dp, pi
         real(dp), intent(in) :: x, a

         if ((x .le. -0.5) .or. (x .ge. 0.5)) then
            TukeyWindow = 0.
         else if (x .le. -0.5 + a) then
            TukeyWindow = 0.5 - 0.5*cos(pi*(x+0.5)/a)
         else if (x .ge. 0.5 - a) then
            TukeyWindow = 0.5 - 0.5*cos(-pi*(x-0.5)/a)
         else
            TukeyWindow = 1.
         endif
      end function TukeyWindow


! ***********************************************************************
      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow(m1/m2,one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p(q))
      end function eval_rlobe


      end module CE_orbit
