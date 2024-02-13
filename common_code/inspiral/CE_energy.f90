! ***********************************************************************
!
!   This file is part of a mesa extension.
!   Authors of this file: Tassos Fragos, Jeff J. Andrews, Matthias U. Kruckow, Jaime Roman-Garza
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

      module CE_energy

      ! NOTE: if you'd like to have some inlist controls for your routine,
      ! you can use the x_ctrl array of real(dp) variables that is in &controls
      ! e.g., in the &controls inlist, you can set
      !     x_ctrl(1) = <my_special_param>
      ! where <my_special_param> is a real value such as 0d0 or 3.59d0
      ! Then in your routine, you can access that by
      !     s% x_ctrl(1)
      ! of course before you can use s, you need to get it using the id argument.
      ! here's an example of how to do that -- add these lines at the start of your routine:
      !         use star_lib, only: star_ptr
      !         type (star_info), pointer :: s
      !         call star_ptr(id, s, ierr)
      !         if (ierr /= 0) then ! OOPS
      !            return
      !         end if
      !
      ! for integer control values, you can use x_integer_ctrl
      ! for logical control values, you can use x_logical_ctrl

      use star_def
      use const_def
      use CE_orbit, only: AtoP, TukeyWindow, calc_quantities_at_comp_position

      implicit none

      contains

! ***********************************************************************
      subroutine CE_inject_energy(id, ierr)
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: CE_test_case, k
         real(dp) :: CE_companion_position

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_test_case = s% x_integer_ctrl(1)
         CE_companion_position = s% xtra(2)

         ! If the star is in the initial relaxation phase, skip energy calculations
         if (s% doing_relax) return
         ! If companion is outside star, skip energy calculations
         if (CE_companion_position*Rsun > s% r(1)) return

         ! If system merged, skip energy deposition
         if (s% lxtra(1)) then
            s% xtra(1) = 0.0d0
            s% xtra(20) = 0.0d0
            s% xtra(22) = 0.0d0
            return
         endif

         ! Call functions to calculate test cases
         if (CE_test_case == 1) then   ! Heat whole hydrogen envelope

            call CE_inject_case1(id, ierr)

         else if (CE_test_case == 2) then   ! Heat just outside helium core

            call CE_inject_case2(id, ierr)

         else if (CE_test_case == 3) then   ! (Outdated) Energy added based on Ostriker (1999)

            call CE_inject_case3(id, ierr)

         else if (CE_test_case == 4) then   ! Energy based on MacLeod & Ramirez-Ruiz

            call CE_inject_case4(id, ierr)

         else if (CE_test_case == 5) then   ! Energy based on Lee & Stahler (2011)

            call CE_inject_case5(id, ierr)
         

         endif
      end subroutine CE_inject_energy


! ***********************************************************************
      subroutine CE_inject_case1(id, ierr)

         use const_def, only: Msun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: ff, mass_to_be_heated, m_bot
         real(dp) :: CE_energy_rate

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_energy_rate = s% x_ctrl(1)

         ! mass (g) of the bottom of the (outer) convective envelope
           ! Based on the inner edge of the convective envelope
!          m_bot = s% conv_mx1_bot * s% mstar
           ! Based on the helium core
         m_bot = s% he_core_mass * Msun

         !First calculate the mass in which the energy will be deposited
         mass_to_be_heated = 0.0
         do k = 1, s% nz
            ff = EnvelopeWindow(s% m(k), m_bot)
            mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
         end do

         !Now redo the loop and add the extra specific heat
         do k = 1, s% nz
            s% extra_heat(k)%val = CE_energy_rate / mass_to_be_heated * EnvelopeWindow(s% m(k), m_bot)
         end do

         ! Save the total erg/second added in this time step
         s% xtra(1) = CE_energy_rate

         contains

         real(dp) function EnvelopeWindow(m_interior, m_bot)
            use const_def, only: pi, Msun
            real(dp), intent(in) :: m_interior, m_bot

            ! arctan function causes a smooth transition from 0 to unity around m_bot.
            ! The 0.002 in the denominator sets the width of the transition, in this case
            ! calibrated so the transition region is roughly 0.01 Msun.
            EnvelopeWindow = 1./pi * atan((s% m(k) - m_bot) / (0.002*Msun)) + 0.5

         end function EnvelopeWindow
      end subroutine CE_inject_case1


! ***********************************************************************
      subroutine CE_inject_case2(id, ierr)

         use const_def, only: Msun, Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_energy_rate

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Get input controls

         CE_energy_rate = s% x_ctrl(1)

         call CE_set_extra_heat(id, CE_energy_rate, ierr)

         ! Save the total erg/second added in this time step
         s% xtra(1) = CE_energy_rate

      end subroutine CE_inject_case2


! ***********************************************************************
      subroutine CE_inject_case3(id, ierr)

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: time, M2
         real(dp) :: I, F_drag, F_coef
         real(dp) :: a_tukey, mass_to_be_heated, ff
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp) :: v_rel, v_rel_div_csound, M_encl, rho_at_companion, scale_height_at_companion

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Alternative energy source here

         ! Get input controls
         CE_energy_rate = s% xtra(1)
         CE_companion_position = s% xtra(2)
         CE_companion_radius = s% xtra(3)
         CE_companion_mass = s% xtra(4)
         CE_n_acc_radii = s% xtra(5)

         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)
         M_encl = s% xtra(15)
         v_rel = s% xtra(16)
         v_rel_div_csound = s% xtra(17)
         rho_at_companion = s% xtra(18)
         scale_height_at_companion = s% xtra(19)

!         ! This is incorrect, but for now, not completely crazy
!         R_acc = (R_acc_low + R_acc_high) / 2.0

         ! Determine drag force
         ! Equations from Ostriker (1999) ApJ, 513, 252
         time = 1.0 ! FIX THIS: What is time here?
         if (v_rel_div_csound .lt. 1.0) then
            I = 0.5 * log((1.0+v_rel_div_csound)/(1.0-v_rel_div_csound)) - v_rel_div_csound
         else
            I = 0.5 * log((v_rel_div_csound+1.0)/(v_rel_div_csound-1.0)) + log(v_rel*time / R_acc)
         end if

         M2 = CE_companion_mass * Msun

         F_coef = 4.0 * pi * standard_cgrav * standard_cgrav * M2 * M2 * rho_at_companion / (v_rel*v_rel)
         F_drag = -F_coef * I

         ! Add hydrodynamic drag
         F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2.0 * rho_at_companion * v_rel**2

         ! Total energy rate= drag force * velocity
         CE_energy_rate = F_drag * v_rel

         call CE_set_extra_heat(id, CE_energy_rate, ierr)

         ! Save the total erg/second added in this time step
         s% xtra(1) = CE_energy_rate
      end subroutine CE_inject_case3


! ***********************************************************************
      subroutine CE_inject_case4(id, ierr)

         use const_def, only: Rsun, Msun, pi, standard_cgrav
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k_bottom
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: M2, R2
         real(dp) :: I, F_drag
         real(dp) :: a_tukey, mass_to_be_heated, ff
         real(dp) :: F_DHL, f1, f2, f3, e_rho
         real(dp) :: mdot_macleod, mdot_HL, log_mdot_factor, a1, a2, a3, a4, L_acc
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp) :: v_rel, v_rel_div_csound, csound
         real(dp) :: M_encl, rho_at_companion, scale_height_at_companion
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Alternative energy source here

         ! Get input controls
         CE_energy_rate = s% xtra(1)
         CE_companion_position = s% xtra(2)
         CE_companion_radius = s% xtra(3)
         CE_companion_mass = s% xtra(4)
         CE_n_acc_radii = s% xtra(5)

         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)
         M_encl = s% xtra(15)
         v_rel = s% xtra(16)
         v_rel_div_csound = s% xtra(17)
         rho_at_companion = s% xtra(18)
         scale_height_at_companion = s% xtra(19)
         csound = v_rel / v_rel_div_csound

!         ! For a first approximation, let's use the average R_acc
!         R_acc = (R_acc_low + R_acc_high) / 2.0
         F_DHL = pi * R_acc**2 * rho_at_companion * v_rel**2

         ! Hydrodynamic drag from MacLeod & Ramirez-Ruiz (2014)
         f1 = 1.91791946d0
         f2 = -1.52814698d0
         f3 = 0.75992092
         e_rho = R_acc / scale_height_at_companion
         F_drag = F_DHL*(f1 + f2*e_rho +f3*e_rho**2)
         F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2 * rho_at_companion * v_rel**2

         ! Add hydrodynamic drag
         F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2.0 * rho_at_companion * v_rel**2

         ! Mass accretion from MacLeod & Ramirez-Ruiz (2014)
         a1 = -2.14034214
         a2 = 1.94694764
         a3 = 1.19007536
         a4 = 1.05762477

         M2 = CE_companion_mass * Msun
         R2 = CE_companion_radius * Rsun    ! NS radius is 10 km

         log_mdot_factor = a1 + a2 / (1.0 + a3*e_rho + a4*e_rho**2)
         mdot_HL = pi * R2**2 * rho_at_companion * v_rel
         mdot_macleod = mdot_HL * 10.0**log_mdot_factor
         s% xtra(22) = mdot_HL
         s% xtra(23) = mdot_macleod

         ! Accretion luminosity luminosity
         L_acc = standard_cgrav * M2 / R2 * mdot_macleod

         ! Total energy rate = drag force * velocity
         if (s% x_logical_ctrl(2)) then
            CE_energy_rate = F_drag * max(v_rel,0.0d0) + L_acc ! Include accretion luminosity depending on inlist input
         else
            CE_energy_rate = F_drag * max(v_rel,0.0d0)
         end if

         call CE_set_extra_heat(id, CE_energy_rate, ierr)

         ! Save the total erg/second added in this time step
         s% xtra(1) = CE_energy_rate
         s% xtra(20) = L_acc

      end subroutine CE_inject_case4


! ***********************************************************************
      subroutine CE_inject_case5(id, ierr)

         use const_def, only: Rsun, Msun, pi, standard_cgrav
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k_bottom
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: M2, R2
         real(dp) :: F_drag
         real(dp) :: F_DHL, f1, f2, f3, e_rho
         real(dp) :: mdot, mdot_macleod, mdot_HL, L_acc, a1, a2, a3, a4
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp) :: v_rel, beta, M_encl, csound
         real(dp) :: rho_at_companion, scale_height_at_companion
         real(dp) :: drag_factor, log_mdot_factor, lambda_squared
         real(dp) :: F_drag_subsonic, F_drag_supersonic
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Alternative energy source here

         ! Get input controls
         CE_energy_rate = s% xtra(1)
         CE_companion_position = s% xtra(2)
         CE_companion_radius = s% xtra(3)
         CE_companion_mass = s% xtra(4)
         CE_n_acc_radii = s% xtra(5)

         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)
         M_encl = s% xtra(15)
         v_rel = s% xtra(16)
         beta = s% xtra(17)
         rho_at_companion = s% xtra(18)
         scale_height_at_companion = s% xtra(19)
         csound = v_rel / beta

         M2 = CE_companion_mass * Msun
         R2 = CE_companion_radius * Rsun    ! NS radius is 10 km

         mdot_HL = pi * R_acc**2 * rho_at_companion * v_rel
         s% xtra(22) = mdot_HL

         if (s% x_integer_ctrl(2) == 2) then

            ! Dynamical drag from MacLeod & Ramirez-Ruiz (2014)
            f1 = 1.91791946d0
            f2 = -1.52814698d0
            f3 = 0.75992092
            e_rho = R_acc / scale_height_at_companion
            drag_factor = (f1 + f2*e_rho +f3*e_rho**2)

            ! Mass accretion from MacLeod & Ramirez-Ruiz (2014)
            a1 = -2.14034214
            a2 = 1.94694764
            a3 = 1.19007536
            a4 = 1.05762477

            log_mdot_factor = a1 + a2 / (1.0 + a3*e_rho + a4*e_rho**2)

            !!Drag at the Supesonic regime based on Macleod & Ramirez-Ruiz (2015)
            ! R_acc = (R_acc_low + R_acc_high) / 2.0
            F_drag_supersonic = pi * R_acc**2 * rho_at_companion * v_rel**2
            F_drag_supersonic = F_drag_supersonic * drag_factor

            !!Drag at the subsonic regime based on Ostriker, E. (1999)
            F_drag_subsonic = (1./2. * log((1+beta)/(1-beta))-beta)
            F_drag_subsonic = F_drag_subsonic * 4.*pi*rho_at_companion*(standard_cgrav*CE_companion_mass*Msun)**2. /v_rel**2.

            if (beta < 0.9 .and. beta > 0.001) then
               F_drag =  F_drag_subsonic
            else if (beta > 0.99) then
               F_drag = F_drag_supersonic
            else if (beta <= 0.001) Then
               F_drag = 0.
            else
               !smooth between the two regimes
               F_drag = ((beta-0.9)*F_drag_supersonic + (0.99-beta)*F_drag_subsonic)/0.09
            endif

            ! Add hydrodynamical drag
            F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2 * rho_at_companion * v_rel**2

            mdot = mdot_HL * 10.0**log_mdot_factor
            s% xtra(23) = mdot

            ! Accretion luminosity luminosity: 10% efficiency
            L_acc = 0.1 * standard_cgrav * M2 / R2 * mdot_HL
         
         else if (s% x_integer_ctrl(2) == 4) then

            if (s% x_logical_ctrl(7)) then 
               ! Accretion luminosity luminosity: 10% efficiency
               L_acc = (s% xtra(30) * Msun * Rsun / (secyer**2.0) )/(s% xtra(16) * Rsun / secyer)
            else
               stop "The value of x_integer_ctrl(2) = 4 is only available for x_logical_ctrl(7) = .true."
            endif
         else if (s% x_integer_ctrl(2) == 3) then

            lambda_squared = exp(3.0) / 16.0

            ! Dimensionless
            mdot = 2.0 * sqrt(lambda_squared + beta*beta) / (1.0 + beta*beta)**2
            ! Add in dimensions
            mdot = mdot * 2.0 * pi * rho_at_companion * standard_cgrav**2 * M2**2 / csound**3
            s% xtra(23) = mdot

            ! Drag force
            F_drag =  beta * csound * mdot
            ! Add hydrodynamic drag
            F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2.0 * rho_at_companion * v_rel**2
            ! Accretion luminosity luminosity: 10% efficiency
            L_acc = 0.1 * standard_cgrav * M2 / R2 * mdot

         else if (s% x_integer_ctrl(2) == 1) then

            F_drag = pi * R_acc**2 * rho_at_companion * v_rel**2

            ! Add hydrodynamic drag
            F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2.0 * rho_at_companion * v_rel**2

            mdot = pi * R_acc**2 * rho_at_companion * v_rel
            s% xtra(23) = mdot

            ! Accretion luminosity luminosity: 10% efficiency
            L_acc = 0.1 * standard_cgrav * M2 / R2 * mdot

         else
            stop "x_integer_ctrl(2) is not defined"
         endif

         ! Limit accretion luminosity to Eddington rate
         L_acc = min(L_acc, 1.26e38 * CE_companion_mass)

         ! Total energy rate = drag force * velocity
         if (.not.(s% x_logical_ctrl(7))) then 
            CE_energy_rate = F_drag * max(v_rel,0.0d0)
         else 
            CE_energy_rate = s% xtra(1)
         end if 

         if (s% x_logical_ctrl(2)) then
            if (s% x_integer_ctrl(3) == 1) then
               call CE_set_extra_heat(id, CE_energy_rate + L_acc, ierr)
            else if (s% x_integer_ctrl(3) == 2) then
               call CE_set_extra_heat2(id, CE_energy_rate + L_acc, ierr)
            else
               stop "s% x_integer_ctrl(3) is not defined"
            endif
         else
            if (s% x_integer_ctrl(3) == 1) then
               call CE_set_extra_heat(id, CE_energy_rate, ierr)
            else if (s% x_integer_ctrl(3) == 2) then
               call CE_set_extra_heat2(id, CE_energy_rate, ierr)
            else
               stop "s% x_integer_ctrl(3) is not defined"
            endif
         end if

         ! Save the total erg/second added in this time step
         s% xtra(1) = CE_energy_rate
         s% xtra(20) = L_acc
      end subroutine CE_inject_case5


! ***********************************************************************
      subroutine CE_set_extra_heat(id, CE_energy_rate, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: CE_energy_rate
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: a_tukey, mass_to_be_heated, ff
         real(dp) :: R_acc, R_acc_low, R_acc_high
         integer :: k, k_bottom
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_companion_position = s% xtra(2)
         CE_n_acc_radii = s% xtra(5)

         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)

         ! Tukey window scale
         a_tukey = 0.5

         ! First calculate the mass in which the energy will be deposited
         mass_to_be_heated = 0.0
         do k = 1, s% nz
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if

            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0 * R_acc), a_tukey)
            mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
            !Energy should be deposited only on the envelope of the star and not in the core
            !When we reach the core boundary we exit the loop
            if (s% m(k) < s% he_core_mass * Msun) exit
         end do
         !this is the limit in k of the boundary between core and envelope
         k_bottom = k-1
         ! If companion is outside star, set mass_to_be_heated arbitrarily low
         if (mass_to_be_heated == 0.) mass_to_be_heated = 1.0

         ! Now redo the loop and add the extra specific heat
         do k = 1, k_bottom
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if

            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0 * R_acc), a_tukey)
            s% extra_heat(k)%val = CE_energy_rate / mass_to_be_heated * ff
            s% xtra6_array(k) = CE_energy_rate / mass_to_be_heated * ff
         end do
      end subroutine CE_set_extra_heat


! ***********************************************************************
      subroutine CE_set_extra_heat2(id, CE_energy_rate, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: CE_energy_rate
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: a_tukey, volume_to_be_heated, ff
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp), DIMENSION(:), ALLOCATABLE :: cell_dr
         integer :: k, k_bottom
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_companion_position = s% xtra(2)
         CE_n_acc_radii = s% xtra(5)

         call calc_quantities_at_comp_position(id, ierr)

         do k=1, s% nz 
            s% xtra6_array(k) = 0.d0
         end do 

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)

         ! Tukey window scale
         a_tukey = 0.5

         ! define cell width
         allocate(cell_dr(s% nz))
         do k = 2, s% nz
            cell_dr(k-1) = s% rmid(k-1) - s% rmid(k)
         end do
         cell_dr(s% nz) = s% rmid(s% nz) - s% R_center

         ! First calculate the mass in which the energy will be deposited
         volume_to_be_heated = 0.0
         do k = 1, s% nz
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if

            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0 * R_acc), a_tukey)
            volume_to_be_heated = volume_to_be_heated + 4.0d0 * pi * s% r(k) * s% r(k) * cell_dr(k) * ff
            !Energy should be deposited only on the envelope of the star and not in the core
            !When we reach the core boundary we exit the loop
            if (s% m(k) < s% he_core_mass * Msun) exit
         end do
         !this is the limit in k of the boundary between core and envelope
         k_bottom = k-1
         ! If companion is outside star, set mass_to_be_heated arbitrarily low
         if (volume_to_be_heated == 0.) volume_to_be_heated = 1.0

         ! Now redo the loop and add the extra specific heat
         do k = 1, k_bottom
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if

            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0 * R_acc), a_tukey)
            s% extra_heat(k)%val = CE_energy_rate * (4.0d0 * pi * s% r(k) * s% r(k) * cell_dr(k) * ff / volume_to_be_heated) / s% dm(k)
            s% xtra6_array(k) = CE_energy_rate * (4.0d0 * pi * s% r(k) * s% r(k) * cell_dr(k) * ff / volume_to_be_heated) / s% dm(k)
         end do
         deallocate(cell_dr)
      end subroutine CE_set_extra_heat2


      end module CE_energy
