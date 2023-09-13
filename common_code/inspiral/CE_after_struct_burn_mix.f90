! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module CE_after_struct_burn_mix

      use const_def, only: dp,ln10, Msun, secyer
      use star_def
      use ionization_def

      implicit none


      contains

      ! set use_other_after_struct_burn_mix = .true. to enable this.
      ! your routine will be called after the standard struct_burn_mix routine


      subroutine CE_other_after_struct_burn_mix(id, dt, res)
         use const_def, only: dp
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: dt
         integer, intent(out) :: res ! keep_going, redo, retry, backup, terminate

         call CE_remove_unbound_envelope(id, dt, res)
         call calc_recombination_after_struct_burn_mix(id, dt, res)

         res = keep_going
      end subroutine CE_other_after_struct_burn_mix



      subroutine CE_remove_unbound_envelope(id, dt, res)

        use const_def, only: Rsun, Msun
        integer, intent(in) :: id
        real(dp), intent(in) :: dt
        integer, intent(out) :: res ! keep_going, redo, retry, backup, terminate
         type (star_info), pointer :: s
         integer :: ierr, k
         real(dp) :: mass_to_remove, CE_mdot, v_rad, vrot, n_tau_to_remove
         real(dp) :: total_envelope_binding_energy

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         k=1
         mass_to_remove = 0.0d0
         n_tau_to_remove = s% x_ctrl(18)
         do while ((k < s% nz) .and. (s% tau(k) < n_tau_to_remove))
            if (.not. is_bound(k)) then
               mass_to_remove = mass_to_remove + s% dm(k)
            else
               exit
            endif
            k=k+1
         enddo

         k=1

         ! Diagnostic to determine envelope binding energy
         ! Includes internal energy
         k=1
         total_envelope_binding_energy = 0.0
         do while ((k < s% nz) .and. (s% m(k) > s% he_core_mass * Msun))

            ! This is necessary depending on the velocity version used
            if (s% u_flag) then
               v_rad = s% u(k)
            else
               v_rad = s% v(k)
            endif

            if (s% rotation_flag) then
               vrot = s% omega(k) * s% r(k)
               total_envelope_binding_energy = total_envelope_binding_energy + &
                              (s% energy(k) - s% cgrav(k)*s% m_grav(k)/s% r(k) + &
                              0.5d0*v_rad*v_rad + 0.5d0*vrot*vrot) * s% dm(k)
            else
               total_envelope_binding_energy = total_envelope_binding_energy + &
                              (s% energy(k) - s% cgrav(k)*s% m_grav(k)/s% r(k) + &
                              0.5d0*v_rad*v_rad) * s% dm(k)
            endif

            k=k+1
         enddo
         s% xtra11 = total_envelope_binding_energy ! In erg


         CE_mdot = - (mass_to_remove) / dt !In gr/s

         if (s% dt .le. s% mass_change_full_off_dt) then
            CE_mdot = 0.0
            write (*,*) "*", s% dt, s% mass_change_full_off_dt
         else if ((s% dt .ge. s% mass_change_full_off_dt) .and. (s% dt .le. s% mass_change_full_on_dt)) then
            CE_mdot = CE_mdot * (s% dt - s% mass_change_full_off_dt)/ &
                (s% mass_change_full_on_dt - s% mass_change_full_off_dt)
            write (*,*)"*", s%dt, s% mass_change_full_off_dt, s% mass_change_full_on_dt
         endif



         s% xtra7 = CE_mdot



         res = keep_going


         contains

            logical function is_bound(k)
               integer, intent(in) :: k
               real(dp) :: val, f_energy, vrot
               logical :: include_internal_energy


               include_internal_energy = s% x_logical_ctrl(1)
               f_energy = logic2dbl(include_internal_energy)

               !When include_internal_energy = true, we shoud include the internal energy
               !only when the shell is mechanically unstable, i.e. Gamma1<4./3.
               if (s% gamma1(k) >= 4./3.) f_energy = 0.d0

               ! This is necessary depending on the velocity version used
               if (s% u_flag) then
                  v_rad = s% u(k)
               else
                  v_rad = s% v(k)
               endif

               !If roation is on, then we add the rotational kinetic energy into the equation
               if (s% rotation_flag) then
                  vrot = s% omega(k) * s% r(k)
                  val = f_energy * s% energy(k) - s% cgrav(k)*s% m_grav(k)/s% r(k) + &
                              0.5d0*v_rad*v_rad + 0.5d0*vrot*vrot
               else
                  val = f_energy * s% energy(k) - s% cgrav(k)*s% m_grav(k)/s% r(k) + &
                              0.5d0*v_rad*v_rad
               endif

               if (val > 0.0d0) then
                  is_bound = .false.
               else
                  is_bound = .true.
               endif

               if (s% x_logical_ctrl(4) .and. (s% v(k)/s% csound(k) .gt. 1.0d0)) is_bound = .false.

               !In order to remove material, the shell should have energetically unbound AND supersonic
               if (s% x_logical_ctrl(5) .and. (s% v(k)/s% csound(k) .lt. 1.0d0)) is_bound = .true.

            end function is_bound


            real(dp) function logic2dbl(a)
               logical, intent(in) :: a

               if (a) then
                  logic2dbl = 1.d0
               else
                  logic2dbl = 0.d0
               end if
            end function logic2dbl

      end subroutine CE_remove_unbound_envelope



      subroutine calc_recombination_after_struct_burn_mix(id, dt, res)
         ! # Here we will consider only H and He in the calculation of ionazion energy. Let us that N_H, N_HI, N_HII are
         ! # the total number of H atoms, the  number of neutral H atoms and the number of ionized H atoms in a specific shell.
         ! # For the Helium, the corresponding numbers would be N_He (total number of He atoms), N_HeI (number of neutral He atoms),
         ! # N_HeII (number of singly ionized He atoms), and N_HeIII (number of doubly ionized He atoms). Also, lets define as
         ! # Q_H the average charge per hydrogen particle (in units of electron charge) in a shell and f_HI the fraction of neutral
         ! # H in a shell. For Helium the respective numbers are Q_He and F_HeI. Given these definitions we can write the equations:
         ! # For hydrogen
         ! # N_HI + N_HII = N_H
         ! # N_HI = f_HI * N_H
         ! # N_HII / NH = Q_H
         ! # For Helium these equations become:
         ! # N_HeI + N_HeII + N_HeIII = N_He
         ! # N_HeI = f_HeI * N_He
         ! # (N_HeII + 2*N_HeIII)/N_He = Q_He
         ! # Solving this simple system of equations can give us the number of H and He atoms at each ionization state:
         ! # N_HII = Q_H * NH
         ! # N_HeII = (2.-2.*f_HeI -Q_He) * N_He
         ! # N_HeIII = (Q_He + f_HeI -1.) * N_He

         integer, intent(in) :: id
         real(dp), intent(in) :: dt
         integer, intent(out) :: res ! keep_going, redo, retry, backup, terminate
         type (star_info), pointer :: s
         integer :: k, ierr
         real(dp) :: erg_in_ev, Eion_HII_pp, Eion_HeII_pp, Eion_HeIII_pp, e_charge, He_mass, H_mass
         real(dp) :: N_H, N_He, N_HII, N_HeII, N_HeIII
         real(dp) :: avg_charge_H, avg_charge_He, neutral_fraction_H, neutral_fraction_He

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         erg_in_ev =1.60217657d-12 !# ergs in 1 eV
         Eion_HII_pp = 13.5924 * erg_in_ev !# Ionization energy HI -> HII of one atom in ergs
         Eion_HeII_pp = 24.5874 * erg_in_ev !# Ionization energy HI -> HII of one atom in ergs
         Eion_HeIII_pp = 54.417760 * erg_in_ev !# Ionization energy HI -> HII of one atom in ergs
         e_charge = 1.60217657d-19 !# electron chanrge in coulomb
         He_mass = 6.64648d-24 !# Mass of He atom in gr
         H_mass = 1.67372d-24 !# Mass of H atom in gr

         s% xtra2_array = 0.d0

         do k = 1, s% nz
            N_H = s% dm(k) * s% X(k) / H_mass
            N_He = s% dm(k) * s% Y(k) / He_mass
            avg_charge_H = get_ion_info(ion_iZ_H,k)
            avg_charge_He = get_ion_info(ion_iZ_He,k)
            neutral_fraction_H = get_ion_info(ion_ifneut_H,k)
            neutral_fraction_He = get_ion_info(ion_ifneut_He,k)
            N_HII = avg_charge_H * N_H
            N_HeII = (2.d0-2.d0* neutral_fraction_He - avg_charge_He) * N_He
            N_HeIII = (avg_charge_He + neutral_fraction_He -1.d0) * N_He

            ! We save the specific energy stored in ionized H and He at the end of the timestep in xtra2_array
            s% xtra2_array(k) = (N_HII*Eion_HII_pp + N_HeII*Eion_HeII_pp + N_HeIII*(Eion_HeII_pp+Eion_HeIII_pp)) / s% dm(k)
            s% xtra2_array(k) = (s% xtra1_array(k) - s% xtra2_array(k))/ dt

            ! We save the difference in ionization energy for each species in xtra_arrays 3, 4, and 5
            s% xtra3_array(k) = (N_HII*Eion_HII_pp / s% dm(k) - s% xtra3_array(k)) / dt
            s% xtra4_array(k) = (N_HeII*Eion_HeII_pp / s% dm(k) - s% xtra4_array(k)) / dt
            s% xtra5_array(k) = (N_HeIII*(Eion_HeII_pp+Eion_HeIII_pp) / s% dm(k) - s% xtra5_array(k)) / dt

         end do
         !And now calculate the specific energy release from recombination and store it again in xtra2_array.

         res = keep_going

         contains

         real(dp) function get_ion_info(id,k)
            use ionization_lib, only: eval_ionization
            integer, intent(in) :: id, k
            integer :: ierr
            real(dp) :: ionization_res(num_ion_vals)

            ierr = 0

            call eval_ionization( &
               1d0 - (s% X(k) + s% Y(k)), s% X(k), s% Rho(k), s% lnd(k)/ln10, &
               s% T(k), s% lnT(k)/ln10, ionization_res, ierr)
            if (ierr /= 0) ionization_res = 0
            get_ion_info = ionization_res(id)

            if (ierr /= 0) stop "Error returned from subroutine eval_ionization"
         end function get_ion_info



      end subroutine calc_recombination_after_struct_burn_mix






      end module CE_after_struct_burn_mix
