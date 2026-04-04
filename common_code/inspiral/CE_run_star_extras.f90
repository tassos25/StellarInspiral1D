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

! copied the content from star/job/standard_run_star_extras.inc here and made the necessary changes

      module CE_run_star_extras

      use star_lib
      use run_star_support, only: failed
      use star_def
      use const_def
      use const_lib
      use colors_lib
      use colors_def
      
      use ionization_lib

      ! Add here all the external modules for CE_mesa here
      use CE_orbit
      use CE_energy
      use CE_torque
      use CE_after_struct_burn_mix
      use CE_before_struct_burn_mix
      use CE_adjust_mdot
      use CE_timestep
      use eos_lib, only: eosDT_get, eosDT_get_T

      implicit none

      contains

! ***********************************************************************
      subroutine CE_extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => CE_extras_startup
         !s% extras_start_step => CE_extras_start_step
         s% extras_check_model => CE_extras_check_model
         s% extras_finish_step => CE_extras_finish_step
         s% extras_after_evolve => CE_extras_after_evolve
         s% how_many_extra_history_columns => CE_how_many_extra_history_columns
         s% data_for_extra_history_columns => CE_data_for_extra_history_columns
         s% how_many_extra_profile_columns => CE_how_many_extra_profile_columns
         s% data_for_extra_profile_columns => CE_data_for_extra_profile_columns

         !s% how_many_extra_history_header_items => CE_how_many_extra_history_header_items
         !s% data_for_extra_history_header_items => CE_data_for_extra_history_header_items
         !s% how_many_extra_profile_header_items => CE_how_many_extra_profile_header_items
         !s% data_for_extra_profile_header_items => CE_data_for_extra_profile_header_items

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         s% job% warn_run_star_extras =.true.

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy (see star_data.inc)

         ! Here we should point to the names of the "other_" functions to be used
         s% other_energy => CE_inject_energy
         if (s% x_integer_ctrl(3)==1) then
            s% other_torque => CE_inject_am
         else if (s% x_integer_ctrl(3)==2) then
            s% other_torque => CE_inject_am2
         else
            stop "s% x_integer_ctrl(3) is not defined"
         endif
         s% other_before_struct_burn_mix => calc_recombination_before_struct_burn_mix
         s% other_after_struct_burn_mix => CE_other_after_struct_burn_mix
         s% other_adjust_mdot => CE_other_adjust_mdot
         !s% other_eosDT_get => my_other_eosDT_get      !does not exist anymore, see old_hooks.f90
         !s% other_eosDT_get_T => my_other_eosDT_get_T  !does not exist anymore, see old_hooks.f90

         ! Reading values of parameters from the extra controls that we are using
         ! Note that "extra_heat" is the specific energy added to the the cell in units of erg/s/gr

         !s% xtra(1) -> CE_energy_rate. It is initially set to 0. It will be calculated when CE_energy is called
         s% xtra(1) = 0.0d0
         !s% xtra(3) -> CE_companion_radius
         s% xtra(3) = s% x_ctrl(3)
         !s% xtra(4) -> CE_companion_mass
         s% xtra(4) = s% x_ctrl(4)
         !s% xtra(5) -> CE_n_acc_radii
         s% xtra(5) = s% x_ctrl(5)
         !s% xtra(6) -> CE_torque. It is initially set to 0. It will be calculated when CE_torque is called
         s% xtra(6) = 0.0d0
         !s% xtra(7) -> CE_mdot. It is initially set to 0. It will be calculated when CE_adjust_mdot is called
         s% xtra(7) = 0.0d0
         !s% xtra(20) -> L_acc. Accretion luminosity, calculated in CE_energy
         s% xtra(20) = 0.0
         !s% xtra(21) -> eta_pulse_wind: Wind enhancement term from Yoon & Cantiello (2010)
         s% xtra(21) = 1.
         !s% xtra(22) -> mdot_HL. Hoyle-Littleton accretion rate onto compact object. Calculated in CE_energy_rate
         s% xtra(22) = 0.0d0
         !s% xtra(23) -> mdot_macleod. Macleod & Ramirez-Ruiz accretion rate onto compact object. Calculated in CE_energy_rate
         s% xtra(23) = 0.0d0
         !s% xtra(24) -> temperature of the star at the position of the companion
         s% xtra(24) = 0.0d0
         !s% xtra(25) -> angular position of the companion (rad)
         s% xtra(25) = 0.0d0
         !s% xtra(26) -> radial component of the companion's velocity (Rsun/yr)
         s% xtra(26) = 0.0d0
         !s% xtra(27) -> angular component of the companion's velocity (rad/yr)
         s% xtra(27) = 0.0d0
         !s% xtra(28) -> radial component of the companion's acceleration (Rsun/yr^2)
         s% xtra(28) = 0.0d0
         !s% xtra(29) -> angular component of the companion's acceleration (rad/yr^2)
         s% xtra(29) = 0.0d0
         !s% xtra(30) -> magnitude of the drag force acting on the companion (Msun Rsun / yr^2)
         s% xtra(30) = 0.0d0


         !s% xtra(7) -> CE_test_case
         s% ixtra(1) = s% x_integer_ctrl(1)

         ! s% job% relax_omega = .true.
         ! s% job% new_omega = s% x_ctrl(15) * 2.*pi/AtoP(1.496112*Msun,s% xtra(4)*Msun,s% xtra(2)*Rsun)
         ! write(*,*) s% job% new_omega, s%x_ctrl(15), 1.496112*Msun,s% xtra(4)*Msun,s% xtra(2)*Rsun
         ! ! ! We set a very small timestep during the relaxation phase, so that the star does not evolve significantly
         ! ! s% job% relax_omega_max_yrs_dt = 1d-8
         !  s% job% set_initial_dt = .True.
         !  s% job% years_for_initial_dt = 1d-8

         
      end subroutine CE_extras_controls


! ***********************************************************************
      subroutine CE_extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_companion_position, R_acc, CE_n_acc_radii
         integer :: CE_test_case
         integer, parameter :: n_colors=11
         logical, parameter :: use_cache=.true.
         ierr = 0
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call ionization_init('ion', '', &
            '../common_code/ionization_data/cache', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'ionization_init failed during initialization'
            return
         end if

         call colors_init(1,(/trim(mesa_dir)//'/colors/data/lcb98cor.dat'/),(/n_colors/),ierr)

         if (ierr /= 0) then
            write(*,*) 'colors_init failed during initialization'
            return
         end if

         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if

         !If we are restarting from a photo, the rest of the synchronization and relaxing steps should be skipped
         if (restart) then

            s% job% set_initial_model_number = .false.

            s% job% change_v_flag = .true.
            s% job% change_initial_v_flag = .false.
            s% job% new_v_flag = .true.

            s% job% new_rotation_flag = .false.
            s% job% change_rotation_flag = .false.

            s% job% set_initial_age = .false.
            s% job% set_initial_model_number = .false.

            return
         endif

         ! Reading values of parameters from the extra controls that we are using
         ! Note that "extra_heat" is the specific energy added to the the  cell in units of erg/s/gr

         !s% xtra(1) -> CE_energy_rate. It is initially set to 0. It will be calculated when CE_energy is called
         s% xtra(1) = 0.0d0
         !s% xtra(3) -> CE_companion_radius
         s% xtra(3) = s% x_ctrl(3)
         !s% xtra(4) -> CE_companion_mass
         s% xtra(4) = s% x_ctrl(4)
         !s% xtra(5) -> CE_n_acc_radii
         s% xtra(5) = s% x_ctrl(5)
         !s% xtra(6) -> CE_torque. It is initially set to 0. It will be calculated when CE_torque is called
         s% xtra(6) = 0.0d0
         !s% xtra(7) -> CE_mdot. It is initially set to 0. It will be calculated when CE_adjust_mdot is called
         s% xtra(7) = 0.0d0

         !s% xtra(7) -> CE_test_case
         s% ixtra(1) = s% x_integer_ctrl(1)
         !s% xtra(2) -> CE_companion_position = CE_companion_initial_position * Rsatr
         s% xtra(2) = s% x_ctrl(2) * s% r(1) / Rsun

         ! s% lxtra(1) -> has the binary merged
         s% lxtra(1) = .false.

            !s% job% relax_omega = .true. ! don't change omega
            s% job% new_omega = s% x_ctrl(15) * 2.*pi/AtoP(s% m(1),s% xtra(4)*Msun,s% xtra(2) * Rsun)
            ! We set a very small timestep during the relaxation phase, so that the star does not evolve significantly
            s% job% relax_omega_max_yrs_dt = 1d-8
            !s% job% set_initial_dt = .True. ! moved to inlist
            !s% job% years_for_initial_dt = 3.78d-8 ! moved to inlist

         ! We are calling here the relax_omega, because we want to first have loaded the model so that we know its radius, and mass.
         if (s% rotation_flag .and. s% job% relax_omega) then
            write(*,*) 'new_omega =', s% job% new_omega
            call star_relax_uniform_omega( &
               id, 0, s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            s% job% relax_omega = .false.
         else
            !call star_relax_num_steps(id, 100, 1d-8 * secyer, ierr)
         endif

         !After relaxation is done, the timestep automatically increases to a "large" timestep. Here we are trying to make this
         !transition smoother
         if (s% job% years_for_initial_dt>0) then
            s% dt_next = s% job% years_for_initial_dt * secyer ! copy from initial step
         endif

         CE_companion_position = s% xtra(2)
         CE_test_case = s% ixtra(1)
         CE_n_acc_radii = s% x_ctrl(5)
         call calc_quantities_at_comp_position(id, ierr)
         R_acc = s% xtra(12)

         ! We need to increase the resolution around the area where the extra heat is deposited
         ! We will do this at the startup and also in the extra_check model, since the position
         ! of the companion will be changing
         if (CE_test_case == 2 .or. CE_test_case == 3 .or. CE_test_case == 4) then
            s% R_function2_param1 = CE_companion_position/(s%r(1)/Rsun) + 2.0* CE_n_acc_radii * R_acc/s%r(1)
            s% R_function2_param2 = CE_companion_position/(s%r(1)/Rsun) - 2.0* CE_n_acc_radii * R_acc/s%r(1)
         endif
      end subroutine CE_extras_startup


! ***********************************************************************
      integer function CE_extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         CE_extras_start_step = 0
      end function CE_extras_start_step


! ***********************************************************************
      integer function CE_extras_check_model(id)
      ! returns either keep_going, retry, backup, or terminate.
         integer, intent(in) :: id
         integer :: ierr, result
         type (star_info), pointer :: s
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_ang_mom_transferred, R_acc, CE_n_acc_radii
         integer :: CE_test_case

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         result = keep_going

         ! Reading initial values of parameters from the extra controls that we are using
         ! Note that "extra_heat" is the specific energy added to the the  cell in units of erg/s/gr
         CE_energy_rate = s% xtra(1)
         CE_companion_position = s% xtra(2)
         CE_companion_radius = s% xtra(3)
         CE_companion_mass = s% xtra(4)
         CE_ang_mom_transferred = s% xtra(6)
         CE_test_case = s% ixtra(1)

         call calc_quantities_at_comp_position(id, ierr)
         R_acc = s% xtra(12)

         ! We need to increase the resolution around the area where the extra heat is deposited
         ! We will do this at the startup and also in the extra_check model, since the position
         ! of the companion will be changing
         CE_n_acc_radii = s% x_ctrl(5)
         if (CE_test_case == 2 .or. CE_test_case == 3 .or. CE_test_case == 4) then
            s% R_function2_param1 = CE_companion_position/(s%r(1)/Rsun) + 2.0* CE_n_acc_radii * R_acc/s%r(1)
            s% R_function2_param2 = CE_companion_position/(s%r(1)/Rsun) - 2.0* CE_n_acc_radii * R_acc/s%r(1)
         endif

         ! For test cases 1 and 2 (heating of the whole envelope and of the base of the envelope) the code below must be skipped
         if (s% x_integer_ctrl(1) .ne. 1 .and. s% x_integer_ctrl(1) .ne. 2) then
            ! Adjust orbital separation based on energy deposited
            ! unless system has merged
            if (.not. s% lxtra(1)) call CE_orbit_adjust(id, ierr)
            ! Added timestep controls
            result = worst_result(result, CE_pick_next_timestep(s))
         endif

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (result == terminate) s% termination_code = t_extras_check_model

         CE_extras_check_model = result
      end function CE_extras_check_model


! ***********************************************************************
      integer function CE_extras_finish_step(id)
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         CE_extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (CE_extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function CE_extras_finish_step


! ***********************************************************************
      subroutine CE_extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !Saving those values back to xtra_controls so that restarts work.
         s% x_ctrl(2) = s% xtra(2) / (s% r(1) / Rsun)
         s% x_ctrl(3) = s% xtra(3)
         s% x_ctrl(4) = s% xtra(4)
      end subroutine CE_extras_after_evolve


! ***********************************************************************
      integer function CE_how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         CE_how_many_extra_history_columns = 26
      end function CE_how_many_extra_history_columns


! ***********************************************************************
      subroutine CE_data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: b_mag, u_mag, v_mag, r_mag, i_mag

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Call function to get rs
         !call get_magnitudes(id, b_mag, u_mag, v_mag, r_mag, i_mag, ierr)

         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         names(1) = 'CE_energy_rate'
         vals(1) = s% xtra(1)
         names(2) = 'L_accretion'
         vals(2) = s% xtra(20)
         names(3) = 'Mdot_macleod'
         vals(3) = s% xtra(23)
         names(4) = 'Mdot_HL'
         vals(4) = s% xtra(22)
         names(5) = 'CE_torque'
         vals(5) = s% xtra(6)
         names(6) = 'CE_companion_position_r'
         vals(6) = s% xtra(2)
         names(7) = 'CE_companion_position_m'
         vals(7) = s% xtra(9)
         names(8) = 'CE_ang_mom_transferred'
         vals(8) = s% xtra(6)
         names(9) = 'envelope_binding_energy'
         vals(9) = s% xtra(11)
         names(10) = 'R_acc'
         vals(10) = s% xtra(12)
         names(11) = 'R_acc_low'
         vals(11) = s% xtra(13)
         names(12) = 'R_acc_high'
         vals(12) = s% xtra(14)
         names(13) = 'v_rel'
         vals(13) = s% xtra(16)
         names(14) = 'v_over_c_sound'
         vals(14) = s% xtra(17)
         names(15) = 'eta_pulse_wind' ! From Yoon & Cantiello (2010)
         vals(15) = s% xtra(21)
         names(16) = 'b_mag'
         vals(16) = b_mag
         names(17) = 'u_mag'
         vals(17) = u_mag
         names(18) = 'v_mag'
         vals(18) = v_mag
         names(19) = 'r_mag'
         vals(19) = r_mag
         names(20) = 'i_mag'
         vals(20) = i_mag

         names(21) = 'CE_companion_position_nu'
         vals(21) = s% xtra(25)
         names(22) = 'CE_companion_velocity_r'
         vals(22) = s% xtra(26)
         names(23) = 'CE_companion_velocity_nu'
         vals(23) = s% xtra(27)
         names(24) = 'CE_companion_acceleration_r'
         vals(24) = s% xtra(28)
         names(25) = 'CE_companion_acceleration_nu'
         vals(25) = s% xtra(29)
         names(26) = 'CE_companion_drag_force_magnitude'
         vals(26) = s% xtra(30)

         ! If a distance provided, adjust from absolute to apparent magnitude
         if (s% x_ctrl(17) .ne. -1) then
            vals(16) = vals(16) + 5.0*(log10(s% x_ctrl(17) * 1000.0) - 1.0)
            vals(17) = vals(17) + 5.0*(log10(s% x_ctrl(17) * 1000.0) - 1.0)
            vals(18) = vals(18) + 5.0*(log10(s% x_ctrl(17) * 1000.0) - 1.0)
            vals(19) = vals(19) + 5.0*(log10(s% x_ctrl(17) * 1000.0) - 1.0)
            vals(20) = vals(20) + 5.0*(log10(s% x_ctrl(17) * 1000.0) - 1.0)
         endif
      end subroutine CE_data_for_extra_history_columns


! ***********************************************************************
      integer function CE_how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         CE_how_many_extra_profile_columns = 3
     ! previously have been 5, but the last three don't exist anymore
      end function CE_how_many_extra_profile_columns


! ***********************************************************************
      subroutine CE_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'CE_data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

         names(1) = 'ionization_energy'
         names(2) = 'eps_recombination'
         names(3) = 'CE_extra_heat'
      !   names(3) = 'eps_visc'
      !   names(4) = 'eta_visc'
      !   names(5) = 'Qvisc'
         do k = 1, nz
           vals(k,1) = s% xtra1_array(k)
           vals(k,2) = s% xtra2_array(k)
           vals(k,3) = s% xtra6_array(k)
      !     vals(k,3) = s% eps_visc(k)  ! does not exist anymore
      !     vals(k,4) = s% eta_visc(k)  ! does not exist anymore
      !     vals(k,5) = s% Qvisc(k)     ! does not exist anymore
         end do
      end subroutine CE_data_for_extra_profile_columns


! ***********************************************************************
      integer function CE_how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         CE_how_many_extra_history_header_items = 0
      end function CE_how_many_extra_history_header_items


! ***********************************************************************
      subroutine CE_data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine CE_data_for_extra_history_header_items


! ***********************************************************************
      integer function CE_how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         CE_how_many_extra_profile_header_items = 0
      end function CE_how_many_extra_profile_header_items


! ***********************************************************************
      subroutine CE_data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine CE_data_for_extra_profile_header_items


! ***********************************************************************
      subroutine get_magnitudes(id,b_Mag,u_Mag,v_Mag,r_Mag,i_Mag,ierr)
         use chem_def, only: zsol
         integer, intent(in) :: id
         type (star_info), pointer :: s
         real(dp), intent(out) :: b_Mag,u_Mag,v_Mag,r_Mag,i_Mag
         integer, intent(out) :: ierr
         real(dp)  :: log_Teff ! log10 of surface temp
         real(dp)  :: Fe_H ! [Fe/H]
         real(dp) :: log_g
         real(dp) :: lum ! in solar luminosities
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         log_Teff = log10(s% Teff)
         Fe_H = safe_log10(get_current_z_at_point(id, 1, ierr) / zsol)
         log_g = safe_log10(s% grav(1))

         if (ierr /= 0) then
            write(*,*) 'failed in get_current_z_at_point'
            return
         end if

         ! Absolute magnitudes from bolometric with corrections
         lum = 10.d0 ** s% log_surface_luminosity

         b_Mag = get_abs_mag_by_name('b',log_Teff,log_g, Fe_H,lum, ierr)
         u_Mag = get_abs_mag_by_name('u',log_Teff,log_g, Fe_H,lum, ierr)
         v_Mag = get_abs_mag_by_name('v',log_Teff,log_g, Fe_H,lum, ierr)
         r_Mag = get_abs_mag_by_name('r',log_Teff,log_g, Fe_H,lum, ierr)
         i_Mag = get_abs_mag_by_name('i',log_Teff,log_g, Fe_H,lum, ierr)

         if (ierr /= 0) then
            write(*,*) 'failed in colors_get'
            return
         end if
      end subroutine get_magnitudes


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3

! ***********************************************************************
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


! ***********************************************************************
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


! ***********************************************************************
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


! ***********************************************************************
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      end subroutine move_extra_info


      end module CE_run_star_extras
