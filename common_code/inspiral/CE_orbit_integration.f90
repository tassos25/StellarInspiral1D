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

    module orbit_cee

    use dop853_module, wp => dop853_wp
    use iso_fortran_env, only: output_unit
    use star_def
    use const_def
    use math_lib
    use CE_orbit, only: AtoP, TukeyWindow, calc_quantities_at_comp_position

    implicit none

    !real(8), DIMENSION(:), ALLOCATABLE :: logrho, logr, mass, csound
    !integer :: num_rows
    !real(8), parameter :: pi = 3.141592653589793d0
    !real(8), parameter :: G = 392512559.8496094d0   ! grav  cnt.        [ Rsun^3 / Msun / yr^2] 
    !real(8), parameter :: c0 = 13598865.132357053d0 ! light speed       [ Rsun / yr]

    integer :: id_dop ! these variables cannot be arguments of this function
    integer :: ierr_dop

    !real(8), parameter :: mcomp = 1.0d0             ! companion mass    [ Msun ]

    contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_mesa_id(id,ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr

        id_dop = id
        ierr_dop = ierr
    end subroutine read_mesa_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine kepler_polar(me,x,y,f)

        !! Right-hand side of van der Pol's equation

        !implicit none

        

        

        class(dop853_class),intent(inout) :: me
        real(8),intent(in)               :: x
        real(8),dimension(:),intent(in)  :: y
        real(8),dimension(:),intent(out) :: f
        real(8), parameter :: G = 392512559.8496094d0   ! grav  cnt.        [ Rsun^3 / Msun / yr^2] 
        real(8), parameter :: c0 = 13598865.132357053d0 ! light speed       [ Rsun / yr]

        
        
        type (star_info), pointer :: s

        

        
        real(8) :: M                                    ! donnor mass       [ Msun ]

        real(8) :: mu                                   ! grav. parameter
        real(8) :: PN25                                 ! 2.5 Post-Newt. terms
        real(8) :: Fext, Fr, Fnu                           ! External force, and its components 

        real(8) :: r,nu

        real(8) :: magv, renv, machn, bmax, bmin, lambda
        real(8) :: tmprho, tmp1, tmp2, tmp3

        real(8) :: CE_companion_mass

        
        
        !real(8), DIMENSION(:), ALLOCATABLE :: logrho, logr, mass, 
        real(8), DIMENSION(:), ALLOCATABLE :: arr
        real(8) :: rho, msunrsun3 = 0.16934021222434983d0, rsunyr = 0.0004536093143596379d0

        real(8) :: value
        INTEGER :: index
        
        integer :: k, k_bottom
        real(8) :: CE_energy_rate, CE_companion_position, CE_companion_radius
        real(8) :: CE_n_acc_radii
        real(8) :: M2, R2
        real(8) :: F_drag
        real(8) :: F_DHL, f1, f2, f3, e_rho
        real(8) :: mdot, mdot_macleod, mdot_HL, L_acc, a1, a2, a3, a4
        real(8) :: R_acc, R_acc_low, R_acc_high
        real(8) :: v_rel, beta, M_encl, csound
        real(8) :: rho_at_companion, scale_height_at_companion
        real(8) :: drag_factor, log_mdot_factor, lambda_squared
        real(8) :: F_drag_subsonic, F_drag_supersonic

         ierr_dop = 0
        call star_ptr(id_dop, s, ierr_dop)
        if (ierr_dop /= 0) return

         ! Alternative energy source here

         ! Get input controls
         CE_energy_rate = s% xtra(1)
         CE_companion_position = s% xtra(2)
         CE_companion_radius = s% xtra(3)
         CE_companion_mass = s% xtra(4)
         CE_n_acc_radii = s% xtra(5)

         !call calc_quantities_at_comp_position(id_dop, ierr_dop)

         R_acc =2.0 * standard_cgrav * M2 / &
             ((v_rel_at_companion*v_rel_at_companion) + csound_at_companion*csound_at_companion) 

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


         else if (s% x_integer_ctrl(2) == 1) then

            F_drag = pi * R_acc**2 * rho_at_companion * v_rel**2

            ! Add hydrodynamic drag
            F_drag = F_drag + pi * (CE_companion_radius * Rsun)**2.0 * rho_at_companion * v_rel**2


         else
            stop "x_integer_ctrl(2) is not defined"
         endif

        

        ! Companion mass 
        CE_companion_mass = s% xtra(4)

        ! r and nu
        r =  y(1)
        nu = y(3)

        ! dr/dt and dnu/dt
        f(1) = y(2)             !dr/dt
        f(3) = y(4)             !dnu/dt

        

        arr = s% r
        value = r * Rsun
        ! Call the subroutine to find the index at the location of the companion
        CALL FindApproximateValueIndex(arr, value, index)
        DEALLOCATE(arr)

        
        ! Eclosed mass at the position of the companion
        M = s% m(index) / Msun

        ! Gravitational parameter 
        mu = G * ( CE_companion_mass + M )

        ! External force and its components
        renv = s% r(1) / Rsun
        magv = sqrt( f(1)*f(1) + r*r*f(3)*f(3) )

        if (r.gt.renv) then
            Fext = 0.d0
        else
            if (.false.) then
                tmprho =  (s% rho(index)) * msunrsun3
                tmp1 = CE_companion_mass * M / (CE_companion_mass + M)
                tmp2 = 2.d0 * pi * G*G * CE_companion_mass * tmp1 * tmprho
                machn = magv / (s% csound(index) * rsunyr) ! => change this, add it to profile
                if (machn.lt. 1.d0)  then
                    tmp3 = log( (1.d0 + machn)/(1.d0 - machn) * exp(-2.d0 * machn) )
                else
                    bmax = 2.d0 * r
                    bmin = G*CE_companion_mass / ((magv - s% omega(index) * s% r(index) / Rsun * secyer)**2.d0)
                    lambda = bmax / bmin
                    tmp3 = log(lambda*lambda - (lambda*lambda / (machn*machn)))
                endif
                ! Dynamical drag
                Fext = - (tmp2 * tmp3 ) / (magv*magv)

                ! Adding hydrodynamical drag
                !Fext = Fext + 
            end if 
            Fext = - F_drag
        endif


        
        Fr  = Fext * (f(1) / magv)
        Fnu = Fext * (r * f(3) / magv) ! be careful with the units

        ! 2.5 Post Newtonian terms
        PN25 = (magv/c0)*(magv/c0) + 1.5d0 * (magv/c0)*(magv/c0)*(magv/c0)*(magv/c0)





        ! d2r/dt2 and d2nu/dt2
        f(2) = - mu / ( y(1)*y(1) ) * (1.0 + PN25) + Fr + r * f(3)*f(3) ! d2r/dt2

        f(4) = Fnu / r - (2.d0 * f(1) * f(3) / r)                             ! d2nu/dt2

    

    end subroutine kepler_polar



    

    SUBROUTINE FindApproximateValueIndex(arr, value, index)
        IMPLICIT NONE
        real(8), DIMENSION(:), INTENT(IN) :: arr
        real(8), INTENT(IN) :: value
        INTEGER, INTENT(OUT) :: index

        INTEGER :: low, high, mid

        ! Initialize the search boundaries
        low = 1
        high = SIZE(arr)
        index = 0

        ! Perform binary search
        DO WHILE (low <= high)
            mid = (low + high) / 2
            IF (arr(mid) >= value) THEN
            index = mid  ! Store the current index
            low = mid + 1
            ELSE
            high = mid - 1
            END IF
        END DO

        ! Handle the case when the value is not found
        IF (index == 0) THEN
            index = 1
        END IF
        END SUBROUTINE FindApproximateValueIndex



    end module orbit_cee


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



