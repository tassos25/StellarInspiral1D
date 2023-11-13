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

    implicit none

    !real(8), DIMENSION(:), ALLOCATABLE :: logrho, logr, mass, csound
    !integer :: num_rows
    !real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: G = 392512559.8496094d0   ! grav  cnt.        [ Rsun^3 / Msun / yr^2] 
    real(8), parameter :: c0 = 13598865.132357053d0 ! light speed       [ Rsun / yr]

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
        !real(8), DIMENSION(:) :: arr

        ierr_dop = 0
        call star_ptr(id_dop, s, ierr_dop)
        if (ierr_dop /= 0) return

        ! Companion mass 
        CE_companion_mass = s% xtra(4)

        ! r and nu
        r =  y(1)
        nu = y(3)

        ! dr/dt and dnu/dt
        f(1) = y(2)             !dr/dt
        f(3) = y(4)             !dnu/dt

        

        arr = s% r
        value = r
        ! Call the subroutine to find the index at the location of the companion
        CALL FindApproximateValueIndex(arr, value, index)
        DEALLOCATE(arr)

        
        ! Eclosed mass at the position of the companion
        M = s% m(index) / Msun

        ! Gravitational parameter 
        mu = G * ( CE_companion_mass + M )

        ! External force and its components
        renv = s% r(1)
        magv = sqrt( f(1)*f(1) + r*r*f(3)*f(3) )

        if (r.gt.renv) then
            Fext = 0.d0
        else
            tmprho =  (s% rho(index)) * msunrsun3
            tmp1 = CE_companion_mass * M / (CE_companion_mass + M)
            tmp2 = 2.d0 * pi * G*G * CE_companion_mass * tmp1 * tmprho
            machn = magv / (s% csound(index) * rsunyr) ! => change this, add it to profile
            if (machn.lt. 1.d0)  then
                tmp3 = log( (1.d0 + machn)/(1.d0 - machn) * exp(-2.d0 * machn) )
            else
                bmax = 2.d0 * r
                bmin = G*CE_companion_mass / ((magv - s% omega(index) * s% r(index) / Rsun)**2.d0)
                lambda = bmax / bmin
                tmp3 = log(lambda*lambda - (lambda*lambda / (machn*machn)))
            endif
            ! Dynamical drag
            Fext = - (tmp2 * tmp3 ) / (magv*magv)

            ! Adding hydrodynamical drag
            !Fext = Fext + 
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



