!*****************************************************************************************
!>
!  2-body problem at 2.5 Post-Newtonian terms + external force 
!
!  edited by Jaime Roman-Garza
!

    module orbit_cee

    use dop853_module, wp => dop853_wp
    use iso_fortran_env, only: output_unit

    implicit none

    real(8), DIMENSION(:), ALLOCATABLE :: logrho, logr, mass, csound
    integer :: num_rows
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: G = 392512559.8496094d0   ! grav  cnt.        [ Rsun^3 / Msun / yr^2] 
    real(8), parameter :: c0 = 13598865.132357053d0 ! light speed       [ Rsun / yr]

    real(8), parameter :: mcomp = 1.0d0             ! companion mass    [ Msun ]

    contains 
    
    SUBROUTINE ReadProfile(filename, logrho, logr, mass, csound)

    real(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: logrho, logr, mass, csound
    CHARACTER(*), INTENT(IN) :: filename
    

    INTEGER :: i, num_rows
    

    OPEN(UNIT=11, FILE=filename, STATUS='old')
    
    ! Determine the number of rows in the file (assuming space-separated columns)
    num_rows = 0
    DO
        READ(11, *, IOSTAT=i) logrho, logr, mass, csound
        IF (i /= 0) EXIT
        num_rows = num_rows + 1
    END DO
    
    ! Allocate memory for the arrays
    ALLOCATE(logrho(num_rows), logr(num_rows), mass(num_rows) , csound(num_rows))
    
    ! Rewind the file to the beginning
    REWIND(11)
    
    ! Read the data and assign to arrays
    DO i = 1, num_rows
        READ(11, *) logrho(i), logr(i), mass(i) , csound(i)
    END DO
    
    CLOSE(11)
    
    END SUBROUTINE ReadProfile

        subroutine fvpol(me,x,y,f)

    ! Test, second order equation as a system of 
    !       two 1st order eqs.
    !       y'' = 1  

    implicit none

    class(dop853_class),intent(inout) :: me
    real(8),intent(in)               :: x
    real(8),dimension(:),intent(in)  :: y
    real(8),dimension(:),intent(out) :: f

    real(8), parameter :: a = 1.0

    !real(wp),parameter :: eps = 1.0d-3
    
    f(1) = y(2) !((1-y(1)**2)*y(2)-y(1))/eps  ! y'(x)
    f(2) = a !y(2)                            ! y''(x)
    

    end subroutine fvpol



    subroutine kepler_polar(me,x,y,f)

    !! Right-hand side of van der Pol's equation

    implicit none

    class(dop853_class),intent(inout) :: me
    real(8),intent(in)               :: x
    real(8),dimension(:),intent(in)  :: y
    real(8),dimension(:),intent(out) :: f



    
    real(8) :: M                                    ! donnor mass       [ Msun ]

    real(8) :: mu                                   ! grav. parameter
    real(8) :: PN25                                 ! 2.5 Post-Newt. terms
    real(8) :: Fext, Fr, Fnu                           ! External force, and its components 

    real(8) :: r,nu

    real(8) :: magv, renv, machn, bmax, bmin, lambda
    real(8) :: tmprho, tmp1, tmp2, tmp3

    
    
    !real(8), DIMENSION(:), ALLOCATABLE :: logrho, logr, mass, 
    real(8), DIMENSION(:), ALLOCATABLE :: arr
    real(8) :: rho, msuncm3 = 0.16934021222434983d0, rsunyr = 0.0004536093143596379d0

    real(8) :: value
    INTEGER :: index
    !real(8), DIMENSION(:) :: arr

    ! r and nu
    r =  y(1)
    nu = y(3)

    ! dr/dt and dnu/dt
    f(1) = y(2)             !dr/dt
    f(3) = y(4)             !dnu/dt

    

    arr = logr
    value = log10(r)
    ! Call the subroutine to find the index at the location of the companion
    CALL FindApproximateValueIndex(arr, value, index)
    DEALLOCATE(arr)

    
    ! Eclosed mass at the position of the companion
    M = mass(index)

    ! Gravitational parameter 
    mu = G * ( mcomp + M )

    ! External force and its components
    renv = 10.d0**(logr(1))
    magv = sqrt( f(1)*f(1) + r*r*f(3)*f(3) )

    if (r.gt.renv) then
        Fext = 0.d0
    else
        tmprho = (10.d0 ** (logrho(index))) * msuncm3
        tmp1 = mcomp * M / (mcomp + M)
        tmp2 = 2.d0 * pi * G*G * mcomp * tmp1 * tmprho
        machn = magv / (csound(index) * rsunyr) ! => change this, add it to profile
        if (machn.lt. 1.d0)  then
            tmp3 = log( (1.d0 + machn)/(1.d0 - machn) * exp(-2.d0 * machn) )
        else
            bmax = 2.d0 * r
            bmin = G*mcomp / (magv*magv)
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



