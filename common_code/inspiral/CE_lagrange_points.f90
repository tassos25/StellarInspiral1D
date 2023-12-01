! *********************************************************************************************************************************
!
!   This file is part of a mesa extension.
!   Authors of this file: Matthias U. Kruckow
!
! *********************************************************************************************************************************
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
! *********************************************************************************************************************************

      module CE_lagrange_points

      use const_def
      use math_lib

      implicit none

      contains

! *********************************************************************************************************************************
      logical function LP_inputs_valid(m1, m2, a, idx, method)
         ! input parameters:
         ! m1: mass of the more massive stellar component
         ! m2: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: m1, m2, a
         ! idx: number of index
         ! method: number of method
         integer, intent(in) :: idx, method
         
         if (m1 .le. 0) then
            write(*,*) "Error in LP_inputs_valid: negative mass: m1=", m1
            LP_inputs_valid = .false.
         else if (m2 .le. 0) then
            write(*,*) "Error in LP_inputs_valid: negative mass: m2=", m2
            LP_inputs_valid = .false.
         else if (a .le. 0) then
            write(*,*) "Error in LP_inputs_valid: negative separation: a=", a
            LP_inputs_valid = .false.
         else if (idx .le. 0) then
            write(*,*) "Error in LP_inputs_valid: unknown index: idx=", idx
            LP_inputs_valid = .false.
         else if (method .lt. 0) then
            write(*,*) "Error in LP_inputs_valid: unknown method: method=", method
            LP_inputs_valid = .false.
         else
            LP_inputs_valid = .true.
         endif
         
      end function LP_inputs_valid


! *********************************************************************************************************************************
      real(dp) function get_roche_lobe(m1, m2, a, idx, method)
         ! Roche lobe size for star of mass m1 with a companion of mass m2 at separation a
         ! input parameters:
         ! m1: mass of the more massive stellar component
         ! m2: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: m1, m2, a
         ! idx: index of the donor component
         ! method: method to use to get the value (0: bisectioning, >0: approximations)
         integer, intent(in) :: idx, method
         ! varibales:
         real(dp) :: q, rl
         
         if (.not. LP_inputs_valid(m1, m2, a, idx, method)) then
            get_roche_lobe = -1.0d0
            return
         else if (idx .gt. 2) then
            write(*,*) "Error in get_roche_lobe: only two stars in a binary: idx=", idx
            get_roche_lobe = -1.0d0
            return
         endif
         
         select case (idx)
         case (1)
            ! first component is donor
            q = m1/m2
         case (2)
            ! second component is donor
            q = m2/m1
         case default
            write(*,*) "Error in get_roche_lobe: unknown case of idx=", idx
            get_roche_lobe = -1.0d0
            return
         end select
         
         select case (method)
         case (0)
            ! estimate the volume with boxes (use 100^3 boxes)
            if (m1 .ge. m2) then
               get_roche_lobe = pow(get_LP_volume(m1, m2, a, 1, idx, method, 100)*3.0d0/(4.0d0*pi),one_third)
            else
               get_roche_lobe = pow(get_LP_volume(m2, m1, a, 1, 3-idx, method, 100)*3.0d0/(4.0d0*pi),one_third)
            endif
            return
         case (1)
            ! approximation of Eggleton, P. P., The Astrophysical Journal, vol. 268, pp. 368–369, 1983. doi:10.1086/160960.
            q = pow(q,one_third)
            rl = 0.49d0*q*q/(0.6d0*q*q + log1p(q))
         case (2)
            ! approximation of Paczyński, B., Annual Review of Astronomy and Astrophysics, vol. 9, p. 183, 1971.
            ! doi:10.1146/annurev.aa.09.090171.001151.
            if (q .ge. 20.0d0) then
               write(*,*) "Warning: use approximation of Paczyński outside of validity M1/M2=", q
            endif
            rl = max(0.46224d0*pow(q/(1.0d0+q),one_third), 0.38d0+0.2d0*log10(q))
         case default
            write(*,*) "Error in get_roche_lobe: unknown case of method=", method
            get_roche_lobe = -1.0d0
            return
         end select
         
         get_roche_lobe = rl*a
         
      end function get_roche_lobe


! *********************************************************************************************************************************
      real(dp) function get_L2_surface(m1, m2, a, idx, method)
         ! Roche lobe size for star of mass m1 with a companion of mass m2 at separation a
         ! input parameters:
         ! m1: mass of the more massive stellar component
         ! m2: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: m1, m2, a
         ! idx: index of the donor component
         ! method: method to use to get the value (0: bisectioning, >0: approximations)
         integer, intent(in) :: idx, method
         ! varibales:
         real(dp) :: q, rl2
         
         if (.not. LP_inputs_valid(m1, m2, a, idx, method)) then
            get_L2_surface = -1.0d0
            return
         else if (idx .gt. 2) then
            write(*,*) "Error in get_L2_surface: only two stars in a binary: idx=", idx
            get_L2_surface = -1.0d0
            return
         endif
         
         select case (idx)
         case (1)
            ! first component is donor
            q = m1/m2
         case (2)
            ! second component is donor
            q = m2/m1
         case default
            write(*,*) "Error in get_L2_surface: unknown case of idx=", idx
            get_L2_surface = -1.0d0
            return
         end select
         
         select case (method)
         case (0)
            ! estimate the volume with boxes (use 100^3 boxes)
            get_L2_surface = pow(get_LP_volume(m1, m2, a, 2, 0, method, 100)*3.0d0/(4.0d0*pi),one_third)
         case (1)
            ! approximation of Misra, D., Fragos, T., Tauris, T. M., Zapartas, E., and Aguilera-Dena, D. R., Astronomy and
            ! Astrophysics, vol. 642, id. A174, 2020. doi:10.1051/0004-6361/202038070.
            q = 1.0d0/q
            if (q .lt. 1.0d0) then
               rl2 = (0.784d0 * pow(q, 1.050d0) * exp(-0.188d0*q) + 1.004d0)
            else
               ! from paper
               !rl2 = (0.290d0 * pow(q, 0.829d0) * exp(-0.016d0*q) + 1.362d0)
               ! from POSYDON
               rl2 = (0.29066811d0 * pow(q, 0.82788069d0) * exp(-0.01572339d0*q) + 1.36176161d0)
            endif
            get_L2_surface = rl2*get_roche_lobe(m1, m2, a, idx, 1)
         case (2)
            ! approximation of Ge, H., Webbink, R. F., and Han, Z., The Astrophysical Journal Supplement Series, vol. 249, no. 1,
            ! id.9, 2020. doi:10.3847/1538-4365/ab98f6.
            if (q .le. 1.0d0) then
               rl2 = (0.179d0 + 0.01d0*q/(1.0d0+q)) * pow(q/(1.0d0+q),0.625d0)
            else
               rl2 = (0.179d0 + 0.01d0*q/(1.0d0+q) - 0.025*(q-1.0d0)/q) * pow(q/(1.0d0+q),0.625d0) * pow(q,-0.74d0)
            endif
            get_L2_surface = get_roche_lobe(m1, m2, a, idx, 1) + rl2*a
         case default
            write(*,*) "Error in get_L2_surface: unknown case of method=", method
            get_L2_surface = -1.0d0
            return
         end select
         
      end function get_L2_surface


! *********************************************************************************************************************************
      real(dp) function get_distance_CoM_comp(m1, m2, a, idx)
         ! input parameters:
         ! m1: mass of the more massive stellar component
         ! m2: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: m1, m2, a
         ! idx: index of the component
         integer, intent(in) :: idx
         
         if (.not. LP_inputs_valid(m1, m2, a, idx, 0)) then
            get_distance_CoM_comp = -1.0d0
            return
         else if (idx .gt. 2) then
            write(*,*) "Error in get_distance_CoM_comp: only two stars in a binary: idx=", idx
            get_distance_CoM_comp = -1.0d0
            return
         endif
         
         select case (idx)
         case (1)
            ! distance between center of mass and the more massive component
            ! a1 = m2/(m1+m2)*a
            get_distance_CoM_comp = m2/(m1+m2)*a
         case (2)
            ! distance between center of mass and the less massive component
            ! a2 = m1/(m1+m2)*a
            get_distance_CoM_comp = m1/(m1+m2)*a
         case default
            write(*,*) "Error in get_distance_CoM_comp: unknown case of idx=", idx
            get_distance_CoM_comp = -1.0d0
            return
         end select
         
      end function get_distance_CoM_comp


! *********************************************************************************************************************************
      real(dp) function get_distance_CoM_LP(mgtr, mlow, a, idx, method)
         ! input parameters:
         ! mgtr: mass of the more massive stellar component
         ! mlow: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: mgtr, mlow, a
         ! idx: number of the lagrange point of interest
         ! method: method to use to get the value (0: bisectioning, >0: approximations)
         integer, intent(in) :: idx, method
         ! varibales:
         real(dp) :: m1, m2, q, d_L2_m1, d_L2_m2
         
         if (.not. LP_inputs_valid(mgtr, mlow, a, idx, method)) then
            get_distance_CoM_LP = -1.0d0
            return
         else if (idx .gt. 3) then
            write(*,*) "Error in get_distance_CoM_LP: only the first 3 lagrange points are supported: idx=", idx
            get_distance_CoM_LP = -1.0d0
            return
         endif
         
         if (mlow .gt. mgtr) then
            write(*,*) "Warning: mlow=", mlow, ">mgtr=", mgtr, " => swap masses"
            m1 = mlow
            m2 = mgtr
         else
            m1 = mgtr
            m2 = mlow
         endif
         
         select case (idx)
         case (1)
            ! distance between center of mass and L1, which is between the components
            ! G*m1/x1^2 = G*m2/x2^2 + G*(m1+m2)*x/a^3, with x1=x+a1, x2=a2-x, a2>x, a1|2 = m2|1/(m1+m2)*a
            select case (method)
            case (0)
               ! solve equation with bisectioning
               ! minimum=0.0: L1 coincide with center of mass (in case of equal masses)
               ! maximum=1.0: L1 at position of less massive component, while center of mass coincide with more massive component
               ! (in case of extrem mass ratios)
               get_distance_CoM_LP = bisec(m1,m2,a,idx,0.0d0,1.0d0,1.0d-9)
            case default
               write(*,*) "Error in get_distance_CoM_LP: unknown case of method=", method, "for idx=", idx
               get_distance_CoM_LP = -1.0d0
               return
            end select
         case (2)
            ! distance between center of mass and L2, which is behind the less massive component
            ! G*(m1+m2)*x/a^3 = G*m1/x1^2 + G*m2/x2^2, with x1=x+a1, x2=x-a2, x>a2, a1|2 = m2|1/(m1+m2)*a
            select case (method)
            case (0)
               ! solve equation with bisectioning
               ! minimum=1.0: L2 at position of less massive component, while center of mass coincide with more massive component
               ! (in case of extrem mass ratios)
               ! maximum=1.5: 2.0 away would be unbound and center of mass is displaced by 0.5 (in case of equal masses)
               get_distance_CoM_LP = bisec(m1,m2,a,idx,1.0d0,1.5d0,1.0d-9)
            case (1)
               ! approximation (q<1) of Misra, D., Fragos, T., Tauris, T. M., Zapartas, E., and Aguilera-Dena, D. R., Astronomy and
               ! Astrophysics, vol. 642, id. A174, 2020. doi:10.1051/0004-6361/202038070.
               q = m2/m1
               d_L2_m1 = (3.334d0 * pow(q, 0.514d0) * exp(-0.052d0*q) + 1.308d0)
               d_L2_m1 = d_L2_m1 * get_roche_lobe(m1, m2, a, 1, 1)
               get_distance_CoM_LP = d_L2_m1 - get_distance_CoM_comp(m1, m2, a, 1)
            case (2)
               ! approximation (q>1) of Misra, D., Fragos, T., Tauris, T. M., Zapartas, E., and Aguilera-Dena, D. R., Astronomy and
               ! Astrophysics, vol. 642, id. A174, 2020. doi:10.1051/0004-6361/202038070.
               q = m1/m2
               ! from paper
               !d_L2_m2 = (-0.040d0 * pow(q, 0.866d0) * exp(-0.040d0*q) + 1.883d0)
               ! from POSYDON
               d_L2_m2 = (-0.04029713d0 * pow(q, 0.862143d0) * exp(-0.04049814d0*q) + 1.88325644d0)
               d_L2_m2 = d_L2_m2 * get_roche_lobe(m1, m2, a, 2, 1)
               get_distance_CoM_LP = d_L2_m2 + get_distance_CoM_comp(m1, m2, a, 2)
            case default
               write(*,*) "Error in get_distance_CoM_LP: unknown case of method=", method, "for idx=", idx
               get_distance_CoM_LP = -1.0d0
               return
            end select
         case (3)
            ! distance between center of mass and L3, which is behind the more massive component
            ! G*(m1+m2)*x/a^3 = G*m1/x1^2 + G*m2/x2^2, with x1=x-a1, x2=x+a2, x>a1, a1|2 = m2|1/(m1+m2)*a
            select case (method)
            case (0)
               ! solve equation with bisectioning
               ! minimum=1.0: L3 at opposite position of less massive component, while center of mass coincide with more massive
               ! component (in case of extrem mass ratios)
               ! maximum=1.5: 2.0 away would be unbound and center of mass is displaced by 0.5 (in case of equal masses)
               get_distance_CoM_LP = bisec(m1,m2,a,idx,1.0d0,1.5d0,1.0d-9)
            case default
               write(*,*) "Error in get_distance_CoM_LP: unknown case of method=", method, "for idx=", idx
               get_distance_CoM_LP = -1.0d0
               return
            end select
         case default
            write(*,*) "Error in get_distance_CoM_LP: unknown case of idx=", idx
            get_distance_CoM_LP = -1.0d0
            return
         end select
         
         contains
         
         real(dp) function resL1(x, m)
            real(dp), intent(in) :: x, m
            real(dp) :: x2, x3, x4, x5, m2, m3, m4
            if ((m .gt. 1.0d0).or.(m .lt. 0.5d0)) then
               ! m not in range
               write (*,*) "Error: m=", m, " out of range [0.5,1.0]"
            endif
            ! get powers of x and m
            x2 = x*x
            x3 = x2*x
            x4 = x3*x
            x5 = x4*x
            m2 = m*m
            m3 = m2*m
            m4 = m3*m
            ! residual of equation G*m1/x1^2 = G*m2/x2^2 + G*(m1+m2)*x/a^3, with x1=x+a1, x2=a2-x, a2>x, a1|2 = m2|1/(m1+m2)*a
            resL1 = x5 + (2.0d0-4.0d0*m)*x4 + (1.0d0-6.0d0*m+6.0d0*m2)*x3 + (1.0d0-4.0d0*m+6.0d0*m2-4.0d0*m3)*x2 &
                       + (2.0d0-4.0d0*m+5.0d0*m2-2.0d0*m3+m4)*x + (1.0d0-3.0d0*m+3.0d0*m2-2.0d0*m3)
            
         end function resL1
         
         real(dp) function resL2(x, m)
            real(dp), intent(in) :: x, m
            real(dp) :: x2, x3, x4, x5, m2, m3, m4
            if ((m .gt. 1.0d0).or.(m .lt. 0.5d0)) then
               ! m not in range
               write (*,*) "Error: m=", m, " out of range [0.5,1.0]"
            endif
            ! get powers of x and m
            x2 = x*x
            x3 = x2*x
            x4 = x3*x
            x5 = x4*x
            m2 = m*m
            m3 = m2*m
            m4 = m3*m
            ! residual of equation G*(m1+m2)*x/a^3 = G*m1/x1^2 + G*m2/x2^2, with x1=x+a1, x2=x-a2, x>a2, a1|2 = m2|1/(m1+m2)*a
            resL2 = x5 + (2.0d0-4.0d0*m)*x4 + (1.0d0-6.0d0*m+6.0d0*m2)*x3 + (-1.0d0-2.0d0*m+6.0d0*m2-4.0d0*m3)*x2 &
                       + (-2.0d0+4.0d0*m+m2-2.0d0*m3+m4)*x + (-1.0d0+3.0d0*m-3.0d0*m2)
            
         end function resL2
         
         real(dp) function resL3(x, m)
            real(dp), intent(in) :: x, m
            real(dp) :: x2, x3, x4, x5, m2, m3, m4
            if ((m .gt. 1.0d0).or.(m .lt. 0.5d0)) then
               ! m not in range
               write (*,*) "Error: m=", m, " out of range [0.5,1.0]"
            endif
            ! get powers of x and m
            x2 = x*x
            x3 = x2*x
            x4 = x3*x
            x5 = x4*x
            m2 = m*m
            m3 = m2*m
            m4 = m3*m
            ! residual of equation G*(m1+m2)*x/a^3 = G*m1/x1^2 + G*m2/x2^2, with x1=x-a1, x2=x+a2, x>a1, a1|2 = m2|1/(m1+m2)*a
            resL3 = x5 + (-2.0d0+4.0d0*m)*x4 + (1.0d0-6.0d0*m+6.0d0*m2)*x3 + (-1.0d0+2.0d0*m-6.0d0*m2+4.0d0*m3)*x2 &
                       + (2.0d0-4.0d0*m+m2-2.0d0*m3+m4)*x + (-1.0d0+3.0d0*m-3.0d0*m2)
            
         end function resL3

         real(dp) function bisec(m1, m2, a, idx, x1_in, x2_in, eps)
            real(dp), intent(in) :: m1, m2, a, x1_in, x2_in, eps
            integer, intent(in) :: idx
            real(dp) :: x, x1, x2, m, res, res1, res2
            if (m1 .ge. m2) then
               m = m1/(m1+m2)
            else
               m = m2/(m1+m2)
            endif
            x1 = x1_in
            x2 = x2_in
            select case (idx)
            case (1)
               res1 = resL1(x1, m)
               res2 = resL1(x2, m)
            case (2)
               res1 = resL2(x1, m)
               res2 = resL2(x2, m)
            case (3)
               res1 = resL3(x1, m)
               res2 = resL3(x2, m)
            case default
               res1 = 0.0d0
               res2 = 0.0d0
            end select
            if (res1*res2 .gt. 0.0d0) then
               write(*,*) "Warning: no bisectioning possible res(", x1, m, ")=", res1, " res(", x2, m, ")=", res2
               ! bisectioning not possible
               if (res1 .lt. res2) then
                  ! x1 has lower residual and is therefore assumed to be closer to the solution
                  x=x1
               else
                  ! x2 has lower residual and is therefore assumed to be closer to the solution
                  x=x2
               endif
            else
               ! do bisectioning: get new point in middle
               x = 0.5d0 * (x1+x2)
               do while (abs(x1-x2)>eps)
                  ! as long as precission isn't reached: get new residual
                  select case (idx)
                  case (1)
                     res = resL1(x, m)
                  case (2)
                     res = resL2(x, m)
                  case (3)
                     res = resL3(x, m)
                  case default
                     res = 0.0d0
                  end select
                  if (res1*res .gt. 0.0d0) then
                     ! solution between x and x2, hence update x1 and res1
                     x1 = x
                     res1 = res
                  else
                     ! solution between x1 and x, hence update x2 and res2
                     x2 = x
                     res2 = res
                  endif
                  ! get next point in middle
                  x = 0.5d0 * (x1+x2)
               enddo
            endif
            if (res1==0) then
               ! x1 is exact solution
               bisec = x1*a
            else if (res2==0) then
               ! x2 is exact solution
               bisec = x2*a
            else
               ! x is surely less then 0.5*eps away from solution, while x1 or x2 could be closer but unknown if so, hence x is
               ! best estimate
               bisec = x*a
            endif
            
         end function bisec
         
      end function get_distance_CoM_LP


! *********************************************************************************************************************************
      real(dp) function get_distance_LP_comp(m1, m2, a, idxLP, idxComp, method)
         ! input parameters:
         ! m1: mass of the more massive stellar component
         ! m2: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: m1, m2, a
         ! idxLP: number of the lagrange point of interest
         ! idxComp: index of the component
         integer, intent(in) :: idxLP, idxComp, method
         
         if (.not. LP_inputs_valid(m1, m2, a, idxLP, method)) then
            get_distance_LP_comp = -1.0d0
            return
         else if (idxLP .gt. 3) then
            write(*,*) "Error in get_distance_LP_comp: only the first 3 lagrange points are supported: idxLP=", idxLP
            get_distance_LP_comp = -1.0d0
            return
         else if ((idxComp .gt. 2).or.(idxComp .lt. 1)) then
            write(*,*) "Error in get_distance_LP_comp: only two stars in a binary: idxComp=", idxComp
            get_distance_LP_comp = -1.0d0
            return
         endif
         
         select case (idxComp)
         case (1)
            ! distance to m1
            if (m1 .ge. m2) then
               get_distance_LP_comp = get_distance_CoM_LP(m1, m2, a, idxLP, method)
               if (get_distance_LP_comp .lt. 0.0d0) then
                  ! Negative distance from get_distance_CoM_LP
                  return
               endif
               select case (idxLP)
               case (1)
                  ! distance to L1 (CoM between m1 and L1)
                  get_distance_LP_comp = get_distance_CoM_comp(m1, m2, a, idxComp) + get_distance_LP_comp
               case (2)
                  ! distance to L2 (CoM between m1 and L2)
                  get_distance_LP_comp = get_distance_LP_comp + get_distance_CoM_comp(m1, m2, a, idxComp)
               case (3)
                  ! distance to L3 (m1 between CoM and L3)
                  get_distance_LP_comp = get_distance_LP_comp - get_distance_CoM_comp(m1, m2, a, idxComp)
               case default
                  write(*,*) "Error in get_distance_LP_comp: unknown case of idxLP=", idxLP
                  get_distance_LP_comp = -1.0d0
                  return
               end select
            else
               get_distance_LP_comp = get_distance_CoM_LP(m2, m1, a, idxLP, method)
               if (get_distance_LP_comp .lt. 0.0d0) then
                  ! Negative distance from get_distance_CoM_LP
                  return
               endif
               select case (idxLP)
               case (1)
                  ! distance to L1 (L1 between CoM and m1)
                  get_distance_LP_comp = get_distance_CoM_comp(m1, m2, a, idxComp) - get_distance_LP_comp
               case (2)
                  ! distance to L2 (m1 between CoM and L2)
                  get_distance_LP_comp = get_distance_LP_comp - get_distance_CoM_comp(m1, m2, a, idxComp)
               case (3)
                  ! distance to L3 (CoM between m1 and L3)
                  get_distance_LP_comp = get_distance_LP_comp + get_distance_CoM_comp(m1, m2, a, idxComp)
               case default
                  write(*,*) "Error in get_distance_LP_comp: unknown case of idxLP=", idxLP
                  get_distance_LP_comp = -1.0d0
                  return
               end select
            endif
         case (2)
            ! distance to m2
            if (m2 .gt. m1) then
               get_distance_LP_comp = get_distance_CoM_LP(m2, m1, a, idxLP, method)
               if (get_distance_LP_comp .lt. 0.0d0) then
                  ! Negative distance from get_distance_CoM_LP
                  return
               endif
               select case (idxLP)
               case (1)
                  ! distance to L1 (CoM between m2 and L1)
                  get_distance_LP_comp = get_distance_CoM_comp(m1, m2, a, idxComp) + get_distance_LP_comp
               case (2)
                  ! distance to L2 (CoM between m2 and L2)
                  get_distance_LP_comp = get_distance_LP_comp + get_distance_CoM_comp(m1, m2, a, idxComp)
               case (3)
                  ! distance to L3 (m2 between CoM and L3)
                  get_distance_LP_comp = get_distance_LP_comp - get_distance_CoM_comp(m1, m2, a, idxComp)
               case default
                  write(*,*) "Error in get_distance_LP_comp: unknown case of idxLP=", idxLP
                  get_distance_LP_comp = -1.0d0
                  return
               end select
            else
               get_distance_LP_comp = get_distance_CoM_LP(m1, m2, a, idxLP, method)
               if (get_distance_LP_comp .lt. 0.0d0) then
                  ! Negative distance from get_distance_CoM_LP
                  return
               endif
               select case (idxLP)
               case (1)
                  ! distance to L1 (L1 between CoM and m2)
                  get_distance_LP_comp = get_distance_CoM_comp(m1, m2, a, idxComp) - get_distance_LP_comp
               case (2)
                  ! distance to L2 (m2 between CoM and L2)
                  get_distance_LP_comp = get_distance_LP_comp - get_distance_CoM_comp(m1, m2, a, idxComp)
               case (3)
                  ! distance to L3 (CoM between m2 and L3)
                  get_distance_LP_comp = get_distance_LP_comp + get_distance_CoM_comp(m1, m2, a, idxComp)
               case default
                  write(*,*) "Error in get_distance_LP_comp: unknown case of idxLP=", idxLP
                  get_distance_LP_comp = -1.0d0
                  return
               end select
            endif
         case default
            write(*,*) "Error in get_distance_LP_comp: unknown case of idxComp=", idxComp
            get_distance_LP_comp = -1.0d0
            return
         end select
         
      end function get_distance_LP_comp


! *********************************************************************************************************************************
      real(dp) function get_j_at_comp(mgtr, mlow, a, idx)
         ! input parameters:
         ! mgtr: mass of the more massive stellar component
         ! mlow: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: mgtr, mlow, a
         ! idx: index of the component
         integer, intent(in) :: idx
         ! varibales:
         real(dp) :: m1, m2, d_comp_CoM, M, omega
         
         if (.not. LP_inputs_valid(mgtr, mlow, a, idx, 0)) then
            get_j_at_comp = -1.0d0
            return
         else if (idx .gt. 2) then
            write(*,*) "Error in get_j_at_comp: only two stars in a binary: idx=", idx
            get_j_at_comp = -1.0d0
            return
         endif
         
         if (mlow .gt. mgtr) then
            write(*,*) "Warning: mlow=", mlow, ">mgtr=", mgtr, " => swap masses"
            m1 = mlow
            m2 = mgtr
         else
            m1 = mgtr
            m2 = mlow
         endif
         ! get total mass
         M = m1+m2
         ! get distance component to center of mass
         d_comp_CoM = get_distance_CoM_comp(m1, m2, a, idx)
         if (d_comp_CoM .le. 0) then
            get_j_at_comp = d_comp_CoM
            return
         endif
         ! get angular velocity (in circular orbit)
         omega = sqrt(standard_cgrav*M/(a*a*a))
         ! result specific angular momentum
         get_j_at_comp = d_comp_CoM * d_comp_CoM * omega
         
      end function get_j_at_comp


! *********************************************************************************************************************************
      real(dp) function get_j_at_LP(mgtr, mlow, a, idx, method)
         ! input parameters:
         ! mgtr: mass of the more massive stellar component
         ! mlow: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: mgtr, mlow, a
         ! idx: number of the lagrange point of interest
         ! method: method to use to get the value (0: bisectioning, >0: approximations)
         integer, intent(in) :: idx, method
         ! varibales:
         real(dp) :: m1, m2, d_LP_CoM, M, omega
         
         if (.not. LP_inputs_valid(mgtr, mlow, a, idx, method)) then
            get_j_at_LP = -1.0d0
            return
         else if (idx .gt. 3) then
            write(*,*) "Error in get_j_at_LP: only the first 3 lagrange points are supported: idx=", idx
            get_j_at_LP = -1.0d0
            return
         endif
         
         if (mlow .gt. mgtr) then
            write(*,*) "Warning: mlow=", mlow, ">mgtr=", mgtr, " => swap masses"
            m1 = mlow
            m2 = mgtr
         else
            m1 = mgtr
            m2 = mlow
         endif
         ! get total mass
         M = m1+m2
         ! get distance lagrange point to center of mass
         d_LP_CoM = get_distance_CoM_LP(m1, m2, a, idx, method)
         if (d_LP_CoM .le. 0) then
            get_j_at_LP = d_LP_CoM
            return
         endif
         ! get angular velocity (in circular orbit)
         omega = sqrt(standard_cgrav*M/(a*a*a))
         ! result specific angular momentum
         get_j_at_LP = d_LP_CoM * d_LP_CoM * omega
         
      end function get_j_at_LP


! *********************************************************************************************************************************
      real(dp) function get_LP_volume(mgtr, mlow, a, idxLP, idxComp, method, n)
         ! input parameters:
         ! mgtr: mass of the more massive stellar component
         ! mlow: mass of the less massive stellar component
         ! a: separation between the two binary components (in circular orbit)
         real(dp), intent(in) :: mgtr, mlow, a
         ! idxLP: number of the lagrange point of interest
         ! idxComp: index of the component
         ! method: method to use to get the value (0: bisectioning, >0: approximations)
         ! n: number of supporting points on each positive axis
         integer, intent(in) :: idxLP, idxComp, method, n
         ! varibales:
         real(dp) :: m1, m2, d_L1_comp, d_L2_CoM, m, phi_at_LP, x, y, z, radius, offset, radius_cell
         integer :: ix, iy, iz, ix_min, ix_max, iy_min, iy_max, iz_min, iz_max
         integer :: c_sphere, c_LP
         
         if (.not. LP_inputs_valid(mgtr, mlow, a, idxLP, method)) then
            get_LP_volume = -1.0d0
            return
         else if (idxLP .gt. 2) then
            write(*,*) "Error in get_LP_volume: only the first 2 lagrange points are supported: idxLP=", idxLP
            get_LP_volume = -1.0d0
            return
         else if ((idxComp .gt. 2).or.(idxComp .lt. 0)) then
            write(*,*) "Error in get_LP_volume: only two stars in a binary: idxComp=", idxComp
            get_LP_volume = -1.0d0
            return
         endif
         
         if (mlow .gt. mgtr) then
            write(*,*) "Warning: mlow=", mlow, ">mgtr=", mgtr, " => swap masses"
            m1 = mlow
            m2 = mgtr
         else
            m1 = mgtr
            m2 = mlow
         endif
         ! get total mass
         m = m1/(m1+m2)
         
         select case (idxLP)
         case (1)
            ! L1
            d_L1_comp = get_distance_CoM_LP(m1, m2, 1.0d0, idxLP, method)
            phi_at_LP = phi_a_over_G_M(d_L1_comp, 0.0d0, 0.0d0, m)
            select case (idxComp)
            case (0)
               d_L1_comp = get_distance_CoM_LP(m1, m2, 1.0d0, idxLP, method)
               offset = 0.0d0
            case (1)
               d_L1_comp = get_distance_LP_comp(m1, m2, 1.0d0, idxLP, idxComp, method)
               offset = m-1.0d0
            case (2)
               d_L1_comp = get_distance_LP_comp(m1, m2, 1.0d0, idxLP, idxComp, method)
               offset = m
            case default
               write(*,*) "Error in get_LP_volume: unknown case of idxComp=", idxComp
               get_LP_volume = -1.0d0
               return
            end select
            radius = d_L1_comp
         case (2)
            ! L2
            d_L2_CoM = get_distance_CoM_LP(m1, m2, 1.0d0, idxLP, method)
            offset = 0.0d0
            phi_at_LP = phi_a_over_G_M(d_L2_CoM, 0.0d0, 0.0d0, m)
            radius = d_L2_CoM
         case default
            write(*,*) "Error in get_LP_volume: unknown case of idxLP=", idxLP
            get_LP_volume = -1.0d0
            return
         end select

         if (radius .le. 0) then
            get_LP_volume = radius
            return
         endif

         c_sphere = 0
         c_LP = 0
         ix_max = n
         ix_min = -ix_max
         radius_cell = radius/(real(n,8)+0.5d0)
         
         do ix=ix_min,ix_max
            iy_max = int(sqrt(real(ix_max*ix_max-ix*ix,8)))
            iy_min = -iy_max
            do iy=iy_min,iy_max
               iz_max = int(sqrt(real(ix_max*ix_max-ix*ix-iy*iy,8)))
               iz_min = -iz_max
               do iz=iz_min,iz_max
                  c_sphere = c_sphere+1
                  x = radius_cell*real(ix,8) + offset
                  y = radius_cell*real(iy,8)
                  z = radius_cell*real(iz,8)
                  if (phi_a_over_G_M(x, y, z, m) .le. phi_at_LP) then
                     c_LP = c_LP+1
                     if (idxLP==1) then
                        select case (idxComp)
                        case (1)
                           if (m/sqrt((x+1.0d0-m)*(x+1.0d0-m)+y*y+z*z) .lt. (1.0d0-m)/sqrt((x-m)*(x-m)+y*y+z*z)) then
                              c_LP = c_LP-1
                           endif
                        case (2)
                           if (m/sqrt((x+1.0d0-m)*(x+1.0d0-m)+y*y+z*z) .gt. (1.0d0-m)/sqrt((x-m)*(x-m)+y*y+z*z)) then
                              c_LP = c_LP-1
                           endif
                        end select
                     endif
                  endif
               enddo
            enddo
         enddo
         ! result fraction of sphere
         get_LP_volume = real(c_LP,8)*pow(radius_cell*a,3.0d0)
         
         contains
         
         real(dp) function phi_a_over_G_M(x, y, z, m_1)
            ! unit free potential from unit free coordinates and mass, where center of mass is at 0,0,0 and two masses are on x axis
            real(dp), intent(in) :: x, y, z, m_1
            real(dp) :: x2, y2, z2, m_2, x_1, x_2, x2_1, x2_2
            if ((m .gt. 1.0d0).or.(m .lt. 0.0d0)) then
               ! m not in range
               write (*,*) "Error: m=", m, " out of range [0.0,1.0]"
            endif
            ! get auxiliary variables
            x2 = x*x
            y2 = y*y
            z2 = z*z
            m_2 = 1.0d0-m_1
            x_1 = (x+m_2)
            x2_1 = x_1*x_1
            x_2 = (x-m_1)
            x2_2 = x_2*x_2
            ! get potential phi(x,y,z)*a/(G*M) = -0.5*([x/a]^2+[y/a]^2) - M1/M/sqrt([x/a+M2/M]^2+[y/a]^2+[z/a]^2) - M2/M/sqrt([x/a-M1/M]^2+[y/a]^2+[z/a]^2)
            phi_a_over_G_M = -0.5*(x2+y2) - m_1/sqrt(x2_1+y2+z2) - m_2/sqrt(x2_2+y2+z2)
            
         end function phi_a_over_G_M
         
      end function get_LP_volume


! *********************************************************************************************************************************


      end module CE_lagrange_points
