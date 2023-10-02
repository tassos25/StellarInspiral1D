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

      module CE_adjust_mdot

      use star_def
      use const_def

      implicit none

      contains

! ***********************************************************************
      subroutine CE_other_adjust_mdot(id, ierr)
      ! set use_other_adjust_mdot = .true. to enable this.
      ! your routine will be called after winds and before mass adjustment
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: CE_mdot, CE_mdot_limit, CE_mdot_factor_increase, CE_mdot_factor_decrease
         real(dp) :: CE_mdot_smooth_limit, CE_mdot_max

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_mdot_factor_increase = s% x_ctrl(11)
         CE_mdot_factor_decrease = s% x_ctrl(12)
         CE_mdot_smooth_limit = s% x_ctrl(13)
         CE_mdot_max = s% x_ctrl(14)

         !CE_mdot = s% xtra(7)

         if (-s% mstar_dot_old < CE_mdot_smooth_limit * Msun/secyer .and. -s% xtra(7) > CE_mdot_smooth_limit * Msun/secyer) then
            CE_mdot = -1.01*CE_mdot_smooth_limit * Msun/secyer
         else if (-s% mstar_dot_old > CE_mdot_smooth_limit * Msun/secyer .and. -s% xtra(7) > CE_mdot_smooth_limit * Msun/secyer) then
            if (-s% xtra(7) > -CE_mdot_factor_increase * s% mstar_dot_old  ) then
               CE_mdot = CE_mdot_factor_increase * s% mstar_dot_old
            else if (-s% xtra(7) < -1./CE_mdot_factor_decrease * s% mstar_dot_old  ) then
               CE_mdot = 1./CE_mdot_factor_decrease* s% mstar_dot_old
            else
               CE_mdot = s% xtra(7)
            endif
         else
            CE_mdot = s% xtra(7)
         endif

         if (CE_mdot < -CE_mdot_max * Msun/secyer) CE_mdot = -CE_mdot_max* Msun/secyer
         if (CE_mdot/(Msun/secyer) < -1d-20) write(*,*) "** CEmdot ** ", CE_mdot/(Msun/secyer)

         s% mstar_dot = s% mstar_dot + CE_mdot !In gr/s

         if (s%x_logical_ctrl(3) .and. (.not. s% doing_relax)) then
            s% Dutch_scaling_factor = s% xtra(21) ** s% x_ctrl(16)
            write(*,*) "**Pulsational Winds** ", s% xtra(21), s% Dutch_scaling_factor
         endif
      end subroutine CE_other_adjust_mdot


      end module CE_adjust_mdot
