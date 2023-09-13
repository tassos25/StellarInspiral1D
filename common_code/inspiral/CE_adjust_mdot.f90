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

module CE_adjust_mdot

   use star_def
   use const_def

   implicit none


contains

   ! set use_other_adjust_mdot = .true. to enable this.
   ! your routine will be called after winds and before mass adjustment

   subroutine CE_other_adjust_mdot(id, ierr)

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

      ! CE_mdot = s% xtra7

      if (-s% mstar_dot_old < CE_mdot_smooth_limit * Msun/secyer .and. -s% xtra7 > CE_mdot_smooth_limit * Msun/secyer) then
         CE_mdot = -1.01*CE_mdot_smooth_limit * Msun/secyer
      else if (-s% mstar_dot_old > CE_mdot_smooth_limit * Msun/secyer .and. -s% xtra7 > CE_mdot_smooth_limit * Msun/secyer) then
         if (-s% xtra7 > -CE_mdot_factor_increase * s% mstar_dot_old  ) then
            CE_mdot = CE_mdot_factor_increase * s% mstar_dot_old
         else if (-s% xtra7 < -1./CE_mdot_factor_decrease * s% mstar_dot_old  ) then
            CE_mdot = 1./CE_mdot_factor_decrease* s% mstar_dot_old
         else
            CE_mdot = s% xtra7
         endif
      else
         CE_mdot = s% xtra7
      endif

      if (CE_mdot < -CE_mdot_max * Msun/secyer) CE_mdot = -CE_mdot_max* Msun/secyer
      if (CE_mdot/(Msun/secyer) < -1d-20) write(*,*) "** CEmdot ** ", CE_mdot/(Msun/secyer)

      s% mstar_dot = s% mstar_dot + CE_mdot !In gr/s


      if (s%x_logical_ctrl(3) .and. (.not. s% doing_relax)) then
         s% Dutch_scaling_factor = s% xtra21 ** s% x_ctrl(16)
         write(*,*) "**Pulsational Winds** ", s% xtra21, s% Dutch_scaling_factor
      endif

   end subroutine CE_other_adjust_mdot




end module CE_adjust_mdot
