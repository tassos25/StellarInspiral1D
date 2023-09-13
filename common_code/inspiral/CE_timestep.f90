! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module CE_timestep

      use const_def
      use star_def

      implicit none


      contains


      integer function CE_pick_next_timestep(s)
         use const_def, only: secyer
         type (star_info), pointer :: s


         CE_pick_next_timestep = keep_going

         !Do the additional timestep controls only before the merger.
         if (.not. s% lxtra1) then
            CE_pick_next_timestep = worst_result(CE_pick_next_timestep, CE_check_energy(s))
            CE_pick_next_timestep = worst_result(CE_pick_next_timestep, CE_check_separation(s))
            CE_pick_next_timestep = worst_result(CE_pick_next_timestep, CE_check_ang_mom(s))
         endif


      end function CE_pick_next_timestep



      integer function CE_check_energy(s)
         type (star_info), pointer :: s
         real(dp) :: dE_fraction, dE_limit


         ! Fractional energy change
         dE_fraction = abs((s% xtra1 * s% dt) / s% xtra8)

         ! Limit from inlist
         dE_limit = s% x_ctrl(8)

         ! Calculate next timestep
         if ((dE_fraction > dE_limit)) then
            CE_check_energy = retry
         else
            CE_check_energy = keep_going
         end if


      end function CE_check_energy




      integer function CE_check_separation(s)
         type (star_info), pointer :: s
         real(dp) :: dA_fraction, dA_limit


         ! Fractional separation change
         dA_fraction = abs((s% xtra2 - s% xtra2_old) / s% xtra2_old)

         ! Limit from inlist
         dA_limit = s% x_ctrl(9)

         ! Calculate next timestep
         if (dA_fraction > dA_limit) then
            CE_check_separation = retry
         else
            CE_check_separation = keep_going
         end if

      end function CE_check_separation




      integer function CE_check_ang_mom(s)
         type (star_info), pointer :: s
         real(dp) :: dJ_fraction, dJ_limit

         ! Fractional angular momentum change
         dJ_fraction = abs(s% xtra6 * s% dt / s% xtra10)

         ! Limit from inlist
         dJ_limit = s% x_ctrl(10)

         ! Calculate next timestep
         if ((dJ_fraction > dJ_limit)) then
            CE_check_ang_mom = retry
         else
            CE_check_ang_mom = keep_going
         end if


      end function CE_check_ang_mom


      integer function worst_result(result1, result2)
         integer, intent(in) :: result1, result2

         if(result1 == terminate .or. result2 == terminate) then
            worst_result = terminate
            return
         end if

         if(result1 == backup .or. result2 == backup) then
            worst_result = backup
            return
         end if

         if(result1 == retry .or. result2 == retry) then
            worst_result = retry
            return
         end if

         if(result1 == redo .or. result2 == redo) then
            worst_result = redo
            return
         end if

         worst_result = keep_going
         return

      end function worst_result


      end module CE_timestep
