! ***********************************************************************
!
!   This file is part of a mesa extension.
!   Authors of this file: Matthias U. Kruckow
!
! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
 
      module CE_d_mix

      ! consult star/other/README for general usage instructions
      ! control name: use_other_D_mix = .true.
      ! procedure pointer: s% other_D_mix => my_routine


      use star_def
      use const_def
      use CE_orbit, only: AtoP, TukeyWindow, calc_quantities_at_comp_position

      implicit none
      
      contains
      
      
      subroutine CE_enhance_D_mix(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_companion_position
         real(dp) :: CE_n_acc_radii
         real(dp) :: a_tukey, enhanced_D_mix, maxenhancement_factor, ff
         real(dp) :: R_acc, R_acc_low, R_acc_high
         integer :: k
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_companion_position = s% xtra(2)

         ! If the star is in the initial relaxation phase, skip additional mixing
         if (s% doing_relax) return
         ! If companion is outside star, skip additional mixing
         if (CE_companion_position*Rsun > s% r(1)) return

         ! If system merged, skip additional mixing
         if (s% lxtra(1)) then
            return
         endif
         
         CE_n_acc_radii = s% xtra(5)

         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra(12)
         R_acc_low = s% xtra(13)
         R_acc_high = s% xtra(14)
         
         ! factor from inlist
         maxenhancement_factor = s% x_ctrl(22)

         ! Tukey window scale
         a_tukey = 0.5d0

         ! Loop through the layers inside the accretion radius
         do k = 1, s% nz
            ! The mixing should only apply to the envelope of the star and not in the core
            ! When we reach the core boundary we exit the loop
            if (s% m(k) < s% he_core_mass * Msun) exit
            ! Get inner or outer accretion radius
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if
            ! Mix through out the accretion radius in the current time step
            enhanced_D_mix = maxenhancement_factor*0.5d0*R_acc*R_acc/s% dt
            ! Use TukeyWindow to determine where D_mix should be between the usual on and the enhanced_D_mix
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0d0 * R_acc), a_tukey)
            s% D_mix(k) = s% D_mix(k) * (1.0d0 - ff) + enhanced_D_mix * ff
         end do
         
      end subroutine CE_enhance_D_mix
      

      end module CE_d_mix
      
      
