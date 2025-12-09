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
         logical :: CE_accretion_scheme
         real(dp) :: f, w, log_mdot_out, mdot_out, SMS_mdot

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_mdot_factor_increase = s% x_ctrl(11)
         CE_mdot_factor_decrease = s% x_ctrl(12)
         CE_mdot_smooth_limit = s% x_ctrl(13)
         CE_mdot_max = s% x_ctrl(14)

         CE_accretion_scheme = s% x_logical_ctrl(9)

         !CE_mdot = s% xtra(7)

         

         !CE wind
         if (-s% mstar_dot_old < CE_mdot_smooth_limit * Msun/secyer .and. -s% xtra(7) > CE_mdot_smooth_limit * Msun/secyer) then
            CE_mdot = -1.01*CE_mdot_smooth_limit * Msun/secyer
         else if (-s% mstar_dot_old > CE_mdot_smooth_limit * Msun/secyer .and. -s% xtra(7) > CE_mdot_smooth_limit * Msun/secyer) then
            if (-s% xtra(7) > -CE_mdot_factor_increase * s% mstar_dot_old  ) then
               CE_mdot = CE_mdot_factor_increase * s% mstar_dot_old
            else if (-s% xtra(7) < -1./CE_mdot_factor_decrease * s% mstar_dot_old  ) then
               CE_mdot = 1./CE_mdot_factor_decrease* s% mstar_dot_old
            else
               CE_mdot = s% xtra(7) !/ s% dt
            endif
         else
            CE_mdot = s% xtra(7) !/ s% dt
         endif

         if (CE_mdot < -CE_mdot_max * Msun/secyer) CE_mdot = -CE_mdot_max* Msun/secyer
         if (CE_mdot/(Msun/secyer) < -1d-20) write(*,*) "** CEmdot ** ", CE_mdot/(Msun/secyer)

         s% mstar_dot = s% mstar_dot + CE_mdot !In gr/s

         if (s%x_logical_ctrl(3) .and. (.not. s% doing_relax)) then
            s% Dutch_scaling_factor = s% xtra(21) ** s% x_ctrl(16)
            write(*,*) "**Pulsational Winds** ", s% xtra(21), s% Dutch_scaling_factor
         endif

         ! Accretion scheme 
         ! Subroutine for accretion of mass based in Haemmerle et al. 2016
         ! used to build extremly massive and supermassive stars in Ramirez-Galeano et al.2025
         SMS_mdot = 0.d0
         if (CE_accretion_scheme) then 
            !s% mstar_dot = 0.0d0
            w = 0.0d0
            !s% explicit_mstar_dot = s% mstar_dot
            ! Mass to reach 

            if (CE_accretion_scheme) then
              if (s% star_mass <= 5.0d0) then
              
                f = 1.0/3.0
                
              else if (s% star_mass > 5.0) then
              
                f = 1.0/11.0
              end if

              log_mdot_out = -5.28d0 + s% log_surface_luminosity *(0.752d0 - 0.0278d0*s% log_surface_luminosity) ![M_sun/yr]
              mdot_out = 10.d0**(log_mdot_out) * (Msun/secyer) ![gr/s]  10 for high accretion
    
              w = f/(1-f)*mdot_out 
                 
            endif
            

            s% mstar_dot = (s% mstar_dot + w)
            s% explicit_mstar_dot = s% mstar_dot 
            SMS_mdot = w

         endif
         
         
      end subroutine CE_other_adjust_mdot

      end module CE_adjust_mdot
