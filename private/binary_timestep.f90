! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton and Pablo Marchant
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


      module binary_timestep

      use const_def
      use crlibm_lib
      use star_lib
      use star_def
      use binary_def

      implicit none

      contains

      subroutine set_star_timesteps(b) ! sets the smallest next timestep for all stars
         type (binary_info), pointer :: b
         integer :: i, l
         real(dp) :: dt_min, rel_overlap
         type (star_info), pointer :: s
         integer :: ierr
         ierr = 0
         dt_min = 1d99
         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            if (s% dt_next < dt_min) then
               dt_min = s% dt_next
            end if
         end do
         if (b% max_timestep <= dt_min) then
            dt_min = b% max_timestep
         else
            b% dt_why_reason = b_Tlim_comp
         end if
         ! just to be sure we dont cause a segfault
         if (b% dt_why_reason < 1 .or. b% dt_why_reason > b_numTlim) then
            dt_why_str(Tlim_binary) = " "
         else
            dt_why_str(Tlim_binary) = binary_dt_why_str(b% dt_why_reason)
         end if
         do l = 1, num_stars
            if (l == 1 .and. b% point_mass_i == 1) then
               i = 2
            else
               i = l
            end if
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            if (s% dt_next > dt_min) then
               s% dt_next = dt_min
               s% why_Tlim = Tlim_binary
            end if
         end do
         
      end subroutine set_star_timesteps

      integer function binary_pick_next_timestep(b)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         
         real(dp) :: &
            env_change, dtm, dtj, dta, dtr, dte, dtdm, &
            j_change, sep_change, rel_gap_change, e_change, set_dt

         include 'formats.inc'

         dtm = 1d99
         dtj = 1d99
         dta = 1d99
         dtr = 1d99
         dte = 1d99
         dtdm = 1d99

         binary_pick_next_timestep = keep_going

         s => b% s_donor

         if (b% max_timestep < 0) b% max_timestep = b% s_donor% dt

         b% env = s% star_mass - s% he_core_mass 
         if (b% env_old /= 0) then
            env_change = b% env - b% env_old
         else
            env_change = 0
         end if
         
         if (b% rl_relative_gap_old(b% d_i) /= 0) then
            rel_gap_change = b% rl_relative_gap_old(b% d_i) - b% rl_relative_gap(b% d_i)
         else
            rel_gap_change = 0
         end if
         
         if (b% angular_momentum_j_old /= 0) then
            j_change = b% angular_momentum_j - b% angular_momentum_j_old
         else
            j_change = 0
         end if
         
         if (b% separation_old /= 0) then
            sep_change = b% separation - b% separation_old
         else
            sep_change = 0
         end if
         if (b% eccentricity_old /= 0) then
             e_change = b% eccentricity - b% eccentricity_old
         else
             e_change = 0
         end if
   
         ! get limits for dt based on relative changes
         if (b% fm > 0) then
            dtm = s% time_step/(abs(env_change/max(b% env, b% fm_limit))/b% fm+1d-99)
         end if
         
         if (b% fr > 0) then
            dtr = s% time_step/ &
                (abs(rel_gap_change/max(-b% rl_relative_gap(b% d_i), b% fr_limit))/b% fr+1d-99)
            if (b% rl_relative_gap_old(b% a_i) /= 0) then
               rel_gap_change = b% rl_relative_gap_old(b% a_i) - b% rl_relative_gap(b% a_i)
            else
               rel_gap_change = 0
            end if
            dtr = min(dtr, s% time_step/ &
                (abs(rel_gap_change/max(-b% rl_relative_gap(b% a_i), b% fr_limit))/b% fr+1d-99))
         end if
         if (dtr < b% fr_dt_limit) dtr = b% fr_dt_limit

         if (b% fj > 0) then
            dtj = s% time_step/(abs(j_change/b% angular_momentum_j)/b% fj+1d-99)
         end if

         if (b% fa > 0) then
            dta = s% time_step/(abs(sep_change/b% separation)/b% fa+1d-99)
         end if

         if (b% fe > 0) then
            dte = s% time_step/(abs(e_change/ max( b% eccentricity, b% fe_limit ))/b% fe+1d-99)
         end if

         if (abs(b% s_donor% mstar_dot) > 0) then
            dtdm = b% fdm * b% s_donor% mstar / abs(b% s_donor% star_mdot) / secyer
         end if

         set_dt = min(dtm, dtr, dtj, dta, dte, dtdm)
         
         if (set_dt == dtm) then
            b% dt_why_reason = b_Tlim_env
         else if (set_dt == dtr) then
            b% dt_why_reason = b_Tlim_roche
         else if (set_dt == dtj) then
            b% dt_why_reason = b_Tlim_jorb
         else if (set_dt == dta) then
            b% dt_why_reason = b_Tlim_sep
         else if (set_dt == dte) then
            b% dt_why_reason = b_Tlim_ecc
         else if (set_dt == dtdm) then
            b% dt_why_reason = b_Tlim_dm
         else
            stop 'Something wrong in binary timestep'
         end if

         if (set_dt < 1d-7) set_dt = 1d-7 ! there's a limit to everything

         b% max_timestep = exp10_cr(b% dt_softening_factor*log10_cr(set_dt*secyer) + &
             (1-b% dt_softening_factor)*log10_cr(b% max_timestep))

         ! use variable varcontrols for different phases of evolution
         if (abs(b% mtransfer_rate)/Msun*secyer > 1d-20) then
            if (b% s_donor% center_h1 > 1d-12 .and. b% varcontrol_case_a > 0d0) then
               b% s_donor% varcontrol_target = b% varcontrol_case_a
               if (b% point_mass_i == 0) &
                   b% s_accretor% varcontrol_target = b% varcontrol_case_a
            else if (b% s_donor% center_h1 < 1d-12 .and. b% varcontrol_case_b > 0d0) then
               b% s_donor% varcontrol_target = b% varcontrol_case_b
               if (b% point_mass_i == 0) &
                   b% s_accretor% varcontrol_target = b% varcontrol_case_b
            end if
         else
            if (b% s_donor% center_h1 > 1d-12) then
               if (b% varcontrol_ms > 0d0) &
                   b% s_donor% varcontrol_target = b% varcontrol_ms
            else
               if (b% varcontrol_post_ms > 0d0) &
                   b% s_donor% varcontrol_target = b% varcontrol_post_ms
            end if

            if (b% point_mass_i == 0) then
               if (b% s_accretor% center_h1 > 1d-12) then
                  if (b% varcontrol_ms > 0d0) &
                      b% s_accretor% varcontrol_target = b% varcontrol_ms
               else
                  if (b% varcontrol_post_ms > 0d0) &
                      b% s_accretor% varcontrol_target = b% varcontrol_post_ms
               end if
            end if
         end if
         
      end function binary_pick_next_timestep
      

      end module binary_timestep
