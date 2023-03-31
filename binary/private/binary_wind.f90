 
! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton and Joris Vos
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

   module binary_wind
   
   use star_lib
   use star_def
   use crlibm_lib
   use binary_def
   
   implicit none

   contains
   
   subroutine eval_wind_xfer_fractions(b)
      type (binary_info), pointer :: b
      integer :: ierr
      ierr = 0
      
      ! for the primary
      if (b% point_mass_i /= 1) then
         if (.not. b% do_wind_mass_transfer_1) then
            b% wind_xfer_fraction(1) = 0d0
         else if(.not. b% use_other_binary_wind_transfer) then
            call Bondy_Hoyle_wind_transfer(b% binary_id, 1, ierr)
         else
            call b% other_binary_wind_transfer(b% binary_id, 1, ierr)
         end if
      end if
      
      ! check if secondary needs wind transfer
      if (b% point_mass_i /= 2) then
         
         if (.not. b% do_wind_mass_transfer_2) then
            b% wind_xfer_fraction(2) = 0d0
         else if(.not. b% use_other_binary_wind_transfer) then
            call Bondy_Hoyle_wind_transfer(b% binary_id, 2, ierr)
         else
            call b% other_binary_wind_transfer(b% binary_id, 2, ierr)
         end if
      end if
         
   end subroutine eval_wind_xfer_fractions
   
   subroutine Bondy_Hoyle_wind_transfer(binary_id, s_i, ierr)
      integer, intent(in) :: binary_id, s_i ! s_i is index of the wind mass losing star
      integer, intent(out) :: ierr

      ! wind transfer fraction based on Bondy-Hoyle mechanism as described in
      ! Hurley et al. 2002, MNRAS, 329, 897-928

      type(binary_info), pointer :: b
      type (star_info), pointer :: s
      real(dp) :: v_orb, v_wind, b_BH
      real(dp) :: alpha  ! Bondy-Hoyle alpha for that star
      real(dp) :: max_xfer  ! Maximum transfer fraction

      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in binary_ptr'
         return
      end if
      
      if (s_i == 1) then
         s => b% s1
         alpha = b% wind_BH_alpha_1
         max_xfer = b% max_wind_transfer_fraction_1
      else
         s => b% s2
         alpha = b% wind_BH_alpha_2
         max_xfer = b% max_wind_transfer_fraction_2
      end if
      
      ! orbital speed Hurley et al 2002 eq. 8
      v_orb = sqrt(s% cgrav(1) * b% m(s_i) / b% separation) !cm/s
      
      ! windspeed from Hurley et al 2002 eq. 9
      v_wind = sqrt( 2d0 / 8d0 *  s% cgrav(1) * b% m(s_i) / b% r(s_i) )
      
      ! Bondy-Hoyle transfer fraction Hurley et al. 2002 eq. 6
      b% wind_xfer_fraction(s_i) = alpha / b% separation**2d0 /&
                  (2d0 * sqrt(1d0 - b% eccentricity**2d0)) *&
                  (s% cgrav(1) * b% m(3-s_i) / v_wind**2d0)**2d0 *&
                  pow_cr(1d0 + (v_orb/v_wind)**2d0,-1.5d0)
                  
      ! limit to provided maximum
      b% wind_xfer_fraction(s_i) = min(max_xfer, b% wind_xfer_fraction(s_i))
      
   end subroutine Bondy_Hoyle_wind_transfer
   
   subroutine Tout_enhance_wind(b, s)
      type (binary_info), pointer :: b
      type (star_info), pointer :: s

      ! Tidaly enhance wind mass loss as described by
      ! Tout & Eggleton 1988,MNRAS,231,823 (eq. 2)
      real(dp) :: B_wind  ! enhancement parameter, B in eq. 2
      integer :: i, s_i
      real(dp) :: dm
      real(dp), DIMENSION(b% anomaly_steps):: rl_d, r_rl, mdot

      if (s% id == b% s1% id) then
         if (.not. b% do_enhance_wind_1) return
         B_wind = b% tout_B_wind_1
         s_i = 1
      else
         if (.not. b% do_enhance_wind_2) return
         B_wind = b% tout_B_wind_2
         s_i = 2
      end if
      
      ! phase dependent roche lobe radius
      rl_d = (1-b%eccentricity**2) / (1+b%eccentricity*cos(b% theta_co)) * b% rl(s_i)
      do i = 1,b% anomaly_steps !limit radius / roche lobe
         r_rl(i) = min(pow6(b% r(s_i) / rl_d(i)), pow6(0.5d0))
      end do
      
      ! actual enhancement
      mdot = s% mstar_dot * (1 + B_wind * r_rl)
      
      dm = 0d0
      do i = 2,b% anomaly_steps ! trapezoidal integration
         dm = dm + 0.5 * (mdot(i-1) + mdot(i)) * (b% time_co(i) - b% time_co(i-1)) 
      end do
      
      ! remember mass-loss is negative!
      !b% mdot_wind_theta = b% mdot_wind_theta + mdot ! store theta dependance for edot
      s% mstar_dot = dm ! return enhanced wind mass loss
      
   end subroutine Tout_enhance_wind
   
   end module binary_wind
   
