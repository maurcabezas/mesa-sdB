 
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

   module binary_rlobe
   
   use star_lib
   use star_def
   use crlibm_lib
   use const_def
   use binary_def
   
   implicit none

   contains
   
   subroutine eval_rlobe(b, s_i, nu)
      integer, intent(in) ::  s_i 
      real(dp), intent(in) :: nu
      type (binary_info), pointer :: b
      real(dp) :: rlobe, q
      
!       ierr = 0
!       call binary_ptr(binary_id, b, ierr)
!       if (ierr /= 0) then
!          write(*,*) 'failed in binary_ptr'
!          return
!       end if
      
      if ( b% rlobe_scheme == "alexey" ) then
      
         if (s_i == 1) then
            q = b% m( 1 )/ b% m( 2 ) 
         else
            q = b% m( 2 )/ b% m( 1 ) 
         end if
      
         rlobe = b% separation * eval_frl1_alexey(q)
      
      else if ( b% rlobe_scheme == "eggleton" ) then
         rlobe = eval_rlobe_Eggleton(b, s_i)
         
      else if ( b% rlobe_scheme == "sepinsky" ) then
         rlobe = eval_rlobe_Sepinsky(b, s_i, nu)
         
      else
         write(*,*) 'Roche lobe scheme not recognized'
         rlobe = 0
         return
         
      end if
      
      ! store the RL size in the binary pointer
      b% rl(s_i) = rlobe
      
   end subroutine eval_rlobe
   
   subroutine eval_rlobe_L3(b, s_i, nu)
      integer, intent(in) ::  s_i 
      real(dp), intent(in) :: nu
      type (binary_info), pointer :: b
      real(dp) :: rlobe, q
      
      if (s_i == 1) then
         q = b% m( 1 )/ b% m( 2 ) 
      else
         q = b% m( 2 )/ b% m( 1 ) 
      end if
      
      rlobe = b% separation * eval_frl3_alexey(q)
      
      ! store the RL size in the binary pointer
      b% rl3(s_i) = rlobe
      
   end subroutine eval_rlobe_L3
   
!============================================================================================
! Location of the L1 and L3 points
!============================================================================================
   
   real(dp) function eval_xrl1_alexey(q) result(xrl1)
      real(dp) :: q, lq ! mass ratio and log10 mass ratio
      ! evaluate the location of the L1 point
         
      lq = safe_log10_cr( q )
      
      xrl1 = -0.301014 + 0.205896*lq - 0.0488536*lq**2 - 0.00395175*lq**3 &
               + 0.00386047*lq**4 + 0.000406739*lq**5 - 0.000604208*lq**6 &
               + 0.00010853*lq**7
      
      xrl1 = exp10_cr(xrl1)
      
   end function eval_xrl1_alexey
   
   
   real(dp) function eval_xrl3_alexey(q) result(xrl3)
      real(dp) :: q, lq ! mass ratio and log10 mass ratio
      ! evaluate the location of the L3 point
         
      lq = safe_log10_cr( q )
      
      xrl3 = - 0.155869 + 0.233746*lq - 0.103023*lq**2 - 0.0226828*lq**3 &
             + 0.0232574*lq**4 + 0.00863963*lq**5 - 0.00807559*lq**6 &
             - 0.00128099*lq**7 + 0.00201497*lq**8 - 0.000398868*lq**9
      
      xrl3 = - exp10_cr(xrl3)
      
   end function eval_xrl3_alexey

!============================================================================================
! Radii of the volume equivalend spheres corresponding to RL1 and RL3 in units of separation
!============================================================================================
   
   real(dp) function eval_frl1_alexey(q) result(frl1)
      real(dp) :: q, lq ! mass ratio and log10 mass ratio
         
      lq = safe_log10_cr( q )
      
      frl1 = -0.420297 + 0.232069*lq - 0.0438153*lq**2 - 0.00567072*lq**3 &
               + 0.00870908*lq**4 - 0.0205853*lq**5 - 0.0169055*lq**6 + 0.0876934*lq**7 &
               - 0.0227041*lq**8 - 0.13918*lq**9 + 0.118513*lq**10 + 0.0627232*lq**11 &
               - 0.122257*lq**12 + 0.0345071*lq**13 + 0.0297013*lq**14 - 0.0253245*lq**15 &
               + 0.00734239*lq**16 - 0.000780009*lq**17
      
      frl1 = exp10_cr(frl1)
      
   end function eval_frl1_alexey
   
   
   real(dp) function eval_frl3_alexey(q) result(frl3)
      real(dp) :: q, lq ! mass ratio and log10 mass ratio
         
      lq = safe_log10_cr( q )
      
      frl3 = -0.300573 + 0.254168*lq - 0.0805788*lq**2 - 0.0288495*lq**3 &
               + 0.0340875*lq**4 + 0.06476*lq**5 - 0.317589*lq**6 - 0.403948*lq**7 &
               + 1.99877*lq**8 + 1.49487*lq**9 - 6.78449*lq**10 -3.00681*lq**11 &
               + 12.9994*lq**12 + 3.29047*lq**13 - 14.0427*lq**14 - 1.845*lq**15 &
               + 7.97442*lq**16 + 0.416024*lq**17 - 1.85215*lq**18
      
      frl3 = exp10_cr(frl3)
      
   end function eval_frl3_alexey
   
   
!============================================================================================
! Radii of the volume equivalend spheres corresponding to RL1 and RL3
!============================================================================================

   
   real(dp) function eval_rlobe_Eggleton(b, s_i) result(rlobe)
      integer :: s_i ! s_i is index of the star for which to calc the roche lobe
      type (binary_info), pointer :: b
      real(dp) :: q
      ! Roche lobe size for star according to
      ! the approximation of Eggleton 1983, apj 268:368-369
      
      if (s_i == 1) then
         q = pow_cr(b% m( 1 )/ b% m( 2 ),one_third)
      else
         q = pow_cr(b% m( 2 )/ b% m( 1 ),one_third)
      end if
      
      rlobe = b% separation*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
   
   end function eval_rlobe_Eggleton
   
   
   real(dp) function A_sepinsky(b, s_i, nu) result(A)
      integer :: s_i
      type (binary_info), pointer :: b
      type (star_info), pointer :: s
      real(dp) :: f, omega_p, nu
      
      ! select the correct star
      if (s_i == 1) then
         s => b% s1
      else
         s => b% s2
      end if
      
      omega_p = 2 * pi / b% period * sqrt( 1 + b% eccentricity ) / pow_cr(1 - b% eccentricity, two_thirds)
      
      f = abs( s% omega_avg_surf ) / omega_p
      
      A = f**2 * pow_cr(1 + b% eccentricity, 4d0) / pow_cr( 1 + b% eccentricity * cos( nu ), 3d0 )
   
   end function A_sepinsky
   
   
   real(dp) function eval_rlobe_Sepinsky(b, s_i, nu) result(rlobe)
      integer :: s_i ! s_i is index of the star for which to calc the roche lobe
      type (binary_info), pointer :: b
      real(dp) :: nu, logq, A, logA, Rl_egg, factor
      real(dp) :: g0, g1, g2, g3
      ! Roche lobe size for star according to
      ! the approximation of Sepinsky et al. 2007 ApJ 660:1624S
      ! Taking into account the rotation of the star and the eccentricity of the orbit
      
      
      if (s_i == 1) then
         logq = log( b% m( 1 )/ b% m( 2 ) )
      else
         logq = log( b% m( 2 )/ b% m( 1 ) )
      end if
      
      A = A_sepinsky(b, s_i, nu)
      logA = log( A )
      
      Rl_egg = eval_rlobe_Eggleton(b, s_i)
      
      
      if ( logq >= 0 .and. logA <= -1 ) then
         factor =  1.226d0 - 0.21d0 * A - 0.15d0 * (1d0 - A) * exp( ( 0.25d0*A - 0.3d0 ) * pow_cr( logq, 1.55d0 ) )
         
      else if ( logq < 0 .and. logA <= -0.1 ) then
         factor = 1 + 0.11d0 * ( 1d0 - A ) - 0.05d0 * ( 1d0 - A ) * exp( -( 0.5d0 * (1d0 + A) + logq )**2d0 )
         
      else if ( logq < 0 .and. logA <= 0.2 .and. logA > -0.1 ) then
         g0 = 0.9978d0 - 0.1229d0 * logA - 0.1273d0 * logA**2d0
         g1 = 0.001d0 + 0.02556d0 * logA
         g2 = 0.0004d0 + 0.0021d0 * logA
         factor = g0 + g1 * logq + g2 * logq**2 
         
      else if (logq >= 0 .and. logA <= 0.2 .and. logA > -0.1 ) then
         g0 = 1.0071d0 - 0.0907d0 * logA - 0.0495d0 * logA**2d0
         g1 = -0.004d0 - 0.163d0 * logA -0.214d0 * logA**2d0
         g2 = 0.00022d0 - 0.0108d0 * logA - 0.02718d0 * logA**2d0
         factor = g0 + g1 * logq + g2 * logq**2
         
      else if (logq < 0 .and. logA > 0.2 ) then
         g0 = 6.3014d0 * pow_cr( logA, 1.3643d0 ) / ( exp(2.3644d0 * pow_cr( logA, 0.70748d0) ) &
         - 1.4413 * exp(-0.0000184d0 * pow_cr( logA, -4.5693d0) ) )
         g1 = logA / ( 0.0015d0 * exp(8.84d0 * pow_cr( logA, 0.282d0 ) ) + 15.78d0 )
         g2 = ( 1d0 + 0.036d0 * exp( 8.01d0 * pow_cr( logA, 0.879d0 ) ) ) / ( 0.105d0 * exp(7.91d0 * pow_cr( logA, 0.879d0 ) ) )
         g3 = 0.991d0 / (1.38d0 * exp( -0.035d0 * pow_cr( logA, 0.76d0 ) ) + 23.0d0 * exp( -2.89d0 * pow_cr( logA, 0.76d0 ) ) )
         factor = g0 + g1 * exp( -g2 * ( logq + g3 )**2d0 )
         
      else
         g0 = 1.895d0 * pow_cr( logA, 0.837d0 ) / ( exp( 1.636d0 * pow_cr( logA, 0.789d0 ) ) - 1d0 )
         g1 = 4.3d0 * pow_cr( logA, 0.98d0 ) / ( exp( 2.5d0 * pow_cr( logA, 0.66d0 ) ) + 4.7d0 )
         g2 = 1d0 / ( 8.8d0 * exp( -2.95d0 * pow_cr( logA, 0.76d0 ) ) + 1.64d0 * exp( -0.03d0 * pow_cr( logA, 0.76d0 ) ) )
         g3 = 0.256d0 * exp( -1.33d0 * pow_cr( logA, 2.9d0 ) ) * ( 5.5d0 * exp( 1.33d0 * pow_cr( logA, 2.9d0 ) ) + 1d0 )
         factor =  g0 + g1 * exp( -g2 * pow_cr( logq, g3 ) )
      
      end if
         
      rlobe = factor * Rl_egg
      
   
   end function eval_rlobe_Sepinsky
   
   
   
   end module binary_rlobe
   
