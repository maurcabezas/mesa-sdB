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


      module binary_evolve

      use const_def
      use crlibm_lib
      use star_lib
      use star_def
      use crlibm_lib
      use binary_def
      use binary_rlobe

      implicit none

      contains

      subroutine binarydata_init(b)
         use utils_lib, only: is_bad
         type (binary_info), pointer :: b
         integer :: finish_step_result
         integer :: i
         real(dp) :: r_isco, Z1, Z2
         include 'formats.inc'

         b% doing_first_model_of_run = .true.
         b% open_new_history_file = .true.

         b% s_donor => b% s1
         b% s_accretor => b% s2
         b% d_i = 1
         b% a_i = 2

         b% max_timestep = 1d99
         b% change_factor = b% initial_change_factor

         if (b% point_mass_i /= 1) then
            initial_mass(1) = b% s1% mstar / Msun
         else
            initial_mass(1) = b% m1
         end if
         if (b% point_mass_i /= 2) then
            initial_mass(2) = b% s2% mstar / Msun
         else
            initial_mass(2) = b% m2
         end if

         b% m(1) = initial_mass(1)*Msun
         b% m(2) = initial_mass(2)*Msun
         if (b% point_mass_i /= 1) then
            b% r(1) = Rsun*b% s1% photosphere_r
         else
            b% r(1) = 0
         end if
         if (b% point_mass_i /= 2) then
            b% r(2) = Rsun*b% s2% photosphere_r
         else
            b% r(2) = 0
         end if
         
         if (b% initial_period_in_days <= 0) then ! calculate from initial_separation_in_Rsuns
            call set_separation_eccentricity(b% binary_id, &
               b% initial_separation_in_Rsuns*Rsun, b% initial_eccentricity)
         else
            call set_period_eccentricity(b% binary_id, &
               b% initial_period_in_days*(24d0*60d0*60d0), b% initial_eccentricity)
         end if

         if (is_bad(b% rl_relative_gap(1))) stop 'binarydata_init'
         if (is_bad(b% rl_relative_gap(2))) stop 'binarydata_init'
         ! these will be adjusted properly by check_radiative_core
         b% have_radiative_core(1) = .true.
         b% have_radiative_core(2) = .true.
         
         ! Set all parameters nessessary for integration over the binary orbit
         ! 1) true anomaly = polar angle from periastron 0 -> 2pi
         do i = 1,b% anomaly_steps 
            b% theta_co(i) = (i-1) * (2 * pi) / b% anomaly_steps
         end do
         
         ! 2) time between periastron and polar angle theta 0 -> 1 (fraction of the
         !    orbital period)
         do i = 1,b% anomaly_steps ! time between periastron and polar angle theta
            b% time_co(i) = ( 2 * atan_cr( sqrt( (1-b% eccentricity)/(1 + b% eccentricity) ) * &
                            tan_cr(b% theta_co(i)/2d0) ) - b% eccentricity * &
                            sqrt(1 - b% eccentricity**2) * sin_cr(b% theta_co(i)) / &
                            (1 + b% eccentricity * cos_cr(b% theta_co(i)) ) ) /2. /pi
            if (i > b% anomaly_steps/2+1) then
               b% time_co(i) = b% time_co(i) + b% time_co(b% anomaly_steps/2+1) * 2
            end if
         end do

         if (b% point_mass_i /= 0) then
            ! this part is only relevant for BH accretors
            if (b% initial_bh_spin < 0d0) then
               b% initial_bh_spin = 0d0
               write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
            else if (b% initial_bh_spin > 1d0) then
               b% initial_bh_spin = 1d0
               write(*,*) "initial_bh_spin is larger than one. It has been set to one."
            end if
            ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
            Z1 = 1d0 + pow_cr(1d0 - b% initial_bh_spin**2,one_third) &
               * (pow_cr(1d0 + b% initial_bh_spin,one_third) + pow_cr(1d0 - b% initial_bh_spin,one_third))
            Z2 = sqrt(3d0*b% initial_bh_spin**2 + Z1**2)
            r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
            ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
            ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
            ! zero spin (r_isco = sqrt(6)).
            b% eq_initial_bh_mass = b% m(b% point_mass_i) * sqrt(r_isco/6d0)
         end if
         
         write(*,*)
         write(*,1) 'm2', b% m2
         write(*,1) 'm1', b% m1
         write(*,1) 'initial_period_in_days', b% initial_period_in_days
         write(*,1) 'initial_separation_in_Rsun', b% separation/Rsun
         write(*,1) 'jdot_multiplier', b% jdot_multiplier
         write(*,1) 'fr', b% fr
         write(*,*)

         b% have_radiative_core = .false.
         b% s_donor% mesh_delta_coeff_pre_ms = 1
         min_binary_period = b% period
         b% min_binary_separation = b% separation
         initial_binary_period = b% period

         b% num_tries = 0

         finish_step_result = binary_finish_step(b)
          
      end subroutine
      
      subroutine set_period_eccentricity(binary_id, period, eccentricity)
         integer, intent(in) :: binary_id
         real(dp) :: period ! in seconds
         real(dp) :: eccentricity
         type (binary_info), pointer :: b
         integer :: ierr
         call binary_ptr(binary_id,b,ierr)
         if (ierr /= 0) return

         b% eccentricity = eccentricity
         b% period = period
         b% separation = &
            pow_cr((b% s_donor% cgrav(1)*(b% m(1)+b% m(2)))*(b% period/(2*pi))**2,one_third)
         call set_angular_momentum_j(binary_id)

      end subroutine set_period_eccentricity
      
      subroutine set_separation_eccentricity(binary_id, separation, eccentricity)
         integer, intent(in) :: binary_id
         real(dp) :: separation ! in cm
         real(dp) :: eccentricity
         type (binary_info), pointer :: b
         integer :: ierr
         call binary_ptr(binary_id,b,ierr)
         if (ierr /= 0) return

         b% eccentricity = eccentricity
         b% separation = separation
         b% period = &
            (2*pi)*sqrt(b% separation*b% separation*b% separation/&
                  (standard_cgrav*(b% m(1)+b% m(2))))
         call set_angular_momentum_j(binary_id)

      end subroutine set_separation_eccentricity
      
      subroutine set_angular_momentum_j(binary_id)
         ! Sets b% angular_momentum_j in terms of the masses, separation and eccentricity
         ! also sets the Roche lobe sizes and relative overflows
         integer, intent(in) :: binary_id
         type (binary_info), pointer :: b
         integer :: ierr
         call binary_ptr(binary_id,b,ierr)
         if (ierr /= 0) return

         b% angular_momentum_j = b% m(1) * b% m(2) * sqrt( b% s_donor% cgrav(1) *&
            b% separation * (1 - b% eccentricity**2) / (b% m(1) + b% m(2)) )

!          b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
!          b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
         call eval_rlobe(b, 1, 0d0)
         call eval_rlobe(b, 2, 0d0)
         call eval_rlobe_L3(b, 1, 0d0)
         call eval_rlobe_L3(b, 2, 0d0)
         
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1) * (1 - b% eccentricity) ) / &
             b% rl(1) / (1 - b% eccentricity) ! gap < 0 means out of contact 
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2) * (1 - b% eccentricity) ) / &
             b% rl(2) / (1 - b% eccentricity) ! gap < 0 means out of contact

      end subroutine set_angular_momentum_j

      subroutine set_donor_star(b)
         type (binary_info), pointer :: b
         logical :: switch_donor
         real(dp) :: mdot_hi_temp
         include 'formats.inc'

         switch_donor = .false.

         if (b% keep_donor_fixed .and. b% mdot_scheme /= "contact") return

         if (b% mdot_scheme == "roche_lobe" .and. &
            abs(b% mtransfer_rate/(Msun/secyer)) < b% mdot_limit_donor_switch .and. &
            b% rl_relative_gap_old(b% a_i) > b% rl_relative_gap_old(b% d_i)) then
            switch_donor = .true.
         else if (b% mtransfer_rate > 0d0) then
            switch_donor = .true.
            b% mtransfer_rate = - b% mtransfer_rate
            mdot_hi_temp = b% mdot_hi
            b% mdot_hi = - b% mdot_lo
            b% mdot_lo = - mdot_hi_temp
            if (.not. b% have_mdot_lo) then
               b% have_mdot_hi = .false.  
            end if
            b% have_mdot_lo = .true.
            b% fixed_delta_mdot = b% fixed_delta_mdot / 2.0
         else if (b% mdot_scheme == "contact" .and. &
            b% rl_relative_gap_old(b% a_i) > b% rl_relative_gap_old(b% d_i) .and. &
            b% rl_relative_gap_old(b% a_i) < - b% implicit_scheme_tolerance .and. &
            b% rl_relative_gap_old(b% d_i) < - b% implicit_scheme_tolerance .and. &
            abs(b% mtransfer_rate/(Msun/secyer)) < b% mdot_limit_donor_switch) then
            switch_donor = .true.
         end if

         if (switch_donor) then
            if (b% report_rlo_solver_progress) write(*,*) "switching donor"
            if (b% d_i == 2) then
               b% d_i = 1
               b% a_i = 2
               b% s_donor => b% s1
               b% s_accretor => b% s2
            else
               b% d_i = 2
               b% a_i = 1
               b% s_donor => b% s2
               b% s_accretor => b% s1
            end if
         end if
      end subroutine

      subroutine binary_evolve_step(b)
         use utils_lib, only: is_bad
         use binary_jdot, only: get_jdot
         use binary_edot, only: get_edot
         type(binary_info), pointer :: b
         integer :: i
         
         include 'formats.inc'
         b% m(b% d_i) = b% s_donor% mstar
         b% time_step = b% s_donor% time_step
         if (b% point_mass_i == 0) then
            b% m(b% a_i) = b% s_accretor% mstar
         else
            b% m(b% a_i) = b% m(b% a_i) - (b% xfer_fraction*b% mtransfer_rate + &
               b% mdot_wind_transfer(b% d_i)) * b% s_donor% dt
         end if
         
         if (b% point_mass_i /= 1) then
            b% r(1) = Rsun*b% s1% photosphere_r
         else
            b% r(1) = 0
         end if
         if (b% point_mass_i /= 2) then
            b% r(2) = Rsun*b% s2% photosphere_r
         else
            b% r(2) = 0
         end if

         ! solve the winds in the system for jdot calculation,
         ! these don't include mass lost due to mass_transfer_efficiency < 1.0
         ! Since s% mstar_dot is not just mass loss, but includes the contribution from
         ! RLO mass transfer and wind mass transfer from the other component, these
         ! need to be removed to get the actual wind. Also, the fraction of the wind
         ! that is fransferred to the other component does not leave the system, and
         ! needs to be removed as well.
         b% mdot_system_wind(b% d_i) = b% s_donor% mstar_dot - b% mtransfer_rate &
            + b% mdot_wind_transfer(b% a_i) - b% mdot_wind_transfer(b% d_i)
         if (b% point_mass_i == 0) then
            b% mdot_system_wind(b% a_i) = b% s_accretor% mstar_dot &
                + b% mtransfer_rate * b% xfer_fraction + b% mdot_wind_transfer(b% d_i) &
                - b% mdot_wind_transfer(b% a_i)
         else
            b% mdot_system_wind(b% a_i) = 0.0d0
         end if

         ! get jdot and update orbital J
         b% jdot = get_jdot(b, b% mtransfer_rate, b% xfer_fraction)
         b% angular_momentum_j = b% angular_momentum_j + b% jdot*b% time_step*secyer

         if (b% angular_momentum_j <= 0) then
            stop 'bad angular_momentum_j'
         end if
         
         ! update the eccentricity (ignore in first step)
         if (.not. b% doing_first_model_of_run) then
            b% eccentricity = b% eccentricity + get_edot(b) *b% time_step*secyer
            if (b% eccentricity < b% min_eccentricity) b% eccentricity = b% min_eccentricity
            if (b% eccentricity > b% max_eccentricity) b% eccentricity = b% max_eccentricity
         end if
         
         !use new eccentricity to calculate new time coordinate
         do i = 1,b% anomaly_steps ! time between periastron and polar angle theta
            b% time_co(i) = ( 2 * atan_cr( sqrt( (1-b% eccentricity)/(1 + b% eccentricity) ) * &
                            tan_cr(b% theta_co(i)/2d0) ) - b% eccentricity * &
                            sqrt(1 - b% eccentricity**2) * sin_cr(b% theta_co(i)) / &
                            (1 + b% eccentricity * cos_cr(b% theta_co(i)) ) ) /2. /pi
            if (i > b% anomaly_steps/2+1) then
               b% time_co(i) = b% time_co(i) + b% time_co(b% anomaly_steps/2+1) * 2
            end if
         end do
         
         ! use the new j to calculate new separation
         b% separation = ((b% angular_momentum_j/(b% m(1)*b% m(2)))**2) *&
             (b% m(1)+b% m(2)) / b% s_donor% cgrav(1) * 1 / (1 - b% eccentricity**2)
         if (b% separation < b% min_binary_separation) &
            b% min_binary_separation = b% separation
         
         b% period = 2*pi*sqrt(pow3(b% separation)/&
               (b% s_donor% cgrav(1)*(b% m(1)+b% m(2)))) 
         if (b% period < min_binary_period) min_binary_period = b% period
         
         ! use the new separation to calculate the new roche lobe radius
         
!          b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
!          b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
         call eval_rlobe(b, 1, 0d0)
         call eval_rlobe(b, 2, 0d0)
         call eval_rlobe_L3(b, 1, 0d0)
         call eval_rlobe_L3(b, 2, 0d0)
         
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1) * (1 - b% eccentricity) ) / &
             b% rl(1) / (1 - b% eccentricity) ! gap < 0 means out of contact 
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2) * (1 - b% eccentricity) ) / &
             b% rl(2) / (1 - b% eccentricity) ! gap < 0 means out of contact

         if (is_bad(b% rl_relative_gap(1)) .or. is_bad(b% rl_relative_gap(2))) then
            stop 'error solving rl_rel_gap'
         end if

         b% model_number = b% model_number + 1
         b% binary_age = b% binary_age + b% time_step

      end subroutine

      integer function binary_check_model(b)
         use binary_mdot, only: rlo_mdot, check_implicit_rlo
         use binary_irradiation
         type (binary_info), pointer :: b

         integer :: i, j, ierr, id
         logical :: implicit_rlo
         real(dp) :: new_mdot, q


         include 'formats.inc'

         binary_check_model = retry
         ierr = 0
         
         implicit_rlo = (b% max_tries_to_achieve > 0 .and. b% implicit_scheme_tolerance > 0d0)
         
         binary_check_model = keep_going
                  
         if (.not. b% ignore_rlof) then
            if (implicit_rlo) then ! check agreement between new r and new rl
               b% s_donor% min_abs_mdot_for_change_limits = 1d99
               binary_check_model = check_implicit_rlo(b, new_mdot)
               if (binary_check_model == keep_going) then
                  b% donor_started_implicit_wind = .false.
               end if
               b% donor_started_implicit_wind = b% donor_started_implicit_wind .or. &
                  b% s_donor% was_in_implicit_wind_limit
            else
               if (.not. b% use_other_rlo_mdot) then
                  call rlo_mdot(b% binary_id, new_mdot, ierr) ! grams per second
                  if (ierr /= 0) then
                     write(*,*) 'failed in rlo_mdot'
                     binary_check_model = retry
                     return
                  end if
               else
                  call b% other_rlo_mdot(b% binary_id, new_mdot, ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in other rlo_mdot'
                     binary_check_model = retry
                     return
                  end if
               end if
               if (new_mdot > 0) then
                  new_mdot = 0.0d0
                  write(*,*) "WARNING: explicit computation of mass transfer results in accreting donor"
                  write(*,*) "Not transfering mass"
               end if
               ! smooth out the changes in mdot
               write (*,*) 'new_mdot original', new_mdot
               new_mdot = b% cur_mdot_frac*b% mtransfer_rate + (1-b% cur_mdot_frac)*new_mdot
               if (-new_mdot/(Msun/secyer) > b% max_explicit_abs_mdot) new_mdot = -b% max_explicit_abs_mdot*Msun/secyer 
            end if
            b% mtransfer_rate = new_mdot
         else
            b% mtransfer_rate = 0
         end if
         call adjust_irradiation(b, b% mtransfer_rate, b% xfer_fraction)

         if (b% point_mass_i == 0 .and. b% rl_relative_gap(b% a_i) >= 0.0d0) then
            if (b% rl_relative_gap(b% a_i) >= b% accretor_overflow_terminate) then
               binary_check_model = terminate
               b% s_donor% termination_code = t_xtra1
               termination_code_str(t_xtra1) = &
                   "Terminate because accretor (r-rl)/rl > accretor_overflow_terminate"
            end if
         end if
         if (b% doing_first_model_of_run .and. b% terminate_if_initial_overflow .and. .not. b% ignore_rlof) then
            if (b% rl_relative_gap(b% d_i) >= 0.0d0 .or. (b% point_mass_i == 0 .and. b% rl_relative_gap(b% a_i) >= 0.0d0)) then
               binary_check_model = terminate
               b% s_donor% termination_code = t_xtra1
               termination_code_str(t_xtra1) = &
                   "Terminate because of overflowing initial model"
            end if
         end if
         if (b% point_mass_i == 0 .and. b% terminate_if_L2_overflow) then
            if (b% m(1) > b% m(2)) then
               q = b% m(2) / b% m(1)
               id = 2
            else
               q = b% m(1) / b% m(2)
               id = 1
            end if
            if (b% rl_relative_gap(id) > 0.29858997d0*atan_cr(1.83530121d0*pow_cr(q,0.39661426d0))) then
               binary_check_model = terminate
               b% s_donor% termination_code = t_xtra1
               termination_code_str(t_xtra1) = &
                   "Terminate because of L2 overflow"
            end if
         end if

      end function binary_check_model

!       real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
!          real(dp), intent(in) :: m1, m2, a
!          real(dp) :: q
!          q = pow_cr(m1/m2,one_third)
!       ! Roche lobe size for star of mass m1 with a
!       ! companion of mass m2 at separation a, according to
!       ! the approximation of Eggleton 1983, apj 268:368-369
!          rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
!       end function eval_rlobe

      integer function binary_finish_step(b)
         type (binary_info), pointer :: b
         real(dp) :: spin_period

         binary_finish_step = keep_going
         ! update change factor in case mtransfer_rate has changed
         if(b% mtransfer_rate_old /= b% mtransfer_rate .and. &
             b% mtransfer_rate /= 0 .and. b% mtransfer_rate_old/=0.0) then
            if(b% mtransfer_rate < b% mtransfer_rate_old) then
               b% change_factor = b% change_factor*(1.0-b% implicit_lambda) + b% implicit_lambda* &
                  (1+b% change_factor_fraction*(b% mtransfer_rate/b% mtransfer_rate_old-1))
            else
               b% change_factor = b% change_factor*(1.0-b% implicit_lambda) + b% implicit_lambda* &
                   (1+b% change_factor_fraction*(b% mtransfer_rate_old/b% mtransfer_rate-1))
            end if
            if(b% change_factor > b% max_change_factor) b% change_factor = b% max_change_factor
            if(b% change_factor < b% min_change_factor) b% change_factor = b% min_change_factor
         end if

         ! store all variables into "old" and "older"
         b% model_number_older = b% model_number_old
         b% binary_age_older = b% binary_age_old
         b% mtransfer_rate_older = b% mtransfer_rate_old
         b% angular_momentum_j_older = b% angular_momentum_j_old
         b% separation_older = b% separation_old
         b% eccentricity_older = b% eccentricity_old
         b% dt_older = b% dt_old
         b% env_older = b% env_old
         b% xfer_fraction_older = b% xfer_fraction_old
         b% sum_div_qloc_older(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc_older(2) = b% sum_div_qloc_old(2)
         b% period_older = b% period_old
         b% rl_relative_gap_older(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap_older(2) = b% rl_relative_gap_old(2)
         b% r_older(1) = b% r_old(1)
         b% r_older(2) = b% r_old(2)
         b% rl_older(1) = b% rl_old(1)
         b% rl_older(2) = b% rl_old(2)
         b% rl3_older(1) = b% rl3_old(1)
         b% rl3_older(2) = b% rl3_old(2)
         b% m_older(1) = b% m_old(1)
         b% m_older(2) = b% m_old(2)
         b% have_radiative_core_older = b% have_radiative_core_old
         b% max_timestep_older = b% max_timestep_old
         b% change_factor_older = b% change_factor_old

         b% d_i_older = b% d_i_old
         b% a_i_older = b% a_i_old
         b% point_mass_i_older = b% point_mass_i_old

         b% dt_why_reason_older = b% dt_why_reason_old

         b% model_number_old = b% model_number
         b% binary_age_old = b% binary_age
         b% mtransfer_rate_old = b% mtransfer_rate
         b% angular_momentum_j_old = b% angular_momentum_j
         b% separation_old = b% separation
         b% eccentricity_old = b% eccentricity
         b% dt_old = b% dt
         b% env_old = b% env
         b% xfer_fraction_old = b% xfer_fraction
         b% sum_div_qloc_old(1) = b% sum_div_qloc(1)
         b% sum_div_qloc_old(2) = b% sum_div_qloc(2)
         b% period_old = b% period
         b% rl_relative_gap_old(1) = b% rl_relative_gap(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap(2)
         b% r_old(1) = b% r(1)
         b% r_old(2) = b% r(2)
         b% rl_old(1) = b% rl(1)
         b% rl_old(2) = b% rl(2)
         b% rl3_old(1) = b% rl3(1)
         b% rl3_old(2) = b% rl3(2)
         b% m_old(1) = b% m(1)
         b% m_old(2) = b% m(2)
         b% have_radiative_core_old = b% have_radiative_core
         b% max_timestep_old = b% max_timestep
         b% change_factor_old = b% change_factor

         b% d_i_old = b% d_i
         b% a_i_old = b% a_i
         b% point_mass_i_old = b% point_mass_i

         b% dt_why_reason_old = b% dt_why_reason

      end function binary_finish_step

      integer function binary_prepare_to_redo(b)
         type (binary_info), pointer :: b

         binary_prepare_to_redo = redo
         call binary_set_current_to_old(b)

      end function binary_prepare_to_redo

      integer function binary_prepare_to_retry(b)
         type (binary_info), pointer :: b

         b% num_tries = 0
         ! this call takes care of restoring variables
         call binary_set_current_to_old(b)
         binary_prepare_to_retry = retry

      end function binary_prepare_to_retry

      integer function binary_do1_backup(b)
         type (binary_info), pointer :: b

         binary_do1_backup = backup
         if (b% s_donor% generations > 2) then
            call binary_restore_older(b)
         end if
         call binary_set_current_to_old(b)

         b% num_tries = 0

      end function binary_do1_backup

      subroutine binary_set_current_to_old(b)
         type (binary_info), pointer :: b
         ! restore variables
         b% model_number = b% model_number_old
         b% binary_age = b% binary_age_old
         ! do not restore mtransfer_rate during implicit rlo
         if (b% num_tries == 0) b% mtransfer_rate = b% mtransfer_rate_old
         b% angular_momentum_j = b% angular_momentum_j_old
         b% separation = b% separation_old
         b% eccentricity = b% eccentricity_old
         b% dt = b% dt_old
         b% env = b% env_old
         b% xfer_fraction = b% xfer_fraction_old
         b% sum_div_qloc(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc(2) = b% sum_div_qloc_old(2)
         b% period = b% period_old
         b% rl_relative_gap(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap(2) = b% rl_relative_gap_old(2)
         b% r(1) = b% r_old(1)
         b% r(2) = b% r_old(2)
         b% rl(1) = b% rl_old(1)
         b% rl(2) = b% rl_old(2)
         b% rl3(1) = b% rl3_old(1)
         b% rl3(2) = b% rl3_old(2)
         b% m(1) = b% m_old(1)
         b% m(2) = b% m_old(2)
         b% have_radiative_core = b% have_radiative_core_old
         b% max_timestep = b% max_timestep_old
         ! do not restore change_factor during implicit mdot
         if (b% num_tries == 0) b% change_factor = b% change_factor_old

         ! do not restore donor and accretor ids during implicit mdot
         if (b% num_tries == 0) b% d_i = b% d_i_old
         if (b% num_tries == 0) b% a_i = b% a_i_old

         b% point_mass_i = b% point_mass_i_old

         b% dt_why_reason = b% dt_why_reason_old
      end subroutine binary_set_current_to_old

      subroutine binary_restore_older(b)
         type (binary_info), pointer :: b

         b% model_number_old = b% model_number_older
         b% mtransfer_rate_old = b% mtransfer_rate_older
         b% angular_momentum_j_old = b% angular_momentum_j_older
         b% separation_old = b% separation_older
         b% eccentricity_old = b% eccentricity_older
         b% dt_old = b% dt_older
         b% env_old = b% env_older
         b% xfer_fraction_old = b% xfer_fraction_older
         b% sum_div_qloc_old(1) = b% sum_div_qloc_older(1)
         b% sum_div_qloc_old(2) = b% sum_div_qloc_older(2)
         b% period_old = b% period_older
         b% rl_relative_gap_old(1) = b% rl_relative_gap_older(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap_older(2)
         b% r_old(1) = b% r_older(1)
         b% r_old(2) = b% r_older(2)
         b% rl_old(1) = b% rl_older(1)
         b% rl_old(2) = b% rl_older(2)
         b% rl3_old(1) = b% rl3_older(1)
         b% rl3_old(2) = b% rl3_older(2)
         b% m_old(1) = b% m_older(1)
         b% m_old(2) = b% m_older(2)
         b% have_radiative_core_old = b% have_radiative_core_older
         b% max_timestep_old = b% max_timestep_older
         b% change_factor_old = b% change_factor_older

         b% d_i_old = b% d_i_older
         b% a_i_old = b% a_i_older
         b% point_mass_i_old = b% point_mass_i_older

         b% dt_why_reason_old = b% dt_why_reason_older

      end subroutine binary_restore_older

      end module binary_evolve
