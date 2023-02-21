! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton and Pablo Marchant
!
!   this file is part of mesa.
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
 
      module run_binary_support

      use star_lib
      use star_def
      use const_def
      use utils_lib
      use binary_def
      use binary_private_def
      use binary_ctrls_io, only: do_one_binary_setup
      use binary_do_one_utils

      
      implicit none

      contains

      subroutine do_run1_binary(tst, &
            ! star extras
            extras_controls, &
            ! binary extras
            extras_binary_controls, &
            ierr, &
            inlist_fname_arg)

         use binary_job_ctrls_io
         use binary_mdot, only: eval_xfer_fraction, adjust_mdots, set_accretion_composition
         use binary_wind, only: eval_wind_xfer_fractions
         use binary_tides, only: sync_spin_orbit_torque
         use binary_evolve
         use mod_other_rlo_mdot
         use mod_other_tsync
         use mod_other_sync_spin_to_orbit
         use mod_other_mdot_edd
         use mod_other_accreted_material_j
         use mod_other_binary_jdot
         use mod_other_binary_wind_transfer
         use mod_other_binary_edot
         use mod_other_binary_extras
         use binary_timestep
         use binary_history
         use binary_history_specs
         use run_star_support
         
         logical, intent(in) :: tst
         
         interface

            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      

            subroutine extras_binary_controls(binary_id, ierr)
               integer :: binary_id
               integer, intent(out) :: ierr
            end subroutine extras_binary_controls      

         end interface
         
         integer, intent(out) :: ierr
         character (len=*) :: inlist_fname_arg
         optional inlist_fname_arg
         

         integer :: id, id_extra, i, j, k, l, i_prev, result, partial_result, result_reason, model_number, iounit,binary_startup
         type (star_info), pointer :: s
         character (len=256) :: restart_filename, photo_filename
         integer :: total, time0, time1, clock_rate
         integer :: model
         logical :: doing_restart, first_try, continue_evolve_loop, just_did_backup
         logical :: get_history_info, write_history, write_terminal
         real(dp) :: sum_times
         real(dp) :: dt
         real(dp) :: timestep_factor
         integer :: binary_id
         type (binary_info), pointer :: b
         character (len=strlen) :: inlist_fname

         include 'formats.inc'

         ierr = 0
         id_extra = 0
         call system_clock(time0,clock_rate)

         call resolve_inlist_fname(inlist_fname,inlist_fname_arg)
         MESA_INLIST_RESOLVED=.true. ! Now any call to resolve_inlist_fname will only return inlist_fname_arg or 'inlist'

         ! Find out if this is a restart
         iounit=alloc_iounit(ierr)
         open(unit=iounit, file='.restart', status='old', action='read',iostat=ierr)
         doing_restart = (ierr == 0)
         if (doing_restart) then
             read(iounit,'(a)', iostat=ierr) photo_filename ! same for both stars
             if (ierr /= 0) then
                 stop "Problem while reading restart info"
             end if
         else
             ierr = 0
         end if
         call free_iounit(iounit)

         binary_id = alloc_binary(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in alloc_binary'
            return
         end if

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (.not. doing_restart) then
             b% model_number = 0
             b% model_number_old = 0
             b% model_number_older = 0
             b% binary_age = 0
             b% binary_age_old = 0
             b% binary_age_older = 0
         end if

         call do_read_binary_job(b, inlist_fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_binary_job'
            return
         end if

         if (.not. b% job% evolve_both_stars) num_stars=1
         
         write(*,*)
         write(*,*)

         result_reason = 0

         ! Setup null hooks
         b% other_rlo_mdot => null_other_rlo_mdot
         b% other_tsync => null_other_tsync
         b% other_sync_spin_to_orbit => null_other_sync_spin_to_orbit
         b% other_mdot_edd => null_other_mdot_edd
         b% other_accreted_material_j => null_other_accreted_material_j
         b% other_jdot_gr => null_other_jdot_gr
         b% other_jdot_ml => null_other_jdot_ml
         b% other_jdot_ls => null_other_jdot_ls
         b% other_jdot_missing_wind => null_other_jdot_missing_wind
         b% other_jdot_mb => null_other_jdot_mb
         b% other_extra_jdot => null_other_extra_jdot
         b% other_binary_wind_transfer => null_other_binary_wind_transfer
         b% other_edot_tidal => null_other_edot_tidal
         b% other_edot_enhance => null_other_edot_enhance
         b% other_extra_edot => null_other_extra_edot

         b% extras_binary_startup => null_extras_binary_startup
         b% extras_binary_check_model => null_extras_binary_check_model
         b% extras_binary_finish_step => null_extras_binary_finish_step
         b% extras_binary_after_evolve => null_extras_binary_after_evolve
         b% how_many_extra_binary_history_columns => null_how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => null_data_for_extra_binary_history_columns
         
         
         call do_one_binary_setup(b, inlist_fname, ierr)
         ! extras_binary_controls is defined in run_binary_extras.f and hooks can
         ! be specified there
         call extras_binary_controls(b% binary_id, ierr)

         b% donor_id = -1
         b% accretor_id = -1
         do i = 1, num_stars
         
            call do_read_star_job(b% job% inlist_names(i), ierr)
            if (failed('do_read_star_job',ierr)) return

            !id = id_from_read_star_job ! star allocated by do_read_star_job
            !id_from_read_star_job = 0

            restart_filename  = trim(photo_filename)

            ! the star is initialized in this call
            call before_evolve_loop(.true., doing_restart, doing_restart, &
               binary_controls, extras_controls, &
               id_from_read_star_job, b% job% inlist_names(i), restart_filename, &
               .false., binary_id, id, id_extra, ierr)

            if (ierr /= 0) return

            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            b% star_ids(i) = id
            b% star_extra_ids(i) = id_extra

            s% include_binary_history_in_log_file = b% append_to_star_history

            ! fix photo output for both stars to the one defined by binary
            s% photo_interval = b% photo_interval
            s% photo_digits = b% photo_digits

            s% how_many_binary_history_columns => how_many_binary_history_columns
            s% data_for_binary_history_columns => data_for_binary_history_columns

            ! additional settings for mass transfer and tides
            if (b% do_j_accretion) then
               s% use_accreted_material_j = .true.
            end if
            s% accrete_given_mass_fractions = .true.
            s% accrete_same_as_surface = .false.
            s% binary_other_torque => sync_spin_orbit_torque
            
            s% doing_timing = .false.
            
            write(*,*)
            write(*,*)

         end do

         b% evolve_both_stars = b% job% evolve_both_stars
         b% warn_binary_extra = b% job% warn_binary_extra
         if (.not. b% evolve_both_stars) then
            b% point_mass_i = 2
         else
            b% point_mass_i = 0
         end if

         ! Initialize arrays for phase dependent calculations
         allocate(b% theta_co(b% anomaly_steps), b% time_co(b% anomaly_steps), &
            b% mdot_donor_theta(b% anomaly_steps))
         allocate(b% edot_theta(b% anomaly_steps), b% e1(b% anomaly_steps), &
            b% e2(b% anomaly_steps), b% e3(b% anomaly_steps))
         ! binary data must be initiated after stars, such that masses are available
         ! if using saved models
         if (.not. doing_restart) then
            call binarydata_init(b)
         else
            if (b% d_i == b% point_mass_i) then
               write(*,*) "WARNING: restart has donor star set as point mass"
               write(*,*) "switching donor"
               b% d_i = b% a_i
               b% a_i = b% point_mass_i
            end if
            if (b% d_i == 1) then
               b% s_donor => b% s1
               b% s_accretor => b% s2
            else
               b% s_donor => b% s2
               b% s_accretor => b% s1
            end if
         end if
         call binary_private_def_init
         call binary_history_column_names_init(ierr)
         call set_binary_history_columns(b, b% job% binary_history_columns_file, ierr)

         if (b% job% show_binary_log_description_at_start .and. .not. doing_restart) then
            write(*,*)
            call do_show_binary_log_description(id, ierr)
            if (failed('show_log_description',ierr)) return
         end if

         if (b% point_mass_i /= 1 .and. b% do_initial_orbit_sync_1 .and. &
             .not. doing_restart) then
            call star_relax_uniform_omega( &
               b% s1% id, 0, (2*pi) / b% period, b% s1% job% num_steps_to_relax_rotation,&
               b% s1% job% relax_omega_max_yrs_dt, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in initial orbital sync'
               return
            end if
         end if

         if (b% point_mass_i /= 2 .and. b% do_initial_orbit_sync_2 .and. &
             .not. doing_restart) then
            call star_relax_uniform_omega( &
               b% s2% id, 0, (2*pi) / b% period, b% s2% job% num_steps_to_relax_rotation,&
               b% s2% job% relax_omega_max_yrs_dt, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in initial orbital sync'
               return
            end if
         end if

         continue_evolve_loop = .true.
         s% doing_timing = .false.
         i_prev = 0
         
         binary_startup = b% extras_binary_startup(b% binary_id,doing_restart,ierr)
         if (ierr /= 0) return

         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            if (b% point_mass_i /= 0) then
               num_stars = 1
            else
               num_stars = 2
            end if
            do l = 1, num_stars
               if (l == 1 .and. b% point_mass_i == 1) then
                  i = 2
               else
                  i = l
               end if

               id = b% star_ids(i)
               call star_ptr(id, s, ierr)
               if (failed('star_ptr',ierr)) return
               call before_step_loop(s, ierr)
               if (ierr /= 0) return
            end do
             
            first_try = .true.
            just_did_backup = .false.
            b% donor_started_implicit_wind = .false.
            b% num_tries = 0
            step_loop: do ! may need to repeat this loop

               call set_donor_star(b)
               call set_star_timesteps(b)
               result = keep_going
               
               !get transfer fraction
               call eval_xfer_fraction(b)
               call eval_wind_xfer_fractions(b)

               ! Store mtransfer_rate used in a step, as it is rewritten by binary_check_model and
               ! that can produce inconsistent output.
               b% step_mtransfer_rate = b% mtransfer_rate

               ! if donor reaches implicit wind limit during mass transfer, set was_in_implicit_wind_limit = .true.
               ! so that further iterations of the implicit wind do not start far off from the required mdot.
               ! system will likely detach suddenly, so set change_factor to double its maximum
               if (b% donor_started_implicit_wind) then
                  b% s_donor% was_in_implicit_wind_limit = .true.
                  b% change_factor = 2*b% max_change_factor 
               end if

               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if

                  ! evolve donor first
                  if (i == 1 .or. b% point_mass_i /= 0) then
                     j = b% d_i
                  else
                     j = b% a_i
                  end if
                  id = b% star_ids(j)

                  ! Avoid repeting the accretor when using the implicit scheme plus
                  ! rotation and implicit winds. When this happens the accretor won't
                  ! usually care about the result of the evolution of the donor.
                  ! EXPERIMENTAL
                  if (j == b% a_i .and. b% num_tries >0 .and. s% was_in_implicit_wind_limit) &
                      cycle

                  ! fix sync timescales to zero. If synching stars these will be
                  ! updated at each star_evolve_step
                  if (j == 1) b% t_sync_1 = 0
                  if (j == 2) b% t_sync_2 = 0
                  
                  result = worst_result(result, star_evolve_step_part1(id, first_try, just_did_backup))

               end do

               ! modify mdots to account for mass transfer
               call adjust_mdots(b)

               ! if both stars are accreting mass at this point, then there is
               ! something very wrong! If one star loses and the other gains mass,
               ! then the mass losing star must be evolved first
               k = b% d_i
               if (b% point_mass_i == 0 .and. b% s_donor% mstar_dot > 0 .and. &
                   b% s_accretor% mstar_dot > 0) then
                  write(*,*) "ERROR: both stars accreting, terminating evolution"
                  result = terminate
                  exit step_loop
               else if (b% point_mass_i /= 0 .and. b% s_donor% mstar_dot > 0) then
                  write(*,*) "ERROR: donor accreting, terminating evolution"
                  result = terminate
                  exit step_loop
               else
                  ! donor can end up accreting due to winds
                  ! NOTE: this is a bit confusing, need to modify things so the
                  ! donor is simply the accreting star.
                  if (b% point_mass_i == 0 .and. b% s_donor% mstar_dot > 0) &
                     k = b% a_i
               end if

               if (result == keep_going) then
                  do l = 1, num_stars
                     if (l == 1 .and. b% point_mass_i == 1) then
                        i = 2
                     else
                        i = l
                     end if
                     if (i == 1 .and. b% point_mass_i == 1) i = 2

                     ! evolve donor first
                     if (i == 1 .or. b% point_mass_i /= 0) then
                        j = k
                     else
                        ! ids are 1 and 2, so 3-k represents the other star
                        j = 3-k
                     end if
                     id = b% star_ids(j)

                     if (j == b% a_i .and. b% num_tries >0 .and. s% was_in_implicit_wind_limit) &
                         cycle

                     ! set accretion composition
                     if (i == 2) call set_accretion_composition(b, j)
                     result = worst_result(result, star_evolve_step_part2(id, first_try, just_did_backup))

                  end do
               end if

               ! do not evolve binary step on failure, its unnecesary and some variables are not properly
               ! setted when the newton solver fails.
               if (result == keep_going) then
                  call binary_evolve_step(b)
               end if

               if (result == keep_going) then
                  result = worst_result(result, binary_check_model(b))
               end if

               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if
                  if (i == 1 .and. b% point_mass_i == 1) i = 2
                  if (result == keep_going) then
                     id = b% star_ids(i)
                     id_extra = b% star_extra_ids(i)
                     call star_ptr(id, s, ierr)
                     result = worst_result(result, s% extras_check_model(id, id_extra))
                     result = worst_result(result, star_check_model(id))
                  end if
               end do
               
               partial_result = b% extras_binary_check_model(b% binary_id)
               result = worst_result(result, partial_result)

               ! solve first binary timestep limit because star_pick_next_timestep needs it
               result = worst_result(result, binary_pick_next_timestep(b))

               if (result == keep_going) then
                  do l = 1, num_stars
                     if (l == 1 .and. b% point_mass_i == 1) then
                        i = 2
                     else
                        i = l
                     end if
                     id = b% star_ids(i)
                     call star_ptr(id, s, ierr)
                     if (failed('star_ptr',ierr)) return
                     result = worst_result(result, star_pick_next_timestep(id))
                  end do
               end if
               if (result == keep_going) then
                  exit step_loop
               end if
               
               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if

                  id = b% star_ids(i)
                  model_number = get_model_number(id, ierr)
                  if (failed('get_model_number',ierr)) return
                  
                  result_reason = get_result_reason(id, ierr)
                  if (result == retry) then
                     if (failed('get_result_reason',ierr)) return
                     if (s% job% report_retries) &
                        write(*,'(i6,3x,a,/)') model_number, &
                           'retry reason ' // trim(result_reason_str(result_reason))
                  else if (result == backup) then
                     if (failed('get_result_reason',ierr)) return
                     if (s% job% report_backups) &
                        write(*,'(i6,3x,a,/)') model_number, &
                           'backup reason ' // trim(result_reason_str(result_reason))
                  end if

               end do
               
               if (result == redo) then
                  result = worst_result(result, binary_prepare_to_redo(b))
                  do l = 1, num_stars
                     if (l == 1 .and. b% point_mass_i == 1) then
                        i = 2
                     else
                        i = l
                     end if
                     id = b% star_ids(i)

                     ! Avoid repeting the accretor when using the implicit scheme plus
                     ! rotation and implicit winds. When this happens the accretor won't
                     ! usually care about the result of the evolution of the donor.
                     ! EXPERIMENTAL

                     if (i == b% a_i .and. b% num_tries >0 .and. s% was_in_implicit_wind_limit) &
                         cycle
                     result = worst_result(result, star_prepare_to_redo(id))
                  end do
               end if
               if (result == retry) then
                  result = worst_result(result, binary_prepare_to_retry(b))
                  do l = 1, num_stars
                     if (l == 1 .and. b% point_mass_i == 1) then
                        i = 2
                     else
                        i = l
                     end if
                     id = b% star_ids(i)
                     result = worst_result(result, star_prepare_to_retry(id))
                  end do
               end if
               if (result == backup) then
                  result = worst_result(result, binary_do1_backup(b))
                  do l = 1, num_stars
                     if (l == 1 .and. b% point_mass_i == 1) then
                        i = 2
                     else
                        i = l
                     end if
                     id = b% star_ids(i)
                     result = worst_result(result, star_do1_backup(id))
                  end do
                  just_did_backup = .true.
               else
                  just_did_backup = .false.
               end if
               if (result == terminate) then
                  continue_evolve_loop = .false.
                  exit step_loop
               end if
               first_try = .false.
               
            end do step_loop

            if(result == keep_going) result = binary_finish_step(b)
            
            partial_result=b% extras_binary_finish_step(b% binary_id)
            result=worst_result(result,partial_result)

            if (result == keep_going) then
               ! write terminal info
               model = b% model_number
               if (b% history_interval > 0) then
                  write_history = (mod(model, b% history_interval) == 0)
               else
                  write_history = .false.
               end if
               if (s% terminal_interval > 0) then
                  write_terminal = (mod(model, b% terminal_interval) == 0)
               else
                  write_terminal = .false.
               end if
               if (write_history) b% need_to_update_binary_history_now = .true.
               get_history_info = b% need_to_update_binary_history_now .or. write_terminal
               if (get_history_info) then
                  if (b% write_header_frequency*b% terminal_interval > 0) then
                     if ( mod(model, b% write_header_frequency*b% terminal_interval) .eq. 0 &
                          .and. .not. b% doing_first_model_of_run) then
                        write(*,*)
                        call write_binary_terminal_header(b)
                     end if
                  end if         
                  if (write_terminal) call do_binary_terminal_summary(b)  
                  if (b% need_to_update_binary_history_now) call write_binary_history_info(b, ierr)
                  b% need_to_update_binary_history_now = .false.
               end if
            end if

            do l = 1, num_stars
               if (l == 1 .and. b% point_mass_i == 1) then
                  i = 2
               else
                  i = l
               end if
               id = b% star_ids(i)
               id_extra = b% star_extra_ids(i)
               call star_ptr(id, s, ierr)
               partial_result = result
               call after_step_loop(s, b% job% inlist_names(i), &
                  id_extra, .false., partial_result, ierr)
               if (ierr /= 0) return
               result =  worst_result(result, partial_result)
            end do
               
            if (result /= keep_going) then
               if (result /= terminate) then
                  write(*,2) 'ERROR in result value in run_star_extras: model', &
                     s% model_number
                  write(*,2) 'result', result
                  exit evolve_loop
               end if
               do l = 1, num_stars
                  if (l == 1 .and. b% point_mass_i == 1) then
                     i = 2
                  else
                     i = l
                  end if
                  id = b% star_ids(i)
                  id_extra = b% star_extra_ids(i)
                  call star_ptr(id, s, ierr)
                  if (s% result_reason == result_reason_normal) then

                     partial_result = result
                     call terminate_normal_evolve_loop(s, &
                         id_extra, .false., partial_result, ierr)
                     if (ierr /= 0) return
                     result =  worst_result(result, partial_result)

                  end if
               end do
               call write_binary_history_info(b, ierr)
               call do_binary_terminal_summary(b)
               exit evolve_loop
            end if
                        
            do l = 1, num_stars
               if (l == 1 .and. b% point_mass_i == 1) then
                  i = 2
               else
                  i = l
               end if
               id = b% star_ids(i)
               id_extra = b% star_extra_ids(i)
               call star_ptr(id, s, ierr)
         
               call do_saves( &
                  id, id_extra, s)

               if (s% doing_timing) then
                  call system_clock(s% job% time1_extra,s% job% clock_rate)
                  s% job% after_step_timing = s% job% after_step_timing + &
                     dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
                  s% job% check_time_end = eval_total_times(s% id, ierr)
                  s% job% check_after_step_timing = s% job% check_after_step_timing + &
                     (s% job% check_time_end - s% job% check_time_start)
               end if
            end do
            call do_saves_for_binary_rlo
      

            if (b% doing_first_model_of_run) b% doing_first_model_of_run = .false.
            
         end do evolve_loop

         ! deallocate arrays used for calculation of phase dependent variables
         deallocate(b% theta_co, b% time_co, b% mdot_donor_theta)
         deallocate(b% edot_theta, b% e1, b% e2, b% e3)
         
         call b% extras_binary_after_evolve(b% binary_id,ierr)
         if (ierr /= 0) then
            return
         end if 
         
         do i = 1, 2
         
            id = b% star_ids(i)
            id_extra = b% star_extra_ids(i)
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) then
               ierr = 0
               cycle
            end if

            call after_evolve_loop(s, id_extra, .true., ierr)

            if (s% doing_timing) then
               call system_clock(s% job% time1_extra,s% job% clock_rate)
               s% job% after_step_timing = s% job% after_step_timing + &
                  dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
               s% job% check_time_end = eval_total_times(s% id, ierr)
               s% job% check_after_step_timing = s% job% check_after_step_timing + &
                  (s% job% check_time_end - s% job% check_time_start)
            end if
            
         end do
         
         call starlib_shutdown


         contains
                  
         
         subroutine do_saves_for_binary_rlo

            integer :: io, i
            character (len=strlen) :: str, str1, str2, filename

            if (b% point_mass_i /= 1) then
               i = index(b% s1% most_recent_photo_name,"/",.true.)
               if (i /= 0) then
                  str1 = b% s1% most_recent_photo_name(i:)
               else
                  str1 = b% s1% most_recent_photo_name
               end if
            end if
            if (b% point_mass_i /= 2) then
               i = index(b% s2% most_recent_photo_name,"/",.true.)
               if (i /= 0) then
                  str2 = b% s2% most_recent_photo_name(i:)
               else
                  str2 = b% s2% most_recent_photo_name
               end if
            end if

            if (b% point_mass_i == 0) then
               str = str1
               if (str1 /= str2) write(*,*) "WARNING: photos off sync"
            else if (b% point_mass_i == 1) then
               str = str2
            else if (b% point_mass_i == 2) then
               str = str1
            end if
                  
            if (b% last_photo_filename /= str) then
         
               b% last_photo_filename = str
               io = alloc_iounit(ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in alloc_iounit'
                  return
               end if
               filename = '.restart'
               open(unit=io, file=trim(filename), action='write', iostat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to open ' // trim(filename)
                  return
               end if
               write(io,'(a)') str
               close(io)
               call free_iounit(io)
               
            end if
            
         end subroutine do_saves_for_binary_rlo


      end subroutine do_run1_binary   

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
      
      
      subroutine binary_controls(id, binary_id, ierr)
         integer, intent(in) :: id, binary_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         s% binary_id = binary_id

         if (b% donor_id == -1) then
            b% donor_id = id
            b% s1 => s
            s% initial_mass = b% m1
            s% other_photo_read => binary_photo_read
            s% other_photo_write => binary_photo_write
         else
            b% accretor_id = id
            b% s2 => s
            s% initial_mass = b% m2
         end if
      end subroutine binary_controls

      subroutine binary_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         type(binary_info), pointer :: b

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         read(iounit, iostat=ierr) &
             b% model_number, b% model_number_old, b% model_number_older, &
             b% mtransfer_rate, b% mtransfer_rate_old, b% mtransfer_rate_older, &
             b% angular_momentum_j, b% angular_momentum_j_old, b% angular_momentum_j_older, & 
             b% separation, b% separation_old, b% separation_older, &
             b% eccentricity, b% eccentricity_old, b% eccentricity_older, &
             b% rl_relative_gap(1), b% rl_relative_gap_old(1), b% rl_relative_gap_older(1), &
             b% rl_relative_gap(2), b% rl_relative_gap_old(2), b% rl_relative_gap_older(2), &
             b% r(1), b% r_old(1), b% r_older(1), &
             b% r(2), b% r_old(2), b% r_older(2), &
             b% rl(1), b% rl_old(1), b% rl_older(1), &
             b% rl(2), b% rl_old(2), b% rl_older(2), &
             b% m(1), b% m_old(1), b% m_older(1), &
             b% m(2), b% m_old(2), b% m_older(2), &
             b% sum_div_qloc(1), b% sum_div_qloc_old(1), b% sum_div_qloc_older(1), &
             b% sum_div_qloc(2), b% sum_div_qloc_old(2), b% sum_div_qloc_older(2), &
             b% dt, b% dt_old, b% dt_older, &
             b% env, b% env_old, b% env_older, &
             b% xfer_fraction, b% xfer_fraction_old, b% xfer_fraction_older, &
             b% eq_initial_bh_mass, &
             b% period, b% period_old, b% period_older, & 
             b% max_timestep, b% max_timestep_old, b% max_timestep_older, &
             b% change_factor, b% change_factor_old, b% change_factor_older, &
             b% min_binary_separation, &
             b% have_radiative_core(1), b% have_radiative_core_old(1), b% have_radiative_core_older(1), &
             b% have_radiative_core(2), b% have_radiative_core_old(2), b% have_radiative_core_older(2), &
             b% d_i, b% d_i_old, b% d_i_older, b% a_i, b% a_i_old, b% a_i_older, &
             b% point_mass_i, b% point_mass_i_old, b% point_mass_i_older, &
             b% dt_why_reason, b% dt_why_reason_old, b% dt_why_reason_older
         if (ierr /= 0) stop "error in binary_photo_read"
      end subroutine binary_photo_read

      subroutine binary_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         type(star_info), pointer :: s
         type(binary_info), pointer :: b

         integer :: ierr

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         write(iounit) &
             b% model_number, b% model_number_old, b% model_number_older, &
             b% mtransfer_rate, b% mtransfer_rate_old, b% mtransfer_rate_older, &
             b% angular_momentum_j, b% angular_momentum_j_old, b% angular_momentum_j_older, & 
             b% separation, b% separation_old, b% separation_older, &
             b% eccentricity, b% eccentricity_old, b% eccentricity_older, &
             b% rl_relative_gap(1), b% rl_relative_gap_old(1), b% rl_relative_gap_older(1), &
             b% rl_relative_gap(2), b% rl_relative_gap_old(2), b% rl_relative_gap_older(2), &
             b% r(1), b% r_old(1), b% r_older(1), &
             b% r(2), b% r_old(2), b% r_older(2), &
             b% rl(1), b% rl_old(1), b% rl_older(1), &
             b% rl(2), b% rl_old(2), b% rl_older(2), &
             b% m(1), b% m_old(1), b% m_older(1), &
             b% m(2), b% m_old(2), b% m_older(2), &
             b% sum_div_qloc(1), b% sum_div_qloc_old(1), b% sum_div_qloc_older(1), &
             b% sum_div_qloc(2), b% sum_div_qloc_old(2), b% sum_div_qloc_older(2), &
             b% dt, b% dt_old, b% dt_older, &
             b% env, b% env_old, b% env_older, &
             b% xfer_fraction, b% xfer_fraction_old, b% xfer_fraction_older, &
             b% eq_initial_bh_mass, &
             b% period, b% period_old, b% period_older, & 
             b% max_timestep, b% max_timestep_old, b% max_timestep_older, &
             b% change_factor, b% change_factor_old, b% change_factor_older, &
             b% min_binary_separation, &
             b% have_radiative_core(1), b% have_radiative_core_old(1), b% have_radiative_core_older(1), &
             b% have_radiative_core(2), b% have_radiative_core_old(2), b% have_radiative_core_older(2), &
             b% d_i, b% d_i_old, b% d_i_older, b% a_i, b% a_i_old, b% a_i_older, &
             b% point_mass_i, b% point_mass_i_old, b% point_mass_i_older, &
             b% dt_why_reason, b% dt_why_reason_old, b% dt_why_reason_older

      end subroutine binary_photo_write

      end module run_binary_support
