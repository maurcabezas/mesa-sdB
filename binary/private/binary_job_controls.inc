      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! binary_job controls

      logical :: show_binary_log_description_at_start
      character (len=1000) :: binary_history_columns_file
      logical :: warn_binary_extra
      character (len=256), dimension(2) :: inlist_names
      logical :: read_extra_binary_job_inlist1
      character (len=256) :: extra_binary_job_inlist1_name 
      logical :: read_extra_binary_job_inlist2
      character (len=256) :: extra_binary_job_inlist2_name 
      logical :: read_extra_binary_job_inlist3
      character (len=256) :: extra_binary_job_inlist3_name 
      logical :: read_extra_binary_job_inlist4
      character (len=256) :: extra_binary_job_inlist4_name 
      logical :: read_extra_binary_job_inlist5
      character (len=256) :: extra_binary_job_inlist5_name 

      logical :: evolve_both_stars
      logical :: relax_primary_to_th_eq
      real(dp) :: log_Lnuc_div_L_for_relax_primary_to_th_eq
      real(dp) :: min_age_for_relax_primary_to_th_eq
      integer :: max_steps_for_relax_primary_to_th_eq
      logical :: no_history_during_relax_primary_to_th_eq
      logical :: reset_age_for_relax_primary_to_th_eq
      real(dp) :: tsync_for_relax_primary_to_th_eq


