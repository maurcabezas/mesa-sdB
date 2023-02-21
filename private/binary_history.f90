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


      module binary_history

      use const_def
      use chem_def
      use star_lib
      use star_def
      use crlibm_lib
      use binary_def
      use binary_private_def
      use binary_history_specs

      implicit none

      contains

      integer function how_many_binary_history_columns(binary_id)
         integer, intent(in) :: binary_id
         integer :: numcols, ierr
         type (binary_info), pointer :: b

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            numcols = 0
            return
         end if

         if (.not. associated(b% history_column_spec)) then
            numcols = 0
         else
            numcols = size(b% history_column_spec, dim=1)
         end if

         how_many_binary_history_columns = numcols
      end function how_many_binary_history_columns
   
   
      subroutine data_for_binary_history_columns( &
            binary_id, n, names, vals, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id, n
         character (len=80) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr

         type (binary_info), pointer :: b
         integer :: c, int_val, i ,j
         logical :: is_int_val
         real(dp) :: val

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         do j=1,n
            c = b% history_column_spec(j)
            if (c > lg_mdot_system_offset) then
               i = c - lg_mdot_system_offset
               names(j) = 'lg_mdot_system_' // trim(chem_isos% name(i))
            else
               names(j) = trim(binary_history_column_name(c))
            end if
            call binary_history_getval( &
               b, c, val, int_val, is_int_val, ierr)
            if (ierr /= 0) then
               write(*,*) "Unknown binary_history_columns.list column"
               return
            end if
            if (is_int_val) then
               vals(j) = int_val
            else
               vals(j) = val
            end if
         end do
      end subroutine data_for_binary_history_columns

      subroutine write_binary_history_info ( b, ierr)
         use utils_lib, only:alloc_iounit, free_iounit
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr
         character (len=maxlen_profile_column_name), pointer :: names(:) ! (num_history_columns)
         real(dp), pointer :: vals(:) ! (num_history_columns)
         logical, pointer :: is_int(:)
         logical, parameter :: write_flag = .true.
         names => null()
         vals => null()
         is_int => null()
         call do_binary_history_info( &
            b, &
            write_flag, names, vals, is_int, ierr)
      end subroutine write_binary_history_info


      subroutine do_binary_history_info( &
            b, &
            write_flag, names, vals, is_int, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         type (binary_info), pointer :: b
         logical, intent(in) :: write_flag
         character (len=maxlen_profile_column_name), pointer :: names(:) ! (num_history_columns)
         real(dp), pointer :: vals(:) ! (num_history_columns)
         logical, pointer :: is_int(:)
         integer, intent(out) :: ierr
         
         character (len=strlen) :: fname, dbl_fmt, int_fmt, txt_fmt
         integer :: numcols, io, i, nz, col, j, i0, num_extra_cols
         character (len=maxlen_binary_history_column_name), pointer :: extra_col_names(:) => null()
         real(dp), pointer :: extra_col_vals(:) => null()
        
         include 'formats'
         
         dbl_fmt = b% history_dbl_format
         int_fmt = b% history_int_format
         txt_fmt = b% history_txt_format
         
         ierr = 0
         
         if (write_flag) then
            io = alloc_iounit(ierr)
            if (ierr /= 0) return
         else
            io = -1 ! set to invalid value to trigger complaint if use it by mistake
         end if
         
         if (.not. associated(b% history_column_spec)) then
            numcols = 0
         else
            numcols = size(b% history_column_spec, dim=1)
         end if
         
         if (numcols == 0) then
            write(*,*) 'WARNING: do not have any output specified for binary logs.'
            if (io > 0) call free_iounit(io)
            return
         end if
         
         num_extra_cols = b% how_many_extra_binary_history_columns(b% binary_id)
         if (num_extra_cols > 0) then
            allocate( &
               extra_col_names(num_extra_cols), extra_col_vals(num_extra_cols), stat=ierr)
            if (ierr /= 0) then
               if (io > 0) call free_iounit(io)
               return
            end if
            call b% data_for_extra_binary_history_columns( &
               b% binary_id, num_extra_cols, extra_col_names, extra_col_vals, ierr)
            if (ierr /= 0) then
               deallocate(extra_col_names, extra_col_vals)
               if (io > 0) call free_iounit(io)
               return
            end if
         end if
         
         i0 = 1
         if (write_flag .and. (open_close_log .or. b% s_donor% model_number == -100)) then
            fname = trim(b% log_directory) // '/' // trim(b% history_name)
            if ((.not. history_file_exists(fname,io)) .or. b% open_new_history_file) then
               ierr = 0
               open(unit=io, file=trim(fname), action='write', iostat=ierr)
               b% open_new_history_file = .false.
            else
               i0 = 3            
               open(unit=io, file=trim(fname), action='write', position='append', iostat=ierr)
            end if
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(fname)
               if (io > 0) call free_iounit(io)
               return
            end if
         end if
         
         if (write_flag .and. i0 == 1) then ! write parameters at start of log
            !call b% other_binary_history_data_initialize(b, ierr)
            !if (ierr /= 0) return
            do i=1,3
               col = 0
               call write_integer(io, col, i, 'version_number', version_number)
               call write_val(io, col, i, 'initial_don_mass', initial_mass(1))
               call write_val(io, col, i, 'initial_acc_mass', initial_mass(2))
               call write_val(io, col, i, 'initial_period_days', &
                   initial_binary_period/(3600*24))
               write(io,*)
            end do
            write(io,*)
         end if

         do i=i0,3 ! add a row to the log
            !call b% other_binary_history_data_add_model(s% id, ierr)
            !if (ierr /= 0) return
            col = 0
            if (i==3) then
            
               !if (write_flag .and. i0 == 1) then
                  !close(io)
                  !stop "enough"
                  !fname = trim(b% log_directory) // '/' // trim(b% history_name)
                  !open(unit=io, file=trim(fname), action='write',status='replace', iostat=ierr)
                  !if (ierr /= 0) then
                  !   call free_iounit(io); return
                  !end if
               !end if
               
            end if
            do j=1,numcols
               call do_col(i, j)
            end do
            do j=1,num_extra_cols
               call do_extra_col(i, j)
            end do
            if (write_flag) write(io, *)
         end do
         
         if (open_close_log) close(io)
         
         if (io > 0) call free_iounit(io)

         if (associated(extra_col_names)) deallocate(extra_col_names)
         if (associated(extra_col_vals)) deallocate(extra_col_vals)
         
         
         contains


         subroutine do_extra_col(pass, j)
            integer, intent(in) :: pass, j
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') j + numcols
            else if (pass == 2) then
               call do_name(j + numcols, extra_col_names(j))
            else if (pass == 3) then
               call do_val(j + numcols, extra_col_vals(j))
            end if
         end subroutine do_extra_col


         subroutine do_name(j, col_name)
            integer, intent(in) :: j
            character (len=*), intent(in) :: col_name
            if (write_flag) then
               write(io, fmt=txt_fmt, advance='no') trim(col_name)
            else
               names(j) = trim(col_name)
            end if
         end subroutine do_name
         

         subroutine do_col(pass, j)
            integer, intent(in) :: pass, j
            if (pass == 1) then
               call do_col_pass1
            else if (pass == 2) then
               call do_col_pass2(j)
            else if (pass == 3) then
               call do_col_pass3(b% history_column_spec(j))
            end if
         end subroutine do_col
         
         
         subroutine do_col_pass1 ! write the column number
            col = col+1
            if (write_flag) write(io, fmt=int_fmt, advance='no') col
         end subroutine do_col_pass1
         
         
         subroutine do_col_pass2(j) ! get the column name
            integer, intent(in) :: j
            character (len=100) :: col_name
            character (len=10) :: str
            integer :: c, i, ii
            c = b% history_column_spec(j)
            if (c > lg_mdot_system_offset) then
               i = c - lg_mdot_system_offset
               col_name = 'lg_mdot_system_' // trim(chem_isos% name(i))
            else
               col_name = trim(binary_history_column_name(c))
            end if
            call do_name(j, col_name)
         end subroutine do_col_pass2
         
         
         subroutine do_col_pass3(c) ! get the column value
            integer, intent(in) :: c
            integer :: i, ii, k, int_val
            logical :: is_int_val
            real(dp) :: val, val1, Ledd, power_photo, frac
            int_val = 0; val = 0; is_int_val = .false.
            call binary_history_getval( &
               b, c, val, int_val, is_int_val, ierr)
            if (ierr /= 0) then
               write(*,*) 'missing log info for ' // trim(binary_history_column_name(c)), j, k
               return
            end if
            if (is_int_val) then
               call do_int_val(j,int_val)
            else
               call do_val(j,val)
            end if
         end subroutine do_col_pass3
         
         
         subroutine do_val(j, val)
            use utils_lib, only: is_bad
            integer, intent(in) :: j
            real(dp), intent(in) :: val
            if (write_flag) then
               if (is_bad(val)) then
                  write(io, fmt=dbl_fmt, advance='no') -1d99
               else
                  write(io, fmt=dbl_fmt, advance='no') val
               end if
            else
               vals(j) = val
               is_int(j) = .false.
            end if
         end subroutine do_val
         
         
         subroutine do_int_val(j, val)
            integer, intent(in) :: j
            integer, intent(in) :: val
            if (write_flag) then
               write(io, fmt=int_fmt, advance='no') val
            else
               vals(j) = dble(val)
               is_int(j) = .true.
            end if
         end subroutine do_int_val
                  
         
      end subroutine do_binary_history_info
      


      subroutine write_integer(io, col, pass, name, val)
         integer, intent(in) :: io, pass
         integer, intent(inout) :: col
         character (len=*), intent(in) :: name
         integer, intent(in) :: val
         if (pass == 1) then
            col = col+1
            write(io, fmt='(i28, 1x)', advance='no') col
         else if (pass == 2) then
            write(io, fmt='(a28, 1x)', advance='no') trim(name)
         else if (pass == 3) then
            write(io, fmt='(i28, 1x)', advance='no') val
         end if
      end subroutine write_integer
      
      
      subroutine write_val(io, col, pass, name, val)
         integer, intent(in) :: io, pass
         integer, intent(inout) :: col
         character (len=*), intent(in) :: name
         real(dp), intent(in) :: val
         if (pass == 1) then
            col = col+1
            write(io, fmt='(i28, 1x)', advance='no') col
         else if (pass == 2) then
            write(io, fmt='(a28, 1x)', advance='no') trim(name)
         else if (pass == 3) then
            write(io, fmt='(1pe28.16e3, 1x)', advance='no') val
         end if
      end subroutine write_val
      
      
      subroutine binary_history_getval( &
            b, c, val, int_val, is_int_val, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: c
         real(dp), intent(out) :: val
         integer, intent(out) :: int_val
         logical, intent(out) :: is_int_val
         integer, intent(out) :: ierr
         integer :: k, i
         
         include 'formats'
         
         ierr = 0
         is_int_val = .false.
         int_val = 0
         val = 0

         if (c > lg_mdot_system_offset) then
            i = c - lg_mdot_system_offset
            k = b% s_donor% net_iso(i)
            val = 0
            if (k > 0) then
               if (b% s_donor% mstar_dot < 0) then
                  val = b% s_donor% xa_removed(k) * b% s_donor% mstar_dot
               end if
            end if
            if (b% point_mass_i == 0) then
               if (b% s_accretor% mstar_dot > 0) then
                  val = val + b% s_donor% xa_removed(k) * b% s_accretor% mstar_dot
               else
                   k = b% s_accretor% net_iso(i)
                   if (k > 0) &
                       val = val + b% s_accretor% xa_removed(k) * b% s_accretor% mstar_dot
               end if
            else
               val = val - b% mtransfer_rate * b% xfer_fraction
            end if
            val = safe_log10_cr(-val / Msun * secyer)
         else
               
            select case(c)
            
            case(bh_model_number)
               int_val = b% model_number
               is_int_val = .true.
            case(bh_age)
               val = b% binary_age
            case(bh_donor_index)
               int_val = b% d_i
               is_int_val = .true.
            case(bh_period_days)
               val = b% period/(60d0*60d0*24d0)
            case(bh_period_hr)
               val = b% period/(60d0*60d0)
            case(bh_period_minutes)
               val = b% period/60d0
            case(bh_lg_separation)
               val = safe_log10_cr(b% separation)
            case(bh_binary_separation)
               val = b% separation/Rsun
            case(bh_eccentricity)
               val = b% eccentricity
            case(bh_star_1_radius)
               val = b% r(1)/Rsun
            case(bh_star_2_radius)
               val = b% r(2)/Rsun
            case(bh_rl_1)
               val = b% rl(1)/Rsun
            case(bh_rl_2)
               val = b% rl(2)/Rsun
            case(bh_rl3_1)
               val = b% rl3(1)/Rsun
            case(bh_rl3_2)
               val = b% rl3(2)/Rsun  
            case(bh_rl_overflow_1)
               val = (b% r(1)-b% rl(1))/Rsun
            case(bh_rl_overflow_2)
               val = (b% r(2)-b% rl(2))/Rsun
            case(bh_rl_relative_overflow_1)
               val = b% rl_relative_gap(1)
            case(bh_rl_relative_overflow_2)
               val = b% rl_relative_gap(2)
            case(bh_P_rot_div_P_orb_1)
               if (b% point_mass_i /= 1) then
                  val = 2 * pi / b% s1% omega_avg_surf / b% period
               else
                  val = 0.0d0
               end if
            case(bh_P_rot_div_P_orb_2)
               if (b% point_mass_i /= 2) then
                  val = 2 * pi / b% s2% omega_avg_surf / b% period
               else
                  val = 0.0d0
               end if
            case(bh_lg_t_sync_1)
               val = safe_log10_cr(abs(b% t_sync_1)/secyer)
            case(bh_lg_t_sync_2)
               val = safe_log10_cr(abs(b% t_sync_2)/secyer)
            case(bh_star_1_mass)
               val = b% m(1)/Msun
            case(bh_lg_star_1_mass)
               val = safe_log10_cr(b% m(1)/Msun)
            case(bh_star_2_mass)
               val = b% m(2)/Msun
            case(bh_lg_star_2_mass)
               val = safe_log10_cr(b% m(2)/Msun)
            case(bh_sum_of_masses)
               val = (b% m(1) + b% m(2))/Msun
            case(bh_lg_mtransfer_rate)
               val = safe_log10_cr(abs(b% mtransfer_rate)/Msun*secyer)
            case(bh_lg_mstar_dot_1)
               if (b% point_mass_i /= 1) then
                  val = safe_log10_cr(abs(b% s1% mstar_dot)/Msun*secyer)
               else
                  val = safe_log10_cr(abs(b% step_mtransfer_rate)*b% xfer_fraction/Msun*secyer)
               end if
            case(bh_lg_mstar_dot_2)
               if (b% point_mass_i /= 2) then
                  val = safe_log10_cr(abs(b% s2% mstar_dot)/Msun*secyer)
               else
                  val = safe_log10_cr(abs(b% step_mtransfer_rate)*b% xfer_fraction/Msun*secyer)
               end if
            case(bh_lg_system_mdot_1)
               val = safe_log10_cr(abs(b% mdot_system_transfer(1))/Msun*secyer)
            case(bh_lg_system_mdot_2)
               val = safe_log10_cr(abs(b% mdot_system_transfer(2))/Msun*secyer)
            case(bh_lg_wind_mdot_1)
               val = safe_log10_cr(abs(b% mdot_system_wind(1))/Msun*secyer)
            case(bh_lg_wind_mdot_2)
               val = safe_log10_cr(abs(b% mdot_system_wind(2))/Msun*secyer)
            case(bh_lg_mdot_thin)
               val = safe_log10_cr(abs(b% mdot_thin)/Msun*secyer)
            case(bh_lg_mdot_thin_L3)
               val = safe_log10_cr(abs(b% mdot_thin_L3)/Msun*secyer)
            case(bh_lg_mdot_thick)
               val = safe_log10_cr(abs(b% mdot_thick)/Msun*secyer)
            case(bh_lg_mdot_thick_L3)
               val = safe_log10_cr(abs(b% mdot_thick_L3)/Msun*secyer)
               
            case(bh_mass_transfer_alpha)
               val = b% mass_transfer_alpha
            case(bh_mass_transfer_beta)
               val = b% mass_transfer_beta
            case(bh_mass_transfer_delta)
               val = b% mass_transfer_delta
            case(bh_mass_transfer_gamma)
               val = b% mass_transfer_gamma
               
            case(bh_star_1_div_star_2_mass)
               val = b% m(1) / b% m(2)
            case(bh_delta_star_1_mass)
               val = b% m(1) - initial_mass(1)
            case(bh_delta_star_2_mass)
               val = b% m(2) - initial_mass(2)
            case(bh_lg_F_irr)
               val = safe_log10_cr(b% s_donor% irradiation_flux)
            case(bh_xfer_fraction)
               val = b% xfer_fraction
            case(bh_lg_mdot_edd)
               if (b% limit_retention_by_mdot_edd) then
                  val = safe_log10_cr(b% mdot_edd/Msun*secyer)
               else
                  val = safe_log10_cr(0d0)
               end if
            case(bh_mdot_edd_eta)
               if (b% limit_retention_by_mdot_edd) then
                  val = b% mdot_edd_eta
               else
                  val = 0d0
               end if
            case(bh_lg_accretion_luminosity)
               if (b% limit_retention_by_mdot_edd) then
                  val = safe_log10_cr (b% mdot_edd_accretion_luminosity / Lsun) 
               else
                  val = safe_log10_cr(0d0)
               end if
            case(bh_bh_spin)
               if (b% point_mass_i /= 0) then
                  val = sqrt(two_thirds) &
                     *(b% eq_initial_bh_mass/min(b% m(b% point_mass_i),sqrt(6d0)*b% eq_initial_bh_mass)) &
                     *(4 - sqrt(18*(b% eq_initial_bh_mass/min(b% m(b% point_mass_i),sqrt(6d0)*b% eq_initial_bh_mass))**2 - 2))
               else
                  val = 0
               end if
            case(bh_v_orb_1)
               val = 2.0d0 * pi * b% m(2)/(b% m(1) + b% m(2)) * b% separation / b% period / 1.0d5
            case(bh_v_orb_2)
               val = 2.0d0 * pi * b% m(1)/(b% m(1) + b% m(2)) * b% separation / b% period / 1.0d5
            case(bh_J_orb)
               val = b% angular_momentum_j
            case(bh_J_spin_1)
               if (b% point_mass_i /= 1) then
                  val = b% s1% total_angular_momentum
               else
                  val = 0d0
               end if
            case(bh_J_spin_2)
               if (b% point_mass_i /= 2) then
                  val = b% s2% total_angular_momentum
               else
                  val = 0d0
               end if
            case(bh_J_total)
               val = b% angular_momentum_j
               if (b% point_mass_i /= 1) &
                  val = val + b% s1% total_angular_momentum
               if (b% point_mass_i /= 2) &
                  val = val + b% s2% total_angular_momentum
               val = val
            case(bh_Jdot)
               val = b% jdot
            case(bh_jdot_mb)
               val = b% jdot_mb
            case(bh_jdot_gr)
               val = b% jdot_gr
            case(bh_jdot_ml)
               val = b% jdot_ml
            case(bh_jdot_ls)
               val = b% jdot_ls
            case(bh_jdot_missing_wind)
               val = b% jdot_missing_wind
            case(bh_extra_jdot)
               val = b% extra_jdot
            case(bh_accretion_mode)
               int_val = b% accretion_mode
               is_int_val = .true.
            case(bh_acc_am_div_kep_am)
               val = b% acc_am_div_kep_am
            case(bh_edot)
               val = b% edot
            case(bh_edot_tidal)
               val = b% edot_tidal
            case(bh_edot_enhance)
               val = b% edot_enhance
            case(bh_extra_edot)
               val = b% extra_edot
            case(bh_point_mass_index)
               is_int_val = .true.
               int_val = b% point_mass_i
            
            case default
               ierr = -1
            
            end select
         end if
         
                  
      end subroutine binary_history_getval

      end module binary_history
