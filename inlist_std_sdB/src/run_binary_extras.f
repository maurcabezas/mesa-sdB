! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use crlibm_lib
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pinters to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% other_jdot_ml => jdot_ml_new
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls
      
! routine other_jdot_ml ----------------------new Jdot(e)

      subroutine jdot_ml_new(binary_id, ierr)
         implicit none
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: jdot_alpha, jdot_beta_new, jdot_delta
         real(dp) :: ml1,msf,mcf,lf,f3,f1,f2
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
         jdot_alpha = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             (b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period *&
             sqrt(1 - b% eccentricity**2)

         !mass lost from circumbinary coplanar toroid
         jdot_delta = b% mdot_system_cct * b% mass_transfer_gamma * &
             sqrt(b% s_donor% cgrav(1) * (b% m(1) + b% m(2)) * b% separation)

         ml1=(b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i)) !zero at begginng 
         msf=(b% m(b% d_i)/(b% m(b% d_i)+b% m(b% a_i)))**2 !MS factor 
         mcf=(b% m(b% d_i)/(b% m(b% d_i)+b% m(1)))**2 !Mcomp factor
         lf=(2*pi/b% period) * (b% separation)**2 * sqrt(1 - b% eccentricity**2)
         f3=(0.500-0.227*log10(b% m(b% d_i)/b% m(b% a_i)))**2
         
         !jdot_beta_min
         f1 = ml1 * msf * lf
         !jdot_beta_max
         f2 = ml1*(mcf+f3)*lf
         !jdot beta
         jdot_beta_new=((1- b% eps_deg_mix)*f1)+ (b% eps_deg_mix*f2)
         b% extra_jdot =jdot_alpha+jdot_beta_new+jdot_delta

      end subroutine jdot_ml_new     
      
!---------------------------------------------     
      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 2
      end function how_many_extra_binary_history_columns
 !---------------------------------------------        
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n),jdot_alpha, jdot_beta, jdot_delta,ml1,msf,mcf,f3,lf,f1,f2
         integer, intent(out) :: ierr
         !real(dp) :: beta
         ! integer :: k
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         names(1) = 'jdot_beta_new'
         names(2) = 'jdot_ml_new'

         ml1=(b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i)) !zero at begginng 
         msf=(b% m(b% d_i)/(b% m(b% d_i)+b% m(b% a_i)))**2 !MS factor 
         mcf=(b% m(b% d_i)/(b% m(b% d_i)+b% m(1)))**2 !Mcomp factor
         lf=(2*pi/b% period) * (b% separation)**2 * sqrt(1 - b% eccentricity**2)
         f3=(0.500-0.227*log10(b% m(b% d_i)/b% m(b% a_i)))**2

         !jdot_beta_min
         f1 = ml1 * msf * lf
         !jdot_beta_max
         f2 = ml1*(mcf+f3)*lf
         !mass lost from circumbinary coplanar toroid
         vals(1) = f1
         vals(2) = f2
         
      end subroutine data_for_extra_binary_history_columns
      
      
 !---------------------------------------------        
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
          extras_binary_startup = keep_going
      end function  extras_binary_startup
 !---------------------------------------------        
      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
 !---------------------------------------------        
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
 !---------------------------------------------        
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
