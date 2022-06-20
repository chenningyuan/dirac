!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module fde_input

   use fde_types
   use fde_cfg
   use fde_data
   use fde_nadd_derv
   use fde_evaluators_dirac
   
   implicit none

   public fde_input_init
   public fde_input_validate_and_postprocess
   public fde_input_show_setup

   public skip_all_imports

   
   logical, save :: skip_all_imports = .false.

   contains

!   ----------------------------------------------------------------------------
      subroutine fde_input_init(log_unit,log_name)
!   ----------------------------------------------------------------------------
         integer          :: log_unit
         character(len=*) :: log_name
         type(fde_files)  :: ftmp

         call fde_initialize_cfg
         
         call fde_get_files_info(ftmp)
         ftmp%logfile%unit = log_unit
         ftmp%logfile%name = log_name
         call fde_set_files_info(ftmp)
         
         call fde_qccode_data_interface
      end subroutine fde_input_init

      
!   ----------------------------------------------------------------------------
      subroutine fde_input_validate_and_postprocess
!   ----------------------------------------------------------------------------
         type(fde_export)  :: etmp
         type(fde_import)  :: itmp
         type(fde_files)   :: ftmp
         integer :: pl
         integer :: lunit

      call fde_get_files_info(ftmp)
      call fde_get_export_info(etmp)
      call fde_get_import_info(itmp)

      lunit = ftmp%logfile%unit

      if (.not.skip_all_imports) then

      if (itmp%im_vemb) then
         if (itmp%im_update_vemb) then
            write (lunit,*) '.UPDATE and .EMBPOT cannot be used together!'
            call fde_quit('conflicting keywords selected')
         endif
      else
         itmp%im_update_vemb = .true.
         itmp%im_frozen      = .true.
      endif
 
! we should arrive here if the update of the embedding potential during the
! scf is asked for, or if fde response        
      if (itmp%im_frozen) then
         itmp%im_n        = .true.
         itmp%im_gn       = .true.
         itmp%im_coulomb  = .true.
      end if
      if (itmp%im_frozen_pert_direct .and. .not. fde_cfg_no_sdft) then
         itmp%im_s_b_direct      = .true.
         itmp%im_gs_b_direct     = .true.
      endif
      if (itmp%im_frozen_pert_reorth) then
         itmp%im_n_b_reorth      = .true.
         itmp%im_gn_b_reorth     = .true.
         if (.not. fde_cfg_no_sdft) then
            itmp%im_s_b_reorth      = .true.
            itmp%im_gs_b_reorth     = .true.
         end if
      endif
!g      write(*, *) 'test1 itmp%im_n_b =', itmp%im_n_b 

      endif ! .not.skip_all_imports

      if (etmp%do_grid_out) then
         etmp%ex_n        = .true.
         etmp%ex_coulomb  = .true.
         etmp%ex_gn       = .true.
!gosia: not tested:
!         if (fde_rsp_mag_lao_export) then
!            etmp%ex_n_b      = .true.
!            etmp%ex_gn_b     = .true.
!            if (.not. fde_cfg_no_sdft) then
!               etmp%ex_s_b      = .true.
!               etmp%ex_gs_b     = .true.
!            end if
!         end if
      endif

      call fde_initialize_nadd_functionals         

      if (itmp%im_vemb) then
         call fde_import_static      
      endif

      if (itmp%im_frozen .or. itmp%im_frozen_pert_direct .or. itmp%im_frozen_pert_reorth) then
         call fde_import_frozen
      endif

      if (etmp%do_grid_out) then
         call fde_import_gridout
      endif
      
! saving information on 
      call fde_set_files_info(ftmp)
      call fde_set_export_info(etmp)
      call fde_set_import_info(itmp)

     end subroutine  


      subroutine fde_input_show_setup
         integer :: unit
         character(len=60) :: string
         type(fde_export)  :: etmp
         type(fde_import)  :: itmp
         type(fde_files)   :: ftmp
         integer :: pl

         call fde_get_print_level(pl)

         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
         call fde_get_import_info(itmp)

         unit = ftmp%logfile%unit
         
         write(unit,*)
         write(unit,'(A,I5)')  ' * Print level                  : ',pl

         if (itmp%im_frozen) call fde_print_nadd_functionals(unit)

         if (itmp%im_vemb) then
            string = ftmp%embpot%name
            
            write(unit,'(A,A60)')  ' * FDE Potential read from file : ',string
            write(unit,*)
            write(unit,'(3X,2A,/3X,2A)') &
            'Enviroment effects included via the ',          &
            'fixed potential method described in:',          &
            'A.S.P. Gomes, C. R. Jacob and L. Visscher, ', &
            'PCCP 10 (2008) 5353-5362.'
            write(unit,*)
            if (itmp%im_remove_vemb_e_dcount) then
               write(unit,'(3X,3A)') &
               'Double-counting of fixed potential will be ',  &
               'removed from active subsystem energy, according to ',  &
               'usual convention in FDE/subsystem DFT theory.'
            else
               write(unit,'(3X,4A)') &
               'Double-counting of fixed potential will NOT be', &
               'removed from active subsystem energy. This is not consistent ', &
               'with FDE/subsystem DFT theory conventions, but do not affect ', &
               'orbitals or orbital energies.' 
            end if
            write(unit,*)
         else
            write(unit,'(A,A60)')  ' * FDE Potential generated from frozen, active densities'
            write(unit,*)
            write(unit,'(3X,2A,/3X,2A,/3X,2A)') &
            'Embedding potential included according to:',    &
            'S. Hofener, A.S.P. Gomes and L. Visscher,',&
            ' J. Chem. Phys. 136 (2012) 044104.',       &
            'S. Hofener, A.S.P. Gomes and L. Visscher,',&
            ' J. Chem. Phys. 139 (2013) 104106.' 
            write(unit,*)
            call fde_test_frozen_density
         endif

         if (itmp%im_frozen) then
            string = ftmp%frozen%name
            write(unit,'(A,A60)') &
                ' * Density (and gradient) from frozen subsystems read from file: ',string
            write(unit,*)
            write(unit,'(3X,2A,/3X,2A,/3X,2A)') &
            'Response contributions to the active subsystem,', &
            'when included, are calculated according to:', &
            'S. Hofener, A.S.P. Gomes and L. Visscher,',&
            ' J. Chem. Phys. 136 (2012) 044104.',       &
            'S. Hofener, A.S.P. Gomes and L. Visscher,',&
            ' J. Chem. Phys. 139 (2013) 104106.' 
            write(unit,*)
         endif

         if (itmp%im_frozen_pert_direct) then
            string = ftmp%frozen_pert_direct%name
            write(unit,'(A,A60)') &
                ' * Direct-LAO part of perturbed density of frozen subsystems read from file: ',string
            write(unit,*)
         endif

         if (itmp%im_frozen_pert_reorth) then
            string = ftmp%frozen_pert_reorth%name
            write(unit,'(A,A60)') &
                ' * Reorthonormalization part of perturbed density of frozen subsystems read from file: ',string
            write(unit,*)
         endif

         if (etmp%do_grid_out) then
            string = ftmp%export%name
            write(unit,'(2A)') &
            ' * FDE Gridfile with updated density written to:',string

            select case(etmp%ex_level)
               case('DHF','MP2','CCSD')
                  write(unit,*)                               &
                  'Outputted density will be taken from: ',trim(etmp%ex_level)
                  
               case default
                  write(unit,*) 'Input for outputted density not a &
     &  x  valid calculation, HF will be used instead. Given:',trim(etmp%ex_level)
                  etmp%ex_level = 'DHF'
                  
            end select
            if (fde_rsp_mag_lao_export) then
               write(unit, *) 'magnetic derivatives of density and density gradient also written to ', string
            end if
         endif
                  
      end subroutine fde_input_show_setup

end module fde_input
