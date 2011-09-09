!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module fde_cfg
! this module holds all options for the embedding code, and routines to set/access
! this information

   use fde_types

   public fde_initialize_cfg
   
   public fde_set_print_level
   public fde_set_import_info
   public fde_set_export_info
   public fde_set_files_info

   public fde_get_print_level
   public fde_get_import_info
   public fde_get_export_info
   public fde_get_files_info

#ifdef FOOBAR   
   public fde_what_to_do_with_vemb
   public fde_parse_export_level
   public fde_set_frozen_import 
   public fde_set_export_file
   public fde_get_export_level
   public fde_get_export_format
   public fde_get_embpot_filename
   public fde_get_frozen_filename
   public fde_get_export_filename
   public fde_parse_vemb
   public fde_parse_frozen
   public fde_set_vemb_import
   public fde_set_vemb_update
   public fde_set_have_frozen, fde_get_have_frozen
   public fde_get_log_unit,  fde_set_log_unit
   public fde_get_embpot_unit
   public fde_get_export_unit
   public fde_get_testing_unit
   public fde_get_frozen_unit
   public fde_get_what_to_export
   public fde_get_what_to_import
   public DO_GRID_OUT
   public fde_export
#endif

   private

   integer, save :: fde_print_level = 0
   real(kind=8), save :: fde_cfg_screening = 1e-8
   ! options controlling other things
   logical,save :: fde_update_vemb = .false.
   logical,save :: fde_have_frozen_density = .false.

   type(fde_export), save :: fde_ex_control
   type(fde_import), save :: fde_im_control
   type(fde_files),  save :: fde_f


   contains
   
!   ----------------------------------------------------------------------------
      subroutine fde_set_print_level(level)
!   ----------------------------------------------------------------------------
         integer :: level
         fde_print_level = level
      end subroutine fde_set_print_level


!   ----------------------------------------------------------------------------
      subroutine fde_get_print_level(level)
!   ----------------------------------------------------------------------------
         integer :: level
         level = fde_print_level 
      end subroutine fde_get_print_level


!   ----------------------------------------------------------------------------
      subroutine fde_initialize_cfg
!   ----------------------------------------------------------------------------
         call fde_initialize_files(fde_f)
         call fde_initialize_export(fde_ex_control)
         call fde_initialize_import(fde_im_control)
      end subroutine fde_initialize_cfg
         
!   ----------------------------------------------------------------------------
      subroutine fde_initialize_files(f)
!   ----------------------------------------------------------------------------
         type(fde_files) :: f
         
         f%embpot%name = 'EMBPOT'
         f%embpot%unit = 45

         f%frozen%name = 'FRZDNS'
         f%frozen%unit = 48

         f%export%name = 'GRIDOUT'
         f%export%unit = 47

         f%test%name = 'TSTFIL'
         f%test%unit = 46

         f%logfile%name = 'DIRAC.OUT'
         f%logfile%unit = 6         
      end subroutine fde_initialize_files
      
      
!   ----------------------------------------------------------------------------
      subroutine fde_initialize_export(ec)
!   ----------------------------------------------------------------------------
         type(fde_export) :: ec
         
         ec%ex_level     = 'DHF'
         ec%ex_format    = 'XML'
         ec%ex_n         = .false.
         ec%ex_gn        = .false.
         ec%ex_coulomb   = .false.
         ec%do_grid_out  = .false.
      end subroutine fde_initialize_export


!   ----------------------------------------------------------------------------
      subroutine fde_initialize_import(ic)
!   ----------------------------------------------------------------------------
         type(fde_import) :: ic
        
         ic%im_level     = 'DHF'
         ic%im_format    = 'TXT'
         ic%im_n         = .false.
         ic%im_gn        = .false.
         ic%im_coulomb   = .false.
         ic%im_vemb      = .false.
         ic%im_frozen    = .false.
         ic%im_update_vemb=.false.
      end subroutine fde_initialize_import


!   ----------------------------------------------------------------------------
      subroutine fde_set_export_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_export), intent(in) :: tmp
         
         fde_ex_control%ex_level     = tmp%ex_level
         fde_ex_control%ex_format    = tmp%ex_format
         fde_ex_control%ex_n         = tmp%ex_n
         fde_ex_control%ex_gn        = tmp%ex_gn 
         fde_ex_control%ex_coulomb   = tmp%ex_coulomb
         fde_ex_control%do_grid_out  = tmp%do_grid_out
      end subroutine fde_set_export_info


!   ----------------------------------------------------------------------------
      subroutine fde_set_import_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_import), intent(in) :: tmp
         
         fde_im_control%im_level     = tmp%im_level
         fde_im_control%im_format    = tmp%im_format
         fde_im_control%im_n         = tmp%im_n
         fde_im_control%im_gn        = tmp%im_gn 
         fde_im_control%im_coulomb   = tmp%im_coulomb
         fde_im_control%im_vemb      = tmp%im_vemb
         fde_im_control%im_frozen    = tmp%im_frozen
         fde_im_control%im_update_vemb=tmp%im_update_vemb
      end subroutine fde_set_import_info


!   ----------------------------------------------------------------------------
      subroutine fde_set_files_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_files), intent(in) :: tmp
         
         fde_f%embpot%name = tmp%embpot%name
         fde_f%embpot%unit = tmp%embpot%unit

         fde_f%frozen%name = tmp%frozen%name
         fde_f%frozen%unit = tmp%frozen%unit

         fde_f%export%name = tmp%export%name
         fde_f%export%unit = tmp%export%unit

         fde_f%test%name = tmp%test%name
         fde_f%test%unit = tmp%test%unit 

         fde_f%logfile%name = tmp%logfile%name
         fde_f%logfile%unit = tmp%logfile%unit         
      end subroutine fde_set_files_info
      
      
!   ----------------------------------------------------------------------------
      subroutine fde_get_export_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_export), intent(out) :: tmp
         
         tmp%ex_level   = fde_ex_control%ex_level
         tmp%ex_format  = fde_ex_control%ex_format
         tmp%ex_n       = fde_ex_control%ex_n         
         tmp%ex_gn      = fde_ex_control%ex_gn           
         tmp%ex_coulomb = fde_ex_control%ex_coulomb
         tmp%do_grid_out= fde_ex_control%do_grid_out    
      end subroutine fde_get_export_info


!   ----------------------------------------------------------------------------
      subroutine fde_get_import_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_import), intent(out) :: tmp
         
         tmp%im_level       = fde_im_control%im_level       
         tmp%im_format      = fde_im_control%im_format      
         tmp%im_n           = fde_im_control%im_n           
         tmp%im_gn          = fde_im_control%im_gn           
         tmp%im_coulomb     = fde_im_control%im_coulomb     
         tmp%im_vemb        = fde_im_control%im_vemb        
         tmp%im_frozen      = fde_im_control%im_frozen      
         tmp%im_update_vemb = fde_im_control%im_update_vemb 
      end subroutine fde_get_import_info

      
!   ----------------------------------------------------------------------------
      subroutine fde_get_files_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_files), intent(out) :: tmp
         
         tmp%embpot%name = fde_f%embpot%name
         tmp%embpot%unit = fde_f%embpot%unit

         tmp%frozen%name = fde_f%frozen%name
         tmp%frozen%unit = fde_f%frozen%unit

         tmp%export%name = fde_f%export%name
         tmp%export%unit = fde_f%export%unit

         tmp%test%name = fde_f%test%name
         tmp%test%unit = fde_f%test%unit 

         tmp%logfile%name = fde_f%logfile%name
         tmp%logfile%unit = fde_f%logfile%unit         
      end subroutine fde_get_files_info
      
      








#ifdef FOOBAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! import/export settings
!
      subroutine fde_parse_export_level(unit)
         integer :: unit

         READ(unit,'(A4)') fde_export_level
         select case(fde_export_level)
         case('DHF','MP2')
            WRITE(fde_log_unit,*)                                  &
            'Outputted density will be taken from: ',trim(fde_export_level)
         case default
            WRITE(fde_log_unit,*) 'Input for outputted density not a &
     &  x  valid calculation, HF will be used instead. Given:',fde_export_level
         end select
      end subroutine fde_parse_export_level


      subroutine fde_set_export_file(unit)
         integer :: unit

         READ(unit,'(A60)') export_file
         DO_GRID_OUT=.TRUE.
      end subroutine fde_set_export_file

      subroutine fde_get_export_level(level)
         CHARACTER(LEN=4), intent(OUT) :: level
         level = fde_export_level
      end subroutine fde_get_export_level

      subroutine fde_get_export_format(fmat)
         CHARACTER(LEN=4), intent(OUT) :: fmat
         fmat = fde_export_format
      end subroutine fde_get_export_format


      subroutine fde_get_embpot_filename(string)
         character(len=60), intent(out) :: string
         string = FDEPOT
      end subroutine fde_get_embpot_filename
         
      subroutine fde_get_frozen_filename(string)
         character(len=60), intent(out) :: string
         string = FDEDNS
      end subroutine

      subroutine fde_get_export_filename(string)
         character(len=60), intent(out) :: string
         string = export_file 
      end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
      subroutine fde_parse_vemb(unit)
         integer :: unit
         read(unit,'(A60)') FDEPOT
         fde_import_vemb = .true.
         fde_update_vemb = .false.
      end subroutine fde_parse_vemb

      subroutine fde_parse_frozen(lucmd)
         integer :: unit
         read(unit,'(A60)') FDEDNS
         fde_have_frozen_density = .true.
      end subroutine fde_parse_frozen

      subroutine fde_set_vemb_import(value)
         logical :: value
         fde_import_vemb = value
      end subroutine fde_set_vemb_import

      subroutine fde_set_vemb_update(value)
         logical :: value
         fde_update_vemb = value
      end subroutine fde_set_vemb_update

      subroutine fde_set_frozen_import(value)
         logical :: value
         fde_have_frozen_density = value
      end subroutine fde_set_frozen_import


      subroutine fde_set_have_frozen(value)
         logical :: value
         fde_have_frozen_density = value
      end subroutine fde_set_have_frozen

      subroutine fde_get_have_frozen(value)
         logical, intent(out) :: value
         value = fde_have_frozen_density 
      end subroutine fde_get_have_frozen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine fde_get_log_unit(unit)
         integer :: unit
         unit = fde_log_unit
      end subroutine fde_get_log_unit

      subroutine fde_set_log_unit(unit)
         integer :: unit
         fde_log_unit = unit
      end subroutine fde_set_log_unit


      subroutine fde_get_embpot_unit(unit)
         integer, intent(out) :: unit
         unit = fde_vemb_unit
      end subroutine fde_get_embpot_unit


      subroutine fde_get_export_unit(unit)
         integer, intent(out) :: unit
         unit = fde_gridout_unit
      end subroutine fde_get_export_unit


      subroutine fde_get_testing_unit(unit)
         integer, intent(out) :: unit
         unit = testing_unit
      end subroutine fde_get_testing_unit


      subroutine fde_get_frozen_unit(unit)
         integer, intent(out) :: unit
         unit = fde_frzdns_unit
      end subroutine fde_get_frozen_unit


      subroutine fde_get_what_to_export(n,gn,vc)
         logical, intent(out) :: n, gn, vc

         vc = fde_export_coulomb  
         gn = fde_export_gn
         n  = fde_export_n
      end  subroutine fde_get_what_to_export

      subroutine fde_get_what_to_import(n,gn,vc,v)
         logical, intent(out) :: n, gn, vc,v

         vc = fde_import_coulomb  
         gn = fde_import_gn
         n  = fde_import_n
         v  = fde_import_vemb
      end  subroutine fde_get_what_to_import

      subroutine fde_what_to_do_with_vemb(import,update)  
         logical, intent(out) :: update, import
         update = fde_update_vemb
         import = fde_import_vemb
      end subroutine fde_what_to_do_with_vemb
#endif

end module fde_cfg
