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

   INTEGER,SAVE :: fde_grid_np=0, fde_nprop=1

   integer, save :: fde_print_level = 0

   real(kind=8), save :: fde_cfg_screening = 1e-8
! io control variables
   character(len=60), save :: FDEPOT='EMBPOT'
   character(len=60), save :: FDEDNS='FRZDNS'
   character(len=60), save :: export_file='GRIDOUT'

   integer, save :: fde_log_unit     =  6
   integer, save :: fde_vemb_unit    = 45
   integer, save :: testing_unit     = 46
   integer, save :: fde_gridout_unit = 47
   integer, save :: fde_frzdns_unit  = 48

! options for importing/exporting data 
   character(len=4), save :: fde_export_level='DHF',  fde_import_level='DHF'
   character(len=3), save :: fde_export_format='XML', fde_import_format='TXT'
   
   logical, save :: do_grid_out    = .false.

   logical, save :: fde_export_n       = .false.
   logical, save :: fde_export_gn      = .false.
   logical, save :: fde_export_coulomb = .false.
   
   logical, save :: fde_import_n       = .false.
   logical, save :: fde_import_gn      = .false.
   logical, save :: fde_import_coulomb = .false.
   logical, save :: fde_import_vemb    = .false.

! options controlling other things
   logical,save :: fde_update_vemb = .false.
   logical,save :: fde_have_frozen_density = .false.

! options controllinh the calculation of the embedding energy/potential/response/...

   logical,save :: fde_skip_kin_contribution = .false.
   logical,save :: fde_skip_xc_contribution  = .false.
   logical,save :: fde_is_rsp_alda           = .false.

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! xcfun-related settings
!
      subroutine fde_set_is_alda(value)
         logical :: value
         fde_is_rsp_alda = value
      end subroutine fde_set_is_alda

!     logical function fde_is_alda
!        fde_is_alda = fde_is_rsp_alda
!     end function fde_is_alda

      subroutine fde_set_skip_xc(value)
         logical :: value
         fde_skip_kin_contribution = value
      end subroutine fde_set_skip_xc

      subroutine fde_set_skip_kin(value)
         logical :: value
         fde_skip_kin_contribution = value
      end subroutine fde_set_skip_kin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
      subroutine fde_set_print(level)
         integer :: level
         fde_printlevel = level
      end subroutine fde_set_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! import/export settings
!
      subroutine parse_export_level(unit)
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
      end subroutine parse_export_level

      subroutine set_export_file(unit)
         integer :: unit

         READ(unit,'(A60)') export_file
         DO_GRID_OUT=.TRUE.
      end subroutine set_export_file

      subroutine get_export_level(level)
         CHARACTER(LEN=4), intent(OUT) :: level
         level = fde_export_level
      end subroutine get_export_level

      subroutine get_export_format(fmat)
         CHARACTER(LEN=4), intent(OUT) :: fmat
         fmat = fde_export_format
      end subroutine get_export_format

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

      subroutine fde_set_have_frozen(value)
         logical :: value
         fde_have_frozen_density = value
      end subroutine fde_set_have_frozen

!     logical function fde_have_frozen
!        fde_have_frozen = fde_have_frozen_density
!     end function fde_have_frozen

!     logical function fde_is_vemb_import
!        fde_is_vemb_import = fde_import_vemb 
!     end function fde_is_vemb_import

!     logical function fde_is_vemb_update
!        fde_is_vemb_update = fde_update_vemb 
!     end function fde_is_vemb_update

!     logical function fde_is_export
!        fde_is_export = do_grid_out
!     end function fde_is_export
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
      subroutine fde_get_log_unit(unit)
         integer :: unit
         unit = fde_log_unit
      end subroutine fde_get_log_unit

      subroutine fde_set_log_unit(unit)
         integer :: unit
         fde_log_unit = unit
      end subroutine fde_set_log_unit


end module fde_cfg
