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

module fde_inputDirac

   use fde_types
   use fde_cfg
   use fde_data
   use fde_nadd_derv
   use fde_eval_dirac
   
   private

   public fdeini
   public fdeinp
   
   logical, save :: skip_all_imports = .false.

   contains

!   ----------------------------------------------------------------------------
      subroutine fdeini(log_unit,log_name)
!   ----------------------------------------------------------------------------
         integer          :: log_unit
         character(len=*) :: log_name
         type(fde_files)  :: ftmp

         call fde_initialize_cfg
         
         call fde_get_files_info(ftmp)
         ftmp%logfile%unit = log_unit
         ftmp%logfile%name = log_name
         call fde_set_files_info(ftmp)
         
         call dirac_data_interface
      end subroutine fdeini


!   ----------------------------------------------------------------------------
      subroutine fdeinp(word,echo_input)
!   ----------------------------------------------------------------------------
         use dirac_cfg
#include "implicit.h"
#include "priunit.h" 
#include "dcbgen.h"
      PARAMETER (NTABLE = 8)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, line*80
 
      logical :: Echo_Input,newdef
      integer :: pl
      integer :: lunit
      type(fde_files)  :: ftmp
      type(fde_export) :: etmp
      type(fde_import) :: itmp
      character(len=80):: nadd_fun
!
      call fde_get_files_info(ftmp)
      call fde_get_export_info(etmp)
      call fde_get_import_info(itmp)
      
      lunit = ftmp%logfile%unit
      
      ICHANG = 0
      IF (WORD(1:4).eq.'*FDE') THEN
        WORD1 = WORD
        DO
          READ (LUCMD, '(A7)') WORD
          SELECT CASE (WORD(1:1))
          CASE('!','#')
            CYCLE
          CASE('.')
            ICHANG = ICHANG + 1
            SELECT CASE(WORD(1:7))
            CASE ('.OPTION')
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',lunit)
              CYCLE

            CASE ('.EMBPOT')
               read(lucmd,'(A60)') ftmp%embpot%name
               itmp%im_vemb        = .true.
               itmp%im_update_vemb = .false.
               
            CASE ('.PRINT')
              read(lucmd,*) pl
              call fde_set_print_level(pl)

#ifdef MOD_UNRELEASED
            CASE ('.FRDENS')
               read(lucmd,'(A60)') ftmp%frozen%name
               itmp%im_frozen   = .true.

            CASE ('.UPDATE')
               itmp%im_update_vemb = .true.

            CASE ('.SKIPK')
               call fde_set_skip_nadd_ke(.true.)
            CASE ('.SKIPX')
               call fde_set_skip_nadd_xc(.true.)

! read kinetic energy functional for non-additive term
            CASE ('.NAXCF')
               read(lucmd,'(A80)') nadd_fun
               call fde_set_nadd_xcf(nadd_fun)

! read kinetic energy functional for non-additive term
            CASE ('.NAKEF')
               read(lucmd,*) nadd_fun
               call fde_set_nadd_kef(nadd_fun)

! flags for response read kinetic energy functional for non-additive term
            CASE ('.RSP')
               itmp%im_frozen = .true.
               dirac_cfg_fde_response = .true.

            CASE ('.RSPLDA')
               itmp%im_frozen = .true.
               call fde_set_nadd_all_alda
               dirac_cfg_fde_response = .true.

            CASE ('.LEVEL')
               read(lucmd,*) etmp%ex_level 

            CASE ('.GRIDOU')
               read(lucmd,*) ftmp%export%name
               etmp%do_grid_out = .true.
! doprp is part of dcbgen.h
               doprp = .true.

            CASE ('.EXONLY')
               read(lucmd,*) ftmp%export%name
               etmp%do_grid_out = .true.
! doprp is part of dcbgen.h
               skip_all_imports = .true.

#endif /* MOD_UNRELEASED */ 

            CASE DEFAULT
              WRITE (lunit,'(/,3A,/)') ' Keyword "',WORD, &
                 '" not recognized in FDEINP.'
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',lunit)
              CALL QUIT('Illegal keyword in FDEINP.')
            END SELECT

          CASE('*')
            EXIT

          CASE DEFAULT
            WRITE (lunit,'(/,3A,/)') ' Prompt "',WORD, &
               '" not recognized in FDEINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',lunit)
            CALL QUIT('Illegal prompt in FDEINP.')

          END SELECT
        END DO
      END IF
      
#ifdef MOD_UNRELEASED
      if (.not.skip_all_imports) then

      if (itmp%im_vemb) then
         if (itmp%im_update_vemb) then
            write (lunit,*) '.UPDATE and .EMBPOT cannot be used together!'
            call quit('conflicting keywords selected')
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
      endif

      if (etmp%do_grid_out) then
         etmp%ex_n        = .true.
         etmp%ex_coulomb  = .true.
         etmp%ex_gn       = .true.
      endif

      endif ! .not.skip_all_imports

#endif /* MOD_UNRELEASED */

! saving information on 
      call fde_set_files_info(ftmp)
      call fde_set_export_info(etmp)
      call fde_set_import_info(itmp)


      if (itmp%im_vemb) then
         call fde_import_static      
      endif

#ifdef MOD_UNRELEASED
      if (itmp%im_frozen) then
         call fde_import_frozen
         call fde_initialize_nadd_functionals         
      endif
#endif
      
!     options are echoed
         if (Echo_Input) then

           write (lunit,*) ' '
           CALL PRSYMB(lunit,'=',75,0)
           write (lunit,*) ' '

           write (lunit,'(1X,A)') &
              ' Frozen Density Embedding (FDE) Calculation setup'  

           call fde_show_setup 

           CALL PRSYMB(lunit,'=',75,0)
           write (lunit,*) ' '
         endif 

      END SUBROUTINE


      subroutine fde_show_setup
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
         
         CALL PRSYMB(unit,'-',48,2)
         write(unit,*)
         write(unit,'(A,I5)')  ' * Print level                  : ',pl

#ifdef MOD_UNRELEASED
         if (itmp%im_frozen) call fde_print_nadd_functionals(unit)
#endif

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

#ifndef MOD_UNRELEASED
         else
            call quit('inconsistent FDE input!')
#else
         else
            write(unit,'(A,A60)')  ' * FDE Potential generated from frozen, active densities'
            write(unit,'(A)') ' ' 
            call fde_test_frozen_density
         endif

         if (itmp%im_frozen) then
            string = ftmp%frozen%name
            write(unit,'(A,A60)') &
                ' * Density (and gradient) from frozen subsystems read from file: ',string
         endif

         if (etmp%do_grid_out) then
            string = ftmp%export%name
            write(unit,'(2A)') &
            ' * FDE Gridfile with updated density written to:',string

            select case(etmp%ex_level)
               case('DHF','MP2')
                  write(unit,*)                               &
                  'Outputted density will be taken from: ',trim(etmp%ex_level)
                  
               case default
                  write(unit,*) 'Input for outputted density not a &
     &  x  valid calculation, HF will be used instead. Given:',trim(etmp%ex_level)
                  etmp%ex_level = 'DHF'
                  
            end select
#endif /* MOD_UNRELEASED */
         endif
                  
      end subroutine fde_show_setup

end module fde_inputDirac
