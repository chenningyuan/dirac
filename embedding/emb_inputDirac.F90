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
         
      end subroutine fdeini


!   ----------------------------------------------------------------------------
      subroutine fdeinp(word,echo_input)
!   ----------------------------------------------------------------------------
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
      
      itmp%im_vemb     = .false.
      itmp%im_frozen   = .true.
      itmp%im_n        = .true.
      itmp%im_gn       = .true.
      itmp%im_coulomb  = .true.
      itmp%im_update_vemb = .true.


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
               itmp%im_vemb     = .true.
               itmp%im_frozen   = .false.
               itmp%im_n        = .false.
               itmp%im_gn       = .false.
               itmp%im_coulomb  = .false.
               itmp%im_update_vemb = .false.
!              call fde_parse_vemb(lucmd)
            CASE ('.FRDENS')
               read(lucmd,'(A60)') ftmp%frozen%name
               itmp%im_vemb     = .false.
               itmp%im_frozen   = .true.
               itmp%im_n        = .true.
               itmp%im_gn       = .true.
               itmp%im_coulomb  = .true.
               itmp%im_update_vemb = .true.
!              call fde_parse_frozen(lucmd)

            CASE ('.PRINT')
              read(lucmd,*) pl 
              call fde_set_print_level(pl)

            CASE ('.SKIPK')
               call fde_set_skip_nadd_ke(.true.)
            CASE ('.SKIPX')
               call fde_set_skip_nadd_xc(.true.)

! read kinetic energy functional for non-additive term
            CASE ('.NAXCF')
               read(lucmd,'(A80)') nadd_fun
               write (*,*) 'debug, nadd_fun is ',nadd_fun
               
               call fde_set_nadd_xcf(nadd_fun)
! read kinetic energy functional for non-additive term
            CASE ('.NAKEF')
               read(lucmd,*) nadd_fun
               call fde_set_nadd_kef(nadd_fun)
! read kinetic energy functional for non-additive term
            CASE ('.RSPLDA')
               call fde_set_nadd_all_alda

            CASE ('.LEVEL')
               read(lucmd,*) etmp%ex_level 
!               call fde_set_export_level(lucmd)
            CASE ('.GRIDOU')
               read(lucmd,*) ftmp%export%name
               etmp%do_grid_out = .true.
               etmp%ex_n        = .true.
               etmp%ex_coulomb  = .true.
               etmp%ex_gn       = .true.
!               call fde_set_export_file(lucmd)
! doprp is part of dcbgen.h
               doprp = .true.

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

! saving information on 
      call fde_set_files_info(ftmp)
      call fde_set_export_info(etmp)
      call fde_set_import_info(itmp)
      
      if (itmp%im_frozen) call fde_import_frozen

      call fde_initialize_nadd_functionals         

!     options are echoed
         if (Echo_Input) then
           write (lunit,*) ' '
           CALL PRSYMB(lunit,'=',101,1)
           write (lunit,*) ' '
           write (lunit,'(1X,A)') &
              ' Frozen Density Embedding Calculation, user-defined setup:'  
           call FDE_show_setup 
           CALL PRSYMB(lunit,'=',101,1)
           write (lunit,*) ' '
         endif 


      END SUBROUTINE


      subroutine FDE_show_setup
         integer :: unit
         character(len=60) :: string
         type(fde_export)  :: etmp
         type(fde_import)  :: itmp
         type(fde_files)   :: ftmp
         logical :: have_frozen
         integer :: pl

         call fde_get_print_level(pl)

         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
         call fde_get_import_info(itmp)

         unit = ftmp%logfile%unit
         
         CALL PRSYMB(unit,'-',44,2)
         write(unit,*)
         write(unit,'(A,I5)')  ' * Print level                : ',pl

         call fde_print_nadd_functionals(unit)

         if (itmp%im_vemb) then
            string = ftmp%embpot%name
            
            write(unit,'(A,A60)')  ' * FDE Potential to be read from file: ',string
            write(unit,*)
            write(unit,'(1X,2A,/1X,2A)') &
            'Frozen Density Embedding to include enviroment effects', &
            ' using the fixed potential method described in:',        &
            'A.S.P. Gomes, C. R. Jacob and L. Visscher, ',            &
            'Phys. Chem. Chem. Phys. 10 (2008) 5353â5362.'

         else
            write(unit,'(A,A60)')  ' * FDE Potential generated from frozen, active densities'
            itmp%im_frozen = .true.
            itmp%im_update_vemb = .true.
            call fde_test_frozen_density
         endif

!         call fde_get_have_frozen(have_frozen)

         if (itmp%im_frozen) then
            string = ftmp%frozen%name
            
!           call fde_get_frozen_filename(string)
            write(unit,'(A,A60)') &
                ' * Density (and gradient) from frozen subsystems read from file: ',string
         endif

         if (etmp%do_grid_out) then
            string = ftmp%export%name
            
!            call fde_get_export_filename(string)
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
         endif
         
         if ((fde_xc_fun_is_lda().and.fde_ke_fun_is_gga()).or. &
             (fde_xc_fun_is_gga().and.fde_ke_fun_is_lda())) then
             call quit('inconsistent selection of KE+XC functionals')
         else 
             write (*,*) 'debug, fde_xc_fun_is_lda() ',fde_xc_fun_is_lda()
             write (*,*) '       fde_ke_fun_is_gga() ',fde_ke_fun_is_gga()
             write (*,*) '       fde_xc_fun_is_gga() ',fde_xc_fun_is_gga()
             write (*,*) '       fde_ke_fun_is_lda() ',fde_ke_fun_is_lda()
         endif
         
! saving information on 
         call fde_set_files_info(ftmp)
         call fde_set_export_info(etmp)
         call fde_set_import_info(itmp)
         
      end subroutine FDE_show_setup

end module fde_inputDirac
