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

   use fde_cfg
   use fde_xcfun_interface
   
   private

   public FDEINI,FDEINP

   CONTAINS

      SUBROUTINE FDEINI()
         integer :: unit

         call fde_get_log_unit(unit)
         call fde_set_have_frozen(.true.)
         call fde_set_vemb_update(.true.)

         call fde_initializeFunctionals

         write (unit,*) ' '
         CALL PRSYMB(unit,'=',101,1)
         write (unit,*) ' '
         write(unit,'(1X,A)') &
            ' Frozen Density Embedding Calculation defaults:'
         call FDE_setupInfo(unit)    
         CALL PRSYMB(unit,'=',101,1)
         write (unit,*) ' '
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDEINP: dirac-specific input reader for the module
!
      SUBROUTINE FDEINP(WORD,Echo_Input)
#include "implicit.h"
#include "priunit.h" 
#include "dcbgen.h"
      PARAMETER (NTABLE = 8)
!
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, line*80
!
!
      logical :: Echo_Input,newdef
      integer :: printlevel
!

      call fde_set_vemb_import(.false.)
      call fde_set_frzdns_import(.false.)

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
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',fde_log_unit)
              CYCLE

            CASE ('.EMBPOT')
              call fde_parse_vemb(lucmd)
            CASE ('.FRDENS')
              call fde_parse_frozen(lucmd)

            CASE ('.PRINT')
              READ(LUCMD,*) printlevel 
              call fde_set_print(printlevel)

            CASE ('.SKIPK')
               call set_skip_kin(.true.)
            CASE ('.SKIPX')
               call set_skip_xc(.true.)

! read kinetic energy functional for non-additive term
            CASE ('.NAXCF')
               call parse_xcf(lucmd)
! read kinetic energy functional for non-additive term
            CASE ('.NAKEF')
               call parse_kef(lucmd)
! read kinetic energy functional for non-additive term
            CASE ('.RSPLDA')
               call set_fde_alda(.true.)

            CASE ('.LEVEL')
               call parse_export_level(lucmd)
            CASE ('.GRIDOU')
               call set_export_file(lucmd)
! doprp is part of dcbgen.h
               doprp = .true.

            CASE DEFAULT
              WRITE (fde_log_unit,'(/,3A,/)') ' Keyword "',WORD, &
                 '" not recognized in FDEINP.'
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',fde_log_unit)
              CALL QUIT('Illegal keyword in FDEINP.')
            END SELECT
          CASE('*')
            EXIT
          CASE DEFAULT
            WRITE (fde_log_unit,'(/,3A,/)') ' Prompt "',WORD, &
               '" not recognized in FDEINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',fde_log_unit)
            CALL QUIT('Illegal prompt in FDEINP.')
          END SELECT
        END DO
      END IF

!     options are echoed
         if (Echo_Input) then
           write (fde_log_unit,*) ' '
           CALL PRSYMB(fde_log_unit,'=',101,1)
           write (fde_log_unit,*) ' '
           write(fde_log_unit,'(1X,A)') &
              ' Frozen Density Embedding Calculation, user-defined setup:'  
           call FDE_setupInfo(fde_log_unit) 
           CALL PRSYMB(fde_log_unit,'=',101,1)
           write (fde_log_unit,*) ' '
         endif 

      END SUBROUTINE


      subroutine FDE_setupInfo(unit)
         integer :: unit
         character(len=60) :: string

         call fde_get_log_unit(unit)

         CALL PRSYMB(unit,'-',44,2)
         write(unit,*)
         write(unit,'(A,I5)')  ' * Print level                : ',print_fde

         call fde_printFunctionals

         if (fde_import_vemb) then
            call fde_get_vemb_filename(string)
            write(unit,'(A,A60)')  ' * FDE Potential to be read from file: ',string
            write(unit,*)
            write(unit,'(1X,2A,/1X,2A)') &
            'Frozen Density Embedding to include enviroment effects', &
            ' using the fixed potential method described in:',        &
            'A.S.P. Gomes, C. R. Jacob and L. Visscher, ',            &
            'Phys. Chem. Chem. Phys. 10 (2008) 5353â5362.'

         else
            write(unit,'(A,A60)')  ' * FDE Potential generated from frozen, active densities'
            call fde_set_vemb_update(.true.)
            call fde_set_have_frozen(.true.)
         endif

         if (fde_have_frozen_density) then
            call fde_get_frozen_filename(string)
            write(unit,'(A,A60)') &
                ' * Density (and gradient) from frozen subsystems read from file: ',string
         endif

         IF (fde_export_n) then
            call fde_get_export_filename(string)
            write(unit,'(2A)') &
            ' * FDE Gridfile with updated density written to:',string

         endif
      end subroutine FDE_setupInfo


end module fde_inputDirac
