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

! $Id$
!
! description : basic fde module, handling i/o of imported and exported
!               data (potentials, densities and density gradients etc) on
!               a grid, and also for setting options (do embedding response,
!               regenerate the embedding potential, set kinetic energy 
!               functional and so on).
!
!               the module also keeps some internal variables (like the grid
!               and the imported quantities) that can be made available to
!               different parts of the code.
!
      MODULE FDE_MOD
      USE GRID_FILE
      use dft_input
      use dft_cfg

      PRIVATE
      PUBLIC FDE_ImportStaticPotential_FileOpen, FDE_ImportStaticPotential_FileClose
      public FDE_setupFDEgrid, FDE_testFrozendensity, FDE_cleanupFDEgrid
      public FDEINI,FDEINP,FDE_GRIDOUT,have_frozen_density
      PUBLIC DO_GRID_OUT
      public fde_grid, fde_frozen_prp, fde_grid_np, fde_nprop
      public skip_kin_contribution, skip_xc_contribution

      LOGICAL,SAVE :: DO_GRID_OUT=.FALSE.
      logical,save :: have_frozen_density = .false.
      logical,save :: import_vemb = .false.
      logical,save :: update_vemb = .false.
      logical,save :: skip_kin_contribution = .false.
      logical,save :: skip_xc_contribution  = .false.
      INTEGER,SAVE :: LUFDEPOT=45,LUTESTING=46,LUGRIDOUT=47, LUFDEFRZDNS=48
      INTEGER,SAVE :: print_fde=0, fde_grid_np=0, fde_nprop=1
      CHARACTER(LEN=60),SAVE :: FDEPOT='EMBPOT',FDEDNS='FRZDNS',GRIDOUT_NAME=''
      CHARACTER(LEN=60),SAVE :: out_level='HF'
      real(kind=8), pointer :: fde_grid(:,:), fde_frozen_prp(:,:) 

      
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_ImportStaticPotential_FileOpen: opens file containing the embedding potential
!                                    (static) exported on a grid
!
         subroutine FDE_ImportStaticPotential_FileOpen(embpot_unit)
            integer, intent(out) :: embpot_unit
            logical              :: embpot_found

            inquire(FILE=FDEPOT,EXIST=embpot_found)

            if (embpot_found) then
                open(lufdepot,               &
                     FORM   = 'FORMATTED',   &
                     STATUS = 'UNKNOWN',     &
                     ACCESS = 'SEQUENTIAL',  &
                     FILE   =  FDEPOT)
                embpot_unit = lufdepot
             else
                call quit("File with embedding potential not found!")
             endif

         end subroutine FDE_ImportStaticPotential_FileOpen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_ImportStaticPotential_FileClose: closes file containing the embedding potential
!                                      (static) exported on a grid
!
         subroutine FDE_ImportStaticPotential_FileClose()
            close(lufdepot) 
         end subroutine FDE_ImportStaticPotential_FileClose


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDEINI: performs the initialization in dirac; at the moment, it simply sets the
!         kinetic energy function to the default value
!
      SUBROUTINE FDEINI()
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
!
#include "memint.h"

      print_fde=0
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
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              CYCLE
            CASE ('.EMBPOT')
              import_vemb = .true.
!             read filename for the embedding potential
              READ(LUCMD,'(A60)') FDEPOT 

            CASE ('.PRINT')
!             set printlevel 
              READ(LUCMD,*) print_fde 

            CASE DEFAULT
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD, &
                 '" not recognized in FDEINP.'
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              CALL QUIT('Illegal keyword in FDEINP.')
            END SELECT
          CASE('*')
            EXIT
          CASE DEFAULT
            WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD, &
               '" not recognized in FDEINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in FDEINP.')
          END SELECT
        END DO
      END IF

!     options are echoed
         IF (Echo_Input) THEN
           write (lupri,*) ' '
           CALL PRSYMB(LUPRI,'=',101,1)

           WRITE(LUPRI,'(1X,2A,/1X,2A)') &
           'Frozen Density Embedding to include enviroment effects', &
           ' using the fixed potential method described in:',        &
           'A.S.P. Gomes, C. R. Jacob and L. Visscher, ',            &
           'Phys. Chem. Chem. Phys. 10 (2008) 5353–5362.'

           CALL PRSYMB(LUPRI,'=',101,1)

           WRITE(LUPRI,*)
           WRITE(LUPRI,'(1X,A)') &
              ' Frozen Density Embedding Calculation set-up:'  
           CALL PRSYMB(LUPRI,'-',44,2)
           WRITE(LUPRI,*)
           WRITE(LUPRI,'(A,I5)')  ' * Print level               : ',print_fde

           if (import_vemb) then
              write(LUPRI,'(A,A60)')  ' * FDE Potential read from file: ',FDEPOT
           endif

         END IF

      END SUBROUTINE


      END MODULE FDE_MOD
