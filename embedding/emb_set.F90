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
!      DIRAC program, see: http://wiki.chem.vu.nl/dirac
!dirac_copyright_end

MODULE emb_settings !MI: module name can not be identical with the file name

  IMPLICIT NONE

  LOGICAL, SAVE :: Do_Embedding         ! Is embedding turned on ?

  ! possible vaules for Emb_Mode
  INTEGER, PARAMETER :: EMB_MODE_OFF = 0                   ! embedding is off
  INTEGER, PARAMETER :: EMB_MODE_FROZEN_ORBITALS = 1       ! use the frozen orbitals (from *DHFCAL .FROZEN)
  ! more to come ...

  INTEGER, SAVE :: Emb_Mode             ! The mode to use for the embedding module 

CONTAINS

!---------------------------------------------------------------------
! init_embedding
!---------------------------------------------------------------------
! 
! Initialize the settings for the embedding module  
! (turn everything of)
!
! This is called at the beginning of HAMINP
!
!---------------------------------------------------------------------
 
    SUBROUTINE init_embedding ()
      IMPLICIT NONE

      Do_Embedding = .FALSE.
      Emb_Mode = EMB_MODE_OFF
      
    END SUBROUTINE init_embedding


    SUBROUTINE input_embedding (word)
      IMPLICIT NONE

#include <priunit.h>

      CHARACTER(len=7), INTENT(INOUT)  :: word

      LOGICAL           :: done
      INTEGER           :: inperr
      CHARACTER(len=70) :: line

      done = .false.

      ! turn embedding on and set defaults
      Do_Embedding = .TRUE.
      Emb_Mode = EMB_MODE_FROZEN_ORBITALS

      DO WHILE (.NOT.done)

         read(LUCMD, '(A7)') word

         if (word(1:1) == '*') then
            done = .true.
         else if ((word(1:1) /= '!') .and. (word(1:1) /= '#')) then
            
            SELECT CASE (word)
               CASE ('.PRORBS')  ! use the projected-out orbitals as frozen density
                                 ! (Since this is the only optin right now, it is also 
                                 !  the default)
                  emb_Mode = EMB_MODE_FROZEN_ORBITALS
               CASE ('.XCFUN ')  ! non-additive XC-functional
                  read (lucmd, '(A70)') line
                  inperr = 0
                  call dft_emb_set_xcfunc (line, inperr)  
                  if (inperr > 0) then
                     call quit("Unknown functional in *EMBEDDING .XCFUNC")
                  end if

               CASE ('.KINFUN')  ! non-additive kinetic energy functional

                  read (lucmd, '(A70)') line
                  inperr = 0
                  call dft_emb_set_kinfunc (line, inperr) 
                  if (inperr > 0) then
                     call quit("Unknown functional in *EMBEDDING .KINFUNC")
                  end if

               CASE DEFAULT
                  write (lupri, '(/,3A,/)') 'Prompt "', word, '" not recognized in emb_input'
                  call quit("Illegal prompt in emb_input")
            END SELECT

         end if

      END DO
 
      ! Check if we have got everything we need, set defaults if needed

    END SUBROUTINE input_embedding


!---------------------------------------------------------------------
! report_embedding
!---------------------------------------------------------------------
!
! Print some information on the embedding scheme used 
!
! I also (mis)use it to make some consistency checks
!
! Has to be called AFTER DHFINP and HAMINP since it relies
! on input from these parts (DFT and frozen orbitals) 
!
!---------------------------------------------------------------------

    SUBROUTINE report_embedding ()
      use dirac_cfg
      IMPLICIT NONE

#include "priunit.h"
#include "dcbham.h"
#include "dcbdhf.h"

      write (lupri, *) " "
      write (lupri, *) "****** Embedding *******"

      IF (Do_Embedding) THEN
          WRITE (LUPRI, *) "Embedding is switched ON"
          WRITE (LUPRI, *) "Embedding mode : ", Emb_Mode

          ! some consistency checks

          IF (.not. dirac_cfg_dft_calculation) THEN
             CALL quit ("embedding only works with DFT at the moment")
          end IF

          IF (Emb_Mode == EMB_MODE_FROZEN_ORBITALS) THEN
             IF (NFRO_DHF(1) + NFRO_DHF(2) == 0) THEN
                call quit ("EMB_MODE_FROZEN_ORBITALS selected, but no frozen orbitals")
             end IF
          END IF

!         call DftEmbReport ()   !MI: this routine is nowhere defined....

      ELSE
          WRITE (LUPRI, *) "Embedding is switched OFF"
      END IF

    END SUBROUTINE report_embedding


END MODULE emb_settings
  
