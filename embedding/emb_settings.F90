MODULE emb_settings

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
!
! Has to be called AFTER DHFINP and HAMINP since it relies
! on input from these parts (DFT and frozen orbitals) 
!
! What it does for now is:
! Turn embedding (using the frozen orbitals) on if there are
! frozen orbitals and if the "embedding" functional is
! selected in the DFT part. 
!
!---------------------------------------------------------------------
 
    SUBROUTINE init_embedding ()
#include "priunit.h"

#include "dcbham.h"
#include "dcbdhf.h"

      LOGICAL, EXTERNAL :: DFTISEMB         ! external function from dft/functionals.c

      IF (DFTHAM .AND. DFTISEMB) THEN       ! DFT is used and an embedding functional is selected
          Do_Embedding = .TRUE. 
      ELSE
          Do_Embedding = .FALSE.
          Emb_Mode = EMB_MODE_OFF
      END IF

      IF (Do_Embedding) THEN                                   ! Decide which mode to use
          IF (NFRO_DHF(1) + NFRO_DHF(2) > 0) THEN 
              Emb_Mode = EMB_MODE_FROZEN_ORBITALS
          ELSE
              Do_Embedding = .FALSE.                    ! No suitable mode found, turn embedding off
          END IF
      END IF

    END SUBROUTINE init_embedding


!---------------------------------------------------------------------
! report_embedding
!---------------------------------------------------------------------
!
! Print some information on the embedding scheme used 
!
!---------------------------------------------------------------------

    SUBROUTINE report_embedding ()
#include "priunit.h"

      IF (Do_Embedding) THEN
          WRITE (LUPRI, *) "Embedding is switched ON"
          WRITE (LUPRI, *) "Embedding mode : ", Emb_Mode
      ELSE
          WRITE (LUPRI, *) "Embedding is switched OFF"
      END IF

    END SUBROUTINE report_embedding


END MODULE emb_settings
  
