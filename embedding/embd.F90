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

MODULE Embedding

  USE Emb_Settings
  USE Emb_Frozen_Orbitals

  IMPLICIT NONE
  PRIVATE 

  PUBLIC :: Do_Embedding, init_embedding, input_embedding, report_embedding  ! From MODULE Emb_Settings 

  PUBLIC :: start_Embedding, end_Embedding
  PUBLIC :: Get_Frozen_Rho, Get_Frozen_Grad_Rho

CONTAINS

  SUBROUTINE start_Embedding ()
    IMPLICIT NONE

    IF (Do_Embedding) THEN
        SELECT CASE (Emb_Mode)
        CASE (EMB_MODE_FROZEN_ORBITALS)
            CALL init_Emb_Frozen_Orbitals ()
        CASE DEFAULT
            CALL QUIT ("Unknown Emb_Mode in start_Embedding") 
        END SELECT
    END IF

  END SUBROUTINE start_Embedding


  SUBROUTINE end_Embedding ()
    IMPLICIT NONE

    IF (Do_Embedding) THEN
        SELECT CASE (Emb_Mode)
        CASE (EMB_MODE_FROZEN_ORBITALS)
            CALL dealloc_Emb_Frozen_Orbitals ()
        CASE DEFAULT
            CALL QUIT ("Unknown Emb_Mode in end_Embedding") 
        END SELECT 
    END IF

  END SUBROUTINE end_Embedding


  
  SUBROUTINE Get_Frozen_Rho (Frozen_Rho, GAO)
    IMPLICIT NONE

#include "dcbbas.h"
    
    REAL*8,                      Intent(OUT) :: Frozen_Rho
    REAL*8, Dimension(NTBAS(0)), Intent(IN)  :: GAO

     SELECT CASE (Emb_Mode)
    CASE (EMB_MODE_FROZEN_ORBITALS)
        CALL Get_Rho_Frozen_Orbitals (Frozen_Rho, GAO)
    CASE (EMB_MODE_OFF)
        CALL QUIT ("Get_Frozen_Rho called in Emb_Mode = EMB_OFF - This should never happen")
    CASE DEFAULT
        CALL QUIT ("Get_Frozen-Rho called with unknown Emb_Mode - This should never happen")
    END SELECT

  END SUBROUTINE Get_Frozen_Rho



  SUBROUTINE Get_Frozen_Grad_Rho (Frozen_Grad_Rho, GAO, GAO1)
    IMPLICIT NONE

#include "dcbbas.h"

    REAL*8, Dimension(3),           Intent(OUT)   :: Frozen_Grad_Rho
    REAL*8, Dimension(NTBAS(0)),    Intent(IN)    :: GAO
    REAL*8, Dimension(NTBAS(0), 3), Intent(IN)    :: GAO1

    SELECT CASE (Emb_Mode)
    CASE (EMB_MODE_FROZEN_ORBITALS)
        CALL Get_Grad_Rho_Frozen_Orbitals (Frozen_Grad_Rho, GAO, GAO1)
    CASE (EMB_MODE_OFF)
        CALL QUIT ("Get_Frozen_Rho called in Emb_Mode = EMB_OFF - This should never happen")
    CASE DEFAULT
        CALL QUIT ("Get_Frozen-Rho called with unknown Emb_Mode - This should never happen")
    END SELECT

  END SUBROUTINE Get_Frozen_Grad_Rho


END MODULE Embedding
