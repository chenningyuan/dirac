C  /* Deck drive_fit */
      SUBROUTINE DRIVE_FIT(DENS,ADENS,LDM,INDAO,INDAOS,INDNAO,IPRINT,
     &                     CFIT)
C*****************************************************************************
C    
C     DRIVE_FIT: Driver for density fitting in Dirac.
C        
C     Written by Luuk Visscher, november 2003
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      DIMENSION DENS(LDM,*),ADENS(LDM,*),CFIT(*)
      DIMENSION INDAO(*),INDAOS(*)
      DIMENSION INDNAO(2,*)
#include "nuclei.h"
#include "symmet.h"
C
#include "ibtfun.h"
C
C     ***************************************
C     ***** Sort the Density matrix      ****
C     ***************************************
C
C     The Density matrix is already in AO-basis but we need it
C     with all functions on one center consecutively.
C
      CALL REORMT (DENS,ADENS,LDM,LDM,INDAO,INDAO)
C
C     **************************************************
C     ***** Loop over symmetry unique center pairs  ****
C     **************************************************
C
      print*,"in driver of density fitting "
C     Loop over symmetry unique centers
      DO IB = 1, NUCIND
         MULB = ISTBNU(IB)
C        Loop over symmetry unique centers
         DO IA = 1, NUCIND
            MULA = ISTBNU(IA)
            MAB  = IBTOR (MULA,MULB)
C           Run over symmetry independent charge distributions
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MAB) .NE. 0) CYCLE
               ICENTA = NUCNUM(IA,1)
               ICENTB = NUCNUM(IB,ISYMOP+1)
C              The LL part of the density matrix
               INDA = INDNAO(1,ICENTA)+1
               INDB = INDNAO(1,ICENTB)+1
C              We need a pointer array that gives the start of the fit functions for a center
               IFTA = (ICENTA-1)*2 + 1
               IFTB = (ICENTB-1)*2 + 1
               print*,"fitting atom pair :",icenta,icentb
               CALL FitPairDensity (ICENTA,ICENTB,ADENS(INDA,INDB),LDM,
     &                              CFIT(IFTA), CFIT(IFTB), IERR)
C              The SS part of the density matrix
               INDA = INDNAO(2,ICENTA)+1
               INDB = INDNAO(2,ICENTB)+1
               CALL FitPairDensity (ICENTA,ICENTB,ADENS(INDA,INDB),LDM,
     &                              CFIT(IFTA), CFIT(IFTB), IERR)
            END DO  
         END DO 
      END DO 
C
      RETURN
      END
