C 
!
!      Copyright (c) 2008 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
C 
C
C
C VERSION : $Revision: 7732 $
C DATE    : $Date: 2008-09-16 22:30:19 +0200 (Tue, 16 Sep 2008) $
C FILE    : dirgeo.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamgeo */
      SUBROUTINE PAMGEO(WORK,LWORK)
C*****************************************************************************
C
C     Module for geometry optimization
C
C     Written by Joern Thyssen 1997/06/17
C     Last revision : 1997/06/17 jth
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbham.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbgrd.h"
C
      LOGICAL NOGRAD
C
      CALL QENTER('PAMGEO')
C
C     Determine whether analytical or numerical gradient is to 
C     be used. Analytical gradient is at the moment only available
C     for the HF module with the Levy-Leblond or DC Hamiltonian.
C
      NOGRAD= BSS.OR.ZORA.OR.(SPINFR .AND. .NOT.LEVYLE)
C
      IF(.NOT.DONGRD) THEN
C       We want the analytical gradient, check whether it is implemented
        DONGRD = NOGRAD.OR.DODFT.OR.DOMP2.OR.DORES.OR.DOCCM.OR.DOCIM.OR.
     &           DOLUCT.OR.DOKRMC.OR.DOLUCIAR
CTROND        DONGRD = .NOT. (DODHF .AND. HAMGRD)
        IF(DONGRD) THEN
          WRITE(LUPRI,'(/A/A)')
     &     ' * PAMGEO INFO: will use numerical molecular gradient',
     &     '         (analytical gradient not implemented for this wave'
     *     //' function type).'
        ENDIF
      ENDIF
      WRKDLM = 0.0D00
      IWOFF  = 0
      CALL OPTMIN(WORK,LWORK,IWOFF,WRKDLM)
C
      CALL QEXIT('PAMGEO')
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck exedir */
      SUBROUTINE EXEDIR(WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/07/09 by J. Thyssen
C
C     EXEcute DIRac; Based on PAMDRV with a few modifications.
C
C     Last revision: 1997/07/08 jth
C
C*****************************************************************************
#include "implicit.h"
      DIMENSION WORK(LWORK)
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C 
#include "priunit.h"
#include "dcbgen.h"
#include "dcbdhf.h"
C
      CALL QENTER('EXEDIR')
C
C     Delete DFCYCL so we don't restart DIIS
C     also delete MCCRES to avoid restarting in relccsd.
C     
      CALL DELFILE(LUCYCL,'DFCYCL')
      CALL DELFILE(16,'MCCRES')
C
C     To avoid readin twice in the 0. iteration skip
C     the readin/setup.
C
C     Input and initialization
C
      IF ( FIRST ) THEN
         FIRST = .FALSE.
      ELSE
C
C***********************************************************************
C*****   I  N  P  U  T    S  E  C  T  I  O  N   ************************
C***********************************************************************
C
         CALL PAMINP(WORK,LWORK)
         IF ( INPTES ) RETURN
C
C*****************************************************************************
C*****  S E T U P  -  M O D U L E  *******************************************
C*****************************************************************************
C
         CALL PAMSET(WORK,LWORK)
      END IF
C
C
C******************************************************************************
C*****  G E T    W A V E    F U N C T I O N   *********************************
C******************************************************************************
C
C     Always do wave function
      DOPSI = .TRUE.
C
C
      IF(DOPIN) CALL PCMOIN(WORK,LWORK)
C
C     Since there is no fun in doing a geometry optimization
C     without a wave function, PAMPSI is always called.
C
      CALL PAMPSI(WORK,LWORK)
      IF(DOACUT) CALL ACMOUT(WORK,LWORK)
      IF(DOPUT) CALL PCMOUT(WORK,LWORK)
C
C
C*****************************************************************************
C*****  A N A L Y S I S    M O D U L E  **************************************
C*****************************************************************************
C
      IF(DOANA) CALL PAMANA(WORK,LWORK)
C
C
C*****************************************************************************
C*****  R E S P O N S E    M O D U L E  **************************************
C*****************************************************************************
C
C
      IF (DOPRP) CALL PAMPRP(WORK,LWORK)
C
C*****************************************************************************
C***** H E I N Z *************************************************************
C*****************************************************************************
C
      IF(DOHEINZ) CALL HEINZ(WORK,LWORK)
C
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('EXEDIR')
      RETURN
      END
C  /* Deck rgetccsd */
      FUNCTION RGETCCSD()
#include "implicit.h"
#include "priunit.h"
C
C     Get energy from RELCCSD
C
      CHARACTER CENERGY*40
#include "../relccsd/results.inc" 
#include "../relccsd/inpt.inc" 
      IF (DOCCSDT) THEN
         CENERGY = 'CCSD(T) energy from RELCCSD'
      ELSE IF (DOCCSD) THEN
         CENERGY = 'CCSD energy from RELCCSD'
      ELSE IF (DOMP2) THEN
         CENERGY = 'MP2 energy from RELCCSD'
      ELSE
         call quit('no energy from RELCCSD')
      END IF
      WRITE(LUPRI,'(2A,1P,E20.10)')
     &   CENERGY,' = ',ETOT
      RGETCCSD = ETOT
      RETURN
      END
C  /* Deck rgetccsd */
      FUNCTION RGETCI()
C***********************************************************************
C
C     Get energy from DIRRCI
C     Taking the lowest root....
C     
C     Written by T. Saue Dec 9 2004
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "../dirrci/param.inc" 
#include "../dirrci/iterat.inc" 
      RGETCI = EIGEN2(1)
      RETURN
      END
