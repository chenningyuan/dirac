C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirgeo.F
C SOURCE  : $Source$
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
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <dcbham.h>
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbgrd.h>
C
      CALL QENTER('PAMGEO')
C
C     Determine whether analytical or numerical gradient is to 
C     be used...
C
      IF(.NOT.DONGRD) THEN
        DONGRD = DOCCM.OR.DOMP2.OR.(DODHF.AND.SPINFR)
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
#include <implicit.h>
      DIMENSION WORK(LWORK)
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C 
#include <priunit.h>
#include <dcbgen.h>
#include <dcbdhf.h>
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
#include <implicit.h>
#include <priunit.h>
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
