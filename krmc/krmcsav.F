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
C...   group of the licence holder may not take any part of Dalton,
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
C FILE    : krmcsav.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmcsav */
      SUBROUTINE KRMCSAV(IOPT,CMO,CREF,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Save info on LUKRMC
C
C     Input:
C        IOPT    - option (current only 0 supported)
C     
C     Written by J. Thyssen - Nov 26 1998
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcbgen.h"
C
      DIMENSION CMO(*), CREF(*)
      DIMENSION WORK(*)
      LOGICAL   FNDLAB, FND, RTRACI
C     
      CHARACTER CNOLBL*6
C
      CALL QENTER('KRMCSAV')
      KFRSAV = KFREE
C
C     Write file KRMCSCF:
C
C     (1) OLDORB  : old MO orbitals
C     (2) RESTART : info needed for restart
C     (3) LREDUCED: info needed for backsteps
C     (4) CREF    : new CI vector
C     (5) NEWORB  : new MO orbitals
C 
C     Clean LUKRMC
C     ------------
C     
      REWIND(LUKRMC)
      FND = FNDLAB('SODLABEL',LUKRMC)
      BACKSPACE LUKRMC
      CALL NEWLAB('SODLABEL',LUKRMC,LUPRI)
      REWIND(LUKRMC)
C
C     (1) OLDORB & OLDCREF: write old MO orbitals
C     -------------------------------------------
C
      CALL WRTKRMC(LUKRMC,'OLDORB  ',CMO,NCMOTQ)
      CALL WRTKRMC(LUKRMC,'OLDCREF ',CREF,NZCONFQ)
C
C     (2) RESTART: write restart info
C     -------------------------------
C
      CALL LABKRMC(LUKRMC,'RESTART ')
      WRITE(LUKRMC)
     &     EMCOLD, DEPRED, BETA, GAMMA, RTRUST, ITMAC,
     &     EMY, EACTIV, EMCOLD, EMCSCF
C
C     (3) LREDUCED: info nededed for backsteps
C     ----------------------------------------
C
      REWIND(LUKRM4)
      READ(LUKRM4) NREDM
C
      CALL MEMGET('REAL',KREDE,NREDM*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS,NREDM*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDGP,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVECR,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEVC,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBPVC,NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBCVC,NREDM,WORK,KFREE,LFREE)
C
C
      REWIND(LUKRM4)
      CALL REALRED(LUKRM4,NZRED,NREDM,WORK(KREDE),WORK(KREDGP),
     &     WORK(KEVECR),WORK(KIBTYP),
     &     WORK(KIBEVC),WORK(KIBPVC),WORK(KIBCVC),NPRED)
C
      CALL LABKRMC(LUKRMC,'LREDUCED')
      CALL WRTLRED(LUKRMC,NZRED,NREDM,WORK(KREDE),WORK(KREDGP),
     &     WORK(KEVECR),WORK(KIBTYP),
     &     WORK(KIBEVC),WORK(KIBPVC),WORK(KIBCVC),NPRED)
C
      CALL LABKRMC(LUKRMC,'EOFLABEL')
C
C     (4) - (5) CREF and NEWORB
C     -------------------------
C
C     Rotate to get new CREF and CMO
C     
      EVALR = D0
      CALL GTNORB(CMO,CREF,WORK(KIBTYP),
     &     WORK(KIBCVC),WORK(KIBEVC),WORK(KIBPVC),
     &     EVALR,WORK(KEVECR),WORK,KFREE,LFREE)
C
C     Rotate to Fock type and/or natural orbitals.
C
      IF (IKRMCCNO .GE. -1) THEN
C        
C
         IF (NATONL) THEN
            CNOLBL = 'ONLYNO'
         ELSE IF (FOCKON) THEN
            CNOLBL = 'ONLYFD'
         ELSE
            CNOLBL = 'FD+NO'
         END IF
C        
         IF (IKRMCCNO .EQ. -1) THEN
            RTRACI = .FALSE.
         ELSE
            RTRACI = .TRUE.
         END IF
C        
         CALL KRMCCNO(CNOLBL,RTRACI,CREF,CMO,IPROPT,WORK,KFREE,LFREE)
C
      END IF
C
      CALL WRTKRMC(LUKRMC,'CREF    ',CREF,NZCONFQ)
      CALL WRTKRMC(LUKRMC,'NEWORB  ',CMO,NCMOTQ)
C
      REWIND(LUKRMC)
C
      CALL MEMREL('KRMCSAV.1',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('KRMCSAV')
      RETURN
      END
