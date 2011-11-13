C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
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
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C
C
C VERSION : $Revision: 1.15.2.2 $
C DATE    : $Date: 2004/09/01 08:08:46 $
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
      LOGICAL HAMGRAD_IMPLEMENTED
C
      CALL QENTER('PAMGEO')
C
C     Determine whether analytical or numerical gradient is to 
C     be used. Analytical gradient is at the moment only available
C     for the HF module with the Levy-Leblond or DC Hamiltonian.
C
      HAMGRAD_IMPLEMENTED = .NOT. (ZORA.OR.(SPINFR .AND. .NOT.LEVYLE))
C
      IF(.NOT.DONGRD) THEN
C       We want the analytical gradient, check whether it is implemented
        DONGRD = .NOT. (DODHF .AND. HAMGRAD_IMPLEMENTED)
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
