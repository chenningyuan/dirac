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
C
C
      CALL QENTER('PAMGEO')
C
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
