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

