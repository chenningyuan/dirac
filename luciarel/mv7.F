      SUBROUTINE MV7(C,HC,LUC,LUHC,INONREL)
*
* Outer routine for sigma vector generation
* GAS version !!!!
*
* Written in terms of RASG3/SBLOCK, May 1997
* Works only for ICISTR .ne. 1, i.e. pure memory version
* has been eliminated.
*
* Generalized to allow for relativistic sigma, T.F. + J.O., march 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      common/glbref/ISPC,ISM,ITERGLB
*
CTF  Variable INONREL added to list. It is passed from CMICDV
*
      CALL QENTER('MV7  ')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MV7   ')
*
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' MV7 does not work for ICISTR = 1'
        WRITE(6,*) ' SWITCH to ICISTR = 2,3 or program'
        Call Abend2( ' MV7 does not work for ICISTR = 1' )
      END IF
*. Arrays for partitioning of sigma
* (*_MS2 arrays are in CSTATE )
      CALL Z_BLKFO(IDC,NMS2VAL,
     &             ISPC,ISM,KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,KLSBLTP,
     &             NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2)
*
      CALL RASSG3(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),LUC,LUHC)
C          RASSG3(CB,SB,LBATS,LEBATS,I1BATS,IBATS,
C    &     LUC,LUHC)
*. Eliminate local memory
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'MV7   ')
*
      CALL QEXIT('MV7  ')
*
      RETURN
      END
