      SUBROUTINE TRP_CITT_BLK(CIN,COUT,ISM,NASTR,NBSTR,NSMST,ICOPY)
*
* Transpose all symmetry blocks of a TT block of a CI expansion
* C(IA,IB) => C(IB,IA)
*
* All blocks are assumed to be present in complete form in CIN
*
*. Order of blocks is dictated by symmetry of first index
*
* Jeppe Olsen, HNIE, July 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "symm.inc"
*. Input
      DIMENSION CIN(*) 
      INTEGER NASTR(NSMST),NBSTR(NSMST)
*. Output
       DIMENSION COUT(*)
*. Local scratch
      INTEGER IOUT_OFF(MXPNSMST)
*. Set up offsets for output block
      IOFF = 1
      DO IBSM = 1, NSMST
        IASM = IDBGMULT(ISM,IBSM)
        IOUT_OFF(IBSM) = IOFF
        IOFF = IOFF + NBSTR(IBSM)*NASTR(IASM)
      END DO
*
      IOFF_IN = 1
      DO IASM = 1, NSMST
        IBSM = IDBGMULT(IASM,ISM)
        NA = NASTR(IASM)
        NB = NBSTR(IBSM)
        CALL TRPMT3(CIN(IOFF_IN),NA,NB,COUT(IOUT_OFF(IBSM)))
        IOFF_IN = IOFF_IN + NA*NB
      END DO
*
      IF(ICOPY.EQ.1) THEN
        LEN = IOFF_IN - 1 
        CALL COPVEC(COUT,CIN,LEN)
      END IF
*
      RETURN
      END
