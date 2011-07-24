      SUBROUTINE LEN_FOR_TF_REL(NSPOBEX,ISPOBEX,LSPOBEX,NGAS,LEN_TF,
     &                          WORK,KFREE,LFREE)
*
* Find largest Length of a T block and the associated F-vector
*
* Jeppe Olsen, November 2001
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Input
      INTEGER ISPOBEX(4*NGAS,NSPOBEX), LSPOBEX(NSPOBEX)
*. Local scratch
      INTEGER ITF_OCC(4*MXPNGAS)
*
      DIMENSION WORK(*)
*
      NTEST = 000
      LEN_TF = 0
      DO IFTP = 1, NSPOBEX+1
        LEN_F = 0
        DO ITTP = 1, NSPOBEX
*. Occupation of ITTP*IFTP
          CALL OP_T_OCC_REL(ISPOBEX(1,ITTP),ISPOBEX(1,IFTP),
     &                  ITF_OCC,IMZERO)
          IADD = 0
          IF(IMZERO.EQ.0) THEN
*. Address of ITF, or number of excitations required to 
*. bring it into space
            CALL INUM_FOR_OCC2_REL(ITF_OCC,INUM,NDIFF,
     &                             WORK,KFREE,LFREE)
*. Space will be allocated only if TF is not in space, 
*. and can be reached by atmost 4 operators
            IF(0.LT.NDIFF.AND.NDIFF.LE.4) IADD = 1
          END IF
          IF(IADD.EQ.1) LEN_F = LEN_F + LSPOBEX(ITTP)
        END DO
        LEN_TF = MAX(LEN_TF,LEN_F+LSPOBEX(IFTP))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Largest length of an T-block and F-vector',
     &               LEN_TF
      END IF
*
      RETURN
      END
