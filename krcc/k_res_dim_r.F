      SUBROUTINE K_RES_DIM_REL(ISPOBEX,NSPOBEX,MAXOP,NELMNT_MAX)
*
* A set of spinorbital excitations, ISPOBEX_TP are given 
*
* Find largest dimension of operators that can be obtained 
* by removing 0-MAXOP operators from these 
*
* Jeppe Olsen, Nov. 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
*. Input
      INTEGER ISPOBEX(NGAS,4,NSPOBEX)
*. Local scratch
      INTEGER KSPOBEX(NGAS,4)
*
      NELMNT_MAX = 0
      DO JSPOBEX = 1, NSPOBEX
*. Number of elements in this TCC block
        NELMNT = 1
        DO ICAAB = 1,4
          DO IGAS = 1, NGAS
            NELMNT = 
     &      NELMNT*IBION(NOBPT(IGAS),ISPOBEX(IGAS,ICAAB,JSPOBEX))
          END DO
        END DO
        NELMNT_MAX = MAX(NELMNT,NELMNT_MAX)
*. This was the number of element in the unmodified type. 
*. Larger dimensions can be obtained only by removing operators from 
*. spaces that are more than half filled. 
*. Use this to locate the spaces where the largest gains are obtained 
*. from removing electrons, done one electron at a time 
        KELMNT_MAX = NELMNT
        CALL ICOPVE(ISPOBEX(1,1,JSPOBEX),KSPOBEX(1,1),4*NGAS)
        DO KOP = 1, MAXOP
*. Place where most will be gained by removing one operator
          KELMNTP = KELMNT_MAX
          JGAS = 0
          JCAAB = 0
          DO ICAAB = 1, 4
          DO IGAS = 1, NGAS
            IF(2*KSPOBEX(IGAS,ICAAB).GT.NOBPT(IGAS)) THEN
*. more than halffiled with operatos
              KELMNT = (KELMNTP*KSPOBEX(IGAS,ICAAB))
     &                /(NOBPT(IGAS)-KSPOBEX(IGAS,ICAAB)+1)
C             KELMNT = (KELMNTP*(NOBPT(IGAS)-KSPOBEX(IGAS,ICAAB)+1))/
C    &                  KSPOBEX(IGAS,ICAAB)
              IF(KELMNT.GT.KELMNT_MAX) THEN
                JGAS = IGAS
                JCAAB = ICAAB
                KELMNT_MAX  = KELMNT
              END IF
            END IF
          END DO
          END DO
          IF(JGAS.GT.0) THEN
            KSPOBEX(JGAS,JCAAB) = KSPOBEX(JGAS,JCAAB)-1
          ELSE
*. No more spaces more than halffilled with operators so
            GOTO 1001
          END IF
        END DO
*       ^ End of loop over KOP
 1001   CONTINUE
        NELMNT_MAX = MAX(NELMNT_MAX,KELMNT_MAX)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max number of operators to be removed ', MAXOP
        WRITE(6,*) ' Largest T block = ', NELMNT_MAX
      END IF
*
      RETURN
      END
