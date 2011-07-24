      SUBROUTINE STRING_IN_CCEXCIT_REL(ISTRING,ICAAB,INSPC,
     &                             NSPOBEX_TP,ISPOBEX_TP)
*
* Is String ISTRING ( given as occ in each gas space)
* include as CAAB string ICAAB in any of the NSPOBEX_TP TCC excitations 
* blocks given in ISPOBEX_TP 
*
* Jeppe Olsen, May 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISTRING(NGAS)
*
      INSPC = 0
      IB = 1 + (ICAAB-1)*NGAS
      DO KTP = 1, NSPOBEX_TP
        IDEL = 0
        DO IGAS = 1, NGAS
          IDEL = IDEL + ABS(ISPOBEX_TP(IB-1+IGAS,KTP)-ISTRING(IGAS))
        END DO
        IF(IDEL.EQ.0) THEN
          INSPC = 1
          GOTO 1001 
        END IF
      END DO
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ICAAB of string = ', ICAAB
        WRITE(6,*) ' Occupation : '
        CALL IWRTMA(ISTRING,1,NGAS,1,NGAS)
        IF(INSPC.EQ.1) THEN
           WRITE(6,*) ' String is in space '
        ELSE
           WRITE(6,*) ' String is out of space '
        END IF
      END IF
*
      RETURN 
      END
