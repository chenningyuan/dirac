      SUBROUTINE DIAG_EXC_CC(NCA,NCB,NAA,NAB,NGAS,IDIAG)
*
* Check of alpha and betaexcitation parts are identical
*
*. Output
* IDIAG = 0 : Not diagonal in alpha, beta
* IDIAG = 1 : Diag in alpha, beta
*
*
#include "implicit.inc"
*. Input
      INTEGER NCA(NGAS),NCB(NGAS),NAA(NGAS),NAB(NGAS)
*
      IDIAG = 1 
      DO IGAS = 1, NGAS
        IF(NCA(IGAS).NE.NCB(IGAS)) IDIAG = 0
        IF(NAA(IGAS).NE.NAB(IGAS)) IDIAG = 0
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCA, NCB, NAA, NAB = '
        CALL IWRTMA(NCA,1,NGAS,1,NGAS)
        CALL IWRTMA(NCB,1,NGAS,1,NGAS)
        CALL IWRTMA(NAA,1,NGAS,1,NGAS)
        CALL IWRTMA(NAB,1,NGAS,1,NGAS)
        WRITE(6,*) ' IDIAG = ', IDIAG
      END IF
*
      RETURN 
      END
