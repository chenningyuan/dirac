      SUBROUTINE SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,
     &                  ISPGRPCON,IPRNT)
*
* FInd connection matrix for string types
*
* ISPGRPCON(ISPGP,JSPGRP) = 0 => spgrps are identical
*                         = 1 => spgrps are connected by single excitation
*      .                  = 2 => spgrps are connected by double excitation
*       .              . ge.3 => spgrps are connected by triple or
*        .                       higher excitation
*
* Jeppe Olsen, September 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IELFSPGRP(MXPNGAS,*)
*. output
      DIMENSION ISPGRPCON(NSPGRP,NSPGRP)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
      DO ISPGRP = 1, NSPGRP
        ISPGRPA = IOFSPGRP-1+ISPGRP
        DO JSPGRP = 1, ISPGRP
          JSPGRPA = IOFSPGRP-1+JSPGRP
          IDIF = 0
          DO IGAS = 1, NGAS
            IDIF = IDIF
     &    + ABS(IELFSPGRP(IGAS,ISPGRPA)-IELFSPGRP(IGAS,JSPGRPA))
          END DO
          NEXC = IDIF/2
          ISPGRPCON(ISPGRP,JSPGRP) = NEXC
          ISPGRPCON(JSPGRP,ISPGRP) = NEXC
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) '==================== '
        WRITE(6,*) 'output from SPGRPCON '
        WRITE(6,*) '==================== '
        WRITE(6,*)
        NEXC1 = 0
        NEXC2 = 0
        DO ISPGRP=1, NSPGRP
          DO JSPGRP = 1, NSPGRP
            IF(ISPGRPCON(ISPGRP,JSPGRP).EQ.1) THEN
              NEXC1 = NEXC1 + 1
            ELSE IF(ISPGRPCON(ISPGRP,JSPGRP).EQ.2) THEN
              NEXC2 = NEXC2 + 1
            END IF
          END DO
        END DO
*
        WRITE(6,*)
     &  ' single excitation interactions',NEXC1,
     &   '( ',float(NEXC1)*100.0/NSPGRP**2,' % ) '
        WRITE(6,*)
     &  ' double excitation interactions',NEXC2,
     &   '( ',float(NEXC2)*100.0/NSPGRP**2,' % ) '
*
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Supergroup connection matrix '
         CALL IWRTMA(ISPGRPCON,NSPGRP,NSPGRP,NSPGRP,NSPGRP)
      END IF
*
      RETURN
      END
