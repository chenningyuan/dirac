      FUNCTION GTIJKL(I,J,K,L)
*
* Obtain  integral (I J ! K L )
* where I,J,K and l refers to active orbitals in
* Type ordering
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"

*.GLobal pointers
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
      IF(INTIMP .EQ. 2 ) THEN
*. LUCAS ordering
        I12S = 0
        I34S = 0
        I1234S = 1
CINA    GTIJKL = GIJKLL(IREOTS(1+NINOB),WORK(KPINT2),WORK(KLSM2),
CINA &           WORK(KINT2),
CINA &           ISMFTO,IBSO,NACOB,NSMOB,NOCOBS,I,J,K,L)
        GTIJKL = GIJKLL(IREOTS(1),WORK(KPINT2),WORK(KLSM2),
     &           WORK(KINT2),
     &           ISMFTO,IBSO,NACOB,NSMOB,NOCOBS,I,J,K,L)
       ELSE IF (INTIMP.EQ.1.OR.INTIMP.EQ.5) THEN
*. MOLCAS OR SIRIUS IMPORT ( I hope integrals are in core !! )
          GTIJKL = GMIJKL(I,J,K,L,WORK(KINT2),WORK(KPINT2))
       END IF
*
      RETURN
      END
