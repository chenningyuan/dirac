      FUNCTION ISYMS1_REL(STRING,NEL)
*
* Symmmetry of string, D2H version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
C     PARAMETER(MXPORB = 500 )
C     PARAMETER(MXPR4T = 10)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPIRR = 20 )
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
*. Specific input
      INTEGER STRING(*)
*
      ISYM = 1
      DO 100 IEL = 1, NEL
        ISYM = SYMPRO(ISYM,ISMFTO(STRING(IEL)))
  100 CONTINUE
*
      ISYMS1_REL = ISYM
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
*
      RETURN
      END
