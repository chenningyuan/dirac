      FUNCTION ISYMS2(STRING,NEL)
*
* Symmetry of string STRING, D inf h, C inf v, O3 version
*
*
      IMPLICIT REAL*8(A-H,O-Z)
C     PARAMETER(MXPORB = 500 )
c     PARAMETER(MXPR4T = 10)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPIRR = 20)
c     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONAB/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*
      INTEGER   STRING(NEL)
*. ML and parity of string
      MLSTR = 0
      IPARI = 1
      DO 10 IEL = 1, NEL
        IF(ISMFTO(STRING(IEL)).LE.NMLOB) THEN
          MLSTR = MLSTR + ISMFTO(STRING(IEL))-1+MNMLOB
        ELSE
          MLSTR = MLSTR + ISMFTO(STRING(IEL))-1+MNMLOB-NMLOB
          IPARI = - IPARI
        END IF
   10 CONTINUE
*
      IF(IPARI.EQ.-1) IPARI = 2
      ISYM  = (IPARI-1) * NMLST+ MLSTR - MNMLST + 1
      ISYMS2 = ISYM
*
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' STRING '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,'(A,3I3)') ' MLSTR, IPARI ISYMS2 ', MLSTR,IPARI,ISYMS2
      END IF
*
      RETURN
      END
