      SUBROUTINE GETINL(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM)
*
* Obtain integrals XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
*.ORBINP
C     PARAMETER(MXPORB=500)
C     PARAMETER(MXPIRR = 20)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPR4T = 10)
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*.Output
      DIMENSION XINT(*)
*
      CALL QENTER('GETIN')
      NI = NTSOB(ITP,ISM)
      NJ = NTSOB(JTP,JSM)
      NK = NTSOB(KTP,KSM)
      NL = NTSOB(LTP,LSM)
*
      IF(IKSM.NE.0) THEN
        NIK = NI*(NI+1)/2
      ELSE
        NIK = NI*NK
      END IF
      IF(JLSM.NE.0) THEN
        NJL = NJ*(NJ+1)/2
      ELSE
        NJL = NJ*NL
      END IF
*
      I = 0
      K = 1
      DO 100 IK = 1, NIK
C       CALL NXTIJ(L,J,NL,NJ,JLSM,NONEW)
        CALL NXTIJ(I,K,NI,NK,IKSM,NONIK)
        IABS_T = I + IBTSOB(ITP,ISM)-1
        KABS_T = K + IBTSOB(KTP,KSM)-1
        J = 0
        L = 1
        DO 50 JL = 1, NJL
          CALL NXTIJ(J,L,NJ,NL,JLSM,NONJL)
          LABS_T = L + IBTSOB(LTP,LSM)-1
          JABS_T = J + IBTSOB(JTP,JSM)-1
C         XINT((IK-1)*NJL+JL) = GTIJKL(IABS_T,JABS_T,KABS_T,LABS_T)
C         IF(IXCHNG.EQ.1)
C    &    XINT((IK-1)*NJL+JL) = XINT((IK-1)*NJL+JL)
C    &    - GTIJKL(IABS_T,LABS_T,KABS_T,JABS_T)
          XINT((JL-1)*NIK+IK) = GTIJKL(IABS_T,JABS_T,KABS_T,LABS_T)
          IF(IXCHNG.EQ.1)
     &    XINT((JL-1)*NIK+IK) = XINT((JL-1)*NIK+IK)
     &    - GTIJKL(IABS_T,LABS_T,KABS_T,JABS_T)
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
      RETURN
      END
