      SUBROUTINE GETINT_REL(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                      IXCHNG,IKSM,JLSM,ICOUL)
*
* Interface to DIRAC relativistic integrals
* (in spinfree formalism) implemented
*
*---------------------------------------------
* Timo Fleig, October 2000
*---------------------------------------------
*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*./CRUN/ : INTIMP used
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,NP1,NP2,NQ,INCORE,MXCIV,ICISTR,
C    &            NOCSF,IDIAG
#include "crun.inc"
*./ORBINP/  : NOBPTS used
*.Memory
#include "wrkspc.inc"
*
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cgas.inc"

C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"
*
      CALL QENTER('GETIN')
      NTEST = 000
*
      IF(INTIMP.EQ.1.OR.INTIMP.EQ.5.or.INTIMP.eq.6) THEN
*. Interface to MOTRA or SIRIUS integrals
*  or to DIRAC (spinfree) integrals
        IF(INCORE.EQ.0) THEN
        CALL GETINM(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM)
        ELSE IF (INCORE.NE.0) THEN
*. Old
COLD      CALL GETINC(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
COLD &                  IXCHNG,IKSM,JLSM,WORK(KINT2),
COLD &                  WORK(KPINT2),NSMOB)
*
COLD       WRITE(6,*) ' Integrals with old fetcher : '
COLD    NI = NOBPTS(ITP,ISM)
COLD    NK = NOBPTS(KTP,KSM)
COLD    IF(IKSM.EQ.0) THEN
COLD      NIK = NI * NK
COLD    ELSE
COLD      NIK = NI*(NI+1)/2
COLD    END IF
COLD    NJ = NOBPTS(JTP,JSM)
COLD    NL = NOBPTS(LTP,LSM)
COLD    IF(JLSM.EQ.0) THEN
COLD      NJL = NJ * NL
COLD    ELSE
COLD      NJL = NJ*(NJ+1)/2
COLD    END IF
COLD    WRITE(6,*) ' 2 electron integral block for TS blocks '
COLD    WRITE(6,*) ' Ixchng :', IXCHNG
COLD    WRITE(6,*) ' After GETINC '
COLD    WRITE(6,'(1H ,4(A,I2,A,I2,A))')
COLD &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
COLD &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
COLD    CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
C
*. New
        CALL GETINCN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &               IXCHNG,IKSM,JLSM,WORK(KINT2),
     &               WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
        END IF
      else if (INTIMP.EQ.2) THEN
C. Interface to LUCAS integrals
        CALL GETINL(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM)
      else if (INTIMP.eq.8) then
*
* Molecule Sweden interface
*
        call GETINCN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &               IXCHNG,IKSM,JLSM,WORK(KINT2),
     &               WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
      END IF
*
      IF(NTEST.NE.0) THEN
        NI = NGSSH(ISM,ITP)
        NK = NGSSH(KSM,KTP)
        IF(IKSM.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
        NJ = NGSSH(JSM,JTP)
        NL = NGSSH(LSM,LTP)
        IF(JLSM.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,*) ' Ixchng :', IXCHNG
        WRITE(6,*) ' After GETINC '
        WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
C     STOP ' Jeppe forced me to stop in GETINT '
      RETURN
      END
