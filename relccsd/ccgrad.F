C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXPVAL (DOO,DVV,DVO,DOV,EOO,EVV,EVO,EOV,EPHASE,CEXPA)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates expectation value of an operator.
C
C---------------Routines called----------------------------------------
C
C     BLAS routines
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 DOO(*),DVV(*),DVO(*),DOV(*)
      REAL*8 EOO(*),EVV(*),EVO(*),EOV(*)
      COMPLEX*16 CEXPA,EPHASE
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 CDOTC,ZDOTC
      REAL*8 DDOT,ZDOT
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CEXPA = ZDOTU (NFVO,EOV,1,DVO,1)
      ELSE
         EXPA = DDOT (NFVO,EOV,1,DVO,1)
      ENDIF
c     print*,' Relaxation contribution 1',expa
C
      IF (CARITH) THEN
         CEXPA = CEXPA + ZDOTU (NFVO,EVO,1,DOV,1)
      ELSE
         EXPA = EXPA + DDOT (NFVO,EVO,1,DOV,1)
      ENDIF
c     print*,' Relaxation contribution ',expa
C
      IF (CARITH) THEN
         CEXPA = CEXPA + ZDOTU (NFOO,EOO,1,DOO,1)
      ELSE
         EXPA = EXPA + DDOT (NFOO,EOO,1,DOO,1)
      ENDIF
c     print*,' OO contribution ',DDOT (NFOO,EOO,1,DOO,1)
C
      IF (CARITH) THEN
         CEXPA = CEXPA + ZDOTU (NFVV,EVV,1,DVV,1)
      ELSE
         EXPA = EXPA + DDOT (NFVV,EVV,1,DVV,1)
      ENDIF
c     print*,' VV contribution ',DDOT (NFVV,EVV,1,DVV,1)
c     print*,' Total MP2 contribution ',expa
C
C     This is the correlation contribution, add the reference value
C
      II = 1
      DO IREP = 1, NREP
         IF (CARITH) THEN
            CEXPA = CEXPA + ZDOTU (NO(IREP),EOO(II),NO(IREP)+1,A1,0)
         ELSE
            EXPA = EXPA + DDOT (NO(IREP),EOO(II),NO(IREP)+1,A1,0)
         ENDIF
         II = II + NO(IREP)*NO(IREP)*RCW
      ENDDO
C
      IF (CARITH) THEN
         CEXPA = EPHASE * CEXPA
      ELSE
         CEXPA = EPHASE * EXPA
      ENDIF
C
      RETURN
      END
