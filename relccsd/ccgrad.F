C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccgrad.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDENSC(IOPT,T1,S1,T2,S2,
     &                   BUF1,BUF2,BUF3,DOO,DVV,DVO,DOV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Menu driven routine. Depending on value of IOPT
C
C     Calculate unperturbed density matrix at the
C       - SCF level (0)
C       - MP2 level (1)
C       - CCSD level (2)
C       - CCSD(T) level (3)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IOPT
      REAL*8 T1(*),S1(*),T2(*),S2(*)
      REAL*8 DOO(*),DVV(*),DVO(*),DOV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NFOO,A0,0,DOO,1)
      CALL XCOPY (NFVV,A0,0,DVV,1)
      IF (IOPT.LT.1) RETURN
C----------------------------------------------------------------------
C D(I,J) = D(I,J) - T(CD,IL) * L(CD,JL) = D(I,J) - T(CDL,I) * L(CDL,J)
C----------------------------------------------------------------------
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF1)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,-A1,BUF1,BUF2,A1,DOO,NREP)
C----------------------------------------------------------------------
C D(I,J) = D(I,J) - T(C,I) * L(C,J) 
C----------------------------------------------------------------------
      CALL CNTRCT ('C','N',NO,NO,NV,-A1,T1,S1,A1,DOO,NREP)
C----------------------------------------------------------------------
C D(A,B) = D(A,B) + L(AD,KL) * T(BD,KL) = H(A,C) + L(A,DKL) * T(B,DKL)
C----------------------------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,A1,BUF1,BUF2,A1,DVV,NREP)
C----------------------------------------------------------------------
C D(A,B) = D(A,B) + L(A,K) * T(B,K)
C----------------------------------------------------------------------
      CALL CNTRCT ('N','C',NV,NV,NO,A1,T1,S1,A1,DVV,NREP)
C
      IF (CARITH) THEN
         CALL CONJUGA (NFOO,DOO,1)
         CALL CONJUGA (NFVV,DVV,1)
      ENDIF
C
      IF (IOPT.LT.2) RETURN
C----------------------------------------------------------------------
C DVO(A,I) = L(A,I)
C----------------------------------------------------------------------
      CALL XCOPY (NDIMT1,S1,1,DVO,1)         
C----------------------------------------------------------------------
C DOV(A,I) = T(A,I)
C----------------------------------------------------------------------
      CALL XCOPY (NDIMT1,T1,1,DOV,1)         
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) + T(AC,IK) * L(C,K)
C----------------------------------------------------------------------
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)
CTROND      CALL XGEMV('N',NDIMT1,NDIMT1,NDIMT1,A1,BUF1,M,S1,1,A1,DOV,1)
      CALL XGEMV ('N',NDIMT1,NDIMT1,A1,BUF1,M,S1,1,A1,DOV,1)
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - T(A,K) * L(C,K) * T(C,I)
C----------------------------------------------------------------------
      CALL CNTRCT ('C','N',NO,NO,NV,-A1,S1,T1,A0,BUF1,NREP)
      CALL CNTRCT ('N','N',NV,NO,NO,A1,T1,BUF1,A1,DOV,NREP)
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - T(A,K) * L(CD,LK) * T(CD,LI)
C----------------------------------------------------------------------
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF1)
      CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,-A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL CNTRCT ('N','N',NV,NO,NO,A1,T1,BUF3,A1,DOV,NREP)
C----------------------------------------------------------------------
C DOV(A,I) = DOV(A,I) - L(CD,KL) * T(AD,KL) * T(D,I)
C----------------------------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL CNTRCT ('N','N',NV,NO,NV,A1,BUF3,T1,A1,DOV,NREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCDENSZ(IOPT,NSP,IREPSPI,EPS,T1,S1,T2,S2,
     &                   BUF1,BUF2,BUF3,NBUF2,NBUF3,CBUF,
     &                   XVO,XOV,DOO,DVV,DVO,DOV,
     &                   EVO,EOV,D2VO,D2OV,DG,CC1,BB1,BB2,CONV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Menu driven routine. Depending on value of IOPT
C
C     Solve Z-vector equations at the
C       - SCF level (0)
C       - MP2 level (1)
C       - CCSD level (2)
C       - CCSD(T) level (3)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER IOPT
      REAL*8 T1(*),S1(*),T2(*),S2(*)
      REAL*8 DOO(*),DVO(*),DOV(*),DVV(*)
      REAL*8 D2VO(*),D2OV(*),DG(*)
      REAL*8 XVO(*),XOV(*)
      REAL*8 EVO(*),EOV(*)
      REAL*8 EPS(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*),CBUF(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "inpt.inc"
      INCLUDE "symm.inc"
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*8 NAME
      COMPLEX*16 EPHASE(MAXOP)
      LOGICAL BDONE,USEDZ
C
C---------------Executable code--------------------------------------
C
      IF (IOPT.LT.1) RETURN
C
C     Calculate X-intermediate
C
      CALL XCOPY (NDIMX,A0,0,XVO,1)
C--------------------------------------------
C X(A,I) = X(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      CALL XCOPY (NDIMT2,S2,1,BUF3,1)
      CALL T1EQN1 (BUF1,BUF2,BUF3,NBUF3,XVO)
C--------------------------------------------
C X(A,I) = X(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      CALL XCOPY (NDIMT2,S2,1,BUF1,1)
      CALL T1EQN2 (BUF1,BUF2,BUF3,NBUF3,XVO)
C---------------------------------------------------------------------
C X(A,I) = X(A,I) + V(IK,AL) * D(K,L)
C---------------------------------------------------------------------
      CALL GETVOOO (BUF1)
      CALL SRT1ST4 (NREP,MULTB,LTR,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              BUF1,BUF2)
      M = NFVO
      K = NFOO
      CALL XGEMV ('N',M,K,A1,BUF2,M,DOO,1,A1,XVO,1)
C------------------------------------
C X(A,I) = X(A,I) + W(CI,DA) * D(C,D)
C------------------------------------
      M = NFVV
      N = NFVO
      IF (N.EQ.0) RETURN
      IF (NBUF3.LT.JVVVO(2)) THEN
         CALL QUIT(" OUT OF CORE STEP IN GV NOT IMPLEMENTED")
      ENDIF
      USEDZ = .FALSE.
      CALL SRT20 (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,BUF3,USEDZ)
      CALL XGEMV ('C',M,N,A1,BUF3,M,DVV,1,A1,XVO,1)
C
      CALL XSCAL(NDIMX,-A1,XVO,1)
      CALL XCOPY (NDIMX,XVO,1,XOV,1)
      IF (CARITH) CALL CONJUGA (NDIMX,XOV,1)
C
C     Solve Z-vector equations by calling the RPA routine
C
      NAME = 'MP2-GRAD'
      DO EOPER = 1, NEOPER
         EOFF = (EOPER-1)*NDIMX*RCW+1
         CALL RDPROP(IPRNT,NAMEE(EOPER),NSP,LFA,CBUF,IREPSPI,
     &               EPHASE(EOPER),EVO(EOFF),EOV(EOFF),DUM,DUM)
      ENDDO
C
C     Silence RPA routine when the print level is smaller than one
C
      IF (IPRNT.GE.1) THEN
         JPRNT = IPRNT
      ELSE
         JPRNT = -1
      ENDIF
C
      BDONE = .FALSE.
      CALL RPAA (NEOPER,NAMEE,NAME,EPHASE,A1,BDONE,EPS,A0,
     &           EVO,EOV,XVO,XOV,DVO,D2VO,DOV,D2OV,DG,BUF1,BUF2,
     &           CC1,BB1,BB2,CONV,JPRNT,DEBUG,TIMING)
C
      RETURN
      END
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
