C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2EQN(EPS,FVO,T1,S1,T2,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves MP2 amplitude or lambda equations.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.
C
      REAL*8 EPS(*),FVO(*),T1(*),S1(*),T2(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "eqns.inc"
      INCLUDE "inpt.inc"
      INCLUDE "results.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CPUDELTA,CPU0,CPUE,CPUVAR,CPUTOT
      REAL*8 DIAG
      LOGICAL TEQ
C
C---------------Executable code--------------------------------------
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSE
         STOP 'Wrong keyword in call to MP2EQNS'
      ENDIF
C
C     Initialize the timing.
C
      CPU0 = CPUDELTA()
      CALL CPUUSED(CPU0)
      CPUE = A0
C
C  INITIALIZE T2 AND T1 AMPLITUDES
C
      CALL XCOPY (NDIMT1,FVO,1,T1,1)
      CALL GETVVOO (T2)
      CALL DENOM (EPS,T1,T2,T1,T2)
      CALL ZCORE (T1,T2)
C
      IF (TEQ) THEN
C
C  EVALUATE CORRELATION ENERGY 
C
         CALL MP2ENG (T1,T2,FVO,S2,EMP2)
         CALL T1DIAG (DIAG,T1)
         CPUE = CPUE + CPUDELTA()
C
         WRITE(IW,6071)
         WRITE(IW,6072) ESCF
         WRITE(IW,6073) EMP2
         ETOT = ESCF + EMP2
         WRITE(IW,6075) ESCF + EMP2
         WRITE(IW,6077) DIAG
C
      ELSE
C
C  INITIALIZE LAMDA AMPLITUDES
C
         CALL XCOPY (NDIMT1,T1,1,S1,1)
         CALL XCOPY (NDIMT2,T2,1,S2,1)
C
      ENDIF
C
      CALL CPUUSED(CPUTOT)
      CPUTOT = CPUTOT-CPU0
      CPUVAR = CPUTOT-CPUE
      IF (TIMING) WRITE(IW,7000) CPU0,CPUE,CPUVAR,CPUTOT
C
      RETURN
 6071 FORMAT(//,'  MP2 results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 7000 FORMAT(//' Timing of routine MP2EQN :'
     &/' Before MP2EQN :',T30,F12.3,' seconds'
     &/' Energy, T1-Diag :',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in MP2EQN :',T30,F12.3,' seconds')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCEQN(EPS,FOO,FVO,FVV,T1,S1,T2,S2,AZ,H,
     &                 HO,HV,HOV,GOO,GVV,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,CONV,
     &                 NBUF2,NBUF3)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves CCSD amplitude or lambda equations.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.
C
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),T1(*),S1(*),T2(*),S2(*)
      REAL*8 AZ(*),H(*),HO(*),HV(*),HOV(*),BUF1(*),BUF2(*),BUF3(*)
      REAL*8 GOO(*),GVV(*)
C
C     for diis ...
C
      REAL*8 CC1(*)
      REAL*8 BB1(*),BB2(*)
C
      REAL*8 RMS,CONV
      REAL*8 DIAG
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "diis.inc"
      INCLUDE "eqns.inc"
      INCLUDE "inpt.inc"
      INCLUDE "results.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CPUDELTA,CPU0,CPUT1,CPUT2,CPUINT,CPU4V,CPUE,CPUVAR,CPUTOT
      LOGICAL TEQ
      real*8 ddot,sdot
      complex*16 zdotc,cdotc
C
C---------------Executable code--------------------------------------
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSE
         STOP 'Wrong keyword in call to CCEQNS'
      ENDIF
C
C     Initialize the timing.
C
      CPU0 = CPUDELTA()
      CALL CPUUSED(CPU0)
      CPUE = A0
      CPUINT = A0
      CPUT1 = A0
      CPUT2 = A0
      CPU4V = A0
C
C     Calculate record size for DIIS files.
C
      ITC = 0
      IT = 0
      NST1 = INTOWP(RCW*NDIMT1)/SEC2I(1)
      IF(NST1*SEC2I(1).NE.INTOWP(RCW*NDIMT1)) NST1 = NST1 + 1
      NST2 = INTOWP(RCW*NDIMT2)/SEC2I(1)
      IF(NST2*SEC2I(1).NE.INTOWP(RCW*NDIMT2)) NST2 = NST2 + 1
C
C     Open DIIS scratch files
C
      NINTSB = 20
      CALL RFILE(ITAMPL,NINTSB)
      CALL RFILE(ITAMPD,NINTSB)
C
      NIT = 0
      RMS = 1.D0
C
      IF (IPRNT.GE.1) THEN
         WRITE(IW,6003)
         WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      ENDIF
C
C     Calculate fixed intermediaries for Lambda equations
C
      IF (.NOT.TEQ) THEN
         CALL GETAMPT (S1,S2)
         CALL HOINTM (FOO,S1,S2,BUF1,BUF2,BUF3,HO)
         CALL HVINTM (FVV,S1,S2,BUF1,BUF2,BUF3,HV)
         CALL HOVINTM (FVO,S1,BUF1,BUF2,HOV)
         CALL GOINTM (FVO,S1,BUF1,BUF2,HO)
         CALL GVINTM (FVO,S1,BUF3,NBUF3,BUF2,NBUF2,HV)
         CALL AINTM (S1,S2,BUF1,BUF2,AZ)
         CALL HINTM (S1,S2,BUF1,BUF2,BUF3,NBUF3,H)
         CPUINT = CPUINT + CPUDELTA()
      ENDIF
C
C  START ITERATIVE PROCEDURE
C
 1000 CONTINUE
      NIT = NIT + 1
      IF (NIT.GT.MAXIT) THEN
         WRITE(IW,6004) RMS
         GO TO 2000
      ENDIF
C
      IF (.NOT.TEQ) THEN
C        
C        Make Gpq intermediates for lambda equations
C
         CALL GETAMPT (S1,S2)
         CALL MKGPP (S2,T2,BUF1,BUF2,GOO,GVV)
      ENDIF
C
C  EVALUATE T1 EQUATIONS
C  CHECK FOR CASES WITH NO NON_ZERO T1'S TO AVOID CRASHES
C
      IF (TEQ) CALL HOINTM (FOO,T1,T2,BUF1,BUF2,BUF3,HO)
      IF (TEQ) CALL HVINTM (FVV,T1,T2,BUF1,BUF2,BUF3,HV)
      IF (NDIMT1.GT.0) THEN
         IF (TEQ) THEN
          CALL HOVINTM (FVO,T1,BUF1,BUF2,HOV)
          CALL T1EQNS (FVO,T1,T2,HO,HV,HOV,H,GOO,GVV,
     &                 BUF1,BUF2,BUF3,NBUF3,S1)
         ELSE
          CALL T1EQNS (HOV,T1,T2,HO,HV,HOV,H,GOO,GVV,
     &                 BUF1,BUF2,BUF3,NBUF3,S1)
         ENDIF
      ENDIF
      CPUT1 = CPUT1 + CPUDELTA()
 111  format (a,2f30.15)
      if (carith.and.debug) then
C64B  write (iw,111)"ho",cdotc(nfoo,ho,1,ho,1)
C64B  write (iw,111)"hv",cdotc(nfvv,hv,1,hv,1)
C64B  write (iw,111)"hov",cdotc(ndimt1,hov,1,hov,1)
C64B  write (iw,111)"s1",cdotc(ndimt1,s1,1,s1,1)
      write (iw,111)"ho",zdotc(nfoo,ho,1,ho,1)
      write (iw,111)"hv",zdotc(nfvv,hv,1,hv,1)
      write (iw,111)"hov",zdotc(ndimt1,hov,1,hov,1)
      write (iw,111)"s1",zdotc(ndimt1,s1,1,s1,1)
      elseif (debug) then
C64B  write (iw,111)"ho",sdot(nfoo,ho,1,ho,1)
C64B  write (iw,111)"hv",sdot(nfvv,hv,1,hv,1)
C64B  write (iw,111)"hov",sdot(ndimt1,hov,1,hov,1)
C64B  write (iw,111)"s1",sdot(ndimt1,s1,1,s1,1)
      write (iw,111)"ho",ddot(nfoo,ho,1,ho,1)
      write (iw,111)"hv",ddot(nfvv,hv,1,hv,1)
      write (iw,111)"hov",ddot(ndimt1,hov,1,hov,1)
      write (iw,111)"s1",ddot(ndimt1,s1,1,s1,1)
      endif
C
C  EVALUATE T2 EQUATIONS
C
      IF (TEQ) THEN
         CALL GOINTM (FVO,T1,BUF1,BUF2,HO)
         CALL GVINTM (FVO,T1,BUF3,NBUF3,BUF2,NBUF2,HV)
         CALL AINTM (T1,T2,BUF1,BUF2,AZ)
         CALL HINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
         CPUINT = CPUINT + CPUDELTA()
      ENDIF
C
      if (carith.and.debug) then
      write (iw,111)"go",zdotc(nfoo,ho,1,ho,1)
      write (iw,111)"gv",zdotc(nfvv,hv,1,hv,1)
      write (iw,111)"az",a2*a2*zdotc(ioooott(nrep+1),az,1,az,1)
      write (iw,111)"h ",zdotc(ivovo(nrep+1),h,1,h,1)
      elseif (debug) then
      write (iw,111)"go",ddot(nfoo,ho,1,ho,1)
      write (iw,111)"gv",ddot(nfvv,hv,1,hv,1)
      write (iw,111)"az",4.d0*ddot(ioooott(nrep+1),az,1,az,1)
      write (iw,111)"h ",ddot(ivovo(nrep+1),h,1,h,1)
      endif
C
C  INTERMEDIATES FORMED --- T2 EQUATIONS
C
      CALL T2EQNS (T1,T2,HO,HV,HOV,H,AZ,GOO,GVV,BUF1,BUF2,BUF3,NBUF3,S2)
      CPUT2 = CPUT2 + CPUDELTA()
C
C
C  MUST NOW ADD IN B INTERMEDIATE CONTRIBUTION
C
      if (carith.and.debug) then
      write (iw,111)"s2",zdotc(ndimt2,s2,1,s2,1)
      elseif (debug) then
      write (iw,111)"s2",ddot(ndimt2,s2,1,s2,1)
      endif
C      
C     Use HOV array as buffer to store T1 amplitudes, when we
C     are solving the lambda equations
C
      CALL BINTM (T1,T2,HOV,BUF1,BUF2,NBUF2,BUF3,NBUF3,S2)
      if (carith.and.debug) then
      write (iw,111)"s2",zdotc(ndimt2,s2,1,s2,1)
      elseif (debug) then
      write (iw,111)"s2",ddot(ndimt2,s2,1,s2,1)
      endif
      CPU4V = CPU4V + CPUDELTA()
C
C  DIVIDE NEW T1s AND T2s BY EIGENVALUES
C
      CALL DENOM (EPS,S1,S2,S1,S2)
      CALL ZCORE (S1,S2)
C
C  EXTRAPOLATE THE AMPLITUDES
C
      CALL DIISX(ITAMPL,ITAMPD,NDIMT1,NDIMT2,T1,S1,T2,S2,IT,
     &           CC1,BB1,ITC,BB2,RMS,NST1,NST2)
C
C  EVALUATE CORRELATION ENERGY AND WRITE OUT RESULTS
C
      IF (TEQ) THEN
         CALL CCENG  (S1,S2,FVO,BUF1,BUF2,ECCSD)
      ENDIF
C
      CALL T1DIAG (DIAG,S1)
      IF (IPRNT.GE.1) WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      CALL FLSHFO(IW)
C
C  PUT S1,S2 INTO T1,T2 AND DO NEXT ITERATION
C
      CALL XCOPY(NDIMT1,S1,1,T1,1)
      CALL XCOPY(NDIMT2,S2,1,T2,1)
C
      CPUE = CPUE + CPUDELTA()
C
      IF(RMS.LT.CONV) GO TO 2000
      GO TO 1000
C
 2000 CONTINUE
      CALL T1DIAG(DIAG,T1)
      CPUE = CPUE + CPUDELTA()
      IF (TEQ) THEN
         ETOT = ESCF + ECCSD
         WRITE(IW,6071)
         WRITE(IW,6072) ESCF
         WRITE(IW,6074) ECCSD
         WRITE(IW,6076) ETOT
      ENDIF
C
      WRITE(IW,6077) DIAG
      WRITE(IW,6078) RMS
      WRITE(IW,6079) NIT
      CALL CPUUSED(CPUTOT)
      CPUTOT = CPUTOT-CPU0
      CPUVAR = CPUTOT-CPUT1-CPUT2-CPUINT-CPU4V-CPUE
      IF (TIMING)
     &WRITE(IW,7000) CPU0,CPUINT,CPUT1,CPUT2,CPU4V,CPUE,CPUVAR,CPUTOT
C
C     Close and delete DIIS scratch files
C
      CALL RCLOSE(ITAMPL,4)
      CALL RCLOSE(ITAMPD,4)
C
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6002 FORMAT(3X,'MP2   =',F20.15/3X,'EMP2  =',F20.15/)
 6003 FORMAT(//3X,'NIT',10X,'ENERGY',20X,'RMS',7X,'T1-DIAGN')
 6004 FORMAT(//' Requested convergence not reached in CCSD : RMS = ',
     &F20.15)
 6071 FORMAT(//,'  CCSD results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6074 FORMAT(' CCSD correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6076 FORMAT(' Total CCSD energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 6078 FORMAT(' Convergence :',T40,F25.15)
 6079 FORMAT(' Number or iterations used :',T60,I5)
 7000 FORMAT(//' Timing of routine CCEQN :'
     &/' Before CCEQN :',T30,F12.3,' seconds'
     &/' Intermediates :',T30,F12.3,' seconds'
     &/' T1 equations :',T30,F12.3,' seconds'
     &/' T2 equations :',T30,F12.3,' seconds'
     &/' 4 virtual contr. :',T30,F12.3,' seconds'
     &/' Energy, DIIS, T1-Diag, etc.:',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in CCEQN :',T30,F12.3,' seconds')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1EQNS (FVO,T1,T2,HO,HV,HOV,H,GOO,GVV,BUF1,BUF2,BUF3,
     &                   NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate T1 or Lambda1 equations
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 FVO(*)
      REAL*8 T1(*),T2(*)
      REAL*8 S1(*)
      REAL*8 HO(*),HV(*),HOV(*),H(*)
      REAL*8 GOO(*),GVV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "eqns.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ
      REAL*8 SDOT,DDOT
      COMPLEX*16 CDOTU,ZDOTU,CT1
C
C---------------Executable code--------------------------------------
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSE
         STOP 'Wrong keyword in call to L1EQNS'
      ENDIF
C-------------------------------
C S(A,I) = F(A,I)
C-------------------------------
      CALL XCOPY (NDIMT1,FVO,1,S1,1)
C--------------------------------
C Terms specific for T1 equations
C--------------------------------
      IF (TEQ) THEN
C------------------------------------------------------------
C S(A,I) = S(A,I) + T(A,K) * [HOV(K,C) - 2 * F(K,C)] * T(C,I)
C------------------------------------------------------------
         CALL XCOPY (NFVO,HOV,1,BUF1,1)
         IF (CARITH) CALL CONJUGA (NFVO,BUF1,1)
         CALL XAXPY (NFVO,-A2,FVO,1,BUF1,1)
         OFF1 = 1
         DO IRP = 1, NREP
            M = NV(IRP)
            N = NO(IRP)
            K1 = NO(IRP)
            K2 = NV(IRP)
            CALL XGEMM ('C','N',K1,N,K2,A1,BUF1(OFF1),K2,T1(OFF1),K2,
     &                   A0,BUF2,K1)
            CALL XGEMM ('N','N',M,N,K1,A1,T1(OFF1),M,BUF2,K1,
     &                   A1,S1(OFF1),M)
            OFF1 = OFF1 + M * N * RCW
         ENDDO
C-----------------------------------------------------------------
C S(A,I) = S(A,I) + HOV(CK) * T2(A,C,I,K)
C-----------------------------------------------------------------
         CALL SRT1 (NREP,MULTB,NV,NV,NO,NO,MVO,JVOVO,JJVO,JJVO,
     &              T2,BUF1,T2,BUF1)
         M = MVO(1)
         K = MVO(1)
         CALL XGEMV ('N',M,K,A1,BUF1,M,HOV,1,A1,S1,1)
      ENDIF
C-----------------------------------------------
C S(A,I) = S(A,I) + HV(A,C) * T(C,I)
C-----------------------------------------------
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NO(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,HV(OFF1),M,T1(OFF2),K,
     &                A1,S1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C-----------------------------------------------
C S(A,I) = S(A,I) - T(A,K) * HO(K,I)
C-----------------------------------------------
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NO(IRP)
         K = NO(IRP)
         CALL XGEMM ('N','N',M,N,K,-A1,T1(OFF1),M,HO(OFF2),K,
     &                A1,S1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF3)
      ELSE
         CALL XCOPY (NDIMT2,T2,1,BUF3,1)
      ENDIF
      CALL T1EQN1 (BUF1,BUF2,BUF3,NBUF3,S1)
C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF1)
      ELSE
         CALL XCOPY (NDIMT2,T2,1,BUF1,1)
      ENDIF
      CALL T1EQN2 (BUF1,BUF2,BUF3,NBUF3,S1)
C--------------------------------------
C S(A,I) = S(A,I) - W(A,K,C,I) * T(C,K)
C--------------------------------------
      IF (TEQ) THEN
         CALL GETVOVO (BUF1)
         CALL SRT7 (NREP,MULTB,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &              BUF1,BUF2,BUF1,BUF2)
      ELSE
         CALL XCOPY (JVOVO(2),H,1,BUF2,1)
      ENDIF
      M = MVO(1)
      K = MVO(1)
      CALL XGEMV ('N',M,K,-A1,BUF2,M,T1,1,A1,S1,1)
C------------------------------------
C Terms specific for Lamba1 equations
C------------------------------------
      IF (.NOT.TEQ) THEN
C--------------------------------------
C S(A,I) = S(A,I) - TV(C,D,L,K) * G(C,D)
C--------------------------------------
c and more see notes
C--------------------------------------
C S(A,I) = S(A,I) - W(C,A,D,I) * G(C,D)
C--------------------------------------
      M = NFVV
      N = NFVO
      IF (NBUF3.LT.JVVVO(2)) THEN
         STOP " OUT OF CORE STEP IN GV NOT IMPLEMENTED"
      ENDIF
      CALL SRT20 (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,BUF3)
      CALL XGEMV ('C',M,N,-A1,BUF3,M,GVV,1,A1,S1,1)
C--------------------------------------
C S(A,I) = S(A,I) + W(A,K,L,I) * G(K,L)
C--------------------------------------
      CALL GETVOOO (BUF1)
      CALL SRT2 (NREP,MULTB,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &           BUF1,BUF2)
      M = NFVO
      K = NFOO
      CALL XGEMV ('N',M,K,A1,BUF2,M,GOO,1,A1,S1,1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T2EQNS (T1,T2,GO,GV,HOV,H,AZ,GOO,GVV,
     &                   BUF1,BUF2,BUF3,NBUF3,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate T2 equations
C
C---------------Routines called----------------------------------------
C
C     BLAS routines :                 XGEMM, XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 T1(*),T2(*)
      REAL*8 S2(*)
      REAL*8 GO(*),GV(*),HOV(*)
      REAL*8 GOO(*),GVV(*)
      REAL*8 AZ(*),H(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "eqns.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ
      REAL*8 T1VAL(2)
      complex*16 zdotc
      real*8 ddot
C
C---------------Executable code--------------------------------------
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.
      ELSE
         STOP 'Wrong keyword in call to L1EQNS'
      ENDIF
C-------------------------------------------
C S(AB,IJ) = W(AB,IJ)
C-------------------------------------------
      CALL GETVVOO (S2)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * AZ(KL,IJ)
C---------------------------------------------
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF1)
      ELSE
         CALL XCOPY (NDIMT2,T2,1,BUF1,1)
      ENDIF
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NVVT(IRP)
         N = NOOT(IRP)
         K = NOOT(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,AZ(OFF2),K,
     &               A1,S2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C----------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PIJ [ T2(AB,IK) * GO(K,J) ]
C----------------------------------------------------------------
      CALL SRT4 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NVVOT(IRP)
         N = NO(IRP)
         K = NO(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,GO(OFF2),K,
     &               A0,BUF2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT10 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF2,BUF1)
      CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)
C-------------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PIJ [ W(AB,CJ) * T1(C,I) ]
C-------------------------------------------------------------------------
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)
  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3)
         IF (CARITH) CALL CONJUGA (NVO(IRP)*NINT,BUF3,1)
         IJ = 0
         DO 20 JREP = 1, NREP
           CREP = MULTB(JREP,IRP+NREP,2)
           IREP = CREP
           IF (IREP.LT.JREP) GOTO 20
           DO J = 1, NO(JREP)
              IMIN = 1
              IF (JREP.EQ.IREP) IMIN = J + 1
              DO I = IMIN, NO(IREP)
                 IJ = IJ + 1
                 ABIJ = (IVVOOTT(IRP)+(IJ-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO C = 1, NV(CREP)
                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = T1(T1CI)
                    IF (CARITH) T1VAL(2) = T1(T1CI+1)
                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABIJ),1)
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
         JI = 0
         DO 21 IREP = 1, NREP
           JREP = MULTB(IREP,IRP+NREP,2)
           CREP = IREP
           IF (JREP.LT.IREP) GOTO 21
           DO I = 1, NO(IREP)
              JMIN = 1
              IF (JREP.EQ.IREP) JMIN = I + 1
              DO J = JMIN, NO(JREP)
                 JI = JI + 1
                 ABJI = (IVVOOTT(IRP)+(JI-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO C = 1, NV(CREP)
                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = - T1(T1CI)
                    IF (CARITH) T1VAL(2) = - T1(T1CI+1)
                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABJI),1)
                 ENDDO
              ENDDO
           ENDDO
 21      CONTINUE
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30    CONTINUE
C--------------------------------
C Terms specific for T2 equations
C--------------------------------
      IF (TEQ) THEN
C----------------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [T1(A,K)*W(BK,IC)*T1(C,J)]
C----------------------------------------------------------------------------
         CALL GETVOVO (BUF1)
         CALL SRT19 (NREP,MULTB,NVO,NV,NO,NVOO,KVOOV,KKVOO,
     &               BUF1,BUF2)
         OFF1 = 1
         OFF2 = 1
         OFF3 = 1
         DO IRP =1 , NREP
            M = NVOO(IRP)
            N = NO(IRP)
            K = NV(IRP)
            CALL XGEMM ('N','N',M,N,K,A1,BUF2(OFF1),M,T1(OFF2),K,
     &                  A0,BUF1(OFF3),M)
            OFF1 = OFF1 + M * K * RCW
            OFF2 = OFF2 + K * N * RCW
            OFF3 = OFF3 + M * N * RCW
         ENDDO
         CALL SRT10 (NREP,MULTB,NVO,NO,NO,NVOO,KVOOO,KKVOO,BUF1,BUF3)
C-------------------------------------------------------------------
C        Integrals are now in VO,OOT format : combine with next term
C-------------------------------------------------------------------
      ENDIF
C----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB [ GV(A,C) * T2(CB,IJ) - T1(A,K) *W(KB,IJ) ]
C----------------------------------------------------------------------
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NV(IRP)
         N = NVOOT(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,GV(OFF1),M,BUF2(OFF2),K,
     &               A0,BUF1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL GETVOOO (BUF2)
      IF (TEQ) CALL XAXPY (NV2,-A1,BUF3,1,BUF2,1)
      CALL SRT11 (NREP,MULTB,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF3)
      CALL SRT13 (NREP,MULTB,NVO,NO,NV,NOOT,LOVOOT,LLVOOT,BUF3,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NV(IRP)
         N = NVOOT(IRP)
         K = NO(IRP)
         CALL XGEMM ('N','N',M,N,K,-A1,T1(OFF1),M,BUF2(OFF2),K,
     &               A1,BUF1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT12 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,BUF1,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
C----------------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [H(AI,CK)*T2(CK,BJ)]
C----------------------------------------------------------------------------
      CALL SRT17 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,T2,BUF1,
     &            T2,BUF1)
      CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF2,1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         JRP = MULTB(IRP+NREP,1+NREP,2)
         M = MVO(IRP)
         N = MVO(JRP)
         K = MVO(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,H(OFF1),M,BUF1(OFF2),K,
     &               A0,BUF2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT18 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,BUF2,BUF1,
     &            BUF2,BUF1)
      CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)
C------------------------------------
C Terms specific for Lamba2 equations
C------------------------------------
      IF (.NOT.TEQ) THEN
C-----------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ. [T(A,I)*H(B,J)]
C-----------------------------------------------
      CALL LAMBD21 (T1,HOV,S2)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB. [G(A,C)*W(CB,IJ)]
C---------------------------------------------
      CALL GETVVOO (BUF1)
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,BUF1,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NVOOT(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,GVV(OFF1),M,BUF2(OFF2),K,
     &               A0,BUF1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT12 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,BUF1,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PIJ. [W(AB,IK)*G(K,J)]
C---------------------------------------------
      CALL GETVVOO(BUF2)
      CALL SRT4 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF2,BUF1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NVVOT(IRP)
         N = NO(IRP)
         K = NO(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,GOO(OFF2),K,
     &               A0,BUF2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT10 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF2,BUF1)
      CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)
C
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T3CORR (FVO,EPS,T1,T2,T2S,T2T,BUF1,BUF2,BUF3,BUF4,
     &                   VVVO,VVOO,VOOO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate non-iterative T3 corrections DE(+T), DE((T)), DE(-T)
C     Use formulation of Deegan and Knowles, CPL 227 (1994), 321.
C     Note that their W is called in the description because we use
C     W for the anti-symmetrized integrals.
C     This version keeps everything in core.
C
C---------------Routines called----------------------------------------
C
C     BLAS routines :                 XGEMM, XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*),EPS(*)
      REAL*8 T1(*),T2(*),T2S(*),T2T(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*),BUF4(*)
      REAL*8 VVVO(*)
      REAL*8 VOOO(*),VVOO(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "files.inc"
      INCLUDE "symm.inc"
      INCLUDE "inpt.inc"
      INCLUDE "results.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPSIJK,EIJKABC,VMIN(2)
      REAL*8 CPUDELTA,CPU0,CPU1,CPU2,CPU3,CPUTOT
C
C---------------Executable code--------------------------------------
C
C     Initialize the timing.
C
      CPU0 = CPUDELTA()
      CALL CPUUSED(CPU0)
      CPU1 = A0
      CPU2 = A0
      CPU3 = A0
C

      ET1 = A0
      ET2 = A0
      ET3 = A0
C
      CALL GETVOOO(VOOO)
      CALL GETVVOO(VVOO)
C
C     Bring arrays in correct order
C
      CALL SRT21 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,T2,T2S)
      CALL SRT22 (NREP,MULTB,NVVT,NO,NO,NVVOT,LOVVOT,LLOVVT,T2,T2T)
      CALL SRT23 (NREP,MULTB,NVO,NV,NO,NVVT,NVVVT,KVVVOT,LLVVVT,
     &            BUF4,VVVO)
      IF (CARITH) CALL CONJUGA (KVVVOT(NREP+1),VVVO,1)
C
C     Put TAU3 in block BUF4 for use in evaluating -T correction
C
      CALL GETTAU3 (T1,T2,BUF4)
C
      CPU1 = CPUDELTA()
C
      DO 300 IJKREP = 1, NREP
      DO 200 KREP = 1, NREP
      IJREP = MULTB(KREP,IJKREP,2)
      DO 100 JREP = KREP, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 100
      IKREP = MULTB(IREP,KREP,1)
      JKREP = MULTB(JREP,KREP,1)
      DO K = 1, NO(KREP)
      KK = IO(KREP) + K
      JMIN = 1
      IF (KREP.EQ.JREP) JMIN = K + 1
      DO J = JMIN, NO(JREP)
      JJ = IO(JREP) + J
      IMIN = 1
      IF (JREP.EQ.IREP) IMIN = J + 1
      DO I = IMIN, NO(IREP)
      II = IO(IREP) + I
      EPSIJK = - EPS(II) - EPS(JJ) - EPS(KK)
C
      IJ = IIOOT(IREP,JREP) + (J-1) * NO(IREP) + I
      IK = IIOOT(IREP,KREP) + (K-1) * NO(IREP) + I
      JK = IIOOT(JREP,KREP) + (K-1) * NO(JREP) + J
      IF (IREP.EQ.JREP) IJ = IJ - J*(J+1)/2
      IF (IREP.EQ.KREP) IK = IK - K*(K+1)/2
      IF (JREP.EQ.KREP) JK = JK - K*(K+1)/2
      IJOFF1 = IVVOOT(IJREP) + (IJ-1) * NVV(IJREP)
      IKOFF1 = IVVOOT(IKREP) + (IK-1) * NVV(IKREP)
      JKOFF1 = IVVOOT(JKREP) + (JK-1) * NVV(JKREP)
      IJOFF2 = IVOOOT(IJREP) + (IJ-1) * NVO(IJREP)
      IKOFF2 = IVOOOT(IKREP) + (IK-1) * NVO(IKREP)
      JKOFF2 = IVOOOT(JKREP) + (JK-1) * NVO(JKREP)
      IJOFF3 = IVVOOTT(IJREP) + (IJ-1) * NVVT(IJREP)
      IKOFF3 = IVVOOTT(IKREP) + (IK-1) * NVVT(IKREP)
      JKOFF3 = IVVOOTT(JKREP) + (JK-1) * NVVT(JKREP)
      IOFF1  = KVVVOT(IREP) + (I-1) * NVVVT(IREP)
      JOFF1  = KVVVOT(JREP) + (J-1) * NVVVT(JREP)
      KOFF1  = KVVVOT(KREP) + (K-1) * NVVVT(KREP)
      IOFF2  = LOVVOT(IREP) + (I-1) * NOVVT(IREP)
      JOFF2  = LOVVOT(JREP) + (J-1) * NOVVT(JREP)
      KOFF2  = LOVVOT(KREP) + (K-1) * NOVVT(KREP)
C
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF1,1)
      OFF1 = 1
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      M = NV(AREP)
      N = NVVT(BCREP)
C-------------------------------------------
C U'(ABC) = T(A,E:IJ) * W(E,BC:K)
C-------------------------------------------
      EREP = MULTB(KREP,BCREP+NREP,2)
      O = NV(EREP)
      AEIJ = (IJOFF1 + IIVV(AREP,EREP)) * RCW + 1
      EBCK = (KOFF1 + LLVVVT(EREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,A1,T2S(AEIJ),M,
     &            VVVO(EBCK),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) + T(A,E:JK) * W(E,BC:I)
C-------------------------------------------
      EREP = MULTB(IREP,BCREP+NREP,2)
      O = NV(EREP)
      AEJK = (JKOFF1 + IIVV(AREP,EREP)) * RCW + 1
      EBCI = (IOFF1 + LLVVVT(EREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,A1,T2S(AEJK),M,
     &            VVVO(EBCI),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) - T(A,E:IK) * W(E,BC:J)
C-------------------------------------------
      EREP = MULTB(JREP,BCREP+NREP,2)
      O = NV(EREP)
      AEIK = (IKOFF1 + IIVV(AREP,EREP)) * RCW + 1
      EBCJ = (JOFF1 + LLVVVT(EREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,-A1,T2S(AEIK),M,
     &            VVVO(EBCJ),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) + W(A,L:IJ) * T(L,BC:K)
C-------------------------------------------
      LREP = MULTB(KREP,BCREP+NREP,2)
      O = NO(LREP)
      ALIJ = (IJOFF2 + IIVO(AREP,LREP)) * RCW + 1
      LBCK = (KOFF2 + LLOVVT(LREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,A1,VOOO(ALIJ),M,
     &            T2T(LBCK),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) + W(A,L:JK) * T(L,BC:I)
C-------------------------------------------
      LREP = MULTB(IREP,BCREP+NREP,2)
      O = NO(LREP)
      ALJK = (JKOFF2 + IIVO(AREP,LREP)) * RCW + 1
      LBCI = (IOFF2 + LLOVVT(LREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,A1,VOOO(ALJK),M,
     &            T2T(LBCI),O,A1,BUF1(OFF1),M)
C-------------------------------------------
C U'(ABC) = U'(ABC) - W(A,L:IK) * T(L,BC:J)
C-------------------------------------------
      LREP = MULTB(JREP,BCREP+NREP,2)
      O = NO(LREP)
      ALIK = (IKOFF2 + IIVO(AREP,LREP)) * RCW + 1
      LBCJ = (JOFF2 + LLOVVT(LREP,BCREP)) * RCW + 1
      CALL XGEMM ('N','N',M,N,O,-A1,VOOO(ALIK),M,
     &            T2T(LBCJ),O,A1,BUF1(OFF1),M)
C
      OFF1 = OFF1 + M * N * RCW
      ENDDO
C
      CPU2 = CPU2 + CPUDELTA()
C
C-------------------------------------------------------------------------
C U(ABC) = U'(ABC) - U'(BAC) + U'(CAB)
C-------------------------------------------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 10 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 10
      ABREP = MULTB(AREP,BREP,1)
      ACREP = MULTB(AREP,CREP,1)
      DO C = 1, NV(CREP)
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      AB = IIVVT(AREP,BREP) + (B-1) * NV(AREP) + AMIN
      AC = IIVVT(AREP,CREP) + (C-1) * NV(AREP) + AMIN
      IF (AREP.EQ.BREP) AB = AB - B*(B+1)/2
      IF (AREP.EQ.CREP) AC = AC - C*(C+1)/2
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC = (IVVVT(BCREP,IJKREP) + (BC-1)*NV(AREP) + AMIN-1) * RCW + 1
      BAC = (IVVVT(ACREP,IJKREP) + (AC-1)*NV(BREP) + B - 1) * RCW + 1
      CAB = (IVVVT(ABREP,IJKREP) + (AB-1)*NV(CREP) + C - 1) * RCW + 1
      DO A = AMIN, NV(AREP)
       BUF3(ABC) = BUF1(ABC) - BUF1(BAC) + BUF1(CAB)
       IF (CARITH) BUF3(ABC+1) = BUF1(ABC+1) - BUF1(BAC+1) + BUF1(CAB+1)
       ABC = ABC + RCW
       BAC = BAC + NV(BREP) * RCW
       CAB = CAB + NV(CREP) * RCW
      ENDDO
      ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C-------------------------------------------------------------------------
C Y2'(ABC) = V(AB:IJ) * T(C:K) + T(AB:IJ) * F(C:K)
C          + V(AB:JK) * T(C:I) + T(AB:JK) * F(C:I)
C          - V(AB:IK) * T(C:J) - T(AB:IK) * F(C:J)
C-------------------------------------------------------------------------
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF2,1)
      ABIJ = IJOFF3 * RCW + 1
      ABJK = JKOFF3 * RCW + 1
      ABIK = IKOFF3 * RCW + 1
      DO C = 1, NV(KREP)
         CK = (IVO(KREP) + (K-1) * NV(KREP) + C - 1) * RCW + 1
         ABC = (JVVVT(KREP,IJKREP)+(C-1)*NVVT(IJREP)) * RCW + 1
         CALL XAXPY (NVVT(IJREP),T1 (CK),VVOO(ABIJ),1,BUF2(ABC),1)
         CALL XAXPY (NVVT(IJREP),FVO(CK),T2  (ABIJ),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(IREP)
         CI = (IVO(IREP) + (I-1) * NV(IREP) + C - 1) * RCW + 1
         ABC = (JVVVT(IREP,IJKREP)+(C-1)*NVVT(JKREP)) * RCW + 1
         CALL XAXPY (NVVT(JKREP),T1 (CI),VVOO(ABJK),1,BUF2(ABC),1)
         CALL XAXPY (NVVT(JKREP),FVO(CI),T2  (ABJK),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(JREP)
         CJ = (IVO(JREP) + (J-1) * NV(JREP) + C - 1) * RCW + 1
         ABC = (JVVVT(JREP,IJKREP)+(C-1)*NVVT(IKREP)) * RCW + 1
         VMIN(1) = - T1(CJ)
         IF (CARITH) VMIN(2) = - T1(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,VVOO(ABIK),1,BUF2(ABC),1)
         VMIN(1) = - FVO(CJ)
         IF (CARITH) VMIN(2) = - FVO(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,T2  (ABIK),1,BUF2(ABC),1)
      ENDDO
C-------------------------------------------------------------------------
C Y2(ABC) = Y2'(ABC) + Y2'(BCA) - Y2'(ACB)
C-------------------------------------------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 20 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 20
      ABREP = MULTB(AREP,BREP,1)
      ACREP = MULTB(AREP,CREP,1)
      DO C = 1, NV(CREP)
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      AB = IIVVT(AREP,BREP) + (B-1) * NV(AREP) + AMIN
      AC = IIVVT(AREP,CREP) + (C-1) * NV(AREP) + AMIN
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (AREP.EQ.BREP) AB = AB - B*(B+1)/2
      IF (AREP.EQ.CREP) AC = AC - C*(C+1)/2
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC1 = (IVVVT(BCREP,IJKREP) + (BC-1)*NV(AREP)+AMIN-1)  * RCW + 1
      ABC2 = (JVVVT(CREP,IJKREP) + (C-1)*NVVT(ABREP)+AB-1)   * RCW + 1
      BCA  = (JVVVT(AREP,IJKREP) + (AMIN-1)*NVVT(BCREP)+BC-1)* RCW + 1
      ACB  = (JVVVT(BREP,IJKREP) + (B-1)*NVVT(ACREP)+AC-1)   * RCW + 1
      DO A = AMIN, NV(AREP)
       BUF1(ABC1) = BUF2(ABC2) + BUF2(BCA) - BUF2(ACB)
       IF (CARITH) BUF1(ABC1+1) = BUF2(ABC2+1)+BUF2(BCA+1)-BUF2(ACB+1)
       ABC1 = ABC1 + RCW
       ABC2 = ABC2 + RCW
       ACB = ACB + RCW
       BCA = BCA + NVVT(BCREP) * RCW
      ENDDO
      ENDDO
      ENDDO
 20   CONTINUE
      ENDDO
C-------------------------------------------------------------------------
C DE(+T)   = - U(ABC) * U(ABC) / EPS(ABC)
C DE'((T)) = - U(ABC) * Y2(ABC) / EPS(ABC)
C-------------------------------------------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 30 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 30
      DO C = 1, NV(CREP)
      CC = IO(NREP+1) + IV(CREP) + C
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      BB = IO(NREP+1) + IV(BREP) + B
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC = (IVVVT(BCREP,IJKREP)+(BC-1)*NV(AREP)+AMIN-1) * RCW + 1
      DO A = AMIN, NV(AREP)
       AA = IO(NREP+1) + IV(AREP) + A
       EIJKABC = EPSIJK + EPS(AA) + EPS(BB) + EPS(CC)
       ET1 = ET1 - BUF3(ABC) * BUF3(ABC) / EIJKABC
       ET2 = ET2 - BUF3(ABC) * BUF1(ABC) / EIJKABC
       IF (CARITH) THEN
         ET1 = ET1 - BUF3(ABC+1) * BUF3(ABC+1) / EIJKABC
         ET2 = ET2 - BUF3(ABC+1) * BUF1(ABC+1) / EIJKABC
       ENDIF
       ABC = ABC + RCW
      ENDDO
      ENDDO
      ENDDO
 30   CONTINUE
      ENDDO
C-------------------------------------------------------------------------
C Y3'(ABC) = TAU3(AB:IJ) * T(C:K)
C          + TAU3(AB:JK) * T(C:I)
C          - TAU3(AB:IK) * T(C:J)
C-------------------------------------------------------------------------
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF2,1)
      ABIJ = IJOFF3 * RCW + 1
      ABJK = JKOFF3 * RCW + 1
      ABIK = IKOFF3 * RCW + 1
      DO C = 1, NV(KREP)
         CK = (IVO(KREP) + (K-1) * NV(KREP) + C - 1) * RCW + 1
         ABC = (JVVVT(KREP,IJKREP)+(C-1)*NVVT(IJREP)) * RCW + 1
         CALL XAXPY (NVVT(IJREP),T1 (CK),BUF4(ABIJ),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(IREP)
         CI = (IVO(IREP) + (I-1) * NV(IREP) + C - 1) * RCW + 1
         ABC = (JVVVT(IREP,IJKREP)+(C-1)*NVVT(JKREP)) * RCW + 1
         CALL XAXPY (NVVT(JKREP),T1 (CI),BUF4(ABJK),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(JREP)
         CJ = (IVO(JREP) + (J-1) * NV(JREP) + C - 1) * RCW + 1
         ABC = (JVVVT(JREP,IJKREP)+(C-1)*NVVT(IKREP)) * RCW + 1
         VMIN(1) = - T1(CJ)
         IF (CARITH) VMIN(2) = - T1(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,BUF4(ABIK),1,BUF2(ABC),1)
      ENDDO
C-------------------------------------------------------------------------
C Y3(ABC) = Y3'(ABC) + Y3'(BCA) - Y3'(ACB)
C-------------------------------------------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 40 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 40
      ABREP = MULTB(AREP,BREP,1)
      ACREP = MULTB(AREP,CREP,1)
      DO C = 1, NV(CREP)
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      AB = IIVVT(AREP,BREP) + (B-1) * NV(AREP) + AMIN
      AC = IIVVT(AREP,CREP) + (C-1) * NV(AREP) + AMIN
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (AREP.EQ.BREP) AB = AB - B*(B+1)/2
      IF (AREP.EQ.CREP) AC = AC - C*(C+1)/2
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC1 = (IVVVT(BCREP,IJKREP) + (BC-1)*NV(AREP)+AMIN-1)  * RCW + 1
      ABC2 = (JVVVT(CREP,IJKREP) + (C-1)*NVVT(ABREP)+AB-1)   * RCW + 1
      BCA  = (JVVVT(AREP,IJKREP) + (AMIN-1)*NVVT(BCREP)+BC-1)* RCW + 1
      ACB  = (JVVVT(BREP,IJKREP) + (B-1)*NVVT(ACREP)+AC-1)   * RCW + 1
      DO A = AMIN, NV(AREP)
       BUF1(ABC1) = BUF2(ABC2) + BUF2(BCA) - BUF2(ACB)
       IF (CARITH) BUF1(ABC1+1) = BUF2(ABC2+1)+BUF2(BCA+1)-BUF2(ACB+1)
       ABC1 = ABC1 + RCW
       ABC2 = ABC2 + RCW
       ACB = ACB + RCW
       BCA = BCA + NVVT(BCREP) * RCW
      ENDDO
      ENDDO
      ENDDO
 40   CONTINUE
      ENDDO
C-------------------------------------------------------------------------
C DE'(-T) = + U(ABC) * Y3(ABC)
C-------------------------------------------------------------------------
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      DO 50 CREP = 1, NREP
      BREP = MULTB(CREP,BCREP+NREP,2)
      IF (BREP.LT.CREP.OR.AREP.LT.BREP) GOTO 50
      DO C = 1, NV(CREP)
      CC = IO(NREP+1) + IV(CREP) + C
      BMIN = 1
      IF (BREP.EQ.CREP) BMIN = C + 1
      DO B = BMIN, NV(BREP)
      BB = IO(NREP+1) + IV(BREP) + B
      AMIN = 1
      IF (AREP.EQ.BREP) AMIN = B + 1
      BC = IIVVT(BREP,CREP) + (C-1) * NV(BREP) + B
      IF (BREP.EQ.CREP) BC = BC - C*(C+1)/2
      ABC = (IVVVT(BCREP,IJKREP)+(BC-1)*NV(AREP)+AMIN-1) * RCW + 1
      DO A = AMIN, NV(AREP)
       AA = IO(NREP+1) + IV(AREP) + A
       ET3 = ET3 + BUF3(ABC) * BUF1(ABC) 
       IF (CARITH) THEN
         ET3 = ET3 + BUF3(ABC+1) * BUF1(ABC+1)
       ENDIF
       ABC = ABC + RCW
      ENDDO
      ENDDO
      ENDDO
 50   CONTINUE
      ENDDO
C
      CPU3 = CPU3 + CPUDELTA()
C
      ENDDO
      ENDDO
      ENDDO
  100 CONTINUE
  200 CONTINUE
  300 CONTINUE
C
      WRITE (IW,1071)
      WRITE (IW,1072) ESCF
      WRITE (IW,1073) ECCSD
      WRITE (IW,1074) ET1
      WRITE (IW,1075) ET2
      WRITE (IW,1076) ET3
      WRITE (IW,1077) ESCF+ECCSD+ET1
      WRITE (IW,1078) ESCF+ECCSD+ET1+ET2
      WRITE (IW,1079) ESCF+ECCSD+ET1+ET3
      ETOT = ESCF+ECCSD+ET1+ET2
C
      CALL CPUUSED (CPUTOT)
      CPUTOT = CPUTOT - CPU0
      IF (TIMING)
     &WRITE (IW,7000) CPU0,CPU1,CPU2,CPU3,CPUTOT-CPU1-CPU2-CPU3,CPUTOT
C
 1071 FORMAT(//,'  Triple excitation estimates',/)
 1072 FORMAT(' SCF energy :',T40,F25.15)
 1073 FORMAT(' CCSD correlation energy :',T40,F25.15)
 1074 FORMAT(' 4th order triples correction :',T40,F25.15)
 1075 FORMAT(' 5th order triples (T) correction :',T40,F25.15)
 1076 FORMAT(' 5th order triples -T  correction :',T40,F25.15)
 1077 FORMAT(' Total CCSD+T  energy :',T40,F25.15)
 1078 FORMAT(' Total CCSD(T) energy :',T40,F25.15)
 1079 FORMAT(' Total CCSD-T  energy :',T40,F25.15)
C
 7000 FORMAT(//' Timing of routine T3CORR :'
     &/' Before T3CORR :',T30,F12.3,' seconds'
     &/' Sorting :',T30,F12.3,' seconds'
     &/' N^7 steps :',T30,F12.3,' seconds'
     &/' Additions :',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in T3CORR :',T30,F12.3,' seconds')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZCORE(T1,T2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Zeroes out core contributions (frozen spinors)
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "freeze.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      OFF = 1
      DO IRP = 1, NREP
         DO I = 1, NO(IRP)
            IF (IFROZ(I,IRP).NE.0)
     &      CALL XCOPY (NV(IRP),A0,0,T1(OFF),1)
            OFF = OFF + NV(IRP) * RCW
         ENDDO
      ENDDO
C
      OFF = 1
      DO IJRP = 1, NREP
      N = NVVT(IJRP)
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            IF (IFROZ(I,IRP).NE.0.OR.IFROZ(J,JRP).NE.0)
     &         CALL XCOPY (N,A0,0,T2(OFF),1)
            OFF = OFF + N * RCW
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIISX(ITAPL,ITAPD,NDIMT1,NDIMT2,T1O,T1N,T2O,T2N,IT,
     &                 CC,BB,ITC,BB2,RMS,NST1,NST2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     DIIS extrapolation of the t1 and t2 amplitudes
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher (based on a routine by T.J. Lee)
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T2O(*),T2N(*),T1O(*),T1N(*)
      REAL*8 CC(MAXDIM),BB(MAXDIM+1,MAXDIM+2),BB2(MAXDIM+1,MAXDIM+2)
      REAL*8 RMS
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc" 
      INCLUDE "complex.inc" 
      INCLUDE "diis.inc"
C
C---------------Local variables--------------------------------------
C 
      REAL*8 XM,XMAX,XEN,XFAC,XADD,DET
      REAL*8 SUM1,SUM2,SSUM,DSUM,SDOT,DDOT
      COMPLEX*16 CSUM1,ZSUM,CSUM,ZDOTC,CDOTC,CFAC
C
C---------------Executable code--------------------------------------
C
      IT = IT + 1
      IF(IT.GT.MAXDIM) IT = MAXDIM
      ITC = ITC + 1
      IF(ITC.GT.MAXDIM) ITC = 1
C
      CALL XAXPY (NDIMT1,-A1,T1N,1,T1O,1)
      CALL XAXPY (NDIMT2,-A1,T2N,1,T2O,1)
      CALL XSCAL (NDIMT1,-A1,T1O,1)
      CALL XSCAL (NDIMT2,-A1,T2O,1)
      IF (CARITH) THEN
         CSUM1 = ZSUM (NDIMT1,T1O,1)
C64B     CSUM1 = CSUM (NDIMT1,T1O,1)
         SUM2 = ZDOTC (NDIMT1,T1O,1,T1O,1)
C64B     SUM2 = CDOTC (NDIMT1,T1O,1,T1O,1)
         CSUM1 = CSUM1 + ZSUM (NDIMT2,T2O,1)
C64B     CSUM1 = CSUM1 + CSUM (NDIMT2,T2O,1)
         SUM2 = SUM2 + ZDOTC (NDIMT2,T2O,1,T2O,1)
C64B     SUM2 = SUM2 + CDOTC (NDIMT2,T2O,1,T2O,1)
         RMS = DCONJG(CSUM1)*CSUM1
C64B     RMS = CONJG(CSUM1)*CSUM1
      ELSE
         SUM1 = DSUM (NDIMT1,T1O,1)
C64B     SUM1 = SSUM (NDIMT1,T1O,1)
         SUM2 = DDOT (NDIMT1,T1O,1,T1O,1)
C64B     SUM2 = SDOT (NDIMT1,T1O,1,T1O,1)
         SUM1 = SUM1 + DSUM (NDIMT2,T2O,1)
C64B     SUM1 = SUM1 + SSUM (NDIMT2,T2O,1)
         SUM2 = SUM2 + DDOT (NDIMT2,T2O,1,T2O,1)
C64B     SUM2 = SUM2 + SDOT (NDIMT2,T2O,1,T2O,1)
         RMS = SUM1 * SUM1
      ENDIF
C
      RMS = DSQRT(DABS(SUM2-RMS))
C64B  RMS = SQRT(ABS(SUM2-RMS))
C
      SECT1 = (ITC-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RWRIT(ITAPD,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAPL,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RWRIT(ITAPD,T2O,INTOWP(NDIMT2*RCW),SECT2)
      CALL RWRIT(ITAPL,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      DO 30 N = 1,IT
C
      SECT1 = (N-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAPD,T1N,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAPD,T2N,INTOWP(NDIMT2*RCW),SECT2)
C
      IF (CARITH) THEN
         CSUM1 = ZDOTC (NDIMT1,T1N,1,T1O,1)
C64B     CSUM1 = CDOTC (NDIMT1,T1N,1,T1O,1)
         CSUM1 = CSUM1 + ZDOTC (NDIMT2,T2N,1,T2O,1)
C64B     CSUM1 = CSUM1 + CDOTC (NDIMT2,T2N,1,T2O,1)
         SUM2 = DBLE(CSUM1) + DIMAG(CSUM1)
C64B     SUM2 = REAL(CSUM1) + AIMAG(CSUM1)
      ELSE
         SUM2 = DDOT (NDIMT1,T1N,1,T1O,1)
C64B     SUM2 = SDOT (NDIMT1,T1N,1,T1O,1)
         SUM2 = SUM2 + DDOT (NDIMT2,T2N,1,T2O,1)
C64B     SUM2 = SUM2 + SDOT (NDIMT2,T2N,1,T2O,1)
      ENDIF
C
      BB2(N,ITC) = SUM2
      BB2(ITC,N) = BB2(N,ITC)
C
   30 CONTINUE
C
C  TRANSFER INTO THE BB ARRAY
C
      DO 60 N = 1,IT
      DO 70 M = 1,N-1
      BB(N,M) = BB2(N,M)
      BB(M,N) = BB2(M,N)
  70  CONTINUE
      BB(N,N) = BB2(N,N)
  60  CONTINUE
C
C  FIND THE MAXIMUM AND SCALE
C
      XM = DABS(BB(1,1))
C64B  XM =  ABS(BB(1,1))
      DO 80 N = 1,IT
      DO 90 M = 1,N
      XFAC = DABS(BB(N,M))
C64B  XFAC =  ABS(BB(N,M))
      XMAX = DMAX1(XM,XFAC)
C64B  XMAX =  MAX (XM,XFAC)
      XM = XMAX
   90 CONTINUE
   80 CONTINUE
C
      XM = AR1/XM
      DO 100 N = 1,IT
      DO 110 M = 1,IT
      BB(N,M) = BB(N,M)*XM
  110 CONTINUE
  100 CONTINUE
C
      IT1 = IT+1
      IT2 = IT+2
      DO 120 N = 1,IT
      BB(N,IT1) = -AR1
      BB(IT1,N) = -AR1
  120 CONTINUE
C
      BB(IT1,IT1) = AR0
      DO 130 N = 1,IT
      BB(N,IT2) = AR0
  130 CONTINUE
      BB(IT1,IT2) = -AR1
C
      CALL FLIN(BB,MAXDIM+1,IT1,1,DET)
      XADD = AR0
      DO 140 N = 1,IT
      CC(N) = BB(N,IT2)
      XADD = XADD + CC(N)
  140 CONTINUE
C
      CALL XCOPY(NDIMT1,A0,0,T1N,1)
      CALL XCOPY(NDIMT2,A0,0,T2N,1)
C
      DO 150 N = 1,IT
C
      SECT1 = (N-1)*(NST1+NST2) + 1
      SECT2 = SECT1 + NST1
      CALL RREAD(ITAPL,T1O,INTOWP(NDIMT1*RCW),SECT1)
      CALL RREAD(ITAPL,T2O,INTOWP(NDIMT2*RCW),SECT2)
C
      CFAC = DCMPLX(CC(N),AR0)
C64B  CFAC = CMPLX(CC(N),AR0)
      CALL XAXPY (NDIMT1,CFAC,T1O,1,T1N,1)
      CALL XAXPY (NDIMT2,CFAC,T2O,1,T2N,1)
C
  150 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCENG(T1,T2,FVO,BUF1,BUF2,ECC)        
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates CC energy. Checks result (should be real).
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
      REAL*8 T1(*),T2(*)
      REAL*8 FVO(*)
      REAL*8 BUF1(*),BUF2(*)
      REAL*8 ECC
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 CECC1,CECC2
      COMPLEX*16 ZDOTC
C64B  COMPLEX*16 CDOTC
      REAL*8 DDOT,SDOT
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CECC1 = ZDOTC (NFVO,FVO,1,T1,1)
C64B     CECC1 = CDOTC (NFVO,FVO,1,T1,1)
         IF (DABS(DIMAG(CECC1)).GT.ACCUR) 
     &      WRITE (6,1000) 'T1',DIMAG(CECC1)
C64B     IF (ABS(AIMAG(CECC1)).GT.ACCUR) 
C64B &      WRITE (6,1000) 'T1',AIMAG(CECC1)
      ELSE
         ECC = DDOT (NFVO,FVO,1,T1,1)
C64B     ECC = SDOT (NFVO,FVO,1,T1,1)
      ENDIF
C
      CALL GETVVOO (BUF1)
      CALL GETTAU  (T1,T2,BUF2)
      IF (CARITH) THEN
         CECC2 = ZDOTC (NDIMT2,BUF1,1,BUF2,1)
C64B     CECC2 = CDOTC (NDIMT2,BUF1,1,BUF2,1)
         IF (DABS(DIMAG(CECC2)).GT.ACCUR) 
     &      WRITE (6,1000) 'T2',DIMAG(CECC2)
C64B     IF (ABS(AIMAG(CECC2)).GT.ACCUR) 
C64B &      WRITE (6,1000) 'T2',AIMAG(CECC2)
         ECC = DBLE(CECC1 + CECC2)
C64B     ECC = REAL(CECC1 + CECC2)
      ELSE
         ECC = ECC +  DDOT (NDIMT2,BUF1,1,BUF2,1)
C64B     ECC = ECC +  SDOT (NDIMT2,BUF1,1,BUF2,1)
      ENDIF
C
      RETURN
1000  FORMAT (/' WARNING : Imaginary contribution to ',A2,
     &' part of CC energy :',F20.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MP2ENG(T1,T2,FVO,BUF1,EMP2)        
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates MP2 energy. Checks result (should be real).
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
      REAL*8 T1(*),T2(*)
      REAL*8 FVO(*)
      REAL*8 BUF1(*)
      REAL*8 EMP2
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 CEMP21,CEMP22
      COMPLEX*16 CDOTC,ZDOTC
      REAL*8 DDOT,ZDOT
C
C---------------Executable code--------------------------------------
C
      IF (CARITH) THEN
         CEMP21 = ZDOTC (NFVO,FVO,1,T1,1)
C64B     CEMP21 = CDOTC (NFVO,FVO,1,T1,1)
         IF (DABS(DIMAG(CEMP21)).GT.ACCUR) 
     &      WRITE (6,1000) 'T1',DIMAG(CEMP21)
C64B     IF (ABS(AIMAG(CEMP21)).GT.ACCUR) 
C64B &      WRITE (6,1000) 'T1',AIMAG(CEMP21)
      ELSE
         EMP2 = DDOT (NFVO,FVO,1,T1,1)
C64B     EMP2 = SDOT (NFVO,FVO,1,T1,1)
      ENDIF
C
      CALL GETVVOO(BUF1)
      IF (CARITH) THEN
         CEMP22 = ZDOTC (NDIMT2,BUF1,1,T2,1)
C64B     CEMP22 = CDOTC (NDIMT2,BUF1,1,T2,1)
         IF (DABS(DIMAG(CEMP22)).GT.ACCUR) 
     &      WRITE (6,1000) 'T2',DIMAG(CEMP22)
C64B     IF (ABS(AIMAG(CEMP22)).GT.ACCUR) 
C64B &      WRITE (6,1000) 'T2',AIMAG(CEMP22)
         EMP2 = DBLE(CEMP21 + CEMP22)
C64B     EMP2 = REAL(CEMP21 + CEMP22)
      ELSE
         EMP2 = EMP2 + DDOT (NDIMT2,BUF1,1,T2,1)
C64B     EMP2 = EMP2 + SDOT (NDIMT2,BUF1,1,T2,1)
      ENDIF
C
      RETURN
1000  FORMAT (/' WARNING : Imaginary contribution to ',A2,
     &' part of MP2 energy :',F20.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOM (EPS,T1,T2,CT1,CT2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Divide T1 and T2 by denominators
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T1(*),T2(*)
      COMPLEX*16 CT1(*),CT2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3
C
C---------------Executable code--------------------------------------
C
      II = 0
      AI = 0
      DO IRP = 1, NREP
      DO I = 1, NO(IRP)
         II = II + 1
         FAC1 = EPS(II)
         AOFF = IO(NREP+1) + IV(IRP)
         DO A = 1, NV(IRP)
             AA = AOFF + A
             FAC = FAC1 - EPS(AA)
             AI = AI + 1
             IF (CARITH) THEN
                CT1(AI) = CT1(AI)/FAC
             ELSE
                T1(AI) = T1(AI)/FAC
             ENDIF
         ENDDO
      ENDDO
      ENDDO
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP) 
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  ABIJ = ABIJ + 1
                  IF (CARITH) THEN
                     CT2(ABIJ) = CT2(ABIJ)/FAC
                  ELSE
                     T2(ABIJ) = T2(ABIJ)/FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1DIAG(DIAG,T1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculate T1-diagnostic
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),DIAG
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER SUM
      REAL*8 RSUM,DDOT,SDOT
      COMPLEX*16 ZDOTC,CDOTC
C
C---------------Executable code--------------------------------------
C
      SUM = 0
      DO I = 1, NREP
         SUM = SUM + NO(I)
      ENDDO
      RSUM = DBLE(SUM)
C64B  RSUM = REAL(SUM)
C
      IF (CARITH) THEN
         DIAG = DBLE(ZDOTC (NDIMT1,T1,1,T1,1))
C64B     DIAG = REAL(CDOTC (NDIMT1,T1,1,T1,1))
      ELSE
         DIAG = DDOT (NDIMT1,T1,1,T1,1)
C64B     DIAG = SDOT (NDIMT1,T1,1,T1,1)
      ENDIF
C
      DIAG = DSQRT(DIAG/RSUM)
C64B  DIAG =  SQRT(DIAG/RSUM)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LINTM(EPS,FOO,FVO,FVV,T1,S1,T2,S2,AZ,H,
     &                 HO,HV,HOV,BUF1,BUF2,BUF3,NBUF2,NBUF3)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Makes intermediates for Lambda equations and writes them to file.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.
C
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),T1(*),S1(*),T2(*),S2(*)
      REAL*8 AZ(*),H(*),HO(*),HV(*),HOV(*),BUF1(*),BUF2(*),BUF3(*)
C
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "files.inc"
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "inpt.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CPUDELTA,CPU0,CPUINT,CPU4V,CPUTOT
      real*8 ddot,sdot
      complex*16 zdotc,cdotc
C
C---------------Executable code--------------------------------------
C
C     Initialize the timing.
C
      CPU0 = CPUDELTA()
      CALL CPUUSED(CPU0)
      CPUINT = A0
      CPU4V = A0
C
C  Add extra term to H-intermediate
C
      CALL HLINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
      CPU4INT = CPU4INT + CPUDELTA()
C
C  Form BL intermediate
C
      CALL BLINTM (T1,T2,BUF1,BUF2,NBUF2,BUF3,NBUF3,S2)
      CPU4V = CPU4V + CPUDELTA()
C
      CALL CPUUSED (CPUTOT)
      CPUTOT = CPUTOT - CPU0
      IF (TIMING)  WRITE(IW,7000) CPU0,CPUINT,CPU4V,CPUTOT-CPUINT-CPU4V,
     &                            CPUTOT
C
 7000 FORMAT(//' Timing of routine LINTM :'
     &/' Before LINTM :',T30,F12.3,' seconds'
     &/' Intermediates :',T30,F12.3,' seconds'
     &/' 4 virtual contr. :',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in LINTM :',T30,F12.3,' seconds')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1EQN1 (BUF1,BUF2,BUF3,NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate part of T1 or Lambda1 equations
C     The same contraction is used in the construction of the
C     X-intermediate.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 S1(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,BUF3,BUF1)
      CALL GETVOOO (BUF3)
C     --------------------------------
C     Sort from 12,34 to 134,2 (VOO,O)
C     --------------------------------
      CALL SRT6 (NREP,MULTB,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           BUF3,BUF2)
C     ------------------------
C     V(K,L,C,I) = V*(C,I,K,L)
C     ------------------------
      IF (CARITH) CALL CONJUGA (LOVOOT(NREP+1),BUF2,1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NO(IRP)
         K = NVOOT(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),K,
     &                A1,S1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE T1EQN2 (BUF1,BUF2,BUF3,NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate part of T1 or Lambda1 equations
C     The same contraction is used in the construction of the
C     X-intermediate.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF3
      REAL*8 S1(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      REAL*8 SDOT,DDOT
      COMPLEX*16 CDOTU,ZDOTU,CT1
C
C---------------Executable code--------------------------------------
C
C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)
  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3)
         IK = 0
         DO 20 KREP = 1, NREP
           IREP = MULTB(KREP,IRP+NREP,2)
           IF (KREP.GT.IREP) GOTO 20
           DO K = 1, NO(KREP)
              IMIN = 1
              IF (KREP.EQ.IREP) IMIN = K + 1
              DO I = IMIN, NO(IREP)
                 IK = IK + 1
                 CDIK = (IVVOOTT(IRP)+(IK-1)*NVVT(IRP)+ISTART)*RCW+1
                 DO A = 1, NV(IREP)
                    AK = (IIVO(IREP,KREP)+(K-1)*NV(IREP)+A-1) * RCW + 1
                    T1AI = (IVO(IREP)+(I-1)*NV(IREP)+A-1) * RCW + 1
                    IF (CARITH) THEN
                       CT1 = ZDOTU(NINT,BUF1(CDIK),1,BUF3(AK),NVO(IRP))
                       S1(T1AI) = S1(T1AI) + DBLE(CT1)
                       S1(T1AI+1) = S1(T1AI+1) + DIMAG(CT1)
C64B                   CT1 = CDOTU(NINT,BUF1(CDIK),1,BUF3(AK),NVO(IRP))
C64B                   S1(T1AI) = S1(T1AI) + REAL(CT1)
C64B                   S1(T1AI+1) = S1(T1AI+1) + AIMAG(CT1)
                    ELSE
                       S1(T1AI) = S1(T1AI)
     &                     + DDOT(NINT,BUF1(CDIK),1,BUF3(AK),NVO(IRP))
C64B &                     + SDOT(NINT,BUF1(CDIK),1,BUF3(AK),NVO(IRP))
                    ENDIF
                 ENDDO
                 DO A = 1, NV(KREP)
                    AI = (IIVO(KREP,IREP)+(I-1)*NV(KREP)+A-1) * RCW + 1
                    T1AK = (IVO(KREP)+(K-1)*NV(KREP)+A-1) * RCW + 1
                    IF (CARITH) THEN
                       CT1 = ZDOTU(NINT,BUF1(CDIK),1,BUF3(AI),NVO(IRP))
                       S1(T1AK) = S1(T1AK) - DBLE(CT1)
                       S1(T1AK+1) = S1(T1AK+1) - DIMAG(CT1)
C64B                   CT1 = CDOTU(NINT,BUF1(CDIK),1,BUF3(AI),NVO(IRP))
C64B                   S1(T1AK) = S1(T1AK) - REAL(CT1)
C64B                   S1(T1AK+1) = S1(T1AK+1) - AIMAG(CT1)
                    ELSE
                       S1(T1AK) = S1(T1AK)
     &                     - DDOT(NINT,BUF1(CDIK),1,BUF3(AI),NVO(IRP))
C64B &                     - SDOT(NINT,BUF1(CDIK),1,BUF3(AI),NVO(IRP))
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30   CONTINUE
C
      RETURN
      END
