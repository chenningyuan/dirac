C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE AINTM(T1,T2,BUF1,BUF2,AZ)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates A intermediate
C
C---------------Routines called----------------------------------------
C
C     GETTAU
C     GETOOOO, GETVVOO, GETVOOO
C     XGEMM, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
      REAL*8 AZ(*),BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      real*8 ddot
C---------------Executable code--------------------------------------
C
C--------------------
C A(KL,IJ) = W(KL,IJ)
C--------------------
      CALL GETOOOO (AZ)
C------------------------------------------------
C A(KL,IJ) = A(KL,IJ) + PIJ ( W(KL,IC) * T(C,J) )
C------------------------------------------------
      CALL GETVOOO (BUF1)
      CALL SRT9 (NREP,MULTB,NVO,NV,NO,NOOT,NOOOT,KOOOVT,KKOOOT,
     &           BUF1,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO JRP = 1, NREP
         M = NOOOT(JRP)
         N = NO   (JRP)
         K = NV   (JRP)
         CALL XCOPY (M*N,A0,0,BUF1(OFF3),1)
         CALL XGEMM ('N','N',M,N,K,A1,BUF2(OFF1),M,T1(OFF2),K,
     &               A0,BUF1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT10 (NREP,MULTB,NOOT,NO,NO,NOOOT,KOOOOT,KKOOOT,BUF1,BUF2)
      CALL XAXPY (IOOOOTT(NREP+1),A1,BUF2,1,AZ,1)
C---------------------------------------------
C A(KL,IJ) = A(KL,IJ) +  W(KL,CD) * TAU(CD,IJ)
C---------------------------------------------
      CALL GETVVOO (BUF1)
      CALL GETTAU (T1,T2,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IJRP = 1, NREP
         M = NOOT(IJRP)
         N = NOOT(IJRP)
         K = NVVT(IJRP)
         CALL XGEMM ('C','N',M,N,K,A1,BUF1(OFF1),K,BUF2(OFF2),K,
     &                A1,AZ(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GOINTM (FVO,T1,BUF1,BUF2,GO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates GO intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*)
      REAL*8 T1(*)
      REAL*8 GO(*)
      REAL*8 BUF1(*),BUF2(*)
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
C----------------------------------
C G(K,I) = H(K,I) + F(K,C) * T(C,I)
C----------------------------------
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NO(IRP)
         N = NO(IRP)
         K = NV(IRP)
         CALL XGEMM ('C','N',M,N,K,A1,FVO(OFF1),K,T1(OFF2),K,
     &               A1,GO(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C-----------------------------------------------------
C G(KI) = G(KI) + W(KL,IC) * T(CL) = G(KI) - W*(CI,KL)
C Order W from CI,KL to CL,KI
C-----------------------------------------------------
      M = MOO(1)
      N = MVO(1)
      IF (N.EQ.0) RETURN
      CALL GETVOOO (BUF1)
      CALL SRT8 (NREP,MULTB,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &           BUF1,BUF2,BUF1,BUF2)
      CALL XGEMV ('C',N,M,-A1,BUF2,N,T1,1,A1,GO,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GVINTM(FVO,T1,BUF3,NBUF3,BUF2,NBUF2,GV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates GV intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*)
      REAL*8 T1(*)
      REAL*8 GV(*)
      REAL*8 BUF3(*),BUF2(*)
      INTEGER NBUF3,NBUF2
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
C
C---------------Executable code--------------------------------------
C
C----------------------------------
C G(A,C) = H(A,C) - T(A,K) * F(K,C)
C----------------------------------
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NV(IRP)
         K = NO(IRP)
         CALL XGEMM ('N','C',M,N,K,-A1,T1(OFF1),M,FVO(OFF2),N,
     &               A1,GV(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C------------------------------------
C G(A,C) = G(A,C) + W(AK,CD) * T(D,K)
C------------------------------------
      M = MVV(1)
      N = MVO(1)
      IF (N.EQ.0) RETURN
      IF (NBUF3.LT.JVVVO(2)) THEN
         STOP " OUT OF CORE STEP IN GV NOT IMPLEMENTED"
      ENDIF
      CALL SRT20 (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,BUF3)
      CALL XGEMV ('N',M,N,A1,BUF3,M,T1,1,A1,GV,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HINTM(T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates H intermediate
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
      REAL*8 T1(*),T2(*),H(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1VAL(2)
      LOGICAL DONE
      real*8 ddot,xsum
C
C---------------Executable code--------------------------------------
C
C---------------------------------------------------------------------
C H(AK,IC) = - W(AK,CI)
C We store h initially as H(AK,CI), ordening and minus sign after next
C step !
C---------------------------------------------------------------------
      CALL GETVOVO (H)
C----------------------------------------
C H(AK,CI) = H(AK,CI) + T(A,L) * W(LK,CI)
C----------------------------------------
      CALL GETVOOO (BUF1)
      CALL SRT14 (NREP,NVO,NOOT,BUF1,BUF2)
      CALL SRT5 (NREP,MULTB,NOOT,NO,NO,NVO,LOOVO,LLOVO,BUF2,BUF1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         M = NV(IREP)
         N = NOVO(IREP)
         K = NO(IREP)
         CALL XCOPY (M*N,A0,0,BUF3(OFF3),1)
         CALL XGEMM ('N','N',M,N,K,A1,T1(OFF1),M,BUF1(OFF2),K,
     &               A0,BUF3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT15 (NREP,MULTB,NVO,NV,NO,NVO,LVOVO,LLOVO,BUF3,BUF1)
      CALL XAXPY (IVOVO(NREP+1),-A1,H,1,BUF1,1)
C------------------------------------------------------------------
C H(AK,CI) = H(AK,CI) - W(AK,CD) * T(D,I)
C H is kept in BUF1, still ordered VOVO
C------------------------------------------------------------------
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)
  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3)
         CD = 0
         DO 20 DREP = 1, NREP
           CREP = MULTB(DREP,IRP+NREP,2)
           IF (CREP.LT.DREP) GOTO 20
           DO D = 1, NV(DREP)
              CMIN = 1
              IF (CREP.EQ.DREP) CMIN = D + 1
              DO C = CMIN, NV(CREP)
                 CD = CD + 1
                 IF (CD.GT.NINT) THEN
C                   -----------------------
C                   We need the next buffer
C                   -----------------------
                    ISTART = ISTART + NINT
                    CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3)
                    CD = 1
                 ENDIF
                 AKCD = (CD-1)*NVO(IRP)*RCW+1
                 DO I = 1, NO(DREP)
                    CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    AKCI = (IVOVO(IRP)+(CI-1)*NVO(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) = - T1(T1DI)
                    IF (CARITH) T1VAL(2) = - T1(T1DI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKCI),1)
                 ENDDO
                 DO I = 1, NO(CREP)
                    DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    AKDI = (IVOVO(IRP)+(DI-1)*NVO(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = T1(T1CI)
                    IF (CARITH) T1VAL(2) = T1(T1CI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKDI),1)
                 ENDDO
              ENDDO
           ENDDO
 20      CONTINUE
 30    CONTINUE
C------------------------------------------------------------------
C Now sort H(AK,CI) to H(AI,CK)
C------------------------------------------------------------------
      CALL SRT16 (NREP,MULTB,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            BUF1,H,BUF1,H)
C------------------------------------------------------------------
C H(AI,CK) = H(AI,CK) + [ 0.5*T(AD,IL) - T(A,L)*T(D,I) ] * W(KL,CD)
C------------------------------------------------------------------
      CALL XCOPY (NDIMT2,T2,1,BUF1,1)
      CALL XSCAL (NDIMT2,AP5,BUF1,1)
      CALL GETTAU (T1,BUF1,BUF2)
      CALL SRT17 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &            BUF2,BUF1,BUF2,BUF1)
C---------------------------------------------------------------
C BUF1 now contains 0.5*T(AD,IL) - T(A,L)*T(D,I) + T(A,I)*T(D,L)
C Subtract T(A,I)*T(D,L)
C---------------------------------------------------------------
      CALL XGEMM ('N','N',NDIMT1,NDIMT1,1,-A1,T1,NDIMT1,T1,1,
     &               A1,BUF1,NDIMT1)
C
      CALL GETVVOO (BUF3)
      IF (CARITH) CALL CONJUGA (NV3,BUF3,1)
      CALL SRT17 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &            BUF3,BUF2,BUF3,BUF2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),K,
     &               A1,H(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HOINTM(FOO,T1,T2,BUF1,BUF2,BUF3,HO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates HO intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FOO(*)
      REAL*8 T1(*),T2(*)
      REAL*8 HO(*)
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
C-----------------------------
C H(K,I) = F(K,I) * DELTA(K,I)
C-----------------------------
      CALL XCOPY (NFOO,FOO,1,HO,1)
      II = 1
      DO IRP = 1, NREP
         N = NO(IRP)
         CALL XCOPY(N,A0,0,HO(II),N+1)
         II = II + N * N * RCW
      ENDDO
C----------------------------------------------------------------------
C H(K,I) = H(K,I) + W(KL,CD) * T(CD,IL) = H(K,I) + W*(CDL,K) * T(CDL,I)
C----------------------------------------------------------------------
      CALL GETVVOO (BUF3)
      CALL SRT4 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF3,BUF1)
      CALL GETTAU  (T1,T2,BUF3)
      CALL SRT4 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF3,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NO(IRP)
         N = NO(IRP)
         K = NVVOT(IRP)
         CALL XGEMM ('C','N',M,N,K,A1,BUF1(OFF1),K,BUF2(OFF2),K,
     &               A1,HO(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HOVINTM(FVO,T1,BUF1,BUF2,HOV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates HOV intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*)
      REAL*8 T1(*)
      REAL*8 HOV(*)
      REAL*8 BUF1(*),BUF2(*)
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
C--------------------------------
C HOV(C,K) = FOV(K,C) = FVO*(C,K)
C--------------------------------
      CALL XCOPY (NFVO,FVO,1,HOV,1)
      IF (CARITH) CALL CONJUGA (NFVO,HOV,1)
C-------------------------------------------------------------------
C HOV(CK) = HOV(CK) + W(KL,CD) * T(DL) = HOV(CK) + W*(CD,KL) * T(DL)
C-------------------------------------------------------------------
      M = MVO(1)
      K = MVO(1)
      CALL GETVVOO (BUF1)
      CALL XCOPY (JVOVO(2),A0,0,BUF2,1)
      CALL SRT1 (NREP,MULTB,NV,NV,NO,NO,MVO,JVOVO,JJVO,JJVO,
     &           BUF1,BUF2,BUF1,BUF2)
      IF (CARITH) CALL CONJUGA (JVOVO(2),BUF2,1)
      CALL XGEMV ('N',M,K,A1,BUF2,M,T1,1,A1,HOV,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HVINTM(FVV,T1,T2,BUF1,BUF2,BUF3,HV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates HV intermediate
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVV(*)
      REAL*8 T1(*),T2(*)
      REAL*8 HV(*)
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
C-----------------------------
C H(A,C) = F(A,C) * DELTA(A,C)
C-----------------------------
      CALL XCOPY (NFVV,FVV,1,HV,1)
      AA = 1
      DO IRP = 1, NREP
         N = NV(IRP)
         CALL XCOPY(N,A0,0,HV(AA),N+1)
         AA = AA + N * N * RCW
      ENDDO
C----------------------------------------------------------------------
C H(A,C) = H(A,C) - W(KL,CD) * T(AD,KL) = H(A,C) - T(A,DKL) * W*(C,DKL)
C----------------------------------------------------------------------
      CALL GETTAU (T1,T2,BUF3)
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,BUF3,BUF1)
      CALL GETVVOO (BUF3)
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,BUF3,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NV(IRP)
         K = NVOOT(IRP)
         CALL XGEMM ('N','C',M,N,K,-A1,BUF1(OFF1),M,BUF2(OFF2),N,
     &               A1,HV(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BINTM(T1,T2,L1,BUF1,BUF2,NBUF2,BUF3,NBUF3,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates B intermediate and adds it to new T2's
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF2,NBUF3
      REAL*8 T1(*),T2(*),L1(*)
      REAL*8 S2(*)
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
      REAL*8 T1MIN(2)
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
C-------------------------------------------------
C B(AB,CD) =  W(AB,CD) - PAB [ T(B,K) * W(AK,CD) ]
C-------------------------------------------------
      IF (TEQ) THEN
         CALL GETTAU (T1,T2,BUF1)
      ELSE
         CALL GETAMPT (L1,BUF1)
         CALL XCOPY   (NDIMT2,T2,1,BUF1,1)
         CALL LAMBD21 (T1,L1,BUF1)
      ENDIF
      OFF1 = 1
      OFF2 = 1
      DO 30 IRP = 1, NREP
         ISTART = 0
         IF (NVVT(IRP).EQ.0) GOTO 30
         NINT = NBUF3/NVVT(IRP)
         IF (TEQ.AND.NVO(IRP).NE.0) NINT = MIN0(NINT,NBUF2/NVO(IRP))
  10     CALL GETVVVV (IRP,ISTART,NINT,DONE,BUF3)
C------------------------------------------------------------------
C We only add this contribution for the amplitude equations, in the
C lambda equations we work with modified integrals
C------------------------------------------------------------------
         IF ((.NOT.TEQ).OR.NVO(IRP).EQ.0) GOTO 21
C-------------------------------------------------------------------
C Make this batch of B : add 3-virtual integrals contracted with T1
C-------------------------------------------------------------------
         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF2)
         AK = 1
         DO 20 KREP = 1, NREP
           BREP = KREP
           AREP = MULTB(KREP,IRP+NREP,2)
           DO K = 1, NO(KREP)
              DO A = 1, NV(AREP)
                 BK = (IVO(KREP)+(K-1)*NV(BREP)) * RCW + 1
                 IF (AREP.GT.BREP) THEN
                    AB = (IIVVT(AREP,BREP)+A-1) * RCW + 1
                    DO B = 1, NV(BREP)
                       T1MIN(1) = - T1(BK)
                       IF (CARITH) T1MIN(2) = - T1(BK+1)
                       CALL XAXPY(NINT,T1MIN,BUF2(AK),NVO(IRP),
     &                            BUF3(AB),NVVT(IRP))
                       AB = AB + NV(AREP) * RCW
                       BK = BK + RCW
                    ENDDO
                 ELSEIF (AREP.EQ.BREP) THEN
                    AB = (IIVVT(AREP,BREP)+A-2) * RCW + 1
                    DO B = 1, A - 1
                       T1MIN(1) = - T1(BK)
                       IF (CARITH) T1MIN(2) = - T1(BK+1)
                       CALL XAXPY(NINT,T1MIN,BUF2(AK),NVO(IRP),
     &                            BUF3(AB),NVVT(IRP))
                       AB = AB + (NV(BREP)-B-1) * RCW
                       BK = BK + RCW
                    ENDDO
                    BK = BK + RCW
                    BA = (IIVVT(BREP,AREP)+(A-1)*NV(BREP)-
     &                    A*(A-1)/2) * RCW + 1
                    DO B =  A + 1, NV(BREP)
                       CALL XAXPY(NINT,T1(BK),BUF2(AK),NVO(IRP),
     &                            BUF3(BA),NVVT(IRP))
                       BA = BA + RCW
                       BK = BK + RCW
                    ENDDO
                 ELSEIF (AREP.LT.BREP) THEN
                    BA = (IIVVT(BREP,AREP)+(A-1)*NV(BREP)) * RCW + 1
                    DO B = 1, NV(BREP)
                       CALL XAXPY(NINT,T1(BK),BUF2(AK),NVO(IRP),
     &                            BUF3(BA),NVVT(IRP))
                       BA = BA + RCW
                       BK = BK + RCW
                    ENDDO
                 ENDIF
                 AK = AK + RCW
              ENDDO
           ENDDO
 20      CONTINUE
 21      CONTINUE
C----------------------------------------
C Calculate contribution of B to new T2's
C----------------------------------------
         M = NVVT(IRP)
         N = NOOT(IRP)
         K = NINT
         CALL XGEMM ('N','N',M,N,K,A1,BUF3,M,BUF1(OFF1+ISTART*RCW),M,
     &               A1,S2(OFF2),M)
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
         OFF1 = OFF1 + M * N * RCW
         OFF2 = OFF2 + M * N * RCW
 30   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETTAU (T1,T2,TAU)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates Tau from T1 and T2 amplitudes
C
C---------------Routines called----------------------------------------
C
C     XCOPY, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),TAU(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1MIN(2)
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT2,T2,1,TAU,1)
C
      ABIJOFF = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) + T1(A,I) * T1(B,J)
C--------------------------------------------------------
            IF (BRP.EQ.JRP) THEN
               AB = 1
               BJ = (IVO(JRP) + (J-1) * NV(JRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AI = (IVO(IRP) + (I-1) * NV(IRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  CALL XAXPY (NA,T1(BJ),T1(AI),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BJ = BJ + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) - T1(A,J) * T1(B,I)
C--------------------------------------------------------
            IF (ARP.EQ.JRP) THEN
               AB = 1
               BI = (IVO(IRP) + (I-1) * NV(IRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AJ = (IVO(JRP) + (J-1) * NV(JRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1MIN(1) = -T1(BI)
                  IF (CARITH) T1MIN(2) = -T1(BI+1)
                  CALL XAXPY (NA,T1MIN,T1(AJ),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BI = BI + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (ARP.NE.BRP) THEN
               ABIJOFF = ABIJOFF + NV(ARP) * NV(BRP) * RCW
            ELSE
               ABIJOFF = ABIJOFF + NV(ARP) * (NV(ARP)-1) * RCW / 2
            ENDIF
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETTAU3 (T1,T2,TAU)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates Tau3 from T1 and T2 amplitudes
C     Tau3 is used in the -T correction : T3 = T2 + 1/3 PIJ T1.T1
C
C---------------Routines called----------------------------------------
C
C     XCOPY, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*),TAU(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
      INCLUDE "param.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1VAL(2)
C
C---------------Executable code--------------------------------------
C
      CALL XCOPY (NDIMT2,T2,1,TAU,1)
C
      ABIJOFF = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) + T1(A,I) * T1(B,J)
C--------------------------------------------------------
            IF (BRP.EQ.JRP) THEN
               AB = 1
               BJ = (IVO(JRP) + (J-1) * NV(JRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AI = (IVO(IRP) + (I-1) * NV(IRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1VAL(1) = T1(BJ) * ARP3
                  IF (CARITH) T1VAL(2) = T1(BJ+1) * ARP3
                  CALL XAXPY (NA,T1VAL,T1(AI),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BJ = BJ + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C TAU(AB,IJ) = T2(AB,IJ) - T1(A,J) * T1(B,I)
C--------------------------------------------------------
            IF (ARP.EQ.JRP) THEN
               AB = 1
               BI = (IVO(IRP) + (I-1) * NV(IRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AJ = (IVO(JRP) + (J-1) * NV(JRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1VAL(1) = -T1(BI) * ARP3
                  IF (CARITH) T1VAL(2) = -T1(BI+1) * ARP3
                  CALL XAXPY (NA,T1VAL,T1(AJ),1,TAU(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BI = BI + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (ARP.NE.BRP) THEN
               ABIJOFF = ABIJOFF + NV(ARP) * NV(BRP) * RCW
            ELSE
               ABIJOFF = ABIJOFF + NV(ARP) * (NV(ARP)-1) * RCW / 2
            ENDIF
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LAMBD21(T1,L1,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates contribution to L2 equations
C
C---------------Routines called----------------------------------------
C
C     XCOPY, XAXPY
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 L1(*),T1(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "symm.inc"
      INCLUDE "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1MIN(2)
C
C---------------Executable code--------------------------------------
C
      ABIJOFF = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
C--------------------------------------------------------------
C L2(AB,IJ) = L2(AB,IJ) + L1(A,I) * T1(B,J) + T1(A,I) * L1(B,J)
C--------------------------------------------------------------
            IF (BRP.EQ.JRP) THEN
               AB = 1
               BJ = (IVO(JRP) + (J-1) * NV(JRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AI = (IVO(IRP) + (I-1) * NV(IRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  CALL XAXPY (NA,T1(BJ),L1(AI),1,S2(ABIJOFF+AB),1)
                  CALL XAXPY (NA,L1(BJ),T1(AI),1,S2(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BJ = BJ + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------------
C S2(AB,IJ) = S2(AB,IJ) - L1(A,J) * T1(B,I) - T1(A,J) * L1(B,I)
C--------------------------------------------------------------
            IF (ARP.EQ.JRP) THEN
               AB = 1
               BI = (IVO(IRP) + (I-1) * NV(IRP)) * RCW + 1
               DO B = 1, NV(BRP)
                  AMIN = 1
                  IF (ARP.EQ.BRP) AMIN = B + 1
                  AJ = (IVO(JRP) + (J-1) * NV(JRP) + AMIN - 1) * RCW + 1
                  NA = NV(ARP) - AMIN + 1
                  T1MIN(1) = -T1(BI)
                  IF (CARITH) T1MIN(2) = -T1(BI+1)
                  CALL XAXPY (NA,T1MIN,L1(AJ),1,S2(ABIJOFF+AB),1)
                  T1MIN(1) = -L1(BI)
                  IF (CARITH) T1MIN(2) = -L1(BI+1)
                  CALL XAXPY (NA,T1MIN,T1(AJ),1,S2(ABIJOFF+AB),1)
                  AB = AB + NA * RCW
                  BI = BI + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (ARP.NE.BRP) THEN
               ABIJOFF = ABIJOFF + NV(ARP) * NV(BRP) * RCW
            ELSE
               ABIJOFF = ABIJOFF + NV(ARP) * (NV(ARP)-1) * RCW / 2
            ENDIF
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE HLINTM(T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Subtracts extra term to make H intermediate for lambda equations
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
      REAL*8 T1(*),T2(*),H(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      real*8 ddot,xsum
C
C---------------Executable code--------------------------------------
C
C------------------------------------------------------------------
C H(AI,CK) = H(AI,CK) + 0.5*T(AD,IL) * W(KL,CD)
C------------------------------------------------------------------
      CALL XCOPY (NDIMT2,T2,1,BUF2,1)
      CALL XSCAL (NDIMT2,AP5,BUF2,1)
      CALL SRT17 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &            BUF2,BUF1,BUF2,BUF1)
      CALL GETVVOO (BUF3)
      IF (CARITH) CALL CONJUGA (NV3,BUF3,1)
      CALL SRT17 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &            BUF3,BUF2,BUF3,BUF2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),K,
     &               A1,H(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLINTM(T1,T2,BUF1,BUF2,NBUF2,BUF3,NBUF3)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates B intermediate for lambda equations and write it to
C     file.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF2,NBUF3
      REAL*8 T1(*),T2(*)
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
      REAL*8 T1MIN(2)
      LOGICAL DONE
C
C---------------Executable code--------------------------------------
C
C-------------------------------------------------
C B(AB,CD) =  W(AB,CD) - PAB [ T(B,K) * W(AK,CD) ]
C-------------------------------------------------
      CALL GETVVOO (BUF2)
      CALL SRT14 (NREP,NVVT,NOOT,BUF2,BUF1)
      OFF1 = 1
      OFF2 = 1
      DO 30 IRP = 1, NREP
         ISTART = 0
         IF (NVVT(IRP).EQ.0) GOTO 30
         NINT = NBUF3/NVVT(IRP)
         IF (NVO(IRP).NE.0) NINT = MIN0(NINT,NBUF2/NVO(IRP))
  10     CALL GETVVVV (IRP,ISTART,NINT,DONE,BUF3)
         IF (NVO(IRP).EQ.0) GOTO 21
C-------------------------------------------------------------------
C Make this batch of B : add 3-virtual integrals contracted with T1
C After addition : make modified 3-virtual integrals
C-------------------------------------------------------------------
         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF2)
         AK = 1
         DO 20 KREP = 1, NREP
           BREP = KREP
           AREP = MULTB(KREP,IRP+NREP,2)
           LREP = AREP
           DO K = 1, NO(KREP)
              DO A = 1, NV(AREP)
                 BK = (IVO(KREP)+(K-1)*NV(BREP)) * RCW + 1
                 AL = (IVO(LREP)+A-1) * RCW + 1
C                --------------------------
                 IF (AREP.GT.BREP) THEN
C                --------------------------
                    AB = (IIVVT(AREP,BREP)+A-1) * RCW + 1
                    DO B = 1, NV(BREP)
                       T1MIN(1) = - T1(BK)
                       IF (CARITH) T1MIN(2) = - T1(BK+1)
                       CALL XAXPY(NINT,T1MIN,BUF2(AK),NVO(IRP),
     &                            BUF3(AB),NVVT(IRP))
                       AB = AB + NV(AREP) * RCW
                       BK = BK + RCW
                    ENDDO
C                   -----> make modified VOOO
                    LK = (IIOOT(LREP,KREP)+(K-1)*NO(LREP)) * RCW + 1
                    DO L = 1, NO(LREP)
                       T1MIN(1) = - T1(AL)
                       IF (CARITH) T1MIN(2) = - T1(AL+1)
                       CALL XAXPY(NINT,T1MIN,BUF1(LK),NOOT(IRP),
     &                            BUF2(AK),NVO(IRP))
                       AL = AL + NV(AREP) * RCW
                       LK = LK + RCW
                    ENDDO
C                --------------------------
                 ELSEIF (AREP.EQ.BREP) THEN
C                --------------------------
                    AB = (IIVVT(AREP,BREP)+A-2) * RCW + 1
                    DO B = 1, A - 1
                       T1MIN(1) = - T1(BK)
                       IF (CARITH) T1MIN(2) = - T1(BK+1)
                       CALL XAXPY(NINT,T1MIN,BUF2(AK),NVO(IRP),
     &                            BUF3(AB),NVVT(IRP))
                       AB = AB + (NV(BREP)-B-1) * RCW
                       BK = BK + RCW
                    ENDDO
                    BK = BK + RCW
                    BA = (IIVVT(BREP,AREP)+(A-1)*NV(BREP)-
     &                    A*(A-1)/2) * RCW + 1
                    DO B =  A + 1, NV(BREP)
                       CALL XAXPY(NINT,T1(BK),BUF2(AK),NVO(IRP),
     &                            BUF3(BA),NVVT(IRP))
                       BA = BA + RCW
                       BK = BK + RCW
                    ENDDO
C                   -----> make modified VOOO
                    KL = (IIOOT(KREP,LREP)+K-2) * RCW + 1
                    DO L = 1, K - 1
                       CALL XAXPY(NINT,T1(AL),BUF1(KL),NOOT(IRP),
     &                            BUF2(AK),NVO(IRP))
                       KL = KL + (NO(KREP)-L-1) * RCW
                       AL = AL + NV(AREP) * RCW
                    ENDDO
                    AL = AL + NV(AREP) * RCW
                    LK = (IIOOT(LREP,KREP)+(K-1)*NO(LREP)-
     &                    K*(K-1)/2) * RCW + 1
                    DO L = K + 1, NO(LREP)
                       T1MIN(1) = - T1(AL)
                       IF (CARITH) T1MIN(2) = - T1(AL+1)
                       CALL XAXPY(NINT,T1MIN,BUF1(LK),NOOT(IRP),
     &                            BUF2(AK),NVO(IRP))
                       AL = AL + NV(AREP) * RCW
                       LK = LK + RCW
                    ENDDO
C                --------------------------
                 ELSEIF (AREP.LT.BREP) THEN
C                --------------------------
                    BA = (IIVVT(BREP,AREP)+(A-1)*NV(BREP)) * RCW + 1
                    DO B = 1, NV(BREP)
                       CALL XAXPY(NINT,T1(BK),BUF2(AK),NVO(IRP),
     &                            BUF3(BA),NVVT(IRP))
                       BA = BA + RCW
                       BK = BK + RCW
                    ENDDO
C                   -----> make modified VOOO
                    KL = (IIOOT(KREP,LREP)+K-1) * RCW + 1
                    DO L = 1, NO(LREP)
                       CALL XAXPY(NINT,T1(AL),BUF1(KL),NOOT(IRP),
     &                            BUF2(AK),NVO(IRP))
                       AL = AL + NV(AREP) * RCW
                       KL = KL + NO(KREP) * RCW
                    ENDDO
                 ENDIF
                 AK = AK + RCW
              ENDDO
           ENDDO
 20      CONTINUE
         CALL PUTVOVV (IRP,ISTART,NINT,DONE,BUF2)
 21      CONTINUE
C-------------------------------------------
C B(AB,CD) =  B(AB,CD) - T(AB,KL) * W(KL,CD)
C-------------------------------------------
         M = NVVT(IRP)
         N = NOOT(IRP)
         K = NINT
         KLCD = OFF1+ISTART*NOOT(IRP)*RCW
         CALL XGEMM ('N','N',M,K,N,A1,T2,M,BUF1(KLCD),M,
     &               -A1,BUF3,M)
         CALL PUTVVVV (IRP,ISTART,NINT,DONE,BUF3)
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
         OFF1 = OFF1 + M * N * RCW
         OFF2 = OFF2 + M * N * RCW
 30   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MLVOOO1(T1,BUF1,BUF2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Make modified VOOO integrals for lambda equations
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*)
      REAL*8 BUF1(*),BUF2(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      real*8 ddot,xsum
      dimension noovt(8) ! DON'T FORGET TO DEFINE in SYMTAB
C
C---------------Executable code--------------------------------------
C
C-------------------------------------------------------
C LVOOO1(CI,KL) = VVOOO(CI,KL) + VVVOO (CD,KL) * T1(D,I)
C-------------------------------------------------------
      CALL GETVVOO (BUF1)
      CALL SRT14 (NREP,NVVT,NOOT,BUF1,BUF2)
      CALL SRT4  (NREP,MULTB,NOOT,NV,NV,NOOVT,KOOVVT,KKOOVT,BUF2,BUF1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NOOVT(IRP)
         N = NO(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,T1(OFF2),K,
     &               A0,BUF2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT24 (NREP,MULTB,NOOT,NV,NO,NOOOT,KOOVOT,KKOOVT,BUF2,BUF1)
      CALL SRT14 (NREP,NVVT,NVO,BUF2,BUF1)
      CALL GETVOOO (BUF2)
      CALL XAXPY (IVOOOT(NREP+1),A1,BUF1,1,BUF2,1)
cluuk CALL PUTVOOO (BUF2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MLVOOO2(T1,T2,HOV,AZ,BUF1,BUF2,BUF3,NBUF3,LVOOO2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C
C     Make modified VOOO integrals for lambda equations
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
      REAL*8 T1(*),T2(*),AZ(*),HOV(*),LVOOO2(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "param.inc"
      INCLUDE "complex.inc"
      INCLUDE "symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 T1VAL(2)
      LOGICAL DONE
      real*8 ddot,xsum
C
C---------------Executable code--------------------------------------
C
C---------------------------------------------------------------------
C LVOOO2(AL,IJ) = W(AL,IJ)
C---------------------------------------------------------------------
      CALL GETVOOO (LVOOO2)
C--------------------------------------------------
C LVOOO2(AL,IJ) = LVOOO2(AL,IJ) - T(A,K) * A(KL,IJ)
C--------------------------------------------------
      CALL SRT5 (NREP,MULTB,NOOT,NO,NO,NOOT,LOOOOT,LLOOOT,AZ,BUF1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         M = NV(IREP)
         N = NOOOT(IREP)
         K = NO(IREP)
         CALL XCOPY (M*N,A0,0,BUF2(OFF3),1)
         CALL XGEMM ('N','N',M,N,K,A1,T1(OFF1),M,BUF1(OFF2),K,
     &               A0,BUF2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT15 (NREP,MULTB,NVO,NV,NO,NOOT,LVOOOT,LLOOOT,BUF2,BUF1)
      CALL XAXPY (NV2,-A1,BUF2,1,LVOOO2,1)
C------------------------------------------------------------
C LVOOO2(AL,IJ) = LVOOO2(AL,IJ) + 0.5 * V(AL,CD) * TAU(CD,IJ)
C------------------------------------------------------------
      CALL GETTAU (T1,T2,BUF1)
      CALL XSCAL (NDIMT2,AP5,BUF1,1)
      OFF1 = 1
      OFF2 = 1
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)
  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3)
         M = NVO (IRP)
         N = NOOT(IRP)
         K = NINT
         CALL XGEMM ('N','N',M,N,K,A1,BUF3,M,BUF1(OFF1+ISTART*RCW),
     &               NVVT(IRP),A1,LVOOO2(OFF2),M)
         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
         OFF1 = OFF1 + M * N * RCW
         OFF2 = OFF2 + M * N * RCW
 30    CONTINUE
C------------------------------------------------------------------
C LVOOO2(AL,IJ) = LVOOO2(AL,IJ) + PIJ.[ T(AC,IK) * LVOOO1(KL,CJ) ]
C------------------------------------------------------------------
      CALL GETVOOO (BUF1)
      IF (CARITH) CALL CONJUGA (NV2,BUF1,1)
      CALL SRT25 (NREP,MULTB,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &            BUF1,BUF2)
      CALL SRT17 (NREP,MULTB,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &            T2,BUF1,T2,BUF1)
      OFF1 = 1
C     ---------------------------------------------------------
C     Ordered as T(AI,CK) and V(CK,JL) : do the multiplication.
C     ---------------------------------------------------------
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MOO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOOO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),K,
     &               A0,BUF3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C     -----------------------------------------------------------------
C     Sort from order (AI,JL) back to (AL,JI) and add with a minus sign
C     -----------------------------------------------------------------
      CALL SRT26 (NREP,MULTB,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &            BUF3,BUF1,BUF3,BUF1)
      CALL XAXPY (NV2,-A1,BUF1,1,LVOOO2,1)
C------------------------------------------------------------------
C LVOOO2(AL,IJ) = LVOOO2(AL,IJ) + PIJ.[ V(AL,ID) * T(D,J) ]
C------------------------------------------------------------------
      CALL GETVOVO (BUF1)
      CALL SRT19 (NREP,MULTB,NVO,NV,NO,NVOO,KVOOV,KKVOO,BUF1,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NVOO(IRP)
         N = NO(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,A1,BUF2(OFF1),M,T1(OFF2),K,
     &               A0,BUF1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT10 (NREP,MULTB,NVO,NO,NO,NVOO,KVOOO,KKVOO,BUF1,BUF2)
      CALL XAXPY (NV2,A1,BUF2,1,LVOOO2,1)
C------------------------------------------------------
C LVOOO2(AL,IJ) = LVOOO2(AL,IJ) + HOV(L,C) * T(AC,IJ) ]
C------------------------------------------------------
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF1)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NO(IRP)
         N = NVOOT(IRP)
         K = NV(IRP)
         CALL XGEMM ('C','N',M,N,K,A1,HOV(OFF1),K,BUF1(OFF2),K,
     &                A0,BUF2(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
      CALL SRT15 (NREP,MULTB,NOV,NO,NV,NOOT,LOVOOT,LLVOOT,BUF2,BUF1)
      CALL SRT11 (NREP,MULTB,NOV,NO,NV,NOOT,IOVOOT,IIOV,BUF1,BUF2)
      CALL XAXPY (NV2,A1,BUF2,1,LVOOO2,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MLVOVV(T1,T2,BUF1,BUF2,NBUF2,BUF3,NBUF3)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates LVOVV intermediate for lambda equations and write it to
C     file.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF2,NBUF3
      REAL*8 T1(*),T2(*)
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
      REAL*8 T1MIN(2)
      LOGICAL DONEAB,DONECD
C
C---------------Executable code--------------------------------------
C
C----------------------------------------------------
C LVOVV(AB,CJ) =  VVOVV(AB,CJ) + BL(AB,CD) * T(D,J)
C Store this intermediate as LVOVV*(CJ,AB)
C----------------------------------------------------
      DO 40 IRP = 1, NREP
         IF (NOV(IRP).EQ.0.OR.NVVT(IRP).EQ.0) GOTO 40
         ISTART = 0
         NINTCD = NBUF3/NVVT(IRP)
 10      CALL GETVVVV (IRP,ISTART,NINTCD,DONECD,BUF3)
         JSTART = 0
         NINTAB = NBUF2/NVO(IRP)
 20      CALL GETVOVV (IRP,JSTART,NINTAB,DONEAB,BUF2)
         ABCD = JSTART*RCW + 1
         CD = 0
         DO 30 DREP = 1, NREP
            CREP = MULTB(DREP,IRP+NREP,2)
            IF (CREP.LT.DREP) GOTO 30
            JREP = DREP
            DO D = 1, NV(DREP)
               CMIN = 1
               IF (CREP.EQ.DREP) CMIN = D + 1
               DO C = CMIN, NV(CREP)
C                 -------------------------------------------------
C                 If the integralblock (**,CD) is in memory we make
C                 the contributions to (AB,CJ)
C                 -------------------------------------------------
                  IF (CD.GE.ISTART.AND.CD.LT.ISTART+NINTCD) THEN
                     DO J = 1, NO(JREP)
                        DJ = (IVO(JREP)+(J-1)*NV(DREP)+D-1) * RCW + 1
                        CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                        CALL XAXPY(NINTAB,T1(DJ),BUF3(ABCD),1,
     &                             BUF2(CJ),NVO(IRP))
                     ENDDO
                     ABCD = ABCD + NVVT(IRP) * RCW
                  ENDIF
                  CD = CD + 1
               ENDDO
            ENDDO
 30      CONTINUE
         CALL PUTVOVV (IRP,JSTART,NINT,DONE,BUF2)
         IF (.NOT.DONEAB) THEN
            JSTART = JSTART + NINTAB
            GOTO 20
         ENDIF
         IF (.NOT.DONECD) THEN
            ISTART = ISTART + NINTCD
            GOTO 10
         ENDIF
 40   CONTINUE
C
C-----------------------------------------------------
C LVOVV(AB,CJ) =  VVOVV(AB,CJ) - TAU(AB,KL) * V(KL,CJ)
C Store this intermediate as LVOVV*(CJ,AB)
C-----------------------------------------------------
      CALL GETTAU (T1,T2,BUF1)
      CALL GETVOOO (BUF2)
      OFF1 = 0
      DO 60 IRP = 1, NREP
         IF (NOV(IRP).EQ.0) GOTO 59
         JSTART = 0
         NINTAB = NBUF3/NVO(IRP)
 50      CALL GETVOVV (IRP,JSTART,NINTAB,DONEAB,BUF3)
         M = NVO(IRP)
         N = NINTAB
         K = NOOT(IRP)
         KLAB = OFF1+JSTART*RCW
         CALL XGEMM ('N','C',M,N,K,A1,BUF2,M,BUF1(KLAB),NVVT(IRP),
     &               -A1,BUF3,M)
         CALL PUTVOVV (IRP,JSTART,NINT,DONE,BUF3)
         IF (.NOT.DONEAB) THEN
            JSTART = JSTART + NINTAB
            GOTO 50
         ENDIF
 59      OFF1 = OFF1 + NVVT(IRP)*NOOT(IRP)*RCW
 60   CONTINUE
c
c   half-finished part : some terms missing !!
c
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MKGPP (T2,S2,BUF1,BUF2,GOO,GVV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Calculates GOO and GVV intermediates for lambda equations.
C     Note that the definition differs from Gauss & Stanton
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 GOO(*),GVV(*)
      REAL*8 T2(*),S2(*)
      REAL*8 BUF1(*),BUF2(*)
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
C---------------Executable code--------------------------------------
C
C----------------------------------------------------------------------
C G(I,J) = G(I,J) - L(CD,IL) * T(CD,JL) = - L(CDL,I) * T(CDL,J)
C----------------------------------------------------------------------
      CALL SRT4 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,S2,BUF1)
      CALL SRT4 (NREP,MULTB,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,T2,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NO(IRP)
         N = NO(IRP)
         K = NVVOT(IRP)
         CALL XGEMM ('C','N',M,N,K,-A1,BUF1(OFF1),K,BUF2(OFF2),K,
     &               A0,GOO(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C----------------------------------------------------------------------
C G(A,B) = L(AD,KL) * T(BD,KL) = L(A,DKL) * T(B,DKL)
C----------------------------------------------------------------------
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,S2,BUF1)
      CALL SRT5 (NREP,MULTB,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP = 1, NREP
         M = NV(IRP)
         N = NV(IRP)
         K = NVOOT(IRP)
         CALL XGEMM ('N','C',M,N,K,A1,BUF1(OFF1),M,BUF2(OFF2),N,
     &               A1,GVV(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C
      IF (CARITH) THEN
         CALL CONJUGA (NFVV,GVV,1)
      ENDIF
C
      RETURN
      END
