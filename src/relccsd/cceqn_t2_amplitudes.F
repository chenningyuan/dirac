      SUBROUTINE amplitude_equation_T2 (T1,T2,
     &          GO,GV,HOV,H,AZ,BUF1,BUF2,BUF3,NBUF3,S2)
C
CPi use contraction module to include tensor contraction routiones
      use contraction
C
      implicit none
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
CPi      INTEGER NBUF3
      INTEGER NBUF2,NBUF3
CPi      REAL*8 T1(*),T2(*)
      REAL*8 T1(NDIMT1*RCW), T2(NDIMT2*RCW)
CPi      REAL*8 S2(*)
      REAL*8 S2(NDIMT2*RCW)
CPi      REAL*8 GO(*),GV(*),HOV(*)
CPi NFOO=MOO(1). See also ccdriv.F:2426
      REAL*8 GO(NFOO*RCW),GV(NFVV*RCW),HOV(NFVO*RCW)
      REAL*8 AZ(*),H(*)
CPi      REAL*8 BUF1(*),BUF2(*),BUF3(*)
      REAL*8 BUF1(NBUF2*RCW),BUF2(NBUF2*RCW),BUF3(NBUF3*RCW)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "complex.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "symm.inc"
#include "ccpar.inc"
C
      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ
      REAL*8 T1VAL(2)
      DATA T1VAL /2*0.0D0/
      complex*16 zdotc
      real*8 ddot
      integer abji,c,cj,crep,i,imin,irep,irp,istart,j,ij,jmin,jrep,k,kld
      integer m,mint,n,nint,off1,off2,off3,t1ci,abij,ji
C
C---------------Executable code--------------------------------------
C
C---- Initialize to zero...
      CALL XCOPY (NDIMT2,A0,0,S2,1)
C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * AZ(KL,IJ)
C---------------------------------------------
      CALL XTIME(0,1,'--- T2EQNS: TAU*AINTM contraction')
      CALL GETTAU (T1,T2,BUF1)
          
      OFF3 = 1
      DO IREP = 1, NREP
        M    = NVVT(IREP)
        N    = NOOT(IREP)
        KLD  = NOOT(IREP)
        K    = IDIST(1,5,IREP)
        OFF1 = IDIST(2,5,IREP)
        OFF2 = IDIST(3,5,IREP)
        CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,AZ(OFF2),KLD,
     &              A1,S2(OFF3),M)
        OFF3 = OFF3 + M * N * RCW
      ENDDO

C ***
C *** original CNTRCT call in serial version
C ***
C *** CALL CNTRCT ('N','N',NVVT,NOOT,NOOT,A1,BUF1,AZ,A1,S2,NREP)
C ***

      CALL XTIME(0,2,'--- T2EQNS: TAU*AINTM contraction')

C ** we will not parallelize a summation over just one occupied index
C ** if the system is small this can lead to zero batches. Only master
C ** does the T2(AB,IK) * GO(K,J) contraction

C ********************
C ** TASK SPLITTING **
C ********************

      IF(MYPROC.EQ.MASTER) THEN
C----------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PIJ [ T2(AB,IK) * GO(K,J) ]
CPi Either the minus sign is wrong or the IK should be KI.
C----------------------------------------------------------------
CPi-old      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
CPi-old     &             T2,BUF1)
CPi-old      CALL CNTRCT ('N','N',NVVOT,NO,NO,A1,BUF1,GO,A0,BUF2,NREP)
CPi-old      CALL SRT1T3 (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
CPi-old     &             BUF2,BUF1)
CPi-old      CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)
CPi-begin
      call contraction_244 ((/"o3","o1"/),(/"p1","p2","o3","o2"/),
     &  (/"p1","p2","o1","o2"/),S2,-1.0d0,1.0d0,nrep,LeftTensor=GO,
     &   RightTensor=T2)
CPi-end

C ***************************
C ** END OF TASK SPLITTING **
C ***************************

      ENDIF
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PIJ [ W(AB,CJ) * T1(C,I) ]
CPi Again problem with sign. C and J are switched compared to reference,
CPi but sign is the same...
C-----------------------------------------------------------------------
CPi Modify this label maybe:
      CALL XTIME(0,1,'--- T2EQNS: VOVV*T1              ')
CPi-begin
      CALL GETVVOV_incore (BUF3)
      CALL DELINT ('VVOV','KKKD',BUF3,0,0,0)
      call contraction_424 ((/"p1","p2","o1","p3"/),(/"p3","o2"/),
     &    (/"p1","p2","o1","o2"/),S2,1.0d0,1.0d0,nrep,LeftTensor=BUF3,
     &    RightTensor=T1)
CPi-end
CPi-old      DO 30 IRP = 1, NREP
CPi-old         IF (NVO(IRP).EQ.0) GOTO 30
CPi-old         DONE = .FALSE.
CPi-old         ISTART = 0
CPi-old         MINT = NBUF3/NVO(IRP)  ! max. # of NVO batches that fit in mem.
CPi-old  10     CONTINUE
CPi-old         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
CPi-old         CALL DELINT ('VOVV','DKKK',BUF3,IRP,ISTART,NINT)
CPi-old         IF (CARITH) CALL CONJUGA (NVO(IRP)*NINT,BUF3,1)
CPi-old         IJ = 0
CPi-old         DO 20 JREP = 1, NREP
CPi-old           CREP = MULTB(JREP,IRP+NREP,2)
CPi-old           IREP = CREP
CPi-old           IF (IREP.LT.JREP) GOTO 20
CPi-old           DO J = 1, NO(JREP)
CPi-old              IMIN = 1
CPi-old              IF (JREP.EQ.IREP) IMIN = J + 1
CPi-old              DO I = IMIN, NO(IREP)
CPi-old                 IJ = IJ + 1
CPi-old                 ABIJ = (IVVOOTT(IRP)+(IJ-1)*NVVT(IRP)+ISTART)*RCW+1
CPi-old                 DO C = 1, NV(CREP)
CPi-old                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
CPi-old                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
CPi-old                    T1VAL(1) = T1(T1CI)
CPi-old                    IF (CARITH) T1VAL(2) = T1(T1CI+1)
CPi-old                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABIJ),1)
CPi-old                 ENDDO
CPi-old              ENDDO
CPi-old           ENDDO
CPi-old 20      CONTINUE
CPi-old         JI = 0
CPi-old         DO 21 IREP = 1, NREP
CPi-old           JREP = MULTB(IREP,IRP+NREP,2)
CPi-old           CREP = IREP
CPi-old           IF (JREP.LT.IREP) GOTO 21
CPi-old           DO I = 1, NO(IREP)
CPi-old              JMIN = 1
CPi-old              IF (JREP.EQ.IREP) JMIN = I + 1
CPi-old              DO J = JMIN, NO(JREP)
CPi-old                 JI = JI + 1
CPi-old                 ABJI = (IVVOOTT(IRP)+(JI-1)*NVVT(IRP)+ISTART)*RCW+1
CPi-old                 DO C = 1, NV(CREP)
CPi-old                    CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
CPi-old                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
CPi-old                    T1VAL(1) = - T1(T1CI)
CPi-old                    IF (CARITH) T1VAL(2) = - T1(T1CI+1)
CPi-old                    CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABJI),1)
CPi-old                 ENDDO
CPi-old              ENDDO
CPi-old           ENDDO
CPi-old 21      CONTINUE
CPi-old         IF (.NOT.DONE) THEN
CPi-old            ISTART = ISTART + NINT
CPi-old            GOTO 10
CPi-old         ENDIF
CPi-old 30    CONTINUE
      CALL XTIME(0,2,'--- T2EQNS: VOVV*T1              ')

C ********************
C ** TASK SPLITTING **
C ********************
C ** here a VOOV dependent contribution (Term 7 in the paper) is to be
C ** calculated only by the master

      IF(MYPROC.EQ.MASTER) THEN

C--------------------------------
C Terms specific for T2 equations
C--------------------------------
                   
C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [T1(A,K)*W(BK,IC)*T1(C,J)]
C-----------------------------------------------------------------------
         CALL GETVOVO (BUF1)
         CALL DELINT ('VOVO','KDDK',BUF1,0,0,0)
         CALL SRT19 (NREP,MULTB,LFA,NVO,NV,NO,NVOO,KVOOV,KKVOO,
     &               BUF1,BUF2)
         CALL CNTRCT ('N','N',NVOO,NO,NV,A1,BUF2,T1,A0,BUF1,NREP)
         CALL SRT1T3 (NREP,MULTB,LTR,NVO,NO,NO,NVOO,KVOOO,KKVOO,
     &                BUF1,BUF3)
C-------------------------------------------------------------------
C        Integrals are now in VO,OOT format : combine with next term
CPi Which next term?
C-------------------------------------------------------------------
           

      ENDIF

C ****************************
C ** END OF TASK SPLITTING. **
C ****************************
C ** this contribution coming from the intermediate G(A,C) has to be
C ** calculated on each node. The fact that G(A,C) contains LOCAL +
C ** DISTRIBUTED parts as well is accounted for by the original calculation
C ** of the G(A,C) array before entering T2EQNS

C-------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB [ GV(A,C) * T2(CB,IJ)]
C-------------------------------------------------
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GV,BUF2,A0,BUF1,NREP)

C ** the contraction of G(A,C) with the T2 amplitudes now is in BUF1 !
C ** these are the partial contributions available on EACH node


C ***
C ****** task splitting
C ***
      IF(MYPROC.EQ.MASTER) THEN

C------------------------------------------------
C S(AB,IJ) = S(AB,IJ) - PAB [ T1(A,K) *W(KB,IJ) ]
C------------------------------------------------
                   
         CALL GETVOOO (BUF2)
         CALL DELINT ('VOOO','KDKK',BUF2,0,0,0)
         CALL XAXPY (NV2,-A1,BUF3,1,BUF2,1)
          
           
      CALL SRT1L1 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF3)
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NO,NV,NOOT,LOVOOT,LLVOOT,
     &             BUF3,BUF2)
      CALL CNTRCT ('N','N',NV,NVOOT,NO,-A1,T1,BUF2,A1,BUF1,NREP)

C ***
C *** end of task splitting
C ***
      ENDIF

      CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF1,BUF2)
      CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)

C ** now the final H*T2 contraction. H is complete since we communicated
C ** it already.

C-----------------------------------------------------------------------
C S(AB,IJ) = S(AB,IJ) + PAB.PIJ [H(AI,CK)*T2(CK,BJ)]
C-----------------------------------------------------------------------

      CALL XTIME(0,1,'--- T2EQNS: HINTM*T2           ')
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              T2,BUF1)
      CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF2,1)

      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M    = MVO(IREP)
         N    = MVO(JREP)
         KLD  = MVO(IREP)
         K    = IDIST(1,10,JREP)
         OFF1 = IDIST(2,10,JREP)
         OFF2 = IDIST(3,10,JREP)
         CALL XGEMM ('N','N',M,N,K,A1,H(OFF1),M,BUF1(OFF2),KLD,
     &                A0,BUF2(OFF3),M)
         OFF3 = OFF3 + M * N * RCW
      ENDDO

      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
     &              BUF2,BUF1)
      CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)

      CALL XTIME(0,2,'--- T2EQNS: HINTM*T2           ')
C
C
C  Finally we ADD the B INTERMEDIATE CONTRIBUTION
C
C
      CALL XTIME(5,1,'-- BINTM                       ')
      CALL BINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,S2)
      CALL XTIME(5,2,'-- BINTM                       ')

      END
