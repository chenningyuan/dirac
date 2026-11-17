      SUBROUTINE T1EQN1 (TAU,S1)
C
      implicit none
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
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"

      real*8,intent(inout) :: S1(ndimt1*rcw)
      real*8,intent(in)    :: TAU(ndimt2*rcw)

C---------------Local variables--------------------------------------
C
      real*8, allocatable ::  tau_sorted(:), vooo(:), vooo_sorted(:)
      
C---------------Executable code--------------------------------------
C
      allocate (tau_sorted(LVVOOT(NREP+1)*RCW))
      allocate (vooo(IVOOOT(NREP+1)*RCW))
      allocate (vooo_sorted(LOVOOT(NREP+1)*RCW))

C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------

C ** TAU now reorganised in BUF1: NV*NVOOT

      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             TAU,tau_sorted)
      CALL GETVOOO (vooo)
      CALL DELINT ('VOOO','DKDD',vooo,0,0,0)

C     --------------------------------
C     Sort from 12,34 to 134,2 (VOO,O)
C     --------------------------------
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo,vooo_sorted)

C ** VOOO integrals now reorganised in BUF2: NVOOT * NO

C     ------------------------
C     V(K,L,C,I) = V*(C,I,K,L)
C     ------------------------
      IF (CARITH) CALL CONJUGA (LOVOOT(NREP+1),vooo_sorted,1)

C
      CALL CNTRCT 
     &   ('N','N',NV,NO,NVOOT,A1,tau_sorted,vooo_sorted,A1,S1,NREP)
C

      deallocate (tau_sorted)
      deallocate (vooo)
      deallocate (vooo_sorted)

      RETURN
      END
      SUBROUTINE T1EQN2 (TAU,BUF3,NBUF3,S1)
C
      implicit none
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
      REAL*8 S1(*),TAU(*)
      REAL*8 BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE
      integer i,ik,ikb,irep,irp,istart,k,k0,kb,kbmax,kk,krep,nint
      integer cdik,vak,vai,t1ai,t1ak,ncd,noi,nok,nvi,nvk,lvo,lvvt
C     The TAU gather buffer is limited to the size NBUF3 of the
C     integral buffer budgeted by the caller, with MXGATH (REAL*8
C     words) as an absolute upper limit
      integer, parameter :: MXGATH = 4194304
      real*8, allocatable :: tbuf(:)
C
C---------------Executable code--------------------------------------
C
C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)
C--------------------------------------------
C
C     V(AK,CD) is held in BUF3 with leading dimension NVO(IRP), for a
C     batch of NINT CD pairs starting at ISTART+1. TAU(CD,IK) has
C     leading dimension NVVT(IRP) with IK running over I>K pairs.
C     Both terms are evaluated as matrix-matrix products:
C       S(A,I) += V(A,K;CD) * TAU(CD;I,K)   (one GEMM per K)
C       S(A,K) -= V(A,I;CD) * TAU(CD;I,K)   (one GEMM per I)
C     For the second term TAU columns for fixed I are equidistant when
C     irep(K) < irep(I); in the diagonal block (K<I within one irrep)
C     they are gathered into a contiguous buffer first (at most
C     MIN(NBUF3,MXGATH) words, see above).
C
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         NINT = NBUF3/NVO(IRP)  ! maximum # of NVO batches
         LVO  = NVO(IRP)
         LVVT = NVVT(IRP)

C ** NINT is updated on exit of GETVOVV !

  10     CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,NINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)
         IF (NINT.LE.0) GOTO 25
         NCD = NINT*RCW
         IKB = 0
         DO 20 KREP = 1, NREP
           IREP = MULTB(KREP,IRP+NREP,2)
           IF (KREP.GT.IREP) GOTO 20
           NOI = NO(IREP)
           NOK = NO(KREP)
           NVI = NV(IREP)
           NVK = NV(KREP)
           IF (KREP.LT.IREP) THEN
C
C             Off-diagonal block: IK = IKB + (K-1)*NOI + I
C
              DO K = 1, NOK
                 IK   = IKB + (K-1)*NOI + 1
                 CDIK = (IVVOOTT(IRP)+(IK-1)*LVVT+ISTART)*RCW+1
                 VAK  = (IIVO(IREP,KREP)+(K-1)*NVI)*RCW+1
                 T1AI = IVO(IREP)*RCW+1
                 CALL XGEMM ('N','N',NVI,NOI,NINT,A1,BUF3(VAK),LVO,
     &                       TAU(CDIK),LVVT,A1,S1(T1AI),NVI)
              ENDDO
              DO I = 1, NOI
                 IK   = IKB + I
                 CDIK = (IVVOOTT(IRP)+(IK-1)*LVVT+ISTART)*RCW+1
                 VAI  = (IIVO(KREP,IREP)+(I-1)*NVK)*RCW+1
                 T1AK = IVO(KREP)*RCW+1
                 CALL XGEMM ('N','N',NVK,NOK,NINT,-A1,BUF3(VAI),LVO,
     &                       TAU(CDIK),LVVT*NOI,A1,S1(T1AK),NVK)
              ENDDO
              IKB = IKB + NOK*NOI
           ELSE
C
C             Diagonal block: IK = IKB + (K-1)*NOI - K*(K-1)/2 + I - K
C
              DO K = 1, NOK-1
                 IK   = IKB + (K-1)*NOI - K*(K-1)/2 + 1
                 CDIK = (IVVOOTT(IRP)+(IK-1)*LVVT+ISTART)*RCW+1
                 VAK  = (IIVO(IREP,KREP)+(K-1)*NVI)*RCW+1
                 T1AI = (IVO(IREP)+K*NVI)*RCW+1
                 CALL XGEMM ('N','N',NVI,NOI-K,NINT,A1,BUF3(VAK),LVO,
     &                       TAU(CDIK),LVVT,A1,S1(T1AI),NVI)
              ENDDO
              IF (NOI.GT.1.AND.NVK.GT.0) THEN
                 KBMAX = MAX(1,MIN(NOI-1,MIN(NBUF3,MXGATH)/NCD))
                 ALLOCATE (TBUF(NCD*KBMAX))
                 DO I = 2, NOI
                    VAI = (IIVO(KREP,IREP)+(I-1)*NVK)*RCW+1
                    DO K0 = 1, I-1, KBMAX
                       KB = MIN(KBMAX,I-K0)
                       DO KK = 1, KB
                          K    = K0 + KK - 1
                          IK   = IKB + (K-1)*NOI - K*(K-1)/2 + I - K
                          CDIK = (IVVOOTT(IRP)+(IK-1)*LVVT+ISTART)*RCW
                          TBUF((KK-1)*NCD+1:KK*NCD) =
     &                         TAU(CDIK+1:CDIK+NCD)
                       ENDDO
                       T1AK = (IVO(KREP)+(K0-1)*NVK)*RCW+1
                       CALL XGEMM ('N','N',NVK,KB,NINT,-A1,BUF3(VAI),
     &                             LVO,TBUF,NINT,A1,S1(T1AK),NVK)
                    ENDDO
                 ENDDO
                 DEALLOCATE (TBUF)
              ENDIF
              IKB = IKB + NOI*(NOI-1)/2
           ENDIF
 20      CONTINUE
 25      IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
 30   CONTINUE
C
      RETURN
      END
