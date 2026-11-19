      SUBROUTINE T3CORR (FVO,EPS,T1,T2,T2S,T2T,BUF1,BUF2,BUF3,TAU3,
     &                   VVVO,VVOO,VOOO)
C
      use interface_to_mpi
#if defined (VAR_MPI)
      use integer_kind_mpilib
#endif
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate non-iterative T3 corrections DE(+T), DE((T)), DE(-T)
C     Use formulation of Deegan and Knowles, CPL 227 (1994), 321.
C     Note that their W is called in the description because we use
C     W for the anti-symmetrized integrals.
C     This version keeps everything in core.
C
C     Parallelization: the VVVO integrals are distributed over the
C     processes (columns BC, see DISLIST), so each process computes
C     its slice of the VVVO contribution to U'(ABC) for every triple
C     IJK. The triples are dealt out round-robin: the slices of
C     triple number N are gathered on process MOD(N,NMPROC). After
C     each round of NMPROC triples every process completes its own
C     triple (VOOO contribution and energy, see T3REST). The energy
C     contributions are summed over the processes at the end.
C     During the gathers BUF3 holds the local slice and BUF2 the
C     gathered slices of the triple owned by this process.
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
      REAL*8 BUF1(*),BUF2(*),BUF3(*),TAU3(*)
      REAL*8 VVVO(*)
      REAL*8 VOOO(*),VVOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER*8 I8,IT3OFF1,JT3OFF1,KT3OFF1,EBCI,EBCJ,EBCK
      integer aeij,aeik,aejk,arep,bcrep,erep,i,ij,ijkrep,ijoff1
      integer ijrep,ik,ikoff1,irep,itrip,iroot,j,jk,jkoff1,jrep,k
      integer krep,m,nn,o,offx
#if defined (VAR_MPI)
      LOGICAL HAVEP
      integer ipend(6),ip,nloc
      integer(kind=integer_kind) :: nwords
      integer(kind=integer_kind), allocatable :: mydist(:),alldist(:),
     &                                           icol(:,:),ncol(:,:),
     &                                           cnt(:,:),dsp(:,:)
      real*8 ebuf(3),esum(3)
#endif
C---------------Executable code--------------------------------------
C
C     Initialize the timing.
C
      CALL XTIME(0,1,'-- T3CORR: Integral resorting      ')
C
      ET1 = D0
      ET2 = D0
      ET3 = D0

      CALL GETVOOO(VOOO)
      CALL GETVVOO(VVOO)
C
C     Bring arrays in correct order
C
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,T2,T2S)
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,LOVVOT,LLOVVT,T2,T2T)

C call SRT23 with appropriate arguments. We use the TAU3 array as scratch space for the integral
C sort. This array is dimensioned like T2.

      CALL SRT23 (NREP,MULTB,NVO,NV,NO,NT3VVT,NT3VVVT,KT3VVVOT,
     &            LLT3VVVT,TAU3,VVVO,NDIMT2)
      IF (CARITH) THEN
C
C      Take complex conjugate. We can not use CONJUGA as the length of this
C      array may exceed the I4 adressing space.
C
         DO I8 = 2, 2*KT3VVVOT(NREP+1), 2
            VVVO(I8) = - VVVO(I8)
         ENDDO
      ENDIF
C
C     Compute TAU3 that is needed in evaluating the -T correction
C
      CALL GETTAU3 (T1,T2,TAU3)

#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
C
C        Columns BC of the VVVO integrals held by each process: first
C        column (offset) ICOL and number of columns NCOL per irrep.
C        CNT and DSP: size and position (in REAL*8 words) of the slice
C        of U'(ABC) computed by each process, for each irrep IJK.
C
         ALLOCATE (MYDIST(2*NREP),ALLDIST(2*NREP*NMPROC))
         ALLOCATE (ICOL(0:NMPROC-1,NREP),NCOL(0:NMPROC-1,NREP))
         ALLOCATE (CNT(NMPROC,NREP),DSP(NMPROC,NREP))
         DO BCREP = 1, NREP
            MYDIST(BCREP)      = IDIST(1,2,BCREP)
            MYDIST(NREP+BCREP) = NT3VVT(BCREP)
         ENDDO
         NWORDS = 2*NREP
         CALL INTERFACE_MPI_ALLGATHER (MYDIST,NWORDS,ALLDIST,NWORDS,
     &                                 GLOBAL_COMMUNICATOR)
         DO IP = 0, NMPROC-1
            DO BCREP = 1, NREP
               ICOL(IP,BCREP) = ALLDIST(2*NREP*IP+BCREP)
               NCOL(IP,BCREP) = ALLDIST(2*NREP*IP+NREP+BCREP)
            ENDDO
         ENDDO
         DO IJKREP = 1, NREP
            DO IP = 0, NMPROC-1
               NLOC = 0
               DO BCREP = 1, NREP
                  AREP = MULTB(BCREP+NREP,IJKREP,2)
                  NLOC = NLOC + NV(AREP) * NCOL(IP,BCREP) * RCW
               ENDDO
               CNT(IP+1,IJKREP) = NLOC
               IF (IP.EQ.0) THEN
                  DSP(IP+1,IJKREP) = 0
               ELSE
                  DSP(IP+1,IJKREP) = DSP(IP,IJKREP) + CNT(IP,IJKREP)
               ENDIF
            ENDDO
         ENDDO
         DEALLOCATE (MYDIST,ALLDIST)
         HAVEP = .FALSE.
      ENDIF
#endif
C
      CALL XTIME(0,2,'-- T3CORR: Integral resorting      ')
C
      ITRIP = 0
      DO 300 IJKREP = 1, NREP
      DO 200 KREP = 1, NREP
      IJREP = MULTB(KREP,IJKREP,2)
      DO 100 JREP = KREP, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 100
      DO K = 1, NO(KREP)
      DO J = 1+MERGE(K,0,KREP.EQ.JREP), NO(JREP)
      DO I = 1+MERGE(J,0,JREP.EQ.IREP), NO(IREP)
C
      IJ = IIOOT(IREP,JREP) + (J-1) * NO(IREP) + I
      IK = IIOOT(IREP,KREP) + (K-1) * NO(IREP) + I
      JK = IIOOT(JREP,KREP) + (K-1) * NO(JREP) + J
      IF (IREP.EQ.JREP) IJ = IJ - J*(J+1)/2
      IF (IREP.EQ.KREP) IK = IK - K*(K+1)/2
      IF (JREP.EQ.KREP) JK = JK - K*(K+1)/2
      IJOFF1 = IVVOOT(IJREP) + (IJ-1) * NVV(IJREP)
      IKOFF1 = IVVOOT(MULTB(IREP,KREP,1)) + (IK-1) *
     &         NVV(MULTB(IREP,KREP,1))
      JKOFF1 = IVVOOT(MULTB(JREP,KREP,1)) + (JK-1) *
     &         NVV(MULTB(JREP,KREP,1))

      IT3OFF1  = KT3VVVOT(IREP) + (I-1) * NT3VVVT(IREP)
      JT3OFF1  = KT3VVVOT(JREP) + (J-1) * NT3VVVT(JREP)
      KT3OFF1  = KT3VVVOT(KREP) + (K-1) * NT3VVVT(KREP)
C
      CALL XTIME(0,1,'-- T3CORR: VOVV contraction        ')
C
C     The local slice of U'(ABC) is built in BUF1 (serial run, where it
C     is the complete U') or in BUF3 (parallel run). For each irrep BC
C     the slice holds the NT3VVT(BCREP) columns BC of this process.
C
      IF (NMPROC.EQ.1) THEN
         CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF1,1)
         CALL ADDVVVO (BUF1)
#if defined (VAR_MPI)
      ELSE
         CALL XCOPY (CNT(MYPROC+1,IJKREP)/RCW,A0,0,BUF3,1)
         CALL ADDVVVO (BUF3)
#endif
      ENDIF
C
      CALL XTIME(0,2,'-- T3CORR: VOVV contraction        ')
C
      IF (NMPROC.EQ.1) THEN
C
C        Serial run: U' is complete, finish this triple right away
C
         CALL T3REST (IJKREP,KREP,JREP,K,J,I,FVO,EPS,T1,T2,T2T,
     &                BUF1,BUF2,BUF3,TAU3,VVOO,VOOO)
#if defined (VAR_MPI)
      ELSE
C
C        Gather the slices of this triple on its owner IROOT
C
         CALL XTIME(0,1,'-- T3CORR: addition of partial U   ')
         IROOT = MOD(ITRIP,NMPROC)
         NWORDS = CNT(MYPROC+1,IJKREP)
         CALL INTERFACE_MPI_GATHERV (BUF3,NWORDS,BUF2,
     &        CNT(:,IJKREP),DSP(:,IJKREP),IROOT,GLOBAL_COMMUNICATOR)
         CALL XTIME(0,2,'-- T3CORR: addition of partial U   ')
         IF (IROOT.EQ.MYPROC) THEN
            IPEND(1) = IJKREP
            IPEND(2) = KREP
            IPEND(3) = JREP
            IPEND(4) = K
            IPEND(5) = J
            IPEND(6) = I
            HAVEP = .TRUE.
         ENDIF
C
C        End of a round: every process completes the triple it owns
C
         IF (IROOT.EQ.NMPROC-1) THEN
            CALL T3UNPK (IPEND(1),ICOL,NCOL,BUF2,BUF1)
            CALL T3REST (IPEND(1),IPEND(2),IPEND(3),IPEND(4),
     &                   IPEND(5),IPEND(6),FVO,EPS,T1,T2,T2T,
     &                   BUF1,BUF2,BUF3,TAU3,VVOO,VOOO)
            HAVEP = .FALSE.
         ENDIF
#endif
      ENDIF
      ITRIP = ITRIP + 1
      ENDDO
      ENDDO
      ENDDO
  100 CONTINUE
  200 CONTINUE
  300 CONTINUE
C
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
C
C        Complete the triple of the last, incomplete, round
C
         IF (HAVEP) THEN
            CALL T3UNPK (IPEND(1),ICOL,NCOL,BUF2,BUF1)
            CALL T3REST (IPEND(1),IPEND(2),IPEND(3),IPEND(4),
     &                   IPEND(5),IPEND(6),FVO,EPS,T1,T2,T2T,
     &                   BUF1,BUF2,BUF3,TAU3,VVOO,VOOO)
         ENDIF
C
C        Sum the energy contributions of all processes
C
         EBUF(1) = ET1
         EBUF(2) = ET2
         EBUF(3) = ET3
         NWORDS = 3
         CALL INTERFACE_MPI_ALLREDUCE_R1_WORK_F77 (EBUF,ESUM,NWORDS,
     &        OP_MPI_SUM,GLOBAL_COMMUNICATOR)
         ET1 = ESUM(1)
         ET2 = ESUM(2)
         ET3 = ESUM(3)
         DEALLOCATE (ICOL,NCOL,CNT,DSP)
      ENDIF
#endif
C
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
 1071 FORMAT(//,'  Perturbative treatment of triple excitations',/)
 1072 FORMAT(' SCF energy :',T40,F25.15)
 1073 FORMAT(' CCSD correlation energy :',T40,F25.15)
 1074 FORMAT(' 4th order triples correction :',T40,F25.15)
 1075 FORMAT(' 5th order triples (T) correction :',T40,F25.15)
 1076 FORMAT(' 5th order triples -T  correction :',T40,F25.15)
 1077 FORMAT(' Total CCSD+T  energy :',T40,F25.15)
 1078 FORMAT(' Total CCSD(T) energy :',T40,F25.15)
 1079 FORMAT(' Total CCSD-T  energy :',T40,F25.15)
C
      RETURN
C
      CONTAINS
C
      SUBROUTINE ADDVVVO (XB)
C
C     Add the VVVO contributions of the local columns BC to U'(ABC):
C       U'(ABC) =  T(A,E:IJ) * W(E,BC:K) + T(A,E:JK) * W(E,BC:I)
C                - T(A,E:IK) * W(E,BC:J)
C     XB holds, for each irrep BC, the NT3VVT(BCREP) local columns.
C
      REAL*8 XB(*)
      OFFX = 1
      DO BCREP = 1, NREP
         AREP = MULTB(BCREP+NREP,IJKREP,2)
         M = NV(AREP)
         NN = NT3VVT(BCREP)
C        U'(ABC) = T(A,E:IJ) * W(E,BC:K)
         EREP = MULTB(KREP,BCREP+NREP,2)
         O = NV(EREP)
         AEIJ = (IJOFF1 + IIVV(AREP,EREP)) * RCW + 1
         EBCK = (KT3OFF1 + LLT3VVVT(EREP,BCREP)) * RCW + 1
         CALL XGEMM ('N','N',M,NN,O,A1,T2S(AEIJ),M,
     &               VVVO(EBCK),O,A1,XB(OFFX),M)
C        U'(ABC) = U'(ABC) + T(A,E:JK) * W(E,BC:I)
         EREP = MULTB(IREP,BCREP+NREP,2)
         O = NV(EREP)
         AEJK = (JKOFF1 + IIVV(AREP,EREP)) * RCW + 1
         EBCI = (IT3OFF1 + LLT3VVVT(EREP,BCREP)) * RCW + 1
         CALL XGEMM ('N','N',M,NN,O,A1,T2S(AEJK),M,
     &               VVVO(EBCI),O,A1,XB(OFFX),M)
C        U'(ABC) = U'(ABC) - T(A,E:IK) * W(E,BC:J)
         EREP = MULTB(JREP,BCREP+NREP,2)
         O = NV(EREP)
         AEIK = (IKOFF1 + IIVV(AREP,EREP)) * RCW + 1
         EBCJ = (JT3OFF1 + LLT3VVVT(EREP,BCREP)) * RCW + 1
         CALL XGEMM ('N','N',M,NN,O,-A1,T2S(AEIK),M,
     &               VVVO(EBCJ),O,A1,XB(OFFX),M)
         OFFX = OFFX + M * NN * RCW
      ENDDO
      END SUBROUTINE ADDVVVO
C
      END
#if defined (VAR_MPI)
      SUBROUTINE T3UNPK (IJKREP,ICOL,NCOL,XG,U)
C
      use integer_kind_mpilib
      implicit none
C
C---------------Description--------------------------------------------
C
C     Copy the slices of U'(ABC) gathered from all processes (stored
C     one process after the other in XG) to their place in U.
C
C---------------Calling variables--------------------------------------
C
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"
      INTEGER IJKREP
      integer(kind=integer_kind) :: ICOL(0:NMPROC-1,NREP),
     &                              NCOL(0:NMPROC-1,NREP)
      REAL*8 XG(*),U(*)
C
C---------------Local variables--------------------------------------
C
      integer arep,bcrep,ip,ipos,iu,nw
C
C---------------Executable code--------------------------------------
C
      IPOS = 1
      DO IP = 0, NMPROC-1
         DO BCREP = 1, NREP
            AREP = MULTB(BCREP+NREP,IJKREP,2)
            NW = NV(AREP) * NCOL(IP,BCREP) * RCW
            IU = (IVVVT(BCREP,IJKREP) + NV(AREP)*ICOL(IP,BCREP))*RCW
            U(IU+1:IU+NW) = XG(IPOS:IPOS+NW-1)
            IPOS = IPOS + NW
         ENDDO
      ENDDO
C
      RETURN
      END
#endif
      SUBROUTINE T3REST (IJKREP,KREP,JREP,K,J,I,FVO,EPS,T1,T2,T2T,
     &                   BUF1,BUF2,BUF3,TAU3,VVOO,VOOO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Complete the triples correction for one triple IJK: add the VOOO
C     contribution to U'(ABC) (held in BUF1 on entry) and add the
C     contributions of this triple to DE(+T), DE((T)) and DE(-T).
C     Called by T3CORR. BUF2 and BUF3 are used as scratch.
C
C---------------Calling variables--------------------------------------
C
      INTEGER IJKREP,KREP,JREP,K,J,I
      REAL*8 FVO(*),EPS(*)
      REAL*8 T1(*),T2(*),T2T(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*),TAU3(*)
      REAL*8 VOOO(*),VVOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "results.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 EPSIJK,EIJKABC,VMIN(2)
      integer a,aa,ab,abc,abc1,abc2,abij,abik,abjk,abrep,aljk
      integer ac,acb,acrep,alij,alik,amin,arep,b,bac
      integer bb,bc,bca,bcrep,bmin,brep,c,cab,cc,ci,cj,ck,crep
      integer ii,ij,ijoff2,ijoff3,ijrep,ik,ikoff2,ikoff3,ikrep
      integer ioff2,irep,jj,jk,jkoff2,jkoff3,jkrep,joff2,kk,koff2
      integer lbci,lbcj,lbck,lrep,m,n,o,off1
C
C---------------Executable code--------------------------------------
C
      CALL XTIME(0,1,'-- T3CORR: energy calculation      ')
C
      IJREP = MULTB(KREP,IJKREP,2)
      IREP = MULTB(JREP,IJREP+NREP,2)
      IKREP = MULTB(IREP,KREP,1)
      JKREP = MULTB(JREP,KREP,1)
      II = IO(IREP) + I
      JJ = IO(JREP) + J
      KK = IO(KREP) + K
      EPSIJK = - EPS(II) - EPS(JJ) - EPS(KK)
C
      IJ = IIOOT(IREP,JREP) + (J-1) * NO(IREP) + I
      IK = IIOOT(IREP,KREP) + (K-1) * NO(IREP) + I
      JK = IIOOT(JREP,KREP) + (K-1) * NO(JREP) + J
      IF (IREP.EQ.JREP) IJ = IJ - J*(J+1)/2
      IF (IREP.EQ.KREP) IK = IK - K*(K+1)/2
      IF (JREP.EQ.KREP) JK = JK - K*(K+1)/2
      IJOFF2 = IVOOOT(IJREP) + (IJ-1) * NVO(IJREP)
      IKOFF2 = IVOOOT(IKREP) + (IK-1) * NVO(IKREP)
      JKOFF2 = IVOOOT(JKREP) + (JK-1) * NVO(JKREP)
      IJOFF3 = IVVOOTT(IJREP) + (IJ-1) * NVVT(IJREP)
      IKOFF3 = IVVOOTT(IKREP) + (IK-1) * NVVT(IKREP)
      JKOFF3 = IVVOOTT(JKREP) + (JK-1) * NVVT(JKREP)
      IOFF2  = LOVVOT(IREP) + (I-1) * NOVVT(IREP)
      JOFF2  = LOVVOT(JREP) + (J-1) * NOVVT(JREP)
      KOFF2  = LOVVOT(KREP) + (K-1) * NOVVT(KREP)
C
      OFF1 = 1
      DO BCREP = 1, NREP
      AREP = MULTB(BCREP+NREP,IJKREP,2)
      M = NV(AREP)
      N = NVVT(BCREP)
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
      OFF1 = OFF1 + M * N * RCW
      ENDDO
C
C-------------------------------------
C U(ABC) = U'(ABC) - U'(BAC) + U'(CAB)
C-------------------------------------
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
C-------------------------------------------------
C Y2'(ABC) = V(AB:IJ) * T(C:K) + T(AB:IJ) * F(C:K)
C          + V(AB:JK) * T(C:I) + T(AB:JK) * F(C:I)
C          - V(AB:IK) * T(C:J) - T(AB:IK) * F(C:J)
C-------------------------------------------------
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
C-----------------------------------------
C Y2(ABC) = Y2'(ABC) + Y2'(BCA) - Y2'(ACB)
C-----------------------------------------
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
C-----------------------------------------
C DE(+T)   = - U(ABC) * U(ABC) / EPS(ABC)
C DE'((T)) = - U(ABC) * Y2(ABC) / EPS(ABC)
C-----------------------------------------
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
C--------------------------------
C Y3'(ABC) = TAU3(AB:IJ) * T(C:K)
C          + TAU3(AB:JK) * T(C:I)
C          - TAU3(AB:IK) * T(C:J)
C--------------------------------
      CALL XCOPY (IVVVT(NREP+1,IJKREP),A0,0,BUF2,1)
      ABIJ = IJOFF3 * RCW + 1
      ABJK = JKOFF3 * RCW + 1
      ABIK = IKOFF3 * RCW + 1
      DO C = 1, NV(KREP)
         CK = (IVO(KREP) + (K-1) * NV(KREP) + C - 1) * RCW + 1
         ABC = (JVVVT(KREP,IJKREP)+(C-1)*NVVT(IJREP)) * RCW + 1
         CALL XAXPY (NVVT(IJREP),T1 (CK),TAU3(ABIJ),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(IREP)
         CI = (IVO(IREP) + (I-1) * NV(IREP) + C - 1) * RCW + 1
         ABC = (JVVVT(IREP,IJKREP)+(C-1)*NVVT(JKREP)) * RCW + 1
         CALL XAXPY (NVVT(JKREP),T1 (CI),TAU3(ABJK),1,BUF2(ABC),1)
      ENDDO
      DO C = 1, NV(JREP)
         CJ = (IVO(JREP) + (J-1) * NV(JREP) + C - 1) * RCW + 1
         ABC = (JVVVT(JREP,IJKREP)+(C-1)*NVVT(IKREP)) * RCW + 1
         VMIN(1) = - T1(CJ)
         IF (CARITH) VMIN(2) = - T1(CJ+1)
         CALL XAXPY (NVVT(IKREP),VMIN,TAU3(ABIK),1,BUF2(ABC),1)
      ENDDO
C-----------------------------------------
C Y3(ABC) = Y3'(ABC) + Y3'(BCA) - Y3'(ACB)
C-----------------------------------------
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
C-----------------------------
C DE'(-T) = + U(ABC) * Y3(ABC)
C-----------------------------
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
      CALL XTIME(0,2,'-- T3CORR: energy calculation      ')
C
      RETURN
      END
