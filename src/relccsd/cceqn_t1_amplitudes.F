      SUBROUTINE amplitude_equation_T1 (FVO,T1,T2,
     &        HO,HV,HOV,BUF1,BUF2,BUF3,NBUF2,NBUF3,S1)
C
CPi use contraction module to include tensor contraction routiones
      use contraction
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Evaluate T1 or Lambda1 equations
C     Here we have to take care of the fact that certain integral classes
C     are now split over individual nodes. The corresponding contributions
C     which are linear in the integrals are added at the end of the T1
C     routines. Therefore we have to take care of calculating the
C     non-distributed contributions to the T1 amplitudes only on the
C     master node in order to avoid double-counting in the PARTS routine !
C     all direct VOVV/VVVV contributions and VOVV/VVVV contributions stemming
C     from intermediates have to be calculated on EACH node.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      external t1eqn1
      external t1eqn2

      INTEGER NBUF2,NBUF3
      REAL*8 FVO(*)
CPi Specify dimensions for general contractions routines to work
CPi      REAL*8 T1(*),T2(*)
      REAL*8 T1(NDIMT1*RCW), T2(NDIMT2*RCW)
CPi      REAL*8 S1(*)
      REAL*8 S1(NDIMT1*RCW)
CPi      REAL*8 HO(*),HV(*),HOV(*)
      REAL*8 HV(NFVV*RCW), HO(NFVV*RCW), HOV(NFVO*RCW)
      REAL*8 BUF1(NBUF2*RCW),BUF2(NBUF2*RCW),BUF3(NBUF3*RCW)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Calling variables--------------------------------------
CPi Allocate space of ooov integrals
CPi      real*8, allocatable :: ooov(IOOOVT(NREP+1)*RCW)
      real*8, allocatable :: vooo(:)

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ,USEDZ,RIGHT
      REAL*8 SDOT,DDOT
      INTEGER irp,k,k1,k2,m,n,off1,off2
C
C---------------Executable code--------------------------------------
C
C-------------------
C Initialize to zero
C-------------------
      CALL XCOPY (NDIMT1,A0,0,S1,1)
C--------------------------------
C Terms specific for T1 equations
C--------------------------------
                   
C ********************
C ** TASK SPLITTING **
C ********************
        IF (MYPROC.EQ.MASTER) THEN
C------------------------------------------------------------
C S(A,I) = S(A,I) + T(A,K) * [HOV(K,C) - 2 * F(K,C)] * T(C,I)
C------------------------------------------------------------
         CALL XTIME(0,1,'--- T1EQNS T*[HOV - F]*T      ')
         CALL XCOPY (NFVO,HOV,1,BUF1,1)
CPi In order to add with -2*FVO we should complex conjugate HOV which is
CPi actually HVO*
         IF (CARITH) CALL CONJUGA (NFVO,BUF1,1)
         CALL XAXPY (NFVO,-A2,FVO,1,BUF1,1)
CPi This is stupid so the transpose of [HOV(K,C) - 2 * F(K,C)] should be
CPi made otherwise
         IF (CARITH) CALL CONJUGA (NFVO,BUF1,1)
CPi DELFCK only applies to FSCC and we should be able to leave it as it is
         CALL DELFCK ('VO','DD',BUF1)
CPi Use new contraction routine for [HOV(K,C) - 2 * F(K,C)] * T(C,I)
         call contraction_222 ((/"p2","o1"/),(/"p2","o2"/),(/"o1","o2"/)&
     &    ,BUF1,T1,BUF2,1.0d0,0.0d0,NREP)
CPi Use new contraction routine for T(A,K) * BUF2 and add to S1
         call contraction_222 ((/"p1","o2"/),(/"o2","o1"/),(/"p1","o1"/)&
     &    ,T1,BUF2,S1,1.0d0,1.0d0,NREP)
CPi Keep the old code just for now. PiRM
CPi         OFF1 = 1
CPi         DO IRP = 1, NREP
CPi            M = NV(IRP)
CPi            N = NO(IRP)
CPi            K1 = NO(IRP)
CPi            K2 = NV(IRP)
CPiCPi         BUF1(OFF1) * T1(OFF1) + 0*BUF2 (empty)
CPi            CALL XGEMM ('C','N',K1,N,K2,A1,BUF1(OFF1),K2,T1(OFF1),K2,
CPi     &                   A0,BUF2,K1)
CPiCPi         T1(OFF1) * BUF2 + S1(OFF1)
CPi            CALL XGEMM ('N','N',M,N,K1,A1,T1(OFF1),M,BUF2,K1,
CPi     &                   A1,S1(OFF1),M)
CPi            OFF1 = OFF1 + M * N * RCW
CPi         ENDDO
         CALL XTIME(0,2,'--- T1EQNS T*[HOV - F]*T      ')
        ENDIF
C ****************************************************
C ** END OF TASK SPLITTING. Parallel execution      **
C ****************************************************
C-----------------------------------------------------------------
C S(A,I) = S(A,I) + HOV(CK) * T2(A,C,I,K)
C-----------------------------------------------------------------
         CALL XTIME(0,1,'--- T1EQNS HOV*T2(A,C,I,K       ')
CPi use new contraction routine - SEQUENTIAL!
CPi I think what we have is HOV(CK)* and that we might have to conjugate
         call contraction_242 ((/"p2","o2"/),(/"p1","p2","o1","o2"/),   &
     &   (/"o1","p1"/),HOV,S1,1.0d0,1.0d0,nrep,RightTensor=T2)
CPi         CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
CPi     &                 MVO,JVOVO,JJVO,JJVO,T2,BUF1)
CPi
CPi         K    = IDIST(1,3,1)
CPi         OFF1 = IDIST(2,3,1)
CPi         OFF2 = IDIST(3,3,1)
CPi         M    = MVO(1)
CPi         CALL XGEMV ('N',M,K,A1,BUF1(OFF1),M,HOV(OFF2),1,A1,S1,1)
         CALL XTIME(0,2,'--- T1EQNS HOV*T2(A,C,I,K       ')
           

C ********************
C ** TASK SPLITTING **
C ********************

      IF (MYPROC.EQ.MASTER) THEN

      CALL XTIME(0,1,'--- T1EQNS   HV*T / T*HO              ')
C-----------------------------------------------
C S(A,I) = S(A,I) + HV(A,C) * T(C,I)
C-----------------------------------------------
!      CALL CNTRCT ('N','N',NV,NO,NV,A1,HV,T1,A1,S1,NREP)
CPi Use new contraction routine:
      call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),  &
     &    HV,T1,S1,1.0d0,1.0d0,NREP)

C-----------------------------------------------
C S(A,I) = S(A,I) - T(A,K) * HO(K,I)
C-----------------------------------------------
C      CALL CNTRCT ('N','N',NV,NO,NO,-A1,T1,HO,A1,S1,NREP)
CPi Use new contraction routine:
      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),  &
     &    HO,T1,S1,-1.0d0,1.0d0,NREP)
      CALL XTIME(0,2,'--- T1EQNS   HV*T / T*HO            ')

C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
CPi allocate space
      allocate (vooo(IVOOOT(NREP+1)*RCW))
C
      CALL XTIME(0,1,'--- T1EQNS VOOO*TAU                ')
      CALL GETTAU (T1,T2,BUF3)
CPi      CALL T1EQN1 (BUF3,S1)
CPi Extract following three calls from T1EQN1
CPi Get two electron integrals
      CALL GETVOOO (vooo)
CPi The following only applies to FSCC and should be kept
      CALL DELINT ('VOOO','DKDD',vooo,0,0,0)
CPi     ------------------------
CPi     V(K,L,C,I) = V*(C,I,K,L)
CPi     ------------------------
      IF (CARITH) CALL CONJUGA (LVOOOT(NREP+1),vooo,1)
CPi Use new contraction routine
      call contraction_442((/"p2","o1","o2","o3"/),(/"p1","p2","o2","o3"&
CPi Should the following give the same? Crashes due to no available sorting routine.
CPi      call contraction_442((/"o2","o3","p2","o1"/),(/"p1","p2","o2","o3"&
     &          /),(/"p1","o1"/),S1,1.0d0,1.0d0,nrep,LeftTensor=vooo,   &
     &          RightTensor=BUF3)
C
CPi Let contraction routine fetch oovo integrals. This doesn't work
CPi because the routine doesn't know that it should fetch vooo and then
CPi complex conjugate.
CPi      call contraction_442((/"o2","o3","p2","o1"/),(/"p1","p2","o2","o3"&
CPi     &          /),(/"p1","o1"/),S1,1.0d0,1.0d0,nrep,                   &
CPi     &          RightTensor=BUF3)
C
CPi Let contraction routine fetch vooo integrals. Works, but CA test
CPi fails due to the lack of complex conjugation
CPi      call contraction_442((/"p2","o1","o2","o3"/),(/"p1","p2","o2","o3"&
CPi     &          /),(/"p1","o1"/),S1,1.0d0,1.0d0,nrep,                   &
CPi     &          RightTensor=BUF3)
      deallocate (vooo)
      CALL XTIME(0,2,'--- T1EQNS VOOO*TAU               ')

C ****************************************************
C ** END OF TASK SPLITTING. THE COMMON STUFF STARTS **
C ****************************************************
      ENDIF

C--------------------------------------------
C S(A,I) = S(A,I) + V(A,K,C,D) * TAU(C,D,I,K)   ! VOVV contributions !
C--------------------------------------------
      CALL GETTAU (T1,T2,BUF1)
      CALL XTIME(0,1,'--- T1EQNS VOVV contribution  ')
      CALL T1EQN2 (BUF1,BUF3,NBUF3,S1)   ! VOVV contractions, term 6 !
      CALL XTIME(0,2,'--- T1EQNS VOVV contribution  ')

C--------------------------------------
C S(A,I) = S(A,I) - W(A,K,C,I) * T(C,K)    ! VOVO contributions
C--------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOVO * T(C,K)         ')
      CALL GETVOVO (BUF1)
      CALL DELINT ('VOVO','KDDK',BUF1,0,0,0)
CPi old      CALL SRT16 (NREP,MULTB,LTR,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
CPi old     &            BUF1,BUF2)
CPi old      M    = MVO(1)
CPi old      K    = IDIST(1,3,1)
CPi old      OFF1 = IDIST(2,3,1)
CPi old      OFF2 = IDIST(3,3,1)
CPi old      CALL XGEMV ('N',M,K,-A1,BUF2(OFF1),M,T1(OFF2),1,A1,S1,1)
      call contraction_422 ((/"p1","o2","p2","o1"/),(/"p2","o2"/),      &
     & (/"p1","o1"/),T1,S1,-1.0d0,1.0d0,nrep,LeftTensor=BUF1)
CPi test     & (/"p1","o1"/),T1,S1,-1.0d0,1.0d0,nrep)
      CALL XTIME(0,2,'--- T1EQNS VOVO * T(C,K)        ')
            
      END
