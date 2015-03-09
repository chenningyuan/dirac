  module intermediates_1b_2b

! Module to make intermediates needed in the CC and Lambda equations and for EOMCC
! Written by Avijit Shee and Lucas Visscher, 2014

  use contraction
  use wrapper_sorting

    implicit none

    public f_bar_mi
    public f_bar_me
    public f_bar_ae   ! these names are subject to be changed if better alternative is found out
    public w_ijmn
    public w_ejmb
    public w_iemn
    public wbar_mbej
    public w_ejmb_out_of_core
    public f_bar_ae_out_of_core
    public w_iemn_out_of_core
     private

    external gettau

     contains

   subroutine f_bar_mi(foo,fvo,t1,t2,ho)

#include "symm.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in) :: foo(rcw,nfoo), fvo(rcw,nfvo)   ! Fock matrix elements
     real(8), intent(in) :: t1(rcw,ndimt1),t2(rcw,ndimt2)  ! Amplitudes
     real(8), intent(out):: ho(rcw,nfoo)                   ! Intermediate
!+++++++++++++++Local variables++++++++++++++++++++++++++++++++++++++++
     real(8), allocatable :: tau(:) !, vooo(:),vooo_sorted(:),vvoo(:)

!    Initialize  f_bar_mi with oo block of Fock matrix
        ho = foo

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_mi = F_bar_mi - V(nm,ef) * T(ef,in) 
!+++++++++++++++++++++++++++++++++++++++++
   allocate(tau(ndimt2))
   call gettau  (t1,t2,tau)
   call contraction_442 ((/"o3","o1","p1","p2"/),(/"p1","p2","o2","o3"/),(/"o1","o2"/),tau,HO,-1.0d0,1.0d0,nrep)
   deallocate(tau) 

!++++++++++++++++++++++++TERM-2+++++++++++++++++++++
!F_bar_mi(m,i) = F_bar_mi(m,i) - V(mn, ei) * T1(e,n) 
!+++++++++++++++++++++++++++++++++++++++++++++++++++

     call contraction_422 ((/"o1","o3","p1","o2"/),(/"p1","o3"/),(/"o1","o2"/),t1,ho,-1.0d0,1.0d0,nrep)

!***********************TERM-3*********
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!**************************************
   
    if (carith) call conjuga (nfvo,fvo,1) 
   call contraction_222 ((/"o1","p1"/),(/"p1","o2"/),(/"o1","o2"/),fvo,t1,ho,1.0d0,1.0d0,nrep)

!   call contraction_222 (fvo, nfvo, t1, ndimt1, ho, 'C', 'N', 'N', &
!                &        no, no, nv, 1.0d0, 1.0d0)

   end subroutine

   subroutine f_bar_ae(fvv,fvo,t1,t2,hv)

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
      real(8), intent(in)  :: fvv(nfvv*rcw), fvo(nfvo*rcw)
      real(8), intent(in)  :: t1(rcw,ndimt1),t2(rcw,ndimt2)  ! Amplitudes.
      real(8), intent(out) :: hv(nfvv*rcw)
!---------------Local variables--------------------------------------
     real(8), allocatable :: tau(:),hv_temp(:),hv_temp_temp(:) ! vooo(:),vooo_sorted(:),vvoo(:)
     integer              :: i 
!---------------------------------------
! H(A,C) = F(A,C)  
!---------------------------------------

     call xcopy (nfvv,fvv,1,hv,1)

!----------------------------------------------------------------------
! H(A,C) = H(A,C) - W(KL,CD) * T(AD,KL) = H(A,C) - T(A,DKL) * W*(C,DKL)
!----------------------------------------------------------------------

   allocate(tau(ndimt2))
   call gettau (t1,t2,tau)

  allocate (hv_temp(nfvv*rcw))

  call contraction_442 ((/"o1","o2","p2","p3"/),(/"p1","p3","o1","o2"/),(/"p1","p2"/),tau,hv_temp,-1.0d0,0.0d0,nrep)
  
  deallocate(tau) 

  allocate (hv_temp_temp(nfvv*rcw))

  CALL SRT1C1 (NREP,NV,NV,hv_temp,hv_temp_temp) 

  deallocate (hv_temp)

  CALL XAXPY (NFVV,A1,hv_temp_temp,1,hv,1) 

  deallocate (hv_temp_temp)


!---------------------------------------------------------------
!F_bar_ae = F_bar_ae + V(am, ef) * T1(f,m) --> V(ae,fm) * T1(fm) 
!---------------------------------------------------------------



!-----------------------------------------------------------------
!F_bar_ae = F_bar_mi + F(m,e) * T1(e,i)
!-----------------------------------------------------------------

   if (carith) call conjuga (NFVO,FVO,1) 

   call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p1","p2"/),fvo,t1,hv,-1.0d0,1.0d0,nrep)


  end subroutine

   subroutine f_bar_me(fvo,t1,hov)

#include "symm.inc"
#include "param.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in),target  :: fvo(rcw,nfvo)                  ! Fock matrix elements
     real(8), intent(in) :: t1(rcw,ndimt1)                 ! Amplitudes
     real(8), intent(out):: hov(rcw,nfvo)                  ! Intermediate
!+++++++++++++++Local variables++++++++++++++++++++++++++++++++++++++++
     integer          :: i 
!+++++initialization of hov+++++++++++++
 
 call srt1c1(nrep,nv,no,fvo,hov)

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_me = F_bar_me + V(mn,ef) * T1(f,n) 
!+++++++++++++++++++++++++++++++++++++++++

 call contraction_422 ((/"o1","o2","p1","p2"/),(/"p2","o2"/),(/"o1","p1"/),t1,hov,-1.0d0,1.0d0,nrep)

  end subroutine

  subroutine w_ijmn(t1,t2,w_oooo)

!---------------Description--------------------------------------------

!     Calculates w_oooo intermediate

!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw,ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_oooo(rcw,nv1)                ! Intermediate.
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
!--------------------------------------------------

!--------------------
! initialize w_oooo
!--------------------

   call getoooo (w_oooo)

!*************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) + (V(ij,en)*T(e,m))
!*************************************************
!whether the anti-symmetrization is missing? 

   call contraction_424 ((/"o1","o2","p1","o3"/),(/"p1","o4"/),(/"o1","o2","o3","o4"/),w_oooo, &
           &      -1.0d0,0.0d0,nrep,RightTensor=t1)

!******************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) +  V(ij,ef) * Tau(ef,mn)
!******************************************************

   allocate(tau(ndimt2))

   call gettau (t1,t2,tau)

   call contraction_444 ((/"o1","o2","p1","p2"/),(/"p1","p2","o3","o4"/),(/"o1","o2","o3","o4"/),tau,w_oooo,1.0d0,1.0d0,nrep)

   deallocate(tau)

  end subroutine

 subroutine w_ejmb(t1,t2,w_vovo)

 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw,ndimt2)  ! Amplitudes.
   real(8), intent(out) :: w_vovo(rcw,nv4)                ! Intermediate.
 !---------------Local Variables----------------------------------------
   real(8), allocatable :: tau(:)
 !--------------------
 ! initialize w_vovo
 !--------------------

  write(*,*)'intermediate w_ejmb starts here'    

  call getvovo (w_vovo)

!---------------------------------------------
!   w_vovo(p2o1, p1o2) = V(o3o1,p1o2)*T1(p2o3)
!---------------------------------------------
    call contraction_424 ((/"o3","o1","p1","o2"/),(/"p2","o3"/),(/"p2","o1","p1","o2"/),w_vovo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t1)


!  call contraction_424 ((/"o3","o1","o2","p1"/),(/"p2","o3"/),(/"p2","o1","o2","p1"/),t1,w_vovo, &
!           &      -1.0d0,0.0d0,nrep)
!---------------------------------------------
!   w_vovo(p2o1, p1o2) = V(o3o1, p1p3) * T2(p2p3, o3o2)
!---------------------------------------------

  allocate (tau(ndimt2))

  call gettau (t1,t2,tau,'approximate')

!   call contraction_444 ((/"o3","o1","p1","p3"/),(/"p2","p3","o3","o2"/),(/"p2","o1","p1","o2"/),tau,w_vovo, &
!              &      -1.0d0,0.0d0,nrep)

! call contraction_444 ((/"o3","o1","p3","p1"/),(/"p2","p3","o2","o3"/),(/"p2","o1","o2","p1"/),tau,w_vovo, &
!            &      -1.0d0,0.0d0,nrep)

  deallocate (tau)

  end subroutine

 subroutine w_iemn(t1,t2,w_ovoo,wbar_ovvo,w_oooo,hov)

 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(out) :: w_ovoo(rcw,nv2)                ! Intermediate.
   real(8), intent(in)  :: wbar_ovvo(rcw*nv4),w_oooo(rcw*nv1),hov(nfvo*rcw)         ! reuse of intermediates.
 !---------------Local Variables----------------------------------------

 real*8, allocatable   :: w_temp(:) 

 !--------------------
 ! initialize wbar_ovvo
 !--------------------
 allocate (w_temp(nv2*rcw))

  call getvooo (w_temp) 
  w_ovoo = 0.0d0

 call srt1l1 (nrep,multb,.false.,nvo,nv,no,noot,iiov,IVOOOT,w_temp,w_ovoo)

 deallocate (w_temp)


 call xscal(nv2*rcw,-A1,w_ovoo,1)

!---------------
! w_ovoo(o1p1,o3o2) = F_me(o1p2) * T2(p2p1o2o3)  
!---------------

  call contraction_244 ((/"o1","p2"/),(/"p2","p1","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
           &      -1.0d0,0.0d0,nrep,LeftTensor=hov,RightTensor=t2)

!---------------------------------------------
!   w_ovoo(o1p1,o2o3) = V(o1o4,p2o2) * T2(p1p2, o3o4)
!---------------------------------------------

  call contraction_444 ((/"o1","o4","p2","o2"/),(/"p1","p2","o3","o4"/),(/"o1","p1","o2","o3"/),t2,w_ovoo, &
             &      -1.0d0,0.0d0,nrep)

!-------------------------------------------
!   w_ovoo(o1p1,o3o2) = V(o1p1,p2o3)*T1(p2o2)
!--------------------------------------------
!careful. the following diagram has its permuted label in ket. this is not the right way of constructing this intermediate. i will come back to it soon.
   call contraction_424 ((/"o1","p1","p2","o3"/),(/"p2","o2"/),(/"o1","p1","o3","o2"/),w_ovoo, &
              &      1.0d0,0.0d0,nrep,LeftTensor=wbar_ovvo,RightTensor=t1)


! write(*,*)'just before the problem'    

!  call contraction_424 ((/"o1","o4","o2","o3"/),(/"p1","o4"/),(/"o1","p1","o2","o3"/),w_ovoo, &
!            &      1.0d0,0.0d0,nrep,LeftTensor=w_oooo,RightTensor=t1)


! call contraction_424 ((/"o2","o3","o1","o4"/),(/"p1","o4"/),(/"o2","o3","o1","p1"/),w_ovoo, &
!           &      1.0d0,0.0d0,nrep,LeftTensor=w_oooo,RightTensor=t1)

  end subroutine

 subroutine wbar_mbej(t2,wbar_ovvo)

 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(out) :: wbar_ovvo(rcw,nv4)                ! Intermediate.
 !---------------Local Variables----------------------------------------

  real*8, allocatable   :: w_vovo(:) 

  write(*,*)'mbej starts here'

 !--------------------
 ! initialize wbar_ovvo
 !--------------------
  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call SRT1L1 (NREP,MULTB,.FALSE.,NVO,NV,NO,NVO, &
     &                   IVOVO,IIOV,w_vovo,wbar_ovvo) 

  deallocate (w_vovo)

  call xscal(nv4*rcw,-A1,wbar_ovvo,1)
!---------------------------------------------
!   w_vovo(p2o1, p1o2) = V(o3o1, p1p3) * T2(p2p3, o3o2)
!---------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),t2,wbar_ovvo, &
             &      -1.0d0,0.0d0,nrep)

  end subroutine


 subroutine w_ejmb_out_of_core(t1,t2,buf1,buf2,buf3,nbuf3,w_vovo)


!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vovo(rcw,nv4)              ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*),buf3(*)        ! reuse of intermediates.
   integer, intent(in)    :: nbuf3 

!---------------Common Blocks--------------------------------------   
                                                                     
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
!#include "eqns.inc"                                                   
#include "ccpar.inc"                                                  
!
!---------------Local variables--------------------------------------
!                                                                    
      REAL*8 T1VAL(2)                                                
      DATA T1VAL /2*0.0D0/                                            
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ci,akci,t1di,di,akdi,t1ci,cmin,c,d,i
      integer :: crep,drep,irp,mint
!------------------------------------------------------------------
! H(AK,CI) = H(AK,CI) - W(AK,CD) * T(D,I)
! H is kept in BUF1, still ordered VOVO
!------------------------------------------------------------------

      IF(MYPROC.EQ.MASTER) THEN
        CALL xcopy(nv4,w_vovo,0,BUF1,1)
      else 
        CALL XCOPY(NV4,A0,0,BUF1,1)
      endif

      CALL XTIME(0,1,'--- HINTM: VOVV*T              ')
      DO 30 IRP = 1, NREP
         IF (NVO(IRP).EQ.0) GOTO 30
         DONE = .FALSE.
         ISTART = 0
         CCO=0
         MINT = NBUF3/NVO(IRP)

! ** reading from GETVOVV: ISTART is updated with the actual batch #
! ** available on the local node, NINT is updated with the number of
! ** actually read integrals from this batch and is a pure output parameter !
! ** We therefore have to count upwards in case of ISTART.NE.0 !!

         CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
         CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)

! ** now start SORTING and contraction loop

         CD = 0
         DO 20 DREP = 1, NREP
           CREP = MULTB(DREP,IRP+NREP,2)
           IF (CREP.LT.DREP) GOTO 20
           DO 15 D = 1, NV(DREP)
              CMIN = 1
              IF (CREP.EQ.DREP) CMIN = D + 1
              DO 10 C = CMIN, NV(CREP)

                 IF(CCO.LT.ISTART) THEN
                   CCO=CCO+1
                   GOTO 10
                 ENDIF

                 CD = CD + 1
                 IF (CD.GT.NINT) THEN
!                   ---------------------------------------
!                   We need the next buffer in this IRREP !
!                   ---------------------------------------
                    ISTART = ISTART + NINT
                    IF(.NOT.DONE) THEN
                      CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF3,MINT)
                      CALL DELINT ('VOVV','KDDD',BUF3,IRP,ISTART,NINT)
                      CD = 1
                      CCO=ISTART
                    ELSE
                      GOTO 30
                    ENDIF
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
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE
!------------------------------------------------------------------
! Now sort H(AK,CI) to H(AI,CK)  is to be done in any case !!
! The H buffer is not changed any more up to the final XGEMM !
! Every other calculation is done in BUF1 !!
!------------------------------------------------------------------
   CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO, &
   &            BUF1,w_vovo)
   CALL XTIME(0,2,'--- HINTM: VOVV*T              ')

 end subroutine


   subroutine f_bar_ae_out_of_core (t1,buf2,nbuf2,hv)

!---------------Calling variables--------------------------------------
     real*8, intent(in)    :: t1(rcw,nfvo),buf2(*)
     real*8, intent(inout) :: hv(rcw*nfvv)
     integer,intent(in)    :: nbuf2 
!---------------Common Blocks--------------------------------------

#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ccpar.inc"

!---------------Local variables----
     logical :: USEDZ,RIGHT
     integer :: m,n
!----------------------------------
      if(myproc==master) then
        call xcopy(nfvv*rcw,hv,0,buf2,1)
      else 
        call xcopy(nfvv,A0,0,buf2,1)
      endif
!------------------------------------
! G(A,C) = G(A,C) + W(AK,CD) * T(D,K)
!------------------------------------
      M = MVV(1)
      N = MVO(1)
      IF (N.EQ.0) RETURN
      USEDZ = .FALSE.
      RIGHT = .TRUE.
      CALL SRT20D(NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV, &
     &            JJVO,BUF2,NBUF2,T1,hv,USEDZ,RIGHT)

   end subroutine

   subroutine w_iemn_out_of_core (t1,t2,buf1,buf2,nbuf3,w_ovoo)

!
!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_ovoo(rcw,nv2)              ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*)                      ! reuse of intermediates.
   integer, intent(in)    :: nbuf3 

!---------------Common Blocks--------------------------------------   
                                                                     
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
!#include "eqns.inc"                                                   
#include "ccpar.inc"                                                  

!---------------Local variables--------------------------------------
                                                                    
      logical :: done,teq
      integer :: kd,istart,nint,m,n,k
      integer :: irp,mint,off1,off2
      real*8,allocatable :: tau(:)


    call xcopy(nv2,A0,0,buf2,1)


   allocate(tau(ndimt2))
   call gettau  (t1,t2,tau)

    OFF1 = 1
    OFF2 = 1
   DO 130 IRP = 1, NREP
    ISTART = 0
     IF (NVVT(IRP).EQ.0.OR.NVO(IRP).EQ.0) GOTO 130
      M = NVO (IRP)
         N = NOOT(IRP)
          KD = NVVT(IRP)
           MINT = NBUF3/NVO(IRP)

  110      CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF1,MINT)

           K = NINT

 ! ** now we construct the full NVO x NOOT matrix
 ! ** but the contributions are only from the restricted CD range
 ! ** the method is additive since ZGEMM works according to
 ! **      C=A*B + C  ==> all the already computed contributions
 ! ** will be added in this IRREP.

           CALL XGEMM ('N','N',M,N,K,A1,BUF1,M,tau(OFF1+ISTART*RCW), &
            &               KD,A1,BUF2(OFF2),M)
           IF (.NOT.DONE) THEN
           ISTART = ISTART + NINT
           GOTO 110
           ENDIF
           OFF1 = OFF1 + KD * N * RCW
           OFF2 = OFF2 + M * N * RCW
 130  CONTINUE

    deallocate (tau)

      CALL SRT1L1 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF1)

      call xaxpy (nv2,A1,buf1,1,w_ovoo,1) 

  end subroutine

  end module
