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

    call contraction_444 ((/"o3","o1","p1","p3"/),(/"p2","p3","o3","o2"/),(/"p2","o1","p1","o2"/),tau,w_vovo, &
               &      -1.0d0,0.0d0,nrep)

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

 call srt1l1 (nrep,multb,.false.,nvo,nv,no,noot, &
    &                  iiov,IVOOOT,w_temp,w_ovoo)

 deallocate (w_temp)

 call xscal(nv2*rcw,-A1,w_ovoo,1)

!---------------
! 
!---------------

  call contraction_244 ((/"o1","p2"/),(/"p2","p1","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
           &      -1.0d0,0.0d0,nrep,LeftTensor=hov,RightTensor=t2)

!---------------------------------------------
!   w_vovo(p2o1, p1o2) = V(o3o1, p1p3) * T2(p2p3, o3o2)
!---------------------------------------------

  call contraction_444 ((/"o1","o4","p2","o2"/),(/"p1","p2","o3","o4"/),(/"o1","p1","o2","o3"/),t2,w_ovoo, &
             &      -1.0d0,0.0d0,nrep)

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

  end module
