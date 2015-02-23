  module intermediates_1b_2b

! Module to make intermediates needed in the CC and Lambda equations and for EOMCC
! Written by Avijit Shee and Lucas Visscher, 2014

  use contraction
  use wrapper_sorting

    implicit none

    public f_bar_mi
    public f_bar_ae   ! these names are subject to be changed if better alternative is found out
!   public f_bar_ia
    public w_ijmn
    public w_ejmb
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

   subroutine f_bar_ae(fvv,t1,t2,hv)

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
      real(8), intent(in)  :: fvv(nfvv*rcw)
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

!CALL GETVVOO (BUF3)

!CALL SRT1C1 (NREP,NVVT,NOOT,BUF3,BUF2) 

!CALL SRT22 (NREP,MULTB,LFA,NOOT,NV,NV,NVOOT,LVOOVT,LLVOOT,BUF2,BUF1)

!CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
!&             tau,BUF2)

!CALL CNTRCT ('C','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A0,hv,NREP)

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

!!!   call getvooo (BUF3)
!!!   call delint ('VOOO','DDDD',BUF3,0,0,0)
!!!   call SRT1C1 (NREP,NVO,NOO,BUF3,BUF2)
!!!   call SRT16 (NREP,MULTB,LFA,LFA,NO,NO,NV,NO,MVO,JVOOO,JJOO,JJVO,&
!!!  &            BUF2,BUF1)
!!!   call XGEMV ('N',NOO,NVO,A1,BUF1,NOO,T1,1,A1,HO,1)
!-----------------------------------------------------------------
!F_bar_ae = F_bar_mi + F(m,e) * T1(e,i)
!-----------------------------------------------------------------

!      if (carith) call conjuga (NFVO,FVO,1) 
!      call CNTRCT ('C','N',NO,NO,NV,A1,FVO,T1,A1,HO,NREP)

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

!  call contraction_444 ((/"o3","o1","p3","p1"/),(/"p2","p3","o2","o3"/),(/"p2","o1","o2","p1"/),tau,w_vovo, &
!             &      -1.0d0,0.0d0,nrep)

  deallocate (tau)

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

 !--------------------
 ! initialize wbar_ovvo
 !--------------------
  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call SRT1L1 (NREP,MULTB,.FALSE.,NVO,NV,NO,NVO, &
     &                   NV4,JJOV,w_vovo,wbar_ovvo)

  deallocate (w_vovo)
!---------------------------------------------
!   w_vovo(p2o1, p1o2) = V(o3o1, p1p3) * T2(p2p3, o3o2)
!---------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),t2,wbar_ovvo, &
             &      -1.0d0,0.0d0,nrep)

  end subroutine

  end module
