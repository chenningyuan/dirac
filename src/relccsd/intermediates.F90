  module intermediates_1b_2b

! Module to make intermediates needed in the CC and Lambda equations and for EOMCC
! Written by Avijit Shee and Lucas Visscher, 2014

  use contraction
  use wrapper_sorting

    implicit none

    public f_bar_mi
    public f_bar_ae   ! these names are subject to be changed if better alternative is found out
!    public f_bar_ia
    public w_ijmn
     private

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

    if (carith) call conjuga (nfvo,fvo,1) 
    call contraction_422 ((/"o1","o3","p1","o2"/),(/"p1","o3"/),(/"o1","o2"/),t1,ho,-1.0d0,1.0d0,nrep)

!***********************TERM-3*********
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!**************************************
   
   call contraction_222 ((/"o1","p1"/),(/"p1","o2"/),(/"o1","o2"/),fvo,t1,ho,1.0d0,1.0d0,nrep)

!   call contraction_222 (fvo, nfvo, t1, ndimt1, ho, 'C', 'N', 'N', &
!                &        no, no, nv, 1.0d0, 1.0d0)
  
   end subroutine

    subroutine f_bar_ae(fvv,buf1,buf2,buf3,t1,t2,hv)

!---------------Calling variables--------------------------------------
      real(8), intent(in) :: FVV(*)
      real(8), intent(in) :: T1(*),T2(*)
      real(8), intent(out):: HV(*)
      real(8), intent(in) :: BUF1(*),BUF2(*),BUF3(*)
!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!---------------Local variables--------------------------------------

     real(8), allocatable :: tau(:) !, vooo(:),vooo_sorted(:),vvoo(:)
!---------------------------------------
! H(A,C) = F(A,C)  
!---------------------------------------
      CALL XCOPY (NFVV,FVV,1,HV,1)
!----------------------------------------------------------------------
! H(A,C) = H(A,C) - W(KL,CD) * T(AD,KL) = H(A,C) - T(A,DKL) * W*(C,DKL)
!----------------------------------------------------------------------
!     CALL GETTAU (T1,T2,BUF3)
!     CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
!    &             BUF3,BUF1)
!     CALL GETVVOO (BUF3)
!     CALL DELINT ('VVOO','DDDD',BUF3,0,0,0)
!     CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
!    &             BUF3,BUF2)
!     CALL CNTRCT ('N','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A1,HV,NREP)

   allocate(tau(ndimt2))
   call gettau  (t1,t2,tau)
   call contraction_442 ((/"o2","o1","p3","p2"/),(/"p3","p1","o2","o1"/),(/"p1","p2"/),tau,HV,-1.0d0,1.0d0,nrep)
   deallocate(tau) 

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
  real(8), intent(in)  :: t1(rcw,ndimt1),t2(rcw,ndimt2)  ! Amplitudes.
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

   call contraction_424 ((/"o1","o2","p1","o3"/),(/"p1","o4"/),(/"o1","o2","o3","o4"/),t1,w_oooo, &
           &      -1.0d0,0.0d0,nrep,nv1)

!******************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) +  V(ij,ef) * Tau(ef,mn)
!******************************************************

   allocate(tau(ndimt2))

   call gettau (t1,t2,tau)

   call contraction_444 ((/"o1","o2","p1","p2"/),(/"p1","p2","o3","o4"/),(/"o1","o2","o3","o4"/),tau,w_oooo,1.0d0,1.0d0,nrep)

   deallocate(tau)

  end subroutine

  end module
