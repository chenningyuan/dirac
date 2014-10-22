  module intermediates_1b_2b

! Module to make intermediates needed in the CC and Lambda equations and for EOMCC
! Written by Avijit Shee and Lucas Visscher, 2014

  use contraction
  use wrapper_sorting

    implicit none

    public f_bar_mi
!    public f_bar_ae   ! these names are subject to be changed if better alternative is found out
!    public f_bar_ia

     private


     contains


   subroutine f_bar_mi(foo,fvo,t1,t2,ho)
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
     real(8), intent(in) :: foo(rcw,nfoo), fvo(rcw,nfvo)   ! Fock matrix elements
     real(8), intent(in) :: t1(rcw,ndimt1),t2(rcw,ndimt2)  ! Amplitudes
     real(8), intent(out):: ho(rcw,nfoo)                   ! Intermediate
!---------------Local variables---------------------------------------
     real(8), allocatable :: tau(:), vooo(:),vooo_sorted(:),vvoo(:) 
!---------------Common Blocks--------------------------------------

!     Initialize  f_bar_mi with oo block of Fock matrix
      ho = foo

!-----------------------------------------------------
! F_bar_mi(m,i) = F_bar_mi(m,i) + W(mn,cd) * T(cd,in)
!-----------------------------------------------------

   allocate (tau(ndimt2)) 
   allocate (vvoo(ndimt2)) 

   call getvvoo (vvoo)
   call gettau (t1,t2,tau)

   call contraction_442 (vvoo, tau, ho, (/"h3","h1","p1","p2"/),(/"p1","p2","h2","h3"/),(/"h1","h2"/), &
           &             no, no, nvvot, 1.0d0, 1.0d0)

   deallocate (tau) 
   deallocate (vvoo) 

!---------------------------------------------------
!F_bar_mi(m,i) = F_bar_mi(m,i) - V(mn, ie) * T1(e,n) 
!---------------------------------------------------

    call contraction_422 ((/"o1","o3","o2","p1"/), (/"p1","o3"/), (/"o1","o2"/), t1, ho, -1.0d0, 1.0d0, nrep)

!-----------------------------------------------------------------
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!-----------------------------------------------------------------

   call contraction_222 (fvo, nfvo, t1, ndimt1, ho, 'C', 'N', 'N', &
                &        no, no, nv, 1.0d0, 1.0d0)
  
   end subroutine

    subroutine f_bar_ae(FVV,FVO,EPS,T1,T2,BUF1,BUF2,BUF3,HV)

!---------------Description--------------------------------------------
!     Calculates HV intermediate
!---------------Routines called----------------------------------------
!---------------Last modified------------------------------------------
!     Author : Luuk Visscher
!---------------Calling variables--------------------------------------

      real(8), intent(in) :: FVV(*),FVO(*),EPS(*)
      real(8), intent(in) :: T1(*),T2(*)
      real(8), intent(out):: HV(*)
      real(8), intent(in) :: BUF1(*),BUF2(*),BUF3(*)
      real(8) :: ddot
      integer :: IRP, ND, AA, A, N 
!---------------Common Blocks--------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Local variables--------------------------------------

!---------------Executable code--------------------------------------

!---------------------------------------
! H(A,C) = F(A,C)  - DELTA(A,C) * EPS(A)
!---------------------------------------
      CALL XCOPY (NFVV,FVV,1,HV,1)
      AA = 1
!     Calculate start of virtual orbital energies in eps
      A = 1
      DO IRP = 1, NREP
         A = A + NO(IRP)
      ENDDO

      DO IRP = 1, NREP
         N = NV(IRP)
         ND = (N+1) * RCW
         CALL DAXPY(N,-A1,EPS(A),1,HV(AA),ND)
         A = A + N
         AA = AA + N * N * RCW
      ENDDO
      CALL DELFCK ('VV','KD',HV)
!----------------------------------------------------------------------
! H(A,C) = H(A,C) - W(KL,CD) * T(AD,KL) = H(A,C) - T(A,DKL) * W*(C,DKL)
!----------------------------------------------------------------------
      CALL GETTAU (T1,T2,BUF3)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
     &             BUF3,BUF1)
      CALL GETVVOO (BUF3)
      CALL DELINT ('VVOO','DDDD',BUF3,0,0,0)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT, &
     &             BUF3,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,-A1,BUF1,BUF2,A1,HV,NREP)


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

  end module
