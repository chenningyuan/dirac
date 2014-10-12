  module intermediates_1b_2b
  use contraction
  use wrapper_sorting

    implicit none

    public f_bar_mi
!    public f_bar_ae   ! these names are subject to be changed if better alternative is found out
!    public f_bar_ia

     private


     contains


    subroutine f_bar_mi(FOO,FVO,EPS,T1,T2,BUF1,BUF2,BUF3,HO)

!---------------Calling variables--------------------------------------
     real(8), intent(in) :: FOO(*),EPS(*), FVO(*)
     real(8), intent(in) :: T1(*),T2(*)
     real(8), intent(out):: HO(*)
     real(8), intent(in) :: BUF1(*),BUF2(*),BUF3(*)
     real(8)  :: v_sort(100000)
     real(8) :: ddot
     integer :: II, I, ND, IRP, N, X, IMIN
     real(8), allocatable :: tau(:), vooo(:),vooo_sorted(:) 
     character            :: vtrans
!     character(*),dimension(*):: v_array, t_array, w_array      
     integer              :: IJKL,J,K,L,JRP,KRP,LRP,KLRP  
!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!#include "ccpar.inc"

!---------------------------------------
! H(K,I) = F(K,I)  - DELTA(K,I) * EPS(I)
!---------------------------------------
      call dzero(ho,nfoo)
      CALL XCOPY (nfoo,FOO,1,HO,1)
!     II = 1
!     I = 1
!     do IRP = 1, NREP
!        N = NO(IRP)
!        ND = (N+1) * RCW
!        CALL DAXPY(N,-A1,EPS(I),1,HO(II),ND)
!        I = I + N
!        II = II + N * N * RCW
!     enddo
      CALL DELFCK ('OO','DK',HO)
!----------------------------------------------------------------------
! H(K,I) = H(K,I) + W(KL,CD) * T(CD,IL) = H(K,I) + W*(CDL,K) * T(CDL,I)
!----------------------------------------------------------------------

   allocate(tau(ndimt2)) 
   CALL GETVVOO (BUF3)
   CALL GETTAU  (T1,T2,tau)

   call contraction_442 (buf3, tau, ho, (/"h3","h1","p1","p2"/),(/"p1","p2","h2","h3"/),(/"h1","h2"/), &
           &        no, no, nvvot, A1, A1,buf2,buf1)

   deallocate(tau) 


!     CALL GETVVOO (BUF3)
!     CALL DELINT ('VVOO','DDDD',BUF3,0,0,0)
!     CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF3,BUF1)
!     CALL GETTAU  (T1,T2,BUF3)
!     CALL SRT1T3(NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,BUF3,BUF2)
!     CALL CNTRCT ('C','N',NO,NO,NVVOT,A1,BUF1,BUF2,A1,HO,NREP)


!---------------------------------------------------------------
!F_bar_mi = F_bar_mi + V(mn, ei) * T1(e,n) --> V(mi,en) * T1(en) 
!---------------------------------------------------------------


   call contraction_422 ((/"o1","o3","o2","p1"/), (/"p1","o3"/), (/"o1","o2"/), T1, HO, -1.0d0, 1.0d0, NREP)

!-----------------------------------------------------------------
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!-----------------------------------------------------------------

  call contraction_222 (FVO, NFVO, T1, ndimt1, HO, 'C', 'N', 'N', &
           &        NO, NO, NV, A1, A1)
  
!      call sort_22_unrestricted((/"o1","o3","o2","p1"/),(/"p1","o3"/),(/"o1","o2"/),vtrans)

!      if (carith) call conjuga (NFVO,FVO,1) 
!      call CNTRCT ('C','N',NO,NO,NV,A1,FVO,T1,A1,HO,NREP)

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
