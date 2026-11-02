  module intermediates_1b_2b

! Module to make intermediates needed in the CC and Lambda equations and for EOMCC
! Written by Avijit Shee and Lucas Visscher, 2014

  use interface_to_mpi
  use contraction
  use interface_to_mpi
  use wrapper_sorting

    implicit none

    public f_bar_mi
    public f_bar_mi_MBPT2
    public f_bar_me
    public f_bar_me_MBPT2
    public f_bar_ae   
    public f_bar_ae_MBPT2
    public w_ijmn
    public w_ijmn_MBPT2
    public w_ejmb 
    public w_ejmb_MBPT2 
    public w_mbej  
    public w_mbej_MBPT2
    public w_iemn
    public w_iemn_MBPT2
    public w_efam
    public w_efam_MBPT2
    public w_amef
    public w_amef_MBPT2
    public w_mnie
    public w_mnie_MBPT2
    public wbar_mbej
    public wbar_mbej_MBPT2 
    public w_ejmb_out_of_core  
    public f_bar_ae_out_of_core
    public w_iemn_out_of_core
    public w_efam_diagram2
    public w_efam_diagram2_MBPT2 
    public w_efab_out_of_core
    public w_efab_out_of_core_MBPT2 
    public w_efab_oocore_wefam_d2_MBPT2  
    public G_intm
     private

     contains

   subroutine f_bar_mi_MBPT2(foo,fvo,t1,t2,ho_MBPT2)

!      use symm_mod, only: NDIMT1, NDIMT2, NFOO, NFVO, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!--------------Calling variables---------------------------------------
     real(8), intent(in) :: foo(nfoo*rcw), fvo(rcw*nfvo)   ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes
     real(8), intent(inout):: ho_MBPT2(nfoo*rcw)           ! Intermediate
!---------------Local variables----------------------------------------
!!!!!!!!!!!!     real(8), allocatable :: tau(:) !, vooo(:),vooo_sorted(:),vvoo(:)

!    Initialize  f_bar_mi with oo block of Fock matrix
       ! Fbar_mi=foo
       ho_MBPT2(1:nfoo*rcw) = foo(1:nfoo*rcw)

!-------------Term-1----------------------
! F_bar_mi_MBPT2 = F_bar_mi_MBPT2 - V(nm,ef) * t2(ef,in) 
!-----------------------------------------
      ! sum einstein nef  ! Vnmef: imef         t2 : efmj               F_bar_mi : ij
      ! sum einstein nef  ! Vnmef: nmef         t2 : efin               F_bar_mi : mi
   call contraction_442 ((/"o3","o2","p1","p2"/),(/"p1","p2","o1","o3"/),(/"o2","o1"/),ho_MBPT2,&
  & -1.0d0,1.0d0,nrep,RightTensor=t2)

!------------------------TERM-2---------------------
!F_bar_mi_MBPT2(m,i) = F_bar_mi_MBPT2(m,i) - V(mn, ei) * t1(e,n) 
!---------------------------------------------------
     ! sum einstein ne    ! Vmnei: imje(miej)       t1: em         F_bar_mi : ij
     ! sum einstein ne    ! Vmnei: mnei             t1: en         F_bar_mi : mi
    call contraction_422 ((/"o2","o3","p1","o1"/),(/"p1","o3"/),(/"o2","o1"/),t1,ho_MBPT2,&
   & -1.0d0,1.0d0,nrep)
!!!!!!!!!!!! Le 1609 : je remrque l'inversion p1 o1 de V qui explique le signe...
!-----------------------TERM-3---------
!F_bar_mi = F_bar_mi + f(m,e) * t1(e,i)
!--------------------------------------
    ! sum einstein e     ! fov: ie        t1: je      Fbar_mi : ij
    ! sum einstein e     ! fov: me        t1: ei      Fbar_mi : mi     
   call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),fvo,t1,ho_MBPT2,&
  & 1.0d0,1.0d0,nrep)

   end subroutine

   subroutine f_bar_mi(foo,fvo,t1,t2,ho)

!      use symm_mod, only: NDIMT1, NDIMT2, NFOO, NFVO, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!--------------Calling variables---------------------------------------
     real(8), intent(in) :: foo(nfoo*rcw), fvo(rcw*nfvo)   ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes
     real(8), intent(inout):: ho(nfoo*rcw)                 ! Intermediate
!---------------Local variables----------------------------------------
     real(8), allocatable :: tau(:) !, vooo(:),vooo_sorted(:),vvoo(:)

!    Initialize  f_bar_mi with oo block of Fock matrix
       ho(1:nfoo*rcw) = foo(1:nfoo*rcw)
!-------------Term-1----------------------
! F_bar_mi = F_bar_mi - V(nm,ef) * T(ef,in) 
!-----------------------------------------
   allocate(tau(ndimt2*rcw))
   call gettau  (t1,t2,tau)
   call contraction_442 ((/"o3","o2","p1","p2"/),(/"p1","p2","o1","o3"/),(/"o2","o1"/),ho,-1.0d0,1.0d0,nrep,RightTensor=tau)
   deallocate(tau) 
!------------------------TERM-2---------------------
!F_bar_mi(m,i) = F_bar_mi(m,i) - V(mn, ei) * T1(e,n) 
!---------------------------------------------------
    call contraction_422 ((/"o2","o3","p1","o1"/),(/"p1","o3"/),(/"o2","o1"/),t1,ho,-1.0d0,1.0d0,nrep)
!-----------------------TERM-3---------
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!--------------------------------------
   call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),fvo,t1,ho,1.0d0,1.0d0,nrep)
   end subroutine

   subroutine f_bar_ae_MBPT2(fvv,fvo,t1,t2,hv_MBPT2)

!---------------Common Blocks--------------------------------------
!      use symm_mod, only: NDIMT1, NDIMT2, NFVO, NFVV, NO, NREP, NV
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
      real(8), intent(in)  :: fvv(nfvv*rcw), fvo(nfvo*rcw)
      real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
      real(8), intent(inout) :: hv_MBPT2(nfvv*rcw)
!---------------Local variables--------------------------------------
      real(8), allocatable :: fov(:)
!---------------------------------------
! hv_MBPT2(a,e) = fvv(a,e)  
!---------------------------------------
    !Fbar_ab = fab
    !Fbar_ae = fae
    hv_MBPT2 = fvv

!----------------------------------------------------------------------
! hv_MBPT2(a,e) = hv(a,e) - V(mn,ef) * t2(af,mn)
!----------------------------------------------------------------------
      !sum einsteim mne  !  Vmnbe: mnbe           t2afmn: aemn          Fbar_ae: ab
      !sum einstein fmn  !  Vmnef: mnef           t2afmn: afmn          Fbar_ae: ae 
   call contraction_442 ((/"o1","o2","p2","p3"/),(/"p1","p3","o1","o2"/),(/"p1","p2"/),hv_MBPT2,-1.0d0,1.0d0,nrep,RightTensor=t2)
  
!---------------------------------------------------------------
!F_bar_ae_MBPT2 = F_bar_ae_MBPT2 + V(am, ef) * t1(f,m) -> V(ae,fm) * t1(fm) 
!---------------------------------------------------------------
     !sum einstein mf    !  Vambf:  ambf         t1: fm       Fbar_ae: ab
     !sum einstein mf    !  Vmaef:  amef         t1: fm       Fbar_ae: ae
  call contraction_422 ((/"p1","o1","p2","p3"/),(/"p3","o1"/),(/"p1","p2"/),t1,hv_MBPT2,1.0d0,1.0d0,nrep)

!-----------------------------------------------------------------
!F_bar_ae_MBPT2 = F_bar_ae_MBPT2 + F(m,a) * T1(e,m)
!-----------------------------------------------------------------

   allocate (fov(nfvo*rcw))
   call srt1c1(nrep,nv,no,fvo,fov) !SUBROUTINE SRT1C1 (NREP,NPAIR1,NPAIR2,BUF1,BUF2) in src/relccsd/ccsort.F
     !sum einstein       !  fov: ma       t1: bm      Fbar_ae: ba
     !sum einstein       !  fov: ma       t1: em      Fbar_ae: ea   !!!!!!!! INVERSION  !!!!!!!   
   call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p1","p2"/),fov,t1,hv_MBPT2,-1.0d0,1.0d0,nrep)

   deallocate (fov)

  end subroutine

   subroutine f_bar_ae(fvv,fvo,t1,t2,hv)

!---------------Common Blocks--------------------------------------
!      use symm_mod, only: NDIMT1, NDIMT2, NFVO, NFVV, NO, NREP, NV
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
      real(8), intent(in)  :: fvv(nfvv*rcw), fvo(nfvo*rcw)
      real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
      real(8), intent(inout) :: hv(nfvv*rcw)
!---------------Local variables--------------------------------------
      real(8), allocatable :: tau(:)
      real(8), allocatable :: fov(:)
!---------------------------------------
! hv(a,e) = fvv(a,e)  
!---------------------------------------
    hv = fvv
      ! hv = 0.0
    

!----------------------------------------------------------------------
! hv(a,e) = hv(a,e) - V(mn,ef) * Tau(af,mn)
!----------------------------------------------------------------------
   allocate(tau(ndimt2*rcw))
   call gettau (t1,t2,tau)
   call contraction_442 ((/"o1","o2","p2","p3"/),(/"p1","p3","o1","o2"/),(/"p1","p2"/),hv,-1.0d0,1.0d0,nrep,RightTensor=tau)
   deallocate(tau) 
!---------------------------------------------------------------
!F_bar_ae = F_bar_ae + V(am, ef) * T1(f,m) -> V(ae,fm) * T1(fm) 
!---------------------------------------------------------------
  call contraction_422 ((/"p1","o1","p2","p3"/),(/"p3","o1"/),(/"p1","p2"/),t1,hv,1.0d0,1.0d0,nrep)
!-----------------------------------------------------------------
!F_bar_ae = F_bar_ae + F(m,a) * T1(e,m)
!-----------------------------------------------------------------
   allocate (fov(nfvo*rcw))
   call srt1c1(nrep,nv,no,fvo,fov)
   call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p1","p2"/),fov,t1,hv,-1.0d0,1.0d0,nrep)
   deallocate (fov)
  end subroutine

  subroutine f_bar_me_MBPT2(fvo,t1,hov_MBPT2)

!      use symm_mod, only: NDIMT1, NFVO, NO, NREP, NV
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in),target  :: fvo(rcw*nfvo)          ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1)                 ! Amplitudes
     real(8), intent(out):: hov_MBPT2(rcw*nfvo)                  ! Intermediate
!+++++initialization of hov+++++++++++++
 
  call srt1c1(nrep,nv,no,fvo,hov_MBPT2)
!  hov_MBPT2 = fia
!  hov_MBPT2 = fvo

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_me_MBPT2 = F_bar_me_MBPT2 + V(mn,ef) * T1(f,n) 
!+++++++++++++++++++++++++++++++++++++++++
   !sum einstein em     ! V mnef : imae         t2: em      Fbar_me : ia
   !sum einstein nf     ! V mnef : mnef         t1: fn      Fbar_me : me
 call contraction_422 ((/"o1","o2","p1","p2"/),(/"p2","o2"/),(/"o1","p1"/),t1,hov_MBPT2,1.0d0,1.0d0,nrep)

  end subroutine


  subroutine f_bar_me(fvo,t1,hov)

!      use symm_mod, only: NDIMT1, NFVO, NO, NREP, NV
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in),target  :: fvo(rcw*nfvo)          ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1)                 ! Amplitudes
     real(8), intent(out):: hov(rcw*nfvo)                  ! Intermediate
!+++++initialization of hov+++++++++++++
 
  call srt1c1(nrep,nv,no,fvo,hov)
!  hov = fvo

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_me = F_bar_me + V(mn,ef) * T1(f,n) 
!+++++++++++++++++++++++++++++++++++++++++
 call contraction_422 ((/"o1","o2","p1","p2"/),(/"p2","o2"/),(/"o1","p1"/),t1,hov,1.0d0,1.0d0,nrep)
  end subroutine

  subroutine w_ijmn_MBPT2(t1,t2,w_oooo_MBPT2)

!---------------Description--------------------------------------------
!     Calculates w_oooo intermediate
!-----------------------------------------------------------------------

!      use symm_mod, only: NDIMT1, NDIMT2, NREP, NV1
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_oooo_MBPT2(rcw*nv1)                ! Intermediate.
!---------------Local Variables--------------------
!--------------------
! initialize w_oooo_MBPT2
!--------------------
        ! Wijmn Wijkl 
   call getoooo (w_oooo_MBPT2)
!*************************************************
! w_oooo_MBPT2(ij,mn) = w_oooo_MBPT2(ij,mn) + (V(ij,en)*t1(e,m))
!*************************************************
    !sum einstein e (29)! V : ijel                t1: ek         W_ijke : ijlk   !!!! INVERSION Attention au signe !!!!
    !sum einstein e     ! V : ijen                t1: em         W_ijmn : ijnm   !!!! INVERSION des INDICES !!!!
   call contraction_424 ((/"o1","o2","p1","o4"/),(/"p1","o3"/),(/"o1","o2","o4","o3"/),w_oooo_MBPT2, &
  &      1.0d0,-1.0d0,nrep,RightTensor=t1)
!******************************************************
! w_oooo_MBPT2(ij,mn) = w_oooo_MBPT2(ij,mn) +  V(ij,ef) * t2(ef,mn)
!******************************************************
    ! sum einstein ef (29) ! Voopp:  ijef        t2: efkl                Wijkl : ijkl
    ! sum einstein         ! Vijef:  ijef        (tau): efkl             Wijmn : ijmn      
   call contraction_444 ((/"o1","o2","p1","p2"/),(/"p1","p2","o3","o4"/),(/"o1","o2","o3","o4"/),w_oooo_MBPT2,1.0d0,1.0d0,nrep, &
  &                      RightTensor=t2)
  end subroutine
 
  subroutine w_ijmn(t1,t2,w_oooo)

!---------------Description--------------------------------------------
!     Calculates w_oooo intermediate
!-----------------------------------------------------------------------

!      use symm_mod, only: NDIMT1, NDIMT2, NREP, NV1
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_oooo(rcw*nv1)                ! Intermediate.
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
!--------------------------------------------------

!--------------------
! initialize w_oooo
!--------------------

   call getoooo (w_oooo)

!*************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) + (V(ij,en)*t1(e,m))
!*************************************************
    !sum einstein  e    ! V : ijen          t1: em      W_ijmn : ijmn   !!!! INVERSION des INDICES !!!!
   call contraction_424 ((/"o1","o2","p1","o4"/),(/"p1","o3"/),(/"o1","o2","o4","o3"/),w_oooo, &
  &      1.0d0,-1.0d0,nrep,RightTensor=t1)

!******************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) +  V(ij,ef) * Tau(ef,mn)
!******************************************************

   allocate(tau(ndimt2*rcw))
   call gettau (t1,t2,tau)
   call contraction_444 ((/"o1","o2","p1","p2"/),(/"p1","p2","o3","o4"/),(/"o1","o2","o3","o4"/),w_oooo,1.0d0,1.0d0,nrep, &
  &                      RightTensor=tau)

   deallocate(tau)

  end subroutine

!===============================================================================================
 subroutine w_ejmb_MBPT2(t1,t2,w_voov_MBPT2)
!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------
!      use symm_mod, only: IOVOOT, IVOOV, NDIMT1, NDIMT2, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_voov_MBPT2(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real*8,allocatable      :: array_temp(:)
!--------------------
! initialize w_voov
!--------------------
 call getvoov (w_voov_MBPT2)
!---------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) - V(nj,mb)*T1(e,n)
!---------------------------------------------------
   call contraction_244 ((/"p2","o3"/),(/"o3","o1","o2","p1"/),(/"p2","o1","o2","p1"/),w_voov_MBPT2, &
            &      -1.0d0,1.0d0,nrep,LeftTensor=t1)
!----------------------------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) + V(nj,fb)*(T2(ef,mn)-T1(f,m)*T1(e,n)
!----------------------------------------------------------------------
  call contraction_444 ((/"p2","p3","o2","o3"/),(/"o3","o1","p3","p1"/),(/"p2","o1","o2","p1"/),w_voov_MBPT2, &
   &      1.0d0,0.0d0,nrep,LeftTensor=t2)
  allocate(array_temp(iovoot(nrep+1)*rcw))
   array_temp = 0.0d0
  call contraction_424 ((/"o3","o1","p1","p3"/),(/"p3","o2"/),(/"o3","o1","o2","p1"/),array_temp, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)
  call contraction_244 ((/"p2","o3"/),(/"o3","o1","o2","p1"/),(/"p2","o1","o2","p1"/),w_voov_MBPT2, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=t1,RightTensor=array_temp)
  deallocate(array_temp)
!--------------------------------------------------
!  w_voov(ej,mb) = w(ej,mb) + V(ej,fb)*t(f,m) 
!--------------------------------------------------
  call contraction_424 ((/"p2","o1","p1","p3"/),(/"p3","o2"/),(/"p2","o1","o2","p1"/),w_voov_MBPT2, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)
 end subroutine
!===========================================================================================================


 
!===============================================================================================
 subroutine w_ejmb(t1,t2,w_voov)
!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------
!      use symm_mod, only: IOVOOT, IVOOV, NDIMT1, NDIMT2, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_voov(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real*8,allocatable      :: array_temp(:)
!--------------------
! initialize w_voov
!--------------------
 call getvoov (w_voov)
!---------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) - V(nj,mb)*T1(e,n)
!---------------------------------------------------
   call contraction_244 ((/"p2","o3"/),(/"o3","o1","o2","p1"/),(/"p2","o1","o2","p1"/),w_voov, &
            &      -1.0d0,1.0d0,nrep,LeftTensor=t1)
!----------------------------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) + V(nj,fb)*(T2(ef,mn)-T1(f,m)*T1(e,n)
!----------------------------------------------------------------------
  call contraction_444 ((/"p2","p3","o2","o3"/),(/"o3","o1","p3","p1"/),(/"p2","o1","o2","p1"/),w_voov, &
   &      1.0d0,0.0d0,nrep,LeftTensor=t2)
  allocate(array_temp(iovoot(nrep+1)*rcw))
   array_temp = 0.0d0
  call contraction_424 ((/"o3","o1","p1","p3"/),(/"p3","o2"/),(/"o3","o1","o2","p1"/),array_temp, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)
  call contraction_244 ((/"p2","o3"/),(/"o3","o1","o2","p1"/),(/"p2","o1","o2","p1"/),w_voov, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=t1,RightTensor=array_temp)
  deallocate(array_temp)
!--------------------------------------------------
!  w_voov(ej,mb) = w(ej,mb) + V(ej,fb)*t(f,m) 
!--------------------------------------------------
  call contraction_424 ((/"p2","o1","p1","p3"/),(/"p3","o2"/),(/"p2","o1","o2","p1"/),w_voov, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)
 end subroutine
!===========================================================================================================

subroutine w_mbej_MBPT2(t1,t2,w_ovvo_MBPT2)
!      use symm_mod, only: IIOV, IOVVO, IVOOV, MULTB, NDIMT1, NDIMT2, NO, NREP, NV, NV4, NVO
!      use complex_mod, only: RCW
  external srt1l1
!---------------Description--------------------------------------------
!     Calculates w_voov_MBPT2 intermediate : Gwaltney eq.34
!-----------------------------------------------------------------------
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_ovvo_MBPT2(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real(8), allocatable :: w_vovo(:),array_temp(:)
!--------------------
! initialize w_voov  : ! De la forme <ai||jb>
!--------------------
  allocate (w_vovo(nv4*rcw))
  call getvovo (w_vovo)
  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,w_ovvo_MBPT2) 
  deallocate (w_vovo)
!---------------------------------------------------
!   w_ovvo_MBPT2(mb,ej) = w_voov(mb,ej) - V(mn,ej)*T1(b,n)  : Probleme d ecriture??? Si non alors aijb==iabj
! G eq 34 : 3eme terme, de la forme Voovo*t1vo 
!---------------------------------------------------
    !sum einstein    !  Vimbj             t1: am           iabj   11 Attention G :aijb  !!!
    !sum einstein n  !  Vmnej             t1: bn           mbej 
   call contraction_424 ((/"o2","o3","p2","o1"/),(/"p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo_MBPT2, &
            &      -1.0d0,1.0d0,nrep,RightTensor=t1)
!-------------------------------------------------------------------------------------
!  w_ovvo_MBPT2(mb,ej) = w_ovvo_MBPT2(mb,ej) - V(mn,ef)*(T2(bf,nj)+T1(f,j)*T1(b,n))
! G eq 34 4eme terme : t2_aemj V_mieb   !! t1*t1 disparait
!-------------------------------------------------------------------------------------
   !sum einstein  nf !      V    mieb                  aemj          iabj !!! Attention au indice a gauche/droite
   !sum einstein  nf !      Vmnef: mnef            t2: bfnj          mbej
  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo_MBPT2, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)
!--------------------------------------------------
!  w_ovvo(mb,ej) = w(mb,ej) + V(mb,ef)*t1(f,j) 
! G eq 34 2 eme terme : tej Vaieb=tvo Vvovv
!--------------------------------------------------
   !sum einstein    !  V iabe                   t1 : ej          iabj  
   !sum einstein f  !  Vmbef vovv               t1 : fj          mbej
  call contraction_424 ((/"o2","p1","p2","p3"/),(/"p3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo_MBPT2, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)
 end subroutine
!-----=========================--------------

 subroutine w_mbej(t1,t2,w_ovvo)
!      use symm_mod, only: IIOV, IOVVO, IVOOV, MULTB, NDIMT1, NDIMT2, NO, NREP, NV, NV2, NV4, NVO
!      use complex_mod, only: RCW
  external srt1l1
!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_ovvo(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real(8), allocatable :: w_vovo(:),array_temp(:)
!--------------------
! initialize w_voov
!--------------------
  allocate (w_vovo(nv4*rcw))
  call getvovo (w_vovo)
  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,w_ovvo) 
  deallocate (w_vovo)
!---------------------------------------------------
!   w_ovvo(mb,ej) = w_voov(mb,ej) - V(mn,ej)*T1(b,n)
!---------------------------------------------------
   call contraction_424 ((/"o2","o3","p2","o1"/),(/"p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
            &      -1.0d0,1.0d0,nrep,RightTensor=t1)
!-------------------------------------------------------------------------------------
!  w_ovvo(mb,ej) = w_ovvo(mb,ej) - V(mn,ef)*(T2(bf,nj)+T1(f,j)*T1(b,n))
!-------------------------------------------------------------------------------------
  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)
  allocate(array_temp(nv2*rcw))
  array_temp = 0.0d0
  call contraction_424 ((/"o2","o3","p2","p3"/),(/"p3","o1"/),(/"o2","o3","p2","o1"/),array_temp, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)
  call contraction_424 ((/"o2","o3","p2","o1"/),(/"p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=array_temp,RightTensor=t1)
  deallocate(array_temp)
!--------------------------------------------------
!  w_ovvo(mb,ej) = w(mb,ej) + V(mb,ef)*t(f,j) 
!--------------------------------------------------
  call contraction_424 ((/"o2","p1","p2","p3"/),(/"p3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)
 end subroutine
!-----=========================--------------


 subroutine w_iemn_MBPT2(t1,t2,fov,w_ovoo_MBPT2)

!---------------Description--------------------------------------------
!     Calculates w_voov intermediate Gwaltney eq.36
!-----------------------------------------------------------------------

!      use symm_mod, only: IOVOOT, NDIMT1, NDIMT2, NFVO, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)                                    ! Amplitudes.
   real(8), intent(out) :: w_ovoo_MBPT2(iovoot(nrep+1)*rcw)                                       ! Intermediate.
!   real(8), intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),w_oooo(rcw*nv1) !,hov(nfvo*rcw)       ! reuse of intermediates.
   real(8), intent(in)  :: fov(nfvo*rcw)   !!!!!!!! Attention verifier les bornes les memes que hov
!---------------Local Variables----------------------------------------
!--------------------
! initialize wbar_ovvo_MBPT2
!--------------------
  call getovoo (w_ovoo_MBPT2) 
!-----------------------------------------------------
! w_ovoo(ie,mn) = w_ovoo(ie,mn) + F_me(if) * T2(ef,mn)  
! w_ovoo(ia,jk) = w_ovoo(ia,jk) + f_me(ie) * T2(ae,jk)  
!-----------------------------------------------------
   ! sum einstein  e   ! fme ie       T2 eajk                  W; iajk  
   ! sum einstein  f   ! Fme if       T2 famn                  W: iamn
  call contraction_244 ((/"o1","p2"/),(/"p2","p1","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo_MBPT2, &
           &       1.0d0,1.0d0,nrep,LeftTensor=fov,RightTensor=t2) !! Attention il faut aller chercher fov


!-------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) +  V(io,mf) * T2(ef,no)
!   w_ovoo(ia,jk) = w_ovoo(ia,jk) +  V(im,je) * T2(ae,mk)   !!!!!!! Besoin d un P ij !!!!!!!
!-------------------------------------------------------
   ! sum einstein  mk  !    Voooe: imje         T2: aekm              Wovoo:  iajk
   ! sum einstein  of  !   Voooe: iomf          T2: efno               Wovoo : iemn       
  call contraction_444 ((/"o1","o4","o2","p2"/),(/"p1","p2","o3","o4"/),(/"o1","p1","o2","o3"/),w_ovoo_MBPT2, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)
!-------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) +  w_bar_mbej(ie,fn)*T1(fm)  !!!!!!!!!!!!!!!!!
!   w_ovoo(ia,jk) = w_ovoo(ia,jk) +  V(ia,ek)*t1(ej)
!--------------------------------------------
   !sum einstein   ! W iefn                       t1: fm            W: ienm  !!!!!!!!! INVERSION des INDICES 
!   call contraction_424 ((/"o1","p1","p2","o3"/),(/"p2","o2"/),(/"o1","p1","o3","o2"/),w_ovoo_MBPT2, &
!             &      1.0d0,-1.0d0,nrep,LeftTensor=wbar_ovvo,RightTensor=t1)
   !sum einstein e   !  V iaek                   t1 : ej          iabj       !!!!!! INVERSION des INDICES
  call contraction_424 ((/"o1","p1","p2","o3"/),(/"p2","o2"/),(/"o1","p1","o3","o2"/),w_ovoo_MBPT2, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)
!-----------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + w_ijmn(io,mn) * t(e,o)   !!!!!!!!!!!!!!!!!!!
!   w_ovoo(ia,jk) = w_ovoo(ia,jk) +  V(ia,ek)*t1(ej)   
!-----------------------------------------------------------
   ! sum einstein        ! Wijmn iomn              t1: eo         W: iemn
!   call contraction_424 ((/"o1","o4","o2","o3"/),(/"p1","o4"/),(/"o1","p1","o2","o3"/),w_ovoo_MBPT2, &
!             &      -1.0d0,1.0d0,nrep,LeftTensor=w_oooo,RightTensor=t1)
   !sum einstein occ  e ! V imjk                  t1 : am         W: aijk    !!!!!!!!!!!! INVERSION des INDICES
  call contraction_424 ((/"o1","o4","o2","o3"/),(/"p1","o4"/),(/"o1","p1","o2","o3"/),w_ovoo_MBPT2, &
            &      1.0d0,1.0d0,nrep,RightTensor=t1) !!!! Changement de signe 1er terme pour addition
!-----------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + V(ie,fg) * tau(fg,mn) 
!  w_ovoo(ia,jk) = w_ovoo(ia,jk) +  V(ia,ek)*t2(ef,jk) 
!-----------------------------------------------------------
!**out of core term..
!  allocate (tau(ndimt2*rcw))
!  call gettau (t1,t2,tau)
!    !sum einstein      ! V iefg                tau fgmn                    W iemn
!  call contraction_444 ((/"o1","p1","p2","p3"/),(/"p2","p3","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo_MBPT2, &
!             &      1.0d0,1.0d0,nrep,RightTensor=tau)
!  deallocate(tau) 
        !sum einstein   ! V iaef                   t2: efjk                 W: iajk       
   call contraction_444 ((/"o1","p1","p2","p3"/),(/"p2","p3","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo_MBPT2,&
  &1.0d0,1.0d0,nrep,RightTensor=t2)


  end subroutine


 subroutine w_iemn(t1,t2,w_ovoo,wbar_ovvo,w_oooo,hov)

!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------

!      use symm_mod, only: IOVOOT, IOVVO, NDIMT1, NDIMT2, NFVO, NREP, NV1
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)                                    ! Amplitudes.
   real(8), intent(out) :: w_ovoo(iovoot(nrep+1)*rcw)                                       ! Intermediate.
   real(8), intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),w_oooo(rcw*nv1),hov(nfvo*rcw)       ! reuse of intermediates.
!---------------Local Variables----------------------------------------
  real*8, allocatable   :: tau(:) 
!--------------------
! initialize wbar_ovvo
!--------------------
  call getovoo (w_ovoo) 
!-----------------------------------------------------
! w_ovoo(ie,mn) = w_ovoo(ie,mn) + F_me(if) * T2(ef,mn)  
!-----------------------------------------------------
  call contraction_244 ((/"o1","p2"/),(/"p2","p1","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
           &       1.0d0,1.0d0,nrep,LeftTensor=hov,RightTensor=t2)
!-------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) +  V(io,mf) * T2(ef,no)
!-------------------------------------------------------
  call contraction_444 ((/"o1","o4","o2","p2"/),(/"p1","p2","o3","o4"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)
!-------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) +  w_bar_mbej(ie,fn)*T1(fm)
!--------------------------------------------
   call contraction_424 ((/"o1","p1","p2","o3"/),(/"p2","o2"/),(/"o1","p1","o3","o2"/),w_ovoo, &
             &      1.0d0,-1.0d0,nrep,LeftTensor=wbar_ovvo,RightTensor=t1)
!-----------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + w_ijmn(io,mn) * t(e,o)
!-----------------------------------------------------------
   call contraction_424 ((/"o1","o4","o2","o3"/),(/"p1","o4"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      -1.0d0,1.0d0,nrep,LeftTensor=w_oooo,RightTensor=t1)
!-----------------------------------------------------------
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + V(ie,fg) * tau(fg,mn) 
!-----------------------------------------------------------
!**out of core term..
  allocate (tau(ndimt2*rcw))
  call gettau (t1,t2,tau)
  call contraction_444 ((/"o1","p1","p2","p3"/),(/"p2","p3","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      1.0d0,1.0d0,nrep,RightTensor=tau)
  deallocate(tau)
  end subroutine

 subroutine wbar_mbej(t2,wbar_ovvo)

!      use symm_mod, only: IIOV, IOVVO, MULTB, NDIMT2, NO, NREP, NV, NV4, NVO
!      use complex_mod, only: RCW
  external srt1l1
 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t2(rcw*ndimt2)                              ! Amplitudes.
   real(8), intent(out) :: wbar_ovvo(rcw*iovvo(nrep+1))                ! Intermediate.
 !---------------Local Variables----------------------------------------
  real*8, allocatable   :: w_vovo(:) 
 !--------------------
 ! initialize wbar_ovvo
 !--------------------

  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,wbar_ovvo) 

  deallocate (w_vovo)

!------------------------------------------------------------
!   wbar_ovvo(mb,ej) = wbar_ovvo(mb,ej) - V(mn,ef) * T2(bf,nj)
!------------------------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o1","o3"/),(/"o2","p1","p2","o1"/),wbar_ovvo, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)

 end subroutine

 subroutine wbar_mbej_MBPT2(t2,wbar_ovvo)  !!!!! Identique a eom_ee mais ne sert pas ????

!      use symm_mod, only: IIOV, IOVVO, MULTB, NDIMT2, NO, NREP, NV, NV4, NVO
!      use complex_mod, only: RCW
  external srt1l1
 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t2(rcw*ndimt2)                              ! Amplitudes.
   real(8), intent(out) :: wbar_ovvo(rcw*iovvo(nrep+1))                ! Intermediate.
 !---------------Local Variables----------------------------------------
  real*8, allocatable   :: w_vovo(:) 
 !--------------------
 ! initialize wbar_ovvo
 !--------------------

  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,wbar_ovvo) 

  deallocate (w_vovo)

!------------------------------------------------------------
!   wbar_ovvo(mb,ej) = wbar_ovvo(mb,ej) - V(mn,ef) * T2(bf,nj)
!------------------------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o1","o3"/),(/"o2","p1","p2","o1"/),wbar_ovvo, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)

 end subroutine




 subroutine w_mnie_MBPT2(t1,int_wooop_MBPT2)
!---------------Description--------------------------------------------
!     Calculates w_mnie intermediate
!-----------------------------------------------------------------------
!      use symm_mod, only: NDIMT1, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables----------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1)                                    ! Amplitudes.
  real(8), intent(out) :: int_wooop_MBPT2(:)                                      ! Intermediate.
!------------------------------------------------------------
! initialize int_wooop
!------------------------------------------------------------
  call getooov(int_wooop_MBPT2)
!----------------------------------------------------------------------------
!     w(mn,ie) = V(mn,ie) + t(f,i) * V(mn,fe)
!     wmnie_MBPT2(ijka) = V(jkia) + t1(e,i)*V(jk,ea) Gwaltney eq.32 
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Pourquoi pas d' INVERSION des INDICES cette fois-ci !!!!!!!!!!!!!!!!!!!!   
  call contraction_424 ((/"o1","o2","p3","p2"/),(/"p3","o3"/),(/"o1","o2","o3","p2"/), &
 &int_wooop_MBPT2,1.0d0,1.0d0,nrep,RightTensor=t1)
 end subroutine


 subroutine w_mnie(t1,int_wooop)
!---------------Description--------------------------------------------
!     Calculates w_mnie intermediate
!-----------------------------------------------------------------------
!      use symm_mod, only: NDIMT1, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables----------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1)                                    ! Amplitudes.
  real(8), intent(out) :: int_wooop(:)                                      ! Intermediate.
!------------------------------------------------------------
! initialize int_wooop
!------------------------------------------------------------
  call getooov(int_wooop)
!----------------------------------------------------------------------------
!     w(mn,ie) = V(mn,ie) + t(f,i) * V(mn,fe)
!----------------------------------------------------------------------------
  call contraction_424 ((/"o1","o2","p3","p2"/),(/"p3","o3"/), &
 &        (/"o1","o2","o3","p2"/),int_wooop,1.0d0,1.0d0,nrep,RightTensor=t1)
 end subroutine

subroutine w_amef_MBPT2(t1,w_v_temp_MBPT2)
!---------------Description--------------------------------------------
!     Calculates w_amef intermediate
!     Gwaltney eq.31
!-----------------------------------------------------------------------
!      use symm_mod, only: NDIMT1, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables----------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1)                                    ! Amplitudes.
  real(8), intent(out) :: w_v_temp_MBPT2(:)                                      ! Intermediate.
!--------------------------------------------------
! initialize w_v_temp_MBPT2
! <ai||bc>
!--------------------------------------------------
    call getvovv_incore(w_v_temp_MBPT2)
!---------------------------------------------------
! w_amef(am,ef) = w_amef(am,ef) - V(nm,ef) * t(a,n) 
! w_aibc = w_aibc - t1(a,m)*V(mi,bc)
!---------------------------------------------------
    ! sum einstein  m    ! t1 am        V: mibc            W: aibc
    !sum einstein   n    ! t1: an       V: nmef            W: amef
    call contraction_244 ((/"p2","o2"/),(/"o2","o1","p3","p1"/),(/"p2","o1","p3","p1"/),&
   & w_v_temp_MBPT2,-1.0d0,1.0d0,nrep,LeftTensor=t1)
 end subroutine


 subroutine w_amef(t1,w_v_temp)
!---------------Description--------------------------------------------
!     Calculates w_amef intermediate
!-----------------------------------------------------------------------
!      use symm_mod, only: NDIMT1, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables----------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1)                                    ! Amplitudes.
  real(8), intent(out) :: w_v_temp(:)                                      ! Intermediate.
!--------------------------------------------------
! initialize w_v_temp
!--------------------------------------------------
    call getvovv_incore(w_v_temp)
!---------------------------------------------------
! w_amef(am,ef) = w_amef(am,ef) - V(nm,ef) * t(a,n) 
!---------------------------------------------------
    call contraction_244 ((/"p2","o2"/),(/"o2","o1","p3","p1"/), &
   &    (/"p2","o1","p3","p1"/),w_v_temp,-1.0d0,1.0d0,nrep,LeftTensor=t1)
 end subroutine

!!!!!!!!!!!!!!!!!!!!1 Verifier que w_efam existe dans sigma_eom.F90  Il y en a 2... ;-)
subroutine w_efam_MBPT2(t1,t2,hov,w_vvvo_MBPT2)
!--------------------DESCRIPTION---------------------------
!        calculates vvvo type intermediate in-core fashion.      
! Gwaltney eq.35 ordre des termes 1-4-3-2-5-6
!----------------------------------------------------------
!----------calling variables----------------------
!      use symm_mod, only: NDIMT1, NDIMT2, NFVO, NREP, NV5
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
   real*8, intent(in)  :: t2(rcw*ndimt2),t1(rcw*ndimt1)         ! Amplitudes.
!   real*8, intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),hov(nfvo*rcw)     ! Intermediate.
   real*8, intent(in)  :: hov(nfvo*rcw)
   real*8, intent(out) :: w_vvvo_MBPT2(rcw*nv5)                ! Intermediate.
!   real(8), intent(in)  :: fov(nfov*rcw)   !!!!!!!! Attention verifier les bornes les memes que hvv
!---------------Local Variables----------------------------------------
   integer :: AllocateStatus 
   real*8, allocatable  :: tau(:)                
   real*8, allocatable  :: vovv1(:),temp(:)                
   type(indices) :: e,f,g
!----------------------------------------------------------
!        Initialize
!    W abci V(ab,ci) 
!----------------------------------------------------------
   call getvvvo_incore(w_vvvo_MBPT2)
!-----------------------------------------------------------
!   diagram1 : w_vvvo(ef,am)  = w_vvvo(ef,am) + p(ef) V(en,ag)*t2(fg,mn)  

!   Wefam Wefam - p(ab)Vcimn t2(abmn)
!   abci
!   Je cherche une contraction entre Vamce et t2_ebim : Wabci
!-----------------------------------------------------------
   ! sum einstein       ! V_vovv ajdc              t2_vvoo: bdij       W_vvvo: abci
   ! sum einstein   ng  ! V_vovv enga              t2_vvoo: fgmn       W_vvvo: efam
  call contraction_444 ((/"p2","o2","p4","p1"/),(/"p3","p4","o1","o2"/),(/"p2","p3","p1","o1"/),w_vvvo_MBPT2, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)
!------------------------------------------------------------
!  diagram2 : w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m) 
! Gwaltney  : w_vvvo(abci)  = w_vvvo(abci) + V(abce)*t(e,i)
!------------------------------------------------------------

! Loic 2509 : ici je fais appel à w_efab special qui ne contient que get vvvv
! Je ne peux pas utiliser le w_efab car il sert pour H_DD je ne peux donc pas le modifier
! Comme je ne sais pas pour l'instant programer un get vvvv out of core seul je crée des sub routines
! qui seront certainement inutiles par la suite


!-------------------------------------------------------------
!  diagram3 : w_vvvo(ef,am) =  w_vvvo(ef,am) + f_bar_me(na) * t(ef,nm)
! Rq pourquoi ecrire hov et parler de Fba_me dans l'equation ?????????????????????????????????
!  sum einstein      ! t2vvoo:efnm          !fbar_me(ov):na           w_vvvo : efamm  !!!!!!! INVERSION des INDICES
!      call contraction_424 ((/"p2","p3","o2","o1"/),(/"o2","p1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
!     &      -1.0d0,1.0d0,nrep,RightTensor=fov,LeftTensor=t2)

!  \td abim \fme mc
!-------------------------------------------------------------
!  sum einstein      ! t2vvoo:abmi          !f_me(ov): mc           w_vvvo : abci
      call contraction_424 ((/"p2","p3","o2","o1"/),(/"o2","p1"/),(/"p2","p3","p1","o1"/),w_vvvo_MBPT2, &
     &      -1.0d0,1.0d0,nrep,RightTensor=hov,LeftTensor=t2)

!-----------------------------------------------------------
!   diagram4 : w_vvvo(ef,am) =  w_vvvo(ef,am) + V(no,am) * tau(ef,no) 
!   Gwaltney ce n est pas tau mais t2 sinon ce sont bien les 2 indices occuppes de t2 qui partent   
!-----------------------------------------------------------
!   sum einstein  !  t2_vvoo : abmn           V mnci               W abci
! sum einstein    ! (tau)vvoo (ef,no)           V_oovo noam                W amef
   call contraction_444 ((/"p2","p3","o2","o3"/),(/"o2","o3","p1","o1"/),(/"p2","p3","p1","o1"/),w_vvvo_MBPT2, &
  &      1.0d0,1.0d0,nrep,LeftTensor=t2)


!---------------------------------------------------------------------------
!  diagram5 : w_vvvo(ef,am) =  w_vvvo(ef,am) - w_bar_ovvo(nf,am) * t(e,n)
! Gwalltney 6eme terme :  -P(ab) t(a,m) V(mbci)  
!---------------------------------------------------------------------------
! Je prends amef du dessus ici sommes sur l occupee n
!sum einstein   n    ! t1: an       V: nmef            W: amef
!    call contraction_244 ((/"p2","o2"/),(/"o2","o1","p3","p1"/),(/"p2","o1","p3","p1"/),&
!   & w_v_temp_MBPT2,-1.0d0,1.0d0,nrep,LeftTensor=t1)
  ! sum einstein sur o2  ! t1 am          V mbci       W: abci
    call contraction_244 ((/"p1","o2"/),(/"o2","p2","p3","o1"/),(/"p1","p2","p3","o1"/),w_vvvo_MBPT2,&
   &-1.0d0,1.0d0,nrep,LeftTensor=t1)


 end subroutine


 subroutine w_efam(t1,t2,wbar_ovvo,hov,w_vvvo)
!--------------------DESCRIPTION---------------------------
!        calculates vvvo type intermediate in-core fashion.      
!----------------------------------------------------------
!----------calling variables----------------------
!      use symm_mod, only: IOVVO, NDIMT1, NDIMT2, NFVO, NREP, NV5
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
   real*8, intent(in)  :: t2(rcw*ndimt2),t1(rcw*ndimt1)         ! Amplitudes.
   real*8, intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),hov(nfvo*rcw)     ! Intermediate.
   real*8, intent(out) :: w_vvvo(rcw*nv5)                ! Intermediate.
!---------------Local Variables----------------------------------------
   integer :: AllocateStatus 
   real*8, allocatable  :: tau(:)                
   real*8, allocatable  :: vovv1(:),temp(:)                
   type(indices) :: e,f,g
!----------------------------------------------------------
!        Initialize
!----------------------------------------------------------
   call getvvvo_incore(w_vvvo)
!-----------------------------------------------------------
!   diagram1 : w_vvvo(ef,am)  = w_vvvo(ef,am) + p(ef) V(en,ag)*t2(fg,mn)  
!-----------------------------------------------------------
  call contraction_444 ((/"p2","o2","p4","p1"/),(/"p3","p4","o1","o2"/),(/"p2","p3","p1","o1"/),w_vvvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)
!------------------------------------------------------------
!  diagram2 : w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m) 
!------------------------------------------------------------


!-------------------------------------------------------------
!  diagram3 : w_vvvo(ef,am) =  w_vvvo(ef,am) + f_bar_me(na) * t(ef,nm)
!-------------------------------------------------------------
      call contraction_424 ((/"p2","p3","o2","o1"/),(/"o2","p1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
     &      -1.0d0,1.0d0,nrep,RightTensor=hov,LeftTensor=t2)
!-----------------------------------------------------------
!   diagram4 : w_vvvo(ef,am) =  w_vvvo(ef,am) + V(no,am) * tau(ef,no) 
!-----------------------------------------------------------
  allocate (tau(ndimt2*rcw))
  call gettau (t1,t2,tau)
   call contraction_444 ((/"p2","p3","o2","o3"/),(/"o2","o3","p1","o1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
  &      1.0d0,1.0d0,nrep,LeftTensor=tau)
  deallocate (tau)
!---------------------------------------------------------------------------
!  diagram5 : w_vvvo(ef,am) =  w_vvvo(ef,am) - w_bar_ovvo(nf,am) * t(e,n)
!---------------------------------------------------------------------------
  call contraction_244 ((/"p2","o2"/),(/"o2","p3","p1","o1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
            &      -1.0d0,1.0d0,nrep,RightTensor=wbar_ovvo,LeftTensor=t1)
 end subroutine

 subroutine w_efam_diagram2_MBPT2(t1,t2,buf2,buf3,nbuf3,w_vvvo_MBPT2)

! Loic 2509 : cette routine fait appel à une fonction w_efab speciale qui ne contient
! que get_vvvv. En effet pour Gwaltney il ne me faut pes W_efab modifié mais uniquement
! \V vvvv

!      use complex_mod, only: CARITH, RCW
!      use symm_mod, only: IIVO, IVO, IVOVVT, MULTB, NDIMT1, NDIMT2, NO, NREP, NV, NV5, NVO, NVVT
#include "complex.inc"
#include "symm.inc"   
!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vvvo_MBPT2(rcw*nv5)                ! Intermediate.
   real(8), intent(in)    :: buf2(*)                        ! buffer arrays
   real(8), intent(inout) :: buf3(*)                        ! buffer arrays
   integer, intent(in)    :: nbuf3 

!---------------Local variables--------------------------------------
                                                                    
      real*8  :: t1val(2)                                                
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ci,akci,t1di,di,akdi,t1ci,cmin,c,d,i,j
      integer :: crep,drep,irp,mint
      real*8,allocatable  :: tau (:)

!------------------------------------------------------------------
!  w_vvvo(ef,am) =  w_vvvo(ef,am) + \V efag \ts gm à la place de :  W (ef,ag)*t(g,m)
!------------------------------------------------------------------

      t1val = (/0.0d0,0.0d0/)

     allocate (tau(ndimt2*rcw))
     call gettau (t1,t2,tau)

      do 30 irp = 1, nrep
         if ((nvvt(irp).eq.0).or.(nvo(irp).eq.0)) goto 30
         done = .false.
         istart = 0
         cco=0

! ** w_efab is constructed over the batches.
! ** ISTART is updated with the actual batch #
! ** available on the local node, NINT is updated with the number of
! ** actually read integrals from this batch and is a pure output parameter !
! ** We therefore have to count upwards in case of ISTART.NE.0 !!

! Dans la ligne ci-dessous ; tau est remplacée par t2 
         call w_efab_oocore_wefam_d2_MBPT2 (irp,istart,done,t2,buf2,t1,nbuf3,buf3,nint)

         call putdz (irp,istart,nint,buf3)

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
! dans la ligne ci-dessous tau est remplacé par t2
                    call w_efab_oocore_wefam_d2_MBPT2 (irp,istart,done,t2,buf2,t1,nbuf3,buf3,nint)
                     
                    call putdz (irp,istart,nint,buf3)
                      CD = 1
                      CCO=ISTART
                    ELSE
                      GOTO 30
                    ENDIF
                 ENDIF
                 AKCD = (CD-1)*NVVT(IRP)*RCW+1
                 DO I = 1, NO(DREP)
                    CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    AKCI = (IVOVVT(IRP)+(CI-1)*NVVT(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) =  T1(T1DI)
                    IF (CARITH) T1VAL(2) =  T1(T1DI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,w_vvvo_MBPT2(AKCI),1)
                 ENDDO
                 DO I = 1, NO(CREP)
                    DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    AKDI = (IVOVVT(IRP)+(DI-1)*NVVT(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = -T1(T1CI)
                    IF (CARITH) T1VAL(2) = -T1(T1CI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,w_vvvo_MBPT2(AKDI),1)
                 ENDDO
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE

   deallocate (tau)

 end subroutine


 subroutine w_efam_diagram2(t1,t2,buf2,buf3,nbuf3,w_vvvo)

!      use complex_mod, only: CARITH, RCW
!      use symm_mod, only: IIVO, IVO, IVOVVT, MULTB, NDIMT1, NDIMT2, NO, NREP, NV, NV5, NVO, NVVT
#include "complex.inc"
#include "symm.inc"   
!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vvvo(rcw*nv5)                ! Intermediate.
   real(8), intent(in)    :: buf2(*)                        ! buffer arrays
   real(8), intent(inout) :: buf3(*)                        ! buffer arrays
   integer, intent(in)    :: nbuf3 

!---------------Local variables--------------------------------------
                                                                    
      real*8  :: t1val(2)                                                
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ci,akci,t1di,di,akdi,t1ci,cmin,c,d,i,j
      integer :: crep,drep,irp,mint
      real*8,allocatable  :: tau (:)

!------------------------------------------------------------------
!  w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m)
!------------------------------------------------------------------

      t1val = (/0.0d0,0.0d0/)

     allocate (tau(ndimt2*rcw))
     call gettau (t1,t2,tau)

      do 30 irp = 1, nrep
         if ((nvvt(irp).eq.0).or.(nvo(irp).eq.0)) goto 30
         done = .false.
         istart = 0
         cco=0

! ** w_efab is constructed over the batches.
! ** ISTART is updated with the actual batch #
! ** available on the local node, NINT is updated with the number of
! ** actually read integrals from this batch and is a pure output parameter !
! ** We therefore have to count upwards in case of ISTART.NE.0 !!

         call w_efab_out_of_core (irp,istart,done,tau,buf2,t1,nbuf3,buf3,nint)

         call putdz (irp,istart,nint,buf3)

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

                    call w_efab_out_of_core (irp,istart,done,tau,buf2,t1,nbuf3,buf3,nint)
                     
                    call putdz (irp,istart,nint,buf3)
                      CD = 1
                      CCO=ISTART
                    ELSE
                      GOTO 30
                    ENDIF
                 ENDIF
                 AKCD = (CD-1)*NVVT(IRP)*RCW+1
                 DO I = 1, NO(DREP)
                    CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    AKCI = (IVOVVT(IRP)+(CI-1)*NVVT(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) =  T1(T1DI)
                    IF (CARITH) T1VAL(2) =  T1(T1DI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,w_vvvo(AKCI),1)
                 ENDDO
                 DO I = 1, NO(CREP)
                    DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    AKDI = (IVOVVT(IRP)+(DI-1)*NVVT(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = -T1(T1CI)
                    IF (CARITH) T1VAL(2) = -T1(T1CI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,w_vvvo(AKDI),1)
                 ENDDO
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE

   deallocate (tau)

 end subroutine

 subroutine w_ejmb_out_of_core(t1,t2,buf1,buf2,buf3,nbuf3,w_vovo)

!      use param_mod, only: A0, D0
!      use complex_mod, only: CARITH, RCW
!      use symm_mod, only: IIOV, IVO, IVOOV, MULTB, NDIMT1, NDIMT2, NO, NREP, NV, NV4, NVO
!      use ccpar_mod, only: MASTER, MYPROC
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vovo(rcw,nv4)                ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*),buf3(*)        ! buffer arrays
   integer, intent(in)    :: nbuf3 

!---------------Local variables--------------------------------------
!                                                                    
      real*8  :: T1VAL(2)                                                
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ic,akic,t1di,id,akid,t1ci,cmin,c,d,i
      integer :: crep,drep,irp,mint
!------------------------------------------------------------------
! H(AK,CI) = H(AK,CI) - W(AK,CD) * T(D,I)
! H is kept in BUF1, still ordered VOVO
!------------------------------------------------------------------

      t1val = (/0.0d0,0.0d0/)

      if(myproc.eq.master) then
        CALL xcopy(nv4,w_vovo,1,BUF1,1)
      else 
        call xcopy(nv4,a0,0,buf1,1)
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
!                   CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    IC = IIOV(CREP,DREP)+(C-1)*NO(DREP)+I
!                   AKCI = (IVOVO(IRP)+(CI-1)*NVO(IRP))*RCW+1
                    AKIC = (IVOOV(IRP)+(IC-1)*NVO(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) =  T1(T1DI)
                    IF (CARITH) T1VAL(2) =  T1(T1DI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKIC),1)
                 ENDDO
                 DO I = 1, NO(CREP)
!                   DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    ID = IIOV(DREP,CREP)+(D-1)*NO(CREP)+I
!                   AKDI = (IVOVO(IRP)+(DI-1)*NVO(IRP))*RCW+1
                    AKID = (IVOOV(IRP)+(ID-1)*NVO(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = -T1(T1CI)
                    IF (CARITH) T1VAL(2) = -T1(T1CI+1)
                    CALL XAXPY(NVO(IRP),T1VAL,BUF3(AKCD),1,BUF1(AKID),1)
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
!   CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO, &
!   &            BUF1,w_vovo)

!  call srt1r1 (nrep,multb,lfa,nvo,nv,no,nvo,ivovo,iiov,buf1,buf2)
!  
!  call xaxpy (nv4,-A1,buf2,1,w_vovo,1) 

   CALL XTIME(0,2,'--- HINTM: VOVV*T              ')

 end subroutine

 subroutine f_bar_ae_out_of_core (t1,buf2,nbuf2,hv)

!      use param_mod, only: A0
!      use complex_mod, only: RCW
!      use symm_mod, only: JJVO, JJVV, JVVVO, MULTB, MVO, MVV, NFVO, NFVV, NO, NREP, NV, NVO
!      use ccpar_mod, only: MASTER, MYPROC
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ccpar.inc"

!---------------Calling variables--------------------------------------
     real*8, intent(in)    :: t1(rcw,nfvo),buf2(*)
     real*8, intent(inout) :: hv(rcw*nfvv)
     integer,intent(in)    :: nbuf2 
!---------------Common Blocks--------------------------------------

!---------------Local variables----
     logical :: usedz,right
     integer :: m,n
!----------------------------------
      if(myproc==master) then
        call xcopy(nfvv*rcw,A0,0,buf2,1)
      else 
        call xcopy(nfvv*rcw,A0,0,buf2,1)
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

!---------------Common Blocks--------------------------------------   
                                                                     
!      use param_mod, only: A0, A1, D0, LFA
!      use complex_mod, only: RCW
!      use symm_mod, only: IIOV, IOVOOT, MULTB, NDIMT1, NDIMT2, NO, NOOT, NREP, NV, NVO, NVVT
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   

!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)! Amplitudes.
   real(8), intent(inout) :: w_ovoo(rcw,iovoot(nrep+1))   ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*)              ! buffer arrays.
   integer, intent(in)    :: nbuf3 

!---------------Local variables----------------
      logical :: done,teq
      integer :: kd,istart,nint,m,n,k
      integer :: irp,mint,off1,off2
      real*8,allocatable :: tau(:)
!----------------------------------------------

!----------------------------------------------
!w(ie,mn) = w(ie,mn) + 1/2 V(ie,fg)*tau(fg,mn)
!----------------------------------------------
   call xcopy(iovoot(nrep+1),A0,0,buf2,1)

   allocate(tau(ndimt2*rcw))
   tau = 0.0d0
   call gettau (t1,t2,tau)

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

           CALL XGEMM ('N','N',M,N,K,A1,BUF1,M,tau(off1+istart*rcw), &
            &               KD,A1,BUF2(OFF2),M)
           IF (.NOT.DONE) THEN
           ISTART = ISTART + NINT
           GOTO 110
           ENDIF
           OFF1 = OFF1 + KD * N * RCW
           OFF2 = OFF2 + M * N * RCW
 130  CONTINUE

    deallocate (tau)

      call srt1l1 (nrep,multb,lfa,nvo,nv,no,noot,iovoot,iiov,buf2,buf1)

      call xaxpy (iovoot(nrep+1),A1,buf1,1,w_ovoo,1) 

  end subroutine


   subroutine w_efab_out_of_core_MBPT2 (irp,istart,done,t2,buf1,t1,nbuf3,w_vvvv_MBPT2,nint)

! Loic 2509 : Je crée cette subroutine cqr Gwaltney en a besoin dans le calcul de H_DD


!-----------------------description------------------------------------------------
!   W_vvvv intermediates are formed such a manner that it could be used as fetching 
!   VVVV integrals, representationwise. However, we will call them only once to construct
!   term2 of W_efam integrals and transfer it to a file, subsequently, inside that routine.
!   batch size is fixed to the distributed batches of VVVV integrals. I assumed that the batch 
!   size for VOVV and VVVV integrals are the same.
!-----------------------------------------------------------------------------------

!      use param_mod, only: A0, A1, D0
!      use complex_mod, only: RCW
!      use symm_mod, only: IVVOOTT, MULTB, NDIMT1, NDIMT2, NO, NOOT, NREP, NT3VVT, NV, NVO, NVVT
!      use ccpar_mod, only: MYPROC, NMPROC
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!--------------------calling variables-------------
   real(8), intent(in)    :: t2(ndimt2*rcw),t1(ndimt1*rcw)  ! Amplitudes.
   real(8), intent(out)   :: w_vvvv_MBPT2(*)              ! Intermediate.
   real(8), intent(in)    :: buf1(*)                ! buffer arrays
   integer, intent(in)    :: nbuf3,irp,istart
   integer, intent(out)   :: nint
   logical, intent(inout) :: done
!--------------------------------------------------

!--------------------local variables----------------
   integer :: mint,off1,off2,off3,off4,off5
   integer :: abmn,m,n,k,i,irrep,irep,mm(nrep)
   integer :: xxt3vvvt(nrep,nrep),yt3vvvt(nrep),zt3vvvot(nrep+1),zt3vvvvt(nrep+1)
   integer, allocatable :: counter(:), displ(:)
   real*8,allocatable :: w_vvvv_temp_temp(:),w_vvvv_temp(:),w_v(:)  
!---------------------------------------------------

         off1 = 1
         off2 = 1
         off3 = 0

       allocate(counter(0:nmproc-1))
       allocate(displ(0:nmproc-1))

       counter = 0
       displ = 0 

         mint = nbuf3/nvvt(irp)
         call getvvvv (irp,istart,nint,done,w_vvvv_MBPT2,mint)

!-----------------------------------------------------------------------
! generate offsets to distribute the intermediate on different nodes
!-----------------------------------------------------------------------
#if defined (VAR_MPI)
       call interface_mpi_allgather(NINT,1, &
     &   counter,1,global_communicator)

      do i = 0, (nmproc-1)

       counter(i) = counter(i)*noot(irp)            

      enddo

         displ(0) = 0
         do i = 1, (nmproc-1) 
             displ(i) = displ(i-1) + counter(i-1)
         enddo
#endif

!----------------------------------------------
! w(ef,ab) =  V(ef,ab) + V(mn,ab) * TAU(ef,mn)
!          =  V(ef,ab) + TAU(ef,mn) * V'*(ab,mn)        
!
!Remplacer par : 
! w_MBPT2(ef,ab) =  V(ef,ab) + V(mn,ab) * T2(ef,mn)
!          =  V(ef,ab) + T2(ef,mn) * V'*(ab,mn)        
!    
!----------------------------------------------
       call getoovv (buf1)
       m = nvvt(irp)
       n = nint
       k = noot(irp)
       abmn = (ivvoott(irp)+displ(myproc))*rcw+1
       off1 = off1 + ivvoott(irp)*rcw  

     if (noot(irp) == 0) then
        call xscal (m*n,a1,w_vvvv_MBPT2,1)
     else 
        call xgemm ('n','n',m,n,k,a1,t2(off1),nvvt(irp),buf1(abmn),k, &
      &               a1,w_vvvv_MBPT2,m)
     endif

!-----------------------------------------------------
! w(ef,ab) = w(ef,ab) - V(em,ab)*t(f,m)
!-----------------------------------------------------

      call offset_outofcore_sorter(irp,xxt3vvvt,yt3vvvt,zt3vvvot,zt3vvvvt)

      if (.not.done) off3 = off3+nint

       allocate(w_v(zt3vvvot(nrep+1)*rcw))  

       w_v = 0.0d0

      call srt33 (irp,istart,nrep,multb,nvo,nv,no,nt3vvt, & 
     &            yt3vvvt,zt3vvvot,xxt3vvvt,w_v,nbuf3,off3,done)

       mm(1:nrep) = 0
       do irrep = 1,nrep
         irep = multb(irrep,irp+nrep,2)
         mm(irrep) = mm(irrep)+ nint*nv(irep)
       enddo

       allocate(w_vvvv_temp_temp(dot_product(mm(1:nrep),nv(1:nrep))*rcw))  

       w_vvvv_temp_temp = 0.0d0

       call cntrct ('N','T',mm,nv,no,-a1,w_v,t1,a0,w_vvvv_temp_temp,nrep)

       deallocate(w_v)

       allocate(w_vvvv_temp(nvvt(irp)*nint*rcw)) 

      w_vvvv_temp = 0.0d0 
      call srt43 (irp,istart,nrep,multb,.true.,nvvt,nv,nv,  &
     &    nt3vvt,yt3vvvt,zt3vvvvt,xxt3vvvt,w_vvvv_temp_temp,&
     &    w_vvvv_temp,nbuf3,off3)

       deallocate(w_vvvv_temp_temp)

       call xaxpy(nint*nvvt(irp),a1,w_vvvv_temp,1,w_vvvv_MBPT2,1)

       deallocate(w_vvvv_temp)

   end subroutine

   subroutine w_efab_oocore_wefam_d2_MBPT2 (irp,istart,done,t2,buf1,t1,nbuf3,w_vvvv_d2_MBPT2,nint)

! Loic 2509 : Je crée cette subroutine qui ne contient que get vvvv pour W_efam diagram2 MBPT2


!-----------------------description------------------------------------------------
!   W_vvvv intermediates are formed such a manner that it could be used as fetching 
!   VVVV integrals, representationwise. However, we will call them only once to construct
!   term2 of W_efam integrals and transfer it to a file, subsequently, inside that routine.
!   batch size is fixed to the distributed batches of VVVV integrals. I assumed that the batch 
!   size for VOVV and VVVV integrals are the same.
!-----------------------------------------------------------------------------------

!      use complex_mod, only: RCW
!      use symm_mod, only: NDIMT1, NDIMT2, NOOT, NREP, NVVT
!      use ccpar_mod, only: NMPROC
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!--------------------calling variables-------------
   real(8), intent(in)    :: t2(ndimt2*rcw),t1(ndimt1*rcw)  ! Amplitudes.
   real(8), intent(out)   :: w_vvvv_d2_MBPT2(*)              ! Intermediate.
   real(8), intent(in)    :: buf1(*)                ! buffer arrays
   integer, intent(in)    :: nbuf3,irp,istart
   integer, intent(out)   :: nint
   logical, intent(inout) :: done
!--------------------------------------------------

!--------------------local variables----------------
   integer :: mint,off1,off2,off3,off4,off5
   integer :: abmn,m,n,k,i,irrep,irep,mm(nrep)
   integer :: xxt3vvvt(nrep,nrep),yt3vvvt(nrep),zt3vvvot(nrep+1),zt3vvvvt(nrep+1)
   integer, allocatable :: counter(:), displ(:)
   real*8,allocatable :: w_vvvv_temp_temp(:),w_vvvv_temp(:),w_v(:)  
!---------------------------------------------------

         off1 = 1
         off2 = 1
         off3 = 0

       allocate(counter(0:nmproc-1))
       allocate(displ(0:nmproc-1))

       counter = 0
       displ = 0 

         mint = nbuf3/nvvt(irp)
         call getvvvv (irp,istart,nint,done,w_vvvv_d2_MBPT2,mint)

!-----------------------------------------------------------------------
! generate offsets to distribute the intermediate on different nodes
!-----------------------------------------------------------------------
#if defined (VAR_MPI)
       call interface_mpi_allgather(NINT,1, &
     &   counter,1,global_communicator)

      do i = 0, (nmproc-1)

       counter(i) = counter(i)*noot(irp)            

      enddo

         displ(0) = 0
         do i = 1, (nmproc-1) 
             displ(i) = displ(i-1) + counter(i-1)
         enddo
#endif
end subroutine


   subroutine w_efab_out_of_core (irp,istart,done,tau,buf1,t1,nbuf3,w_vvvv,nint)

!-----------------------description------------------------------------------------
!   W_vvvv intermediates are formed such a manner that it could be used as fetching 
!   VVVV integrals, representationwise. However, we will call them only once to construct
!   term2 of W_efam integrals and transfer it to a file, subsequently, inside that routine.
!   batch size is fixed to the distributed batches of VVVV integrals. I assumed that the batch 
!   size for VOVV and VVVV integrals are the same.
!-----------------------------------------------------------------------------------

!      use param_mod, only: A0, A1, D0
!      use complex_mod, only: RCW
!      use symm_mod, only: IVVOOTT, MULTB, NDIMT1, NDIMT2, NO, NOOT, NREP, NT3VVT, NV, NVO, NVVT
!      use ccpar_mod, only: MYPROC, NMPROC
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!--------------------calling variables-------------
   real(8), intent(in)    :: tau(ndimt2*rcw),t1(ndimt1*rcw)  ! Amplitudes.
   real(8), intent(out)   :: w_vvvv(*)              ! Intermediate.
   real(8), intent(in)    :: buf1(*)                ! buffer arrays
   integer, intent(in)    :: nbuf3,irp,istart
   integer, intent(out)   :: nint
   logical, intent(inout) :: done
!--------------------------------------------------

!--------------------local variables----------------
   integer :: mint,off1,off2,off3,off4,off5
   integer :: abmn,m,n,k,i,irrep,irep,mm(nrep)
   integer :: xxt3vvvt(nrep,nrep),yt3vvvt(nrep),zt3vvvot(nrep+1),zt3vvvvt(nrep+1)
   integer, allocatable :: counter(:), displ(:)
   real*8,allocatable :: w_vvvv_temp_temp(:),w_vvvv_temp(:),w_v(:)  
!---------------------------------------------------

         off1 = 1
         off2 = 1
         off3 = 0

       allocate(counter(0:nmproc-1))
       allocate(displ(0:nmproc-1))

       counter = 0
       displ = 0 

         mint = nbuf3/nvvt(irp)
         call getvvvv (irp,istart,nint,done,w_vvvv,mint)

!-----------------------------------------------------------------------
! generate offsets to distribute the intermediate on different nodes
!-----------------------------------------------------------------------
#if defined (VAR_MPI)
       call interface_mpi_allgather(NINT,1, &
     &   counter,1,global_communicator)

      do i = 0, (nmproc-1)

       counter(i) = counter(i)*noot(irp)            

      enddo

         displ(0) = 0
         do i = 1, (nmproc-1) 
             displ(i) = displ(i-1) + counter(i-1)
         enddo
#endif

!----------------------------------------------
! w(ef,ab) =  V(ef,ab) + V(mn,ab) * TAU(ef,mn)
!          =  V(ef,ab) + TAU(ef,mn) * V'*(ab,mn)            
!----------------------------------------------
       call getoovv (buf1)
       m = nvvt(irp)
       n = nint
       k = noot(irp)
       abmn = (ivvoott(irp)+displ(myproc))*rcw+1
       off1 = off1 + ivvoott(irp)*rcw  

     if (noot(irp) == 0) then
        call xscal (m*n,a1,w_vvvv,1)
     else 
        call xgemm ('n','n',m,n,k,a1,tau(off1),nvvt(irp),buf1(abmn),k, &
      &               a1,w_vvvv,m)
     endif

!-----------------------------------------------------
! w(ef,ab) = w(ef,ab) - V(em,ab)*t(f,m)
!-----------------------------------------------------

      call offset_outofcore_sorter(irp,xxt3vvvt,yt3vvvt,zt3vvvot,zt3vvvvt)

      if (.not.done) off3 = off3+nint

       allocate(w_v(zt3vvvot(nrep+1)*rcw))  

       w_v = 0.0d0

      call srt33 (irp,istart,nrep,multb,nvo,nv,no,nt3vvt, & 
     &            yt3vvvt,zt3vvvot,xxt3vvvt,w_v,nbuf3,off3,done)

       mm(1:nrep) = 0
       do irrep = 1,nrep
         irep = multb(irrep,irp+nrep,2)
         mm(irrep) = mm(irrep)+ nint*nv(irep)
       enddo

       allocate(w_vvvv_temp_temp(dot_product(mm(1:nrep),nv(1:nrep))*rcw))  

       w_vvvv_temp_temp = 0.0d0

       call cntrct ('N','T',mm,nv,no,-a1,w_v,t1,a0,w_vvvv_temp_temp,nrep)

       deallocate(w_v)

       allocate(w_vvvv_temp(nvvt(irp)*nint*rcw)) 

      w_vvvv_temp = 0.0d0 
      call srt43 (irp,istart,nrep,multb,.true.,nvvt,nv,nv,  &
     &    nt3vvt,yt3vvvt,zt3vvvvt,xxt3vvvt,w_vvvv_temp_temp,&
     &    w_vvvv_temp,nbuf3,off3)

       deallocate(w_vvvv_temp_temp)

       call xaxpy(nint*nvvt(irp),a1,w_vvvv_temp,1,w_vvvv,1)

       deallocate(w_vvvv_temp)

   end subroutine

  subroutine G_intm (t2,l2,G_vv,G_oo)

!---------------Common Blocks--------------------------------------
!      use symm_mod, only: NDIMT2, NREP
!      use complex_mod, only: RCW
#include "symm.inc"
#include "complex.inc"
!     Calling variables
!     -----------------  
      real(8), intent(in)  :: l2(rcw*ndimt2),t2(rcw*ndimt2)  ! Amplitudes.
      real(8), intent(inout) :: G_vv(:),G_oo(:)

!----------initialization--------------
! G_vv = 0.0d0
! G_oo = 0.0d0
!--------------------------------------

!----------------------------------------------
!    G_ea(e,a)  = - l2(mn,af) * t2(ef,mn)
!----------------------------------------------

  call contraction_442 ((/"p2","p3","o1","o2"/),(/"o1","o2","p1","p3"/),(/"p2","p1"/),G_vv,-1.0d0,1.0d0,nrep,&
                        &  LeftTensor=t2,RightTensor=l2 )

!---------------------------------------------
!     G_im(i,m) = l2(in,ef) * t2(ef,mn)
!---------------------------------------------

  call contraction_442 ((/"o1","o3","p1","p2"/),(/"p1","p2","o2","o3"/),(/"o1","o2"/),G_oo,1.0d0,1.0d0,nrep,&
                        & LeftTensor=l2,RightTensor=t2)

  end subroutine


      subroutine offset_outofcore_sorter(irp,xxt3vvvt,yt3vvvt,zt3vvvot,zt3vvvvt)

!      use symm_mod, only: MULTB, NO, NREP, NT3VVT, NV
#include "symm.inc"
!------------------------------------------------
        integer, intent(in)  :: irp
        integer, intent(out) :: xxt3vvvt(nrep,nrep),yt3vvvt(nrep),zt3vvvot(nrep+1) 
        integer, intent(out) :: zt3vvvvt(nrep+1) 
!------------------------------------------------

!------------------------------------------------
        integer :: irep,krep,ijkrep,jrep
!------------------------------------------------

        YT3VVVT(1:nrep) = 0
        XXT3VVVT(1:nrep,1:nrep) = 0

        DO KREP = 1,NREP
          IJKREP = MULTB(KREP,IRP+NREP,2)
          XXT3VVVT(KREP,IRP)=YT3VVVT(IJKREP)
          YT3VVVT(IJKREP)=YT3VVVT(IJKREP) + NT3VVT(IRP) * NV(KREP)
        ENDDO

      ZT3VVVOT(1)=0
      DO IREP=1,NREP
        ZT3VVVOT(IREP+1) = ZT3VVVOT(IREP) + YT3VVVT(IREP)*NO(IREP)
      ENDDO

      ZT3VVVVT(1)=0
      DO IREP=1,NREP
        ZT3VVVVT(IREP+1) = ZT3VVVVT(IREP) + YT3VVVT(IREP)*NV(IREP)
      ENDDO

      end subroutine


  end module
