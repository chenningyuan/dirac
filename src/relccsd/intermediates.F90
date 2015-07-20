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
    public w_mbej
    public w_iemn
    public w_efam
    public wbar_mbej
    public w_ejmb_out_of_core
    public f_bar_ae_out_of_core
    public w_iemn_out_of_core
    public w_efam_out_of_core
    public w_efam_diagram2
    public w_efab_out_of_core
    public G_intm
     private

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface


     contains

   subroutine f_bar_mi(foo,fvo,t1,t2,ho,buf1,buf2,buf3)

#include "symm.inc"
#include "param.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in) :: foo(rcw*nfoo), fvo(rcw*nfvo)   ! Fock matrix elements
     real(8), intent(in) :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes
     real(8), intent(in) :: buf1(*),buf2(*)  ! Amplitudes
     real(8), intent(inout) :: buf3(*)
     real(8), intent(out):: ho(rcw*nfoo)                   ! Intermediate
!+++++++++++++++Local variables++++++++++++++++++++++++++++++++++++++++
     real(8), allocatable :: tau(:) !, vooo(:),vooo_sorted(:),vvoo(:)
     real*8, external     :: ddot
     integer     :: i 

!    Initialize  f_bar_mi with oo block of Fock matrix
       ho = foo

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_mi = F_bar_mi - V(nm,ef) * T(ef,in) 
!+++++++++++++++++++++++++++++++++++++++++

   allocate(tau(ndimt2*rcw))
   call gettau  (t1,t2,tau)
   call contraction_442 ((/"o3","o2","p1","p2"/),(/"p1","p2","o1","o3"/),(/"o2","o1"/),ho,-1.0d0,1.0d0,nrep,RightTensor=tau)
   deallocate(tau) 

!++++++++++++++++++++++++TERM-2+++++++++++++++++++++
!F_bar_mi(m,i) = F_bar_mi(m,i) - V(mn, ei) * T1(e,n) 
!+++++++++++++++++++++++++++++++++++++++++++++++++++

    call contraction_422 ((/"o2","o3","p1","o1"/),(/"p1","o3"/),(/"o2","o1"/),t1,ho,-1.0d0,1.0d0,nrep)

!***********************TERM-3*********
!F_bar_mi = F_bar_mi + F(m,e) * T1(e,i)
!**************************************
   
!   call contraction_222 ((/"o2","p1"/),(/"p1","o1"/),(/"o2","o1"/),fvo,t1,ho,1.0d0,1.0d0,nrep)


   do i = 1, nfoo
   if (abs(ho(i)) <= 1.0d-14) then  
   write(*,*)'value of weird hoo', ho(i)
   ho(i) = 0.0d0
   endif
   enddo

   end subroutine

   subroutine f_bar_ae(fvv,fvo,t1,t2,hv)

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
      real(8), intent(in)  :: fvv(nfvv*rcw), fvo(nfvo*rcw)
      real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
      real(8), intent(out) :: hv(nfvv*rcw)
!---------------Local variables--------------------------------------
      real(8), allocatable :: tau(:)
      real(8), allocatable :: fov(:)
      real*8      :: ddot
!---------------------------------------
! hv(a,e) = fvv(a,e)  
!---------------------------------------

!!  call xcopy (nfvv,fvv,1,hv,1)

    hv = fvv

!----------------------------------------------------------------------
! hv(a,e) = hv(a,e) - V(mn,ef) * Tau(af,mn)
!----------------------------------------------------------------------

   allocate(tau(ndimt2*rcw))
   call gettau (t1,t2,tau)

   call contraction_442 ((/"o1","o2","p2","p3"/),(/"p1","p3","o1","o2"/),(/"p1","p2"/),hv,-1.0d0,1.0d0,nrep,RightTensor=tau)
  
   deallocate(tau) 

    write(*,*)'wpp 1',ddot(nfvv,hv,1,hv,1)
!---------------------------------------------------------------
!F_bar_ae = F_bar_ae + V(am, ef) * T1(f,m) -> V(ae,fm) * T1(fm) 
!---------------------------------------------------------------

  call contraction_422 ((/"p1","o1","p2","p3"/),(/"p3","o1"/),(/"p1","p2"/),t1,hv,1.0d0,1.0d0,nrep)
 
  write(*,*)'wpp 2',ddot(nfvv,hv,1,hv,1)

!-----------------------------------------------------------------
!F_bar_ae = F_bar_ae + F(m,a) * T1(e,m)
!-----------------------------------------------------------------

   allocate (fov(nfvo*rcw))
   call srt1c1(nrep,nv,no,fvo,fov)

   call contraction_222 ((/"o1","p1"/),(/"p2","o1"/),(/"p1","p2"/),fov,t1,hv,-1.0d0,1.0d0,nrep)

   deallocate (fov)

  end subroutine

  subroutine f_bar_me(fvo,t1,hov)

#include "symm.inc"
#include "param.inc"
#include "complex.inc"
!++++++++++++++Calling variables+++++++++++++++++++++++++++++++++++++++
     real(8), intent(in),target  :: fvo(rcw,nfvo)          ! Fock matrix elements
     real(8), intent(in) :: t1(rcw,ndimt1)                 ! Amplitudes
     real(8), intent(out):: hov(rcw,nfvo)                  ! Intermediate
!+++++initialization of hov+++++++++++++
 
! call srt1c1(nrep,nv,no,fvo,hov)

  hov = fvo

!+++++++++++++Term-1++++++++++++++++++++++
! F_bar_me = F_bar_me + V(mn,ef) * T1(f,n) 
!+++++++++++++++++++++++++++++++++++++++++

 call contraction_422 ((/"o1","o2","p1","p2"/),(/"p2","o2"/),(/"o1","p1"/),t1,hov,1.0d0,1.0d0,nrep)

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
  real(8), intent(out) :: w_oooo(rcw*nv1)                ! Intermediate.
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real*8      :: ddot
  integer     :: i
!--------------------------------------------------

!--------------------
! initialize w_oooo
!--------------------

   call getoooo (w_oooo)

!*************************************************
! w_oooo(ij,mn) = w_oooo(ij,mn) + (V(ij,en)*T(e,m))
!*************************************************

   call contraction_424 ((/"o1","o2","p1","o4"/),(/"p1","o3"/),(/"o1","o2","o4","o3"/),w_oooo, &
           &      1.0d0,-1.0d0,nrep,RightTensor=t1)

   do i = 1, nv1
   if (abs(w_oooo(i)) <= 1.0d-14) then  
   w_oooo(i) = 0.0d0
   endif
   enddo

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

 subroutine w_ejmb(t1,t2,w_voov)

!!---------------Description--------------------------------------------
!!     Calculates w_voov intermediate
!!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_voov(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real(8), allocatable :: tau(:)
 real*8      :: ddot
!--------------------
! initialize w_voov
!--------------------


 call getvoov (w_voov)

!---------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) - V(nj,mb)*T1(e,n)
!---------------------------------------------------

   call contraction_424 ((/"o3","o1","o2","p1"/),(/"p2","o3"/),(/"p2","o1","o2","p1"/),w_voov, &
            &      -1.0d0,1.0d0,nrep,RightTensor=t1)

!----------------------------------------------------------------------
!   w_voov(ej,mb) = w_voov(ej,mb) + V(nj,fb)*(T2(ef,mn)-T1(f,m)*T1(e,n)
!----------------------------------------------------------------------

  allocate (tau(ndimt2*rcw))

  tau = 0.0d0

!  call gettau_approx (t1,t2,tau,'approximate')
!  call gettau_approx (t1,t2,tau)

  call contraction_444 ((/"o3","o1","p3","p1"/),(/"p2","p3","o2","o3"/),(/"p2","o1","o2","p1"/),w_voov, &
   &      1.0d0,0.0d0,nrep,RightTensor=tau)

  deallocate (tau)

!--------------------------------------------------
!  w_voov(ej,mb) = w(ej,mb) + V(ej,fb)*t(f,m) 
!--------------------------------------------------

  call contraction_424 ((/"p2","o1","p1","p3"/),(/"p3","o2"/),(/"p2","o1","o2","p1"/),w_voov, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1)

 end subroutine


!===========================================================================================================

 subroutine w_mbej(t1,t2,w_ovvo)

 external srt1l1

!---------------Description--------------------------------------------
!     Calculates w_voov intermediate
!-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Calling variables--------------------------------------
  real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
  real(8), intent(out) :: w_ovvo(rcw*ivoov(nrep+1))                ! Intermediate.
!---------------Local Variables----------------------------------------
 real(8), allocatable :: tau(:)
 real(8), allocatable :: w_vovo(:),array_temp(:)
 real*8      :: ddot
 integer :: i
!--------------------
! initialize w_voov
!--------------------


  w_ovvo = 0.0d0

  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,w_ovvo) 

  deallocate (w_vovo)

!---------------------------------------------------
!   w_ovvo(mb,ej) = w_voov(mb,ej) - V(mn,ej)*T1(b,n)
!---------------------------------------------------

   call contraction_424 ((/"o2","o3","p2","o1"/),(/"p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
            &      -1.0d0,1.0d0,nrep,RightTensor=t1)

!----------------------------------------------------------------------
!   w_ovvo(mb,ej) = w_ovvo(mb,ej) + V(mn,ef)*(T2(bf,jn)-T1(f,j)*T1(b,n)
!----------------------------------------------------------------------

! allocate (tau(ndimt2*rcw))

!  call gettau (t1,t2,tau,'approximate')
! call gettau_approx (t1,t2,tau)

! call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
!  &      1.0d0,0.0d0,nrep,RightTensor=tau)

! deallocate (tau)

!-------------------------------------------------------------------------------------
!  w_ovvo(mb,ej) = w_ovvo(mb,ej) - V(mn,ef)*(T2(fb,jn)+T1(f,j)*T1(b,n)
!-------------------------------------------------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)

  allocate(array_temp(nv4*rcw))

  array_temp = 0.0d0

  call xgemm ('n','n',ndimt1,ndimt1,1,a1,t1,ndimt1,t1,1, &
 &               a1,array_temp,ndimt1)


  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p3","o1","p1","o3"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=array_temp)

  deallocate(array_temp)
!--------------------------------------------------
!  w_ovvo(mb,ej) = w(mb,ej) + V(mb,ef)*t(f,j) 
!--------------------------------------------------

  call contraction_424 ((/"o2","p1","p2","p3"/),(/"p3","o1"/),(/"o2","p1","p2","o1"/),w_ovvo, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)


   do i = 1, ivoov(nrep+1)
   if (abs(w_ovvo(i)) <= 1.0d-14) then  
!   write(*,*)'zero of mbej',w_ovvo(i)
   w_ovvo(i) = 0.0d0
   endif
   enddo


 end subroutine

!-----=========================--------------

 subroutine w_iemn(t1,t2,w_ovoo,wbar_ovvo,w_oooo,hov)

 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t1(rcw*ndimt1),t2(rcw*ndimt2)                         ! Amplitudes.
   real(8), intent(out) :: w_ovoo(iovoot(nrep+1)*rcw)                                       ! Intermediate.
   real(8), intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),w_oooo(rcw*nv1),hov(nfvo*rcw)      ! reuse of intermediates.
 !---------------Local Variables----------------------------------------
  real*8      :: ddot
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
!   w_ovoo(ie,mn) = w_ovoo(ie,mn) + 1/2 V(ie,fg) * tau(fg,mn) 
!-----------------------------------------------------------
!**out of core term..

  allocate (tau(ndimt2*rcw))
  call gettau (t1,t2,tau)

  call contraction_444 ((/"o1","p1","p2","p3"/),(/"p2","p3","o2","o3"/),(/"o1","p1","o2","o3"/),w_ovoo, &
             &      1.0d0,1.0d0,nrep,RightTensor=tau)

  deallocate(tau)

  end subroutine

 subroutine wbar_mbej(t2,wbar_ovvo)

 !---------------Description--------------------------------------------
 !     Calculates w_voov intermediate
 !-----------------------------------------------------------------------

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real(8), intent(in)  :: t2(rcw*ndimt2)                    ! Amplitudes.
   real(8), intent(out) :: wbar_ovvo(rcw*iovvo(nrep+1))                ! Intermediate.
 !---------------Local Variables----------------------------------------
  real*8      :: ddot
  real*8, allocatable   :: w_vovo(:) 

 !--------------------
 ! initialize wbar_ovvo
 !--------------------
  wbar_ovvo = 0.0d0

  allocate (w_vovo(nv4*rcw))

  call getvovo (w_vovo)

  call srt1l1 (nrep,multb,.false.,nvo,nv,no,nvo,iovvo,iiov,w_vovo,wbar_ovvo) 

  deallocate (w_vovo)

!------------------------------------------------------------
!   wbar_ovvo(mb,ej) = wbar_ovvo(mb,ej) - V(mn,ef) * T2(bf,nj)
!------------------------------------------------------------

  call contraction_444 ((/"o2","o3","p2","p3"/),(/"p1","p3","o3","o1"/),(/"o2","p1","p2","o1"/),wbar_ovvo, &
             &      -1.0d0,0.0d0,nrep,RightTensor=t2)

 end subroutine

 subroutine w_efam(t1,t2,wbar_ovvo,hov,w_vvvo,buf2,buf3)

!--------------------DESCRIPTION---------------------------
!        calculates vvvo type intermediate in-core fashion.      
!----------------------------------------------------------

!----------calling variables----------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"

 !---------------Calling variables--------------------------------------
   real*8, intent(in)  :: t2(rcw*ndimt2),t1(rcw*ndimt1)         ! Amplitudes.
   real*8, intent(in)  :: buf2(*),buf3(*)         
   real*8, intent(in)  :: wbar_ovvo(rcw*iovvo(nrep+1)),hov(nfvo*rcw)     ! Intermediate.
   real*8, intent(out) :: w_vvvo(rcw,nv5)                ! Intermediate.
 !---------------Local Variables----------------------------------------

   real*8      :: ddot
   real*8,allocatable  :: tau(:)                ! Intermediate.
   real*8,allocatable  :: huge_temp(:),array_temp(:)                ! Intermediate.

!----------------------------------------------------------
!        Initialize
!----------------------------------------------------------

   call getvvvo_incore(w_vvvo)

       write(*,*)'w~efam: init',ddot(nv5,w_vvvo,1,w_vvvo,1)
!-----------------------------------------------------------
!   diagram1 : w_vvvo(ef,am)  = w_vvvo(ef,am) + p(ef) V(en,ag)*tau(fg,mn)  
!-----------------------------------------------------------

  allocate (tau(ndimt2*rcw))
  call gettau (t1,t2,tau)
  call contraction_444 ((/"p2","o2","p1","p4"/),(/"p3","p4","o1","o2"/),(/"p2","p3","p1","o1"/),w_vvvo, &
             &      1.0d0,0.0d0,nrep,RightTensor=t2)
  deallocate (tau)

!-----------------------------------------------------------------------
!  w_temp(en,am) = w_temp(en,am) + V(en,ag)*t1(g,m)
!  w (ef,am) = w(ef,am) - w_temp(en,am) * t(f,n) 
!-----------------------------------------------------------------------

  allocate(array_temp(nv4*rcw))

   array_temp = 0.0d0

  call contraction_424 ((/"p2","o2","p1","p4"/),(/"p4","o1"/),(/"p2","o2","p1","o1"/),array_temp, &
             &      1.0d0,1.0d0,nrep,RightTensor=t1)

  call contraction_424 ((/"p2","o2","p1","o1"/),(/"p3","o2"/),(/"p2","p3","p1","o1"/),w_vvvo, &
             &      -1.0d0,1.0d0,nrep,RightTensor=t1,LeftTensor=array_temp)

  deallocate(array_temp)
!allocate(array_temp(nv4*rcw))

! array_temp = 0.0d0

! call xgemm ('n','n',ndimt1,ndimt1,1,a1,t1,ndimt1,t1,1, &
!&               a1,array_temp,ndimt1)

! call contraction_444 ((/"p2","o2","p1","p4"/),(/"p3","o1","p4","o2"/),(/"p2","p3","p1","o1"/),w_vvvo, &
!            &      -1.0d0,0.0d0,nrep,RightTensor=array_temp)

!deallocate(array_temp)
  write(*,*)'write me inside w~efam,dia2',ddot(nv5,w_vvvo,1,w_vvvo,1)
!------------------------------------------------------------
!  diagram2 : w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m) 
!------------------------------------------------------------


!-------------------------------------------------------------
!  diagram3 : w_vvvo(ef,am) =  w_vvvo(ef,am) + f_bar_me(na) * t(ef,nm)
!-------------------------------------------------------------

      allocate(huge_temp(kvvovt(nrep+1)*rcw))

      huge_temp = 0.0d0
      call srt1t3 (nrep,multb,.false.,nvvt,no,no,nvvot,kvvoot,kkvvot, &
     &             t2,buf2)
  
      call xcopy(kvvovt(nrep+1),a0,0,buf3,1)

      call cntrct ('c','c',nv,nvvot,no,a1,hov,buf2,a0,buf3,nrep)

      call srt1c1 (nrep,nv,nvvot,buf3,huge_temp)

      call srt19 (nrep,multb,.true.,nvvt,nv,no,nvvot,kvvovt,kkvvot,huge_temp,buf3)

      call xaxpy (nv5,a1,buf3,1,w_vvvo,1)

      deallocate(huge_temp)

!     call contraction_244 ((/"o2","p1"/),(/"p2","p3","o2","o1"/),(/"p2","p3","p1","o1"/),w_vvvo, &
!    &      -1.0d0,1.0d0,nrep,LeftTensor=hov,RightTensor=t2)

  write(*,*)'write me inside w~efam,dia3',ddot(nv5,w_vvvo,1,w_vvvo,1)
!---------------------------------------------------------------------------
!  diagram4 : w_vvvo(ef,am) =  w_vvvo(ef,am) + 0.5d0 * V(no,am) * tau(ef,no) 
!---------------------------------------------------------------------------

  allocate (tau(ndimt2*rcw))
  call gettau (t1,t2,tau)
  call contraction_444 ((/"o2","o3","p1","o1"/),(/"p2","p3","o2","o3"/),(/"p2","p3","p1","o1"/),w_vvvo, &
             &      1.0d0,0.0d0,nrep,RightTensor=tau)
  deallocate (tau)

   write(*,*)'w~efam: vhhph',ddot(nv5,w_vvvo,1,w_vvvo,1)
!---------------------------------------------------------------------------
!  diagram5 : w_vvvo(ef,am) =  w_vvvo(ef,am) - w_bar_ovvo(nf,am) * t(e,n)
!---------------------------------------------------------------------------

  call contraction_424 ((/"o2","p3","p1","o1"/),(/"p2","o2"/),(/"p2","p3","p1","o1"/),w_vvvo, &
            &      -1.0d0,1.0d0,nrep,LeftTensor=wbar_ovvo,RightTensor=t1)

      write(*,*)'w~efam: whpph',ddot(nv5,w_vvvo,1,w_vvvo,1)

 end subroutine

 subroutine w_efam_diagram2(t1,t2,buf1,buf2,buf3,nbuf3,w_vvvo)

#include "param.inc"  
#include "complex.inc"
#include "symm.inc"   
#include "ccpar.inc"
!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)  ! Amplitudes.
   real(8), intent(inout) :: w_vvvo(rcw*nv5)                ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*)                ! buffer arrays
   real(8), intent(inout) :: buf3(*)                        ! buffer arrays
   integer, intent(in)    :: nbuf3 

!---------------Local variables--------------------------------------
                                                                    
      real*8  :: t1val(2)                                                
      logical :: done
      integer :: cd,cco,istart,nint,akcd,ci,akci,t1di,di,akdi,t1ci,cmin,c,d,i
      integer :: crep,drep,irp,mint
      real*8,allocatable  :: tau (:)
      real*8,allocatable  :: array_temp(:) !huge
!------------------------------------------------------------------
!  w_vvvo(ef,am) =  w_vvvo(ef,am) + w (ef,ag)*t(g,m)
!------------------------------------------------------------------

      t1val = (/0.0d0,0.0d0/)

     allocate(array_temp(nv5*rcw))

     allocate (tau(ndimt2*rcw))
     call gettau (t1,t2,tau)

      if(myproc.eq.master) then
        call xcopy(nv5,w_vvvo,1,array_temp,1)
      else 
        call xcopy(nv5,a0,0,array_temp,1)
      endif

      do 30 irp = 1, nrep
         if (nvv(irp).eq.0) goto 30
         done = .false.
         istart = 0
         cco=0

! ** w_efab are constructed over the batches.
! ** ISTART is updated with the actual batch #
! ** available on the local node, NINT is updated with the number of
! ** actually read integrals from this batch and is a pure output parameter !
! ** We therefore have to count upwards in case of ISTART.NE.0 !!

         call w_efab_out_of_core (irp,istart,done,tau,buf1,nbuf3,buf3,nint)

         mint = nbuf3/nvo(irp)
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
                    call w_efab_out_of_core (irp,istart,done,tau,buf1,nbuf3,buf3,nint)
                      CD = 1
                      CCO=ISTART
                    ELSE
                      GOTO 30
                    ENDIF
                 ENDIF
                 AKCD = (CD-1)*NVVT(IRP)*RCW+1
                 DO I = 1, NO(DREP)
                    CI = IIVO(CREP,DREP)+(I-1)*NV(CREP)+C
                    AKCI = (IOVVVT(IRP)+(CI-1)*NVVT(IRP))*RCW+1
                    T1DI = (IVO(DREP)+(I-1)*NV(DREP)+D-1) * RCW + 1
                    T1VAL(1) =  T1(T1DI)
                    IF (CARITH) T1VAL(2) =  T1(T1DI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,array_temp(AKCI),1)
                 ENDDO
                 DO I = 1, NO(CREP)
                    DI = IIVO(DREP,CREP)+(I-1)*NV(DREP)+D
                    AKDI = (IOVVVT(IRP)+(DI-1)*NVVT(IRP))*RCW+1
                    T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
                    T1VAL(1) = -T1(T1CI)
                    IF (CARITH) T1VAL(2) = -T1(T1CI+1)
                    CALL XAXPY(NVVT(IRP),T1VAL,BUF3(AKCD),1,array_temp(AKDI),1)
                 ENDDO
 10           CONTINUE
 15        CONTINUE
 20      CONTINUE
 30    CONTINUE

   deallocate (tau)

   call xcopy(nv5,array_temp,1,w_vvvo,1)

   deallocate(array_temp)
 end subroutine


 subroutine w_ejmb_out_of_core(t1,t2,buf1,buf2,buf3,nbuf3,w_vovo)


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
                                                                     
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2)! Amplitudes.
   real(8), intent(inout) :: w_ovoo(rcw,iovoot(nrep+1))              ! Intermediate.
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

   subroutine w_efam_out_of_core (t1,t2,hov,w_vovo,wbar_ovvo,buf1,buf2,buf3,nbuf3,nbuf2,w_vvvo)

!---------------Common Blocks--------------------------------------   
                                                                     
#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
!#include "eqns.inc"                                                   
#include "ccpar.inc"                                                  

!---------------Calling variables-------------------------------------

   real(8), intent(in)    :: t1(rcw*ndimt1),t2(rcw*ndimt2) ! Amplitudes.
   real(8), intent(in)    :: hov(nfvo*rcw)       
   real(8), intent(in)    :: w_vovo(rcw*nv4),wbar_ovvo(rcw*nv4)       
   real(8), intent(out)   :: w_vvvo(rcw*nv5)             ! Intermediate.
   real(8), intent(in)    :: buf1(*),buf2(*)             ! buffer arrays.
   real(8), intent(in),target    :: buf3(*)              ! buffer arrays.
   integer, intent(in)    :: nbuf3,nbuf2 

!---------------Local variables----------------
      logical :: done,teq
      integer :: kd,istart,nint,m,n,k
      integer :: irp,mint,off1,off2,off3
      integer :: a,ab,abkj,abkl,ak,akcj,amin,arep,b,bk,bkcj,c,cd
      integer :: cj,ck,cmin,crep,brep,d,dj,drep,j,i,jrep,kmin,krep 
      integer :: AllocateStatus 
      real*8,allocatable :: tau(:)
      real*8,allocatable :: t_sort(:,:)
      real*8,allocatable :: w_vvvv(:)
      real*8             :: t1con(2)
      real*8, pointer    :: xx(:,:) 
!----------------------------------------------
      integer :: ijkl, ik0,il0,imin,ip0,iq0,irep,kco,kj,kl,klstar,lj,ljrep, &
       &         lq0,lrep,p,pl0,pq,pq0,pqrep,prep,l,q,qrep,kjrep,lqrep,lqip, & 
       &         ilrep,iqrep,ilpq,pliq 
      real*8  :: tval(2) 


!-------------Initialization---------------

      t1con = (/0.0d0,0.0d0/)
      tval = (/0.0d0,0.0d0/)

      do i = 1,nrep
       xx(1:mvo(i),1:mvv(i)) => buf3 (1:mvo(i)*mvv(i)) 
      enddo

      call xcopy(nv5,A0,0,w_vvvo,1)

!     allocate (t_sort(j2vovo(nrep+1)*rcw))
     allocate (t_sort(sum(mvo,nrep),sum(mvo,nrep)))

     t_sort = 0.0d0

      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,mvo,j2vovo,jjvo,jjvo, &
      &              t2,t_sort)

      allocate (tau(ndimt2),stat = AllocateStatus)
      if (AllocateStatus /= 0) stop "*** not enough memory ***"

      call gettau (t1,t2,tau)

      off3 = 1
      off1 = 1

      done = .false.

      DO 50 IRP = 1, NREP
         IF (NOV(IRP).EQ.0.OR.NVVT(IRP).EQ.0) GOTO 50
         MINT = MIN0(NBUF3/NVVT(IRP),NBUF2/NVO(IRP))
         ISTART = 0
 10      CALL GETVOVV (IRP,ISTART,NINT,DONE,BUF2,MINT)

         allocate (w_vvvv(nint*nvvt(irp)),stat = AllocateStatus)!hopefully ok. not 100% sure.
         if (AllocateStatus /= 0) stop "*** not enough memory ***"


         w_vvvv = 0.0d0

         call w_efab_out_of_core (irp,istart,done,tau,buf1,nbuf3,w_vvvv,nint)

!-----------------------------------------------------------------------
!        w(ef,am) = w(ef,am) + w_vvvv(ef,ag)*t(gm)
!-----------------------------------------------------------------------
         CD = 1
         DO DREP = 1, NREP
            CREP = MULTB(DREP,IRP+NREP,2)
            IF (CREP.LT.DREP) cycle
            DO D = 1, NV(DREP)
               CMIN = 1
               IF (CREP.EQ.DREP) CMIN = D + 1
               DO C = CMIN, NV(CREP)
                  JREP = DREP
                  DO J = 1, NO(JREP)
                     DJ = (IVO(JREP)+(J-1)*NV(DREP)+D-1) * RCW + 1
                     CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                     T1CON(1) =  T1(DJ)
                     IF (CARITH) T1CON(2) = - T1(DJ+1)
                     CALL XAXPY(NINT,T1CON,w_vvvv(CD),1, &
                      &         BUF2(CJ),NVO(IRP))
                  ENDDO
                  JREP = CREP
                  DO J = 1, NO(JREP)
                     CJ = (IVO(JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                     DJ = (IIVO(DREP,JREP)+(J-1)*NV(DREP)+D-1) * RCW + 1
                     T1CON(1) =  - T1(CJ)
                     IF (CARITH) T1CON(2) = T1(CJ+1)
                     CALL XAXPY(NINT,T1CON,w_vvvv(CD),1, &
                       &       BUF2(DJ),NVO(IRP))
                  ENDDO
                  CD = CD + RCW
               ENDDO
            ENDDO
         ENDDO
         
         deallocate (w_vvvv)


!!=================================term to investigate==========================
!-------------------------------------------------------
!new term:  w(ef,am) = w(ef,am) + p(ef) V(en,ag)*t(gf,mn)
!-------------------------------------------------------


   !   IJKL = 1
   !   KL = 0
   !   KCO = 0
   !   KLSTAR = 0

   !   DO LREP = 1, NREP
   !    KREP = MULTB(LREP,IRP+NREP,2)
   !    IF (KREP.LT.LREP) cycle
   !    DO L = 1, NV(LREP)
   !     KMIN = 1
   !     IF (KREP.EQ.LREP) KMIN = L + 1
   !     DO K = KMIN, NV(KREP)

   !      IF(KCO.LT.KLSTAR) THEN
   !        KCO = KCO + 1
   !        cycle
   !      ENDIF

   !      KL = KL + 1
 
   !      IF (KL.EQ.NINT+1) THEN
!!-!---------------------
!!W! need the next buffer
!!-!---------------------
   !         KLSTAR = KLSTAR + NINT
   !         NINT = MINT
   !            IF(.NOT.DONE) THEN
   !              CALL GETVOVV (IRP,KLSTAR,NINT,DONE,BUF1,MINT)
   !            ELSE
   !              GOTO 50
   !            ENDIF

 ! !         KL = 1
 ! !         IJKL = 1
 ! !         KCO = KLSTAR
   !      ENDIF

   !         DO  JREP = 1, NREP
   !          IREP = MULTB(JREP,IRP+NREP,2)
   !          LJREP = MULTB(LREP,JREP,2)
   !          KJREP = MULTB(KREP,JREP,2)
   !          LJ    = JJVO(LREP,JREP) + L - 1
   !          LJ    = LJ * RCW + 1
   !          KJ    = JJVO(KREP,JREP) + K - 1
   !          KJ    = KJ * RCW + 1

   !          IK0   = JJVV(IREP,KREP) + (K-1)*NV(IREP)
   !          IK0   = IK0 * RCW + 1
   !          IL0   = JJVV(IREP,LREP) + (L-1)*NV(IREP)
   !          IL0   = IL0 * RCW + 1

   !          PQ = NVO(LJREP)

   !        DO  J = 1, NO(JREP)

!!-!----------new loop over (P,Q) to get all the elements of sorted t2 array------------

   !           CALL XGEMM ('N','N',NV(IREP),PQ,1,A1,BUF1(IJKL),NVV(IRP),t_sort(LJ,OFF3), &
   !           &               NVO(IRP),A1,xx(IK0,OFF3),NVV(IRP))

   !              CALL XCOPY (PQ,t_sort(KJ,OFF3),1,TVAL,1)
   !              CALL XSCAL (PQ,-A1,TVAL,1)

   !              CALL XGEMM ('N','N',NV(IREP),PQ,1,A1,BUF1(IJKL),NVV(IRP),TVAL, &
   !              &               NVO(IRP),A1,xx(IL0,OFF3),NVV(IRP))

   !           IJKL = IJKL + NV(IREP) * RCW
   !           LJ   = LJ   + NV(LREP) * RCW
   !           KJ   = KJ   + NV(KREP) * RCW

   !          enddo

   !           OFF3 = OFF3 + PQ * RCW
   !         enddo

   !     enddo
   !    enddo
   !   enddo

!!-!------------------------------------------------
!! ! sorting step. major caveat: buf3 is not triangular.
!  !IL,PQ) -> (LQ,IP)-(LQ,PI) -> (LQ,I>P)  
!!-!------------------------------------------------

   !     DO PQREP = 1,NREP
   !      DO  PREP = 1,NREP
   !          QREP = MULTB(PREP,PQREP+NREP,2) 
   !       DO Q = 1, NO (QREP) 

   !         DO P = 1, NV (PREP)

   !          DO IREP = 1, NREP 
   !             LREP = MULTB(IREP,PQREP+NREP,2) 

   !        
   !           LQREP = MULTB(LREP,QREP,2)
   !           ILREP = MULTB(IREP,LREP,2)
   !           IQREP = MULTB(IREP,QREP,2)
   !           PQ0   = JJVV(PREP,QREP) + (Q-1)*NV(PREP)+P
   !           PQ0   = PQ0 * RCW + 1 

   !           PL0   = JJVV(PREP,LREP) +  P
   !           PL0   = PL0 * RCW + 1


   !           IMIN = 1
   !           IF (IREP.EQ.PREP) IMIN = P + 1
   !           DO I = IMIN, NV(IREP)


   !           IP0   = JJVV(IREP,PREP) + (P-1)*NV(IREP)+I
   !           IP0   = IP0 * RCW + 1 

   !           IL0   = JJVV(IREP,KREP) +  I
   !           IL0   = IK0 * RCW + 1

   !           IQ0   = JJVO(KREP,QREP) + (Q-1)*NV(IREP) + I
   !           IQ0   = IQ0 * RCW + 1

   !           LQ0   = JJVO(LREP,QREP) + (Q-1)*NV(LREP)
   !           LQ0   = LQ0 * RCW + 1


   !           LQIP = (IP0-1)*MVO(LQREP)  + LQ0
!  !            ILPQ = (PQ0-1)*MVV(ILREP)  + IL0
!  !            PLIQ = (PL0-1)*MVV(IQREP)  + IQ0

   !           CALL XCOPY (NV(LREP),xx(IL0,PQ0),NV(IREP),buf2(LQIP),1)

   !           CALL XAXPY (NV(LREP),-A1,xx(PL0,IQ0),NV(PREP),buf2(LQIP),1)


   !           LQ0 = LQ0 + RCW*NV(LREP)
   !           enddo
   !           IL0 = IL0 + RCW*NV(LREP) 
   !           PL0 = PL0 + RCW*NV(LREP) 
   !          enddo
   !         enddo

   !       enddo
   !      enddo
   !     enddo

!             DO I = 1, NV (PREP) 
!              PMIN = 1
!              IF (IREP.EQ.PREP) PMIN = I + 1
!              DO P = PMIN, NV(PREP)

!              IK0   = JJVV(IREP,KREP) + (K-1)*NV(IREP) + I
!              IK0   = IK0 * RCW + 1

!              PI0   = JJVV(IREP,PREP) + (I-1)*NV(IREP)+P
!              PI0   = PI0 * RCW + 1 

!              IF (LJREP.EQ.1) THEN
!!               CALL XCOPY (1,BUF2(IK0,OFF3),1,buf3(KQ0,PI0),1)
!              CALL XAXPY (1,-A1,BUF2(IK0,OFF3),1,buf3(KQ0,PI0),1)
!              ENDIF

!              IF (KJREP.EQ.1) THEN
!!               CALL XCOPY (1,BUF2(IL0,OFF3),1,buf3(LQ0,PI0),1)
!              CALL XAXPY (1,-A1,BUF2(IL0,OFF3),1,buf3(LQ0,PI0),1)
!              ENDIF

!              enddo
!             enddo


!!================================ term to investigate ==================================

!-------------------------------------------------
!    w(ef,am) = w(ef,am) + 1/2 V(no,am)*tau(ef,no)
!-------------------------------------------------

         CALL GETVOOO (BUF3)
         M = NVO(IRP)
         N = NINT
         K = NOOT(IRP)
         ABKL = (IVVOOTT(IRP)+ISTART)*RCW+1
         CALL XGEMM ('N','C',M,N,K,A1,BUF3,M,tau(ABKL),NVVT(IRP), &
     &               A1,BUF2,M)

!-------------------------------------------
!   w(ef,am) = w(ef,am) + f_bar(na)*t(ef,mn) 
!-------------------------------------------

         ABKJ = (IVVOOTT(IRP)+ISTART)*RCW+1
         DO 30 JREP = 1, NREP
            KREP = MULTB(JREP,IRP+NREP,2)
            IF (KREP.LT.JREP) GOTO 30
            DO J = 1, NO(JREP)
               KMIN = 1
               IF (KREP.EQ.JREP) KMIN = J + 1
               DO K = KMIN, NO(KREP)
                  CREP = KREP
                  DO C = 1, NV(CREP)
                     CK = (IVO(KREP)+(K-1)*NV(CREP)+C-1) * RCW + 1
                     CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                     T1CON(1) =  - HOV(CK)
                     IF (CARITH) T1CON(2) = HOV(CK+1)
                     CALL XAXPY(NINT,T1CON,T2(ABKJ),1, &
                       &        BUF2(CJ),NVO(IRP))
                  ENDDO
                  CREP = JREP
                  DO C = 1, NV(CREP)
                     CJ = (IVO(JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
                     CK = (IIVO(CREP,KREP)+(K-1)*NV(CREP)+C-1) * RCW + 1
                     T1CON(1) =  HOV(CJ)
                     IF (CARITH) T1CON(2) = - HOV(CJ+1)
                     CALL XAXPY(NINT,T1CON,T2(ABKJ),1, &
                       &        BUF2(CK),NVO(IRP))
                  ENDDO
                  ABKJ = ABKJ + NVVT(IRP) * RCW
               ENDDO
            ENDDO
 30      CONTINUE

!------------------------------------------------
!    w(ef,am) = w(ef,am) + w_bar(nf,am)*t(en) 
!------------------------------------------------

         AB = 0
         DO 40 BREP = 1, NREP
            AREP = MULTB(BREP,IRP+NREP,2)
            IF (AREP.LT.BREP) GOTO 40
            DO B = 1, NV(BREP)
               AMIN = 1
               IF (AREP.EQ.BREP) AMIN = B + 1
               DO A = AMIN, NV(AREP)
                  IF (AB.LT.ISTART) GOTO 40
                  IF (AB.GE.ISTART+NINT) GOTO 41
                  KREP = BREP
                  AKCJ = (IVOVO(IRP)+IIVO(AREP,KREP)+A-1)
                  DO K = 1, NO(KREP)
                     BK = (IVO(KREP)+(K-1)*NV(BREP)+B-1) * RCW + 1
                     T1CON(1) = - T1(BK)
                     IF (CARITH) T1CON(2) = T1(BK+1)
                     CALL XAXPY (NVO(IRP),T1CON,wbar_ovvo(AKCJ),NVO(IRP), &
                       &         BUF2(AB*RCW+1),1)
                     AKCJ = AKCJ + NV(AREP) * RCW
                  ENDDO
                  KREP = AREP
                  BKCJ = (IVOVO(IRP)+IIVO(BREP,KREP)+B-1)
                  DO K = 1, NO(KREP)
                     AK = (IVO(KREP)+(K-1)*NV(BREP)+A-1) * RCW + 1
                     T1CON(1) = T1(AK)
                     IF (CARITH) T1CON(2) = - T1(AK+1)
                     CALL XAXPY (NVO(IRP),T1CON,wbar_ovvo(BKCJ),NVO(IRP), &
                      &          BUF2(AB*RCW+1),1)
                     BKCJ = BKCJ + NV(BREP) * RCW
                  ENDDO
                  AB = AB + 1
               ENDDO
            ENDDO
   40    CONTINUE
   41    CONTINUE

        call xcopy(nint*nvo(irp),buf2,1,w_vvvo(istart*nvo(irp)*rcw+off1),1)

         IF (.NOT.DONE) THEN
            ISTART = ISTART + NINT
            GOTO 10
         ENDIF
         off1 = off1 + nvo(irp)*nvvt(irp)*rcw
 50   CONTINUE

    deallocate(tau)
    deallocate(t_sort)

   end subroutine

   subroutine w_efab_out_of_core (irp,istart,done,tau,buf1,nbuf3,w_vvvv,nint)

#include "param.inc"                                                  
#include "complex.inc"                                                
#include "symm.inc"                                                   
#include "ccpar.inc"                                                  

!--------------------calling variables-------------
   real(8), intent(in)    :: tau(:)  ! Amplitudes.
   real(8), intent(out)   :: w_vvvv(*)                ! Intermediate.
   real(8), intent(in)    :: buf1(*)                  ! buffer arrays
   integer, intent(in)    :: nbuf3,irp,istart
   integer, intent(out)   :: nint
   logical, intent(inout) :: done
!--------------------------------------------------


!--------------------local variables----------------
   integer :: mint,off1
   integer :: abkl,m,n,k
   real*8,allocatable :: w_vvvv_temp(:,:)  
!---------------------------------------------------

         off1 = 1

         mint = nbuf3/nvvt(irp)
         call getvvvv (irp,istart,nint,done,w_vvvv,mint)
         if (carith) then
            call conjuga (nint*nvvt(irp),w_vvvv,1)
         endif
   
!----------------------------------------------
! w(ef,ab) =  V(ef,ab) + V(mn,ab) * TAU*(ef,mn)
!----------------------------------------------

       allocate(w_vvvv_temp(nint,nvvt(irp))) !be careful about rcw. probably this is OK. 

       w_vvvv_temp = 0.0d0

       call getvvoo (buf1)
       m = nvvt(irp)
       n = nint
       k = noot(irp)
       abkl = (ivvoott(irp)+istart)*rcw+1
       off1 = off1 + ivvoott(irp)*rcw  

       call xgemm ('n','c',n,m,k,a1,buf1(abkl),m,tau(off1),nvvt(irp), &
     &               a0,w_vvvv_temp,m)

       call xaxpy(nint*nvvt(irp),a1,transpose(w_vvvv_temp),1,w_vvvv,1)

       deallocate(w_vvvv_temp)

   end subroutine

  subroutine G_intm (t2,l2,G_vv,G_oo)

!---------------Common Blocks--------------------------------------
#include "param.inc"
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
!    G_ea(e,a)  = -1/2 l2(mn,af) * t2(ef,mn)
!----------------------------------------------

  call contraction_442 ((/"o1","o2","p1","p3"/),(/"p2","p3","o1","o2"/),(/"p2","p1"/),G_vv,-1.0d0,1.0d0,nrep,&
                        &  LeftTensor=l2,RightTensor=t2 )

!---------------------------------------------
!     G_im(i,m) = 1/2 l2(in,ef) * t2(ef,mn)
!---------------------------------------------

  call contraction_442 ((/"o1","o3","p1","p2"/),(/"p1","p2","o2","o3"/),(/"o1","o2"/),G_oo,1.0d0,1.0d0,nrep,&
                        & LeftTensor=l2,RightTensor=t2)

  end subroutine

  end module
