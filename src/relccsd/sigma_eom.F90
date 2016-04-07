  module sigma_eom

!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

  use contraction
  use intermediates_1b_2b

!    public eom_ip  
!    public eom_ea
    public sigma_ee_1h1p
    public IM_container

    private
    public intermediates
    type intermediates
      real(8),allocatable,dimension(:) :: Fbar_mi, Fbar_me, Fbar_ae                                         ! all 1-b intermediates 
      real(8),allocatable,dimension(:) :: W_ijmn, W_ejmb, W_mbej, W_iemn, W_efam, W_amef, W_mnie, Wbar_mbej ! all 2-b intermediates   
    endtype intermediates


   contains


   subroutine IM_container(eomtype,t1,t2,A,buf1,buf2,buf3)

#include "symm.inc"
#include "complex.inc"
!-----------------------------------------------------------------------------   
!  In this routine we calculate all the intermediates required for a specific
!  problem (IP/EA/EE). 
!-----------------------------------------------------------------------------

!-------------------------------Variables--------------------------------------

!Calling Variables
!----------------- 
character(8), intent (in) :: eomtype  
type(intermediates), intent(out) :: A
real(8), intent(in) :: t1(:), t2(:)
real(8), intent(in) :: buf1(*),buf2(*)
real(8), intent(inout) :: buf3(*)

!Local Variables
!----------------
real(8), allocatable :: foo(:),fvo(:),fvv(:) 

!-----------------------------------------------------------------------------

     allocate (foo(nfoo*rcw))
     allocate (fvo(nfvo*rcw))
     allocate (fvv(nfvv*rcw))


     call fmtofile (.false.,fvo,foo,fvv) ! read one-body integrals from file.

     call f_bar_mi(foo,fvo,t1,t2,A%Fbar_mi)

     call f_bar_me(fvo,t1,A%Fbar_me)

     call f_bar_ae(fvv,fvo,t1,t2,A%Fbar_ae)

     call w_ijmn(t1,t2,A%W_ijmn)

     call wbar_mbej(t2,A%Wbar_mbej)

     call w_mbej(t1,t2,A%W_mbej)

     call w_iemn(t1,t2,A%W_iemn,A%Wbar_mbej,A%W_ijmn,A%Fbar_me)

     if (eomtype /= "IP") then

      call w_efam(t1,t2,A%Wbar_mbej,A%Fbar_me,A%W_efam,buf2,buf3)

      call w_efam_diagram2 (t1,t2,buf1,buf2,buf3,nbuf3,A%W_efam)

     endif 

     call w_mnie(t1,A%W_mnie)

     if (eomtype /= "IP") then

      call w_amef(t1,A%W_amef)

     endif 

   deallocate(foo)
   deallocate(fvo)
   deallocate(fVV)


   end subroutine

   subroutine sigma_ee_1h1p(r1,r2,t1,t2,B,sigma1,sigma2,buf2,buf3)
#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Description--------------------------------------------
!     Calculates 1h1p type sigma vectors for excitation energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),target   :: r1(:,:), r2(:,:)           ! R-Amplitudes
  real(8), intent(in)   :: t1(:),t2(:)           ! R-Amplitudes
  real(8), intent(in)  :: buf2(*),buf3(*)           ! T-Amplitudes
  type(intermediates), intent(in)  :: B          ! T-Amplitudes
  real(8), intent(out),target  :: sigma1(:,:), sigma2(:,:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), pointer :: r1_ptr(:),r2_ptr(:),sigma1_ptr(:),sigma2_ptr(:)

  real(8), allocatable :: int_1b(:),sigma2_temp(:) 
  integer :: off1,off2,irp,istart,mint,m,n,k
  logical :: done 
!--------------------------------------------------

  r1_ptr(1:ndimt1*rcw)     => r1(1:ndimt1*rcw,1)
  r2_ptr(1:ndimt2*rcw)     => r2(1:ndimt2*rcw,1)
  sigma1_ptr(1:ndimt1*rcw) => sigma1(1:ndimt1*rcw,1)
  sigma2_ptr(1:ndimt2*rcw) => sigma2(1:ndimt2*rcw,1)


! class : (HBar_SS*R)(a,i)
!=========================

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),B%Fbar_ae,r1_ptr,sigma1_ptr,1.0d0,1.0d0,nrep)

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

    call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),B%Fbar_mi,r1,sigma1,-1.0d0,1.0d0,nrep)

!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

    call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),r1,sigma1,1.0d0,1.0d0,nrep,B%W_mbej)


! class : [HBar_SD*R](a,i)
!=========================

!S_ai = S_ai + F_bar_me(m,e)*r2(ae,im)
!--------------------------------------

    call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), &
   &   (/"p1","o1"/),B%Fbar_me,sigma1,1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

!S_ai = S_ai + w_amef(am,ef) * r2(ef,im)
!----------------------------------------
 
    call contraction_442((/"p1","o2","p2","p3"/),(/"p2","p3","o1","o2"/),(/"p1","o1"/),sigma1,1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_amef,RightTensor=r2_ptr)


!S_ai = S_ai - W_mnie(mn,ie)*r2(ae,mn)
!---------------------------------------

    call contraction_442((/"o2","o3","o1","p2"/),(/"p1","p2","o2","o3"/),(/"p1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_mnie,RightTensor=r2_ptr)


! class : [HBar_DS*R](ab,ij)
!===========================

!S(ab,ij) = S(ab,ij) + P(ab) w_iemn(ma,ij)*r1(b,m)
!--------------------------------------------------

    call contraction_424 ((/"o2","p1","o1","o2"/),(/"p2","o2"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%W_iemn, &
   &    RightTensor=r1_ptr) 

!S(ab,ij) = S(ab,ij) + P(ij) w_efam(ab,ej)*r1(e,i)
!--------------------------------------------------

    call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","o1"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%W_efam, &
   &    RightTensor=r1_ptr) 

!int_1b(b,f) = int_1b(b,f) + w_amef (bm,fe)*r1(e,m)
!S(ab,ij) = S(ab,ij) + P(ab) int_1b(b,f)*t2(af,ij)
!---------------------------------------------------

    allocate(int_1b(nfvv*rcw))


    int_1b = 0.0d0

    call contraction_422 ((/"p2","o3","p4","p3"/),(/"p3","o3"/),(/"p2","p4"/),r1_ptr,int_1b,1.0d0,1.0d0,nrep,B%W_amef)

    call contraction_244 ((/"p2","p4"/),(/"p1","p4","o1","o2"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)

    deallocate(int_1b)


!int_1b(n,j) = int_1b(n,j) + w_mnie(nm,je)*r1(e,m)
!S(ab,ij) = S(ab,ij) + int_1b(n,j) * t2(ab,in)
!---------------------------------------------------

    allocate(int_1b(nfoo*rcw))

    int_1b = 0.0d0

    call contraction_422 ((/"o4","o3","o2","p3"/),(/"p3","o3"/),(/"o4","o2"/),r1,int_1b,1.0d0,1.0d0,nrep,LeftTensor=B%W_mnie)

    call contraction_244 ((/"o4","o2"/),(/"p1","p2","o1","o4"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)

    deallocate(int_1b)

! class : [HBar_DD*R](ab,ij) 
!===========================

! S(ab,ij) = S(ab,ij) + P(ab) f_bar_ae(b,e) * r2(ae,ij)
!-------------------------------------------------------

    call contraction_244 ((/"p2","p3"/),(/"p1","p3","o1","o2"/), &
   & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_ae,RightTensor=r2_ptr)
     
! S(ab,ij) = S(ab,ij) - P(ij) f_bar_mi(m,j)*r2(ab,im) 
!-----------------------------------------------------

    call contraction_244 ((/"o3","o2"/),(/"p1","p2","o1","o3"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_mi,RightTensor=r2_ptr)

! S(ab,ij) = S(ab,ij) + w_vvvv(ab,ef)*r2(ef,ij)
!------------------------------------------------

      off1 = 1
      off2 = 1
      do irp = 1, nrep
         istart = 0
         if (nvvt(irp).eq.0) cycle
         mint = nbuf3/nvvt(irp)
  10     call getdz (irp,istart,nint,done,buf3,mint)
         m = nvvt(irp)
         n = noot(irp)
         k = nint
         call xgemm ('n','n',m,n,k,a1,buf3,m,r2_ptr(off1+istart*rcw),m, &
     &               a1,sigma2_ptr(off2),m)
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off1 = off1 + m * n * rcw
         off2 = off2 + m * n * rcw
      enddo

! S(ab,ij) = S(ab,ij) + w_ijmn(mn,ij) * r2(ab,mn)
!--------------------------------------------------

      call contraction_444((/"o3","o4","o1","o2"/),(/"p1","p2","o1","o2"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_ijmn)

! S(ab,ij) = S(ab,ij) + P(ab)P(ij) w_mbej(mb,ej)*r2(ae,im)
!----------------------------------------------------------

      call srt1tt5 (nrep,multb,.false.,.false.,nv,nv,no,no,mvo,j2voov,jjvo,jjov, &
     &              r2,buf2)


      call srt1ss4 (nrep,multb,.false.,.false.,no,nv,nv,no,mov,j2ovvo,jjov,jjvo, &
     &              B%W_mbej,buf3)


      allocate (sigma2_temp(j2vovo(nrep+1)*rcw))

      resl2_temp = 0.0d0


      do i = 1, nrep
         j = multb(i+nrep,1+nrep,2)

          m = mvo(i)
          n = mvo(i)
          k = mov(j)

          nd = mvo(i)
          off1 = 1 + j2voov(i)*rcw  
          off2 = 1 + j2ovvo(i)*rcw  
          off3 = 1 + j2vovo(i)*rcw  


            call xgemm ('n','t',m,n,k,a1,buf2(off1),m, &
      &                  buf3(off2),nd,a0,sigma2_temp(off3),m)

      enddo

      call srt1tt4 (nrep,multb,.false.,.true.,nv,nv,no,no,mvo,j2vovo,jjvo,jjvo, &
     &              sigma2_temp,buf2)

      deallocate(sigma2_temp)

      call xaxpy (ndimt2,a1,buf2,1,sigma2,1)


!int_1b(a,f) = int_1b(a,f) - V(nm,fe)*r2(ea,mn)  
!S(ab,ij) = S(ab,ij) + int_1b(a,f)*t2(fb,ij)  
!-----------------------------------------------
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p3","p1","o3","o4"/),(/"p1","p4"/),int_1b,-1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    call contraction_244 ((/"p1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)

!int_1b(i,n) = int_1b(i,n) + V(nm,fe)*r2(fe,im)
!S(ab,ij) = S(ab,ij) - int_1b(i,n)*t(ab,jn)
!------------------------------------------------

    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p4","p3","o1","o3"/),(/"o1","o4"/),int_1b,1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    call contraction_244 ((/"o1","o4"/),(/"p1","p2","o2","o4"/), &
    & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)


  end subroutine

   end module 
