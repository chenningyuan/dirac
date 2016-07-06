  module sigma_eom

  use contraction
  use intermediates_1b_2b
  implicit none

!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

!    public eom_ip  
!    public eom_ea
    public sigma_ee_1h1p
    public sigma_IP
    public IM_container
    public create_sigma_vector
    public create_diagonal
    public initialize_buffers
    public free_buffers
    public get_matrix_size

    private
    public intermediates
    type intermediates
      character(8) :: eomtype
      real(8),allocatable,dimension(:) :: Fbar_mi, Fbar_me, Fbar_ae                                         ! all 1-b intermediates 
      real(8),allocatable,dimension(:) :: W_ijmn, W_ejmb, W_mbej, W_iemn, W_efam, W_amef, W_mnie, Wbar_mbej ! all 2-b intermediates   
      real(8),allocatable,dimension(:) :: t1, t2
    endtype intermediates

    real(kind=8), pointer :: local_buf1(:) 
    real(kind=8), pointer :: local_buf2(:) 
    real(kind=8), pointer :: local_buf3(:) 

   contains

   function get_matrix_size(A)
      integer :: get_matrix_size
      type(intermediates), intent(in) :: A

      get_matrix_size = size(A%t1) + size(A%t2)
   end function

   subroutine initialize_buffers(buf1, buf2, buf3)
       real(kind=8), allocatable, target ::  buf1(:), buf2(:), buf3(:)
       local_buf1 => buf1 
       local_buf2 => buf2 
       local_buf3 => buf3 
   end subroutine

   subroutine free_buffers()
       if (associated(local_buf1)) then
           deallocate(local_buf1)
           nullify(local_buf1)
       end if
       if (associated(local_buf2)) then
           deallocate(local_buf2)
           nullify(local_buf2)
       end if
       if (associated(local_buf3)) then
           deallocate(local_buf3)
           nullify(local_buf3)
       end if
   end subroutine

   function create_diagonal(A)
      real(kind=8), pointer :: create_diagonal(:)
      type(intermediates) :: A
      integer :: N
      integer :: nbuf3, ndimt1, ndimt2

      ndimt1 = size(A%t1)
      ndimt2 = size(A%t2)
      nbuf3 = size(local_buf3)

      N = ndimt1 + ndimt2
      allocate(create_diagonal(N))
      create_diagonal = 1.0d0

      select case (A%eomtype)
         case ("IP")
         case ("EA")
         case ("EE")
             create_diagonal = 1.0d0
      end select

   end function

   function create_sigma_vector(A, r)
      real(kind=8), pointer :: create_sigma_vector(:,:)
      type(intermediates) :: A
      real(kind=8), intent(in) :: r(:,:)
      integer :: L, N
      integer :: nbuf3, ndimt1, ndimt2

      N = size(r,1)
      L = size(r,2)
      allocate(create_sigma_vector(N,L))

      nbuf3 = size(local_buf3)
      ndimt1 = size(A%t1)
      ndimt2 = size(A%t2)

      select case (A%eomtype)
         case ("IP")
         case ("EA")
         case ("EE")
            call sigma_ee_1h1p(r(1:ndimt1,:),r(ndimt1+1:N,:),A, &
                               create_sigma_vector(1:ndimt1,:), &
                               create_sigma_vector(ndimt1+1:N,:),&
                               local_buf2,&
                               local_buf3,&
                               nbuf3)
      end select

   end function

   subroutine IM_container(eomtype,A)

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
type(intermediates), intent(inout) :: A

!Local Variables
!----------------
real(8), allocatable :: foo(:),fvo(:),fvv(:) 
integer :: nbuf3

!-----------------------------------------------------------------------------

     nbuf3 = int(size(local_buf3),4)

     allocate (foo(nfoo*rcw))
     allocate (fvo(nfvo*rcw))
     allocate (fvv(nfvv*rcw))

     A%eomtype = eomtype

     call fmtofile (.false.,fvo,foo,fvv) ! read one-body integrals from file.

     call f_bar_mi(foo,fvo,A%t1,A%t2,A%Fbar_mi)

     call f_bar_me(fvo,A%t1,A%Fbar_me)

     call f_bar_ae(fvv,fvo,A%t1,A%t2,A%Fbar_ae)

     call w_ijmn(A%t1,A%t2,A%W_ijmn)

     call wbar_mbej(A%t2,A%Wbar_mbej)

     call w_mbej(A%t1,A%t2,A%W_mbej)

     call w_iemn(A%t1,A%t2,A%W_iemn,A%Wbar_mbej,A%W_ijmn,A%Fbar_me)

     if (eomtype /= "IP") then

      call w_efam(A%t1,A%t2,A%Wbar_mbej,A%Fbar_me,A%W_efam,local_buf2,local_buf3)

      call w_efam_diagram2 (A%t1,A%t2,local_buf1,local_buf2,local_buf3,nbuf3,A%W_efam)

     endif 

     call w_mnie(A%t1,A%W_mnie)

     if (eomtype /= "IP") then

      call w_amef(A%t1,A%W_amef)

     endif 

   deallocate(foo)
   deallocate(fvo)
   deallocate(fVV)


   end subroutine

   subroutine sigma_ee_1h1p(r1,r2,B,sigma1,sigma2,buf2,buf3,nbuf3)
#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Description--------------------------------------------
!     Calculates 1h1p & 2h2p type sigma vectors for excitation energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),target   :: r1(:,:), r2(:,:) ! R-Amplitudes
  integer, intent(in)     :: nbuf3    
  real(8), intent(in)     :: buf2(*)           
  real(8), intent(inout)  :: buf3(*)           
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout)           :: sigma1(:,:), sigma2(:,:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), pointer     :: r1_ptr(:),r2_ptr(:)
  real(8), allocatable :: int_1b(:),sigma2_temp(:) 
  integer :: off1,off2,irp,istart,mint,m,n,k,nint
  integer :: off3,nd,i,j
  logical :: done 
!--------------------------------------------------

  r1_ptr(1:ndimt1*rcw)     => r1(1:ndimt1*rcw,1)
  r2_ptr(1:ndimt2*rcw)     => r2(1:ndimt2*rcw,1)

! class : (HBar_SS*R)(a,i)
!=========================

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),B%Fbar_ae,r1,sigma1,1.0d0,1.0d0,nrep)

    write(*,*)'inside sigma 1', dot_product(sigma1(ndimt1*rcw,:),sigma1(ndimt1*rcw,:))

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

    write(*,*)'dot product fbar_pp', dot_product(B%Fbar_mi(1:nfoo),B%Fbar_mi(1:nfoo))
    call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),B%Fbar_mi,r1,sigma1,-1.0d0,1.0d0,nrep)

    write(*,*)'inside sigma 2', dot_product(sigma1(ndimt1*rcw,:),sigma1(ndimt1*rcw,:))

!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

    call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),r1,sigma1,1.0d0,1.0d0,nrep,B%W_mbej)


    write(*,*)'inside sigma 3', dot_product(sigma1(ndimt1*rcw,:),sigma1(ndimt1*rcw,:))

! class : [HBar_SD*R](a,i)
!=========================

!S_ai = S_ai + F_bar_me(m,e)*r2(ae,im)
!--------------------------------------

    call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), &
   &   (/"p1","o1"/),B%Fbar_me,sigma1,1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    write(*,*)'inside sigma 4', dot_product(sigma1(ndimt1*rcw,:),sigma1(ndimt1*rcw,:))

!S_ai = S_ai + w_amef(am,ef) * r2(ef,im)
!----------------------------------------
 
    call contraction_442((/"p1","o2","p2","p3"/),(/"p2","p3","o1","o2"/),(/"p1","o1"/),sigma1,1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_amef,RightTensor=r2_ptr)


    write(*,*)'inside sigma 5', dot_product(sigma1(ndimt1*rcw,:),sigma1(ndimt1*rcw,:))
!S_ai = S_ai - W_mnie(mn,ie)*r2(ae,mn)
!---------------------------------------

    call contraction_442((/"o2","o3","o1","p2"/),(/"p1","p2","o2","o3"/),(/"p1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_mnie,RightTensor=r2_ptr)

    write(*,*)'inside sigma 6', dot_product(sigma1(ndimt1*rcw,:),sigma1(ndimt1*rcw,:))

! class : [HBar_DS*R](ab,ij)
!===========================

!S(ab,ij) = S(ab,ij) + P(ab) w_iemn(ma,ij)*r1(b,m)
!--------------------------------------------------

    call contraction_424 ((/"o3","p1","o1","o2"/),(/"p2","o3"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%W_iemn, &
   &    RightTensor=r1_ptr) 


    write(*,*)'inside sigma 7', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))

!S(ab,ij) = S(ab,ij) + P(ij) w_efam(ab,ej)*r1(e,i)
!--------------------------------------------------

    call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","o1"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%W_efam, &
   &    RightTensor=r1_ptr) 


    write(*,*)'inside sigma 8', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))

!int_1b(b,f) = int_1b(b,f) + w_amef (bm,fe)*r1(e,m)
!S(ab,ij) = S(ab,ij) + P(ab) int_1b(b,f)*t2(af,ij)
!---------------------------------------------------

    allocate(int_1b(nfvv*rcw))


    int_1b = 0.0d0

    call contraction_422 ((/"p2","o3","p4","p3"/),(/"p3","o3"/),(/"p2","p4"/),r1_ptr,int_1b,1.0d0,1.0d0,nrep,B%W_amef)

    call contraction_244 ((/"p2","p4"/),(/"p1","p4","o1","o2"/),(/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0, &
   &                       nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)


    write(*,*)'inside sigma 8', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))
!int_1b(n,j) = int_1b(n,j) + w_mnie(nm,je)*r1(e,m)
!S(ab,ij) = S(ab,ij) + int_1b(n,j) * t2(ab,in)
!---------------------------------------------------

    allocate(int_1b(nfoo*rcw))

    int_1b = 0.0d0

    call contraction_422 ((/"o4","o3","o2","p3"/),(/"p3","o3"/),(/"o4","o2"/),r1,int_1b,1.0d0,1.0d0,nrep,LeftTensor=B%W_mnie)

    call contraction_244 ((/"o4","o2"/),(/"p1","p2","o1","o4"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)


    write(*,*)'inside sigma 9', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))


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


    write(*,*)'inside sigma 10', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))

! S(ab,ij) = S(ab,ij) + w_vvvv(ab,ef)*r2(ef,ij)
!------------------------------------------------

     allocate(int_1b(ndimt2*rcw))

     call gettau (B%t1,B%t2,int_1b)
      off1 = 1
      off2 = 1
      do irp = 1, nrep
         istart = 0
         done =.false.
         if (nvvt(irp).eq.0) cycle
   10    call w_efab_out_of_core (irp,istart,done,int_1b,buf2,B%t1,nbuf3,buf3,nint)
         m = nvvt(irp)
         n = noot(irp)
         k = nint
         call xgemm ('n','n',m,n,k,a1,buf3,m,r2_ptr(off1+istart*rcw),m, &
     &               a1,sigma2(off2,1),m)
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off1 = off1 + m * n * rcw
         off2 = off2 + m * n * rcw
      enddo

      deallocate(int_1b)

! S(ab,ij) = S(ab,ij) + w_ijmn(mn,ij) * r2(ab,mn)
!--------------------------------------------------

      call contraction_444((/"o3","o4","o1","o2"/),(/"p1","p2","o3","o4"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_ijmn)


! S(ab,ij) = S(ab,ij) + P(ab)P(ij) w_mbej(mb,ej)*r2(ae,im)
!----------------------------------------------------------

      call contraction_444((/"o3","p2","p3","o2"/),(/"p1","p3","o1","o3"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_mbej)

!int_1b(a,f) = int_1b(a,f) - V(nm,fe)*r2(ea,mn)  
!S(ab,ij) = S(ab,ij) + int_1b(a,f)*t2(fb,ij)  
!-----------------------------------------------

    allocate(int_1b(nfvv*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p3","p1","o3","o4"/),(/"p1","p4"/),int_1b,-1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    call contraction_244 ((/"p1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)

    write(*,*)'inside sigma 12', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))

!int_1b(i,n) = int_1b(i,n) + V(nm,fe)*r2(fe,im)
!S(ab,ij) = S(ab,ij) - int_1b(i,n)*t(ab,jn)
!------------------------------------------------

    allocate(int_1b(nfoo*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p4","p3","o1","o3"/),(/"o1","o4"/),int_1b,1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    call contraction_244 ((/"o1","o4"/),(/"p1","p2","o2","o4"/), &
    & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)
    write(*,*)'inside sigma 13', dot_product(sigma2(ndimt2*rcw,:),sigma2(ndimt2*rcw,:))

  end subroutine

   subroutine sigma_IP(nvx,r1,r2,B,sigma1,sigma2,ndimr1,ndimr2)
#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Description--------------------------------------------
!     Calculates 2h1p and 1h type sigma vectors for ionization energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),target   :: r1(:,:), r2(:,:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout)           :: sigma1(:,:),sigma2(:,:)  ! sigma-vectors
  integer, intent(in) :: ndimr1, ndimr2,nvx(nrep)
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), pointer     :: r1_ptr(:),r2_ptr(:)
  real(8), allocatable :: int_1b(:) 
  real(8)              :: ddot
  integer :: i
!--------------------------------------------------
  r1_ptr => r1(:,1)
  r2_ptr => r2(:,1)
!------------------------------------------
! S(a*,i) = S(a*,i) - Fbar_mi(m,i)*r(a*,m)
!------------------------------------------

  call contraction_222 ((/"o2","o1"/),(/"c1","o2"/),(/"c1","o1"/),B%Fbar_mi,r1,sigma1,-1.0d0,1.0d0,nrep)

!  write(*,*)'diagram 1',dot_product(sigma1(ndimr1*rcw,:),sigma1(ndimr1*rcw,:))

  write(*,*)'diagram 1 ', ddot(ndimr1,sigma1(:,1:1),1,sigma1(:,1:1),1)

!--------------------------------------------
! S(a*,i) = S(a*,i) + Fbar_me(m,e)*r2(a*e,im)
!--------------------------------------------

    call contraction_242 ((/"o2","p2"/),(/"c1","p2","o1","o2"/), &
   &   (/"c1","o1"/),B%Fbar_me,sigma1,1.0d0,1.0d0,nrep,RightTensor=r2_ptr)


  write(*,*)'diagram 2 ', ddot(ndimr1,sigma1(:,1:1),1,sigma1(:,1:1),1)
!-------------------------------------
!S(a*,i) = S(a*,i) - W_mnie(mn,ie)*r2(a*e,mn)
!-------------------------------------

    call contraction_442((/"o2","o3","o1","p2"/),(/"c1","p2","o2","o3"/),(/"c1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_mnie,RightTensor=r2_ptr)


  write(*,*)'diagram 3 ', ddot(ndimr1,sigma1(:,1:1),1,sigma1(:,1:1),1)

!--------------------------------------------
!  S(a*b,ij) = S(a*b,ij)-W(mb,ij)*r(a*,m)
!--------------------------------------------

    call contraction_424 ((/"o3","p2","o1","o2"/),(/"c1","o3"/), &
   &    (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=B%W_iemn, &
   &    RightTensor=r1_ptr) 


  write(*,*)'diagram 4 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!----------------------------------------------------
! S(a*b,ij) = S(a*b,ij) - P(ij) f_bar_mi(m,j)*r2(a*b,im) 
!-----------------------------------------------------

    call contraction_244 ((/"o3","o2"/),(/"c1","p2","o1","o3"/), &
   & (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_mi,RightTensor=r2_ptr)


  write(*,*)'diagram 5 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!---------------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + P(ij) w_mbej(mb,ej)*r2(a*e,im)
!----------------------------------------------------------

      call contraction_444((/"o3","p2","p3","o2"/),(/"c1","p3","o1","o3"/), &
     &    (/"c1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_mbej)


  write(*,*)'diagram 6 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!--------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + w_ijmn(mn,ij) * r2(a*b,mn)
!--------------------------------------------------

      call contraction_444((/"o3","o4","o1","o2"/),(/"c1","p2","o3","o4"/), &
     &    (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_ijmn)


  write(*,*)'diagram 7 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!----------------------------------------------------
! S(a*b,ij) = S(a*b,ij) - P(ij) f_bar_mi(m,j)*r2(a*b,im) 
!-----------------------------------------------------

    call contraction_244 ((/"o3","o2"/),(/"c1","p2","o1","o3"/), &
   & (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_mi,RightTensor=r2_ptr)

  write(*,*)'diagram 8 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)
!------------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + f_bar_ae(b,e) * r2(a*e,ij)
!------------------------------------------------------

    call contraction_244 ((/"p2","p3"/),(/"c1","p3","o1","o2"/), &
   & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_ae,RightTensor=r2_ptr)


  write(*,*)'diagram 9 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!int_1b(a*,f) = int_1b(a*,f) - V(nm,fe)*r2(ea*,mn)  
!S(a*b,ij) = S(a*b,ij) + int_1b(a*,f)*t2(fb,ij)  
!-----------------------------------------------

    allocate(int_1b(nvx(1)*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"c1","p3","o3","o4"/),(/"c1","p4"/),int_1b,-1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    call contraction_244 ((/"c1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)

  write(*,*)'diagram 10 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

   end subroutine

   end module 
