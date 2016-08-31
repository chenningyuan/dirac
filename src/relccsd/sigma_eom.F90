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
    public create_eom_ndet
    public create_eom_nvx
    public initialize_buffers
    public free_buffers
    public free_intermediates
    public free_eom_ndet
    public free_eom_nvx
    public get_matrix_size
    public set_eom_ndet
    public set_eom_nvx
    public get_eom_nvx
    public trial_vectors_ee
    public trial_vectors_ip

    private
    public intermediates
    type intermediates
      character(8) :: eomtype
      real(8),allocatable,dimension(:) :: Fbar_mi, Fbar_me, Fbar_ae                                         ! all 1-b intermediates 
      real(8),allocatable,dimension(:) :: W_ijmn, W_ejmb, W_mbej, W_iemn, W_efam, W_amef, W_mnie, Wbar_mbej ! all 2-b intermediates   
      real(8),allocatable,dimension(:) :: t1, t2
      real(8),allocatable,dimension(:) :: H_II
    endtype intermediates

    real(kind=8), pointer :: local_buf1(:) 
    real(kind=8), pointer :: local_buf2(:) 
    real(kind=8), pointer :: local_buf3(:) 
    integer, pointer :: eom_ndet(:)
    integer, pointer :: eom_nvx(:)

   contains


   subroutine create_eom_nvx(nrep)
       integer, intent(in) :: nrep

       if (.not.associated(eom_nvx)) allocate(eom_nvx(nrep))
       eom_nvx = 0
   end subroutine

   subroutine set_eom_nvx(item,value)
       integer, intent(in) :: item, value

       eom_nvx(item) = value
   end subroutine

   function get_eom_nvx(item)
       integer :: get_eom_nvx
       integer, intent(in) :: item

       get_eom_nvx = eom_nvx(item)  
   end function

   subroutine free_eom_nvx()
       deallocate(eom_nvx)
       nullify(eom_nvx)
   end subroutine

   subroutine create_eom_ndet(max_excitation_level)
       integer, intent(in) :: max_excitation_level

       if (.not.associated(eom_ndet)) allocate(eom_ndet(max_excitation_level))
       eom_ndet = 0
   end subroutine

   subroutine free_eom_ndet()
       if (associated(eom_ndet)) deallocate(eom_ndet)
       nullify(eom_ndet)
   end subroutine

   subroutine set_eom_ndet(excitation_level, nr_det)
       integer, intent(in) :: excitation_level
       integer, intent(in) :: nr_det
       eom_ndet(excitation_level) = nr_det
   end subroutine 

   function get_matrix_size(A)
      integer :: get_matrix_size
      type(intermediates), intent(in) :: A

      get_matrix_size = eom_ndet(1) + eom_ndet(2)
   end function

   subroutine free_intermediates(A)
      type(intermediates), intent(inout) :: A

      if (allocated(A%Fbar_mi)) deallocate(A%Fbar_mi)
      if (allocated(A%Fbar_me)) deallocate(A%Fbar_me)
      if (allocated(A%Fbar_ae)) deallocate(A%Fbar_ae)
      if (allocated(A%W_ijmn)) deallocate(A%W_ijmn)
      if (allocated(A%W_ejmb)) deallocate(A%W_ejmb)
      if (allocated(A%W_mbej)) deallocate(A%W_mbej)
      if (allocated(A%W_iemn)) deallocate(A%W_iemn)
      if (allocated(A%W_efam)) deallocate(A%W_efam)
      if (allocated(A%W_amef)) deallocate(A%W_amef)
      if (allocated(A%W_mnie)) deallocate(A%W_mnie)
      if (allocated(A%Wbar_mbej)) deallocate(A%Wbar_mbej)
      if (allocated(A%t1)) deallocate(A%t1)
      if (allocated(A%t2)) deallocate(A%t2)
      if (allocated(A%H_II)) deallocate(A%H_II)

   end subroutine

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

      ndimt1 = eom_ndet(1)
      ndimt2 = eom_ndet(2)
      nbuf3 = size(local_buf3)

      N = ndimt1 + ndimt2

      allocate(create_diagonal(N))
      create_diagonal = 1.0d0

      select case (A%eomtype)
         case ("IP")
             create_diagonal = A%H_II 
         case ("EA")
         case ("EE")
             create_diagonal = A%H_II 
      end select

   end function

   function create_sigma_vector(A, r, old_sigma_vector)
      real(kind=8), pointer :: create_sigma_vector(:,:)
      real(kind=8), pointer :: old_sigma_vector(:,:)
      type(intermediates) :: A
      real(kind=8), intent(in) :: r(:,:)
      integer :: L, N, L_already_done 
      integer :: nbuf3, ndimr1, ndimr2
      integer :: start_r1, end_r1, start_r2, end_r2
      integer :: iroot


      nbuf3 = size(local_buf3)

      N = size(r,1)
      L = size(r,2)

      allocate(create_sigma_vector(N,L))

      if (associated(old_sigma_vector)) then 
          L_already_done = size(old_sigma_vector,2) + 1
          do iroot = 1, L_already_done
              create_sigma_vector(:,iroot) = old_sigma_vector(:,iroot)
          end do
          deallocate(old_sigma_vector)
          nullify(old_sigma_vector)
      else
          L_already_done = 1
      endif

      ndimr1 = eom_ndet(1)
      ndimr2 = eom_ndet(2)

      start_r1 = 1
      end_r1   = ndimr1
      start_r2 = end_r1 + 1
      end_r2   = end_r1 + ndimr2

      select case (A%eomtype)
         case ("IP")

             print *, "debug: eom_nvx:",eom_nvx
             do iroot = L_already_done, L 

                create_sigma_vector(:,iroot:iroot) = 0.0d0
                 print "(2x,a,i6,a,i6)", "debug, constructing EOM-IP sigma vector ",iroot," out of ",L
                 call sigma_IP(eom_nvx,  &
                 &             r(start_r1:end_r1, iroot:iroot), &
                 &             r(start_r2:end_r2, iroot:iroot), &
                 &             A,                               &
                 &             create_sigma_vector(start_r1:end_r1, iroot:iroot), &
                 &             create_sigma_vector(start_r2:end_r2, iroot:iroot), &
                 &             ndimr1, &
                 &             ndimr2) 
             enddo

         case ("EA")
         case ("EE")

             do iroot = L_already_done, L

                create_sigma_vector(:,iroot:iroot) = 0.0d0
                 print "(2x,a,i6,a,i6)", "debug, constructing EOM-EE sigma vector ",iroot," out of ",L
                 call sigma_ee_1h1p(r(start_r1:end_r1, iroot:iroot), &
                 &                  r(start_r2:end_r2, iroot:iroot), & 
                 &                  A,                               &
                 &                  create_sigma_vector(start_r1:end_r1, iroot:iroot), &
                 &                  create_sigma_vector(start_r2:end_r2, iroot:iroot), &
                 &                  local_buf2, &
                 &                  local_buf3, &
                 &                  int(nbuf3,4))
             enddo
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

     call w_ejmb(A%t1,A%t2,A%W_ejmb)

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
  real(8), intent(in),contiguous,target   :: r1(:,:), r2(:,:) ! R-Amplitudes
  integer, intent(in)     :: nbuf3    
  real(8), intent(in)     :: buf2(*)           
  real(8), intent(inout)  :: buf3(*)           
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout)           :: sigma1(:,:), sigma2(:,:)  ! sigma-vectors
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), pointer     :: r1_ptr(:),r2_ptr(:)
  real(8), allocatable :: int_1b(:),sigma2_temp(:),ResL2_temp(:)
  integer :: off1,off2,irp,istart,mint,m,n,k,nint
  integer :: off3,nd,i,j
  logical :: done 
  real(8)              :: ddot
!--------------------------------------------------

  r1_ptr(1:ndimt1*rcw)  => r1
  r2_ptr(1:ndimt2*rcw)  => r2

! class : (HBar_SS*R)(a,i)
!=========================

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),B%Fbar_ae,r1,sigma1,1.0d0,1.0d0,nrep)

    write(*,*)'diagram 1 ', ddot(ndimt1,sigma1(:,1:1),1,sigma1(:,1:1),1)

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

!    call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),B%Fbar_mi,r1,sigma1,-1.0d0,1.0d0,nrep)


    call contraction_222 ((/"p1","o2"/),(/"o2","o1"/),(/"p1","o1"/),r1,B%Fbar_mi,sigma1,-1.0d0,1.0d0,nrep)

    write(*,*)'diagram 2 ', ddot(ndimt1,sigma1(:,1:1),1,sigma1(:,1:1),1)

!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

    call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),r1_ptr,sigma1,1.0d0,1.0d0,nrep,B%W_mbej)

    write(*,*)'diagram 3 ', ddot(ndimt1,sigma1(:,1:1),1,sigma1(:,1:1),1)

! class : [HBar_SD*R](a,i)
!=========================

!S_ai = S_ai + F_bar_me(m,e)*r2(ae,im)
!--------------------------------------

    call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), &
   &   (/"p1","o1"/),B%Fbar_me,sigma1(:,1:1),1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    write(*,*)'diagram 4 ', ddot(ndimt1,sigma1(:,1:1),1,sigma1(:,1:1),1)

!S_ai = S_ai + w_amef(am,ef) * r2(ef,im)
!----------------------------------------
 
    call contraction_442((/"p1","o2","p2","p3"/),(/"p2","p3","o1","o2"/),(/"p1","o1"/),sigma1,1.0d0,1.0d0,nrep, &
   &         LeftTensor=B%W_amef,RightTensor=r2_ptr)

    write(*,*)'diagram 5 ', ddot(ndimt1,sigma1(:,1:1),1,sigma1(:,1:1),1)
!S_ai = S_ai - W_mnie(mn,ie)*r2(ae,mn)
!---------------------------------------

!   call contraction_442((/"o2","o3","o1","p2"/),(/"p1","p2","o2","o3"/),(/"p1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
!  &         LeftTensor=B%W_mnie,RightTensor=r2_ptr)


     call contraction_442((/"p1","p2","o2","o3"/),(/"o2","o3","o1","p2"/),(/"p1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
    &         RightTensor=B%W_mnie,LeftTensor=r2_ptr)


    write(*,*)'diagram 6 ', ddot(ndimt1,sigma1(:,1:1),1,sigma1(:,1:1),1)

! class : [HBar_DS*R](ab,ij)
!===========================

!S(ab,ij) = S(ab,ij) - P(ab) w_iemn(mb,ij)*r1(a,m)
!--------------------------------------------------

!   call contraction_424 ((/"o3","p2","o1","o2"/),(/"p1","o3"/), &
!  &    (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=B%W_iemn, &
!  &    RightTensor=r1_ptr) 


    call contraction_244 ((/"p1","o3"/),(/"o3","p2","o1","o2"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%W_iemn, &
   &   LeftTensor=r1_ptr) 


    write(*,*)'diagram 7 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!S(ab,ij) = S(ab,ij) + P(ij) w_efam(ab,ej)*r1(e,i)
!--------------------------------------------------

!   call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","o1"/), &
!  &    (/"p1","p2","o2","o1"/),sigma2,1.0d0,-1.0d0,nrep,LeftTensor=B%W_efam, &
!  &    RightTensor=r1_ptr) 



    call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","o1"/), &
   &    (/"p1","p2","o2","o1"/),sigma2,1.0d0,-1.0d0,nrep, &
   &    RightTensor=r1_ptr) 




    write(*,*)'diagram 8 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!int_1b(b,f) = int_1b(b,f) + w_amef (bm,fe)*r1(e,m)
!S(ab,ij) = S(ab,ij) + P(ab) int_1b(b,f)*t2(af,ij)
!---------------------------------------------------

    allocate(int_1b(nfvv*rcw))


    int_1b = 0.0d0

    call contraction_422 ((/"p2","o3","p4","p3"/),(/"p3","o3"/),(/"p2","p4"/),r1_ptr,int_1b,1.0d0,1.0d0,nrep,B%W_amef)

!   call contraction_244 ((/"p2","p4"/),(/"p1","p4","o1","o2"/),(/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0, &
!  &                       nrep,LeftTensor=int_1b,RightTensor=B%t2)


    call contraction_424 ((/"p1","p4","o1","o2"/),(/"p2","p4"/),(/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0, &
   &                       nrep,RightTensor=int_1b,LeftTensor=B%t2)

    deallocate(int_1b)

    write(*,*)'diagram 9 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!int_1b(n,j) = int_1b(n,j) + w_mnie(nm,je)*r1(e,m)
!S(ab,ij) = S(ab,ij) - int_1b(n,j) * t2(ab,in)
!---------------------------------------------------

    allocate(int_1b(nfoo*rcw))

    int_1b = 0.0d0

    call contraction_422 ((/"o4","o3","o2","p3"/),(/"p3","o3"/),(/"o4","o2"/),r1_ptr,int_1b,1.0d0,1.0d0,nrep,LeftTensor=B%W_mnie)

!   call contraction_244 ((/"o4","o2"/),(/"p1","p2","o1","o4"/), &
!  & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)


    call contraction_424 ((/"p1","p2","o1","o4"/),(/"o4","o2"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=int_1b,LeftTensor=B%t2)


    deallocate(int_1b)

    write(*,*)'diagram 10 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

! class : [HBar_DD*R](ab,ij) 
!===========================

! S(ab,ij) = S(ab,ij) + P(ab) f_bar_ae(b,e) * r2(ae,ij)
!-------------------------------------------------------

!   call contraction_244 ((/"p2","p3"/),(/"p1","p3","o1","o2"/), &
!  & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_ae,RightTensor=r2_ptr)


    call contraction_424 ((/"p1","p3","o1","o2"/),(/"p2","p3"/), &
   & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=B%Fbar_ae,LeftTensor=r2_ptr)


    write(*,*)'diagram 11 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)
     
! S(ab,ij) = S(ab,ij) - P(ij) f_bar_mi(m,j)*r2(ab,im) 
!-----------------------------------------------------

!   call contraction_244 ((/"o3","o2"/),(/"p1","p2","o1","o3"/), &
!  & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=B%Fbar_mi,RightTensor=r2_ptr)


    call contraction_424 ((/"p1","p2","o1","o3"/),(/"o3","o2"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%Fbar_mi,LeftTensor=r2_ptr)



    write(*,*)'diagram 12 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

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

    write(*,*)'diagram 13 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)


! S(ab,ij) = S(ab,ij) + w_ijmn(mn,ij) * r2(ab,mn)
!--------------------------------------------------

      call contraction_444((/"p1","p2","o3","o4"/),(/"o3","o4","o1","o2"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=r2_ptr,RightTensor=B%W_ijmn)

    write(*,*)'diagram 14 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

! S(ab,ij) = S(ab,ij) + P(ab)P(ij) w_mbej(mb,ej)*r2(ae,im)
!----------------------------------------------------------

      call contraction_444((/"p1","p3","o1","o3"/),(/"o3","p2","p3","o2"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,LeftTensor=r2_ptr,RightTensor=B%W_mbej)


!      call contraction_444((/"o3","p2","p3","o2"/),(/"p1","p3","o1","o3"/), &
!     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_mbej)


!     call contraction_444((/"p1","p3","o1","o3"/),(/"p2","o3","o2","p3"/), &
!    &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,LeftTensor=r2_ptr,RightTensor=B%W_ejmb)


!     call contraction_444((/"p2","o3","o2","p3"/),(/"p1","p3","o1","o3"/), &
!    &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,0.0d0,nrep,RightTensor=r2_ptr,LeftTensor=B%W_ejmb)


    write(*,*)'diagram 15 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!int_1b(a,f) = int_1b(a,f) - V(nm,fe)*r2(ea,mn)  
!S(ab,ij) = S(ab,ij) + int_1b(a,f)*t2(fb,ij)  
!-----------------------------------------------

    allocate(int_1b(nfvv*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p1","p3","o4","o3"/),(/"p1","p4"/),int_1b,-1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

    call contraction_244 ((/"p1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)

    write(*,*)'diagram 16 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!int_1b(n,i) = int_1b(n,i) + V(nm,fe)*r2(fe,im)
!S(ab,ij) = S(ab,ij) - int_1b(n,i)*t(ba,jn)
!------------------------------------------------

    allocate(int_1b(nfoo*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p4","p3","o1","o3"/),(/"o4","o1"/),int_1b,1.0d0,1.0d0,nrep,RightTensor=r2_ptr)

!   call contraction_244 ((/"o4","o1"/),(/"p1","p2","o2","o4"/), &
!   & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)


    call contraction_424 ((/"p1","p2","o2","o4"/),(/"o4","o1"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=int_1b,LeftTensor=B%t2)

    deallocate(int_1b)

    write(*,*)'diagram 17 ', ddot(ndimt2,sigma2(:,1:1),1,sigma2(:,1:1),1)

  end subroutine

   subroutine sigma_IP(nvx,r1,r2,B,sigma1,sigma2,ndimr1,ndimr2)

#include "param.inc"
#include "symm.inc"
#include "complex.inc"

!---------------Description--------------------------------------------
!     Calculates 2h1p and 1h type sigma vectors for ionization energy
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------
  real(8), intent(in),contiguous,target   :: r1(:,:), r2(:,:) ! R-Amplitudes
  type(intermediates), intent(in)  :: B      
  real(8), intent(inout),contiguous :: sigma1(:,:),sigma2(:,:)  ! sigma-vectors
  integer, intent(in) :: ndimr1, ndimr2,nvx(nrep)
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
  real(8), pointer,contiguous :: r1_ptr(:),r2_ptr(:)
  real(8), allocatable :: int_1b(:) 
  real(8)              :: ddot
!--------------------------------------------------
  r1_ptr(1:ndimr1) => r1
  r2_ptr(1:ndimr2) => r2
!------------------------------------------
! S(a*,i) = S(a*,i) - Fbar_mi(m,i)*r(a*,m)
!------------------------------------------

  call contraction_222 ((/"c1","o2"/),(/"o2","o1"/),(/"c1","o1"/),r1_ptr,B%Fbar_mi,sigma1,-1.0d0,1.0d0,nrep)

  write(*,*)'diagram 1 ', ddot(ndimr1,sigma1,1,sigma1,1)
!--------------------------------------------
! S(a*,i) = S(a*,i) + Fbar_me(m,e)*r2(a*e,im)
!--------------------------------------------

    call contraction_422 ((/"c1","p2","o1","o2"/),(/"o2","p2"/),(/"c1","o1"/),B%Fbar_me, &
  &                      sigma1,1.0d0,1.0d0,nrep,r2_ptr)

  write(*,*)'diagram 2 ', ddot(ndimr1,sigma1(:,1:1),1,sigma1(:,1:1),1)
!-------------------------------------
!S(a*,i) = S(a*,i) - W_mnie(mn,ie)*r2(a*e,mn)
!-------------------------------------

    call contraction_442((/"c1","p2","o2","o3"/),(/"o2","o3","o1","p2"/),(/"c1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
   &         RightTensor=B%W_mnie,LeftTensor=r2_ptr)

  write(*,*)'diagram 3 ', ddot(ndimr1,sigma1(:,1:1),1,sigma1(:,1:1),1)

!--------------------------------------------
!  S(a*b,ij) = S(a*b,ij)-W(mb,ij)*r(a*,m)
!--------------------------------------------

    call contraction_244 ((/"c1","o3"/),(/"o3","p2","o1","o2"/), &
   &    (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%W_iemn, &
   &    LeftTensor=r1_ptr) 

  write(*,*)'diagram 4 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!----------------------------------------------------
! S(a*b,ij) = S(a*b,ij) - P(ij) f_bar_mi(m,j)*r2(a*b,im) 
!-----------------------------------------------------

    call contraction_424 ((/"c1","p2","o1","o3"/),(/"o3","o2"/), &
   & (/"c1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,RightTensor=B%Fbar_mi,LeftTensor=r2_ptr)

  write(*,*)'diagram 5 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!---------------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + P(ij) w_mbej(mb,ej)*r2(a*e,im)
!----------------------------------------------------------

      call contraction_444((/"c1","p3","o3","o1"/),(/"o3","p2","p3","o2"/), &
     &    (/"c1","p2","o1","o2"/),sigma2,-1.0d0,0.0d0,nrep,LeftTensor=r2_ptr,RightTensor=B%W_mbej)

  write(*,*)'diagram 6 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!--------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + w_ijmn(mn,ij) * r2(a*b,mn)
!--------------------------------------------------

      call contraction_444((/"c1","p2","o3","o4"/), (/"o3","o4","o1","o2"/),&
     &    (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=r2_ptr,RightTensor=B%W_ijmn)

  write(*,*)'diagram 7 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!------------------------------------------------------
! S(a*b,ij) = S(a*b,ij) + f_bar_ae(b,e) * r2(a*e,ij)
!------------------------------------------------------

    call contraction_424 ((/"c1","p3","o1","o2"/),(/"p2","p3"/), &
   & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=B%Fbar_ae,LeftTensor=r2_ptr)

  write(*,*)'diagram 8 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

!int_1b(a*,f) = int_1b(a*,f) - V(nm,fe)*r2(a*e,nm)  
!S(a*b,ij) = S(a*b,ij) + int_1b(a*,f)*t2(fb,ij)  
!-----------------------------------------------

    allocate(int_1b(nvx(1)*rcw))
    int_1b = 0.0d0

    call contraction_442 ((/"c1","p3","o4","o3"/),(/"o4","o3","p4","p3"/),(/"c1","p4"/),int_1b,-1.0d0,1.0d0,nrep,LeftTensor=r2_ptr)

    call contraction_244 ((/"c1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"c1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=B%t2)

    deallocate(int_1b)

  write(*,*)'diagram 9 ',ddot(ndimr2,sigma2(:,1:1),1,sigma2(:,1:1),1)

 end subroutine

   subroutine trial_vectors_ee (r1,Fbar_mi,Fbar_ae,W_mbej)

#include "symm.inc"
!---------------Description--------------------------------------------
!     generating R1 trial vectors by diagonalizing HBar_SS matrix.
!-----------------------------------------------------------------------

!---------------Calling variables
  real(8), intent(inout)   :: r1(:,:)       ! R-Amplitudes
  real(8), intent(in)  :: Fbar_ae(:),Fbar_mi(:),W_mbej(:)      

!---------------Local ------------------

  real(8), allocatable :: sigma1(:,:),ev(:),evec(:,:) 
  integer              :: i,iroot, nroot,jroot

!--------------------------------------------------

    nroot = size(r1,1)

    allocate(sigma1(nroot,nroot))

    sigma1 = 0.0d0

  do iroot = 1, nroot
! class : (HBar_SS*R)(a,i)
!=========================

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/),Fbar_ae,r1(:,iroot:iroot),sigma1(:,iroot:iroot), &
  &                    1.0d0,1.0d0,nrep)

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

    call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),Fbar_mi,r1(:,iroot:iroot), &
  &                      sigma1(:,iroot:iroot),-1.0d0,1.0d0,nrep)


!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

    call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),r1(:,iroot:iroot), &
  &                      sigma1(:,iroot:iroot),1.0d0,1.0d0,nrep,W_mbej)

  enddo

       allocate (ev(nroot))
       allocate (evec(nroot,nroot))

       evec = 0.0d0
       ev = 0.0d0
       call solve_cg(nroot,sigma1,evec,ev)
       write(*,*)'H eigenvalues', (ev(i),i=1,min(10,nroot))

     deallocate(sigma1)
! copy the eigenvectors to r1
         r1 = evec

       deallocate(ev)
       deallocate(evec)

 end subroutine

   subroutine trial_vectors_ip (r1,Fbar_mi)

#include "symm.inc"
!---------------Description--------------------------------------------
!     generating R1 trial vectors by diagonalizing HBar_SS matrix.
!-----------------------------------------------------------------------

!---------------Calling variables--------------------------------------

  real(8), intent(inout)   :: r1(:,:) ! R-Amplitudes
  real(8), intent(in)  :: Fbar_mi(:)      

!---------------Local Variables--------------------

  real(8), allocatable :: sigma1(:,:),ev(:),evec(:,:) 
  integer              :: i,iroot,nroot,jroot

!--------------------------------------------------

    nroot = size(r1,1)

    allocate(sigma1(nroot,nroot))

    sigma1 = 0.0d0

   do iroot = 1, nroot

     call contraction_222 ((/"o2","o1"/),(/"c1","o2"/),(/"c1","o1"/),Fbar_mi,r1(:,iroot:iroot),sigma1(:,iroot:iroot), &
    &                       -1.0d0,1.0d0,nrep)

   enddo

       allocate (ev(nroot))
       allocate (evec(nroot,nroot))

       evec = 0.0d0
       ev = 0.0d0
       call solve_cg(nroot,sigma1,evec,ev)
         write(*,*)'H eigenvalues', (ev(i),i=1,nroot)

     deallocate(sigma1)
! copy the eigenvectors to r1
         r1 = evec

       deallocate(ev)
       deallocate(evec)

 end subroutine

   end module 
