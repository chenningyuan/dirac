  module eom_density
  use contraction
  implicit none

  public density_IP

    private

    character(8) :: eomtype
    integer :: nroot,state_sym

   contains

    subroutine density_IP (l1,l2,r1,r2,t1,t2,dmo)

  use symmetry_offset
#include "symm.inc"
#include "complex.inc"
!----------------calling variables-------------------------
   real*8, intent(in),contiguous,target :: l1(:,:),l2(:,:),r1(:,:),r2(:,:)
   real*8, intent(in) :: t1(:),t2(:)
   real*8, intent(out) :: dmo(*)
!----------------local variables---------------------------
   real*8, allocatable :: doo(:),dvv(:),dvo(:),dov(:)  
   real*8, allocatable :: work1(:),work2(:) 
   real*8, allocatable :: woo_temp(:),wpp_temp(:),wop_temp(:)
   real(8), pointer,contiguous :: l1_ptr(:),l2_ptr(:),r1_ptr(:),r2_ptr(:)
   type(Offset)         :: e
   integer              :: ndiml1,ndiml2
!-----------------------------------------------------------

  ndiml1 = size(l1,1) 
  ndiml2 = size(l2,1) 

  l1_ptr(1:ndiml1) => l1
  l2_ptr(1:ndiml2) => l2
  r1_ptr(1:ndiml1) => r1
  r2_ptr(1:ndiml2) => r2

    allocate(doo(nfoo*rcw))
    allocate(dvv(nfvv*rcw))
    allocate(dvo(nfvo*rcw))
    allocate(dov(nfvo*rcw))

     doo = 0.0d0
     dvv = 0.0d0
     dvo = 0.0d0
     dov = 0.0d0

!-------------------------------------------------------------------
!   doo(j,i) = doo(j,i) - l2(jm,e*f) * r2(e*f,im)
!-------------------------------------------------------------------

  call contraction_442 ((/"o2","o3","c1","p2"/),(/"c1","p2","o1","o3"/),(/"o2","o1"/),doo,-1.0d0,1.0d0,nrep, &
 &  LeftTensor=l2_ptr,RightTensor=r2_ptr)

   write(*,*)'diagram1'

!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l1(j,e*)*r1(e*,i)
!---------------------------------------------------------------------

  call contraction_222 ((/"o2","c1"/),(/"c1","o1"/),(/"o2","o1"/),l1,r1,doo,-1.0d0,1.0d0,nrep)

   write(*,*)'diagram2'
!---------------------------------------------------------------------
!    doo(j,i) = doo(j,i) - l2(im,a*e)*r1(a*,m)*t1(e,j) 

!    wop_temp(i,e) = wop_temp(i,e) + l2(im,a*e)*r1(a*,m) 
!    doo(j,i) = doo(j,i) + wop_temp(i,e) * t1(e,j)
!---------------------------------------------------------------------

  allocate(wop_temp(nfvo*rcw)) 

  call contraction_422 ((/"o1","o3","c1","p2"/),(/"c1","o3"/),(/"o1","p2"/),r1,wop_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)

  call contraction_222 ((/"o1","p2"/),(/"p2","o2"/),(/"o1","o2"/),wop_temp,t1,doo,1.0d0,1.0d0,nrep)

 deallocate(wop_temp)

   write(*,*)'diagram3'
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,ea)*t2(eb,mn) 
!----------------------------------------------------------------------

! call contraction_442 ((/"o1","o2","c1","p1"/),(/"c1","p2","o1","o2"/),(/"p2","p1"/),dvv,1.0d0,1.0d0,nrep, &
!&     LeftTensor=l2_ptr,RightTensor=t2 )

   write(*,*)'diagram4'
!----------------------------------------------------------------------
!     dvv(b,a) = dvv(b,a) + l2(mn,cb)*r1(c,n)*t1(a,m)
!     wop_temp(m,b) = wop_temp(m,b) - l2(mn,cb)*r1(c,n) 
!     dvv(b,a) = dvv(b,a) + wop_temp(m,b)*t1(a,m)
!----------------------------------------------------------------------

  allocate(wop_temp(nfvo*rcw))

  call contraction_422 ((/"o1","o2","c1","p2"/),(/"c1","o2"/),(/"o1","p2"/),r1,wop_temp,-1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)

  call contraction_222 ((/"o1","p2"/),(/"p1","o1"/),(/"p2","p1"/),wop_temp,t1,dvv,-1.0d0,1.0d0,nrep)

  deallocate(wop_temp)

   write(*,*)'diagram5'
!----------------------------------------------------------------------------
!     dvo(a,i) = t1(a,i)
!----------------------------------------------------------------------------

      call xcopy (nfvo,t1,1,dvo,1)

   write(*,*)'diagram6'
!--------------------------------------------------------------------------
!     dvo(a,i) = dvo(a,i) + l1(m,b) * r(ba,im)
!--------------------------------------------------------------------------

!     call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), 
!    & (/"p1","o1"/),l1,dvo,1.0d0,1.0d0,nrep,RightTensor=t2)

      call contraction_422 ((/"c1","p1","o1","o2"/),(/"o2","c1"/), &
     & (/"p1","o1"/),l1_ptr,dvo,-1.0d0,1.0d0,nrep,LeftTensor=r2_ptr)       

   write(*,*)'diagram6'
!----------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l1(m,e) * r(e,i)*t(a,m)
!----------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_222 ((/"o2","c1"/),(/"c1","o1"/),(/"o2","o1"/), &
     &              l1_ptr,r1_ptr,woo_temp,1.0d0,1.0d0,nrep)                   


      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/), &
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)                 

      deallocate(woo_temp)

   write(*,*)'diagram7'
!-------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*r2(ef,in)*t1(a,m) 
!-------------------------------------------------------------

      allocate(woo_temp(nfoo*rcw))

      woo_temp = 0.0d0

      call contraction_442((/"o2","o3","c1","p3"/),(/"c1","p3","o1","o3" &
     &       /),(/"o2","o1"/),woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr,   &
     &          RightTensor=r2_ptr)

      call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/),&
     &              woo_temp,t1,dvo,-1.0d0,1.0d0,nrep)                

      deallocate(woo_temp)

   write(*,*)'diagram8'
!-----------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*r1(e,i)*t2(af,mn)
!-----------------------------------------------------------------------


    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,nv,ncont,.false.,.false.)

      allocate(wpp_temp(e%oneNonDirac(1)*rcw))

      wpp_temp = 0.0d0

      call contraction_442((/"o2","o3","c1","p3"/),(/"p1","p3","o2","o3" &
     &       /),(/"p1","c1"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=l2_ptr, & 
     &          RightTensor=t2)

      call contraction_222 ((/"p1","c1"/),(/"c1","o1"/),(/"p1","o1"/), &
     &              wpp_temp,r1_ptr,dvo,1.0d0,1.0d0,nrep)                  

      deallocate(wpp_temp)
  call dealloc_array(e)

   write(*,*)'diagram9'
!------------------------------------------------------------------
!      dvo(a,i) = dvo(a,i) - l2(mn,ef)*t1(f,i)*r2(ea,mn)
!------------------------------------------------------------------

      allocate(wpp_temp(nfvv*rcw))

      wpp_temp = 0.0d0

!     call contraction_442((/"c1","p1","o2","o3"/),(/"o2","o3","c1","p3"/), &
!    &       (/"p1","p3"/),wpp_temp,-1.0d0,1.0d0,nrep,LeftTensor=r2_ptr, & 
!    &          RightTensor=l2_ptr)

!     call contraction_222 ((/"p1","p3"/),(/"p3","o1"/),(/"p1","o1"/), &
!    &              wpp_temp,t1,dvo,1.0d0,1.0d0,nrep)                  

      deallocate(wpp_temp)


   write(*,*)'diagram10'
!--------------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) +  l2(mn,ef)*t2(af,in)*r1(e,m) 
!--------------------------------------------------------------------

      allocate(woo_temp(nfvo*rcw))

      woo_temp = 0.0d0

      call contraction_422((/"o2","o3","c1","p3"/),(/"c1","o2"/), &
     &       (/"o3","p3"/),r1_ptr,woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)   

      call contraction_422 ((/"p1","p3","o1","o3"/),(/"o3","p3"/),(/"p1","o1"/),&
     &              t1,dvo,1.0d0,1.0d0,nrep,woo_temp)                

      deallocate(woo_temp)

   write(*,*)'diagram11'
!--------------------------------------------------------------------
!       dvo(a,i) = dvo(a,i) -  l2(mn,ef)*t1(f,i)*t1(a,n)*r1(e,m) 
!--------------------------------------------------------------------

      allocate(woo_temp(nfvo*rcw))

      woo_temp = 0.0d0

      call contraction_422((/"o2","o3","c1","p3"/),(/"c1","o2"/), &
     &       (/"o3","p3"/),r1_ptr,woo_temp,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)   

     allocate(work1(nfoo*rcw)) 

      call contraction_222 ((/"p3","o1"/),(/"o3","p3"/),(/"o3","o1"/),&
     &              t1,woo_temp,work1,1.0d0,1.0d0,nrep)                

      deallocate(woo_temp)

      call contraction_222 ((/"p1","o3"/),(/"o3","o1"/),(/"p1","o1"/),&
     &              t1,work1,dvo,-1.0d0,1.0d0,nrep)                

     deallocate(work1)

   write(*,*)'diagram12'
!----------------------------------------------------------
!       dov(i,a) = dov(i,a) + l2(im,ba) * r1(b,m)
!----------------------------------------------------------

      call contraction_422 ((/"o1","o2","c1","p1"/),(/"c1","o2"/), & 
     & (/"o1","p1"/),r1_ptr,dov,1.0d0,1.0d0,nrep,LeftTensor=l2_ptr)

   write(*,*)'diagram13'
      if (carith) then
         call conjuga (nfoo,doo,1)
         call conjuga (nfvv,dvv,1)
         call conjuga (nfvo,dvo,1)
         call conjuga (nfvo,dov,1)
      endif


     allocate(work1(nfvv*rcw))
     allocate(work2(nfvv*rcw))

      call convdm(work1,work2,doo,dvv,dov,dvo,dmo)

     deallocate(work1,work2)  

    deallocate(dov)
    deallocate(dvo)
    deallocate(doo)
    deallocate(dvv)

      end subroutine


   end module 
