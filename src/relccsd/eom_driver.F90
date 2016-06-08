  module eom_driver
  use sigma_eom
  use relcc_cfg
  implicit none
!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

!    public eom_ip  
!    public eom_ea
    public eom_ee

    private

    character(8) :: eomtype
    integer :: nroot,iarep

   contains

  subroutine read_input()

    iarep = nr_roots(1)
    nroot = nr_roots(2)

    if (relcc_do_eomip) eomtype = 'IP'
    if (relcc_do_eomea) eomtype = 'EA'
    if (relcc_do_eomee) eomtype = 'EE'

  end subroutine

   subroutine input_processor(Fbar_HH,Fbar_PP,Vbar_SS,R1,R2)

#include "symm.inc"
#include "complex.inc"
!description
!===========
!1)  this subroutine processes the guess vectors for eomcc.
!2)  it also defines the denominator for the error vectors, which is relevant in
!    the Davidson step.

!calling variables
!================

   real*8, intent(in)  :: Fbar_HH(:),Fbar_PP(:),Vbar_SS(:)
   real*8, intent(out) :: R1(:,:), R2(:,:) ! (ndimt1,nroot) and (ndimt2,nroot)

!Local Variables
!==============

   integer :: i,ii,irep,j,iaai,aoff,a,aa,ai,arep,iroot,k(1)
   real*8  :: fac1
   real*8,allocatable  :: H_II(:)
   real*8  :: eps((iv(nrep+1)+io(nrep+1))) 
   integer :: icount, acount

!   iarep = nr_roots(1)
!   nroot = nr_roots(2)


!     CALL PRSYMB(LUPRI,'-',75,0)
!     WRITE(LUPRI,'(A,8(/10X,A,I2,3A,I5))')
!    &   ' * No. of excitations evaluated in each boson symmetry:',
!    &   (' boson symmetry',J,' (',REP(J-1),') :',KEXCNV(J),J=1,NBSYM)

!initialize the R-vector
!=======================
 
    R1(1:ndimt1*rcw,1:nroot) = 0.0d0 
    R2(1:ndimt2*rcw,1:nroot) = 0.0d0 

! 1> Determine the minimum energy determinants for IP, EA & EE. The number of them is equal to the number of roots. 

     i = 1
     ii = 1
     do irep = 1, nrep
        call dcopy (no(irep),Fbar_HH(ii),((no(irep)+1)*rcw),eps(i),1)
        i  =  i + no(irep)
        ii =  ii + NO(irep)* no(irep) * rcw
     enddo

     ii = 1
     do irep = 1, nrep
        call dcopy (nv(irep),Fbar_PP(ii),((nv(irep)+1)*rcw),eps(i),1)
        i  =  i + nv(irep)
        ii =  ii + NV(irep)* nv(irep) * rcw
     enddo

   select case (eomtype)

   case ('EE')

    allocate(H_II(nfvo*rcw))

    ai = 1
    ii = 0
    do irep = 1, nrep
       arep = multb (irep,iarep+nrep,2)
        do i = 1, no(irep)
         ii = ii + 1
         fac1 = eps(ii)
         aoff = io(nrep+1) + iv(arep)
          do a = 1, nv(arep)
             aa = aoff + a
             iaai = (iovvo(iarep) + (iivo(arep,irep)+(i-1)*nv(arep)+a-1)*nov(iarep) + &
       &      iiov(irep,arep) + (a-1)*no(irep))*rcw+i ! careful about the complex algebra  
             H_II(ai) = -fac1 + eps(aa) - vbar_ss (iaai)
             ai = ai + 1
          enddo
        enddo
    enddo

    do iroot = 1, nroot
     k = minloc(H_II) 
     write(*,*)'roots',iroot,H_II(k(1)) 
     i = k(1)
     H_II(i) = 8000.0d0
     R1(k,iroot) = 1.0d0
    enddo

     deallocate(H_II)

   case ('IP')

! As of now, I am targetting the lowest energy determinants across the symmetries.


     allocate(H_II(io(nrep+1)*rcw))    

    H_II = eps(1:(io(nrep+1)*rcw))

    do iroot = 1, nroot
     k = minloc(H_II) 
     write(*,*)'roots',iroot,H_II(k(1)) 
     i = k(1)
     H_II(i) = 8000.0d0
     R1(k,iroot) = 1.0d0
    enddo

     deallocate(H_II)

   case ('EA')

     allocate(H_II(nv(iarep)*rcw))    

      H_II = eps((1+no(nrep+1)):(no(nrep+1)+nv(nrep+1)))

      do iroot = 1, nroot
        k = minloc(H_II) 
        write(*,*)'roots',iroot,H_II(k(1)) 
        i = k(1)
        H_II(i) = 8000.0d0
        R1(k,iroot) = 1.0d0
      enddo
     deallocate(H_II)

  end select 


   end subroutine


   subroutine eom_ee(icalc,mxcorr)

   use memory_allocator
   use allocator_parameters, only : klongint

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
!--------------------------------------------------
 integer, intent(in)                :: icalc
 integer(kind=klongint), intent(inout) :: mxcorr
!--------------------------------------------------

   type(intermediates) :: B
   real*8, allocatable :: t1(:),t2(:) 
   real*8, allocatable :: sigma1(:,:),sigma2(:,:) 
   real*8, allocatable :: R1(:,:),R2(:,:) 
   real*8, allocatable :: ev(:),evec(:,:) 
   real*8, allocatable :: buf1(:),buf2(:),buf3(:) ! may be not all of them are necessary. 
!  real*8 :: sub_H(nr_roots(2),nr_roots(2))! sub_H(nroot,nroot)
   real*8, allocatable:: sub_H(:,:) ! sub_H(nroot,nroot)
   integer :: irp,iroot,jroot,i,j,irep
   integer :: ndimr1,ndimr2
   integer :: AllocateStatus

!memory allocation
!------------------

  integer(kind=8)        :: current_mem_use = 0
  integer(kind=8)        :: start_mem_use   = 0
  integer                :: freemem = 0
  integer                :: work_memsize
  integer                :: nbuf1,nbuf2
  integer(kind=8)        :: nbuf3
  integer*8              :: nv5irp,nv6irp,nvirp

    nbuf1 = max(joooo(nrep+1),jvooo(nrep+1),jvvoo(nrep+1),jvovo(nrep+1))
    nbuf2 = nbuf1

!     call alloc( buf1  , nbuf1*rcw             , id="buf1" )
!     call alloc( buf2  , nbuf2*rcw             , id="buf2" )

    call allocator_get_words_inuse(start_mem_use)
    call legacy_lwork_get(work_memsize)

    call allocator_get_maxbuff(freemem,kind(BUF3))

       NBUF3 = freemem / RCW

         NV5IRP = MAX(JVOVV(NREP+1),JOVVV(NREP+1))
         NV6IRP = 1
         DO IRP = 1, NREP
            NV6IRP = MAX(NV6IRP,INT(NVVT(IRP),8)*INT(NVVT(IRP),8))
         ENDDO

!        Check the minimum size for the general use of this array

         NVIRP = MAX(NV5IRP,NV6IRP)

         NBUF3 = MIN(NVIRP,INT(NBUF3,8))
         NBUF3 = MAX(NBUF1,NBUF3)

!     Allocate now also this remaining work array
        allocate(buf3(nbuf3*rcw))
        allocate(buf1(nbuf1*rcw))
        allocate(buf2(nbuf2*rcw))

!!        call alloc( BUF3  , NBUF3*RCW             , id="buf3" )


!     For the moment I am using as much in-core memory as possible. However we
!     may think of writing them in files, replacing the integral ones.  


      allocate (B%Fbar_mi(nfoo*rcw))  
      allocate (B%Fbar_me(nfvo*rcw))  
      allocate (B%Fbar_ae(nfvv*rcw))  

      allocate (B%W_ijmn(nv1*rcw)) 
      allocate (B%W_ejmb(iovvo(nrep+1)*rcw))
      allocate (B%W_mbej(iovvo(nrep+1)*rcw))
      allocate (B%W_iemn(iovoot(nrep+1)*rcw))
      allocate (B%W_efam(nv5*rcw))
      allocate (B%W_amef(nv5*rcw))
      allocate (B%W_mnie(iovoot(nrep+1)*rcw))
      allocate (B%Wbar_mbej(iovvo(nrep+1)*rcw))
       
! allocate the amplitudes

      call read_input()

      allocate (sub_H(nroot,nroot))

      allocate (t1(ndimt1*rcw))
      allocate (t2(ndimt2*rcw))

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use

      if (icalc.eq.1) then
         write (iw,1001) "EOMCC calculation",mxcorr
         goto 451
      endif

    call ccdump

    call getampt(t1,t2)

    call IM_container(eomtype,t1,t2,B,buf1,buf2,buf3,int(nbuf3,4))

     deallocate(B%Wbar_mbej) 

   select case (eomtype)

    case ('EE')

      allocate (r1(ndimt1*rcw,nroot))
      allocate (r2(ndimt2*rcw,nroot))

      allocate (sigma1(ndimt1*rcw,nroot),stat=AllocateStatus)

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      allocate (sigma2(ndimt2*rcw,nroot),stat = AllocateStatus)

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

!  call alloc (sigma1(ndimt1*rcw,nroot))
!  call alloc (sigma2(ndimt2*rcw,nroot))

    call input_processor(B%Fbar_mi,B%Fbar_ae,B%W_mbej,R1,R2)

    do iroot = 1, nroot
      sigma1(1:ndimt1*rcw,iroot:iroot) = 0.0d0 
      sigma2(1:ndimt2*rcw,iroot:iroot) = 0.0d0 

      call sigma_ee_1h1p(r1(1:ndimt1*rcw,iroot:iroot),r2(1:ndimt2*rcw,iroot:iroot),t1,t2,B, &
     &                  sigma1(1:ndimt1*rcw,iroot:iroot),sigma2(1:ndimt2*rcw,iroot:iroot),buf2,buf3,int(nbuf3,4))
    enddo

    case ('IP')

      do irep = 1, nrep
         ndimr2 = ndimr2 + nv(irep)*noot(irep)
      enddo

       allocate (r1(io(nrep+1)*rcw,nroot))
       allocate (r2(ndimr2*rcw,nroot))

       allocate (sigma1(io(nrep+1)*rcw,nroot))
       allocate (sigma2(ndimr2*rcw,nroot))


    case('EA')

      do irep = 1, nrep
         ndimr2 = ndimr2 + nvvt(irep)*no(irep)
      enddo

       allocate (r1(iv(nrep+1)*rcw,nroot))
       allocate (r2(ndimr2*rcw,nroot))

       allocate (sigma1(iv(nrep+1)*rcw,nroot))
       allocate (sigma2(ndimr2*rcw,nroot))


   end select 

    sub_H = 0.0d0

    do jroot = 1, nroot
     do iroot = 1, nroot
      sub_H(iroot,jroot) = dot_product(r1(:,iroot),sigma1(:,jroot)) + dot_product(r2(:,iroot),sigma2(:,jroot)) 
     enddo
    enddo

     write(*,*)'eigenvalues', (sub_H(i,i),i=1,nroot)

     allocate (ev(nroot))
     allocate (evec(nroot,nroot))

    call solve_cg(nr_roots(2),sub_H,evec,ev)

     write(*,*)'eigenvalues', (ev(i),i=1,nroot)
     write(*,*)'eigenvectors', ((evec(i,j),i=1,nroot),j=1,nroot)

     deallocate(ev)
     deallocate(evec)

!    call davidson 

     deallocate (t1,t2,sigma1,sigma2,r1,r2)
 451  continue

!    call dealloc(buf1)  
!    call dealloc(buf2)  
!    call dealloc(buf3)  
     deallocate(buf3)
     deallocate(buf1)
     deallocate(buf2)

     deallocate (B%Fbar_mi)  
     deallocate (B%Fbar_me)  
     deallocate (B%Fbar_ae)  

     deallocate (B%W_ijmn) 
     deallocate (B%W_ejmb)
     deallocate (B%W_mbej)
     deallocate (B%W_iemn)
     deallocate (B%W_efam)
     deallocate (B%W_amef)
     deallocate (B%W_mnie)
     deallocate (sub_H)

1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")

   end subroutine

   end module 
