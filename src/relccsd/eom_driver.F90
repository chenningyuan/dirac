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
    integer :: nroot,state_sym

   contains

  subroutine read_input()
#include "symm.inc"
#include "dcbham.h"

    state_sym = nr_roots(1)
    nroot = nr_roots(2)

    if (relcc_do_eomip) eomtype = 'IP'
    if (relcc_do_eomea) eomtype = 'EA'
    if (relcc_do_eomee) eomtype = 'EE'

!Assigning Dimension to the continuous orbital
!===============================================

     ncont(1:nrep) = 0

     ncont(state_sym) = 1
    if (spinfr) then
       ncont(state_sym+(nrep/4)) = 1 
    else
       ncont(state_sym+(nrep/2)) = 1 
    endif

  end subroutine

   subroutine input_processor(Fbar_HH,Fbar_PP,Vbar_SS,R1,R2,H_II)

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
   real*8, intent(inout) :: R1(:,:), R2(:,:) ! (ndimt1,nroot) and (ndimt2,nroot)
   real*8, intent(inout),target :: H_II(:)

!Local Variables
!==============

   integer :: i,ii,irep,j,iaai,aoff,a,aa,ai,arep,iroot,k(1)
   real*8  :: fac1
   real*8  :: eps((iv(nrep+1)+io(nrep+1))) 
   real*8, pointer :: diag_H(:)
   integer :: icount, acount

   integer :: abij, ijrp, abrp, jrp, jj, irp, ioff, imin, bb, amin, arp, brp, b
   real(kind=8) :: fac, fac3, fac2
   integer :: ndimr1, ndimr2, search_dimension, xroot
   real*8,allocatable :: R1_dum(:,:) 


!initialize the R-vector
!=======================

    ndimr1 = size(R1,1)
    ndimr2 = size(R2,1)
    print *, "dimensions of R1, R2", ndimr1, ndimr2
    R1(1:ndimr1,1:nroot) = 0.0d0 
    R2(1:ndimr2,1:nroot) = 0.0d0 

! 1> Determine the minimum energy determinants for IP, EA & EE. The number of them is equal to the number of roots. 

     i = 1
     ii = 1
     do irep = 1, nrep
        if (no(irep) .gt. 0) call dcopy (no(irep),Fbar_HH(ii),((no(irep)+1)*rcw),eps(i),1)
        i  =  i + no(irep)
        ii =  ii + NO(irep)* no(irep) * rcw
     enddo

     ii = 1
     do irep = 1, nrep
        if (nv(irep) .gt. 0) call dcopy (nv(irep),Fbar_PP(ii),((nv(irep)+1)*rcw),eps(i),1)
        i  =  i + nv(irep)
        ii =  ii + NV(irep)* nv(irep) * rcw
     enddo

   select case (eomtype)

   case ('EE')

    ai = 1
    ii = 0
    do irep = 1, nrep
       arep = multb (irep,state_sym+nrep,1)
        do i = 1, no(irep)
         ii = ii + 1
         fac1 = eps(ii)
         aoff = io(nrep+1) + iv(arep)
          do a = 1, nv(arep)
             aa = aoff + a
             iaai = (iovvo(state_sym) + (iivo(arep,irep)+(i-1)*nv(arep)+a-1)*nov(state_sym) + &
       &      iiov(irep,arep) + (a-1)*no(irep))*rcw+i ! careful about the complex algebra  
             H_II(ai) = -fac1 + eps(aa) - vbar_ss (iaai)
             ai = ai + 1
          enddo
        enddo
    enddo

!aspg+as
      ABIJ = ai 
      DO IJRP = 1, NREP
       abrp = multb (ijrp+nrep,state_sym+nrep,1)
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + 1
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
!aspg+as

! New trial vectors for R1 
!if (nroot.eq.(ndimr1+ndimr2)) then
!    do iroot = 1, ndimr1
!        R1(iroot,iroot) = 1.0d0
!    end do
!    do iroot = 1, ndimr2
!        R2(iroot,ndimr1+iroot) = 1.0d0
!    end do
!else
!  allocate (R1_dum(ndimr1,ndimr1))
!    R1_dum = 0.0d0

!    do iroot = 1, ndimr1
!        R1_dum(iroot,iroot) = 1.0d0
!    end do

!     call trial_vectors_ee(R1_dum,Fbar_HH,Fbar_PP,Vbar_SS)
!     
!     R1 = R1_dum(:,1:nroot) 

!   deallocate(R1_dum)
!end if

     allocate(diag_H(ndimr1))
     diag_H = H_II(1:ndimr1)

       do iroot = 1, nroot
        k = minloc(diag_H) 
        write(*,*)'roots',iroot,H_II(k(1)),k(1) 
        i = k(1)
        diag_H(i) = 8000.0d0
        print *, "debug: i, iroot, ndimr1, ndimr2,size(H_II)",i,iroot,ndimr1,ndimr2,size(H_II)
        R1(i,iroot) = 1.0d0
       enddo
     deallocate(diag_H)

   case ('IP')

!As of now, I am targetting the lowest energy determinants across the symmetries.

    ai = 1
    ii = 0
    do irep = 1, nrep
       arep = multb (irep,1+nrep,1)
        do i = 1, no(irep)
         ii = ii + 1
         fac1 = eps(ii)
          do a = 1, ncont(arep)
             H_II(ai) = -fac1 
             ai = ai + 1
          enddo
        enddo
    enddo

!      H_II = -eps(1:(ndimr1*rcw)) !+ 100.0d0

      ABIJ = ai 
      DO IJRP = 1, NREP
       abrp = multb (ijrp+nrep,1+nrep,1)
      DO 30 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 30
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO 40 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               DO A = 1, NCONT(ARP)
                  FAC = FAC3! - 100.0d0
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + 1
               ENDDO
            ENDDO
 40         CONTINUE
         ENDDO
      ENDDO
 30   CONTINUE
      ENDDO

!  allocate (R1_dum(ndimr1,ndimr1))
!    R1_dum = 0.0d0

!    do iroot = 1, ndimr1
!        R1_dum(iroot,iroot) = 1.0d0
!    end do

!     call trial_vectors_ip(R1_dum,Fbar_HH)
!     
!     R1 = R1_dum(:,1:nroot) 

!   deallocate(R1_dum)

 !  if (nroot.eq.(ndimr1+ndimr2)) then
 !   do iroot = 1, ndimr1
 !       R1(iroot,iroot) = 1.0d0
 !   end do
 !   do iroot = 1, ndimr2
 !       R2(iroot,ndimr1+iroot) = 1.0d0
 !   end do
 !  else

     search_dimension = ncont(state_sym)*no(state_sym) 

     allocate(diag_H(ndimr1))
     diag_H = H_II(1:ndimr1)

        xroot = 1 

       do iroot = 1, 2*nroot
        k = minloc(diag_H) 
        write(*,*)'roots',iroot,H_II(k(1)),k(1) 
        i = k(1)
        diag_H(i) = 8000.0d0

        if (i <= search_dimension) then
         R1(i,xroot) = 1.0d0
         xroot = xroot + 1
        endif

       enddo
     deallocate(diag_H)
 ! end if 

   case ('EA')

    ai = 1
    do irep = 1, nrep
       arep = multb (irep,1+nrep,1)
        do i = 1, ncont(irep)
         aoff = io(nrep+1) + iv(arep)
          do a = 1, nv(arep)
             aa = aoff + a
             H_II(ai) =  eps(aa) 
             ai = ai + 1
          enddo
        enddo
    enddo

      ABIJ = ai 
      DO IJRP = 1, NREP
       abrp = multb (ijrp+nrep,1+nrep,1)
      DO 50 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         DO I = 1, NCONT(IRP)
            FAC2 =  FAC1
            DO 60 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 60
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + 1
               ENDDO
            ENDDO
 60         CONTINUE
         ENDDO
      ENDDO
 50   CONTINUE
      ENDDO


! calculate the orbital dimension within which our desired root should belong..

  search_dimension = ncont(state_sym)*nv(state_sym) 

!! allocate (r1_dum(ndimr1,ndimr1))
!!   r1_dum = 0.0d0

!!   do iroot = 1, ndimr1
!!       r1_dum(iroot,iroot) = 1.0d0
!!   end do

!!    call trial_vectors_EA(R1_dum,Fbar_PP)
!!    
!!    r1 = r1_dum(:,1:nroot) 

!!  deallocate(R1_dum)

 !  if (nroot.eq.(ndimr1+ndimr2)) then
 !   do iroot = 1, ndimr1
 !       R1(iroot,iroot) = 1.0d0
 !   end do
 !   do iroot = 1, ndimr2
 !       R2(iroot,ndimr1+iroot) = 1.0d0
 !   end do
 !  else

     allocate(diag_H(ndimr1))
     diag_H = H_II((no(nrep+1)+1):(ndimr1+no(nrep+1)))

        xroot = 1

       do iroot = 1, 2*nroot
        k = minloc(diag_H) 
        write(*,*)'roots',iroot,H_II(k(1)),k(1) 
        i = k(1)
        diag_H(i) = 8000.0d0
        if (i <= search_dimension) then
         R1(i,xroot) = 1.0d0
         xroot = xroot + 1
        endif
       enddo
     deallocate(diag_H)
!! end if 

  end select 

   end subroutine

   subroutine eom_ee(icalc,mxcorr,eps)
   use davidson
   use memory_allocator
   use allocator_parameters, only : klongint

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
!--------------------------------------------------
 integer, intent(in)                   :: icalc
 integer(kind=klongint), intent(inout) :: mxcorr
 real*8                                :: eps(*)
!--------------------------------------------------

   type(intermediates) :: B
   real*8, allocatable :: sigma1(:,:),sigma2(:,:),sigma(:,:) 
   real*8, allocatable :: R1(:,:),R2(:,:)
   real*8, pointer :: r(:,:) 
   real*8, allocatable :: ev(:),evec(:,:) 
   real*8, allocatable, target :: buf1(:),buf2(:),buf3(:) ! may be not all of them are necessary. 
   real*8, allocatable:: sub_H(:,:) ! sub_H(nroot,nroot)
   integer :: irp,iroot,jroot,i,j,irep,jrep,ijrep,nvx(nrep),nox(nrep)
   integer :: ndimr1,ndimr2
   integer :: AllocateStatus

   type(davidson_results) :: davResults
   logical :: left_vectors, right_vectors
   logical :: symmetric = .false.
   logical :: verbose = .false.
   integer :: N, N1, N2 ! size of determinantal space 

!memory allocation
!------------------

  integer(kind=8)        :: current_mem_use = 0
  integer(kind=8)        :: start_mem_use   = 0
  integer                :: freemem = 0
  integer                :: work_memsize
  integer                :: nbuf1,nbuf2
  integer                :: nbuf3
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

!!        call alloc( BUF3  , NBUF3*RCW             , id="buf3" )

!     For the moment I am using as much in-core memory as possible. However we
!     may think of writing them in files, replacing the integral ones...  


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
      allocate (B%t1(ndimt1*rcw))
      allocate (B%t2(ndimt2*rcw))
       
! allocate the amplitudes

      call read_input()

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use

      if (icalc.eq.1) then
         write (iw,1001) "EOMCC calculation",mxcorr
         goto 451
      endif

        allocate(buf3(nbuf3*rcw))
        allocate(buf1(nbuf1*rcw))
        allocate(buf2(nbuf2*rcw))

    call initialize_buffers(buf1, buf2, buf3)

    call ccdump

    call getampt(B%t1,B%t2)

    call IM_container(eomtype,B)

     call free_buffers(buf1,buf2,buf3)

    call create_eom_ndet(2)

   select case (eomtype)

    case ('EE')

      ndimr1 = ndimt1
      ndimr2 = ndimt2

      call set_eom_ndet(1, ndimr1)
      call set_eom_ndet(2, ndimr2)

      N1 = ndimr1*rcw
      N2 = ndimr2*rcw
      N = ndimr1 + ndimr2

!     nroot = N

      allocate (r (N*rcw,nroot), stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      allocate (B%H_II(N*rcw),  stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      call input_processor(B%Fbar_mi, &
      &                    B%Fbar_ae, &
      &                    B%W_mbej,  &
      &                    r(1:ndimr1*rcw,:), &
      &                    r((ndimr1*rcw)+1:N*rcw,:),  &
      &                    B%H_II)

      if (nroot .eq. N) then
! this is to test the full diagonalizaion

          allocate (sigma(N*rcw,nroot),stat = AllocateStatus)
          IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

          do iroot = 1, nroot
              print *, "debug, constructing sigma vector ",iroot
              sigma(:,iroot:iroot) = 0.0d0
              call sigma_ee_1h1p(r(1:ndimr1*rcw,iroot:iroot),         &
              &                  r((ndimr1*rcw)+1:N*rcw,iroot:iroot), &
              &                  B,                                   &
              &                  sigma(1:ndimr1*rcw,iroot:iroot),     &
              &                  sigma((ndimr1*rcw)+1:N*rcw,iroot:iroot), &
              &                  nbuf3)
          enddo
      endif

    case ('IP')

       ndimr2 = 0
       ndimr1 = 0
       nvx = 0
       
       call create_eom_nvx(nrep)

      do jrep = 1, nrep
        do irep = 1, nrep
        ijrep = multb(irep,jrep,1)
        nvx(ijrep) = nvx(ijrep) + ncont(irep)*nv(jrep)
        call set_eom_nvx(ijrep, (get_eom_nvx(ijrep) + ncont(irep)*nv(jrep)))
        enddo
      enddo

        do irep = 1, nrep
         ndimr1 = ndimr1 + ncont(irep)*no(irep)
        enddo

      do irep = 1, nrep
          ndimr2 = ndimr2 + get_eom_nvx(irep)*noot(irep)
      enddo

       N = ndimr1 + ndimr2
!       nroot = N

       call set_eom_ndet(1, ndimr1*rcw)
       call set_eom_ndet(2, ndimr2*rcw)

      allocate (r (N*rcw,nroot), stat = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
       allocate (B%H_II(N*rcw)) 
       
       r = 0.0d0

       B%H_II = 1.0d0

      call input_processor(B%Fbar_mi, &
      &                    B%Fbar_ae, &
      &                    B%W_mbej,  &
      &                    r(1:ndimr1*rcw,:), &
      &                    r((ndimr1*rcw)+1:N*rcw,:),  &
      &                    B%H_II)


      if (nroot .eq. N) then
! this is to test the full diagonalizaion
          allocate (sigma(N*rcw,nroot))
          do iroot = 1, nroot
              print *, "debug, constructing sigma vector ",iroot
              sigma(:,iroot:iroot) = 0.0d0
              call sigma_IP(nvx,   &
              &             r(1:ndimr1*rcw,iroot:iroot),         &
              &             r((ndimr1*rcw)+1:N*rcw,iroot:iroot), &
              &             B,                                   &
              &             sigma(1:ndimr1*rcw,iroot:iroot),     &
              &             sigma((ndimr1*rcw)+1:N*rcw,iroot:iroot), &
              &             ndimr1,                              &
              &             ndimr2)
          enddo
       end if


    case('EA')

       ndimr2 = 0
       ndimr1 = 0
       nox(1:nrep) = 0
       
       call create_eom_nvx(nrep)

      do jrep = 1, nrep
        do irep = 1, nrep
        ijrep = multb(irep,jrep,1)
        nox(ijrep) = nox(ijrep) + ncont(irep)*no(jrep)
        call set_eom_nvx(ijrep, (get_eom_nvx(ijrep) + ncont(irep)*no(jrep)))
        enddo
      enddo

        do irep = 1, nrep
         ndimr1 = ndimr1 + ncont(irep)*nv(irep)
        enddo

      do irep = 1, nrep
          ndimr2 = ndimr2 + get_eom_nvx(irep)*nvvt(irep)
      enddo

       N = ndimr1 + ndimr2
!       nroot = N

       call set_eom_ndet(1, ndimr1)
       call set_eom_ndet(2, ndimr2)

       allocate (r(N*rcw,nroot))
       allocate (B%H_II(N*rcw)) 
       
       r = 0.0d0

       B%H_II = 1.0d0

      call input_processor(B%Fbar_mi, &
      &                    B%Fbar_ae, &
      &                    B%W_mbej,  &
      &                    r(1:ndimr1*rcw,:), &
      &                    r((ndimr1*rcw)+1:N*rcw,:),  &
      &                    B%H_II)

!      do iroot = 1, nroot

!         sigma1(1:ndimr1*rcw,iroot:iroot) = 0.0d0 
!         sigma2(1:ndimr2*rcw,iroot:iroot) = 0.0d0 
!         call sigma_IP(nvx,r1(:,iroot:iroot),r2(:,iroot:iroot),B, &
!    &                  sigma1(:,iroot:iroot),sigma2(:,iroot:iroot),ndimr1,ndimr2)

!         sigma(1:N*rcw,iroot:iroot) = 0.0d0 
!         call sigma_IP(nvx, &
!                       r(1:ndimr1*rcw,iroot:iroot), &
!         &             r((ndimr1*rcw)+1:N*rcw,iroot:iroot), &
!         &             B,                            &
!         &             sigma(1:ndimr1*rcw,iroot:iroot), &
!         &             sigma((ndimr1*rcw)+1:N*rcw,iroot:iroot), &
!         &             ndimr1, &
!         &             ndimr2)

!       end do

   end select 

   if (nroot.eq.N) then
! this is to test the full diagonalizaion
      allocate (sub_H(nroot,nroot))
      sub_H = 0.0d0

      do jroot = 1, nroot
         do iroot = 1, nroot
            sub_H(iroot,jroot) = dot_product(r(:,iroot),sigma(:,jroot))
         enddo
      enddo

      call dump_to_file(42,"CCjacobianH.bin",sub_H)

       allocate (ev(nroot))
       allocate (evec(nroot,nroot))

       call solve_cg(nroot,sigma,evec,ev)
       write(*,*)'Hp eigenvalues', (ev(i),i=1,min(20,nroot))

       deallocate(ev)
       deallocate(evec)
       deallocate (sub_H)
       if (allocated(sigma)) deallocate (sigma)

    end if

    if (.true.) then
!       call davidson_setup(1.0e-2, subspace_size, nroot)
        symmetric = .false.
        verbose   = .false.
        left_vectors = .false.
        right_vectors = .true.

        call initialize_davidson_results(davResults, N, nroot, left_vectors, right_vectors)
        call davidson_driver(nroot, B, symmetric, davResults, verbose, r, left_vectors, right_vectors)

!      Fock-Space like output for EA, EE and IP

       call eom_analysis(davResults%eValuesRe,davResults%eVectorsR,eps)

       call free_davidson_results(davResults)
        
    end if
    call free_eom_ndet()
    call free_eom_nvx()

      if (allocated(B%H_II)) deallocate(B%H_II)


! todo:
! allocations/deallocations of intermediates would be more simply handled via
!    call initialize_intermediates() 

 451 continue

     call free_intermediates(B)

1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
   end subroutine

   function merge_arrays(dimension_to_merge, a,b)
       real(kind=8), pointer :: merge_arrays(:,:)
       real(kind=8), intent(in) :: a(:,:), b(:,:)
       integer, intent(in) :: dimension_to_merge

       integer :: Na, Nb, La, Lb, Ns, Ls

       Na = size(a,1)
       La = size(a,2)

       Nb = size(b,1)
       Lb = size(b,2)

       if ((dimension_to_merge .eq. 1) .and. (La .eq. Lb)) then
           Ns = Na + Nb
           Ls = La
           allocate(merge_arrays(Ns,Ls))
           merge_arrays(1:Na, :) = a(1:Na, :) 
           merge_arrays(Na+1:Ns, :) = b(1:Nb, :)
       else if ((dimension_to_merge .eq. 2) .and. (Na .eq. Nb)) then
           Ls = La + Lb
           Ns = Na
           allocate(merge_arrays(Ns,Ls))
           merge_arrays(:, 1:La) = a(:, 1:La) 
           merge_arrays(:, La+1:Ls) = b(:, 1:Lb)
       else
           stop('mismatch in array dimensions')
       end if

   end function

   subroutine dump_to_file(fileno, filename, matrix)
       integer, intent(in) :: fileno
       character(len=14), intent(in) :: filename
       real(kind=8) :: matrix(:,:)
       integer :: irec
  
       inquire (iolength = irec) matrix
       open (unit=fileno, file=filename, form="unformatted", access="direct", recl=irec)
       write(unit=fileno, rec=1) matrix
       close(unit=fileno)
   end subroutine


   subroutine eom_analysis(eigenvalues,eigenvectors,eps) 

! Analyse the R-vectors of EOM-IP,EE,EA. Based on FS_ANALYSIS1 routine by Luuk

#include "param.inc"
#include "inpt.inc"
#include "complex.inc"
#include "symm.inc"
#include "results.inc"
#include "ccpar.inc"
#include "files.inc"

!---------------Calling variables--------------------------------------

     real*8,intent(inout) :: eigenvalues(rcw,1:*),eigenvectors(*)
     real*8,intent(in)    :: eps(*)

!---------------Local variables -----------------------------------
      INTEGER , ALLOCATABLE :: INDX_ALL(:), IRPS_ALL(:), nstate_sym(:,:)  
      INTEGER :: IERR, ISIZE, select_eom,dum2
      REAL*8 :: IND(MXREP),dum1(rcw),eigenvalues_temp(1:rcw,0:nroot)
      REAL*8 :: E_OFF
!---------------Executable code--------------------------------------

! Generate an array denoting number of states per symmetry...

    dum1(1:rcw) = 0.0d0
    dum2 = 0

   allocate(nstate_sym(nrep,3))  

      nstate_sym = 0

   select case(eomtype)

    case('IP')

      nstate_sym(state_sym,3) = nroot
      select_eom = 3

    case('EA')

      nstate_sym(state_sym,2) = nroot

      select_eom = 2

    case('EE')

     select_eom = 4

  end select  

      ISIZE = nroot + 1 ! possibly correct
      ALLOCATE (INDX_ALL(ISIZE),STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in INDX_ALL allocation !')
      ENDIF
      ALLOCATE (IRPS_ALL(ISIZE),STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in IRPS_ALL allocation !')
      ENDIF

!     Initialize the reference energy

      E_OFF = ESCF + ECCSD

!     Print an ordered list of eigenvalues

      IF ((select_eom==2).OR.(select_eom==3)) THEN
!        Odd electron case : fermion irreps
         CALL PRT_EV (IW,NREP,nstate_sym(1,select_eom),REPNA(1), &
     &                E_OFF,ECCSDIM,RCW, &
     &                EIGENVALUES(1:RCW,1:ISIZE),ISIZE, &
     &                INDX_ALL,IRPS_ALL,IND,dum1,dum2)
      ELSE 
!        Add energy of (0,0) sector to list of printed energies for the (1,1) sector
!        This makes the routine print excitation energies
         nstate_sym(1,select_eom) = nstate_sym(1,select_eom) + 1
         EIGENVALUES(1:RCW,0) = 0.D0
!        Even electron case : boson irreps
         CALL PRT_EV (IW,NREP,nstate_sym(1,select_eom),REPNA(NREP+1), &
     &                E_OFF,ECCSDIM,RCW, &
     &                EIGENVALUES(1:RCW,0:ISIZE),ISIZE, &
     &                INDX_ALL,IRPS_ALL,IND,dum1,dum2)
!        One should not forget this line : the real dimension of the first irrep should be retained !
         nstate_sym(1,select_eom) = nstate_sym(1,select_eom) - 1
      ENDIF

      call evector_analysis (select_eom,nstate_sym,eps,iprnt,e_off,eigenvalues(1:rcw,1:isize),eigenvectors)

      DEALLOCATE (INDX_ALL,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in INDX_ALL deallocation !')
      ENDIF
      DEALLOCATE (IRPS_ALL,STAT=IERR)
      IF (IERR.NE.0) THEN
       CALL QUIT('Error in IRPS_ALL deallocation !')
      ENDIF

     deallocate(nstate_sym)

   end subroutine


   subroutine evector_analysis (select_eom,nstate_sym,eps,ipr,e_off,ev_all,eigenvectors)

!---------------Description--------------------------------------------
!
!     Analyze the excited states and calculate transition moments
!           SELECT_EOM : selection of eomtype
!                  EPS : orbital energies
!                  IPR : print level
!                E_OFF : CCSD energy to be added to eigenvalue
!               EV_ALL : eigenvalues for all irreps
!         eigenvectors : R-Vectors for all irreps (packed)
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee from FS_ANALYSIS routine by Luuk.
!
!
!---------------Calling variables--------------------------------------

#include "complex.inc"
#include "symm.inc"

   real*8, intent(in) :: eps(*),e_off
   real*8, intent(in) :: ev_all(rcw,nroot),eigenvectors(rcw,*) 
   integer            :: select_eom,ipr,nstate_sym(nrep,nroot)

!---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "ccpar.inc"

!---------------Local variables--------------------------------------

   integer       :: i_eigenvalue, j_eigenvalue,ndt, irpoff
   character*100 :: label
   real*8        :: value,tres
   integer       :: i,ii,irep,j,aoff,a,aa,ai,arep,ssym
   integer       :: ijrp, abrp, jrp, jj, kk, irp, ioff, imin, bb, amin, arp, brp, b

!---------------Executable code--------------------------------------
 
!     Now print an analysis of the vectors
!     Initialize and get labels for the excited state coefficiens

      if (ipr.le.0) then
         tres = 0.1d0
      elseif (ipr.le.2) then
         tres = 0.01d0
      else
         tres = 0.0d0
      end if

      i_eigenvalue = 0

      write (iw,1000) tres
      do ssym = 1,nrep
 
!     get the dimension of the space
 
         ndt = nstate_sym(ssym,select_eom) 
         if (ndt.eq.0) cycle

!     add an offset for irrep names for bosonic irreps

         if (select_eom.le.3) then
          irpoff = 0
         else
          irpoff = nrep
         endif

!   Loop over the eigenvalues

     ii = 1
     do j_eigenvalue = 1, nstate_sym(ssym,select_eom) 
          i_eigenvalue = i_eigenvalue + 1 
!    I have seen crazy ouputs with large fock space calcs, put a limit
          if (j_eigenvalue.lt.100) then
          write (iw,1010) repna(irpoff+ssym),j_eigenvalue,  &
    &             e_off+ev_all(1,i_eigenvalue),ev_all(1,i_eigenvalue)
!    Loop over the R-vectors and print the interesting ones

    select case(eomtype)

      case('EA')

      do irep=1,nrep
        arep = multb (irep,1+nrep,1)
        do i = 1, ncont(irep)
          aoff = io(nrep+1) + iv(arep)
         do a = 1, nv(arep)
          aa = aoff + a
          write(label,1001) repna(arep),(a+no(arep)),eps(aa)
            value = abs(eigenvectors(1,ii))
            if (carith) value = value + abs(eigenvectors(2,ii))
            if (carith .and. value.gt.tres) then
                write (iw,'(2f15.5,3x,a50)') & 
  &             eigenvectors(1,ii),eigenvectors(2,ii),label
            else if (value.gt.tres) then              
            write (iw,'(f15.5,3x,a50)') eigenvectors(1,ii),label
            end if
          ii = ii + 1
         enddo
        enddo
      enddo

      do irep=1,nrep
      ijrp = irep
      abrp = multb (ijrp+nrep,1+nrep,1)
      do jrp = 1, nrep
      jj = io(jrp)
      irp = multb(jrp,ijrp+nrep,2)
      do j = 1, no(jrp)
         jj = jj + 1
         do i = 1, ncont(irp)
            do brp = 1, nrep
            bb = iv(brp) + io(nrep+1)
            arp = multb(brp,abrp+nrep,2)
            if (arp.lt.brp) cycle
            aoff = iv(arp) + io(nrep+1)
            do b = 1, nv(brp)
               bb = bb + 1
               amin = 1
               if (arp.eq.brp) amin = b + 1
               do a = amin, nv(arp)
                  aa = aoff + a
                  write(label,1002) repna(jrp),j,eps(jj), &
                &          repna(arp),(a+no(arp)),eps(aa),&
                &          repna(brp),(b+no(brp)),eps(bb)

                  value = abs(eigenvectors(1,ii))
                  if (carith) value = value + abs(eigenvectors(2,ii))
                  if (carith .and. value.gt.tres) then
                      write (iw,'(2f15.5,3x,a70)') & 
            &             eigenvectors(1,ii),eigenvectors(2,ii),label
                  else if (value.gt.tres) then              
                  write (iw,'(f15.5,3x,a70)') eigenvectors(1,ii),label
                  end if
                  ii = ii + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo
      case('IP')

    jj = 0
  do irep = 1, nrep
     arep = multb (irep,1+nrep,1)
      do i = 1, no(irep)
       jj = jj + 1
        do a = 1, ncont(arep)
          write(label,1001) repna(irep),i,eps(jj)
            value = abs(eigenvectors(1,ii))
            if (carith) value = value + abs(eigenvectors(2,ii))
            if (carith .and. value.gt.tres) then
                write (iw,'(2f15.5,3x,a50)') & 
  &             eigenvectors(1,ii),eigenvectors(2,ii),label
            else if (value.gt.tres) then              
            write (iw,'(f15.5,3x,a50)') eigenvectors(1,ii),label
            end if
          ii = ii + 1
        enddo
      enddo
  enddo


  do ijrp = 1, nrep
       abrp = multb (ijrp+nrep,1+nrep,1)
    do jrp = 1, nrep
      jj = io(jrp)
      irp = multb(jrp,ijrp+nrep,2)
      if (irp.lt.jrp) cycle
      ioff = io(irp)
      do j = 1, no(jrp)
         jj = jj + 1
         imin = 1
         if (irp.eq.jrp) imin = j + 1
         do i = imin, no(irp)
            kk = ioff + i
            do brp = 1, nrep
            bb = iv(brp) + io(nrep+1)
            arp = multb(brp,abrp+nrep,2)
            aoff = iv(arp) + io(nrep+1)
            do b = 1, nv(brp)
               bb = bb + 1
               do a = 1, ncont(arp)

                  write(label,1003) repna(irp),i,eps(kk), &
                &          repna(jrp),j,eps(jj),&
                &          repna(brp),(b+no(brp)),eps(bb)

                  value = abs(eigenvectors(1,ii))
                  if (carith) value = value + abs(eigenvectors(2,ii))
                  if (carith .and. value.gt.tres) then
                      write (iw,'(2f15.5,3x,a70)') & 
            &             eigenvectors(1,ii),eigenvectors(2,ii),label
                  else if (value.gt.tres) then              
                  write (iw,'(f15.5,3x,a70)') eigenvectors(1,ii),label
                  end if
                  ii = ii + 1
               enddo
            enddo
            enddo
         enddo
      enddo
    enddo
  enddo

      case('EE')

          write(*,*)'under preparation'  

      end select
      endif
         end do
      end do
 
 1000 FORMAT (//' Analysis of eigenvectors of the effective Hamiltonian' &
     &/' First line  :  Energy, Eigenvalue' &
     &/' other lines : Coefficients of contributing determinants',&
     & ' (above a treshold of ',E6.1,')')
 1010 FORMAT (/' Irrep ',A4,' State ',I5,2F16.8)

!for R1 of IP and EA
 1001 format ('| ',a4,' #',i4,' (',f8.3,') |',25x)

!for R2 of EA
 1002 format (1x,a4,' #',i4,' (',f8.3,') -> ', &
     &           a4,' #',i4,' (',f8.3,'), ',   &
     &           a4,' #',i4,' (',f8.3,')  ')


!for R2 of IP
 1003 format (1x,a4,' #',i4,' (',f8.3,'), ',   &
                 a4,' #',i4,' (',f8.3,') -> ', &
     &           a4,' #',i4,' (',f8.3,')  ')

!for R1 of EE
 1004 format (1x,a4,' #',i4,' (',f8.3,') -> ', &
     &             a4,' #',i4,' (',f8.3,')  ')

!for R2 of EE
 1005 format (1x,a4,' #',i4,' (',f8.3,'), ',   &
                 a4,' #',i4,' (',f8.3,') -> ', &
     &           a4,' #',i4,' (',f8.3,'), ',   &
     &           a4,' #',i4,' (',f8.3,')  ')

  end subroutine

   end module 
