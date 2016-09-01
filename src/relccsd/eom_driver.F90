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

   integer :: abij, ijrp, jrp, jj, irp, ioff, imin, bb, amin, arp, brp, b
   real(kind=8) :: fac, fac3, fac2
   integer :: ndimr1, ndimr2
   real*8,allocatable :: R1_dum(:,:) 


!   iarep = nr_roots(1)
!   nroot = nr_roots(2)


!     CALL PRSYMB(LUPRI,'-',75,0)
!     WRITE(LUPRI,'(A,8(/10X,A,I2,3A,I5))')
!    &   ' * No. of excitations evaluated in each boson symmetry:',
!    &   (' boson symmetry',J,' (',REP(J-1),') :',KEXCNV(J),J=1,NBSYM)

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

!   allocate(H_II(nfvo*rcw))

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

!aspg+as
      ABIJ = ai 
      DO IJRP = 1, NREP
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
            ARP = MULTB(BRP,IJRP+NREP,2)
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

if (nroot.eq.(ndimr1+ndimr2)) then
     do iroot = 1, ndimr1
         R1(iroot,iroot) = 1.0d0
     end do
     do iroot = 1, ndimr2
         R2(iroot,ndimr1+iroot) = 1.0d0
     end do
else
   allocate (R1_dum(ndimr1,ndimr1))
     R1_dum = 0.0d0

     do iroot = 1, ndimr1
         R1_dum(iroot,iroot) = 1.0d0
     end do

      call trial_vectors_ee(R1_dum,Fbar_HH,Fbar_PP,Vbar_SS)
      
      R1 = R1_dum(:,1:nroot) 

    deallocate(R1_dum)
end if


   case ('IP')

!As of now, I am targetting the lowest energy determinants across the symmetries.

      H_II = -eps(1:(ndimr1*rcw)) !+ 100.0d0

      ABIJ = ndimr1+1 
      DO IJRP = 1, NREP
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
            ARP = MULTB(BRP,IJRP+NREP,2)
!            IF (ARP.LT.BRP) GOTO 40
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
!              AMIN = 1
!              IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = 1, NCONT(ARP)
!                 AA = AOFF + A
                  FAC = FAC3 !- 100.0d0
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + 1
               ENDDO
            ENDDO
 40         CONTINUE
         ENDDO
      ENDDO
 30   CONTINUE
      ENDDO

   allocate (R1_dum(ndimr1,ndimr1))
     R1_dum = 0.0d0

     do iroot = 1, ndimr1
         R1_dum(iroot,iroot) = 1.0d0
     end do

      call trial_vectors_ip(R1_dum,Fbar_HH)
      
      R1 = R1_dum(:,1:nroot) 

    deallocate(R1_dum)

 !  if (nroot.eq.(ndimr1+ndimr2)) then
 !   do iroot = 1, ndimr1
 !       R1(iroot,iroot) = 1.0d0
 !   end do
 !   do iroot = 1, ndimr2
 !       R2(iroot,ndimr1+iroot) = 1.0d0
 !   end do
 !  else

!!   allocate(diag_H(ndimr1))
!!   diag_H = H_II(1:ndimr1)

!!     do iroot = 1, nroot
!!      k = minloc(diag_H) 
!!      write(*,*)'roots',iroot,H_II(k(1)),k(1) 
!!      i = k(1)
!!      diag_H(i) = 8000.0d0
!!      print *, "debug: i, iroot, ndimr1, ndimr2, size(H_II)",i,iroot,ndimr1,ndimr2,size(H_II)
!!      R1(i,iroot) = 1.0d0
!!     enddo
!!   deallocate(diag_H)
 ! end if 

   case ('EA')

      H_II = eps((1+no(nrep+1)):(no(nrep+1)+nv(nrep+1)))

      do iroot = 1, nroot
        k = minloc(H_II) 
        write(*,*)'roots',iroot,H_II(k(1)) 
        i = k(1)
        H_II(i) = 8000.0d0
        R1(i,iroot) = 1.0d0
      enddo

  end select 


   end subroutine


   subroutine eom_ee(icalc,mxcorr)
   use davidson
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
   real*8, allocatable :: sigma1(:,:),sigma2(:,:),sigma(:,:) 
   real*8, allocatable :: R1(:,:),R2(:,:)
   real*8, pointer :: r(:,:) 
   real*8, allocatable :: ev(:),evec(:,:) 
   real*8, allocatable, target :: buf1(:),buf2(:),buf3(:) ! may be not all of them are necessary. 
   real*8, allocatable:: sub_H(:,:) ! sub_H(nroot,nroot)
   integer :: irp,iroot,jroot,i,j,irep,jrep,ijrep,nvx(nrep)
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
        allocate(buf3(nbuf3*rcw))
        allocate(buf1(nbuf1*rcw))
        allocate(buf2(nbuf2*rcw))

!!        call alloc( BUF3  , NBUF3*RCW             , id="buf3" )
        call initialize_buffers(buf1, buf2, buf3)

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
       
! allocate the amplitudes

      call read_input()

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use

      if (icalc.eq.1) then
         write (iw,1001) "EOMCC calculation",mxcorr
         goto 451
      endif

      allocate (B%t1(ndimt1*rcw))
      allocate (B%t2(ndimt2*rcw))

    call ccdump

    call getampt(B%t1,B%t2)

    call IM_container(eomtype,B)
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
              &                  buf2, &
              &                  buf3, &
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
      print *, "debug, in eom_driver, nvx:",nvx
      do jrep = 1, nrep
         print *, "debug, in eom_driver, get_nvx(jrep):",jrep,get_eom_nvx(jrep)
      end do

        do irep = 1, nrep
         ndimr1 = ndimr1 + ncont(irep)*no(irep)
        enddo

      do irep = 1, nrep
!       ndimr2 = ndimr2 + nvx(irep)*noot(irep)
          ndimr2 = ndimr2 + get_eom_nvx(irep)*noot(irep)
         print *, "debug, ndimr2: irep",irep,"nvx(irep)",nvx(irep),"noot(irep)",noot(irep)
      enddo
       
       print *, "debug, ndimr1",ndimr1,"ndimr2",ndimr2

       N = ndimr1 + ndimr2
!       nroot = N

       call set_eom_ndet(1, ndimr1)
       call set_eom_ndet(2, ndimr2)

       allocate (sigma(N*rcw,nroot))
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

      do irep = 1, nrep
         ndimr2 = ndimr2 + nvvt(irep)*no(irep)
      enddo

       allocate (r1(iv(nrep+1)*rcw,nroot))
       allocate (r2(ndimr2*rcw,nroot))

!      allocate (sigma1(iv(nrep+1)*rcw,nroot))
!      allocate (sigma2(ndimr2*rcw,nroot))

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
!      write(*,*)'matrix to diagonalize', sub_H
      call dump_to_file(42,"CCjacobianH.bin",sub_H)

!      write(*,*)'matrix to diagonalize', sigma 
      call dump_to_file(42,"CCjacobianHp.bin", sigma)

       allocate (ev(nroot))
       allocate (evec(nroot,nroot))

       call solve_cg(nroot,sigma,evec,ev)
       write(*,*)'Hp eigenvalues', (ev(i),i=1,min(20,nroot))

!      evec = 0.0d0
!      ev = 0.0d0
!      call solve_cg(nroot,sub_H,evec,ev)
!      write(*,*)'H eigenvalues', (ev(i),i=1,min(20,nroot))

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
        call davidson_driver(nroot, B, symmetric, davResults, verbose, r)
        print *, "debug, after davidson !"
        call free_davidson_results(davResults)
        print *, "debug, after free_davidson_results !"
        
    end if

! this seems to be the offending code, ucomment
!   if (associated(r)) then
!       deallocate (r)
!       r => NULL()
!   end if

    call free_buffers()
        print *, "debug, after free_buffers !"
    call free_eom_ndet()
        print *, "debug, after free_eom_ndet !"
    call free_eom_nvx()
        print *, "debug, after free_eom_nvx !"

 451 continue
! todo:
! allocations/deallocations of intermediates would be more simply handled via
!    call initialize_intermediates() 
     call free_intermediates(B)
        print *, "debug, after free_intermediates !"

1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
        print *, "debug, before the end of eom_driver !"
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

   end module 
