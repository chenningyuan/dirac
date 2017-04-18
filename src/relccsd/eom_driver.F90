  module eom_driver
  use sigma_eom
  use trial_vectors_eom
  use relcc_cfg
  implicit none
!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

!    public eom_ip  
!    public eom_ea
    public eom_ee
    public number_of_states

    private

    character(8) :: eomtype
    integer :: nroot,state_sym 

   contains

  subroutine set_continuous_orbital(active_irrep)
#include "symm.inc"
#include "dcbham.h"
#include "dgroup.h"
     integer, intent(in) :: active_irrep

!Assigning Dimension to the continuous orbital
!===============================================
     ncont(1:nrep) = 0
     ncont(active_irrep) = 1
     if (spinfr) then
        ncont(active_irrep+(nrep/4)) = 1 
     elseif (linear) then
        ncont(active_irrep+1) = 1 
     else
        ncont(active_irrep+(nrep/2)) = 1 
     endif
  end subroutine

  subroutine verify_input()
#include "symm.inc"
#include "dcbham.h"
#include "dgroup.h"
#include "files.inc"
    logical :: reset_roots_in_symmetry = .false.

    if (relcc_do_eomip) eomtype = 'IP'
    if (relcc_do_eomea) eomtype = 'EA'
    if (relcc_do_eomee) eomtype = 'EE'

    if (    (relcc_do_eomip .and. relcc_do_eomea) &
        .or.(relcc_do_eomip .and. relcc_do_eomee) &
        .or.(relcc_do_eomea .and. relcc_do_eomee)) &
          call quit('Specifying more than one EOM kind of calculation at a time (IP/EE/EA)  not allowed !') 

     do state_sym = 1, nrep 
        if (relcc_eom_nroots(state_sym).ne.0) then
! checking whether or not, for ip/ea, we can actually ask for this many roots 
!
! we choose to restrict the number of roots one cas ask for in the eom-ee/ip/ea calculationsnos to
! avoid setting up guess trial vectors which would have to sample the 2h1p (ip, ea) or 2h2p subblocks:
!
! 1. IP: maximum number of roots for irrepN is reset back to the number of occupied spinors (nocc) in irrepN
! 2. EA: maximum number of roots for irrepN is reset back to the number of virtual  spinors (nvirt), per symmetry 
! 3. EE: maximum number of roots for irrepN is reset back to  nocc(irrepM)*nvirt(irrepO) for irrepN = irrepM * irrepO

          write (iw,'(2x,A,I6,2x,A,2x,I3)') "Requested",relcc_eom_nroots(state_sym),"root(s) for symmetry",state_sym
          select case (eomtype)
             case ('IP')
                if (relcc_eom_nroots(state_sym).gt.no(state_sym)) then
                   relcc_eom_nroots(state_sym) = no(state_sym)
                   reset_roots_in_symmetry = .true.
                else if (relcc_eom_nroots(state_sym).eq.-1) then
                   reset_roots_in_symmetry = .true.
                else
                   reset_roots_in_symmetry = .false.
                end if
             case ('EA')
                if (relcc_eom_nroots(state_sym).gt.nv(state_sym)) then
                    relcc_eom_nroots(state_sym) = nv(state_sym)
                   reset_roots_in_symmetry = .true.
                else if (relcc_eom_nroots(state_sym).eq.-1) then
                   reset_roots_in_symmetry = .true.
                else
                   reset_roots_in_symmetry = .false.
                end if
             case ('EE')
                if (relcc_eom_nroots(state_sym).gt.nvo(state_sym)) then
                    relcc_eom_nroots(state_sym) = nvo(state_sym)
                   reset_roots_in_symmetry = .true.
                else if (relcc_eom_nroots(state_sym).eq.-1) then
                   reset_roots_in_symmetry = .true.
                else
                   reset_roots_in_symmetry = .false.
                end if
          end select

          if (reset_roots_in_symmetry) then
             write (iw,'(4x,A)') "Warning: more EOM-"//eomtype(1:2)//" roots than there are occupied/virtual orbitals/pairs."
             write (iw,'(13x,A,I3,2x,A,2x,I3)') "Number or roots for symmetry",state_sym,"set to",relcc_eom_nroots(state_sym)
          end if
       end if
    end do

  end subroutine

   subroutine input_processor_right(Fbar_HH,Fbar_PP,Vbar_SS,R1,R2,H_II)
      use davidson 

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"

!description
!===========
!1)  this subroutine processes the guess vectors for eomcc.
!2)  it also defines the Preconditioner for the error vectors, which is relevant at
!    the Davidson step.

!calling variables
!================

   real*8, intent(in)  :: Fbar_HH(:),Fbar_PP(:),Vbar_SS(:)
   real*8, intent(inout) :: R1(:,:), R2(:,:) ! (ndimt1,nroot) and (ndimt2,nroot)
   real*8, intent(inout),target :: H_II(:)

!Local Variables
!==============

   integer :: i,ii,irep,j,iaai,aoff,a,aa,ai,arep,iroot,k(1),jroot,newarep
   real*8  :: fac1
   real*8  :: eps((iv(nrep+1)+io(nrep+1))) 
   real*8, pointer :: diag_H(:)
   integer :: icount, acount

   integer :: abij, ijrp, abrp, jrp, jj, irp, ioff, imin, bb, amin, arp, brp, b
   real(kind=8) :: fac, fac3, fac2
   integer :: ndimr1, ndimr2, search_dimension, xroot, N1, N2, N
   real*8,allocatable :: R1_dum(:,:) 
   integer :: fudge_factor


!initialize the R-vector
!=======================

    N1 = size(R1,1)
    N2 = size(R2,1)
    ndimr1 = N1 / rcw
    ndimr2 = N2 / rcw

    R1(1:N1,1:nroot) = 0.0d0 
    R2(1:N2,1:nroot) = 0.0d0 

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

   case ("EE")

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
      DO JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) CYCLE
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            IF (ARP.LT.BRP) CYCLE
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
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
!aspg+as

      if (relcc_mfd_trial_full_matrix) then
! New trial vectors for R1 
          do i = 1, min(ndimr1,nroot)
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1(iroot,jroot) = 1.0d0
          end do
          if (nroot.gt.ndimr1) then
              do i = 1, min(ndimr2,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
                  R2(iroot,ndimr1+jroot) = 1.0d0
             end do
          end if

      else if (relcc_mfd_trial_ccs) then
          write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
          allocate (R1_dum(N1,ndimr1))
          R1_dum = 0.0d0

          do i = 1, ndimr1
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1_dum(iroot,jroot) = 1.0d0
          end do
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector         ")

          call trial_vectors_ee(R1_dum,Fbar_HH,Fbar_PP,Vbar_SS)
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector  CD1    ",complex_dimension=1)

          R1 = R1_dum(:,1:nroot)
          if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")

          deallocate(R1_dum)


      else if (relcc_mfd_trial_diagonal) then
          write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a
! particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

          allocate(diag_H(N1))
          diag_H = H_II(1:N1)

          xroot = 1

          fudge_factor = 2
! here we set the factor to 1 because for C1, with only one fermion irrep,
! we can't get around to working with the degenerate roots
          if (nrep.eq.1) fudge_factor = 1

          do iroot = 1, fudge_factor*nroot
              k = minloc(diag_H)
              if (iprnt.ge.3) write(iw,*)'roots',iroot,H_II(k(1)),k(1)
              i = k(1)
              diag_H(i) = 8000.0d0
              jroot = i
              if (rcw.eq.2) jroot = 2*i - 1
              R1(jroot,xroot) = 1.0d0
              xroot = xroot + 1 
          enddo
          deallocate(diag_H)

      else if (relcc_mfd_trial_restart) then
          stop('restart from prior RHS vectors not implemented')
      else
          stop('unknown RHS trial vectors generation scheme')
      end if ! if chain for trial vectors 


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
      DO JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) CYCLE
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO BRP = 1, NREP
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
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      if (relcc_mfd_trial_full_matrix) then
! New trial vectors for R1 
          do i = 1, min(ndimr1,nroot)
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1 
              jroot = i 
              R1(iroot,jroot) = 1.0d0
          end do
          if (nroot.gt.ndimr1) then
              do i = 1, min(ndimr2,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1 
                  jroot = i 
                  R2(iroot,ndimr1+jroot) = 1.0d0
             end do
          end if

      else if (relcc_mfd_trial_ccs) then
          write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
          allocate (R1_dum(N1,ndimr1))
          R1_dum = 0.0d0

          do i = 1, ndimr1
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1_dum(iroot,jroot) = 1.0d0
          end do
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector         ")

          call trial_vectors_ip(R1_dum,Fbar_HH)
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector  CD1    ",complex_dimension=1)
      
          R1 = R1_dum(:,1:nroot) 

          if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")

          deallocate(R1_dum)

      else if (relcc_mfd_trial_diagonal) then
          write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

          search_dimension = ncont(state_sym)*no(state_sym) 

          allocate(diag_H(N1))
          diag_H = H_II(1:N1)

          xroot = 1 

          fudge_factor = 2
! here we set the factor to 1 because for C1, with only one fermion irrep,
! we can't get around to working with the degenerate roots
          if (nrep.eq.1) fudge_factor = 1

          do iroot = 1, fudge_factor*nroot
              k = minloc(diag_H) 
              if (iprnt.ge.3) write(iw,*)'roots',iroot,H_II(k(1)),k(1) 
              i = k(1)
              diag_H(i) = 8000.0d0
              if (i <= search_dimension) then
                  jroot = i
                  if (rcw.eq.2) jroot = 2*i - 1
                  R1(jroot,xroot) = 1.0d0
                  if (rcw.eq.2) then
                      xroot = xroot + 1
                  else
                      xroot = xroot + 1
                  endif
              end if
          enddo
          deallocate(diag_H)

      else if (relcc_mfd_trial_restart) then
          stop('restart from prior RHS vectors not implemented')
      else
          stop('unknown RHS trial vectors generation scheme')
      end if ! if chain for trial vectors 


   case ("EA")

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
      DO JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         DO I = 1, NCONT(IRP)
            FAC2 =  FAC1
            DO BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            IF (ARP.LT.BRP) CYCLE
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
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      if (relcc_mfd_trial_full_matrix) then
! New trial vectors for R1 
          do i = 1, min(ndimr1,nroot)
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1(iroot,jroot) = 1.0d0
          end do
          if (nroot.gt.ndimr1) then
              do i = 1, min(ndimr2,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
                  R2(iroot,ndimr1+jroot) = 1.0d0
             end do
          end if

      else if (relcc_mfd_trial_ccs) then
          write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
          allocate (R1_dum(N1,ndimr1))
          R1_dum = 0.0d0

          do i = 1, ndimr1
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1_dum(iroot,jroot) = 1.0d0
          end do
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector         ")

          call trial_vectors_EA(R1_dum,Fbar_PP)
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector  CD1    ",complex_dimension=1)

          R1 = R1_dum(:,1:nroot)

          if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")

          deallocate(R1_dum)

! calculate the orbital dimension within which our desired root should belong..

      else if (relcc_mfd_trial_diagonal) then
          write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

          search_dimension = ncont(state_sym)*nv(state_sym)

          allocate(diag_H(N1))
          diag_H = H_II(1:N1)

          xroot = 1

          fudge_factor = 2
! here we set the factor to 1 because for C1, with only one fermion irrep,
! we can't get around to working with the degenerate roots
          if (nrep.eq.1) fudge_factor = 1

          do iroot = 1, fudge_factor*nroot
              k = minloc(diag_H)
              if (iprnt.ge.3) write(iw,*)'roots',iroot,H_II(k(1)),k(1)
              i = k(1)
              diag_H(i) = 8000.0d0
              if (i <= search_dimension) then
                  jroot = i
                  if (rcw.eq.2) jroot = 2*i - 1
                  if (jroot.le.N1) then
                      R1(jroot,xroot) = 1.0d0
                  else if (jroot.le.N2) then
                      R2(jroot,xroot) = 1.0d0
                  end if
                  if (rcw.eq.2) then
                      xroot = xroot + 1
                  else
                      xroot = xroot + 1
                  endif
              end if
          enddo
          deallocate(diag_H)

      else if (relcc_mfd_trial_restart) then
          stop('restart from prior RHS vectors not implemented')
      else
          stop('unknown RHS trial vectors generation scheme')
      end if ! if chain for trial vectors 

  end select 

   end subroutine


   subroutine input_processor_left(Fbar_HH,Fbar_PP,Vbar_SS,eVectorsR,L1,L2,H_II)

  use symmetry_offset
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"


!description
!===========
!1)  this subroutine processes the left hand guess vectors for eomcc.
!2)  it also defines the Preconditioner for the error vectors, which is relevant at
!    the Davidson step.

!calling variables
!================

   real*8, intent(in)  :: Fbar_HH(:),Fbar_PP(:),Vbar_SS(:),eVectorsR(:,:) 
   real*8, intent(inout) :: L1(:,:), L2(:,:) ! (ndimt1,nroot) and (ndimt2,nroot)
   real*8, intent(inout),target :: H_II(:)

!Local Variables
!==============

   integer :: i,ii,irep,j,iaai,aoff,a,aa,ia,arep,iroot,k(1)
   real*8  :: fac1
   real*8  :: eps((iv(nrep+1)+io(nrep+1))) 
   real*8, pointer :: diag_H(:)
   integer :: icount, acount

   integer :: ijab, ijrp, abrp, jrp, jj, irp, ioff, imin, bb, amin, arp, brp, b
   real(kind=8) :: fac, fac3, fac2
   integer :: ndimL1, ndimL2, search_dimension, xroot, N, N1, N2, fudge_factor, jroot
   type(Offset)         :: e
   real*8,allocatable :: L1_dum(:,:) 


!initialize the R-vector
!=======================

    N1 = size(L1,1)
    N2 = size(L2,1)
    N  = N1 + N2

    ndimL1 = N1 / rcw
    ndimL2 = N2 / rcw

!   L1(1:N1,1:nroot) = 0.0d0 
!   L2(1:N2,1:nroot) = 0.0d0 
    L1 = 0.0d0 
    L2 = 0.0d0 

! 1> store the diagonal elements of the transformed Hamiltonian..  

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

   case ("EE")

    ia = 1
    ii = 0

      do arp = 1, nrep
         irep = multb (arep,state_sym+nrep,1)
         aa = io(nrep+1) + iv(arp)
      do a = 1, nv(arp)
         aa = aa + 1
         fac1 = eps(aa)
         ioff = io(arp)
         do i = 1, no(arp)
             ii = ioff + i
             h_ii(ia) = fac1 - eps(ii) ! I will add the missing two-electron term later 
             ia = ia + 1
         enddo
      enddo
      enddo


      ijab = ia
      do abrp = 1, nrep
      do brp = 1, nrep
      bb = iv(brp) + io(nrep+1)
      arp = multb(brp,abrp+nrep,2)
      if (arp.lt.brp) cycle
      do b = 1, nv(brp)
         bb = bb + 1
         fac1 = eps(bb)
            aoff = iv(arp) + io(nrep+1)
           amin = 1
           if (arp.eq.brp) amin = b + 1
         do a = amin, nv(arp)
            aa = aoff + a
            fac2 = eps(aa) + fac1
            do jrp = 1, nrep
            jj = io(jrp)
            irp = multb(jrp,abrp+nrep,2)
            if (irp.lt.jrp) cycle
            ioff = io(irp)
            do j = 1, no(jrp)
               jj = jj + 1
               fac3 = fac2 - eps(jj)
               imin = 1
               if (irp.eq.jrp) imin = j + 1
               do i = imin, no(irp)
                  ii = ioff + i
                  fac = fac3 - eps(ii)
                  h_ii(ijab) = fac
                  ijab = ijab + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo


    do iroot = 1, nroot

      call srt1c1 (nrep,nv,no,eVectorsR(1:N1,iroot:iroot), &
     &             L1(1:N1,iroot:iroot)) 
      call srt1c1 (nrep,nvvt,noot,eVectorsR(N1+1:N,iroot:iroot), &
     &             L2(1:N2,iroot:iroot)) 

    enddo

   case ('IP')

    ia = 1
    ii = 0

      do arp = 1, nrep
         irep = multb (arp,1+nrep,1)
         aa = io(nrep+1) + iv(arp)
      do a = 1, ncont(arp)
         ioff = io(irep)
         do i = 1, no(irep)
             ii = ioff + i
             h_ii(ia) = -eps(ii)  
             ia = ia + 1
         enddo
      enddo
      enddo


      ijab = ia
      do abrp = 1, nrep
      do brp = 1, nrep
      bb = iv(brp) + io(nrep+1)
      arp = multb(brp,abrp+nrep,2)
      do b = 1, nv(brp)
         bb = bb + 1
         fac1 = eps(bb)
         do a = 1, ncont(arp)
            fac2 =  fac1
            do jrp = 1, nrep
            jj = io(jrp)
            irp = multb(jrp,abrp+nrep,2)
            if (irp.lt.jrp) cycle
            ioff = io(irp)
            do j = 1, no(jrp)
               jj = jj + 1
               fac3 = fac2 - eps(jj)
               imin = 1
               if (irp.eq.jrp) imin = j + 1
               do i = imin, no(irp)
                  ii = ioff + i
                  fac = fac3 - eps(ii)
                  h_ii(ijab) = fac
                  ijab = ijab + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo

! In the following section we generate the guess vectors belonging to a particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,nv,ncont,.false.,.false.)

         
 !   search_dimension = ncont(state_sym)*no(state_sym) 

 !   allocate(diag_H(ndiml1))
 !   diag_H = H_II(1:ndiml1)

 !      xroot = 1 

 !     do iroot = 1, 2*nroot
 !      k = minloc(diag_H) 
!!ashee        write(iw,*)'roots',iroot,H_II(k(1)),k(1) 
 !      i = k(1)
 !      diag_H(i) = 8000.0d0

 !      if (i <= search_dimension) then
 !       L1(i,xroot) = 1.0d0
 !       xroot = xroot + 1
 !      endif

 !     enddo
 !   deallocate(diag_H)

      if (relcc_mfd_trial_lhs_use_rhs) then
         write (iw,*) '  Using RHS eigenvectors as LHS trial vectors'
! if we are reusing the right-hand side eigenvectors, we resort them

          do iroot = 1, nroot
              call srt1c1 (nrep,ncont,no,eVectorsR(1:N1,iroot:iroot), &
     &                     L1(1:N1,iroot:iroot))
              call srt1c1 (nrep,e%oneNonDirac,noot,eVectorsR(N1+1:N,iroot:iroot), &
     &                     L2(1:N2,iroot:iroot))
          enddo

      else
! otherwise, we are going to generate trial vectors for the lhs
          if (relcc_mfd_trial_full_matrix) then
! New trial vectors for L1 
              do i = 1, min(ndiml1,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
!                 print *, 'iroot',iroot,' jroot',jroot
                  L1(iroot,jroot) = 1.0d0
              end do
              if (nroot.gt.ndiml1) then
                  do i = 1, min(ndiml2,nroot)
                      iroot = i
                      if (rcw.eq.2) iroot = 2*i - 1
                      jroot = i
                      L2(iroot,ndiml1+jroot) = 1.0d0
                 end do
              end if

          else if (relcc_mfd_trial_ccs) then
              write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
              allocate (L1_dum(N1,ndiml1))
              L1_dum = 0.0d0

              do i = 1, ndiml1
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
                  L1_dum(iroot,jroot) = 1.0d0
              end do
              if (iprnt.ge.3) call davidson_print_vectors(L1_dum, label="debug: L1_dum vector         ")

! this is a temporary hack (as of 31/03/2017), since the routine below only returns the right-hand vectors
! to generalize
              call trial_vectors_ip(L1_dum,Fbar_HH) 
              if (iprnt.ge.3) call davidson_print_vectors(L1_dum, label="debug: L1_dum vector  CD1    ",complex_dimension=1)

              L1 = L1_dum(:,1:nroot)

              if (iprnt.ge.3) call davidson_print_vectors(L1, label="debug:  R1 vector             ")

              deallocate(L1_dum)

          else if (relcc_mfd_trial_diagonal) then
              write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a
! particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

              search_dimension = ncont(state_sym)*no(state_sym)

              allocate(diag_H(N1))
              diag_H = H_II(1:N1)

              xroot = 1

              fudge_factor = 2
! here we set the factor to 1 because for C1, with only one fermion irrep,
! we can't get around to working with the degenerate roots
              if (nrep.eq.1) fudge_factor = 1

              do iroot = 1, fudge_factor*nroot
                  k = minloc(diag_H)
                  if (iprnt.ge.3) write(iw,*)'roots',iroot,H_II(k(1)),k(1)
                  i = k(1)
                  diag_H(i) = 8000.0d0
                  if (i <= search_dimension) then
                      jroot = i
                      if (rcw.eq.2) jroot = 2*i - 1
                      L1(jroot,xroot) = 1.0d0
                      if (rcw.eq.2) then
                          xroot = xroot + 1
                      else
                          xroot = xroot + 1
                      endif
                  end if
              enddo
              deallocate(diag_H)
          else if (relcc_mfd_trial_restart) then
              stop('restart from prior LHS vectors not implemented')
          else
              stop('unknown LHS trial vectors generation scheme')
          end if ! if chain for trial vectors 

      end if ! relcc_mfd_trial_lhs_use_rhs


  call dealloc_array(e)

   case ("EA")

    ia = 1
    aa = 0

      do arp = 1, nrep
         irep = multb (arep,state_sym+nrep,1)
         aa = io(nrep+1) + iv(arp)
      do a = 1, nv(arp)
         aa = aa + 1
         fac1 = eps(aa)
         do i = 1, ncont(arp)
             h_ii(ia) = fac1  
             ia = ia + 1
         enddo
      enddo
      enddo

      ijab = ia
      do abrp = 1, nrep
      do brp = 1, nrep
      bb = iv(brp) + io(nrep+1)
      arp = multb(brp,abrp+nrep,2)
      if (arp.lt.brp) cycle
      do b = 1, nv(brp)
         bb = bb + 1
         fac1 = eps(bb)
            aoff = iv(arp) + io(nrep+1)
           amin = 1
           if (arp.eq.brp) amin = b + 1
         do a = amin, nv(arp)
            aa = aoff + a
            fac2 = eps(aa) + fac1
            do jrp = 1, nrep
            jj = io(jrp)
            irp = multb(jrp,abrp+nrep,2)
            do j = 1, no(jrp)
               jj = jj + 1
               fac3 = fac2 - eps(jj)
               do i = 1, ncont(irp)
                  fac = fac3  
                  h_ii(ijab) = -fac
                  ijab = ijab + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo

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

     allocate(diag_H(N1))
     diag_H = H_II(1:N1)
        xroot = 1
       do iroot = 1, 2*nroot
        k = minloc(diag_H) 
        i = k(1)
        diag_H(i) = 8000.0d0
        if (i <= search_dimension) then
         L1(i,xroot) = 1.0d0
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
   use eom_density
   use symmetry_offset 

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dgroup.h"
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
   integer :: ndimr1,ndimr2,sym
   integer :: AllocateStatus

   type(davidson_results) :: davResultsR(nrep), davResultsL(nrep)
   logical :: left_vectors, right_vectors
   integer :: N, N1, N2 ! size of determinantal space 
   type(Offset)         :: e
   real*8               :: overlap 

!memory allocation
!------------------

  integer(kind=8)        :: current_mem_use = 0
  integer(kind=8)        :: start_mem_use   = 0
  integer                :: freemem = 0
  integer                :: work_memsize
  integer                :: nbuf1,nbuf2
  integer                :: nbuf3
  integer                :: state
  integer*8              :: nv5irp,nv6irp,nvirp
  real*8,allocatable     :: dmo(:) 
  real*8,allocatable     :: l1_temp(:),l2_temp(:) 

  real*8 :: difference_in_evalue
  real*8 :: rhs_lhs_evalue(2)
  real*8 :: threshold
  double precision dznrm2

    nbuf1 = max(joooo(nrep+1),jvooo(nrep+1),jvvoo(nrep+1),jvovo(nrep+1))
    nbuf2 = nbuf1

    call alloc( buf1  , nbuf1*rcw             , id="buf1" )
    call alloc( buf2  , nbuf2*rcw             , id="buf2" )

    call allocator_get_words_inuse(start_mem_use)
    call legacy_lwork_get(work_memsize)

    call allocator_get_maxbuff(freemem,kind(BUF3))

    NBUF3 = freemem / RCW

    NV5IRP = MAX(JVOVV(NREP+1),JOVVV(NREP+1))
    NV6IRP = 1
    DO IRP = 1, NREP
       NV6IRP = MAX(NV6IRP,INT(NVVT(IRP),8)*INT(NVVT(IRP),8))
    ENDDO

!   Check the minimum size for the general use of this array

    NVIRP = MAX(NV5IRP,NV6IRP)

    NBUF3 = MIN(NVIRP,INT(NBUF3,8))
    NBUF3 = MAX(NBUF1,NBUF3)

!     Allocate now also this remaining work array

    call alloc( BUF3  , NBUF3*RCW             , id="buf3" )

!   For the moment I am using as much in-core memory as possible. However we
!   may think of writing them in files, replacing the integral ones...  
! aspg: we should collect these actions into a routine, like
!    call initialize_intermediates() 
    call alloc(B%Fbar_mi,   (nfoo*rcw), id="Fbar_mi")  
    call alloc(B%Fbar_me,   (nfvo*rcw), id="Fbar_me")  
    call alloc(B%Fbar_ae,   (nfvv*rcw), id="Fbar_ae")  

    call alloc(B%W_ijmn,    (nv1*rcw),            id="W_ijmn") 
    call alloc(B%W_ejmb,    (iovvo(nrep+1)*rcw),  id="W_ejmb")
    call alloc(B%W_mbej,    (iovvo(nrep+1)*rcw),  id="W_mbej")
    call alloc(B%W_iemn,    (iovoot(nrep+1)*rcw), id="W_iemn")

! todo: get special handling for these in the case of i4 ? we keep using deallocate here for the time being
!   call alloc(B%W_efam,    (nv5*rcw),            id="W_efam")
!   call alloc(B%W_amef,    (nv5*rcw),            id="W_amef")
    allocate (B%W_efam(nv5*rcw))
    allocate (B%W_amef(nv5*rcw))

    call alloc(B%W_mnie,    (iovoot(nrep+1)*rcw), id="W_mnie")
    call alloc(B%Wbar_mbej, (iovvo(nrep+1)*rcw),  id="Wbar_mbej")

    call alloc(B%t1, (ndimt1*rcw), id="t1")
    call alloc(B%t2, (ndimt2*rcw), id="t2")
       
! allocate the amplitudes

    call allocator_get_words_inuse(current_mem_use)
    mxcorr = current_mem_use - start_mem_use

    if (icalc.eq.1) then
       write (iw,1001) "EOMCC calculation",mxcorr
       goto 451
    endif

    write (iw,*) " "
    write (iw,*) "Equation of Motion CC module setup"
    write (iw,*) " "

    call verify_input()
    call initialize_buffers(buf1, buf2, buf3)
    call ccdump
    call getampt(B%t1,B%t2)
    call IM_container(eomtype,B)
    call free_buffers(buf1,buf2,buf3)

    right_vectors = .true.
    left_vectors =  .false.
    call eom_solve_for_left_or_right(left_vectors, right_vectors, davResultsL, davResultsR, B, eps)

    if (relcc_do_eomprop) then

       right_vectors = .false.
       left_vectors =  .true.
       call eom_solve_for_left_or_right(left_vectors, right_vectors, davResultsL, davResultsR, B, eps)

!!!!!
! turn the loops below into a separare routine ?
!      call eom_check_left_right_roots(davResultsL,davResultsR)
       do sym = 1, nrep  
           if (relcc_eom_nroots(sym) <= 0) cycle
           write (iw,*) " "
           do iroot = 1, relcc_eom_nroots(sym)
               rhs_lhs_evalue(1) = davResultsR(sym)%eValues(1,iroot) - davResultsL(sym)%eValues(1,iroot)
               rhs_lhs_evalue(2) = davResultsR(sym)%eValues(2,iroot) - davResultsL(sym)%eValues(2,iroot)
               difference_in_evalue = dznrm2(1,rhs_lhs_evalue,1)
               write (iw,'(2X,A,I4,A,I3,2X,A,2X,E16.8)') &
                      "INFO: difference between RHS and LHS eigenvalues for root ",iroot,&
     &                " irrep ",sym,":",difference_in_evalue
               write (iw,'(2X,A)') "INFO:              Re                  Im"
               write (iw,'(2X,A,2X,E16.8,4X,E16.8)') "INFO: E(RHS)   ",davResultsR(sym)%eValues(:,iroot)
               write (iw,'(2X,A,2X,E16.8,4X,E16.8)') "INFO: E(LHS)   ",davResultsL(sym)%eValues(:,iroot)
           write (iw,*) " "
           end do  
       end do  
!!!!!
       write (iw,*) "Determining EOM density matrices"

       allocate(dmo(norbt*norbt*nz))

       dmo =0.0d0
       state =1
       do sym = 1, nrep  

           if (relcc_eom_nroots(sym) <= 0) cycle

           do iroot = 1, relcc_eom_nroots(sym)
               select case(eomtype)
                  case ('IP')
!                    call set_continuous_orbital(sym)
                     ndimr1 = 0     

                     do irep = 1, nrep
                         ndimr1 = ndimr1 + ncont(irep)*no(irep)
                     enddo

                     call alloc_array(e,nrep)
                     call auto_symmetry_offset(e,ncont,nv,.false.,.false.)
     
                     ndimr2 = 0     

                     do irep = 1, nrep
                         ndimr2 = ndimr2 + e%oneNonDirac(irep)*noot(irep)
                     enddo

                     N1 = ndimr1*rcw
                     N2 = ndimr2*rcw
                     N  = N1 + N2

                     allocate(l1_temp(N1))
                     allocate(l2_temp(N2))

                     call srt1c1(nrep,no,ncont,davResultsL(sym)%eVectorsL(1:N1,iroot:iroot),l1_temp)
                     call srt1c1(nrep,noot,e%oneNonDirac,davResultsL(sym)%eVectorsL(1+N1:N,iroot:iroot),l2_temp)

                     call dealloc_array(e)

                     overlap = dot_product(l1_temp,davResultsR(sym)%eVectorsR(1:N1,iroot)) + &
    &                          dot_product(l2_temp,davResultsR(sym)%eVectorsR(1+N1:N,iroot))

                     write(*,*)'overlap', overlap 

                     overlap = 1.0d0/overlap

                     call dscal(N,overlap,davResultsL(sym)%eVectorsL(1:N,iroot:iroot),1)  
                     call dscal(N1,overlap,l1_temp,1)  
                     call dscal(N2,overlap,l2_temp,1)  

                     write(*,*) dot_product(l1_temp,davResultsR(sym)%eVectorsR(1:N1,iroot)) + &
    &                           dot_product(l2_temp,davResultsR(sym)%eVectorsR(1+N1:N,iroot))

                     deallocate(l1_temp,l2_temp)

                     call density_IP(davResultsL(sym)%eVectorsL(1:N1,iroot),&
  &                                  davResultsL(sym)%eVectorsL(1+N1:N,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1:N1,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1+N1:N,iroot),B%t1,B%t2,dmo)

                  case ('EE')

                     ndimr1 = ndimt1
                     ndimr2 = ndimt2

                     N1 = ndimr1*rcw
                     N2 = ndimr2*rcw
                     N  = N1 + N2

                     allocate(l1_temp(N1))
                     allocate(l2_temp(N2))

                     call srt1c1(nrep,no,nv,davResultsL(sym)%eVectorsL(1:N1,iroot),l1_temp)
                     call srt1c1(nrep,noot,nvvt,davResultsL(sym)%eVectorsL(1+N1:N,iroot),l2_temp)

                     overlap = dot_product(l1_temp,davResultsR(sym)%eVectorsR(1:N1,iroot)) + &
    &                          dot_product(l2_temp,davResultsR(sym)%eVectorsR(1+N1:N,iroot))

                     write(*,*)'overlap', overlap 

                     overlap = 1.0d0/overlap

                     call dscal(N,overlap,davResultsL(sym)%eVectorsL(1:N,iroot:iroot),1)  
                     call dscal(N1,overlap,l1_temp,1)  
                     call dscal(N2,overlap,l2_temp,1)  

                     write(*,*) dot_product(l1_temp,davResultsR(sym)%eVectorsR(1:N1,iroot)) + &
    &                           dot_product(l2_temp,davResultsR(sym)%eVectorsR(1+N1:N,iroot))
 
!                    write(*,*)dot_product(davResultsR(sym)%eVectorsR(1:,iroot),davResultsR(sym)%eVectorsR(1:,iroot))

                     deallocate(l1_temp,l2_temp)

                     call density_EE(davResultsL(sym)%eVectorsL(1:N1,iroot),&
  &                                  davResultsL(sym)%eVectorsL(1+N1:N,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1:N1,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1+N1:N,iroot),B%t1,B%t2,dmo)

                 case ('EA')
                   write(*,*)'Not yet implemented'

               end select

               write (iw,*) "debug, before store_eom_density for state",state
               call store_eom_density ('EOM ',state,dmo)
               write (iw,*) "debug,  after store_eom_density for state",state
               state = state + 1
           enddo ! loop over roots 
       enddo ! loop over symmetries
       deallocate(dmo)
   endif ! relcc_do_eomprop

     do sym = 1, nrep 

         if (relcc_eom_nroots(sym) <= 0) cycle

          if (right_vectors) call davidson_cleanup(davResultsR(sym))
          if (left_vectors)  call davidson_cleanup(davResultsL(sym))
     end do

 451 continue

!
! final cleanup
!
     call free_intermediates(B)
     if (icalc.eq.1) then
        call free_buffers(buf1,buf2,buf3)
     end if

1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
   end subroutine


   subroutine number_of_states(state)
#include "symm.inc"

   integer, intent(out) :: state 
   integer              :: sym

   state = 0
   do sym = 1, nrep  

     state = state + relcc_eom_nroots(sym) 

   enddo

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
      INTEGER , ALLOCATABLE :: INDX_ALL(:), IRPS_ALL(:)  
      INTEGER :: IERR, ISIZE, select_eom,dum2
      REAL*8 :: IND(MXREP),dum1(rcw),eigenvalues_temp(1:rcw,0:nroot)
      REAL*8 :: E_OFF
      INTEGER :: nstate_sym(nrep,4)
!---------------Executable code--------------------------------------

    dum1(1:rcw) = 0.0d0
    dum2 = 0

! Generate an array denoting number of states per symmetry...

   nstate_sym(1:nrep,1:4) = 0

   select case(eomtype)

    case('IP')

      nstate_sym(state_sym,3) = nroot
      select_eom = 3
      ISIZE = nroot

    case('EA')

      nstate_sym(state_sym,2) = nroot
      select_eom = 2
      ISIZE = nroot 

    case('EE')

     nstate_sym(state_sym,4) = nroot
     select_eom = 4
     ISIZE = nroot + 1 

  end select  

!     ISIZE = nroot + 1 ! possibly correct
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
      ELSE IF ((select_eom==4)) THEN 
!        Add energy of (0,0) sector to list of printed energies for the (1,1) sector
!        This makes the routine print excitation energies
         nstate_sym(1,select_eom) = nstate_sym(1,select_eom) + 1

         EIGENVALUES_temp(1:RCW,1:nroot)=eigenvalues(1:RCW,1:nroot)
         EIGENVALUES_temp(1:RCW,0) = 0.D0

!        Even electron case : boson irreps
         CALL PRT_EV (IW,NREP,nstate_sym(1,select_eom),REPNA(NREP+1), &
     &                E_OFF,ECCSDIM,RCW, &
     &                EIGENVALUES_temp(1:RCW,0:ISIZE),ISIZE, &
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

!     deallocate(nstate_sym)

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

      if (ipr.le.1) then
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

    jj = 0
  do irep = 1, nrep
     arep = multb (irep,1+nrep,1)
      do i = 1, no(irep)
       jj = jj + 1
          aoff = io(nrep+1) + iv(arep)
        do a = 1, nv(arep)
          aa = aoff + a
          write(label,1004)repna(irep),i,eps(jj), &
            repna(arep),(a+no(arep)),eps(aa)
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
            if (arp.lt.brp) cycle
            aoff = iv(arp) + io(nrep+1)
            do b = 1, nv(brp)
               bb = bb + 1
               amin = 1
               if (arp.eq.brp) amin = b + 1
               do a = amin, nv(arp)
                  aa = aoff + a
                  write(label,1005) repna(irp),i,eps(kk), &
                &          repna(jrp),j,eps(jj),&
                &          repna(arp),(a+no(arp)),eps(aa),&
                &          repna(brp),(b+no(brp)),eps(bb)

                  value = abs(eigenvectors(1,ii))
                  if (carith) value = value + abs(eigenvectors(2,ii))
                  if (carith .and. value.gt.tres) then
                      write (iw,'(2f15.5,3x,a70)') & 
            &             eigenvectors(1,ii),eigenvectors(2,ii),label
                  else if (value.gt.tres) then              
                  write (iw,'(f15.5,3x,a100)') eigenvectors(1,ii),label
                  end if
                  ii = ii + 1
               enddo
            enddo
            enddo
         enddo
      enddo
    enddo
  enddo

      end select
      endif
         end do
      end do
 
 1000 FORMAT (//' Analysis of EOM eigenvectors' &
     &/' First line  :  State Energy, Eigenvalue' &
     &/' other lines : Largest R1 and R2 Vectors',&
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



subroutine eom_solve_for_left_or_right(left_vectors, right_vectors, davResultsL, davResultsR, B, eps)
   use davidson
   use allocator_parameters, only : klongint
   use symmetry_offset

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"
#include "ccpar.inc"
!--------------------------------------------------
   type(davidson_results), intent(inout) :: davResultsL(nrep)
   type(davidson_results), intent(inout) :: davResultsR(nrep)
   type(intermediates), intent (inout)   :: B
   real(kind=8)                          :: eps(*)
   logical, intent(in)                   :: left_vectors
   logical, intent(in)                   :: right_vectors
!--------------------------------------------------
   real*8, pointer :: trial_vector(:,:)
   integer :: irp,iroot,jroot,i,j,irep,jrep,ijrep,nvx(nrep),nox(nrep)
   integer :: ndimr1,ndimr2
   integer :: AllocateStatus
   integer :: N, N1, N2 ! size of determinantal space 
   character(len=3) :: left_or_right
   type(Offset) :: e,f,g,h

! 
!  we set the state for the diagonalizer, as this should not change within a given type of calculation,
! e.g. excitation energies
!
    call davidson_setup(unit_output= iw,   &
                        solve_right=right_vectors, &
                        solve_left=left_vectors,   &
                        symmetric=.false., &
                        verbose=relcc_mfd_verbose,   &
                        complex_mode=(rcw.eq.2), &
                        convergence_threshold = relcc_mfd_convergence_threshold, &
                        max_subspace_size = relcc_mfd_max_subspace_size, &
                        max_iterations = relcc_mfd_max_iterations, &
                        overlap_sorting = relcc_mfd_overlap_sorting)

    if (left_vectors .and. right_vectors) call quit('EOM does left OR right vectors for now')

    if (right_vectors) then
        left_or_right = "RHS"
    else if (left_vectors) then
        left_or_right = "LHS"
    end if
!
! be very careful with state_sym and nroot, now these are global variables in this module !!
!
    do state_sym = 1, nrep ! 16 
       call create_eom_ndet(2)
       nroot = relcc_eom_nroots(state_sym)

       if (nroot.ne.0) then
          write (iw,*) ""
          write (iw,'(2x,A42,2x,I3)') "Solving "//left_or_right(1:3)//" EOM-"//eomtype(1:2)//" equations for symmetry",state_sym
!         write (iw,'(2x,A,2x,I3)') "Solving EOM-"//eomtype(1:2)//" equations for symmetry",state_sym
          write (iw,*) ""

! setting up information for continuous orbital. this is a trick to perform
! EOM-IP and EOM-EA calculations within the same framework as EOM-EE 
       call set_continuous_orbital(state_sym)

       select case (eomtype)

          case ('EE')

      call alloc_array(e,nrep,f,g,h)
      call auto_symmetry_offset(e,nv,no,.false.,.false.)
      call auto_symmetry_offset_triangular(f,nv,nv)
      call auto_symmetry_offset_triangular(g,no,no)
      call auto_symmetry_offset(h,f%oneNonDirac,g%oneNonDirac,.true.,.true.)

         ndimr1 = e%oneNonDirac(state_sym)
         ndimr2 = h%oneNonDirac(state_sym)
 
      call dealloc_array(e,f,g,h)

             if (nroot.eq.-1) nroot = ndimr1 + ndimr2

             N1 = ndimr1*rcw
             N2 = ndimr2*rcw
             N =  N1 + N2 

             call set_eom_ndet(1, N1)
             call set_eom_ndet(2, N2)

             allocate (trial_vector(N,nroot), stat = AllocateStatus)
             IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
             allocate (B%H_II(N),  stat = AllocateStatus)
             IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

             trial_vector = 0.0d0
             call eom_initialize_HII(B%H_II)
!            B%H_II = 1.0d0

             if (right_vectors.and..not.left_vectors) then
                 call input_processor_right(B%Fbar_mi, &
      &                               B%Fbar_ae, &
      &                               B%W_mbej,  &
      &                               trial_vector(1:N1,:), &
      &                               trial_vector(N1+1:N,:),  &
      &                               B%H_II)
             else if (left_vectors.and..not.right_vectors) then    
                 call input_processor_left(B%Fbar_mi, &
      &                               B%Fbar_ae, &
      &                               B%W_mbej,davResultsR(state_sym)%eVectorsR(1:N,:),  &
      &                               trial_vector(1:N1,:), &
      &                               trial_vector(N1+1:N,:),  &
      &                               B%H_II)
 
             end if


          case ('IP')

             ndimr2 = 0
             ndimr1 = 0
             nvx = 0
       
             call create_eom_nvx(nrep)

! this need a bit of explaining
             do jrep = 1, nrep
                do irep = 1, nrep
                   ijrep = multb(irep,jrep,1)
                   nvx(ijrep) = nvx(ijrep) + ncont(irep)*nv(jrep)
                   call set_eom_nvx(ijrep, (get_eom_nvx(ijrep) + ncont(irep)*nv(jrep)))
                enddo
             enddo

             do irep = 1, nrep
                if (iprnt.ge.3) write(*,*)'continuous orbital ',irep,ncont(irep),nroot
                ndimr1 = ndimr1 + ncont(irep)*no(irep)
             enddo

             do irep = 1, nrep
                 ndimr2 = ndimr2 + get_eom_nvx(irep)*noot(irep)
             enddo

             N1 = ndimr1*rcw
             N2 = ndimr2*rcw
             N = N1 + N2 

             call set_eom_ndet(1, N1)
             call set_eom_ndet(2, N2)

             if (nroot.eq.-1) nroot = ndimr1 + ndimr2

             allocate (trial_vector(N,nroot), stat = AllocateStatus)
             IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
             allocate (B%H_II(N)) 
       
             trial_vector = 0.0d0
             call eom_initialize_HII(B%H_II)
 !           B%H_II = 1.0d0

             if (right_vectors.and..not.left_vectors) then
                 call input_processor_right(B%Fbar_mi, &
      &                               B%Fbar_ae, &
      &                               B%W_mbej,  &
      &                               trial_vector(1:N1,:), &
      &                               trial_vector(N1+1:N,:),  &
      &                               B%H_II)
             else if (left_vectors.and..not.right_vectors) then    
                 call input_processor_left(B%Fbar_mi, &
      &                               B%Fbar_ae, &
      &                               B%W_mbej,davResultsR(state_sym)%eVectorsR(1:N,:),  &
      &                               trial_vector(1:N1,:), &
      &                               trial_vector(N1+1:N,:),  &
      &                               B%H_II)

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

             N1 = ndimr1*rcw
             N2 = ndimr2*rcw
             N = N1 + N2

             call set_eom_ndet(1, N1)
             call set_eom_ndet(2, N2)

             if (nroot.eq.-1) nroot = ndimr1 + ndimr2

             allocate (trial_vector(N,nroot))
             allocate (B%H_II(N)) 
       
             trial_vector = 0.0d0

             call eom_initialize_HII(B%H_II)
!            B%H_II = 1.0d0

             if (right_vectors.and..not.left_vectors) then 
                 call input_processor_right(B%Fbar_mi, &
      &                               B%Fbar_ae, &
      &                               B%W_mbej,  &
      &                               trial_vector(1:N1,:), &
      &                               trial_vector(N1+1:N,:),  &
      &                               B%H_II)
             else if (left_vectors.and..not.right_vectors) then 
                 call input_processor_left(B%Fbar_mi, &
      &                               B%Fbar_ae, &
      &                               B%W_mbej,davResultsR(state_sym)%eVectorsR,  &
      &                               trial_vector(1:N1,:), &
      &                               trial_vector(N1+1:N,:),  &
      &                               B%H_II)
 
             end if

         end select 

         if (right_vectors) then
             call davidson_driver(nroot, B, davResultsR(state_sym), trial_vector,state_sym)
             call eom_analysis(davResultsR(state_sym)%eValues(1:rcw,:),davResultsR(state_sym)%eVectorsR,eps) ! Fock-Space like output for EA, EE and IP

         else if (left_vectors) then

             call davidson_driver(nroot, B, davResultsL(state_sym), trial_vector,state_sym)

         end if

         call free_eom_ndet()
         call free_eom_nvx()

         if (allocated(B%H_II)) deallocate(B%H_II)
! aspg: trial_vector should have been nullified inside davidson_driver, check
! why this crashes
!        if (associated(trial_vector)) deallocate(trial_vector)

      end if ! this ends the if checking whether there are any roots required for this symmetry
   end do    ! this ends the loop over symmetries

   if (right_vectors) then
       write (iw,*) ""
       write (iw,'(2x,A)') "Final results for "//left_or_right//" EOM-"//eomtype(1:2)
       write (iw,*) ""
       call eom_analysis_allirreps(davResultsR, eps)
   else if (left_vectors) then
       write (iw,*) ""
       write (iw,'(2x,A)') "Final results for "//left_or_right//" EOM-"//eomtype(1:2)
       write (iw,*) ""
       call eom_analysis_allirreps(davResultsL, eps)
   end if
   
end subroutine

subroutine eom_initialize_HII(H_II)
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"

    real(kind=8), intent(inout) :: H_II(:)
    integer :: i

    H_II = 0.0d0
    do i = 1, size(H_II,1), rcw
        H_II(i) = 1.0d0
    end do

end subroutine



subroutine eom_analysis_allirreps(davResults, eps)
    use davidson
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
!--------------------------------------------------
   type(davidson_results), intent(inout) :: davResults(nrep)
   real(kind=8)                          :: eps(*)
!--------------------------------------------------
   write (iw,*) "Analysis of all EOM roots not implemented yet !"
end subroutine

   end module 
