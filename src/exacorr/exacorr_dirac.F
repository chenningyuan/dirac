!
!    Collection of routines needed for exacorr that use the DIRAC common
!    block infrastructure.
!
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_MOLTRA()
C***********************************************************************
C
C     Driver routine for exacorr MOLTRA implementation
C
C     Written by L.Visscher Jan 2017
C
C***********************************************************************

      use exacorr_global,only   : get_talsh_only
      use exacorr_ao_to_mo,only : exacorr_compute_mo_integrals
      use talsh_ao_to_mo,only   : talsh_compute_mo_integrals
#if defined (VAR_MPI)
      use interface_to_mpi
#endif

#include "implicit.h"! because common blocks may assume this
#include "../moltra/dcbtra.h"  ! to access tra4_indstr
#include "infpar.h"  ! we need to broadcast tra4_indstr
#include "dcborb.h"  ! to access nesh (number of positive energy "electronlike" orbitals)
#include "dgroup.h"  ! to access nfsym (should be 1)

      character (len=72)   :: orbital_string(4) ! String with range of MO's to be transformed for each index
      integer              :: kvec              ! Dummy variable needed in first call to numlst
      integer, allocatable :: mo_list(:),mo_occ(:),mo_vir(:)
      integer              :: nmo(4),noccu,nvir

      if (nfsym.ne.1) then
         call quit ('exacorr can not handle symmetry')
      end if
#if defined (VAR_MPI)
      ! parallel version of DIRAC
C     Transfer necessary information from common blocks
      if (mytid == 0) then
         orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
         ! instruct slaves to go to the current routine (exa_moltra) as well
         call dirac_parctl( 18 )
      end if
      ! slaves get input on a need-to-know basis, this is needed to determine nmo and mo_list
      ! TODO: better is to broadcast nmo and mo_list instead 
      call interface_mpi_bcast(orbital_string,72*4,mparid,
     &                         global_communicator)
      call interface_mpi_bcast(nesh,1,mparid,global_communicator)
      if (get_talsh_only()) then
         if (mytid /= 0) return ! Needed for TALSH, keep going for exatensor
      end if
#else
      ! serial version of DIRAC
      if (get_talsh_only()) then
         ! TALSH version works fine with serial version of DIRAC
         orbital_string = tra4_indstr(:,1) ! This is to be changed once we have an exacorr-specific input driver
      else
         ! EXATENSOR version does not work with serial version of DIRAC
         print*, "Attempting to use exatensor with a non MPI version"
         call quit ("error in setup: exatensor requires MPI")
      end if
#endif

C
C     Loop over indices to convert string information to list of orbitals
C
C     Count number of orbitals to be transformed for each index 
      ntmo = 0
      do i = 1, 4
        nmo(i) = 0
        call numlst (orbital_string(i),kvec,1,1,nesh(1),1,nmo(i))
        ntmo = ntmo + nmo(i)
      end do

C     Fill the array with orbital indices
      allocate (mo_list(ntmo))
      imo = 1
      do i = 1, 4
       call numlst
     & (orbital_string(i),mo_list(imo),nmo(i),1,nesh(1),1,nmo(i))
       imo = imo + nmo(i)
      end do

C
      if (get_talsh_only()) then
         call talsh_compute_mo_integrals (nmo, mo_list)
      else
         call exacorr_compute_mo_integrals (nmo, mo_list)
      end if
C
      deallocate (mo_list)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXA_CC()
C***********************************************************************
C
C     Driver routine for exacorr CC implementation
C
C     Written by L.Visscher Jul 2017
C
C***********************************************************************

      Use exacc_cfg
      Use talsh_cc

#include "implicit.h"! because common blocks may assume this
#include "dcborb.h"  ! to access nesh (number of positive energy "electronlike" orbitals)
#include "dgroup.h"  ! to access nfsym (should be 1)

      integer              :: kvec              ! Dummy variable needed in first call to numlst
      integer, allocatable :: mo_list(:),mo_occ(:),mo_vir(:)
      integer              :: nmo(4),noccu,nvir

      if (nfsym.ne.1) then
         call quit ('exacorr can not handle symmetry')
      end if
C
      noccu = 0
      call numlst (string_occupied,kvec,1,1,nesh(1),1,noccu)
      allocate (mo_occ(noccu))
      call numlst (string_occupied,mo_occ,noccu,1,nesh(1),1,noccu)

      nvir = 0
      call numlst (string_virtual,kvec,1,1,nesh(1),1,nvir)
      allocate (mo_vir(nvir))
      call numlst (string_virtual,mo_vir,nvir,1,nesh(1),1,nvir)

      !TODO: make TALSH/EXACORR switch similar to the one in exa_moltra
      call cc_driver (noccu,nvir,mo_occ,mo_vir)

      deallocate (mo_occ)
      deallocate (mo_vir)

      end

       subroutine Initialize_Shells (shell,ierr)

!      get information from common blocks about the basis set
!      raise an error for contracted sets that can not yet be handled by Interest

!      NB: this routine is outdated, replaced by read_from_dfcoef. Just keeping it in case we still need it.

       use exacorr_datatypes

#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"

       integer ishell, iPrimStart, iPrimEnd, k, ierr
       type(basis_func_info_t) :: shell(*)

       do ishell = 1, nlrgsh
          shell(ishell)%orb_momentum = nhkt(ishell)
          shell(ishell)%atom_number  = ncent(ishell)
          shell(ishell)%n_primitives = nuco(ishell)
          nullify  (shell(ishell)%exponent)
          allocate (shell(ishell)%exponent(nuco(ishell)))
          nullify  (shell(ishell)%coefficient)
          allocate (shell(ishell)%coefficient(nuco(ishell)))
          iPrimStart = jstrt(ishell) + 1
          iPrimEnd   = jstrt(ishell) + nuco(ishell)
          k          = numcf(ishell)
         shell(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
         shell(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
         shell(ishell)%coord = cent(ishell,1:3,1)
!        In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
         call compress_unc_shell (shell(ishell),ierr)
       end do

       return 
       end

       subroutine compress_unc_shell (shell,ierr)
!If basis functions were grouped together in the input, they appear as quasi-contracted functions
!with a lot of zero coefficients. Remove these zeroes and reduce the length to just one function.
!Routine will produce an error, if the function is a true contraced function (more than one non-zero coefficient)

          use exacorr_datatypes

          integer, intent(out)     :: ierr
          type(basis_func_info_t), intent(inout) :: shell

          integer i, j, non_zero
          real(8), allocatable :: exponents(:), coefficients(:)
          real(8), parameter   :: treshold=1.D-12

          allocate(exponents(shell%n_primitives))
          allocate(coefficients(shell%n_primitives))
          exponents    = shell%exponent
          coefficients = shell%coefficient

          ierr = 0
          non_zero = 0
          do i = 1, shell%n_primitives
             if (abs(coefficients(i)) > treshold) then
                non_zero = non_zero + 1
                j = i
             end if
          end do

          if (non_zero > 1) then
             ierr = 1
             return   ! leave shell untouched and return with error condition
          else if (non_zero == 0) then
             print*, " WARNING: compress_unc_shell found empty shell"
             ierr = 2
             j = 1    ! set j to the first function, to avoid complications with zero length arrays
          end if

          deallocate (shell%exponent)
          deallocate (shell%coefficient)
          allocate   (shell%exponent(1))
          allocate   (shell%coefficient(1))

          shell%n_primitives = 1
          shell%exponent     = exponents(j)
          shell%coefficient  = coefficients(j)

          deallocate (exponents)
          deallocate (coefficients)

       end subroutine compress_unc_shell
  
       subroutine read_from_dfcoef (ao_basis,qmo,ierr)
  
       use dirac_mo
       use exacorr_datatypes
       use exacorr_global,only: ncartesian

       implicit none

!      read all positive energy mo coefficients present on DFCOEF file
!      check for imcompatibiklities (use of symmetry) as well and continue
!      by reading also the basis set information into gtos

       type(basis_set_info_t), intent(inout) :: ao_basis
       type(MO), intent(inout):: qmo
       integer, intent(out)   :: ierr 
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
       integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old DFCOEF?)
       integer, parameter     :: CONTRAC  = 5 ! basis set is contracted (Interest can not handle this)
       integer, parameter     :: INCONSIS = 6 ! DFCOEF does not contain consistent information wrt number of basis functions
       character              :: text*74      ! holds text read from DFCOEF, only used for printing this
       integer                :: nsym
       integer                :: idim(3,2)
       real(8)                :: total_energy
       logical                :: tobe
       integer                :: lucoef
       integer                :: i, j, k, ishell, iprimstart, iprimend
       integer                :: nao, nao_bas, nshells, iao, ibas
       integer                :: nmo, npo, nto, npriexp, ncnt
       integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
       real(8), allocatable   :: coeff(:,:,:),eig(:)
       integer, allocatable   :: ibeig(:),index(:)
       integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
       integer, allocatable   :: jstrt(:)
       real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)

       ierr = NOERROR

       ! open file and read dimension info
       inquire (file='DFCOEF',exist=tobe)
       if (.not.tobe) then
          print*, "error: MO coefficient file DFCOEF not found"
          ierr = NOTFOUND
          return
       else
          call get_free_fileunit(lucoef)
          open (lucoef,file='DFCOEF',status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
          read (lucoef,end=10,err=10) text,nsym,
     &         ((idim(i,j),i=1,3),j=1,nsym),total_energy
       end if

       ! check symmetry and get dimensions of the coefficient array
       if (nsym > 1) then
          ierr = SYMMETRY
          return
       else
          nao = idim(3,1) ! number of basis functions
          npo = idim(1,1) ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
          nmo = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
          nto = npo+nmo   ! total number of solutions on this file
          if (npo > 0) print*, "warning: reading from 4C DFCOEF file" ! this warning will go once we have the full index transform implementedd
       end if

       allocate(coeff(nao,nto,nz))
       allocate(eig(nto))
       allocate(ibeig(nto))
       allocate(index(nto))
        
       ! read the mo coefficients and orbital eigenvalues
       read (lucoef,end=10,err=10) coeff
       read (lucoef,end=10,err=10) eig
       read (lucoef,end=10,err=10) ibeig
       ! fill in index values
       do i = 1, nmo
          index(i) = i ! we set the index compatible with the fact that we store all mo's
       end do

       ! put in qmo
       call create_qmo (qmo,nao,nmo,nz)
       qmo%total_energy = total_energy
       qmo%index        = index
       qmo%coeff        = coeff(:,npo+1:nto,1:nz)
       qmo%energy       = eig(npo+1:nto)
       qmo%boson_irrep  = ibeig(npo+1:nto)

       deallocate(coeff)
       deallocate(eig)
       deallocate(ibeig)
       deallocate(index)

       ! read the basis set information and check for consistency
       read (LUCOEF,err=20,end=20) nao_bas,nshells,npriexp,ncnt 
       if (nao /= nao_bas) then
          ierr = INCONSIS
          return
       end if

       ! allocate the arrays to read the basis set information in DIRAC format
       allocate(ncent(nshells))
       allocate(nhkt(nshells))
       allocate(nuco(nshells))
       allocate(numcf(nshells))
       allocate(jstrt(nshells))
       allocate(cent(nshells,3))
       allocate(priexp(npriexp))
       allocate(priccf(npriexp,ncnt))

!      orb_momentum, atom_number, n_primives for all shells
       read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
!      start of exponent values, number of contracted for all shells
       read (LUCOEF,err=20,end=20) jstrt,numcf
!      coordinates (x, y, z) of center for all shells
       read (LUCOEF,err=20,end=20) cent
!      values of the exponents for all shells (packed into a 1-d array)
       read (LUCOEF,err=20,end=20) priexp
!      values of the contraction coefficients for all shells (packed into a 1-d array)
       read (LUCOEF,err=20,end=20) priccf

       ! put in gtos
       ao_basis%nshells = nshells
       ao_basis%nao     = nao
       nullify(ao_basis%gtos)
       allocate(ao_basis%gtos(nshells))
       do ishell = 1, nshells
          ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
          ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
          ao_basis%gtos(ishell)%n_primitives = nuco(ishell)
          nullify  (ao_basis%gtos(ishell)%exponent) 
          allocate (ao_basis%gtos(ishell)%exponent(nuco(ishell)))
          nullify  (ao_basis%gtos(ishell)%coefficient)
          allocate (ao_basis%gtos(ishell)%coefficient(nuco(ishell)))
          iPrimStart = jstrt(ishell) + 1
          iPrimEnd   = jstrt(ishell) + nuco(ishell)
          k          = numcf(ishell)
          ao_basis%gtos(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
          ao_basis%gtos(ishell)%coefficient = 
     &                                   priccf(iPrimStart:iPrimEnd,k)
          ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
!         In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
          call compress_unc_shell (ao_basis%gtos(ishell),ierr)
          if (ierr == 1) ierr = CONTRAC
       end do

       deallocate(ncent)
       deallocate(nhkt)
       deallocate(nuco)
       deallocate(numcf)
       deallocate(jstrt)
       deallocate(cent)
       deallocate(priexp)
       deallocate(priccf)

       ! make array with pointers from a basis function index to its shell
       nullify(ao_basis%shell_indices)
       allocate(ao_basis%shell_indices(nao))
       iao = 0
       do ishell = 1, nshells
          do ibas = 1, ncartesian(ao_basis%gtos(ishell)%orb_momentum)
             iao = iao + 1
             if (iao.gt.nao) then
                ierr = INCONSIS
                call quit ('Error in read_from_dfcoef')
             endif
             ao_basis%shell_indices(iao) = ishell
          end do
       end do
       if (iao.gt.nao) then
          ierr = INCONSIS
          call quit ('Error in read_from_dfcoef')
       endif

       close (lucoef,status='KEEP')

       call print_date("read from DFCOEF:"//text)
       return

  10   ierr = CORRUPT
       return

  20   ierr = NOBASIS
       return

       end subroutine read_from_dfcoef

       subroutine Write_To_Mdcint (nmo,rcw,mo_tens)

!      Utility routine to write to the DIRAC moltra format, makes it easier to test the index transform
 
!      Note that the integral matrix is transposed in moltra, this is why the 3 and 4th indices are written to the left below
!      This routine should therefore be rewritten if the matrix transpose I(ij,kl) .ne. I(kl,ij)

#include "implicit.h"

       integer rcw    ! this is used to toggle between complex (so, rcw=2) and real (scalar, rcw=1) tensors
       integer iclass, nclass
       integer nmo(4)
       real(8) mo_tens(rcw,nmo(3)*rcw,nmo(4)*rcw,nmo(1)*rcw,nmo(2)*rcw) ! the first rcw is real/complex, the other toggles orbitals/spinors
       real(8), allocatable :: fmat(:,:),rkl(:,:)
       integer, allocatable :: indk(:),indl(:)
       integer :: nmo_start(8,4),nmo_end(8,4)

       iprint = 0

!      Open the MDCINT file
       call mdinti (0) 

!      We use no symmetry, but may have both scalar and spin-orbit tensors as input
       irepij = 1
       if (rcw == 1) then
          nclass = 1  ! this implies spinfree
       else
          nclass = 8
       end if

!      We only need to write a subset of integrals, the rest will be generated by permutation and/or Kramers' symmetry
 
!      Classes are defined to have a standard subset of integrals, make these consistent with the definition in mdintb (in moltra)
!      Class 1 : (k    l    | i    j   )
!      Class 2 : (kbar lbar | i    j   )
!      Class 3 : (k    lbar | i    jbar)
!      Class 4 : (kbar l    | i    jbar)
!      Class 9 : (kbar l    | i    j   )
!      Class 10: (k    lbar | i    j   )
!      Class 11: (k    l    | ibar j   )
!      Class 12: (k    l    | i    jbar)

!      NB: classes 5 through 8 are not used, we therefore shift the final 4 classes in the array below to get a contiguous index

       do i = 1, 4
          do iclass = 1, 8
             nmo_start(iclass,i) = 1
             nmo_end  (iclass,i) = nmo(i)
          end do
       end do

!      define bars for class 2
       nmo_start(2,3) = nmo(3) + 1
       nmo_end  (2,3) = 2 * nmo(3)
       nmo_start(2,4) = nmo(4) + 1
       nmo_end  (2,4) = 2 * nmo(4)
!      define bars for class 3
       nmo_start(3,2) = nmo(2) + 1
       nmo_end  (3,2) = 2 * nmo(2)
       nmo_start(3,4) = nmo(4) + 1
       nmo_end  (3,4) = 2 * nmo(4)
!      define bars for class 4
       nmo_start(4,2) = nmo(2) + 1
       nmo_end  (4,2) = 2 * nmo(2)
       nmo_start(4,3) = nmo(3) + 1
       nmo_end  (4,3) = 2 * nmo(3)
!      define bars for class 5=9
       nmo_start(5,3) = nmo(3) + 1
       nmo_end  (5,3) = 2 * nmo(3)
!      define bars for class 6=10
       nmo_start(6,4) = nmo(4) + 1
       nmo_end  (6,4) = 2 * nmo(4)
!      define bars for class 7=11
       nmo_start(7,1) = nmo(1) + 1
       nmo_end  (7,1) = 2 * nmo(1)
!      define bars for class 8=12
       nmo_start(8,2) = nmo(2) + 1
       nmo_end  (8,2) = 2 * nmo(2)

!      Allocate memory needed by mdintb
       nkl = nmo(3) * nmo(4)
       allocate(rkl(rcw,nkl))
       allocate(indk(nkl))
       allocate(indl(nkl)) 
       allocate(fmat(nkl,rcw))

!      loop over classes to write out subsets of integrals
       do iclass = 1, nclass
          if (iclass.le.4) then
             jclass = iclass
          else
             jclass = iclass + 4  ! shift iclass to obey the convention used in mdintb (given above)
          end if

!         loop over blocks of integrals, taking i and j as outermost indices like is done in moltra
          do j = nmo_start(iclass,2), nmo_end(iclass,2)
             if (j.le.nmo(2)) then ! mdintb takes as index Kramers pairs not spinors
                jj = j
             else
                jj = j - nmo(2)
             endif
             do i = nmo_start(iclass,1),nmo_end(iclass,1)
               if (i.le.nmo(1)) then ! mdintb takes as index Kramers pairs not spinors
                  ii = i
               else
                  ii = i - nmo(1)
               endif

!              In moltra the real/imaginary index within a block is last rather than first, so we need to reorganize

!              If the particle/particle (index 12 <-> index 34) symmetry is not fulfilled we should reorganize further
!              this can be done here, but will be quite inefficient and not worth the trouble for this utility routine that
!              is to discarded later anyhow.
     
               rkl = reshape(mo_tens(1:rcw,
     &                       nmo_start(iclass,3):nmo_end(iclass,3),
     &                       nmo_start(iclass,4):nmo_end(iclass,4),i,j),
     &                       (/rcw, nmo(3)*nmo(4)/))
               fmat = transpose(rkl)

!              write a block of integrals to file. mdintb takes care of thresholding and labelling
               call mdintb(iprint,rcw,jclass,irepij,nmo(3),nmo(4),ii,jj,
     &                     fmat,nkl,indk,indl,rkl)

             end do
          end do
       end do

!      Close the MDCINT file
       call mdintf 
       deallocate(fmat)
       deallocate(rkl)
       deallocate(indk)
       deallocate(indl)

       return
       end
