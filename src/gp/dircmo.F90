
! Module to handle file writing (modernized versions of WRICMO and REACMO)
! Defined as a module to able to use interfaces and optional arguments

! Written by Lucas Visscher, spring 2021, one of the Covid years

module dircmo
  implicit none

  private

  public wricmo
  public reacmo_new

  interface wricmo
      module procedure wricmo
  end interface

  interface reacmo_new
      module procedure reacmo
  end interface

  contains

      SUBROUTINE WRICMO(IUNIT,CMO,EIG,IBEIG,TOTERG,TITLE,IDIMS)
!*****************************************************************************
!
!     Write to file:
!       CMO    - coefficients
!       EIG    - eigenvalues
!       IBEIG  - supersymmetry information
!       TOTERG - total energy
!       IDIMS  - dimension information (simple integer array to facilitate use
!                in ancient Fortran pieces of the code)
!
!     SSYM = T activates writing of second array of supersymmetry information
!
!     IDIM(1,J) - number of positronic solutions
!     IDIM(2,J) - number of electronic solutions
!     IDIM(3,J) - number of AO-basis functions
!
!     Written by T.Saue Sept 1 1995
!     Modernized by L.Visscher in 2021
!
!*****************************************************************************
      use labeled_storage
      use checkpoint
      implicit none
      integer, intent(in) :: iunit,ibeig(*)
      real(8), intent(in) :: toterg,cmo(*),eig(*)
      character*(*),intent(in), optional :: title
      integer,intent(in), optional :: idims(:)

      INTEGER   CMO_DIMS(12),IDIM(3,2),NCMOTQ,NORBT,NZ,NFSYM,NSSYM
      integer   ibeig_dim, ifsym, ierr
      integer   n_basis(2),n_mo(2),n_po(2)
      LOGICAL   SSYM
      character(len=74) text

!     Check whether we should get some information from the common blocks
      if (present(title)) then
         text(1:50) = title
      else
         call get_title (text)
      end if
      CALL GTINFO(TEXT(51:74)) ! Add date & time info to title

      if (present(idims)) then
         cmo_dims = idims
      else
         call get_dimension_info(cmo_dims)
      end if

!     Unpack the information stored in cmo_dims
      NZ     = cmo_dims(1) ! Number of active quaternion units (1, 2 or 4)
      NFSYM  = cmo_dims(2) ! Number of fermion symmetries (no inversion 1, inversion 2)
      NCMOTQ = cmo_dims(3) ! Size of MO coefficient array
      NORBT  = cmo_dims(5) ! Total number of orbitals
      NSSYM  = cmo_dims(6) ! Number of blocks with supersymmetry information
      ! unpack NPSH, NESH and NFBAS(:,0) for each fermion irrep
      IDIM(:,1) = cmo_dims(7:9)
      if (nfsym == 2) IDIM(:,2) = cmo_dims(10:12)

      SSYM   = NSSYM == 2

      !REWIND IUNIT
      !CALL LAB_WRITE(IUNIT,'INFO    ')
      !WRITE(IUNIT) TEXT,NFSYM,NZ,IDIM(:,1:NFSYM),TOTERG
      !CALL LAB_WRITE(IUNIT,'COEFS   ',CMO(1:NCMOTQ))
      !CALL LAB_WRITE(IUNIT,'EVALS   ',EIG(1:NORBT))
      !CALL LAB_WRITE(IUNIT,'SUPERSYM',idata=IBEIG(1:NORBT))
!     !Add additional information about the basis set
      !CALL LAB_WRITE(IUNIT,'BASINFO ')
!     !This call introduces a dependency on common blocks, should be moved as soon as the checkpoint file is fully functional
      !call append_basisinfo_to_dfcoef (IUNIT)
      !CALL LAB_WRITE(IUNIT,'EOFLABEL')

      ! write all data also to an hdf5 file
      ! simplify the handling of fermion symmetry by always writing the same size dimension arrays (just padding with zeroes if nfsym==1)
      n_mo = 0
      n_po = 0
      n_basis = 0
      do ifsym = 1, nfsym
         n_mo(ifsym)    = cmo_dims(4+3*ifsym)+cmo_dims(5+3*ifsym)
         n_po(ifsym)    = cmo_dims(4+3*ifsym)
         n_basis(ifsym) = cmo_dims(6+3*ifsym)
      end do
      call checkpoint_write ('/result/wavefunctions/scf/energy',rdata=toterg)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/mobasis_id',idata=1)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/eigenvalues',rdata=eig(1:norbt))
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/nz',idata=nz)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/n_fsym',idata=nfsym)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/n_basis',idata=n_basis)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/n_mo',idata=n_mo)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/n_po',idata=n_po)
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/symmetry',idata=ibeig(1:norbt))
      call checkpoint_write ('/result/wavefunctions/scf/mobasis/orbitals',rdata=cmo(1:ncmotq))

      RETURN
      END

      SUBROUTINE REACMO(TEXT,CMO,EIG,IBEIG,TOTERG,IDIMS)
!***********************************************************************
!
!     Read SCF-coefficients and other information (eigenvalues) from CHECKPOINT file
!
!     Written by T.Saue Sept 1 1995
!     Modernized by L.Visscher in 2021
!
!***********************************************************************
      use checkpoint
      implicit none
      character(len=74), intent(out), optional :: text
      integer, intent(out), optional :: ibeig(*),idims(*)
      real(8), intent(out), optional :: toterg,cmo(*),eig(*)

      integer              :: nfsym, nz, idim(3,2), ncoef, neig
      integer              :: n_mo(2),n_po(2),n_basis(2)
      real(8)              :: energy
      real(8), allocatable :: coefficients(:), eigenvalues(:)
      integer, allocatable :: boson_irreps(:), atom_irreps(:)
      character(len=74)    :: infotext
      character*(:), allocatable :: group, mobasis
      logical              :: natural_orbitals, x2c_mobasis

      n_mo = 0
      n_po = 0
      n_basis = 0
      infotext = ' Data read from the checkpoint file '

      call checkpoint_query('/result/wavefunctions/mbpt/mobasis/nz',natural_orbitals)
      call checkpoint_query('/result/wavefunctions/scf/x2c_mobasis/mobasis/nz',x2c_mobasis)
      if (natural_orbitals) then ! This is a hack as we currently have only two options, needs to be generalized !
         mobasis = '/result/wavefunctions/mbpt/mobasis/'
         group = '/result/wavefunctions/mbpt/'
      else if (x2c_mobasis) then
         mobasis = '/result/wavefunctions/scf/x2c_mobasis/mobasis/'
         group = '/result/wavefunctions/scf/'
         infotext = ' X2C-transformed data read from the checkpoint file '
      else
         mobasis = '/result/wavefunctions/scf/mobasis/'
         group = '/result/wavefunctions/scf/'
      endif

      call checkpoint_read (mobasis//'nz',idata=nz)
      call checkpoint_read (mobasis//'n_fsym',idata=nfsym)
      call checkpoint_read (mobasis//'n_basis',idata=n_basis)
      call checkpoint_read (mobasis//'n_mo',idata=n_mo)
      call checkpoint_read (mobasis//'n_po',idata=n_po)
      call checkpoint_read (group//'energy',rdata=energy)

      neig  = sum(n_mo)
      ncoef = nz * sum(n_basis*n_mo)
      idim(1,:) = n_po
      idim(2,:) = n_mo-n_po
      idim(3,:) = n_basis

      if (present(text))   text   = infotext
      if (present(toterg)) toterg = energy

      if (present(idims))  then
       ! Pack dimension information (can be used to allocate and/or check consistency before reading)
       ! NB: in the old format no information on sym was stored on dfcoef, this is necessarily
       ! constrained to be 1 to keep backwards compatibility.
         idims(1) = nz        ! Number of active quaternion units (1, 2 or 4)
         idims(2) = nfsym     ! Number of fermion symmetries (no inversion 1, inversion 2)
         idims(3) = ncoef     ! Size of MO coefficient array
         idims(4) = 1         ! Used to be NKRBLK, now obsolete.
         idims(5) = neig      ! Total number of orbitals
         idims(6) = 1         ! Number of blocks with supersymmetry information (atomic symmetry unsupported, see above)
         idims(7:9)   = idim(:,1) ! npsh, nesh, nfbas for fermion irrep 1
         idims(10:12) = idim(:,2) ! npsh, nssh, nfbas for fermion irrep 2 (zero if nfsym is 1)
      end if

      if (present(cmo)) then
         allocate (coefficients(ncoef))
         call checkpoint_read (mobasis//'orbitals',rdata=coefficients)
         cmo(1:ncoef) = coefficients
         deallocate(coefficients)
      end if

      if (present(eig)) then
         allocate (eigenvalues(neig))
         call checkpoint_read (mobasis//'eigenvalues',rdata=eigenvalues)
         eig(1:neig) = eigenvalues
         deallocate(eigenvalues)
      end if

      if (present(ibeig)) then
         allocate (boson_irreps(neig))
         call checkpoint_read (mobasis//'symmetry',idata=boson_irreps)
         ibeig(1:neig) = boson_irreps
         deallocate(boson_irreps)
      end if

      END

end module
