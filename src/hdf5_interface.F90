module hdf5_interface

    private 

    public read_from_hdf5, write_mos_to_hdf5

    interface read_from_hdf5
      module procedure read_hdf5
    end interface

    interface write_mos_to_hdf5
      module procedure replace_mobasis
      module procedure write_to_checkpoint
    end interface

    interface
      module subroutine write_to_checkpoint (spinors)
        use rose_mo, only: mo_basis
        implicit none
        type(mo_basis), intent(inout):: spinors
      end subroutine write_to_checkpoint
    end interface
contains

    subroutine read_hdf5 (datafile_path,ierr,ao_basis,spinors,molecule)

       use rose_mo
       use rose_ao
       use rose_utils
       use labeled_storage
       use read_xyzfile

       implicit none

!      read information from datafile formatted like DIRACs CHECKPOINT.h5

       character(len=*),intent(in) :: datafile_path
       integer, intent(out)   :: ierr
       type(basis_set_info_t), intent(inout), optional :: ao_basis
       type(basis_set_info_t) :: ao_basis_local ! local variable in case we only need to return the mo_basis
       type(mo_basis), intent(inout), optional :: spinors
       type(molecule_info), intent(inout), optional :: molecule
       type(file_info_t)      :: datafile
       character*(:), allocatable :: ao_label, mo_label
       integer                :: natoms
       real(8), allocatable   :: geometry(:),nuc_charge(:),bas_charge(:)
       logical                :: label_present

       ! Initialize the file_info object and sections to read the AO and MO data from.
       datafile%name = trim(datafile_path)//'.h5'
       datafile%type = 2

       ! Generic labels (used in the Psi4 and PysSCF interface)
       ao_label = '/aobasis/' 
       mo_label = '/mobasis/'
       ! Check whether we have DIRAC (same data structure, but label is different) 
       call lab_query(datafile,'/input/aobasis/3',label_present)
       if (label_present) then
         ao_label = '/input/aobasis/3/'
         mo_label = '/result/wavefunctions/scf/mobasis/'
         call lab_query(datafile,'/result/wavefunctions/scf/x2c_mobasis/mobasis/nz',label_present)
         if (label_present) mo_label = '/result/wavefunctions/scf/x2c_mobasis/mobasis/'
       end if

       if (present(ao_basis)) then
         if (present(spinors)) then
            ! Read AO and MO basis information
            call read_aomobasis (datafile,ierr,ao_label,ao_basis,mo_label,spinors)
         else
            ! Read only AO basis information
            call read_aomobasis (datafile,ierr,ao_label,ao_basis)
         end if
       elseif (present(spinors)) then
            ! Read only MO basis information
            call read_aomobasis (datafile,ierr,ao_label,ao_basis_local,mo_label,spinors)
       end if

       ! Read molecule information and compute total nuclear charge (znumber)
       if (present(molecule)) then
	      call lab_read(datafile, '/input/molecule/n_atoms',idata=natoms)
	      allocate(geometry(3*natoms))
	      allocate(nuc_charge(natoms))
	      allocate(bas_charge(natoms))
	      call lab_read(datafile, '/input/molecule/nuc_charge',rdata=nuc_charge)
	      call lab_query(datafile, '/input/molecule/bas_charge', label_present)
	      if (label_present) then
	         call lab_read(datafile, '/input/molecule/bas_charge',rdata=bas_charge)
	      else
	         bas_charge = nuc_charge
	      end if
	      call lab_read(datafile, '/input/molecule/geometry',rdata=geometry)
	      geometry = geometry/xtang  ! xyz files have coordinates in Angstroms, convert to au
	      call get_atoms (natoms,reshape(geometry,(/3,natoms/)),molecule%firstatom, &
	                      nuc_charges=nuc_charge,bas_charges=bas_charge)
	      molecule%natoms  = natoms
	      molecule%znumber = sum(nuc_charge)
       end if

!      datafile is local to this routine, so close it here rather than leaving
!      the handle for HDF5 to clean up from its atexit handler.
       call lab_close (datafile)

    end subroutine read_hdf5

    subroutine read_aomobasis (datafile,ierr,ao_label,ao_basis,mo_label,spinors)

       use rose_mo
       use rose_ao
       use rose_utils
       use labeled_storage

       implicit none

!      read a mo basis from datafile formatted like DIRACs CHECKPOINT.h5
!      this subroutine is to kept in sync with the DIRAC subroutine with the same name
!      differs mainly in:
!         - using type(mo_basis) instead of DIRACs type(cmo) for storing MO information
!         - not having the argument string_active to select a subset of MOs

       type(file_info_t),intent(inout):: datafile
       character(len=*),intent(in) :: ao_label
       type(basis_set_info_t), intent(inout) :: ao_basis
       character(len=*),intent(in),optional :: mo_label
       type(mo_basis), intent(inout), optional:: spinors
       integer, intent(out)   :: ierr
       integer, parameter     :: NOERROR  = 0 ! no error
       integer, parameter     :: NOTFOUND = 1 ! file not found
       integer, parameter     :: CORRUPT  = 2 ! file unreadable
       integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
       integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old CHECKPOINT?)
       integer, parameter     :: INCONSIS = 6 ! CHECKPOINT does not contain consistent information wrt number of basis functions
       integer                :: nao, nmo, nmo_all, npo
       integer                :: nz = 4 ! we use no symmetry, so always have this
       real(8), allocatable   :: coeff(:)
       real(8), allocatable   :: coeff_all(:,:,:)
       real(8), allocatable   :: energy_all(:)
       integer, allocatable   :: irrep_all(:), shell_id_all(:)
       integer                :: n_mo(2),n_po(2),n_basis(2)
       logical                :: label_present, small_basis_present
       character*(:), allocatable :: C1

       ierr = NOERROR
       n_mo = 0
       n_po = 0
       n_basis = 0

       ! read AO basis information
       call read_ao_basis(datafile,ierr,ao_label,ao_basis,small_component=.false.)
       if (ierr /= NOERROR) return

       call check_ao_overlap(datafile,ao_basis)

       call lab_query(datafile,'/input/aobasis/4/n_ao',small_basis_present)
       if (small_basis_present) then
          call read_ao_basis(datafile,ierr,'/input/aobasis/4/',ao_basis,small_component=.true.)
          if (ierr /= NOERROR) return
       end if

       if (.not.present(spinors)) return

       ! read the dimensions of the mo basis and check for gerade/ungerade symmetry-blocking (not supported)
       call lab_read (datafile,mo_label//'n_basis',idata=n_basis)
       call lab_read (datafile,mo_label//'n_mo',idata=n_mo)
       call lab_read (datafile,mo_label//'nz',idata=nz)
       call lab_query(datafile,mo_label//'n_po',label_present)
       if (label_present) then
          call lab_read (datafile,mo_label//'n_po',idata=n_po)
       else
          n_po = 0
       end if
       if (n_mo(2) /= 0) then
          nao = sum(n_basis)
          nmo_all = sum(n_mo)
          npo = sum(n_po)
       else
          nao = n_basis(1)
          nmo_all = n_mo(1)
          npo = n_po(1)
       end if
       nmo = nmo_all - npo

       call lab_query(datafile,mo_label//'orbitals_C1',label_present)

       if (label_present) then
          C1 = '_C1' ! DIRAC, always taking coefficients without symmetry
          nz = 4
       else
          C1 = ''    ! Generic code, no extension in the label
       end if
       ! we do not know whether we have restricted or unrestricted, assume for now always restricted
       ! in that case nz equals the type (1: real (nonrelativistic) restricted, 4: quaternion (relativistic) restricted

       call alloc_mo_basis (spinors,nao,nmo,nz)

       ! read the mo coefficients
       allocate(coeff(nao*nmo_all*nz))
       call lab_read (datafile,mo_label//'orbitals'//C1,rdata=coeff)
       allocate(coeff_all(nao,nmo_all,nz))
       coeff_all = reshape(coeff,(/nao,nmo_all,nz/))
       ! put in spinors
       if (nz==1) then
          spinors%coeff_r      = coeff_all(:,npo+1:nmo_all,:)
       elseif (nz==2) then
          stop 'reading complex spinors is not yet implemented'
       elseif (nz==4) then
          spinors%coeff_q      = coeff_all(:,npo+1:nmo_all,:)
          if (npo > 0) then
             spinors%has_parent    = .true.
             spinors%parent_read   = .true.
             spinors%parent_active = .true.
             spinors%nao_large     = ao_basis%nao
             spinors%nao_small     = nao - ao_basis%nao
             allocate(spinors%coeff_parent_q(nao,nmo,4))
             spinors%coeff_parent_q = spinors%coeff_q
          end if
       end if
       deallocate(coeff_all)

       ! read the orbital eigenvalues if they are available
       call lab_query(datafile,mo_label//'eigenvalues'//C1,label_present)
       if (label_present) then
          allocate(energy_all(nmo_all))
          call lab_read (datafile,mo_label//'eigenvalues'//C1,rdata=energy_all)
          spinors%energy = energy_all(npo+1:nmo_all)
          deallocate(energy_all)
       else
          spinors%energy = 0.D0
       end if

       ! read orbital symmetry information if this is available
       call lab_query(datafile,mo_label//'symmetry',label_present)
       if (label_present) then
          allocate(irrep_all(nmo_all))
          call lab_read (datafile,mo_label//'symmetry',idata=irrep_all)
          spinors%irrep = irrep_all(npo+1:nmo_all)
          deallocate(irrep_all)
       else
          spinors%irrep        = 0 ! because nosym
       end if

       ! shell_id is only safe to import when it follows the same orbital ordering.
       ! DIRAC's C1 orbitals may be reordered relative to symmetry-blocked shell_id.
       call lab_query(datafile,mo_label//'shell_id',label_present)
       if (label_present .and. C1 .eq. '') then
          allocate(shell_id_all(nmo_all))
          call lab_read (datafile,mo_label//'shell_id',idata=shell_id_all)
          spinors%shell_id = shell_id_all(npo+1:nmo_all)
          deallocate(shell_id_all)
       end if

       call lab_query(datafile,'/result/wavefunctions/scf/energy',label_present)
       if (label_present) then
          call lab_read (datafile,'/result/wavefunctions/scf/energy',rdata=spinors%total_energy)
       else
          spinors%total_energy = 0.0
       end if

       ! check consistency with the size of the ao basis
       if (spinors%has_parent) then
          if (.not. ao_basis%has_small_component) then
             ierr = INCONSIS
             return
          end if
          if (spinors%nao_large /= ao_basis%nao) then
             ierr = INCONSIS
             return
          end if
          if (spinors%nao_small /= ao_basis%nao_small) then
             ierr = INCONSIS
             return
          end if
          if (spinors%nao /= ao_basis%nao + ao_basis%nao_small) then
             ierr = INCONSIS
             return
          end if
       elseif (nao .ne. ao_basis%nao) then
          ierr = INCONSIS
          return
       end if

    end subroutine read_aomobasis

    subroutine check_ao_overlap (datafile,ao_basis)

!      If the file carries a reference AO overlap, rebuild it from the basis
!      that was just read and report the largest deviation.
!
!      rose_kf_bridge stores that reference when it converts an ADF
!      calculation, taking it from Matrices%Smat. The comparison is the gate on
!      the conversion: a basis that reproduces ADF's own overlap is in ADF's AO
!      order, in ADF's normalisation, with the l convention read correctly, and
!      the overlap routine that produced it works on a real basis rather than
!      on the analytic cases of its self-test. One number settles all four, and
!      it keeps settling them after the ADF route is deleted, because the
!      reference travels with the fixture rather than with ADF.
!
!      The reference is the lower triangle, row by row, exactly as it sits on
!      TAPE15: element (i,j) with j <= i at i*(i-1)/2 + j.
!
!      This is a check, not a computation: nothing downstream uses the result,
!      and a file without the reference passes through untouched.

       use rose_ao
       use labeled_storage
       use slater_overlap, only: sto_overlap_matrix

       implicit none

       type(file_info_t), intent(inout)   :: datafile
       type(basis_set_info_t), intent(in) :: ao_basis

       real(8), allocatable :: reference(:), computed(:,:)
       real(8)              :: worst, scale
       integer              :: nao, i, j, idx
       logical              :: label_present
       character(len=*), parameter :: ref_label = '/validation/ao_overlap_tri'

       call lab_query (datafile,ref_label,label_present)
       if (.not. label_present) return

       if (ao_basis%basis_type /= 2 .or. .not. associated(ao_basis%stos)) then
          write (*,'(/A)') ' AO overlap check skipped: the reference '// &
                                    'is only interpretable for a Slater basis'
          return
       end if

       nao = ao_basis%nao
       allocate (reference(nao*(nao+1)/2))
       call lab_read (datafile,ref_label,rdata=reference)

       allocate (computed(nao,nao))
       call sto_overlap_matrix (ao_basis%stos,ao_basis%stos,computed, &
                                angular=ao_basis%basis_angular)

       worst = 0.0d0
       scale = 0.0d0
       idx   = 0
       do i = 1, nao
          do j = 1, i
             idx   = idx + 1
             worst = max(worst,abs(computed(i,j)-reference(idx)))
             scale = max(scale,abs(reference(idx)))
          end do
       end do

       write (*,'(/A)')      ' AO overlap check against the reference on file'
       write (*,'(A,I8)')    '   dimension                      ',nao
       write (*,'(A,ES12.4)') '   largest absolute deviation    ',worst
       write (*,'(A,ES12.4)') '   largest reference element     ',scale
       if (worst > 1.0d-6) then
          write (*,'(A)') '   **WARNING** the computed overlap does not match the'
          write (*,'(A)') '   reference. The AO ordering, the normalisation or the'
          write (*,'(A)') '   angular convention of the converted basis is wrong.'
       end if

       deallocate (reference,computed)

    end subroutine check_ao_overlap


    subroutine read_ao_basis(datafile,ierr,ao_label,ao_basis,small_component)

       use rose_ao
       use rose_utils
       use labeled_storage

       implicit none

       type(file_info_t),intent(inout) :: datafile
       integer, intent(out)            :: ierr
       character(len=*),intent(in)     :: ao_label
       type(basis_set_info_t), intent(inout) :: ao_basis
       logical, intent(in)             :: small_component
       type(gaussian_func_info_t), pointer :: gtos(:)
       type(slater_func_info_t), pointer   :: stos(:)
       integer, pointer                :: shell_indices(:)
       integer                         :: nao, nshells, basis_angular, basis_type
       integer, parameter              :: NOERROR  = 0
       integer, parameter              :: CONTRAC  = 5
       integer, parameter              :: INCONSIS = 6
       integer                         :: ishell, iprimstart, iprimend
       integer                         :: iao, ibas, npriexp
       integer, allocatable            :: orbmom(:), n_prim(:), n_cont(:), prinqn(:)
       real(8), allocatable            :: priexp(:), priccf(:), coord(:)
       logical                         :: label_present

       ierr = NOERROR

       call lab_read (datafile,ao_label//'angular',idata=basis_angular)
       call lab_read (datafile,ao_label//'n_ao',idata=nao)
       call lab_read (datafile,ao_label//'n_shells',idata=nshells)

!      The kind of basis function. No file written before Slater support existed
!      carries this label, so its absence means Gaussian and every such file keeps
!      reading exactly as before.
       basis_type = 1
       call lab_query (datafile,ao_label//'basis_type',label_present)
       if (label_present) call lab_read (datafile,ao_label//'basis_type',idata=basis_type)

       allocate (orbmom(nshells))
       allocate (n_prim(nshells))
       allocate (n_cont(nshells))
       allocate (coord(3*nshells))
!      read information that is defined for each shell
       call lab_read(datafile,ao_label//'orbmom',idata=orbmom)
       call lab_read(datafile,ao_label//'n_prim',idata=n_prim)
       call lab_read(datafile,ao_label//'n_cont',idata=n_cont)
       call lab_read(datafile,ao_label//'center',rdata=coord)
!      values of the exponents and coefficients for all shells (packed into a 1-d arrays)
       npriexp = sum(n_prim)
       allocate (priexp(npriexp))
       allocate (priccf(npriexp))
       call lab_read(datafile,ao_label//'exponents',rdata=priexp)
       call lab_read(datafile,ao_label//'contractions',rdata=priccf)
!      a Slater primitive carries a principal quantum number alongside its exponent
       if (basis_type == 2) then
          allocate (prinqn(npriexp))
          call lab_read(datafile,ao_label//'n',idata=prinqn)
       end if

       nullify(gtos)
       nullify(stos)
       select case (basis_type)

       case (1)

          allocate(gtos(nshells))
          iPrimStart = 1
          do ishell = 1, nshells
             gtos(ishell)%orb_momentum = orbmom(ishell)
             gtos(ishell)%n_primitives = n_prim(ishell)
             nullify  (gtos(ishell)%exponent)
             allocate (gtos(ishell)%exponent(n_prim(ishell)))
             nullify  (gtos(ishell)%coefficient)
             allocate (gtos(ishell)%coefficient(n_prim(ishell)))
             iPrimEnd   = iPrimStart+n_prim(ishell)-1
             gtos(ishell)%exponent   = priexp(iPrimStart:iPrimEnd)
             gtos(ishell)%coefficient= priccf(iPrimStart:iPrimEnd)
             gtos(ishell)%coord = coord(ishell*3-2:3*ishell)
!            In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
             call compress_unc_shell (gtos(ishell),ierr)
             if (ierr == 1) ierr = CONTRAC
             iPrimStart = iPrimStart + n_prim(ishell)
          end do

       case (2)

!         Slater shells. There is no quasi-contraction to undo here: that is a DIRAC
!         convention for Gaussians, so compress_unc_shell has no counterpart.
          if (small_component) stop 'read_ao_basis: a Slater basis has no small component'
          allocate(stos(nshells))
          iPrimStart = 1
          do ishell = 1, nshells
             stos(ishell)%orb_momentum = orbmom(ishell)
             stos(ishell)%n_primitives = n_prim(ishell)
             allocate (stos(ishell)%n(n_prim(ishell)))
             allocate (stos(ishell)%exponent(n_prim(ishell)))
             allocate (stos(ishell)%coefficient(n_prim(ishell)))
             iPrimEnd   = iPrimStart+n_prim(ishell)-1
             stos(ishell)%n           = prinqn(iPrimStart:iPrimEnd)
             stos(ishell)%exponent    = priexp(iPrimStart:iPrimEnd)
             stos(ishell)%coefficient = priccf(iPrimStart:iPrimEnd)
             stos(ishell)%coord = coord(ishell*3-2:3*ishell)
             iPrimStart = iPrimStart + n_prim(ishell)
          end do

       case default

          write(*,*) '**ERROR** unknown basis_type on file: ',basis_type
          stop 'read_ao_basis: unknown basis_type'

       end select

!      The shell to AO map only needs the angular momenta, which compress_unc_shell
!      leaves alone, so orbmom serves both kinds of basis.
       nullify(shell_indices)
       allocate(shell_indices(nao))
       iao = 0
       do ishell = 1, nshells
          do ibas = 1, nfunctions(orbmom(ishell),basis_angular)
             iao = iao + 1
             if (iao.gt.nao) then
                ierr = INCONSIS
                return
             endif
             shell_indices(iao) = ishell
          end do
       end do

       deallocate(orbmom)
       deallocate(n_prim)
       deallocate(n_cont)
       deallocate(coord)
       deallocate(priexp)
       deallocate(priccf)
       if (allocated(prinqn)) deallocate(prinqn)

       if (small_component) then
          ao_basis%has_small_component = .true.
          ao_basis%nao_small = nao
          ao_basis%nshells_small = nshells
          ao_basis%basis_angular_small = basis_angular
          ao_basis%gtos_small => gtos
          ao_basis%shell_indices_small => shell_indices
       else
          ao_basis%nao           = nao
          ao_basis%nshells       = nshells
          ao_basis%basis_angular = basis_angular
          ao_basis%basis_type    = basis_type
          ao_basis%has_small_component = .false.
          ao_basis%nao_small = 0
          ao_basis%nshells_small = 0
          ao_basis%basis_angular_small = 1
          ao_basis%gtos => gtos
          ao_basis%stos => stos
          ao_basis%shell_indices => shell_indices
       end if

    end subroutine read_ao_basis

    subroutine replace_mobasis (spinors,datafile_path1,datafile_path2, &
                                ribo_occ_per_frag,ribo_vir_per_frag,RIBO_Fock,ncore)

       use rose_mo
       use labeled_storage

       implicit none

!      replace mos on a checkpoint file by localized ones
!      should ideally be a new section but this requires a change of DIRAC schema

       character(len=*),intent(in) :: datafile_path1,datafile_path2
       type(mo_basis), intent(inout):: spinors
!      The RIBO metadata, optional because only the ribo call supplies it. Until
!      now these four quantities reached the ADF and AMS writers only, so a
!      calculation whose results live in h5 could not be written back to an rkf.
!      Written into /mobasis/ of the standalone schema only; see below.
       integer, intent(in), optional :: ribo_occ_per_frag(:)  ! per recanonisation group
       integer, intent(in), optional :: ribo_vir_per_frag(:)  ! idem, not per fragment
       real(8), intent(in), optional :: RIBO_Fock(:,:,:)
       integer, intent(in), optional :: ncore
       type(file_info_t)      :: datafile
       character*(:), allocatable :: mo_label
       integer                :: nz, nribo
       logical                :: label_present

!      first make the copy of the hdf5 file
       call execute_command_line ('cp '//trim(datafile_path1)//'.h5 '//trim(datafile_path2)//'.h5')

       ! Initialize the file_info object and sections to read the AO and MO data from.
       datafile%name = trim(datafile_path2)//'.h5'
       datafile%type = 2
       datafile%status = 0
       mo_label = '/mobasis/'
       ! Check whether we have DIRAC 
       call lab_query(datafile,'/result/wavefunctions/scf/mobasis/',label_present)
       if (label_present) then
!        switch to Kramers-restricted format
         nz = 4
         call convert_mo_basis(spinors,nz)
         mo_label = '/result/wavefunctions/scf/mobasis/'
         ! write the mo coefficients and eigenvalues (will not work when symmetry is used, need better labeling in DIRAC)
         call lab_write (datafile,mo_label//'orbitals',rdata=reshape(spinors%coeff_q,(/spinors%nao*spinors%nmo*nz/)))
         call lab_write (datafile,mo_label//'orbitals_C1',rdata=reshape(spinors%coeff_q,(/spinors%nao*spinors%nmo*nz/)))
         call lab_write (datafile,mo_label//'eigenvalues',rdata=spinors%energy)
         call lab_write (datafile,mo_label//'eigenvalues_C1',rdata=spinors%energy)
         if (present(ribo_occ_per_frag) .or. present(RIBO_Fock)) &
            stop 'replace_mobasis: the RIBO metadata is not written to a DIRAC checkpoint'
       else
          select case (spinors%type)
            case (1,2) ! real (un)restricted
               nz = 1
               call lab_write (datafile,mo_label//'orbitals',rdata=reshape(spinors%coeff_r,(/spinors%nao*spinors%nmo/)))
            case (3) ! complex unrestricted
               nz = 2
               stop 'writing complex spinors is not yet implemented'
            case (4) ! quaternion restricted
               nz = 4
               call lab_write (datafile,mo_label//'orbitals',rdata=reshape(spinors%coeff_q,(/spinors%nao*spinors%nmo*nz/)))
          end select
          ! write the orbital energies
          call lab_write (datafile,mo_label//'eigenvalues',rdata=spinors%energy)

          ! write dimension info (especially nmo may have changed)
          call lab_write (datafile,mo_label//'n_basis',idata=spinors%nao)
          call lab_write (datafile,mo_label//'n_mo',idata=spinors%nmo)
          call lab_write (datafile,mo_label//'nz',idata=nz)

!         The RIBO metadata. These are the four quantities that belong in the
!         Rose Data section of an rkf, which is what ADF's TDDFT consumes;
!         write_ROSE_results puts them there inside AMS, and rose_kf_bridge.py
!         does it from outside. Storing them here is what makes an h5 result a
!         complete record of the calculation, and so what makes that conversion
!         possible at all.
!
!         They live in /mobasis/, which they belong to and which already exists:
!         create_group_hdf5 deliberately cannot create a top-level group, so
!         every label written here has to land in a group the producing code
!         already made. This is the schema pyscf, psi4 and the KF converter
!         share. DIRAC's checkpoint has a layout of its own and a small
!         component besides, and is left alone -- it does not reach this branch,
!         and the branch above refuses the arguments rather than guessing.
!
!         Counts are per recanonisation group, not per fragment; the two coincide
!         only when no grouping is used, as the AMS README records.
!
!         ribo_fock is stored with the SECOND index running fastest, which is the
!         order the Rose Data section uses, so a converter can copy it straight
!         across. The matrix is symmetric, so the choice is invisible in practice
!         and is made only to keep the two paths comparable. Its
!         dimension is spinors%nmo - ncore: the frozen orbitals have been
!         re-inserted into spinors by this point but take no part in it.
          if (present(ribo_occ_per_frag) .and. present(ribo_vir_per_frag)) then
             if (size(ribo_occ_per_frag) /= size(ribo_vir_per_frag)) &
                stop 'replace_mobasis: occupied and virtual group counts differ in length'
             call lab_write (datafile,'/mobasis/ribo_nfrag',idata=size(ribo_occ_per_frag))
             call lab_write (datafile,'/mobasis/ribo_occ_per_frag',idata=ribo_occ_per_frag)
             call lab_write (datafile,'/mobasis/ribo_vir_per_frag',idata=ribo_vir_per_frag)
          end if
          if (present(ncore)) call lab_write (datafile,'/mobasis/ribo_n_core',idata=ncore)
          if (present(RIBO_Fock)) then
             nribo = size(RIBO_Fock,1)
             if (size(RIBO_Fock,2) /= nribo) stop 'replace_mobasis: RIBO_Fock is not square'
             if (present(ncore)) then
                if (nribo /= spinors%nmo - ncore) &
                   stop 'replace_mobasis: RIBO_Fock dimension is not nmo - ncore'
             end if
             call lab_write (datafile,'/mobasis/ribo_n_active',idata=nribo)
             call lab_write (datafile,'/mobasis/ribo_fock', &
                             rdata=reshape(transpose(RIBO_Fock(:,:,1)),(/nribo*nribo/)))
          end if
       end if

       call lab_close (datafile)

    end subroutine replace_mobasis

       subroutine compress_unc_shell (shell,ierr)
!If basis functions were grouped together in the input, they appear as quasi-contracted functions
!with a lot of zero coefficients. Remove these zeroes and reduce the length to just one function.
!Routine will produce an error, if the function is a true contraced function (more than one non-zero coefficient)

          use rose_ao

          integer, intent(out)     :: ierr
          type(gaussian_func_info_t), intent(inout) :: shell

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
             return   ! leave shell untouched and return
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

end module hdf5_interface
