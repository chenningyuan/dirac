 module rose_intrinsic_mos

   private
   public :: rose_count_mospace
   public :: rose_count_occupied
   public :: rose_split_mospace
   public :: generate_iao
   public :: rose_generate_avas
   public :: generate_ibos
   public :: generate_fragvirs
   public :: copy_original_occ
   public :: copy_original_vir
   public :: rose_write_mos
   public :: rose_merge_virtuals

   interface rose_count_mospace
      module procedure rose_count_mospace
   end interface

   interface rose_count_occupied
      module procedure rose_count_occupied
   end interface

   interface rose_split_mospace
      module procedure rose_split_mospace
   end interface

   interface generate_iao
      module procedure generate_iao
   end interface

   interface generate_ibos
      module procedure generate_ibos
   end interface

   interface generate_fragvirs
      module procedure generate_fragvirs
   end interface

   interface copy_original_occ
      module procedure copy_original_occ
   end interface

   interface copy_original_vir
      module procedure copy_original_vir
   end interface

   interface rose_write_mos
      module procedure rose_write_mos
   end interface

   interface rose_merge_virtuals
      module procedure rose_merge_virtuals
   end interface

contains

   subroutine rose_count_occupied (molecule,nocc)
        use read_xyzfile
        use rose_input_data

        implicit none

        type(molecule_info),intent(in)     :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        integer,intent(out)                :: nocc ! number of occupied spinors / orbitals

        nocc = molecule%znumber - charge
        if (restricted) nocc = nocc/2

   end subroutine rose_count_occupied

   subroutine rose_count_mospace (molecule,fragments,spinor_original,nocc_active,ncore,nmo_frag)
        use rose_mo
        use read_xyzfile
        use rose_input_data
        use rose_global

        implicit none

        type(molecule_info),intent(in)     :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info),intent(in)     :: fragments(:) ! same information for the fragments
        type(mo_basis),intent(in)          :: spinor_original      ! Canonical MO-coefficients
        integer,intent(out)                :: nocc_active ! number of occupied spinors / orbitals in the active MO space
        integer,intent(out)                :: ncore ! number of core spinors / orbitals
        integer, allocatable, intent(out)  :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        integer                            :: nao1,nmo1,nmo2
        integer                            :: i
        integer, allocatable               :: nco_frag(:) ! Number of core orbitals for each fragment.
        
        ! Initialize dimensions for the supermolecule
        nao1 = spinor_original%nao
        nmo1 = spinor_original%nmo
        nocc_active = molecule%znumber - charge
        if (restricted) nocc_active = nocc_active/2
 
        ! Figure out how many fragment spinors are to be kept
        allocate (nco_frag(n_fragments))
        allocate (nmo_frag(n_fragments))       
        do i = 1, n_fragments
           nmo_frag(i) = fragments(i)%nvale
           nco_frag(i) = fragments(i)%ncore
           if (restricted) nco_frag(i) = nco_frag(i) / 2
           if (restricted) nmo_frag(i) = nmo_frag(i) / 2
        end do

        ncore = sum(nco_frag)
        nmo1 = nmo1 - ncore   ! Number of supermolecular orbitals to manipulate
        nmo2 = sum(nmo_frag)  ! Number of valence fragment orbitals
        if (nmo1 .lt. nmo2) stop "More fragment valence than supermolecular orbitals, can not proceed"
        nocc_active = nocc_active - ncore   ! Adjust number of occupied orbitals after the removal of the core

        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#         Basis set size           #"
        write(luout,*) "#               and                #"
        write(luout,*) "#      Canonical MO Energies       #"
        write(luout,*) "####################################"
        write(luout,*)
        write(luout,'(A,I6)') " No. basis functions    B1: ",nao1
        write(luout,'(A,I6)') " No. molecular orbitals B1: ",nmo1
        write(luout,'(A,I6)') " No. core      orbitals B1: ",ncore
        write(luout,'(A,I6)') " No. occupied  orbitals B1: ",nocc_active
        write(luout,'(A,I6)') " No. molecular orbitals B2: ",nmo2
        write(luout,'(A,999I6)') " No. core frag orbitals   : ",nco_frag
        write(luout,'(A,999I6)') " No. MO   frag orbitals   : ",nmo_frag
        write(luout,*)

   end subroutine rose_count_mospace

   subroutine rose_split_mospace (nfrozen,spinor_original,spinor_frozen)
        use rose_mo

        implicit none

        integer,intent(in)                 :: nfrozen
        type(mo_basis),intent(inout)       :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis),intent(out)         :: spinor_frozen ! Frozen MO-coefficients (left untouched)
        integer                            :: nmo1

        ! Split off the leading frozen orbitals, so that we can manipulate or just re-insert these later.
        call copy_mo_basis (spinor_original,spinor_frozen)
        nmo1 = spinor_original%nmo - nfrozen
        call reduce_mo_basis (spinor_frozen,nfrozen)
        call reduce_mo_basis (spinor_original,nmo1,nfrozen+1)
        spinor_frozen%shell_id = -1
        spinor_original%shell_id = 0

   end subroutine rose_split_mospace

   subroutine generate_iao (nocc,nmo_frag,fragments,ao_basis_1,spinor_original, &
                            IAO,ICO,Cocc_IAO,Cvir_IAO)
        use rose_ao
        use rose_mo
        use read_xyzfile
        use rose_import_fragments
        use properties
        use rose_input_data
        use rose_global

        implicit none

        integer,intent(in)                 :: nocc ! number of occupied spinors / orbitals
        integer, intent(in)                :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        type(molecule_info),intent(in)     :: fragments(:) ! molecule information for the fragments
        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(in)          :: spinor_original      ! Canonical MO-coefficients
        real(8), allocatable,intent(out)   :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), allocatable,intent(out)   :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable,intent(out)   :: Cocc_IAO(:,:,:) ! occupied MOs in terms of IAOs
        real(8), allocatable,intent(out)   :: Cvir_IAO(:,:,:) ! valence virtual MOs in terms of IAOs

        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis), allocatable        :: spinor_fragment(:) ! MO-coefficients of the fragments
        integer                            :: i
        integer                            :: nmo1,nmo2,nmo3

        ! Initialize dimensions
        nmo1 = spinor_original%nmo  ! Number of original orbitals / spinors
        nmo2 = sum(nmo_frag)        ! Number of valence fragment orbitals / spinors
        nmo3 = nmo1 - nmo2          ! Number of complementary orbitals / spinors

        ! Read the MOs from the fragments
        allocate (ao_basis_2(n_fragments))
        allocate (spinor_fragment(n_fragments))
        call rose_read_fragments (1,fragments,ao_basis_2,spinor_fragment)

        ! Allocate arrays to store transformation coefficients to intrinsic and complementary orbitals
        allocate (IAO(nmo1,nmo2,rcq))
        allocate (ICO(nmo1,nmo3,rcq))

        call get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,nocc,IAO,ICO)

        ! Calculate partial charges based on the localized occupied MOs expressed in terms of IAOs.
        call express_in_intrinsic_mos (IAO,nocc,Cocc_IAO,Cvir_IAO)
        call molecular_partial_charges(Cocc_IAO,spinor_fragment,fragments,restricted) 

        do i = 1, n_fragments
           call delete_adf_objects(ao_basis_2(i))
        end do

   end subroutine generate_iao

   subroutine rose_generate_avas (ao_basis_1,spinor_original,spinor_frozen,IAO,nmo_frag,nocc)
        use rose_ao
        use rose_mo
        use avas
        use hdf5_interface
        use fchk_interface
        use rose_input_data

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(in)         :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis),intent(in)         :: spinor_frozen ! Frozen MO-coefficients (left untouched)
        real(8), intent(in)               :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        integer,intent(in)                :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        integer,intent(in)                :: nocc ! number of occupied spinors / orbitals

        type(mo_basis)                    :: spinor_avas ! MO coefficients for the AVAS approach to define an active space
        type(mo_basis)                    :: spinor_iao ! MO-coefficients for the iao basis

        call IMO_to_mobasis(IAO,spinor_original,spinor_iao)
        call get_iao_avas (ao_basis_1,spinor_original,spinor_iao,nmo_frag,nocc,spinor_avas)
        call insert_mo_basis (spinor_avas,spinor_frozen)
        if (trim(file_for_MOs) .eq. "dirac") then
           call write_mos_to_hdf5(spinor_avas,"CHECKPOINT","avas")
        else if (trim(file_for_MOs) .eq. "h5") then
           call write_mos_to_hdf5(spinor_avas,"CHECKPOINT","avas")
        else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
           call write_mos_to_fchk(spinor_avas,"CHECKPOINT","avas")
        else
           stop "**ERROR** Wrong extension, I cannot write the coefficients into file."
        end if

   end subroutine rose_generate_avas

   subroutine generate_ibos (first_mo,CMO_IAO,IAO,energies,nmo_frag, &
                             LMO,energies_lmo,nlo_frag)
        use linear_algebra
        use make_ibo
        use properties
        use rose_input_data
        use rose_global

        implicit none

        integer,intent(in)                 :: first_mo ! first canonical orbital in the subspace to localize
        integer,intent(in)                 :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        real(8), intent(in)                :: CMO_IAO(:,:,:) ! MO coefficient matrix in terms of IAOs
        real(8), intent(in)                :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), intent(in)                :: energies(:) ! Energies of the core-stripped canonical orbitals
        real(8), allocatable,intent(out)   :: LMO(:,:,:)! Localized MOs in terms of original MOs of B1
        real(8), allocatable,intent(out)   :: energies_lmo(:) ! Localized/recanonicalized energies
        integer, allocatable,intent(out)   :: nlo_frag(:)! Number of localized orbitals assigned to each fragment

        integer                            :: nmo1,nmo2,nlo
        real(8), allocatable               :: CMO_IBO(:,:,:)  ! Localized MOs in terms of IAOs.

        ! Look up dimensions
        nmo1 = size(IAO,1)
        nmo2 = size(IAO,2)
        nlo  = size(CMO_IAO,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Calculate the intrinsic bonding orbital (IBO) basis                                        !
!         Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !-----Compute IBOs (Localized MO in terms of IAOs)----- (nmo2 x nlo)
        allocate(CMO_IBO(nmo2,nlo,rcq))
        !-----Compute IBOs (Localized MO in terms of original MOs)----- (nmo1 x nlo)
        allocate(LMO(nmo1,nlo,rcq))
        allocate(energies_lmo(nlo))
        write(luout,*)
        write(luout,*)
        write(luout,*) "###########################"
        write(luout,*) "# Intrinsic Localized MOs #"
        write(luout,*) "###########################"
        write(luout,*)
        call print_orbital_delocalization(CMO_IAO,nmo_frag,"before localization",&
                                          energies(first_mo:first_mo+nlo-1))
        call Localization (CMO_IAO,CMO_IBO,nmo_frag,Expnt)

        ! Express the localized MOs in the CMO basis (LMO) and sort them by ascending order:
        call matrix_mul (IAO,CMO_IBO,LMO)
        call sort_LMO_by_energy(LMO,energies,energies_lmo,restricted)
        call print_orbital_delocalization(LMO,IAO,energies_lmo,nmo_frag,"after localization")

        ! Sort localized MOs by fragments
        allocate(nlo_frag(n_fragments))
        call sort_on_fragments (CMO_IBO,nmo_frag,nlo_frag,frag_bias)
        ! Express the localized MOs in the CMO basis.
        call matrix_mul (IAO,CMO_IBO,LMO)

        ! Recanonization: block diagonalize the Fock matrix within each fragment.
        call recanonize (LMO,energies,energies_lmo,nlo_frag)
        call print_orbital_delocalization(LMO,IAO,energies_lmo,nmo_frag,"after recanonization")

        deallocate(CMO_IBO)

   end subroutine generate_ibos

   subroutine copy_original_occ (nocc,energies,Locc,energies_occ,nocc_frag)
        use rose_input_data

        implicit none

        integer,intent(in)                 :: nocc
        real(8),intent(in)                 :: energies(:)
        real(8),allocatable,intent(out)    :: Locc(:,:,:)
        real(8),allocatable,intent(out)    :: energies_occ(:)
        integer,allocatable,intent(out)    :: nocc_frag(:)
        integer                            :: i,nmo1

        nmo1 = size(energies)
        allocate(Locc(nmo1,nocc,rcq))
        allocate(energies_occ(nocc))
        allocate(nocc_frag(n_fragments))
        Locc = 0.0D0
        do i = 1, nocc
           Locc(i,i,1) = 1.0D0
        end do
        energies_occ = energies(1:nocc)
        nocc_frag = 0
        if (n_fragments > 0) nocc_frag(1) = nocc

   end subroutine copy_original_occ

   subroutine copy_original_vir (nocc,energies,Lvir,energies_vir,nvir_frag)
        use rose_input_data

        implicit none

        integer,intent(in)                 :: nocc
        real(8),intent(in)                 :: energies(:)
        real(8),allocatable,intent(out)    :: Lvir(:,:,:)
        real(8),allocatable,intent(out)    :: energies_vir(:)
        integer,allocatable,intent(out)    :: nvir_frag(:)
        integer                            :: i,nmo1,nvir

        nmo1 = size(energies)
        nvir = nmo1 - nocc
        allocate(Lvir(nmo1,nvir,rcq))
        allocate(energies_vir(nvir))
        allocate(nvir_frag(n_fragments))
        Lvir = 0.0D0
        do i = 1, nvir
           Lvir(nocc+i,i,1) = 1.0D0
        end do
        energies_vir = energies(nocc+1:nmo1)
        nvir_frag = 0
        if (n_fragments > 0) nvir_frag(1) = nvir

   end subroutine copy_original_vir

   subroutine generate_fragvirs (fragments,ao_basis_1,spinor_original,ICO,Ladd,nadd_frag)
        use rose_ao
        use rose_mo
        use linear_algebra
        use read_xyzfile
        use rose_import_fragments
        use make_ibo
        use rose_input_data
        use rose_global

        implicit none

        type(molecule_info),intent(in)     :: fragments(:) ! molecule information for the fragments
        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(in)          :: spinor_original      ! Canonical MO-coefficients
        real(8), intent(inout)             :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable,intent(out)   :: Ladd(:,:,:)! Localized MOs in terms of original MOs of B1
        integer, allocatable,intent(out)   :: nadd_frag(:) ! Number of molecular orbitals (minimal complementary basis) for each fragment.

        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis), allocatable        :: spinor_fragment(:) ! MO-coefficients of the fragments
        type(mo_basis)                     :: spinor_ico ! MO-coefficients for the space complementary to the iaos (and the core)
        integer                            :: i
        integer                            :: nmo1,nmo3,nmo4,nvir_add
        real(8), allocatable               :: IVO(:,:,:) ! Transformation to intrinsic virtual orbitals (similar to IAOs, but in compl. space)
        real(8), allocatable               :: energies_ico(:) ! Recanonicalized energies of the complementary space
        integer                            :: nhard(1) ! Number of hard virtuals

        nmo1 = spinor_original%nmo
        nmo3 = size(ICO,2)

        ! The complementary space was obtained by projection. Recanonicalize it here,
        ! close to the additional-virtual selection that needs the resulting energies.
        nhard(1) = nmo3
        allocate(energies_ico(nmo3))
        call recanonize (ICO,spinor_original%energy,energies_ico,nhard)
        call IMO_to_mobasis(ICO,spinor_original,spinor_ico,energies_ico)

        ! We will now construct intrinsic "hard" virtuals on basis of a set of reference virtuals
        ! Reread the fragment spinors to prepare for creating additional intrinsic virtuals
        allocate (ao_basis_2(n_fragments))
        allocate (spinor_fragment(n_fragments))
        call rose_read_fragments (2,fragments,ao_basis_2,spinor_fragment)

        allocate (nadd_frag(n_fragments))
        nadd_frag = spinor_fragment(:)%nmo

        ! For the supermolecule the "core" is already removed as we have only the complementary space in spinor_ico
        nmo4 = sum(nadd_frag) ! the total number of reference virtuals
        ! Determine the number of important (low-lying) virtuals that should always be included
        nvir_add = 0
        do i = 1, nmo3
           if ( energies_ico(i) > additional_virtuals_cutoff) exit
           nvir_add = nvir_add + 1
        end do
        ! For numerical stability we should use fewer low-lying virtuals than we have reference virtuals, check this and adjust if needed.
        nvir_add = min0(nmo4-1,nvir_add) ! the number of important virtuals (plays the role that nocc has in the IAO generation)

        if (nmo4 < 2) then
           write(luout,*) " WARNING : too few reference virtuals to make additional localized virtuals"
           nvir_add = 0
           nmo4 = 0
           nadd_frag = 0
        else if (nvir_add == 0) then
           write(luout,*) "Number of selected virtuals is equal to 0, no additional virtuals added."
           nmo4 = 0
           nadd_frag = 0
        else
           write(luout,*)
           write(luout,*) "####################################"
           write(luout,*) "#      Additional virtuals         #"
           write(luout,*) "####################################"
           write(luout,*)
           write(luout,'(A,I6)') " No. molecular virtuals B3: ",nmo3
           write(luout,'(A,I6)') " No. important virtuals B3: ",nvir_add
           write(luout,'(A,I6)') " No. reference virtuals B4: ",nmo4
           write(luout,*)

           ! Some sanity checks before proceeding
           if (nmo3 .lt. nmo4) then
              write(luout,*) "Too many reference virtuals, please use a smaller AO basis for the fragments"
              error stop "more reference virtuals than supermolecular virtuals, can not proceed"
           end if

           allocate (IVO(nmo3,nmo4,rcq))
           call get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_ico,spinor_fragment,nvir_add,IVO)

           ! Express the localized MOs in the basis of the original CMOs.
           allocate (Ladd(nmo1,nmo4,rcq))
           call matrix_mul (ICO,IVO,Ladd)
           deallocate(IVO)
        end if

        do i = 1, n_fragments
           call delete_adf_objects(ao_basis_2(i))
        end do

   end subroutine generate_fragvirs

   subroutine rose_merge_virtuals (nvir_frag,nadd_frag,Lvir,Ladd,energies_original,energies_vir)
        use make_ibo
        use rose_input_data
        use rose_global

        implicit none

        integer,intent(in)                 :: nadd_frag(:)! Number of additional virtual orbitals assigned to each fragment
        integer, allocatable,intent(inout) :: nvir_frag(:)! Number of localized virtual orbitals assigned to each fragment
        real(8), allocatable,intent(inout) :: Lvir(:,:,:),Ladd(:,:,:)! Localized virtual MOs in terms of original MOs of B1
        real(8),intent(in)                 :: energies_original(:) ! Energies of the core-stripped canonical orbitals
        real(8), allocatable,intent(out)   :: energies_vir(:) ! Recanonicalized virtual energies

        integer                            :: i,j,k,off1,off2,off3
        integer                            :: nmo1,nvir_add,nvir_val
        integer, allocatable               :: nvir_all(:)! Number of localized virtual orbitals assigned to each fragment
        integer                            :: nvir_global(1)
        real(8), allocatable               :: Lvirt(:,:,:)! Merged localized virtual MOs in terms of original MOs of B1
        real(8)                            :: pnorm

        ! Look up the relevant dimensions for later convenience
        if (allocated(Lvir)) then
           nmo1 = size(Lvir,1)
           nvir_val = size(Lvir,2)
        elseif (allocated(Ladd)) then
           nmo1 = size(Ladd,1)
           nvir_val = 0
           if (.not. allocated(nvir_frag)) then
              allocate(nvir_frag(n_fragments))
              nvir_frag = 0
           end if
        else
           error stop "No generated virtual orbitals available for merging"
        end if
        nvir_add = 0
        if (allocated(Ladd)) nvir_add = size(Ladd,2)

        ! Gather and recanonize valence and additional virtuals
        allocate(Lvirt(nmo1,nvir_val+nvir_add,rcq))
        allocate(nvir_all(n_fragments))
        nvir_all = nvir_frag + nadd_frag
 
        ! Loop over fragments to merge the two kinds of virtuals
        if (nvir_add > 0 .or. sum(nvir_frag) > 0) then
           ! Combine Lvir and Ladd for each fragment
           off1 = 0
           off2 = 0
           off3 = 0
           do i = 1, n_fragments
              if (nvir_frag(i) > 0) then
                 Lvirt(:,off1+1:off1+nvir_frag(i),:) = Lvir(:,off2+1:off2+nvir_frag(i),:)
                 off1 = off1 + nvir_frag(i)
                 off2 = off2 + nvir_frag(i)
              end if
              if (nadd_frag(i) > 0) then
                 Lvirt(:,off1+1:off1+nadd_frag(i),:) = Ladd(:,off3+1:off3+nadd_frag(i),:)
                 off1 = off1 + nadd_frag(i)
                 off3 = off3 + nadd_frag(i)
              end if
           enddo
           if (allocated(Lvir)) deallocate(Lvir)
           if (allocated(Ladd)) deallocate(Ladd)
        end if

        ! Recanonization:    
        allocate(energies_vir(nvir_val+nvir_add))
        select case (trim(virtual_mode))
        case ("canonical")
           nvir_global(1) = size(Lvirt,2)
           call recanonize (Lvirt,energies_original,energies_vir,nvir_global)
           nvir_all = 0
           if (size(nvir_all) > 0) nvir_all(1) = nvir_global(1)
        case ("localized")
           call recanonize (Lvirt,energies_original,energies_vir,nvir_all)
        case default
           error stop "Generated virtual orbitals require localized or canonical virtual mode"
        end select
        call move_alloc(Lvirt,Lvir)
        nvir_frag = nvir_all
        
        write(luout,*)
        if (trim(virtual_mode) .eq. "canonical") then
           write(luout,*)"----- Globally recanonicalized virtual energies -----"
        else
           write(luout,*)"----- Recanonicalized virtual energies of fragments -----"
        end if
        write(luout,*)
        write(luout,'(2A12,A20)')"Fragment","MO","Energy"
        off1=0
        do i=1,n_fragments
           do j=1,nvir_frag(i)
             write(luout,'(2I12,F20.9)') i, j, energies_vir(off1+j)
           end do
           off1 = off1 + nvir_frag(i)
        enddo

        if (.false.) then ! Useful while testing, but too much output for production use
        write (luout,'(/A)') " Contribution of original MOs: energies and norms."
        do j = 1, nmo1
           pnorm = 0.0D0
           do k = 1, rcq
              pnorm = pnorm + dot_product(Lvir(j,:,k),Lvir(j,:,k))
           end do
           if (pnorm < 1.0D-8) cycle
           write (luout,'(I4,5X,2F15.6)') j, energies_original(j), pnorm
        end do
        end if

   end subroutine rose_merge_virtuals

   subroutine rose_write_mos (nocc_frag,nvir_frag,Locc,Lvir,energies_occ,energies_vir, &
                              spinor_frozen, spinor_original, IAO)
        use rose_mo
        use rose_utils
        use linear_algebra
        use make_ibo
        use hdf5_interface
        use fchk_interface
        use ADF_interface
        use rose_input_data
        use rose_global

        implicit none

        integer,intent(in)                 :: nocc_frag(:)! Number of localized occupied orbitals assigned to each fragment
        integer,intent(in)                 :: nvir_frag(:)! Number of localized virtual orbitals assigned to each fragment
        real(8), intent(in)                :: Locc(:,:,:) ! Localized occupied MOs in terms of original MOs of B1
        real(8), intent(in)                :: Lvir(:,:,:)! Localized virtual MOs in terms of original MOs of B1
        real(8), intent(in)                :: energies_occ(:) ! Recanonicalized occupied energies
        real(8), intent(in)                :: energies_vir(:) ! Recanonicalized virtual energies
        real(8), intent(in),optional       :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        type(mo_basis),intent(in)          :: spinor_frozen ! Frozen MO-coefficients (left untouched)
        type(mo_basis),intent(inout)       :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis)                     :: spinor_iao ! MO-coefficients for the iao basis
        type(mo_basis)                     :: spinor_ribo ! occupied MO-coefficients after localization and recanonization
        type(mo_basis)                     :: spinor_rido ! valence(+hard) virtual MO-coefficients after localization and recanonization
        integer                            :: ncore ! number of core spinors / orbitals
        integer                            :: nocc  ! number of occupied spinors / orbitals
        integer                            :: nmo1,nvirt
        ! additional options not controlled by the input file:
        logical                            :: extract_RIBO_Fock = .true.
        real(8), allocatable               :: RIBO_Fock(:,:,:),RIBO(:,:,:)

        ! Look up the relevant dimensions for later convenience
        nmo1  = spinor_original%nmo
        ncore = spinor_frozen%nmo
        nocc  = size(Locc,2)
        nvirt = size(Lvir,2)

        ! Express everything in terms of the original AOs.
        call IMO_to_mobasis(Lvir,spinor_original,spinor_rido,energies_vir,nvirt)
        ! Get the Fock matrix in this basis, should be block diagonal.
        if (extract_RIBO_Fock) then
            ! Gather all occupied and virtual (valence + hard)
            allocate(ribo(nmo1,nocc+nvirt,rcq))
            ribo      = 0.D0
            ribo(:,1:nocc,:)            = Locc
            ribo(:,nocc+1:nocc+nvirt,:) = Lvir
            call make_RIBO_Fock(RIBO, spinor_original%energy, RIBO_Fock)
            deallocate(ribo)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         WRITE THE COEFFICIENTS INTO FILES           ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (present(IAO)) then
           call IMO_to_mobasis(IAO,spinor_original,spinor_iao)
           if (trim(file_for_MOs) .eq. "dirac") then
              call write_mos_to_hdf5(spinor_iao,"CHECKPOINT","iao")
           elseif (trim(file_for_MOs) .eq. "h5") then
              call write_mos_to_hdf5(spinor_iao,"CHECKPOINT","iao")
           else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
              call write_mos_to_fchk(spinor_iao,"CHECKPOINT","iao")
           else if (trim(file_for_MOs) .eq. "adf") then
              call write_mos_to_ADF ("iao",spinor_iao)
           else
              print*,"**WARNING** error in the extension of file? I cannot write the resulting IAO to file."
           endif
        end if

        ! Re-insert the original frozen orbitals
        call IMO_to_mobasis(Locc,spinor_original,spinor_ribo,energies_occ,nocc)
        spinor_ribo%shell_id = -1
        call insert_mo_basis (spinor_original,spinor_frozen)
        call insert_mo_basis (spinor_ribo,spinor_frozen)
        ! Add the selected virtuals after the occupied orbitals
        if (nvirt>0) call insert_mo_basis (spinor_ribo,spinor_rido,ncore+nocc)
        if (trim(file_for_MOs) .eq. "dirac") then
           call write_mos_to_hdf5(spinor_ribo)
        elseif (trim(file_for_MOs) .eq. "h5") then
           call write_mos_to_hdf5(spinor_ribo,"CHECKPOINT","ribo")
        else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
           call write_mos_to_fchk(spinor_ribo,"CHECKPOINT","ribo")
        else if (trim(file_for_molfrag) .eq. "ams") then
           call write_ROSE_results (spinor_ribo,nocc_frag,nvir_frag, RIBO_Fock,ncore)
        else if (trim(file_for_MOs) .eq. "adf") then
           call write_mos_to_ADF ("ribo",spinor_ribo) 
           if (extract_RIBO_Fock) call write_data_rose ("ribo",spinor_ribo,nocc_frag,nvir_frag,RIBO_Fock,ncore)
        else
           print*,"**WARNING** error in the extension of file? I cannot write the resulting IAO/IBO to file."
        endif

   end subroutine rose_write_mos
   
   subroutine get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,nocc,IAO,ICO)

        use rose_ao
        use rose_mo
        use rose_input_data, only: version
        use linear_algebra
        use rose_global, only: luout
        use make_iao
        use get_overlap

        implicit none

        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis),intent(in)          :: spinor_original    ! Canonical MO-coefficients
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO-coefficients of the fragments
        integer, intent(in)                :: nocc
        real(8), intent(out)               :: IAO(:,:,:) ! IAO coefficient matrices to transform from B1 to IAO space
        real(8), intent(out), optional     :: ICO(:,:,:) ! Space complementary to the IAO space (also expressed in B1)
  
        integer                            :: nmo1,nmo2
        integer                            :: i, rcq
        real(8), allocatable               :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        real(8), allocatable               :: P12(:,:,:) ! Projection matrix
        real(8), allocatable               :: coeff1_occ(:,:,:) ! MO coefficient matrices.
        real(8), allocatable               :: eig(:),U(:,:,:),VT(:,:,:) ! for Singular Value Decomposition subroutine.

        nmo1 = size(IAO,1)
        nmo2 = size(IAO,2)
        rcq  = size(IAO,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate (S11(nmo1,nmo1,rcq))
        allocate (S22(nmo2,nmo2,rcq))
        allocate (S12(nmo1,nmo2,rcq))
        allocate (S21(nmo2,nmo1,rcq))
        call compute_all_overlaps (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,S11,S12,S21,S22)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Projection Matrix     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! P12 = inv(S11) S12 ---> solve S11 P12 = S12
        ! call solve_lineq(S11,S12,P12) ! solve_lineq(A,B,X) <--> AX=B where A and B are known.
        ! S11 is identity so one can simply set P12 = S12
        ! TODO: check integrity of S11
        allocate (P12(nmo1,nmo2,rcq))
        P12 = S12
        
        !-----Compute Cocc----- (nmo1 x nocc)
        allocate(coeff1_occ(nmo1,nocc,rcq))
        coeff1_occ = 0.D0
        forall (i=1:nocc) coeff1_occ(i,i,1) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Calculate the intrinsic atomic orbital (IAO) basis                                         !
!         Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! -----Compute IAO that span the occupied space----- (nmo1 x nmo2)

        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#   Intrinsic Fragments Orbitals   #"
        write(luout,*) "####################################"
        write(luout,*)

        select case (version)
           case ('Stndrd_2013')
             write(luout,*)" Version: Standard 2013"
             call MakeIAO_Standard_2013(IAO,coeff1_occ,S11,S21,S22,P12)
           case ('Simple_2013')
             write(luout,*)" Version: Simple 2013"
             call MakeIAO_Simple_2013(IAO,coeff1_occ,S11,S12,S21,S22,P12)
           case ('Simple_2014')
             write(luout,*)" Version: Simple 2014"
             call MakeIAO_Simple_2014(IAO,coeff1_occ,S21,S22,P12)
           case default
             write(luout,*)" Default scheme chosen."
             write(luout,*)" Version: Standard 2013"
             call MakeIAO_Standard_2013(IAO,coeff1_occ,S11,S21,S22,P12)
        end select
        write(luout,*)
        deallocate(S22,P12,S21,S12)

        !-----Symmetrically orthonormalize IAOs with S11-----
        call SymOrth(IAO,S11)

        if (present(ICO)) then
           ! Perform SVD to get the complementary space (i.e. the virtual space not contained in the IAO space)
           allocate(eig(min(nmo1,nmo2)))
           allocate(U(nmo1,nmo1,rcq))
           allocate(VT(nmo2,nmo2,rcq))
           U = 0.D0
           VT = 0.D0
           call SVD(IAO,eig,U,VT)
           ICO = U(:,nmo2+1:nmo1,:)
           deallocate(U,VT,eig)
        end if

   end subroutine get_intrinsic_mos

   subroutine express_in_intrinsic_mos (IAO,nocc,Cocc_IAO,Cvir_IAO)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Express occupied and valence virtuals in terms of IAOs     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use linear_algebra
        use rose_global, only: luout

        implicit none

        real(8), intent(in)                :: IAO(:,:,:) ! IAO coefficient matrices to transform from B1 to IAO space
        integer, intent(in)                :: nocc
        real(8), intent(out), allocatable  :: Cocc_IAO(:,:,:) ! orbital coefficient matrices in terms of IAOs.
        real(8), intent(out), optional, allocatable:: Cvir_IAO(:,:,:) ! orbital coefficient matrices in terms of IAOs.
  
        real(8), allocatable               :: IAO_bar(:,:,:), Cvir_IAO_raw(:,:,:)
        integer                            :: nmo1,nmo2,nvir,nvir_val
        integer                            :: i, k, rcq
        real(8), allocatable               :: coeff1_occ(:,:,:),coeff1_vir(:,:,:) ! MO coefficient matrices.
        real(8), allocatable               :: eig(:),U(:,:,:),VT(:,:,:) ! for Singular Value Decomposition subroutine.

        nmo1 = size(IAO,1)
        nmo2 = size(IAO,2)
        rcq  = size(IAO,3)

        !-----Compute Cocc----- (nmo1 x nocc)
        allocate(coeff1_occ(nmo1,nocc,rcq))
        coeff1_occ = 0.D0
        forall (i=1:nocc) coeff1_occ(i,i,1) = 1.D0

        ! Express the occupied MO exactly in terms of IAO
        allocate(Cocc_IAO(nmo2,nocc,rcq))
        allocate(IAO_bar(nmo2,nmo1,rcq))
        call transpose_conjg(IAO,IAO_bar)
        call matrix_mul(IAO_bar,coeff1_occ,Cocc_IAO)
        deallocate(coeff1_occ)

        if (.not.present(Cvir_IAO)) return

        !-----Compute Cvir----- (nmo1 x nvir)
        nvir = nmo1-nocc
        nvir_val = nmo2-nocc
        allocate(coeff1_vir(nmo1,nvir,rcq))
        coeff1_vir = 0.D0
        forall (i=1:nvir) coeff1_vir(i+nocc,i,1) = 1.D0

        ! Express the valence orbitals in terms of IAO
        allocate(Cvir_IAO_raw(nmo2,nvir,rcq))
        call matrix_mul(IAO_bar,coeff1_vir,Cvir_IAO_raw)

        ! Perform SVD to express the valence orbitals in terms of IAO, because Cvir_IAO_raw is not normalized.
        allocate(eig(min(nmo2,nvir)))
        allocate(U(nmo2,nmo2,rcq))
        allocate(VT(nvir,nvir,rcq))
        U = 0.D0
        VT = 0.D0
        call SVD(Cvir_IAO_raw,eig,U,VT)
        allocate(Cvir_IAO(nmo2,nvir_val,rcq))
        Cvir_IAO = 0.D0
        k = 0
        do i=1,min(nmo2,nvir) ! left singular vectors (U) are stored columnwise.
           if (eig(i) .gt. 0.9) then
              k = k + 1
              if (k .le. nvir_val) Cvir_IAO(:,k,:) = U(:,i,:)
           endif
        enddo

        ! check that the number of non-zero singular values are equal to the number of valence virtuals
        if (k .ne. nvir_val) then
           write(luout,*) "ERROR: SVD of the virtual MOs in terms of IFO led to ", &
                         k, " virtuals instead of ", nvir_val
           stop
        end if
        deallocate(U,VT,eig,coeff1_vir,Cvir_IAO_raw)

   end subroutine express_in_intrinsic_mos

   end module rose_intrinsic_mos
