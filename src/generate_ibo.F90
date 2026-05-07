   subroutine generate_ibo

        use rose_ao
        use rose_mo
        use rose_utils
        use linear_algebra
        use get_overlap
        use read_xyz
        use import_fragments
        use make_iao
        use make_ibo
        use check
        use properties
        use hdf5_interface
        use fchk_interface
        use ADF_interface
        use rose_input_data
        use rose_global

        implicit none

        complex(8), parameter              :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        real(8), parameter                 :: bohr_to_angs = 0.5291772083D0
        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info), allocatable   :: fragments(:) ! same information for the fragments
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis)                     :: spinor_original_core ! Canonical core MO-coefficients (left untouched)
        type(mo_basis)                     :: spinor_original      ! Canonical MO-coefficients
        type(mo_basis), allocatable        :: spinor_fragment(:) ! MO-coefficients of the fragments
        type(mo_basis)                     :: spinor_iao ! MO-coefficients for the iao basis 
        type(mo_basis)                     :: spinor_ibo ! occupied MO-coefficients after localization 
        type(mo_basis)                     :: spinor_ico ! MO-coefficients for the space complementary to the iaos (and the core)
        type(mo_basis)                     :: spinor_ido ! valence(+hard) virtual MO-coefficients after localization
        type(mo_basis)                     :: spinor_ribo ! occupied MO-coefficients after localization and recanonization
        type(mo_basis)                     :: spinor_rido ! valence(+hard) virtual MO-coefficients after localization and recanonization
        type(mo_basis)                     :: spinor_avas ! MO coefficients for the AVAS approach to define an active space
        type(atom), pointer                :: current_atom ! linked list to store coordinate information
        integer                            :: luxyz,ierr,i,j,k,off1,off2,off3
        integer                            :: nao1,nmo1,nmo2,nmo3, nmo4, nocc,nvir,nvir_val,nvir_add
        integer                            :: imo,jmo,nmoi,nmoj,ncore,nvale,rcq
        integer, allocatable               :: nco_frag(:) ! Number of core orbitals for each fragment.
        integer, allocatable               :: nmo_frag(:) ! Number of molecular orbitals (minimal valence basis) for each fragment.
        integer, allocatable               :: nmo_frag2(:) ! Number of molecular orbitals (minimal complementary basis) for each fragment.
        integer, allocatable               :: nocc_frag(:)! Number of localized occupied orbitals assigned to each fragment
        integer, allocatable               :: nvir_frag(:)! Number of localized valence virtual orbitals assigned to each fragment
        integer, allocatable               :: nadd_frag(:)! Number of localized additional virtual orbitals assigned to each fragment
        integer, allocatable               :: nvall_frag(:)! Number of localized virtual orbitals assigned to each fragment
        real(8), allocatable               :: IAO(:,:,:) ! Transformation from canonical to valence intrinsic fragment orbitals (IAOs)
        real(8), allocatable               :: IBO(:,:,:) ! Transformation to intrinsic bonding orbitals (localized orbitals)
        real(8), allocatable               :: ICO(:,:,:) ! Transformation the space complementary (orthogonal) to the IAOs
        real(8), allocatable               :: IVO(:,:,:) ! Transformation to intrinsic virtual orbitals (similar to IAOs, but in compl. space)
        real(8), allocatable               :: IAO_frag(:,:,:) ! IAO coefficient matrix of one fragment.
        real(8), allocatable               :: Cocc_IAO(:,:,:),Cvir_IAO_val(:,:,:) ! orbital coefficient matrices in terms of IAOs
        real(8), allocatable               :: Cadd_IVO(:,:,:),Cdel_IVO(:,:,:)     ! orbital coefficient matrices in terms of IVOs
        real(8), allocatable               :: Cocc_IBO(:,:,:),Cvir_IBO(:,:,:),Cadd_IBO(:,:,:)  ! Localized MOs in terms of IAOs or IVOs.
        real(8), allocatable               :: Ladd3(:,:,:) ! Localized MOs in terms of the complementary space (dimension NMO3)
        real(8), allocatable               :: Locc(:,:,:),Lvir(:,:,:),Ladd1(:,:,:)! Localized MOs in terms of original MOs of B1
        real(8), allocatable               :: HF_Fock(:,:,:),IAO_Fock(:,:,:) ! Fock matrices
        real(8), allocatable               :: coeff1_occ(:,:,:),coeff1_vir(:,:,:) ! MO coefficient matrices.
        real(8), allocatable               :: energies(:),energies_occ(:),energies_vir(:) ! energies of canonical or localized orbitals.
        real(8), allocatable               :: energies_com(:),energies_add(:) ! complementary orbitals and additional virtuals.
        real(8), allocatable               :: energies_iao(:),energies_ivo(:) ! (pseudo-)energies of IAOs and IVOs
        real(8), allocatable               :: energies_iao_frag(:) ! energies of the IAOs for fragments.
        logical                            :: tobe
        ! additional options not controlled by the input file:
        logical                            :: extract_IAO_Fock = .false.
        logical                            :: extract_IBO_Fock = .false.
        logical                            :: extract_RIBO_Fock = .true.
        real(8), allocatable               :: RIBO_Fock(:,:,:),RIBO(:,:,:)
     
        ! We first define the type of algebra (real/complex/quaternion) that will be used
        call define_algebra_type(spatial_orbitals,restricted,rcq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ INFORMATION FROM FILES (MO COEFF, BASIS SET)         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! For the full molecule:
        if (file_for_MOs .eq. "dirac") then
           call read_from_hdf5 ('MOLECULE',ierr,ao_basis_1,spinor_original,molecule)
           call convert_mo_basis(spinor_original,mo_type(restricted,spatial_orbitals))
        elseif (file_for_MOs .eq. "h5") then
           call read_from_hdf5 ('MOLECULE',ierr,ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "fchk") then
           call read_from_fchk ('MOLECULE',ao_basis_1,spinor_original,restricted,molecule)
        else if (file_for_MOs .eq. "ams") then
           molecule%name = 'adf'
           call read_in_ADF (ao_basis_1,spinor_original,molecule)
        else if (file_for_MOs .eq. "adf") then
           call read_from_ADF ("Molecule_B1_MO1",ao_basis_1,spinor_original,molecule)
        end if

        nmo1 = spinor_original%nmo
        nao1 = spinor_original%nao
        nocc = molecule%znumber - charge
        if (restricted) nocc = nocc/2
        call set_coreval_info (molecule,include_core)

        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "#      Molecular  Coordinates      #"
        write(luout,*) "####################################"
        write(luout,*)
        write(luout,'(A4,A5,3A12)') "Atom","Z","X","Y","Z"

        current_atom => molecule%firstatom
        do while(associated(current_atom))
           write(luout,'(A4,F5.2,3F12.8)') current_atom%name,current_atom%bascharge,current_atom%x
           current_atom => current_atom%next
        end do

        ! Define the list of fragments
        if (n_fragments == 0) then
            n_fragments = molecule%natoms
            if (.not. file_for_MOs .eq. "adf" .and. .not. file_for_MOs .eq. "ams") atom_fragments = .true.
            allocate(frag_ncor(n_fragments))
            allocate(frag_nval(n_fragments))
            allocate(frag_bias(n_fragments))
            frag_ncor = -1 ! signal that this is not set by the user
            frag_nval = -1 ! signal that this is not set by the user
            frag_bias = 0.0 ! by default we do not bias (except the first fragment, but this is done later)
        end if
        if (maxval(frag_bias) < 0.001 ) then
           ! no bias applied, bias fragment 1 with a default value of 0.1 to have a tie-breaker for non-polar bonds
           frag_bias(1) = 0.1
        end if
        allocate (fragments(n_fragments))
        allocate (nco_frag(n_fragments))
        allocate (nmo_frag(n_fragments))
        allocate (ao_basis_2(n_fragments))
        allocate (spinor_fragment(n_fragments))

        if (atom_fragments) then
           ! We make the fragments using information from the supermolecule
           current_atom => molecule%firstatom
           do i = 1, n_fragments
              fragments(i)%natoms  = 1
              fragments(i)%name    = current_atom%name
              fragments(i)%Znumber = current_atom%bascharge
              allocate(fragments(i)%firstatom)
              fragments(i)%firstatom = current_atom
              nullify(fragments(i)%firstatom%next)
              current_atom => current_atom%next
           end do
        else
           do i = 1, n_fragments
              write(fragments(i)%name,'(a4,i0)') 'frag',i ! This name is used to find the fragment checkpoint file
           enddo
        end if

        ! Read the basis sets and relevant spinors of all fragments
        do i = 1, n_fragments
           ! Read the basis set and the full set of spinors for each fragment
           ! TODO: make it possible to mix atomic and molecular fragments
           if (atom_fragments) then
              Call get_atomic_fragment (nint(fragments(i)%Znumber),fragments(i)%firstatom%x,ao_basis_2(i),spinor_fragment(i),&
                                        mo_type(restricted,spatial_orbitals),file_for_molfrag)
           else
              call get_molecular_fragment (fragments(i),ao_basis_2(i),spinor_fragment(i),&
                                           mo_type(restricted,spatial_orbitals),file_for_molfrag)
           end if
           ! Set the predefined information about core and valence
           call set_coreval_info (fragments(i),include_core)
           nmo_frag(i) = fragments(i)%nvale
           nco_frag(i) = fragments(i)%ncore
           if (restricted) nco_frag(i) = nco_frag(i) / 2
           if (restricted) nmo_frag(i) = nmo_frag(i) / 2
           ! Overwrite with user-defined information (if not set, this has a value of -1)
           ! We have no sanity checks, user should make sure ncor and nval make sense for these fragments
           if (allocated(frag_nval)) then
              if (frag_nval(i) .ne. -1) nmo_frag(i) = frag_nval(i)
           end if
           if (allocated(frag_ncor)) then
              if (frag_ncor(i) .ne. -1) nco_frag(i) = frag_ncor(i)
           end if
           ! Keep only the spinors that we need
           call reduce_mo_basis (spinor_fragment(i),nmo_frag(i),nco_frag(i)+1)
        end do

        ! Remove core from the supermolecule. Done here so that we are always consistent with the number of core orbitals from the fragments
        ! First make a copy of the core, so that we can manipulate or just re-insert this later
        call copy_mo_basis (spinor_original,spinor_original_core)
        ! Now remove the core orbitals
        ncore = sum(nco_frag)
        nmo1 = nmo1 - ncore
        nmo2 = sum(nmo_frag)
        if (nmo1 .lt. nmo2) stop "nmo1 shouldn't be inferior to nmo2"
        nmo3 = nmo1 - nmo2
        call reduce_mo_basis (spinor_original_core,ncore)
        call reduce_mo_basis (spinor_original,nmo1,ncore+1)
        nocc = nocc - ncore

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
        write(luout,'(A,I6)') " No. occupied  orbitals B1: ",nocc
        write(luout,'(A,I6)') " No. molecular orbitals B2: ",nmo2
        write(luout,'(A,999I6)') " No. core frag orbitals   : ",nco_frag
        write(luout,'(A,999I6)') " No. MO   frag orbitals   : ",nmo_frag
        write(luout,*)

        allocate (IAO(nmo1,nmo2,rcq))
        allocate (ICO(nmo1,nmo3,rcq))
        call get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,nocc, &
                                IAO,energies_iao,Cocc_IAO,Cvir_IAO_val,ICO)

        ! Write fragment files with IAO's. Only the ones specified in the input.
        ! Separate the IAO coefficient matrix for different fragments defined by list_frag_avas.
        ! For instance, it will be interesting to get the IAOs of the metallic atom in a transition metal complex. 
        ! This can be used in the AVAS theory (https://arxiv.org/abs/1701.07862)
        if (avas_orbitals%array_type /= ' ') then
           ! count orbitals to be put in the active space
           jmo = 0
           do i = 1, n_fragments
              jmo = jmo + avas_orbitals%column_length(i)
           enddo
           allocate(IAO_frag(nmo1,jmo,rcq))
           allocate(energies_iao_frag(jmo))
           jmo = 0
           nmoi = 0
           do j = 1, n_fragments
              do i = 1, avas_orbitals%column_length(j)
                 jmo = jmo + 1
                 call JAGetElement(avas_orbitals,i,j,imo)
                 IAO_frag(:,jmo,:) = IAO(:,nmoi+imo,:)
                 energies_iao_frag(jmo) = energies_iao(nmoi+imo)
              enddo
              nmoi = nmoi + nmo_frag(j)
           end do
           call IMO_to_mobasis(IAO_frag,spinor_original,spinor_iao,energies_iao_frag)
           deallocate(IAO_frag,energies_iao_frag)

           call make_avas(ao_basis_1,ao_basis_1,spinor_original,spinor_iao,nocc,avas_threshold,spinor_avas)
           call insert_mo_basis (spinor_avas,spinor_original_core)
           if (trim(file_for_MOs) .eq. "dirac" .or. trim(file_for_MOs) .eq. "h5") then
              call write_mos_to_hdf5(spinor_avas,"MOLECULE","avas")
           else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
              call write_mos_to_fchk(spinor_avas,ao_basis_1,molecule,"avas")
           else
              stop "**ERROR** Wrong extension, I cannot write the coefficients into file."
           end if
        end if

        ! Calculate partial charges based on the localized occupied MOs expressed in terms of IAOs.
        call molecular_partial_charges(Cocc_IAO,spinor_fragment,fragments,restricted)
        ! Express everything in terms of the original AOs.
        call IMO_to_mobasis(IAO,spinor_original,spinor_iao,energies_iao)

        ! The complementary space was obtained by projection, recanonize to get proper energies.
        ! For this purpose we can use the general canonization routine and define everything as one fragment.
        allocate(nadd_frag(1))
        nadd_frag(1) = nmo3
        allocate(energies_com(nmo3))
        call recanonize (ICO,spinor_original%energy,energies_com,nadd_frag)
        call IMO_to_mobasis(ICO,spinor_original,spinor_ico,energies_com)              

        deallocate(nadd_frag)
        do i = 1, n_fragments
           call delete_adf_objects(ao_basis_2(i))
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Calculate the intrinsic bonding orbital (IBO) basis                                        !
!         Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !-----Compute IBOs (Localized occupied MO in terms of IAOs)----- (nmo2 x nocc)
        allocate(Cocc_IBO(nmo2,nocc,rcq))
        !-----Compute IBOs (Localized occupied MO in terms of original MOs)----- (nmo1 x nocc)
        allocate(Locc(nmo1,nocc,rcq))
        allocate(energies_occ(nocc))
        write(luout,*)
        write(luout,*)
        write(luout,*) "####################################"
        write(luout,*) "# Intrinsic Localized occupied MOs #"
        write(luout,*) "####################################"
        write(luout,*)
        call print_orbital_delocalization(Cocc_IAO,nmo_frag,"before the localization procedure",spinor_original%energy(1:nocc))
        call Localization (Cocc_IAO,Cocc_IBO,nmo_frag,Expnt)
        ! Express the localized MOs in the CMO basis (Locc) and sort them by ascending order:
        call matrix_mul (IAO,Cocc_IBO,Locc)
        call sort_LMO_by_energy(Locc,spinor_original%energy,energies_occ,restricted)
        call print_orbital_delocalization(Locc,IAO,energies_occ,nmo_frag,"after localization")
        ! Express the occupied ILMOs in terms of the original AOs (sorted by ascending-ordered energies).
        call IMO_to_mobasis(Locc,spinor_original,spinor_ibo,energies_occ,nocc)

        ! Recanonization: Sort by fragments and block diagonalize the Fock matrix to get "recanonized" energies.
        allocate(nocc_frag(n_fragments))
        call sort_on_fragments (Cocc_IBO,nmo_frag,nocc_frag,frag_bias)
        ! Express the localized MOs in the CMO basis.
        call matrix_mul (IAO,Cocc_IBO,Locc)
        call recanonize (Locc,spinor_original%energy,energies_occ,nocc_frag)
        call print_orbital_delocalization(Locc,IAO,energies_occ,nmo_frag,"after recanonization")
        ! Express the occupied RILMOs in terms of the original AOs.
        call IMO_to_mobasis(Locc,spinor_original,spinor_ribo,energies_occ,nocc)

        nvir_val = size(Cvir_IAO_val,2)
        allocate(Cvir_IBO(nmo2,nvir_val,rcq))
        allocate(Lvir(nmo1,nvir_val,rcq))
        allocate(energies_vir(nvir_val))
        write(luout,*)
        write(luout,*)
        write(luout,*) "############################################"
        write(luout,*) "# Intrinsic Localized valence virtual MOs  #"
        write(luout,*) "############################################"
        write(luout,*)
        call print_orbital_delocalization(Cvir_IAO_val,nmo_frag,"before the localization procedure",&
                                          spinor_original%energy(nocc+1:nocc+nvir_val))
        call Localization (Cvir_IAO_val,Cvir_IBO,nmo_frag,Expnt)
        ! Express the localized MOs in the CMO basis (Lvir) and sort them by ascending order:
        call matrix_mul (IAO,Cvir_IBO,Lvir)
        call sort_LMO_by_energy(Lvir,spinor_original%energy,energies_vir,restricted)
        call print_orbital_delocalization(Lvir,IAO,energies_vir,nmo_frag,"after localization")

        ! At this stage we are done with the occupied and valence virtual space.
        ! We will now construct intrinsic "hard" virtuals on basis of a set of reference virtuals

        ! Reread the original spinors to prepare for creating these additional intrinsic virtuals
        allocate (nmo_frag2(n_fragments)) ! Need to keep nmo_frag for later.
        do i = 1, n_fragments
           ! Read the basis set and the full set of spinors for each fragment
           if (atom_fragments) then
              call get_atomic_fragment (nint(fragments(i)%Znumber),fragments(i)%firstatom%x,ao_basis_2(i),spinor_fragment(i),&
                                        mo_type(restricted,spatial_orbitals),file_for_molfrag)
           else
              call get_molecular_fragment (fragments(i),ao_basis_2(i),spinor_fragment(i),&
                                           mo_type(restricted,spatial_orbitals),file_for_molfrag)
           end if
           ! Keep only the spinors that we need, now "core" is everything up to the valence space
           nco_frag(i) = nco_frag(i) + nmo_frag(i)
           nmo_frag2(i) = 0
           do j = nco_frag(i)+1, spinor_fragment(i)%nmo
              if (spinor_fragment(i)%energy(j) > fragment_cutoff) exit
              nmo_frag2(i) = j - nco_frag(i)
           end do
           call reduce_mo_basis (spinor_fragment(i),nmo_frag2(i),nco_frag(i)+1)
        end do
     
        ! For the supermolecule the "core" is already removed as we have only the complementary space in spinor_ico
        nmo4 = sum(nmo_frag2) ! the total number of reference virtuals
        ! Determine the number of important (low-lying) virtuals that should always be included
        nvir_add = 0
        do i = 1, nmo3
              if ( energies_com(i) > additional_virtuals_cutoff) exit
              nvir_add = nvir_add + 1
        end do
        ! For numerical stability we should use fewer low-lying virtuals than we have reference virtuals, check this and adjust if needed.
        nvir_add = min0(nmo4-1,nvir_add) ! the number of important virtuals (plays the role that nocc has in the IAO generation)
        allocate(nadd_frag(n_fragments))
        nadd_frag = 0 ! initialize to zero, will become number of additional virtuals for each fragment

        if (nmo4 < 2) then
           write(luout,*) " WARNING : too few reference virtuals to make additional localized virtuals"
           nvir_add = 0
        else if (nvir_add == 0) then
           write(luout,*) "Number of selected virtuals is equal to 0, no additional virtuals added."

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
           call get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_ico,spinor_fragment,nvir_add, &
                                   IVO,energies_ivo)

           ! Express the localized MOs in the basis of the original CMOs.
           allocate (Ladd1(nmo1,nmo4,rcq))
           call matrix_mul (ICO,IVO,Ladd1)

        end if

        ! Array to gather and recanonize valence and additional virtuals
        allocate(IBO(nmo1,nvir_val+nmo4,rcq))
        allocate(energies(nvir_val+nmo4))

        ! Sort valence virtuals by fragments
        allocate(nvir_frag(n_fragments))
        allocate(nvall_frag(n_fragments))
        call sort_on_fragments (Cvir_IBO,nmo_frag,nvir_frag,frag_bias)
        ! Express the localized valence virtual MOs in the CMO basis.
        call matrix_mul (IAO,Cvir_IBO,Lvir)

        ! Add the already sorted additional virtuals
        if (nvir_add > 0 .or. sum(nvir_frag) > 0) then
           ! Combine Lvir and Ladd1 for each fragment
           nvall_frag = nvir_frag + nmo_frag2
           off1 = 0
           off2 = 0
           off3 = 0
           do i = 1, n_fragments
              if (nvir_frag(i) > 0) then
                 IBO(:,off1+1:off1+nvir_frag(i),:) = Lvir(:,off2+1:off2+nvir_frag(i),:)
                 off1 = off1 + nvir_frag(i)
                 off2 = off2 + nvir_frag(i)
              end if
              if (nmo_frag2(i) > 0) then
                 IBO(:,off1+1:off1+nmo_frag2(i),:) = Ladd1(:,off3+1:off3+nmo_frag2(i),:)
                 off1 = off1 + nmo_frag2(i)
                 off3 = off3 + nmo_frag2(i)
              end if
           enddo
           if (nvir_add > 0) deallocate(ICO,IVO,Ladd1)
        end if

        ! Recanonization:
        call recanonize (IBO,spinor_original%energy,energies,nvall_frag)
        deallocate(Cvir_IBO,energies_vir,Lvir)
        call delete_adf_objects(ao_basis_1)
        do i = 1, n_fragments
           call delete_adf_objects(ao_basis_2(i))
        end do

        write(luout,*)
        write(luout,*)"----- Recanonicalized virtual energies of fragments -----"
        write(luout,*)
        write(luout,'(2A12,A20)')"Fragment","MO","Energy"
        off1=0
        do i=1,n_fragments
           do j=1,nvall_frag(i)
             write(luout,'(2I12,F20.9)') i, j, energies(off1+j)
           end do
           off1 = off1 + nvall_frag(i)
        enddo
         
        ! Express everything in terms of the original AOs.
        call IMO_to_mobasis(IBO,spinor_original,spinor_rido,energies,nvir_val+nmo4)
        ! Get the Fock matrix in this basis, should be block diagonal.
        if (extract_RIBO_Fock) then
            ! Gather all occupied and virtual (valence + hard)
            allocate(ribo(nmo1,nmo2+nmo4,rcq))
            ribo      = 0.D0
            ribo(:,1:nocc,:)            = Locc
            ribo(:,nocc+1:nmo2+nmo4,:)  = IBO 
            call make_RIBO_Fock(RIBO, spinor_original%energy, RIBO_Fock)
            deallocate(ribo)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         WRITE THE COEFFICIENTS INTO FILES           ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Re-insert the original core orbitals
        ! Note that one could also add core orbitals of the FRAGMENTS to the IAOs 
        ! (from spinor_fragment(:) before removing the core, and expressed in the basis B1 using change_aobasis subroutine)
        ! However, in contrast to the IAOs which have been orthogonalized, the core orbitals of the fragments are not. 
        ! Hence, one should orthogonalize them with the IAOs if wanted.
        call insert_mo_basis (spinor_original,spinor_original_core)
        call insert_mo_basis (spinor_ribo,spinor_original_core)
        ! Add the selected virtuals after the occupied orbitals
        if (nvir_val+nmo4>0) call insert_mo_basis (spinor_ribo,spinor_rido,ncore+nocc)
        if (trim(file_for_MOs) .eq. "dirac" .or. trim(file_for_MOs) .eq. "h5") then
           call write_mos_to_hdf5(spinor_iao,"MOLECULE","iao")
           call write_mos_to_hdf5(spinor_ribo,"MOLECULE","ribo")
        else if (trim(file_for_MOs) .eq. "fchk" .and. rcq .eq. 1) then
           call write_mos_to_fchk(spinor_iao,ao_basis_1,molecule,"iao")
           call write_mos_to_fchk(spinor_ribo,ao_basis_1,molecule,"ribo")
        else if (trim(file_for_molfrag) .eq. "ams") then
           call write_ROSE_results (spinor_ribo,nocc_frag,nvall_frag, RIBO_Fock,ncore)
        else if (trim(file_for_MOs) .eq. "adf") then
           call write_mos_to_ADF ("iao",spinor_iao)  
           call write_mos_to_ADF ("ribo",spinor_ribo) 
           if (extract_RIBO_Fock) call write_data_rose ("ribo",spinor_ribo,nocc_frag,nvall_frag,RIBO_Fock,ncore)
        else
           print*,"**WARNING** error in the extension of file? I cannot write the resulting IAO/IBO to file."
        endif
        deallocate(energies,IBO)

!!!!!!!!!!!!!!!!!!!!
!      TEST        !
!!!!!!!!!!!!!!!!!!!!

        if (test) then
          call test_IAO_and_IBO(IAO,Cocc_IAO,Locc,spinor_original,ncore,restricted)
        endif
        deallocate(Cocc_IAO,Locc)
        deallocate(IAO)

   ! Fortran Notes:
   !  - routine below needs to be inside the program definition to be able to query sizes, hence the contains
   !  - in the routine below the Sxx blocks are intent(inout) because intent(out) will nullify their allocations

   contains
   
   subroutine get_intrinsic_mos (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,nocc, &
                                 IAO,energies_iao,Cocc_IAO,Cvir_IAO_val,ICO)

        use rose_ao
        use rose_mo
        use linear_algebra
        use make_iao
        use rose_input_data, only: version

        implicit none

        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis),intent(in)          :: spinor_original    ! Canonical MO-coefficients
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO-coefficients of the fragments
        integer, intent(in)                :: nocc


        real(8), intent(inout)             :: IAO(:,:,:) ! IAO coefficient matrices to transform from B1 to IAO space
        real(8), intent(inout), allocatable:: energies_iao(:) ! IAO pseudo-energies (diagonal of Fock matrix in the IAO basis)
        real(8), intent(inout), allocatable, optional:: Cocc_IAO(:,:,:),Cvir_IAO_val(:,:,:) ! orbital coefficient matrices in terms of IAOs.
        real(8), intent(inout), optional   :: ICO(:,:,:) ! Space complementary to the IAO space (also expressed in B1)
  
        real(8), allocatable               :: IAO_bar(:,:,:), Cvir_IAO(:,:,:)
        integer                            :: luxyz,i,k
        integer                            :: nmo1,nmo2,nvir,nvir_val
        integer                            :: rcq
        real(8), allocatable               :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        real(8), allocatable               :: P12(:,:,:) ! Projection matrix
        real(8), allocatable               :: HF_Fock(:,:,:),IAO_Fock(:,:,:) ! Fock matrices
        real(8), allocatable               :: coeff1_occ(:,:,:),coeff1_vir(:,:,:) ! MO coefficient matrices.
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
        !-----Compute Cvir----- (nmo1 x nvir)
        nvir = nmo1-nocc
        nvir_val = nmo2-nocc
        allocate(coeff1_vir(nmo1,nvir,rcq))
        coeff1_vir = 0.D0
        forall (i=1:nvir) coeff1_vir(i+nocc,i,1) = 1.D0

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

        ! Compute the Fock matrix in the IAO basis to get the IAO pseudo-energies.
        allocate(IAO_bar(nmo2,nmo1,rcq))
        call transpose_conjg(IAO,IAO_bar)
        allocate(HF_Fock(nmo1,nmo1,rcq))
        allocate(energies_iao(nmo2))
        HF_Fock = 0.D0
        forall (i=1:nmo1) HF_Fock(i,i,1) = spinor_original%energy(i)
        ! Fock matrix in the IAO basis
        allocate(IAO_Fock(nmo2,nmo2,rcq))
        call matrix_mul(IAO_bar,HF_Fock,IAO,IAO_Fock)
        forall (i=1:nmo2) energies_iao(i) = IAO_Fock(i,i,1)
        deallocate(HF_Fock,IAO_Fock)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Express occupied and valence virtuals in terms of IAOs     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (present(Cocc_IAO)) then
           ! Express the occupied MOs exactly in terms of IAOs
           allocate(Cocc_IAO(nmo2,nocc,rcq))
           call matrix_mul(IAO_bar,S11,coeff1_occ,Cocc_IAO)
        end if
        deallocate(coeff1_occ)

        if (present(Cvir_IAO_val)) then
           ! Express the valence orbitals in terms of IAOs
           allocate(Cvir_IAO(nmo2,nvir,rcq))
           call matrix_mul(IAO_bar,S11,coeff1_vir,Cvir_IAO)
        
           ! Perform SVD to express the valence orbitals in terms of IAOs, because Cvir_IAO is not normalized.
           allocate(eig(min(nmo2,nvir)))
           allocate(U(nmo2,nmo2,rcq))
           allocate(VT(nvir,nvir,rcq))
           U = 0.D0
           VT = 0.D0
           call SVD(Cvir_IAO,eig,U,VT)
           allocate(Cvir_IAO_val(nmo2,nvir_val,rcq))
           k = 0
           do i=1,nvir_val ! left singular vectors (U) are stored columnwise.
              if (eig(i) .gt. 0.9) then
                 k = k + 1
                 Cvir_IAO_val(:,i,:) = U(:,i,:)
              endif
           enddo
    
           ! check that the number of non-zero singular values are equal to the number of valence virtuals
           if (min(nmo2,nvir) .gt. nvir_val) then
              do i=nvir_val+1,min(nmo2,nvir)
                 if (eig(i) .gt. 0.9) k = k+1
              enddo
              if (k .ne. nvir_val) print*,"WARNING: SVD of the virtual MOs in terms of &
                                           IFO led to ", k, " virtuals instead of ", nvir_val
           end if
           deallocate(U,VT,eig,Cvir_IAO)
        end if
        deallocate(coeff1_vir)

        if (present(ICO)) then
           ! Perform SVD to get the complementary space (i.e. unused for the IAOs).
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

   subroutine  compute_all_overlaps (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,S11,S12,S21,S22)

        use rose_ao
        use rose_mo
        use linear_algebra
        use get_overlap

        implicit none
        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis of the supermolecule
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis),intent(in)          :: spinor_original    ! MO-coefficients of the supermolecule
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO-coefficients of the fragments
        real(8),intent(inout)              :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        real(8), allocatable               :: S_block(:,:,:),transconjg_S_block(:,:,:)    ! Subblocks of the overlap matrices
        integer                            :: imo, jmo, i, j, nmo1, nmoi, nmoj, rcq, n_fragments

        if (ao_basis_1%basis_type == -2) then
           ! With only KF access to ADF, overlap matrices are read from file, otherwise we compute them ourselves.
           call read_overlap_from_adf (spinor_original,spinor_fragment,S11,S12,S21,S22)
           return
        end if

        nmo1        = size(S11,1)
        rcq         = size(S11,3)
        n_fragments = size(spinor_fragment)

        !-----Compute S11----- (nmo1 x nmo1) (Should be equal to Identity)
        call get_mo_overlap (ao_basis_1,spinor_original,ao_basis_1,spinor_original,S11)

        !-----Compute S22----- (nmo2 x nmo2)
        imo = 0
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           jmo = imo
           do j = i, n_fragments
              nmoj = spinor_fragment(j)%nmo
              allocate (S_block(nmoi,nmoj,rcq))
              call get_mo_overlap (ao_basis_2(i),spinor_fragment(i),ao_basis_2(j),spinor_fragment(j),S_block)
              S22(imo+1:imo+nmoi,jmo+1:jmo+nmoj,:) = S_block
              if (i /= j) then
                 allocate (transconjg_S_block(nmoj,nmoi,rcq))
                 call transpose_conjg(S_block,transconjg_S_block)
                 S22(jmo+1:jmo+nmoj,imo+1:imo+nmoi,:) = transconjg_S_block
                 deallocate (transconjg_S_block)
              end if
              jmo = jmo + nmoj
              deallocate (S_block)
           enddo
           imo = imo + nmoi
        enddo

        !-----Compute S21----- (nmo2 x nmo1)
        imo = 0
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           allocate (S_block(nmoi,nmo1,rcq))
           call get_mo_overlap (ao_basis_2(i),spinor_fragment(i),ao_basis_1,spinor_original,S_block)
           S21(imo+1:imo+nmoi,:,:) = S_block
           imo = imo + nmoi
           deallocate (S_block)
        enddo

        !-----Compute S12----- (nmo1 x nmo2)
        call transpose_conjg(S21,S12)

   end subroutine compute_all_overlaps

   subroutine change_aobasis(ao_basis_new,ao_basis_old,spinor_old,spinor_new)

        ! Projects spinors expressed in ao_basis_old to ao_basis_new.
        ! Spinors are no longer normalized (unless ao_basis_new exactly spans ao_basis_old).

        use rose_ao
        use rose_mo
        use linear_algebra

        implicit none
        type(basis_set_info_t), intent(in) :: ao_basis_new, ao_basis_old
        type(mo_basis), intent(in)         :: spinor_old
        type(mo_basis), intent(out)        :: spinor_new
        real(8), allocatable               :: Snew(:,:), Snew_inv(:,:) ! the AO overlap matrix of ao_basis_new
        real(8), allocatable               :: Snew_old(:,:) ! the AO overlap matrix between ao_basis_1 and ao_basis_2
        complex(8), allocatable            :: SC(:,:) ! Overlap times MO coefficients
        real(8),    allocatable            :: Snew_tmp(:,:,:), Snew_old_tmp(:,:,:)

        call alloc_mo_basis(spinor_new,ao_basis_new%nao,spinor_old%nmo,spinor_old%type)
        allocate (Snew(ao_basis_new%nao,ao_basis_new%nao))
        allocate (Snew_old(ao_basis_new%nao,ao_basis_old%nao))
        call get_ao_overlap (ao_basis_new,ao_basis_new,Snew)
        call get_ao_overlap (ao_basis_new,ao_basis_old,Snew_old)

        spinor_new%nao = ao_basis_new%nao

        ! C_new = Snew^-1 Snew_old C_old
        allocate(Snew_inv(spinor_new%nao,spinor_new%nao))
        allocate(Snew_tmp(spinor_new%nao,spinor_new%nao,1))
        Snew_tmp = 0.D0
        Snew_tmp(:,:,1) = Snew
        call matrix_inv(Snew_tmp)
        Snew_inv = Snew_tmp(:,:,1)
        deallocate(Snew_tmp)

        allocate(SC(spinor_new%nao,spinor_old%nmo))
        if (spinor_new%type .eq. 1) then
           SC = matmul(Snew_old,spinor_old%coeff_r(:,:,1))
           spinor_new%coeff_r(:,:,1) = matmul(Snew_inv,SC)
        else if (spinor_new%type .eq. 2) then
           SC = matmul(Snew_old,spinor_old%coeff_r(:,:,1))
           spinor_new%coeff_r(:,:,1) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_r(:,:,2))
           spinor_new%coeff_r(:,:,2) = matmul(Snew_inv,SC)
        else if (spinor_new%type .eq. 3) then
           SC = matmul(Snew_old,spinor_old%coeff_c(:,:,1))
           spinor_new%coeff_c(:,:,1) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_c(:,:,2))
           spinor_new%coeff_c(:,:,2) = matmul(Snew_inv,SC)
        else if (spinor_new%type .eq. 4) then
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,1))
           spinor_new%coeff_q(:,:,1) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,2))
           spinor_new%coeff_q(:,:,2) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,3))
           spinor_new%coeff_q(:,:,3) = matmul(Snew_inv,SC)
           SC = matmul(Snew_old,spinor_old%coeff_q(:,:,4))
           spinor_new%coeff_q(:,:,4) = matmul(Snew_inv,SC)
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if
        deallocate(SC)
        deallocate(Snew,Snew_inv,Snew_old)

   end subroutine change_aobasis

   subroutine make_avas(ao_basis_1,ao_basis_2,spinor_B1,spinor_B2,nocc,threshold,spinor_avas)
   ! Ref: E.R. Sayfutyarova, Q. Sun, G.K.-L. Chan, G. Knizia, Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals, J Chem Theory Comput. 13 (2017) 4063–4078. doi:10.1021/acs.jctc.7b00128.
   ! Basis B2 can be the IAOs, and one wants to compute the overlap between this "minimal" basis and a larger basis B1. Especially, one could take only a target set of AOs in B2, like the one corresponding to the d orbitals of a transition metal. Then, a complete active space is automatically constructed according to the calculated overlap matrices for occupied and virtual orbitals.

        use rose_ao
        use rose_mo
        use rose_utils
        use read_xyz
        use linear_algebra
        use get_overlap
        use rose_input_data

        implicit none
        type(basis_set_info_t), intent(in) :: ao_basis_1 ! the ao basis B1
        type(basis_set_info_t), intent(in) :: ao_basis_2 ! the ao basis B2
        type(mo_basis), intent(out)        :: spinor_avas ! the AVAS spinors: inactive occupied, active occupied, active virtual, inactive virtual, all expressed in B1
        type(mo_basis),intent(in)          :: spinor_B1  ! MO-coefficients of basis B1
        type(mo_basis),intent(in)          :: spinor_B2  ! MO-coefficients of basis B2 (only valence)
        integer, intent(in)                :: nocc
        real(8), intent(in)                :: threshold

        integer                 :: i,j,k,ierr
        integer                 :: nao1,nmo1,nao2,nmo,nmo2,rcq,nact_occ,nact_vir,nvir
        real(8), allocatable    :: Proj(:,:,:),Socc(:,:,:),Cocc(:,:,:),Svir(:,:,:),Cvir(:,:,:)
        real(8), allocatable    :: Cocc_tilde(:,:,:),Cvir_tilde(:,:,:),Ctilde(:,:,:)
        real(8), allocatable    :: eigvec_occ(:,:,:),eigvec_vir(:,:,:)
        real(8), allocatable    :: eigval_occ(:),eigval_vir(:),eigval(:)
        real(8), allocatable    :: S21(:,:,:),S12(:,:,:),S22(:,:,:)
        character(len=13)       :: fmt_tmp
       
        write(*,*)
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(*,*) "@                                                @"
        write(*,*) "@        Automated Valence Active Space          @"
        write(*,*) "@                    (AVAS)                      @"
        write(*,*) "@                                                @"
        write(*,*) "@               By Bruno Senjean                 @"
        write(*,*) "@                                                @"
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(*,*)

        write(*,'(A,E10.2)') " Threshold: ",threshold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CHECK IF OPTIONS ARE CORRECT         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call define_algebra_type(spatial_orbitals,restricted,rcq)

        nao2 = spinor_B2%nao
        nmo2 = spinor_B2%nmo
        if (mod(nmo2,2) .ne. 0 .and. .not. restricted) stop "**ERROR** Number of spinors to keep should be an even number"

        nao1 = spinor_B1%nao
        nmo1 = spinor_B1%nmo

        if (nmo1 .lt. nmo2) stop "nmo1 shouldn't be inferior to nmo2"

        write(*,*)
        write(*,*) "####################################"
        write(*,*) "#         Basis set size           #"
        write(*,*) "####################################"
        write(*,*)
        write(*,'(A,I6)') " No. basis functions    B1: ",nao1
        write(*,'(A,I6)') " No. molecular orbitals B1: ",nmo1
        write(*,'(A,I6)') " No. occupied  orbitals B1: ",nocc
        write(*,'(A,I6)') " No. basis functions    B2: ",nao1
        write(*,'(A,I6)') " No. molecular orbitals B2: ",nmo2
        write(*,*)      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(S21(nmo2,nmo1,rcq))
        allocate(S12(nmo1,nmo2,rcq))
        allocate(S22(nmo2,nmo2,rcq))
        call get_mo_overlap (ao_basis_2,spinor_B2,ao_basis_2,spinor_B2,S22)
        call get_mo_overlap (ao_basis_2,spinor_B2,ao_basis_1,spinor_B1,S21)
        call transpose_conjg(S21,S12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Projection Matrices    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! B2 form an orthonormalized MO basis, so sigma^(-1) in Eq.13 of the article is the identity, so not considered here.
        allocate(Proj(nmo1,nmo1,rcq))
        call matrix_mul(S12,S21,Proj)

        ! Compute the projected overlap matrices Socc = C†occ P Cocc and Svir = C†vir P Cvir
        ! In the MO basis, Cocc and Cvir are diagonal matrices.
        nvir = nmo1 - nocc
        allocate(Socc(nocc,nocc,rcq))
        allocate(Svir(nvir,nvir,rcq))
        Socc = Proj(1:nocc,1:nocc,:)
        Svir = Proj(nocc+1:nocc+nvir,nocc+1:nocc+nvir,:)

        ! deallocate matrices not used anymore:
        deallocate(S12,S21,S22,Proj)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DIAGONALIZE PROJECTED OVERLAP MATRICES TO DEFINE THE AVAS     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(eigvec_occ(nocc,nocc,rcq))
        allocate(eigval_occ(nocc))
        allocate(eigvec_vir(nvir,nvir,rcq))
        allocate(eigval_vir(nvir))
        eigvec_occ = Socc
        eigvec_vir = Svir
        call diag_matrix(eigvec_occ,eigval_occ)
        call diag_matrix(eigvec_vir,eigval_vir)

        ! Count the number of (spatial) active orbitals to be considered.
        ! There will be at most 2*nmo2 of them.
        nact_occ = 0
        nact_vir = 0
        do i=1,nocc
         if (eigval_occ(i) > threshold) nact_occ = nact_occ + 1
        enddo
        do i=1,nvir
         if (eigval_vir(i) > threshold) nact_vir = nact_vir + 1
        enddo
        if ((nact_occ+nact_vir) .gt. 2*nmo2) print*,"**WARNING**: nact is larger than 2*nmo2 which shouldn't be possible."

        ! Now we sort the coefficients to have frozen-occupied / actives / unoccupied-virtuals:
        ! Eigenvectors are sorted in an ascending order (lowest eigenvalues first) such that occupied orbitals
        ! are already sorted, but the virtuals need to be in reverse order.
        allocate(Ctilde(nmo1,nmo1,rcq))
        allocate(eigval(nmo1))
        Ctilde = 0.0
        Ctilde(1:nocc,1:nocc,:) = eigvec_occ
        eigval(1:nocc)          = eigval_occ
        do j=1,nvir
          Ctilde(nocc+1:nocc+nvir,nocc+j,:) = eigvec_vir(:,nvir-j+1,:)
          eigval(nocc+j)                    = eigval_vir(nvir-j+1)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       WRITE INFORMATION ABOUT AVAS IN THE OUTPUT     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(*,*)
        write(*,*) "####################################"
        write(*,*) "#    Active space construction     #"
        write(*,*) "####################################"
        write(*,*)
        write(*,'(A,I6)') " No. of frozen occupied: ",nocc-nact_occ
        write(*,'(A,I6)') " No. of active occupied: ",nact_occ
        write(*,'(A,I6)') " No. of active virtuals: ",nact_vir
        write(*,'(A,I6)') " No. of frozen virtuals: ",nvir-nact_vir
        write(*,*)
        write(*,*) "--------- Eigenvalues of the Projected overlap matrices (active MOs) ---------"
        write(*,*)
        write(*,'(A4,A1,A24)') "MO"," ","eigval"
        do i=nocc-nact_occ+1,nocc+nact_vir
          write(*,'(I4,A1,F24.16)') i," ",eigval(i)
        enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       WRITE COEFFICIENT MATRICES THAT WILL BE USED IN POST-HF CALCULATION     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call IMO_to_mobasis(Ctilde,spinor_B1,spinor_avas)

        write(*,*) 
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(*,*) "@                                                @"
        write(*,*) "@        Automated Valence Active Space          @"
        write(*,*) "@                    (AVAS)                      @"
        write(*,*) "@                                                @"
        write(*,*) "@                     END                        @"
        write(*,*) "@                                                @"
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

   end subroutine make_avas

   end subroutine generate_ibo
