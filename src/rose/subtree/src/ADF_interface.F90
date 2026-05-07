module ADF_interface


contains

!
!    Interface with ADF. We have two situations to deal with:
!    1. only access to KF modules
!       - read only geometry and MO coefficients using direct KF read commands
!       - read AO overlap from precomputed TAPE15
!    2. full access to AMS modules
!       - use FragmentTypeModule to initialize and read fragments
!       - compute AO overlap with the FragOverlapMatrix subroutine from the FDE code
!
!    This module deals with the first situation

    subroutine read_from_ADF (filename,ao_basis,spinors,molecule)

       use read_xyz, only : molecule_info, get_atoms
       use rose_ao
       use rose_mo
       use CoreKF

       implicit none
       character(*), intent(in)            :: filename
       type(mo_basis), intent(inout)       :: spinors
       type(molecule_info), intent(inout)  :: molecule
       type(basis_set_info_t), intent(inout) :: ao_basis

       integer                :: natoms, nao, nmo, nmo_alpha, nmo_beta, nelec, charge, nspin
       integer                :: i
       real(8), allocatable   :: energy(:),coeff(:,:)
       real(8)                :: total_energy
       real(8), allocatable   :: energy_A(:), coeff_A(:,:), energy_B(:), coeff_B(:,:)                          
       real(8), allocatable   :: coord(:,:)
       integer, allocatable   :: atomic_numbers(:)
       logical                :: restricted=.True.

       ! Signal in the basis set type that we have no explicit AO information
       ao_basis%basis_type = -2

       ! Extract, process and store general info from the rkf file
       call extract_other_data(filename, nelec, charge, natoms, nspin, total_energy, atomic_numbers, coord)
       molecule%nelec        = nelec
       molecule%natoms       = natoms
       molecule%charge       = charge
       molecule%nalpha       = (molecule%nelec + (molecule%multiplicity - 1))/2
       molecule%nbeta        = (molecule%nelec - (molecule%multiplicity - 1))/2
       ! process coordinates and atomic_numbers
       call get_atoms (natoms,coord,molecule%firstatom,atomic_numbers)
       molecule%znumber = sum(atomic_numbers)

       ! Extract MO_coeff from the rkf file
       if (restricted) then
        call extract_MOs(filename, 'A', nao, nmo_alpha, energy_A, coeff_A)
        nmo = nmo_alpha
        call alloc_mo_basis(spinors,nao,nmo,mo_type(restricted,.true.))
        spinors%coeff_r(:,:,1) = coeff_A(:,:)
        spinors%energy(:) = energy_A(:)
       else
        !LV: this code needs to be revised for truly unrestricted calculations (with different nmo_alpha and nmo_beta)
        call extract_MOs(filename, 'A', nao, nmo_alpha, energy_A, coeff_A)
        call extract_MOs(filename, 'B', nao, nmo_beta , energy_B, coeff_B)
        nmo = nmo_alpha + nmo_beta
        call alloc_mo_basis(spinors,nao,nmo,mo_type(restricted,.true.))
        do i = 1, nmo_alpha
           spinors%coeff_r(:,2*i-1,1) = coeff_A(:,i)
           spinors%energy(2*i-1) = energy_A(i)
        enddo
        do i = 1, nmo_beta
           spinors%coeff_r(:,2*i,2) = coeff_B(:,i)
           spinors%energy(2*i) = energy_B(i)
        enddo
       endif

       spinors%total_energy = total_energy
       deallocate(energy_A)
       deallocate(coeff_A)
       if (allocated(energy_B)) deallocate(energy_B)
       if (allocated(coeff_B))  deallocate(coeff_B)
              
    end subroutine read_from_ADF

    subroutine write_mos_to_ADF (filename,spinors)

      use read_xyz
      use rose_utils
      use rose_mo
      use CoreKF

      implicit none
      character(len=*), intent(in)  :: filename
      type(mo_basis), intent(inout) :: spinors
      integer                       :: i, j, k, l, nmo_A, lucoef, idx, iu21, ius
      character(len=20)             :: FMT
      logical                       :: tobe
      real(8), allocatable          :: mocoef_A(:), mocoef_B(:), energy_A(:), energy_B(:)
      character(len=128)            :: filepath
      integer                       :: varint(1),nel, chg,nocc   
      real(8),allocatable           :: froc(:)  
      real(8)                       :: varreal(1)                           

      filepath="Molecule_B1_MO1_"//trim(filename)//".rkf"

!     make the copy of the original rkf file if the file does not yet exist
      tobe = kfFileExists(trim(filepath))
      if (.not.tobe) call execute_command_line ('cp Molecule_B1_MO1.rkf '//trim(filepath))

      iu21 = kfOpenFile(trim(filepath))

      if (spinors%type .eq. 1) then
        if (kfQuerySection(iu21, 'Rose Data') == 1) then
          ius = kfDeleteSection(iu21, 'Rose Data')
        endif
        ius = kfCreateSection(iu21, 'Rose Data')
        nmo_A = spinors%nmo
        allocate(mocoef_A(spinors%nao*nmo_A))
        mocoef_A = reshape(spinors%coeff_r,(/nmo_A*spinors%nao/))

        varint(1) = nmo_A
        ius = kfWriteInt(iu21,'Rose Data%nmo_A', varint, 1, 1)
        ius = kfWriteReal(iu21,'Rose Data%Eigen-Bas_A', mocoef_A,nmo_A*spinors%nao, 1)
        ius = kfWriteReal(iu21,'Rose Data%eps_A', spinors%energy(:), nmo_A, 1)
        deallocate(mocoef_A)
        allocate(froc(nmo_A))
        froc = 0.D0
        ius = kfReadReal(iu21, 'General%electrons', varreal, 1, 1)
        nel = INT(varreal(1))
        ius = kfReadReal(iu21, 'Molecule%Charge', varreal, 1, 1)
        chg = INT(varreal(1))  
        nocc = nel - chg
        nocc = nocc/2
        froc(1:nocc) = 2.D0   
        ius = kfWriteReal(iu21,'Rose Data%froc_A', froc(:), nmo_A, 1) 
        deallocate(froc)

      else if (spinors%type .eq. 2) then 

        if (kfQuerySection(iu21, 'Rose Data') == 1) then
          ius = kfDeleteSection(iu21, 'Rose Data')
        endif
        ius = kfCreateSection(iu21, 'Rose Data')

        nmo_A = spinors%nmo/2

        varint(1) = nmo_A

        ius = kfWriteInt(iu21,'Rose Data%nmo_A', varint, 1, 1)

        ius = kfWriteInt(iu21,'Rose Data%nmo_B', varint, 1, 1)

        allocate(mocoef_A(spinors%nao*nmo_A))

        !LV: This only works for perfectly interleaved alpha and beta spinors, should be improved !
        idx = 0
        do i = 1, nmo_A
          mocoef_A(idx+1:idx+spinors%nao) = spinors%coeff_r(:,2*i-1,1)
          idx = idx + spinors%nao
        end do

        ius = kfWriteReal(iu21,'Rose Data%Eigen-Bas_A', reshape(mocoef_A,(/nmo_A*spinors%nao/)), nmo_A*spinors%nao, 1)

        allocate(froc(nmo_A))
        froc = 0.D0
        ius = kfReadReal(iu21, 'General%electrons', varreal, 1, 1)
        nel = INT(varreal(1))
        ius = kfReadReal(iu21, 'Molecule%Charge', varreal, 1, 1)
        chg = INT(varreal(1))  
        nocc = nel - chg
        nocc = nocc/2
        froc(1:nocc) = 1.D0   
        ius = kfWriteReal(iu21,'Rose Data%froc_A', froc(:), nmo_A, 1)

        allocate(energy_A(nmo_A))
        do i = 1, nmo_A
          energy_A(i) = spinors%energy(2*i-1)
        enddo  
        ius = kfWriteReal(iu21,'Rose Data%eps_A', energy_A(:), nmo_A, 1)

        ! Same for the Beta MOs !!!

        allocate(mocoef_B(spinors%nao*nmo_A))

        idx = 0
        do i = 1, nmo_A
          mocoef_B(idx+1:idx+spinors%nao) = spinors%coeff_r(:,2*i,1)
          idx = idx + spinors%nao
        end do

        ius = kfWriteReal(iu21,'Rose Data%Eigen-Bas_B', mocoef_B, nmo_A*spinors%nao, 1)
        ius = kfWriteReal(iu21,'Rose Data%froc_B', froc(:), nmo_A, 1)
        deallocate(froc)

        allocate(energy_B(nmo_A))
        do i = 1, nmo_A
          energy_B(i) = spinors%energy(2*i)
        enddo  
        ius = kfWriteReal(iu21,'Rose Data%eps_B', energy_B(:), nmo_A, 1)

        deallocate(mocoef_A, energy_A, mocoef_B, energy_B)

      endif  

      ius = kfCloseFile(iu21)


    end subroutine write_mos_to_ADF

    subroutine read_overlap_from_ADF (spinors,spinors_frag,S11,S12,S21,S22)

       use CoreKF
       use rose_mo
       use read_xyz
       use linear_algebra

       implicit none
       type(mo_basis), intent(in) :: spinors,spinors_frag(:)

       real(8), intent(inout) :: S11(:,:,:), S22(:,:,:), S12(:,:,:), S21(:,:,:)
       real(8), allocatable   :: smat_ao_triangular(:),smat_ao_full(:,:)
       integer                :: nfragments, naos, naosx, nao1, nao2, nmo1, nmo2
       integer                :: idx, i, j, iao, imo, jao, jmo, naoi, naoj, nmoi, nmoj
       logical                :: bool=.false.
       integer, allocatable   :: nao_frag(:),nmo_frag(:),order_ao(:)
       real(8), allocatable   :: S11_AO(:,:), S22_AO(:,:), S12_AO(:,:)
       real(8), allocatable   :: S_block(:,:,:), transconjg_S_block(:,:,:), S_block_AO(:,:)
       integer                :: iu15, iu21, ius, varint(1)
       character(len=12)      :: name

       ! Read full AO overlap matrix (molecule + all fragments) from file
       bool = kfFileExists("Molecule_B2_S12.t15")
       iu15 = kfOpenFile("Molecule_B2_S12.t15")

       ius = kfReadInt(iu15,'Basis%naos', varint, 1, 1)
       naos = varint(1)
       naosx = (naos*(naos + 1))/2
       allocate(smat_ao_triangular(naosx))
       allocate(smat_ao_full(naos,naos))
       ius = kfReadReal(iu15,'Matrices%Smat', smat_ao_triangular(:), naosx, 1)

       idx = 0
       do i =1,naos
        do j =1, i
          idx = idx + 1
          smat_ao_full(i,j) = smat_ao_triangular(idx)
          smat_ao_full(j,i) = smat_ao_triangular(idx)
        enddo
       enddo

       deallocate(smat_ao_triangular)
       ius = kfCloseFile(iu15)

       ! Compute and check dimensions
       nmo1  = spinors%nmo
       nao1  = spinors%nao
       nfragments = size(spinors_frag(:))
       nao2 = 0
       nmo2 = 0
       do i = 1, nfragments
          nao2 = nao2 + spinors_frag(i)%nao
          nmo2 = nmo2 + spinors_frag(i)%nmo
       end do

       if (nao1+nao2 /= naos) stop "Molecule_B2_S12.t15 has incompatible number of aos"

       ! read reorder information (stored overlap matrix is in different order as the MO coefficients !)
       ! simplification possible if we also have the .rkf of the B2_S12 calculation
       if(allocated(order_ao)) deallocate (order_ao)
       allocate(order_ao(naos))

       iu21 = kfOpenFile("Molecule_B1_MO1.rkf")
       ius = kfReadInt(iu21,'A%npart', order_ao(1:nao1), nao1, 1)
       ius = kfCloseFile(iu21)
       iao = 0
       do i = 1, nfragments
          naoi = spinors_frag(i)%nao
          name = '            '
          write(name,'(a4,i0)') 'frag',i ! This name is used to find the fragment checkpoint file
          iu21 = kfOpenFile(trim(name)//".rkf")
          ius = kfReadInt(iu21,'A%npart', order_ao(nao1+iao+1:nao1+iao+naoi), naoi, 1)
          order_ao(nao1+iao+1:nao1+iao+naoi) = order_ao(nao1+iao+1:nao1+iao+naoi) + iao
          ius = kfCloseFile(iu21)
          iao = iao + naoi
       end do

       ! Distribute over the blocks 1 (basis of the supermolecule) and 2 (concatenated bases of the fragments)
       allocate(S11_AO(nao1,nao1)) 
       allocate(S22_AO(nao2,nao2))
       allocate(S12_AO(nao1,nao2)) 

       do j = 1, nao1
        do i = 1, nao1
          S11_AO(i,j) = smat_ao_full(order_ao(i),order_ao(j))
        enddo
       enddo

       do j = 1, nao2
        do i = 1, nao1
          S12_AO(i,j) = smat_ao_full(order_ao(i),nao1+order_ao(nao1+j))
        enddo
       enddo

       do j = 1, nao2
        do i = 1,  nao2
          S22_AO(i,j) = smat_ao_full(nao1+order_ao(nao1+i),nao1+order_ao(nao1+j))
        enddo
       enddo 

       deallocate (smat_ao_full)
       deallocate (order_ao)
       
      ! Compute the overlap matrices in MO basis:

      ! S11
      call get_overlap_ADF(spinors,spinors,S11_AO,S11)
      deallocate (S11_AO)
       
      ! S22
      ! need to divide S22_AO into blocks and compute S_block in the MO basis.
      imo = 0
      iao = 0
      do i = 1, nfragments
          nmoi = spinors_frag(i)%nmo
          naoi = spinors_frag(i)%nao
          jmo = imo
          jao = iao
          do j = i, nfragments
             nmoj = spinors_frag(j)%nmo
             naoj = spinors_frag(j)%nao
             allocate (S_block(nmoi,nmoj,1))
             allocate (S_block_AO(naoi,naoj))
             S_block_AO(:,:) = S22_AO(iao+1:iao+naoi,jao+1:jao+naoj)
             call get_overlap_ADF(spinors_frag(i),spinors_frag(j),S_block_AO,S_block)
             S22(imo+1:imo+nmoi,jmo+1:jmo+nmoj,1) = S_block(:,:,1)
             if (i /= j) then
                allocate (transconjg_S_block(nmoj,nmoi,1))
                call transpose_conjg(S_block,transconjg_S_block)
                S22(jmo+1:jmo+nmoj,imo+1:imo+nmoi,1) = transconjg_S_block(:,:,1)
                deallocate (transconjg_S_block)
             end if
             jmo = jmo + nmoj
             jao = jao + naoj
             deallocate (S_block)
             deallocate (S_block_AO)
          enddo
          imo = imo + nmoi
          iao = iao + naoi
      enddo
      deallocate (S22_AO)
             
      ! S12
      ! need to divide S12_AO into blocks and compute S_block in the MO basis.
      imo = 0
      iao = 0
      do i = 1, nfragments
          nmoi = spinors_frag(i)%nmo
          naoi = spinors_frag(i)%nao
          allocate (S_block(nmo1,nmoi,1))
          allocate (S_block_AO(nao1,naoi))
          S_block_AO(:,:) = S12_AO(:,iao+1:iao+naoi)
          call get_overlap_ADF(spinors,spinors_frag(i),S_block_AO,S_block)
          S12(:,imo+1:imo+nmoi,1) = S_block(:,:,1)
          deallocate (S_block)
          deallocate (S_block_AO)
          imo = imo + nmoi
          iao = iao + naoi
      enddo
      deallocate (S12_AO)

      call transpose_conjg(S12,S21)

    end subroutine read_overlap_from_ADF

    subroutine get_overlap_ADF(spinor_bra,spinor_ket,SAO,SMO)

       use linear_algebra
       use rose_mo
       use rose_utils
       implicit none
       type(mo_basis), intent(in)  :: spinor_bra
       type(mo_basis), intent(in)  :: spinor_ket
       real(8), intent(in)         :: SAO(:,:)
       real(8), intent(inout)      :: SMO(:,:,:)
       real(8), allocatable        :: spinor_bra_coeff_bar(:,:,:), SAO_tmp(:,:,:), SC(:,:)

       if (spinor_bra%type .eq. 1) then  ! real restricted
          allocate(SAO_tmp(spinor_bra%nao,spinor_ket%nao,1))
          allocate(spinor_bra_coeff_bar(spinor_bra%nmo,spinor_bra%nao,1))
          SAO_tmp = 0.D0
          SAO_tmp(:,:,1) = SAO
          call transpose_conjg(spinor_bra%coeff_r,spinor_bra_coeff_bar)
          call matrix_mul(spinor_bra_coeff_bar,SAO_tmp,spinor_ket%coeff_r,SMO)
          deallocate(spinor_bra_coeff_bar)
          deallocate(SAO_tmp)
       elseif (spinor_bra%type .eq. 2) then ! real unrestricted
          allocate (SC(spinor_bra%nao,spinor_ket%nmo))
          SC  = matmul(SAO,spinor_ket%coeff_r(:,:,1))                               ! multiply with alpha coefficients
          SMO(:,:,1) = matmul(transpose(spinor_bra%coeff_r(:,:,1)),SC)              ! alpha-alpha contribution
          SC  = matmul(SAO,spinor_ket%coeff_r(:,:,2))                               ! multiply with beta coefficients
          SMO(:,:,1) = SMO(:,:,1) + matmul(transpose(spinor_bra%coeff_r(:,:,2)),SC) ! add the beta-beta contribution
          deallocate(SC)
       end if

    end subroutine get_overlap_ADF

    subroutine extract_MOs(filename, spin, naos, nmos, eps_A, MO_coeff)

       !! Routine to extract the MO co-efficients from an adf.rkf file.

       use CoreKF

       implicit none

       logical                :: bool=.false.
       integer                :: iu21, ius, i, j
       
       character, intent(in)    :: spin
       character(*), intent(in) :: filename
       integer, intent(out)   :: naos, nmos
       real(8), allocatable, intent(out)   :: eps_A(:), MO_coeff(:,:)
       real(8), allocatable   :: MO_coeff_tmp(:)
       integer                :: varint(1)

       bool = kfFileExists(trim(filename)//".rkf")
       iu21 = kfOpenFile(trim(filename)//".rkf")

       ius = kfReadInt(iu21,'Basis%naos', varint, 1, 1)
       naos = varint(1)

       ius = kfReadInt(iu21,'A%nmo_'//spin, varint, 1, 1)
       nmos = varint(1)

       if(allocated(MO_coeff)) deallocate (MO_coeff)
       allocate(MO_coeff(naos,nmos))
       if(allocated(MO_coeff_tmp)) deallocate (MO_coeff_tmp)
       allocate(MO_coeff_tmp(naos*nmos))
       ius = kfReadReal(iu21,'A%Eigen-Bas_'//spin, MO_coeff_tmp(:), naos*nmos, 1)

       MO_coeff = reshape(MO_coeff_tmp,(/naos, nmos/))
       deallocate(MO_coeff_tmp)

       if(allocated(eps_A)) deallocate (eps_A)
       allocate(eps_A(nmos))
       ius = kfReadReal(iu21,'A%eps_'//spin, eps_A(:), nmos, 1) 

       ius = kfCloseFile(iu21)

    end subroutine extract_MOs

    subroutine extract_overlap(filename, naos, smatao)

        ! extract overlap matrix from TAPE15

       use CoreKF

       implicit none

       logical                             :: bool
       integer                             :: iu15, ius,  idx, iao, jao
       
       character(*), intent(in)            :: filename
       integer, intent(in)                 :: naos
       real(8), allocatable, intent(out)   :: smatao(:,:)
       real(8), allocatable                :: smatao_tmp(:)
       integer                             :: naosx

       bool = .false.
       bool = kfFileExists(trim(filename))

       iu15 = kfOpenFile(trim(filename))
       if(allocated(smatao)) deallocate(smatao)
       allocate(smatao(naos,naos))
       naosx = (naos*(naos + 1))/2
       if(allocated(smatao_tmp)) deallocate(smatao_tmp)
       allocate(smatao_tmp(naosx))
       ius = kfReadReal(iu15,'Matrices%Smat', smatao_tmp(:), naosx, 1)

       idx = 0 
       do iao =1,naos 
        do jao =1, iao
          idx = idx + 1 
          smatao(iao,jao) = smatao_tmp(idx)
          smatao(jao,iao) = smatao(iao,jao)
        enddo
       enddo

       deallocate(smatao_tmp)

       ius = kfCloseFile(iu15)

    end subroutine extract_overlap

    subroutine extract_other_data(filename, nel, charge, natoms, nspin, tot_enrgy, atomic_numbers, coord)

        !! Routine to retrieve other data from ADF rkf file

       use CoreKF

       implicit none  
       
       real(8), allocatable, intent(inout) :: coord(:,:)
       integer, allocatable, intent(inout) :: atomic_numbers(:)
       character(*), intent(in) :: filename
       logical                :: bool=.false.
       integer                :: iu21, ius,  idx, varint(1), i,j
       real(8),intent(out)    :: tot_enrgy
       integer,intent(out)    :: nel, charge, natoms, nspin
       real(8)                :: varreal(1)

       bool = kfFileExists(trim(filename)//".rkf")
       if (.not.bool) then
           print*," ERROR, missing ",trim(filename)//".rkf"
           stop 'error: missing file'
       end if
       iu21 = kfOpenFile(trim(filename)//".rkf")

       ius = kfReadReal(iu21, 'General%electrons', varreal, 1, 1)
       nel = INT(varreal(1))

       ius = kfReadInt(iu21,'Geometry%nr of atoms', varint, 1, 1)
       natoms = varint(1)

       allocate(atomic_numbers(natoms))
       ius = kfReadInt(iu21,'Molecule%AtomicNumbers',atomic_numbers, natoms, 1)

       allocate(coord(3,natoms))
       ius = kfReadReal(iu21,'Molecule%coords',coord, 3*natoms, 1)
       ! Name of variable has been Capatalized in latest version, check
       if (ius .ne. 0) then
          ius = kfReadReal(iu21,'Molecule%Coords',coord, 3*natoms, 1)
       end if

       ius = kfReadInt(iu21,'General%nspin', varint, 1, 1)
       nspin = varint(1)

       ius = kfReadReal(iu21, 'Molecule%Charge', varreal, 1, 1)
       charge = INT(varreal(1))

       ius = kfReadInt(iu21,'Molecule%nspin', varint, 1, 1)
       nspin = varint(1)

       ius = kfReadReal(iu21,'Energy%Bond Energy', varreal, 1, 1)
       tot_enrgy = varreal(1)
       
    end subroutine extract_other_data

    subroutine write_data_rose(filename,spinors,ribo_occ_per_frag, ribo_vir_per_frag, RIBO_Fock,ncor)

      use read_xyz
      use rose_utils
      use rose_mo
      use CoreKF

      implicit none
      character(len=*), intent(in)  :: filename
      type(mo_basis), intent(in)    :: spinors
      integer, intent(in)           :: ribo_occ_per_frag(:), ribo_vir_per_frag(:)  ! Added for TDDFT calculation   
      real(8), intent(in)           :: RIBO_Fock(:,:,:)                            ! Added for TDDFT calculation   
      real(8),allocatable           :: RIBO_Fock_tmp(:)                            ! Added for TDDFT calculation 
      logical                       :: bool
      character(len=1028)           :: filepath
      integer                       :: i, j, nmo_A, iu21, ius, idx
      integer,intent(in)            :: ncor
      integer                       :: varint(1)

      filepath="Molecule_B1_MO1_"//trim(filename)//".rkf"

      bool = .false.
      bool = kfFileExists(trim(filepath))
 
      iu21 = kfOpenFile(trim(filepath))

        ! This section is required for the TDDFT calculations later with the different orbital spaces (SS, July2021)
      if (spinors%type .eq. 1) then

          nmo_A = spinors%nmo

      else if (spinors%type .eq. 2) then  ! S.S. : Need to re-write this part properly for unrestricted calculation

          nmo_A = spinors%nmo/2  

      endif            

      nmo_A = nmo_A - ncor               !a clumsy correction, to be removed soon


      ius = kfWriteInt(iu21, 'Rose Data%ribo_occ_per_frag_A', ribo_occ_per_frag, size(ribo_occ_per_frag), 1)
      ius = kfWriteInt(iu21, 'Rose Data%ribo_vir_per_frag_A', ribo_vir_per_frag, size(ribo_vir_per_frag), 1)
      varint(1) = size(ribo_occ_per_frag)
      ius = kfWriteInt(iu21,'Rose Data%nfrag', varint, 1, 1)


      allocate(RIBO_Fock_tmp(nmo_A*nmo_A))
      idx = 0
      do i = 1, nmo_A
        do j = 1, nmo_A
          idx = idx + 1
          RIBO_Fock_tmp(idx) = RIBO_Fock(i,j,1) 
        enddo 
      enddo 

      ius = kfWriteReal(iu21,'Rose Data%Fock_RIBO_A', RIBO_Fock_tmp(:), nmo_A*nmo_A, 1)

      deallocate(RIBO_Fock_tmp)
      ius = kfCloseFile(iu21)

    end subroutine write_data_rose    

end module ADF_interface
