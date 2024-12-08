module visual_dirac

  use visual_datatypes
  use visual_cfg
  use visual_utils
  use visual_mpi
  use dirac_dmat_factory

  !use dirac_rsp_utils

! DIRAC modules
  use checkpoint
  use basis_set_datatypes
  use dirac_rsp_utils
  use dirac_conmat_lao_utils

#ifdef VAR_MPI
   use interface_to_mpi
   use integer_kind_mpilib
#endif

  implicit none

  public visual_input_show_setup

  public visual_initialize_and_sync_nodes
  public visual_initialize_data_from_dirac
  public visual_update_data

  public dirac_dmat_init
  public get_xyzcomp_from_rspvec


#ifdef VAR_MPI
  public synchronize_grid_data
  public synchronize_grid_function_data
  public synchronize_data_from_dirac
  public synchronize_dmat
#endif


  private

contains

! -----------------------------------------------------------------------------
  subroutine dirac_dmat_init(mode, cmo, &
                             dmat_ao, dmat_mo, &
                             grid_function)
    character(*), intent(in) :: mode
    type(cmo_t), intent(in) :: cmo
    type(dmat_t), optional, intent(out) :: dmat_ao
    type(dmat_t), optional, intent(out) :: dmat_mo
    type(grid_function_t), optional, intent(inout):: grid_function

    type(cmo_t) :: cmo_i
    type(response_info_t) :: rsp_info
    type(conmat_info_t) :: conmat_info
    type(dmat_t) :: tmpmat_ao, tmpmat_mo
    integer :: nao, nmo, nz
    integer :: i, ivec, ivec_rsp, ivec_con, rspvec_ndim, rspvec_nz, n1

    character(len=:),allocatable :: f_name
    integer :: iq, im, iz, ipqtoq(4), ihqmat(4)
    real(8) :: factor(4)
!
!   this subroutine constructs density matrices for one specific grid function
!
!   it is based on how density matrices were set in old VISUAL code;
!   I group the old code here, so that it is in one place:
!   * it should now be easier to modify it, e.g., to read density matrices from outside,
!   * but still, the code is too coupled to other structures (cmo, response, etc.)
!

    nao = sum(cmo%n_ao)
    nmo = sum(cmo%n_mo)
    nz  = cmo%nz

    select case (trim(mode))
      case ('unperturbed')

!       construct unperturbed density matrix from occupation numbers
    
        tmpmat_mo%mat_dim = nmo
        tmpmat_mo%mat_nz  = nz
        call tmpmat_mo%alloc()
        do i = 1, nmo
          tmpmat_mo%elements((i-1)*nmo + i) = cmo%occupation(i)
        end do
    
        if (present(dmat_mo)) dmat_mo = tmpmat_mo
    
        if (present(dmat_ao)) then
          dmat_ao%mat_dim = nao
          dmat_ao%mat_nz  = nz
!         fixme symmetries hardcoded! assuming C1
          dmat_ao%irep = 0
          dmat_ao%ih_sym = 1
          dmat_ao%tr_sym = 1
          call dmat_ao%alloc()
          call dmat_transform_mo2ao(cmo, tmpmat_mo, dmat_ao)
          call dmat_ao%scale_elements(2.0d0)
        end if
    
        call tmpmat_mo%erase()

      case ('xlr')

!       construct perturbed density matrix from rsp vectors on PAMXVC
!       D = \sum_{\mu\nu} C^\dagger_{\mu p} \kappa_{pq} C_q\nu, 
!       where \kappa are elements of the LR solution vector

!       read response data from 'RESPONSE.h5' file
        if (.not.present(grid_function)) stop
        f_name = grid_function%rspvec_file
        call rsp_info%read_from_checkpoint(trim(f_name))
        ivec_rsp = grid_function%rspvec_id

!       scatter response vector to matrix:
        tmpmat_mo%mat_dim = nmo
        tmpmat_mo%mat_nz = nz
        call tmpmat_mo%alloc()
        call rsp_info%vectors(ivec_rsp)%scatter_to_matrix_q(tmpmat_mo%elements, nmo, nz)

!       dmat inherits symmetry of the response vector:
        !call tmpmat_mo%assign_ih(rsp_info%vectors(ivec_rsp)%get_ih_sym())
        !call tmpmat_mo%assign_tr(rsp_info%vectors(ivec_rsp)%get_tr_sym())
        !call tmpmat_mo%assign_irep(rsp_info%vectors(ivec_rsp)%get_irep())
 
        if (present(dmat_mo)) dmat_mo = tmpmat_mo

        if (present(dmat_ao)) then
          dmat_ao%mat_dim = nao
          dmat_ao%mat_nz  = nz
          call dmat_ao%alloc()
!         construct non-symmetrized AO dmat:
          call construct_perturbed_dmat(dmat_ao%elements,&
                                        tmpmat_mo%elements,&
                                        ivec_rsp, rsp_info, cmo%coeff)
          !call symmetrize_mat(dmat_ao%elements)
          do iz = 1, nz
            iq = rsp_info%vectors(ivec_rsp)%get_ipqtoq(iz)
            im = rsp_info%vectors(ivec_rsp)%get_ihqmat(iq)
            if (im .eq. 1) then
              factor(iz) = 1.0d0
            else if (im .eq. 2) then
              factor(iz) =-1.0d0
            end if
          end do
          call pertdenmat_sym(nao, nz, dmat_ao%elements, factor)
          call dmat_ao%scale_elements(2.0d0)
        end if

        call tmpmat_mo%erase()
        call rsp_info%erase()

      case ('lao-connection-1')

!       gosia fixme:
!       now, i'm reusing old code here

!       read response data from 'RESPONSE.h5' file
        if (.not.present(grid_function)) stop

        call get_xyzcomp_from_rspvec(grid_function, ivec_con)

        if (present(dmat_ao)) then
          dmat_ao%mat_dim = nao
          dmat_ao%mat_nz  = nz
          call dmat_ao%alloc()

!
!         fixme list:
!         * reading again rspvec_file, only to get the label and set ivec_con;
!           this is to ensure that the same vector component is read from "TBMO" file,
!           but it should be changed once TBMO is refactored (as now is PAMXVC)
!

          call conmat_info%read_from_checkpoint('TBMO.h5')
          do ivec = 1, conmat_info%nr_comp

            if (conmat_info%vectors(ivec)%id == ivec_con) then
              tmpmat_mo%mat_dim = nmo
              tmpmat_mo%mat_nz = nz
              call tmpmat_mo%alloc()
              call conmat_info%vectors(ivec)%get_copy_of_values(tmpmat_mo%elements,nmo*nmo*nz)
              if (rs_wf%hamiltonian/='LEVYLE') then
                call tmpmat_mo%scale_elements(0.5d0)
              end if
              
              call tmpmat_mo%get_dmat_q_mo_sector('pi ', cmo%ioffs, cmo%n_mo, cmo%n_fsym)
              
              call construct_perturbed_dmat_lao_e(dmat_ao%elements,&
                            tmpmat_mo%elements,&
                            cmo%coeff)
              call construct_perturbed_dmat_lao_p(dmat_ao%elements,&
                            tmpmat_mo%elements,&
                            cmo%coeff)
              
              factor(1) = -1.0d0
              factor(2:4) = 1.0d0
              call pertdenmat_sym(nao, nz, dmat_ao%elements, factor)
  

            end if
          end do

          call dmat_ao%scale_elements(2.0d0)
          call tmpmat_mo%erase()
          call conmat_info%erase()
        end if

      case default
        stop
    end select
  end subroutine dirac_dmat_init

! ------------------------------------------------------------------------------
  subroutine visual_initialize_and_sync_nodes(hello)

!   hello = .true.:
!   * make sure that all workers have all information to start
!   * wake up nodes, pass initial information,
!     and leave the rest to master 
!   * happens when called from dirnod; on workers
!
!   hello = .false.:
!   * proceed to calculations
!   * happens when called from dirac; on master
!
    use visual_mpi
#ifdef VAR_MPI
    use interface_to_mpi
    use integer_kind_mpilib
#endif
#include "dirac_partask.h"  
! legacy
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h" ! cvalue
#include "orgcom.h" ! gagorg

    logical, intent(in), optional :: hello

#ifdef VAR_MPI
    ! wake up workers
    if (is_parallel()) then
      call dirac_parctl( VISUAL_PAR )
      if (get_my_rank() /= irank_master) then
        call sync_visual_cfg()
        return
      end if
    end if
#endif

    ! then, the rest happens on master only
    if (.not. rs_input%set_cvalue) then
      rs_input%cvalue = cval  ! speed of light used in Dirac
    end if

    rs_input%gauge_origin(1:3) = gagorg(1:3) ! magnetic gauge origin defined globally in Dirac

#ifdef VAR_MPI
    if (is_parallel()) call sync_visual_cfg()
#endif

  end subroutine visual_initialize_and_sync_nodes


! ------------------------------------------------------------------------------
  subroutine visual_update_data(cmo, print_level)
    type(cmo_t), intent(inout):: cmo
    integer, intent(in), optional :: print_level
    logical:: small_component
    integer:: i

!   update MO data that was not read from CHECKPOINT:
    call setup_orbital_string(cmo)
    call setup_mo_offsets(cmo)

!   update grid function data with information specific to this job
    do i = 1, rs_input%nr_grid_functions
      call modify_grid_function_due_to_method(rs_input%grid_functions(i))
    end do

  end subroutine


! ------------------------------------------------------------------------------
  subroutine modify_grid_function_due_to_method(grid_function)
    class(grid_function_t), intent(inout):: grid_function

    grid_function%gf_type = "unperturbed"
    if (grid_function%wf_order > 0) then
      grid_function%gf_type = "response"
    end if

    if (grid_function%lao > 0) then
      if (grid_function%max_ao_m_order < 1) grid_function%max_ao_m_order = 1
      if (grid_function%lao_skip(2) == 0) then
!       do not skip LAO-direct contributions
        grid_function%gf_type = "composite"
      end if
    end if

    if (rs_wf%hamiltonian=='LEVYLE') then
      if ((grid_function%name == 'j') .or. &
          (grid_function%name == 'ndip') .or. &
          (grid_function%name == 'bdip')) then
        if (grid_function%max_ao_g_order < 1) then
          grid_function%max_ao_g_order = 1
        end if
      end if
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine visual_initialize_data_from_dirac(mol, mo, ao_basis)
    type(mol_t), intent(out):: mol
    type(cmo_t), intent(out):: mo
    type(basis_set_info_t), intent(out) :: ao_basis
    character*(:),allocatable:: group, section
    character(len=80) :: wf_data
    logical:: dirac_mol, dirac_mo, dirac_ao, tobe
!
!   read data from DIRAC's CHECKPOINT file(s) into DTs
!
    group = '/input/molecule/'
    call checkpoint_query(group, dirac_mol)
    if (dirac_mol) then
      call read_mol_from_dirac(mol, group)
      if ( rs_input%print_level > 2 ) then
        write(*, *) 'Molecule data read from DIRAC CHECKPOINT'
        call mol%display
      end if
    end if

    group = '/result/wavefunctions/'
    call checkpoint_query(group, dirac_mo)
    if (dirac_mo) then
!     todo: generalize to other wfs
      section = 'scf'
      call read_mo_from_dirac(section, mo)
      if ( rs_input%print_level > 2 ) then
        write(*, *) 'MO data read from DIRAC CHECKPOINT'
        call mo%display
      end if

!     read other data stored in this section:
      call checkpoint_query(group//trim(section)//'/hamiltonian', tobe)
      if (tobe) then
        wf_data=''
        call checkpoint_read(group//trim(section)//'/hamiltonian',sdata=wf_data)
        if (len_trim(wf_data) > 0) rs_wf%hamiltonian = trim(wf_data)
      else
        print*, "WARNING: no hamiltonian on CHECKPOINT. Assuming Dirac-Coulomb"
      end if

      call checkpoint_query(group//trim(section)//'/xc_functional', tobe)
      if (tobe) then
        wf_data=''
        call checkpoint_read(group//trim(section)//'/xc_functional',sdata=wf_data)
        if (len_trim(wf_data) > 0) rs_wf%xc_functional = trim(wf_data)
      else
        print*, "WARNING: no xc_functional on CHECKPOINT. Assuming non-DFT job"
      end if
      if ( rs_input%print_level > 2 ) then
        write(*, *) 'WF data read from DIRAC CHECKPOINT'
        call rs_wf%display
      end if
    end if

!
!   fixme list:
!   * ao_basis is kept as a placeholder (unused now)
!
    !group = '/input/aobasis/'
    !call checkpoint_query(group, dirac_ao)
    !if (dirac_ao) then
    !  call read_ao_from_dirac(ao_basis)
    !  if ( rs_input%print_level > 2 ) then
    !    write(*, *) 'AO basis data read from DIRAC CHECKPOINT'
    !    call ao_basis%display
    !  end if
    !end if

    !if (dirac_mo .and. dirac_ao) then
    !  call check_data_consistency(ao_basis, mo)
    !end if

!
!   read data from DIRAC's common blocks
!

  end subroutine visual_initialize_data_from_dirac


! ------------------------------------------------------------------------------
  subroutine read_mol_from_dirac(mol, group)
#include "codata.h" !xtang
    type(mol_t), intent(out) :: mol
    character(len=*),intent(in) :: group
    real(8), allocatable :: temp(:)

    call checkpoint_read(group//'n_atoms', idata=mol%n_atoms)

    allocate(mol%symbols(mol%n_atoms))
    allocate(mol%geometry(3, mol%n_atoms))
    allocate(mol%nuc_charge(mol%n_atoms))
    allocate(temp(3*mol%n_atoms))

    call checkpoint_read(group//'/symbols',   sarray=mol%symbols, slen=4)
    call checkpoint_read(group//'geometry',   rdata=temp)
    call checkpoint_read(group//'nuc_charge', rdata=mol%nuc_charge)

!   Note: on CHECKPOINT.h5, atomic coordinates are stored in Angstroms;
!   here, they are rescaled to AU, because in VISUAL module we work in AU
    temp = temp/xtang
    mol%geometry = reshape(temp, (/3, mol%n_atoms/))

    deallocate(temp)

  end subroutine read_mol_from_dirac


! ------------------------------------------------------------------------------
  subroutine read_mo_from_dirac(section,cmo)
    type(cmo_t), intent(inout) :: cmo
    character(len=*),intent(in) :: section
    character*(:),allocatable:: group
    logical:: tobe_dirac
    integer:: nao,npo,nmo,nto,ncoef,nfsym,nz
    integer:: n_basis(2),n_po(2),n_mo(2)
    real(8):: total_energy
    integer, allocatable   :: ibeig(:), index(:)
    real(8), allocatable   :: eig(:), mo_coeff(:), mo_coeff_q(:,:,:)
    integer:: maxopr,inversion
    logical:: sym_info
    integer:: i
!
!   this subroutine is based on how things are done in dircmo.F90 and exacorr
!

    call checkpoint_query ('/result/symmetry',sym_info)
    if (.not.sym_info) then
       print*, "WARNING: old CHECKPOINT, assuming NOSYM"
       maxopr = 0
       inversion = 1
    else
       call checkpoint_read ('/result/symmetry/maxopr',idata=maxopr)
       call checkpoint_read ('/result/symmetry/inversion',&
            idata=inversion)
    end if

    group = '/result/wavefunctions/'//trim(section)//'/'
    call checkpoint_query(group//'mobasis', tobe_dirac)
    if (.not.tobe_dirac) call visual_quit('No MO data available on CHECKPOINT. Stop.')

    n_basis = 0
    n_mo = 0
    n_po = 0
    total_energy = 0.0d0
    call checkpoint_read(group//'mobasis/nz',idata=nz)
    call checkpoint_read(group//'mobasis/n_basis',idata=n_basis)
    call checkpoint_read(group//'mobasis/n_mo',idata=n_mo)
    call checkpoint_read(group//'mobasis/n_po',idata=n_po)
    call checkpoint_read(group//'energy',rdata=total_energy)

    ! check symmetry and get dimensions of the coefficient array
    if (inversion == 2) then
       nao = sum(n_basis)
       npo = sum(n_po)
       nmo = sum(n_mo)
    else
       nao = n_basis(1)
       npo = n_po(1)
       nmo = n_mo(1)
    end if
    nto = npo+nmo
    ncoef = nz*nao*nmo

    cmo%n_ao = n_basis
    cmo%n_mo = n_mo
    cmo%n_po = n_po
    cmo%nz   = nz
    cmo%total_energy = total_energy

    call cmo%alloc_all()

    call checkpoint_query(group//'mobasis/orbitals',tobe_dirac)
    if (tobe_dirac) then
      allocate(mo_coeff(ncoef))
      mo_coeff = 0.0d0
      call checkpoint_read(group//'/mobasis/orbitals',rdata=mo_coeff)
      cmo%coeff  = mo_coeff(1:ncoef)
      deallocate(mo_coeff)
    end if

    call checkpoint_query(group//'mobasis/eigenvalues',tobe_dirac)
    if (tobe_dirac) then
      allocate(eig(nmo))
      eig = 0.0d0
      call checkpoint_read(group//'/mobasis/eigenvalues',rdata=eig)
      cmo%energy   = eig(npo+1:nmo)
      deallocate(eig)
    end if

    call checkpoint_query(group//'mobasis/symmetry',tobe_dirac)
    if (tobe_dirac) then
      allocate(ibeig(nmo))
      ibeig = 0
      call checkpoint_read(group//'/mobasis/symmetry',idata=ibeig)
      cmo%symmetry = ibeig(npo+1:nmo)
      deallocate(ibeig)
    end if

    allocate(index(nmo))
    do i = 1, nmo
       index(i) = i ! we set the index compatible with the fact that we store all mo's
    end do
    cmo%index    = index(npo+1:nmo)
    deallocate(index)

    cmo%data_is_read = .true.

  end subroutine

! ------------------------------------------------------------------------------
  subroutine read_ao_from_dirac(ao_basis)
    type(basis_set_info_t), intent(inout) :: ao_basis
    logical ::  small_component
    character*(:),allocatable:: group
    integer, allocatable:: orbmom(:),n_prim(:),n_cont(:)
    real(8), allocatable:: priexp(:),priccf(:),coord(:)
    integer:: nao, nshells, npriexp, iPrimStart, iPrimEnd
    integer:: iao, ibas, ishell, ierr
    integer:: maxopr,inversion
    logical:: sym_info
!
!   this subroutine is based on how things are done in dirgp.F and exacorr
!   gosia fixme unused now (WIP)
!


    call checkpoint_query ('/result/symmetry',sym_info)
    if (.not.sym_info) then
       print*, "WARNING: old CHECKPOINT, assuming NOSYM"
       maxopr = 0
       inversion = 1
    else
       call checkpoint_read ('/result/symmetry/maxopr',idata=maxopr)
       call checkpoint_read ('/result/symmetry/inversion',&
            idata=inversion)
    end if

    small_component = rs_input%force_small_ao

!   read general AO basis information
    if (maxopr == 0) then
!     SCF was done in C1
      if (.not. small_component) then
!       read the LC basis
        group = '/input/aobasis/1/'
      else
!       read the SC basis
        group = '/input/aobasis/2/'
      end if
    else
!     we need the AO basis expanded from the symmetry-packed form
      if (.not. small_component) then
!       read the LC basis
        group = '/input/aobasis/3/'
      else
!       read the SC basis
        group = '/input/aobasis/4/'
      end if
    end if

    call checkpoint_read(group//'angular',idata=ao_basis%basis_angular)
    call checkpoint_read(group//'n_ao',idata=nao)
    call checkpoint_read(group//'n_shells',idata=nshells)
    allocate (orbmom(nshells))
    allocate (n_prim(nshells))
    allocate (n_cont(nshells))
    allocate (coord(3*nshells))
!   read information that is defined for each shell
    call checkpoint_read (group//'orbmom',idata=orbmom)
    call checkpoint_read (group//'n_prim',idata=n_prim)
    call checkpoint_read (group//'n_cont',idata=n_cont)
    call checkpoint_read (group//'center',rdata=coord)
!   values of the exponents and coefficients for all shells (packed into a 1-d arrays)
    npriexp = sum(n_prim)
    allocate (priexp(npriexp))
    allocate (priccf(npriexp))
    call checkpoint_read (group//'exponents',rdata=priexp)
    call checkpoint_read (group//'contractions',rdata=priccf)

    ao_basis%nshells       = nshells
    ao_basis%nao           = nao
    nullify(ao_basis%gtos)
    allocate(ao_basis%gtos(nshells))
    iPrimStart = 1
    do ishell = 1, nshells
      ao_basis%gtos(ishell)%orb_momentum = orbmom(ishell)
      ao_basis%gtos(ishell)%n_primitives = n_prim(ishell)
      nullify  (ao_basis%gtos(ishell)%exponent) 
      allocate (ao_basis%gtos(ishell)%exponent(n_prim(ishell)))
      nullify  (ao_basis%gtos(ishell)%coefficient)
      allocate (ao_basis%gtos(ishell)%coefficient(n_prim(ishell)))
      iPrimEnd   = iPrimStart+n_prim(ishell)-1
      ao_basis%gtos(ishell)%exponent   = priexp(iPrimStart:iPrimEnd)
      ao_basis%gtos(ishell)%coefficient= priccf(iPrimStart:iPrimEnd)
      ao_basis%gtos(ishell)%coord = coord(ishell*3-2:3*ishell)
!     In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
      call compress_unc_shell (ao_basis%gtos(ishell),ierr)
      iPrimStart = iPrimStart + n_prim(ishell)
    end do

    deallocate(orbmom)
    deallocate(n_prim)
    deallocate(n_cont)
    deallocate(coord)
    deallocate(priexp)
    deallocate(priccf)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine check_data_consistency(ao_basis,mo)
    type(basis_set_info_t), intent(in) :: ao_basis
    type(cmo_t), intent(in) :: mo
    integer:: ierr
    integer:: inconsis

!   gosia fixme, after excorr_dirac.F
    if (sum(mo%n_ao) < ao_basis%nao) then
       ierr = INCONSIS
       return
    else if (ao_basis%nao > sum(mo%n_ao)) then
       print*,"WARNING: assuming Levy-Leblond Hamiltonian"
       print*,"Reading fewer basis functions than MO-coefficients"
       !allocate (ao_map(2,nao_bas))
       !forall (i=1:nao_bas) ao_map(:,i)=i
       !call compress_mo (cspinor,ao_map)
    end if


  end subroutine


#ifdef VAR_MPI
! ------------------------------------------------------------------------------
  subroutine synchronize_data_from_dirac(mol, cmo)
    type(mol_t), intent(inout):: mol
    type(cmo_t), intent(inout):: cmo
    integer:: len_hamiltonian, len_xc_functional
    integer:: nmo,nao,nz,ncoeff
!
!   this subroutine syncs data stored in cmo, mol, and rs_wf DTs, 
!   if any new variable has been added to these DTs, please also add it here
!

    if (mpi_is_master()) then
      if (allocated(rs_wf%hamiltonian)) then
        len_hamiltonian = len_trim(rs_wf%hamiltonian)
      else
        len_hamiltonian = 0
      end if
      if (allocated(rs_wf%xc_functional)) then
        len_xc_functional = len_trim(rs_wf%xc_functional)
      else
        len_xc_functional = 0
      end if
    end if
    call interface_mpi_bcast(len_hamiltonian,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_xc_functional,1,irank_master,global_communicator)

    call interface_mpi_bcast(mol%n_atoms,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      allocate(mol%symbols(mol%n_atoms))
      allocate(mol%geometry(3, mol%n_atoms))
      allocate(mol%nuc_charge(mol%n_atoms))
      if (len_hamiltonian > 0) allocate(character(len_hamiltonian):: rs_wf%hamiltonian)
      if (len_xc_functional > 0) allocate(character(len_xc_functional):: rs_wf%xc_functional)
    end if

    call interface_mpi_bcast(rs_wf%hamiltonian,len_hamiltonian,irank_master,global_communicator)
    call interface_mpi_bcast(rs_wf%xc_functional,len_xc_functional,irank_master,global_communicator)

    call interface_mpi_bcast(mol%symbols,mol%n_atoms,irank_master,global_communicator)
    call interface_mpi_bcast(mol%geometry,3*mol%n_atoms,irank_master,global_communicator)
    call interface_mpi_bcast(mol%nuc_charge,mol%n_atoms,irank_master,global_communicator)

    call interface_mpi_bcast(cmo%n_ao,2,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%n_mo,2,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%n_fsym,1,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%nz,1,irank_master,global_communicator)
    nmo = sum(cmo%n_mo)
    nao = sum(cmo%n_ao)
    nz  = cmo%nz
    ncoeff = nmo*nao*nz

    if (.not.mpi_is_master()) then
      call cmo%alloc_all()
    end if
    call interface_mpi_bcast(cmo%energy,nmo,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%coeff,ncoeff,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%coeff_q,ncoeff,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%symmetry,nmo,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%index,nmo,irank_master,global_communicator)

  end subroutine

! ------------------------------------------------------------------------------
  subroutine synchronize_dmat(dmat)
    type(dmat_t), intent(inout):: dmat
    integer:: n

    call interface_mpi_bcast(dmat%mat_dim,1,irank_master,global_communicator)
    call interface_mpi_bcast(dmat%mat_nz,1,irank_master,global_communicator)
    call interface_mpi_bcast(dmat%ih_sym,1,irank_master,global_communicator)
    call interface_mpi_bcast(dmat%tr_sym,1,irank_master,global_communicator)
    call interface_mpi_bcast(dmat%irep,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      call dmat%alloc()
    end if

    n = dmat%mat_dim*dmat%mat_dim*dmat%mat_nz
    call interface_mpi_bcast(dmat%elements,n,irank_master,global_communicator)

  end subroutine

! ------------------------------------------------------------------------------
  subroutine synchronize_grid_data(grid)
    type(grid_t), intent(inout):: grid
    integer:: len_gfn_inp, len_gfn_out

    if (mpi_is_master()) then
      if (allocated(grid%file_inp)) then
        len_gfn_inp = len_trim(grid%file_inp)
      else
        len_gfn_inp = 0
      end if
      if (allocated(grid%file_out)) then
        len_gfn_out = len_trim(grid%file_out)
      else
        len_gfn_out = 0
      end if
    end if
    call interface_mpi_bcast(len_gfn_inp,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gfn_out,1,irank_master,global_communicator)

    call interface_mpi_bcast(grid%id,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%ndim,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%nr_points_total,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%nr_points,3,irank_master,global_communicator)
    call interface_mpi_bcast(grid%origin,3,irank_master,global_communicator)
    call interface_mpi_bcast(grid%extents,6,irank_master,global_communicator)
    call interface_mpi_bcast(grid%spacing,3,irank_master,global_communicator)
    call interface_mpi_bcast(grid%status,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%input,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%export,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%typ,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%weights,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%margin,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%origin_known,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%extents_known,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%nr_points_total_known,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%nr_points_known,3,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%spacing_known,3,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      if (len_gfn_inp > 0) allocate(character(len_gfn_inp):: grid%file_inp)
      if (len_gfn_out > 0) allocate(character(len_gfn_out):: grid%file_out)
    end if

    call interface_mpi_bcast(grid%file_inp,len_gfn_inp,irank_master,global_communicator)
    call interface_mpi_bcast(grid%file_out,len_gfn_out,irank_master,global_communicator)

  end subroutine synchronize_grid_data

! ------------------------------------------------------------------------------
  subroutine synchronize_grid_function_data(grid_function)
    type(grid_function_t), intent(inout):: grid_function
    integer:: len_gffn_rsp, len_gffn_inp, len_gffn_out, len_gf_name, len_gf_typ

    if (mpi_is_master()) then
      if (allocated(grid_function%rspvec_file)) then
        len_gffn_rsp = len_trim(grid_function%rspvec_file)
      else
        len_gffn_rsp = 0
      end if
      if (allocated(grid_function%file_inp)) then
        len_gffn_inp = len_trim(grid_function%file_inp)
      else
        len_gffn_inp = 0
      end if
      if (allocated(grid_function%file_out)) then
        len_gffn_out = len_trim(grid_function%file_out)
      else
        len_gffn_out = 0
      end if
      if (allocated(grid_function%name)) then
        len_gf_name = len_trim(grid_function%name)
      else
        len_gf_name = 0
      end if
      if (allocated(grid_function%gf_type)) then
        len_gf_typ = len_trim(grid_function%gf_type)
      else
        len_gf_typ = 0
      end if
    end if
    call interface_mpi_bcast(len_gffn_rsp,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gffn_inp,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gffn_out,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gf_name,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gf_typ,1,irank_master,global_communicator)

    call interface_mpi_bcast(grid_function%id,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%status,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%input,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid_function%export,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%nuc,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%wf_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%rspvec_id,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%lao,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%lao_skip,3,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%tensor_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%max_ao_g_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%max_ao_m_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%id_assigned_grid,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%purpose,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%modify,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%cartesian_power,3,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%up_scale,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%down_scale,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      if (len_gffn_rsp > 0) allocate(character(len_gffn_rsp):: grid_function%rspvec_file)
      if (len_gffn_inp > 0) allocate(character(len_gffn_inp):: grid_function%file_inp)
      if (len_gffn_out > 0) allocate(character(len_gffn_out):: grid_function%file_out)
      if (len_gf_name > 0) allocate(character(len_gf_name):: grid_function%name)
      if (len_gf_typ > 0) allocate(character(len_gf_typ):: grid_function%gf_type)
    end if

    call interface_mpi_bcast(grid_function%rspvec_file,len_gffn_rsp,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%file_inp,len_gffn_inp,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%file_out,len_gffn_out,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%name,len_gf_name,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%gf_type,len_gf_typ,irank_master,global_communicator)

  end subroutine synchronize_grid_function_data

#endif

! ------------------------------------------------------------------------------
  subroutine dmat_transform_mo2ao(cmo, dmat_mo, dmat_ao)
    type(cmo_t), intent(in):: cmo
    type(dmat_t), intent(inout):: dmat_ao
    type(dmat_t), intent(in):: dmat_mo
    integer:: iprint
    integer:: nao, nmo, nz, nrao, ncao, nrmo, ncmo, irep
    integer:: i2basx(2, 2), i2orbx(2, 2), icmoq(2), ipqtoq(2, 0:7)
    integer:: i,j
!
!   dmat_ao(mu,nu)=cmo(mu,p)*dmat_mo(p,q)*cmo^{-1}(q,nu)
!   gosia fixme: no symmetry
! ------------------------------------------------------------------------------

    ! gosia fixme - this causes "catastrophic error: **Internal compiler error: internal abort** Please report this error along with the circumstances in which it occurred in a Software Problem Report.  Note: File and line given may not be explicit cause of this error."
    !if (   .not.dmat_ao%associated_elements  &
    !  .and..not.dmat_mo%associated_elements  &
    !  .and..not.cmo%associated_cmo_coeff) return

    iprint = 0

    nao  = sum(cmo%n_ao)
    nmo  = sum(cmo%n_mo)
    nz   = cmo%nz

    i = 1
    j = 1
      nrao = cmo%n_ao(i)
      ncao = cmo%n_ao(j)
      nrmo = cmo%n_mo(i)
      ncmo = cmo%n_mo(j)
      irep = 0
      call qtrans90('MOAO','S',0.D0,nrao,ncao,nrmo,ncmo,&
                    dmat_ao%elements,nao,nao,nz,(/1,2,3,4/),& 
                    dmat_mo%elements,nmo,nmo,nz,(/1,2,3,4/),& 
                    cmo%coeff,nrao,nrmo,nz,(/1,2,3,4/), &
                    cmo%coeff,ncao,ncmo,nz,(/1,2,3,4/), &
                    iprint)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine construct_perturbed_dmat(dmat,bmat,  &
                                      ivec, &
                                      rsp_info,&
                                      cmo)
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbham.h"

    integer, intent(in)  :: ivec
    type(response_info_t), intent(in) :: rsp_info
    real(kind=8), intent(out) :: dmat(ntbas(0)*ntbas(0), nz)
    real(kind=8), intent(in) :: bmat(norbt, norbt, nz)
    real(kind=8), intent(in) :: cmo(*)

    integer :: irep,iopsy,itim
    integer :: iprint
    integer :: i1, i2, i, j, iz, ns, is, js, ni, ii, ji
    character(2) ::typ
    logical :: debug_print
!

    iprint = 0  ! needed for qtrans

!   construct half dmat = Ci*bmat*dag(Ca)

    dmat(1:ntbas(0)*ntbas(0), 1:nz) = 0.0d0

    call rsp_info%vectors(ivec)%get_symmetry_info(irep,&
            iopsy,itim,typ)

    write(*,*) 'in construct_perturbed_dmat: typ  = ', typ
    write(*,*) 'in construct_perturbed_dmat: itim = ', itim
    write(*,*) 'in construct_perturbed_dmat: irep = ', irep
    write(*,*) 'in construct_perturbed_dmat: iopsy = ', iopsy

    do i1 = 1, nfsym
      i2 = mod(i1+iopsy,2) + 1

!     second index: inactive index:
      ni = nish(i2)
      ii = iorb(i2) + npsh(i2) + 1
      ji = icmoq(i2) + npsh(i2)*nfbas(i2,0)+1

!     first index: secondary index:
      if (typ(1:1) == 'E') then
!        .. secondary index is electronic:
         ns = nesh(i1)
         is = iorb(i1) + npsh(i1) + 1
         js = icmoq(i1) + npsh(i1)*nfbas(i1, 0) + 1
      else if (typ(1:1) == 'P') then
!        .. secondary index is positronic:
         ns = npsh(i1)
         is = iorb(i1) + 1
         js = icmoq(i1) + 1
      end if

      if (ni > 0 .and. ns > 0) then

        call qtrans90('MOAO', &
                      'S',    &
                      0.0d0,  &
                      nfbas(i1,0), nfbas(i2,0), &
                      ns, ni, &
!                     ------------------------------------------------------
                      dmat(i2basx(i1,i2)+1, 1),   &
                      ntbas(0), ntbas(0), nz,     &
                      ipqtoq(1, irep),                  &
!                     ------------------------------------------------------
                      bmat(is, ii, 1),                    &
                      norbt, norbt, nz, &
                      ipqtoq(1, irep),                  &
!                     ------------------------------------------------------
                      cmo(js),                                  &
                      nfbas(i1,0), norb(i1), nz,&
                      ipqtoq(1, 0),                     &
!                     ------------------------------------------------------
                      cmo(ji),                                  &
                      nfbas(i2,0), norb(i2), nz,&
                      ipqtoq(1, 0),                     &
!                     ------------------------------------------------------
                      iprint)
      end if

     end do

  end subroutine





! ------------------------------------------------------------------------------
  subroutine construct_perturbed_dmat_lao_e(dmat,bmat,  &
                                      cmo)
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbham.h"

    real(kind=8), intent(inout) :: dmat(ntbas(0)*ntbas(0), nz)
    real(kind=8), intent(inout) :: bmat(norbt, norbt, nz)
    real(kind=8), intent(in) :: cmo(*)

    integer :: irep,iopsy,itim
    integer :: iprint
    integer :: i1, i2, i, j, iz, ns, is, js, ni, ii, ji
    character(2) ::typ
    logical :: debug_print
!

    iprint = 0
    irep = 0
    typ = 'P-'

    do i1 = 1, 1
      i2 = 1

!     second index: inactive index:
      ni = nish(i2)
      ii = iorb(i2) + npsh(i2) + 1
      ji = icmoq(i2) + npsh(i2)*nfbas(i2,0)+1

!     first index: secondary index:
!        .. secondary index is electronic:
         ns = nesh(i1)
         is = iorb(i1) + npsh(i1) + 1
         js = icmoq(i1) + npsh(i1)*nfbas(i1, 0) + 1

      if (ni > 0 .and. ns > 0) then

        call qtrans90('MOAO', &
                      'S',    &
                      1.0d0,  &
                      nfbas(i1,0), nfbas(i2,0), &
                      ns, ni, &
!                     ------------------------------------------------------
                      dmat(i2basx(i1,i2)+1, 1),   &
                      ntbas(0), ntbas(0), nz,     &
                      ipqtoq(1, irep),                  &
!                     ------------------------------------------------------
                      bmat(is, ii, 1),                    &
                      norbt, norbt, nz, &
                      ipqtoq(1, irep),                  &
!                     ------------------------------------------------------
                      cmo(js),                                  &
                      nfbas(i1,0), norb(i1), nz,&
                      ipqtoq(1, 0),                     &
!                     ------------------------------------------------------
                      cmo(ji),                                  &
                      nfbas(i2,0), norb(i2), nz,&
                      ipqtoq(1, 0),                     &
!                     ------------------------------------------------------
                      iprint)
      end if

     end do

  end subroutine
  subroutine construct_perturbed_dmat_lao_p(dmat,bmat,  &
                                      cmo)
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbham.h"

    real(kind=8), intent(inout) :: dmat(ntbas(0)*ntbas(0), nz)
    real(kind=8), intent(inout) :: bmat(norbt, norbt, nz)
    real(kind=8), intent(in) :: cmo(*)

    integer :: irep,iopsy,itim
    integer :: iprint
    integer :: i1, i2, i, j, iz, ns, is, js, ni, ii, ji
    character(2) ::typ
    logical :: debug_print
!

    iprint = 0
    irep = 0
    typ = 'P-'

    do i1 = 1, 1
      i2 = 1

!     second index: inactive index:
      ni = nish(i2)
      ii = iorb(i2) + npsh(i2) + 1
      ji = icmoq(i2) + npsh(i2)*nfbas(i2,0)+1

!     first index:
!        .. secondary index is positronic:
         ns = npsh(i1)
         is = iorb(i1) + 1
         js = icmoq(i1) + 1

      if (ni > 0 .and. ns > 0) then

        call qtrans90('MOAO', &
                      'S',    &
                      1.0d0,  &
                      nfbas(i1,0), nfbas(i2,0), &
                      ns, ni, &
!                     ------------------------------------------------------
                      dmat(i2basx(i1,i2)+1, 1),   &
                      ntbas(0), ntbas(0), nz,     &
                      ipqtoq(1, irep),                  &
!                     ------------------------------------------------------
                      bmat(is, ii, 1),                    &
                      norbt, norbt, nz, &
                      ipqtoq(1, irep),                  &
!                     ------------------------------------------------------
                      cmo(js),                                  &
                      nfbas(i1,0), norb(i1), nz,&
                      ipqtoq(1, 0),                     &
!                     ------------------------------------------------------
                      cmo(ji),                                  &
                      nfbas(i2,0), norb(i2), nz,&
                      ipqtoq(1, 0),                     &
!                     ------------------------------------------------------
                      iprint)
      end if

     end do

  end subroutine


!
! ------------------------------------------------------------------------------
  subroutine pertdenmat_sym(mat_dim,mat_nz,mat, f)

    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: mat_nz
    real(8), intent(inout) :: mat(mat_dim, mat_dim, mat_nz)
    real(8), intent(in)    :: f(4)
    integer :: i, j, iz
    real(8) :: fullmat
!
!   after dmoli3 in prp/pamrvc.F and fulmat in gp/gpsaue.F
!
    do iz = 1, mat_nz
      do j = 1, mat_dim
        do i = j, mat_dim
          fullmat = mat(i, j, iz) + f(iz)*mat(j, i, iz)
          mat(i, j, iz) =  fullmat
          mat(j, i, iz) =  f(iz)*fullmat
        end do
      end do
    end do

  end subroutine



! ------------------------------------------------------------------------------
  subroutine setup_orbital_string(cmo)
#include "dcborb.h" !nish, nash
#include "maxorb.h" !mxprim
    type(cmo_t), intent(inout) :: cmo
    integer:: i, ifsym, ioff

!   this routine will set up occupation vector
!   based on**VISUAL
!            .OCCUPATION
!   if this is not given, then it will occupy all
!   "occupied" electrons with 1.0

    if (.not. cmo%data_is_read) return
    if (.not. cmo%allocated_occupation()) return

    cmo%occupation = 0.0d0

    ioff = 0
    do ifsym = 1, cmo%n_fsym
      ioff = ioff+cmo%n_po(ifsym)
      if (rs_input%set_occupation) then
        if (.not. rs_input%allocated_occupation()) return
        print *, 'Orbital occupations will be modified according to .OCCUPATION'
        do i = 1, (cmo%n_mo(ifsym) - cmo%n_po(ifsym))
          cmo%occupation(ioff+i) = rs_input%occupation(i, ifsym)
        end do
      else
        if (nash(ifsym) > 0) then
          print *, 'ERROR: with open shells you should use .OCCUPATION'
          print *, '       otherwise occupation will be wrong'
          print *, '       so we better stop here'
          stop 1
        end if
        do i = 1, nish(ifsym)
          cmo%occupation(ioff+i) = 1.0d0
        end do
      end if
      ioff = cmo%n_mo(ifsym)
      if (cmo%n_ip(ifsym) == 0) then
        cmo%n_ip(ifsym) = nish(ifsym)
      end if
    end do

    ! todo: most of cmo data is first read from checkpoint,
    ! but it lacks information on n_ip and occupation, which are set here;
    ! this is to catch any inconsistency: 
    !call check_cmo_consistency(cmo)
    cmo%data_is_verified = .true.

  end subroutine


! ------------------------------------------------------------------------------
  subroutine setup_mo_offsets(cmo)
    type(cmo_t), intent(inout) :: cmo
    integer:: i, ifsym, ioff
!
!   order of offsets for MOs:
!   ioffs(1) - start of gerade secondary negative energy
!   ioffs(2) - start of gerade inactive positive energy
!   ioffs(3) - start of gerade secondary positive energy
!   ioffs(4) - start of ungerade secondary negative energy
!   ioffs(5) - start of ungerade inactive positive energy
!   ioffs(6) - start of ungerade secondary positive energy
!

    if (.not. cmo%data_is_verified) return

    cmo%ioffs = 0

    cmo%ioffs(1) = 1
    cmo%ioffs(2) = cmo%ioffs(1) + cmo%n_po(1)
    cmo%ioffs(3) = cmo%ioffs(2) + cmo%n_ip(1)

    if (cmo%n_fsym == 2) then
      cmo%ioffs(4) = 1 + cmo%n_mo(1)
      cmo%ioffs(5) = cmo%ioffs(4) + cmo%n_po(2)
      cmo%ioffs(6) = cmo%ioffs(5) + cmo%n_ip(2)
    end if

  end subroutine

  subroutine get_xyzcomp_from_rspvec(grid_function, ivec)
    type(grid_function_t), intent(inout):: grid_function
    integer,intent(out),optional::ivec
    character(len=:),allocatable :: f_name
    type(response_info_t):: rsp_info_this
    integer:: rspvec_id

    f_name = grid_function%rspvec_file
    call rsp_info_this%read_from_checkpoint(f_name)
    if ( rs_input%print_level > 2 ) call rsp_info_this%display
    rspvec_id = grid_function%rspvec_id
 
    if (rsp_info_this%vectors(rspvec_id)%get_label() == 'LAO-XRM1H1') then
      if (present(ivec)) ivec = 1
      grid_function%lao_comp = 'X'
    else if (rsp_info_this%vectors(rspvec_id)%get_label() == 'LAO-YRM1H1') then
      if (present(ivec)) ivec = 2
      grid_function%lao_comp = 'Y'
    else if (rsp_info_this%vectors(rspvec_id)%get_label() == 'LAO-ZRM1H1') then
      if (present(ivec)) ivec = 3
      grid_function%lao_comp = 'Z'
    end if

    call rsp_info_this%erase()
  end

! ------------------------------------------------------------------------------
  subroutine visual_input_show_setup()

    write(*, '(/,2x, a)') '*** Visualization module ***' 
    if ( rs_input%print_level > 2 ) then
      write(*, *) "DEBUG: rs_input:"
#ifdef VAR_MPI
      call rs_input%display(my_rank=get_my_rank())
#else
      call rs_input%display()
#endif
    end if
 
  end subroutine


! ------------------------------------------------------------------------------
  subroutine verify_input_keywords()

  end subroutine verify_input_keywords


end module
