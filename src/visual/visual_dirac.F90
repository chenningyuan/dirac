module visual_dirac

  use checkpoint
  use visual_datatypes, only: mol_t, ao_t, cmo_t, grid_t, grid_function_t
  use visual_cfg

  use visual_mpi
#ifdef VAR_MPI
   use interface_to_mpi
   use integer_kind_mpilib
#endif


  implicit none

  public visual_initialize_and_sync_nodes
  public visual_initialize_data_from_dirac

  public dirac_dmat_init

#ifdef VAR_MPI
  public synchronize_grid_data
  public synchronize_grid_function_data
  public synchronize_data_from_dirac
#endif


  private


contains

  subroutine dirac_dmat_init(wf_order, cmo, dmat_ao, dmat_mo)
!   ----------------------------------------------------------------------------
    integer, intent(in) :: wf_order
    type(cmo_t), intent(in) :: cmo
    real(8), optional, pointer, intent(out) :: dmat_ao(:,:,:)
    real(8), optional, pointer, intent(out) :: dmat_mo(:,:,:)
    integer :: nto
    real(8), allocatable :: occupation(:)
!   ----------------------------------------------------------------------------

!    setup density matrix in AO basis from occupation numbers (user-controlled)
!    get orbital string

     nto = sum(cmo%n_mo) ! total number of orbitals
     allocate(occupation(nto))
     call setup_orbital_string(cmo,occupation)

!    get CMOs

!    !construct dmats

     !dmat = 0.0d0

     !do m = 1, norbt
     !   D_mo%elms((m-1)*norbt+m) = occupation(m)
     !end do

     !D_mo%ih_sym = 1
     !D_mo%tr_sym = 1
     !D_mo%irep   = 0

     !if (wf_order == 0) dmat = C*(D_mo*dag(C))
!
!    !clean
!
     deallocate(occupation)

  end subroutine dirac_dmat_init


  subroutine setup_orbital_string(cmo,occupation)
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
    type(cmo_t), intent(in) :: cmo
!   ----------------------------------------------------------------------------

    ! this routine will set up occupation vector
    ! based on**VISUAL
    !          .OCCUPATION
    ! if this is not given, then it will occupy all
    ! "occupied" electrons with 1.0
    real(8), intent(out):: occupation(:)
    integer:: i, ifsym, ioff
!   ----------------------------------------------------------------------------

    occupation = 0.0d0

    ioff = 0
    do ifsym = 1, cmo%n_fsym
      ioff = ioff+npsh(ifsym)
      !ioff = ioff+cmo%n_po(ifsym)
      if (visual_cfg_use_orbital_string) then
        !do i = 1, cmo%n_mo(ifsym) - cmo%n_po(ifsym)
        do i = 1, nesh(ifsym)
          occupation(ioff+i) = visual_cfg_occupation(i, ifsym)
        end do
      else
        !if (cmo%n_mo(ifsym)-cmo%n_po(ifsym) > 0) then
        if (nash(ifsym) > 0) then
          print *, 'ERROR: with open shells you should use .OCCUPATION'
          print *, '       otherwise occupation will be wrong'
          print *, '       so we better stop here'
          stop 1
        end if
        !do i = 1, cmo%n_mo(ifsym)
        do i = 1, nish(ifsym)
          occupation(ioff+i) = 1.0d0
        end do
      end if
      ioff = cmo%n_mo(ifsym)
    end do
    write(*,'(a,i6)') 'NEW OCCUPATION SIZE: ', sum(cmo%n_mo)
    !do i = 1, sum(cmo%n_mo)
    do i = 1, norbt
      write(*,'(a,i6,e16.8)') 'NEW OCCUPATION: ', i, occupation(i)
    end do

!!!!!!!!!!!!!!!
    write(*, *)                                                                                          
    write(*, *) 'OLD CMO:'
    write(*, *) 'OLD n_ao =         ', nfbas
    write(*, *) 'OLD n_mo =         ', NORB
    write(*, *) 'OLD n_po =         ', NPSH
    write(*, *) 'OLD n_fsym =       ', nfsym
    write(*, *) 'OLD nz   =         ', nz
    write(*, *)                                                                                          



!!!!!!!!!!!!!!


  end subroutine

  subroutine visual_initialize_and_sync_nodes(hello)
    use dirac_cfg, only: dirac_cfg_visual
    use visual_mpi
#ifdef VAR_MPI
    use interface_to_mpi
    use integer_kind_mpilib
#endif
#include "dirac_partask.h"  
! legacy
#include "dcbbas.h"
#include "dgroup.h"

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
    ! legacy
    mat_ao_dim = ntbas(0)
    mat_ao_nz = nz
    mat_ao_irep = 0

#ifdef VAR_MPI
    if (is_parallel()) call sync_visual_cfg()
#endif

  end subroutine visual_initialize_and_sync_nodes


  subroutine visual_initialize_data_from_dirac(mol, ao, mo, print_level)

    type(mol_t), intent(out):: mol
    type(ao_t),  intent(out):: ao
    type(cmo_t), intent(out):: mo
    integer, optional :: print_level
    character*(:),allocatable:: group
    logical:: tobe_dirac

    group = '/input/molecule/'
    call checkpoint_query(group, tobe_dirac)
    if (tobe_dirac) then
      call read_mol_from_dirac(mol, group)
      if (present(print_level)) then
        if (print_level > 0) then
          write(*, *) 'Molecule data read from DIRAC CHECKPOINT'
          call mol%display
        end if
      end if

      !call write_to_rs_checkpoint('/input/molecule', mol)
      !if (print_level > 0) &
      !  write(*, *) 'Molecule data written to DIRAC RS-CHECKPOINT'
    end if

!   gosia FIXME - the code is validated for scf only
    group = '/result/wavefunctions/scf/'
    call checkpoint_query(group//'mobasis/nz', tobe_dirac)
    if (tobe_dirac) then
      call read_ao_mo_from_dirac(ao, mo, group)
      !if (present(print_level)) then
      !  if (print_level > 0) then
!g          write(*, *) 'AO data read from DIRAC CHECKPOINT'
!g          !call ao%display
          write(*, *) 'MO data read from DIRAC CHECKPOINT'
          call mo%display
      !  end if
      !end if
    end if

  end subroutine visual_initialize_data_from_dirac


  subroutine read_mol_from_dirac(mol, group)
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

    mol%geometry = reshape(temp, (/3, mol%n_atoms/))

    deallocate(temp)

  end subroutine read_mol_from_dirac


  subroutine read_ao_mo_from_dirac(ao, cmo, group)
!   ----------------------------------------------------------------------------
    type(ao_t), intent(out):: ao
    type(cmo_t), intent(inout) :: cmo
    character(len=*),intent(in) :: group
    logical:: tobe_dirac
    integer:: nao,npo,nmo,nto,ncoef,nfsym,nz
    integer:: n_basis(2),n_po(2),n_mo(2)
    real(8):: total_energy
    integer, allocatable   :: ibeig(:), index(:)
    real(8), allocatable   :: eig(:), mo_coeff(:), mo_coeff_q(:,:,:)
    integer:: i
!
!   this subroutine reads data from DIRAC CHECKPOINT files
!   and fills the "ao" and "mo" derived type
!   now - based on how things are done in dircmo.F90 and exacorr
!   gosia, summer 2023
!   ----------------------------------------------------------------------------
!
!   read MO data
!
    n_basis = 0
    n_mo = 0
    n_po = 0
 
    call checkpoint_read(group//'mobasis/nz',idata=nz)
    call checkpoint_read(group//'mobasis/n_fsym',idata=nfsym)
    call checkpoint_read(group//'mobasis/n_basis',idata=n_basis)
    call checkpoint_read(group//'mobasis/n_mo',idata=n_mo)
    call checkpoint_read(group//'mobasis/n_po',idata=n_po)
    call checkpoint_read(group//'energy',rdata=total_energy)

    nao = sum(n_basis)
    npo = sum(n_po)
    nmo = sum(n_mo)
    nto = npo+nmo
    ncoef = nz * sum(n_basis*n_mo)

    cmo%n_ao = n_basis
    cmo%n_mo = n_mo
    cmo%nz   = nz
    cmo%total_energy = total_energy

    call cmo%allocate_all()

    call checkpoint_query(group//'mobasis/orbitals',tobe_dirac)
    if (tobe_dirac) then
      allocate(mo_coeff(ncoef))
      mo_coeff = 0.0d0
      call checkpoint_read(group//'/mobasis/orbitals',rdata=mo_coeff)
    ! put in cmo in quaternion format
      allocate(mo_coeff_q(nao,nto,nz))
      mo_coeff_q = 0.0d0
      mo_coeff_q = reshape(mo_coeff,(/nao,nto,nz/))
      cmo%coeff  = mo_coeff_q(:,npo+1:nto,:)
      deallocate(mo_coeff)
      deallocate(mo_coeff_q)
    end if

    call checkpoint_query(group//'mobasis/eigenvalues',tobe_dirac)
    if (tobe_dirac) then
      allocate(eig(nto))
      eig = 0.0d0
      call checkpoint_read(group//'/mobasis/eigenvalues',rdata=eig)
      cmo%energy   = eig(npo+1:nto)
      deallocate(eig)
    end if

    call checkpoint_query(group//'mobasis/symmetry',tobe_dirac)
    if (tobe_dirac) then
      allocate(ibeig(nto))
      ibeig = 0
      call checkpoint_read(group//'/mobasis/symmetry',idata=ibeig)
      cmo%symmetry = ibeig(npo+1:nto)
      deallocate(ibeig)
    end if

    ! fill in index values
    allocate(index(nmo))
    do i = 1, nmo
       index(i) = i ! we set the index compatible with the fact that we store all mo's
    end do
    cmo%index    = index(npo+1:nto)
    deallocate(index)

!   switch to unrestricted spinor picture
    !call convert_mo_to_complex(cmo)

  end subroutine read_ao_mo_from_dirac


#ifdef VAR_MPI
  subroutine synchronize_data_from_dirac(mol, ao, cmo)
    !
    ! this subroutine syncs data stored in ao and mol dts, 
    ! which are read from DIRAC CHECKPOINT
    ! if any new variable has been added to these dts, please also add it here
    !
    type(mol_t), intent(inout):: mol
    type(ao_t), intent(inout):: ao
    type(cmo_t), intent(inout):: cmo
    integer:: n,m

    call interface_mpi_bcast(mol%n_atoms,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      allocate(mol%symbols(mol%n_atoms))
      allocate(mol%geometry(3, mol%n_atoms))
      allocate(mol%nuc_charge(mol%n_atoms))
    end if

    call interface_mpi_bcast(mol%symbols,mol%n_atoms,irank_master,global_communicator)
    call interface_mpi_bcast(mol%geometry,3*mol%n_atoms,irank_master,global_communicator)
    call interface_mpi_bcast(mol%nuc_charge,mol%n_atoms,irank_master,global_communicator)

    ! gosia todo/fixme what ao data do you need??
    call interface_mpi_bcast(ao%n_ao,1,irank_master,global_communicator)

    !write(*, *) 'in synchronize_data_from_dirac, mol on node: ', get_my_rank()
    !call mol%display()

    !cmo data
    call interface_mpi_bcast(cmo%n_ao,2,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%n_mo,2,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%n_fsym,1,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%nz,1,irank_master,global_communicator)
    n = sum(cmo%n_mo)
    m = sum(cmo%n_ao)
    if (.not.mpi_is_master()) then
      nullify(cmo%energy)
      allocate(cmo%energy(n))
      cmo%energy(1:n) = 0.0d0
  
      nullify(cmo%coeff)
      allocate(cmo%coeff(m,n,cmo%nz))
      cmo%coeff(1:m,1:n,1:cmo%nz) = 0.0d0
  
      nullify(cmo%symmetry)
      allocate(cmo%symmetry(n))
      cmo%symmetry(1:n) = 0
  
      nullify(cmo%index)
      allocate(cmo%index(n))
      cmo%index(1:n) = 0
    end if
    call interface_mpi_bcast(cmo%energy,n,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%coeff,n*m*cmo%nz,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%symmetry,n,irank_master,global_communicator)
    call interface_mpi_bcast(cmo%index,n,irank_master,global_communicator)

  end subroutine synchronize_data_from_dirac


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
    call interface_mpi_bcast(grid%status,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%input,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%typ,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%weights,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%margin,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%export,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%origin_known,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%nr_points_total_known,1,irank_master,global_communicator)

    call interface_mpi_bcast(grid%nr_points,grid%ndim,irank_master,global_communicator)
    call interface_mpi_bcast(grid%origin,grid%ndim,irank_master,global_communicator)
    call interface_mpi_bcast(grid%spacing,grid%ndim,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%nr_points_known,grid%ndim,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid%spacing_known,grid%ndim,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
    !!  allocate(grid%points(grid%ndim, grid%nr_points_total))
    !!  grid%points = 0.0d0
      if (len_gfn_inp > 0) allocate(character(len_gfn_inp):: grid%file_inp)
      if (len_gfn_out > 0) allocate(character(len_gfn_out):: grid%file_out)
    end if

    call interface_mpi_bcast(grid%file_inp,len_gfn_inp,irank_master,global_communicator)
    call interface_mpi_bcast(grid%file_out,len_gfn_out,irank_master,global_communicator)

  end subroutine synchronize_grid_data


  subroutine synchronize_grid_function_data(grid_function)
    type(grid_function_t), intent(inout):: grid_function
    integer:: len_gffn_inp, len_gffn_out, len_gf_name

    if (mpi_is_master()) then
      !write(*,*) 'IN synchronize_grid_function_data: ', grid_function%file_out
      !write(*,*) 'IN synchronize_grid_function_data: ', len_trim(grid_function%file_out)
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
    end if
    call interface_mpi_bcast(len_gffn_inp,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gffn_out,1,irank_master,global_communicator)
    call interface_mpi_bcast(len_gf_name,1,irank_master,global_communicator)

    call interface_mpi_bcast(grid_function%id,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%id_assigned_grid,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%status,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%input,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%wf_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%tensor_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%max_ao_g_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%max_ao_m_order,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%purpose,1,irank_master,global_communicator)
    call interface_mpi_bcast_l0(grid_function%export,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      if (len_gffn_inp > 0) allocate(character(len_gffn_inp):: grid_function%file_inp)
      if (len_gffn_out > 0) allocate(character(len_gffn_out):: grid_function%file_out)
      if (len_gf_name > 0) allocate(character(len_gf_name):: grid_function%name)
    end if

    call interface_mpi_bcast(grid_function%file_inp,len_gffn_inp,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%file_out,len_gffn_out,irank_master,global_communicator)
    call interface_mpi_bcast(grid_function%name,len_gf_name,irank_master,global_communicator)

  end subroutine synchronize_grid_function_data

#endif

end module
