module visual_dirac

  use checkpoint
  use visual_datatypes, only: mol_t, ao_t, grid_t
  use visual_cfg

  use visual_mpi
#ifdef VAR_MPI
   use interface_to_mpi
   use integer_kind_mpilib
#endif


  implicit none

  public visual_initialize_nodes
#ifdef VAR_MPI
  public synchronize_grid_data
#endif

  private

  interface read_from_dirac
    module procedure read_from_dirac
  end interface

contains

  subroutine visual_initialize_nodes(mol, aobasis, hello)
    use dirac_cfg, only: dirac_cfg_visual
    use visual_mpi
#ifdef VAR_MPI
    use interface_to_mpi
    use integer_kind_mpilib
#endif
#include "dirac_partask.h"  
    type(mol_t), intent(out):: mol
    type(ao_t), intent(out):: aobasis
    logical, intent(in), optional :: hello

#ifdef VAR_MPI
    ! wake up workers
    if (is_parallel()) call dirac_parctl( VISUAL_PAR )
#endif

    ! initialize and sync data:
    ! * global data: "rs_input" + legacy (dirac_cfg* variables);
    !   the only public data is going to be "rs_input",
    !   and this should be known now (from reading DIRAC input)
    call initialize_visual_data(mol, aobasis)
    if (present(hello)) then
      if (hello) return
    end if

    ! * data read from DIRAC CHECKPOINT (mol, ao)
    if (mpi_is_master()) then
      call read_from_dirac(mol, aobasis, rs_input%print_level)
    end if

#ifdef VAR_MPI
    if (is_parallel()) call synchronize_data_from_dirac(mol, aobasis)
#endif

  end subroutine visual_initialize_nodes


  subroutine initialize_visual_data(mol, aobasis)
    use dirac_cfg, only: dirac_cfg_visual
    type(mol_t), intent(out):: mol
    type(ao_t), intent(out):: aobasis

#ifdef VAR_MPI
    if (is_parallel()) then
      if (get_my_rank() /= irank_master) then
        call sync_visual_cfg
        return
      end if
    end if
#endif

    ! then, the rest happens on master only

#ifdef VAR_MPI
    call sync_visual_cfg
#endif
  end subroutine initialize_visual_data


  subroutine read_from_dirac(mol, ao, print_level)
    integer, intent(in):: print_level

    type(mol_t), intent(out):: mol
    type(ao_t), intent(out):: ao
    character*(:),allocatable:: group
    logical:: tobe_dirac

    group = '/input/molecule/'
    call checkpoint_query(group, tobe_dirac)
    if (tobe_dirac) then
      call read_mol_from_dirac(mol, group)
      write(*, *) 'Molecule data read from DIRAC CHECKPOINT'
      !if (print_level > 0) call mol%display

      !call write_to_rs_checkpoint('/input/molecule', mol)
      !if (print_level > 0) &
      !  write(*, *) 'Molecule data written to DIRAC RS-CHECKPOINT'
    end if

    !group = '/input/aobasis/'
    !call checkpoint_query(group, tobe_dirac)
    !if (tobe_dirac) then
    !  call read_ao_from_dirac(ao, group)
    !  write(*, *) 'AO data read from DIRAC CHECKPOINT'
    !end if

  end subroutine read_from_dirac


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


  subroutine read_ao_from_dirac(ao, group)
    type(ao_t), intent(out) :: ao
    character(len=*),intent(in) :: group

  end subroutine read_ao_from_dirac


#ifdef VAR_MPI
  subroutine synchronize_data_from_dirac(mol, ao)
    !
    ! this subroutine syncs data stored in ao and mol, read from DIRAC CHECKPOINT
    ! if any new variable has been added to these dts, please also add it here
    !
    type(mol_t), intent(inout):: mol
    type(ao_t), intent(inout):: ao

    call interface_mpi_bcast(mol%n_atoms,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      allocate(mol%symbols(mol%n_atoms))
      allocate(mol%geometry(3, mol%n_atoms))
      allocate(mol%nuc_charge(mol%n_atoms))
    end if

    call interface_mpi_bcast(mol%symbols,mol%n_atoms,irank_master,global_communicator)
    call interface_mpi_bcast(mol%geometry,3*mol%n_atoms,irank_master,global_communicator)
    call interface_mpi_bcast(mol%nuc_charge,mol%n_atoms,irank_master,global_communicator)

    call interface_mpi_bcast(ao%n_ao,1,irank_master,global_communicator)

    write(*, *) 'in synchronize_data_from_dirac, mol on node: ', get_my_rank()
    call mol%display()

  end subroutine synchronize_data_from_dirac


  subroutine synchronize_grid_data(grid)
    type(grid_t), intent(inout):: grid
    integer:: len_gl, len_fn_inp, len_fn_out


    !if (mpi_is_master()) then
    !  !len_gl = len_trim(grid_labels%x_label%name)
    !  len_fn_inp = len_trim(grid%file_inp)
    !  len_fn_out = len_trim(grid%file_out)
    !end if
    !call interface_mpi_bcast(len_fn_inp,1,irank_master,global_communicator)
    !call interface_mpi_bcast(len_fn_out,1,irank_master,global_communicator)

    call interface_mpi_bcast(grid%id,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%ndim,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%nr_points_total,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%status,1,irank_master,global_communicator)
    call interface_mpi_bcast(grid%input,1,irank_master,global_communicator)
!    call interface_mpi_bcast(grid%typ,1,irank_master,global_communicator)
!    call interface_mpi_bcast(grid%weights,1,irank_master,global_communicator)
!    call interface_mpi_bcast(grid%margin,1,irank_master,global_communicator)
!    call interface_mpi_bcast_l0(grid%export,1,irank_master,global_communicator)
!    call interface_mpi_bcast_l0(grid%origin_known,1,irank_master,global_communicator)
!    call interface_mpi_bcast_l0(grid%nr_points_total_known,1,irank_master,global_communicator)
!
!    !call interface_mpi_bcast(grid%nr_points,grid%ndim,irank_master,global_communicator)
!    !call interface_mpi_bcast(grid%origin,grid%ndim,irank_master,global_communicator)
!    !call interface_mpi_bcast(grid%spacing,grid%ndim,irank_master,global_communicator)
!    !call interface_mpi_bcast_l0(grid%nr_points_known,grid%ndim,irank_master,global_communicator)
!    !call interface_mpi_bcast_l0(grid%spacing_known,grid%ndim,irank_master,global_communicator)
!
!    !if (.not.mpi_is_master()) then
!    !  allocate(grid%nr_points(grid%ndim))
!    !  allocate(grid%origin(grid%ndim))
!    !  allocate(grid%spacing(grid%ndim))
!    !  allocate(grid%nr_points_known(grid%ndim))
!    !  allocate(grid%spacing_known(grid%ndim))
!    !!  allocate(grid%points(grid%ndim, grid%nr_points_total))
!    !!  grid%points = 0.0d0
!    !!  if (len_fn_inp > 0) allocate(character(len_fn_inp):: grid%file_inp)
!    !!  if (len_fn_out > 0) allocate(character(len_fn_out):: grid%file_out)
!    !end if
!
!    !call interface_mpi_bcast(grid%file_inp,len_fn_inp,irank_master,global_communicator)
!    !call interface_mpi_bcast(grid%file_out,len_fn_out,irank_master,global_communicator)

    print *, 'grid%ndim on node ', grid%id, grid%ndim, get_my_rank()

  end subroutine synchronize_grid_data

#endif

end module
