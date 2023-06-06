module visual_dirac

  use checkpoint
  use visual_datatypes, only: mol_t, ao_t, grid_t, grid_function_t
  use visual_cfg

  use visual_mpi
#ifdef VAR_MPI
   use interface_to_mpi
   use integer_kind_mpilib
#endif


  implicit none

  public visual_initialize_and_sync_nodes
  public visual_initialize_data_from_dirac
#ifdef VAR_MPI
  public synchronize_grid_data
  public synchronize_grid_function_data
  public synchronize_data_from_dirac
#endif


  private


contains

  subroutine visual_initialize_and_sync_nodes(hello)
    use dirac_cfg, only: dirac_cfg_visual
    use visual_mpi
#ifdef VAR_MPI
    use interface_to_mpi
    use integer_kind_mpilib
#endif
#include "dirac_partask.h"  
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

#ifdef VAR_MPI
    if (is_parallel()) call sync_visual_cfg()
#endif

  end subroutine visual_initialize_and_sync_nodes


  subroutine visual_initialize_data_from_dirac(mol, aobasis)

    type(mol_t), intent(out):: mol
    type(ao_t), intent(out):: aobasis
    character*(:),allocatable:: group
    logical:: tobe_dirac

    group = '/input/molecule/'
    call checkpoint_query(group, tobe_dirac)
    if (tobe_dirac) then
      call read_mol_from_dirac(mol, group)
      !write(*, *) 'Molecule data read from DIRAC CHECKPOINT'
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

    !write(*, *) 'in synchronize_data_from_dirac, mol on node: ', get_my_rank()
    !call mol%display()

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
