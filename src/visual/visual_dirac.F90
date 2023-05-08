module visual_dirac

  use checkpoint
  use visual_datatypes, only: mol_t, ao_t

  implicit none
  public read_from_dirac


  private

  interface read_from_dirac
    module procedure read_from_dirac
  end interface

contains

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
      if (print_level > 0) call mol%display

      !call write_to_rs_checkpoint('/input/molecule', mol)
      if (print_level > 0) &
        write(*, *) 'Molecule data written to DIRAC RS-CHECKPOINT'
    end if

    group = '/input/aobasis/'
    call checkpoint_query(group, tobe_dirac)
    if (tobe_dirac) then
      call read_ao_from_dirac(ao, group)
      write(*, *) 'AO data read from DIRAC CHECKPOINT'
    end if


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

end module
