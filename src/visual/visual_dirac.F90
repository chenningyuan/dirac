module visual_dirac

  use visual_datatypes
  use visual_cfg
  use visual_utils
  use visual_mpi

! DIRAC modules
  use checkpoint

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
  public synchronize_dmat
#endif

  private

contains

  subroutine dirac_dmat_init(wf_order, cmo, dmat_ao, dmat_mo)
!   ----------------------------------------------------------------------------
    integer, intent(in) :: wf_order
    type(cmo_t), intent(in) :: cmo
    type(dmat_t), optional, intent(inout) :: dmat_ao
    type(dmat_t), optional, intent(inout) :: dmat_mo
    type(dmat_t) :: dmat_mo_from_occupation
    real(8), allocatable :: occupation(:)
    integer :: nao, nmo, nz
    integer :: i
!   ----------------------------------------------------------------------------

    if (wf_order == 0) then
!     initialize unperturbed density matrix from occupation numbers

      nao = sum(cmo%n_ao)
      nmo = sum(cmo%n_mo)
      nz  = cmo%nz

!     get orbital string
      allocate(occupation(nmo))
      call setup_orbital_string(cmo,occupation)

      dmat_mo_from_occupation%mat_dim = nmo
      dmat_mo_from_occupation%mat_nz  = nz
      call dmat_mo_from_occupation%alloc()
      do i = 1, nmo
        dmat_mo_from_occupation%elements((i-1)*nmo + i) = occupation(i)
      end do

      if (present(dmat_mo)) then
        dmat_mo = dmat_mo_from_occupation
      end if

      if (present(dmat_ao)) then
        dmat_ao%mat_dim = nao
        dmat_ao%mat_nz  = nz
        call dmat_ao%alloc()
        call dmat_transform_mo2ao(cmo, dmat_mo_from_occupation, dmat_ao)
        call dmat_ao%scale_elements(2.0d0)
      end if

      deallocate(occupation)

    else if (wf_order == 1) then
!     use PAMXVC
    end if
  end subroutine dirac_dmat_init


! -----------------------------------------------------------------------------
  subroutine setup_orbital_string(cmo,occupation)
!   "dcborb.h" has missing nish(ifsym) and nash(ifsym) information
#include "dcborb.h"
    type(cmo_t), intent(in) :: cmo
    real(8), intent(out):: occupation(:)
    integer:: i, ifsym, ioff

!   this routine will set up occupation vector
!   based on**VISUAL
!            .OCCUPATION
!   if this is not given, then it will occupy all
!   "occupied" electrons with 1.0
! -----------------------------------------------------------------------------

    occupation = 0.0d0

    ioff = 0
    do ifsym = 1, cmo%n_fsym
      ioff = ioff+cmo%n_po(ifsym)
      if (visual_cfg_use_orbital_string) then
        do i = 1, (cmo%n_mo(ifsym) - cmo%n_po(ifsym))
          occupation(ioff+i) = visual_cfg_occupation(i, ifsym)
        end do
      else
        if (nash(ifsym) > 0) then
          print *, 'ERROR: with open shells you should use .OCCUPATION'
          print *, '       otherwise occupation will be wrong'
          print *, '       so we better stop here'
          stop 1
        end if
        do i = 1, nish(ifsym)
          occupation(ioff+i) = 1.0d0
        end do
      end if
      ioff = cmo%n_mo(ifsym)
    end do

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

#ifdef VAR_MPI
    if (is_parallel()) call sync_visual_cfg()
#endif

  end subroutine visual_initialize_and_sync_nodes


! ------------------------------------------------------------------------------
  subroutine visual_initialize_data_from_dirac(mol, mo, print_level)
    type(mol_t), intent(out):: mol
    type(cmo_t), intent(out):: mo
    integer, optional :: print_level
    character*(:),allocatable:: group,section
    logical:: tobe_dirac
!
!   read data from DIRAC's CHECKPOINT file(s) into DTs
! ------------------------------------------------------------------------------

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

      ! todo - move this to a separate subroutine + decide what to save on checkpoint
      !call write_to_rs_checkpoint('/input/molecule', mol)
      !if (print_level > 0) &
      !  write(*, *) 'Molecule data written to DIRAC RS-CHECKPOINT'
    end if

!   gosia FIXME - the code is validated for scf only
    section = 'scf'
    call read_mo_from_dirac(section, mo)
      if (present(print_level)) then
        if (print_level > 0) then
          write(*, *) 'MO data read from DIRAC CHECKPOINT'
          call mo%display
        end if
      end if

  end subroutine visual_initialize_data_from_dirac


! ------------------------------------------------------------------------------
  subroutine read_mol_from_dirac(mol, group)
    type(mol_t), intent(out) :: mol
    character(len=*),intent(in) :: group
    real(8), allocatable :: temp(:)
! ------------------------------------------------------------------------------

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
    integer:: i
!
!   this subroutine is based on how things are done in dircmo.F90 and exacorr
!   gosia, summer 2023
! ------------------------------------------------------------------------------

    group = '/result/wavefunctions/'//trim(section)//'/'
    call checkpoint_query(group//'mobasis', tobe_dirac)
    if (.not.tobe_dirac) call visual_quit('No MO data available on CHECKPOINT. Stop.')

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
    ncoef = nz * sum(n_basis*n_mo)

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

    ! fill in index values
    allocate(index(nmo))
    do i = 1, nmo
       index(i) = i ! we set the index compatible with the fact that we store all mo's
    end do
    cmo%index    = index(npo+1:nmo)
    deallocate(index)

!   switch to unrestricted spinor picture
    !call convert_mo_to_complex(cmo)

  end subroutine read_mo_from_dirac


#ifdef VAR_MPI
! ------------------------------------------------------------------------------
  subroutine synchronize_data_from_dirac(mol, cmo)
    type(mol_t), intent(inout):: mol
    type(cmo_t), intent(inout):: cmo
    integer:: nmo,nao,nz,ncoeff
!
!   this subroutine syncs data stored in cmo and mol DTs, 
!   if any new variable has been added to these DTs, please also add it here
! ------------------------------------------------------------------------------

    call interface_mpi_bcast(mol%n_atoms,1,irank_master,global_communicator)

    if (.not.mpi_is_master()) then
      allocate(mol%symbols(mol%n_atoms))
      allocate(mol%geometry(3, mol%n_atoms))
      allocate(mol%nuc_charge(mol%n_atoms))
    end if

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

  end subroutine synchronize_data_from_dirac


! ------------------------------------------------------------------------------
  subroutine synchronize_dmat(dmat)
    type(dmat_t), intent(inout):: dmat
    integer:: n
! ------------------------------------------------------------------------------

    call interface_mpi_bcast(dmat%mat_dim,1,irank_master,global_communicator)
    call interface_mpi_bcast(dmat%mat_nz,1,irank_master,global_communicator)
    n = dmat%mat_dim*dmat%mat_dim*dmat%mat_nz

    if (.not.mpi_is_master()) then
      call dmat%alloc()
!gosia fixme
      !call dmat%alloc_q()
    end if
    call interface_mpi_bcast(dmat%elements,n,irank_master,global_communicator)
    !call interface_mpi_bcast(dmat%elements_q,n,irank_master,global_communicator)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine synchronize_grid_data(grid)
    type(grid_t), intent(inout):: grid
    integer:: len_gfn_inp, len_gfn_out
! ------------------------------------------------------------------------------

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
      if (len_gfn_inp > 0) allocate(character(len_gfn_inp):: grid%file_inp)
      if (len_gfn_out > 0) allocate(character(len_gfn_out):: grid%file_out)
    end if

    call interface_mpi_bcast(grid%file_inp,len_gfn_inp,irank_master,global_communicator)
    call interface_mpi_bcast(grid%file_out,len_gfn_out,irank_master,global_communicator)

  end subroutine synchronize_grid_data


! ------------------------------------------------------------------------------
  subroutine synchronize_grid_function_data(grid_function)
    type(grid_function_t), intent(inout):: grid_function
    integer:: len_gffn_inp, len_gffn_out, len_gf_name
! ------------------------------------------------------------------------------

    if (mpi_is_master()) then
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

! ------------------------------------------------------------------------------
  subroutine dmat_transform_mo2ao(cmo, dmat_mo, dmat_ao)

    type(cmo_t), intent(in):: cmo
    type(dmat_t), intent(inout):: dmat_ao
    type(dmat_t), intent(inout):: dmat_mo
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

end module
