module exacc_cfg

  use exacorr_datatypes, only : exacc_input

  implicit none

#ifdef VAR_MPI
  public exacc_sync_cw
#endif

  integer, parameter :: LENGTH=100
  character(LENGTH), public :: string_occupied, string_virtual, string_occ_beta, string_vir_beta
  type(exacc_input), public :: exa_input

#ifdef VAR_MPI

contains
  subroutine exacc_sync_cw(exa_MPI_rank)
!
! please add bcast here if a new variable has been added to the list above
!
  use interface_to_mpi

  integer :: exa_MPI_rank

  call interface_mpi_bcast_l0(exa_input%ccd,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%cc2,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%cis,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%lambda,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%noints,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%beta_occ,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%beta_vir,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast_l0(exa_input%do_triples,1,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%exa_blocksize,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%talsh_buff,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%moint_scheme,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%ncycles,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%print_level,1,exa_MPI_RANK,global_communicator)
  
  call interface_mpi_bcast(exa_input%mokr_occ,exa_input%nkr_occ,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mokr_vir,exa_input%nkr_vir,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mo_occ,exa_input%nocc,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%mo_vir,exa_input%nvir,exa_MPI_RANK,global_communicator)

  call interface_mpi_bcast(exa_input%t_econv,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%level_shift,1,exa_MPI_RANK,global_communicator)
  call interface_mpi_bcast(exa_input%t_cholesky,1,exa_MPI_RANK,global_communicator)
  
  end subroutine exacc_sync_cw
#endif

end module
