module exacc_cfg

  use exacorr_datatypes, only : exacc_input

  implicit none

#ifdef VAR_MPI
  public exacc_sync_cw
#endif

  integer, parameter :: LENGTH=100
  character(LENGTH), public :: string_occupied, string_virtual
  type(exacc_input), public :: exa_input

#ifdef VAR_MPI

contains
  subroutine exacc_sync_cw()
!
! please add bcast here if a new variable has been added to the list above
!
  use interface_to_mpi
 
  call interface_mpi_bcast_l0(exa_input%talsh,1,0,global_communicator)
  call interface_mpi_bcast_l0(exa_input%ccd,1,0,global_communicator)
  call interface_mpi_bcast(exa_input%nkr_occ,1,0,global_communicator)
  call interface_mpi_bcast(exa_input%nkr_vir,1,0,global_communicator)
  call interface_mpi_bcast(exa_input%mokr_occ,exa_input%nkr_occ,0,global_communicator)
  call interface_mpi_bcast(exa_input%mokr_vir,exa_input%nkr_vir,0,global_communicator)
  call interface_mpi_bcast(exa_input%t_econv,1,0,global_communicator)

  end subroutine exacc_sync_cw
#endif

end module
