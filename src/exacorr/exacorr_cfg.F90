module exacc_cfg

  implicit none

#ifdef VAR_MPI
  public exacc_sync_cw
#endif

  integer, parameter        :: LENGTH=100
  character(LENGTH), public :: string_occupied, string_virtual
  logical, public           :: exacc_use_talsh, exacc_do_ccd
  real(8), public           :: econv

#ifdef VAR_MPI

contains
  subroutine exacc_sync_cw()
!
! please add bcast here if a new variable has been added to the list above
!
  use interface_to_mpi

  call interface_mpi_bcast(string_occupied,LENGTH,0,global_communicator)
  call interface_mpi_bcast(string_virtual,LENGTH,0,global_communicator)
  call interface_mpi_bcast_l0(exacc_use_talsh,1,0,global_communicator)
  call interface_mpi_bcast_l0(exacc_do_ccd,1,0,global_communicator)
  call interface_mpi_bcast(econv,1,0,global_communicator)

  end subroutine exacc_sync_cw
#endif

end module
