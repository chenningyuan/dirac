subroutine RELCASPT2_MASTER
   use module_global_variables, only: rank, nprocs
#ifdef HAVE_MPI
   use interface_to_mpi
#endif
   implicit none
#ifdef HAVE_MPI
   call interface_mpi_comm_rank(global_communicator, rank)
   call interface_mpi_comm_size(global_communicator, nprocs)
   if (nprocs > 1) call DIRAC_PARCTL(22)
#else
   rank = 0; nprocs = 1
#endif

   ! call print_head
   print '(2(a,1x,i0))', 'initialization of mpi, rank : ', rank, ' nprocs : ', nprocs
   print *, ''
   print *, ' START RELATIVISIC CASPT2 PROGRAM'
   print *, ''

   call CASPT2_MAIN
end subroutine RELCASPT2_MASTER

subroutine RELCASPT2_NODE
   use module_global_variables, only: rank, nprocs
   use module_time
#ifdef HAVE_MPI
   use interface_to_mpi
#endif
   implicit none
#ifdef HAVE_MPI
   call interface_mpi_comm_rank(global_communicator, rank)
   call interface_mpi_comm_size(global_communicator, nprocs)
#else
   rank = 0; nprocs = 1
#endif

   call CASPT2_MAIN
end subroutine RELCASPT2_NODE

subroutine CASPT2_MAIN
    use module_time
    use module_index_utils, only: set_global_index

   call sync_input_variables

   call write_allocated_memory_size
   call get_current_time(init_time); call print_time(init_time); start_time = init_time
   ! call set_global_index

   call dcaspt2_run_subprograms

   call pam_caspt2_mpi_end
contains
   subroutine sync_input_variables
#ifdef HAVE_MPI
      use module_mpi, only: bcast_wrapper
#endif
      use module_global_variables
      implicit none
#ifdef HAVE_MPI
      ! Essential variables
      call bcast_wrapper(ninact, 0)
      call bcast_wrapper(nact, 0)
      call bcast_wrapper(nsec, 0)
      call bcast_wrapper(nelec, 0)
      call bcast_wrapper(dirac_version, 0)
      ! Optional variables
      call bcast_wrapper(eshift, 0)
      call bcast_wrapper(debug, 0)
      call bcast_wrapper(mdcint_scheme, 0)
      call bcast_wrapper(ncore, 0)
      ! RAS specific
      call bcast_wrapper(ras1_max_hole, 0)
      call bcast_wrapper(min_hole_ras1, 0)
      call bcast_wrapper(ras3_max_elec, 0)
      call bcast_wrapper(ras1_size, 0)
      call bcast_wrapper(ras2_size, 0)
      call bcast_wrapper(ras3_size, 0)
      ! worker nodes rasx_list are not allocated, if rasx_list is not allocated, allocate it with rasx_size
      if (.not. allocated(ras1_list)) allocate (ras1_list(ras1_size))
      if (.not. allocated(ras2_list)) allocate (ras2_list(ras2_size))
      if (.not. allocated(ras3_list)) allocate (ras3_list(ras3_size))
      call bcast_wrapper(ras1_list, 0)
      call bcast_wrapper(ras2_list, 0)
      call bcast_wrapper(ras3_list, 0)
#endif
      call set_global_index

   end subroutine sync_input_variables

   subroutine pam_caspt2_mpi_end
#ifdef HAVE_MPI
      use module_global_variables, only: nprocs
      implicit none
#include "dirac_partask.h"
      if (nprocs > 1) call DIRAC_PARCTL(RELEASE_NODES)
#endif
   end subroutine pam_caspt2_mpi_end
end subroutine CASPT2_MAIN
