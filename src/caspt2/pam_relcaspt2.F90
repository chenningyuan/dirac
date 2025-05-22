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

   call print_title("CASPT2")

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
   use module_global_variables, only: tmem
   use module_index_utils, only: set_global_index
   use module_time

   tmem = 0.0d+00
   call sync_input_variables
   call set_global_index

   call write_allocated_memory_size
   call get_current_time(init_time); call print_time(init_time); start_time = init_time

   call dcaspt2_run_subprograms(.false.)

   call pam_caspt2_mpi_end
contains
   subroutine pam_caspt2_mpi_end
#ifdef HAVE_MPI
      use module_global_variables, only: nprocs
      implicit none
#include "dirac_partask.h"
      if (nprocs > 1) call DIRAC_PARCTL(RELEASE_NODES)
#endif
   end subroutine pam_caspt2_mpi_end
end subroutine CASPT2_MAIN

subroutine RELIVO_MASTER
   use module_global_variables, only: rank, nprocs
#ifdef HAVE_MPI
   use interface_to_mpi
#endif
   implicit none
#ifdef HAVE_MPI
   call interface_mpi_comm_rank(global_communicator, rank)
   call interface_mpi_comm_size(global_communicator, nprocs)
   if (nprocs > 1) call DIRAC_PARCTL(23)
#else
   rank = 0; nprocs = 1
#endif

   call print_title("IVO")

   call IVO_MAIN
end subroutine RELIVO_MASTER

subroutine RELIVO_NODE
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

   call IVO_MAIN
end subroutine RELIVO_NODE

subroutine IVO_MAIN
   use module_index_utils, only: set_global_index
   use module_time
   use module_global_variables
   implicit none
   integer :: ninact_orig, nact_orig, nsec_orig, nelec_orig

   tmem = 0.0d+00
   call sync_input_variables
   ! define inact,act,sec,elec values for IVO calculation since input values are for CASPT2
   ninact_orig = ninact; nact_orig = nact; nsec_orig = nsec; nelec_orig = nelec
   ! nact and nelec is the number of occupied orbitals, nsec is the number of virtual orbitals
   ninact = 0; nact = ninact_orig + nelec_orig; nelec = nact; nsec = nsec_orig + nact_orig - nelec_orig
   call set_global_index

   call dcaspt2_run_subprograms(.true.)

   if (rank == 0) then
      ! PCMOIN reads DFPCMO and write it to CHECKPOINT
      ! therefore we need to rename DFPCMONEW to DFPCMO in order to write the IVO result to CHECKPOINT
#if defined SYS_WINDOWS
      call system('move DFPCMO DFPCMO.orig')
      call system('move DFPCMONEW DFPCMO')
#else
      call system('mv DFPCMO DFPCMO.orig')
      call system('mv DFPCMONEW DFPCMO')
#endif
      call PCMOIN ! Write DFPCMO to CHECKPOINT
   end if

   ! restore original input values
   ninact = ninact_orig; nact = nact_orig; nsec = nsec_orig; nelec = nelec_orig

   call pam_ivo_mpi_end
contains
   subroutine pam_ivo_mpi_end
#ifdef HAVE_MPI
      use module_global_variables, only: nprocs
      implicit none
#include "dirac_partask.h"
      if (nprocs > 1) call DIRAC_PARCTL(RELEASE_NODES)
#endif
   end subroutine pam_ivo_mpi_end
end subroutine IVO_MAIN


subroutine sync_input_variables
#ifdef HAVE_MPI
   use module_mpi, only: bcast_wrapper
#endif
   use module_global_variables
   implicit none
#include "dcbpsi.h"
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
   ! run CASCI/CASPT2?
   call bcast_wrapper(docasci, 0)
   call bcast_wrapper(docaspt2, 0)
   ! IVO input variables
   call bcast_wrapper(doivo, 0)
   call bcast_wrapper(occ_mo_num, 0)
   call bcast_wrapper(vcut_mo_num, 0)
#endif
end subroutine sync_input_variables

subroutine print_title(program_name)
   use module_global_variables, only: rank
   implicit none
   character(len=*), intent(in) :: program_name
   character(:), allocatable :: title
   integer, parameter :: title_length = 118

   title = "RELATIVISTIC " // program_name // " PROGRAM"
   if (rank == 0) then
      call TITLER(title, "*", title_length)
      print '(A)', "Written by: Kohei Noda, Minori Abe, Sumika Iwamuro & Yasuto Masuda"
      print '(A)', "Reference: Y. Masuda, K. Noda, S. Iwamuro, M. Hada, N. Nakatani, M. Abe. &
                   J. Chem. Theory Comput. 2025, 21, 3, 1249–1258"
      print '(A)', ""
   end if
end subroutine print_title