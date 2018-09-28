module exacorr_global

!This module controls storage and interaction with the global data used in ExaCorr

        use exacorr_datatypes

        implicit none
        private

!       Global data (should limit this to the minimum required to function)
        integer :: my_MPI_master=0 ! Define the MPI master proces (is always zero in the DIRAC MPI environment)
        integer :: number_of_aos=0 ! This is set by master and communicated in the initialization step. 
        integer, allocatable :: shell_indices(:)  ! for each basis function the index of the shell to which it belongs
        type(basis_func_info_t), allocatable :: gtos(:)  ! list of all shells
        logical :: TALSH_ONLY = .TRUE. ! use single-node TALSH library instead of the full EXATENSOR library

!       Functions and routines to interact with these data
        public get_nao
        public get_talsh_only
        public set_talsh_only
        public initialize_global_data
        public delete_global_data
        public ncartesian
        public get_gto
        public fetch_gtos
        public get_shell_index

       contains
         
         subroutine set_talsh_only(new_value)
             logical, intent(in) :: new_value
             TALSH_ONLY = new_value
         end subroutine set_talsh_only

         logical function get_talsh_only()
             get_talsh_only = TALSH_ONLY
         end function get_talsh_only

         subroutine set_nao(new_value)
             integer, intent(in) :: new_value
             number_of_aos = new_value
         end subroutine set_nao

         integer function get_nao()
             get_nao = number_of_aos
         end function get_nao

         subroutine initialize_global_data

#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: my_MPI_rank
#endif
             integer :: number_aos, number_aoshells, nao, nshells
             integer :: ish, iao, ibas, ierr

#if defined (VAR_MPI)
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
             if (my_MPI_rank /= my_MPI_master) then
!               Not all information in common blocks is initialized for the slaves, directly go to the sync and get all that is necessary from master
                call synchronize_global_data
                return
             end if 
#endif
             ! store basis set definitions
             nshells = number_aoshells()  ! basis functions are organized in shells, this is their number
             allocate(gtos(nshells))
             call initialize_shells(gtos,ierr)
             if (ierr == 1) call quit ('This version can not handle contracted functions')
             nao     = number_aos(1)      ! retrieve number of aos from the host code (e.g. DIRAC common block)
             call set_nao(nao)            ! store in global data
             print*, " Retrieved basis set information from DIRAC common blocks"
             
             ! make array with pointers from a basis function index to its shell
             allocate(shell_indices(nao))
             iao = 0
             do ish = 1, nshells
                do ibas = 1, ncartesian(gtos(ish)%orb_momentum)
                   iao = iao + 1
                   if (iao.gt.nao) call quit ('Error in initialize_global_data routine')
                   shell_indices(iao) = ish
                end do
             end do

             if (iao.ne.nao) call quit ('Error in initialize_global_data routine')
             call synchronize_global_data
 
         end subroutine initialize_global_data

         subroutine synchronize_global_data
 
#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: number_aos, number_aoshells, nao, nshells
             integer :: ish
             integer :: my_MPI_rank

             if (TALSH_ONLY) return
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
             if (my_MPI_rank == my_MPI_master) nao = get_nao()
             if (my_MPI_rank == my_MPI_master) nshells = number_aoshells() 

             call interface_mpi_bcast (nshells,1,my_MPI_master,global_communicator)
             call interface_mpi_bcast (nao,1,my_MPI_master,global_communicator)
   
             if (my_MPI_rank /= my_MPI_master) call set_nao(nao)

             if (my_MPI_rank /= my_MPI_master) allocate(gtos(nshells))
             if (my_MPI_rank /= my_MPI_master) allocate(shell_indices(nao))
             
!            Note that this communication step could be speeded up by packing and unpacking the data for all shells
!            Assuming that the time for this step is insignificant, a simpler option is chosen here.
             do ish = 1, nshells
                call interface_mpi_bcast (gtos(ish)%orb_momentum,1,my_MPI_master,global_communicator)
                call interface_mpi_bcast (gtos(ish)%atom_number,1,my_MPI_master,global_communicator)
                call interface_mpi_bcast (gtos(ish)%n_primitives,1,my_MPI_master,global_communicator)
             end do

!            We now know the size of the shells, the exponents and coefficients are already allocated on master,
!            but not on the slaves, so do this first before starting the next communication loop
             if (my_MPI_rank /= my_MPI_master) then
                do ish = 1, nshells
                   nullify  (gtos(ish)%exponent)
                   allocate (gtos(ish)%exponent(gtos(ish)%n_primitives))
                   nullify  (gtos(ish)%coefficient)
                   allocate (gtos(ish)%coefficient(gtos(ish)%n_primitives))
                end do
             end if
 
             do ish = 1, nshells
                call interface_mpi_bcast (gtos(ish)%coord,3,my_MPI_master,global_communicator)
                call interface_mpi_bcast (gtos(ish)%exponent,gtos(ish)%n_primitives,my_MPI_master,global_communicator)
                call interface_mpi_bcast (gtos(ish)%coefficient,gtos(ish)%n_primitives,my_MPI_master,global_communicator)
             end do

!            We finally broadcast the pointers from a basis function index to its shell that were determined on master
             call interface_mpi_bcast (shell_indices,nao,my_MPI_master,global_communicator)
#else
             return ! For running a serial version no action is required.
#endif
 
         end subroutine synchronize_global_data

         subroutine delete_global_data
             deallocate(shell_indices)
             deallocate(gtos)
         end subroutine delete_global_data

         integer function ncartesian(l_value)
         ! Computes number of Cartesian functions for a given l-value
             integer, intent(in) :: l_value
             integer :: i
             ncartesian = 1
             do i = 2, l_value
                ncartesian = ncartesian + i
             end do
         end function ncartesian

         integer function get_shell_index(ifun)
             integer, intent(in) :: ifun
             get_shell_index = shell_indices(ifun)
         end function get_shell_index

         type(basis_func_info_t) function get_gto(ish)
             integer, intent(in) :: ish
             get_gto = gtos(ish)
         end  function get_gto

         subroutine fetch_gtos(first_ao,nao,gto,nshells)

          integer, intent(in):: first_ao, nao
          type(basis_func_info_t), intent(inout), allocatable :: gto(:)
          integer, intent(out) :: nshells
          integer :: ish, first_shell, last_shell

          first_shell = get_shell_index (first_ao)
          last_shell  = get_shell_index (first_ao+nao-1)
          nshells = last_shell - first_shell + 1

          allocate (gto(nshells))
          do ish = 1, nshells
             gto(ish) = get_gto(first_shell+ish-1)
          end do

         end subroutine fetch_gtos
!-------------------------------------------------------------------------------------------
        subroutine generate_ao_segments(num_large_ao,num_small_ao,desired_seg_size,ao_info,&
                                       &num_large_segs,num_small_segs,segs_large,segs_small,ierr)
!This subroutine partitions the range of large/small basis functions into segments
!such that the size of each segment is close to the <desired_seg_size> as much as possible.
!Different segments may have different sizes in order to align the segments to atomic shell boundaries.
         implicit none
         integer, intent(in):: num_large_ao                !in: total number of large component basis functions
         integer, intent(in):: num_small_ao                !in: total number of small component basis functions
         integer, intent(in):: desired_seg_size(0:2)       !in: desired size of a segment (0: target value, 1: minimal value, 2: maximum value)
         type(basis_func_info_t), intent(in):: ao_info(1:) !in: information on each basis function
         integer, intent(out):: num_large_segs             !out: number of large component segments
         integer, intent(out):: num_small_segs             !out: number of small component segments
         integer, intent(out), allocatable:: segs_large(:) !out: size of each segment in the large component range
         integer, intent(out), allocatable:: segs_small(:) !out: size of each segment in the small component range
         integer, intent(out):: ierr                       !out: error code (0:success)

!        place holder initializations for testing compilation and cmake inclusion of of this module
         num_large_segs = 1 
         allocate (segs_large(num_large_segs))
         segs_large(1) = num_large_ao
         num_small_segs = 0 
         ierr=0
         return
        end subroutine generate_ao_segments
!----------------------------------------------------------------------------------
        subroutine generate_mo_segments(num_occ,num_virt,desired_seg_size,mo_info,&
                                       &num_occ_segs,num_virt_segs,segs_occ,segs_virt,ierr)
!This subroutine partitions the range of occupied/virtual molecular spinors into segments
!such that the size of each segement is close to the <desired_seg_size> as much as possible.
!Different segements may have different sizes in order to align the segements by symmetry IR.
         implicit none

         integer, intent(in):: num_occ                    !in: total number of occupied spinors
         integer, intent(in):: num_virt                   !in: total number of virtual spinors
         integer, intent(in):: desired_seg_size(0:2)      !in: desired size of a segment (0: target value, 1: minimal value, 2: maximum value)
         type(mo_spinor_info_t), intent(in):: mo_info(1:) !in: information on each molecular spinor
         integer, intent(out):: num_occ_segs              !out: number of segments in the occupied range
         integer, intent(out):: num_virt_segs             !out: number of segments in the virtual range
         integer, intent(out), allocatable:: segs_occ(:)  !out: size of each segment in the occupied range
         integer, intent(out), allocatable:: segs_virt(:) !out: size of each segment in the virtual range
         integer, intent(out):: ierr                      !out: error code (0:success)

!        place holder initializations for testing compilation and cmake inclusion of of this module
         num_occ_segs  = 1 
         allocate (segs_occ(num_occ_segs))
         segs_occ(1) = num_occ
         num_virt_segs = 1 
         allocate (segs_virt(num_virt_segs))
         segs_virt(1) = num_virt
         ierr=0
         return
        end subroutine generate_mo_segments

end module exacorr_global
