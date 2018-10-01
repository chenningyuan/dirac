module exacorr_global

!This module controls storage and interaction with the global data used in ExaCorr

        use exacorr_datatypes
        use dirac_mo

        implicit none
        private

!       Global data (should limit this to the minimum required to function)
        integer :: my_MPI_master=0 ! Define the MPI master proces (is always zero in the DIRAC MPI environment)
        integer :: number_of_aos=0 ! This is set by master and communicated in the initialization step. 
        integer, allocatable :: shell_indices(:)  ! for each basis function the index of the shell to which it belongs
        type(basis_func_info_t), allocatable :: gtos(:)  ! list of all shells
        type(MO):: qmo_all ! Stores all MO-coefficients and energies obtained from the parent code 
        logical :: TALSH_ONLY = .TRUE. ! use single-node TALSH library instead of the full EXATENSOR library

!       Functions and routines to interact with these data
        public ncartesian
        public get_nao
        public get_talsh_only
        public set_talsh_only
        public initialize_global_data
        public delete_global_data
        public get_gtos
        public get_mo_coefficients

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
             integer :: Number_AOShells, Number_AOs
             integer :: nao, nshells
             integer :: ish, iao, ibas, ierr

#if defined (VAR_MPI)
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
             if (my_MPI_rank /= my_MPI_master) then
!               Not all information in common blocks is initialized for the slaves, directly go to the sync and get all that is necessary from master
                call synchronize_global_data
                return
             end if 
#endif
             ! Get dimension information and allocate gtos. 
             ! TODO: figure out how to allocate inside read_from_dfcoef
             ! TODO: neither passing an allocatable array nor a pointer array seems to work
             nshells = Number_AOShells()
             nao = Number_AOs(1)
             call set_nao(nao)
             allocate(gtos(nshells))

             ! read MO coefficients and all other necessary data from DFCOEF
             call read_from_dfcoef (nao,nshells,gtos,qmo_all,ierr)

             if (ierr == 0) then
                call print_date ("retrieved basis set information from DIRAC common blocks")
             else
                print*, "Error code",ierr," while reading DFCOEF"
                call quit ('Could not read DFCOEF')
             end if
             
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
             call delete_qmo (qmo_all)
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

         subroutine get_gtos(first_ao,nao,gto,nshells)

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

         end subroutine get_gtos

     subroutine get_MO_coefficients (qmo,mo_list,nmo)
!***********************************************************************
!
!     Fill array with subset of MO coefficients and related information
!     Retrives this from the qmo_all array that should be synced
!
!     Input
!        nmo     : number of selected MOs
!        mo_list : indices of selected MOs
!
!     Output
!        QMO  : MOs (coefficients, energies, etc.)
!
!
!     Written by L.Visscher Sep 2018
!
!
!***********************************************************************

      use dirac_mo

      integer, intent(in)     :: nmo,mo_list(nmo)
      type(MO), intent(inout) :: qmo

      integer  :: i, nao, nmo_all, nz

      nao     = qmo_all%nao
      nmo_all = qmo_all%nmo
      nz      = qmo_all%nz 

      ! allocate arrays in qmo to hold the active MO's specified in mo_list
      call create_qmo (qmo,nao,nmo,nz)

      ! copy the requested mo's
      do i = 1, nmo
         qmo%index(i)   = mo_list(i)
         qmo%coeff(:,i,:) = qmo_all%coeff(:,mo_list(i),:)
         qmo%energy(i) = qmo_all%energy(mo_list(i))
         qmo%boson_irrep(i) = qmo_all%boson_irrep(mo_list(i))
      enddo

      end

end module exacorr_global
