module exacorr_global

!This module controls storage and interaction with the global data used in ExaCorr

        use exacorr_datatypes
        use dirac_mo

        implicit none
        private

!       Global data (should limit this to the minimum required to function)
        integer :: my_MPI_master=0 ! Define the MPI master proces (is always zero in the DIRAC MPI environment)
        type(basis_set_info_t) :: ao_basis ! information about the ao basis 
        type(MO):: qmo_all ! Stores all MO-coefficients and energies obtained from the parent code 
        logical, public, parameter :: DEBUG = .TRUE. ! switch on and off debug

!       Functions and routines to interact with these data
        public ncartesian
        public get_nao
        public initialize_global_data
        public delete_global_data
        public get_gtos
        public get_mo_coefficients

       contains

         integer function get_nao()
             get_nao = ao_basis%nao
         end function get_nao

         subroutine initialize_global_data

#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: my_MPI_rank
#endif
             integer :: ierr

#if defined (VAR_MPI)
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
             if (my_MPI_rank /= my_MPI_master) then
!               Slaves go to the sync immediately and get all that is necessary from master
                call synchronize_global_data
                return
             end if 
#endif
             ! read MO coefficients and ao basis set information from DFCOEF
             call read_from_dfcoef (ao_basis,qmo_all,ierr)

             if (ierr == 0) then
                call print_date ("retrieved basis set information")
             else
                print*, "Error code",ierr," while reading DFCOEF"
                call quit ('Could not read DFCOEF')
             end if
             
             call synchronize_global_data
 
         end subroutine initialize_global_data

         subroutine synchronize_global_data
 
#if defined (VAR_MPI)
             use interface_to_mpi
             use exacc_cfg

             integer            :: nao, nshells, n_primitives
             integer            :: ish
             integer            :: my_MPI_rank

!            The next line is to be commented out when running the stand-alone code
             call exacc_sync_cw()
             call print_date("Slave was here")
             if (exa_input%talsh) return
             call print_date("Slave passed the talsh check")
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)

             call sync_qmo(qmo_all,my_MPI_master)

             call interface_mpi_bcast (ao_basis%nshells,1,my_MPI_master,global_communicator)
             call interface_mpi_bcast (ao_basis%nao,1,my_MPI_master,global_communicator)
             nshells = ao_basis%nshells
             nao     = ao_basis%nao
   
             if (my_MPI_rank /= my_MPI_master) nullify (ao_basis%gtos)
             if (my_MPI_rank /= my_MPI_master) allocate(ao_basis%gtos(nshells))
             if (my_MPI_rank /= my_MPI_master) nullify (ao_basis%shell_indices)
             if (my_MPI_rank /= my_MPI_master) allocate(ao_basis%shell_indices(nao))
             
!            Note that this communication step could be speeded up by packing and unpacking the data for all shells
!            Assuming that the time for this step is insignificant, a simpler option is chosen here.
             do ish = 1, nshells
                call interface_mpi_bcast (ao_basis%gtos(ish)%orb_momentum,1,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%atom_number,1,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%n_primitives,1,my_MPI_master,global_communicator)
             end do

!            We now know the size of the shells, the exponents and coefficients are already allocated on master,
!            but not on the slaves, so do this first before starting the next communication loop
             if (my_MPI_rank /= my_MPI_master) then
                do ish = 1, nshells
                   n_primitives = ao_basis%gtos(ish)%n_primitives
                   nullify  (ao_basis%gtos(ish)%exponent)
                   allocate (ao_basis%gtos(ish)%exponent(n_primitives))
                   nullify  (ao_basis%gtos(ish)%coefficient)
                   allocate (ao_basis%gtos(ish)%coefficient(n_primitives))
                end do
             end if
 
             do ish = 1, nshells
                n_primitives = ao_basis%gtos(ish)%n_primitives
                call interface_mpi_bcast (ao_basis%gtos(ish)%coord,3,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%exponent,n_primitives,my_MPI_master,global_communicator)
                call interface_mpi_bcast (ao_basis%gtos(ish)%coefficient,n_primitives,my_MPI_master,global_communicator)
             end do

!            We finally broadcast the pointers from a basis function index to its shell that were determined on master
             call interface_mpi_bcast (ao_basis%shell_indices,nao,my_MPI_master,global_communicator)
#else
             return ! For running a serial version no action is required.
#endif
 
         end subroutine synchronize_global_data

         subroutine delete_global_data
         !TODO: check for memory leaks here. No so important as we usually stop right after and there is not much data anyhow
             deallocate(ao_basis%shell_indices)
             deallocate(ao_basis%gtos)
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
             get_shell_index = ao_basis%shell_indices(ifun)
         end function get_shell_index

         type(basis_func_info_t) function get_gto(ish)
             integer, intent(in) :: ish
             get_gto = ao_basis%gtos(ish)
         end  function get_gto

         subroutine get_gtos(first_ao,nao,gtos,nshells)

          integer, intent(in):: first_ao, nao
          type(basis_func_info_t), intent(inout), allocatable :: gtos(:)
          integer, intent(out) :: nshells
          integer :: ish, first_shell, last_shell
          
          first_shell = get_shell_index (first_ao)
          last_shell  = get_shell_index (first_ao+nao-1)
          nshells = last_shell - first_shell + 1

          allocate (gtos(nshells))
          do ish = 1, nshells
             gtos(ish) = get_gto(first_shell+ish-1)
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
