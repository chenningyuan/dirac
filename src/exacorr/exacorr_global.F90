module exacorr_global

!This module controls storage and interaction with the global data used in ExaCorr

        use exacorr_datatypes
        use exacorr_mo

        implicit none
        private

!       Global data (should limit this to the minimum required to function)
        integer :: my_MPI_master=0 ! Define the MPI master proces (is always zero in the DIRAC MPI environment)
        type(basis_set_info_t) :: ao_basis ! information about the ao basis 
        type(cmo):: cspinor_all ! Stores all MO-coefficients and energies obtained from the parent code 
        type(one_el_t)  :: int_1el ! stores 1 electron integrals
        logical :: TALSH_ONLY = .TRUE. ! use single-node TALSH library instead of the full EXATENSOR library
        logical, public, parameter :: debug_exa = .false. ! switch on and off debug
        logical :: one_el_exist=.FALSE. ! logical if MRCONEE is available or not
        integer :: m_cholesky=-666  ! length of cholesky vector

!       Functions and routines to interact with these data
        public nfunctions
        public get_nao
        public get_basis_angular
        public get_nsolutions
        public get_talsh_only
        public set_talsh_only
        public initialize_global_data
        public delete_global_data
        public get_gtos
        public get_mo_coefficients
        public exist_one_el
        public get_one_el
        public set_m_cholesky
        public get_m_cholesky

       contains
         
         subroutine set_talsh_only(new_value)
             logical, intent(in) :: new_value
             TALSH_ONLY = new_value
         end subroutine set_talsh_only

         logical function get_talsh_only()
             get_talsh_only = TALSH_ONLY
         end function get_talsh_only

         integer function get_nao()
             get_nao = ao_basis%nao
         end function get_nao

         integer function get_basis_angular()
             get_basis_angular = ao_basis%basis_angular
         end function get_basis_angular

         integer function get_nsolutions()
             get_nsolutions = cspinor_all%nmo
         end function get_nsolutions

         subroutine initialize_global_data

             use exacorr_respect
             use exacc_cfg
             use exacorr_utils, only : print_date

#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: my_MPI_rank
#endif
             integer :: ierr
             logical :: tobe

#if defined (VAR_MPI)
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)
             if (my_MPI_rank /= my_MPI_master) then
!               Slaves go to the sync immediately and get all that is necessary from master
                call synchronize_global_data
                return
             end if 
#endif
             ! read MO coefficients and ao basis set information
             inquire (file='DFCOEF',exist=tobe)
             if (tobe) then
                call read_from_dirac (ao_basis,cspinor_all,'DFCOEF',ierr)
                if (exa_input%beta_occ .and. (.not.TALSH_ONLY)) then
                  if (exa_input%level_shift.lt.1.0D-10) then
                    call quit ('Open shell computations with DIRAC and EXACORR require level shift (LSHIFT)')
                  end if
               end if
             else
                inquire (file='RSD_MOS',exist=tobe)
                if (tobe) then
                   call read_from_respect (ao_basis,cspinor_all,ierr)
                end if
             end if

             if (ierr == 0) then
                call print_date ("retrieved basis set information")
             else
                print*, "Error code",ierr," while reading restart file"
                call quit ('Could not read restart file')
             end if

             ! read one electron integrals
             inquire (file='MRCONEE',exist=tobe)
             if (tobe) then
                call read_from_mrconee (int_1el, ierr)
                if (ierr == 0) then
                  call print_date("read one el. integrals from MRCONEE")
                  one_el_exist = .TRUE.
                else
                  print*, "Error code",ierr," while reading MRCONEE"
                  call quit ('Could not read MRCONEE')
                end if
             else
               one_el_exist = .FALSE.
               if (exa_input%beta_occ) then
                 call quit ('One electron integrals required in the open shell case (MOLTRA)')
               end if
             end if
             
             call synchronize_global_data

         end subroutine initialize_global_data

         subroutine synchronize_global_data

#if defined (VAR_MPI)
             use interface_to_mpi
             integer :: nao, nshells, n_primitives, nsp
             integer :: ish
             integer :: my_MPI_rank

!            The next line is to be commented out when running the stand-alone code
             if (TALSH_ONLY) return
             call interface_mpi_comm_rank (global_communicator,my_MPI_rank)

             call sync_mo(cspinor_all,my_MPI_master)

             call interface_mpi_bcast (ao_basis%nshells,1,my_MPI_master,global_communicator)
             call interface_mpi_bcast (ao_basis%nao,1,my_MPI_master,global_communicator)
             call interface_mpi_bcast (ao_basis%basis_angular,1,my_MPI_master,global_communicator)
             nshells       = ao_basis%nshells
             nao           = ao_basis%nao
   
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
             
!            One electron integrals
             call interface_mpi_bcast_l0 (one_el_exist,1,my_MPI_master,global_communicator)
             if(one_el_exist) then
               call interface_mpi_bcast (int_1el%n_spinor,1,my_MPI_master,global_communicator)
               nsp=int_1el%n_spinor
               call interface_mpi_bcast (int_1el%e_core,1,my_MPI_master,global_communicator)

               if (my_MPI_rank /= my_MPI_master) nullify (int_1el%h_core)
               if (my_MPI_rank /= my_MPI_master) allocate(int_1el%h_core(nsp,nsp))

               call interface_mpi_bcast (int_1el%h_core,2*nsp*nsp,my_MPI_master,global_communicator)
             end if
#else
             return ! For running a serial version no action is required.
#endif
 
         end subroutine synchronize_global_data

         subroutine delete_global_data
         !TODO: check for memory leaks here. No so important as we usually stop right after and there is not much data anyhow
             deallocate(ao_basis%shell_indices)
             deallocate(ao_basis%gtos)
             call dealloc_mo (cspinor_all)
             if(one_el_exist) then
               deallocate(int_1el%h_core)
             end if
         end subroutine delete_global_data

         integer function nfunctions(l_value,switch)
         ! Computes number of Cartesian/spherical functions for a given l-value
             integer, intent(in) :: l_value
             integer, intent(in) :: switch
             if (switch == 1) then
!               Cartesian
                nfunctions = l_value * ( l_value + 1 ) / 2
             else if (switch == 2) then
!               spherical
                nfunctions = 2 * (l_value - 1) + 1
             end if
         end function nfunctions

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

     subroutine get_MO_coefficients (cspinor,mo_list,nmo)
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

      use exacorr_mo

      integer, intent(in)     :: nmo,mo_list(nmo)
      type(cmo), intent(inout) :: cspinor

      integer  :: i, nao, nmo_all

      nao     = cspinor_all%nao
      nmo_all = cspinor_all%nmo !WHY IS THIS HERE?

      ! allocate arrays in cmo to hold the active MO's specified in mo_list
      call alloc_mo (cspinor,nao,nmo)

      ! copy the requested mo's
      do i = 1, nmo
         cspinor%index(i)       = mo_list(i)
         cspinor%coeff(:,i,:)   = cspinor_all%coeff(:,mo_list(i),:)
         cspinor%energy(i)      = cspinor_all%energy(mo_list(i))
         cspinor%boson_irrep(i) = cspinor_all%boson_irrep(mo_list(i))
      enddo

      end subroutine get_MO_coefficients

      logical function exist_one_el(e_core)

        real(8), intent(out) :: e_core

        exist_one_el=one_el_exist
        if (one_el_exist) then
           e_core=int_1el%e_core
        else
           e_core = 0.D0
        end if

      end function

      subroutine get_one_el(integrals,mo1_list,nmo1,mo2_list,nmo2,min_occ)

        complex(8), intent(out) :: integrals(:,:) !one lectron integrals
        integer, intent(in)     :: nmo1          ! the length of the mo basis
        integer, intent(in)     :: mo1_list(:)   ! and their indices
        integer, intent(in)     :: nmo2          ! the length of the mo basis
        integer, intent(in)     :: mo2_list(:)   ! and their indices
        integer, intent(in)     :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

        integer :: pmo, qmo

        if((mo1_list(nmo1)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MRCONEE (MOLTRA)')
        if((mo2_list(nmo2)-min_occ+1).gt.int_1el%n_spinor) call quit ('Not enough orbitals in MRCONEE (MOLTRA)')

        do qmo = 1, nmo2
          do pmo = 1, nmo1
            integrals(pmo, qmo) = int_1el%h_core(mo1_list(pmo)-min_occ+1,mo2_list(qmo)-min_occ+1)
          end do
        end do

     end subroutine get_one_el

     subroutine get_m_cholesky(m)

      integer, intent(out)  ::  m

      m=m_cholesky

     end subroutine get_m_cholesky

     subroutine set_m_cholesky(m)
       use interface_to_mpi
       
       integer, intent(in)  :: m
       integer :: my_MPI_rank

       call interface_mpi_comm_rank (global_communicator,my_MPI_rank)

       m_cholesky=m

       call interface_mpi_bcast (m_cholesky,1,my_MPI_rank,global_communicator)

    end subroutine set_m_cholesky


end module exacorr_global
