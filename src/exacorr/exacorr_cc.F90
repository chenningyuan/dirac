module exacorr_cc

!This module contains the routines needed to calculate a coupled cluster wave function
!almost everything requires MPI, but we provide a stub function to make compilation of 
!non-MPI versions of DIRAC possible.

#if defined (VAR_MPI)

        use exatensor
        use exacorr_datatypes
        use exacorr_utils
        use intermediates

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0),ONE_QUARTER=(0.25D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0)!, THREE=(3.0,0.0), SIX=(6.0,0.0)
        type(tens_rcrsv_t),public :: one_tensor
        real(8) :: e_core
        logical :: one_el_exist

        private

        public exacorr_cc_driver

        interface diis
           module procedure diis
        end interface

       contains

!-------------------------------------------------------------------------------------------
        subroutine exacorr_cc_driver (exa_input)

!        This subroutine drives the ccsd energy calculation
!        The implementation is based on the parallel ExaTensor library

!        Written by Lucas Visscher, winter 2016/2017

         use interface_to_mpi
         use exacorr_global
         use exacorr_tensor_methods
         use exacorr_ao_to_mo
         use exacorr_gradient
         

         implicit none

         type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir           ! the size of the mo basis for occupied and virtual spinors

         integer(INTD) :: my_role
         integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode

         integer(INTD)                  :: ierr
         complex(8)                     :: result_val

!        Common tensors
         type(exatns_comm_tens_t) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_intg_tens_t) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t) :: l1_tensor, l2_tensor
!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t)         :: result_tensor
!        space ids and roots
         type(space_dims) :: dims
!        CCSD control and result variables
         real(8)            :: mp2_energy, ccsd_energy !, t_energy(3) ! triples not implemented
         real(8)            :: t_target_precision

         ! classes for the definition of spaces (AO and 2 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         class(h_space_t), pointer                :: occ_space, vir_space

         ! variables for the definition of spaces (AO and 2 different MO spaces)
         integer(INTL)              :: ao_space_root, occ_space_root, vir_space_root
         integer(INTL),dimension(2) :: oo_root, ov_root, vo_root, vv_root
         integer(INTL),dimension(4) :: oooo_root, ooov_root, oovv_root, vovo_root, &
                                       vvoo_root, vovv_root, vvvv_root
         integer(INTD),dimension(2) :: oo_id, ov_id, vo_id, vv_id 
         integer(INTD),dimension(4) :: oooo_id, ooov_id, oovv_id, vovo_id, &
                                       vvoo_id, vovv_id, vvvv_id
         integer(INTD)              :: ao_space_id, occ_space_id, vir_space_id

         integer :: ncycles

         !CCD
         logical :: CCD, lambda, noints

         !controlling print
         integer :: print_level

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) call print_date(" Starting CC calculation with exatensor")

!        Initialize and register AO and MO spaces and get their IDs for further reference
         call exacorr_cc_spaces (exa_input,ao_space,occ_space,vir_space,ao_space_id,occ_space_id,vir_space_id)
         if (my_MPI_rank == 0) call print_date(" Registered AO and MO spaces")

!        Register and initialize user defined methods needed in the CC calculation
         call exacorr_cc_methods (exa_input,ao_space_id,occ_space_id,vir_space_id)
         if (my_MPI_rank == 0) call print_date(" Initialized and registered tensor methods")

!        copy input variables to local ones for ease of reference
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         t_target_precision = exa_input%t_econv
         CCD = exa_input%ccd
         lambda = exa_input%lambda
         noints = exa_input%noints
         ncycles = exa_input%ncycles
         print_level = exa_input%print_level

         if (my_MPI_rank == 0) then
            write(*,*) ""
            write(*,'(A)') " Characteristics of this calculation"
            write(*,'(A)') " -----------------------------------"
            write(*,'(A,L1)') " - Solve CCD equations :              ", exa_input%ccd
            write(*,'(A,L1)') " - Solve Lambda equations :           ", exa_input%lambda
            write(*,'(A,I3)') " - Number of atomic orbitals :      ",ao_space%get_root_id(ierr)
            write(*,'(A,I3)') " - Number of occupied spinors :     ", exa_input%nocc
            write(*,'(A,I3)') " - Number of virtual spinors :      ", exa_input%nvir
            write(*,'(A,ES10.1)') " - Target precision of amplitudes :", exa_input%t_econv
            write(*,'(A,I3)') " - Maximum number of cycles :       ", exa_input%ncycles
            write(*,'(A,I3)') " - ExaTensor block size :           ",exa_input%exa_blocksize
            write(*,'(A,I3)') " - Integral transformation scheme : ",exa_input%moint_scheme
            if (exa_input%level_shift.gt.1.0D-14) then
              write(*,'(A,ES10.1)') " - Level Shift :                    ",exa_input%level_shift
            end if
            write(*,'(A,I3)') " - Print Level :                    ",exa_input%print_level
            if(exa_input%lambda) then
              write(*,*) ""
              write(*,*) "Attention: For properties the density matrix is saved on a node."
              write(*,*) "           buffer size: ",exa_input%talsh_buff,"GB"
            end if
            write(*,*) ""
         end if

         !Create all id and root arrays, neccessary for exatns_tensor_create
         ao_space_root=ao_space%get_root_id(ierr)
         occ_space_root=occ_space%get_root_id(ierr)
         vir_space_root=vir_space%get_root_id(ierr)
         dims%occ_space_id=occ_space_id
         dims%vir_space_id=vir_space_id
         dims%occ_space_root=occ_space_root
         dims%vir_space_root=vir_space_root
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vo_id = (/vir_space_id, occ_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vo_root = (/vir_space_root, occ_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oooo_id = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
         ooov_id = (/occ_space_id, occ_space_id, occ_space_id, vir_space_id/)
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvoo_id = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do.
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(int(global_communicator,INTD))

         if(ierr.eq.EXA_SUCCESS) then
            ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
            if(my_role.eq.EXA_DRIVER) then

!           Create and initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(result_tensor,"result_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(one_tensor,"another_one_tensor",EXA_DATA_KIND_C8)
            ! Note:  also get_integral_tensors will define a one_tensor and call it "one_tensor", so we needed another name
            ierr=exatns_tensor_init(one_tensor,ONE)
            
!           Get fixed two-body tensors
            ierr=exatns_tensor_create(int_t%oooo,"oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%ooov,"ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%oovv,"oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovo,"vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovv,"vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vvvv,"vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
            call get_CC_integrals (ao_space_id,ao_space,int_t,exa_input%moint_scheme)
            call print_date('initialized two electron integrals')

            ! If we have the 1-body integrals we rebuild the Fock matrix and check consistency, otherwise we use the orbital energies
            ierr=exatns_tensor_create(comm_t%foo,"foo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fov,"fov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fvv,"fvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            call init_fock (comm_t,int_t,oo_id,oo_root,ov_id,ov_root,vv_id,vv_root, print_level)
            call print_date('initialized Fock matrix')

!           Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
            if (.not.CCD) then
               ierr=exatns_tensor_create(comm_t%t1,"t1",vo_id,vo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_init(comm_t%t1,ZERO)
            end if

!           Initialize t2 amplitudes with the MP2 values
            ierr=exatns_tensor_create(comm_t%t2,"t2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(comm_t%t2,ZERO)
!           the tensor_add routine is not implemented in exatensor, work around this by using a contraction
            ierr=exatns_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)*K()",comm_t%t2,int_t%oovv,one_tensor) ! o(2)v(2)
            ierr=exatns_tensor_transform(comm_t%t2,'Denominate')

            ierr=exatns_tensor_init(result_tensor,ZERO)
!           Calculate and print the MP2 energy
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv,ONE_QUARTER) ! o(2)v(2)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=real(result_val,8) ! Casting into real(8)

            print*, ""
            print*," MP2 energy = ", mp2_energy
            print*, ""

!           Start CC(S)D iterations

            ierr=exatns_tensor_create(comm_t%goo,"goo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%gvv,"gvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hoo,"hoo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hvv,"hvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            if (.not.CCD) ierr=exatns_tensor_create(comm_t%hov,"hov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%a_int,"aint",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%h_int,"hint",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%tau,"tau",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            if (.not.CCD) then
                call diis (comm_t%t2,vvoo_id,vvoo_root,comm_t%t1,vo_id,vo_root,mode=1)
            else
                call diis (comm_t%t2,vvoo_id,vvoo_root,mode=1)
            end if
            if (ierr.eq.EXA_SUCCESS) then
               print*," ccdriver: allocated all tensors"
            else
               print*," ccdriver: error in tensor allocation, code=",ierr
               return
            end if

            call solve_t_equations(comm_t,int_t,dims,CCD,t_target_precision,ccsd_energy,ncycles)

!           Calculate (T) energy correction, not implemented yet
!           call t3_energy (t_energy,eps_occ,eps_vir,t1_tensor,t2_tensor,fov_tensor,ooov_tensor,vovv_tensor,oovv_tensor)
!           write (*,'(A,F20.15)') "  4th order triples correction :    ",t_energy(1)
!           write (*,'(A,F20.15)') "  5th order triples (T) correction :",t_energy(2)
!           write (*,'(A,F20.15)') "  5th order triples -T  correction :",t_energy(3)
!           write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",ccsd_energy+t_energy(1)
!           write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(2)
!           write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(3)

            if (.not.CCD) then
                print*," Final CCSD energy = ",ccsd_energy
            else
                print*," Final CCD energy = ",ccsd_energy
            end if
!           print*," Final CCSD(T) energy = ",ccsd_energy+t_energy(1)+t_energy(2)  ! no (T) implementation yet

!           Destroy tensors and shut down library
            if (.not.CCD) then
                call diis (comm_t%t2,t1_tensor=comm_t%t1,mode=2)
            else
                call diis(comm_t%t2,mode=2)
            end if

            if (lambda ) then

!            construct lambda tensors + initialize with t-amplitudes
             if (.not.CCD) then
                ierr=exatns_tensor_create(l1_tensor,"l1",ov_id,ov_root,EXA_DATA_KIND_C8)
                ierr=exatns_tensor_init(l1_tensor,ZERO)
                ierr=exatns_tensor_contract("L(i,a)+=T+(a,i)*K()",l1_tensor,comm_t%t1,one_tensor)
             end if

             ierr=exatns_tensor_create(l2_tensor,"l2",oovv_id,oovv_root,EXA_DATA_KIND_C8)
             ierr=exatns_tensor_init(l2_tensor,ZERO)
             ierr=exatns_tensor_contract("L(i,j,a,b)+=T+(a,b,i,j)*K()",l2_tensor,comm_t%t2,one_tensor)

!            solve lambda-equations and caluclate density matrix in MO basis
             if (.not.CCD) then
               call solve_lambda_equations(comm_t,int_t,l2_tensor,dims,t_target_precision,ncycles,print_level,l1_tensor)
               call exacorr_dm_driver(exa_input,dims,comm_t%t2,l2_tensor,comm_t%t1,l1_tensor)
             else
               call solve_lambda_equations(comm_t,int_t,l2_tensor,dims,t_target_precision,ncycles,print_level)
               call exacorr_dm_driver(exa_input,dims,comm_t%t2,l2_tensor)
             end if

            end if

            ierr=exatns_tensor_destroy(int_t%oooo)
            ierr=exatns_tensor_destroy(int_t%ooov)
            ierr=exatns_tensor_destroy(int_t%oovv)
            ierr=exatns_tensor_destroy(int_t%vovo)
            ierr=exatns_tensor_destroy(int_t%vovv)
            ierr=exatns_tensor_destroy(int_t%vvvv)
            ierr=exatns_tensor_destroy(comm_t%foo)
            ierr=exatns_tensor_destroy(comm_t%fov)
            ierr=exatns_tensor_destroy(comm_t%fvv)
            ierr=exatns_tensor_destroy(comm_t%tau)
            ierr=exatns_tensor_destroy(comm_t%hoo)
            if (.not.CCD) ierr=exatns_tensor_destroy(comm_t%hov)
            ierr=exatns_tensor_destroy(comm_t%hvv)
            ierr=exatns_tensor_destroy(comm_t%goo)
            ierr=exatns_tensor_destroy(comm_t%gvv)
            ierr=exatns_tensor_destroy(comm_t%a_int)
            ierr=exatns_tensor_destroy(comm_t%h_int)
            ierr=exatns_tensor_destroy(result_tensor)
            ierr=exatns_tensor_destroy(one_tensor)
            if (.not.CCD) then
              ierr=exatns_tensor_destroy(comm_t%t1)
              ierr=exatns_tensor_destroy(l1_tensor)
            end if
            ierr=exatns_tensor_destroy(comm_t%t2)
            ierr=exatns_tensor_destroy(l2_tensor)

            !Stop ExaTENSOR runtime
            ierr=exatns_stop()

            end if

         else
            write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         if (my_MPI_rank == 0) then
            print*, ""
            print*, " Leaving cc_driver routine"
         end if
         return

        end subroutine exacorr_cc_driver

        subroutine exacorr_cc_spaces (exa_input,ao_space,occ_space,vir_space,ao_space_id,occ_space_id,vir_space_id)

!        This subroutine defines and registers all spaces for ExaTensor

!        Written by Lucas Visscher, January 2019

         use interface_to_mpi
         use exacorr_global
         use exacorr_ao_to_mo

         implicit none

         type(exacc_input), intent(in) :: exa_input
         ! classes for the definition of spaces (AO and 2 different MO spaces)
         class(h_space_t), pointer     :: ao_space, occ_space, vir_space
         ! ids of these spaces (are not stored in the space class itself ?)
         integer(INTD), intent(out)    :: ao_space_id, occ_space_id, vir_space_id

         integer                :: nocc,nvir       ! the size of the mo basis for occupied and virtual ORBITALS
         integer(INTL)          :: nocc_l, nvir_l  ! the size of the mo basis for occupied and virtual SPINORS

         integer                                :: nao      ! the number of atomic orbitals)
         integer                                :: nshells  ! the number of shells (orbitals with the same center and l-value)
         integer                                :: basis_angular    ! 1=cartesian, 2=spherical
         type(basis_func_info_t), allocatable   :: gto(:)   ! array with information about the shells

         ! variables for the definition of spaces (AO and 2 different MO spaces)
         type(subspace_basis_t)     :: basis_ao, basis_occ, basis_vir
         integer(INTL)              :: nao_l ! (same as nao, as long integer)
         integer(INTD)              :: branch_factor, tavp_mng_depth
         integer                    :: my_MPI_rank

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable    :: color(:)
         integer(INTL)                          :: labs

         ! Loop counters and error codes
         integer(INTD) :: ierr, lsh

         integer:: l, nl

         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)

         nao           = get_nao()   ! This is the total number of ao's
         nao_l         = int(nao,8)  ! we need this as a long int for exatensor
         basis_angular = get_basis_angular() ! needed for setting basis functions

!        Get the basis set information
         call get_gtos(1,nao,gto,nshells)

         !Create a basis for the ao vector space:
         ierr = 0
         call basis_ao%subspace_basis_ctor(nao_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         !Assign the colors for the ao space (this indicates to exatensor where it may break up into subspaces)
         allocate(color(nshells))
         labs = 0
         do lsh = 1, nshells !set basis functions
            nl   =  nfunctions(gto(lsh)%orb_momentum,basis_angular)
            call color(lsh)%color_symmetry_ctor(lsh) ! take one shell as a color
            do  l = 1, nl
                labs = labs + 1
                call basis_ao%set_basis_func(labs,BASIS_ABSTRACT,ierr,symm=color(lsh))
                if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
            end do
         enddo

         !Finalize the ao space basis
         call basis_ao%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

         !Get ExaTENSOR TAVP-MNG depth:
         ierr = exatns_virtual_depth(tavp_mng_depth)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_virtual_depth() failed!')
         if(debug_exa .and. my_MPI_rank == 0) write(*,*) 'ExaCorr: Inferred TAVP-MNG depth = ',tavp_mng_depth

         !Register the ao vector space:
         branch_factor = max(2, (int(get_num_segments(int(nao,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('AO_space',basis_ao,ao_space_id,ao_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')
         if(debug_exa .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') ' ExaCorr: Registered AO vector space of dimension   ',nao, &
         ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1)) !debug
         if(debug_exa .and. my_MPI_rank == 0) call ao_space%print_it()

         !Copy user input to local arrays for the definition of MO spaces (for standard CC we only distinguish occupied and virtual spinors)
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         !Get the MO basis for the occupied and virtual spaces
         ierr = 0
         nocc_l    = int(nocc,8)  ! NB: imo, jmo and nmo count orbitals, nocc_l counts spinors
         nvir_l    = int(nvir,8)  ! NB: imo, jmo and nmo count orbitals, nvir_l counts spinors
         call basis_occ%subspace_basis_ctor(nocc_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')
         call basis_vir%subspace_basis_ctor(nvir_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         !Set the functions inside the occupied basis (without colors, not needed for MOs)
         do labs = 1, nocc_l
            call basis_occ%set_basis_func(labs,BASIS_ABSTRACT,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
         enddo

         !Set the functions inside the virtual basis (without colors, not needed for MOs)
         do labs = 1, nvir_l
            call basis_vir%set_basis_func(labs,BASIS_ABSTRACT,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
         enddo

         !Finalize the mo space basis
         call basis_occ%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')
         call basis_vir%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

         !Register occ mo vector space:
         branch_factor = max(2, (int(get_num_segments(int(nocc,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('Occupied',basis_occ,occ_space_id,occ_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(occ) failed!')
         if(debug_exa .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') ' ExaCorr: Registered OCC vector space of dimension  ',nocc, &
         ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1))

         !Register vir mo vector space:
         branch_factor = max(2, (int(get_num_segments(int(nvir,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('Virtual',basis_vir,vir_space_id,vir_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(vir) failed!')
         if(debug_exa .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') ' ExaCorr: Registered VIRT vector space of dimension ',nvir, &
         ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1)) !debug

        end subroutine exacorr_cc_spaces

        subroutine exacorr_cc_methods (exa_input,ao_space_id,occ_space_id,vir_space_id)

!        This subroutine initalizes and registers all user-defined methods for ExaTensor

!        Written by Lucas Visscher, January 2019

         use exacorr_global
         use exacorr_tensor_methods
         use exacorr_ao_to_mo

         implicit none

         ! container with user input (we only need the indices of the occupied and virtual MOs here)
         type(exacc_input), intent(in) :: exa_input
         ! identifiers of the spaces that we will use in the methods
         integer(INTD), intent(in)     :: ao_space_id, occ_space_id, vir_space_id

         integer                :: nocc,nvir           ! the size of the mo basis for occupied and virtual orbitals (Kramers pairs)
         integer, allocatable   :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         integer, allocatable   :: mo_oo(:),mo_vo(:),mo_ov(:),mo_vv(:) ! joined lists needed for dm init (ugly code, should be refactored)

         integer(INTD)                :: ierr, spin

         ! types of the initializer methods
         type(tens_printer_t)         :: tens_printer         ! generic method, but needs to be defined here in specific form
         type(compute_2e_ao_tensor_t) :: aoint_calculator     ! method to initialize the ao-integrals tensor
         type(init_mocoef_tensor_t)   :: occupied_init(2), virtual_init(2) ! methods to initialize MO coeffcient tensors
         type(init_dm_tensor_t)       :: dm_init(4)           ! array of initializers for MO transformation algorithm 2
         type(scale_with_denom_t)     :: denominate           ! method to scale CC amplitudes by the inverse of orbital energy differences
         type(set_diagonal_t)         :: one_diag             ! identity matrix
         type(init_one_el_t)          :: one_el(3)            ! one electron integrals
         type(set_energy_t)           :: orb_set              ! place orbital energies on diagonal

         integer:: min_occ

         !Copy user input to local arrays for the definition of MO spaces (for standard CC we only distinguish occupied and virtual spinors)
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir

         min_occ=minval(mo_occ)

         !Define separate initializers for the alpha(A) and beta(B) occupied and virtual MO-coefficients
         do spin = 1, 2
            !Register this initializer for occupied as a method
            call occupied_init(spin)%init_mocoef_ctor(mo_occ,nocc,spin)
            ierr=exatns_method_register(mocoef_initlabel(ao_space_id,occ_space_id,spin),occupied_init(spin))
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            !Same procedure for the virtual space
            call virtual_init(spin)%init_mocoef_ctor(mo_vir,nvir,spin)
            ierr=exatns_method_register(mocoef_initlabel(ao_space_id,vir_space_id,spin),virtual_init(spin))
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         end do

         if (exa_input%moint_scheme == 2) then
           !Define names for the initializers for density matrices
           allocate(mo_oo(nocc+nocc))
           allocate(mo_vo(nvir+nocc))
           allocate(mo_ov(nocc+nvir))
           allocate(mo_vv(nvir+nvir))
           mo_oo(1:nocc) = mo_occ
           mo_ov(1:nocc) = mo_occ
           mo_vo(1:nvir) = mo_vir
           mo_vv(1:nvir) = mo_vir
           mo_oo(nocc+1:nocc+nocc) = mo_occ
           mo_ov(nocc+1:nocc+nvir) = mo_vir
           mo_vo(nvir+1:nvir+nocc) = mo_occ
           mo_vv(nvir+1:nvir+nvir) = mo_vir
           !Register the dm initializers as methods
           call dm_init(1)%init_dm_ctor(mo_oo,(/nocc,nocc/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/occ_space_id,occ_space_id/)),dm_init(1))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
           call dm_init(2)%init_dm_ctor(mo_vo,(/nvir,nocc/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/vir_space_id,occ_space_id/)),dm_init(2))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
           call dm_init(3)%init_dm_ctor(mo_ov,(/nocc,nvir/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/occ_space_id,vir_space_id/)),dm_init(3))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
           call dm_init(4)%init_dm_ctor(mo_vv,(/nvir,nvir/))
           ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),(/vir_space_id,vir_space_id/)),dm_init(4))
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         end if

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(aoint) failed!')

         !Register the denominator initializer as a method
         call denominate%scale_with_denom_ctor(mo_occ,nocc,mo_vir,nvir,exa_input%level_shift)
         ierr=exatns_method_register('Denominate',denominate)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(Denominate) failed!')

         !Register the standard print method with treshold such that it only prints interesting elements
         call tens_printer%reset_thresh(1.0D-6)
         ierr=exatns_method_register('PrintTensor',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(printer) failed!')

         !Register method to set diagonal to one
         call one_diag%init_const_ctor(ONE)
         ierr=exatns_method_register('OneDiag',one_diag)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OneDiag) failed!')

         one_el_exist = exist_one_el(e_core)
         if (one_el_exist) then
            !Register methods to read one electron integrals
            call one_el(1)%init_one_el_int(mo_occ,nocc,mo_occ,nocc,min_occ)
            ierr=exatns_method_register('OE_oo',one_el(1))
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_oo) failed!')
            call one_el(2)%init_one_el_int(mo_occ,nocc,mo_vir,nvir,min_occ)
            ierr=exatns_method_register('OE_ov',one_el(2))
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_ov) failed!')
            call one_el(3)%init_one_el_int(mo_vir,nvir,mo_vir,nvir,min_occ)
            ierr=exatns_method_register('OE_vv',one_el(3))
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_vv) failed!')

            !Register methods to make matrix with orbital energies on the diagonal
            call orb_set%init_eps(mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,exa_input%level_shift)
            ierr=exatns_method_register('OrbSet',orb_set)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(OrbSet) failed!')
         end if

        end subroutine exacorr_cc_methods

        subroutine solve_t_equations(comm_t,int_t,dims,CCD,t_target_precision,ccsd_energy,ncycles)

!        Routine for solving the t amplitude equations

         use exacorr_ao_to_mo, only : print_tensornorm2
         use exacorr_global,   only : debug_exa

!        Input
         type(exatns_comm_tens_t), intent(inout) :: comm_t
         type(exatns_intg_tens_t), intent(inout) :: int_t
         type(space_dims),intent(in) :: dims
         integer :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         logical, intent(in) :: CCD
         real(8), intent(in) :: t_target_precision
         real(8), intent(out) :: ccsd_energy
!        solution tensors
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!        intermediate
         type(tens_rcrsv_t) :: tau2_tensor, cint_tensor
!        auxiliary tensors
         type(tens_rcrsv_t) :: hov_aux_tensor, hoo_aux_tensor
!        result tensor
         type(tens_rcrsv_t) :: result_t_tensor
         complex(8)         :: result_val
!        Ids and roots for tensors
         integer(INTD),dimension(2) :: oo_id,ov_id,vo_id
         integer(INTD),dimension(4) :: ovoo_id,vvoo_id
         integer(INTL),dimension(2) :: oo_root,ov_root,vo_root
         integer(INTL),dimension(4) :: ovoo_root,vvoo_root

         real(8) :: t_convergence

         integer :: iccsd
         integer,intent(in) :: ncycles
         integer(INTD) :: ierr

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         if (.not.CCD) then
            oo_id     = (/occ_space_id, occ_space_id/)
            ov_id     = (/occ_space_id, vir_space_id/)
            vo_id     = (/vir_space_id, occ_space_id/)
            oo_root   = (/occ_space_root, occ_space_root/)
            ov_root   = (/occ_space_root, vir_space_root/)
            vo_root   = (/vir_space_root, occ_space_root/)
         end if
         ovoo_id   = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vvoo_id   = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

!        Initialize scalar that are to be used as tensors in contractions
         ierr=exatns_tensor_create(result_t_tensor,"result_t",EXA_DATA_KIND_C8)

!        Construct intermediate tensors
         if (.not.CCD) then
            ierr=exatns_tensor_create(hoo_aux_tensor,"hoo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hov_aux_tensor,"hov_aux",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(s1_tensor,"s1",vo_id,vo_root,EXA_DATA_KIND_C8)
         end if
         ierr=exatns_tensor_create(s2_tensor,"s2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(tau2_tensor,"tau2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(cint_tensor,"cint",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)

         do iccsd = 1, ncycles

          call get_tau(comm_t,CCD,one_tensor)

          ierr=exatns_tensor_init(s2_tensor)

!         The vvvv contribution
          ierr=exatns_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,int_t%vvvv,comm_t%tau,ONE_HALF) !o(2)v(4)

!         vv
          call get_hvv(comm_t,int_t%oovv,one_tensor)
          call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)

!         vovo
          ierr=exatns_tensor_init(tau2_tensor)
          ierr=exatns_tensor_contract("T(a,d,l,i)+=R(a,d,l,i)*K()",tau2_tensor,comm_t%t2,one_tensor,ONE_HALF)
          if (.not.CCD) ierr=exatns_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,comm_t%t1,comm_t%t1) ! o(2)v(2)
          call get_h(comm_t,int_t,CCD,one_tensor)
          ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",comm_t%h_int,int_t%oovv,tau2_tensor)! o(3)v(3)

!         This term appears four times as we need to permute both i and j and a and b
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,comm_t%h_int,comm_t%t2)  ! o(3)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,comm_t%h_int,comm_t%t2)  ! o(3)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,comm_t%h_int,comm_t%t2)  ! o(3)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,comm_t%h_int,comm_t%t2)  ! o(3)v(3)

!         ovoo

          ierr=exatns_tensor_init(cint_tensor)
          ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,int_t%vovv,comm_t%tau,MINUS_ONE_HALF) ! o(3)v(3)
          if (.not.CCD) then
              ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,int_t%vovo,comm_t%t1) ! o(3)v(2)
              ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,int_t%vovo,comm_t%t1,MINUS_ONE) ! o(3)v(2)
          end if
          ierr=exatns_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,int_t%ooov,one_tensor) ! o(3)v(1)

!         oo
          call get_hoo(comm_t,int_t%oovv,one_tensor)
          call get_goo(comm_t,int_t%ooov,CCD,one_tensor)

!         ov
          if (.not.CCD) then
            call get_hov(comm_t,int_t%oovv,one_tensor)
            ierr=exatns_tensor_init(hov_aux_tensor)
            ierr=exatns_tensor_contract("L(k,c)+=H(k,c)*K()",hov_aux_tensor,comm_t%hov,one_tensor)        ! copy common term
            ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hov_aux_tensor,comm_t%fov,one_tensor,MINUS_TWO) ! subtract fock matrix
          end if

!         oooo
          call get_a(comm_t,int_t,CCD,one_tensor)

          if (ierr.ne.EXA_SUCCESS) then
             print*," ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

          if (debug_exa) then
             if (.not.CCD) write(*,'(A7,F20.15)') "T1amp = ",print_tensornorm2 (comm_t%t1)
             write(*,'(A7,F20.15)') "T2amp = ",print_tensornorm2 (comm_t%t2)
             write(*,'(A7,F20.15)') "Tau   = ",print_tensornorm2 (comm_t%tau)
             write(*,'(A7,F20.15)') "hoint = ",print_tensornorm2 (comm_t%hoo)
             write(*,'(A7,F20.15)') "goint = ",print_tensornorm2 (comm_t%goo)
             write(*,'(A7,F20.15)') "hvint = ",print_tensornorm2 (comm_t%hvv)
             write(*,'(A7,F20.15)') "gvint = ",print_tensornorm2 (comm_t%gvv)
             if (.not.CCD) write(*,'(A7,F20.15)') "hovint= ",print_tensornorm2 (comm_t%hov)
             write(*,'(A7,F20.15)') "aint  = ",print_tensornorm2 (comm_t%a_int)
             write(*,'(A7,F20.15)') "cint  = ",print_tensornorm2 (cint_tensor)
             write(*,'(A7,F20.15)') "hint  = ",print_tensornorm2 (comm_t%h_int)
          end if

!         This term appears twice as we need to permute i and j
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,comm_t%goo,comm_t%t2)                   ! o(3)v(2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,comm_t%goo,comm_t%t2)                   ! o(3)v(2)
!         Again two times, also note the complex conjugation as we use vovv instead of vvvo
          if (.not.CCD) then
              ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,int_t%vovv,comm_t%t1)                   ! o(2)v(3)
              ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,int_t%vovv,comm_t%t1,MINUS_ONE)         ! o(2)v(3)
          end if

!         This term appears twice as we need to permute a and b
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! o(2)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! o(2)v(3)
!         Again two times, with the C intermediate (note: in RELCCSD this is done differently)
          if (.not.CCD) then
              ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,comm_t%t1,MINUS_ONE) ! o(3)v(2)
              ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,comm_t%t1)           ! o(3)v(2)
          end if

!         The oooo contribution
          ierr=exatns_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,comm_t%a_int,comm_t%tau,ONE_HALF) !o(4)v(2)

          ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)    ! o(2)v(2)

!         calculate s1 tensor
          if (.not.CCD) then
            ierr=exatns_tensor_init(s1_tensor)
            ierr=exatns_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,ONE_HALF) !o(2)v(3)
            ierr=exatns_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                   !o(1)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                   !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)               !o(2)v(2)
            ierr=exatns_tensor_init(hoo_aux_tensor)
            ierr=exatns_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)          !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)               !o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,ONE_HALF) !o(3)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,MINUS_ONE)     !o(2)v(2)
          end if

          if (ierr.ne.EXA_SUCCESS) then
             print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

!         prepare for next iteration
          if (.not.CCD) then
              ierr=exatns_tensor_init(comm_t%t1,ZERO)
              ierr=exatns_tensor_contract("T1(a,i)+=S1(a,i)*K()",comm_t%t1,s1_tensor,one_tensor)
              ierr=exatns_tensor_transform(comm_t%t1,'Denominate')
          end if
          ierr=exatns_tensor_init(comm_t%t2,ZERO)
          ierr=exatns_tensor_contract("T2(a,b,i,j)+=S2(a,b,i,j)*K()",comm_t%t2,s2_tensor,one_tensor)
          ierr=exatns_tensor_transform(comm_t%t2,'Denominate')

          if (.not.CCD) then
              call diis (comm_t%t2,vvoo_id,vvoo_root,comm_t%t1,vo_id,vo_root,t_convergence)
          else
              call diis (comm_t%t2,vvoo_id,vvoo_root,t_convergence=t_convergence)
          end if

          ierr=exatns_tensor_init(result_t_tensor)
          if (.not.CCD) then
!             Will always be 0 contribution for now, since fov=0 (canonical orbitals)
              ierr=exatns_tensor_contract("E()+=T(a,i)*F(i,a)",result_t_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
              ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
              ccsd_energy=real(result_val,8)
              ierr=exatns_tensor_init(result_t_tensor)
          else
              ccsd_energy=real(ZERO,8)
          end if

!         Update tau with converged values
          call get_tau(comm_t,CCD,one_tensor)

          ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_t_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
          ccsd_energy = ccsd_energy + real(result_val*ONE_QUARTER,8)

          call print_iteration(iccsd,t_convergence,ccsd_energy)
          if (t_convergence .lt. t_target_precision) exit
         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
         write(*,*) ""

         if (.not.CCD) ierr=exatns_tensor_destroy(s1_tensor)
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(tau2_tensor)
         ierr=exatns_tensor_destroy(cint_tensor)
         if (.not.CCD) then
            ierr=exatns_tensor_destroy(hoo_aux_tensor)
            ierr=exatns_tensor_destroy(hov_aux_tensor)
         end if
         ierr=exatns_tensor_destroy(result_t_tensor)

        end subroutine solve_t_equations

        subroutine solve_lambda_equations(comm_t,int_t,l2_tensor,dims,t_target_precision,ncycles,print_level,l1_tensor)

!        Routine for solving the lambda equations

         use exacorr_ao_to_mo, only : print_tensornorm2
         use exacorr_global,   only : debug_exa

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_comm_tens_t), intent(inout) :: comm_t
         type(exatns_intg_tens_t), intent(inout) :: int_t
         type(exatns_lambda_tens_t)              :: lambda_t
!        solution tensors
         type(tens_rcrsv_t), intent(inout)           :: l2_tensor
         type(tens_rcrsv_t), intent(inout), optional :: l1_tensor
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!        Tensor dimensions
         type(space_dims),intent(in) :: dims
         integer(INTD) :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         integer(INTD), dimension(2) :: oo_id, ov_id, vv_id
         integer(INTD), dimension(4) :: oooo_id,ovoo_id,ooov_id,vovo_id,oovv_id, &
                                  vvvo_id,vovv_id,vvvv_id
         integer(INTL), dimension(2) :: oo_root, ov_root, vv_root
         integer(INTL), dimension(4) :: oooo_root,ovoo_root,ooov_root,vovo_root,oovv_root, &
                                        vvvo_root,vovv_root,vvvv_root
         integer :: ilambda
         integer,intent(in) :: ncycles
         integer,intent(in) :: print_level
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         logical :: CCD=.false.
         integer(INTD) :: ierr

         write(*,*) ""
         write(*,'(A)') "*******************************************************************************"
         write(*,'(A)') "********************** Solving Left Eigenvector Equation **********************"
         write(*,'(A)') "*******************************************************************************"
         write(*,*) ""

         if (.not.present(l1_tensor)) CCD = .true.

!**********************************
!        create all ids and roots *
!**********************************

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oooo_id = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
         ovoo_id = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         ooov_id = (/occ_space_id, occ_space_id, occ_space_id, vir_space_id/)
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvvo_id = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!***************************************
!        construct fixed intermediates *
!***************************************

         ierr=exatns_tensor_create(lambda_t%fbar_oo,"fbar_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%fbar_ov,"fbar_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%fbar_vv,"fbar_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_oooo,"w_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_ovoo,"w_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_ooov,"w_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vovo,"w_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%wbar_vovo,"wbar_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vvvo,"w_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vovv,"w_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vvvv,"w_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)

!******************************************
!        obtain intermediates from CC(S)D *
!******************************************

         call get_hvv(comm_t,int_t%oovv,one_tensor)
         call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)
         call get_h(comm_t,int_t,CCD,one_tensor)
         call get_hoo(comm_t,int_t%oovv,one_tensor)
         call get_goo(comm_t,int_t%ooov,CCD,one_tensor)
         call get_hov(comm_t,int_t%oovv,one_tensor)
         call get_a(comm_t,int_t,CCD,one_tensor)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD,one_tensor)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t,one_tensor)
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov,one_tensor)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD,one_tensor)
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,occ_space_id,occ_space_root, &
                          vir_space_id,vir_space_root,CCD,one_tensor)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int,one_tensor)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD,one_tensor)
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv,one_tensor)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD,one_tensor)
         call print_date('Finished calculating fixed intermediates')

         if (print_level > 0) then
           if (.not.CCD) print*, "T1amp      = ", print_tensornorm2(comm_t%t1)
           print*, "T2amp      = ", print_tensornorm2(comm_t%t2)*0.25D0
           ! print*, "Tau        = ", print_tensornorm2(comm_t%tau)*0.25D0
           ! print*, "Fbar_oo    = ", print_tensornorm2(lambda_t%fbar_oo)*0.5D0
           print*, "Fbar_ov    = ", print_tensornorm2(lambda_t%fbar_ov)
           ! print*, "Fbar_vv    = ", print_tensornorm2(lambda_t%fbar_vv)*0.5D0
           ! print*, "w_vvvv     = ", print_tensornorm2(lambda_t%w_vvvv)*0.25D0
           print*, "wbar_vovo  = ", print_tensornorm2(lambda_t%wbar_vovo)
           print*, "w_vvvo     = ", print_tensornorm2(lambda_t%w_vvvo)*0.5D0
           print*, "w_vovv     = ", print_tensornorm2(lambda_t%w_vovv)*0.5D0
           print*, "w_vovo     = ", print_tensornorm2(lambda_t%w_vovo)
           print*, "w_ovoo     = ", print_tensornorm2(lambda_t%w_ovoo)*0.5D0
           print*, "w_ooov     = ", print_tensornorm2(lambda_t%w_ooov)*0.5D0
           print*, "w_oooo     = ", print_tensornorm2(lambda_t%w_oooo)*0.25D0
         end if

!************************************
!        construct solution tensors *
!************************************

         if (.not.CCD) then
            ierr=exatns_tensor_create(s1_tensor,"s1",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(s1_tensor,ZERO)
         end if 
         ierr=exatns_tensor_create(s2_tensor,"s2",oovv_id,oovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_init(s2_tensor,ZERO)

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=exatns_tensor_create(lambda_t%goo,"goo_mixed",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%gvv,"gvv_mixed",vv_id,vv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_init(lambda_t%goo,ZERO)
         ierr=exatns_tensor_init(lambda_t%gvv,ZERO)

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (l2_tensor,oovv_id,oovv_root,l1_tensor,ov_id,ov_root,mode=1)
         else
           call diis (l2_tensor,oovv_id,oovv_root,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)

!         debug
          if (debug_exa) then
            if (.not.CCD) print*, "L1  = ", print_tensornorm2(l1_tensor)
            print*, "L2  = ", print_tensornorm2(l2_tensor)*0.25D0
            print*, "goo = ", print_tensornorm2(lambda_t%goo)
            print*, "gvv = ", print_tensornorm2(lambda_t%gvv)
          end if

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=exatns_tensor_init(s1_tensor,ZERO)

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=F(i,a)*K()",s1_tensor,lambda_t%fbar_ov,one_tensor)
!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=exatns_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=exatns_tensor_init(s2_tensor,ZERO)

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |
!----------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |
!----------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |
!------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |
!----------------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!-------------------------------------------------------
!         term 6: S(ij,ab) += P(a,b) V(ij,ae) * G(e,b) |
!-------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,lambda_t%gvv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,lambda_t%gvv)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(,j,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)
          end if

!-------------------------------------------------------
!         term 8: S(ij,ab) -= P(i,j) V(im,ab) * G(j,m) |
!-------------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,lambda_t%goo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,lambda_t%goo)

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
            ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |
!-----------------------------------------------------

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,ONE_HALF)

          if (ierr.ne.EXA_SUCCESS) then
            print*," ccdriver: error in evaluation of lambda equations, code=",ierr
            return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=exatns_tensor_init(l1_tensor)
            ierr=exatns_tensor_contract("L1(i,a)+=S1(i,a)*K()",l1_tensor,s1_tensor,one_tensor)
            ierr=exatns_tensor_transform(l1_tensor,"Denominate")
          end if
          ierr=exatns_tensor_init(l2_tensor)
          ierr=exatns_tensor_contract("L2(i,j,a,b)+=S2(i,j,a,b)*K()",l2_tensor,s2_tensor,one_tensor)
          ierr=exatns_tensor_transform(l2_tensor,"Denominate")

          if (.not.CCD) then
            call diis (l2_tensor,t1_tensor=l1_tensor,t_convergence=t_convergence)
          else
            call diis (l2_tensor,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", print_tensornorm2(l1_tensor)
         write(*,'(A,ES22.16)') "  Final L2amp = ", print_tensornorm2(l2_tensor)*0.25D0
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=exatns_tensor_destroy(s1_tensor)
           call diis(l2_tensor,t1_tensor=l1_tensor,mode=2)
         else
           call diis(l2_tensor,mode=2)
         end if
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(one_tensor)
         ierr=exatns_tensor_destroy(lambda_t%fbar_oo)
         ierr=exatns_tensor_destroy(lambda_t%fbar_ov)
         ierr=exatns_tensor_destroy(lambda_t%fbar_vv)
         ierr=exatns_tensor_destroy(lambda_t%w_oooo)
         ierr=exatns_tensor_destroy(lambda_t%w_ooov)
         ierr=exatns_tensor_destroy(lambda_t%w_ovoo)
         ierr=exatns_tensor_destroy(lambda_t%w_vovo)
         ierr=exatns_tensor_destroy(lambda_t%w_vvvo)
         ierr=exatns_tensor_destroy(lambda_t%w_vovv)
         ierr=exatns_tensor_destroy(lambda_t%w_vvvv)
         ierr=exatns_tensor_destroy(lambda_t%wbar_vovo)
         ierr=exatns_tensor_destroy(lambda_t%goo)
         ierr=exatns_tensor_destroy(lambda_t%gvv)

        end subroutine solve_lambda_equations

        subroutine get_CC_integrals (ao_space_id,ao_space,int_t,moint_scheme)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017

         use exacorr_ao_to_mo
         use exacorr_global

!        information needed to compute the ao-tensor
         integer(INTD),intent(in)              :: ao_space_id
         class(h_space_t), pointer, intent(in) :: ao_space
!        target integral tensors to be filled and given back to the caller
         type(exatns_intg_tens_t), intent(inout) :: int_t 

         integer(INTL)              :: ao_space_root,subspace_id(4),num_subspaces
         integer(INTL),allocatable  :: ao_subspaces(:)
         type(tens_rcrsv_t)         :: a2int_tensor, ovvo_tensor ! Auxilliary tensors
         integer(INTL),dimension(4) :: a2int_root, vovo_root, ovvo_root
         integer(INTD),dimension(4) :: a2int_id, vovo_id, ovvo_id
         
         integer(INTD)                     :: ierr
         integer                           :: i, sub3, sub4
         integer                           :: moint_scheme ! (0: put zeroes for testing. 1: n^5 algorithm, 2: n^6 algorithm
                                                           !  3: n^5 with re-use of ints 4: as 3, but with sliced AO and HT)
         character(len=15)                 :: text

!        auxilliary integral tensor (only needed inside this routine)
         type(tens_rcrsv_t) :: aoint_tensor

!        get the root space, covers all ao's
         ao_space_root=ao_space%get_root_id(ierr)

!        Determine space ids and root for auxilliary voov tensor from vovo tensor
         call int_t%vovo%get_space_ids(vovo_id,vovo_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('vovo tensor corrupted')
         ovvo_id(1)   = vovo_id(2)
         ovvo_id(2)   = vovo_id(1)
         ovvo_id(3:4) = vovo_id(3:4)
         ovvo_root(1)   = vovo_root(2)
         ovvo_root(2)   = vovo_root(1)
         ovvo_root(3:4) = vovo_root(3:4)

!        Get tensor with AO integrals
         if (moint_scheme /= 4) then
            ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            call print_date('allocated AO integral tensor')
         end if

         !Initialize ao-tensor:
         if (moint_scheme == 0) then
!        test code : initialize integrals to zero to be able to time the integral generation
           ierr=exatns_tensor_init(aoint_tensor)
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('initalized all ao integrals to zero')
         elseif (moint_scheme == 1 .or. moint_scheme == 2) then
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('finished ao integral calculation')

!          Get vvvv tensor
           call get_integral_tensor (aoint_tensor,int_t%vvvv,12,algorithm=moint_scheme)
           call print_date('finished vvvv integrals')

!          Get vovv tensor
           call get_integral_tensor (aoint_tensor,int_t%vovv,34,algorithm=moint_scheme)
           call print_date('finished vovv integrals')

!          Get oovv tensor
           call get_integral_tensor (aoint_tensor,int_t%oovv,12,algorithm=moint_scheme)
           call print_date('finished oovv integrals')

!          Get ooov tensor
           call get_integral_tensor (aoint_tensor,int_t%ooov,12,algorithm=moint_scheme)
           call print_date('finished ooov integrals')

!          Get oooo tensor
           call get_integral_tensor (aoint_tensor,int_t%oooo,12,algorithm=moint_scheme)
           call print_date('finished oooo integrals')

!          Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!          Get vovo tensor
           call get_integral_tensor (aoint_tensor,int_t%vovo,0,algorithm=moint_scheme)
           ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
           call get_integral_tensor (aoint_tensor,ovvo_tensor,0,algorithm=moint_scheme)
!          Make anti-symmetrized vovo tensor
           ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
           ierr=exatns_tensor_destroy(ovvo_tensor)
           call print_date('finished vovo integrals')

         elseif (moint_scheme == 3) then
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('finished ao integral calculation')
           call print_date('choosing moint_scheme 3')

           ! Setting dimensions for the half transformed integral tensor, the last two are always AO spaces
           a2int_id(3:4)   = ao_space_id
           a2int_root(3:4) = ao_space_root

           ! get the integrals with vv as first transform
           a2int_id(1:2)   = vovo_id(1)
           a2int_root(1:2) = vovo_root(1)
           ierr=exatns_tensor_create(a2int_tensor,"a2vv",a2int_id,a2int_root,EXA_DATA_KIND_C8) 
           call get_ht1_integral_tensor (aoint_tensor,a2int_tensor)
           call print_date('vv halftransformation completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vvvv,12)
           call print_date('vvvv integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vovv,34)
           call print_date('vovv integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vovo,0)
           call print_date('vovo integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

           ! get the integrals with ov as first transform
           a2int_id(1:2)   = vovo_id(2:3)
           a2int_root(1:2) = vovo_root(2:3)
           ierr=exatns_tensor_create(a2int_tensor,"a2ov",a2int_id,a2int_root,EXA_DATA_KIND_C8) 
           call get_ht1_integral_tensor (aoint_tensor,a2int_tensor)

           call print_date('ov halftransformation completed')
           ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
           call get_ht2_integral_tensor (a2int_tensor,ovvo_tensor,0)
           call print_date('ovvo integrals completed')
!          Make anti-symmetrized vovo tensor
           ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
           ierr=exatns_tensor_destroy(ovvo_tensor)
           call print_date('vovo integrals anti-symmetrized')
           call get_ht2_integral_tensor (a2int_tensor,int_t%oovv,12)
           call print_date('oovv integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

           ! get the integrals with oo as first transform
           a2int_id(1:2)   = vovo_id(2)
           a2int_root(1:2) = vovo_root(2)
           ierr=exatns_tensor_create(a2int_tensor,"a2oo",a2int_id,a2int_root,EXA_DATA_KIND_C8) 
           call get_ht1_integral_tensor (aoint_tensor,a2int_tensor)
           call print_date('oo halftransformation completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%ooov,12)
           call print_date('ooov integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%oooo,12)
           call print_date('oooo integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

         elseif (moint_scheme == 4) then
           
           call print_date('choosing moint_scheme 4')
           ! in this algorithm we only keep part of the AO and HT integrals in memory
           subspace_id(1:2) = ao_space_root
           call ao_space%get_level_composition(int(1,INTD),ao_subspaces,num_subspaces,ierr)
           do sub4 = 1, num_subspaces
              subspace_id(4)=ao_subspaces(sub4)
              do sub3 = 1, num_subspaces
                 subspace_id(3)=ao_subspaces(sub3)
                 ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),subspace_id,EXA_DATA_KIND_C8)
                 if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
                 ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
                 if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')

                 ! Setting dimensions for the half transformed integral tensor, the last two are always AO spaces
                 a2int_id(3:4)   = ao_space_id
                 a2int_root(3:4) = subspace_id(3:4)

                 ! get the integrals with vv as first transform
                 a2int_id(1:2)   = vovo_id(1)
                 a2int_root(1:2) = vovo_root(1)
                 ierr=exatns_tensor_create(a2int_tensor,"a2vv",a2int_id,a2int_root,EXA_DATA_KIND_C8) 
                 call get_ht1_integral_tensor (aoint_tensor,a2int_tensor)
                 call get_ht2_integral_tensor (a2int_tensor,int_t%vvvv,12)
                 call get_ht2_integral_tensor (a2int_tensor,int_t%vovv,34)
                 call get_ht2_integral_tensor (a2int_tensor,int_t%vovo,0)
                 ierr=exatns_tensor_destroy(a2int_tensor)

                 ! get the integrals with ov as first transform
                 a2int_id(1:2)   = vovo_id(2:3)
                 a2int_root(1:2) = vovo_root(2:3)
                 ierr=exatns_tensor_create(a2int_tensor,"a2ov",a2int_id,a2int_root,EXA_DATA_KIND_C8) 
                 call get_ht1_integral_tensor (aoint_tensor,a2int_tensor)
                 ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
                 call get_ht2_integral_tensor (a2int_tensor,ovvo_tensor,0)
!                Make anti-symmetrized vovo tensor
                 ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
                 ierr=exatns_tensor_destroy(ovvo_tensor)
                 call get_ht2_integral_tensor (a2int_tensor,int_t%oovv,12)
                 ierr=exatns_tensor_destroy(a2int_tensor)

                 ! get the integrals with oo as first transform
                 a2int_id(1:2)   = vovo_id(2)
                 a2int_root(1:2) = vovo_root(2)
                 ierr=exatns_tensor_create(a2int_tensor,"a2oo",a2int_id,a2int_root,EXA_DATA_KIND_C8) 
                 call get_ht1_integral_tensor (aoint_tensor,a2int_tensor)
                 call get_ht2_integral_tensor (a2int_tensor,int_t%ooov,12)
                 call get_ht2_integral_tensor (a2int_tensor,int_t%oooo,12)
                 ierr=exatns_tensor_destroy(a2int_tensor)
                 ierr=exatns_tensor_destroy(aoint_tensor)
              end do
              write (text,'(I6,A,I6)') sub4,' of',num_subspaces
              call print_date('finished transforming subspace'//text)
           end do

         end if

!        Clean-up
         if (moint_scheme /= 4) then
            ierr=exatns_tensor_destroy(aoint_tensor)
         end if

        end subroutine get_CC_integrals

        subroutine diis (t2_tensor,vvoo_id,vvoo_root,t1_tensor,vo_id,vo_root,t_convergence,mode)

         integer, parameter             :: maxdim = 5 ! maximum number of tensors that will be stored
         integer, parameter             :: mindim = 3 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(tens_rcrsv_t), intent(inout),optional :: t1_tensor
         type(tens_rcrsv_t), intent(inout)          :: t2_tensor
         integer(INTD)                              :: ierr, tens_rank
         integer                                    :: i
         integer(INTD),dimension(2),optional        :: vo_id
         integer(INTD),dimension(4),optional        :: vvoo_id
         integer(INTL),dimension(2),optional        :: vo_root
         integer(INTL),dimension(4),optional        :: vvoo_root
         logical                                    :: init, clean_up, CCD
         integer, save                              :: itens, ntens
         complex(8), save                           :: b_diis(maxdim,maxdim)
         complex(8)                                 :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         complex(8)                                 :: scalar
         real(8)                                    :: b_max, b_scale
         type(tens_rcrsv_t), save                   :: t1_inp_tensor, t2_inp_tensor
         type(tens_rcrsv_t), save                   :: t1_res_tensor(maxdim), t2_res_tensor(maxdim)
         type(tens_rcrsv_t), save                   :: t1_err_tensor(maxdim), t2_err_tensor(maxdim)
         integer                                    :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)                                    :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(tens_rcrsv_t), save                   :: scalar_tensor

         character(len=1) :: vs_str

         init       = .false.
         clean_up   = .false.
         CCD        = .false.

         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if

         if (.not.present(t1_tensor)) CCD = .true.

         if (init) then
!           get dimensions of the tensors
            if (.not.CCD) then
                tens_rank=t1_tensor%get_rank(ierr)
                if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
            end if 
            tens_rank=t2_tensor%get_rank(ierr)
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            if (.not.CCD) ierr=exatns_tensor_create(t1_inp_tensor,"t1_inp_tensor",vo_id,vo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(t2_inp_tensor,"t2_inp_tensor",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            do i = 1, maxdim
               write (vs_str, '(I1)') i
               if (.not.CCD) then
                   ierr=exatns_tensor_create(t1_res_tensor(i),"t1_res_tensor_" // vs_str,vo_id,vo_root,EXA_DATA_KIND_C8)
                   ierr=exatns_tensor_create(t1_err_tensor(i),"t1_err_tensor_" // vs_str,vo_id,vo_root,EXA_DATA_KIND_C8)
               end if
               ierr=exatns_tensor_create(t2_res_tensor(i),"t2_res_tensor_" // vs_str,vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_create(t2_err_tensor(i),"t2_err_tensor_" // vs_str,vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            end do
            if(ierr.ne.EXA_SUCCESS) call quit('failed to create storage for input')
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            if (.not.CCD) then
                ierr=exatns_tensor_init(t1_inp_tensor)
                ierr=exatns_tensor_contract("T1O(a,i)+=T1N(a,i)*K()",t1_inp_tensor,t1_tensor,one_tensor)
            end if
            ierr=exatns_tensor_init(t2_inp_tensor)
            ierr=exatns_tensor_contract("T2O(a,b,i,j)+=T2N(a,b,i,j)*K()",t2_inp_tensor,t2_tensor,one_tensor)
            ierr=exatns_tensor_create(scalar_tensor,"scalar_tensor",EXA_DATA_KIND_C8)
            return ! no extrapolation possible
         end if

         if (clean_up) then
!           de-allocate tensors
            if (.not.CCD) ierr=exatns_tensor_destroy(t1_inp_tensor)
            ierr=exatns_tensor_destroy(t2_inp_tensor)
            do i = 1, maxdim
                if (.not.CCD) then
                    ierr=exatns_tensor_destroy(t1_res_tensor(i))
                    ierr=exatns_tensor_destroy(t1_err_tensor(i))
                end if
                ierr=exatns_tensor_destroy(t2_res_tensor(i))
                ierr=exatns_tensor_destroy(t2_err_tensor(i))
            end do
            ierr=exatns_tensor_destroy(scalar_tensor)
            return
         end if

!        update counters
         itens = itens + 1
         if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
         ntens = ntens + 1                    ! number of tensors stored after returning
         if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!        store the result tensors
         if (.not.CCD) then
            ierr=exatns_tensor_init(t1_res_tensor(itens))
            ierr=exatns_tensor_contract("T1R(a,i)+=T1(a,i)*K()",t1_res_tensor(itens),t1_tensor,one_tensor)
         end if
         ierr=exatns_tensor_init(t2_res_tensor(itens))
         ierr=exatns_tensor_contract("T2R(a,b,i,j)+=T2(a,b,i,j)*K()",t2_res_tensor(itens),t2_tensor,one_tensor)

!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         if (.not.CCD) then
            ierr=exatns_tensor_init(t1_err_tensor(itens))
            ierr=exatns_tensor_contract("T1E(a,i)+=T1R(a,i)*K()",t1_err_tensor(itens),t1_tensor,one_tensor)
            ierr=exatns_tensor_contract("T1E(a,i)+=T1I(a,i)*K()",t1_err_tensor(itens),t1_inp_tensor,one_tensor,prefactor=MINUS_ONE)
         end if
         ierr=exatns_tensor_init(t2_err_tensor(itens))
         ierr=exatns_tensor_contract("T2E(a,b,i,j)+=T2R(a,b,i,j)*K()",t2_err_tensor(itens),t2_tensor,one_tensor)
         ierr=exatns_tensor_contract("T2E(a,b,i,j)+=T2I(a,b,i,j)*K()",t2_err_tensor(itens),t2_inp_tensor,one_tensor, &
                                                                                                prefactor=MINUS_ONE)

!        update the b_diis array that holds the scalar products of the error tensors
         do i = 1, ntens
            ierr=exatns_tensor_init(scalar_tensor)
            if (.not.CCD) ierr=exatns_tensor_contract("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
            ierr=exatns_tensor_contract("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
            ierr=exatns_tensor_get_scalar(scalar_tensor,scalar)
            b_diis(itens,i)=scalar
            b_diis(i,itens)=dconjg(scalar)
         end do

         if (present(t_convergence)) then
            t_convergence = dsqrt(dble(b_diis(itens,itens)))
         end if 

!        transfer to a work array in the way requested by the solver and normalize the matrix
         b_max = maxval (abs(b_diis))
         if (b_max .gt. 1.D-30) then
             b_scale = 1.D0 / b_max
         else
             print*, "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
             return
         end if
         bb = ZERO
         bb(1:ntens,0) = MINUS_ONE
         bb(0,1:ntens) = MINUS_ONE
         bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!        solve the linear system of equations using the LAPACK solver ZHESV
         coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
         coeff(1:ntens) = ZERO 
         call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!        make the extrapolated amplitudes if the space is large enough
         if (ntens.ge.mindim) then
            if (.not.CCD) ierr=exatns_tensor_init(t1_tensor)
            ierr=exatns_tensor_init(t2_tensor)
            do i = 1, ntens
               if (.not.CCD) ierr=exatns_tensor_contract("T1(a,i)+=T1O(a,i)*K()",t1_tensor,t1_res_tensor(i), &
                                                                                one_tensor,prefactor=coeff(i))
               ierr=exatns_tensor_contract("T2(a,b,i,j)+=T2O(a,b,i,j)*K()",t2_tensor,t2_res_tensor(i),one_tensor,prefactor=coeff(i))
            end do
         end if
        
!        the extrapolated amplitudes are the new input tensors
         if (.not.CCD) then
            ierr=exatns_tensor_init(t1_inp_tensor)
            ierr=exatns_tensor_contract("T1O(a,i)+=T1(a,i)*K()",t1_inp_tensor,t1_tensor,one_tensor)
         end if
         ierr=exatns_tensor_init(t2_inp_tensor)
         ierr=exatns_tensor_contract("T2O(a,b,i,j)+=T2(a,b,i,j)*K()",t2_inp_tensor,t2_tensor,one_tensor)

        end subroutine diis

        subroutine init_fock (comm_t,int_t,oo_id,oo_root,ov_id,ov_root,vv_id,vv_root, print_level)

          use exacorr_global
          use exacorr_tensor_methods

          implicit none
  
          type(exatns_comm_tens_t), intent(inout)  :: comm_t
          type(exatns_intg_tens_t), intent(inout)  :: int_t
          integer(INTL),dimension(2), intent(in) :: oo_root, ov_root, vv_root
          integer(INTD),dimension(2), intent(in) :: oo_id, ov_id, vv_id 
          integer, intent(in)                    :: print_level

          real(8)                    :: ev_energy, off_diag, hf_energy
          !        scalars (need to be defined as tensor types)
          type(tens_rcrsv_t)         :: r_tensor, o_tensor, orb_tensor, h_tensor
          !        identity matrix for contraction 
          type(tens_rcrsv_t)         :: id_mat_occ
          integer(INTD)              :: ierr
          complex(8)                 :: result_val

          if (one_el_exist) then

            ! Allocate two auxilliary tensors needed to extract results or multiply by one
            ierr=exatns_tensor_create(r_tensor,"r_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(o_tensor,"o_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(o_tensor, ONE)

            ! Make the density matrix (which is an identity matrix in the MO basis)
            ierr=exatns_tensor_create(id_mat_occ,"id_mat_occ",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(id_mat_occ,ZERO)
            ierr=exatns_tensor_transform(id_mat_occ,'OneDiag')
            if (print_level >0) call print_date('Density matrix initialized')

            ! Initialize the sub blocks of the Fock matrix with the (effective) one-electron integrals
            ierr=exatns_tensor_init(comm_t%foo,'OE_oo')
            ierr=exatns_tensor_init(comm_t%fov,'OE_ov')
            ierr=exatns_tensor_init(comm_t%fvv,'OE_vv')
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'Creating one electron integrals failed!')
            if (print_level >0) call print_date('one electron integrals added to Fock Matrix')

            ! Calculate the Hartree-Fock energy
            ierr=exatns_tensor_init(r_tensor, ZERO)
            ierr=exatns_tensor_contract("R()=F(i,j)*I(i,j)",r_tensor,comm_t%foo,id_mat_occ)
            ierr=exatns_tensor_create(h_tensor,"h_tensor",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_get_scalar(r_tensor,result_val)
            ierr=exatns_tensor_init(h_tensor, ZERO)
            ierr=exatns_tensor_contract("H(i,j)+=V(i,k,j,l)*I(k,l)",h_tensor,int_t%oooo,id_mat_occ, ONE_HALF)
            ierr=exatns_tensor_contract("R()+=H(i,j)*I(i,j)",r_tensor,h_tensor,id_mat_occ)
            ierr=exatns_tensor_get_scalar(r_tensor,result_val)
            hf_energy=real(result_val,8) ! Casting into real(8)
            ierr=exatns_tensor_destroy(h_tensor)
            write (*,'(A,F23.15)') " Energy of reference wave function (HF single determinant) = ", hf_energy+e_core

            ! Make the complete Fock matrix by addding the valence electron contributions
            ierr=exatns_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",comm_t%foo,int_t%oooo,id_mat_occ, ONE)
            ierr=exatns_tensor_contract("F(i,a)+=V(k,i,l,a)*I(k,l)",comm_t%fov,int_t%ooov,id_mat_occ, ONE)
            ierr=exatns_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",comm_t%fvv,int_t%vovo,id_mat_occ, ONE)
            if (print_level >0) call print_date('two electron integrals added to Fock Matrix')

            if (print_level >3) ierr=exatns_tensor_traverse(comm_t%foo,'PrintTensor')
            if (print_level >3) ierr=exatns_tensor_traverse(comm_t%fov,'PrintTensor')
            if (print_level >3) ierr=exatns_tensor_traverse(comm_t%fvv,'PrintTensor')

            ! Subtract orbital energies from the diagonal (they are moved to the other side of the CC equations)
            ierr=exatns_tensor_create(orb_tensor,"orb_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(orb_tensor,'OrbSet')
            ierr=exatns_tensor_contract("F(i,j)+=E(i,j)*O()",comm_t%foo,orb_tensor,o_tensor, MINUS_ONE)
            ierr=exatns_tensor_destroy(orb_tensor)
            ierr=exatns_tensor_create(orb_tensor,"orb_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(orb_tensor,'OrbSet')
            ierr=exatns_tensor_contract("F(a,b)+=E(a,b)*O()",comm_t%fvv,orb_tensor,o_tensor, MINUS_ONE)
            ierr=exatns_tensor_destroy(orb_tensor)
            if (print_level >0) call print_date('orbital energies subtracted from the diagonal')

            ! Compute norms (check whether this single determinant is contructed from optimized orbitals)
            ierr=exatns_tensor_init(r_tensor, ZERO)
            ierr=exatns_tensor_contract("R()=F(i,j)*F(i,j)",r_tensor,comm_t%foo,comm_t%foo)
            ierr=exatns_tensor_get_scalar(r_tensor,result_val)
            off_diag=real(result_val,8) ! Casting into real(8)
            if (print_level >0) print*," Norm of OO part Fock matrix (without diag) = ", off_diag
            
            ierr=exatns_tensor_init(r_tensor, ZERO)
            ierr=exatns_tensor_contract("R()=F(i,a)*F(i,a)",r_tensor,comm_t%fov,comm_t%fov)
            ierr=exatns_tensor_get_scalar(r_tensor,result_val)
            off_diag=real(result_val,8) ! Casting into real(8)
            if (print_level >0) print*," Norm of OV part Fock matrix =                ", off_diag

            ierr=exatns_tensor_init(r_tensor, ZERO)
            ierr=exatns_tensor_contract("R()=F(a,b)*F(a,b)",r_tensor,comm_t%fvv,comm_t%fvv)
            ierr=exatns_tensor_get_scalar(r_tensor,result_val)
            off_diag=real(result_val,8) ! Casting into real(8)
            if (print_level >0) print*," Norm of VV part Fock matrix (without diag) = ", off_diag

            ierr=exatns_tensor_destroy(h_tensor)
            ierr=exatns_tensor_destroy(id_mat_occ)
            ierr=exatns_tensor_destroy(r_tensor)
            ierr=exatns_tensor_destroy(o_tensor)
          else
            ierr=exatns_tensor_init(comm_t%foo,ZERO)
            ierr=exatns_tensor_init(comm_t%fov,ZERO)
            ierr=exatns_tensor_init(comm_t%fvv,ZERO)

          end if

        end subroutine init_fock

#else
        implicit none

        public exacorr_cc_driver

       contains

        subroutine exacorr_cc_driver (exa_input)
         use exacorr_datatypes
         implicit none
         type(exacc_input), intent(in) :: exa_input
         print*,  'called exatensor version of cc_driver in a serial version'
         print*,  'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_cc_driver

#endif

end module exacorr_cc
