module exacorr_cc

!This module contains the routines needed to calculate a coupled cluster wave function
!almost everything requires MPI, but we provide a stub function to make compilation of 
!non-MPI versions of DIRAC possible.

#if defined (VAR_MPI)

        use exatensor
        use exacorr_datatypes
        use intermediates

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0),ONE_QUARTER=(0.25D0,0.D0)!, &
!                                THREE=(3.0,0.0), SIX=(6.0,0.0)
        type(tens_rcrsv_t),public :: one_tensor

        private

        public exacorr_cc_driver

        interface diis
           module procedure diis
        end interface

       contains

!-------------------------------------------------------------------------------------------
        subroutine exacorr_cc_driver (exa_input)

!        This subroutine drives the ccsd energy calculation
!        Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code

!        Written by Lucas Visscher, winter 2016/2017

         use interface_to_mpi
         use exacorr_global
         use exacorr_tensor_methods
         use exacorr_ao_to_mo
         use exacorr_gradient
         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir           ! the size of the mo basis for occupied and virtual spinors
         integer(INTL)          :: nocc_l, nvir_l
         integer, allocatable   :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer :: my_role, my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode
         integer                                :: nshells  ! the number of shells (orbitals with the same center and l-value)
         integer                                :: nao      ! the number of atomic orbitals)
         type(basis_func_info_t), allocatable   :: gto(:)   ! array with information about the shells

         integer                        :: ierr, i
         complex(8)                     :: result_val

!        Common tensors
         type(exatns_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_intg_tens) :: int_t
         type(tens_rcrsv_t) :: voov_tensor
!        Lambda tensors
         type(tens_rcrsv_t) :: l1_tensor, l2_tensor
!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t)         :: result_tensor
!        space ids and roots
         type(space_dims) :: dims
!        CCSD control and result variables
         real(8)            :: mp2_energy, ccsd_energy!, t_energy(3) ! triples not implemented
         real(8)            :: t_target_precision

         ! classes for the definition of spaces (AO and 4 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         class(h_space_t), pointer                :: occ_space, vir_space

         ! variables for the definition of spaces (AO and 2 different MO spaces)
         type(subspace_basis_t)     :: basis_ao, basis_occ, basis_vir
         integer(INTL)              :: nao_l, ao_space_root
         integer(INTL)              :: occ_space_root, vir_space_root
         integer(INTL),dimension(2) :: oo_root, ov_root, vo_root, vv_root
         integer(INTL),dimension(4) :: oooo_root, ooov_root, ovoo_root, oovv_root, vovo_root, &
                                       vvoo_root, vovv_root, vvvv_root, voov_root
         integer,dimension(2)       :: oo_id, ov_id, vo_id, vv_id 
         integer,dimension(4)       :: oooo_id, ooov_id, ovoo_id, oovv_id, vovo_id, &
                                       vvoo_id, vovv_id, vvvv_id, voov_id
         integer                    :: ao_space_id, occ_space_id, vir_space_id

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable    :: color(:)
         integer(INTL)                          :: labs

         ! types of the initializer methods
         type(compute_2e_ao_tensor_t) :: aoint_calculator ! method to initialize the ao-integrals tensor
         type(init_mocoef_tensor_t)   :: occupied_init(2), virtual_init(2)
         type(tens_printer_t)         :: tens_printer     ! generic method, but needs to be defined here in specific form
         type(scale_with_denom_t)     :: denominate

         ! variables for the definition of the initializer methods
         character(len=5)  :: initlabel
         character(len=2)  :: vs_str

         integer:: spin
         integer:: lsh, l, nl

         !CCD
         logical :: CCD, lambda

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) print*, " Starting CC calculation with exatensor"

!        Initialize global data (read from DFCOEF or similar by master)
         call set_talsh_only(exa_input%talsh)
         call initialize_global_data

         nao     = get_nao()   ! This is the total number of ao's
         nao_l   = int(nao,8)  ! we need this as a long int for exatensor

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
            nl   =  ncartesian(gto(lsh)%orb_momentum)
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

         !Register the ao vector space:
         ierr=exatns_space_register('AO_space',basis_ao,ao_space_id,ao_space)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')
         if (my_MPI_rank == 0) print*," registered AO space with nao, nshells:",nao,nshells

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         CCD = exa_input%ccd
         lambda = exa_input%lambda

         !Get the MO-coeffcient basis for the occupied and virtual spaces
         ierr = 0
         nocc_l   = 2*int(nocc,8)  ! NB: imo, jmo and nmo count orbitals, nocc_l counts spinors
         nvir_l   = 2*int(nvir,8)  ! NB: imo, jmo and nmo count orbitals, nvir_l counts spinors
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

         !Register this mo vector space :
         ierr=exatns_space_register('Occupied',basis_occ,occ_space_id,occ_space)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(occ) failed!')
         ierr=exatns_space_register('Virtual',basis_vir,vir_space_id,vir_space)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(vir) failed!')

         !Define separate initializers for the alpha(A) and beta(B) occupied and virtual MO-coefficients
         do spin = 1, 2
            !Define a unique name for the initializer for this spin (A/B)
            if (spin == 1) write (vs_str,'(A2)') '_A'
            if (spin == 2) write (vs_str,'(A2)') '_B'
            initlabel = 'occ' // vs_str
            !Register this initializer for occupied as a method
            call occupied_init(spin)%init_mocoef_ctor(mo_occ(1:nocc),nocc,spin)
            ierr=exatns_method_register(initlabel,occupied_init(spin))
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            !Same procedure for the virtual space
            initlabel = 'vir' // vs_str
            !Register this initializer for virtuals as a method
            call virtual_init(spin)%init_mocoef_ctor(mo_vir(1:nvir),nvir,spin)
            ierr=exatns_method_register(initlabel,virtual_init(spin))
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         end do

         if (my_MPI_rank == 0) print*," registered MO spaces and initializers with nocc and nvir:",nocc, nvir

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(aoint) failed!')

         !Register the denominator initializer aa a method
         call denominate%scale_with_denom_ctor(mo_occ,nocc,mo_vir,nvir)
         ierr=exatns_method_register('Denominate',denominate)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the standard print method with a low treshold to always print result_tensor
         call tens_printer%reset_thresh(1.0D-15)
         ierr=exatns_method_register('PrintTensor',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(printer) failed!')

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
         voov_id = (/vir_space_id, occ_space_id, occ_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)
         voov_root = (/vir_space_root, occ_space_root, occ_space_root, vir_space_root/)

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do.
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(global_communicator)

         if(ierr.eq.EXA_SUCCESS) then
            ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
            if(my_role.eq.EXA_DRIVER) then

            if (.not.CCD) then
                call print_date('starting ccsd calculation')
            else
                call print_date('starting ccd calculation')
            end if

!           Create and initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(result_tensor,"result_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(one_tensor,"another_one_tensor",EXA_DATA_KIND_C8)
            ! Note:  also get_integral_tensors will define a one_tensor and call it "one_tensor", so we needed another name
            ierr=exatns_tensor_init(one_tensor,ONE)

!           Get fixed one-body tensors => set to 0 for the moment (canonical orbitals)
            ierr=exatns_tensor_create(comm_t%foo,"foo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fov,"fov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fvv,"fvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            call print_date('initialized fock matrix')

!           Get fixed two-body tensors
            ierr=exatns_tensor_create(int_t%oooo,"oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%ooov,"ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%oovv,"oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovo,"vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovv,"vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vvvv,"vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(voov_tensor,"voov",voov_id,voov_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
            call get_CC_integrals (ao_space_id,ao_space_root,voov_tensor,int_t)
            ierr=exatns_tensor_destroy(voov_tensor)

            call print_date('initialized two electron integrals')

!           Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
            if (.not.CCD) ierr=exatns_tensor_create(comm_t%t1,"t1",vo_id,vo_root,EXA_DATA_KIND_C8)

!           Initialize t2 amplitudes with the MP2 values
            ierr=exatns_tensor_create(comm_t%t2,"t2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
!           the tensor_add routine is not implemented in exatensor, work around this by using a contraction
            ierr=exatns_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)*K()",comm_t%t2,int_t%oovv,one_tensor) ! o(2)v(2)
            ierr=exatns_tensor_transform(comm_t%t2,'Denominate')

!           Calculate and print the MP2 energy
            ierr=exatns_tensor_init(result_tensor)
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv,ONE_QUARTER) ! o(2)v(2)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=real(result_val,8) ! Casting into real(8)
            print*," MP2 energy = ", mp2_energy

!           Start CC(S)D iterations

            ierr=exatns_tensor_create(comm_t%goo,"goo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%gvv,"gvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hoo,"hoo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hvv,"hvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%hov,"hov",ov_id,ov_root,EXA_DATA_KIND_C8)
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

            call solve_t_equations(comm_t,int_t,dims,CCD,t_target_precision,ccsd_energy)

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

            if (lambda) then
             print*, ""
             call print_date("Begin solving Lambda equations")

!            Construct lambda tensors
             ierr=exatns_tensor_create(l1_tensor,"l1",ov_id,ov_root,EXA_DATA_KIND_C8)
             ierr=exatns_tensor_create(l2_tensor,"l2",oovv_id,oovv_root,EXA_DATA_KIND_C8)

!            Initialize l1 and l2 with converged t amplitudes
             ierr=exatns_tensor_contract("L(i,a)+=T+(a,i)*K()",l1_tensor,comm_t%t1,one_tensor)
             ierr=exatns_tensor_contract("L(i,j,a,b)+=T+(a,b,i,j)*K()",l2_tensor,comm_t%t2,one_tensor)

             call solve_lambda_equations(comm_t,int_t,l1_tensor,l2_tensor,dims,CCD,t_target_precision)            

             if (.not.CCD) call get_dm_mo(comm_t%t1,comm_t%t2,l1_tensor,l2_tensor,dims)
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
            if (.not.CCD) ierr=exatns_tensor_destroy(comm_t%t1)
            ierr=exatns_tensor_destroy(comm_t%t2)
            ierr=exatns_tensor_destroy(comm_t%tau)
            ierr=exatns_tensor_destroy(comm_t%hoo)
            ierr=exatns_tensor_destroy(comm_t%hov)
            ierr=exatns_tensor_destroy(comm_t%hvv)
            ierr=exatns_tensor_destroy(comm_t%goo)
            ierr=exatns_tensor_destroy(comm_t%gvv)
            ierr=exatns_tensor_destroy(comm_t%a_int)
            ierr=exatns_tensor_destroy(comm_t%h_int)
            ierr=exatns_tensor_destroy(result_tensor)
            ierr=exatns_tensor_destroy(one_tensor)
            ierr=exatns_tensor_destroy(l1_tensor)
            ierr=exatns_tensor_destroy(l2_tensor)

!           Get statistics 
            !ierr=exatns_status()    !Need to find correct exatns alternatives. Probably exatns_status

            !Stop ExaTENSOR runtime
            ierr=exatns_stop()

            end if

         else
            write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

         deallocate(mo_occ)
         deallocate(mo_vir)

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         if (my_MPI_rank == 0) print*, " Leaving cc_driver routine"
         return

        end subroutine exacorr_cc_driver

        subroutine solve_t_equations(comm_t,int_t,dims,CCD,t_target_precision,ccsd_energy)

!        Routine for solving the t amplitude equations

         use exacorr_ao_to_mo, only : print_tensornorm2

!        Input
         type(exatns_comm_tens), intent(inout) :: comm_t
         type(exatns_intg_tens), intent(inout) :: int_t
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
         integer,dimension(2) :: oo_id,ov_id,vo_id
         integer,dimension(4) :: ovoo_id,vvoo_id
         integer(INTL),dimension(2) :: oo_root,ov_root,vo_root
         integer(INTL),dimension(4) :: ovoo_root,vvoo_root

         real(8) :: t_convergence

         integer :: iccsd, nccsd=30
         integer :: ierr
         
         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         oo_id     = (/occ_space_id, occ_space_id/)
         ov_id     = (/occ_space_id, vir_space_id/)
         vo_id     = (/vir_space_id, occ_space_id/)
         oo_root   = (/occ_space_root, occ_space_root/)
         ov_root   = (/occ_space_root, vir_space_root/)
         vo_root   = (/vir_space_root, occ_space_root/)
         ovoo_id   = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vvoo_id   = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

!        Initialize scalar that are to be used as tensors in contractions
         ierr=exatns_tensor_create(result_t_tensor,"result_t",EXA_DATA_KIND_C8)

!        Construct intermediate tensors
         ierr=exatns_tensor_create(hoo_aux_tensor,"hoo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(hov_aux_tensor,"hov_aux",ov_id,ov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(s1_tensor,"s1",vo_id,vo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(s2_tensor,"s2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(tau2_tensor,"tau2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(cint_tensor,"cint",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)

         do iccsd = 1, nccsd

          call get_tau(comm_t,CCD,one_tensor)

!         NOTE: contractions marked with "validated" gave numerically significant contributions and could be validated by comparison
!         to the reference implementation, the ones marked "moretests" were too small to decide whether it was correct or not.
!         calculate intermediaries

          ierr=exatns_tensor_init(s2_tensor)

!         The vvvv contribution
          ierr=exatns_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,int_t%vvvv,comm_t%tau,ONE_HALF) !validated o(2)v(4)

!         vv ! validated
          call get_hvv(comm_t,int_t%oovv,one_tensor)
          call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)

!         vovo ! moretests
          ierr=exatns_tensor_init(tau2_tensor)
          ierr=exatns_tensor_contract("T(a,d,l,i)+=R(a,d,l,i)*K()",tau2_tensor,comm_t%t2,one_tensor,ONE_HALF)
          if (.not.CCD) ierr=exatns_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,comm_t%t1,comm_t%t1) ! o(2)v(2)
          call get_h(comm_t,int_t,CCD,one_tensor)
          ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",comm_t%h_int,int_t%oovv,tau2_tensor)! o(3)v(3)

!         This term appears four times as we need to permute both i and j and a and b
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,comm_t%h_int,comm_t%t2)  ! validated o(3)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,comm_t%h_int,comm_t%t2)  ! validated o(3)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,comm_t%h_int,comm_t%t2)  ! validated o(3)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,comm_t%h_int,comm_t%t2)  ! validated o(3)v(3)

!         ovoo ! validated
          ierr=exatns_tensor_init(cint_tensor)
          ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,int_t%vovv,comm_t%tau,MINUS_ONE_HALF) ! o(3)v(3)
          if (.not.CCD) then
              ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,int_t%vovo,comm_t%t1) ! o(3)v(2)
              ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,int_t%vovo,comm_t%t1,MINUS_ONE) ! o(3)v(2)
          end if
          ierr=exatns_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,int_t%ooov,one_tensor) ! o(3)v(1)

!         oo ! validated
          call get_hoo(comm_t,int_t%oovv,one_tensor)
          call get_goo(comm_t,int_t%ooov,CCD,one_tensor)

!         ov ! validated
          call get_hov(comm_t,int_t%oovv,CCD,one_tensor)
          ierr=exatns_tensor_init(hov_aux_tensor)
          ierr=exatns_tensor_contract("L(k,c)+=H(k,c)*K()",hov_aux_tensor,comm_t%hov,one_tensor)        ! copy common term
          ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hov_aux_tensor,comm_t%fov,one_tensor,MINUS_ONE) ! subtract fock matrix

!         oooo ! validated
          call get_a(comm_t,int_t,CCD,one_tensor)

          if (ierr.ne.EXA_SUCCESS) then
             print*," ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

          if (.false.) then
             if (.not.CCD) write(*,'(A7,F20.15)'), "T1amp = ",print_tensornorm2 (comm_t%t1)
             write(*,'(A7,F20.15)'), "T2amp = ",print_tensornorm2 (comm_t%t2)
             write(*,'(A7,F20.15)'), "Tau   = ",print_tensornorm2 (comm_t%tau)
             write(*,'(A7,F20.15)'), "hoint = ",print_tensornorm2 (comm_t%hoo)
             write(*,'(A7,F20.15)'), "goint = ",print_tensornorm2 (comm_t%goo)
             write(*,'(A7,F20.15)'), "hvint = ",print_tensornorm2 (comm_t%hvv)
             write(*,'(A7,F20.15)'), "gvint = ",print_tensornorm2 (comm_t%gvv)
             write(*,'(A7,F20.15)'), "hovint= ",print_tensornorm2 (comm_t%hov)
             write(*,'(A7,F20.15)'), "aint  = ",print_tensornorm2 (comm_t%a_int)
             write(*,'(A7,F20.15)'), "cint  = ",print_tensornorm2 (cint_tensor)
             write(*,'(A7,F20.15)'), "hint  = ",print_tensornorm2 (comm_t%h_int)
          end if

!         This term appears twice as we need to permute i and j
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,comm_t%goo,comm_t%t2)                   ! validated o(3)v(2)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,comm_t%goo,comm_t%t2)                   ! validated o(3)v(2)
!         Again two times, also note the complex conjugation as we use vovv instead of vvvo
          if (.not.CCD) then
              ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,int_t%vovv,comm_t%t1)                   ! validated o(2)v(3)
              ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,int_t%vovv,comm_t%t1,MINUS_ONE)         ! validated o(2)v(3)
          end if

!         This term appears twice as we need to permute a and b
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! validated o(2)v(3)
          ierr=exatns_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! validated o(2)v(3)
!         Again two times, with the C intermediate (note: in RELCCSD this is done differently)
          if (.not.CCD) then
              ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,comm_t%t1,MINUS_ONE) ! o(3)v(2)
              ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,comm_t%t1)           ! o(3)v(2)
          end if

!         The oooo contribution
          ierr=exatns_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,comm_t%a_int,comm_t%tau,ONE_HALF) !validated o(4)v(2)

          ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)    ! validated o(2)v(2)

!         calculate s1 tensor
          ierr=exatns_tensor_init(s1_tensor)
          ierr=exatns_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
          ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,ONE_HALF) !validated o(2)v(3)
          if (.not.CCD) then    
              ierr=exatns_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                  !validated o(1)v(2)
              ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                        !validated o(2)v(1)
          end if
          ierr=exatns_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)                   !validated o(2)v(2)
          ierr=exatns_tensor_init(hoo_aux_tensor)                                                                        !moretests
          if (.not.CCD) then    
              ierr=exatns_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)                    !moretests o(2)v(1)
              ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)                        !moretests o(2)v(1)
          end if
          ierr=exatns_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,ONE_HALF) !validated o(3)v(2)
          if (.not.CCD) ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,MINUS_ONE)     !validated o(2)v(2)

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
          end if

!         Update tau with converged values
          call get_tau(comm_t,CCD,one_tensor)

          ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_t_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ierr=exatns_tensor_get_scalar(result_t_tensor,result_val)
          ccsd_energy = ccsd_energy + real(result_val*ONE_QUARTER,8)
          if (.not.CCD) then
              write (*,'(A,I4,A,F20.15,A,ES10.1)') "  Iteration ",iccsd," CCSD energy = ", &
                                                  ccsd_energy, " convergence =",t_convergence
          else
              write (*,'(A,I4,A,F20.15,A,ES10.1)') "  Iteration ",iccsd," CCD energy = ", &
                                                  ccsd_energy, " convergence =",t_convergence
          end if
          if (t_convergence .lt. t_target_precision) exit
         end do

         ierr=exatns_tensor_destroy(s1_tensor)
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(tau2_tensor)
         ierr=exatns_tensor_destroy(cint_tensor)
         ierr=exatns_tensor_destroy(hoo_aux_tensor)
         ierr=exatns_tensor_destroy(hov_aux_tensor)
         ierr=exatns_tensor_destroy(result_t_tensor)

        end subroutine solve_t_equations

        subroutine solve_lambda_equations(comm_t,int_t,l1_tensor,l2_tensor,dims,CCD,t_target_precision)

!        Routine for solving the lambda equations

         use exacorr_ao_to_mo, only       : print_tensornorm2

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_comm_tens), intent(inout) :: comm_t
         type(exatns_intg_tens), intent(inout) :: int_t
         type(exatns_lambda_tens)              :: lambda_t
!        Mixed type intermediates
         type(tens_rcrsv_t) :: goo_tensor, gvv_tensor
!        solution tensors
         type(tens_rcrsv_t) :: l1_tensor, l2_tensor
         type(tens_rcrsv_t) :: s1_tensor, s2_tensor
!        Aux tensor
         type(tens_rcrsv_t) :: oo_aux_tensor
!        Tensor dimensions
         type(space_dims),intent(in) :: dims
         integer :: occ_space_id, vir_space_id
         integer(INTL) :: occ_space_root, vir_space_root
         integer, dimension(2) :: oo_id, ov_id, vv_id
         integer, dimension(4) :: oooo_id,ovoo_id,ooov_id,vovo_id,oovv_id, &
             vvvo_id,vovv_id,vvvv_id
         integer(INTL), dimension(2) :: oo_root, ov_root, vv_root
         integer(INTL), dimension(4) :: oooo_root,ovoo_root,ooov_root,vovo_root,oovv_root, &
                   vvvo_root,vovv_root,vvvv_root
         integer :: iccsd, nccsd=30
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence

         logical, intent(in) :: CCD
         integer :: ierr, prints=0

         call print_date('Entered lambda solver routine')

!        Create all ids and roots
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

!        Construct fixed intermediates
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

!        Construct auxiliary tensor
         ierr=exatns_tensor_create(oo_aux_tensor,"oo_aux",oo_id,oo_root,EXA_DATA_KIND_C8)

!        Obtain intermediates from CCSD
         call print_date("Calculating common blocks")
         call get_tau(comm_t,CCD,one_tensor)
         call get_hvv(comm_t,int_t%oovv,one_tensor)
         call get_gvv(comm_t,int_t%vovv,CCD,one_tensor)
         call get_h(comm_t,int_t,CCD,one_tensor)
         call get_hoo(comm_t,int_t%oovv,one_tensor)
         call get_goo(comm_t,int_t%ooov,CCD,one_tensor)
         call get_hov(comm_t,int_t%oovv,CCD,one_tensor)
         call get_a(comm_t,int_t,CCD,one_tensor)
         call print_date("Finished calculating common blocks")

!        Get fixed intermediates
         call print_date('Calculating fixed intermediates')
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

         !Debug
         print*, "T1amp      = ", print_tensornorm2(comm_t%t1)
         print*, "T2amp      = ", print_tensornorm2(comm_t%t2)*0.25D0
!        print*, "Tau        = ", print_tensornorm2(comm_t%tau)*0.25D0
!        print*, "Fbar_oo    = ", print_tensornorm2(lambda_t%fbar_oo)*0.5D0
         print*, "Fbar_ov    = ", print_tensornorm2(lambda_t%fbar_ov)
!        print*, "Fbar_vv    = ", print_tensornorm2(lambda_t%fbar_vv)*0.5D0
!        print*, "w_vvvv     = ", print_tensornorm2(lambda_t%w_vvvv)*0.25D0
         print*, "wbar_vovo  = ", print_tensornorm2(lambda_t%wbar_vovo)
         print*, "w_vvvo     = ", print_tensornorm2(lambda_t%w_vvvo)*0.5D0
         print*, "w_vovv     = ", print_tensornorm2(lambda_t%w_vovv)*0.5D0
         print*, "w_vovo     = ", print_tensornorm2(lambda_t%w_vovo)
         print*, "w_ovoo     = ", print_tensornorm2(lambda_t%w_ovoo)*0.5D0
         print*, "w_ooov     = ", print_tensornorm2(lambda_t%w_ooov)*0.5D0
         print*, "w_oooo     = ", print_tensornorm2(lambda_t%w_oooo)*0.25D0

!        Construct solution tensors
         ierr=exatns_tensor_create(s2_tensor,"s2",oovv_id,oovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(s1_tensor,"s1",ov_id,ov_root,EXA_DATA_KIND_C8)

!        Construct mixed type intermediates
         ierr=exatns_tensor_create(goo_tensor,"goo_mixed",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(gvv_tensor,"gvv_mixed",vv_id,vv_root,EXA_DATA_KIND_C8)

!        Reinitialize diis, so it can be used in lambda routine
         call diis (l2_tensor,oovv_id,oovv_root,l1_tensor,ov_id,ov_root,mode=1)

         call print_date('Begin iterating')
         do iccsd = 1, nccsd

!         Form mixed intermediates
          ierr=exatns_tensor_init(gvv_tensor)
          ierr=exatns_tensor_contract("G(e,a)+=L(m,n,a,f)*T(e,f,m,n)",gvv_tensor,l2_tensor,comm_t%t2,MINUS_ONE_HALF)
          ierr=exatns_tensor_init(goo_tensor)
          ierr=exatns_tensor_contract("G(i,m)+=L(i,n,e,f)*T(e,f,m,n)",goo_tensor,l2_tensor,comm_t%t2,ONE_HALF)

          if (prints .ge. 1) then
            print*, "L1  = ", print_tensornorm2(l1_tensor)
            print*, "L2  = ", print_tensornorm2(l2_tensor)*0.25D0
            print*, "goo = ", print_tensornorm2(goo_tensor)
            print*, "gvv = ", print_tensornorm2(gvv_tensor)
          end if

!         Calculate s1 tensor
          ierr=exatns_tensor_init(s1_tensor)
          ierr=exatns_tensor_contract("S(i,a)+=F(i,a)*K()",s1_tensor,lambda_t%fbar_ov,one_tensor)
          ierr=exatns_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)
          ierr=exatns_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,goo_tensor,lambda_t%w_ooov)
          ierr=exatns_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,gvv_tensor,lambda_t%w_vovv)

!         Calculate s2 tensor
          ierr=exatns_tensor_init(s2_tensor)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)

!         Need to permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)

!         Need to permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,ONE_HALF)

!         Need to permute i and j, as well as a and b, so term appears four times
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!         Need to permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,gvv_tensor)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,gvv_tensor)

!         Need to permute a and b, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)

!         Need to permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,goo_tensor)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,goo_tensor)

!         Need to permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)

!         Need to permute i and j, as well as a and b, so term appears four times
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)

          ierr=exatns_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,ONE_HALF)


          if (ierr.ne.EXA_SUCCESS) then
            print*," ccdriver: error in evaluation of lambda equations, code=",ierr
            return
          end if

!         prepare for next iteration
          ierr=exatns_tensor_init(l1_tensor)
          ierr=exatns_tensor_contract("L1(i,a)+=S1(i,a)*K()",l1_tensor,s1_tensor,one_tensor)
          !print*, "L1 before = ", print_tensornorm2(l1_tensor)
          ierr=exatns_tensor_transform(l1_tensor,"Denominate")
          !print*, "L1 after = ", print_tensornorm2(l1_tensor)
          ierr=exatns_tensor_init(l2_tensor)
          ierr=exatns_tensor_contract("L2(i,j,a,b)+=S2(i,j,a,b)*K()",l2_tensor,s2_tensor,one_tensor)
          !print*, "L2 before = ", print_tensornorm2(l2_tensor)*ONE_QUARTER
          ierr=exatns_tensor_transform(l2_tensor,"Denominate")
          !print*, "L2 after = ", print_tensornorm2(l2_tensor)*ONE_QUARTER

          call diis (l2_tensor,t1_tensor=l1_tensor,t_convergence=t_convergence)

          write (*,'(A,I4,A,ES10.1)') "  Iteration ",iccsd," convergence =",t_convergence
          if (t_convergence .lt. t_target_precision) exit

         end do

         call print_date("Finished iterating")

         if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
         write(*,'(A,ES22.16)'),"  Final L1amp = ", print_tensornorm2(l1_tensor)
         write(*,'(A,ES22.16)'),"  Final L2amp = ", print_tensornorm2(l2_tensor)*ONE_QUARTER

!        Destroy diis tensors
         call diis(l2_tensor,t1_tensor=l1_tensor,mode=2)
         ierr=exatns_tensor_destroy(s1_tensor)
         ierr=exatns_tensor_destroy(s2_tensor)
         ierr=exatns_tensor_destroy(goo_tensor)
         ierr=exatns_tensor_destroy(gvv_tensor)
         ierr=exatns_tensor_destroy(oo_aux_tensor)
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

         call print_date('Leaving lambda solver routine')

        end subroutine solve_lambda_equations

        subroutine get_CC_integrals (ao_space_id,ao_space_root,voov_tensor,int_t)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017

         use exacorr_ao_to_mo
         use exacorr_global

!        information needed to compute the ao-tensor
         integer,intent(in)       :: ao_space_id
         integer(INTL),intent(in) :: ao_space_root
!        target integral tensors to be filled and given back to the caller
         type(exatns_intg_tens), intent(inout) :: int_t 
         type(tens_rcrsv_t),intent(inout) :: voov_tensor ! Auxilliary tensor that is not anti-symmetrized

         integer                           :: i,ierr

         ! variables for the definition of the initializer methods
         character(len=5)  :: initlabel(2,4),initocc(2),initvir(2)

!        auxilliary integral tensor (only needed inside this routine)
         type(tens_rcrsv_t) :: aoint_tensor

!        scalars (need to be defined as tensor types)

         initocc(1) = 'occ_A'
         initocc(2) = 'occ_B'
         initvir(1) = 'vir_A'
         initvir(2) = 'vir_B'

!        Get tensor with AO integrals
         ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

         !Initialize ao-tensor:
         ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
         if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
         call print_date('finished ao integral calculation')

!        Note that the spaces for the tensor are given in physicists notation <12||34>
!        whereas the initializer labels needed in the index transform are in mullikens notation (13||24)

!        Get oovv tensor
         initlabel(:,1)=initocc
         initlabel(:,2)=initvir
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,int_t%oovv,12,initlabel)

!        Get oooo tensor
         initlabel(:,1)=initocc
         initlabel(:,2)=initocc
         initlabel(:,3)=initocc
         initlabel(:,4)=initocc
         call get_integral_tensor (aoint_tensor,int_t%oooo,12,initlabel)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initvir
         initlabel(:,3)=initocc
         initlabel(:,4)=initocc
         call get_integral_tensor (aoint_tensor,int_t%vovo,0,initlabel)
         initlabel(:,1)=initvir
         initlabel(:,2)=initocc
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,voov_tensor,0,initlabel)
!        Make anti-symmetrized vovo tensor
         ierr=exatns_tensor_contract("V(a,i,b,j)+=V(a,i,j,b)*K()",int_t%vovo,voov_tensor,one_tensor,MINUS_ONE)

!        Get vvvv tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initvir
         initlabel(:,3)=initvir
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,int_t%vvvv,12,initlabel)

!        Get ooov tensor
         initlabel(:,1)=initocc
         initlabel(:,2)=initocc
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,int_t%ooov,12,initlabel)

!        Get vovv tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initvir
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,int_t%vovv,34,initlabel)

!        Clean-up
         ierr=exatns_tensor_destroy(aoint_tensor)

        end subroutine get_CC_integrals

        subroutine diis (t2_tensor,vvoo_id,vvoo_root,t1_tensor,vo_id,vo_root,t_convergence,mode)

         integer, parameter             :: maxdim = 5 ! maximum number of tensors that will be stored
         integer, parameter             :: mindim = 3 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(tens_rcrsv_t), intent(inout),optional :: t1_tensor
         type(tens_rcrsv_t), intent(inout)          :: t2_tensor
         integer                                    :: ierr, i, tens_rank, j
         integer,dimension(2),optional              :: vo_id
         integer,dimension(4),optional              :: vvoo_id
         integer(INTL),dimension(2),optional        :: vo_root
         integer(INTL),dimension(4),optional        :: vvoo_root
         logical                                    :: init, clean_up, CCD
         integer, save                              :: itens, jtens, ntens
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
