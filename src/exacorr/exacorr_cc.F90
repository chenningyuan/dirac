module exacorr_cc
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra  ! still necessary for type definitions, like C_SIZE_T, etc
        use exatensor

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0),ONE_QUARTER=(0.25D0,0.D0)!, &
!                                THREE=(3.0,0.0), SIX=(6.0,0.0)
        private

        public cc_driver
        
       contains
!-------------------------------------------------------------------------------------------
        subroutine cc_driver (nocc,nvir,mo_occ,mo_vir)

!        This subroutine drives the ccsd energy calculation
!        Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
 
!        Written by Lucas Visscher, winter 2016/2017
!        ExaTensor conversion by Anastasios Papadopoulos, October 2018

!        COMMENTED OUT (IN MY OPINION) UNNECESSARY INIT STATEMENTS
         
         use interface_to_mpi
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use, intrinsic:: ISO_C_BINDING

         implicit none

         integer, intent(in) :: nocc,nvir           ! the size of the mo basis for occupied and virtual spinors
         integer(INTL)       :: nocc_l, nvir_l
         integer, intent(in) :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer :: my_role,mpi_th_provided, my_MPI_rank    ! needed to start exatensor, MPI should already be started in threaded mode
         integer                                :: nshells  ! the number of shells (orbitals with the same center and l-value)
         integer                                :: nao      ! the number of atomic orbitals)
         type(basis_func_info_t), allocatable   :: gto(:)   ! array with information about the shells
         
         integer                        :: ierr, i
         complex(8), pointer            :: energy(:),result_tens(:)
         integer(C_SIZE_T), parameter   :: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_SIZE_T)              :: host_buf_size
         integer(C_INT)                 :: host_arg_max

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(tens_rcrsv_t) :: foo_tensor, fvv_tensor, fov_tensor
         type(tens_rcrsv_t) :: oooo_tensor, ooov_tensor, oovv_tensor, vovo_tensor, vovv_tensor, vvvv_tensor
!        intermediate tensors (same notation as in relccsd)
         type(tens_rcrsv_t) :: goint_tensor, gvint_tensor, hoint_tensor, hvint_tensor, hovint_tensor
         type(tens_rcrsv_t) :: aint_tensor, bint_tensor, cint_tensor, hint_tensor
         type(tens_rcrsv_t) :: tau_tensor, tau2_tensor
!        auxilliary tensors for t1 equation
         type(tens_rcrsv_t) :: hovaux_tensor, hoaux_tensor
!        solution tensors
         type(tens_rcrsv_t) :: t1_tensor, t2_tensor, s1_tensor, s2_tensor
!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t) :: result_tensor, one_tensor
!        CCSD control and result variables
         integer            :: iccsd, nccsd
         real(8)            :: mp2_energy, ccsd_energy!, t_energy(3), t_convergence  ! triples not implemented
!        real(8), parameter :: t_target_precision = 1.D-13

         ! classes for the definition of spaces (AO and 4 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         class(h_space_t), pointer                :: occ_space, vir_space

         ! variables for the definition of spaces (AO and 4 different MO spaces)
         type(subspace_basis_t) :: basis_ao, basis_occ, basis_vir
         integer(INTL)          :: nao_l, ao_space_root
         integer(INTL)          :: occ_space_root, vir_space_root
         integer                :: ao_space_id, occ_space_id, vir_space_id

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable    :: color(:)
         integer(INTL)                          :: labs
         integer                                :: l_col=1

         ! types of the initializer methods
         type(compute_2e_ao_tensor_t) :: aoint_calculator ! method to initialize the ao-integrals tensor
         type(init_mocoef_tensor_t)   :: occupied_init(2,nocc), virtual_init(2,nvir)
         type(tens_printer_t)         :: tens_printer     ! generic method, but needs to be defined here in specific form
         type(scale_with_denom_t)     :: denominate

         ! variables for the definition of the initializer methods
         character(len=5)  :: initlabel(2)
         character(len=2)  :: vs_str

         integer:: spin
         integer:: lsh, l, nl

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) print*, " Starting CC calculation with exatensor"

!        Initialize global data (read from DFCOEF or similar by master)
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
!           call color(lsh)%color_symmetry_ctor(l_col) ! give them all the same color
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
         if (my_MPI_rank == 0) print*," registered AO vector space with nao, nshells:",nao,nshells
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Get the MO-coeffcient basis for the occupied and virtual spaces
         ierr = 0
         nocc_l   = 2*int(nocc,8)  ! NB: imo, jmo and nmo count orbitals, nocc_l counts spinors
         nvir_l   = 2*int(nvir,8)  ! NB: imo, jmo and nmo count orbitals, nvir_l counts spinors
         call basis_occ%subspace_basis_ctor(nocc_l,ierr)
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
         call basis_vir%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

         !Register this mo vector space :
         write (vs_str,'(A,I1)') '_',i
         !ierr=exatns_space_register('Occupied',basis_occ(i),occ_space_id,occ_space%ptr) => LUUK WROTE BASIS_OCC(I) => WHY?
         !ierr=exatns_space_register('Virtual',basis_cir(i),vir_space_id,vir_space%ptr)
         ierr=exatns_space_register('Occupied',basis_occ,occ_space_id,occ_space)
         ierr=exatns_space_register('Virtual',basis_vir,vir_space_id,vir_space)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Might be faster if spin is outer loop => ASK
         !Define separate initializers for the alpha(A) and beta(B) occupied MO-coefficients for each index 
         do i = 1, nocc
            do spin = 1, 2
               !Define a unique name for the initializer for this spin (A/B) of space i
               if (spin == 1) write (vs_str,'(A2)') '_A'
               if (spin == 2) write (vs_str,'(A2)') '_B'
               initlabel(spin) = 'occ' // vs_str
               !Register this initializer as a method
               call occupied_init(spin,i)%init_mocoef_ctor(mo_occ(1:nocc),nocc,spin)
               ierr=exatns_method_register(initlabel(spin),occupied_init(spin,i))
            end do
         end do
 
         do i = 1, nvir   
            do spin = 1, 2
               initlabel(spin) = 'vir' // vs_str
               !Register this initializer as a method
               call virtual_init(spin,i)%init_mocoef_ctor(mo_vir(1:nvir),nvir,spin)
               ierr=exatns_method_register(initlabel(spin),virtual_init(spin,i))
               if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            end do         
         end do

         if (my_MPI_rank == 0) print*," registered MO vectors spaces and initializers with nocc and nvir:",nocc, nvir

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the denominator initializer aa a method (testing only, registration will move to exacorr_cc)
         !call denominate%scale_with_denom_ctor(mo_list(1:2),2,mo_list(3:nmo(1)),nmo(1)-2)
         !ierr=exatns_method_register('Denominate',denominate)
         !if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
   
         !Register the standard print method with a high treshold to print only important tensor elements
         call tens_printer%reset_thresh(0.01D0)
         ierr=exatns_method_register('PrintTensor',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do. 
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(global_communicator)

         if(ierr.eq.EXA_SUCCESS) then
            ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
            if(my_role.eq.EXA_DRIVER) then

            !call print_date('starting ao integral calculation')

!           Get fixed one-body tensors => set to 0 
            call get_CC_Fockmatrix (occ_space,vir_space,occ_space_id,vir_space_id,foo_tensor,fov_tensor,fvv_tensor)

!           Get fixed two-body tensors
!           call get_CC_integrals (nocc,nvir,mo_occ,mo_vir, &
!                               oooo_tensor,ooov_tensor,oovv_tensor, &
!                               vovo_tensor,vovv_tensor,vvvv_tensor)

!           Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
            !t1_dims(1) = nvir * 2
            !t1_dims(2) = nocc * 2
            !ierr=exatns_tensor_create(t1_tensor,"t1_tensor",(),(),EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(t1_tensor,ZERO)

!           Create and initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(one_tensor,"one_tensor",EXA_DATA_KIND_C8)   
            ierr=exatns_tensor_create(result_tensor,"result_tensor",EXA_DATA_KIND_C8)   
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         
            ierr=exatns_tensor_init(one_tensor,ONE)
            ierr=exatns_tensor_init(result_tensor,ZERO)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

!           Initialize t2 amplitudes with the MP2 values
            !t2_dims(1:2) = nvir * 2
            !t2_dims(3:4) = nocc * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(t2_tensor,ZERO)
!           ierr=exatns_tensor_add("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor)
!           work-around because complex conjugation does not work in the above statement
            ierr=exatns_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor,one_tensor) ! o(2)v(2)
!           end of work-around

!           Calculate and print the MP2 energy
            !ierr=exatns_tensor_init(result_tensor) => already initialized. Why is it here again in the original implementation?
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor) ! o(2)v(2)
            mp2_energy=result_tens(1) * ONE_QUARTER
            print*," MP2 energy = ", mp2_energy

!           Start CCSD iterations

            nccsd = 30 ! hard wired for the moment

            !goint_dims = nocc  * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(goint_tensor,ZERO)
            !gvint_dims = nvir * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(gvint_tensor,ZERO)
            !hoint_dims = nocc  * 2
            !ierr=exatns_tensor_create()
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(hoint_tensor,ZERO)
            ierr=exatns_tensor_init(hoaux_tensor,ZERO)
            !hvint_dims = nvir * 2!
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(hvint_tensor,ZERO)
            !hovint_dims(1) = nocc * 2
            !hovint_dims(2) = nvir * 2
            !ierr=exatns_tensor_create()
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(hovint_tensor,ZERO)
            ierr=exatns_tensor_init(hovaux_tensor,ZERO)
            !aint_dims  = nocc  * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(aint_tensor,ZERO)
            !cint_dims(1) = nocc * 2
            !cint_dims(2) = nvir * 2
            !cint_dims(3) = nocc * 2
            !cint_dims(4) = nocc * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(cint_tensor,ZERO)
            !hint_dims(1) = nvir * 2
            !hint_dims(2) = nocc * 2
            !hint_dims(3) = nvir * 2
            !hint_dims(4) = nocc * 2
            !ierr=exatns_tensor_create()
            !ierr=exatns_tensor_create()
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(hint_tensor,ZERO)
            ierr=exatns_tensor_init(tau_tensor,ZERO)
            ierr=exatns_tensor_init(tau2_tensor,ZERO)
            !s1_dims(1) = nvir * 2
            !s1_dims(2) = nocc * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(s1_tensor,ZERO)
            !s2_dims(1:2) = nvir * 2
            !s2_dims(3:4) = nocc * 2
            !ierr=exatns_tensor_create()
            ierr=exatns_tensor_init(s2_tensor,ZERO)
!           call diis (t1_tensor,t2_tensor,mode=1)
            if (ierr.eq.EXA_SUCCESS) then
               print*," ccdriver: allocated all tensors"
            else
               print*," ccdriver: error in tensor allocation, code=",ierr
               return
            end if

            do iccsd = 1, nccsd

!           make tau ! validated
            !ierr=exatns_tensor_init(tau_tensor)
            ierr=exatns_tensor_add("T(a,b,i,j)+=R(a,b,i,j)",tau_tensor,t2_tensor)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",tau_tensor,t1_tensor,t1_tensor) ! o(2)v(2)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",tau_tensor,t1_tensor,t1_tensor,MINUS_ONE) ! o(2)v(2)

!           NOTE: contractions marked with "validated" gave numerically significant contributions and could be validated by comparison
!           to the reference implementation, the ones marked "moretests" were too small to decide whether it was correct or not.
!           calculate intermediaries

            !ierr=exatns_tensor_init(s2_tensor)

!           The vvvv contribution
            ierr=exatns_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,vvvv_tensor,tau_tensor,ONE_HALF) !validated o(2)v(4)

!           vv ! validated
            !ierr=exatns_tensor_init(hvint_tensor)
            ierr=exatns_tensor_add("H(a,c)+=F(a,c)",hvint_tensor,fvv_tensor)
            ierr=exatns_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",hvint_tensor,oovv_tensor,tau_tensor,ONE_HALF) ! o(2)v(3)
            !ierr=exatns_tensor_init(gvint_tensor)
            ierr=exatns_tensor_add("G(a,c)+=H(a,c)",gvint_tensor,hvint_tensor)
            ierr=exatns_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",gvint_tensor,vovv_tensor,t1_tensor) ! o(1)v(3)
            ierr=exatns_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",gvint_tensor,fov_tensor,t1_tensor,MINUS_ONE) ! o(1)v(2)

!           vovo ! moretests
            !ierr=exatns_tensor_init(tau2_tensor)
            ierr=exatns_tensor_add("T(a,d,l,i)+=R(a,d,l,i)",tau2_tensor,t2_tensor,ONE_HALF)
            ierr=exatns_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,t1_tensor,t1_tensor) ! o(2)v(2)
            !ierr=exatns_tensor_init(hint_tensor)
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",hint_tensor,oovv_tensor,tau2_tensor)! o(3)v(3)
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",hint_tensor,vovv_tensor,t1_tensor) ! o(2)v(3)
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",hint_tensor,ooov_tensor,t1_tensor,MINUS_ONE) ! o(3)v(2)
            ierr=exatns_tensor_add("H(a,k,c,i)+=V(a,k,c,i)",hint_tensor,vovo_tensor)

!           This term appears four times as we need to permute both i and j and a and b
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)

!           ovoo ! validated
            !ierr=exatns_tensor_init(cint_tensor)
            ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,vovv_tensor,tau_tensor,MINUS_ONE_HALF) ! o(3)v(3)
            ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,vovo_tensor,t1_tensor) ! o(3)v(2)
            ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,vovo_tensor,t1_tensor,MINUS_ONE) ! o(3)v(2)
            ierr=exatns_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,ooov_tensor,one_tensor) ! o(3)v(1)

!           oo ! validated
            !ierr=exatns_tensor_init(hoint_tensor)
            ierr=exatns_tensor_add("H(k,j)+=F(k,j)",hoint_tensor,foo_tensor)
            ierr=exatns_tensor_contract("H(k,j)+=V(k,l,c,d)*T(c,d,l,j)",hoint_tensor,oovv_tensor,tau_tensor,ONE_HALF) ! o(3)v(2)
            !ierr=exatns_tensor_init(goint_tensor)
            ierr=exatns_tensor_add("G(k,j)+=H(k,j)",goint_tensor,hoint_tensor)
            ierr=exatns_tensor_contract("G(k,j)+=V(k,l,j,c)*T(c,l)",goint_tensor,ooov_tensor,t1_tensor,MINUS_ONE) ! o(3)v(1)
            ierr=exatns_tensor_contract("G(k,j)+=F(k,c)*T(c,j)",goint_tensor,fov_tensor,t1_tensor) ! o(2)v(1)

!           ov ! validated
            !ierr=exatns_tensor_init(hovint_tensor)
            !ierr=exatns_tensor_init(hovaux_tensor)
            ierr=exatns_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",hovint_tensor,oovv_tensor,t1_tensor) ! o(2)v(2)
            ierr=exatns_tensor_add("L(k,c)+=H(k,c)",hovaux_tensor,hovint_tensor)                                  ! copy common term
            ierr=exatns_tensor_add("H(k,c)+=F(k,c)",hovint_tensor,fov_tensor)                                     ! add fock matrix
            ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hovaux_tensor,fov_tensor,one_tensor,MINUS_ONE) ! subtract fock matrix

!           oooo ! validated
            !ierr=exatns_tensor_init(aint_tensor)
            ierr=exatns_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",aint_tensor,oooo_tensor)
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",aint_tensor,ooov_tensor,t1_tensor) ! o(4)v(1)
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",aint_tensor,ooov_tensor,t1_tensor,MINUS_ONE) ! o(4)v(1)
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",aint_tensor,oovv_tensor,tau_tensor,ONE_HALF) ! o(4)v(2)

            if (ierr.ne.EXA_SUCCESS) then
               print*," ccdriver: error in calculation of intermediates, code=",ierr
               return
            end if

            if (debug) then
!           For debugging: print quasi-norm (NOT using complex conjugation so that we can also see whether tensors are imaginary)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=T(a,i)*T(a,i)",result_tensor,t1_tensor,t1_tensor)
               print*, "T1amp = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=T(a,b,i,j)*T(a,b,i,j)",result_tensor,t2_tensor,t2_tensor,ONE_QUARTER)
               print*, "T2amp = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=T(a,b,i,j)*T(a,b,i,j)",result_tensor,tau_tensor,tau_tensor,ONE_QUARTER)
               print*, "Tau   = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=H(k,j)*H(k,j)",result_tensor,hoint_tensor,hoint_tensor)
               print*, "hoint = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=G(k,j)*G(k,j)",result_tensor,goint_tensor,goint_tensor)
               print*, "goint = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=H(a,c)*H(a,c)",result_tensor,hvint_tensor,hvint_tensor)
               print*, "hvint = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=G(a,c)*G(a,c)",result_tensor,gvint_tensor,gvint_tensor)
               print*, "gvint = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=H(k,c)*H(k,c)",result_tensor,hovint_tensor,hovint_tensor)
               print*, "hovint = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=A(k,l,i,j)*A(k,l,i,j)",result_tensor,aint_tensor,aint_tensor)
               print*, "a int = ",result_tens(1)/4.D0
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=H(a,i,c,k)*H(a,i,c,k)",result_tensor,hint_tensor,hint_tensor)
               print*, "h int = ",result_tens(1)
               !ierr=exatns_tensor_init(result_tensor)
               ierr=exatns_tensor_contract("S()+=C(k,b,i,j)*C(k,b,i,j)",result_tensor,cint_tensor,cint_tensor)
               print*, "c int = ",result_tens(1)
            end if

!           This term appears twice as we need to permute i and j
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,goint_tensor,t2_tensor)                   ! validated o(3)v(2)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,goint_tensor,t2_tensor)                   ! validated o(3)v(2)
!           Again two times, also note the complex conjugation as we use vovv instead of vvvo
            ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,vovv_tensor,t1_tensor)                   ! validated o(2)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,vovv_tensor,t1_tensor,MINUS_ONE)         ! validated o(2)v(3)

!           This term appears twice as we need to permute a and b
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,gvint_tensor,t2_tensor)  ! validated o(2)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,gvint_tensor,t2_tensor)  ! validated o(2)v(3)
!           Again two times, with the C intermediate (note: in RELCCSD this is done differently)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,t1_tensor,MINUS_ONE) ! o(3)v(2)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,t1_tensor)           ! o(3)v(2)

!           The oooo contribution
            ierr=exatns_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,aint_tensor,tau_tensor,ONE_HALF) !validated o(4)v(2)

!           ierr=exatns_tensor_add("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!           work-around because complex conjugation does not work in the above statement
            ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,oovv_tensor,one_tensor)    ! validated o(2)v(2)
!           end of work-around

!           calculate s1 tensor
            !ierr=exatns_tensor_init(s1_tensor)
!           ierr=exatns_tensor_add("S(a,i)+=F+(i,a)",s1_tensor,fov_tensor)
!           work-around because complex conjugation does not work in the above statement
            ierr=exatns_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,fov_tensor,one_tensor)
!           end of work-around
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,vovv_tensor,tau_tensor,ONE_HALF) !validated o(2)v(3)
            ierr=exatns_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,hvint_tensor,t1_tensor)                  !validated o(1)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoint_tensor,t1_tensor)                        !validated o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,hovint_tensor,t2_tensor)                   !validated o(2)v(2)
            !ierr=exatns_tensor_init(hoaux_tensor)                                                                        !moretests
            ierr=exatns_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoaux_tensor,hovaux_tensor,t1_tensor)                    !moretests o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoaux_tensor,t1_tensor)                        !moretests o(2)v(1)
            ierr=exatns_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,ooov_tensor,tau_tensor,ONE_HALF) !validated o(3)v(2)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,vovo_tensor,t1_tensor,MINUS_ONE)     !validated o(2)v(2)

            if (ierr.ne.EXA_SUCCESS) then
               print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
               return
            end if

!           prepare for next iteration
            !ierr=exatns_tensor_init(t1_tensor)
            ierr=exatns_tensor_add("T1(a,i)+=S1(a,i)",t1_tensor,s1_tensor)
            !ierr=exatns_tensor_init(t2_tensor)
            ierr=exatns_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",t2_tensor,s2_tensor)

!           call diis (t1_tensor,t2_tensor,t_convergence=t_convergence)

            !ierr=exatns_tensor_init(result_tensor)
            ierr=exatns_tensor_contract("E()+=T(a,i)*F(i,a)",result_tensor,t1_tensor,fov_tensor) ! o(1)v(1)
            ccsd_energy=result_tens(1)
            !ierr=exatns_tensor_init(result_tensor)
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,tau_tensor,oovv_tensor) ! o(2)v(2)
            ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER
            write (*,'(A,I4,A,F20.15,A,E10.1)') "  Iteration ",iccsd," CCSD energy = ",ccsd_energy!, " convergence !=",t_convergence
!           if (t_convergence .lt. t_target_precision) exit
            exit   ! for testing first one cycle   
            end do

!           Calculate (T) energy correction, not implemented yet
!           call t3_energy (t_energy,eps_occ,eps_vir,t1_tensor,t2_tensor,fov_tensor,ooov_tensor,vovv_tensor,oovv_tensor)
!           write (*,'(A,F20.15)') "  4th order triples correction :    ",t_energy(1)
!           write (*,'(A,F20.15)') "  5th order triples (T) correction :",t_energy(2)
!           write (*,'(A,F20.15)') "  5th order triples -T  correction :",t_energy(3)
!           write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",ccsd_energy+t_energy(1)
!           write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(2)
!           write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(3)

!           Destroy tensors and shut down library
!           call diis (t1_tensor,t2_tensor,mode=2)
            ierr=exatns_tensor_destroy(oooo_tensor)
            ierr=exatns_tensor_destroy(oovv_tensor)
            ierr=exatns_tensor_destroy(vovo_tensor)
            ierr=exatns_tensor_destroy(vvvv_tensor)
            ierr=exatns_tensor_destroy(ooov_tensor)
            ierr=exatns_tensor_destroy(vovv_tensor)
            ierr=exatns_tensor_destroy(s1_tensor)
            ierr=exatns_tensor_destroy(s2_tensor)
            ierr=exatns_tensor_destroy(t1_tensor)
            ierr=exatns_tensor_destroy(t2_tensor)
            ierr=exatns_tensor_destroy(hoint_tensor)
            ierr=exatns_tensor_destroy(hoaux_tensor)
            ierr=exatns_tensor_destroy(hvint_tensor)
            ierr=exatns_tensor_destroy(hovint_tensor)
            ierr=exatns_tensor_destroy(hovaux_tensor)
            ierr=exatns_tensor_destroy(goint_tensor)
            ierr=exatns_tensor_destroy(gvint_tensor)
            ierr=exatns_tensor_destroy(aint_tensor)
            ierr=exatns_tensor_destroy(cint_tensor)
            ierr=exatns_tensor_destroy(hint_tensor)
!           Get statistics 
!           ierr=talsh_stats()    !Need to find correct exatns alternatives. Probably exatns_status

            !Stop ExaTENSOR runtime
            ierr=exatns_stop()

            end if

         else
            write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*," Final CCSD energy    = ",ccsd_energy
!        print*," Final CCSD(T) energy = ",ccsd_energy+t_energy(1)+t_energy(2)  ! no (T) implementation yet
         print*, " Leaving cc_driver routine"
         return

        end subroutine cc_driver

        subroutine get_CC_Fockmatrix (occ_space,vir_space,occ_space_id,vir_space_id,foo_tensor,fov_tensor,fvv_tensor)

!        Routine to get tensors with Fock matrix elements.
!        In this implementation we assume canonical orbitals, so we may simply construct zero tensors.
!        We also assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.
     
!        Written by Anastasios Papadopoulos, October 2018        
            
!        target tensors to be filled and given back to the caller
         type(tens_rcrsv_t), intent(inout)  :: foo_tensor, fov_tensor, fvv_tensor

         class(h_space_t), pointer  :: occ_space, vir_space
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
        
         integer                    :: ierr

         integer,dimension(2)       :: foo_id_array, fov_id_array, fvv_id_array
         integer(INTL), dimension(2)    :: foo_root_array, fov_root_array, fvv_root_array

         occ_space_root=occ_space%get_root_id(ierr)
         vir_space_root=vir_space%get_root_id(ierr)

         !For some reason these have to be defined. if directly inserted, won't work
         foo_id_array = (/occ_space_id,occ_space_id/)
         fov_id_array = (/occ_space_id,vir_space_id/)
         fvv_id_array = (/vir_space_id,vir_space_id/)
         foo_root_array = (/occ_space_root,occ_space_root/)
         fov_root_array = (/occ_space_root,vir_space_root/)
         fvv_root_array = (/vir_space_root,vir_space_root/)
 
         ierr=exatns_tensor_create(foo_tensor,"foo_tensor",foo_id_array,foo_root_array,EXA_DATA_KIND_C8) 
         ierr=exatns_tensor_create(fov_tensor,"fov_tensor",fov_id_array,fov_root_array,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(fvv_tensor,"fvv_tensor",fvv_id_array,fvv_root_array,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

         ierr=exatns_tensor_init(foo_tensor,ZERO)
         ierr=exatns_tensor_init(fov_tensor,ZERO)
         ierr=exatns_tensor_init(fvv_tensor,ZERO)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init(fock_matrix) failed!')

        end subroutine get_CC_Fockmatrix

        subroutine construct_4d_tensor (tensor,occ_space_id,vir_space_id,occ_space_root,vir_space_root,d1,d2,d3,d4)

!        Written by Anastasios Papadopoulos, October 2018

         type(tens_rcrsv_t), intent(inout)   :: tensor
         integer                             :: occ_space_id, vir_space_id
         integer(INTL)                       :: occ_space_root, vir_space_root
         integer, dimension(4)               :: tensor_id_array
         integer(INTL), dimension(4)         :: tensor_root_array
         character(len=3), intent(in)        :: d1, d2, d3, d4  ! "occ" or "vir"

         integer                             :: ierr

         !Create correct arrays => Need to figure out how to get the right names, cannot find dynamic variable names
         !tensor_id_array = (/d1_space_id, d2_space_id, d3_space_id, d4_space_id/)
         !tensor_root_array = (/d1_space_root, d2_space_root, d3_space_root, d4_space_root/)

         !Create 4d tensor
         ierr=exatns_tensor_create(tensor,"tensor",tensor_id_array,tensor_root_array,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
       
         !Initialize tensor to 0 
         ierr=exatns_tensor_init(tensor,ZERO)        
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init(tensor) failed!')

        end subroutine

        subroutine get_CC_integrals (nocc,nvir,mo_occ,mo_vir, &
                                    oooo_tensor,ooov_tensor,oovv_tensor, &
                                    vovo_tensor,vovv_tensor,vvvv_tensor)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017
!        ExaTensor conversion by Anastasios Papadopoulos, October 2018   

!        DISABLED ALL talsh_tensor_constructs AND DIM CALCS TO SUPPRESS ALL COMPILATIO ERRORS
         use exacorr_datatypes
         use exacorr_ao_to_mo
         use exacorr_global

         integer, intent(in)               :: nocc,nvir             ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)               :: mo_occ(:),mo_vir(:)   ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(tens_rcrsv_t), intent(inout) :: oooo_tensor, ooov_tensor, oovv_tensor, vovo_tensor, vovv_tensor, vvvv_tensor
         !integer                           :: oooo_dims(1:4),ooov_dims(1:4),oovv_dims(1:4),vovo_dims(1:4),vovv_dims(1:4),vvvv_dims(1:4)

         integer                           :: ierr, i
         complex(8), pointer, contiguous   :: ao_tens(:,:,:,:)

         integer    :: nao

!        auxilliary integral tensors (only needed inside this routine)
         type(tens_rcrsv_t) :: aoint_tensor, ovvo_tensor
         integer            :: aoint_dims(1:4), ovvo_dims(1:4)

!        arrays needed to communicate with the MO integral generator
         integer                :: nmo(4)
         integer, allocatable   :: mo_list(:)

!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t)    :: minusone_tensor
         integer               :: minusone_dims(1)

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get tensor with AO integrals
         aoint_dims      = nao
         !ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         call compute_ao_integrals (nao,ao_tens)

!        Get antisymmetric oovv tensor
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         !oovv_dims(1:2) = nocc * 2
         !oovv_dims(3:4) = nvir * 2
         !ierr=talsh_tensor_construct(oovv_tensor,C8,oovv_dims,init_val=ZERO)
!        whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
         allocate (mo_list(2*nocc+2*nvir))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nocc
         nmo(4) = nvir
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nocc+nvir+1:2*nocc+nvir)     = mo_occ(1:nocc)
         mo_list(2*nocc+nvir+1:2*nocc+2*nvir) = mo_vir(1:nvir)
         !call get_integral_tensor (aoint_tensor,oovv_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get antisymmetric oooo tensor
         !oooo_dims = nocc * 2
         !ierr=talsh_tensor_construct(oooo_tensor,C8,oooo_dims,init_val=ZERO)
         allocate (mo_list(4*nocc))
         nmo = nocc 
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         !call get_integral_tensor (aoint_tensor,oooo_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         !vovo_dims(1) = nvir * 2
         !vovo_dims(2) = nocc * 2
         !vovo_dims(3) = nvir * 2
         !vovo_dims(4) = nocc * 2
         !ierr=talsh_tensor_construct(vovo_tensor,C8,vovo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         !call get_integral_tensor (aoint_tensor,vovo_tensor,0,nmo,mo_list)
         deallocate (mo_list)
!        Get ovvo tensor
         !ovvo_dims(1) = nocc * 2
         !ovvo_dims(2) = nvir * 2
         !ovvo_dims(3) = nvir * 2
         !ovvo_dims(4) = nocc * 2
         !ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nvir
         nmo(4) = nocc
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nvir+nocc+1:2*nvir+nocc)     = mo_vir(1:nvir)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         !call get_integral_tensor (aoint_tensor,ovvo_tensor,0,nmo,mo_list)
         deallocate (mo_list)
!        Make anti-symmetrized vovo tensor
!        minusone_dims(1) = 1
         !ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
         ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",vovo_tensor,ovvo_tensor,minusone_tensor)
         ierr=exatns_tensor_destroy(ovvo_tensor)

!        Get vvvv tensor
         !vvvv_dims = nvir * 2
         !ierr=talsh_tensor_construct(vvvv_tensor,C8,vvvv_dims,init_val=ZERO)
         allocate (mo_list(4*nvir))
         nmo = nvir
         mo_list(1:nvir)          = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
         mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
         mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
         !call get_integral_tensor (aoint_tensor,vvvv_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get ooov tensor
         !ooov_dims(1:3) = nocc * 2
         !ooov_dims(4)   = nvir * 2
         !ierr=talsh_tensor_construct(ooov_tensor,C8,ooov_dims,init_val=ZERO)
         allocate (mo_list(3*nocc+nvir))
         nmo(1:3) = nocc
         nmo(4)   = nvir
         mo_list(1:nocc)               = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)        = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc)      = mo_occ(1:nocc)
         mo_list(3*nocc+1:3*nocc+nvir) = mo_vir(1:nvir)
         !call get_integral_tensor (aoint_tensor,ooov_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get vovv tensor
         !vovv_dims(1)   = nvir * 2
         !vovv_dims(2)   = nocc * 2
         !vovv_dims(3:4) = nvir * 2
         !ierr=talsh_tensor_construct(vovv_tensor,C8,vovv_dims,init_val=ZERO)
         allocate (mo_list(3*nvir+nocc))
         nmo(1)   = nvir
         nmo(2)   = nvir
         nmo(3)   = nocc
         nmo(4)   = nvir
         mo_list(1:nvir)                    = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)             = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)      = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:3*nvir+nocc) = mo_vir(1:nvir)
         !call get_integral_tensor (aoint_tensor,vovv_tensor,34,nmo,mo_list)
         deallocate (mo_list)

!        Clean-up
         ierr=exatns_tensor_destroy(aoint_tensor)

        end subroutine get_CC_integrals

end module exacorr_cc
