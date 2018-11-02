module exacorr_cc
!This module contains the routines needed to calculate ae coupled cluster wave function

        use exatensor

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

#if defined (VAR_MPI)
!-------------------------------------------------------------------------------------------
        !subroutine exacorr_cc_driver (nocc,nvir,mo_occ,mo_vir)
        subroutine exacorr_cc_driver (exa_input)

!        This subroutine drives the ccsd energy calculation
!        Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
 
!        Written by Lucas Visscher, winter 2016/2017

         use interface_to_mpi
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use exacorr_ao_to_mo
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

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(tens_rcrsv_t) :: foo_tensor, fvv_tensor, fov_tensor
         type(tens_rcrsv_t) :: oooo_tensor, ooov_tensor, oovv_tensor, vovo_tensor, vovv_tensor, vvvv_tensor
         type(tens_rcrsv_t) :: voov_tensor
!        intermediate tensors (same notation as in relccsd)
         type(tens_rcrsv_t) :: goint_tensor, gvint_tensor, hoint_tensor, hvint_tensor, hovint_tensor
         type(tens_rcrsv_t) :: aint_tensor, bint_tensor, cint_tensor, hint_tensor
         type(tens_rcrsv_t) :: tau_tensor, tau2_tensor
!        auxilliary tensors for t1 equation
         type(tens_rcrsv_t) :: hovaux_tensor, hoaux_tensor
!        solution tensors
         type(tens_rcrsv_t) :: t1_tensor, t2_tensor, s1_tensor, s2_tensor
!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t)         :: result_tensor!, one_tensor
!        CCSD control and result variables
         integer            :: iccsd, nccsd
         real(8)            :: mp2_energy, ccsd_energy, t_convergence!, t_energy(3) ! triples not implemented
         real(8), parameter :: t_target_precision = 1.D-12

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
         logical :: CCD=.false.

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
         ovoo_id = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvoo_id = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         voov_id = (/vir_space_id, occ_space_id, occ_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
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

            call print_date('starting ccsd calculation')

!           Create and initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(result_tensor,"result_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(one_tensor,"another_one_tensor",EXA_DATA_KIND_C8)   
            ! Note:  also get_integral_tensors will define a one_tensor and call it "one_tensor", so we needed another name
            ierr=exatns_tensor_init(one_tensor,ONE)

!           Get fixed one-body tensors => set to 0 for the moment (canonical orbitals)
            ierr=exatns_tensor_create(foo_tensor,"foo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(fov_tensor,"fov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(fvv_tensor,"fvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            call print_date('initialized fock matrix')

!           Get fixed two-body tensors
            ierr=exatns_tensor_create(oooo_tensor,"oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(ooov_tensor,"ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(oovv_tensor,"oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(vovo_tensor,"vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(vovv_tensor,"vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(vvvv_tensor,"vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(voov_tensor,"voov",voov_id,voov_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals 
            call get_CC_integrals (ao_space_id,ao_space_root, &
                                   oooo_tensor,ooov_tensor,oovv_tensor, &
                                   vovo_tensor,vovv_tensor,vvvv_tensor,voov_tensor)
            ierr=exatns_tensor_destroy(voov_tensor) ! no longer needed 

            call print_date('initialized two electron integrals')

!           Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
            if (.not.CCD) ierr=exatns_tensor_create(t1_tensor,"t1",vo_id,vo_root,EXA_DATA_KIND_C8)

!           Initialize t2 amplitudes with the MP2 values
            ierr=exatns_tensor_create(t2_tensor,"t2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
!           the tensor_add routine is not implemented in exatensor, work around this by using a contraction
            ierr=exatns_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)*K()",t2_tensor,oovv_tensor,one_tensor) ! o(2)v(2)
            ierr=exatns_tensor_transform(t2_tensor,'Denominate')

!           Calculate and print the MP2 energy
            ierr=exatns_tensor_init(result_tensor)
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor,ONE_QUARTER) ! o(2)v(2)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=real(result_val,8) ! Casting into real(8)
            print*," MP2 energy = ", mp2_energy

!           Start CCSD iterations

            nccsd = 30 ! hard wired for the moment

            ierr=exatns_tensor_create(goint_tensor,"goint",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(gvint_tensor,"gvint",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hoint_tensor,"hoint",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hoaux_tensor,"hoaux",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hvint_tensor,"hvint",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hovint_tensor,"hovint",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hovaux_tensor,"hovaux",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(aint_tensor,"aint",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(cint_tensor,"cint",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(hint_tensor,"hint",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(tau_tensor,"tau",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(tau2_tensor,"tau2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(s1_tensor,"s1",vo_id,vo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(s2_tensor,"s2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            if (.not.CCD) call diis (t1_tensor,t2_tensor,vo_id,vo_root,vvoo_id,vvoo_root,mode=1)
            if (ierr.eq.EXA_SUCCESS) then
               print*," ccdriver: allocated all tensors"
            else
               print*," ccdriver: error in tensor allocation, code=",ierr
               return
            end if

            do iccsd = 1, nccsd

!           make tau ! validated
            ierr=exatns_tensor_init(tau_tensor)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,b,i,j)*K()",tau_tensor,t2_tensor,one_tensor)
            if (.not.CCD) then
                ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",tau_tensor,t1_tensor,t1_tensor) ! o(2)v(2)
                ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",tau_tensor,t1_tensor,t1_tensor,MINUS_ONE) ! o(2)v(2)
            end if

!           NOTE: contractions marked with "validated" gave numerically significant contributions and could be validated by comparison
!           to the reference implementation, the ones marked "moretests" were too small to decide whether it was correct or not.
!           calculate intermediaries

            ierr=exatns_tensor_init(s2_tensor)

!           The vvvv contribution
            ierr=exatns_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,vvvv_tensor,tau_tensor,ONE_HALF) !validated o(2)v(4)

!           vv ! validated
            ierr=exatns_tensor_init(hvint_tensor)
            ierr=exatns_tensor_contract("H(a,c)+=F(a,c)*K()",hvint_tensor,fvv_tensor,one_tensor)
            ierr=exatns_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",hvint_tensor,oovv_tensor,tau_tensor,ONE_HALF) ! o(2)v(3)
            ierr=exatns_tensor_init(gvint_tensor)
            ierr=exatns_tensor_contract("G(a,c)+=H(a,c)*K()",gvint_tensor,hvint_tensor,one_tensor)
            if (.not.CCD) then
                ierr=exatns_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",gvint_tensor,vovv_tensor,t1_tensor) ! o(1)v(3)
                ierr=exatns_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",gvint_tensor,fov_tensor,t1_tensor,MINUS_ONE) ! o(1)v(2)
            end if

!           vovo ! moretests
            ierr=exatns_tensor_init(tau2_tensor)
            ierr=exatns_tensor_contract("T(a,d,l,i)+=R(a,d,l,i)*K()",tau2_tensor,t2_tensor,one_tensor,ONE_HALF)
            if (.not.CCD) ierr=exatns_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,t1_tensor,t1_tensor) ! o(2)v(2)
            ierr=exatns_tensor_init(hint_tensor)
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",hint_tensor,oovv_tensor,tau2_tensor)! o(3)v(3)
            if (.not.CCD) then
                ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",hint_tensor,vovv_tensor,t1_tensor) ! o(2)v(3)
                ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",hint_tensor,ooov_tensor,t1_tensor,MINUS_ONE) ! o(3)v(2)
            end if
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,i)*K()",hint_tensor,vovo_tensor,one_tensor)

!           This term appears four times as we need to permute both i and j and a and b
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)

!           ovoo ! validated
            ierr=exatns_tensor_init(cint_tensor)
            ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,vovv_tensor,tau_tensor,MINUS_ONE_HALF) ! o(3)v(3)
            if (.not.CCD) then
                ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,vovo_tensor,t1_tensor) ! o(3)v(2)
                ierr=exatns_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,vovo_tensor,t1_tensor,MINUS_ONE) ! o(3)v(2)
            end if
            ierr=exatns_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,ooov_tensor,one_tensor) ! o(3)v(1)

!           oo ! validated
            ierr=exatns_tensor_init(hoint_tensor)
            ierr=exatns_tensor_contract("H(k,j)+=F(k,j)*K()",hoint_tensor,foo_tensor,one_tensor)
            ierr=exatns_tensor_contract("H(k,j)+=V(k,l,c,d)*T(c,d,l,j)",hoint_tensor,oovv_tensor,tau_tensor,ONE_HALF) ! o(3)v(2)
            ierr=exatns_tensor_init(goint_tensor)
            ierr=exatns_tensor_contract("G(k,j)+=H(k,j)*K()",goint_tensor,hoint_tensor,one_tensor)
            if (.not.CCD) then 
                ierr=exatns_tensor_contract("G(k,j)+=V(k,l,j,c)*T(c,l)",goint_tensor,ooov_tensor,t1_tensor,MINUS_ONE) ! o(3)v(1)
                ierr=exatns_tensor_contract("G(k,j)+=F(k,c)*T(c,j)",goint_tensor,fov_tensor,t1_tensor) ! o(2)v(1)
            end if

!           ov ! validated
            ierr=exatns_tensor_init(hovint_tensor)
            ierr=exatns_tensor_init(hovaux_tensor)
            if (.not.CCD) ierr=exatns_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",hovint_tensor,oovv_tensor,t1_tensor) ! o(2)v(2)
            ierr=exatns_tensor_contract("L(k,c)+=H(k,c)*K()",hovaux_tensor,hovint_tensor,one_tensor)        ! copy common term
            ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hovint_tensor,fov_tensor,one_tensor)           ! add fock matrix
            ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",hovaux_tensor,fov_tensor,one_tensor,MINUS_ONE) ! subtract fock matrix

!           oooo ! validated
            ierr=exatns_tensor_init(aint_tensor)
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,j)*K()",aint_tensor,oooo_tensor,one_tensor)
            if (.not.CCD) then 
                ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",aint_tensor,ooov_tensor,t1_tensor) ! o(4)v(1)
                ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",aint_tensor,ooov_tensor,t1_tensor,MINUS_ONE) ! o(4)v(1)
            end if
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",aint_tensor,oovv_tensor,tau_tensor,ONE_HALF) ! o(4)v(2)

            if (.true.) then
               if (.not.CCD) then   
                write(*,'(A7,F20.15)'), "T1amp = ",print_tensornorm2 (t1_tensor)
               end if
               write(*,'(A7,F20.15)'), "T2amp = ",print_tensornorm2 (t2_tensor)
               write(*,'(A7,F20.15)'), "Tau   = ",print_tensornorm2 (tau_tensor)
               write(*,'(A7,F20.15)'), "hoint = ",print_tensornorm2 (hoint_tensor)
               write(*,'(A7,F20.15)'), "goint = ",print_tensornorm2 (goint_tensor)
               write(*,'(A7,F20.15)'), "hvint = ",print_tensornorm2 (hvint_tensor)
               write(*,'(A7,F20.15)'), "gvint = ",print_tensornorm2 (gvint_tensor)
               write(*,'(A7,F20.15)'), "hovint= ",print_tensornorm2 (hovint_tensor)
               write(*,'(A7,F20.15)'), "aint  = ",print_tensornorm2 (aint_tensor)
               write(*,'(A7,F20.15)'), "hint  = ",print_tensornorm2 (hint_tensor)
               write(*,'(A7,F20.15)'), "cint  = ",print_tensornorm2 (cint_tensor)
            end if

!           This term appears twice as we need to permute i and j
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,goint_tensor,t2_tensor)                   ! validated o(3)v(2)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,goint_tensor,t2_tensor)                   ! validated o(3)v(2)
!           Again two times, also note the complex conjugation as we use vovv instead of vvvo
            if (.not.CCD) then    
                ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,vovv_tensor,t1_tensor)                   ! validated o(2)v(3)
                ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,vovv_tensor,t1_tensor,MINUS_ONE)         ! validated o(2)v(3)
            end if

!           This term appears twice as we need to permute a and b
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,gvint_tensor,t2_tensor)  ! validated o(2)v(3)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,gvint_tensor,t2_tensor)  ! validated o(2)v(3)
!           Again two times, with the C intermediate (note: in RELCCSD this is done differently)
            if (.not.CCD) then    
                ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,t1_tensor,MINUS_ONE) ! o(3)v(2)
                ierr=exatns_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,t1_tensor)           ! o(3)v(2)
            end if

!           The oooo contribution
            ierr=exatns_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,aint_tensor,tau_tensor,ONE_HALF) !validated o(4)v(2)
            ierr=exatns_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,oovv_tensor,one_tensor)    ! validated o(2)v(2)

!           calculate s1 tensor
            ierr=exatns_tensor_init(s1_tensor)
            ierr=exatns_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,fov_tensor,one_tensor)
            ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,vovv_tensor,tau_tensor,ONE_HALF) !validated o(2)v(3)
            if (.not.CCD) then    
                ierr=exatns_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,hvint_tensor,t1_tensor)                  !validated o(1)v(2)
                ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoint_tensor,t1_tensor)                        !validated o(2)v(1)
            end if
            ierr=exatns_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,hovint_tensor,t2_tensor)                   !validated o(2)v(2)
            ierr=exatns_tensor_init(hoaux_tensor)                                                                        !moretests
            if (.not.CCD) then    
                ierr=exatns_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoaux_tensor,hovaux_tensor,t1_tensor)                    !moretests o(2)v(1)
                ierr=exatns_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoaux_tensor,t1_tensor)                        !moretests o(2)v(1)
            end if
            ierr=exatns_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,ooov_tensor,tau_tensor,ONE_HALF) !validated o(3)v(2)
            if (.not.CCD) ierr=exatns_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,vovo_tensor,t1_tensor,MINUS_ONE)     !validated o(2)v(2)

            if (ierr.ne.EXA_SUCCESS) then
               print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
               return
            end if

!           prepare for next iteration
            if (.not.CCD) then
                ierr=exatns_tensor_init(t1_tensor,ZERO)
                ierr=exatns_tensor_contract("T1(a,i)+=S1(a,i)*K()",t1_tensor,s1_tensor,one_tensor)
                ierr=exatns_tensor_transform(t1_tensor,'Denominate')
            end if
            ierr=exatns_tensor_init(t2_tensor,ZERO)
            ierr=exatns_tensor_contract("T2(a,b,i,j)+=S2(a,b,i,j)*K()",t2_tensor,s2_tensor,one_tensor)
            ierr=exatns_tensor_transform(t2_tensor,'Denominate')

            if (.not.CCD) call diis (t1_tensor,t2_tensor,vo_id,vo_root,vvoo_id,vvoo_root,t_convergence=t_convergence)

            ierr=exatns_tensor_init(result_tensor)
            if (.not.CCD) then
!               Will always be 0 contribution for now, since fov=0 (canonical orbitals)                
                ierr=exatns_tensor_contract("E()+=T(a,i)*F(i,a)",result_tensor,t1_tensor,fov_tensor) ! o(1)v(1)
                ierr=exatns_tensor_get_scalar(result_tensor,result_val)
                ccsd_energy=real(result_val,8)
                ierr=exatns_tensor_init(result_tensor)
            end if
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,tau_tensor,oovv_tensor) ! o(2)v(2)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            ccsd_energy = ccsd_energy + real(result_val*ONE_QUARTER,8)
            write (*,'(A,I4,A,F20.15,A,E10.1)') "  Iteration ",iccsd," CCSD energy = ",ccsd_energy, " convergence =",t_convergence
            if (t_convergence .lt. t_target_precision) exit
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
            if (.not.CCD) call diis (t1_tensor,t2_tensor,mode=2)
            ierr=exatns_tensor_destroy(one_tensor)
            ierr=exatns_tensor_destroy(result_tensor)
            ierr=exatns_tensor_destroy(foo_tensor)
            ierr=exatns_tensor_destroy(fov_tensor)
            ierr=exatns_tensor_destroy(fvv_tensor)
            ierr=exatns_tensor_destroy(oooo_tensor)
            ierr=exatns_tensor_destroy(oovv_tensor)
            ierr=exatns_tensor_destroy(vovo_tensor)
            ierr=exatns_tensor_destroy(vvvv_tensor)
            ierr=exatns_tensor_destroy(ooov_tensor)
            ierr=exatns_tensor_destroy(vovv_tensor)
            ierr=exatns_tensor_destroy(s1_tensor)
            ierr=exatns_tensor_destroy(s2_tensor)
            if (.not.CCD) ierr=exatns_tensor_destroy(t1_tensor)
            ierr=exatns_tensor_destroy(t2_tensor)
            ierr=exatns_tensor_destroy(tau_tensor)
            ierr=exatns_tensor_destroy(tau2_tensor)
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
            !ierr=exatns_status()    !Need to find correct exatns alternatives. Probably exatns_status

            print*," Final CCSD energy = ",ccsd_energy
!           print*," Final CCSD(T) energy = ",ccsd_energy+t_energy(1)+t_energy(2)  ! no (T) implementation yet
            
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

        subroutine get_CC_integrals (ao_space_id,ao_space_root, &
                                    oooo_tensor,ooov_tensor,oovv_tensor, &
                                    vovo_tensor,vovv_tensor,vvvv_tensor,voov_tensor)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017

         use exacorr_datatypes
         use exacorr_ao_to_mo
         use exacorr_global

!        information needed to compute the ao-tensor
         integer,intent(in)                :: ao_space_id
         integer(INTL),intent(in)          :: ao_space_root
!        target integral tensors to be filled and given back to the caller
         type(tens_rcrsv_t), intent(inout) :: oooo_tensor, ooov_tensor, oovv_tensor, vovo_tensor, vovv_tensor, vvvv_tensor
         type(tens_rcrsv_t), intent(inout) :: voov_tensor ! Auxilliary tensor that is not anti-symmetrized

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
         call get_integral_tensor (aoint_tensor,oovv_tensor,12,initlabel)

!        Get oooo tensor
         initlabel(:,1)=initocc
         initlabel(:,2)=initocc
         initlabel(:,3)=initocc
         initlabel(:,4)=initocc
         call get_integral_tensor (aoint_tensor,oooo_tensor,12,initlabel)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initvir
         initlabel(:,3)=initocc
         initlabel(:,4)=initocc
         call get_integral_tensor (aoint_tensor,vovo_tensor,0,initlabel)
!        Get voov tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initocc
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,voov_tensor,0,initlabel)
!        Make anti-symmetrized vovo tensor
         ierr=exatns_tensor_contract("V(a,i,b,j)+=V(a,i,j,b)*K()",vovo_tensor,voov_tensor,one_tensor,MINUS_ONE)

!        Get vvvv tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initvir
         initlabel(:,3)=initvir
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,vvvv_tensor,12,initlabel)

!        Get ooov tensor
         initlabel(:,1)=initocc
         initlabel(:,2)=initocc
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,ooov_tensor,12,initlabel)

!        Get vovv tensor
         initlabel(:,1)=initvir
         initlabel(:,2)=initvir
         initlabel(:,3)=initocc
         initlabel(:,4)=initvir
         call get_integral_tensor (aoint_tensor,vovv_tensor,34,initlabel)

!        Clean-up
         ierr=exatns_tensor_destroy(aoint_tensor)

        end subroutine get_CC_integrals

        subroutine diis (t1_tensor,t2_tensor,vo_id,vo_root,vvoo_id,vvoo_root,t_convergence,mode)
        
         integer, parameter             :: maxdim = 8 ! maximum number of tensors that will be stored 
         integer, parameter             :: mindim = 4 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes
        
         type(tens_rcrsv_t), intent(inout)   :: t1_tensor, t2_tensor
         integer                             :: ierr, i, tens_rank, j
         integer,dimension(2),optional       :: vo_id
         integer,dimension(4),optional       :: vvoo_id
         integer(INTL),dimension(2),optional :: vo_root
         integer(INTL),dimension(4),optional :: vvoo_root
         logical                             :: init, clean_up
         integer, save                       :: itens, jtens, ntens
         complex(8), save                    :: b_diis(maxdim,maxdim)
         complex(8)                          :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         complex(8)                          :: scalar
         real(8)                             :: b_max, b_scale
         type(tens_rcrsv_t), save            :: t1_inp_tensor, t2_inp_tensor
         type(tens_rcrsv_t), save            :: t1_res_tensor(maxdim), t2_res_tensor(maxdim)
         type(tens_rcrsv_t), save            :: t1_err_tensor(maxdim), t2_err_tensor(maxdim)
         integer                             :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)                             :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(tens_rcrsv_t), save            :: scalar_tensor
       
         character(len=1) :: vs_str
 
         init       = .false.
         clean_up   = .false.
        
         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if 

         if (init) then
!           get dimensions of the tensors
            tens_rank=t1_tensor%get_rank(ierr)
            if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
            !ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
            tens_rank=t2_tensor%get_rank(ierr)
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            ierr=exatns_tensor_create(t1_inp_tensor,"t1_inp_tensor",vo_id,vo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(t2_inp_tensor,"t2_inp_tensor",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            do i = 1, maxdim
               write (vs_str, '(I1)') i
               ierr=exatns_tensor_create(t1_res_tensor(i),"t1_res_tensor_" // vs_str,vo_id,vo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_create(t1_err_tensor(i),"t1_err_tensor_" // vs_str,vo_id,vo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_create(t2_res_tensor(i),"t2_res_tensor_" // vs_str,vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
               ierr=exatns_tensor_create(t2_err_tensor(i),"t2_err_tensor_" // vs_str,vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            end do
            if(ierr.ne.EXA_SUCCESS) call quit('failed to create storage for input')
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            ierr=exatns_tensor_init(t1_inp_tensor)
            ierr=exatns_tensor_contract("T1O(a,i)+=T1N(a,i)*K()",t1_inp_tensor,t1_tensor,one_tensor)
            ierr=exatns_tensor_init(t2_inp_tensor)
            ierr=exatns_tensor_contract("T2O(a,b,i,j)+=T2N(a,b,i,j)*K()",t2_inp_tensor,t2_tensor,one_tensor)
            ierr=exatns_tensor_create(scalar_tensor,"scalar_tensor",EXA_DATA_KIND_C8)
            return ! no extrapolation possible
         end if
       
         if (clean_up) then
!           de-allocate tensors
            ierr=exatns_tensor_destroy(t1_inp_tensor)
            ierr=exatns_tensor_destroy(t2_inp_tensor)
            do i = 1, maxdim
                ierr=exatns_tensor_destroy(t1_res_tensor(i))
                ierr=exatns_tensor_destroy(t1_err_tensor(i))
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
         ierr=exatns_tensor_init(t1_res_tensor(itens))
         ierr=exatns_tensor_contract("T1R(a,i)+=T1(a,i)*K()",t1_res_tensor(itens),t1_tensor,one_tensor)
         ierr=exatns_tensor_init(t2_res_tensor(itens))
         ierr=exatns_tensor_contract("T2R(a,b,i,j)+=T2(a,b,i,j)*K()",t2_res_tensor(itens),t2_tensor,one_tensor)
     
!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         ierr=exatns_tensor_init(t1_err_tensor(itens))
         ierr=exatns_tensor_contract("T1E(a,i)+=T1R(a,i)*K()",t1_err_tensor(itens),t1_tensor,one_tensor)
         ierr=exatns_tensor_contract("T1E(a,i)+=T1I(a,i)*K()",t1_err_tensor(itens),t1_inp_tensor,one_tensor,prefactor=MINUS_ONE)
         ierr=exatns_tensor_init(t2_err_tensor(itens))
         ierr=exatns_tensor_contract("T2E(a,b,i,j)+=T2R(a,b,i,j)*K()",t2_err_tensor(itens),t2_tensor,one_tensor)
         ierr=exatns_tensor_contract("T2E(a,b,i,j)+=T2I(a,b,i,j)*K()",t2_err_tensor(itens),t2_inp_tensor,one_tensor, &
                                                                                                prefactor=MINUS_ONE)
 
!        update the b_diis array that holds the scalar products of the error tensors
         do i = 1, ntens
            ierr=exatns_tensor_init(scalar_tensor)
            ierr=exatns_tensor_contract("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
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
            ierr=exatns_tensor_init(t1_tensor)
            ierr=exatns_tensor_init(t2_tensor)
            do i = 1, ntens
               ierr=exatns_tensor_contract("T1(a,i)+=T1O(a,i)*K()",t1_tensor,t1_res_tensor(i),one_tensor,prefactor=coeff(i))
               ierr=exatns_tensor_contract("T2(a,b,i,j)+=T2O(a,b,i,j)*K()",t2_tensor,t2_res_tensor(i),one_tensor,prefactor=coeff(i))
            end do
         end if
        
!        the extrapolated amplitudes are the new input tensors
         ierr=exatns_tensor_init(t1_inp_tensor)
         ierr=exatns_tensor_contract("T1O(a,i)+=T1(a,i)*K()",t1_inp_tensor,t1_tensor,one_tensor)
         ierr=exatns_tensor_init(t2_inp_tensor)
         ierr=exatns_tensor_contract("T2O(a,b,i,j)+=T2(a,b,i,j)*K()",t2_inp_tensor,t2_tensor,one_tensor)

        end subroutine diis

#else
        subroutine exacorr_cc_driver (nocc,nvir,mo_occ,mo_vir)
         implicit none
         integer, intent(in) :: nocc,nvir           ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in) :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         print*,  'called exatensor version of cc_driver in a serial version'
         print*,  'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_cc_driver

#endif

end module exacorr_cc
