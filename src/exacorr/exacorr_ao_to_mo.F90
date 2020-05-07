module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exatensor
        use exacorr_utils
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private
        public exacorr_compute_mo_integrals
#if defined (VAR_MPI)

        interface get_integral_tensor
            module procedure get_integral_tensor
        end interface get_integral_tensor

        public mocoef_initlabel
        public dm_initlabel
        public get_integral_tensor
        public get_ht2_integral_tensor
        public print_tensornorm2
        public ao2mo_exat
        public mulliken_to_dirac_sort

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         ! It is not possible to make an array of pointers, so we need a container
         type,public ::  h_spaces_container_t
              class(h_space_t), pointer :: ptr
         end type h_spaces_container_t

       contains

        subroutine exacorr_compute_mo_integrals (nmo, mo_list,print_level)

!        This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!        It currently tests two different algorithms as well as the get_integral_tensor routine used by the CC code.
!        Intended final purpose is two run the most efficient algorithm only and write integrals out to a file that can be used by correlation methods.
!        See talsh_compute_mo_integrals to see how this works for the single node implementation.

!        Written by Lucas Visscher, march/september 2018

         use interface_to_mpi
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use exatensor
         

         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         integer, intent(in) :: print_level

         integer(INTD)       :: my_role
         integer             :: my_MPI_rank    ! needed to start exatensor, MPI should already be started in threaded mode
         integer             :: nshells                 ! the number of shells (orbitals with the same center and l-value)
         integer             :: nao                     ! the number of atomic orbitals)
         integer             :: basis_angular           ! 1=cartesian, 2=spherical
         type(basis_func_info_t), allocatable :: gto(:) ! array with information about the shells

         ! classes for the definition of spaces (AO and 4 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         type(h_spaces_container_t), dimension(4) :: mo_space

         ! variables for the definition of spaces (AO and 4 different MO spaces)
         type(subspace_basis_t):: basis_ao, basis_mo(4)
         integer(INTL)::          nao_l, ao_space_root
         integer(INTL)::          nmo_l, mo_space_root(4)
         integer(INTD)::          ao_space_id, mo_space_id(4)

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL)::                        labs

         ! types of the initializer methods
         type(compute_2e_ao_tensor_t) :: aoint_calculator ! method to initialize the ao-integrals tensor
         type(init_mocoef_tensor_t)   :: mocoef_init(2,4) ! array of initializers for algorithm 1, first index is alpha/beta, second is 1,2,3,4
         type(init_dm_tensor_t)       :: dm_init(2)       ! array of initializers for algorithm 2, index indicates first/second halftransformation
         type(init_ht_tensor_t)       :: ht_init          ! initializers for algorithm 3 in which the HT matrix is made without ao-integral tensor
         type(tens_printer_t)         :: tens_printer     ! generic method, but needs to be defined here in specific form

         ! the two tensors that we need
         type(tens_rcrsv_t):: aoint_tensor
         type(tens_rcrsv_t):: moint_tensor

         ! loop and other auxilliary variables
         integer(INTD) ::  spin,ierr, lsh
         integer:: i,imo,jmo,imo12,jmo12,imo34,jmo34
         integer:: l, nl
         character(len=2)  :: vs_str
         integer(INTD) branch_factor

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) print*, " Starting AO-to-MO transform with exatensor"

!        Initialize global data (read from DFCOEF or similar by master)
         call set_talsh_only(.false.)
         call initialize_global_data

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

         !Register the ao vector space (block size hardwired to 75 for the moment):
         branch_factor = nao / 75 + 1
         if (branch_factor < 2) branch_factor = 2
         ierr=exatns_space_register('AOspace',basis_ao,ao_space_id,ao_space,branch_factor)
         if (my_MPI_rank == 0) print*," registered AO vector space with nao, nshells:",nao,nshells
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Get the MO-coefficient basis for the four indices that are to be transformed
         ierr = 0
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            !nmo_l   = 2*int(nmo(i),8)  ! NB: imo, jmo and nmo count orbitals, nmo_l counts spinors
            nmo_l   = int(nmo(i),8)  ! NB: imo, jmo and nmo count orbitals, nmo_l counts spinors
            call basis_mo(i)%subspace_basis_ctor(nmo_l,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

            !Set the functions inside this basis (without colors, not needed for MOs)
            do labs = 1, nmo_l
               call basis_mo(i)%set_basis_func(labs,BASIS_ABSTRACT,ierr)
               if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
            enddo

            !Finalize the mo space basis
            call basis_mo(i)%finalize(ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

            !Register this mo vector space :
            branch_factor = nmo(i) / 75 + 1
            if (branch_factor < 2) branch_factor = 2
            write (vs_str,'(A,I1)') '_',i
            ierr=exatns_space_register('MOspace'//vs_str,basis_mo(i),mo_space_id(i),mo_space(i)%ptr,branch_factor)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

            imo = imo + nmo(i) ! make imo point to the start of the list for the next index
         end do

         !Define separate initializers for the alpha(A) and beta(B) MO-coefficients for each index
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            do spin = 1, 2
               !Register this initializer as a method
               call mocoef_init(spin,i)%init_mocoef_ctor(mo_list(imo:jmo),nmo(i),spin)
               ierr=exatns_method_register(mocoef_initlabel(ao_space_id,mo_space_id(i),spin),mocoef_init(spin,i))
               if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            end do
         end do

         !Define separate initializers for the first and second density matrix
         imo12 = 1
         jmo12 = nmo(1)+nmo(2)
         imo34 = jmo12 + 1
         jmo34 = jmo12 + nmo(3)+nmo(4)
         !Register the dm initializers as methods
         call dm_init(1)%init_dm_ctor(mo_list(imo12:jmo12),nmo(1:2))
         ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),mo_space_id(1:2)),dm_init(1))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         call dm_init(2)%init_dm_ctor(mo_list(imo34:jmo34),nmo(3:4))
         ierr=exatns_method_register(dm_initlabel((/ao_space_id,ao_space_id/),mo_space_id(3:4)),dm_init(2))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         if (my_MPI_rank == 0) print*," registered MO vectors spaces and initializers with nmo:",nmo

         !Define separate initializer for the construction of the half-transformed integrals (algorithm 3)
         imo12 = 1
         jmo12 = nmo(1)+nmo(2)
         !Register the ht initializer as method
         call ht_init%init_ht_ctor(mo_list(imo12:jmo12),nmo(1:2))
         ierr=exatns_method_register('HalfTrans',ht_init)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the ao integral generation as a method
         ierr=exatns_method_register('AointCalculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the standard print method with a high treshold to print only important tensor elements
         call tens_printer%reset_thresh(0.01D0)
         ierr=exatns_method_register('TensorPrinter',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do. 
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(int(global_communicator,INTD))

         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
          if(my_role.eq.EXA_DRIVER) then

           call print_date('starting ao integral calculation')

           !Create ao-integral-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'AointCalculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('finished ao integral calculation')

           ! the  next three lines are for debugging,
           print*, "aoint = ", print_tensornorm2 (aoint_tensor) ! this is usually sufficent for checking
           ! when uncommenting the next line: keep in mind that printing an entire tensor yields an enormous output
           ! ierr=exatns_tensor_traverse(aoint_tensor,'TensorPrinter')

           !Create mo-integral-tensor:
           do i = 1, 4
              mo_space_root(i)=mo_space(i)%ptr%get_root_id(ierr)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           call print_date('created mo integral tensor')

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call aomo_transform_1 (aoint_tensor,moint_tensor,print_level)
           call print_date('completed mo integral tensor: algorithm 1')
           print*, "moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call aomo_transform_2 (aoint_tensor,moint_tensor,print_level)
           call print_date('completed mo integral tensor: algorithm 2')
           print*, "moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call aomo_transform_3 (ao_space_id,ao_space_root,moint_tensor,print_level)
           call print_date('completed mo integral tensor: algorithm 3')
           print*,"moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call get_integral_tensor (aoint_tensor,moint_tensor,0,print_level)
           call print_date('completed dirac-sorted mo integral tensor - not anti-symmetrized')
           print*,"moint = ", print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_destroy(moint_tensor)

           !Recreate mo-integral-tensor with equal index spaces so that we can anti-symmetrize:
           do i = 2, 4
              mo_space_id(i)=mo_space_id(1)
              mo_space_root(i)=mo_space_root(1)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           
           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call get_integral_tensor (aoint_tensor,moint_tensor,12,print_level)
           call print_date('completed dirac-sorted mo integral tensor - 12 anti-symmetrized')
           print*,"moint = ", print_tensornorm2 (moint_tensor)
           
           ierr=exatns_tensor_init(moint_tensor,ZERO)
           call get_integral_tensor(aoint_tensor,moint_tensor,34,print_level)
           call print_date('completed dirac-sorted mo integral tensor - 34 anti-symmetrized')
           print*, "moint = ", print_tensornorm2 (moint_tensor)

           !Destroy aoint-tensor:
           ierr=exatns_tensor_destroy(aoint_tensor)

           !Destroy moint-tensor:
           ierr=exatns_tensor_destroy(moint_tensor)

           !Stop ExaTENSOR runtime:
           ierr=exatns_stop()

          end if

         else
          write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Process ",my_MPI_rank, " leaving compute_mo_integrals_exa routine"
         return

        end subroutine exacorr_compute_mo_integrals

        subroutine aomo_transform_1 (aoint_tensor,moint_tensor,print_level)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 1: n^5 scaling with index-by-index transform

         implicit none

         ! AO integrals are to be provided on input and are not changed, MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor   ! for technical reasons intent(inout) needs to be used for all tensors
         integer, intent(in)                :: print_level
         ! Auxilliary tensors that are created and destroyed inside this subroutine
         type(tens_rcrsv_t):: a2int_tensor

         integer(INTD)     :: ierr
         integer(INTD)     :: aoint_id(4),moint_id(4),a2int_id(4)
         integer(INTL)     :: aoint_root(4),moint_root(4),a2int_root(4)

         if (print_level >5) call print_date('entered aomo_transform_1')

!        Determine space ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the auxilliary one and allocate it
         a2int_id(1:2) = moint_id(1:2)
         a2int_id(3:4) = aoint_id(3:4)
         a2int_root(1:2) = moint_root(1:2)
         a2int_root(3:4) = aoint_root(3:4)
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(a2int_tensor,ZERO)

!        Do 1st and 2nd quarter-transformation
         call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
         if (print_level >5) call print_date('first ht done')
!        Do 3rd and 4th quarter-transformation
         call ao2mo_exat (a2int_tensor,moint_tensor,print_level)
         if (print_level >5) call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=35; return; endif

         if (print_level >5) call print_date('leaving aomo_transform_1')

        end subroutine aomo_transform_1

        subroutine aomo_transform_2 (aoint_tensor,moint_tensor,print_level)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 2: n^6 scaling with density matrices

         implicit none

         ! AO integrals are to be provided on input and are not changed, MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor
         integer,intent(in)                 :: print_level

         ! Auxilliary tensor needed for this algorithm
         type(tens_rcrsv_t):: dm_tensor
         type(tens_rcrsv_t):: a2int_tensor

         integer(INTD)     :: ierr
         integer(INTD)     :: aoint_id(4),moint_id(4),a2int_id(4),dm_id(4)
         integer(INTL)     :: aoint_root(4),moint_root(4),a2int_root(4),dm_root(4)

         if (print_level >5) call print_date('entered aomo_transform_2')

!        Determine ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Create and initialize the density matrix tensor for the first half-transformation
         dm_id(1:2)    = aoint_id(1:2)
         dm_id(3:4)    = moint_id(1:2)
         dm_root(1:2)  = aoint_root(1:2)
         dm_root(3:4)  = moint_root(1:2)
         ierr=exatns_tensor_create(dm_tensor,'dm12tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm1) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,dm_initlabel(dm_id(1:2),dm_id(3:4)))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         if (print_level >5) call print_date('initialized first density matrix')

!        Copy the information about spaces from the starting and target tensors to the ht integral tensor
         a2int_id(1:2) = moint_id(1:2)
         a2int_id(3:4) = aoint_id(3:4)
         a2int_root(1:2) = moint_root(1:2)
         a2int_root(3:4) = aoint_root(3:4)

!        Do 1st half-transformation
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('ERROR in aomo_transform_2: int. tens. too large ')
         ierr=exatns_tensor_init(a2int_tensor,ZERO)
         if(ierr.ne.EXA_SUCCESS) call quit('ERROR in aomo_transform_2: initialization failed ')
         if (print_level >5) call print_date('constructed a2')
         ierr=exatns_tensor_contract(a2int_tensor,aoint_tensor,dm_tensor,"D(i,j,r,s)+=L(p,q,r,s)*R(p,q,i,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=1; return; endif

!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         if (print_level >5) call print_date('first ht done')

!        Create and initialize the density matrix tensor for the second half-transformation
         dm_id(1:2)    = aoint_id(3:4)
         dm_id(3:4)    = moint_id(3:4)
         dm_root(1:2)  = aoint_root(3:4)
         dm_root(3:4)   = moint_root(3:4)
         ierr=exatns_tensor_create(dm_tensor,'dm34tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,dm_initlabel(dm_id(1:2),dm_id(3:4)))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         if (print_level >5) call print_date('initialized second density matrix')

!        Do 2nd half-transformation
         ierr=exatns_tensor_contract(moint_tensor,a2int_tensor,dm_tensor,"D(i,j,k,l)+=L(i,j,r,s)*R(r,s,k,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=2; return; endif

         if (print_level >5) call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif
!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         if (print_level >5) call print_date('leaving aomo_transform_2')

        end subroutine aomo_transform_2

        subroutine aomo_transform_3 (ao_space_id,ao_space_root,moint_tensor,print_level)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 3: n^6 scaling with fock matrix builds

         implicit none

         ! MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: moint_tensor

         ! Auxilliary tensor needed for this algorithm
         type(tens_rcrsv_t):: dm_tensor
         type(tens_rcrsv_t):: a2int_tensor

         integer,intent(in) :: print_level
         integer(INTD)     :: ierr
         integer(INTD)     :: ao_space_id, moint_id(4),a2int_id(4),dm_id(4)
         integer(INTL)     :: ao_space_root,moint_root(4),a2int_root(4),dm_root(4)

         if (print_level >5) call print_date('entered aomo_transform_3')

!        Determine ids and root
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the ht integral tensor
         a2int_id(1:2) = ao_space_id
         a2int_id(3:4) = moint_id(3:4)
         a2int_root(1:2) = ao_space_root
         a2int_root(3:4) = moint_root(3:4)

!        Initialize half-transformed integrals
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(a2int_tensor,'HalfTrans')
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=1; return; endif

         if (print_level >5) call print_date('first ht done')

!        Create the density matrix tensor for the second half-transformation
         dm_id(1:2)    = ao_space_id
         dm_id(3:4)    = moint_id(1:2)
         dm_root(1:2)  = ao_space_root
         dm_root(3:4)   = moint_root(1:2)
         ierr=exatns_tensor_create(dm_tensor,'dm12tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,dm_initlabel(dm_id(1:2),dm_id(3:4)))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         if (print_level >5) call print_date('initialized density matrix')

!        Do 2nd half-transformation
         ierr=exatns_tensor_contract(moint_tensor,a2int_tensor,dm_tensor,"D(i,j,k,l)+=L(p,q,k,l)*R(p,q,i,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=2; return; endif

         if (print_level >5) call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif
!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         if (print_level >5) call print_date('leaving aomo_transform_3')

        end subroutine aomo_transform_3

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,print_level,algorithm)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, integral_tensor
         type(tens_rcrsv_t)                 :: moint_tensor

         integer, intent(in)                :: antisymmetrize
         integer, optional                  :: algorithm
         integer, intent(in)                :: print_level

         integer(INTD)     :: ierr, tens_rank
         integer(INTD)     :: aoint_id(4),integral_id(4),moint_id(4)
         integer(INTL)     :: aoint_root(4),integral_root(4),moint_root(4)

         if (.not.present(algorithm)) algorithm = 1

!        Check integrity of tensors
         tens_rank = aoint_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('input tensor corrupted in aomo_transform')
         tens_rank = integral_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('output tensor corrupted in aomo_transform')
         
!        Determine ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call integral_tensor%get_space_ids(integral_id,integral_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('integral tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the auxilliary ones
!        Use relation (12|34) = <13|24>
         moint_id(1) = integral_id(1)
         moint_id(2) = integral_id(3)
         moint_id(3) = integral_id(2)
         moint_id(4) = integral_id(4)

         moint_root(1) = integral_root(1)
         moint_root(2) = integral_root(3)
         moint_root(3) = integral_root(2)
         moint_root(4) = integral_root(4)

!        Allocate auxilliary tensor to get integrals (without anti-symmetrization) in Mulliken order
         ierr=exatns_tensor_create(moint_tensor,'mulliken',moint_id,moint_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(moint_tensor,ZERO)
        
!        Transform integrals
         select case (algorithm)
         case (1)
            call aomo_transform_1 (aoint_tensor,moint_tensor,print_level)
         case (2)
            call aomo_transform_2 (aoint_tensor,moint_tensor,print_level)
         case default
            call quit ("asked for non-existing algorithm in get_integral_tensor")
         end select

!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=exatns_tensor_destroy(moint_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine get_integral_tensor

        subroutine get_ht2_integral_tensor (a2int_tensor,integral_tensor,antisymmetrize,print_level)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: a2int_tensor, integral_tensor
         type(tens_rcrsv_t)                 :: moint_tensor

         integer, intent(in)                :: antisymmetrize, print_level

         integer(INTD)     :: ierr, tens_rank
         integer(INTD)     :: a2int_id(4),integral_id(4),moint_id(4)
         integer(INTL)     :: a2int_root(4),integral_root(4),moint_root(4)

!        Check integrity of tensors
         tens_rank = a2int_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('input tensor corrupted in aomo_transform')
         tens_rank = integral_tensor%get_rank(ierr)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('output tensor corrupted in aomo_transform')

!        Determine ids and root
         call a2int_tensor%get_space_ids(a2int_id,a2int_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('a2int tensor corrupted')
         call integral_tensor%get_space_ids(integral_id,integral_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('integral tensor corrupted')

!        Copy the information about spaces from the starting and target tensors to the auxilliary ones
!        Use relation (12|34) = <13|24>
         moint_id(1) = integral_id(1)
         moint_id(2) = integral_id(3)
         moint_id(3) = integral_id(2)
         moint_id(4) = integral_id(4)

         moint_root(1) = integral_root(1)
         moint_root(2) = integral_root(3)
         moint_root(3) = integral_root(2)
         moint_root(4) = integral_root(4)

!        Allocate auxilliary tensor to get integrals (without anti-symmetrization) in Mulliken order
         ierr=exatns_tensor_create(moint_tensor,'mulliken',moint_id,moint_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
         ierr=exatns_tensor_init(moint_tensor,ZERO)

!        Transform integrals
         call ao2mo_exat (a2int_tensor,moint_tensor,print_level)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=exatns_tensor_destroy(moint_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine get_ht2_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         integer, intent(in)                :: antisymmetrize


         integer(INTD)   :: mulliken_id(4),dirac_id(4)
         integer(INTL)   :: mulliken_root(4),dirac_root(4)

         integer(INTD)   :: ierr, mismatch

         type(tens_rcrsv_t) :: one_tensor

!        get dimensions of the input tensors and check consistency
         call mulliken_tensor%get_space_ids(mulliken_id,mulliken_root,ierr)
         call dirac_tensor%get_space_ids(dirac_id,dirac_root,ierr)
         mismatch = 0
         if (mulliken_id(1) .ne. dirac_id(1)) mismatch = mismatch + 1
         if (mulliken_id(2) .ne. dirac_id(3)) mismatch = mismatch + 1
         if (mulliken_id(3) .ne. dirac_id(2)) mismatch = mismatch + 1
         if (mulliken_id(4) .ne. dirac_id(4)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 12 .and. dirac_id(1) .ne. dirac_id(2)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 34 .and. dirac_id(3) .ne. dirac_id(4)) mismatch = mismatch + 1
         if (mismatch.gt.0) stop 'spaces mismatch in mulliken_to_dirac_sort'

         ! We first need to define a scalar as a tensor to work around the missing tensor addition
         ierr=exatns_tensor_create(one_tensor,'one_tensor',EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(one_tens) failed!')
         ierr=exatns_tensor_init(one_tensor,ONE)

         ! Now do the actual sort, note that we do not initialize as we may be accumulating in the target tensor
         ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(p,r,q,s)*C()")

         if (antisymmetrize .eq. 12) then
             ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(q,r,p,s)*C()",prefactor=MINUS_ONE)
         elseif (antisymmetrize .eq. 34) then
             ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(p,s,q,r)*C()",prefactor=MINUS_ONE)
         end if

         ierr=exatns_tensor_destroy(one_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine mulliken_to_dirac_sort

        function print_tensornorm2 (tens) result(scalar_real)

          type(tens_rcrsv_t):: tens, scalar
          integer(INTD) :: ierr,rank
          complex(8) :: scalar_complex
          real(8)    :: scalar_real

          ! create a rank-0 tensor
          ierr=exatns_tensor_create(scalar,'scalar',EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')

          ! get the rank and do the contraction
          rank = tens%get_rank(ierr)
          select case (rank)
          case(0)
              ierr=exatns_tensor_contract(scalar,tens,tens,'E()+=D+()*D()')
          case(1)
              ierr=exatns_tensor_contract(scalar,tens,tens,'E()+=D+(a)*D(a)')
          case(2)
              ierr=exatns_tensor_contract(scalar,tens,tens,'E()+=D+(a,b)*D(a,b)')
          case(3)
              ierr=exatns_tensor_contract(scalar,tens,tens,'E()+=D+(a,b,c)*D(a,b,c)')
          case(4)
              ierr=exatns_tensor_contract(scalar,tens,tens,'E()+=D+(a,b,c,d)*D(a,b,c,d)')
          case default
              ierr=1
          end select
          if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_contract() failed!')
          ierr=exatns_tensor_get_scalar(scalar,scalar_complex)
          scalar_real=real(scalar_complex,8)
          ierr=exatns_tensor_destroy(scalar); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif

        end function print_tensornorm2
 
        function mocoef_initlabel (ao_space_id,mo_space_id,spin) result(initlabel)

         !Define a unique name for the initializer of mo coeffiecients for a particular spin
         !we assume that we can use a single digit to indicate the ao and mo spaces
         !if we start using more spaces, we can simply define a longer string here

         integer(INTD),intent(in):: ao_space_id, mo_space_id, spin
         character(len=9)  :: initlabel
         character(len=2)  :: vs_str

         if (ao_space_id>9.or.mo_space_id>9) then
             call quit ('Revise naming scheme for MO coef initializers')
         else
              write (vs_str,'(2I1)') ao_space_id,mo_space_id
         end if
         if (spin == 1) initlabel = 'MOCoef' // vs_str // 'A'
         if (spin == 2) initlabel = 'MOCoef' // vs_str // 'B'

        end function mocoef_initlabel

        function dm_initlabel (ao_space_id,mo_space_id) result(initlabel)

         !Define a unique name for the initializer of density matrices
         !we assume that we can use a single digit to indicate the ao and mo spaces
         !if we start using more spaces, we can simply define a longer string here

         integer(INTD),intent(in):: ao_space_id(2), mo_space_id(2)
         character(len=9)        :: initlabel
         character(len=2)        :: vs_str(2)
         integer                 :: i

         do i = 1, 2
            if (ao_space_id(i)>9.or.mo_space_id(i)>9) then
               call quit ('Revise naming scheme for DM initializers')
            else
              write (vs_str(i),'(2I1)') ao_space_id(i),mo_space_id(i)
           end if
         end do
         initlabel = 'DensM' // vs_str(1) // vs_str(2)

        end function dm_initlabel

      subroutine ao2mo_exat(ao_tensor,mo_tensor,print_level)

!       Transform integral tensor from ao-basis to mo-basis

        implicit none

        type(tens_rcrsv_t), intent(inout) :: ao_tensor, mo_tensor 
        integer,intent(in)                :: print_level

        ! Auxilliary tensors that are created and destroyed inside this subroutine
        type(tens_rcrsv_t):: coef_tensor
        type(tens_rcrsv_t):: t_tensor

        integer           :: i, ind_a, ind_b
        integer(INTD)     :: ierr,spin, aorank, morank
        integer(INTD)     :: ao4_id(4),ao3_id(3),ao2_id(2)
        integer(INTD)     :: mo4_id(4),mo3_id(3),mo2_id(2)
        integer(INTD)     :: t4_id(4), t3_id(3), t2_id(2)
        integer(INTD)     :: coef_id(2)
        integer(INTL)     :: ao4_root(4),ao3_root(3),ao2_root(2)
        integer(INTL)     :: mo4_root(4),mo3_root(3),mo2_root(2)
        integer(INTL)     :: t4_root(4), t3_root(3), t2_root(2)
        integer(INTL)     :: coef_root(2)
        character(len=9)  :: coef_name
        character(len=38) :: info_str

        if (print_level >5) call print_date('=== entered ao2mo_exat')

!       Determine ranks
        aorank = ao_tensor%get_rank(ierr)
        if (ierr.ne.0) call quit('ao2mo_exat: ao tensor corrupted (rank)')
        morank = mo_tensor%get_rank(ierr)
        if (ierr.ne.0 .or. aorank.ne.morank) call quit('ao2mo_exat: mo tensor corrupted (rank) ')

        if (aorank.eq.2) then
          call ao_tensor%get_space_ids(ao2_id,ao2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (2) ')
          call mo_tensor%get_space_ids(mo2_id,mo2_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: mo tensor corrupted (2)')

          !create temporary tensor
          t2_id(1)   = mo2_id(1)
          t2_id(2)   = ao2_id(2)
          t2_root(1) = mo2_root(1)
          t2_root(2) = ao2_root(2)

          ierr=exatns_tensor_create(t_tensor,'temp',t2_id,t2_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp2 tensor not created')

          do spin = 1, 2
            ierr=exatns_tensor_init(t_tensor,ZERO)

            i=1 ! first index
            coef_id(1)   = ao2_id(i)
            coef_root(1) = ao2_root(i)
            coef_id(2)   = mo2_id(i)
            coef_root(2) = mo2_root(i)

            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            ierr=exatns_tensor_contract(t_tensor,ao_tensor,coef_tensor,"B(p,s)+=A(r,s)*M+(r,p)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 1 failed')

            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'half(2) transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if

            i=2  ! second index
            coef_id(1)   = ao2_id(i)
            coef_root(1) = ao2_root(i)
            coef_id(2)   = mo2_id(i)
            coef_root(2) = mo2_root(i)

            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            ierr=exatns_tensor_contract(mo_tensor,t_tensor,coef_tensor,"C(p,q)+=B(p,s)*M(s,q)")
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 2 failed')

            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'quarter transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if

          end do

          ierr=exatns_tensor_destroy(t_tensor)

        else if (aorank.eq.3) then
          call quit('ao2mo_exat: for now rank 3 not implemented')

        else if (aorank.eq.4) then
          ! Algorithm 1: n^5 scaling with index-by-index transform
          call ao_tensor%get_space_ids(ao4_id,ao4_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: ao tensor corrupted (4)')
          call mo_tensor%get_space_ids(mo4_id,mo4_root,ierr)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: mo tensor corrupted (4)')

          ! create temporary tensor
          if (ao4_id(1).eq.ao4_id(3) .and. ao4_id(2).eq.ao4_id(4)) then
            !first halftransformation
            t4_id(1)     = mo4_id(1)
            t4_id(2:4)   = ao4_id(2:4)
            t4_root(1)   = mo4_root(1)
            t4_root(2:4) = ao4_root(2:4)
            ind_a=1
            ind_b=2
          else
            t4_id(1:3)   = mo4_id(1:3)
            t4_id(4)     = ao4_id(4)
            t4_root(1:3) = mo4_root(1:3)
            t4_root(4)   = ao4_root(4)
            ind_a=3
            ind_b=4
          end if

          ierr=exatns_tensor_create(t_tensor,'temp',t4_id,t4_root,EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: temp4 tensor not created')

          do spin = 1, 2
            ierr=exatns_tensor_init(t_tensor,ZERO)

            i=ind_a ! first index
            coef_id(1)   = ao4_id(i)
            coef_root(1) = ao4_root(i)
            coef_id(2)   = mo4_id(i)
            coef_root(2) = mo4_root(i)
          
            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            if (ind_a.eq.1) then
              ierr=exatns_tensor_contract(t_tensor,coef_tensor,ao_tensor,"D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 1 failed')
            else
              ierr=exatns_tensor_contract(t_tensor,coef_tensor,ao_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 3 failed')
            end if
            
            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'quarter transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if

            i=ind_b  ! second index
            coef_id(1)   = ao4_id(i)
            coef_root(1) = ao4_root(i)
            coef_id(2)   = mo4_id(i)
            coef_root(2) = mo4_root(i)

            write (coef_name,'(A7,2I1)') 'mocoef_',spin,i
            ierr=exatns_tensor_create(coef_tensor,coef_name,coef_id,coef_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not created')
            ierr=exatns_tensor_init(coef_tensor,mocoef_initlabel(coef_id(1),coef_id(2),spin))
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: coef_tensor not set')
            if (print_level>10) write(*,*) "created and initialized coef_tensor ",i,spin, &
              mocoef_initlabel(coef_id(1),coef_id(2),spin),print_tensornorm2(coef_tensor)

            if (ind_b.eq.2) then
              ierr=exatns_tensor_contract(mo_tensor,t_tensor,coef_tensor,"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 2 failed')
            else
              ierr=exatns_tensor_contract(mo_tensor,t_tensor,coef_tensor,"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
              if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: contraction 4 failed')
            end if

            ierr=exatns_tensor_destroy(coef_tensor)
            if(ierr.ne.EXA_SUCCESS) call quit('ao2mo_exat: tensor not removed')

            if (print_level >5) then
              write (info_str,'(A31,I2,A3,I2)') 'quarter transformation done: s=',spin,',i=',i
              call print_date(info_str)
            end if
          end do
          
          ierr=exatns_tensor_destroy(t_tensor)

        else
          call quit('ao2mo_exat: wrong rank or mode')
        end if

        if (print_level >5) call print_date('=== leaving ao2mo_exat')

        end subroutine ao2mo_exat

! stub non-MPI routine to satisfy references

#else

       contains

        subroutine exacorr_compute_mo_integrals (nmo, mo_list)
         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         print*,  'called exatensor version of compute_mo_integrals in a serial version'
         print*,  'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_compute_mo_integrals

#endif

end module exacorr_ao_to_mo
