module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exatensor

        implicit none
        private
        public exacorr_compute_mo_integrals
#if defined (VAR_MPI)
        public get_integral_tensor
        public print_tensornorm2
#endif
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         ! It is not possible to make an array of pointers, so we need a container
         type,public ::  h_spaces_container_t
              class(h_space_t), pointer :: ptr
         end type h_spaces_container_t

       contains

#if defined (VAR_MPI)
        subroutine exacorr_compute_mo_integrals (nmo, mo_list)

!        This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!        It currently tests two different algorithms as well as the get_integral_tensor routine used by the CC code.
!        Intended final purpose is two run the most efficient algorithm only and write integrals out to a file that can be used by correlation methods.
!        See talsh_compute_mo_integrals to see how this works for the single node implementation.

!        Written by Lucas Visscher, march/september 2018

         use interface_to_mpi
         use dirac_MO
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use exatensor
         use, intrinsic:: ISO_C_BINDING

         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer:: my_role,mpi_th_provided, my_MPI_rank ! needed to start exatensor, MPI should already be started in threaded mode
         integer             :: nshells                 ! the number of shells (orbitals with the same center and l-value)
         integer             :: nao                     ! the number of atomic orbitals)
         type(basis_func_info_t), allocatable :: gto(:) ! array with information about the shells

         ! classes for the definition of spaces (AO and 4 different MO spaces)
         class(h_space_t), pointer                :: ao_space
         type(h_spaces_container_t), dimension(4) :: mo_space

         ! variables for the definition of spaces (AO and 4 different MO spaces)
         type(subspace_basis_t):: basis_ao, basis_mo(4)
         integer(INTL)::          nao_l, ao_space_root
         integer(INTL)::          nmo_l, mo_space_root(4)
         integer::                ao_space_id, mo_space_id(4)

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL)::                        labs
         integer::                              l_col=1

         ! types of the initializer methods
         type(compute_2e_ao_tensor_t) :: aoint_calculator ! method to initialize the ao-integrals tensor
         type(init_mocoef_tensor_t)   :: mocoef_init(2,4) ! array of initializers for algorithm 1, first index is alpha/beta, second is 1,2,3,4
         type(init_dm_tensor_t)       :: dm_init(2)       ! array of initializers for algorithm 2, index indicates first/second halftransformation
         type(tens_printer_t)         :: tens_printer     ! generic method, but needs to be defined here in specific form
         type(scale_with_denom_t)     :: denominate

         ! variables for the definition of the initializer methods
         character(len=9)  :: initlabel(2,4),initlabel2(2)
         character(len=7)  :: MOCoef_str = "MOCoef_"
         character(len=7)  :: DM_str     = "DM_mat_"
         character(len=2)  :: vs_str

         ! the two tensors that we need
         type(tens_rcrsv_t):: aoint_tensor
         type(tens_rcrsv_t):: moint_tensor

         ! loop variables
         integer:: i, ierr
         integer:: spin,imo,jmo,imo12,jmo12,imo34,jmo34
         integer:: lsh, l, nl

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) print*, " Starting AO-to-MO transform with exatensor"

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

         !Get the MO-coeffcient basis for the four indices that are to be transformed
         ierr = 0
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            nmo_l   = 2*int(nmo(i),8)  ! NB: imo, jmo and nmo count orbitals, nmo_l counts spinors
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
            write (vs_str,'(A,I1)') '_',i
            ierr=exatns_space_register('MOspace'//vs_str,basis_mo(i),mo_space_id(i),mo_space(i)%ptr)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

            imo = imo + nmo(i) ! make imo point to the start of the list for the next index
         end do

         !Define separate initializers for the alpha(A) and beta(B) MO-coefficients for each index 
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            do spin = 1, 2
               !Define a unique name for the initializer for this spin (A/B) of space i
               if (spin == 1) write (vs_str,'(I1,A)') i,'A'
               if (spin == 2) write (vs_str,'(I1,A)') i,'B'
               initlabel(spin,i) = MOCoef_str // vs_str
               !Register this initializer as a method
               call mocoef_init(spin,i)%init_mocoef_ctor(mo_list(imo:jmo),nmo(i),spin)
               ierr=exatns_method_register(initlabel(spin,i),mocoef_init(spin,i))
               if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            end do
         end do

         !Define separate initializers for the first and second density matrix
         imo12 = 1
         jmo12 = nmo(1)+nmo(2)
         imo34 = jmo12 + 1
         jmo34 = jmo12 + nmo(3)+nmo(4)
         !Define names for the initializers for density matrices
         initlabel2(1) = DM_str // '12'
         initlabel2(2) = DM_str // '34'
         !Register the dm initializers as methods
         call dm_init(1)%init_dm_ctor(mo_list(imo12:jmo12),nmo(1:2))
         ierr=exatns_method_register(initlabel2(1),dm_init(1))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
         call dm_init(2)%init_dm_ctor(mo_list(imo34:jmo34),nmo(3:4))
         ierr=exatns_method_register(initlabel2(2),dm_init(2))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         if (my_MPI_rank == 0) print*," registered MO vectors spaces and initializers with nmo:",nmo

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the denominator initializer aa a method (testing only, registration will move to exacorr_cc)
         call denominate%scale_with_denom_ctor(mo_list(1:2),2,mo_list(3:nmo(1)),nmo(1)-2)
         ierr=exatns_method_register('Denominate',denominate)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
   
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

           call print_date('starting ao integral calculation')

           !CALC GOES WRONG IN THIS PART
           !Create ao-integral-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           call print_date('finished ao integral calculation')

           ! the  next three lines are for debugging, 
           call print_tensornorm2 (aoint_tensor) ! this is usually sufficent for checking
           ! when uncommenting the next line: keep in mind that printing an entire tensor yields an enormous output
           ! ierr=exatns_tensor_traverse(aoint_tensor,'PrintTensor')

           !Create mo-integral-tensor:
           do i = 1, 4
              mo_space_root(i)=mo_space(i)%ptr%get_root_id(ierr)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           call print_date('created mo integral tensor')

           call aomo_transform_1 (aoint_tensor,moint_tensor,initlabel)
           call print_date('completed mo integral tensor: algorithm 1')
           call print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_init(moint_tensor,zero)
           call aomo_transform_2 (aoint_tensor,moint_tensor,initlabel2)
           call print_date('completed mo integral tensor: algorithm 2')
           call print_tensornorm2 (moint_tensor)

           call get_integral_tensor (aoint_tensor,moint_tensor,0,initlabel)
           call print_date('completed dirac-sorted mo integral tensor - not anti-symmetrized')
           call print_tensornorm2 (moint_tensor)

           ierr=exatns_tensor_destroy(moint_tensor)

           !Recreate mo-integral-tensor with equal index spaces so that we can anti-symmetrize:
           do i = 2, 4
              mo_space_id(i)=mo_space_id(1)
              mo_space_root(i)=mo_space_root(1)
              initlabel(1:2,i)=initlabel(1:2,1)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')

           call get_integral_tensor (aoint_tensor,moint_tensor,12,initlabel)
           call print_date('completed dirac-sorted mo integral tensor - 12 anti-symmetrized')
           call print_tensornorm2 (moint_tensor)

           call get_integral_tensor (aoint_tensor,moint_tensor,34,initlabel)
           call print_date('completed dirac-sorted mo integral tensor - 34 anti-symmetrized')
           call print_tensornorm2 (moint_tensor)

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

        subroutine aomo_transform_1 (aoint_tensor,moint_tensor,initlabel)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 1: n^5 scaling with index-by-index transform

         implicit none

         ! The strings initlabel are used to find the appropriate initializer for the mo-coefficients
         ! This initializer has the coefficients stored (on all nodes) and just picks out the proper subset.
         character(len=*), intent(in)       :: initlabel(2,4)
         
         ! AO integrals are to be provided on input and are not changed, MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor   ! for technical reasons intent(inout) needs to be used for all tensors

         ! Auxilliary tensors that are created and destroyed inside this subroutine
         type(tens_rcrsv_t):: mocoef_tensor(2,4)
         type(tens_rcrsv_t):: a1int_tensor,a2int_tensor, a3int_tensor

         integer           :: ierr,i,spin
         integer           :: aoint_id(4),moint_id(4),a1int_id(4),a2int_id(4),a3int_id(4),mocoef_id(2)
         integer(INTL)     :: aoint_root(4),moint_root(4),a1int_root(4),a2int_root(4),a3int_root(4),mocoef_root(2)
         character(len=9)  :: mocoef_name

         call print_date('entered aomo_transform_1')

!        Determine space ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Initialize tensors with MO coefficients, one for each index and alpha/beta
         mocoef_id(1)   = aoint_id(1)
         mocoef_root(1) = aoint_root(1)
         do i = 1, 4
           mocoef_id(2)   = moint_id(i)
           mocoef_root(2) = moint_root(i)
           do spin = 1, 2
!             The name of the tensor needs to be unique, define a name that makes sense
              write (mocoef_name,'(A7,2I1)') 'mocoef_',spin,i
!             Create mo-coefficient-tensor
              ierr=exatns_tensor_create(mocoef_tensor(spin,i),mocoef_name, &
               mocoef_id,mocoef_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
!             Initialize mocoef-tensor 
              ierr=exatns_tensor_init(mocoef_tensor(spin,i),initlabel(spin,i))
              if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')
!             next two lines can be used for debugging
!             print*,"created and initialized mocoef_tensor ",spin,i
!             call print_tensornorm2 (mocoef_tensor(spin,i))
           end do
         end do
!        next two lines can be used for debugging
!        ierr=exatns_tensor_traverse(mocoef_tensor(1,1),'PrintTensor')
!        ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

         call print_date('initialized mo-coefficients')

!        Copy the information about spaces from the starting and target tensors to the auxilliary ones
         a1int_id(1)   = moint_id(1)
         a1int_id(2:4) = aoint_id(2:4)
         a2int_id(1:2) = moint_id(1:2)
         a2int_id(3:4) = aoint_id(3:4)
         a3int_id(1:3) = moint_id(1:3)
         a3int_id(4)   = aoint_id(4)

         a1int_root(1)   = moint_root(1)
         a1int_root(2:4) = aoint_root(2:4)
         a2int_root(1:2) = moint_root(1:2)
         a2int_root(3:4) = aoint_root(3:4)
         a3int_root(1:3) = moint_root(1:3)
         a3int_root(4)   = aoint_root(4)

!        Do 1st and 2nd quarter-transformation
         ierr=exatns_tensor_create(a1int_tensor,'a1int',a1int_id,a1int_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_contract(a1int_tensor,mocoef_tensor(1,1),aoint_tensor,"D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=exatns_tensor_contract(a2int_tensor,a1int_tensor,mocoef_tensor(1,2),"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         ierr=exatns_tensor_init(a1int_tensor,zero) ! alpha-alpha part done, continue with beta-beta part
         ierr=exatns_tensor_contract(a1int_tensor,mocoef_tensor(2,1),aoint_tensor,"D(i,q,r,s)=R+(p,i)*L(p,q,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=exatns_tensor_contract(a2int_tensor,a1int_tensor,mocoef_tensor(2,2),"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=22; return; endif

         call print_date('first ht done')

!        Destruct a1 integrals:
         ierr=exatns_tensor_destroy(a1int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=25; return; endif

!        Do 3rd and 4th quarter-transformation
         ierr=exatns_tensor_create(a3int_tensor,'a3int',a3int_id,a3int_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_init(a3int_tensor,zero)
         ierr=exatns_tensor_contract(a3int_tensor,mocoef_tensor(1,3),a2int_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=exatns_tensor_contract(moint_tensor,a3int_tensor,mocoef_tensor(1,4),"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=exatns_tensor_init(a3int_tensor,zero) ! alpha-alpha part done, continue with beta-beta part
         ierr=exatns_tensor_contract(a3int_tensor,mocoef_tensor(2,3),a2int_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=32; return; endif
         ierr=exatns_tensor_contract(moint_tensor,a3int_tensor,mocoef_tensor(2,4),"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=42; return; endif

         call print_date('second ht done')

!        Destruct a2 and a3 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=35; return; endif
         ierr=exatns_tensor_destroy(a3int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif
!        Destruct mocoef tensors:
         do i = 1, 4
           do spin = 1, 2
            ierr=exatns_tensor_destroy(mocoef_tensor(spin,i)); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif
            end do
         end do

         call print_date('leaving aomo_transform_1')

        end subroutine aomo_transform_1

        subroutine aomo_transform_2 (aoint_tensor,moint_tensor,initlabel2)

!        Transform integral tensor from ao-basis to mo-basis
!        Algorithm 1: n^6 scaling with density matrices

         implicit none

         ! The strings in initlabel2 are used to find the appropriate initializer for the density matrices.
         ! This initializer has the MO-coefficients stored (on all nodes) and just picks out the proper subsets.
         character(len=*), intent(in)       :: initlabel2(2)

         ! AO integrals are to be provided on input and are not changed, MO integrals are returned on output
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor

         ! Auxilliary tensor needed for this algorithm
         type(tens_rcrsv_t):: dm_tensor
         type(tens_rcrsv_t):: a2int_tensor

         integer           :: ierr,i
         integer           :: aoint_id(4),moint_id(4),a2int_id(4),dm_id(4)
         integer(INTL)     :: aoint_root(4),moint_root(4),a2int_root(4),dm_root(4)

         call print_date('entered aomo_transform_2')

!        Determine ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Create the density matrix tensor for the first half-transformation
         dm_id(1:2)    = aoint_id(1:2)
         dm_id(3:4)    = moint_id(1:2)
         dm_root(1:2)  = aoint_root(1:2)
         dm_root(3:4)  = moint_root(1:2)
         ierr=exatns_tensor_create(dm_tensor,'dm12tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm1) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,initlabel2(1))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         call print_date('initialized first density matrix')

!        Copy the information about spaces from the starting and target tensors to the ht integral tensor
         a2int_id(1:2) = moint_id(1:2)
         a2int_id(3:4) = aoint_id(3:4)
         a2int_root(1:2) = moint_root(1:2)
         a2int_root(3:4) = aoint_root(3:4)

!        Do 1st half-transformation
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         call print_date('constructed a2')
         ierr=exatns_tensor_contract(a2int_tensor,aoint_tensor,dm_tensor,"D(i,j,r,s)+=L(p,q,r,s)*R(p,q,i,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=1; return; endif

!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         call print_date('first ht done')

!        Create the density matrix tensor for the second half-transformation
         dm_id(1:2)    = aoint_id(3:4)
         dm_id(3:4)    = moint_id(3:4)
         dm_root(1:2)  = aoint_root(3:4)
         dm_root(3:4)   = moint_root(3:4)
         ierr=exatns_tensor_create(dm_tensor,'dm34tensor',dm_id,dm_root,EXA_DATA_KIND_C8)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create(dm) failed!')
!        Initialize dm-tensor 
         ierr=exatns_tensor_init(dm_tensor,initlabel2(2))
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')

         call print_date('initialized second density matrix')

!        Do 2nd half-transformation
         ierr=exatns_tensor_contract(moint_tensor,a2int_tensor,dm_tensor,"D(i,j,k,l)+=L(i,j,r,s)*R(r,s,k,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=2; return; endif

         call print_date('second ht done')

!        Destruct a2 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif
!        Destruct dm tensor:
         ierr=exatns_tensor_destroy(dm_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=3; return; endif

         call print_date('leaving aomo_transform_2')

        end subroutine aomo_transform_2

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,initlabel)

         implicit none

         ! The strings initlabel are used to find the appropriate initializer for the mo-coefficients
         ! This initializer has the coefficients stored (on all nodes) and just picks out the proper subset.
         character(len=*), intent(in)       :: initlabel(2,4) ! Note that these are already given in Mulliken order !
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, integral_tensor
         type(tens_rcrsv_t)                 :: moint_tensor

         integer, intent(in)                :: antisymmetrize

         integer           :: ierr, tens_rank
         integer           :: aoint_id(4),integral_id(4),moint_id(4)
         integer(INTL)     :: aoint_root(4),integral_root(4),moint_root(4)

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

!        Transform integrals
         call aomo_transform_1 (aoint_tensor,moint_tensor,initlabel)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=exatns_tensor_destroy(moint_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine get_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         implicit none

         type(tens_rcrsv_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         integer, intent(in)                :: antisymmetrize


         integer           :: mulliken_id(4),dirac_id(4)
         integer(INTL)     :: mulliken_root(4),dirac_root(4)

         integer   :: ierr, tens_rank, mismatch

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

         ! Now do the actual sort
         ierr=exatns_tensor_init(dirac_tensor,ZERO)
         ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(p,r,q,s)*C()")

         if (antisymmetrize .eq. 12) then
             ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(q,r,p,s)*C()",prefactor=MINUS_ONE)
         elseif (antisymmetrize .eq. 34) then
             ierr=exatns_tensor_contract(dirac_tensor,mulliken_tensor,one_tensor,"D(p,q,r,s)+=M(p,s,q,r)*C()",prefactor=MINUS_ONE)
         end if

         ierr=exatns_tensor_destroy(one_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=1; return; endif

        end subroutine mulliken_to_dirac_sort

        subroutine print_tensornorm2 (tens)

          type(tens_rcrsv_t):: tens, scalar
          integer ierr,rank
 
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
          ierr=exatns_tensor_traverse(scalar,'PrintTensor')
          if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_traverse() failed!')
          ierr=exatns_tensor_destroy(scalar); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif

        end subroutine print_tensornorm2

! stub non-MPI routine to satisfy references

#else
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
