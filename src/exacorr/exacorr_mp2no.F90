module exacorr_mp2no

!This module contains the routines needed to calculate a coupled cluster wave function
!almost everything requires MPI, but we provide a stub function to make compilation of 
!non-MPI versions of DIRAC possible.

#if defined (VAR_MPI)

        use exatensor
        use exacorr_datatypes
        use exacorr_utils
        use intermediates
        use exacorr_ao_to_mo,  only : print_tensornorm2
        use, intrinsic:: ISO_C_BINDING

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0), &
                                 ONE_QUARTER=(0.25D0,0.D0), MINUS_ONE_QUARTER=(-0.25D0,0.D0),  &
                                 MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0), &
                                 THREE=(3.D0,0.D0),SIX=(6.D0,0.D0)
        type(tens_rcrsv_t),public :: one_tensor
        real(8) :: e_core
        logical :: one_el_exist

        private

        public exacorr_mp2no_driver


       contains

!-----------------------------------------------------------------------------------------------------
      subroutine exacorr_mp2no_driver (exa_input,N_cano)

!       This module is used to construct the MP2 density matrix by talsh.
!       The complex density matrix are saved on the file 'DM_complex'.  
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new reduced canonical orbitals are written into 'MP2NOs_AO' by Module3.
!       Written by Xiang Yuan
            use exacorr_global
            use interface_to_mpi

            type(exacc_input), intent(in) :: exa_input
            integer, intent(inout)        :: N_cano                  ! number of new canonical orbital
            integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode
    
    
    !       Point connected with tensors 
            integer                :: ierr
            real(8),pointer        :: NOs_AO(:,:,:)           ! new quaternion orbital in AO basis
            real(kind=8), pointer  :: N_OCCUPATION(:)         ! occupation number from qdia of DM
            real(kind=8), pointer  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
            logical                :: alive
            real(kind=8)           :: threshold
            integer                :: i
    
            
            threshold = exa_input%t_mp2no
    
            inquire(file="NOs_MO", exist=alive)
            if (alive) then
                call print_date('File NOs_MO existed')
            else
                call print_date('File NOs_MO is not existed. Start constructing MP2 density matrix')
                call exacorr_generate_MP2_VNO_MObasis(exa_input)   
                call get_transform_Fock(exa_input,'DM_complex')
            end if
            call get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION, &
                                  New_orbital_energy,'NOs_MO',N_cano)
            call write_to_newfile(exa_input,NOs_AO,N_OCCUPATION,New_orbital_energy,'DFCOEF','MP2NOs_AO')
    
            call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
            if (my_MPI_rank == 0) call print_date(" Leave MP2_NO Module")

            ierr=exatns_method_unregister('ZERO')
            print *, "unregister method ierr", ierr
            ierr=exatns_method_unregister('Aoint_calculator_mp2')
            print *, "unregister method ierr Aoint_calculator", ierr
      !      ierr=exatns_method_unregister('AointSegment_mp2')
      !      print *, "unregister method ierr AointSegment", ierr
            ierr=exatns_method_unregister('Denominate_mp2')
            print *, "unregister method ierr Denominate", ierr
      !      ierr=exatns_method_unregister('PrintTensor_mp2')
      !      print *, "unregister method ierr PrintTensor", ierr
            ierr=exatns_method_unregister('Unity_mp2')
            print *, "unregister method ierr Unity", ierr

            


      end subroutine exacorr_mp2no_driver



!-------------------------------------------------------------------------------------------
        subroutine exacorr_generate_MP2_VNO_MObasis (exa_input)

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
!         integer, intent(inout)        :: N_cano                  ! number of new canonical orbital

         integer(INTD) :: my_role
         integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode

         integer(INTD)                  :: ierr
         complex(8)                     :: result_val

!        Common tensors
         type(exatns_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(exatns_intg_tens) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t) :: l1_tensor, l2_tensor
!        scalars (need to be defined as tensor types)
         type(tens_rcrsv_t)         :: result_tensor
         type(exatns_dm_tens)       :: MP2_density_tensor
!        space ids and roots
         type(space_dims) :: dims
!        CCSD control and result variables
         real(8)            :: scf_energy, mp2_energy, cc_energy, t_energy(3) ! triples not implemented
         real(8)            :: t1diag
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


         !methods
         type(delta_t)      :: f_delta
         type(denom_t)      :: denom
         type(set_energy_t) :: orb_set
         type(denom3_t)     :: denom3
         type(chol_J)       :: j_chol
         type(set_ff_t)     :: ff_set

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         if (my_MPI_rank == 0) call print_date(" Starting MP2NO calculation with exatensor")

!        Initialize and register AO and MO spaces and get their IDs for further reference
         call exacorr_mp2no_spaces (exa_input,ao_space,occ_space,vir_space,ao_space_id,occ_space_id,vir_space_id)
         if (my_MPI_rank == 0) call print_date(" Registered AO and MO spaces")

!        Register and initialize user defined methods needed in the CC calculation

         call exacorr_mp2no_methods (exa_input,ao_space_id,occ_space_id,vir_space_id, &
                                  f_delta,orb_set,denom,denom3,j_chol,ff_set)
         if (my_MPI_rank == 0) call print_date(" Initialized and registered tensor methods")

!        copy input variables to local ones for ease of reference
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         t_target_precision = exa_input%t_econv



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
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

!        All preliminary work is now done on all MPI-nodes, we may proceed to start ExaTensor.
!        After the start of ExaTensor, the driver process (the last MPI process) will be executing
!        the commands and tell the others (that stay inside exatns_process_role) what to do.
!        These other nodes return only after exatns_stop() is called by the driver.

!        Start ExaTENSOR within MPI_COMM_WORLD (called global_communicator in the MPI interface that DIRAC uses):
         !call exatns_ctrl_reset_logging(1,2) !logging for debugging purposes
         ierr=exatns_start(int(global_communicator,INT_MPI))

         if(ierr.eq.EXA_SUCCESS) then
            ierr=exatns_process_role(my_role) ! only EXA_DRIVER will return immediately
            if(my_role.eq.EXA_DRIVER) then

!           Create and initialize scalars that are to be used as tensors in contractions
            ierr=exatns_tensor_create(MP2_density_tensor%prime_vv,"gprime_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(result_tensor,"result_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(one_tensor,"another_one_tensor",EXA_DATA_KIND_C8)
            ! Note:  also get_integral_tensors will define a one_tensor and call it "one_tensor", so we needed another name
            ierr=exatns_tensor_init(one_tensor,ONE)
            
!           Get fixed two-body tensors
            call print_date('-Start- Integral Transformation')
            ierr=exatns_tensor_create(int_t%oooo,"oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%ooov,"ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%oovv,"oovv",oovv_id,oovv_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(int_t%vovo,"vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
            call get_MP2_integrals (ao_space_id,ao_space,int_t,exa_input%moint_scheme,exa_input%print_level)
            call print_date('-End- Integral Transformation')

            ! If we have the 1-body integrals we rebuild the Fock matrix and check consistency, otherwise we use the orbital energies
            ierr=exatns_tensor_create(comm_t%foo,"foo",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fov,"fov",ov_id,ov_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(comm_t%fvv,"fvv",vv_id,vv_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            call init_fock_exat (comm_t,int_t,oo_id,oo_root,ov_id,ov_root,vv_id,vv_root, &
                            f_delta,denom,orb_set,ff_set,exa_input,scf_energy)
            call print_date('initialized Fock matrix')


!           Initialize t2 amplitudes with the MP2 values
            ierr=exatns_tensor_create(comm_t%t2,"t2",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            ierr=exatns_tensor_init(comm_t%t2,ZERO)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)*K()",comm_t%t2,int_t%oovv,one_tensor)
            ierr=exatns_tensor_transform(comm_t%t2,denom)

!           Calculate and print the MP2 energy
            ierr=exatns_tensor_init(result_tensor,ZERO)
            ierr=exatns_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv,ONE_QUARTER)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=real(result_val,8)
            ierr=exatns_tensor_init(result_tensor,ZERO)
            ierr=exatns_tensor_contract("E()+=T(a,i)*V(i,a)",result_tensor,comm_t%t1,comm_t%fov)
            ierr=exatns_tensor_get_scalar(result_tensor,result_val)
            mp2_energy=mp2_energy+real(result_val,8)
         
            write(*,*) ""
            write(*,*) "    MP2 energy = ", mp2_energy
            write(*,*) ""


            ierr=exatns_tensor_transform(int_t%oovv,denom)
            ierr=exatns_tensor_init(MP2_density_tensor%prime_vv,ZERO)
            ierr=exatns_tensor_contract("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_density_tensor%prime_vv, &
                                    int_t%oovv,int_t%oovv,ONE_HALF)               
!        Make some noise so that we know we are proceeding 
            
            call write_matrix(exa_input,MP2_density_tensor,'DM_complex')
            call print_date('Save complex non-digonal DM to file DM_complex')

            ierr=exatns_tensor_destroy(int_t%oooo)
            ierr=exatns_tensor_destroy(int_t%ooov)
            ierr=exatns_tensor_destroy(int_t%oovv)
            ierr=exatns_tensor_destroy(int_t%vovo)
            ierr=exatns_tensor_destroy(comm_t%foo)
            ierr=exatns_tensor_destroy(comm_t%fov)
            ierr=exatns_tensor_destroy(comm_t%fvv)
            ierr=exatns_tensor_destroy(comm_t%tau)
            ierr=exatns_tensor_destroy(comm_t%hoo)
            ierr=exatns_tensor_destroy(comm_t%hov)
            ierr=exatns_tensor_destroy(comm_t%hvv)
            ierr=exatns_tensor_destroy(comm_t%goo)
            ierr=exatns_tensor_destroy(comm_t%gvv)
            ierr=exatns_tensor_destroy(comm_t%a_int)
            ierr=exatns_tensor_destroy(comm_t%h_int)
            ierr=exatns_tensor_destroy(result_tensor)
            ierr=exatns_tensor_destroy(MP2_density_tensor%prime_vv)
            ierr=exatns_tensor_destroy(one_tensor)
            ierr=exatns_tensor_destroy(comm_t%t1)
            ierr=exatns_tensor_destroy(comm_t%t2)

            !Stop ExaTENSOR runtime
            

            !ierr=exatns_space_unregister('Virtual')
            !print *, "unregister space ierr Virtual", ierr
            ierr=exatns_stop()

            end if

         else
            write(*,*) ' Process ',my_MPI_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         
         call delete_global_data

!        Make some noise so that we know we are leaving
         if (my_MPI_rank == 0) then
            write(*,*) ""
            write(*,*) " Leaving mp2no DM construction routine"
         end if
         return
                              

        end subroutine exacorr_generate_MP2_VNO_MObasis

!------------------------------------------------------------------------------------------------------------------

        subroutine write_matrix(exa_input,gamma,file_out)
            ! Routine for writing complex density matrix into the file'DM_complex'. 
            
            use exacorr_global,         only : get_nsolutions
            use talsh
            use tensor_algebra
            use talsh_common_routines,  only : print_tensor
#if defined (VAR_MPI)
            use interface_to_mpi
      
      !     variables for mpi
            integer(C_INT)    :: read_node
            integer           :: my_MPI_rank
#endif
            
            ! input variables
                  type(exacc_input),   intent(in   )  :: exa_input
                  type(exatns_dm_tens)                :: gamma
                  character(len=*),intent(in)         :: file_out
            ! talsh tensor to get local copy
                  type(talsh_tens_t)           :: tensor_slice
                  integer(INTD), dimension(2)  :: h_dims
            
            ! matrix dimensions
                  integer :: nocc, nvir
            ! body pointers
                  type(C_PTR)                     :: body_p2
                  complex(8), pointer, contiguous :: dm_block(:,:)
                  complex(8), allocatable         :: dm_mo_fort(:,:)
            ! error code
                  integer :: ierr
            ! loop variables
                  integer :: i, j, k
                  integer :: mp2dm
            ! debug
                  logical :: debug=.false.
            ! variables for talsh, density matrix needs to be smaller than buffer size 
                  integer(C_SIZE_T)               :: buf_size=1024_8*1024_8*1024_8 
                  integer(C_INT)                  :: host_arg_max
            
            
                  read_node=0
            !***************
            !     set dims *
            !***************
                  nocc    = 0
                  nvir    = exa_input%nvir
                     
            !***************************************
            !     copy DM blocks into one large DM *
            !***************************************
                     
                  allocate (dm_mo_fort(nvir,nvir))
                  dm_mo_fort = ZERO
                     
                  buf_size=exa_input%talsh_buff*buf_size
                  ierr=talsh_init(buf_size,host_arg_max)
                  call print_date('Initialized talsh library')
                  write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

                  !-----------
                  !| OO | OV |
                  !-----------
                  !| VO | VV |
                  !-----------
                              
                  !     vv
                  h_dims=nvir
                  ierr=talsh_tensor_construct(tensor_slice,C8,h_dims,init_val=ZERO)
                  ierr=exatns_tensor_get_slice(gamma%prime_vv,tensor_slice)
                  ierr=talsh_tensor_get_body_access(tensor_slice,body_p2,C8,read_node,DEV_HOST)
                  call c_f_pointer(body_p2,dm_block,h_dims) ! to use <dm_mo_block> as a regular Fortran 2d array
                  dm_mo_fort(1:nvir,1:nvir) = dm_block

                  
                  call get_free_fileunit(mp2dm)
                  open (mp2dm,file=file_out,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
                  ! write dimension of the matrix 
                  write(mp2dm) nocc+nvir
                  ! write matrix element 
                  write(mp2dm) dm_mo_fort

                  ierr=talsh_tensor_destruct(tensor_slice)
      
                  ierr = talsh_shutdown()
            
!                  if (debug) then
!                  print *,'<<<<< DM, complete, before formatting to DM_complex'
!                        do j=1,(nocc+nvir)
!                           do i=1,(nocc+nvir)
!                              write(mp2dm,'(2x,2F28.16,2x,2I4)') dm_mo_fort(i,j),i,j
!                           end do
!                        end do
!                        print *,'>>>>> DM, complete, after formatting to DM_complex'
!                  end if
            
                  close(mp2dm,status='KEEP')
            
            end subroutine write_matrix

!-------------------------------------------------------------------------------------------------------------------

        subroutine exacorr_mp2no_spaces (exa_input,ao_space,occ_space,vir_space,ao_space_id,occ_space_id,vir_space_id)

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
         integer                    :: my_MPI_rank, MPI_size

         ! variables for the definition of colors inside the spaces
         type(color_symmetry_t), allocatable    :: color(:)
         integer(INTL)                          :: labs

         ! Loop counters and error codes
         integer(INTD) :: ierr, lsh
         integer :: l, nl
         integer :: num_blocks

         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         call interface_mpi_comm_size(global_communicator,MPI_size)

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
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write(*,*) 'ExaCorr: Inferred TAVP-MNG depth = ',tavp_mng_depth

         !Register the ao vector space:
         branch_factor = max(2, (int(get_num_segments(int(nao,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('AO_space_mp2',basis_ao,ao_space_id,ao_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') & 
           ' ExaCorr: Registered AO vector space of dimension   ',nao, &
           ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1)) !debug
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) call ao_space%print_it()

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
         ierr=exatns_space_register('Occupied_mp2',basis_occ,occ_space_id,occ_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(occ) failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') &
           ' ExaCorr: Registered OCC vector space of dimension  ',nocc, &
           ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1))
         num_blocks=(branch_factor*(2**(tavp_mng_depth-1)))**2

         !Register vir mo vector space:
         branch_factor = max(2, (int(get_num_segments(int(nvir,8),int(exa_input%exa_blocksize,8)),4) - 1)/tavp_mng_depth + 1)
         ierr=exatns_space_register('Virtual_mp2',basis_vir,vir_space_id,vir_space,branch_factor)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register(vir) failed!')
         if(exa_input%print_level.gt.8 .and. my_MPI_rank == 0) write (*,'(A,I3,A,I3)') &
           ' ExaCorr: Registered VIRT vector space of dimension ',nvir, &
           ': Number of segments = ', branch_factor*(2**(tavp_mng_depth-1)) !debug
         num_blocks=num_blocks*(branch_factor*(2**(tavp_mng_depth-1)))**2
         
         if(exa_input%print_level.gt.-1 .and. my_MPI_rank == 0) write (*,'(A,I9,A,I9)') &
           'number of blocks = ',num_blocks,'; number of MPI processes =',MPI_size
         
         if(num_blocks<MPI_size) then
            write (*,'(A)') 'In the current implementation there is a scaling limit.'
            write (*,'(A)') 'Using your configuration there would be nodes doing no work.'
            write (*,'(A)') 'Therefore decrease the number of MPI processes or .EXA_BLOCKSIZE'
            call quit('Not enough work for MPI processes')
         end if

        end subroutine exacorr_mp2no_spaces


        subroutine exacorr_mp2no_methods (exa_input,ao_space_id,occ_space_id,vir_space_id, &
                                       f_delta,orb_set,denom,denom3,j_chol, ff_set)

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
         type(delta_t)                 :: f_delta  ! fill delta matrix
         type(set_energy_t)            :: orb_set  ! place orbital energies on diagonal
         type(denom_t)                 :: denom    ! method to scale CC amplitudes by the inverse of orbital energy differences
         type(denom3_t)                :: denom3   ! scale with fixed occ value
         type(chol_J)                  :: j_chol   ! Cholesky: compute diagonal elements
         type(set_ff_t)                :: ff_set   ! method to add property integrals

         integer                :: nocc,nvir           ! the size of the mo basis for occupied and virtual orbitals (Kramers pairs)
         integer, allocatable   :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
         integer, allocatable   :: mo_oo(:),mo_vo(:),mo_ov(:),mo_vv(:) ! joined lists needed for dm init (ugly code, should be refactored)

         integer(INTD)                :: ierr, spin

         ! types of the initializer methods
         type(tens_printer_t)         :: tens_printer         ! generic method, but needs to be defined here in specific form
         type(compute_2e_ao_tensor_t) :: aoint_calculator     ! method to initialize the ao-integrals tensor
         type(compute_2e_ao_seg_t)    :: aoint_seg            ! only for segment
         type(init_mocoef_tensor_t)   :: occupied_init(2), virtual_init(2) ! methods to initialize MO coeffcient tensors
         type(init_dm_tensor_t)       :: dm_init(4)           ! array of initializers for MO transformation algorithm 2
         type(set_diagonal_t)         :: one_diag             ! identity matrix
         type(set_1el_t)              :: one_el(3)            ! one electron integrals
         type(square_t)               :: t_square             ! square matrix elements
         type(set_projection_t)       :: set_proj             ! projection matrix for triples
         type(set_zero_t)             :: cmplx_zero,imag_zero ! zero out tensors
         type(chol_diag)              :: diag_chol            ! Cholesky: compute diagonal elements

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
         ierr=exatns_method_register('Aoint_calculator_mp2',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(aoint) failed!')

         !Register the ao integral generation for a segment
         !ierr=exatns_method_register('AointSegment_mp2',aoint_seg)
         !if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(AointSegment) failed!')

         !Register the denominator initializer as a method
         call denom%denom_t_init(mo_occ,nocc,mo_vir,nvir,exa_input%level_shift,occ_space_id)
         ierr=exatns_method_register('Denominate_mp2',denom)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(Denominate) failed!')


         !Register the standard print method with treshold such that it only prints interesting elements
         call tens_printer%reset_thresh(1.0D-6)
         ierr=exatns_method_register('PrintTensor_mp2',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(printer) failed!')

         !Register method to set diagonal to one
         call one_diag%set_diagonal_init(ONE,.TRUE.)
         ierr=exatns_method_register('Unity_mp2',one_diag)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(Unity) failed!')


         one_el_exist = exist_one_el(e_core)
!         if (one_el_exist) then
            !Register methods to read one electron integrals
!            call one_el(1)%set_1el_init(mo_occ,nocc,mo_occ,nocc,min_occ)
!            ierr=exatns_method_register('OE_oo_mp2',one_el(1))
!            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_oo_mp2) failed!')
!            call one_el(2)%set_1el_init(mo_occ,nocc,mo_vir,nvir,min_occ)
!            ierr=exatns_method_register('OE_ov_mp2',one_el(2))
!            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_ov_mp2) failed!')
!            call one_el(3)%set_1el_init(mo_vir,nvir,mo_vir,nvir,min_occ)
!            ierr=exatns_method_register('OE_vv_mp2',one_el(3))
!            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(OE_vv_mp2) failed!')

            !Register methods to make matrix with orbital energies on the diagonal
!            call orb_set%set_energy_init(mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,exa_input%level_shift)
!            ierr=exatns_method_register('OrbSet_mp2',orb_set)
!            if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register(OrbSet_mp2) failed!')
!         end if

         !Register method to add property integrals
!         call ff_set%set_ff_init(mo_occ,nocc,occ_space_id,mo_vir,nvir,vir_space_id,min_occ)
!         ierr=exatns_method_register('FF_set',ff_set)
!         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(FF_set_mp2) failed!')

         
         !Register method to set to ZERO
         call cmplx_zero%set_zero_init(3)
         ierr=exatns_method_register('ZERO',cmplx_zero)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(cmplx_zero) failed!')
         
         !Register method to set imaginary part to zero
!         call imag_zero%set_zero_init(1)
!         ierr=exatns_method_register('ZERO_I',imag_zero)
!         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register(imag_zero) failed!')


        end subroutine exacorr_mp2no_methods




        subroutine get_MP2_integrals (ao_space_id,ao_space,int_t,moint_scheme,print_level)

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
         type(exatns_intg_tens), intent(inout) :: int_t
         integer, intent(in)                   :: print_level 

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

!        initialise tensors to zero
         ierr=exatns_tensor_init(int_t%oooo,ZERO)
         ierr=exatns_tensor_init(int_t%ooov,ZERO)
         ierr=exatns_tensor_init(int_t%oovv,ZERO)
         ierr=exatns_tensor_init(int_t%vovo,ZERO)


!        Get tensor with AO integrals
         if (moint_scheme < 4 ) then
            ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            call print_date('allocated AO integral tensor')
         end if

         !Initialize ao-tensor:
         if (moint_scheme == 0) then
!        test code : initialize integrals to zero to be able to time the integral generation
           ierr=exatns_tensor_init(aoint_tensor)
           call print_date('initalized all ao integrals to zero')
         elseif (moint_scheme == 1 .or. moint_scheme == 2) then
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator_mp2')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           if(print_level > 6) write(*,*) " N (aoint) = ", print_tensornorm2(aoint_tensor)
           call print_date('finished ao integral calculation')


!          Get oovv tensor
           call get_integral_tensor (aoint_tensor,int_t%oovv,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (oovv) =",print_tensornorm2(int_t%oovv)
           call print_date('finished oovv integrals')

!          Get ooov tensor
           call get_integral_tensor (aoint_tensor,int_t%ooov,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (ooov) =",print_tensornorm2(int_t%ooov)
           call print_date('finished ooov integrals')

!          Get oooo tensor
           call get_integral_tensor (aoint_tensor,int_t%oooo,12,print_level,algorithm=moint_scheme)
           if(print_level > 6) write(*,*) " N (oooo) =",print_tensornorm2(int_t%oooo)
           call print_date('finished oooo integrals')

!          Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!          Get vovo tensor
           call get_integral_tensor (aoint_tensor,int_t%vovo,0,print_level,algorithm=moint_scheme)
           ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tens. too large ')
           ierr=exatns_tensor_init(ovvo_tensor,ZERO)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: initialization failed ')
           call get_integral_tensor (aoint_tensor,ovvo_tensor,0,print_level,algorithm=moint_scheme)
!          Make anti-symmetrized vovo tensor
           ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
           ierr=exatns_tensor_destroy(ovvo_tensor)
           if(print_level > 6) write(*,*) " N (vovo) =",print_tensornorm2(int_t%vovo)
           call print_date('finished vovo integrals')

         elseif (moint_scheme == 3) then
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator_mp2')
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
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ') 
           ierr=exatns_tensor_init(a2int_tensor,ZERO)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: init. failed ')
           call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
           call print_date('vv halftransformation completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%vovo,0,print_level)
           call print_date('vovo integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

           ! get the integrals with ov as first transform
           a2int_id(1:2)   = vovo_id(2:3)
           a2int_root(1:2) = vovo_root(2:3)
           ierr=exatns_tensor_create(a2int_tensor,"a2ov",a2int_id,a2int_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
           ierr=exatns_tensor_init(a2int_tensor,ZERO)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: init. failed ')
           call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
           call print_date('ov halftransformation completed')
           ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           ierr=exatns_tensor_init(ovvo_tensor,ZERO)
           call get_ht2_integral_tensor (a2int_tensor,ovvo_tensor,0,print_level)
           call print_date('ovvo integrals completed')
!          Make anti-symmetrized vovo tensor
           ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
           ierr=exatns_tensor_destroy(ovvo_tensor)
           if(print_level > 6) write(*,*) " N (vovo) =",print_tensornorm2(int_t%vovo)
           call print_date('vovo integrals anti-symmetrized')
           call get_ht2_integral_tensor (a2int_tensor,int_t%oovv,12,print_level)
           if(print_level > 6) write(*,*) " N (oovv) =",print_tensornorm2(int_t%oovv)
           call print_date('oovv integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

           ! get the integrals with oo as first transform
           a2int_id(1:2)   = vovo_id(2)
           a2int_root(1:2) = vovo_root(2)
           ierr=exatns_tensor_create(a2int_tensor,"a2oo",a2int_id,a2int_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
           ierr=exatns_tensor_init(a2int_tensor,ZERO)
           if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: init. failed ')
           call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
           call print_date('oo halftransformation completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%ooov,12,print_level)
           if(print_level > 6) write(*,*) " N (ooov) =",print_tensornorm2(int_t%ooov)
           call print_date('ooov integrals completed')
           call get_ht2_integral_tensor (a2int_tensor,int_t%oooo,12,print_level)
           if(print_level > 6) write(*,*) " N (oooo) =",print_tensornorm2(int_t%oooo)
           call print_date('oooo integrals completed')
           ierr=exatns_tensor_destroy(a2int_tensor)

         elseif (moint_scheme == 4) then
           
           call print_date('choosing moint_scheme 4')
           ! in this algorithm we only keep part of the AO and HT integrals in memory
           subspace_id(1:3) = ao_space_root
           call ao_space%get_level_composition(int(1,INTD),ao_subspaces,num_subspaces,ierr)
           do sub4 = 1, num_subspaces
              subspace_id(4)=ao_subspaces(sub4)
              ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),subspace_id,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
              ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator_mp2')
              if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')

              ! Setting dimensions for the half transformed integral tensor, the last two are always AO spaces
              a2int_id(3:4)   = ao_space_id
              a2int_root(3:4) = subspace_id(3:4)

              ! get the integrals with vv as first transform
              a2int_id(1:2)   = vovo_id(1)
              a2int_root(1:2) = vovo_root(1)
              ierr=exatns_tensor_create(a2int_tensor,"a2vv",a2int_id,a2int_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
              ierr=exatns_tensor_init(a2int_tensor,ZERO)
              call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%vovo,0,print_level)
              ierr=exatns_tensor_destroy(a2int_tensor)

              ! get the integrals with ov as first transform
              a2int_id(1:2)   = vovo_id(2:3)
              a2int_root(1:2) = vovo_root(2:3)
              ierr=exatns_tensor_create(a2int_tensor,"a2ov",a2int_id,a2int_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
              ierr=exatns_tensor_init(a2int_tensor,ZERO)
              call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
              ierr=exatns_tensor_create(ovvo_tensor,"ovvo",ovvo_id,ovvo_root,EXA_DATA_KIND_C8) ! auxilliary for get_CC_integrals
              ierr=exatns_tensor_init(ovvo_tensor,ZERO)
              call get_ht2_integral_tensor (a2int_tensor,ovvo_tensor,0,print_level)
!             Make anti-symmetrized vovo tensor
              ierr=exatns_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*K()",int_t%vovo,ovvo_tensor,one_tensor,MINUS_ONE)
              ierr=exatns_tensor_destroy(ovvo_tensor)
              call get_ht2_integral_tensor (a2int_tensor,int_t%oovv,12,print_level)
              ierr=exatns_tensor_destroy(a2int_tensor)

              ! get the integrals with oo as first transform
              a2int_id(1:2)   = vovo_id(2)
              a2int_root(1:2) = vovo_root(2)
              ierr=exatns_tensor_create(a2int_tensor,"a2oo",a2int_id,a2int_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('ERROR in get_CC_integrals: tensor too large ')
              ierr=exatns_tensor_init(a2int_tensor,ZERO)
              call ao2mo_exat (aoint_tensor,a2int_tensor,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%ooov,12,print_level)
              call get_ht2_integral_tensor (a2int_tensor,int_t%oooo,12,print_level)
              ierr=exatns_tensor_destroy(a2int_tensor)
              ierr=exatns_tensor_destroy(aoint_tensor)
              write (text,'(I6,A,I6)') sub4,' of',num_subspaces
              call print_date('finished transforming subspace'//text)
           end do

           if(print_level > 6) write(*,*) " N (vovo) =",print_tensornorm2(int_t%vovo)
           if(print_level > 6) write(*,*) " N (oovv) =",print_tensornorm2(int_t%oovv)
           if(print_level > 6) write(*,*) " N (ooov) =",print_tensornorm2(int_t%ooov)
           if(print_level > 6) write(*,*) " N (oooo) =",print_tensornorm2(int_t%oooo)

         end if

!        Clean-up
         if (moint_scheme < 4) then
            ierr=exatns_tensor_destroy(aoint_tensor)
         end if

        end subroutine get_MP2_integrals


        subroutine init_fock_exat (comm_t,int_t,oo_id,oo_root,ov_id,ov_root,vv_id,vv_root,&
                              f_delta,denom,orb_set,ff_set,exa_input,scf_energy)

          use exacorr_global
          use exacorr_tensor_methods
          use exacorr_ao_to_mo, only : print_tensornorm2

          implicit none
  
          type(exatns_comm_tens), intent(inout)  :: comm_t
          type(exatns_intg_tens), intent(inout)  :: int_t
          integer(INTL),dimension(2), intent(in) :: oo_root, ov_root, vv_root
          integer(INTD),dimension(2), intent(in) :: oo_id, ov_id, vv_id 
          type(delta_t)                          :: f_delta
          type(denom_t)                          :: denom
          type(set_energy_t)                     :: orb_set
          type(set_ff_t)                         :: ff_set
          type(exacc_input), intent(in)          :: exa_input
          real(8), intent(out)                   :: scf_energy
          integer                                :: iff=1   !TODO: make a loop for field strngths

          real(8)                    :: ev_energy, off_diag
          !        scalars (need to be defined as tensor types)
          type(tens_rcrsv_t)         :: val_tensor, tensor_one, tensor_h
          !        identity matrix for contraction 
          type(tens_rcrsv_t)         :: id_mat_occ, orb_tensor
          integer(INTD)              :: ierr
          complex(8)                 :: result_val
          !        finite field stuff
          integer                    :: i,i_prop

          if (one_el_exist) then

            ! Allocate two auxilliary tensors needed to extract results or multiply by one
            ierr=exatns_tensor_create(val_tensor,"val_tensor",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_create(tensor_one,"tensor_one",EXA_DATA_KIND_C8)
            ierr=exatns_tensor_init(tensor_one, ONE)

            ! Make the density matrix (which is an identity matrix in the MO basis)
            ierr=exatns_tensor_create(id_mat_occ,"id_mat_occ",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_transform(id_mat_occ,'Unity_mp2')
            if (exa_input%print_level >0) call print_date('Density matrix initialized')

            ! Initialize the sub blocks of the Fock matrix with the (effective) one-electron integrals
            ierr=exatns_tensor_init(comm_t%foo,'OE_oo_mp2')
            ierr=exatns_tensor_init(comm_t%fov,'OE_ov_mp2')
            ierr=exatns_tensor_init(comm_t%fvv,'OE_vv_mp2')
            if(ierr.ne.EXA_SUCCESS) call quit(ierr,'Creating one electron integrals failed!')
            if (exa_input%print_level >0) call print_date('one electron integrals added to Fock Matrix')

            !adding property integrals
            do i=1,exa_input%nff(1)
              call get_ff_ind (i_prop,exa_input%ff_names(i))
              call ff_set%reset(i_prop)

              ierr=exatns_tensor_create(tensor_h,"tensor_h",oo_id,oo_root,EXA_DATA_KIND_C8)
              write(*,*) 'ff_set create error', ierr
              ierr=exatns_tensor_init(tensor_h, ff_set)
              write(*,*) 'ff_set init error', ierr
              ierr=exatns_tensor_contract("F(p,q)+=P(p,q)*O()",comm_t%foo,tensor_h,tensor_one,exa_input%ff(i,iff))
              write(*,*) 'ff_set contract error', ierr
              write(*,*) " prop =",print_tensornorm2(tensor_h)
              write(*,*) " field =",exa_input%ff(i,iff)
              ierr=exatns_tensor_traverse(tensor_h,'PrintTensor_mp2')
              ierr=exatns_tensor_destroy(tensor_h)
              write(*,*) " F_oo+prop =",print_tensornorm2(comm_t%foo)
              write(*,*) 'ff_set destroy error', ierr
              ierr=exatns_tensor_traverse(comm_t%foo,'PrintTensor_mp2')
            
            end do


            ! Calculate the Hartree-Fock energy
            ierr=exatns_tensor_init(val_tensor, ZERO)
            ierr=exatns_tensor_contract("R()=F(i,j)*I(i,j)",val_tensor,comm_t%foo,id_mat_occ)
            ierr=exatns_tensor_create(tensor_h,"tensor_h",oo_id,oo_root,EXA_DATA_KIND_C8)
            ierr=exatns_tensor_get_scalar(val_tensor,result_val)
            ierr=exatns_tensor_init(tensor_h, ZERO)
            ierr=exatns_tensor_contract("H(i,j)+=V(i,k,j,l)*I(k,l)",tensor_h,int_t%oooo,id_mat_occ, ONE_HALF)
            ierr=exatns_tensor_contract("R()+=H(i,j)*I(i,j)",val_tensor,tensor_h,id_mat_occ)
            ierr=exatns_tensor_get_scalar(val_tensor,result_val)
            scf_energy=real(result_val,8) ! Casting into real(8)
            ierr=exatns_tensor_destroy(tensor_h)
            scf_energy=scf_energy+e_core
            write(*,*) ""
            write (*,'(A,F23.15)') "  SCF energy (single determinant, calculated) = ", scf_energy
            write(*,*) ""
            ! Make the complete Fock matrix by addding the valence electron contributions
            ierr=exatns_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",comm_t%foo,int_t%oooo,id_mat_occ, ONE)
            ierr=exatns_tensor_contract("F(i,a)+=V(k,i,l,a)*I(k,l)",comm_t%fov,int_t%ooov,id_mat_occ, ONE)
            ierr=exatns_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",comm_t%fvv,int_t%vovo,id_mat_occ, ONE)
            if (exa_input%print_level >0) call print_date('two electron integrals added to Fock Matrix')

            if (exa_input%print_level >10) then
              ierr=exatns_tensor_traverse(comm_t%foo,'PrintTensor_mp2')
              ierr=exatns_tensor_traverse(comm_t%fov,'PrintTensor_mp2')
              ierr=exatns_tensor_traverse(comm_t%fvv,'PrintTensor_mp2')
            end if

            call update_eps(comm_t, f_delta, denom, orb_set, oo_root, vv_root, oo_id, vv_id, exa_input%print_level)

            ! Subtract orbital energies from the diagonal (they are moved to the other side of the CC equations)
            ierr=exatns_tensor_create(orb_tensor,"orb_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            ierr=exatns_tensor_init(orb_tensor,orb_set)
            if(ierr.ne.EXA_SUCCESS) call quit('init_fock: set orbitals failed!')
            ierr=exatns_tensor_contract("F(i,j)+=E(i,j)*O()",comm_t%foo,orb_tensor,tensor_one, MINUS_ONE)
            ierr=exatns_tensor_destroy(orb_tensor)
            ierr=exatns_tensor_create(orb_tensor,"orb_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
            ierr=exatns_tensor_init(orb_tensor,orb_set)
            if(ierr.ne.EXA_SUCCESS) call quit('init_fock: set orbitals failed!')
            ierr=exatns_tensor_contract("F(a,b)+=E(a,b)*O()",comm_t%fvv,orb_tensor,tensor_one, MINUS_ONE)
            ierr=exatns_tensor_destroy(orb_tensor)
            if (exa_input%print_level >0) call print_date('orbital energies subtracted from the diagonal')
            
            if (exa_input%print_level >0) then
              ! Compute norms (check whether this single determinant is contructed from optimized orbitals)
              write(*,*) ""
              ierr=exatns_tensor_init(val_tensor, ZERO)
              ierr=exatns_tensor_contract("R()=F(i,j)*F(i,j)",val_tensor,comm_t%foo,comm_t%foo)
              ierr=exatns_tensor_get_scalar(val_tensor,result_val)
              off_diag=real(result_val,8) ! Casting into real(8)
              write(*,*) "    Norm of Fock matrix (OO, no diagonal) = ", off_diag
            
              ierr=exatns_tensor_init(val_tensor, ZERO)
              ierr=exatns_tensor_contract("R()=F(i,a)*F(i,a)",val_tensor,comm_t%fov,comm_t%fov)
              ierr=exatns_tensor_get_scalar(val_tensor,result_val)
              off_diag=real(result_val,8) ! Casting into real(8)
              write(*,*) "    Norm of Fock matrix (OV, all        ) = ", off_diag

              ierr=exatns_tensor_init(val_tensor, ZERO)
              ierr=exatns_tensor_contract("R()=F(a,b)*F(a,b)",val_tensor,comm_t%fvv,comm_t%fvv)
              ierr=exatns_tensor_get_scalar(val_tensor,result_val)
              off_diag=real(result_val,8) ! Casting into real(8)
              write(*,*) "    Norm of Fock matrix (VV, no diagonal) = ", off_diag
              write(*,*) ""
            end if 
            ierr=exatns_tensor_destroy(tensor_h)
            ierr=exatns_tensor_destroy(id_mat_occ)
            ierr=exatns_tensor_destroy(val_tensor)
            ierr=exatns_tensor_destroy(tensor_one)
          else
            ierr=exatns_tensor_init(comm_t%foo,ZERO)
            ierr=exatns_tensor_init(comm_t%fov,ZERO)
            ierr=exatns_tensor_init(comm_t%fvv,ZERO)

            call get_scf ( scf_energy )
            write(*,*) ""
            write (*,'(A,F23.15)') "    SCF energy ( read ) = ", scf_energy
            write(*,*) ""

          end if

        end subroutine init_fock_exat

        subroutine update_eps(comm_t, f_delta, denom, orb_set, oo_root, vv_root, oo_id, vv_id, print_level)
          
          use exacorr_tensor_methods
          use exacorr_utils, only : print_orbital_energy

          implicit none

          type(exatns_comm_tens)                 :: comm_t
          type(delta_t)                          :: f_delta 
          type(denom_t)                          :: denom
          type(set_energy_t)                     :: orb_set 
          integer(INTL),dimension(2), intent(in) :: oo_root, vv_root
          integer(INTD),dimension(2), intent(in) :: oo_id, vv_id 
          integer, intent(in)                    :: print_level

          real(8),allocatable        :: eps_occ(:),eps_vir(:)
          integer                    :: nocc, nvir, i
          integer(INTL),dimension(2) :: delta_set
          type(tens_rcrsv_t)         :: d_mat, res_val
          complex(8)                 :: res
          integer(INTD)              :: ierr

          nocc=oo_root(1)
          nvir=vv_root(1)

          allocate(eps_occ(nocc))
          allocate(eps_vir(nvir))

          ierr=exatns_tensor_create(res_val,"res_val",EXA_DATA_KIND_C8)

          ierr=exatns_tensor_create(d_mat,"d_mat",oo_id,oo_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(d_mat,ZERO)
          do i = 1, nocc
              delta_set=i
              call f_delta%reset(delta_set,2)
              ierr=exatns_tensor_transform(d_mat,f_delta)

              ierr=exatns_tensor_init(res_val,ZERO)
              ierr=exatns_tensor_contract("R()+=W(p,q)*C(p,q)",res_val,comm_t%foo,d_mat)
              ierr=exatns_tensor_get_scalar(res_val,res)
              eps_occ(i)=real(res,8)

          end do
          ierr=exatns_tensor_destroy(d_mat)

          ierr=exatns_tensor_create(d_mat,"d_mat",vv_id,vv_root,EXA_DATA_KIND_C8)
          ierr=exatns_tensor_init(d_mat,ZERO)
          do i = 1, nvir
              delta_set=i
              call f_delta%reset(delta_set,2)
              ierr=exatns_tensor_transform(d_mat,f_delta)

              ierr=exatns_tensor_init(res_val,ZERO)
              ierr=exatns_tensor_contract("R()+=W(p,q)*C(p,q)",res_val,comm_t%fvv,d_mat)
              ierr=exatns_tensor_get_scalar(res_val,res)
              eps_vir(i)=real(res,8)

          end do
          ierr=exatns_tensor_destroy(d_mat)
          print *, "eps_occ", eps_occ
          print *, "eps_occ", eps_vir
          
          if (print_level.gt.3) call print_orbital_energy(eps_occ,nocc,eps_vir,nvir,2)

          call denom%reset(eps_occ,eps_vir,nocc,nvir)
          call orb_set%reset(eps_occ,eps_vir,nocc,nvir)

          deallocate(eps_occ)
          deallocate(eps_vir)
          ierr=exatns_tensor_destroy(res_val)

        end subroutine update_eps

!------------------------------------------------------------------------------------------------------------------
        subroutine get_transform_Fock (exa_input,Matrix_file)
            ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
            ! The new basis set comes from the digonalization of density matrix.
            ! The all natural orbitals and occupation numbers are saved in file NOs_MO 
            ! The transformed Fock matrix are saved in file FM_in_NatOrb
            ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix. 
            ! This routine doesn't truncate the new basis set.    
            ! Written by Xiang Yuan   
   
               use matrix_complex_to_quaternion 
               use write_orbital_to_file
               use interface_to_mpi
   
                     character(len=*),intent(in)        :: Matrix_file
                     type(exacc_input), intent(in)      :: exa_input
   
                     real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                     complex(kind=8), pointer           :: MP2_density_comp(:,:)   ! MP2 density matrix from DM_complex
                     real(8), pointer                   :: MP2_density_quat(:,:,:) ! quaternion MP2 density matrix 
                     real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                     real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
                     real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
                     integer                            :: sz,flb                  ! size of matrix and flag for block
                     character                          :: mat                     ! flag used for matrix basis transform
                     integer                            :: i
                     integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode
   
                                                           
                  call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
                  if (my_MPI_rank == 0) then

                     call read_matrix_comp(MP2_density_comp,Matrix_file)
                     call print_date('Read Density Matrix from DM_complex')
   
                     call complex_to_quaternion(exa_input,MP2_density_comp,MP2_density_quat)
                     call print_date('Transform DM from complex into quaternion')
   
                     call diagoanl_matrix_quater(MP2_density_quat,N_OCCUPATION,NOs_MO)    
                     call print_date('Qdiag Density matrix')   
   
                     call sort_values_vectors (N_OCCUPATION,NOs_MO)
                     call print_date('Sort quaterbion natural orbital')
   
                     call quater_write(NOs_MO,N_OCCUPATION,'NOs_MO')
                     call print_date('Write the all natural orbital and occupation number in file NOs_MO')
                       
                     flb = 0                     ! we read vv block of Fock matrix
                     sz = exa_input%nkr_vir      ! the size of the FM is the same as DM
                     call quater_read_from_dirac_fock(Fock_HF,sz,flb)
                     call print_date('Read Fock matrix from DFCOEF')
   
                     mat = 'Matrix'
                     call basis_transform(Fock_HF,NOs_MO,Fock_NO,mat)
                     call print_date('Get full vv Fock matrix in natural orbital basis')
   
                     call quater_write(Fock_NO,"FM_in_NatOrb")
                     call print_date('Write the Full Fock matrix of natural orbital in file FM_in_NatOrb')
                  end if 
                           
         end subroutine get_transform_Fock
   
            
   !------------------------------------------------------------------------------------------------------------------------------
         subroutine get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION_reduced,New_orbital_energy,file_orb,N_reduced)
         ! This routine is used for getting new truncated canonical orbital in AO basis
         ! and the corresponding orbital energy.
   
               use matrix_complex_to_quaternion 
               use write_orbital_to_file
               use interface_to_mpi
   
   
               type(exacc_input), intent(in)      :: exa_input
               real(kind=8),intent(in)            :: threshold
               real(kind=8), pointer,intent(out)  :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis 
               real(kind=8), pointer,intent(out)  :: N_OCCUPATION_reduced(:) ! occupation number from qdia of DM
               real(kind=8), pointer,intent(out)  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
               character(len=*),intent(in)        :: file_orb
               integer,intent(inout)              :: N_reduced               ! number of reduced natural orbital 
   
   
               real(kind=8), pointer              :: NOs_MO_reduced(:,:,:)   ! natural orbitals which have occupation number larger than threshold
               real(kind=8), pointer              :: New_cano_MO(:,:,:)      ! new canonical orbital in AO basis
               real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
               real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
               real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital 
               integer                            :: sz,flb, fl              ! size of matrix and flag for block
               character                          :: mat                     ! flag used for matrix basis transform
               integer                            :: i,j,k
               integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode
               
            
            call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
            if (my_MPI_rank == 0) then   
               call quater_read(NOs_MO_reduced,N_OCCUPATION_reduced,threshold,N_reduced,file_orb)
               call print_date('Get the truncated the natural orbital and occupation number from NOs_MO')
               write(*,"(A45,2X,D8.3)") 'Current MP2 occupation number Threshold: ', threshold
               write(*,"(A15,2X,I4,2X,A15)") 'Now we have', N_reduced, 'New Orbitals'
   
   
               flb = 0                       ! we only read virtual orbital 
               sz = exa_input%nkr_vir        ! the number of virtual orbital 
               call quater_read_from_dirac_orbital(HF_AO,sz,flb)
               call print_date('Read HF orbital coefficient from DFCOEF')
   
               
               flb = 0                         ! we read vv block of Fock matrix
               sz = exa_input%nkr_vir          ! the size of the FM is the same as DM
               call quater_read_from_dirac_fock(Fock_HF,sz,flb)
               call print_date('Read Fock matrix from DFCOEF')
   
   
               mat = 'Matrix'
               call basis_transform(Fock_HF,NOs_MO_reduced,Fock_NO,mat)
               call print_date('Get reduced vv Fock Matrix in Natural Orbital Basis')
   
               
               call get_Cano_Orbital(NOs_MO_reduced,New_cano_MO,Fock_NO,New_orbital_energy)
               call print_date('Get new Canonical orbital in MO basis')
   
   
   !      Debug part, write HF orbital into the new file to test the code.
   !*********************************************** 
   !            New_cano_MO = 0.0                
   !            do i=1, size(New_cano_MO,2)      
   !                  New_cano_MO(i,i,1) = 1.0   
   !            end do                            
   !***********************************************
               
               fl = 0 
               call basis_transform(New_cano_MO,HF_AO,NOs_AO,fl)
               call print_date('Get new orbital in AO basis')
               
   
               print *, " Print several lowest New orbital Energy"
               if (N_reduced < 10) then
                     do i=1, N_reduced
                           !write(*,"(A8,2X,I2,2X,A8,2X,F20.11)") ' Orbital',i,':Energy:', New_orbital_energy(i)
                           print *, ' Orbital',i,'  :Energy: ', New_orbital_energy(i)
                     end do
               else
                     do i=1, 10
                           !write(*,"(A8,2X,I2,2X,A8,2X,F20.11)") ' Orbital',i,':Energy:', New_orbital_energy(i)
                           print *, ' Orbital',i,'  :Energy: ', New_orbital_energy(i)
                     end do
               end if

            end if
         end subroutine
   
   
   !-----------------------------------------------------------------------------------------------------------------------------
         subroutine write_to_newfile (exa_input,NOs_AO,N_OCCUPATION,New_orbital_energy,file_old,file_new)
     !      write new orbital information into new file 'MP2NOs_AO', which is similar to DFCOEF 
     !      copy the basis information from DFCOEF and use new orbital energy and coefficient to replace the old orbital.   
     !      revise subroutine write_mos_to_dirac of exacorr_dirac.F
   
            use exacorr_mo
            use exacorr_datatypes
            use exacorr_utils
            use write_orbital_to_file
            use interface_to_mpi
     
            implicit none
       
            type(exacc_input), intent(in)         :: exa_input
            character(len=*),intent(in)           :: file_old,file_new
            real(8),pointer,intent(in)            :: NOs_AO(:,:,:)          ! new canonical orbital in AO basis
            real(kind=8),pointer,intent(in)       :: N_OCCUPATION(:)        ! occupation number
            real(kind=8),pointer,intent(in)       :: New_orbital_energy(:)  ! new orbital energy 
   
   
            type(basis_set_info_t)   :: ao_basis
            type(qumo)               :: qorbital        ! save the new canonical orbital 
            type(qumo)               :: qorbital_hf     ! save the original HF canonical orbital 
   
            integer                  :: nocc            ! occupied orbital 
            integer                  :: nao, nmo, nno   ! number of atomic orbital; molecular orbital; natural orbital
            integer,allocatable      :: index(:)
            integer                  :: i,j
            integer                  :: ierr
            real(kind=8)             :: diff_1, diff_2, diff_3, diff_4
            integer       :: my_MPI_rank                    ! needed to start exatensor, MPI should already be started in threaded mode
   

      call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
      if (my_MPI_rank == 0) then 
   !        read information of from DFCOEF
            call quater_read_from_dirac_all(ao_basis,qorbital_hf,file_old,ierr)
   
            nocc = exa_input%nkr_occ
            nao = qorbital_hf%nao
            nno = size(NOs_AO,2)
            nmo = nocc+nno
            allocate(index(nmo))
            do i = 1, nmo
               index(i) = i ! we set the index compatible with the fact that we store all mo's
            end do
   
            call alloc_mo(qorbital,nao,nmo,4)
            qorbital%allocated = qorbital_hf%allocated
            qorbital%restricted = qorbital_hf%restricted
            qorbital%complex_algebra = qorbital_hf%complex_algebra
            qorbital%total_energy = qorbital_hf%total_energy
            qorbital%fermion_irrep = qorbital_hf%fermion_irrep
            qorbital%boson_irrep = 0
            qorbital%index = index
   
   
            do i = 1, nocc
               qorbital%energy(i)=qorbital_hf%energy(i)
            end do
   
            do i = nocc+1, qorbital%nmo
               qorbital%energy(i)=New_orbital_energy(i-nocc)
   !            qorbital%energy(i)=qorbital_hf%energy(i)
            end do
   
            do j = 1, nocc
               do  i = 1, qorbital%nao
                     qorbital%coeff(i,j,1)=qorbital_hf%coeff(i,j,1)
                     qorbital%coeff(i,j,2)=qorbital_hf%coeff(i,j,2)
                     qorbital%coeff(i,j,3)=qorbital_hf%coeff(i,j,3)
                     qorbital%coeff(i,j,4)=qorbital_hf%coeff(i,j,4)
               enddo
            enddo
   
   
            do  j = nocc+1, qorbital%nmo
               do  i = 1, qorbital%nao
                     qorbital%coeff(i,j,1)=NOs_AO(i,j-nocc,1)
                     qorbital%coeff(i,j,2)=NOs_AO(i,j-nocc,2)
                     qorbital%coeff(i,j,3)=NOs_AO(i,j-nocc,3)
                     qorbital%coeff(i,j,4)=NOs_AO(i,j-nocc,4)
               enddo
            enddo
   
   
            call quater_write_mos_to_dirac(ao_basis,qorbital,file_new)
            call print_date('Write new orbital into MP2NOs_AO file like DFCOEF')
   
      end if 
      end subroutine write_to_newfile

!-----------------------------------------------------------------------------------------------------
#else
        implicit none

        public exacorr_cc_driver

       contains

        subroutine exacorr_mp2no_driver (exa_input,N_cano)
         use exacorr_datatypes
         implicit none
         type(exacc_input), intent(in) :: exa_input
         integer, intent(inout)        :: N_cano                  ! number of new canonical orbital
         write(*,*) 'called exatensor version of cc_driver in a serial version'
         write(*,*) 'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_mp2no_driver

#endif

end module exacorr_mp2no
