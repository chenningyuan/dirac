module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exatensor

        implicit none
        private
        public exacorr_compute_mo_integrals
        public test_integral_generation

       contains

#if defined (VAR_MPI)
        subroutine exacorr_compute_mo_integrals (nmo, mo_list)

!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor

!        Written by Lucas Visscher, march 2018

         use interface_to_mpi
         use dirac_MO
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use module_interest_interface
         use module_interest_eri
         use exatensor
         use, intrinsic:: ISO_C_BINDING


         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: i, ierr
         integer             :: nshells                    ! the number of shells (orbitals with the same center and l-value)
         integer             :: nao                        ! the number of atomic orbitals)
         type(basis_func_info_t), allocatable :: gto(:)    ! array with information about the shells

         type(subspace_basis_t):: basis_ao, basis_mo
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL):: nao_l,labs, ao_space_root
         integer(INTL):: nmo_l, mo_space_root
         class(h_space_t), pointer:: ao_space, mo_space
         integer:: mpi_th_provided, my_rank,my_role, ao_space_id, mo_space_id
         integer:: lsh, l, nl
         integer:: l_col=1

         type(tens_rcrsv_t):: aoint_tensor
         type(tens_rcrsv_t):: mocoef_tensor
!later   type(tens_rcrsv_t):: moint_tensor
         type(compute_2e_ao_tensor_t) :: aoint_calculator
         type(init_mocoef_tensor_t) :: mocoefa_init
         type(init_mocoef_tensor_t) :: mocoefb_init
         type(tens_init_test_t):: init1369

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_rank)
         print*, " Process", my_rank," arrived in compute_mo_integrals_exa routine"

!        Initialize libraries
         call interest_initialize
         call initialize_global_data

         nao     = get_nao()   ! This is the otal number of ao's
         nao_l   = int(nao,8)  ! we need this a long int for exatensor

!        Get the basis set information
         call get_gtos(1,nao,gto,nshells)
         print*,"my_rank,nao,nshells",my_rank,nao,nshells

         !Create a basis for the ao vector space:
         ierr = 0
         call basis_ao%subspace_basis_ctor(nao_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         !Assign the colors for the ao space (this is where exatensor may break up into subspaces)
         allocate(color(nshells))
         labs = 0
         do lsh = 1, nshells !set basis functions
            nl   =  ncartesian(gto(lsh)%orb_momentum)
!           call color(lsh)%color_symmetry_ctor(lsh) ! for the moment take one shell as a color
            call color(lsh)%color_symmetry_ctor(l_col) ! give the all the same color
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
         print*,"my_rank registered vector space, ierr=",my_rank,ierr
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Get the MO-coeffcient tensor for the first index that is to be transformed
         ierr = 0
         nmo_l   = 2*int(nmo(1),8)
         call basis_mo%subspace_basis_ctor(nmo_l,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         !Set the functions inside this basis
         do labs = 1, nmo_l
            call basis_mo%set_basis_func(labs,BASIS_ABSTRACT,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.set_basis_func() failed!')
         enddo

         !Finalize the mo space basis
         call basis_mo%finalize(ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.finalize() failed!')

         !Register the mo vector space:
         ierr=exatns_space_register('MO_space',basis_mo,mo_space_id,mo_space)
         print*,"my_rank registered vector space, ierr=",my_rank,ierr
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Register the template initialization routine (for testing purposes only)
         call init1369%tens_init_test_ctor(13.69d0)
         ierr=exatns_method_register('SetTo13.69',init1369)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the mo coefficient initializer aa a method
         call mocoefa_init%init_mocoef_ctor(mo_list,nmo(1),'a')
         ierr=exatns_method_register('MOCoeff_initializer',mocoefa_init)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

!        Run ExaTENSOR within MPI_COMM_WORLD (called golbal_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(global_communicator)
         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role)
          print*,"i am ",my_rank,my_role
          if(my_role.eq.EXA_DRIVER) then

           print*,"i am ",my_rank,my_role," and will create two tensors"

           !Create ao-integral-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           !Create mo-coefficient-tensors:
           mo_space_root=mo_space%get_root_id(ierr)
           ierr=exatns_tensor_create(mocoef_tensor,'mocoeftensor', &
            (/ao_space_id,mo_space_id/),(/ao_space_root,mo_space_root/),EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

!          Initialize mocoef-tensor:
!          ierr=exatns_tensor_init(mocoef_tensor,'SetTo13.69')
           ierr=exatns_tensor_init(mocoef_tensor,'MOCoeff_initializer')
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           !Destroy aoint-tensor:
           ierr=exatns_tensor_destroy(aoint_tensor)
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           !Destroy mocoef-tensor:
           ierr=exatns_tensor_destroy(mocoef_tensor)
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           !Stop ExaTENSOR runtime:
           ierr=exatns_stop()

          end if

         else
          write(*,*) ' Process ',my_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Process ",my_rank, " leaving compute_mo_integrals_exa routine"
         return
        end subroutine exacorr_compute_mo_integrals

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

        subroutine test_integral_generation (print_level)

!This subroutine tests the parallel AO integral generation

!        Written by Lucas Visscher, march 2018

         use dirac_MO
         use exacorr_datatypes
         use exacorr_global
         use module_interest_interface
         use module_interest_eri

         implicit none
         integer, intent(in) :: print_level  ! 0: only show that we running this test
                                             ! 1: print only overall summary
                                             ! 2: print information for each batch

         integer             :: nao           ! the size of the atomic orbital basis

         integer :: ierr
         integer :: i, j, k, l, ish, jsh, ksh, lsh, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         real(8), parameter :: fijkl = 1.0
         real(8) :: batch_value, ddot
         real(8), allocatable :: batch_values(:,:,:,:)
         integer            :: nshells ! number of AO shells
         type(basis_func_info_t), allocatable :: gto(:)   ! array with basis function information
         real(8) :: relative_error

!        Make some noise so that we know we are here
         print*, " Testing AO integral generation by calculating all integrals once"

!        Initialize libraries
         call interest_initialize(print_level.gt.1)
         if (print_level > 0) print*, " Initialized interest library"
         call initialize_global_data
         if (print_level > 0) print*, " Initialized global data"

!        Retrieve basis set information
         nao     = get_nao()

!        Get AOs and allocate array to store batch values
         call get_gtos(1,nao,gto,nshells)
         allocate(batch_values(nshells,nshells,nshells,nshells))

         do lsh = 1, nshells
            ll   =  gto(lsh)%orb_momentum
            el   =  gto(lsh)%exponent(1)
            xl   =  gto(lsh)%coord(1)
            yl   =  gto(lsh)%coord(2)
            zl   =  gto(lsh)%coord(3)
            cl   =  gto(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            do ksh = 1, nshells
               lk   =  gto(ksh)%orb_momentum
               ek   =  gto(ksh)%exponent(1)
               xk   =  gto(ksh)%coord(1)
               yk   =  gto(ksh)%coord(2)
               zk   =  gto(ksh)%coord(3)
               ck   =  gto(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               do jsh = 1, nshells
                  lj   =  gto(jsh)%orb_momentum
                  ej   =  gto(jsh)%exponent(1)
                  xj   =  gto(jsh)%coord(1)
                  yj   =  gto(jsh)%coord(2)
                  zj   =  gto(jsh)%coord(3)
                  cj   =  gto(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  do ish = 1, nshells
                     li   =  gto(ish)%orb_momentum
                     ei   =  gto(ish)%exponent(1)
                     xi   =  gto(ish)%coord(1)
                     yi   =  gto(ish)%coord(2)
                     zi   =  gto(ish)%coord(3)
                     ci   =  gto(ish)%coefficient(1)
                     ni   =  ncartesian(li)
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                     batch_value = ddot(ni*nj*nk*nl,gout,1,gout,1)
                     batch_values(ish,jsh,ksh,lsh) = batch_value
                     if (print_level > 1) then
                        write (*,'(a,4i6,e20.8)') 'batch',ish,jsh,ksh,lsh,batch_value 
                     endif
                  end do
               end do
            end do
         end do

!        Test correctness by checking permutational symmetry
         ierr = 0
         do lsh = 1, nshells
            do ksh = 1, lsh
               do jsh = 1, lsh
                  do ish = 1, jsh
                     if (relative_error(batch_values(ish,jsh,ksh,lsh),batch_values(jsh,ish,ksh,lsh)).gt.1.e-10) then
                        ierr = ierr + 1
                        write (*,'(a,3e20.10,4i4)') 'error in perm. symmetry #1',  &
                        batch_values(ish,jsh,ksh,lsh),batch_values(jsh,ish,ksh,lsh), &
                        batch_values(ish,jsh,ksh,lsh)-batch_values(jsh,ish,ksh,lsh), &
                        ish,jsh,ksh,lsh
                     end if
                     if (relative_error(batch_values(ish,jsh,ksh,lsh),batch_values(ish,jsh,lsh,ksh)).gt.1.e-10) then
                        ierr = ierr + 1
                        write (*,'(a,3e20.10,4i4)') 'error in perm. symmetry #2',  &
                        batch_values(ish,jsh,ksh,lsh),batch_values(ish,jsh,lsh,ksh), &
                        batch_values(ish,jsh,ksh,lsh)-batch_values(ish,jsh,lsh,ksh), &
                        ish,jsh,ksh,lsh
                     end if
                     if (relative_error(batch_values(ish,jsh,ksh,lsh),batch_values(ksh,lsh,ish,jsh)).gt.1.e-10) then
                        ierr = ierr + 1
                        write (*,'(a,3e20.10,4i4)') 'error in perm. symmetry #3',  &
                        batch_values(ish,jsh,ksh,lsh),batch_values(ksh,lsh,ish,jsh), &
                        batch_values(ish,jsh,ksh,lsh)-batch_values(ksh,lsh,ish,jsh), &
                        ish,jsh,ksh,lsh
                     end if
                  end do
               end do
            end do
         end do

         if (ierr == 0) Then
            print*," passed permutational symmetry test"
         else
            print*,ierr," failures in permutational symmetry test"
         end if

         deallocate(gto)
         deallocate(batch_values)

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Leaving test_integral_generation routine"
         return

        end subroutine test_integral_generation

end module exacorr_ao_to_mo
