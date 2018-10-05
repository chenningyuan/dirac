module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exatensor

        implicit none
        private
        public exacorr_compute_mo_integrals
        public test_integral_generation

        complex(8), parameter              :: zero=(0.D0,0.D0)

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

         class(h_space_t), pointer:: ao_space
         ! It is not possible to make an array of pointers, so we need a container
         type mo_spaces_container
              class(h_space_t), pointer :: ptr
         end type mo_spaces_container
         type(mo_spaces_container),  dimension(4) :: mo_space

         type(subspace_basis_t):: basis_ao, basis_mo(4)
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL):: nao_l,labs, ao_space_root
         integer(INTL):: nmo_l, mo_space_root(4)
         integer:: mpi_th_provided, my_rank,my_role, ao_space_id, mo_space_id(4)
         integer:: lsh, l, nl
         integer:: l_col=1
         integer:: spin,imo,jmo
         character(len=9)  :: initlabel(2,4)
         character(len=7)  :: MOCoef_str = "MOCoef_"
         character(len=2)  :: vs_str

         type(tens_rcrsv_t):: aoint_tensor
         type(tens_rcrsv_t):: moint_tensor
         type(compute_2e_ao_tensor_t) :: aoint_calculator
         type(init_mocoef_tensor_t) :: mocoef_init(2,4) ! array of initializers, first index is alpha/beta
         type(scale_with_denom_t):: denominate
         type(tens_printer_t):: tens_printer

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
            call color(lsh)%color_symmetry_ctor(lsh) ! for the moment take one shell as a color
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
         print*," process ",my_rank," registered AO vector space with ierr=",ierr
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

         !Get the MO-coeffcient basis for the four indices that are to be transformed
         ierr = 0
         imo = 1
         do i = 1, 4
            jmo     = imo + nmo(i) - 1 ! last MO for index i, whereas imo is the first
            nmo_l   = 2*int(nmo(i),8)  ! NB: imo, jmo and nmo count orbitals, nmo_l counts spinors
            call basis_mo(i)%subspace_basis_ctor(nmo_l,ierr)
            if(ierr.ne.EXA_SUCCESS) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

            !Set the functions inside this basis
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
            print*," process ",my_rank," registered MO vector space",i," with ierr=",ierr
            if(ierr.ne.EXA_SUCCESS) call quit('exatns_space_register() failed!')

            !Define separate initializers for the alpha and beta MO-coefficients
            do spin = 1, 2
               !Define a name for this MO vector space and its initializer for MOcoefficients
               if (spin == 1) write (vs_str,'(I1,A)') i,'A'
               if (spin == 2) write (vs_str,'(I1,A)') i,'B'
               initlabel(spin,i) = MOCoef_str // vs_str
               !Register the mo coefficient initializer for this space as a method
               call mocoef_init(spin,i)%init_mocoef_ctor(mo_list(imo:jmo),nmo(i),spin)
               ierr=exatns_method_register(initlabel(spin,i),mocoef_init(spin,i))
               if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
            end do
            imo = imo + nmo(i) ! make imo point to the start of the list for the next index
         end do

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')

         !Register the denominator initializer aa a method (testing only, registration will move to exacorr_cc)
         call denominate%scale_with_denom_ctor(mo_list(1:2),2,mo_list(3:nmo(1)),nmo(1)-2)
         ierr=exatns_method_register('Denominate',denominate)
         if(ierr.ne.EXA_SUCCESS) call quit('exatns_method_register() failed!')
   
         !Register the standard print method (needs to be done at the user side)
         ierr=exatns_method_register('PrintTensor',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')

!        Run ExaTENSOR within MPI_COMM_WORLD (called golbal_communicator in the MPI interface that DIRAC uses):
         ierr=exatns_start(global_communicator)
         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role)
          if(my_role.eq.EXA_DRIVER) then

           call print_date('starting ao integral calculation')

           !Create ao-integral-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'aotensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')
           if(ierr.ne.EXA_SUCCESS) call quit('Aoint_calculator failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           call print_date('finished ao integral calculation')
           call print_tensornorm (aoint_tensor)

           !Create mo-integral-tensor:
           do i = 1, 4
              mo_space_root(i)=mo_space(i)%ptr%get_root_id(ierr)
           end do
           ierr=exatns_tensor_create(moint_tensor,'motensor',mo_space_id,mo_space_root,EXA_DATA_KIND_C8)
           if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           call print_date('created mo integral tensor')
           call print_tensornorm (moint_tensor)

           call aomo_transform (aoint_tensor,moint_tensor,initlabel)

           !Destroy aoint-tensor:
           ierr=exatns_tensor_destroy(aoint_tensor)
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

           call print_date('completed mo integral tensor')

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

        subroutine aomo_transform (aoint_tensor,moint_tensor,initlabel)

!        Transform integral tensor from ao-basis to mo-basis

         implicit none

         ! The strings initlabel are used to find the appropriate initializer for the mo-coefficients
         ! This initializer has the coefficients stored (on all nodes) and just picks out the proper subset.
         character(len=*), intent(in)       :: initlabel(2,4)
         type(tens_rcrsv_t), intent(inout)  :: aoint_tensor, moint_tensor

         type(tens_rcrsv_t):: mocoef_tensor(2,4)
         type(tens_rcrsv_t):: a1int_tensor,a2int_tensor, a3int_tensor

         integer           :: ierr,i,spin
         integer           :: aoint_id(4),moint_id(4),a1int_id(4),a2int_id(4),a3int_id(4),mocoef_id(2)
         integer(INTL)     :: aoint_root(4),moint_root(4),a1int_root(4),a2int_root(4),a3int_root(4),mocoef_root(2)
         character(len=9)  :: mocoef_name

         call print_date('entered aomo_transform')

!        Determine ids and root
         call aoint_tensor%get_space_ids(aoint_id,aoint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('aoint tensor corrupted')
         call moint_tensor%get_space_ids(moint_id,moint_root,ierr)
         if(ierr.ne.EXA_SUCCESS) call quit('moint tensor corrupted')

!        Initialize tensors with MO coefficients, one for each index and alpha/beta
         mocoef_id(1)    = aoint_id(1)
         mocoef_root(1)  = aoint_root(1)
         do i = 1, 4
           mocoef_id(2)    = moint_id(i)
           mocoef_root(2)  = moint_root(i)
           do spin = 1, 2
              write (mocoef_name,'(A7,2I1)') 'mocoef_',spin,i
!             Create mo-coefficient-tensor
              ierr=exatns_tensor_create(mocoef_tensor(spin,i),mocoef_name, &
               mocoef_id,mocoef_root,EXA_DATA_KIND_C8)
              if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_create() failed!')
              ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
!             Initialize mocoef-tensor 
              ierr=exatns_tensor_init(mocoef_tensor(spin,i),initlabel(spin,i))
              if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_init() failed!')
              ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
              print*,"created and initialized mocoef_tensor ",spin,i
              call print_tensornorm (mocoef_tensor(spin,i))
           end do
         end do

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
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit('exatns_sync() failed!')
         ierr=exatns_tensor_create(a2int_tensor,'a2int',a2int_id,a2int_root,EXA_DATA_KIND_C8)
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         call print_date('constructed a1 and a2')
         ierr=exatns_tensor_contract(a1int_tensor,mocoef_tensor(1,1),aoint_tensor,"D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)")
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=exatns_tensor_contract(a2int_tensor,a1int_tensor,mocoef_tensor(1,2),"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_init(a1int_tensor,zero)
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_contract(a1int_tensor,mocoef_tensor(2,1),aoint_tensor,"D(i,q,r,s)=R+(p,i)*L(p,q,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_contract(a2int_tensor,a1int_tensor,mocoef_tensor(2,2),"D(i,j,r,s)+=L(i,q,r,s)*R(q,j)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=22; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         call print_date('first ht done')

!        Destruct a1 integrals:
         ierr=exatns_tensor_destroy(a1int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=25; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

!        Do 3rd and 4th quarter-transformation
         ierr=exatns_tensor_create(a3int_tensor,'a3int',a3int_id,a3int_root,EXA_DATA_KIND_C8)
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_init(a3int_tensor,zero)
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         call print_date('constructed a3')
         ierr=exatns_tensor_contract(a3int_tensor,mocoef_tensor(1,3),a2int_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_contract(moint_tensor,a3int_tensor,mocoef_tensor(1,4),"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_init(a3int_tensor,zero)
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_contract(a3int_tensor,mocoef_tensor(2,3),a2int_tensor,"D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=32; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_contract(moint_tensor,a3int_tensor,mocoef_tensor(2,4),"D(i,j,k,l)+=L(i,j,k,s)*R(s,l)")
         if(ierr.ne.EXA_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=42; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         call print_date('second ht done')

!        Destruct a2 and a3 integrals:
         ierr=exatns_tensor_destroy(a2int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=35; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
         ierr=exatns_tensor_destroy(a3int_tensor); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif
         ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
!        Destruct mocoef tensors:
         do i = 1, 4
           do spin = 1, 2
            ierr=exatns_tensor_destroy(mocoef_tensor(spin,i)); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif
            ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
            end do
         end do
         call print_date('leaving aomo_transform')

        end subroutine aomo_transform
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

        subroutine print_tensornorm (tens)

          type(tens_rcrsv_t):: tens, scalar
          integer ierr,rank
 
          ! create a rank-0 tensor
          ierr=exatns_tensor_create(scalar,'scalar',EXA_DATA_KIND_C8)
          if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
          ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')

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
          ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit('exatns_sync() failed!')
          ierr=exatns_tensor_traverse(scalar,'PrintTensor')
          if(ierr.ne.EXA_SUCCESS) call quit('exatns_tensor_traverse() failed!')
          ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit('exatns_sync() failed!')
          ierr=exatns_tensor_destroy(scalar); if(ierr.ne.EXA_SUCCESS) then; ierr=45; return; endif
          ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit('exatns_sync() failed!')

        end subroutine print_tensornorm


end module exacorr_ao_to_mo
