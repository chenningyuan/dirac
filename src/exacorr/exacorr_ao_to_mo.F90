module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exatensor

        implicit none
        private
        public compute_ao_integrals
        public compute_mo_integrals
        public get_integral_tensor
        public test_integral_generation

       contains

        subroutine get_mocoef_tensor (nao, nmo, mo_list, mocoefa_tensor, mocoefb_tensor)

!This subroutine reads the desired subset of mo coefficients from DIRAC

!        Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 30 C)

         use dirac_MO
         use tensor_algebra
         use talsh

         implicit none
         integer, intent(in) :: nao          ! the length of the ao basis
         integer, intent(in) :: nmo          ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         type(talsh_tens_t), intent(inout) :: mocoefa_tensor,mocoefb_tensor

         integer             :: j, k, ierr
         type(MO)            :: qmo
         integer             :: ifrp = 1     ! we do not use symmetry, so ifrp (g/u) is always 1
         integer             :: nz = 4       ! we do not use symmetry, so nz is always 4
         integer             :: tens_rank, mocoef_dims(1:2)

         complex(8), pointer, contiguous:: camo_tens(:,:),cbmo_tens(:,:)
         type(C_PTR):: body_p


         ierr=talsh_tensor_dimensions(mocoefa_tensor,tens_rank,mocoef_dims)
         ierr=talsh_tensor_get_body_access(mocoefa_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,camo_tens,mocoef_dims(1:2)) ! to use camo_tens as a regular Fortran 2d array
         ierr=talsh_tensor_get_body_access(mocoefb_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,cbmo_tens,mocoef_dims(1:2)) ! to use cbmo_tens as a regular Fortran 2d array

         call get_dirac_mo_coefficients (ifrp,qmo,mo_list,nmo)

!        Copy the coefficients into the alpha and beta tensors, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         do k = 1, nmo
            do j = 1, nao
                camo_tens(j,k) = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
                cbmo_tens(j,k) = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
                camo_tens(j,k+nmo) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
                cbmo_tens(j,k+nmo) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
            end do
         end do
         call delete (qmo)

         return
        end subroutine get_mocoef_tensor

        subroutine compute_ao_integrals (nao,ao_tens)

!This subroutine drives the generation of AO integrals by InteRest
!
         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri

         implicit none
         integer, intent(in) :: nao ! number of basis functions
         complex(8),intent(out) :: ao_tens(nao,nao,nao,nao)

         integer :: ierr
         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl, ijkl2
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
!lv      code marked lv is tested and functional but not needed yet as we don't have the code to transform the additional integral classes
!lv      real(8) :: gout2(21*21*5*21*21*5) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         real(8), parameter :: fijkl = 1.0
         integer            :: nshells  ! number of AO shells
         type(basis_func_info_t), allocatable :: gto(:)   ! array with basis function information

         call fetch_gtos(1,nao,gto,nshells)

         ![eric] OpenMP task-based parallelism
         !$omp parallel
         !$omp single
         loff = 0
         do lsh = 1, nshells
            ll   =  gto(lsh)%orb_momentum
            el   =  gto(lsh)%exponent(1)
            xl   =  gto(lsh)%coord(1)
            yl   =  gto(lsh)%coord(2)
            zl   =  gto(lsh)%coord(3)
            cl   =  gto(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            koff = 0
            do ksh = 1, nshells
               lk   =  gto(ksh)%orb_momentum
               ek   =  gto(ksh)%exponent(1)
               xk   =  gto(ksh)%coord(1)
               yk   =  gto(ksh)%coord(2)
               zk   =  gto(ksh)%coord(3)
               ck   =  gto(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               joff = 0
               do jsh = 1, nshells
                  lj   =  gto(jsh)%orb_momentum
                  ej   =  gto(jsh)%exponent(1)
                  xj   =  gto(jsh)%coord(1)
                  yj   =  gto(jsh)%coord(2)
                  zj   =  gto(jsh)%coord(3)
                  cj   =  gto(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  ioff = 0
                  do ish = 1, nshells
                     li   =  gto(ish)%orb_momentum
                     ei   =  gto(ish)%exponent(1)
                     xi   =  gto(ish)%coord(1)
                     yi   =  gto(ish)%coord(2)
                     zi   =  gto(ish)%coord(3)
                     ci   =  gto(ish)%coefficient(1)
                     ni   =  ncartesian(li)
!                    call interest and get only the (LL|LL) type integrals
!                    output order of eri is (c,d,a,b), so input k,l first and then i,j to get the order that we want

                     ![eric]
                     !$omp task
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )

!lv                  call interest and get (LL|LL), (SS|LL) and (SS|SS)  type integrals
!lv                  output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
!lv                  call interest_eri_rkb('ssss',fijkl,gout2,&
!lv                               lk,ek,xk,yk,zk,ck,&
!lv                               ll,el,xl,yl,zl,cl,&
!lv                               li,ei,xi,yi,zi,ci,&
!lv                               lj,ej,xj,yj,zj,cj )
                     ijkl = 0
                     ijkl2 = 1
                     do l = 1, nl
                        do k = 1, nk
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = ijkl + 1
                                 ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
!lv                              ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout2(ijkl2),0.D0)
!lv                              ijkl2 = ijkl2 + 5
                              end do
                           end do
!lv                        ijkl2 = ijkl2 + ni*nj*5*4
                        end do
                     end do
                     !$omp end task

                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do
         !$omp end single nowait
         !$omp end parallel

         deallocate(gto)

        end subroutine compute_ao_integrals

        subroutine aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)

!        Transform integral tensor from ao-basis to mo-basis

         use tensor_algebra
         use talsh

         implicit none

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         type(talsh_tens_t), intent(inout)  :: aoint_tensor, moint_tensor

         type(talsh_tens_t):: mocoefa_tensor(4),mocoefb_tensor(4)
         type(talsh_tens_t):: a1int_tensor,a2int_tensor, a3int_tensor
         integer:: aoint_dims(1:4),mocoef_dims(1:2),a1int_dims(1:4)
         integer:: a2int_dims(1:4),a3int_dims(1:4),moint_dims(1:4)

         complex(8):: zero=(0.D0,0.D0)
         integer   :: ierr, tens_rank, i, imo, nao

         call print_date('entered aomo_transform')
!        Get tensor dimension of input and output tensors
         ierr = talsh_tensor_dimensions(aoint_tensor,tens_rank,aoint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: aoint tensor corrupted'
         ierr = talsh_tensor_dimensions(moint_tensor,tens_rank,moint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: moint tensor corrupted'

!        Initialize tensors with MO coefficients
         nao = aoint_dims(1)
         imo = 0
         do i = 1, 4
            mocoef_dims(1) = aoint_dims(i)
            if (mocoef_dims(1) .ne. nao)      stop 'program error: inconsistent dimensions'
            mocoef_dims(2) = moint_dims(i)
            if (mocoef_dims(2) .ne. 2*nmo(i)) stop 'program error: inconsistent dimensions'
            ierr=talsh_tensor_construct(mocoefa_tensor(i),C8,mocoef_dims,init_val=zero)
            ierr=talsh_tensor_construct(mocoefb_tensor(i),C8,mocoef_dims,init_val=zero)
            call get_mocoef_tensor (nao, nmo(i), mo_list(imo+1:imo+nmo(i)), mocoefa_tensor(i), mocoefb_tensor(i))
            imo = imo + nmo(i)
         end do
         call print_date('initialized mo-coefficients')

!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         a1int_dims(1)   = moint_dims(1)
         a1int_dims(2:4) = aoint_dims(2:4)
         a2int_dims(1:2) = moint_dims(1:2)
         a2int_dims(3:4) = aoint_dims(3:4)
         a3int_dims(1:3) = moint_dims(1:3)
         a3int_dims(4)   = aoint_dims(4)

!        Do 1st and 2nd quarter-transformation
         ierr=talsh_tensor_construct(a1int_tensor,C8,a1int_dims,init_val=zero)
         ierr=talsh_tensor_construct(a2int_tensor,C8,a2int_dims,init_val=zero)
         call print_date('constructed a1 and a2')
         ierr=talsh_tensor_contract("D(i,q,r,s)+=R+(p,i)*L(p,q,r,s)",a1int_tensor,mocoefa_tensor(1),aoint_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=11; return; endif
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefa_tensor(2))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         ierr=talsh_tensor_init(a1int_tensor,zero)
         ierr=talsh_tensor_contract("D(i,q,r,s)=R+(p,i)*L(p,q,r,s)",a1int_tensor,mocoefb_tensor(1),aoint_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=12; return; endif
         ierr=talsh_tensor_contract("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefb_tensor(2))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=22; return; endif
         call print_date('first ht done')

!        Destruct a1 and a1 integrals:
         ierr=talsh_tensor_destruct(a1int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=25; return; endif

!        Do 3rd and 4th quarter-transformation
         ierr=talsh_tensor_construct(a3int_tensor,C8,a3int_dims,init_val=zero)
         call print_date('constructed a3')
         ierr=talsh_tensor_contract("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefa_tensor(3),a2int_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefa_tensor(4))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=talsh_tensor_init(a3int_tensor,zero)
         ierr=talsh_tensor_contract("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefb_tensor(3),a2int_tensor)
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=32; return; endif
         ierr=talsh_tensor_contract("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefb_tensor(4))
         if(ierr.ne.TALSH_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=42; return; endif
         call print_date('second ht done')

!        Destruct a2 and a3 integrals:
         ierr=talsh_tensor_destruct(a2int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=35; return; endif
         ierr=talsh_tensor_destruct(a3int_tensor); if(ierr.ne.TALSH_SUCCESS) then; ierr=45; return; endif
         call print_date('leaving aomo_transform')

         return

        end subroutine aomo_transform

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,nmo,mo_list)


         use tensor_algebra
         use talsh

         implicit none

         complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         type(talsh_tens_t), intent(inout)  :: aoint_tensor,integral_tensor
         integer, intent(in)                :: nmo(4), mo_list(:)
         integer, intent(in)                :: antisymmetrize

         integer                            :: ierr, moint_dims(4)
         type(talsh_tens_t)                 :: moint_tensor

         moint_dims = 2 * nmo
         ierr=talsh_tensor_construct(moint_tensor,C8,moint_dims,init_val=ZERO)
!        Transform integrals
         call aomo_transform (aoint_tensor,moint_tensor,nmo,mo_list)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=talsh_tensor_destruct(moint_tensor)

        end subroutine get_integral_tensor

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         use tensor_algebra
         use talsh

         implicit none

         complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

         type(talsh_tens_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         integer, intent(in)                :: antisymmetrize

         integer:: mulliken_dims(1:4),dirac_dims(1:4)

         integer   :: ierr, tens_rank, mismatch

!        work-around: addition with permutation is not implemented
         type(talsh_tens_t):: one_tensor, minusone_tensor
         integer:: scalar_dims(1)
         scalar_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,scalar_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(minusone_tensor,C8,scalar_dims(1:0),init_val=MINUS_ONE)
!        end work-round

!        get dimensions of the tensors and check consistency
         ierr = talsh_tensor_dimensions(mulliken_tensor,tens_rank,mulliken_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input tensor corrupted in mulliken_to_dirac_sort'
         ierr = talsh_tensor_dimensions(dirac_tensor,tens_rank,dirac_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'output tensor corrupted in mulliken_to_dirac_sort'
         mismatch = 0
         if (mulliken_dims(1) .ne. dirac_dims(1)) mismatch = mismatch + 1
         if (mulliken_dims(2) .ne. dirac_dims(3)) mismatch = mismatch + 1
         if (mulliken_dims(3) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (mulliken_dims(4) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 12 .and. dirac_dims(1) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 34 .and. dirac_dims(3) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (mismatch.gt.0) stop 'dimension mismatch in mulliken_to_dirac_sort'

         ierr=talsh_tensor_init(dirac_tensor)
!shorter ierr=talsh_tensor_add("D(p,q,r,s)+=M(p,r,q,s)",dirac_tensor,mulliken_tensor)
         ierr=talsh_tensor_contract("D(p,q,r,s)+=M(p,r,q,s)*C()",dirac_tensor,mulliken_tensor,one_tensor)

         if (antisymmetrize .eq. 12) then
!shorter     ierr=talsh_tensor_add("D(p,q,r,s)+=M(q,r,p,s)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
             ierr=talsh_tensor_contract("D(p,q,r,s)+=M(q,r,p,s)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
         elseif (antisymmetrize .eq. 34) then
!shorter     ierr=talsh_tensor_add("D(p,q,r,s)+=M(p,s,q,r)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
              ierr=talsh_tensor_contract("D(p,q,r,s)+=M(p,s,q,r)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
         end if

!        work-around: addition with permutation is not implemented
         ierr=talsh_tensor_destruct(one_tensor)
         ierr=talsh_tensor_destruct(minusone_tensor)
!        end work-round

        end subroutine mulliken_to_dirac_sort


        subroutine fetch_gtos(first_ao,nao,gto,nshells)

         use exacorr_datatypes
         use exacorr_global

         implicit none

         integer, intent(in):: first_ao, nao
         type(basis_func_info_t), intent(inout), allocatable :: gto(:)
         integer, intent(out) :: nshells

         integer :: ish, first_shell, last_shell

         first_shell = get_shell_index (first_ao)
         last_shell  = get_shell_index (first_ao+nao-1)
         nshells = last_shell - first_shell + 1

         allocate (gto(nshells))
         do ish = 1, nshells
            gto(ish) = get_gto(first_shell+ish-1)
         end do

        end subroutine fetch_gtos

#if defined (VAR_MPI)
        subroutine compute_mo_integrals (nmo, mo_list)

!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor

!        Written by Lucas Visscher, march 2018

         use dirac_MO
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use module_interest_interface
         use module_interest_eri
         use exatensor
         use mpi
         use, intrinsic:: ISO_C_BINDING

         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: i, ierr
         integer             :: number_aos                 ! the size of the atomic orbital basis
         integer             :: nshells                    ! the number of shells (orbitals with the same center and l-value)
         type(basis_func_info_t), allocatable :: gto(:)    ! array with information about the shells

         type(subspace_basis_t):: basis_ao
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL):: nao_l,labs, ao_space_root
         class(h_space_t), pointer:: ao_space
         integer:: mpi_th_provided, my_rank,my_role, ao_space_id
         integer:: nao, lsh, l, nl

         type(tens_rcrsv_t):: aoint_tensor
!later   type(tens_rcrsv_t):: moint_tensor
         type(compute_2e_ao_tensor_t) :: aoint_calculator

!        Make some noise so that we know we are here
         call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
         print*, " Process", my_rank," arrived in compute_mo_integrals_exa routine"

!        Initialize libraries
         call interest_initialize
         call initialize_global_data

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)
         nao_l   = nao ! copy of nao to a long integer

!        Get the basis set information
         call fetch_gtos(1,nao,gto,nshells)

!        Temporary: leave here because the exatensor library is not yet working
         print*, " WARNING: we quit the index transform because exatensor is not yet operational"
         call delete_global_data
         print*, " Process ",my_rank, " leaving compute_mo_integrals_exa routine"
         return

         !Create a basis for a vector space:
         call basis_ao%subspace_basis_ctor(nao_l,ierr)
         !This needs to be an array that is persistent on all nodes
         allocate(color(nao))
         if(ierr.ne.0) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         labs = 0
         do lsh = 1, nshells !set basis functions
            nl   =  ncartesian(gto(lsh)%orb_momentum)
            do  l = 1, nl
                labs = labs + 1
                call color(labs)%color_symmetry_ctor(lsh) ! for the moment take one shell as a color
                call basis_ao%set_basis_func(labs,BASIS_ABSTRACT,ierr,symm=color(labs))
                if(ierr.ne.0) call quit('subspace_basis_t.set_basis_func() failed!')
            end do
         enddo

         call basis_ao%finalize(ierr)
         if(ierr.ne.0) call quit('subspace_basis_t.finalize() failed!')

         !Register a vector space:
         ierr=exatns_space_register('AO_space',basis_ao,ao_space_id,ao_space)
         if(ierr.ne.0) call quit('exatns_space_register() failed!')

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)

!        Run ExaTENSOR within MPI_COMM_WORLD:
         ierr=exatns_start(MPI_COMM_WORLD)
         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role)
          if(my_role.eq.EXA_DRIVER) then

           !Create ao-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'ao_tensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')

           !Destroy ao-tensor:
           ierr=exatns_tensor_destroy(aoint_tensor)
          end if

          !Stop ExaTENSOR runtime:
          ierr=exatns_stop()

         else
          write(*,*) ' Process ',my_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Process ",my_rank, " leaving compute_mo_integrals_exa routine"
         return
        end subroutine compute_mo_integrals

#else
        subroutine compute_mo_integrals (nmo, mo_list)
         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         print*,  'called exatensor version of compute_mo_integrals in a serial version'
         print*,  'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine compute_mo_integrals

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

         integer             :: nao, number_aos            ! the size of the atomic orbital basis

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

!        Make some noise so that we know we are here
         print*, " Testing AO integral generation by calculating all integrals once"

!        Initialize libraries
         call interest_initialize(print_level.gt.1)
         if (print_level > 0) print*, " Initialized interest library"
         call initialize_global_data
         if (print_level > 0) print*, " Initialized global data"

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)

!        Get AOs and allocate array to store batch values
         call fetch_gtos(1,nao,gto,nshells)
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

        real function relative_error(quantity_1,quantity_2)

        real(8),intent(in):: quantity_1,quantity_2
        real(8)           :: error

        error = quantity_1 - quantity_2

        if (abs(quantity_1).gt.1.0) then
           relative_error = abs (error / quantity_1)
        else
           relative_error = abs (error)
        end if

        end function relative_error
           
        subroutine print_date(text)
         character(len=10) :: datex
         character(len=8)  :: timex
         character(len=*)  :: text
         integer           :: lentext
         lentext = len(text)
         call daytime(datex,timex)
         write (*,'(1x,a10,1x,a8,1x,a)') datex,timex,text(1:lentext)
        end subroutine print_date

end module exacorr_ao_to_mo
