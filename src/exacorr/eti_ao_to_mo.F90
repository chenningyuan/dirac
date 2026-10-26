module eti_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use tensor_algebra
        use exacorr_tensor_interface
        use eti_common_routines, only : tensor_norm2
        use exacorr_utils
        use exacorr_scaling_parameters
        use, intrinsic:: ISO_C_BINDING

        implicit none

        private
        public eti_compute_mo_integrals
        public compute_ao_integrals
        public get_integral_tensor
        public do_ht1
        public do_ht2
        public ao2mo_ind
        public ao2mo_vec
        public mulliken_to_dirac_sort

        interface compute_ao_integrals
           module procedure compute_ao_integrals
           module procedure compute_selected_ao_integrals
        end interface

        type, public  :: ao2mo_info
            ! type to pack information needed in the ao to mo transformation steps
            integer             :: &
               nao, &               ! number of aos active in this transformation
               nmo, &               ! number of mos active in this transformation
               first_ao, &          ! index of first active AO
               last_ao              ! index of last active AO

            integer, pointer :: &
               mo_list (:)          ! indices of the active MOs
         
         end type ao2mo_info

       contains
!-----------------------------------------------------------

        subroutine eti_compute_mo_integrals (nmo, mo_list, buff, executor, print_level, rcw, th_cholesky)

!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
!The integrals are written to MDCINT, to allow for testing with unmodified DIRAC correlation modules

!        Written by Lucas Visscher, winter 2016/2017

         use exacorr_mo
         use exacorr_datatypes
         use exacorr_global
         use module_interest_interface
         use eti_CD_routines
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
         use interface_to_mpi
#endif

         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         integer, intent(in) :: buff   ! buffer size for talsh in MB
         integer, intent(in) :: executor   ! executor
         integer, intent(in) :: print_level  
         integer, intent(in) :: rcw          ! real or complex 
         real(8), intent(in) :: th_cholesky

         integer(INTD)       :: ierr
         integer             :: nao, my_MPI_rank
         type(tensor_t)  :: aoint_tensor
         type(tensor_t)  :: moint_tensor
         integer(C_INT)      :: aoint_dims(1:4),moint_dims(1:4)
         integer             :: m, i, imo
         complex(8), pointer, contiguous :: ao_tens(:,:,:,:),mo_tens(:,:,:,:)
         type(C_PTR)                     :: body_p

         type(tensor_t), allocatable :: chol_tensor(:), temp(:)
         type(tensor_t), allocatable :: upq(:), urs(:)
         integer(C_INT)                  :: dim2(1:2)

         type(tensor_t)              :: chol_tens, temp_tens
         type(tensor_t)              :: upq_tens, urs_tens
         integer(C_INT)              :: dim3(1:3), dim_chol(1:3)

         type(ao2mo_info)            :: aomo_spaces(4)

!        Make some noise so that we know we are here
         call print_date('Entered compute_mo_integrals routine')

!        Initialize global data
         call initialize_global_data
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
#else
         my_MPI_rank = 0
#endif
         if (my_MPI_rank /= 0) return ! TALSH is not parallelized, leave the work to master
         call print_date('Initialized global data')

!        The experimental Cholesky code works only with TALSH as it needs slicing and reshaping of tensors
         if (th_cholesky > 0) then
             write(*,*) ' Cholesky feature requested, thus always using executor 2 (TALSH)'
             ierr = set_default_executor(2)
         end if

         call eti_executor_initialize(buff,get_default_executor())

         call interest_initialize(.true.)
         call print_date('Initialized InteRest library')

!        Retrieve basis set information
         nao = get_nao()

!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         aoint_dims = nao
         moint_dims = nmo

         if (th_cholesky.le.0) then
!          Get tensor for AO integrals
           ierr=construct_tensor_t(aoint_tensor,C8,aoint_dims,init_val=ZERO)
!          Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
           ierr=direct_data_pointer(aoint_tensor,body_p)
           call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
           call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
           call print_date('Calculated AO integrals')

!          Transform to MO basis
           ierr=construct_tensor_t(moint_tensor,C8,moint_dims,init_val=ZERO)
           aomo_spaces(:)%nao = nao
           aomo_spaces(:)%first_ao = 1
           aomo_spaces(:)%last_ao = nao
           aomo_spaces(:)%nmo = nmo
           imo = 0
           do i = 1,4
              allocate(aomo_spaces(i)%mo_list(aomo_spaces(i)%nmo))
              aomo_spaces(i)%mo_list = mo_list(imo+1:imo+nmo(i))
              imo = imo + nmo(i)
           enddo
           call aomo_transform (aoint_tensor,moint_tensor,aomo_spaces)
           if (print_level.gt.2) then
             print*, "AOint = ",tensor_norm2(aoint_tensor)
             print*, "MOint = ",tensor_norm2(moint_tensor)
             print*, "additional transform to get anti-symmetrized integrals"
           end if

         else
            call print_date('Cholesky method 42')
            !determine cholesky vectors for ao integrals
            call decompose_cholesky_eti (chol_tens, th_cholesky, print_level)
!           ierr = tensor_dimensions(chol_tens,i,dim_chol)
            ierr = get_nmodes(chol_tens,i)
            ierr = get_extents_static(chol_tens,dim_chol)

            !  transform upq
            dim3=dim_chol
            dim3(1)=nmo(1)
            ierr=construct_tensor_t(temp_tens,C8,dim3,init_val=ZERO)
            dim3(1:2)=nmo(1:2)
            ierr=construct_tensor_t(upq_tens,C8,dim3,init_val=ZERO)

            call ao2mo_ind(chol_tens,temp_tens,nmo(1),mo_list(1:nmo(1)),11)
            call ao2mo_ind(temp_tens,upq_tens,nmo(2),mo_list(nmo(1)+1:sum(nmo(1:2))),21)
            ierr=initialize_tensor_t(temp_tens)
            call ao2mo_ind(chol_tens,temp_tens,nmo(1),mo_list(1:nmo(1)),12)
            call ao2mo_ind(temp_tens,upq_tens,nmo(2),mo_list(nmo(1)+1:sum(nmo(1:2))),22)

            ierr=destroy_tensor_t(temp_tens)

            if (print_level.gt.2) call print_date(' upq done ')

            !  transform urs
            dim3=dim_chol
            dim3(1)=nmo(3)
            ierr=construct_tensor_t(temp_tens,C8,dim3,init_val=ZERO)
            dim3(1:2)=nmo(3:4)
            ierr=construct_tensor_t(urs_tens,C8,dim3,init_val=ZERO)

            call ao2mo_ind(chol_tens,temp_tens,nmo(3),mo_list(sum(nmo(1:2))+1:sum(nmo(1:3))),11)
            call ao2mo_ind(temp_tens,urs_tens,nmo(4),mo_list(sum(nmo(1:3))+1:sum(nmo)),21)
            ierr=initialize_tensor_t(temp_tens)
            call ao2mo_ind(chol_tens,temp_tens,nmo(3),mo_list(sum(nmo(1:2))+1:sum(nmo(1:3))),12)
            call ao2mo_ind(temp_tens,urs_tens,nmo(4),mo_list(sum(nmo(1:3))+1:sum(nmo)),22)

            ierr=destroy_tensor_t(temp_tens)
            ierr=destroy_tensor_t(chol_tens)

            if (print_level.gt.2) call print_date(' urs done ')

            ierr=construct_tensor_t(moint_tensor,C8,moint_dims,init_val=ZERO)
            ierr=t_einsum("T(p,q,r,s)+=K(p,q,u)*L(r,s,u)",moint_tensor,upq_tens,urs_tens)
            if (ierr.ne.0) stop 'program error: cholesky not working'

            ierr=destroy_tensor_t(upq_tens)
            ierr=destroy_tensor_t(urs_tens)

         end if
         call print_date(' 2-el Integrals done ')

!        Write to old format MO-integrals file to allow for easy debugging of the integral transform
         ierr=direct_data_pointer(moint_tensor,body_p)
         call c_f_pointer(body_p,mo_tens,moint_dims(1:4))
         call write_to_mdcint (nmo,mo_tens,rcw)

         if (print_level.gt.8) then
            print*, "additional transform to get anti-symmetrized integralss"
            call get_integral_tensor (aoint_tensor,moint_tensor,12,aomo_spaces)
            print*, "MOint-antisym = ",tensor_norm2(moint_tensor)
         end if

!        Destruct remaining tensors and shut down library
         ierr=destroy_tensor_t(aoint_tensor)
         ierr=destroy_tensor_t(moint_tensor)

!        shutdown tensor libraries
         call eti_executor_shutdown

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         call print_date('Leaving compute_mo_integrals routine')

        end subroutine eti_compute_mo_integrals

        subroutine get_mocoef_tensor (first_ao, last_ao, nmo, mo_list, mocoefa_tensor, mocoefb_tensor)

!This subroutine reads the desired subset of mo coefficients from DIRAC

!        Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 30 C)

         use exacorr_mo
         use exacorr_global

         implicit none
         integer, intent(in) :: first_ao, last_ao ! range of the ao basis
         integer, intent(in) :: nmo          ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         type(tensor_t), intent(inout) :: mocoefa_tensor,mocoefb_tensor

         integer(INTD)       :: ierr, tens_rank
         type(cmo)           :: cspinor
         integer(INTD)       :: mocoef_dims(1:2)

         complex(8), pointer, contiguous:: camo_tens(:,:),cbmo_tens(:,:)
         type(C_PTR):: body_p

!        ierr=tensor_dimensions(mocoefa_tensor,tens_rank,mocoef_dims)
         ierr = get_nmodes(mocoefa_tensor, tens_rank)
         ierr = get_extents_static(mocoefa_tensor, mocoef_dims)

         ierr=direct_data_pointer(mocoefa_tensor,body_p)
         call c_f_pointer(body_p,camo_tens,mocoef_dims(1:2)) ! to use camo_tens as a regular Fortran 2d array
         ierr=direct_data_pointer(mocoefb_tensor,body_p)
         call c_f_pointer(body_p,cbmo_tens,mocoef_dims(1:2)) ! to use cbmo_tens as a regular Fortran 2d array

         call get_mo_coefficients (cspinor,mo_list,nmo)

!        Copy the coefficients into the alpha and beta tensors, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         camo_tens = cspinor%coeff(first_ao:last_ao,:,1)
         cbmo_tens = cspinor%coeff(first_ao:last_ao,:,2)

        end subroutine get_mocoef_tensor

        subroutine compute_ao_integrals (nao,ao_tens_complex,ao_tens_real)

!This subroutine drives the generation of AO integrals by InteRest
!
         implicit none
         integer, intent(in) :: nao           ! number of basis functions
         complex(8),intent(out),optional :: ao_tens_complex(nao,nao,nao,nao)
         real(8),intent(out),optional    :: ao_tens_real(nao,nao,nao,nao)

         call compute_selected_ao_integrals ((/1,1,1,1/),(/nao,nao,nao,nao/), &
                                             ao_tens_complex,ao_tens_real)

        end subroutine compute_ao_integrals

        subroutine compute_selected_ao_integrals (first_ao,nao,ao_tens_complex,ao_tens_real)

!This subroutine drives the generation of AO integrals by InteRest.
!The shell quartets are distributed over OpenMP threads: each quartet
!writes its own block of the output, InteRest is initialized per thread.
!
         use exacorr_ao
         use exacorr_global
         use exacorr_eri

         implicit none
         integer, intent(in) :: first_ao(4)      ! first ao for each index
         integer, intent(in) :: nao(4)           ! number of basis functions for each index
         integer             :: basis_angular    ! 1=cartesian, 2=spherical
         complex(8),intent(out),optional :: ao_tens_complex(nao(1),nao(2),nao(3),nao(4))
         real(8),intent(out),optional    :: ao_tens_real(nao(1),nao(2),nao(3),nao(4))

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ijkl, nmax
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8), allocatable :: gout(:)
         real(8), parameter :: fijkl = 1.0
         integer            :: nshells(4)                   ! number of AO shells
         integer, allocatable :: off1(:), off2(:), off3(:), off4(:) ! first function of each shell - 1
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)   ! arrays with basis function information
         logical complex_output

         if (present(ao_tens_complex)) then
            complex_output = .true.
         else if (present(ao_tens_real)) then
            complex_output = .false.
         else
            call quit ('no output array given in compute_selected_ao_integrals')
         end if

         basis_angular = get_basis_angular()
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))
         call shell_offsets(gto1,nshells(1),off1)
         call shell_offsets(gto2,nshells(2),off2)
         call shell_offsets(gto3,nshells(3),off3)
         call shell_offsets(gto4,nshells(4),off4)
         nmax = 1
         do ish = 1, nshells(1)
            nmax = max(nmax,off1(ish+1)-off1(ish))
         end do
         do ish = 1, nshells(2)
            nmax = max(nmax,off2(ish+1)-off2(ish))
         end do
         do ish = 1, nshells(3)
            nmax = max(nmax,off3(ish+1)-off3(ish))
         end do
         do ish = 1, nshells(4)
            nmax = max(nmax,off4(ish+1)-off4(ish))
         end do

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k,l,ish,jsh,ksh,lsh,ijkl,li,lj,lk,ll,ni,nj,nk,nl, &
!$OMP         xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,gout)
         call interest_initialize(.false.)
         allocate(gout(nmax**4))
!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
         do lsh = 1, nshells(4)
            do ksh = 1, nshells(3)
               ll   =  gto4(lsh)%orb_momentum
               xl   =  gto4(lsh)%coord(1)
               yl   =  gto4(lsh)%coord(2)
               zl   =  gto4(lsh)%coord(3)
               nl   =  off4(lsh+1) - off4(lsh)
               lk   =  gto3(ksh)%orb_momentum
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               nk   =  off3(ksh+1) - off3(ksh)
               do jsh = 1, nshells(2)
                  lj   =  gto2(jsh)%orb_momentum
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  nj   =  off2(jsh+1) - off2(jsh)
                  do ish = 1, nshells(1)
                     li   =  gto1(ish)%orb_momentum
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ni   =  off1(ish+1) - off1(ish)
!                    call interest and get the (LL|LL) type integrals
!                    output order of eri is (c,d,a,b), so input k,l first and then i,j to get the order that we want
                     call compute_eri('llll',fijkl,gout(1:ni*nj*nk*nl),&
                                  lk,gto3(ksh)%exponent,xk,yk,zk,gto3(ksh)%coefficient,&
                                  ll,gto4(lsh)%exponent,xl,yl,zl,gto4(lsh)%coefficient,&
                                  li,gto1(ish)%exponent,xi,yi,zi,gto1(ish)%coefficient,&
                                  lj,gto2(jsh)%exponent,xj,yj,zj,gto2(jsh)%coefficient )

                     ijkl = 0
                     do l = 1, nl
                        do k = 1, nk
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = ijkl + 1
                                 if (complex_output) then
                                    ao_tens_complex(i+off1(ish),j+off2(jsh),k+off3(ksh),l+off4(lsh)) = &
                                       dcmplx(gout(ijkl),0.D0)
                                 else
                                    ao_tens_real(i+off1(ish),j+off2(jsh),k+off3(ksh),l+off4(lsh)) = gout(ijkl)
                                 end if
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
!$OMP END DO
         deallocate(gout)
!$OMP END PARALLEL

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)

        contains

         subroutine shell_offsets(gto,nsh,off)
         type(basis_func_info_t), intent(in) :: gto(:)
         integer, intent(in)                 :: nsh
         integer, allocatable, intent(out)   :: off(:)
         integer                             :: s

         allocate(off(nsh+1))
         off(1) = 0
         do s = 1, nsh
            off(s+1) = off(s) + nfunctions(gto(s)%orb_momentum,basis_angular)
         end do
         end subroutine shell_offsets

        end subroutine compute_selected_ao_integrals

        subroutine get_integral_tensor (aoint_tensor,integral_tensor,antisymmetrize,aomo_spaces)

         implicit none

         type(tensor_t), intent(inout)  :: aoint_tensor,integral_tensor
         type(ao2mo_info), intent(in)   :: aomo_spaces(4)
         integer, intent(in)            :: antisymmetrize

         integer(INTD)                  :: ierr, moint_dims(4)
         type(tensor_t)                 :: moint_tensor

         moint_dims = aomo_spaces(1:4)%nmo
         ierr=construct_tensor_t(moint_tensor,C8,moint_dims,init_val=ZERO)
         if (ierr.ne.0) stop 'get_integral_tensor: tensor construction failed'
!        Transform integrals
         call aomo_transform (aoint_tensor,moint_tensor,aomo_spaces)
!        Switch to Dirac order and antisymmetrize
         call mulliken_to_dirac_sort (moint_tensor,integral_tensor,antisymmetrize)
         ierr=destroy_tensor_t(moint_tensor)

        end subroutine get_integral_tensor

        subroutine aomo_transform (aoint_tensor,moint_tensor,aomo_spaces)

!        Transform integral tensor from ao-basis to mo-basis(all 4 indices)

         use exacorr_global

         implicit none

         type(ao2mo_info), intent(in)   :: aomo_spaces(4)
         type(tensor_t), intent(inout)  :: aoint_tensor, moint_tensor

         type(tensor_t):: a2int_tensor
         integer(INTD) :: aoint_dims(4),moint_dims(4),a2int_dims(4)
 
         integer(INTD) :: ierr, tens_rank, first_ao(2), last_ao(2)
 
!        Get tensor dimension of input and output tensors
         ierr = get_nmodes(aoint_tensor,tens_rank)
         ierr = get_extents_static(aoint_tensor, aoint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: aoint tensor corrupted'

         ierr = get_nmodes(moint_tensor,tens_rank)
         ierr = get_extents_static(moint_tensor, moint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: moint tensor corrupted'

 
         ! Make tensor to hold half-transformed integrals
         a2int_dims(1:2) = moint_dims(1:2)
         a2int_dims(3:4) = aoint_dims(3:4)
         ierr=construct_tensor_t(a2int_tensor,C8,a2int_dims,init_val=ZERO)

         ! Transform the first two indices
         call do_ht1 (aoint_tensor,a2int_tensor,aomo_spaces(1:2))

         ! Transform the last two indices, we need to indicate that we have the full AO ranges for these
         call do_ht2 (a2int_tensor,moint_tensor,aomo_spaces(3:4))

!        Destruct half-transformed integrals:
         ierr=destroy_tensor_t(a2int_tensor); if(ierr.ne.ETI_SUCCESS) then; ierr=35; return; endif

        end subroutine aomo_transform

        subroutine do_ht1 (aoint_tensor,a2int_tensor,aomo_spaces)         

!        Perform first half-transformation (indices 1 and 2)

         implicit none

         type(ao2mo_info), intent(in)   :: aomo_spaces(2)
         type(tensor_t), intent(inout)  :: aoint_tensor, a2int_tensor

         type(tensor_t):: mocoefa_tensor(2),mocoefb_tensor(2)
         type(tensor_t):: a1int_tensor
         integer(INTD) :: aoint_dims(1:4),mocoef_dims(1:2),a1int_dims(1:4)
         integer(INTD) :: a2int_dims(1:4)

         integer(INTD) :: ierr, tens_rank
         integer       :: i, nao 
         type(C_PTR)   :: body_p

!        Get tensor dimension of input and output tensors
         ierr = get_nmodes(aoint_tensor,tens_rank)
         ierr = get_extents_static(aoint_tensor, aoint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: aoint tensor corrupted'

         ierr = get_nmodes(a2int_tensor,tens_rank)
         ierr = get_extents_static(a2int_tensor, a2int_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: a2int tensor corrupted'

!        Initialize tensors with MO coefficients
         nao = aoint_dims(1) ! Note: We assume here that all AOs are active
         do i = 1, 2
            mocoef_dims(1) = aoint_dims(i)
            if (mocoef_dims(1) .ne. aomo_spaces(i)%nao) stop 'program error: inconsistent dimensions'
            mocoef_dims(2) = a2int_dims(i)
            if (mocoef_dims(2) .ne. aomo_spaces(i)%nmo) stop 'program error: inconsistent dimensions'
            ierr=construct_tensor_t(mocoefa_tensor(i),C8,mocoef_dims,init_val=ZERO)
            ierr=construct_tensor_t(mocoefb_tensor(i),C8,mocoef_dims,init_val=ZERO)
            call get_mocoef_tensor (1,nao, aomo_spaces(i)%nmo, aomo_spaces(i)%mo_list, mocoefa_tensor(i), mocoefb_tensor(i))
          end do

!        Set tensor dimensions for the index transformation (simple brute-force algorithm)
         a1int_dims(1)   = a2int_dims(1)
         a1int_dims(2:4) = aoint_dims(2:4)

!        Do 1st and 2nd quarter-transformation
         ierr=construct_tensor_t(a1int_tensor,C8,a1int_dims,init_val=ZERO)
         ! call print_date('constructed a1 and a2')
         call first_quarter_transform(mocoefa_tensor(1),aoint_tensor,a1int_tensor,aoint_dims,a1int_dims(1))
         ierr=t_einsum("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefa_tensor(2))
         if(ierr.ne.ETI_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=21; return; endif
         call first_quarter_transform(mocoefb_tensor(1),aoint_tensor,a1int_tensor,aoint_dims,a1int_dims(1))
         ierr=t_einsum("D(i,j,r,s)+=L(i,q,r,s)*R(q,j)",a2int_tensor,a1int_tensor,mocoefb_tensor(2))
         if(ierr.ne.ETI_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=22; return; endif
         ! call print_date('first ht done')

!        Destruct a1 integrals:
         ierr=destroy_tensor_t(a1int_tensor); if(ierr.ne.ETI_SUCCESS) then; ierr=25; return; endif

!        call print_date('leaving ao2mo_eti')

        end subroutine do_ht1

        subroutine first_quarter_transform (mocoef_tensor,l_tensor,d_tensor,l_dims,nmo)

!        D(i,q,r,s) = conj(C(p,i)) * L(p,q,r,s): with the last three indices taken together
!        this is D = C^H L, done as one ZGEMM on the (host, column-major) tensor bodies.
!        The executors do this contraction with conjugation much slower and without threads.

         implicit none

         type(tensor_t), intent(inout) :: mocoef_tensor, l_tensor, d_tensor
         integer(INTD), intent(in)     :: l_dims(4), nmo

         type(C_PTR)          :: pc, pl, pd
         complex(8), pointer  :: c(:,:), l(:,:), d(:,:)
         integer              :: nao, ncol
         integer(INTD)        :: ierr

         nao  = l_dims(1)
         ncol = l_dims(2)*l_dims(3)*l_dims(4)
         ierr = direct_data_pointer(mocoef_tensor,pc)
         ierr = direct_data_pointer(l_tensor,pl)
         ierr = direct_data_pointer(d_tensor,pd)
         call c_f_pointer(pc,c,(/nao,int(nmo)/))
         call c_f_pointer(pl,l,(/nao,ncol/))
         call c_f_pointer(pd,d,(/int(nmo),ncol/))
         call zgemm('C','N',int(nmo),ncol,nao,(1.0d0,0.0d0),c,nao,l,nao,(0.0d0,0.0d0),d,int(nmo))

        end subroutine first_quarter_transform

        subroutine do_ht2 (a2int_tensor,moint_tensor,aomo_spaces)

         ! Perform second half-transformation (indices 3 and 4)

         implicit none

         type(ao2mo_info), intent(in)   :: aomo_spaces(2)
         type(tensor_t), intent(inout)  :: a2int_tensor, moint_tensor

         type(tensor_t):: mocoefa_tensor(2),mocoefb_tensor(2)
         type(tensor_t):: a3int_tensor

         integer(INTD) :: a2int_dims(1:4),mocoef_dims(1:2),a3int_dims(1:4)
         integer(INTD) :: moint_dims(1:4)

 
         integer(INTD)     :: ierr, nao, tens_rank

         ! Check integrity of tensors
         ierr = get_nmodes(a2int_tensor,tens_rank)
         ierr = get_extents_static(a2int_tensor,a2int_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('input tensor corrupted in do_ht2')

         ierr = get_nmodes(moint_tensor,tens_rank)
         ierr = get_extents_static(moint_tensor,moint_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) call quit('output tensor corrupted in do_ht2')

!        Initialize tensors with MO coefficients for index 3
         nao = a2int_dims(3)
         mocoef_dims(1) = a2int_dims(3)
         mocoef_dims(2) = moint_dims(3)
         ierr=construct_tensor_t(mocoefa_tensor(1),C8,mocoef_dims,init_val=ZERO)
         ierr=construct_tensor_t(mocoefb_tensor(1),C8,mocoef_dims,init_val=ZERO)
         call get_mocoef_tensor (aomo_spaces(1)%first_ao, aomo_spaces(1)%last_ao, &
            aomo_spaces(1)%nmo, aomo_spaces(1)%mo_list, mocoefa_tensor(1), mocoefb_tensor(1))

!        Initialize tensors with MO coefficients for index 4
         nao = a2int_dims(3)
         mocoef_dims(1) = a2int_dims(4)
         mocoef_dims(2) = moint_dims(4)
         ierr=construct_tensor_t(mocoefa_tensor(2),C8,mocoef_dims,init_val=ZERO)
         ierr=construct_tensor_t(mocoefb_tensor(2),C8,mocoef_dims,init_val=ZERO)
         call get_mocoef_tensor (aomo_spaces(2)%first_ao, aomo_spaces(2)%last_ao, &
            aomo_spaces(2)%nmo, aomo_spaces(2)%mo_list, mocoefa_tensor(2), mocoefb_tensor(2))

!        Do 3rd and 4th quarter-transformation
         a3int_dims(1:3) = moint_dims(1:3)
         a3int_dims(4)   = a2int_dims(4)
         ierr=construct_tensor_t(a3int_tensor,C8,a3int_dims,init_val=ZERO)
         ierr=t_einsum("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefa_tensor(1),a2int_tensor)
         if(ierr.ne.ETI_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=31; return; endif
         ierr=t_einsum("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefa_tensor(2))
         if(ierr.ne.ETI_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=41; return; endif
         ierr=initialize_tensor_t(a3int_tensor,ZERO)
         ierr=t_einsum("D(i,j,k,s)+=R+(r,k)*L(i,j,r,s)",a3int_tensor,mocoefb_tensor(1),a2int_tensor)
         if(ierr.ne.ETI_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=32; return; endif
         ierr=t_einsum("D(i,j,k,l)+=L(i,j,k,s)*R(s,l)",moint_tensor,a3int_tensor,mocoefb_tensor(2))
         if(ierr.ne.ETI_SUCCESS) then; write(*,'("Error ",i11)') ierr; ierr=42; return; endif
         ! call print_date('second ht done')

!        Destruct a3 integrals:
         ierr=destroy_tensor_t(a3int_tensor); if(ierr.ne.ETI_SUCCESS) then; ierr=45; return; endif

        end subroutine ! do_ht2

        subroutine mulliken_to_dirac_sort  (mulliken_tensor,dirac_tensor,antisymmetrize)

         implicit none

         type(tensor_t), intent(inout)  :: mulliken_tensor,dirac_tensor
         integer, intent(in)                :: antisymmetrize

         integer(INTD) :: mulliken_dims(1:4),dirac_dims(1:4)

         integer(INTD) :: ierr, tens_rank, mismatch

!        work-around: addition with permutation is not implemented
         type(tensor_t) :: one_tensor, minusone_tensor
         integer(C_INT)     :: scalar_dims(1)
         scalar_dims(1) = 1
         ierr=construct_tensor_t(one_tensor,C8,scalar_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(minusone_tensor,C8,scalar_dims(1:0),init_val=MINUS_ONE)
!        end work-round

!        get dimensions of the tensors and check consistency
         ierr = get_nmodes(mulliken_tensor,tens_rank)
         ierr = get_extents_static(mulliken_tensor,mulliken_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input tensor corrupted in mulliken_to_dirac_sort'

         ierr = get_nmodes(dirac_tensor,tens_rank)
         ierr = get_extents_static(dirac_tensor,dirac_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'output tensor corrupted in mulliken_to_dirac_sort'

         mismatch = 0
         if (mulliken_dims(1) .ne. dirac_dims(1)) mismatch = mismatch + 1
         if (mulliken_dims(2) .ne. dirac_dims(3)) mismatch = mismatch + 1
         if (mulliken_dims(3) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (mulliken_dims(4) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 12 .and. dirac_dims(1) .ne. dirac_dims(2)) mismatch = mismatch + 1
         if (antisymmetrize .eq. 34 .and. dirac_dims(3) .ne. dirac_dims(4)) mismatch = mismatch + 1
         if (mismatch.gt.0) stop 'dimension mismatch in mulliken_to_dirac_sort'

         ierr=initialize_tensor_t(dirac_tensor)
!shorter ierr=tensor_add("D(p,q,r,s)+=M(p,r,q,s)",dirac_tensor,mulliken_tensor)
         ierr=t_einsum("D(p,q,r,s)+=M(p,r,q,s)*C()",dirac_tensor,mulliken_tensor,one_tensor)

         if (antisymmetrize .eq. 12) then
!shorter     ierr=tensor_add("D(p,q,r,s)+=M(q,r,p,s)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
             ierr=t_einsum("D(p,q,r,s)+=M(q,r,p,s)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
         elseif (antisymmetrize .eq. 34) then
!shorter     ierr=tensor_add("D(p,q,r,s)+=M(p,s,q,r)",dirac_tensor,mulliken_tensor,scale=MINUS_ONE)
              ierr=t_einsum("D(p,q,r,s)+=M(p,s,q,r)*C()",dirac_tensor,mulliken_tensor,minusone_tensor)
         end if

!        work-around: addition with permutation is not implemented
         ierr=destroy_tensor_t(one_tensor)
         ierr=destroy_tensor_t(minusone_tensor)
!        end work-round

        end subroutine mulliken_to_dirac_sort

        subroutine ao2mo_ind(ao_tensor,mo_tensor,nmo,mo_list,ind)
!         Transform integral tensor from ao-basis to mo-basis

          use eti_common_routines

          implicit none

          type(tensor_t), intent(inout) :: ao_tensor
          type(tensor_t), intent(inout) :: mo_tensor
          integer(INTD), intent(in) :: nmo          ! the size of the mo basis for each of the 2 indices (mulliken ordering)
          integer, intent(in)       :: mo_list(:)   ! and their indices (2 lists consecutively)
          integer, intent(in)       :: ind

          type(tensor_t):: mocoef_alpha, mocoef_beta
          integer(INTD)     :: mocoef_dims(1:2)
          integer(INTD)     :: ao_dim2(1:2), mo_dim2(1:2)
          integer(INTD)     :: ao_dim3(1:3), mo_dim3(1:3)
          integer(INTD)     :: ao_dim4(1:4), mo_dim4(1:4)
          integer(INTD)     :: ierr, rank, i

!         Get rank of ao tensor
          ierr = get_nmodes(ao_tensor,rank)

!         for one electron integrals / matrices
          if (rank.eq.2) then
            ierr = get_nmodes(ao_tensor,rank)
            ierr = get_extents_static(ao_tensor,ao_dim2)
            if (ierr.ne.0 .or. rank.ne.2) stop 'error in ao2mo_ind: ao tensor corrupted'

!           check output tensor
            ierr = get_nmodes(mo_tensor,rank)
            ierr = get_extents_static(mo_tensor,mo_dim2)
            if (ierr.ne.0 .or. rank.ne.2) stop 'error in ao2mo_ind: mo tensor corrupted'

            if (mo_dim2(ind).ne.nmo) stop 'error in ao2mo_ind: mo tensor wrong nmo1'

!           Initialize tensors with MO coefficients
            mocoef_dims(1) = ao_dim2(ind)
            mocoef_dims(2) = mo_dim2(ind)
            ierr=construct_tensor_t(mocoef_alpha,C8,mocoef_dims,init_val=ZERO)
            ierr=construct_tensor_t(mocoef_beta,C8,mocoef_dims,init_val=ZERO)
            call get_mocoef_tensor (1,mocoef_dims(1),mocoef_dims(2), & 
                        mo_list(1:mocoef_dims(2)), mocoef_alpha, mocoef_beta)

!           transform the index
            if (ind.eq.1) then
              ierr=t_einsum("B(p,q)+=M+(i,p)*A(i,q)",mo_tensor,mocoef_alpha,ao_tensor)
              if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 2a1'
              ierr=t_einsum("B(p,q)+=M+(i,p)*A(i,q)",mo_tensor,mocoef_beta,ao_tensor)
              if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 2b1'
            else if (ind.eq.2) then
                ierr=t_einsum("C(p,q)+=B(p,i)*M(i,q)",mo_tensor,ao_tensor,mocoef_alpha)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 2a2'
                ierr=t_einsum("C(p,q)+=B(p,i)*M(i,q)",mo_tensor,ao_tensor,mocoef_beta)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 2b2'
            else
              stop 'error in ao2mo_ind: index not available'
            end if 

!           clean up
            ierr=destroy_tensor_t(mocoef_alpha)
            ierr=destroy_tensor_t(mocoef_beta)
            
!         cholesky vectors
          else if (rank.eq.3) then
!           ierr = tensor_dimensions(ao_tensor,rank,ao_dim3)
            ierr = get_nmodes(ao_tensor,rank)
            ierr = get_extents_static(ao_tensor,ao_dim3)
            if (ierr.ne.0 .or. rank.ne.3) stop 'error in ao2mo_ind: ao tensor corrupted'

!           check output tensor
!           ierr = tensor_dimensions(mo_tensor,rank,mo_dim3)
            ierr = get_nmodes(mo_tensor,rank)
            ierr = get_extents_static(mo_tensor,mo_dim3)
            if (ierr.ne.0 .or. rank.ne.3) stop 'error in ao2mo_ind: mo tensor corrupted'

            if (mo_dim3(3).ne.ao_dim3(3)) stop 'error in ao2mo_ind: mo tensor wrong chol. dim.'
            if (ind.eq.11 .or. ind.eq.12) then
                i=1
            else if (ind.eq.21 .or. ind.eq.22) then
                i=2
            else
              stop 'error in ao2mo_ind: index not available'
            end if
            if (mo_dim3(i).ne.nmo) stop 'error in ao2mo_ind: mo tensor wrong nmo'

!           Initialize tensors with MO coefficients
            mocoef_dims(1) = ao_dim3(i)
            mocoef_dims(2) = mo_dim3(i)
            ierr=construct_tensor_t(mocoef_alpha,C8,mocoef_dims,init_val=ZERO)
            ierr=construct_tensor_t(mocoef_beta,C8,mocoef_dims,init_val=ZERO)
            call get_mocoef_tensor (1,mocoef_dims(1),mocoef_dims(2), & 
                       mo_list(1:mocoef_dims(2)), mocoef_alpha, mocoef_beta)

!           transform index
            if (ind.eq.11) then
                ierr=t_einsum("A(p,s,u)+=M+(r,p)*C(r,s,u)",mo_tensor,mocoef_alpha,ao_tensor)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 3a1'
            else if (ind.eq.12) then
                ierr=t_einsum("A(p,s,u)+=M+(r,p)*C(r,s,u)",mo_tensor,mocoef_beta,ao_tensor)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 3b1'
            else if (ind.eq.21) then
                ierr=t_einsum("D(p,q,u)+=A(p,s,u)*M(s,q)",mo_tensor,ao_tensor,mocoef_alpha)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 3a2'
            else if (ind.eq.22) then
                ierr=t_einsum("D(p,q,u)+=A(p,s,u)*M(s,q)",mo_tensor,ao_tensor,mocoef_beta)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_ind: 3b2'
            end if

!           clean up
            ierr=destroy_tensor_t(mocoef_alpha)
            ierr=destroy_tensor_t(mocoef_beta)

          else 
            stop 'error in ao2mo_ind: wrong rank'
          end if

        end subroutine ao2mo_ind

        subroutine ao2mo_vec(ao_vec,mo_vec,m,nmo,mo_list,ind)
!         Transform integral tensor vector from ao-basis to mo-basis

          use eti_common_routines

          implicit none
         
          type(tensor_t), allocatable :: ao_vec(:)
          type(tensor_t), allocatable :: mo_vec(:)
          integer(INTD), intent(in)       :: nmo          ! the size of the mo basis for each of the 2 indices (mulliken ordering)
          integer, intent(in)             :: mo_list(:)   ! and their indices (2 lists consecutively)
          integer, intent(in)             :: ind, m

          type(tensor_t) :: mocoef_alpha, mocoef_beta
          integer(INTD)      :: mocoef_dims(1:2)
          integer(INTD)      :: ao_dim2(1:2), mo_dim2(1:2)
          integer(INTD)      :: ierr, rank
          integer            :: i,ii

!         Get info about ao tensor
          ierr = get_nmodes(ao_vec(1), rank)
          ierr = get_extents_static(ao_vec(1),ao_dim2)

          if (ierr.ne.0 .or. rank.ne.2) stop 'error in ao2mo_vec: ao tensor corrupted'

          if (ind.eq.11 .or. ind.eq.12) then
            ii=1
          else if (ind.eq.21 .or. ind.eq.22) then
            ii=2
          else
            stop 'error in ao2mo_ind: index not available'
          end if
          
          mo_dim2=ao_dim2
          mo_dim2(ii)=nmo

!         Initialize tensors with MO coefficients
          mocoef_dims(1) = ao_dim2(ii)
          mocoef_dims(2) = mo_dim2(ii)
          ierr=construct_tensor_t(mocoef_alpha,C8,mocoef_dims,init_val=ZERO)
          ierr=construct_tensor_t(mocoef_beta,C8,mocoef_dims,init_val=ZERO)
          call get_mocoef_tensor (1,mocoef_dims(1),mocoef_dims(2), & 
                        mo_list(1:mocoef_dims(2)), mocoef_alpha, mocoef_beta)

          do i=1,m
!           transform the index
            if (ind.eq.11) then
              ierr=t_einsum("B(p,q)+=M+(i,p)*A(i,q)",mo_vec(i),mocoef_alpha,ao_vec(i))
              if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_vec: 2a1'
            else if (ind.eq.12) then
              ierr=t_einsum("B(p,q)+=M+(i,p)*A(i,q)",mo_vec(i),mocoef_beta,ao_vec(i))
              if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_vec: 2b1'
            else if (ind.eq.21) then
                ierr=t_einsum("C(p,q)+=B(p,i)*M(i,q)",mo_vec(i),ao_vec(i),mocoef_alpha)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_vec: 2a2'
            else if (ind.eq.22) then
                ierr=t_einsum("C(p,q)+=B(p,i)*M(i,q)",mo_vec(i),ao_vec(i),mocoef_beta)
                if(ierr.ne.ETI_SUCCESS) stop 'error in ao2mo_vec: 2b2'
            end if 
          end do

!         clean up
          ierr=destroy_tensor_t(mocoef_alpha)
          ierr=destroy_tensor_t(mocoef_beta)

        end subroutine ao2mo_vec

end module eti_ao_to_mo
