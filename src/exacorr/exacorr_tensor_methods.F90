module exacorr_tensor_methods

!This module contains routines to modify (initialize, scale with denominators, etc.) that can be registered with exatensor

        use exatensor
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

        type, extends(tens_method_uni_t), public:: tens_init_test_t
         complex(8), private:: init_val=(0d0,0d0)
         contains
          procedure, private:: TensInitTestCtorReal
          procedure, private:: TensInitTestCtorComplex
          generic, public:: tens_init_test_ctor=>TensInitTestCtorReal,TensInitTestCtorComplex
          procedure, public:: apply=>TensInitTestApply
        end type tens_init_test_t

        !Initializer for the 2e AO integral tensor:
        type, extends(tens_method_uni_t), public:: compute_2e_ao_tensor_t
         contains
          procedure, public:: apply=>compute_2e_ao_tensor
        end type compute_2e_ao_tensor_t

        !Initializer for the mo coefficients tensor:
        type, extends(tens_method_uni_t), public:: init_mocoef_tensor_t
         complex(8), allocatable, private:: mocoef_all(:,:) ! stores mo-coefficients for both indices, both spins
         contains
          procedure, public:: init_mocoef_ctor=>init_mocoef
          procedure, public:: apply=>init_mocoef_tensor
        end type init_mocoef_tensor_t

        !Initializer for the density matrix tensor:
        type, extends(tens_method_uni_t), public:: init_dm_tensor_t
         complex(8), allocatable, private:: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)
         contains
          procedure, public:: init_dm_ctor=>init_dm
          procedure, public:: apply=>init_dm_tensor
        end type init_dm_tensor_t

        !Initializer for the half-transformed integrals:
        type, extends(tens_method_uni_t), public:: init_ht_tensor_t
         complex(8), allocatable, private:: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)
         contains
          procedure, public:: init_ht_ctor=>init_ht
          procedure, public:: apply=>init_ht_tensor
        end type init_ht_tensor_t

        !Scale by denominators:
        type, extends(tens_method_uni_t), public:: scale_with_denom_t
         real(8), allocatable, private:: eps_occ(:),eps_vir(:)
         contains
          procedure, public:: scale_with_denom_ctor=>init_epsilons
          procedure, public:: apply=>scale_with_denom
        end type scale_with_denom_t

        !set diagonal
        type, extends(tens_method_uni_t), public:: set_diagonal_t
         real(8), private :: constant
         contains
          procedure, public:: init_const_ctor=>set_const
          procedure, public:: apply=>diag_set
        end type set_diagonal_t

       contains

       subroutine TensInitTestCtorReal(this,val)
         implicit none
         class(tens_init_test_t), intent(out):: this
         real(8), intent(in):: val

         this%init_val=dcmplx(val,0d0)
        end subroutine TensInitTestCtorReal

        subroutine TensInitTestCtorComplex(this,val)
         implicit none
         class(tens_init_test_t), intent(out):: this
         complex(8), intent(in):: val

         this%init_val=val
        end subroutine TensInitTestCtorComplex

        function TensInitTestApply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(tens_init_test_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         integer(INTL):: vol
         type(tens_dense_t):: tens
         real(4), pointer:: arr_r4(:)
         real(8), pointer:: arr_r8(:)
         complex(4), pointer:: arr_c4(:)
         complex(8), pointer:: arr_c8(:)

         tens=tensor%get_dense_adapter(ierr)
         if(ierr.eq.EXA_SUCCESS) then
          vol=tens_dense_volume(tens)
          if(vol.gt.0) then
           select case(tens%data_kind)
           case(EXA_DATA_KIND_R4)
            call c_f_pointer(tens%body_ptr,arr_r4,(/vol/))
            arr_r4(:)=real(this%init_val,4)
           case(EXA_DATA_KIND_R8)
            call c_f_pointer(tens%body_ptr,arr_r8,(/vol/))
            arr_r8(:)=real(this%init_val,8)
           case(EXA_DATA_KIND_C4)
            call c_f_pointer(tens%body_ptr,arr_c4,(/vol/))
            arr_c4(:)=cmplx(real(this%init_val,4),real(imag(this%init_val),4),4)
           case(EXA_DATA_KIND_C8)
            call c_f_pointer(tens%body_ptr,arr_c8,(/vol/))
            arr_c8(:)=this%init_val
           case default
            ierr=2
           end select
          else
           ierr=1
          endif
         endif
        end function TensInitTestApply

        function compute_2e_ao_tensor(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr,n
         class(compute_2e_ao_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.4) call quit('Error: dimension not 4 in compute_2e_ao_tensor')
         ierr=tens_init_ao_integrals(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n))
         return
        end function compute_2e_ao_tensor

        integer function tens_init_ao_integrals(tens_body, data_kind, tens_dims, tens_bases) result (ierr)
        !Initialization routine that is used to fill a tensor with AO integrals

         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri
         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(C_PTR), value:: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(4), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)

         integer:: nao(4)
         complex(8),pointer :: ao_tens(:,:,:,:)
         integer :: first_ao(1:4)

         integer :: nshells(4)
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         !$OMP threadprivate(gout)
         real(8), parameter :: fijkl = 1.0

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

!        Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
         call interest_initialize(.false.)

         ! map the arguments to names used in the context of integral evaluation
         first_ao(1:4) = tens_bases(1:4)
         nao(1:4) = tens_dims(1:4)
         call c_f_pointer(tens_body,ao_tens,nao(1:4)) !map C pointer to 4-dimensional array

!        Allocate and get shell information for all shells that we will treat
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))

!        Intranode parallelization should be possible in the loops below because Interest is thread safe

         loff = 0
         do lsh = 1, nshells(4)
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            koff = 0
            do ksh = 1, nshells(3)
               lk   =  gto3(ksh)%orb_momentum
               ek   =  gto3(ksh)%exponent(1)
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               ck   =  gto3(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               joff = 0
               do jsh = 1, nshells(2)
                  lj   =  gto2(jsh)%orb_momentum
                  ej   =  gto2(jsh)%exponent(1)
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  cj   =  gto2(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  ioff = 0
                  do ish = 1, nshells(1)
                     li   =  gto1(ish)%orb_momentum
                     ei   =  gto1(ish)%exponent(1)
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ci   =  gto1(ish)%coefficient(1)
                     ni   =  ncartesian(li)
!                    output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                     ijkl = 0
                     do l = 1, nl
                        do k = 1, nk
                           do j = 1, nj
                              do i = 1, ni
                                 ijkl = ijkl + 1
                                 ao_tens(i+ioff,j+joff,k+koff,l+loff) = dcmplx(gout(ijkl),0.D0)
                              end do
                           end do
                        end do
                     end do
                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)

        end function tens_init_ao_integrals


        subroutine init_epsilons(this,mo_occ,nocc,mo_vir,nvir)

         use dirac_MO
         use exacorr_global

         implicit none
         integer, intent(in) :: nocc, nvir
         integer, intent(in) :: mo_occ(:),mo_vir(:)
         class(scale_with_denom_t), intent(out):: this

         type(MO) :: qmo
         integer  :: k

         ! we are assuming closed shell (Kramers-restricted) cases here
         call get_mo_coefficients (qmo,mo_occ,nocc)
         allocate(this%eps_occ(2*nocc))
         do k = 1, nocc
            this%eps_occ(k)      = qmo%energy(k)
            this%eps_occ(k+nocc) = qmo%energy(k)
         end do
         call delete_qmo (qmo)

         call get_mo_coefficients (qmo,mo_vir,nvir)
         allocate(this%eps_vir(2*nvir))
         do k = 1, nvir
            this%eps_vir(k)      = qmo%energy(k)
            this%eps_vir(k+nvir) = qmo%energy(k)
         end do
         call delete_qmo (qmo)

        end subroutine init_epsilons

        subroutine init_mocoef(this,mo_list,nmo,spin)

         use dirac_MO
         use exacorr_global

         implicit none
         integer, intent(in) :: nmo          ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         integer, intent(in) :: spin
         class(init_mocoef_tensor_t), intent(out):: this

         integer  :: nao,j,k
         type(MO) :: qmo

!        Copy the coefficients into the alpha or beta array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         nao = get_nao()
         allocate (this%mocoef_all(nao,2*nmo))
         call get_mo_coefficients (qmo,mo_list,nmo)

         if (spin == 1) then
            do k = 1, nmo
               do j = 1, nao
                this%mocoef_all(j,k)     = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
                this%mocoef_all(j,k+nmo) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
               end do
            end do
         elseif (spin == 2) then
            do k = 1, nmo
               do j = 1, nao
                this%mocoef_all(j,k)     = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
                this%mocoef_all(j,k+nmo) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
               end do
            end do
         else
            call quit ('init_mocoef called with unknown spin value')
         end if
         call delete_qmo (qmo)

        end subroutine init_mocoef

        subroutine init_dm(this,mo_list,nmo)

         use dirac_MO
         use exacorr_global

         implicit none
         integer, intent(in) :: nmo(2)       ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         class(init_dm_tensor_t), intent(out):: this

         integer  :: nao,j,k
         type(MO) :: qmo

!        Copy the coefficients into the alpha or beta array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         nao = get_nao()

         allocate (this%mocoef_all_1(nao,2*nmo(1),2))
         call get_mo_coefficients (qmo,mo_list(1:nmo(1)),nmo(1))
         do k = 1, nmo(1)
            do j = 1, nao
             this%mocoef_all_1(j,k,1)     = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
             this%mocoef_all_1(j,k,2)     = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_1(j,k+nmo(1),1) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_1(j,k+nmo(1),2) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
            end do
         end do
         call delete_qmo (qmo)

         allocate (this%mocoef_all_2(nao,2*nmo(2),2))
         call get_mo_coefficients (qmo,mo_list(nmo(1)+1:nmo(1)+nmo(2)),nmo(2))
         do k = 1, nmo(2)
            do j = 1, nao
             this%mocoef_all_2(j,k,1)     = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
             this%mocoef_all_2(j,k,2)     = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_2(j,k+nmo(2),1) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_2(j,k+nmo(2),2) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
            end do
         end do
         call delete_qmo (qmo)

        end subroutine init_dm

        subroutine init_ht(this,mo_list,nmo)

         use dirac_MO
         use exacorr_global

         implicit none
         integer, intent(in) :: nmo(2)       ! the length of the mo basis
         integer, intent(in) :: mo_list(:)   ! and their indices
         class(init_ht_tensor_t), intent(out):: this

         integer  :: nao,j,k
         type(MO) :: qmo

!        Copy the coefficients into the alpha or beta array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.
!        The coefficients for the barred spinors (that start after the nmo(i) unbarred spinors) are generated
!        using Kramers' symmetry.

         nao = get_nao()

         allocate (this%mocoef_all_1(nao,2*nmo(1),2))
         call get_mo_coefficients (qmo,mo_list(1:nmo(1)),nmo(1))
         do k = 1, nmo(1)
            do j = 1, nao
             this%mocoef_all_1(j,k,1)     = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
             this%mocoef_all_1(j,k,2)     = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_1(j,k+nmo,1) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_1(j,k+nmo,2) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
            end do
         end do
         call delete_qmo (qmo)

         allocate (this%mocoef_all_2(nao,2*nmo(2),2))
         call get_mo_coefficients (qmo,mo_list(nmo(1)+1:nmo(1)+nmo(2)),nmo(2))
         do k = 1, nmo(2)
            do j = 1, nao
             this%mocoef_all_2(j,k,1)     = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))
             this%mocoef_all_2(j,k,2)     = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_2(j,k+nmo,1) = dcmplx(qmo%coeff(j,k,3),qmo%coeff(j,k,4))
             this%mocoef_all_2(j,k+nmo,2) = dcmplx(qmo%coeff(j,k,1),-qmo%coeff(j,k,2))
            end do
         end do
         call delete_qmo (qmo)

        end subroutine init_ht

        function scale_with_denom(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr,n
         class(scale_with_denom_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         ierr=denominate(tens%num_dims,tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                         this%eps_occ,this%eps_vir)

        end function scale_with_denom

        integer function denominate(num_dims,tens_body,data_kind,tens_dims,tens_bases, &
                         eps_occ,eps_vir) result (ierr)

        !Scale with denominators for 2-d or 4-d tensor for regular CC and for lambda

         use, intrinsic:: ISO_C_BINDING

         implicit none

         integer,intent(in):: num_dims            !in: number of dimensions of the tensor
         type(C_PTR), value:: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(4), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)
         real(8), intent(in):: eps_occ(:),eps_vir(:)

         integer :: i,j,a,b
         real(8) :: denominator
         complex(8),pointer :: t1_tens(:,:),t2_tens(:,:,:,:)

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         select case (num_dims)
         case (2)
            call c_f_pointer(tens_body,t1_tens,tens_dims(1:2)) !map C pointer to 2-dimensional array
!           if (tens_dims(1)*2.eq.size(eps_vir)) then 
               do i = 1, tens_dims(2)
                  do a = 1, tens_dims(1)
                     denominator = eps_occ(i+tens_bases(2)-1) - eps_vir(a+tens_bases(1)-1)
                     t1_tens(a,i) = t1_tens(a,i) / denominator
                  end do
               end do
!           else
!              do i = 1, tens_dims(1)
!                 do a = 1, tens_dims(2)
!                    denominator = eps_occ(i+tens_bases(1)-1) - eps_vir(a+tens_bases(2)-1)
!                    t1_tens(i,a) = t1_tens(i,a) / denominator
!                 end do
!              end do
!           end if
         case (4)
            call c_f_pointer(tens_body,t2_tens,tens_dims(1:4)) !map C pointer to 2-dimensional array
!           if (tens_dims(1)*2.eq.size(eps_vir)) then
               do j = 1, tens_dims(4)
                  do i = 1, tens_dims(3)
                     do b = 1, tens_dims(2)
                        do a = 1, tens_dims(1)
                           denominator = eps_occ(i+tens_bases(3)-1) + eps_occ(j+tens_bases(4)-1) &
                                       - eps_vir(a+tens_bases(1)-1) - eps_vir(b+tens_bases(2)-1)
                           t2_tens(a,b,i,j) = t2_tens(a,b,i,j) / denominator
                        end do
                     end do
                  end do
               end do
!           else
!              do j = 1, tens_dims(2)
!                 do i = 1, tens_dims(1)
!                    do b = 1, tens_dims(4)
!                       do a = 1, tens_dims(3)
!                          denominator = eps_occ(i+tens_bases(1)-1) + eps_occ(j+tens_bases(2)-1) &
!                                      - eps_vir(a+tens_bases(3)-1) - eps_vir(b+tens_bases(4)-1)
!                          t2_tens(i,j,a,b) = t2_tens(i,j,a,b) / denominator
!                       end do
!                    end do
!                 end do
!              end do
!           end if
         case default
            call quit ('Error: wrong dimension of tensor in denominate')
         end select

        end function denominate

        function init_mocoef_tensor(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr,n
         class(init_mocoef_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.2) call quit('Error: dimension not 2 in init_mocoef_tensor')
         ierr=tens_init_mocoef(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                               this%mocoef_all)

        end function init_mocoef_tensor

        integer function tens_init_mocoef(tens_body,data_kind,tens_dims,tens_bases, &
                         mocoef_all) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(C_PTR), value:: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(4), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)
         complex(8), intent(in):: mocoef_all(:,:)

         integer :: nao, nmo
         integer :: first_ao, last_ao, first_mo, last_mo
         complex(8),pointer :: mocoef_tens(:,:)

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         nao = tens_dims(1)
         nmo = tens_dims(2)
         call c_f_pointer(tens_body,mocoef_tens,(/nao,nmo/)) !map C pointer to 2-dimensional array

         first_ao = tens_bases(1)
         first_mo = tens_bases(2)
         last_ao  = first_ao + nao - 1
         last_mo  = first_mo + nmo - 1

         mocoef_tens(1:nao,1:nmo) = mocoef_all(first_ao:last_ao,first_mo:last_mo)

        end function tens_init_mocoef

        function init_dm_tensor(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr,n
         class(init_dm_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.4) call quit('Error: dimension not 4 in init_dm_tensor')
         ierr=tens_init_dm(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                               this%mocoef_all_1,this%mocoef_all_2)

        end function init_dm_tensor

        integer function tens_init_dm(tens_body,data_kind,tens_dims,tens_bases, &
                         mocoef_all_1,mocoef_all_2) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(C_PTR), value:: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(4), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)
         complex(8), intent(in):: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)

         integer :: p,q,i,j,ia,ja,pa,qa,nao(2),nmo(2)
         complex(8),pointer :: dm_tens(:,:,:,:)

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         nao = tens_dims(1:2)
         nmo = tens_dims(3:4)
         call c_f_pointer(tens_body,dm_tens,(/nao(1),nao(2),nmo(1),nmo(2)/)) !map C pointer to 4-dimensional array

         ! carry out the spin integration and get the tensor
         do j = 1, nmo(2)
            ja = tens_bases(4) + j - 1 
            do i = 1, nmo(1)
               ia = tens_bases(3) + i - 1 
               do q = 1, nao(2)
                  qa = tens_bases(2) + q - 1 
                  do p = 1, nao(1)
                     pa = tens_bases(1) + p - 1 
                     dm_tens(p,q,i,j) = &
                       dconjg(mocoef_all_1(pa,ia,1))*mocoef_all_2(qa,ja,1) & ! alpha-alpha
                     + dconjg(mocoef_all_1(pa,ia,2))*mocoef_all_2(qa,ja,2)   ! beta-beta
                  end do
               end do
            end do
         end do

        end function tens_init_dm

        function init_ht_tensor(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr,n
         class(init_ht_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         n=tens%num_dims
         if (n.ne.4) call quit('Error: dimension not 4 in init_dm_tensor')
         ierr=tens_init_ht(tens%body_ptr,tens%data_kind,tens%dims(1:n),tens%bases(1:n), &
                               this%mocoef_all_1,this%mocoef_all_2)

        end function init_ht_tensor

        integer function tens_init_ht(tens_body,data_kind,tens_dims,tens_bases, &
                         mocoef_all_1,mocoef_all_2) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri
         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(C_PTR), value:: tens_body           !in: C pointer to the tensor body (tensor elements)
         integer(4), value:: data_kind            !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in):: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in):: tens_bases(1:*) !in: base offsets for each tensor dimension)
         complex(8), intent(in):: mocoef_all_1(:,:,:),mocoef_all_2(:,:,:)

         complex(8),pointer :: ht_tens(:,:,:,:)

         integer :: first_ao(1:4)
         integer :: nshells(4)
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: ioff,joff,koff,loff,kmo,lmo,ka,la,pa,qa,ra,sa,nao(4),nmo(2),ijkl
         integer :: i, j, k, l, ish, jsh, ksh, lsh
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         !$OMP threadprivate(gout)
         complex(8) :: dP
         real(8), parameter :: fijkl = 1.0

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

!        Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
         call interest_initialize(.false.)

         ! map the arguments to names used in the context of integral evaluation
         first_ao(1:2) = tens_bases(1:2)
         first_ao(3:4) = 1
         nao(1:2)      = tens_dims(1:2)
         nmo(1:2)      = tens_dims(3:4)
         nao(3:4)      = get_nao()
         call c_f_pointer(tens_body,ht_tens,(/nao(1),nao(2),nmo(1),nmo(2)/)) !map C pointer to 4-dimensional array
         ht_tens = ZERO

!        Allocate and get shell information for all shells that we will treat
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))

         loff = 0
         do lsh = 1, nshells(4) ! loop over all ao's, this index is transformed to mo
            ll   =  gto4(lsh)%orb_momentum
            el   =  gto4(lsh)%exponent(1)
            xl   =  gto4(lsh)%coord(1)
            yl   =  gto4(lsh)%coord(2)
            zl   =  gto4(lsh)%coord(3)
            cl   =  gto4(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            koff = 0
            do ksh = 1, nshells(3) ! loop over all ao's, this index is transformed to mo
               lk   =  gto3(ksh)%orb_momentum
               ek   =  gto3(ksh)%exponent(1)
               xk   =  gto3(ksh)%coord(1)
               yk   =  gto3(ksh)%coord(2)
               zk   =  gto3(ksh)%coord(3)
               ck   =  gto3(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               joff = 0
               do jsh = 1, nshells(2) ! loop over the ao's active in the HT integral tensor block
                  lj   =  gto2(jsh)%orb_momentum
                  ej   =  gto2(jsh)%exponent(1)
                  xj   =  gto2(jsh)%coord(1)
                  yj   =  gto2(jsh)%coord(2)
                  zj   =  gto2(jsh)%coord(3)
                  cj   =  gto2(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  ioff = 0
                  do ish = 1, nshells(1) ! loop over the ao's active in the HT integral tensor block
                     li   =  gto1(ish)%orb_momentum
                     ei   =  gto1(ish)%exponent(1)
                     xi   =  gto1(ish)%coord(1)
                     yi   =  gto1(ish)%coord(2)
                     zi   =  gto1(ish)%coord(3)
                     ci   =  gto1(ish)%coefficient(1)
                     ni   =  ncartesian(li)
                     ! output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )

                     do lmo = 1, nmo(2) ! loop over the mo's active in the HT integral tensor block
                        la = tens_bases(4) + lmo - 1 
                        do kmo = 1, nmo(1) ! loop over the mo's active in the HT integral tensor block
                           ka = tens_bases(3) + kmo - 1 
                           ijkl = 0
                           do l = 1, nl
                              sa = loff + l
                              do k = 1, nk
                                 ra = koff + k
                                 dP = dconjg(mocoef_all_1(ra,ka,1))*mocoef_all_2(sa,la,1) & ! alpha-alpha
                                    + dconjg(mocoef_all_1(ra,ka,2))*mocoef_all_2(sa,la,2)   ! beta-beta
                                 do j = 1, nj
                                    qa = joff + j
                                    do i = 1, ni
                                       pa = ioff + i
                                       ijkl = ijkl + 1
                                       ht_tens(pa,qa,kmo,lmo) = ht_tens(pa,qa,kmo,lmo) &
                                                              + dP * gout(ijkl)
                                    end do ! i
                                 end do ! j
                              end do ! k
                           end do ! l
                        end do ! kmo
                     end do ! lmo
                     ioff = ioff + ni
                  end do
                  joff = joff + nj
               end do
               koff = koff + nk
            end do
            loff = loff + nl
         end do

         deallocate(gto1)
         deallocate(gto2)
         deallocate(gto3)
         deallocate(gto4)

        end function tens_init_ht

        !subroutine init_diag(this,diag)
        ! use dirac_MO

        ! class(set_diagonal_t), intent(inout)   :: this
        ! real(8), intent(in) :: diag(:)
        ! integer :: ddim, k


        ! ddim=SIZE(diag)

        ! allocate(this%diag_vector(ddim))

        ! do k = 1, ddim
        !     this%diag_vector(k)=diag(k)
        ! end do

        !end subroutine init_diag

        subroutine set_const(this, const)

          class(set_diagonal_t), intent(out)   :: this
          real(8), intent(in) :: const

          this%constant=const

        end subroutine set_const

        function diag_set(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr,rank
         class(set_diagonal_t), intent(in)    :: this
         class(tens_rcrsv_t), intent(inout)   :: tensor
         complex(8), intent(inout), optional  :: scalar
         type(tens_dense_t):: tens

         ierr=0
         tens=tensor%get_dense_adapter(ierr)
         if (ierr.ne.EXA_SUCCESS) call quit('Error: getting dense_adapter')
         rank=tens%num_dims
         if (rank.ne.2) call quit('Error: dimension not 2 in init_dm_tensor')

         ierr=make_diagonal(tens%body_ptr,tens%data_kind,tens%dims(1:rank),tens%bases(1:rank), this%constant)

        end function diag_set

        function make_diagonal(tens_body,data_kind,tens_dims,tens_bases, constant) result (ierr)
         ! routine to fill the diagonal of a matrix

         use, intrinsic:: ISO_C_BINDING

         implicit none

         integer(INTD)          :: ierr
         type(C_PTR), value     :: tens_body       !in: C pointer to the tensor body (tensor elements)
         integer(4), value      :: data_kind       !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(8), intent(in) :: tens_dims(1:*)  !in: tensor dims (extent of each tensor dimension)
         integer(8), intent(in) :: tens_bases(1:*) !in: base offsets for each tensor dimension)
         real(8), intent(in)    :: constant


         integer            :: lo, la, ro, ra
         complex(8),pointer :: h_tens(:,:)
         complex(8)         :: const


          const=dcmplx(constant,0d0)

          ierr=0

          call c_f_pointer(tens_body,h_tens,tens_dims(1:2)) !map C pointer to 2-dimensional array

          !print *, "dim1=", tens_dims(1)
          !print *, "dim2=", tens_dims(2)

          do lo = 1, tens_dims(1) ! loop over lines
             la = tens_bases(1) + lo - 1
                  do ro = 1, tens_dims(2) ! loop over rows
                     ra = tens_bases(2) + ro - 1
                     if (la.eq.ra) h_tens(lo, ro) = const
                  end do
          end do

        end function make_diagonal

end module exacorr_tensor_methods
