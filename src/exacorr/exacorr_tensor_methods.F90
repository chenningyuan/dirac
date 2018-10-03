module exacorr_tensor_methods

!This module contains routines to modify (initialize, scale with denominators, etc.) that can be registered with exatensor

        use exatensor
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private

        type, extends(tens_method_uni_t), public:: tens_init_test_t
         complex(8), private:: init_val=(0d0,0d0)
         contains
          procedure, private:: TensInitTestCtorReal
          procedure, private:: TensInitTestCtorComplex
          generic, public:: tens_init_test_ctor=>TensInitTestCtorReal,TensInitTestCtorComplex
          procedure, public:: apply=>TensInitTestApply
        end type tens_init_test_t

        !Block-wise initializer for the 2e AO integral tensor:
        type, extends(tens_method_uni_t), public:: compute_2e_ao_tensor_t
         contains
          procedure, public:: apply=>compute_2e_ao_tensor
        end type compute_2e_ao_tensor_t

        !Initializer for the mo coefficients tensor:
        type, extends(tens_method_uni_t), public:: initialize_mocoef_tensor_t
         contains
          procedure, public:: apply=>initialize_mocoef_tensor
        end type initialize_mocoef_tensor_t

       contains

       subroutine TensInitTestCtorReal(this,val)
         implicit none
         class(tens_init_test_t), intent(out):: this
         real(8), intent(in):: val

         this%init_val=cmplx(val,0d0)
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
         integer(INTD):: ierr
         class(compute_2e_ao_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         integer(INTL):: tens_shape(1:MAX_TENSOR_RANK)
         integer(INTL):: tens_signa(1:MAX_TENSOR_RANK)
         integer(INTD):: data_kind,n
         type(C_PTR):: tens_body_p

         ierr=0
         data_kind=tensor%get_data_type(ierr)
         if(ierr.ne.0) then; ierr=-6; return; endif
         call tensor%get_dims(tens_shape,n,ierr)
         if(ierr.ne.0) then; ierr=-5; return; endif
         if(n.eq.4) then !order-4 tensor
          call tensor%get_bases(tens_signa,n,ierr)
          if(ierr.ne.0) then; ierr=-4; return; endif
          tens_body_p=tensor%get_body_ptr(ierr)
          if(ierr.ne.0) then; ierr=-3; return; endif
          ierr=tens_init_ao_integrals(tens_body_p,data_kind,tens_shape(1:n),tens_signa(1:n))
          if(ierr.ne.0) then; ierr=-2; return; endif
         else
          ierr=-1
         endif
         return
        end function compute_2e_ao_tensor

        integer function tens_init_ao_integrals(tens_body, data_kind, tens_shape, tens_base) result (ierr)
        !Initialization routine that is used to fill a tensor with AO integrals

         use exacorr_datatypes
         use exacorr_global
         use module_interest_eri
         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(C_PTR), value:: tens_body                     !in: C pointer to the tensor body (tensor elements)
         integer(C_INT), value:: data_kind                  !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(C_LONG_LONG), intent(in):: tens_shape(1:*) !in: tensor shape (extent of each tensor dimension)
         integer(C_LONG_LONG), intent(in):: tens_base(1:*)  !in: base offsets for each tensor dimension)

         integer:: nao(4)
         complex(8),pointer :: ao_tens(:,:,:,:)
         integer :: first_ao(1:4)

         integer :: nshells(4)
         type(basis_func_info_t), allocatable :: gto1(:), gto2(:), gto3(:), gto4(:)  ! arrays with basis function information

         integer :: i, j, k, l, ish, jsh, ksh, lsh, ioff, joff, koff, loff, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         real(8), parameter :: fijkl = 1.0

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         ! map the arguments to names used in the context of integral evaluation
         first_ao(1:4) = tens_base(1:4)
         nao(1:4) = tens_shape(1:4)
         call c_f_pointer(tens_body,ao_tens,nao(1:4)) !map C pointer to 4-dimensional array

!        Allocate and get shell information for all shells that we will treat
         call get_gtos(first_ao(1),nao(1),gto1,nshells(1))
         call get_gtos(first_ao(2),nao(2),gto2,nshells(2))
         call get_gtos(first_ao(3),nao(3),gto3,nshells(3))
         call get_gtos(first_ao(4),nao(4),gto4,nshells(4))

!        Intranode parallelization should be possible in the loops below because Interest is thread safe

         print*,"hi I'm inside tens_init_ao_integrals",nshells
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

        function initialize_mocoef_tensor(this,tensor,scalar) result(ierr)
         use, intrinsic:: ISO_C_BINDING
         implicit none
         integer(INTD):: ierr
         class(initialize_mocoef_tensor_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         integer(INTL):: tens_shape(1:MAX_TENSOR_RANK)
         integer(INTL):: tens_signa(1:MAX_TENSOR_RANK)
         integer(INTD):: data_kind,n
         type(C_PTR):: tens_body_p

         ierr=0
         data_kind=tensor%get_data_type(ierr)
         if(ierr.ne.0) then; ierr=-6; return; endif
         call tensor%get_dims(tens_shape,n,ierr)
         if(ierr.ne.0) then; ierr=-5; return; endif
         if(n.eq.2) then !order-2 tensor
          call tensor%get_bases(tens_signa,n,ierr)
          if(ierr.ne.0) then; ierr=-4; return; endif
          tens_body_p=tensor%get_body_ptr(ierr)
          if(ierr.ne.0) then; ierr=-3; return; endif
          ierr=tens_init_mocoef(tens_body_p,data_kind,tens_shape(1:n),tens_signa(1:n))
          if(ierr.ne.0) then; ierr=-2; return; endif
         else
          ierr=-1
         endif
         return
        end function initialize_mocoef_tensor

        integer function tens_init_mocoef(tens_body, data_kind, tens_shape, tens_base) result (ierr)
        !Initialization routine that is used to fill a tensor with mo coefficients

         use exacorr_datatypes
         use exacorr_global
         use, intrinsic:: ISO_C_BINDING

         implicit none

         type(C_PTR), value:: tens_body                     !in: C pointer to the tensor body (tensor elements)
         integer(C_INT), value:: data_kind                  !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
         integer(C_LONG_LONG), intent(in):: tens_shape(1:*) !in: tensor shape (extent of each tensor dimension)
         integer(C_LONG_LONG), intent(in):: tens_base(1:*)  !in: base offsets for each tensor dimension)

         integer:: dimensions(2),nao,nmo
         integer :: first_ao, first_mo
         complex(8),pointer :: mocoef_tens(:,:)
         complex(8),parameter :: ZERO=(0.D0,0.D0)

         ierr=0 !set error code to success
         if(data_kind /= EXA_DATA_KIND_C8) then
           ierr=-1 !invalid data kind requested
           return
         endif

         ! map the arguments to names used in the context of mocoef retrieval
         first_ao = tens_base(1)
         nao = tens_shape(1)
         first_mo = tens_base(2)
         nmo = tens_shape(2)
         dimensions(1) = nao
         dimensions(2) = nmo
         call c_f_pointer(tens_body,mocoef_tens,dimensions) !map C pointer to 2-dimensional array

!        Placeholder: just initialize to zero
         mocoef_tens = ZERO

        end function tens_init_mocoef

end module exacorr_tensor_methods
