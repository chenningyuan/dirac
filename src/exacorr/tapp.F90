module tapp
! This module contains the fortran bindings of TAPP standard tensor interface
! Written by Jan Brandejs, Dec 2024

    use, intrinsic:: ISO_C_BINDING
    ! use TAPP_create_tensor_info
    
    implicit none
    

    ! TAPP_create_tensor_info
    interface
    function TAPP_create_tensor_info(info_A_, dtype_, nmode_A_, extents_A_, strides_A_) BIND(C, name='TAPP_create_tensor_info') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR) :: info_A_ 
      type(C_PTR), value :: extents_A_, strides_A_ 
      integer(c_int), intent(in), value :: dtype_ ! TAPP_datatype
      integer(c_int), intent(in), value :: nmode_A_ 
      integer(c_int) :: TAPP_create_tensor_info !, intent(out)
    END function TAPP_create_tensor_info
    END interface

    ! TAPP_create_tensor
    interface
    function TAPP_create_tensor(info_A_, dtype_, nmode_A_, extents_A_, strides_A_, init_val_) BIND(C, name='TAPP_create_tensor') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none

      type(C_PTR) :: info_A_
      type(C_PTR), value :: extents_A_, strides_A_
      integer(c_int), intent(in), value :: dtype_ ! TAPP_datatype
      integer(c_int), intent(in), value :: nmode_A_
      type(C_PTR), value  :: init_val_
      integer(c_int) :: TAPP_create_tensor !, intent(out)
    END function TAPP_create_tensor
    END interface

    ! TAPP_initialize_tensor
    interface
    function TAPP_initialize_tensor(info_A_, init_val_) BIND(C, name='TAPP_initialize_tensor') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none

      type(C_PTR), value  :: info_A_
      type(C_PTR), value  :: init_val_
      integer(c_int) :: TAPP_initialize_tensor !, intent(out)
    END function TAPP_initialize_tensor
    END interface

    ! TAPP_copy_tensor
    interface
    function TAPP_copy_tensor(info_dest_, info_src_) BIND(C, name='TAPP_copy_tensor') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none

      type(C_PTR), value  :: info_dest_
      type(C_PTR), value  :: info_src_
      integer(c_int) :: TAPP_copy_tensor !, intent(out)
    END function TAPP_copy_tensor
    END interface

    ! TAPP_copy_tensor
    interface
    function TAPP_reshape_tensor(info_dest_, info_src_, nmode, extents) BIND(C, name='TAPP_reshape_tensor') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none

      type(C_PTR), value  :: info_dest_
      type(C_PTR), value  :: info_src_
      type(C_PTR), value :: extents 
      integer(c_int), intent(in), value :: nmode 
      integer(c_int) :: TAPP_reshape_tensor !, intent(out)
    END function TAPP_reshape_tensor
    END interface

    ! TAPP_broadcast_tensor_body_data
    interface
    function TAPP_broadcast_tensor_body_data(info_A_, data_) BIND(C, name='TAPP_broadcast_tensor_body_data') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none

      type(C_PTR), value  :: info_A_
      type(C_PTR), value  :: data_
      integer(c_int) :: TAPP_broadcast_tensor_body_data !, intent(out)
    END function TAPP_broadcast_tensor_body_data
    END interface

    ! TAPP_gather_tensor_body_data
    interface
    function TAPP_gather_tensor_body_data(info_A_, data_) BIND(C, name='TAPP_gather_tensor_body_data') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none

      type(C_PTR), value  :: info_A_
      type(C_PTR), value  :: data_
      integer(c_int) :: TAPP_gather_tensor_body_data !, intent(out)
    END function TAPP_gather_tensor_body_data
    END interface

    ! TAPP_create_tensor_product 
    interface
    function TAPP_create_tensor_product(plan_, handle_, op_A_, info_A_, idx_A_, op_B_, info_B_, idx_B_, &
      op_C_, info_C_, idx_C_, op_D_, info_D_, idx_D_, prec_) BIND(C, name='TAPP_create_tensor_product') 
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: info_A_, info_B_, info_C_, info_D_
      type(C_PTR), value :: handle_ ! TAPP_handle
      type(C_PTR) :: plan_ ! TAPP_tensor_product
      integer(c_int), value :: op_A_, op_B_, op_C_, op_D_ ! 0 .. identitity
      type(C_PTR), value:: idx_A_, idx_B_, idx_C_, idx_D_
      integer(c_int), value :: prec_ !  TAPP_prectype -1 .. default
      integer(c_int) :: TAPP_create_tensor_product 
    END function TAPP_create_tensor_product
    END interface


    ! TAPP_execute_product
    interface
    function TAPP_execute_product(plan_, exec_, status_, alpha_, A_, B_, beta_, C_, D_) BIND(C, name='TAPP_execute_product') 
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: plan_ ! TAPP_tensor_product
      type(C_PTR), value :: exec_ ! TAPP_executor
      type(C_PTR) :: status_ ! TAPP_status
      type(C_PTR), value :: alpha_, beta_ ! scalar void*  ..  (void*)&alpha_
      type(C_PTR), value :: A_, B_, C_ ! data arr's .. (void*)A_
      type(C_PTR), value :: D_ ! data arr's .. (void*)A_
      integer(c_int) :: TAPP_execute_product
    END function TAPP_execute_product
    END interface


    ! TAPP_destroy_tensor_product
    interface
    function TAPP_destroy_tensor_product(plan_) BIND(C, name='TAPP_destroy_tensor_product') 
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: plan_ ! TAPP_tensor_product
      integer(c_int) :: TAPP_destroy_tensor_product 
    END function TAPP_destroy_tensor_product
    END interface


    ! TAPP_destroy_tensor_info
    interface
    function TAPP_destroy_tensor_info(info_) BIND(C, name='TAPP_destroy_tensor_info') 
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: info_ ! TAPP_tensor_info
      integer(c_int) :: TAPP_destroy_tensor_info 
    END function TAPP_destroy_tensor_info
    END interface


    ! TAPP_get_nmodes
    interface
    function TAPP_get_nmodes(info_) BIND(C, name='TAPP_get_nmodes') 
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: info_ ! TAPP_tensor_info
      integer(c_int) :: TAPP_get_nmodes
    END function TAPP_get_nmodes
    END interface


    ! TAPP_get_extents
    interface
    subroutine TAPP_get_extents(info_, extents_) BIND(C, name='TAPP_get_extents') 
      use, intrinsic :: ISO_C_BINDING, only: C_PTR
      implicit none 
      
      type(C_PTR), value :: info_ ! TAPP_tensor_info. Already exists, passed by value.
      type(C_PTR), value :: extents_ 
    END subroutine TAPP_get_extents
    END interface

    ! TAPP_scale_with_denominators
    interface
    function TAPP_scale_with_denominators(info_, n_occ_, n_vir_, eps_occ_, eps_vir_, &
        eps_ijk_) BIND(C, name='TAPP_scale_with_denominators')
      use, intrinsic :: ISO_C_BINDING, only: c_int, C_PTR
      implicit none

      type(C_PTR), value :: info_
      integer(c_int), intent(in), value :: n_occ_
      integer(c_int), intent(in), value :: n_vir_
      type(C_PTR), value :: eps_occ_, eps_vir_
      type(C_PTR), value  :: eps_ijk_
      integer(c_int) :: TAPP_scale_with_denominators
    END function TAPP_scale_with_denominators
    END interface

    ! TAPP_finalize
    interface
    function TAPP_finalize() BIND(C, name='TAPP_finalize')
      use, intrinsic :: ISO_C_BINDING, only: c_int
      implicit none

      integer(c_int) :: TAPP_finalize
    END function TAPP_finalize
    END interface
end module tapp
