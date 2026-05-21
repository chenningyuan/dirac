module tapp
! This module contains the fortran bindings of TAPP standard tensor interface
! Written by Jan Brandejs, Dec 2024

    use, intrinsic:: ISO_C_BINDING
    ! use TAPP_create_tensor_info
    
    implicit none
    
    ! TAPP_create_executor
    interface
    function TAPP_create_executor(executor_) BIND(C, name='TAPP_create_executor') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR) :: executor_ 
      integer(c_int) :: TAPP_create_executor !, intent(out)
    END function TAPP_create_executor
    END interface
    
    ! TAPP_destroy_executor
    interface
    function TAPP_destroy_executor(executor_) BIND(C, name='TAPP_destroy_executor') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: executor_ 
      integer(c_int) :: TAPP_destroy_executor !, intent(out)
    END function TAPP_destroy_executor
    END interface
    
    ! TAPP_create_handle
    interface
    function TAPP_create_handle(handle_) BIND(C, name='TAPP_create_handle') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR) :: handle_ 
      integer(c_int) :: TAPP_create_handle !, intent(out)
    END function TAPP_create_handle
    END interface
    
    ! TAPP_destroy_handle
    interface
    function TAPP_destroy_handle(handle_) BIND(C, name='TAPP_destroy_handle') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR), value :: handle_ 
      integer(c_int) :: TAPP_destroy_handle !, intent(out)
    END function TAPP_destroy_handle
    END interface

    ! TAPP_create_tensor_info
    interface
    function TAPP_create_tensor_info(info_A_, handle, dtype_, nmode_A_, extents_A_, strides_A_) BIND(C, name='TAPP_create_tensor_info') ! result(outx)
      use, intrinsic :: ISO_C_BINDING, only: C_PTR, c_int
      implicit none 
      
      type(C_PTR) :: info_A_ 
      type(C_PTR), value :: handle, extents_A_, strides_A_ 
      integer(c_int), intent(in), value :: dtype_ ! TAPP_datatype
      integer(c_int), intent(in), value :: nmode_A_ 
      integer(c_int) :: TAPP_create_tensor_info !, intent(out)
    END function TAPP_create_tensor_info
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

end module tapp
