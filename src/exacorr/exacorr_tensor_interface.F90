!This module defines the unified tensor interface
       module exacorr_tensor_interface

#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        use exatensor
#endif
        use talsh
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private
!PARAMETERS:
!Generic:
        !integer, parameter, private:: REALD=8 ! default real for this module
!DERIVED TYPES:
!Unified tensor type
        type, public  :: tensor
         !library_tensor_instance_pointer
         type(tens_rcrsv_t) :: tens_exatensor  
         type(talsh_tens_t) :: tens_talsh
         type(c_ptr) :: tens_standard
         
         !context, metadata

        end type tensor

        !contract
        !initialize
        !copy
        !tensor methods
        !mapping from one library to another


       end module exacorr_tensor_interface
