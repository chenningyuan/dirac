module exacorr_tensor_interface

! This module defines the unified tensor interface
! Written by Jan Brandejs and Luuk Visscher, June 2024

        use exatensor
        use talsh
!        use standard_iterface

        use, intrinsic:: ISO_C_BINDING

        implicit none
        private
!PARAMETERS:
!Generic:
        !integer, parameter, private:: REALD=8 ! default real for this module
        integer(4) INTL=8
        integer(4) INTD=4
        
!DERIVED TYPES:
        public :: tensor_data_type_real_8, tensor_data_type_complex_8
        enum, bind(c)
                enumerator :: tensor_data_type = 0 ! name of the enum
                enumerator :: tensor_data_type_real_8 = 1
                enumerator :: tensor_data_type_complex_8 = 2
        end enum

!Unified tensor type
        type, public  :: tensor_t
         !library_tensor_instance_pointer
         type(c_ptr) :: tens_standard
         type(tens_rcrsv_t) :: tens_exatensor  
         type(talsh_tens_t) :: tens_talsh

         integer(kind(tensor_data_type)) :: data_type

         !tens_rcrsv_t has to have a placeholder too

         integer, public:: which_library = 1 !set this default value by CMAKE
         
         !context, metadata
                !dimensions?
                !indices?

        end type tensor_t

        

        type, public  :: index_label
         !tensor index label
         character :: val = 'a'
        end type index_label

        type, public  :: index_symmetry_structure
         !index space
         !defines the symmetry sructure of a index
         !like AO/MOspace. E.g. MOspace(virtual), MOspace(occupied)
         ! h_space_t in exatensor
         class(h_space_t), pointer :: struct_exatensor
         type(c_ptr) :: struct_talsh
         type(c_ptr) :: struct_standard

        end type index_symmetry_structure

        contains

        !performs D += alpha A * B (+ beta C ), with possible complex conjugates
        !implement +beta C as part of this function?
        INTERFACE contract
                SUBROUTINE contract_r(einsum_str, d, alpha, a, b) ! D += alpha A * B, with alpha real
                        Import :: tensor_t
                        character(len = *), intent(in)   :: einsum_str
                        type(tensor_t), intent(inout)      :: d
                        real(8), intent(in)              :: alpha
                        type(tensor_t), intent(in)         :: a,b
                        integer(INTD)                    :: ierr
                        character(len = *)               :: d_indices, a_indices, b_indices  
                        logical                          :: a_conj, b_conj

                        SELECT CASE(a%which_library)
                        CASE(1)
                                ierr=standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj)
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta)
                        CASE(2)
                                ierr=exatns_tensor_contract(einsum_str, d, a, b, prefactor=alpha)
                        CASE(3)
                                ierr=talsh_tensor_contract(einsum_str,d,a,b, scale=alpha)
                        END SELECT

                END SUBROUTINE

                SUBROUTINE contract_c(einsum_str, d, alpha, a, b) ! D += alpha A * B, with alpha complex
                        Import :: tensor_t
                        character(len = *), intent(in)   :: einsum_str
                        type(tensor_t), intent(inout)      :: d
                        complex(8), intent(in)           :: alpha
                        type(tensor_t), intent(in)         :: a,b
                        integer(INTD)                    :: ierr
                        character(len = *)               :: d_indices, a_indices, b_indices  
                        logical                          :: a_conj, b_conj

                        SELECT CASE(a%which_library)
                        CASE(1)
                                ierr=standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj)
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta)
                        CASE(2)
                                ierr=exatns_tensor_contract(einsum_str, d, a, b, prefactor=alpha)
                        CASE(3)
                                ierr=talsh_tensor_contract(einsum_str,d,a,b, scale=alpha)
                        END SELECT

                END SUBROUTINE

                !the following routine is present four times in this file, with different data types of alpha and beta
                SUBROUTINE contract_r_r(einsum_str, d, alpha, a, b, beta, c) ! D += alpha A * B + beta C, with alpha real, beta real
                        Import :: tensor_t
                        character(len = *), intent(in)   :: einsum_str
                        type(tensor_t), intent(inout)      :: d
                        real(8), intent(in)              :: alpha
                        real(8), intent(in)              :: beta
                        type(tensor_t), intent(in)         :: a, b, c
                        integer(INTD)                    :: ierr
                        type(tensor_t)                     :: beta_tensor
                        character(len = *)               :: addition_str, d_indices, a_indices, b_indices, c_indices  
                        logical                          :: a_conj, b_conj, c_conj  
                        
                        SELECT CASE(a%which_library)
                        CASE(1)
                                ierr=standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj)
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta, c, c_indices, c_conj)
                        CASE(2)
                                ierr=exatns_tensor_contract(einsum_str, d, a, b, prefactor=alpha)
                                
                                ! addition
                                ! use add through the interface instead?
                                ! beta_tensor scalar is just a workaround that replaces exatensor's problematic addition by a contraction
                                ierr=exatns_tensor_create(beta_tensor,"aux_beta_tensor",a%data_type)
                                ierr=exatns_tensor_init(beta_tensor,beta)
                                ! here the einsum_str has to be parsed to get the indices of D (addition_str). e.g. "D(m1,e1)+=C(m1,e1)*beta()"
                                ! parse using the standard to avoid code duplicate? (standard api will have to do this anyway)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=exatns_tensor_contract(addition_str, d, a, beta_tensor) 
                        CASE(3)
                                ierr=talsh_tensor_contract(einsum_str,d,a,b, scale=alpha)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=talsh_tensor_add(addition_str, d,c, scale=beta)
                END SELECT

                END SUBROUTINE

                SUBROUTINE contract_c_c(einsum_str, d, alpha, a, b, beta, c) ! D += alpha A * B + beta C, with alpha complex, beta complex 
                        Import :: tensor_t
                        character(len = *), intent(in)   :: einsum_str
                        type(tensor_t), intent(inout)      :: d
                        complex(8), intent(in)           :: alpha
                        complex(8), intent(in)           :: beta
                        type(tensor_t), intent(in)         :: a, b, c
                        integer(INTD)                    :: ierr
                        type(tensor_t)                     :: beta_tensor
                        character(len = *)               :: addition_str, d_indices, a_indices, b_indices, c_indices  
                        logical                          :: a_conj, b_conj, c_conj  

                        SELECT CASE(a%which_library)
                        CASE(1)
                                ierr=standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj)
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta, c, c_indices, c_conj)
                        CASE(2)
                                ierr=exatns_tensor_contract(einsum_str, d, a, b, prefactor=alpha)
                                
                                ! addition
                                ierr=exatns_tensor_create(beta_tensor,"aux_beta_tensor",a%data_type)
                                ierr=exatns_tensor_init(beta_tensor,beta)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=exatns_tensor_contract(addition_str, d, a, beta_tensor)
                        CASE(3)
                                ierr=talsh_tensor_contract(einsum_str,d,a,b, scale=alpha)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=talsh_tensor_add(addition_str, d,c, scale=beta)
                END SELECT
                END SUBROUTINE

                SUBROUTINE contract_r_c(einsum_str, d, alpha, a, b, beta, c) ! D += alpha A * B + beta C, with alpha real, beta complex  
                        Import :: tensor_t
                        character(len = *), intent(in)   :: einsum_str         
                        type(tensor_t), intent(inout)      :: d
                        real(8), intent(in)             :: alpha
                        complex(8), intent(in)           :: beta
                        type(tensor_t), intent(in)         :: a, b, c
                        integer(INTD)                    :: ierr
                        type(tensor_t)                     :: beta_tensor
                        character(len = *)               :: addition_str, d_indices, a_indices, b_indices, c_indices  
                        logical                          :: a_conj, b_conj, c_conj

                        SELECT CASE(a%which_library)
                        CASE(1)
                                ierr=standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj)
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta, c, c_indices, c_conj)
                                
                        CASE(2)
                                ierr=exatns_tensor_contract(einsum_str, d, a, b, prefactor=alpha)
                                
                                ! addition
                                ierr=exatns_tensor_create(beta_tensor,"aux_beta_tensor",a%data_type)
                                ierr=exatns_tensor_init(beta_tensor,beta)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=exatns_tensor_contract(addition_str, d, a, beta_tensor)
                        CASE(3)
                                ierr=talsh_tensor_contract(einsum_str,d,a,b, scale=alpha)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=talsh_tensor_add(addition_str, d,c, scale=beta)
                END SELECT
                END SUBROUTINE

                SUBROUTINE contract_c_r(einsum_str, d, alpha, a, b, beta, c) ! D += alpha A * B + beta C, with alpha complex, beta real 
                        Import :: tensor_t
                        character(len = *), intent(in)   :: einsum_str         
                        type(tensor_t), intent(inout)      :: d
                        complex(8), intent(in)           :: alpha
                        real(8), intent(in)              :: beta
                        type(tensor_t), intent(in)         :: a, b, c
                        integer(INTD)                    :: ierr
                        type(tensor_t)                     :: beta_tensor
                        character(len = *)               :: addition_str, d_indices, a_indices, b_indices, c_indices  
                        logical                          :: a_conj, b_conj, c_conj

                        SELECT CASE(a%which_library)
                        CASE(1)
                                ierr=standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj)
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta, c, c_indices, c_conj)
                                
                        CASE(2)
                                ierr=exatns_tensor_contract(einsum_str, d, a, b, prefactor=alpha)

                                ! addition
                                ierr=exatns_tensor_create(beta_tensor,"aux_beta_tensor",a%data_type)
                                ierr=exatns_tensor_init(beta_tensor,beta)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=exatns_tensor_contract(addition_str, d, a, beta_tensor)
                        CASE(3)
                                ierr=talsh_tensor_contract(einsum_str,d,a,b, scale=alpha)
                                ierr=standard_parse_einsum(einsum_str, d_indices)
                                ierr=indices2addition_str(d_indices, addition_str)
                                ierr=talsh_tensor_add(addition_str, d,c, scale=beta)
                END SELECT
                END SUBROUTINE
        END INTERFACE
        
        SUBROUTINE standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj)
                character(len = *), intent(in)   :: einsum_str     
                character(len = *), intent(out)  :: d_indices, a_indices, b_indices 
                logical , intent(out)            :: a_conj, b_conj
                character(len = *), intent(out), optional  :: c_indices  
                logical , intent(out), optional  :: c_conj

                !should this function be a feature of the standard c-api instead?
     
        END SUBROUTINE

        SUBROUTINE indices2addition_str(d_indices, addition_str)
                character(len = *), intent(in)   :: d_indices     
                character(len = *), intent(out)  :: addition_str
     
        END SUBROUTINE

        !todo:
        !add
                !add standard
                !add exatensor
                !add talsh
        !initialize
                !contract standard
                !contract exatensor
                !contract talsh
        !copy
                !contract standard
                !contract exatensor
                !contract talsh

                !copy from talsh to standard
                !tensor type mapping from one library to another
                
        !optional tensor method register
                !method register standard
                !method register exatensor
                
        !tensor method apply
                !tensor method apply standard
                !tensor method apply exatensor
                !tensor method apply talsh

                !type, extends(tens_method_uni_t), public:: compute_2e_ao_tensor_t
                ! contains
                !  procedure, public:: apply=>compute_2e_ao_tensor
                !end type compute_2e_ao_tensor_t

end module exacorr_tensor_interface
