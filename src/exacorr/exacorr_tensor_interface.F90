module exacorr_tensor_interface

! This module defines the unified tensor interface
! Written by Jan Brandejs, Andre S.P. Gomes and Lucas Visscher, June 2024

        use exatensor
        use talsh
!        use tapp

        use, intrinsic:: ISO_C_BINDING

        implicit none
        private
!PARAMETERS:
!Generic:
        !integer, parameter, private:: REALD=8 ! default real for this module
        ! aspg: to check whether definitions such as integer(4) :: T_INTL=8 do not clash with imported modules
        ! jan: this is defined in exatensor_placeholders
        
!DERIVED TYPES:
        public :: tensor_data_type_real_8, tensor_data_type_complex_8
        enum, bind(c)
                enumerator :: tensor_data_type = 0 ! name of the enum
                enumerator :: tensor_data_type_real_double = 1
                enumerator :: tensor_data_type_complex_double = 2
                enumerator :: tensor_data_type_int64 = 3
        end enum

!Unified tensor type
        type, public  :: tensor_t
         !library_tensor_instance_pointer
         type(C_PTR) :: tens_info_tapp        ! executor 1       !Do not access in your code directly outside this file (and tapp_utils)!
         type(talsh_tens_t) :: tens_talsh     ! executor 2       !Do not access in your code directly outside this file. Will be deprecated soon.
         type(tens_rcrsv_t) :: tens_exatensor ! executor 3       !Do not access in your code directly outside this file. Will be deprecated soon.

         ! TAPP reference implementation does not store data
         complex(c_double_complex), allocatable, target :: tapp_tens_data_complex(:) !Do not access in your code outside this file (and tapp_utils)!
         real(c_double),            allocatable, target :: tapp_tens_data_real(:)    !Do not access in your code outside this file (and tapp_utils)!
         integer(c_int64_t),        allocatable, target :: tapp_tens_data_int64(:)   !Do not access in your code outside this file (and tapp_utils)!


         integer(kind(tensor_data_type)) :: data_type

         integer, public:: which_executor = 1 !set this default value by CMAKE. 1=tapp, 2=talsh, 3=exatensor
         
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
         type(c_ptr) :: struct_tapp

        end type index_symmetry_structure

        public contract

        contains

        ! old: contract performs D = alpha A * B (+ beta C ) in einsum notation with possible complex conjugates
        ! old: function einsum(expresssion, d, a, b, alpha, beta, c) result(ierr) ! The case C = D is supported. 

        ! Performs D += scale A * B, or D += scale A in talsh einsum notation with possible 
        ! Complex conjugate is denoted by + like A+(i,j,k). 
        function t_einsum(expresssion, d, a, b, scale) result(ierr)  
        !Import :: tensor_t
                character(len = *), intent(in)   :: expresssion
                type(tensor_t), intent(inout)    :: d, a, b
                complex(8), intent(in)           :: scale

                integer(INTD)                    :: ierr, t_rank
                integer(c_int64_t)                    :: tapp_err
                type(tensor_t)                   :: beta_tensor
                character(len = *)               :: addition_str, d_indices, a_indices, b_indices, c_indices  
                logical                          :: a_conj, b_conj, c_conj 
                integer(INTL),allocatable :: t_dim(:),subspace_id(:)
                integer(INTD),allocatable :: space_id(:)

                integer(INTD) :: t_einsum
                ! here the einsum_str has to be parsed to get the indices of D (addition_str). e.g. "D(m1,e1)+=C(m1,e1)*beta()"
                ierr=standard_parse_einsum(expresssion, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj)
                ierr=indices2addition_str(d_indices, c_conj, addition_str)


                SELECT CASE(a%which_executor)
                CASE(1)
                        IF (alpha%im .ne. 0.D0 .and. beta%im .ne. 0.D0) THEN
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta, c, c_indices, c_conj)
                        ELSE IF (alpha%im .ne. 0.D0 .and. beta%im .eq. 0.D0) THEN
                                ierr=standard_contract(d, d_indices, &
                                        alpha, a, a_indices, a_conj, b, b_indices, b_conj, beta%re, c, c_indices, c_conj)
                        ELSE IF (alpha%im .eq. 0.D0 .and. beta%im .ne. 0.D0) THEN
                                ierr=standard_contract(d, d_indices, &
                                        alpha%re, a, a_indices, a_conj, b, b_indices, b_conj, beta, c, c_indices, c_conj)
                        ELSE IF (alpha%im .eq. 0.D0 .and. beta%im .eq. 0.D0) THEN
                                ierr=standard_contract(d, d_indices, &
                                        alpha%re, a, a_indices, a_conj, b, b_indices, b_conj, beta%re, c, c_indices, c_conj)
                        END IF
                        ierr = tapp_err
                CASE(2)
                        ierr=talsh_tensor_contract(einsum_str,d%tens_talsh,a%tens_talsh,b%tens_talsh, scale=alpha)

                        IF (beta .ne. (0.D0,0.D0)) THEN
                                ierr=talsh_tensor_add(addition_str, d%tens_talsh,c%tens_talsh, scale=beta)
                        END IF
                CASE(3)
                        ierr=exatns_tensor_contract(einsum_str, d%tens_exatensor, a%tens_exatensor, b%tens_exatensor, prefactor=alpha)
                              
                        IF (beta .ne. (0.D0,0.D0)) THEN
                                t_rank=tensor%get_rank(ierr)
                                allocate(t_dim(t_rank))
                                call c%get_dims(t_dim,t_rank,ierr)
                                ! get space id
                                allocate(space_id(t_rank),subspace_id(t_rank))
                                call c%get_space_ids(space_id,subspace_id,ierr)

                                ! addition
                                ! use add through the interface instead? why not.
                                ! beta_tensor scalar is just a workaround that replaces exatensor's problematic addition by a contraction
                                ierr=exatns_tensor_create(beta_tensor%tens_exatensor,"aux_beta_tensor",space_id,t_dim,a%data_type)
                                ierr=exatns_tensor_init(beta_tensor%tens_exatensor,beta)
                                ierr=exatns_tensor_contract(addition_str, d%tens_exatensor, a%tens_exatensor, beta_tensor%tens_exatensor) 
                                ierr=exatns_tensor_destroy(beta_tensor%tens_exatensor)
                        END IF
                END SELECT

        t_einsum = ierr
        end function t_einsum

        function indices2addition_str(d_indices, c_conj, addition_str) result(ierr)
                integer(INTD)  :: ierr
                character(len = *), intent(in)   :: d_indices
                logical , intent(in)             :: c_conj     
                character(len = *), intent(out)  :: addition_str

                ierr = -1
     
        end function

        ! external tensor body initialization/transformation (unary) method
        INTERFACE tensor_method_register
                FUNCTION tensor_method_register_tapp(method_name, method) result(ierr)
                        character(*), intent(in):: method_name         
                        integer(T_INTD)                    :: ierr
                        interface 
                                integer(T_INTD) function method(a,b,i,j)  !some function
                                        integer, intent(in) :: a,b,i,j
                                end function
                        end interface

                        ierr = tapp_method_register(method_name,method)
                END FUNCTION

                FUNCTION tensor_method_register_exa(method_name, method) result(ierr)
                        character(*), intent(in):: method_name            
                        class(tens_method_uni_t), intent(inout):: method  !exatensor
                        integer(T_INTD)                    :: ierr
                        ! procedure(tens_method_uni_i), public, deferred:: apply ! possible implementation

                        ierr = exatns_method_register(method_name,method)
                END FUNCTION
        END INTERFACE

        ! this will be a part of the standard interface, implemented inside the c code of the c-api, but also here for practical reasons in case of problems with the compilation of the c part
        ! FUNCTION standard_parse_einsum(einsum_str, d_indices, a_indices, a_conj, b_indices, b_conj, c_indices, c_conj) result (ierr)
        !         integer :: ierr
        !         character(len = *), intent(in)   :: einsum_str     
        !         character(len = *), intent(out)  :: d_indices, a_indices, b_indices 
        !         logical , intent(out)            :: a_conj, b_conj
        !         character(len = *), intent(out), optional  :: c_indices  
        !         logical , intent(out), optional  :: c_conj
        !
        ! should this function be a feature of the standard c-api instead?
        !        ierr = -1
     
        ! END FUNCTION 

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

        !tensor method reset

end module exacorr_tensor_interface
