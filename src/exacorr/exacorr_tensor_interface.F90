module exacorr_tensor_interface

! This module defines the unified tensor interface
! Written by Jan Brandejs, Andre S.P. Gomes and Lucas Visscher, June 2024

        use exatensor
        use talsh
        use tapp
        use tapp_top

        use, intrinsic:: ISO_C_BINDING

        implicit none
        private
!PARAMETERS:
!Generic:
        !integer, parameter, private:: REALD=8 ! default real for this module
        ! aspg: to check whether definitions such as integer(4) :: T_INTL=8 do not clash with imported modules
        ! jan: this is defined in exatensor_placeholders
        
!DERIVED TYPES:
        ! we follow TAPP_datatype enumeration

!Unified tensor type
        type, public  :: tensor_t
         !library_tensor_instance_pointer
         type(tapp_tensor_t) :: tens_tapp     ! executor 1       !Do not access in your code directly outside this file (and tapp_top)!
         type(talsh_tens_t) :: tens_talsh     ! executor 2       !Do not access in your code directly outside this file. Will be deprecated soon.
         type(tens_rcrsv_t) :: tens_exatensor ! executor 3       !Do not access in your code directly outside this file. Will be deprecated soon.

         integer(kind(TAPP_datatype)) :: datatype

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

        ! Most commmon use case:
        ! Performs D += scale A * B, or D += scale A in talsh einsum notation
        ! Complex conjugate is denoted by + like A+(i,j,k). 
        ! if C is defined (present): this function performs D = scale A * B + scale_c C.
        ! if C not defined, this function assumes C=D and scale_c=1
        ! if B and C not defined, performs D += scale A
        ! if neither of A,B,C is defined, sets each element of D to scale value
        function t_einsum(expresssion, d, a, b, c, scale, scale_c)  
          !Import :: tensor_t
                character(len = *), intent(in)   :: expression
                type(tensor_t), intent(inout)    :: d
                type(tensor_t), intent(inout), optional    :: a, b, c
                complex(8), intent(in), optional          :: scale, scale_c
                integer(INTD) :: t_einsum ! return value
                ! internal
                integer(INTD)                    :: ierr
                integer(c_int64_t)               :: tapp_err

                ! tapp
                integer(c_int) :: op_A, op_B, op_C, op_D ! 0 .. identitity, 1 .. element-wise complex conjugate
                integer(c_int64_t), allocatable, target :: arr_idx_A(:), arr_idx_B(:), arr_idx_C(:), arr_idx_D(:)
                complex(c_double_complex)        :: alpha, beta
                character(len=(len_trim(expression)+6))               :: contraction_str,addition_str

                ! exatensor
                type(tensor_t)                   :: one_tensor


                ! here the einsum_str has to be parsed to get indices and throw errors if any
                ierr=parse_t_einsum(expression, arr_idx_A, op_A, arr_idx_B, op_B, arr_idx_C, op_C, arr_idx_D, op_D)
              
                ! check number of indices
                 ! to.do

                ! check extents
                 ! to.do
                
                alpha = 1.0
                beta = 1.0
                if (present(scale)) alpha = scale 
                if (present(scale_c)) beta = scale_c
                
                SELECT CASE(a%which_executor)
                CASE(1) ! tapp
                        ! to.do: what if A/B/C is scalar? empty idx arrays and nmodes=0?
                        if     (present(a) .and. present(b) .and. present(c)) then
                          tapp_err = tapp_einsum(d%tens_tapp, arr_idx_D, op_D, a%tens_tapp, arr_idx_A, op_A, &
                            b%tens_tapp, arr_idx_B, op_B, c%tens_tapp, arr_idx_C, op_C, alpha=alpha, beta=beta)
                        elseif (present(a) .and. present(b)) then
                          tapp_err = tapp_einsum(d%tens_tapp, arr_idx_D, op_D, a%tens_tapp, arr_idx_A, op_A, &
                            b%tens_tapp, arr_idx_B, op_B, alpha=alpha, beta=beta)
                        elseif (present(a)) then
                          tapp_err = tapp_einsum(d%tens_tapp, arr_idx_D, op_D, a%tens_tapp, arr_idx_A, op_A, &
                            alpha=alpha, beta=beta)
                        else
                          tapp_err = tapp_einsum(d%tens_tapp, arr_idx_D, op_D, alpha=alpha, beta=beta)
                        endif
                        ierr = tapp_err
                CASE(2) ! talsh
                        if (present(a) .and. present(b) .and. present(c)) then
                          ! the case D=C has to be supported 
                          ! to.do check if complex conjugate works in copy 
                          ierr=split_DABC_expression_str(expression, op_C, contraction_str, addition_str)
                          ierr=talsh_tensor_copy(trim(addition_str), d%tens_talsh,c%tens_talsh)
                          ierr=talsh_tensor_scale(d%tens_talsh, beta)
                          ierr=talsh_tensor_contract(trim(contraction_str),d%tens_talsh,a%tens_talsh,b%tens_talsh, scale=alpha)
                        elseif (present(a) .and. present(b)) then
                          ierr=talsh_tensor_contract(trim(expression),d%tens_talsh,a%tens_talsh,b%tens_talsh, scale=alpha)
                        elseif (present(a)) then
                          ierr=talsh_tensor_add(trim(expression), d%tens_talsh,a%tens_talsh, scale=alpha)
                        else
                          ierr=talsh_tensor_init(d%tens_talsh, alpha)
                        endif

                        ! if (a%datatype == 1) then
                        ! IF (alpha%im .ne. 0.D0 .and. beta%im .ne. 0.D0) THEN

                CASE(3) ! exatensor
                        if (present(a) .and. present(b) .and. present(c)) then
                          ! the case D=C has to be supported 
                          ! to.do check if complex conjugate works in copy 
                          ierr=split_DABC_expression_str(expression, op_C, contraction_str, addition_str)
                          ierr=exatns_tensor_copy(trim(addition_str), d%tens_talsh,c%tens_talsh)
                          ierr=exatns_tensor_scale(d%tens_talsh, beta)
                          ierr=exatns_tensor_contract(trim(contraction_str),d%tens_exatensor,a%tens_exatensor,b%tens_exatensor, prefactor=alpha)
                        elseif (present(a) .and. present(b)) then
                          ierr=exatns_tensor_contract(trim(expression),d%tens_exatensor,a%tens_exatensor,b%tens_exatensor, prefactor=alpha)
                        elseif (present(a)) then
                          ierr=exatns_tensor_create(one_tensor%tens_exatensor,"autokern_one_tensor",EXA_DATA_KIND_C8)
                          ierr=exatns_tensor_init(one_tensor,1.0d0)
                          addition_str = trim(expression) // "*O01()"
                          addition_str=exatns_tensor_contract(trim(addition_str), d%tens_exatensor, a%tens_exatensor, &
                                         one_tensor%tens_exatensor, prefactor=alpha) 
                          ierr=exatns_tensor_destroy(one_tensor%tens_exatensor)
                        else
                          ierr==exatns_tensor_init(d%tens_exatensor, alpha)
                        endif
                              
                END SELECT

        if (allocated(arr_idx_A)) deallocate(arr_idx_A)
        if (allocated(arr_idx_B)) deallocate(arr_idx_B)
        if (allocated(arr_idx_C)) deallocate(arr_idx_C)
        if (allocated(arr_idx_D)) deallocate(arr_idx_D)

        t_einsum = ierr
        end function t_einsum

        function split_DABC_expression_str(expression, op_C, contraction_str, addition_str)
                character(len = *), intent(in)   :: expression
                integer(c_int), intent(in) :: op_C ! 0 .. identitity, 1 .. element-wise complex conjugate
                character(len = len_trim(expression)), intent(out)  :: contraction_str, addition_str
                integer(INTD)  :: split_DABC_expression_str
                integer :: i, pos, curr, curr2, eq_found
                
                pos = -1
                do i=1:len_trim(expression)
                  curr = len_trim(expression)-i+1
                  if (expression(curr:curr) == '+') then
                    if (op_C == 1 .and. pos /= -1) then
                      pos = curr
                      exit
                    elseif (op_C == 1 .and. pos == -1) then
                      pos = curr
                    else
                      pos = curr
                      exit
                    endif
                enddo
                if (pos == -1) then ! assumes that + is used for addition
                  write(*,'("  Error. + not found in expression D = A*B+C.")')
                  stop 45
                end if
                contraction_str = ""
                eq_found = 0
                do i=1:(pos-1)
                  contraction_str(i:i)=expression(i:i)
                  if (eq_found == 0) addition_str(i:i)=expression(i:i)
                  if (expression(i:i) == '=') eq_found = i
                enddo
                addition_str = ""
                do i=1:(len_trim(expression)-pos)
                  curr = eq_found + i
                  curr2 = pos + i
                  addition_str(curr:curr)=expression(curr2:curr2)
                  if (addition_str(curr:curr) == '*') then ! assumes that + operation is on the right of *
                    write(*,'("  Error. in expression D = A*B+C we assume + on the right of *.")')
                    stop 46
                  end if
                enddo

                split_DABC_expression_str = 0
        end function
 
        ! copies data and destroys the source talsh tensor, unless staying in talsh
        function talsh_to_tensor_t(tens_talsh, out_tens, out_executor)
          implicit none
          type(talsh_tens_t), intent(inout)    :: tens_talsh
          type(tensor_t), intent(inout)    :: out_tens
          integer, intent(in), optional :: out_executor ! If not set value set by Cmake is used. 1=tapp, 2=talsh, 3=exatensor.

          integer             :: DataKind(1),nData
          integer(INTD)       :: ierr, talsh_to_tensor_t
          type(tensor_t)    :: in_tens

          ierr=talsh_tensor_data_kind(tens_talsh,nData,DataKind)
          if (DataKind(1) == R8) in_tens%datatype = 1
          if (DataKind(1) == C8) in_tens%datatype = 3

          in_tens%which_executor = 2 
          in_tens%tens_talsh = tens_talsh

          if (present(out_executor)) then
            ierr=convert(in_tens, out_tens, out_executor)
          else
            ierr=convert(in_tens, out_tens)
          endif

          talsh_to_tensor_t = ierr
        end function talsh_to_tensor_t

        ! copies data and destroys the source talsh tensor, unless staying in talsh
        function tensor_t_to_talsh(in_tens, tens_talsh)
          implicit none
          type(tensor_t), intent(inout)    :: in_tens
          type(talsh_tens_t), intent(inout)    :: tens_talsh

          integer(INTD)       :: ierr, tensor_t_to_talsh
          type(tensor_t)    :: out_tens
          
          ierr=convert(in_tens, out_tens, out_executor=2)
          tens_talsh = out_tens%tens_talsh

          tensor_t_to_talsh = ierr
        end function tensor_t_to_talsh

        ! copies data and destroys the source tensor, unless staying in the same executor
        function convert(in_tens, out_tens, out_executor) ! If out_executor not set value set by Cmake is used.
          implicit none
          type(tensor_t), intent(inout)    :: in_tens
          type(tensor_t), intent(inout)    :: out_tens
          integer, intent(in), optional:: out_executor 
          integer(INTD)       :: ierr, convert
          integer             :: in_out_exec
          ! TAPP_create_tensor_info
          type(C_PTR) :: info
          integer(c_int) :: dtype, nmodes, tapp_err
          integer(c_int64_t), allocatable, target :: extents(:), strides(:)
          integer(8) :: i, k, stride, numel
          integer(INTD), allocatable, target :: out_dims(:)

          complex(c_double_complex), allocatable, target :: arr_c(:)
          real(c_double), allocatable, target :: arr_r(:)
          integer(INTD), allocatable, target :: lin_dim(:)
          complex(8), pointer  :: lin_tens_c(:)
          real(8), pointer  :: lin_tens_r(:)
          type(C_PTR) :: talshsptr
           
          if (present(out_executor)) then
            out_tens%which_executor = out_executor
          endif

          out_tens%datatype = in_tens%datatype
          
          in_out_exec = 100 * in_tens%which_executor + out_executor 
          SELECT CASE( in_out_exec) ! e.g. 201 means convert from exec 2 to 1 
            CASE(102) 
              nmodes = TAPP_get_nmodes(in_tens%tens_tapp%tensor_info)
              
              ! to.do replace with unified get_dims
              allocate(extents(nmodes))
              call TAPP_get_extents(in_tens%tens_tapp%tensor_info, c_loc(extents(1)))
              numel = 1
              do i=1,nmodes
                numel = numel * extents(i)
              enddo
              
              allocate(out_dims(nmodes))
              out_dims = int(extents, kind=INTD)
              allocate(lin_dim(1))
              lin_dim(1) = numel

              if (out_tens%datatype == 1 ) then
                ierr=talsh_tensor_construct(out_tens%tens_talsh, R8, out_dims, init_val=ZERO)
                ierr=talsh_tensor_get_body_access(out_tens%tens_talsh,talshsptr,R8,int(0,C_INT),DEV_HOST)
                call c_f_pointer(talshsptr,lin_tens_r,lin_dim)
                call dcopy(numel, in_tens%tens_tapp%data_real, 1, lin_tens_r, 1)
                deallocate(in_tens%tens_tapp%data_real)

              elseif (out_tens%datatype == 3 ) then
                ierr=talsh_tensor_construct(out_tens%tens_talsh, C8, out_dims, init_val=ZERO)
                ierr=talsh_tensor_get_body_access(out_tens%tens_talsh,talshsptr,C8,int(0,C_INT),DEV_HOST)
                call c_f_pointer(talshsptr,lin_tens_c,lin_dim)
                call zcopy(numel, in_tens%tens_tapp%data_complex, 1, lin_tens_c, 1)
                deallocate(in_tens%tens_tapp%data_complex)
              endif

              tapp_err = TAPP_destory_tensor_info(in_tens%tens_tapp%tensor_info)

            CASE(201) 
              nmodes = talsh_tensor_rank(in_tens%tens_talsh)
              allocate(out_dims(nmodes))
              ierr = talsh_tensor_dimensions(in_tens%tens_talsh,nmodes,out_dims)
              
              allocate(extents(nmodes), strides(nmodes))
              extents = int(out_dims, kind=8)
              
              ! to.do replace with unified create_tensor
              numel = 1
              do i=1,nmodes
                stride = 1
                do k=1,i-1
                  stride = stride * extents(k)
                enddo
                strides(i) = stride
                numel = numel * extents(i)
              enddo
              
              dtype = out_tens%datatype
              
              tapp_err = TAPP_create_tensor_info(out_tens%tens_tapp%tensor_info, dtype, nmodes, &
                            c_loc(extents(1)), c_loc(strides(1)))
              if (tapp_err.ne.0) then
                write(*,'("  Error. Creating TAPP tensor info failed with error code: ",i11,".")') tapp_err
                stop 40
              end if

              ! copy data
              allocate(lin_dim(1))
              lin_dim(1) = numel
  
              if (out_tens%datatype == 1 ) then
                ierr=talsh_tensor_get_body_access(in_tens%tens_talsh,talshsptr,R8,int(0,C_INT),DEV_HOST)
                call c_f_pointer(talshsptr,lin_tens_r,lin_dim)
                allocate(arr_r(numel)) 
                call dcopy(numel, lin_tens_r, 1, arr_r, 1)
                out_tens%tens_tapp%data_real = arr_r
              elseif (out_tens%datatype == 3 ) then
                ierr=talsh_tensor_get_body_access(in_tens%tens_talsh,talshsptr,C8,int(0,C_INT),DEV_HOST)
                call c_f_pointer(talshsptr,lin_tens_c,lin_dim)
                allocate(arr_c(numel)) 
                call zcopy(numel, lin_tens_c, 1, arr_c, 1)
                out_tens%tens_tapp%data_complex = arr_c
              endif

              ierr=talsh_tensor_destruct(in_tens%tens_talsh)

            CASE(202) 
              out_tens%tens_talsh = tens_talsh
            CASE default 
              write(*,'("  Error. Conversion from executor ",i11," to ",i11," not yet implemented.")') &
                in_tens%which_executor, out_executor
              stop 55
          END SELECT

          if (allocated(extents)) deallocate(extents)
          if (allocated(strides)) deallocate(strides)
          if (allocated(out_dims)) deallocate(out_dims)
          if (allocated(lin_dim)) deallocate(lin_dim)
          if (allocated(lin_tens_r)) deallocate(lin_tens_r)
          if (allocated(lin_tens_c)) deallocate(lin_tens_c)

          convert = ierr
        end function convert


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
