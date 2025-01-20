module exacorr_tensor_interface

! This module defines the unified tensor interface
! Written by Jan Brandejs, Andre S.P. Gomes and Lucas Visscher, June 2024

        use exatensor
        use talsh
        use tapp
        use tapp_top

        use, intrinsic:: ISO_C_BINDING

        implicit none
!         private
!PARAMETERS:
!Generic:
        ! aspg: to check whether definitions such as integer(4) :: T_INTL=8 do not clash with imported modules
        ! jan: this is defined in exatensor_placeholders
        
!DERIVED TYPES:

  !Unified tensor type
        type, public  :: tensor_t
         !library_tensor_instance_pointer
         type(tapp_tensor_t) :: tens_tapp     ! executor 1       !Do not access in your code directly outside this file (and tapp_top)!
         type(talsh_tens_t) :: tens_talsh     ! executor 2       !Do not access in your code directly outside this file. Will be deprecated soon.
         type(tens_rcrsv_t) :: tens_exatensor ! executor 3       !Do not access in your code directly outside this file. Will be deprecated soon.

         ! we follow TAPP_datatype enumeration
         integer(kind(TAPP_datatype)) :: datatype

         integer, public:: which_executor = 2 !set this default value by CMAKE. 1=tapp, 2=talsh, 3=exatensor
         
         ! contains
         !   procedure :: copy_tensor_t
         !   generic :: assignment(=) => copy_tensor_t
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

        ! external tensor body initialization/transformation (unary) method
        ! INTERFACE tensor_method_register

!Functions
        public construct_tensor_t
        public t_einsum
        public talsh_to_tensor_t
        public tensor_t_to_talsh
        public convert
        public destroy_tensor_t
        public direct_data_pointer
        public get_extents
        public get_nmodes
        public get_numel
        public is_scalar_tensor_t
        public is_empty_tensor_t
        public initialize_tensor_t
        ! public copy_tensor_t
        ! public get_scalar
        ! public get_real_scalar
        ! public tensor_t_equal

        contains

        ! Most commmon use case:
        ! Performs D += scale A * B, or D += scale A in talsh einsum notation
        ! Complex conjugate is denoted by + like A+(i,j,k). 
        ! if C is defined (present): this function performs D = scale A * B + scale_c C.
        ! if C not defined, this function assumes C=D and scale_c=1
        ! if B and C not defined, performs D += scale A
        ! if neither of A,B,C is defined, sets each element of D to scale value
        function t_einsum(expression, d, a, b, c, scale, scale_c)  
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
                complex(8), parameter :: ONE_c=(1.D0,0.D0)

                ! exatensor
                type(tens_rcrsv_t)                   :: one_tensor


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
                          ! tapp_err = tapp_einsum(d%tens_tapp, arr_idx_D, op_D, alpha=alpha, beta=beta)
                          tapp_err = 0
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
                        ! else
                        !   ierr=talsh_tensor_init(d%tens_talsh, alpha)
                        endif

                        ! IF (alpha%im .ne. 0.D0 .and. beta%im .ne. 0.D0) THEN

                CASE(3) ! exatensor
                        if (present(a) .and. present(b) .and. present(c)) then
                          ! the case D=C has to be supported 
                          ! to.do check if complex conjugate works in copy 
                          ierr=split_DABC_expression_str(expression, op_C, contraction_str, addition_str)
                          ierr=exatns_tensor_copy(trim(addition_str), d%tens_exatensor,c%tens_exatensor)
                          ierr=exatns_tensor_scale(d%tens_exatensor, beta)
                          ierr=exatns_tensor_contract(trim(contraction_str),d%tens_exatensor,a%tens_exatensor, &
                                 b%tens_exatensor, prefactor=alpha)
                        elseif (present(a) .and. present(b)) then
                          ierr=exatns_tensor_contract(trim(expression),d%tens_exatensor,a%tens_exatensor, &
                                 b%tens_exatensor, prefactor=alpha)
                        elseif (present(a)) then
                          ierr=exatns_tensor_create(one_tensor,"autokern_one_tensor",EXA_DATA_KIND_C8)
                          ierr=exatns_tensor_init(one_tensor,ONE_c)
                          addition_str = trim(expression) // "*O01()"
                          ierr=exatns_tensor_contract(trim(addition_str), d%tens_exatensor, a%tens_exatensor, &
                                         one_tensor, prefactor=alpha) 
                          ierr=exatns_tensor_destroy(one_tensor)
                        ! else
                        !   ierr=exatns_tensor_init(d%tens_exatensor, alpha)
                        endif
                              
                END SELECT
        
        if (.not. (present(a) .or. present(b) .or. present(c))) then
          ierr = initialize_tensor_t(d, alpha)
        endif

        if (allocated(arr_idx_A)) deallocate(arr_idx_A)
        if (allocated(arr_idx_B)) deallocate(arr_idx_B)
        if (allocated(arr_idx_C)) deallocate(arr_idx_C)
        if (allocated(arr_idx_D)) deallocate(arr_idx_D)

        t_einsum = ierr
        end function t_einsum

 
        ! copies data and destroys the source talsh tensor, unless staying in talsh
        function talsh_to_tensor_t(tens_talsh, out_tens, out_executor, destroy)
          implicit none
          type(talsh_tens_t), intent(inout)    :: tens_talsh
          type(tensor_t), intent(inout)    :: out_tens
          integer, intent(in), optional :: out_executor ! If not set value set by Cmake is used. 1=tapp, 2=talsh, 3=exatensor.
          integer, intent(in), optional:: destroy ! if set to 0 source tensor not destroyed, 1 it's destroyed (default) 

          integer             :: DataKind(1),nData
          integer(INTD)       :: ierr, talsh_to_tensor_t
          integer             :: destroy_
          type(tensor_t)    :: in_tens

          destroy_ = 1
          if (present(destroy)) then
            destroy_ = destroy
          endif

          ierr=talsh_tensor_data_kind(tens_talsh,nData,DataKind)
          if (DataKind(1) == R8) in_tens%datatype = 1
          if (DataKind(1) == C8) in_tens%datatype = 3

          in_tens%which_executor = 2 
          in_tens%tens_talsh = tens_talsh

          if (present(out_executor)) then
            ierr=convert(in_tens, out_tens, out_executor, destroy=destroy_)
          else
            ierr=convert(in_tens, out_tens, destroy=destroy_)
          endif

          talsh_to_tensor_t = ierr
        end function talsh_to_tensor_t

        ! copies data and destroys the source talsh tensor, unless staying in talsh
        function tensor_t_to_talsh(in_tens, tens_talsh, destroy)
          implicit none
          type(tensor_t), intent(inout)    :: in_tens
          type(talsh_tens_t), intent(inout)    :: tens_talsh
          integer, intent(in), optional:: destroy ! if set to 0 source tensor not destroyed, 1 it's destroyed (default) 

          integer(INTD)       :: ierr, tensor_t_to_talsh
          integer             :: destroy_
          type(tensor_t)    :: out_tens
          
          destroy_ = 1
          if (present(destroy)) then
            destroy_ = destroy
          endif
          ierr=convert(in_tens, out_tens, out_executor=2, destroy=destroy_)
          tens_talsh = out_tens%tens_talsh

          tensor_t_to_talsh = ierr
        end function tensor_t_to_talsh

        ! copies data and destroys the source tensor, unless staying in the same executor
        function convert(in_tens, out_tens, out_executor, destroy) ! If out_executor not set value set by Cmake is used.
          implicit none
          type(tensor_t), intent(inout)    :: in_tens
          type(tensor_t), intent(inout)    :: out_tens
          integer, intent(in), optional:: out_executor 
          integer, intent(in), optional:: destroy ! if set to 0 source tensor not destroyed, 1 it's destroyed (default) 
          integer(INTD)       :: ierr, convert
          integer             :: in_out_exec, destroy_

          integer(c_int) :: tapp_err, nmodes
          integer(c_int64_t), allocatable, target :: extents(:)
          integer(8) :: numel
          integer(INTD), allocatable :: out_dims(:)

          complex(c_double_complex), allocatable, target :: arr_c(:)
          real(c_double), allocatable, target :: arr_r(:)
          integer(INTD), allocatable, target :: lin_dim(:)
          complex(8), pointer, contiguous  :: lin_tens_c(:)
          real(8), pointer, contiguous  :: lin_tens_r(:)
          type(C_PTR) :: talshsptr
          complex(8), parameter :: ZERO_c=(0.D0,0.D0)
           
          if ( is_empty_tensor_t(in_tens) ) then
            write(*,'("  Error. Input tensor empty tensor in convert.")') 
            stop 98
          endif

          if (present(out_executor)) then
            out_tens%which_executor = out_executor
          endif
          destroy_ = 1
          if (present(destroy)) then
            destroy_ = destroy
          endif

          out_tens%datatype = in_tens%datatype
          
          in_out_exec = 100 * in_tens%which_executor + out_executor 
          SELECT CASE( in_out_exec) ! e.g. 201 means convert from exec 2 to 1 
            CASE(101, 202) 
              ierr = copy_tensor_t(in_tens, out_tens)
            CASE(102) 
              ierr = get_extents(in_tens, out_dims)
              ierr = construct_tensor_t(out_tens, out_tens%datatype, out_dims, executor=out_tens%which_executor)

              ierr = get_numel(in_tens, numel)
              allocate(lin_dim(1))
              lin_dim(1) = numel

              ierr=direct_data_pointer(out_tens,talshsptr)

              if (out_tens%datatype == 1 ) then
                call c_f_pointer(talshsptr,lin_tens_r,lin_dim)
                call dcopy(numel, in_tens%tens_tapp%data_real, 1, lin_tens_r, 1)
              elseif (out_tens%datatype == 3 ) then
                call c_f_pointer(talshsptr,lin_tens_c,lin_dim)
                call zcopy(numel, in_tens%tens_tapp%data_complex, 1, lin_tens_c, 1)
              endif

            CASE(201) 
              ierr = get_extents(in_tens, out_dims)
              ierr = construct_tensor_t(out_tens, out_tens%datatype, out_dims, executor=out_tens%which_executor)
              
              if (ierr.ne.0) then
                write(*,'("  Error. Converting tensor_t failed with error code: ",i11,".")') ierr
                stop 53
              end if

              ! copy data
              ierr = get_numel(in_tens, numel)
              allocate(lin_dim(1))
              lin_dim(1) = numel
  
              ierr=direct_data_pointer(in_tens,talshsptr)

              if (out_tens%datatype == 1 ) then
                call c_f_pointer(talshsptr,lin_tens_r,lin_dim)
                allocate(arr_r(numel)) 
                call dcopy(numel, lin_tens_r, 1, arr_r, 1)
                out_tens%tens_tapp%data_real = arr_r
              elseif (out_tens%datatype == 3 ) then
                call c_f_pointer(talshsptr,lin_tens_c,lin_dim)
                allocate(arr_c(numel)) 
                call zcopy(numel, lin_tens_c, 1, arr_c, 1)
                out_tens%tens_tapp%data_complex = arr_c
              endif

            CASE default 
              write(*,'("  Error. Conversion from executor ",i11," to ",i11," not yet implemented.")') &
                in_tens%which_executor, out_executor
              stop 55
          END SELECT

          if (destroy_ == 1) then
            ierr = destroy_tensor_t(in_tens)
          endif

          if (allocated(extents)) deallocate(extents)
          if (allocated(out_dims)) deallocate(out_dims)
          if (allocated(lin_dim)) deallocate(lin_dim)
          convert = ierr
        end function convert

        ! constructor for tensor_t
        function construct_tensor_t(tensor, datatype, extents, init_val, executor)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          integer, intent(in), optional    :: datatype ! 1=real double, 3 = complex double
          integer(INTD), intent(in), optional :: extents(:)
          ! integer(INTD), pointer, intent(in), optional :: extents(:)     
          integer(INTD), allocatable :: out_dims(:)
          complex(c_double_complex), intent(in), optional    :: init_val
          complex(c_double_complex)     :: init_val_
          integer, intent(in), optional    :: executor 

          integer(INTD)       :: ierr, construct_tensor_t
          integer(c_int)       :: tapp_err
          integer(c_int64_t), allocatable :: extents_(:)
          integer(INTD) :: extents_scalar(1)

          if ( .not. is_empty_tensor_t(tensor) ) then
            write(*,'("  Error. Tensor already constructed in construct_tensor_t.")') 
            stop 98
          endif

          if (present(executor)) tensor%which_executor = executor

          init_val_ = 0.0
          if (present(init_val)) init_val_ = init_val
          
          if (present(datatype)) tensor%datatype = datatype
          if (tensor%datatype == R8) then
            tensor%datatype = 1
          elseif (tensor%datatype == C8) then
            tensor%datatype = 3
          endif

          SELECT CASE(tensor%which_executor)
          CASE(1)
            if(.not. present(extents) .or. size(extents) == 0) then ! scalar
              tapp_err=construct_tapp_tensor_t(tensor%tens_tapp, tensor%datatype, init_val=init_val_)
            else
              allocate(extents_(size(extents)))
              extents_ = int(extents, kind=8)
              tapp_err=construct_tapp_tensor_t(tensor%tens_tapp, tensor%datatype, extents_, init_val_)
              deallocate(extents_)
            endif
            ierr = tapp_err
          CASE(2)
            extents_scalar = 1
            if (tensor%datatype == 1) then
              if(.not. present(extents) .or. size(extents) == 0) then ! scalar
                ierr=talsh_tensor_construct(tensor%tens_talsh, R8, extents_scalar(1:0), init_val=init_val_)
              else
                ierr=talsh_tensor_construct(tensor%tens_talsh, R8, extents, init_val=init_val_)
              endif
            elseif (tensor%datatype == 3) then
              if(.not. present(extents) .or. size(extents) == 0) then ! scalar
                ierr=talsh_tensor_construct(tensor%tens_talsh, C8, extents_scalar(1:0), init_val=init_val_)
              else
                ierr=talsh_tensor_construct(tensor%tens_talsh, C8, extents, init_val=init_val_)
              endif
            endif
          CASE default 
            write(*,'("  Error. construct_tensor_t for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
             
          END SELECT

          construct_tensor_t = ierr
        end function construct_tensor_t

        ! destructor for tensor_t. Safe to call for an empty tensor too.
        function destroy_tensor_t(tensor)
          implicit none
          type(tensor_t), intent(inout)    :: tensor

          integer(INTD)       :: ierr, destroy_tensor_t
          integer(c_int)       :: tapp_err

          if ( is_empty_tensor_t(tensor) ) then
            destroy_tensor_t = 0
            return
          endif

          SELECT CASE(tensor%which_executor)
          CASE(1)
             tapp_err=destroy_tapp_tensor_t(tensor%tens_tapp)
             ierr = tapp_err
          CASE(2)
             ierr=talsh_tensor_destruct(tensor%tens_talsh)
          CASE default 
            write(*,'("  Error. destroy_tensor_t for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
          END SELECT

          if ( .not. is_empty_tensor_t(tensor) ) then
            write(*,'("  Error. Failed sanity check in destroy_tensor_t.")') 
            stop 98
          endif

          destroy_tensor_t = ierr
        end function destroy_tensor_t


        ! get direct access to the tensor's elements
        function direct_data_pointer(tensor, ptr)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          type(C_PTR), intent(out) :: ptr

          integer(INTD)       :: ierr, direct_data_pointer   

          if ( is_empty_tensor_t(tensor) ) then
            write(*,'("  Error. Attempting to access data of an empty tensor in direct_data_pointer.")') 
            stop 98
          endif

          SELECT CASE(tensor%which_executor)
          CASE(1)
             ptr = get_data_ptr(tensor%tens_tapp)
             ierr = 0
          CASE(2)
              if (tensor%datatype == 1 ) then
                ierr=talsh_tensor_get_body_access(tensor%tens_talsh,ptr,R8,int(0,C_INT),DEV_HOST)
              elseif (tensor%datatype == 3 ) then
                ierr=talsh_tensor_get_body_access(tensor%tens_talsh,ptr,C8,int(0,C_INT),DEV_HOST)
              endif
          CASE default 
            write(*,'("  Error. direct_data_pointer for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
          END SELECT

          if ( .not. c_associated(ptr) ) then
            write(*,'("  Error. Null pointer exception in direct_data_pointer.")') 
            stop 99
          endif
          direct_data_pointer = ierr
        end function direct_data_pointer

        ! get nmodes (number of indices) of a tensor. Returns 0 for a scalar.
        function get_nmodes(tensor, nmodes)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          integer(c_int), intent(out) :: nmodes

          integer(INTD)       :: ierr, get_nmodes   

          if ( is_empty_tensor_t(tensor) ) then
            nmodes = 0
            get_nmodes = 0
            return
            ! write(*,'("  error. attempting to access data of an empty tensor in get_nmodes.")') 
            ! stop 98
          endif

          SELECT CASE(tensor%which_executor)
          CASE(1)
             ierr =  get_nmodes_tapp_tensor_t(tensor%tens_tapp, nmodes)
          CASE(2)
             nmodes = talsh_tensor_rank(tensor%tens_talsh)
             ierr = 0 
          CASE default 
            write(*,'("  Error. get_nmodes for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
          END SELECT

          get_nmodes = ierr
        end function get_nmodes

        ! get extents (dimensions) of a tensor. In case of scalar, it returns empty array
        function get_extents(tensor, extents)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          integer(INTD), allocatable, intent(inout) :: extents(:)
          integer(c_int64_t), allocatable :: extents_(:)
          integer(c_int)          :: nmodes

          integer(INTD)       :: ierr, get_extents 

          if ( is_empty_tensor_t(tensor) ) then
            write(*,'("  Error. Attempting to access an empty tensor in get_extents.")') 
            stop 98
          endif

          ierr =  get_nmodes(tensor, nmodes)
          if (.not. allocated(extents)) then
            allocate(extents(nmodes))
          endif
          if (.not. size(extents) == nmodes) then 
            write(*,'("  Error. An array of wrong size passed to get_extents.")') 
            stop 58
          endif

          select case(tensor%which_executor)
          case(1)
                allocate(extents_(nmodes))
                ierr = get_extents_tapp_tensor_t(tensor%tens_tapp, extents_)

                if (size(extents_) .ne. 0) then 
                  extents = int(extents_, kind=INTD)
                endif
                deallocate(extents_) 
          case(2)
             ierr = talsh_tensor_dimensions(tensor%tens_talsh,nmodes,extents)
          case default 
            write(*,'("  Error. get_extents for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
          end select
          
          if ( .not. size(extents) == 0 .and. is_scalar_tensor_t(tensor)) then
            write(*,'("  Error. Failed sanity check in get_extents.")') 
            stop 98
          endif

          get_extents = ierr
        end function get_extents

        ! get the total number of elements (scalars) which the tensor occupies in memory
        function get_numel(tensor, numel)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          integer(c_int)       :: i, nmodes
          integer(INTD), allocatable :: extents(:)

          integer(8)       :: numel   
          integer(INTD)       :: ierr, get_numel   
          
          if (is_empty_tensor_t(tensor)) then
            numel = 0
            ierr = 0
          else
            ierr = get_nmodes(tensor, nmodes)
            numel = 1
            if (nmodes .ne. 0) then !not scalar
              ierr = get_extents(tensor, extents)
              do i=1,nmodes ! problem: assumes a dense tensor
                numel = numel * extents(i)
              enddo
            endif
          endif 

          if (allocated(extents)) deallocate(extents)
          get_numel = ierr
        end function get_numel

        ! check if a tensor is empty (constructor not used)
        function is_empty_tensor_t(tensor)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          logical             :: ret, is_empty_tensor_t   
          
          SELECT CASE(tensor%which_executor)
          CASE(1)
             ret=is_empty_tapp_tensor_t(tensor%tens_tapp)
          CASE(2)
             ret=talsh_tensor_is_empty(tensor%tens_talsh)
          CASE default 
            write(*,'("  Error.  is_empty_tensor_t for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
          END SELECT
 
          is_empty_tensor_t = ret
        end function is_empty_tensor_t

        ! check if tensor_t is a scalar (0 indices, i.e. can be different to a vector with 1 index and 1 element)
        function is_scalar_tensor_t(tensor)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          integer(INTD)       :: ierr   
          integer(c_int) ::  nmodes
          logical             :: ret, is_scalar_tensor_t  
          
          if (is_empty_tensor_t(tensor)) then
            ret = .FALSE.
          else
            ierr = get_nmodes(tensor, nmodes)
            if (nmodes == 0) then !scalar
              ret = .TRUE.
            else
              ret = .FALSE.
            endif
          endif 
          is_scalar_tensor_t = ret
        end function is_scalar_tensor_t

        ! initialize tensor_t to a given scalar value
        function initialize_tensor_t(tensor, init_val)
          implicit none
          type(tensor_t), intent(inout)    :: tensor
          complex(c_double_complex), intent(in), optional    :: init_val
          complex(c_double_complex)     :: init_val_

          integer(INTD), allocatable, target :: lin_dim(:)
          complex(8), pointer, contiguous  :: lin_tens_c(:)
          real(8), pointer, contiguous  :: lin_tens_r(:)
          integer(INTD), allocatable :: extents(:)
          integer(8) :: numel

          integer(INTD)       :: ierr, initialize_tensor_t
          type(C_PTR) :: ptr
          integer(c_int)       :: tapp_err, nmodes
          
          init_val_ = 0.0
          if (present(init_val)) init_val_ = init_val

          ierr = get_numel(tensor, numel)

          allocate(lin_dim(1))
          lin_dim(1) = numel
          ierr = direct_data_pointer(tensor, ptr)
          SELECT CASE(tensor%which_executor)
          CASE(1)
              if (tensor%datatype == 1 ) then
                call c_f_pointer(ptr,lin_tens_r,lin_dim)
                lin_tens_r = init_val_%re
              elseif (tensor%datatype == 3 ) then
                call c_f_pointer(ptr,lin_tens_c,lin_dim)
                lin_tens_c = init_val_
              endif
          CASE(2)
            ierr=talsh_tensor_init(tensor%tens_talsh,init_val_)
          CASE default 
            write(*,'("  Error. initialize_tensor_t for executor ",i11," not yet implemented.")') &
              tensor%which_executor
            stop 55
          END SELECT

          if (allocated(lin_dim)) deallocate(lin_dim)
          initialize_tensor_t = ierr
        end function initialize_tensor_t

        ! copies a tensor (copy constructor, deep copy)
        function copy_tensor_t(src_tens, dest_tens)
          implicit none
          type(tensor_t), intent(inout)    :: src_tens
          type(tensor_t), intent(inout)    :: dest_tens

          integer(INTD), allocatable, target :: lin_dim(:)
          complex(8), pointer, contiguous  :: lin_src_tens_c(:)
          complex(8), pointer, contiguous  :: lin_dest_tens_c(:)
          real(8), pointer, contiguous  :: lin_src_tens_r(:)
          real(8), pointer, contiguous  :: lin_dest_tens_r(:)
          integer(INTD), allocatable :: out_dims(:)
          integer(8) :: numel
          type(C_PTR) :: src_ptr, dest_ptr
          integer(INTD)       :: ierr, copy_tensor_t

          if ( is_empty_tensor_t(src_tens) ) then
            write(*,'("  Error. Source tensor empty tensor in copy_tensor_t.")') 
            stop 98
          endif

          ierr = get_extents(src_tens, out_dims)
          ierr = get_numel(src_tens, numel)
          
          allocate(lin_dim(1))
          lin_dim(1) = numel

          ! if the destination tensor already exists, destroy it
          ierr = destroy_tensor_t(dest_tens)

          ierr = construct_tensor_t(dest_tens, src_tens%datatype, out_dims, executor=src_tens%which_executor)

          ! copy data
          ierr = direct_data_pointer(src_tens, src_ptr)
          ierr = direct_data_pointer(dest_tens, dest_ptr)

          if (src_tens%datatype == 1 ) then
            call c_f_pointer(src_ptr,lin_src_tens_r,lin_dim)
            call c_f_pointer(dest_ptr,lin_dest_tens_r,lin_dim)
            call dcopy(numel, lin_src_tens_r, 1, lin_dest_tens_r, 1)
          elseif (src_tens%datatype == 3 ) then
            call c_f_pointer(src_ptr,lin_src_tens_c,lin_dim)
            call c_f_pointer(dest_ptr,lin_dest_tens_c,lin_dim)
            call zcopy(numel, lin_src_tens_c, 1, lin_dest_tens_c, 1)
          endif

          if (allocated(out_dims)) deallocate(out_dims)
          copy_tensor_t = ierr
        end function copy_tensor_t
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
