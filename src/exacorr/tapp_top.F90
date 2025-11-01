module tapp_top
! This module is a high-level generalization of TAPP standard tensor interface. See github.com/TAPPorg
! Written by Jan Brandejs, Dec 2024

 use talsh
 use tensor_algebra
 use, intrinsic:: ISO_C_BINDING
 use tapp
    
 implicit none

#ifndef DEFAULT_EXECUTOR_WITHIN_TAPP
#  define DEFAULT_EXECUTOR_WITHIN_TAPP 1
#endif
integer, parameter, public :: TAPP_DEFAULT_EXECUTOR = DEFAULT_EXECUTOR_WITHIN_TAPP ! 1 bruteforece, 2 tblis, 3 cyclops


!DERIVED TYPES:
        public :: TAPP_F64, TAPP_C64
        enum, bind(c)
                ! we follow TAPP_datatype enumeration
                enumerator :: TAPP_datatype = 0 ! name of the enum
                enumerator :: TAPP_F64 = 1 ! double precision real
                enumerator :: TAPP_C64 = 3 ! double precision complex
        end enum
 
  !tapp_top's own tensor type (not part of low level TAPP)
        type, public  :: tapp_tensor_t
          !library_tensor_instance_pointer
          type(C_PTR) :: tensor_info  ! TAPP_tensor_info       
          logical     :: constructed = .false.

          ! TAPP reference implementation does not store data
          complex(c_double_complex), allocatable :: data_complex(:) !Do not access in your code outside this file (and exacorr_tensor_interface)!
          real(c_double),            allocatable :: data_real(:)    !Do not access in your code outside this file (and exacorr_tensor_interface)!
          ! integer(c_int64_t),        

          integer(kind(TAPP_datatype)) :: datatype
        end type tapp_tensor_t

 contains

 function tapp_einsum(D, arr_idx_D, op_D, A, arr_idx_A, op_A, B, arr_idx_B, op_B, C, arr_idx_C, op_C, &
                                       alpha, beta)
    implicit none
    type(tapp_tensor_t), target, intent(inout)    :: D
    type(tapp_tensor_t), target, intent(inout), optional    :: A, B, C
    integer(c_int64_t), allocatable, target, intent(in) :: arr_idx_D(:)
    integer(c_int64_t), allocatable, target, intent(in), optional :: arr_idx_A(:), &
      arr_idx_B(:), arr_idx_C(:)
    integer(c_int), intent(in), optional :: op_A, op_B, op_C, op_D ! 0 .. identitity, 1 .. conjugate
    complex(c_double_complex), intent(in), optional :: alpha, beta
    integer(c_int) :: tapp_einsum
    
    ! TAPP_create_tensor_product 
    integer(c_int) :: tapp_err, prec !  TAPP_prectype -1 .. default auto
    type(C_PTR) :: handle, plan ! TAPP_handle, plan=TAPP_tensor_product
    type(C_PTR) :: info_A, info_B, info_C, info_D ! TAPP_tensor_info
    type(C_PTR):: idx_A, idx_B, idx_C, idx_D
    integer(c_int) :: op_A_, op_B_, op_C_, op_D_ ! 0 .. identitity, 1 .. conjugate

    ! TAPP_execute_product
    type(C_PTR) :: exec, status ! TAPP_executor, TAPP_status
    type(C_PTR):: A_, B_, C_, D_
    complex(c_double_complex), target :: alpha_, beta_
    complex(c_double), target :: alpha_r, beta_r
    integer(c_int), target  :: executor_id ! 1 = bruteforce, 2 = tblis, 12 = tblis + bruteforce check
    type(C_PTR):: alpha_ptr, beta_ptr
    
    ! TAPP_create_tensor_info
    integer(c_int) :: datatype, nmodes ! datatype = TAPP_datatype
    type(C_PTR):: extents, strides
    real(c_double), target :: tmpr           ! just a temporary to make a scalar tensor
    complex(c_double_complex), target :: tmpc ! just a temporary to make a scalar tensor
    
    integer :: i

    alpha_ = 1.0 ! default values of alpha and beta are 1.0
    beta_ = 1.0
    if (present(alpha)) alpha_ = alpha
    if (present(beta)) beta_ = beta
    ! to.do throw error when alpha or beta complex but the tensors are real?

    prec = -1 ! default prec: computational precision is determined as *at least* the lowest precision 
              ! of the storage data types of all input and output operands. See TAPP_prectype for further options

    handle = c_null_ptr ! tapp library handle will need to be stored somewhere. to.do

    op_A_ = 0 
    op_B_ = 0
    op_C_ = 0 
    op_D_ = 0
    if (present(op_A)) op_A_ = op_A
    if (present(op_B)) op_B_ = op_B
    if (present(op_C)) op_C_ = op_C
    if (present(op_D)) op_D_ = op_D

    if (TAPP_DEFAULT_EXECUTOR.ne.3) then
      D_ = get_data_ptr(D)
    else
      D_ = c_null_ptr
    endif 
    info_D = D%tensor_info 
    idx_D = c_loc(arr_idx_D(1)) 
    if (present(A)) then
      if (TAPP_DEFAULT_EXECUTOR.ne.3) then
        A_ = get_data_ptr(A)
      else
        A_ = c_null_ptr
      endif 
      info_A = A%tensor_info 
      idx_A = c_loc(arr_idx_A(1)) 
    endif
    if (present(B)) then
      if (TAPP_DEFAULT_EXECUTOR.ne.3) then
        B_ = get_data_ptr(B)
      else
        B_ = c_null_ptr
      endif 
      info_B = B%tensor_info 
      idx_B = c_loc(arr_idx_B(1)) 
    endif
    if (present(C)) then
      if (TAPP_DEFAULT_EXECUTOR.ne.3) then
        C_ = get_data_ptr(C)
      else
        C_ = c_null_ptr
      endif 
      info_C = C%tensor_info 
      idx_C = c_loc(arr_idx_C(1)) 
    endif

    if (present(A) .and. present(A) .and. present(C)) then ! D = alpha A * B + beta C
      ! all set
    elseif (present(A) .and. present(B)) then ! this defaults to D = alpha A * B + beta D
      C_ = D_
      info_C = info_D
      idx_C = idx_D
      op_C_ = op_D_
    elseif (present(A)) then ! this defaults to D = alpha A + beta D
      C_ = D_
      info_C = info_D
      idx_C = idx_D
      op_C_ = op_D_
      
      ! B has to be created as a scalar. 
      ! to.do : check that contracting with scalar tensor stays a part of the standard
      ! to.do : make this code part of create_tensor routine
      nmodes = 0
      extents = c_null_ptr
      strides = c_null_ptr
      if (TAPP_DEFAULT_EXECUTOR.ne.3) then
        tapp_err = TAPP_create_tensor_info(info_B, A%datatype, nmodes, extents, strides)
        if (tapp_err.ne.0) then
          write(*,'("  Error. Creating temporary TAPP tensor info in tapp_einsum failed with error code: ",i11,".")') &
          tapp_err
          stop 40
        end if
      end if
      if (A%datatype == 1 ) then
        tmpr   = 1.0
        tmp = c_loc(tmpr) 
      elseif (A%datatype == 3 ) then
        tmpc   = 1.0
        tmp = c_loc(tmpc) 
      endif
      if (TAPP_DEFAULT_EXECUTOR == 3) then
        tapp_err = TAPP_create_tensor(info_B, A%datatype, nmodes, extents, strides, tmp)
        if (tapp_err.ne.0) then
          write(*,'("  Error. Creating temporary TAPP tensor info in tapp_einsum failed with error code: ",i11,".")') &
            tapp_err
          stop 40
        end if
      endif
      idx_B = c_null_ptr 
      
    else ! this is interpreted elementwise as D = alpha and requires alpha defined
      if (.not. present(alpha)) then
        write(*,'("  Error. At least a scalar alpha has to be given in tapp_einsum when A,B,C not specified.")')
        stop 48
      end if

      if (TAPP_DEFAULT_EXECUTOR.ne.3 .and. D%datatype == 1 .and. allocated(D%data_real)) then
        do i=1,size(D%data_real)
          D%data_real(i) = alpha_%re
        enddo
        tapp_einsum = tapp_err
        return
      elseif (TAPP_DEFAULT_EXECUTOR == 3 .and. D%datatype == 1) then
         alpha_r=alpha_%re
         alpha_ptr = c_loc(alpha_r)
         tapp_err = TAPP_initialize_tensor(info_D, alpha_ptr)
         tapp_einsum = tapp_err
         return
      elseif (TAPP_DEFAULT_EXECUTOR.ne.3 .and. D%datatype == 3 .and. allocated(D%data_complex)) then
        do i=1,size(D%data_complex)
          D%data_complex(i) = alpha_
        enddo
        tapp_einsum = tapp_err
        return
      elseif (TAPP_DEFAULT_EXECUTOR == 3 .and. D%datatype == 3) then
         alpha_ptr = c_loc(alpha_)
         tapp_err = TAPP_initialize_tensor(info_D, alpha_ptr)
         tapp_einsum = tapp_err
         return
      else
        ! to.do create a new tensor
        write(*,'("  Error. The case D = alpha without D allocated is not implemented yet in tapp_einsum.")')
        stop 47
      endif
    endif

    ! prepare a contraction plan
    tapp_err = TAPP_create_tensor_product(plan, handle, op_A_, info_A, idx_A, op_B_, info_B, idx_B, &
      op_C_, info_C, idx_C, op_D_, info_D, idx_D, prec)

    if (tapp_err.ne.0) then
      write(*,'("  Error. Creating TAPP contraction plan failed with error code: ",i11,".")') tapp_err
      stop 41
    end if

    if (A%datatype == 1 .and. (.not. present(B) .or. B%datatype == 1) ) then ! A is present at this point
        alpha_r=alpha_%re
        beta_r=beta_%re ! todo separate condition for beta?
        alpha_ptr = c_loc(alpha_r) 
        beta_ptr = c_loc(beta_r) 
    else 
        alpha_ptr = c_loc(alpha_) 
        beta_ptr = c_loc(beta_) 
    endif
    
    executor_id = TAPP_DEFAULT_EXECUTOR ! 1 = bruteforce, 2 = tblis, 12 = tblis + bruteforce check
    exec = c_loc(executor_id)

    ! execute the contraction plan
    tapp_err = TAPP_execute_product(plan, exec, status, alpha_ptr, A_, B_, beta_ptr, C_, D_) ! diff between plan par here and elsewhere?

    if (tapp_err.ne.0) then
      write(*,'("  Error. TAPP tensor product failed with error code: ",i11,".")') tapp_err
      call quit(" failed product")
      stop 42
    end if

    tapp_err = TAPP_destroy_tensor_product(plan)
    tapp_einsum = tapp_err
 end function tapp_einsum


 ! parses the talsh format einsum string into 
 function parse_t_einsum(expression, arr_idx_A, op_A, arr_idx_B, op_B, arr_idx_C, op_C, arr_idx_D, op_D)
  implicit none
  character(len=*), intent(in) :: expression
  integer(c_int64_t), allocatable, target, intent(inout) :: arr_idx_A(:), &
   arr_idx_B(:), arr_idx_C(:), arr_idx_D(:)
  integer(c_int), intent(out) :: op_A, op_B, op_C, op_D ! 0 .. identitity, 1 .. conjugate

  character(len=50), allocatable :: tensor_names(:), algebraic_operations(:)
  character(len=10), allocatable :: indices(:,:)
  integer, allocatable :: index_count(:), tensor_conjugate(:)
  character(len=50), allocatable :: distinct_indices(:)
  integer(c_int64_t), allocatable :: map_to(:)
  
  integer(c_int64_t) :: i8
  integer :: ierr, i, j, k
  integer :: parse_t_einsum

  ierr = pre_parse_t_einsum(expression, tensor_names, indices, index_count, tensor_conjugate, algebraic_operations)
  do i = 1, size(indices, 1)
    do j = 1, index_count(i)
      if (allocated(distinct_indices)) then
        do k = 1, size(distinct_indices)
          if (trim(indices(i,j)) == trim(distinct_indices(k))) then
            exit
          end if
          if (k == size(distinct_indices)) then
            call add_str(distinct_indices)
            distinct_indices(k+1) = trim(indices(i,j))
          end if 
        end do
      else
        call add_str(distinct_indices)
        distinct_indices(1) = trim(indices(i,j))
      end if
    end do
  end do
  
  call bubble_sort(distinct_indices) 
  allocate(map_to(size(distinct_indices)))
  do i = 1, size(distinct_indices) 
     map_to(i) = 96+i ! 97 is the ichar value of 'a'
  end do
  
  ! now create idx_A,B, .. index arrays
  if (size(indices, 1) > 0) then
    allocate(arr_idx_D(index_count(1)))
    do j = 1, index_count(1)
      arr_idx_D(j) = map_str_to_i8(indices(1,j), distinct_indices, map_to)
    end do
    op_D = tensor_conjugate(1)
  end if 
  if (size(indices, 1) > 1) then
    allocate(arr_idx_A(index_count(2)))
    do j = 1, index_count(2)
      arr_idx_A(j) = map_str_to_i8(indices(2,j), distinct_indices, map_to)
    end do
    op_A = tensor_conjugate(2)
  end if 
  if (size(indices, 1) > 2) then
    allocate(arr_idx_B(index_count(3)))
    do j = 1, index_count(3)
      arr_idx_B(j) = map_str_to_i8(indices(3,j), distinct_indices, map_to)
    end do
    op_B = tensor_conjugate(3)
  end if 
  if (size(indices, 1) > 3) then
    allocate(arr_idx_C(index_count(4)))
    do j = 1, index_count(4)
      arr_idx_C(j) = map_str_to_i8(indices(4,j), distinct_indices, map_to)
    end do
    op_C = tensor_conjugate(4)
    if (trim(algebraic_operations(3)) /= "+" .and. trim(algebraic_operations(3)) /= "-") then
      write(*,'("  Error after parsing t_einsum expression: tesor C can follow only after + or -.")')
      stop 49
    end if
  end if 
  
  deallocate(tensor_names, indices, index_count, distinct_indices, map_to, tensor_conjugate, algebraic_operations)
  parse_t_einsum = ierr
 end function parse_t_einsum 

 function pre_parse_t_einsum(expression, tensor_names, indices, index_count, tensor_conjugate, algebraic_operations)
  implicit none

  character(len=*), intent(in) :: expression
  character(len=50), allocatable, intent(out) :: tensor_names(:)
  character(len=10), allocatable, intent(out) :: indices(:,:)
  integer, allocatable, intent(out) :: index_count(:)
  integer, allocatable, intent(out) :: tensor_conjugate(:)
  character(len=50), allocatable, intent(out) :: algebraic_operations(:)
  integer :: pre_parse_t_einsum

  ! Internal variables
  character(len=50) :: current_tensor
  character(len=10) :: current_index
  character(len=50) :: current_operator
  integer :: i, j, k, l, parenthesis_level, ierr, ntens, current_conjugate

  ! Initialize
  ierr = -1
  ! tensor_names = []
  ! indices = reshape([], [0, 0])

  ! Simple parsing (assuming well-formed expression)
  parenthesis_level = 0
  j = 1
  k = 1
  ntens = 0
  current_tensor = ""
  current_index = ""
  current_conjugate = 0
  current_operator = ""
  do i = 1, len_trim(expression)
    ! Extract tensor names
    if (expression(i:i) == '(') then
      parenthesis_level = parenthesis_level + 1
      if (parenthesis_level == 1) then
        ntens = ntens + 1
        call add_str(tensor_names)
        tensor_names(ntens) = trim(current_tensor)
        
        call add_str_row(indices)
        call push_back_int(index_count, 0)
        call push_back_int(tensor_conjugate, current_conjugate)
        call add_str(algebraic_operations)
        algebraic_operations(ntens) = ""
        if (ntens > 1) then ! otherwise current_operator should be empty
          algebraic_operations(ntens-1) = trim(current_operator)
        endif
        
        current_tensor = ""
        current_conjugate = 0
        current_operator = ""
        k = 1
        l = 1
      end if
    elseif (expression(i:i) == ')' .or. expression(i:i) == ',') then
      if (parenthesis_level == 1 .and. j>1) then
        index_count(ntens) = index_count(ntens) + 1
        if (size(indices, 2) < index_count(ntens)) then
          call add_str_col(indices)
        endif
        indices(ntens, index_count(ntens)) = trim(current_index)
        current_index = ""
        j = 1
      end if
      if (expression(i:i) == ')') then
        parenthesis_level = parenthesis_level - 1
      endif
    elseif (expression(i:i) == '+' .and. k > 1 .and. parenthesis_level == 0) then
        current_conjugate = 1
    elseif (expression(i:i) == '=' .or. &
            expression(i:i) == '+' .or. &
            expression(i:i) == '-' .or. &
            expression(i:i) == '*') then
      if (k == 1 .and. parenthesis_level == 0) then
        current_operator(l:l) = expression(i:i)
        l = l + 1
      endif 
    elseif (expression(i:i) /= ' ') then
      if (parenthesis_level == 0) then
        current_tensor(k:k) = expression(i:i)
        k = k + 1
      else
        current_index(j:j) = expression(i:i)
        j = j + 1
      end if
    end if
  end do

  ! current_operator should be empty at this point

  ierr = 0
  pre_parse_t_einsum = ierr
 end function pre_parse_t_einsum 

 subroutine add_str(list)
  IMPLICIT NONE
  integer :: i, isize
  character(len=50), allocatable, intent(inout) :: list(:)
  character(len=50), allocatable :: clist(:)
  
  if(allocated(list)) then
      isize = size(list)
      allocate(clist(isize+1))
      do i=1,isize          
      clist(i) = list(i)
      end do
      clist(isize+1) = ""
   
      deallocate(list)
      call move_alloc(clist, list)
  
  else
      allocate(list(1))
      list(1) = ""
  end if
 end subroutine add_str

 subroutine push_back_int(list, element)
  IMPLICIT NONE
  integer :: i, isize
  integer, intent(in) :: element
  integer, allocatable, intent(inout) :: list(:)
  integer, allocatable :: clist(:)
  
  if(allocated(list)) then
      isize = size(list)
      allocate(clist(isize+1))
      do i=1,isize          
      clist(i) = list(i)
      end do
      clist(isize+1) = element
   
      deallocate(list)
        call move_alloc(clist, list)
  
  else
      allocate(list(1))
      list(1) = element
  end if
 end subroutine push_back_int

 subroutine add_str_row(list)
  IMPLICIT NONE
  integer :: i,j, isize
  character(len=10), allocatable, intent(inout) :: list(:,:)
  character(len=10), allocatable :: clist(:,:)
  
  if(allocated(list)) then
      isize = size(list, 1)
      allocate(clist(isize+1, size(list, 2)))
      do i=1,isize          
        do j=1,size(list, 2)          
          clist(i,j) = list(i,j)
        end do
      end do
      do j=1,size(list, 2)          
        clist(isize+1, j) = ""
      end do
   
      deallocate(list)
      call move_alloc(clist, list)
  
  else
      allocate(list(1,16))
      do j=1,size(list, 2)          
        list(1, j) = ''
      end do
  end if
 end subroutine add_str_row

 subroutine add_str_col(list)
  IMPLICIT NONE
  integer :: i,j, jsize
  character(len=10), allocatable, intent(inout) :: list(:,:)
  character(len=10), allocatable :: clist(:,:)
  
  if(allocated(list)) then
      jsize = size(list, 2)
      allocate(clist(size(list, 1), jsize+1))
      do i=1,size(list, 1)          
        do j=1,jsize          
          clist(i,j) = list(i,j)
        end do
      end do
      do i=1,size(list, 1)          
        clist(i, jsize+1) = ""
      end do
   
      deallocate(list)
      call move_alloc(clist, list)
  
  else
      allocate(list(1,16))
      do j=1,size(list, 2)          
        list(1, j) = ""
      end do
  end if
 end subroutine add_str_col

 subroutine bubble_sort(arr)
    implicit none
    character(len=50), allocatable, intent(inout) :: arr(:)
    integer :: n
    integer :: i, j
    character(len=50) :: temp
    n = size(arr)
    do i = 1, n - 1
      do j = i + 1, n
        if (arr(i) > arr(j)) then
          temp = arr(i)
          arr(i) = arr(j)
          arr(j) = temp
        end if
      end do
    end do
 end subroutine bubble_sort

 function map_str_to_i8(element, list, map_to)
    implicit none
    character(len=10), intent(in) :: element
    character(len=50), allocatable, intent(in) :: list(:)
    integer(c_int64_t), allocatable, intent(in) ::  map_to(:)
    integer :: i
    integer(c_int64_t) :: i8
    integer(c_int64_t) :: map_str_to_i8
    
    i8 = -1
    do i = 1, size(list)
      if (trim(element) == trim(list(i))) then
        i8 = map_to(i)
        exit
      end if
    end do
    map_str_to_i8 = i8
 end function map_str_to_i8


 function broadcast_tapp_tensor_t_body_data(tensor, custom_ptr)
    implicit none
    type(tapp_tensor_t), target, intent(inout)    :: tensor
    type(C_PTR), intent(inout), optional :: custom_ptr
    integer(c_int)       :: broadcast_tapp_tensor_t_body_data, tapp_err

    if(present(custom_ptr)) then
      tapp_err = TAPP_broadcast_tensor_body_data(tensor%tensor_info, custom_ptr)
    else
      ! allocate and initialize data
      if (tensor%datatype == 1 ) then
        if (.not. allocated(tensor%data_real)) then
             write(*,'("  Error. Tensor nor allocated in broadcast_tapp_tensor_t_body_data.")')
             stop 62
        endif
        tapp_err = TAPP_broadcast_tensor_body_data(tensor%tensor_info, c_loc(tensor%data_real(1)))

      elseif (tensor%datatype == 3 ) then
        if (.not. allocated(tensor%data_complex)) then
             write(*,'("  Error. Tensor nor allocated in broadcast_tapp_tensor_t_body_data.")')
             stop 62
        endif
        tapp_err = TAPP_broadcast_tensor_body_data(tensor%tensor_info, c_loc(tensor%data_complex(1)))
      endif
    endif

    broadcast_tapp_tensor_t_body_data = tapp_err
 end function broadcast_tapp_tensor_t_body_data


 function gather_tapp_tensor_t_body_data(tensor)
    implicit none
    type(tapp_tensor_t), target, intent(inout)    :: tensor
    integer(8) :: i, numel

    complex(c_double_complex), allocatable, target :: arr_c(:)
    real(c_double), allocatable, target :: arr_r(:)
    integer(c_int64_t), allocatable, target  :: extents(:)
    integer(c_int)       :: gather_tapp_tensor_t_body_data, tapp_err, nmodes


    tapp_err = get_nmodes_tapp_tensor_t(tensor, nmodes)
    tapp_err = get_extents_tapp_tensor_t(tensor, extents)
    numel = 1
    do i=1,nmodes
      numel = numel * extents(i)
    enddo
    ! allocate and initialize data
    if (tensor%datatype == 1 ) then
      if (.not. allocated(tensor%data_real)) then
        allocate(arr_r(numel))
        tensor%data_real = arr_r
      endif
      tapp_err = TAPP_gather_tensor_body_data(tensor%tensor_info, c_loc(tensor%data_real(1)))

    elseif (tensor%datatype == 3 ) then
      if (.not. allocated(tensor%data_complex)) then
        allocate(arr_c(numel))
        tensor%data_complex = arr_c
      endif
      tapp_err = TAPP_gather_tensor_body_data(tensor%tensor_info, c_loc(tensor%data_complex(1)))
    endif

    if (allocated(extents)) deallocate(extents)

    gather_tapp_tensor_t_body_data = tapp_err
 end function gather_tapp_tensor_t_body_data

 function get_data_ptr(tensor)
    implicit none
    type(tapp_tensor_t), target, intent(inout)    :: tensor
    type(C_PTR):: get_data_ptr
    
    get_data_ptr = c_null_ptr

    if (TAPP_DEFAULT_EXECUTOR.ne.3) then
      ! Warning: creating FULL tensor array locally (likely just to get a data pointer). Not memory-scalable!
      tapp_err = gather_tapp_tensor_t_body_data(tensor)
    end if

    if (tensor%datatype == 1 ) then
      get_data_ptr = c_loc(tensor%data_real(1)) 
    elseif (tensor%datatype == 3 ) then
      get_data_ptr = c_loc(tensor%data_complex(1)) 
    endif
 end function get_data_ptr

 ! splits talsh-formatter einsum expression into contraction part and addition part for D = AB + C
 function split_DABC_expression_str(expression, op_C, contraction_str, addition_str)
         character(len = *), intent(in)   :: expression
         integer(c_int), intent(in) :: op_C ! 0 .. identitity, 1 .. element-wise complex conjugate
         character(len = len_trim(expression)), intent(out)  :: contraction_str, addition_str
         integer(INTD)  :: split_DABC_expression_str
         integer :: i, pos, curr, curr2, eq_found
         
         pos = -1
         do i=1,len_trim(expression)
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
           endif
         enddo
         if (pos == -1) then ! assumes that + is used for addition
           write(*,'("  Error. + not found in expression D = A*B+C.")')
           stop 45
         end if
         contraction_str = ""
         eq_found = 0
         do i=1,(pos-1)
           contraction_str(i:i)=expression(i:i)
           if (eq_found == 0) addition_str(i:i)=expression(i:i)
           if (expression(i:i) == '=') eq_found = i
         enddo
         addition_str = ""
         do i=1,(len_trim(expression)-pos)
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

 ! constructor for tapp_tensor_t
 function construct_tapp_tensor_t(tensor, datatype, extents, init_val)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor
   integer(c_int), intent(in), optional    :: datatype ! 1=real double, 3 = complex double
   integer(c_int64_t), allocatable, target, intent(in), optional :: extents(:)
   complex(c_double_complex), intent(in), optional    :: init_val
   complex(c_double_complex), target     :: init_val_
   real(c_double), target     :: init_val_r

   integer(c_int)       :: tapp_err, nmodes, datatype_,  construct_tapp_tensor_t
   integer(c_int64_t), allocatable, target :: strides(:)
   integer(8) :: i, k, stride, numel
   type(c_ptr)     :: tmp

   complex(c_double_complex), allocatable, target :: arr_c(:)
   real(c_double), allocatable, target :: arr_r(:)

   if (present(datatype)) tensor%datatype = datatype
   init_val_ = 0.0
   if (present(init_val)) init_val_ = init_val

   if (tensor%datatype == 1) then
     init_val_r = init_val_%re
     tmp = c_loc(init_val_r)
   elseif (tensor%datatype == 3) then
     tmp = c_loc(init_val_)
   endif

   if (.not. is_empty_tapp_tensor_t(tensor)) then
     write(*,'("  Error. tensor already constructed in construct_tapp_tensor_t.")')
     call quit(" Error. tensor already constructed in construct_tapp_tensor_t.")
     stop 30
   endif

   tensor%constructed = .TRUE.

   nmodes = 0
   if (present(extents) .and. size(extents) .gt. 0) then
     nmodes = size(extents)
   else ! scalar
     if (TAPP_DEFAULT_EXECUTOR.ne.3) then
      tapp_err = TAPP_create_tensor_info(tensor%tensor_info, tensor%datatype, nmodes, &
                 c_null_ptr, c_null_ptr)
      if (tensor%datatype == 1) then
       allocate(arr_r(1))
       arr_r(1) = init_val_%re 
       tensor%data_real = arr_r
      elseif (tensor%datatype == 3) then
       allocate(arr_c(1))
       arr_c(1) = init_val_ 
       tensor%data_complex = arr_c
      endif 

     else
      tapp_err = TAPP_create_tensor(tensor%tensor_info, tensor%datatype, nmodes, &
                 c_null_ptr, c_null_ptr, tmp)
     endif 
     
     construct_tapp_tensor_t = tapp_err
     return
   endif

   allocate(strides(nmodes))
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
   
   datatype_ = tensor%datatype
   
   if (TAPP_DEFAULT_EXECUTOR.ne.3) then
    tapp_err = TAPP_create_tensor_info(tensor%tensor_info, datatype_, nmodes, &
                 c_loc(extents(1)), c_loc(strides(1)))
   else
    tapp_err = TAPP_create_tensor(tensor%tensor_info, datatype_, nmodes, &
                 c_loc(extents(1)), c_loc(strides(1)), tmp)
   endif
   if (tapp_err.ne.0) then
     write(*,'("  Error. Creating TAPP tensor info failed with error code: ",i11,".")') tapp_err
     stop 40
   end if

   if (TAPP_DEFAULT_EXECUTOR.ne.3) then
   ! allocate and initialize data
    if (tensor%datatype == 1 ) then
     allocate(arr_r(numel)) 
     call dscal(numel, init_val_%re, arr_r, 1)
     tensor%data_real = arr_r
    elseif (tensor%datatype == 3 ) then
     allocate(arr_c(numel)) 
     call zscal(numel, init_val_, arr_c, 1)
     tensor%data_complex = arr_c
    endif
   end if

   construct_tapp_tensor_t = tapp_err
 end function construct_tapp_tensor_t

 ! initialize to a given value or zero by default. Data are already allocated previously in any case.
 function initialize_tapp_tensor_t(tensor, init_val)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor
   complex(c_double_complex), intent(in), optional    :: init_val
   complex(c_double_complex), target     :: init_val_
   real(c_double), target     :: init_val_r
   type(c_ptr)     :: tmp

   integer(c_int)       :: tapp_err, initialize_tapp_tensor_t

   init_val_ = 0.0
   if (present(init_val)) init_val_ = init_val

   tapp_err = 1
   if (.not. is_empty_tapp_tensor_t(tensor)) then
      if (tensor%datatype == 1 ) then
        init_val_r = init_val_%re
        tmp = c_loc(init_val_r)
      elseif (tensor%datatype == 3 ) then
        tmp = c_loc(init_val_)
      endif
     tapp_err = TAPP_initialize_tensor(tensor%tensor_info, tmp)
  else
     write(*,'("  Error. non-construced tensor in initialize_tapp_tensor_t.")')
     call quit(" Error. non-construced tensor in initialize_tapp_tensor_t.")
     stop 31
   endif

   initialize_tapp_tensor_t = tapp_err
 end function initialize_tapp_tensor_t

 ! copy_tapp_tensor_t
 function copy_tapp_tensor_t(dest, src)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: dest, src
   integer(c_int)       :: tapp_err, copy_tapp_tensor_t

   tapp_err = 1
   if (.not. is_empty_tapp_tensor_t(dest) .and. .not. is_empty_tapp_tensor_t(src)) then
     tapp_err = TAPP_copy_tensor(dest%tensor_info, src%tensor_info)
   else
     write(*,'("  Error. copying non-construced tensor in copy_tapp_tensor_t.")')
     call quit(" Error. copying non-construced tensor in copy_tapp_tensor_t.")
     stop 31
   endif

   copy_tapp_tensor_t = tapp_err
 end function copy_tapp_tensor_t

 ! destructor for tapp_tensor_t. Safe to call for an empty tensor as well.
 function destroy_tapp_tensor_t(tensor)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor

   integer(c_int)       :: tapp_err, destroy_tapp_tensor_t
   
   tapp_err = 0
   if (c_associated(tensor%tensor_info)) then
     tapp_err = TAPP_destroy_tensor_info(tensor%tensor_info)
   endif

   if (allocated(tensor%data_real)) deallocate(tensor%data_real)
   if (allocated(tensor%data_complex)) deallocate(tensor%data_complex)
    
   tensor%constructed = .FALSE.
   destroy_tapp_tensor_t = tapp_err
 end function destroy_tapp_tensor_t

 ! deallocates data stored in tensor_t temporarily after direct_data_pointer is called in case full tensor direct data access is needed
 function deallocate_temporary_full_data_tapp_tensor_t(tensor)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor
   integer(c_int)       :: tapp_err, deallocate_temporary_full_data_tapp_tensor_t

   tapp_err = 0

   if (allocated(tensor%data_real)) deallocate(tensor%data_real)
   if (allocated(tensor%data_complex)) deallocate(tensor%data_complex)

   deallocate_temporary_full_data_tapp_tensor_t = tapp_err
 end function deallocate_temporary_full_data_tapp_tensor_t



 ! get nmodes (rank) from tapp_tensor_t
 function get_nmodes_tapp_tensor_t(tensor, nmodes)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor
   integer(c_int), intent(out) :: nmodes

   integer(c_int)       :: tapp_err, get_nmodes_tapp_tensor_t
   
   tapp_err = 0
   nmodes = TAPP_get_nmodes(tensor%tensor_info)
   get_nmodes_tapp_tensor_t = tapp_err
 end function get_nmodes_tapp_tensor_t

 ! get extents (dimensions) of a tensor. Returns an empty array in case of a scalar
 function get_extents_tapp_tensor_t(tensor, extents)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor
   integer(c_int64_t), allocatable, target, intent(inout) :: extents(:)
   integer(c_int)          :: nmodes

   integer(c_int)       :: tapp_err, get_extents_tapp_tensor_t 

   tapp_err =  get_nmodes_tapp_tensor_t(tensor, nmodes)
   if (.not. allocated(extents)) then
     allocate(extents(nmodes))
   endif
   if (.not. size(extents) == nmodes) then 
     write(*,'("  Error. An array of wrong size passed to get_extents_tapp_tensor_t.")') 
     stop 58
   endif

   if (size(extents) .ne. 0) then 
     call TAPP_get_extents(tensor%tensor_info, c_loc(extents(1)))
   endif

   get_extents_tapp_tensor_t = tapp_err
 end function get_extents_tapp_tensor_t

 ! check if a tensor is empty (constructor not used)
 function is_empty_tapp_tensor_t(tensor)
   implicit none
   type(tapp_tensor_t), intent(inout)    :: tensor
   logical             :: ret, is_empty_tapp_tensor_t   

   is_empty_tapp_tensor_t = .not. tensor%constructed
 end function is_empty_tapp_tensor_t
   
 function scale_with_denominators_tapp_tensor_t(tensor, eps_occ, eps_vir, eps_ijk)
    implicit none
    type(tapp_tensor_t), target, intent(inout)    :: tensor
    real(8), intent(in), target                   :: eps_occ(:), eps_vir(:)
    real(8), intent(in), optional                 :: eps_ijk
    integer(c_int)       :: tapp_err
    complex(c_double_complex), allocatable, target               :: eps_occ_(:), eps_vir_(:)
    real(c_double), target                               :: eps_ijk_r
    complex(c_double_complex), target                            :: eps_ijk_
    integer(c_int) :: i, n_occ, n_vir
    integer(c_int) :: scale_with_denominators_tapp_tensor_t

    eps_ijk_r = 0.0
    if(present(eps_ijk)) eps_ijk_r = eps_ijk
    eps_ijk_ = eps_ijk_r

   if (is_empty_tapp_tensor_t(tensor)) then
     ! write(*,'("  Error. non-constructed tensor in scale_with_denominators_tapp_tensor_t.")')
     call quit(" Error. non-constructed tensor in scale_with_denominators_tapp_tensor_t.")
     stop 31
   endif

    n_occ = size(eps_occ)
    n_vir = size(eps_vir)
    tapp_err =  -1
    if (tensor%datatype == 1) then
      tapp_err =TAPP_scale_with_denominators(tensor%tensor_info, n_occ, n_vir, &
                  c_loc(eps_occ(1)), c_loc(eps_vir(1)),  c_loc(eps_ijk_r))
    elseif (tensor%datatype == 3) then
      allocate(eps_occ_(n_occ))
      allocate(eps_vir_(n_vir))
      do i = 1, n_occ
        eps_occ_(i) = eps_occ(i)
      end do
      do i = 1, n_vir
        eps_vir_(i) = eps_vir(i)
      end do
      tapp_err =TAPP_scale_with_denominators(tensor%tensor_info, n_occ, n_vir, &
                  c_loc(eps_occ_(1)), c_loc(eps_vir_(1)),  c_loc(eps_ijk_))
    endif

    if (allocated(eps_occ_)) deallocate(eps_occ_)
    if (allocated(eps_vir_)) deallocate(eps_vir_)
    scale_with_denominators_tapp_tensor_t = tapp_err
 end function scale_with_denominators_tapp_tensor_t

! stop workers and call mpi finalize
 function finalize_tapp_tensor_t()
   implicit none
   integer(c_int)       :: finalize_tapp_tensor_t

   finalize_tapp_tensor_t = TAPP_finalize()
 end function finalize_tapp_tensor_t


end module tapp_top
