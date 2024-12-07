module tapp_utils
! This module contains utilities for tha usage of TAPP standard tensor interface
! Written by Jan Brandejs, Dec 2024

 use talsh
 use tensor_algebra
 use, intrinsic:: ISO_C_BINDING
 use tapp
    
 implicit none
 
 contains
 function tapp_einsum(expression, D_talsh, A_talsh, B_talsh, scale)
    implicit none
    character(len=*), intent(in) :: expression
    type(talsh_tens_t), intent(inout) :: D_talsh, A_talsh, B_talsh
    COMPLEX(c_double_complex), optional, intent(in) :: scale

    ! TAPP_create_tensor_info
    type(C_PTR) :: info_A,info_B,info_C,info_D ! integer(c_intptr_t)
    integer(c_int) :: dtype ! data  type 0.. tapp f32
    integer(c_int) :: nmode_A, nmode_B, nmode_C, nmode_D
    integer(c_int) :: tapp_err
    integer(c_int64_t), allocatable, target :: arr_extents_A(:), arr_extents_B(:), arr_extents_D(:)
    integer(c_int64_t), allocatable, target :: arr_strides_A(:), arr_strides_B(:), arr_strides_D(:)
    type(C_PTR):: extents_A, extents_B, extents_C, extents_D
    type(C_PTR):: strides_A, strides_B, strides_C, strides_D


    ! TAPP_create_tensor_product 
    type(C_PTR) :: handle ! TAPP_handle
    type(C_PTR) :: plan ! TAPP_tensor_product
    integer(c_int) :: op_A, op_B, op_C, op_D ! 0 .. identitity
    integer(c_int64_t), allocatable, target :: arr_idx_A(:), arr_idx_B(:), arr_idx_C(:), arr_idx_D(:)
    type(C_PTR):: idx_A, idx_B, idx_C, idx_D
    integer(c_int) :: prec !  TAPP_prectype -1 .. default


    ! TAPP_execute_product
    type(C_PTR) :: exec ! TAPP_executor
    type(C_PTR) :: status ! TAPP_status
    COMPLEX(c_double_complex), target :: alpha, beta
    type(C_PTR) :: alpha_ptr, beta_ptr 
    ! COMPLEX(c_double_complex), allocatable, target :: arr_A(:), arr_B(:), arr_C(:), arr_D(:)
    type(C_PTR):: A, B, C, D

    ! TAPP_get_nmodes
    integer(c_int) :: nmodes

    ! TAPP_get_extents
    integer(c_int64_t), allocatable, target :: arr_extents(:)
    type(C_PTR):: extents

    ! internal + talsh
    integer :: i, k, stride, numel
    integer(INTD), allocatable, target :: out_dimsA(:), out_dimsB(:), out_dimsD(:)
    integer(INTD), allocatable, target :: lin_dim(:)
    complex(8), pointer  :: lin_tens(:)

    ! err
    integer :: ierr
    integer :: tapp_einsum


    ierr = parse_einsum(expression, arr_idx_A, arr_idx_B, arr_idx_C, arr_idx_D)

    nmode_A = talsh_tensor_rank(A_talsh)
    allocate(out_dimsA(nmode_A))
    ierr = talsh_tensor_dimensions(A_talsh,nmode_A,out_dimsA)

    allocate(arr_extents_A(nmode_A), arr_strides_A(nmode_A))
    arr_extents_A = int(out_dimsA, kind=8)

    do i=1,nmode_A
      ! fortran memory order, so e.g. for 
      ! arr_extents_A   = (/4, 3, 5/)
      ! we want
      ! arr_strides_A   = (/12, 4, 1/)
      stride = 1
      !do k=1,nmode_A-i
      do k=1,i-1
        stride = stride * arr_extents_A(k)
      enddo
      arr_strides_A(i) = stride
    enddo
    
    dtype = 3 ! 0 TAPP_F32, 1 TAPP_F64, 3 TAPP_C64

    extents_A = c_loc(arr_extents_A(1)) 
    strides_A = c_loc(arr_strides_A(1)) 
    tapp_err = TAPP_create_tensor_info(info_A, dtype, nmode_A, extents_A, strides_A)
    if (tapp_err.ne.0) then
      write(*,'("  Error. Creating TAPP tensor info failed with error code: ",i11,".")') tapp_err
      stop 40
    end if



    nmode_B = talsh_tensor_rank(B_talsh)
    allocate(out_dimsB(nmode_B))
    ierr = talsh_tensor_dimensions(B_talsh,nmode_B,out_dimsB)

    allocate(arr_extents_B(nmode_B), arr_strides_B(nmode_B))
    arr_extents_B = int(out_dimsB, kind=8)
    do i=1,nmode_B
      stride = 1
      !do k=1,nmode_B-i
      do k=1,i-1
        stride = stride * arr_extents_B(k)
      enddo
      arr_strides_B(i) = stride
    enddo

    extents_B = c_loc(arr_extents_B(1)) 
    strides_B = c_loc(arr_strides_B(1)) 
    tapp_err = TAPP_create_tensor_info(info_B, dtype, nmode_B, extents_B, strides_B)



    nmode_D = talsh_tensor_rank(D_talsh)
    allocate(out_dimsD(nmode_D))
    ierr = talsh_tensor_dimensions(D_talsh,nmode_D,out_dimsD)

    allocate(arr_extents_D(nmode_D), arr_strides_D(nmode_D))
    arr_extents_D = int(out_dimsD, kind=8)
    do i=1,nmode_D
      stride = 1
      !do k=1,nmode_D-i
      do k=1,i-1
        stride = stride * arr_extents_D(k)
      enddo
      arr_strides_D(i) = stride
    enddo

    extents_D = c_loc(arr_extents_D(1)) 
    strides_D = c_loc(arr_strides_D(1)) 
    tapp_err = TAPP_create_tensor_info(info_D, dtype, nmode_D, extents_D, strides_D)


    nmodes = TAPP_get_nmodes(info_A)
    allocate(arr_extents(nmodes))
    arr_extents = -1
    extents = c_loc(arr_extents(1)) 

    call TAPP_get_extents(info_A, extents)

    nmode_C = nmode_D

    handle = c_null_ptr
    op_A = 0 ! identity
    op_B = 0
    op_C = 0
    op_D = 0

    idx_A = c_loc(arr_idx_A(1)) 
    idx_B = c_loc(arr_idx_B(1)) 
    idx_D = c_loc(arr_idx_D(1)) 
    prec = -1 ! default prec
    tapp_err = TAPP_create_tensor_product(plan, handle, op_A, info_A, idx_A, op_B, info_B, idx_B, &
      op_C, info_D, idx_D, op_D, info_D, idx_D, prec)

    if (tapp_err.ne.0) then
      write(*,'("  Error. Creating TAPP contraction plan failed with error code: ",i11,".")') tapp_err
      stop 41
    end if
    
    exec = c_null_ptr
    alpha = 1.0
    if (present(scale)) then
      alpha = scale
    end if

    beta = 1.0
    
    alpha_ptr = c_loc(alpha)
    beta_ptr = c_loc(beta)

    ierr=talsh_tensor_get_body_access(A_talsh,A,C8,int(0,C_INT),DEV_HOST)
    ierr=talsh_tensor_get_body_access(B_talsh,B,C8,int(0,C_INT),DEV_HOST)
    ierr=talsh_tensor_get_body_access(D_talsh,D,C8,int(0,C_INT),DEV_HOST)

    tapp_err = TAPP_execute_product(plan, exec, status, alpha_ptr, A, B, beta_ptr, D, D) ! diff between plan par here and elsewhere?

    ierr = tapp_err
    if (tapp_err.ne.0) then
      write(*,'("  Error. TAPP tensor product failed with error code: ",i11,".")') tapp_err
      stop 42
    end if

    ! debug print
    ! numel = 1
    ! do i=1,nmode_D
    !   numel = numel * arr_extents_D(i)
    ! enddo
    ! allocate(lin_dim(1))
    ! lin_dim(1) = numel
    ! call c_f_pointer(D,lin_tens,lin_dim)
    ! write(*,*) " result array D:"
    ! do i=1,numel
    !   write (*,'(" (",F10.4,",",F10.4,")")') lin_tens(i)
    ! enddo

    tapp_err = TAPP_destory_tensor_product(plan)

    tapp_err = TAPP_destory_tensor_info(info_A)
    tapp_err = TAPP_destory_tensor_info(info_B)
    ! tapp_err = TAPP_destory_tensor_info(info_C)
    ! write(*,'("  destroy info_C ",i11,".")') tapp_err
    tapp_err = TAPP_destory_tensor_info(info_D)

    deallocate(arr_extents_A, arr_extents_B, arr_extents_D)
    deallocate(arr_strides_A, arr_strides_B, arr_strides_D)
    if (allocated(arr_idx_A)) deallocate(arr_idx_A)
    if (allocated(arr_idx_B)) deallocate(arr_idx_B)
    if (allocated(arr_idx_C)) deallocate(arr_idx_C)
    if (allocated(arr_idx_D)) deallocate(arr_idx_D)

    if (allocated(lin_dim)) deallocate(lin_dim)
    deallocate(arr_extents, out_dimsA, out_dimsB, out_dimsD)

  tapp_einsum = ierr
 end function tapp_einsum
 function parse_einsum(expression, arr_idx_A, arr_idx_B, arr_idx_C, arr_idx_D)
  implicit none
  character(len=*), intent(in) :: expression
  character(len=50), allocatable :: tensor_names(:)
  character(len=10), allocatable :: indices(:,:)
  integer, allocatable :: index_count(:)
  character(len=50), allocatable :: distinct_indices(:)
  integer(c_int64_t), allocatable, target, intent(inout) :: arr_idx_A(:), &
   arr_idx_B(:), arr_idx_C(:), arr_idx_D(:)
  integer(c_int64_t), allocatable :: map_to(:)
  
  integer(c_int64_t) :: i8
  integer :: ierr, i, j, k
  integer :: parse_einsum

  ierr = pre_parse_einsum(expression, tensor_names, indices, index_count)
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
  end if 
  if (size(indices, 1) > 1) then
    allocate(arr_idx_A(index_count(2)))
    do j = 1, index_count(2)
      arr_idx_A(j) = map_str_to_i8(indices(2,j), distinct_indices, map_to)
    end do
  end if 
  if (size(indices, 1) > 2) then
    allocate(arr_idx_B(index_count(3)))
    do j = 1, index_count(3)
      arr_idx_B(j) = map_str_to_i8(indices(3,j), distinct_indices, map_to)
    end do
  end if 
  if (size(indices, 1) > 3) then
    allocate(arr_idx_C(index_count(4)))
    do j = 1, index_count(4)
      arr_idx_C(j) = map_str_to_i8(indices(4,j), distinct_indices, map_to)
    end do
  end if 

  deallocate(tensor_names, indices, index_count, distinct_indices, map_to)
  parse_einsum = ierr
 end function parse_einsum 

 function pre_parse_einsum(expression, tensor_names, indices, index_count)
  implicit none

  character(len=*), intent(in) :: expression
  character(len=50), allocatable, intent(out) :: tensor_names(:)
  character(len=10), allocatable, intent(out) :: indices(:,:)
  integer, allocatable, intent(out) :: index_count(:)
  integer :: pre_parse_einsum

  ! Internal variables
  character(len=50) :: current_tensor
  character(len=10) :: current_index
  integer :: i, j, k, l, parenthesis_level, ierr, ntens

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
        current_tensor = ""
        k = 1
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
    elseif (expression(i:i) /= '+' .and. &
            expression(i:i) /= '=' .and. &
            expression(i:i) /= '*' .and. &
            expression(i:i) /= ' ') then
      if (parenthesis_level == 0) then
        current_tensor(k:k) = expression(i:i)
        k = k + 1
      else
        current_index(j:j) = expression(i:i)
        j = j + 1
      end if
    end if
  end do

  ierr = 0
  pre_parse_einsum = ierr
 end function pre_parse_einsum 

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

end module tapp_utils
