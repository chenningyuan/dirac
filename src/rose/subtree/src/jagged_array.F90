module jagged_array
  implicit none
  save
  public NewJA, DeleteJA, JASetColumn, JAGetColumn, JAGetElement

  type IntegerArrayType
       integer, allocatable :: column(:)
  end type

  type RealArrayType
       real(8), allocatable :: column(:)
  end type

  type JaggedArray
          private        
          character, public  :: array_type=' '
          integer, allocatable, public :: column_length(:)
          type(IntegerArrayType), pointer :: integer_array(:)
          type(RealArrayType), pointer :: real_array(:)
  end type
 
  ! Overloaded procedure interfaces
  interface NewJA
          module procedure New
  end interface

  interface DeleteJA
          module procedure Delete
  end interface

  interface JASetColumn
          module procedure SetColumn
  end interface

  interface JAGetColumn
          module procedure GetColumn
  end interface

  interface JAGetElement
          module procedure GetElement
  end interface

contains

  ! Constructor, allocates integer or real pointer array
  subroutine New(self, array_length, array_type)
    type(JaggedArray), intent(out) :: self
    character, intent(in) :: array_type
    integer, intent(in) :: array_length
    allocate(self%column_length(array_length))
    self%column_length = 0
    if (array_type == 'R' .or. array_type == 'r') then
       self%array_type = 'r'
       allocate (self%real_array(array_length))
    else if (array_type == 'I' .or. array_type == 'i') then
       self%array_type = 'i'
       allocate (self%integer_array(array_length))
    else
       error stop 'JaggedArray called with unknown type'
    end if
  end subroutine

  ! Destructor, deallocates array
  subroutine Delete(self)
    type (JaggedArray), intent(inout) :: self
    integer :: i
    deallocate(self%column_length)
    if (self%array_type == 'r') then
       do i = 1, size(self%real_array)
          deallocate (self%real_array(i)%column)
       end do
       deallocate (self%real_array)
    else if (self%array_type == 'i') then
       do i = 1, size(self%integer_array)
          deallocate (self%integer_array(i)%column)
       end do
       deallocate (self%integer_array)
    end if
  end subroutine

  ! Setter, initializes column
  subroutine SetColumn(self, i, integer_column, real_column)
    type (JaggedArray), intent(inout) :: self
    integer, intent(in) :: i
    integer, intent(in), optional :: integer_column(:)
    real(8), intent(in), optional :: real_column(:)
    integer :: length
    if (present(integer_column)) then
       length = size(integer_column)
       if (allocated(self%integer_array(i)%column)) deallocate(self%integer_array(i)%column)
       allocate(self%integer_array(i)%column(length))
       self%integer_array(i)%column = integer_column
    else if (present(real_column)) then
       length = size(real_column)
       if (allocated(self%real_array(i)%column)) deallocate(self%real_array(i)%column)
       allocate(self%real_array(i)%column(length))
       self%real_array(i)%column = real_column
    end if 
    self%column_length(i) = length
  end subroutine

  ! Getter, gets entire column
  subroutine GetColumn(self,i,integer_column,real_column)
    type (JaggedArray), intent(in) :: self
    integer, intent(in) :: i
    integer, intent(out), optional :: integer_column(:)
    real(8), intent(out), optional :: real_column(:)
    if (present(integer_column)) then
       integer_column = self%integer_array(i)%column
    else if (present(real_column)) then
       real_column = self%real_array(i)%column
    end if
  end subroutine

  ! Getter, gets element value
  subroutine GetElement(self,i,j,integer_value,real_value)
    type (JaggedArray), intent(in) :: self
    integer, intent(in) :: i,j
    integer, intent(out), optional :: integer_value
    real(8), intent(out), optional :: real_value
    if (present(integer_value)) then
       integer_value = self%integer_array(j)%column(i)
    else if (present(real_value)) then
       real_value = self%real_array(j)%column(i)
    end if
  end subroutine

end module
