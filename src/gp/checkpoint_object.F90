! Object-oriented helpers for working with composite checkpoint data.
!
! This module deliberately keeps all checkpoint I/O in the checkpoint module.
! The abstract base type only stores the root label and per-field state, and
! forwards typed reads/writes to checkpoint_read/checkpoint_write.

module checkpoint_object

  use checkpoint, only: checkpoint_query, checkpoint_read, checkpoint_write

  implicit none

  private

  public :: checkpoint_field_t
  public :: checkpoint_object_t

  type :: checkpoint_field_t
     character(len=:), allocatable :: name
     character(len=:), allocatable :: data_type
     character(len=:), allocatable :: rank
     logical                      :: required = .false.
     logical                      :: exists = .false.
     logical                      :: dirty = .false.
     integer                      :: data_size = -1
  end type checkpoint_field_t

  type, abstract :: checkpoint_object_t
     character(len=:), allocatable :: root
     type(checkpoint_field_t), allocatable :: fields(:)
  contains
     procedure :: set_root
     procedure :: register_field
     procedure :: label
     procedure :: query_field
     procedure :: mark_dirty
     procedure :: mark_clean
     procedure :: field_is_dirty
     procedure :: field_exists
     procedure :: field_size
     procedure :: read_integer
     procedure :: read_integer_array
     procedure :: read_real
     procedure :: read_real_array
     procedure :: read_string
     procedure :: read_string_array
     procedure :: write_integer
     procedure :: write_integer_array
     procedure :: write_real
     procedure :: write_real_array
     procedure :: write_string
     procedure :: write_string_array
  end type checkpoint_object_t

contains

  subroutine set_root(self, root)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: root

    self%root = trim(root)

  end subroutine set_root

  subroutine register_field(self, name, data_type, rank, required)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: name
    character(len=*), intent(in)              :: data_type
    character(len=*), intent(in)              :: rank
    logical, intent(in), optional             :: required

    type(checkpoint_field_t), allocatable     :: new_fields(:)
    integer                                   :: n

    n = field_index(self, name)
    if (n == 0) then
       if (allocated(self%fields)) then
          n = size(self%fields) + 1
          allocate(new_fields(n))
          new_fields(1:n-1) = self%fields
          call move_alloc(new_fields, self%fields)
       else
          n = 1
          allocate(self%fields(1))
       end if
    end if

    self%fields(n)%name = trim(name)
    self%fields(n)%data_type = trim(data_type)
    self%fields(n)%rank = trim(rank)
    self%fields(n)%required = .false.
    if (present(required)) self%fields(n)%required = required

  end subroutine register_field

  function label(self, field) result(full_label)

    class(checkpoint_object_t), intent(in) :: self
    character(len=*), intent(in)           :: field
    character(len=:), allocatable          :: full_label

    if (.not.allocated(self%root)) then
       stop 'checkpoint_object: root label is not set'
    end if
    full_label = trim(self%root)//'/'//trim(field)

  end function label

  subroutine query_field(self, field, exists, data_size)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    logical, intent(out)                      :: exists
    integer, intent(out), optional            :: data_size

    integer                                   :: i, local_size

    call checkpoint_query(self%label(field), exists, local_size)

    i = field_index(self, field)
    if (i > 0) then
       self%fields(i)%exists = exists
       self%fields(i)%data_size = local_size
    end if

    if (present(data_size)) data_size = local_size

  end subroutine query_field

  subroutine mark_dirty(self, field)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field

    call set_dirty_state(self, field, .true.)

  end subroutine mark_dirty

  subroutine mark_clean(self, field)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field

    call set_dirty_state(self, field, .false.)

  end subroutine mark_clean

  logical function field_is_dirty(self, field)

    class(checkpoint_object_t), intent(in) :: self
    character(len=*), intent(in)           :: field
    integer                                :: i

    field_is_dirty = .false.
    i = field_index(self, field)
    if (i > 0) field_is_dirty = self%fields(i)%dirty

  end function field_is_dirty

  logical function field_exists(self, field)

    class(checkpoint_object_t), intent(in) :: self
    character(len=*), intent(in)           :: field
    integer                                :: i

    field_exists = .false.
    i = field_index(self, field)
    if (i > 0) field_exists = self%fields(i)%exists

  end function field_exists

  integer function field_size(self, field)

    class(checkpoint_object_t), intent(in) :: self
    character(len=*), intent(in)           :: field
    integer                                :: i

    field_size = -1
    i = field_index(self, field)
    if (i > 0) field_size = self%fields(i)%data_size

  end function field_size

  subroutine read_integer(self, field, value)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(out)                      :: value

    call checkpoint_read(self%label(field), idata=value)
    call set_read_state(self, field, 1)

  end subroutine read_integer

  subroutine read_integer_array(self, field, values)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(out)                      :: values(:)

    call checkpoint_read(self%label(field), idata=values)
    call set_read_state(self, field, size(values))

  end subroutine read_integer_array

  subroutine read_real(self, field, value)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    real(8), intent(out)                      :: value

    call checkpoint_read(self%label(field), rdata=value)
    call set_read_state(self, field, 1)

  end subroutine read_real

  subroutine read_real_array(self, field, values)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    real(8), intent(out)                      :: values(:)

    call checkpoint_read(self%label(field), rdata=values)
    call set_read_state(self, field, size(values))

  end subroutine read_real_array

  subroutine read_string(self, field, value)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    character(len=*), intent(out)             :: value

    call checkpoint_read(self%label(field), sdata=value)
    call set_read_state(self, field, 1)

  end subroutine read_string

  subroutine read_string_array(self, field, values, string_length)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(in)                       :: string_length
    character(len=string_length), intent(out) :: values(:)

    call checkpoint_read(self%label(field), sarray=values, slen=string_length)
    call set_read_state(self, field, size(values))

  end subroutine read_string_array

  subroutine write_integer(self, field, value)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(in)                       :: value

    call checkpoint_write(self%label(field), idata=value)
    call set_written_state(self, field, 1)

  end subroutine write_integer

  subroutine write_integer_array(self, field, values)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(in)                       :: values(:)

    call checkpoint_write(self%label(field), idata=values)
    call set_written_state(self, field, size(values))

  end subroutine write_integer_array

  subroutine write_real(self, field, value)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    real(8), intent(in)                       :: value

    call checkpoint_write(self%label(field), rdata=value)
    call set_written_state(self, field, 1)

  end subroutine write_real

  subroutine write_real_array(self, field, values)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    real(8), intent(in)                       :: values(:)

    call checkpoint_write(self%label(field), rdata=values)
    call set_written_state(self, field, size(values))

  end subroutine write_real_array

  subroutine write_string(self, field, value)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    character(len=*), intent(in)              :: value

    call checkpoint_write(self%label(field), sdata=value)
    call set_written_state(self, field, 1)

  end subroutine write_string

  subroutine write_string_array(self, field, values, string_length)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(in)                       :: string_length
    character(len=string_length), intent(in)  :: values(:)

    call checkpoint_write(self%label(field), sarray=values, slen=string_length)
    call set_written_state(self, field, size(values))

  end subroutine write_string_array

  integer function field_index(self, field)

    class(checkpoint_object_t), intent(in) :: self
    character(len=*), intent(in)           :: field
    integer                                :: i

    field_index = 0
    if (.not.allocated(self%fields)) return

    do i = 1, size(self%fields)
       if (allocated(self%fields(i)%name)) then
          if (trim(self%fields(i)%name) == trim(field)) then
             field_index = i
             return
          end if
       end if
    end do

  end function field_index

  subroutine set_dirty_state(self, field, dirty)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    logical, intent(in)                       :: dirty
    integer                                   :: i

    i = field_index(self, field)
    if (i == 0) return
    self%fields(i)%dirty = dirty

  end subroutine set_dirty_state

  subroutine set_read_state(self, field, data_size)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(in)                       :: data_size
    integer                                   :: i

    i = field_index(self, field)
    if (i == 0) return
    self%fields(i)%exists = .true.
    self%fields(i)%dirty = .false.
    self%fields(i)%data_size = data_size

  end subroutine set_read_state

  subroutine set_written_state(self, field, data_size)

    class(checkpoint_object_t), intent(inout) :: self
    character(len=*), intent(in)              :: field
    integer, intent(in)                       :: data_size
    integer                                   :: i

    i = field_index(self, field)
    if (i == 0) return
    self%fields(i)%exists = .true.
    self%fields(i)%dirty = .false.
    self%fields(i)%data_size = data_size

  end subroutine set_written_state

end module checkpoint_object
