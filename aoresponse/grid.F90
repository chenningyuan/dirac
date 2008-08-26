module grid_m

  use kinds_m
  use memory_allocator
  use matrix_oper
  use matrix_mod

  implicit none

  type grid_t
    integer                    :: nr_points
    real(kind_dp), allocatable :: x(:)
    real(kind_dp), allocatable :: y(:)
    real(kind_dp), allocatable :: z(:)
    real(kind_dp), allocatable :: w(:)
  end type

contains

  subroutine init_dft_grid(D, self)

    type(matrix),  intent(inout) :: D
    type(grid_t),  intent(inout) :: self
!   ----------------------------------------------------------------------------
    real(kind_dp), allocatable   :: work(:)
    integer                      :: kfree, lfree
    integer                      :: file_unit
    integer                      :: nr_points
    integer                      :: counter
    integer                      :: i

    file_unit = 44

    kfree = 1
    lfree = 10000000
    call alloc(work, lfree)

    call dftgrd(D%elms, work(kfree), lfree, 1)

    call dealloc(work)


!   first round: get the total nr_points
!   ====================================

    rewind file_unit
    nr_points = 137 ! just set it to something > 0

    counter = 0
    do while (nr_points > 0)

      read(file_unit) nr_points

      if (nr_points > 0) then

        read(file_unit)
        counter = counter + nr_points

      end if
    end do


!   second round: allocate and read
!   ===============================

    self%nr_points = counter

    call alloc(self%x, self%nr_points)
    call alloc(self%y, self%nr_points)
    call alloc(self%z, self%nr_points)
    call alloc(self%w, self%nr_points)

    rewind file_unit
    nr_points = 137 ! just set it to something > 0

    counter = 0
    do while (nr_points > 0)

      read(file_unit) nr_points

      if (nr_points > 0) then

        call reaqua(self%x(1+counter:1+counter+nr_points), &
                    self%y(1+counter:1+counter+nr_points), &
                    self%z(1+counter:1+counter+nr_points), &
                    self%w(1+counter:1+counter+nr_points), &
                    file_unit,           &
                    nr_points)

        counter = counter + nr_points

      end if
    end do

    rewind file_unit
!   close(file_unit, status = 'keep')

  end subroutine

  subroutine del_grid(self)

    type(grid_t), intent(inout) :: self

    call dealloc(self%x)
    call dealloc(self%y)
    call dealloc(self%z)
    call dealloc(self%w)

  end subroutine

end module
