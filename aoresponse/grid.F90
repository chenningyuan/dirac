module grid_m

  use kinds_m
  use memory_allocator
  use matrix_genop
  use matrix_defop
  use xc_driver

  implicit none

  type grid_t
    integer                    :: nr_points
    real(kind_dp), allocatable :: x(:)
    real(kind_dp), allocatable :: y(:)
    real(kind_dp), allocatable :: z(:)
    real(kind_dp), allocatable :: w(:)
  end type

  type cube_t
    integer                    :: nr_atoms
    real(kind_dp), allocatable :: atom_xyz(:, :)
    integer,       allocatable :: atom_charge(:)
    integer                    :: nr_points
    real(kind_dp)              :: initial_xyz(3)
    real(kind_dp)              :: step_vector(3, 3)
    real(kind_dp)              :: add_space = 6.0d0
  end type

contains

  subroutine init_cube(self, nr_points)

    type(cube_t), intent(inout) :: self
    integer,      intent(in)    :: nr_points
!   ----------------------------------------------------------------------------
    integer                     :: i, j, k, l, m

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"

    self%nr_atoms = nucdep

    call alloc(self%atom_xyz, self%nr_atoms, 3)
    call alloc(self%atom_charge, self%nr_atoms)

    self%nr_points = nr_points

    self%initial_xyz = 0.0d0
    self%step_vector = 0.0d0

    do i = 1, 3
      self%initial_xyz(i)    = cord(i, 1)
      self%step_vector(i, i) = cord(i, 1)
    end do

    i = 0
    do j = 1, nucind
      k = istbnu(j)
      do l = 0, maxopr
        if (iand(l, k) == 0) then
          i = i + 1
          do m = 1, 3
            self%atom_xyz(i, m)    = pt(iand(isymax(m, 1), l))*cord(m, j)
            self%initial_xyz(m)    = min(self%initial_xyz(m),    self%atom_xyz(i, m))
            self%step_vector(m, m) = max(self%step_vector(m, m), self%atom_xyz(i, m))
          end do
          self%atom_charge(i) = charge(j)
        end if
      end do
    end do

    do i = 1, 3
      self%initial_xyz(i)    = self%initial_xyz(i)    - self%add_space
      self%step_vector(i, i) = self%step_vector(i, i) + self%add_space*2.0d0

      self%step_vector(i, i) = self%step_vector(i, i)*0.5291772083d0
      self%step_vector(i, i) = (self%step_vector(i, i) - self%initial_xyz(i))/(self%nr_points - 1)
    end do

  end subroutine

  subroutine del_cube(self)

    type(cube_t), intent(inout) :: self

    call dealloc(self%atom_xyz)
    call dealloc(self%atom_charge)

  end subroutine

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

    call init_grid(self, counter)

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

  subroutine init_grid(self, nr_points)

    type(grid_t), intent(inout) :: self
    integer,      intent(in)    :: nr_points

    self%nr_points = nr_points

    call alloc(self%x, self%nr_points)
    call alloc(self%y, self%nr_points)
    call alloc(self%z, self%nr_points)
    call alloc(self%w, self%nr_points)

  end subroutine

  subroutine del_grid(self)

    type(grid_t), intent(inout) :: self

    call dealloc(self%x)
    call dealloc(self%y)
    call dealloc(self%z)
    call dealloc(self%w)

  end subroutine

end module
