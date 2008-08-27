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

  type cube_grid_t
    integer                    :: nr_atoms
    real(kind_dp), allocatable :: atom_xyz(:)
    real(kind_dp), allocatable :: atom_charge(:)
    integer                    :: nr_points_along_side
    type(grid_t)               :: grid
    real(kind_dp)              :: step_vector(3, 3)
    real(kind_dp)              :: initial_xyz(3)
  end type

contains

  subroutine init_even_cube_grid(self, nr_points_along_side)

    type(cube_grid_t), intent(inout) :: self
    integer,           intent(in)    :: nr_points_along_side
!   ----------------------------------------------------------------------------
!   real(kind_dp), allocatable   :: work(:)
!   integer                      :: kfree, lfree
!   integer                      :: file_unit
!   integer                      :: nr_points
!   integer                      :: counter
    integer                          :: i, j, k, l, m

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"

!   integer                    :: nr_atoms
!   real(kind_dp), allocatable :: atom_xyz(:, :)
!   integer,       allocatable :: atom_charge(:)
!   integer                    :: nr_points_along_side
!   type(grid_t)               :: grid
!   real(kind_dp)              :: step_vector(3, 3)
!   real(kind_dp)              :: initial_xyz(3)

    self%nr_atoms = nucdep

    call alloc(self%atom_xyz, self%nr_atoms*3)
    call alloc(self%atom_charge, self%nr_atoms)

    self%nr_points_along_side = nr_points_along_side

    call init_grid(self%grid, nr_points_along_side**3)






!       NCUBE(0)      number of atoms
!       NCUBE(1)      number of points in the "x"-direction
!       NCUBE(2)      number of points in the "y"-direction
!       NCUBE(3)      number of points in the "z"-direction

!       STEP_3D(*,0)  coordinates of initial point
!       STEP_3D(*,1)  step vector for the "x"-direction
!       STEP_3D(*,2)  step vector for the "y"-direction
!       STEP_3D(*,3)  step vector for the "z"-direction

!       ATOM_XYZ(1,I) x-coordinate of atom I
!       ATOM_XYZ(2,I) y-coordinate of atom I
!       ATOM_XYZ(3,I) z-coordinate of atom I
!       ATOM_XYZ(4,I) charge of atom I

!       CALL DZERO(STEP_3D,12)
!       DO I = 1,3
!         STEP_3D(I,0) = CORD(I,1)
!         STEP_3D(I,I) = CORD(I,1)
!       ENDDO




        i = 0
        do j = 1, nucind
          k = istbnu(j)
          do l = 0, maxopr
            if (iand(l, k) == 0) then
              i = i + 1
              do m = 1, 3
!               atom_xyz(m,i) = pt(ibtand(isymax(m,1),l))*cord(m,j)

!               step_3d(m,0) = min(step_3d(m,0),atom_xyz(m,i))
!               step_3d(m,m) = max(step_3d(m,m),atom_xyz(m,i))
              end do
!             atom_xyz(4,i) = charge(j)
            end if
          end do
        end do






!       add some space, calculate step sizes and convert step vectors to angstrom

!       DO I = 1,3
!         STEP_3D(I,0) = STEP_3D(I,0) - ADD_3D
!         STEP_3D(I,I) = (STEP_3D(I,I) + (ADD_3D*D2))*XTANG
!         STEP_3D(I,I) = (STEP_3D(I,I) - STEP_3D(I,0))/(NCUBE(I) - 1)
!       ENDDO

!       write header

!       WRITE(LUNIT_3D,'(A)') 'some'
!       WRITE(LUNIT_3D,'(A)') 'text'
!       DO I = 0,3
!         WRITE(LUNIT_3D,'(I5,4F12.6)') NCUBE(I), (STEP_3D(J,I),J=1,3)
!       ENDDO
!       DO I = 1,NUCDEP
!         WRITE(LUNIT_3D,'(I5,4F12.6)')
!    &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
!       ENDDO

!       grid walk

!       DO ICUBE = 1,NCUBE(1)
!         DO JCUBE = 1,NCUBE(2)
!           DO KCUBE = 1,NCUBE(3)

!             PX = STEP_3D(1,0) + (ICUBE-1)*STEP_3D(1,1)
!    &                          + (JCUBE-1)*STEP_3D(1,2)
!    &                          + (KCUBE-1)*STEP_3D(1,3)
!             PY = STEP_3D(2,0) + (ICUBE-1)*STEP_3D(2,1)
!    &                          + (JCUBE-1)*STEP_3D(2,2)
!    &                          + (KCUBE-1)*STEP_3D(2,3)
!             PZ = STEP_3D(3,0) + (ICUBE-1)*STEP_3D(3,1)
!    &                          + (JCUBE-1)*STEP_3D(3,2)
!    &                          + (KCUBE-1)*STEP_3D(3,3)





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
