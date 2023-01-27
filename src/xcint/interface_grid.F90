!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module interface_grid

   implicit none

   public num_grid_write
   public num_grid_read

   private

contains

   subroutine num_grid_write(          &
                             x,        &
                             y,        &
                             z,        &
                             w,        &
                             lunit,    &
                             nr_points &
                            )

!     --------------------------------------------------------------------------
      real(8), intent(in) :: x(*)
      real(8), intent(in) :: y(*)
      real(8), intent(in) :: z(*)
      real(8), intent(in) :: w(*)
      integer, intent(in) :: lunit
      integer, intent(in) :: nr_points
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      write(lunit, *) nr_points
      if (nr_points > 0) then
         do i = 1, nr_points
            write(lunit, '(4d20.12)') x(i), y(i), z(i), w(i)
         enddo
      end if

   end subroutine

   subroutine num_grid_read(          &
                            x,        &
                            y,        &
                            z,        &
                            w,        &
                            lunit,    &
                            nr_points &
                           )

!     --------------------------------------------------------------------------
      real(8), intent(out) :: x(*)
      real(8), intent(out) :: y(*)
      real(8), intent(out) :: z(*)
      real(8), intent(out) :: w(*)
      integer, intent(in)  :: lunit
      integer, intent(in)  :: nr_points
!     --------------------------------------------------------------------------
      integer              :: i
!     --------------------------------------------------------------------------

      do i = 1, nr_points
         read(lunit, *) x(i), y(i), z(i), w(i)
      enddo

   end subroutine

end module
