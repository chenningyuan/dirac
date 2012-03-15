module unit_testing

   implicit none

   public compare_matrix
   public matrix_to_file
   public matrix_from_file

   private

   integer, parameter :: file_unit = 137
   logical            :: file_exists

contains

   subroutine compare_matrix(file_name, matrix, tolerance)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: file_name
      real(8)                  :: matrix(:, :)
      real(8),      intent(in) :: tolerance
!     --------------------------------------------------------------------------
      integer                  :: i, j, l, n
      real(8)                  :: f, diff, max_diff
!     --------------------------------------------------------------------------

      l = size(matrix)**0.5

      inquire(file = file_name, exist = file_exists)
      if (file_exists) then
         open(file_unit,             &
              file   = file_name,    &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
      else
         print *, 'error: file '//file_name//' not found'
         stop 1
      end if
      rewind(file_unit)

      read(file_unit, *) n
      if (n /= l*l) then
         print *, 'error: matrix dimensions do not match in file '//file_name
         stop 1
      end if

      max_diff = 0.0d0
      do i = 1, l
         do j = 1, l
            read(file_unit, *) f
            diff = abs(f - matrix(j, i))
            if (diff > max_diff) then
               max_diff = diff
            end if
         end do
      end do

      close(file_unit, status = 'keep')

      if (max_diff > tolerance) then
         stop 1
      end if

   end subroutine

   subroutine matrix_to_file(file_name, matrix)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: file_name
      real(8)                  :: matrix(:, :)
!     --------------------------------------------------------------------------
      integer                  :: i, j, l
!     --------------------------------------------------------------------------

      inquire(file = file_name, exist = file_exists)
      if (file_exists) then
         open(file_unit,             &
              file   = file_name,    &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
         close(file_unit, status = 'delete')
      end if
      open(file_unit,             &
           file   = file_name,    &
           status = 'new',        &
           form   = 'formatted',  &
           access = 'sequential')
      rewind(file_unit)

      l = size(matrix)**0.5

      write(file_unit, *) l*l
      do i = 1, l
         do j = 1, l
            write(file_unit, *) matrix(j, i)
         end do
      end do

      close(file_unit, status = 'keep')

   end subroutine

   subroutine matrix_from_file(file_name, matrix)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: file_name
      real(8)                  :: matrix(:, :)
!     --------------------------------------------------------------------------
      integer                  :: i, j, l, n
!     --------------------------------------------------------------------------

      l = size(matrix)**0.5

      inquire(file = file_name, exist = file_exists)
      if (file_exists) then
         open(file_unit,             &
              file   = file_name,    &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
      else
         print *, 'error: file '//file_name//' not found'
         stop 1
      end if
      rewind(file_unit)

      read(file_unit, *) n
      if (n /= l*l) then
         print *, 'error: matrix dimensions do not match in file '//file_name
         stop 1
      end if

      do i = 1, l
         do j = 1, l
            read(file_unit, *) matrix(j, i)
         end do
      end do

      close(file_unit, status = 'keep')

   end subroutine

end module
