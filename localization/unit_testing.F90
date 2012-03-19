module unit_testing

   implicit none

   public compare_matrix
   public compare_files
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

   subroutine compare_files(iunit, file1, file2, thresh)
   implicit none

   ! External variables
   character(*), intent(in) :: file1, file2
   integer,      intent(in) :: iunit
   real*8,       intent(in) :: thresh

   ! Internal variables
   integer :: iu1, iu2, ier, ier1, ier2
   real*8  :: dnum1, dnum2, dmax

   iu1 = iunit
   iu2 = iunit + 1

   ! open file1
   open(iu1,file=file1,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_files: ERROR while opening file ',file1
     write(0,*) 'compare_files: ERROR while opening file ',file1
     stop 1
   endif

   ! open file2
   open(iu2,file=file2,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_files: ERROR while opening file ',file2
     write(0,*) 'compare_files: ERROR while opening file ',file2
     stop 1
   endif

   ! compare two files
   dmax = 0.0d0
   do
     read(iu1,iostat=ier1) dnum1
     read(iu2,iostat=ier2) dnum2
     if( (ier1 == -1) .or. (ier2 == -1) ) exit
     if( (ier1 /=  0) .or. (ier2 /=  0) )then
       write(6,*) 'compare_files: ERROR reading one of the files ',file1,file2
       write(0,*) 'compare_files: ERROR reading one of the files ',file1,file2
       stop 1
     endif

     if( (dabs(dnum1-dnum2)) > dmax ) dmax = dabs(dnum1-dnum2)

   enddo

   ! close file1
   close(iu1,status='keep')

   ! close file2
   close(iu2,status='keep')

   ! check result
   if(dmax > thresh) stop 1

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
