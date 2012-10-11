module unit_testing

   implicit none

   public compare_matrix
   public compare_vector
   public compare_2vectors
   public compare_files
   public compare_text_files
   public matrix_to_file
   public matrix_from_file

interface compare_2vectors
   module procedure compare_1d_real_vectors
   module procedure compare_2d_real_vectors
   module procedure compare_3d_real_vectors
   module procedure compare_4d_real_vectors
   module procedure compare_5d_real_vectors
   module procedure compare_1d_complex_vectors
   module procedure compare_2d_complex_vectors
end interface compare_2vectors

   private

   integer, parameter :: file_unit = 137
   logical            :: file_exists

   private :: strip_blank, str_size, print_err_info
   private :: compare_real_vectors, compare_complex_vectors

interface print_err_info
   module procedure print_real_err_info
   module procedure print_complex_err_info
end interface print_err_info

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
         stop
      end if
      rewind(file_unit)

      read(file_unit, *) n
      if (n /= l*l) then
         print *, 'error: matrix dimensions do not match in file '//file_name
         stop
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
         stop
      end if

   end subroutine compare_matrix

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_vector(iunit, file_name, vector, thresh)
   implicit none

   ! External variables
   character(*), intent(in) :: file_name
   integer,      intent(in) :: iunit
   real*8,       intent(in) :: thresh, vector(*)

   ! Internal variables
   integer :: ier, ic, icmax
   real*8  :: dnum, errmax, sdev, dval1, dval2, diff

   ! open file_name
   open(iunit,file=file_name,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_vector: ERROR while opening file ',file_name
     write(0,*) 'compare_vector: ERROR while opening file ',file_name
     stop
   endif
   rewind iunit

   ! compare data from file and vector
   errmax = 0.0d0
   sdev   = 0.0d0
   dval1  = 0.0d0
   dval2  = 0.0d0
   icmax  = 0
   ic     = 0
   do

     ic = ic + 1
     read(iunit,*,iostat=ier) dnum
     if(ier == -1)exit
     if(ier /=  0)then
       write(6,*) 'compare_vector: ERROR reading file ',file_name
       write(0,*) 'compare_vector: ERROR reading file ',file_name
       write(6,*) 'ier: ',ier
       stop
     endif

     diff = dabs(dnum-vector(ic))
     sdev = sdev + diff*diff

     if( diff > errmax )then
       errmax = diff
       dval1  = dnum
       dval2  = vector(ic)
       icmax  = ic
     endif

   enddo

   ! close file_name
   close(iunit,status='keep')

   ! info
   ic = ic - 1
   sdev = dsqrt(sdev/dfloat(ic))
   call print_err_info(ic,thresh,errmax,dval1,dval2,icmax,sdev)

   ! check result
   if(errmax > thresh)then
     write(6,'(a)')
     write(6,'(a)') ' ERROR in compare_vector routine: Examined file and vector differ'
     write(0,'(a)') ' ERROR in compare_vector routine: Examined file and vector differ'
     stop
   endif

   end subroutine compare_vector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_1d_real_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer, intent(in) :: ndim
   real*8,  intent(in) :: thresh, vector1(:), vector2(:)
   call compare_real_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_1d_real_vectors

   subroutine compare_2d_real_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer, intent(in) :: ndim
   real*8,  intent(in) :: thresh, vector1(:,:), vector2(:,:)
   call compare_real_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_2d_real_vectors

   subroutine compare_3d_real_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer, intent(in) :: ndim
   real*8,  intent(in) :: thresh, vector1(:,:,:), vector2(:,:,:)
   call compare_real_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_3d_real_vectors

   subroutine compare_4d_real_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer, intent(in) :: ndim
   real*8,  intent(in) :: thresh, vector1(:,:,:,:), vector2(:,:,:,:)
   call compare_real_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_4d_real_vectors

   subroutine compare_5d_real_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer, intent(in) :: ndim
   real*8,  intent(in) :: thresh, vector1(:,:,:,:,:), vector2(:,:,:,:,:)
   call compare_real_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_5d_real_vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_real_vectors(ndim, vector1, vector2, thresh)
   implicit none

   ! External variables
   integer, intent(in) :: ndim
   real*8,  intent(in) :: thresh, vector1(*), vector2(*)

   ! Internal variables
   integer :: ic, icmax
   real*8  :: sdev, errmax, diff, dval1, dval2

   ! compare data from file and vector
   errmax = 0.0d0
   sdev   = 0.0d0
   dval1  = 0.0d0
   dval2  = 0.0d0
   icmax  = 0
   do ic=1, ndim

     diff = dabs(vector1(ic)-vector2(ic))
     sdev = sdev + diff*diff

     if( diff > errmax )then
       errmax = diff
       dval1  = vector1(ic)
       dval2  = vector2(ic)
       icmax  = ic
     endif

   enddo

   ! info
   sdev = dsqrt(sdev/dfloat(ndim))
   call print_err_info(ndim,thresh,errmax,dval1,dval2,icmax,sdev)

   ! check result
   if(errmax > thresh)then
     write(6,'(a)')
     write(6,'(a)') ' ERROR in compare_real_vectors routine: Examined arrays differ too much'
     write(0,'(a)') ' ERROR in compare_real_vectors routine: Examined arrays differ too much'
     stop
   endif

   end subroutine compare_real_vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_1d_complex_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer,    intent(in) :: ndim
   real*8,     intent(in) :: thresh
   complex*16, intent(in) :: vector1(:), vector2(:)
   call compare_complex_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_1d_complex_vectors

   subroutine compare_2d_complex_vectors(ndim, vector1, vector2, thresh)
   implicit none
   integer,    intent(in) :: ndim
   real*8,     intent(in) :: thresh
   complex*16, intent(in) :: vector1(:,:), vector2(:,:)
   call compare_complex_vectors(ndim, vector1, vector2, thresh)
   end subroutine compare_2d_complex_vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_complex_vectors(ndim, vector1, vector2, thresh)
   implicit none

   ! External variables
   integer,    intent(in) :: ndim
   real*8,     intent(in) :: thresh
   complex*16, intent(in) :: vector1(*), vector2(*)

   ! Internal variables
   integer :: ic, icmax
   real*8  :: sdev, errmax, diff
   complex*16 :: dval1, dval2

   ! compare data from file and vector
   errmax = 0.0d0
   sdev   = 0.0d0
   dval1  = 0.0d0
   dval2  = 0.0d0
   icmax  = 0
   do ic=1, ndim

     diff = (dreal(vector1(ic))-dreal(vector2(ic))) * (dreal(vector1(ic))-dreal(vector2(ic))) &
          + (dimag(vector1(ic))-dimag(vector2(ic))) * (dimag(vector1(ic))-dimag(vector2(ic)))
     diff = dsqrt(diff)
     sdev = sdev + diff*diff

     if( diff > errmax )then
       errmax = diff
       dval1  = vector1(ic)
       dval2  = vector2(ic)
       icmax  = ic
     endif

   enddo

   ! info
   sdev = dsqrt(sdev/dfloat(ndim))
   call print_err_info(ndim,thresh,errmax,dval1,dval2,icmax,sdev)

   ! check result
   if(errmax > thresh)then
     write(6,'(a)')
     write(6,'(a)')' ERROR in compare_complex_vectors routine: Examined arrays differ too much'
     write(0,'(a)')' ERROR in compare_complex_vectors routine: Examined arrays differ too much'
     stop
   endif

   end subroutine compare_complex_vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_files(iunit, file1, file2, thresh)
   implicit none

   ! External variables
   character(*), intent(in) :: file1, file2
   integer,      intent(in) :: iunit
   real*8,       intent(in) :: thresh

   ! Internal variables
   integer :: iu1, iu2, ier, ier1, ier2, ic, icmax
   real*8  :: dnum1, dnum2, errmax, sdev, dval1, dval2, diff

   iu1 = iunit
   iu2 = iunit + 1

   ! open file1
   open(iu1,file=file1,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_files: ERROR while opening file ',file1
     write(0,*) 'compare_files: ERROR while opening file ',file1
     stop
   endif
   rewind iu1

   ! open file2
   open(iu2,file=file2,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_files: ERROR while opening file ',file2
     write(0,*) 'compare_files: ERROR while opening file ',file2
     stop
   endif
   rewind iu2

   ! compare two files
   errmax = 0.0d0
   sdev   = 0.0d0
   dval1  = 0.0d0
   dval2  = 0.0d0
   icmax  = 0
   ic     = 0
   do

     ic = ic + 1
     read(iu1,*,iostat=ier1) dnum1
     read(iu2,*,iostat=ier2) dnum2
     if( (ier1 == -1) .and. (ier2 == -1) ) exit
     if( (ier1 == -1) .or.  (ier2 == -1) )then
       write(6,*) ' ERROR in compare_files: Examined files have different length '
       write(0,*) ' ERROR in compare_files: Examined files have different length '
       stop
     endif
     if( (ier1 /=  0) .or. (ier2 /=  0) )then
       write(6,*) ' ERROR in compare_files: reading one of the files ',file1,' ',file2
       write(0,*) ' ERROR in compare_files: reading one of the files ',file1,' ',file2
       write(6,*) ' ier1, ier2: ',ier1, ier2
       stop
     endif

     diff = dabs(dnum1-dnum2)
     sdev = sdev + diff*diff

     if( diff > errmax )then
       errmax = diff
       dval1  = dnum2  ! reference file
       dval2  = dnum1
       icmax  = ic
     endif

   enddo

   ! close file1
   close(iu1,status='keep')

   ! close file2
   close(iu2,status='keep')

   ! info
   ic = ic - 1
   sdev = dsqrt(sdev/dfloat(ic))
   call print_err_info(ic,thresh,errmax,dval1,dval2,icmax,sdev)

   ! check result
   if(errmax > thresh)then
     write(6,'(a)')
     write(6,'(a)') ' ERROR in compare_files routine: Examined files differ'
     write(0,'(a)') ' ERROR in compare_files routine: Examined files differ'
     stop
   endif

   end subroutine compare_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

   subroutine compare_text_files(iunit, file1, file2)
   implicit none

   ! External variables
   character(*), intent(in) :: file1, file2
   integer,      intent(in) :: iunit

   ! Internal variables
   character*(250) :: line1, line2
   integer :: ic, iu1, iu2, ier, ier1, ier2

   iu1 = iunit
   iu2 = iunit + 1

   ! open file1
   open(iu1,file=file1,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_text_files: ERROR while opening file ',file1
     write(0,*) 'compare_text_files: ERROR while opening file ',file1
     stop
   endif
   rewind iu1

   ! open file2
   open(iu2,file=file2,status='old',form='formatted',access='sequential',iostat=ier)
   if(ier /= 0)then
     write(6,*) 'compare_text_files: ERROR while opening file ',file2
     write(0,*) 'compare_text_files: ERROR while opening file ',file2
     stop
   endif
   rewind iu2

   ! compare two files
   ic = 0
   do

     ic = ic + 1
     read(iu1,'(a)',iostat=ier1) line1
     read(iu2,'(a)',iostat=ier2) line2
     if( (ier1 == -1) .and. (ier2 == -1) ) exit
     if( (ier1 == -1) .or.  (ier2 == -1) )then
       write(6,*) ' ERROR: Examined files have different length '
       write(0,*) ' ERROR: Examined files have different length '
       stop
     endif
     if( (ier1 /=  0) .or. (ier2 /=  0) )then
       write(6,*) 'compare_text_files: ERROR reading one of the files ',file1,' ',file2
       write(0,*) 'compare_text_files: ERROR reading one of the files ',file1,' ',file2
       write(6,*) 'ier1, ier2: ',ier1, ier2
       stop
     endif

     if( line1 /= line2 )then
       write(6,*) ' ERROR in compare_text_files routine'
       write(6,*) '       files (',file1,' and ',file2,') differ on line ',ic
       write(6,'(a)') strip_blank(line1)
       write(6,'(a)') strip_blank(line2)
       write(0,*) ' ERROR: files differ (see output for more details) '
       stop
     endif

   enddo

   ! close file1
   close(iu1,status='keep')

   ! close file2
   close(iu2,status='keep')

   end subroutine compare_text_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

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

   end subroutine matrix_to_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

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
         stop
      end if
      rewind(file_unit)

      read(file_unit, *) n
      if (n /= l*l) then
         print *, 'error: matrix dimensions do not match in file '//file_name
         stop
      end if

      do i = 1, l
         do j = 1, l
            read(file_unit, *) matrix(j, i)
         end do
      end do

      close(file_unit, status = 'keep')

   end subroutine matrix_from_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function strip_blank(string) result(str)
  implicit none
  character*(*), intent(in) :: string
  character(len=str_size(string)) :: str
  integer :: i
  do i=1,len(str)
    str(i:i) = string(i:i)
  enddo
  end function strip_blank

  pure function str_size(string) result(is)
  implicit none
  character*(*), intent(in) :: string
  integer :: i, is
  do i=len(string),1,-1
    if(string(i:i) /= ' ') exit
  enddo
  is = i
  end function str_size

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine print_real_err_info(ndim,thresh,errmax,dval1,dval2,icmax,sdev)
  implicit none

  integer, intent(in) :: ndim, icmax
  real*8,  intent(in) :: thresh, errmax, dval1, dval2, sdev

  write(6,'(a)')
  write(6,'(a)') ' --------------------------------'
  write(6,'(a)') ' Compare two sets of real numbers'
  write(6,'(a)') '   a --> reference  set'
  write(6,'(a)') '   b --> calculated set'
  write(6,'(a,d12.4)')  ' Threshold  = ',thresh
  write(6,'(a,i12)')    ' Number of compared elements  = ',ndim
  write(6,'(a,d12.4)')  ' The largest difference (max|a-b|) = ',errmax
  write(6,'(a,d25.16)') ' a = ',dval1
  write(6,'(a,d25.16)') ' b = ',dval2
  write(6,'(a,i12)')    ' Position of these elements = ',icmax
  write(6,'(a,d12.4)')  ' Standard deviation: sqrt( (a-b)^2 / n ) = ',dsqrt(sdev/dfloat(ndim))

  end subroutine print_real_err_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine print_complex_err_info(ndim,thresh,errmax,dval1,dval2,icmax,sdev)
  implicit none

  integer,    intent(in) :: ndim, icmax
  real*8,     intent(in) :: thresh, errmax, sdev
  complex*16, intent(in) :: dval1, dval2

  write(6,'(a)')
  write(6,'(a)') ' --------------------------------'
  write(6,'(a)') ' Compare two sets of real numbers'
  write(6,'(a)') '   a --> reference  set'
  write(6,'(a)') '   b --> calculated set'
  write(6,'(a,d12.4)')  ' Threshold  = ',thresh
  write(6,'(a,i12)')    ' Number of compared elements  = ',ndim
  write(6,'(a,d12.4)')  ' The largest difference (max|a-b|) = ',errmax
  write(6,'(a,2d25.16)')' a = ',dval1
  write(6,'(a,2d25.16)')' b = ',dval2
  write(6,'(a,i12)')    ' Position of these elements = ',icmax
  write(6,'(a,d12.4)')  ' Standard deviation: sqrt( (a-b)^2 / n ) = ',dsqrt(sdev/dfloat(ndim))

  end subroutine print_complex_err_info

end module
