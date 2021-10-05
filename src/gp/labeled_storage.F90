
! Interface to read/write data with a label
! Contains the simple label writer used for DIRACs Fortran unformatted files
! DFCOEF and AOPROPER to allow for easy backwards compatibility

! Written by Lucas Visscher, spring 2021, one of the Covid years

module labeled_storage
  implicit none

  private

  public lab_read
  public lab_write

! Define type to hold information for an active file
  type, public:: file_info_t
     integer:: type=1                     ! Future use: which type of file type we use (hpf,kf, dirac, ..)
     integer:: status=-1                  ! -1: never used, 0: closed, 1: open
     integer:: unit=0                     ! unit used (for Fortran type writes)
     character(len=:),allocatable :: name ! file name
  end type file_info_t

  interface lab_read
      module procedure read_dirac
      module procedure read_general
  end interface

  interface lab_write
      module procedure write_dirac
      module procedure write_general
  end interface

  contains

    subroutine read_dirac (file_unit,label,rdata,idata)
      integer, intent(in)                    :: file_unit
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)

      integer :: ierr

      ! We only can read one type of data, check for real and integer.
      ! Planned extension is to allow for a call-back with a reader for a datatype defined elsewhere.
      ! At the moment we just position the file at the right position for such cases

      ierr = dirac_labsearch (file_unit,label)
      if (present(rdata)) then
         if (ierr==0) read (file_unit) rdata
      elseif (present(idata)) then
         if (ierr==0) read (file_unit) idata
      end if
    end subroutine read_dirac

    subroutine read_general (file_info,label,rdata,idata)
      type(file_info_t),intent(in)           :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)

      integer :: ierr

      if (file_info%type == 1) then
         ! Interfaces to file types different from dirac still to be implemented.
         if (present(rdata)) then
            call read_dirac(file_info%unit,label,rdata=rdata)
         elseif (present(idata)) then
            call read_dirac(file_info%unit,label,idata=idata)
         end if
      end if
      
    end subroutine read_general

    subroutine write_dirac (file_unit,label,rdata,idata)
      integer, intent(in)                    :: file_unit
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)

      integer :: ierr

      ! We only can write one type of data, check for real and integer.
      ! Planned extension is to allow for a call-back with a reader for a datatype defined elsewhere.
      ! At the moment we just position the file at the right position for such cases

      call dirac_newlab (file_unit,label)
      if (present(rdata)) then
         write (file_unit) rdata
      elseif (present(idata)) then
         write (file_unit) idata
      end if
    end subroutine write_dirac

    subroutine write_general(file_info,label,rdata,idata)
      type(file_info_t),intent(in)           :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)

      integer :: ierr

      if (file_info%type == 1) then
         ! Interfaces to file types different from dirac still to be implemented.
         if (present(rdata)) then
            call write_dirac(file_info%unit,label,rdata=rdata)
         elseif (present(idata)) then
            call write_dirac(file_info%unit,label,idata=idata)
         end if
      end if
    end subroutine write_general

    integer function dirac_labsearch (file_unit,label)
      integer,intent(in)        :: file_unit
      character*(*), intent(in) :: label
      character(len=8)          :: stars,label_record(4)

      logical :: found

      stars = '********'
      found = .false.
      rewind (file_unit)
      do while (.not. found)
         read (file_unit,end=1,err=2) label_record
         if (label_record(1) .ne. stars) go to 2
         if (label_record(4) .eq. label) found = .true.
         if ((.not.found).or.label_record(1).eq.'INFO    ') read(file_unit)
      end do
     
      dirac_labsearch = 0 ! label found, file is now positioned for a read
      return
   1  dirac_labsearch = 1 ! label is not found
      return
   2  dirac_labsearch = 2 ! file is corrupted
      return

    end function dirac_labsearch
      
    subroutine dirac_newlab (file_unit,label)
      integer,intent(in)        :: file_unit
      character*(*), intent(in) :: label
      character(len=8)          :: dirac_label,time_stamp,date_stamp,stars
      character(len=10)         :: time10

      stars = '********'
      call date_and_time(date=date_stamp,time=time10)
      time_stamp = time10(1:8) ! we need to truncate to be compatible with the DIRAC formar, chop off the ms.
      write (file_unit) stars,date_stamp,time_stamp,label

    end subroutine dirac_newlab

    integer function get_unit()
    ! find a free unit number (go fortran!)
      implicit none
      integer :: unit
      logical :: isOpen

      integer, parameter :: MIN_UNIT_NUMBER = 10
      integer, parameter :: MAX_UNIT_NUMBER = 399

      do unit = MAX_UNIT_NUMBER, MIN_UNIT_NUMBER,-1
         inquire(unit = unit, opened = isOpen)
         if (.not. isOpen) then
            get_unit = unit
            return
         end if
      end do
    end function get_unit

end module labeled_storage
