
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
     integer(kind=8):: id=0               ! id (for hdf5 type writes)
     character(len=:),allocatable :: name ! file name
  end type file_info_t

  interface lab_read
      module procedure read_general
      module procedure read_single_real
      module procedure read_single_int
      module procedure read_dirac
  end interface

  interface lab_write
      module procedure write_general
      module procedure write_single_real
      module procedure write_single_int
      module procedure write_dirac
  end interface

  contains

!!!!!! Definition of the general read / write routines  !!!!!!!!

    subroutine read_general (file_info,label,rdata,idata)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)

      integer :: ierr

      select case (file_info%type)
      case (1)
         ! We use DIRAC files (FORTRAN unformatted sequential)
         if (present(rdata)) then
            call read_dirac(file_info%unit,label,rdata=rdata)
         elseif (present(idata)) then
            call read_dirac(file_info%unit,label,idata=idata)
         end if
      case (2)
         ! We use HDF5 files
         if (present(rdata)) then
            call read_hdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call read_hdf5(file_info,label,idata=idata)
         end if
      end select

    end subroutine read_general

    subroutine write_general(file_info,label,rdata,idata)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)

      integer :: ierr

      select case (file_info%type)
      case (1)
         ! We use DIRAC files (FORTRAN unformatted sequential)
         if (present(rdata)) then
            call write_dirac(file_info%unit,label,rdata=rdata)
         elseif (present(idata)) then
            call write_dirac(file_info%unit,label,idata=idata)
         end if
      case (2)
         ! We use HDF5 files
         if (present(rdata)) then
            call write_hdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call write_hdf5(file_info,label,idata=idata)
          end if
      end select
    end subroutine write_general

    subroutine write_single_real (file_info,label,rdata)
      ! Wrapper around write_general to allow writing of single numbers
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in)                     :: rdata
      real(8)                                :: rdata_array(1)

      rdata_array(1) = rdata
      call write_general (file_info,label,rdata=rdata_array)

    end subroutine write_single_real

    subroutine read_single_real (file_info,label,rdata)
      ! Wrapper around read_general to allow reading of single numbers
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out)                    :: rdata
      real(8)                                :: rdata_array(1)

      call read_general (file_info,label,rdata=rdata_array)
      rdata = rdata_array(1)

    end subroutine read_single_real

    subroutine write_single_int (file_info,label,idata)
      ! Wrapper around write_general to allow writing of single numbers
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      integer,intent(in)                     :: idata
      integer                                :: idata_array(1)

      idata_array(1) = idata
      call write_general (file_info,label,idata=idata_array)

    end subroutine write_single_int

    subroutine read_single_int (file_info,label,idata)
      ! Wrapper around read_general to allow reading of single numbers
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      integer,intent(out)                    :: idata
      integer                                :: idata_array(1)

      call read_general (file_info,label,idata=idata_array)
      idata = idata_array(1)

    end subroutine read_single_int

!!!!!! Selfcontained implementation of The DIRAC format, used for backwards compatibility !!!!!!

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

!!!!!! Wrapper around the HDF5 interface !!!!!!

    subroutine read_hdf5 (file_info,label,rdata,idata)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)

      integer :: ierr
      integer(kind=8), allocatable :: i8buf(:)

      if (file_info%status == 0) then
         file_info%id = mh5_open_file_rw(file_info%name)
         file_info%status = 1
      end if

      if (present(rdata)) then
         call mh5_fetch_dset(file_info%id,label,rdata)
      elseif (present(idata)) then
         allocate (i8buf(size(idata)))
         call mh5_fetch_dset(file_info%id,label,i8buf)
         idata = int(i8buf)
         deallocate (i8buf)
      end if

      call mh5_close_file(file_info%id)
      file_info%status = 0

    end subroutine read_hdf5

    subroutine write_hdf5 (file_info,label,rdata,idata)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)

      integer :: ierr
      integer(kind=8) :: file_id, dset_id, group_id

      if (file_info%status == -1) then
         file_info%id = mh5_create_file(file_info%name)
         file_info%status = 1
         group_id = mh5_open_group(file_info%id,'/result')
         group_id = mh5_open_group(file_info%id,'/result/wavefunctions')
         group_id = mh5_open_group(file_info%id,'/result/wavefunctions/scf')
         group_id = mh5_open_group(file_info%id,'/result/wavefunctions/scf/mobasis')
      elseif (file_info%status == 0) then
         file_info%id = mh5_open_file_rw(file_info%name)
         file_info%status = 1
      end if

      file_id = file_info%id

      if (present(rdata)) then
         dset_id = mh5_exists_dset(file_id, label)
         if (dset_id == 0) then
            dset_id = mh5_create_dset_real(file_id, label, int8(1), [size(rdata,kind=8)])
         else
            dset_id = mh5_open_dset(file_id, label)
         end if
         call mh5_put_dset(dset_id, rdata)
         call mh5_close_dset(dset_id)
      elseif (present(idata)) then
         dset_id = mh5_exists_dset(file_id, label)
         if (dset_id == 0) then
            dset_id = mh5_create_dset_int(file_id, label, int8(1), [size(idata,kind=8)])
         else
            dset_id = mh5_open_dset(file_id, label)
         end if
         call mh5_put_dset(dset_id, int8(idata))
         call mh5_close_dset(dset_id)
      end if

      call mh5_close_file(file_id)
      file_info%status = 0

    end subroutine write_hdf5

end module labeled_storage
