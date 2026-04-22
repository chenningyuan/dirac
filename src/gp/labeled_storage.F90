
! Interface to hdf5 functionality to read/write data with a label. 
! Contains also the simple label writer used for DIRACs Fortran unformatted files (for backward compatibility).
! Contains also a simple unix-files based alternative in case hdf5 is not available.

! Written by Lucas Visscher, spring 2021, one of the Covid years

module labeled_storage

  implicit none

  private

  public lab_init
  public lab_close
  public lab_read
  public lab_write
  public lab_query
  public group_from_label

! Define type to hold information for an active file
  type, public:: file_info_t
     integer:: type=1                     ! which type of file type we use (old dirac(1), hdf5(2), hdf5-fallback(3), ..)
     integer:: status=-1                  ! -1: never used, 0: closed, 1: open
     integer:: unit=0                     ! unit used (for Fortran type writes)
     integer(kind=8):: id=0               ! id (for hdf5 type writes)
     logical:: keep_open=.false.          ! default is to close hdf5 files after each operation; opt in only with explicit lab_close
     character(len=:),allocatable :: name ! file name
  end type file_info_t

  type(file_info_t), save, public   :: dirac_datafiles(3) ! keep a limited number of file_info's here for convenience.

  interface lab_init
      module procedure init_general
  end interface

  interface lab_close
      module procedure close_general
  end interface

  interface lab_read
      module procedure read_general
      module procedure read_single_real
      module procedure read_3d_real
      module procedure read_single_int
      module procedure read_string_array
      module procedure read_dirac
  end interface

  interface lab_write
      module procedure write_general
      module procedure write_single_real
      module procedure write_single_int
      module procedure write_string_array
      module procedure write_dirac
  end interface

  interface lab_query
      module procedure query_general
  end interface

  interface group_from_label
      module procedure group_from_label_general
  end interface

  contains

!!!!!! Definition of the general init / read / write / query routines  !!!!!!!!

    subroutine init_general (file_info,top_groups)
      type(file_info_t),intent(inout)        :: file_info
      character(len=80), intent(in)          :: top_groups(:)

      select case (file_info%type)
      case (2)
         call init_hdf5(file_info,top_groups)
      case (3)
         call init_nohdf5(file_info,top_groups)
      end select

    end subroutine init_general

    subroutine close_general (file_info)
      type(file_info_t),intent(inout)        :: file_info

      select case (file_info%type)
      case (2)
         call close_hdf5(file_info)
      case default
         file_info%status = 0
      end select

    end subroutine close_general

    subroutine read_general (file_info,label,rdata,idata,sdata)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)
      character*(*),intent(out), optional    :: sdata

      integer :: ndata

      ndata = count_data_args(present(rdata),present(idata),present(sdata))
      if (ndata /= 1) then
         write (*,*) " ERROR: lab_read needs exactly one data argument for label ",trim(label)
         error stop ' error detected by read_general '
      end if

      select case (file_info%type)
      case (1)
         ! We use DIRAC files (FORTRAN unformatted sequential)
         if (present(rdata)) then
            call read_dirac(file_info%unit,label,rdata=rdata)
         elseif (present(idata)) then
            call read_dirac(file_info%unit,label,idata=idata)
         elseif (present(sdata)) then
            error stop "reading strings is not supported in DIRAC format"
         end if
      case (2)
         ! We use HDF5 files
         if (present(rdata)) then
            call read_hdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call read_hdf5(file_info,label,idata=idata)
         elseif (present(sdata)) then
            call read_hdf5(file_info,label,sdata=sdata)
         end if
      case (3)
         ! We use a UNIX fallback for HDF5 files
         if (present(rdata)) then
            call read_nohdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call read_nohdf5(file_info,label,idata=idata)
         elseif (present(sdata)) then
            call read_nohdf5(file_info,label,sdata=sdata)
         end if
      end select

    end subroutine read_general

    subroutine write_general(file_info,label,rdata,idata,sdata)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)
      character*(*),intent(in), optional     :: sdata

      integer :: ndata

      ndata = count_data_args(present(rdata),present(idata),present(sdata))
      if (ndata > 1) then
         write (*,*) " ERROR: lab_write accepts at most one data argument for label ",trim(label)
         error stop ' error detected by write_general '
      end if
      if (ndata == 0 .and. file_info%type == 2) then
         write (*,*) " ERROR: lab_write needs a data argument for HDF5 label ",trim(label)
         error stop ' error detected by write_general '
      end if

      select case (file_info%type)
      case (1)
         ! We use DIRAC files (FORTRAN unformatted sequential)
         if (present(rdata)) then
            call write_dirac(file_info%unit,label,rdata=rdata)
         elseif (present(idata)) then
            call write_dirac(file_info%unit,label,idata=idata)
         elseif (present(sdata)) then
            error stop "writing strings is not supported in DIRAC format"
         end if
      case (2)
         ! We use HDF5 files
         if (present(rdata)) then
            call write_hdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call write_hdf5(file_info,label,idata=idata)
         elseif (present(sdata)) then
            call write_hdf5(file_info,label,sdata=sdata)
          end if
      case (3)
         ! We use a UNIX fallback for HDF5 files
         if (present(rdata)) then
            call write_nohdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call write_nohdf5(file_info,label,idata=idata)
         elseif (present(sdata)) then
            call write_nohdf5(file_info,label,sdata=sdata)
          end if
      end select
    end subroutine write_general

    subroutine query_general (file_info,label,exist,data_size)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      logical,intent(out), optional          :: exist
      integer,intent(out), optional          :: data_size

      integer :: ierr
      logical :: tobe

      select case (file_info%type)
      case (1)
         ! We use DIRAC files (FORTRAN unformatted sequential)
         if (present(exist)) then
            error stop "query functionality is not supported in DIRAC format"
         end if
      case (2)
         ! We use HDF5 files
         if (present(data_size)) then
            call query_hdf5(file_info,label,exist=tobe,data_size=data_size)
         else
            call query_hdf5(file_info,label,exist=tobe)
         end if
         if (present(exist))     exist = tobe
      case (3)
         ! We use a UNIX fallback for HDF5 files
         if (present(data_size)) then
            call query_nohdf5(file_info,label,exist=tobe,data_size=data_size)
         else
            call query_nohdf5(file_info,label,exist=tobe)
         end if
         if (present(exist))     exist = tobe
      end select

    end subroutine query_general

    function group_from_label_general(label) result(group)

       character(len=*)            :: label
       character, parameter        :: delimiter='/'
       character(:),allocatable    :: group
       integer                     :: i

       ! find rightmost delimiter
       i = scan(label, delimiter, .true.)

       ! the first part is the group
       if (i > 1) then
          group = label(1:i-1)
       elseif (i == 1) then
          group = label ! This is one of the top groups
       else
          group = 'nogroups' ! No / found, no hierarchical structure
       end if

     end function group_from_label_general

    integer function count_data_args(has_rdata,has_idata,has_sdata) result(ndata)
      logical,intent(in) :: has_rdata, has_idata, has_sdata

      ndata = 0
      if (has_rdata) ndata = ndata + 1
      if (has_idata) ndata = ndata + 1
      if (has_sdata) ndata = ndata + 1

    end function count_data_args

!!!!!! Wrapper for the general read and write routines  !!!!!!!!

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

    subroutine read_3d_real (file_info,label,rdata)
      ! Wrapper around read_general to allow reading of 3d arrays
      use, intrinsic :: ISO_C_BINDING
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out),target             :: rdata(:,:,:)
      real(8),            pointer            :: rdata1d(:)

      nullify(rdata1d)
      call C_F_POINTER (C_LOC(rdata), rdata1d, [size(rdata)])
      !write(6,*) '#0 read_3d     '; call flush(6)
      call read_general (file_info,label,rdata=rdata1d)
      if(associated(rdata1d)) nullify(rdata1d)

    end subroutine read_3d_real

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

    subroutine write_string_array (file_info,label,sarray,slen)
      ! Wrapper around write_general to allow writing of string arrays
      ! Fixed length strings only (variable length string arrays are non-existent in Fortran)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      integer,intent(in)                     :: slen
      character(len=slen),intent(in)         :: sarray(:)
      character(len=:), allocatable          :: sdata
      integer                                :: i, sdim

      select case (file_info%type)
      case (1)
         error stop "writing strings is not supported in DIRAC format"
      case (2,3)
         ! There might be better ways to flatten an array of strings in Fortran, 
         ! but I could not find one and wasted already enough time on this.
         sdim = size(sarray)
         sdata = sarray(1)
         do i = 2, sdim
            sdata = sdata // sarray(i)
         end do
         if (file_info%type==2) then
             call write_hdf5(file_info,label,sdata=sdata,slen=slen,sdim=sdim)
         elseif (file_info%type==3) then 
             call write_nohdf5(file_info,label,sdata=sdata,slen=slen,sdim=sdim)
         end if
      end select

    end subroutine write_string_array

    subroutine read_string_array (file_info,label,sarray,slen)
      ! Wrapper around read_hdf5 to allow reading of string arrays
      ! Fixed length strings only (variable length string arrays are non-existent in Fortran)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      integer,intent(in)                     :: slen
      character(len=slen),intent(out)        :: sarray(:)
      character(len=:), allocatable          :: sdata
      integer                                :: i,j, sdim

      select case (file_info%type)
      case (1)
         error stop "reading strings is not supported in DIRAC format"
      case (2,3)
         if (file_info%type==2) then
             call read_hdf5(file_info,label,sdata=sdata)
         elseif (file_info%type==3) then 
             call read_nohdf5(file_info,label,sdata=sdata)
         end if
         sdim = len(sdata)/slen
         do i = 1, sdim
            sarray(i) = sdata((i-1)*slen+1:i*slen)
         end do
      end select

    end subroutine read_string_array

!!!!!! Selfcontained implementation of the old DIRAC format, used for backwards compatibility !!!!!!

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

!!!!!! Wrapper around the HDF5 interface mh5 !!!!!!

    subroutine init_hdf5 (file_info,top_groups)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character(len=80), intent(in)          :: top_groups(:)
      integer(kind=8)                        :: group_id
      integer                                :: i

      !  check first whether we can actually use mh5
      if (.not.mh5_is_available()) then
         print*," HDF5 interface is not available, using fallback Unix implementation"
         call toggle_fallback(file_info)
         call init_nohdf5(file_info,top_groups)
         return
      end if

      ! create the hdf5 file
      file_info%id = mh5_create_file(file_info%name)
      file_info%status = 1

      ! create the top level groups
      do  i = 1, size(top_groups)
         group_id = mh5_create_group(file_info%id,top_groups(i))
         call  mh5_close_group(group_id)
      end do

    end subroutine init_hdf5

    subroutine read_hdf5 (file_info,label,rdata,idata,sdata)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)
      character*(*),intent(out), optional    :: sdata

      logical :: exist
      integer :: ierr, data_size
      integer(kind=8), allocatable :: i8buf(:)

      !  check first whether we can actually use mh5
      inquire (file=file_info%name,exist=exist)
      if (.not.mh5_is_available().or..not.exist) then
         ! check whether the fallback can be used instead
         call toggle_fallback(file_info)
         inquire(file=file_info%name//'/.initialized',exist=exist)
         if (.not.exist) then
           ! can not be used, toggle back and stop
           if (mh5_is_available()) call toggle_fallback(file_info)
           write (*,*) " ERROR: trying to read data set ",trim(label)," from non-existing file ",trim(file_info%name)
           error stop ' error detected by read_hdf5 '
         end if
         if (present(rdata)) then
            call read_nohdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call read_nohdf5(file_info,label,idata=idata)
         elseif (present(sdata)) then
            call read_nohdf5(file_info,label,sdata=sdata)
         end if
         return
      end if

      if (file_info%status == 0) then
         file_info%id = mh5_open_file_rw(file_info%name)
         file_info%status = 1
      end if

      call query_hdf5 (file_info,label,exist,data_size)
      if (.not.exist) then
        write (*,*) " ERROR: trying to read non-existent data set ",trim(label)," from ",trim(file_info%name)
        error stop ' error detected by read_hdf5 '
      end if

      if (present(rdata)) then
         if (data_size > size(rdata)) then
           write (*,*) " ERROR: array too small for data set ",trim(label)," from ",trim(file_info%name) 
           write (*,'(2(A,I8/))') " Required ",data_size," Available ",size(rdata) 
           error stop ' error detected by read_hdf5 '
         end if
         call mh5_fetch_dset(file_info%id,label,rdata)
      elseif (present(idata)) then
         if (data_size > size(idata)) then
           write (*,*) " ERROR: array too small for data set ",trim(label)," from ",trim(file_info%name) 
           write (*,'(2(A,I8/))') " Required ",data_size," Available ",size(idata) 
           error stop ' error detected by read_hdf5 '
         end if
         allocate (i8buf(size(idata)))
         i8buf = 0
         call mh5_fetch_dset(file_info%id,label,i8buf)
         idata = int(i8buf)
         deallocate (i8buf)
      elseif (present(sdata)) then
         call mh5_fetch_dset(file_info%id,label,sdata)
      end if

      if (.not.file_info%keep_open) call close_hdf5(file_info)

    end subroutine read_hdf5

    subroutine write_hdf5 (file_info,label,rdata,idata,sdata,slen,sdim)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)
      character*(*),intent(in), optional     :: sdata
      integer,intent(in),optional            :: slen, sdim

      integer, parameter                     :: MAX_RANK=7
      integer(kind=8)                        :: dims(MAX_RANK),newdim
      integer :: ierr
      integer(kind=8) :: file_id, dset_id
      logical(kind=8) :: dyn
      logical :: group_exists, dset_log

      !  check first whether we can actually use mh5
      if (.not.mh5_is_available()) then
         ! use the fallback
         call toggle_fallback(file_info)
         if (present(rdata)) then
            call write_nohdf5(file_info,label,rdata=rdata)
         elseif (present(idata)) then
            call write_nohdf5(file_info,label,idata=idata)
         elseif (present(sdata).and.present(slen).and.present(sdim)) then
            call write_nohdf5(file_info,label,sdata=sdata,slen=slen,sdim=sdim)
         elseif (present(sdata)) then
            call write_nohdf5(file_info,label,sdata=sdata)
         end if
         return
      end if

      if (file_info%status == -1) then
         file_info%id = mh5_create_file(file_info%name)
         file_info%status = 1
      elseif (file_info%status == 0) then
         file_info%id = mh5_open_file_rw(file_info%name)
         file_info%status = 1
      end if

      file_id = file_info%id

      group_exists = create_group_hdf5(file_id,group_from_label(label)) ! check whether group is present and create it if not
      dset_log = mh5_exists_dset(file_id, label)                         ! check whether we already have data

      if (dset_log) then
          dset_id = mh5_open_dset(file_id, label)
          call mh5_get_dset_dims(dset_id,dims)
          ! if we are writing an array, we need to check whether a resize of the existing dataset is required
          if (present(rdata)) then
             newdim = size(rdata,kind=8)
          elseif (present(idata)) then
             newdim = size(idata,kind=8)
          elseif (present(sdata).and.present(slen).and.present(sdim)) then
             newdim = int8(sdim)
          else
             newdim = int8(1) ! writing a scalar
          end if
          if (newdim > 1) then
             ! check and resize the data set if size has changed
             if (newdim /= dims(1)) call mh5_resize_dset(dset_id,[newdim])
          end if
      end if

      if (present(rdata)) then
         dyn = size(rdata) > 1 ! Do not create resizable (chunked) arrays for single numbers
         if (.not.dset_log) dset_id = mh5_create_dset_real(file_id, label, int8(1), [size(rdata,kind=8)],dyn)
         call mh5_put_dset(dset_id, rdata)
         if (file_info%keep_open) call mh5_close_dset(dset_id)
      elseif (present(idata)) then
         dyn = size(idata) > 1 ! Do not create resizable (chunked) arrays for single numbers
         if (.not.dset_log) dset_id = mh5_create_dset_int(file_id, label, int8(1), [size(idata,kind=8)],dyn)
         call mh5_put_dset(dset_id, int8(idata))
         if (file_info%keep_open) call mh5_close_dset(dset_id)
      elseif (present(sdata).and.present(slen).and.present(sdim)) then
         dyn = .true.
         if (.not.dset_log) dset_id = mh5_create_dset_str(file_id, label, int8(1), [int8(sdim)], int8(slen),dyn)
         call mh5_put_dset(dset_id, sdata)
         if (file_info%keep_open) call mh5_close_dset(dset_id)
      elseif (present(sdata)) then
         if (.not.dset_log) dset_id = mh5_create_dset_str(file_id, label, int8(len(sdata)))
         call mh5_put_dset(dset_id, sdata)
         if (file_info%keep_open) call mh5_close_dset(dset_id)
      end if

      if (.not.file_info%keep_open) call close_hdf5(file_info)

    end subroutine write_hdf5

    subroutine close_hdf5 (file_info)
      use mh5
      type(file_info_t),intent(inout)        :: file_info

      if (file_info%status == 1) then
         call mh5_close_file(file_info%id)
         file_info%status = 0
         file_info%id = 0
      end if

    end subroutine close_hdf5

    subroutine query_hdf5 (file_info,label,exist,data_size)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      logical,intent(out)                    :: exist
      integer, parameter                     :: MAX_RANK=7
      integer(kind=8)                        :: group_id,dset_id,dims(MAX_RANK)
      integer,intent(out), optional          :: data_size
      character*(:),allocatable              :: group

      !  check first whether we can actually use mh5
      inquire (file=file_info%name,exist=exist)
      if (.not.mh5_is_available().or..not.exist) then
         call toggle_fallback(file_info)
         if (present(data_size)) then
            call query_nohdf5(file_info,label,exist,data_size)
         else
            call query_nohdf5(file_info,label,exist)
         end if
         ! restore original state and go back
         call toggle_fallback(file_info)
         return
      end if

      if (file_info%status == -1) then
         ! unknown file, check whether it exists and do not open (and create) if it does not
         inquire(file=file_info%name,exist=exist)
         if (.not.exist) return
      end if

      if (file_info%status /= 1) then
          file_info%id = mh5_open_file_rw(file_info%name)
          if (file_info%id < 0) then
             ! invalid identifier, error when opening this file
             print*,"  - Could not open file ",file_info%name
             exist = .false.
             return
          end if
          file_info%status = 1
      end if

      group = group_from_label(label)
      if (group=='nogroups') group = label

      if (group == label) then
         ! a top group that can be checked immediately
         exist = mh5_exists_dset(file_info%id,label)
      else
         ! check the existence of the group first
         exist = mh5_exists_group(file_info%id,group)
         ! check the existence of the data set itself
         if (exist) exist = mh5_exists_dset(file_info%id,label)
      end if

      if (present(data_size)) then
         if (exist) then
            dset_id = mh5_open_dset(file_info%id, label)
            call mh5_get_dset_dims(dset_id,dims)
            data_size = dims(1)
            call mh5_close_dset(dset_id)
         else
            data_size = 0
         end if
      end if

    end subroutine query_hdf5

    recursive function create_group_hdf5 (file_id,label) result (group_exists)

     ! Check whether the group does exist and create it if it does not
     ! All groups should be subgroups of an existing top group
     use mh5
     character(len=*),intent(in) :: label
     logical                     :: group_exists,parent_exists
     character(:),allocatable    :: parent_label
     integer(kind=8)             :: file_id,group_id

     ! only in case there is no hierarchical structure we need not create any groups
     group_exists = label == 'nogroups'
     if (group_exists) return

     ! otherwise we keep unpeeling the groups till we hit a parent that already exists
     group_exists = mh5_exists_group(file_id,label)
     if (group_exists) return

     parent_label  = group_from_label(label)
     parent_exists = mh5_exists_group(file_id,parent_label)

     if (parent_exists) then
        group_id = mh5_create_group(file_id,label)
        call mh5_close_group(group_id)
        group_exists = .true.
     else
        group_exists = create_group_hdf5  (file_id,parent_label)
        group_id = mh5_create_group (file_id,label)
        call mh5_close_group(group_id)
     end if

    end function create_group_hdf5

!!!!!! Fallback if HDF5 is not available !!!!!!

    subroutine init_nohdf5 (file_info,top_groups)
      type(file_info_t),intent(inout)        :: file_info
      character(len=80), intent(in)          :: top_groups(:)
      integer                                :: i, ierr
      character*(:), allocatable             :: topdir,subdir

      topdir = './'//trim(file_info%name)
      call execute_command_line ('mkdir -p '//topdir,exitstat=ierr)
      if (ierr.ne.0) then
         print*, " ERROR: could not create directory ",topdir
         error stop " error in init_nohdf5 "
      end if
      call execute_command_line ('date >'//topdir//'/.initialized')
      file_info%status = 1

      ! create the top level groups
      do  i = 1, size(top_groups)
          subdir = topdir//trim(top_groups(i))
          call execute_command_line ('mkdir -p '//subdir,exitstat=ierr)
          if (ierr.ne.0) then
             print*, " ERROR: could not create directory ",subdir
             error stop " error in create_group_nohdf5 "
          end if
          call execute_command_line ('date >'//subdir//'/.initialized')
      end do

    end subroutine init_nohdf5

    subroutine read_nohdf5 (file_info,label,rdata,idata,sdata)
      use mh5
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)
      character*(*),intent(out), optional    :: sdata

      integer, parameter          :: file_unit=66623
      logical                     :: file_exists, group_exists
      character*(:), allocatable  :: data_file
      integer                     :: data_size
      integer(8)                  :: data_size8,data_len8
      character(len=4)            :: data_type

      call query_nohdf5 (file_info,label,file_exists,data_size)
      if (.not.file_exists) then
        write (*,*) " ERROR: trying to read non-existent data set ",trim(label)," from ",trim(file_info%name)
        error stop ' error detected by read_hdf5 '
      end if

      data_file = './'//trim(file_info%name)//'/'//label
      open (file_unit,file=data_file,form='unformatted')
      read (file_unit) data_type,data_size8
      data_size = int(data_size8)

      if (present(rdata)) then
         if (data_size > size(rdata)) then
           write (*,*) " ERROR: array too small for data set ",trim(label)," from ",trim(file_info%name) 
           write (*,'(2(A,I8/))') " Required ",data_size," Available ",size(rdata) 
           error stop ' error detected by read_hdf5 '
         end if
         read (file_unit) rdata(1:data_size)
      elseif (present(idata)) then
         if (data_size > size(idata)) then
           write (*,*) " ERROR: array too small for data set ",trim(label)," from ",trim(file_info%name) 
           write (*,'(2(A,I8/))') " Required ",data_size," Available ",size(idata) 
           error stop ' error detected by read_hdf5 '
         end if
         read (file_unit) idata(1:data_size)
      elseif (present(sdata)) then
         read (file_unit) sdata
      end if

      close (file_unit,status='keep')

    end subroutine read_nohdf5

    subroutine write_nohdf5 (file_info,label,rdata,idata,sdata,slen,sdim)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)
      character*(*),intent(in), optional     :: sdata
      integer,intent(in),optional            :: slen, sdim

      integer, parameter          :: file_unit=66623
      logical                     :: file_exists, group_exists
      character*(:), allocatable  :: data_file

      if (file_info%status == -1) then
         ! the file has not yet been opened in this run, create the top directory (if needed)
         file_exists = create_group_nohdf5(file_info%name,'')
         file_info%status = 1
      end if

      group_exists = create_group_nohdf5(file_info%name,group_from_label(label)) ! check whether group is present and create it if not
      data_file = './'//trim(file_info%name)//'/'//label
      open (file_unit,file=data_file,form='unformatted')

      if (present(rdata)) then
         write (file_unit) 'real',size(rdata,kind=8),int8(1)
         write (file_unit) rdata
      elseif (present(idata)) then
        if (kind(idata)==4) write (file_unit) 'int4',size(idata,kind=8),int8(1)
        if (kind(idata)==8) write (file_unit) 'int8',size(idata,kind=8),int8(1)
         write (file_unit) idata
      elseif (present(sdata).and.present(slen).and.present(sdim)) then
         write (file_unit) 'strn',int8(slen),int8(sdim)
         write (file_unit) sdata
      elseif (present(sdata)) then
         write (file_unit) 'str1',int8(len(sdata)),int8(1)
         write (file_unit) sdata
      end if

      close (file_unit,status='keep')
      file_info%status = 0

    end subroutine write_nohdf5

    subroutine query_nohdf5 (file_info,label,exist,data_size)
      type(file_info_t),intent(inout)        :: file_info
      character*(*),intent(in)               :: label
      logical,intent(out)                    :: exist
      integer,intent(out), optional          :: data_size

      integer, parameter          :: file_unit=66623
      logical                     :: file_exists
      character*(:), allocatable  :: data_file
      integer(8)                  :: data_size8
      character(len=4)            :: data_type

      data_file = './'//trim(file_info%name)//'/'//label
      inquire (file=data_file,exist=exist)

      ! could be an inquiry about a directory, so check this possibility as well
      if (.not.exist) inquire(file=data_file//'/.initialized',exist=exist)

      if (present(data_size)) then
         if (exist) then
            open (file_unit,file=data_file,form='unformatted')
            read (file_unit) data_type,data_size8
            close (file_unit,status='keep')
            data_size = data_size8
         else
            data_size = 0
         end if
      end if

    end subroutine query_nohdf5

    function exists_group_nohdf5 (file_name,label) result (group_exists)
     ! Check whether a group does exist
     character(len=*),intent(in) :: file_name,label
     logical                     :: group_exists
     character*(:), allocatable  :: subdir
     if (len_trim(label) == 0) then
        subdir = './'//trim(file_name)
     else
        subdir = './'//trim(file_name)//'/'//trim(label)
     end if
     inquire (file=subdir//'/.initialized',exist=group_exists)
    end function exists_group_nohdf5

    recursive function create_group_nohdf5 (file_name,label) result (group_exists)
     ! Check whether the group does exist and create it if it does not
     character(len=*),intent(in) :: file_name,label
     character(:),allocatable    :: parent_label
     logical                     :: group_exists,parent_exists
     integer                     :: ierr
     character*(:), allocatable  :: subdir

     ! only in case there is no hierarchical structure we need not create any groups
     group_exists = label == 'nogroups'
     if (group_exists) return

     if (len_trim(label) == 0) then
        subdir = './'//trim(file_name)
        call execute_command_line ('mkdir -p '//subdir,exitstat=ierr)
        if (ierr.ne.0) then
           print*, " ERROR: could not create directory ",trim(file_name)
           error stop " error in create_group_nohdf5 "
        end if
        call execute_command_line ('date >'//subdir//'/.initialized')
        group_exists = .true.
        return
     end if

     group_exists = exists_group_nohdf5(file_name,label)
     if (group_exists) return

     parent_label  = group_from_label(label)
     if (parent_label == 'nogroups') then
        parent_exists = create_group_nohdf5(file_name,'')
     else
        parent_exists = exists_group_nohdf5(file_name,parent_label)
     end if

     subdir = './'//trim(file_name)//'/'//trim(label)
     if (parent_exists) then
        call execute_command_line ('mkdir -p '//subdir,exitstat=ierr)
        if (ierr.ne.0) then
           print*, " ERROR: could not create directory ",trim(file_name)
           error stop " error in create_group_nohdf5 "
        end if
        call execute_command_line ('date >'//subdir//'/.initialized')
        group_exists = .true.
     else
        group_exists = create_group_nohdf5(file_name,parent_label)
        call execute_command_line ('mkdir -p '//subdir,exitstat=ierr)
        if (ierr.ne.0) then
           print*, " ERROR: could not create directory ",trim(file_name)
           error stop " error in create_group_nohdf5 "
        end if
        call execute_command_line ('date >'//subdir//'/.initialized')
        group_exists = .true.
     end if

    end function create_group_nohdf5

    subroutine toggle_fallback(targetfile)
         type(file_info_t), intent(inout)    :: targetfile
         integer                             :: len_name

         len_name = len_trim(targetfile%name)
         if (targetfile%type == 2 .and. targetfile%name(len_name-1:len_name) == 'h5') then
           targetfile%name = targetfile%name(1:len_name-2)//'noh5'
           targetfile%type = 3
         elseif (targetfile%type == 3 .and. targetfile%name(len_name-3:len_name) == 'noh5') then
           targetfile%name = targetfile%name(1:len_name-4)//'h5'
           targetfile%type = 2
         endif

    end subroutine toggle_fallback


end module labeled_storage
