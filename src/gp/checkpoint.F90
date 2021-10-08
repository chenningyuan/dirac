
! Module to work with a checkpoint file that contains all essential and
! some optional data. Only data defined in the schema can be read or written.
! For non-standard data one may define other (hdf5 or simple Fortran) files.

! Written by Lucas Visscher, fall 2021

module checkpoint

  use labeled_storage

  implicit none

  private

  ! Constants
  logical,parameter            :: VALID = .true.                   ! defined for convenience
  character(len=13), parameter :: CHECKPOINT_NAME="CHECKPOINT.h5"  ! name of the checkpoint file
  character(len=17), parameter :: SCHEMA_NAME="schema_labels.txt"  ! pathname of the file that contains the labels (created from DIRACschema.txt by a python function)

  ! The file object that is stored in global memory
  type(file_info_t), save      :: checkpoint_file

  ! These 4 functions should provide all functionality needed to write/read data defined in the schema.
  ! Note that the file is closed after each read/write, so no explicit close is needed here.

  public checkpoint_open
  public checkpoint_valid
  public checkpoint_write
  public checkpoint_read

  interface checkpoint_write
      module procedure checkpoint_write
      module procedure write_single_int
      module procedure write_single_real
  end interface

  interface checkpoint_read
      module procedure checkpoint_read
      module procedure read_single_int
      module procedure read_single_real
  end interface

  contains

!------- Public procedures ----

  integer function checkpoint_open()

      ! Open a new or existing file and initialize or check all groups defined in the schema
      ! After this step individual data sets can be written

      logical           :: restart
      integer,parameter :: SUCCESS=0
      integer,parameter :: FAILURE=1

      ! Status should be -1 (checkpoint_file never used) when calling this routine, check this
      if (checkpoint_file%status /= -1) then
          checkpoint_open = FAILURE
          return
      end if

      ! Proceed by initializing the file object as a hdf5 file with the correct name
      checkpoint_file%type = 2
      checkpoint_file%name = CHECKPOINT_NAME

      ! Inquire whether (a yet unopened) file with this name is present and attempt a restart if it is
      inquire (file=CHECKPOINT_NAME,exist=restart)
      if (restart) then
         print*," Opening existing HDF5 checkpoint file and checking its contents"
         if (checkpoint_valid(checkpoint_file)) then
            print*, " Checkpoint file is valid, will be used for restart"
            checkpoint_open = SUCCESS
            return
         else 
            print*, " Checkpoint file is not valid, will be overwritten"
            call unlink(CHECKPOINT_NAME)
            checkpoint_file%status = -1
         endif
      endif 

      ! If we come this far we are dealing with a pristine directory and need to create and init the hdf5 file
      call create_all_groups

      checkpoint_open = SUCCESS
      return
      
  end function checkpoint_open

  logical function checkpoint_valid(file_info)

      ! Check whether the file given as argument is a valid checkpoint file
      ! Criteria used for checking:
      ! - file should be readable (hdf5 formatted)
      ! - file contains all allowed groups (data inside may be empty)
      ! - file contains all required datasets (this data should be present)

      ! function is called at the end of a run to check the integrity of the file
      ! and can also be called at the beginning of a run to see if restart is possible

      type(file_info_t) :: file_info

      !TODO: implement this check, this is just a placeholder at the moment
      checkpoint_valid = VALID
      return

  end function checkpoint_valid

  subroutine checkpoint_write(label,rdata,idata,sdata)

      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)
      character*(*),intent(in), optional     :: sdata

      logical valid

      if (present(rdata)) then
         valid = check_label(label,'real')
         if (valid) then
             call lab_write(checkpoint_file,label,rdata=rdata)
         else
             print*,"WARNING: invalid real data ",label
         end if
      elseif (present(idata)) then
         valid =  check_label(label,'integer')
         if (valid) then
             call lab_write(checkpoint_file,label,idata=idata)
         else
             print*,"WARNING: invalid integer data ",label
         end if
      elseif (present(sdata)) then
         valid = check_label(label,'string')
         if (valid) then
             call lab_write(checkpoint_file,label,sdata=sdata)
         else
             print*,"WARNING: invalid string data ",label
         end if
      end if

  end subroutine checkpoint_write

  subroutine checkpoint_read(label,rdata,idata,sdata)

      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)
      character*(*),intent(out), optional    :: sdata

      !TODO: implement error handling
      if (present(rdata)) then
         call lab_read(checkpoint_file,label,rdata=rdata)
      elseif (present(idata)) then
         call lab_read(checkpoint_file,label,idata=idata)
      elseif (present(sdata)) then
         call lab_read(checkpoint_file,label,sdata=sdata)
      end if

  end subroutine checkpoint_read

  subroutine write_single_real (label,rdata)
      ! Wrapper to allow writing of single numbers
      character*(*),intent(in)               :: label
      real(8),intent(in)                     :: rdata
      real(8)                                :: rdata_array(1)

      rdata_array(1) = rdata
      call checkpoint_write (label,rdata=rdata_array)

  end subroutine write_single_real

  subroutine read_single_real (label,rdata)
      ! Wrapper to allow reading of single numbers
      character*(*),intent(in)               :: label
      real(8),intent(out)                    :: rdata
      real(8)                                :: rdata_array(1)

      call checkpoint_read (label,rdata=rdata_array)
      rdata = rdata_array(1)

  end subroutine read_single_real

  subroutine write_single_int (label,idata)
      ! Wrapper to allow writing of single numbers
      character*(*),intent(in)               :: label
      integer,intent(in)                     :: idata
      integer                                :: idata_array(1)

      idata_array(1) = idata
      call checkpoint_write (label,idata=idata_array)

  end subroutine write_single_int

  subroutine read_single_int (label,idata)
      ! Wrapper to allow reading of single numbers
      character*(*),intent(in)               :: label
      integer,intent(out)                    :: idata
      integer                                :: idata_array(1)

      call checkpoint_read (label,idata=idata_array)
      idata = idata_array(1)

  end subroutine read_single_int

!------- Private procedures ----

  subroutine create_all_groups

      ! Part of initialization procedure, hdf5 requires that a group is created
      ! before attempting to write data that is part of this group. Do this for
      ! all groups defined in the schema.

      use mh5
      integer :: unit
      logical :: composite_data = .true.
      character(len=:),allocatable:: fformat,label,type,rank,use,description 
      integer(kind=8) :: group_id

      !open(newunit=unit,form='formatted',file=SCHEMA_NAME)
      unit = 237
      open(unit,form='formatted',file=SCHEMA_NAME)
      read(unit,'(A)') fformat
      fformat='(A80,4A10,A60)'
      print*," read format from ",SCHEMA_NAME," :",fformat
      do while (composite_data)
         read(unit,fformat) label,type,rank,use,description
         print*,label,type,rank,use,description
         composite_data = index(type,'composite').ne.0
         if (composite_data) print*, "create group:",label
         if (composite_data) group_id = mh5_open_group(checkpoint_file%id,label)
      end do
      
      close(unit)

  end subroutine create_all_groups

  logical function check_label(mylabel,mytype)

      ! Check whether the label is known and refers to the right type of data

      character*(*), intent(in) :: mylabel, mytype
      logical :: found = .false.
      integer :: unit, iostatus=0
      character(len=:),allocatable:: format,label,type,rank,use,description 

      ! Check existence of mylabel
      ! Should be made more efficient than implemented below.
      open(newunit=unit,form='formatted',file=SCHEMA_NAME)
      read(unit,*) format
      format='(A80,4A10,A60)'
      do while (iostatus > -1)
         read(unit,format,iostat=iostatus) label,type,rank,use,description
         found = index(label,mylabel).ne.0
         print*,iostatus,found,label,type,rank,use,description
         if (found) exit
      end do
      close(unit)

      ! Check whether the type matches
      if (found) then
         if (index(type,mytype).ne.0) then
            check_label = VALID
            return
         else
            check_label = .not.VALID
            return
         end if
      else
         check_label = .not.VALID
         return
      end if

  end function check_label

end module
