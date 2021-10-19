
! Module to work with a checkpoint file that contains all essential and
! some optional data. Only data defined in the schema can be read or written.
! For non-standard data one may define other (hdf5 or simple Fortran) files.

! Written by Lucas Visscher, fall 2021

module checkpoint

  use labeled_storage

  implicit none

  private

  ! A type to store the schema
  type :: label_t
      character(len=:),allocatable:: label,type,rank,use,description
  end type label_t

  ! Constants
  logical,parameter            :: VALID = .true.                   ! defined for convenience
  character(len=13), parameter :: CHECKPOINT_NAME="CHECKPOINT.h5"  ! name of the checkpoint file
  character(len=17), parameter :: SCHEMA_NAME="schema_labels.txt"  ! pathname of the file that contains the labels (created from DIRACschema.txt by a python function)

  ! The file object that is stored in global memory
  type(file_info_t), save      :: checkpoint_file

  ! Keep also the schema in the modules (private) global memory
  type(label_t), allocatable, save :: labels(:)

  ! These 4 functions should provide all functionality needed to write/read data defined in the schema.

  public checkpoint_open
  public checkpoint_close
  public checkpoint_write
  public checkpoint_read

  interface checkpoint_write
      module procedure checkpoint_write
      module procedure write_single_int
      module procedure write_single_real
      module procedure write_string_array
  end interface

  interface checkpoint_read
      module procedure checkpoint_read
      module procedure read_single_int
      module procedure read_single_real
      module procedure read_string_array
  end interface

  contains

!------- Public procedures ----

  subroutine checkpoint_open

      ! Open a new or existing file and initialize or check data sets defined in the schema
      ! After this step individual data sets can be written (or read in case of restart)

      logical           :: restart

      ! Read the schema to define all allowed labels
      call get_labels

      ! Checkpoint_file should not be used before when calling this routine, check this
      if (checkpoint_file%status /= -1) then
          stop " Error in checkpoint_open: file is already in use"
      end if

      ! Proceed by initializing the file object as a hdf5 file with the correct name
      checkpoint_file%type = 2
      checkpoint_file%name = CHECKPOINT_NAME

      ! Inquire whether (a yet unopened) file with this name is present and attempt a restart if it is
      inquire (file=CHECKPOINT_NAME,exist=restart)
      if (restart) then
         print*,"   - Existing checkpoint file: check contents"
         if (checkpoint_valid()) then
            print*, "   - Checkpoint file is valid, will be used for restart"
            return
         else 
            print*, "   - Checkpoint file is not valid, will be overwritten"
            call unlink(CHECKPOINT_NAME)
            checkpoint_file%status = -1
         endif
      endif 

      ! If we come this far we are dealing with a pristine directory and need to create and init the hdf5 file
      call create_new_checkpoint
      print*, "   - New checkpoint file created"

      call checkpoint_write ('/result/execution/status',idata=0)

  end subroutine checkpoint_open

  subroutine checkpoint_close

      use mh5

      ! Close the checkpoint file after having checked that all required data is written

      logical           :: restart

      call checkpoint_write ('/result/execution/status',idata=2)

      if (checkpoint_valid()) then
            print*, "   - Checkpoint file contains all required data, can be used for restarts"
      else
            print*, "   - Checkpoint file does not contains all required data, can not be used for restarts"
      end if

      call mh5_close_file(checkpoint_file%id)
      checkpoint_file%status = 0
      print*, "   - Checkpoint file closed"

  end subroutine checkpoint_close

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
             print*,"WARNING: skipped writing invalid real data ",label
         end if
      elseif (present(idata)) then
         valid =  check_label(label,'integer')
         if (valid) then
             call lab_write(checkpoint_file,label,idata=idata)
         else
             print*,"WARNING: skipped writing invalid integer data ",label
         end if
      elseif (present(sdata)) then
         valid = check_label(label,'string')
         if (valid) then
             call lab_write(checkpoint_file,label,sdata=sdata)
         else
             print*,"WARNING: skipped writing invalid string data ",label
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

  subroutine write_string_array (label,sarray,slen)
      ! Wrapper to write string arrays
      character*(*),intent(in)               :: label
      integer                                :: slen
      character(len=slen),intent(in)         :: sarray(:)
      logical                                :: valid

      valid = check_label(label,'string')
      if (valid) then
         call lab_write(checkpoint_file,label,sarray=sarray,slen=slen)
      else
         print*,"WARNING: skipped writing invalid string array ",label
      end if

  end subroutine write_string_array

  subroutine read_string_array (label,sarray,slen)
      ! Wrapper to read string arrays
      character*(*),intent(in)               :: label
      integer                                :: slen
      character(len=slen),intent(out)        :: sarray(:)

      call lab_read(checkpoint_file,label,sarray=sarray,slen=slen)

  end subroutine read_string_array

!------- Private procedures ----

  logical function checkpoint_valid()

      use mh5

      ! Check whether the file given as argument is a valid checkpoint file
      ! Criteria used for checking:
      ! - file should be readable (hdf5 formatted)
      ! - file contains all required datasets (this data should be present)

      ! function is called at the end of a run to check the integrity of the file
      ! and can also be called at the beginning of a run to see if restart is possible

      integer         :: i
      integer(kind=8) :: group_id

      checkpoint_valid = VALID

      ! open the file (if not already open)
      if (checkpoint_file%status /= 1) then
          checkpoint_file%id = mh5_open_file_rw(checkpoint_file%name)
          if (checkpoint_file%id < 0) then
             ! invalid identifier, error when opening this file
             print*,"  - Could not open checkpoint file"
             checkpoint_valid = .NOT. VALID
             return
          end if
          checkpoint_file%status = 1
      end if

      ! check the existence of all required data
      do i = 1, size(labels)
         if (index(labels(i)%use,'required').ne.0) then
             if (.not.mh5_exists_dset(checkpoint_file%id,labels(i)%label)) then
                print*, "   - Required data set ",trim(labels(i)%label)," is missing"
                ! dataset does not exist
                checkpoint_valid = .NOT. VALID
             end if
         end if
      end do

  end function checkpoint_valid

  subroutine create_new_checkpoint

      ! Part of initialization procedure, hdf5 requires that a group is created
      ! before attempting to write data that is part of this group. Do this for
      ! all groups defined in the schema.

      use mh5
      integer :: i
      logical :: composite_data
      integer(kind=8) :: group_id

      ! create the hdf5 file
      checkpoint_file%id = mh5_create_file(checkpoint_file%name)
      checkpoint_file%status = 1

      ! Composite data corresponds to groups, find these labels and initialize the file by creating (empty) groups
      do i = 1, size(labels)
         if (index(labels(i)%type,'composite').ne.0) then
             group_id = mh5_create_group(checkpoint_file%id,labels(i)%label)
             call  mh5_close_group(group_id)
         end if
      end do

  end subroutine create_new_checkpoint

  subroutine get_labels

      ! Read the allowed labels from the schema file and store them.

      integer :: unit, lengths(5), i, n_labels=0, iostatus=0
      character(len=300):: line

      open(newunit=unit,form='formatted',file=SCHEMA_NAME)
      read(unit,*) lengths
      if (sum(lengths)>300) stop " fix module checkpoint to read long lines"

      ! Do a first pass to determine the number of labels
      do
         read(unit,'(A)',iostat=iostatus) line(1:sum(lengths))
         if (iostatus < 0) exit
         n_labels=n_labels+1
      end do

      ! Allocate the labels and read them in
      if (allocated(labels)) deallocate(labels) ! If the routine was called before, we simply delete the old set of labels
      allocate(labels(n_labels))
      iostatus = 0
      rewind(unit)
      read(unit,*) lengths
      do i = 1, n_labels
         read(unit,'(A)') line
         labels(i)%label        = line(1                  :lengths(1))
         labels(i)%type         = line(lengths(1)+1       :sum(lengths(1:2)))
         labels(i)%rank         = line(sum(lengths(1:2))+1:sum(lengths(1:3)))
         labels(i)%use          = line(sum(lengths(1:3))+1:sum(lengths(1:4)))
         labels(i)%description  = line(sum(lengths(1:4))+1:sum(lengths(1:5)))
      end do
      close(unit)

  end subroutine get_labels

  logical function check_label(mylabel,mytype)

      ! Check whether the label is known and refers to the right type of data

      character*(*), intent(in) :: mylabel, mytype
      logical :: found = .false.
      integer :: i

      ! Try to find the label in the list of allowed labels
      do i = 1, size(labels)
         ! if this label is generic, only the groups need to match
         if (trim(labels(i)%use) == 'generic') then
            found = index(group_from_label(mylabel),group_from_label(labels(i)%label)).ne.0
            if (found) exit
         end if
         ! normal case: a complete match is required
         found = index(labels(i)%label,mylabel).ne.0
         if (found) exit
      end do

      ! Check whether the type matches
      if (found) then
         if (index(labels(i)%type,mytype).ne.0) then
            check_label = VALID
         else
            check_label = .not.VALID
         end if
      else
         check_label = .not.VALID
      end if

  end function check_label

  function group_from_label(mylabel) result(group)

       character(len=*)            :: mylabel
       character, parameter        :: delimiter='/'
       character(:),allocatable    :: group
       integer                     :: i

       ! find rightmost delimiter
       i = scan(mylabel, delimiter, .true.)

       ! the first part is the group
       if (i > 1) then
          group = mylabel(1:i-1)
       else
          group = 'invalid_label'
       end if

  end function group_from_label

end module
