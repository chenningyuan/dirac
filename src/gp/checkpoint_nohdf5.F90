
! Stub for handling checkpoint calls in case no hdf5 library was found.
! No checkpoint file is written / read in this case, obviously.

! Written by Lucas Visscher, fall 2021

module checkpoint

  use labeled_storage

  implicit none

  private

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

      print*, "   - DIRAC was compiled without hdf5, checkpointing is not possible"

  end subroutine checkpoint_open

  subroutine checkpoint_close

      print*, "   - DIRAC was compiled without hdf5, checkpointing was not possible"

  end subroutine checkpoint_close

  subroutine checkpoint_write(label,rdata,idata,sdata)

      character*(*),intent(in)               :: label
      real(8),intent(in), optional           :: rdata(:)
      integer,intent(in), optional           :: idata(:)
      character*(*),intent(in), optional     :: sdata

      return

  end subroutine checkpoint_write

  subroutine checkpoint_read(label,rdata,idata,sdata)

      character*(*),intent(in)               :: label
      real(8),intent(out), optional          :: rdata(:)
      integer,intent(out), optional          :: idata(:)
      character*(*), optional                :: sdata

      return

  end subroutine checkpoint_read

  subroutine write_single_real (label,rdata)
      ! Wrapper to allow writing of single numbers
      character*(*),intent(in)               :: label
      real(8),intent(in)                     :: rdata

      return

  end subroutine write_single_real

  subroutine read_single_real (label,rdata)
      ! Wrapper to allow reading of single numbers
      character*(*),intent(in)               :: label
      real(8)                                :: rdata

      return

  end subroutine read_single_real

  subroutine write_single_int (label,idata)
      ! Wrapper to allow writing of single numbers
      character*(*),intent(in)               :: label
      integer,intent(in)                     :: idata

      return

  end subroutine write_single_int

  subroutine read_single_int (label,idata)
      ! Wrapper to allow reading of single numbers
      character*(*),intent(in)               :: label
      integer                                :: idata

      return

  end subroutine read_single_int

  subroutine write_string_array (label,sarray,slen)
      ! Wrapper to write string arrays
      character*(*),intent(in)               :: label
      integer                                :: slen
      character(len=slen),intent(in)         :: sarray(:)
      logical                                :: valid

      return

  end subroutine write_string_array

  subroutine read_string_array (label,sarray,slen)
      ! Wrapper to read string arrays
      character*(*),intent(in)               :: label
      integer                                :: slen
      character(len=slen)                    :: sarray(:)

      return

  end subroutine read_string_array

end module
