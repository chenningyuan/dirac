module test_generator

implicit none

public  :: add_input_file, add_output_file, print_input_files, print_output_files

private :: init_string, add_string
private :: errall

private

type :: string
  character, pointer :: str(:)
end type string

type ( string ), pointer, save :: input_files(:), output_files(:)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_input_file(file_name)
  implicit none

  ! External variables
  character(*), intent(in) :: file_name

  ! add file_name to the list of input files (input_files)
  call add_string(file_name,input_files)

  end subroutine add_input_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_output_file(file_name)
  implicit none

  ! External variables
  character(*), intent(in) :: file_name

  ! add file_name to the list of output files (output_files)
  call add_string(file_name,output_files)

  end subroutine add_output_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine print_input_files
  implicit none
  integer :: i

  write(6,'(/a)') ' LIST OF INPUT FILES:'
  do i=1,size(input_files)
    write(6,'(a,i2)') ' File #',i
    write(6,*) '  ',input_files(i)%str
  enddo

  end subroutine print_input_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine print_output_files
  implicit none
  integer :: i

  write(6,'(/a)') ' LIST OF OUTPUT FILES:'
  do i=1,size(output_files)
    write(6,'(a,i2)') ' File #',i
    write(6,*) '  ',output_files(i)%str
  enddo

  end subroutine print_output_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine init_string(type_string,ndim)
  !==============================================================================!
  !     allocate unnamed array "type_string" with dimension of the pointer       !
  !   "input_files" plus one and nullify the character pointer of string type    !
  !==============================================================================!
  implicit none

  ! External variables
  integer, intent(in) :: ndim
  type ( string ), pointer, intent(out) :: type_string(:)

  ! Internal variables
  integer :: i, ier
  
  allocate(type_string(ndim), stat = ier); call errall(ier,'init_string','type_string')

  do i=1,ndim
    nullify(type_string(i)%str)
  enddo

  end subroutine init_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_string(string_name,type_string)
  !==============================================================================!
  !      Add string "string_name" at the end of the "type_string" array          !
  !==============================================================================!
  implicit none

  ! External variables
  character(*), intent(in) :: string_name
  type ( string ), pointer, intent(inout) :: type_string(:)

  ! Internal variables
  integer :: i, ier, ndim, nstr
  type ( string ), pointer :: tmp_list(:)

  ! initialize tmp_list
  ndim = 0
  if(associated(type_string))then
    ndim = size(type_string)
    call init_string(tmp_list,ndim)
  endif

  ! save type_string in tmp_list
  if(associated(type_string)) tmp_list => type_string

  ! reallocate type_string with greater dimension
  nullify(type_string)
  allocate( type_string(ndim+1), stat = ier ); call errall(ier,'add_string','type_string')

  ! copy tmp_list back to type_string
  if(ndim > 0)then
    do i=1,ndim
      type_string(i)%str => tmp_list(i)%str
    enddo
  endif

  ! allocate unnamed array str character array
  nstr = len(string_name)
  allocate( type_string(ndim+1)%str( nstr ), stat = ier )
  call errall(ier,'add_string','type_string()%str')

  ! copy string_name to newly allocated str array
  do i=1,nstr
    type_string(ndim+1)%str(i) = string_name(i:i)
  enddo

  ! nullify tmp_list
  nullify( tmp_list )

  end subroutine add_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine errall(ier,chsub,charr)
  implicit none

  integer      :: ier
  character(*) :: chsub, charr

  if(ier /= 0)then
    write(6,*)
    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
    write(6,'(a,i4)'   )' iostat =',ier
    stop 'ERROR in allocation (see output form more information)'
  endif

  end subroutine errall

end module test_generator
