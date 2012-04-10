module unit_test_generator

implicit none

public :: add_inp_int_array, add_inp_real_array, add_out_int_array, add_out_real_array
public :: set_driver_routine, add_input_file, add_output_file, add_module, set_keys
public :: redefine_variable, set_suffix, get_output_name
public :: print_input_files, print_output_files, print_variable_info
public :: generate_test, unit_test_generator_destructor

public :: add_inp_variable, add_out_variable, add_opt_variable

interface add_inp_variable
  module procedure add_inp_logical_variable
  module procedure add_inp_integer_variable
  module procedure add_inp_integer_1d_array
  module procedure add_inp_integer_2d_array
  module procedure add_inp_real_variable
  module procedure add_inp_real_1d_array
  module procedure add_inp_real_2d_array
  module procedure add_inp_real_3d_array
  module procedure add_inp_real_4d_array
  module procedure add_inp_real_5d_array
end interface add_inp_variable

interface add_out_variable
  module procedure add_out_integer_variable
  module procedure add_out_integer_1d_array
  module procedure add_out_integer_2d_array
  module procedure add_out_real_variable
  module procedure add_out_real_1d_array
  module procedure add_out_real_2d_array
  module procedure add_out_real_3d_array
  module procedure add_out_real_4d_array
  module procedure add_out_real_5d_array
end interface add_out_variable

interface set_suffix
  module procedure set_suffix_string
  module procedure set_suffix_integer
  module procedure set_suffix_logical
end interface set_suffix

private :: init_string, add_string, errall, erropen, i2c, l2c, gets, csize, change_string
private :: add_data_routine, find_name, strip_blank, str_size, add_compare_info
private :: write_data, write_integer_array, write_real_array

interface write_data
  module procedure write_logical_variable
  module procedure write_integer_variable
  module procedure write_real_variable
end interface write_data

private

type :: string
  integer :: key
  character, pointer :: str(:)
end type string

type ( string ), pointer, save :: input_files(:), output_files(:), driver(:), test_name(:)
type ( string ), pointer, save :: var_name(:), var_type(:), var_dim(:)
type ( string ), pointer, save :: var_info(:), var_value(:)
type ( string ), pointer, save :: modules(:), suffix(:)
type ( string ), pointer, save :: compare_technique(:), thresholds(:), compare_array(:)
type ( string ), pointer, save :: compare_info(:)

type :: keys
  logical :: mod_call = .true.
end type keys

type ( keys ), save :: key

integer, parameter :: iunit = 60

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_int_array(ivar,vname,ndim)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim, ivar(*)

  ! Internal variables
  integer :: ier
  integer, allocatable :: itmp(:)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    allocate( itmp(ndim), stat=ier ); call errall(ier,'add_inp_int_array','itmp')
    itmp(1:ndim) = ivar(1:ndim)
  
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_int_array',gets(input_files(1)%str))
  
    write(iunit,'(a)')    vname
    write(iunit,'(6i12)') itmp
  
    close(iunit,status='keep')

    deallocate(itmp)
  else
    call write_integer_array(ivar,ndim,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer, allocatable, dimension(:)',var_type)
  call add_string('(' // i2c(ndim) // ')',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_int_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_array(dvar,vname,ndim)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim
  real*8,       intent(in) :: dvar(*)

  ! Internal variables
  integer :: ier
  real*8, allocatable :: dtmp(:)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    allocate( dtmp(ndim), stat=ier ); call errall(ier,'add_inp_real_array','dtmp')
    dtmp(1:ndim) = dvar(1:ndim)

    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_array',gets(input_files(1)%str))

    write(iunit,'(a)')       vname
    write(iunit,'(3D24.16)') dtmp

    close(iunit,status='keep')

    deallocate(dtmp)
  else
    call write_real_array(dvar,ndim,vname)
  endif
  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:)',var_type)
  call add_string('(' // i2c(ndim) // ')',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_int_array(vname,ndim)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer,  allocatable, dimension(:)',var_type)
  call add_string('(' // i2c(ndim) // ')',var_dim)
  call add_string('output',var_info)
  call add_string('own_name',var_value)

  end subroutine add_out_int_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_array(vname,ndim)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:)',var_type)
  call add_string('(' // i2c(ndim) // ')',var_dim)
  call add_string('output',var_info)
  call add_string('own_name',var_value)

  end subroutine add_out_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_logical_variable(lvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  logical,      intent(in) :: lvar

  ! Internal variables
  integer :: ier

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_logical_variable',gets(input_files(1)%str))

    write(iunit,'(a)') vname
    write(iunit,'(l)') lvar

    close(iunit,status='keep')
  else
    call write_data(lvar,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('logical',var_type)
  call add_string('0',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_logical_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_variable(ivar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar

  ! Internal variables
  integer :: ier

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_integer_variable',gets(input_files(1)%str))

    write(iunit,'(a)')   vname
    write(iunit,'(i12)') ivar

    close(iunit,status='keep')
  else
    call write_data(ivar,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer',var_type)
  call add_string('0',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_1d_array(ivar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar(:)

  ! Internal variables
  integer :: n1, ier

  n1 = size(ivar,1)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_integer_1d_array',gets(input_files(1)%str))

    write(iunit,'(a)')    vname
    write(iunit,'(6i12)') ivar

    close(iunit,status='keep')
  else
    call write_integer_array(ivar,n1,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer, allocatable, dimension(:)',var_type)
  call add_string( '(' // i2c(n1) // ')' , var_dim )
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_integer_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_2d_array(ivar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar(:,:)

  ! Internal variables
  integer :: n1, n2, ier

  n1 = size(ivar,1)
  n2 = size(ivar,2)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_integer_2d_array',gets(input_files(1)%str))
  
    write(iunit,'(a)')    vname
    write(iunit,'(6i12)') ivar
  
    close(iunit,status='keep')
  else
    call write_integer_array(ivar,n1*n2,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer, allocatable, dimension(:,:)',var_type)
  call add_string( '(' // i2c(n1) // ',' // i2c(n2) // ')' , var_dim )
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_integer_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_variable(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar

  ! Internal variables
  integer :: ier

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_variable',gets(input_files(1)%str))
  
    write(iunit,'(a)')      vname
    write(iunit,'(D24.16)') dvar
  
    close(iunit,status='keep')
  else
    call write_data(dvar,vname)
  endif


  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8',var_type)
  call add_string('0',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_1d_array(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:)

  ! Internal variables
  integer :: n1, ier

  n1 = size(dvar,1)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_1d_array',gets(input_files(1)%str))

    write(iunit,'(a)')       vname
    write(iunit,'(3D24.16)') dvar

    close(iunit,status='keep')
  else
    call write_real_array(dvar,n1,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:)',var_type)
  call add_string( '(' // i2c(n1) // ')' , var_dim )
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_2d_array(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:)

  ! Internal variables
  integer :: n1, n2, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_2d_array',gets(input_files(1)%str))

    write(iunit,'(a)')       vname
    write(iunit,'(3D24.16)') dvar

    close(iunit,status='keep')
  else
    call write_real_array(dvar,n1*n2,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:)',var_type)
  call add_string('('//i2c(n1)//','//i2c(n2)//')',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_3d_array(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_3d_array',gets(input_files(1)%str))

    write(iunit,'(a)')       vname
    write(iunit,'(3D24.16)') dvar

    close(iunit,status='keep')
  else
    call write_real_array(dvar,n1*n2*n3,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:,:)',var_type)
  call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//')',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_4d_array(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_4d_array',gets(input_files(1)%str))
  
    write(iunit,'(a)')       vname
    write(iunit,'(3D24.16)') dvar
  
    close(iunit,status='keep')
  else
    call write_real_array(dvar,n1*n2*n3*n4,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:,:,:)',var_type)
  call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//','//i2c(n4)//')',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_5d_array(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4, n5, ier

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)
  n5 = size(dvar,5)

  ! store input data to the file
  !   1.) input file which will be read at the beginning of the test routine
  !   2.) temporary file from which will be created subroutine which will return the input data
  if(associated(input_files))then
    open(iunit,file=gets(input_files(1)%str),status='unknown',form='formatted', &
         access='sequential',position='append',iostat=ier)
    call erropen(ier,'add_inp_real_5d_array',gets(input_files(1)%str))

    write(iunit,'(a)')       vname
    write(iunit,'(3D24.16)') dvar

    close(iunit,status='keep')
  else
    call write_real_array(dvar,n1*n2*n3*n4*n5,vname)
  endif

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:,:,:,:)',var_type)
  call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//','//i2c(n4)//','//i2c(n5)//')',var_dim)
  call add_string('input',var_info)
  call add_string('own_name',var_value)

  end subroutine add_inp_real_5d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_logical_variable(lvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  logical,      intent(in) :: lvar

  ! Internal variables
  integer :: ier

  open(iunit,file=vname,status='new',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'write_logical_variable',vname)

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', vname, '(', vname, ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  logical :: ',vname
  write(iunit,'(a)')
  if(lvar)then
    write(iunit,'(a,a,a)')  '  ', vname, ' = .True.'
  else
    write(iunit,'(a,a,a)')  '  ', vname, ' = .False.'
  endif
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine write_logical_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_integer_variable(ivar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar

  ! Internal variables
  integer :: ier

  open(iunit,file=vname,status='new',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'write_integer_variable',vname)
 
  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', vname, '(', vname, ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  integer :: ',vname
  write(iunit,'(a)')
  write(iunit,'(2x,a,a,i12)') vname, ' = ', ivar
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'
 
  close(iunit,status='keep')

  end subroutine write_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_real_variable(dvar,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar

  ! Internal variables
  integer :: ier

  open(iunit,file=vname,status='new',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'write_real_variable',vname)
 
  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', vname, '(', vname, ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  real*8 :: ',vname
  write(iunit,'(a)')
  write(iunit,'(2x,a,a,D24.16)') vname, ' = ', dvar
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'
 
  close(iunit,status='keep')

  end subroutine write_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_integer_array(ivar,ndim,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim, ivar(*)

  ! Internal variables
  integer :: i, ier

  open(iunit,file=vname,status='new',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'write_integer_array',vname)

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', vname, '(', vname, ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  integer, dimension(*) :: ',vname
  write(iunit,'(a)')
  do i=1,ndim
    write(iunit,'(2x,a,a,i12,a,6i12)') vname, '(', i, ') = ', ivar(i)
  enddo
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine write_integer_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine write_real_array(dvar,ndim,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim
  real*8,       intent(in) :: dvar(*)

  ! Internal variables
  integer :: i, ier

  open(iunit,file=vname,status='new',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'write_real_array',vname)

  write(iunit,'(a,a,a,a,a)')'  subroutine get_data_', vname, '(', vname, ')'
  write(iunit,'(a)')        '  implicit none'
  write(iunit,'(a,a)')      '  real*8, dimension(*) :: ',vname
  write(iunit,'(a)')
  do i=1,ndim
    write(iunit,'(2x,a,a,i12,a,3D24.16)') vname, '(', i, ') = ', dvar(i)
  enddo
  write(iunit,'(a)')
  write(iunit,'(a)')        '  end subroutine'

  close(iunit,status='keep')

  end subroutine write_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_data_routine(iw,vname)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: iw

  ! Internal variables
  integer :: i, io, ier
  character*(250) :: line

  io = iw + 1
  open(io,file=vname,status='old',form='formatted',access='sequential',iostat=ier)
  call erropen(ier,'add_data_routine',vname)
  rewind io

  do

  read(io,'(a)',iostat=ier) line
  if(ier == -1) exit
  write(iw,'(a)') strip_blank(line)

  enddo

  close(io,status='delete')

  end subroutine add_data_routine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  integer function find_name(string) result(ipos)
  implicit none

  ! External variables
  character(*), intent(in) :: string

  ! Internal variables
  integer :: i

  ipos = 0
  do i=1,size(var_name)
    if(string == gets(var_name(i)%str))then
      ipos = i
      exit
    endif
  enddo

  end function find_name

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

  subroutine add_compare_info(vname,ndim,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  integer,      intent(in) :: ndim

  call add_string(compare,compare_info)
  select case(compare)
    case('array-array')
      call add_string('  call compare_2vectors(' // i2c(ndim), compare_technique)
      call add_string(strip_blank(vname), compare_array)

    case('file-array')
      call add_string('  call compare_vector', compare_technique)
      call add_string(strip_blank(vname), compare_array)
      if(.not.present(file))then
        write(6,*) 'ERROR: You forgot to define "file" name add_out_variable call'
        stop
      endif
      ! add file (with suffix) to the list of output files (output_files)
      if(associated(suffix))then
        call add_string(file // '_' // gets(suffix(1)%str), output_files)
      else
        call add_string(file, output_files)
      endif

    case('file-file')
      call add_string('  call compare_files', compare_technique)
      call add_string('0', compare_array)
      if(.not.present(file))then
        write(6,*) 'ERROR: You forgot to define "file" name in add_out_variable call'
        stop
      endif
      ! add file (with suffix) to the list of output files (output_files)
      if(associated(suffix))then
        call add_string(file // '_' // gets(suffix(1)%str), output_files)
      else
        call add_string(file, output_files)
      endif

    case default
      write(6,*) 'ERROR: Wrong optional parameter "compare" in add_out_variable call'
      write(6,*) '       Use string "file-array" or "file-file"'
      stop
  end select

  ! set threshold
  if(present(threshold))then
    call add_string(threshold, thresholds)
  else
    call add_string('1.0d-15', thresholds)
  endif

  end subroutine add_compare_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_integer_variable(ivar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('integer',var_type)
    call add_string('0',var_dim)
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,0,compare,threshold,file)
  else
    call add_compare_info(vname,0,'array-array',threshold,file)
    call write_data(ivar,strip_blank(vname) // '_ref')
    write(6,*) ' STOP: output integer data are not yet implemented'
    stop
  endif

  end subroutine add_out_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_integer_1d_array(ivar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar(:)

  ! Internal variables
  integer :: n1

  n1 = size(ivar,1)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('integer, allocatable, dimension(:)',var_type)
    call add_string( '(' // i2c(n1) // ')' , var_dim )
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1,compare,threshold,file)
  else
    call add_compare_info(vname,n1,'array-array',threshold,file)
    call write_integer_array(ivar,n1,strip_blank(vname) // '_ref')
    write(6,*) ' STOP: output integer data are not yet implemented'
    stop
  endif

  end subroutine add_out_integer_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_integer_2d_array(ivar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar(:,:)

  ! Internal variables
  integer :: n1, n2

  n1 = size(ivar,1)
  n2 = size(ivar,2)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('integer, allocatable, dimension(:,:)',var_type)
    call add_string( '(' // i2c(n1) // ',' // i2c(n2) // ')' , var_dim )
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1*n2,compare,threshold,file)
  else
    call add_compare_info(vname,n1*n2,'array-array',threshold,file)
    call write_integer_array(ivar,n1*n2,strip_blank(vname) // '_ref')
    write(6,*) ' STOP: output integer data are not yet implemented'
    stop
  endif

  end subroutine add_out_integer_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_variable(dvar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('real*8',var_type)
    call add_string('0',var_dim)
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,0,compare,threshold,file)
  else
    call add_compare_info(vname,0,'array-array',threshold,file)
    call write_data(dvar,strip_blank(vname) // '_ref')
    write(6,*) ' STOP: output real variable is not yet implemented'
    stop
  endif

  end subroutine add_out_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_1d_array(dvar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:)

  ! Internal variables
  integer :: n1

  n1 = size(dvar,1)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('real*8,  allocatable, dimension(:)',var_type)
    call add_string( '(' // i2c(n1) // ')' , var_dim )
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1,compare,threshold,file)
  else
    call add_compare_info(vname,n1,'array-array',threshold,file)
    call write_real_array(dvar,n1,strip_blank(vname) // '_ref')
  endif

  end subroutine add_out_real_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_2d_array(dvar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:)

  ! Internal variables
  integer :: n1, n2

  n1 = size(dvar,1)
  n2 = size(dvar,2)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('real*8,  allocatable, dimension(:,:)',var_type)
    call add_string('('//i2c(n1)//','//i2c(n2)//')',var_dim)
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1*n2,compare,threshold,file)
  else
    call add_compare_info(vname,n1*n2,'array-array',threshold,file)
    call write_real_array(dvar,n1*n2,strip_blank(vname) // '_ref')
  endif

  end subroutine add_out_real_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_3d_array(dvar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:)

  ! Internal variables
  integer :: n1, n2, n3

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('real*8,  allocatable, dimension(:,:,:)',var_type)
    call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//')',var_dim)
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1*n2*n3,compare,threshold,file)
  else
    call add_compare_info(vname,n1*n2*n3,'array-array',threshold,file)
    call write_real_array(dvar,n1*n2*n3,strip_blank(vname) // '_ref')
  endif

  end subroutine add_out_real_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_4d_array(dvar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('real*8,  allocatable, dimension(:,:,:,:)',var_type)
    call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//','//i2c(n4)//')',var_dim)
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1*n2*n3*n4,compare,threshold,file)
  else
    call add_compare_info(vname,n1*n2*n3*n4,'array-array',threshold,file)
    call write_real_array(dvar,n1*n2*n3*n4,strip_blank(vname) // '_ref')
  endif

  end subroutine add_out_real_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_5d_array(dvar,vname,compare,threshold,file)
  implicit none

  ! External variables
  character(*), intent(in), optional :: compare, threshold, file
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:,:,:)

  ! Internal variables
  integer :: n1, n2, n3, n4, n5

  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)
  n5 = size(dvar,5)

  ! if vname is already defined just change its status to "inout" variable and write the data
  if(find_name(vname) == 0)then
    ! add name, type, dimension, info and order in driver call
    !   of the variable to the corresponding list
    call add_string(vname,var_name)
    call add_string('real*8,  allocatable, dimension(:,:,:,:,:)',var_type)
    call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//','//i2c(n4)//','//i2c(n5)//')',var_dim)
    call add_string('output',var_info)
    call add_string('own_name',var_value)
  else
    call change_string(var_info,'inout',find_name(vname))
  endif

  ! set what you want to compare
  if(present(compare))then
    call add_compare_info(vname,n1*n2*n3*n4*n5,compare,threshold,file)
  else
    call add_compare_info(vname,n1*n2*n3*n4*n5,'array-array',threshold,file)
    call write_real_array(dvar,n1*n2*n3*n4*n5,strip_blank(vname) // '_ref')
  endif

  end subroutine add_out_real_5d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_opt_variable(vname,vvalue,vtype)
  implicit none
  character(*), intent(in) :: vname, vvalue, vtype

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string(vtype,var_type)
  call add_string('0',var_dim)
  call add_string('optional',var_info)
  call add_string(vvalue,var_value)

  end subroutine add_opt_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_driver_routine(drv_name)
  implicit none
  character(*), intent(in) :: drv_name

  ! add name of driver routine to the list
  call add_string(drv_name,driver)

  end subroutine set_driver_routine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_input_file(file_name)
  implicit none
  character(*), intent(in) :: file_name

  ! add file_name (with suffix) to the list of input files (input_files)
  if(associated(suffix))then
    call add_string(file_name // '_' // gets(suffix(1)%str), input_files)
  else
    call add_string(file_name, input_files)
  endif

  end subroutine add_input_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_output_file(file_name)
  implicit none
  character(*), intent(in) :: file_name

  ! add file_name (with suffix) to the list of output files (output_files)
  if(associated(suffix))then
    call add_string(file_name // '_' // gets(suffix(1)%str), output_files)
  else
    call add_string(file_name, output_files)
  endif

  end subroutine add_output_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_module(string)
  implicit none
  character(*), intent(in) :: string

  ! add module to the test
  call add_string(string,modules)

  end subroutine add_module

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_keys(mod_call)
  implicit none
  logical, intent(in), optional :: mod_call

  ! set keys
  key%mod_call = mod_call

  end subroutine set_keys

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine redefine_variable(variable_type,variable_size,nvar)
  implicit none

  ! External variables
  character(*), intent(in), optional :: variable_type, variable_size
  integer,      intent(in), optional :: nvar

  ! Internal variables
  integer :: i

  ! position of the string
  i = size(var_name)
  if(present(nvar)) i = nvar

  ! reset
  if(present(variable_type)) call change_string(var_type,variable_type,i)
  if(present(variable_size)) call change_string(var_dim, variable_size,i)

  end subroutine redefine_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_suffix_string(string)
  implicit none
  character(*), intent(in) :: string

  ! add string to the name of the test source file
  if(associated(suffix))then
    call change_string(suffix, gets(suffix(1)%str) // '_' // string, 1)
  else
    call add_string(string, suffix)
  endif

  end subroutine set_suffix_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_suffix_integer(ivar)
  implicit none
  integer, intent(in) :: ivar

  ! add integer to the name of the test source file
  if(associated(suffix))then
    call change_string(suffix, gets(suffix(1)%str) // '_' // i2c(ivar), 1)
  else
    call add_string(i2c(ivar), suffix)
  endif

  end subroutine set_suffix_integer

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_suffix_logical(lvar)
  implicit none
  logical, intent(in) :: lvar

  ! add integer to the name of the test source file
  if(associated(suffix))then
    call change_string(suffix, gets(suffix(1)%str) // '_' // l2c(lvar), 1)
  else
    call add_string(l2c(lvar), suffix)
  endif

  end subroutine set_suffix_logical

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function get_output_name(nout) result(out_name)
  implicit none
  integer, intent(in) :: nout
  character(len=size(output_files(nout)%str)) :: out_name
  integer :: i

  do i=1,size(output_files(nout)%str)
    out_name(i:i)=output_files(nout)%str(i)
  enddo

  end function get_output_name

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine print_input_files()
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

  subroutine print_output_files()
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

  subroutine print_variable_info()
  implicit none
  integer :: i

  write(6,'(/a)') ' VARIABLE INFO:'
  do i=1,size(var_name)
    write(6,'(a,i2)') ' Variable #',i
    write(6,*) '  name     : ',var_name(i)%str
    write(6,*) '  dimension: ',var_dim(i)%str
    write(6,*) '  type     : ',var_type(i)%str
    write(6,*) '  info     : ',var_info(i)%str
    write(6,*) '  key      : ',var_name(i)%key
  enddo

  end subroutine print_variable_info

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
      type_string(i)%key =  tmp_list(i)%key
    enddo
  endif

  ! allocate unnamed array str character array
  nstr = len(string_name)
  if(nstr > 70)then
    write(6,*)
    write(6,*) 'ERROR: Input string is too long. It would cause problems in test.F90 generation.'
    write(6,*) string_name
    stop
  endif
  allocate( type_string(ndim+1)%str( nstr ), stat = ier )
  call errall(ier,'add_string','type_string()%str')

  ! copy string_name to newly allocated str array
  do i=1,nstr
    type_string(ndim+1)%str(i) = string_name(i:i)
    type_string(ndim+1)%key    = 0
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
    stop 'ERROR in allocation (see output for more information)'
  endif

  end subroutine errall

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine erropen(ier,chsub,chfile)
  implicit none

  integer      :: ier
  character(*) :: chsub, chfile

  if(ier /= 0)then
    write(6,*)
    write(6,'(a,a,a,a)')' ERROR while opening ',chfile,' in subroutine ',chsub
    write(6,'(a,i4)'   )' iostat =',ier
    stop 'ERROR in opening file (see output for more information)'
  endif

  end subroutine erropen

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function i2c(i) result(c)
  implicit none
  integer, intent(in) :: i
  character(len=csize(i)) :: c
  write(c,'(i0)')i
  end function i2c

  function l2c(l) result(c)
  implicit none
  logical, intent(in) :: l
  character(len=1) :: c
  if(l) c='T'
  if(.not.l) c='F'
  end function l2c

  function gets(str) result(c)
  implicit none
  character, pointer, intent(in) :: str(:)
  character(len=size(str)) :: c
  integer :: i
  do i=1,size(str)
    c(i:i)=str(i)
  enddo
  end function gets

  pure function csize(i) result(is)
  implicit none
  integer, intent(in) :: i
  integer :: is
  is = floor(log10(real(i,kind(1.0d0)))) + 1
  end function csize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine change_string(type_string,new_string,npos)
  implicit none

  ! External variables
  character(*), intent(in) :: new_string
  integer,      intent(in) :: npos
  type ( string ), pointer, intent(inout) :: type_string(:)

  ! Internal variables
  integer :: i, ier, nstr

  ! get length of new string
  nstr = len(new_string)

  ! nullify old string
  nullify( type_string(npos)%str )

  ! allocate string array with new length
  allocate( type_string(npos)%str( nstr ), stat = ier )
  call errall(ier,'change_string','type_string(npos)%str')

  ! reset
  do i=1,nstr
    type_string(npos)%str(i) = new_string(i:i)
  enddo

  end subroutine change_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine generate_test
  implicit none

  ! Internal variables
  character :: coma*1
  integer   :: i, io, ier

  ! set name of the test source file
  if(associated(suffix))then
    call add_string('test.F90_' // gets(suffix(1)%str), test_name)
  else
    call add_string('test.F90', test_name)
  endif

  open(iunit,file=gets(test_name(1)%str),status='new',form='formatted',access='sequential',&
       iostat=ier); call erropen(ier,'generate_test',gets(test_name(1)%str))

  ! header
  write(iunit,'(a)')'program test'
  write(iunit,'(a)')
  write(iunit,'(a)')'  use unit_testing'
  do i=1,size(modules)
    write(iunit,'(a,100a)')'  use ' , gets(modules(i)%str)
  enddo
  write(iunit,'(a)')'  implicit none'
  write(iunit,'(a)')
  write(iunit,'(a)')'  integer :: ier'
  ! declaration of variables
  io = 0
  do i=1,size(var_name)
    write(iunit,'(2x,70a,a,70a)') gets(var_type(i)%str),' :: ',gets(var_name(i)%str)
    ! declare reference array if we want to use array-array comparison
    if( gets(var_info(i)%str) == 'output' .or. gets(var_info(i)%str) == 'inout' )then
      io = io + 1
      if(gets(compare_info(io)%str) == 'array-array')then
        write(iunit,'(2x,70a,a,70a)') gets(var_type(i)%str),' :: ', &
                                      strip_blank(gets(var_name(i)%str)) // '_ref'
      endif
    endif
  enddo
  ! allocate arrays
  io = 0
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( var_dim(i)%str(1) /= '0')then
      write(iunit,'(a,70a,70a,a,70a,a)') &
            '  allocate( ' , gets(var_name(i)%str) , gets(var_dim(i)%str) , &
            " , stat=ier ); call errall(ier,'test','" , gets(var_name(i)%str) ,"')"
      ! allocate reference array if we want to use array-array comparison
      if( gets(var_info(i)%str) == 'output' .or. gets(var_info(i)%str) == 'inout' )then
        io = io + 1
        if(gets(compare_info(io)%str) == 'array-array')then
          write(iunit,'(a,70a,70a,a,70a,a)') '  allocate( ' , &
                strip_blank(gets(var_name(i)%str)) // '_ref' , gets(var_dim(i)%str) , &
                " , stat=ier ); call errall(ier,'test','" , gets(var_name(i)%str) ,"')"
        endif
      endif
    endif
  enddo
  ! open input file (if data are in the file)
  if(associated(input_files))then
    if(size(input_files) > 1)then
      write(6,*) 'ERROR: For now the program works only with one input file.'
      stop
    endif
    write(iunit,'(a)')
    write(iunit,'(a,70a,a)') "  open(60,file='" , gets(input_files(1)%str)  , &
                               "',status='old',form='formatted',access='sequential',iostat=ier)"
    write(iunit,'(a)')       '  if(ier /= 0)then'
    write(iunit,'(a,70a,a)') "    write(6,*) 'test: ERROR while creating " , &
                               gets(input_files(1)%str) , " file'"
    write(iunit,'(a)')       '    stop'
    write(iunit,'(a)')       '  endif'
  endif
  ! get input data
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( gets(var_info(i)%str) == 'input' .or. gets(var_info(i)%str) == 'inout' )then
      if(associated(input_files))then
        write(iunit,'(a)')     '  read(60,*)'
        write(iunit,'(a,70a)') '  read(60,*) ' , gets(var_name(i)%str)
      else
        write(iunit,'(a,70a,a,70a,a)') &
              '  call get_data_', gets(var_name(i)%str), '(', gets(var_name(i)%str), ')'
      endif
    endif
  enddo
  ! close input file
  if(associated(input_files))then
    write(iunit,'(a)')
    write(iunit,'(a)')"  close(60,status='keep')"
  endif
  ! set optional parameters
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( gets(var_info(i)%str) == 'optional' )then
      write(iunit,'(2x,70a)') gets(var_name(i)%str), ' = ', gets(var_value(i)%str)
    endif
  enddo
  ! call driver routine
  write(iunit,'(a,70a,a)')"  call ", gets(driver(1)%str) , '( &'
  coma = ' '
  if(key%mod_call)then
    do i=1,size(var_name)
      write(iunit,'(9x,a,70a,a)')coma,gets(var_name(i)%str),' = ',gets(var_name(i)%str),' &'
      coma = ','
    enddo
  else
    do i=1,size(var_name)
      write(iunit,'(9x,a,70a,a)')coma,gets(var_name(i)%str),' &'
      coma = ','
    enddo
  endif
  write(iunit,'(9x,a)') ')'
  ! get output data
  io = 0
  do i=1,size(var_name)
    if( gets(var_info(i)%str) == 'output' .or. gets(var_info(i)%str) == 'inout' )then
      io = io + 1
      if(gets(compare_info(io)%str) == 'array-array')then
        write(iunit,'(a)')
        write(iunit,'(a,70a,a,70a,a)') '  call get_data_', &
              strip_blank(gets(var_name(i)%str)), '_ref(', & 
              strip_blank(gets(var_name(i)%str)), '_ref)'
      endif
    endif
  enddo
  ! compare results
  write(iunit,'(a)')
  do i=1,size(compare_technique)
    if(gets(compare_info(i)%str) == 'file-file')then
      write(iunit,'(70a,a,70a,a,70a,a,70a,a)') &
           gets(compare_technique(i)%str), "(60,'", gets(output_files(i)%str), "','", &
           gets(output_files(i)%str), "_ref',", gets(thresholds(i)%str), ')'
    elseif(gets(compare_info(i)%str) == 'file-array')then
      write(iunit,'(70a,a,70a,a,70a,a,70a,a)') &
           gets(compare_technique(i)%str), "(60,'", gets(output_files(i)%str), "',", &
           gets(compare_array(i)%str), ",", gets(thresholds(i)%str), ')'
    elseif(gets(compare_info(i)%str) == 'array-array')then
      write(iunit,'(70a,a,70a,a,70a,a,70a,a)') &
           gets(compare_technique(i)%str), ',', gets(compare_array(i)%str), '_ref,', &
           gets(compare_array(i)%str), ',', gets(thresholds(i)%str), ')'
    endif
  enddo
  ! deallocation
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( var_dim(i)%str(1) /= '0')then
      write(iunit,'(a,70a,a)')"  deallocate( " , gets(var_name(i)%str) , " )"
    endif
  enddo
  ! write error routine
  write(iunit,'(a)')
  write(iunit,'(a)')"contains"
  write(iunit,'(a)')
  write(iunit,'(a)')"  subroutine errall(ier,chsub,charr)"
  write(iunit,'(a)')"  implicit none"
  write(iunit,'(a)')
  write(iunit,'(a)')"  integer      :: ier"
  write(iunit,'(a)')"  character(*) :: chsub, charr"
  write(iunit,'(a)')
  write(iunit,'(a)')"  if(ier /= 0)then"
  write(iunit,'(a)')"    write(6,*)"
  write(iunit,'(a)')"    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub"
  write(iunit,'(a)')"    write(6,'(a,i4)'   )' iostat =',ier"
  write(iunit,'(a)')"    stop 'ERROR in allocation (see output form more information)'"
  write(iunit,'(a)')"  endif"
  write(iunit,'(a)')
  write(iunit,'(a)')"  end subroutine"
  ! add input data routines
  if(.not. associated(input_files))then
    do i=1,size(var_name)
      if( gets(var_info(i)%str) == 'input' .or. gets(var_info(i)%str) == 'inout' )then
        write(iunit,'(a)')
        call add_data_routine(iunit,gets(var_name(i)%str))
      endif
    enddo
  endif
  ! add output data routines
  io = 0
  do i=1,size(var_name)
    if( gets(var_info(i)%str) == 'output' .or. gets(var_info(i)%str) == 'inout' )then
      io = io + 1
      if(gets(compare_info(io)%str) == 'array-array')then
        write(iunit,'(a)')
        call add_data_routine(iunit,strip_blank(gets(var_name(i)%str)) // '_ref')
      endif
    endif
  enddo
  ! end of the test program
  write(iunit,'(a)')"end program"

  close(iunit,status='keep')

  call unit_test_generator_destructor

  end subroutine generate_test

  subroutine unit_test_generator_destructor
  implicit none

  ! clean all information about just generated source
  ! so in one run can be generated multiple files
  nullify( input_files ); nullify( output_files ); nullify( driver ); nullify( test_name )
  nullify( var_name ); nullify( var_type ); nullify( var_dim )
  nullify( var_info ); nullify( var_value )
  nullify( modules ); nullify( suffix )
  nullify( compare_technique ); nullify( thresholds ); nullify( compare_array )
  nullify( compare_info )

  end subroutine unit_test_generator_destructor

end module unit_test_generator
