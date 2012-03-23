module test_generator

implicit none

public :: add_inp_int_array, add_inp_real_array, add_out_int_array, add_out_real_array
public :: set_driver_routine, add_input_file, add_output_file
public :: print_input_files, print_output_files, print_variable_info
public :: generate_test

public :: add_inp_variable

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

private :: init_string, add_string, errall, erropen, get_name, i2c, gets, csize, check_idrv

private

type :: string
  integer :: idrv
  character, pointer :: str(:)
end type string

type ( string ), pointer, save :: input_files(:), output_files(:), driver(:)
type ( string ), pointer, save :: var_name(:), var_type(:), var_dim(:), var_info(:)

integer, parameter :: iunit = 60, nchar = 100

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_int_array(ivar,vname,ndim,vsize,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname, vsize
  integer,      intent(in) :: ndim, ivar(*)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, ier
  character :: namef*(nchar)
  integer, allocatable :: itmp(:)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write array to the file
  allocate( itmp(ndim), stat=ier ); call errall(ier,'add_inp_int_array','itmp')
  itmp(1:ndim) = ivar(1:ndim)

  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_int_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) itmp

  close(iunit,status='keep')

  deallocate(itmp)
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer, allocatable, dimension(:)',var_type)
  call add_string(vsize,var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_int_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_array(dvar,vname,ndim,vsize,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname, vsize
  integer,      intent(in) :: ndim
  real*8,       intent(in) :: dvar(*)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, ier
  character :: namef*(nchar)
  real*8, allocatable :: dtmp(:)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write array to the file
  allocate( dtmp(ndim), stat=ier ); call errall(ier,'add_inp_real_array','dtmp')
  dtmp(1:ndim) = dvar(1:ndim)

  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dtmp

  close(iunit,status='keep')

  deallocate(dtmp)
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:)',var_type)
  call add_string(vsize,var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_int_array(vname,vsize,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname, vsize
  integer, intent(in), optional :: idrv

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer,  allocatable, dimension(:)',var_type)
  call add_string(vsize,var_dim)
  call add_string('output',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_out_int_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_out_real_array(vname,vsize,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname, vsize
  integer, intent(in), optional :: idrv

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:)',var_type)
  call add_string(vsize,var_dim)
  call add_string('output',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_out_real_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_logical_variable(lvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  logical,      intent(in) :: lvar
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write variable to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_logical_variable',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) lvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('logical',var_type)
  call add_string('0',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_logical_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_variable(ivar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write variable to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_integer_variable',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) ivar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer',var_type)
  call add_string('0',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_integer_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_1d_array(ivar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar(:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_integer_1d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) ivar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer, allocatable, dimension(:)',var_type)
  n1 = size(ivar,1)
  call add_string( '(' // i2c(n1) // ')' , var_dim )
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_integer_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_integer_2d_array(ivar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  integer,      intent(in) :: ivar(:,:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, n2, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_integer_2d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) ivar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('integer, allocatable, dimension(:,:)',var_type)
  n1 = size(ivar,1)
  n2 = size(ivar,2)
  call add_string( '(' // i2c(n1) // ',' // i2c(n2) // ')' , var_dim )
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_integer_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_variable(dvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write variable to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_variable',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8',var_type)
  call add_string('0',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_variable

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_1d_array(dvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_1d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:)',var_type)
  n1 = size(dvar,1)
  call add_string( '(' // i2c(n1) // ')' , var_dim )
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_1d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_2d_array(dvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, n2, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_2d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:)',var_type)
  n1 = size(dvar,1)
  n2 = size(dvar,2)
  call add_string('('//i2c(n1)//','//i2c(n2)//')',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_2d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_3d_array(dvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, n2, n3, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_3d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:,:)',var_type)
  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//')',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_3d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_4d_array(dvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:,:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, n2, n3, n4, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_4d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:,:,:)',var_type)
  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)
  call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//','//i2c(n4)//')',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_4d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine add_inp_real_5d_array(dvar,vname,nfile,idrv)
  implicit none

  ! External variables
  character(*), intent(in) :: vname
  real*8,       intent(in) :: dvar(:,:,:,:,:)
  integer, intent(in), optional :: nfile, idrv

  ! Internal variables
  integer   :: n, n1, n2, n3, n4, n5, ier
  character :: namef*(nchar)

  ! get name of the input file
  call get_name(namef,n,input_files,nfile)

  ! --\-- write 2d array to the input file
  open(iunit,file=namef(1:n),status='unknown',form='formatted',access='sequential',&
       position='append',iostat=ier); call erropen(ier,'add_inp_real_5d_array',namef(1:n))

  write(iunit,*) vname
  write(iunit,*) dvar

  close(iunit,status='keep')
  ! --/--

  ! add name, type, dimension, info and order in driver call
  !   of the variable to the corresponding list
  call add_string(vname,var_name)
  call add_string('real*8,  allocatable, dimension(:,:,:,:,:)',var_type)
  n1 = size(dvar,1)
  n2 = size(dvar,2)
  n3 = size(dvar,3)
  n4 = size(dvar,4)
  n5 = size(dvar,5)
  call add_string('('//i2c(n1)//','//i2c(n2)//','//i2c(n3)//','//i2c(n4)//','//i2c(n5)//')',var_dim)
  call add_string('input',var_info)
  if(present(idrv)) var_name(size(var_name))%idrv = idrv

  end subroutine add_inp_real_5d_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine set_driver_routine(drv_name)
  implicit none

  ! External variables
  character(*), intent(in) :: drv_name

  ! add name of driver routine to the list
  call add_string(drv_name,driver)

  end subroutine set_driver_routine

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

  subroutine print_variable_info
  implicit none
  integer :: i

  write(6,'(/a)') ' VARIABLE INFO:'
  do i=1,size(var_name)
    write(6,'(a,i2)') ' Variable #',i
    write(6,*) '  name     : ',var_name(i)%str
    write(6,*) '  dimension: ',var_dim(i)%str
    write(6,*) '  type     : ',var_type(i)%str
    write(6,*) '  info     : ',var_info(i)%str
    write(6,*) '  idrv     : ',var_name(i)%idrv
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
      type_string(i)%idrv = tmp_list(i)%idrv
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
    type_string(ndim+1)%idrv   = 0
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

  subroutine get_name(namef,n,type_string,nfile)
  implicit none

  ! External variables
  character, intent(out) :: namef*(nchar)
  integer,   intent(out) :: n
  integer,   intent( in), optional :: nfile
  type ( string ), pointer, intent(in) :: type_string(:)
  
  ! Internal variables
  integer :: i, ifile

  ifile = 1
  if(present(nfile)) ifile = nfile

  n = size(type_string(ifile)%str)

  if(n > nchar)then
    write(6,*)
    write(6,*) 'ERROR in test_generator module (get_name routine)'
    write(6,*) '  the size of the name of some file is too big; enlarge nchar'
    write(6,*) 'Problematic file:'
    write(6,*) type_string(ifile)%str
    stop
  endif

  do i=1,n
    namef(i:i) = type_string(ifile)%str(i)
  enddo

  end subroutine get_name

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  function i2c(i) result(c)
  implicit none
  integer, intent(in) :: i
  character(len=csize(i)) :: c
  write(c,'(i0)')i
  end function i2c

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

  function check_idrv result(lg)
  implicit none
  logical :: lg
  integer :: i
  lg = .false.
  do i=1,size(var_name)
    if(var_name(i)%idrv /= 0) lg = .true.
  enddo
  end function check_idrv

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine generate_test
  implicit none

  ! Internal variables
  character :: coma*1
  integer   :: i, j, k, ier

  open(iunit,file='test.F90',status='new',form='formatted',access='sequential',&
       iostat=ier); call erropen(ier,'generate_test','test.F90')

  ! header
  write(iunit,'(a)')'program test'
  write(iunit,'(a)')
  write(iunit,'(a)')'  use unit_testing'
  write(iunit,'(a)')'  implicit none'
  write(iunit,'(a)')
  write(iunit,'(a)')'  integer :: ier'
  write(iunit,'(a)')'  integer :: key'
  ! declaration of variables
  do i=1,size(var_name)
    write(iunit,'(2x,70a,a,70a)') gets(var_type(i)%str),' :: ',gets(var_name(i)%str)
  enddo
  ! open input file
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
  ! first read dimensions and parameters
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( var_info(i)%str(1) == 'i')then ! input variables
      if( var_dim(i)%str(1) == '0')then
        write(iunit,'(a)')     '  read(60,*)'
        write(iunit,'(a,70a)') '  read(60,*) ' , gets(var_name(i)%str)
      endif
    endif
  enddo
  ! allocation
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( var_dim(i)%str(1) /= '0')then
      write(iunit,'(a,70a,70a,a,70a,a)') &
            '  allocate( ' , gets(var_name(i)%str) , gets(var_dim(i)%str) , &
            " , stat=ier ); call errall(ier,'test','" , gets(var_name(i)%str) ,"')"
    endif
  enddo
  ! read arrays
  write(iunit,'(a)')
  do i=1,size(var_name)
    if( var_info(i)%str(1) == 'i')then ! input variables
      if( var_dim(i)%str(1) /= '0')then
        write(iunit,'(a)')    '  read(60,*)'
        write(iunit,'(a,70a)')'  read(60,*) ' , gets(var_name(i)%str)
      endif
    endif
  enddo
  ! close input file
  write(iunit,'(a)')
  write(iunit,'(a)')"  close(60,status='keep')"
  ! call driver routine
  write(iunit,'(a)')
  write(iunit,'(a)') '  key = 1'
  write(iunit,'(a,70a,a)')"  call ", gets(driver(1)%str) , '( &'
  coma = ' '
  if(check_idrv())then
    do k=1,size(var_name)
      do i=1,size(var_name)
        if(var_name(i)%idrv == k)then
          write(iunit,'(9x,a,70a,a)') coma,gets(var_name(i)%str) , ' &'
          coma = ','
        endif
      enddo
    enddo
    write(iunit,'(9x,a)') ',key    &'
  else
    do i=1,size(var_name)
      write(iunit,'(9x,a,70a,a)')coma,gets(var_name(i)%str),' = ',gets(var_name(i)%str),' &'
      coma = ','
    enddo
    write(iunit,'(9x,a)') ',key    = key    &'
  endif
  write(iunit,'(9x,a)') ')'
  ! compare results
  write(iunit,'(a)')
  do i=1,size(output_files)
    write(iunit,'(a,70a,a,70a,a,a)')"  call compare_files(60,'",       &
                                    gets(output_files(i)%str) , "','", &
                                    gets(output_files(i)%str) , "_ref'", ",1.0d-15)"
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
  write(iunit,'(a)')
  write(iunit,'(a)')"end program"

  close(iunit,status='keep')

  end subroutine generate_test

end module test_generator
