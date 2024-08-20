module dirac_conmat_lao_utils

  use labeled_storage
  use checkpoint

  implicit none

  public conmat_info_t

  private

  type conmat_vector_t
    private
    ! --- variables for a single conmat vector --- !
    real(8), allocatable :: values(:,:)
    integer, public   :: id    = 0
    integer, public   :: ind   = 0 ! index on AOPROPER ("ipcon")
    integer, public   :: irep  = 0
    integer, public   :: trsym = 0
    integer, public   :: iopsy = 0
    integer, public   :: ndim1 = 0 !temp
    integer, public   :: ndim2 = 0 !temp
    character(len=16), public :: lab
    logical           :: data_is_initiated = .false.
    logical           :: data_is_read = .false.
  contains
    private
    procedure, pass, public :: erase => erase_conmatvec_data
    procedure, pass   :: allocated_values
    procedure, pass, public :: alloc_values => alloc_conmatvec_values
    procedure, pass, public :: fill_values
    procedure, pass, public :: get_label
    procedure, pass, public :: get_copy_of_values
    procedure, pass   :: read_values => read_conmatvec_values_from_checkpoint
    procedure, pass, public :: write_to_checkpoint => write_conmatvec_to_checkpoint
    procedure, pass, public :: display => display_conmatvec_data
  end type

  type conmat_info_t
    private
    integer, public :: nr_comp = 0  ! total number of components; for T(1) there will be 3 components (X/Y/Z), for T(2) - 6 components (XX,XY,XZ,YY,YZ,ZZ)
    integer :: epreorth = 0            ! legacy; can be controlled from user input
    logical :: data_is_read= .false.
    type(conmat_vector_t), allocatable, public :: vectors(:)
  contains
    private
    procedure, pass         :: allocated_vectors
    procedure, pass         :: alloc => alloc_conmat_vectors
    procedure, pass, public :: erase => erase_conmat_data
    procedure, pass, public :: display => display_conmat_data
    procedure, pass, public :: init => init_conmat_data_from_DIRAC
    procedure, pass, public :: read_from_DIRAC => read_conmat_data_from_TBMO
    procedure, pass, public :: read_from_checkpoint => read_conmat_data_from_checkpoint
    procedure, pass, public :: write_to_checkpoint => write_conmat_data_to_checkpoint
  end type

  interface flatten
    module procedure flatten_real_2d_to_1d
    module procedure flatten_int_2d_to_1d
  end interface

contains

! ------------------------------------------------------------------------------
  pure function allocated_vectors(self)
    class(conmat_info_t), intent(in):: self
    logical:: allocated_vectors
    allocated_vectors = allocated(self%vectors)
  end function

! ------------------------------------------------------------------------------
  subroutine alloc_conmat_vectors(self)
    class(conmat_info_t), intent(inout):: self
    integer :: nr_comp
    character(len=:),allocatable :: errmsg
    if (self%allocated_vectors()) return
    nr_comp = self%nr_comp
    if (nr_comp > 0) then
      allocate(self%vectors(nr_comp))
    else
      write(errmsg, '(a)') 'error in allocation of conmat vectors'
      call quit(errmsg)
    end if
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_conmat_data(self)
    class(conmat_info_t), intent(out) :: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine fill_values(self, mat, mat_dim)
    class(conmat_vector_t), intent(inout) :: self
    integer, intent(in):: mat_dim
    real(8), intent(in):: mat(mat_dim)
    integer :: n, n1, n2, i, i1, i2

    if (.not.self%allocated_values()) return

    n1=self%ndim1
    n2=self%ndim2
    n=n1*n2
    if (mat_dim < n) stop
    self%values = 0.0d0
    i = 1
    do i2 = 1, n2
      do i1 = 1, n1
        self%values(i1,i2) = mat(i)
        i = i + 1
      end do
    end do
  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_copy_of_values(self, mat, mat_ndim)
    class(conmat_vector_t), intent(in) :: self
    integer, intent(in):: mat_ndim
    real(8), intent(out):: mat(mat_ndim)
    integer :: n, n1, n2, i, i1, i2

    n1=self%ndim1
    n2=self%ndim2
    n=n1*n2
    if (mat_ndim < n) stop
    mat = 0.0d0
    i = 1
    do i2 = 1, n2
      do i1 = 1, n1
        mat(i) = self%values(i1,i2)
        i = i + 1
      end do
    end do
  end subroutine


! ------------------------------------------------------------------------------
  function get_label(self)
    class(conmat_vector_t), intent(in) :: self
    character(len=16) :: get_label
    if (self%data_is_read) then
      get_label = self%lab
    else
      print *, 'conmat_vector data is not read. get_label unknown. stop'
      stop
    end if
  end function


! ------------------------------------------------------------------------------
  subroutine read_conmat_data_from_TBMO(self, icomp)
    class(conmat_info_t), intent(inout) :: self
    integer, intent(in) :: icomp
    integer :: f_unit

    real(8), allocatable    :: D_ai(:, :, :), D_ji(:, :, :)
    integer                 :: file_unit, nr_xyz_comp, which_icomp
    logical                 :: debug_me, dummy , sum_ep_rn


!   read T^B matrix in MO basis to D_con_mo
!   ---------------------------------------

    file_unit = 71

    open(file_unit,                 &
         file   = 'TBMO',           &
         status = 'old',            &
         form   = 'unformatted',    &
         access = 'sequential',     &
         position = 'rewind',       &
         action = 'read')

1   continue

    read(file_unit) nr_xyz_comp
    read(file_unit) which_icomp
    read(file_unit) self%vectors(which_icomp)%irep
    read(file_unit) sum_ep_rn
    read(file_unit) self%vectors(which_icomp)%values

    if (icomp == which_icomp) then
      go to 2
    else
      go to 1
    end if

2   continue

    close(file_unit, status = 'keep')

    write(*, *) 'visual: connection matrix (T^B in MO basis) read from TBMO file'
    write(*, *) 'nr_xyz_comp = ', nr_xyz_comp
    write(*, *) 'component = ', icomp
    write(*, *) 'irep = ', self%vectors(icomp)%irep
    write(*, *) 'sum over e+p shells = ', sum_ep_rn
    write(*, *) 'mat = ', self%vectors(icomp)%values

  end subroutine

! ------------------------------------------------------------------------------
  function open_conmat_file(fname)
    character(*), intent(in) :: fname
    integer :: open_conmat_file
    integer :: funit
    logical :: exists, is_opened
    character(len=:),allocatable :: errmsg

    inquire(file=trim(fname), exist=exists, opened=is_opened)
    if (exists) then
      if (.not.is_opened) then
        open(file = trim(fname), &
             newunit = funit, &
             status = 'old',  &
             access = 'sequential',&
             form   = 'unformatted')
      else
        inquire(file=trim(fname), number=funit)
      end if
    else
      write(errmsg, '(a)') 'file ',trim(fname),' does not exist. Stopping now!'
      call quit(errmsg)
    end if
    open_conmat_file = funit
  end function

! ------------------------------------------------------------------------------
  subroutine flatten_real_2d_to_1d(arr_2d, arr_1d)
    real(8), intent(in) :: arr_2d(:,:)
    real(8), pointer, intent(out) :: arr_1d(:)
    integer :: n

    n = size(arr_2d)
    nullify(arr_1d)
    allocate(arr_1d(n))
    arr_1d(1:n) = 0.0d0
    arr_1d(1:n) = reshape(arr_2d, (/n/)) 

  end subroutine

! ------------------------------------------------------------------------------
  subroutine flatten_int_2d_to_1d(arr_2d, arr_1d)
    integer, intent(in) :: arr_2d(:,:)
    integer, pointer, intent(out) :: arr_1d(:)
    integer :: n

    n = size(arr_2d)
    nullify(arr_1d)
    allocate(arr_1d(n))
    arr_1d(1:n) = 0
    arr_1d(1:n) = reshape(arr_2d, (/n/)) 

  end subroutine

! ------------------------------------------------------------------------------
  pure function allocated_values(self)
    class(conmat_vector_t), intent(in):: self
    logical:: allocated_values
    allocated_values = allocated(self%values)
  end function

! ------------------------------------------------------------------------------
  subroutine alloc_conmatvec_values(self, ndim1, ndim2)
    class(conmat_vector_t), intent(inout):: self
    integer, intent(in) :: ndim1, ndim2
    if (self%allocated_values()) return
    allocate(self%values(ndim1, ndim2))
    self%values = 0.0d0
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_conmatvec_data(self)
    class(conmat_vector_t), intent(out) :: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine display_conmat_data(self)
    class(conmat_info_t), intent(in) :: self
    integer :: i

    write(*,*) 'data on conmat_info_t:'
    write(*,*) ' nr_comp    =', self%nr_comp
    write(*,*) ' epreorth   =', self%epreorth

    if (self%nr_comp > 0) then
      if (self%allocated_vectors()) then
        do i = 1, self%nr_comp
          write(*,*) 'data on conmat_vector_t ', i
          call self%vectors(i)%display()
        end do
      else
        write(*,*) 'conmat vectors NOT ALLOCATED '
      end if
    end if

  end subroutine

! ------------------------------------------------------------------------------
  subroutine display_conmatvec_data(self)
    class(conmat_vector_t), intent(in) :: self
    integer :: i, j, k, n1, n2

    if (self%data_is_read) then
      i = self%id
      write(*,*) 'data on conmat_vector_t ', i
      write(*,*) ' id       =', self%id
      write(*,*) ' ind      =', self%ind
      write(*,*) ' irep     =', self%irep
      write(*,*) ' trsym    =', self%trsym
      write(*,*) ' iopsy    =', self%iopsy
      write(*,*) ' ndim1    =', self%ndim1
      write(*,*) ' ndim2    =', self%ndim2
      write(*,*) ' lab      =', self%lab
      if (self%allocated_values()) then
        do j = 1, n2
          do k = 1, n1
            write(*,*) ' conmat_vector value =', k, j, self%values(k, j)
          end do
        end do
      else
        write(*,*) ' conmat_vector value NOT ALLOCATED'
      end if
    end if

  end subroutine


! ----------------------------------------------------------------------------
  subroutine init_conmat_data_from_DIRAC(self, nr_comp, print_level)
#include "implicit.h"
#include "dcbnmr.h"  ! epreorth
    class(conmat_info_t), intent(inout) :: self
    integer, intent(in):: nr_comp
    integer, intent(in):: print_level

!   non-allocatable data
    if (epreorth) then
      self%epreorth = 1
    else
      self%epreorth = 0
    end if
    self%nr_comp = nr_comp
    self%data_is_read = .true.

!   allocatable data
    call self%alloc()

    if (print_level > 2) then
      call self%display()
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine write_conmat_data_to_checkpoint(self, f_name)
    class(conmat_info_t), intent(in) :: self
    character(*), intent(in) :: f_name
    type(file_info_t) :: conmat_file

    call init_conmat_checkpoint(f_name, conmat_file)
    call write_conmat_nonallocatables_to_checkpoint(self, conmat_file)
    call conmat_checkpoint_close(conmat_file)

  end subroutine

! ------------------------------------------------------------------------------
  subroutine write_conmatvec_to_checkpoint(self, f_name, ivec)
    class(conmat_vector_t), intent(in) :: self
    character(*), intent(in) :: f_name
    integer, intent(in) :: ivec
    type(file_info_t) :: conmat_file

    call init_conmat_checkpoint(f_name, conmat_file)
    call write_conmat_allocatables_to_checkpoint(self, conmat_file, ivec)
    call conmat_checkpoint_close(conmat_file)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_conmat_data_from_checkpoint(self, f_name)
    class(conmat_info_t), intent(out) :: self
    character(*), intent(in) :: f_name
    type(file_info_t) :: conmat_file

    call init_conmat_checkpoint(f_name, conmat_file, valid=.true.)
    call read_conmat_nonallocatables_from_checkpoint(self, conmat_file)
    call read_conmat_allocatables_from_checkpoint(self, conmat_file)
    call conmat_checkpoint_close(conmat_file)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_conmat_nonallocatables_from_checkpoint(conmat_info, conmat_file)
    type(conmat_info_t), intent(inout) :: conmat_info
    type(file_info_t), intent(inout) :: conmat_file
    character(len=80):: group

    if (conmat_info%data_is_read) return

    group = '/result/conmat'
    call lab_read(conmat_file, trim(group)//'/nr_comp', &
                  idata=conmat_info%nr_comp)
    call lab_read(conmat_file, trim(group)//'/epreorth',&
                  idata = conmat_info%epreorth)

    conmat_info%data_is_read = .true.
  end subroutine

! ------------------------------------------------------------------------------
  subroutine read_conmat_allocatables_from_checkpoint(conmat_info, conmat_file)
    type(conmat_info_t), intent(inout) :: conmat_info
    type(file_info_t), intent(inout) :: conmat_file
    character(len=80):: group
    character(len=4):: ich
    real(8), pointer:: tmp_rdata(:)
    integer, pointer:: tmp_idata(:)
    integer:: ivec, nr_vec

    if (.not. conmat_info%data_is_read) return

    if (.not. conmat_info%allocated_vectors()) then
      call conmat_info%alloc()
    end if

    nr_vec = conmat_info%nr_comp
    do ivec = 1, nr_vec
      if (conmat_info%vectors(ivec)%data_is_read) return
      write(ich, '(i4)') ivec
      group = '/result/conmat/vectors/vec_'//trim(adjustl(ich))
      call lab_read(conmat_file, trim(group)//'/id', &
                   idata=conmat_info%vectors(ivec)%id)
      call lab_read(conmat_file, trim(group)//'/ind', &
                   idata=conmat_info%vectors(ivec)%ind)
      call lab_read(conmat_file, trim(group)//'/irep', &
                   idata=conmat_info%vectors(ivec)%irep)
      call lab_read(conmat_file, trim(group)//'/trsym', &
                   idata=conmat_info%vectors(ivec)%trsym)
      call lab_read(conmat_file, trim(group)//'/iopsy', &
                   idata=conmat_info%vectors(ivec)%iopsy)
      call lab_read(conmat_file, trim(group)//'/ndim1', &
                   idata=conmat_info%vectors(ivec)%ndim1)
      call lab_read(conmat_file, trim(group)//'/ndim2', &
                   idata=conmat_info%vectors(ivec)%ndim2)
      call lab_read(conmat_file, trim(group)//'/lab', &
                   sdata=conmat_info%vectors(ivec)%lab)

!     all non-allocatables of vector ivec are read, so alloc and read allocatables
      conmat_info%vectors(ivec)%data_is_read = .true.

      call conmat_info%vectors(ivec)%read_values(conmat_file, group)
    end do

  end subroutine


! ------------------------------------------------------------------------------
  subroutine init_conmat_checkpoint(f_name, conmat_file, valid)
    character(*), intent(in) :: f_name
    type(file_info_t), intent(inout) :: conmat_file
    logical, optional, intent(in) :: valid
    character(len=80) :: top_groups(1)
    logical :: exists
    character(100) :: errmsg

    inquire(file=trim(f_name), exist=exists)
    if (.not. exists) then

      if (present(valid)) then
        if (valid) then
          write(errmsg, *) 'file ',trim(f_name),' does not exist. Stopping now!'
          call quit(errmsg)
        end if
      end if

      write(*,*) " Opening new HDF5 file for conmat data export: ", conmat_file%name
      top_groups(1) = '/result'
      conmat_file%name = trim(f_name)
      conmat_file%type = 2
      conmat_file%status = -1
      call lab_init(conmat_file, top_groups)
    else
      write(*,*) " Opening existing HDF5 file for conmat data export: ", conmat_file%name
      conmat_file%name = trim(f_name)
      conmat_file%type = 2
      conmat_file%status = 0
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine write_conmat_nonallocatables_to_checkpoint(conmat_info, conmat_file)
    type(conmat_info_t), intent(in) :: conmat_info
    type(file_info_t), intent(inout) :: conmat_file
    character(len=80):: group

    group = '/result/conmat'
    call lab_write(conmat_file, trim(group)//'/nr_comp', &
                   idata=conmat_info%nr_comp)
    call lab_write(conmat_file, trim(group)//'/epreorth',&
                   idata = conmat_info%epreorth)

  end subroutine

! ------------------------------------------------------------------------------
  subroutine write_conmat_allocatables_to_checkpoint(conmat_vector, conmat_file, &
                                                     ivec)
    type(conmat_vector_t), intent(in) :: conmat_vector
    type(file_info_t), intent(inout) :: conmat_file
    integer, intent(in) :: ivec
    character(len=80):: group
    character(len=4):: ich
    real(8), pointer:: tmp_rdata(:)

    !if (.not. conmat_info%allocated_vectors()) return

    write(ich, '(i4)') ivec
    group = '/result/conmat/vectors/vec_'//trim(adjustl(ich))
    call lab_write(conmat_file, trim(group)//'/id', &
                 idata=conmat_vector%id)
    call lab_write(conmat_file, trim(group)//'/ind', &
                 idata=conmat_vector%ind)
    call lab_write(conmat_file, trim(group)//'/irep', &
                 idata=conmat_vector%irep)
    call lab_write(conmat_file, trim(group)//'/trsym', &
                 idata=conmat_vector%trsym)
    call lab_write(conmat_file, trim(group)//'/iopsy', &
                 idata=conmat_vector%iopsy)
    call lab_write(conmat_file, trim(group)//'/ndim1', &
                 idata=conmat_vector%ndim1)
    call lab_write(conmat_file, trim(group)//'/ndim2', &
                 idata=conmat_vector%ndim2)
    call lab_write(conmat_file, trim(group)//'/lab', &
                 sdata=conmat_vector%lab)
    if (conmat_vector%allocated_values()) then
      call flatten(conmat_vector%values, tmp_rdata)
      call lab_write(conmat_file, trim(group)//'/values', &
                 rdata=tmp_rdata)
      deallocate(tmp_rdata)
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_conmatvec_values_from_checkpoint(self, conmat_file, group)
    class(conmat_vector_t), intent(inout) :: self
    type(file_info_t), intent(inout) :: conmat_file
    character(len=80), intent(in):: group
    integer :: ndim, ndim1, ndim2
    real(8), pointer:: tmp_rdata(:)

    ndim1=self%ndim1
    ndim2=self%ndim2
    ndim=ndim1*ndim2

    nullify(tmp_rdata)
    allocate(tmp_rdata(ndim))
    tmp_rdata = 0.0d0
    call lab_read(conmat_file, trim(group)//'/values', &
                   rdata=tmp_rdata)

    if (.not.self%allocated_values()) then
      call self%alloc_values(ndim1, ndim2)
    end if

    self%values = reshape(tmp_rdata, (/ndim1,ndim2/))

    deallocate(tmp_rdata)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine conmat_checkpoint_close(conmat_file)
    type(file_info_t), intent(inout) :: conmat_file
    conmat_file%status = 0
    print*, " conmat Checkpoint file closed"
  end subroutine

end module
