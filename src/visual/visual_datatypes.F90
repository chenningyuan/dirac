module visual_datatypes

  implicit none
  private
  public grid_t
  public grid_function_t
  !public grid_function_dirac_t

  !integer, parameter, private:: kind_real = 8

  type:: grid_t
    private
    integer, public:: id = 0
    integer, public:: ndim = 0
    integer, public:: nr_points_total = 0
    integer, public:: action = 0                            ! 0 = generate, 1 = read from file
    integer, public:: typ = 0                               ! 0 = unknown, 1 = regular cubic
    character(len=:), public, allocatable:: file            ! name of a file for grid I/O
    !type(grid_point_t), public, allocatable:: points(:,:)
    real(8), public, allocatable:: points(:,:)
  contains
    private
    procedure, pass, public:: allocated_points
    procedure, pass, public:: erase   => erase_grid
    procedure, pass, public:: display => display_grid
  end type grid_t

  type:: grid_function_t
    private
    integer, public:: id = 0
    character(len=:), public, allocatable:: name
    integer, public:: status = -1                           ! -1 = never used, 0 = initialized, but not assigned, 1 = assigned (ready to work on)
    integer, public:: wf_order = 0                          ! 0 = depends on unperturbed wave function only, >0 requires perturbed wave functions (1 = linear response, 2 = quadratic response, etc.)
    integer, public:: tensor_order = 0                      ! 0 = a scalar field, 1 = a vector field, 2 = a tensor field of rank 2, etc.
    integer, public:: id_assigned_grid = 0                  ! id of assigned grid
    integer, public:: purpose = 0                           ! 0 = visualization, 1 = integration
    character(len=:), public, allocatable:: file            ! name of a file for reading gridfunction from/writing gridfunction to
    real(8), public, allocatable:: values(:,:)
  contains
    private
    procedure, pass, public:: assigned_to_grid
    procedure, pass, public:: allocated_values
    procedure, pass, public:: get_ndim
    procedure, pass, public:: erase   => erase_grid_function
    procedure, pass, public:: display => display_grid_function
  end type grid_function_t

!  type, extends(grid_function_t):: grid_function_dirac_t
!!   this type allows to compose a property density 
!!   in a similar way as is done with 1-el operators in DIRAC
!!   see: 
!    character(len=:), allocatable:: operator_type
!    character(len=:), allocatable:: operator_label
!    real(8):: factor
!    real(8):: comfactor
!  end type grid_function_dirac_t

contains

  pure function allocated_points(self)
    class(grid_t), intent(in):: self
    logical:: allocated_points
    allocated_points = allocated(self%points)
  end function

  pure function allocated_values(self)
    class(grid_function_t), intent(in):: self
    logical:: allocated_values
    allocated_values = allocated(self%values)
  end function

  pure function get_ndim(self)
    class(grid_function_t), intent(in):: self
    integer:: get_ndim
    get_ndim = 3**(self%tensor_order)
  end function

  pure function assigned_to_grid(self)
    class(grid_function_t), intent(in):: self
    logical:: assigned_to_grid
    assigned_to_grid = (self%id_assigned_grid > 0)
  end function

  subroutine erase_grid(self)
    class(grid_t), intent(out):: self
  end subroutine

  subroutine erase_grid_function(self)
    class(grid_function_t), intent(out):: self
  end subroutine

  function get_coord(self, id)
    class(grid_t), target, intent(in):: self
    integer, intent(in):: id
    real(8), pointer:: get_coord(:)

    if (id > self%ndim) stop
    if (.not. self%allocated_points()) return

    nullify(get_coord)
    get_coord(1:self%nr_points_total)  => self%points(id, 1:self%nr_points_total)
  end function

  function get_coord_label(self, id)
    class(grid_t), intent(in):: self
    integer, intent(in):: id
    character(len = 7):: get_coord_label

    if (id > self%ndim) stop
    if (.not. self%allocated_points()) return

    if (id < 10)  then
      write(get_coord_label, "(a6, i1)") "coord_", id
    else
      stop
    end if
  end function

  subroutine display_grid(self)
    class(grid_t), intent(in):: self
    if (self%id == 0) return
    write(*, *) 'id =              ', self%id
    write(*, *) 'ndim =            ', self%ndim
    write(*, *) 'nr_points_total = ', self%nr_points_total
    write(*, *) 'action =          ', self%action
    write(*, *) 'typ =             ', self%typ
    write(*, *) 'file =            ', self%file
   end subroutine

  subroutine display_grid_function(self)
    class(grid_function_t), intent(in):: self
    if (self%id == 0) return
    write(*, *) 'id =              ', self%id
    write(*, *) 'name =            ', self%name
    write(*, *) 'status =          ', self%status
    write(*, *) 'wf_order =        ', self%wf_order
    write(*, *) 'tensor_order      ', self%tensor_order
    write(*, *) 'id_assigned_grid =', self%id_assigned_grid
    write(*, *) 'purpose =         ', self%purpose
    write(*, *) 'file =            ', self%file
   end subroutine

end module
