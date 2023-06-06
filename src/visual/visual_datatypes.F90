module visual_datatypes

  implicit none
  private

  public rs_input_t
  public rs_label_t
  public grid_t
  public grid_function_t
  !public grid_function_dirac_t
  public mol_t
  public ao_t

  !integer, parameter, private:: kind_real = 8

  type mol_t
    private
    integer, public :: n_atoms = 0
    real(8), allocatable, public :: geometry(:,:)
    real(8), allocatable, public :: nuc_charge(:)
    character*4, allocatable, public :: symbols(:)
  contains
    private
    procedure, pass, public:: in_use => mol_in_use
    procedure, pass, public:: erase  => erase_mol
    procedure, public:: display => display_mol
  end type mol_t

  type ao_t
    private
    integer, public:: n_ao = 0
  end type ao_t

  type:: rs_input_t
    private
    integer, public:: print_level = 0
    logical, public:: outpri = .false.   ! if true, print to DIRAC output
    !logical:: small_ao = .false. ! if true, force small ao
    !logical:: occup_ao = .false. ! if true, use orbital strings
    integer, public:: nr_grids = 0
    integer, public:: nr_grid_functions = 0
    type(grid_t), allocatable, public:: grids(:)
    type(grid_function_t), allocatable, public:: grid_functions(:)
!   default values, thresholds, convergence criteria
  contains
    private
!    procedure, public:: setup_grids
    procedure, public:: display => display_rs_setup
    procedure, pass, public:: allocated_grids
    procedure, pass, public:: allocated_grid_functions
    procedure, pass, public:: erase   => erase_rs_input
  end type rs_input_t

  type:: rs_label_t
    private
    character(len=:), allocatable, public:: name
  contains
    procedure, pass, public:: in_use => rs_label_in_use
  end type rs_label_t

  type:: grid_t
    ! gosia todo:
    ! - store units; now this is messed up: 
    !   molecule data (from DIRAC CHECKPOINT) 
    !   and all that the user can provide in input (margin and spacing) is originally in Angstroms, 
    !   but then everything is converted to a.u. (in grid_input_cfg) and stored in a.u.
    private
    integer, public:: id = 0
    integer, public:: ndim = 0                              ! grid dimensionality (0= a single vertex, 1=line, 2=2D, 3=3D)
    integer, public:: nr_points_total = 0
    integer, public:: nr_points(1:3) = 0
    real(8), public:: origin(1:3) = 0                       ! grid origin
    real(8), public:: spacing(1:3) = 0                      ! grid spacing
    integer, public:: status = -1                           ! -1 = never used, 0 = ready to work on
    integer, public:: input = 0                             ! 0 = create, 1 = import
    logical, public:: export = .false.                      ! also export to a separate file
    integer, public:: typ = 0                               ! 0 = unknown, 1 = rectilinear
    integer, public:: weights = 0                           ! 0/1 = points w/wo weights
    real(8), public:: margin = 0.0d0                        ! margin (in Angstrom) to be added around molecule when creating the grid
    character(len=:), public, allocatable:: file_inp        ! name of a file for grid import
    character(len=:), public, allocatable:: file_out        ! name of a file for grid export
    !type(grid_point_t), public, allocatable:: points(:,:)
    real(8), public, allocatable:: points(:,:)
    logical, public:: ndim_known = .false.
    logical, public:: origin_known = .false.
    logical, public:: nr_points_total_known = .false.
    logical, public:: nr_points_known(1:3) = .false.
    logical, public:: spacing_known(1:3) = .false.
  contains
    private
    procedure, pass, public:: allocated_points
    procedure, pass, public:: erase   => erase_grid
    procedure, pass, public:: display => display_grid
    procedure, pass, public:: display_points => display_grid_points
  end type grid_t

  type:: grid_function_t
    private
    integer, public:: id = 0
    character(len=:), public, allocatable:: name
    integer, public:: status = -1                           ! -1 = never used, 0 = initialized, but not assigned, 1 = assigned (ready to work on)
    integer, public:: input = 0                             ! 0 = calculate, 1 = import
    logical, public:: export = .true.                       ! export to file
    integer, public:: wf_order = 0                          ! 0 = depends on unperturbed wave function only, >0 requires perturbed wave functions (1 = linear response, 2 = quadratic response, etc.)
    integer, public:: tensor_order = 0                      ! 0 = a scalar field, 1 = a vector field, 2 = a tensor field of rank 2, etc.
    !integer, public:: max_ao_g_order
    !integer, public:: max_ao_m_order
    integer, public:: id_assigned_grid = 0                  ! id of assigned grid
    integer, public:: purpose = 0                           ! 0 = visualization, 1 = integration
    character(len=:), public, allocatable:: file_inp        ! name of a file for grid function import
    character(len=:), public, allocatable:: file_out        ! name of a file for grid function export
    real(8), public, allocatable:: values(:,:)
  contains
    private
    !procedure, pass, public:: in_use  => grid_function_in_use
    procedure, pass, public:: assigned_to_grid
    procedure, pass, public:: allocated_values
    procedure, pass, public:: get_ndim
    procedure, pass, public:: get_npoints
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

  pure function allocated_grids(self)
    class(rs_input_t), intent(in):: self
    logical:: allocated_grids
    allocated_grids = allocated(self%grids)
  end function

  pure function allocated_grid_functions(self)
    class(rs_input_t), intent(in):: self
    logical:: allocated_grid_functions
    allocated_grid_functions = allocated(self%grid_functions)
  end function

  pure function get_ndim(self)
    class(grid_function_t), intent(in):: self
    integer:: get_ndim
    get_ndim = 3**(self%tensor_order)
  end function

  function get_npoints(self, grid)
    class(grid_function_t), intent(in):: self
    type(grid_t), intent(in):: grid
    integer:: get_npoints
    if (self%id_assigned_grid /= grid%id) return
    get_npoints = grid%nr_points_total
  end function

  pure function rs_label_in_use(self)
    class(rs_label_t), intent(in):: self
    logical:: rs_label_in_use
    rs_label_in_use = (len(self%name) > 0)
  end function

  pure function mol_in_use(self)
    class(mol_t), intent(in):: self
    logical:: mol_in_use
    mol_in_use = ((self%n_atoms) > 0)
  end function

  pure function assigned_to_grid(self)
    class(grid_function_t), intent(in):: self
    logical:: assigned_to_grid
    assigned_to_grid = (self%id_assigned_grid > 0)
  end function

  subroutine erase_rs_input(self)
    class(rs_input_t), intent(out):: self
  end subroutine

  subroutine erase_mol(self)
    class(mol_t), intent(out):: self
  end subroutine

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


  subroutine display_mol(self)
    class(mol_t), intent(in):: self
    integer:: i, j
    write(*, *)
    write(*, *) ' MOL:'
    write(*, *) 'n_atoms =    ', self%n_atoms
    do i = 1, self%n_atoms
      write(*, *) 'position of    ', i, (self%geometry(j,i), j = 1, 3)
      write(*, *) 'nuc_charge of  ', i, self%nuc_charge(i)
      !gosia fixme write(*, *) 'symbol of      ', i, self%symbols(i)
    end do
    write(*, *)
   end subroutine


  subroutine display_rs_setup(self)
    class(rs_input_t), intent(in):: self
    write(*, *)
    write(*, *) ' RS SETUP:'
    write(*, *) 'print_level =         ', self%print_level
    write(*, *) 'outpri =              ', self%outpri
    write(*, *) 'nr_grids =            ', self%nr_grids
    write(*, *) 'nr_grid_functions =   ', self%nr_grid_functions
    write(*, *)
   end subroutine

  subroutine display_grid(self,my_rank)
    class(grid_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    logical:: with_rank_info
    integer::i, j

    if (self%id == 0) return

    with_rank_info = .false.
    if (present(my_rank)) with_rank_info = .true.

    write(*, *)
    write(*, *) ' GRID:'
    if (with_rank_info) then
      write(*, *) 'id =              ', self%id, my_rank
      write(*, *) 'ndim =            ', self%ndim, my_rank 
      write(*, *) 'nr_points       = ', (self%nr_points(i), i=1,self%ndim), my_rank
      write(*, *) 'nr_points_total = ', self%nr_points_total, my_rank
      write(*, *) 'origin [a.u.]     = ', (self%origin(i), i=1,self%ndim), my_rank
      write(*, *) 'spacing [a.u.]    = ', (self%spacing(i), i=1,self%ndim), my_rank
      write(*, *) 'status =          ', self%status, my_rank
      write(*, *) 'input =           ', self%input, my_rank
      write(*, *) 'export =          ', self%export, my_rank
      write(*, *) 'typ =             ', self%typ, my_rank
      write(*, *) 'weights =         ', self%weights, my_rank
      write(*, *) 'margin =          ', self%margin, my_rank
      write(*, *) 'file_inp =        ', self%file_inp, my_rank
      write(*, *) 'file_out =        ', self%file_out, my_rank
      write(*, *) 'ndim_known =      ', self%ndim_known, my_rank
      write(*, *) 'origin_known =    ', self%origin_known, my_rank
      write(*, *) 'nr_points_total_known = ', self%nr_points_total_known, my_rank
      write(*, *) 'nr_points_known(1:3) = ', (self%nr_points_known(i), i=1,self%ndim), my_rank
      write(*, *) 'spacing_known(1:3) = ', (self%spacing_known(i), i=1,self%ndim), my_rank
    else
      write(*, *) 'id =              ', self%id
      write(*, *) 'ndim =            ', self%ndim
      write(*, *) 'nr_points       = ', (self%nr_points(i), i=1,self%ndim)
      write(*, *) 'nr_points_total = ', self%nr_points_total
      write(*, *) 'origin [a.u.]     = ', (self%origin(i), i=1,self%ndim)
      write(*, *) 'spacing [a.u.]    = ', (self%spacing(i), i=1,self%ndim)
      write(*, *) 'status =          ', self%status
      write(*, *) 'input =           ', self%input
      write(*, *) 'export =          ', self%export
      write(*, *) 'typ =             ', self%typ
      write(*, *) 'weights =         ', self%weights
      write(*, *) 'margin =          ', self%margin
      write(*, *) 'file_inp =        ', self%file_inp
      write(*, *) 'file_out =        ', self%file_out
      write(*, *) 'ndim_known =      ', self%ndim_known
      write(*, *) 'origin_known =    ', self%origin_known
      write(*, *) 'nr_points_total_known = ', self%nr_points_total_known
      write(*, *) 'nr_points_known(1:3) = ', (self%nr_points_known(i), i=1,self%ndim)
      write(*, *) 'spacing_known(1:3) = ', (self%spacing_known(i), i=1,self%ndim)
    end if
    write(*, *)
   end subroutine

  subroutine display_grid_function(self,my_rank)
    class(grid_function_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    logical:: with_rank_info

    if (self%id == 0) return

    with_rank_info = .false.
    if (present(my_rank)) with_rank_info = .true.

    write(*, *)
    write(*, *) ' GRID FUNCTION:'
    if (with_rank_info) then
      write(*, *) 'id =              ', self%id, my_rank
      if (allocated(self%name)) write(*, *) 'name =            ', self%name, my_rank
      write(*, *) 'status =          ', self%status, my_rank
      write(*, *) 'input =           ', self%input, my_rank
      write(*, *) 'export =          ', self%export, my_rank
      write(*, *) 'wf_order =        ', self%wf_order, my_rank
      write(*, *) 'tensor_order      ', self%tensor_order, my_rank
      write(*, *) 'id_assigned_grid =', self%id_assigned_grid, my_rank
      write(*, *) 'purpose =         ', self%purpose, my_rank
      if (allocated(self%file_inp)) write(*, *) 'file_inp =        ', self%file_inp, my_rank
      if (allocated(self%file_out)) write(*, *) 'file_out =        ', self%file_out, my_rank
    else
      write(*, *) 'id =              ', self%id
      if (allocated(self%name)) write(*, *) 'name =            ', self%name
      write(*, *) 'status =          ', self%status
      write(*, *) 'input =           ', self%input
      write(*, *) 'export =          ', self%export
      write(*, *) 'wf_order =        ', self%wf_order
      write(*, *) 'tensor_order      ', self%tensor_order
      write(*, *) 'id_assigned_grid =', self%id_assigned_grid
      write(*, *) 'purpose =         ', self%purpose
      if (allocated(self%file_inp)) write(*, *) 'file_inp =        ', self%file_inp
      if (allocated(self%file_out)) write(*, *) 'file_out =        ', self%file_out
    end if
    write(*, *)
  end subroutine

  subroutine display_grid_points(self,my_rank)
    class(grid_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    logical:: with_rank_info
    integer::i

    if (self%id == 0) return
    if (.not. self%allocated_points()) return

    with_rank_info = .false.
    if (present(my_rank)) with_rank_info = .true.

    write(*, *)
    if (with_rank_info) then
      write(*, *) ' GRID points (grid id = ):', self%id, my_rank
      call print_rs_2d(self%points, self%ndim, self%nr_points_total, my_rank)
    else
      write(*, *) ' GRID points (grid id = ):', self%id
      call print_rs_2d(self%points, self%ndim, self%nr_points_total)
    end if
  end subroutine

  subroutine print_rs_2d(array, n, m, my_rank)
    real(8), intent(in):: array(n, m)
    integer, intent(in):: n, m
    integer, intent(in), optional:: my_rank
    logical:: with_rank_info
    integer :: i, j

    with_rank_info = .false.
    if (present(my_rank)) with_rank_info = .true.

    if (with_rank_info) then
      do j = 1, m
        write(*, *) my_rank, (array(i, j), i = 1, n)
      end do
    else
      do j = 1, m
        write(*, *) (array(i, j), i = 1, n)
      end do
    end if
  end subroutine

end module
