module visual_datatypes

  implicit none
  private

  public rs_input_t
  public rs_label_t
  public rs_wf_t
  public grid_t
  public grid_function_t
  public mol_t
  public cmo_t

  type mol_t
    private
    integer, public :: n_atoms = 0
    real(8), allocatable, public :: geometry(:,:)
    real(8), allocatable, public :: nuc_charge(:)
    character*4, allocatable, public :: symbols(:)
!   temp (fix for esp)
    integer, public :: nucdep = 0
    real(8), allocatable, public :: geometry_dep(:,:)
    real(8), allocatable, public :: nuc_charge_dep(:)
  contains
    private
    procedure, pass, public:: in_use => mol_in_use
    procedure, pass, public:: erase  => erase_mol
    procedure, public:: display => display_mol
  end type mol_t

  type cmo_t
    private
    integer, public:: n_ao(2) = 0                   ! number of AOs in each fermion symmetry
    integer, public:: n_mo(2) = 0                   ! number of all (positronic+electronic) MOs in each fermion symmetry
    integer, public:: n_ip(2) = 0                   ! number of inactive positive energy MOs in each fermion symmetry
    integer, public:: n_po(2) = 0                   ! number of positronic MOs in each fermion symmetry
    integer, public:: n_fsym  = 1                   ! Number of fermion symmetries (no inversion 1, inversion 2); default=1
    integer, public:: nz   = 4                      ! number of quaternion units; default=4
    real(8), allocatable, public:: coeffs(:)        ! MO coefficients
    real(8), allocatable, public:: eigenvalues(:)   ! orbital energies
    real(8), allocatable, public:: occupations(:)   ! occupation vector: 1.0 fully occupied, 0.0 empty
    integer, allocatable, public:: symmetry(:)      ! symmetry labels
    integer, allocatable, public:: shell_id(:)      ! shell -1: closed, 0: virtual, 1: open shell
    integer, public :: ioffs(6) = 0                 ! MO offsets
    logical, public :: data_is_read  = .false.
    logical, public :: data_is_verified  = .false.
  contains
    private
    procedure, pass, public:: allocated_coeffs
    procedure, pass, public:: allocated_eigenvalues
    procedure, pass, public:: allocated_occupations
    procedure, pass, public:: allocated_symmetry
    procedure, pass, public:: allocated_shell_id
    procedure, pass, public:: alloc_all => allocate_cmo_all
    procedure, pass, public:: erase   => erase_cmo
    procedure, pass, public:: display => display_cmo
  end type cmo_t

  type:: grid_t
    ! todo
    ! - add information about units; now:
    !   molecule data (from DIRAC CHECKPOINT) and all that the user can provide 
    !   in input (margin and spacing) is originally in Angstroms;
    !   then, everything in VISUAL is converted to a.u. (in grid_input_cfg) and stored in a.u.
    private
    integer, public:: ndim = 3 ! grid dimensionality; all grids are embedded in 3D space
    integer, public:: id = 0
    integer, public:: nr_points_total = 0
    integer, public:: nr_points(1:3) = 0
    real(8), public:: origin(1:3) = 0.0d0              ! grid origin [a.u.]
    real(8), public:: spacing(1:3) = 0.0d0             ! grid spacing [a.u.]
    real(8), public:: extents(0:2,1:3) = 0.0d0         ! grid extents [a.u.]
    real(8), public:: radial_length = 0.0d0            ! radial length [a.u.] ->todo fixme make sure it is in a.u....
    integer, public:: radial_nr_steps = 0              ! number of steps for radial grids
    integer, public:: radial_nr_points_angular = 0     ! number of angular points for radial grids
    integer, public:: radial_lebedev_l = 0             ! Lebedev grid L (now, set up to = num_grid_cfg_angint)
    integer, public:: status = -1                      ! -1 = never used, 0 = ready to work on
    integer, public:: input = 0                        ! 0=create, 1=import from hdf5 file, -1=import from 'numerical_grid', 2='list', 3='line', 4='2d', -4='2d_int', 5='3d'
    integer, public:: typ = 0                          ! 0 = unknown, 1 = rectilinear uniform, 2 = radial
    real(8), public:: margin = 4.0d0                   ! margin [a.u.] to be added around molecule when creating the grid
    integer, public:: integration_order_2d = 0         ! order of integration in 2d
    character(len=:), public,allocatable:: import      ! import from file; specify file format: 'hdf5', 'csv', 'txt'
    character(len=:), public,allocatable:: export      ! also export to a separate file; specify file format: 'hdf5', 'csv', 'txt'
    character(len=:), public, allocatable:: file_inp   ! name of a file for grid import
    character(len=:), public, allocatable:: file_out   ! name of a file for grid export
    real(8), public, allocatable:: radial_coord(:)
    real(8), public, allocatable:: points(:,:)
    real(8), public, allocatable:: weights(:)
    logical, public:: origin_known = .false.
    logical, public:: extents_known = .false.
    logical, public:: nr_points_total_known = .false.
    logical, public:: nr_points_known(1:3) = .false.
    logical, public:: spacing_known(1:3) = .false.
  contains
    private
    procedure, pass, public:: allocated_points
    procedure, pass, public:: allocated_radial_coord
    procedure, pass, public:: allocated_weights
    procedure, pass, public:: requires_weights
    !procedure, pass, public:: get_grid_size
    procedure, pass, public:: get_ndim => get_grid_ndim
    procedure, pass, public:: erase   => erase_grid
    procedure, pass, public:: printout => printout_grid_info
    procedure, pass, public:: display => display_grid
    procedure, pass, public:: display_points => display_grid_points
    procedure, pass, public:: display_weights => display_grid_weights
  end type grid_t

  type:: grid_function_t
    private
    integer, public:: id = 0
    character(len=:), public, allocatable:: name
    integer, public:: status = -1                           ! -1 = never used, 0 = initialized, but not assigned, 1 = assigned (ready to work on)
    logical, public:: outpri = .false.                      ! if true, print to DIRAC output
    integer, public:: input = 0                             ! 0 = calculate, 1 = import
    integer, public:: nuc = 0                               ! id of a nucleus (for grid functions specific to atomic centers)
    integer, public:: wf_order = 0                          ! 0 = depends on unperturbed wave function only, >0 requires perturbed wave functions (1 = linear response, 2 = quadratic response, etc.)
    integer, public:: rspvec_id = 0                         ! now: record ID on PAMXVC file
    integer, public:: lao = 0                               ! 0 = no LAOs, 1 = LAOs (all LAO-contributions)
    integer, public:: lao_skip(1:3) = 0                     ! 0/1 = do not/do skip LAO contributions; ('response', 'direct', 'reorthonormalization')
    character(len=1), public:: lao_comp = ' '               ! temp: 'X', 'Y', 'Z'
    logical, public:: nics = .false.                        ! NICS-like shift for magnetic property densities
    real(8), public:: nics_origin(1:3) = 0.0d0              ! NICS-like shift for magnetic property densities
    integer, public:: tensor_order = 0                      ! 0 = a scalar field, 1 = a vector field, 2 = a tensor field of rank 2, etc.
    integer, public:: max_ao_g_order = 0 !legacy
    integer, public:: max_ao_m_order = 0 !legacy
    integer, public:: id_assigned_grid = 0                  ! id of assigned grid
    integer, public:: purpose = 0                           ! 0 = visualization, 1 = integration
    integer, public:: modify = 0                            ! 0/1 do not/do modify grid function
    integer, public:: radial_power = 0
    integer, public:: cartesian_power(1:3) = 0.0d0
    real(8), public:: up_scale = 1.0d0
    real(8), public:: down_scale = 1.0d0
    real(8), public:: bed_wave_vector(1:3) = 0.0d0
    real(8), public:: bed_pol_vector(1:3) = 0.0d0
    real(8), public:: bed_freq = 0.0d0
    logical, public:: bedfix_known= .false.
    character(len=:), public, allocatable:: rspvec_file     ! name of a file with response parameters
    character(len=:), public, allocatable:: gf_type         ! must be one of: "unperturbed", "response", "composite"
    character(len=:), public, allocatable:: export          ! export to file; specify file format: 'hdf5', 'csv', 'txt', 'cube'
    character(len=:), public, allocatable:: file_inp        ! name of a file for grid function import
    character(len=:), public, allocatable:: file_out        ! name of a file for grid function export
    real(8), public, allocatable:: integrated_values(:)
    real(8), public, allocatable:: values(:,:)
    logical, public:: composite = .false. ! is a sum of a few contributions (e.g., depends on a few density matrices)
  contains
    private
    procedure, pass, public:: assigned_to_grid
    procedure, pass, public:: allocated_integrated_values
    procedure, pass, public:: allocated_values
    procedure, pass, public:: get_ndim => get_grid_function_ndim
    procedure, pass, public:: get_npoints
    procedure, pass, public:: should_depend_on_LAOs
    procedure, pass, public:: erase   => erase_grid_function
    procedure, pass, public:: printout => printout_grid_function_info
    procedure, pass, public:: display => display_grid_function
  end type grid_function_t

  type:: rs_input_t
    private
    integer, public:: print_level = 0
    integer, public:: print_unit  = 6
    logical, public:: set_occupation = .false.      ! if true, modify orbital occupations
    logical, public:: force_small_ao = .false.      ! if true, force evaluation of small component basis functions
    logical, public:: set_cvalue     = .false.      ! if true, modify cvalue
    real(8), public:: cvalue = 0.0d0                ! cvalue in VISUAL
    real(8), public:: gauge_origin(1:3) = 0.0d0     ! gauge origin in VISUAL
    real(8), allocatable, public:: occupation(:,:)  ! occupation vector
    integer, public:: nr_grids = 0
    integer, public:: nr_grid_functions = 0
    type(grid_t), allocatable, public:: grids(:)
    type(grid_function_t), allocatable, public:: grid_functions(:)
  contains
    private
    procedure, public:: display => display_rs_setup
    procedure, pass, public:: allocated_grids
    procedure, pass, public:: allocated_grid_functions
    procedure, pass, public:: allocated_occupation => allocated_occupation_from_input
    procedure, pass, public:: erase   => erase_rs_input
  end type rs_input_t

  type:: rs_label_t
    private
    character(len=:), allocatable, public:: name
  contains
    procedure, pass, public:: in_use => rs_label_in_use
  end type rs_label_t

  type:: rs_wf_t
    ! store information on wavefunction
    ! (read from 'wavefunctions' sections on CHECKPOINT.h5)
    private
    character(len=:), public, allocatable:: hamiltonian
    character(len=:), public, allocatable:: xc_functional
  contains
    procedure, public:: display => display_rs_wf
  end type rs_wf_t

contains

! ------------------------------------------------------------------------------
  pure function allocated_eigenvalues(self)
    class(cmo_t), intent(in):: self
    logical:: allocated_eigenvalues
    allocated_eigenvalues = allocated(self%eigenvalues)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_coeffs(self)
    class(cmo_t), intent(in):: self
    logical:: allocated_coeffs
    allocated_coeffs = allocated(self%coeffs)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_symmetry(self)
    class(cmo_t), intent(in):: self
    logical:: allocated_symmetry
    allocated_symmetry = allocated(self%symmetry)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_shell_id(self)
    class(cmo_t), intent(in):: self
    logical:: allocated_shell_id
    allocated_shell_id = allocated(self%shell_id)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_occupations(self)
    class(cmo_t), intent(in):: self
    logical:: allocated_occupations
    allocated_occupations = allocated(self%occupations)
  end function

! ------------------------------------------------------------------------------
  subroutine erase_cmo(self)
    class(cmo_t), intent(out):: self
  end subroutine

! ------------------------------------------------------------------------------
  pure function allocated_points(self)
    class(grid_t), intent(in):: self
    logical:: allocated_points
    allocated_points = allocated(self%points)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_radial_coord(self)
    class(grid_t), intent(in):: self
    logical:: allocated_radial_coord
    allocated_radial_coord = allocated(self%radial_coord)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_weights(self)
    class(grid_t), intent(in):: self
    logical:: allocated_weights
    allocated_weights = allocated(self%weights)
  end function

! ------------------------------------------------------------------------------
  pure function requires_weights(self)
    class(grid_t), intent(in):: self
    logical:: requires_weights
    requires_weights = (self%typ == 2    & ! radial grid 
                  .or.  self%input == -1 & ! integration on dft grid ("numerical_grid")
                  .or.  self%input == -4 & ! integration on 2d grid ("2d_int")
                  .or.  self%input == -5)  ! integration on 3d grid ("3d_int")
  end function

! ------------------------------------------------------------------------------
  pure function get_grid_ndim(self)
    class(grid_t), intent(in):: self
    integer:: get_grid_ndim
    get_grid_ndim = self%ndim
  end function

! ------------------------------------------------------------------------------
  pure function allocated_values(self)
    class(grid_function_t), intent(in):: self
    logical:: allocated_values
    allocated_values = allocated(self%values)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_integrated_values(self)
    class(grid_function_t), intent(in):: self
    logical:: allocated_integrated_values
    allocated_integrated_values = allocated(self%integrated_values)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_grids(self)
    class(rs_input_t), intent(in):: self
    logical:: allocated_grids
    allocated_grids = allocated(self%grids)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_grid_functions(self)
    class(rs_input_t), intent(in):: self
    logical:: allocated_grid_functions
    allocated_grid_functions = allocated(self%grid_functions)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_occupation_from_input(self)
    class(rs_input_t), intent(in):: self
    logical:: allocated_occupation_from_input
    allocated_occupation_from_input = allocated(self%occupation)
  end function

! ------------------------------------------------------------------------------
  pure function get_grid_function_ndim(self)
    class(grid_function_t), intent(in):: self
    integer:: get_grid_function_ndim
    get_grid_function_ndim = 3**(self%tensor_order)
  end function

! ------------------------------------------------------------------------------
  function get_npoints(self, grid)
    class(grid_function_t), intent(in):: self
    type(grid_t), intent(in):: grid
    integer:: get_npoints
    if (self%id_assigned_grid /= grid%id) return
    get_npoints = grid%nr_points_total
  end function

! ------------------------------------------------------------------------------
  pure function rs_label_in_use(self)
    class(rs_label_t), intent(in):: self
    logical:: rs_label_in_use
    rs_label_in_use = (len(self%name) > 0)
  end function

! ------------------------------------------------------------------------------
  pure function mol_in_use(self)
    class(mol_t), intent(in):: self
    logical:: mol_in_use
    mol_in_use = ((self%n_atoms) > 0)
  end function

! ------------------------------------------------------------------------------
  pure function assigned_to_grid(self)
    class(grid_function_t), intent(in):: self
    logical:: assigned_to_grid
    assigned_to_grid = (self%id_assigned_grid > 0)
  end function

! ------------------------------------------------------------------------------
  pure function should_depend_on_LAOs(self, lao_mode)
    class(grid_function_t), intent(in):: self
    logical, intent(in):: lao_mode
    logical:: should_depend_on_LAOs
    should_depend_on_LAOs = (self%lao > 0 .and. (lao_mode))
  end function

! ------------------------------------------------------------------------------
  subroutine erase_rs_input(self)
    class(rs_input_t), intent(out):: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_mol(self)
    class(mol_t), intent(out):: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_grid(self)
    class(grid_t), intent(out):: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_grid_function(self)
    class(grid_function_t), intent(out):: self
  end subroutine

! ------------------------------------------------------------------------------
  function get_coord(self, id)
    class(grid_t), target, intent(in):: self
    integer, intent(in):: id
    real(8), pointer:: get_coord(:)

    if (id > self%get_ndim()) stop
    if (.not. self%allocated_points()) return

    nullify(get_coord)
    get_coord(1:self%nr_points_total)  => self%points(id, 1:self%nr_points_total)
  end function

! ------------------------------------------------------------------------------
  function get_coord_label(self, id)
    class(grid_t), intent(in):: self
    integer, intent(in):: id
    character(len = 7):: get_coord_label

    if (id > self%get_ndim()) stop
    if (.not. self%allocated_points()) return

    if (id < 10)  then
      write(get_coord_label, "(a6, i1)") "coord_", id
    else
      stop
    end if
  end function

! ------------------------------------------------------------------------------
  subroutine allocate_cmo_all(self)
    class(cmo_t), intent(inout):: self
    integer :: nmo,ncoeff

    if (self%allocated_coeffs() .or.&
        self%allocated_eigenvalues() .or.&
        self%allocated_occupations() .or.&
        self%allocated_shell_id() .or.&
        self%allocated_symmetry()) return

    nmo = sum(self%n_mo)
    ncoeff = self%nz*sum(self%n_ao*self%n_mo)

    allocate(self%coeffs(ncoeff))
    self%coeffs(1:ncoeff) = 0.0d0

    allocate(self%eigenvalues(nmo))
    self%eigenvalues(1:nmo) = 0.0d0

    allocate(self%occupations(nmo))
    self%occupations(1:nmo) = 0.0d0

    allocate(self%shell_id(nmo))
    self%shell_id(1:nmo) = 0

    allocate(self%symmetry(nmo))
    self%symmetry(1:nmo) = 0
  end subroutine


! ------------------------------------------------------------------------------
  subroutine display_cmo(self)
    class(cmo_t), intent(in):: self
    integer:: i,nmo,ncoeff

    nmo = sum(self%n_mo)
    ncoeff = self%nz*sum(self%n_ao*self%n_mo)

    write(*, *)
    write(*, *) 'CMO in VISUAL:'
    write(*, *) 'n_ao =         ', self%n_ao
    write(*, *) 'n_mo =         ', self%n_mo
    write(*, *) 'n_ip =         ', self%n_ip
    write(*, *) 'n_po =         ', self%n_po
    write(*, *) 'n_fsym =       ', self%n_fsym
    write(*, *) 'nz   =         ', self%nz
    if (self%allocated_coeffs()) write(*,'(a,es16.8e2)') 'coeffs = ', (self%coeffs(i), i=1,ncoeff)
    if (self%allocated_eigenvalues()) write(*,'(a,es16.8e2)') 'eigenvalues = ', (self%eigenvalues(i), i=1,nmo)
    if (self%allocated_occupations()) write(*,*) 'occupations = ', (self%occupations(i), i=1,nmo)
    if (self%allocated_symmetry()) write(*,*) 'symmetry = ', (self%symmetry(i), i=1,nmo)
    if (self%allocated_shell_id()) write(*,*) 'shell_id = ', (self%shell_id(i), i=1,nmo)
    write(*, *)
  end subroutine


! ------------------------------------------------------------------------------
  subroutine display_mol(self)
    class(mol_t), intent(in):: self
    integer:: i, j
    write(*, *)
    write(*, *) ' MOL:'
    write(*, *) 'n_atoms =    ', self%n_atoms
    write(*, *) 'nucdep  =    ', self%nucdep
    do i = 1, self%n_atoms
      write(*, *) 'position of    ', i, (self%geometry(j,i), j = 1, 3)
      write(*, *) 'nuc_charge of  ', i, self%nuc_charge(i)
      !write(*, *) 'symbol of      ', i, self%symbols(i)
    end do
    do i = 1, self%nucdep
      write(*, *) '(symmetry-dep) position of    ', i, (self%geometry_dep(j,i), j = 1, 3)
      write(*, *) '(symmetry-dep) nuc_charge of  ', i, self%nuc_charge_dep(i)
    end do
    write(*, *)
   end subroutine


! ------------------------------------------------------------------------------
  subroutine display_rs_setup(self,my_rank)
    class(rs_input_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    logical:: with_rank_info
    integer :: i, j

    with_rank_info = .false.
    if (present(my_rank)) with_rank_info = .true.

    write(*, *)
    if (with_rank_info) then
      write(*, *) ' VISUAL: rs_input on node ', my_rank
    else
      write(*, *) ' VISUAL: rs_input'
    end if
    write(*, *) 'print_level =         ', self%print_level
    write(*, *) 'print_unit  =         ', self%print_unit
    write(*, *) 'set_occupation =      ', self%set_occupation
    write(*, *) 'force_small_ao =      ', self%force_small_ao
    write(*, *) 'cvalue =              ', self%cvalue
    write(*, *) 'gauge_origin =        ', self%gauge_origin
    write(*, *) 'set_cvalue     =      ', self%set_cvalue
    write(*, *) 'nr_grids =            ', self%nr_grids
    write(*, *) 'nr_grid_functions =   ', self%nr_grid_functions
    if (self%allocated_occupation()) then  
      do j = 1, size(self%occupation,2)
        do i = 1, size(self%occupation,1)
          write(*,*) 'occupation   = ', self%occupation(i,j)
        end do
      end do
    end if
    write(*, *)
   end subroutine


! ------------------------------------------------------------------------------
  subroutine display_rs_wf(self,my_rank)
    class(rs_wf_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    logical:: with_rank_info

    with_rank_info = .false.
    if (present(my_rank)) with_rank_info = .true.

    write(*, *)

    write(*, *)
    if (with_rank_info) then
      write(*, *) ' VISUAL: rs_wf on node ', my_rank
    else
      write(*, *) ' VISUAL: rs_wf'
    end if
    write(*, *) ' Hamiltonian and method used for the wavefunction optimization:'
    if (allocated(self%hamiltonian)) &
      write(*, *) 'hamiltonian =   ', trim(self%hamiltonian)
    if (allocated(self%xc_functional)) &
      write(*, *) 'xc_functional = ', self%xc_functional
    write(*, *)
   end subroutine


! ------------------------------------------------------------------------------
  subroutine printout_grid_info(self)
    class(grid_t), intent(in):: self
    integer::i, j

    if (self%id == 0) return

    select case (self%input)

      case (1)
        write(*, '(2x,a)') 'Grid is read from file: ', self%file_inp

      case (-1)
        write(*, '(2x,a)') 'Grid points are read from the "numerical_grid" file '//&
                           'from DFT calculations.'

      case (2)
        write(*, '(2x,a)') 'Coordinates of grid points '//&
                           'are provided in the user input (.LIST)'
        write(*, '(4x,a,i8)') '* Total number of points: ', self%nr_points_total

      case (3)
        write(*, '(2x,a)') 'Grid points are generated as equidistant points, '//&
                           'on line defined in the user input (.LINE): '
        write(*, '(4x,a,i8)') '* The line is defined by two points: '
        write(*, '(6x,a,3f13.4)') 'start: ', self%extents(0,1:3)
        write(*, '(6x,a,3f13.4)') 'end:   ', self%extents(1,1:3)
        write(*, '(4x,a,i8)') '* Total number of points: ', self%nr_points_total

      case (4, -4)
        write(*, '(2x,a)') 'Grid points are generated to span a 2D plane, '//&
                           'defined in the user input (.2D or .2D_INT): '
        write(*, '(4x,a)') '* The plane is spanned by three points: '
        write(*, '(6x,a,3f13.4)') '"origin": ', self%extents(0,1:3)
        write(*, '(6x,a,3f13.4)') '"right" : ', self%extents(1,1:3)
        write(*, '(6x,a,3f13.4)') '"top"   : ', self%extents(2,1:3)
        write(*, '(4x,a,i8)') '* The number of points to "right": ', self%nr_points(1)
        write(*, '(4x,a,i8)') '* The number of points to "top"  : ', self%nr_points(2)
        write(*, '(4x,a,i8)') '* Total number of points:       ', self%nr_points_total

      case (0, 5, -5)
        if (self%input == 0) then
          write(*, '(2x,a)') 'Grid points are generated to span a 3D rectilinear cube, '//&
                             'defined in the user input (.3D or .3D_INT): '
        else
          write(*, '(2x,a)') 'Grid points are generated to span a 3D rectilinear cube, '//&
                             'defined in the grid specification line: '
        end if
        write(*, '(4x,a)') '* Cube parameters: '
        write(*, '(6x,a)') '* "x" direction: '
        write(*,'(8x,a,i8)') 'nr points: ', self%nr_points(1)
        write(*,'(8x,a,3f13.4)') 'spacing: ', self%spacing(1)
        write(*,'(8x,a,f13.4,a,f13.4,a)') 'extend from: ', self%extents(0,1), ' to ', self%extents(1,1), ' on x axis'
        write(*, '(6x,a)') '* "y" direction: '
        write(*,'(8x,a,i8)') 'nr points: ', self%nr_points(2)
        write(*,'(8x,a,3f13.4)') 'spacing: ', self%spacing(2)
        write(*,'(8x,a,f13.4,a,f13.4,a)') 'extend from: ', self%extents(0,2), ' to ', self%extents(1,2), ' on y axis'
        write(*, '(6x,a)') '* "z" direction: '
        write(*,'(8x,a,i8)') 'nr points: ', self%nr_points(3)
        write(*,'(8x,a,3f13.4)') 'spacing: ', self%spacing(3)
        write(*,'(8x,a,f13.4,a,f13.4,a)') 'extend from: ', self%extents(0,3), ' to ', self%extents(1,3), ' on z axis'
        write(*, '(4x,a,i8)') '* Total number of points:       ', self%nr_points_total

      case (6)
        write(*, '(2x,a)') 'Grid points are generated to span radial grid, '//&
                           'defined in the user input (.RADIAL): '
        write(*, '(4x,a,3f13.4)') '* Grid origin:                  ', self%origin(1:3)
        write(*, '(4x,a,i8)') '* Lebedev angular grid L:       ', self%radial_lebedev_l 
        write(*, '(4x,a,i8)') '* The number of angular points: ', self%radial_nr_points_angular
        write(*, '(4x,a,i8)') '* Total number of points:       ', self%nr_points_total

    end select

      if (allocated(self%export)) write(*, '(2x,a)') 'Grid will be exported; export format: ', self%export
      if (allocated(self%file_out)) write(*, '(2x,a)') 'Grid will be exported to file: ', self%file_out

    !if (self%purpose == 1) then
    !  if (self%typ == 1) then
    !    write(*, *) 'Integration on rectilinear grid: '
    !  else if (self%typ == 2) then
    !    write(*, *) 'Integration on radial grid: '
    !  end if
    !end if
    write(*, *)
   end subroutine


! ------------------------------------------------------------------------------
  subroutine printout_grid_function_info(self)
    class(grid_function_t), intent(in):: self
    integer::i, j

    if (self%id == 0) return

    write(*, '(2x,a,i8)') 'Grid function is a tensor field of order ', self%tensor_order

    if (self%purpose == 1) then
      write(*, '(2x,a,i8)') 'Grid function will be integrated on grid ', self%id_assigned_grid
    else if (self%purpose == 0) then
      write(*, '(2x,a,i8)') 'Grid function will be calculated in vertices of grid', self%id_assigned_grid
    end if

    if (self%modify == 1) then
      write(*, '(2x,a)') 'Grid function values will be modified (pointwise): '
      write(*, '(4x,a,f13.4)') '* upscaled by: ', self%up_scale
      write(*, '(4x,a,f13.4)') '* downscaled by: ', self%down_scale
      write(*, '(4x,a,3i3)') '* multiplied with cartesian powers: ', self%cartesian_power(1:3)
      write(*, '(4x,a,i3)') '* multiplied with radial power: ', self%radial_power
    end if

    if (self%wf_order > 0) then
      write(*, '(2x,a)') 'Grid function depends on the perturbed electron density:'
      write(*, '(4x,a,i3)') '* perturbation order: ', self%wf_order
      write(*, '(4x,a)') '* File with response vectors data: ', self%rspvec_file
      write(*, '(4x,a,i3)') '* ID of a response vector: ', self%rspvec_id
    else
      write(*, '(2x,a)') 'Grid function depends on the unperturbed electron density'
    end if

    if (self%lao == 1) then
      write(*, '(2x,a)') 'London atomic orbitals (LAOs) are used in the formulation of grid function'
      if (self%lao_skip(1) == 1) &
        write(*, '(4x,a)') '* LAOs: response contribution to the grid function is neglected'
      if (self%lao_skip(2) == 1) &
        write(*, '(4x,a)') '* LAOs: "direct" contribution to the grid function is neglected'
      if (self%lao_skip(3) == 1) &
        write(*, '(4x,a)') '* LAOs: reorthonormalization contribution to the grid function is neglected'
    end if

      if (allocated(self%export)) write(*, '(2x,a)') 'Grid function will be exported; export format: ', self%export
      if (allocated(self%file_out)) write(*, '(2x,a)') 'Grid function will be exported to file: ', self%file_out
    write(*, *)
   end subroutine


! ------------------------------------------------------------------------------
  subroutine display_grid(self,my_rank)
    class(grid_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    integer::i, j, my_rank_local

    if (self%id == 0) return

    if (present(my_rank)) then
      my_rank_local = my_rank
    else
      my_rank_local = 0
    end if

    write(*, *)
    write(*, *) ' VISUAL: grid data on proc' 
    write(*, *) '  id = ', self%id, my_rank_local
    write(*, *) '  ndim = ', self%get_ndim(), my_rank_local 
    write(*, *) '  nr_points = ', (self%nr_points(i), i=1,3), my_rank_local
    write(*, *) '  nr_points_total = ', self%nr_points_total, my_rank_local
    write(*, *) '  origin [a.u.] = ', (self%origin(i), i=1,3), my_rank_local
    write(*, *) '  spacing [a.u.] = ', (self%spacing(i), i=1,3), my_rank_local
    write(*, *) '  extents [a.u.] = ', (self%extents(0:2,i), i=1,3), my_rank_local
    write(*, *) '  radial_length [a.u.] = ', self%radial_length, my_rank_local
    write(*, *) '  radial_nr_steps  = ', self%radial_nr_steps, my_rank_local
    write(*, *) '  radial_nr_points_angular  = ', self%radial_nr_points_angular, my_rank_local
    write(*, *) '  radial_lebedev_l  = ', self%radial_lebedev_l, my_rank_local
    write(*, *) '  status = ', self%status, my_rank_local
    write(*, *) '  input = ', self%input, my_rank_local
    write(*, *) '  typ = ', self%typ, my_rank_local
    write(*, *) '  margin = ', self%margin, my_rank_local
    write(*, *) '  integration_order_2d = ', self%integration_order_2d, my_rank_local
    if (allocated(self%import)) &
      write(*, *) '  import = ', self%import, my_rank_local
    if (allocated(self%export)) &
      write(*, *) '  export = ', self%export, my_rank_local
    if (allocated(self%file_inp)) &
      write(*, *) '  file_inp = ', self%file_inp, my_rank_local
    if (allocated(self%file_out)) &
      write(*, *) '  file_out = ', self%file_out, my_rank_local
    write(*, *) '  origin_known = ', self%origin_known, my_rank_local
    write(*, *) '  extents_known = ', self%extents_known, my_rank_local
    write(*, *) '  nr_points_total_known = ', self%nr_points_total_known, my_rank_local
    write(*, *) '  nr_points_known(1:3) = ', (self%nr_points_known(i), i=1,3), my_rank_local
    write(*, *) '  spacing_known(1:3) = ', (self%spacing_known(i), i=1,3), my_rank_local
    write(*, *)
   end subroutine

! ------------------------------------------------------------------------------
  subroutine display_grid_function(self,my_rank)
    class(grid_function_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    integer:: my_rank_local

    if (self%id == 0) return

    if (present(my_rank)) then
      my_rank_local = my_rank
    else
      my_rank_local = 0
    end if

    write(*, *)
    write(*, *) ' VISUAL: grid function on proc'
    write(*, *) '  id = ', self%id, my_rank_local
    if (allocated(self%name)) write(*, *) 'name = ', self%name, my_rank_local
    write(*, *) '  status = ', self%status, my_rank_local
    write(*, *) '  outpri = ', self%outpri, my_rank_local
    write(*, *) '  input = ', self%input, my_rank_local
    write(*, *) '  wf_order = ', self%wf_order, my_rank_local
    write(*, *) '  rspvec_id = ', self%rspvec_id, my_rank_local
    write(*, *) '  nuc = ', self%nuc, my_rank_local
    write(*, *) '  lao = ', self%lao, my_rank_local
    write(*, *) '  lao_skip = ', self%lao_skip, my_rank_local
    write(*, *) '  lao_comp = ', self%lao_comp, my_rank_local
    write(*, *) '  nics = ', self%nics, my_rank_local
    write(*, *) '  nics_origin = ', self%nics_origin, my_rank_local
    write(*, *) '  tensor_order = ', self%tensor_order, my_rank_local
    write(*, *) '  max_ao_g_order = ', self%max_ao_g_order, my_rank_local
    write(*, *) '  max_ao_m_order = ', self%max_ao_m_order, my_rank_local
    write(*, *) '  modify = ', self%modify, my_rank_local
    write(*, *) '  radial_power = ', self%radial_power, my_rank_local
    write(*, *) '  cartesian_power = ', self%cartesian_power, my_rank_local
    write(*, *) '  up_scale = ', self%up_scale, my_rank_local
    write(*, *) '  down_scale = ', self%down_scale, my_rank_local
    write(*, *) '  bed_wave_vector = ', self%bed_wave_vector, my_rank_local
    write(*, *) '  bed_pol_vector = ', self%bed_pol_vector, my_rank_local
    write(*, *) '  bed_freq = ', self%bed_freq, my_rank_local
    write(*, *) '  bedfix_known = ', self%bedfix_known, my_rank_local
    write(*, *) '  id_assigned_grid =', self%id_assigned_grid, my_rank_local
    write(*, *) '  purpose = ', self%purpose, my_rank_local
    if (allocated(self%rspvec_file)) &
      write(*, *) '  rspvec_file = ', self%rspvec_file, my_rank_local
    if (allocated(self%gf_type)) &
      write(*, *) '  gf_type = ', self%gf_type, my_rank_local
    if (allocated(self%export)) &
      write(*, *) '  export = ', self%export, my_rank_local
    if (allocated(self%file_inp)) &
      write(*, *) '  file_inp = ', self%file_inp, my_rank_local
    if (allocated(self%file_out)) &
      write(*, *) '  file_out = ', self%file_out, my_rank_local
    write(*, *)
  end subroutine

! ------------------------------------------------------------------------------
  subroutine display_grid_points(self,my_rank)
    class(grid_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    integer::i, my_rank_local

    if (self%id == 0) return
    if (.not. self%allocated_points()) return

    if (present(my_rank)) then
      my_rank_local = my_rank
    else
      my_rank_local = 0
    end if

    write(*, *)
    write(*, *) ' GRID points (grid id, proc):', self%id, my_rank_local
    call print_rs_2d(self%points, &
                     self%get_ndim(), self%nr_points_total, my_rank_local)
  end subroutine

! ------------------------------------------------------------------------------
  subroutine display_grid_weights(self,my_rank)
    class(grid_t), intent(in):: self
    integer, intent(in), optional:: my_rank
    integer::i, my_rank_local

    if (self%id == 0) return
    if (.not. self%allocated_weights()) return

    if (present(my_rank)) then
      my_rank_local = my_rank
    else
      my_rank_local = 0
    end if

    write(*, *)
    write(*, *) ' GRID weights (grid id, proc):', self%id, my_rank_local
    call print_rs_2d(self%weights, &
                     1, self%nr_points_total, my_rank_local)
  end subroutine

! ------------------------------------------------------------------------------
  subroutine print_rs_2d(array, n, m, my_rank)
    real(8), intent(in):: array(n, m)
    integer, intent(in):: n, m
    integer, intent(in), optional:: my_rank
    integer :: i, j, my_rank_local

    if (present(my_rank)) then
      my_rank_local = my_rank
    else
      my_rank_local = 0
    end if

    do j = 1, m
      write(*, *) my_rank_local, (array(i, j), i = 1, n)
    end do
  end subroutine

end module
