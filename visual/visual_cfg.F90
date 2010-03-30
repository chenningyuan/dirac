module visual_cfg

  implicit none

  save

  integer, parameter, public :: max_nr_xvectors       = 6
  integer, parameter, public :: max_nr_points_in_list = 100

  character(6),  public :: visual_cfg_xvector_file(max_nr_xvectors)
  character(72), public :: visual_cfg_orbital_string(2)
  character(6),  public :: visual_cfg_xvector_file(max_nr_xvectors)

  logical,       public :: visual_cfg_list               = .false.
  logical,       public :: visual_cfg_line               = .false.
  logical,       public :: visual_cfg_2d                 = .false.
  logical,       public :: visual_cfg_3d                 = .false.
  logical,       public :: visual_cfg_use_orbital_string = .false.
  logical,       public :: visual_cfg_integrate          = .false.
  logical,       public :: visual_cfg_radial             = .false.
  logical,       public :: visual_cfg_force_small_ao     = .false.
  logical,       public :: visual_cfg_london             = .false.

  integer,       public :: visual_cfg_ncube(0:3)                       = 80
  integer,       public :: visual_cfg_nr_points_in_list                = 0
  integer,       public :: visual_cfg_nr_xvectors                      = 0
  integer,       public :: visual_cfg_iq_xvector(max_nr_xvectors)      = 0
  integer,       public :: visual_cfg_irep_xvector(max_nr_xvectors)    = 0
  integer,       public :: visual_cfg_iwhich_xvector(max_nr_xvectors)  = 1
  integer,       public :: visual_cfg_ioption_xvector(max_nr_xvectors) = 0
  integer,       public :: visual_cfg_line_nr_steps                    = 0
  integer,       public :: visual_cfg_orientation                      = 1
  integer,       public :: visual_cfg_connection                       = 0

  real(8),       public :: visual_cfg_add_3d                             = 4.0d0
  real(8),       public :: visual_cfg_common_factor                      = 1.0d0
  real(8),       public :: visual_cfg_xyz_list(max_nr_points_in_list, 3) = 0.0d0
  real(8),       public :: visual_cfg_origin_2d(3)                       = 0.0d0
  real(8),       public :: visual_cfg_gauge_origin(3)                    = 0.0d0
  real(8),       public :: visual_cfg_add_2d                             = 4.0d0
  real(8),       public :: visual_cfg_step_2d                            = 0.1d0
  real(8),       public :: visual_cfg_line_from(3)                       = 0.0d0
  real(8),       public :: visual_cfg_line_to(3)                         = 0.0d0

  private

end module
