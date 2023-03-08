!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful, 
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module visual_cfg

   implicit none

   save

   integer,      public, parameter:: max_batch_length = 500
   integer,      public, parameter:: max_para = 100
   integer,      public            :: visual_cfg_nr_dmat = 0
   character(6), public            :: visual_cfg_dmat_file(max_para)
   integer,      public            :: visual_cfg_dmat_file_record(max_para)
   integer,      public            :: visual_cfg_ref_nucleus(max_para)
   integer,      public            :: visual_cfg_property(max_para)
   logical,      public            :: visual_cfg_skip(max_para) = .false.
   integer,      public            :: visual_cfg_property_to_matrix(max_para) = 0
   integer,      public            :: visual_cfg_irep_xvector(max_para)
   logical,      public            :: visual_cfg_out_file_custom = .false.
   logical,      public            :: visual_cfg_out_scalar_file_custom = .false.
   logical,      public            :: visual_cfg_out_vector_file_custom = .false.
   logical,      public            :: visual_cfg_out_2ndtensor_file_custom = .false.
   logical,      public            :: visual_cfg_add_grid_points_grad_ed = .false.
   logical,      public            :: visual_cfg_add_grid_points_grad_ed_scaled = .false.
   character(60), public           :: visual_cfg_out_scalar_file_name
   character(60), public           :: visual_cfg_out_vector_file_name
   character(60), public           :: visual_cfg_out_2ndtensor_file_name
   logical,      public            :: visual_cfg_london_skip_cube = .false.
   logical,      public            :: visual_cfg_london_skip_vtk  = .false.


!
!  options for exporting real-space data to files
!
! WIP
!  character(len = :), public               :: visual_cfg_plot_file_formats(2) = (/'txt', 'cube', 'csv', 'hdf5'/)
!  character(len = 6), public               :: visual_cfg_txt_plot_file_domains(4) = &
!                                              [character(len = 6):: '2D', '3D', 'line', 'radial']
!
   logical,      public                     :: visual_cfg_plot_file_custom = .false.
   character(len = :), allocatable, public  :: visual_cfg_plot_file_format
   character(len = :), allocatable, public  :: visual_cfg_plot_domain
   character(len = :), allocatable, public  :: visual_cfg_plot_property

!  options for exporting to txt files
!  by default, export real-space data to txt files and use legacy file names ("plot.3d.scalar", etc.)
   logical,      public                     :: visual_cfg_txt_plot_file_use_legacy = .true.
   character(len = :), allocatable, public  :: visual_cfg_txt_plot_file_name
   character(len = :), allocatable, public  :: visual_cfg_txt_plot_file_float_format

!  options for exporting to cube files (applies to scalar fields on 3D cube grids)
!  by default, use old export to cube files
   logical,      public                     :: visual_cfg_cube_plot_file_use_legacy = .true.
   character(len = :), allocatable, public  :: visual_cfg_cube_plot_file_name
   character(len = :), allocatable, public  :: visual_cfg_cube_plot_file_float_format



!
!  option for input data files
!
   logical,      public                     :: visual_cfg_inp_data_files_use_legacy= .false.
   integer,      public                     :: visual_cfg_nr_plot_files = 0

   integer, public:: iq_gamma5   =  1
   integer, public:: iq_density  =  2
   integer, public:: iq_elf      =  3
   integer, public:: iq_j        =  4
   integer, public:: iq_rotj     =  5
   integer, public:: iq_divj     =  6
   integer, public:: iq_s        =  7
   integer, public:: iq_rots     =  8
   integer, public:: iq_divs     =  9
   integer, public:: iq_jdia     = 10
   integer, public:: iq_edipx    = 11
   integer, public:: iq_edipy    = 12
   integer, public:: iq_edipz    = 13
   integer, public:: iq_bdipx    = 14
   integer, public:: iq_bdipy    = 15
   integer, public:: iq_bdipz    = 16
   integer, public:: iq_ndipx    = 17
   integer, public:: iq_ndipy    = 18
   integer, public:: iq_ndipz    = 19
   integer, public:: iq_bdipxdia = 20
   integer, public:: iq_bdipydia = 21
   integer, public:: iq_bdipzdia = 22
   integer, public:: iq_ndipxdia = 23
   integer, public:: iq_ndipydia = 24
   integer, public:: iq_ndipzdia = 25
   integer, public:: iq_esp      = 26
   integer, public:: iq_espe     = 27
   integer, public:: iq_espn     = 28
   integer, public:: iq_esprho   = 29
   integer, public:: iq_esperho  = 30
   integer, public:: iq_espnrho  = 31
   integer, public:: iq_kin      = 32
   integer, public:: iq_kin_ls   = 33
   integer, public:: iq_kin_sl   = 34
   integer, public:: iq_kin_tau  = 35
   integer, public:: iq_kin_lap  = 36
   integer, public:: iq_kin_nr   = 37
   integer, public:: iq_jx       = 38
   integer, public:: iq_jy       = 39
   integer, public:: iq_jz       = 40
   integer, public:: iq_tcos     = 41
   integer, public:: iq_tsin     = 42   
   integer, public:: iq_density_gradient     =  43
   integer, public:: iq_density_hessian      =  44
   integer, public:: iq_density_laplacian    =  45
   integer, public:: iq_density_sign_lambda2 =  46
   integer, public:: iq_rdg      = 47
   integer, public:: iq_gradj    = 48
   

   logical, public:: visual_cfg_print_data_to_dirac_output = .false.

!                     visualization in 2d
   logical, public:: visual_cfg_2d             = .false.
   real(8), public:: visual_cfg_2d_p_origin(3) = 0.0d0
   real(8), public:: visual_cfg_2d_p_right(3)  = 0.0d0
   real(8), public:: visual_cfg_2d_p_top(3)    = 0.0d0
   integer, public:: visual_cfg_2d_nr_right    = 0
   integer, public:: visual_cfg_2d_nr_top      = 0
   character(len = 6), public:: visual_cfg_2d_data_types(4) = &
                       [character(len = 6):: 'scalar', 'vector', 'tensor', 'flow']

!                     integration in 2d
   logical, public:: visual_cfg_2d_integration             = .false.
   real(8), public:: visual_cfg_2d_integration_p_origin(3) = 0.0d0
   real(8), public:: visual_cfg_2d_integration_p_right(3)  = 0.0d0
   real(8), public:: visual_cfg_2d_integration_p_top(3)    = 0.0d0
   integer, public:: visual_cfg_2d_integration_nr_right    = 0
   integer, public:: visual_cfg_2d_integration_nr_top      = 0
   integer, public:: visual_cfg_2d_integration_order       = 0

!                     visualization in 3d
   logical, public:: visual_cfg_3d         = .false.
   logical, public:: visual_cfg_3d_fast    = .false.
   integer, public:: visual_cfg_ncube(1:3) = 0
   real(8), public:: visual_cfg_add_3d     = 4.0d0
   logical, public:: visual_cfg_3d_gridimp = .false.
   logical, public:: visual_cfg_3d_readjb  = .false.
   character(80), public:: visual_cfg_3d_gridfil
   character(80), public:: visual_cfg_3d_grid_units
   character(80), public:: visual_cfg_3d_jbfile
   character(len = 10), public:: visual_cfg_3d_data_types(5) = &
                       [character(len = 10)::'scalar', 'vector', 'tensor', 'cube', 'vector_vtk']

!                     integration in 3d
   logical, public:: visual_cfg_3d_integration = .false.

!                     visualization along line
   logical, public:: visual_cfg_line          = .false.
   real(8), public:: visual_cfg_line_from(3)  = 0.0d0
   real(8), public:: visual_cfg_line_to(3)    = 0.0d0
   integer, public:: visual_cfg_line_nr_steps = 0
   character(len = 6), public:: visual_cfg_line_data_types(3) = &
                       [character(len = 6)::'scalar', 'vector', 'tensor']

!                                  densities for a list of points
   logical, public              :: visual_cfg_list              = .false.
   integer, public              :: visual_cfg_nr_points_in_list = 0
   real(8), public, allocatable:: visual_cfg_xyz_list(:, :)
   character(len = 6), public:: visual_cfg_list_data_types(3) = &
                       [character(len = 6)::'scalar', 'vector', 'tensor']


!                     radial distributions
   logical, public:: visual_cfg_radial          = .false.
   real(8), public:: visual_cfg_radial_from(3)  = 0.0d0
   integer, public:: visual_cfg_radial_nr_steps = 0
   real(8), public:: visual_cfg_radial_length   = 0
   character(len = 6), public:: visual_cfg_radial_data_types(3) = &
                       [character(len = 6)::'scalar', 'vector', 'tensor']
   
   real(8),       public:: visual_cfg_scale              = 1.0d0
   real(8),       public:: visual_cfg_cartesian_power(3) = 0.0   
   real(8),       public:: visual_cfg_gauge_origin(3)    = 0.0d0
   real(8),       public:: visual_cfg_wave_vector(3)     = 0.0d0
   real(8),       public:: visual_cfg_freq               = 0.0d0   
   real(8),       public:: visual_cfg_pol_vector(3)      = 0.0d0   
   logical,       public:: visual_cfg_force_small_ao     = .false.

   logical,       public:: visual_cfg_nics              = .false.
   real(8),       public:: visual_cfg_nics_origin(3)    = 0.0d0

   logical,      public:: visual_cfg_london               = .false.
   logical,      public:: visual_cfg_london_skip_kappa    = .false.
   logical,      public:: visual_cfg_london_skip_direct   = .false.
   logical,      public:: visual_cfg_london_skip_ro       = .false.
   logical,      public:: visual_cfg_london_none          = .false.
   character(1), public:: visual_cfg_london_component

   integer,      public:: visual_cfg_irep_conmat = 0
   integer,      public:: visual_cfg_itim_conmat = 0

   logical, public              :: visual_cfg_use_orbital_string = .false.
   real(8), public, allocatable:: visual_cfg_occupation(:, :)

   logical,      public:: visual_cfg_tensor          = .false.

   private

end module
