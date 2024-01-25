module visual_dirac_input

   use visual_cfg

   implicit none

   public verify_input_keywords

   private

   integer, parameter:: lupri = 6

   contains

   subroutine verify_input_keywords()

      if (visual_cfg_3d .and. visual_cfg_3d_gridimp) then
          call quit('Mutually exclusive keywords found in the input: .3D and .3D_IMP, please choose one')
      end if
      if (visual_cfg_3d_readjb .and. visual_cfg_3d_gridimp) then
          call quit('Mutually exclusive keywords found in the input: .READJB and .3D_IMP, please choose one.')
      end if

      if (visual_cfg_3d_gridimp) then
          if (visual_cfg_3d_grid_units .eq. 'angstrom') then
              write(LUPRI, *) 'imported grid is in Angstrom, it will be transformed to a.u.'
          else if (visual_cfg_3d_grid_units .eq. 'au') then
              write(LUPRI, *) 'imported grid is in atomic units'
          else
              call quit('unrecognized unit in 3D_IMP keyword; choose "angstrom" or "au"')
          end if
      end if

      if (.not. visual_cfg_txt_plot_file_use_legacy &
        .and..not.visual_cfg_cube_plot_file_use_legacy) then
         visual_cfg_plot_file_custom = .true.
      end if

      call write_visual_cfg()

   end subroutine verify_input_keywords


   subroutine write_visual_cfg()
 
      write(lupri, *)    '--------------------------------------------------------------------------------------'
 
      write(lupri, '(/,2x, a)') '*** Visualisation module ***'

      if (visual_cfg_3d) then
         write(lupri, *) '***** Visualization in 3d done with the following settings *****'
         write(lupri, *) 'Number of points in x, y, z direction: ', visual_cfg_ncube(1), &
                                                                    visual_cfg_ncube(2), &
                                                                    visual_cfg_ncube(3) 
         write(lupri, *) 'Space around the cube file is ', visual_cfg_add_3d
         write(lupri, *) 'Gauge origin is ', visual_cfg_gauge_origin(1), &
                                             visual_cfg_gauge_origin(2), &
                                             visual_cfg_gauge_origin(3)
      end if

      if (visual_cfg_3d_fast) then 
         write(lupri, *) 'Visualization in 3d with fast evaluation of the molecular electrostatic potential', &
                         ' done with the following settings:'
         write(lupri, *) 'Number of points in x, y, z direction: ', visual_cfg_ncube(1), &
                                                                    visual_cfg_ncube(2), &
                                                                    visual_cfg_ncube(3) 
         write(lupri, *) 'Space around the cube file is ', visual_cfg_add_3d
      end if
 
      if (visual_cfg_3d_gridimp) then
         write(lupri, *) 'Visualization in 3d done on the grid imported from file ', visual_cfg_3d_gridfil
      end if
 
      if (visual_cfg_3d_readjb) then
         write(lupri, *) 'Visualization in 3d done on the grid imported from file ', visual_cfg_3d_jbfile
         write(lupri, *) 'And the magnetically-induced current density written on that file is used'
      end if

      write(lupri, *) '***** Grids used for data export *****'
      write(lupri, *) 'Number of grid files in use: ', visual_cfg_nr_grid_files

      write(lupri, *) '***** Visualization data exported to files *****'
      if (visual_cfg_plot_file_format .eq. 'txt') then
         if (visual_cfg_txt_plot_file_use_legacy) then
            write(lupri, *) 'using legacy files'
         else
            write(lupri, *) 'file name: ', visual_cfg_txt_plot_file_name
            write(lupri, *) 'domain:    ', visual_cfg_plot_domain
            write(lupri, *) 'property:  ', visual_cfg_plot_property
         end if
      end if

      write(lupri, *)    '--------------------------------------------------------------------------------------'
 
   end subroutine


end module
