module visual_dirac_input

   use visual_cfg

   implicit none

   public verify_input_keywords

   private

   integer, parameter:: lupri = 6

   contains

   subroutine verify_input_keywords()

      call write_visual_cfg()

   end subroutine verify_input_keywords


   subroutine write_visual_cfg()
 
      write(lupri, *)    '--------------------------------------------------------------------------------------'
 
      write(lupri, '(/,2x, a)') '*** Visualisation module ***'

      if (visual_cfg_3d) then
         write(lupri, *) 'Visualization in 3d done with the following settings:'
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

      write(lupri, *)    '--------------------------------------------------------------------------------------'
 
   end subroutine


end module
