module visual_grid_functions


   use visual_io
   use visual_cfg
   use visual_utils
   use visual_datatypes

   implicit none

   public get_grid_function_info
   public get_grid_function_values

contains

  subroutine get_grid_function_info(grid_function_id, grid_function)
    integer, intent(in):: grid_function_id !fixme
    type(grid_function_t), intent(inout):: grid_function

    if (grid_function%input  ==  1) then
      !call get_grid_function_info_from_file(grid_function)
      !if (grid_function%id == 0) grid_function%id = grid_function_id
      print *, 'INFO from get_grid_function_info: grid_function%id = ', grid_function%id
    !else
    !  call get_grid_function_essentials(grid_function)
    end if

    if (grid_function%id == 0) grid_function%id = grid_function_id

  end subroutine get_grid_function_info


  subroutine get_grid_function_values(grid_function)
    type(grid_function_t), intent(inout):: grid_function

    if (.not. grid_function%allocated_values()) return

    ! gosia todo
    !if (grid_function%input == 1) then
    !  call get_grid_function_values_from_file(grid_function)
    !  write(*, *) 'grid_function values for grid_function id = ', grid_function%id, ' are read from file ', grid_function%file_inp
    !else if (grid_function%input == 0) then
    !  call calculate_grid_function_values(grid_function)
    !  write(*, *) 'grid_function values for grid_function id = ', grid_function%id, ' are calculated'
    !end if
  end subroutine

end module
