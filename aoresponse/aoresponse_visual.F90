module aoresponse_visual_m

  use kinds_m
  use density_m
  use grid_m
! use memory_allocator
  use matrix_oper
  use matrix_mod

  implicit none

contains

  subroutine aoresponse_visual_driver(D_0, D)

    type(matrix)     :: D_0, D(:)
!   ----------------------------------------------------------------------------
    type(grid_t)     :: grid
    type(gao_data_t) :: gao_data
    type(density_t)  :: n(100)
    real(kind_dp)    :: integral(100)
    integer          :: i, ipoint

    do i = 1, size(D)
      call init_density(n(i), '1')
    end do

    call init_dft_grid(D_0, grid)
    call init_gao_data(gao_data, 0)

    integral = 0.0d0

    do ipoint = 1, grid%nr_points

      call get_gao_data(gao_data,       &
                        grid%x(ipoint), &
                        grid%y(ipoint), &
                        grid%z(ipoint))

      do i = 1, size(D)
        call get_or_put_density(n(i), D(i), .true., gao_data)
        integral(i) = integral(i) - n(i)%scalar*grid%w(ipoint)*grid%z(ipoint)
      end do

    end do

    do i = 1, size(D)
      write(*,*) 'blub', i, integral(i)*4.0d0
    end do

    call del_grid(grid)
    call del_gao_data(gao_data)

  end subroutine

end module
