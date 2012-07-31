module openrsp_visual

  use kinds_m
  use scf_density_m
! use grid_m
! use memory_allocator
! use matrix_oper
  use matrix_defop_old

  implicit none

  public openrsp_visual_driver

  private

  type cube_file_t
    integer      :: matrix_pointer
    integer      :: file_unit
    character(9) :: file_name
    integer      :: op
  end type

  real(kind_dp), parameter :: small = 1.0d-10

contains

  subroutine init_cube_file(self, matrix_pointer, file_unit, file_name, op)

    type(cube_file_t) :: self
!   ----------------------------------------------------------------------------
    integer           :: matrix_pointer, file_unit, op
    character(*)      :: file_name

    self%matrix_pointer = matrix_pointer
    self%file_unit      = file_unit
    self%file_name      = file_name
    self%op             = op

    file_unit = file_unit + 1

  end subroutine

  subroutine openrsp_visual_driver(D)

    type(matrix)                   :: D(:)
!   ----------------------------------------------------------------------------
    type(grid_t)                   :: dft_grid
    type(cube_t)                   :: cube
    type(gao_data_t)               :: gao_data
    type(density_t)                :: n
    real(kind_dp),     allocatable :: integral(:)
    real(kind_dp),     allocatable :: cube_slice(:, :)
    type(cube_file_t), allocatable :: cube_file(:)
    real(kind_dp)                  :: r(3), s
    integer                        :: ip, ic, id, ix, ia, iz, iq
    integer                        :: k, l, m
    integer,           parameter   :: ncube = 80
    integer                        :: file_unit = 113

#include "dgroup.h"


!   scale matrices
!   ==============

    do id = 1, size(D)
      D(id) = 4.0d0*D(id)
    end do


!   half-phase
!   ==========

    if (nz < 4) then
      do iz = 1, nz
        do id = 1, size(D)
          iq = ipqtoq(iz, D(id)%irep)
          call q2bphase('D', iq, 1, D(id)%elms((iz-1)*D(id)%nrow*D(id)%ncol + 1))
        end do
      end do
    end if


    call init_density(n, '1')

    call init_gao_data(gao_data, 0)

    allocate(cube_file(size(D) + 3*(size(D) - 1)))
    allocate(cube_slice(size(cube_file), ncube))
    allocate(integral(size(cube_file)))

    if (size(D) == 4) then
      call init_cube_file(cube_file( 1),  1, file_unit, '3d_n_____', 0)
      call init_cube_file(cube_file( 2),  2, file_unit, '3d_n_x___', 0)
      call init_cube_file(cube_file( 3),  3, file_unit, '3d_n_y___', 0)
      call init_cube_file(cube_file( 4),  4, file_unit, '3d_n_z___', 0)
      call init_cube_file(cube_file( 5),  2, file_unit, '3d_n_x__x', 1)
      call init_cube_file(cube_file( 6),  3, file_unit, '3d_n_y__x', 1)
      call init_cube_file(cube_file( 7),  4, file_unit, '3d_n_z__x', 1)
      call init_cube_file(cube_file( 8),  2, file_unit, '3d_n_x__y', 2)
      call init_cube_file(cube_file( 9),  3, file_unit, '3d_n_y__y', 2)
      call init_cube_file(cube_file(10),  4, file_unit, '3d_n_z__y', 2)
      call init_cube_file(cube_file(11),  2, file_unit, '3d_n_x__z', 3)
      call init_cube_file(cube_file(12),  3, file_unit, '3d_n_y__z', 3)
      call init_cube_file(cube_file(13),  4, file_unit, '3d_n_z__z', 3)
    elseif (size(D) == 13) then
      call init_cube_file(cube_file( 1),  1, file_unit, '3d_n_____', 0)
      call init_cube_file(cube_file( 2),  2, file_unit, '3d_n_x___', 0)
      call init_cube_file(cube_file( 3),  3, file_unit, '3d_n_y___', 0)
      call init_cube_file(cube_file( 4),  4, file_unit, '3d_n_z___', 0)
      call init_cube_file(cube_file( 5),  5, file_unit, '3d_n_xx__', 0)
      call init_cube_file(cube_file( 6),  6, file_unit, '3d_n_xy__', 0)
      call init_cube_file(cube_file( 7),  7, file_unit, '3d_n_xz__', 0)
      call init_cube_file(cube_file( 8),  8, file_unit, '3d_n_yx__', 0)
      call init_cube_file(cube_file( 9),  9, file_unit, '3d_n_yy__', 0)
      call init_cube_file(cube_file(10), 10, file_unit, '3d_n_yz__', 0)
      call init_cube_file(cube_file(11), 11, file_unit, '3d_n_zx__', 0)
      call init_cube_file(cube_file(12), 12, file_unit, '3d_n_zy__', 0)
      call init_cube_file(cube_file(13), 13, file_unit, '3d_n_zz__', 0)
      call init_cube_file(cube_file(14),  2, file_unit, '3d_n_x__x', 1)
      call init_cube_file(cube_file(15),  3, file_unit, '3d_n_y__x', 1)
      call init_cube_file(cube_file(16),  4, file_unit, '3d_n_z__x', 1)
      call init_cube_file(cube_file(17),  5, file_unit, '3d_n_xx_x', 1)
      call init_cube_file(cube_file(18),  6, file_unit, '3d_n_xy_x', 1)
      call init_cube_file(cube_file(19),  7, file_unit, '3d_n_xz_x', 1)
      call init_cube_file(cube_file(20),  8, file_unit, '3d_n_yx_x', 1)
      call init_cube_file(cube_file(21),  9, file_unit, '3d_n_yy_x', 1)
      call init_cube_file(cube_file(22), 10, file_unit, '3d_n_yz_x', 1)
      call init_cube_file(cube_file(23), 11, file_unit, '3d_n_zx_x', 1)
      call init_cube_file(cube_file(24), 12, file_unit, '3d_n_zy_x', 1)
      call init_cube_file(cube_file(25), 13, file_unit, '3d_n_zz_x', 1)
      call init_cube_file(cube_file(26),  2, file_unit, '3d_n_x__y', 2)
      call init_cube_file(cube_file(27),  3, file_unit, '3d_n_y__y', 2)
      call init_cube_file(cube_file(28),  4, file_unit, '3d_n_z__y', 2)
      call init_cube_file(cube_file(29),  5, file_unit, '3d_n_xx_y', 2)
      call init_cube_file(cube_file(30),  6, file_unit, '3d_n_xy_y', 2)
      call init_cube_file(cube_file(31),  7, file_unit, '3d_n_xz_y', 2)
      call init_cube_file(cube_file(32),  8, file_unit, '3d_n_yx_y', 2)
      call init_cube_file(cube_file(33),  9, file_unit, '3d_n_yy_y', 2)
      call init_cube_file(cube_file(34), 10, file_unit, '3d_n_yz_y', 2)
      call init_cube_file(cube_file(35), 11, file_unit, '3d_n_zx_y', 2)
      call init_cube_file(cube_file(36), 12, file_unit, '3d_n_zy_y', 2)
      call init_cube_file(cube_file(37), 13, file_unit, '3d_n_zz_y', 2)
      call init_cube_file(cube_file(38),  2, file_unit, '3d_n_x__z', 3)
      call init_cube_file(cube_file(39),  3, file_unit, '3d_n_y__z', 3)
      call init_cube_file(cube_file(40),  4, file_unit, '3d_n_z__z', 3)
      call init_cube_file(cube_file(41),  5, file_unit, '3d_n_xx_z', 3)
      call init_cube_file(cube_file(42),  6, file_unit, '3d_n_xy_z', 3)
      call init_cube_file(cube_file(43),  7, file_unit, '3d_n_xz_z', 3)
      call init_cube_file(cube_file(44),  8, file_unit, '3d_n_yx_z', 3)
      call init_cube_file(cube_file(45),  9, file_unit, '3d_n_yy_z', 3)
      call init_cube_file(cube_file(46), 10, file_unit, '3d_n_yz_z', 3)
      call init_cube_file(cube_file(47), 11, file_unit, '3d_n_zx_z', 3)
      call init_cube_file(cube_file(48), 12, file_unit, '3d_n_zy_z', 3)
      call init_cube_file(cube_file(49), 13, file_unit, '3d_n_zz_z', 3)
    else
      call quit('error in openrsp_visual')
    end if

    call init_dft_grid(D(1), dft_grid)

    integral = 0.0d0

    do ip = 1, dft_grid%nr_points

      call get_gao_data(gao_data,       &
                        dft_grid%x(ip), &
                        dft_grid%y(ip), &
                        dft_grid%z(ip))

      do ic = 1, size(cube_file)

        call get_or_put_density(n, D(cube_file(ic)%matrix_pointer), .true., gao_data)
        s = n%scalar*dft_grid%w(ip)

        select case (cube_file(ic)%op)

          case (1)
            s = s*dft_grid%x(ip)
          case (2)
            s = s*dft_grid%y(ip)
          case (3)
            s = s*dft_grid%z(ip)

        end select

        integral(ic) = integral(ic) + s

      end do

    end do

    do ic = 1, size(cube_file)
      if (abs(integral(ic)) > small) then
        write(*,*) 'integral', cube_file(ic)%file_name, integral(ic)
      end if
    end do

    call del_grid(dft_grid)

    call init_cube(cube, ncube)

    do ic = 1, size(cube_file)

      open(cube_file(ic)%file_unit,          &
           file   = cube_file(ic)%file_name, &
           form   = 'formatted',             &
           access = 'sequential',            &
           status = 'new')

      write(cube_file(ic)%file_unit, '(a)') 'some'
      write(cube_file(ic)%file_unit, '(a)') 'text'

      write(cube_file(ic)%file_unit, '(i5, 3f12.6)') cube%nr_atoms,       &
                                                     cube%initial_xyz(1), &
                                                     cube%initial_xyz(2), &
                                                     cube%initial_xyz(3)

      do ix = 1, 3
        write(cube_file(ic)%file_unit, '(i5, 3f12.6)') cube%nr_points,          &
                                                       cube%step_vector(ix, 1), &
                                                       cube%step_vector(ix, 2), &
                                                       cube%step_vector(ix, 3)
      end do

      do ia = 1, cube%nr_atoms
        write(cube_file(ic)%file_unit, '(2i5, a7, 3f12.6)') cube%atom_charge(ia), &
                                                            cube%atom_charge(ia), &
                                                            '.000000',            &
                                                            cube%atom_xyz(ia, 1), &
                                                            cube%atom_xyz(ia, 2), &
                                                            cube%atom_xyz(ia, 3)
      end do

    end do

    do k = 1, cube%nr_points
      do l = 1, cube%nr_points

        cube_slice = 0.0d0
        do m = 1, cube%nr_points

           do ix = 1, 3
             r(ix) = cube%initial_xyz(ix) + (k - 1)*cube%step_vector(ix, 1) &
                                          + (l - 1)*cube%step_vector(ix, 2) &
                                          + (m - 1)*cube%step_vector(ix, 3)
           end do

           call get_gao_data(gao_data, r(1), r(2), r(3))

           do ic = 1, size(cube_file)

             call get_or_put_density(n, D(cube_file(ic)%matrix_pointer), .true., gao_data)
             cube_slice(ic, m) = n%scalar

             select case (cube_file(ic)%op)

               case (1)
                 cube_slice(ic, m) = cube_slice(ic, m)*r(1)
               case (2)
                 cube_slice(ic, m) = cube_slice(ic, m)*r(2)
               case (3)
                 cube_slice(ic, m) = cube_slice(ic, m)*r(3)

             end select

           end do

        end do

        do ic = 1, size(cube_file)
          write(cube_file(ic)%file_unit, '(6e13.5)') (cube_slice(ic, m), m = 1, cube%nr_points)
        end do

      end do
    end do

    do ic = 1, size(cube_file)
      close(cube_file(ic)%file_unit, status = 'keep')
    end do

    call del_cube(cube)
    call del_gao_data(gao_data)

    deallocate(integral)
    deallocate(cube_slice)
    deallocate(cube_file)


!   scale matrices
!   ==============

    do id = 1, size(D)
      D(id) = 0.25d0*D(id)
    end do


!   half-phase
!   ==========

    if (nz < 4) then
      do iz = 1, nz
        do id = 1, size(D)
          iq = ipqtoq(iz, D(id)%irep)
          call q2bphase('D', iq, 1, D(id)%elms((iz-1)*D(id)%nrow*D(id)%ncol + 1))
        end do
      end do
    end if

  end subroutine

end module
