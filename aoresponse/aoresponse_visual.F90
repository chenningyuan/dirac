module aoresponse_visual_m

  use kinds_m
  use scf_density_m
  use grid_m
! use memory_allocator
  use matrix_oper
  use matrix_mod

  implicit none

contains

  subroutine aoresponse_visual_driver(D_0, D)

    type(matrix)     :: D_0, D(:)
!   ----------------------------------------------------------------------------
    type(grid_t)     :: dft_grid
    type(cube_t)     :: cube
    type(gao_data_t) :: gao_data
    type(density_t)  :: n(100)
    real(kind_dp)    :: integral(100)
    integer          :: i, ipoint
    integer          :: file_unit
    character(6)     :: file_name

    call init_gao_data(gao_data, 0)

    do i = 1, size(D)
      call init_density(n(i), '1')
    end do

    call init_dft_grid(D_0, dft_grid)

    integral = 0.0d0

    do ipoint = 1, dft_grid%nr_points

      call get_gao_data(gao_data,           &
                        dft_grid%x(ipoint), &
                        dft_grid%y(ipoint), &
                        dft_grid%z(ipoint))

      do i = 1, size(D)
        call get_or_put_density(n(i), D(i), .true., gao_data)
        integral(i) = integral(i) - n(i)%scalar*dft_grid%w(ipoint)*dft_grid%z(ipoint)
      end do

    end do

    do i = 1, size(D)
      write(*,*) 'blub', i, integral(i)*4.0d0
    end do

    call del_grid(dft_grid)








    call init_cube(cube, 80)

    file_unit = 13
    file_name = 'PLOTbl'

    open(file_unit,             &
         file   = file_name,    &
         form   = 'formatted',  &
         access = 'sequential', &
         status = 'new')

    write(file_unit, '(a)') 'some'
    write(file_unit, '(a)') 'text'

    write(file_unit, '(i5, 3f12.6)') cube%nr_atoms,       &
                                     cube%initial_xyz(1), &
                                     cube%initial_xyz(2), &
                                     cube%initial_xyz(3)

    do i = 1, 3
      write(file_unit, '(i5, 3f12.6)') cube%nr_points,         &
                                       cube%step_vector(i, 1), &
                                       cube%step_vector(i, 2), &
                                       cube%step_vector(i, 3)
    end do

    do i = 1, cube%nr_atoms
      write(file_unit, '(2i5, a7, 3f12.6)') cube%atom_charge(i), &
                                            cube%atom_charge(i), &
                                            '.000000',           &
                                            cube%atom_xyz(i, 1), &
                                            cube%atom_xyz(i, 2), &
                                            cube%atom_xyz(i, 3)
    end do


!       DO I = 1,NUCDEP
!         WRITE(LUNIT_3D,'(I5,4F12.6)')
!    &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
!       ENDDO


!       calculate cube parameters

!       NCUBE(0)      number of atoms
!       NCUBE(1)      number of points in the "x"-direction
!       NCUBE(2)      number of points in the "y"-direction
!       NCUBE(3)      number of points in the "z"-direction

!       STEP_3D(*,0)  coordinates of initial point
!       STEP_3D(*,1)  step vector for the "x"-direction
!       STEP_3D(*,2)  step vector for the "y"-direction
!       STEP_3D(*,3)  step vector for the "z"-direction

!       ATOM_XYZ(1,I) x-coordinate of atom I
!       ATOM_XYZ(2,I) y-coordinate of atom I
!       ATOM_XYZ(3,I) z-coordinate of atom I
!       ATOM_XYZ(4,I) charge of atom I





    close(file_unit, status = 'keep')


















    call del_cube(cube)








    call del_gao_data(gao_data)

  end subroutine

end module
