program main

   use interface_ao
   use interface_file_io
   use xcint_main
   use xc_derv
   use fake_input

   implicit none

   real(8), allocatable :: dmat(:, :)
   real(8), allocatable :: fmat(:, :, :, :)
   real(8), allocatable :: fmat_reference(:, :)
   integer, parameter   :: unit_grid = 44

   call read_input()
   call interface_ao_read(.false.)

   fun_is_automatic = .true.

   allocate(dmat(nr_ao, nr_ao))
   dmat = 0.0d0

   allocate(fmat(nr_ao, nr_ao, nr_quaternion_blocks, 1))
   fmat = 0.0d0

   call matrix_from_file('xc_dmat0', nr_ao, dmat)

   open(unit_grid,                 &
        file   = 'numerical_grid', &
        status = 'unknown',        &
        form   = 'formatted',      &
        access = 'sequential')
   rewind(unit_grid)

   call xcint_potential_rks(nr_ao, &
                            dmat,  &
                            fmat)
!  call integrate_xc(                                        &
!                     xc_mat_dim      = nr_ao,                &
!                     xc_nz           = nr_quaternion_blocks, &
!                     xc_dmat_0       = dmat,                 &
!                     xc_nr_dmat      = 0,                    &
!                     xc_nr_fmat      = 1,                    &
!                     xc_fmat         = fmat,                 &
!                     xc_do_potential = .true.                &
!                    )

   deallocate(dmat)

   allocate(fmat_reference(nr_ao, nr_ao))
   call matrix_from_file('xc_fmat', nr_ao, fmat_reference)
   call compare_matrices(nr_ao, fmat, fmat_reference)

   deallocate(fmat)
   deallocate(fmat_reference)

contains

   subroutine compare_matrices(nr_ao, mat, mat_reference)

!     --------------------------------------------------------------------------
      integer, intent(in) :: nr_ao
      real(8), intent(in) :: mat(*)
      real(8), intent(in) :: mat_reference(*)
!     --------------------------------------------------------------------------
      integer             :: i
      real(8)             :: diff
      real(8)             :: max_diff
!     --------------------------------------------------------------------------

      max_diff = 0.0d0
      do i = 1, nr_ao*nr_ao
         diff = abs(mat(i) - mat_reference(i))
         if (diff > max_diff) then
            max_diff = diff
         end if
      end do
      print *, 'compare matrices: largest deviation =', max_diff

   end subroutine

end program
