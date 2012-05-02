program main

   use interface_ao
   use interface_file_io

   implicit none

   real(8), allocatable :: dmat(:, :)
   real(8), allocatable :: fmat(:, :, :, :)
   real(8), allocatable :: fmat_reference(:, :)
   integer, parameter   :: unit_grid = 44

   character(80) :: line

   line = 'lda'
   call xcint_parse_functional(line)

   call interface_ao_read(.false.)

   allocate(dmat(nr_ao, nr_ao))
   dmat = 0.0d0

   allocate(fmat(nr_ao, nr_ao, nr_quaternion_blocks, 1))
   fmat = 0.0d0

   call matrix_from_file('xc_dmat0', nr_ao, dmat)
   dmat = 2.0d0*dmat

   open(unit_grid,                 &
        file   = 'numerical_grid', &
        status = 'unknown',        &
        form   = 'formatted',      &
        access = 'sequential')
   rewind(unit_grid)

   call xcint_potential_rks(nr_ao, &
                            dmat,  &
                            fmat)

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
