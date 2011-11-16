subroutine integrate_xc_potential_real(mat_dim, &
                                       dmat_0,  &
                                       fmat)

   use xcint_main

   implicit none

   integer, intent(in) :: mat_dim
   real(8), intent(in) :: dmat_0(*)
   real(8), intent(in) :: fmat(*)

   call integrate_xc(xc_mat_dim      = mat_dim, &
                     xc_nz           = 1,       &
                     xc_dmat_0       = dmat_0,  &
                     xc_nr_dmat      = 0,       &
                     xc_nr_fmat      = 1,       &
                     xc_fmat         = fmat,    &
                     xc_do_potential = .true.)

end subroutine
