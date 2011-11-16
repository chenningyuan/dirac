subroutine xcint_potential_rks(mat_dim, &
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

subroutine xcint_launch_slave_process()

  use xcint_main

  implicit none

  call integrate_xc(xc_mat_dim = 1,         &
                    xc_nz      = 1,         &
                    xc_dmat_0  = (/0.0d0/), &
                    xc_nr_dmat = 0,         &
                    xc_nr_fmat = 0)

end subroutine

function xcint_get_xc_energy()

  use xcint_main

  implicit none

  real(8) :: xcint_get_xc_energy

  xcint_get_xc_energy = get_xc_energy()

end function

function xcint_get_xc_mat_energy()

  use xcint_main

  implicit none

  real(8) :: xcint_get_xc_mat_energy

  xcint_get_xc_mat_energy = get_xc_mat_energy()

end function
