module test_derv_m

  use kinds_m
  use vwn_m
  use vwn_derv_m

  implicit none

contains

  subroutine test_vwn_derv()

    real(kind_dp) :: n, s, fn, fs, w
    real(kind_dp) :: d(0:60), d_p(0:60), d_m(0:60)

    n  = 1.0d0
    fn = 1.0d-5
    s  = 0.0d0
    fs = 0.0d0
    w  = 1.0d0

    call vwn_derv(0.5d0*(n + fn + s + fs), 0.5d0*(n + fn - s - fs), w, d)
    call dcopy(61, d, 1, d_p, 1)

    call vwn_derv(0.5d0*(n - fn + s - fs), 0.5d0*(n - fn - s + fs), w, d)
    call dcopy(61, d, 1, d_m, 1)

    call vwn_derv(0.5d0*n, 0.5d0*n, w, d)

    d_p(0) = vwn_energy(0.5d0*(n + fn + s + fs), 0.5d0*(n + fn - s - fs), w)
    d_m(0) = vwn_energy(0.5d0*(n - fn + s - fs), 0.5d0*(n - fn - s + fs), w)

    write(0,*) 'd(10)', (d_p( 0) - d_m( 0))/(2.0d0*fn), d(10)
    write(0,*) 'd(20)', (d_p(10) - d_m(10))/(2.0d0*fn), d(20)
    write(0,*) 'd(30)', (d_p(20) - d_m(20))/(2.0d0*fn), d(30)
    write(0,*) 'd(40)', (d_p(30) - d_m(30))/(2.0d0*fn), d(40)
    write(0,*) 'd(50)', (d_p(40) - d_m(40))/(2.0d0*fn), d(50)
    write(0,*) 'd(60)', (d_p(50) - d_m(50))/(2.0d0*fn), d(60)

  end subroutine

end module
