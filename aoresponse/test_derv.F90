module test_derv_m

  use kinds_m
  use vwn_m
  use vwn_derv_m

  implicit none

contains

  subroutine test_vwn_derv()

    real(kind_dp) :: n, s, fn, fs, w
    real(kind_dp) :: d(10)
    real(kind_dp) :: d_0,   d_r,   d_rr,   d_rrr,   d_rrrr
    real(kind_dp) :: d_0_p, d_r_p, d_rr_p, d_rrr_p, d_rrrr_p
    real(kind_dp) :: d_0_m, d_r_m, d_rr_m, d_rrr_m, d_rrrr_m

    n  = 1.0d0
    fn = 1.0d-5
    s  = 0.0d0
    fs = 0.0d0
    w  = 1.0d0

    d_0_p = vwn_energy(0.5d0*(n + fn + s + fs), 0.5d0*(n + fn - s - fs), w)
    d_0_m = vwn_energy(0.5d0*(n - fn + s - fs), 0.5d0*(n - fn - s + fs), w)

    call vwn_derv(0.5d0*n, 0.5d0*n, w, d)

    d_r    = d(1)
    d_rr   = d(2)
    d_rrr  = d(3)
    d_rrrr = d(4)

    call vwn_derv(0.5d0*(n + fn + s + fs), 0.5d0*(n + fn - s - fs), w, d)

    d_r_p    = d(1)
    d_rr_p   = d(2)
    d_rrr_p  = d(3)
    d_rrrr_p = d(4)

    call vwn_derv(0.5d0*(n - fn + s - fs), 0.5d0*(n - fn - s + fs), w, d)

    d_r_m    = d(1)
    d_rr_m   = d(2)
    d_rrr_m  = d(3)
    d_rrrr_m = d(4)

    write(0,*) 'd_r',    (d_0_p   - d_0_m)/(2.0d0*fn),   d_r
    write(0,*) 'd_rr',   (d_r_p   - d_r_m)/(2.0d0*fn),   d_rr
    write(0,*) 'd_rrr',  (d_rr_p  - d_rr_m)/(2.0d0*fn),  d_rrr
    write(0,*) 'd_rrrr', (d_rrr_p - d_rrr_m)/(2.0d0*fn), d_rrrr

  end subroutine

end module
