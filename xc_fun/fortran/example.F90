
! Show how to set up the PBE functional and calculate it at a point
program xc_example
  use xc_fun
  ! There are 21 (=bin(2+5,2)) 0:th, 1:st and 2:nd order coefficients
  double precision fun(XC_PARAMS_LEN), d(5), res(21)
  do i=1,XC_PARAMS_LEN
     fun(i) = 0
  enddo
  fun(XC_PBE_CORRELATION) = 1.0
  fun(XC_PBE_EXCHANGE) = 1.0
  call set_functional(fun)
  d(XC_R) = 1.0 !Total density
  d(XC_S) = 0.0 !rho_a - rho_b
  d(XC_Z) = 2.0 !|grad R|^2
  d(XC_X) = 1.0 !|grad S|^2
  d(XC_Y) = 1.0 !(grad R).(grad S)
  call eval_gga(res,2,d)
  print *,'XC energy:',res(1)
  print *,'dE_xc/dZZ:',res(GGA_ZZ)
end program xc_example
