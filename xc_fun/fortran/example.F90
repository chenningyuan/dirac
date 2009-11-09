
! Show how to set up the PBE functional and calculate it at a point
program xc_example
  use xc_fun
  ! There are 21 (=bin(2+5,2)) 0:th, 1:st and 2:nd order coefficients
  double precision fun(XC_PARAMS_LEN), d(5), res(21)
  fun = 0
  fun(XC_PBE_CORRELATION) = 1.0
!  fun(XC_PBE_EXCHANGE) = 1.0
  call xc_set_functional(XC_ABFGH,fun)
  d(1) = 1.7 !rho_a
  d(2) = 1.7 !rho_b
  d(3) = 1.7 !|grad rho_a|^2
  d(4) = 1.7 !|grad rho_b|^2
  d(5) = 1.7 !(grad rho_a).(grad rho_b)
  call xc_eval(res,2,d)
  print *,'Calculated',xc_len(XC_ABFGH,2),'terms, for example:'
  print *,'XC energy:',res(1)
  print *,'d^2/drhoa^2:',res(xc_index(XC_ABFGH,(/2,0,0,0,0/)))
  print *,'d^2/drhoa drhob:',res(xc_index(XC_ABFGH,(/1,1,0,0,0/)))
end program xc_example
