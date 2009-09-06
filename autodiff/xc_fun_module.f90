module xc_fun
  implicit none
  integer, parameter ::  XC_PARAMS_LEN = 6
  integer, parameter ::  XC_MAX_ORDER = 5 

  integer, parameter ::  XC_SLATER_EXCHANGE = 1
  integer, parameter ::  XC_VWN_CORRELATION = 2
  integer, parameter ::  XC_BECKE_EXCHANGE  = 3
  integer, parameter ::  XC_LYP_CORRELATION = 4
  integer, parameter ::  XC_PBE_EXCHANGE    = 5
  integer, parameter ::  XC_PBE_CORRELATION = 6

  integer, parameter :: XC_R = 1
  integer, parameter :: XC_S = 2
  integer, parameter :: XC_Z = 3
  integer, parameter :: XC_X = 4
  integer, parameter :: XC_Y = 5



contains
  subroutine set_functional(params)
    double precision, intent(in) :: params(XC_PARAMS_LEN)
    call xc_set_functional(params)
  end subroutine set_functional

  subroutine eval_lda(res, order, densvars)
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(2)
    call xc_eval_lda(res,order,densvars)
  end subroutine eval_lda

  subroutine eval_gga(res, order, densvars)
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(5)
    call xc_eval_gga(res,order,densvars)
  end subroutine eval_gga
end module
