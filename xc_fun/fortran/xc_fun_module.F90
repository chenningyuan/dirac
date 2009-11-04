module xc_fun
  implicit none
  integer, parameter ::  XC_PARAMS_LEN = 9
  integer, parameter ::  XC_MAX_ORDER = 5 

! Index of each functional parameter in the argument
! of set_functional
  integer, parameter ::  XC_SLATER_EXCHANGE = 1
  integer, parameter ::  XC_VWN5_CORRELATION = 2
  integer, parameter ::  XC_BECKE_EXCHANGE  = 3
  integer, parameter ::  XC_LYP_CORRELATION = 4
  integer, parameter ::  XC_PBE_EXCHANGE    = 5
  integer, parameter ::  XC_PBE_CORRELATION = 6
  integer, parameter ::  XC_RPBE_EXCHANGE    = 7
  integer, parameter ::  XC_OLD_PBE_CORRELATION = 8
  integer, parameter ::  XC_PW92_CORRELATION = 9

! Indices into the densvars array (lda_eval and gga_eval)
  integer, parameter :: XC_R = 1
  integer, parameter :: XC_S = 2
  integer, parameter :: XC_Z = 3
  integer, parameter :: XC_X = 4
  integer, parameter :: XC_Y = 5

! Indices for the LDA first and second derivatives in the 
! output from lda_eval()
  integer, parameter :: LDA_R = 2
  integer, parameter :: LDA_S = 3
  integer, parameter :: LDA_RR = 4
  integer, parameter :: LDA_RS = 5
  integer, parameter :: LDA_SS = 6

! Indices for the GGA first and second derivatives in the 
! output from gga_eval()
  integer, parameter :: GGA_R = 2
  integer, parameter :: GGA_S = 3
  integer, parameter :: GGA_Z = 4
  integer, parameter :: GGA_X = 5
  integer, parameter :: GGA_Y = 6

  integer, parameter :: GGA_RR = 7
  integer, parameter :: GGA_RS = 8
  integer, parameter :: GGA_RZ = 9
  integer, parameter :: GGA_RX = 10
  integer, parameter :: GGA_RY = 11
  integer, parameter :: GGA_SS = 12
  integer, parameter :: GGA_SZ = 13
  integer, parameter :: GGA_SX = 14
  integer, parameter :: GGA_SY = 15
  integer, parameter :: GGA_ZZ = 16
  integer, parameter :: GGA_ZX = 17
  integer, parameter :: GGA_ZY = 18
  integer, parameter :: GGA_XX = 19
  integer, parameter :: GGA_XY = 20
  integer, parameter :: GGA_YY = 21


contains
  subroutine set_functional(params)
    double precision, intent(in) :: params(XC_PARAMS_LEN)
    call xc_set_functional(params,XC_PARAMS_LEN)
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
