
#if 0
Fortran version of xc_fun. Note that the indices XC_* are increased
by one to match with fortran indexing.

Fortran interface:

subroutine xc_set_functional(params)
  double precision, dimension(XC_PARAMS_LEN), intent(in) :: densvars

subroutine xc_eval_gga(result, order, densvars)
  double precision, dimension(:), intent(out) :: result
  integer, intent(in) :: order
  double precision, dimension(5), intent(in) :: densvars

subroutine xc_eval_lda(result, order, densvars)
  Same as xc_eval_gga, but dimension of densvars is 2 instead of 5

This file must be kept in sync with xc_fun.h
#endif

#define XC_PARAMS_LEN      6
#define XC_MAX_ORDER       5 

#define XC_SLATER_EXCHANGE 1
#define XC_VWN_CORRELATION 2
#define XC_BECKE_EXCHANGE  3
#define XC_LYP_CORRELATION 4
#define XC_PBE_EXCHANGE    5
#define XC_PBE_CORRELATION 6
