
#if 0
Fortran version of xc_fun. Note that the indices XC_* are increased
by one to match with fortran indexing.

Fortran interface (order is integer, the others 
		   double precision):

subroutine xc_set_functional(params)
subroutine xc_eval_lda(result, order, densvars)
subroutine xc_eval_gga(result, order, densvars)

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
