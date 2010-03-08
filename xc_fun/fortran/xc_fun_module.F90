!The following variables are used in the code:
! A = \rho_\alpha
! B = \rho_\beta
! F = \nabla A \cdot \nabla A (usually called \sigma_\alpha\alpha)
! G = \nabla B \cdot \nabla B (\sigma_\beta\beta)
! H = \nabla A \cdot \nabla B (\sigma_\alpha\beta)
!Alternatively you can use;
! R = \rho_\alpha + \rho_\beta
! S = \rho_\alpha - \rho_\beta
! Z = \nabla R \cdot \nabla R (sometimes called \sigma)
! X = \nabla S \cdot \nabla S
! Y = \nabla R \cdot \nabla S

module xc_fun
  implicit none

#ifndef XC_FUN_INTERFACE_F90
#define XC_FUN_INTERFACE_F90
#endif
#include "../include/xc_fun_interface.h"

contains
! Set up the functional used by xc_eval. The mode is one of 
! XC_A,XC_R, .. defined above. It defines the variables that
! will be used as input and output to xc_eval. For example
! using XC_ABFGH will require densvars to be an array of
! the five variables A,B,F,G,H described above. The output
! derivatives are given wrt these same variables. 
  subroutine xc_set_functional(mode,params)
    integer, intent(in) :: mode
    double precision, intent(in) :: params(XC_PARAMS_LEN)
    call xc_set_functional_fortran(mode,params,XC_PARAMS_LEN)
  end subroutine xc_set_functional

! Eval the functional previously set by xc_set_functional.
! The size and layout of res and densvars depends on the
! mode given to xc_set_functional.
  subroutine xc_eval(res, order, densvars)
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
    call xc_eval_fortran(res,order,densvars)
  end subroutine xc_eval

#ifdef WITH_QD
! Use quad-double internally but round the output to double.
  subroutine xc_eval_qd(res, order, densvars)
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
    call xc_eval_fortran_qd(res,order,densvars)
  end subroutine xc_eval_qd
#endif
  function xc_index(mode, exponents)
    integer, intent(in) :: mode
    integer, intent(in) :: exponents(*)
    integer xc_index, xc_index_fortran
    xc_index = xc_index_fortran(mode,exponents)
  end function xc_index

  function xc_len(mode, order)
    integer, intent(in) :: mode
    integer, intent(in) :: order
    integer xc_len, xc_len_fortran
    xc_len = xc_len_fortran(mode,order)
  end function xc_len
end module
