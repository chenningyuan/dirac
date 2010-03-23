
module xc_fun
  implicit none
! These parameters should mirror those in xcfun.h
  integer, parameter ::  XC_VARS_A  = 0
  integer, parameter ::  XC_VARS_R  = 1
  integer, parameter ::  XC_VARS_AB = 2
  integer, parameter ::  XC_VARS_RS = 3
  integer, parameter ::  XC_LDA  = 0
  integer, parameter ::  XC_GGA  = 0
  integer, parameter ::  XC_MGGA = 0
contains
  function xcfun_version()
    double precision xcfun_version, xcfunversionfortran
    xcfun_version = xcfunversionfortran()
  end function xcfun_version

  subroutine xcfun_splash(text)
    character, intent(out) :: text*(*)
    ! FIXME
  end subroutine xcfun_splash

  ! Create a new, "empty" functional and return its id
  function xc_new_functional()
    integer xc_new_functional, xcnewfunctionalfortran
    xc_new_functional = xcnewfunctionalfortran()
  end function xc_new_functional

  subroutine xc_set_setting(funid, setting, val, stat)
    integer, intent(in) :: funid
    integer, intent(out) :: stat
    double precision, intent(in) :: val
    character, intent in :: setting*(*)
    !FIXME
  end subroutine xc_set_setting

  function xc_get_setting(funid,setting)
    integer, intent(in) :: funid
    double precision, intent(out) :: xc_get_setting
    character, intent in :: setting*(*)
    !FIXME
  end function xc_get_setting

  function xc_is_set(funid,setting)
    integer, intent(in) :: funid
    logical, intent(out) :: xc_get_setting
    character, intent in :: setting*(*)
    !FIXME
  end function

  subroutine xc_short_description(funid,setting,description)
    integer, intent(in) :: funid
    character, intent(out) :: description*(*)
    character, intent in :: setting*(*)
    !FIXME
  end subroutine

  subroutine xc_long_description(funid,setting,description)
    integer, intent(in) :: funid
    character, intent(out) :: description*(*)
    character, intent in :: setting*(*)
    !FIXME
  end subroutine

  subroutine xc_set_mode(funid, mode)
    integer, intent(in) :: funid, mode
    !FIXME
  end subroutine xc_set_mode

  function xc_get_type(funid)
    integer, intent(in) :: funid
    integer xc_get_type
    !FIXME
  end subroutine xc_get_type

  subroutine xc_eval(funid, res, order, densvars)
    integer, intent(in) :: funid
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
    !FIXME
  end subroutine xc_eval

  function xc_index(funid, exponents)    integer, intent(in) :: 
    integer, intent(in) :: exponents(*)
    integer xc_index, xc_index_fortran
    xc_index = xcindexfortran(mode,exponents)
  end function xc_index

  function xc_input_length(funid)
    integer, intent(in) :: funid
    integer xc_input_length, xcinlenfortran
    xc_input_length = xcinlenfortran(funid)
  end function xc_input_length

  function xc_output_length(funid, order)
    integer, intent(in) :: funid, order
    integer xc_output_length, xcinlenfortran
    xc_input_length = xcinlenfortran(funid,order)
  end function xc_output_length
end module
