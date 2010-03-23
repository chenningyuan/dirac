
module xc_fun
  implicit none
contains
  function xc_new_functional()
    integer xc_new_functional, xcnewfunctional
    xc_new_functional = xcnewfunctionalfortran()
  end function xc_new_functional

  subroutine xc_set_setting(funid, setting, val, stat)
    integer, intent(in) :: funid
    integer, intent(out) :: stat
    double precision, intent(in) :: val
    character, intent in :: setting*(*)
  end subroutine xc_set_setting

  function xc_get_setting(funid,setting)
    integer, intent(in) :: funid
    double precision, intent(out) :: xc_get_setting
    character, intent in :: setting*(*)
  end function xc_get_setting

  function xc_is_set(funid,setting)
    integer, intent(in) :: funid
    logical, intent(out) :: xc_get_setting
    character, intent in :: setting*(*)
  end function

  subroutine xc_short_description(funid,setting,description)
    integer, intent(in) :: funid
    character, intent(out) :: description*(*)
    character, intent in :: setting*(*)
  end subroutine

  subroutine xc_long_description(funid,setting,description)
    integer, intent(in) :: funid
    character, intent(out) :: description*(*)
    character, intent in :: setting*(*)
  end subroutine


  subroutine xc_set_mode(funid, mode)
    integer, intent(in) :: funid, mode
  end subroutine xc_set_mode

  function xc_get_type(funid)
    integer, intent(in) :: funid
    integer xc_get_type
  end subroutine xc_get_type

  subroutine xc_eval(funid, res, order, densvars)
    integer, intent(in) :: funid
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
  end subroutine xc_eval

  function xc_index(funid, exponents)    integer, intent(in) :: 
    integer, intent(in) :: exponents(*)
    integer xc_index, xc_index_fortran
    xc_index = xcindexfortran(mode,exponents)
  end function xc_index

  function xc_input_length(mode, order)
    integer, intent(in) :: mode
    integer, intent(in) :: order
    integer xc_len, xc_len_fortran
    xc_len = xclenfortran(mode,order)
  end function xc_input_length

  function xc_output_length(mode, order)
    integer, intent(in) :: mode
    integer, intent(in) :: order
    integer xc_len, xc_len_fortran
  end function xc_output_length
end module
