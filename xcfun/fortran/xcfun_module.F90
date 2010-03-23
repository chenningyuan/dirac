
module xc_fun
  implicit none
! These parameters should mirror those in xcfun.h
  integer, parameter ::  XC_VARS_A  = 0
  integer, parameter ::  XC_VARS_R  = 1
  integer, parameter ::  XC_VARS_AB = 2
  integer, parameter ::  XC_VARS_RS = 3
  integer, parameter ::  XC_LDA  = 0
  integer, parameter ::  XC_GGA  = 1
  integer, parameter ::  XC_MGGA = 2
contains
  ! We pass strings as null terminated integer arrays to C, this
  ! should be portable if C and Fortran uses the same character set.
  subroutine str2ints(str,ints)
! ints must have length at least len(str)+1
    character, intent(in) :: str*(*)
    integer, intent(out) :: ints(*)
    integer i
    do i=1,len(str)
       ints(i) = ichar(str(i:i))
    enddo
    ints(len(str)+1) = 0
  end subroutine str2ints

  subroutine ints2str(ints,str)
    character, intent(out) :: str*(*)
    integer, intent(in) :: ints(*)
    integer i,j
    i = 1
    do while (ints(i).ne.0)
       str(i:i) = char(ints(i))
       i = i + 1
    enddo
    do j=i,len(str)
       str(j:j) = ' '
    enddo
  end subroutine ints2str


  function xcfun_version()
    double precision xcfun_version, xcfunversionfortran
    xcfun_version = xcfunversionfortran()
  end function xcfun_version

  subroutine xcfun_splash(text)
    character, intent(out) :: text*(*)
    integer :: ibuf(len(text)+1)
    integer le
    le = len(text)
    allocate(ibuf(le+1))
    call xcsplashfortran(ibuf,le)
    call ints2str(ibuf,text)
  end subroutine xcfun_splash

  ! Create a new, "empty" functional and return its id
  function xc_new_functional()
    integer xc_new_functional, xcnewf
    xc_new_functional = xcnewf()
  end function xc_new_functional

  subroutine xc_set_setting(funid, setting, val, stat)
    integer, intent(in) :: funid
    integer, intent(out) :: stat
    double precision, intent(in) :: val
    character, intent(in) :: setting*(*)
    integer :: ibuf(len(setting)+1)
    str2ints(setting,ibuf)
    call xcsets(funid,ibuf,val)
  end subroutine xc_set_setting

  function xc_get_setting(funid,setting)
    integer, intent(in) :: funid
    double precision, intent(out) :: xc_get_setting
    double precision xcgets
    character, intent(in) :: setting*(*)
    integer :: ibuf(len(setting)+1)
    str2ints(setting,ibuf)
    xc_get_setting = xcgets(funid,ibuf)
  end function xc_get_setting

  function xc_is_set(funid,setting)
    integer, intent(in) :: funid
    logical, intent(out) :: xc_is_set
    character, intent(in) :: setting*(*)
    integer :: ibuf(len(setting)+1)
    integer xcisse,res
    str2ints(setting,ibuf)
    res = xcisse(funid,ibuf)
    if (res.eq.0) then
       xc_is_set = .true.
    else
       xc_is_set = .false.
    endif
  end function

  function xc_is_functional(funid,setting)
    integer, intent(in) :: funid
    logical, intent(out) :: xc_is_functional
    character, intent(in) :: setting*(*)
    integer :: ibuf(len(setting)+1)
    integer xcisfu,res
    str2ints(setting,ibuf)
    res = xcisfu(funid,ibuf)
    if (res.eq.0) then
       xc_is_functional = .true.
    else
       xc_is_functional = .false.
    endif
  end function

  subroutine xc_short_description(funid,setting,description)
    integer, intent(in) :: funid
    character, intent(out) :: description*(*)
    character, intent(in) :: setting*(*)
    integer :: ibuf(len(setting)+1)
    integer :: idescr(len(description)+1)
    str2ints(setting,ibuf)
    
  end subroutine

  subroutine xc_long_description(funid,setting,description)
    integer, intent(in) :: funid
    character, intent(out) :: description*(*)
    character, intent(in) :: setting*(*)
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
  end function xc_get_type

  subroutine xc_eval(funid, res, order, densvars)
    integer, intent(in) :: funid
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
    !FIXME
  end subroutine xc_eval

  function xc_index(funid, exponents)
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
