!-----------------------------------------
! Module with parameters previously defined
! at src/gp/xml/xml_settings.h include file
!-----------------------------------------
module xml_settings_mod
  implicit none
  integer, parameter :: name_len = 16
  integer, parameter :: attr_len = 16
  integer, parameter :: data_len = 16

  integer, parameter :: real_kind = 8

! For layout
! CJ: the e3 at the end of the format is to use 3 digits for the exponent
  character(len=80), parameter :: real_fmt = '(e28.16e3)'
  integer, parameter :: tab_default = 3
end module xml_settings_mod
