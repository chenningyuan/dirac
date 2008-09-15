!dirac_copyright_start
!dirac_copyright_end
!
!      here we can collect variable type definitions
!      i'm not so sure about the decimals and exponent ranges so please edit
!
!      radovan bast (2008-04-05)
!      inspiration: discussions with jonas juselius
!-------------------------------------------------------------------------------

module kinds_m

  implicit none

  integer, public, parameter :: kind_i2 = selected_int_kind(4)
  integer, public, parameter :: kind_i4 = selected_int_kind(8)

  integer, public, parameter :: kind_sp = selected_real_kind(6, 30)
  integer, public, parameter :: kind_dp = selected_real_kind(14, 200)

end module kinds_m
