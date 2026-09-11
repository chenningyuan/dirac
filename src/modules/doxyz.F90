!-----------------------------------------
!
! Module with variables previously defined
!
! at doxyz.h include file
!
!-----------------------------------------

module doxyz_mod

   implicit none

   logical :: DOXYZ(3)

   logical :: DOX, DOY, DOZ

   equivalence (DOX, DOXYZ(1)), (DOY, DOXYZ(2)), (DOZ, DOXYZ(3))

end module doxyz_mod
