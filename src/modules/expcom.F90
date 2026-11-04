!-----------------------------------------
! Module with variables previously defined
! at expcom.h include file (COMMON /EXPCOM/)
!-----------------------------------------
   module expcom_mod
   implicit none
   real*8 :: SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z, SIGN3X, SIGN3Y, SIGN3Z, &
             SIGN4X, SIGN4Y, SIGN4Z
   integer :: NCENT1, NCENT2, NCENT3, NCENT4, ICENT1, ICENT2, ICENT3, ICENT4, ISO1, ISO2, &
              ISO3, ISO4
   logical :: DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN
   end module expcom_mod
