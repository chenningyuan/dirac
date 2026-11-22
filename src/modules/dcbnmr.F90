!-----------------------------------------
! Module with variables previously defined
! at src/include/dcbnmr.h include file (COMMON /CBLNMR/, /CBINMR/)
!-----------------------------------------
! ******************************************
! ***            NMR module              ***
! ***
!  NOTWOL - do not include two-electron London contribution
!  NOOONEI - do not include reorthonormalization contributions
!  INTNMR - what kind of integrals for LAO 2-el.contribution
! ******************************************
   module dcbnmr_mod
   implicit none
   integer, parameter :: LUTWOL = 69
   integer, parameter :: LUFCKL = 70
   integer :: INTNMR
   logical :: NOTWOL, NOONEI, EPREORTH
   end module dcbnmr_mod
