!-----------------------------------------
! Module with variables previously defined
! at dcbnmr.h include file
!-----------------------------------------
   module dcbnmr_mod
!  ******************************************
!  ***            NMR module              ***
!  ***
!  ***   LUTWOL - file number for two-ele.contrib. of LAO
!  ***   LUFCKL - file number  of LAO Fock matrix
! 
!   NOTWOL - do not include two-electron London contribution
!   NOONEI - do not include reorthonormalization contributions
!   INTNMR - what kind of integrals for LAO 2-el.contribution
! 
!  ******************************************
   implicit none
!
   !PARAMETERS
   INTEGER, PARAMETER :: LUTWOL = 69
   INTEGER, PARAMETER :: LUFCKL = 70
!
   !VARIABLES
   INTEGER :: INTNMR
   LOGICAL :: NOTWOL, NOONEI, EPREORTH
!
   end module dcbnmr_mod