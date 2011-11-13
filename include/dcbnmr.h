! ******************************************
! ***            NMR module              ***
! ***
! ***   LUTWOL - file number for two-ele.contrib. of LAO
! ***   LUFCKL - file number  of LAO Fock matrix
! 
!  NOTWOL - do not include two-electron London contribution
!  NOOONEI - do not include reorthonormalization contributions
!  DOEPRN  - do include positronic shells into 
!            the summation for reorthonormalization
!  INTNMR - what kind of integrals for LAO 2-el.contribution
! 
! ******************************************
      PARAMETER (LUTWOL = 69,LUFCKL = 70)
      LOGICAL NOTWOL, NOONEI, DOEPRN
      COMMON/CBLNMR/NOTWOL, NOONEI, DOEPRN
      INTEGER INTNMR
      COMMON/CBINMR/INTNMR
