!-----------------------------------------
! Module with variables previously defined
! at src/krcc/cintfo.inc include file (COMMON /KRCC_CINTFO/)
!-----------------------------------------
!. At most 10 double group integral lists
   module cintfo_krcc_mod
   implicit none
   integer :: I12S, I34S, I1234S, NINT1, NINT2, NBINT1, NBINT2, NL1D, IUB1D(2,10), L1D(10), &
              NL2D, IUB2D(4,10), ISM2D(3,10), L2D(10)
   end module cintfo_krcc_mod
