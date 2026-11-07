!-----------------------------------------
! Module with variables previously defined
! at molinp.h include file (COMMON /MOLINP/, /MOLINP_C/)
!-----------------------------------------
!  input .mol file is saved in MLINE(1:NMLINE)
   module molinp_mod
   use mxcent_mod, only: MXCENT
   implicit none
   integer, parameter :: KMLINE = 1000
   integer, parameter :: LEN_MLINE = 80
   character(len=LEN_MLINE) :: MLINE(KMLINE)
   integer :: NMLINE, NCLINE(MXCENT), NMLAU
   end module molinp_mod
