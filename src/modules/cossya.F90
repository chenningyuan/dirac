!-----------------------------------------
! Module with variables previously defined
! at cossya.h include file (COMMON /COSSYA0/, /COSSYA1/)
!-----------------------------------------
   module cossya_mod
   implicit none
   integer, parameter :: N_NFSYM = 16
   integer, parameter :: N_NREP = 32
   integer, parameter :: N_NOPEN = 16
   integer :: SYA_PRINT, SYA_INACT(N_NFSYM), SYA_ACTIV(N_NFSYM), SYA_NOPEN, &
              SYA_GASO(N_NREP,N_NOPEN), SYA_IELC(N_NOPEN), SYA_MAXE(N_NOPEN), &
              SYA_MINE(N_NOPEN), SYA_NTRDM(N_NREP), SYA_NREP, MSTATE
   logical :: KCOSINP, SYA_KTRDMAT
   end module cossya_mod
