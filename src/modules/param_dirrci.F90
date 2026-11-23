!-----------------------------------------
! Module with variables previously defined
! at src/dirrci/param.inc include file (COMMON /MAXMEM/)
!-----------------------------------------
! SOURCE  : $Source$
!   DEFINITION OF PARAMETERS DIRRCI
!   N2  MAXIMUM NUMBER OF VECTORS
!   N3  MAXIMUM NUMBER OF ORBITALS
!   N4  MAXIMUM NUMBER OF ELECTRONS
!   N8  MAXIMUM NUMBER OF ROOTS OPTIMIZED SIMULTANEOUSLY
   module param_dirrci_mod
   implicit none
   integer, parameter :: N2 = 100
   integer, parameter :: N3 = 300
   integer, parameter :: N4 = 80
   integer, parameter :: N8 = 36
   integer :: MAXCORE
   end module param_dirrci_mod
