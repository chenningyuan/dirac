!-----------------------------------------
! Module with variables previously defined
! at mxcent.h include file
!-----------------------------------------
!
!     MXATOM  = max number of different atomic types in the molecule
!     MXCENT  = max number of nuclei + point charges + ghost orbital centers in the molecule
!     MXBSETS = max number of basis sets (p.t. 5 : large comp., small comp. +
!               density fit sets for both + extended Huckel basis)
!     MAX_NONT= the maximum number of centers that can be defined within an atomic type (MI+PB)
!
   module mxcent_mod

   implicit none

   integer, parameter :: MXATOM   = 20
   integer, parameter :: MXCENT   = 300
   integer, parameter :: MXCOOR   = 3 * MXCENT
   integer, parameter :: MXBSETS  = 5
   integer, parameter :: MAX_NONT = 150

   end module mxcent_mod
