!-----------------------------------------
! Module with variables previously defined
! at trkoor.h include file (COMMON /TRKOOR/)
!-----------------------------------------
!     NCOOR  - total number of Cartesian coordinates
!     NCOORS - number of Cartesian coordinates of each symmetry
!     DEPEND - true for dependent symmetry coordinates
!     NCDEP  - total number of dependent/trarot coordinates
!     NCDEPS - number of dependent/trarot coordinates of each symmetry
!     IPTTRO - identifies symmetry-ordered trarot coordinate as
!              Tx, Ty, Tz, Rx, Ry, or Rz
   module trkoor_mod
   use mxcent_mod, only: MXCOOR
   implicit none
   integer :: NCOOR, NCOORS(8), NCDEP, NCDEPS(8), NCIND, NCINDS(8), IPTTRO(6,8), NTRREP(0:7)
   logical :: DEPEND(MXCOOR)
   end module trkoor_mod
