! orgcom.h :
!   CMXYZ  : center-of-mass
!   NCHBC  : nuclei charge barycenter
!   ORIGIN : general one-electron operator origin
!   DIPORG : origin for multipole operators, including dipole
!   GAGORG : magnetic gauge origin
!     GAGORG_SET true: GAGORG has been set by user
!   CAVORG : cavity origin for spherical cavity model
      REAL*8 CMXYZ(3), ORIGIN(3), DIPORG(3), GAGORG(3), CAVORG(3)
      REAL*8 NCHBC(3)
      LOGICAL GAGORG_SET
      COMMON /ORGCOM/ CMXYZ, ORIGIN, DIPORG, GAGORG, CAVORG,            &
     &   GAGORG_SET, NCHBC
