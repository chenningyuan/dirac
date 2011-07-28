!
! FILE    : dcbrkb.h
!
! ** COMMON block RKB pre-projection
!    IRKB - for large shells, pointer to L+1 small partner
!         - for small shells, pointer to L-2 small partner
      INTEGER NRKBSF, IRKB(MXPRIM)
      COMMON /RKBCOM/ NRKBSF, IRKB


