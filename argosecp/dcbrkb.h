!
! VERSION : $Revision: 1.5 $
! DATE    : $Date: 2004/05/26 13:06:54 $
! FILE    : dcbrkb.h
!
! ** COMMON block RKB pre-projection
!    IRKB - for large shells, pointer to L+1 small partner
!         - for small shells, pointer to L-2 small partner
      INTEGER NRKBSF, IRKB(MXPRIM)
      COMMON /RKBCOM/ NRKBSF, IRKB


