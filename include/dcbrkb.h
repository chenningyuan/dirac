C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbrkb.h
C
C ** COMMON block RKB pre-projection
C    IRKB - for large shells, pointer to L+1 small partner
C         - for small shells, pointer to L-2 small partner
      COMMON /RKBCOM/ NRKBSF,IRKB(MXPRIM)


