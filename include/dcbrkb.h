C ** COMMON block RKB pre-projection
C    NRKBSF - Number of scalar RKB S-component functions
C    CRKB   - Transformation matrix Cartesian -> RKB
C    CRKBBK - Transformation matrix RKB -> Cartesian
C    ISRKB  - Pointer to Transformation matrix (CRKB) for given L
C    IRKBLS - Array of S-comp functions that never are
C             L, only L-2 partner in pair. Lenght is IRKBLE
C    IRKBPT - Pointer array to L-2 partner for S function
C    LRKBPT - Pointer array to L-1 partner for L function
C    
#include <sphtrm.h>
      COMMON /RKBCOM/ CRKB(NCSP),CRKBBK(NCSP),ISRKB(MXQN),
     &  NRKBSF, ISTRT(MXSHEL),IRKBLS(MXSHEL),IRKBLE,
     &  IRKBPT(MXSHEL,2),LRKBPT(MXSHEL,2)

