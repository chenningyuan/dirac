      SUBROUTINE FIFAM(FIFA)
*
*. Construct inactive + active fock matrix
*
*. On input FIFAM Should be the inactive Fock matrix, in symmetry packed form
*
* Jeppe Olsen
      ImplICIT REAL*8(A-H,O-Z)
*
      DIMENSION FIFA(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "lucinp.inc"
*
*
      CALL FIFAMS(FIFA,WORK(KRHO1),IBSO,NSMOB,
     &            NTOOBS,NACOB,NTOOB,IREOST)
*
      RETURN
      END
