      SUBROUTINE GETH1(H,ISM,ITP,JSM,JTP)
*
* One-electron integrals over orbitals belonging to
* given OS class
*
*
* The orbital symmetries  are used to obtain the total
* symmetry of the one-electron integrals.
* It is therefore assumed that ISM, JSM represents a correct symmetry block
* of the integrals
*
* Jeppe Olsen, Version of fall 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
*.Global pointers
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
*.Output
      DIMENSION H(*)
*
C     IF
      NI = NOBPTS(ITP,ISM)
      NJ = NOBPTS(JTP,JSM)
      IJ = 0
      DO 100 J = 1, NJ
        DO 50 I = 1, NI
          IJ = IJ+1
          H(IJ) = GETH1E(I,ITP,ISM,J,JTP,JSM)
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
        CALL WRTMAT(H,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END
