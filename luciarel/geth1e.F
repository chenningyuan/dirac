      FUNCTION GETH1E(IORB,ITP,ISM,JORB,JTP,JSM)
*
* One-electron integral for active
* orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
*
* The orbital symmetries are used to obtain the
* total symmetry of the operator
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
*
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "multd2h.inc"
C     INCLUDE 'intform.inc'
*. Normal integrals
      IH1FORM = 1
      IJSM = MULTD2H(ISM,JSM)
      IF(IH1FORM.EQ.1) THEN
*. Normal integrals, lower triangular packed
        IF(IJSM.EQ.1) THEN
          GETH1E =
     &    GTH1ES(IREOTS,WORK(KPINT1),WORK(KINT1),IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        ELSE
C         GETH1E =
C    &    GTH1ES(IREOTS,WORK(KPGINT1(IJSM)),WORK(KINT1),IBSO,MXPNGAS,
C    &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        END IF
      ELSE
*. Integrals are in full blocked form
C       GETH1E =
C    &  GTH1ES(IREOTS,WORK(KPGINT1A(IJSM)),WORK(KINT1),IBSO,MXPNGAS,
C    &         IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,0)
      END IF
*
      RETURN
      END
