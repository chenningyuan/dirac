      Function I2EAD2(ISM,IORB,JSM,JORB,KSM,KORB,LSM,LORB)
*
* Find adress of integral in LUCIA order
*
* Indeces IORB,JORB,KORB,LORB are relative to start of symmetry block of orbs.
*
      IMPLICIT REAL*8           (A-H,O-Z)
*
#include "mxpdim.inc"
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR

#include "glbbas.inc"
*
#include "wrkspc.inc"
*
      I2EAD2 = I2EAD2S(ISM,IORB,JSM,JORB,KSM,KORB,LSM,LORB,WORK(KPINT2))
*
      RETURN
      END
