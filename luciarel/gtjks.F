      SUBROUTINE GTJKS(J,K,NORB)
*
* Obtain Coulomb and Exchange integrals
* from complete integral list stored in core
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL*8           J(NORB,NORB),K(NORB,NORB)
*
      DO 200 IORB = 1, NORB
	DO 100 JORB = 1, NORB
	  J(IORB,JORB) = GTIJKL(IORB,IORB,JORB,JORB)
	  K(IORB,JORB) = GTIJKL(IORB,JORB,JORB,IORB)
  100   CONTINUE
  200 CONTINUE
*
      RETURN
      END
