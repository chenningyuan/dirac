      SUBROUTINE NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
*
* an input supergroup is given.
* apply a string of elementary operators to this supergroup and
* obtain supergroup number of new string
*
* Jeppe Olsen, October 1993
* GAS-version : July 95
*
* ------
* Input
* ------
*
* INSPGP  : input super group ( given occupation in each AS )
* IACOP(I) = 1 : operator I is an annihilation operator
*          = 2 : operator I is a  creation   operator
* ITPOP(I) : orbitals space of operator I
* NOP : Number of operators
*
* Output
* ------
* OUTSPGP  : supergroup of resulting string
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER ITPOP(*),IACOP(*)
*. Number of active spaces  (NGAS )
#include "mxpdim.inc"
#include "cgas.inc"
#include "strbas.inc"
#include "wrkspc.inc"
*. Local scratch
      DIMENSION IEL(MXPNGAS)
*. output
      INTEGER OUTSPGP
*
      INEW_OR_OLD = 1
      IF(INEW_OR_OLD.EQ.1) THEN
        CALL NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &       NGAS,WORK(KSPGPAN),WORK(KSPGPCR),OUTSPGP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
      ELSE IF(INEW_OR_OLD.EQ.2) THEN
*. Number of electrons in AS1,AS2, .... for input supergroup
        CALL  GTSPGP(IEL,INSPGP,2)
*
        IDELTA = 0
        DO IOP = 1, NOP
*. Change in number of orbitals
          IF(IACOP(IOP).EQ.1) THEN
            IDELTA = IDELTA - 1
            IDELTAI = -1
          ELSE
            IDELTA = IDELTA + 1
            IDELTAI = +1
          END IF
          IEL(ITPOP(IOP)) = IEL(ITPOP(IOP)) + IDELTAI
        END DO
*. output supergroup
        CALL  GTSPGP(IEL,OUTSPGP,1)
      END IF
*     /\ End  of NEW_OR_OLD switch
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NEWTYP ,  OUTSPGP ', OUTSPGP
      END IF
*
      RETURN
      END
