      SUBROUTINE NST_SPGP(IOCC,NSTFSM)
* Number of strings for given supergroup. 
*.Input supergroup is defined by  occupation in each orb space
*
* Jeppe Olsen , March 99   
*
* Does not work for double group symmetries !
* Use nst_spgrp_dbg and loop over symmetries instead.
*   Timo Fleig
* 
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "cgas.inc"
*. Input
      DIMENSION IOCC(NGAS)
*. Output
      DIMENSION NSTFSM(NSMST)        
*. Scratch 
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER IGRP(MXPNGAS)
*
      NTEST = 100
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NST_SPGP is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Occupation of supergroup '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
      END IF
*. Occupation => groups
C  OCC_TO_GRP_(NSTR,IOCC,IGRP,IWAY)
      CALL OCC_TO_GRP(IOCC,IGRP,1)
      IZERO = 0
      CALL ISETVC(NSTFSM,IZERO,NSMST)
*
*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. First symmetry combination
      DO IGAS = 1, NGAS
         ISM(IGAS) = MNSM(IGAS)
      END DO
*. Loop over symmetries in each gas space
      IFIRST = 1
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
      ELSE
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in 
*  this supergroup
        ISMSPGP = ISM(1)
        NST = NSTFSMGP(ISM(1),IGRP(1))
        DO JGRP = 2, NGAS
          CALL SYMCOM(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTFSMGP(ISM(JGRP),IGRP(JGRP))
        END DO
        NSTFSM(ISMSPGP) =   NSTFSM(ISMSPGP) + NST
        GOTO 1000
      END IF         
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of strings per symmetry for supergroup'
        CALL IWRTMA10(NSTFSM(1),1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
