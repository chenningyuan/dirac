      SUBROUTINE NSTPTP_GAS(NGAS,ISPGRP,NSTSGP,NSMST,
     &                      NSTSSPGP,IGRP,MXNSTR)
*
* From number of strings per group and sym to number of strings
* per supergroup and sym for given super group.
*
* Jeppe Olsen , Fall of 94
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Output
      DIMENSION NSTSSPGP(NSMST,IGRP)
*. Scratch
#include "mxpdim.inc"
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NSTPTP_GAS is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' NSMST = ' , NSMST
        write(6,*) 'First 3 columns of NSTSGP'
        CALL IWRTMA(NSTSGP,NSMST,3,NSMST,3)
      END IF
*
      CALL ISETVC(NSTSSPGP(1,IGRP),0,NSMST)
*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
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
C       NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        ISMSPGP = ISM(1)
        NST = NSTSGP(ISM(1),ISPGRP(1))
        DO JGRP = 2, NGAS
          CALL SYMCOM(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTSGP(ISM(JGRP),ISPGRP(JGRP))
        END DO
C?      WRITE(6,*) ' Symmetry of groups '
C?      CALL IWRTMA(ISM,1,NGAS,1,NGAS)
C?      WRITE(6,*) ' Symmetry of super group and number of strings '
C?      WRITE(6,*) ISMSPGP , NST
        NSTSSPGP(ISMSPGP,IGRP) = NSTSSPGP(ISMSPGP,IGRP) + NST
        GOTO 1000
      END IF
*
      MXNSTR = -1
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
      END DO
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
      END IF
*
      RETURN
      END
