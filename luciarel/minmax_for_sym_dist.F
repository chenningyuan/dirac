      SUBROUTINE MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Include blocks
#include "mxpdim.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "wrkspc.inc"
*. Input
      DIMENSION IGRP(NIGRP)
*.Output
      DIMENSION MNVAL(NIGRP),MXVAL(NIGRP)
*. Local scratch
      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*
      NTEST = 0000
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
*. Number of strings per sym and group
      DO JGRP = 1, NIGRP
        CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(JGRP)-1)*NSMST+1,
     &               NSMST,LSMGP(1,JGRP))
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' LSMGP '
        CALL IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
      END IF
*. Max and min sym in each group
      DO JGRP = 1, NIGRP
*
        IMAX = 1
        DO ISM=1, NSMST
          IF(LSMGP(ISM,JGRP).GT.0) IMAX = ISM
        END DO
        MXVAL(JGRP) = IMAX
*
        IMIN = NSMST
        DO ISM = NSMST,1,-1
         IF(LSMGP(ISM,JGRP).GT.0) IMIN = ISM
        END DO
        MNVAL(JGRP) = IMIN
      END DO
*. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Group combination : '
        WRITE(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1, NIGRP)
        WRITE(6,*)
        WRITE(6,*) ' Group Minsym Maxsym'
        WRITE(6,*) ' ==================='
        DO JGRP = 1, NIGRP
          WRITE(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Total number of distributions', NDIST
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
*
      RETURN
      END
