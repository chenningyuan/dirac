      SUBROUTINE MINMAX_FOR_SYM_DIST_DBG(NIGRP,IGRP,MNVAL,MXVAL,NDIST,
     &                                   IAB)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*              April 1998     From  MINMAX_SM_GP
*
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
*
      DO JGRP = 1, NIGRP
C Requires setting up of  MINMAX_SM_GP and  MINMAX_SM_GP2 
        IF(IAB.EQ.1) THEN
          MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
          MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
        ELSE
          MNVAL(JGRP) = MINMAX_SM_GP2(1,IGRP(JGRP))
          MXVAL(JGRP) = MINMAX_SM_GP2(2,IGRP(JGRP))
        END IF
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
