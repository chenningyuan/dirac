      SUBROUTINE NST_SPGRP_DBG(NIGRP,IGRP,ISM_TOT,NSMST,NSTRIN,NDIST,
     &                         IAB)
*
* Number of strings for given combination of groups and 
* symmetry.
*
*. Input
*        
*
*   NIGRP : Number of active groups 
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSMST   : Number of string symmetries
*   IAB : = 1 => unbarred operators
*         = 2 => barred operators 
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
*. Double group version
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*. Specific Input
      DIMENSION IGRP(NIGRP)
*. General input
#include "gasstr.inc"
C     DIMENSION NSTSGP(NSMST,*)
*. Scratch 
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' NST_SPGRP_DBG is speaking '
        WRITE(6,*) ' =========================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
      END IF
*
      IF(NIGRP.EQ.0) THEN
        IF(ISM_TOT.EQ.1) THEN 
          LENGTH = 1
        ELSE
          LENGTH = 0
        END IF
        GOTO 1001
      END IF
*
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST_DBG(NIGRP,IGRP,MNSM,MXSM,NDISTX,IAB)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0 
      NDIST = 0
*. Last group with symmetry differing from total symmetric
      NIGRPL = 1
      DO JGRP = 1, NIGRP
        IF(MXSM(JGRP).GT.1) NIGRPL = JGRP
      END DO
*. Number of strings in groups after NIGRPL
      NSTRL = 1
      DO JGRP = NIGRPL+1,NIGRP
        IF(IAB.EQ.1) THEN 
          NSTRL = NSTRL*NSTFSMGP(1,IGRP(JGRP))
        ELSE  
          NSTRL = NSTRL*NSTFSMGP2(1,IGRP(JGRP))
        END IF
      END DO
 1000 CONTINUE
*. Next symmetry distribution
* \/  NEXT_SYM_DISTR should work for double group symmetry
        CALL NEXT_SYM_DISTR(NIGRPL,MNSM,MXSM,ISM,ISM_TOT,IFIRST,NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = NSTRL
          DO JGRP = 1, NIGRPL
            IF(IAB.EQ.1) THEN
              LDIST = LDIST*NSTFSMGP(ISM(JGRP),IGRP(JGRP))
            ELSE 
              LDIST = LDIST*NSTFSMGP2(ISM(JGRP),IGRP(JGRP))
            END IF
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
 1001 CONTINUE
      NSTRIN = LENGTH
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', LENGTH
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
*
      RETURN
      END
