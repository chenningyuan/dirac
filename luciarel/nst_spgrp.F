      SUBROUTINE NST_SPGRP(NGRP,IGRP,ISM_TOT,NSTSGP,NSMST,
     &                     NSTRIN,NDIST)
*
* Number of strings for given combination of groups and
* symmetry.
*
*. Input
*
*
*   NGRP : Number of active groups
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSTSGP  : Number of strings per symmetry and supergroup
*   NSMST   : Number of string symmetries
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific Input
      DIMENSION IGRP(NGRP)
*. General input
      DIMENSION NSTSGP(NSMST,*)
*. Scratch
#include "mxpdim.inc"
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' NST_SPGP is speaking '
        WRITE(6,*) ' ====================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NGRP)
      END IF
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST(NGRP,IGRP,MNSM,MXSM,NDISTX)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0
      NDIST = 0
 1000 CONTINUE
*. Next symmetry distribution
        CALL NEXT_SYM_DISTR(NGRP,MNSM,MXSM,ISM,ISM_TOT,IFIRST,NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = 1
          DO JGRP = 1, NGRP
            LDIST = LDIST*NSTSGP(ISM(JGRP),IGRP(JGRP))
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
      NSTRIN = LENGTH
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', LENGTH
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
*
      RETURN
      END
