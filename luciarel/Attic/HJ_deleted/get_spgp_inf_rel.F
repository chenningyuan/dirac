      SUBROUTINE GET_SPGP_INF_REL(ISPGP,ITP,IGRP)
*
* Obtain groups defining supergroup ISPGP,ITP
*
* Jeppe Olsen, November 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Output
      DIMENSION IGRP(*)
*
      NTEST = 000
*. Absolute group number
      ISPGPABS = ISPGP + IBSPGPFTP(ITP) -1
      CALL ICOPVE(ISPGPFTP(1,ISPGPABS),IGRP,NGAS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_SPGP_INF_REL : ISPGP ITP ISPGPABS',
     &              ISPGP, ITP, ISPGPABS
        WRITE(6,*) ' Groups of supergroups'
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
