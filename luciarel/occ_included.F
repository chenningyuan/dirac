      SUBROUTINE OCC_INCLUDED(IOCC,INCLUDED)
*
* A supergroup is given as occupation for each gas space
*
* Check if all these groups are included, and return 
* answer in INCLUDED
*
* Jeppe and Timo, Oct 2001 
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
*. Specific input
      INTEGER IOCC(NGAS)
*
      NTEST = 000
*
      INCLUDED = 1
      DO IOBTP = 1, NGAS
        JJGRP = 0
        DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
          IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
        END DO
        IF(JJGRP.EQ.0 ) INCLUDED = 0
      END DO
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation  '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Included = ',INCLUDED
      END IF
*
      RETURN
      END 
