      SUBROUTINE OCC_TO_GRP(IOCC,IGRP,IWAY)
*. Translate between occupation and group labels of supergroup
*
* IWAY = 1 => OCC to group
*      = 2 => Group to occ
*
* Number of strings in supergroup with given symmetry and number of 
* electrons in each space
*
* Jeppe Olsen, March 1999
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
C     INTEGER JGRP(MXPNGAS)
*. Specific input/output
      INTEGER IOCC(NGAS),IGRP(NGAS)
*
      NTEST = 000
*
      IF(IWAY.EQ.1) THEN 
*. Occupation => Group number 
        DO IOBTP = 1, NGAS
          JJGRP = 0
          DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
            IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
          END DO
          IGRP(IOBTP) = JJGRP
*
          IF(JJGRP.EQ.0.or.NTEST.ge.10) THEN
           WRITE(6,*) ' Group not included in list '
           WRITE(6,*) ' Input occupations : '
           CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
           write(6,*) 'Current orbital type :',IOBTP
           write(6,*) 'IBGPSTR array :'
           call iwrtma(IBGPSTR,1,NGAS,1,NGAS)
           write(6,*) 'NGPSTR array :'
           call iwrtma(NGPSTR,1,NGAS,1,NGAS)
           write(6,*) 'Group range ',
     &                 IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
           write(6,*) 'List of electron occupations: '
           do K=IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
             write(6,*) 'NELFGP(',K,') =', NELFGP(K)
           end do
           if (JJGRP.eq.0) STOP ' Group not included in list '
          END IF
        END DO
*
      ELSE
*. Group => Occupation 
        DO IOBTP = 1, NGAS
          IOCC(IOBTP) = NELFGP(IGRP(IOBTP))
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation and corresponding group array '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS) 
      END IF
*
      RETURN
      END 
