      SUBROUTINE REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY,MXNOP)
* An operatorstring may be specifed as
*
* IOPGAS : Number of operators per GASspace
* IOP    : GASpace of each operator
*
* Transform between these two form
*
* Iway = 1 : IOPGAS => IOP
* Iway = 2 : IOP    => IOPGAS
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
*. Input/Output
      INTEGER IOPGAS(NGAS),IOP(MXNOP)
*
      IF(IWAY.EQ.1) THEN
        JOP = 0
        DO JGAS = 1, NGAS
          LJGAS = IOPGAS(JGAS)
          DO JJOP = 1, LJGAS
            JOP = JOP + 1
            IOP(JOP) = JGAS
          END DO
        END DO
        NOP = JOP
        if (NOP.gt.MXNOP) then
          write(6,*) 'NOP too large in REF_OP. NOP =',NOP
          write(6,*) 'Increase MXNOP or find bug !!'
        end if
      ELSE
        DO JGAS = 0, NGAS
          JOP = 0
          DO JJOP = 1, NOP
            IF(IOP(JJOP).EQ.JGAS) JOP = JOP +1 
          END DO
          IOPGAS(JGAS) = JOP
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' IOPGAS => IOP '
        ELSE
         WRITE(6,*) ' IOP => IOPGAS '
        END IF
        WRITE(6,*) ' IOPGAS and IOP '
        CALL IWRTMA(IOPGAS,1,NGAS,1,NGAS)
        CALL IWRTMA(IOP,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
