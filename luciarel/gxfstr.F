      SUBROUTINE GXFSTR(ILA,ILB,IRA,IRB,ICA,IAA,ICB,IAB,NGAS,ICON)
*
* Left type defined by ILA, ILB
* Right type defined by IRA, IRB
* 
* Excitation operator defined by ICA, IAA, ICB,IAB
*
* Does this operator connect the two types ?
*
* ICON = 1 => Connection
* ICON = 0 => No connection
*
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
*. Input
      INTEGER ILA(NGAS),ILB(NGAS),IRA(NGAS),IRB(NGAS)
      INTEGER ICA(NGAS),ICB(NGAS),IAA(NGAS),IAB(NGAS)
*
      ICON = 1
      DO IGAS = 1, NGAS
        IF(.NOT.( ICA(IGAS) - IAA(IGAS) + IRA(IGAS).EQ.ILA(IGAS)))
     &  ICON = 0
        IF(.NOT.( ICB(IGAS) - IAB(IGAS) + IRB(IGAS).EQ.ILB(IGAS)))
     &  ICON = 0
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'GXFSTR : IRA, IRB, ILA, ILB '
        CALL IWRTMA(IRA,1,NGAS,1,NGAS)
        CALL IWRTMA(IRB,1,NGAS,1,NGAS)
        CALL IWRTMA(ILA,1,NGAS,1,NGAS)
        CALL IWRTMA(ILB,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' ICA, IAA, ICB, IAB '
        CALL IWRTMA(ICA,1,NGAS,1,NGAS)
        CALL IWRTMA(IAA,1,NGAS,1,NGAS)
        CALL IWRTMA(ICB,1,NGAS,1,NGAS)
        CALL IWRTMA(IAB,1,NGAS,1,NGAS)
        IF(ICON.EQ.1) THEN
          WRITE(6,*) ' We have contact '
        ELSE
          WRITE(6,*) ' No contact '
        END IF
      END IF
*
      RETURN
      END
