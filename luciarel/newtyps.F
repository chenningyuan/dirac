      SUBROUTINE NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &           NGAS,ISPGPAN,ISPGPCR,OUTSPGP)
*
* Strinf of operators X supergroup => new supergroup
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ISPGPAN(NGAS,*),ISPGPCR(NGAS,*)
*. Specific input
      DIMENSION IACOP(NOP),ITPOP(NOP)
*. Output
      INTEGER OUTSPGP
*
      OUTSPGP = INSPGP
      DO IOP = 1, NOP
        IF(IACOP(IOP).EQ.1) THEN
          OUTSPGP = ISPGPAN(ITPOP(IOP),OUTSPGP)
        ELSE
          OUTSPGP = ISPGPCR(ITPOP(IOP),OUTSPGP)
        END IF
*
        IF(OUTSPGP.EQ.0) THEN
C         WRITE(6,*) ' NEWTYPS, cul de sac'
C         WRITE(6,*) ' undefined supergroup type '
C         WRITE(6,*) ' String of operator : IAC,ITPOP'
C         CALL IWRTMA(IACOP,1,NOP,1,NOP)
C         CALL IWRTMA(ITPOP,1,NOP,1,NOP)
          GOTO 1001
        END IF
      END DO
*
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NEWTYPS : IAC and ITPOP'
        CALL IWRTMA(IACOP,1,NOP,1,NOP)
        CALL IWRTMA(ITPOP,1,NOP,1,NOP)
        WRITE(6,*) ' Input and Output group ',INSPGP,OUTSPGP
      END IF
*
      RETURN
      END
