      FUNCTION ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
*
* Adress of string IOCC
*
* version of Winter 1990 , Jeppe Olsen
*
      implicit real*8 (A-H,O-Z)
*
      INTEGER Z,IOCC(*)
      DIMENSION NEWORD(*),Z(NORB,*)
*
      IZ = 1
      DO 100 I = 1,NEL
        IZ = IZ + Z(IOCC(I),I)
  100 CONTINUE
*
      IF(IREORD.EQ.0) THEN
        ISTRNM = IZ
      ELSE
        ISTRNM = NEWORD(IZ)
      END IF
*
      NTEST = 0
      IF ( NTEST .GT. 1 ) THEN
        WRITE(6,*) ' STRING'
        CALL IWRTMA(IOCC,1,NEL,1,NEL)
C       WRITE(6,*) ' Z matrix '
C       CALL IWRTMA(Z,NORB,NEL,NORB,NEL)
        WRITE(6,*) ' ADDRESS OF STRING ',ISTRNM
        WRITE(6,*) ' REV LEX number : ', IZ
      END IF
*
      RETURN
      END
