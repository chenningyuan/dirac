      SUBROUTINE GRAPW(W,Y,MINEL,MAXEL,NORB,NEL,NTEST)
*
* A graph of strings has been defined from
*
*      MINEL(I) is the smallest allowed number of electrons in
*      orbitals 1 through I
*
*      MAXEL(I) is the largest allowed number of electrons in
*      orbitals 1 through I
*
* Set up vertex weights W
* Set up arc weights    Y
*
* Reverse lexical ordering is used with
* weights of unoccupied orbitals set to 0
*
* Jeppe Olsen
*
       IMPLICIT REAL*8(A-H,O-Z)
       INTEGER W(NORB+1,NEL+1)
       INTEGER Y(NORB,NEL)
       INTEGER MAXEL(NORB),MINEL(NORB)
*
C      NTEST = 0
       CALL ISETVC(W,0,(NEL+1)*(NORB+1) )
       CALL ISETVC(Y,0,NEL*NORB)
*
*================
*  Vertex weights
*================
*
*. (Weight for vertex(IEL,IORB) is stored in W(IORB+1,IEL+1) )
      W(1,1) = 1
      DO 300 IEL = 0, NEL
        DO 200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            IF( IEL .GT. 0 ) THEN
              W(IORB+1,IEL+1) = W(IORB-1+1,IEL+1)
     &                        + W(IORB-1+1,IEL-1+1)
            ELSE
              W(IORB+1,1) = W(IORB-1+1,1)
            END IF
          END IF
  200   CONTINUE
  300 CONTINUE
*
*=============
* Arc weights
*=============
*
*. Weight for arc connecting vertices (IORB-1,IEL-1) and(IORB,IEL)
*. is stored in Y(IORB,IEL)
*. Y(IORB,IEL) = W(IORB-1,IEL)
      DO 1300 IEL = 1, NEL
        DO 1200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            Y(IORB,IEL) = W(IORB-1+1,IEL+1)
          END IF
 1200   CONTINUE
 1300 CONTINUE
*
      IF( NTEST .GE.10 ) THEN
C       WRITE(6,'(A)') ' Matrix of vertex weights '
C       WRITE(6,'(A)') ' ========================'
C       CALL IWRTMA(W,NORB+1,NEL+1,NORB+1,NEL+1)
        WRITE(6,'(A)') '  Matrix for arc weights  '
        WRITE(6,'(A)') '  ======================'
        CALL IWRTMA(Y,NORB,NEL,NORB,NEL)
      END IF
*
      RETURN
      END
