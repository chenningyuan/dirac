      SUBROUTINE HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
*
* Multiply P1P2Q preconditioner with a vector
* Jeppe Olsen , July 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General Input
      DIMENSION PHP(*),PHQ(NP1,NQ),QHQ(NQ)
*. PHP is in lower triangular form
*. Specific input
      DIMENSION VECIN(*)
*.Output
      DIMENSION VECUT(*)
*
      NP = NP1+NP2
      NPQ = NP + NQ
*
      CALL SETVEC(VECUT,0.0D0,NPQ)
*. PHQ * VECIN
      CALL MATVCC(PHQ,VECIN(1+NP),VECUT,NP1,NQ,0)
*. PHP * VECIN
      IJ = 0
      DO 60 I = 1, NP
      DO 50 J = 1, I
        IJ = IJ + 1
        VECUT(I) = VECUT(I)+PHP(IJ)*VECIN(J)
        VECUT(J) = VECUT(J)+PHP(IJ)*VECIN(I)
   50 CONTINUE
      VECUT(I) = VECUT(I)-PHP(IJ)*VECIN(I)
   60 CONTINUE
* QHP * VECIN
      CALL MATVCC(PHQ,VECIN,VECUT(1+NP),NP1,NQ,1)
*. QHQ * VECIN
      DO 100 I = 1, NQ
        VECUT(I+NP) = VECUT(I+NP)+QHQ(I)*VECIN(NP+I)
  100 CONTINUE
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTV '
        CALL WRTMAT(VECIN,1,NPQ,1,NPQ)
        CALL WRTMAT(VECUT,1,NPQ,1,NPQ)
      END IF
      RETURN
      END
