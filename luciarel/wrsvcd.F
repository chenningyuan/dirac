      SUBROUTINE WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM)
*
* Write scattered vector to disc
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION IPLAC(*),VAL(*)
*.Scratch
      DIMENSION VEC1(*)
*
      IF(LBLK.GT.0) THEN
*. Write the vector without markers in one block
         CALL SETVEC(VEC1,0.0D0,NDIM)
         DO 100 IEFF = 1, NSCAT
           VEC1(IPLAC(IEFF)) = VAL(IEFF)
  100    CONTINUE
         CALL TODSC(VEC1,NDIM,-1,LU)
      ELSE
*. Write the vector with the block format of file LUFORM
        CALL REWINO(LUFORM,LBLK)
        IBOT = 1
*. Loop over records
 1000   CONTINUE
        CALL IFRMDS(LBL,1,-1,LUFORM)
        CALL ITODS (LBL,1,-1,LU)
        IF(LBL.GE.0) THEN
C?        write(6,*) ' IBOT,IBOT + LBL -1 ',IBOT,
C?   &    IBOT + LBL -1
          CALL SETVEC(VEC1,0.0D0,LBL)
          DO 200 IEFF = 1, NSCAT
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
     &      VEC1(IPLAC(IEFF)-IBOT+1) = VAL(IEFF)
C?          IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
C?   &      write(6,*) ' Catch : IPLAC(IEFF) VAL(IEFF) ',
C?   &      IPLAC(IEFF),VAL(IEFF)
  200     CONTINUE
          CALL TODSC(VEC1,LBL,-1,LU)
*.Skip record on LUFORM
          CALL FRMDSC(VEC1,LBL,-1,LUFORM,IMZERO,IAMPACK)
          IBOT = IBOT + LBL
          GOTO 1000
        END IF
*
C?      CALL WRTVCD(VEC1,LU,1,LBLK)
      END IF
*
      RETURN
      END
