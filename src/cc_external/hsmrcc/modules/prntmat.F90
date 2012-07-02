!     THIS ROUTINE IS ADOPTED FROM GAMESS

      SUBROUTINE PRNTMAT(MATNAME,UNITNO,M,N,V,LDV,ROWIND,COLIND)

      implicit none
      CHARACTER*(*) MATNAME
      INTEGER UNITNO
      INTEGER M, N, LDV
      DOUBLE PRECISION V(LDV,N)
      INTEGER ROWIND(M)
      INTEGER COLIND(N)

      INTEGER IU
      INTEGER IMIN, IMAX
      INTEGER MAXCOL
      INTEGER I, J
      CHARACTER*104 MATN
      INTEGER LENNM, LR, LB

      IU = UNITNO
      IF (IU.LE.0) IU = 6

      IF ((LDV.LT.M).OR.(M.LE.0).OR.(N.LE.0)) THEN
         WRITE(IU,8001)
         WRITE(IU,*) (MATNAME(I:I),I=1,LEN(MATNAME))
         WRITE(IU,8003)
         RETURN
      ENDIF

 8001 FORMAT("INCORRECT PARAMETERS FOR PRNTMAT")
 8003 FORMAT("THE MATRIX IS NOT PRINTED")

      LENNM = len_trim(MATNAME)
      if (LENNM.gt.84) LENNM = 84

      LR = 104 - LENNM
      LB = LR/2

      DO I = 1, 104
         MATN(I:I) = ' '
      ENDDO

      DO I = 1, LENNM
         MATN(LB+I:LB+I) = MATNAME(I:I)
      ENDDO

      WRITE(IU,9001)
      WRITE(IU,9002)
      WRITE(IU,9003) (MATN(I:I),I=1,104)
      WRITE(IU,9002)
      WRITE(IU,9001)

 9001 FORMAT(1X)
 9002 FORMAT(6X,104('@'))
 9003 FORMAT(6X,104(A1))

      MAXCOL = 8

      IMAX = 0
  100 IMIN = IMAX + 1
      IMAX = IMAX + MAXCOL
      IF (IMAX .GT. N) IMAX = N                          

      WRITE (IU,9031)
      WRITE (IU,9032) (COLIND(I),I = IMIN,IMAX)
      WRITE (IU,9031)
      DO J = 1, M
         WRITE (IU,9033) ROWIND(J),(V(J,I),I = IMIN,IMAX)
      ENDDO
      IF (IMAX .LT. N) GO TO 100
      WRITE (IU,9031)                                                   

 9031 FORMAT(1X)
 9032 FORMAT(6X,8(4X,I5,4X))
 9033 FORMAT(I5,1X,8(F12.7,1X)) 

      CALL FLUSH(IU)

      RETURN                                                            

      END                                                               
