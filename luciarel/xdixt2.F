      SUBROUTINE XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
*
* Obtain XDX = X * (DIA+Shift)-1 * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION XDX(NXRDM,NXRDM)
      DIMENSION X(NXRDM,NXCDM),DIA(NXCDM)
      DIMENSION SCR(NXCDM)
*
      CALL SETVEC(XDX,0.0D0,NXRDM ** 2 )
      THRES = 1.0D-9
      DO 100 J=1,NXRDM
*.Scr(k) = X(J,K)/(SHIFT+DIA(K)
        DO 10 K = 1, NXCDM
          IF(ABS(SHIFT+DIA(K)) .GT. THRES ) THEN
            SCR(K) = X(J,K)/(SHIFT+DIA(K))
          ELSE
            SCR(K) = X(J,K)/THRES
          END IF
   10   CONTINUE
*
        DO 20 K = 1, NXCDM
          FACTOR = SCR(K)
          CALL VECSUM(XDX(1,J),XDX(1,J),X(1,K),1.0D0,FACTOR,NXRDM)
   20   CONTINUE
*
  100 CONTINUE
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' x (Dia + shift ) - 1 X(T) from XDIXT2 '
        CALL WRTMAT(XDX,NXRDM,NXRDM,NXRDM,NXRDM)
      END IF
*
      RETURN
      END
