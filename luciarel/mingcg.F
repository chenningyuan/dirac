      SUBROUTINE MINGCG(MV8,LU1,LU2,LU3,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPROD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CONVER = .FALSE.
      ITER = 1
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X
      CALL REWINO(LU1)
      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINO(LU2)
      CALL FRMDSC(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
      CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-1.0D0,NVAR)
*
      RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
      ERROR(1) = RNORM
      CALL REWINO(LU2)
      CALL TODSC(VEC1,NVAR,-1  ,LU2)
*. Preconditioner H times initial vector , H * R
*.H * R
      CALL REWINO(LUDIA)
      CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
      CALL DIAVC2(VEC2,VEC1,VEC2,W,NVAR)
      IF(LUPROJ.NE.0) THEN
        CALL REWINO(LUPROJ)
        CALL FRMDSC(VEC1,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL REWINO(LU2)
        CALL FRMDSC(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPROD(VEC1,VEC2,NVAR)
*. P = RHO * H*R
      RHO = 1.0D0
      CALL SCALVE(VEC2,RHO,NVAR)
      CALL REWINO(LU3)
      CALL TODSC(VEC2,NVAR,-1  ,LU3)
      CALL COPVEC(VEC2,VEC1,NVAR)
*.S = AP
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINO (LU3)
      CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*.    P is assumed in VEC1 and S = A*P in VEC2

        NITER = NITER + 1
       IF ( NTEST .GE. 10 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPROD(VEC1,VEC2,NVAR)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S
        CALL REWINO(LU2)
        CALL FRMDSC(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-A,NVAR)
        CALL REWINO(LU2)
        CALL TODSC(VEC1,NVAR,-1  ,LU2)
*.    new residual has been obtained , check for convergence
        RNORM = INPROD(VEC1,VEC1,NVAR)
        ERROR(ITER+1) = SQRT(RNORM)
*.    X = X + A * P
        CALL REWINO(LU1)
        CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        CALL REWINO(LU3)
        CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC2,VEC1,1.0D0,A,NVAR)
        CALL REWINO(LU1)
        CALL TODSC(VEC1,NVAR,-1  ,LU1)
*
        IF( SQRT(RNORM) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.       H * R
           CALL REWINO(LU2)
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           CALL REWINO(LUDIA)
           CALL FRMDSC(VEC1,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL DIAVC2(VEC1,VEC2,VEC1 ,W,NVAR)
           IF(LUPROJ.NE.0) THEN
             CALL REWINO(LUPROJ)
             CALL FRMDSC(VEC2,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
             CALL REWINO(LU2)
             CALL FRMDSC(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
           END IF
           GAMMA = INPROD(VEC1,VEC2,NVAR)
           B = GAMMA/C
*.       P = RHO*(H*R + B*P)
           CALL REWINO(LU3)
           CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,B,NVAR)
*.       Define next RHO
           RHO = 1.0D0
           CALL SCALVE(VEC1,RHO,NVAR)
           CALL REWINO(LU3)
           CALL TODSC(VEC1,NVAR,-1  ,LU3)
*.       S = MATRIX * P
           CALL MV8(VEC1,VEC2,0,0)
           CALL REWINO(LU3)
           CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
*.End of prepations for next iteration
        END IF
*
*
 1000 CONTINUE
 1001 CONTINUE
      IF(NTEST .GT. 0 ) THEN
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(1H0,'  convergence was obtained in...',I3,' iterations',/,
     +        1H ,'  norm of residual..............',F13.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT +1 )
 1020  FORMAT(1H0,' convergence was not obtained in',I3,'iterations',/,
     +        1H ,' norm of residual...............',F13.8)
      END IF
      END IF
C
      IF(NTEST.GT. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(1H0,' solution to set of linear equations')
       CALL REWINO(LU1)
       CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
       CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call wrtmat(vec2,1,nvar,1,nvar)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(1H0,10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(1H ,12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
C
      RETURN
      END
