      SUBROUTINE CEIGEN(HR,HI,NDIM,EIGVAL,EIGVEC,SCR)
*
* A hermitian matrix is given in the form of a real part HR
* and an imaginary part HI.
*
* Find the eigenvalues and eigenvectors
* using a routine for real symmetric problems
*
* Output
* =======
*
* EIGVAL : Eigenvalues
* EIGVEC : Eigenvectors:EIGVEC(1,IVEC)     -EIGVEC(NDIM,IVEC) is real coefs
*                       EIGVEC(NDIM+1,IVEC)-EIGVEC(2*NDIM,IVEC) is imag coefs
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION HR(NDIM*(NDIM+1)/2),HI(NDIM*(NDIM+1)/2)
*. Output
      DIMENSION EIGVAL(2*NDIM)
      DIMENSION EIGVEC(2*NDIM,2*NDIM)
*. Note that EIGVEC is dimensioned as having 2*NDIM columns !
*. and that EIGVAL has 2NDIM elements
*. Scratch
      DIMENSION SCR(*)
*. Required length of SCR : (2*NDIM+1)*2*NDIM/2
      NTEST = 00
      IF(NTEST.GT.00) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' Welcome to CEIGEN '
        WRITE(6,*) ' =================='
        WRITE(6,*)
       END IF
*
*                                                   (Hr    -Hi)
* 1 : Assemble lower part of the symmetric matrix   (         )
*                                                   (Hi     Hr)
      CALL COPVEC(HR,SCR(1),NDIM*(NDIM+1)/2)
      DO IROW = NDIM+1,2*NDIM
        IROWOF = (IROW-1)*IROW/2+1
        IROWI = IROW-NDIM
        IROWIOF = IROWI*(IROWI-1)/2+1
*
        DO ICOL = 1, IROWI
          SCR(IROWOF-1+ICOL) = HI(IROWIOF-1+ICOL)
        END DO
*
        DO ICOL = IROWI+1,NDIM
          SCR(IROWOF-1+ICOL) = -HI(ICOL*(ICOL-1)/2+IROWI)
        END DO
*
        CALL COPVEC(HR(IROWI*(IROWI-1)/2+1),SCR(IROWOF+NDIM),IROWI)
      END DO
*
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Expanded complex matrix '
         WRITE(6,*) ' ======================= '
         CALL PRSYM(SCR,2*NDIM)
      END IF
*
* 2 : Diagonalize expanded matrix, and extract unique solutions
*
      CALL EIGEN(SCR,EIGVEC,2*NDIM,0,1)
      DO I = 1, 2*NDIM
        EIGVAL(I) = SCR(I*(I+1)/2)
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     & ' Eigenvalues and eigenvectors of expanded matrix'
        WRITE(6,*)
        CALL WRTMAT(EIGVAL,1,2*NDIM,1,2*NDIM)
CTF     do IVAL=1,2*NDIM,1
CTF       write(6,*) IVAL,EIGVAL(IVAL)
CTF     end do
        WRITE(6,*)
        CALL WRTMAT(EIGVEC,2*NDIM,2*NDIM,2*NDIM,2*NDIM)
      END IF
*
*. Each eigensolution occurs twice, eliminate one of each pair.
*. A bit of caution is required to distinguish between this
*. additional degeneracy and true degeneracies.
*  note that THRES is a sensitive thing which may have to be adjusted
*  put ntest to 10 to see what is going on
*
      ONE = 1.0D0
      ZERO = 0.0D0
      THRES = 1.0D-7
      IOFF = 1
 1000 CONTINUE
*. Find next set of degenerate eigenvalues
        EVAL = EIGVAL(IOFF)
        LDEG = 1
        DO I = IOFF+1,2*NDIM
          IF(ABS(EIGVAL(I)-EIGVAL(IOFF)).LE.THRES) THEN
            LDEG = LDEG + 1
          END IF
        END DO
*
        LDEG = 2
        IF(LDEG.EQ.2) THEN
*. No problem, just the additional degenerency is present,
*. zap eigenvector IOFF + 1
          ZERO = 0.0D0
          CALL SETVEC(EIGVEC(1,IOFF+1),ZERO,2*NDIM)
        ELSE
*. more than 2 degenerate eigenvectors, eliminate by Gram orthogonalization
*. (Gram-Schmidt)
          DO IVEC = IOFF,IOFF+LDEG-1
*. | Vec I > = Scalar * ( | Vec I > - Sum(J<I) <Vec J | Vec I > | Vec J > )
*. Orthogonalize
            DO JVEC = IOFF, IVEC-1
*. < Vec J | Vec I >
               XR = INPROD(EIGVEC(1,JVEC),EIGVEC(1,IVEC),NDIM)
     &            + INPROD(EIGVEC(1+NDIM,JVEC),
     &              EIGVEC(1+NDIM,IVEC),NDIM)
               XI = INPROD(EIGVEC(1,JVEC),
     &              EIGVEC(1+NDIM,IVEC),NDIM)
     &            - INPROD(EIGVEC(1+NDIM,JVEC),
     &              EIGVEC(1,IVEC),NDIM)
              CALL VECSUMC(EIGVEC(1,IVEC),EIGVEC(1,IVEC),
     &                     EIGVEC(1,JVEC),
     &                     ONE,ZERO,-XR,-XI,NDIM)
             END DO
*. Normalize
             XNORM = INPROD(EIGVEC(1,IVEC),EIGVEC(1,IVEC),
     &                      2*NDIM)
        if(ntest.ge.10) WRITE(6,*) ' ldeg XNORM = ',ldeg, XNORM
             IF(ABS(XNORM).LE.THRES) THEN
*. This is a zero vector
               CALL SETVEC(EIGVEC(1,IVEC),ZERO,2*NDIM)
             ELSE
*. This is not a zero vector
               SCALE = ONE/SQRT(XNORM)
               CALL SCALVE(EIGVEC(1,IVEC),SCALE,2*NDIM)
             END IF
           END DO
*. End of not so simple case
         END IF
*. End of this batch of degenerate eigenvalues
         IOFF = IOFF + LDEG
       IF(IOFF.LT.2*NDIM) GOTO 1000
*
*
* Only unique eigenvectors have survived, collect these
*
      IEFF = 0
      DO IVEC = 1, 2*NDIM
        XNORM = INPROD(EIGVEC(1,IVEC),
     &                 EIGVEC(1,IVEC),2*NDIM)
        if(ntest.ge.10) WRITE(6,*) ' XNORM = ',XNORM
        IF(XNORM.NE.0.0D0) THEN
          IEFF = IEFF + 1
          EIGVAL(IEFF) = EIGVAL(IVEC)
          CALL COPVEC(EIGVEC(1,IVEC),EIGVEC(1,IEFF),
     &                2*NDIM)
        END IF
      END DO
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     & ' Unique Eigenvalues and eigenvectors  matrix'
        WRITE(6,*)
        CALL WRTMAT(EIGVAL,1,2*NDIM,1,2*NDIM)
        WRITE(6,*)
        CALL WRTMAT(EIGVEC,2*NDIM,2*NDIM,2*NDIM,2*NDIM)
      END IF
      IF(IEFF.NE.NDIM) THEN
        WRITE(6,*) ' CEIGEN in trouble again '
        WRITE(6,*) ' Obtained number of unique eigenpairs',IEFF
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Unique eigenvalues and eigenvectors'
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        CALL WRTMAT(EIGVAL,1,NDIM,1,NDIM)
        WRITE(6,*)
        CALL WRTMAT(EIGVEC,2*NDIM,NDIM,2*NDIM,2*NDIM)
      END IF
*
      RETURN
      END
