      SUBROUTINE MINDV4(MV7,
     &                  VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &                  IOLSEN,IPICO)
*
* Davidson algorithm , requires two vectors in core
* Multi root version
*
* Allows updating of preconditioning matrix so this is
* the current eigenvector approximation
* is an eigenvector for the preconditioner
*
* Jeppe Olsen Sept 89
*             Jan  92 : MV7 entry
*
* Input :
* =======
*        MV7 : Name of routine performing matrix*vector calculation
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    complete vector
*        LU2,LU3   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
* On input LU1 is supposed to hold initial guess to eigenvectors
*
* IOLSEN : Use inverse iteration modified Davidson
* IPICO  : Use perturbation estimate of new vector instead of
*          variational method
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
*
       EXTERNAL MV7
*
       IOLSTM = IOLSEN
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Variational'
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Perturbational'
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Normal Davidson, Variational '
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Normal Davidson, Perturbational'
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' SORRY MINDV2 WOUNDED , MAXVEC .LT. 2*NROOT '
         Call Abend2( ' ENFORCED STOP IN MINDV2' )
       END IF
*
       IF(IPICO.NE.0) THEN
         MAXVEC = 2*NROOT
       END IF
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-7
       CONVER = .FALSE.
       DO 1234 MACRO = 1,1
*
*.   INITIAL ITERATION
       ITER = 1
       CALL REWINO( LU1 )
       CALL REWINO( LU2 )
       DO 10 IVEC = 1,NINVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL MV7(VEC1,VEC2,0,0)
         CALL TODSC(VEC2,NVAR,-1  ,LU2)
*        PROJECTED MATRIX
         CALL REWINO( LU1)
         DO 8 JVEC = 1, IVEC
           CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*  DIAGONALIZE INITIAL PROJECTED MATRIX
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF( IPRT  .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Initial set of eigenvalues '
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF),IROOT=1,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
** LOOP OVER ITERATIONS
*
 1000 CONTINUE
      IF(IPRT  .GE. 10 ) THEN
       WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
      END IF


        ITER = ITER + 1
*
** 1          NEW DIRECTION TO BE INCLUDED
*
*   1.1 : R = H*X - EIGAPR*X
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         CALL SETVEC(VEC1,0.0D0,NVAR)
*
         CALL REWINO( LU2)
         DO 60 IVEC = 1, NVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   60    CONTINUE
*
         EIGAPR = EIG(ITER-1,IROOT)
         CALL REWINO( LU1)
         DO 50 IVEC = 1, NVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   50    CONTINUE
           IF ( IPRT  .GE.600 ) THEN
             WRITE(6,*) '  ( HX - EX ) '
             CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
         RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
         RNRM(ITER-1,IROOT) = RNORM
         IF(RNORM.LT. TEST ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
*.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
           CALL REWINO( LUDIA)
           CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL H0M1TV(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
              CALL REWINO(LU3)
              CALL TODSC(VEC1,NVAR,-1,LU3)
* Current eigen vector
              CALL REWINO( LU1)
              CALL SETVEC(VEC1,0.0D0,NVAR)
              DO 59 IVEC = 1, NVEC
                CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
                FACTOR = AVEC((IROOT-1)*NVEC+ IVEC)
                CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   59         CONTINUE
              IF ( IPRT  .GE. 600 ) THEN
                WRITE(6,*) ' And X  '
                CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
              END IF
              CALL TODSC(VEC1,NVAR,-1,LU3)
* (H0 - E )-1  * X
              CALL REWINO( LUDIA)
              CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
              CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
     &                   H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
              CALL TODSC(VEC2,NVAR,-1,LU3)
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPROD(VEC2,VEC1,NVAR)
* is X an eigen vector for (H0 - 1 ) - 1
              CALL VECSUM(VEC2,VEC1,VEC2,GAMMA,-1.0D0,NVAR)
              VNORM = SQRT(MAX(0.0D0,INPROD(VEC2,VEC2,NVAR)))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                CALL REWINO(LU3)
                CALL FRMDSC(VEC2,NVAR,-1,LU3,IMZERO,IAMPACK)
                DELTA = INPROD(VEC1,VEC2,NVAR)
                CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSUM(VEC1,VEC1,VEC2,FACTOR,1.0D0,NVAR)
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) '  Modified new trial vector '
                  CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
                END IF
              ELSE
                IF(IPRT.GT.0) WRITE(6,*)
     &          ' Inverse correction switched of'
                CALL REWINO(LU3)
                CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
              END IF
            END IF
*. 1.3 ORTHOGONALIZE TO ALL PREVIOUS VECTORS
           XNRMI =    INPROD(VEC1,VEC1,NVAR)
           CALL REWINO( LU1 )

           DO 80 IVEC = 1,NVEC+IADD-1
             CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
   80      CONTINUE
*. 1.4 Normalize vector and check for linear dependency
           SCALE = INPROD(VEC1,VEC1,NVAR)
           IF(ABS(SCALE)/XNRMI .LT. 1.0D-10) THEN
*. Linear dependency
             IADD = IADD - 1
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(6,*) '  Trial vector linear dependent so OUT !!! '
             END IF
           ELSE
             C1NRM = SQRT(SCALE)
             FACTOR = 1.0D0/SQRT(SCALE)
             CALL SCALVE(VEC1,FACTOR,NVAR)
*
             CALL TODSC(VEC1,NVAR,-1  ,LU1)
             IF ( IPRT  .GE.600 ) THEN
               WRITE(6,*) 'ORTHONORMALIZED (D-E)-1 *( HX - EX ) '
               CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
             END IF
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
*
*  2.1: MULTIPLY NEW DIRECTION WITH MATRIX
       CALL REWINO( LU1)
       CALL REWINO( LU2)
       DO 110 IVEC = 1, NVEC
         CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
         CALL FRMDSC(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
  110  CONTINUE
*
      DO 150 IVEC = 1, IADD
        CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        CALL MV7(VEC1,VEC2,0,0)
        CALL TODSC(VEC2,NVAR,-1  ,LU2)
*   AUGMENT PROJECTED MATRIX
        CALL REWINO( LU1)
        DO 140 JVEC = 1, NVEC+IVEC
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
  140   CONTINUE
  150 CONTINUE
*  DIAGONALIZE PROJECTED MATRIX
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
       END IF
*
      IF(IPRT .GE. 5) THEN
        WRITE(6,*) ' PROJECTED MATRIX AND EIGEN PAIRS '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  PERHAPS RESET OR ASSEMBLE CONVERGED EIGENVECTORS
*
  901 CONTINUE
*
      IPULAY = 1
      IF(IPULAY.EQ.1 .AND. MAXVEC.EQ.3 .AND.NVEC.GE.2.
     &   .AND. .NOT.CONVER) THEN
* Save trial vectors : 1 -- current trial vector
*                      2 -- previous trial vector orthogonalized
        CALL REWINO( LU3)
        CALL REWINO( LU1)
*. Current trial vector
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2200 IVEC = 1, NVEC
          CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
         CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2200   CONTINUE
        SCALE = INPROD(VEC1,VEC1,NVAR)
        SCALE  = 1.0D0/SQRT(SCALE)
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC(VEC1,NVAR,-1  ,LU3)
* Previous trial vector orthonormalized
        CALL REWINO(LU1)
        CALL FRMDSC(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        SCALE2 = INPROD(VEC2,VEC2,NVAR)
        SCALE2 = 1.0D0/SQRT(SCALE2)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC(VEC2,NVAR,-1,LU3)
*
        CALL REWINO( LU1)
        CALL REWINO( LU3)
        DO 2411 IVEC = 1,2
          CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC1,NVAR,-1,  LU1)
 2411   CONTINUE
*. Corresponding sigma vectors
        CALL REWINO ( LU3)
        CALL REWINO( LU2)
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2250 IVEC = 1, NVEC
          CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
          CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2250   CONTINUE
*
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC(VEC1,NVAR,-1,  LU3)
* Sigma vector corresponding to second vector on LU1
        CALL REWINO(LU2)
        CALL FRMDSC(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC(VEC2,NVAR,-1,LU3)
*
        CALL REWINO( LU2)
        CALL REWINO( LU3)
        DO 2400 IVEC = 1,2
          CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC2,NVAR,-1  ,LU2)
 2400   CONTINUE
        NVEC = 2
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 2410 IROOT = 1,NVEC
          AVEC((IROOT-1)*NVEC+IROOT) = 1.0D0
 2410   CONTINUE
*.Projected hamiltonian
       CALL REWINO( LU1 )
       DO 2010 IVEC = 1,NVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL REWINO( LU2)
         DO 2008 JVEC = 1, IVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
 2008    CONTINUE
 2010  CONTINUE
      END IF
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINO( LU3)
        DO 320 IROOT = 1, NROOT
          CALL REWINO( LU1)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 200 IVEC = 1, NVEC
            CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  200     CONTINUE
*
          SCALE = INPROD(VEC1,VEC1,NVAR)
          SCALE  = 1.0D0/SQRT(SCALE)
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC(VEC1,NVAR,-1  ,LU3)
  320   CONTINUE
        CALL REWINO( LU1)
        CALL REWINO( LU3)
        DO 411 IVEC = 1,NROOT
          CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC1,NVAR,-1,  LU1)
  411   CONTINUE
* CORRESPONDING SIGMA VECTOR
        CALL REWINO ( LU3)
        DO 329 IROOT = 1, NROOT
          CALL REWINO( LU2)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 250 IVEC = 1, NVEC
            CALL FRMDSC(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  250     CONTINUE
*
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC(VEC1,NVAR,-1,  LU3)
  329   CONTINUE
* PLACE C IN LU1 AND HC IN LU2
        CALL REWINO( LU2)
        CALL REWINO( LU3)
        DO 400 IVEC = 1,NROOT
          CALL FRMDSC(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL TODSC (VEC2,NVAR,-1  ,LU2)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
C
      END IF
C
C     IF( ITER .LT. MAXIT .AND. .NOT. CONVER) GOTO 1000
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*. Place first eigenvector in vec1
      CALL REWINO(LU1)
      CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)



* ( End of loop over iterations )
*
*
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*. Final eigenvalues
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
*
      IF ( IPRT .GT. 1 ) THEN
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL REWINO( LU1)
            CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),
     &                 RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
 1234 CONTINUE
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
