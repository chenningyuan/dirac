      SUBROUTINE MICDV6(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version
*
* From MICDV4, Jeppe Olsen, April 1997
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
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
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
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
       REAL*8 INPRDD
       DIMENSION XJEP(1000)
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
       IPRT = 600
C      IPRT = IPRTXX
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV6 '
         Call Abend1( 20 )
       END IF
*
       WRITE(6,*) ' LU1 LU2 LU3 LU4 LU5 LU6',
     &              LU1,LU2,LU3,LU4,LU5,LU6
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       DO IVEC = 1,NINVEC
         CALL REWINO(LU5)
         CALL REWINO(LU6)
         CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU5,LU6)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU6)
         CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
         DO JVEC = 1, IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
         END DO
*        ^ End of loop over JVEC
       END DO
*      ^ End of loop over IVEC
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(6,*) ' Initial set of eigen values '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
C     TRAVCD(VEC1,VEC2,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
       CALL REWINO(LU3)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINO(LU3)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)
*
       IF(IPRT.GE.600) THEN
         WRITE(6,*) ' Initial set of eigenvectors '
         CALL REWINO(LU1)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(6,*) ' Initial set of sigma vectors '
         CALL REWINO(LU2)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*
       NVEC = NROOT
       IF (MAXIT .EQ. 1 ) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        IF(IPRT  .GE. 10 ) THEN
         WRITE(6,*) ' Info from iteration .... ', ITER
        END IF
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       CALL REWINO(LU1)
       CALL REWINO(LU2)
*
       DO 100 IROOT = 1, NROOT
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINO(LU7)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
         CALL REWINO(LU7)
         CALL REWINO(LU5)
         FACTOR = - EIGAPR
         CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU2,LU7,LU5,0,LBLK)
*
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(6,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU5,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
         IF(IPRT.GE.1000) WRITE(6,*) ' RNORM = ', RNORM
         IF(RNORM.LT. TEST .OR.
     &      (ITER.GT.2.AND.
     &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU6
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(6,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU5)
           CALL REWINO(LU6)
           CALL H0M1TD(LU6,LUDIA,LU5,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
C     H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU6,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU5
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             GAMMA = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
             IF(IPRT.GE.1000) WRITE(6,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINO(LU5)
             CALL REWINO(LU7)
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU7,LU5,0,LBLK))
              IF(IPRT.GE.1000) write(6,*) ' VNORM = ', VNORM
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                FACTOR = -(DELTA/GAMMA)
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU6,LU5,LU7,1,LBLK)
                CALL COPVCD(LU7,LU6,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU6,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
*.. Vectors on LU1

           CALL REWINE( LU1 ,LBLK)
           DO IVEC = 1, NROOT
               CALL REWINO(LU6)
               WORK(IVEC) = -INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
           END DO
           IF(IPRT.GE.1000) THEN
             Write(6,*) ' Overlap with vectors on LU1'
             call wrtmat(work,1,nroot,1,nroot)
           END IF
           ONE = 1.0D0
           CALL REWINO(LU1)
           CALL MVCSMD(LU1,WORK,LU7,LU5,VEC1,VEC2,NROOT,1,LBLK)
           CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU6,LU5,1,LBLK)
           IF(IPRT.GE.1000) THEN
             write(6,*) ' orthogonalized to vectors on LU1'
             CALL WRTVCD(VEC1,LU5,1,LBLK)
           END IF
*.. Vectors on LU3
           IF(NVEC+IADD-1-NROOT.GT.0) THEN
             CALL REWINE( LU3 ,LBLK)
             DO IVEC = 1, NVEC+IADD-1-NROOT
                 CALL REWINO(LU6)
                 WORK(IVEC) = -INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
             END DO
             ONE = 1.0D0
             CALL REWINO(LU3)
             CALL MVCSMD(LU3,WORK,LU7,LU6,VEC1,VEC2,NVEC+IADD-1-NROOT,
     &            1,LBLK)
             CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU5,LU6,1,LBLK)
           ELSE
             CALL REWINO(LU3)
             CALL COPVCD(LU5,LU6,VEC1,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU6,LBLK)
           CALL SCLVCD(LU6,LU3,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU6,LU7,FACTOR,VEC1,1,LBLK)
             WRITE(6,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU7,1,LBLK)
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
* ====================================================
*  2 : Optimal combination of new and old directions
* ====================================================
*
*  2.1: Multiply new directions with matrix
*
      IF(IPRT.GE.1000) THEN
        WRITE(6,*) ' Vectors on LU3'
        WRITE(6,*) 'NVEC-NROOT+IADD =',NVEC-NROOT+IADD
        CALL REWINO(LU3)
        DO IVEC = 1, NVEC-NROOT+IADD
          CALL WRTVCD(VEC1,LU3,0,LBLK)
        END DO
      END IF
*
      CALL SKPVCD(LU3,NVEC-NROOT,VEC1,1,LBLK)
      CALL SKPVCD(LU4,NVEC-NROOT,VEC1,1,LBLK)
      DO IVEC = 1, IADD
        CALL REWINE(LU5,LBLK)
        CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
        CALL MV7(VEC1,VEC2,LU5,LU6)
        CALL REWINE(LU6,LBLK)
        CALL COPVCD(LU6,LU4,VEC1,0,LBLK)
*.2.2 Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO JVEC = 1, NROOT
          CALL REWINE(LU6,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
        END DO
*
        CALL REWINE(LU3,LBLK)
C       DO JVEC = NROOT+1, NVEC+IADD
        DO JVEC = NROOT+1, NVEC+IVEC
         CALL REWINE(LU6,LBLK)
         IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
         APROJ(IJ) = INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
        END DO
      END DO
*     /\ End do over new trial vectors
*. 2.3 Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
*
      DO IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
      END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(6,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset
*
*
* case 1 : Only NROOT vectors can be stored
*          save current eigenvector approximations
* Case 2 : Atleast 2*NROOT eigenvectors can be saved
*          Current eigenvactor approximations+
*          vectors allowing generation of previous approxs.
*
*
C     IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
*
*. 3.1 : Current wave function approximations, collect on LU7
*
        IF(IPRT.GE.1000) THEN
        write(6,*) ' I am going to reset '
        write(6,*) ' nroot, nvec ', nroot,nvec
        END IF
        IF(NVEC.GT.NROOT) THEN
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU1 to LU5
          CALL REWINO(LU1)
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of c reset, part 1'
*. and add LU3 to LU5
          CALL REWINO(LU3)
          ONE = 1.0D0
          CALL MVCSMD2(LU3,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of c reset, part 2'
C              MVCSMD2(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
C    &         VEC1,VEC2,NVEC,IREW,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
          CALL REWINE(LU5,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
*. scale LU5 => LU7
          CALL REWINO(LU5)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of c reset, part 3'
        END DO
*. Transfer C vectors from LU7 to LU1
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU1,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU1,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 4'
*. 3.2 : corresponding sigma vectors
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU2 to LU5
          CALL REWINO(LU2)
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 1'
*. and add LU4 to LU5
          CALL REWINO(LU4)
          CALL MVCSMD2(LU4,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 2'
          SCALE  = WORK(IROOT)
*. scale LU5 => LU7
          CALL REWINO(LU5)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 3'
        END DO
*. Transfer sigma  vectors from LU7 to LU2
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU2,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU2,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of s reset, part 4'
        NNVEC = NROOT
        IF(3*NROOT.LE.MAXVEC) THEN
*
*. Orthogonalize the
*. last set of correction vectors to the current
*. eigenvectors on LU1, and save on LU2
*. Overlap with root approximations
*. Start of last set of trial vectors
          ISTART = NVEC-NROOT-IADD+1
          CALL SKPVCD(LU3,ISTART-1,VEC1,1,LBLK)
*
          CALL REWINO(LU5)
          DO JVEC = 1, IADD
*. Orthogonalize to vectors on LU1
            CALL REWINO(LU7)
            CALL COPVCD(LU3,LU7,VEC1,0,LBLK)
            IF(IPRT.GE.1000) THEN
            write(6,*) ' Initial vector on LU7'
            CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
            CALL REWINO(LU1)
            DO IROOT = 1, NROOT
              CALL REWINO(LU7)
              WORK(IROOT+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU1,0,LBLK)
            END DO
            IF(IPRT.GE.1000)
     &      write(6,*) ' c, part1 finito'
*. And to trial vectors on LU5
            CALL REWINO(LU5)
            DO KVEC = 1, JVEC-1
              CALL REWINO(LU7)
              WORK(NROOT+KVEC+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU5,0,LBLK)
            END DO
            WORK(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
            IF(IPRT.GE.1000)
     &      write(6,*) ' c, part2 finito'
*
            ONE = 1.0D0
            CALL MVCSMD2(LU1,WORK(1+(JVEC-1)*2*NROOT),ONE ,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            ONE = 1.0D0
            CALL MVCSMD2(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000) THEN
              write(6,*) ' c, part4 finito'
              write(6,*) ' Vector after sec ort'
              CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
*
            FACTOR = INPRDD(VEC1,VEC1,LU7,LU7,1,LBLK)
            SCALE = 1.0D0/SQRT(FACTOR)
            CALL SCALVE(WORK((JVEC-1)*2*NROOT+1),SCALE,
     &           NROOT+JVEC)
            CALL REWINO(LU7)
            CALL SCLVCD(LU7,LU5,SCALE,VEC1,0,LBLK)
C                SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
           IF(IPRT.GE.1000)
     &       write(6,*) ' c, part6 finito'
         END DO
*        /\ End of loop over orthogonalized directions
*. Transfer modified directions to LU3
         CALL REWINO(LU3)
         CALL REWINO(LU5)
         DO JVEC =1, IADD
           CALL COPVCD(LU5,LU3,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000) THEN
           write(6,*) ' c, part7 finito'
           WRITE(6,*) ' Additional trial vectors on LU3'
           CALL REWINO(LU3)
           DO JVEC = 1, IADD
            CALL WRTVCD(VEC1,LU3,0,LBLK)
           END DO
         END IF
* Sigma vectors corresponding to orthogonalized directions
         CALL SKPVCD(LU4,ISTART-1,VEC1,1,LBLK)
         CALL REWINO(LU5)
         DO JVEC = 1, IADD
            CALL REWINO(LU7)
            CALL COPVCD(LU4,LU7,VEC1,0,LBLK)
*
            FACT = WORK(NROOT+JVEC+(JVEC-1)*2*NROOT)
            CALL MVCSMD2(LU2,WORK(1+(JVEC-1)*2*NROOT),FACT,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 1 finito '
            ONE = 1.0D0
            CALL MVCSMD2(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 2 finito '
            CALL REWINO(LU7)
            CALL COPVCD(LU7,LU5,VEC1,0,LBLK)
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 3 finito '
         END DO
*        /\ End of loop over orthogonalized directions
*. Copy back to LU4
         CALL REWINO(LU4)
         CALL REWINO(LU5)
         DO JVEC = 1, IADD
            CALL COPVCD(LU5,LU4,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000)
     &   write(6,*) ' s, part 4 finito '
         NNVEC = NROOT + IADD
       END IF
       IF(IPRT.GE.1000) THEN
         WRITE(6,*) ' Additional sigma vectors on LU4'
         CALL REWINO(LU4)
         DO JVEC = 1, IADD
           CALL WRTVCD(VEC1,LU4,0,LBLK)
         END DO
       END IF
*
*      /\ End if more than NROOT vectors to be reset
       NVEC = NNVEC
       END IF
*.     ^ End if reset
*. New subspace
*. Calculate subspace Hamiltonian from actual vectors
*. on disc
       IF(IPRT.GE.1000) write(6,*) ' Subspace hamiltonian'
       CALL REWINO(LU1)
       CALL REWINO(LU3)
       DO IVEC = 1, NVEC
*
         CALL REWINO(LU5)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINO(LU2)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
         END DO
         CALL REWINO(LU4)
         DO JVEC = NROOT+1,IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
         END DO
       END DO
       write(6,*) ' Reset hamiltonian'
       call prsym(aproj,nvec)
*. Test : Orthogonality of new vectors
       CALL REWINO(LU1)
       CALL REWINO(LU3)
       DO IVEC = 1, NVEC
*
         CALL REWINO(LU5)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINO(LU1)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU1,0,LBLK)
         END DO
         CALL REWINO(LU3)
         DO JVEC = NROOT+1,IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU3,0,LBLK)
         END DO
       END DO
       write(6,*) ' Overlap matrix    '
       call prsym(xjep,nvec)
*. End of resetting business
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE

* ( End of loop over iterations )
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
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
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
      ELSE
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
