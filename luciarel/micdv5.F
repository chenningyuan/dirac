      SUBROUTINE MICDV5(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  THRES_E)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Special version for NROOT = 1, MAXVEC = 2 !!
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
*   THRES_E   : Convergence criteria for eigenvalue
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
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
       REAL*8 INPRDD
*
       IF(NROOT.NE.1) THEN
         WRITE(6,*)
     &   ' MICDV5 is wrong path, since NROOT .ne. 1, NROOT =',
     &     NROOT
         Call Abend1( 20 )
       END IF
       IF(MAXVEC.NE.2) THEN
         WRITE(6,*)
     &   ' MICDV5 is wrong path, since MAXVEC.ne. 2, MAXVEC=',
     &     MAXVEC
         Call Abend1( 20 )
       END IF
       IF(NINVEC.NE.1) THEN
         WRITE(6,*) ' Several input vectors ', NINVEC
         WRITE(6,*) ' Only the first will be used '
       END IF
*
       IPICO = 0
       IOLSTM = 1
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
*

       IPRT = 3
*
       IF(IPRT.GE.1) THEN
         WRITE(6,*) ' Convergence threshold for eigenvalues', THRES_E
       END IF
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
       IROOT = 1
*
* ===================
*.Initial iteration
* ===================
*
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       CALL MV7(VEC1,VEC2,LU1,LU2)
*
       IF(IPRT.GE.600) THEN
         WRITE(6,*) ' Initial sigma-vector'
         CALL WRTVCD(VEC1,LU2,1,LBLK)
       END IF
*. Projected matrix
       APROJ(1) = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,1     )
       END IF
*. Diagonalize initial projected matrix : Dimension one :simple
       EIG(1,IROOT) = APROJ(1)
       AVEC(1) = 1.0D0
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       NVEC = 1
       ITERP = 1
*
* ======================
*. Loop over iterations
* ======================
*
      DO 1000 ITER = 2, MAXIT+1
       ITERP = ITER -1
       IF(IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Info from iteration .... ', ITER
       END IF
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IROOT = 1
*
       EIGAPR = EIG(ITER-1,IROOT)
       FACHC = 1.0D0
       FACC  = -EIGAPR
       CALL VECSMD(VEC1,VEC2,FACC,FACHC,LU1,LU2,LU4,1,LBLK)
C           VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
       IF ( IPRT  .GE. 600 ) THEN
         WRITE(6,*) '  ( HX - EX ) '
         CALL WRTVCD(VEC1,LU4,1,LBLK)
       END IF
*  Strange place to put convergence but ....
       RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
       RNRM(ITER-1,IROOT) = RNORM
       IF(RNORM.LT. TEST ) THEN
          CONVER = .TRUE.
          RTCNV(IROOT) = .TRUE.
       ELSE
          RTCNV(IROOT) = .FALSE.
          CONVER = .FALSE.
       END IF
       IF(ITER.GT.2.AND.
     & EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT).LT.THRES_E) CONVER = .TRUE.
       IF( ITER .EQ. MAXIT+1 .OR. CONVER ) GOTO 1001
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU3
       IF( .NOT. RTCNV(IROOT) ) THEN
         IADD = IADD + 1
         CALL REWINO(LUDIA)
         CALL REWINO(LU3)
         CALL REWINO(LU4)
         CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &               H0,-EIGAPR,H0SCR,XH0IX,
     &               NP1,NP2,NQ,VEC1,VEC2,IPRT)
         IF ( IPRT  .GE. 600) THEN
           WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
           CALL WRTVCD(VEC1,LU3,1,LBLK)
         END IF
*
         IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU4
           CALL REWINO(LU1)
           CALL REWINO(LU4)
           CALL REWINO(LUDIA)
*
           CALL H0M1TD(LU4,LUDIA,LU1,LBLK,Np1+Np2+NQ,
     &                 IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)

* Gamma = X(T) * (H0 - E) ** -1 * X
           GAMMA = INPRDD(VEC1,VEC2,LU1,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
           VNORM =
     &     SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU1,LU4,1,LBLK))
           IF(VNORM .GT. 1.0D-7 ) THEN
             IOLSAC = 1
           ELSE
             IOLSAC = 0
           END IF
           IF(IOLSAC .EQ. 1 ) THEN
             IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
             DELTA = INPRDD(VEC1,VEC2,LU1,LU3,1,LBLK)
             FACTOR = -DELTA/GAMMA
             IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
             IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
             CALL COPVCD(LU5,LU3,VEC1,1,LBLK)

             IF(IPRT.GE.600) THEN
               WRITE(6,*) ' Modified trial vector '
               CALL WRTVCD(VEC1,LU3,1,LBLK)
             END IF
*
            END IF
         END IF
*
*. 1.3 Orthogonalize to current vector
*
         OVLAP  = INPRDD(VEC1,VEC2,LU1,LU3,1,LBLK)
         ONE = 1.0D0
         CALL VECSMD(VEC1,VEC2,-OVLAP,ONE,LU1,LU3,
     &                LU4,1,LBLK)
*
         IF ( IPRT  .GE. 600 ) THEN
           WRITE(6,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
         END IF
*
*. 1.4 Normalize vector
*
         SCALE = INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK)
         FACTOR = 1.0D0/SQRT(SCALE)
         CALL SCLVCD(LU4,LU3,FACTOR,VEC1,1,LBLK)
         IF(IPRT.GE.600) THEN
           WRITE(6,*) '   normalized     (D-E)-1 *( HX - EX ) '
           CALL WRTVCD(VEC1,LU3,1,LBLK)
         END IF
*
       END IF
*
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
        CALL REWINO(LU3)
        CALL MV7(VEC1,VEC2,LU3,LU4)
*. Augment projected matrix
        CALL REWINO(LU1)
* <X!H! Delta>
         APROJ(2) = INPRDD(VEC1,VEC2,LU1,LU4,1,LBLK)
*<Delta!H!Delta>
         APROJ(3) = INPRDD(VEC1,VEC2,LU3,LU4,1,LBLK)
*. Diagonalize projected matrix
      CALL COPVEC(APROJ,WORK(KAPROJ),2*(2+1)/2)
C     write(6,*) ' work(aproj) '
C     call prsym(work(kaproj),2)
      CALL EIGEN(WORK(KAPROJ),AVEC,2,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR **2 + C1VAR **2 )
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
        EIG(ITER,IROOT) = WORK(KAPROJ)
*
C?     WRITE(6,*) ' APROJ(2),APROJ(3)',APROJ(2),APROJ(3)
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)')
     &   ' Eigenvalue and residual of iteration ..', ITER
         WRITE(6,'(2F18.13)') EIG(ITER,1)+EIGSHF, RNORM
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,2)
        WRITE(6,'(2X,E13.7)') EIG(ITER,1)
        CALL WRTMAT(AVEC,2,1,2,1)
      END IF
*
*. Reset
*
      CX = AVEC(1)
      CDELTA = AVEC(2)
*. Eigenvector
      CALL VECSMD(VEC1,VEC2,CX,CDELTA,LU1,LU3,LU5,1,LBLK)
      XNORM = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
      SCALE = 1.0D0/SQRT(XNORM)
      CALL SCLVCD(LU5,LU1,SCALE,VEC1,1,LBLK)
*. Sigma vector
      CXS = CX*SCALE
      CDELTAS = CDELTA*SCALE
      CALL VECSMD(VEC1,VEC2,CXS,CDELTAS,LU2,LU4,LU5,1,LBLK)
      CALL COPVCD(LU5,LU2,VEC1,1,LBLK)
*
      APROJ(1) = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
*
      IF(CONVER) GOTO 1001
 1000 CONTINUE
* ( End of loop over iterations )
 1001 CONTINUE
      ITER = ITERP
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINO(LU1)
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
