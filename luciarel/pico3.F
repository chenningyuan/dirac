      SUBROUTINE PICO3(VEC1,VEC2,LU1,LU2,LU3,LU4,RNRM,EIG,FINEIG,MAXIT,
     &                 NBATCH,LLBATCHB,LLBATCHE,LBLOCK,IBLOCK,IPRTXX,
     &                 NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                 THRES_ET,THRES_EC,THRES_CC,
     &                 E_CONV,C_CONV,ICLSSEL,
     &                 IBLK_TO_CLS,NCLS,CLS_C,CLS_E,CLS_CT,CLS_ET,
     &                 ICLS_A,ICLS_L,IBLKS_A,CLS_DEL,CLS_DELT,ISKIPEI,
     &                 I2BLK,VEC3,ICLS_A2,MXLNG)
C Jeppe REMEMBER to add MXLNG,ICLSA2 to calling sequence
*
* Davidson algorithm , requires three blocks in core
*
* Only three vectors in on DISC
*
* Lu4 should only hold a batch of coefficients
*
*
*
* Block processing version
*
* Jeppe Olsen Winter of 1996
*
* Revision of june 97, PICO2 => PICO3
* Last modifications : Dynamic construction of of  batches
*                      Only relevant sigma blocks constructed
*
* Initial version - Only Diagonal preconditioner,
*
* Special version for NROOT = 1, MAXVEC = 2 !!
*
* Input :
* =======
*        LU1 : Initial  vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU2,LU3   : Scatch files
*        MAXIT     : Largest allowed number of iterations
*        NBATCH    : Number of batches of vector
*        LBATCHB   : Number of blocks in each batch
*        LBATCHE   : Number of elements  in each batch
*        IBLOCK    : Some additional informaition about the blocking
*                    that this routine does not care about !!!!
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
*   THRES_ET   : Convergence criteria for eigenvalue
*
*   THRES_EC   : Threshold for second order energies for individual terms
*   THRES_CC   : Threshold for first  order wavefunction  for individual terms
*
*
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL * 8   INPROD, INPRDD, INPRODB
*
      DIMENSION LLBATCHB(NBATCH),LLBATCHE(NBATCH),IBLOCK(*)
      DIMENSION LBLOCK(*)
*. I2BLK should atleast have length of number of
*. blocks
      DIMENSION I2BLK(*)
*
      DIMENSION RNRM(MAXIT,1),EIG(MAXIT,1),FINEIG(1)
      DIMENSION H0(*),IPNTR(1)
      DIMENSION H0SCR(*)
*. Class and block information
C     IBLK_TO_CLS,NCLS,CLS_C,CLS_E,CLS_A,ICLS_L,IBLKS_A)
      DIMENSION IBLK_TO_CLS(*)
      DIMENSION CLS_C(NCLS),CLS_E(NCLS),CLS_CT(NCLS),CLS_ET(NCLS)
      DIMENSION CLS_DEL(*), CLS_DELT(*)
*.Initial VEC3
      DIMENSION VEC3(*)
      INTEGER ICLS_A(NCLS), ICLS_L(NCLS),IBLKS_A(*),ICLS_A2(NCLS)
*
*     H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
      LOGICAL CONVER
*
      write(6,*) ' PICO2 first 4 rows of IBLOCK'
      call iwrtma(iblock,8,4,8,4)
      IPICO = 0
      IF(IPICO.NE.0) THEN
C?      WRITE(6,*) ' Perturbative solver '
        MAXVEC = MIN(MAXVEC,2)
      ELSE IF(IPICO.EQ.0) THEN
C?      WRITE(6,*) ' Variational  solver '
      END IF
*
      IF(ICLSSEL.EQ.1) THEN
        WRITE(6,*) ' Class selection will be performed '
        WRITE(6,*) ' Number of classes ', NCLS
C?      WRITE(6,*) ' Dimension of each class'
C?      CALL IWRTMA(ICLS_L,1,NCLS,1,NCLS)
      END IF

      IPRT = 10
      IOLSTM = 1
      IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     &WRITE(6,*) ' Inverse iteration modified Davidson '
      IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     &WRITE(6,*) ' Normal Davidson method '
*
      WRITE(6,*) ' LU1 LU2 LU3 LU4 = ', LU1,LU2,LU3,LU4
      IF(IPRT.GE.1) THEN
        WRITE(6,*) ' Convergence threshold for eigenvalues', THRES_ET
        WRITE(6,*)
        WRITE(6,*) ' Elements of trial vectors discarded if '
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*)
     &  '    Estimate of contribution to wavefunction is less than ',
     &  THRES_CC
        WRITE(6,*)
     &  '    Estimate of contribution to Energy is less than ',
     &  THRES_EC
      END IF
      IF(IPRT.GE.10)
     &WRITE(6,*) ' Number of batches of vector ', NBATCH
*
*. Total number of blocks
      NBLOCKT = 0
      DO IBATCH = 1, NBATCH
        NBLOCKT = NBLOCKT + LLBATCHB(IBATCH)
      END DO
      IF(IPRT.GE.10)
     &WRITE(6,*) ' Total number of BLOCKS ', NBLOCKT
      TEST = 1.0D-6
      CONVER = .FALSE.
      IROOT = 1
      NROOT = 1
      ZERO = 0.0D0
*. Play around with dynamic alloxation of batches ...
      IDYNBATCH = 1
*
* ===================
*.Initial iteration
* ===================
*
      ITER = 1
*. Obtain energy of initial vector
      IF(ISKIPEI.EQ.0) THEN
*. active classses of initial vector
*
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU1,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
C              FIND_ACTIVE_CLASSES(LUIN,LBLK,I_BLK_TO_CLS,ICLS_A,NCLS,SEGMNT)
*. Mark active blocks and find required number of batches
C     REPART_CIV(IBATCH,NBATCH,LBATCH,LEBATCH,
C    &                  MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLKS)
          write(6,*) ' MXLNG = ',MXLNG
          CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
        ELSE
          NBATCHL = NBATCH
        END IF
        E = 0.0D0
        IRESTRICT = 1
        IOFF = 1
        DO IBATCH = 1, NBATCHL
          LBATCHB = LLBATCHB(IBATCH)
          LBATCHE = LLBATCHE(IBATCH)
          IF(IPRT.GE.10)
     &    WRITE(6,*) '  <Delta ! H ! Delta >, Batch : ',IBATCH
          CALL SETVEC(VEC1,ZERO,LBATCHE)
C         SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC)
          CALL SBLOCK(LBATCHB,IBLOCK,IOFF,VEC2,VEC1,LU1,IRESTRICT,0)
          IF(IPRT.GE.200) THEN
            WRITE(6,*) ' Initial batch of S, number', IBATCH
            CALL WRTBLKN(VEC1,LBATCHB,LBLOCK(IOFF))
          END IF
*. Obtain corresponding C blocks
          CALL GET_BLOCKS_FROM_DISC
     &    (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
          E = E + INPROD(VEC1,VEC2,LBATCHE)
          IF(IPRT.GE.200) THEN
            WRITE(6,*) ' Initial batch of C, number', IBATCH
            CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
          END IF
          IOFF = IOFF + LBATCHB
        END DO
        IF(IRESTRICT.EQ.1) E = 2*E
        EIG(1,IROOT) = E
      ELSE IF(ISKIPEI.EQ.1) THEN
        E = FINEIG(IROOT)
        WRITE(6,*) ' Initial energy obtained from FINEIG as ',E
        EIG(1,IROOT) = E
      END IF

*
      IF(IPRT .GE. 3 ) THEN
        WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
        WRITE(6,'(5F18.13)')
     &  ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
      END IF
      ITERX = 1
*
* ======================
*. Loop over iterations
* ======================
*
      DO 1000 ITER = 2, MAXIT
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*) ' Info from iteration .... ', ITER
       END IF
       ZERO = 0.0D0
       IF(ICLSSEL.EQ.1) THEN
         CALL SETVEC(CLS_CT,ZERO,NCLS)
         CALL SETVEC(CLS_ET,ZERO,NCLS)
         CALL SETVEC(CLS_DELT,ZERO,NCLS)
         CALL SETVEC(CLS_C ,ZERO,NCLS)
         CALL SETVEC(CLS_E ,ZERO,NCLS)
         CALL SETVEC(CLS_DEL,ZERO,NCLS)
       END IF
*. Active classes
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU1,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
C              FIND_ACTIVE_CLASSES(LUIN,LBLK,I_BLK_TO_CLS,ICLS_A,NCLS,SEGMNT)
*. Mark active blocks and find required number of batches
C     REPART_CIV(IBATCH,NBATCH,LBATCH,LEBATCH,
C    &                  MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLKS)
          CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
        ELSE
          NBATCHL = NBATCH
        END IF
       EIGAPR = EIG(ITER-1,1)
* ===============================
*. Obtain C(T) (H0-E)**-1 C
* ===============================
       GAMMA = 0.0D0
       IOFF = 1
       CALL REWINO(LU1)
       DO IBATCH = 1, NBATCHL
         LBATCHB = LLBATCHB(IBATCH)
         LBATCHE = LLBATCHE(IBATCH)
*. Retrieve Batch of C
         NO_ZEROING = 0
         NO_ZEROING1= 1
         CALL FRMDSCN3(VEC1,LBATCHB,LBLK,LU1,NO_ZEROING1,I2BLK(IOFF),
     &                 LBLOCK(IOFF))
         CALL COPVEC(VEC1,VEC2,LBATCHE)
*. Multiply with (H0-E)** -1
         FACTOR = -EIGAPR
         ITASK = 1
         CALL  DIATERM_GAS(FACTOR,ITASK,VEC2,LBATCHB,IBLOCK,IOFF,0,
     &         NO_ZEROING1,I2BLK(IOFF))
         IF(NO_ZEROING1.EQ.0) THEN
           GAMMA = GAMMA + INPROD(VEC1,VEC2,LBATCHE)
         ELSE
C                           INPRODB(VEC1,VEC2,NBLK,LBLK,I0BLK)
            GAMMA = GAMMA + INPRODB(VEC1,VEC2,LBATCHB,LBLOCK(IOFF),
     &                      I2BLK(IOFF))
         END IF
         IOFF = IOFF + LBATCHB
       END DO
       IF(IPRT.GE.10)
     & WRITE(6,*) ' Gamma  calculated ',GAMMA

*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = (H0-E)-1 (H*C - EIGAPR*C) : Obtain in batches and save on DISC
*
       EIGAPR = EIG(ITER-1,1)
       NELMNT = 0
       NZERO = 0.0D0
*
       RNORM = 0.0D0
       DELTA = 0.0D0
       CHEDEL =0.0D0
       DELNORM = 0.0D0
*
       DELTAT = 0.0D0
       CHEDELT =0.0D0
       DELNORMT = 0.0D0
       ECC = 0.0D0
*
       CALL REWINO(LU2)
       IOFF = 1
*. Find partitioning of sigma
*. BLocks obttained by double excitations from classes in C
       NEXC  = 2
       CALL EXCCLS(NCLS,ICLS_A,ICLS_A2,NEXC)
C           EXCCLS(NCLS,IACTIN,IACTOUT,IEXC)
*. Partitioning of sigma vector
        CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A2,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)

       DO IBATCH = 1, NBATCHL
         LBATCHB = LLBATCHB(IBATCH)
         LBATCHE = LLBATCHE(IBATCH)
         NELMNT = NELMNT + LBATCHE
         IF(IPRT.GE.10)
     &   WRITE(6,*) '  Delta, Batch : ',IBATCH
* Batch of HC in VEC1
         ZERO = 0.0D0
         CALL SETVEC(VEC1,ZERO,LBATCHE)
         CALL SBLOCK(LBATCHB,IBLOCK,IOFF,VEC2,VEC1,LU1,0,0)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' Batch of H C '
           CALL WRTBLKN(VEC1,LBATCHB,LBLOCK(IOFF))
         END IF
*. Retrieve Batch of C
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' C batch read in '
           CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
         END IF
*. Update energy
         ECC = ECC + INPROD(VEC1,VEC2,LBATCHE)
* Batch of (H-E)C in VEC1
         ONE = 1.0D0
         FACTOR = -EIGAPR
         CALL VECSUM(VEC1,VEC1,VEC2,ONE,FACTOR,LBATCHE)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' Batch of (H - E ) C '
           CALL WRTBLKN(VEC1,LBATCHB,LBLOCK(IOFF))
         END IF
*. Norm of residual
         RNORM = RNORM + INPROD(VEC1,VEC1,LBATCHE)
*. Batch of (H0-E)-1(H-E)C  in VEC2
         CALL COPVEC(VEC1,VEC2,LBATCHE)
         FACTOR = -EIGAPR
         ITASK = 1
         I12 = 1
         CALL  DIATERM_GAS(FACTOR,ITASK,VEC2,LBATCHB,IBLOCK,IOFF,0,0,0)
         DELNORMT = DELNORMT + INPROD(VEC2,VEC2,LBATCHE)
* C(H-E)(H0-E0)-1(H-E)C
         CHEDELT = CHEDELT + INPROD(VEC1,VEC2,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
*. Contributions divided into occupation classes, complete expansion
*. Wave function correction
           CALL CLASS_PROD3(VEC2,VEC2,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_CT)
*. Energy correction
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_ET)
         END IF
*.[(H0-E0)-1(H-E)C]_{truncated} and
*.C(H-E) [(H0-E0)-1(H-E)C]_{truncated}
         IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
           ZERO = 0.0D0
           CALL SETVEC(VEC3,ZERO,LBATCHE)
           NZERO = NZERO + LBATCHE
           DO I = 1, LBATCHE
             IF(ABS(VEC2(I)*VEC1(I)).GE.THRES_EC.OR.
     &          ABS(VEC2(I)).GE.THRES_CC           ) THEN
               CHEDEL = CHEDEL + VEC1(I)*VEC2(I)
               DELNORM = DELNORM + VEC2(I)*VEC2(I)
               VEC3(I)=VEC2(I)
               NZERO = NZERO - 1
             END IF
           END DO
*
           IF(ICLSSEL.EQ.1) THEN
*. Contributions divided into occupation classes, truncated expansion
*. Wave function correction
             CALL CLASS_PROD3(VEC3,VEC3,IOFF,LBATCHB,IBLOCK,
     &                        IBLK_TO_CLS,NCLS,CLS_C)
*. Energy correction
             CALL CLASS_PROD3(VEC1,VEC3,IOFF,LBATCHB,IBLOCK,
     &                        IBLK_TO_CLS,NCLS,CLS_E)
           END IF
         END IF
*. retrieve c batch from disc
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC1,1)
* C(H0-E0)-1(H-E)C
         DELTAT = DELTAT + INPROD(VEC1,VEC2,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_DELT)
         END IF
         IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
* C[(H0-E0)-1(H-E)C]{truncated}
           DELTA = DELTA + INPROD(VEC1,VEC3,LBATCHE)
           IF(ICLSSEL.EQ.1) THEN
             CALL CLASS_PROD3(VEC1,VEC3,IOFF,LBATCHB,IBLOCK,
     &                        IBLK_TO_CLS,NCLS,CLS_DEL)
           END IF
         END IF
*
         IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
*. Write packed version to DISC
*. Pack out so zero blocks are given zero entries
*JEPPE!!!!!
           CALL TODSCNP(VEC3,LBATCHB,LBLOCK(IOFF),LBLK,LU2)
         ELSE
*. Write complete batch of delta to disc
           CALL TODSCN(VEC2,LBATCHB,LBLOCK(IOFF),LBLK,LU2)
           CHEDEL = CHEDELT
           DELTA = DELTAT
           DELNORM = DELNORMT
         END IF
         IF(IPRT.GE.200) THEN
           WRITE(6,*) ' Batch of blocks of trial vector '
           CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
         END IF
         IOFF = IOFF + LBATCHB
       END DO
*      /\ End of loop over batches of correction vector
       WRITE(6,*) ' ECC = ',ECC
       CALL ITODS(-1,1,LBLK,LU2)
       RNRM(ITER-1,1) = SQRT(RNORM)
*. (End of loop over batches of (H0-E)-1(H-E)C)
*. Predicted energy
       IF(ICLSSEL.EQ.1) THEN
         IF(THRES_EC.EQ.0.0D0.AND.THRES_CC.EQ.0.0D0) THEN
*. No truncation of small terms were performed, so the energy of
*. term truncated expansion equals energy of untruncated expansion
           CALL COPVEC(CLS_CT,CLS_C,NCLS)
           CALL COPVEC(CLS_ET,CLS_E,NCLS)
         END IF
*. decide which classes should be truncated
         CALL CLASS_TRUNC(NCLS,ICLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                    E_CONV,ICLS_A,N_CLS_TRN,E_CLS_TRN,W_TRN)
         WRITE(6,*) ' Home from class truncation'
         WRITE(6,*) 'N_CLS_TRN,E_CLS_TRN',N_CLS_TRN,E_CLS_TRN
*. Update corrections for class elimination
         N_ACT_CLS = 0
         CHEDEL = CHEDEL - (-1.0D0) * E_CLS_TRN
         DO JCLS = 1, NCLS
           IF(ICLS_A(JCLS).EQ.0) THEN
             DELTA = DELTA-CLS_DEL(JCLS)
           ELSE
             N_ACT_CLS = N_ACT_CLS + 1
           END IF
         END DO

*. And do the truncation
*. Truncation of classes => truncation of blocks
         IF(N_CLS_TRN.NE.0) THEN
           CALL CLS_TO_BLK(NBLOCKT,IBLK_TO_CLS,ICLS_A,IBLKS_A)
*. from LU2 to LU3 and back to LU2
           CALL ZAP_BLOCK_VEC(LU2,LBLK,IBLKS_A,VEC2,LU3)
C               ZAP_BLOCK_VEC(LUIN,LBLK,IBLKS_A,SEGMNT,LUSCR)
         END IF
       ELSE
         N_CLS_TRN = 0
         E_CLS_TRN = 0.0D0
         W_CLS_TRN = 0.0D0
       END IF
*. Predicted energy
       IF(GAMMA.NE.0.0D0) THEN
         E2PREDIT = - CHEDELT + DELTAT**2/GAMMA
         E2PREDI  = - CHEDEL  + DELTA * DELTAT /GAMMA
       ELSE
         E2PREDIT = - CHEDELT
         E2PREDI  = - CHEDEL
         IF(ICLSSEL.EQ.1) THEN
           CALL ICOPVE(CLS_CT,CLS_C,NCLS)
           CALL ICOPVE(CLS_ET,CLS_E,NCLS)
         END IF
       END IF
*
*.
       WRITE(6,*) ' Information for untruncated expansion:'
       WRITE(6,*) ' ======================================'
       WRITE(6,*)
     & ' CHEDELT DELTAT GAMMA ', CHEDELT,DELTAT,GAMMA
       IF(GAMMA.NE.0.0D0) THEN
         WRITE(6,*)
     & '   Orthogonalization(Total) term to E2 ', DELTAT**2/GAMMA
       END IF
       WRITE(6,*)
     & ' Predicted energy(no truncation), change and total ',
     &              E2PREDIT,EIGAPR+E2PREDIT+EIGSHF
       IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
         WRITE(6,*)
         WRITE(6,*) ' Information for truncated expansion:'
         WRITE(6,*) ' ======================================'
         WRITE(6,*)
     &   ' CHEDEL DELTA GAMMA ', CHEDEL,DELTA,GAMMA
         IF(GAMMA.NE.0.0D0) THEN
           WRITE(6,*)
     &   '   Orthogonalization term to E2 ', DELTA**2/GAMMA
         END IF
         WRITE(6,*)
     &   ' Predicted energy (truncated), change and total ',
     &                E2PREDI,EIGAPR+E2PREDI+EIGSHF
         WRITE(6,*) ' Estimated norm of eliminated terms ',
     &   DELNORMT-DELNORM
         WRITE(6,*)
     &   ' Estimated energy contributions of eliminated terms',
     &   E2PREDIT-E2PREDI
         WRITE(6,*) ' Number of zero elements in delta ',
     &    NZERO
         WRITE(6,*) ' Number of nonzero terms in delta ',
     &    NELMNT-NZERO
       ELSE
       END IF
*
       IF(ICLSSEL.EQ.1) THEN
         IF(N_ACT_CLS .EQ. 0   ) THEN
*. All classes were eliminated so we are home -and hopefully dry
           WRITE(6,*) ' No active classes  '
           WRITE(6,*) ' I will therefore end the diagonalization'
           GOTO 1001
         END IF
       END IF
*
* ============================================
* 1.5 : Inverse Iteration Correction to Delta
* ============================================
*
*
* Update delta to
* -(H0-E0)-1(H-E)|0> + delta/gamma * (H0-E0)-1 |0>
*
* ( was +(H0-E0)-1(H-E)|0>)
*
       CALL REWINO(LU3)
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       IOFF = 1
       IF(IOLSTM.EQ.1.AND.ABS(GAMMA).GT.1.0D-6) THEN
         WRITE(6,*) ' Inverse iteration correction will be added '
         DO IBATCH = 1, NBATCH
           LBATCHB = LLBATCHB(IBATCH)
           LBATCHE = LLBATCHE(IBATCH)
*. Retrieve Batch of C
           NO_ZEROING = 0
           CALL FRMDSCN3(VEC1,LBATCHB,LBLK,LU1,NO_ZEROING,I2BLK(IOFF),
     &                   LBLOCK(IOFF))
C          CALL GET_BLOCKS_FROM_DISC
C    &     (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC1,1)
*. Multiply with (H0-E)** -1
           FACTOR = -EIGAPR
           ITASK = 1
           CALL  DIATERM_GAS(FACTOR,ITASK,VEC1,LBATCHB,IBLOCK,IOFF,0
     &                       0,0)
*. Retrieve Batch of Delta
           NO_ZEROING = 0
           CALL FRMDSCN3(VEC2,LBATCHB,LBLK,LU2,NO_ZEROING,I2BLK(IOFF),
     &                   LBLOCK(IOFF))
C          CALL GET_BLOCKS_FROM_DISC
C    &     (LU2,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
*. And add
           FAC1 = -1.0D0
           FAC2 = DELTA/GAMMA
           CALL VECSUM(VEC2,VEC2,VEC1,FAC1,FAC2,LBATCHE)
*. Transfer to Disc
           IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
             CALL TODSCNP(VEC2,LBATCHB,LBLOCK(IOFF),LBLK,LU3)
           ELSE
             CALL TODSCN(VEC2,LBATCHB,LBLOCK(IOFF),LBLK,LU3)
           END IF
           IOFF = IOFF + LBATCHB
         END DO
*        ^ End of loop over batches
         CALL ITODS(-1,1,LBLK,LU3)
*. Well, it is nice to have the correction vector on LU2 so
         IREW = 1
         IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
            CALL COPVCDP(LU3,LU2,VEC1,IREW,LBLK)
         ELSE
            CALL COPVCD(LU3,LU2,VEC1,IREW,LBLK)
         END IF
       END IF
*
* ===================================
* 2 : Calculate <Delta ! H ! Delta >
* ===================================
*
*. Active classes
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU2,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
          CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
        ELSE
          NBATCHL = NBATCH
        END IF
* Loop over batches of H !delta>
       IOFF = 1
       DELHDEL = 0.0D0
       IRESTRICT = 1
       DO IBATCH = 1, NBATCHL
         IF(IPRT.GE.10)
     &   WRITE(6,*) '  <Delta ! H ! Delta >, Batch : ',IBATCH
         LBATCHB = LLBATCHB(IBATCH)
         LBATCHE = LLBATCHE(IBATCH)
         CALL SETVEC(VEC1,ZERO,LBATCHE)
         CALL SBLOCK(LBATCHB,IBLOCK,IOFF,VEC2,VEC1,LU2,IRESTRICT,0)
*. Retrieve Batch of Delta
         CALL GET_BLOCKS_FROM_DISC
     &   (LU2,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
         IF(IPRT.GE.1000) THEN
           WRITE(6,*) ' Blocks of delta retrieved '
           CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
         END IF
         DELHDEL = DELHDEL + INPROD(VEC1,VEC2,LBATCHE)
         IOFF = IOFF + LBATCHB
       END DO
       IF(IRESTRICT.EQ.1) DELHDEL = 2.0D0*DELHDEL
*
* ===========================================
* 3 : Solve 2 by 2 problem : Nonorthogonal !!
* ===========================================)
*
*      Norm of delta and overlap between 0 and delta
       S12 = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
       S22 = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
C      H11 = EIGAPR
       H11 = ECC
       IF(IOLSTM.EQ.1.AND.ABS(GAMMA).GT.1.0D-6) THEN
        H12 = -CHEDEL + DELTA*DELTAT/GAMMA
       ELSE
        H12 = CHEDEL + EIGAPR*DELTA
       END IF
       H22 = DELHDEL
*
       S11 = 1.0D0
*
*.( H11  H12 ) (X1)       (S11    S12 )(X1)
* (          ) (  )   = E (           )( )
* ( H12  H22 ) (X2)       (S12    S22 )(X2)
*
* The eigenvalues
*
        A = S11*S22 -S12 **2
        B = 2*S12*H12-S11*H22-S22*H11
        C = H11*H22 -H12**2
*
        EA = -B/(2*A) - SQRT(B**2 - 4*A*C)/(2*A)
        EB = -B/(2*A) + SQRT(B**2 - 4*A*C)/(2*A)
*. And the lowest eigenvalue is
        E = MIN(EA,EB)
*. The corresponding eigenvector
*. Intermediate normalization
        X1 = 1.0D0
        X2 = -(H11-E*S11)/(H12-E*S12)
*. Normalized
        XNORM2 = S11*X1**2 + S22*X2**2 + 2.0*S12*X1*X2
        XNORM = SQRT(XNORM2)
        X1 = X1/XNORM
        X2 = X2/XNORM
*
        IF(IPRT.GE.10) THEN
          WRITE(6,*) ' 2 X 2 Generalized eigenvalue problem, H and S '
          WRITE(6,*)
          WRITE(6,'(4X,E20.10)') H11
          WRITE(6,'(4X,2E20.10)') H12,H22
          WRITE(6,*)
          WRITE(6,'(4X,E20.10)') S11
          WRITE(6,'(4X,2E20.10)') S12,S22
*
          WRITE(6,*) ' Lowest eigenvalue (with shift)', E + EIGSHF
          WRITE(6,*) ' Corresponding eigenvector ', X1,X2
        END IF

*. Save corresponding eigenvector on file LU1 ( first LU3, then COPY)
* VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        IREW = 1
        IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
          CALL VECSMDP(VEC1,VEC2,X1,X2,LU1,LU2,LU3,IREW,LBLK)
          CALL COPVCDP(LU3,LU1,VEC1,IREW,LBLK)
        ELSE
          CALL VECSMD(VEC1,VEC2,X1,X2,LU1,LU2,LU3,IREW,LBLK)
          CALL COPVCD(LU3,LU1,VEC1,IREW,LBLK)
        END IF
        EIG(ITER,1) = E
        ITERX = ITER
*. Convergence ??
      IF(ABS(EIG(ITER,1) - EIG(ITER-1,1)).LE.THRES_ET)
     &   CONVER = .TRUE.
      IF(CONVER) GOTO 1001
 1000 CONTINUE
* ( End of loop over iterations )
 1001 CONTINUE
      ITER = ITERX
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
          DO 1330 I=1,ITER-1
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          WRITE(6,1340) ITER,EIG(ITER,IROOT)+EIGSHF
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        DO 1601 IROOT = 1, NROOT
COLD       FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           FINEIG(IROOT) = EIG(ITER,IROOT)
 1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER-1,IROOT)
 1607   CONTINUE
      END IF
C
*. Clean up : IBLOCK in original form
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &     ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
