C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T.
CFleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V.
CLarsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T.
CWinther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the
Cresearch
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
#if defined (VAR_MPI)
***********************************************************************
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE MICDV4_ENLMD_REL_PAR(T_CC,VEC1,VEC2,LU1,LU2,
     &                  RNRM,EIG,FINEIG,MAXIT,
     &                  LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,IBLOCKL,NPARBLOCK,NBLOCKS,
     &                  LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LUCLIST,NBATV,
     &                  LBATV,LEBATV,I1BATV,IBATV)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Version using H0 + Lambda V as Sigma routine
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
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
#endif
#include "lucipar_rel.inc"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1), IBLOCKL(*), NPARBLOCK(*)
       DIMENSION H0SCR(*),T_CC(*)
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
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*
       call qenter('MICDV4')
*
       WRITE(LUWRT,*) '             ::::::::::::::::::::::::::::   '
       WRITE(LUWRT,*) '               Entering par_MICDV4 (rel)    '
       WRITE(LUWRT,*) '             ::::::::::::::::::::::::::::   '
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(LUWRT,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(LUWRT,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(LUWRT,*) ' Sorry MICDV4 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(LUWRT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(LUWRT,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(LUWRT,*) ' Enforced stop on MICDV4 '
         STOP 20
       END IF
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
*
       DO 10 IVEC = 1,NINVEC
*        
*        copy c-vector to working-file ILUC
*
         CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
CSK         CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,LBATV,LEBATV,I1BATV,IBATV,
CSK     &                       MY_LU1_OFF,LU1LIST,LUCLIST,IVEC-1,1)
*
*        start calculation: sigma = H x C
*
************************************************************************
         call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,1)
************************************************************************
*
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU4)
         CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
         DO 8 JVEC = 1, IVEC
           CALL REWINO(LU3)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU3,0,LBLK)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.3 ) THEN
           WRITE(LUWRT,*) ' INITIAL PROJECTED MATRIX  '
           CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL DCOPY(NINVEC*(NINVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
       CALL EIGEN(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*

       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRT,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(LUWRT,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(LUWRT,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        CALL GETTIM(CPUITR1,WALLITR1)
        write(LUWRT,*)
        write(LUWRT,'(A21,3X,I3)') ' Info from iteration ',ITER
        write(LUWRT,*) '_______________________'
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
       DO 100 IROOT = 1, NROOT
         EIGAPR = EIG(ITER-1,IROOT)
*
         CALL REWINO(LU1)
         CALL REWINO(LU2)
         EIGAPR = EIG(ITER-1,IROOT)
         DO 60 IVEC = 1, NVEC
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           IF(IVEC.EQ.1) THEN
             CALL REWINO( LU3 )
*                 SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
             CALL SCLVCD(LU2,LU3,FACTOR,VEC1,0,LBLK)
           ELSE
             CALL REWINO(LU3)
             CALL REWINO(LU4)
C                 VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU2,LU3,0,LBLK)
           END IF
C
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(LUWRT,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
         END IF
*  Strange place to put convergence but ....
C                      INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
*
         WRITE(LUWRT,'(A19,7X,I3,3X,1E18.13,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
         CALL FLSHFO(LUWRT)
*
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
*. (H0-E) -1 *(HX-EX) on LU3
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(LUWRT,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
*. Assuming diagonal preconditioner
           IPRECOND = 1
           CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
C               H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND)
           IF ( IPRT  .GE. 600) THEN
             WRITE(LUWRT,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* Current eigen-vector on LU5
             CALL REWINO(LU1)
             DO 66 IVEC = 1, NVEC
               FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
               IF(IVEC.EQ.1) THEN
                 IF(NVEC.EQ.1) THEN
                   CALL REWINO( LU5 )
                   CALL SCLVCD(LU1,LU5,FACTOR,VEC1,0,LBLK)
                 ELSE
                   CALL REWINO( LU4 )
                   CALL SCLVCD(LU1,LU4,FACTOR,VEC1,0,LBLK)
                 END IF
               ELSE
                 CALL REWINO(LU5)
                 CALL REWINO(LU4)
                 CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                 CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
               END IF
   66        CONTINUE
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(LUWRT,*) '  (current  X ) '
               CALL WRTVCD(VEC1,LU5,1,LBLK)
             END IF
* (H0 - E )-1  * X on LU4
             CALL REWINO(LU5)
             CALL REWINO(LU4)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD(LU4,LUDIA,LU5,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPRDD(VEC1,VEC2,LU5,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU5,LU4,1,LBLK))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(LUWRT,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU5,LU3,1,LBLK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(LUWRT,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(LUWRT,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
                CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(LUWRT,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU3,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
           CALL REWINE( LU1 ,LBLK)
           DO 80 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             WORK(IVEC) = INPRDD(VEC1,VEC2,LU1,LU3,0,LBLK)
C?       WRITE(6,*) ' MICDV4 : Overlap ', WORK(IVEC)
   80      CONTINUE
*
           CALL REWINE(LU1,LBLK)
           DO 82 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             CALL REWINE(LU4,LBLK)
             CALL VECSMD(VEC1,VEC2,-WORK(IVEC),1.0D0,LU1,LU3,
     &                   LU4,0,LBLK)
             CALL COPVCD(LU4,LU3,VEC1,1,LBLK)
   82      CONTINUE
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(LUWRT,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU3,LU4,FACTOR,VEC1,1,LBLK)
             WRITE(LUWRT,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
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
      IF( MYPROC .eq. MASTER ) WRITE(LUWRT,*) 'THIS IS NVEC',NVEC,MYPROC
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,NVEC,VEC1,1,LBLK)
      DO 150 IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
CTF     CALL ENLMD(VEC1,VEC2,LU3,LU4)
        call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,1)
        CALL REWINE(LU4,LBLK)
        CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO 140 JVEC = 1, NVEC+IVEC
          CALL REWINE(LU4,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU4,0,LBLK)
  140   CONTINUE
  150 CONTINUE
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      write(LUWRT,*) 'Dimension of projected matrix:',NVEC
      CALL DCOPY(NVEC*(NVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
      CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(LUWRT,*) ' Var and Pert solution, energy and coefficients'
        WRITE(LUWRT,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(LUWRT,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(LUWRT,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(LUWRT,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(LUWRT,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRT,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset      
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
      IF( MYPROC .eq. MASTER ) WRITE(LUWRT,*) 
     & 'I am here',NVEC,NROOT,MAXVEC,CONVER
        
        CALL REWINE( LU5,LBLK)
        DO 320 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  320   CONTINUE
*. Transfer C vectors to LU1
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 411 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  411   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 329 IROOT = 1, NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  329   CONTINUE
*
* Transfer HC's to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 400 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  400   CONTINUE
        NVEC = NROOT
*
        CALL DZERO(AVEC,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
*
      END IF
*
*  Timing of this iteration
      CALL GETTIM(CPUITR2,WALLITR2)
      CPUTID = SECTID(CPUITR2-CPUITR1)
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,9300) CPUTID,WALLTID
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
 
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(LUWRT,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(LUWRT,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 0 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(LUWRT,*)
          WRITE(LUWRT,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(LUWRT,*)
     &    '============================================'
          WRITE(LUWRT,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(LUWRT,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(LUWRT,1200)
 1200       FORMAT(' The final approximation to eigenvector ')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(LUWRT,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(LUWRT,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(LUWRT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        CALL REWINE(LU1,LBLK)
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
        DO 1601 IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
          write(LUWRT,9000) ' Final Eigenvalue ',IROOT,' : ',
     & FINEIG(IROOT)
          IF(IPRT.GE.400) THEN
            WRITE(LUWRT,1200)
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
 1601   CONTINUE
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
      END IF
*
 9000 format(A18,1I3,A3,F20.13)
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(LUWRT,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
*
      call qexit('MICDV4')
*
* Delete scratch units
      close(unit=LU3,status='DELETE')
      close(unit=LU4,status='DELETE')
      close(unit=LU5,status='DELETE') ! think about that!!!!!!!!
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
 9300 FORMAT(' >>>  CPU (WALL) TIME IN ITERATION: ',A,'(',A,')')
      END
*
#else
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_DUMMY_REL2
*     dummy routine for normal compilation
      END
#endif
