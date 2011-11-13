!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

#if defined (VAR_MPI2)
***********************************************************************

      SUBROUTINE MICDV4_ENLMD_REL_PAR(T_CC,T_BUFF,VEC1,VEC2,
     &                  RNRM,EIG,FINEIG,MAXIT,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,IPNTR,LBLK,EIGSHF,
     &                  THRES_E,IBLOCKL,NPARBLOCK,NBLOCKS,
     &                  LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LUCLIST,NBATV,
     &                  LBATV,LEBATV,I1BATV,IBATV,RCCTOS,SCRRED,
     &                  IGROUPLIST,IPROCLIST,IT_TTPL,IT_TTOL)
*
*     Davidson algorithm , requires two blocks in core
*     parallel Multi root version
*
*     orig. written by Jeppe Olsen Winter of 1991
*
*     updated to allow general preconditioner, October 1993
*
*     Version using H0 + Lambda V as Sigma routine
*
* Input :
* =======
*        ILU1      : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        IDIA      : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
* For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
*           LBLK : Defines block structure of matrices
* On input ILU1 is supposed to hold initial guesses to eigenvectors
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
#include "parluci.h"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION IPNTR(1), IBLOCKL(*), NPARBLOCK(*)
       DIMENSION T_CC(*), SCRRED(*), T_BUFF(*)
       INTEGER   RCCTOS(*)
       DIMENSION IT_TTPL(*) 
       INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2                               
*      SCRRED : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
*
       DIMENSION FINEIG(NROOT)
       LOGICAL CONVER,RTCNV(NROOT)
       REAL*8 INPRDD
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12, WALLTSTEP*12
       CHARACTER WPART22*12, FILELAB*8
       INTEGER NZERO,LZERO_SCRRED
       INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR
*
       call qenter('MICDV4')
*
       WRITE(LUWRT,*) '                                            '
       WRITE(LUWRT,*) '             ============================   '
       WRITE(LUWRT,*) '               Entering rel. MICDV4_par     '
       WRITE(LUWRT,*) '             ============================   '
       WRITE(LUWRT,*) '                                            '
C
C      no screening 
       CSCREEN = .FALSE.
C
C      do more explicit timing
       TIMING = .false.
C
C      set check point parameters - approx after every 5th iteration 
C      of MAXIT
       ITER_CHECKP = 0
       ICHPARAM    = 0
       IF( CHECKPOINT_LUCIX )THEN
         ICHPARAM    = 5
       END IF
C
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
*      initialize LZERO_SCRRED
       LZERO_SCRRED = 0
       LZERO_SCRRED = MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2 )
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(LUWRT,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(LUWRT,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(LUWRT,*) ' Sorry MICDV4 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(LUWRT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(LUWRT,*) ' Raise MXCIV_CI to be at least 2 * Nroot '
         WRITE(LUWRT,*) ' Enforced stop on MICDV4 '
         STOP 20
       END IF
*
       NVEC = 0
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
       NZERO = 0
*
*      ===================
*       Initial iteration
*      ===================
*
*      start timing of initial iteration
       WALLITR1 = MPI_WTIME()
        
       CALL DZERO(APROJ, MAXVEC*(MAXVEC+1)/2 )
*
       ITER        = 1
       ITER_CHECKP = 1
*
       DO 10 IVEC = 1,NINVEC
*        
*        copy c-vector to working-file ILUC
*
         CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
         CALL IZERO(LUCLIST,NUM_BLOCKS2)
*
*         WRITE(LUWRT,*) 'LUCLIST for the 1st time:'
*         CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUWRT)
*         WRITE(LUWRT,*) 'LU1LIST for the 1st time:'
*         CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
*
*        reset LUCLIST; MY_ACT_BLK_ALL = NUMBLOCKS!
*
         CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                       IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                       LUCLIST,IBLOCKL,IVEC-1)
CSK         CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
*         WRITE(LUWRT,*) 'LUCLIST for the 2nd time:'
*         CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUWRT)
*         WRITE(LUWRT,*) 'LU1LIST for the 2nd time:'
*         CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
*
*        set offset for sigma-file
*
         JVEC_SF = IVEC - 1
*
*        start calculation: sigma = H x C
*
*        timing this sigma-vector computation
*
         sigmatime = MPI_WTIME()
*
*=======================================================================
         call sigden_ctrl(VEC1,VEC2,ILUC,ILU2,T_CC,T_BUFF,1
#if defined (VAR_MPI2)
     &                   ,LUCLIST,LU2LIST,IBLOCKL,NPARBLOCK,
     &                    IGROUPLIST,IPROCLIST,RCCTOS,
     &                    IT_TTPL,IT_TTOL
     &                   ,IBLOCKL,NPARBLOCK
#endif
     &                   )
*=======================================================================
*
*        end of timing
         sigmatime2 = MPI_WTIME()
         WALLTID = SECTID(sigmatime2-sigmatime)
         WRITE(LUWRT,9777) WALLTID
         WRITE(LUWRT,*) '   '
*
*        projected matrix using batch structure of CI vector(s)
*
         CALL INPROD_B_PAR_RL(ILU1,ILU2,VEC1,VEC2,APROJ,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                        LU2LIST,IVEC)
*
   10  CONTINUE
*
*      timing of initial iteration
*
       WALLITR2 = MPI_WTIME()
       WALLTID = SECTID(WALLITR2-WALLITR1)
       WRITE(LUWRT,9888) WALLTID
       WRITE(LUWRT,*) '   '
*
*      synchronize MPI_COMM_WORLD
* 
       IREDL = 0
       IREDL = NINVEC*(NINVEC-1)/2 + NINVEC
       CALL DZERO(SCRRED,IREDL)
       CALL REDVEC_REL(APROJ,SCRRED,IREDL,2,MPI_SUM,MPI_COMM_WORLD,-1)
       CALL DCOPY(IREDL,SCRRED,1,APROJ,1)
*
*
       IF( IPRT .GE.3 ) THEN
         WRITE(LUWRT,*) ' INITIAL PROJECTED MATRIX  '
         CALL WRTMATMN(APROJ,1,IREDL,1,IREDL,LUWRT)
       END IF
*
*      diagonalize initial projected matrix
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
         CALL WRTMATMN(EIG(1,1),1,NROOT,MAXIT,NROOT,LUWRT)
       END IF
CSK       WRITE(LUWRT,*) 'NVEC, NINVEC' , NVEC, NINVEC
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
*      assuming diagonal preconditioner
*
       IPRECOND = 1
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
*       start timing of iteration
        WALLITR1 = MPI_WTIME()
        starttime = MPI_WTIME()
*
        write(LUWRT,*)
        write(LUWRT,'(A21,3X,I3)') ' Info from iteration ',ITER
        write(LUWRT,*) '_______________________'
        ITER        = ITER + 1
        ITER_CHECKP = ITER_CHECKP + 1
*
*
*===========================================================
*                       PART 1                             =
*                                                          =
*              New directions to be included               =
*                                                          =
*===========================================================
*
*      1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       DO 100 IROOT = 1, NROOT
*
*        reset scratch file lists ...
*
         CALL IZERO(LU3LIST,IALL_LU3)
         CALL IZERO(LU4LIST,IALL_LU4)
         CALL IZERO(LU5LIST,IALL_LU5)
*
         EIGAPR = EIG(ITER-1,IROOT)
*
*        calculate residues ...
*
         CALL P1_B_PAR_RL_1(VEC1,VEC2,AVEC,EIGAPR,RNRM,EIGSHF,
     &                      EIG,TEST,THRES_E,RTCNV,CONVER,ITER,MAXIT,
     &                      IROOT,LU2LIST,LU1LIST,LU5LIST,NBATV,LBATV,
     &                      LEBATV,I1BATV,IBATV,MY_LU2_OFF,MY_LU1_OFF,
     &                      MY_LU5_OFF,SCRRED,NVEC,ILU2,ILU1,ILU5)
*
         IF( ITER .GT. MAXIT) GOTO 100
*
*        new direction needed?
*
*        1.2 : multiply with inverse Hessian approximation 
*              to get new direction
*
         IF( .NOT. RTCNV(IROOT) ) THEN
*
*          (D-E)-1 *( HX - EX )
*
           IADD = IADD + 1
*
           CALL H0M1TD_REL_PAR(ILU3,IDIA,ILU5,-EIGAPR,VEC1,VEC2,
     &                         LU3LIST,LU5LIST,NBATV,LBATV,LEBATV,
     &                         I1BATV,IBATV,MY_LU3_OFF,MY_DIA_OFF,
     &                         MY_LU5_OFF,1,THRES_E)
*               H0M1TD_REL_PAR(LUOUT,LUDIA,LUIN,SHIFT,VEC1,VEC2,LISTOUT,
*    &                         LISTIN,NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
*    &                         OFFSET_OUT,OFFSET_DIAG,OFFSET_IN,INV,THRES_E)
*
*
*          hardwired to 1!
*         
           IF(IOLSTM .NE. 0 ) THEN
*
*            add Olsen correction if neccessary
*            current eigen-vector on LU5
*
             CALL IZERO(LU4LIST,IALL_LU4)
             CALL IZERO(LU5LIST,IALL_LU5)
*
*
             CALL P1_B_PAR_RL_2(VEC1,VEC2,AVEC,-EIGAPR,IROOT,LU1LIST,
     &                          LU4LIST,LU5LIST,LU3LIST,
     &                          NBATV,LBATV,LEBATV,
     &                          I1BATV,IBATV,MY_LU1_OFF,MY_LU4_OFF,
     &                          MY_LU5_OFF,MY_LU3_OFF,MY_DIA_OFF,
     &                          NVEC,ILU1,ILU4,ILU5,ILU3,IDIA,1)
*
           END IF
*
*          1.3 orthogonalize to all previous vectors
*          1.4 normalize vector
*
           CALL P1_B_PAR_RL_3(VEC1,VEC2,WORK,LU1LIST,LU3LIST,
     &                        LU1LIST,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU1_OFF,MY_LU3_OFF,MY_LU1_OFF,
     &                        SCRRED,NVEC,IADD,ILU1,ILU3,ILU1)
*
*
         END IF
*        ^ converged?
  100 CONTINUE
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      IF( TIMING )
     &WRITE(LUWRT,9250) WALLTID
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
*===========================================================
*                       PART 2                             = 
*                                                          =
*         Optimal combination of new and old directions    =
*                                                          =
*===========================================================
*
*     2.1: multiply new directions with matrix
*
      starttime = MPI_WTIME()
      xixidletime = 0.0D0
*

      CALL DZERO(SCRRED,LZERO_SCRRED)
      IMUSTRED = 0
*
      DO 150 IVEC = 1, IADD
*
*       copy c-vector to working-file ILUC
*
*
        xidletime = MPI_WTIME()
        CALL MPI_BARRIER(MYNEW_COMM,IERR)
        xixidletime = xixidletime - xidletime + MPI_WTIME()
*
CSK        WRITE(LUWRT,*) 'LU1LIST again:'
CSK        CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
*
*       reset LUCLIST
*
        CALL IZERO(LUCLIST,MY_ACT_BLK_ALL)
*
        CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                      IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                      LUCLIST,IBLOCKL,NVEC+IVEC-1)
*
CSK         CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
CSK        WRITE(LUWRT,*) 'LUCLIST for the 2nd time:'
CSK        CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUWRT)
CSK        WRITE(LUWRT,*) 'LU1LIST for the 2nd time:'
CSK        CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
*
*       set offset for sigma-file
*
        JVEC_SF = NVEC + IVEC - 1
*
*       start calculation: sigma = H x C
*
*       timing this sigma-vector computation
        sigmatime = MPI_WTIME()
*
*=======================================================================
        call sigden_ctrl(VEC1,VEC2,ILUC,ILU2,T_CC,T_BUFF,1
#if defined (VAR_MPI2)
     &                  ,LUCLIST,LU2LIST,IBLOCKL,NPARBLOCK,
     &                   IGROUPLIST,IPROCLIST,RCCTOS,
     &                   IT_TTPL,IT_TTOL
     &                  ,IBLOCKL,NPARBLOCK
#endif
     &                  )
*=======================================================================
*
C       end of timing
        sigmatime2 = MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
        WRITE(LUWRT,9400) WALLTID
*
*       augment projected matrix using batch structure of CI vector(s)
*
        CALL INPROD_B_PAR_RL2(ILU2,ILU1,VEC1,VEC2,SCRRED,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU2_OFF,MY_LU1_OFF,LU2LIST,
     &                        LU1LIST,IVEC,NVEC,IMUSTRED,ISTRED)
*
  150 CONTINUE
*
*      synchronize MPI_COMM_WORLD
*
       xidletime = MPI_WTIME()
       CALL DZERO(APROJ(ISTRED),IMUSTRED)
       CALL REDVEC_REL(SCRRED(ISTRED),APROJ(ISTRED),IMUSTRED,2,
     &                 MPI_SUM,MPI_COMM_WORLD,-1)
       xixidletime = xixidletime - xidletime + MPI_WTIME()
*
*     2.2: diagonalize projected matrix
*
      NVEC = NVEC + IADD
      CALL DCOPY(NVEC*(NVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
      CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
*
      endtime      = MPI_WTIME()
      tottime      = 0.0D0
      tottime_save = 0.0D0
      tottime      = endtime - starttime
C     TIMING FOR PARTS 2.1 - 2.2
      tottime_save = tottime
      WPART22 = SECTID(tottime)
      WRITE(LUWRT,9350) WPART22
*
      IF( TIMING )THEN
*
        xixidletime_save = 0.0D0
        xixidletime_save = xixidletime
*
        WALLTID = SECTID(xixidletime)
*
*       print idle time
*
        WRITE(LUWRT,'(/A,1X,A)')
     &  ' accumulated idle time in part 2                 :',WALLTID
        xpercent = (xixidletime_save/tottime_save) * 100
        WRITE(LUWRT,'(A,F14.9,A/)')
     &  ' ratio (idle time)/(time part 2) =',xpercent,' %'
*
      END IF
*
*     store new eigenvalues
*
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
*      debug printings
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
        CALL WRTMATMN(AVEC,NVEC,NROOT,MAXVEC,NROOT,LUWRT)
      END IF
*
*
  901 CONTINUE
*
*===========================================================
*                         PART 3                           =
*                                                          =
*      perhaps reset or assemble converged eigenvectors    =
*                                                          =
*===========================================================
*
      IF( NVEC+NROOT .gt. MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER )THEN
*
*     check timing
      timer3 = 0.0D0
      starttimer = MPI_WTIME()
*
*
CSK        CALL DZERO(SCRRED,NROOT)
        CALL IZERO(LU3LIST,IALL_LU3)
*
*       c vectors to ILU1
*
        DO IROOT = 1, NROOT
*
          CALL P3_B_PAR_RL_1(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
     &                       LU1LIST,LU3LIST,
     &                       NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                       MY_LU1_OFF,MY_LU3_OFF,
     &                       NVEC,IROOT,ILU1,ILU3)
*
        END DO
*
*       update WORK array to get correct scaling factor
*
        CALL DZERO(WORK,NROOT)
*
*       no scaling, we should already work in a normalized basis
*
        CALL SETVEC(WORK,1.0D0,NROOT)
        CALL IZERO(LU1LIST,NROOT*MY_ACT_BLK1)
*
        DO IROOT = 1, NROOT
*
          CALL COPVCD_PP_B_RL(VEC1,LU3LIST,LU1LIST,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU3_OFF,MY_LU1_OFF,IROOT,ILU3,ILU1)
*
        END DO
*
        CALL IZERO(LU3LIST,IALL_LU3)
*
*
*       corresponding sigma vectors to ILU2
*
        DO IROOT = 1, NROOT
*
          CALL P3_B_PAR_RL_1(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
     &                       LU2LIST,LU3LIST,
     &                       NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                       MY_LU2_OFF,MY_LU3_OFF,
     &                       NVEC,IROOT,ILU2,ILU3)
*
        END DO
*
!       fix for memory crash on "idle" co-workers
!       CALL IZERO(LU2LIST,NROOT)
        CALL IZERO(LU2LIST,NROOT*MY_ACT_BLK2)
*
        DO IROOT = 1, NROOT
*
          CALL COPVCD_PP_B_RL(VEC1,LU3LIST,LU2LIST,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU3_OFF,MY_LU2_OFF,IROOT,ILU3,ILU2)
*
        END DO
*
*       reset NVEC
*
        NVEC = NROOT
*
*       reset subspace matrices ...
*
        CALL DZERO(AVEC,NVEC**2)
        CALL DZERO(APROJ,NVEC*(NVEC+1)/2)
*
        DO IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
        END DO
C
C       write to check point file KRCI_CVECS.x (unit 61 see io_r.F)
        IF( (.NOT. CONVER) .or. (ITER .lt. MAXIT))THEN
          IF( CHECKPOINT_LUCIX .and. (ITER_CHECKP.ge.ICHPARAM))THEN
C
C           reset ITER_CHECKP
            ITER_CHECKP = 0
C           scratch file LUSC41 (unit 79 see io_r.F)
            CALL REWINE(79,-1)
            DO JROOT = 1, NROOT
              CALL COPVCD_PAR_BDRIV_REL(ILU1,79,VEC1,
     &                                  NPARBLOCK,NUM_BLOCKS,
     &                                  IBLOCKL,MPI_COMM_WORLD,-1,
     &                                  JROOT,LU1LIST,MY_LU1_OFF,1)
            END DO
            IF( MYPROC .eq. MASTER )THEN
              CALL REWINE(79,-1)
              CALL REWINE(61,-1)
              DO JROOT = 1, NROOT
                CALL COPVCDC(79,61,VEC1,0,1,-1)
              ENDDO
              CALL REWINE(61,-1)
            END IF
          END IF
        END IF
C
C       finish timing part3
        timer3 = timer3 + MPI_WTIME() - starttimer
        WALLTSTEP = SECTID(timer3)
        WRITE(LUWRT,9600) WALLTSTEP
      END IF
*     ^ NVEC+NROOT > MAXVEC or CONVER == .TRUE. or ITER == MAXIT
*
*     timing of this iteration
*
      WALLITR2 = MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,9300) WALLTID
*
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
C
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         WRITE(LUWRT,1170) MAXIT
 1170    FORMAT(/,'  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         WRITE(LUWRT,1180) ITER
 1180    FORMAT(/,'  Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 0 ) THEN
        DO 1600 IROOT = 1, NROOT
          WRITE(LUWRT,*)
          WRITE(LUWRT,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(LUWRT,*)
     &    '============================================'
          WRITE(LUWRT,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          IF (RTCNV(IROOT)) THEN
             WRITE(LUWRT,1190) IROOT,FINEIG(IROOT)+EIGSHF
          ELSE
             WRITE(LUWRT,1191) IROOT,FINEIG(IROOT)+EIGSHF
          END IF
 1190 FORMAT(' The final eigenvalue',I5,F22.10,' (converged)')
 1191 FORMAT(' The final eigenvalue',I5,F22.10,' (NOT converged)')
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
        write(LUWRT,'(/A/)')
     &  '  *++++++++++++++++++++++++++++++++++++++++++++*'
        DO 1601 IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
          IF (RTCNV(IROOT)) THEN
             WRITE(LUWRT,1190) IROOT,FINEIG(IROOT)
          ELSE
             WRITE(LUWRT,1191) IROOT,FINEIG(IROOT)
          END IF
 1601   CONTINUE
        write(LUWRT,'(/A/)')
     &  '  *++++++++++++++++++++++++++++++++++++++++++++*'
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(LUWRT,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
C     store information about energies for each root on file
C     KRCI_CVECS.INFO (fh LU_INFO) - used in property calculations
      IF( MYPROC .eq. MASTER ) THEN
        DO I = 1, NROOT
           FINEIG(I) = EIG(ITER,I)+EIGSHF
        END DO
        WRITE(FILELAB,'(A5,A3)') "eroot",SYMFLABEL
        CALL KRCI_PRPFILE(LU_INFO,FILELAB,FINEIG,NROOT,0)
      END IF
*
      call qexit('MICDV4')
*
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))

 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION            : ',A)
 9250 FORMAT(' >>>  WALL TIME IN STEP 1 OF CURRENT ITERATION   : ',A)
 9400 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A/)
 9600 FORMAT(' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A)
 9350 FORMAT(' >>>  WALL TIME FOR PART 2.1 - 2.2               : ',A)
 9777 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9888 FORMAT(' >>>  WALL TIME FOR INITIAL ITERATION            : ',A)
*
      END
***********************************************************************

      SUBROUTINE CMICDV_PAR(VEC1,VEC2,RNRM,EIG,FINEIG,
     &                      MAXIT,NVAR,NROOT,MAXVEC,NINVEC,T_CC,T_BUFF,
     &                      APROJI,APROJR,AVEC,WORK,IPRT,
     &                      NPRDIM,IPNTR,LBLK,EIGSHF,
     &                      IHAM12,SIGDEN_ROUTE,
     &                      CONVLIM,IBLOCKL,NPARBLOCK,NBLOCKS,
     &                      LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                      LUCLIST,NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      RCCTOS,SCRRED_R,SCRRED_I,IGROUPLIST,
     &                      IPROCLIST,IT_TTPL,IT_TTOL,AMATC,CWORK,
     &                      RWORK,EVL,IWORK)
*
* Multi root Davidson routine for complex Hermitian
* matrices, requires two blocks in core
*
* May 1994, Jeppe Olsen
*
* Initial version without general preconditioner
* or inverse iteration corrections
*
* General revision for parallel implementation - CMICDV_PAR
*
* June 2007, Stefan Knecht
*
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
*
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
*
* The vectors are in general multi record vectors.
* The real and imaginary parts are separated by end of vector marks
* (-1) in my simple world
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      REAL * 8   INPROD,INPRDD
      character*6 SIGDEN_ROUTE
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12, WALLTSTEP*12
      CHARACTER WPART22*12, FILELAB*8
      LOGICAL CONVER,RTCNV(100),EXSTOP
*
#include "krmcluci_inf.h"
*
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
      DIMENSION SCRRED_R(*), SCRRED_I(*)
      INTEGER   RCCTOS(*)
      DIMENSION IGROUPLIST(*), IPROCLIST(*)
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION RNRM(NROOT,MAXIT),EIG(2*MAXVEC,MAXIT)
      DIMENSION APROJI(*),APROJR(*),AVEC(*),WORK(*)
      DIMENSION IPNTR(1), IBLOCKL(*), NPARBLOCK(*)
      DIMENSION IT_TTPL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL(*)
      COMPLEX*16 AMATC(MAXVEC,*), CWORK(*)
      DIMENSION RWORK(*),EVL(*)
      INTEGER   IWORK(*)
C
C     min. DIMENSION of matrices and arrys required for complex 
C     hermitian matrix diagonalizer
C
C     AMATC: MAXVEC**2
C     CWORK: 2*MAXVEC + MAXVEC**2
C     RWORK: 1 + 5*MAXVEC + 2*MAXVEC**2    
C     EVL  : 1*MAXVEC
C     IWORK: 3 + 5*MAXVEC
C
*
* Dimensioning required of local vectors
*      APROJR : MAXVEC*(MAXVEC+1)/2
*      APROJI : MAXVEC*(MAXVEC+1)/2
*      AVEC   : (2*MAXVEC) ** 2
*      WORK   : 2*MAXVEC*(2*MAXVEC+1)/2
*
      DIMENSION FINEIG(*)
*
      call qenter('CMICDV_PAR')
*
      NTEST = 0000
      IPRT = MAX(NTEST,IPRT)
C
C     no screening 
      CSCREEN = .FALSE.
C
C     timing
      TIMING = .FALSE.
C
C     set check point parameters - approx after every 5th iteration
      ITER_CHECKP = 0
      ICHPARAM    = 0
      IF( CHECKPOINT_LUCIX )THEN
        ICHPARAM    = 5
      END IF
C
C     initialize LZERO_SCRRED
      LZERO_SCRRED = 0
      LZERO_SCRRED = MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2 )
C
C
      WRITE(LUWRT,'(/7X,A)') 
     & ' ==========================================='
      WRITE(LUWRT,'(7X, A)') 
     & ' start of complex Davidson solver CMICDV_PAR'
      WRITE(LUWRT,'(7X,A/)') 
     & ' ==========================================='
C
      WRITE(LUWRT,'(/7X,A)') 
     &              '=============================================='
      WRITE(LUWRT,'(7X,A)') ' I will take you to a place '
      WRITE(LUWRT,'(7X,A)') ' Where we shall find our'
      WRITE(LUWRT,'(7X,A)') '  ROOTS, BLOODY ROOTS '
      WRITE(LUWRT,'(7X,A)')
      WRITE(LUWRT,'(7X,A)') '               Sepultura '
      WRITE(LUWRT,'(7X,A/)') 
     &              '=============================================='
*
      WRITE(LUWRT,'(/A,1F25.15/)') 
     &              '  Nuclear repulsion + core energy: ',EIGSHF
*
      IF( MAXVEC .LT. 2 * NROOT ) THEN
        WRITE(LUWRT,*) ' Sorry CMICDV wounded , MAXVEC .LT. 2*NROOT '
        WRITE(LUWRT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
        WRITE(LUWRT,*) ' Increase MXCIV_CI to be at least 2 * Nroot '
        Call Abend1(' Enforced stop in CMICDV ')
      END IF
*
*     Complex diagonalizer, therefore
      IRC = 2
*
      CONVER = .FALSE.
*
      IF( MAXIT .le. 0 ) GOTO 9999
*
      NZERO = 0
*
*     ===================
*      Initial iteration
*     ===================
      ITER_CHECKP = 1
      ITER = 1
      NVEC = 0
*     start timing of initial iteration
      xinit_time = MPI_WTIME()
*              
      CALL DZERO(APROJR, MAXVEC*(MAXVEC+1)/2 )
      CALL DZERO(APROJI, MAXVEC*(MAXVEC+1)/2 )
      CALL DZERO(SCRRED_R, MAXVEC*(MAXVEC+1)/2 )
      CALL DZERO(SCRRED_I, MAXVEC*(MAXVEC+1)/2 )
*
      DO 10 IVEC = 1,NINVEC
*
*       copy c-vector to working-file ILUC
        CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
*       reset LUCLIST; MY_ACT_BLK_ALL = NUMBLOCKS!
*
        CALL IZERO(LUCLIST,NUM_BLOCKS2)
*
        CALL COPVCD_PP_CC_B_C(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                        IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                        LUCLIST,IBLOCKL,IVEC-1)
*
*       set offset for sigma-file
*
        JVEC_SF = IVEC - 1
*
*       start calculation: sigma = H x C
*
*       timing this sigma-vector computation
*
        sigmatime = MPI_WTIME()
*=======================================================================
        call sigden_ctrl(VEC1,VEC2,ILUC,ILU2,T_CC,T_BUFF,1
#if defined (VAR_MPI2)
     &                  ,LUCLIST,LU2LIST,IBLOCKL,NPARBLOCK,
     &                   IGROUPLIST,IPROCLIST,RCCTOS,
     &                   IT_TTPL,IT_TTOL
     &                  ,IBLOCKL,NPARBLOCK
#endif
     &                  )
*=======================================================================
*       end of timing
        sigmatime2 = MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
        WRITE(LUWRT,9777) WALLTID
        WRITE(LUWRT,*) '   '
*
*       projected matrix using batch structure of CI vector(s)
*
        CALL INPROD_B_PAR_CPLX(ILU1,ILU2,VEC1,VEC2,APROJR,APROJI,
     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                         MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                         LU2LIST,IVEC)
*
 10   CONTINUE

*     timing of initial iteration
*
      xinit_time2 = MPI_WTIME()
      WALLTID = SECTID(xinit_time2-xinit_time)
      WRITE(LUWRT,9888) WALLTID
*
*     synchronize MPI_COMM_WORLD
*
      CONE = - 1.0D0
      IREDL = 0
      IREDL = NINVEC*(NINVEC-1)/2 + NINVEC
      CALL DZERO(SCRRED_R,IREDL)
      CALL DZERO(SCRRED_I,IREDL)
*     ... sum up
      CALL REDVEC_REL(APROJR,SCRRED_R,IREDL,2,MPI_SUM,
     &                MPI_COMM_WORLD,-1)
      CALL REDVEC_REL(APROJI,SCRRED_I,IREDL,2,MPI_SUM,
     &                MPI_COMM_WORLD,-1)
*     ... update APROJR and APROJI 
      CALL DCOPY(IREDL,SCRRED_R,1,APROJR,1)
      CALL DSCAL(IREDL,CONE,SCRRED_I,1)
      CALL DCOPY(IREDL,SCRRED_I,1,APROJI,1)
*
      IPRT = 00
*
      IF( IPRT .GE.10) THEN
        WRITE(LUWRT,*) ' Initial projected matrix, real part'
        CALL PRSYMMN(APROJR,NINVEC,LUWRT)
        WRITE(LUWRT,*) ' Initial projected matrix, imag part'
        CALL PRSYMMN(APROJI,NINVEC,LUWRT)
      END IF
C
C     diagonalize initial projected matrix
C
C     S. Knecht - Nov 2008
C     disabled old diagonalizer routine which gives wrong eigenvectors 
C     when we are dealing with (two) states of exact degeneracy.
C     I have now implemented an interface to the complex hermitian 
C     matrix diagonalizer ZHEEVD (see e.g. pdpack/gp_zlapack.F) that 
C     is nowadays part of Intels MKL /Lapack libraries.
CSK   
CSK   CALL CEIGEN(APROJR,APROJI,NINVEC,EIG(1,1),AVEC,WORK)
CSK
      CALL CEIGEN_NEW(APROJR,APROJI,NINVEC,EIG(1,1),AVEC,AMATC,CWORK,
     &                RWORK,EVL,IWORK)
C
      IF(IPRT .GE. 3) THEN
        WRITE(LUWRT,'(A,I4)') ' Eigenvalues of initial iteration '
        WRITE(LUWRT,'(F22.13)')
     &  ( EIG(IROOT,1)+EIGSHF,IROOT=1,NROOT)
      END IF
      IPRT = 0
      NVEC = NINVEC
*
      IF ( MAXIT .eq. 1) GOTO  901
*
*     ======================
*     Loop over iterations
*     ======================
*
 1000 CONTINUE
*       start timing of iteration
        WALLITR1 = MPI_WTIME()
*
        write(LUWRT,'(/A21,3X,I3)') ' Info from iteration ',ITER
        write(LUWRT,*) '_______________________'
*
*
*===========================================================
*                       PART 1                             =
*                                                          =
*              New directions to be included               =
*                                                          =
*===========================================================
*
*      1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       DO 100 IROOT = 1, NROOT
*
*        reset scratch file lists ...
*
         CALL IZERO(LU3LIST,IALL_LU3)
         CALL IZERO(LU4LIST,IALL_LU4)
         CALL IZERO(LU5LIST,IALL_LU5)
*
         EIGAPR = EIG(IROOT,ITER)
C        EIGAPR = EIG(IROOT,ITER-1)
*
*        calculate residues ...
*
         CALL P1_B_PAR_CX_1(VEC1,VEC2,AVEC,WORK,EIGAPR,RNRM,EIGSHF,
     &                      RTCNV,CONVER,CONVLIM,ITER,NROOT,MAXIT,
     &                      IROOT,LU1LIST,LU2LIST,LU5LIST,NBATV,LBATV,
     &                      LEBATV,I1BATV,IBATV,MY_LU1_OFF,MY_LU2_OFF,
     &                      MY_LU5_OFF,SCRRED_R,SCRRED_I,NVEC,2*NVEC,
     &                      ILU1,ILU2,ILU5)
*
*
*
*        I suppose this is only a check for the necessity of defining
*        a new search direction. Using intermediate ITERIM instead ...
CTF      ITER = ITER + 1
*
         ITERIM = ITER + 1
*
         IF( ITERIM .GT. MAXIT) GOTO 100
*
*        new direction needed?
*
*        1.2 : multiply with inverse Hessian approximation
*              to get new direction
*
         IF( .NOT. RTCNV(IROOT) ) THEN
*
*          (D-E)-1 *( HX - EX )
*
           IADD = IADD + 1
*
csk           WRITE(LUWRT,*) 'THIS IS MY LU5LIST after part1'
csk           CALL IWRTMAMN(LU5LIST,1,IALL_LU5,1,IALL_LU5,LUWRT)
*
           CALL H0M1TD_REL_PAR_CX(ILU3,IDIA,ILU5,-EIGAPR,VEC1,VEC2,
     &                            LU3LIST,LU5LIST,NBATV,LBATV,LEBATV,
     &                            I1BATV,IBATV,MY_LU3_OFF,MY_DIA_OFF,
     &                            MY_LU5_OFF,1,IRC,XNORM)
*
           XNORMO = XNORM
CSK        WRITE(LUWRT,*) 'THIS IS XNORMO',XNORMO
*
*          check for linear dependency...
           CALL DZERO(SCRRED_R,2*(NVEC+IADD))
           CALL DZERO(SCRRED_I,2*(NVEC+IADD))
*
           CALL INPROD_B_PAR_CPLX2(ILU1,ILU3,VEC1,VEC2,SCRRED_R,
     &                             SCRRED_I,
     &                             NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                             MY_LU1_OFF,MY_LU3_OFF,LU1LIST,
     &                             LU3LIST,NVEC,IADD)
*
*          communication is always essential ...
           CALL DZERO(VEC1,NVEC+IADD-1)
           CALL DZERO(VEC2,NVEC+IADD-1)
*
           CALL REDVEC_REL(SCRRED_R,VEC1,NVEC+IADD-1,2,
     &                     MPI_SUM,MPI_COMM_WORLD,-1)
           CALL REDVEC_REL(SCRRED_I,VEC2,NVEC+IADD-1,2,
     &                     MPI_SUM,MPI_COMM_WORLD,-1)
*
*          transfer information to WORK array...
           CALL DCOPY(NVEC+IADD-1,VEC1,1,WORK(1),2)
           CALL DCOPY(NVEC+IADD-1,VEC2,1,WORK(2),2)
*
           XR = 0.0D0
           XI = 0.0D0
*
           DO JVEC = 1,NVEC+IADD-1
             XR = WORK(2*(JVEC-1)+1)
             XI = WORK(2*(JVEC-1)+2)
             XNORM = XNORM - XR**2 - XI**2
           END DO
*
           THRSLDP = 1.0D-10
           IF( XNORM / XNORMO .le. THRSLDP )THEN
             WRITE(LUWRT,*)
     &       ' CMICDV : Vector eliminated due to linear dependence '
             IADD = IADD -1
             GOTO 100
           END IF
*
           SCALE = 1.0D0/SQRT(XNORM)
CSK        WRITE(LUWRT,*) ' Scale, scale**2*xnorm',scale,
CSK  &                      scale**2*xnorm
*
           CALL DSCAL( 2*(NVEC+IADD-1), SCALE, WORK, 1 )
C           
C          1.3 and 1.4 orthogonalize to all previous vectors 
C          and normalize vector
           ONEM = -1.0D0
           CALL DSCAL(NVEC+IADD-1,ONEM,WORK(2),2)
C
           CALL P1_B_PAR_CX_2(VEC1,VEC2,WORK,SCALE,LU1LIST,LU3LIST,
     &                        LU1LIST,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU1_OFF,MY_LU3_OFF,MY_LU1_OFF,
     &                        NVEC,IADD,ILU1,ILU3,ILU1)
C
         END IF
C        ^ not converged...
  100 CONTINUE
CTF
      ITER        = ITER + 1
      ITER_CHECKP = ITER_CHECKP + 1
CTF
      IF( CONVER ) THEN
        ITER = ITER - 1
        GOTO  901
      END IF
*
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
*===========================================================
*                       PART 2                             =
*                                                          =
*         Optimal combination of new and old directions    =
*                                                          =
*===========================================================
*
*     2.1: multiply new directions with matrix
*
      starttime = MPI_WTIME()
*
      CALL DZERO(SCRRED_R,LZERO_SCRRED)
      CALL DZERO(SCRRED_I,LZERO_SCRRED)
      IMUSTRED = 0
      ISTRED   = 0
*
      DO 150 IVEC = 1, IADD
*
*       copy c-vector to working-file ILUC
*
*
        CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
csk     WRITE(LUWRT,*) 'LU1LIST again:'
csk     CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
*
*       reset LUCLIST
*
        CALL IZERO(LUCLIST,NUM_BLOCKS2)
*
        CALL COPVCD_PP_CC_B_C(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                        IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                        LUCLIST,IBLOCKL, NVEC + IVEC - 1 )
*
csk     WRITE(LUWRT,*) 'LUCLIST for the 2nd time:'
csk     CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUWRT)
csk     WRITE(LUWRT,*) 'LU1LIST for the 2nd time:'
csk     CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
*
*       set offset for sigma-file
        JVEC_SF = NVEC + IVEC - 1
*
*       start calculation: sigma = H x C
*
*       timing this sigma-vector computation
        sigmatime = MPI_WTIME()
*=======================================================================
        call sigden_ctrl(VEC1,VEC2,ILUC,ILU2,T_CC,T_BUFF,1
#if defined (VAR_MPI2)
     &                  ,LUCLIST,LU2LIST,IBLOCKL,NPARBLOCK,
     &                   IGROUPLIST,IPROCLIST,RCCTOS,
     &                   IT_TTPL,IT_TTOL
     &                  ,IBLOCKL,NPARBLOCK
#endif
     &                  )
*=======================================================================
*       end of timing
        sigmatime2 = MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
        WRITE(LUWRT,9400) WALLTID
*
*       projected matrix using batch structure of CI vector(s)
        CALL INPROD_B_PAR_CPLX3(ILU1,ILU2,VEC1,VEC2,SCRRED_R,SCRRED_I,
     &                          NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                          MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                          LU2LIST,NVEC,IVEC,IMUSTRED,ISTRED)
*
 150  CONTINUE 
*
*      synchronize MPI_COMM_WORLD
*
       CALL DZERO(APROJR(ISTRED),IMUSTRED)
       CALL DZERO(APROJI(ISTRED),IMUSTRED)
*
       CALL REDVEC_REL(SCRRED_R(ISTRED),APROJR(ISTRED),IMUSTRED,2,
     &                 MPI_SUM,MPI_COMM_WORLD,-1)
       CALL REDVEC_REL(SCRRED_I(ISTRED),APROJI(ISTRED),IMUSTRED,2,
     &                 MPI_SUM,MPI_COMM_WORLD,-1)
*
      endtime = MPI_WTIME()
      tottime = 0.D0
      tottime = endtime - starttime
*     TIMING FOR PARTS 2.1 - 2.2
      WPART22 = SECTID(tottime)
      WRITE(LUWRT,9350) WPART22
*
      IPRT = 00
*
*     diagonalize projected matrix
      NVEC = NVEC + IADD
      IF( IPRT  .GE. 10 ) THEN
        WRITE(LUWRT,*) ' Projected matrices before CEIGEN_new '
        CALL PRSYM(APROJR,NVEC,LUWRT)
        CALL PRSYM(APROJI,NVEC,LUWRT)
      END IF
C
C     S. Knecht - Nov 2008
C     disabled old diagonalizer routine which gives wrong eigenvectors
C     when we are dealing with (two) states of exact degeneracy.
C     I have now implemented an interface to the complex hermitian
C     matrix diagonalizer ZHEEVD (see e.g. pdpack/gp_zlapack.F) that
C     is nowadays part of Intels MKL /Lapack libraries.
CSK
CSK   CALL CEIGEN(APROJR,APROJI,NVEC,EIG(1,ITER),AVEC,WORK)
CSK
      CALL CEIGEN_NEW(APROJR,APROJI,NVEC,EIG(1,ITER),AVEC,AMATC,CWORK,
     &                RWORK,EVL,IWORK)
C
      IF(IPRT .GE. 3 ) THEN
        WRITE(LUWRT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(LUWRT,'(5F22.13)')
     &  ( EIG(IROOT,ITER)+EIGSHF,IROOT=1,NROOT)
      END IF
C
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRT,*) ' Projected matrix and eigen pairs '
        CALL PRSYMMN(APROJR,NVEC,LUWRT)
        CALL PRSYMMN(APROJI,NVEC,LUWRT)
        CALL WRTMATMN(AVEC,2*NVEC,NROOT,2*NVEC,NROOT,LUWRT)
      END IF
      IPRT = 0
C
  901 CONTINUE
C
C===========================================================
C                         PART 3                           =
C                                                          =
C      perhaps reset or assemble converged eigenvectors    =
C                                                          =
C===========================================================
C
      IF( NVEC + NROOT .gt. MAXVEC .or. CONVER .or. MAXIT .eq. ITER)THEN
*
*       check timing
        timer3 = 0.0D0
        starttimer = MPI_WTIME()
*
        CALL IZERO(LU3LIST,IALL_LU3)
*
*       c vectors from ILU1 to ILU1
        DO IROOT = 1, NROOT
          DO JVEC = 1, NVEC
            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
          END DO
          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU1LIST,LU3LIST,
     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      MY_LU1_OFF,MY_LU3_OFF,
     &                      NVEC,2*NVEC,IROOT,ILU1,ILU3)
        END DO
*
*       no scaling, we should already work in a normalized basis
        CALL IZERO(LU1LIST,NROOT*MY_ACT_BLK2)
*
        DO IROOT = 1, NROOT
          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU1LIST,
     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                         MY_LU3_OFF,MY_LU1_OFF,IROOT,ILU3,ILU1)
          WORK(2*NROOT*NVEC+IROOT) = 1.0D0
        END DO
*
        CALL IZERO(LU3LIST,IALL_LU3)
*
*       corresponding sigma vectors from ILU2 to ILU2
        DO IROOT = 1, NROOT
          DO JVEC = 1, NVEC
            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
          END DO
          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU2LIST,LU3LIST,
     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      MY_LU2_OFF,MY_LU3_OFF,
     &                      NVEC,2*NVEC,IROOT,ILU2,ILU3)
        END DO
*
*       no scaling, we should already work in a normalized basis
        CALL IZERO(LU2LIST,NROOT*MY_ACT_BLK2)
*
        DO IROOT = 1, NROOT
          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU2LIST,
     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                         MY_LU3_OFF,MY_LU2_OFF,IROOT,ILU3,ILU2)
        END DO
*
        CALL IZERO(LU3LIST,IALL_LU3)
*
        NVEC = NROOT
*
*       reset subspace matrices ...
*
        CALL DZERO(AVEC,2*NVEC**2)
        CALL DZERO(APROJR,NVEC*(NVEC+1)/2)
        CALL DZERO(APROJI,NVEC*(NVEC+1)/2)
*
        DO IROOT = 1,NROOT
          AVEC((IROOT-1)*2*NROOT+IROOT) = 1.0D0
          APROJR(IROOT*(IROOT+1)/2 )    = EIG(IROOT,ITER)
        END DO
C
C       write to check point file KRCI_CVECS.x (unit 61 see io_r.F)
        IF( (.NOT. CONVER) .or. (ITER .lt. MAXIT))THEN
          IF( CHECKPOINT_LUCIX .and. (ITER_CHECKP.ge.ICHPARAM))THEN
C
C           reset ITER_CHECKP
            ITER_CHECKP = 0
C           scratch file LUSC41 (unit 79 see io_r.F)
            CALL REWINE(79,-1)
            DO JROOT = 1, NROOT
              CALL COPVCD_PAR_BDRIV_REL(ILU1,79,VEC1,
     &                                  NPARBLOCK,NUM_BLOCKS,
     &                                  IBLOCKL,MPI_COMM_WORLD,LBLK,
     &                                  JROOT,LU1LIST,MY_LU1_OFF,IRC)
            END DO
            IF(MYPROC.EQ.MASTER) THEN
              CALL REWINE(79,-1)
              CALL REWINE(61,-1)
              DO JROOT = 1, NROOT
                CALL COPVCDC(79,61,VEC1,0,IRC,LBLK)
              ENDDO
              CALL REWINE(61,-1)
            END IF
          END IF
        END IF
C
C     finish timing part3
      timer3 = timer3 + MPI_WTIME() - starttimer
      WALLTSTEP = SECTID(timer3)
      WRITE(LUWRT,9600) WALLTSTEP
      END IF
*     ^ NVEC + NROOT > MAXVEC or CONVER == .TRUE. or ITER == MAXIT
*
*     timing of this iteration
      WALLITR2 = MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,9300) WALLTID
*
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 0 )
     &   WRITE(LUWRT,1170) MAXIT
 1170    FORMAT(/'   Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
C!?      ITER = ITER - 1
         IF (IPRT .GE. 0 )
     &   WRITE(LUWRT,1180) ITER
 1180    FORMAT(/'  Convergence was obtained in ',I3,' iterations')
      END IF
*
      if (IPRT.gt.0) then
        DO 1600 IROOT = 1, NROOT
          WRITE(LUWRT,'(/A,I3/A/)')
     &    ' Information about convergence for root... ' ,IROOT,
     &    ' ============================================'
          FINEIG(IROOT) = EIG(IROOT,ITER)
          IF (RTCNV(IROOT)) THEN
             WRITE(LUWRT,1190) IROOT,FINEIG(IROOT)+EIGSHF
          ELSE
             WRITE(LUWRT,1191) IROOT,FINEIG(IROOT)+EIGSHF
          END IF
 1190 FORMAT(' The final eigenvalue',I5,F22.10,' (converged)')
 1191 FORMAT(' The final eigenvalue',I5,F22.10,' (NOT converged)')
          WRITE(LUWRT,1300)
 1300     FORMAT(/'  Summary of iterations '
     +           /'  ----------------------')
          WRITE(LUWRT,1310)
 1310     FORMAT
     &    (/'  Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(LUWRT,1340) I,EIG(IROOT,I)+EIGSHF,RNRM(IROOT,I)
 1340     FORMAT(7X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      else
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
        do 1601 IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(IROOT,ITER)+EIGSHF
          IF (RTCNV(IROOT)) THEN
             WRITE(LUWRT,1190) IROOT,FINEIG(IROOT)
          ELSE
             WRITE(LUWRT,1191) IROOT,FINEIG(IROOT)
          END IF
 1601   continue
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
      end if
*
      ECI = EIG(1,ITER)
*
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(LUWRT,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(IROOT,ITER),
     &                 RNRM(IROOT,1)/RNRM(IROOT,ITER)
 1607   CONTINUE
      END IF
C
C     store information about energies for each root on file
C     KRCI_CVECS.INFO (fh LU_INFO) - used in property calculations
      IF( MYPROC .eq. MASTER ) THEN
        DO I = 1, NROOT
           FINEIG(I) = EIG(I,ITER)+EIGSHF
        END DO
        WRITE(FILELAB,'(A5,A3)') "eroot",SYMFLABEL
        CALL KRCI_PRPFILE(LU_INFO,FILELAB,FINEIG,NROOT,0)
      END IF
C
 9999 CALL QEXIT('CMICDV_PAR')
*
      RETURN
 1030 FORMAT(/3X,7F15.8,/,(3X,7F15.8))
 1120 FORMAT(/3X,I3,7F15.8,/,(6X,7F15.8))
*
 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION            : ',A/)
 9400 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A/)
 9600 FORMAT(' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A/)
 9350 FORMAT(' >>>  WALL TIME FOR PART 2.1 - 2.2               : ',A/)
 9777 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A/)
 9888 FORMAT(' >>>  WALL TIME FOR INITIAL ITERATION            : ',A/)
*
      END
#else
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PAR_DUMMY_REL2
C**********************************************************************
C     dummy routine for normal compilation
C**********************************************************************
      END
#endif
