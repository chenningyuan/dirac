!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

#if defined (VAR_MPI2)
***********************************************************************

      SUBROUTINE MICDV4_ENLMD_REL_PAR(T_CC,T_BUFF,VEC1,VEC2,
     &                  RNRM,EIG,FINEIG,MAXIT,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,AVECO,WORK,IPRT,
     &                  NPRDIM,IPNTR,LBLK,EIGSHF,
     &                  thres_G,thres_E,IBLOCKL,NPARBLOCK,NBLOCKS,
     &                  LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LU6LIST,LU7LIST,
     &                  LUCLIST,NBATV,
     &                  LBATV,LEBATV,I1BATV,IBATV,RCCTOS,SCRRED,
     &                  IGROUPLIST,IPROCLIST,IT_TTPL,IT_TTOL,
     &                  AMATC,RWORK,EVL,eciold,ITERSEOUT)

*     Davidson algorithm , requires two blocks in core
*     parallel Multi root version

*     orig. written by Jeppe Olsen Winter of 1991

*     updated to allow general preconditioner, October 1993

*     Version using H0 + Lambda V as Sigma routine

*     Modified version of MICDV4_ENLMD_REL_PAR as a 
*     generalized version of MICDV6. It accepts an (within memory
*     limits) arbitrarily large Davidson subspace, then reduces to
*     2*NRoot vectors when the space has been filled. Andreas Nyvang and
*     Jeppe Olsen, September 2019         

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

*           LBLK : Defines block structure of matrices
* On input ILU1 is supposed to hold initial guesses to eigenvectors


      use interface_to_mpi
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION LU6LIST(*), LU7LIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
#include "parluci.h"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXVEC,MAXIT)
       DIMENSION APROJ(*),AVEC(*),AVECO(*),WORK(*)
       DIMENSION IPNTR(1), IBLOCKL(*), NPARBLOCK(*)
       DIMENSION T_CC(*), SCRRED(*), T_BUFF(*)
       INTEGER   RCCTOS(*)
       DIMENSION IT_TTPL(*) 
       INTEGER(KIND=df_MPI_OFFSET_KIND) IT_TTOL(*)
       DIMENSION AMATC(MAXVEC,*), RWORK(*), EVL(*), eciold(*)

* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2                               
*      SCRRED : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)

       DIMENSION FINEIG(NROOT)
       LOGICAL CONVER,RTCNV(NROOT)
       REAL*8 INPRDD
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12, WALLTSTEP*12
       CHARACTER WPART22*12, FILELAB*8
       INTEGER NZERO,LZERO_SCRRED
       INTEGER(KIND=df_MPI_OFFSET_KIND) IOFF_SCR
       Integer(Kind=8) IError

       real*8, allocatable, dimension(:)  :: scrcno
       real*8, allocatable, dimension(:)  :: CNO
       real*8, allocatable, dimension(:)  :: XJEP
       integer*8, allocatable, dimension(:)  :: IXJEP
       logical :: lex,orth_twice
       character(len=80) :: dumchar

       call qenter('MICDV4')

      write(LUWRT,*)
      write(LUWRT,*) ' ----------------------------------------------'
      write(LUWRT,*) " I'll take you to a place"
      write(LUWRT,*) ' Where we shall find our'
      write(LUWRT,*) '  ROOTS, BLOODY ROOTS'
      write(LUWRT,*)
      write(LUWRT,*) '               Sepultura'
      write(LUWRT,*) ' ----------------------------------------------'
      write(LUWRT,*)

C      no screening 
       CSCREEN = .FALSE.

C      do more explicit timing
       TIMING = .false.

C      set check point parameters - approx after every 5th iteration 
C      of MAXIT
       ITER_CHECKP = 0
       ICHPARAM    = 0
       IF( CHECKPOINT_LUCIX )THEN
         ICHPARAM    = 5
       END IF


!tb disable IPICO and IRESET
       IPICO = 0
       IRESET = 0
       IROOTHOMING = 0
       IDONEWBAS = 1
      IF ( MYPROC .EQ. MASTER) THEN
       inquire(file='luciarelparam.dat',exist=lex)
       if ( lex ) then
        write(luwrt,*) ' Master found file luciarelparam.dat'
        open(12234,file='luciarelparam.dat') 
        read(12234,*) dumchar,IROOTHOMING
        read(12234,*) dumchar,IDONEWBAS
       close(12234)
       call interface_mpi_bcast(IROOTHOMING,1,MASTER,
     &                             global_communicator)
       call interface_mpi_bcast(IDONEWBAS,1,MASTER,
     &                             global_communicator)
       endif
      ENDIF
!       if ( IROOTHOMING.EQ.1 ) 
!     & write(LUWRT,*) ' IROOTHOMING switched on !'
!       if ( IDONEWBAS.EQ.1 ) 
!     & write(LUWRT,*) ' IDONEWBAS switched on !'






*      initialize LZERO_SCRRED
       LZERO_SCRRED = 0
       LZERO_SCRRED = MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2 )

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

       NVEC = 0

       KAPROJ = 1
       KFREE  = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       CONVER = .FALSE.
       NZERO  = 0

*      ===================
*       Initial iteration
*      ===================

*      start timing of initial iteration
       WALLITR1 = interface_MPI_WTIME()
        
       CALL DZERO(APROJ, MAXVEC*(MAXVEC+1)/2 )

       ITER        = 1
       ITER_CHECKP = 1

       DO 10 IVEC = 1,NINVEC

*        copy c-vector to working-file ILUC

         call interface_mpi_BARRIER(MYNEW_COMM)

*        reset LUCLIST; MY_ACT_BLK_ALL = NUMBLOCKS!
         CALL IZERO(LUCLIST,NUM_BLOCKS2)

         CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                       IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                       LUCLIST,IBLOCKL,IVEC-1)

*        set offset for sigma-file

         JVEC_SF = IVEC - 1

*        start calculation: sigma = H x C

*        timing this sigma-vector computation

         sigmatime = interface_MPI_WTIME()

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

*        end of timing
         sigmatime2 = interface_MPI_WTIME()
         WALLTID = SECTID(sigmatime2-sigmatime)

         if(timing)then
           WRITE(LUWRT,9777) WALLTID
           WRITE(LUWRT,*) '   '
         end if

*        projected matrix using batch structure of CI vector(s)

         CALL INPROD_B_PAR_RL(ILU1,ILU2,VEC1,VEC2,APROJ,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                        LU2LIST,IVEC)

   10  CONTINUE

*      timing of initial iteration

       WALLITR2 = interface_MPI_WTIME()
       WALLTID = SECTID(WALLITR2-WALLITR1)
       if(timing)then
         WRITE(LUWRT,9888) WALLTID
         WRITE(LUWRT,*) '   '
       end if

*      synchronize global_communicator

       IREDL = 0
       IREDL = NINVEC*(NINVEC-1)/2 + NINVEC
       CALL DZERO(SCRRED,IREDL)
       CAll redvec(APROJ,SCRRED,IREDL,2,op_MPI_SUM,
     &             global_communicator,-1)
       CALL DCOPY(IREDL,SCRRED,1,APROJ,1)


       IF( IPRT .GE.3 ) THEN
         WRITE(LUWRT,*) ' INITIAL PROJECTED MATRIX  '
         CALL WRTMATMN(APROJ,1,IREDL,1,IREDL,LUWRT)
       END IF

*      diagonalize initial projected matrix
       CALL DCOPY(NINVEC*(NINVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
!      CALL EIGEN(WORK(KAPROJ),AVEC,NINVEC,0,1)
!       CALL EIGEN_NEW(WORK(KAPROJ),NINVEC,EIG(1,1),AVEC,AMATC,
!     &                RWORK,EVL)
!tb  only master should copy in order to avoid phase-problems
      IF ( MYPROC .EQ. MASTER) THEN 
         
      CALL EIGEN_NEW(WORK(KAPROJ),NINVEC,EIG(1,1),AVEC,AMATC,
     &               RWORK,EVL)
      ENDIF
       call interface_mpi_bcast(EIG(1,1),NINVEC,MASTER,
     &                             global_communicator)
       call interface_mpi_bcast(AVEC(1),NINVEC**2,MASTER,
     &                             global_communicator)
!      DO 20 IROOT = 1, NROOT
!        EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
!  20  CONTINUE

       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRT,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(LUWRT,'(5F18.13)')
     &   ( EIG(IROOT,1)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(LUWRT,*) ' Initial set of eigen values (no shift) '
         CALL WRTMATMN(EIG(1,1),1,NROOT,NROOT,MAXIT,LUWRT)
       END IF

!       WRITE(LUWRT,*) 'NVEC, NINVEC' , NVEC, NINVEC
       NVEC = NINVEC
C       NVEC = NROOT

       IF (MAXIT .EQ. 1 ) GOTO  901

*      assuming diagonal preconditioner

       IPRECOND = 1

* ======================
*. Loop over iterations
* ======================

 1000 CONTINUE
*       start timing of iteration
        WALLITR1 = interface_MPI_WTIME()
        starttime = interface_MPI_WTIME()

        if (ITERSEOUT.eq.1) then
          write(LUWRT,'(/A,I6/A)') ' Iteration',ITER,
     &                             '_________________'
        else
        WRITE(LUWRT,'(//A,I4)')
     &        ' (MICDV4_ENLMD_REL_PAR) CI microiteration no.',ITER
        end if

        ITER        = ITER + 1
        ITER_CHECKP = ITER_CHECKP + 1


*===========================================================
*                       PART 1                             =
*                                                          =
*              New directions to be included               =
*                                                          =
*===========================================================

*      1.1 : R = H*X - EIGAPR*X

       IADD = 0
       CONVER = .TRUE.

       DO IROOT = 1, NROOT

*        reset scratch file lists ...

         CALL IZERO(LU3LIST,IALL_LU3)
         CALL IZERO(LU4LIST,IALL_LU4)
         CALL IZERO(LU5LIST,IALL_LU5)

         EIGAPR = EIG(IROOT,ITER-1)

*        calculate residues ...

         CALL P1_B_PAR_RL_1(VEC1,VEC2,AVEC,EIGAPR,RNRM,EIGSHF,
     &                      EIG,thres_G,thres_E,RTCNV,CONVER,ITER,MAXIT,
     &                      IROOT,LU2LIST,LU1LIST,LU5LIST,NBATV,LBATV,
     &                      LEBATV,I1BATV,IBATV,MY_LU2_OFF,MY_LU1_OFF,
     &                      MY_LU5_OFF,SCRRED,NVEC,ILU2,ILU1,ILU5,
     &                      MAXVEC)

         IF( ITER .GT. MAXIT) Cycle !IRoot

*        new direction needed?

*        1.2 : multiply with inverse Hessian approximation 
*              to get new direction

         IF( .NOT. RTCNV(IROOT) ) THEN

*          (D-E)-1 *( HX - EX )

           IADD = IADD + 1

           CALL H0M1TD_REL_PAR(ILU3,IDIA,ILU5,-EIGAPR,VEC1,VEC2,
     &                         LU3LIST,LU5LIST,NBATV,LBATV,LEBATV,
     &                         I1BATV,IBATV,MY_LU3_OFF,MY_DIA_OFF,
     &                         MY_LU5_OFF,1,thres_E)
*               H0M1TD_REL_PAR(LUOUT,LUDIA,LUIN,SHIFT,VEC1,VEC2,LISTOUT,
*    &                         LISTIN,NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
*    &                         OFFSET_OUT,OFFSET_DIAG,OFFSET_IN,INV,thres_E)


*          IOLSTM hardwired to 1!

           IF(IOLSTM .NE. 0 ) THEN

*            add Olsen correction if neccessary
*            current eigen-vector on LU5

             CALL IZERO(LU4LIST,IALL_LU4)
             CALL IZERO(LU5LIST,IALL_LU5)


             CALL P1_B_PAR_RL_2(VEC1,VEC2,AVEC,-EIGAPR,IROOT,LU1LIST,
     &                          LU4LIST,LU5LIST,LU3LIST,
     &                          NBATV,LBATV,LEBATV,
     &                          I1BATV,IBATV,MY_LU1_OFF,MY_LU4_OFF,
     &                          MY_LU5_OFF,MY_LU3_OFF,MY_DIA_OFF,
     &                          NVEC,ILU1,ILU4,ILU5,ILU3,IDIA,1)

           END IF

*          1.3 orthogonalize to all previous vectors
CAN        September 2019: Now always orthogonalizing twice
*          1.4 normalize vector

        start_orth = interface_MPI_WTIME()
           CALL P1_B_PAR_RL_3(VEC1,VEC2,WORK,LU1LIST,LU3LIST,
     &                        LU1LIST,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU1_OFF,MY_LU3_OFF,MY_LU1_OFF,
     &                        SCRRED,NVEC,IADD,ILU1,ILU3,ILU1)

        end_orth = interface_MPI_WTIME()

!      write(6,*) 'Time on orthogonalisation ',
!     &           SECTID(end_orth-start_orth)
         END IF
*        ^ converged?

      End do !Loop over IRoot

      endtime = interface_MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      IF( TIMING )WRITE(LUWRT,9250) WALLTID

      IF (IPRT .GE. 10) THEN
          WRITE(LUWRT,'(/A,F20.8)')
     &      '                        Core energy        =',EIGSHF
          WRITE(LUWRT,'(A,4F20.8,/,(37X,4F20.8))')
     &      '                        Active energy      =',
     &                        EIG(1:NROOT,ITER-1)
       END IF
      if (ITERSEOUT.eq.1) then
        do I=1,NROOT
          WRITE(LUWRT,'(A19,7X,I3,3X,1E18.13,3X,1F19.10)')
     &    ' Iter RNORM EIGAPR ',
     &      ITER-1,RNRM(ITER-1,I),EIG(I,ITER-1)+EIGSHF
        end do
      else
       WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &      '                        Total energy       =',
     &                        (EIG(I,ITER-1)+EIGSHF,i=1,NROOT)
       WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &      '                        Lowering of energy =',
     &      (eciold(i)-EIG(I,ITER-1)-EIGSHF,I=1,nroot)
       end if

       do i = 1, nroot
         eciold(i) = EIG(i,ITER-1) + EIGSHF
       end do

       if (ITERSEOUT.ne.1) then
       WRITE(LUWRT,'(/A,1P,D8.2,A/,(I15,D20.5))')
     &   '                        Norm of CI residuals'//
     &   ' (thr = ',thres_G,')',
     &   (i,RNRM(ITER-1,i),i=1,NROOT)
       end if
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF

*===========================================================
*                       PART 2                             = 
*                                                          =
*         Optimal combination of new and old directions    =
*                                                          =
*===========================================================

*     2.1: multiply new directions with matrix

      starttime = interface_MPI_WTIME()
      xixidletime = 0.0D0


      CALL DZERO(SCRRED,LZERO_SCRRED)
      IMUSTRED = 0

      DO IVEC = 1, IADD

*       copy c-vector to working-file ILUC


        xidletime = interface_MPI_WTIME()
        call interface_mpi_BARRIER(MYNEW_COMM)
        xixidletime = xixidletime - xidletime + interface_MPI_WTIME()

*       reset LUCLIST

        CALL IZERO(LUCLIST,MY_ACT_BLK_ALL)

        CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                      IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                      LUCLIST,IBLOCKL,NVEC+IVEC-1)

*       set offset for sigma-file

        JVEC_SF = NVEC + IVEC - 1

*       start calculation: sigma = H x C

*       timing this sigma-vector computation
        sigmatime = interface_MPI_WTIME()

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

C       end of timing
        sigmatime2 = interface_MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
!        write(6,*) 'Time in sigden_ctrl',WALLTID
        if (ITERSEOUT.eq.1) then
          if (IVEC.eq.1) WRITE(LUWRT,*)
          WRITE(LUWRT,9401) WALLTID
        else
        if(timing)then
          WRITE(LUWRT,9400) WALLTID
        end if
        end if

*       augment projected matrix using batch structure of CI vector(s)

        CALL INPROD_B_PAR_RL2(ILU2,ILU1,VEC1,VEC2,SCRRED,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU2_OFF,MY_LU1_OFF,LU2LIST,
     &                        LU1LIST,IVEC,NVEC,IMUSTRED,ISTRED)

      End do !Loop over IAdd new directions

*      synchronize global_communicator

       xidletime = interface_MPI_WTIME()
       CALL DZERO(APROJ(ISTRED),IMUSTRED)
       CAll redvec(SCRRED(ISTRED),APROJ(ISTRED),IMUSTRED,2,
     &                 op_MPI_SUM,global_communicator,-1)
       xixidletime = xixidletime - xidletime + interface_MPI_WTIME()

       nveco = nvec
       call dcopy(nveco*nveco,avec,1,aveco,1)

      !Sending eigenvector approximation AVEC from previous iteration
      Call Interface_mpi_bcast(NVeco,1,Master,global_communicator)
      call interface_mpi_bcast(AVECo(1:NVeco*NVeco),NVECO*NVECO,MASTER,
     &                             global_communicator)

*     2.2: diagonalize projected matrix

!      Write(LUWRT,*) ' NVEC changed', NVEC,' -> ', NVEC+IADD

      NVEC = NVEC + IADD
!        if(NVEC.GT.MAXVEC) NVEC = MAXVEC
      CALL DCOPY(NVEC*(NVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
!      CALL EIGEN_NEW(WORK(KAPROJ),NVEC,EIG(1,ITER),AVEC,AMATC,
!     &               RWORK,EVL)
!     CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
! only master should diagonalize in order to avoind phase-problems
      IF ( MYPROC .EQ. MASTER) THEN
         CALL EIGEN_NEW(WORK(KAPROJ),NVEC,EIG(1,ITER),AVEC,AMATC,
     &                  RWORK,EVL)
      ENDIF
       call interface_mpi_bcast(EIG(1,ITER),NVEC,MASTER,
     &                             global_communicator)
       call interface_mpi_bcast(AVEC(1),NVEC**2,MASTER,
     &                             global_communicator)

C      IF(IPICO.NE.0) THEN
C        E0VAR = WORK(KAPROJ)
C        C0VAR = AVEC(1)
C        C1VAR = AVEC(2)
C        C1NRM = SQRT(C0VAR**2+C1VAR**2)
C*. overwrite with pert solution
C        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
C        AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
C        E0PERT = AVEC(1)**2*APROJ(1)
C     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
C     &         + AVEC(2)**2*APROJ(3)
C        WORK(KAPROJ) = E0PERT
C        WRITE(LUWRT,*) ' Var and Pert solution, energy and coefficients'
C        WRITE(LUWRT,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
C        WRITE(LUWRT,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
C      END IF
!
C      IF(IROOTHOMING.EQ.1) THEN
C*
C*. Reorder roots so the NROOT with the largest overlap with
C*  the original roots become the first
C*
C       NNVEC=max(nvec,nroot)
C       allocate(xjep(nnvec**2),ixjep(nnvec))
C*. Norm of wavefunction in previous space
C       DO IVEC = 1, NVEC
C         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
C     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
C       END DO
C       WRITE(LUWRT,*)
C     & ' Norm of projections to previous vector space '
C       CALL WRTMAT(XJEP,1,NVEC,1,NVEC)
C*. My sorter arranges in increasing order, multiply with minus 1
C*  so the eigenvectors with largest overlap comes out first
C       ONEM = -1.0D0
C       CALL SCALVE(XJEP,ONEM,NVEC)
C       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
C       IF(NSORT.LT.NVEC) THEN
C         WRITE(LUWRT,*) ' Warning : Some elements lost in sorting '
C         WRITE(LUWRT,*) ' NVEC,NSORT = ', NSORT,NVEC
C       END IF
C       IF(IPRT.GE.3) THEN
C         WRITE(LUWRT,*) ' New roots choosen as vectors '
C         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
C       END IF
C*. Reorder
C       DO INEW = 1, NVEC
C         IOLD = IXJEP(INEW)
C         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
C       END DO
C       CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
C       DO INEW = 1, NVEC
C         IOLD = IXJEP(INEW)
C         XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
C       END DO
C       DO INEW = 1, NVEC
C         WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
C       END DO
C*
C       IF(IPRT.GE.3) THEN
C         WRITE(LUWRT,*) ' Reordered WORK and AVEC arrays '
C         CALL PRSYM(WORK,NVEC)
C         CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
C       END IF
C*
C       deallocate(xjep,ixjep)
C      END IF
*     ^ End of root homing procedure

      endtime      = interface_MPI_WTIME()
      tottime      = 0.0D0
      tottime_save = 0.0D0
      tottime      = endtime - starttime
C     TIMING FOR PARTS 2.1 - 2.2
      tottime_save = tottime
      WPART22 = SECTID(tottime)

      IF( TIMING )THEN
        WRITE(LUWRT,9350) WPART22

        xixidletime_save = 0.0D0
        xixidletime_save = xixidletime

        WALLTID = SECTID(xixidletime)

*       print idle time

        WRITE(LUWRT,'(/A,1X,A)')
     &  ' accumulated idle time in part 2                 :',WALLTID
        xpercent = (xixidletime_save/tottime_save) * 100
        WRITE(LUWRT,'(A,F14.9,A/)')
     &  ' ratio (idle time)/(time part 2) =',xpercent,' %'

      END IF

       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(LUWRT,'(5F18.13)')
     &   ( EIG(IROOT,ITER)+EIGSHF,IROOT=1,NROOT)
         WRITE(LUWRT,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(LUWRT,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF

      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRT,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRT,'(2X,E13.7)') (EIG(IROOT,ITER),IROOT = 1, NROOT)
        CALL WRTMATMN(AVEC,NVEC,NROOT,MAXVEC,NROOT,LUWRT)
      END IF


  901 CONTINUE !Jumps here in case of convergence
*
*===========================================================
*                         PART 3                           =
*                                                          =
*      perhaps reset or assemble converged eigenvectors    =
*                                                          =
*===========================================================
*

!      NADD = NROOT

      IF( NVEC+NROOT .gt. MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER )
     &THEN
!      IF( NVEC+NROOT .gt. MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER )THEN
!       IDONEWBAS set in top of subroutine



      If((NVEC+NROOT .gt. MAXVEC).and.(.not.CONVER)) then

      IF ( MYPROC .EQ. MASTER) THEN

         nnvec1=max(nroot,nvec)  
         allocate(cno(2*nnvec1**2))
         nnvec2=3*max(2*nroot,nvec)+1
         allocate(scrcno(MAX(4*2*NROOT*NVEC,4*(2*NROOT)**2)))
!         allocate(scrcno(3*nvec*nvec))
C         NVECUD=2*NROOT
         THR_ORTBAS = 1.0d-11
C         THR_ORTBAS = 1.0d-7
!         call get_cnewcold_bas_2(avec,nvec,nroot,nvecud,thr_ortbas)
!          write(6,*) 'nvec',nvec,'nroot',nroot,'nnvec1',nnvec1

          
        CALL GET_CNEWCOLD_BAS2(AVEC,AVECO,cno,
!        CALL GET_CNEWCOLD_BAS(AVEC(1:nvec*nroot),CNO(1:2*nnvec1**2),
     &                      NVEC,NVECO,NROOT,
     &                       scrcno(1),NVECUD,
     &                      IPRT,THR_ORTBAS,RTCNV)
!          write(6,*) 'nvec',nvec
!          write(6,*) 'NVECUD',NVECUD,'NADD',NVECUD-NROOT

         CALL COPVEC(CNO,AVEC,NVECUD*NVEC)

         deallocate(cno,scrcno)
      END IF
      Call Interface_mpi_bcast(NVecud,1,Master,global_communicator)
      call interface_mpi_bcast(AVEC(1:NVecUd*NVec),NVECUD*NVEC,MASTER,
     &                             global_communicator)
      Else
         NVecUd = NRoot
         Call Interface_mpi_bcast(NVecud,1,Master,global_communicator)


      end if



           
      !Sending new AVEC to all processes

C       call interface_mpi_BARRIER(global_communicator)




*     check timingG
      timer3 = 0.0D0
      starttimer = interface_MPI_WTIME()


CSK        CALL DZERO(SCRRED,NROOT)
C        CALL IZERO(LU3LIST,IALL_LU3)
C
C*       c vectors to ILU1
C
C        DO IROOT = 1, NROOT
C
C          CALL P3_B_PAR_RL_1(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
C     &                       LU1LIST,LU3LIST,
C     &                       NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                       MY_LU1_OFF,MY_LU3_OFF,
C     &                       NVEC,IROOT,ILU1,ILU3)
C
C        END DO
C
C*       update WORK array to get correct scaling factor
C
        CALL DZERO(WORK,NROOT)
C
C*       no scaling, we should already work in a normalized basis
C
        CALL SETVEC(WORK,1.0D0,NROOT)
C        CALL IZERO(LU1LIST,NROOT*MY_ACT_BLK1)
C
C        DO IROOT = 1, NROOT
C
C          CALL COPVCD_PP_B_RL(VEC1,LU3LIST,LU1LIST,
C     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                        MY_LU3_OFF,MY_LU1_OFF,IROOT,ILU3,ILU1)
C
C        END DO
C
C
C        CALL IZERO(LU3LIST,IALL_LU3)
C
C
C*       corresponding sigma vectors to ILU2
C
C        DO IROOT = 1, NROOT
C
C          CALL P3_B_PAR_RL_1(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
C     &                       LU2LIST,LU3LIST,
C     &                       NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                       MY_LU2_OFF,MY_LU3_OFF,
C     &                       NVEC,IROOT,ILU2,ILU3)
C
C        END DO
C
C!       fix for memory crash on "idle" co-workers
C!       CALL IZERO(LU2LIST,NROOT)
C        CALL IZERO(LU2LIST,NROOT*MY_ACT_BLK2)
C
C        DO IROOT = 1, NROOT
C
C          CALL COPVCD_PP_B_RL(VEC1,LU3LIST,LU2LIST,
C     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                        MY_LU3_OFF,MY_LU2_OFF,IROOT,ILU3,ILU2)
C
C        END DO


C     transform vectors to the actual eigenvector approximations: C and sigma


      CALL IZERO(LU3LIST,IALL_LU3)

      CALL TRAVC_B_RL_DRV(VEC1,VEC2,AVEC,LU1LIST,LU3LIST,
     &                    NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                    MY_LU1_OFF,MY_LU3_OFF,
C     &                    NVEC,NVECUD,ILU1,ILU3,IALL_LU1)
     &                    NVEC,NVECUD,ILU1,ILU3,NVecUd*MY_ACT_BLK1)

      CALL IZERO(LU3LIST,IALL_LU3)

      CALL TRAVC_B_RL_DRV(VEC1,VEC2,AVEC,LU2LIST,LU3LIST,
     &                    NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                    MY_LU2_OFF,MY_LU3_OFF,
C     &                    NVEC,NVECUD,ILU2,ILU3,IALL_LU2)
     &                    NVEC,NVECUD,ILU2,ILU3,NVecUd*MY_ACT_BLK2)

      CALL IZERO(LU3LIST,IALL_LU3)
      CALL IZERO(LU4LIST,IALL_LU4)

C          write(LUWRT,*) 'NVEC ', NVEC,' reset to ',nvecud
          NVEC = NVECUD

C       call interface_mpi_BARRIER(global_communicator)
      call interface_mpi_bcast(NVec,1,MASTER,
     &                             global_communicator)
!         NVEC=NROOT


*       reset subspace matrices ...

        CALL DZERO(AVEC,NVEC**2)
        CALL DZERO(APROJ,NVEC*(NVEC+1)/2)

C         CALL DZERO(AVEC,MAXVEC**2)
C         CALL DZERO(APROJ,MAXVEC*(MAXVEC+1)/2)


        DO JVEC = 1,NVEC
          AVEC(NVEC*(JVEC-1)+JVEC) = 1.0D0
C          APROJ(IROOT*(IROOT+1)/2 ) = EIG(IROOT,ITER)
        END DO

C      Call Interface_mpi_comm_rank(global_communicator,Iirank)
C      Call Interface_mpi_comm_size(global_communicator,Iisize)


       IREDL = 0
       IREDL = NVEC*(NVEC+1)/2
C       IREDL = NVEC*(NVEC-1)/2 + NVEC
C       write(6,*) 'IREDL ',IREDL
       CALL DZERO(SCRRED,IREDL)


         DO IVEC = 1,NVEC

         JVEC_SF = IVEC - 1


C         call interface_mpi_BARRIER(global_communicator)
         CALL INPROD_B_PAR_RL(ILU1,ILU2,VEC1,VEC2,APROJ,
     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                        MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                        LU2LIST,IVEC)

         ! computes APROJ=<ILU1,ILU2>

         ! check which communicator is working here


         END DO







C       call interface_mpi_BARRIER(global_communicator)
C         write(6,*) 'After barrier'
C         write(6,*) 'iredl ',iredl
C       call MPI_REDUCE(APROJ(1:IREDL), SCRRED(1:IREDL), IREDL, MPI_REAL8
C     &                 ,op_MPI_Sum, Master, global_communicator, IERROR)
C         write(6,*) 'After reduce'
C       call MPI_BCAST(ScrRed, IREDL, MPI_REAL8, Master,
C     &                global_communicator, IERROR)
C       write(6,*) 'After broadcast'

C         Call Interface_Mpi_AllReduce_R1_Work_F77
C     &       (AProj,ScrRed,IRedL,op_MPI_SUM,global_communicator)
       CAll redvec(APROJ(1),SCRRED(1),IREDL,2,op_MPI_SUM,
     &             global_communicator,-1)

       CALL DCOPY(IREDL,SCRRED,1,APROJ,1)

       IF( IPRT .GE.3 ) THEN
         WRITE(LUWRT,*) ' Reset projected matrix  '
         CALL WRTMATMN(APROJ,1,IREDL,1,IREDL,LUWRT)
       END IF


C       write to check point file KRCI_CVECS.x (unit 61 see io_r.F)
        IF( (.NOT. CONVER) .or. (ITER .lt. MAXIT))THEN
          IF( CHECKPOINT_LUCIX .and. (ITER_CHECKP.ge.ICHPARAM))THEN

C           reset ITER_CHECKP
            ITER_CHECKP = 0
C           scratch file LUSC41 (unit 79 see io_r.F)
            CALL REWINE(79,-1)
            DO JROOT = 1, NROOT
              CALL COPVCD_PAR_BDRIV_REL(ILU1,79,VEC1,
     &                                  NPARBLOCK,NUM_BLOCKS,
     &                                  IBLOCKL,global_communicator,-1,
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

C       finish timing part3
        timer3 = timer3 + interface_MPI_WTIME() - starttimer
        WALLTSTEP = SECTID(timer3)
        if(timing)then
          WRITE(LUWRT,9600) WALLTSTEP
        end if
      END IF
*     ^ NVEC+NROOT > MAXVEC or CONVER == .TRUE. or ITER == MAXIT

*     timing of this iteration

      WALLITR2 = interface_MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,'(/A,2A)')
     &         '                        Wall time for microiteration',
     &         '      ',WALLTID

      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
C
* ( End of loop over iterations )


      if (ITERSEOUT.eq.1) then

      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(LUWRT,1170) MAXIT
 1170    FORMAT(/'  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(LUWRT,1180) ITER
 1180    FORMAT(/,'  Convergence was obtained in ',I3,' iterations')
        END IF

      IF ( IPRT .GT. 0 ) THEN
        CALL REWINE(ILU1,LBLK)
        DO IROOT = 1, NROOT
          WRITE(LUWRT,'(/A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(LUWRT,*)
     &    '============================================'
          WRITE(LUWRT,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          IF (RTCNV(IROOT)) THEN
             WRITE(6,2190) IROOT,FINEIG(IROOT)+EIGSHF
          ELSE
             WRITE(6,2191) IROOT,FINEIG(IROOT)+EIGSHF
          END IF
 2190 FORMAT(' The final eigenvalue',I5,F22.10,' (converged)')
 2191 FORMAT(' The final eigenvalue',I5,F22.10,' (NOT converged)')
          IF(IPRT.GE.400) THEN
            WRITE(LUWRT,2200)
 2200       FORMAT(' The final approximation to eigenvector ')
            CALL WRTVCD(VEC1,ILU1,0,LBLK)
          END IF
          WRITE(LUWRT,2300)
 2300     FORMAT(/,'  Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(LUWRT,2310)
 2310     FORMAT
     &    (/,'  Iteration point        Eigenvalue         Residual ')
          DO I=1,ITER
 2330     WRITE(LUWRT,2340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 2340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
        END DO
        END DO
      ELSE
        CALL REWINE(ILU1,LBLK)
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
        DO IROOT = 1, NROOT
C          FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
          FINEIG(IROOT) = EIG(IROOT,ITER)+EIGSHF
          IF (RTCNV(IROOT)) THEN
             WRITE(6,2190) IROOT,FINEIG(IROOT)
          ELSE
             WRITE(6,2191) IROOT,FINEIG(IROOT)
          END IF
        END DO
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
      END IF


      IF(IPRT .EQ. 1 ) THEN
        DO IROOT = 1, NROOT
          WRITE(LUWRT,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
        END DO
      END IF

      else   ! new output

      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         WRITE(LUWRT,'(//A,I4,A)') ' (MICDV4_ENLMD_REL_PAR) '//
     &           'WARNING: maximum number of micro iterations,',
     &           MAXIT, ', is reached, CI aborted.'
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         WRITE(LUWRT,'(//A)')
     &        ' (MICDV4_ENLMD_REL_PAR) Micro iterations converged.'
      END IF

      do iroot = 1, nroot
        fineig(iroot) = eig(iroot,iter)+eigshf
      end do
      WRITE(LUWRT,'(/A,F20.8)')
     &   '                        Core energy        = ',eigshf
      WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &   '                        Active energies    = ',
     &                            EIG(1:NROOT,ITER)
      WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &   '                        Final CI energies  = ',
     &     fineig(1:nroot)

      do iroot = 1, nroot
        IF (RTCNV(IROOT)) THEN
          WRITE(LUWRT,'(/A,i4,a)')
     &    '                        root ',iroot,' ...... converged!'
        else
          WRITE(LUWRT,'(/A,i4,a)')
     &    '                        root ',iroot,' did not converge!'
        end if
      end do
C      IF ( IPRT .GT. 0 ) THEN
        WRITE(LUWRT,'(/A,1P,D8.2/A,D8.2)')
     &     ' Gradient convergence threshold ',thres_G,
     &     ' Energy   convergence threshold ',thres_E

        DO 1600 IROOT = 1, NROOT
          WRITE(LUWRT,'(/A,I3/A/)')
     &  ' Information about convergence for root... ' ,IROOT,
     &  ' ============================================'
          FINEIG(IROOT) = EIG(IROOT,ITER)
          WRITE(LUWRT,1300)
 1300     FORMAT(/,'  Summary of iterations',
     +           /,' -----------------------')
          WRITE(LUWRT,1310)
 1310     FORMAT
     &    (/'  Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(LUWRT,1340) I,EIG(IROOT,I)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(7X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
C      END IF

      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(LUWRT,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
      end if

C     store information about energies for each root on file
C     KRCI_CVECS.INFO (fh LU_INFO) - used in property calculations
      IF( MYPROC .eq. MASTER ) THEN
        DO I = 1, NROOT
           FINEIG(I) = EIG(I,ITER)+EIGSHF
        END DO
        WRITE(FILELAB,'(A5,A3)') "eroot",SYMFLABEL
        CALL KRCI_PRPFILE(LU_INFO,FILELAB,FINEIG,NROOT,0)

        ! for +Q correction
        open(file="energies.CI",unit=10,status="unknown",
     &       form="unformatted",access="sequential")
!       TODO: check proper reference energies for states other than the ground state and MCSCF reference wave functions!!!
        write(10) nroot,eig(1,1)+eigshf
        write(10) ((eig(i,iter)+eigshf),i=1,nroot)
        close(10,status="keep")
      END IF

      call qexit('MICDV4')

 1030 FORMAT(/3X,7F15.8,/,(3X,7F15.8))
 1120 FORMAT(/3X,I3,7F15.8,/,(6X,7F15.8))

 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION            : ',A)
 9250 FORMAT(' >>>  WALL TIME IN STEP 1 OF CURRENT ITERATION   : ',A)
 9400 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A/)
 9401 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9600 FORMAT(' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A)
 9350 FORMAT(' >>>  WALL TIME FOR PART 2.1 - 2.2               : ',A)
 9777 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9888 FORMAT(' >>>  WALL TIME FOR INITIAL ITERATION            : ',A)
*
      END
***********************************************************************

      SUBROUTINE GET_CNEWCOLD_BAS2(CN,CNP,CNO,NVEC,NVECP,NROOT,SCR,
     &           NVECUT,IPRT,THRES,RTCNV)
*
* A subspace expansion of eigenvectors is given by CN
* The previous set of eigenvectors is given in CNP
* Obtain orthogonal expansion of basis for New + previous NROOT eigenvector expansions
*
*. Jeppe Olsen, August 2019 - from GET_CNEWCOLD_BAS 
*
*
C      INCLUDE 'implicit.inc'
      Implicit None
      REAL*8 INPROD
*. Input
      Real(Kind=8) CN(NVEC,NROOT), CNP(NVECP,NROOT), Thres
      LOGICAL RTCNV(NROOT)
      Integer NVec, NVecP, NRoot, IPrt
*. Output
      Real(Kind=8) CNO(NVEC,*)
      Integer NVecUt
*. Scratch: Min size: MAX(4*2*NROOT*NVEC,4*(2*NROOT)**2)
      Real(Kind=8) SCR(3*NVEC**2)
*. Local scratch
      Integer, PARAMETER :: MXL_NROOT = 1000
      INTEGER ISCRL(MXL_NROOT)
      !Local variables
      Real(Kind=8) Zero, One, OneM, Thres2, XOrt, XNorm
      Real(Kind=8) Scale, XMax_Dia, XMax_Ofd
      !Functions
      Real(Kind=8) Sqrt
      !Scratch offsets, print flags
      Integer KFree, KS, KX, KScr, KVec, NTest
      Integer Irc, NCorr, IRoot, NDO, MaxIt, Iter
      

*
*. Corrections as second set of vectors
      IRC = 2
*. Partitioning of SCR
      KS = 1
      KFREE = KS + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KX = KFREE
      KFREE = KX + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KSCR = KFREE
      KFREE = KSCR + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KVEC = KFREE
      KFREE = KVEC + 2*NROOT
*
      ZERO = 0.0D0
      ONE = 1.0D0
      ONEM = -1.0D0
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GET_CNEWCOLD_BAS2 '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,'(A,3I4)') ' NROOT, NVEC, NVECP = ', NROOT, NVEC, NVECP
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' The initial CN matrix '
       CALL WRTMAT(CN,NVEC,NROOT,NVEC,NROOT)
       WRITE(6,*) ' The CNP matrix '
       CALL WRTMAT(CNP,NVECP,NROOT,NVECP,NROOT)
      END IF
*. Find the previous eigenvectors that has the largest overlap with new eigenvectors
COLD  IZERO = 0
COLD  CALL ISETVC(ISCRL,IZERO,NROOT)
COLD  DO IROOT = 1, NROOT
COLD   OVLAP_MAX = 0.0D0
COLD   DO IROOTP = 1, NROOT
COLD     IF(ISCRL(IROOTP).EQ.0) THEN
COLD       OVLAP = INPROD(CN(1,IROOT),CNP(1,IROOTP),NVECP)
COLD       IF(ABS(OVLAP).GT.
COLD     END IF
COLD    END DO
COLD  END DO
*
*. CNO: First Nroot vectors are the new vectors, 
*. the last NROOT are the original first roots(IRC = 1) or corrections
*. IRC = 2
      THRES2 = 1.0D-11
      CALL SETVEC(CNO,ZERO,2*NROOT*NVEC)
      CALL COPVEC(CN,CNO,NROOT*NVEC)
      NCORR = 0
*. Change the last NROOT vectors to corrections
      DO IROOT = 1, NROOT
*. Part of vector that is not in previous eigenvector-space
        XORT = INPROD(CN(NVECP+1,IROOT),CN(NVECP+1,IROOT),NVEC-NVECP)
        XORT = SQRT(XORT)
        IF(XORT.GT.THRES2.AND. .NOT.RTCNV(IROOT)) THEN
*. There is a significant part of root IROOT that is outside the first NVECP vectors, so include also a vector for previous vector.
          NCORR = NCORR + 1
*
*. We will assume that there is no root-flipping, so old and new root IROOT describe the same state. Could and 
*. should be improved. Also, it is assumed that the the converged root do not change from their 
*. original position - hold nearly always, but checks shoul be added.
* 
*. In CNO(*, IROOT) = New Iroot - Previous IROOT
          CALL COPVEC(CNO(1,IROOT),CNO(1,NROOT+NCORR),NVEC)
          CALL VECSUM(CNO(1,NROOT+NCORR),CN(1,IROOT),CNP(1,IROOT),
     &    ONE,ONEM,NVECP)
C?        IF(NTEST.GE.100) THEN
C?          WRITE(6,'(A)') ' New vector for old eigenvector before ort.'
C?          CALL WRTMAT(CNO(1,NROOT+NCORR),1,NVEC,1,NVEC)
C?        END IF
*
          XNORM = INPROD(CNO(1,NROOT+NCORR),CNO(1,NROOT+NCORR),
     &                   NVEC)
          SCALE = 1.0D0/SQRT(XNORM)
          CALL SCALVE(CNO(1,NROOT+NCORR),SCALE,NVEC)
        END IF
      END DO
        
      NDO = NROOT+NCORR
*
* Iterate over orthogonalizations
*
      MAXIT = 2
      DO ITER = 1, MAXIT
        IF(NTEST.GE.100) WRITE(6,*) ' Info from orth. iter. ', ITER
        IF(ITER.GT.1) THEN
          NDO = NVECUT
        END IF
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The CNO matrix '
          CALL WRTMAT(CNO,NVEC,NDO,NVEC,NDO)
        END IF
*. Overlap matrix S = CNO^T CNO
C         MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                    NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        CALL MATML7(SCR(KS),CNO,CNO,NDO,NDO,NVEC,NDO,NVEC,NDO,ZERO,ONE,
     &              1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Overlap of NDO basis '
          CALL WRTMAT(SCR(KS),NDO,NDO,NDO,NDO)
        END IF
*. Perform Gram-Schmidt orthogonalization
C       MGS4(X,S,NDIM,SCR1,THRES,NVECUT)
        CALL MGS4(SCR(KX),SCR(KS),NDO,SCR(KVEC),THRES,NVECUT)
C   TRNMAD(A,X,SCR,NDIMI,NDIMO)
*.   Check that delivered X orthonormalizes
        CALL TRNMAD(SCR(KS),SCR(KX),SCR(KSCR),NDO,NVECUT)
        IF(NTEST.GE.100) THEN
             WRITE(6,*) ' New overlap matrix '
           CALL WRTMAT(SCR(KS),NVECUT,NVECUT,NVECUT)
        END IF
        IF(IPRT.GE.10)
     &  CALL CHECK_UNIT_MAT(SCR(KS),NVECUT,XMAX_DIA,XMAX_OFD,0)
*. Transformation from original to orthonormal basis
        CALL MATML7(SCR(KS),CNO,SCR(KX),NVEC,NVECUT,
     &              NVEC,NDO,NDO,NVECUT,ZERO,ONE,0)
        CALL COPVEC(SCR(KS),CNO,NVEC*NVECUT)
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,'(A,2I5)') 
     &  ' Output from GET_CNEW.., NROOT and NVECUT = ',
     &  NROOT, NVECUT
      END IF
*
      IF(NTEST.GE.1.AND.NVECUT.LT.2*NROOT) THEN
        WRITE(6,'(A,2I5)') 
     &  ' GET_CNEW.., Reduced number of vectors: 2*NROOT and NVECUT = ',
     &  2*NROOT, NVECUT
      END IF
        
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_CNEWCOLD_BAS2_CPX(CN,CNP,CNO,NVEC,NVECP,NROOT,SCR,
     &           NVECUT,IPRT,THRES,RTCNV)
*
* A subspace expansion of eigenvectors is given by CN
* The previous set of eigenvectors is given in CNP
* Obtain orthogonal expansion of basis for New + previous NROOT eigenvector expansions
*
*. Jeppe Olsen, August 2019 - from GET_CNEWCOLD_BAS 
*  Andreas Nyvang, September 2019 - complex version of GET_CNEWCOLD_BAS2
*
*
C      INCLUDE 'implicit.inc'
      Implicit none
      REAL*8 INPROD
*. Input
      Integer NVec, NVecP, NRoot, IPrt
      Real(Kind=8) CN(2*NVEC,NROOT), CNP(2*NVECP,NROOT)
      Real(Kind=8) Thres
      LOGICAL RTCNV(NROOT)
*. Output
      Integer NVecUt
      Real(Kind=8) CNO(2*NVEC,*)
*. Scratch: Min size: MAX(4*2*NROOT*NVEC,4*(2*NROOT)**2)
      Real(Kind=8) SCR(6*NVEC**2)
*. Local scratch
      Integer, PARAMETER :: MXL_NROOT = 1000
      INTEGER ISCRL(MXL_NROOT) 

      !Local variables
      Real(Kind=8) Zero, One, OneM, Thres2, XOrt, XNorm
      Real(Kind=8) Scale, XMax_Dia, XMax_Ofd
      !Functions
      Real(Kind=8) Sqrt, DDot
      !Scratch offsets, print flags
      Integer KFree, KS, KX, KScr, KVec, NTest
      Integer Irc, NCorr, IRoot, NDO, MaxIt, Iter
      Integer LDCNO, LDCNO1, LDCNO2, LDKS, LDKX

*
*. Corrections as second set of vectors
      IRC = 2
*. Partitioning of SCR
      KS = 1
      KFREE = KS + 2*(2*NROOT)*MAX(2*NROOT,NVEC)
*
      KX = KFREE
      KFREE = KX + 2*(2*NROOT)*MAX(2*NROOT,NVEC)
*
      KSCR = KFREE
      KFREE = KSCR + 2*(2*NROOT)*MAX(2*NROOT,NVEC)
*
      KVEC = KFREE
      KFREE = KVEC + 4*NROOT
*
      ZERO = 0.0D0
      ONE = 1.0D0
      ONEM = -1.0D0
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GET_CNEWCOLD_BAS2_CPX '
        WRITE(6,*) ' =============================== '
        WRITE(6,*)
        WRITE(6,'(A,3I4)') ' NROOT, NVEC, NVECP = ', NROOT, NVEC, NVECP
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' The initial CN matrix '
       CALL WRTMAT(CN,2*NVEC,NROOT,2*NVEC,NROOT)
       WRITE(6,*) ' The CNP matrix '
       CALL WRTMAT(CNP,2*NVECP,NROOT,2*NVECP,NROOT)
      END IF
*. Find the previous eigenvectors that has the largest overlap with new eigenvectors
COLD  IZERO = 0
COLD  CALL ISETVC(ISCRL,IZERO,NROOT)
COLD  DO IROOT = 1, NROOT
COLD   OVLAP_MAX = 0.0D0
COLD   DO IROOTP = 1, NROOT
COLD     IF(ISCRL(IROOTP).EQ.0) THEN
COLD       OVLAP = INPROD(CN(1,IROOT),CNP(1,IROOTP),NVECP)
COLD       IF(ABS(OVLAP).GT.
COLD     END IF
COLD    END DO
COLD  END DO
*
*. CNO: First Nroot vectors are the new vectors, 
*. the last NROOT are the original first roots(IRC = 1) or corrections
*. IRC = 2

      THRES2 = 1.0D-11
      CALL SETVEC(CNO,ZERO,4*NROOT*NVEC)
      CALL COPVEC(CN,CNO,2*NROOT*NVEC)
      NCORR = 0

*. Change the last NROOT vectors to corrections
      DO IROOT = 1, NROOT
*. Part of vector that is not in previous eigenvector-space
C Real part of overlap     
        XORT = DDot(NVec-NVecP,CN(NVECP+1,IROOT),1,CN(NVECP+1,IROOT),1)

C Imaginary part
        XORT = XORT + 
     &         DDot(NVec-NVecP,CN(NVECP+NVEC+1,IROOT),1,
     &         CN(NVECP+NVEC+1,IROOT),1)

        XORT = SQRT(XORT)

        IF(XORT.GT.THRES2.AND. .NOT.RTCNV(IROOT)) THEN
*. There is a significant part of root IROOT that is outside the first NVECP vectors, so include also a vector for previous vector.
          NCORR = NCORR + 1

*. We will assume that there is no root-flipping, so old and new root IROOT describe the same state. Could and 
*. should be improved. Also, it is assumed that the the converged root do not change from their 
*. original position - hold nearly always, but checks shoul be added.

*. In CNO(*, IROOT) = New Iroot - Previous IROOT
          CALL DCopy(2*NVec,CNO(1,IROOT),1,CNO(1,NROOT+NCORR),1)
C          CALL COPVEC(CNO(1,IROOT),CNO(1,NROOT+NCORR),2*NVEC)
C Real part
          CALL VECSUM(CNO(1,NROOT+NCORR),CN(1,IROOT),CNP(1,IROOT),
     &    ONE,ONEM,NVECP)
C Imaginary part
          CALL VECSUM(CNO(1+nvec,NROOT+NCORR),CN(1+nvec,IROOT),
     &    CNP(1+nvecp,IROOT),
     &    ONE,ONEM,NVECP)

        IF(NTEST.GE.100) THEN
          WRITE(6,'(A)') ' New vector for old eigenvector before ort.'
          CALL WRTMAT(CNO(1,NROOT+NCORR),1,NVEC,1,NVEC)
        END IF
*
C Real part
          XNORM = INPROD(CNO(1,NROOT+NCORR),CNO(1,NROOT+NCORR),
     &                   NVEC)
C Imaginary part
          XNORM = XNORM + INPROD(CNO(1+nvec,NROOT+NCORR),
     &    CNO(1+nvec,NROOT+NCORR),
     &                   NVEC)
          SCALE = 1.0D0/SQRT(XNORM)
          CALL SCALVE(CNO(1,NROOT+NCORR),SCALE,2*NVEC)
        END IF
      END DO
        
      NDO = NROOT+NCORR
*
* Iterate over orthogonalizations
*
      MAXIT = 2
      DO ITER = 1, MAXIT
        IF(NTEST.GE.100) WRITE(6,*) ' Info from orth. iter. ', ITER

        IF(ITER.GT.1) THEN
          NDO = NVECUT
        END IF

        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The CNO matrix (real first nvec rows) '
          CALL WRTMAT(CNO,2*NVEC,NDO,2*NVEC,NDO)
        END IF

*. Overlap matrix S = CNO^H CNO
        LDCNO1 = max(1,2*NVec)
        LDCNO2 = max(1,2*NVec)
        LDKS = max(1,2*NDO)

C Real part
        CALL DGEMM('T','N',NDO,NDO,2*NVEC,ONE,CNO,LDCNO1,
     &         CNO,LDCNO2,ZERO,SCR(KS),LDKS)

C Imaginary part
        CALL DGEMM('T','N',NDO,NDO,NVEC,ONE,CNO,LDCNO1,
     &         CNO(nvec+1,1),LDCNO2,ZERO,SCR(KS+NDO),LDKS)
        CALL DGEMM('T','N',NDO,NDO,NVEC,-1.0d0,CNO(nvec+1,1),LDCNO1,
     &         CNO,LDCNO2,ONE,SCR(KS+NDO),LDKS)

        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Overlap of NDO basis '
          CALL WRTMAT(SCR(KS),2*NDO,NDO,2*NDO,NDO)
        END IF

*. Perform Gram-Schmidt orthogonalization
        CALL MGS4_CPX(SCR(KX),SCR(KS),NDO,SCR(KVEC),THRES,NVECUT)

*.   Check that delivered X orthonormalizes
        Call DZero(Scr(KScr),2*NDO*NVECUT)

        CALL TRNMAD_CPX(SCR(KS),SCR(KX),SCR(KSCR),NDO,NVECUT)

        IF(NTEST.GE.100) THEN
             WRITE(6,*) ' New overlap matrix '
           CALL WRTMAT(SCR(KS),2*NVECUT,NVECUT,2*NVECUT,NVECUT)
        END IF

        IF(IPRT.GE.10)
     &  CALL CHECK_UNIT_MAT(SCR(KS),NVECUT,XMAX_DIA,XMAX_OFD,0)

*. Transformation from original to orthonormal basis

        LDCNO = max(1,2*NVec)
        LDKX = max(1,2*NDo)
        LDKS = max(1,2*NVec)

        !Real part
        Call DGemm('N','N',NVec,NVecUt,NDo,1.0D0,CNO,LDCNO,
     &             SCR(KX),LDKX,0.0D0,Scr(KS),LDKS)   

        Call DGemm('N','N',NVec,NVecUt,NDo,-1.0D0,CNO(NVec+1,1),LDCNO,
     &             SCR(KX+NDO),LDKX,1.0D0,Scr(KS),LDKS)   

        !Imaginary part
        Call DGemm('N','N',NVec,NVecUt,NDo,1.0D0,CNO(NVec+1,1),LDCNO,
     &             SCR(KX),LDKX,0.0D0,Scr(KS+NVec),LDKS)   
        
        Call DGemm('N','N',NVec,NVecUt,NDo,1.0D0,CNO,LDCNO,
     &             SCR(KX+NDO),LDKX,1.0D0,Scr(KS+NVec),LDKS)   

        CALL COPVEC(SCR(KS),CNO,2*NVEC*NVECUT)
         
C        write(6,*) 'New CNO matrix'
C        call wrtmat(CNO,2*NVec,NVecut,2*NVec,NVecut)
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,'(A,2I5)') 
     &  ' Output from GET_CNEW.., NROOT and NVECUT = ',
     &  NROOT, NVECUT
      END IF
*
      IF(NTEST.GE.1.AND.NVECUT.LT.2*NROOT) THEN
        WRITE(6,'(A,2I5)') 
     &  ' GET_CNEW.., Reduced number of vectors: 2*NROOT and NVECUT = ',
     &  2*NROOT, NVECUT
      END IF
        
*
      END 

      SUBROUTINE MGS4_CPX(X,S,NDIM,SCR1,THRES,NVECUT)

* Modified Gram-Schmidt procedure by forward orthogonalization

*  watch out for zero columns indicating linear dependency

* Jeppe Olsen, March 2013, added thres to MGS3

* Andreas Nyvang, September 2019: Complex version of MGS4

* S is input overlap matrix, X is output set of orthonormalized vectors
* First NDim entries in a column of S and X are the real parts,
* the latter NDim entries are the imaginary parts

* Thres is min norm of linear independent vector- only meaningfull if all
* initial vectors have identical norm


      Implicit none
      REAL*8 INPROD
*. input
      Integer NDim
      Real(kind=8) S(2*NDIM,NDIM), Thres
*. Output
      Real(kind=8) X(2*NDIM,*)
      Integer NVecUt
*. Scratch : vector of length NDIM
      Real(kind=8) SCR1(*)
      
      Real(kind=8) Zero, One, Ddot, XNorm, Factor
      Real(kind=8) XSX_Re, XSX_Im
      Integer NTest, I, IVec, JVec

      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from MGS4 '
        WRITE(6,*) ' ==============='
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Initial overlap matrix (real first ndim rows)'
       CALL WRTMAT(S,2*NDIM,NDIM,2*NDIM,NDIM)
      END IF

*. Initialize X to unit matrix

      ZERO = 0.0D0
      ONE = 1.0D0
      CALL SETVEC(X,ZERO,2*NDIM**2)

      Do I =1,NDim
         X(I,I) = 1.0d0
      End do
*
      DO IVEC = 1, NDIM
*. Normalize vector IVEC - using Scr1 for, well, scratch...

        Call DCopy(2*NDim,S(1,IVec),1,Scr1,1)
        !Jeppe's procedure originally calculated the inner product of S*X(1,IVec),
        !but there is no real numerical difference from just calculating
        !the inner product between X(1,IVec) and S(1,IVec)
        !because the initialisation of X should already make the other
        !elements zero
C        Call DGemv('N',NDim,NDim,1.0D0,S,max(1,2*NDim),
C     &             X(1,IVec),1,0.0D0,Scr1,1)
C
C        Call DGemv('N',NDim,NDim,-1.0D0,S(NDim+1,1),max(1,2*NDim),
C     &             X(NDim+1,IVec),1,1.0D0,Scr1,1)
C
C        Call DGemv('N',NDim,NDim,1.0D0,S,max(1,2*NDim),
C     &             X(NDim+1,IVec),1,0.0D0,Scr1(NDim+1),1)
C
C        Call DGemv('N',NDim,NDim,1.0D0,S(NDim+1,1),max(1,2*NDim),
C     &             X(1,IVec),1,1.0D0,Scr1(NDim+1),1)

C        write(6,*) 'Scr1 after inprod'
C        Call WrtMat(Scr1,2*NDim,1,2*NDim,1)

        !Procedure is made such that only real elements are left in the
        !inner product between X(1,IVec) and S(1,IVec), so
        XNORM = InProd(X(1,IVec),Scr1,2*NDim)

        IF(NTEST.GE.100) THEN
          WRITE(6,'(A,I4,E15.7)') ' IVEC, XNORM = ', IVEC, XNORM
        END IF
*
        IF (XNorm.LE.THRES) THEN
          Factor = 0.0D0
        ELSE
          FACTOR = 1.0D0/SQRT(XNORM)
        END IF

C        write(6,*) 'Factor ',Factor

        CALL SCALVE(X(1,IVEC), FACTOR, 2*NDIM)
        CALL SCALVE(SCR1,FACTOR,2*NDIM)

C        write(6,*) 'IVec ',IVec
C        write(6,*) 'X after scaling'

C        Call Wrtmat(X(1,IVec),2*NDim,1,2*NDim,1)

C        write(6,*) 'Scr1 after scaling'

C        Call Wrtmat(Scr1,2*NDim,1,2*NDim,1)

*. Subtract X(1,IVEC) from all remaining vectors
        DO JVEC = IVEC+1,NDIM

          XSX_Re = INPROD(Scr1,X(1,JVEC),2*NDIM)

          XSX_Im = INPROD(Scr1,X(NDim+1,JVEC),NDIM)

C          XSX_Re = XSX_Re + INPROD(Scr1(NDim+1),X(NDim+1,JVEC),NDIM)

          XSX_Im = XSX_Im - INPROD(Scr1(NDim+1),X(1,JVEC),NDIM)

C          write(6,*) 'XSX_Re ',XSX_Re
C          write(6,*) 'XSX_Im ',XSX_Im

C          CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,IVEC),ONE,-XSX,NDIM)

C         write(6,*) 'Old X for JVec = ',JVec
C        Call wrtmat(X(1,JVec),2*NDim,1,2*NDim,1)

          !Subtraction, real parts
          CALL VECSUM(X(1,JVEC),X(1,JVEC),
     &                X(1,IVEC),ONE,-XSX_Re,NDIM)

          CALL VECSUM(X(1,JVEC),X(1,JVEC),
     &                X(1+NDim,IVEC),ONE,XSX_Im,NDIM)

          !Imaginary parts
          CALL VECSUM(X(1+NDim,JVEC),X(1+NDim,JVEC),
     &                X(1+NDim,IVEC),ONE,-XSX_Re,NDIM)
          
          CALL VECSUM(X(1+NDim,JVEC),X(1+NDim,JVEC),
     &                X(1,IVEC),ONE,-XSX_Im,NDIM)

C          write(6,*) 'New X for JVec = ',JVec
C          Call wrtmat(X(1,JVec),2*NDim,1,2*NDim,1)
        END DO
      END DO
*
*. And remove zero vectors
*
      NVECUT = 0
      DO IVEC = 1, NDIM
        XNORM = INPROD(X(1,IVEC),X(1,IVEC),2*NDIM)
        IF(XNORM.GT.0.0D0) THEN
          NVECUT = NVECUT + 1
          IF(NVECUT.NE.IVEC) CALL COPVEC(X(1,IVEC),X(1,NVECUT),2*NDIM)
        END IF
      END DO
*
      IF(NTEST.GE.1.AND.NVECUT.NE.NDIM) THEN
        WRITE(6,*)' MGS4 reduced dim, from and to ', NDIM,NVECUT
      ELSE IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization information:'
        WRITE(6,*) ' Number of linear independent vectors ', NVECUT
        WRITE(6,*) ' Orthonormalized vectors '
        CALL WRTMAT(X,2*NDIM,NVECUT,2*NDIM,NVECUT)
      END IF
*

      END
***********************************************************************
      SUBROUTINE CMICDV_PAR(VEC1,VEC2,RNRM,EIG,FINEIG,
     &                      MAXIT,NVAR,NROOT,MAXVEC,NINVEC,T_CC,T_BUFF,
     &                      APROJI,APROJR,AVEC,AVECO,WORK,IPRT,
     &                      NPRDIM,IPNTR,LBLK,EIGSHF,
     &                      IHAM12,thres_G,thres_E,
     &                      IBLOCKL,NPARBLOCK,NBLOCKS,
     &                      LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                      LUCLIST,NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      RCCTOS,SCRRED_R,SCRRED_I,IGROUPLIST,
     &                      IPROCLIST,IT_TTPL,IT_TTOL,AMATC,CWORK,
     &                      RWORK,EVL,eciold,IWORK)
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
* Revised 2020, Andreas Nyvang, use by default 2 vectors instead of 1 for Davidson space
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
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#include "parluci.h"
      REAL * 8   INPROD,INPRDD
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12, WALLTSTEP*12
      CHARACTER WPART22*12, FILELAB*8
      LOGICAL CONVER,RTCNV(100),EXSTOP
*
#include "krmcluci_inf.h"
*
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
C      DIMENSION LU6LIST(*), LU7LIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
      DIMENSION SCRRED_R(*), SCRRED_I(*)
      INTEGER   RCCTOS(*)
      DIMENSION IGROUPLIST(*), IPROCLIST(*)
      DIMENSION VEC1(*),VEC2(*),eciold(*)
      DIMENSION RNRM(NROOT,MAXIT),EIG(2*MAXVEC,MAXIT)
      DIMENSION APROJI(*),APROJR(*),AVEC(*),AVECO(*),WORK(*)
      DIMENSION IPNTR(1), IBLOCKL(*), NPARBLOCK(*)
      DIMENSION IT_TTPL(*)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IT_TTOL(*)
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
      real*8, allocatable, dimension(:)  :: scrcno, AVec_Test
      real*8, allocatable, dimension(:)  :: CNO, Eig_test
*
      call qenter('CMICDV_PAR')
*
      IDONEWBAS = 1
      NTEST = 0000
      IPRT = MAX(NTEST,IPRT)
C
!     !> (re-)set convergence threshold
       if(abs(thres_G) < 1.0d-07) then
          write(luwrt,'(/A,1P,E10.2,A)')
     &    'INFO: CI gradient convergence threshold too tight,'//
     &    ' changed from',thres_G,' to 1.0e-7'
          thres_G = 1.0d-07
       end if
C      write(6,*) 'thres_e ',thres_e

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
      write(LUWRT,*)
      write(LUWRT,*) ' ----------------------------------------------'
      write(LUWRT,*) " I'll take you to a place"
      write(LUWRT,*) ' Where we shall find our'
      write(LUWRT,*) '  ROOTS, BLOODY ROOTS'
      write(LUWRT,*)
      write(LUWRT,*) '               Sepultura'
      write(LUWRT,*) ' ----------------------------------------------'
      write(LUWRT,*)

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
      xinit_time = interface_MPI_WTIME()
*              
      CALL DZERO(APROJR, MAXVEC*(MAXVEC+1)/2 )
      CALL DZERO(APROJI, MAXVEC*(MAXVEC+1)/2 )
      CALL DZERO(SCRRED_R, MAXVEC*(MAXVEC+1)/2 )
      CALL DZERO(SCRRED_I, MAXVEC*(MAXVEC+1)/2 )
*
      DO 10 IVEC = 1,NINVEC
*
*       copy c-vector to working-file ILUC
        call interface_mpi_BARRIER(MYNEW_COMM)
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
        sigmatime = interface_MPI_WTIME()
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
        sigmatime2 = interface_MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
!       WRITE(LUWRT,9777) WALLTID
!       WRITE(LUWRT,*) '   '
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
      xinit_time2 = interface_MPI_WTIME()
      WALLTID = SECTID(xinit_time2-xinit_time)
!     WRITE(LUWRT,9888) WALLTID
*
*     synchronize global_communicator
*
      CONE = - 1.0D0
      IREDL = 0
      IREDL = NINVEC*(NINVEC-1)/2 + NINVEC
      CALL DZERO(SCRRED_R,IREDL)
      CALL DZERO(SCRRED_I,IREDL)
*     ... sum up
      CAll redvec(APROJR,SCRRED_R,IREDL,2,op_MPI_SUM,
     &                global_communicator,-1)
      CAll redvec(APROJI,SCRRED_I,IREDL,2,op_MPI_SUM,
     &                global_communicator,-1)
*     ... update APROJR and APROJI 
      CALL DCOPY(IREDL,SCRRED_R,1,APROJR,1)
      CALL DSCAL(IREDL,-1.0D0,SCRRED_I,1)
C      CALL DSCAL(IREDL,CONE,SCRRED_I,1)
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
      IPRT = 00
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
      If(MyProc == Master) Then
      CALL CEIGEN_NEW(APROJR,APROJI,NINVEC,EIG(1,1),AVEC,AMATC,CWORK,
     &                RWORK,EVL,IWORK)
      End if

      call interface_mpi_bcast(EIG(1,1),2*NINVEC,MASTER,
     &                             global_communicator)
      call interface_mpi_bcast(AVEC(1),2*NINVEC**2,MASTER,
     &                             global_communicator)


C
      IF(IPRT .GE. 3) THEN
        WRITE(LUWRT,'(A,I4)') ' Eigenvalues of initial iteration '
        WRITE(LUWRT,'(F22.13)')
     &  ( EIG(IROOT,1)+EIGSHF,IROOT=1,NROOT)
      END IF

      IPRT = 0
      NVEC = NINVEC
C      call interface_mpi_bcast(NVec,1,Master,global_communicator)
*
      IF ( MAXIT .eq. 1) GOTO  901
*
*     ======================
*     Loop over iterations
*     ======================
*
 1000 CONTINUE
*       start timing of iteration
        WALLITR1 = interface_MPI_WTIME()
*
        WRITE(LUWRT,'(//A,I4)')
     &        ' (CMICDV_PAR)           CI microiteration no.',ITER
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

C         write(6,*) 'thres_e ',thres_e
         CALL P1_B_PAR_CX_1(VEC1,VEC2,AVEC,WORK,EIGAPR,RNRM,EIGSHF,
     &                      EIG,thres_g,thres_e,
     &                      RTCNV,CONVER,ITER,NROOT,MAXIT,
     &                      IROOT,LU1LIST,LU2LIST,LU5LIST,NBATV,LBATV,
     &                      LEBATV,I1BATV,IBATV,MY_LU1_OFF,MY_LU2_OFF,
     &                      MY_LU5_OFF,SCRRED_R,SCRRED_I,NVEC,2*NVEC,
     &                      ILU1,ILU2,ILU5,MaxVec)
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
C           XNORMO = XNORM
CSK        WRITE(LUWRT,*) 'THIS IS XNORMO',XNORMO
           IOlsTm = 0
           if(IOlsTm .ne. 0) then
C Now calculating Olsen correction
         
            CALL IZERO(LU4LIST,IALL_LU4)
            CALL IZERO(LU5LIST,IALL_LU5)
*

      Call Olsen_Corr_Cpx(Vec1,Vec2,AVec,Work,-EigAPr,IRoot,
     &                             LU1LIST,LU4LIST,LU5LIST,
     &                             LU3LIST,
     &                             NBatv,LBatv,LeBatv,I1Batv,IBatv,
     &                             MY_LU1_OFF,MY_LU4_OFF,
     &                             MY_LU5_OFF,MY_LU3_OFF, 
     &                             MY_DIA_OFF,NVEC,ILU1,
     &                             ILU4,ILU5,ILU3,IDIA,1)

         end if !Use Olsen correction
*
*          check for linear dependency...
*
         !Orthonormalise new direction wrt prior directions
C         write(6,*) myproc,' before orthonormalisation'
      Call P1_B_PAR_CX_3(VEC1,VEC2,Work,
     &                         LU1LIST,LU3LIST,LU1LIST,
     &                         NBatv,LBatv,LEBatv,I1Batv,IBatv,
     &                         MY_LU1_OFF,MY_LU3_OFF,
     &                         MY_LU1_OFF,SCRRED_R,SCRRED_I,
     &                         NVEC,IADD,ILU1,ILU3,ILU1)

C           CALL INPROD_B_PAR_CPLX2(ILU1,ILU3,VEC1,VEC2,SCRRED_R,
C     &                             SCRRED_I,
C     &                             NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                             MY_LU1_OFF,MY_LU3_OFF,LU1LIST,
C     &                             LU3LIST,NVEC,IADD)
C*
C*          communication is always essential ...
C           CALL DZERO(VEC1,NVEC+IADD-1)
C           CALL DZERO(VEC2,NVEC+IADD-1)
C*
C           CAll redvec(SCRRED_R,VEC1,NVEC+IADD-1,2,
C     &                     op_MPI_SUM,global_communicator,-1)
C           CAll redvec(SCRRED_I,VEC2,NVEC+IADD-1,2,
C     &                     op_MPI_SUM,global_communicator,-1)
C*
C*          transfer information to WORK array...
C           CALL DCOPY(NVEC+IADD-1,VEC1,1,WORK(1),2)
C           CALL DCOPY(NVEC+IADD-1,VEC2,1,WORK(2),2)
C*
C           XR = 0.0D0
C           XI = 0.0D0
C*
C           DO JVEC = 1,NVEC+IADD-1
C             XR = WORK(2*(JVEC-1)+1)
C             XI = WORK(2*(JVEC-1)+2)
C             XNORM = XNORM - XR**2 - XI**2
C           END DO
C*
C           THRSLDP = 1.0D-10
C           IF( XNORM / XNORMO .le. THRSLDP )THEN
C             WRITE(LUWRT,*)
C     &       ' CMICDV : Vector eliminated due to linear dependence '
C             IADD = IADD -1
C             GOTO 100
C           END IF
C*
C           SCALE = 1.0D0/SQRT(XNORM)
CCSK        WRITE(LUWRT,*) ' Scale, scale**2*xnorm',scale,
CCSK  &                      scale**2*xnorm
C*
C           CALL DSCAL( 2*(NVEC+IADD-1), SCALE, WORK, 1 )
CC           
CC          1.3 and 1.4 orthogonalize to all previous vectors 
CC          and normalize vector
C           ONEM = -1.0D0
C           CALL DSCAL(NVEC+IADD-1,ONEM,WORK(2),2)
CC
C           CALL P1_B_PAR_CX_2(VEC1,VEC2,WORK,SCALE,LU1LIST,LU3LIST,
C     &                        LU1LIST,
C     &                        NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                        MY_LU1_OFF,MY_LU3_OFF,MY_LU1_OFF,
C     &                        NVEC,IADD,ILU1,ILU3,ILU1)
C
         END IF
C        ^ not converged...
  100 CONTINUE
      IF (IPRT .GE. 10) THEN
          WRITE(LUWRT,'(/A,F20.8)')
     &      '                        Core energy        =',EIGSHF
          WRITE(LUWRT,'(A,4F20.8,/,(37X,4F20.8))')
     &      '                        Active energy      =',
     &                        EIG(1:NROOT,ITER)
       END IF
       WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &      '                        Total energy       =',
     &                        (EIG(I,ITER)+EIGSHF,i=1,NROOT)
       WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &      '                        Lowering of energy =',
     &      (eciold(i)-EIG(I,ITER)-EIGSHF,I=1,nroot)

       do i = 1, nroot
         eciold(i) = EIG(i,ITER) + EIGSHF
       end do

       WRITE(LUWRT,'(/A,1P,D8.2,A/,(I15,D20.5))')
     &   '                        Norm of CI residuals'//
     &   ' (thr = ',thres_G,')',
     &   (i,RNRM(i,ITER),i=1,NROOT)
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
      starttime = interface_MPI_WTIME()
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
        call interface_mpi_BARRIER(MYNEW_COMM)
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
        sigmatime = interface_MPI_WTIME()
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
        sigmatime2 = interface_MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
!       WRITE(LUWRT,9400) WALLTID
*
*       projected matrix using batch structure of CI vector(s)
        CALL INPROD_B_PAR_CPLX3(ILU1,ILU2,VEC1,VEC2,SCRRED_R,SCRRED_I,
     &                          NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                          MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                          LU2LIST,NVEC,IVEC,IMUSTRED,ISTRED)
*
 150  CONTINUE 
*
*      synchronize global_communicator
*
       CALL DZERO(APROJR(ISTRED),IMUSTRED)
       CALL DZERO(APROJI(ISTRED),IMUSTRED)
*
       CAll redvec(SCRRED_R(ISTRED),APROJR(ISTRED),IMUSTRED,2,
     &                 op_MPI_SUM,global_communicator,-1)
       CAll redvec(SCRRED_I(ISTRED),APROJI(ISTRED),IMUSTRED,2,
     &                 op_MPI_SUM,global_communicator,-1)
*
      endtime = interface_MPI_WTIME()
      tottime = 0.D0
      tottime = endtime - starttime
*     TIMING FOR PARTS 2.1 - 2.2
      WPART22 = SECTID(tottime)
!     WRITE(LUWRT,9350) WPART22
*
      IPRT = 00
*
      NVecO = NVec
      call dcopy(2*nveco*nveco,avec,1,aveco,1)
      !Sending previous eigenvector approximation to all processes
      Call interface_mpi_bcast(NVecO,1,Master,global_communicator)
      call interface_mpi_bcast(AVecO(1:2*NVECO**2),2*NVECO**2,MASTER,
     &                             global_communicator)

*     diagonalize projected matrix
C      write(6,*) 'Aveco new:'
C      Call Wrtmat(Aveco,2*Nveco,Nveco,2*NVeco,Nveco)
C      Write(LUWRT,*) ' NVEC changed', NVEC,' -> ', NVEC+IADD
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
      If (MyProc == Master) then
      CALL CEIGEN_NEW(APROJR,APROJI,NVEC,EIG(1,ITER),AVEC,AMATC,CWORK,
     &                RWORK,EVL,IWORK)
      End if

      call interface_mpi_bcast(EIG(1,ITER),2*NVEC,MASTER,
     &                             global_communicator)
       call interface_mpi_bcast(AVEC(1),2*NVEC**2,MASTER,
     &                             global_communicator)

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



      if((NVEC+NROOT .gt. MAXVEC).and.(.not.Conver)) then

         If(Myproc==Master) then

         nnvec1=max(2*nroot,nvec)
         allocate(cno(4*nnvec1**2))

         allocate(scrcno(6*nvec*nvec))

         THR_ORTBAS = 1.0d-7

C          write(6,*) 'nvec',nvec,'nroot',nroot,'nnvec1',nnvec1

        CALL GET_CNEWCOLD_BAS2_Cpx(AVEC,AVECO,cno,
     &                      NVEC,NVECO,NROOT,
     &                       scrcno(1),NVECUD,
     &                      IPRT,THR_ORTBAS,RTCNV)
C          write(6,*) 'nvec',nvec
C          write(6,*) 'NVECUD',NVECUD,'NADD',NVECUD-NROOT

C         NADD=NVECUD-NROOT
!        IF (IPRT.GE.3 ) THEN
C           WRITE(LUWRT,*) 'NADD SET to:',NROOT,' -> ', NADD
!         ENDIF
         CALL COPVEC(CNO,AVEC,2*NVECUD*NVEC)
C         write(6,*) 'New Avec'
C         Call wrtmat(Avec,2*NvecUd,NVec,2*NVecUd,NVec)

         deallocate(cno,scrcno)
         !Sending reduced space of eigenvector approximation to all
         !processes

         End If
         Call interface_mpi_bcast(NVecUd,1,Master,global_communicator)
         call interface_mpi_bcast(AVEC(1:2*NVecUd*NVec),2*NVECUD*NVEC,
     &                            MASTER,global_communicator)
      Else
         NVecUd = NRoot
         Call interface_mpi_bcast(NVecUd,1,Master,global_communicator)

      End If


          

C      write(6,*) 'myproc ',myproc,'nvec ',nvec
C      write(6,*) 'myproc ',myproc,'nvecud ',nvecud

C      Call Wrtmat(Avec,2*NVecUd,NVec,2*NVecUd,NVec)

        Call DZero(Work,NRoot)

        Call SetVec(Work,1.0D0,NRoot)

      Call IZero(LU3LIST,IAll_LU3)
      
      Call TRAVC_B_CX_DRV(VEC1,VEC2,AVec,LU1LIST,
     &                          LU3LIST,NBatv,LBatv,LeBatv,I1Batv,
     &                          IBatv,MY_LU1_OFF,MY_LU3_OFF,
C     &                          NVEC,NVECUD,ILU1,ILU3,IALL_LU1)
     &                     NVEC,NVECUD,ILU1,ILU3,NVecUd*My_Act_Blk1)

      Call IZero(LU3LIST,IAll_LU3)

      Call TRAVC_B_CX_DRV(VEC1,VEC2,AVec,LU2LIST,
     &                          LU3LIST,NBatv,LBatv,LeBatv,I1Batv,
     &                          IBatv,MY_LU2_OFF,MY_LU3_OFF,
C     &                          NVEC,NVECUD,ILU2,ILU3,IALL_LU2)
     &                     NVEC,NVECUD,ILU2,ILU3,NVecUd*My_Act_Blk2)

            CALL IZERO(LU3LIST,IALL_LU3)
      Call Interface_mpi_barrier(global_communicator)


C*
C*       c vectors from ILU1 to ILU1
C        DO IVec = 1, NVecUd
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IVec-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = -AVEC((IVec-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU1LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU1_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IVec,ILU1,ILU3)
C        END DO
C*
C*       no scaling, we should already work in a normalized basis
C        CALL IZERO(LU1LIST,NVecUd*MY_ACT_BLK2)
C*
C        DO IVec = 1, NVecUd
C          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU1LIST,
C     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                         MY_LU3_OFF,MY_LU1_OFF,IVec,ILU3,ILU1)
CC          WORK(2*NROOT*NVEC+IROOT) = 1.0D0
C        END DO
C*
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C*
C*       corresponding sigma vectors from ILU2 to ILU2
C        DO IVec = 1, NVecUd
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IVec-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = - AVEC((IVec-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU2LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU2_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IVec,ILU2,ILU3)
C        END DO
C*
C*       no scaling, we should already work in a normalized basis
C        CALL IZERO(LU2LIST,NVecUd*MY_ACT_BLK2)
C*
C        DO IVec = 1, NVecUd
C          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU2LIST,
C     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                         MY_LU3_OFF,MY_LU2_OFF,IVec,ILU3,ILU2)
C        END DO
C*
C        CALL IZERO(LU3LIST,IALL_LU3)

      CALL IZERO(LU3LIST,IALL_LU3)
C      CALL IZERO(LU4LIST,IALL_LU4)

C      write(LUWRT,*) 'NVEC ', NVEC,' reset to ',nvecud
      NVEC = NVECUD
      call interface_mpi_bcast(nvec,1,Master,global_communicator)


        CALL DZERO(AVEC,2*NVEC**2)
        CALL DZERO(APROJR,NVEC*(NVEC+1)/2)
        CALL DZERO(APROJI,NVEC*(NVEC+1)/2)
C        CALL DZERO(AVEC,2*MAXVEC**2)
C        CALL DZERO(APROJR,MAXVEC*(MAXVEC+1)/2)
C        CALL DZERO(APROJI,MAXVEC*(MAXVEC+1)/2)
*
C        DO IROOT = 1,NROOT
C          AVEC((IROOT-1)*2*NROOT+IROOT) = 1.0D0
C          APROJR(IROOT*(IROOT+1)/2 )    = EIG(IROOT,ITER)
C        END DO

        DO JVEC = 1,NVEC
          AVEC( 2*NVec*(JVec-1) + JVec ) = 1.0D0
C          APROJR(IROOT*(IROOT+1)/2 )    = EIG(IROOT,ITER)
        END DO

C      call interface_mpi_bcast(AVec(1:2*NVec**2),2*NVec**2,
C     &                         Master,global_communicator)

      IREDL = 0
      IREDL = NVEC*(NVEC+1)/2

      CALL DZERO(SCRRED_R,IREDL)
      CALL DZERO(SCRRED_I,IREDL)

        Do IVec =1,NVec
          JVec_SF = IVec - 1
C      write(6,*) 'myproc ',myproc,'iredl ',iredl
          call interface_mpi_barrier(global_communicator)
          CALL INPROD_B_PAR_CPLX(ILU1,ILU2,VEC1,VEC2,APROJR,APROJI,
     &                           NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                           MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                           LU2LIST,IVEC)
C        write(6,*) 'myproc,dim,normr',myproc,iredl,dnrm2(iredl,aprojr,1)
C        write(6,*) 'myproc,dim,normi',myproc,iredl,dnrm2(iredl,aproji,1)
C        write(6,*) 'myproc ',myproc
C        CALL PRSYM(APROJR,NVEC,LUWRT)
C        CALL PRSYM(APROJI,NVEC,LUWRT)
        End do



*     ... sum up
      CAll redvec(APROJR,SCRRED_R,IREDL,2,op_MPI_SUM,
     &                global_communicator,-1)
      CAll redvec(APROJI,SCRRED_I,IREDL,2,op_MPI_SUM,
     &                global_communicator,-1)

*     ... update APROJR and APROJI 
      CALL DCOPY(IREDL,SCRRED_R,1,APROJR,1)

      CALL DSCAL(IREDL, -1.0D0, SCRRED_I,1)

      CALL DCOPY(IREDL,SCRRED_I,1,APROJI,1)
C        write(6,*) 'Projected matrix after reset'
C        write(6,*) 'Real part'
C        CALL PRSYM(APROJR,NVEC,LUWRT)
C        write(6,*) 'Imaginary part'
C        CALL PRSYM(APROJI,NVEC,LUWRT)

C#ifdef LUCI_DEBUG
C         Allocate(Eig_TEST(2*NVEC))
C         Call DZero(EIG_TEST,2*NVEC)
C         Allocate(AVec_Test(4*NVec**2))
C
C      CALL CEIGEN_NEW(APROJR,APROJI,NVEC,EIG_TEST,AVEC_TEST,
C     &                 AMATC,CWORK,
C     &                RWORK,EVL,IWORK)
C         write(6,*) 'New eigenvalues:'
C         call wrtmat(Eig_Test,1,2*nvec,1,2*nvec)
C         write(6,*) 'New eigenvectors :'
C         call wrtmat(Avec_Test,2*nvec,nvec,2*nvec,nvec)
C         Deallocate(Eig_Test,AVec_Test)
C#endif LUCI_DEBUG


C      Else !Converged or reached maxit
C            CALL IZERO(LU3LIST,IALL_LU3)
C*
C*       c vectors from ILU1 to ILU1
C        DO IROOT = 1, NROOT
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU1LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU1_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IROOT,ILU1,ILU3)
C        END DO
C*
C*       no scaling, we should already work in a normalized basis
C        CALL IZERO(LU1LIST,NROOT*MY_ACT_BLK2)
C*
C        DO IROOT = 1, NROOT
C          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU1LIST,
C     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                         MY_LU3_OFF,MY_LU1_OFF,IROOT,ILU3,ILU1)
C          WORK(2*NROOT*NVEC+IROOT) = 1.0D0
C        END DO
C*
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C*       corresponding sigma vectors from ILU2 to ILU2
C        DO IROOT = 1, NROOT
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU2LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU2_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IROOT,ILU2,ILU3)
C        END DO
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C*       corresponding sigma vectors from ILU2 to ILU2
C        DO IROOT = 1, NROOT
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU2LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU2_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IROOT,ILU2,ILU3)
C        END DO
C*
C*       no scaling, we should already work in a normalized basis
C        CALL IZERO(LU2LIST,NROOT*MY_ACT_BLK2)
C*
C        DO IROOT = 1, NROOT
C          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU2LIST,
C     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                         MY_LU3_OFF,MY_LU2_OFF,IROOT,ILU3,ILU2)
C        END DO
C*
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C        NVEC = NROOT
C*
C*       reset subspace matrices ...
C*
C        CALL DZERO(AVEC,2*NVEC**2)
C        CALL DZERO(APROJR,NVEC*(NVEC+1)/2)
C        CALL DZERO(APROJI,NVEC*(NVEC+1)/2)
C*
C        DO IROOT = 1,NROOT
C          AVEC((IROOT-1)*2*NROOT+IROOT) = 1.0D0
C          APROJR(IROOT*(IROOT+1)/2 )    = EIG(IROOT,ITER)
C        END DO
C
C        End if !Converged or reached maxit

        timer3 = 0.0D0
        starttimer = interface_MPI_WTIME()
C        end if !Subspace filled, and not converged



C      End if



*       check timing

*
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C*       c vectors from ILU1 to ILU1
C        DO IROOT = 1, NROOT
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU1LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU1_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IROOT,ILU1,ILU3)
C        END DO
C*
C*       no scaling, we should already work in a normalized basis
C        CALL IZERO(LU1LIST,NROOT*MY_ACT_BLK2)
C*
C        DO IROOT = 1, NROOT
C          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU1LIST,
C     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                         MY_LU3_OFF,MY_LU1_OFF,IROOT,ILU3,ILU1)
C          WORK(2*NROOT*NVEC+IROOT) = 1.0D0
C        END DO
C*
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C*       corresponding sigma vectors from ILU2 to ILU2
C        DO IROOT = 1, NROOT
C          DO JVEC = 1, NVEC
C            WORK((JVEC-1)*2+1) =   AVEC((IROOT-1)*2* NVEC + JVEC )
C            WORK((JVEC-1)*2+2) = - AVEC((IROOT-1)*2* NVEC + NVEC + JVEC)
C          END DO
C          CALL P3_B_PAR_CPX(VEC1,VEC2,WORK,AVEC,LU2LIST,LU3LIST,
C     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                      MY_LU2_OFF,MY_LU3_OFF,
C     &                      NVEC,2*NVEC,IROOT,ILU2,ILU3)
C        END DO
C*
C*       no scaling, we should already work in a normalized basis
C        CALL IZERO(LU2LIST,NROOT*MY_ACT_BLK2)
C*
C        DO IROOT = 1, NROOT
C          CALL COPVCD_PP_B_CPX(VEC1,LU3LIST,LU2LIST,
C     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
C     &                         MY_LU3_OFF,MY_LU2_OFF,IROOT,ILU3,ILU2)
C        END DO
C*
C        CALL IZERO(LU3LIST,IALL_LU3)
C*
C        NVEC = NROOT
*
*       reset subspace matrices ...
*

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
     &                                 NPARBLOCK,NUM_BLOCKS,
     &                                 IBLOCKL,global_communicator,LBLK,
     &                                 JROOT,LU1LIST,MY_LU1_OFF,IRC)
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
      timer3 = timer3 + interface_MPI_WTIME() - starttimer
      WALLTSTEP = SECTID(timer3)
!     WRITE(LUWRT,9600) WALLTSTEP
      END IF
*     ^ NVEC + NROOT > MAXVEC or CONVER == .TRUE. or ITER == MAXIT
*
*     timing of this iteration
      WALLITR2 = interface_MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,'(/A,2A)')
     &         '                        Wall time for microiteration',
     &         '      ',WALLTID
*
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*
* ( End of loop over iterations )
*
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         WRITE(LUWRT,'(//A,I4,A)') ' (CMICDV_PAR)           '//
     &           'WARNING: maximum number of micro iterations,',
     &           MAXIT, ', is reached, CI aborted.'
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         WRITE(LUWRT,'(//A)')
     &        ' (CMICDV_PAR)           Micro iterations converged.'
      END IF
*
      do iroot = 1, nroot
        fineig(iroot) = eig(iroot,iter+1)+eigshf
      end do
      WRITE(LUWRT,'(/A,F20.8)')
     &   '                        Core energy        = ',eigshf
      WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &   '                        Active energies    = ',
     &                            EIG(1:NROOT,ITER)
      WRITE(LUWRT,'(A,4F20.8/,(37X,4F20.8))')
     &   '                        Final CI energies  = ',
     &     fineig(1:nroot)

      do iroot = 1, nroot
        IF (RTCNV(IROOT)) THEN
          WRITE(LUWRT,'(/A,i4,a)')
     &    '                        root ',iroot,' ...... converged!'
        else
          WRITE(LUWRT,'(/A,i4,a)')
     &    '                        root ',iroot,' did not converge!'
        end if
      end do
      IPRT = 1
      if (IPRT.gt.0) then
        DO 1600 IROOT = 1, NROOT
          WRITE(LUWRT,'(/A,I3/A/)')
     &    ' Information about convergence for root... ' ,IROOT,
     &    ' ============================================'
          FINEIG(IROOT) = EIG(IROOT,ITER)
          WRITE(LUWRT,1300)
 1300     FORMAT(/'  Summary of iterations '
     +           /'  ----------------------')
          WRITE(LUWRT,1310)
 1310     FORMAT
     &    (/'  Iteration point        Eigenvalue         Residual ')
          If(Conver) then
            DO 1330 I=1,ITER+1
 1330       WRITE(LUWRT,1340) I,EIG(IROOT,I)+EIGSHF,RNRM(IROOT,I)
          Else
            DO 1331 I=1,ITER
 1331       WRITE(LUWRT,1340) I,EIG(IROOT,I)+EIGSHF,RNRM(IROOT,I)
          End if
 1340     FORMAT(7X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      end if
      IPRT = 0
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

        ! for +Q correction
        open(file="energies.CI",unit=10,status="unknown",
     &       form="unformatted",access="sequential")
!       TODO: check proper reference energies for states other than the ground state and MCSCF reference wave functions!!!
        write(10) nroot,eig(1,1)+eigshf
        write(10) ((eig(i,iter)+eigshf),i=1,nroot)
        close(10,status="keep")
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
      SUBROUTINE TRNMAD(A,X,SCR,NDIMI,NDIMO)
*
* Obtain X(T) A X and store it in A
* Allows different dimensions in input and output matrices
*     Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION A(*), X(NDIMI,NDIMO)
*. Scratch
      DIMENSION SCR(NDIMI*NDIMO)
      NTEST = 000
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from TRNMAD '
        WRITE(6,*) '   NDIMI, NDIMO = ', NDIMI,NDIMO
        WRITE(6,*) ' Input X matrix '
        CALL WRTMAT(X,NDIMI,NDIMO,NDIMI,NDIMO)
        WRITE(6,*) ' Input A matrix '
        CALL WRTMAT(A,NDIMI,NDIMI,NDIMI,NDIMI)
       END IF

*
*. 1 : X(T) A in SCR
      ZERO = 0.D0
      CALL SETVEC(SCR,ZERO,NDIMI*NDIMO)
      CALL MATML7(SCR,X,A,NDIMO,NDIMI,NDIMI,NDIMO,NDIMI,NDIMI,
     &              0.0D0,1.0D0,1)
*. X(T) A X in A
      CALL MATML7(A,SCR,X,NDIMO,NDIMO,NDIMO,NDIMI,NDIMI,NDIMO,
     &            0.0D0,1.0D0,0)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformed matrix : '
        CALL WRTMAT(A,NDIMO,NDIMO,NDIMO,NDIMO)
      END IF
*
      RETURN
      END
      SUBROUTINE TRNMAD_CPX(A,X,SCR,NDIMI,NDIMO)
*
* Obtain X(H) A X and store it in A
* Allows different dimensions in input and output matrices
*     Complex version of TrnMad. September 2019, Andreas Nyvang
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION A(*), X(2*NDIMI,NDIMO)
*. Scratch
      DIMENSION SCR(2*NDIMO,NDIMI)
      NTEST = 0000
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from TRNMAD '
        WRITE(6,*) '   NDIMI, NDIMO = ', NDIMI,NDIMO
        WRITE(6,*) ' Input X matrix '
        CALL WRTMAT(X,2*NDIMI,NDIMO,2*NDIMI,NDIMO)
        WRITE(6,*) ' Input A matrix '
        CALL WRTMAT(A,2*NDIMI,NDIMI,2*NDIMI,NDIMI)
       END IF

*
*. 1 : X(H) A in SCR
      ZERO = 0.D0

      LDX = max(1,2*NDimi)
      LDA = max(1,2*NDimi)
      LDScr = max(1,2*NDimo)

      !Real part
      Call DGemm ('T','N',NDimo,NDimi, 2*NDimi, 1.0D0, X, LDX, A, LDA,
     &            0.0D0, Scr, LDScr )

      !Imaginary part
      Call DGemm ( 'T', 'N', NDimo, NDimi, NDimi, 1.0D0, X, LDX,
     &  A(NDimi+1), LDA, 0.0D0, Scr(NDimo+1,1), LDScr )


      Call DGemm ( 'T', 'N', NDimo, NDimi, NDimi, -1.0D0, X(NDimi+1,1), 
     &             LDX, A, LDA, 1.0D0, Scr(NDimo+1,1), LDScr )


*. (X^H A)X in A
      LDA = max(1,2*NDimo) !Note the difference from prior LDA
C      write(6,*) 'NDimi ',NDimi
C      write(6,*) 'NDimo ',NDimo
C      write(6,*) 'LDScr ',LDScr
C      write(6,*) 'LDX ',LDX
C      write(6,*) 'LDA ',LDA

      !Real part
      Call DGemm ('N','N',NDimo,NDimo,NDimi,1.0D0,Scr(1,1),LDScr,
     &            X,LDX,0.0D0,A,LDA)

      Call DGemm ('N','N',NDimo,NDimo,NDimi,-1.0D0,Scr(NDimo+1,1),
     &            LDScr,X(NDimi+1,1),LDX,1.0D0,A,LDA)

      !Imaginary part
      Call DGemm ('N','N',NDimo,NDimo,NDimi,1.0D0,Scr(1,1),LDScr,
     &            X(NDimi+1,1),LDX,0.0D0,A(NDimo+1),LDA)

      Call DGemm ('N','N',NDimo,NDimo,NDimi,1.0D0,Scr(NDimo+1,1),
     &            LDScr,X,LDX,1.0D0,A(NDimo+1),LDA)


      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformed matrix : '
        CALL WRTMAT(A,2*NDIMO,NDIMO,2*NDIMO,NDIMO)
      END IF

      RETURN
      END
      SUBROUTINE CHECK_UNIT_MAT(UNI,NDIM,XMAX_DIFF_DIAG,XMAX_DIFF_OFFD,
     &           ISYM)
*
* A matrix UNI is given. Check difference between UNI and UNIT matrix
* and report in:
*     XMAX_DIFF_DIAG: Max absolute difference between between diagonal
*                     element and 1
*     XMAX_DIFF_OFFD: Max absolute difference between off diagonal and zero
*
*. Jeppe Olsen, July 2011 (Thought I had written this routine before...)
*  Last modification; Feb 27, 2013; Jeppe Olsen; ISYM added
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION UNI(*)
*. Diagonal element
      XMAX_DIFF_DIAG = 0.0D0
      DO I = 1, NDIM
         IF(ISYM.EQ.0) THEN
           II = (I-1)*NDIM + I
         ELSE
          II = I*(I-1)/2 + I
        END IF
        XMAX_DIFF_DIAG = MAX(XMAX_DIFF_DIAG,ABS(UNI(II)-1.0D0))
      END DO
*. Off diagonal elements
      XMAX_DIFF_OFFD = 0.0D0
      DO I = 1, NDIM
       DO J = 1, I-1
        IF(ISYM.EQ.0) THEN
          JI = (I-1)*NDIM + J
          IJ = (J-1)*NDIM + I
          XMAX_DIFF_OFFD = MAX(XMAX_DIFF_OFFD,ABS(UNI(IJ)),ABS(UNI(JI)))
        ELSE
         IJ = I*(I-1)/2 + J
         XMAX_DIFF_OFFD = MAX(XMAX_DIFF_OFFD,ABS(UNI(IJ)))
        END IF
       END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Comparison of matrix with unit matrix: '
       WRITE(6,*) '   Largest deviation of diagonal elements ',
     & XMAX_DIFF_DIAG
       WRITE(6,*) '   Largest deviation of of-diagonal elements ',
     & XMAX_DIFF_OFFD
      END IF
      NTEST = 00
*
      RETURN
      END
      SUBROUTINE MGS4(X,S,NDIM,SCR1,THRES,NVECUT)
*
* Modified Gram-Schmidt procedure by forward orthogonalization
*
*  watch out for zero columns indicating linear dependency
*
* Jeppe Olsen, March 2013, added thres to MGS3
*
* S is input overlap matrix, X is output set of orthonormalized vectors
*
* Thres is min norm of linear independent vector- only meaningfull if all
* initial vectors have identical norm
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. input
      DIMENSION S(NDIM,NDIM)
*. Output
      DIMENSION X(NDIM,*)
*. Scratch : vector of length NDIM
      DIMENSION SCR1(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from MGS4 '
        WRITE(6,*) ' ==============='
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Initial overlap matrix'
       CALL WRTMAT(S,NDIM,NDIM,NDIM,NDIM)
      END IF
*
*. Initialize X to unit matrix
*
      ZERO = 0.0D0
      ONE = 1.0D0
      CALL SETVEC(X,ZERO,NDIM**2)
      CALL SETDIA(X,ONE,NDIM,0)
C          SETDIA(MATRIX,VALUE,NDIM,IPACK)
*
      DO IVEC = 1, NDIM
*. Normalize vector IVEC
C        CALL MATVCB(S,X(1,IVEC),SCR1,NDIM,NDIM,0)
        Call DCopy(NDim,S(1,IVec),1,Scr1,1)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
*. avoid NaN's by putting norm to at least zero
        XNORM = INPROD(X(1,IVEC),SCR1,NDIM)
        IF(NTEST.GE.100) THEN
          WRITE(6,'(A,I4,E15.7)') ' IVEC, XNORM = ', IVEC, XNORM
        END IF
*
        IF (XNORM.LE.THRES) THEN
          FACTOR = 0.0D0
        ELSE
          FACTOR = 1.0D0/SQRT(XNORM)
        END IF
        CALL SCALVE(X(1,IVEC), FACTOR, NDIM)
        CALL SCALVE(SCR1,FACTOR,NDIM)
*. Subtract X(1,IVEC) from all remaining vectors
        DO JVEC = IVEC+1,NDIM
          XSX = INPROD(SCR1,X(1,JVEC),NDIM)
          CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,IVEC),ONE,-XSX,NDIM)
        END DO
      END DO
*
*. And remove zero vectors
*
      NVECUT = 0
      DO IVEC = 1, NDIM
        XNORM = INPROD(X(1,IVEC),X(1,IVEC),NDIM)
        IF(XNORM.GT.0.0D0) THEN
          NVECUT = NVECUT + 1
          IF(NVECUT.NE.IVEC) CALL COPVEC(X(1,IVEC),X(1,NVECUT),NDIM)
        END IF
      END DO
*
      IF(NTEST.GE.1.AND.NVECUT.NE.NDIM) THEN
        WRITE(6,*)' MGS4 reduced dim, from and to ', NDIM,NVECUT
      ELSE IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization information:'
        WRITE(6,*) ' Number of linear independent vectors ', NVECUT
        WRITE(6,*) ' Orthonormalized vectors '
        CALL WRTMAT(X,NDIM,NVECUT,NDIM,NVECUT)
      END IF
*
      NTEST = 00
      RETURN
      END
      SUBROUTINE GET_CNEWCOLD_BAS(CN,CNO,NVEC,NROOT,SCR,NVECUT,
     &           IPRT,THRES)
*
* A subspace expansion is given by CN
* Obtain orthogonal expansion of basis for New + old vectors
*
*. Jeppe Olsen, March 2013- reducing orthogonalization errors
*
*. Last modification; Jeppe Olsen; Aug 6, 2015; Checking that
*                     eigenvectors contains corrections.
*
      INCLUDE 'implicit.inc'
!      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 INPROD
      INTEGER*8 :: NVEC,NROOT,IRC,KS,KX,KFREE,NTEST,IPRT
*. Input
      DIMENSION CN(NVEC,NROOT)
*. Output
      DIMENSION CNO(NVEC,*)
*. Scratch: Min size: MAX(4*2*NROOT*NVEC,4*(2*NROOT)**2)
      DIMENSION SCR(*)
*
*. Old roots or corrections as second set of vectors
      IRC = 2
*. Partitioning of SCR
      KS = 1
      KFREE = KS + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KX = KFREE
      KFREE = KX + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KSCR = KFREE
      KFREE = KSCR + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KVEC = KFREE
      KFREE = KVEC + 2*NROOT
*
      ZERO = 0.0D0
      ONE = 1.0D0
*
      NTEST = 100
      NTEST = MAX(NTEST,IPRT)

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GET_CNEWCOLD_BAS '
        WRITE(6,*) ' ==========================='
        WRITE(6,*)
        WRITE(6,'(A,2I4)') ' NROOT, NVEC = ', NROOT, NVEC
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' The initial CN matrix '
       CALL WRTMAT(CN,NVEC,NROOT,NVEC,NROOT)
      END IF
*
*. CNO: First Nroot vectors are the new vectors,
*. the last NROOT are the original first roots(IRC = 1) or corrections
*. IRC = 2
      THRES2 = 1.0D-11
      ZERO = 0.0D0
      CALL SETVEC(CNO,ZERO,2*NROOT*NVEC)
      CALL COPVEC(CN,CNO,NROOT*NVEC)
      NCORR = 0
      IF(IRC.EQ.1) THEN
*. The original first roots
        DO IROOT = 1, NROOT
          CNO(IROOT,NROOT+IROOT) = 1.0D0
        END DO
      ELSE
*. Change the last NROOT vectors to corrections
COLD    CALL COPVEC(CN,CNO(1,NROOT + 1),NROOT*NVEC)
        DO IROOT = 1, NROOT
          XORT = INPROD(CN(NROOT+1,IROOT),CN(NROOT+1,IROOT),NVEC-NROOT)
          XORT = SQRT(XORT)
          IF(XORT.GT.THRES2) THEN
            NCORR = NCORR + 1
            CALL COPVEC(CN(1,IROOT),CNO(1,NROOT+NCORR),NVEC)
            CNO(IROOT,NROOT+NCORR) = CNO(IROOT,NROOT+NCORR) - 1.0D0
            XNORM = INPROD(CNO(1,NROOT+NCORR),CNO(1,NROOT+NCORR),
     &                     NVEC)
            SCALE = 1.0D0/SQRT(XNORM)
            CALL SCALVE(CNO(1,NROOT+IROOT),SCALE,NVEC)
          END IF
        END DO
C        DO IROOT = 1, NROOT
C          XORT = INPROD(CN(NROOT+1,IROOT),CN(NROOT+1,IROOT),NVEC-NROOT)
C          XORT = SQRT(XORT)
C          IF(XORT.GT.THRES2) THEN
C            NCORR = NCORR + 1
C            CALL COPVEC(CN(1,IROOT),CNO(1,IROOT),NVEC)
C            CNO(IROOT,IROOT) = CNO(IROOT,IROOT) - 1.0D0
C            XNORM = INPROD(CNO(1,NROOT+NCORR),CNO(1,NROOT+NCORR),
C     &                     NVEC)
C            SCALE = 1.0D0/SQRT(XNORM)
C            CALL SCALVE(CNO(1,NROOT+IROOT),SCALE,NVEC)
C          END IF
C        END DO
      END IF
      NDO = NROOT+NCORR
*
* Iterate over orthogonalizations
*
      MAXIT = 2
      DO ITER = 1, MAXIT
        IF(NTEST.GE.100) WRITE(6,*) ' Info from orth. iter. ', ITER
        IF(ITER.GT.1) THEN
          NDO = NVECUT
        END IF
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The CNO matrix '
          CALL WRTMAT(CNO,NVEC,NDO,NVEC,NDO)
        END IF
*. Overlap matrix S = CNO^T CNO
C         MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                    NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        CALL MATML7(SCR(KS),CNO,CNO,NDO,NDO,NVEC,NDO,NVEC,NDO,ZERO,ONE,
     &              1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Overlap of NDO basis '
          CALL WRTMAT(SCR(KS),NDO,NDO,NDO,NDO)
        END IF
*. Perform Gram-Schmidt orthogonalization
C       MGS4(X,S,NDIM,SCR1,THRES,NVECUT)
        CALL MGS4(SCR(KX),SCR(KS),NDO,SCR(KVEC),THRES,NVECUT)
C   TRNMAD(A,X,SCR,NDIMI,NDIMO)
*.   Check that delivered X orthonormalizes
        CALL TRNMAD(SCR(KS),SCR(KX),SCR(KSCR),NDO,NVECUT)
        IF(NTEST.GE.100) THEN
             WRITE(6,*) ' New overlap matrix '
           CALL WRTMAT(SCR(KS),NVECUT,NVECUT,NVECUT)
        END IF
        IF(IPRT.GE.10)
     &  CALL CHECK_UNIT_MAT(SCR(KS),NVECUT,XMAX_DIA,XMAX_OFD,0)
*. Transformation from original to orthonormal basis
        CALL MATML7(SCR(KS),CNO,SCR(KX),NVEC,NVECUT,
     &              NVEC,NDO,NDO,NVECUT,ZERO,ONE,0)
        CALL COPVEC(SCR(KS),CNO,NVEC*NVECUT)
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,'(A,2I5)')
     &  ' Output from GET_CNEW.., NROOT and NVECUT = ',
     &  NROOT, NVECUT
      END IF
*
      IF(NTEST.GE.1.AND.NVECUT.LT.2*NROOT) THEN
        WRITE(6,'(A,2I5)')
     &  ' GET_CNEW.., Reduced number of vectors: 2*NROOT and NVECUT = ',
     &  2*NROOT, NVECUT
      END IF

*
      RETURN
      END
      SUBROUTINE SETDIA(MATRIX,VALUE,NDIM,IPACK)
*
* Set diagonal elements of matrix MATRIX to VALUE
*
* IPACK = 0 => full quadratic matrix
* IPACK = 1 => lower triangular matrix, row packed
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MATRIX(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I=1,NDIM
100     MATRIX((I-1)*NDIM+I) = VALUE
      ELSE IF (IPACK .EQ. 1 ) THEN
        DO 200 I = 1, NDIM
 200    MATRIX(I*(I+1)/2) = VALUE
      ELSE
        WRITE(6,*) ' IPACK called with IPACK = ', IPACK
        STOP ' SETDIA ,IPACK out of range '
      END IF
*
      RETURN
      END
