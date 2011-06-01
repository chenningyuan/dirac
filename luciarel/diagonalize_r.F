!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!dirac_copyright_end

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE CIEIG5_REL(INICI,EROOTS,VEC1,
     &           VEC2,MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           NDIM,NBLK,NROOTS,MAXVEC,MXCIIT,MXROOTS,LUINCI,
     &           IPRT,NPRDET,IPNTR,NP1,NP2,NQ,
     &           SCR1,ISCR1,SCR2,T_CC,T_BUFF,
     &           EIGSHF,ICISTR,LBLK,IDIAG,THRES_E,
     &           NBATCH,LBATCHB,LBATCHE,
     &           IBLOCK,INIDEG,E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,
     &           NCLS,
     &           IRC,IDCOMH,IHAM12,SIGDEN_ROUTE,
     &           IBLOCKL,NPARBLOCK,IPROCLIST,IGROUPLIST
#if defined (VAR_MPI2)
     &           ,LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,LUCLIST,
     &           RCCTOS,IT_TTPL,IT_TTOL
#endif
     &     )
*
* Master routine for CI diagonalization
*
* Modified to handle PQ - preconditioner , May 1990
* PICO,MICDV4 added spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
      DIMENSION VEC1(*), VEC2(*)
      DIMENSION INIDET(100)
      DIMENSION IPNTR(*), ISCR1(*)
      DIMENSION SCR1(*), SCR2(*)
      DIMENSION IBLOCKL(*), NPARBLOCK(*) 
      DIMENSION IPROCLIST(*),IGROUPLIST(*)
#include "infpar.h"
#if defined (VAR_MPI2)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION IT_TTPL(*)
      INTEGER   RCCTOS(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL(*)
#include "cstate.inc"
#endif
#include "parluci.h"
*
*     output from subspace diagonalization
*
*
      DIMENSION EROOTS(NROOTS)
*
      IF( IPRT.GT. 1 )  WRITE(LUWRT,'(/A)')
     &'          *** information from ci diagonalization  ***'
      WRITE(LUWRT,*) ' '
*
      NTEST=0
*
CSK#if defined (VAR_MPI2)
CSK      WRITE(LUWRT,*)'  THIS IS MY LUCLIST'
CSK      CALL IWRTMAMN( LUCLIST, 1, IALL_LUC, 1, IALL_LUC, LUWRT)
CSK      WRITE(LUWRT,*)'  '
CSK#endif
*               ====================================
** 1 :               INITIAL VARIATIONAL SUBSPACE
*               ====================================
*
      IF( INICI .EQ. 0 ) THEN
        IF(NPRDET .EQ. 0 ) THEN
* ==================================================
*. Initial guess from lowest elements of CI diagonal
* ==================================================
* In order to treat degeneracies, the lowest 4 * NROOT elements are
*.obtained

          NFINDM = MIN(NDIM,4*NROOTS)
#if defined (VAR_MPI2)
          CALL FNDMND_PAR_REL(LUDIA,LBLK,VEC1,NFINDM,NFINDA,
     &                        ISCR1(1+2*NFINDM),SCR1(1+2*NFINDM),ISCR1,
     &                        SCR1,IBLOCKL,NPARBLOCK,NUM_BLOCKS,IPRT)
C
#else
          CALL FNDMND(LUDIA,LBLK,VEC1,NFINDM,NFINDA,ISCR1(1+2*NFINDM),
     &                SCR1(1+2*NFINDM),ISCR1,SCR1,IPRT)
C              FNDMND(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,SCR,ISCAT,SUBVAL
#endif
          CALL REWINE(LU1,-1)
          IBASE = 1
          TEST = 1.0D-10
*
          DO 100 IROOTS = 1, NROOTS
*. Number of degenerate elements
            NDEG = 1
            XVAL = SCR1(IBASE)
   90       CONTINUE
            IF(IBASE-1+NDEG+1.LE.NFINDA) THEN
              IF (ABS(SCR1(IBASE-1+NDEG+1)-XVAL).LE.TEST) THEN
                NDEG = NDEG + 1
                GOTO 90
              END IF
            END IF
            WRITE(LUWRT,*) ' IROOT NDEG ', IROOTS,NDEG
*
            IF (INIDEG.EQ.0.AND.NDEG.GT.1) THEN
              WRITE(LUWRT,*) ' WARNING WARNING WARNING WARNING ! '
              WRITE(LUWRT,*) ' DEGENERATE INITIAL VECTORS FOR CI '
              NDEG = 1
            END IF

*. Initial guess in compressed form in SCR1
            SCALE = 1.0D0/SQRT(DFLOAT(NDEG))
            DO 250 II = 1,NDEG
*. Anti symmetric combination
              IF(INIDEG.EQ.-1) THEN
                SCR1(II) = (-1.0D0)**II * SCALE
*. Symmetric combination
              ELSE IF (INIDEG.EQ.1.OR.INIDEG.EQ.0) THEN
                SCR1(II) =  SCALE
              END IF
  250       CONTINUE
C                WRSVCD_REL(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM)
* Real part
            CALL REWINE(LUDIA,LBLK)
c           ISCR1(1) = 16
            IAM_REAL = 1
            CALL WRSVCD_REL(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
     &                      LUDIA,NUM_BLOCKS,IBLOCKL,NPARBLOCK,IAM_REAL,
     &                      IROOTS
#if defined (VAR_MPI2)
     &                      ,LU1LIST
#endif
     &                      )
            IAM_REAL = 2
*. Imaginary part
            if (IRC.eq.2) then
              CALL WRSVCD_REL(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,0,NDIM,
     &                        LUDIA,NUM_BLOCKS,IBLOCKL,NPARBLOCK,
     &                        IAM_REAL,IROOTS
#if defined (VAR_MPI2)
     &                        ,LU1LIST
#endif
     &                        )
          
            end if
            IBASE = IBASE + NDEG
  100     CONTINUE
        END IF
      END IF
*
*
*                 ========================
* 2 :                  Diagonalization
*                 ========================
*
      CALL QENTER('CIEIG5')
* Inverse iteration modified Davidson with 2 vectors in core
      IF(IPRT .GE. 2 ) THEN
         WRITE(LUWRT,*)
         WRITE(LUWRT,'(A,I3)')
     &   '  Number of roots to be converged..  ',NROOTS
         WRITE(LUWRT,*)
         WRITE(LUWRT,'(A,I3)')
     &   '  Largest allowed number of vectors..',MAXVEC
         WRITE(LUWRT,*)
         WRITE(LUWRT,'(A,I3)')
     &   '  Allowed number of CI iterations  ..',MXCIIT
      END IF
*
      MRNRM = 1
      IFAC = 1
      IF(IRC.EQ.2) IFAC = 2
C
      MEIG = MRNRM + (MXCIIT+1)*NROOTS
      MAPROJ = MEIG  + (MXCIIT+1)*2*(MAXVEC + NROOTS)
      IF(IRC.EQ.2) THEN
        MAPROJI = MAPROJ +  MAXVEC*(MAXVEC+1)/2
      ELSE
        MAPROJI = MAPROJ
      END IF
      MAVEC = MAPROJI + MAXVEC*(MAXVEC+1)/2
      MWORK = MAVEC + MAXVEC ** 2*IFAC*IFAC
      MLFRSC = MWORK + MAXVEC*(MAXVEC+1)*IFAC*IFAC
C
C     new complex hermitian matrix diagonalizer
      IF( IRC .eq. 2)THEN
C     +1 added for safety reasons
C     double complex!
      LCWORK  = IRC * (2*MAXVEC + MAXVEC**2 + 1)
C     +1 added for safety reasons
      LRWORK = 1 + 5*MAXVEC + 2*MAXVEC**2 + 1
C     +1 added for safety reasons
      LIWORK = 3 + 5*MAXVEC + 1
C     double complex!
      LAMATC = IRC * MAXVEC * MAXVEC
      LEIGVL = MAXVEC
      END IF
C
*
      IF( IPRT .GE. 5 ) THEN
         WRITE(LUWRT,*) ' MRNRM MEIG MAPROJ MAVEC MWORK MLFRSC '
         WRITE(LUWRT,'(6I4)')  MRNRM,MEIG,MAPROJ,MAVEC,MWORK,MLFRSC
      END IF
*
*     proper allocation of scratch arrays:
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'CIEIG5')
      call memmar(KRNRM,MEIG-1,'ADDL  ',2,'KRNRMS')
      call memmar(KEIG,MAPROJ-MEIG,'ADDL  ',2,'KEIGS ')
      if (IRC.eq.1) then
        call memmar(KAPROJ,(MAXVEC*(MAXVEC+1)/2),'ADDL  ',
     &              2,'KAPROJ')
      else
        call memmar(KAPROJ,MAPROJI-MAPROJ,'ADDL  ',2,'KAPROJ')
        call memmar(KAPROJI,MAVEC-MAPROJI,'ADDL  ',2,'KAPROI')
C
C       new stuff for complex hermitian matrix diagonalizer
        call memmar(KAMATC,LAMATC,'ADDL  ',2,'KAMATC')
        call memmar(KCWORK,LCWORK,'ADDL  ',2,'KCWORK')
        call memmar(KRWORK,LRWORK,'ADDL  ',2,'KRWORK')
        call memmar(KEIGVL,LEIGVL,'ADDL  ',2,'KEIGVL')
        call memmar(KIWORK,LIWORK,'ADDL  ',1,'KIWORK')
C
C       ... initialize
        CALL DZERO(WORK(KAMATC),LAMATC)
        CALL DZERO(WORK(KCWORK),LCWORK)
        CALL DZERO(WORK(KRWORK),LRWORK)
        CALL DZERO(WORK(KEIGVL),LEIGVL)
        CALL IZERO(WORK(KIWORK),LIWORK)
      end if
      call memmar(KAVEC,MWORK-MAVEC,'ADDL  ',2,'KAVECS')
      call memmar(KWORK,MLFRSC-MWORK,'ADDL  ',2,'KWORKS')
      
#if defined (VAR_MPI2)
*
*     scratch array for reducing operations
      MAXSCRRED = MAX( (MAXVEC*(MAXVEC+1)/2) , MAXVEC**2 )
      CALL MEMMAR(KSCRRED,MAXSCRRED,'ADDL  ',2,'KREDUC')
      IF( IRC .eq. 2 )THEN
        CALL MEMMAR(KSCRRED2,MAXSCRRED,'ADDL  ',2,'KREDU2')
      END IF
*
*.    Arrays for partitioning of the CI vector
*     (*_MS2 arrays are in CSTATE )
      JCMBSPC = 1
      call z_blkfo_rel_par(IDC,NMS2VAL,JCMBSPC,IREFSM,
     &                     KLVLBT,KLVLEBT,KLVI1BT,KLVIBT,KDMXXX,
     &                     NBATCH,NUM_BLOCKS,NBLK_MS2,IBLK_MS2,
     &                     NBAT_MS2,IBAT_MS2,0,1,NPARBLOCK)
CSK      call z_blkfo_rel(IDC,NMS2VAL,JCMBSPC,IREFSM,
CSK     &                 KLVLBT,KLVLEBT,KLVI1BT,KLVIBT,KDUM,
CSK     &                 NBATCH,NUM_BLOCKS,NBLK_MS2,IBLK_MS2,
CSK     &                 NBAT_MS2,IBAT_MS2,0,0,NPARBLOCK)
*
#endif
*
*     Construct and diagonalize complete Hamiltonian 
      IF( IDCOMH .eq. 1 )THEN
#if defined(VAR_MPI2)
        IF(NMPROC .ge. 1 ) THEN
          WRITE(LUWRT,*)  ' I am sorry, but this route does not work'
          WRITE(LUWRT,*)  ' at the moment'
          CALL Abend2('STOP for IDCOMH == 1 ')
        END IF
#endif
#if defined CODE_NOT_WORKING
        CALL COMHAM(VEC1,VEC2,NDIM,T_CC,SIGDEN_ROUTE)
#endif
      END IF
*
      IF (IRC.EQ.1) THEN
C
C       REAL DAVIDSON solver
#if defined (VAR_MPI2)
*
        call micdv4_enlmd_rel_par(T_CC,T_BUFF,VEC1,VEC2,
     &                  WORK(KRNRM),WORK(KEIG),
     &                  EROOTS,MXCIIT,NROOTS,
     &                  MAXVEC,NROOTS,WORK(KAPROJ),WORK(KAVEC),
     &                  WORK(KWORK),IPRT,NPRDIM,
     &                  IPNTR,LBLK,EIGSHF,THRES_E,
     &                  IBLOCKL,NPARBLOCK,NUM_BLOCKS,
     &                  LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LUCLIST,NBATCH,
     &                  WORK(KLVLBT),WORK(KLVLEBT),WORK(KLVI1BT),
     &                  WORK(KLVIBT),RCCTOS,
     &                  WORK(KSCRRED),IGROUPLIST,IPROCLIST,
     &                  IT_TTPL,IT_TTOL)
*
#else
*
        call micdv4_enlmd_rel(T_CC,T_BUFF,VEC1,VEC2,LU1,LU2,
     &                  WORK(KRNRM),WORK(KEIG),
     &                  EROOTS,MXCIIT,LU3,LU4,LU5,LUDIA,NROOTS,
     &                  MAXVEC,NROOTS,WORK(KAPROJ),WORK(KAVEC),
     &                  WORK(KWORK),IPRT,NPRDIM,
     &                  H0DUM,IPNTR,NP1,NP2,NQ,H0SCRDUM,
     &                  LBLK,EIGSHF,THRES_E)
*
#endif
        IROOTHOMING = 0
      ELSE
C
C     Complex DAVIDSON solver 
#if defined (VAR_MPI2)
C
        CALL CMICDV_PAR(VEC1,VEC2,WORK(KRNRM),
     &                  WORK(KEIG),SCR2(1),MXCIIT,NDIM,
     &                  NROOTS,MAXVEC,NROOTS,
     &                  T_CC,T_BUFF,WORK(KAPROJI),WORK(KAPROJ),
     &                  WORK(KAVEC),
     &                  WORK(KWORK),IPRT,NPRDIM,IPNTR,
     &                  LBLK,EIGSHF,IHAM12,SIGDEN_ROUTE,
     &                  THRES_E,IBLOCKL,NPARBLOCK,NUM_BLOCKS,
     &                  LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LUCLIST,NBATCH,
     &                  WORK(KLVLBT),WORK(KLVLEBT),WORK(KLVI1BT),
     &                  WORK(KLVIBT),RCCTOS,WORK(KSCRRED),
     &                  WORK(KSCRRED2),IGROUPLIST,IPROCLIST,
     &                  IT_TTPL,IT_TTOL,WORK(KAMATC),WORK(KCWORK),
     &                  WORK(KRWORK),WORK(KEIGVL),WORK(KIWORK))
*
#else
*
        CALL CMICDV(VEC1,VEC2,LU1,LU2,WORK(KRNRM),WORK(KEIG),
     &              SCR2(1),MXCIIT,
     &              NDIM,LU3,LU4,LU5,LU6,LUDIA,NROOTS,MAXVEC,NROOTS,
     &              T_CC,T_BUFF,
     &              WORK(KAPROJI),WORK(KAPROJ),WORK(KAVEC),WORK(KWORK),
     &              IPRT,NPRDIM,IPNTR,NP1,NP2,NQ,LBLK,EIGSHF,
     &              IHAM12,SIGDEN_ROUTE,THRES_E,WORK(KAMATC),
     &              WORK(KCWORK),WORK(KRWORK),WORK(KEIGVL),WORK(KIWORK))
*
#endif
      END IF
*     ^ End if switch real/complex numbers
      CALL MEMMAR(KDUM ,IDUM,'FLUSM ',2,'CIEIG5')
*
      CALL QEXIT('CIEIG5')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CEIGEN_NEW(HR,HI,NDIM,EIGVAL,EIGVEC,AMATC,CWORK,
     &                      RWORK,EVL,IWORK)
C***********************************************************************
C
C     find the eigenvalues and eigenvectors using 
C     the general complex hermitian solver routine ZHEEVD
C
C     a complex hermitian matrix (double prec) is given in the form of 
C     a real part HR and an imaginary part HI
C
C
C     OUTPUT: 
C             - eigenvalues in EIGVAL
C             - eigenvectors in EIGVEC where
C               real coefficients are stored from EIGVEC(1   ,  IVEC) 
C                                              to EIGVEC(NDIM,  IVEC)
C
C               and imaginary coefficients   from EIGVEC(NDIM+1,IVEC)
C                                                 EIGVEC(2*NDIM,IVEC)
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C     INPUT
      DIMENSION HR(NDIM*(NDIM+1)/2),HI(NDIM*(NDIM+1)/2)
C     OUTPUT
      DIMENSION EIGVAL(2*NDIM)
      DIMENSION EIGVEC(2*NDIM,2*NDIM)
C     scratch for LAPACK solver
      INTEGER            IWORK(*)
      DOUBLE PRECISION   RWORK(*), EVL(*)
      COMPLEX*16         AMATC(NDIM,*), CWORK(*), ABUFF, A0
      CHARACTER*1 JOBZ, UPLO

      NTEST = 00
      IF(NTEST.GT.00) THEN
        WRITE(LUPRI,'(/A )') ' ======================='
        WRITE(LUPRI,'( A )') '  WELCOME to CEIGEN_new '
        WRITE(LUPRI,'( A/)') ' ======================='
      END IF
C
      A0 = (0.0D0,0.0D0)
C     initialize
      DO ICOL = 1, NDIM
        DO JROW = 1, NDIM
          AMATC(JROW,ICOL) = A0
        END DO
      END DO
      AREALP = 0.0D0
      AIMAGP = 0.0D0
C
C     fill the lower triangular part of matrix AMATC
      DO IROW = 1, NDIM
        DO ICOL = 1, IROW
          AREALP           = HR(IROW*(IROW-1)/2 + ICOL)
          AIMAGP           = HI(IROW*(IROW-1)/2 + ICOL)
          ABUFF            = DCMPLX(AREALP,AIMAGP)
          AMATC(IROW,ICOL) = ABUFF
#ifdef MOD_DEBUG
          WRITE(LUPRI,*) ' IROW, ICOL, IROW*(IROW-1)/2 + ICOL, HR(..)',
     &    IROW, ICOL, IROW*(IROW-1)/2 + ICOL, HR(IROW*(IROW-1)/2 + ICOL)
          WRITE(LUPRI,*) ' IROW, ICOL, IROW*(IROW-1)/2 + ICOL, HI(..)',
     &    IROW, ICOL, IROW*(IROW-1)/2 + ICOL, HI(IROW*(IROW-1)/2 + ICOL)
          WRITE(LUPRI,*) ' IROW, ICOL, AMATC(IROW,ICOL)',
     &                     IROW, ICOL, AMATC(IROW,ICOL)
#endif
        END DO
      END DO
C
C     diagonalize complex hermitian matrix AMATC(NDIM,NDIM)
C     - get eigenvalues and eigenvectors
      JOBZ   = 'V'
      UPLO   = 'L'
      INFO   =  0
      LADIM  = NDIM
C     +1 added for safety reasons
      LCWORK  = 2*NDIM + NDIM**2 + 1
C     +1 added for safety reasons
      LRWORK = 1 + 5*NDIM + 2*NDIM**2 + 1
C     +1 added for safety reasons
      LIWORK = 3 + 5*NDIM + 1
C
C     call library routine
      CALL ZHEEVD( JOBZ, UPLO, LADIM, AMATC, LADIM, EVL, CWORK, LCWORK,
     &             RWORK, LRWORK, IWORK, LIWORK, INFO )
      IF(INFO.NE.0) THEN
         WRITE(LUPRI,'(/A,I4)') ' *** ERROR in CEIGEN_NEW: info code'//
     &                  ' NOT zero from ZHEEVD:',INFO
         CALL QUIT ('*** Error in CEIGEN_NEW (lapack zheevd) ***')
      ENDIF
C
C     eigenvalues in EVL(NDIM) - copy to EIGVAL
      N2DIM = 2*NDIM
      CALL DZERO(EIGVAL,N2DIM)
      CALL DCOPY(NDIM,EVL,1,EIGVAL,1)
C
C     eigenvectors in AMATC - insert into EIGVEC(N2DIM,N2DIM)
      CALL DZERO(EIGVEC,N2DIM**2)
      AREALP = 0.0D0
      AIMAGP = 0.0D0
      DO I = 1, NDIM
        DO J = 1, NDIM
C          real and imag part
           AREALP           = DBLE(AMATC(J,I))
           AIMAGP           = AIMAG(AMATC(J,I))
           EIGVEC(J,I)      = AREALP
           EIGVEC(NDIM+J,I) = AIMAGP
        END DO
      END DO
C
C     debug print section 
      IF( NTEST .ge. 10 )THEN
        WRITE(LUPRI,'(/A/)') ' *** CEIGEN_NEW ***: final eigenvalues'
        CALL WRTMATMN(EIGVAL,1,N2DIM,1,N2DIM,LUPRI)
        WRITE(LUPRI,'(/A/)') ' *** CEIGEN_NEW ***: final eigenvectors'
        CALL WRTMATMN(EIGVEC,N2DIM,N2DIM,N2DIM,N2DIM,LUPRI)
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
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
#include "ipoist8.inc"
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
      CALL DCOPY(NDIM*(NDIM+1)/2,HR,1,SCR(1),1)
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
        CALL DCOPY(IROWI,HR(IROWI*(IROWI-1)/2+1),1,SCR(IROWOF+NDIM),1)
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
          CALL DZERO(EIGVEC(1,IOFF+1),2*NDIM)
        ELSE
*. more than 2 degenerate eigenvectors, eliminate by Gram orthogonalization
*. (Gram-Schmidt)
          DO IVEC = IOFF,IOFF+LDEG-1
*. | Vec I > = Scalar * ( | Vec I > - Sum(J<I) <Vec J | Vec I > | Vec J > )
*. Orthogonalize
            DO JVEC = IOFF, IVEC-1
*. < Vec J | Vec I >
               XR = DDOT(NDIM,EIGVEC(1,JVEC),1,EIGVEC(1,IVEC),1)
     &            + DDOT(NDIM,EIGVEC(1+NDIM,JVEC),1,EIGVEC(1+NDIM,IVEC)
     &              ,1)
               XI = DDOT(NDIM,EIGVEC(1,JVEC),1,EIGVEC(1+NDIM,IVEC),1)
     &            - DDOT(NDIM,EIGVEC(1+NDIM,JVEC),1,EIGVEC(1,IVEC),1)
              CALL VECSUMC(EIGVEC(1,IVEC),EIGVEC(1,IVEC),
     &                     EIGVEC(1,JVEC),
     &                     ONE,ZERO,-XR,-XI,NDIM)
             END DO
*. Normalize
             XNORM = DDOT(2*NDIM,EIGVEC(1,IVEC),1,EIGVEC(1,IVEC),1)
        if(ntest.ge.10) WRITE(6,*) ' ldeg XNORM = ',ldeg, XNORM
             IF(ABS(XNORM).LE.THRES) THEN
*. This is a zero vector
               CALL DZERO(EIGVEC(1,IVEC),2*NDIM)
             ELSE
*. This is not a zero vector
               SCALE = ONE/SQRT(XNORM)
               CALL DSCAL(2*NDIM,SCALE,EIGVEC(1,IVEC),1)
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
        XNORM = DDOT(2*NDIM,EIGVEC(1,IVEC),1,EIGVEC(1,IVEC),1)
        if(ntest.ge.10) WRITE(6,*) ' XNORM = ',XNORM
        IF(XNORM.NE.0.0D0) THEN
          IEFF = IEFF + 1
          EIGVAL(IEFF) = EIGVAL(IVEC)
          CALL DCOPY(2*NDIM,EIGVEC(1,IVEC),1,EIGVEC(1,IEFF),1)
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
      END
***********************************************************************

      SUBROUTINE CMICDV(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  T_CC,T_BUFF,
     &                  APROJI,APROJR,AVEC,WORK,IPRT,
     &                  NPRDIM,IPNTR,NP1,NP2,NQ,LBLK,EIGSHF,
     &                  IHAM12,SIGDEN_ROUTE,
     &                  CONVLIM,AMATC,CWORK,RWORK,EVL,IWORK)
*
* Multi root Davidson routine for complex Hermitian
* matrices, requires two blocks in core
*
* May 1994, Jeppe Olsen
*
* Initial version without general preconditioner
* or inverse iteration corrections
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
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
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
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
      REAL * 8   INPROD,INPRDD
      character*6 SIGDEN_ROUTE
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12, FILELAB*8
      LOGICAL CONVER,RTCNV(NROOT),EXSTOP
*
#include "krmcluci_inf.h"
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION RNRM(NROOT,MAXIT),EIG(2*MAXVEC,MAXIT)
      DIMENSION APROJI(*),APROJR(*),AVEC(*),WORK(*)
      DIMENSION IPNTR(1)
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
C
* Dimensioning required of local vectors
*      APROJR : MAXVEC*(MAXVEC+1)/2
*      APROJI : MAXVEC*(MAXVEC+1)/2
*      AVEC   : (2*MAXVEC) ** 2
*      WORK   : 2*MAXVEC*(2*MAXVEC+1)/2
*
      DIMENSION FINEIG(*)
*
      call qenter('CMICDV')
*
      NTEST = 0000
      IPRT = MAX(NTEST,IPRT)
*
       WRITE(LUWRT,*) '                  :::::::::::::::::::'
       WRITE(LUWRT,*) '                    Entering CMICDV  '
       WRITE(LUWRT,*) '                  :::::::::::::::::::'
*
*
      if (NTEST.ge.1) then
        write(LUWRT,*)
        write(LUWRT,*) '----------------------------------------------'
        write(LUWRT,*) ' I´ll take you to a place '
        write(LUWRT,*) ' Where we shall find our'
        write(LUWRT,*) '  ROOTS, BLOODY ROOTS '
        write(LUWRT,*)
        write(LUWRT,*) '               Sepultura '
        write(LUWRT,*) '----------------------------------------------'
        write(LUWRT,*)
      end if
*
      write(LUWRT,*)
      write(LUWRT,*) ' Nuclear repulsion + core energy: ',EIGSHF
*
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(LUWRT,*) ' Sorry CMICDV wounded , MAXVEC .LT. 2*NROOT '
         WRITE(LUWRT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(LUWRT,*) ' Increase MXCIV_CI to be at least 2 * Nroot '
         Call Abend1(' Enforced stop in CMICDV ')
       END IF
*
*  Complex diagonalizer, therefore
       IRC = 2
*
       CONVER = .FALSE.
*
       IF(MAXIT.EQ.0) THEN
         WRITE(LUWRT,*) ' MAXIT = 0, so just a fast visit '
         RETURN
       END IF
*
       if (IPRT.ge.1) then
         IF( MYPROC .eq. MASTER ) THEN
           WRITE(LUWRT,*) ' Diagonal in CMICDV'
           WRITE(LUWRT,*) ' =================='
           CALL REWINE(LUDIA,LBLK)
           CALL WRTVCDC(VEC1,LUDIA,1,IRC,LBLK)
         END IF
       end if
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       NVEC = 0
       IF( MAXIT .GE. 1 ) THEN
*
       IF(IPRT.GE.1000) THEN
         IF( MYPROC .eq. MASTER ) THEN
           WRITE(LUWRT,*) ' Initial vectors '
           WRITE(LUWRT,*) ' ================'
           CALL REWINE(LU1,LBLK)
           DO IVEC = 1, NINVEC
             WRITE(LUWRT,*) ' Initial vector ...  ', IVEC
             CALL WRTVCDC(VEC1,LU1,0,IRC,LBLK)
           END DO
         END IF
       END IF
*
         CALL REWINE(LU1,LBLK)
         CALL REWINE(LU2,LBLK)
         DO 10 IVEC = 1,NINVEC
           CALL REWINE(LU3,LBLK)
           CALL REWINE(LU4,LBLK)
*. Copy to scratch unit LU3
           CALL COPVCDC(LU1,LU3,VEC1,0,IRC,LBLK)
*
*--------------------------------------------------------------
* Generalization to SIGDEN routine block implemented.
*--------------------------------------------------------------
           call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,T_BUFF,1)
*
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
           CALL REWINE(LU4,LBLK)
           CALL COPVCDC(LU4,LU2,VEC1,0,IRC,LBLK)
*. Projected matrix
           DO 8 JVEC = 1, IVEC
             IJ = IVEC*(IVEC-1)/2 + JVEC
             CALL INPRDDC(VEC1,VEC2,LU2,LU3,JVEC,1,1,LBLK,XR,XI)
C                 INPRDDC(VEC1,VEC2,LU1,LU2,IOFF1,IOFF2,IREW,LBLK,XR,XI)
             APROJR(IJ) = XR
             APROJI(IJ) = -XI
    8      CONTINUE
   10    CONTINUE
*
         IF(IPRT.GE.1000) THEN
           WRITE(LUWRT,*) ' Initial sigma vectors '
           WRITE(LUWRT,*) ' ======================'
           CALL REWINE(LU2,LBLK)
           DO IVEC = 1, NINVEC
             WRITE(LUWRT,*) ' Initial vector ...  ', IVEC
             CALL WRTVCDC(VEC1,LU2,0,IRC,LBLK)
           END DO
         END IF
*
         IF( IPRT .GE.10) THEN
           WRITE(LUWRT,*) ' Initial projected matrix, real part'
           CALL PRSYM(APROJR,NINVEC)
           WRITE(LUWRT,*) ' Initial projected matrix, imag part'
           CALL PRSYM(APROJI,NINVEC)
          END IF
C
C     Diagonalize initial projected matrix
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
C
         IF(IPRT .GE. 3) THEN
           WRITE(LUWRT,'(A,I4)') ' Eigenvalues of initial iteration '
           WRITE(LUWRT,'(F22.13)')
     &     ( EIG(IROOT,1)+EIGSHF,IROOT=1,NROOT)
         END IF
         NVEC = NINVEC
       END IF
*
       IF (MAXIT.EQ.1) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
       CALL GETTIM(CPUITR1,WALLITR1)
       write(LUWRT,*)
       WRITE(LUWRT,*) ' Info from iteration .... ', ITER
       WRITE(LUWRT,*) ' _________________________'
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
         EIGAPR = EIG(IROOT,ITER)
C        EIGAPR = EIG(IROOT,ITER-1)
*. ==========================
*. Real part of HX - EIGAPR*X
*. ==========================
         DO IVEC = 1, NVEC
           WORK(2*(IVEC-1)+1) = AVEC((IROOT-1)*2*NVEC+IVEC)
           WORK(2*(IVEC-1)+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+IVEC)
         END DO
*. Real part, Hx
         CALL MVCSMD(LU2,WORK,LU3,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
C             MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
*. Real part X
         CALL MVCSMD(LU1,WORK,LU5,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
*
C?       WRITE(6,*) 
C?   &   ' Current approximation to real part of eigenvector'
C?       CALL WRTVCD(VEC1,LU5,1,LBLK)
*. Real part of Hx -Ex on LU4
         EIGAPRM = - EIGAPR
         ONE = 1.0D0
         CALL REWINE(LU3,LBLK)
         CALL REWINE(LU5,LBLK)
         CALL REWINE(LU4,LBLK)
         CALL VECSMD(VEC1,VEC2,ONE,EIGAPRM,LU3,LU5,LU4,0,LBLK)
*. ==========================
*. Imag part of HX - EIGAPR*X
*. ==========================
         DO IVEC = 1, NVEC
           WORK(2*(IVEC-1)+1) = AVEC((IROOT-1)*2*NVEC+NVEC+IVEC)
           WORK(2*(IVEC-1)+2) = AVEC((IROOT-1)*2*NVEC     +IVEC)
         END DO
*. Imag part of HX
         CALL MVCSMD(LU2,WORK,LU3,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
C               MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
*. Imag part X
         CALL MVCSMD(LU1,WORK,LU5,LU6,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imag part of Hx -Ex as second vector on LU4
         EIGAPRM = - EIGAPR
         ONE = 1.0D0
         CALL REWINE(LU3,LBLK)
         CALL REWINE(LU5,LBLK)
         CALL VECSMD(VEC1,VEC2,ONE,EIGAPRM,LU3,LU5,LU4,0,LBLK)
         IF ( IPRT  .GE. 600 ) THEN
           WRITE(LUWRT,*) '  ( HX - EX )'
           CALL WRTVCDC(VEC1,LU4,1,IRC,LBLK)
         END IF
*  Strange place to put convergence but ....
         CALL REWINE(LU4,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,0,LBLK)
     &         +       INPRDD(VEC1,VEC1,LU4,LU4,0,LBLK) )
         RNRM(IROOT,ITER) = RNORM
*
* Flush output buffer
         CALL FLSHFO(6)
*
* Set logical for convergence results
         IF(RNORM.LT. CONVLIM ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         WRITE(LUWRT,'(A19,5X,I5,3X,1E18.13,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER,RNORM,EIGAPR+EIGSHF
* I suppose this is only a check for the necessity of defining
* a new search direction. Using intermediate ITERIM instead ...
CTF      ITER = ITER + 1
         ITERIM = ITER + 1
         IF( ITERIM .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
*. Inverted diagonal times (HX-EX) on LU3
*. Real part
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU4,LU3,-EIGAPR,1,1,LBLK)
*. Imag part
           CALL REWINE(LUDIA,LBLK)
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU4,LU3,-EIGAPR,0,1,LBLK)
           IF ( IPRT  .GE. 600) THEN
             WRITE(LUWRT,*) '(D-E)-1 *( HX - EX ) '
             CALL WRTVCDC(VEC1,LU3,1,IRC,LBLK)
           END IF
* =====================================================================
*. 1.3 Orthogonalize to all previous vectors
* =====================================================================
           CALL REWINE(LU3,LBLK)
           XNORM =  INPRDD(VEC1,VEC1,LU3,LU3,0,LBLK)
     &           +  INPRDD(VEC1,VEC1,LU3,LU3,0,LBLK)
           XNORMO = XNORM
           DO JVEC = 1,NVEC+IADD-1
C    9            INPRDDC(VEC1,VEC2,LU1,LU2,IOFF1,IOFF2,IREW,LBLK,XR,XI)
             CALL INPRDDC(VEC1,VEC2,LU1,LU3,JVEC,1,1,LBLK,
     &                    XR,XI)
             WORK(2*(JVEC-1)+1) = XR
             WORK(2*(JVEC-1)+2) = XI
             XNORM = XNORM - XR**2 - XI**2
           END DO
           THRSLDP = 1.0D-10
           IF(XNORM/XNORMO .LE. THRSLDP) THEN
*. Linear dependency
             WRITE(LUWRT,*)
     &       ' CMICDV : Vector eliminated due to linear dependence '
             IADD = IADD -1
             GOTO 100
           END IF
           SCALE = 1.0D0/SQRT(XNORM)
C          write(6,*) ' Scale, scale**2*xnorm',scale,
C    &                scale**2*xnorm
           CALL DSCAL( 2*(NVEC+IADD-1), SCALE, WORK, 1 )
*. Real part of orthogonalized vector on LU5
           DO JVEC = 1, NVEC+IADD-1
             WORK(2*(JVEC-1)+2) =
     &      -WORK(2*(JVEC-1)+2)
           END DO
           CALL MVCSMD(LU1,WORK,LU4,LU5,VEC1,VEC2,2*(NVEC+IADD-1),
     &                 1,LBLK)
C               MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
           CALL REWINE(LU3,LBLK)
           CALL REWINE(LU4,LBLK)
           CALL REWINE(LU5,LBLK)
           ONEM = -1.0D0
           CALL VECSMD(VEC1,VEC2,SCALE,ONEM,LU3,LU4,LU5,0,LBLK)
*. Imag part of orthogonal vector
           DO JVEC = 1, NVEC+IADD-1
             BUF = WORK(2*(JVEC-1)+1)
             WORK(2*(JVEC-1)+1) = - WORK(2*(JVEC-1)+2)
             WORK(2*(JVEC-1)+2) = BUF
           END DO
           CALL MVCSMD(LU1,WORK,LU4,LU6,VEC1,VEC2,2*(NVEC+IADD-1),
     &                 1,LBLK)
C               MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
           CALL SKPVCD(LU3,1,VEC1,1,LBLK)
           CALL REWINE(LU4,LBLK)
           ONEM = -1.0D0
           CALL VECSMD(VEC1,VEC2,SCALE,ONEM,LU3,LU4,LU5,0,LBLK)
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(LUWRT,*) '   Orthonormalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCDC(VEC1,LU5,1,IRC,LBLK)
           END IF
C               SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
           CALL SKPVCD(LU1,2*(NVEC+IADD-1),VEC1,1,LBLK)
           CALL REWINE(LU5,LBLK)
           CALL COPVCDC(LU5,LU1,VEC1,0,IRC,LBLK)
         END IF
  100 CONTINUE
CTF
      ITER = ITER + 1
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
* =======================================================
**  2 : Optimal combination of new and old directions
* =======================================================
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,2*NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,2*NVEC,VEC1,1,LBLK)
      DO IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCDC(LU1,LU3,VEC1,0,IRC,LBLK)
        CALL REWINE(LU4,LBLK)
        IF(IPRT.GE.10000) THEN
          WRITE(LUWRT,*) ' Input vector to mV7 '
          CALL WRTVCDC(VEC1,LU3,1,IRC,LBLK)
        END IF
*
* Generalization to SIGDEN routine block implemented.
        call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,T_BUFF,1)
*
        IF(IPRT.GE.10000) THEN
          WRITE(LUWRT,*) ' Output vector from sigden_ctrl '
          CALL WRTVCDC(VEC1,LU4,1,IRC,LBLK)
        END IF
        CALL REWINE(LU4,LBLK)
        CALL COPVCDC(LU4,LU2,VEC1,0,IRC,LBLK)
*. Augment projected matrix
        DO  JVEC = 1, NVEC+IVEC
          IJ = IVEC*(IVEC-1)/2 + JVEC
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          CALL INPRDDC(VEC1,VEC2,LU3,LU2,1,JVEC,1,LBLK,XR,XI)
C               INPRDDC(VEC1,VEC2,LU1,LU2,IOFF1,IOFF2,IREW,LBLK,XR,XI)
          APROJR(IJ) = XR
          APROJI(IJ) = XI
        END DO
      END DO
*
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      IF( IPRT  .GE. 10 ) THEN
        WRITE(LUWRT,*) ' Projected matrices before CEIGEN '
        CALL PRSYM(APROJR,NVEC)
        CALL PRSYM(APROJI,NVEC)
      END IF
      WRITE(LUWRT,*) 'Dimension of projected matrix:',NVEC
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
CJO     call towave(nroot,eig(1,iter))
      END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRT,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJR,NVEC)
        CALL PRSYM(APROJI,NVEC)
        CALL WRTMAT(AVEC,2*NVEC,NROOT,2*NVEC,NROOT)
      END IF
  901 CONTINUE
*
* =============================================
**  reset or assemble eigenvectors
* =============================================
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
*. Reset C vectors
        CALL REWINE( LU5,LBLK)
        DO IROOT = 1, NROOT
*. Real part of vector IROOT on LU3
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+JVEC)
           WORK((JVEC-1)*2+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
          END DO
          CALL MVCSMD(LU1,WORK,LU3,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imaginary part of vector on LU6
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
           WORK((JVEC-1)*2+2) = AVEC((IROOT-1)*2*NVEC+JVEC)
          END DO
          CALL MVCSMD(LU1,WORK,LU6,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
     &          + INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
C          write(LUWRT,*) ' Scale, scale**2*xnorm',scale,
C    &                scale**2*xnorm
          WORK(2*NROOT*NVEC+IROOT) = SCALE
*. Copy current vector to LU5
C?        WRITE(LUWRT,*) ' Current vector on LU3, IROOT = ', IROOT
C?        CALL WRTVCD(VEC1,LU3,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
*
C?        WRITE(LUWRT,*) ' Current vector on LU6, IROOT = ', IROOT
C?        CALL WRTVCD(VEC1,LU6,1,LBLK)
*
          CALL REWINE(LU6,LBLK)
          CALL SCLVCD(LU6,LU5,SCALE,VEC1,0,LBLK)
CE        IF(IPRT.GT.100) THEN
CE           WRITE(LUWRT,*) ' Normalized vector '
CE           call rewine(LU5,LBLK)
CE           CALL WRTVCDC(VEC1,LU5,1,IRC,LBLK)
CE        END IF
        END DO
*. Transfer C vectors to LU1
*
        if (IPRT.ge.1000) then
          WRITE(LUWRT,*) ' Initial vectors on LU1 '
          CALL REWINE(LU1,LBLK)
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) ' Vector ', IVEC
            CALL WRTVCDC(VEC1,LU1,0,IRC,LBLK)
            WRITE(LUWRT,*) ' Writing of vector completed '
          END DO
*
          WRITE(LUWRT,*) ' Vectors on LU5 '
          CALL REWINE(LU5,LBLK)
          DO IVEC = 1, NROOT
            WRITE(LUWRT,*) ' Vector ', IVEC
            CALL WRTVCDC(VEC1,LU5,0,IRC,LBLK)
            WRITE(LUWRT,*) ' Writing of vector completed '
          END DO
        end if
*
        CALL REWINE(LU1,LBLK)
        CALL REWINE(LU5,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCDC(LU5,LU1,VEC1,0,IRC,LBLK)
C?        WRITE(6,*) ' Vector copied : ', IVEC
        END DO
        IF(IPRT.GE.1000) THEN
          WRITE(LUWRT,*) ' New reset C-vectors '
          WRITE(LUWRT,*) ' =================== '
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCDC(VEC1,LU1,0,IRC,LBLK)
          END DO
        END IF
*. corresponding sigma vectors
        CALL REWINE(LU5,LBLK)
        DO IROOT = 1, NROOT
*. Real part of vector IROOT on LU3
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+JVEC)
           WORK((JVEC-1)*2+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
          END DO
          CALL MVCSMD(LU2,WORK,
     &    LU3,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imaginary part of vector on LU6
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
           WORK((JVEC-1)*2+2) = AVEC((IROOT-1)*2*NVEC+JVEC)
          END DO
          CALL MVCSMD(LU2,WORK,
     &    LU6,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
          SCALE = WORK(2*NROOT*NVEC+IROOT)
*. Copy to LU5
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
          CALL REWINE(LU6,LBLK)
          CALL SCLVCD(LU6,LU5,SCALE,VEC1,0,LBLK)
        END DO
*. Transfer HC vectors to LU2
        CALL REWINE(LU2,LBLK)
        CALL REWINE(LU5,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCDC(LU5,LU2,VEC1,0,IRC,LBLK)
        END DO
        IF(IPRT.GE.1000) THEN
          WRITE(LUWRT,*) ' New reset sigma -vectors '
          WRITE(LUWRT,*) ' =================== '
          CALL REWINE(LU2,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCDC(VEC1,LU2,0,IRC,LBLK)
          END DO
        END IF
        NVEC = NROOT
*. Reset projected matrices
        CALL SETVEC(AVEC,0.0D0,2*NVEC**2)
        CALL SETVEC(APROJR,0.0D0,NVEC*(NVEC+1)/2)
        CALL SETVEC(APROJI,0.0D0,NVEC*(NVEC+1)/2)
        DO IROOT = 1,NROOT
          AVEC((IROOT-1)*2*NROOT+IROOT) = 1.0D0
          APROJR(IROOT*(IROOT+1)/2 ) = EIG(IROOT,ITER)
        END DO
*
      END IF
*
      CALL GETTIM(CPUITR2,WALLITR2)
      CPUTID = SECTID(CPUITR2-CPUITR1)
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,9300) CPUTID,WALLTID
*
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE

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
        CALL REWINE(LU1,LBLK)
        CALL REWINE(LU2,LBLK)
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
          IF(IPRT.GE.10) THEN
            WRITE(LUWRT,1201) IROOT
 1201       FORMAT(/' The final approximation to eigenvector',I5)
            CALL WRTVCDC(VEC1,LU1,0,IRC,LBLK)
 1200       FORMAT(/' The final sigma vector')
            CALL WRTVCDC(VEC1,LU2,0,IRC,LBLK)
          END IF
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
      call qexit('CMICDV')
*
* Delete scratch units
      close(unit=LU3,status='DELETE')
      close(unit=LU4,status='DELETE')
      close(unit=LU5,status='DELETE')
      close(unit=LU6,status='DELETE')
*
      RETURN
 1030 FORMAT(/3X,7F15.8,/,(3X,7F15.8))
 1120 FORMAT(/3X,I3,7F15.8,/,(6X,7F15.8))
 9300 FORMAT(' >>>  CPU (WALL) TIME IN ITERATION: ',A,'(',A,')')
      END
***********************************************************************

      SUBROUTINE MICDV6_REL(
     &                  T_CC,T_BUFF,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,
     &                  MAXIT,
     &                  LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,NPRDIM,
     &                  H0,IPNTR,NP1,NP2,NQ,H0SCR,
     &                  LBLK,EIGSHF,E_CONV,IROOTHOMING)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version
*
* From MICDV4, Jeppe Olsen, April 1997
*              Modified : Oct 97 : root homing added
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
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
CTF
CTF  Preliminary version without H0 
CTF  Simple inverse diagonal Hessian update
CTF
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
       dimension T_CC(*), T_BUFF(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
*      IROOTHOMING : Do roothoming, i.e. select the
*      eigenvectors in iteration n+1 as the approximations
*      with largest overlap with the previous space
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(NROOT)
       REAL*8 INPRDD, INPROD
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*. Notice XJEP is also used for ROOTHOMING, should be allocated
*  outside (for roothoming :dim = 3*MAXVEC )
       DIMENSION XJEP(10000)
       INTEGER   IXJEP(10000)
*
       WRITE(6,*) '             ::::::::::::::::::::::::   '
       WRITE(6,*) '               Entering MICDV6 (rel)    '
       WRITE(6,*) '             ::::::::::::::::::::::::   '
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
c      IPRT = 6000
       IPRT = IPRTXX
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV_CI to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV6 '
         Call Abend1( 20 )
       END IF
*
C?     WRITE(6,*) ' LU1 LU2 LU3 LU4 LU5 LU6 LU7',
C?   &              LU1,LU2,LU3,LU4,LU5,LU6,LU7
       IF(IROOTHOMING.EQ.1) THEN
         WRITE(6,*) ' Root homing performed '
       END IF
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINE(LU1,-1)
       CALL REWINE(LU2,-1)
       DO IVEC = 1,NINVEC
         CALL REWINE(LU5,-1)
         CALL REWINE(LU6,-1)
         CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
CTF      CALL MV7(VEC1,VEC2,LU5,LU6)
         call sigden_ctrl(VEC1,VEC2,LU5,LU6,T_CC,T_BUFF,1)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINE(LU6,-1)
         CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINE(LU2,-1)
         DO JVEC = 1, IVEC
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
         END DO
*        ^ End of loop over JVEC
       END DO
*      ^ End of loop over IVEC
*
       IF( IPRT .GE.5 ) THEN
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
         WRITE(6,*) ' Eigenvalues of initial iteration (with shift)'
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(6,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
C     TRAVCD(VEC1,VEC2,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)
*
       IF(IPRT.GE.600) THEN
         WRITE(6,*) ' Initial set of eigenvectors '
         CALL REWINE(LU1,-1)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(6,*) ' Initial set of sigma vectors '
         CALL REWINE(LU2,-1)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       ZERO = 0.0D0
       CALL SETVEC(APROJ,ZERO,NROOT*(NROOT+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO
*
       NVEC = NROOT
       IF (MAXIT .EQ. 1 ) GOTO  1001
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        CALL GETTIM(CPUITR1,WALLITR1)
        write(6,*)
        write(6,'(A21,3X,I3)') ' Info from iteration ',ITER
        write(6,*) '_______________________'
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
       CALL REWINE(LU1,-1)
       CALL REWINE(LU2,-1)
*
       DO 100 IROOT = 1, NROOT
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINE(LU7,-1)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
*
         CALL REWINE(LU7,-1)
         CALL REWINE(LU5,-1)
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
         WRITE(6,'(A19,7X,I3,3X,1E18.13,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
*
         CALL FLSHFO(6)
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
*. (H0-E)-1 * (HX-EX) on LU6
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(6,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINE(LUDIA,-1)
           CALL REWINE(LU5,-1)
           CALL REWINE(LU6,-1)
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
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
             CALL REWINE(LUDIA,-1)
*
             CALL H0M1TD(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
             GAMMA = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
             IF(IPRT.GE.1000) WRITE(6,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
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
*
*. 1.3 Orthogonalize to all previous vectors
*.. Vectors on LU1
           CALL REWINE( LU1 ,LBLK)
           DO IVEC = 1, NROOT
               CALL REWINE(LU6,-1)
               WORK(IVEC) = -INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
           END DO
           IF(IPRT.GE.1000) THEN
             Write(6,*) ' Overlap with vectors on LU1'
             call wrtmat(work,1,nroot,1,nroot)
           END IF
           ONE = 1.0D0
           CALL REWINE(LU1,-1)
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
                 CALL REWINE(LU6,-1)
                 WORK(IVEC) = -INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
             END DO
             ONE = 1.0D0
             CALL REWINE(LU3,-1)
             CALL MVCSMD(LU3,WORK,LU7,LU6,VEC1,VEC2,NVEC+IADD-1-NROOT,
     &            1,LBLK)
             CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU5,LU6,1,LBLK)
           ELSE
             CALL REWINE(LU3,-1)
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
*
CJAN20IF( CONVER ) GOTO  901
      IF( CONVER ) THEN
*. Well, we converged
C        ITER = ITER-1
         GOTO  1001
      END IF
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
        CALL REWINE(LU3,-1)
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
* Check timings for sigma vector generation
C?        CALL GETTIM(CPUSIG1,WALLSIG1)
CTF     CALL MV7(VEC1,VEC2,LU5,LU6)
        call sigden_ctrl(VEC1,VEC2,LU5,LU6,T_CC,T_BUFF,1)
C?        CALL GETTIM(CPUSIG2,WALLSIG2)
C?        CPUTID = SECTID(CPUSIG2-CPUSIG1)
C?        WALLTID = SECTID(WALLSIG2-WALLSIG1)
C?        WRITE(6,9400) CPUTID,WALLTID
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
      write(6,*) 'Dimension of projected matrix:',NVEC
      CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
*. Test : transform projected matrix
C TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
C     CALL TRAN_SYM_BLOC_MAT(APROJ,AVEC,1,NVEC,XJEP(1000),XJEP(1))
C     WRITE(6,*) ' Explicitly transformed matrix'
C     CALL PRSYM(XJEP(1000),NVEC)

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
      IF(IROOTHOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first
*
*. Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
       END DO
       WRITE(6,*)
     & ' Norm of projections to previous vector space '
       CALL WRTMAT(XJEP,1,NVEC,1,NVEC)
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL SCALVE(XJEP,ONEM,NVEC)
       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
       IF(NSORT.LT.NVEC) THEN
         WRITE(6,*) ' Warning : Some elements lost in sorting '
         WRITE(6,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
       IF(IPRT.GE.3) THEN
         WRITE(6,*) ' New roots choosen as vectors '
         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
       END IF
*. Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
       END DO
       CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
       END DO
       DO INEW = 1, NVEC
         WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
       END DO
*
       IF(IPRT.GE.3) THEN
         WRITE(6,*) ' Reordered WORK and AVEC arrays '
         CALL PRSYM(WORK,NVEC)
         CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
       END IF
*
      END IF
*     ^ End of root homing procedure
      DO IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
      END DO
*
      IF(IPRT .GE. 3 ) THEN
        WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(6,'(5F18.13)')
     &  ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
        WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
        WRITE(6,'(5F18.13)')
     &  ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
      END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
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
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' Initial vectors on LU1'
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
          WRITE(6,*) ' Initial vectors on LU3'
          CALL REWINE(LU3,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF

C       IF(NVEC.GT.NROOT) THEN
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU1 to LU5
          CALL REWINE(LU1,-1)
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of c reset, part 1'
*. and add LU3 to LU5
          CALL REWINE(LU3,-1)
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
C         WRITE(6,*) ' SCALE = ', SCALE
          WORK(IROOT) = SCALE
*. scale LU5 => LU7
          CALL REWINE(LU5,-1)
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
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' Reset C-vectors on LU1 '
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
        END IF
*. 3.2 : corresponding sigma vectors
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU2 to LU5
          CALL REWINE(LU2,-1)
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 1'
*. and add LU4 to LU5
          CALL REWINE(LU4,-1)
          CALL MVCSMD2(LU4,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 2'
          SCALE  = WORK(IROOT)
*. scale LU5 => LU7
          CALL REWINE(LU5,-1)
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
*. Overlap between first vectors on LU1 and LU2
        E11 = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
C?      WRITE(6,*) ' <Lu1!Lu2> = ', E11
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
          CALL REWINE(LU5,-1)
          DO JVEC = 1, IADD
*. Orthogonalize to vectors on LU1
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU3,LU7,VEC1,0,LBLK)
            IF(IPRT.GE.1000) THEN
            write(6,*) ' Initial vector on LU7'
            CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
            CALL REWINE(LU1,-1)
            DO IROOT = 1, NROOT
              CALL REWINE(LU7,-1)
              WORK(IROOT+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU1,0,LBLK)
            END DO
            IF(IPRT.GE.1000)
     &      write(6,*) ' c, part1 finito'
*. And to trial vectors on LU5
            CALL REWINE(LU5,-1)
            DO KVEC = 1, JVEC-1
              CALL REWINE(LU7,-1)
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
            CALL REWINE(LU7,-1)
            CALL SCLVCD(LU7,LU5,SCALE,VEC1,0,LBLK)
C                SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
           IF(IPRT.GE.1000)
     &       write(6,*) ' c, part6 finito'
         END DO
*        /\ End of loop over orthogonalized directions
*. Transfer modified directions to LU3
         CALL REWINE(LU3,-1)
         CALL REWINE(LU5,-1)
         DO JVEC =1, IADD
           CALL COPVCD(LU5,LU3,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000) THEN
           write(6,*) ' c, part7 finito'
           WRITE(6,*) ' Additional trial vectors on LU3'
           CALL REWINE(LU3,-1)
           DO JVEC = 1, IADD
            CALL WRTVCD(VEC1,LU3,0,LBLK)
           END DO
         END IF
* Sigma vectors corresponding to orthogonalized directions
         CALL SKPVCD(LU4,ISTART-1,VEC1,1,LBLK)
         CALL REWINE(LU5,-1)
         DO JVEC = 1, IADD
            CALL REWINE(LU7,-1)
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
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU7,LU5,VEC1,0,LBLK)
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 3 finito '
         END DO
*        /\ End of loop over orthogonalized directions
*. Copy back to LU4
         CALL REWINE(LU4,-1)
         CALL REWINE(LU5,-1)
         DO JVEC = 1, IADD
            CALL COPVCD(LU5,LU4,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000)
     &   write(6,*) ' s, part 4 finito '
         NNVEC = NROOT + IADD
       END IF
       IF(IPRT.GE.1000) THEN
         WRITE(6,*) ' Additional sigma vectors on LU4'
         CALL REWINE(LU4,-1)
         DO JVEC = 1, IADD
           CALL WRTVCD(VEC1,LU4,0,LBLK)
         END DO
       END IF
*
*      /\ End if more than NROOT vectors to be reset
       NVEC = NNVEC
C      END IF
*.     ^ End if reset
*. New subspace
*. Calculate subspace Hamiltonian from actual vectors
*. on disc
       IF(IPRT.GE.1000) write(6,*) ' Subspace hamiltonian'
       CALL REWINE(LU1,-1)
       CALL REWINE(LU3,-1)
       DO IVEC = 1, NVEC
*
         CALL REWINE(LU5,-1)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINE(LU2,-1)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
         END DO
         CALL REWINE(LU4,-1)
         DO JVEC = NROOT+1,IVEC
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
         END DO
       END DO
       if (IPRT.ge.2) then
         write(6,*) ' Reset hamiltonian'
         call prsym(aproj,nvec)
       end if
*. Test : Orthogonality of new vectors
       CALL REWINE(LU1,-1)
       CALL REWINE(LU3,-1)
       DO IVEC = 1, NVEC
*
         CALL REWINE(LU5,-1)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINE(LU1,-1)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU1,0,LBLK)
         END DO
         CALL REWINE(LU3,-1)
         DO JVEC = NROOT+1,IVEC
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU3,0,LBLK)
         END DO
       END DO
       if (IPRT.ge.2) then
         write(6,*) ' Overlap matrix    '
         call prsym(xjep,nvec)
       end if
*
*  Timing of this iteration
       CALL GETTIM(CPUITR2,WALLITR2)
       CPUTID = SECTID(CPUITR2-CPUITR1)
       WALLTID = SECTID(WALLITR2-WALLITR1)
       WRITE(6,9300) CPUTID,WALLTID
*. End of resetting business
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE

* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT(/'  Convergence was NOT obtained in',I4,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(/'  Convergence was obtained in',I4,' iterations')
        END IF
*
      IF ( IPRT .GT. 0 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          write(6,*)
          write(6,*) '------------------------'
          write(6,*) 'Root number  ',IROOT
          write(6,*) '------------------------'
          DO I=1,ITER
            WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          end do
          WRITE(6,'(/A,I4/A/)')
     &    ' Information about convergence for root...' ,IROOT,
     &    '==============================================='
          FINEIG(IROOT) = EIG(ITER,IROOT)
          IF (RTCNV(IROOT)) THEN
             WRITE(6,1190) IROOT,FINEIG(IROOT)+EIGSHF
          ELSE
             WRITE(6,1191) IROOT,FINEIG(IROOT)+EIGSHF
          END IF
 1190 FORMAT(' The final eigenvalue',I5,F22.10,' (converged)')
 1191 FORMAT(' The final eigenvalue',I5,F22.10,' (NOT converged)')
          IF(IPRT.GE.400) THEN
            WRITE(6,1200) IROOT
 1200       FORMAT(/' The final approximation to eigenvector',I5)
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
 1340     FORMAT(I10,F28.13,E14.5)
 1600   CONTINUE
      ELSE
        write(6,*)
 1310   FORMAT
     &  (/'  Iteration point        Eigenvalue         Residual ')
        DO IROOT = 1, NROOT
          write(6,*) '------------------------'
          write(6,*) 'Root number  ',IROOT
          write(6,*) '------------------------'
          DO I=1,ITER
            WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          end do
        end do
        write(6,*)
        write(6,*) '**********************************'//
     &             '**************************'
        write(6,*) '   Iter  Root       Energy        '//
     &             'RESIDUAL     RESRATIO '
        write(6,*) '**********************************'//
     &             '**************************'
        write(6,*)
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           WRITE(6,'(3X,I3,2X,I3,3X,F18.10,2X,E10.3,3X,E10.3)')
     &               ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &               RNRM(1,IROOT)/RNRM(ITER,IROOT)
           write(6,'(A,I5,F20.10)') 
     &               ' Final energy root',IROOT,FINEIG(IROOT)
 1601   CONTINUE
      END IF
*
      CALL FLSHFO(6)
*
      Call Add_Info('E_CI',FINEIG,NROOT)
*
      RETURN
 1030 FORMAT(/,3X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(/,3X,I3,7F15.8,/,(1H ,5X,7F15.8))
 9300 FORMAT(' >>>  CPU (WALL) TIME IN ITERATION: ',A,'(',A,')')
 9400 FORMAT(' >>>  CPU (WALL) TIME IN SIGMA VECTOR CALL: ',A,'(',A,')')
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
*
      SUBROUTINE MICDV4_ENLMD_REL(T_CC,T_BUFF,VEC1,VEC2,LU1,LU2,
     &                  RNRM,EIG,FINEIG,MAXIT,
     &                  LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV)
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
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*),T_CC(*), T_BUFF(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2                               
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(NROOT)
       REAL*8 INPRDD
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12,FILELAB*8
*
       call qenter('MICDV4')
*
       WRITE(LUWRT,*) '             ::::::::::::::::::::::::   '
       WRITE(LUWRT,*) '               Entering MICDV4 (rel)    '
       WRITE(LUWRT,*) '             ::::::::::::::::::::::::   '
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
!      IPRT = 5000 ! debug
!      IPRT = 0000

       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINE(LU1,-1)
       CALL REWINE(LU2,-1)
       DO 10 IVEC = 1,NINVEC
         CALL REWINE(LU3,-1)
         CALL REWINE(LU4,-1)
         CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
CTF      CALL ENLMD(VEC1,VEC2,LU3,LU4)
         call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,T_BUFF,1)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINE(LU4,-1)
         CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINE(LU2,-1)
         DO 8 JVEC = 1, IVEC
           CALL REWINE(LU3,-1)
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
        write(LUWRT,'(/A,I6/A)') ' Info from iteration',ITER,
     &                             '___________________________'
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
         CALL REWINE(LU1,-1)
         CALL REWINE(LU2,-1)
         EIGAPR = EIG(ITER-1,IROOT)
         DO 60 IVEC = 1, NVEC
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           IF(IVEC.EQ.1) THEN
             CALL REWINE( LU3 ,-1)
*                 SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
             CALL SCLVCD(LU2,LU3,FACTOR,VEC1,0,LBLK)
           ELSE
             CALL REWINE(LU3,-1)
             CALL REWINE(LU4,-1)
C                 VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU2,LU3,0,LBLK)
           END IF
C
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINE(LU3,-1)
           CALL REWINE(LU4,-1)
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
*. 1.2 : Multiply with inverse Hessian approximation to get new
*        direction
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU3
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(LUWRT,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINE(LUDIA,-1)
           CALL REWINE(LU3,-1)
           CALL REWINE(LU4,-1)
*. Assuming diagonal preconditioner
           IPRECOND = 1
           CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
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
             CALL REWINE(LU1,-1)
             DO 66 IVEC = 1, NVEC
               FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
               IF(IVEC.EQ.1) THEN
                 IF(NVEC.EQ.1) THEN
                   CALL REWINE( LU5 ,-1)
                   CALL SCLVCD(LU1,LU5,FACTOR,VEC1,0,LBLK)
                 ELSE
                   CALL REWINE( LU4 ,-1)
                   CALL SCLVCD(LU1,LU4,FACTOR,VEC1,0,LBLK)
                 END IF
               ELSE
                 CALL REWINE(LU5,-1)
                 CALL REWINE(LU4,-1)
                 CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                 CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
               END IF
   66        CONTINUE
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(LUWRT,*) '  (current  X ) '
               CALL WRTVCD(VEC1,LU5,1,LBLK)
             END IF
* (H0 - E )-1  * X on LU4
             CALL REWINE(LU5,-1)
             CALL REWINE(LU4,-1)
             CALL REWINE(LUDIA,-1)
*
             CALL H0M1TD(LU4,LUDIA,LU5,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
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
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,NVEC,VEC1,1,LBLK)
      DO 150 IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
CTF     CALL ENLMD(VEC1,VEC2,LU3,LU4)
        call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,T_BUFF,1)
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
        CALL REWINE( LU5,LBLK)
        CALL REWINE( LU1,LBLK)
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
 1170    FORMAT(/'  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(LUWRT,1180) ITER
 1180    FORMAT(/,'  Convergence was obtained in ',I3,' iterations')
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
          IF (RTCNV(IROOT)) THEN
             WRITE(6,1190) IROOT,FINEIG(IROOT)+EIGSHF
          ELSE
             WRITE(6,1191) IROOT,FINEIG(IROOT)+EIGSHF
          END IF
 1190 FORMAT(' The final eigenvalue',I5,F22.10,' (converged)')
 1191 FORMAT(' The final eigenvalue',I5,F22.10,' (NOT converged)')
          IF(IPRT.GE.400) THEN
            WRITE(LUWRT,1200)
 1200       FORMAT(' The final approximation to eigenvector ')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(LUWRT,1300)
 1300     FORMAT(/,'  Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(LUWRT,1310)
 1310     FORMAT
     &    (/,'  Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(LUWRT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        CALL REWINE(LU1,LBLK)
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
        DO IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
          IF (RTCNV(IROOT)) THEN
             WRITE(6,1190) IROOT,FINEIG(IROOT)
          ELSE
             WRITE(6,1191) IROOT,FINEIG(IROOT)
          END IF
        END DO
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
      END IF
*
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
C
C     Delete scratch units
      close(unit=LU3,status='DELETE')
      close(unit=LU4,status='DELETE')
      close(unit=LU5,status='DELETE')
C
      call qexit('MICDV4')
C
 1030 FORMAT(/,3X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(/,I6,7F15.8,/,(1H ,5X,7F15.8))
 9300 FORMAT(' >>>  CPU (WALL) TIME IN ITERATION: ',A,'(',A,')')
      END
