C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
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
C...   This means also that persons (e.g. post-docs) leaving the research
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
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE CIEIG5_REL(MV7_REL,INICI,EROOT,VEC1,
     &           VEC2,MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           NDIM,NBLK,NROOT,MAXVEC,MXCIIT,MXROOT,LUINCI,
     &           IPRT,PEIGVC,NPRDET,IPNTR,NP1,NP2,NQ,
     &           SCR1,ISCR1,SCR2,T_CC,
     &           EIGSHF,ICISTR,LBLK,IDIAG,THRES_E,
     &           NBATCH,LBATCHB,LBATCHE,
     &           IBLOCK,INIDEG,E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,
     &           NCLS,
     &           MXLNG,IRC,IDCOMH,IHAM12,SIGDEN_ROUTE)
*
* Master routine for CI diagonalization
*
* Modified to handle PQ - preconditioner , May 1990
* PICO,MICDV4 added spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "wrkspc.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION INIDET(100)
      DIMENSION IPNTR(*),ISCR1(*)
      dimension PEIGVC(*),SCR1(*),SCR2(*)
*
*. Output from Subspace dagonalization
*.
      DIMENSION EROOT(NROOT)
*
      IF( IPRT.GT. 1 )  WRITE(6,'(/A)')
     &'          *** information from ci diagonalization  ***'
      WRITE(6,*)
*
      NTEST=0
*
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

          NFINDM = MIN(NDIM,4*NROOT)
          CALL FNDMND(LUDIA,LBLK,VEC1,NFINDM,NFINDA,ISCR1(1+2*NFINDM),
     &                SCR1(1+2*NFINDM),ISCR1,SCR1,IPRT)
C              FNDMND(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,SCR,ISCAT,SUBVAL
          CALL REWINO(LU1)
          IBASE = 1
          TEST = 1.0D-10
*
          DO 100 IROOT = 1, NROOT
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
            WRITE(6,*) ' IROOT NDEG ', IROOT,NDEG
*
            IF (INIDEG.EQ.0.AND.NDEG.GT.1) THEN
              WRITE(6,*) ' WARNING WARNING WARNING WARNING ! '
              WRITE(6,*) ' DEGENERATE INITIAL VECTORS FOR CI '
C             WRITE(6,*) ' I AM NOT ALLOWED TO TAKE THIS INTO '
C             WRITE(6,*) ' CONSIDERATION SINCE  INIDEG = 0 '
C             WRITE(6,*)
C             WRITE(6,*) ' I hope you know what you are doing '
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
            CALL WRSVCD_REL(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
     &                      LUDIA)
*. Imaginary part
            CALL WRSVCD_REL(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,0,NDIM,
     &                      LUDIA)
            IBASE = IBASE + NDEG
  100     CONTINUE
        ELSE
* =====================================
*. Initial approximations are in PEIGVC
* =====================================
CTF  Currently never used !!!
          CALL REWINO(LU1)
          DO 1984 IROOT = 1, NROOT
C                WRSVCD_REL(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM)
            CALL WRSVCD_REL(LU1,LBLK,VEC1,IPNTR,
     &           PEIGVC((IROOT-1)*NPRDET+1),NPRDET,NDIM,LUDIA)
 1984     CONTINUE
        END IF
      END IF
*
* ======================================================
* Initial CI vectors are already on file LU1, do nothing
* ======================================================
*
      IF( INICI .LT. 0 ) THEN
*. Vectors assumed already in LU1
        IF(IPRT.GT. 1 )
     &  WRITE(6,*)' Initial CI vector assumed in place '
      END IF
*
*                 ========================
* 2 :                  Diagonalization
*                 ========================
*
      CALL QENTER('CIEIG5')
* Inverse iteration modified Davidson with 2 vectors in core
      IF(IPRT .GE. 2 ) THEN
         WRITE(6,*)
         WRITE(6,'(A,I3)')
     &   '  Number of roots to be converged..  ',NROOT
         WRITE(6,*)
         WRITE(6,'(A,I3)')
     &   '  Largest allowed number of vectors..',MAXVEC
         WRITE(6,*)
         WRITE(6,'(A,I3)')
     &   '  Allowed number of CI iterations  ..',MXCIIT
      END IF
*
      MRNRM = 1
      IF(IRC.EQ.1) THEN
       IFAC = 1
      ELSE
       IFAC = 2
      END IF
      MEIG = MRNRM + (MXCIIT+1)*NROOT
      MAPROJ = MEIG  + (MXCIIT+1)*2*MAXVEC + NROOT
      IF(IRC.EQ.2) THEN
        MAPROJI = MAPROJ +  MAXVEC*(MAXVEC+1)/2
      ELSE
        MAPROJI = MAPROJ
      END IF
      MAVEC = MAPROJI + MAXVEC*(MAXVEC+1)/2
      MWORK = MAVEC + MAXVEC ** 2*IFAC*IFAC
      MLFRSC = MWORK + MAXVEC*(MAXVEC+1)*IFAC*IFAC
*
      IF( IPRT .GE. 5 ) THEN
         WRITE(6,*) ' MRNRM MEIG MAPROJ MAVEC MWORK MLFRSC '
         WRITE(6,'(6I4)')  MRNRM,MEIG,MAPROJ,MAVEC,MWORK,MLFRSC
      END IF
*
*     IF( MLFRSC-1 .GT. LLWRK) THEN
*        WRITE(6,'(A,2I7)' )
*    &   ' Not enough memory in CIEIG5 : needed and available ',
*    &   MLFRSC-1, LLWRK
*        Call Abend2( ' insufficient memory in cieig5 ' )
*     END IF
*
* Proper allocation of scratch arrays:
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'CIEIG5')
      call memmar(KRNRM,MEIG-1,'ADDL  ',2,'KRNRMS')
      call memmar(KEIG,MAPROJ-MEIG,'ADDL  ',2,'KEIGS ')
      call memmar(KAPROJ,MAPROJI-MAPROJ,'ADDL  ',2,'KAPROJ')
      call memmar(KAPROJI,MAVEC-MAPROJI,'ADDL  ',2,'KAPROI')
      call memmar(KAVEC,MWORK-MAVEC,'ADDL  ',2,'KAVECS')
      call memmar(KWORK,MLFRSC-MWORK,'ADDL  ',2,'KWORKS')
*
      IF(IRC.EQ.1) THEN
        write(6,*) ' Real algebra not implemented.'
	write(6,*) ' Use LUCITA or some other non-relativistic code.'
	stop ' Quitting in CIEIG5_REL.'
*
      ELSE
*. Complex numbers may be encountered
*
        CALL CMICDV(MV7_REL,VEC1,VEC2,LU1,LU2,WORK(KRNRM),WORK(KEIG),
     &              SCR2(1),MXCIIT,
     &              NDIM,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NROOT,
     &              T_CC,
     &              WORK(KAPROJI),WORK(KAPROJ),WORK(KAVEC),WORK(KWORK),
     &              IPRT,NPRDIM,IPNTR,NP1,NP2,NQ,LBLK,EIGSHF,
     &              0,0,IDCOMH,IHAM12,SIGDEN_ROUTE)
*
C      SUBROUTINE CMICDV(MV7_REL,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
C     &                  NVAR,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NINVEC,
C     &                  T_CC,
C     &                  APROJI,APROJR,AVEC,WORK,IPRT,
C     &                  NPRDIM,IPNTR,NP1,NP2,NQ,LBLK,EIGSHF,
C     &                  IPREDZN,INONREL,IDCOMH,IHAM12,SIGDEN_ROUTE)
C
        call memmar(KDUM ,IDUM,'FLUSM ',2,'CIEIG5')
      END IF
*     ^ End if switch real/complex numbers
      ENOT = SCR1(KFIN)
*
      call qexit('CIEIG5')
      RETURN
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
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE CLASS_TRUNC_REL(NCLS,ICLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                       E_CONV,ICLS_A,N_TRN_CLS,E_TRUNC,W_TRUNC)
*
* Decide which classes of parameters that can be eliminated
*
* Jeppe Olsen, Jan 97
*              March '97 updated
*
*. Note in current version all energy contributions are
*. positive.
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION ICLS_L(NCLS)
      DIMENSION CLS_CT(NCLS),CLS_ET(NCLS),CLS_C(NCLS),CLS_E(NCLS)
*. Output
      INTEGER ICLS_A(NCLS)
*. Giving  the truncated classes
*  Additional output is
*              N_TRN_CLS : Number of truncated classes
*              E_TRUNC   : Estimated error of eliminating these classes
*                          (compared to expansion with only largest
*                           coefficients included).
*                          E_TRUNC is thus the error that should be
*                          added to the estmated error arising from
*                          eliminating small terms
*
      WRITE(6,*)
      WRITE(6,*) ' Welcome to TRUNC_CLASS '
      WRITE(6,*) ' ======================='
CM    WRITE(6,*)
CM    WRITE(6,*) ' Information about classes '
CM    WRITE(6,*) ' =========================='
CM    WRITE(6,*)
CM    WRITE(6,*)
CM   & ' Class    Dimension      E          E(Trunc)      C',
CM   & '         C(Trunc)'
CM    WRITE(6,*)
CM   & ' =============================================================',
CM   & '========='
CM     DO JCLS = 1, NCLS
CM       IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
CM   &   THEN
CM         WRITE(6,'(3X,I4,1X,I10,3X,4(E12.6,1X))')
CM   &     JCLS,ICLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
CM   &     CLS_CT(JCLS),CLS_C(JCLS)
CM       END IF
CM     END DO
*
       WRITE(6,*)
     & ' Required threshold for convergence of energy',E_CONV
*
*. First find total energy correction and
*. largest class contribution
*
      E_TOT = 0.0D0
      E_CLS_MAX= 0.0D0
      DO JCLS = 1, NCLS
       E_TOT = E_TOT + CLS_ET(JCLS)
       E_CLS_MAX = MAX( E_CLS_MAX,CLS_ET(JCLS))
      END DO
      WRITE(6,*) ' E_CLS_MAX = ', E_CLS_MAX
*
*. The truncation is done in two steps, first an overall
*  threshold for deleting classes is constructed,
*  and all classes with energy contributions less than this
*  threshold are excluded.
*. Then one takes a pass through remaining classes and
*  eliminates untill the quota is filled
*
*. Loop over thesholds, increase by a given amount (1.5) in each it
*. Until total energy threshold is passed
*
* Another strategy : start with E_CLS_MAX and decrease from there
      ILOOP = 0
      STEP = 1.2
 1000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
C         E_THRES = E_CONV/NCLS
          E_THRES =  E_CLS_MAX
        ELSE
C         E_THRES = E_THRES*STEP
          E_THRES = E_THRES/STEP
        END IF
        E_ELI = 0.0D0
        N_ELI = 0
        DO JCLS = 1, NCLS
          IF(ABS(CLS_ET(JCLS)).LE.E_THRES) THEN
            E_ELI = E_ELI + CLS_ET(JCLS)
            N_ELI = N_ELI + 1
          END IF
        END DO
        IF(ILOOP.EQ.10000) THEN
          WRITE(6,*) ' Loop count exceeded 10000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          Call Abend2( 'CLASS_TRUNC: Iloop.eq.10000' )
        END If
      IF(E_ELI.GT. E_CONV) GOTO 1000
*. We stepped out when the energy of the eliminated classes
*  became to large, go back to previous threshold
COLD  E_THRES = E_THRES/STEP
      WRITE(6,*)
     & ' Overall Threshold for eliminating classes',E_THRES
*
      N_PAS_CLS = 0
      L_PAS_CLS = 0
      N_TRN_CLS = 0
      L_TRN_CLS = 0
      N_ACT_CLS = 0
      L_ACT_CLS = 0
*
      E_TRUNC = 0.0D0
      E_TRUNCT = 0.0D0
      W_TRUNCT= 0.0D0
      W_TRUNC = 0.0D0
*
*. Eliminate classes with energy contribution less than E_THRES
*
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).NE.0.0D0.AND.CLS_ET(JCLS).LE. E_THRES) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        END IF
      END DO
      WRITE(6,*)
     &  ' Number of classes with contributions less than E_THRES',
     &   N_TRN_CLS
      WRITE(6,*) ' Energy contributions from these classes ',
     &   E_TRUNCT
* Eliminate remaining classes until thres hold is obtained
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).EQ.0.0D0) THEN
*. Passive class, no contribution before truncation
          N_PAS_CLS = N_PAS_CLS + 1
          L_PAS_CLS = L_PAS_CLS + ICLS_L(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1
     &  .AND.ABS(E_TRUNCT+CLS_ET(JCLS)).LT.E_CONV) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1) THEN
*. Class is active
          N_ACT_CLS = N_ACT_CLS + 1
          L_ACT_CLS = L_ACT_CLS + ICLS_L(JCLS)
           ICLS_A(JCLS) = 1
        END IF
      END DO
*. Correct for missing minus in first order correction
      E_TRUNC  = -E_TRUNC
      E_TRUNCT = -E_TRUNCT
*
C?    WRITE(6,*) ' ICLS_A giving active classes '
C?    WRITE(6,*) ' ============================='
C?    WRITE(6,*)
C?    CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
      WRITE(6,*)
      WRITE(6,*) ' Estimated complete energy contribution of ',
     &            '  eliminated classes ', E_TRUNCT
      WRITE(6,*) ' Estimated truncated energy contribution of ',
     &            '  eliminated classes ', E_TRUNC
      WRITE(6,*) ' Truncation error in eliminated classes was ',
     &              E_TRUNC-E_TRUNCT
      WRITE(6,*) ' Estimated truncated weight contribution of ',
     &            '  eliminated classes ', W_TRUNC
      WRITE(6,*) ' Estimated energy contribution without trunc',
     &             -E_TOT
      WRITE(6,*) ' Energy contribution of of active classes ',
     &             -E_TOT-E_TRUNCT

      WRITE(6,*)
      WRITE(6,*) '  Class      Number     Dimension  '
      WRITE(6,*) ' ================================= '
      WRITE(6,'(A,5X,I5,5X,I10)') ' Passive  ', N_PAS_CLS,L_PAS_CLS
      WRITE(6,'(A,5X,I5,5X,I10)') ' Truncated', N_TRN_CLS,L_TRN_CLS
      WRITE(6,'(A,5X,I5,5X,I10)') ' Active   ', N_ACT_CLS,L_ACT_CLS
*
      WRITE(6,*)
      WRITE(6,*) ' Information about classes '
      WRITE(6,*) ' =========================='
      WRITE(6,*)
      WRITE(6,*)
     & ' Class    Dimension      E          E(Trunc)      C',
     & '         C(Trunc)   Active'
      WRITE(6,*)
     & ' =============================================================',
     & '================'
       DO JCLS = 1, NCLS
         IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
     &   THEN
           WRITE(6,'(3X,I4,1X,I10,3X,4(E12.6,1X),1X,I2)')
     &     JCLS,ICLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
     &     CLS_CT(JCLS),CLS_C(JCLS),ICLS_A(JCLS)
         END IF
       END DO
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE CMICDV(mv7_rel,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LU6,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  T_CC,
     &                  APROJI,APROJR,AVEC,WORK,IPRT,
     &                  NPRDIM,IPNTR,NP1,NP2,NQ,LBLK,EIGSHF,
     &                  IPREDZN,INONREL,IDCOMH,IHAM12,SIGDEN_ROUTE)
*
* Multi root Davidson routine for complex Hermitian
* matrices , requires two blocks in core
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
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
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
      REAL * 8   INPROD,INPRDD
      character*6 SIGDEN_ROUTE
      LOGICAL CONVER,RTCNV(100),EXSTOP
*
#include "krmcluci_inf.h"
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION RNRM(NROOT,MAXIT),EIG(2*MAXVEC,MAXIT)
      DIMENSION APROJI(*),APROJR(*),AVEC(*),WORK(*)
      DIMENSION IPNTR(1)
*
* Dimensioning required of local vectors
*      APROJR : MAXVEC*(MAXVEC+1)/2
*      APROJI : MAXVEC*(MAXVEC+1)/2
*      AVEC   : (2*MAXVEC) ** 2
*      WORK   : 2*MAXVEC*(2*MAXVEC+1)/2
*
      DIMENSION FINEIG(*)
      EXTERNAL mv7_rel
*
      call qenter('CMICDV')
*
      NTEST = 0000
      IPRT = MAX(NTEST,IPRT)
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) '----------------------------------------------'
        write(6,*) ' I´ll take you to a place '
        write(6,*) ' Where we shall find our'
        write(6,*) '  ROOTS, BLOODY ROOTS '
        write(6,*)
        write(6,*) '               Sepultura '
        write(6,*) '----------------------------------------------'
        write(6,*)
      end if
*
      write(6,*)
      write(6,*) ' Nuclear repulsion + core energy: ',EIGSHF
*
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry CMICDV wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Increase MXCIV to be at least 2 * Nroot '
         Call Abend1(' Enforced stop in CMICDV ')
       END IF
*
       CONVLIM = 1.0D-5
       CONVER = .FALSE.
*
       IF(MAXIT.EQ.0) THEN
         WRITE(6,*) ' MAXIT = 0, so just a fast visit '
         RETURN
       END IF
*
       if (IPRT.ge.1) then
         WRITE(6,*) ' Diagonal in CMICDV'
         WRITE(6,*) ' =================='
         CALL REWINE(LUDIA,LBLK)
         CALL WRTVCDC(VEC1,LUDIA,1,LBLK)
       end if
*
       if (IDCOMH.eq.0) then
          I_DO_COMHAM = 0
       else if (IDCOMH.eq.1) then
          I_DO_COMHAM = 1
       else
          write(6,*) 'Problemo in CMICDV.'
          Call Abend2( 'Quitting.' )
       end if
*
       IF(I_DO_COMHAM .EQ.1) THEN
         write(6,*)
         write(6,*) ' ===================================== '
         WRITE(6,*) ' Control will be transferred to COMHAM '
         write(6,*) ' ===================================== '
         WRITE(6,*) ' NVAR = ', NVAR
         CALL COMHAM(NVAR,T_CC,SIGDEN_ROUTE)
         Call Abend2( ' After COMHAM  ' )
       END IF
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       NVEC = 0
       IF( MAXIT .GE. 1 ) THEN
*
         IF(IPRT.GE.1000) THEN
           WRITE(6,*) ' Initial vectors '
           WRITE(6,*) ' ================'
           CALL REWINE(LU1,LBLK)
           DO IVEC = 1, NINVEC
             WRITE(6,*) ' Initial vector ...  ', IVEC
             CALL WRTVCDC(VEC1,LU1,0,LBLK)
           END DO
         END IF
*
         CALL REWINE(LU1,LBLK)
         CALL REWINE(LU2,LBLK)
         DO 10 IVEC = 1,NINVEC
           CALL REWINE(LU3,LBLK)
           CALL REWINE(LU4,LBLK)
*. Copy to scratch unit LU3
           CALL COPVCDC(LU1,LU3,VEC1,0,LBLK)
*
*--------------------------------------------------------------
* Generalization to SIGDEN routine block implemented.
*--------------------------------------------------------------
           if (SIGDEN_ROUTE.eq.'TRADIT') then
             call mv7_rel(VEC1,VEC2,LU3,LU4,INONREL)
           else if (SIGDEN_ROUTE.eq.'SIGDEN') then
             call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,1)
           else
             write(6,*) 'Route for sigma and densities not allowed.'
             write(6,*) 'SIGDEN_ROUTE is ',SIGDEN_ROUTE
             call abend ('Quitting.')
           end if
*
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
           CALL REWINE(LU4,LBLK)
           CALL COPVCDC(LU4,LU2,VEC1,0,LBLK)
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
           WRITE(6,*) ' Initial sigma vectors '
           WRITE(6,*) ' ======================'
           CALL REWINE(LU2,LBLK)
           DO IVEC = 1, NINVEC
             WRITE(6,*) ' Initial vector ...  ', IVEC
             CALL WRTVCDC(VEC1,LU2,0,LBLK)
           END DO
         END IF
*
         IF( IPRT .GE.10.OR.IPREDZN.NE.0 ) THEN
           WRITE(6,*) ' Initial projected matrix, real part'
           CALL PRSYM(APROJR,NINVEC)
           WRITE(6,*) ' Initial projected matrix, imag part'
           CALL PRSYM(APROJI,NINVEC)
          END IF
*. Diagonalize initial projected matrix
C             CEIGEN(HR,HI,NDIM,EIG,C,SCR
C?       write(6,*) 'NINVEC FOR CEIGEN ',NINVEC
         CALL CEIGEN(APROJR,APROJI,NINVEC,EIG(1,1),AVEC,WORK)
*
         IF(IPRT .GE. 3.OR.IPREDZN.NE.0 ) THEN
           WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
           WRITE(6,'(F22.13)')
     &     ( EIG(IROOT,1)+EIGSHF,IROOT=1,NROOT)
         END IF
         NVEC = NINVEC
       END IF
*
       IF (MAXIT.EQ.1.OR.IPREDZN.EQ.1) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
       write(6,*)
       WRITE(6,*) ' Info from iteration .... ', ITER
       WRITE(6,*) ' _________________________'
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
           WRITE(6,*) '  ( HX - EX )'
           CALL WRTVCDC(VEC1,LU4,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         CALL REWINE(LU4,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,0,LBLK)
     &         +       INPRDD(VEC1,VEC1,LU4,LU4,0,LBLK) )
         RNRM(IROOT,ITER) = RNORM
         WRITE(6,'(A19,7X,I3,3X,1E18.13,3X,1F19.12)')
     &     ' Iter RNORM EIGAPR ', ITER,RNORM,EIGAPR+EIGSHF
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
             WRITE(6,*) '(D-E)-1 *( HX - EX ) '
             CALL WRTVCDC(VEC1,LU3,1,LBLK)
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
             WRITE(6,*)
     &       ' CMICDV : Vector eliminated due to linear dependence '
             IADD = IADD -1
             GOTO 100
           END IF
           SCALE = 1.0D0/SQRT(XNORM)
C          write(6,*) ' Scale, scale**2*xnorm',scale,
C    &                scale**2*xnorm
           CALL SCALVE(WORK,SCALE,2*(NVEC+IADD-1) )
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
*. Imag part of othogonal vector
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
             WRITE(6,*) '   Orthonormalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCDC(VEC1,LU5,1,LBLK)
           END IF
C               SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
           CALL SKPVCD(LU1,2*(NVEC+IADD-1),VEC1,1,LBLK)
           CALL REWINE(LU5,LBLK)
           CALL COPVCDC(LU5,LU1,VEC1,0,LBLK)
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
        CALL COPVCDC(LU1,LU3,VEC1,0,LBLK)
        CALL REWINE(LU4,LBLK)
        IF(IPRT.GE.10000) THEN
          WRITE(6,*) ' Input vector to mV7 '
          CALL WRTVCDC(VEC1,LU3,1,LBLK)
        END IF
*
* Generalization to SIGDEN routine block implemented.
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          call mv7_rel(VEC1,VEC2,LU3,LU4,INONREL)
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
          call sigden_ctrl(VEC1,VEC2,LU3,LU4,T_CC,1)
        end if
*
        IF(IPRT.GE.10000) THEN
          WRITE(6,*) ' Output vector from mv7_rel '
          CALL WRTVCDC(VEC1,LU4,1,LBLK)
        END IF
        CALL REWINE(LU4,LBLK)
        CALL COPVCDC(LU4,LU2,VEC1,0,LBLK)
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
        WRITE(6,*) ' Projected matrices before CEIGEN '
        CALL PRSYM(APROJR,NVEC)
        CALL PRSYM(APROJI,NVEC)
      END IF
      write(6,*) 'Dimension of projected matrix:',NVEC
      CALL CEIGEN(APROJR,APROJI,NVEC,EIG(1,ITER),AVEC,WORK)
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F22.13)')
     &   ( EIG(IROOT,ITER)+EIGSHF,IROOT=1,NROOT)
CJO      call towave(nroot,eig(1,iter))
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
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
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER.
     &   OR.IPREDZN.EQ.1)THEN
*. Reset C vectors
        CALL REWINE( LU5,LBLK)
        DO IROOT = 1, NROOT
*. Real part of vector IROOT on LU3
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+JVEC)
           WORK((JVEC-1)*2+2) =-AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
          END DO
          CALL MVCSMD(LU1,WORK,
     &    LU3,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
*. Imaginary part of vector on LU6
          DO JVEC = 1, NVEC
           WORK((JVEC-1)*2+1) = AVEC((IROOT-1)*2*NVEC+NVEC+JVEC)
           WORK((JVEC-1)*2+2) = AVEC((IROOT-1)*2*NVEC+JVEC)
          END DO
          CALL MVCSMD(LU1,WORK,
     &    LU6,LU4,VEC1,VEC2,2*NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
     &          + INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
C          write(6,*) ' Scale, scale**2*xnorm',scale,
C    &                scale**2*xnorm
          WORK(2*NROOT*NVEC+IROOT) = SCALE
*. Copy current vector to LU5
C?        WRITE(6,*) ' Current vector on LU3, IROOT = ', IROOT
C?        CALL WRTVCD(VEC1,LU3,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
*
C?        WRITE(6,*) ' Current vector on LU6, IROOT = ', IROOT
C?        CALL WRTVCD(VEC1,LU6,1,LBLK)
*
          CALL REWINE(LU6,LBLK)
          CALL SCLVCD(LU6,LU5,SCALE,VEC1,0,LBLK)
CE        IF(IPRT.GT.100) THEN
CE           WRITE(6,*) ' Normalized vector '
CE           call rewine(LU5,LBLK)
CE           CALL WRTVCDC(VEC1,LU5,1,LBLK)
CE        END IF
        END DO
*. Transfer C vectors to LU1
*
        if (IPRT.ge.1000) then
          WRITE(6,*) ' Initial vectors on LU1 '
          CALL REWINE(LU1,LBLK)
          DO IVEC = 1, NROOT
            WRITE(6,*) ' Vector ', IVEC
            CALL WRTVCDC(VEC1,LU1,0,LBLK)
            WRITE(6,*) ' Writing of vector completed '
          END DO
*
          WRITE(6,*) ' Vectors on LU5 '
          CALL REWINE(LU5,LBLK)
          DO IVEC = 1, NROOT
            WRITE(6,*) ' Vector ', IVEC
            CALL WRTVCDC(VEC1,LU5,0,LBLK)
            WRITE(6,*) ' Writing of vector completed '
          END DO
        end if
*
        CALL REWINE(LU1,LBLK)
        CALL REWINE(LU5,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCDC(LU5,LU1,VEC1,0,LBLK)
C?        WRITE(6,*) ' Vector copied : ', IVEC
        END DO
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' New reset C-vectors '
          WRITE(6,*) ' =================== '
          CALL REWINO(LU1)
          DO IVEC = 1, NROOT
             CALL WRTVCDC(VEC1,LU1,0,LBLK)
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
          CALL COPVCDC(LU5,LU2,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000) THEN
          WRITE(6,*) ' New reset sigma -vectors '
          WRITE(6,*) ' =================== '
          CALL REWINO(LU2)
          DO IVEC = 1, NROOT
             CALL WRTVCDC(VEC1,LU2,0,LBLK)
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
cbs
      EXSTOP =.FALSE. 
      INQUIRE(FILE='STOPPER',EXIST=EXSTOP)
      if (EXSTOP) then
      write(6,*) 'I will stop now, like requested '
      goto 1001
      endif
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER.AND.IPREDZN.EQ.0) THEN
        GOTO 1000
      ELSE IF(IPREDZN.EQ.1) THEN
        WRITE(6,*) ' PREDIA vectors written to file LUCVECT'
        RETURN
      END IF
 1001 CONTINUE

* ( End of loop over iterations )
*
      IF(IPREDZN.EQ.0) THEN
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 0 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT(/'   Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
C!?      ITER = ITER - 1
         IF (IPRT .GE. 0 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(/'  Convergence was obtained in ',I3,' iterations')
      END IF
      END IF
*
      if (IPRT.gt.0) then
        CALL REWINE(LU1,LBLK)
        CALL REWINE(LU2,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,'(/A,I3/A/)')
     &    ' Information about convergence for root... ' ,IROOT,
     &    ' ============================================'
          FINEIG(IROOT) = EIG(IROOT,ITER)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F22.10)
          IF(IPRT.GE.10) THEN
            WRITE(6,1201)
 1201       FORMAT(/' The final approximation to eigenvector')
            CALL WRTVCDC(VEC1,LU1,0,LBLK)
 1200       FORMAT(/' The final sigma vector')
            CALL WRTVCDC(VEC1,LU2,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(/'  Summary of iterations '
     +           /'  ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (/'  Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(IROOT,I)+EIGSHF,RNRM(IROOT,I)
 1340     FORMAT(7X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      else
        write(6,*)
        write(6,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(6,*)
        do 1601 IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(IROOT,ITER)+EIGSHF
          write(6,9000) ' Final Eigenvalue ',IROOT,' : ',FINEIG(IROOT)
 1601   continue
        write(6,*)
        write(6,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(6,*)
      end if
*
      ECI = EIG(1,ITER)
*
 9000 format(A18,1I3,A3,F20.13)
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(IROOT,ITER),
     &                 RNRM(IROOT,1)/RNRM(IROOT,ITER)
 1607   CONTINUE
      END IF
*
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
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      SUBROUTINE COMHAM(NVAR,T_CC,SIGDEN_ROUTE)
*
* Construct and diagonalize complete Hamiltonian
*
*
* LUDIA needed for giving form of vectors
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "clunit.inc"
*
      character*6 SIGDEN_ROUTE
*
      common/CECORER/ECORE
*
      PARAMETER(MXDIM = 700)
      DIMENSION VEC1(MXDIM),VEC2(MXDIM)
      DIMENSION SCR1((2*MXDIM)**2),SCR2((2*MXDIM)**2)
*
      NTEST = 000
*
      write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(6,*) 'MXDIM in COMHAM set to ',MXDIM
      write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
*
      LUIN = LUSC92
      LUOUT = LUSC93
      LBLK = -1
*
      LENH = NVAR**2
      call memmar(KHR,LENH,'ADDS  ',2,'HREAL  ')
      call memmar(KHI,LENH,'ADDS  ',2,'HIMAG  ')
*
*. Loop over elements in real part
      DO I = 1, NVAR
        WRITE(6,*) ' COMHAM : I = ', I
        CALL FLSHFO(6)
* Real part
        CALL REWINE(LUDIA,LBLK)
        CALL REWINE(LUIN,LBLK)
        ONE = 1.0D0
        CALL WRSVCD_REL(LUIN,LBLK,VEC1,I,ONE,1,NVAR,LUDIA)
* Imaginary part
        CALL WRSVCD_REL(LUIN,LBLK,VEC1,I,ONE,0,NVAR,LUDIA)
* Do the sigma
        CALL REWINE(LUOUT,LBLK)
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          call mv7_rel(VEC1,VEC2,LUIN,LUOUT,0)
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
          call sigden_ctrl(VEC1,VEC2,LUIN,LUOUT,T_CC,1)
        else
          write(6,*) 'Route for sigma and densities not allowed.'
          write(6,*) 'SIGDEN_ROUTE is ',SIGDEN_ROUTE
          call abend ('Quitting.')
        end if
*.
        CALL REWINE(LUOUT,LBLK)
        CALL MERGE_VEC(LUOUT,WORK(KHR+(I-1)*NVAR))
        CALL MERGE_VEC(LUOUT,WORK(KHI+(I-1)*NVAR))
C            MERGE_VEC(LU,VEC)
C       WRITE(6,*) I,' th Sigma vector read in in COMHAM '
C       CALL WRTMAT(WORK(KHR+(I-1)*NVAR),1,NVAR,1,NVAR)
      END DO
*
      if (NTEST.ge.50) then
         WRITE(6,*) ' Real and Imaginary parts of Hamiltonian '
         CALL WRTMAT(WORK(KHR),NVAR,NVAR,NVAR,NVAR)
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)
         CALL WRTMAT(WORK(KHI),NVAR,NVAR,NVAR,NVAR)
      end if
*
*. Check imag part > Antihermitian
      XMXDIFF = 1.0D-11
      DO I = 1, NVAR
        DO J = 1, I
          IF(ABS(WORK(KHI+(I-1)*NVAR+J-1)+WORK(KHI+(J-1)*NVAR+I-1)).
     &           GT.XMXDIFF) THEN
            Write(6,*)
            WRITE(6,*)
     &      ' Imag matrix not antihermitean, I,J, H(I,J),H(J,I), DIFF'
            WRITE(6,*) I,J,WORK(KHI+(I-1)*NVAR+J-1),
     &                     WORK(KHI+(J-1)*NVAR+I-1),
     &                     WORK(KHI+(I-1)*NVAR+J-1)
     &                     +WORK(KHI+(J-1)*NVAR+I-1)
          END IF
        END DO
      END DO
*. Check real part : Hermitian ?
      XMXDIFF = 1.0D-11
      DO I = 1, NVAR
        DO J = 1, I
          IF(ABS(WORK(KHR+(I-1)*NVAR+J-1)-WORK(KHR+(J-1)*NVAR+I-1))
     &                .GT.XMXDIFF) THEN
            Write(6,*)
            WRITE(6,*)
     &      ' Real matrix not hermitean, I,J, H(I,J),H(J,I), DIFF'
            WRITE(6,*) I,J,WORK(KHR+(I-1)*NVAR+J-1),
     &                     WORK(KHR+(J-1)*NVAR+I-1),
     &                     WORK(KHR+(I-1)*NVAR+J-1)
     &                     -WORK(KHR+(J-1)*NVAR+I-1)
          END IF
        END DO
      END DO
*
*. Diagonalize
*
*. But first reform to lower half form
      SIGNTP = 1.0
      CALL TRIPAK(WORK(KHR),SCR1,1,NVAR,NVAR,SIGNTP)
      CALL COPVEC(SCR1,WORK(KHR),NVAR*(NVAR+1)/2)
*. Read in two-electron terms from unit 30
*. Kept for historical (sentimental) reasons
C     IDIR2E = 1
C     IF(IDIR2E.EQ.1) THEN
C       LU30 =
C       REWIND (LU30)
C       DO IDIM = 1, NVAR
C         READ(LU30,*) (SCR1(IDIM*(IDIM-1)/2+JCOL),JCOL=1,IDIM)
C       END DO
C       WRITE(6,*) ' 2 el  terms read in '
C       CALL PRSYM(SCR1,NVAR)
*. And add
C       call VECSUM(HR,HR,SCR1,1.0D0,1.0D0,NVAR*(NVAR+1)/2)
C       write(6,*) ' Updated HR '
C       CALL PRSYM(HR,NVAR)
C     END IF
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
      CALL TRIPAK(WORK(KHI),SCR1,1,NVAR,NVAR,SIGNTP)
      CALL COPVEC(SCR1,WORK(KHI),NVAR*(NVAR+1)/2)
*
      write(6,*) 'NVAR : ',NVAR
      NVAR_SAVE = NVAR
      CALL CEIGEN(WORK(KHR),WORK(KHI),NVAR,VEC1,SCR1,SCR2)
C          CEIGEN(HR,HI,NDIM,EIGVAL,EIGVEC,SCR)
      write(6,*) 'NVAR_SAVE : ',NVAR_SAVE
*
*  Add nuclear repulsion energy
      DO I = 1, NVAR_SAVE
         VEC1(I) = VEC1(I) + ECORE
      END DO
*
      WRITE(6,*) ' Lowest Eigenvalue ', VEC1(1)
      WRITE(6,*) ' All eigenvalues '
      DO I = 1, NVAR_SAVE
         WRITE(6,'(X,I4,2X,1F18.10)') I,VEC1(I)
      END DO
*
*  Print MXANAVEC eigenvectors
      if (IPRWFC.ge.1) then
        do IANA=1,MXANAVEC,1
          write(6,*)
          write(6,'(A,X,I4)') ' Eigenvector number ... ',IANA
        end do
      end if
*
      close(unit=LUIN,status='DELETE')
      close(unit=LUOUT,status='DELETE')
*
      RETURN
      END
