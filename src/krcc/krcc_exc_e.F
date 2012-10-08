      SUBROUTINE KRCC_EXC_E(CCAMP,CCVEC1,CCVEC2,
     &                        CCDIA,ECC,
     &                       WORK,KFREE,LFREE)
*
* Master routine for coupled cluster linear response calculation 
* of excitation energies
*
* Output : CC-excitation operators are stored on LU_CCEXC_OP
*
* Based on outer setup from Jeppe Olsen, May 2000
*
* New algorithm from Lasse 2012
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "csm.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "clunit.inc"
#include "glbbas.inc"
!#include "cc_exc.inc" Check if there is something on this common block
!#include "lorr.inc"
#include "symm.inc"
*
*. Scratch vectors ( To show Jesper that I am on the right track )
      DIMENSION CCAMP(*),CCVEC1(*),WORK(*)
      DIMENSION CCVEC2(*),CCDIA(*)
*. Local scratch, assuming atmost 1000 roots per sym 
      INTEGER IREO(1000)
*
      IDUM = 0
*
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) ' **********************************************'
      WRITE(6,*) ' *                                            *'
      WRITE(6,*) ' * CCLR calculation of excitation energies    *'
      WRITE(6,*) ' *                                            *'
      WRITE(6,*) ' * Based on KRCC algorithm by Lasse Soerensen *'
      WRITE(6,*) ' *                                            *'
      WRITE(6,*) ' * Outer setup from Jeppe Olsens CI-driven CC *'
      WRITE(6,*) ' *                                            *'
      WRITE(6,*) ' **********************************************'
      DO i=1,100
      WRITE(6,*) ' WHAT IS MV8 '
      END DO
*
* Check if there is a vector to start from
*
      CALL INI_CCLR(II_RES_EXC)
*
*  Analyze the CC expansion in terms of leading amplitudes
* Why hardwire this?
      I_ANA_CC = 0
*
*  For controlling some dimensionings in sigden_ctrl
* will not be needed but check what it will do
      I_AM_IN_CC_EXC_E = 1
*
*. Copy initial sets of excitation operators from file 
      IF(II_RES_EXC.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' Reading excitation operators from file ',
     &             LU_CCEXC_OP
        CALL REWINO(LU_CCEXC_OP)
        CALL REWINO(LU_CCVECFL)
*
        DO ISM = 1,NSMST
          IF(NEXC_PER_SYM(ISM).NE.0) THEN 
C          IMSCOMB_CC = 0
           CALL IDIM_TCC_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TPE,ISM,
     &          MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &          WORK(KLLSOBEX_CC),WORK(KLIBSOBEX_CC),LEN_T_VEC,
     &          MX_SBSTR,
     &          IFTONE,IFTTWO,-1,N1ELINT,N2ELINT)
C               FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
           NROOT_CC = NEXC_PER_SYM(ISM)
           CALL IFRMDS(LEN_T_VEC_PREV,1,-1,LU_CCEXC_OP)
           DO IROOT_CC = 1, NROOT_CC
             ZERO = 0.0D0
             CALL SETVEC(CCVEC1,ZERO,LEN_T_VEC)
             CALL FRMDSC(CCVEC1,LEN_T_VEC_PREV,-1,LU_CCEXC_OP,
     &                   IMZERO,I_AM_PACKED)
             CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_CCVECFL)
           END DO
          END IF 
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Writing information to file ',LU_CCVECFL
      END IF
*     ^ End of this run is a restart 
*
      CALL REWINO(LU_CCEXC_OP)
      IF (II_RES_EXC.NE.0) CALL REWINO(LU_CCVECFL)


        print*
        print*,' SYMMETRY HANDLING NEEDS CLEANUP. '
        print*,' USING LUCIAREL ORDERING ON NONREL. INPUT.'
        print*,' FIXME !!'
        print*,' FIXME !!'
        print*,' FIXME !!'
        print*,' NSMST : ',NSMST
        print*,' NIRR_PN : ',NIRR_PN
        print*,' NEXC_PER_SYM : ',(NEXC_PER_SYM(J),J=1,NIRR_PN,1)
        print*,' NOT SURE IF I SHOULD USE NIRR_PN '
*
* Loop over either boson or fermion double group irreps
      DO ISM = 1, NSMST
        IF(NEXC_PER_SYM(ISM).NE.0) THEN 
          WRITE(6,*)
          WRITE(6,*) ' ============================================='
          WRITE(6,*) ' Information for excitations of symmetry',ISM
          WRITE(6,*) ' ============================================='
          WRITE(6,*)
          WRITE(6,*) ' Number of roots required : ', NEXC_PER_SYM(ISM)
          NROOT_CC = NEXC_PER_SYM(ISM)
          ITEX_SM = ISM
*
*. Information needed for diagonal of this operator symmetry
          CALL IDIM_TCC_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TPE,ISM,
     &                      MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &                      WORK(KLLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &                      LEN_T_VEC,MX_SBSTR,
     &                      IFTONE,IFTTWO,-1,N1ELINT,N2ELINT)
          N_CC_AMP = LEN_T_VEC
! No point in finding roots without amplitudes. Lasse
          IF(N_CC_AMP.EQ.0) CYCLE 
* Reset due to symmetry dependence
          WRITE(6,*) ' Number of CC amplitudes ',N_CC_AMP
*
*. Set up Diagonal for this symmetry
          CALL GENCC_F_DIAG_M_KRCC(WORK(KLSOBEX_CC),NSPOBEX_TPE,CCDIA,
     &                        ISM,KDUM,KDUM,MX_ST_TSOSO,MX_ST_TSOSO_BLK,
     &                        N1ELINT,N2ELINT,WORK(KT_CC),
     &                        WORK,KFREE,LFREE)
*
*. Save on LUDIA for future use
          CALL REWINO(LUDIA)
          CALL TODSC(CCDIA,LEN_T_VEC,-1,LUDIA)
*
*. Initialization : Restart or lowest diagonal elements 
          IF (II_RES_EXC.EQ.0) THEN
*. Lowest elements of diagonal
            CALL DUMSORT(CCDIA,LEN_T_VEC,NROOT_CC,IREO,CCVEC2)
*. Save corresponding initial guesses
            CALL REWINO(LU_CCVECF)
            DO IROOT = 1, NROOT_CC
              ZERO = 0.0D0
              CALL SETVEC(CCVEC1,ZERO,LEN_T_VEC)
              CCVEC1(IREO(IROOT)) = 1.0D0
              CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_CCVECF)
C?            WRITE(6,*) ' Initial vector for IROOT = ', IROOT
C?            CALL WRTMAT(CCVEC1,1,LEN_T_VEC,1,LEN_T_VEC) 
            END DO
          ELSE IF (II_RES_EXC.NE.0) THEN
            WRITE(6,*)
            WRITE(6,*) ' Reading initial CCLR vectors from file ',
     &                 LU_CCVECFL
            CALL REWINO(LU_CCVECF)
*. Read initial vectors from LU_CCVECFL
            DO IROOT =1, NROOT_CC
             CALL FRMDSC(CCVEC1,LEN_T_VEC,-1,LU_CCVECFL,
     &                   IMZERO,I_AM_PACKED)
             CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_CCVECF)
            END DO
          END IF
*         ^ End of shift between different forms of initialization
*. Allocate scratch space for the diagonalization 
*. Length of APROJ : MAXVEC**2
*. Length of AVEC  : MAXVEC**2
*  Length of WORK  : 4*MAXVEC**2 + 7*MAXVEC
*  Length of H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
          MAXVEC = MXCIV*NROOT_CC
          CALL MEMGET('REAL',KLAPROJ,MAXVEC**2,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KLAVEC ,MAXVEC**2,WORK,KFREE,LFREE)
          LWORK2 = 4*MAXVEC**2 + 7*MAXVEC
          CALL MEMGET('REAL',KLWORK,LWORK2,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KLSCR ,0,WORK,KFREE,LFREE)
          LEN = NROOT_CC*MXITLR
          CALL MEMGET('REAL',KLRNRM,LEN,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KLEIG ,LEN,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KLFINEIG,NROOT_CC,WORK,KFREE,LFREE)
          print*,'MXCIV,NROOT_CC,LWORK2,LEN,MXITLR',
     &            MXCIV,NROOT_CC,LWORK2,LEN,MXITLR
*
*. Performing Jacobian times right vectors for excitation energies
          L_OR_R = 2
          IPRDIA = 1
          CALL GENMAT_DIAG_REL(ECC,CCVEC1,CCVEC2,
     &         LU_CCVECF,LU_CCVECL,LU_CCEXC_OP,WORK(KLRNRM),
     &         WORK(KLEIG),WORK(KLFINEIG),MXITLR,LEN_T_VEC,LU_CCVEC,
     &         LUDIA,
     &         NROOT_CC,MAXVEC,NROOT_CC,WORK(KLAPROJ),WORK(KLAVEC),
     &         WORK(KLWORK),IPRDIA,0,0.0D0,0,0,0,0,0.0D0,0.0D0,0,0,
     &         CCCONV_EX)
*
* Analyze excitation vectors - and copy to LU_CCEXC
          if (I_ANA_CC.eq.1) then
            CALL REWINO(LU_CCVECF)
            CALL ITODS(LEN_T_VEC,1,-1,LU_CCEXC_OP)
*
            DO IROOT = 1, NROOT_CC
              WRITE(6,*)
              WRITE(6,'(2X ,A,I3)')
     &        ' * * * * * * * * * * * * * * * * * * * * * * * * * *'
              WRITE(6,'(2X ,A,I3)')
     &        ' Analysis of occupation for ROOT = ',IROOT
              WRITE(6,'(2X ,A,I3)')
     &        ' * * * * * * * * * * * * * * * * * * * * * * * * * *'
              WRITE(6,*)
*
              CALL FRMDSC(CCVEC1,LEN_T_VEC,-1,LU_CCVECF,IMZERO,IAMPACK)
              CALL ANA_GENCC_CC_KRCC(CCVEC1,ISM,
     &                               WORK,KFREE,LFREE)
              CALL TODSC(CCVEC1,LEN_T_VEC,-1,LU_CCEXC_OP)
            END DO
          end if
          CALL MEMCHK_KRCC(WORK)
          CALL MEMREL('SYM',WORK,KLAPROJ,KLAPROJ,KFREE,LFREE)
*
        END IF
*       ^ End if NEXC_PER_SYM(ISM).NE.0
        CALL MEMCHK_KRCC(WORK)
      END DO
      CALL MEMCHK_KRCC(WORK)
*     ^ End loop over symmetries
*
      RETURN
      END
*
      SUBROUTINE INI_CCLR(II_RES_EXC)
*
* Find out if there is a vector to start from
*
* Initialize Coupled Cluster amplitudes
*
* Made automated restart from LU_CCVECFL if it is there
* Also checks for formatted or unformatted file
*
* Lasse, 2012
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
*(    ^contains number of single and double excitations and N_CC_AMP)
#include "clunit.inc"
*. Amplitudes to be defined
      LOGICAL YN
      CHARACTER(LEN=2) UNFORM
*
* Automated restart if LU_CCVECFL and LU_CCEXC_OP are present. 
* If only one is present we quit.
*
      INQUIRE(FILE='LU_CCVECFL',EXIST=YN)
      IF(YN) THEN
        II_RES_EXC = 1
        INQUIRE(FILE='LU_CCVECFL',UNFORMATTED=UNFORM)
        IF(UNFORM(1:2).EQ.'UN') THEN
          I_FORMATTED = 0
        ELSE
          I_FORMATTED = 1
        END IF
      ELSE
        II_RES_EXC = 0
      END IF
*
      INQUIRE(FILE='LU_CCEXC_OP',EXIST=YN)
      IF(YN) THEN
        IF(II_RES_EXC.EQ.0) THEN 
          CALL QUIT(' Only LU_CCEXC_OP and not LU_CCVECFL present ') 
        END IF
        INQUIRE(FILE='LU_CCEXC_OP',UNFORMATTED=UNFORM)
        IF(UNFORM(1:2).EQ.'UN') THEN
          I_FORMATTED = 0
        ELSE
          I_FORMATTED = 1
        END IF
      ELSE
        IF(II_RES_EXC.EQ.1) THEN 
          CALL QUIT(' Only LU_CCVECFL and not LU_CCEXC_OP present ') 
        END IF
      END IF 
*
      IF(II_RES_EXC.EQ.1) THEN
        WRITE(6,*) ' A set of amplitudes have been found '
        WRITE(6,*) ' This is a restarted calculation '
      END IF
*
      RETURN
      END
*
      SUBROUTINE DUMSORT(VEC,NDIM,NELMNT,IREO,ISCR)      
*
* Extremely stupid routine for finding the ordering of 
* elements in VEC. Only the lowest NELMNT elements are obtained
*
* On output : IREO: ordered => unordered index
*
* Author : Prefers to be anonymous, May 2000
*
*
#include "implicit.inc"
*. Input
      DIMENSION VEC(NDIM)
*. Output 
      INTEGER IREO(NELMNT)                
*. Scratch through parameter list
      INTEGER ISCR(NDIM)
*
      IZERO = 0
      CALL ISETVC(ISCR,IZERO,NDIM)   
      XMAX = FNDMNX(VEC,NDIM,2)
C?    WRITE(6,*) ' XMAX = ', XMAX
*
      DO I = 1, NELMNT
*. Find the I'th lowest element 
        IELMNT = 0
        XVAL = XMAX
        DO J = 1, NDIM
C?       WRITE(6,*) ' I,J,XVAL,VEC(J) = ', I,J,XVAL,VEC(J)
         IF(VEC(J).LE.XVAL.AND.ISCR(J).EQ.0) THEN
C?         WRITE(6,*) ' Update of lowest element I and J = ', I,J
           XVAL = VEC(J)
           IELMNT = J
         END IF
        END DO
        IREO(I) = IELMNT
        ISCR(IELMNT) = I
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Elements to be sorted : '
        CALL WRTMAT(VEC,1,NDIM,1,NDIM)
      END  IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Reorder array : new => old order '
        CALL IWRTMA(IREO,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
*
      REAL*8 FUNCTION FNDMNX(VECTOR,NDIM,MINMAX)
C
C     FIND SMALLEST(MINMAX=1) OR LARGEST(MINMAX=2)
C     ABSOLUTE VALUE OF ELEMENTS IN VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(1)
C
      IF(MINMAX.EQ.1) THEN
       RESULT=ABS(VECTOR(1))
       DO 100 I=2,NDIM
        RESULT=MIN(RESULT,ABS(VECTOR(I)))
  100  CONTINUE
      END IF
C
      IF(MINMAX.EQ.2) THEN
       RESULT=ABS(VECTOR(1))
       DO 200 I=2,NDIM
        RESULT=MAX(RESULT,ABS(VECTOR(I)))
  200  CONTINUE
       END IF
C
       FNDMNX=RESULT
      RETURN
      END
*
      SUBROUTINE GENMAT_DIAG_REL(ECC,VEC1,VEC2,
     &                  LU1,LU2,LU_CCEXC_OP,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &                  IOLSEN,IPICO,CCCONV_EX)
*
* p.t. matrix vector routine is hardwired
*
* Iterative solution of eigenvalue problem for general matrix
*
* Final and intermediate eigenvalues and eigenvectors are assumed real
*
* MIN version requiring two vectors in core
*
*
* Jeppe Olsen May 2000, from MINDV4
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors, each must be dimensioned to hold
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
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       LOGICAL CONVER
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION FINEIG(1)
       DIMENSION H0(*),IPNTR(*)
*
*. Scratch through argument list
*
*. Length of APROJ : MAXVEC**2
*. Length of AVEC  : MAXVEC**2
*  Length of WORK  : 4*MAXVEC**2 + 7*MAXVEC
*  Length of H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0SCR(*)
*. Local scratch - atmost 1000 roots
       LOGICAL RTCNV(1000)
*. Ad hoc arrays for root hooming
       DIMENSION XJEP(3000),IXJEP(3000)
*
*
       iprt = 0
       IOLSTM = 0
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Variational'
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Perturbational'
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Normal Davidson, Variational '
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Normal Davidson, Perturbational'
C      IF( MAXVEC .LT. 2 * NROOT ) THEN
C        WRITE(6,*) ' SORRY MINDV2 WOUNDED , MAXVEC .LT. 2*NROOT '
C        STOP ' ENFORCED STOP IN MINDV2'
C      END IF
       WRITE(6,*) ' Convergence threshold for residual : ', CCCONV_EX
*
       I_DO_ROOTHOOMING = 0
       IF(I_DO_ROOTHOOMING.EQ.1) THEN
         WRITE(6,*) ' Root hooming active '
       END IF
*
       IF(IPICO.NE.0) THEN
         MAXVEC = 2*NROOT
       END IF
*. Division of scratch memory
      KAPROJ = 1
      KFREE = KAPROJ + MAXVEC**2
*
      KARVAL = KFREE
      KFREE = KARVAL + MAXVEC
*
      KAIVAL = KFREE
      KFREE = KAIVAL + MAXVEC
*
      KARVEC = KFREE
      KFREE = KARVEC + MAXVEC**2
*
      KAIVEC = KFREE
      KFREE = KAIVEC + MAXVEC**2
*
      KZ = KFREE
      KFREE = KZ    + MAXVEC**2
*
      KW = KFREE
      KFREE = KW    + MAXVEC
*
      KSCR1 = KFREE
      KFREE = KSCR1 + 4*MAXVEC
*
C     CCCONV_EX = 1.0D-10
      CONVER = .FALSE.
      DO 1234 MACRO = 1,1
*
*.   INITAL ITERATION
*
* ===========================================================
*. The initial vectors does not neccessarily constitute an
*. orthonormal basis. Start by orthogonalizing
* ===========================================================
*
       DO IROOT = 1, NINVEC
*. Read vector IROOT in
         CALL REWINO(LU1)
         DO ISKP = 1, IROOT
           CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
         END DO
*. Diagonal element
         WORK(KAPROJ-1+(IROOT-1)*NINVEC+IROOT) =
     &   INPROD(VEC1,VEC1,NVAR)
         DO JROOT = IROOT+1, NINVEC
           CALL FRMDSC(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
           WORK(KAPROJ-1+(IROOT-1)*NINVEC+JROOT) =
     &     INPROD(VEC1,VEC2,NVAR)
           WORK(KAPROJ-1+(JROOT-1)*NINVEC+IROOT) =
     &     WORK(KAPROJ-1+(IROOT-1)*NINVEC+JROOT)
         END DO
       END DO
*
       IF(IPRT.GE.1) THEN
        WRITE(6,*) ' Overlap matrix of initial basis '
        CALL WRTMAT(WORK(KAPROJ),NROOT,NROOT,NROOT,NROOT)
       END IF
*. Orthonormal basis for the NROOT lowest roots
       CALL MGS3(WORK(KARVEC),WORK(KAPROJ),NINVEC,WORK(KAIVEC))
*. Orthogonalize, save orthogonlized vectors on LU3
       CALL REWINO( LU3)
       DO IROOT = 1, NINVEC
         CALL REWINO( LU1)
         CALL SETVEC(VEC1,0.0D0,NVAR)
         DO IVEC = 1, NINVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR =  WORK(KARVEC-1+(IROOT-1)*NINVEC+IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
         END DO
         CALL TODSC(VEC1,NVAR,-1  ,LU3)
       END DO
*
       CALL REWINO( LU1)
       CALL REWINO( LU3)
       DO IVEC = 1,NINVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
         CALL TODSC (VEC1,NVAR,-1,  LU1)
       END DO
*
       ITER = 1
       CALL REWINO( LU1 )
       CALL REWINO( LU2 )
* Mat * initial  vectors
       DO JVEC = 1,NINVEC
         CALL FRMDSC(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
*
         CALL MV10_REL(VEC1,VEC2,WORK,KFREE,LFREE)
         CALL TODSC(VEC2,NVAR,-1  ,LU2)
       END DO
*. Projected matrix
       CALL GET_NONSYM_SUBMAT(LU1,LU2,NINVEC,0,APROJ,NVAR,
     &                        VEC1,VEC2,0.0D0)
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL WRTMAT(APROJ,NINVEC,NINVEC,NINVEC,NINVEC)
       END IF
*  Diagonalize initial subspace matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*NINVEC)
       CALL EIGGMT3(WORK(KAPROJ),NINVEC,WORK(KARVAL),WORK(KAIVAL),
     &              WORK(KARVEC),WORK(KAIVEC),WORK(KZ),WORK(KW),
     &              WORK(KSCR1),1)
C          EIGGMT3(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
C    &                   Z,W,SCR,IORD)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KARVAL-1+IROOT)
   20  CONTINUE
       CALL COPVEC(WORK(KARVEC),AVEC,NINVEC**2)
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
      IF(IPRT  .GE. 0 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Running CCLR iteration   ', ITER
      END IF
*
* Save a copy of excitation vectors to LU_CCEXC
c      call rewino(LU1)
c      call itods(NVAR,1,-1,LU_CCEXC_OP)
c      do IROOT = 1, NROOT_CC
c        call frmdsc(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
c        call todsc(VEC1,NVAR,-1,LU_CCEXC_OP)
c      end do
*
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
         CALL REWINO(LU1)
         DO 50 IVEC = 1, NVEC
           CALL FRMDSC(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR = (-EIGAPR)*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   50    CONTINUE
           IF ( IPRT  .GE.600 ) THEN
             WRITE(6,*) '  ( HX - EX ) '
             CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
           END IF
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
         RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
         RNRM(ITER-1,IROOT) = RNORM
C         write(6,*) ' Current norm : ',RNORM
         WRITE(6,'(A19,7X,I3,3X,1E18.13,3X,1F19.12)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
         IF(RNORM.LT. CCCONV_EX ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
*
         IF( ITER .GT. MAXIT) GOTO 100
*.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
           CALL REWINO( LUDIA)
           CALL FRMDSC(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL H0M1TV_REL(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
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
              CALL H0M1TV_REL(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
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
                FACTOR = (-DELTA)/GAMMA
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
        CALL MV10_REL(VEC1,VEC2,WORK,KFREE,LFREE)
        CALL TODSC(VEC2,NVAR,-1  ,LU2)
  150 CONTINUE
*.Augment projected matrix
C     GET_NONSYM_SUBMAT(LUC,LUHC,NVEC,NVECP,SUBMAT,NDIM,
C    &                             VEC1,VEC2,SUBMATP)
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC**2)
      CALL GET_NONSYM_SUBMAT(LU1,LU2,NVEC+IADD,NVEC,APROJ,NVAR,
     &                       VEC1,VEC2,WORK(KAPROJ))
*  DIAGONALIZE PROJECTED MATRIX
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*NVEC)
      CALL EIGGMT3(WORK(KAPROJ),NVEC,WORK(KARVAL),WORK(KAIVAL),
     &             WORK(KARVEC),WORK(KAIVEC),WORK(KZ),WORK(KW),
     &             WORK(KSCR1),1)

      DO  IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KARVAL-1+IROOT)
      END DO
      CALL COPVEC(WORK(KARVEC),AVEC,NROOT*NVEC)
      IF(I_DO_ROOTHOOMING.EQ.1) THEN
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
       IF(IPRT.GE.0) THEN
         WRITE(6,*) ' New roots choosen as vectors '
         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
       END IF
*. Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
       END DO
       CALL COPVEC(XJEP,AVEC,NVEC*NVEC)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW) = WORK(KARVAL-1+IOLD)
       END DO
       DO INEW = 1, NROOT
         EIG(ITER,INEW)=XJEP(INEW)
       END DO
*
       IF(IPRT.GE.3) THEN
         WRITE(6,*) ' Reordered AVEC arrays '
         CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
       END IF
*
      END IF
*     ^ End of root homing procedure
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
           WRITE(6,*) ' Norms of residuals '
           WRITE(6,'(10E13.5)') (RNRM(ITER-1,KROOT),KROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' PROJECTED MATRIX AND EIGEN PAIRS '
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
      END IF
*
**  PERHAPS RESET OR ASSEMBLE CONVERGED EIGENVECTORS
*
  901 CONTINUE
*
      IPULAY = 0
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
*. Select space spanning lowest NROOT as new subspace
*. Note that subspace is required to be orthogonal although
*. the eigenvectors in general not are orthogonal
*
*. Overlap matrix of lowest NROOT eigenvectors
        DO IROOT = 1, NROOT
          DO JROOT = 1, NROOT
            WORK(KAPROJ-1+(IROOT-1)*NROOT+JROOT) =
     &      INPROD(AVEC((IROOT-1)*NVEC+1),
     &             AVEC((JROOT-1)*NVEC+1),NVEC)
          END DO
        END DO
*
        IF(IPRT.GE.10) THEN
         WRITE(6,*) ' Overlap matrix in new basis '
         CALL WRTMAT(WORK(KAPROJ),NROOT,NROOT,NROOT,NROOT)
        END IF
*. Orthonormal basis for the NROOT lowest roots
C       MGS3(X,S,NDIM,SCR1)
        CALL MGS3(WORK(KARVEC),WORK(KAPROJ),NROOT,WORK(KAIVEC))
*. Orthogonalization matrix is now in WORK(KARVEC)
*  Obtain new basis  - if iteration procedure continues
        IF(ITER.LT.MAXIT .AND. (.NOT.CONVER)) THEN
          CALL MATML4(WORK(KAIVEC),AVEC,WORK(KARVEC),
     &    NVEC,NROOT,NVEC,NROOT,NROOT,NROOT,0)
          CALL COPVEC(WORK(KAIVEC),AVEC,NROOT*NVEC)
        END IF
*
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
        CALL GET_NONSYM_SUBMAT(LU1,LU2,NROOT,0,APROJ,NVAR,
     &                       VEC1,VEC2,WORK(KAPROJ))
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
      iprt = 2
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
*
      write(6,*)
      write(6,*) '**********************************'//
     &           '********************'
      write(6,*) '               Root       Energy  '//
     &           '       RESIDUAL     '
      write(6,*) '**********************************'//
     &           '********************'
      write(6,*)
      write(6,'(1X,A,F18.10)') ' ground state       ',ECC
      do IROOT = 1, NROOT
         FINEIG(IROOT) = FINEIG(IROOT) + ECC
         write(6,'(2X,A,I3,3X,F18.10,2X,E10.3)')
     &             'Final energy ',
     &             IROOT,FINEIG(IROOT),RNRM(ITER,IROOT)
      end do
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
*
      SUBROUTINE MGS3(X,S,NDIM,SCR1)
*
* Modified Gram-Schmidt procedure by forward orthogonalization
*
* Jeppe Olsen, Summer of 99
*
* S is input overlap matrix, X is output set of orthonormalized vectors
*
#include "implicit.inc"
      REAL*8 INPROD
*. input 
      DIMENSION S(NDIM,NDIM)
*. Output
      DIMENSION X(NDIM,NDIM)
*. Scratch : vector of length NDIM
      DIMENSION SCR1(*)
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
        CALL MATVCB(S,X(1,IVEC),SCR1,NDIM,NDIM,0)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
        XNORM = INPROD(X(1,IVEC),SCR1,NDIM)
        FACTOR = 1.0D0/SQRT(XNORM)
        CALL SCALVE(X(1,IVEC), FACTOR, NDIM)
        CALL SCALVE(SCR1,FACTOR,NDIM)
*. Subtract X(1,IVEC) from all remaining vectors
        DO JVEC = IVEC+1,NDIM
          XSX = INPROD(SCR1,X(1,JVEC),NDIM)
          CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,IVEC),ONE,-XSX,NDIM) 
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization information '
        WRITE(6,*) ' Overlap matrix '
        CALL WRTMAT(S,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Orthonormalized vectors '
        CALL WRTMAT(X,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
*
      SUBROUTINE SETDIA(MATRIX,VALUE,NDIM,IPACK)
*
* Set diagonal elements of matrix MATRIX to VALUE
*
* IPACK = 0 => full quadratic matrix
* IPACK = 1 => lower triangular matrix, row packed
*
#include "implicit.inc"
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
*
      SUBROUTINE GET_NONSYM_SUBMAT(LUC,LUHC,NVEC,NVECP,SUBMAT,NDIM,
     &                             VEC1,VEC2,SUBMATP)
*
* Obtain subspace matrix from trialvectors and matrix times 
* trial vectors. No assumption about symmetry of matrix
*
* Jeppe Olsen, May 2000
*
* Version assuming the ability to hold two vectors in core
*
#include "implicit.inc"
      REAL*8 INPROD
*. Input : Previous subspace matrix
      DIMENSION SUBMATP(*)
*. Output : New subspace matrix 
      DIMENSION SUBMAT(*)
*. Scratch 
      DIMENSION VEC1(*),VEC2(*)
*
      LBLK = -1
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Input vectors to GET_NONSYM_SUBMAT '
        WRITE(6,*) ' ==================================='
*
        WRITE(6,*) ' C vectors '
        CALL REWINO(LUC)
        DO IVEC = 1, NVEC
          CALL FRMDSC(VEC1,NDIM,LBLK,LUC,IMZERO,IAMPACK)
          WRITE(6,*) ' Input C vector number = ', IVEC
          CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END DO
*
        WRITE(6,*) ' HC vectors '
        CALL REWINO(LUHC)
        DO IVEC = 1, NVEC
          CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
          WRITE(6,*) ' Input HC vector number = ', IVEC
          CALL WRTMAT(VEC1,1,NDIM,1,NDIM)
        END DO
*
      END IF
*
      IMZERO = 0
      IAMPACK = 0
*.Reform previous matrix from (NVECP,NVECP) to (NVEC,NVEC) format
      IF(NVECP.GT.0) THEN
C            COPMAT(AIN,AOUT,NIN,NOUT)
        CALL COPMAT(SUBMATP,SUBMAT,NVECP,NVEC)
      END IF
*
*. Add new columns of subspace matrix
*
C?    WRITE(6,*) ' NVEC, NVECP = ', NVEC,NVECP
      CALL REWINO(LUHC)
      DO IVEC = 1, NVECP
        CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
      END DO
      DO J = 1, NVEC-NVECP
*. Read in Hc(nvecp+j)
C?      WRITE(6,*) ' J = ', J
        CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
        CALL REWINO(LUC)
        DO I = 1, NVEC
C?        WRITE(6,*) ' I = ', I
*. Read in c(i)
          CALL FRMDSC(VEC2,NDIM,LBLK,LUC,IMZERO,IAMPACK)
* Submat(i,nvecp+j) = c(i)t Hc(nvecp+j)
          SUBMAT((NVECP+J-1)*NVEC+I) = INPROD(VEC2,VEC1,NDIM)
        END DO
      END DO
*
*. Add new rows of subspace matrix 
*
C?    WRITE(6,*) ' Part 2 ' 
      CALL REWINO(LUC)
      DO I = 1, NVECP
        CALL FRMDSC(VEC2,NDIM,LBLK,LUC,IMZERO,IAMPACK)
      END DO
C     CALL SKPRC3(NVECP,LUC)
      DO I = 1, NVEC-NVECP
C?      WRITE(6,*) ' I = ', I
*. Read in c(nvecp+i)
        CALL FRMDSC(VEC2,NDIM,LBLK,LUC,IMZERO,IAMPACK)
        CALL REWINO(LUHC)
        DO J = 1, NVECP
C?        WRITE(6,*) ' J = ', J
*. Read in Hc(j)
          CALL FRMDSC(VEC1,NDIM,LBLK,LUHC,IMZERO,IAMPACK)
* Submat(i,nvecp+j) = c(i)t Hc(nvecp+j)
          SUBMAT((J-1)*NVEC+I+NVECP) = INPROD(VEC2,VEC1,NDIM)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated subspace matrix '
        CALL WRTMAT(SUBMAT,NVEC,NVEC,NVEC,NVEC)
      END IF
*
      RETURN
      END
*
      SUBROUTINE COPMAT(AIN,AOUT,NIN,NOUT)
C
C COPY MATRIX AIN OF DIMENSION NIN,NIN INTO
C      MATRIX AOUT OF DIMENSAION NOUT,NOUT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AIN(NIN,NIN)
      DIMENSION AOUT(NOUT,NOUT)
C
      DO 100 J = 1, NOUT
       CALL COPVEC(AIN(1,J),AOUT(1,J),NOUT)
  100 CONTINUE
C
      RETURN
      END
*
      SUBROUTINE EIGGMT3(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
     &                   Z,W,SCR,IORD)
*
* Outer routine for calculating eigenvectors and eigenvalues
* of a general real matrix
*
* Version employing EISPACK path RG
*
* Current implementation is rather wastefull with respect to
* memory but at allows one to work with real arithmetic
* outside this routine
*
* If IORD.EQ.1, the eigenvalues are oredered according to the 
* size of the real part of the eigenvalues
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8 INPROD
      DIMENSION AMAT(NDIM,NDIM),SCR(*)
      DIMENSION ARVAL(NDIM),AIVAL(NDIM)
      DIMENSION ARVEC(NDIM,NDIM),AIVEC(NDIM,NDIM)
      DIMENSION Z(NDIM,NDIM),W(NDIM)
*
* Diagonalize
*
      NSCR = 2*NDIM
      CALL RG(NDIM,NDIM,AMAT,ARVAL,AIVAL,1,Z,SCR(1),SCR(1+NDIM),IERR)
      IF( IERR.NE.0) THEN
        WRITE(6,*) ' Problem in EIGGMTN, no convergence '
        WRITE(6,*) ' I have to stop '
        STOP ' No convergence in EIGGMTN '
      END IF
*
* Extract real and imaginary parts according to Eispack manual p.89
*
      DO 150 K = 1, NDIM
*
        IF(AIVAL(K).NE.0.0D0) GOTO 110
        CALL COPVEC(Z(1,K),ARVEC(1,K),NDIM)
        CALL SETVEC(AIVEC(1,K),0.0D0,NDIM)
        GOTO 150
*
  110   CONTINUE
        IF(AIVAL(K).LT.0.0D0) GOTO 130
        CALL COPVEC(Z(1,K),ARVEC(1,K),NDIM)
        CALL COPVEC(Z(1,K+1),AIVEC(1,K),NDIM)
        GOTO 150
*
  130   CONTINUE
        CALL COPVEC(ARVEC(1,K-1),ARVEC(1,K),NDIM)
        CALL VECSUM(AIVEC(1,K),AIVEC(1,K),AIVEC(1,K-1),
     &              0.0D0,-1.0D0,NDIM)
*
  150 CONTINUE
*
* explicit orthogonalization of eigenvectors with
* (degenerate eigenvalues are not orthogonalized by DGEEV)
*
      GOTO 201
      TEST = 1.0D-11
      DO 200 IVEC = 1, NDIM
         RNORM = INPROD(ARVEC(1,IVEC),ARVEC(1,IVEC),NDIM)
     &         + INPROD(AIVEC(1,IVEC),AIVEC(1,IVEC),NDIM)
         FACTOR = 1.0d0/SQRT(RNORM)
         CALL SCALVE(ARVEC(1,IVEC),FACTOR,NDIM)
         CALL SCALVE(AIVEC(1,IVEC),FACTOR,NDIM)
         DO 190 JVEC = IVEC+1,NDIM
           IF(ARVAL(IVEC)-ARVAL(JVEC).LE.TEST) THEN
* orthogonalize jvec to ivec
           OVLAPR = INPROD(ARVEC(1,IVEC),ARVEC(1,JVEC),NDIM)
     &            + INPROD(AIVEC(1,JVEC),AIVEC(1,IVEC),NDIM)
           OVLAPI = INPROD(ARVEC(1,IVEC),AIVEC(1,JVEC),NDIM)
     &            - INPROD(AIVEC(1,IVEC),ARVEC(1,JVEC),NDIM)
           CALL VECSUM(ARVEC(1,JVEC),ARVEC(1,JVEC),ARVEC(1,IVEC),
     &                 1.0D0,-OVLAPR,NDIM )
           CALL VECSUM(AIVEC(1,JVEC),AIVEC(1,JVEC),AIVEC(1,IVEC),
     &                 1.0D0,-OVLAPR,NDIM )
           CALL VECSUM(ARVEC(1,JVEC),ARVEC(1,JVEC),AIVEC(1,IVEC),
     &                 1.0D0,OVLAPI,NDIM )
           CALL VECSUM(AIVEC(1,JVEC),AIVEC(1,JVEC),ARVEC(1,IVEC),
     &                 1.0D0,-OVLAPI,NDIM )
         END IF
  190    CONTINUE
  200 CONTINUE
  201 CONTINUE
*
* Normalize eigenvectors
*
      DO 300 L = 1, NDIM
        XNORM = INPROD(ARVEC(1,L),ARVEC(1,L),NDIM)
     &        + INPROD(AIVEC(1,L),AIVEC(1,L),NDIM)
        FACTOR = 1.0D0/SQRT(XNORM)
        CALL SCALVE(ARVEC(1,L),FACTOR,NDIM)
        CALL SCALVE(AIVEC(1,L),FACTOR,NDIM)
  300 CONTINUE
*
* Order eigensolutions after size of real part of A
*
      IF(IORD.EQ.1) THEN
*. Get reorder array
        CALL DUMSORT(ARVAL,NDIM,NDIM,W,SCR)
C            DUMSORT(VEC,NDIM,NELMNT,IREO,ISCR)      
*. Reorder eigenvalues
C  REO_VEC(IREO,VECIN,NDIM,VECOUT,IWAY)
        CALL REO_VEC(W,ARVAL,NDIM,SCR,1)
        CALL COPVEC(SCR,ARVAL,NDIM)
        CALL REO_VEC(W,AIVAL,NDIM,SCR,1)
        CALL COPVEC(SCR,AIVAL,NDIM)
*. Reorder eigenvectors
C  REO_COL_MAT(IREO,AIN,AOUT,NR,NC,IWAY)
        CALL REO_COL_MAT(W,ARVEC,Z,NDIM,NDIM,1)
        CALL COPVEC(Z,ARVEC,NDIM**2)
        CALL REO_COL_MAT(W,AIVEC,Z,NDIM,NDIM,1)
        CALL COPVEC(Z,AIVEC,NDIM**2)
*
      END IF
*
      NTEST = 0
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Output from EIGGMT '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Real and imaginary parts of eigenvalues '
        CALL WRTMAT_EP(ARVAL,1,NDIM,1,NDIM)
        CALL WRTMAT_EP(AIVAL,1,NDIM,1,NDIM)
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' real part of eigenvectors '
        CALL WRTMAT(ARVEC,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' imaginary part of eigenvectors '
        CALL WRTMAT(AIVEC,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
*
      SUBROUTINE WRTMAT_EP(A,NROW,NCOL,NMROW,NMCOL)
*
* Print matrix, extended precision (E25.15)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(A(I,J),J=1,NCOL)
 1010 FORMAT(1H0,I3,2X,2(1X,E25.15),/,(1H ,5X,2(1X,E25.15)))
  100 CONTINUE
      RETURN
      END
*
      SUBROUTINE REO_VEC(IREO,VECIN,NDIM,VECOUT,IWAY)
*
* Reorder vector
*
* IWAY = 1 : VECOUT(I) = VECIN(IREO(I))
* IWAY = 2 : VECOUT(IREO(I)) = VECIN(I)
*
* Jeppe Olsen, May 2000
*
#include "implicit.inc"
*. Input
      DIMENSION VECIN(NDIM)
      INTEGER IREO(NDIM)
*. Output
      DIMENSION VECOUT(NDIM)
*
      IF(IWAY.EQ.1) THEN
        DO I = 1, NDIM
          VECOUT(I) = VECIN(IREO(I))
        END DO
      ELSE
        DO I = 1, NDIM
          VECOUT(IREO(I)) = VECIN(I)
        END DO 
      END IF
*
      RETURN
      END
*
      SUBROUTINE REO_COL_MAT(IREO,AIN,AOUT,NR,NC,IWAY)
*
* Reorder columns of matrix
*
* IWAY = 1 : AOUT(I,J) = AIN(I,IREO(J))
* IWAY = 2 : AOUT(I,IREO(J)) = AIN(I,J)
*
* Jeppe Olsen, May of 2000
*
#include "implicit.inc"
*. Input
      DIMENSION AIN(NR,NC) 
      INTEGER IREO(*)
*. Output
      DIMENSION AOUT(NR,NC)
*
      IF(IWAY.EQ.1) THEN
        DO J = 1, NC
          CALL COPVEC(AIN(1,IREO(J)),AOUT(1,J),NR)
        END DO
      ELSE
        DO J = 1, NC
          CALL COPVEC(AIN(1,J),AOUT(1,IREO(J)),NR)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE H0M1TV_REL(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal Diag
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements in DIAG corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , May 1990

*
*
* =====
* Input
* =====
* DIAG : Diagonal of matrix
* VECIN : Input vector
* NVAR : Dimension of full matrix
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
* ======
* Output
* ======
* VECUT : Output vector (you guessed ?? ), can occupy same space
*         as VECIN or DIAG
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      REAL * 8  INPROD
*
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
C?    write(6,*) ' H0M1TV , NPQDM = ', NPQDM
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
        CALL GATVEC(WORK(KLV1),VECIN,IPNTR,NPQDM)
* X(T)(DIAG+SHIFT)-1 X in subspace, for later subtraction
        CALL GATVEC(WORK(KLV2),DIAG,IPNTR,NPQDM)
        CALL DIAVC3(WORK(KLV2),WORK(KLV1),
     &       WORK(KLV2),SHIFT,NPQDM,X1)
       ELSE
         X1 = 0.0D0
       END IF
*
      CALL DIAVC3(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
*
      IF(NPQDM .NE. 0 ) THEN
C                H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
C    &           X,RHS,S,SCR)
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
C?     write(6,*) ' KLPHP KLPHQ KLQHQ ',KLPHP,KLPHQ,KLQHQ
*
         CALL H0LNSL_REL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
         X3 = INPROD(WORK(KLV1),WORK(KLV2),NPQDM)
         CALL SCAVEC(VECUT,WORK(KLV2),IPNTR,NPQDM)
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(6,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3


*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TV '
        WRITE(6,*) ' ========================= '
        CALL WRTMAT(VECUT,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
*
      SUBROUTINE H0LNSL_REL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS

*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVER
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM_REL
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL QENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          CALL TRIPAK(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,1.0D0)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM,1.0D0)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMAT(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM,ISING)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================

         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LUDIA = 37
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LUDIA)
*. Initial Guess
        CALL REWINE(LU1,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LU1)
*. Right hand side
        CALL REWINE(LU2,-1)
        CALL TODSC(RHS,NPQDM,-1,LU2)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM_REL,LU1,LU2,LU3,LUDIA,SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LU1,-1)
        CALL FRMDSC(SCR(KLV1),NPQDM,-1,LU1,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      CALL QEXIT('H0LNS')
      RETURN
      END
*
      SUBROUTINE HPQTVM_REL(VECIN,VECUT,IDUM,JDUM,
     &                      WORK,KFREE,LFREE)
*
* Outer routine for Preconditioner times vector
*
*. LUCIA version
#include "implicit.inc"
#include "mxpdim.inc"
#include "glbbas.inc"
#include "crun.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
C#include "wrkspc.inc"
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,
C    &            ICISTR,NOCSF,IDIAG

      DIMENSION WORK(*)


      COMMON/SHFT/SHIFT
*
      NPDM = MXP1+MXP2
      NQDM = MXQ
      KLPHP = KH0
      KLPHQ = KH0 + NPDM*(NPDM+1)/2
      KLQHQ = KLPHQ + MXP1*MXQ
*
C          HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
      CALL HPQTV(MXP1,MXP2,MXQ,
     &           WORK(KLPHP),WORK(KLPHQ),WORK(KLQHQ),
     &           VECIN,VECUT)
*
      IF(SHIFT.NE.0.0D0) THEN
       CALL VECSUM(VECUT,VECUT,VECIN,1.0D0,SHIFT,NPDM+NQDM)
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTVM'
        CALL WRTMAT(VECIN,1,NPDM+NQDM,1,NPDM+NQDM)
        CALL WRTMAT(VECUT,1,NPDM+NQDM,1,NPDM+NQDM)
      END IF
*
      RETURN
      END
*
      SUBROUTINE MV10_REL(CCIN,CCOUT,WORK,KFREE,LFREE)
*
* Initial outer routine for cc Jacobian times vector
*
* Jeppe Olsen, June 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "glbbas.inc"
C#include "wrkspc.inc" !pass down work
#include "ctcc.inc"
#include "ctccp.inc"
C#include "lorr.inc"
#include "crun.inc"
C     COMMON/LORR/L_OR_R
*
      DIMENSION CCIN(LEN_T_VEC_MX),CCOUT(LEN_T_VEC_MX),WORK(*)
*
      ZERO = 0.0D0
      NTEST = 00
C
      CALL SETVEC(CCOUT,ZERO,LEN_T_VEC_MX)
C
      CALL KRCC_VEC_FNC_LR(WORK(KCC1),CCIN,CCOUT,WORK,KFREE,LFREE)
*WORK(KCC1) are the CC amplitudes 
*CCIN is the right vector VR 
          IF(NTEST.GE.100) THEN
            WRITE(6,*)
            WRITE(6,*) ' ====================================== '
            WRITE(6,*) ' Input and output KRCC LRCC' 
            WRITE(6,*) ' ====================================== '
            WRITE(6,*) 
            CALL WRTMAT(CCIN,1,LEN_T_VEC_MX,1,LEN_T_VEC_MX)
            WRITE(6,*)
            CALL WRTMAT(CCOUT,1,LEN_T_VEC_MX,1,LEN_T_VEC_MX)
          END IF
*
      RETURN
      END 
*
      SUBROUTINE KRCC_VEC_FNC_LR(CC_AMP,CC_VEC,
     &                          WORK,KFREE,LFREE)
*
* Calculate energy and CC vector function for 
* a set of CC amplitudes defined by CC_AMP
*
* New version by Lasse 2010
*
* Now for all possible contractions in the (p,h) scheme
* Shown as Operator --Type of contarction with T --> New operator
* H for Hamiltonian
* M for intermediate (easier to see)
* T for Cluster operator (no index needed)
* Indices are free of the ph definition
*
* This version is programmed as shown below. This is not yet the optimum
* but can "easily" be changed. These are in execution order without
* scaling considerations.
*
* First H00 -- OK
*
* H00 -- 00 -- T 
*
* Second all single lines -- OK
*
* H22 -- 22 -- T
* H21 -- 21 -- T
* H12 -- 12 -- T
* H20 -- 20 -- T
*
* Third through M10
*
* H20 -- 10 -- M10 (multiply by 1/2)
* H22 -- 12 -- M10
* H21 -- 11 -- M10
*
* M10 + H10 -- 10 -- T
*
* Fourth through M02
*
* H22 -- 10 -- M12 (multiply by 1/2)
* M12 + H12 -- 10 -- M02 
* H22 -- 20 --  M02 
* M02 + H02 -- 02 -- T (save M02 + H02)
*
* Fifth through M11
*
* H22 -- 11 -- M11 (multiply by 1/2, save M11 before multiplication)
* M11 + H11 -- 11 -- T 
*
* Sixth through M01
*
* H21 -- 10 -- M11 (multiply by 1/2 on additional array(At the moment))
* H22 -- 11 -- M11 (reuse the above M11, add to above)
* M11 + H11 -- 10 -- M01
* H02 + M02 -- 01 -- M01 (reuse the above M02 + H02, multiply by 1/2)
* H22 -- 21 -- M01
* H21 -- 20 -- M01
* H12 -- 11 -- M01
* M01 + H01 -- 01 -- T
*
* I guess that is all we need then!
*
*      
      IMPLICIT REAL*8(A-H,O-Z)
*. Input   
      DIMENSION CC_AMP(*)
*. Output   
      DIMENSION CC_VEC(*)
*
      DIMENSION WORK(*)
*.
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "clunit.inc"
#include "glbbas.inc"
#include "crun.inc"
#include "cands.inc"
#include "cstate.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "cgas.inc"
#include "cecore.inc"
#include "interm.inc"
*
      LBLK = -1
      NTEST = 10
*
* Will determine the number of commutators (Of course one one with zero)
      IONE = 1
      IZERO = 0
      ZERO = 0.0D0
      HALF = 0.5D0
      ONE = 1.0D0
      DOUBLE = 2.0D0 !!! forgot it?
* Debug
      ISKIP = 5 ! from which commutator to skip
*         
* Real or complex algebra
      IF(IRECOM.EQ.1) THEN
        IMULTFAC = 1
      ELSE
        IMULTFAC = 2
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' CC_VEC_FNC in action '
        WRITE(6,*) ' ====================='
      END IF
*
* Hardwire symmetry
      IHSM = 1
      ITSM = 1
*      
*. May want to allocate some temporary arrays
*
      CALL EXP_MT_H_EXP_TC_MEM_KRCC(1,WORK,KFREE,LFREE)
*
*. Start by zeroing the CC-vector function
      CALL SETVEC(CC_VEC,ZERO,N_CC_AMP+1)
*
* Start with zero commutators! (If first iteration and no restart this
* should be the only one done!)
*
* H00 --00--> T SHOULD BE THE ONLY ONE DONE FOR IT=1 AND NO RESTART
* 
* First H00
*
* H00 -- 00 -- T
      IF(ISKIP.GE.1) THEN 
      IHM = 1 !Still the Hamiltonian that is contracted from 
      IFHM = 1 !integrals are fetched
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IZERO,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TP_CC,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH00,WORK(KNH00),
     &     NH00,WORK(KM00TOH00),WORK(KM00TOH00),
     &     IDUM,IDUM,
     &     IDUM,IDUM,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
! First iteration with no restart only do H00
      CALL MEMCHK_KRCC(WORK)
      IF(IFORM.EQ.1.AND.ITER.EQ.1) GO TO 22
      END IF
*
* Second all single lines
*
* H22 -- 22 -- T
*
      IF(ISKIP.GE.2) THEN
      print*,'starting single contractions'
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOT,WORK(KNH22TOT),WORK(KNH22TT),WORK(KCOMBNH22TOT),
     &     WORK(KPERMNH22TOT),WORK(KFACNH22TOT),WORK(KCOMBFACNH22TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
*
* H21 -- 21 -- T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      print*,'doing H21'
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOT,WORK(KNH21TOT),WORK(KNH21TT),WORK(KCOMBNH21TOT),
     &     WORK(KPERMNH21TOT),WORK(KFACNH21TOT),WORK(KCOMBFACNH21TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H12 -- 12 -- T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH12,WORK(KNH12),
     &     NH12TOT,WORK(KNH12TOT),WORK(KNH12TT),WORK(KCOMBNH12TOT),
     &     WORK(KPERMNH12TOT),WORK(KFACNH12TOT),WORK(KCOMBFACNH12TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H20 -- 20 -- T
*
      END IF !ISKIP
      IF(ISKIP.GE.2) THEN
      print*,'single contraction'
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH20,WORK(KNH20),
     &     NH20TOT,WORK(KNH20TOT),WORK(KNH20TT),WORK(KCOMBNH20TOT),
     &     WORK(KPERMNH20TOT),WORK(KFACNH20TOT),WORK(KCOMBFACNH20TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
*
* Third through M10
*
* H20 -- 10 -- M10 (multiply by 1/2) Should be done first!!!
*
* Allocate the M10 intermediate
*
      NPART = 1
      NHOLE = 0 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
*
      IF(ISKIP.GE.3) THEN
*
      print*,'to intermediate'
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH20,WORK(KNH20),
     &     NH20TOM10,WORK(KNH20TOM10),WORK(KNH20T10),
     &     WORK(KCOMBNH20TOM10),WORK(KPERMNH20TOM10),
     &     WORK(KFACNH20TOM10),WORK(KCOMBFACNH20TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
      print*,'after H20 -- 10 -- M10'
*
* Multiply by 1/2
*
      CALL SCALVE(WORK(KM10),HALF,LEN_M10_VEC_CC)
      CALL MEMCHK_KRCC(WORK)
      END IF !ISKIP
      IF(ISKIP.GE.3) THEN
*
* H22 -- 12 -- M10
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM10,WORK(KNH22TOM10),WORK(KNH22T10),
     &     WORK(KCOMBNH22TOM10),WORK(KPERMNH22TOM10),
     &     WORK(KFACNH22TOM10),WORK(KCOMBFACNH22TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
      IF(ISKIP.GE.3) THEN
*
* H21 -- 11 -- M10
*
      print*,'H21 to M10'
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOM10,WORK(KNH21TOM10),WORK(KNH21T10),
     &     WORK(KCOMBNH21TOM10),WORK(KPERMNH21TOM10),
     &     WORK(KFACNH21TOM10),WORK(KCOMBFACNH21TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* M10 + H10 -- 10 -- T
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      END IF !ISKIP
      IF(ISKIP.GE.2) THEN
      print*,' adding M10 and H10'
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH10,WORK(KNH10),
     &                                   NINTER10,WORK(KINTM10),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM10),WORK(KM10TOH10),
     &                                   WORK(KLIHIND),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      print*,'second contraction M10'
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM10), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &    NINTER10,WORK(KINTM10),WORK(KLIBSOBEX_M10),WORK(KLIHIND),IHSM,
     &     NINTER10,WORK(KLLSOBEX_M10),
     &     NM10TOT,WORK(KNM10TOT),WORK(KNM10TT),
     &     WORK(KCOMBNM10TOT),WORK(KPERMNM10TOT),
     &     WORK(KFACNM10TOT),WORK(KCOMBFACNM10TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
      END IF !ISKIP
*
* Deallocate the M10 intermediate
*
      NPART = 1
      NHOLE = 0 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,0)
*
* Fourth through M02
*
* H22 -- 10 -- M12 (multiply by 1/2)
*
* Allocate the M12 intermediate
*
      NPART = 1
      NHOLE = 2 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
*
      IF(ISKIP.GE.5) THEN
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM12),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M12),WORK(KLLSOBEX_M12),NINTER12,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM12,WORK(KNH22TOM12),WORK(KNH22T12),
     &     WORK(KCOMBNH22TOM12),WORK(KPERMNH22TOM12),
     &     WORK(KFACNH22TOM12),WORK(KCOMBFACNH22TOM12),
     &     WORK(KINTM12),NINTER12,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Multiply by 1/2
*
      CALL SCALVE(WORK(KM12),HALF,LEN_M12_VEC_CC)
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
*
* M12 + H12 -- 10 -- M02 
*
* Add amplitudes and integrals 
*
* Allocate the M02 intermediate
*
      NPART = 0
      NHOLE = 2 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
*
      IF(ISKIP.GE.3) THEN
*
      print*,'adding M12 to H12'
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH12,WORK(KNH12),
     &                                   NINTER12,WORK(KINTM12),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM12),WORK(KM12TOH12),
     &                                   WORK(KLIHIND),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      print*,'doing the contraction'
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM02),IONE,IHM,
     &     IFHM,WORK(KM12), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M02),WORK(KLLSOBEX_M02),NINTER02,
     &    NINTER12,WORK(KINTM12),WORK(KLIBSOBEX_M12),WORK(KLIHIND),IHSM,
     &     NINTER12,WORK(KLLSOBEX_M12),
     &     NM12TOM02,WORK(KNM12TOM02),WORK(KNM12T02),
     &     WORK(KCOMBNM12TOM02),WORK(KPERMNM12TOM02),
     &     WORK(KFACNM12TOM02),WORK(KCOMBFACNM12TOM02),
     &     WORK(KINTM02),NINTER02,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
      END IF !ISKIP
      IF(ISKIP.GE.3) THEN
*
* H22 -- 20 -- M02
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM02),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M02),WORK(KLLSOBEX_M02),NINTER02,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM02,WORK(KNH22TOM02),WORK(KNH22T02),
     &     WORK(KCOMBNH22TOM02),WORK(KPERMNH22TOM02),
     &     WORK(KFACNH22TOM02),WORK(KCOMBFACNH22TOM02),
     &     WORK(KINTM02),NINTER02,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* M02 + H02 -- 02 -- T (save M02 + H02)
*
* Add amplitudes and integrals 
*
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
      IF(ISKIP.GE.2) THEN
      print*,' adding M02 and H02'
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH02,WORK(KNH02),
     &                                   NINTER02,WORK(KINTM02),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM02),WORK(KM02TOH02),
     &                                   WORK(KLIHIND),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* Save vector to disc since it will be needed later
*
      CALL VEC_TO_DISC(WORK(KM02),LEN_M02_VEC_CC,1,LBLK,LUSC5)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM02), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &    NINTER02,WORK(KINTM02),WORK(KLIBSOBEX_M02),WORK(KLIHIND),IHSM,
     &     NINTER02,WORK(KLLSOBEX_M02),
     &     NM02TOT,WORK(KNM02TOT),WORK(KNM02TT),
     &     WORK(KCOMBNM02TOT),WORK(KPERMNM02TOT),
     &     WORK(KFACNM02TOT),WORK(KCOMBFACNM02TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
      END IF !ISKIP
*
* Deallocate the M12 intermediate (Will also deallocate M02)
*
      NPART = 1
      NHOLE = 2
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,0)
*
* Fifth through M11
*
* H22 -- 11 -- M11 (multiply by 1/2)
*
* Allocate the M11 intermediate
*
      NPART = 1
      NHOLE = 1 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
*
      IF(ISKIP.GE.3) THEN
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
      print*,'H22 to M11 is now in session'
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM11),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M11),WORK(KLLSOBEX_M11),NINTER11,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM11,WORK(KNH22TOM11),WORK(KNH22T11),
     &     WORK(KCOMBNH22TOM11),WORK(KPERMNH22TOM11),
     &     WORK(KFACNH22TOM11),WORK(KCOMBFACNH22TOM11),
     &     WORK(KINTM11),NINTER11,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Save vector to disc since it will be needed later
*
      CALL VEC_TO_DISC_KRCC(WORK(KM11),LEN_M11_VEC_CC,1,
     &                      I_READ_BLOCK,LUSC4)
*
* Multiply by 1/2
*
      CALL SCALVE(WORK(KM11),HALF,LEN_M11_VEC_CC)
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
      IF(ISKIP.GE.2) THEN
*
* M11 + H11 -- 11 -- T 
*
* Add amplitudes and integrals 
*
      print*,' adding M11 and H11'
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH11,WORK(KNH11),
     &                                   NINTER11,WORK(KINTM11),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM11),WORK(KM11TOH11),
     &                                   WORK(KLIHIND),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
      print*,'after add for M11'
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM11), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &    NINTER11,WORK(KINTM11),WORK(KLIBSOBEX_M11),WORK(KLIHIND),IHSM,
     &     NINTER11,WORK(KLLSOBEX_M11),
     &     NM11TOT,WORK(KNM11TOT),WORK(KNM11TT),
     &     WORK(KCOMBNM11TOT),WORK(KPERMNM11TOT),
     &     WORK(KFACNM11TOT),WORK(KCOMBFACNM11TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      END IF !ISKIP
*
* Deallocate the M11 intermediate 
*
      NPART = 1
      NHOLE = 1
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,0)
*
      CALL MEMCHK_KRCC(WORK)
*
* Allocate the M01 intermediate. Since in this way we will never need
* more than two intermediates in memory
*
      NPART = 0
      NHOLE = 1 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
*
* Allocate the M11 intermediate
*
      NPART = 1
      NHOLE = 1 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
      IF(ISKIP.GE.4) THEN
*
* Sixth through M01
*
* H21 -- 10 -- M11
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM11),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M11),WORK(KLLSOBEX_M11),NINTER11,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOM11,WORK(KNH21TOM11),WORK(KNH21T11),
     &     WORK(KCOMBNH21TOM11),WORK(KPERMNH21TOM11),
     &     WORK(KFACNH21TOM11),WORK(KCOMBFACNH21TOM11),
     &     WORK(KINTM11),NINTER11,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Multiply by 1/2
*
      CALL SCALVE(WORK(KM11),HALF,LEN_M11_VEC_CC)
*
* H22 -- 11 -- M11 (reuse the above M11). Add KM11 and LUSC4
*
      CALL ADD_VEC_FROM_DISC_TO_VEC(WORK(KM11),LEN_M11_VEC_CC,1,
     &                              I_READ_BLOCK,LUSC4)
      CALL MEMCHK_KRCC(WORK)
C
      END IF !ISKIP
      IF(ISKIP.GE.3) THEN
*
* M11 + H11 -- 10 -- M01
*
* Add amplitudes and integrals 
      print*,' adding M11 and H11 Second time'
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH11,WORK(KNH11),
     &                                   NINTER11,WORK(KINTM11),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM11),WORK(KM11TOH11),
     &                                   WORK(KLIHIND),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
      print*,'does the 11 -- 10 -- M01 contraction'
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KM11), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &    NINTER11,WORK(KINTM11),WORK(KLIBSOBEX_M11),WORK(KLIHIND),IHSM,
     &     NINTER11,WORK(KLLSOBEX_M11),
     &     NM11TOM01,WORK(KNM11TOM01),WORK(KNM11T01),
     &     WORK(KCOMBNM11TOM01),WORK(KPERMNM11TOM01),
     &     WORK(KFACNM11TOM01),WORK(KCOMBFACNM11TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Deallocate the M11 intermediate 
*
      NPART = 1
      NHOLE = 1
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,0)
*
      END IF !ISKIP
      IF(ISKIP.GE.3) THEN
*
* Allocate the M02 intermediate
*
      NPART = 0
      NHOLE = 2 
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,1)
*
* H02 + M02 -- 01 -- M01 (reuse the above M02 + H02, multiply by 1/2)
*
      CALL VEC_FROM_DISC(WORK(KM02),LEN_M02_VEC_CC,1,LBLK,LUSC5)
*
* Multiply by 1/2
*
      CALL SCALVE(WORK(KM02),HALF,LEN_M02_VEC_CC)
* 
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      print*,'intermediate H02 is being fetched'
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KM02), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &    NINTER02,WORK(KINTM02),WORK(KLIBSOBEX_M02),WORK(KLIHIND),IHSM,
     &     NINTER02,WORK(KLLSOBEX_M02),
     &     NM02TOM01,WORK(KNM02TOM01),WORK(KNM02T01),
     &     WORK(KCOMBNM02TOM01),WORK(KPERMNM02TOM01),
     &     WORK(KFACNM02TOM01),WORK(KCOMBFACNM02TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Deallocate the M02 intermediate 
*
      NPART = 0
      NHOLE = 2
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,0)
*
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
      IF(ISKIP.GE.3) THEN
*
* H22 -- 21 -- M01
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM01,WORK(KNH22TOM01),WORK(KNH22T01),
     &     WORK(KCOMBNH22TOM01),WORK(KPERMNH22TOM01),
     &     WORK(KFACNH22TOM01),WORK(KCOMBFACNH22TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      END IF !ISKIP
      IF(ISKIP.GE.3) THEN
*
* H21 -- 20 -- M01
*
      print*,' H21 to M01'
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOM01,WORK(KNH21TOM01),WORK(KNH21T01),
     &     WORK(KCOMBNH21TOM01),WORK(KPERMNH21TOM01),
     &     WORK(KFACNH21TOM01),WORK(KCOMBFACNH21TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      END IF !ISKIP
      IF(ISKIP.GE.3) THEN
*
* H12 -- 11 -- M01
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),WORK(KLIHIND),IHSM,
     &     NH12,WORK(KNH12),
     &     NH12TOM01,WORK(KNH12TOM01),WORK(KNH12T01),
     &     WORK(KCOMBNH12TOM01),WORK(KPERMNH12TOM01),
     &     WORK(KFACNH12TOM01),WORK(KCOMBFACNH12TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* M01 + H01 -- 01 -- T
*
* Add amplitudes and integrals 
*
      END IF !ISKIP
      CALL MEMCHK_KRCC(WORK)
      IF(ISKIP.GE.2) THEN
      print*,' adding M01 and H01'
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH01,WORK(KNH01),
     &                                   NINTER01,WORK(KINTM01),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM01),WORK(KM01TOH01),
     &                                   WORK(KLIHIND),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      print*,'second contraction for M01'
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM01), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &    NINTER01,WORK(KINTM01),WORK(KLIBSOBEX_M01),WORK(KLIHIND),IHSM,
     &     NINTER01,WORK(KLLSOBEX_M01),
     &     NM01TOT,WORK(KNM01TOT),WORK(KNM01TT),
     &     WORK(KCOMBNM01TOT),WORK(KPERMNM01TOT),
     &     WORK(KFACNM01TOT),WORK(KCOMBFACNM01TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Deallocate the M01 intermediate 
*
      NPART = 0
      NHOLE = 1
      CALL IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,0)
*
      END IF !ISKIP
22    CONTINUE
      CALL MEMCHK_KRCC(WORK)
* End of the line
*
C     E_CC = CC_VEC(N_CC_AMP+1) IF IDENTITY WAS LAST BUT MINE IS FIRST
      E_CC = CC_VEC(1)
*
*. Mission over, print 
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' CC-vector function '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' DO NOT USE WRT_CC_VEC2 YET '
        STOP ' DO NOT USE WRT_CC_VEC2 YET '
        CALL WRT_CC_VEC2(CC_VEC,6,CCTYPE,
     &                   WORK,KFREE,LFREE)
        WRITE(6,*) ' Initial amplitudes '
        CALL WRT_CC_VEC2(CC_AMP,6,CCTYPE,
     &                   WORK,KFREE,LFREE)
      END IF
*
* Release temp mem
*
      CALL EXP_MT_H_EXP_TC_MEM_KRCC(0,WORK,KFREE,LFREE)
* 
      RETURN
      END
*
****************************************************************
*                                                              *
* Silly subroutine by Lasse to find largest difference between *
* neighbouring elements in an array.                           *
*                                                              *
****************************************************************
      SUBROUTINE FIND_LARGEST_DIFF_FOR_NEIGHBOURS(NOUT,NARR,NLENGHT)
      IMPLICIT NONE
      INTEGER NOUT,NLENGHT,NARR(NLENGHT),I,J
      NOUT = 0
      DO I=2,NLENGHT
        J = NARR(I)-NARR(I-1)
        IF(J.GT.NOUT) NOUT = J
      END DO
      END
