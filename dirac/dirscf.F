C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psidhf */
      SUBROUTINE PSIDHF(WORK,LWORK)
C*****************************************************************************
C
C     Solve DHF - equations; print results
C
C     Written by T.Saue November 1994
C     Last revision : Sep 11 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#endif
#if defined (VAR_PVM)
#include <fpvm3h.h>
#endif
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <mxorb.h>
C
#include <ccom.h>
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <blocks.h>
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PSIDHF')
#include <memint.h>
C
      CALL TITLER('Dirac-Fock calculation','*',125)
C
C     Get hold of the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = 1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
      CALL FLSHFO(LUPRI)
C
C     Initialize
C     ==========
C
      CALL SETDHF(IPRSCF)
C
C     Memory allocation
C
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,      WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFOCK,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2TMT         ,WORK,KFREE,LFREE)
C
C     Solve DHF-equations
C     ===================
C
      CALL DHFSCF(WORK(KCMO),WORK(KEIG),WORK(KFOCK),
     &            WORK(KTMAT),WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('PSIDHF.DHFSCF',WORK,KWORK,KFOCK,KFREE,LFREE)
C
C     Output section
C
      CALL DHFOUT(WORK(KCMO),WORK(KEIG),WORK(KFREE),LFREE,IPRINT)
C
C     Release the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = -1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
C
C     Memory deallocation
C
      CALL MEMREL('PSIDHF.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('PSIDHF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck predhf */
      SUBROUTINE PREDHF(CMO,EIG,FOCK,TMAT,DMAT,BMAT,
     &                  NPOS,WORK,LWORK)
C*****************************************************************************
C
C     Obtain orthonormal trial vectors for Dirac-Fock
C     either by:
C       i) bare nucleus approximation
C      ii) starting from trial vectors
C     iii) starting from two-electron Fock matrix
C     Also serves as solver module for the Dirac-equation in the
C     finite basis approximation (ONESYS)
C
C     Numerous possible schemes for restart has been implemtented.
C     They are based on the program looking for necessary files.
C     Three modes should be discerned:
C
C     1) One-electron system (ONESYS)
C     2) Many-electron system
C         a) DFCYCL present: 
C            The program will restart on as much as possible.
C            (coefficients, convergence acceleration, 
C             differential density etc.)
C         b) DFCYCL abscent: 
C            Restart is only possible on coefficients/Fock matrix.
C
C     Written by T.Saue November 1994
C     Last revision : Jan 27 1997 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0,D0 = 0.0D0)
C
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbham.h>
#include <dcbdhf.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
      LOGICAL TOBE,CHKCNV,DOSEL,DSCFON,BARNUC
      CHARACTER CPUTID*12,DAYTID*24,SECTID*12
      DIMENSION FOCK(*),EIG(*),CMO(*),TMAT(*),
     &          BMAT(MXDIIS,MXDIIS),DMAT(*),
     &          NPOS(*),WORK(LWORK)
C
C     ****************************************************
C     ***  I N I T I A L I Z E  a n d   R E S T A R T ****
C     ****************************************************
C
      BARNUC  = .TRUE.
      CACC    = '        '
      DOSEL   = .FALSE.
      DSCFON  = .FALSE.
      CHKCNV  = .FALSE.
      KITER   = 1
C
C     1. One-electron system
C     ======================
C
      IF(ONESYS) THEN
        TRIVEC = .FALSE.
        TRIFCK = .FALSE.
        NITER  = 1
      ELSE
C
C     2. Many-electron system
C     =======================
C
C     * Check for trial function
C     --------------------------
C
        IF(TRIFCK) THEN
          INQUIRE(FILE='DFFCK2',EXIST=TOBE)
          IF(TOBE) THEN
            INTTYP = 'Old FOCK2'
            BARNUC = .FALSE.
          ELSE    
            WRITE(LUPRI,'(/A/)')
     &      '*** WARNING *** No trial 2-electron Fock matrix found.'//
     &      'Looking for trial vectors instead'
            TRIFCK = .FALSE.
            TRIVEC = .TRUE.
          ENDIF
        ENDIF
        IF(TRIVEC) THEN
          INQUIRE(FILE='DFCOEF',EXIST=TOBE)
          IF(TOBE) THEN
            BARNUC = .FALSE.
          ELSE
            WRITE(LUPRI,'(/A/)')
     &      '*** WARNING *** No trial vectors found.'//
     &      'Using bare nucleus approximation instead.'
            TRIVEC = .FALSE.
          ENDIF
        ENDIF
C
C     * Check for cycle history on DFCYCL
C     -----------------------------------
C
        INQUIRE(FILE='DFCYCL',EXIST=TOBE)
        IF(TOBE) THEN
C
C         a) DFCYCL is present
C         --------------------
C
          OPEN(LUCYCL,FILE ='DFCYCL',STATUS='OLD',
     &                ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
          REWIND LUCYCL
 10       CONTINUE
          READ(LUCYCL,1010,END=20) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,INTTYP,DAYTID
          GOTO 10
 20       CONTINUE
#if defined (SYS_IBM) || defined (SYS_CRAY) || defined (VAR_MFDS)
C        backspace on multifile systems for correct positioning
C        at the end of the file
          BACKSPACE LUCYCL
#endif
          WRITE(LUPRI,'(A)') '*** WARNING *** : SCFCYC restart info:'
          WRITE(LUPRI,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,INTTYP,DAYTID
          ILLINT = 0
          IF(INTTYP(1:2).EQ.'LL') ILLINT = 1
          ISLINT = 0
          IF(INTTYP(4:5).EQ.'SL') ISLINT = 1
          ISSINT = 0
          IF(INTTYP(7:8).EQ.'SS') ISSINT = 1
          INTBUF = ILLINT + 2*ISLINT + 4*ISSINT
          NITER = NITER + 1
          CALL INTCON(INTFLG,INTBUF,INTDEF,
     &                ERGVAL,CNVINT,NITER,ITRINT,INTTYP)
C
C         * Check on convergence acceleration
C         -----------------------------------
C
          INQUIRE(FILE='DFFOCK',EXIST=CHKCNV)
          CHKCNV = CHKCNV.AND.(DODIIS.OR.DODAMP)
          IF(DODIIS) THEN
            DIISON = CHKCNV.AND.INTFLG.EQ.INTBUF.AND.NITER.GT.1
            INQUIRE(FILE='DFDIIS',EXIST=TOBE)
            DIISON = DIISON.AND.TOBE
            INQUIRE(FILE='DFCMOS',EXIST=TOBE)
            DIISON = DIISON.AND.TOBE
            INQUIRE(FILE='DFEVEC',EXIST=TOBE)
            DIISON = DIISON.AND.TOBE
            IF(DIISON) THEN
              CALL OPNFIL(LUDIIS,'DFDIIS','OLD','SCFINT')
              READ(LUDIIS) ITDIIS,MDIIS2,NELMBM,BMAT
              IF(MDIIS2.EQ.MXDIIS) THEN
                CHKCNV = .TRUE.
                WRITE(LUPRI,'(A)') '*** WARNING *** : DIIS RESTART'
                WRITE(LUPRI,'(3X,A,I5)')
     &          '* Previous DIIS iteration        :', ITDIIS,
     &          '* Number of elements in B-matrix :',NELMBM
              ELSE
                ITDIIS = 0
              ENDIF
            ELSE
              ITDIIS = 0
            ENDIF
          ENDIF
C
C         Check for differential density matrix
C
          IF(DODSCF.AND.(INTBUF.EQ.INTFLG)) THEN
            INQUIRE(FILE='DFFCK2',EXIST=DSCFON)
            INQUIRE(FILE='DFDENS',EXIST=TOBE)
            DSCFON = DSCFON.AND.TOBE
            IF(DSCFON) THEN
              WRITE(LUPRI,'(A)') 
     &        '*** WARNING *** : Differential density restart'
            ENDIF
          ENDIF
        ELSE
C
C         b) DFCYCL is abscent
C         --------------------
C
          NITER  = 1
          OPEN(LUCYCL,FILE ='DFCYCL',STATUS='NEW',
     &                ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
          IF(TRIVEC) CALL INTCON(INTFLG,INTBUF,INTDEF,
     &                    DUMMY,CNVINT,NITER,ITRINT,INTTYP)
        ENDIF
C
C       Initiate SCF cycle
C
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +  '########## START ITERATION NO.',NITER,'##########',
     +  DAYTID
C
C       Open files necessary for Dirac-Fock calculation
C
        OPEN (LUEVEC, FILE='DFEVEC', FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*N2EVECTQ,STATUS='UNKNOWN')
        CALL OPNFIL(LUFCK2,'DFFCK2','UNKNOWN','PREDHF')
        CALL OPNFIL(LUDENS,'DFDENS','UNKNOWN','PREDHF')
        CALL OPNFIL(LUCMOS,'DFCMOS','UNKNOWN','PREDHF')
      ENDIF
      OPEN (LUFOCK, FILE='DFFOCK', FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*N2ORBTQ,STATUS='UNKNOWN')
      CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','PREDHF')
      CALL OPNFIL(LUFCK1,'DFFCK1','OLD','PREDHF')
C
C     *****************************************************************
C     *****  G E N E R A T E   F U L L    F O C K    M A T R I X  *****
C     *****************************************************************
C
      TID   = SECOND()
C
C     Get MO-transformation matrix
C     ============================
C
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','PREDHF')
      CALL READT(LUTMAT,N2TMT,TMAT)
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Get 2-electron Fock matrix (save in CMO)
C     ========================================
C
C     1) - using trial vectors
C     -------------------------
C
      IF(TRIVEC) THEN
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,TOTERG,2)
C
C       Generate density matrix
C
        CALL DENSTY(DMAT,CMO,IPRSCF)
C
C       Differential density matrix approach
C
        TDF2 = SECOND()
        IF(DSCFON) THEN
          WRITE(LUPRI,'(I5,A)') NITER,
     &       ' *** Differential density matrix ***'
          CALL READNS(LUDENS,FOCK)
          IF(FIXDIF) THEN
            DCOVLP = D1
          ELSE
            DCOVLP = DDOT(N2BBASXQ,FOCK,1,FOCK,1)
            DCOVLP = DDOT(N2BBASXQ,DMAT,1,FOCK,1) / DCOVLP
          ENDIF
          CALL DSCAL(N2BBASXQ,(-DCOVLP),FOCK,1)
          CALL DAXPY(N2BBASXQ,D1,DMAT,1,FOCK,1)
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,CMO,FOCK,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK,LWORK)
          CALL REAFCK(LUFCK2,FOCK,.FALSE.)
          CALL DAXPY(N2BBASXQ,DCOVLP,FOCK,1,CMO,1)
        ELSE
C
C       Standard approach
C
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,CMO,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK,LWORK)
        ENDIF
        IF(DOTSC) CALL TSCORR(CMO,WORK,LWORK,IPRSCF)
        TDF2 = SECOND() - TDF2
C
C       If necessary, generate matrix for overlap selection
C
        IF(OVLSEL) THEN
          DOSEL = .TRUE.
          INQUIRE(FILE='DFSMOS',EXIST=TOBE)
          IF(.NOT.TOBE) CALL SMOGEN(CMO,TMAT,WORK,LWORK)
        ENDIF
        CALL WRIDNS(LUDENS,DMAT)
        CALL WRIFCK(LUFCK2,CMO)
      ENDIF
C
C     2) - using old 2-electron Fock matrix
C     --------------------------------------
C
      IF(TRIFCK) CALL REAFCK(LUFCK2,CMO,.TRUE.)
C
C     Get one-electron Fock matrix (save in FOCK)
C     ===========================================
C
      CALL REAFCK(LUFCK1,FOCK,.FALSE.)
      IF(TRIVEC) CALL ERGCAL(FOCK,CMO,DMAT)
C
C     If necessary, add one- and two- electron Fock matrices
C     ======================================================
C
      IF(TRIVEC.OR.TRIFCK) THEN
        CALL DAXPY(N2BBASXQ,D1,CMO,1,FOCK,1)
      ENDIF
C
C     **************************************************************************
C     *****   S O L V E    D F    E I G E N V A L U E    P R O B L E M    ******
C     **************************************************************************
C
C
      CALL DFSOLV(CMO,EIG,FOCK,DMAT,TMAT,BMAT,CHKCNV,DOSEL,
     &            WORK,LWORK)
      IF(BARNUC.OR.TRIFCK) CALL ONEERG(EIG)
C
      TID = SECOND() - TID
      IF(IPRSCF.GE.1) THEN
        CPUTID = SECTID(TDDG)
        IF(NZ.EQ.1) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Real diagonalization      (CPU):',CPUTID
        IF(NZ.EQ.2) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Complex diagonalization   (CPU):',CPUTID
        IF(NZ.EQ.4) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Quaternion diagonalization(CPU):',CPUTID
      ENDIF
      IF(ONESYS) THEN
        IF(PHCOEF) THEN
          DO I = 1,NFSYM
            CALL PHATRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ)
          ENDDO
          CALL WRICMO(LUCOEF,CMO,EIG,DHFERG)
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
        CLOSE(LUFCK1,STATUS='KEEP')
        CLOSE(LUFOCK,STATUS='KEEP')
        CPUTID = SECTID(TID)
        WRITE(LUPRI,'(A,A12)') 
     &    '* Total time for one-electron system:',CPUTID
      ELSE
        ITRSCF(INTFLG) = ITRSCF(INTFLG) + 1
        SCFTID(INTFLG) = SCFTID(INTFLG) + TID
        IF(IPRSCF.GE.1) THEN
          CPUTID = SECTID(TDF2)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* Two-electron Fock matrix  (CPU):',CPUTID
        ENDIF
        CPUTID = SECTID(TID)
        CALL GTINFO(DAYTID)
        WRITE(LUCYCL,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,INTTYP,DAYTID
        WRITE(LUPRI,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,INTTYP,DAYTID
      ENDIF
      CALL FLSHFO(LUPRI)
C
 666  CONTINUE
      RETURN
 1010 FORMAT('It. ',I3,1P,E20.10,3X,3E10.2,3X,A8,3X,A12,3X,A9,
     &       3X,A10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dfsolv */
      SUBROUTINE DFSOLV(FMO,EIG,FAO,DMAT,TMAT,BMAT,CHKCNV,DOSEL,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Solve Dirac-Fock eigenvalue equation
C     On input:
C       CMO - MO-coefficients from previous iteration
C       FAO - Fock-matrix in AO-basis
C     On output:
C       CMO - MO-coefficients from current iteration
C       EIG - eigenvalues
C
C     Input/output:
C       
C     Written by T.Saue Mar 09 1995
C     Last revision: Mar 09 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      LOGICAL CHKCNV,DOSEL
      DIMENSION FMO(*),EIG(*),FAO(*),DMAT(*),TMAT(*),BMAT(*),
     &          WORK(LWORK),IQBF(4)
C
C     Transform to MO-basis
C     =====================
C
      CALL DZERO(FMO,N2ORBTQ)
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NORB(I),NORB(I),
     &            FAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,0),
     &            FMO(I2ORBT(I)+1),NORB(I),NORB(I),NBF,IQBF,
     &            TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &            TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &            WORK,LWORK,IPRSCF)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** DFSOLV: NBF.NE.NZ !'
 10   CONTINUE
C
C     Convergence control/acceleration
C     ================================
C
      CALL DHFCNV(FMO,DMAT,FAO,TMAT,BMAT,CHKCNV,WORK,LWORK)
C
      TDDG = SECOND()
      DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
C
C       Diagonalize Fock-matrix in MO-basis
C         - store coefficients in FAO
C       ===================================
C
        IF(DOJACO) THEN
          CALL RSJACO(NORB(I),NORB(I),NORB(I),FMO(I2ORBT(I)+1),
     &                EIG(IORB(I)+1),1,1,0,FAO(I2ORBT(I)+1),
     &                WORK,LWORK)
        ELSE
          CALL QDIAG(NZ,NORB(I),FMO(I2ORBT(I)+1),NORB(I),NORB(I),
     &             EIG(IORB(I)+1),1,FAO(I2ORBT(I)+1),NORB(I),
     &             NORB(I),WORK,LWORK,IERR)
        ENDIF
        IF(IPRSCF.GE.2) CALL EIGCHK(EIG(IORB(I)+1),NPSH(I))
        IF(IPRSCF.GE.3) THEN
          CALL HEADER('ELECTRON STATES - eigenvalues',-1)
          IVAL = NPSH(I)
          CALL PREIGV(EIG(IORB(I)+1),NESH(I),IVAL,WORK,LWORK)
          CALL HEADER('POSITRON STATES - eigenvalues',-1)
          IVAL = 0
          CALL PREIGV(EIG(IORB(I)+1),NPSH(I),IVAL,WORK,LWORK)
        ENDIF
 20   CONTINUE
      TDDG = SECOND() - TDDG
C
C
C     Overlap selection
C
      IF(DOSEL) CALL SMOSEL(FAO,WORK,LWORK)
C
      DO 30 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 30
C
C       Backtransform all spinors
C         - temporary storage in FMO
C       ==============================
C
        CALL BCKTRA(FMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              FAO(I2ORBT(I)+1),NORB(I),NORB(I),
     &              NORB(I),NZ,IPQTOQ(1,0),NORB(I),1,NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,
     &              IPQTOQ(1,0),IPRSCF)
 30   CONTINUE
C
C     Write coefficients to file
C     ==========================
C
C     Write untransformed coefficients
      REWIND LUCMOS
      CALL WRITT(LUCMOS,N2ORBTQ,FAO)
C     Write transformed coefficients
      CALL WRICMO(LUCOEF,FMO,EIG,DHFERG)
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfscf */
      SUBROUTINE DHFSCF(CMO,EIG,FOCK,TMAT,WORK,LWORK)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <dcbgen.h>
#include <dcbham.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <cbihr1.h>
      DIMENSION FOCK(*),EIG(*),CMO(*),TMAT(*),
     &          WORK(LWORK)
C
C     Generate one-electron Fock matrix and write to file
C     ===================================================
C
      CALL ONEFCK(FOCK,IPRONE,WORK,LWORK)
      CALL OPNFIL(LUFCK1,'DFFCK1','UNKNOWN','ONEFC1')
      CALL WRIFCK(LUFCK1,FOCK)
      CLOSE(LUFCK1,STATUS='KEEP')
      IF(ONESYS) THEN
C
C     Solve one-electron system ...
C     =============================
C 
        CALL PREDHF(CMO,EIG,FOCK,TMAT,DUMMY,DUMMY,
     &              IDUMMY,WORK,LWORK)
      ELSE
C
C     ... or perform SCF-cycle
C     ========================
C 
        CALL SCFCYC(CMO,EIG,FOCK,TMAT,WORK,LWORK)
      ENDIF                                                        
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densty */
      SUBROUTINE DENSTY(DMAT,CMO,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),CMO(*)
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),
     &              CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &              NPSH(I)+1,NISH(I),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL TITLER('Output from DENSTY','*',103)
        DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck denst1 */
      SUBROUTINE DENST1(DMAT,LRD,LCD,CMO,LRC,LCC,JVEC,NVEC,NBAS)
C*****************************************************************************
C
C     Construct density matrix for fermion corep IFRP
C
C     Written by T.Saue Mar 3 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0,D0=0.0D0)
#include <dgroup.h>
C
      DIMENSION DMAT(LRD,LCD,NZ),CMO(LRC,LCC,NZ)
      CALL QGEMM(NBAS,NBAS,NVEC,D1,
     &     'N','N',IPQTOQ(1,0),CMO(1,JVEC,1),LRC,LCC,NZ,
     &     'H','N',IPQTOQ(1,0),CMO(1,JVEC,1),LRC,LCC,NZ,
     &          D0,IPQTOQ(1,0),DMAT,LRD,LCD,NZ)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck scfcyc */
      SUBROUTINE SCFCYC(CMO,EIG,FOCK,TMAT,WORK,LWORK)
C*****************************************************************************
C
C     Perform Dirac-Fock SCF-cycle
C
C     Written by T.Saue Mar 3 1995
C     Last revision: tsaue Mar 3 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
      DIMENSION FOCK(*),EIG(*),CMO(*),
     &          TMAT(*),WORK(LWORK)
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <cbihr2.h>
#include <blocks.h>
      CALL QENTER('SCFCYC')
#include <memint.h>
C
C     Memory allocation
C
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT,(MXDIIS*(MXDIIS+1))/2,
     &            WORK,KFREE,LFREE)
      NPOS = 0
      IF(PARCAL) THEN
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
      ENDIF
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
C     Initiate/restart SCF-cycle
C 
      CALL SCFCY1(CMO,EIG,FOCK,TMAT,WORK(KDMAT),WORK(KBMAT),
     &            WORK(KPOS),WORK(KFREE),LFREE)
      CALL MEMREL('SCFCYC.SCFCY1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT ('SCFCYC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck scfcy1 */
      SUBROUTINE SCFCY1(CMO,EIG,FOCK,TMAT,DMAT,BMAT,
     &                  NPOS,WORK,LWORK)
C*****************************************************************************
C
C     Perform Dirac-Fock SCF-cycle
C
C     Written by T.Saue Mar 3 1995
C     Last revision: tsaue Mar 3 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,DM1=-1.0D0,D1=1.0D0)
C
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbham.h>
#include <dcbint.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      CHARACTER CPUTID*12,SECTID*12,DAYTID*24
      DIMENSION FOCK(*),EIG(*),CMO(*),DMAT(*),
     &          BMAT(*),TMAT(*),NPOS(*),WORK(LWORK)
      CALL QENTER('SCFCY1')
#include <memint.h>
C
C
C     ***********************************
C     ***** Perform SCF-iterations ******
C     ***********************************
C
C     Initiate/restart SCF cycle
C
      CALL PREDHF(CMO,EIG,FOCK,TMAT,DMAT,BMAT,
     &            NPOS,WORK,LWORK)
C
C     Start cycle
C
      DO 10 ITER = 2,MAXITR
        NITER  = NITER + 1
        KITER  = KITER + 1        
        CALL GTINFO(DAYTID)
        CALL GETTIM(CTIDSTR,WTIDSTR)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +    '########## START ITERATION NO.',NITER,'##########',
     &    DAYTID
C
C	Check what integral types contribute
C       ====================================
C
        CALL INTCON(INTFLG,INTBUF,INTDEF,
     &              ERGVAL,CNVINT,NITER,ITRINT,INTTYP)
C
C       Construct density matrix
C       ========================
C
        CALL DENSTY(DMAT,CMO,IPRSCF)
C
C       Construct two-electron Fock matrix
C       (stored in FOCK)
C       ==================================
C
        TDF2 = SECOND()
        IF(DODSCF.AND.INTFLG.EQ.INTBUF) THEN
          CALL READNS(LUDENS,CMO)
          IF(FIXDIF) THEN
            DCOVLP = D1
          ELSE
            DCOVLP = DDOT(N2BBASXQ,CMO,1,CMO,1)
            DCOVLP = DDOT(N2BBASXQ,DMAT,1,CMO,1) / DCOVLP
          ENDIF
          WRITE(LUPRI,'(I5,A,F6.4)') NITER,
     &       ' *** Differential density matrix. DCOVLP =',
     &       DCOVLP
          CALL DSCAL(N2BBASXQ,(-DCOVLP),CMO,1)
          CALL DAXPY(N2BBASXQ,D1,DMAT,1,CMO,1)
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,CMO,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK,LWORK)
          CALL REAFCK(LUFCK2,CMO,.FALSE.)
          CALL DAXPY(N2BBASXQ,DCOVLP,CMO,1,FOCK,1)
        ELSE
          CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FOCK,DMAT,NFMAT,
     &                NPOS,INTFLG,IPRSCF,WORK,LWORK)
        ENDIF
        IF(DOTSC) CALL TSCORR(FOCK,WORK,LWORK,IPRSCF)
        TDF2 = SECOND() - TDF2
C
C       Write current two-electron Fock matrix to file
C       ==============================================
C
        CALL WRIDNS(LUDENS,DMAT)
        CALL WRIFCK(LUFCK2,FOCK)
C
C       Read one-electron Fock matrix from file
C       =======================================
C
        CALL REAFCK(LUFCK1,CMO,.FALSE.)
C
C       Calculate total energy
C       ======================
C
        CALL ERGCAL(CMO,FOCK,DMAT)
C
C       Add one-electron Fock matrix to two-electron matrix
C       ===================================================
C
        CALL DAXPY(N2BBASXQ,D1,CMO,1,FOCK,1)
        IF(IPRSCF.GE.5) THEN
          CALL HEADER('Total Fock matrix',-1)
          CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C
C       Solve DF eigenvalue problem
C       ===========================
C       Input:  FOCK - Total Fock matrix on AO-basis
C               DMAT - Density matrix
C               TMAT - MO transformation matrix
C               BMAT - DIIS matrix
C       Output: CMO  - MO-coefficients for current iteration
C
        CALL DFSOLV(CMO,EIG,FOCK,DMAT,TMAT,BMAT,.TRUE.,
     &              OVLSEL,WORK,LWORK)
 20     CONTINUE
        IF(IPRSCF.GE.1) THEN
          CPUTID = SECTID(TDF2)
          WRITE(LUPRI,'(3X,A,A12)')
     &      '* Two-electron Fock matrix  (CPU):',CPUTID
          CPUTID = SECTID(TDDG)
          IF(NZ.EQ.1) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Real diagonalization      (CPU):',CPUTID
          IF(NZ.EQ.2) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Complex diagonalization   (CPU):',CPUTID
          IF(NZ.EQ.4) WRITE(LUPRI,'(3X,A,A12)')
     &      '* Quaternion diagonalization(CPU):',CPUTID
        ENDIF
        CALL GETTIM(CTIDEND,WTIDEND)
        CTIDTOT = CTIDEND - CTIDSTR
        WTIDTOT = WTIDEND - WTIDSTR
        ITRSCF(INTFLG) = ITRSCF(INTFLG) + 1
        SCFTID(INTFLG) = SCFTID(INTFLG) + WTIDTOT
        CPUTID = SECTID(CTIDTOT)
        WRITE(LUPRI,'(A,A12)') 
     &    '>>> Total CPU-time :',CPUTID
        CPUTID = SECTID(WTIDTOT)
        WRITE(LUPRI,'(A,A12)') 
     &    '>>> Total wall time:',CPUTID
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)')
     +    '########## END ITERATION NO.',NITER,'##########',
     &    DAYTID
        WRITE(LUCYCL,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,INTTYP,DAYTID(1:10)
        WRITE(LUPRI,1010) NITER,DHFERG,ERGVAL,FCKVAL,EVCVAL,
     &                CACC,CPUTID,INTTYP,DAYTID(1:10)
        CALL FLSHFO(LUPRI)
        IF(EXIT) GOTO 30
        INTBUF = INTFLG
 10   CONTINUE
C
C     **********************************
C     ***** End of SCF-iterations ******
C     **********************************
C
 30   CONTINUE
C
C     Phase-adjustment of coefficients, if requested...
C
      IF(PHCOEF) THEN
        DO I = 1,NFSYM
          CALL PHATRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ)
        ENDDO
        CALL WRICMO(LUCOEF,CMO,EIG,DHFERG)
      ENDIF
C
C     Close necessary files
C
      CLOSE(LUCOEF,STATUS='KEEP')
      CLOSE(LUCMOS,STATUS='KEEP')
      CLOSE(LUDENS,STATUS='KEEP')
      CLOSE(LUFOCK,STATUS='KEEP')
      CLOSE(LUFCK1,STATUS='KEEP')
      CLOSE(LUFCK2,STATUS='KEEP')
      CLOSE(LUCYCL,STATUS='KEEP')
      CLOSE(LUEVEC,STATUS='KEEP')
C
      CALL QEXIT('SCFCY1')
      RETURN
 1010 FORMAT('It. ',I3,1P,E20.10,3X,3E10.2,3X,A8,3X,A12,3X,A9,
     &       3X,A10)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ergcal */
      SUBROUTINE ERGCAL(F1,F2,DMAT)
C****************************************************************************
C     PURPOSE:
C        Calculate Dirac-Hartree-Fock energy
C
C	E = D(IJ)(*)(F1(IJ)+0.5F2(IJ))
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D2=2.0D0)
#include <mxcent.h>
C
#include <nuclei.h>
#include <frame.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbdhf.h>
#include <dcbpsi.h>
      LOGICAL LVCON
      CHARACTER M(4)*2,FMT*6,MXFORM*6
      DIMENSION F1(N2BBASX,NZ),F2(N2BBASX,NZ),DMAT(N2BBASX,NZ),
     &          E(4,2)
      DATA M/'AR','AI','BR','BI'/
C
      CALL QENTER('ERGCAL')
C
      LVCON = DOLVC.AND.(INTFLG.EQ.3)
      CALL DZERO(E,8)
      ERGBUF = DHFERG
      E1PART = D0
      E2PART = D0
      DO 10 IZ = 1,NZ
        E(IZ,1) = DDOT(N2BBASX,DMAT(1,IZ),1,F1(1,IZ),1)
        E1PART  = E1PART + E(IZ,1)
        E(IZ,2) = DDOT(N2BBASX,DMAT(1,IZ),1,F2(1,IZ),1)
        E2PART  = E2PART + E(IZ,2)
   10 CONTINUE
      E1PART = D2*E1PART
      ELERGY = E1PART + E2PART
      DHFERG = ELERGY + POTNUC
      IF(LVCON) DHFERG = DHFERG + CORRLV
C
      IF(IPRSCF.GE.1) THEN
        Call HEADER('Output from ERGCAL',-1)
        ERGMAX = MAX(ABS(E1PART),ABS(E2PART))*10
        FMT = MXFORM(ERGMAX,20)
        IF(LVCON) THEN
          WRITE(LUPRI,'(A,'//FMT//')')
     +       'Electronic energy         :  ',ELERGY,
     +       'Nuclear potential energy  :  ',POTNUC,
     +       'SS Coulombic correction   :  ',CORRLV,
     +       'Total energy              :  ',DHFERG
        ELSE
          WRITE(LUPRI,'(A,'//FMT//')')
     +       'Electronic energy         :  ',ELERGY,
     +       'Nuclear potential energy  :  ',POTNUC,
     +       'Total energy              :  ',DHFERG
        ENDIF
        IF(IPRSCF.GE.2) THEN
          WRITE(LUPRI,'(6X,A,'//FMT//')')
     +     '* 1-electron energy :  ',E1PART
          IF(NZ.GT.1) THEN
            DO 20 IZ = 1,NZ
              WRITE(LUPRI,'(6X,A,A2,A,'//FMT//')')
     +        '            from ',M(IZ),' :  ',D2*E(IZ,1)
   20       CONTINUE
          ENDIF
          WRITE(LUPRI,'(6X,A,'//FMT//')')
     +     '* 2-electron energy :  ',E2PART
          IF(NZ.GT.1) THEN
            DO 30 IZ = 1,NZ
              WRITE(LUPRI,'(6X,A,A2,A,'//FMT//')')
     +        '            from ',M(IZ),' :  ',E(IZ,2)
   30       CONTINUE
          ENDIF
        ENDIF
      ENDIF
C
      CALL QEXIT('ERGCAL')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intcon */
      SUBROUTINE INTCON(INTFLG,INTOLD,INTDEF,
     &                  CNV,CNVINT,ITR,ITRINT,INTTYP)
C*****************************************************************************
C
C     Check what integrals contribute in a given SCF-iteration
C       CNVINT has precedence over ITRINT !!!
C       The information is bitpacked in INTDEF:
C         001   LL-contribution
C         010   SL-contribution
C         100   SS-contribution
C       and all combinations.....
C
C     Written by T.Saue Mar 09 1995
C     Last revision: Aug 24 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      LOGICAL LBIT
      CHARACTER INTTYP*9
      DIMENSION CNVINT(2),ITRINT(2)
      INTTYP = '         '
C
C     LL-integrals
C     ============
C
      IF(LBIT(INTDEF,1)) THEN
        ILLINT = 1
      ELSE
        ILLINT = 0
      ENDIF
C
C     SL-integrals
C     ============
C
      IF(LBIT(INTDEF,2)) THEN
        ISLINT = 1
        IF(.NOT.LBIT(INTOLD,2)) THEN
          IF(ABS(CNV).GT.CNVINT(1)) ISLINT = 0
          IF(ITR     .LT.ITRINT(1)) ISLINT = 0
        ENDIF
      ELSE
        ISLINT = 0
      ENDIF      
C
C     SS-integrals
C     ============
C
      IF(LBIT(INTDEF,3)) THEN
        ISSINT = 1
        IF(.NOT.LBIT(INTOLD,3)) THEN
          IF(ABS(CNV).GT.CNVINT(2)) ISSINT = 0
          IF(ITR     .LT.ITRINT(2)) ISSINT = 0
        ENDIF
      ELSE
        ISSINT = 0
      ENDIF      
C
      IF(ILLINT.EQ.1) INTTYP(1:3) = 'LL '
      IF(ISLINT.EQ.1) INTTYP(4:6) = 'SL '
      IF(ISSINT.EQ.1) INTTYP(7:9) = 'SS '
      INTFLG= ILLINT + 2*ISLINT + 4*ISSINT
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smogen */
      SUBROUTINE SMOGEN(CMO,TMAT,WORK,LWORK)
C*****************************************************************************
C
C     Multiply overlap matrix onto occupied vectors and write to file;
C     this forms the basis for overlap selection
C
C     Written by T.Saue Nov 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C Used from COMMON blocks:
C  DCBBAS: N2BBASX,NFBAS
C  DCBORB: NOCC
C  DHFCAL: NSMOTQ
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION CMO(*),TMAT(*),WORK(LWORK)
C
      CALL QENTER('SMOGEN')
#include <memint.h>
      NBFDIM = 0
      DO I = 1,NFSYM
        NBFDIM = MAX(NBFDIM,NOCC(I)*NFBAS(I,0))
      ENDDO
      NBFDIM = NBFDIM*NZ
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMO,NSMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF,NBFDIM,WORK,KFREE,LFREE)
      CALL SMOGE1(WORK(KSMO),CMO,TMAT,WORK(KSMAT),WORK(KBUF))
      CALL MEMREL('SMOGEN.SMOGE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SMOGEN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smoge1 */
      SUBROUTINE SMOGE1(SMO,CMO,TMAT,SMAT,BUF)
C*****************************************************************************
C
C     Multiply overlap matrix onto occupied vectors and write to file;
C     this forms the basis for overlap selection
C
C     Overlap selection is based on a previous set CPRE of vectors in AO-basis
C     We form the matrix
C          SMO  =  T^(+) * S * CPRE
C     where S is the overlap matrix and T the MO-transformation matrix.
C     In each SCF-iteration we form the matrix
C          SEL  = CMO^(+) * SMO
C     using the current CMO-coefficients in MO-basis.
C
C     Written by T.Saue Nov 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
C Used from COMMON blocks:
C  PAMIOU: LUOVLP
C  DCBBAS:
C  DCBORB:
C  DGROUP: NFSYM
C
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbdhf.h>
      DIMENSION SMO(*),CMO(*),TMAT(*),SMAT(*),BUF(*)
C
C     Get overlap matrix
C
      CALL OPNFIL(LUOVLP,'DFOVLP','OLD','SMOGE1')
      CALL READT(LUOVLP,N2BBASX,SMAT)
      CLOSE(LUOVLP,STATUS='KEEP')
      DO I = 1,NFSYM
      IF(NOCC(I).NE.0) THEN
        IOFF = ICMOQ(I) + NPSH(I)*NFBAS(I,0)
        CALL QGEMM(NFBAS(I,0),NOCC(I),NFBAS(I,0),D1,
     &             'N','N',IPQTOQ(1,0),SMAT(I2BASX(I,I)+1),
     &                     NTBAS(0),NTBAS(0),1,
     &             'N','N',IPQTOQ(1,0),CMO(IOFF+1),
     &                     NFBAS(I,0),NORB(I),NZ,
     &        D0,IPQTOQ(1,0),BUF,NFBAS(I,0),NOCC(I),NZ)
        CALL QGEMM(NORB(I),NOCC(I),NFBAS(I,0),D1,
     &             'H','N',IPQTOQ(1,0),TMAT(I2TMT(I)+1),
     &                     NFBAS(I,0),NORB(I),NZT,
     &             'N','N',IPQTOQ(1,0),BUF,NFBAS(I,0),NOCC(I),NZ,
     &        D0,IPQTOQ(1,0),SMO(ISMOQ(I)+1),NORB(I),NOCC(I),NZ)
        IF(IPRSCF.GE.4) THEN
          CALL HEADER('SMOGE1: SMO-matrix',-1)
          CALL PRQMAT(SMO(ISMOQ(I)+1),NORB(I),NOCC(I),NORB(I),NOCC(I),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
      ENDIF
      ENDDO
C
C     Write to file
C
      CALL OPNFIL(LUSMOS,'DFSMOS','UNKNOWN','SMOGE1')
      CALL WRITT(LUSMOS,NSMOTQ,SMO)
      CLOSE(LUSMOS,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smosel */
      SUBROUTINE SMOSEL(CFU,WORK,LWORK)
C*****************************************************************************
C
C     Select vectors based on overlap.
C     Input is coefficients in MO-basis.
C
C     Written by T.Saue Nov 11 1995
C     Last revision Nov 11 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION CFU(*),WORK(LWORK)
C Used from COMMON blocks:
C  DHFCAL: NSMOTQ
C
#include <dcbdhf.h>
      CALL QENTER('SMOSEL')
#include <memint.h>
C     Memory allocation
      CALL MEMGET('REAL',KSMO,NSMOTQ,WORK,KFREE,LFREE)
      CALL SMOSE1(CFU,WORK(KSMO),WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('SMOSEL.SMOSE1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('SMOSEL')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smose1 */
      SUBROUTINE SMOSE1(CFU,SMO,WORK,LWORK)
C*****************************************************************************
C
C     Select vectors based on overlap
C
C     Written by T.Saue Nov 11 1995
C     Last revision Nov 11 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION CFU(*),SMO(*),WORK(LWORK)
C Used from COMMON blocks:
C  DCBORB: NSMOTQ,NOCC
C  DGROUP: NFSYM
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
      CALL QENTER('SMOSE1')
#include <memint.h>
C
C     Read vectors to select from
C
      CALL OPNFIL(LUSMOS,'DFSMOS','OLD','SMOSE1')
      CALL READT(LUSMOS,NSMOTQ,SMO)
      CLOSE(LUSMOS,STATUS='KEEP')
C
      DO I = 1,NFSYM
      IF(NOCC(I).NE.0) THEN
C       Memory allocation
        CALL MEMGET('REAL',KSBUF,NESH(I)*NOCC(I)*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KCBUF,NESH(I)*NORB(I)*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KCSEL,NESH(I),WORK,KFREE,LFREE)
        CALL SMOSE2(I,CFU(I2ORBT(I)+1),SMO(ISMOQ(I)+1),WORK(KSBUF),
     &              WORK(KCBUF),WORK(KCSEL),ISELM,
     ^              NORB(IFRP),NOCC(IFRP),NESH(IFRP),
     &              WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('SMOSE1.SMOSE2',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
      ENDDO
      CALL QEXIT('SMOSE1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck smose2 */
      SUBROUTINE SMOSE2(IFRP,CFU,SMO,SBUF,CBUF,ICSEL,ISELM,
     &                  KNORB,KNOCC,KNESH,WORK,LWORK)
C*****************************************************************************
C
C     Select vectors based on overlap
C
C     Written by T.Saue Nov 11 1995
C     Last revision: Jan 7 1998 - jth NORB(IFRP) etc. is now parameter
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
C Used from COMMON blocks:
C
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
#include <dcbdhf.h>
C
C     KNORB = NORB(IFRP)
C     KNOCC = NOCC(IFRP)
C     KNESH = NESH(IFRP)
C
      CHARACTER FSYM*2
      DIMENSION CFU(KNORB,KNORB,NZ),
     &          SMO(KNORB,KNOCC,NZ),
     &          SBUF(KNESH,KNOCC,NZ),
     &          CBUF(KNORB,KNESH,NZ),
     &          ICSEL(KNESH),WORK(LWORK)
C
C     Form product CFU * (SMAT * COCC)
C
      CALL QGEMM(NESH(IFRP),NOCC(IFRP),NORB(IFRP),D1,
     &       'H','N',IPQTOQ(1,0),CFU(1,NPSH(IFRP)+1,1),
     &                NORB(IFRP),NORB(IFRP),NZ,
     &       'N','N',IPQTOQ(1,0),SMO,NORB(IFRP),NOCC(IFRP),NZ,
     &       D0,IPQTOQ(1,0),SBUF,NESH(IFRP),NOCC(IFRP),NZ)
      IF(IPRSCF.GE.4) THEN
        CALL HEADER('SMOSE1: SBUF-matrix',-1)
        CALL PRQMAT(SBUF,NESH(IFRP),NOCC(IFRP),NESH(IFRP),NOCC(IFRP),
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      IF(NZ.EQ.1) THEN
        DO J = 1,NOCC(IFRP)
          DO I = 1,NESH(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
          ENDDO
        ENDDO
      ELSEIF(NZ.EQ.2) THEN
        DO J = 1,NOCC(IFRP)
          DO I = 1,NESH(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
     &                  + SBUF(I,J,2)*SBUF(I,J,2)
          ENDDO
        ENDDO
      ELSEIF(NZ.EQ.4) THEN
        DO J = 1,NOCC(IFRP)
          DO I = 1,NESH(IFRP)
            SBUF(I,J,1) = SBUF(I,J,1)*SBUF(I,J,1)
     &                  + SBUF(I,J,2)*SBUF(I,J,2)
     &                  + SBUF(I,J,3)*SBUF(I,J,3)
     &                  + SBUF(I,J,4)*SBUF(I,J,4)
          ENDDO
        ENDDO
      ENDIF
      IF(IPRSCF.GE.4) THEN
        CALL HEADER('SMOSE1: Squared-matrix',-1)
        CALL PRQMAT(SBUF,NESH(IFRP),NOCC(IFRP),NESH(IFRP),NOCC(IFRP),
     &              1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Perform selection
C
      N2DIM = NORB(IFRP)*NESH(IFRP)
      DO IZ = 1,NZ
        CALL DCOPY(N2DIM,CFU(1,NPSH(IFRP)+1,IZ),1,CBUF(1,1,IZ),1)
      ENDDO
      ISELM = 0
      DO J = 1,NOCC(IFRP)
        JJ = NPSH(IFRP) + J
        ISEL = IDAMAX(NESH(IFRP),SBUF(1,J,1),1)
        ICSEL(J) = ISEL
        ISELM = MAX(ISELM,ISEL)
C       Zero out the selected overlap
        CALL DCOPY(NOCC(IFRP),D0,0,SBUF(ISEL,1,1),NESH(IFRP))
C       Copy buffer into vector
        DO IZ = 1,NZ
          CALL DCOPY(NORB(IFRP),CBUF(1,ISEL,IZ),1,CFU(1,JJ,IZ),1)
        ENDDO
      ENDDO
C
C     Copy in missing vectors
C
      IF(ISELM.GT.NOCC(IFRP)) THEN
        JSTR = NOCC(IFRP)+1
        DO J = JSTR,ISELM
          JJ = NPSH(IFRP) + J
          ISEL = IDAMAX(ISELM,SBUF(1,J,1),1)
          ICSEL(J) = ISEL
C         Zero out the selected overlap
          CALL DCOPY(NOCC(IFRP),D0,0,SBUF(ISEL,1,1),NESH(IFRP))
C         Copy buffer into vector
          DO IZ = 1,NZ
            CALL DCOPY(NORB(IFRP),CBUF(1,ISEL,IZ),1,CFU(1,JJ,IZ),1)
          ENDDO
        ENDDO
      ENDIF
C
C     Print section
C
      WRITE(FSYM,'(I2)') IFRP
      CALL HEADER('SMOSEL: Reordered vectors in fermion ircop '
     &            //FSYM,-1)
      NVEC = 0
      DO J = 1,NOCC(IFRP)
      IF(ICSEL(J).NE.J) THEN
        NVEC = NVEC + 1
        WRITE(LUPRI,'(I5,A,I5)') J,' --> ',ICSEL(J)
      ENDIF
      ENDDO
      IF(NVEC.EQ.0) WRITE(LUPRI,'(A)') 'No reordering'
      IF(IPRSCF.GE.2) THEN
        WRITE(LUPRI,'(20I6)') (ICSEL(I),I=1,NOCC(IFRP))
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck denorb */
      SUBROUTINE DENORB(DMAT,JVEC,NVEC,IFRP,CMO,IPRINT)
C*****************************************************************************
C
C     Construct density matrix for NVEC orbitals of fermion corep IFRP
C     starting with orbital JVEC (with respect to NPSH(IFRP))
C
C     Written by T.Saue Oct 22 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0,D0=0.0D0)
      DIMENSION DMAT(*),CMO(*)
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
C
      JOFF = NPSH(IFRP) + JVEC 
      CALL DZERO(DMAT,N2BBASXQ)
      CALL DENST1(DMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),
     &              CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &              JOFF,NVEC,NFBAS(IFRP,0))
C
C     Print section
C
      IF(IPRINT.GE.4) THEN
        CALL TITLER('Output from DENORB','*',103)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',IFRP
        CALL PRQMAT(DMAT(I2BASX(IFRP,IFRP)+1),NFBAS(IFRP,0),
     &              NFBAS(IFRP,0),NTBAS(0),NTBAS(0),
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck oneerg */
      SUBROUTINE ONEERG(EIG)
C****************************************************************************
C     PURPOSE:
C        Calculate one-electron energy as sum of eigenvalues
C
C        Written by T.Saue Jan 28 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D2=2.0D0)
C
#include <mxcent.h>
#include <frame.h>
#include <nuclei.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbdhf.h>
      CHARACTER FMT*6,MXFORM*6
      DIMENSION EIG(*)
C
      CALL QENTER('ONEERG')
C
      ERGBUF = DHFERG
      E1PART = D0
      E2PART = D0
      DO IFRP = 1,NFSYM
        IF(NOCC(IFRP).GT.0) THEN
          E1PART = E1PART + D2*DSUM(NOCC(IFRP),
     &                       EIG(IORB(IFRP)+NPSH(IFRP)+1),1)
        ENDIF
      ENDDO
      ELERGY = E1PART
      DHFERG = ELERGY + POTNUC
C
      IF(IPRSCF.GE.1) THEN
        Call HEADER('Output from ONEERG',-1)
        ERGMAX = ABS(E1PART)*10
        FMT = MXFORM(ERGMAX,20)
        WRITE(LUPRI,'(A,'//FMT//')')
     +     'Electronic energy         :  ',ELERGY,
     +     'Nuclear potential energy  :  ',POTNUC,
     +     'Total energy              :  ',DHFERG
      ENDIF
C
      CALL QEXIT('ONEERG')
      RETURN
C
      END
