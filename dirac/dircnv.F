C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfcnv */
      SUBROUTINE DHFCNV(FOCK,EVEC,FBUF,TMAT,BMAT,CHKCNV,WORK,LWORK)
C*****************************************************************************
C     Driver routine for checking convergence etc.
C
C     Written by T.Saue October 1995
C     Last revision March 12 1996  - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
C Used from COMMON blocks
C   DCBBAS: N2BBASX
C   DCBORB: N2ORBT
C   DGROUP: NZ
C
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      LOGICAl CHKCNV
      DIMENSION FOCK(*),EVEC(*),FBUF(*),BMAT(*),TMAT(*),WORK(LWORK)
C
C     Calculate error vector for this iteration
C     =========================================
C             FOCK  - Fock matrix in MO-basis
C             TMAT - MO-transformation matrix
C             FBUF  - work area
C
      IF(CHKCNV) CALL DIISER(EVEC,FOCK,FBUF,TMAT,WORK,LWORK)
C
C
C     Read Fock matrix from previous cycle into FBUF...
C     =================================================
C
      IF(CHKCNV) THEN
        IF(MXDIIS.EQ.0) THEN
          IREC = 1
        ELSE
          IREC = MOD(NITER-2,MXDIIS-1) + 1
        ENDIF
        CALL READAC(LUFOCK,N2ORBTQ,FBUF,IREC)
      ENDIF

C     ...and write current Fock matrix to file
C     ==============================================
C
      IF(MXDIIS.EQ.0) THEN
        IREC = 1
        CALL WRTDAC(LUFOCK,N2ORBTQ,FOCK,IREC)
      ELSE
        IREC = MOD(NITER-1,MXDIIS-1) + 1
        CALL WRTDAC(LUFOCK,N2ORBTQ,FOCK,IREC)
      ENDIF
C
C     If no check on convergence, then return here...
C
      IF(.NOT.CHKCNV) RETURN
C
C     Calculate largest change in Fock matrix for this iteration
C     ==========================================================
C
      CALL FCKCMP(FOCK,FBUF,WORK,LWORK)
C
C     Check convergence
C     =================
C
      ERGVAL  = ERGBUF - DHFERG
      IF    (ERGCNV) THEN
        CONVRG = ERGVAL
      ELSEIF(FCKCNV) THEN
        CONVRG = FCKVAL
      ELSEIF(EVCCNV) THEN
        CONVRG = EVCVAL
      ENDIF
      IF(ABS(CONVRG).LT.SCFCNV) SUCCES = .TRUE.
      IF(SUCCES.OR.(KITER.GT.MAXITR)) EXIT = .TRUE.
C
C     Convergence acceleration
C     ========================
C
      CACC = '        '
      DIISON = DODIIS.AND.ABS(EVCVAL).LE.DIISTH.AND.INTFLG.EQ.INTBUF
C
C     DIIS --- Direct Inversion of Iterative Subspace
C     ===============================================
C
      IF(DIISON) THEN
        ITDIIS = ITDIIS + 1
        CALL DIISDR(BMAT,EVEC,FOCK,FBUF,IBSTAT,WORK,LWORK)
        IF(IBSTAT.EQ.1) THEN
          DIISON = .FALSE.
C
C         If DODAMP Fock matrix from previous cycle has to be reread from
C         previous cycle after being overwritten in DIISDR
C
          IF(DODAMP) THEN
            IREC = MOD(NITER-2,MXDIIS-1) + 1
            CALL READAC(LUFOCK,N2ORBTQ,FBUF,IREC)
          ENDIF
        ELSE
          CACC(1:4) = 'DIIS'
          WRITE(CACC(5:8),'(I4)') NELMBM
          IF(IPRSCF.GE.5) THEN
            CALL HEADER('DHFCNV: Fock matrix after DIIS:',-1)
            CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
        ENDIF
      ELSE
        IF(ITDIIS.GT.0) CLOSE(LUDIIS,STATUS='DELETE')
        ITDIIS = 0
      ENDIF
C
C     Damping of total Fock matrix
C     ============================
C
      DAMPON = DODAMP.AND.(.NOT.DIISON).AND.INTFLG.EQ.INTBUF
      IF(DAMPON) THEN
        CACC(1:4) = 'DAMP'
        WRITE(CACC(5:8),'(F4.2)') DAMPFC
        DNEW = D1 - DAMPFC
        CALL DSCAL(N2ORBTQ,DNEW  ,FOCK,1)
        CALL DAXPY(N2ORBTQ,DAMPFC,FBUF,1,FOCK,1)
        IF(IPRSCF.GE.5) THEN
          DO I = 1,NFSYM
            WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
            CALL HEADER('DHFCNV: Fock matrix after DAMP:',-1)
            CALL PRQMAT(FOCK(I2ORBT(I)+1),NORB(I),NORB(I),
     &                  NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diiser */
      SUBROUTINE DIISER(EVEC,FOCK,EBUF,TMAT,WORK,LWORK)
C****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0,D1=1.0D0,D2 = 2.0D0)
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbdhf.h>
C
C Used from COMMON blocks:
C    DCBBAS: N2BBASX,NTBAS(0)
C    DCBORB: NFORB,
C    DGROUP: NZ,NFSYM
C
      DIMENSION FOCK(*),EVEC(*),EBUF(*),TMAT(*),EMAX(2),
     &          WORK(LWORK)
C
      CALL QENTER('DIISER')
C
C       Construct density matrix in MO-basis
C       ====================================
C
      REWIND LUCMOS
      CALL READT(LUCMOS,N2ORBTQ,EVEC)
      DO 10 I = 1,NFSYM
        IF(NOCC(I).EQ.0) GOTO 10
        IOFF = I2ORBT(I)+NORB(I)*NPSH(I)+1
        CALL QGEMM(NORB(I),NORB(I),NISH(I),D1,
     &     'N','N',IPQTOQ(1,0),EVEC(IOFF),NORB(I),NORB(I),NZ,
     &     'H','N',IPQTOQ(1,0),EVEC(IOFF),NORB(I),NORB(I),NZ,
     &          D0,IPQTOQ(1,0),EBUF(I2ORBT(I)+1),NORB(I),NORB(I),NZ)
      IF(IPRSCF.GE.5) THEN
        WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
        CALL HEADER('DIISER:Density matrix in MO-basis',-1)
        CALL PRQMAT(EBUF(I2ORBT(I)+1),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
 10   CONTINUE
C
C     Form error vector and find max element EVCVAL
C     =============================================
C
      IOFF = 0
      EMAX(1) = D0
      EMAX(2) = D0
      EVCVAL = D0
      DO 20 I = 1,NFSYM
        IF(NOCC(I).EQ.0) GOTO 20
        CALL QGEMM(NORB(I),NORB(I),NORB(I),D1,
     &      'N','N',IPQTOQ(1,0),FOCK(I2ORBT(I)+1),NORB(I),NORB(I),NZ,
     &      'N','N',IPQTOQ(1,0),EBUF(I2ORBT(I)+1),NORB(I),NORB(I),NZ,
     &           D0,IPQTOQ(1,0),EVEC(IOFF+1),NORB(I),NORB(I),NZ)
        IF(IPRSCF.GE.5) THEN
          CALL HEADER('DIISER:FD matrix',-1)
          WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
          CALL PRQMAT(EVEC(IOFF+1),NORB(I),NORB(I),NORB(I),NORB(I),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        CALL DIFMMH(EVEC(IOFF+1),NORB(I),NZ,NORB(I),NORB(I))
        IF(IPRSCF.GE.5) THEN
          CALL HEADER('DIISER:Error matrix',-1)
          WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
          CALL PRQMAT(EVEC(IOFF+1),NORB(I),NORB(I),NORB(I),NORB(I),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IOFF = IOFF + NORB(I)*NORB(I)*NZ
   20 CONTINUE
C
C     Find total max element in EVEC
C     ===============================
C
      EVCVAL = D2*DNRM2(N2EVECTQ,EVEC,1)
C
      CALL QEXIT('DIISER')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckcmp */
      SUBROUTINE FCKCMP(FNEW,FOLD,WORK,LWORK)
C*****************************************************************************
C
C     Find largest change in Fock matrix between two iterations
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C Used from COMMON blocks:
C   DCBORB: N2ORBTQ
C   DGROUP: NZ
C
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION FNEW(*),FOLD(*),WORK(LWORK)
C
      CALL QENTER('FCKCMP')
#include <memint.h>
C
C     Memory allocation
      CALL MEMGET('REAL',KDIFF,N2ORBTQ,WORK,KFREE,LFREE)
C
      FCKVAL =  DIFFMX(FNEW,FOLD,WORK(KDIFF),N2ORBTQ)
C
C     Memory deallocation
      CALL MEMREL('FCKCMP.DIFFMX',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT ('FCKCMP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diffmx */
      FUNCTION DIFFMX(AVEC,BVEC,DIFF,NELM)
C*****************************************************************************
C
C     Find maximum differnce between two vectors
C
C     Written by T.Saue Aug 25 1995
C     Last revision Aug 25 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION AVEC(NELM),BVEC(NELM),DIFF(NELM)
C
      DO I = 1,NELM
        DIFF(I) = AVEC(I) - BVEC(I)
      ENDDO
      IMAX = IDAMAX(NELM,DIFF,1)
      DIFFMX = DIFF(IMAX)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diisdr */
      SUBROUTINE DIISDR(BMAT,EVEC,FOCK,EBUF,IBSTAT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for DIIS - Direct Inversion of Iterative Subspaces
C
C     Written by Jon Laerdahl and T. Saue 1995
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C Used from COMMON blocks:
C   DHFCAL: MXDIIS
C
#include <dcbdhf.h>
      DIMENSION BMAT(*),EVEC(*),FOCK(*),EBUF(*),WORK(LWORK)
C
      CALL QENTER('DIISDR')
#include <memint.h>
      CALL MEMGET('REAL',KCVEC,MXDIIS,WORK,KFREE,LFREE)
      CALL DIISD1(BMAT,EVEC,FOCK,EBUF,WORK(KCVEC),IBSTAT,
     &            WORK(KFREE),LFREE)
      CALL MEMREL('DIISDR.DIISD1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('DIISDR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diisd1 */
      SUBROUTINE DIISD1(BMAT,EVEC,FOCK,EBUF,CVEC,IBSTAT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for DIIS - Direct Inversion of Iterative Subspaces
C
C     Written by T.Saue Aug 28 1995
C     Last revision Aug 28 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DM1 = -1.0D0,D0 = 0.0D0)
C
C Used from COMMON blocks:
C   DHFCAL: LUEVEC,IPRSCF,ITDIIS,MXDIIS
C
#include <dcbdhf.h>
      DIMENSION BMAT(MXDIIS*MXDIIS),EVEC(*),FOCK(*),EBUF(*),CVEC(*),
     &          WORK(LWORK)
      IPRINT = IPRSCF
C
C     Write current error vector to direct access unit
C     ================================================
C
      IREC = MOD(ITDIIS-1,MXDIIS-1) + 1
      CALL WRTDAC(LUEVEC,N2EVECTQ,EVEC,IREC)
C
C     First DIIS iteration; initialize
C     ================================
C
      IF(ITDIIS.EQ.1) THEN
         NELMBM = 1
        BMAT(1) = D0
        II = 1
        DO I = 1,MXDIIS-1
          II = II + I
          BMAT(II) = DM1
        ENDDO
      ENDIF
      IF(ITDIIS.EQ.2) CALL OPNFIL(LUDIIS,'DFDIIS','UNKNOWN','DIISDR')
C
C     Solve DIIS equation for current iteration
C     ========================================
C
      CALL EQDIIS(BMAT,CVEC,EVEC,EBUF,IBSTAT,IPRINT,
     &            WORK,LWORK)
C
C     Construct total Fock matrix using coefficients of CVEC
C
      IF(IBSTAT.EQ.0) THEN
        CALL FMDIIS(FOCK,EBUF,CVEC,IPRINT)
        REWIND LUDIIS
        WRITE (LUDIIS) ITDIIS,MXDIIS,NELMBM,BMAT
      ELSEIF(IBSTAT.EQ.1) THEN
        CLOSE(LUDIIS,STATUS='DELETE')
      ENDIF
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Output from DIISD1:',-1)
        WRITE(LUPRI,'(A,I10)') '* DIIS iteration no.:',ITDIIS
        WRITE(LUPRI,'(A,I10)') '* B matrix dimension:',NELMBM
        WRITE(LUPRI,'(A,1P,D10.3)') '* B matrix condition:',BMCOND
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(A)') '* Current B matrix:'
          CALL PRMUTC(BMAT,NELMBM,1,LUPRI)
        ENDIF
        IF(IBSTAT.EQ.0) THEN
          CFSUM = D0
          DO I = 2,NELMBM
            CFSUM = CFSUM + CVEC(I)
          ENDDO
          WRITE(LUPRI,'(/A,1P,D9.3/)') '* DIIS coefficients.Sum: ',CFSUM
          IOFF = NITER - NELMBM
          WRITE(LUPRI,'(5X,I5,5X,1P,D10.3)')
     &       (IOFF+I,CVEC(I),I=2,NELMBM)
        ENDIF
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eqdiis */
      SUBROUTINE EQDIIS(BMAT,CVEC,EVEC,EBUF,IBSTAT,IPRINT,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Solve DIIS equations
C
C     Written by Jon Laerdahl & Trond Saue
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C Used from COMMON blocks:
C   DHFCAL: MXDIIS
#include <dcbdhf.h>
      DIMENSION BMAT(*),CVEC(*),EVEC(*),EBUF(*),WORK(LWORK)
C
      CALL QENTER('EQDIIS')
#include <memint.h>
C
C     Memory allocation
      CALL MEMGET('REAL',KBFAC,(MXDIIS*(MXDIIS+1))/2,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KKPVT,MXDIIS,WORK,KFREE,LFREE)
      CALL EQDII1(BMAT,WORK(KBFAC),EVEC,EBUF,CVEC,IBSTAT,
     &            WORK(KKPVT),IPRINT,WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('EQDIIS.EQDII1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('EQDIIS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eqdii1 */
      SUBROUTINE EQDII1(BMAT,BFAC,EVEC,EBUF,CVEC,IBSTAT,KPVT,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Solve DIIS equations
C
C     Written by Jon Laerdahl & Trond Saue
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DM1 = -1.0D0)
      DIMENSION BMAT(*),BFAC(*),EVEC(*),EBUF(*),CVEC(*),
     &          KPVT(*),WORK(LWORK)
C Used from COMMON blocks
C  DHFCAL: NELMBM
C
#include <dcbdhf.h>
C
C     Construct and factorize B - matrix
C     ==================================
C
      CALL BMDIIS(BMAT,BFAC,EVEC,EBUF,CVEC,KPVT,
     &            IBSTAT,IPRINT)
C
C     Solve DIIS equation
C     ===================
C
      IF(IBSTAT.EQ.0) THEN
        CALL DZERO(CVEC,NELMBM)
        CVEC(1) = DM1
        CALL DSPSL(BFAC,NELMBM,KPVT,CVEC)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bmdiis */
      SUBROUTINE BMDIIS(BMAT,BFAC,EVEC,EBUF,CVEC,KPVT,IBSTAT,
     &                  IPRINT)
C*****************************************************************************
C
C     Calculate and factor DIIS B matrix.
C     If the matrix is ill-conditioned
C     Elements are of form B_ij = <e_i|e_j> where e_i,e_j are error vectors.
C     The B matrix is symmetric and is therefore upper triangular
C     columnwise packed.
C
C     Written by Jon Laerdahl & T.Saue
C     Last revision: Aug 27 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
C Used from COMMON blocks:
C  DHFCAL: LUEVEC,NELMBM,ITDIIS,MXDIIS
C  DGROUP: NZ
C  DCBBAS: N2VCTQ
C
#include <dcbdhf.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
      DIMENSION BMAT(*),BFAC(*),EVEC(*),EBUF(*),CVEC(*),KPVT(*)
      CALL QENTER('BMDIIS')
C
      NELMBM  = NELMBM + 1
C
C     If necessary, make room for more elements in B matrix
C
      IF(NELMBM.GT.MXDIIS) THEN
        NELMBM = NELMBM - 1
        CALL BSHIFT(BMAT,NELMBM)
      ENDIF
C
C     Add new column to B matrix; IOFF is offset to
C     first element of new column
C
      IOFF  = (NELMBM*(NELMBM-1))/2 + 1
      NPREV = NELMBM - 2
      JREC  = ITDIIS - NPREV - 2
      DO I = 1, NPREV
        IREC = MOD(JREC+I,MXDIIS-1) + 1
        CALL READAC(LUEVEC,N2EVECTQ,EBUF,IREC)
        BMAT(IOFF+I) = BIJELM(EBUF,EVEC)
      ENDDO
      BMAT(IOFF+NPREV+1) = BIJELM(EVEC,EVEC)
C
C     Factor and check condition of B matrix
C     ======================================
C
   10 CONTINUE
      IF(NELMBM.LE.2) THEN
        IBSTAT = 1
        GOTO 999
      ELSE
        NXBM = (NELMBM*(NELMBM+1))/2
        CALL DCOPY(NXBM,BMAT,1,BFAC,1)
        CALL DSPCO(BFAC,NELMBM,KPVT,BMCOND,CVEC)
        IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(A,I5)') 'ITDIIS: ',ITDIIS,
     & 'NELMBM: ',NELMBM
          WRITE(LUPRI,'(A,1P,E9.3)') 'BMCOND : ',BMCOND
        ENDIF
C
C       Ill-conditioned matrix; delete oldest elements
C
        IF (D1+BMCOND.EQ.D1) THEN
          IF(IPRINT.GE.3) THEN
            WRITE(LUPRI,'(A,1P,E9.3)') '*** WARNING *** B-cond:',BMCOND
            WRITE(LUPRI,'(/A/)') 'B matrix:'
            CALL PRMUTC(BMAT,NELMBM,1,LUPRI)
          ENDIF
          CALL BSHIFT(BMAT,NELMBM)
          NELMBM = NELMBM - 1
          GOTO 10
        ENDIF
        IBSTAT = 0
      ENDIF
C
 999  CONTINUE
      CALL QEXIT ('BMDIIS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bijelm */
      FUNCTION BIJELM(EVEC1,EVEC2)
C*****************************************************************************
C
C     Calculate an element of the DIIS B matrix
C
C     Written by T.Saue Aug 28 1995
C     Last revision: Oct 18 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0)
C
C Used from COMMON blocks:
C   DHFCAL: N2EVECQ,I2EVECTQ
C   DGROUP: NFSYM
C   DCBORB: NOCC
C
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION EVEC1(*),EVEC2(*)
C
      BIJELM = DDOT(N2EVECTQ,EVEC1,1,EVEC2,1)
      BIJELM = BIJELM + BIJELM
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bshift */
      SUBROUTINE BSHIFT(BMAT,NELMBM)
C*****************************************************************************
C
C     Delete row 2 and column 2 of B matrix and shift all the
C     other elements accordingly diagonally upwards, that is
C       B(i,j) --> B(i-1,j-1)
C     The B matrix is upper triangularly columnwise packed, so
C     the transformation is
C       B([ij]) --> B([ij] - j)
C     The formula for a general shift
C       B(i,j) --> B(i-n,j-n)
C     is
C       B([ij]) --> B([ij]-N), N = [nj - n(n-1)/2]
C
C     Written by T.Saue Aug 28 1995
C     Last revision Aug 28 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION BMAT(*)
C
C     J runs over columns, I over rows;
C     ID gives the supervectorindex of diagonal element of previous
C     column
C
      ID = 3
      DO J = 3,NELMBM
        DO I = 3,NELMBM
          IJ = ID + I
          BMAT(IJ-J) = BMAT(IJ)
        ENDDO
        ID = ID + J
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fmdiis */
      SUBROUTINE FMDIIS(FOCK,FBUF,CVEC,IPRINT)
C*****************************************************************************
C
C     Construct total Fock matrix using coefficients obtained
C     by solving DIIS equations
C
C     Written by T.Saue Aug 28 1995
C     Last revision : Aug 28 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C Used from COMMON blocks:
C  DHFCAL: LUFOCK,NITER,NELMBM
C  DCBBAS: N2BBASX
C  DGROUP: NZ
C
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION FOCK(*),FBUF(*),CVEC(*)
C
      CALL DSCAL(N2ORBTQ,CVEC(NELMBM),FOCK,1)
      IOFF  = 1
      NPREV = NELMBM - 2
      JREC  = NITER - NPREV - 2
      DO I = 1,NPREV
        IREC = MOD(JREC+I,MXDIIS-1) + 1
        CALL READAC(LUFOCK,N2ORBTQ,FBUF,IREC)
        CALL DAXPY(N2ORBTQ,CVEC(IOFF+I),FBUF,1,FOCK,1)
      ENDDO
C
      RETURN
      END
