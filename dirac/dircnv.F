C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dircnv.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfcnv */
      SUBROUTINE DHFCNV(FMO,FAO,DMAT,EVEC,FBUF,TMAT,BMAT,
     &                  CHKCNV,SMOAO,WORK,LWORK)
C*****************************************************************************
C     Driver routine for checking convergence etc.
C
C     Written by T.Saue October 1995
C     Last revision 971028-jth: Lot of changes to make one open shell work.
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
      DIMENSION FMO(*), EVEC(*),FBUF(*),BMAT(*),TMAT(*),
     &          DMAT(*),SMOAO(*),WORK(LWORK), FAO(*)
C
C     Calculate error vector for this iteration
C     =========================================
C             FMO   - Fock matrix in MO-basis 
C             FAO   - Fock matrix in AO-basis (FDao and FVao)
C             DMAT  - Density matrix in AO-basis (DCao and DVao)
C             TMAT  - MO-transformation matrix
C             SMOAO - Left index transformed overlap matrix
C             FBUF  - work area
C
      IF(CHKCNV) CALL DIISER(EVEC,FMO,FAO,DMAT,SMOAO,
     &                       FBUF,TMAT,WORK,LWORK)
C
C
C     Read MO Fock matrix from previous cycle into FBUF...
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

C     ...and write current MO Fock matrix to file
C     ==============================================
C
      IF(MXDIIS.EQ.0) THEN
        IREC = 1
        CALL WRTDAC(LUFOCK,N2ORBTQ,FMO,IREC)
      ELSE
        IREC = MOD(NITER-1,MXDIIS-1) + 1
        CALL WRTDAC(LUFOCK,N2ORBTQ,FMO,IREC)
      ENDIF
C
C     If no check on convergence, then return here...
C
      IF(.NOT.CHKCNV) RETURN
C
C     Calculate largest change in AO Fock matrix for this iteration
C     =============================================================
C
      CALL FCKCMP(FMO,FBUF,WORK,LWORK)
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
C     I (jth) want DIIS from iteration no. 2 and no damping.
C     hjj has done non-relativistic tests of DIIS, and it was
C     never found that DIIS slowed things down; in the worst
C     cases there was the same number of iterations.
C
      DIISON = DODIIS.AND.ABS(EVCVAL).LE.DIISTH
C
C     If change in which integral classes included (LL, LS, SS),
C     then we must restart DIIS
C
      IF (INTFLG.NE.INTBUF) THEN
        IF(ITDIIS.GT.0) CLOSE(LUDIIS,STATUS='DELETE')
        ITDIIS = 0
      END IF
C
C     DIIS --- Direct Inversion of Iterative Subspace
C     ===============================================
C
      IF(DIISON) THEN
        ITDIIS = ITDIIS + 1
        CALL DIISDR(BMAT,EVEC,FMO,FBUF,IBSTAT,WORK,LWORK)
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
             DO I = 1,NFSYM
                WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
                CALL HEADER('DHFCNV: MO Fock matrix after DIIS:',-1)
                CALL PRQMAT(FMO(I2ORBT(I)+1),NORB(I),NORB(I),
     &                   NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
             END DO
          END IF 
        ENDIF
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
        CALL DSCAL(N2ORBTQ,DNEW  ,FMO,1)
        CALL DAXPY(N2ORBTQ,DAMPFC,FBUF,1,FMO,1)
        IF(IPRSCF.GE.5) THEN
          DO I = 1,NFSYM
             WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
             CALL HEADER('DHFCNV: MO Fock matrix after DAMP:',-1)
             CALL PRQMAT(FMO(I2ORBT(I)+1),NORB(I),NORB(I),
     &                   NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          END DO
        ENDIF
      ENDIF
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diiser */
      SUBROUTINE DIISER(EVEC,FMO,FAO,DMAT,SMOAO,EBUF,TMAT,WORK,LWORK)
C****************************************************************************
C
C     Last revision 971029-jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0,D1=1.0D0,D2 = 2.0D0, DM1 = -1.00D00, 
     &          D4 = 4.00D00, DP5 = 0.50D00)
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
C     Note that EBUF has dimension NFMAT*N2BBASXQ
C 
      DIMENSION FAO(*),DMAT(*),EVEC(*),EBUF(*),TMAT(*),
     &          SMOAO(*),WORK(LWORK), FMO(*)
      CHARACTER*1 CTMP
C
#include <memint.h>
C
C
      CALL QENTER('DIISER')
C
C
      IF (DIISAO.OR.(NASHT.GT.0)) THEN
C
         CALL DZERO(EBUF,NFMAT*N2BBASXQ)
         CALL DZERO(EVEC,N2EVECTQ)
C
         IOFF = 0
         IOFF2 = 1
         EVCVAL = D0
C
C        The gradient is
C           g = DCao * FDao + f_o DVao * (FDao + (a_o-1) QVao) + H.C.
C             = (DCao+f_o DVao) FDao + f_o (a_o-1) DVao QVao + H.C.
C
C   
         IF (NASHT.GT.0) THEN
C
C           Construct total density matrix
C
            DO IOPEN = 1,NOPEN
               CALL DAXPY(N2BBASXQ,DF(IOPEN),DMAT(1+N2BBASXQ*IOPEN),1,
     &                                       DMAT,1)
            END DO
         END IF
         DO 10 I = 1,NFSYM
            IF (NOCC(I).EQ.0) GOTO 10
C
C              (DCao+f_o DVao) FDao
C
            CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),D1,
     &                 'N','N',IPQTOQ(1,0),DMAT(1+I2BASX(I,I)),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 'N','N',IPQTOQ(1,0),FAO(1+I2BASX(I,I)),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 D0,IPQTOQ(1,0),EBUF(1+I2BASX(I,I)),
     &                         NTBAS(0),NTBAS(0),NZ)
            IF ( IPRSCF .GE. 10 ) THEN
               WRITE(CTMP,'(I1)') I
               CALL HEADER('DIISER: (DCao+DVao) FDao, corep '//CTMP,-1)
               CALL PRQMAT(EBUF(1+I2BASX(I,I)),
     &                     NFBAS(I,0),NFBAS(I,0),
     &                     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IF (NASHT.GT.0) THEN
               DO 15 IOPEN = 1,NOPEN
                  IF (NACSH(I,IOPEN).EQ.0) GOTO 15
C
C                 (a_o-1) f_o DVao QVao
C                 = (a_o-1) f_o DVao f FVao
C                 where FVao is the Fock matrix from TWOINT
C
                  DSCALE = DF(IOPEN)*(DA(IOPEN)-D1)*DF(IOPEN)
                  CALL QGEMM(NFBAS(I,0),NFBAS(I,0),NFBAS(I,0),DSCALE,
     &                 'N','N',IPQTOQ(1,0),
     &                         DMAT(1+I2BASX(I,I)+N2BBASXQ*IOPEN),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 'N','N',IPQTOQ(1,0),
     &                         FAO(1+I2BASX(I,I)+N2BBASXQ*IOPEN),
     &                         NTBAS(0),NTBAS(0),NZ,
     &                 D1,IPQTOQ(1,0),EBUF(1+I2BASX(I,I)),
     &                         NTBAS(0),NTBAS(0),NZ)
                  IF ( IPRSCF .GE. 10 ) THEN
                     WRITE(CTMP,'(I1)') I
                     CALL HEADER('DIISER: (DCao+f_o DVao) FDao + '//
     &                           'f_o (a_o-1) DVao QVao, corep '//
     &                           CTMP,-1)
                     CALL PRQMAT(EBUF(1+I2BASX(I,I)),
     &                           NFBAS(I,0),NFBAS(I,0),
     &                           NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
                  END IF
 15            CONTINUE
            END IF
C
C           Transform to MO-basis.: Ct Sao A C
C           store in FAO
C
            CALL QGEMM(NORB(I),NFBAS(I,0),NFBAS(I,0),D1,
     &                 'N','N',IPQTOQ(1,0),SMOAO(IOFF2),
     &                   NORB(I),NFBAS(I,0),NZ,
     &                 'N','N',IPQTOQ(1,0),EBUF(1+I2BASX(I,I)),
     &                   NTBAS(0),NTBAS(0),NZ,
     &                 D0,IPQTOQ(1,0),FAO,
     &                   NORB(I),NFBAS(I,0),NZ)
C
            CALL QGEMM(NORB(I),NORB(I),NFBAS(I,0),D1,
     &                 'N','N',IPQTOQ(1,0),FAO,
     &                    NORB(I),NFBAS(I,0),NZ,
     &                 'N','N',IPQTOQ(1,0),TMAT(I2TMT(I)+1),
     &                    NFBAS(I,0),NORB(I),NZT,
     &                 D0,IPQTOQ(1,0),EVEC(1+IOFF),
     &                    NORB(I),NORB(I),NZ)
C
            IF ( IPRSCF .GE. 10 ) THEN
               WRITE(CTMP,'(I1)') I
               CALL HEADER('DIISER: Error vector in MO-basis, corep '
     &                     //CTMP,-1)
               CALL PRQMAT(EVEC(1+IOFF),NORB(I),NORB(I),
     &                     NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
C
C           Calculate ErrVec - ErrVec(Hermitian conj.)
C
            CALL DIFMMH(EVEC(IOFF+1),NORB(I),NZ,NORB(I),NORB(I))
C
            IF ( IPRSCF .GE. 10 ) THEN
               WRITE(CTMP,'(I1)') I
               CALL HEADER('DIISER: ErrVec - ErrVec(Hermitian conj.)'//
     &                     ', corep '//CTMP,-1)
               CALL PRQMAT(EVEC(1+IOFF),NORB(I),NORB(I),
     &                  NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
C
            IOFF2 = IOFF2 + NZ*NFBAS(I,0)*NORB(I)
            IOFF = IOFF + NORB(I)*NORB(I)*NZ
 10      CONTINUE     
C
C        The error vector is now 0.25 as in SIRIUS, so multiply with 4
C        But in SIRIUS the norm is taken of the triangular matrix
C           so we divide by SQRT(2). 
C
         EVCVAL = D4*DNRM2(N2EVECTQ,EVEC,1)/SQRT(D2)
C
C
      ELSE
C
C        Old DIIS.
C           Only for closed shell molecules.
C
C
C        Construct density matrix in MO-basis
C        ====================================
C
         REWIND LUCMOS
         CALL READT(LUCMOS,N2ORBTQ,EVEC)
         DO 20 I = 1,NFSYM
            IF(NOCC(I).EQ.0) GOTO 20
            IOFF = I2ORBT(I)+NORB(I)*NPSH(I)+1
            CALL QGEMM(NORB(I),NORB(I),NISH(I),D1,
     &           'N','N',IPQTOQ(1,0),EVEC(IOFF),NORB(I),NORB(I),NZ,
     &           'H','N',IPQTOQ(1,0),EVEC(IOFF),NORB(I),NORB(I),NZ,
     &             D0,IPQTOQ(1,0),EBUF(I2ORBT(I)+1),NORB(I),NORB(I),NZ)
            IF(IPRSCF.GE.5) THEN
               WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
               CALL HEADER('DIISER:Density matrix in MO-basis',-1)
               CALL PRQMAT(EBUF(I2ORBT(I)+1),NORB(I),NORB(I),
     &                   NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
 20      CONTINUE
C
C        Form error vector and find max element EVCVAL
C        =============================================
C
         IOFF = 0
         EVCVAL = D0
         DO 30 I = 1,NFSYM
            IF(NOCC(I).EQ.0) GOTO 30
            CALL QGEMM(NORB(I),NORB(I),NORB(I),D1,
     &         'N','N',IPQTOQ(1,0),FMO(I2ORBT(I)+1),NORB(I),NORB(I),NZ,
     &         'N','N',IPQTOQ(1,0),EBUF(I2ORBT(I)+1),NORB(I),NORB(I),NZ,
     &           D0,IPQTOQ(1,0),EVEC(IOFF+1),NORB(I),NORB(I),NZ)
            IF(IPRSCF.GE.5) THEN
               CALL HEADER('DIISER:FD matrix',-1)
               WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
               CALL PRQMAT(EVEC(IOFF+1),NORB(I),NORB(I),NORB(I),NORB(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
            CALL DIFMMH(EVEC(IOFF+1),NORB(I),NZ,NORB(I),NORB(I))
            IF(IPRSCF.GE.5) THEN
               CALL HEADER('DIISER:Error matrix',-1)
               WRITE(LUPRI,'(3X,A,I2)') '*** Fermion ircop ',I
               CALL PRQMAT(EVEC(IOFF+1),NORB(I),NORB(I),NORB(I),NORB(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
            IOFF = IOFF + NORB(I)*NORB(I)*NZ
   30    CONTINUE
C
C        Find total max element in EVEC
C        ===============================
C
         EVCVAL = D2*DNRM2(N2EVECTQ,EVEC,1)
C
      END IF
C
      CALL FLSHFO(LUPRI)
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
C   DCIDHF: NFMAT
C   DGROUP: NZ
C
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbbas.h>
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
          WRITE(LUPRI,'(5X,A6,4X,1P,D10.3)')
     &       'lambda',CVEC(1)
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
      PARAMETER(DM1 = -1.0D0, D0 = 0.00D00)
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
      PARAMETER(DP5=0.5D00,D8=8.00D00)
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
COld Code
C
      IF (DIISMO) THEN
         BIJELM = DDOT(N2EVECTQ,EVEC1,1,EVEC2,1)
         BIJELM = BIJELM + BIJELM
      ELSE
C   
C     jth: The new errorvector is half as big as in DALTON
C       => the trace of two matrices is one fourth as big
C       => we multiply by four (old code multiplied by two, so
C                               a total of 8)
C
         BIJELM = D8*DDOT(N2EVECTQ,EVEC1,1,EVEC2,1)
      END IF

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
C  DCIDHF: NFMAT
C  DGROUP: NZ
C
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbbas.h>
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
