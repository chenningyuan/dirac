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
C
C
C VERSION : $Revision: 1.36.4.3 $
C DATE    : $Date: 2004/09/01 08:16:06 $
C FILE    : pamexp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpexp */
      SUBROUTINE PRPEXP(WF,EXPVAL,WORK,LWORK)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*),WORK(LWORK)
      CHARACTER WF*4
C
#include "dcbexp.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('PRPEXP')
#include "memint.h"
C
c     CALL MEMCHK('begin prpexp',WORK,1)
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KEXPST ,2*NEXPP   ,WORK,KFREE,LFREE)
CHJ TODO : check if 2 DMAT or NOPEN+1 DMAT is the correct ...
      CALL MEMGET('REAL',KDMAT  ,NZ*N2BBASX*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT  ,N2BBASX   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTRI  ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST ,NZ        ,WORK,KFREE,LFREE)
C
      CALL PRPEX1(WF,EXPVAL,WORK(KEXPST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),WORK(KPTRI),
     &            WORK(KFIRST),WORK(KFREE),LFREE)
C
C
C     Memory deallocation
C
      CALL MEMREL('PRPEXP',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PRPEXP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex1 */
      SUBROUTINE PRPEX1(WF,EXPVAL,LEXPST,CMO,EIG,DMAT,PMAT,PTRI,
     &                  FIRST,WORK,LWORK)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0,DTOL=1.0D-16)
C
      LOGICAL FIRST(NZ),LEXPST(2,*),FILEX
      DIMENSION EXPVAL(*),CMO(*),DMAT(*),PMAT(*),PTRI(NNBBASX,NZ),
     &          EIG(*),WORK(LWORK),NP(2,0:2)
      CHARACTER MXFORM*6,FMT*6
      CHARACTER WF*4
      INTEGER   GETIOPEN
      LOGICAL FNDLAB
C
#include "cbihr1.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
#include "dcbprj.h"
C
C
C     Initialize
C
#include "memint.h"
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Initialize
C     ==========
C
      CALL DZERO(EXPVAL,NEXPP)
      N2EXP = 2*NEXPP
      CALL LSET(N2EXP,.FALSE.,LEXPST)
C
C     Generate density matrix
C     =======================
C
      IF (WF .EQ. 'DHF' .OR. WF .EQ. 'DFT') THEN
         INQUIRE(FILE='DFCOEF',EXIST=FILEX)
         IF(.NOT.FILEX) GOTO 1000
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP')
         CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6,
     &               WORK,KFREE,LFREE)
         CLOSE(LUCOEF,STATUS='KEEP')
         IF (.NOT. ORBANA) THEN
            CALL GENDEN(DMAT,CMO,1,IPREXP)
            DO IOPEN = 1, NOPEN
               CALL ACTDEN(DMAT(1+N2BBASXQ),CMO,IOPEN,IPREXP)
               CALL DAXPY(N2BBASXQ,D1,DMAT(1+N2BBASXQ),1,DMAT,1)
            END DO
         END IF
#ifdef MOD_KRMC
      ELSE IF (WF .EQ. 'KRMC') THEN
C
         IF (ORBANA) 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'ORBANA not implemented for MCSCF')
C
         INQUIRE(FILE='KRMCSCF', EXIST=FILEX)
C
         IF (.NOT. FILEX) 
     &      CALL QUIT('*** ERROR in PRPEX1 *** ' //
     &        'file "KRMCSCF" not found!')
C
         CALL OPNFIL(LUCOEF,'KRMCSCF','OLD','PRPEXP')
         CALL REAKRMC(LUCOEF,'NEWORB  ',CMO,NCMOTQ)
         CLOSE(LUCOEF,STATUS='KEEP')
C
         INQUIRE(FILE='KRMC_FOCK', EXIST=FILEX)
C
         IF (.NOT. FILEX) 
     &      CALL QUIT('*** ERROR in PRPEX1 *** ' //
     &        'file "KRMC_FOCK" not found!')
C
         CALL OPNFIL(LUCOEF,'KRMC_FOCK','OLD','PRPEXP')
         CALL REAKRMC(LUCOEF,'DTMO    ',DMAT(1+N2BBASXQ),N2ORBXQ)
         CLOSE(LUCOEF,STATUS='KEEP')
C
C        Get AO density matrix.
C
         DO I = 1, NFSYM
            CALL QTRANS('MOAO','S',D0,
     &           NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &           DMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DMAT(1+N2BBASXQ+I2ORBX(I,I)),
     &           NORBT,NORBT,NZ,IPQTOQ(1,0),
     &           CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK,LWORK,0)
         END DO
C
C        Result is scaled by 2 in PRPEX2...
C
         CALL DSCAL(N2BBASXQ,DP5,DMAT,1)
C
#endif
      ELSE
         WRITE(LUPRI,'(1X,2A)')
     &        '*** ERROR in PRPEX1 *** Unknown wave function: ',WF
         CALL QUIT('*** ERROR in PRPEX1 *** Unknown wave function')
      END IF
C
C     Projection analysis
C     ===================
C
      IF(PRJANA) THEN
        INQUIRE(FILE='DFPRJC',EXIST=FILEX)
        IF(.NOT.FILEX) THEN
          WRITE(LUPRI,'(A)')
     &      '** WARNING ! No projection coefficients found !',
     &      '    - Skipping expectation value projection analysis !'
          PRJANA = .FALSE.
          KRVEC = KFREE
          KBVEC = KFREE
          KPRMT = KFREE
        ELSE
          CALL OPNFIL(LUCOEF,'DFPRJC','OLD','PRPEXP')
          CALL MEMGET('INTE',KNR,NFSYM*NREFS,WORK,KFREE,LFREE)
          CALL REAPRJ(LUCOEF,NP,WORK(KNR),KRVEC,KBVEC,
     &                WORK,KFREE,LFREE)
          CLOSE(LUCOEF,STATUS='KEEP')
          NBMAX = 0
          NPFMT = 0
          DO IFRP = 1,NFSYM
            NBMAX = MAX(NBMAX,NP(IFRP,0))
            NPFMT = NPFMT + NP(IFRP,0)*NP(IFRP,0)*NZ
          ENDDO
          NN = NBMAX*(NBMAX+1)/2
          CALL MEMGET('REAL',KPFMT,NPFMT,WORK,KFREE,LFREE)       
          CALL MEMGET('REAL',KDFMT,NPFMT,WORK,KFREE,LFREE)       
          CALL MEMGET('INTE',KIBUF,NN,WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KINDX,NBMAX*2,WORK,KFREE,LFREE)
          NTOT = NREFS + 1
          NFDIM = NTOT*NTOT
          CALL MEMGET('REAL',KFMAT,NFDIM,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KPM  ,NTOT,WORK,KFREE,LFREE)
C.........Calculate density matrix in fragment basis
          CALL PRJEX4(WORK(KDFMT),WORK(KBVEC),NP(1,0),NP(1,2),IPREXP)
        ENDIF
      ELSE
        KRVEC = KFREE
        KBVEC = KFREE
        KPFMT = KFREE
      ENDIF
C
C
C     Expectation values normal procedure
C     ===================================
C
      IF(.NOT.ORBANA) THEN
        CALL TITLER('Expectation values','*',127)
        WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
        CALL PRSYMB(LUPRI,'-',76,0)
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
C.........generate property matrix in SO-basis
          CALL PRPEX2(LU1INT,INDXPR,EXPVAL(I),LEXPST(1,I),
     &                IREP,ITIM,PTRI,PMAT,FIRST,
     &                WORK(KFREE),LFREE,IPREXP)
          IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 30
C.............calculate expectation value
          CALL PRPEX3(EXPVAL(I),IREP,ITIM,
     &        DMAT,PMAT,PTRI,FIRST,ORBANA,
     &        WORK(KFREE),LFREE,IPREXP)
          IF(PRJANA) THEN
C.........Transform property matrix to fragment orbital basis
            CALL PRJEX1(WORK(KPFMT),PTRI,PMAT,WORK(KRVEC),
     &                  INDXPR,IREP,ITIM,NP(1,0),FIRST,IPREXP,
     &                  WORK,KFREE,LFREE)
          ENDIF
 30       CONTINUE
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
     &       PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.',
     &       LEXPST(1,I),LEXPST(2,I)
          CALL FLSHFO(LUPRI)
          IF(PRJANA.AND.ABS(EXPVAL(I)).GT.DTOL) THEN
            CALL PRJEX2(EXPVAL(I),PTRI,WORK(KPFMT),WORK(KDFMT),
     &                  WORK(KBVEC),
     &                  WORK(KIBUF),INDXPR,NP(1,0),NP(1,2),
     &                  WORK(KNR),NTOT,WORK(KFMAT),WORK(KPM),
     &                  WORK(KINDX),
     &                  IPREXP,PMAT,WORK,KFREE,LFREE)
            CALL PRSYMB(LUPRI,'-',76,0)
          ENDIF
        ENDDO
        CALL PRSYMB(LUPRI,'-',76,0)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
C
C     Expectation values by orbital
C    ************************************
C
      ELSE
        CALL TITLER('Expectation values','*',127)
        WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
        CALL PRSYMB(LUPRI,'-',76,0)
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          CALL PRSYMB(LUPRI,'-',60,4)
          CALL HEADER('Operator '//PRPNAM(INDXPR)//':',0)            
          CALL HEADER('Expectation value for individual orbitals',-1)
          EXPVAL(I) = D0
          DO 10 IFRP = 1,NFSYM
            IF(NORB(IFRP).EQ.0) GOTO 10
C.........generate property matrix in SO-basis
            CALL PRPEX2(LU1INT,INDXPR,EXPBUF,LEXPST(1,I),
     &                  IREP,ITIM,PTRI,PMAT,FIRST,
     &                  WORK(KFREE),LFREE,IPREXP)
            IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 20
            DO J = 1,NOCC(IFRP)
              CALL DENORB(DMAT,J,1,IFRP,CMO,IPREXP)
              EXPBUF = D0
C.............calculate expectation value
              CALL PRPEX3(EXPBUF,IREP,ITIM,
     &          DMAT,PMAT,PTRI,FIRST,ORBANA,
     &          WORK(KFREE),LFREE,IPREXP)
              IOPEN = GETIOPEN(J,IFRP)
              EXPBUF = EXPBUF*DF(IOPEN)
              EXPVAL(I) = EXPVAL(I) + EXPBUF
              FMT = MXFORM(EXPBUF,15)
              WRITE(LUPRI,'(4X,A3,I5,3X,'//FMT//',1X,A4)')
     &          FREP(IFRP),J,EXPBUF,'a.u.'
              CALL FLSHFO(LUPRI)
            ENDDO
            IF(PRJANA) THEN
C...........Transform property matrix to fragment orbital basis
            CALL PRJEX1(WORK(KPFMT),PTRI,PMAT,WORK(KRVEC),
     &                  INDXPR,IREP,ITIM,NP(1,0),FIRST,IPREXP,
     &                  WORK,KFREE,LFREE)
            ENDIF
 10       CONTINUE
 20       CONTINUE
          CALL PRSYMB(LUPRI,'-',60,4)
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A8,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
     &     'Total : ',EXPVAL(I),'a.u.',
     &     LEXPST(1,I),LEXPST(2,I)
          CALL PRSYMB(LUPRI,'-',60,4)
          IF(PRJANA.AND.ABS(EXPVAL(I)).GT.DTOL) THEN
            CALL PRJEX2(EXPVAL(I),PTRI,WORK(KPFMT),WORK(KDFMT),
     &                  WORK(KBVEC),
     &                  WORK(KIBUF),INDXPR,NP(1,0),NP(1,2),
     &                  WORK(KNR),NTOT,WORK(KFMAT),WORK(KPM),
     &                  WORK(KINDX),IPREXP,PMAT,WORK,KFREE,LFREE)
            CALL PRSYMB(LUPRI,'-',76,0)
          ENDIF
        ENDDO
C        CALL PRSYMB(LUPRI,'-',76,0)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
      ENDIF
C
      CLOSE(LU1INT,STATUS = 'KEEP')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPEXP: Coefficient file not found !'
      CALL QUIT('PRPEXP: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex2 */
      SUBROUTINE PRPEX2(LU1INT,INDXPR,EXPVAL,LEXPST,IREP,ITIM,
     &                  PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     For totally symmetric operators return property matrix;
C     for other symmetries return .FALSE. in
C         LEXPST(1) - molecular point group symmetry
C     or  LEXPST(2) - time reversal symmetry
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      LOGICAL FIRST(NZ),LEXPST(2)
      DIMENSION PTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
      EXPVAL = D0
C
C     Check symmetry
C     ==============
C
C     1. Molecular point group symmetry
C
      IREP = IPRPSYM(INDXPR)-1
      IF(IREP.NE.0) LEXPST(1) = .TRUE.
C
C     2. Time reversal symmetry
C
      ITIM = IPRPTIM(INDXPR)
      IF(ITIM.NE.1) LEXPST(2) = .TRUE.
C
C     If not totally symmetric, return with zero value
C
      IF (LEXPST(1).OR.LEXPST(2)) RETURN
C
C     Totally symmetric operator ==> proceed
C
      CALL LSET(NZ,.TRUE.,FIRST)
      CALL PRPMAO(LU1INT,INDXPR,PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C
C     
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex3 */
      SUBROUTINE PRPEX3(EXPVAL,IREP,ITIM,
     &           DMAT,PMAT,PTRI,FIRST,ORBFLG,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      LOGICAL FIRST(NZ),LEXPST(2),ORBFLG
      DIMENSION DMAT(N2BBASX,NZ,*),PMAT(N2BBASX),
     &          PTRI(NNBBASX,NZ),WORK(LWORK)
C
C     Get expectation value, using only the matrices containing
C     integrals
C
      CALL DZERO(PMAT,N2BBASX)
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        IQ = IPQTOQ(IZ,IREP)
        IH = IHQMAT(IQ,ITIM)
        IF(IH.EQ.1) THEN
          CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT)
        ELSEIF(IH.EQ.2) THEN
          CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT)
        ENDIF
        CALL BUTOBS(PMAT,1,WORK,LWORK)
        IF (.NOT. ORBFLG) THEN
           DO IOPEN = 0,NOPEN
              EXPVAL = EXPVAL + DF(IOPEN)
     &           *DDOT(N2BBASX,DMAT(1,IZ,IOPEN+1),1,PMAT,1)
            END DO
         ELSE
            EXPVAL = EXPVAL + DDOT(N2BBASX,DMAT(1,IZ,1),1,PMAT,1)
         END IF
      END IF
      END DO
      EXPVAL = EXPVAL + EXPVAL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjex1 */
      SUBROUTINE PRJEX1(PFO,PTRI,PAO,RVEC,INDXPR,IREP,ITIM,NCNT,
     &                  FIRST,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Transform property matrix to fragment orbital basis
C
C     Written by Trond Saue Apr 27 2003 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION PFO(*),PAO(N2BBASX),NCNT(*),
     &          RVEC(*),PTRI(NNBBASX,NZ),WORK(*)
C
C     Transform on the fly
C
      IPF = 1
      DO IFRP = 1,NFSYM
        NDIM = NCNT(IFRP)*NCNT(IFRP)*NZ
        CALL DZERO(PFO(IPF),NDIM)
        IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
      ENDDO
      FADD = D0
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        IQ = IPQTOQ(IZ,IREP)
        IH = IHQMAT(IQ,ITIM)
        IF(IH.EQ.1) THEN
          CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PAO)
        ELSEIF(IH.EQ.2) THEN
          CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PAO)
        ENDIF
        CALL BUTOBS(PAO,1,WORK(KFREE),LFREE)
        IPA = 1
        IPF = 1
        IPR = 1
        DO IFRP = 1,NFSYM
          CALL QTRANS('AOMO','S',FADD,
     &         NFBAS(IFRP,0),NFBAS(IFRP,0),NCNT(IFRP),NCNT(IFRP),
     &         PAO(IPA),NTBAS(0),NTBAS(0),1,IPQTOQ(IZ,IREP),
     &         PFO(IPF),NCNT(IFRP),NCNT(IFRP),NZ,IPQTOQ(1,IREP),
     &         RVEC(IPR),NFBAS(IFRP,0),NCNT(IFRP),NZ,IPQTOQ(1,0),
     &         RVEC(IPR),NFBAS(IFRP,0),NCNT(IFRP),NZ,IPQTOQ(1,0),
     &         WORK(KFREE),LFREE,IPRINT)
          IPA = IPA + NFBAS(IFRP,0)*(NTBAS(0)+1)
          IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
          IPR = IPR + NFBAS(IFRP,0)*NCNT(IFRP)*NZ
        ENDDO
        FADD = D1
      ENDIF
      ENDDO
C
      IF(IPRINT.GE.5) THEN
        IPF = 1
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A16,A,I3)')
     &    'Property matrix of ',PRPNAM(INDXPR),
     &    'in fragment basis..',IFRP
          CALL PRQMAT(PFO(IPF),NCNT(IFRP),NCNT(IFRP),
     &         NCNT(IFRP),NCNT(IFRP),NZ,
     &         IPQTOQ(1,IREP),LUPRI)
          IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
        ENDDO       
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjex2 */
      SUBROUTINE PRJEX2(EXPVAL,PTRI,PFO,DFO,BVEC,IBUF,INDXPR,NCNT,NMO,
     &                  NR,NTOT,PM,FMAT,INDX,
     &                  IPRINT,PBUF,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Make generalized population matrix
C
C     Written by Trond Saue Apr 27 2003 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbprj.h"
C
      DIMENSION PFO(*),DFO(*),NCNT(*),NMO(*),
     &          BVEC(*),PTRI(NNBBASX),IBUF(*),PBUF(*),
     &          NR(NREFS,NFSYM),FMAT(NTOT,NTOT),PM(NTOT),INDX(2,*),
     &          WORK(*)
C
      N2TOT = NTOT*NTOT
      CALL DZERO(FMAT,N2TOT)
      CALL DZERO(PM,NTOT)
C
C     Make density matrix in fragment orbitals
C
      IPD = 1
      DO 10 IFRP = 1,NFSYM
        IF(NMO(IFRP).EQ.0) GOTO 10
        II = 0
        DO I = 1,NREFS
          DO J = 1,NR(I,IFRP)
            II = II + 1
            INDX(1,II) = I
            INDX(2,II) = J
          ENDDO
        ENDDO
        I = NREFS + 1
        DO J = 1,NMO(IFRP)
          II = II + 1
          INDX(1,II) = I
          INDX(2,II) = J
        ENDDO
        CALL PRJEX3(EXPVL2,PTRI,PFO(IPD),DFO(IPD),NCNT(IFRP),
     &              NR(1,IFRP),NMO(IFRP),IBUF,PBUF,FMAT,PM,INDX,
     &              NTOT,INDXPR,IPRINT)
        IPD = IPD + NCNT(IFRP)*NCNT(IFRP)*NZ
 10   CONTINUE
C
C     Give fragment contributions
C
      CALL HEADER('* Total fragment contributions',-1)
      WRITE(LUPRI,'(A)') ' * Intraatomic contributions:'
      TOT = D0
      DO J = 1,NREFS
        TMP = D2*FMAT(J,J)
        TOT = TOT + TMP
        FAC = TMP/EXPVAL
        WRITE(LUPRI,'(I5,2X,A6,16X,E16.8,F10.4)') 
     &    J,REFFIL(J),TMP,FAC
        WRITE(LUPRI,'(9X,A,E16.8,F10.4)')
     &    '- principal moment: ',PM(J),PM(J)/EXPVAL
        TMP = TMP - D2*PM(J)
        WRITE(LUPRI,'(9X,A,E16.8,F10.4)')
     &    '- hybridization   : ',TMP,TMP/EXPVAL
      ENDDO
      WRITE(LUPRI,'(3X,A,19X,E16.8,F10.4)') 'Total: ',TOT,TOT/EXPVAL
      WRITE(LUPRI,'(A)') ' * Interatomic contributions:'
      TOT = D0
      DO J = 1,NREFS
        DO I = 1,(J-1)
          TMP = D2*(FMAT(I,J)+FMAT(J,I))
          TOT = TOT + TMP
          FAC = TMP/EXPVAL
          WRITE(LUPRI,'(2(I5,2X,A6),3X,E16.8,F10.4)') 
     &      J,REFFIL(J),I,REFFIL(I),TMP,FAC
        ENDDO
      ENDDO
      WRITE(LUPRI,'(3X,A,19X,E16.8,F10.4)') 'Total: ',TOT,TOT/EXPVAL
      WRITE(LUPRI,'(A)') ' * Polarization contributions:'
      J = NREFS+1
      TOT = D0
      DO I = 1,NREFS
        TMP = D2*(FMAT(I,J)+FMAT(J,I))
        TOT = TOT + TMP
        FAC = TMP/EXPVAL
        WRITE(LUPRI,'(2(I5,2X,A6),3X,E16.8,F10.4)') 
     &      J,REFFIL(J),I,REFFIL(I),TMP,FAC
      ENDDO
      TMP = D2*FMAT(J,J)
      TOT = TOT + TMP
      FAC = TMP/EXPVAL
      WRITE(LUPRI,'(2(I5,2X,A6),3X,E16.8,F10.4)') 
     &    J,REFFIL(J),J,REFFIL(J),TMP,FAC
      WRITE(LUPRI,'(3X,A,19X,E16.8,F10.4)') 'Total: ',TOT,TOT/EXPVAL
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_efg */
      SUBROUTINE NUC_EFG(IPRINT,DOATOM)
C***********************************************************************
C     
C     Calculates nuclear contribution to electric field gradients
C
C     Written by T.Enevoldsen  - Sep 1996
C     Last revision by W.A. de Jong - April 1999
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D3 = 3.0D0)
      DIMENSION  GEOM(3*8*MXCENT)
      LOGICAL    DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
      CALL DZERO(EFGNUC,9*MXCENT)
      KK = 1
      NATOM = 1
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               LL = 1
               DO 100 I2 = 1, NUCIND
                  IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
                  DO 110 ISYM = 0, MAXOPR
                     IF (IBTAND(ISYM,ISTBNU(I2)) .EQ. 0) THEN
                        IF (LL .NE. KK) THEN
                           XCOOR = GEOM(LL) - GEOM(KK)
                           YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                           ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                           R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                           R5 = R2*R2*SQRT(R2)
                           EFGNUC(1,1,NATOM) = EFGNUC(1,1,NATOM)
     &                             + CHARGE(I2)*(D3*XCOOR*XCOOR - R2)/R5
                           EFGNUC(2,2,NATOM) = EFGNUC(2,2,NATOM)
     &                             + CHARGE(I2)*(D3*YCOOR*YCOOR - R2)/R5
                           EFGNUC(3,3,NATOM) = EFGNUC(3,3,NATOM)
     &                             + CHARGE(I2)*(D3*ZCOOR*ZCOOR - R2)/R5
                           EFGNUC(1,2,NATOM) = EFGNUC(1,2,NATOM)
     &                             + CHARGE(I2)*D3*XCOOR*YCOOR/R5
                           EFGNUC(1,3,NATOM) = EFGNUC(1,3,NATOM)
     &                             + CHARGE(I2)*D3*XCOOR*ZCOOR/R5
                           EFGNUC(2,3,NATOM) = EFGNUC(2,3,NATOM)
     &                             + CHARGE(I2)*D3*YCOOR*ZCOOR/R5
                        END IF
                        LL = LL + 3
                     END IF
 110              CONTINUE
 100           CONTINUE
               EFGNUC(2,1,NATOM) = EFGNUC(1,2,NATOM)
               EFGNUC(3,1,NATOM) = EFGNUC(1,3,NATOM)
               EFGNUC(3,2,NATOM) = EFGNUC(2,3,NATOM)
               IF (IPRINT .GE. 2) THEN
                  WRITE (LUPRI ,'(/2X,2A,/)')
     &               'Nuclear EFG Tensor of nucleus: ',NAMDEP(NATOM)
                  CALL OUTPUT(EFGNUC(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
               END IF
               NATOM = NATOM + 1
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck getiopen */
      INTEGER FUNCTION GETIOPEN(IORBITAL,IFRP)
C***********************************************************************
C     
C     Which shell does electron orbital IORBITAL belong to?
C
C     Written by J. Thyssen - Feb 11 1999
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbdhf.h"
#include "dcborb.h"
C
      IF (IORBITAL .LT. 1)
     &   CALL QUIT('illegal value in GETIOPEN(neg)')
C
      IF (IORBITAL .LE. NISH(IFRP)) THEN
         GETIOPEN = 0
         RETURN
      END IF
C
      ILOW = NISH(IFRP)+1
      DO IOPEN = 1,NOPEN
         IHIGH = ILOW + NACSH(IFRP,IOPEN)
         IF ( (IORBITAL .GE. ILOW) .AND. (IORBITAL .LE. IHIGH) ) THEN
            GETIOPEN = IOPEN
            RETURN
         END IF
         ILOW = IHIGH+1
      END DO
      CALL QUIT('illegal value in GETIOPEN(pos)')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpexp_mp2 */
      SUBROUTINE PRPEXP_MP2(DMO,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate correlation contribution to static first order properties
C     defined in /CBIEXP/
C
C     Based on PRPEXP
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Joost van Stralen, april 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION DMO(*)
      DIMENSION WORK(LFREE)
C
#include "dcbexp.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('PRPEXP_MP2')
C
C     We come from RELCCSD, The total memory available for the WORK array,
C     LWORK is equal to LFREE.
C
      KFRSAV=KFREE
C
C     Memory allocation
C     -----------------
      CALL MEMGET('REAL',KCMO ,N2BBASXQ   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT      ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KEXPST,2*NEXPP   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT ,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
      CALL MEMGET('REAL',KEXP  ,NEXPP     ,WORK,KFREE,LFREE)
C
      CALL MEMCHK('voor prpex1_mp2',WORK,1)
      CALL PRPEX1_MP2(WORK(KEXP),WORK(KEXPST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),WORK(KPTRI),
     &            DMO,WORK(KFIRST),WORK,
     &            KFREE,LFREE)
C
C
C     Memory deallocation
C     -------------------
      CALL MEMREL('PRPEXP_MP2',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPEXP_MP2')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex1_mp2 */
      SUBROUTINE PRPEX1_MP2(EXPVAL,LEXPST,CMO,EIG,DMAT,PMAT,PTRI,
     &                      DMO,FIRST,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate correlation contribution to static first order properties
C     defined in /CBIEXP/
C
C     Based on PRPEXP1
C
C
C     INPUT:   DMO - second order density matrix
C
C     Written by Joost van Stralen, april 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0)
C
      LOGICAL FIRST(NZ),LEXPST(2,*),FILEX
      DIMENSION EXPVAL(*),CMO(*),DMAT(*),PMAT(*),PTRI(*),
     &          EIG(*),DMO(*),WORK(*)
      CHARACTER MXFORM*6,FMT*6
C
#include "dgroup.h"
ctec dcborb included
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
C
      CALL QENTER('PRPEX1_MP2')
C
C     Initialize
C     ==========
C
c#include "memint.h"
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
      LWORK=LFREE
C
      CALL MEMCHK('voor dzero',WORK,1)
      CALL DZERO(EXPVAL,NEXPP)
      CALL MEMCHK('na dzero',WORK,1)
      N2EXP = 2*NEXPP
      CALL LSET(N2EXP,.FALSE.,LEXPST)
C
C     Open coefficient file
C     =====================
C
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP_MP2')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C
C     Transform DMO to 'orthonormal SO' basis
C     =======================================
C
      CALL DMP2AO(DMO,DMAT,CMO,WORK,KFREE,LFREE)
C  
C     Calculate correlation contribution to property
C     ==============================================
C
      CALL TITLER('MP2 correlation contribution to property','*',115)
      WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
      CALL PRSYMB(LUPRI,'-',76,0)
      DO I = 1,NEXPP
        INDXPR = LEXPP(I)
          CALL PRPEX2(LU1INT,INDXPR,EXPVAL(I),LEXPST(1,I),
     &                IREP,ITIM,PTRI,PMAT,FIRST,
     &                WORK(KFREE),LFREE,IPREXP)
          IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 30
          CALL PRPEX3(EXPVAL(I),IREP,ITIM,
     &        DMAT,PMAT,PTRI,FIRST,ORBANA,
     &        WORK(KFREE),LFREE,IPREXP)
 30       CONTINUE
        FMT = MXFORM(EXPVAL(I),15)
        WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')   
     &     PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.',
     &     LEXPST(1,I),LEXPST(2,I)
        CALL FLSHFO(LUPRI)
      ENDDO
      CALL MEMCHK('na verwachtingswaarde',WORK,1)
C
C     Write expectation values to file
C     ================================
C
      CALL WEXPVTOF(EXPVAL)
C
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(4X,A)')
     & 's0 = T : Expectation value zero by point group symmetry.'
      WRITE(LUPRI,'(4X,A)')
     & 't0 = T : Expectation value zero by time reversal symmetry.'
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
      CALL QEXIT('PRPEX1_MP2')
C
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPEXP_MP2: Coefficient file not found !'
      CALL QUIT('PRPEXP_MP2: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wexpvtof */
      SUBROUTINE WEXPVTOF(EXPVAL)
C*****************************************************************************
C
C     Write expectation values to file
C
C     Written by Joost van Stralen, June 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*)
#include "dcbexp.h"
#include "dcbprp.h"
C
      CALL QENTER('WEXPVTOF')
C
      OPEN (LUCORP,FILE='CORPROP',FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE (LUCORP) (EXPVAL(I),I=1,NEXPP)
      CLOSE (LUCORP,STATUS='KEEP')
C
      CALL QEXIT('WEXPVTOF')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readexpv */
      SUBROUTINE READEXPV(EXPVAL,TOBE)
C*****************************************************************************
C
C     read expectation values from file
C
C     Written by Joost van Stralen, June 2002
C     
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*)
      LOGICAL TOBE
#include "dcbexp.h"
#include "dcbprp.h"
C
      CALL QENTER('READEXPV')
C     
      CALL DZERO(EXPVAL,NEXPP)
C
      INQUIRE(FILE='CORPROP',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE (LUPRI,1000)
      ELSE
         OPEN (LUCORP,FILE='CORPROP',FORM='UNFORMATTED',STATUS='OLD')
         REWIND (LUCORP)
         READ (LUCORP) (EXPVAL(I),I=1,NEXPP)
         CLOSE (LUCORP)
      ENDIF
C
      CALL QEXIT('READEXPV')
C
 1000 FORMAT (//" WARNING : correlation contribution was not calculated"
     &/" did you remember to specify DOFOPR=T in RELCCSD ?"//)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck exptot */
      SUBROUTINE EXPTOT(EXPVAL,EXP_CORR,EXP_TOT)
C*****************************************************************************
C
C     put HF and correlation contribution to expectation value together
C
C     Written by Joost van Stralen, June 2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*), EXP_CORR(*), EXP_TOT(*)
#include "dcbexp.h"
#include "dcbprp.h"
C
      CALL QENTER('EXPTOT')
C
      CALL DZERO(EXP_TOT,NEXPP)
C
      DO I = 1,NEXPP
         EXP_TOT(I) = EXPVAL(I) + EXP_CORR(I)
      ENDDO
C
      CALL QEXIT('EXPTOT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_nef */
      SUBROUTINE NUC_EFN(IPRINT,DOATOM,EFNNUC,GEOM)
C*****************************************************************************
C
C     Calculates nuclear contribution to the electric field at 
C     individual nuclei
C
C     Based on NUC_EFG
C     Written by T.Saue  - Feb 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D3 = 3.0D0)
#include "nuclei.h"
#include "symmet.h"
      LOGICAL    DOATOM(*)
      DIMENSION  GEOM(3*NUCDEP),EFNNUC(3,NUCDEP)
#include "ibtfun.h"
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
      CALL DZERO(EFNNUC,3*NUCDEP)
      KK = 1
      NATOM = 1
      DO 10 I = 1, NUCIND
        DO 20 ISYMOP = 0, MAXOPR
        IF (IBTAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
          LL = 1
          DO 100 I2 = 1, NUCIND
          IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
            DO 110 ISYM = 0, MAXOPR
            IF (IBTAND(ISYM,ISTBNU(I2)) .EQ. 0) THEN
              IF (LL .NE. KK) THEN
                XCOOR = GEOM(LL) - GEOM(KK)
                YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                R3 = R2*SQRT(R2)
                EFNNUC(1,NATOM) = EFNNUC(1,NATOM)- CHARGE(I2)*XCOOR/R3
                EFNNUC(2,NATOM) = EFNNUC(2,NATOM)- CHARGE(I2)*YCOOR/R3
                EFNNUC(3,NATOM) = EFNNUC(3,NATOM)- CHARGE(I2)*ZCOOR/R3
              END IF
              LL = LL + 3
            END IF
 110      CONTINUE
 100      CONTINUE
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,2A,/)')
     &        'Nuclear electric field at nucleus: ',NAMDEP(NATOM)
            IOFF = ICHAR('w')
            DO K = 1,3
               WRITE(LUPRI,'(4X,A1,A,E20.10)') CHAR(IOFF+K),
     &              '-component: ',EFNNUC(K,NATOM)
            ENDDO
          END IF
          NATOM = NATOM + 1
          KK = KK + 3
        END IF
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck prjex3 */
      SUBROUTINE PRJEX3(EXPVAL,PTRI,PFO,DFO,NCNT,NR,NMO,IBUF,PBUF,
     &                  FMAT,PM,INDX,NTOT,INDXPR,IPRINT)
C***********************************************************************
C
C     Calculate generalized population matrix
C
C     Written by T. Saue Apr 27 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2=2.0D0,D0=0.0D0,DTOL=1.0D-4)
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbprj.h"
      DIMENSION PTRI(*),PFO(NCNT,NCNT,NZ),DFO(NCNT,NCNT,NZ),
     &          NR(*),PBUF(*),IBUF(*),FMAT(NTOT,NTOT),PM(NTOT),INDX(2,*)
C
#include "sdpre.h"
      NTRI = NCNT*(NCNT+1)/2
C
C     Make population matrix
C
      CALL DZERO(PTRI,NTRI)
      DO IZ = 1,NZ
        IJ = 0
        DO J = 1,NCNT
          JR = INDX(1,J)
          DO I = 1,(J-1)
            IR = INDX(1,I)
            IJ = IJ + 1
            TMP = D2*PFO(I,J,IZ)*DFO(I,J,IZ)
            FMAT(IR,JR) = FMAT(IR,JR) + TMP
            PTRI(IJ) = PTRI(IJ) + TMP
          ENDDO
          IJ = IJ + 1
          TMP = PFO(J,J,IZ)*DFO(J,J,IZ)
          PTRI(IJ) = PTRI(IJ) + TMP
          PM(JR) = PM(JR) + TMP
          FMAT(JR,JR) = FMAT(JR,JR) + TMP
        ENDDO
      ENDDO      
      CALL DSCAL(NTRI,D2,PTRI,1)
      EXPVAL = DSUM(NTRI,PTRI,1)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(A,A16,E16.10)') 
     &    '* Expectation value of property ',PRPNAM(INDXPR),EXPVAL
        CALL PRMUTC(PTRI,NCNT,1,LUPRI)      
      ENDIF
C
C     Make ordered list of contributions
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Ordered list of contributions',-1)            
        DO IJ = 1,NTRI
          PBUF(IJ) = -ABS(PTRI(IJ))
        ENDDO
        CALL INDEXX (NTRI,PBUF,IBUF)
        TMP = D0      
        WRITE(LUPRI,'(2(A4,10X),2A16,A10)')
     &    'Ind1','Ind2','Contribution','Accumulated','Fraction'
        DO M = 1,NTRI
          IJ = IBUF(M)
          IF(ABS(PTRI(IJ)/EXPVAL).LT.DTOL) GOTO 10
          J  = INT(SQRT(SDPRE(2*IJ) + 0.25D0) + 0.4999D0)
          I  = IJ - J*(J-1)/2
          JR = INDX(1,J)
          JJ = INDX(2,J)
          IR = INDX(1,I)
          II = INDX(2,I)
          TMP = TMP + PTRI(IJ)
          FAC = TMP/EXPVAL
          WRITE(LUPRI,'(2(A6,I5,3X),2E16.8,F10.4)') 
     &          REFFIL(JR),JJ,REFFIL(IR),II,PTRI(IJ),TMP,FAC
        ENDDO
 10     CONTINUE
      ENDIF
C      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjex4 */
      SUBROUTINE PRJEX4(DFO,BVEC,NCNT,NMO,IPRINT)
C*****************************************************************************
C
C     Make density matrix in fragment basis
C
C     Written by Trond Saue Oct 21 2003 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dgroup.h"
      DIMENSION DFO(*),NCNT(*),NMO(*),BVEC(*)
C
C     Make density matrix in fragment orbitals
C
      IPD = 1
      IPB = 1
      DO 10 IFRP = 1,NFSYM
        IF(NMO(IFRP).EQ.0) GOTO 10
        CALL DENST1(DFO(IPD),NCNT(IFRP),NCNT(IFRP),
     &              BVEC(IPB),NCNT(IFRP),NMO(IFRP),
     &              1,NMO(IFRP),NCNT(IFRP))  
        IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(A,I3)')
     &    'Density matrix in fragment basis..',IFRP
          CALL PRQMAT(DFO(IPD),NCNT(IFRP),NCNT(IFRP),
     &         NCNT(IFRP),NCNT(IFRP),NZ,
     &         IPQTOQ(1,0),LUPRI)
        ENDIF
        IPD = IPD + NCNT(IFRP)*NCNT(IFRP)*NZ
        IPB = IPB + NCNT(IFRP)*NMO(IFRP)*NZ
 10   CONTINUE
      RETURN
      END
