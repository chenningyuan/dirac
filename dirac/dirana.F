C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirana.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamana */
      SUBROUTINE PAMANA(WORK,LWORK)
C*****************************************************************************
C
C     Module for wave function analysis
C
C     Written by T.Saue Sep 11 1995
C     Last revision : April 27 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
C
      DIMENSION WORK(LWORK)
      CALL QENTER('PAMANA')
#include <memint.h>
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
C
      CALL TITLER('Wave function analysis','*',125)
      CALL PAMAN1(WORK(KCMO),WORK(KEIG),WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('PAMANA.PAMAN1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('PAMANA')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck paman1 */
      SUBROUTINE PAMAN1(CMO,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Analyze Dirac-Fock wavefunction
C
C     Written by T.Saue Sep 11 1995
C     Last revision : KApril 27 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbana.h>
#include <dcbbas.h>
C Used from COMMON blocks:
C  CBIANA: DOPOP,DOANA,DOEXP
C  DCBBAS: N2BBASX
C
      LOGICAL FILEX
      DIMENSION CMO(*),EIG(*),WORK(LWORK)
C
      CALL QENTER('PAMAN1')
C
C*****************************************************************************
C*****  R E A D    W A V E F U N C T I O N  **********************************
C*****************************************************************************
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
C
C     Read coefficients
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,TOTERG,2)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C*****************************************************************************
C*****  V E C T O R    P R I N T *********************************************
C*****************************************************************************
C
      IF(DOVEC) CALL ANAVEC(CMO,EIG,WORK,LWORK)
C
C*****************************************************************************
C*****  M U L L I K E N    P O P U L A T I O N     A N A L Y S I S  **********
C*****************************************************************************
C
      IF(DOPOP) CALL MULPOP(CMO,EIG,WORK,LWORK)
C
C*****************************************************************************
C*****  P R O J E C T I O N   O N T O    A N O T H E R   S O L U T I O N *****
C*****************************************************************************
C
      IF(DOPRJ) CALL ANAPRJ(CMO,EIG,WORK,LWORK)
C
      CALL QEXIT('PAMAN1')
C
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PAMANA: Coefficient file not found !'
      CALL QUIT('PAMANA: Coefficients not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpop */
      SUBROUTINE MULPOP(CMO,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Mulliken population analysis
C
C     Written by T.Saue - Sep 11 1995
C     Last revision: Aoril 27 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxorb.h>
C
C Used from COMMON blocks:
C  DCBBAS: N2BBASX
C  LABELS: NLAB
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbpop.h>
#include <dgroup.h>
C
      DIMENSION CMO(*),EIG(*),WORK(LWORK)
C
#include <dcbibt.h>
      CALL QENTER('MULPOP')
#include <memint.h>
C
C     Get labels for population analysis
C
      CALL GETPLB(IPRPOP)
C
C     Memory allocation
      NNLAB = (NPOPLAB*(NPOPLAB+1))/2
      N2LAB = NPOPLAB*NPOPLAB
      CALL MEMGET('INTE',KIPGR ,NTBAS(0)       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT ,N2LAB*2        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT ,N2BBASX        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOPG ,NPOPLAB*2      ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIEFF ,NPOPLAB        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPGTOT,NPOPLAB*NFSYM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPNTOT,N2LAB*NFSYM  ,WORK,KFREE,LFREE)
C
      CALL MULPO1(WORK(KPMAT),WORK(KSMAT),CMO,EIG,WORK(KIPGR),
     &     WORK(KPOPG),WORK(KIEFF),WORK(KPGTOT),WORK(KPNTOT),
     &     WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('MULPOP.MULPO1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPOP')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpo1 */
      SUBROUTINE MULPO1(PMAT,SMAT,CMO,EIG,IPGRS,POPG,IEFF,
     &                  PGTOT,PNTOT,WORK,LWORK)
C*****************************************************************************
C
C     This routine generates net and gross Mulliken populations
C     separately for each component of a 2-spinor. The array
C     IPLAB determines the grouping of basis functions into labels.
C
C     PNET - net population(triangularly row-packed)
C     PGRS - gross population
C     SMAT - overlap matrix
C     COEF - coefficients
C
C     Written by T.Saue July 24 1994
C     LAST VERSION: Jan 9 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D2 = 2.0D0)
#include <dummy.h>
#include <mxorb.h>
C Used from COMMON blocks:
C  DCBORB: NSPH
C  DCBBAS:
C  DGROUP: NZ
C
#include <dcbgen.h>
#include <dcbpop.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dgroup.h>
C
      DIMENSION PMAT(NPOPLAB,NPOPLAB,2),SMAT(*),CMO(*),EIG(*),IPGRS(*),
     &          POPG(NPOPLAB,2),IEFF(NPOPLAB),
     &          PNTOT(NPOPLAB,NPOPLAB,NFSYM),
     &          PGTOT(NPOPLAB,NFSYM),WORK(LWORK)
C
      CALL QENTER('MULPO1')
#include <memint.h>
C
      CALL TITLER('Mulliken population analysis','*',122)
C
      N2LAB = NPOPLAB*NPOPLAB
      CALL DZERO(PNTOT,N2LAB*NFSYM)
      CALL DZERO(PGTOT,NPOPLAB*NFSYM)
C
C     Make pointers to group labels
C     =============================
C
      DO I = 1,NTBAS(0)
        IPGRS(I) = IPOPLAB(IPLAB(I,ILABDF))
      ENDDO
C
C     Get overlap matrix
C     ==================
C
      CALL GTOVLP(SMAT,SSMTRC,WORK,LWORK)
C
C     If necessary, backtransform overlap matrix to AO-basis
C     ======================================================
C
      IF(ILABDF.EQ.1.AND.NBSYM.GT.1) THEN
        CALL MEMGET('REAL',KBUF ,N2BBASX ,WORK,KFREE,LFREE)
        CALL MTBSBU(SMAT,WORK(KBUF))
        CALL MTSOAO(WORK(KBUF),SMAT,NTBAS(0),0,IPRPOP)
        CALL MEMREL('MULPOP.MULPO1',WORK,KWORK,KWORK,KFREE,LFREE)
      ENDIF
C
C     Loop over fermion ircops
C     ========================
C
      DO 10 IFRP = 1,NFSYM
C
C       Find number of electronic/positronic vectors
C       ============================================
C
        NVEC = 1
        CALL MEMGET('INTE',KJVEC,NORB(IFRP),WORK,KFREE,LFREE)
        CALL  NUMLST(VECPOP(IFRP),WORK(KJVEC),NORB(IFRP),
     &             -NPSH(IFRP),NESH(IFRP),NVEC,WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KJVEC),NVEC,IFRP,NPVEC,NEVEC)
        IF(NVEC.EQ.0) GOTO 20
        CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
C
        IF(ILABDF.EQ.1) THEN
C
C       Extract/expand coefficients in AO-basis
C       =======================================
C
          NCBAS = NTBAS(0)
          IF(NVEC.EQ.NORB(IFRP)) THEN
            IF(NBSYM.EQ.1) THEN
              NCBF1 = NTBAS(0)*NVEC*4
              CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
              CALL CFSORT(WORK(KCBF1),CMO(ICMOQ(IFRP)+1),NVEC,IFRP,
     &                    NFBAS(IFRP,0))
            ELSE
              NCBF1 = NTBAS(0)*NVEC*4
              NCBF2 = NFBAS(IFRP,0)*NVEC*4
              CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
              CALL CFSORT(WORK(KCBF2),CMO(ICMOQ(IFRP)+1),NVEC,IFRP,
     &                    NFBAS(IFRP,0))
              CALL WTSOAO(WORK(KCBF2),WORK(KCBF1),IFRP,
     &                    NFBAS(IFRP,0),NTBAS(0),NVEC*4,IPRPOP)
              CALL MEMREL('MULPO1.AO',WORK,1,KCBF2,KFREE,LFREE)
            ENDIF
          ELSE
            NCBF1 = NTBAS(0)*NVEC*4
            NCBF2 = NTBAS(0)*NVEC*4
            CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
            IF(NBSYM.EQ.1) THEN
              CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF2),
     &                    WORK(KJVEC),NPVEC,NEVEC,
     &                    NFBAS(IFRP,0),NORB(IFRP))
              CALL CFSORT(WORK(KCBF1),WORK(KCBF2),NVEC,IFRP,
     &                    NFBAS(IFRP,0))
            ELSE
              CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF1),
     &                    WORK(KJVEC),NPVEC,NEVEC,
     &                    NFBAS(IFRP,0),NORB(IFRP))
              CALL CFSORT(WORK(KCBF2),WORK(KCBF1),NVEC,IFRP,
     &                    NFBAS(IFRP,0))
              CALL WTSOAO(WORK(KCBF2),WORK(KCBF1),IFRP,
     &                    NFBAS(IFRP,0),NTBAS(0),NVEC*4,IPRPOP)
            ENDIF
          ENDIF
        ELSE
C
C       Extract/expand coefficients in SO-basis
C       =======================================
C
          NCBAS = NFBAS(IFRP,0)
          IF(NVEC.EQ.NORB(IFRP)) THEN
            NCBF  = NFBAS(IFRP,0)*NVEC
            NCBF1 = NCBF*4
            CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
            CALL CFSORT(WORK(KCBF1),CMO(ICMOQ(IFRP)+1),NVEC,IFRP,
     &                    NFBAS(IFRP,0))
          ELSE
            NCBF  = NFBAS(IFRP,0)*NVEC
            NCBF1 = NCBF*4
            NCBF2 = NCBF*NZ
            CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
            CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF2),
     &                  WORK(KJVEC),NPVEC,NEVEC,
     &                  NFBAS(IFRP,0),NORB(IFRP))
            CALL CFSORT(WORK(KCBF1),WORK(KCBF2),NVEC,IFRP,
     &                    NFBAS(IFRP,0))
            CALL MEMREL('MULPO1.SO',WORK,1,KCBF2,KFREE,LFREE)
          ENDIF
        ENDIF
C
C       Perform Mulliken population analysis
C       ====================================
C
        CALL MULPO2(PMAT,IFRP,SMAT,WORK(KCBF1),EIG(IORB(IFRP)+1),
     &              IPGRS,POPG,IEFF,PGTOT,PNTOT,NPVEC,NEVEC,NCBAS,
     &              WORK(KJVEC),WORK(KFREE),LFREE)
 20     CONTINUE
        CALL MEMREL('MULPO1',WORK,1,KWORK,KFREE,LFREE)
 10   CONTINUE
C
C     Write out total poulations
C     ==========================
C
      IF(NFSYM.GT.1) THEN
        WRITE(LUPRI,'(//A)') '*** Total gross poulation ***'
        CALL DAXPY(NPOPLAB,D1,PGTOT(1,2),1,PGTOT(1,1),1)        
        CALL WRPOPT(PGTOT(1,1),POPLAB,NPOPLAB,IEFF)                
        IF(DONETP) THEN
          CALL DAXPY(N2LAB,D1,PNTOT(1,1,2),1,PNTOT(1,1,1),1)
          WRITE(LUPRI,'(//A)') '*** Total net poulation ***'
          CALL WRPOPN(PNTOT(1,1,1),PGTOT(1,1),
     &                POPLAB,NPOPLAB,IEFF,'total')
        ENDIF
      ENDIF
C
C     Memory deallocation
      CALL MEMREL('MULPO1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPO1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpo2 */
      SUBROUTINE MULPO2(PMAT,IFRP,SMAT,CMO,EIG,IPGRS,POPG,IEFF,
     &           PGTOT,PNTOT,NPVEC,NEVEC,NCBAS,JVEC,WORK,LWORK)
C*****************************************************************************
C
C     This routine generates net and gross Mulliken populations
C     separately for each component of a 2-spinor. The array
C     IPLAB determines the grouping of basis functions into labels.
C
C     PNET - net population(triangularly row-packed)
C     PGRS - gross population
C     SMAT - overlap matrix
C     COEF - coefficients
C
C     Written by T.Saue July 24 1994
C     Last revision: Jan 7 1998 - jth LINUX
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D2 = 2.0D0)
#include <dummy.h>
#include <mxorb.h>
C Used from COMMON blocks:
C  DCBORB: NSPH
C  DCBBAS:
C  DGROUP: NZ
C
#include <dcbgen.h>
#include <dcbpop.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
C 
#ifdef VAR_G77
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0))
#endif
      CHARACTER FMT*6,MXFORM*6
      DIMENSION PMAT(NPOPLAB,NPOPLAB,2),JVEC(*),
     &          CMO(NCBAS,NPVEC+NPVEC,NZ),EIG(*),IPGRS(*),
     &          POPG(NPOPLAB,2),IEFF(NPOPLAB),
     &          PNTOT(NPOPLAB,NPOPLAB,NFSYM),
     &          PGTOT(NPOPLAB,NFSYM),WORK(LWORK)
C
      CALL QENTER('MULPO2')
#include <memint.h>
C
      NVEC = NPVEC + NEVEC
      DO I = 1,NVEC
        IF(I.GT.NPVEC) THEN
          II = NPSH(IFRP)+JVEC(I)
          FMT = MXFORM(EIG(II),16)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +       '* Electronic eigenvalue nr.',JVEC(I),': ',
     &       EIG(II)
        ELSE
          II  = NPSH(IFRP)+1+JVEC(I)
          FMT = MXFORM(EIG(II),16)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +       '* Positronic eigenvalue nr.',II,': ',
     &       EIG(II)
        ENDIF
        CALL PRSYMB(LUPRI,'=',48,0)
C
C       Calculate SO population matrix
C       ==============================
C
        CALL POPMAT(IFRP,I,PMAT,SMAT,CMO,NVEC,
     &              IPGRS,NPOPLAB,NCBAS,WORK,LWORK,IPRPOP)
C
C       Calculate gross populations
C       ===========================
C
        CALL GRSPOP(POPG,PMAT,NPOPLAB)
        CALL WRPOPG(POPG,POPLAB,NPOPLAB,IEFF)
        CALL DAXPY(NPOPLAB,D2,POPG(1,1),1,PGTOT(1,IFRP),1)
        CALL DAXPY(NPOPLAB,D2,POPG(1,2),1,PGTOT(1,IFRP),1)
C
C       Calculate net populations
C       ===========================
C
        IF(DONETP) THEN
          CALL WRPOPN(PMAT(1,1,1),POPG(1,1),POPLAB,NPOPLAB,
     &         IEFF,'alpha')
          CALL WRPOPN(PMAT(1,1,2),POPG(1,2),POPLAB,NPOPLAB,
     &         IEFF,'beta ')
          CALL DAXPY(N2LAB,D2,PMAT(1,1,1),1,PNTOT(1,1,IFRP),1)
          CALL DAXPY(N2LAB,D2,PMAT(1,1,2),1,PNTOT(1,1,IFRP),1)
        ENDIF
      ENDDO
C
C     Write out total poulations
C     ==========================
C
      IF(NFSYM.EQ.1) THEN
        WRITE(LUPRI,'(//A)') '*** Total gross poulation ***'
      ELSE
        WRITE(LUPRI,'(//A,A3,A)') 
     &    '** Total gross population of fermion ircop ',
     &    FREP(IFRP),' **'
      ENDIF
      CALL WRPOPT(PGTOT(1,IFRP),POPLAB,NPOPLAB,IEFF)                
      IF(DONETP) THEN
        IF(NFSYM.EQ.1) THEN
          WRITE(LUPRI,'(//A)') '*** Total net poulation ***'
        ELSE
          WRITE(LUPRI,'(//A,A3,A)') 
     &      '** Total net population of fermion ircop ',
     &      FREP(IFRP),' **'
        ENDIF
        CALL WRPOPN(PNTOT(1,1,IFRP),PGTOT(1,IFRP),
     &              POPLAB,NPOPLAB,IEFF,'total')
      ENDIF
      CALL MEMREL('MULPO1.POPMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPO2')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popmat */
      SUBROUTINE POPMAT(IFRP,IVEC,PMAT,SMAT,CMO,NVEC,IGLAB,NGLAB,
     &                  NCBAS,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Construct population matrix
C         M(i,j) = Re(CONJ[CMO(i,IVEC)]*SMAT(i,j)*CMO(j,IVEC))
C     The final population matrix is compressed based on a labeling system.
C
C     Written by T.Saue Oct 9 1995
C     Last revision: Oct 9 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C Used from COMMON blocks
C  DCBAS: N2BBASX
#include <dcbbas.h>
      DIMENSION PMAT(*),SMAT(*),CMO(*),IGLAB(*),WORK(LWORK)
C
      CALL QENTER('POPMAT')
#include <memint.h>
C
C     Memory allocation
      CALL MEMGET('REAL',KPBUF,2*N2BBASX,WORK,KFREE,LFREE)
C
      CALL POPMA1(IFRP,IVEC,PMAT,SMAT,CMO,NVEC,WORK(KPBUF),
     &     IGLAB,NGLAB,NCBAS,WORK(KFREE),LFREE,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('POPMAT.POPMA1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('POPMAT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popma1 */
      SUBROUTINE POPMA1(IFRP,IVEC,PMAT,SMAT,CMO,NVEC,PBUF,
     &                  IGLAB,NGLAB,NCBAS,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Construct population matrix
C         M(i,j) = Re(CONJ[CMO(i,IVEC)]*SMAT(i,j)*CMO(j,IVEC))
C     The final population matrix is compressed based on a labeling system.
C
C     Written by T.Saue Oct 9 1995
C     Last revision: Jan 7 1998 - jth LINUX
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0)
#include <mxorb.h>
C
C Used from COMMON blocks
C  DCBBAS: NTBAS,NFBAS
C  DCBORB: NORB
C  DGROUP: NZ
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbpop.h>
C
#ifdef VAR_G77
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST),PBUF(NBAST,NBAST,*)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0)),PBUF(NTBAS(0),NTBAS(0),*)
#endif
      DIMENSION PMAT(NGLAB,NGLAB,2),CMO(NCBAS,NVEC,4),IGLAB(*),
     &          WORK(LWORK)
C
C     Generate population matrix for basis
C     ====================================
C
      SUMA = D0
      SUMB = D0
      IOFF = IBAS(IFRP)
      CALL DZERO(PBUF,N2BBASX*2)
      IF(ILABDF.EQ.1) THEN
        DO J = 1,NTBAS(0)
          DO I = 1,NTBAS(0)
            PBUF(I,J,1) = (CMO(I,IVEC,1)*CMO(J,IVEC,1)+
     &                    CMO(I,IVEC,2)*CMO(J,IVEC,2))*SMAT(I,J)
            SUMA = SUMA + PBUF(I,J,1)
            PBUF(I,J,2) = (CMO(I,IVEC,3)*CMO(J,IVEC,3)+
     &                    CMO(I,IVEC,4)*CMO(J,IVEC,4))*SMAT(I,J)
            SUMB = SUMB + PBUF(I,J,2)
          ENDDO
        ENDDO
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(3X,A)') 'POPMAT: AO-basis'
          WRITE(LUPRI,'(6X,A,F9.4)') 'alpha-part:',SUMA
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(1,1,1),1,NTBAS(0),
     &                1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
          WRITE(LUPRI,'(6X,A,F9.4)') 'beta-part :',SUMB
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(1,1,2),1,NTBAS(0),
     &                  1,NTBAS(0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF
      ELSE
        DO J = 1,NFBAS(IFRP,0)
          JJ = IOFF + J
          DO I = 1,NFBAS(IFRP,0)
            II = IOFF + I
            PBUF(II,JJ,1) = (CMO(I,IVEC,1)*CMO(J,IVEC,1)+
     &                       CMO(I,IVEC,2)*CMO(J,IVEC,2))*SMAT(II,JJ)
            SUMA = SUMA + PBUF(II,JJ,1)
            PBUF(II,JJ,2) = (CMO(I,IVEC,3)*CMO(J,IVEC,3)+
     &                       CMO(I,IVEC,4)*CMO(J,IVEC,4))*SMAT(II,JJ)
            SUMB = SUMB + PBUF(II,JJ,2)
          ENDDO
        ENDDO
        IF(IPRINT.GE.2) THEN
          IOFF = IBAS(IFRP) + 1
          WRITE(LUPRI,'(3X,A)') 'POPMAT: SO-basis'
          WRITE(LUPRI,'(6X,A,F9.4)') 'alpha-part:',SUMA
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(IOFF,IOFF,1),1,NFBAS(IFRP,0),
     &                1,NFBAS(IFRP,0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
          WRITE(LUPRI,'(6X,A,F9.4)') 'beta-part :',SUMB
          IF(IPRINT.GE.4) THEN
            CALL OUTPUT(PBUF(IOFF,IOFF,2),1,NFBAS(IFRP,0),
     &                  1,NFBAS(IFRP,0),NTBAS(0),NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF
      ENDIF
C
C     Compress population matrix for labels
C     =====================================
C
      N2PMAT = NGLAB*NGLAB*2
      CALL DZERO(PMAT,N2PMAT)
      CALL GATMAT(1,NTBAS(0),NTBAS(0),PBUF(1,1,1),PMAT(1,1,1),
     &            IGLAB,NGLAB)
      CALL GATMAT(1,NTBAS(0),NTBAS(0),PBUF(1,1,2),PMAT(1,1,2),
     &            IGLAB,NGLAB)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(3X,A)') 'POPMAT: labels'
        SUMA = D0
        DO J = 1,NGLAB
          DO I = 1,NGLAB
            SUMA = SUMA + PMAT(I,J,1)
          ENDDO
        ENDDO
        WRITE(LUPRI,'(6X,A,F9.4)') 'alpha-part:',SUMA
        IF(IPRINT.GE.3) THEN
          CALL OUTPUT(PMAT(1,1,1),1,NGLAB,
     &                1,NGLAB,NGLAB,NGLAB,1,LUPRI)
        ENDIF
        SUMB = D0
        DO J = 1,NGLAB
          DO I = 1,NGLAB
            SUMB = SUMB + PMAT(I,J,2)
          ENDDO
        ENDDO
        WRITE(LUPRI,'(6X,A,F9.4)') 'beta-part :',SUMB
        IF(IPRINT.GE.3) THEN
          CALL OUTPUT(PMAT(1,1,2),1,NGLAB,
     &                1,NGLAB,NGLAB,NGLAB,1,LUPRI)
        ENDIF
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grspop */
      SUBROUTINE GRSPOP(POPG,PMAT,NLAB)
C*****************************************************************************
C
C     From Mulliken population matrix generate gross populations for
C     2-spinor.
C
C     Written by T.Saue - July 24 1994
C     LAST VERSION: Oct 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION POPG(NLAB,2),PMAT(NLAB,NLAB,2)
C
      DO 10 J = 1,NLAB
        POPG(J,1) = DSUM(NLAB,PMAT(1,J,1),1)
        POPG(J,2) = DSUM(NLAB,PMAT(1,J,2),1)
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopg */
      SUBROUTINE WRPOPG(POPG,LAB,NLAB,IEFF)
C*****************************************************************************
C
C     Output routine for Mulliken gross population
C
C     Written by T.Saue - Oct 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DTOL=0.0001D0,D0 = 0.0D0)
      CHARACTER LAB(NLAB)*12,FMT*7
      DIMENSION POPG(NLAB,2),IEFF(*)
C
C     Find nonzero contributions
C     ==========================
C
      SUMA = D0
      SUMB = D0
      NEFF = 0
      DO I = 1,NLAB
        SUMA = SUMA + POPG(I,1)
        SUMB = SUMB + POPG(I,2)
        POPT = ABS(POPG(I,1))+ABS(POPG(I,2))
        IF(POPT.GE.DTOL) THEN
          NEFF = NEFF + 1
          IEFF(NEFF) = I
        ENDIF
      ENDDO
C
C     Write out nonzero contributions
C     ===============================
C
      NROW = (NEFF+1)/8
      IF(MOD((NEFF+1),8).GT.0) NROW = NROW + 1
      NPRI = MIN(NEFF+1,8)
      NFMT = 8 + 15*NPRI
      NPRI = NPRI - 1
      WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
      WRITE(LUPRI,'(/A8,2X,A6,2X,A1,4X,7(A12,3X))')
     &     'Gross   ','Total ','|',(LAB(IEFF(I)),I=1,NPRI)
      WRITE(LUPRI,FMT) ('-',I=1,NFMT)
      WRITE(LUPRI,'(A8,F8.4,2X,A1,4X,7(F8.4,7X))')
     &     ' alpha  ',SUMA,'|',(POPG(IEFF(I),1),I=1,NPRI)
      WRITE(LUPRI,'(A8,F8.4,2X,A1,4X,7(F8.4,7X))')
     &     ' beta   ',SUMB,'|',(POPG(IEFF(I),2),I=1,NPRI)
      NOFF = NPRI
      DO J = 2,NROW
        NPRI = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NPRI
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(1X,A12,2X))')
     &       'Gross  |',(LAB(IEFF(NOFF+I)),I=1,NPRI)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        WRITE(LUPRI,'(A8,8(1X,F8.4,6X))')
     &       ' alpha |',(POPG(IEFF(NOFF+I),1),I=1,NPRI)
        WRITE(LUPRI,'(A8,8(1X,F8.4,6X))')
     &       ' beta  |',(POPG(IEFF(NOFF+I),2),I=1,NPRI)
        NOFF = NOFF + NPRI
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck netpop */
      SUBROUTINE NETPOP(POPN,PMAT,NLAB)
C*****************************************************************************
C
C     Mulliken net populations (upper triangular columnwise packed)
C
C     Written  by T.Saue Oct 11 1995
C     Last revision: Oct 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D2 = 2.0D0)
      DIMENSION POPN((NLAB*(NLAB+1)/2),2),PMAT(NLAB,NLAB,2)
C
C     Generate net populations
C
      IND = 0
      DO J = 1,NLAB
        IND = IND + 1
        POPN(IND,1) = PMAT(J,J,1)
        POPN(IND,2) = PMAT(J,J,2)
        DO I = 1,J
          IND = IND + 1
          POPN(IND,1) = D2*PMAT(I,J,1)
          POPN(IND,2) = D2*PMAT(I,J,2)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopn */
      SUBROUTINE WRPOPN(PMAT,POPG,LAB,NLAB,IEFF,COMP)
C*****************************************************************************
C
C     Output routine for Mulliken net populations
C
C     Written by T.Saue Oct 11 1995
C     Last revision Oct 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DTOL=0.0001D0,D2 = 2.0D0)
      CHARACTER LAB(NLAB)*12,FMT*7,PFMT*23,COMP*5
      DIMENSION PMAT(NLAB,NLAB),POPG(NLAB),IEFF(NLAB)
C
C     Find nonzero contributions
C     ==========================
C
      NEFF = 0
      DO I = 1,NLAB
        POPT = ABS(POPG(I))
        IF(POPT.GE.DTOL) THEN
          NEFF = NEFF + 1
          IEFF(NEFF) = I
        ENDIF
      ENDDO
      IF(NEFF.EQ.0) RETURN
C
      CALL HEADER('Net population -- '//COMP,0)
C
C     Write net populations
C     =====================
C
      NBATCH = NEFF/8
      IF(MOD(NEFF,8).GT.0) NBATCH = NBATCH + 1
      NOFF = 0
      DO J = 1,NBATCH
        NEND = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NEND
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(A12,3X))')
     &       'Net.pop.',(LAB(IEFF(NOFF+K)),K=1,NEND)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        NT = NOFF+NEND
        DO I = 1,NT
          NBEG = MAX(NOFF+1,I)
          NPRI = NT - NBEG + 1
          NX   = NEND-NPRI
          IF(NX.GT.0) THEN
            WRITE(PFMT,'(A,I4,A,I2,A)')
     &         '(A12,',15*NX,'X,',NPRI,'(F8.4,7X))'
          ELSE
            WRITE(PFMT,'(A,I2,A)') '(A12,',NPRI,'(F8.4,7X))      '
          ENDIF
          IF(I.EQ.NBEG) THEN
            WRITE(LUPRI,PFMT) LAB(IEFF(I)),
     &         PMAT(IEFF(I),IEFF(I)),
     &         (D2*PMAT(IEFF(NOFF+K),IEFF(I)),K=NX+2,NEND)
          ELSE
            WRITE(LUPRI,PFMT) LAB(IEFF(I)),
     &       (D2*PMAT(IEFF(NOFF+K),IEFF(I)),K=NX+1,NEND)
          ENDIF
        ENDDO
        NOFF = NOFF + NEND
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpprj */
      SUBROUTINE ANAPRJ(CMO,EIG,WORK,LWORK)
C*****************************************************************************
C
C
C     Project molecular solution down on another solution(e.g. atomic)
C     in the same basis
C
C         Y(mol)_i = Y(ref)_jc_ji
C
C     Projection coeffients:  c_ji = <Y(ref)_j|Y(mol)_i>
C     Coefficients for Y(atom) are read from file REFFIL.
C     PROTHR is threshold for print of projection coeffiencients.
C
C     Written by Trond Saue March 13 1996
C     Last revision April 27 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      LOGICAL FILEX
      DIMENSION CMO(*),EIG(*),WORK(LWORK)
C Used from COMMON blocks:
C     DCBBAS: N2BBASX,N2BBASXQ
C     DCBORB: NOCC
C     DGROUP: NZ
#include <dcbprj.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbham.h>
#include <dgroup.h>
C
      CALL QENTER('ANAPRJ')
#include <memint.h>
C
C     Memory allocation
      CALL MEMGET('REAL',KSMAT,N2BBASX ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF,NCMOTQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEREF,NTBAS(0),WORK,KFREE,LFREE)
C
C     Get overlap matrix
C     ===================
C
      CALL GTOVLP(WORK(KSMAT),SSMTRC,WORK(KFREE),LFREE)
C
C     Read reference coefficients
C     ===========================
C
C     Check if coefficients are on file
C
      INQUIRE(FILE=REFFIL,EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
C
C     Read coefficients
C
      LUATOM = 1
      CALL OPNFIL(LUATOM,REFFIL,'OLD','ANAAN1')
      CALL REACMO(LUATOM,REFFIL,WORK(KCREF),WORK(KEREF),TOTERG,2)
      CLOSE(LUATOM,STATUS='KEEP')
C
C     Perform projection
C     ==================     
C
      CALL ANAPR1(CMO,EIG,WORK(KSMAT),WORK(KCREF),WORK(KEREF),
     &            WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('ANAPRJ.ANAAN1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ANAPRJ')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A,A6,A)') 'ANAPRJ: Coefficient file ',REFFIL,
     &                        ' not found !'
      CALL QUIT('ANAPRJ: Reference coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anapr1 */
      SUBROUTINE ANAPR1(CMO,EIG,SMAT,CREF,EREF,WORK,LWORK)
C*****************************************************************************
C
C
C     Project molecular solution down on another solution(e.g. atomic)
C     in the same basis
C
C         Y(mol)_i = Y(ref)_jc_ji
C
C     Projection coeffients:  c_ji = <Y(ref)_j|Y(mol)_i>
C     Coefficients for Y(atom) are read from file PROFIL.
C     PROTHR is threshold for print of projection coeffiencients.
C
C     Written by Trond Saue March 13 1996
C     Last revision Mar 13 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      DIMENSION CMO(*),EIG(*),SMAT(*),CREF(*),EREF(*),WORK(LWORK)
#include <dgroup.h>
#include <dcbprj.h>
#include <dcborb.h>
#include <dcbbas.h>
C
#include <memint.h>
C
C     Loop over fermion ircops
C     ========================
C
      DO 10 IFRP = 1,NFSYM
C
C       Find number of electronic/positronic vectors
C       of coefficients to be projected
C       ============================================
C
        NVEC = 1
        CALL MEMGET('INTE',KJVEC,NORB(IFRP),WORK,KFREE,LFREE)
        CALL  NUMLST(VECPRJ(IFRP),WORK(KJVEC),NORB(IFRP),
     &             -NPSH(IFRP),NESH(IFRP),NVEC,WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KJVEC),NVEC,IFRP,NPVEC,NEVEC)
        IF(NVEC.EQ.0) GOTO 20
C
C       Find number of electronic/positronic vectors
C       of reference coefficents
C       ============================================
C
        NREF = 1
        CALL MEMGET('INTE',KJREF,NORB(IFRP),WORK,KFREE,LFREE)
        CALL  NUMLST(VECREF(IFRP),WORK(KJREF),NORB(IFRP),
     &             -NPSH(IFRP),NESH(IFRP),NREF,WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KJREF),NREF,IFRP,NPREF,NEREF)
        IF(NREF.EQ.0) GOTO 20
C
C       Extract/expand coefficients to be projected
C       ===========================================
C
        CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
        IF(NVEC.EQ.NORB(IFRP)) THEN
          NCBF  = NFBAS(IFRP,0)*NVEC
          NCBF1 = NCBF*4
          CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
          CALL CFSORT(WORK(KCBF1),CMO(ICMOQ(IFRP)+1),NVEC,IFRP,
     &                NFBAS(IFRP,0))
        ELSE
          NCBF  = NFBAS(IFRP,0)*NVEC
          NCBF1 = NCBF*4
          NCBF2 = NCBF*NZ
          CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
          CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF2),
     &                WORK(KJVEC),NPVEC,NEVEC,
     &                NFBAS(IFRP,0),NORB(IFRP))
          CALL CFSORT(WORK(KCBF1),WORK(KCBF2),NVEC,IFRP,
     &                NFBAS(IFRP,0))
          CALL MEMREL('ANAPR1.CMO',WORK,1,KCBF2,KFREE,LFREE)
        ENDIF
C
C       Extract/expand reference coefficients
C       =====================================
C
        IF(NREF.EQ.NORB(IFRP)) THEN
          NRBF  = NFBAS(IFRP,0)*NREF
          NRBF1 = NRBF*4
          CALL MEMGET('REAL',KRBF1,NRBF1,WORK,KFREE,LFREE)
          CALL CFSORT(WORK(KRBF1),CREF(ICMOQ(IFRP)+1),NREF,IFRP,
     &                NFBAS(IFRP,0))
        ELSE
          NRBF  = NFBAS(IFRP,0)*NREF
          NRBF1 = NRBF*4
          NRBF2 = NRBF*NZ
          CALL MEMGET('REAL',KRBF1,NRBF1,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KRBF2,NRBF2,WORK,KFREE,LFREE)
          CALL SELCFS(CREF(ICMOQ(IFRP)+1),IFRP,WORK(KRBF2),
     &                WORK(KJREF),NPREF,NEREF,
     &                NFBAS(IFRP,0),NORB(IFRP))
          CALL CFSORT(WORK(KRBF1),WORK(KRBF2),NREF,IFRP,
     &                NFBAS(IFRP,0))
          CALL MEMREL('ANAPR1.REF',WORK,1,KRBF2,KFREE,LFREE)
        ENDIF
C
C       Perform projection
C       ==================
C
        NPROJ = NVEC*NREF*4
        CALL MEMGET('REAL',KPROJ,NPROJ,WORK,KFREE,LFREE)
        CALL ANAPR2(WORK(KFREE),LFREE,WORK(KPROJ),IFRP,SMAT,
     &       WORK(KCBF1),EIG (IORB(IFRP)+1),NPVEC,NEVEC,WORK(KJVEC),
     &       WORK(KRBF1),EREF(IORB(IFRP)+1),NPREF,NEREF,WORK(KJREF),
     &       NFBAS(IFRP,0))
 20     CONTINUE
        CALL MEMREL('ANAPR1',WORK,1,KWORK,KFREE,LFREE)
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anapr2 */
      SUBROUTINE ANAPR2(WORK,LWORK,CPROJ,IFRP,SMAT,
     &                  CVEC,EVEC,NPVEC,NEVEC,JVEC,
     &                  CREF,EREF,NPREF,NEREF,JREF,
     &                  KNFBAS)
C*****************************************************************************
C
C
C     Project molecular solution down onto reference solution(e.g. atomic)
C     in the same basis
C
C         Y(mol)_i = Y(ref)_jc_ji
C
C     Projection coeffients:  c_ji = <Y(ref)_j|Y(mol)_i>
C     PROTHR is threshold for print of projection coeffiencients.
C
C     Written by Trond Saue March 13 1996
C     Last revision Jan 7 1998 - jth NFBAS(IFRP,0) is now parameter
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1 = -1.0D0)
C
#include <dgroup.h>
#include <dcbprj.h>
#include <dcborb.h>
#include <dcbbas.h>
C
C     KNFBAS = NFBAS(IFRP,0)
C
      LOGICAL FIRST
      DIMENSION CVEC(KNFBAS,NEVEC+NPVEC,4),EVEC(*),JVEC(*),
     &          CREF(KNFBAS,NEREF+NPREF,4),EREF(*),JREF(*),
     &          CPROJ(NEREF+NPREF,NEVEC+NPVEC,4),SMAT(*),WORK(LWORK)
      CHARACTER FMT*6,MXFORM*6,BLINE*105
      DIMENSION IQBF(4),IQCF(4)
      DATA IQCF/1,2,3,4/
      DATA IQBF/1,2,3,4/
      NVEC = NEVEC + NPVEC
      NREF = NEREF + NPREF
C
C     Perform projection
C     ==================
C
      CALL DSCAL(NVEC*NFBAS(IFRP,0),DM1,CVEC(1,1,3),1)
      CALL DSCAL(NREF*NFBAS(IFRP,0),DM1,CREF(1,1,3),1)
      CALL QTRANS('AOMO','S',D0,NFBAS(IFRP,0),NFBAS(IFRP,0),NREF,NVEC,
     &     SMAT(I2BASX(IFRP,IFRP)+1),NTBAS(0),NTBAS(0),1,IQBF,
     &     CPROJ,NREF,NVEC,NBF,IQBF,
     &     CREF,NFBAS(IFRP,0),NREF,4,IQCF,
     &     CVEC,NFBAS(IFRP,0),NVEC,4,IQCF,
     &     WORK,LWORK,IPRINT)
      IF(NBF.NE.4) WRITE(LUPRI,'(A)')
     &   '*** WARNING *** ANAPR2: NBF.NE.4 !'
      CALL DSCAL(NREF*NVEC,DM1,CPROJ(1,1,3),1)
C
C     Print section
C     =============  
C
      DO J = 1,NVEC
        FIRST = .TRUE.
        IF(J.GT.NPVEC) THEN
          JJ = NPSH(IFRP)+JVEC(J)
          FMT = MXFORM(EVEC(JJ),16)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +       '* Electronic eigenvalue nr.',JVEC(J),': ',
     &       EVEC(JJ)
        ELSE
          JJ  = NPSH(IFRP)+1+JVEC(J)
          FMT = MXFORM(EVEC(JJ),16)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//')')
     +       '* Positronic eigenvalue nr.',JJ,': ',
     &       EVEC(JJ)
        ENDIF
        CALL PRSYMB(LUPRI,'=',48,0)
C
C       Positronic orbitals
C       -------------------
C
        CPTOT = D0
        DO 20 I = 1,NPREF
          II = NPSH(IFRP)+JREF(I)+1
          CAABS = PYTHAG(CPROJ(I,J,1),CPROJ(I,J,2))
          CBABS = PYTHAG(CPROJ(I,J,3),CPROJ(I,J,4))
          CPABS = PYTHAG(CAABS,CBABS)
          CPTOT = CPTOT + CPABS*CPABS
          IF(CPABS.LT.PROTHR) GOTO 20
          IF(FIRST) THEN
            WRITE(LUPRI,'(4X,A7,4X,A5,11X,A10,10X,A17,13X,A17)')
     &         'Orbital','Total','Eigenvalue',
     &         'Kramers partner 1',
     &         'Kramers partner 2'
            FIRST = .FALSE.
          ENDIF
          BLINE = ' '
          WRITE(BLINE(1:45),100) 'P',II,CPABS,EREF(II)
          IF(CAABS.GT.PROTHR) THEN
            CP1 = CPROJ(I,J,1)/CAABS
            CP2 = CPROJ(I,J,2)/CAABS
            WRITE(BLINE(51:75),110) CAABS,'(',CP1,',',CP2,')'
          ENDIF
          IF(CBABS.GT.PROTHR) THEN
            CP1 = CPROJ(I,J,3)/CBABS
            CP2 = CPROJ(I,J,4)/CBABS
            WRITE(BLINE(81:105),110) CBABS,'(',CP1,',',CP2,')'
          ENDIF
          WRITE(LUPRI,'(A105)') BLINE
   20   CONTINUE
C
C       Electronic orbitals
C       -------------------
C
        CETOT = D0
        DO 30 I = NPREF+1,NREF
          II  = NPSH(IFRP)+JREF(I)
          CAABS = PYTHAG(CPROJ(I,J,1),CPROJ(I,J,2))
          CBABS = PYTHAG(CPROJ(I,J,3),CPROJ(I,J,4))
          CPABS = PYTHAG(CAABS,CBABS)
          CETOT = CETOT + CPABS*CPABS
          IF(CPABS.LT.PROTHR) GOTO 30
          IF(FIRST) THEN
            WRITE(LUPRI,'(4X,A7,4X,A5,11X,A10,10X,A17,13X,A17)')
     &         'Orbital','Total','Eigenvalue',
     &         'Kramers partner 1',
     &         'Kramers partner 2'
            FIRST = .FALSE.
          ENDIF
          BLINE = ' '
          WRITE(BLINE(1:45),100) 'E',JREF(I),CPABS,EREF(II)
          IF(CAABS.GT.PROTHR) THEN
            CP1 = CPROJ(I,J,1)/CAABS
            CP2 = CPROJ(I,J,2)/CAABS
            WRITE(BLINE(51:75),110) CAABS,'(',CP1,',',CP2,')'
          ENDIF
          IF(CBABS.GT.PROTHR) THEN
            CP1 = CPROJ(I,J,3)/CBABS
            CP2 = CPROJ(I,J,4)/CBABS
            WRITE(BLINE(81:105),110) CBABS,'(',CP1,',',CP2,')'
          ENDIF
          WRITE(LUPRI,'(A105)') BLINE
   30   CONTINUE
        CTOT = CPTOT + CETOT
        IF(FIRST) THEN
          WRITE(LUPRI,'(A)') 'No contributions above threshold.'
        ELSE
          WRITE(LUPRI,'(/)')
          WRITE(LUPRI,'(4X,A,1P,E11.4)')
     &   'Electronic contribution:',CETOT,
     &   'Positronic contribution:',CPTOT,
     &   'Total contribution     :',CTOT
        ENDIF
      ENDDO
C
      RETURN
 100  FORMAT(4X,A1,I6,3X,F11.8,5X,1P,E15.8)
 110  FORMAT(F11.8,1X,A1,2(F5.2,A1))
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getplb */
      SUBROUTINE GETPLB(IPRINT)
C*****************************************************************************
C
C     Generate labels for use in population analysis
C
C     Written by Trond Saue May 1996
C     Last revision: tsaue - May 2 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxorb.h>
#include <dcbpop.h>
C
#include <dcbibt.h>
C
      IF(.NOT.LABDEF) THEN
        DO I = 1,NPLAB(1)
          POPLAB(I) = PLABEL(I,ILABDF)
          IPOPLAB(I)  = I
        ENDDO
        NPOPLAB = NPLAB(1)
      ENDIF
      IF(ADDSML) THEN
        JCENT = 0
        JDEG  = 0
        DO I = NPLAB(1)+1,NPLAB(0)
          ICENT = JGET(IATTR(I,ILABDF))
          IDEG  = KGET(IATTR(I,ILABDF))
          IF(ICENT.NE.JCENT.OR.IDEG.NE.JDEG) THEN
            NPOPLAB = NPOPLAB + 1
            POPLAB(NPOPLAB) = PLABEL(I,ILABDF)(3:8)//'_small'
            JCENT = ICENT
            JDEG  = IDEG
          ENDIF
          IPOPLAB(I) = NPOPLAB
        ENDDO
      ELSE
        IF(.NOT.LABDEF) THEN
          DO I = NPLAB(1)+1,NPLAB(0)
            POPLAB(I) = PLABEL(I,ILABDF)
            IPOPLAB(I)  = I
          ENDDO
          NPOPLAB = NPOPLAB + NPLAB(2)
        ENDIF
      ENDIF
C
C     Print section
C     =============
C
      IF(IPRINT.GE.1) THEN
        CALL HEADER('Output from GETPLB: Group labels',-1)
        WRITE(LUPRI,'(A,I5)') '* Total number : ',NPOPLAB
        WRITE(LUPRI,'(6(3X,I3,2X,A12))') (I,POPLAB(I),I=1,NPOPLAB)
        IF(IPRINT.GE.5) THEN
          CALL PRIVEC('IPOPLAB ',IPOPLAB,NPLAB(0))
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anavec */
      SUBROUTINE ANAVEC(CMO,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Print vectors
C
C     Written by Trond Saue May 1996
C     Last revision: tsaue - Jan 9 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbvec.h>
C
      DIMENSION CMO(*),EIG(*),WORK(LWORK)
C
      IF(IPRCMP.LE.0) RETURN
      CALL TITLER('Vector print','*',130)
      CALL PRCOEF(CMO,EIG,VECPRI,IPRCMP,ILABDF,IPRVEC,WORK,LWORK)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopt */
      SUBROUTINE WRPOPT(POPG,LAB,NLAB,IEFF)
C*****************************************************************************
C
C     Output routine for Mulliken gross population
C
C     Written by T.Saue - Oct 11 1995
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DTOL=0.0001D0,D0 = 0.0D0)
      CHARACTER LAB(NLAB)*12,FMT*7
      DIMENSION POPG(NLAB),IEFF(*)
C
C     Find nonzero contributions
C     ==========================
C
      SUMA = D0
      NEFF = 0
      DO I = 1,NLAB
        SUMA = SUMA + POPG(I)
        POPT = ABS(POPG(I))
        IF(POPT.GE.DTOL) THEN
          NEFF = NEFF + 1
          IEFF(NEFF) = I
        ENDIF
      ENDDO
C
C     Write out nonzero contributions
C     ===============================
C
      NROW = (NEFF+1)/8
      IF(MOD((NEFF+1),8).GT.0) NROW = NROW + 1
      NPRI = MIN(NEFF+1,8)
      NFMT = 8 + 15*NPRI
      NPRI = NPRI - 1
      WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
      WRITE(LUPRI,'(/A8,2X,A6,2X,A1,4X,7(A12,3X))')
     &     'Gross   ','Total ','|',(LAB(IEFF(I)),I=1,NPRI)
      WRITE(LUPRI,FMT) ('-',I=1,NFMT)
      WRITE(LUPRI,'(A8,F8.4,2X,A1,4X,7(F8.4,7X))')
     &     ' total  ',SUMA,'|',(POPG(IEFF(I)),I=1,NPRI)
      NOFF = NPRI
      DO J = 2,NROW
        NPRI = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NPRI
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(1X,A12,2X))')
     &       'Gross  |',(LAB(IEFF(NOFF+I)),I=1,NPRI)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        WRITE(LUPRI,'(A8,8(1X,F8.4,6X))')
     &       ' total |',(POPG(IEFF(NOFF+I)),I=1,NPRI)
        NOFF = NOFF + NPRI
      ENDDO
      RETURN
      END
