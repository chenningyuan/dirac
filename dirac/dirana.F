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
C VERSION : $Revision: 1.63.2.3 $
C DATE    : $Date: 2004/09/01 08:08:44 $
C FILE    : dirana.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck paman */
      SUBROUTINE PAMANA(WORK,LWORK)
C***********************************************************************
C
C     Analyze Dirac-Fock wavefunction
C
C     Written by T.Saue Sep 11 1995
C     Last revision : KApril 27 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbbas.h"
C Used from COMMON blocks:
C  CBIANA: DOPOP,DOANA,DOEXP
C  DCBBAS: N2BBASX
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PAMANA')
C
C***********************************************************************
C*****  V E C T O R    P R I N T ***************************************
C***********************************************************************
C
      IF(DOVEC) CALL ANAVEC(WORK,LWORK)
C
C***********************************************************************
C*****  M U L L I K E N    P O P U L A T I O N     A N A L Y S I S  ****
C***********************************************************************
C
      IF(DOPOP) CALL MULPOP(WORK,LWORK)
C
C***********************************************************************
C*****  P L O T   1 - D I M E N S I O N A L     D E N S I T Y   ********
C***********************************************************************
C
      IF(DO1RHO) CALL GT1RHO(WORK,LWORK)
C
C***********************************************************************
C*****  P L O T   3 - D I M E N S I O N A L     D E N S I T Y   ********
C***********************************************************************
C
      IF(DO3RHO) CALL GT3RHO(WORK,LWORK)
C
C***********************************************************************
C*****  P R O J E C T I O N   O N T O    A N O T H E R   S O L U T I O N
C***********************************************************************
C
      IF(DOPRJ) CALL ANAPRJ(WORK,LWORK)
C
C
C***********************************************************************
C*****  P L O T  I N T E G R A T I O N                   ***************
C*****  W E I G H T S   I N   O N E   D I M E N S I O N  ***************
C***********************************************************************
#if defined (MOD_DFT)
      IF(DO1WT) CALL GT1WT(WORK,LWORK)
#endif
C
C***********************************************************************
C
      CALL QEXIT('PAMANA')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpop */
      SUBROUTINE MULPOP(WORK,LWORK)
C***********************************************************************
C
C     Mulliken population analysis
C
C     Written by T.Saue - Sep 11 1995
C     Last revision: April 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
C Used from COMMON blocks:
C  DCBBAS: N2BBASX
C  LABELS: NLAB
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbpop.h"
#include "dgroup.h"
C
      LOGICAL TOBE
      DIMENSION WORK(LWORK)
C
      CALL QENTER('MULPOP')
#include "memint.h"
C
      CALL TITLER('Mulliken population analysis','*',122)
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBEIG,NORBT   ,WORK,KFREE,LFREE)
C
C     Read coefficients, eigenvalues and irrep identification + restart info
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KBEIG),
     &            TOTERG,14,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
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
      CALL MULPO1(WORK(KPMAT),WORK(KSMAT),WORK(KCMO),WORK(KEIG),
     &            WORK(KIPGR),WORK(KPOPG),WORK(KIEFF),WORK(KPGTOT),
     &            WORK(KPNTOT),WORK(KFREE),LFREE)
C
C     Memory deallocation
      CALL MEMREL('MULPOP',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPOP')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'MULPOP: Coefficient file not found !'
      CALL QUIT('MULPOP: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpo1 */
      SUBROUTINE MULPO1(PMAT,SMAT,CMO,EIG,IPGRS,POPG,IEFF,
     &                  PGTOT,PNTOT,WORK,LWORK)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxorb.h"
C Used from COMMON blocks:
C  DCBORB: NSPH
C  DCBBAS:
C  DGROUP: NZ
C
#include "dcbgen.h"
#include "dcbpop.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
C
      DIMENSION PMAT(NPOPLAB,NPOPLAB,2),SMAT(*),CMO(*),EIG(*),IPGRS(*),
     &          POPG(NPOPLAB,2),IEFF(NPOPLAB),
     &          PNTOT(NPOPLAB,NPOPLAB,NFSYM),
     &          PGTOT(NPOPLAB,NFSYM),WORK(LWORK)
C
      CALL QENTER('MULPO1')
#include "memint.h"
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
      CALL GTOVLX(SMAT,SSMTRC,WORK,LWORK)
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
        NVECS = 1
        CALL MEMGET('INTE',KJVEC,NORB(IFRP),WORK,KFREE,LFREE)
        CALL NUMLST(VECPOP(IFRP),WORK(KJVEC),NORB(IFRP),
     &             -NPSH(IFRP),NESH(IFRP),IFRP,NVECS,WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KJVEC),NVECS,IFRP,NPVECS,NEVECS)
        IF(NVECS.EQ.0) GOTO 20
        CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
C
        IF(ILABDF.EQ.1) THEN
C
C       Extract/expand coefficients in AO-basis
C       =======================================
C
          NCBAS = NTBAS(0)
          IF(NVECS.EQ.NORB(IFRP)) THEN
            IF(NBSYM.EQ.1) THEN
              NCBF1 = NTBAS(0)*NVECS*4
              CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
              CALL CFSORT(WORK(KCBF1),CMO(ICMOQ(IFRP)+1),NVECS,IFRP,
     &                    NFBAS(IFRP,0))
            ELSE
              NCBF1 = NTBAS(0)*NVECS*4
              NCBF2 = NFBAS(IFRP,0)*NVECS*4
              CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
              CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
              CALL CFSORT(WORK(KCBF2),CMO(ICMOQ(IFRP)+1),NVECS,IFRP,
     &                    NFBAS(IFRP,0))
              CALL WTSOAO(WORK(KCBF2),WORK(KCBF1),IFRP,
     &                    NFBAS(IFRP,0),NTBAS(0),NVECS*4,IPRPOP)
              CALL MEMREL('MULPO1.AO',WORK,1,KCBF2,KFREE,LFREE)
            ENDIF
          ELSE
            NCBF1 = NTBAS(0)*NVECS*4
            NCBF2 = NTBAS(0)*NVECS*4
            CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
            LCB = NPVECS+NEVECS
            IF(NBSYM.EQ.1) THEN
              CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF2),
     &                    LCB,WORK(KJVEC),NPVECS,NEVECS,
     &                    NFBAS(IFRP,0),NORB(IFRP))
              CALL CFSORT(WORK(KCBF1),WORK(KCBF2),NVECS,IFRP,
     &                    NFBAS(IFRP,0))
            ELSE
              CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF1),
     &                    LCB,WORK(KJVEC),NPVECS,NEVECS,
     &                    NFBAS(IFRP,0),NORB(IFRP))
              CALL CFSORT(WORK(KCBF2),WORK(KCBF1),NVECS,IFRP,
     &                    NFBAS(IFRP,0))
              CALL WTSOAO(WORK(KCBF2),WORK(KCBF1),IFRP,
     &                    NFBAS(IFRP,0),NTBAS(0),NVECS*4,IPRPOP)
            ENDIF
          ENDIF
        ELSE
C
C       Extract/expand coefficients in SO-basis
C       =======================================
C
          NCBAS = NFBAS(IFRP,0)
          IF(NVECS.EQ.NORB(IFRP)) THEN
            NCBF  = NFBAS(IFRP,0)*NVECS
            NCBF1 = NCBF*4
            CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
            CALL CFSORT(WORK(KCBF1),CMO(ICMOQ(IFRP)+1),NVECS,IFRP,
     &                  NFBAS(IFRP,0))
          ELSE
            NCBF  = NFBAS(IFRP,0)*NVECS
            NCBF1 = NCBF*4
            NCBF2 = NCBF*NZ
            LCB = NPVECS+NEVECS
            CALL MEMGET('REAL',KCBF1,NCBF1,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KCBF2,NCBF2,WORK,KFREE,LFREE)
            CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KCBF2),
     &                  LCB,WORK(KJVEC),NPVECS,NEVECS,
     &                  NFBAS(IFRP,0),NORB(IFRP))
            CALL CFSORT(WORK(KCBF1),WORK(KCBF2),NVECS,IFRP,
     &                  NFBAS(IFRP,0))
            CALL MEMREL('MULPO1.SO',WORK,1,KCBF2,KFREE,LFREE)
          ENDIF
        ENDIF
C
C       Perform Mulliken population analysis
C       ====================================
C
        CALL MEMGET('REAL',KOCC,NVECS,WORK,KFREE,LFREE)
        CALL MULPO2(PMAT,IFRP,SMAT,WORK(KCBF1),EIG(IORB(IFRP)+1),
     &              IPGRS,POPG,IEFF,PGTOT,PNTOT,NPVECS,NEVECS,NCBAS,
     &              WORK(KJVEC),WORK(KOCC),WORK(KFREE),LFREE)
        CALL MEMREL('MULPO1.occ',WORK,KOCC,KOCC,KFREE,LFREE)
 20     CONTINUE
        CALL MEMREL('MULPO1',WORK,1,KWORK,KFREE,LFREE)
 10   CONTINUE
C
C     Write out total populations
C     ==========================
C
      IF(NFSYM.GT.1) THEN
        WRITE(LUPRI,'(//A)') '*** Total gross population ***'
        CALL DAXPY(NPOPLAB,D1,PGTOT(1,2),1,PGTOT(1,1),1)        
        CALL WRPOPT(PGTOT(1,1),POPLAB,NPOPLAB,IEFF)                
        IF(DONETP) THEN
          CALL DAXPY(N2LAB,D1,PNTOT(1,1,2),1,PNTOT(1,1,1),1)
          WRITE(LUPRI,'(//A)') '*** Total net population ***'
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mulpo2 */
      SUBROUTINE MULPO2(PMAT,IFRP,SMAT,CMO,EIG,IPGRS,POPG,IEFF,
     &           PGTOT,PNTOT,NPVECS,NEVECS,NCBAS,JVEC,OCC,WORK,LWORK)
C***********************************************************************
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
C     Revised       : Jan 7 1998 - jth LINUX
C     Last revision : Aug 20 2002 - jkp (take fractional occ. into 
C                                        account)
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2 = 2.0D0)
#include "dummy.h"
#include "maxorb.h"
C Used from COMMON blocks:
C  DCBORB: NSPH
C  DCBBAS:
C  DGROUP: NZ
C  DCBDHF: NOPEN
C
#include "dcbgen.h"
#include "dcbpop.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
C 
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0))
#endif
      CHARACTER FMT*6,MXFORM*6
      DIMENSION PMAT(NPOPLAB,NPOPLAB,2),JVEC(*),OCC(*),
     &          CMO(NCBAS,NPVECS+NPVECS,NZ),EIG(*),IPGRS(*),
     &          POPG(NPOPLAB,2),IEFF(NPOPLAB),
     &          PNTOT(NPOPLAB,NPOPLAB,NFSYM),
     &          PGTOT(NPOPLAB,NFSYM),WORK(LWORK)
C
      CALL QENTER('MULPO2')
#include "memint.h"
C
      N2LAB = NPOPLAB*NPOPLAB
      NVEC = NPVECS + NEVECS
      CALL DZERO(OCC,NVEC)
      DO I = 1,NVEC
c
c       Find fractional occupation for each requested orbital string
c
        NVORB =  NISH(IFRP)
        DO IOPEN=1,NOPEN
c     ... open shell
           NVORB =  NVORB + NACSH(IFRP,IOPEN)
           IF (JVEC(I).GT.NISH(IFRP).AND.JVEC(I).LE.NVORB) THEN
              OCC(I) = DF(IOPEN)
              GO TO 100
           ENDIF
        END DO
        IF (JVEC(I) .GT. 0 .AND. JVEC(I) .LE. NISH(IFRP)) THEN
c     ... closed shell
           OCC(I) = DF(0)
        ENDIF
100     CONTINUE
c
        IF(I.GT.NPVECS) THEN
          II = NPSH(IFRP)+JVEC(I)
          FMT = MXFORM(EIG(II),16)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//',A,F6.4,A)')
     +       '* Electronic eigenvalue nr.',JVEC(I),': ',
     &       EIG(II),'       (Occupation : f = ',OCC(I),')'
        ELSE
          II  = NPSH(IFRP)+1+JVEC(I)
          FMT = MXFORM(EIG(II),16)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//',A,F6.4,A)')
     +       '* Positronic eigenvalue nr.',II,': ',
     &       EIG(II),'       (Occupation : f = ',OCC(I),')'
        ENDIF
        CALL PRSYMB(LUPRI,'=',80,0)
C
C       Calculate SO population matrix
C       ==============================
C
        CALL POPMAT(IFRP,I,PMAT,SMAT,CMO,NVEC,
     &              IPGRS,NPOPLAB,NCBAS,WORK(KFREE),LFREE,IPRPOP)
C
C       Calculate gross populations
C       ===========================
C
        CALL GRSPOP(POPG,PMAT,NPOPLAB)
        CALL WRPOPG(POPG,POPLAB,NPOPLAB,IEFF,IPRPOP)
C
C       Multiply fractional occupation on population.
C       =============================================
C
        FACTOR  = D2*OCC(I)
        CALL DAXPY(NPOPLAB,FACTOR,POPG(1,1),1,PGTOT(1,IFRP),1)
        CALL DAXPY(NPOPLAB,FACTOR,POPG(1,2),1,PGTOT(1,IFRP),1)
C
C       Calculate net populations
C       ===========================
C
        IF(DONETP) THEN
          IF (IPRPOP.GT.0) THEN
            CALL WRPOPN(PMAT(1,1,1),POPG(1,1),POPLAB,NPOPLAB,
     &           IEFF,'alpha')
            CALL WRPOPN(PMAT(1,1,2),POPG(1,2),POPLAB,NPOPLAB,
     &           IEFF,'beta ')
          END IF
          CALL DAXPY(N2LAB,D2,PMAT(1,1,1),1,PNTOT(1,1,IFRP),1)
          CALL DAXPY(N2LAB,D2,PMAT(1,1,2),1,PNTOT(1,1,IFRP),1)
        ENDIF
      ENDDO
C
C     Write out total populations
C     ==========================
C
      IF(NFSYM.EQ.1) THEN
        WRITE(LUPRI,'(//A)') '*** Total gross population ***'
      ELSE
        WRITE(LUPRI,'(//A,A3,A)') 
     &    '** Total gross population of fermion ircop ',
     &    FREP(IFRP),' **'
      ENDIF
      CALL WRPOPT(PGTOT(1,IFRP),POPLAB,NPOPLAB,IEFF)                
      IF(DONETP) THEN
        IF(NFSYM.EQ.1) THEN
          WRITE(LUPRI,'(//A)') '*** Total net population ***'
        ELSE
          WRITE(LUPRI,'(//A,A3,A)') 
     &      '** Total net population of fermion ircop ',
     &      FREP(IFRP),' **'
        ENDIF
        CALL WRPOPN(PNTOT(1,1,IFRP),PGTOT(1,IFRP),
     &              POPLAB,NPOPLAB,IEFF,'total')
      ENDIF
      CALL MEMREL('MULPO2.POPMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MULPO2')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popmat */
      SUBROUTINE POPMAT(IFRP,IVEC,PMAT,SMAT,CMO,NVECS,IGLAB,NGLAB,
     &                  NCBAS,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Construct population matrix
C         M(i,j) = Re(CONJ[CMO(i,IVEC)]*SMAT(i,j)*CMO(j,IVEC))
C     The final population matrix is compressed based on a labeling system.
C
C     Written by T.Saue Oct 9 1995
C     Last revision: Oct 9 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C Used from COMMON blocks
C  DCBAS: N2BBASX
#include "dcbbas.h"
      DIMENSION PMAT(*),SMAT(*),CMO(*),IGLAB(*),WORK(LWORK)
C
      CALL QENTER('POPMAT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KPBUF,2*N2BBASX,WORK,KFREE,LFREE)
C
      CALL POPMA1(IFRP,IVEC,PMAT,SMAT,CMO,NVECS,WORK(KPBUF),
     &     IGLAB,NGLAB,NCBAS,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('POPMAT.POPMA1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('POPMAT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popma1 */
      SUBROUTINE POPMA1(IFRP,IVEC,PMAT,SMAT,CMO,NVECS,PBUF,
     &                  IGLAB,NGLAB,NCBAS,IPRINT)
C***********************************************************************
C
C     Construct population matrix
C         M(i,j) = Re(CONJ[CMO(i,IVEC)]*SMAT(i,j)*CMO(j,IVEC))
C     The final population matrix is compressed based on a labeling system.
C
C     Written by T.Saue Oct 9 1995
C     Last revision: Jan 7 1998 - jth LINUX
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
#include "maxorb.h"
C
C Used from COMMON blocks
C  DCBBAS: NTBAS,NFBAS
C  DCBORB: NORB
C  DGROUP: NZ
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbpop.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST),PBUF(NBAST,NBAST,*)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0)),PBUF(NTBAS(0),NTBAS(0),*)
#endif
      DIMENSION PMAT(NGLAB,NGLAB,2),CMO(NCBAS,NVECS,4),IGLAB(*)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grspop */
      SUBROUTINE GRSPOP(POPG,PMAT,NLAB)
C***********************************************************************
C
C     From Mulliken population matrix generate gross populations for
C     2-spinor.
C
C     Written by T.Saue - July 24 1994
C     LAST VERSION: Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION POPG(NLAB,2),PMAT(NLAB,NLAB,2)
C
      DO 10 J = 1,NLAB
        POPG(J,1) = DSUM(NLAB,PMAT(1,J,1),1)
        POPG(J,2) = DSUM(NLAB,PMAT(1,J,2),1)
   10 CONTINUE
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopg */
      SUBROUTINE WRPOPG(POPG,LAB,NLAB,IEFF,IPRINT)
C***********************************************************************
C
C     Output routine for Mulliken gross population
C
C     Written by T.Saue - Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
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
        WRITE(LUPRI,'(/A,F8.5)') '* Gross populations greater than',DTOL
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
C
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,'(/A)') '* Gross populations (all):'
        DO I = 1,NLAB
          WRITE(LUPRI,'(A12,1P,2E16.8)') LAB(I),POPG(I,1),POPG(I,2)
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck netpop */
      SUBROUTINE NETPOP(POPN,PMAT,NLAB)
C***********************************************************************
C
C     Mulliken net populations (upper triangular columnwise packed)
C
C     Written  by T.Saue Oct 11 1995
C     Last revision: Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopn */
      SUBROUTINE WRPOPN(PMAT,POPG,LAB,NLAB,IEFF,COMP)
C***********************************************************************
C
C     Output routine for Mulliken net populations
C
C     Written by T.Saue Oct 11 1995
C     Last revision Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anaprj */
      SUBROUTINE ANAPRJ(WORK,LWORK)
C***********************************************************************
C
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0,DM1=-1.0D0,D1=1.0D0)
C
#include "dcbgen.h"
#include "dcbprj.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
C
      LOGICAL TOBE
      DIMENSION WORK(LWORK)
      DIMENSION GROSS(2,0:2,MAXREF),POL(2)
      DIMENSION NSTR(2,0:2,0:MAXREF),KVEC(2,0:MAXREF),
     &          KE(2,2),KQ(2,2),NVECS(2,0:2),NQ(2),NE(2),
     &          KBV(2)
C
#include "memint.h"
      KFRSAV = KFREE
      CALL TITLER('Projection analysis','*',122)
C
C     Find total numbers of reference orbitals
C
      DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        DO I = 1,NREFS
          CALL MEMGET('INTE',KVEC(IFRP,I),NORB(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL NUMLST(VECREF(IFRP,I),WORK(KVEC(IFRP,I)),
     &              NORB(IFRP),-NPSH(IFRP),NESH(IFRP),
     &              IFRP,NSTR(IFRP,0,I),WORK(KFREE),LFREE)
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &                NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NVECS(IFRP,1) = NVECS(IFRP,1) + NSTR(IFRP,0,I)
          CALL MEMREL('PRJANA',WORK,KVEC(IFRP,I),KVEC(IFRP,I),
     &                KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO
C
C     Find number of MO-orbitals to analyse
C
      DO IFRP = 1,NFSYM
        CALL MEMGET('INTE',KVEC(IFRP,0),NORB(IFRP),WORK,KFREE,LFREE)
        NSTR(IFRP,0,0) = 1
        CALL NUMLST(VECPRJ(IFRP),WORK(KVEC(IFRP,0)),
     &              NORB(IFRP),-NPSH(IFRP),NESH(IFRP),
     &              IFRP,NSTR(IFRP,0,0),WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KVEC(IFRP,0)),NSTR(IFRP,0,0),IFRP,
     &              NSTR(IFRP,2,0),NSTR(IFRP,1,0))
        NSTR(IFRP,0,0) = NSTR(IFRP,1,0) + NSTR(IFRP,2,0)
        NVECS(IFRP,2) = NSTR(IFRP,0,0)
        CALL MEMREL('PRJANA',WORK,KVEC(IFRP,0),KVEC(IFRP,0),
     &              KFREE,LFREE)
        CALL MEMGET('INTE',KVEC(IFRP,0),NSTR(IFRP,0,0),
     &                WORK,KFREE,LFREE)
      ENDDO
C
C     Selected reference eigenvalues
C
      NDIM = 0
      DO IFRP = 1, NFSYM
        NVECS(IFRP,0) = NVECS(IFRP,1) + NVECS(IFRP,2)
        NDIM = NDIM + NVECS(IFRP,1)
      ENDDO
      CALL MEMGET('REAL',KE(1,1),NDIM,WORK,KFREE,LFREE)
C
C     Selected molecular eigenvalues
C
      NDIM = 0
      DO IFRP = 1, NFSYM
        NDIM = NDIM + NVECS(IFRP,2)
      ENDDO
      CALL MEMGET('REAL',KE(1,2),NDIM,WORK,KFREE,LFREE)
C
C     Calculate dimensions of coefficient arrays:
C      reference orbitals ............... KQ(IFRP,1) 
C      are followed by molecular ones ... KQ(IFRP,2)
C
      NCDIM = 0
      DO IFRP = 1,NFSYM
        NQ(IFRP) = NFBAS(IFRP,0)*NVECS(IFRP,0)*NZ
        NCDIM    = NCDIM + NQ(IFRP)
      ENDDO       
      CALL MEMGET('REAL',KQ(1,1),NCDIM,WORK,KFREE,LFREE)
      KQ(1,2)   = KQ(1,1) + NFBAS(1,0)*NVECS(1,1)
      IF(NFSYM.EQ.2) THEN
        KE(2,1) = KE(1,1) + NVECS(1,1)
        KE(2,2) = KE(1,2) + NVECS(1,2)
        KQ(2,1) = KQ(1,1) + NQ(1)
        KQ(2,2) = KQ(2,1) + NFBAS(2,0)*NVECS(2,1)
      ELSE
        KE(2,1) = KFREE
        KE(2,2) = KFREE
        KQ(2,1) = KFREE
        KQ(2,2) = KFREE
      ENDIF
C
C     COEFFICIENT SELECTION
C     =====================
C
      ICFS = 6
      IF(IPRPRJ.GE.4) ICFS = 7
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT   ,WORK,KFREE,LFREE)
      KIBE = KFREE
C
C     Select reference orbitals and eigenvalues
C
      DO IFRP = 1,NFSYM
        NQ(IFRP) = KQ(IFRP,1)
        NE(IFRP) = KE(IFRP,1)
      ENDDO
      DO I = 1,NREFS
        INQUIRE(FILE=REFFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
           WRITE(LUPRI,'(A,A6,A)') 
     &     'Reference coefficient file ',REFFIL(I),' not present'
           CALL QUIT('PRJANA: No reference coefficients !')
        ENDIF
        CALL OPNFIL(LUCOEF,REFFIL(I),'OLD','PRJANA')
        IF(IPRPRJ.GE.4) WRITE(LUPRI,'(A,A6)') 
     &       '* Coefficients for reference ',REFFIL(I)
        IF(OWNBAS) THEN
          CALL SELOWN(7,WORK(KCMO),WORK(KEIG),WORK(KIBE),
     &                I,LUCOEF,REFFIL(I),NQ,NE,IDUM,
     &                NVECS(1,0),KVEC(1,I),
     &                NSTR(1,0,I),WORK,KFREE,LFREE)
          DO IFRP = 1, NFSYM
            IF(NSTR(IFRP,0,I).GT.0) THEN
              NQ(IFRP) = NQ(IFRP) + NSTR(IFRP,0,I)*NFBAS(IFRP,0)
              NE(IFRP) = NE(IFRP) + NSTR(IFRP,0,I)
            ENDIF
         ENDDO
        ELSE
          CALL REACMO(LUCOEF,REFFIL(I),WORK(KCMO),WORK(KEIG),
     &                IDUM,TOTERG,ICFS,WORK,KFREE,LFREE)
          DO IFRP = 1, NFSYM
            IF(NSTR(IFRP,0,I).GT.0) THEN
              CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(NQ(IFRP)),
     &                       NVECS(IFRP,0),WORK(KVEC(IFRP,I)),
     &                       NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                       NFBAS(IFRP,0),NORB(IFRP))
              CALL SELEIG(WORK(KEIG+IORB(IFRP)),IFRP,WORK(NE(IFRP)),
     &             WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I))
              NQ(IFRP) = NQ(IFRP) + NSTR(IFRP,0,I)*NFBAS(IFRP,0)
              NE(IFRP) = NE(IFRP) + NSTR(IFRP,0,I)
            ENDIF
          ENDDO
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
      ENDDO
C
C     Get molecular coefficients and eigenvalues
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Molecular coefficient file not present'
         CALL QUIT('PRJANA: No coefficients !')
      ENDIF
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRJANA')
      IF(IPRPRJ.GE.4) WRITE(LUPRI,'(A)') '* Molecular coefficients:'
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),IDUM,
     &            TOTERG,ICFS,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
      DO IFRP = 1, NFSYM
        IF(NVECS(IFRP,2).GT.0) THEN
          CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,2)),
     &                   NVECS(IFRP,0),WORK(KVEC(IFRP,0)),
     &                   NSTR(IFRP,2,0),NSTR(IFRP,1,0),
     &                   NFBAS(IFRP,0),NORB(IFRP))
          CALL SELEIG(WORK(KEIG+IORB(IFRP)),IFRP,WORK(KE(IFRP,2)),
     &         WORK(KVEC(IFRP,0)),NSTR(IFRP,2,0),NSTR(IFRP,1,0))
        ENDIF
      ENDDO
C
C     Throw away the full coefficients and eigenvalues
C     ================================================
C
      CALL MEMREL('PRJANA.cfs',WORK,KCMO,KCMO,KFREE,LFREE)
C
C     Do projection analysis
C     ======================
C
C     The polarization contribution will be appended to the solution vector
C
      NBDIM = 0
      DO IFRP = 1,NFSYM
        NBDIM = NBDIM + NVECS(IFRP,0)*NVECS(IFRP,2)*NZ
      ENDDO       
      CALL MEMGET('REAL',KBV(1),NBDIM,WORK,KFREE,LFREE)
      IF(NFSYM.EQ.2) THEN
        KBV(2)   = KBV(1) + NVECS(1,0)*NVECS(1,2)*NZ
      ELSE
        KBV(2)   = KFREE
      ENDIF
      CALL MEMGET('REAL',KSMAT,N2BAST,WORK,KFREE,LFREE)
      CALL GTOVLT(WORK(KSMAT),SSMTRC,WORK(KFREE),LFREE)
      IS = KSMAT
      POLT = D0
      DO IFRP = 1,NFSYM
        CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
        KBUF   = KFREE
        NMOL   = NVECS(IFRP,2)
        NREF   = NVECS(IFRP,1)
        N2REF  = NREF*NREF
        CALL MEMGET('REAL',KGVEC,(NREF+1)*2,WORK,KFREE,LFREE)
        IF(NREFS.GT.1) THEN
          CALL MEMGET('REAL',KAMAT,N2REF*NZ,WORK,KFREE,LFREE)
        ELSE
          KAMAT = KFREE
        ENDIF
        CALL MEMGET('INTE',KJREF,NREF*3,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIPIV,NREF  ,WORK,KFREE,LFREE)
        IVEC = 0
        DO I = 1,NREFS
          CALL PRJINF(WORK(KJREF),IFRP,I,IVEC,NSTR,WORK(KVEC(IFRP,I)))
          IVEC = IVEC + NSTR(IFRP,0,I)
        ENDDO
        CALL PRJAN1(IFRP,NVECS(IFRP,2),NVECS(IFRP,1),NVECS(IFRP,0),
     &       KQ,KE,WORK(KAMAT),WORK(KBV(IFRP)),
     &       WORK(KGVEC),WORK(KVEC(IFRP,0)),WORK(KJREF),NSTR,
     &       WORK(IS),WORK(KIPIV),GROSS,POL(IFRP),WORK,KFREE,LFREE)
        POLT = POLT + POL(IFRP)
C
C       Calculate polarization contribution
C       ===================================
C
        CALL QGEMM(NFBAS(IFRP,0),NVECS(IFRP,2),NVECS(IFRP,1),DM1,
     &             'N','N',IPQTOQ(1,0),
     &             WORK(KQ(IFRP,1)),NFBAS(IFRP,0),NVECS(IFRP,0),NZ,
     &             'N','N',IPQTOQ(1,0),
     &             WORK(KBV(IFRP)),NVECS(IFRP,0),NVECS(IFRP,2),NZ,
     &             D1,IPQTOQ(1,0),
     &             WORK(KQ(IFRP,2)),NFBAS(IFRP,0),NVECS(IFRP,0),NZ)
        IF(IPRPRJ.GE.5) THEN
          CALL HEADER('Polarization vectors',0)
          CALL PRIPOL(IFRP,WORK(KQ(IFRP,2)),WORK(KE(IFRP,2)),
     &                NVECS(IFRP,2),NFBAS(IFRP,0),NVECS(IFRP,0))
        ENDIF      
C       Add polarization contribution to expansion coefficients
        CALL DUNIT2(WORK(KBV(IFRP)+NVECS(IFRP,1)),NVECS(IFRP,2),
     &              NVECS(IFRP,0),NVECS(IFRP,2),NZ)
        IS = IS + N2BAS(IFRP)
        CALL MEMREL('PRJANA.loop',WORK,KBUF,KBUF,KFREE,LFREE)
      ENDDO
C
      IF(NFSYM.EQ.2) THEN
        CALL HEADER('Total gross contributions:',-1)
        GTOT = D0
        DO I = 1,NREFS
          WRITE(LUPRI,'(3X,A6,8X,3X,F7.4,2(3X,A4,F7.4))')
     &        REFFIL(I),(GROSS(1,0,I)+GROSS(2,0,I)),
     &        'E - '   ,(GROSS(1,2,I)+GROSS(2,2,I)),
     &        'P - '   ,(GROSS(1,1,I)+GROSS(2,1,I))
          GTOT = GTOT + GROSS(1,0,I)+GROSS(2,0,I)
        ENDDO
        WRITE(LUPRI,'(3X,A,F10.4)') 'Polarization: ',POLT
      ENDIF
C
C     Dump reference and projection coefficients to file
C     ==================================================
C
      CALL OPNFIL(LUCOEF,'DFPRJC','UNKNOWN','PRJANA')
      WRITE(LUCOEF) NREFS
      DO IFRP = 1,NFSYM
C       Dimensions....NREF..........NMOL.........
        WRITE(LUCOEF) NVECS(IFRP,1),NVECS(IFRP,2),
     &                (NSTR(IFRP,0,I),I=1,NREFS)
      ENDDO
C.....Reference + polarization coefficients
      CALL WRITT(LUCOEF,NCDIM,WORK(KQ(1,1)))
C.....Projection coefficients
      CALL WRITT(LUCOEF,NBDIM,WORK(KBV(1)))
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL MEMREL('PRJANA',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjan1 */
      SUBROUTINE PRJAN1(IFRP,NMOL,NREF,NTOT,KQ,KE,AMAT,BVEC,GVEC,
     &           JMOL,JREF,NSTR,SMAT,IPIVOT,GROSS,POL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Projection analysis for fermion ircop IFRP
C
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,DC=1.0D-2,D2=2.0D0)
C
#include "dcbprj.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
C
      LOGICAL DOPOP
      DIMENSION GROSS(2,0:2,MAXREF)
      DIMENSION NSTR(2,0:2,0:MAXREF),JREF(3,*),JMOL(*),
     &          KE(2,2),KQ(2,2),AMAT(NREF,NREF,NZ),
     &          BVEC(NTOT,NMOL,NZ),SMAT(*),GVEC(*),IPIVOT(*),WORK(*)
C
      KFRSAV = KFREE
      DOPOP = NREFS.GT.1
      PTOL = PROTHR*DC
C
C     Make projection vector  BVEC = <\psi^B_j1\psi^{MO}_i>
C
      CALL QTRANS('AOMO','S',D0,
     &            NFBAS(IFRP,0),NFBAS(IFRP,0),NREF,NMOL,
     &            SMAT,NFBAS(IFRP,0),NFBAS(IFRP,0),1,IPQTOQ(1,0),
     &            BVEC,NTOT,NMOL,NZ,IPQTOQ(1,0),
     &            WORK(KQ(IFRP,1)),NFBAS(IFRP,0),NTOT,NZ,IPQTOQ(1,0),
     &            WORK(KQ(IFRP,2)),NFBAS(IFRP,0),NTOT,NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRPRJ)
      CALL MEMGET('REAL',KOCC,NMOL,WORK,KFREE,LFREE)
      IF(NREFS.GT.1) THEN
C
C       If requested , print raw projection coefficients
C
        IF(IPRPRJ.GE.3) THEN
          WRITE(LUPRI,'(A,A3)') 
     &     '* Raw projection coefficients for fermion ircop ',FREP(IFRP)
          CALL ANAPR2(IFRP,.FALSE.,NREF,NMOL,NTOT,SMAT,GVEC,AMAT,BVEC,
     &            JMOL,WORK(KE(IFRP,2)),NSTR,
     &            JREF,WORK(KE(IFRP,1)),GROSS,WORK(KOCC),POL)
        ENDIF
C
C       Make overlap matrix AMAT = <\psi^B_j1\psi^A_i>
C
        CALL QTRANS('AOMO','S',D0,NFBAS(IFRP,0),NFBAS(IFRP,0),
     &           NREF,NREF,
     &           SMAT,NFBAS(IFRP,0),NFBAS(IFRP,0),1,IPQTOQ(1,0),
     &           AMAT,NREF,NREF,NZ,IPQTOQ(1,0),
     &           WORK(KQ(IFRP,1)),NFBAS(IFRP,0),NTOT,NZ,IPQTOQ(1,0),
     &           WORK(KQ(IFRP,1)),NFBAS(IFRP,0),NTOT,NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRPRJ)
C
C       Print overlaps
C
        JOFF = 1
        DO J = 1,NREFS
          IOFF = 1
          DO I = 1,(J-1)
            WRITE(LUPRI,'(A,2(3X,A6))') '*Overlaps :',
     &           REFFIL(I),REFFIL(J)
            CALL PRQMAT(AMAT(IOFF,JOFF,1),
     &           NSTR(IFRP,0,I),NSTR(IFRP,0,J),NREF,NREF,NZ,
     &           IPQTOQ(1,0),LUPRI)
            IOFF = IOFF + NSTR(IFRP,0,I)
          ENDDO
          JOFF = JOFF + NSTR(IFRP,0,J)
        ENDDO
C
C       Solve linear system by Cholesky decomposition
C       with full pivoting
C        --> notice that overlap matrix SMAT is destroyed !
C
        JOB = 1
        CALL QCHOLD(AMAT,NREF,NZ,NREF,NREF,SMAT,PTOL,NEFF,JOB,IPIVOT)
        CALL QCHOLS(AMAT,NREF,NEFF,NMOL,NZ,NREF,NREF,SMAT,BVEC,
     &       NTOT,NMOL,BVEC,NTOT,NMOL,JOB,IPIVOT,SMAT(1+NREF))
        IF(NEFF.LT.NREF) THEN
          WRITE(LUPRI,'(A,A)') '* WARNING: ',
     &     ' linear dependencies detected in Cholesky decomposition'
          WRITE(LUPRI,'(11X,A,I5,A,I5)') 
     &     ' Reduced dimensionality: ',NREF,' --> ',NEFF
        ENDIF
C
C       Restore overlap matrix
C
        CALL QHMRST(AMAT,NREF,NZ,NREF,NREF)
      ENDIF
C
C     Population analysis
C
      CALL ANAPR2(IFRP,DOPOP,NREF,NMOL,NTOT,SMAT,GVEC,AMAT,BVEC,
     &            JMOL,WORK(KE(IFRP,2)),NSTR,
     &            JREF,WORK(KE(IFRP,1)),GROSS,WORK(KOCC),POL)
C
      CALL MEMREL('PRJAN1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjinf */
      SUBROUTINE PRJINF(JREF,IFRP,I,IVEC,NSTR,KREF)
C***********************************************************************
C
C     Make pointer vector of reference orbital information for use 
C     in projection analysis.
C
C     Wrintten by T.Saue March 31 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbprj.h"
      DIMENSION JREF(3,*),NSTR(2,0:2,0:MAXREF),KREF(*)
C
C     Positrons
C
      IP = 2      
      DO J = 1,NSTR(IFRP,2,I)
        JREF(1,IVEC+J) = NPSH(IFRP)+1+KREF(J)
        JREF(2,IVEC+J) = IP
        JREF(3,IVEC+J) = I
      ENDDO
C
C     Electrons
C
      IP = 1
      DO J = NSTR(IFRP,2,I)+1,NSTR(IFRP,0,I)
        JREF(1,IVEC+J) = KREF(J)
        JREF(2,IVEC+J) = IP
        JREF(3,IVEC+J) = I
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anapr2 */
      SUBROUTINE ANAPR2(IFRP,DOPOP,NREF,NMOL,NTOT,PMAT,GVEC,AMAT,BVEC,
     &                  JMOL,EMOL,NSTR,JREF,EREF,GROSS,OCC,POL)
C***********************************************************************
C
C
C     Project molecular solution down onto reference solution
C     (e.g. atomic) in the same basis
C
C         Y(mol)_i = Y(ref)_jc_ji
C
C     Projection coeffients:  c_ji = <Y(ref)_j|Y(mol)_i>
C     PROTHR is threshold for print of projection coeffiencients.
C
C     Written by Trond Saue March 13 1996
C     Last revision Jan 7 1998 - jth NFBAS(IFRP,0) is now parameter
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbprj.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbdhf.h"
C
      LOGICAL DOPOP
      CHARACTER FMT*6,MXFORM*6,OTYP*1
      DIMENSION GROSS(2,0:2,MAXREF)
      DIMENSION PMAT(NREF,NREF),GVEC(NREF+1,2),AMAT(NREF,NREF,NZ),
     &          BVEC(NTOT,NMOL,NZ),JMOL(*),EMOL(*),
     &          NSTR(2,0:2,0:MAXREF),JREF(3,*),EREF(*),OCC(*)
      DIMENSION P(4),OTYP(2)
      DATA OTYP/'E','P'/
      SAVE OTYP
C
C     Initialization
C
      N2REF = NREF*NREF
      CALL DZERO(GVEC,(NREF+1))
      DO I = 1,4
        P(I) = D0
      ENDDO
C
C     Loop over molecular orbitals
C
      CALL DZERO(OCC,NMOL)
      DO 10 K = 1,NMOL
        FMT = MXFORM(EMOL(K),16)
        IF(K.GT.NSTR(IFRP,2,0)) THEN
C...    POSITIVE-ENERGY (electronic) ORBITAL
c
c         Find fractional occupation for each requested orbital string
c
          NVORB =  NISH(IFRP)
          DO IOPEN=1,NOPEN
c       ... open shell
             NVORB =  NVORB + NACSH(IFRP,IOPEN)
             IF (JMOL(K).GT.NISH(IFRP).AND.JMOL(K).LE.NVORB) THEN
               OCC(K) = DF(IOPEN)
               GO TO 11
             ENDIF
          END DO
          IF (JMOL(K).GT.0.AND.JMOL(K).LE.NISH(IFRP)) THEN
c       ... closed shell
             OCC(K) = DF(0)
          ENDIF
 11       CONTINUE
          WRITE(LUPRI,'(/A,I3,A,'//FMT//',A,F6.4,A)')
     +       '* Electronic eigenvalue nr.',JMOL(K),': ',
     &       EMOL(K),'       (Occupation : f = ',OCC(K),')'
        ELSE
C...    NEGATIVE-ENERGY (positronic) ORBITAL
          KK  = NPSH(IFRP)+1+JMOL(K)
          WRITE(LUPRI,'(/A,I3,A,'//FMT//',A,F6.4,A)')
     +       '* Positronic eigenvalue nr.',KK,': ',
     &       EMOL(K),'       (Occupation : f = ',OCC(K),')'
        ENDIF
        CALL PRSYMB(LUPRI,'=',48,0)
        WRITE(LUPRI,100) 'Orbital','Total','Eigenvalue',
     &       'Kramers partner 1','Kramers partner 2'
C
C       Print coefficients
C
        NPRI = 0
C
        DO 20 I = 1,NREF
          CW = D0
          DO IZ = 1,NZ
            CW = CW + BVEC(I,K,IZ)*BVEC(I,K,IZ)
          ENDDO
          CW = SQRT(CW)
          IF(CW.LT.PROTHR) GOTO 20
          NPRI = NPRI + 1
          KORB = JREF(1,I)
          KTYP = JREF(2,I)
          KREF = JREF(3,I)
          DO IZ = 1,NZ
            P(IZ) = BVEC(I,K,IZ)/CW
          ENDDO
C
          WRITE(LUPRI,110) REFFIL(KREF),OTYP(KTYP),KORB,CW,EREF(I),
     &     '(',P(1),',',P(2),')','(',P(3),',',P(4),')'
   20   CONTINUE
        IF(NPRI.EQ.0) THEN
          WRITE(LUPRI,'(A)') ' --> No contributions !'	
          GOTO 10
        ENDIF
C
C       Population analysis
C
        IF(DOPOP) THEN
          CALL DZERO(PMAT,N2REF)
C
C         Calculate net population matrix
C
          DO JZ = 1,NZ
            DO IZ = 1,NZ
              KZ  = IQMULT(IZ,JZ,1)
              FAC = IQSIGN(IZ,2,1)*IQSIGN(JZ,1,1)*IQPHASE(IZ,JZ,KZ)
              DO J = 1,NREF
                DO I = 1,NREF
                  PMAT(I,J) = PMAT(I,J) + FAC*
     &                        BVEC(I,K,IZ)*AMAT(I,J,KZ)*BVEC(J,K,JZ)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C
C         Calculate gross population vector
C
          IF(LWGPOP) THEN
            DO J = 1,NREF
              GVEC(J,2) = PMAT(J,J)
              DO I = 1,(J-1)
                CW = PMAT(I,I) + PMAT(J,J)
                IF(CW.GT.D0) THEN
                  GVEC(J,2) = GVEC(J,2) + D2*(PMAT(J,J)/CW)*PMAT(I,J)
                  GVEC(I,2) = GVEC(I,2) + D2*(PMAT(I,I)/CW)*PMAT(I,J)
                ENDIF              
              ENDDO
            ENDDO
          ELSE
            DO J = 1,NREF
              GVEC(J,2) = DSUM(NREF,PMAT(1,J),1)
            ENDDO
          ENDIF
          IF(IPRPRJ.GE.2) THEN
            WRITE(LUPRI,'(A)') '* Net population'
            CALL OUTPUT(PMAT,1,NREF,1,NREF,NREF,NREF,1,LUPRI)
            WRITE(LUPRI,'(A)') '* Gross population'
            CALL OUTPUT(GVEC(1,2),1,NREF,1,1,NREF,1,1,LUPRI)
          ENDIF
        ELSE
C
C         Calculate gross population vector
C
          DO J = 1,NREF
            GVEC(J,2) = BVEC(J,K,1)*BVEC(J,K,1)
            DO JZ = 2,NZ
              GVEC(J,2) = GVEC(J,2) + BVEC(J,K,JZ)*BVEC(J,K,JZ)
            ENDDO
          ENDDO
          IF(IPRPRJ.GE.2) THEN
            WRITE(LUPRI,'(A)') '* Gross population'
            CALL OUTPUT(GVEC(1,2),1,NREF,1,1,NREF,1,1,LUPRI)
          ENDIF
        ENDIF
        WRITE(LUPRI,'(A)') '* Gross contributions:'
        CALL GROSSC(GROSS,NREFS,IFRP,NREF,D1,GVEC(1,2),NSTR,REFFIL)
        FACTOR  = D2*OCC(K)
        CALL DAXPY((NREF+1),FACTOR,GVEC(1,2),1,GVEC(1,1),1)
  10  CONTINUE
      CALL PRSYMB(LUPRI,'=',48,0)
      WRITE(LUPRI,'(/A)') 
     &    '* Total reference orbital contributions:'
      DO I = 1,NREF
        KORB = JREF(1,I)
        KTYP = JREF(2,I)
        KREF = JREF(3,I)
        WRITE(LUPRI,110)
     &   REFFIL(KREF),OTYP(KTYP),KORB,GVEC(I,1),EREF(I)
      ENDDO
      WRITE(LUPRI,'(/A)') '* Total gross contributions:'
      TOTM = D2*DSUM(NMOL,OCC,1)
      CALL GROSSC(GROSS,NREFS,IFRP,NREF,TOTM,GVEC(1,1),NSTR,REFFIL)
      POL = GVEC(NMOL+1,1)
C
      RETURN
 100  FORMAT(4X,A7,12X,A5,11X,A10,4X,2(3X,A17))
 110  FORMAT(4X,A6,2X,A1,I6,3X,F11.5,5X,E15.8,
     &       2(3X,A1,F5.2,A1,F5.2,A1,3X))
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getplb */
      SUBROUTINE GETPLB(IPRINT)
C***********************************************************************
C
C     Generate labels for use in population analysis
C
C     Written by Trond Saue May 1996
C     Last revision: tsaue - May 2 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "dcbpop.h"
C
#include "dcbibt.h"
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
      IF(IPRINT.GE.2) THEN
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck anavec */
      SUBROUTINE ANAVEC(WORK,LWORK)
C***********************************************************************
C
C     Print vectors
C
C     Written by Trond Saue May 1996
C     Last revision: tsaue - Jan 9 1997
C     Vector analysis for KRMC added: T. Fleig, Oct 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "orbene.h"
#include "dcbgen.h"
#include "dcbvec.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "consts.h"
C
      LOGICAL TOBE
      DIMENSION WORK(*)
      DIMENSION IOFPFS(2,2)
#include "memint.h"
C
      NTESTL = 000
C
      IF(IPRCMP.LE.0) RETURN
      CALL TITLER('Vector print','*',130)
#ifdef MOD_KRMC
C
C     Check if MCSCF coefficients are on file
C
      INQUIRE (FILE ='KRMCSCF',EXIST=TOBE)
C
      IF(.NOT.TOBE) GOTO 990
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
C
      CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','PAMANA')
      IOPT = 1       ! Read from NEWORB label
      JRDMO = -1
      CALL RREADMO(WORK(KCMO),JRDMO,IOPT,LUKRMC)
      CLOSE(LUKRMC,STATUS='KEEP')
C
      IF (NTESTL.GE.20) THEN
C
C     Print all coefficients (including positronic and S component)
C
        DO I = 1, NFSYM
           WRITE(LUPRI,'(A,I2)') ' (ROPTST)  Coefficients, irrep ',I
           CALL PRQMAT(WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &          NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        END DO
      END IF
C
      CALL HEADER('Coefficients from KRMCSCF',2)
      CALL PRCOEF(WORK(KCMO),ORBENE,VECPRI,IPRCMP,ILABDF,
     &            IPRVEC,WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('ANAVEC',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
#endif /* MOD_KRMC */
C
C     No MCSCF, then check if RHF on DFCOEF :
C
 990  CONTINUE
      INQUIRE (FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBEIG,NORBT   ,WORK,KFREE,LFREE)
C
C     Read coefficients
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KBEIG),
     &            TOTERG,14,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
      IF (NTESTL.GE.10) THEN
C
C     Check orbital energies.
C
        DO J=1,NORBT,1
          WRITE(6,*) 'ORBENE(',J,') = ',WORK(KEIG+J-1)
        END DO
      END IF
C
      CALL HEADER('Coefficients from DFCOEF',2)
      CALL PRCOEF(WORK(KCMO),WORK(KEIG),VECPRI,IPRCMP,ILABDF,IPRVEC,
     &            WORK(KFREE),LFREE)
C     Memory deallocation
      CALL MEMREL('ANAVEC',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(//A)') 'ANAVEC: Coefficient file not found !'
      CALL QUIT('ANAVEC: Coefficient file not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrpopt */
      SUBROUTINE WRPOPT(POPG,LAB,NLAB,IEFF)
C***********************************************************************
C
C     Output routine for Mulliken gross population
C
C     Written by T.Saue - Oct 11 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
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
      WRITE(LUPRI,'(/A8,3X,A6,3X,A1,4X,7(A12,3X))')
     &     'Gross   ','Total ','|',(LAB(IEFF(I)),I=1,NPRI)
      WRITE(LUPRI,FMT) ('-',I=1,NFMT)
      WRITE(LUPRI,'(A8,F10.5,2X,A1,4X,7(F9.5,6X))')
     &     ' total  ',SUMA,'|',(POPG(IEFF(I)),I=1,NPRI)
      NOFF = NPRI
      DO J = 2,NROW
        NPRI = MIN(NEFF-NOFF,8)
        NFMT = 8 + 15*NPRI
        WRITE(FMT,'(A1,I3,A3)') '(',NFMT,'A1)'
        WRITE(LUPRI,'(/A8,8(1X,A12,2X))')
     &       'Gross  |',(LAB(IEFF(NOFF+I)),I=1,NPRI)
        WRITE(LUPRI,FMT) ('-',I=1,NFMT)
        WRITE(LUPRI,'(A8,8(1X,F9.5,5X))')
     &       ' total |',(POPG(IEFF(NOFF+I)),I=1,NPRI)
        NOFF = NOFF + NPRI
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck grossc */
      SUBROUTINE GROSSC(GROSS,NREFS,IFRP,NREF,TOTM,GVEC,NSTR,REFFIL)
C***********************************************************************
C
C     Determine gross contribution in projection analysis
C
C     Written by T.Saue Oct 5 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER REFFIL*6
      DIMENSION GROSS(2,0:2,NREFS)
      DIMENSION GVEC(*),NSTR(2,0:2,0:NREFS),REFFIL(*)
      GTOT = D0
      IVEC = 0
      DO I = 1,NREFS
C...    Negative-energy (positronic) orbitals
        GROSS(IFRP,1,I) = D0
        DO J = 1,NSTR(IFRP,2,I)
          IVEC = IVEC + 1
          GROSS(IFRP,1,I)   = GROSS(IFRP,1,I) + GVEC(IVEC)
        ENDDO
        GTOT = GTOT + GROSS(IFRP,1,I)
C...    Positive-energy (electronic) orbitals
        GROSS(IFRP,2,I) = D0
        DO J = NSTR(IFRP,2,I)+1,NSTR(IFRP,0,I)
          IVEC = IVEC + 1
          GROSS(IFRP,2,I)   = GROSS(IFRP,2,I) + GVEC(IVEC)
        ENDDO
        GROSS(IFRP,0,I) = GROSS(IFRP,1,I) + GROSS(IFRP,2,I)
        GTOT = GTOT + GROSS(IFRP,2,I)
        WRITE(LUPRI,'(3X,A6,8X,3X,F7.4,2(3X,A4,F7.4))')
     &      REFFIL(I),GROSS(IFRP,0,I),
     &      'E - ',GROSS(IFRP,2,I),'P - ',GROSS(IFRP,1,I)
      ENDDO
      GVEC(NREF+1) = TOTM - GTOT
      WRITE(LUPRI,'(3X,A,F10.4)')
     &    'Polarization: ',GVEC(NREF+1)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1rho  */
      SUBROUTINE GT1RHO(WORK,LWORK)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
      PARAMETER (D2=2.0D0)
      LOGICAL TOBE
      DIMENSION WORK(LWORK)
#include "dcbgen.h"
#include "dcbrho1.h"
#include "dcbbas.h"
C
C for NFMAT:
#include "dcbdhf.h"
C
      CALL QENTER('GT1RHO')
#include "memint.h"
C
      CALL TITLER('Density along bonds','*',122)
      IF (NFMAT .GT. 1) THEN
        WRITE(LUPRI,'(/A)')
     &  '--> Not implemented for open shell yet, nothing done.'
        GO TO 9999
      END IF
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
        WRITE(LUPRI,'(/A)')'GT1RHO: Coefficient file DFCOEF not found !'
        CALL QUIT('GT1RHO: Coefficients not found !')
      END IF
C     Memory allocation
CTROND: Consider open shell ! (NFMAT...)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAOS,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOOR,6       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Read coefficients + give info
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GT1RHO')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),TEMP,TEMP,
     &            TOTERG,3,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Generate density matrix
C     Be aware: Both density matrices are normed to 1
C     Coeffients may now be deallocated.....
C
      CALL GENDEN(WORK(KDMAT),WORK(KCMO),1,IPRHO1)
      CALL DSCAL(N2BBASXQ,D2,WORK(KDMAT),1)
CTROND      DO IOPEN = 1,NOPEN
CTROND        CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPRSCF)
CTROND      END DO
      CALL GT1RH1(WORK(KDMAT),WORK(KGAOS),WORK(KCOOR),WORK(KBUF),
     &            IPRHO1,WORK,KFREE,LFREE)
C
C     Memory deallocation
      CALL MEMREL('GT1RHO',WORK,KWORK,KWORK,KFREE,LFREE)
C
 9999 CALL QEXIT('GT1RHO')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1rh1 */
      SUBROUTINE GT1RH1(DMAT,GAO,COORD,BUF,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "nuclei.h"
      DIMENSION DMAT(*),GAO(*),COORD(3,2),BUF(*),WORK(*)
      CHARACTER FILNAM*12
#include "ibtfun.h"
C
      LURHO = 1
C
C     Loop over pairs of nuclear centers
C
      DO ICENTA = 1, NUCIND
         MULA  = ISTBNU(ICENTA)
         COORD(1,1) = CORD(1,ICENTA)
         COORD(2,1) = CORD(2,ICENTA)
         COORD(3,1) = CORD(3,ICENTA)
         DO ICENTB = ICENTA+1, NUCIND
           MULB = ISTBNU(ICENTB)
           IDEGB = 0 
           DO ISYMOP = 0, MAXOPR
             IF (IBTAND(ISYMOP,MULB) .EQ. 0) THEN
               IDEGB = IDEGB + 1
               COORD(1,2)=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENTB)
               COORD(2,2)=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENTB)
               COORD(3,2)=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENTB)
               WRITE(FILNAM,'(A2,2A4,I2)') 
     &           'RH',NAMN(ICENTA),NAMN(ICENTB),IDEGB
               DO I = 1,12
                 IF(FILNAM(I:I).EQ.' ') FILNAM(I:I)='_'
               ENDDO
               OPEN(LURHO,FILE =FILNAM,STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
               CALL GT1RH2(DMAT,GAO,COORD,BUF,LURHO,
     &             IPRINT,WORK,KFREE,LFREE)
               CLOSE(LURHO,STATUS='KEEP')
             ENDIF
           ENDDO
         ENDDO
         IF (MULT(MULA).GT.1) THEN
           IDEGB = 1
           DO ISYMOP = 1, MAXOPR
             IF (IBTAND(ISYMOP,MULA) .EQ. 0) THEN
               IDEGB = IDEGB + 1
               COORD(1,2)=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENTA)
               COORD(2,2)=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENTA)
               COORD(3,2)=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENTA)
               WRITE(FILNAM,'(A2,2(A3,I1)))') 
     &           'RH',NAMN(ICENTA),NAMN(ICENTA),IDEGB
               DO I = 1,10
                 IF(FILNAM(I:I).EQ.' ') FILNAM(I:I)='_'
               ENDDO
               OPEN(LURHO,FILE =FILNAM,STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
               CALL GT1RH2(DMAT,GAO,COORD,BUF,LURHO,
     &                     IPRINT,WORK,KFREE,LFREE)
               CLOSE(LURHO,STATUS='KEEP')
              ENDIF
           ENDDO
         ENDIF 
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Gt1rh2 */
      SUBROUTINE GT1RH2(DMAT,GAO,COORD,BUF,LURHO,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
#include "codata.h"
      DIMENSION DMAT(*),GAO(*),COORD(3,2),RHO(0:2),BUF(*),WORK(*)
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbrho1.h"
C
#include "sdpre.h"
      MC = 2
      IF(LEVYLE) MC = 1      
C
C     Calculate distance between centers
C
      DX = (COORD(1,2)-COORD(1,1))
      DY = (COORD(2,2)-COORD(2,1))
      DZ = (COORD(3,2)-COORD(3,1))
      DISTAN = SQRT(DX*DX+DY*DY+DZ*DZ)
      DX     = DX/DISTAN
      DY     = DY/DISTAN
      DZ     = DZ/DISTAN
      ASTEP  = DSTEP/XTANG
      NPOINT = INT(DISTAN/ASTEP)
      IF(MOD(NPOINT,2).EQ.1) NPOINT = NPOINT+1
      ASTEP  = DISTAN/SDPRE(NPOINT)
C
C     Calculate density along between centers:
C     If bond distance is 2R, then points are calculated between (-R,3R)
C 
      IP1 = -NPOINT/2
      IP2 = NPOINT - IP1 
      DO IP = IP1,IP2
        PP  = SDPRE(IP)*ASTEP
        PX  = COORD(1,1) + PP*DX
        PY  = COORD(2,1) + PP*DY
        PZ  = COORD(3,1) + PP*DZ
C       Get value of SOs at (PX,PY,PZ)
        CALL GETSOS(GAO,TEMP,TEMP,TEMP,TEMP,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),0,.FALSE.,IPRINT)
C       Get density in this point
        RHO(0) = D0
        DO IC = 1,MC
          RHO(IC) = D0
          CALL GETRH1(RHO(IC),IC,0,GAO,DMAT,BUF)
          RHO(0) = RHO(0) + RHO(IC)
        ENDDO
        PP = XTANG*PP
        WRITE(LURHO,'(F10.6,3E20.10)') PP,(RHO(I),I=0,MC)
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt3rho  */
      SUBROUTINE GT3RHO(WORK,LWORK)
C***********************************************************************
C
C     Prepare for 3-dimensional density plot
C
C     Written by Trond Saue Nov 5 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0)
      LOGICAL TOBE,LBIT
      DIMENSION WORK(LWORK),CUBE(3,0:4),LUCUBE(2)
#include "dcbgen.h"
#include "dcbrho.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
C
      CALL QENTER('GT3RHO')
c
c     Start timer
c
      TIMSTR = SECOND()
#include "memint.h"
C
      CALL TITLER('3-dimensional density plot','*',122)
C
      IF(LEVYLE) ILSCUB = 1
      IF(LBIT(ILSCUB,1)) THEN
         WRITE(LUPRI,'(A)') 
     &   ' GT3RHO : Gaussian cube file generated for large component.'
         LUCUBE(1) = 11
         OPEN(LUCUBE(1),FILE='rhoL.cube',STATUS='NEW',
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')   
      ENDIF
      IF (LBIT(ILSCUB,2)) THEN
         WRITE(LUPRI,'(A)') 
     &   ' GT3RHO : Gaussian cube file generated for small component.'
         LUCUBE(2) = 12
         OPEN(LUCUBE(2),FILE='rhoS.cube',STATUS='NEW',
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')   
      ENDIF
      CALL MEMGET('REAL',KATOM,4*NUCDEP,WORK,KFREE,LFREE)
C     Get cube parameters
      CALL CUBPAR(CUBE,WORK(KATOM))
C     Write header for cube file
      DO IC = 1,2
        IF(LBIT(ILSCUB,IC)) CALL CUBHED(LUCUBE(IC),CUBE,WORK(KATOM))
      ENDDO
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
C     Memory allocation
CTROND: Consider open shell ! (NFMAT...)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAOS ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF  ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Read coefficients + give restart info
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),TEMP,TEMP,
     &            TOTERG,3,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Generate density matrix
C     Be aware: Both density matrices are normed to 1
      CALL GENDEN(WORK(KDMAT),WORK(KCMO),1,IPRRHO)
CTROND      DO IOPEN = 1,NOPEN
CTROND        CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPRSCF)
CTROND      END DO
C     Coeffients may now be deallocated.....
      CALL MEMREL('GT3RHO.cmo',WORK,KCMO,KCMO,KFREE,LFREE)
C
      CALL MEMGET('REAL',KRHO,NCUBE(3)*2,WORK,KFREE,LFREE)
      CALL GT3RH1(LUCUBE,CUBE,WORK(KRHO),
     &            WORK(KDMAT),WORK(KGAOS),WORK(KBUF),
     &            WORK,KFREE,LFREE)
C
      DO IC = 1,2
        IF(LBIT(ILSCUB,IC)) CLOSE(LUCUBE(IC),STATUS='KEEP')
      ENDDO
C     Memory deallocation
      CALL MEMREL('GT3RHO',WORK,KWORK,KWORK,KFREE,LFREE)
C
c    
c     Stop timer
c           
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      CALL TIMTXT('>>> Time used in GT3RHO is',TIME,LUPRI)
      CALL QEXIT('GT3RHO')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GT3RHO: Coefficient file not found !'
      CALL QUIT('GT1RHO: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Gt3rh1 */
      SUBROUTINE GT3RH1(LUCUBE,CUBE,RHO,DMAT,GAO,BUF,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     3-dimensional plot of density
C
C     Written by Trond Saue Nov 5 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
      LOGICAL LBIT
      DIMENSION LUCUBE(*),CUBE(3,0:3),RHO(2,*),DMAT(*),GAO(*),
     &          BUF(*),WORK(*)
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbrho.h"
C
#include "sdpre.h"
C     
      DO I1 = 1,NCUBE(1)
        DO I2 = 1,NCUBE(2)
          DO I3 = 1,NCUBE(3)
            PX = CUBE(1,0)+(I1-1)*CUBE(1,1)
     &                    +(I2-1)*CUBE(1,2)
     &                    +(I3-1)*CUBE(1,3)
            PY = CUBE(2,0)+(I1-1)*CUBE(2,1)
     &                    +(I2-1)*CUBE(2,2)
     &                    +(I3-1)*CUBE(2,3)
            PZ = CUBE(3,0)+(I1-1)*CUBE(3,1)
     &                    +(I2-1)*CUBE(3,2)
     &                    +(I3-1)*CUBE(3,3)
C           Get value of SOs at (PX,PY,PZ)
            CALL GETSOS(GAO,TEMP,TEMP,TEMP,TEMP,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),0,.FALSE.,IPRRHO)
C           Get density in this point
            DO IC = 1,2
              IF(LBIT(ILSCUB,IC)) THEN
                RHO(IC,I3) = D0
                CALL GETRH1(RHO(IC,I3),IC,0,GAO,DMAT,BUF)
              ENDIF
            ENDDO
          ENDDO
C         Write buffer
          DO IC = 1,2
            IF(LBIT(ILSCUB,IC))
     &      WRITE(LUCUBE(IC),'(6E13.5)') (RHO(IC,I),I=1,NCUBE(3))
          ENDDO
        ENDDO
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Cubpar */
      SUBROUTINE CUBPAR(CUBE,ATOMS)
C***********************************************************************
C
C     Generate parameters for the cube. All coordinates in Angstroms.
C
C     NCUBE(0) - Number of atoms
C     NCUBE(1) - Number of points in the "X"-direction
C     NCUBE(2) - Number of points in the "Y"-direction
C     NCUBE(3) - Number of points in the "Z"-direction
C
C     CUBE(*,0) - Coordinates of initial point
C     CUBE(*,1) - Step vector for the "X"-direction
C     CUBE(*,2) - Step vector for the "Y"-direction
C     CUBE(*,3) - Step vector for the "Z"-direction
C
C     ATOMS(1,I) - x-coordinate of atom I
C     ATOMS(2,I) - y-coordinate of atom I
C     ATOMS(3,I) - z-coordinate of atom I
C     ATOMS(4,I) - Charge of atom I
C
C     Written by Trond Saue Nov 6 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "dcbrho.h"
      DIMENSION CUBE(3,0:3),ATOMS(4,NUCDEP)
C
#include "ibtfun.h"
C
      NCUBE(0) = NUCDEP
C     Initialize
      CALL DZERO(CUBE,12)
      DO I = 1,3
        TEMP = CORD(I,1)
        CUBE(I,0) = TEMP
        CUBE(I,I) = TEMP
      ENDDO
C     Loop over all centers
      JCENT = 0
      DO ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
C        Loop over all centers
         DO ISYMOP = 0, MAXOPR
         IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
           JCENT = JCENT + 1
           DO I = 1,3
             ATOMS(I,JCENT) = 
     &             PT(IBTAND(ISYMAX(I,1),ISYMOP))*CORD(I,ICENT)
             CUBE(I,0) = MIN(CUBE(I,0),ATOMS(I,JCENT))
             CUBE(I,I) = MAX(CUBE(I,I),ATOMS(I,JCENT))
           ENDDO
           ATOMS(4,JCENT)=CHARGE(ICENT)
         ENDIF
         ENDDO
      ENDDO      
C     Add some space, calculate step sizes and convert step vectors to Angstrom
      DO I = 1,3
        CUBE(I,0) = CUBE(I,0)  - CUBADJ(1)
        CUBE(I,I) = (CUBE(I,I) + CUBADJ(2))*XTANG
        CUBE(I,I) = (CUBE(I,I)-CUBE(I,0))/(NCUBE(I)-1)
      ENDDO
C      CALL DSCAL(12,XTANG,CUBE,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Cubhed */
      SUBROUTINE CUBHED(LUCUBE,CUBE,ATOMS)
C***********************************************************************
C
C     Write header for formatted Gaussian cube file
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "nuclei.h"
#include "dcbrho.h"
C
      DIMENSION CUBE(3,0:3),ATOMS(4,NUCDEP)
      WRITE(LUCUBE,10) TITLE
      WRITE(LUCUBE,10) 'Total density'
      DO I = 0,3
        WRITE(LUCUBE,20) NCUBE(I),(CUBE(J,I),J=1,3)
      ENDDO      
      DO I = 1,NUCDEP
        WRITE(LUCUBE,20) NINT(ATOMS(4,I)),ATOMS(4,I),
     &                   (ATOMS(J,I),J=1,3)
      ENDDO
C
      RETURN
   10 FORMAT(1X,A)
   20 FORMAT(I5,4F12.6)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
