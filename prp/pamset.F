C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : pamset.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxop */
      SUBROUTINE SETXOP(INDSTR,IOPSY,NXOPE,KXOPE,NXOPP,KXOPP,
     &                  WORK,KFREE,LFREE,SKIPEP,SKIPEE,QED,IPRINT)
C*****************************************************************************
C
C     Set orbital rotation index vectors
C
C     Orbitals are classified as:
C       inactive  : doubly occupied in all determinants
C       active    : singly occupied in some determinants
C       secondary : unoccupied in all determinants
C       positronic: positronic orbitals
C
C     Ordering is:
C       e-e rotations: (skipped if SKIPEE = T)
C          IA : inactive-active   
C          IS : inactive-secondary
C          AS : active-secondary
C       e-p rotations: (skipped if SKIPEP = T)
C          IP : inactive-positronic
C          AP : active-positronic
C         (SP : secondary-positronic   -  only for QED=T)
C
C          
C
C     Written by  Hans Joergen Aa. Jensen and T.Saue 27-Jan-1995
C     Last revision: July 21 1997 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C Used from common blocks:
C  DGROUP: NFSYM
C  DCBORB: IORB,NPSH,NFRO,NISH,NASH
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
C
      LOGICAL SKIPEP,SKIPEE,QED
      CHARACTER INDSTR(3,2)*72
      DIMENSION WORK(*),NSTR(0:2,3,2)
C
      CALL QENTER('SETXOP')
      KXOPE = KFREE
      KXOPP = KFREE
      NXOPE = 0
      NXOPP = 0
C
C     Allocate buffer areas
C     =====================
C
      NBUF = N2BBASX
      CALL MEMGET('INTE',KBUF ,NBUF          ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,MXFBAS*3*NFSYM,WORK,KFREE,LFREE)
C
C     Unpack strings
C     ==============
C
      CALL UNPSTR(INDSTR,3,WORK(KIBUF),MXFBAS,NSTR,
     &            WORK,KFREE,LFREE,IPRINT)
C
C     ***********************************
C     *** Electron-electron rotations ***
C     ***********************************
C      
 10   CONTINUE
      IF(SKIPEE) GOTO 20
      CALL SETXEE(WORK(KBUF),IOPSY,NXOPE,WORK(KIBUF),NSTR,IPRINT)
      IF(NXOPE.GT.0) THEN
        CALL MEMREL('SETXOP.e-e',WORK,1,KBUF,KFREE,LFREE)
        CALL MEMGET('INTE',KXOPE,2*NXOPE,WORK,KFREE,LFREE)
        NBUF = N2BBASX - 2*(NXOPE+1)
        CALL MEMGET('INTE',KBUF ,NBUF,WORK,KFREE,LFREE)
      ENDIF
C
C     ***********************************
C     *** Electron-positron rotations ***
C     ***********************************
C      
 20   CONTINUE
      IF(SKIPEP) GOTO 30
      CALL SETXEP(WORK(KBUF),IOPSY,NXOPP,QED,WORK(KIBUF),NSTR,IPRINT)
      CALL MEMREL('SETXOP.e-p',WORK,1,KBUF,KFREE,LFREE)
      CALL MEMGET('INTE',KXOPP,2*NXOPP,WORK,KFREE,LFREE)
C
 30   CONTINUE
      CALL QEXIT('SETXOP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setxoe */
      SUBROUTINE SETXEE(JXOPE,IOPSY,NXOPE,IND,NIND,IPRINT)
C***********************************************************************
C   
C     Generate electron-electron rotations
C
C     Written by T.Saue July 18 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
C
      DIMENSION JXOPE(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
C
      CALL QENTER('SETXEE')

      NXOPE = 0
      DO 10 ISYM = 1,NFSYM
        JSYM  = MOD(ISYM+IOPSY,2) + 1
C
        ISTI  = NIND(1,1,ISYM) + 1
        IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
        ISTA  = NIND(1,2,ISYM) + 1
        IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
C
        JSTA  = NIND(1,2,JSYM) + 1
        JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
        JSTS  = NIND(1,3,JSYM) + 1
        JENDS = NIND(1,3,JSYM) + NIND(2,3,JSYM)
C
        IF(NIND(1,1,ISYM).NE.0) THEN
          WRITE(LUPRI,'(A,I5)')
     &      '*** WARNING(SETXEE): Inactive positrons detected:',
     &      NIND(1,1,ISYM)
        ENDIF
        IF(NIND(2,1,ISYM).EQ.0) GOTO 300
C
C       Inactive - active elements
C       ==========================
C
 100    CONTINUE
        IF(NIND(1,2,JSYM).NE.0) THEN
          WRITE(LUPRI,'(A,I5)')
     &      '*** WARNING(SETXEE): Active positrons detected:',
     &      NIND(1,2,JSYM)
        ENDIF
        IF (NIND(2,1,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 200
        DO IN = ISTI,IENDI
          DO JA = JSTA,JENDA
            NXOPE          = NXOPE + 1
            JXOPE(1,NXOPE) = IND(IN,1,ISYM)
            JXOPE(2,NXOPE) = IND(JA,2,JSYM)
          ENDDO
        ENDDO
C
C       Inactive-secondary elements
C       ===========================
C
 200    CONTINUE
        IF (NIND(2,1,ISYM).EQ.0.OR.NIND(2,3,JSYM).EQ.0) GO TO 300
        DO IN = ISTI,IENDI
          DO JS = JSTS,JENDS
            NXOPE          = NXOPE + 1
            JXOPE(1,NXOPE) = IND(IN,1,ISYM)
            JXOPE(2,NXOPE) = IND(JS,3,JSYM)
          ENDDO
        ENDDO
C
C       Active-secondary elements
C       =========================
C
 300    CONTINUE
        IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,3,JSYM).EQ.0) GO TO 10
        DO IA = ISTA,IENDA
          DO JS = JSTS,JENDS
            NXOPE          = NXOPE+1
            JXOPE(1,NXOPE) = IND(IN,2,ISYM)
            JXOPE(2,NXOPE) = IND(JS,3,JSYM)
          ENDDO
        ENDDO
 10   CONTINUE
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from SETXEE',-1)        
        WRITE(LUPRI,'(A,I10)') 'e-e rotations: ',NXOPE
        WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,NXOPE)
      ENDIF
C
      CALL QEXIT('SETXEE')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setxep */
      SUBROUTINE SETXEP(JXOPP,IOPSY,NXOPP,QED,IND,NIND,IPRINT)
C***********************************************************************
C   
C     Generate electron-positron rotations
C
C     Written by T.Saue July 21 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
C
      LOGICAL QED
      DIMENSION JXOPP(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
C
      CALL QENTER('SETXEP')
C
      NXOPP = 0
      DO 10 ISYM = 1,NFSYM
        JSYM  = MOD(ISYM+IOPSY,2) + 1
C
        ISTI  = NIND(1,1,ISYM) + 1
        IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
        ISTA  = NIND(1,2,ISYM) + 1
        IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
C
C       Inactive-positronic elements
C       ===========================
C
 100    CONTINUE 
        IF (NIND(2,1,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 200
        DO IN = ISTI,IENDI
          DO JP = 1,NIND(1,3,JSYM)
            NXOPP          = NXOPP+1
            JXOPP(1,NXOPP) = IND(IN,1,ISYM)
            JXOPP(2,NXOPP) = IND(JP,3,JSYM)
          ENDDO
        ENDDO
C
C       Active-positronic elements
C       ==========================
C
 200    CONTINUE
        IF (NIND(2,2,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 300
        DO IA = ISTA,IENDA
          DO JP = 1,NIND(1,3,JSYM)
            NXOPP          = NXOPP+1
            JXOPP(1,NXOPP) = IND(IA,1,ISYM)
            JXOPP(2,NXOPP) = IND(JP,3,JSYM)
          ENDDO
        ENDDO
C
C       Secondary-positronic elements
C       =============================
C
 300    CONTINUE
        IF (QED) THEN
          IF (NIND(1,3,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 10
          ISTS  = NIND(1,3,ISYM) + 1
          IENDS = NIND(1,3,ISYM) + NIND(2,3,ISYM)
          DO IS = ISTS,IENDS
            DO JP = 1,NIND(1,3,JSYM)
              NXOPP          = NXOPP+1
              JXOPP(1,NXOPP) = IND(IS,1,ISYM)
              JXOPP(2,NXOPP) = IND(JP,3,JSYM)
            ENDDO
          ENDDO
        ENDIF
C
   10 CONTINUE
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from SETXEP',-1)        
        WRITE(LUPRI,'(A,I10)') 'e-p rotations: ',NXOPP
        WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,NXOPP)
      ENDIF
C
      CALL QEXIT('SETXEP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck getfck */
      SUBROUTINE GETFCK(FMO,IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Get Fock-matrix in MO-basis; use coefficients for MO-transformation
C
C     Written by T.Saue, Sep 16 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
C
      LOGICAL TOBE
      DIMENSION FMO(*),WORK(*)
C
      CALL QENTER('GETFCK')
      KFRSAV = KFREE
C
C     First check if Fock-matrix is on file PAMFCK
C
      LUBUF = 22
      INQUIRE(FILE='PAMFCK',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'PAMFCK','OLD','GETFC1')
        CALL READT(LUBUF,N2ORBXQ,FMO)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        CALL MEMGET('REAL',KFAO,N2BBASXQ,WORK,KFREE,LFREE)      
        CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C
        CALL GETFC1(FMO,WORK(KFAO),WORK(KCMO),
     &            IPRINT,WORK,KFREE,LFREE)
C
        CALL OPNFIL(LUBUF,'PAMFCK','NEW','GETFC1')
        CALL WRITT(LUBUF,N2ORBXQ,FMO)
        CLOSE(LUBUF,STATUS='KEEP')
        CALL MEMREL('GETFCK',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
          CALL HEADER('GETFCK:Total Fock matrix in MO-basis',-1)
          CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      CALL QEXIT('GETFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETFC1 */
      SUBROUTINE GETFC1(FMO,FAO,CMO,IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Get Fock-matrix in MO-basis
C
C     Written by T.Saue, Sep 16 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include <dcbgen.h>
#include <dcbfir.h>
#include <dcbdhf.h>
#include <dcbham.h>
#include <blocks.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <cbihr2.h>
CTOBS: Note the use of dhfcal ! Should be generalized !!!
      LOGICAL TOBE,BF1,BF2,BF3
      DIMENSION FMO(*),FAO(*),CMO(*),WORK(*),IQBF(4)
C
      CALL QENTER('GETFC1')
      LUBUF = 22
C
C     1. Get coefficients
C     ===================
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','GETFC1')
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUM,IDUM,TOTERG,4)
      CLOSE(LUBUF,STATUS='KEEP')
C
C     1. Get two-electron Fock matrix in AO-basis
C     ===========================================
C
      INQUIRE(FILE='DFFCK2',EXIST=TOBE)   
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK2','OLD','GETFC1')
        CALL REAFCK(LUBUF,FMO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GETFCK: No 2-Fock found. Regenerating.'
        NPOS = 0
        IF(PARCAL) THEN
          BF1  = FIRST1
          BF2  = FIRST2
          BF3  = FIRST3
          FIRST1 = .TRUE.
          FIRST2 = .TRUE.
          FIRST3 = .TRUE.
          NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
          NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
          NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
        ENDIF
        CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
C       Closed shell Hartree - Fock
C       ===========================
CTROND: This must be generalized !!!
        NFMAT = 1
C       Totally symmetric operator
        ISYMOP(1) = 1
C       Fock matrix type
        IFCKOP(1) = 1
C       Hermitian operator
        IHRMOP(1) = 1
        INQUIRE(FILE='DFDENS',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFDENS','OLD','GETFC1')
          CALL READNS(LUBUF,FAO)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
          CALL DENSTY(FAO,CMO,IPRINT)        
        ENDIF
        INTFLG = 7
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMO,FAO,NFMAT,
     &       WORK(KPOS),INTFLG,IPRINT,WORK(KFREE),LFREE)
        IF(PARCAL) THEN
          FIRST1 = BF1
          FIRST2 = BF2
          FIRST3 = BF3
        ENDIF
      ENDIF
C
C     3a. Solvent: get solvent contribution
C     =====================================
C
      IF(SOLVEN) THEN
C        CALL DENSTY(FAO,CMO,IPRINT)        
        NF = 1
        CALL SOLFCK(FMO,FAO,NF,ESOLVE,ESOLVN,WORK(KFREE),LFREE,IPRSOL)
      ENDIF           
C
C     3. Get one-electron Fock matrix
C     ===============================
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK1','OLD','GETFC1')
        CALL REAFCK(LUBUF,FAO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GETFCK: No 1-Fock found. Regenerating.'
        CALL ONEFCK(FAO,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
C     4. Add 1- and 2-electron matrix and transform to MO-basis
C     =========================================================
C
      CALL DAXPY(N2BBASXQ,D1,FMO,1,FAO,1) 
      CALL DZERO(FMO,N2ORBXQ)
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &              FAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),
     &              FMO(I2ORBX(I,I)+1),NORBT,NORBT,NBF,IQBF,
     &              CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,IPRINT)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** GETFC1: NBF.NE.NZ !'
 10   CONTINUE
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GETFC1:Total Fock matrix in MO-basis',-1)
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &    IPQTOQ(1,0),LUPRI)
      ENDIF
C
      CALL QEXIT('GETFC1')
      RETURN
C
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GETFC1: Coefficient file not found !'
      CALL QUIT('GETFC1: Coefficients not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck unpstr */
      SUBROUTINE UNPSTR(INDSTR,NTYP,IND,LIND,NIND,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     General routine for unpacking of orbital strings
C
C     Written by T.Saue July 19 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>      
      CHARACTER INDSTR(NTYP,NFSYM)*72
      DIMENSION IND(LIND,NTYP,NFSYM),NIND(2,NTYP,NFSYM),
     &          WORK(*)
C
      DO IFRP = 1,NFSYM
        DO ITYP = 1,NTYP
          NVEC = 1
          CALL  NUMLST(INDSTR(ITYP,IFRP),IND(1,ITYP,IFRP),
     &                 NFBAS(IFRP,0),-NFBAS(IFRP,2),NFBAS(IFRP,1),
     &                 NVEC,WORK(KFREE),LFREE)
          CALL ORBCNT(IND(1,ITYP,IFRP),NVEC,IFRP,
     &                NIND(1,ITYP,IFRP),NIND(2,ITYP,IFRP))
C
C         Absolute indices of positronic orbitals
C
          IOFF = NPSH(IFRP) + 1 + IORB(IFRP)
          DO I = 1,NIND(1,ITYP,IFRP)
            IND(I,ITYP,IFRP) = IOFF + IND(I,ITYP,IFRP)
          ENDDO
C
C         Absolute indices of electronic orbitals
C
          IOFF = NPSH(IFRP) + IORB(IFRP)
          JOFF = NIND(1,ITYP,IFRP)
          DO I = 1,NIND(2,ITYP,IFRP)
            IND(I+JOFF,ITYP,IFRP) = IOFF + IND(I+JOFF,ITYP,IFRP) 
          ENDDO
        ENDDO
      ENDDO
C
C     Print section
C     =============      
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Output from UNPSTR',-1)
        DO ITYP = 1,NTYP
          WRITE(LUPRI,'(A,I5)') '* Orbital string no.',ITYP
          DO IFRP = 1,NFSYM
            WRITE(LUPRI,'(3X,A,A3)') '* Fermion ircop ',FREP(IFRP)
            IF(NIND(1,ITYP,IFRP).EQ.0) THEN
              WRITE(LUPRI,'(5X,A)') 'No positronic orbitals.'
            ELSE
              WRITE(LUPRI,'(5X,A)') 'Positronic orbitals:'
              WRITE(LUPRI,'(8X,12I5)')
     &             (IND(I,ITYP,IFRP),I=1,NIND(1,ITYP,IFRP))
            ENDIF
            IF(NIND(2,ITYP,IFRP).EQ.0) THEN
              WRITE(LUPRI,'(5X,A)') 'No electronic orbitals.'
            ELSE
              JOFF = NIND(1,ITYP,IFRP)
              WRITE(LUPRI,'(5X,A)') 'Electronic orbitals:'
              WRITE(LUPRI,'(8X,12I5)')
     &             (IND(I+JOFF,ITYP,IFRP),I=1,NIND(2,ITYP,IFRP))
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck REDPAR */
      SUBROUTINE REDPAR(GPCI,GPOE,GPOP,JCONF,JXOPE,JXOPP,
     &                  THRCOM,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Reduce number of parameters by deleting parameters
C     for which the corresponding gradient element is
C     below threshold THRCOM
C
C     Written by T.Saue - July 22 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,DC = 100.0D0)
C
#include <dcbxrs.h>
#include <dgroup.h>
      DIMENSION GPCI(*),GPOE(*),GPOP(*),WORK(*),
     &          JCONF(*),JXOPE(2,NZXOPE),JXOPP(2,NZXOPP)
C
      KFRSAV = KFREE
      IF(IPRINT.GE.4) CALL HEADER('Output from REDPAR',-1)        

C
C     Configurational part
C     ====================
C
      IF(NZCONF.GT.0) THEN
        WRITE(LUPRI,'(A)') '* REDPAR: No configurational code.'
        CALL QUIT('REDPAR: No conf.code')
      ENDIF
C
C     e-e orbital rotations
C     =====================
C
      IF(NZXOPE.GT.0) THEN
        CALL MEMGET('REAL',KBUF,NZXOPEQ,WORK,KFREE,LFREE)
        II = 0
        DO I = 1,NZXOPE
          DTST = D0
          IZOFF = 0
          DO IZ = 1,NZ
            DTST = DTST + GPOE(I+IZOFF)*GPOE(I+IZOFF)  
            IZOFF = IZOFF + NZXOPE
          ENDDO
          IF(SQRT(DTST).GE.THRCOM) THEN
            II = II + 1
            JXOPE(1,II) = JXOPE(1,I)
            JXOPE(2,II) = JXOPE(2,I)
            IOFF = KBUF-1
            JOFF = 0
            DO IZ = 1,NZ
              WORK(IOFF+II) = GPOE(JOFF+I)
              IOFF = IOFF + NZXOPE
              JOFF = JOFF + NZXOPE
            ENDDO
          ENDIF
        ENDDO
        IF(II.LT.NZXOPE) THEN
          IOFF = KBUF
          JOFF = 1
          DO IZ = 1,NZ
            CALL DCOPY(II,WORK(IOFF),1,GPOE(JOFF),1)
            IOFF = IOFF + NZXOPE
            JOFF = JOFF + II
          ENDDO
          IF(IPRINT.GE.4) THEN
            DTST = DC*(NZXOPE-II)/NZXOPE
            WRITE(LUPRI,'(A,I10,A,F5.1,A)') 
     &            'e-e rotations: ',II,' (reduction: ',DTST,'%)'
            WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,II)
            CALL PRQMAT(GPOE,II,1,II,1,NZ,1,LUPRI)
          ENDIF
          NZXOPE  = II
          NZXOPEQ = NZXOPE*NZ
        ENDIF
        CALL MEMREL('REDPAR.ee',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF      
C
C     e-p orbital rotations
C     =====================
C
      IF(NZXOPP.GT.0) THEN
        CALL MEMGET('REAL',KBUF,NZXOPPQ,WORK,KFREE,LFREE)
        II = 0
        DO I = 1,NZXOPP
          DTST = D0
          IZOFF = 0
          DO IZ = 1,NZ
            DTST = DTST + GPOP(I+IZOFF)*GPOP(I+IZOFF)  
            IZOFF = IZOFF + NZXOPP
          ENDDO
          IF(SQRT(DTST).GE.THRCOM) THEN
            II = II + 1
            JXOPP(1,II) = JXOPP(1,I)
            JXOPP(2,II) = JXOPP(2,I)
            IOFF = KBUF-1
            JOFF = 0
            DO IZ = 1,NZ
              WORK(IOFF+II) = GPOP(JOFF+I)
              IOFF = IOFF + NZXOPP
              JOFF = JOFF + NZXOPP
            ENDDO
          ENDIF
        ENDDO
        IF(II.LT.NZXOPP) THEN
          IOFF = KBUF
          JOFF = 1
          DO IZ = 1,NZ
            CALL DCOPY(II,WORK(IOFF),1,GPOP(JOFF),1)
            IOFF = IOFF + NZXOPP
            JOFF = JOFF + II
          ENDDO
          IF(IPRINT.GE.4) THEN
            DTST = DC*(NZXOPP-II)/NZXOPP
            WRITE(LUPRI,'(A,I10,A,F5.1,A)') 
     &            'e-p rotations: ',II,'(reduction: ',DTST,'%)'
            WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,II)
          ENDIF
          NZXOPP  = II
          NZXOPPQ = NZXOPP*NZ
        ENDIF
        CALL MEMREL('REDPAR.ep',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF      
C
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT *NZ
C
C     Total number of parameters
C     ==========================
C
      NZVAR    = NZCONF + NZXOPT
      NZVARQ   = NZVAR*NZ
      RETURN
      END
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DEFORB */
      SUBROUTINE DEFORB(IORBCL,NIND,NTYP,INDSTR,WORK,KFREE,LFREE,
     &                  IPRINT)
C***********************************************************************
C
C     This routine reads orbital strings and classifies
C     and indexes orbitals
C
C     Written by T.Saue Oct 3 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
      CHARACTER INDSTR(NTYP,2)*72
      DIMENSION IORBCL(2,NORBT),NIND(2,NTYP,2),WORK(*)
C
      CALL QENTER('DEFORB')
      KFRSAV = KFREE
C
C     Allocate buffer areas
C     =====================
C
      CALL MEMGET('INTE',KIND ,MXFBAS*NTYP*NFSYM,WORK,KFREE,LFREE)
      CALL DEFOR1(IORBCL,NIND,NTYP,INDSTR,WORK(KIND),
     &            WORK,KFREE,LFREE,IPRINT)
      CALL MEMREL('DEFORB',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DEFORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DEFOR1 */
      SUBROUTINE DEFOR1(IORBCL,NIND,NTYP,INDSTR,IND,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     This routine reads orbital strings and classifies
C     and indexes orbitals
C
C     Written by T.Saue Oct 3 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
      CHARACTER INDSTR(NTYP,2)*72
      DIMENSION IORBCL(2,NORBT),IND(MXFBAS,NTYP,NFSYM),
     &          NIND(2,NTYP,2),WORK(KFREE)
C
      NDIM = NORBT*2
      CALL IZERO(IORBCL,NDIM)
C
C     Unpack strings
C     ==============
C
      DO ITYP = 1,NTYP
        KOFF = 0
        DO IFRP = 1,NFSYM
          NVEC = 1
          CALL  NUMLST(INDSTR(ITYP,IFRP),IND(1,ITYP,IFRP),
     &                 NFBAS(IFRP,0),-NFBAS(IFRP,2),NFBAS(IFRP,1),
     &                 NVEC,WORK(KFREE),LFREE)
          CALL ORBCNT(IND(1,ITYP,IFRP),NVEC,IFRP,
     &                NIND(1,ITYP,IFRP),NIND(2,ITYP,IFRP))
C
C         Absolute indices of positronic orbitals
C
          IOFF = NPSH(IFRP) + 1 + IORB(IFRP)
          DO I = 1,NIND(1,ITYP,IFRP)
            J           = IOFF + IND(I,ITYP,IFRP)
            IORBCL(1,J) = ITYP
            IORBCL(2,J) = KOFF + I
          ENDDO
          KOFF = KOFF + NIND(1,ITYP,IFRP)
C
C         Absolute indices of electronic orbitals
C
          ICLS =  ITYP
          IOFF = NPSH(IFRP) + IORB(IFRP)
          JOFF = NIND(1,ITYP,IFRP)
          DO I = 1,NIND(2,ITYP,IFRP)
            J = IOFF + IND(I+JOFF,ITYP,IFRP) 
            IORBCL(1,J) = ITYP
            IORBCL(2,J) = KOFF + I
          ENDDO
          KOFF = KOFF + NIND(2,ITYP,IFRP)
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from DEFORB',-1)
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A3)') '* Fermion ircop:',FREP(IFRP)
          WRITE(LUPRI,'(5X,3I5)')
     &    (I,IORBCL(1,IORB(IFRP)+I),IORBCL(2,IORB(IFRP)+I),
     &     I = 1,NORB(IFRP))
        ENDDO
      ENDIF
      RETURN
      END
