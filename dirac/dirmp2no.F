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
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirmp2no.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mp2noinp */
      SUBROUTINE MP2NOINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for generation of MP2 natural orbitals
C
C     Written by S.Knecht and H.J.Aa. Jensen - March 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmp2no.h"
C
      LOGICAL SET, NEWDEF, RESET,LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
C  FIXME ideas for options: (7-Mar-2008)
C  .MAXVIRT - maximum number of virtuals in MP2
C             if energy criterium gives more then max energy will be
C             reduced
C  .FOCKNO  - transform to selected NO's to canonical Fock orbitals in
C             subspace (?? for better convergence in CI?)
C  .SCHEME  - for MOLTRA (default 4 or 6 ?)
C  .SELECT  - select min NO for active orbitals afterwards (relccsd,
C             lucita, kr-mcscf, luciarel, ?).
C             Default: e.g. "NO-occ 1.99 0.001 10%" in analogy with the
C             "energy -10.0 20.0 0.2"
C
C  
      SAVE SET
      DATA TABLE /'.ACTIVE','.PRINT ','.INTFLG','.F2CONT',
     &            'xxxxxxx','xxxxxxx','xxxxxxx','xXXXXXX'/
      DATA SET/.FALSE./
C
#include "ibtfun.h"
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C
C     Initialize /CBIMP2NO/
C     =====================
C
      IPRMP2NO = 0
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      IFCMP2NO = 1
      IFEMP2NO = 1
      DO I = 1,NFSYM
        VECMP2NO(I) = 'energy -10.0 20.0 0.2 '
        IF(NOCC(I).EQ.0) THEN
          WRITE(VECMP2NO(I),'(I1)') 1
        ELSEIF(NOCC(I).GT.0) THEN
          WRITE(VECMP2NO(I),'(I1,A2,I4)') 1,'..',NOCC(I)
        ENDIF
      ENDDO
C
C     Process input from CBIMP2NO
C     =========================
C
      NEWDEF = (WORD .EQ. '*MP2 NO')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in *MP2 NO.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in *MP2 NO.')
    1          CONTINUE
C&&&& ACTIVE: Spinors to be included in the MP2 calculation
                  READ(LUCMD,'(A72)') (VECMP2NO(I),I=1,NFSYM)
               GO TO 100
    2          CONTINUE
C&&&& PRINT: General print level
                  READ(LUCMD,*) IPRMP2NO
               GO TO 100
    3          CONTINUE
C&&& INTFLG: Specify what two-itegrals should be included in this run
                  IF (IGTINT.EQ.1) THEN
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT,IGTINT
                  ELSE
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  END IF
               GO TO 100
    4          CONTINUE
C&&& F2CONT: Specify Coulomb/exchange contributions
                  READ(LUCMD,*) IFCMP2NO,IFEMP2NO
               GO TO 100
    5          CONTINUE
C&&& 
               GO TO 100
    6          CONTINUE
C&&& 
               GO TO 100
    7          CONTINUE
C&&&
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in *MP2 NO.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *MP2 NO.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      INTMP2NO = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
      IF2MP2NO = IFCMP2NO + 2*IFEMP2NO
C
C     Print section
C     =============
C
      IF(.NOT.DOMP2NO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' MP2NOGEN: MP2 natural orbitals'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I3)') ' * General print level:',IPRMP2NO
      WRITE(LUPRI,'(A)') ' * Vectors defining the DHF potential:'
      NTOT = 0
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECMP2NO(I),IDUMMY,NFBAS(I,0),
     &            -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)')
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          NTOT = NTOT + NVEC
          WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &    '- Orbitals in fermion ircop ',FREP(I),' :',VECMP2NO(I)
        ENDIF
      ENDDO
      IF(ILLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No LL contributions'
      IF(ISLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SL contributions'
      IF(ISSINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SS contributions'
      IF(IGTINT.EQ.0) WRITE(LUPRI,'(A)') ' * No GT contributions'
      IF(IF2MP2NO.NE.3) THEN
        WRITE(LUPRI,'(A)')' * Two-electron contributions restricted to:'
        IF(LBIT(IF2MP2NO,1)) WRITE(LUPRI,'(3X,A)')
     &     '+ direct contribution'
        IF(LBIT(IF2MP2NO,2)) WRITE(LUPRI,'(3X,A)')
     &     '+ exchange contribution'
      ENDIF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mp2nogen */
      SUBROUTINE MP2NOGEN(WORK,LWORK)
C***********************************************************************
C
C     Generate MP2 natural orbitals
C
C     Written by S.Knecht and H.J.Aa. Jensen - March 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "aovec.h"
C
#include "dcbbas.h"
#include "dcborb.h"
C
#include "dcbfir.h"
#include "dcbgen.h"
#include "dcbmp2no.h"
#include "blocks.h"
C
      LOGICAL BF1,BF2,BF3
      DIMENSION WORK(*)
C
#include "memint.h"
      CALL TITLER('MP2 natural orbitals','*',125)
C
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBRP,NTBAS(0),      WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFOCK,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,N2BBASXQ,WORK,KFREE,LFREE)
      NPOS = 0
      IF (PARCAL) THEN
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
      CALL MP2NOGE1(WORK(KCMO),WORK(KEIG),WORK(KIBRP),
     &            WORK(KFOCK),WORK(KDMAT),
     &            WORK(KPOS),WORK(KFREE),LFREE)
      CALL MEMREL('MP2NOGEN',WORK,KWORK,KWORK,KFREE,LFREE)
C
      IF (PARCAL) THEN
         FIRST1 = BF1
         FIRST2 = BF2
         FIRST3 = BF3
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mp2noge1 */
      SUBROUTINE MP2NOGE1(CMO,EIG,IBEIG,FOCK,DMAT,NPOS,WORK,LWORK)
C***********************************************************************
C     Generate MP2 natural orbitals
C
C     Written by S.Knecht and H.J.Aa. Jensen - March 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DM1=-1.0D0,DM2=-2.0D0,D2=2.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbmp2no.h"
#include "dcbham.h"
      LOGICAL TOBE
      DIMENSION CMO(*),EIG(*),IBEIG(*),FOCK(*),DMAT(*),NPOS(*),WORK(*)
      DIMENSION IC(2),IE(2),ID(2)
      CHARACTER CTMP*1
      CALL QENTER('MP2NOGEN')
#include "memint.h"
C
C     Initialization
C
      LUBUF = 22
C     full Fock matrix
      IFCMP2NO = 1
C     only direct (Coulomb) contributions 
      IF(IF2MP2NO.EQ.1) IFCMP2NO = 2
C     only exchange contributions
      IF(IF2MP2NO.EQ.2) IFCMP2NO = 3
      ISYMP2NO = 1
      IHRMP2NO = 1
C
C     Get all coefficients
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('MP2NOGEN: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','MP2NOGEN')
         CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14,
     &               WORK,KFREE,LFREE)
      ENDIF
#ifdef mp2no_fock
      IF(IONMVO) THEN
C
C       Generate MP2NOs by the Bauschlicher method (ionic Fock)
C
        CALL DENMP2NO(DMAT,CMO,FOCK,WORK(KFREE),LFREE)
        CALL TWOFCK(ISYMVO,IHRMVO,IFCMVO,FOCK,DMAT,1,NPOS,
     &              INTMVO,IPRMVO,WORK(KFREE),LFREE)
        INQUIRE(FILE='DFFCK1',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFFCK1','OLD','MP2NOGE1')
          CALL REAFCK(LUBUF,DMAT,.TRUE.,1)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
          CALL ONEFCK(DMAT,IPRINT,WORK(KFREE),LFREE)
        ENDIF
        CALL DAXPY(N2BBASXQ,D1,DMAT,1,FOCK,1)
      ENDIF
#endif
C
C       First transform to virtual MO - basis;
C       this liberates FOCK. Transformed matrix in DMAT
C
      ND = 1
      DO I = 1,NFSYM
      IF (NSSH(I).GT.0 ) THEN
        IC(I) = ICMOQ(I)+NFBAS(I,0)*(NOCC(I)+NPSH(I))+1
        IE(I) = IORB(I) + NOCC(I) + NPSH(I) + 1
        ID(I) = ND
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &           NSSH(I),NSSH(I),
     &           FOCK(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DMAT(ID(I)),NSSH(I),NSSH(I),NZ,IPQTOQ(1,0),
     &           CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRMP2NO)
        IF (IPRMP2NO.GE.10 ) THEN
          WRITE(CTMP,'(I1)') I
          CALL HEADER('MP2NO: Fock matrix in virtual MO-basis'//
     &                 ', corep '//CTMP,-1)
          CALL PRQMAT(DMAT(ID(I)),NSSH(I),NSSH(I),
     &                NSSH(I),NSSH(I),NZ,IPQTOQ(1,0),LUPRI)
        END IF
        ND = ND + NSSH(I)*NFBAS(I,0)*NZ
      ENDIF
      ENDDO
C
C     Diagonalize matrix and backtransform vectors
C
      DO I = 1,NFSYM
      IF (NSSH(I).GT.0 ) THEN
        CALL QDIAG(NZ,NSSH(I),DMAT(ID(I)),NSSH(I),NSSH(I),
     &             EIG(IE(I)),1,FOCK(ID(I)),NSSH(I),
     &             NSSH(I),WORK(KFREE),LFREE,IERR)
        CALL BCKTR1(DMAT(ID(I)),NFBAS(I,0),NSSH(I),
     &              FOCK(ID(I)),NSSH(I),NSSH(I),
     &              NSSH(I),NZ,NSSH(I),1,NFBAS(I,0),
     &              CMO(IC(I)),NFBAS(I,0),NORB(I),NZ,
     &              IPRMP2NO)
        IF(IPRMP2NO.GE.2) THEN
          WRITE(LUPRI,'(A,A3)')
     &     '* MP2 natural orbital occupancies of ircop ',FREP(I)
          WRITE(LUPRI,'(5F15.8)') (EIG(IE(I)+J),J=0,(NSSH(I)-1))
        ENDIF
C
C       Insert new coefficients
C
        NDIM = NSSH(I)*NFBAS(I,0)
        ICC = IC(I)
        IDD = ID(I)
        DO IZ = 1,NZ
          CALL DCOPY(NDIM,DMAT(IDD),1,CMO(ICC),1)
          ICC = ICC + NCMO(I)
          IDD = IDD + NFBAS(I,0)*NSSH(I)
        ENDDO
C
C       Insert dummy eigenvalues...
C       FIXME: natural orbital occupations ? or if MP2-FOCK, then
C       orbital eigenvalues /hjaaj 7-mar-08
C
Chj     CALL DCOPY(NSSH(I),DUMMY,0,EIG(IE(I)),1)
      END IF
      END DO
      CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL QEXIT('MP2NOGEN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck denmp2no */
      SUBROUTINE DENMP2NO(DMAT,CMO,BUF,WORK,LWORK)
C***********************************************************************
C
C     Generate density matrix for chosen MP2NO vectors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmp2no.h"
      DIMENSION DMAT(*),CMO(*),BUF(*),WORK(*)
      DIMENSION NR(2),IR(2,2),KVEC(2,2),NVEC(2,0:2)
      DIMENSION ICMV(2)
#include "memint.h"
      KFRSAV = KFREE
C
C     Set up index arrays for coefficients contributing
C     =================================================
C
      IC   = 1
      DO I = 1,NFSYM
        NR(I)   = NESH(I)
        IR(1,I) = 0
        IR(2,I) = NESH(I)
      ENDDO
      DO IFRP = 1,NFSYM
        CALL MEMGET('INTE',KVEC(IFRP,1),NR(IFRP),WORK,KFREE,LFREE)
        NVEC(IFRP,0) = 1
        CALL  NUMLST(VECMP2NO(IFRP),WORK(KVEC(IFRP,1)),
     &               NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &               IFRP,NVEC(IFRP,0),WORK(KFREE),LFREE)
        CALL ORBCNT(WORK(KVEC(IFRP,1)),NVEC(IFRP,0),IFRP,
     &              NVEC(IFRP,2),NVEC(IFRP,1))
        NVEC(IFRP,0) = NVEC(IFRP,1) + NVEC(IFRP,2)
        CALL MEMREL('DENMP2NO',WORK,KFRSAV,KVEC(IFRP,1),KFREE,LFREE)
        CALL MEMGET('INTE',KVEC(IFRP,1),NVEC(IFRP,0),
     &              WORK,KFREE,LFREE)
      ENDDO
C
C     Calculate dimensions of the coefficient array
C
      NCMV = 0
      DO IFRP = 1,NFSYM
        ICMV(IFRP) = NCMV
        NCMV       = NCMV + NFBAS(IFRP,0)*NVEC(IFRP,0)*NZ
      ENDDO
C
C     Select the set that we need
C
      DO IFRP = 1, NFSYM
         CALL SELCFS (CMO(ICMOQ(IFRP)+1),IFRP,BUF(ICMV(IFRP)+1),
     &                NVEC(IFRP,0),
     &                WORK(KVEC(IFRP,1)),NVEC(IFRP,2),NVEC(IFRP,1),
     &                NFBAS(IFRP,0),NORB(IFRP))
      ENDDO
C
C     Get the density matrix
C
CTROND
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NVEC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),D1,D0,
     &              BUF(ICMV(I)+1),NFBAS(I,0),NVEC(I,0),
     &              1,NVEC(I,0),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF(IPRMP2NO.GE.4) THEN
        CALL TITLER('Output from DENMP2NO','*',103)
        DO 20 I = 1,NFSYM
        IF(NVEC(I,0).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
