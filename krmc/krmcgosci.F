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
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
C
C VERSION : $Revision$
C DATE    : $Date$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_setup */
      SUBROUTINE GOSCI_SETUP()
C***********************************************************************
C
C     Setup and initialize variables in dcbgosci.h
C
C     Input :
C
C     Output:
C
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dcbbos.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgosci.h"
#include "pgroup.h"
C
C
      CALL QENTER('GOSCI_SETUP')
C
C     Setup symmetry mult. table
C     --------------------------
C
C     GMULTSF/GMULTA returns:
C       NREP     : #(boson irreps) = #(fermion irreps);
C       REPNA    : strings identifying the irreps
C       MULTB    : mult. table
C
      IF (OPT_UCIBOS) THEN
C        Use boson symmetry (for spin free and Levy-Leblond)
C        or use approximate boson symmetry
         CALL GMULTSF(NREP,REPNA,MULTB(1,1))
C        hjaaj aug 2003: NOTE that this info is NOT used in this
C        version of GOSCI !!!!!!!!!!!!!
C        (MULTB is not used as it includes spin-symmetry -- see GMULTSF --
C        and we treat the spin Ms separately as in non-rel codes)
      ELSE
         CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB(1,1))
      ENDIF
C     We extend MULTB(1:2*NREP,1:2*NREP) with MULTB(1:2*NREP,0) for use in
C     GOSCI_SYMOK when we add the first spinor to a new
C     string (i.e. symmetry zero in second index means vacuum)
      DO I = 1, 2 * NREP
         MULTB(I,0) = I
      END DO
C
C     Setup symmetry
C     --------------
C
C     IRP  : IOPT_SYMMETRY is user input for CI symmetry
C     IRRP : symmetry of wave function according to:
C            1 .. NREP           : fermion irreps
C            NREP + 1 .. 2 * NREP: boson irreps
C
C             fermion               boson
C     C1   :  A                     a
C     Ci   :  AG, AU                ag, au
C
C     C2   :  1E, 2E                a, b
C     Cs   :  1E, 2E                a, b
C     C2h  :  1Eg, 2Eg, 1Eu, 2Eu    ag, bg, au, bu
C
C     C2v  :  1E, 2E                a, b
C     D2h  :  1Eg, 2Eg, 1Eu, 2Eu    ag, bg, au, bu
C
C
      IRP = IOPT_SYMMETRY
      IF (IRP.LE.0) THEN
C        ... code for include all symmetries
         IRRP = IRP
      ELSE IF ( MOD ( NAELEC, 2 ) .EQ. 0 ) THEN
C        ...boson symmetry
         IRRP = NREP + IRP
      ELSE
C        ...fermion symmetry
         IRRP = IRP
      END IF
C
C     Setup orbital symmetries
C     ------------------------
C
C     IRPAMO is the symmetry of each spinor.
C     IREPA counts fermion irrep of double group (see GMULTA)
C     (use boson irreps without spin for SPINFR .OR. LEVYLE,
C      which implies OPT_UCIBOS true, or use approx. boson symmetry
C      if OPT_UCIBOS true for relativistic calculation)
C
      II    = 0
      IREPA = 0
      DO I = 1, NFSYM
         IREPA = IREPA + 1
         DO ISTR = 1, NASH(I)
            II = II + 1
            IF (OPT_UCIBOS) THEN
               IG = IDXU2G(IASH(I)+ISTR)
               IRPAMO(II) = IBOSYM(IG)
            ELSE
               IRPAMO(II) = IREPA
            END IF
         ENDDO
         IF (NZ.NE.4.OR.SPINFR.OR.LEVYLE) IREPA = IREPA + 1
C        for NZ.eq.4: unbarred and barred in same fermion irrep !
C        (see GMULTA for def. of fermion irreps)
         DO ISTR = 1, NASH(I)
            II = II + 1
            IF (OPT_UCIBOS) THEN
               IG = IDXU2G(IASH(I)+ISTR)
               IRPAMO(II) = IBOSYM(IG)
C              to use spin symmetry here with MULTB from GMULTSF
C              we should have had NBSYM + IBOSYM(IG) here,
C              but we take care of the spin separately in GOSCI,
C              as in non-relativistic CI programs. /hjaaj aug 03
            ELSE
               IRPAMO(II) = IREPA
            END IF
         ENDDO
      END DO
C
C     Setup pointer from active orbital to a- and b-strings
C
      II = 0
      DO J = 1, NFSYM
         DO I = 1, NASH(J)
            II = II + 1
            IPTA2O(I + IASH(J) ) = II
            IPTO2A(II) = I + IASH(J)
            IPTO2B(II) = -1
         END DO
         DO I = 1, NASH(J)
            II = II + 1
            IPTB2O(I + IASH(J) ) = II
            IPTO2A(II) = -1
            IPTO2B(II) = I + IASH(J)
         END DO
      END DO
C
C
C     Print section
C     -------------
C
      IF ( IPROPT .GE. 2 ) THEN
         WRITE(LUPRI,'(A,3I3/A)')
     &      ' (GOSCI_SETUP) NREP, IRP, IRRP = ', NREP, IRP, IRRP,
     &      ' -> IRRP.lt.0 means include all symmetries'
         IF (IRRP.GT.0) WRITE(LUPRI,'(2A)')
     &      ' (GOSCI_SETUP) REPNA(IRRP)     = ', REPNA(IRRP)
         WRITE(LUPRI,'(/A)')
     &      ' Symmetry of active orbitals (unbarred, barred):'
         DO I = 1,NASHT
            IA = IPTA2O(I)
            IB = IPTB2O(I)
            IF (.NOT. OPT_UCIBOS) THEN
               WRITE(LUPRI,'(I15,2(5X,A))')
     &         I, REPNA(IRPAMO(IA)), REPNA(IRPAMO(IB))
            ELSE
               WRITE(LUPRI,'(I15,2(5X,A))')
     &         I, REP(IRPAMO(IA)), REP(IRPAMO(IB))
            END IF
         END DO
      END IF
C
      CALL QEXIT('GOSCI_SETUP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_gendet */
      SUBROUTINE GOSCI_GENDET(NDET,IDET,IASTR,IBSTR,MAXDET)
C***********************************************************************
C
C     Find determinants.
C
C     Input :
C
C     Output:
C        NDET : number of determinants
C        IDET : determinants
C
C     Written by J. Thyssen and H.J.Aa. Jenen - Jan 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgosci.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "pgroup.h"
#include "dgroup.h"
C
      DIMENSION NA(MAXASH_GOSCI), NB(MAXASH_GOSCI)
      DIMENSION IDET(2,*), IASTR(*), IBSTR(*)
      LOGICAL   GOSCI_GASOK, GOSCI_SYMOK
C
      INTEGER MSKGAS(MXMCGAS), IBTAB(0:255)
C
#include "ibtfun.h"
C
      IF ( IPROPT .GE. 3 ) THEN
         WRITE(LUPRI,'(/A)')
     &       ' (GOSCI_GENDET) Output from generation of determinants:'
         IF (OPT_UCIBOS) WRITE(LUPRI,'(15X,A)')
     &    'Using boson spinor symmetry (exact or approximate) in the CI'
         IF (IRRP .LE. 0) THEN
            WRITE(LUPRI,'(15X,A)')
     &        'Generating determinants of all symmetries'
         ELSE IF (OPT_UCIBOS) THEN
            WRITE(LUPRI,'(15X,3A,I4,A,I3)')
     &        'Generating dets. with symmetry ', REP(IRRP-1),
     &        ' and 2*MS interval:',MINMK2,' to',MAXMK2
         ELSE
            WRITE(LUPRI,'(15X,2A)')
     &        'Generating dets. with symmetry ', REPNA(IRRP)
         END IF
         WRITE(LUPRI,'(15X,A,I10)')
     &        'Max no. of determinants',MAXDET
      END IF
C
      IF (NASHT .GT. MAXASH_GOSCI) THEN
         WRITE(LUPRI,'(//A,I3,A,I3)')
     &      ' *** ERROR in GOSCI_GENDET *** Number of active'//
     &      ' orbitals (',NASHT,') must be <',MAXASH_GOSCI+1
         CALL QUIT('*** ERROR in GOSCI_GENDET ***')
      END IF
      CALL FLSHFO(LUPRI)
C
C     Setup bit table
C     ---------------
C
      CALL IBTABINI(IBTAB)
C
C     Setup GAS stuff
C     ---------------
C
      IF (NGAS .GT. 0) THEN
C
         CALL IZERO(MSKGAS,MXMCGAS)
C
         L = 1
         DO I = 1, NFSYM
            DO J = 1, NGAS
               DO K = 1, NGSH(I,J)
                  MSKGAS(J) = IBTOR(MSKGAS(J),L)
                  L = L * 2
               END DO
            END DO
         END DO
C
#ifdef MOD_DEBUG
         DO I = 1, NGAS
            write(lupri,*) 'gas bit mask before: ',i,mskgas(i)
         END DO
#endif
C
         DO I = 2, NGAS
            MSKGAS(I) = MSKGAS(I-1) + MSKGAS(I)
         END DO
C
#ifdef MOD_DEBUG
         DO I = 1, NGAS
            write(lupri,*) 'gas bit mask after: ',i,mskgas(i)
         END DO
#endif
C
      END IF
C
C     Find min and max number of a & b electrons
C     ------------------------------------------
C
D     write(lupri,*) 'maxmk2, minmk2',maxmk2,minmk2
      MAXAEL = (MAXMK2 + NAELEC) / 2
      MAXAEL = MIN(MAXAEL,NASHT)
C
      MAXBEL = (NAELEC - MINMK2) / 2
      MAXBEL = MIN(MAXBEL,NASHT)
C
      MINAEL = NAELEC - MAXBEL
      MINBEL = NAELEC - MAXAEL
C
      MINOCC = MIN(MINAEL,MINBEL)
      MAXOCC = MAX(MAXAEL,MAXBEL)
C
C
C     Generate determinants:
C     -----------------------
C
D     write(lupri,*) 'IAEL   IBEL'
      NDET = 0
      DO 10 IAEL = MINAEL, MAXAEL
         IBEL = NAELEC - IAEL
         MK2 = IAEL - IBEL
         IMK2(MK2) = NDET
D        write(lupri,'(I4,I7)') IAEL,IBEL
D        call flshfo(lupri)
C
         CALL GOSCI_GENSTR(IAEL,IASTR,NASTR,MSKGAS,MAXDET)
         IF (NASTR .EQ. 0) GOTO 10
         CALL GOSCI_GENSTR(IBEL,IBSTR,NBSTR,MSKGAS,MAXDET)
         IF (NBSTR .EQ. 0) GOTO 10
D        write(lupri,'(12X,I8,I11,I11)') NASTR,NBSTR,NDET
D        call flshfo(lupri)
C
C HJ TODO : the next can easily be optimized by not redoing the IA
C           for each IB in GOSCI_GASOK and GOSCI_SYMOK /aug 03
C
         DO IA = 1, NASTR
            DO IB = 1, NBSTR
               IF (GOSCI_GASOK(IASTR(IA),IBSTR(IB),MSKGAS,IBTAB)) THEN
                  IF (GOSCI_SYMOK(IASTR(IA),IBSTR(IB))) THEN
                     NDET = NDET + 1
                     IF (NDET .GE. MAXDET) THEN
                        WRITE(LUPRI,'(2A,I9,A)')
     &                       ' *** ERROR in GOSCI_GENDET ***',
     &                       'Maximum number of determinants reached!',
     &                       MAXDET,
     &                       'Allocate more memory!'
                        CALL QUIT('*** ERROR in GOSCI_GENDET ***')
                     END IF
                     IDET(1,NDET) = IASTR(IA)
                     IDET(2,NDET) = IBSTR(IB)
                  END IF
               END IF
            END DO
         END DO
D        write(LUPRI,'(T32,I11)') NDET
D        call flshfo(LUPRI)
         NMK2(MK2) = NDET - IMK2(MK2)
 10   CONTINUE
C
C
C     Write statistics about determinants
C
      WRITE(LUPRI,9100)
      DO MK2 = MINMK2, MAXMK2
         IF (NMK2(MK2) .GT. 0) THEN
            IF (MOD(MK2,2) .EQ. 0) THEN
               WRITE(LUPRI,9101) MK2/2,NMK2(MK2)
            ELSE
               WRITE(LUPRI,9102) MK2,NMK2(MK2)
            END IF
         END IF
      END DO
      WRITE(LUPRI,9103) NDET
 9100 FORMAT(//' Distribution of determinants wrt. MK'
     &        /5X,'MK',T28,'#dets'/1X,38('-'))
 9101 FORMAT(I8,T23,I10)
 9102 FORMAT(I6,'/2',T23,I10)
 9103 FORMAT(1X,38('-') / 5X,'all',T23,I10)
C
      CALL FLSHFO(LUPRI)
C
      IF (IPROPT .GE. 10 .OR. NDET .LE. 50) THEN
         WRITE(LUPRI,9001)
         DO MK2 = MINMK2, MAXMK2
            IF (NMK2(MK2) .GT. 0) THEN
               IF (MOD(MK2,2) .EQ. 0) THEN
                  WRITE(LUPRI,9004) NMK2(MK2),MK2/2
               ELSE
                  WRITE(LUPRI,9005) NMK2(MK2),MK2
               END IF
            END IF
C
            DO I = 1, NMK2(MK2)
               IDA = IDET(1,I+IMK2(MK2))
               DO J=1,NASHT
                  NA(J)=0
                  IF (IBTAND (IDA, IBTSHL (1, J - 1)) .NE. 0)
     &                 NA(J)=1
               END DO
               IDB = IDET(2,I+IMK2(MK2))
               DO J=1,NASHT
                  NB(J)=0
                  IF (IBTAND (IDB, IBTSHL (1, J - 1)) .NE. 0)
     &                 NB(J)=1
               END DO
               WRITE(LUPRI,9002) I+IMK2(MK2),(NA(J),J=1,NASHT)
               WRITE(LUPRI,9003)             (NB(J),J=1,NASHT)
            END DO
         END DO
         CALL FLSHFO(LUPRI)
      END IF
 9001 FORMAT(//' Det.',T20,'String',/1X,77('-'))
 9002 FORMAT(I7,T20,100I1)
 9003 FORMAT(   T20,100I1)
 9004 FORMAT(I10,' determinants with MK =',I3,' :')
 9005 FORMAT(I10,' determinants with MK =',I3,'/2 :')
C
C
      RETURN
 1010 FORMAT(I6,I10,32I2)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_gasok */
      LOGICAL FUNCTION GOSCI_GASOK(IA,IB,MSKGAS,IBTAB)
C***********************************************************************
C
C     Check that determinant is allowed for GASSCF calculations.
C
C     Input :
C       I : determinant
C
C     Output:
C
C     Written by J. Thyssen - Feb 11 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgosci.h"
#include "dcbopt.h"
#include "dcborb.h"
C
      DIMENSION MSKGAS(*), IBTAB(*)
C
#include "ibtfun.h"
C
      GOSCI_GASOK = .TRUE.
      IF (NGAS .EQ. 0) RETURN
C
      DO IGAS = 1, NGAS
         IGA = IBTAND(MSKGAS(IGAS),IA)
         IGB = IBTAND(MSKGAS(IGAS),IB)
         CALL FNBITS(NA,IGA,IBTAB)
         CALL FNBITS(NB,IGB,IBTAB)
         NG = NA + NB
         IF ( NG .LT. NGSP(1,IGAS) .OR. NG .GT. NGSP(2,IGAS)) THEN
            GOSCI_GASOK = .FALSE.
            RETURN
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ibtabini */
      SUBROUTINE IBTABINI(IBTAB)
C***********************************************************************
C
C     Initialize table with IBTAB(I) = NBITS(I)
C
C     Input : none
C
C     Output:
C       IBTAB
C
C     Written by J. Thyssen - Feb 27 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IBTAB(0:255)
C
#include "ibtfun.h"
C
      DO I = 0, 255
         ID = I
         NB = 0
         DO J = 1, 8
            NB = NB + IBTAND(1,ID)
            ID = ID / 2
         END DO
         IBTAB(I) = NB
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_krlist */
      SUBROUTINE GOSCI_KRLIST(NDET,IDET,KRLIST,IPRINT)
C***********************************************************************
C
C     For each determinant find Kramers partner.
C
C     Input :
C        NDET : number of determinants
C        IDET : determinants
C
C     Output:
C        KRLIST:
C
C     Written by J. Thyssen - Jan 23 2001
C     Revised Sep 2003 hjaaj (was with old, obsolete IDET format)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IDET(2,*), KRLIST(*)
C
      CALL IZERO(KRLIST,NDET)
C
      DO 10 I = 1, NDET
         IF (KRLIST(I) .NE. 0) GOTO 10
         IDA = IDET(1,I)
         IDB = IDET(2,I)
         DO 20 J = NDET, I, -1
            IF (KRLIST(J) .NE. 0) GOTO 20
            IF (IDA .EQ. IDET(2,J) .AND. IDB .EQ. IDET(1,J)) THEN
               KRLIST(I) = J
               KRLIST(J) = I
               GOTO 30
            END IF
 20      CONTINUE
         WRITE(LUPRI,*) 'WARNING: No Kramers partner found for det.',I
 30      CONTINUE
 10   CONTINUE
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Output from GOSCI_KRLIST',-1)
         WRITE(LUPRI,'(T10,A,T30,A)') 'Determinant','Kramers partner'
         DO I = 1, NDET
            WRITE(LUPRI,'(T10,I9,T30,I9)') I,KRLIST(I)
         END DO
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ngosci_nbits */
      INTEGER FUNCTION NGOSCI_NBITS(JDET,IBTAB)
C***********************************************************************
C
C     Find number if bits in integer JDET.
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
C
      DIMENSION IBTAB(0:255)
C
      INTEGER*1 ID1(4)
      INTEGER*4 ID4
      EQUIVALENCE (ID1(1),ID4)
C
C     Anyone who can find with a faster routine?
C
      ID4 = JDET
D     write(lupri,*) 'jdet',jdet,' and id4',id4
D     write(lupri,*) 'id1',id1
      NGOSCI_NBITS
     &     = IBTAB(ID1(1))
     &     + IBTAB(ID1(2))
     &     + IBTAB(ID1(3))
     &     + IBTAB(ID1(4))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_symok */
      LOGICAL FUNCTION GOSCI_SYMOK(IA,IB)
C***********************************************************************
C
C     Find symmetry of determinant (IA,IB).
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
#include "dcbopt.h"
C
#include "ibtfun.h"
C
      IF (IOPT_SYMMETRY .LE. 0) THEN
C        ... include all symmetries
         GOSCI_SYMOK = .TRUE.
      ELSE IF (.NOT. OPT_UCIBOS) THEN
C
C        The normal symmetry scheme.
C
         JRRP = 0
         ID=IA
         DO I=1,NASHT
            IF (MOD(ID,2).EQ.1) THEN
               JRRP = MULTB(IRPAMO(IPTA2O(I)),JRRP)
            ENDIF
            ID=ID/2
         END DO
         ID=IB
         DO I=1,NASHT
            IF (MOD(ID,2).EQ.1) THEN
               JRRP = MULTB(IRPAMO(IPTB2O(I)),JRRP)
            ENDIF
            ID=ID/2
         END DO
         GOSCI_SYMOK = (JRRP .EQ. IRRP)
      ELSE
C
C        Boson symmetry (spin free or Levy-Leblond).
C        (do not use MULTB which includes spin-symmetry -- see GMULTSF --
C         as we treat the spin Ms separately as in non-rel codes /hjaaj aug 03)
C
         JRRP = 0
         ID = IA
         DO I = 1, NASHT
            IF ( MOD(ID,2) .EQ. 1) THEN
               JRRP = IBTXOR(IRPAMO(IPTA2O(I)),JRRP)
            END IF
            ID = ID / 2
         END DO
         ID = IB
         DO I = 1, NASHT
            IF ( MOD(ID,2) .EQ. 1) THEN
               JRRP = IBTXOR(IRPAMO(IPTB2O(I)),JRRP)
            END IF
            ID = ID / 2
         END DO
         GOSCI_SYMOK = (JRRP .EQ. IOPT_SYMMETRY-1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_startvec */
      SUBROUTINE GOSCI_STARTVEC(IDET,CMO,CREF,ECI,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors from GOSCI by diagonalizing the
C     full CI matrix.
C
C     Input : CMO - MO coefficents
C
C     Output: CREF - CI vector
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), CREF(NZCONF,*), CMO(*), IDET(*)
C
      CALL QENTER('GOSCI_STARTVEC')
C
      KFRSAV = KFREE
      MZ = MIN(2,NZ)
      N2ZCONF = NZCONF*NZCONF
C
      CALL MEMGET('REAL',KEIGVL,NZCONF,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,N2ZCONF*MZ,WORK,KFREE,LFREE)
C
C
C     Do full CI
C
      CALL GOSCI_FULLCI(IDET,CMO,WORK(KEIGVL),WORK(KEIGVC),
     &                  ECORE,WORK,KFREE,LFREE)
C
C     Get CI vector:
C     --------------
C
C     The CI vector is stored in colunms; the lowest eigenstate is in
c     the first column.
C
      IOFF = KEIGVC - 1 + ( IOPT_ISTATE - 1 ) * NZCONF
      IF ( NZ .EQ. 1 ) THEN
         DO I = 1,NZCONF
            CREF(I,1) = WORK(IOFF + I)
         END DO
      ELSE
         DO I = 1,NZCONF
            CREF(I,1) = WORK(IOFF + I )
            CREF(I,2) = WORK((IOFF + N2ZCONF) + I )
         END DO
      END IF
C
C DEBUG DEBUG DEBUG
c      call dzero(cref,nzconfq)
c      if (nz .eq. 1) then
c      cref(2,1) = 1.0d0
c      cref(2,2) = 1.0d0/sqrt(2.0d0)
c      else
c         CREF(2,2) = 1.0D00
c      endif
c      DO I = 2, NZCONF
c         CREF(I,1) = 0.0D00
c      END DO
c      IF ( NZCONF .GE. 2) THEN
c         CREF(1,1) = 0.8
c         CREF(2,1) = 0.6
c      END IF
c      CREF(1,1) = 0.8D00
c      IF ( NZCONF .GE. 2) CREF(2,1) = 0.6D00
c      DO I = 3, NZCONF
c         CREF(I,1) = 0.0D00
c      END DO
C
C
C     Get CI energy
C
      ECI = WORK(KEIGVL + IOPT_ISTATE - 1)
      WRITE(LUPRI,'(/A,F20.8)')
     &     ' (GOSCI_STARTVEC) Core energy    = ',ECORE
      WRITE(LUPRI,'(A,F20.8)')
     &     '                  Active energy  = ',ECI-ECORE
      WRITE(LUPRI,'(A,F20.8)')
     &     '                  Total energy   = ',ECI
C
C     Check norm of CI vector
C
      CALL RCINORM(CREF,10)
C
      IF ( IPROPT .GE. 5 ) THEN
         WRITE(LUPRI,'(A)') ' (GOSCI_STARTVEC) CI start vector'
         CALL RPRCI(CREF,NZCONF,NZ,OPT_CIPROGRAM, THRPCI,LUPRI)
      END IF
C
      CALL MEMREL('GOSCI_STARTVEC',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GOSCI_STARTVEC')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_fullci */
      SUBROUTINE GOSCI_FULLCI(IDET,CMO,EIGVAL,EIGVEC,ECORE,
     &                        WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors from GOSCI by diagonalizing the
C     full CI matrix.
C
C     Input : IDET   - determinant information
C             CMO    - MO coefficents
C
C     Output: EIGVAL - CI eigenvalues
C             EIGVEC - CI eigenvectors
C             ECORE  - core energy
C
C     Extracted from GOSCI_STARTVEC July 2002 /hjaaj.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), EIGVAL(*), EIGVEC(*), CMO(*), IDET(*)
C
      CALL QENTER('GOSCI_FULLCI')
C
      KFRSAV = KFREE
      MZ = MIN(2,NZ)
      N2ZCONF = NZCONF * NZCONF
      N2ZCONFQ = N2ZCONF * MZ
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET('REAL',KFCACM,(2*NASHT)*(2*NASHT)*MZ,
     &     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Calculate FCmo (TOFILE false: do not write to LUKRM3) :
C
      CALL RGETFC(.FALSE.,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
      CALL MEMREL('GOSCI_FULLCI after RGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * MZ, WORK,KFREE,LFREE)
C
      LH2AC = NASHT*NASHT*NNASHX*NZ*3
      CALL MEMGET('REAL',KH2AC,LH2AC,WORK,KFREE,LFREE)
C
C     Read integrals from 4IND1XXXX0:
C
      CALL RGETH21(.FALSE.,.TRUE.,DUMMY,DUMMY,WORK(KH2AC),
     &     WORK(KFREE),LFREE,IPROPT)
C
C     Transform integrals to Molfdir format:
C
      CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPROPT)
      CALL MEMREL('GOSCI_FULLCI.H2AC',WORK,1,KH2AC,KFREE,LFREE)
C
C     Allocate memory for matrices:
C
      CALL MEMGET('REAL',KHR,N2ZCONFQ,WORK,KFREE,LFREE)
      KHI = KHR + N2ZCONF
C
C     Construct Hamiltonian:
C     ----------------------
C
      CALL GOSCI_MAKEH(.TRUE.,1,IDET,WORK(KHR),WORK(KHI),
     &                 WORK(KFCACM),WORK(KMUUUU))
C
C
C     Diagonalize Hamiltonian:
C     ------------------------
C
      IF ( IPROPT .GE. 20 ) THEN
         WRITE(LUPRI,'(A)') ' (GOSCI_FULLCI) CI matrix - real part:'
         CALL OUTPUT(WORK(KHR),1,NZCONF,1,NZCONF,NZCONF,NZCONF,1,LUPRI)
        IF (MZ .GT. 1) THEN
         WRITE(LUPRI,'(A)') ' (GOSCI_FULLCI) CI matrix - imag part:'
         CALL OUTPUT(WORK(KHI),1,NZCONF,1,NZCONF,NZCONF,NZCONF,1,LUPRI)
        END IF
      END IF
      CALL FLSHFO(LUPRI)
C
C
      IF (MZ .EQ. 1) THEN
         IJOB = 1
         IORDER = 1
         IPACK = 0
         CALL RSJACO(NZCONF,NZCONF,NZCONF,WORK(KHR),EIGVAL,
     &        IJOB,IORDER,IPACK,EIGVEC,WORK(KFREE),LFREE)
      ELSE
         MATZ = 1
         CALL QDIAG(2,NZCONF,WORK(KHR),NZCONF,NZCONF,
     &        EIGVAL,MATZ,EIGVEC,NZCONF,NZCONF,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in GOSCI_FULLCI ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in GOSCI_FULLCI *** QDIAG failed')
         END IF
      END IF
C
C     Add core energy
C
      DO I = 1,NZCONF
         EIGVAL(I) = EIGVAL(I) + ECORE
      END DO
C
      IF ( IPROPT .GE. 5 ) THEN
         WRITE(LUPRI,'(/A)') ' (GOSCI_FULLCI) CI eigenvalues'
         CALL OUTPUT(EIGVAL,1,1,1,NZCONF,1,NZCONF,1,LUPRI)
       IF ( IPROPT .GE. 20 .OR. NZCONF .LE. 12) THEN
         WRITE(LUPRI,'(/A)')
     &   ' (GOSCI_FULLCI) CI eigenvectors (real part)'
         CALL OUTPUT(EIGVEC,1,NZCONF,1,NZCONF,
     &        NZCONF,NZCONF,1,LUPRI)
         IF (NZ .GT. 1) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (GOSCI_FULLCI) CI eigenvectors (imag part)'
            CALL OUTPUT(EIGVEC(1+N2ZCONF),1,NZCONF,1,NZCONF,
     &           NZCONF,NZCONF,1,LUPRI)
         END IF
       END IF
      END IF
C
      CALL MEMREL('GOSCI_FULLCI',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GOSCI_FULLCI')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gosci_ireps */
      SUBROUTINE GOSCI_IREPS (NZ,NFSYM,NSYMRPA,REPA)
C***********************************************************************
C
C     Rip-off from GMULTA from MOLTRA./jth-03072000
C
C     Hard-wired multiplication table for the abelian subgroups of D2h
C     Defines representation names as well.
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER*4 REPA(8)
C
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
         NSYMRPA = 1
         REPA(1) = '   A'
         REPA(2) = '   a'
      ELSEIF (NFSYM.EQ.1.AND.NZ.LE.2) THEN
         NSYMRPA = 2
         REPA(1) = '  1E'
         REPA(2) = '  2E'
         REPA(3) = '   a'
         REPA(4) = '   b'
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRPA = 2
         REPA(1) = '  AG'
         REPA(2) = '  AU'
         REPA(3) = '  ag'
         REPA(4) = '  au'
      ELSEIF (NFSYM.EQ.2.AND.NZ.LE.2) THEN
         NSYMRPA = 4
         REPA(1) = ' 1Eg'
         REPA(2) = ' 2Eg'
         REPA(3) = ' 1Eu'
         REPA(4) = ' 2Eu'
         REPA(5) = '  ag'
         REPA(6) = '  bg'
         REPA(7) = '  au'
         REPA(8) = '  bu'
      ELSE
         CALL QUIT('Multiplication table is not defined')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_rmakdm */
      SUBROUTINE GOSCI_RMAKDM(CALCDV,CALCPV,IDET,CL,CR,DV,PV,IPRINT)
C***********************************************************************
C
C     Construct one- and two-electron active density matrix
C
C     D_{pq} = < CL | X^+_{pq} | CR >
C
C     P_{pqrs} = < CL | x^+_{pq,rs} | CR >
C
C     Input:
C        CL, CR  - CI vectors
C
C     Output:
C        DV      - one electron active density matrix
C        PV      - two electron active density matrix
C
C     DV has dimension (2*NASHT)**2 * MIN(NZ,2).
C     PV has dimensions (2*NASHT)**4 * MIN(NZ,2). Triangular packing
C     would be better, but whatever...
C     The code is now very straightforward.
C
C     Written by J. Thyssen - Oct 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcbgosci.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
C
      DIMENSION CL(NZCONF,NZ),CR(NZCONF,NZ),DV(2*NASHT,2*NASHT,2)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,2)
      DIMENSION IDET(*)
C
      LOGICAL CALCDV, CALCPV
C
C
      CALL QENTER('GOSCI_RMAKDM')
C
      IF (.NOT. CALCDV) THEN
         WRITE(LUPRI,'(/A)')
     $        ' *** ERROR in GOSCI_RMAKDM *** ' //
     $        'CALCDV false is not implemented.'
         CALL QUIT('*** ERROR in GOSCI_RMAKDM ***')
      END IF
C
C     Calculate density matrices
C     --------------------------
C
      CALL GOSCI_MAKETDM(CALCPV,IDET,CL,CR,DV,PV)
C
C     Print section
C     -------------
C
      IF ( IPRINT .GE. 10 ) THEN
C
C        Print DV
C
         IF (CALCDV) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GOSCI_RMAKDM) real(DV) in Molfdir basis'
            CALL OUTPUT(DV(1,1,1),1,2*NASHT,1,2*NASHT,
     &           2*NASHT,2*NASHT,1,LUPRI)
            IF ( NZ .GE. 2 ) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (GOSCI_RMAKDM) imag(DV) in Molfdir basis'
               CALL OUTPUT(DV(1,1,2),1,2*NASHT,1,2*NASHT,
     &              2*NASHT,2*NASHT,1,LUPRI)
            END IF
         END IF
C
C        Print PV
C
         IF (CALCPV) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GOSCI_RMAKDM) PV in Molfdir basis'
            WRITE(LUPRI,'(4X,A,14X,A,10X,A)')
     $           'Index','Real(PV)','Imag(PV)'
            DO IU = 1,2*NASHT
               DO IV = 1,2*NASHT
                  DO IX = 1,2*NASHT
                     DO IY = 1,2*NASHT
                        IF ( NZ .EQ. 1 ) THEN
                           IF ( ABS ( PV(IU,IV,IX,IY,1) ) .GT. D0 ) THEN
                              WRITE(LUPRI,'(4I3,2X,1P,D17.10)')
     $                             IU,IV,IX,IY,PV(IU,IV,IX,IY,1)
                           END IF
                        ELSE
                           IF ( ABS ( PV(IU,IV,IX,IY,1) +
     &                          PV(IU,IV,IX,IY,2) ) .GT. D0 ) THEN
                              WRITE(LUPRI,
     &                             '(4I3,2X,1P,D17.10,4X,D17.10)')
     $                             IU,IV,IX,IY,
     $                             PV(IU,IV,IX,IY,1),PV(IU,IV,IX,IY,2)
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END IF
C
      CALL QEXIT('GOSCI_RMAKDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_maketdn */
      SUBROUTINE GOSCI_MAKETDM(DOPV,IDET,CL,CR,DV,PV)
C***********************************************************************
C
C     Calculate density matrices: <CL|Pq|CR>
C
C     Input:
C        CL, CR - CI vectors
C        DOPV - calculate 2-electron density matrix
C
C     Output:
C        DV, PV - density matrices
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL   DOPV
      DIMENSION CL(*), CR(*), DV(*), PV(*)
      DIMENSION IDET(*)
C
      CALL QENTER('GOSCI_MAKETDM')
C
      CALL GOSCI_CI(
     &     .FALSE.,DUMMY,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,CL,CR,DV,PV,
     &     DOPV,IDET)
C     CALL GOSCI_CI(
C    &     DODIA,DIAG,
C    &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C    &     DODEN,CL,CR,DV,PV,
C    &     DO2EL,IDET)
C
      CALL QEXIT('GOSCI_MAKETDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_pmat */
      SUBROUTINE GOSCI_PMAT(DOPV,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &     DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C***********************************************************************
C
C     Calculate contributions determinant I and J.
C
C     Input:
C        DOPV   - calculate PV
C        DOIMAG - calculate imaginary contributions
C        IT     - the number of different determinants in I and J
C        IOCCI, JOCCI - the strings for determinants I and J
C        K,L,M,N - IT = 0: not used
C                  IT = 2: K and M is the orbitals where I and J differ
C                  IT = 4: K,L,M,N  ditto.
C        DMRIJ, DMRJI - real part of cl^*_{i} cr_{j} and cl^*_{j} cr_{i}
C        DMIIJ, DMIJI - imag part of ditto.
C
C     Output:
C        DV, PV: density matrices
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgosci.h"
C
      LOGICAL   DOPV, DOIMAG
      DIMENSION IOCCI(2*MAXASH_GOSCI), JOCCI(2*MAXASH_GOSCI)
      DIMENSION DV(2*NASHT,2*NASHT,*)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
#include "ibtfun.h"
C
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
         DO 30 II = 1, NAELEC
            I = IOCCI(II)
            DV(I,I,1) = DV(I,I,1) + DMRIJ
            IF (DOIMAG) DV(I,I,2) = DV(I,I,2) + DMIIJ
            IF (DOPV) THEN
               DO 20 JJ = II + 1, NAELEC
                  J = IOCCI(JJ)
C---------------------
C                    H = H + (I, I, J, J)
C                    H = H - (I, J, J, I)
C---------------------
                  PV(J,J,I,I,1) = PV(J,J,I,I,1) + DMRIJ
                  PV(I,J,J,I,1) = PV(I,J,J,I,1) - DMRIJ
C---------------------
                  IF (DOIMAG) THEN
                     PV(J,J,I,I,2) = PV(J,J,I,I,2) + DMIIJ
                     PV(I,J,J,I,2) = PV(I,J,J,I,2) - DMIIJ
                  END IF
 20            CONTINUE
            END IF
 30      CONTINUE
         RETURN
      END IF
C
      IF(IT.EQ.2) THEN
C
C        ===================================
C        DETERMINANTS DIFFER BY ONE FUNCTION K -> M
C        ===================================
C
         IPHASE = 0
         DO IELEC = 1,NAELEC
            IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
         END DO
         IPHASE = (-1)**IPHASE
C
         DMRIJ1 = IPHASE * DMRIJ
         DMRJI1 = IPHASE * DMRJI
         DV(K,M,1) = DV(K,M,1) + DMRIJ1
         DV(M,K,1) = DV(M,K,1) + DMRJI1
         IF (DOIMAG) THEN
            DMIIJ1 = IPHASE * DMIIJ
            DMIJI1 = IPHASE * DMIJI
            DV(K,M,2) = DV(K,M,2) + DMIIJ1
            DV(M,K,2) = DV(M,K,2) + DMIJI1
         END IF
C
         IF (DOPV) THEN
            DO IELEC = 1, NAELEC
               I = IOCCI(IELEC)
               IF (I .NE. K) THEN
C---------------------
C                 H = H + (I, I, K, M)
C                 H = H - (I, M, K, I)
C---------------------
                  PV(I,I,K,M,1) = PV(I,I,K,M,1) + DMRIJ1
                  PV(I,M,K,I,1) = PV(I,M,K,I,1) - DMRIJ1
                  PV(I,I,M,K,1) = PV(I,I,M,K,1) + DMRJI1
                  PV(I,K,M,I,1) = PV(I,K,M,I,1) - DMRJI1
C---------------------
                  IF (DOIMAG) THEN
                     PV(I,I,K,M,2) = PV(I,I,K,M,2) + DMIIJ1
                     PV(I,M,K,I,2) = PV(I,M,K,I,2) - DMIIJ1
                     PV(I,I,M,K,2) = PV(I,I,M,K,2) + DMIJI1
                     PV(I,K,M,I,2) = PV(I,K,M,I,2) - DMIJI1
                  END IF
               END IF
            END DO
         END IF
         RETURN
      END IF
C
      IF(IT.EQ.4 .AND. DOPV) THEN
C
C        ====================================
C        DETERMINANTS DIFFER BY TWO FUNCTIONS  K,L -> M,N
C        ====================================
C
         IPHASE = 0
         DO IELEC = 1,NAELEC
            IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
            IF (IOCCI(IELEC) .EQ. L) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. N) IPHASE = IPHASE + IELEC
         END DO
         IPHASE = (-1)**IPHASE
C
         DMRIJ1 = IPHASE * DMRIJ
         DMRJI1 = IPHASE * DMRJI
         IF (DOIMAG) THEN
            DMIIJ1 = IPHASE * DMIIJ
            DMIJI1 = IPHASE * DMIJI
         END IF
C
C---------------------
C        H = H + (K, M, L, N)
C---------------------
         PV(K,M,L,N,1) = PV(K,M,L,N,1) + DMRIJ1
         PV(K,N,L,M,1) = PV(K,N,L,M,1) - DMRIJ1
         PV(M,K,N,L,1) = PV(M,K,N,L,1) + DMRJI1
         PV(M,L,N,K,1) = PV(M,L,N,K,1) - DMRJI1
C---------------------
C        H = H - (K, N, L, M)
C---------------------
         IF (DOIMAG) THEN
            PV(K,N,L,M,2) = PV(K,N,L,M,2) + DMIIJ1
            PV(K,M,L,N,2) = PV(K,M,L,N,2) - DMIIJ1
            PV(M,K,N,L,2) = PV(M,K,N,L,2) + DMIJI1
            PV(M,L,N,K,2) = PV(M,L,N,K,2) - DMIJI1
         END IF
         RETURN
      END IF
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C
      CALL QUIT('GOSCI_PMAT called with IT .gt. 4 !')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_maketdn */
      SUBROUTINE GOSCI_MAKEH(DO2EL,IHRM,IDET,HR,HI,HCORE,VMUUUU)
C***********************************************************************
C
C     Construct Hamiltonian for use in gosci.
C
C     Input:
C        DO2EL  - true: both one and two electron integrals
C                 false: only one electron integrals
C        IHRM   - .gt. 0: Hermitian operator
C                 .lt. 0: Anti-hermitian operator
C                 .eq. 0: non-Hermitian operator
C                         (or debug to calculate full CI Hamiltonian)
C        HCORE  - one electron integrals
C        VMUUUU - two electron integrals.
C
C     Output:
C        HR - real part of CI Hamiltonian
C        HI - imag part of CI Hamiltonian
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
C
      DIMENSION HR(NZCONF,NZCONF),HI(NZCONF,NZCONF)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
      DIMENSION HCORE(*), IDET(2,*)
      LOGICAL   DO2EL,DOIMAG
      DIMENSION IBTAB(0:255)
C
      DIMENSION IOCCI(2*MAXASH_GOSCI), JOCCI(2*MAXASH_GOSCI)
      LOGICAL HFULL
C
#include "ibtfun.h"
C
      CALL IBTABINI(IBTAB)
      HFULL = IHRM .EQ. 0
      DOIMAG  = NZ.GT.1
      CALL DZERO(HR,NZCONF*NZCONF)
      IF (DOIMAG) CALL DZERO(HI,NZCONF*NZCONF)
C
      IF (DO2EL) THEN
C        ... one and two electron terms
         ITMAX = 4
      ELSE
C        ... only one electron terms
         ITMAX = 2
      END IF
C
      IDAOLD = -1
      DO I = 1, NZCONF
C
C        Use that outer loop in gosci_gendet is over A-strings,
C        inner loop is over B-strings
C
         IDA = IDET(1,I)
         IF (IDA .NE. IDAOLD) THEN
            IDAOLD = IDA
            IELECA = 0
            DO K = 1, NASHT
               IF (IBTAND (IDA, IBTSHL (1, K - 1)) .NE. 0) THEN
                  IELECA = IELECA + 1
                  IOCCI(IELECA) = IPTA2O(K)
               END IF
            END DO
         END IF
C
         IDB = IDET(2,I)
         IELEC = IELECA
         DO K = 1, NASHT
            IF (IBTAND (IDB, IBTSHL (1, K - 1)) .NE. 0) THEN
               IELEC = IELEC + 1
               IOCCI(IELEC) = IPTB2O(K)
            END IF
         END DO
         IF (IELEC .NE. NAELEC)
     &      CALL QUIT('GOSCI_MAKEH: IELEC.ne.NAELEC')
C
         IF (HFULL) THEN
            JEND = NZCONF
         ELSE
chj aug 2002: code generalized to non-hermitian operators
chj         CALL GOSCI_HMAT(DO2EL,.FALSE.,0,IOCCI,IOCCI,0,0,0,0,
chj  &           HR(I,I),DUMMY,HCORE,VMUUUU)
            CALL GOSCI_HMAT(DO2EL,DOIMAG,0,IOCCI,IOCCI,0,0,0,0,
     &           HR(I,I),HI(I,I),HCORE,VMUUUU)
            JEND = I - 1
         END IF
         JDAOLD = -1
         DO J = 1, JEND
C
C           Use that outer loop in gosci_gendet is over A-strings,
C           inner loop is over B-strings
C
            JDA = IDET(1,J)
            IF (JDA .NE. JDAOLD) THEN
               JDAOLD = JDA
               IDAV = IBTXOR (IDA, JDA)
               CALL FNBITS(NDAV,IDAV,IBTAB)
C              Can we get a contribution from <I| H |J> ?
               IF (NDAV .LE. ITMAX) THEN
                  JELECA = 0
                  IAI = IBTAND (IDAV, IDA)
                  KA = 0
                  IAJ = IBTAND (IDAV, JDA)
                  MA = 0
                  DO II = 1, NASHT
                     IIBIT = IBTSHL (1, II - 1)
                     IF (IBTAND (JDA, IIBIT) .NE. 0) THEN
                        JELECA = JELECA + 1
                        JOCCI(JELECA) = IPTA2O(II)
                     END IF
                     IF (IBTAND (IAI, IBTSHL (1, II - 1)) .NE. 0) THEN
                        L  = KA
                        KA = IPTA2O(II)
                     END IF
                     IF (IBTAND (IAJ, IBTSHL (1, II - 1)) .NE. 0) THEN
                        N  = MA
                        MA = IPTA2O(II)
                     END IF
                  END DO
               END IF
            END IF
C
C           Any reason to check the B-string difference?
C
            IF (NDAV .LE. ITMAX) THEN
C
C           Yes, do it.
C
               JDB = IDET(2,J)
               IDBV = IBTXOR (IDB, JDB)
               CALL FNBITS(NDBV,IDBV,IBTAB)
C
               IT = NDAV + NDBV
C
C           Can this give a non-zero contribution ?
C
            IF ( IT .LE. ITMAX ) THEN
C
C           Yes, calculate the contribution.
C
               IELEC = JELECA
               DO K = 1, NASHT
                  IF (IBTAND (JDB, IBTSHL (1, K - 1)) .NE. 0) THEN
                     IELEC = IELEC + 1
                     JOCCI(IELEC) = IPTB2O(K)
                  END IF
               END DO
C
               IBI = IBTAND (IDBV, IDB)
               IBJ = IBTAND (IDBV, JDB)
               K = KA
               M = MA
               DO II = 1, NASHT
                  IF (IBTAND (IBI, IBTSHL (1, II - 1)) .NE. 0) THEN
                     L = K
                     K = IPTB2O(II)
                  END IF
                  IF (IBTAND (IBJ, IBTSHL (1, II - 1)) .NE. 0) THEN
                     N = M
                     M = IPTB2O(II)
                  END IF
               END DO
C
               CALL GOSCI_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &              HR(I,J),HI(I,J),HCORE,VMUUUU)
               IF( IHRM .GT. 0 ) THEN
                  HR (J,I) =  HR (I, J)
                  IF (DOIMAG) HI (J,I) = -HI (I, J)
               ELSE IF( IHRM .LT. 0 ) THEN
                  HR (J,I) = -HR (I, J)
                  IF (DOIMAG) HI (J,I) =  HI (I, J)
               END IF
            END IF
C           ... for IF ( IT .LE. ITMAX ) THEN
            END IF
C           ... for IF (NDAV .LE. ITMAX) THEN
         ENDDO
      ENDDO
C
c     write(LUPRI,*) 'real part of HMAT in GOSCI_MAKEH ='
c     call output(hr,1,nzconf,1,nzconf,nzconf,nzconf,1,lupri)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_hmat */
      SUBROUTINE GOSCI_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &     HR,HI,HCORE,VMUUUU)
C***********************************************************************
C
C     Calculate Hamiltonian matrix element
C
C     Written by H. J. Aa. Jensen & J. Thyssen - Dec 30 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
C
      LOGICAL   DO2EL,DOIMAG
      DIMENSION IOCCI(2*MAXASH_GOSCI), JOCCI(2*MAXASH_GOSCI)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
      DIMENSION HCORE(2*NASHT,2*NASHT,*)
C
#include "ibtfun.h"
C
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
         HR = D0
         IF (DOIMAG) HI = D0
         DO 30 II = 1, NAELEC
            I = IOCCI(II)
            HR = HR + HCORE(I,I,1)
            IF (DOIMAG) HI = HCORE(I,I,2)
           IF (DO2EL) THEN
            DO 20 JJ = II + 1, NAELEC
               J = IOCCI(JJ)
C---------------------
C              H = H + (I, I, J, J) - (I, J, J, I)
C---------------------
               HR = HR + VMUUUU(J,J,I,I,1) - VMUUUU(I,J,J,I,1)
               IF (DOIMAG)
     &         HI = HI + VMUUUU(J,J,I,I,2) - VMUUUU(I,J,J,I,2)
 20         CONTINUE
           END IF
 30      CONTINUE
         RETURN
      END IF
C
      IF(IT.EQ.2) THEN
C
C        ===================================
C        DETERMINANTS DIFFER BY ONE FUNCTION K -> M
C        ===================================
C
         IPHASE = 0
         DO IELEC = 1,NAELEC
            IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
         END DO
         IPHASE = (-1)**IPHASE
C
         HR = HCORE(K,M,1)
         IF (DOIMAG) HI = HCORE(K,M,2)
C
        IF (DO2EL) THEN
         DO IELEC = 1, NAELEC
            I = IOCCI(IELEC)
            IF (I .NE. K) THEN
C---------------------
C              H = H + (I, I, K, M)
C              H = H - (I, M, K, I)
C---------------------
               HR = HR + VMUUUU(I,I,K,M,1) - VMUUUU(I,M,K,I,1)
C---------------------
               IF (DOIMAG)
     &              HI = HI + VMUUUU(I,I,K,M,2) - VMUUUU(I,M,K,I,2)
            END IF
         END DO
        END IF
         HR = IPHASE * HR
         IF (DOIMAG) HI = IPHASE * HI
         RETURN
      END IF
C
      IF(IT.EQ.4 .AND. DO2EL) THEN
C
C        ====================================
C        DETERMINANTS DIFFER BY TWO FUNCTIONS  K,L -> M,N
C        ====================================
C
         IPHASE = 0
         DO IELEC = 1,NAELEC
            IF (IOCCI(IELEC) .EQ. K) IPHASE = IPHASE + IELEC
            IF (IOCCI(IELEC) .EQ. L) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. M) IPHASE = IPHASE + IELEC
            IF (JOCCI(IELEC) .EQ. N) IPHASE = IPHASE + IELEC
         END DO
         IPHASE = (-1)**IPHASE
C
C---------------------
C        H = H + (K, M, L, N) - (K, N, L, M)
C---------------------
         HR = IPHASE * (VMUUUU(K,M,L,N,1) - VMUUUU(K,N,L,M,1))
         IF (DOIMAG)
     &   HI = IPHASE * (VMUUUU(K,M,L,N,2) - VMUUUU(K,N,L,M,2))
         RETURN
      END IF
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C
      CALL QUIT('GOSCI_HMAT called with IT .gt. 4 !')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_svc */
      SUBROUTINE GOSCI_SVC(SAVDIA,DIAG,SVEC,IDET,BVEC,HCORE,VMUUUU)
C***********************************************************************
C
C     Calculate CI sigma vectors: SVEC(I) = \sum_J H_{IJ} B_{J}
C
C     Input:
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C        BVEC   - trial vector
C        IDET   - determinant strings
C        SAVDIA - save diagonal of Hamiltonian (H_{II})
C
C     Output:
C        SVEC   - sigma vector
C        DIAG   - diagonal of Hamiltonian (if SAVDIA)
C
C     Written by J. Thyssen - Jan 3 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      LOGICAL   SAVDIA
      DIMENSION DIAG(*), SVEC(*), BVEC(*), HCORE(*), VMUUUU(*)
      DIMENSION IDET(*)
C
      CALL QENTER('GOSCI_SVC')
C
#define xxx
#ifdef xxx
C hj 16-jun-2004 new sigma code with screening
      THRSCR = 1.0D-10
      CALL GOSCI_CISIG_SCREEN(SVEC,BVEC,HCORE,VMUUUU,.TRUE.,IDET,THRSCR)
      IF (SAVDIA) CALL GOSCI_CI(SAVDIA,DIAG,
     &     .FALSE.,SVEC,BVEC,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,IDET)
#else
      CALL GOSCI_CI(SAVDIA,DIAG,
     &     .TRUE.,SVEC,BVEC,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,IDET)
#endif
C
C     CALL GOSCI_CI(
C    &     DODIA,DIAG,
C    &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C    &     DODEN,CL,CR,DV,PV,
C    &     DO2EL,IDET)
C
      CALL QEXIT('GOSCI_SVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_hdiag */
      SUBROUTINE GOSCI_HDIAG(DIAG,IDET,HCORE,VMUUUU)
C***********************************************************************
C
C     Calcualte diagonal of CI Hamiltonian.
C
C     Input:
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C        IDET   - determinant strings
C
C     Output:
C        DIAG   - diagonal of Hamiltonian (if SAVDIA)
C
C     Written by J. Thyssen - Jan 3 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DIAG(*), IDET(*), HCORE(*), VMUUUU(*)
C
      CALL QENTER('GOSCI_HDIAG')
C
      CALL GOSCI_CI(.TRUE.,DIAG,
     &     .FALSE.,DUMMY,DUMMY,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,IDET)
C
C     CALL GOSCI_CI(
C    &     DODIA,DIAG,
C    &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C    &     DODEN,CL,CR,DV,PV,
C    &     DO2EL,IDET)
C
      CALL QEXIT('GOSCI_HDIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_rcistd */
      SUBROUTINE GOSCI_RCISTD(IDET,CMO,CREF,ECI,WORK,LWORK)
C***********************************************************************
C
C     Do MAXCIT CI iterations.
C
C     Input :
C        IDET - determinant strings
C        CMO  - MO coeff.
C
C     Output:
C        CREF - CI ref.
C        ECI  - CI energy.
C
C
C     Written by J. Thyssen - Jan 3 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), CREF(NZCONF,*), CMO(*), IDET(2,*)
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
#include "memint.h"
C
      CALL QENTER('GOSCI_RCISTD')
C
      CALL GETTIM(CPU1,WLL1)
C
      IF (NCIROOT .GT. 1) THEN
         CALL QUIT('.CIROOT .gt. 1 has not been implemented yet.')
      END IF
C
C     ************************************
C     *** Calculate necesary integrals ***
C     ************************************
C
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET('REAL',KFCACM,(2*NASHT)*(2*NASHT)*MIN(NZ,2),
     $     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Calculate FCmo (TOFILE false: do not write to LUKRM3) :
C
      CALL RGETFC(.FALSE.,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
      CALL MEMREL('RCISTD after RGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,WORK,KFREE,LFREE)
C
C     Read integrals from 4IND1XXXX0:
C
      CALL RGETH21(.FALSE.,.TRUE.,DUMMY,DUMMY,WORK(KH2AC),
     &     WORK(KFREE),LFREE,IPROPT)
C
C     Make sure imag(VMUUUU) = 0 for NZ=1
C     (DNZ32M only DZERO'es real part for NZ = 1)
C
      IF (NZ .EQ. 1) CALL DZERO(WORK(KMUUUU+(2*NASHT)**4),(2*NASHT)**4)
C
C     Transform integrals to Molfdir format:
C
      CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPROPT)
      CALL MEMREL('GOSCI_RCISTD.H2AC',WORK,KWORK,KH2AC,KFREE,LFREE)
C
C
C     ********************************
C     *** Find start trial vectors ***
C     ********************************
C
C
      CALL MEMGET('REAL',KHDIAG,NZCONF,WORK,KFREE,LFREE)
C
C     Assume we can have all sigma and trial vectors in memory.
C     This could get a little expensive for large CI's.
C
      MAXR = MAX(50,MAXCIT+4)
C     MAXCIT+4: allocating space for up to 5-double degeneracy in diagonal
      CALL MEMGET('REAL',KBCVEC,MAXR * NZCONFQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCVEC,MAXR * NZCONFQ,WORK,KFREE,LFREE)
C
C     Calculate diagonal of CI matrix:
C     --------------------------------
C
c     CALL GOSCI_HDIAG(WORK(KHDIAG),IDET,WORK(KFCACM),WORK(KMUUUU))
      CALL GOSCI_CI(.TRUE.,WORK(KHDIAG),
     &     .FALSE.,DUMMY,DUMMY,WORK(KFCACM),WORK(KMUUUU),
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,IDET)
C     CALL GOSCI_CI(DODIA,DIAG,
C                   DOSIG,SVEC,BVEC,HCORE,VMUUUU,
C                   DODEN,CL,CR,DV,PV,
C                   DO2EL,IDET)
      IF (IPROPT .GE. 10) THEN
         WRITE(LUPRI,'(A)') ' (GOSCI_RCISTD) Diagonal of Hamiltonian'
         CALL OUTPUT(WORK(KHDIAG),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
      END IF
C
C     Find lowest elements in HDIAG:
C
      II = IDMIN(NZCONF,WORK(KHDIAG),1)
C
C     Set first trial vector equal to:
C     BCVEC(I) = DELTA(I,II) / HDIAG(II)
C     add trial vectors for degenerate determinants within 1.0D-5
C
      HDIAGTST = WORK(KHDIAG-1+II) + 1.0D-5
      NCSIM = 0
      DO  II = 1, NZCONF
          IF (WORK(KHDIAG-1+II) .LE. HDIAGTST) THEN
             JBCVEC = KBCVEC +  NCSIM*NZCONFQ
             CALL DZERO(WORK(JBCVEC),NZCONFQ)
             WORK(JBCVEC + (II-1)) = D1
             CALL RCINORM(WORK(JBCVEC),IPROPT)
             NCSIM = NCSIM + 1
          END IF
      END DO
C
C
C
      MZ = MIN(2,NZ)
C
      JCRED = 0
      NCRED = NCSIM
      RTHRCI = 1.0D-5
c     RTHRCI = 1.0D-12
c     RTHRCI = 1.0D-6
      ECI = D0
C     GOSCI screening, a calculation in C2 symmetry showd that
C     THRSCR = RTHRCI**2 was not good enough for the complex group /hjaaj July 04
      THRSCR = RTHRCI**2 * 0.01D0
C
      LHRED = MAXR*MAXR*MZ
      CALL MEMGET('REAL',KHRED,LHRED,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVL,MAXR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,LHRED,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHRED),LHRED)
      KHREDI = KHRED + MAXR*MAXR
      KEIGVCI = KEIGVC + MAXR*MAXR
C
      ITMIC = 0
      ICONV = 0
      CALL GETTIM(CPU2,WLL2)
 100  CONTINUE
C
C        **********************************
C        *** Start micro iteration loop ***
C        **********************************
C
         ITMIC = ITMIC + 1
         CALL GETTIM(CPU3,WLL3)
C
         WRITE(LUPRI,'(//A,I4)')
     &        ' (GOSCI_RCISTD) CI microiteration no.',ITMIC
         CALL FLSHFO(LUPRI)
C
C        Calculate sigma vector(s)
C        -------------------------
C
         DO I = 1, NCSIM
c           CALL GOSCI_SVC(.FALSE.,DUMMY,
c    &           WORK(KSCVEC+(I-1+JCRED)*NZCONFQ),
c    &           IDET,WORK(KBCVEC+(I-1+JCRED)*NZCONFQ),
c    &           WORK(KFCACM),WORK(KMUUUU))
#define xxx
#ifdef xxx
C hj 16-jun-2004 new sigma code with screening
            CALL GOSCI_CISIG_SCREEN(
     &           WORK(KSCVEC+(I-1+JCRED)*NZCONFQ),
     &           WORK(KBCVEC+(I-1+JCRED)*NZCONFQ),
     &           WORK(KFCACM),WORK(KMUUUU),
     &           .TRUE.,IDET,THRSCR)
#else
            CALL GOSCI_CI(.FALSE.,DUMMY,
     &           .TRUE.,WORK(KSCVEC+(I-1+JCRED)*NZCONFQ),
     &           WORK(KBCVEC+(I-1+JCRED)*NZCONFQ),
     &           WORK(KFCACM),WORK(KMUUUU),
     &           .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &           .TRUE.,IDET)
#endif
C           CALL GOSCI_CI(DODIA,DIAG,
C                         DOSIG,SVEC,BVEC,
C                         HCORE,VMUUUU,
C                         DODEN,CL,CR,DV,PV,
C                         DO2EL,IDET)
C
            IF (IPROPT .GE. 10) THEN
C
C              Trial vector
C
               WRITE(LUPRI,'(/A,I3)')
     &              ' (GOSCI_RCISTD) Trial vector no. ',I
               CALL OUTPUT(WORK(KBCVEC+(I-1+JCRED)*NZCONFQ),
     &              1,NZCONF,1,MZ,NZCONF,MZ,1,LUPRI)
C
C              Sigma vector
C
               WRITE(LUPRI,'(/A,I3)')
     &              ' (GOSCI_RCISTD) Sigma vector no. ',I
               CALL OUTPUT(WORK(KSCVEC+(I-1+JCRED)*NZCONFQ),
     &              1,NZCONF,1,MZ,NZCONF,MZ,1,LUPRI)
            END IF
         END DO
C
C        Extend reduced Hamiltonian:
C        ---------------------------
C
C        HRED(I,J) = BCVEC(I)^+ * SCVEC(J)
C
C        New trial vectors * Old sigma vectors
C
         DO J = 1, JCRED
            DO I = JCRED+1, NCRED
               WORK(KHRED+(I-1)*MAXR+(J-1)) =
     &              DDOT(NZCONFQ,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ),1)
               IF (MZ .GT. 1) THEN
                  WORK(KHREDI+(I-1)*MAXR+(J-1)) =
     &              DDOT(NZCONF,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ+NZCONF),1)
     &            - DDOT(NZCONF,WORK(KBCVEC+(I-1)*NZCONFQ+NZCONF),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ),1)
               END IF
            END DO
         END DO
C
C        Old trial vectors * New sigma vectors
C
         DO J = JCRED+1, NCRED
            DO I = 1, JCRED
               WORK(KHRED+(I-1)*MAXR+(J-1)) =
     &              DDOT(NZCONFQ,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ),1)
               IF (MZ .GT. 1) THEN
                  WORK(KHREDI+(I-1)*MAXR+(J-1)) =
     &              DDOT(NZCONF,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ+NZCONF),1)
     &            - DDOT(NZCONF,WORK(KBCVEC+(I-1)*NZCONFQ+NZCONF),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ),1)
               END IF
            END DO
         END DO
C
C        New trial vectors * New sigma vectors
C
         DO J = JCRED+1, NCRED
            DO I = JCRED+1, NCRED
               WORK(KHRED+(I-1)*MAXR+(J-1)) =
     &              DDOT(NZCONFQ,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ),1)
               IF (MZ .GT. 1) THEN
                  WORK(KHREDI+(I-1)*MAXR+(J-1)) =
     &              DDOT(NZCONF,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ+NZCONF),1)
     &            - DDOT(NZCONF,WORK(KBCVEC+(I-1)*NZCONFQ+NZCONF),1,
     &              WORK(KSCVEC+(J-1)*NZCONFQ),1)
               END IF
            END DO
         END DO
C
C        Check for inconsistencies:
C
chj      1.0D-10 can give problems with 1.0D-10 screening in GOSCI
chj      THRDAS = 1.0D-10
         THRDAS = THRSCR*10.0D0
         DAS = D0
         DO I = 1, NCRED
            DO J = 1, I - 1
               DAS = DAS +
     &              ABS( WORK(KHRED + (J-1) + (I-1)*MAXR) -
     &              WORK(KHRED + (I-1) + (J-1)*MAXR) )
               IF (MZ .GT. 1) DAS = DAS +
     &              ABS( WORK(KHREDI + (J-1) + (I-1)*MAXR) +
     &              WORK(KHREDI + (I-1) + (J-1)*MAXR) )
            END DO
         END DO
         IF (IPROPT .GE. 6 .OR. DAS .GT. THRDAS) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GOSCI_RCISTD) Reduced Hamiltonian, real part'
            CALL OUTPUT(WORK(KHRED),1,NCRED,1,NCRED,MAXR,MAXR,1,LUPRI)
           IF (MZ .GT. 1) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GOSCI_RCISTD) Reduced Hamiltonian, imaginary part'
            CALL OUTPUT(WORK(KHREDI),1,NCRED,1,NCRED,MAXR,MAXR,1,LUPRI)
           END IF
         END IF
         IF (IPROPT .GE. 5 .OR. DAS .GT. 0.01D0*THRDAS) THEN
            WRITE(LUPRI,'(A,1P,D10.3)')
     &      ' (GOSCI_RCISTD) Anti-hermiticity of red. Hamiltonian = ',
     &      DAS
            IF (DAS .GT. THRDAS) THEN
               call QUIT('red. Hamiltonian not hermitian')
            END IF
         END IF
C
C        Diagonalize reduced Hamiltonian
C        -------------------------------
C
         IF (MZ.EQ.1) THEN
            IJOB = 1
            IORDER = 1
            IPACK = 0
            CALL RSJACO(MAXR,NCRED,NCRED,WORK(KHRED),WORK(KEIGVL),
     &        IJOB,IORDER,IPACK,WORK(KEIGVC),WORK(KFREE),LFREE)
         ELSE
            MATZ = 1
            CALL MEMGET('REAL',KHREDCP,LHRED,WORK,KFREE,LFREE)
            CALL DCOPY(LHRED,WORK(KHRED),1,WORK(KHREDCP),1)
            CALL QDIAG(2,NCRED,WORK(KHREDCP),MAXR,MAXR,
     &        WORK(KEIGVL),MATZ,WORK(KEIGVC),MAXR,MAXR,
     &        WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in GOSCI_RCISTD ***: ',
     &           'QDIAG failed with error code ',IERR
               CALL QUIT('*** ERROR in GOSCI_RCISTD *** QDIAG failed')
            END IF
            CALL MEMREL('GOSCI_RCISTD.QDIAG',
     &         WORK,KWORK,KHREDCP,KFREE,LFREE)
         END IF
C
         IF (IPROPT .GE. 5) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (GOSCI_RCISTD) Eigenvalues of reduced Hamiltonian'
            CALL OUTPUT(WORK(KEIGVL),1,1,1,NCRED,1,MAXR,1,LUPRI)
         END IF
         IF (IPROPT .GE. 6) THEN
            WRITE(LUPRI,'(/A)') ' (GOSCI_RCISTD)'//
     &      ' Eigenvectors of reduced Hamiltonian, real part'
            CALL OUTPUT(WORK(KEIGVC),1,NCRED,1,NCRED,MAXR,MAXR,1,LUPRI)
           IF (MZ .GT. 1) THEN
            WRITE(LUPRI,'(/A)') ' (GOSCI_RCISTD)'//
     &      ' Eigenvectors of reduced Hamiltonian, imaginary part'
            CALL OUTPUT(WORK(KEIGVCI),1,NCRED,1,NCRED,MAXR,MAXR,1,LUPRI)
           END IF
         END IF
         ECIOLD = ECI
         ECIA = WORK(KEIGVL + IOPT_ISTATE - 1)
         ECI = ECORE + ECIA
         IF (IPROPT .GE. 4) THEN
            WRITE(LUPRI,'(2(/A,F20.10))')
     &        '                Core energy        = ',ECORE,
     &        '                Active energy      = ',ECIA
         END IF
         WRITE(LUPRI,'(A,F20.10)')
     &        ' (GOSCI_RCISTD) Total energy       = ',ECI,
     &        '                Lowering of energy = ',(ECIOLD-ECI)
         CALL FLSHFO(LUPRI)
C
C        Construct residual
C        ------------------
C
C        r = (H-E)x  (store in CREF)
C          = Hx - Ex
C          = SUM_J C(J) SCVEC(J) - ECIA * SUM_J C(J) BCVEC(J)
C
         CALL DZERO(CREF,NZCONFQ)
         DO I = 1, NCRED
            FAC = WORK(KEIGVC+(I-1)+(IOPT_ISTATE-1)*MAXR)
            CALL DAXPY(NZCONFQ,FAC,WORK(KSCVEC+(I-1)*NZCONFQ),1,
     &           CREF,1)
            CALL DAXPY(NZCONFQ,-ECIA*FAC,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &           CREF,1)
         END DO
         DNRES = DNORM2(NZCONFQ,CREF,1)
         WRITE(LUPRI,'(A,1P,D18.3,A,D8.2,A)')
     &        '                Norm of residual   =',DNRES,
     &        '  (thr = ',RTHRCI,')'
         IF (IPROPT .GE. 10) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GOSCI_RCISTD) Residual'
            CALL OUTPUT(CREF,1,NZCONF,1,MZ,NZCONF,MZ,1,LUPRI)
            WRITE(LUPRI,'()')
         END IF
C
         CALL FLSHFO(LUPRI)
C
         IF (DNRES .LE. RTHRCI) THEN
            ICONV = 1
            GOTO 200
         END IF
         IF (ITMIC .GE. MAXCIT) THEN
            ICONV = -2
            GOTO 200
         END IF
         IF (NCRED .GE. MAXR) THEN
            ICONV = -3
            GOTO 200
         END IF
C
C        Construct 1 new trial vector
C        ----------------------------
C
         JCRED = NCRED
         NCRED = NCRED + 1
         NCSIM = 1
C
         DO IZ = 1, MZ
            DO I = 1, NZCONF
               HD = WORK(KHDIAG+I-1)
               IF (ABS(HD-ECIA) .GT. THRZER) THEN
                  WORK(KBCVEC+(I-1)+(NCRED-1)*NZCONFQ + (IZ-1)*NZCONF) =
     &                 - CREF(I,IZ)/(HD - ECIA)
               ELSE
                  WORK(KBCVEC+(I-1)+(NCRED-1)*NZCONFQ + (IZ-1)*NZCONF) =
     &                 - CREF(I,IZ)
               END IF
            END DO
         END DO
         CALL RCINORM(WORK(KBCVEC+(NCRED-1)*NZCONFQ),IPROPT)
C
C        Orthogonalize trial vector(s)
C        -----------------------------
C
         IROUND = 0
         ITURN  = 0
 400     CONTINUE
C
            ITURN = ITURN + 1
C
            IF (NCSIM .GT. 1)
     &           CALL QUIT('*** ERROR in GOSCI_RCISTD ***'//
     &           ' NCSIM.gt.1 not implemented')
C
C
C           Project out all previous trial vectors
C
            DNBEF = DNORM2(NZCONFQ,WORK(KBCVEC+(NCRED-1)*NZCONFQ),1)
            DO I = 1, JCRED
C
C              FIXME: loop over new trial vectors
C
C              Calculate <B(new)|B(old)> / <B(old)|B(old)>
C              = <B(new)|B(old)>  (as trial vectors are normed to 1)
C
               OVLP = DDOT(NZCONFQ,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &              WORK(KBCVEC+(NCRED-1)*NZCONFQ),1)
C
C              Calculate B(new) = B(new) - OVLP * B(old)
C
               IF (ABS(OVLP) .GT. THRZER) THEN
                  CALL DAXPY(NZCONFQ,-OVLP,WORK(KBCVEC+(I-1)*NZCONFQ),1,
     &                 WORK(KBCVEC+(NCRED-1)*NZCONFQ),1)
               END IF
            END DO
C
C           FIXME: Orthogonalize new trial vectors against each other
C
            DNAFT = DNORM2(NZCONFQ,WORK(KBCVEC+(NCRED-1)*NZCONFQ),1)
            IF (ABS(DNAFT) .LE. THRZER) THEN
               ICONV = -1
               GO TO 200
            END IF
            IF (DNAFT .LE. 1.0D-4 * DNBEF .AND. ITURN .EQ. 1)
     &           IROUND = IROUND + 1
C
            CALL DSCAL(NZCONFQ,1/DNAFT,WORK(KBCVEC+(NCRED-1)*NZCONFQ),1)
C hj 16-jun-2004 TODO: screen on elements here ...
C
            IF (IROUND .GT. 0) THEN
               IROUND = 0
               IF (ITURN .EQ. 1) GO TO 400
               CALL QUIT('*** ERROR in GOSCI_RCISTD ***' //
     &              'IROUND .gt. 0 and ITURN .gt. 1')
            END IF
C
C           Print timings
C
            IF (IPROPT .GE. 2) THEN
               CALL GETTIM(CPU4,WLL4)
               CPUTID = SECTID(CPU4-CPU3)
               WALLTID = SECTID(WLL4-WLL3)
               WRITE(LUPRI,'(/A,I3,5A)')
     &              ' (GOSCI_RCISTD) CPU (Wall) time for iteration ',
     &              ITMIC,': ',CPUTID,' (',WALLTID,')'
            END IF
C
         GOTO 100
C
 200  CONTINUE
C
C     Print timings
C
      CALL GETTIM(CPU4,WLL4)
      IF (IPROPT .GE. 2) THEN
         CPUTID = SECTID(CPU4-CPU3)
         WALLTID = SECTID(WLL4-WLL3)
         WRITE(LUPRI,'(/A,I4,5A)')
     &        ' (GOSCI_RCISTD) CPU (Wall) time for iteration',
     &        ITMIC,': ',CPUTID,' (',WALLTID,')'
      END IF
      IF (IPROPT .GE. 0) THEN
         CPUTID = SECTID(CPU4-CPU2)
         WALLTID = SECTID(WLL4-WLL2)
         WRITE(LUPRI,'(/6A)')
     &        ' (GOSCI_RCISTD) Total CPU (Wall) '//
     &        'time for microiterations: ',
     &        CPUTID,' (',WALLTID,')'
      END IF
C
      IF (ICONV .GT. 0) THEN
         WRITE(LUPRI,'(//A)')
     &        ' (GOSCI_RCISTD) Micro iterations converged.'
      ELSE
         IF (ICONV .EQ. -1) THEN
            WRITE(LUPRI,'(//A)') ' (GOSCI_RCISTD) '//
     &           'WARNING: linear dependency in trial vectors.'
         ELSE IF (ICONV .EQ. -2) THEN
            WRITE(LUPRI,'(/A,I4,A)') ' (GOSCI_RCISTD) '//
     &           'WARNING: maximum number of micro iterations',
     &           MAXCIT, ' reached.'
         ELSE IF (ICONV .EQ. -3) THEN
            WRITE(LUPRI,'(//A,I4,A)') ' (GOSCI_RCISTD) '//
     &           'WARNING: maximum reduced space dimension of',
     &           MAXR, ' reached.'
         ELSE
            WRITE(LUPRI,'(//A,I4/A)') ' (GOSCI_RCISTD) '//
     &        'WARNING: micro iterations aborted while not converged.',
     &        ' ICONV value :',ICONV
         END IF
      END IF
C
C
      MAXCIT = ITMIC
      ECI = ECORE + ECIA
      WRITE(LUPRI,'(3(/A,F20.10))')
     &     '                Core energy       = ',ECORE,
     &     '                Active energy     = ',ECIA,
     &     ' (GOSCI_RCISTD) Final CI energy   = ',ECI
C
C     Construct CREF:
C
      CALL DZERO(CREF,NZCONFQ)
      DO I = 1, NCRED
         CALL DAXPY(NZCONFQ,WORK(KEIGVC+(I-1)+(IOPT_ISTATE-1)*MAXR),
     &        WORK(KBCVEC+(I-1)*NZCONFQ),1,CREF,1)
      END DO
C
C     Check norm of CI vector
C
      CALL RCINORM(CREF,10)
C
      IF ( IPROPT .GE. 5 ) THEN
         WRITE(LUPRI,'(/A)') ' (GOSCI_RCISTD) CI start vector'
         CALL RPRCI(CREF,NZCONF,NZ,OPT_CIPROGRAM, THRPCI,LUPRI)
      END IF
C
      CALL GETTIM(CPU2,WLL2)
      IF (IPROPT .GE. 0) THEN
         CPUTID = SECTID(CPU2-CPU1)
         WALLTID = SECTID(WLL2-WLL1)
         WRITE(LUPRI,9000) CPUTID,WALLTID
      END IF
 9000 FORMAT(/' (GOSCI_RCISTD) Total CPU (Wall) time for routine: ',
     &     A,' (',A,')')
C
      CALL MEMREL('GOSCI_RCISTD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GOSCI_RCISTD')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck idet2str */
      SUBROUTINE IDET2STR(IA,IB,IASTR,IBSTR,IAELEC,IBELEC)
C***********************************************************************
C
C     Transform bit packed IA and IB into alpha- and beta string.
C
C     Input :
C        JDET - determinant
C
C     Output:
C        IASTR, IBSTR - strings
C        IAELEC, IBELEC - alpha and beta electrons.
C
C     Written by J. Thyssen - Jan 5 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION IASTR(*), IBSTR(*)
C
      LOGICAL LBIT
C
#include "ibtfun.h"
C
      IAELEC = 0
      IBELEC = 0
C
      DO J = 1, NASHT
         IF (LBIT(IA,J)) THEN
            IAELEC = IAELEC + 1
            IASTR(IAELEC) = IPTA2O(J)
         END IF
         IF (LBIT(IB,J)) THEN
            IBELEC = IBELEC + 1
            IBSTR(IBELEC) = IPTB2O(J)
         END IF
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_anaci */
      SUBROUTINE GOSCI_ANACI(CREF,IDET,WORK,LWORK)
C***********************************************************************
C
C     Print nice analysis of CI vector.
C
C     Input :
C        CREF - CI ref. vector
C        IDET - determintants
C
C     Output:
C
C     Written by J. Thyssen - Jan 5 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
      PARAMETER (D10 = 10.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcbgosci.h"
C
      DIMENSION CREF(NZCONF,*), IDET(2,*), WORK(*)
C
      DIMENSION IASTR(2*MAXASH_GOSCI), IBSTR(2*MAXASH_GOSCI)
C
      CALL QENTER('GOSCI_ANACI')
C
      MZ = MIN(2,NZ)
C
      DMX = D1
      DMN = D1 / D10
 10   CONTINUE
         WRITE(LUPRI,'(/A,1P,D7.1,A,D7.1/)')
     &     ' Print of CI coefficients in range ',DMN,' to ',DMX
         K = 0
         DO I = 1, NZCONF
            DI = D0
            DO IZ = 1, MZ
               DI = DI + CREF(I,IZ) * CREF(I,IZ)
            END DO
            IF (DI .GT. DMN .AND. DI .LE. DMX) THEN
               K = K + 1
               IF (MZ .EQ. 1) THEN
                  WRITE(LUPRI,'(3X,A,I8,A,1P,D20.10,A,0P,F10.8)')
     &                 'Determinant no. ',I,' coeff. ',
     &                 CREF(I,1),'; weight ',DI
               ELSE
                  WRITE(LUPRI,
     &                 '(3X,A,I8,A,1P,D20.10,A,D20.10,A,0P,F10.8)')
     &                 'Determinant no. ',I,' coeff. ',
     &                 CREF(I,1),' + ',CREF(I,2),' i; weight ',DI
               END IF
               CALL IDET2STR(IDET(1,I),IDET(2,I),IASTR,IBSTR,IA,IB)
               WRITE(LUPRI,'(3X,A,40I4)')
     &              'P-string     : ',(IASTR(J),J=1,IA)
               WRITE(LUPRI,'(3X,A,40I4/)')
     &              'P(bar)-string: ',(IBSTR(J),J=1,IB)
            END IF
         END DO
         IF (K .EQ. 0) WRITE(LUPRI,'(3X,A/)') 'None'
         DMX = DMN
         DMN = DMN / D10
         IF (DMX .GT. THRPCI) GO TO 10
C
      CALL QEXIT('GOSCI_ANACI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_chkkr */
      SUBROUTINE GOSCI_CHKKR(DOIMAG,BVEC,WORK,LWORK)
C***********************************************************************
C
C     Check that BVEC is Kramers restricted.
C
C     Input :
C       BVEC: CI vector
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      LOGICAL   DOIMAG
      DIMENSION BVEC(*), WORK(*)
C
#include "memint.h"
C
      IF ( MOD(NAELEC,2) .NE. 0) RETURN
C
      CALL MEMGET('INTE',KRLIST,NZCONF,WORK,KFREE,LFREE)
      CALL IREAKRMC(LUKRMC,'KRLIST  ',WORK(KRLIST),NZCONF)
      CALL GOSCI_CHKKR1(DOIMAG,BVEC,WORK(KRLIST))
      CALL MEMREL('GOSCI_CHKKR',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_chkkr1 */
      SUBROUTINE GOSCI_CHKKR1(DOIMAG,BVEC,KRLIST)
C***********************************************************************
C
C     Check that BVEC is Kramers restricted.
C
C     Input :
C       BVEC: CI vector
C       KRLIST: list of determinants and their Kramers partner.
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dgroup.h"
C
      LOGICAL   DOIMAG
      DIMENSION BVEC(NZCONF,*), KRLIST(*)
C
      DO I = 1, NZCONF
         IK = KRLIST(I)
         IF (IK .NE. I) THEN
C HJAaJ TODO aug 2002: we may divide by zero in next line!
            IF ( ABS(BVEC(I,1) - BVEC(IK,1))/ABS(BVEC(I,1))
     &           .GT. 1.0D-13) THEN
              write(LUPRI,'(A,2I4,3D25.10)') ' GOSCI_CHKKR1: real diff',
     &              I,IK,BVEC(I,1),BVEC(IK,1),
     &              ABS(BVEC(I,1) - BVEC(IK,1))
            END IF
C HJAaJ TODO aug 2002: we may divide by zero in next line!
            IF (DOIMAG) THEN
            IF ( ABS(BVEC(I,2) + BVEC(IK,2))/ABS(BVEC(I,2))
     &           .GT. 1.0D-13) THEN
              write(LUPRI,'(A,2I4,3D25.10)') ' GOSCI_CHKKR1: imag diff',
     &              I,IK,BVEC(I,2),BVEC(IK,2),
     &              ABS(BVEC(I,2) + BVEC(IK,2))
            END IF
            END IF
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_genstr */
      SUBROUTINE GOSCI_GENSTR(IEL,ISTR,NSTR,MSKGAS,MAXDET)
C***********************************************************************
C
C     Generate all strings with IEL electrons and GAS constraints
C     fullfilled.
C
C     Input :
C       IEL: number of electrons
C       MSKGAS: bit mask for GAS spaces
C
C     Output:
C       ISTR: strings
C       NSTR: #strings
C
C     Written by J. Thyssen - Jan 23 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgosci.h"
C
      DIMENSION ISTR(*), MSKGAS(*)
C
      DIMENSION NGSPC(2*MXMCGAS)
C
C     Calculate minimum and maximum number of
C     electrons in gas spaces.
C
      J = 1
      NGTOT = 0
      DO I = 1, NGAS
C
         NGTOT = NGTOT + NGSHT(I)
         NGSPC(J) = MAX(0,NGSP(1,I) - NGTOT)
         NGSPC(J+1) = MIN(NGSP(2,I),NGTOT)
C        write(lupri,*) 'gosci_genstr: ngspc: ',ngspc(J),ngspc(j+1)
         J = J + 2
      END DO
C
      CALL CGENDET(IEL,NASHT,NSTR,ISTR,NGAS,MSKGAS,NGSPC,MAXDET)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_ci */
      SUBROUTINE GOSCI_CI(
     &     DODIA,DIAG,
     &     DOSIG,SVEC,BVEC,HCORE,VMUUUU,
     &     DODEN,CL,CR,DV,PV,
     &     DO2EL,IDET)
C***********************************************************************
C
C     Calculate:
C     (1) diagonal of CI Hamiltonian (DODIA)
C     (2) CI sigma vector (DOSIG)
C     (3) Density matrices (DODEN)
C
C     Input:
C        DODIA, DOSIG, DODEN - control flags
C
C        BVEC - trial vector
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C
C        CL, CR - left and right vector for density matrices
C
C        DO2EL - include 2-electron part
C                -- DODIA and DOSIG: include VMUUUU
C                -- DODEN: calculate PV
C        IDET  - determinant strings
C
C     Output:
C        DIAG   - diagonal of Hamiltonian (if SAVDIA)
C
C        SVEC   - sigma vector
C
C        DV, PV - density matrices
C
C     Written by J. Thyssen and H.J.Aa. Jensen - Jan-Mar 2001
C     Last revision : July 2002/Sep 2003 hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "thrzer.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
C
      DIMENSION DIAG(*)
C
      DIMENSION SVEC(NZCONF,*), BVEC(NZCONF,*)
      DIMENSION VMUUUU(*), HCORE(*)
C
      DIMENSION CL(NZCONF,*), CR(NZCONF,*)
      DIMENSION DV(2*NASHT,2*NASHT,*)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      LOGICAL   DODIA, DOSIG, DODEN, DO2EL, DODENCHK
C
      DIMENSION IDET(2,*)
C
      LOGICAL   DOIMAG, DOIMAG1
C
      DIMENSION IOCCI(2*MAXASH_GOSCI), JOCCI(2*MAXASH_GOSCI)
      DIMENSION IBTAB(0:255)
      INTEGER*8 FNCALLS, CNCALLS
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
#include "ibtfun.h"
C
      CALL GETTIM(CPU1,WLL1)
C
C     Initialization
C     --------------
C
C     Bit table:
C
      CALL IBTABINI(IBTAB)
C
C     Zero counters
C
      CALL GETNBITSCALLS(1,FNCALLS,CNCALLS)
C
C     Zero vectors:
C
      MZ = MIN(NZ,2)
      IF (DODIA) CALL DZERO(DIAG,NZCONF)
      IF (DOSIG) CALL DZERO(SVEC,NZCONFQ)
      IF (DODEN) THEN
         CALL DZERO(DV,(2*NASHT)**2 * MZ)
         IF (DO2EL) CALL DZERO(PV,(2*NASHT)**4 * MZ)
      END IF
C
      IF (DO2EL) THEN
C        max 4 differences between 2 determinants
C        for 2-electron terms
         ITMAX = 4
      ELSE
C        max 2 differences between 2 determinants
C        when only 1-electron terms
         ITMAX = 2
      END IF
C
      DOIMAG  = NZ.GT.1
      DOIMAG1 = DOIMAG
C
C     screening for density matrices:
C
      DODENCHK = DODEN .AND. .NOT.DOSIG .AND. .NOT.DODIA
      THRSCR = 1.0D-10
      NDSKIP = 0
C
      MK2STP = 2
C
      DO 10 MK2I = MINMK2, MAXMK2, MK2STP
C        ... MK2 is NAEL - NBEL
         IF (NMK2(MK2I) .EQ. 0) GOTO 10
D        write(lupri,*) '(GOSCI_CI) MK2I = ',MK2I
         ISTRT = IMK2(MK2I) + 1
         IEND  = IMK2(MK2I) + NMK2(MK2I)
C
Chj 16June2004: I don't understand how the -MAXMK2 in next line can be right ...
Chj      MK2JSTRT = MAX(MINMK2,-MAXMK2,MK2I-ITMAX)
         MK2JSTRT = MAX(MINMK2,MK2I-ITMAX)
         MK2JEND  = MK2I - MK2STP
C
C
C           Loop over determinants
C           (use that outer loop in gosci_gendet is over A-strings,
C           inner loop is over B-strings)
C
            IDAOLD = -1
            DO 100 I = ISTRT, IEND
               IF (DODENCHK) THEN
                  DENCHK = ABS(CL(I,1)) + ABS(CR(I,1))
                  IF (DOIMAG)
     &            DENCHK = DENCHK + ABS(CL(I,2)) + ABS(CR(I,2))
                  IF (DENCHK .LT. THRSCR) THEN
                     NDSKIP = NDSKIP + 1
                     GO TO 100
                  END IF
               ENDIF
C
C
               IDA = IDET(1,I)
               IF (IDA .NE. IDAOLD) THEN
C                 ... outer loop in gosci_gendet is over A-strings
                  IAELEC = 0
                  DO K = 1, NASHT
                     IF (IBTAND (IDA, IBTSHL (1, K - 1)) .NE. 0) THEN
                        IAELEC = IAELEC + 1
                        IOCCI(IAELEC) = IPTA2O(K)
                     END IF
                  END DO
                  IDAOLD = IDA
               END IF
C
               IELEC = IAELEC
               IDB = IDET(2,I)
               DO K = 1, NASHT
                  IF (IBTAND (IDB, IBTSHL (1, K - 1)) .NE. 0) THEN
                     IELEC = IELEC + 1
                     IOCCI(IELEC) = IPTB2O(K)
                  END IF
               END DO
               IF (IELEC .NE. NAELEC)
     &            CALL QUIT('GOSCI_CI: IELEC.ne.NAELEC')
C
C        *************************
C        *** Off-diagonal loop *** (MK2J .lt. MK2I)
C        *************************
C
C
            IF (.NOT. (DODEN .OR. DOSIG)) GOTO 21
C           ...only off-diagonal blocks for dens.mat or sigma-vec.
C
            DO 20 MK2J = MK2JSTRT, MK2JEND, MK2STP
               IF (NMK2(MK2J) .EQ. 0) GOTO 20
D              write(lupri,*) '(GOSCI_CI) MK2I MK2J = ',MK2I,MK2J
D              call flshfo(lupri)
C
               JSTRT = IMK2(MK2J) + 1
               JEND  = IMK2(MK2J) + NMK2(MK2J)
               MXNDAV = ITMAX + (MK2J - MK2I)/2
C              ... we know that NDBV .ge. (MK2I-MK2J)/2
C
               JDAOLD = -1
               DO 200 J = JSTRT, JEND
C
                  JDA = IDET(1,J)
                  IF (JDA .NE. JDAOLD) THEN
C                    ... outer loop in gosci_gendet is over A-strings
                     IDAV = IBTXOR (IDA, JDA)
                     CALL FNBITS(NDAV,IDAV,IBTAB)
                     IF (NDAV .GT. MXNDAV) GOTO 200
C
                     JAELEC = 0
                     IAI = IBTAND (IDAV, IDA)
                     KA  = 0
                     IAJ = IBTAND (IDAV, JDA)
                     MA  = 0
                     DO II = 1, NASHT
                        IIBIT = IBTSHL (1, II - 1)
                        IF (IBTAND (JDA, IIBIT) .NE. 0) THEN
                           JAELEC = JAELEC + 1
                           JOCCI(JAELEC) = IPTA2O(II)
                        END IF
                        IF (IBTAND (IAI, IIBIT) .NE. 0) THEN
                           L  = KA
                           KA = IPTA2O(II)
                        END IF
                        IF (IBTAND (IAJ, IIBIT) .NE. 0) THEN
                           N  = MA
                           MA = IPTA2O(II)
                        END IF
                     END DO
                     JDAOLD = JDA
                  END IF
                  IF (NDAV .GT. MXNDAV) GOTO 200
C                 ... see note above about MXNDAV
C
                  JDB  = IDET(2,J)
                  IDBV = IBTXOR (IDB, JDB)
                  CALL FNBITS(NDBV,IDBV,IBTAB)
C
                  IT = NDAV + NDBV
                  IF (IT .GT. ITMAX) GOTO 200
C                 ... IT = number of different spinors
C                          in det(I) and det(J)
                  IF (DODENCHK) THEN
                     DENCHK = ABS(CL(J,1)) + ABS(CR(J,1))
                     IF (DOIMAG)
     &               DENCHK = DENCHK + ABS(CL(J,2)) + ABS(CR(J,2))
                     IF (DENCHK .LT. THRSCR) GO TO 200
                  ENDIF
C
                  IELEC = JAELEC
                  IBI = IBTAND (IDBV, IDB)
                  K = KA
                  IBJ = IBTAND (IDBV, JDB)
                  M = MA
                  DO II = 1, NASHT
                     IIBIT = IBTSHL (1, II - 1)
                     IF (IBTAND (JDB, IIBIT) .NE. 0) THEN
                        IELEC = IELEC + 1
                        JOCCI(IELEC) = IPTB2O(II)
                     END IF
C
                     IF (IBTAND (IBI, IIBIT) .NE. 0) THEN
                        L = K
                        K = IPTB2O(II)
                     END IF
                     IF (IBTAND (IBJ, IIBIT) .NE. 0) THEN
                        N = M
                        M = IPTB2O(II)
                     END IF
                  END DO
C
C                 Hamiltonian matrix element
C                 --------------------------
C
                  IF (DOSIG) THEN
C
                     CALL GOSCI_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,
     &                    K,L,M,N,HR,HI,HCORE,VMUUUU)
C
C                    S_I += H_{IJ} B_{J}
C
                     SVEC(I,1) = SVEC(I,1) + HR * BVEC(J,1)
                     IF (DOIMAG) THEN
                        SVEC(I,1) = SVEC(I,1) - HI * BVEC(J,2)
                        SVEC(I,2) = SVEC(I,2)
     &                       + HR * BVEC(J,2) + HI * BVEC(J,1)
                     END IF
C
C                    S_J += H_{IJ}^* B_{I}
C
                     SVEC(J,1) = SVEC(J,1) + HR * BVEC(I,1)
                     IF (DOIMAG) THEN
                        SVEC(J,1) = SVEC(J,1) + HI * BVEC(I,2)
                        SVEC(J,2) = SVEC(J,2)
     &                       + HR * BVEC(I,2) - HI * BVEC(I,1)
                     END IF
                  END IF
C
C                 Density matrix elements
C                 -----------------------
C
                  IF (DODEN) THEN
C
                     DMRIJ = CL(I,1)*CR(J,1)
                     DMRJI = CL(J,1)*CR(I,1)
                     IF (DOIMAG) THEN
                        DMRIJ = DMRIJ + CL(I,2)*CR(J,2)
                        DMIIJ = CL(I,1)*CR(J,2) - CL(I,2)*CR(J,1)
                        DMRJI = DMRJI + CL(J,2)*CR(I,2)
                        DMIJI = CL(J,1)*CR(I,2) - CL(J,2)*CR(I,1)
                        DOIMAG1 = (ABS(DMIIJ)+ABS(DMIJI)) .GT. THRZER
                     END IF
                     CALL GOSCI_PMAT(DO2EL,DOIMAG1,IT,IOCCI,JOCCI,
     &                    K,L,M,N,
     &                    DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C
                  END IF
 200           CONTINUE
 20         CONTINUE
 21         CONTINUE
C
C
C
C        ***********************
C        *** Diagonal blocks *** (MK2J .eq. MK2I)
C        ***********************
C
C           Diagonal element (Det_I .eq. DET_J)
C
C           Calculate Hamiltonian matrix element
C           ------------------------------------
C
            IF (DODIA .OR. DOSIG) THEN
C
               CALL GOSCI_HMAT(DO2EL,.FALSE.,0,IOCCI,IOCCI,0,0,0,0,
     &              HR,DUMMY,HCORE,VMUUUU)
C
C              Sigma-vector elements
C
C              SVEC_I += H_{II} B_{I}
C
               IF (DOSIG) THEN
C
                  SVEC(I,1) = SVEC(I,1) + HR * BVEC(I,1)
                  IF (DOIMAG) SVEC(I,2) = SVEC(I,2) + HR * BVEC(I,2)
C
               END IF
C
               IF (DODIA) DIAG(I) = HR
C
            END IF
C
C           Calculate density matrix elements
C           ---------------------------------
C
            IF (DODEN) THEN
C
               IF (.NOT. DOIMAG) THEN
                  DMR = CL(I,1)*CR(I,1)
               ELSE
                  DMR = CL(I,1)*CR(I,1) + CL(I,2)*CR(I,2)
                  DMI = CL(I,1)*CR(I,2) - CL(I,2)*CR(I,1)
                  DOIMAG1 = ABS(DMI) .GT. THRZER
               END IF
               CALL GOSCI_PMAT(DO2EL,DOIMAG1,0,IOCCI,IOCCI,0,0,0,0,
     &              DMR,DMI,DMR,DMI,DV,PV)
C
            END IF
C
C           Off-diagonal elements
C           ---------------------
C
            IF (.NOT. (DOSIG .OR. DODEN)) GOTO 100
            JDAOLD = -1
            DO 400 J = ISTRT, I - 1
C
               JDA = IDET(1,J)
               IF (JDA .NE. JDAOLD) THEN
C                 ... outer loop in gosci_gendet is over A-strings
                  IDAV = IBTXOR (IDA, JDA)
                  CALL FNBITS(NDAV,IDAV,IBTAB)
                  IF (NDAV .GT. ITMAX) GOTO 400
C
                  JAELEC = 0
                  IAI = IBTAND (IDAV, IDA)
                  KA  = 0
                  IAJ = IBTAND (IDAV, JDA)
                  MA  = 0
                  DO II = 1, NASHT
                     IIBIT = IBTSHL(1, II - 1)
                     IF (IBTAND (JDA, IIBIT) .NE. 0) THEN
                        JAELEC = JAELEC + 1
                        JOCCI(JAELEC) = IPTA2O(II)
                     END IF
                     IF (IBTAND (IAI, IIBIT) .NE. 0) THEN
                        L  = KA
                        KA = IPTA2O(II)
                     END IF
                     IF (IBTAND (IAJ, IIBIT) .NE. 0) THEN
                        N  = MA
                        MA = IPTA2O(II)
                     END IF
                  END DO
                  JDAOLD = JDA
               END IF
               IF (NDAV .GT. ITMAX) GOTO 400
C
               IF (DODENCHK) THEN
                  DENCHK = ABS(CL(J,1)) + ABS(CR(J,1))
                  IF (DOIMAG)
     &               DENCHK = DENCHK + ABS(CL(J,2)) + ABS(CR(J,2))
                  IF (DENCHK .LT. THRSCR) GO TO 400
               ENDIF
C
               JDB = IDET(2,J)
               IF (NDAV .EQ. ITMAX) THEN
                  IF (JDB .NE. IDB) GO TO 400
                  IT = NDAV
                  IBI = 0
                  IBJ = 0
                  K = KA
                  M = MA
                  DO IELEC = JAELEC + 1, NAELEC
                     JOCCI(IELEC) = IOCCI(IELEC)
                  END DO
               ELSE
                  IDBV = IBTXOR (IDB, JDB)
                  CALL FNBITS(NDBV,IDBV,IBTAB)
                  IT = NDAV + NDBV
                  IF (IT .GT. ITMAX) GOTO 400
C
                  IELEC = JAELEC
                  IBI = IBTAND (IDBV, IDB)
                  IBJ = IBTAND (IDBV, JDB)
                  K = KA
                  M = MA
                  DO II = 1, NASHT
                     IIBIT = IBTSHL(1, II - 1)
                     IF (IBTAND (JDB, IIBIT) .NE. 0) THEN
                        IELEC = IELEC + 1
                        JOCCI(IELEC) = IPTB2O(II)
                     END IF
C
                     IF (IBTAND (IBI, IIBIT) .NE. 0) THEN
                        L = K
                        K = IPTB2O(II)
                     END IF
                     IF (IBTAND (IBJ, IIBIT) .NE. 0) THEN
                        N = M
                        M = IPTB2O(II)
                     END IF
                  END DO
               END IF
C
C              Hamiltonian matrix element
C              --------------------------
C
               IF (DOSIG) THEN
C
                  CALL GOSCI_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
     &                 HR,HI,HCORE,VMUUUU)
C
C                 S_I += H_{IJ} B_{J}
C
                  SVEC(I,1) = SVEC(I,1) + HR * BVEC(J,1)
                  IF (DOIMAG) THEN
                     SVEC(I,1) = SVEC(I,1) - HI * BVEC(J,2)
                     SVEC(I,2) = SVEC(I,2)
     &                    + HR * BVEC(J,2) + HI * BVEC(J,1)
                  END IF
C
C                 S_J += H_{IJ}^* B_{I}
C
                  SVEC(J,1) = SVEC(J,1) + HR * BVEC(I,1)
                  IF (DOIMAG) THEN
                     SVEC(J,1) = SVEC(J,1) + HI * BVEC(I,2)
                     SVEC(J,2) = SVEC(J,2)
     &                    + HR * BVEC(I,2) - HI * BVEC(I,1)
                  END IF
               END IF
C
C              Density matrix elements
C              -----------------------
C
               IF (DODEN) THEN
C
                  DMRIJ = CL(I,1)*CR(J,1)
                  DMRJI = CL(J,1)*CR(I,1)
                  IF (DOIMAG) THEN
                     DMRIJ = DMRIJ + CL(I,2)*CR(J,2)
                     DMIIJ = CL(I,1)*CR(J,2) - CL(I,2)*CR(J,1)
                     DMRJI = DMRJI + CL(J,2)*CR(I,2)
                     DMIJI = CL(J,1)*CR(I,2) - CL(J,2)*CR(I,1)
                     DOIMAG1 = (ABS(DMIIJ)+ABS(DMIJI)) .GT. THRZER
                  END IF
                  CALL GOSCI_PMAT(DO2EL,DOIMAG1,IT,IOCCI,JOCCI,K,L,M,N,
     &                 DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C
               END IF
 400        CONTINUE
 100     CONTINUE
C
 10   CONTINUE
C
      CALL GETTIM(CPU2,WLL2)
      CPUTID  = SECTID(CPU2-CPU1)
      WALLTID = SECTID(WLL2-WLL1)
      CALL GETNBITSCALLS(1,CNCALLS,FNCALLS)
D     write(lupri,*) '#calls to nbits: ',cncalls,fncalls
      WRITE(LUPRI,'(/5A)')
     &   ' CPU (Wall) time in GOSCI_CI kernel: ',CPUTID,' (',WALLTID,')'
      IF (NDSKIP .GT. 0) WRITE(LUPRI,'(A,I10,A,I10)')
     &   ' Number of configuration elements skipped:',NDSKIP,
     &   ' out of',NZCONF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_cisig_screen */
      SUBROUTINE GOSCI_CISIG_SCREEN(SVEC,BVEC,HCORE,VMUUUU,
     &                              DO2EL,IDET,THRSCR)
C***********************************************************************
C
C     Calculate:
C         CI sigma vector in SVEC with screening
C
C     Input:
C
C        BVEC - trial vector
C        HCORE  - core Fock matrix
C        VMUUUU - two electron integrals.
C
C        DO2EL - include 2-electron part : VMUUUU
C        IDET  - determinant strings
C        THRSCR - skip abs(BVEC(i)) .le. THRSCR
C
C     Output:
C
C        SVEC   - sigma vector
C
C
C     Written by J. Thyssen and H.J.Aa. Jensen - Jan-Mar 2001
C     Last revision : July 2002/Sep 2003 hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgosci.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
C
      DIMENSION SVEC(NZCONF,*), BVEC(NZCONF,*)
      DIMENSION VMUUUU(*), HCORE(*)
C
C
      LOGICAL   DO2EL
C
      DIMENSION IDET(2,*)
C
      LOGICAL   DOIMAG
C
      DIMENSION IOCCI(2*MAXASH_GOSCI), JOCCI(2*MAXASH_GOSCI)
      DIMENSION IBTAB(0:255)
      INTEGER*8 FNCALLS, CNCALLS
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
#include "ibtfun.h"
C
      CALL GETTIM(CPU1,WLL1)
C
C     Initialization
C     --------------
C
C     Bit table:
C
      CALL IBTABINI(IBTAB)
C
C     Zero counters
C
      CALL GETNBITSCALLS(1,FNCALLS,CNCALLS)
C
C     Zero vectors:
C
      MZ = MIN(NZ,2)
      CALL DZERO(SVEC,NZCONFQ)
C
      IF (DO2EL) THEN
C        max 4 differences between 2 determinants
C        for 2-electron terms
         ITMAX = 4
      ELSE
C        max 2 differences between 2 determinants
C        when only 1-electron terms
         ITMAX = 2
      END IF
C
      DOIMAG  = NZ.GT.1
C
      MK2STP = 2
C
      NBSKIP = 0
      DO 10 MK2I = MINMK2, MAXMK2, MK2STP
C        ... MK2 is NAEL - NBEL
         IF (NMK2(MK2I) .EQ. 0) GOTO 10
D        write(lupri,*) '(GOSCI_CISIG_SCREEN) MK2I = ',MK2I
         ISTRT = IMK2(MK2I) + 1
         IEND  = IMK2(MK2I) + NMK2(MK2I)
C
         MK2JSTRT = MAX(MINMK2,MK2I-ITMAX)
         MK2JEND  = MIN(MAXMK2,MK2I+ITMAX)
C
C        Loop over determinants
C        (use that outer loop in gosci_gendet is over A-strings,
C        inner loop is over B-strings)
C
         IDAOLD = -1
         DO 100 I = ISTRT, IEND
C
C        Screening:
C
            IF (ABS(BVEC(I,1)) .GT. THRSCR) GOTO 12
            IF (DOIMAG) THEN
               IF (ABS(BVEC(I,2)) .GT. THRSCR) GOTO 12
            END IF
            NBSKIP = NBSKIP + 1
            GOTO 100
C
C           ---------------
C
   12       IDA = IDET(1,I)
            IF (IDA .NE. IDAOLD) THEN
C              ... outer loop in gosci_gendet is over A-strings
               IAELEC = 0
               DO K = 1, NASHT
                  IF (IBTAND (IDA, IBTSHL (1, K - 1)) .NE. 0) THEN
                     IAELEC = IAELEC + 1
                     IOCCI(IAELEC) = IPTA2O(K)
                  END IF
               END DO
               IDAOLD = IDA
            END IF
C
            IELEC = IAELEC
            IDB = IDET(2,I)
            DO K = 1, NASHT
               IF (IBTAND (IDB, IBTSHL (1, K - 1)) .NE. 0) THEN
                  IELEC = IELEC + 1
                  IOCCI(IELEC) = IPTB2O(K)
               END IF
            END DO
            IF (IELEC .NE. NAELEC)
     &         CALL QUIT('GOSCI_CISIG_SCREEN: IELEC.ne.NAELEC')
C
C
C        *************************
C        *************************
C
C
            DO 20 MK2J = MK2JSTRT, MK2JEND, MK2STP
               IF (NMK2(MK2J) .EQ. 0) GOTO 20
D           write(lupri,*) '(GOSCI_CISIG_SCREEN) MK2I MK2J = ',MK2I,MK2J
D              call flshfo(lupri)
C
               JSTRT = IMK2(MK2J) + 1
               JEND  = IMK2(MK2J) + NMK2(MK2J)
               MXNDAV = ITMAX - ABS(MK2J - MK2I)/2
C              ... we know that NDBV .ge. ABS(MK2I-MK2J)/2
C
               JDAOLD = -1
               DO 200 J = JSTRT, JEND
C
                  JDA = IDET(1,J)
                  IF (JDA .NE. JDAOLD) THEN
C                    ... outer loop in gosci_gendet is over A-strings
                     IDAV = IBTXOR (IDA, JDA)
                     CALL FNBITS(NDAV,IDAV,IBTAB)
                     IF (NDAV .GT. MXNDAV) GOTO 200
C
                     JAELEC = 0
                     IAI = IBTAND (IDAV, IDA)
                     KA  = 0
                     IAJ = IBTAND (IDAV, JDA)
                     MA  = 0
                     DO II = 1, NASHT
                        IIBIT = IBTSHL (1, II - 1)
                        IF (IBTAND (JDA, IIBIT) .NE. 0) THEN
                           JAELEC = JAELEC + 1
                           JOCCI(JAELEC) = IPTA2O(II)
                        END IF
                        IF (IBTAND (IAI, IIBIT) .NE. 0) THEN
                           L  = KA
                           KA = IPTA2O(II)
                        END IF
                        IF (IBTAND (IAJ, IIBIT) .NE. 0) THEN
                           N  = MA
                           MA = IPTA2O(II)
                        END IF
                     END DO
                     JDAOLD = JDA
                  END IF
                  IF (NDAV .GT. MXNDAV) GOTO 200
C                 ... see note above about MXNDAV
C
                  JDB  = IDET(2,J)
                  IF (NDAV .EQ. ITMAX) THEN
                     IF (JDB .NE. IDB) GO TO 200
                     IT = NDAV
                     IBI = 0
                     IBJ = 0
                     K = KA
                     M = MA
                     DO IELEC = JAELEC + 1, NAELEC
                        JOCCI(IELEC) = IOCCI(IELEC)
                     END DO
                  ELSE
                     IDBV = IBTXOR (IDB, JDB)
                     CALL FNBITS(NDBV,IDBV,IBTAB)
C
                     IT = NDAV + NDBV
                     IF (IT .GT. ITMAX) GOTO 200
C                    ... IT = number of different spinors
C                          in det(I) and det(J)
C
                     IELEC = JAELEC
                     IBI = IBTAND (IDBV, IDB)
                     IBJ = IBTAND (IDBV, JDB)
                     K = KA
                     M = MA
                     DO II = 1, NASHT
                        IIBIT = IBTSHL (1, II - 1)
                        IF (IBTAND (JDB, IIBIT) .NE. 0) THEN
                           IELEC = IELEC + 1
                           JOCCI(IELEC) = IPTB2O(II)
                        END IF
C
                        IF (IBTAND (IBI, IIBIT) .NE. 0) THEN
                           L = K
                           K = IPTB2O(II)
                        END IF
                        IF (IBTAND (IBJ, IIBIT) .NE. 0) THEN
                           N = M
                           M = IPTB2O(II)
                        END IF
                     END DO
                  END IF
C
C                 Hamiltonian matrix element
C                 --------------------------
C
C
                  CALL GOSCI_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,
     &                    K,L,M,N,HR,HI,HCORE,VMUUUU)
C
C                 S_J += H_{IJ}^* B_{I}
C
                  SVEC(J,1) = SVEC(J,1) + HR * BVEC(I,1)
                  IF (DOIMAG) THEN
                     SVEC(J,1) = SVEC(J,1) + HI * BVEC(I,2)
                     SVEC(J,2) = SVEC(J,2)
     &                    + HR * BVEC(I,2) - HI * BVEC(I,1)
                  END IF
C              end loop 200 J and end loop 20 MKJ
 200           CONTINUE
 20         CONTINUE
C        end loop 100 I and end loop 10 MKI
 100     CONTINUE
 10   CONTINUE
C
      CALL GETTIM(CPU2,WLL2)
      CPUTID  = SECTID(CPU2-CPU1)
      WALLTID = SECTID(WLL2-WLL1)
D     CALL GETNBITSCALLS(1,CNCALLS,FNCALLS)
D     write(lupri,*) '#calls to nbits: ',cncalls,fncalls
      WRITE(LUPRI,'(/5A/A,I10,A,I10)')
     &   ' CPU (Wall) time in GOSCI_CISIG_SCREEN kernel: ',CPUTID,' (',
     &   WALLTID,')',' Number of BVEC elements skipped:',NBSKIP,
     &   ' out of',NZCONF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
