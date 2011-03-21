!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C VERSION : $Revision: 7375 $
C DATE    : $Date: 2008-08-15 16:30:18 +0200 (Fri, 15 Aug 2008) $
C FILE    : krmcgascip.F
C main authors: H.J.Aa. Jensen & J. Thyssen
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GASCIP_SETUP(UCIBOS,ICIREP,IPRINT)
C***********************************************************************
C
C     Setup and initialize variables in dcbgascip.h
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
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbbos.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbgascip.h"
#include "pgroup.h"
      LOGICAL UCIBOS
C
C
      CALL QENTER('GASCIP_SETUP')
C
C     Setup symmetry mult. table
C     --------------------------
C
C     GMULTSF/GMULTA returns:
C       NREP     : #(boson irreps) = #(fermion irreps);
C       REPNA    : strings identifying the irreps
C       MULTB    : mult. table
C
      IF (UCIBOS) THEN
C        Use boson symmetry (for spin free and Levy-Leblond)
C        or use approximate boson symmetry
         CALL GMULTSF(NREP,REPNA,MULTB(1,1))
C        hjaaj aug 2003: NOTE that this info is NOT used in this
C        version of GASCIP !!!!!!!!!!!!!
C        (MULTB is not used as it includes spin-symmetry -- see GMULTSF --
C        and we treat the spin Ms separately as in non-rel codes)
      ELSE
         CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB(1,1))
      ENDIF
C     We extend MULTB(1:2*NREP,1:2*NREP) with MULTB(1:2*NREP,0) for use in
C     GASCIP_SYMOK when we add the first spinor to a new
C     string (i.e. symmetry zero in second index means vacuum)
      DO I = 1, 2 * NREP
         MULTB(I,0) = I
      END DO
C
C     Setup symmetry
C     --------------
C
C     IRP  : ICIREP is user input for CI symmetry
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
      IRP = ICIREP
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
C      which implies UCIBOS true, or use approx. boson symmetry
C      if UCIBOS true for relativistic calculation)
C
      II    = 0
      IREPA = 0
      DO I = 1, NFSYM
         IREPA = IREPA + 1
         DO ISTR = 1, NASH(I)
            II = II + 1
            IF (UCIBOS) THEN
               IG = IDXU2G(IASH(I)+ISTR)
               IRPAMO(II) = IBOSYM(IG)
            ELSE
               IRPAMO(II) = IREPA
            END IF
         ENDDO
         IF (NZ.NE.4.OR.SPINFR.or.levyle) IREPA = IREPA + 1
C        for NZ.eq.4: unbarred and barred in same fermion irrep !
C        (see GMULTA for def. of fermion irreps)
         DO ISTR = 1, NASH(I)
            II = II + 1
            IF (UCIBOS) THEN
               IG = IDXU2G(IASH(I)+ISTR)
               IRPAMO(II) = IBOSYM(IG)
C              to use spin symmetry here with MULTB from GMULTSF
C              we should have had NBSYM + IBOSYM(IG) here,
C              but we take care of the spin separately in GASCIP,
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
      IF (IPRINT .GE. 2 ) THEN
         WRITE(LUPRI,'(A,3I3/A)')
     &      ' (GASCIP_SETUP) NREP, IRP, IRRP = ', NREP, IRP, IRRP,
     &      ' -> IRRP.lt.0 means include all symmetries'
         IF (IRRP.GT.0) WRITE(LUPRI,'(2A)')
     &      ' (GASCIP_SETUP) REPNA(IRRP)     = ', REPNA(IRRP)
         WRITE(LUPRI,'(/A)')
     &      ' Symmetry of active orbitals (unbarred, barred):'
         DO I = 1,NASHT
            IA = IPTA2O(I)
            IB = IPTB2O(I)
            IF (.NOT. UCIBOS) THEN
               WRITE(LUPRI,'(I15,2(5X,A))')
     &         I, REPNA(IRPAMO(IA)), REPNA(IRPAMO(IB))
            ELSE
               WRITE(LUPRI,'(I15,2(5X,A))')
     &         I, REP(IRPAMO(IA)), REP(IRPAMO(IB))
            END IF
         END DO
      END IF
C
      CALL QEXIT('GASCIP_SETUP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_gendet */
      SUBROUTINE GASCIP_GENDET(NDET,IDET,IASTR,IBSTR,MAXDET,UCIBOS,
     &                        ICIREP,IPRINT)
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
#include "dcbgascip.h"
#include "dcborb.h"
#include "pgroup.h"
#include "dgroup.h"
C
      DIMENSION NA(MAXASH_GASCIP), NB(MAXASH_GASCIP)
      DIMENSION IDET(2,*), IASTR(*), IBSTR(*)
      LOGICAL GASCIP_GASOK,GASCIP_SYMOK,UCIBOS
C
      INTEGER MSKGAS(MXMCGAS), IBTAB(0:255)
C
#include "ibtfun.h"
C
C     Use GASCIP:
C
      CALL GASCIP_SETUP(UCIBOS,ICIREP,IPRINT)
C
      IF ( IPRINT .GE. 3 ) THEN
         WRITE(LUPRI,'(/A)')
     &       ' (GASCIP_GENDET) Output from generation of determinants:'
         IF (UCIBOS) WRITE(LUPRI,'(15X,A)')
     &    'Using boson spinor symmetry (exact or approximate) in the CI'
         IF (IRRP .LE. 0) THEN
            WRITE(LUPRI,'(15X,A)')
     &        'Generating determinants of all symmetries'
         ELSE IF (UCIBOS) THEN
            WRITE(LUPRI,'(15X,3A,I4,A,I3)')
     &        'Generating dets. with spatial symmetry ', REP(IRP-1),
     &        ' and 2*MS interval:',MINMK2,' to',MAXMK2
         ELSE
            WRITE(LUPRI,'(15X,2A)')
     &        'Generating dets. with symmetry ', REPNA(IRRP)
         END IF
         WRITE(LUPRI,'(15X,A,I10)')
     &        'Max no. of determinants',MAXDET
      END IF
C
      IF (NASHT .GT. MAXASH_GASCIP) THEN
         WRITE(LUPRI,'(//A,I3,A,I3)')
     &      ' *** ERROR in GASCIP_GENDET *** Number of active'//
     &      ' orbitals (',NASHT,') must be <',MAXASH_GASCIP+1
         CALL QUIT('*** ERROR in GASCIP_GENDET ***')
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
      IF (NGAS_DC .GT. 0) THEN
C
         CALL IZERO(MSKGAS,MXMCGAS)
C
         L = 1
         DO I = 1, NFSYM
            DO J = 1, NGAS_DC
               DO K = 1, NGSH(I,J)
                  MSKGAS(J) = IBTOR(MSKGAS(J),L)
                  L = L * 2
               END DO
            END DO
         END DO
C
#ifdef MOD_DEBUG
         DO I = 1, NGAS_DC
            write(lupri,*) 'gas bit mask before: ',i,mskgas(i)
         END DO
#endif
C
         DO I = 2, NGAS_DC
            MSKGAS(I) = MSKGAS(I-1) + MSKGAS(I)
         END DO
C
#ifdef MOD_DEBUG
         DO I = 1, NGAS_DC
            write(lupri,*) 'gas bit mask after: ',i,mskgas(i)
         END DO
#endif
C
      END IF
C
C     Find min and max number of a & b electrons
C     ------------------------------------------
C
C     write(lupri,*) 'maxmk2, minmk2',maxmk2,minmk2
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
C     write(lupri,*) 'IAEL   IBEL'
      NDET = 0
      DO 10 IAEL = MINAEL, MAXAEL
         IBEL = NAELEC - IAEL
         MK2 = IAEL - IBEL
         IMK2(MK2) = NDET
C        write(lupri,'(I4,I7)') IAEL,IBEL
C        call flshfo(lupri)
C
         CALL GASCIP_GENSTR(IAEL,IASTR,NASTR,MSKGAS,MAXDET)
         IF (NASTR .EQ. 0) GOTO 10
         CALL GASCIP_GENSTR(IBEL,IBSTR,NBSTR,MSKGAS,MAXDET)
         IF (NBSTR .EQ. 0) GOTO 10
C        write(lupri,'(12X,I8,I11,I11)') NASTR,NBSTR,NDET
C        call flshfo(lupri)
C
C HJ TODO : the next can easily be optimized by not redoing the IA
C           for each IB in GASCIP_GASOK and GASCIP_SYMOK /aug 03
C
         DO IA = 1, NASTR
            DO IB = 1, NBSTR
               IF (GASCIP_GASOK(IASTR(IA),IBSTR(IB),MSKGAS,IBTAB)) THEN
                  IF (GASCIP_SYMOK(IASTR(IA),IBSTR(IB),UCIBOS)) THEN
                     NDET = NDET + 1
                     IF (NDET .GE. MAXDET) THEN
                        WRITE(LUPRI,'(/A,I12/A)')
     &                     ' *** ERROR in GASCIP_GENDET *** '//
     &                       'Maximum number of determinants reached!',
     &                       MAXDET,
     &                     'Allocate more memory!'
                        CALL QUIT('*** ERROR in GASCIP_GENDET ***')
                     END IF
                     IDET(1,NDET) = IASTR(IA)
                     IDET(2,NDET) = IBSTR(IB)
                  END IF
               END IF
            END DO
         END DO
C        write(LUPRI,'(T32,I11)') NDET
C        call flshfo(LUPRI)
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
      IF (IPRINT .GE. 10 .OR. NDET .LE. 50) THEN
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
C  /* Deck gascip_gasok */
      LOGICAL FUNCTION GASCIP_GASOK(IA,IB,MSKGAS,IBTAB)
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
#include "dcbgascip.h"
#include "dcborb.h"
C
      DIMENSION MSKGAS(*), IBTAB(*)
C
#include "ibtfun.h"
C
      GASCIP_GASOK = .TRUE.
      IF (NGAS_DC .EQ. 0) RETURN
C
      DO IGAS = 1, NGAS_DC
         IGA = IBTAND(MSKGAS(IGAS),IA)
         IGB = IBTAND(MSKGAS(IGAS),IB)
         CALL FNBITS(NA,IGA,IBTAB)
         CALL FNBITS(NB,IGB,IBTAB)
         NG = NA + NB
         IF ( NG .LT. NGSP(1,IGAS) .OR. NG .GT. NGSP(2,IGAS)) THEN
            GASCIP_GASOK = .FALSE.
            RETURN
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_krlist */
      SUBROUTINE GASCIP_KRLIST(NDET,IDET,KRLIST,IPRINT)
C***********************************************************************
C
C     For each determinant find its Kramers partner.
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
         CALL HEADER('Output from GASCIP_KRLIST',-1)
         WRITE(LUPRI,'(T10,A,T30,A)') 'Determinant','Kramers partner'
         DO I = 1, NDET
            WRITE(LUPRI,'(T10,I9,T30,I9)') I,KRLIST(I)
         END DO
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ngascip_nbits */
      INTEGER FUNCTION NGASCIP_NBITS(JDET,IBTAB)
C***********************************************************************
C
C     Find number of bits in integer JDET.
C
C     Written by J. Thyssen - Jan 2 2001
C     Last revision : S. Knecht - Nov 17 2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgascip.h"
#include "dcborb.h"
C
      INTEGER   IBTAB(0:255)
C
C     17-11-2007-SK: INTEGER*1 is an illegal value on 64-Bit SX machines
C     Is this solution still working?
#if defined (SYS_SX64)
!     INTEGER*4 ID1(4)
      call qenter('NGASCIP_NBITS')
      call quit('ERROR: INTEGER*1 was not defined on 64-bit SX machines')
#else
      INTEGER*1 ID1(4)
      INTEGER*4 ID4
      EQUIVALENCE (ID1(1),ID4)
#endif /* SYS_SX64 */
C
C     Anyone who can find with a faster routine?
C
      ID4 = JDET
C     write(lupri,*) 'jdet',jdet,' and id4',id4
C     write(lupri,*) 'id1',id1
      NGASCIP_NBITS
     &     = IBTAB(ID1(1))
     &     + IBTAB(ID1(2))
     &     + IBTAB(ID1(3))
     &     + IBTAB(ID1(4))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_symok */
      LOGICAL FUNCTION GASCIP_SYMOK(IA,IB,UCIBOS)
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
#include "dcbgascip.h"
#include "dcborb.h"
C
      LOGICAL UCIBOS
#include "ibtfun.h"
C
      IF (IRP .LE. 0) THEN
C        ... include all symmetries
         GASCIP_SYMOK = .TRUE.
      ELSE IF (.NOT. UCIBOS) THEN
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
         GASCIP_SYMOK = (JRRP .EQ. IRRP)
      ELSE
C
C        Boson symmetry (spin free or Levy-Leblond).
C        (do not use IRRP and MULTB which includes spin-symmetry -- see GMULTSF --
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
         GASCIP_SYMOK = (JRRP .EQ. IRP-1)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_startvec */
      SUBROUTINE GASCIP_STARTVEC(IDET,CMO,CREF,ECI,clvcorr,NDET,ISTATE,
     &                           THRPCI,CIPROG,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors from GASCIP by diagonalizing the
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
#include "dgroup.h"
C
      CHARACTER CIPROG*8
      DIMENSION WORK(*), CREF(NDET,*), CMO(*), IDET(*)
C
      CALL QENTER('GASCIP_STARTVEC')
C
      KFRSAV = KFREE
      MZ = MIN(2,NZ)
      N2DET = NDET*NDET
C
      CALL MEMGET('REAL',KEIGVL,NDET,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,N2DET*MZ,WORK,KFREE,LFREE)
C
C
C     Do full CI
C
      CALL GASCIP_FULLCI(NDET,IDET,CMO,WORK(KEIGVL),WORK(KEIGVC),
     &                  ECORE,IPRINT,WORK,KFREE,LFREE)
C
C     Get CI vector:
C     --------------
C
C     The CI vector is stored in columns; the lowest eigenstate is in
c     the first column.
C
      IOFF = KEIGVC - 1 + (ISTATE-1) * NDET
      IF ( NZ .EQ. 1 ) THEN
         DO I = 1,NDET
            CREF(I,1) = WORK(IOFF + I)
         END DO
      ELSE
         DO I = 1,NDET
            CREF(I,1) = WORK(IOFF + I )
            CREF(I,2) = WORK((IOFF + N2DET) + I )
         END DO
      END IF
C
C DEBUG DEBUG DEBUG
c      call dzero(cref,ndetq)
c      if (nz .eq. 1) then
c      cref(2,1) = 1.0d0
c      cref(2,2) = 1.0d0/sqrt(2.0d0)
c      else
c         CREF(2,2) = 1.0D00
c      endif
c      DO I = 2, NDET
c         CREF(I,1) = 0.0D00
c      END DO
c      IF ( NDET .GE. 2) THEN
c         CREF(1,1) = 0.8
c         CREF(2,1) = 0.6
c      END IF
c      CREF(1,1) = 0.8D00
c      IF ( NDET .GE. 2) CREF(2,1) = 0.6D00
c      DO I = 3, NDET
c         CREF(I,1) = 0.0D00
c      END DO
C
C
C     Get CI energy
C
      ECI = WORK(KEIGVL + ISTATE - 1)
!
!     add LVCORR correction (zero if flag is not active)
      ECORE = ECORE + clvcorr
      WRITE(LUPRI,'(/A,F20.8)')
     &     ' (GASCIP_STARTVEC) Core energy    = ',ECORE
      WRITE(LUPRI,'(A,F20.8)')
     &     '                  Active energy  = ',ECI-ECORE
      WRITE(LUPRI,'(A,F20.8)')
     &     '                  Total energy   = ',ECI
C
C     Check norm of CI vector
C
      CALL RCINORM(CREF,10)
C
      IF (IPRINT .GE. 5 ) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_STARTVEC) CI start vector'
         CALL RPRCI(CREF,NDET,NZ,CIPROG, THRPCI,LUPRI)
      END IF
C
      CALL MEMREL('GASCIP_STARTVEC',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GASCIP_STARTVEC')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_fullci */
      SUBROUTINE GASCIP_FULLCI(NDET,IDET,CMO,EIGVAL,EIGVEC,ECORE,
     &                        IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate CI start vectors from GASCIP by diagonalizing the
C     full CI matrix.
C
C     Input : IDET   - determinant information
C             CMO    - MO coefficents
C
C     Output: EIGVAL - CI eigenvalues
C             EIGVEC - CI eigenvectors
C             ECORE  - core energy
C
C     Extracted from GASCIP_STARTVEC July 2002 /hjaaj.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION WORK(*), EIGVAL(*), EIGVEC(*), CMO(*), IDET(*)
C
      CALL QENTER('GASCIP_FULLCI')
C
      KFRSAV = KFREE
      MZ = MIN(2,NZ)
      N2DET = NDET * NDET
      N2DETQ = N2DET * MZ
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
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPRINT)
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPRINT)
C
      CALL MEMREL('GASCIP_FULLCI after RGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * MZ, WORK,KFREE,LFREE)
C
      LH2AC = NASHT*NASHT*NNASHX*NZ*3
      CALL MEMGET('REAL',KH2AC,LH2AC,WORK,KFREE,LFREE)
C
C     read integrals from 4IND*: out --> work(kh2ac)
C     ==============================================

!     get dimensions for scracth matrices
      call rgeth2_dim(ndim_ij,ndim_kl,ndim_h12_tmp)

      call memget('REAL',khmat,ndim_h12_tmp,work,kfree,lfree)
      call memget('REAL',kh2xy,ndim_h12_tmp,work,kfree,lfree)
      call memget('REAL',kh2yx,ndim_h12_tmp,work,kfree,lfree)

      call memget('INTE',kdxb12, 2*ndim_ij,work,kfree,lfree)
      call memget('INTE',kdxb34, 2*ndim_kl,work,kfree,lfree)
      call memget('INTE',kdxkr12,2*ndim_ij,work,kfree,lfree)

      call memget('INTE',kibeig,norbt,work,kfree,lfree)

      call izero(work(kibeig),norbt)
      if(spinfr.or.levyle)then
        call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
      else if(linear)then
        call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
      end if

      call rgeth21(.false.,.true.,dummy,dummy,work(kh2ac),
     &             work(khmat),work(kh2xy),work(kh2yx),
     &             work(kdxkr12),work(kdxb12),work(kdxb34),
     &             work(kibeig),iprint,lupri)
C
C     Transform integrals to Molfdir format:
C
C     remove potential "garbage" from imaginary part...
      IF (NZ .EQ. 1) CALL DZERO(WORK(KMUUUU+(2*NASHT)**4),(2*NASHT)**4)

      CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPRINT)
      CALL MEMREL('GASCIP_FULLCI.H2AC',WORK,1,KH2AC,KFREE,LFREE)
C
C     Allocate memory for matrices:
C
      CALL MEMGET('REAL',KHR,N2DETQ,WORK,KFREE,LFREE)
      KHI = KHR + N2DET
C
C     Construct Hamiltonian:
C     ----------------------
C
      CALL GASCIP_MAKEH(.TRUE.,1,NDET,IDET,WORK(KHR),WORK(KHI),
     &                 WORK(KFCACM),WORK(KMUUUU))
C
C
C     Diagonalize Hamiltonian:
C     ------------------------
C
      IF ( IPRINT .GE. 20 ) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_FULLCI) CI matrix - real part:'
         CALL OUTPUT(WORK(KHR),1,NDET,1,NDET,NDET,NDET,1,LUPRI)
        IF (MZ .GT. 1) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_FULLCI) CI matrix - imag part:'
         CALL OUTPUT(WORK(KHI),1,NDET,1,NDET,NDET,NDET,1,LUPRI)
        END IF
      END IF
      CALL FLSHFO(LUPRI)
C
C
      CALL GETTIM(CPU1,WALL1)
      IF (MZ .EQ. 1) THEN
         IJOB = 1
         IORDER = 1
         IPACK = 0
         CALL RSJACO(NDET,NDET,NDET,WORK(KHR),EIGVAL,
     &        IJOB,IORDER,IPACK,EIGVEC,WORK(KFREE),LFREE)
      ELSE
         MATZ = 1
         CALL QDIAG(2,NDET,WORK(KHR),NDET,NDET,
     &        EIGVAL,MATZ,EIGVEC,NDET,NDET,
     &        WORK(KFREE),LFREE,IERR)
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in GASCIP_FULLCI ***: ',
     &           'QDIAG failed with error code ',IERR
            CALL QUIT('*** ERROR in GASCIP_FULLCI *** QDIAG failed')
         END IF
      END IF
      CALL GETTIM(CPU2,WALL2)
      WRITE (LUPRI,'(/A,2F20.2)')
     &   "CPU and WALL times for diagonalization:",CPU2-CPU1,WALL2-WALL1
C
C     Add core energy
C
      DO I = 1,NDET
         EIGVAL(I) = EIGVAL(I) + ECORE
      END DO
C
      IF ( IPRINT .GE. 5 ) THEN
         WRITE(LUPRI,'(/A)') ' (GASCIP_FULLCI) CI eigenvalues'
         CALL OUTPUT(EIGVAL,1,1,1,NDET,1,NDET,1,LUPRI)
       IF ( IPRINT .GE. 20 .OR. NDET .LE. 12) THEN
         WRITE(LUPRI,'(/A)')
     &   ' (GASCIP_FULLCI) CI eigenvectors (real part)'
         CALL OUTPUT(EIGVEC,1,NDET,1,NDET,NDET,NDET,1,LUPRI)
         IF (NZ .GT. 1) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (GASCIP_FULLCI) CI eigenvectors (imag part)'
            CALL OUTPUT(EIGVEC(1+N2DET),1,NDET,1,NDET,
     &           NDET,NDET,1,LUPRI)
         END IF
       END IF
      END IF
C
      CALL MEMREL('GASCIP_FULLCI',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GASCIP_FULLCI')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gascip_ireps */
      SUBROUTINE GASCIP_IREPS (NZ,NFSYM,NSYMRPA,REPA)
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
C  /* Deck gascip_rmakdm */
      SUBROUTINE GASCIP_RMAKDM(NDET,CALCDV,CALCPV,IDET,CL,CR,DV,PV,
     &                        IPRINT)
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
#include "dcbgascip.h"
#include "dgroup.h"
#include "dcborb.h"
C
      DIMENSION CL(NDET,NZ),CR(NDET,NZ),DV(2*NASHT,2*NASHT,2)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,2)
      DIMENSION IDET(*)
C
      LOGICAL CALCDV, CALCPV
C
C
      CALL QENTER('GASCIP_RMAKDM')
C
      IF (.NOT. CALCDV) THEN
         WRITE(LUPRI,'(/A)')
     $        ' *** ERROR in GASCIP_RMAKDM *** ' //
     $        'CALCDV false is not implemented.'
         CALL QUIT('*** ERROR in GASCIP_RMAKDM ***')
      END IF
C
C     Calculate density matrices
C     --------------------------
C
      CALL GASCIP_MAKETDM(CALCPV,NDET,IDET,CL,CR,DV,PV,IPRINT)
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
     &           ' (GASCIP_RMAKDM) real(DV) in Molfdir basis'
            CALL OUTPUT(DV(1,1,1),1,2*NASHT,1,2*NASHT,
     &           2*NASHT,2*NASHT,1,LUPRI)
            IF ( NZ .GE. 2 ) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (GASCIP_RMAKDM) imag(DV) in Molfdir basis'
               CALL OUTPUT(DV(1,1,2),1,2*NASHT,1,2*NASHT,
     &              2*NASHT,2*NASHT,1,LUPRI)
            END IF
         END IF
C
C        Print PV
C
         IF (CALCPV) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RMAKDM) PV in Molfdir basis'
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
      CALL QEXIT('GASCIP_RMAKDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_maketdn */
      SUBROUTINE GASCIP_MAKETDM(DOPV,NDET,IDET,CL,CR,DV,PV,IPRINT)
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
      CALL QENTER('GASCIP_MAKETDM')
C
      CALL GASCIP_CI(
     &     .FALSE.,DUMMY,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,CL,CR,DV,PV,
     &     DOPV,NDET,IDET,IPRINT)
C
      CALL QEXIT('GASCIP_MAKETDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_pmat */
      SUBROUTINE GASCIP_PMAT(DOPV,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
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
#include "dcbgascip.h"
C
      LOGICAL   DOPV, DOIMAG
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
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
      CALL QUIT('GASCIP_PMAT called with IT .gt. 4 !')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_maketdn */
      SUBROUTINE GASCIP_MAKEH(DO2EL,IHRM,NDET,IDET,HR,HI,HCORE,VMUUUU)
C***********************************************************************
C
C     Construct Hamiltonian for use in gascip.
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
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
      DIMENSION HCORE(*), IDET(2,*)
      LOGICAL   DO2EL,DOIMAG
      DIMENSION IBTAB(0:255)
C
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      LOGICAL HFULL
C
#include "ibtfun.h"
C
      CALL IBTABINI(IBTAB)
      HFULL = IHRM .EQ. 0
      DOIMAG  = NZ.GT.1
      CALL DZERO(HR,NDET*NDET)
      IF (DOIMAG) CALL DZERO(HI,NDET*NDET)
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
      DO I = 1, NDET
C
C        Use that outer loop in gascip_gendet is over A-strings,
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
         IF (IELEC .NE. NAELEC) THEN
            WRITE(LUPRI,*) 'GASCIP_CI: IELEC.ne.NAELEC',IELEC,NAELEC
            WRITE(LUPRI,'(A,B64)') 'IDA ',IDA
            WRITE(LUPRI,'(A,B64)') 'IDB ',IDB
            CALL QUIT('GASCIP_MAKEH: IELEC.ne.NAELEC')
         END IF
C
         IF (HFULL) THEN
            JEND = NDET
         ELSE
chj aug 2002: code generalized to non-hermitian operators
            CALL GASCIP_HMAT(DO2EL,DOIMAG,0,IOCCI,IOCCI,0,0,0,0,
     &           HR(I,I),HI(I,I),HCORE,VMUUUU)
            JEND = I - 1
         END IF
         JDAOLD = -1
         DO J = 1, JEND
C
C           Use that outer loop in gascip_gendet is over A-strings,
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
               CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
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
c     write(LUPRI,*) 'real part of HMAT in GASCIP_MAKEH ='
c     call output(hr,1,nzconf,1,nzconf,nzconf,nzconf,1,lupri)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_hmat */
      SUBROUTINE GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
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
#include "dcbgascip.h"
#include "dcborb.h"
C
      LOGICAL   DO2EL,DOIMAG
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
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
      CALL QUIT('GASCIP_HMAT called with IT .gt. 4 !')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_svc */
      SUBROUTINE GASCIP_SVC(SAVDIA,DIAG,SVEC,NDET,IDET,BVEC,
     &                      HCORE,VMUUUU,IPRINT)
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
      CALL QENTER('GASCIP_SVC')
C
      IF (NDET .GT. 1000) THEN
C hj 16-jun-2004 new sigma code with screening
         THRSCR = 1.0D-10
         CALL GASCIP_CISIG_SCREEN(SVEC,BVEC,HCORE,VMUUUU,.TRUE.,
     &      NDET,IDET,THRSCR,IPRINT)
         IF (SAVDIA) CALL GASCIP_CI(SAVDIA,DIAG,
     &     .FALSE.,SVEC,BVEC,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IPRINT)
      ELSE
         CALL GASCIP_CI(SAVDIA,DIAG,
     &     .TRUE.,SVEC,BVEC,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IPRINT)
      END IF
C
      CALL QEXIT('GASCIP_SVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_hdiag */
      SUBROUTINE GASCIP_HDIAG(DIAG,NDET,IDET,HCORE,VMUUUU,IPRINT)
C***********************************************************************
C
C     Calculate diagonal of CI Hamiltonian.
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
      CALL QENTER('GASCIP_HDIAG')
C
      CALL GASCIP_CI(.TRUE.,DIAG,
     &     .FALSE.,DUMMY,DUMMY,HCORE,VMUUUU,
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IPRINT)
C
      CALL QEXIT('GASCIP_HDIAG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_rcistd */
      SUBROUTINE GASCIP_RCISTD(NDET,IDET,CMO,CREF,ECI,clvcorr,NCIROOT,
     &                        ISTATE,MAXCIT,THRPCI,CIPROG,
     &                        IPRINT,WORK,LWORK)
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
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CHARACTER CIPROG*8
      DIMENSION WORK(*), CREF(NDET,*), CMO(*), IDET(2,*)
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
#include "memint.h"
C
      CALL QENTER('GASCIP_RCISTD')
C
      CALL GETTIM(CPU1,WALL1)
C
      IF (NCIROOT .GT. 1) THEN
         CALL QUIT('.CIROOT .gt. 1 has not been implemented yet.')
      END IF
C
      MZ = MIN(NZ,2)
C
C     ************************************
C     *** Calculate necesary integrals ***
C     ************************************
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET('REAL',KFCACM,(2*NASHT)*(2*NASHT)*MZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Calculate FCmo (TOFILE false: do not write to LUKRM3) :
C
      CALL RGETFC(.FALSE.,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
C     Get active-active part of FCmo:
C
!
!     add LVCORR correction (zero if flag is not active)
      ECORE = ECORE + clvcorr

      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPRINT)
      IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(A)') 'Active Fock matrix:'
        DO I=1,N2ASHXQ,1
          WRITE(6,'(1I3,1F18.14)') I,WORK(KFCAC+I-1)
        END DO
      END IF
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPRINT)
C
      CALL MEMREL('RCISTD after RGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,WORK,KFREE,LFREE)
C
C     read integrals from 4IND*: out --> work(kh2ac)
C     ==============================================

!     get dimensions for scracth matrices
      call rgeth2_dim(ndim_ij,ndim_kl,ndim_h12_tmp)

      call memget('REAL',khmat,ndim_h12_tmp,work,kfree,lfree)
      call memget('REAL',kh2xy,ndim_h12_tmp,work,kfree,lfree)
      call memget('REAL',kh2yx,ndim_h12_tmp,work,kfree,lfree)

      call memget('INTE',kdxb12, 2*ndim_ij,work,kfree,lfree)
      call memget('INTE',kdxb34, 2*ndim_kl,work,kfree,lfree)
      call memget('INTE',kdxkr12,2*ndim_ij,work,kfree,lfree)

      call memget('INTE',kibeig,norbt,work,kfree,lfree)

      call izero(work(kibeig),norbt)
      if(spinfr.or.levyle)then
        call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
      else if(linear)then
        call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
      end if

      call rgeth21(.false.,.true.,dummy,dummy,work(kh2ac),
     &             work(khmat),work(kh2xy),work(kh2yx),
     &             work(kdxkr12),work(kdxb12),work(kdxb34),
     &             work(kibeig),iprint,lupri)

C
C     Make sure imag(VMUUUU) = 0 for NZ=1
C     (DNZ32M only DZERO'es real part for NZ = 1)
C
      IF (NZ .EQ. 1) CALL DZERO(WORK(KMUUUU+(2*NASHT)**4),(2*NASHT)**4)
C
C     Transform integrals to Molfdir format:
C
      CALL DNZ32M(WORK(KH2AC),WORK(KMUUUU),IPRINT)
      CALL MEMREL('GASCIP_RCISTD.H2AC',WORK,KWORK,KH2AC,KFREE,LFREE)
C
C
C     ********************************
C     *** Find start trial vectors ***
C     ********************************
C
C
      NDETQ = NDET*MZ
      CALL MEMGET('REAL',KHDIAG,NDET,WORK,KFREE,LFREE)
C
C     Assume we can have all sigma and trial vectors in memory.
C     This could get a little expensive for large CI's.
C
      MAX_HRED = MAX(50,MAXCIT+4)
C     MAXCIT+4: allocating space for up to 5-double degeneracy in diagonal
      CALL MEMGET('REAL',KBCVEC,MAX_HRED*NDETQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCVEC,MAX_HRED*NDETQ,WORK,KFREE,LFREE)
C
C     Calculate diagonal of CI matrix:
C     --------------------------------
C
c     CALL GASCIP_HDIAG(WORK(KHDIAG),NDET,IDET,WORK(KFCACM),WORK(KMUUUU))
      CALL GASCIP_CI(.TRUE.,WORK(KHDIAG),
     &     .FALSE.,DUMMY,DUMMY,WORK(KFCACM),WORK(KMUUUU),
     &     .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &     .TRUE.,NDET,IDET,IPRINT)
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(A)') ' (GASCIP_RCISTD) Diagonal of Hamiltonian'
         CALL OUTPUT(WORK(KHDIAG),1,NDET,1,1,NDET,1,1,LUPRI)
      END IF
C
C     Find lowest elements in HDIAG:
C
      II = IDMIN(NDET,WORK(KHDIAG),1)
C
C     Set first trial vector equal to:
C     BCVEC(I) = DELTA(I,II) / HDIAG(II)
C     add trial vectors for degenerate determinants within 1.0D-5
C
      HDIAGTST = WORK(KHDIAG-1+II) + 1.0D-5
      N_CSIM = 0
      DO  II = 1, NDET
          IF (WORK(KHDIAG-1+II) .LE. HDIAGTST) THEN
             JBCVEC = KBCVEC +  N_CSIM*NDETQ
             CALL DZERO(WORK(JBCVEC),NDETQ)
             WORK(JBCVEC + (II-1)) = D1
             CALL RCINORM(WORK(JBCVEC),IPRINT)
             N_CSIM = N_CSIM + 1
          END IF
      END DO
C
C
C
      MZ = MIN(2,NZ)
C
      JCRED = 0
      N_CRED = N_CSIM
      RTHRCI = 1.0D-5
c     RTHRCI = 1.0D-12
c     RTHRCI = 1.0D-6
      ECI = D0
C     GASCIP screening, a calculation in C2 symmetry showed that
C     THRSCR = RTHRCI**2 was not good enough for the complex group /hjaaj July 04
      THRSCR = RTHRCI**2 * 0.01D0
C
      LHRED = MAX_HRED*MAX_HRED*MZ
      CALL MEMGET('REAL',KHRED,LHRED,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVL,MAX_HRED,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,LHRED,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHRED),LHRED)
      KHREDI = KHRED + MAX_HRED*MAX_HRED
      KEIGVCI = KEIGVC + MAX_HRED*MAX_HRED
C
      ITMIC = 0
      ICONV_TOT = 0
      CALL GETTIM(CPU2,WALL2)
C
C        **********************************
C        *** Start micro iteration loop ***
C        **********************************
C
 100  CONTINUE  ! DO WHILE (not converged and no limits reached)
         ITMIC = ITMIC + 1
         CALL GETTIM(CPU3,WALL3)
C
         WRITE(LUPRI,'(//A,I4)')
     &        ' (GASCIP_RCISTD) CI microiteration no.',ITMIC
         CALL FLSHFO(LUPRI)
C
C        Calculate sigma vector(s)
C        -------------------------
C
         DO I = 1, N_CSIM
            IF (NDET .GT. 1000) THEN
C hj 16-jun-2004 new sigma code with screening
               CALL GASCIP_CISIG_SCREEN(
     &           WORK(KSCVEC+(I-1+JCRED)*NDETQ),
     &           WORK(KBCVEC+(I-1+JCRED)*NDETQ),
     &           WORK(KFCACM),WORK(KMUUUU),
     &           .TRUE.,NDET,IDET,THRSCR,IPRINT)
            ELSE
               CALL GASCIP_CI(.FALSE.,DUMMY,
     &           .TRUE.,WORK(KSCVEC+(I-1+JCRED)*NDETQ),
     &           WORK(KBCVEC+(I-1+JCRED)*NDETQ),
     &           WORK(KFCACM),WORK(KMUUUU),
     &           .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &           .TRUE.,NDET,IDET,IPRINT)
            END IF
c
            IF (IPRINT .GE. 10) THEN
C
C              Trial vector
C
               WRITE(LUPRI,'(/A,I3)')
     &              ' (GASCIP_RCISTD) Trial vector no. ',I
               CALL OUTPUT(WORK(KBCVEC+(I-1+JCRED)*NDETQ),
     &              1,NDET,1,MZ,NDET,MZ,1,LUPRI)
C
C              Sigma vector
C
               WRITE(LUPRI,'(/A,I3)')
     &              ' (GASCIP_RCISTD) Sigma vector no. ',I
               CALL OUTPUT(WORK(KSCVEC+(I-1+JCRED)*NDETQ),
     &              1,NDET,1,MZ,NDET,MZ,1,LUPRI)
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
            DO I = JCRED+1, N_CRED
               WORK(KHRED+(I-1)*MAX_HRED+(J-1)) =
     &              DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KSCVEC+(J-1)*NDETQ),1)
               IF (MZ .GT. 1) THEN
                  WORK(KHREDI+(I-1)*MAX_HRED+(J-1)) =
     &              DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KSCVEC+(J-1)*NDETQ+NDET),1)
     &            - DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ+NDET),1,
     &              WORK(KSCVEC+(J-1)*NDETQ),1)
               END IF
            END DO
         END DO
C
C        Old trial vectors * New sigma vectors
C
         DO J = JCRED+1, N_CRED
            DO I = 1, JCRED
               WORK(KHRED+(I-1)*MAX_HRED+(J-1)) =
     &              DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KSCVEC+(J-1)*NDETQ),1)
               IF (MZ .GT. 1) THEN
                  WORK(KHREDI+(I-1)*MAX_HRED+(J-1)) =
     &              DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KSCVEC+(J-1)*NDETQ+NDET),1)
     &            - DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ+NDET),1,
     &              WORK(KSCVEC+(J-1)*NDETQ),1)
               END IF
            END DO
         END DO
C
C        New trial vectors * New sigma vectors
C
         DO J = JCRED+1, N_CRED
            DO I = JCRED+1, N_CRED
               WORK(KHRED+(I-1)*MAX_HRED+(J-1)) =
     &              DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KSCVEC+(J-1)*NDETQ),1)
               IF (MZ .GT. 1) THEN
                  WORK(KHREDI+(I-1)*MAX_HRED+(J-1)) =
     &              DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KSCVEC+(J-1)*NDETQ+NDET),1)
     &            - DDOT(NDET,WORK(KBCVEC+(I-1)*NDETQ+NDET),1,
     &              WORK(KSCVEC+(J-1)*NDETQ),1)
               END IF
            END DO
         END DO
C
C        Check for inconsistencies:
C
chj      1.0D-10 can give problems with 1.0D-10 screening in GASCIP
chj      THRDAS = 1.0D-10
         THRDAS = THRSCR*10.0D0
         DASR = D0
         DASI = D0
         DO I = 1, N_CRED
            DO J = 1, I - 1
               DASR = DASR +
     &            ABS( WORK(KHRED + (J-1) + (I-1)*MAX_HRED) -
     &                 WORK(KHRED + (I-1) + (J-1)*MAX_HRED) )
               IF (MZ .GT. 1) DASI = DASI +
     &            ABS( WORK(KHREDI + (J-1) + (I-1)*MAX_HRED) +
     &                 WORK(KHREDI + (I-1) + (J-1)*MAX_HRED) )
            END DO
         END DO
         DAS = DASR + DASI
         IF (IPRINT .GE. 6 .OR. DAS .GT. THRDAS) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RCISTD) Reduced Hamiltonian, real part'
            CALL OUTPUT(WORK(KHRED),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,1,LUPRI)
           IF (MZ .GT. 1) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RCISTD) Reduced Hamiltonian, imaginary part'
            CALL OUTPUT(WORK(KHREDI),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,1,LUPRI)
           END IF
         END IF
         IF (IPRINT .GE. 5 .OR. DAS .GT. 0.01D0*THRDAS) THEN
            WRITE(LUPRI,'(A,1P,2D10.3)')
     &      ' (GASCIP_RCISTD) Anti-hermiticity of red. Hamiltonian = ',
     &      DASR, DASI
            IF (DAS .GT. THRDAS) THEN
               WRITE(LUPRI,'(/A,1P,D10.3)')
     &          ' FATAL ERROR, this is greater than threshold',THRDAS
               CALL QUIT('Anti-hermiticity of red. Ham. > threshold')
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
            CALL RSJACO(MAX_HRED,N_CRED,N_CRED,WORK(KHRED),WORK(KEIGVL),
     &        IJOB,IORDER,IPACK,WORK(KEIGVC),WORK(KFREE),LFREE)
         ELSE
            MATZ = 1
            CALL MEMGET('REAL',KHREDCP,LHRED,WORK,KFREE,LFREE)
            CALL DCOPY(LHRED,WORK(KHRED),1,WORK(KHREDCP),1)
            CALL QDIAG(2,N_CRED,WORK(KHREDCP),MAX_HRED,MAX_HRED,
     &        WORK(KEIGVL),MATZ,WORK(KEIGVC),MAX_HRED,MAX_HRED,
     &        WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/2A,I4)')
     &           ' *** ERROR in GASCIP_RCISTD ***: ',
     &           'QDIAG failed with error code ',IERR
               CALL QUIT('*** ERROR in GASCIP_RCISTD *** QDIAG failed')
            END IF
            CALL MEMREL('GASCIP_RCISTD.QDIAG',
     &         WORK,KWORK,KHREDCP,KFREE,LFREE)
         END IF
C
         IF (IPRINT .GE. 5) THEN
            WRITE(LUPRI,'(/A)')
     &      ' (GASCIP_RCISTD) Eigenvalues of reduced Hamiltonian'
            CALL OUTPUT(WORK(KEIGVL),1,1,1,N_CRED,1,MAX_HRED,1,LUPRI)
         END IF
         IF (IPRINT .GE. 6) THEN
            WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD)'//
     &      ' Eigenvectors of reduced Hamiltonian, real part'
            CALL OUTPUT(WORK(KEIGVC),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,1,LUPRI)
           IF (MZ .GT. 1) THEN
            WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD)'//
     &      ' Eigenvectors of reduced Hamiltonian, imaginary part'
            CALL OUTPUT(WORK(KEIGVCI),1,N_CRED,1,N_CRED,
     &                  MAX_HRED,MAX_HRED,1,LUPRI)
           END IF
         END IF
         ECIOLD = ECI
         ECIA = WORK(KEIGVL + ISTATE - 1)
         ECI = ECORE + ECIA
         IF (IPRINT .GE. 4) THEN
            WRITE(LUPRI,'(2(/A,F20.10))')
     &        '                Core energy        = ',ECORE,
     &        '                Active energy      = ',ECIA
         END IF
         WRITE(LUPRI,'(A,F20.10)')
     &        ' (GASCIP_RCISTD) Total energy       = ',ECI,
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
         CALL DZERO(CREF,NDETQ)
         DO I = 1, N_CRED
            FAC = WORK(KEIGVC+(I-1)+(ISTATE-1)*MAX_HRED)
            CALL DAXPY(NDETQ,FAC,WORK(KSCVEC+(I-1)*NDETQ),1,
     &           CREF,1)
            CALL DAXPY(NDETQ,-ECIA*FAC,WORK(KBCVEC+(I-1)*NDETQ),1,
     &           CREF,1)
         END DO
         DNRES = DNORM2(NDETQ,CREF,1)
         WRITE(LUPRI,'(A,1P,D18.3,A,D8.2,A)')
     &        '                Norm of residual   =',DNRES,
     &        '  (thr = ',RTHRCI,')'
         IF (IPRINT .GE. 10) THEN
            WRITE(LUPRI,'(/A)')
     &           ' (GASCIP_RCISTD) Residual'
            CALL OUTPUT(CREF,1,NDET,1,MZ,NDET,MZ,1,LUPRI)
            WRITE(LUPRI,'()')
         END IF
C
         CALL FLSHFO(LUPRI)
C
         IF (DNRES .LE. RTHRCI) THEN
            ICONV_TOT = 1
            GOTO 200
         END IF
         IF (ITMIC .GE. MAXCIT) THEN
            ICONV_TOT = -2
            GOTO 200
         END IF
         IF (N_CRED .GE. MAX_HRED) THEN
            ICONV_TOT = -3
            GOTO 200
         END IF
C
C        Construct 1 new trial vector
C        ----------------------------
C
         JCRED = N_CRED
         N_CRED = N_CRED + 1
         N_CSIM = 1
C
         DO IZ = 1, MZ
            DO I = 1, NDET
               HD = WORK(KHDIAG+I-1)
               IF (ABS(HD-ECIA) .GT. THRZER) THEN
                  WORK(KBCVEC+(I-1)+(N_CRED-1)*NDETQ + (IZ-1)*NDET) =
     &                 - CREF(I,IZ)/(HD - ECIA)
               ELSE
                  WORK(KBCVEC+(I-1)+(N_CRED-1)*NDETQ + (IZ-1)*NDET) =
     &                 - CREF(I,IZ)
               END IF
            END DO
         END DO
         CALL RCINORM(WORK(KBCVEC+(N_CRED-1)*NDETQ),IPRINT)
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
            IF (N_CSIM .GT. 1)
     &           CALL QUIT('*** ERROR in GASCIP_RCISTD ***'//
     &           ' N_CSIM .gt. 1 is not implemented')
C
C
C           Project out all previous trial vectors
C
            DNBEF = DNORM2(NDETQ,WORK(KBCVEC+(N_CRED-1)*NDETQ),1)
            DO I = 1, JCRED
C
C              FIXME: loop over new trial vectors
C
C              Calculate <B(new)|B(old)> / <B(old)|B(old)>
C              = <B(new)|B(old)>  (as trial vectors are normed to 1)
C
               OVLP = DDOT(NDETQ,WORK(KBCVEC+(I-1)*NDETQ),1,
     &              WORK(KBCVEC+(N_CRED-1)*NDETQ),1)
C
C              Calculate B(new) = B(new) - OVLP * B(old)
C
               IF (ABS(OVLP) .GT. THRZER) THEN
                  CALL DAXPY(NDETQ,-OVLP,WORK(KBCVEC+(I-1)*NDETQ),1,
     &                 WORK(KBCVEC+(N_CRED-1)*NDETQ),1)
               END IF
            END DO
C
#warning    FIXME: Orthogonalize new trial vectors against each other
C
            DNAFT = DNORM2(NDETQ,WORK(KBCVEC+(N_CRED-1)*NDETQ),1)
            IF (ABS(DNAFT) .LE. THRZER) THEN
               ICONV_TOT = -1
               GO TO 200
            END IF
            IF (DNAFT .LE. 1.0D-4 * DNBEF .AND. ITURN .EQ. 1)
     &           IROUND = IROUND + 1
C
            CALL DSCAL(NDETQ,1/DNAFT,WORK(KBCVEC+(N_CRED-1)*NDETQ),1)
C hj 16-jun-2004 TODO: screen on elements here ...
C
            IF (IROUND .GT. 0) THEN
               IROUND = 0
               IF (ITURN .EQ. 1) GO TO 400
               CALL QUIT('*** ERROR in GASCIP_RCISTD ***' //
     &              'IROUND .gt. 0 and ITURN .gt. 1')
            END IF
C
C           Print timings
C
            IF (IPRINT .GE. 2) THEN
               CALL GETTIM(CPU4,WALL4)
               CPUTID = SECTID(CPU4-CPU3)
               WALLTID = SECTID(WALL4-WALL3)
               WRITE(LUPRI,'(/A,I3,5A)')
     &              ' (GASCIP_RCISTD) CPU (Wall) time for iteration ',
     &              ITMIC,': ',CPUTID,' (',WALLTID,')'
            END IF
C
         GOTO 100
C
 200  CONTINUE
C
C     Print timings
C
      IF (IPRINT .GE. 2) THEN
         CALL GETTIM(CPU4,WALL4)
         CPUTID = SECTID(CPU4-CPU3)
         WALLTID = SECTID(WALL4-WALL3)
         WRITE(LUPRI,'(/A,I4,5A)')
     &        ' (GASCIP_RCISTD) CPU (Wall) time for iteration',
     &        ITMIC,': ',CPUTID,' (',WALLTID,')'
      END IF
      IF (IPRINT .GE. 0) THEN
         CPUTID = SECTID(CPU4-CPU2)
         WALLTID = SECTID(WALL4-WALL2)
         WRITE(LUPRI,'(/6A)')
     &        ' (GASCIP_RCISTD) Total CPU (Wall) '//
     &        'time for microiterations: ',
     &        CPUTID,' (',WALLTID,')'
      END IF
C
      IF (ICONV_TOT .GT. 0) THEN
         WRITE(LUPRI,'(//A)')
     &        ' (GASCIP_RCISTD) Micro iterations converged.'
      ELSE
         IF (ICONV_TOT .EQ. -1) THEN
            WRITE(LUPRI,'(//A)') ' (GASCIP_RCISTD) '//
     &           'WARNING: linear dependency in trial vectors.'
         ELSE IF (ICONV_TOT .EQ. -2) THEN
            WRITE(LUPRI,'(/A,I4,A)') ' (GASCIP_RCISTD) '//
     &           'WARNING: maximum number of micro iterations',
     &           MAXCIT, ' reached.'
         ELSE IF (ICONV_TOT .EQ. -3) THEN
            WRITE(LUPRI,'(//A,I4,A)') ' (GASCIP_RCISTD) '//
     &           'WARNING: maximum reduced space dimension of',
     &           MAX_HRED, ' reached.'
         ELSE
            WRITE(LUPRI,'(//A,I4/A)') ' (GASCIP_RCISTD) '//
     &        'WARNING: micro iterations aborted while not converged.',
     &        ' ICONV_TOT value :',ICONV_TOT
         END IF
      END IF
C
C
      MAXCIT = ITMIC
      ECI = ECORE + ECIA
      WRITE(LUPRI,'(3(/A,F20.10))')
     &     '                Core energy       = ',ECORE,
     &     '                Active energy     = ',ECIA,
     &     ' (GASCIP_RCISTD) Final CI energy   = ',ECI
C
C     Construct CREF:
C
      CALL DZERO(CREF,NDETQ)
      DO I = 1, N_CRED
         CALL DAXPY(NDETQ,WORK(KEIGVC+(I-1)+(ISTATE-1)*MAX_HRED),
     &        WORK(KBCVEC+(I-1)*NDETQ),1,CREF,1)
      END DO
C
C     Check norm of CI vector
C
      CALL RCINORM(CREF,10)
C
      IF ( IPRINT .GE. 5 ) THEN
         WRITE(LUPRI,'(/A)') ' (GASCIP_RCISTD) CI start vector'
         CALL RPRCI(CREF,NDET,NZ,CIPROG, THRPCI,LUPRI)
      END IF
C
      CALL GETTIM(CPU2,WALL2)
      IF (IPRINT .GE. 0) THEN
         CPUTID = SECTID(CPU2-CPU1)
         WALLTID = SECTID(WALL2-WALL1)
         WRITE(LUPRI,9000) CPUTID,WALLTID
      END IF
 9000 FORMAT(/' (GASCIP_RCISTD) Total CPU (Wall) time for routine: ',
     &     A,' (',A,')')
C
      CALL MEMREL('GASCIP_RCISTD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GASCIP_RCISTD')
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
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION IASTR(*), IBSTR(*)
C
      LOGICAL LBIT
C
#include "ibtfun.h"
C
#ifndef OLDCODE
      IELEC = 0
      CALL IZERO(IASTR,NASHT)
      CALL IZERO(IBSTR,NASHT)
      DO J = 1, NASHT
         IF (LBIT(IA,J) .OR. LBIT(IB,J)) THEN
            IELEC = IELEC + 1
            IF (LBIT(IA,J)) IASTR(IELEC) = J
            IF (LBIT(IB,J)) IBSTR(IELEC) = J
         END IF
      END DO
      IAELEC = IELEC
      IBELEC = IELEC
#else
      IAELEC = 0
      IBELEC = 0
C
      DO J = 1, NASHT
         IF (LBIT(IA,J)) THEN
            IAELEC = IAELEC + 1
            IASTR(IAELEC) = J
         END IF
         IF (LBIT(IB,J)) THEN
            IBELEC = IBELEC + 1
            IBSTR(IBELEC) = J
         END IF
      END DO
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_anaci */
      SUBROUTINE GASCIP_ANACI(CREF,NDET,IDET,THRPCI,WORK,LWORK)
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
CMI   PARAMETER (D1=1.0D0,D10 = 10.0D0)
CMI   ...D1 is defined in "consts.h"
      PARAMETER (D10 = 10.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgascip.h"
C
      DIMENSION CREF(NDET,*), IDET(2,*), WORK(*)
C
      DIMENSION IASTR(2*MAXASH_GASCIP), IBSTR(2*MAXASH_GASCIP)
C
      CALL QENTER('GASCIP_ANACI')
C
      MZ = MIN(2,NZ)
C
      DMX = D1
      DMN = D1 / D10
 10   CONTINUE
         WRITE(LUPRI,'(/A,1P,D7.1,A,D7.1/)')
     &     ' Print of CI coefficients in range ',DMN,' to ',DMX
         K = 0
         DO I = 1, NDET
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
               WRITE(LUPRI,'(5X,A,40I4)')
     &              'P-string     : ',(IASTR(J),J=1,IA)
               WRITE(LUPRI,'(5X,A,40I4/)')
     &              'P(bar)-string: ',(IBSTR(J),J=1,IB)
            END IF
         END DO
         IF (K .EQ. 0) WRITE(LUPRI,'(3X,A/)') 'None'
         DMX = DMN
         DMN = DMN / D10
         IF (DMX .GT. THRPCI) GO TO 10
C
      CALL QEXIT('GASCIP_ANACI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_chkkr */
      SUBROUTINE GASCIP_CHKKR(DOIMAG,BVEC,NDET,WORK,LWORK)
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
C
      LOGICAL   DOIMAG
      DIMENSION BVEC(*), WORK(*)
C
#include "memint.h"
C
      IF ( MOD(NAELEC,2) .NE. 0) RETURN
C
      CALL MEMGET('INTE',KRLIST,NDET,WORK,KFREE,LFREE)
      CALL IREAKRMC(LUKRMC,'KRLIST  ',WORK(KRLIST),NDET)
      CALL GASCIP_CHKKR1(DOIMAG,BVEC,NDET,WORK(KRLIST))
      CALL MEMREL('GASCIP_CHKKR',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_chkkr1 */
      SUBROUTINE GASCIP_CHKKR1(DOIMAG,BVEC,NDET,KRLIST)
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
#include "dgroup.h"
C
      LOGICAL   DOIMAG
      DIMENSION BVEC(NDET,*), KRLIST(*)
C
      DO I = 1, NDET
         IK = KRLIST(I)
         IF (IK .NE. I) THEN
            BVEC_DIFF = ABS(BVEC(I,1) - BVEC(IK,1))
ChjAaJ oct 08: first absolute check, if not OK, then relative check
C              (if not both, we might divide by zero!)
            IF ( BVEC_DIFF .GT. 1.0D-13 ) THEN
            IF ( BVEC_DIFF/ABS(BVEC(I,1)) .GT. 1.0D-13 ) THEN
               write(LUPRI,'(A,2I4,3D25.10)')
     &            ' GASCIP_CHKKR1: real diff',
     &            I,IK,BVEC(I,1),BVEC(IK,1),BVEC_DIFF
            END IF
            END IF
            IF (DOIMAG) THEN
               BVEC_DIFF = ABS(BVEC(I,2) + BVEC(IK,2))
ChjAaJ oct 08: first absolute check, if not OK, then relative check
C              (if not both, we might divide by zero!)
               IF ( BVEC_DIFF .GT. 1.0D-13 ) THEN
               IF ( BVEC_DIFF/ABS(BVEC(I,1)) .GT. 1.0D-13 ) THEN
                    write(LUPRI,'(A,2I4,3D25.10)')
     &                 ' GASCIP_CHKKR1: imag diff',
     &                 I,IK,BVEC(I,2),BVEC(IK,2),BVEC_DIFF
               END IF
               END IF
            END IF
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_genstr */
      SUBROUTINE GASCIP_GENSTR(IEL,ISTR,NSTR,MSKGAS,MAXDET)
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
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgascip.h"
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
      DO I = 1, NGAS_DC
C
         NGTOT = NGTOT + NGSHT(I)
         NGSPC(J) = MAX(0,NGSP(1,I) - NGTOT)
         NGSPC(J+1) = MIN(NGSP(2,I),NGTOT)
C        write(lupri,*) 'gascip_genstr: ngspc: ',ngspc(J),ngspc(j+1)
         J = J + 2
      END DO
C
      CALL CGENDET(IEL,NASHT,NSTR,ISTR,NGAS_DC,MSKGAS,NGSPC,MAXDET)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_ci */
      SUBROUTINE GASCIP_CI(DODIA,DIAG,DOSIG,SVEC,BVEC,HCORE,VMUUUU,
     &                    DODEN,CL,CR,DV,PV,DO2EL,NDET,IDET,IPRINT)
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
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION DIAG(*)
C
      DIMENSION SVEC(NDET,*), BVEC(NDET,*)
      DIMENSION VMUUUU(*), HCORE(*)
C
      DIMENSION CL(NDET,*), CR(NDET,*)
      DIMENSION DV(2*NASHT,2*NASHT,*)
      DIMENSION PV(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      LOGICAL   DODIA, DOSIG, DODEN, DO2EL, DODENCHK
C
      DIMENSION IDET(2,*)
C
      LOGICAL   DOIMAG, DOIMAG1
C
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      DIMENSION IBTAB(0:255)
      INTEGER*8 FNCALLS, CNCALLS
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
#include "ibtfun.h"
C
      CALL GETTIM(CPU1,WALL1)
C
C     Initialization
C     --------------
C
      NDETQ = NDET*NZ
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
      IF (DODIA) CALL DZERO(DIAG,NDET)
      IF (DOSIG) CALL DZERO(SVEC,NDETQ)
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
C        write(lupri,*) '(GASCIP_CI) MK2I = ',MK2I
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
C           (use that outer loop in gascip_gendet is over A-strings,
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
C                 ... outer loop in gascip_gendet is over A-strings
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
     &            CALL QUIT('GASCIP_CI: IELEC.ne.NAELEC')
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
C              write(lupri,*) '(GASCIP_CI) MK2I MK2J = ',MK2I,MK2J
C              call flshfo(lupri)
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
C                    ... outer loop in gascip_gendet is over A-strings
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
                     CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,
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
                     CALL GASCIP_PMAT(DO2EL,DOIMAG1,IT,IOCCI,JOCCI,
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
               CALL GASCIP_HMAT(DO2EL,.FALSE.,0,IOCCI,IOCCI,0,0,0,0,
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
               CALL GASCIP_PMAT(DO2EL,DOIMAG1,0,IOCCI,IOCCI,0,0,0,0,
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
C                 ... outer loop in gascip_gendet is over A-strings
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
                  CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,K,L,M,N,
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
                  CALL GASCIP_PMAT(DO2EL,DOIMAG1,IT,IOCCI,JOCCI,K,L,M,N,
     &                 DMRIJ,DMIIJ,DMRJI,DMIJI,DV,PV)
C
               END IF
 400        CONTINUE
 100     CONTINUE
C
 10   CONTINUE
C
      IF (IPRINT .GT. 0) THEN
        CALL GETTIM(CPU2,WALL2)
        IF ((WALL2-WALL1)*IPRINT .GT. 5.D0) THEN
C        ... do not print all the small numbers
          CPUTID  = SECTID(CPU2-CPU1)
          WALLTID = SECTID(WALL2-WALL1)
          WRITE(LUPRI,'(/5A)') ' CPU (Wall) time in GASCIP_CI kernel: ',
     &      CPUTID,' (',WALLTID,')'
          IF (NDSKIP .GT. 0) WRITE(LUPRI,'(A,I10,A,I10)')
     &      ' Number of configuration elements skipped:',NDSKIP,
     &      ' out of',NDET
          IF (IPRINT .GT. 2) THEN
            CALL GETNBITSCALLS(1,CNCALLS,FNCALLS)
            write(lupri,*) '#calls to nbits: ',cncalls,fncalls
          END IF
        END IF
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_cisig_screen */
      SUBROUTINE GASCIP_CISIG_SCREEN(SVEC,BVEC,HCORE,VMUUUU,
     &                              DO2EL,NDET,IDET,THRSCR,IPRINT)
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
#include "dcbgascip.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION SVEC(NDET,*), BVEC(NDET,*)
      DIMENSION VMUUUU(*), HCORE(*)
C
C
      LOGICAL   DO2EL
C
      DIMENSION IDET(2,*)
C
      LOGICAL   DOIMAG
C
      DIMENSION IOCCI(2*MAXASH_GASCIP), JOCCI(2*MAXASH_GASCIP)
      DIMENSION IBTAB(0:255)
      INTEGER*8 FNCALLS, CNCALLS
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
C
#include "ibtfun.h"
C
      CALL GETTIM(CPU1,WALL1)
C
C     Initialization
C     --------------
C
      NDETQ = NDET*NZ
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
      CALL DZERO(SVEC,NDETQ)
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
C        write(lupri,*) '(GASCIP_CISIG_SCREEN) MK2I = ',MK2I
         ISTRT = IMK2(MK2I) + 1
         IEND  = IMK2(MK2I) + NMK2(MK2I)
C
         MK2JSTRT = MAX(MINMK2,MK2I-ITMAX)
         MK2JEND  = MIN(MAXMK2,MK2I+ITMAX)
C
C        Loop over determinants
C        (use that outer loop in gascip_gendet is over A-strings,
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
C              ... outer loop in gascip_gendet is over A-strings
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
     &         CALL QUIT('GASCIP_CISIG_SCREEN: IELEC.ne.NAELEC')
C
C
C        *************************
C        *************************
C
C
            DO 20 MK2J = MK2JSTRT, MK2JEND, MK2STP
               IF (NMK2(MK2J) .EQ. 0) GOTO 20
C           write(lupri,*) '(GASCIP_CISIG_SCREEN) MK2I MK2J = ',MK2I,MK2J
C              call flshfo(lupri)
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
C                    ... outer loop in gascip_gendet is over A-strings
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
                  CALL GASCIP_HMAT(DO2EL,DOIMAG,IT,IOCCI,JOCCI,
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
      IF (IPRINT .GT. 0) THEN
         CALL GETTIM(CPU2,WALL2)
       IF ((WALL2-WALL1)*IPRINT .GT. 5.D0) THEN
C        ... do not print all the small numbers
         CPUTID  = SECTID(CPU2-CPU1)
         WALLTID = SECTID(WALL2-WALL1)
         WRITE(LUPRI,'(/5A/A,I10,A,I10)')
     &   ' CPU (Wall) time in GASCIP_CISIG_SCREEN kernel: ',CPUTID,' (',
     &   WALLTID,')',' Number of BVEC elements skipped:',NBSKIP,
     &   ' out of',NDET
        IF (IPRINT .GT. 2) THEN
         CALL GETNBITSCALLS(1,CNCALLS,FNCALLS)
         write(lupri,*) '#calls to nbits: ',cncalls,fncalls
        END IF
       END IF
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck gascip_inp */
      SUBROUTINE GASCIP_INP(WORD,RESET,WORK,LWORK)
!***********************************************************************
!
!     Input section for GASCIP-module
!
!     Written by: S. Yamamoto + tsaue, Feb 06, 2007
!
!     Copied from KRMCINP and modified.
!     Common block data will be modified in this routine.
!     The data is stored in the common block in the cossya.h file.
!     Called by: DIRCTL->PAMINP->PSIINP->GASCIP_INP
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
!
      PARAMETER (NTABLE = 31)
!
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK),IMKSET(3)
      DIMENSION MULTB_TMP(64,64)
!
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcborb.h"
#include "dcbgascip.h"
!
C.s/sya,2007.02.07
#include "dcbgos2.h"
C.q
!     ...The variable "SET" must have save attribute.
      SAVE SET
      DATA TABLE /'.PRINT ','.GASSH ','.GASSPC','.MK2REF','.MK2DEL',
     &            '.INACTI','.CASSCF','.CI SYM','.RELSYM','.MK2LIM',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX'/
      DATA SET/.FALSE./
!
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
C     Initialize DC_GASCIP_ORB
C     ================
C
      IPRGASCIP = 0
      IGASCIP_REP = -1
      DO IFRP = 1,NFSYM
        IGASCIP_NISH(IFRP) = 0
        IGASCIP_NASH(IFRP) = 0
      ENDDO
      IGASCIP_STATE = 1
      GASCIP_UCIBOS = .FALSE.
      CALL IZERO(IGASCIP_NGSHT,MAXASH_GASCIP)
C
C     Initialize DCBORB
C     =================
C     NOTE that we use the parts of it not used by HF !
C     We thus assume no MCSCF in the same run...
C
      IGASCIP_NGAS = 0
      DO I = 1,2
         IGASCIP_NISH(I) = 0
         IGASCIP_NASH(I) = 0
         DO J = 1, MAXASH_GASCIP 
            IGASCIP_NGSH(I,J) = 0
         END DO
      END DO
      IGASCIP_NAELEC = -1
C
C.s/sya,2007.02.07
C
C     Initialize DCBGOS2
C     ==================
      KDOGASCIP = .FALSE.
C.q
C
C     Process input
C     =============
C
      CALL IZERO(IMKSET,3)
      NEWDEF = (WORD .EQ. '*GASCIP')
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,
     &                     31), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in COSCI_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in COSCI_INP.')
!
    1          CONTINUE
!&&&  PRINT : print level in GASCIP module
                  READ(LUCMD,*) IPRGASCIP
               GO TO 100
!
    2          CONTINUE
!&&&  GASSH : GAS setup
                  IF (IGASCIP_NGAS.GT.0) THEN
                     CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You cannot specify both .GASSH and .CASSCF')
                  END IF
                  READ(LUCMD,*) IGASCIP_NGAS
                  IF (IGASCIP_NGAS .LT. 1 .OR.
     &                IGASCIP_NGAS .GT. MXMCGAS) THEN
                     CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                         'illegal value of IGASCIP_NGAS')
                  END IF
                  DO I = 1, IGASCIP_NGAS
                     READ(LUCMD,*) (IGASCIP_NGSH(J,I),J=1,NFSYM)
                  END DO
               GO TO 100
!
    3          CONTINUE
C&&&  GASSPC: GAS space constraints
                  IF (IGASCIP_NGAS .EQ. 0) THEN
                     WRITE(LUPRI,'(A)') ' *** ERROR in GASCIP_INP ***'
     &                  //' .GASSH must be specified before .GASSP!'
                     CALL QUIT('*** ERROR in GASCIP_INP ***')
                  END IF
                  DO I = 1, IGASCIP_NGAS
                     READ(LUCMD,*) (IGASCIP_NGSP(J,I),J=1,2)
                  END DO
                  IF (IGASCIP_NAELEC.LT.0) IGASCIP_NAELEC = 
     &               IGASCIP_NGSP(2,IGASCIP_NGAS)
               GO TO 100
!
    4          CONTINUE
C&&&  MK2REFERENCE: 2 * M_K reference
                  READ(LUCMD,*) IGASCIP_MK2REF
                  IMKSET(2) = 1
               GO TO 100
    5          CONTINUE
C&&&  MK2DELTA: 2 * DELTA M_K
                  READ(LUCMD,*) IGASCIP_MK2DEL
                  IMKSET(3) = 1
               GO TO 100
    6          CONTINUE
C&&&  INACTIVE ORBITALS
               READ(LUCMD,*) (IGASCIP_NISH(I),I=1,NFSYM)
               GO TO 100
    7          CONTINUE
C&&&  CASSCF -- this defines a CAS space
               IF (IGASCIP_NGAS .GT. 0) THEN
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You cannot specify both .GASSH and .CASSCF')
               END IF
               READ(LUCMD,*) IGASCIP_NAELEC
               READ(LUCMD,*) (IGASCIP_NASH(I),I=1,NFSYM)
               IGASCIP_NGAS = 1
               IGASCIP_NGSP(1,IGASCIP_NGAS) = IGASCIP_NAELEC
               IGASCIP_NGSP(2,IGASCIP_NGAS) = IGASCIP_NAELEC 
    8          CONTINUE
C&&&  SYMMETRY: symmetry of the wave function
               IF (IGASCIP_REP .GT. 0) THEN
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You can only specify symmetry once')
               END IF
               READ(LUCMD,*) IGASCIP_REP
               GO TO 100
!
    9          CONTINUE
C&&&  RELSYM: symmetry of the wave function
!.s/sya,2007.02.06
               IF (IGASCIP_REP .GT. 0) THEN
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                 'You can only specify symmetry once')
               END IF
               READ(LUCMD,*) IGASCIP_REP
!.q
               GO TO 100
!
   10          CONTINUE
C&&&  MK2LIM: min/max value of 2 * MK
               READ(LUCMD,*) IGASCIP_MINMK2,IGASCIP_MAXMK2
               IMKSET(1) = 1
               GO TO 100
!
   11          CONTINUE
   12          CONTINUE
   13          CONTINUE
   14          CONTINUE
   15          CONTINUE
   16          CONTINUE
   17          CONTINUE
   18          CONTINUE
   19          CONTINUE
   20          CONTINUE
   21          CONTINUE
   22          CONTINUE
   23          CONTINUE
   24          CONTINUE
   25          CONTINUE
   26          CONTINUE
   27          CONTINUE
   28          CONTINUE
   29          CONTINUE
   30          CONTINUE
   31          CONTINUE
!&&&  XXXXXX: invalid data
                  CALL QUIT('*** ERROR in GASCIP_INP *** ' //
     &                      'invalid data XXXXXX')
               GO TO 100
!
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in GASCIP_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in GASCIP_INP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process input
C     =============
C
      IF (IGASCIP_NGAS .GT. 0) THEN
         IGASCIP_NASHT = 0
         DO J = 1, NFSYM
            IGASCIP_NASH(J) = 0
            DO I = 1, IGASCIP_NGAS
               IGASCIP_NASH(J) = IGASCIP_NASH(J) + IGASCIP_NGSH(J,I)
            END DO
            IGASCIP_NASHT = IGASCIP_NASHT + IGASCIP_NASH(J)
         END DO
      END IF
      IF (IGASCIP_NASHT.GT.MAXASH_GASCIP) THEN
         WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &      ' *** ERROR in GASCIP_INP *** Too many active orbitals: ',
     &      IGASCIP_NASHT, ' -- max is ',MAXASH_GASCIP
         CALL QUIT(
     &      '*** ERROR in GASCIP_INP *** Too many active orbitals!')
      ENDIF
CTROND: set UCIBOS !!!!
C
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' KRMCINP: Set-up for GASCIP calculation'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I4)') ' * General print level   : ',IPRGASCIP
      WRITE(LUPRI,'(A,2I4)')
     &     ' * Inactive orbitals     : ',(IGASCIP_NISH(I),I=1,NFSYM)
      WRITE(LUPRI,'(A,2I4)')
     &     ' * Active orbitals       : ',(IGASCIP_NASH(I),I=1,NFSYM)
      WRITE(LUPRI,'(A,I4)')
     &     ' * Active electrons      : ',IGASCIP_NAELEC
      IF (IGASCIP_NGAS .GT. 0) THEN
         WRITE(LUPRI,'(A,I3,A)')
     &        ' * GAS space setup for ',IGASCIP_NGAS,' GAS space(s) : '
         DO I = 1, IGASCIP_NGAS
            WRITE(LUPRI,'(A,I3,A,2I4)')
     &           '   - GAS space ',I,'       : ',
     &           (IGASCIP_NGSH(J,I),J=1,NFSYM)
            WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &      '    (constraints: min/max active electrons after space : ',
     &           IGASCIP_NGSP(1,I),'/',IGASCIP_NGSP(2,I),')'
         END DO
         IF (IGASCIP_NGSP(2,IGASCIP_NGAS) .NE. IGASCIP_NAELEC)
     &        CALL QUIT('*** ERROR in GOCSCI_INP *** ' //
     &        'Inconsistency between .GASSP and .CASSCF')
      END IF
C
C     Symmetry of the wave function:
C     ------------------------------
C
C
C     MK constraints in the CI expansion.
C
      IF (IMKSET(1) .EQ. 1) THEN
C
C        MINMK2/MAXMK2 explicitly set
C
      ELSE
         IF (IMKSET(2) .EQ. 1 .AND. IMKSET(3) .EQ. 1) THEN
C
C           MK2REF/MK2DEL set, transform to MINMK2/MAXMK2
C
            IGASCIP_MINMK2 = IGASCIP_MK2REF - IGASCIP_MK2DEL
            IGASCIP_MAXMK2 = IGASCIP_MK2REF + IGASCIP_MK2DEL
C
         ELSE
C
C           Assign default values for MINMK2/MAXMK2
C
            IGASCIP_MAXMK2 = 2*MIN(IGASCIP_NAELEC,IGASCIP_NASHT)
     &                     - IGASCIP_NAELEC
            IGASCIP_MINMK2 = IGASCIP_NAELEC
     &                     - 2*MIN(IGASCIP_NAELEC,IGASCIP_NASHT)
            IGASCIP_MK2REF = 0
            IGASCIP_MK2DEL = IGASCIP_MAXMK2
         END IF
      END IF
C
C     Boson symmetry in the CI expansion
C
      IF(IGASCIP_REP.GT.0) THEN
        IF ((SPINFR.or.levyle).AND.GASCIP_UCIBOS) THEN
          NREP = NBSYM
          WRITE(LUPRI,'(A)')
     &    ' * Boson spatial spinor symmetry used in CI expansion.'
          WRITE(LUPRI,'(2A)')
     &    ' * Spatial symmetry (boson) of wave function : ',
     &    REP( IGASCIP_REP - 1 )
          WRITE(LUPRI,'(A,I4,A,I3)')
     &    ' * Allowed interval of 2 * MS :',
     &    IGASCIP_MINMK2, ' to ', IGASCIP_MAXMK2
        ELSE IF (.NOT.(SPINFR.or.levyle).AND.GASCIP_UCIBOS) THEN
          NREP = NBSYM
          WRITE(LUPRI,'(A)') ' * Approximate boson spatial spinor'//
     &    ' symmetry (as non-rel) used in CI expansion.'
          WRITE(LUPRI,'(2A)')
     &    ' * Approximate spatial symmetry (boson) of wave function : ',
     &    REP ( IGASCIP_REP - 1 )
          WRITE(LUPRI,'(A,I4,A,I3)')
     &    ' * Allowed interval of approximate 2 * MS :',
     &    IGASCIP_MINMK2, ' to ', IGASCIP_MAXMK2
        ELSE
          IF (SPINFR.or.levyle) THEN
            CALL GMULTSF(NREP,REPNA,MULTB_TMP)
          ELSE
            CALL GMULTA(NZ,NFSYM,NREP,REPNA,MULTB_TMP)
          ENDIF
          IF ( MOD ( IGASCIP_NAELEC, 2 ) .EQ. 0 ) THEN
            IRRP = NREP + IGASCIP_REP
            WRITE(LUPRI,'(2A)')
     &        ' * Symmetry (boson) of wave function: ',
     &        REPNA(IRRP)
          ELSE
            WRITE(LUPRI,'(2A)')
     &        ' * Symmetry (fermion) of wave function: ',
     &        REPNA(IGASCIP_REP)
          END IF
          WRITE(LUPRI,'(A,I4,A,I3)')
     &     ' * Allowed interval of 2 * MK :',
     &     IGASCIP_MINMK2, ' to ', IGASCIP_MAXMK2
        END IF
      ELSE
CTROND: imporive on this,...
        WRITE(LUPRI,'(A)') 
     &    '* All symmetries used'
      ENDIF
C
      IF (IGASCIP_REP .GT. NREP) THEN
            WRITE(LUPRI,'(/A//A,2I3)')
     &         ' INPUT ERROR, requested .SYMMETRY is non-existent!',
     &         ' symmetry, max value of symmetry = ', 
     &           IGASCIP_REP, NREP
            CALL QUIT('INPUT ERROR for *GASCIP  .SYMMETRY')
      END IF
C
C.s/sya,2007.02.07
      KDOGASCIP = .TRUE.
C.q
C
      CALL GASCIP_CPY
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gascip_cpy */
      SUBROUTINE GASCIP_CPY
C***********************************************************************
C
C     Transfer information from common block in dcbgascip.h to 
C     common block in dcbkrmc.h.
C
C     Written by: S. Yamamoto - Fri Feb  9 15:26:44 CET 2007
C     Last revision: 
C     Called by: GASCIP_INP
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxmcgas.h"
C
#include "dgroup.h"
#include "dcbgascip.h"
#include "dcbkrmc.h"
C
      CALL QENTER('GASCIP_CPY')
C
      NKRMCGAS = IGASCIP_NGAS
      NKRMC_MK2REF = IGASCIP_MK2REF
      NKRMC_MK2DEL = IGASCIP_MK2DEL
!radovan: IGASCIP_MINMK2 is not initialized
      NKRMC_MINMK2 = IGASCIP_MINMK2
!radovan: IGASCIP_MAXMK2 is not initialized
      NKRMC_MAXMK2 = IGASCIP_MAXMK2
ccc   CALL IZERO(NKRMC_GSHT,MXMCGAS)
      DO I = 1,NFSYM
         NKRMCISH(I) = IGASCIP_NISH(I)
         NKRMCASH(I) = IGASCIP_NASH(I)
         DO J = 1, NKRMCGAS
            NKRMCGSH(I,J) = IGASCIP_NGSH(I,J)
ccc         NGSHT(J) = NGSHT(J) + NGSH(I,J)
         END DO
ccc      NOCC(I) = IGASCIP_NISH(I) + IGASCIP_NASH(I)
      END DO
      NKRMCAELEC = IGASCIP_NAELEC
C     min and max electrons in each gas space:
      DO J = 1, NKRMCGAS
         NKRMCGSP(1,J) = IGASCIP_NGSP(1,J)
         NKRMCGSP(2,J) = IGASCIP_NGSP(2,J)
      END DO
C
!radovan: IPRGASCIP is not initialized
      IPRKRMC = IPRGASCIP
      IKRMC_ISTATE = 1
      IKRMC_SYMMETRY = IGASCIP_REP
C
      KRMC_UCIBOS = GASCIP_UCIBOS
      KRMC_CIPROGRAM = 'GASCIP'
C
      CALL QEXIT('GASCIP_CPY')
      RETURN
      END
