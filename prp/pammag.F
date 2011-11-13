!
!      Copyright (c) 2008 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
C
C "For Zion's sake will I not hold my peace, and for Jerusalem's sake I will not rest,
C  until her righteousness go forth as brightness, and her salvation as a lamp that burneth."
C                                       (Is 62:1)
C

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck NMRINP */
      SUBROUTINE NMRINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for magnetic properties:
C --------------------------------------------
C  NMR shielding, magnetizabilities with/ without London atomic orbitals ...
C 
C     Written by T.Enevoldsen - Nov 1997
C
C     Added: MI & HJAaJ, 2002/2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "orgcom.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbham.h"
#include "nuclei.h"
C

      PARAMETER (NTABLE = 15)
      LOGICAL SET, NEWDEF, RESET, GORG, LBIT, USECM, MODRDSUSLL
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &          PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.GAUGEO','.NOLOND','.SYMCON','.NOORTH',
     &            '.INTFLG','.NOTWO ','.NOONEI','.DOEPRN',
     &            '.NOEPRN','.LONDON','.USECM ','.MLLSUS',
     &            '.XXXXXX','.XXXXXX','.XXXXXX'   /
      DATA SET/.FALSE./
C
      CALL QENTER('NMRINP')
#include "memint.h"
C
C     IF (RESET) GOTO 999

      NEWDEF = (WORD .EQ. '*NMR')
      IF (SET) THEN
         IF (NEWDEF)
     &      CALL QUIT('NMRINP: Only one "*NMR" input section allowed')
C        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         GOTO 999
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C     Initialize 
C     ===================
C
CMI   LONDON = SHIELD .OR. MSUSCP .OR. OPTROT
CMI ... by default we use conventional orbitals to ease running DIRAC tests
      LONDON = .FALSE.
      SYMCON = .FALSE.
      NOORTH = .FALSE.
      NOTWOL = .FALSE.
      NOONEI = .FALSE.
      GORG   = .FALSE.
      USECM  = .FALSE.
      IF (LEVYLE) THEN
       RDSUSLLMOD = .TRUE.
      ELSE
       RDSUSLLMOD = .FALSE.
      ENDIF
      MODRDSUSLL = .FALSE.
C ... For the restricted kinetic balance, do not sum
C     reorthonormalization terms over positronic shells !
      DOEPRN = .FALSE.
C ... for the unrestricted kinetic balance sum by default
C     reorthonormalization terms over positronic shells
C ... when you want to change the DOEPRN default value in URKBAL,
C      use the NOEPRN keyword (see below)
      IF (URKBAL) DOEPRN=.TRUE.
CTROND: I can not see why these should be set to zero here....!!!
CTROND  ..they are initialized already in HERINP...
CTROND      CALL DZERO(GAGORG,3)
CTROND      CALL DZERO(DIPORG,3)

C INTNMR here in NMRINP is read before INTXLR (in XLRINP), 
C  so it has its own default value
C of all 3 integral classed (LL+LS+SS)
      INTNMR = 7  
C
C     Process input 
C     =========================
C
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in NMRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in NMRINP.')
    1          CONTINUE
C ... read user specified gauge origin
                  READ(LUCMD,*) (GAGORG(ICOOR), ICOOR = 1,3)
                  GORG = .TRUE.
               GO TO 100
    2          CONTINUE
C ... use conventional atomic orbitals instead of London atomic orbitals
               LONDON = .FALSE.
               GO TO 100
    3          CONTINUE
C ...  use symmetric connection for reorthonormalization terms
C ... default is to use the natural (non-symmetric) connection
               SYMCON = .TRUE.
               GO TO 100
    4          CONTINUE
C ... do not calculate reorthonormalization terms, {T(1),h(m)},
C     for the exp.value term h(Bm)
               NOORTH = .TRUE.
               GO TO 100
    5          CONTINUE
C INTFLG: specify what two-electron integrals to include
C         for the two-electron London contibution g(B)_mnpq
                  READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  INTNMR = ILLINT + 2*ISLINT + 4*ISSINT
               GO TO 100
    6          CONTINUE
C  ... do not calculate two-electron London contribution, g(B)_mnpq
               NOTWOL = .TRUE.
               GO TO 100
    7          CONTINUE
C ... do not calculate  reorthonormalization terms {T(1),H(0)}
C      for the property gradient H(B)
               NOONEI = .TRUE.
               GO TO 100
    8          CONTINUE
C  ...  do summations in reorthonormalization terms
C  over both electronic and positronic states
               DOEPRN = .TRUE.
               GO TO 100
    9          CONTINUE
C ... NOEPRN: deactivate the DOEPRN keyword which is by default set to true
C     in the URKBAL scheme; in RKBAL the default is .FALSE.
               DOEPRN = .FALSE.
               GO TO 100
   10          CONTINUE
               WRITE(LUPRI,'(1X,A)') '* This version of DIRAC does not'
     & //' support London atomic orbitals !'
     
               GO TO 100
   11          CONTINUE
C&&&& USECM: use center of mass as a default gauge origin
               USECM = .TRUE.
               GO TO 100
   12          CONTINUE
C&&&& MLSUS: Activate the modification of the relativistic RDSUSLL term
               RDSUSLLMOD = .TRUE.
               MODRDSUSLL = .TRUE.
               GO TO 100
   13          CONTINUE
               GO TO 100
   14          CONTINUE
               GO TO 100
   15          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *         '" not recognized in NMRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in NMRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
C ... header:
C     CALL PRSYMB(LUPRI,'=',75,0)
C     WRITE(LUPRI,'(A)')
C    &   ' NMRINP: Magnetic properties'
C     CALL PRSYMB(LUPRI,'=',75,0)

      IF (USECM.AND.GORG) THEN
       WRITE(LUPRI,'(A)')
     & 'You have asked for center of mass to be the gauge origin!'
       WRITE(LUPRI,'(A)')
     & '...your defined gauge origin is then ignored.'
       GORG = .FALSE.
      ENDIF

      IF (GORG) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A)') '*NMRINP: '
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(A,3F12.6,/)') 'User defined gauge '//
     &   'origin used:',(GAGORG(I), I = 1, 3)
         IF (.NOT.SYMCON.AND.LONDON) WRITE(LUPRI,'(1X,A,/)')
     &  ' WARNING: User defined gauge '//
     &  'may give wrong results when using symmetry !'   
CMI      .... Don't forget to copy GAGORG into the DIPORG !!!
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
      ELSE
       IF (LONDON.OR.USECM) THEN
C ... make the center of mass the gauge origin (default in dalton)
         DO I = 1, MXCENT
            ISOTOP(I) = 1
         END DO
         CALL MEMGET('REAL',KGEOM  ,3*NATOMS ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KMASS  ,  NATOMS ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KNAT   ,  NATOMS ,WORK,KFREE,LFREE)
C ... call the routine in abacus/hergp.F
         CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),ISOTOP,IPRDEF)
         CALL MEMREL('NMRINP.CMMASS',WORK,KWORK,KGEOM,KFREE,LFREE)
         CALL DCOPY(3,CMXYZ,1,GAGORG,1)
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
         WRITE(LUPRI,'(A,3F12.6)')
     &  'NMRINP: Center of mass used as the gauge & dipole '//
     &  'origin :',(GAGORG(I), I = 1, 3)
       END IF
      END IF

      IF (LONDON) THEN

       IF (NOONEI) THEN
         WRITE(LUPRI,'(/,1X,A)')
     &': No calculation of {T(1),H(0)} reorthonormalization'//
     &' terms for the H(B) property gradient !'
       ENDIF

       IF (NOORTH.AND. (SHIELD.OR.OPTROT)) THEN
         WRITE(LUPRI,'(/,1X,A)')
     &': No calculation of {T(1),H(m_K)} reorthonormalization'//
     &' terms for the H(B m_K) expectation value terms !'
       ENDIF

       IF ((LEVYLE.OR.BSS).AND.DOEPRN) THEN
        WRITE(LUPRI,'(2X,A)') ': For the LevyLe/BSS '//
     &'Hamiltonian the only choice for DOEPRN is .FALSE.!  '
         DOEPRN = .FALSE.
       ENDIF

       IF ((LEVYLE.OR.BSS).AND.(INTNMR.GT.1)) THEN
         INTNMR=1 ! Include only LL integrals...
         WRITE(LUPRI,'(2X,A)') ': For the LevyLe/BSS '//
     &'Hamiltonian include only LL integrals for the dg(1,2)/dB terms!'
       ENDIF

       IF (.NOT.DOEPRN) THEN
         WRITE(LUPRI,'(/,1X,A)')
     & ': Reorthonormalization terms only over '//
     & 'electronic shells (recommended default for RSKINB,LevyLe,BSS) !'
        IF (URKBAL) WRITE(LUPRI,'(/1X,A)')
     &'...though the URKBAL scheme is on, the DOEPRN is .false.!'
       ELSE
         WRITE(LUPRI,'(/,1X,A)')
     & ': Reorthonormalization terms over both '//
     & 'positronic and electronic shells '//
     &' (recommended as the default for URKBAL)! '
        IF (.NOT.URKBAL) WRITE(LUPRI,'(1X,A,/,1X,A)')
     & 'THIS IS WRONG FOR THE RESTRICTED KINETIC BALANCE SCHEME !',
     & '------  This option is only experimental ! --------'
       ENDIF

       IF (SYMCON) THEN
         WRITE(LUPRI,'(1X,A)')
     & ': Symmetric connection used for LAO-PDBS.'
       ELSE
         WRITE(LUPRI,'(1X,A)')
     & ': Natural connection used for LAO-PDBS.'
       ENDIF

       IF (NOTWOL) THEN
         WRITE(LUPRI,'(1X,A)')
     & ': No two-electron London contribution to magn.prop. gradient!'
       ELSE
         WRITE(LUPRI,'(1X,A)')
     & ': Include London magnet.prop.grad.contributions from '//
     & 'the following two-electron integrals:'
         IF(LBIT(INTNMR,1)) WRITE(LUPRI,'(6X,A)') '- LL-integrals'
         IF(LBIT(INTNMR,2)) WRITE(LUPRI,'(6X,A)') '- SL-integrals'
         IF(LBIT(INTNMR,3)) WRITE(LUPRI,'(6X,A)') '- SS-integrals'
       ENDIF

C =========================================================================
CMI  ... for the magnetic susceptibilities take care the RDSUSLLMOD flag
C =========================================================================
       IF (MSUSCP) THEN
        IF (LEVYLE.AND.MODRDSUSLL) THEN
         WRITE(LUPRI,'(1X,A)') 
     &  '...Though the Levy-Lebl mode the full form of the RDSUSLL '//
     &  'operator used to achieve the compatibility with rel. calcul.'
          RDSUSLLMOD = .FALSE.
        ELSE
          IF (RDSUSLLMOD) THEN
            WRITE(LUPRI,'(1X,A)') 
     &  '...SL block only for the RDSUSLL operator used to achieve'//
     &  ' the non-relativistic limit of the magnetizability !'
          ELSE
            WRITE(LUPRI,'(1X,A)') 
     &  '...Full relativistic form of the RDSUSLL operator used!'
          ENDIF
        ENDIF
       ENDIF
      ELSE 
C   .... defined conventional atomic orbitals (AO), not London AO
       IF (SHIELD.OR.MSUSCP.OR.OPTROT) THEN
        WRITE(LUPRI,'(1X,A)')
     &': Used conventional atomic orbitals, not London orbitals !'
        WRITE(LUPRI,'(1X,A)')
     &'WARNING: Calculated magnetic properties are dependent'//
     &' on the choice of the gauge origin!'   
         WRITE(LUPRI,'(A,3F12.6)')
     &  'NMRINP: The gauge & dipole '//
     &  'origin :',(GAGORG(I), I = 1, 3)
       END IF
      END IF

  999 CONTINUE
      CALL QEXIT('NMRINP')
      RETURN
      END



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getmdm */
      SUBROUTINE GETMDM(EPMIX,ISYMT,DMAT,CMO,IBEIG,TMAT,
     &                   WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C  PURPOSE:
C ----------
C  Generate SYMMETRIZED modified density matrixes which are used for
C the additional two-electron reorthonormalization term 
C when employing London orbitals.
C    
C     We define (based on Helgaker,Jorgensen,JCP 95(1991) 2595):
C   -----------------------------------------------------------------
C             DMAT_{ab} = Sum_{i,o}[ c_{ai} T_{oi} c_{bi}^* ]
C             a,b - AO-indices (L+S)
C             o   - secondary MO-index (electronic or both e+p shells)
C             i   - inactive MO-index (=occupied electronic MO's)
C
C     DMAT_{ab} is general matrix, so do symmetrize it 
C
C   On input: EPMIX=.true.=> o-summation runs over e+p shells, 
C             otherwise only over e-shells 
C             ISYMT - boson symmetry (1-NBSYM) of enetring T operator matrix
C             CMO,TMAT - MO coeff. + transformation matrix T
C
C  The T is expected to be a general matrix (natural connection), 
C  or an antihermitian (symmetric connection) matrix.
C
C  Written by Miro Ilias (based on TEC's routine TECTR2),
C    august 2003, Odense.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0 = 0.0D0,  D2 = 2.0D0)
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),IBEIG(*),
     &          TMAT(NORBT,NORBT,NZ),IQSYM(4),WORK(*)
      LOGICAL EPMIX

      CALL QENTER('GETMDM')
      KFRSAV = KFREE
C
C     Initialization
C     ==============
C
      IREP = ISYMT-1 ! irre[ pf T(1)_icomp (0 till nirrep-1)
      IOPSY = JBTOF(IREP,1)

      CALL DZERO(DMAT,N2BBASXQ)
C
C     Print section
C
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,'(/,A)')
     &  'GETMDM: Entering T(1) MO-connection matrix:'
        CALL PRQMAT(TMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &               IPQTOQ(1,IREP),LUPRI)

       IF(IPRINT.GE.10) THEN
        WRITE(LUPRI,'(/A)')
     &  'GETMDM: Entering MO coefficients:'
        CALL PRQMAT(CMO,NTBAS(0),NORBT,NTBAS(0),NORBT,NZ,
     &               IPQTOQ(1,0),LUPRI)
       ENDIF
      ENDIF

C
C     Generate MDM = modified density matrix of f.symm I1,I2
C   ===========================================================
C
      DO  I1 = 1, NFSYM
          I2 = MOD(I1+IOPSY,2) + 1
C
C  DMAT_{ab} = Sum_{i,o}[ c_{ao} T_{oi} c_{bi}^* ] ... C T(oi) C^+
C
        IF (EPMIX) THEN ! Run 'o' over both (p+e) shell
          IF (NOCC(I2).GT.0.AND.NORB(I1).GT.0) THEN
         CALL QTRANS(
     &    'MOAO','S',D0,NFBAS(I1,0),NFBAS(I2,0),
     &                  NORB(I1),NOCC(I2),
     &     DMAT(I2BASX(I1,I2)+1,1),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
!          T(o,i)
     &     TMAT(IORB(I1)+1,IORB(I2)+NPSH(I2)+1,1),
     &       NORBT,NORBT,NZ,IPQTOQ(1,IREP),
!          C(a,o)^+
     &     CMO(ICMOQ(I1)+1),
     &       NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
!          C(i,b)
     &     CMO(ICMOQ(I2)+(NPSH(I2)*NFBAS(I2,0))+1),
     &       NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)
         ENDIF
        ELSE ! Run o only over electronic states
        IF (NOCC(I2).GT.0.AND.NESH(I1).GT.0) THEN
      IF (IPRINT.GE.10) THEN
        CALL HEADER(
     &'GETMDM: Before QTRANS whole modified density matrix',-1)
        IF (EPMIX) THEN
         WRITE(LUPRI,'(2X,A)') "summation over both e+p shells..."
        ELSE
         WRITE(LUPRI,'(2X,A)') "summation over e shells..."
        ENDIF
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL HEADER(
     &'GETMDM: Before QTRANS BLOCK modified density matrix',-1)
        Write(LUPRI,'(2X,A,2I2)') "Ferm.symmetries I1,I2:",I1,I2
        CALL PRQMAT(DMAT(I2BASX(I1,I2)+1,1),NFBAS(I1,0),
     &              NFBAS(I2,0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
          ENDIF
         CALL QTRANS(
     &    'MOAO','S',D0,NFBAS(I1,0),NFBAS(I2,0),
     &       NESH(I1),NOCC(I2),
     &    DMAT(I2BASX(I1,I2)+1,1),
     &       NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
     &    TMAT(IORB(I1)+NPSH(I1)+1,IORB(I2)+NPSH(I2)+1,1),
     &       NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &    CMO(ICMOQ(I1)+(NPSH(I1)*NFBAS(I1,0))+1),
     &       NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &    CMO(ICMOQ(I2)+(NPSH(I2)*NFBAS(I2,0))+1),
     &       NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &    WORK(KFREE),LFREE,IPRINT)
      IF (IPRINT.GE.10) THEN
        CALL HEADER(
     &'GETMDM: Actual WHOLE modified density matrix',-1)
        IF (EPMIX) THEN
         WRITE(LUPRI,'(2X,A)') "summation over both e+p shells..."
        ELSE
         WRITE(LUPRI,'(2X,A)') "summation over e shells..."
        ENDIF
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         CALL HEADER(
     &'GETMDM: After QTRANS BLOCK modified density matrix',-1)
         Write(LUPRI,'(2X,A,2I2)') "Ferm.symmetries I1,I2:",I1,I2
         CALL PRQMAT(DMAT(I2BASX(I1,I2)+1,1),NFBAS(I1,0),
     &              NFBAS(I2,0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
          ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'GETMDM: NON-HERMITIAN modified density matrix',-1)
        IF (EPMIX) THEN
         WRITE(LUPRI,'(2X,A)') "summation over both e+p shells..."
        ELSE
         WRITE(LUPRI,'(2X,A)') "summation over e shells..."
        ENDIF
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
C
C     Take symmetric/antisymmetric combinations to form MDMAT
C    =========================================================
      IH = -1 ! T is TR antisymm.operator !
      DO IZ = 1,NZ
        IQ = IPQTOQ(IZ,IREP) ! Real quatern.number (1-4)
        IM = IHQMAT(IQ,IH) ! 1-S/2-A
        IQSYM(IZ) = IM
      ENDDO
C  ... in gp/gpsaue.F
      CALL QGETAM(DMAT,NTBAS(0),NZ,NTBAS(0),NTBAS(0),IQSYM)
      IF(IPRINT.GE.5) THEN
        CALL HEADER(
     &'GETMDM: Total SYMMETRIZED modified density matrix',-1)
        IF (EPMIX) THEN
         WRITE(LUPRI,'(2X,A)') "summation over both e+p shells..."
        ELSE
         WRITE(LUPRI,'(2X,A)') "summation over e shells..."
        ENDIF
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
      CALL QEXIT('GETMDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsklond */
      SUBROUTINE RSKLOND(FMAT,SKLTON,FACTOR,INDEX,IPOINT,
     &                   NDMAT,ISYMDM,IPRINT)
C*****************************************************************************
C
C   Purpose: Get two-electron Fock London matrixes in SA-AO basis
C
C   Written by Thomas Enevoldsen(TEC), 1998
C       and taken as a black-box by MI &HJAaJ, 2002 
C   Last modifications:  MI, March/2006 for pure twocomponent mode
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      INTEGER A, B, AP, BP, X
      DIMENSION FMAT(NBASIS,NBASIS,NDMAT,3), SKLTON(NBASIS,NBASIS),
     &          INDEX(KMAX,0:MAXREP,*), FACTOR(NBASIS,0:MAXREP),
     &          IPOINT(NBASIS,0:MAXREP),ISYMDM(NDMAT)
#include "nuclei.h"
#include "pincom.h"
#include "shells.h"
#include "symmet.h"
#include "ibtfun.h"
#include "sdpre.h"

      CALL QENTER('RSKLOND')
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER(
     & 'RSKLOND: Input skeleton Fock matrices ',1)
         DO 100 X = 1, 3
            DO N = 1, NDMAT
             WRITE(LUPRI,'(2X,A,I1,A,I1)')
     & '=== Component(X=1,Y=2,Z=3):',X,' quaternion(1-4):',N 
             CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,
     &            NBASIS,NBASIS,1,LUPRI)
            END DO
  100    CONTINUE
      END IF
C
C     Form anti/symmetric skeleton Fock matrix
C     ========================================
C ... in this case Fock matrixes are imaginary=antisymmetric

      SIGN = -D1
      DO N = 1, NDMAT
         DO 200 X = 1, 3
            DO 210 I = 1, NBASIS
            DO 210 J = 1, I
              FMATIJ = DP5*(FMAT(I,J,N,X)+SIGN*FMAT(J,I,N,X))
              FMAT(I,J,N,X) =      FMATIJ
              FMAT(J,I,N,X) = SIGN*FMATIJ
 210        CONTINUE
 200     CONTINUE
         IF (IPRINT .GT. 5) THEN
            CALL HEADER(
     &  'RSKLOND: Skeleton Fock matrixes',1)
            DO 300 X = 1, 3
               WRITE(LUPRI,'(/,1X,A,I1,A,I1,A,I1)')
     &'N(1-NDMAT)/NDMAT=',N,'/',NDMAT,'  of component(x-1/y-2/z-3)=',X 
               CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,NBASIS,NBASIS
     &              ,1,LUPRI)
 300        CONTINUE
         END IF
         SIGN = D1
      END DO ! End of N

C===============================================================
C
C     Form AO Fock matrix from symmetric skeleton matrix
C     ==================================================
C
C===============================================================
      IF (MAXREP .GT. 0) THEN
C
C        First construct index and symmetry factor arrays
C
         IORB = 0
         DO 400 ISHELL = 1, KMAX
            DO 410 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)).EQ.0) THEN
               DO 420 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  INDEX(ISHELL,ISYMOP,ICOMP) = IORB
  420          CONTINUE
            ELSE
               ISYMF = IBTXOR(ISYMOP,IBTAND(ISYMOP,ISTBAO(ISHELL)))
               DO 430 ICOMP  = 1, KHKT(ISHELL)
                  INDEX(ISHELL,ISYMOP,ICOMP) = INDEX(ISHELL,ISYMF,ICOMP)
  430          CONTINUE
            END IF
  410       CONTINUE
  400    CONTINUE
C
         IORB = 0
         DO 500 ISHELL = 1, NLRGSH
            NHKTA = NHKT(ISHELL)
            DO 510 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 520 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 530 ISYM = 0, MAXREP
                     ISYMPR = IBTXOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) =
     &                        PT(IBTAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  530             CONTINUE
  520          CONTINUE
            END IF
  510       CONTINUE
  500    CONTINUE

C
C        The next paragraph concerns relativistic calculations only !
C        IPARity is invoked since the operation of inversion 
C        receives a minus sign when acting on small components....
C
         DO ISHELL = NLRGSH+1,KMAX
            NHKTA = NHKT(ISHELL)
            DO ISYMOP = 0, MAXREP
               IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
                  DO ICOMP  = 1, KHKT(ISHELL)
                     IORB = IORB + 1
                     DO ISYM = 0, MAXREP
                        ISYMPR = IBTXOR(ISYMOP,ISYM)
                        IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                        FACTOR(IORB,ISYM) = IPAR(ISYM)*
     &                       PT(IBTAND(ISYM,ISYMAO(NHKTA,ICOMP)))
                     END DO
                  END DO
               END IF
            END DO
         END DO
C
C        END of relativistic code .....
C
C
C        Construct full matrix from skeleton
C
         DO N = 1, NDMAT
            DO  X = 1, 3
               CALL DCOPY(NBASIS*NBASIS,FMAT(1,1,N,X),1,SKLTON,1)
               DO 600 ISYMOP = 1, MAXREP
                  FP = PT(IBTAND(ISYMOP,IBTXOR(ISYMDM(N),ISYMAX(X,2))))
                  DO 620 A = 1, NBASIS
                     AP = IPOINT(A,ISYMOP)
                     FA = FACTOR(A,ISYMOP)
                     DO 630 B = 1, NBASIS
                        BP = IPOINT(B,ISYMOP)
                        FB = FACTOR(B,ISYMOP)
                        FMAT(AP,BP,N,X)=FMAT(AP,BP,N,X)
     &                       +FA*FB*FP*SKLTON(A,B)
 630                 CONTINUE
 620              CONTINUE
 600           CONTINUE
            END DO
         END DO
C
         FAC = D1/SDPRE(MAXREP + 1)
         NN = NDMAT*3*NBASIS*NBASIS
         CALL DSCAL(NN,FAC,FMAT,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("AO Fock matrices in RSKLOND",1)
            DO N = 1, NDMAT
               DO 700 X = 1, 3
                WRITE(LUPRI,'(/,1X,A,I1,A,I1,A,I1)')
     &'N(1-NDMAT)/NDMAT=',N,'/',NDMAT,'  of component(x-1/y-2/z-3)=',X 
                  CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,
     &                 NBASIS,NBASIS,1,LUPRI)
 700           CONTINUE
            END DO
         END IF
C
C      Check for relativistic calculation....
C
C     IF(KMAX.GT.NLRGSH) RETURN
CMI   IF(KMAX.GT.NLRGSH) GOTO 134
CMI   IF(KMAX.GE.NLRGSH) GOTO 134
      IF(.NOT.(KMAX.GE.NLRGSH)) THEN
CMI  ... seems that this body is never called....
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(3*(NBASIS**2),FMAT,1,SKLTON,1)
         CALL DZERO(FMAT,3*NBASIS*NBASIS)
         ISTRA = 1
         DO 800 IREPA = 0, MAXREP
            NORBA = NAOS(IREPA+1)
            DO 810 I = ISTRA,ISTRA + NORBA - 1
               IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
               NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
CTROND               IA    = IBTSHR(IPIND(I),8)
CTROND               NA    = IBTSHR(IPIND(I),4) - IA*16
               NHKTA = NHKT(IA)
               KHKTA = KHKT(IA)
               MULA  = ISTBAO(IA)
               INDA  = KSTRT(IA) + NA - KHKTA
               IVARA = IBTXOR(IREPA,ISYMAO(NHKTA,NA))
               DO 820 ISYMA = 0, MAXOPR
               IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
                  FA    = PT(IBTAND(ISYMA,IVARA))
                  INDA  = INDA + KHKTA
                  ISTRB = 1
                  DO 830 IREPB = 0, MAXREP
                     NORBB = NAOS(IREPB+1)
                     DO 840 J = ISTRB,ISTRB + NORBB - 1
                        IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                        NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
CTROND                        IB    = IBTSHR(IPIND(J),8)
CTROND                        NB    = IBTSHR(IPIND(J),4) - IB*16
                        NHKTB = NHKT(IB)
                        KHKTB = KHKT(IB)
                        MULB  = ISTBAO(IB)
                        INDB  = KSTRT(IB) + NB - KHKTB
                        IVARB = IBTXOR(IREPB,ISYMAO(NHKTB,NB))
                        DO 850 ISYMB = 0, MAXOPR
                        IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                           INDB = INDB + KHKTB
                           FAB  = FA*PT(IBTAND(ISYMB,IVARB))
                           DO 860 X = 1, 3
                           IF (IBTXOR(IREPA,IREPB).EQ.ISYMAX(X,2)) THEN
                              FMAT(I,J,N,X) = FMAT(I,J,N,X)
     &                                    + FAB*SKLTON(INDA,INDB)
                           END IF
  860                      CONTINUE
                        END IF
  850                   CONTINUE
  840                CONTINUE
                     ISTRB = ISTRB + NORBB
  830             CONTINUE
               END IF
  820          CONTINUE
  810       CONTINUE
            ISTRA = ISTRA + NORBA
  800    CONTINUE

         IF (IPRINT .GT. 5) THEN
            CALL HEADER("SO Fock matrix in RSKLOND",1)
            DO 900 X = 1, 3
               CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,
     &                       NBASIS,NBASIS,1,LUPRI)
  900       CONTINUE
         END IF
      END IF

      ENDIF

C134  CONTINUE ! Needed this statement...

      CALL QEXIT('RSKLOND')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qmphase */
      SUBROUTINE QMPHASE(TYP,IM,JSYMOP,ONEINT,DOINT)
C*****************************************************************************
C
C     Quaternion phase insertion in lower triangular matrix
C       Phase: (e_K1*)(e_IM)(e_K2)
C
C     Theory: T.Saue &  H.J.Aa.Jensen
C     Written by T.Saue October 1994
C     Last revision: May 27 1998 - tec
C                    Added TYP 'S' for Square matrices, 'L' for lower triangle
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "symmet.h"
C
#include "dgroup.h"
#include "dcbqal.h"
#include "dcbbas.h"
      LOGICAL DOINT(2,2),DIAG,SAME
      DIMENSION ONEINT(*)
      EXTERNAL BDQALGEBRA
      CHARACTER TYP*1
C
#include "ibtfun.h"
C
      IL = 2
      IR = 1
      I1 = 0
      IF (TYP.EQ.'L') THEN
         DO 10 IRP1 = 0,NBSYM-1
            DO 20 IC1 = 1,2
               NB1 = NBBAS(IRP1,IC1)
               K1  = JQBAS(IRP1,IC1)
               I2 = 0
               DO 30 IRP2 = 0,IRP1
                  SAME = IRP1.EQ.IRP2
                  ICMX = 2
                  IF(SAME) ICMX = IC1
                  DO 40 IC2 = 1,ICMX
                     NB2 = NBBAS(IRP2,IC2)
                     K2  = JQBAS(IRP2,IC2)
                     DIAG = SAME.AND.(IC1.EQ.IC2)
                     IF(DOINT(IC1,IC2).AND.(IBTXOR(IRP1,IRP2).EQ.JSYMOP)
     &                    ) THEN
                        KFAC = IQPHASE(K1,K2,IM)*IQPH(K1,IL)*IQPH(K2,IR)
                        IF(KFAC.EQ.-1) THEN
                           DO 50 J1 = 1,NB1
                              JMX = NB2
                              IF(DIAG) JMX = J1
                              JOFF = (I1+J1)*(I1+J1-1)/2 + I2
                              DO 60 J2 = 1,JMX
                                 ONEINT(JOFF+J2) = -ONEINT(JOFF+J2)
 60                           CONTINUE
 50                        CONTINUE
                        ENDIF
                     ENDIF
                     I2 = I2 + NB2
 40               CONTINUE
 30            CONTINUE
               I1   = I1 + NB1
 20         CONTINUE
 10      CONTINUE
      ELSE IF (TYP.EQ.'S') THEN
         DO IRP1 = 0,NBSYM-1
            ISY1 = IRP1 + 1
            DO IC1 = 1,2
               NB1 = NBBAS(IRP1,IC1)
               K1  = JQBAS(IRP1,IC1)
               I2 = 0
               DO IRP2 = 0,NBSYM-1
                  ISY2 = IRP2 + 1
                  DO IC2 = 1,2
                     NB2 = NBBAS(IRP2,IC2)
                     K2  = JQBAS(IRP2,IC2)
                     IF(DOINT(IC1,IC2).AND.
     &                    (IBTXOR(IRP1,IRP2).EQ.JSYMOP)) THEN
                        KFAC = IQPHASE(K1,K2,IM)*IQPH(K1,IL)*IQPH(K2,IR)
                        IF(KFAC.EQ.-1) THEN
                           IOFF = I2COSX(ISY1,ISY2,IC1,IC2)
                           DO ICOL = 1,NB2
                              CALL DSCAL(NB1,DM1,ONEINT(IOFF+1),1)
                              IOFF = IOFF + NTBAS(0)
                           END DO
                        ENDIF
                     ENDIF
                     I2 = I2 + NB2
                  END DO
               END DO
               I1   = I1 + NB1
            END DO
         END DO
      ELSE
         CALL QUIT('QMPHASE: Unknown TYP !')
      ENDIF
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtcnvd */
      SUBROUTINE GTCNVD(DMAT,CMO,IBEIG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/06/09 by J. Thyssen
C
C     Get contravariant denstity matrix without folding
C     Copied from GTCNV1
C     
C     On input: CMO - MO coeficcients
C               IBEIG - info on boson irreps      
C     
C     Last revisions: 1997/06/10 jth
C                     1998  TEC
C                     jan.2003 MI 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WORK(*),DMAT(N2BBASX,4,NFMAT),CMO(*),IBEIG(*)
      LOGICAL FOLD
      PARAMETER ( D1 = 1.0 D00, D2 = 2.0 D00 )
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbqal.h"
#include "dummy.h"
C
#include "ibtfun.h"
      IQOFF(IZ,IMAT) = (IZ-1) * N2BBASX + (IMAT-1) * N2BBASXQ
C
#include "memint.h"
C
C
      CALL QENTER('GTCNVD')
C 
      CALL MEMGET('REAL',KTMP,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Contruct density matrix
C
      DO IMAT = 1, NFMAT ! Run over closed/open-shells
         IBIT = IBTSHL(1,IMAT-1)
C      ... routine in ../dirac/dirden.F
         CALL GENDEN(WORK(KTMP),CMO,IBIT,IPRINT)
C
C        Scale inactive density matrix with D2 (= occupation)
C
         IF (IMAT.EQ.1)
     &     CALL DSCAL(N2BBASXQ,D2,WORK(KTMP),1)
C
C        Transform to unsorted basis
C
         IREP = ISYMOP(IMAT)-1
         CALL BSTOBU(WORK(KTMP),NZ,WORK(KFREE),LFREE)
         IF(IPRINT.GE.5) THEN
           CALL HEADER('GTCNVD: '//
     &     'DMAT After BSTOBU to unsorted basis',-1)
           WRITE(LUPRI,'(A,I3)') '*** Matrix no.: ',IMAT
           CALL PRQMAT(WORK(KTMP),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         ENDIF
C
C        Insert half-phases when NZ < 4 (real and complex groups)
C     =============================================================
C
         IREP = ISYMOP(IMAT)-1
         IF(NZ.LT.4) THEN
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREP)
              CALL Q2PHASE('D',IQ,1,WORK(KTMP+N2BBASX*(IZ-1)))
            ENDDO
            IF(IPRINT.GE.5) THEN
              CALL HEADER('GTCNVD: '//
     &  'Unsorted DMAT with half-phases inserted',-1)
              CALL PRQMAT(WORK(KTMP),NTBAS(0),NTBAS(0),NTBAS(0),
     &             NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            ENDIF
         ENDIF
C
C        Transform density matrix from QO to AO
C        ======================================
C   .... analogous to sr. DIRSCF...
C
         CALL DZERO(DMAT(1,1,IMAT),4*N2BBASX)
         IREP = ISYMOP(IMAT)-1
         IPAR = JBTOF(IREP,1)
         DO IZ = 1,4
           IREPD = IRQMAT(IZ,IREP)
           IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
           IPQ   = IQTOPQ(IQ,IREP)
           IF(IPQ.GT.NZ) CALL QUIT('GTCNVD.DENSITY: IPQ.GT.NZ !')
           IF(IPQ.LE.0 ) CALL QUIT('GTCNVD.DENSITY: IPQ.LE.0 !')
           IADR = (IPQ-1) * N2BBASX

           CALL DTSOAO(WORK(KTMP+IADR),DMAT(1,IZ,IMAT),
     &                 NTBAS(0),IREPD,IPRINT)
         ENDDO
      ENDDO

      IF (NASHT.GT.0) THEN ! When open-shell are included
         CALL DAXPY(N2BBASX*4,D1,DMAT(1,1,2),1,DMAT(1,1,1),1)
      END IF
C
C        Final print section
C
      IF (IPRINT .GE. 4) THEN
       CALL HEADER('GTCNVD: Final quaternion density matrices.',-1)
       DO IMAT = 1,NFMAT
         WRITE(LUPRI,'(2X,A,I1,A,I1,A)')
     &    '**** >>>>> Matrix number(1-cs,2-os1,3-os2...): ',
     &    IMAT,'/',NFMAT,' <<<<<< ****'
         DO IZ = 1, 4
          WRITE(LUPRI,'(/2X,A,I1)')
     &    '=====  quaternion number (1-4): ',IZ
          CALL OUTPUT(DMAT(1,IZ,IMAT),1,NTBAS(0),1,NTBAS(0),
     &                  NTBAS(0),NTBAS(0),1,LUPRI) 
         END DO
       END DO
      END IF

C     Dispose used memory
      CALL MEMREL('GTCNVD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GTCNVD')
      RETURN
      END


      INTEGER FUNCTION ICONVXYZ(C)
C******************************************************
C 
C   Integer function returning 1,2 or 3,
C  depending on input(x,y,z)
C 
C  Written by M.Ilias, Odense, august 2003
C
C******************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*1 C

      CALL QENTER('ICONVXYZ')

      IF (C.EQ.'X'.OR.C.EQ.'x') THEN
        ICONVXYZ=1
      ELSE IF (C.EQ.'Y'.OR.C.EQ.'y') THEN
        ICONVXYZ=2
      ELSE IF (C.EQ.'Z'.OR.C.EQ.'z') THEN
        ICONVXYZ=3
      ELSE
        WRITE(LUPRI,'(A,A)') 'ICONVXYZ: The input character:',C
        WRITE(LUPRI,'(A)') 'ICONVXYZ: and it has to be x/y/z...!!!'
        CALL QUIT('ICONVXYZ: Not x/y/z input!')
      ENDIF

C       write(lupri,*) 'ICONVXYZ: C=',C,' ICONVXYZ=',ICONVXYZ

      CALL QEXIT('ICONVXYZ')
      RETURN
      END
