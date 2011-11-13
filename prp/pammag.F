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
      use codata
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
     &            '.GO_ANG','.XXXXXX','.XXXXXX'   /
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

C INTNMR here in NMRINP is read before INTXLR (in XLRINP), 
C  so it has its own default value
C of all 3 integral classed (LL+LS+SS)
!miro TODO: take the integral default from HAMINP
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

!              .GAUGEO
    1          CONTINUE
C ... read user specified gauge origin
                  READ(LUCMD,*,IOSTAT=IOS) (GAGORG(ICOOR),ICOOR=1,3)
                  IF (IOS.NE.0) THEN
                   WRITE(LUPRI,*) '(GAGORG(ICOOR),ICOOR=1,3:',
     &             (GAGORG(ICOOR),ICOOR=1,3)
                   CALL QUIT(
     & 'NMRINP: error in reading (GAGORG(ICOOR),ICOOR=1,3) !')
                  ENDIF
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
                 READ(LUCMD,*,IOSTAT=IOS) ILLINT,ISLINT,ISSINT
                 IF (IOS.NE.0) THEN
                  WRITE(LUPRI,*)
     &  'ILLINT,ISLINT,ISSINT:',ILLINT,ISLINT,ISSINT
                  CALL QUIT(
     & 'NMRINP: Error in reading ILLINT,ISLINT,ISSINT!')
                 ENDIF
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
#ifdef MOD_LONDON
C&&&& LONDON: use London atomic orbitals rather than conventional atomic orbitals
               LONDON = .TRUE.
#else
               WRITE(LUPRI,'(1X, A)') '* This version of DIRAC does not'
     & //' support London atomic orbitals !'
     
#endif
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

!              .GO_ANG
   13          continue
!              read user specified gauge origin in angstrom
               read(lucmd, *, iostat=ios) gagorg
               if (ios /= 0) then
                 call quit('error in reading .GO_ANG')
               end if
!              convert from angstrom to bohr
               gagorg = gagorg/xtang
               gorg = .true.
               go to 100

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
         WRITE(LUPRI,'(3X,A)')
     &   'NMRINP: Input for London magnetic field gradient'//
     &   '(NMR shieldings, magnetizabilities).'
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(3X,A,3F12.6)') 'User defined gauge '//
     &   'origin used (a.u.):',(GAGORG(I), I = 1, 3)
         IF (.NOT.SYMCON.AND.LONDON) WRITE(LUPRI,'(3X,A,/)')
     &  'WARNING: User defined gauge '//
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
         WRITE(LUPRI,'(1X,A,3F12.6)')
     &  '* Center of mass used for the gauge & dipole '//
     &  'origin (a.u.) :',(GAGORG(I), I = 1, 3)
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

       IF ((LEVYLE.OR.BSS.or.x2c).AND.DOEPRN) THEN
        WRITE(LUPRI,'(2X,A)') ': For the LevyLe/BSS '//
     &'Hamiltonian the only choice for DOEPRN is .FALSE.!  '
         DOEPRN = .FALSE.
       ENDIF

       IF ((LEVYLE.OR.BSS.or.x2c).AND.(INTNMR.GT.1)) THEN
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
     &'  Used conventional atomic orbitals, not London orbitals !'
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

#ifdef MOD_LONDON
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETCON */
      SUBROUTINE GETCON(W,CONMAT,DMAT,ORTO,WORK,
     &                  KFREE,LFREE,IPRINT)
C***********************************************************************
C
C   Get the reortonormalization contribution term,
C   {T(1),h(m_K)}, when using LAO. Calculated as expectation value.
C
C   MI/sept.2003: The .DOEPRN parameter decides whether to do
C   reorthonormalization over both p+e shells, or only over
C   the electronic shells. 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcborb.h"   
#include "dcbbas.h"   
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbxlr.h"
#include "dummy.h"
      PARAMETER(D2 = 2.0D0,D1 = 1.0D0,D0 = 0.0D0)
      PARAMETER(DM1 = -1.0D0,DM2 = - 2.0D0)
!               Array for {T,h(mK)}_ic product
      DIMENSION W(NORBT,NORBT,NZ,3),
!              Array for h(mK) property matrix
     &         CONMAT(NORBT,NORBT,NZ),
!              Density matrix
     &         DMAT(NORBT,NORBT,NZ),
     &         ORTO(3,*),WORK(*)
      LOGICAL AH,EPMIX,DOPCTRA_SAVE
C
#include "ibtfun.h"
      CALL QENTER('GETCON')

      KFRSAV = KFREE

C  ...array for MO coefficients
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C boson irreps identification (read with MO coeff.)
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GETCON')
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident.
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,WORK(KBEIG),
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,IDUM,
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ENDIF

C Open property file
      CALL OPNFIL(LU1INT,'AOPROPER','OLD','GETCON')

C  ... arrays for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KTMAT,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTHMAT,N2ORBXQ,WORK,KFREE,LFREE)

      IF (IPRINT.GE.5) CALL HEADER(
     &'GETCON: Creating transformed matrix; {T_c, h_mK}',-1)
C --------------------------------------------
C -- Calculate expect. value of {T(1),h(m)} --
C --------------------------------------------
      NSCOOR = 3*NUCDEP
      NOP = 0
      DO N = 1, NSCOOR  ! Run over coordinates m_K
         INDOP = IPSHIELD(N)
         IF (INDOP.GT.0) THEN ! Get h(m_K) in MO basis into CONMAT
           NOP = NOP + 1
           INDXPR = LLRAPU(INDOP) ! Index of h(m_K)
           ISYM  = IPRPSYM(INDXPR)
           IREP = ISYM - 1
           IOPSY = JBTOF(ISYM-1,1)
           CALL DZERO(CONMAT,N2ORBXQ)
           CALL PRPMAT(INDXPR,IOPSY,CONMAT,WORK(KCMO),WORK(KBEIG),
     &                  ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
         
           IF (IPRINT.GE.5) THEN
            WRITE(LUPRI,'(2A)')
     &      'GETCON: Get MO matrix of property:', PRPNAM(INDXPR)
            CALL PRQMAT(CONMAT,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
           ENDIF
C --------------------------------------------------------------------
C ...  calculate the {T(1),h(m_K)}_ic reorthonormalization terms 
C    as expectation value
C --------------------------------------------------------------------
C .. do not perform the picture change transformation 
C    of connection matrixes !
           DOPCTRA_SAVE=DOPCT
           DOPCT=.FALSE.

           DO ICOMP  = 1, 3 ! Run over x,y,z components
C .... get the connection matrix (symmetric or natural) into WORK(KTMAT)
             INDXPR2 = IPCON(ICOMP) ! Get index of S(1)/W(1)
             ISYM2   = IPRPSYM(INDXPR2)
             ITIM2   = IPRPTIM(INDXPR2)
             IOPSY2  = JBTOF(ISYM2-1,1)
             IREP2   = ISYM2 - 1 ! 0-7
             CALL DZERO(WORK(KTMAT),N2ORBXQ)
             IF (.NOT.SYMCON) THEN
              IF (ITIM2.NE.0)
     &         CALL QUIT('GETCON: Wrong TR symm. of W(1)!')
              CALL PRPMAT(INDXPR2,IOPSY2,WORK(KTMAT),WORK(KCMO),
     &             WORK(KBEIG),ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
             ELSE
              IF (ITIM2.NE.-1)
     &         CALL QUIT('GETCON: Wrong TR symm. of S(1)!')

              CALL PRPMAT(INDXPR2,IOPSY2,WORK(KTMAT),WORK(KCMO),
     &             WORK(KBEIG),ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
             ENDIF

             IREP12 = IBTXOR(IREP2,IREP) ! resulting irrep of T_ic x m_K
             CALL DZERO(WORK(KTHMAT),N2ORBXQ)
             EPMIX = DOEPRN
             AH  = .TRUE. ! T_icomp is antihermitian ! 
C ... do the reorthonormalization ...
             CALL TTRA(EPMIX,AH,ISYM2,ISYM,WORK(KTMAT),
     &                   CONMAT,WORK(KTHMAT),IPRINT)

             IF (IPRINT.GE.5) THEN
              WRITE(LUPRI,'(4A/)')
     &         'GETCON: Transformed term {}', PRPNAM(INDXPR),';',
     &          PRPNAM(IPCON(ICOMP))
              CALL PRQMAT(WORK(KTHMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
             ENDIF

C ... calculate the expectation value of the {T_ic,h(mK)} matrix
              EXPVAL = D0
              DO IFS = 1, NFSYM
              IF (NOCC(IFS).GT.0) THEN
C ... run only over occupied electronic states
              DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS)
               IADR1  = IORB(IFS)+I
               IADR2  = (IADR1-1)*NORBT+IADR1 ! diagonal element
               EXPVAL = EXPVAL+WORK(KTHMAT+IADR2-1)
              ENDDO
              ENDIF
              ENDDO

C ... store the reorthonormalization contribution into ORTO
              ITMP = IPTAX(ICOMP,2)
              ORTO(ITMP,NOP) = DM2*EXPVAL

              IF (IPRINT.GE.4) THEN
                 WRITE(LUPRI,'(4A,1X,E15.8,/)')
     &           'GETCON: Expect. value {Tc,h(m_K)} for ',
     &            PRPNAM(IPCON(ICOMP)),' and ',
     &            PRPNAM(INDXPR),
     &            ORTO(ITMP,NOP)
              ENDIF
           ENDDO ! Of x,y,z connection matrixes...
           DOPCT=DOPCTRA_SAVE ! restore the pict.change transf.
         ENDIF
      ENDDO
C
C close property file
C
      CLOSE(LU1INT,STATUS = 'KEEP')      

C ... release the memmory...
      CALL MEMREL('GETCON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETCON')
      RETURN
      END
#endif /* MOD_LONDON  */

#ifdef MOD_LONDON
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TWOLON */
      SUBROUTINE TWOLON(CMO,IBEIG,WORK,LWORK)
C***********************************************************************
C
C   PURPOSE: 
C     Construct the two electron contribution
C     to the property gradient dH~(B)/dB|B=0
C     of the differentiated Dirac Hamiltonian,dg(B)/dB|B=0,
C     when using London atomic orbitals.
C
C   On input: CMO - MO coefficients
C             IBEIG - boson info
C
C
C     Written by T. Enevoldsen Feb 4 1998
C     Added by MI & HJAaJ, 2002/2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "shells.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "iratdef.h"
#include "ccom.h"
#include "symmet.h"
C
      DIMENSION CMO(*),IBEIG(*),WORK(LWORK)
      CALL QENTER('TWOLON')
#include "memint.h"
C
C     Memory allocation
C
      NCOMP = 3 ! x.y,z components 
      NDMAT = 4 ! 4 quatern. components
C ... Fock matrices
      CALL MEMGET('REAL',KFMAT,N2BBASX*NDMAT*NCOMP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIREP,NDMAT*NCOMP,WORK,KFREE,LFREE)
C ... quaternionic densities
      CALL MEMGET('REAL',KDMAT,N2BBASX*NDMAT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFCT,NDMAT*NCOMP,WORK,KFREE,LFREE)
C.....Get skeleton AO-Fock matrix with two-electron London integrals
      CALL TWOLO1A(WORK(KFMAT),WORK(KDMAT),CMO,
     &             IBEIG,NCOMP,WORK(KIREP),WORK(KIFCT),
     &             WORK(KFREE),LFREE,IPRPRP)
      CALL MEMREL('TWOLON.TWOLON1A',WORK,KWORK,KDMAT,KFREE,LFREE)
C
      CALL MEMGET('REAL',KSKLTN,NDMAT*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC, NTBAS(0)*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KINDEX,(KMAX*(MAXREP+1)*KHK(NHTYP)+1)/IRAT,
     &             WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPOINT, (NTBAS(0)*(MAXREP + 1) + 1)/IRAT,
     &             WORK,KFREE,LFREE)
C.....Get full AO-Fock matrix with two-electron London integrals
      CALL RSKLOND(WORK(KFMAT),WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &             WORK(KPOINT),NDMAT,WORK(KIREP),IPRPRP)
      CALL MEMREL('TWOLON.RSKLOND',WORK,KWORK,KSKLTN,KFREE,LFREE)

      CALL MEMGET('REAL',KFOCK,N2BBASXQ,WORK,KFREE,LFREE)
      CALL TWOLO1B(WORK(KFMAT),WORK(KFOCK),NCOMP,
     &             WORK(KFREE),LFREE,IPRPRP)
      CALL MEMREL('TWOLON.TWOLO1B',WORK,KWORK,KFOCK,KFREE,LFREE)

!     it has to happen after TWOLO1B
!     because XC is done in SA-AO basis
      call memget('REAL',kdmat0,n2bbasxq,work,kfree,lfree)
      call add_dft_to_twolon(work(kfmat),work(kdmat0),cmo,
     &                       work(kfree),lfree,iprprp)
      call memrel('twolon.dft',work,kwork,kdmat0,kfree,lfree)
!fixme: radovan:
!     XC contribution with nz < 4 does not work
!     i think this is because twolon works with unpacked AO matrices
!     and XC with packed
!     i don't know why we don't work with packed AO matrices here
!     and at the moment don't feel like unpacking the dft contribution
!     because i would have to allocate extra space
!     (dft adds on top), sort to hermit sorting, unpack, etc.

      CALL MEMGET('REAL',KGBMO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL TWOLO1C(WORK(KGBMO),WORK(KFMAT),CMO,NCOMP,
     &             WORK(KFREE),LFREE,IPRPRP)
      CALL MEMREL('TWOLON.TWOLO1C',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('TWOLON')

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine add_dft_to_twolon(fmat,dmat0,cmo,work,lwork,iprint)
!
      use dirac_cfg
      use xc_driver
      use dft_cfg

      implicit none
!
#include "priunit.h"  
#include "dummy.h"   
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "shells.h"
#include "ccom.h"
#include "cbihr2.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "iratdef.h"
C
#include "aovec.h"
#include "blocks.h"
#include "dcbgen.h"
#include "dcbnmr.h"
#include "dcbfir.h"
C
      real(8) :: fmat(*), dmat0(*), cmo(*), work(*)
      integer :: isymop(3), ihrmop(3), lwork, iprint, nr_fmat, i

      CALL QENTER('add_dft_to_twolon')

      if (dirac_cfg_dft_calculation) then
      nr_fmat = 3

      call genden(dmat0, cmo, 1, iprint)

!     CALL DZERO(FMAT,NTBAS(0)*NTBAS(0)*NZ*3)
      DO I=1,3
!        isymop(i) = 1
         isymop(i) = isymax(i, 2) + 1
!        anti hermitian
         IHRMOP(I) = -1
!     write(*, *) 'debug before direct', i
!     CALL PRQMAT(fmat,NTBAS(0),NTBAS(0),NTBAS(0),
!    &             NTBAS(0),NZ,IPQTOQ(1,isymax(i, 2)),LUPRI)
      END DO

      call xc_mpi_wake_up_coworkers()
      if (dft_cfg_no_sdft) then
         call xc_integrator(xc_mat_dim           = ntbas(0),
     &                      xc_nz                = nz,
     &                      xc_dmat_0            = dmat0,
     &                      xc_nr_dmat           = 0,
     &                      xc_nr_fmat           = nr_fmat,
     &                      xc_fmat              = fmat,
     &                      xc_fmat_pg_sym       = isymop,
     &            xc_do_london_rhs_direct = .true.)
      else
         call xc_integrator(xc_mat_dim           = ntbas(0),
     &                      xc_nz                = nz,
     &                      xc_dmat_0            = dmat0,
     &                      xc_nr_dmat           = 0,
     &                      xc_nr_fmat           = nr_fmat,
     &                      xc_fmat              = fmat,
     &                      xc_fmat_pg_sym       = isymop,
     &            xc_response_order_mo = 1,
     &            xc_do_london_rhs_direct = .true.)
      end if
      call xc_mpi_release_coworkers()

!     DO I=1,3
!     write(*, *) 'debug after direct', i
!     CALL PRQMAT(fmat,NTBAS(0),NTBAS(0),NTBAS(0),
!    &             NTBAS(0),NZ,IPQTOQ(1,isymax(i, 2)),LUPRI)
!     END DO

      endif
      call qexit('add_dft_to_twolon')

      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TWOLO1A */
      SUBROUTINE TWOLO1A(FMAT,DMAT,CMO,IBEIG,NCOMP,IREPDM,IFCTYP,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Construct two electron part of differentiated Hamiltonian
C     when using London Orbitals, dg(B)/dB(p,q,r,s)_X/Y/Z|B=0
C
C     Written by T. Enevoldsen Feb 4 1998
C
C     Added by MI&HJAaJ, 2002
C
C     Last modif: Dec 2005, MI - extended for .SPINFREE case
C                 March 2006,MI
C
C***********************************************************************
      use dirac_cfg
      use xc_driver
      use dft_cfg

#include "implicit.h"
#include "priunit.h"  
#include "dummy.h"   
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "shells.h"
#include "ccom.h"
#include "cbihr2.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "iratdef.h"
C
#include "aovec.h"
#include "blocks.h"
#include "dcbgen.h"
#include "dcbnmr.h"
#include "dcbfir.h"
C
      DIMENSION FMAT(N2BBASX,4,NCOMP),
     &          IBEIG(*),
     &          DMAT(N2BBASXQ),CMO(*),WORK(*)
      DIMENSION IFCTYP(4,NCOMP),IREPDM(4,NCOMP)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3),MTOTTK(3)
      LOGICAL LBIT,NODV,NOPV,NOCONT,RETUR,
     &         FCKDDR,DDFOCK,LONDON,FIRST,PARCAL_SAVE

      CALL QENTER('TWOLO1A')
#include "memint.h"
C  ... Initialize

      NDMAT = 4 ! 4 matrixes for each quater.unit. - closed-shell

      ITYPE = -5  ! Integral type - London orbitals !
      MAXDIF = 1  ! First derivative...
      NODV = .TRUE.
      NOPV = .TRUE.
      IATOM = 0
      NOCONT = .FALSE.
      RETUR  = .FALSE.
      RELCAL = .TRUE.
      DDFOCK = .TRUE.
      FCKDDR = .TRUE.
      LONDON = .TRUE.  !
      FIRST = .TRUE.
C
      DO I=1,3
         ISYMOP(I) = 1
         !ISYMOP(I) = isymax(i,2)+1
C anti hermitisk
!radovan: ??? but it's hermitian
         IHRMOP(I) = 1
         !IHRMOP(I) = - 1
         IFCKOP(I) = 1
      END DO

C ...  routine in 'dirac/dirtwo.F'
C      Set up IFCTYP and IREPDM
      CALL SETFCK(IFCTYP,IREPDM,NCOMP,4,ISYMOP,IHRMOP,IFCKOP,IPRINT)
!gosia      write(*, *) 'ifctyp = ', ifctyp
!gosia      write(*, *) 'irepdm = ', irepdm
C
C     Get contravariant density matrix, no folding
C  --------------------------------------------------
       CALL GTCNVD(DMAT,CMO,IBEIG,IPRINT,WORK(KFREE),LWORK)

C
C     Prepare for two-electron contribution
C
      LWRK = LWORK - KFREE + 1
      CALL DZERO(FMAT,NCOMP*NDMAT*N2BBASX)
C
C     Check what type of two-electron integrals to include
C
#if defined (VAR_MPI) || defined (VAR_PVM)
      MTOTTK(1) = NLRGBL*(NLRGBL+1)/2 ! 
      MTOTTK(2) = NLRGBL*(NLRGBL+1)/2 ! MI: lenght OK ...
      MTOTTK(3) = NSMLBL*(NSMLBL+1)/2
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
#endif

! here if dft, then scale HF exchange
!      do i = 1, 3
!      do j = 1, 4
!        write(*, *) 'ifctyp(j, i)=', j, i, ifctyp(j, i)
!      enddo
!      enddo

CMI   ... fix the parallel run...
C     PARCAL_SAVE = PARCAL
C     PARCAL = .FALSE.

      CALL GETTIM(CPU1,WALL1)

! gosia:
!        open(88, file='twolo1a', status='unknown', form='formatted',
!     &       access = 'append')

      DO I2TYP = 1,3  ! LL,LS,SS
        IF (LBIT(INTNMR,I2TYP)) THEN
#if defined (VAR_MPI) || defined (VAR_PVM)
         IF (PARCAL) THEN
            CALL MEMGET('INTE',KNPOS,MTOTTK(I2TYP),WORK,KFREE,LFREE)
CMI      .... replace PARDRV by HER_PARDRV and IFTHRS by SCRFCK
            CALL DIRAC_PARCTL( 1 )
CMI ... 11.8.'03 - replaced SCRFCK by DUMMY 
CMI   - Cauchy-Schwarz inequality does not hold for dg(1,2)/dB integrals !
            CALL HER_PARDRV(WORK(KFREE),LWRK,FMAT,DMAT,NDMAT,IREPDM,
     &            IFCTYP,ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,
     &            TKTIME,.FALSE.,FIRST,WORK(KNPOS),MTOTTK(I2TYP),
     &            I2TYP,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &            DUMMY)
         ELSE
#endif
             CALL TWOINT(WORK(KFREE),LWRK,FMAT,DMAT,NDMAT,IREPDM,IFCTYP
     &        ,DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT
     &        ,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY
     &        ,I2TYP,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,RELCAL,.false.
     &        ,IDUMMY,DUMMY)


#if defined (VAR_MPI) || defined (VAR_PVM)
         ENDIF
#endif
         ENDIF
! gosia:
!        write(88, *) 'i2typ = ', i2typ
!        call prqmat(fmat(1, 1, i2typ), 
!     &              ntbas(0), ntbas(0), ntbas(0), ntbas(0), 
!     &              nz, ipqtoq(1, 0), 88)
      ENDDO
!gosia:
!        close(88)

C     PARCAL = PARCAL_SAVE

      CALL GETTIM(CPU2,WALL2)
C
      CPUTOT = CPU2 - CPU1
      WLLTOT = WALL2 - WALL1
      WRITE(LUPRI,'(/A,F13.4,A,F13.4,A/)') 
     &'Cpu time:',CPUTOT,' Wall time:',WLLTOT,
     &' used in call to TWOINT(in TWOLO1).'
      CALL FLSHFO(LUPRI)
      CALL QEXIT('TWOLO1A')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TWOLO1B */
      SUBROUTINE TWOLO1B(FMAT,FOCK,NCOMP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Construct two electron part of differentiated Hamiltonian
C     when using London Orbitals, dg(B)/dB(p,q,r,s)_X/Y/Z|B=0
C
C     Written by T. Enevoldsen Feb 4 1998
C
C     Added by MI&HJAaJ, 2002
C
C     Last modif: Dec 2005, MI - extended for .SPINFREE case
C                 March 2006,MI
C
C***********************************************************************
      use dirac_cfg
      use xc_driver
      use dft_cfg

#include "implicit.h"
#include "priunit.h"  
#include "dummy.h"   
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "shells.h"
#include "ccom.h"
#include "cbihr2.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "iratdef.h"
C
#include "aovec.h"
#include "blocks.h"
#include "dcbgen.h"
#include "dcbnmr.h"
#include "dcbfir.h"
C
      DIMENSION FMAT(N2BBASX,4,NCOMP),
     &          FOCK(N2BBASX,NZ),WORK(*)
      CALL QENTER('TWOLO1B')
#include "memint.h"

! gosia: check enetering fock matrix after rsklond:
!        open(77, file='after_rsklond', status='unknown', 
!     &       form='formatted', access = 'append')
!        do i2typ = 1, 3
!        write(77, *) 'gosia: i2typ, isymax(i2typ,2) = ', 
!     &                       i2typ, isymax(i2typ,2)
!        IREP = ISYMAX(i2typ,2)       ! irep of x,y and z
!        call prqmat(fmat(1, 1, i2typ), 
!     &              ntbas(0), ntbas(0), ntbas(0), ntbas(0), 
!     &              nz, ipqtoq(1, irep), 77)
!        end do
!        close(77)
!end
!        open(66, file='twolo1b', status='unknown', form='formatted',
!     &       access = 'append')
C----------------------------------------------------------------
C       Transform Fock matrix from AO to QO
C       ===================================
C----------------------------------------------------------------
      DO I = 1, NCOMP  ! x,y,z components
         CALL DZERO(FOCK,N2BBASXQ)
         IREP = ISYMAX(I,2)       ! irep of x,y and z
         IPARI = JBTOF(IREP,1)    !  ??parity ??

!gosia         write(lupri, *) 'twolo1b, i     = ', i    
!gosia         write(lupri, *) 'twolo1b, irep  = ', irep
!gosia         write(lupri, *) 'twolo1b, ipari = ', ipari

         DO IZ = 1,4
            IREPD = IRQMAT(IZ,IREP)
            IQ    = IQMULT(1,JQBAS(IREPD,IPARI),IZ)            
            IPQ   = IQTOPQ(IQ,IREP)

!gosia            write(lupri, *) 'twolo1b, iz    = ', iz
!gosia            write(lupri, *) 'twolo1b, irepd = ', irepd
!gosia            write(lupri, *) 'twolo1b, iq    = ', iq
!gosia            write(lupri, *) 'twolo1b, ipq   = ', ipq

C ...    transform AO - SO
            CALL MTAOSO(FMAT(1,IZ,I),FOCK(1,IPQ),
     &           NTBAS(0),IREPD,IPRINT)
         ENDDO

         IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &'TWOLO1: Transformed g(B) matrix in QO basis',-1)
          WRITE(LUPRI,'(2X,A,I2)') '** component: ',I
          CALL PRQMAT(FOCK(1,1),NTBAS(0),NTBAS(0),NTBAS(0),
     &                 NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         ENDIF  
C
C     Insert half-phases and symmetrize
C     =================================
C
         IF(NZ.LT.4) THEN
           DO IZ = 1,NZ
             IQ = IPQTOQ(IZ,IREP)
             CALL Q2PHASE('D',IQ,1,FOCK(1,IZ))
           ENDDO
           IF(IPRINT.GE.5) THEN
             CALL HEADER(
     & 'TWOLO1: Unsorted g(B) matrix with half-phases inserted',-1)
             WRITE(LUPRI,'(2X,A,I2)') '** component: ',I
             CALL PRQMAT(FOCK(1,1),NTBAS(0),NTBAS(0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
           ENDIF
         ENDIF

C ... finally get g(B) into SA-AO basis - copy to lower triangle of symm.matrix
         DO IZ = 1, NZ
           CALL DSLTSP(NTBAS(0),FOCK(1,IZ),FMAT(1,IZ,I))
         END DO

C ... unpack to the full AO - form
         CALL DZERO(FOCK,N2BBASXQ)
         ITIM = -1 ! g(B)_comp derivatives are TR antisymmetric 
         DO IZ = 1, NZ
          IQ = IPQTOQ(IZ,IREP)
          IH = IHQMAT(IQ,ITIM)
          IF    (IH.EQ.1) THEN
           CALL DSPTSI(NTBAS(0),FMAT(1,IZ,I),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
           CALL DAPTGE(NTBAS(0),FMAT(1,IZ,I),FOCK(1,IZ))
          ELSE
           CALL QUIT('TWOLO1: Not defined IH!')
          ENDIF
         ENDDO

C     ...sort AO to SO-AO, so somehow ...
         CALL BUTOBS(FOCK,NZ,WORK(KFREE),LFREE)
C ...    
         CALL DCOPY(N2BBASXQ,FOCK,1,FMAT(1,1,I),1)
         IF(IPRINT.GE.5) THEN
           CALL HEADER(
     & 'TWOLO1:Total SA-AO two-electron London contribution',-1)
           WRITE(LUPRI,'(2X,A,I2)') '** component: ',I
           CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),
     &            NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         ENDIF
!
!        write(66, *) 'gosia: i, irep = ', i, irep
!        call prqmat(fmat(1, 1, i), 
!     &              ntbas(0), ntbas(0), ntbas(0), ntbas(0), 
!     &              nz, ipqtoq(1, irep), 66)
!
      ENDDO ! Of I = x,y,z
!      close(66)
      CALL QEXIT('TWOLO1B')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TWOLO1C */
      SUBROUTINE TWOLO1C(GBMO,FMAT,CMO,NCOMP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Construct two electron part of differentiated Hamiltonian
C     when using London Orbitals, dg(B)/dB(p,q,r,s)_X/Y/Z|B=0
C
C     Written by T. Enevoldsen Feb 4 1998
C
C     Added by MI&HJAaJ, 2002
C
C     Last modif: Dec 2005, MI - extended for .SPINFREE case
C                 March 2006,MI
C
C***********************************************************************
      use dirac_cfg
      use xc_driver
      use dft_cfg

#include "implicit.h"
#include "priunit.h"  
#include "dummy.h"   
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "shells.h"
#include "ccom.h"
#include "cbihr2.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "iratdef.h"
C
#include "aovec.h"
#include "blocks.h"
#include "dcbgen.h"
#include "dcbnmr.h"
#include "dcbfir.h"
C
      DIMENSION GBMO(N2ORBXQ),FMAT(N2BBASX,4,NCOMP),
     &          CMO(*),WORK(*)
      CALL QENTER('TWOLO1C')
#include "memint.h"
!gosia:
!        open(55, file='twolo1c', status='unknown', form='formatted',
!     &       access = 'append')

C----------------------------------------------------------------
C   Open new file for MO-transformed g(B) x,y,z matrixes
C----------------------------------------------------------------
      OPEN(LUTWOL,FILE='TWOLON',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='UNKNOWN')
      DO I = 1, NCOMP  ! x,y,z components
         IREP = ISYMAX(I,2)       ! irep of x,y and z
C
C     Transform to MO-basis only the blocks of correct fermion symmetry
C     =================================================================
C
C   ....  it would be good to add also QBTRANS for spin-free case
     
      CALL DZERO(GBMO,N2ORBXQ)
      IOPSY = JBTOF(IREP,1)
      DO I1 = 1, NFSYM
         I2 = MOD(I1+IOPSY,2) + 1
         CALL QTRANS('AOMO','S',D0,
     &      NFBAS(I1,0),NFBAS(I2,0),
     &      NORB(I1),NORB(I2),
     &     FMAT(I2BASX(I1,I2)+1,1,I),
     &      NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
     &     GBMO(1+I2ORBX(I1,I2)),
     &      NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &     CMO(ICMOQ(I1)+1),
     &      NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &     CMO(ICMOQ(I2)+1),
     &      NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)
      ENDDO

      IF(IPRINT.GE.5) THEN
        CALL HEADER(
     & 'TWOLO1: Two-electron London contribution in MO',-1)
        WRITE(LUPRI,'(2X,A,I2)') '** component: ',I
        CALL PRQMAT(GBMO,NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP),LUPRI)
      ENDIF

C .... write the two-electron London contribution to the file...
!gosia       write(*, *) 'twolo1c, i=', i
       CALL WRTDAC(LUTWOL,N2ORBXQ,GBMO,I)

!gosia
!        write(55, *) 'gosia: i = ', i
!        call prqmat(gbmo, 
!     &              ntbas(0), ntbas(0), ntbas(0), ntbas(0), 
!     &              nz, ipqtoq(1, irep), 55)

      ENDDO ! Of I = x,y,z
!gosia:
!      close(55)

      CLOSE(LUTWOL,STATUS='KEEP')
      CALL QEXIT('TWOLO1C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck TTRA */
      SUBROUTINE TTRA(EPMIX,AH,ISYMT,ISYMH,TMAT,HMAT,THMAT,IPRINT)
C*****************************************************************************
C
C   PURPOSE:
C ===========
C    Do one-index transformation of entering MO-matrixes, {T,h}.
C    This is expression is used for reorthonormalization terms when
C    employing London atomic orbitals.
C
C Definition (Ruud et al., Chem.Phys. 195(1995) 157):
C -----------------------------------------------------
C  {T,h}= THMAT_{mn} = Sum_{o}(TMAT_{om}^*.HMAT_{on} + HMAT_{mo}.TMAT_{no}) 
C
C           in matrix notation:    T^{+}.H + H.T
C
C COMMENT: Helgaker,Jorgensen, JCP 95(1991) 2595 prefer (1 1^*) definition,
C          however, definition (1^* 1) is correct here because 
C          we have the same type
C          of one- and two-electron (1^* 1|2^* 2) matrix elements.
C          We have to obey the convention carefully since we are
C          working with general (quaternion,complex) quantities.
C
C
C On input:
C-----------
C  EPMIX = .true - run o-summation over both positronic and electronic shells,
C   else (EPMIX=.false.) run only over electronic shells.
C
C  AH = .true. - take -1 for Hermitian conjugate, otherwise take +1
C
C  ISYMT,ISYMH ... boson symmetries on entering operators T,h
C
C  Note: The TMAT might also be a general matrix,
C        the entering h is usually anti-/hermitian
C
C  Written by Miro Ilias (based on TEC's routine TECTRA), august 2003, Odense.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, DM1 = -1.0D0)
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      LOGICAL EPMIX, AH
      DIMENSION  TMAT(NORBT,NORBT,NZ), 
     &           HMAT(NORBT,NORBT,NZ),
     &          THMAT(NORBT,NORBT,NZ)
C
#include "ibtfun.h"
      CALL QENTER('TTRA')

      IREPT = ISYMT-1
      IOPSYT = JBTOF(IREPT,1)

      IREPH = ISYMH-1
      IOPSYH = JBTOF(IREPH,1)

      IREPTH = IBTXOR(IREPT,IREPH) ! Symm.product of G(T)xG(H)
      IOPSYTH = JBTOF(IREPTH,1)

C ... specify sign for herm.conjugate of TMAT 
      IF (AH) THEN ! anti-hermitian
        DSGN = DM1
      ELSE
        DSGN = D1
      ENDIF

      IF(IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: Entering total H MO matrix of boson symmetry ',
     & ISYMH,'/',NBSYM
       CALL PRQMAT(HMAT,NORBT,NORBT,NORBT,
     &             NORBT,NZ,IPQTOQ(1,IREPH),LUPRI)
       WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: Entering total T MO matrix of boson symmetry ',
     & ISYMT,'/',NBSYM
       CALL PRQMAT(TMAT,NORBT,NORBT,NORBT,
     &             NORBT,NZ,IPQTOQ(1,IREPT),LUPRI)
      ENDIF

C ... do matrix multiplications - run over fermion irreps
C ... (I1,I2) - f.irreps of T, (I2,I3) - f.irreps of H

      DO I1 = 1, NFSYM

         I2 = MOD(I1+IOPSYT,2) + 1  ! I1,I2 ... T
         I3 = MOD(I1+IOPSYH,2) + 1  ! I1,I3 ... H

C        Write(lupri,'(2x,a,4i3)') 'I1,I2,I3,I4:',I1,I2,I3,I4
C        Write(lupri,'(2x,a,3i3)')
C    & 'IREPT,IREPH,IREPTH:',IREPT,IREPH,IREPTH

         IF (EPMIX) THEN ! Matrix multiplication in full e+p MO space

C        ... do TH(I2,I3) := dsgn*T^{H}(I2,I1).H(I1,I3)
          CALL QGEMM(NORB(I2),NORB(I3),NORB(I1),DSGN,
     &    'H','N',IPQTOQ(1,IREPT),
     &    TMAT(IORB(I2)+1,IORB(I1)+1,1),NORBT,NORBT,NZ,
     &    'N','N',IPQTOQ(1,IREPH),
     &    HMAT(IORB(I1)+1,IORB(I3)+1,1),NORBT,NORBT,NZ, 
     &    D0,IPQTOQ(1,IREPTH),
     &    THMAT(IORB(I2)+1,IORB(I3)+1,1),NORBT,NORBT,NZ)

          IF(IPRINT.GE.5) THEN
            WRITE(LUPRI,'(/2X,A,I1,A,I1,F6.2)')
     &    'TTRA: 1.part transformed T^{+}.H over e-p :',
     &     IREPTH+1,'/',NBSYM,DSGN
            CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
          ENDIF

C  ...  do TH(I3,I2) := TH(I3,I2) + H(I3,I1).T(I1,I2)
          CALL QGEMM(NORB(I3),NORB(I2),NORB(I1),D1,
     &   'N','N',IPQTOQ(1,IREPH),
     &   HMAT(IORB(I3)+1,IORB(I1)+1,1),NORBT,NORBT,NZ,
     &   'N','N',IPQTOQ(1,IREPT),
     &   TMAT(IORB(I1)+1,IORB(I2)+1,1),NORBT,NORBT,NZ, 
     &   D1,IPQTOQ(1,IREPTH),
     &   THMAT(IORB(I3)+1,IORB(I2)+1,1),NORBT,NORBT,NZ)

         IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: TOTAL transformed {T,H} over e-p MO matrix of symm ',
     & IREPTH+1,'/',NBSYM
          CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &           NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
         ENDIF

         ELSE  ! Transformation in the e-e space only

C        ... do TH(I2,I1) = dsgn*T^{H}(I2,I1).H(I1,I3)
          CALL QGEMM(NESH(I2),NESH(I3),NESH(I1),DSGN,
     &    'H','N',IPQTOQ(1,IREPT),
     &    TMAT(IORB(I2)+NPSH(I2)+1,IORB(I1)+NPSH(I1)+1,1),
     &    NORBT,NORBT,NZ,
     &    'N','N',IPQTOQ(1,IREPH),
     &    HMAT(IORB(I1)+NPSH(I1)+1,IORB(I3)+NPSH(I3)+1,1),
     &    NORBT,NORBT,NZ,D0,IPQTOQ(1,IREPTH),
     &    THMAT(IORB(I2)+NPSH(I2)+1,IORB(I3)+NPSH(I3)+1,1),
     &    NORBT,NORBT,NZ)

          IF(IPRINT.GE.5) THEN
            WRITE(LUPRI,'(/2X,A,I1,A,I1,F6.2)')
     &    'TTRA: 1.part transformed T^{+}.H over e:',
     &     IREPTH+1,'/',NBSYM,DSGN
            CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
          ENDIF

C        ...  do TH(I3,I2) = TH(I3,I2) + H(I3,I1)*T(I1,I2)
          CALL QGEMM(NESH(I3),NESH(I2),NESH(I1),D1,
     &   'N','N',IPQTOQ(1,IREPH),
     &   HMAT(IORB(I3)+NPSH(I3)+1,IORB(I1)+NPSH(I1)+1,1),
     &   NORBT,NORBT,NZ,
     &   'N','N',IPQTOQ(1,IREPT),
     &   TMAT(IORB(I1)+NPSH(I1)+1,IORB(I2)+NPSH(I2)+1,1),
     &   NORBT,NORBT,NZ,D1,IPQTOQ(1,IREPTH),
     &   THMAT(IORB(I3)+NPSH(I3)+1,IORB(I2)+NPSH(I2)+1,1),
     &   NORBT,NORBT,NZ)

         IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: TOTAL transformed {T,H} MO matrix over e of symm ',
     & IREPTH+1,'/',NBSYM
          CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &           NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
         ENDIF

        ENDIF

      ENDDO

      CALL QEXIT('TTRA')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck RNTLON */
      SUBROUTINE RNTLON(CMO,IBEIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C   PURPOSE: Calculate the {T(1),H(0)} contribution
C            to the property gradient matrix dH~(B)/dB|B=0
C            for NMR shieldings and magnetizabilities
C            with London atomic orbitals.
C            Write resulting matrix (in the MO basis)
C            to the DA file 'FCKLON'. 
C
C   Definition (Helgaker et al: JCP  84(1986) 6266)
C -----------------------------------------------------
C     {T(1),H(0)} =  {T(1),F(0)} + W(0)
C
C   Called from:  pamxlr.F/PRPXL1
C
C On input: CMO - MO coefficientsa
C           IBEIG  - info on boson irreps (LevyLe or SpinFREE)
C
C   Written by M.Ilias & HJAaJ
C     (based on T. Enevoldsen's codes), august 2003, Odense 
C
C   Last revisons: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"  
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"

      PARAMETER (D1=1.0D0,D0=0.0D0)

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      DIMENSION ISYMOP(3),IHRMOP(3),IFCKOP(3)
      DIMENSION IREP(3),IOPSY(3),ISYM(3)
      DIMENSION KF1TR(3),KF2TR(3),KMDM(3)

      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE


      CALL QENTER('RNTLON')
      KFRSAV = KFREE
C Open file for the London reorthonormalization contribution, {T(1),H(0)} 
      OPEN(LUFCKL,FILE='FCKLON',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='UNKNOWN') ! Why 8* and not 3*?

C Get the total Fock matrix in MO basis
      CALL MEMGET('REAL',KFMAT,N2BBASXQ,WORK,KFREE,LFREE)
C    .... prp/pamset.F
      CALL GETFCK(WORK(KFMAT),IPRPRP,WORK,KFREE,LFREE) 
      IPRINT = IPRPRP

      IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &'RNTLON: Total Fock matrix in MO basis:',-1)
        CALL PRQMAT(WORK(KFMAT),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

C ... get space for one-index transformed Fock matrixes,
C     {T(1),F(0)}_x/y/z, in MO basis - 1. part of {T(1),H(0)}
      CALL MEMGET('REAL',KF1TR(1),3*N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KF1TRC,N2ORBXQ,WORK,KFREE,LFREE)
      KF1TR(2) = KF1TR(1) + N2ORBXQ
      KF1TR(3) = KF1TR(2) + N2ORBXQ


C ... get space for (1) connection matrix, T(1)_x/y/z, in MO basis
      CALL MEMGET('REAL',KTMATC,N2ORBXQ,WORK,KFREE,LFREE)

C  ... get the space for all three symmetrized 
C      modified density matrixes, MDM_x/y/z, in AO basis
      CALL MEMGET('REAL',KMDM(1),3*N2BBASXQ,WORK,KFREE,LFREE)
      KMDM(2) = KMDM(1) + N2BBASXQ
      KMDM(3) = KMDM(2) + N2BBASXQ

      CALL MEMGET('REAL',KMDMC,N2BBASXQ,WORK,KFREE,LFREE)

C ... define e-p/e shells for reorthonormalization
      EPMIX = DOEPRN
C Open property file (needed connection matrixes, Fock MO matrix)
      CALL OPNFIL(LU1INT,'AOPROPER','OLD','RNTLON')
C
      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.
C ... run over components x/y/z
      DO ICOMP = 1, 3
C ... get connection matrixes (may be SQUARE for natural connection!),
C     T(1)_ICOMP, in MO basis, into WORK(KTMATC)
       INDXPR          =  IPCON(ICOMP)
       ITIM            =  IPRPTIM(INDXPR)
       ISYM(ICOMP)     =  IPRPSYM(INDXPR)
       IOPSY(ICOMP)    =  JBTOF(ISYM(ICOMP)-1,1)
       IREP(ICOMP)     =  ISYM(ICOMP) - 1 
       CALL DZERO(WORK(KTMATC),N2ORBXQ)
C ...      routine generalized for square matrixes too... 
       CALL PRPMAT(INDXPR,IOPSY,WORK(KTMATC),CMO,
     &           IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)
       IF(IPRINT.GE.5) THEN
        CALL HEADER(
     &  'RNTLON: T(1)_ic connection matrix in MO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(KTMATC),NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
       ENDIF
C
C ...  get 1.transformed term, {T(1),F(0)}_ICOMP, in MO basis
C      into WORK(KF1TRC)
C
       ITYP = 1 ! e-e only !
       HERM  = .FALSE.  ! antisymmetric
       CALL DZERO(WORK(KF1TRC),N2ORBXQ)

       ISYMH = 1 ! Do {T,F}
       AH = .TRUE. ! 
       CALL TTRA(EPMIX,AH,ISYM(ICOMP),ISYMH,WORK(KTMATC),
     &           WORK(KFMAT),WORK(KF1TRC),IPRPRP)

C  ... and store them into WORK(KF1TR(1-3))
       CALL DCOPY(N2ORBXQ,WORK(KF1TRC),1,WORK(KF1TR(ICOMP)),1)

       IF(IPRINT.GE.5) THEN
        CALL HEADER(
     &  'RNTLON: Transformed {T(1),F(0)}_icomp in MO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(KF1TRC),NORBT,NORBT,NORBT,NORBT,
     &               NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
       ENDIF
  
C ...  get symmetrized modified density matrix in AO basis
C ...  into WORK(KMDMC)
       ITYP = 1 
       IH = -1 ! Anti-hermitian
       CALL DZERO(WORK(KMDMC),N2BBASXQ)

       CALL GETMDM(EPMIX,ISYM(ICOMP),
     &            WORK(KMDMC),CMO,IBEIG,WORK(KTMATC),
     &            WORK,KFREE,LFREE,IPRPRP)

C ... and store it into WORK(KMDM(ICOMP))
       CALL DCOPY(N2BBASXQ,WORK(KMDMC),1,WORK(KMDM(ICOMP)),1)

       IF(IPRINT.GE.5) THEN
         CALL HEADER(
     &'RNTLON: Total SYMMETRIZED modified density matrix',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
         CALL PRQMAT(WORK(KMDM(ICOMP)),NTBAS(0),NTBAS(0),
     &     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
         CALL FLSHFO(LUPRI)
       ENDIF

      ENDDO ! x/y/z

C===================================================================
C Needed: space for {T(1),g(0)_add} and mod.symm.dens.matrixes
C===================================================================
C  Get space for two-electron contributions, {T(1),g(0)_add}  in AO basis
C===================================================================
      CALL MEMGET('REAL',KF2TR(1),3*N2BBASXQ,WORK,KFREE,LFREE)
      KF2TR(2) = KF2TR(1) + N2BBASXQ 
      KF2TR(3) = KF2TR(2) + N2BBASXQ 

      ISYMOP(1) = ISYM(1)
      ISYMOP(2) = ISYM(2)
      ISYMOP(3) = ISYM(3)
      IHRMOP(1) = -1
      IHRMOP(2) = -1
      IHRMOP(3) = -1
C     IFCKOP(1) = 1
      IFCKOP(1) = JBTOF(ISYM(1)-1,1)
C     IFCKOP(2) = 1
      IFCKOP(2) = JBTOF(ISYM(2)-1,1)
C     IFCKOP(3) = 1
      IFCKOP(3) = JBTOF(ISYM(3)-1,1)

      NFMTX = 3 ! Number of Fock matrixes - 3(x,y,z)

      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /)) 
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif

C ... get all 3 terms, {T(1),g(0)-add}_x/y/z,
C ...  into WORK(KF2TR(1-3)) AT ONCE !
      CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,
     &      WORK(KF2TR(1)),WORK(KMDM(1)),
     &      NFMTX,WORK(KPOS),INTNMR,
     &      IPRINT,WORK(KFREE),LFREE)



!!!!!do transformation AO- MO !!!!!!!!!!!!!11
C  ... run over components
      DO ICOMP = 1, 3

C  copy WORK(KF2TR(ICOMP)) -> WORK(KMDMC)
        CALL DCOPY(N2BBASXQ, WORK(KF2TR(ICOMP)),
     &                    1,WORK(KMDMC),1)

C  Control output
       IF(IPRINT.GE.5) THEN
        CALL HEADER(
     &  'RNTLON: {T(1),g(0)_add}_ic in AO basis ',-1)
        WRITE(LUPRI,'(2X,A,I1)') '** component:',ICOMP
        CALL PRQMAT(WORK(KMDMC),NTBAS(0),NTBAS(0),
C       CALL PRQMAT(WORK(KF2TR(ICOMP)),NTBAS(0),NTBAS(0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP(ICOMP)),LUPRI)
       ENDIF

C  Make transformation of {T(1),g(0)_add}_icomp from AO to the MO basis.
C  Transform to MO-basis only the blocks of correct fermion symmetry
C  WORK(KF2TR(ICOMP)) ->  WORK(KTMATC)
C  =================================================================
C  ... it would be good to extend it using QBTRANS....
       CALL DZERO(WORK(KTMATC),N2ORBXQ)

       DO I1 = 1, NFSYM
          I2 = MOD(I1+IOPSY(ICOMP),2) + 1
          CALL QTRANS('AOMO','S',D0,
     &     NFBAS(I1,0),NFBAS(I2,0),
     &     NORB(I1),NORB(I2),
     &     WORK(KMDMC+I2BASX(I1,I2)),   !result in AO basis
     &     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP(ICOMP)),
     &     WORK(KTMATC+I2ORBX(I1,I2)),  !space for new result in MO basis
     &     NORBT,NORBT,NZ,IPQTOQ(1,IREP(ICOMP)),
     &     CMO(1+ICMOQ(I1)),
     &     NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &     CMO(1+ICMOQ(I2)),
     &     NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)
       ENDDO

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &'RNTLON: Total transf. term, {T(1),g(0)_add}_ic, in MO:',-1)
         WRITE(LUPRI,'(2X,A,I2)') '*** component :',ICOMP
         CALL PRQMAT(WORK(KTMATC),NORBT,NORBT,
     &               NORBT,NORBT,NZ,
     &            IPQTOQ(1,IREP(ICOMP)),LUPRI)
       ENDIF

C Sum {T(1),F(0)}_c + {T(1),g(0)_add}_c into WORK(KF1TR)
C         ****  F1TR = F1TR + F2TR  ****
C  MI(27.11'03) - D2h symmetry does not work correctly for F2TR evaluation!
C      CALL DAXPY(N2ORBXQ,D0,WORK(KTMATC),1,WORK(KF1TR(ICOMP)),1) 
      CALL DAXPY(N2ORBXQ,D1,WORK(KTMATC),1,WORK(KF1TR(ICOMP)),1) 
      IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &'RNTLON: Total transformed term, {T(1),H(0)}_comp, in MO:',-1)
        WRITE(LUPRI,'(2X,A,I1)') '*** Component x-1/y-2/z-3:',ICOMP
        CALL PRQMAT(WORK(KF1TR(ICOMP)),
     &              NORBT,NORBT,NORBT,NORBT,NZ,
     &            IPQTOQ(1,IREP(ICOMP)),LUPRI)
       ENDIF

C  Write the total reorthonormalization term to the DA file
       CALL WRTDAC(LUFCKL,N2ORBXQ,WORK(KF1TR(ICOMP)),ICOMP)

      ENDDO    

C ... restore the picture change transformation keyword
      DOPCT=DOPCTRA_SAVE

C
C close property and FCKLON files
C
      CLOSE(LU1INT,STATUS = 'KEEP')
      CLOSE(LUFCKL,STATUS = 'KEEP') 

C     .... release the whole block used in this subroutine
      CALL MEMREL('RNTLON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('RNTLON')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onelon */
      SUBROUTINE ONELON(CMO,IBEIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C   PURPOSE:
C ==========
C  Get the total one-electron magnetic field property gradient dh~(B)/dB
C  in MO basis when using London orbitals and write it to the DA file ONELON.
C
C  Expression:  RM1H1 + RM1H2 + RM1H3(sq) + RM1RN(sq) (x,y,z)
C
C  On output: CMO exports MO coefficients
C             IBEIG and boson irreps info...
C
C  Written by Miro Ilias, august 2003, Odense  
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"  
      PARAMETER(D0 = 0.0D0, D1= 1.0D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dummy.h"

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      DATA LUONELO / 72 /

      CALL QENTER('ONELON')

      KFRSAV = KFREE

      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','ONELON')
C ... distinguish between spin-free (boson symm.) mode
C     and the full relativistic (fermion symm.) mode
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident. 
       CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IBEIG,
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IDUM,
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')

C Open file for the London one-electron contributions, dh(B)/dB_c
      OPEN(LUONELO,FILE='ONELON',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='UNKNOWN') ! Why 8* and not 3*?

      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')

C Allocate memory for storing components of dh(B)/dB_c

      CALL MEMGET('REAL',KDHB,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDHBIX,N2ORBXQ,WORK,KFREE,LFREE)

      DO ICOMP = 1, 3
       CALL DZERO(WORK(KDHB),N2ORBXQ)
       DO IX = 1, 4
         CALL DZERO(WORK(KDHBIX),N2ORBXQ)
         INDXPR=IPLONDON(IX,ICOMP)
         ITIM   = IPRPTIM(INDXPR)
         IREP   = IPRPSYM(INDXPR) - 1
         CALL PRPMAT(INDXPR,IOPSY,WORK(KDHBIX),CMO,
     &      IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)
         IF (IPRPRP.GE.5) THEN
           WRITE(LUPRI,'(/A,A16)') 
     &     'ONELON: *MO-matrix of property ',
     &      PRPNAM(INDXPR)
            CALL PRQMAT(WORK(KDHBIX),NORBT,NORBT,
     &             NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
         ENDIF

C   ...  accumulate all 1-el. prop. grad. contributions
         CALL DAXPY(N2ORBXQ,D1,WORK(KDHBIX),1,WORK(KDHB),1)
       ENDDO

C ... and write it the file:
       CALL WRTDAC(LUONELO,N2ORBXQ,WORK(KDHB),ICOMP)
       IF (IPRPRP.GE.5) THEN
         CALL HEADER(
     &'ONELON: The total LAO 1 el. mag.field gradient MO matrix:',-1)
         WRITE(LUPRI,'(2X,A,I2)') '** component:',ICOMP
         CALL PRQMAT(WORK(KDHB),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
       ENDIF
      ENDDO

      CLOSE(LUONELO,STATUS='KEEP')
      CLOSE(LU1INT,STATUS='KEEP')

C     .... release the whole block used in this subroutine
      CALL MEMREL('ONELON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('ONELON')
      RETURN
      END
#endif /*  MOD_LONDO  */

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

!gosia fixme
!       write(lupri, *) 'density matrix before twoint'
!       call prqmat(dmat(1, 1, 1),
!     &           ntbas(0), ntbas(0), ntbas(0), ntbas(0), 
!     &           nz, ipqtoq(1, 0), lupri)
!


C     Dispose used memory
      CALL MEMREL('GTCNVD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GTCNVD')
      RETURN
      END

#ifdef MOD_LONDON
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck EXPVSUSLO  */
      SUBROUTINE EXPVSUSLO(WORK,LWORK)
C*****************************************************************************
C
C  Routine that calculates expectation value terms for magnetic
C  susceptibilities - two-electron term and both higher and lower order
C  reorthonormalization terms.
C
C  Written by Miroslav Ilias, Prievidza, june 2003.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dummy.h"
C
      DIMENSION WORK(*)

      CALL QENTER('EXPVSUSLO')
#include "memint.h" 

C     WRITE(LUPRI,*) "IN DIASUSLO routine!"


C 1. Get tensor of two-electron expectation value terms 

C--- Try to get the exp.value over two-electron integrals !!!
C prepare all arrays and density matrixes !
C
      NDMAT = 4 ! Number of dens.matrixes (closed shell)
      NCOMP = 6 ! tensor components 

      CALL MEMGET('REAL',KFMAT,N2BBASX*NDMAT*NCOMP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,N2BBASX*NDMAT,WORK,KFREE,LFREE) ! Contravariant density matrixes
      CALL MEMGET('REAL',KFOCK,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE) ! MO coefficients

      CALL MEMGET('INTE',KIREP,NDMAT*NCOMP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFCT,NDMAT*NCOMP,WORK,KFREE,LFREE)

C ... boson irreps identification
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)

C get MO coefficients and boson irreps info
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','EXVSUSLO')
C ... distinguish between spin-free (boson symm.) mode
C     and the full relativistic (fermion symm.) mode
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident.
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,WORK(KBEIG),
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,IDUM,
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')

      CALL SUSC2EL(WORK(KFMAT), WORK(KFOCK), WORK(KBEIG),
     &     WORK(KCMO),WORK(KDMAT),WORK(KIREP),WORK(KIFCT),
     &     NCOMP,WORK(KFREE),LFREE)

C  2. Compute highest order reorthonormalization term {T(2),h(0)+g(0)}

      CALL HORDRNT(WORK(KCMO),WORK(KBEIG),WORK,KFREE,LFREE)

C  3. Compute lowest order reorthonormalization term  ...many more terms...

      CALL LORDRNT(WORK(KCMO),WORK(KBEIG),WORK,KFREE,LFREE)

      CALL MEMREL('EXPVSUSLO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('EXPVSUSLO')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SUSC2EL  */
      SUBROUTINE SUSC2EL(FMAT,FOCK,IBEIG,
     &             CMO,DMAT,IREPDM,IFCTYP,NCOMP,
     &                      WORK,LWORK)
C*****************************************************************************
C
C   Routine that calculates two-electron expectation value contribution
C to magnetic susceptibility when using London atomic orbitals.
C
C     Written by Miroslav Ilias, Prievidza, june 2003.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"  
#include "dummy.h"   
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "shells.h"
#include "ccom.h"
#include "cbihr2.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "iratdef.h"
C
#include "aovec.h"
#include "blocks.h"
#include "dcbgen.h"
#include "dcbnmr.h"
#include "dcbprp.h"
#include "dcbfir.h"
C
      DIMENSION FMAT(N2BBASX,4,NCOMP),FOCK(N2BBASX,NZ)
      DIMENSION DMAT(N2BBASX*4),CMO(*),IBEIG(*),WORK(*)
      DIMENSION IFCTYP(4,NCOMP),IREPDM(4,NCOMP)
      DIMENSION ISYMOP(6),IHRMOP(6),IFCKOP(6),MTOTTK(3)
      LOGICAL   LBIT,NODV,NOPV,NOCONT,RETUR,
     &          FCKDDR,DDFOCK,FIRST,SUSCEP
      LOGICAL PARCAL_SAVE

      CALL QENTER('SUSC2EL')
#include "memint.h"

      NDMAT = 4
      NCOMP = 6

C     get contavariant density matrixes

!gosia: gtcnvd won't always work: 
!here for dmat we allocate n2bbasx*4 memory, 
!while inside gtcnvd we need n2bbasx*4*nrfmat memory for this dmat
!(here it is ok, as we loop over 1 'nrfmat' only inside gtcnvd)

      CALL GTCNVD(DMAT,CMO,IBEIG,IPRPRP,WORK,LWORK)

C set up for second derivative g12 over B
      DO I=1,NCOMP
         ISYMOP(I) = 1
         IHRMOP(I) = 1
         IFCKOP(I) = 1
      END DO
C prepare IFCTYP and IREPDM
      CALL SETFCK(IFCTYP,IREPDM,NCOMP,4,ISYMOP,IHRMOP,IFCKOP,IPRPRP)

C call TWOINT and modify routines inside her2ave

C ...Integral type - London orbitals and without EXPECT in her2ave !!!
      ITYPE = -5  
CMI   EXPECT is for nucl.coordinate derivates !!!
      MAXDIF = 2  ! Second derivate...=> SUSCEPTIBILITIES !!
      NODV = .TRUE.
      NOPV = .TRUE.
      NDMAT = 4 ! 4 matrixes for each quater.unit.
      IATOM = 0
      NOCONT = .FALSE.
      RETUR  = .FALSE.
      RELCAL = .TRUE.
      DDFOCK = .TRUE.
      FCKDDR = .TRUE.
      LONDON = .TRUE.  
      SUSCEP = .TRUE.
      FIRST = .TRUE.

      LWRK = LWORK - KFREE + 1
      CALL DZERO(FMAT,3*NDMAT*N2BBASX)
C
C     Check what type of two-electron integrals to include
C
#if defined (VAR_MPI) || defined (VAR_PVM)
      MTOTTK(1) = NLRGBL*(NLRGBL+1)/2 ! MI...hm..some lenght...
      MTOTTK(2) = NLRGBL*(NLRGBL+1)/2
      MTOTTK(3) = NSMLBL*(NSMLBL+1)/2
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
#endif
CLV   ... temporary fix for parallel runs...
      PARCAL_SAVE = PARCAL
CMI   IF (PARCAL) PARCAL = .FALSE.
CMI   IF (PARCAL_SAVE) WRITE (LUPRI,*) 
CMI  &   " WARNING: Running in serial, parallel code not ready yet"
      CALL GETTIM(CPU1,WALL1)
      DO I2TYP = 1,3
         IF (LBIT(INTNMR,I2TYP)) THEN
#if defined (VAR_MPI) || defined (VAR_PVM)
         IF (PARCAL) THEN
CMI        CALL QUIT('Parallel suscpeptibilities not ready!')

            CALL MEMGET('INTE',KNPOS,MTOTTK(I2TYP),WORK,KFREE,LFREE)
CMI      .... replace PARDRV by HER_PARDRV and IFTHRS by SCRFCK
            CALL DIRAC_PARCTL( 1 )
            CALL HER_PARDRV(WORK(KFREE),LWRK,FMAT,DMAT,NDMAT,IREPDM,
     &            IFCTYP,ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,
     &            TKTIME,.FALSE.,FIRST,WORK(KNPOS),MTOTTK(I2TYP),
     &            I2TYP,IDUMMY,SCRFCK,DUMMY,DUMMY,DUMMY,
     &            DUMMY)
         ELSE
#endif
             CALL TWOINT(WORK(KFREE),LWRK,FMAT,DMAT,NDMAT,IREPDM,IFCTYP
     &        ,DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT
     &        ,TKTIME,IPRPRP,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY
     &        ,I2TYP,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,RELCAL,.false.
     &        ,IDUMMY,DUMMY)
#if defined (VAR_MPI) || defined (VAR_PVM)
         ENDIF
#endif
         ENDIF
      END DO
      CALL GETTIM(CPU2,WALL2)
C
      PARCAL = PARCAL_SAVE
C
      CPUTOT = CPU2 - CPU1
      WLLTOT = WALL2 - WALL1
      WRITE(LUPRI,'(/A,F13.4,A,F13.4,A/)') 
     &'Cpu time',CPUTOT,'Wall time',WLLTOT,
     &'in call to TWOINT(SUSC2EL).'
      CALL FLSHFO(LUPRI)

      CALL QEXIT('SUSC2EL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck hordrnt  */
      SUBROUTINE HORDRNT(CMO,IBEIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C   PURPOSE:
C  ----------
C   Compute the highest order expectation value
C   reorthonormalization term, {T(2), H(0)},
C   (where {T(2),H(0)} = {T(2),F(0)}) for the magnetic susceptibility 
C   with London atomic orbitals.
C   In the Levy-Leblond case this is identical
C   with the nonrelativistic DALTON counterpart.
C
C   Definition:
C -------------
C   The exact form (Ruud et al, JCP...) is -1/2{S(2),H(0)} 
C   (for both connections beacuse the S(2) operator is unique)
C
C   On input:  CMO - MO coefficients, IBEIG - boson irreps info
C
C   Called from: EXPVSUSLO (prp/pammag.F)
C
C   Written by Miro Ilias,  jun 2003,  Prievidza.
C   Last revisions:
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
      PARAMETER(D0=0.0D0, D2=2.0D0, D1 = 1.0D0)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbsusc.h"

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
      EXTERNAL BDQALGEBRA

      CALL QENTER('HORDRNT')

      KFRSAV = KFREE

C   Get the Fock matrix in MO basis

      CALL MEMGET('REAL',KFMAT,N2ORBXQ,WORK,KFREE,LFREE)
      IPRINT = IPRPRP
      CALL GETFCK(WORK(KFMAT),IPRINT,WORK,KFREE,LFREE)
      IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &'HORDRNT: Total Fock matrix in MO basis',-1)
        CALL PRQMAT(WORK(KFMAT),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')

C ... arrays for transformed F - {T(2),F(0)}_comp 
      CALL MEMGET('REAL',KTFMAT,N2ORBXQ,WORK,KFREE,LFREE)

C  ... array for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KT2MAT,N2ORBXQ,WORK,KFREE,LFREE)

      IPRINT = IPRPRP
      EPMIX = DOEPRN 
      IF (IPRINT.GE.5)  THEN
        WRITE(LUPRI,'(/,3X,A,L2 )')
     & "HORDRNT: The EPMIX variable is set to:", EPMIX
      ENDIF

C ....  DO NOT PICT.CHANGE TRANSFORMATIN OF CONNECTION MATRIXES!
      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.
 
      DO ICOMP= 1, 6 ! XX,XY,XZ,YY,YZ,ZZ
       INDXPR = IPCON23(ICOMP) ! Get indx of T(2)
       ISYM   = IPRPSYM(INDXPR)
       ITIM   = IPRPTIM(INDXPR)
       IOPSY  = JBTOF(ISYM-1,1)
       IREP   = ISYM - 1
       CALL DZERO(WORK(KT2MAT),N2ORBXQ)
C    --- Get the T(2) symm.connection matrix in MO basis ---
       IF (ITIM.EQ.0)
     & CALL QUIT('HORDRNT: The S(2) term must have itim=1!')
       CALL PRPMAT(INDXPR,IOPSY,WORK(KT2MAT),CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

C ... Call TTRA 
        AH = .FALSE. ! S(2) is hermitian ...
        ISYMT = ISYM
        ISYMH = 1 ! It's the H(0) term
        CALL TTRA(EPMIX,AH,ISYMT,ISYMH,WORK(KT2MAT),
     &    WORK(KFMAT),WORK(KTFMAT),IPRINT)

        IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &'HORDRNT: Transformed Fock matrix, {T(2),F(0)}',-1)
         WRITE(LUPRI,'(2X,A,I1)')
     & '*** tensor component(1-6): ',ICOMP
         CALL PRQMAT(WORK(KTFMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
        ENDIF

C----------------------------------------------------------------
C  Calculate the expectation value for the higher order
C reorthonormalization term, SUSFS(3,3).
C----------------------------------------------------------------
       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KTFMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO
       EXPVAL = EXPVAL*D2

C       ... control print-out ...
       IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/2X,A,I2,1X,A,D15.9)')
     &  'HORDRNT: Contribution to the highest order reorthonorm. term'
     &  //' -exp.value {S(2),H(0)} of component ',ICOMP,' is:',EXPVAL
       ENDIF

       IF (ICOMP.EQ.1) THEN
         SUSFS(IPTAX(1,2),IPTAX(1,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.2) THEN
         SUSFS(IPTAX(1,2),IPTAX(2,2)) = EXPVAL
         SUSFS(IPTAX(2,2),IPTAX(1,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.3) THEN
         SUSFS(IPTAX(1,2),IPTAX(3,2)) = EXPVAL
         SUSFS(IPTAX(3,2),IPTAX(1,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.4) THEN
         SUSFS(IPTAX(2,2),IPTAX(2,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.5) THEN
         SUSFS(IPTAX(2,2),IPTAX(3,2)) = EXPVAL
         SUSFS(IPTAX(3,2),IPTAX(2,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.6) THEN
         SUSFS(IPTAX(3,2),IPTAX(3,2)) = EXPVAL
       ELSE
         CALL QUIT('HORDRNT: ICOMP exceeds 6!')
       ENDIF

      ENDDO ! of XX,XY,XZ,YY,YZ,ZZ .....
      DOPCT=DOPCTRA_SAVE

C     .... release the whole memory block used in this subroutine
      CALL MEMREL('HORDRNT',WORK,1,KFRSAV,KFREE,LFREE)

C     ... close the AOPROPER file
      CLOSE(LU1INT,STATUS = 'KEEP')

      CALL QEXIT('HORDRNT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lordrrnt  */
      SUBROUTINE LORDRNT(CMO,IBEIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C   PURPOSE:
C  ----------
C   Compute the lowest order expectation value
C   reorthonormalization terms needed for calculating
C   the magnetic susceptibility 
C   with London atomic orbitals. 
C   In the Levy-Leblond case these are identical
C   with the DALTON counterparts.
C
C   DEFINITIONS OF TRANSFORMED OPERATORS FOR THE EXPECTATION VALUE:
C     (Ruud et. al, Chem.Phys. Vol. 195 (1995) 157,
C                 J.Chem.Phys. Vol. 99  (1993) 3847 )
C
C        Using the symmetric connection :
C -------------------------------------------------------
C   1/2{ S(1)S(1),H(0)} - {S(1),H(1)} + 1/4{S(1),S(1),H(0)}
C
C        Using the natural connection:
C -------------------------------------------------------
C   { S(1)W(1), H(0)} - {W(1),H(1)} + {W(1),W(1),H(0)}
C
C On input:  CMO - MO coefficients, IBEIG - boson irreps info
C
C   Written by Miro Ilias,  august 2003,  Odense.
C
C   Last revisions:
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
      PARAMETER(D0=0.0D0, DI2=0.5D0, DM1=-1.0D0, D2=2.0D0, D1 = 1.0D0)
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbsusc.h"

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      LOGICAL EPMIX,HERM,AH,DOPCTRA_SAVE
#include "ibtfun.h"

      CALL QENTER('LORDRNT')

      KFRSAV = KFREE

C ... allocate space for the F(0) matrix and for the MO coefficients

C   Get the Fock matrix in MO basis
      CALL MEMGET('REAL',KFMAT,N2ORBXQ,WORK,KFREE,LFREE)
      IPRINT = IPRPRP
      CALL GETFCK(WORK(KFMAT),IPRINT,WORK,KFREE,LFREE)
      IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &'LORDRNT: Total Fock matrix in MO basis',-1)
        CALL PRQMAT(WORK(KFMAT),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')

C  ... arrays for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KSWMAT,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT ,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KWMAT ,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXMAT ,N2ORBXQ,WORK,KFREE,LFREE)

C****************************************************************
C. Get  { S(1)W(1), H(0) } = { S(1)W(1), F(0) } expect. value term
C or 1/2{S(1)S(1),H(0)} = 1/2{S(1)S(1),F(0)} for symm.connection
C****************************************************************

C     WRITE(LUPRI,*)"IPTAX:",IPTAX(1,2),IPTAX(2,2),IPTAX(3,2)

      IPRINT = IPRPRP

C ... define e-p/e-shells for reorthonormalization
      EPMIX = DOEPRN
      IF (IPRINT.GE.5)  WRITE(LUPRI,'(3X,A,L2 )')
     & "LORDRNT: The EPMIX is set to:", EPMIX

C DO NOT PERFORM THE PICT.CHANGE TRANSFORMATIN OF CONNECTION MATRIXES!
      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.

      DO IS1 = 1, 3
       INDXPR1 = IPCON21(IS1) ! Get index of S(1)
       ISYM1   = IPRPSYM(INDXPR1)
       ITIM1   = IPRPTIM(INDXPR1)
       IOPSY1  = JBTOF(ISYM1-1,1)
       IREP1   = ISYM1 - 1
       CALL DZERO(WORK(KSMAT),N2ORBXQ)
C    --- Get the S(1) connection matrix in MO basis ---
       IF (ITIM1.NE.-1)
     &  CALL QUIT('LORDRNT:Wrong TRsymm. of S(1)!')
       CALL PRPMAT(INDXPR1,IOPSY1,WORK(KSMAT),CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

       DO IW1 = 1, 3
C       DO IW1 = IS1, 3
        IF (.NOT.SYMCON) THEN
C ... it's the natural connection
         INDXPR2 = IPCON22(IW1) ! Get index of W(1)
         ISYM2   = IPRPSYM(INDXPR2)
         ITIM2   = IPRPTIM(INDXPR2)
         IOPSY2  = JBTOF(ISYM2-1,1)
         IREP2   = ISYM2 - 1 ! 0-7
         CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis ---
         IF (ITIM2.NE.0) CALL QUIT('LORDRNT: Wrong TR symm. of W(1)!')
         CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ELSE
C ... it's the symmetric connection
         INDXPR2 = IPCON21(IW1) ! Get index of S(1)
         ISYM2   = IPRPSYM(INDXPR2)
         ITIM2   = IPRPTIM(INDXPR2)
         IOPSY2  = JBTOF(ISYM2-1,1)
         IREP2   = ISYM2 - 1 ! 0-7
         CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the S(1) connection matrix in MO basis ---
         IF (ITIM2.NE.-1) CALL QUIT('LORDRNT: Wrong TR symm. of S(1)!')
         CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ENDIF
C  ... make a product of operators S(1)_ir1 W(1)/S(1)_ir2 - 
C     to deterimine resulting irrep, use irreps in range 0-7,
C     otherwise use MULTD2H
C ...it's based on order of irreps in DALTON (invented by J.Almlof,1970)
C       so logical function XOR works nice...
!        Direct product of two boson irreps
         IREP12 = IBTXOR(IREP1,IREP2)
         CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &    'N','N',IPQTOQ(1,IREP1),
     &     WORK(KSMAT),NORBT,NORBT,NZ,
     &    'N','N',IPQTOQ(1,IREP2),WORK(KWMAT),NORBT,NORBT,NZ,
     &     D0,IPQTOQ(1,IREP12),WORK(KSWMAT),NORBT,NORBT,NZ)

         IF (IPRINT.GE.5) THEN
          IF (.NOT.SYMCON) THEN 
           CALL HEADER(
     &'LORDRNT: Transformation matrix product, S(1)*W(1)',-1)
          ELSE
           CALL HEADER(
     &'LORDRNT: Transformation matrix product, S(1)*S(1)',-1)
          ENDIF
        WRITE(LUPRI,'(2X,A,2I2)')
     & '*** tensor components: ',IS1,IW1
           CALL PRQMAT(WORK(KSWMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
         ENDIF

C ... call TTRA for {S(1)W(1), H(0)} = {S(1)W(1), F(0)}
C ... call TTRA for {S(1)S(1), H(0)} = {S(1)S(1), F(0)}

        AH  = .FALSE. ! product of S(1)W(1) is hermitian 
        ISYM12 = IREP12 + 1 ! Total symmetry of S(1)*W(1)/S(1) operators
        ISYMF=1 ! F(0) is fully symmetric
        CALL TTRA(EPMIX,AH,ISYM12,ISYMF,WORK(KSWMAT),
     &             WORK(KFMAT),WORK(KXMAT),IPRINT)

        IF (IPRINT.GE.5) THEN
         IF (.NOT.SYMCON) THEN 
          CALL HEADER(
     &'LORDRNT: Transformed Fock matrix, {S(1)W(1),F(0)}',-1)
         ELSE
          CALL HEADER(
     &'LORDRNT: Transformed Fock matrix, {S(1)S(1),F(0)}',-1)
         ENDIF
        WRITE(LUPRI,'(2X,A,2I2)')
     &'*** tensor component: ',IS1,IW1
         CALL PRQMAT(WORK(KXMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
        ENDIF

C  ... save the expectation values for the lower order rn term
C      into SUSFSY(3,3)

       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KXMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO

       IF (.NOT.SYMCON) THEN
        EXPVAL = -EXPVAL*D2 ! Because of -S(1)W(1)
       ELSE
        EXPVAL = -EXPVAL ! Take care of -1/2 S(1)S(1)
       ENDIF

       IF (IPRINT.GE.5) THEN
        IF (.NOT.SYMCON) THEN  
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {S(1)W(1),H(0)} of component ',IS1,IW1,EXPVAL
        ELSE
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {S(1)S(1),H(0)} of component ',IS1,IW1,EXPVAL
        ENDIF
       ENDIF

       SUSFSY(IPTAX(IS1,2),IPTAX(IW1,2)) = EXPVAL

       ENDDO
      ENDDO

C===================================================
C   Proceed to the second term,
C -------------------------------------------------
C    - {W(1)/S(1),H(1)} = {S(1)/W(1),F(1)},
C           where F(1)=h(1)+u(1)...
C===================================================

C     ====   Get h(1) and 2*u(1) MO terms ====

      DO IW1 = 1, 3  ! Components of W(1)
       IF (.NOT.SYMCON) THEN
        INDXPR2 = IPCON22(IW1) ! Get index of +W(1)
        ISYM2   = IPRPSYM(INDXPR2)
        ITIM2   = IPRPTIM(INDXPR2)
        IOPSY2  = JBTOF(ISYM2-1,1)
        IREP2   = ISYM2 - 1
        CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis into WORK(KWMAT) ---
        IF (ITIM2.NE.0) CALL QUIT('LORDRNT:Wrong TR of W(1)!')
        CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),CMO,IBEIG,
     &                 ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
      ELSE
        INDXPR2 = IPCON21(IW1) ! Get index of +S(1)
        ISYM2   = IPRPSYM(INDXPR2)
        ITIM2   = IPRPTIM(INDXPR2)
        IOPSY2  = JBTOF(ISYM2-1,1)
        IREP2   = ISYM2 - 1
        CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis into WORK(KWMAT) ---
        IF (ITIM2.NE.-1) CALL QUIT('LORDRNT:Wron TR of s(1)!')
        CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),CMO,IBEIG,
     &                 ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
      ENDIF
C ... get components of F(1)=h(1)+u(1) in MO basis into WORK(KSMAT)
        DO IH1 = 1, 3 ! Components of H(1)
C       DO IH1 = IW1, 3 ! Components of H(1)
          IB = 3 ! 0011 ... get ONELON + TWOLON
          INDXPR1 = IPLONDON(1,IH1) ! RM1H1_comp
          ISYM1 = IPRPSYM(INDXPR1) 
          IREP1 = ISYM1 - 1
          IREP12 = IBTXOR(IREP2,IREP1) 
          CALL DZERO(WORK(KSMAT),N2ORBXQ)
          CALL GETDHB(IB,D0,WORK(KSMAT),
     &           INDXPR1,IOPSY,WORK,KFREE,LFREE,IPRINT)

C ... call TTRA for {W/S(1),F(1)} .. get all into KSWMAT
          AH  = .TRUE. ! T=W =>evidently antisymm...
          CALL TTRA(EPMIX,AH,ISYM2,ISYM1,WORK(KWMAT),
     &            WORK(KSMAT),WORK(KSWMAT),IPRINT)

          IF (IPRINT.GE.5) THEN
           IF (.NOT.SYMCON) THEN
            CALL HEADER(
     &     'LORDRNT: Transformed Fock-LAO matrix, {W(1),F(1)}',-1)
           ELSE
            CALL HEADER(
     &     'LORDRNT: Transformed Fock-LAO matrix, {S(1),F(1)}',-1)
           ENDIF
           WRITE(LUPRI,'(2X,A,2I2)')
     & '*** tensor component: ',IW1,IH1
           CALL PRQMAT(WORK(KSWMAT),NORBT,NORBT,NORBT,
     &                 NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
          ENDIF

C  ... save the expectation values for the lower order rn term
C      into SUSFSY(3,3)

       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KSWMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO

       IF (.NOT.SYMCON) THEN
         EXPVAL = D2*D2*EXPVAL ! Times 2*
       ELSE
         EXPVAL = D2*EXPVAL ! Standard..
       ENDIF

       IF (IPRINT.GE.5) THEN
        IF (.NOT.SYMCON) THEN
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT:Contribution to the lowest order reorthonorm. term'
     &//' {W(1),H(1)} of component ',IW1,IH1,EXPVAL
        ELSE
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT:Contribution to the lowest order reorthonorm. term'
     &//' {S(1),H(1)} of component ',IW1,IH1,EXPVAL
        ENDIF
       ENDIF

        SUSFSY(IPTAX(IW1,2),IPTAX(IH1,2)) = 
     &  SUSFSY(IPTAX(IW1,2),IPTAX(IH1,2)) + EXPVAL
         
      ENDDO

      ENDDO

C=================================================================
C
C Proceed to the third term, {W(1),W(1),H(0)}={W(1),{W(1),H(0)}}
C
C Get F'=({W(1),F(0)}+W) contribution and  then make {W(1),F'}
C
C The same for S(1)
C
C=================================================================

      DO I1  =   1, 3
       IF (.NOT.SYMCON) THEN
C ...  get W(1)_icomp(I1) into WORK(KSMAT)
        INDXPR1 = IPCON22(I1) ! Get index of +W(1)
        ISYM1   = IPRPSYM(INDXPR1)
        ITIM1   = IPRPTIM(INDXPR1)
        IOPSY1  = JBTOF(ISYM1-1,1)
        IREP1   = ISYM1 - 1 ! 0-7
        CALL DZERO(WORK(KSMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis into WORK(KSMAT) ---
        IF (ITIM1.NE.0) CALL QUIT('LORDRNT: Wrong TR symm. of W(1)!')
        CALL PRPMAT(INDXPR1,IOPSY1,WORK(KSMAT),CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ELSE
C ...  get S(1)_icomp(I1) into WORK(KSMAT)
        INDXPR1 = IPCON21(I1) ! Get index of +S(1)
        ISYM1   = IPRPSYM(INDXPR1)
        ITIM1   = IPRPTIM(INDXPR1)
        IOPSY1  = JBTOF(ISYM1-1,1)
        IREP1   = ISYM1 - 1 ! 0-7
        CALL DZERO(WORK(KSMAT),N2ORBXQ)
C    --- Get the S(1) connection matrix in MO basis into WORK(KSMAT) ---
        IF (ITIM1.NE.-1) CALL QUIT('LORDRNT: Wrong TR symm. of S(1)!')
        CALL PRPMAT(INDXPR1,IOPSY1,WORK(KSMAT),CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ENDIF

C       DO I2  =  I1, 3
        DO I2  =  1, 3
C ...  get FCKLON_icomp(I2) into WORK(KWMAT)
         INDXPR2 = IPLONDON(1,I2) ! Get index of RM1H1_comp
         ISYM2   = IPRPSYM(INDXPR2)
         IREP2   = ISYM2 - 1 ! 0-7
         IREP12  = IBTXOR(IREP1,IREP2) ! Final symm.
         CALL DZERO(WORK(KWMAT),N2ORBXQ)
         IB = 4 ! 0100B -> get F'=FCKLON... {T(1),F(0)}+W into WORK(KWMAT)
         CALL GETDHB(IB,D0,WORK(KWMAT),
     &          INDXPR2,IOPSY,WORK,KFREE,LFREE,IPRINT)

C ...  Do transform  F''={W(1)/S(1)_ic1,F'_ic2} into WORK(KSWMAT)
         AH  = .TRUE. ! T=W=>evidently antiherm...
         CALL DZERO(WORK(KSWMAT),N2ORBXQ)
         CALL TTRA(EPMIX,AH,ISYM1,ISYM2,WORK(KSMAT),
     &          WORK(KWMAT),WORK(KSWMAT),IPRINT)

         IF (IPRINT.GE.5) THEN
          IF (.NOT.SYMCON) THEN
           CALL HEADER(
     & 'LORDRNT: TOTAL transformed matrix {W(1),W(1),H(0)}',-1)
          ELSE
           CALL HEADER(
     & 'LORDRNT: TOTAL transformed matrix {S(1),S(1),H(0)}',-1)
          ENDIF
          WRITE(LUPRI,'(2X,A,2I2)')
     & '*** tensor component: ',I1,I2
          CALL PRQMAT(WORK(KSWMAT),NORBT,NORBT,NORBT,
     &               NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
         ENDIF

C ...  Store expect. values from resulting matrix WORK(KSWMAT)

       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KSWMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO

       IF (.NOT.SYMCON) THEN
        EXPVAL = D2*EXPVAL
       ELSE
        CONTINUE  ! Consider the 1/2
       ENDIF

       IF (IPRINT.GE.5) THEN
        IF (.NOT.SYMCON) THEN
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     &'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {W(1),W(1),H(0)} of component ',I1,I2,EXPVAL
        ELSE
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     &'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {S(1),S(1),H(0)} of component ',I1,I2,EXPVAL
        ENDIF
       ENDIF

       SUSFSY(IPTAX(I1,2),IPTAX(I2,2)) = 
     &  SUSFSY(IPTAX(I1,2),IPTAX(I2,2)) + EXPVAL

      ENDDO
      ENDDO

C Symmetrize the off-diagonal terms of SUSFSY !!

      DO I1 = 1, 2
      DO I2 = I1+1, 3
       PS =
     & (SUSFSY(IPTAX(I1,2),IPTAX(I2,2))+
     &  SUSFSY(IPTAX(I2,2),IPTAX(I1,2)))/D2
       SUSFSY(IPTAX(I1,2),IPTAX(I2,2))=PS
       SUSFSY(IPTAX(I2,2),IPTAX(I1,2))=PS
      ENDDO
      ENDDO

C     .... release the whole block used in this subroutine
      CALL MEMREL('LORDRNT',WORK,1,KFRSAV,KFREE,LFREE)

      CLOSE(LU1INT,STATUS = 'KEEP')

      DOPCT=DOPCTRA_SAVE

      CALL QEXIT('LORDRNT')
C     CALL QUIT('Exit from LORDRNT')
      
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_d1hblond */
       SUBROUTINE DEF_D1HBLOND(IMP,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C   PURPOSE:
C ===============
C     Define all one-electron operators
C     for magnetic property
C     gradient dH/dB|B=0 (A,B,or both)
C     when using London orbitals for the specified magnetic dependent
C     property IMP:
C
C     IMP = 1 .... NMR shielding, optical rotation (TEC)
C     IMP = 2 .... magnetic susceptibility (MI)
C     IMP = 3 ...  perturbation from one component of magnetic field
C                  (gosia: temporarily i need it in visual module)
C     IMP = 4 ...  and more (to be added later)
C
C  Written by:  MI, jan.2003, based on TEC's routines
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER(DI2 = +0.5D00, DMI2 = -0.5D0, DM1 = -1.0D0,
     &          D2 = 2.0D0, D1 = 1.0D0)
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbsusc.h"
#include "dcbnmr.h"
      CHARACTER PNAME*16, PLABEL(3)*8, LABELC*7
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)

      CALL QENTER('DEF_D1HBLOND')
#include "memint.h"
      IPRINT = IPRPRP

      IF (IMP.EQ.1) THEN
        IF(SHIELD) CALL HEADER(
     &'DEF_D1HBLOND: Operators for NMR shieldings with LAO!',-1)
        IF(OPTROT) CALL HEADER(
     &'DEF_D1HBLOND: Operators for optical rotation with LAO!',-1)
      ELSE IF (IMP.EQ.2) THEN
       CALL HEADER(
     &'DEF_D1HBLOND: Operators for magnetic'//
     &' susceptibility with LAO!',-1)
      ELSE IF (IMP.EQ.3) THEN
        CALL HEADER(
     &'DEF_D1HBLOND: Operators for perturbation from '//
     &' one component of magnetic field with LAO!',-1)
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & 'DEF_D1HBLOND: Entering magnetic property IMP '//
     & ' has to be 1 or 2 or 3!'
       CALL QUIT('DEF_D1HBLOND: wrong IMP value!')
      ENDIF

C------------------------------------------------
C          X component of the B-field
C------------------------------------------------
C     RM1H1 operator ... i Q_MN <r beta' c^2 > 
C
      PNAME    = 'LAO-XRM1H1'
      IPTYP    = 1
      NPCOMP   = 1
      PFAC(1)  = -D2*CVAL*CVAL
      PLABEL(1)= 'dS/dBX'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(1,1) = INDXPR

      IF (IMP.EQ.1) THEN
        CALL OP1IND('XSHIELB',IPSHIELD(1+MXCOOR),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF (IMP.EQ.2) THEN
        CALL OP1IND('MSUSCX1',IPMSUSC(1,1),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('MSUSCX2',IPMSUSC(2,1),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF ((IMP.EQ.3) .AND. BXLAO) THEN
        CALL OP1IND('BLAOX1',IP_BLAO(1,1),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('BLAOX2',IP_BLAO(2,1),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE
        CONTINUE
      ENDIF

C--------------------------------------------------------------
C       RM1H2 operator ... i Q_MN < r_vec Vnuc >
C--------------------------------------------------------------
      PNAME    = 'LAO-XRM1H2'
      IPTYP    = 1
      NPCOMP   = 1
      PFAC(1)  = D1
      PLABEL(1)= 'XRM1H2'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(2,1) = INDXPR
C--------------------------------------------------------------
C       RM1H3  operator ... i c Q_MN < r_vec (alp.p) >
C--------------------------------------------------------------
      PNAME     = 'LAO-XRM1H3'
      IPTYP     = 8
      NPCOMP    = 3
      PFAC(1)  = CVAL
      PFAC(2)  = CVAL
      PFAC(3)  = CVAL
      PLABEL(1) = 'XRM1H3X '
      PLABEL(2) = 'XRM1H3Y '
      PLABEL(3) = 'XRM1H3Z '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(3,1) = INDXPR
C
C  RM1RN-X operator of +(c/2).< r_n x alpha >
C
      PNAME    = 'LAO-XRM1RN'
      IPTYP    = 5
      NPCOMP   = 2
      PFAC(1)  = 0.5d0*cval
      PFAC(2)  = 0.5d0*cval
      PLABEL(1) = 'ZRM1RN'
      PLABEL(2) = 'YRM1RN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(4,1) = INDXPR

#if !defined (XH4)
C-------------------------------------------------------------
C            Y component of B-field
C-------------------------------------------------------------
C
C         RM1H1-Y operators
C
      PNAME     = 'LAO-YRM1H1'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = -D2*CVAL*CVAL
      PLABEL(1) = 'dS/dBY'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &       INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(1,2) = INDXPR

      IF (IMP.EQ.1) THEN
       CALL OP1IND('XSHIELB',IPSHIELD(2+MXCOOR),LLRBPU,NLRBPT,INDXPR
     &             ,MAXLLR)
      ELSE IF (IMP.EQ.2) THEN
       CALL OP1IND('MSUSCY1',IPMSUSC(1,2),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
       CALL OP1IND('MSUSCY2',IPMSUSC(2,2),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF ((IMP.EQ.3) .AND. BYLAO) THEN
        CALL OP1IND('BLAOY1',IP_BLAO(1,2),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('BLAOY2',IP_BLAO(2,2),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE
       CONTINUE
      ENDIF
C
C        RM1H2-Y operator
C
      PNAME     = 'LAO-YRM1H2'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YRM1H2'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(2,2) = INDXPR
C
C       RM1H3-Y operator
C
      PNAME    = 'LAO-YRM1H3'
      IPTYP    = 8
      NPCOMP   = 3
      PFAC(1)  = CVAL
      PFAC(2)  = CVAL
      PFAC(3)  = CVAL
      PLABEL(1) = 'YRM1H3X '
      PLABEL(2) = 'YRM1H3Y '
      PLABEL(3) = 'YRM1H3Z '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(3,2) = INDXPR
C
C          RM1RN-Y operator
C
      PNAME     = 'LAO-YRM1RN'
      IPTYP     = 6
      NPCOMP    = 2
      PFAC(1)   = 0.5d0*cval
      PFAC(2)   = 0.5d0*cval
      PLABEL(1) = 'XRM1RN'
      PLABEL(2) = 'ZRM1RN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(4,2) = INDXPR

C-------------------------------------------------------------
C     Z component of B-field
C-------------------------------------------------------------
C
C         RM1H1-Z operator
C
      PNAME     = 'LAO-ZRM1H1'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = -D2*CVAL*CVAL
      PLABEL(1) = 'dS/dBZ'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &      INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(1,3) = INDXPR

      IF (IMP.EQ.1) THEN
        CALL OP1IND('XSHIELB',IPSHIELD(3+MXCOOR),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF (IMP.EQ.2) THEN
        CALL OP1IND('MSUSCZ1',IPMSUSC(1,3),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('MSUSCZ2',IPMSUSC(2,3),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF ((IMP.EQ.3) .AND. BZLAO) THEN
        CALL OP1IND('BLAOZ1',IP_BLAO(1,3),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('BLAOZ2',IP_BLAO(2,3),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE
       CONTINUE
      ENDIF
C
C         RM1H2-Z operator
C
      PNAME     = 'LAO-ZRM1H2'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'ZRM1H2'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(2,3) = INDXPR
C
C         RM1H3-Z operator
C
      PNAME     = 'LAO-ZRM1H3'
      IPTYP     = 8
      NPCOMP    = 3
      PFAC(1)   = CVAL
      PFAC(2)   = CVAL
      PFAC(3)   = CVAL
      PLABEL(1) = 'ZRM1H3X '
      PLABEL(2) = 'ZRM1H3Y '
      PLABEL(3) = 'ZRM1H3Z '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &     INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(3,3) = INDXPR
C
C        RM1RN-Z operator of +(c/2).< r_n x alpha >
C
      PNAME     = 'LAO-ZRM1RN'
      IPTYP     = 7
      NPCOMP    = 2
      PFAC(1)   = 0.5d0*cval
      PFAC(2)   = 0.5d0*cval
      PLABEL(1) = 'YRM1RN'
      PLABEL(2) = 'XRM1RN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPLONDON(4,3) = INDXPR
#endif

C----------------------------------------------------------------
C   Define connection matrix T for the first derivate dH/dB
C----------------------------------------------------------------
! for imp=3 always write the correct indxpr to ipcon(1), it is easier for a moment
      IPTYP  = 1
      NPCOMP = 1
      IF (SYMCON) THEN
        PFAC(1) = DMI2
        LABELC='dS/dB'
        ILABELC = 5
      ELSE
        PFAC(1) = -D1
        LABELC='d|S>/dB'
        ILABELC = 7
      END IF
      PNAME     = LABELC(1:ILABELC)//'X'
      PLABEL(1) = LABELC(1:ILABELC)//'X'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPCON(1) = INDXPR
      if ((imp .eq. 3) .and. bxlao) then
        ipcon(1) = indxpr
      end if
      PNAME     = LABELC(1:ILABELC)//'Y'
      PLABEL(1) = LABELC(1:ILABELC)//'Y'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPCON(2) = INDXPR
      if ((imp .eq. 3) .and. bylao) then
        ipcon(1) = indxpr
      end if
      PNAME     = LABELC(1:ILABELC)//'Z'
      PLABEL(1) = LABELC(1:ILABELC)//'Z'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPCON(3) = INDXPR
      if ((imp .eq. 3) .and. bzlao) then
        ipcon(1) = indxpr
      end if

! if (imp == 3) only ipcon(1) is important, lines below are to catch errors 
      if (imp .eq. 3) then
        ipcon(2) = 0
        ipcon(3) = 0
      end if

C----------------------------------------------------------------
C   Define connection matrixes T(1),T(2) for 
C   the second derivate d2H/dB2 
C----------------------------------------------------------------
      IF (IMP.EQ.2) THEN

C ... Symmetry connection term T(1)=S(1)
       IPTYP  = 1
       NPCOMP = 1
C      PFAC(1) = DMI2
       PFAC(1) = D1
       LABELC='dS/dB'
       ILABELC = 5
       PNAME     = LABELC(1:ILABELC)//'X'
       PLABEL(1) = LABELC(1:ILABELC)//'X'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON21(1) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'Y'
       PLABEL(1) = LABELC(1:ILABELC)//'Y'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON21(2) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'Z'
       PLABEL(1) = LABELC(1:ILABELC)//'Z'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON21(3) = INDXPR

C ... Natural connection term T(1)=W(1)
C      PFAC(1) = -D1
       PFAC(1) =  D1
       LABELC='d|S>/dB'
       ILABELC = 7
       PNAME     = LABELC(1:ILABELC)//'X'
       PLABEL(1) = LABELC(1:ILABELC)//'X'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON22(1) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'Y'
       PLABEL(1) = LABELC(1:ILABELC)//'Y'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON22(2) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'Z'
       PLABEL(1) = LABELC(1:ILABELC)//'Z'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON22(3) = INDXPR

C Define d2S/dBXX,XY,...
C Specify latter PFAC(1)!!!
       IPTYP  = 1
       NPCOMP = 1
       PFAC(1) = DMI2 ! For the  "-1/2 {S(2),H(0)}" term
       LABELC='dS/dB2'
       ILABELC =6 
       PNAME     = LABELC(1:ILABELC)//'XX'
       PLABEL(1) = LABELC(1:ILABELC)//'XX'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON23(1) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'XY'
       PLABEL(1) = LABELC(1:ILABELC)//'XY'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON23(2) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'XZ'
       PLABEL(1) = LABELC(1:ILABELC)//'XZ'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON23(3) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'YY'
       PLABEL(1) = LABELC(1:ILABELC)//'YY'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON23(4) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'YZ'
       PLABEL(1) = LABELC(1:ILABELC)//'YZ'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON23(5) = INDXPR
       PNAME     = LABELC(1:ILABELC)//'ZZ'
       PLABEL(1) = LABELC(1:ILABELC)//'ZZ'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
       IPCON23(6) = INDXPR
      ENDIF

      CALL QEXIT('DEF_D1HBLOND')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sum_d1hblo */
       SUBROUTINE GETDHB(IB,FAC2LO,PRPMO,INDXPR,
     &              IOPSY,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C    PURPOSE:
C  ==============
C   Accumulate selected terms (depends on IB)
C   of the (dH/dB)|B=0 property gradient
C   in MO-basis into the enetring PRPMO MO matrix
C   when using London orbitals.
C
C   The gradient is needed for calculating the NMR shieldings and
C   the magnetizabilities.
C
C *** Terms:
C
C       dH/dB|B=0|_c = ONELON + FCKLON + (FAC2LO*)TWOLON
C
C  where ONELON = (RM1RN(sq) + RM1H1+ RM1H2 + RM1H3(sq))
C
C  Bit flags:
C ---------------
C
C    1.bit(1)  -  ONELON
C    2.bit(2)  -  TWOLON
C    3.bit(4)  -  FCKLON
C    4.bit(8)  -  FAC2LO*TWOLON  (+ 2.bit ON!)
C
C  Examples:
C    IB = 0111 (7)  - get all ONELON + FCKLON + TWOLON (full mag.field.grad.)
C    IB = 1011 (11) - get ONELON + FAC2LO*TWOLON (FAC2LO is entering)
C    IB = 1111 (15) - get ONELON + FCKLON + FAC2LO*TWOLON (FAC2LO is entering)
C    IB = 0001 (1)  - get only ONELON 
C    IB = 0101 (5)  - get only ONELON + FCKLON
C
C  The entering INDXPR for this case corresponds to some component of the
C  RM1H1 operator (X/Y/Z). 
C
C  The resulting PRPMO_c OUTPUT matrix is TR antisymmetric, what is 
C  is checked inside this routine.
C
C   Called from: GOPGET
C                GRDOR1
C
C   Written by Miro Ilias, August 2003, Odense
C 
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,  D1 = 1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
C
#include "mxcent.h"
#include "dcbnmr.h"
C
!gosia, i need bxlao, bylao, bzlao
#include "dcbprp.h"
      DIMENSION PRPMO(NORBT,NORBT,NZ),WORK(*)
      INTEGER ICONVXYZ
      CHARACTER*1 COMP
      DATA THRS / 1.0D-7 /
      DATA LUONELO / 72 /
#include "ibtfun.h"

      CALL QENTER('GETDHB')

      KFRSAV = KFREE

      CALL DZERO(PRPMO, N2ORBXQ)

C     ... PRPNAM(INDXPR) is the antisymmetric 'LAO-cRM1H1'
      IF (PRPNAM(INDXPR)(6:10).NE.'RM1H1') THEN
       WRITE(LUPRI,'(2X,A,I3,A)')
     & 'GETDHB: Entering INDXPR, PRPNAM(INDXPR):',  
     &   INDXPR, PRPNAM(INDXPR)
       CALL QUIT('GETDHB: Entering wrong INDXPR!')
      ENDIF
      READ (PRPNAM(INDXPR),'(4X,A1)') COMP
      ICOMP = ICONVXYZ(COMP)

      IREP = IPRPSYM(INDXPR)-1

      CALL MEMGET('REAL',KPRPMOLO,N2ORBXQ,WORK,KFREE,LFREE)

C     WRITE(LUPRI,*) 'IB=',IB
C     WRITE(LUPRI,*) 'IBTAND(IB,1):',IBTAND(IB,1)
C     WRITE(LUPRI,*) 'IBTAND(IB,2):',IBTAND(IB,2)
C     WRITE(LUPRI,*) 'IBTAND(IB,4):',IBTAND(IB,4)
C     WRITE(LUPRI,*) 'IBTAND(IB,8):',IBTAND(IB,8)

C   ... 1. accumulate ONELON into WORK(KPRPMOLO) if 1.bit is ON
      IF (IBTAND(IB,1).EQ.1) THEN
         OPEN(LUONELO,FILE='ONELON',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='OLD')
         CALL READAC(LUONELO,N2ORBXQ,WORK(KPRPMOLO),ICOMP)
         CLOSE(LUONELO,STATUS = 'KEEP')
         IF (IBTAND(IB,8).NE.8) THEN
          CALL DAXPY(N2ORBXQ,D1,WORK(KPRPMOLO),1,PRPMO,1)
         ELSE
           WRITE(LUPRI,'(2X,A,D12.7)')
     &    'GETDHBLO: Added ONELON with FAC2LO=',FAC2LO
          CALL DAXPY(N2ORBXQ,FAC2LO,WORK(KPRPMOLO),1,PRPMO,1)
         ENDIF
         CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A)')'GETDHBLO: Added ONELON'
         ENDIF
      ENDIF

C   ... 2. accumulate TWOLON into WORK(KPRPMOLO) if 2.bit is ON
      IF (IBTAND(IB,2).EQ.2.AND.(.NOT.NOTWOL)) THEN
         OPEN(LUTWOL,FILE='TWOLON',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='OLD')
         !CALL READAC(LUTWOL,N2ORBXQ,WORK(KPRPMOLO),ICOMP)
         CALL READAC(LUTWOL,N2ORBXQ,WORK(KPRPMOLO),ICOMP)
C   ...   if 4.bit is ON, acc TWOLON with some factor FAC2LO
         IF (IBTAND(IB,8).EQ.8) THEN
          CALL DAXPY(N2ORBXQ,FAC2LO,WORK(KPRPMOLO),1,PRPMO,1)
         ELSE
          CALL DAXPY(N2ORBXQ,D1,WORK(KPRPMOLO),1,PRPMO,1)
         ENDIF
         CLOSE(LUTWOL,STATUS = 'KEEP')
         CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A)')'GETDHBLO: Added TWOLON'
          IF (IBTAND(IB,8).EQ.8) THEN
           WRITE(LUPRI,'(2X,A,D12.7)')
     &    'GETDHBLO: Added TWOLON with FAC2LO=',FAC2LO
          ENDIF
         ENDIF
      ENDIF

C   ... 3. accumulate FCKLON into WORK(KPRPMOLO) if 3.bit is ON
!gosia (i'll fix it later)
      if (bxlao .or. bylao .or. bzlao) then
!       we always write to component "1"
        icomp_reorth = 1
      else
        icomp_reorth = icomp
      end if
      IF (IBTAND(IB,4).EQ.4.AND.(.NOT.NOONEI)) THEN
         OPEN(LUFCKL,FILE='FCKLON',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='OLD')
         CALL READAC(LUFCKL,N2ORBXQ,WORK(KPRPMOLO),icomp_reorth)
         CALL DAXPY(N2ORBXQ,D1,WORK(KPRPMOLO),1,PRPMO,1)
         CLOSE(LUFCKL,STATUS = 'KEEP')
         CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A)')'GETDHBLO: Added FCKLON'
         ENDIF
      ENDIF

C ... print the total matrix and CHECK if it is anti-symmetric !

        IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/A,A,A,I3)')
     & 'GETDHBLO: For *** LAO-PDBS ***'
     &//' the total accumulated MO-matrix of property ',
     &  PRPNAM(INDXPR), 'with flag',IB
        CALL PRQMAT(PRPMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
        ENDIF

C ... check the symmetry of the total property
C      gradient matrix, PRPMO for LAO-PDBS

      DSX=D0
      DO IZ=1,NZ
       DO I=1,NORBT-1
       DO J=I+1,NORBT
C   ... accumulate sum of differences
        DSX = DSX +
     & (DABS(PRPMO(I,J,IZ))-DABS(PRPMO(J,I,IZ)) )
       ENDDO
       ENDDO
       ENDDO
       DSX = DSX/DFLOAT( ((NORBT*NORBT-NORBT)/2)*NZ)

       IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(4X,A,/,A,A,A,1X,D12.5)')
     &'*** GETDHB: The average symmetry of the final dH/dB',
     &' MO-property matrix ',PRPNAM(INDXPR),
     &' with LAO is: ',DSX
       ENDIF

      IF (DABS(DSX).GE.THRS) THEN
        WRITE(LUPRI,'(4X,A,/,A,A,A,1X,D12.5)')
     &'*** GETDHB: The average symmetry of the final dH/dB',
     &' MO-property matrix ',PRPNAM(INDXPR),
     &' with LAO is: ',DSX

        WRITE(LUPRI,'(A)')
     &'GETDHB:: The total propery gradient'
     &//'matrix dH/dB for LAO-PDBS is not (fully?) symmetric !'
CMI     CALL QUIT(
CMI  &'GETDHB: The total MO matrix in not symmetric!')
      ENDIF

      CALL MEMREL('GETDHB',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETDHB')
      RETURN
      END
#endif /* MOD_LONDON */

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

      RETURN
      END

