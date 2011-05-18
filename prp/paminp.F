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

C VERSION : $Revision$
C DATE    : $Date$
C FILE    : paminp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpinp */
      SUBROUTINE PRPINP(WRDSRC,WORK,LWORK)
C***********************************************************************
C
C <<< Property Input for DIRAC >>>
C
C
C    Called from:  PAMINP - first reading of the input
C                  MINPRP - called when optimize is on 
C                  (properties at final geometry)
C
C       Written by Trond Saue and Hans Joergen Aa.Jensen
C       Last revision: June 14 1996 - tsaue
C
C***********************************************************************

         use dirac_cfg
         use input_reader
#ifdef MOD_OPENRSP
         use dirac_interface
         use openrsp_cfg
         use openrsp_interface_1el
#endif

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 10, NTABLE = 40)
#include "maxorb.h"
#include "mxcent.h"
#include "dcbham.h"
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7, WRDSRC*7
      DIMENSION WORK(LWORK)
      LOGICAL CHGO
C
#include "orgcom.h"
#include "dcbprp.h"
#include "dgroup.h"
#include "xmlout.h"
C
      DATA TABDIR /'*END OF','*EXPECT','*LINEAR','*EXCITA',
     &             '*MOLGRD','*NMR   ','*ESR   ','*QUADRA',
     &             '*STEX  ','*OPENRS'/
!                     1         2         3         4
      DATA TABLE  /'.DIPOLE','.QUADRU','.POLARI','.MAGNET',
!                     5         6         7         8
     &             '.NMR   ','.SPIN-S','.SHIELD','.PRINT ',
!                     9        10        11        12
     &             '.ABUNDA','.EFG   ','.NQCC  ','.DSO   ',
!                    13        14        15        16
     &             '.MOLGRD','.NSTDIA','.ESR   ','.PVC   ',
!                    17        18        19        20
     &             '.EFT   ','.EFTNTL','.PVCSHI','.EFN   ',
!                    21        22        23        24
     &             '.WAVE F','.RHONUC','.FIRST ','.TWO-PH',
!                    25        26        27        28
     &             '.VERDET','.ROTG  ','.MSTDIA','.PVCSPN',
!                    29        30        31        32
     &             '.NOPCTR','.EXCDIP','.MAGOUT','.OPTROT', 
!                    33        34        35        36
     &             '.C6    ','.RKBIMP','.BXLAO ','.BYLAO ',
!                    37        38        39        40
     &             '.BZLAO ','.EFFDEN','.XMLOUT','.XXXXXX' /

C
      CALL QENTER('PRPINP')
#include "memint.h"
C
      NJOBS = 0
C
C     Initialize /CBIPRP/
C     ===================
C
      IPRPRP       = 0
      DOEXP        = .FALSE.
      DOESR        = .FALSE.
      DOXLR        = .FALSE.
      DOXPP        = .FALSE.
      DOXQR        = .FALSE.
      DOTPA        = .FALSE.
      DOVER        = .FALSE.
      DOEXCPRP     = .FALSE.
      EXCPRPTYP    = 0
      DIPOLE       = .FALSE.
      QUADRU       = .FALSE.
      EPOLAR       = .FALSE.
      DOVDW        = .FALSE.
      MSUSCP       = .FALSE.
      IMAGOUT      = 0
      SHIELD       = .FALSE.
      PVC_SHIELD   = .FALSE.
      PVC_SPINSPIN = .FALSE.
      SPNSPN       = .FALSE.
      EFG          = .FALSE.
      EFN          = .FALSE.
      EFT          = .FALSE.
      EFTNTL       = .FALSE.
      OPTROT       = .FALSE.
      NQCC         = .FALSE.
      DSO          = .FALSE.
      NSTDIA       = .FALSE.
      MSUSCDIA     = .FALSE.
      ROTG         = .FALSE.
      LONDON       = .FALSE.
      PVC          = .FALSE.
      RHONUC       = .FALSE.
      ABUND        = D1
      NDSOINT      = 0
      MGRAD        = .FALSE.
      ESRGTENS     = .FALSE.
      ESRANISO     = .FALSE.
      NPRPWF       = 1
      PRPWF(1)     = 'DHF'
      DORKBIMP     = .FALSE.
      BXLAO        = .FALSE.
      BYLAO        = .FALSE.
      BZLAO        = .FALSE.
      EFFDEN       = .FALSE.


CMI Activate the  picture change transformations
C   of four-component property matrixes 
C   for the two-component (BSS) infinite order.
C TODO with TS/HJJ: BSS spin-free picture change transformation
      NOPCT = .FALSE.
      IF (BSS.AND.(IBSS.GE.99.OR.IBSS.EQ.009.OR.IBSS.GE.109)) THEN
        DOPCT = .TRUE.
C       ... by default DO the picture change transformation
        NOPCT = .FALSE.
        WRITE(LUPRI,'(/3X,A/)')
     &'PRPINP: Picture change transformation of properties activated !'
      ELSE
        DOPCT = .FALSE.
      ENDIF
C
C     Read menu file
C     ==============
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      IF (WORD .EQ. WRDSRC) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
C        Tell user which property input is missing
C        (**PROPE or **PRP F or ?)
         WRITE(LUPRI,'(/3A/A)')
     &   ' * INFORMATION: No property input "',WRDSRC,
     &   '" found in input file.',
     &   ' - No properties calculated.'
         CALL QEXIT('PRPINP')
         RETURN
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no property input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /CBIPRP/
C     ==================================
C
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,
     &                109,110,111,112,113,114,115,116,
     &                117,118,119,120,121,122,123,124,
     &                125,126,127,128,129,130,131,132,
     &                133,134,135,136,137,138,139,140), I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in PRPINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in PRPINP.')
  101    CONTINUE
C&&&& DIPOLE - dipole moment
            DIPOLE   = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  102    CONTINUE
C&&&& QUADRUpole moment
            DOEXP   = .TRUE.
            QUADRU  = .TRUE.
            GO TO 100
  103    CONTINUE
C&&&& POLAR - electronic dipole polarizabilities
            DOXLR  = .TRUE.
            EPOLAR = .TRUE.
            GO TO 100
  104    CONTINUE
C&&&& MAGNET - magnetic susceptibilities
#ifdef MOD_MAGNETIZ
            DOXLR  = .TRUE.
            MSUSCP = .TRUE.
#endif
            GO TO 100
  105    CONTINUE
C&&&& NMR parameters:
            DOXLR  = .TRUE.
            SHIELD = .TRUE.
            SPNSPN = .TRUE.
            GO TO 100
  106    CONTINUE
C&&&& SPIN-S: NMR indirect spin-spin coupling constants:
            DOXLR  = .TRUE.
            SPNSPN = .TRUE.
            GO TO 100
  107    CONTINUE
C&&&& SHIELD: NMR shielding constants:
            DOXLR = .TRUE.
            SHIELD = .TRUE.
            GO TO 100
  108    CONTINUE
C&&&& PRINT: print level
            READ(LUCMD,*) IPRPRP
            GO TO 100
  109    CONTINUE
C&&&& ABUND: abundance of atomic type
            READ (LUCMD,*) ABUND
            GO TO 100
  110    CONTINUE
C&&&& EFG - Electric field gradients
            EFG      = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  111    CONTINUE
C&&&& NQCC - Nuclear Quadrupole Coupling Constants
            NQCC     = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  112    CONTINUE
C&&&& DSO - diamagnetic spin-orbit 
            DSO   = .TRUE.
            DOEXP = .TRUE.
            GO TO 100
  113    CONTINUE
C&&&& MGRAD - molecular gradient
            MGRAD = .TRUE.
            GO TO 100
  114    CONTINUE
C&&&& NSTDIA diamagnetic contribution to shielding
            NSTDIA   = .TRUE.
            DOEXP = .TRUE.
            GO TO 100
  115    CONTINUE
C&&&& ESR parameters
            DOESR = .TRUE.
            ESRGTENS = .TRUE.
            ESRANISO = .TRUE.
            GO TO 100
  116    CONTINUE
C&&&& PVC Parity violation - chirality
            PVC   = .TRUE.
            DOEXP = .TRUE.
            GO TO 100
  117    CONTINUE
C&&&& EFT - traceless electric field third derivatives
            EFT      = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  118    CONTINUE
C&&&& EFTNTL - non-traceless electric field third derivatives
            EFTNTL   = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
  119    CONTINUE
C&&&& PV contribution to NMR shielding
            DOXLR      = .TRUE.
            PVC_SHIELD = .TRUE.
            IF(NBSYM .GT. 1) CALL QUIT('.PVCSHI only for C1 symmetry')
            GOTO 100
  120    CONTINUE
C&&&& EFN: Electric field at nuclei
            DOEXP = .TRUE.
            EFN   = .TRUE.
            GOTO 100
  121    CONTINUE
C&&&& WAVE FUNCTIONS - define which wave functions we should use
            READ(LUCMD,*) NPRPWF
            IF (NPRPWF .LT. 1 .OR. NPRPWF .GT. MXPRPWF) THEN
               WRITE(LUPRI,'(1X,2A,I10)')
     &              '*** ERROR in PRPINP *** ',
     &              'Illegal value for NPRPWF: ',NPRPWF
               CALL QUIT('*** ERROR in PRPINP ***')
            END IF
            DO I = 1, NPRPWF
               READ(LUCMD,'(1X,A4)') PRPWF(I)
            END DO
            GOTO 100
 122     CONTINUE
C&&&& RHONUC - electronic density at nuclei
            RHONUC   = .TRUE.
            DOEXP    = .TRUE.
            GO TO 100
 123     CONTINUE
C&&&& FIRST HYPERPOLARIZABILITY - compute complete beta tensor
            DOXQR = .TRUE.
            CALL DEF_1STHYP(WORK,LWORK)
            GO TO 100
 124     CONTINUE
C&&&& TWO-PHOTON ABSORPTION - compute all two-photon absorption matrix elements
            CALL DEF_TPA(WORK,LWORK,DOTPA)
            GO TO 100
 125     CONTINUE
C&&&& .VERDET - Calculates the Verdet constant with the QR module
            DOVER = .TRUE.
            DOXQR = .TRUE.
            CALL DEF_VERDET(WORK,LWORK)
            GO TO 100
 126     CONTINUE
C&&&& ROTG rotational g-tensor
CMI ... at present it is not implemented
            ROTG  = .TRUE.
            DOXLR = .TRUE.
            GO TO 100
  127    CONTINUE 
C&&&& MSUSCD diamagnetic contribution to the magnetizabilities
#ifdef MOD_MAGNETIZ
            MSUSCDIA = .TRUE.
            DOEXP    = .TRUE.
#endif
            GO TO 100
  128    CONTINUE 
C&&&& PV contribution to spin-spin coupling
            DOXLR        = .TRUE.
            PVC_SPINSPIN = .TRUE.
            IF(NBSYM .GT. 1) CALL QUIT('.PVCSPN only for C1 symmetry')
            GO TO 100
  129    CONTINUE 
C&&&& NOPCT - even if DOPCT is .true. (because of default settings) 
C             do NOT activate picture change transformation...
            NOPCT = .TRUE.
            GO TO 100
  130    CONTINUE 
C&&&& .EXCDIP - Calculates excited state expectation values of dipole op.
            CALL EXCPRPINP(1,DOEXCPRP,EXCPRPTYP,WORK,LWORK)
            CALL EXCPRPINP(4,DOEXCPRP,EXCPRPTYP,WORK,LWORK)
            GO TO 100
  131    CONTINUE 
C&&&& .MAGOUT - specify the hard-output for the magnetizability
            READ(LUCMD,*) IMAGOUT
            GO TO 100
  132    CONTINUE 
C&&&& OPTROT - do optical rotation
            DOXLR = .TRUE.
            OPTROT = .TRUE.
            GO TO 100
  133    CONTINUE 
C&&&& .C6 - evaluate C6 dispersion coefficients
            DOXLR  = .TRUE.
            EPOLAR = .TRUE.
            DOVDW  = .TRUE.
            GO TO 100
  134    CONTINUE 
C&&&& RKBIMPort; read RKB coefficients and add complementary UKB space
            DORKBIMP = .TRUE.
            GO TO 100
  135    CONTINUE 
C&&&& BXLAO 
            bxlao  = .true.
            doxlr  = .true.
            london = .true.
            GO TO 100
  136    CONTINUE
C&&&& BYLAO 
            bylao  = .true.
            doxlr  = .true.
            london = .true.
            GO TO 100
  137    CONTINUE
C&&&& BZLAO 
            bzlao  = .true.
            doxlr  = .true.
            london = .true.
            GO TO 100
  138    CONTINUE
C&&&& EFFDEN; calculate effective electronic densities at the nuclei
            DOEXP  = .TRUE.
            EFFDEN = .TRUE. 
            GO TO 100
  139    CONTINUE
C&&&& XMLOUT: Turn on xml format output by opening a file

            GO TO 100
  140    CONTINUE
C&&&& XXXXXX
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in PRPINP.')
      END IF
  180 CONTINUE
C
C     Print section
C     =============
C
      CALL TITLER('Property module','*',129)
      WRITE(LUPRI,'(1X,A,I3)') '* Print level:',IPRPRP
      WRITE(LUPRI,'(2A)')   ' * Input label:',WRDSRC
      WRITE(LUPRI,'(1X,A)')
     &     '* Properties calculated for the following wave functions: '
      DO I = 1, NPRPWF
         WRITE(LUPRI,'(5X,I1,A,A)') I,': ',PRPWF(I)
      END DO

CMI ...for magnetic susceptibility, NMR shielding constants and rotational
CMI   g_tensors operator centers might be changed
      WRITE(LUPRI,'(1X,A)') 
     &'These initial settings of center and origins'//
     &' might be changed later:'
      WRITE (LUPRI,'(1X,A,3F18.10)')
     &     '* Operator center (a.u.):', (ORIGIN(I),I=1,3),
     &     '* Gauge origin    (a.u.):', (GAGORG(I),I=1,3),
     &     '* Dipole origin   (a.u.):', (DIPORG(I),I=1,3)
      GAGORGX = GAGORG(1)
      GAGORGY = GAGORG(2)
      GAGORGZ = GAGORG(3)

!     set compulsory default values for X2C module.
      if(x2c)then
        DOPCT = .true.
        NOPCT = .false.
      end if

      IF (DOPCT) THEN
       IF (.NOT.NOPCT) THEN
         WRITE(LUPRI,'(2X,A)') 
     & '...perform the picture change transformation of the'//
     & ' four-comp. property operators'
       ELSE
         WRITE(LUPRI,'(2X,A)') 
     & '...do NOT perform the picture change transformation of the'//
     & ' four-comp. property operators - take only LL block(s)'
       ENDIF
      ELSE
        WRITE (LUPRI,'(1X,A)')
     &  '* do NOT perform the picture change transformation of the'//
     &  ' four-component property operators'
      ENDIF
      IF(DORKBIMP) THEN
        IF(URKBAL) THEN
          WRITE(LUPRI,'(1X,A,A)')
     &    '* Read RKB coefficients ',
     &    'and extended with complementary UKB space'
        ELSE
          DORKBIMP = .FALSE.
          WRITE(LUPRI,'(1X,A)') 
     &      '*** WARNING! RKBIMPort ignored since no URKBAL'
        ENDIF
      ENDIF
C
C     Process input for various program sections
C     ==========================================
C
C     Now read input
C
  200 CONTINUE
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in PRPINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in PRPINP, error in prompt.')
      END IF
    2 CONTINUE
        CALL EXPINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    3 CONTINUE
        CALL XLRINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    4 CONTINUE
        CALL XPPINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    5 CONTINUE
        CALL GRDINP(WORD,.FALSE.)
        GO TO 200

!     *NMR
    6 CONTINUE
CMI ... changed for LAO-based magnetic properties
        CALL NMRINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200

    7 CONTINUE
#ifdef MOD_ESR
        CALL ESRINP(WORD,.FALSE.,WORK,LWORK)
#else
        CALL QUIT('*ESR    not included in this version')
#endif
        GO TO 200
    8 CONTINUE
         CALL XQRINP(WORD,.FALSE.,WORK,LWORK)
         GO TO 200
C     *STEX
    9 CONTINUE
         CALL STEX_INPUT(WORD,.FALSE.,WORK,LWORK)
         GO TO 200

!     *OPENRSP    
   10 continue
#ifdef MOD_OPENRSP
         dirac_cfg_openrsp = .true.
         call init_dirac_interface()
         call move_to_next_star(word)
         if (openrsp_cfg_skip_1el) then
            write(*, *) 'OpenRSP: skipping calculation of 1el integrals'
         else
            call order_1el_integrals()
         end if
#else
         call quit('*OpenRSP not included in this version')
#endif
      go to 200

    1 CONTINUE
      IF(DOXLR) CALL XLRINP(WORD,.TRUE.,WORK,LWORK)
CMI ... read again the input
      CALL NMRINP(WORD,.TRUE.,WORK,LWORK)
!gosia: i have to reset these variables once again...
!or always remember to put *NMR .LONDON in the input...
      if (bxlao .or. bylao .or. bzlao) then
        london = .true.
        doxlr  = .true.
      end if
!
      CALL GRDINP(WORD,.TRUE.)
      IF(DOXPP) CALL XPPINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOXQR) CALL XQRINP(WORD,.TRUE.,WORK,LWORK)
#ifdef MOD_ESR
      IF(DOESR) CALL ESRINP(WORD,.TRUE.,WORK,LWORK)
#endif
C===========================================================================
C Define special operators for the list of properties programmed in DIRAC
C===========================================================================
      CALL PRPDEF(WORK,LWORK)

C .... check the change of the gauge origin ....
      THRS1=1.0D-7
      CHGO= DABS(GAGORGX-GAGORG(1)).GE.THRS1.OR.
     &      DABS(GAGORGY-GAGORG(2)).GE.THRS1.OR.
     &      DABS(GAGORGZ-GAGORG(3)).GE.THRS1
      IF (CHGO) THEN
       WRITE (LUPRI,'(1X,A,3F18.10)')
     &     '* Operator center:', (ORIGIN(I),I=1,3),
     &     '* Gauge origin   :', (GAGORG(I),I=1,3),
     &     '* Dipole origin  :', (DIPORG(I),I=1,3)
      ENDIF

      CALL FLSHFO(LUPRI)
      CALL QEXIT('PRPINP')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck expinp */
      SUBROUTINE EXPINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for expectation values
C
C     Written by T.Saue - May 1996
C     Last revision: May 26  - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 8, DM1 = -1.0D00, D0 = 0.0D00)
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbxpr.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.OPERAT','.ORBANA','.PROJEC',
     &            '.PRPCAN','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
      NEWDEF = (WORD .EQ. '*EXPECT')
      IF (SET) THEN
         IF (NEWDEF)
     &      CALL QUIT('Only one "*EXPECT" input section allowed')
C        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C     Initialize /DCBEXP/
C     ===================
C
      IPREXP = 0
      NEXPP  = 0
      ORBANA = .FALSE.
      PRJANA = .FALSE.
      PRPCAN = .FALSE.
C
C     Process input for DCBEXP
C     =========================
C
      INPERR = 0
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
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in EXPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in EXPINP.')
    1          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPREXP
               GO TO 100
    2          CONTINUE
C&&&& OPERATOR: Define operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPREXP)
                  CALL OP1IND('EXPINP',INDEXP,LEXPP,NEXPP,INDXPR,
     &                        MAXEXP)
                  DOEXP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
C&&&& ORBANA: Analyze expectation value in terms of individual orbitals
               ORBANA = .TRUE.
               GO TO 100
    4          CONTINUE
C&&&& PROJEC: Analyze expectation value in terms of fragment orbitals
               PRJANA = .TRUE.
               GO TO 100
    5          CONTINUE
C&&&& PRPCAN: Generate canonical occupied property orbitals
               PRPCAN = .TRUE.
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in EXPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in EXPINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' EXPINP: Expectation values'
      IF(ORBANA.OR.PRJANA) THEN
        WRITE(LUPRI,'(A)') 
     &   '* Expectation values analyzed in terms of '
        IF(ORBANA) WRITE(LUPRI,'(3X,A)') '- individual orbitals'
        IF(PRJANA) WRITE(LUPRI,'(3X,A)') '- fragment orbitals'
      ENDIF
      IF(PRPCAN) THEN
        WRITE(LUPRI,'(A)') 
     &     '* Generating canonical orbitals for each property.'
      ENDIF
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(ICHANG.GT.0) THEN
        WRITE(LUPRI,'(A)') 
     &  ' * The following expectation values will be calculated:'
      ENDIF
      IF(NEXPP.GT.0) THEN
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF  
      IF(DIPOLE) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Dipole moment'
      ENDIF
      IF(QUADRU) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Quadrupole moment'
      ENDIF
      IF(PVC) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Parity violation - chirality'
      ENDIF
      IF(EFN) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Electric field at nuclei'
      ENDIF
      IF(EFG) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Electric field gradient at nuclei'
      ENDIF
      IF(EFT) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')
     &     '- traceless electric field third derivattives'
      ENDIF
      IF(EFTNTL) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')
     &     '- non-traceless electric field third derivattives'
      ENDIF
      IF(NQCC) THEN    
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(4X,A)')'- Nuclear Quadrupole Coupling Constant'
      ENDIF
C
      CALL PRSYMB(LUPRI,'-',75,0)
  999 CONTINUE
      IF (INPERR.GT.0) CALL QUIT('Input error in *EXPECT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xprinp */
      SUBROUTINE XPRINP(LUCMDX,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
C  Last revision: May 14 1996 - tsaue
C                 Febr. 2011   M Ilias (TRS_FLIP)
C
C  Assign a new property operator based on the input
C  and the /XCBPRP/ list.
C
C  An operator in 4-component theory is a linear combination
C  of products of a 4x4 matrix and some scalar operator P.
C  The 4x4 matrices are the alpha matrices and the identity.
C  XPRINP recognizes eighteen different kinds of operators:
C
C    1. DIAGONAL
C       P             * scalar operator
C
C    2. XALPHA
C       [alpha_x]P    * x-component of alpha times scalar operator
C
C    3. YALPHA
C       [alpha_y]P    * y-component of alpha times scalar operator
C
C    4. ZALPHA
C       [alpha_z]P    * z-component of alpha times scalar operator
C
C    5. XAVECTOR
C       [Alpha x P]_x * vector product of alpha and vector operator,
C                       x-component    
C    6. YAVECTOR
C       [Alpha x P]_y * vector product of alpha and vector operator,
C                       y-component    
C    7. ZAVECTOR
C       [Alpha x P]_z * vector product of alpha and vector operator,
C                       z-component    
C    8. ALPHADOT
C       A.P           * dot-product of alpha and vector operator
C
C    9. GAMMA5        
C       gamma5 P      * The gamma5 matrix times a scalar operator
C                       gamma5 = | 0 I | 
C                                | I 0 |
C   10. XSIGMA
C       [sigma_x]p    * x-component of sigma times scalar operator
C
C   11. YSIGMA
C       [sigma_y]p    * y-component of sigma times scalar operator
C
C   12. ZSIGMA 
C       [sigma_z]p    * z-component of sigma times scalar operator
C
C   13. XBETASIG
C       [betasig_x]p  * x-component of beta sigma times scalar operator
C
C   14. YBETASIG
C       [betasig_y]p  * Y-component of beta sigma times scalar operator
C
C   15. ZBETASIG
C       [betasig_z]p  * z-component of beta sigma times scalar operator
C
C   16. XiBETAAL
C       i[betaalp_x]p * x-component of beta alpha times scalar operator
C 
C   17. YiBETAAL
C       i[betaalp_y]p * y-component of beta alpha times scalar operator
C
C   18. ZiBETAAL
C       i[betaalp_z]p * z-component of beta alpha times scalar operator
C
C   19. BETA   
C       [beta]P * beta times scalar operator
C
C   20 .is occupied - iS.P         * dot-product of sigma and vector operator
C
C   21. iBETAGAM
C       i[beta*gamma5] P           * i beta*gamma5 times scalar operator
C   
C
C  Operator types 5,6 and 7 form the components of a vectorproduct 
C  of the alpha matrices and some scalar operator.
C
C  The scalar operators may be supplied with factors, e.g. field strengths.
C
C  The pointer IKW controls what is to be read
C     0   Property name
C     1   Keywords
C     2   Property labels (to be found on file)
C     3   Class label
C     4   Separate factor for each component
C     5   Common factor for all components
C     6   Keyword to force the time-reversal symmetry flip
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
      CHARACTER*7 WORD
C
      PARAMETER (NTABKW = 32)
      CHARACTER*8 TABKW(NTABKW)
      CHARACTER   LINE*80, KEYWD*8, TOKEN16*16, BLANK8*8, BLANK16*16
      DIMENSION PFAC(3),WORK(LWORK)
      CHARACTER   PLABEL(3)*8, PNAME*16
      LOGICAL     CLABEL 
      LOGICAL     TRS_FLIP
C
      PARAMETER  (D1 = 1.0D0)
      PARAMETER (BLANK8='        ', BLANK16='                ' )
      SAVE TABKW
      DATA TABKW/'CLASS   ', 'FACTORS ', 'COMFACTO', 'CMULT   ',
     &           'TRS_FLIP',
     &           'FOCKMAT ', 'DIAGONAL', 'XALPHA  ', 'YALPHA  ',
     &           'ZALPHA  ', 'XAVECTOR', 'YAVECTOR', 'ZAVECTOR',
     &           'ALPHADOT', 'GAMMA5  ', 'XSIGMA  ', 'YSIGMA  ',
     &           'ZSIGMA  ', 'XBETASIG', 'YBETASIG', 'ZBETASIG',
     &           'XiBETAAL', 'YiBETAAL', 'ZiBETAAL', 'BETA    ',
     &           'XXXXXXXX', 'iBETAGAM', 'XXXXXXXX', 'XXXXXXXX',
     &           'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX'            /
C
      CALL QENTER('XPRINP')
C
C     Initialization
C     ==============
C
      INDXPR = -1
      ISYXPR = -1
      ITRXPR = -1
      PFAC(1) = D1
      PFAC(2) = D1
      PFAC(3) = D1
      CLABEL  = .FALSE.
      TRS_FLIP = .FALSE.
      !mi - position of the Fock operator in the TABKW list, or "border" index
      I_FOCKMAT=6
C
      IKW = 0
C
C     *****************
C     *** Read line ***
C     *****************
C
  100 READ (LUCMDX,'(A)',END=8000) LINE
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') '***** Input line read in XPRINP:'
         WRITE (LUPRI,'(A)') LINE
      END IF
C
C     Line to be skipped
C     ==================
      IF (LINE(1:1) .EQ. '!' .OR. LINE(1:1) .EQ. '#') GO TO 100
C
C     Line is keyword from other section: analyze what you already have
C     =================================================================
C
      IF (LINE(1:1) .EQ. '*' .OR. LINE(1:1) .EQ. '.') THEN
         WORD = LINE(1:7)
         BACKSPACE LUCMDX
         GO TO 500
      END IF
C
C     Line must start with a blank
C     ============================
C
      IF (LINE(1:1) .NE. ' ') THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A/A)')
     *      ' XPRINP: INPUT ERROR following '//WORD,
     *      ' First character is not blank; line in error:', LINE
      END IF
C
C     Line accepted; proceed
C     ======================
C
      ITKN = 2
C
C     XPRTKN returns initial adress(ITKN) and length (LTKN) of next token
C     on line. If ITKN -1, then no more tokens on LINE.
C
  200 CALL XPRTKN(LINE,ITKN,LTKN)
C
C     If no more tokens, read new line
C     ================================
C
      IF (LTKN .LE. 0) GO TO 100
C
C     ************************
C     *** Specify operator ***
C     ************************
C
C
C     Read property name
C     ==================
C
      IF (IKW .EQ. 0) THEN
C        IKW .eq. 0 :  Property name
C
         PNAME     = BLANK16
         LTKNR     = MIN(16,LTKN)
         READ(LINE(ITKN:80),'(A)') PNAME(1:LTKNR)
C        default diagonal operator name = label
         PLABEL(1) = PNAME(1:8)
         IPTYP     = 1
         NPCOMP    = 1
         NFAC      = 0
         IKW       = 1
      ELSE IF (IKW .EQ. 1) THEN
C        IKW .eq. 1 :  Property operator keywords
         LKW = MIN(8,LTKN)
         READ(LINE(ITKN:80),'(A)') KEYWD(1:LKW)
         IF (LKW .LT. 3) THEN
            INPERR = INPERR + 1
             WRITE (LUPRI,'(/A/A)')
     &       ' XPRINP: INPUT ERROR following '//WORD,
     &       ' Keyword "'//KEYWD(1:LKW)//
     &       '" too short, min 3 char required.'
         END IF
C
C        Look for keyword from table
C        ===========================
C
         DO I = 1,NTABKW
         IF (KEYWD(1:LKW) .EQ. TABKW(I)(1:LKW)) THEN
C
C          If keyword defines operator type, then initialize
C          =================================================
C
           IF ( I .GT. I_FOCKMAT ) THEN
C             (the first five entries in TABKW are not operator types)
CMI      ---- assignement of IPTYP needed in further routines -----
              IPTYP = I - I_FOCKMAT
              NLABEL= 0
              NFAC  = 0
              IKW = 2
           END IF

C
C          If keyword Fock defines matrix read, then initialize
C          ====================================================
C          (the matrix is assumed totally symmetric 
C           under point group and time reversal symmetry)
C
           IF (I .EQ. I_FOCKMAT) THEN
              IPTYP = 0
              NLABEL= 0
              NFAC  = 0
              IKW = 2
              ISYXPR = 1
              ITRXPR = 1
           END IF
 
           ! ... initialize ...
           TRS_FLIP=.FALSE.
           ! write(LUPRI,*) 'operatorik:',I,KEYWD(1:LKW)

!mi     ...jump on operator types according to the variable I prepared above, I : 1-26
!           4-8 are operators characteritics, TABKW(1-5)
           GO TO (4,5,6,7,8,1,1,1,1,1,2,2,2,3,
     &            1,1,1,1,1,1,1,1,1,1,1,1,1), I
         END IF
         END DO
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')
     &     ' XPRINP: INPUT ERROR following '//WORD,
     &     ' Keyword "'//KEYWD(1:LKW)//'" not recognized.'
         GO TO 300
C&&&&    One component: DIAGONAL, [XYZ]ALPHA, GAMMA5, [XYZ]SIGMA,
C&&&&                   [XYZ]BETASIG, [XYZ]iBETAAL, BETA, FOCKMAT, iBETAGAMMA5
    1    CONTINUE
           NPCOMP = 1
         GO TO 300
C&&&&    Two components: [XYZ]AVECTOR
    2    CONTINUE
           NPCOMP = 2
C          note that the minus in the second component
C          is introduced later
         GO TO 300
C&&&&    Three components: ALPHADOT
    3    CONTINUE
           NPCOMP = 3
         GO TO 300
C&&&&    Class keyword, CLASS, TABKW(1)
    4    CONTINUE
           CLABEL = .TRUE.
           IKW = 3
         GO TO 300
C&&&&    Read a multiplicative factor for each component, FACTORS, TABKW(2)
    5    CONTINUE
           IKW  = 4
         GO TO 300
C&&&&    Read a common factor for each component, COMFACTOR, TABKW(3)
    6    CONTINUE
           IKW  = 5
         GO TO 300
C&&&     Multiply with current value of the speed of light, CMULT
    7    CONTINUE
           PFAC(1) = PFAC(1)*CVAL          
           PFAC(2) = PFAC(2)*CVAL          
           PFAC(3) = PFAC(3)*CVAL 
         GO TO 300
C&&&     miro: Time-reversal symmetry flip, TRS_FLIP, as newly added operator characteristic
    8    CONTINUE
           TRS_FLIP=.TRUE.
           IKW=6
         GO TO 300
         !GO TO 100
C
C     Read property labels for each component
C     =======================================
C
      ELSE IF (IKW .EQ. 2) THEN
C        IKW .eq. 2 :  Property labels (to be found on file)
         NLABEL = NLABEL + 1
         PLABEL(NLABEL) = BLANK8
         LTKNR = MIN(8,LTKN)
         READ(LINE(ITKN:80),'(A)') PLABEL(NLABEL)(1:LTKNR)
         IF (NLABEL .GE. NPCOMP) IKW = 1
C
C     Read class keyword
C     ==================
C
      ELSE IF (IKW .EQ. 3) THEN
C        IKW .eq. 3 :  Class label
CHJMAERKE: implement interpretation of class keyword
!miro - missing field, HJJ has an idea...
C
C     FACTORS, read factor for each component
C     ==================================================
C
      ELSE IF (IKW .EQ. 4) THEN
C        IKW .eq. 4 :  Separate factor for each component
         NFAC = NFAC + 1
         JTKN = ITKN - 1 + LTKN
         I    = 16 - LTKN + 1
         TOKEN16 = BLANK16
         TOKEN16(I:16) = LINE(ITKN:JTKN)
         READ(TOKEN16,'(F16.0)') TMP
         PFAC(NFAC) = TMP*PFAC(NFAC)
         IF (NFAC .GE. NPCOMP) IKW = 1
C
C     COMFACTOR, read a common factor for all components
C     =======================================
      ELSE IF (IKW .EQ. 5) THEN
C        IKW .eq. 5 :  Common factor for all components
         JTKN = ITKN - 1 + LTKN
         I    = 16 - LTKN + 1
         TOKEN16 = BLANK16
         TOKEN16(I:16) = LINE(ITKN:JTKN)
         READ(TOKEN16,'(F16.0)') COMFAC
         PFAC(1) = COMFAC*PFAC(1)
         PFAC(2) = COMFAC*PFAC(2)
         PFAC(3) = COMFAC*PFAC(3)
         IKW = 1

!         TRS_FLIP
      ELSE IF (IKW .EQ. 6) THEN
         write(LUPRI,*) 'TRS_FLIP !'
         TRS_FLIP=.TRUE.
         IKW=1
C
C     Unknown IKW
C     ===========
C
      ELSE
C        IKW .gt. 5
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &    'Programming error in XPRINP, unknown IKW=',IKW
      END IF
C
  300 ITKN = ITKN + LTKN
      GO TO 200
C
C**********************************************************
C***** Look up in list of operators already defined. ******
C**********************************************************
C
C     Find INDXPR, ISYXPR, ITRXPR
C
  500 CONTINUE
      IF (CLABEL) THEN
CHJMAERKE: implement interpretation of class keyword
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')
     &      ' XPRINP: INPUT ERROR following '//WORD,
     &      ' Keyword "'//KEYWD(1:LKW)//'" not implemented yet.'
      ELSE
      ! miro-output and XPRIND
         IF (IPRINT.GE.5) THEN
           WRITE(LUPRI,'(/,6x,a)') '*** Output from XPRINP ***'
           WRITE(LUPRI,'(2x,a,a16)') 'PNAME=',PNAME
           WRITE(LUPRI,'(2x,a,i2,a,i2)') 
     &     'IPTYP=',IPTYP,' NPCOMP=',NPCOMP
           WRITE(LUPRI,'(2x,a,3a8)') 
     &     'PLABEL(1-NPCOMP)=',(PLABEL(I),I=1,NPCOMP)
           WRITE(LUPRI,'(2x,a,3f10.6)') 
     &     'PFAC(1-NPCOMP)=',(PFAC(I),I=1,NPCOMP)
           WRITE(LUPRI,'(2x,a,i2,a,i2,a,i2)') 
     &     'initial values: INDXPR=',INDXPR,
     &    ' ISYXPR=',ISYXPR,' ITRXPR=',ITRXPR
         ENDIF
CMI      ... assign remaining parameters
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         IF (TRS_FLIP) THEN
           ! miro: perform the time-reversal symmetry flip of the operaror
           CALL OPERATOR_TRS_FLIP(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                            INDXPR,ISYXPR,ITRXPR,IPRINT)
           TRS_FLIP=.FALSE.
         ENDIF
         IF (IPRINT.GE.8) THEN
           WRITE(LUPRI,'(/,6x,a)')
     &      '*** Output from XPRINP after XPRIND ***'
           WRITE(LUPRI,'(2x,a,i2,a,i2,a,i2)') 
     &     'INDXPR=',INDXPR,' ISYXPR=',ISYXPR,' ITRXPR=',ITRXPR
           WRITE(LUPRI,'(2x,a,i2,a,i2)') 
     &     'IPTYP=',IPTYP,' NPCOMP=',NPCOMP
           WRITE(LUPRI,*)
         ENDIF
      END IF
      GO TO 9999
C
 8000 INPERR = INPERR + 1
      WORD = ' INPEOF'
      WRITE (LUPRI,'(/A/A)')
     &' XPRINP: INPUT ERROR following '//WORD,
     &' End of file encountered; input must finish with "**something".'
C
 9999 CALL QEXIT ('XPRINP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xprtkn */
      SUBROUTINE XPRTKN(LINE,ITKN,LTKN)
C
C  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
C
C  Find initial address (ITKN) and length (ITKN) of next token
C  If no more tokens return ITKN = -1
C
      PARAMETER (MAXLL = 80)
      CHARACTER*(MAXLL) LINE
      CHARACTER DELIM*2, DELIM2*1
      PARAMETER (DELIM = '''"')
C
      CALL QENTER('XPRTKN')
C
C     Find first character in next token
C
      DO I = ITKN,MAXLL
         IF (LINE(I:I) .NE. ' ') THEN
            IF (INDEX(DELIM,LINE(I:I)) .GT. 0) THEN
               DELIM2 = LINE(I:I)
               IND1 = I + 1
            ELSE
               DELIM2 = ' '
               IND1 = I
            END IF
            GO TO 100
         END IF
      END DO
      ITKN = -1
      LTKN = 0
      GO TO 9999
C
C     Find last character of this token
C
  100 CONTINUE
      INDL = MAXLL
      DO I = IND1,MAXLL
         IF (LINE(I:I) .EQ. DELIM2) THEN
            INDL = I - 1
            LINE(I:I) = ' '
            GO TO 200
         END IF
      END DO
  200 ITKN = IND1
      LTKN = (INDL + 1 - IND1)
C
 9999 CALL QEXIT('XPRTKN')
      RETURN
      END
      SUBROUTINE OPERATOR_TRS_FLIP(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
! Perform the operator time-reversal symmetry flip (thanks to TRS_FLIP characteritics)
!
! Called by XPRINT
!
! M.Ilias, Banska Bystrica, February 2011
!
#include "implicit.h"
#include "priunit.h"
#include "dcbxpr.h"
#include "dcbprl.h"
      CHARACTER PNAME*16
      DIMENSION PFAC(NPCOMP)
      WRITE(LUPRI,*)
     &'Flipping time reversal symmetry for the selected operator ',
     & PNAME
      IPRPTIM(INDXPR) =  - IPRPTIM(INDXPR)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xprind */
      SUBROUTINE XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
C*********************************************************************** 
C
C  Given a 4-component operator this routine returns an index
C  INDXPR to the list of operators in COMMON block XCBPRP and
C  returns its symmetry under time reversal(ITRXPR) and its
C  spatial symmetry(ISYXPR).
C
C  This routine is called from many places 
C
C  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
C  Revised May 14 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PNAME*16, PLABEL(NPCOMP)*8
      DIMENSION PFAC(NPCOMP),WORK(LWORK)
C
C Used from common blocks:
C  XCBPRP: everything
C  XCBPRL: PRPLBL()
C  XCBPRI: LUPRI,NINFO
C
#include "dcbxpr.h"
#include "dcbprl.h"
C
!     CALL QENTER('XPRIND')
      INDXPR = -1

      IF (IPRINT.GE.15) THEN
        WRITE(LUPRI,'(/,4x,a)') '***** Output from XPRIND *****'
        WRITE(LUPRI,'(2x,a,a16)') 'Entering PNAME=',PNAME
        WRITE(LUPRI,'(2x,a,i2)') 'Entering NPCOMP=',NPCOMP
        WRITE(LUPRI,'(2x,a,6a8)')
     &   'Entering PLABEL=',(PLABEL(I),I=1,NPCOMP)
        WRITE(LUPRI,'(2x,a,6f12.6)') 
     &    'Entering PFAC=',(PFAC(I),I=1,NPCOMP)
        WRITE(LUPRI,'(2x,a,i2,a,i3)') 
     &   'Entering IPTYP=',IPTYP,' INDXPR=',INDXPR
      ENDIF

C================================================================
C       Search for this property operator in existing list
C================================================================
      DO 200 I = 1,NPRPS
C        Check operator type
         IF (IPRPTYP(I) .NE. IPTYP) GO TO 200
C        Check factors and property labels for each component
         DO 100 J = 1,NPCOMP
            IF (PFAC(J) .NE. FACPRP(J,I)) GO TO 200
            K = IPRPLBL(J,I)
            IF (PRPLBL(K) .NE. PLABEL(J)) GO TO 200
  100    CONTINUE
C
         IF (PNAME .NE. PRPNAM(I)) THEN
C        Operator exists, but under another name; rename
            WRITE (LUPRI,'(3(/A))')
     &      ' XPRIND INFO: property operator named '//PNAME,
     &      '   is already existing under the name '//PRPNAM(I),
     &      '   previous operator renamed to       '//PNAME
            PRPNAM(I) = PNAME
         END IF
         INDXPR = I
         GO TO 1000
  200 CONTINUE

C========================================================================
C     This property operator has not been defined previously.
C     Append the new property operator to /XCBPRP/ list.
C=========================================================================
      CALL XPRADD(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,INDXPR,
     &            WORK,LWORK,IPRINT)
C
 1000 CONTINUE
      IF(IPTYP.EQ.0) THEN
C       ... i.e. FOCKMAT 
        IPRPSYM(INDXPR) = ISYXPR
        IPRPTIM(INDXPR) = ITRXPR
      ELSE
        ISYXPR = IPRPSYM(INDXPR)
        ITRXPR = IPRPTIM(INDXPR)
      ENDIF
!     CALL QEXIT('XPRIND')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xpradd */
      SUBROUTINE XPRADD(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,INDXPR,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C  Add property to list of properties in COMMON block XCBPRP
C  Copyright (c) 6-Feb-1995 Hans Joergen Aa. Jensen
C
C  Entering: PNAME
C            IPTYP
C            NPCOMP
C
C  Routine is called only from XPRIND.
C
C  MI (&HAaaJ), 2003,2004: Extended for LAO integral types
C               as well as other integrals
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PNAME*16, PLABEL(NPCOMP)*8
      DIMENSION PFAC(NPCOMP),WORK(LWORK)
C
C Used from common blocks:
C  XCBPRP: everything
C  XCBPRI: LUPRI,NINFO
C
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbham.h"
CMI ... needs CVAL
#include "dcbgen.h"

      CHARACTER*4 PCOMB
C     ... the four characters in PCOMB define inclusion of the
C         LL, SL, LS, and SS blocks of one-electron operators.
C
      CALL QENTER('XPRADD')
C
C     This property operator has not been defined previously.
C     Append the new property operator to /XCBPRP/ list.


      IF (IPRINT .GE. 6) THEN
         WRITE (LUPRI,'(/,4X,A/2A/A,I5/A,I5/A,1P,3(D20.10))')
     &   '*** XPRADD called for specifying PCOMB ***',
     &   '  PNAME   =',PNAME,
     &   '  IPTYP   =',IPTYP,
     &   '  NPCOMP  =',NPCOMP,
     &   '  PFAC    =',(PFAC(I),I=1,NPCOMP)
         WRITE (LUPRI,'(A,3(A8,2X))')
     &   '  PLABEL  =',(PLABEL(I),I=1,NPCOMP)
      END IF
CMI  List of exceptions:
CMI/Dec.2005 I guess this routine could be reprogrammed as exceptions are before and after cycle.
      IF ((PLABEL(1)(1:7).EQ.'BETAMAT').OR.(PNAME(6:10).EQ.'RM1H1')
     &    .OR.(PNAME(7:11).EQ.'RM2H1')) THEN
CMI (&TEC): RM1H1/RM2H1 is first/second London derivative of BETAMAT
         PCOMB = '000+'
CMI/HJAaJ
C ... NEWS for LAO-ccRDSUSLL in Levy-Leblond - cancel the SL(2,1) block
C      to get compatibility with DALTON
C     ELSE IF (PNAME(7:13).EQ.'RDSUSLL'.AND.LEVYLE) THEN
      ELSE IF (PNAME(7:13).EQ.'RDSUSLL'.AND.RDSUSLLMOD) THEN
C ... fix the LS block for the LevyLe option
         PCOMB = '00+0'
        IF (IPRINT.GE.6) WRITE(LUPRI,'(/A)')
     & 'XPRADD: Only the LS block kept for the RDSUSLL operator!'
      ELSE IF ((PNAME(1:5).EQ.'d|S>/'.OR.PNAME(1:3).EQ.'dS/').AND.BSS)
     &     THEN
        PCOMB='+000'
        IF (IPRINT.GE.6) WRITE(LUPRI,'(/A)')
     &  'XPRADD: Only the LL block for connection'//
     &  ' matrixes in two-component mode!'
      ELSE
         IF (IPTYP.EQ.0) THEN
C          ... for FOCKMAT
           PCOMB = 'xxxx'
           GOTO 200           
         ELSEIF (IPTYP.LT.0 .OR. IPTYP.GT.21) THEN
           CALL QUIT('XPRADD: Illegal IPTYP (not 1-20)!')
         ENDIF
         GO TO (1,2,2,2,2,2,2,2,2,
     &          1,1,1,3,3,3,4,4,4,3,4,4),IPTYP
C
C
C        DIAGONAL and SIGMA:
C
    1    PCOMB = '+00+'
CJTH DEBUG!!
#ifdef UNDEF
         write(6,*) 'debug debug: operator ',PNAME
         IF (PNAME(1:2) .EQ. 'EF') THEN
            PCOMB = '+000'
            write(6,*) '           : setting comp to',pcomb
         END IF
#endif
         GO TO 100
C
C        ALPHA and GAMMA5:
C
    2    PCOMB = '0++0'
Chj      PCOMB = 'FTTF'
         GO TO 100
C
C        "BETA" and "BETA SIGMA":
C
    3    PCOMB = '+00-'
         IF (LEVYLE) PCOMB = '+000'
         GO TO 100
C
C        "i BETA ALPHA" and "i BETA GAMMA5" :
C
    4    PCOMB = '0+-0'
         GO TO 100
C
  100    CONTINUE

CMI  ... for the Levy-Leblond mode cancel the SS block!
CMI ... check first if PLABEL(1) is not operator from
CMI ... pilot Dirac bare nucleus (from subroutine OP1INI)
        IF (LEVYLE.AND.PLABEL(1).NE.'OVERLAP '
     &       .AND.PLABEL(1).NE.'MOLFIELD'.AND.
     &       PLABEL(1).NE.'BETAMAT '
     &       .and.PLABEL(1).NE.'FDEVEMB ' ) PCOMB(4:4) = '0'
CLV ... same thing for ZORA!
        IF ((ZORA.AND..NOT.ZORA4).AND.PLABEL(1).NE.'OVERLAP '
     &       .AND.PLABEL(1).NE.'MOLFIELD'.AND.
     &       PLABEL(1).NE.'BETAMAT '
     &       .and.PLABEL(1).NE.'FDEVEMB ' ) PCOMB(4:4) = '0'
      END IF

CMI ... inform about the assigned PCOMB
      IF (IPRINT .GE. 6) THEN
         WRITE(LUPRI,'(A,A4)') '  assigned PCOMB:',PCOMB
      END IF

!     ...  assign proper factor for EDM  ...
      IF (PLABEL(1).EQ.'EDM') THEN
        PFAC(1) = -4.0D0 * CVAL * PFAC(1)
        WRITE(LUPRI,'(2X,A,D12.6)') 
     &  'EDM operator: operator strenght multiplied by -4*CVAL,'
     & //' PFAC(1)=', PFAC(1)
      ENDIF

C
C     Define property
C
 200  CONTINUE
      NPRPS  = NPRPS + 1
      IF(NPRPS.GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &   ' ** ERROR in XPRADD **  Redimension MAXPRPS to ',NPRPS
        CALL QUIT('XPRADD: Too may properties !')
      ENDIF
CMI   ... assign key variables ...
      INDXPR = NPRPS
C     Property name
      PRPNAM(NPRPS)  = PNAME
C     Property type
      IPRPTYP(NPRPS) = IPTYP
C     Pointers to property labels
      IF(IPTYP.EQ.0) THEN
        DO 300 I = 1,NPRPLBL
C         Check property label
          IF (PLABEL(1).NE.PRPLBL(I)) GOTO 300
          IF (PCOMB.NE.PDOINT(I))  GOTO 300
          INDX = I
        GO TO 400
  300   CONTINUE
        NPRPLBL = NPRPLBL + 1
        INDX    = NPRPLBL
  400   CONTINUE      
        FACPRP(1,NPRPS)  = PFAC(1)
        IPRPLBL(1,NPRPS) = INDX
        PRPLBL(INDX) = PLABEL(1)
        PDOINT(INDX) = PCOMB
      ELSE
CMI     ... (IPTYP.NE.0) operator has more components
        DO 600 J = 1,NPCOMP
          FACPRP(J,NPRPS)  = PFAC(J)
          INDX = INDXPL(PLABEL(J),PCOMB,WORK,LWORK,IPRINT)
          IF (INDX .LE. 0) THEN
            INDXPR = -1
            WRITE (LUPRI,'(3A)')
     &      'PLABEL, PCOMB not found in /XCBPRL/ : ', PLABEL(J), PCOMB
          END IF
          IPRPLBL(J,NPRPS) = INDX
  600   CONTINUE
        IF (INDXPR .EQ. -1) CALL QUIT('XPRADD: PLABEL, PCOMB not found')
C       Symmetry under time reversal and point group symmetry
        CALL PRPSYM(INDXPR,IPRINT)
      ENDIF
C
      CALL QEXIT('XPRADD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck indxpl */
      INTEGER FUNCTION INDXPL(PLABEL,PCOMB,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     9-Mar-1995 hjaaj
C     Complete revision : May 14 1996 - tsaue
C
C     Return value is pointer to list of property labels in 
C     COMMON block PRPLBL. The routine will furthermore update
C     list of integrals to calculate found in COMMON block
C     PRPCLB
C
C     If INDXPL .le. 0 then no match was found.
C
C     Last revision: May 15 1996 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PLABEL*8, PCOMB*4,CLASLB*7
C
C Used from common blocks:
C   XCBPRL : PRPFILE
C
#include "dcbprl.h"
#include "dcbcls.h"
C
      CALL QENTER('INDXPL')
C
      INDXPL = -1
      DO 100 I = 1,NPRPLBL
C       Check property label
        IF (PLABEL.NE.PRPLBL(I)) GOTO 100
        IF (PCOMB.NE.PDOINT(I))  GOTO 100
C
C       Label and block specification already defined; return
C
        INDXPL = I
        GO TO 9999
  100 CONTINUE
C
C     Class label not yet defined; first check for class label
C
      IF(PLABEL.EQ.'MOLFIELD') THEN
        CLASLB = 'MOLFLD '
        IORDER  = 0
      ELSEIF(PLABEL(1:7).EQ.'FDEVEMB') THEN
        CLASLB = 'FDEVEMB'
        IORDER   = 0
      ELSEIF(PLABEL(1:7).EQ.'BETAMAT') THEN
        CLASLB = 'BETAMAT'
        IORDER   = 0
      ELSEIF(PLABEL(1:3).EQ.'EDM') THEN
        CLASLB = 'EDM    '
        IORDER   = 0
      ELSEIF(PLABEL(1:2).EQ.'ED') THEN
        CLASLB = 'EFFDEN '
        IORDER   = 0
      ELSE
        CALL GETCLS(PLABEL,CLASLB,IORDER,IPRINT)
        IF(CLASLB.EQ.'Unknown') THEN
          WRITE(LUPRI,'(A,A8)')
     &    '***** INDXPL ERROR ***** No class label found for ',
     &    PLABEL
          CALL QUIT('INDXPL ERROR: No class label found !')
        ENDIF
      ENDIF
      NPRPCLS = NPRPCLS + 1
      IF (NPRPCLS .GT. MAXCLS) THEN
         WRITE(LUPRI,'(/A/A,I5,A/A,1X,A)')
     &      'INDXPL ERROR: MAXCLS in dcbcls.h too small',
     &      'Limit of',MAXCLS,' exceeded for operator:',
     &      PLABEL, PCOMB
         CALL QUIT('INDPXL: MAXCLS in dcbcls.h too small')
      END IF
      CLSINT(NPRPCLS) = CLASLB
      CLSCMB(NPRPCLS) = PCOMB
      CLSCAL(NPRPCLS) = .TRUE.
      IORDCL(NPRPCLS) = IORDER
      CALL XPLDEF(PLABEL,CLASLB,IORDER,PCOMB,INDX,
     &                WORK,LWORK,IPRINT)
      INDXPL = INDX
 9999 CONTINUE
      CALL QEXIT('INDXPL')
      RETURN
C
 1000 CONTINUE
 1010 CONTINUE
!radovan: rest of the routine is unreachable
      WRITE (LUPRI,'(/A/2A/2A)')
     &'INDXPL ERROR: illegal PCOMB or PDOINT',
     &'    PCOMB    = ',PCOMB,
     &'    PDOINT   = ',PDOINT(I)
      STOP 'INDXPL ERROR: illegal PCOMB or PDOINT(I)'
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getcls */
      SUBROUTINE GETCLS(LABINT,CLASLB,IORDER,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Given a file label for a given integral, this subroutine
C     returns the label for the general class of integrals
C     the integral belongs to. For moment integrals the order is returned
C     as well.
C
C     Jan-18-1995 - tsaue
C
C  2003, MI,HJAaJ extended list for integrals over LAO
C        for NMR shielding and magnetic susceptibilities
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      CHARACTER LABINT*8,CLASLB*7

      CALL QENTER('GETCLS')

C     Overlap integral
C
      IF(LABINT.EQ.'OVERLAP') THEN
        CLASLB = 'OVERLAP'
C
C     Dipole length/velocity integrals
C
      ELSE IF(LABINT(2:4).EQ.'DIP') THEN
        IF(LABINT(5:7).EQ.'VEL') THEN
          CLASLB = 'DIPVEL '
        ELSE
          CLASLB = 'DIPLEN '
        ENDIF
C
C     Quadrupole moments integrals
C
      ELSE IF(LABINT(3:8).EQ.'QUADRU') THEN
        CLASLB = 'QUADRUP'
C
C     Spin-orbit integrals
C
      ELSE IF(LABINT(3:8).EQ.'SPNORB') THEN
        CLASLB = 'SPNORB '
C          
C     SOPP integrals for ECP
C          
      ELSE IF(LABINT(3:6) .EQ.'SOPP') THEN
        CLASLB = 'SOPPINT'
C
C     Second moments integrals
C
      ELSE IF(LABINT(3:8).EQ.'SECMOM') THEN
        CLASLB = 'SECMOM '
C
C     Traceless theta quadrupole integrals
C
      ELSE IF(LABINT(3:7).EQ.'THETA') THEN
        CLASLB = 'THETA  '

!requested by Malaya Kumar Nayak <mk.nayak72@gmail.com>
      ELSE IF(LABINT(2:8).EQ.'NUCFIEL') THEN
        CLASLB = 'NUCFIEL'

C
C     Cartesian/spherical/electronic solvent moments integrals
C
      ELSE IF(LABINT(2:2).EQ.'M') THEN
        IF(LABINT(1:1).EQ.'C') THEN
          CLASLB = 'CARMOM '    
          READ(LABINT(3:4),'(I2)') NX
          READ(LABINT(5:6),'(I2)') NY
          READ(LABINT(7:8),'(I2)') NZ
          IORDER = NX + NY + NZ
        ELSE
          CLASLB = 'SPHMOM '
          READ(LABINT(3:4),'(I2)') IORDER
        ENDIF
C
C     One-electron Fermi contact integrals
C
      ELSE IF(LABINT(1:2).EQ.'FC') THEN
        CLASLB = 'FERMI C'
C
C     Paramagnetic spin-orbit integrals
C
      ELSE IF(LABINT(1:3).EQ.'PSO') THEN
        CLASLB = 'PSO    '
C     
C     Spin-dipole integrals
C
      ELSE IF(LABINT(1:2).EQ.'SD') THEN
C
C     Spin-dipole + Fermi contact integrals
C
        IF(LABINT(3:3).EQ.'C') THEN
          CLASLB = 'SDFC   '
        ELSE
C     
C     Spin-dipole integrals
C
          CLASLB = 'SPIN-DI'
        ENDIF
C
C     Diamagnetic spin-orbit integrals
C
      ELSE IF(LABINT(1:3).EQ.'DSO') THEN
        CLASLB = 'DSO    '
C
C     Half-derivative overlap integrals
C
      ELSE IF(LABINT(1:3).EQ.'HDO') THEN
        CLASLB = 'HDO   '
C
C     Contribution from overlap matrix to magnetic properties
C
      ELSE IF(LABINT(1:5).EQ.'dS/dB') THEN
        IF(LABINT(6:6).EQ.'2') THEN
C
C       Contribution from overlap matrix to magnetic properties
C
          CLASLB = 'S2MAG  '
        ELSE
C
C       Second order contribution from overlap matrix to 
C       magnetic properties
C
          CLASLB = 'S1MAG  '
        ENDIF
C
C     Angular momentum around the nuclei
C
      ELSE IF(LABINT(2:7).EQ.'ANGLON') THEN
        CLASLB = 'ANGLON '
C
C     Electronic angular momentum around origin
C
      ELSE IF(LABINT(2:7).EQ.'ANGMOM') THEN
        CLASLB = 'ANGMOM '
C
C     London orbital contribution to angular momentum
C
      ELSE IF(LABINT(2:7).EQ.'LONMOM') THEN
        CLASLB = 'LONMOM '
C
C     One-electron contributions to magnetic moment
C
      ELSE IF(LABINT(1:5).EQ.'dh/dB') THEN
        CLASLB = 'MAGMOM '
C
C     Electronic kinetic energy
C
      ELSE IF(LABINT(1:7).EQ.'KINENER') THEN   
        CLASLB = 'KINENER'
C
C     Diamagnetic susceptibility
C
      ELSE IF(LABINT(4:6).EQ.'SUS') THEN
        IF    (LABINT(7:8).EQ.'NL') THEN
C
C         Diamagnetic susceptibility without London contribution
C
          CLASLB = 'DSUSNOL'
        ELSEIF(LABINT(7:8).EQ.'LL') THEN
C
C         Angular London orbital contribution to diamagnetic
C         susceptibility
C
          CLASLB = 'DSUSLAN'
        ELSE
C
C       Angular London orbital contribution to diamagnetic 
C       susceptibility
C
          CLASLB = 'DSUSLH'
        ENDIF
C
C     Angular London orbital contribution to diamagnetic susceptibility
C
      ELSE IF(LABINT(3:8).EQ.'dh/dB2') THEN
        CLASLB = 'DIASUS '
C
C     Nuclear shielding integrals without London orbital contribution
C
      ELSE IF(LABINT(4:7).EQ.'NSNL') THEN
        CLASLB = 'NUCSNLO'
C
C     London orbital contribution to nuclear shielding tensor integrals
C
      ELSE IF(LABINT(4:7).EQ.'NSLO') THEN
        CLASLB = 'NUCSLO '
C
C     Nuclear shielding tensor integrals
CMI ... use this rather than 'NST'
      ELSE IF(LABINT(4:7).EQ.' NST') THEN
        CLASLB = 'NUCSHI '
C
C     Electric field at the individual nuclei
C
      ELSE IF(LABINT(1:3).EQ.'NEF') THEN
        CLASLB = 'NEFIELD'
C
C     traceless electric field third derivatives
C
      ELSE IF(LABINT(5:5).EQ.'V') THEN
        CLASLB = 'EFT    '
C
C     Parity violation - chirality operator
C
      ELSE IF(LABINT(1:3).EQ.'PVC') THEN
        CLASLB = 'PVC    '
C
C     Electric field gradient at the individual nuclei, spherical
C
      ELSE IF(LABINT(3:5).EQ.'FGS') THEN
        CLASLB = 'ELFGRDS'
C
C     Bra-differentiation of overlap matrix with respect to magnetic 
C     field
C
      ELSEIF(LABINT(1:7).EQ.'d<S|/dB') THEN
        CLASLB = 'S1MAGL'
C
C     Ket-differentiation of overlap matrix with respect to 
C     magnetic field
C
      ELSE IF(LABINT(1:7).EQ.'d|S>/dB') THEN
        CLASLB = 'S1MAGR'
C
C     Ket-differentation of HDO-integrals with respect to magnetic field
C
      ELSE IF(LABINT(4:6).EQ.'HDB') THEN
        CLASLB = 'HDOBR  '
C
C     Potential energy at the nuclei
C
      ELSE IF(LABINT(1:5).EQ.'POT.E') THEN
        CLASLB = 'NUCPOT '
C
C     Half B-differentiated overlap matrix
C
      ELSE IF(LABINT(2:5).EQ.'HBDO') THEN
        CLASLB = 'HBDO   '
C
C     Half-derivative overlap integrals not to be antisymmetrized
C
      ELSE IF(LABINT(1:5).EQ.'SQHDO') THEN
        CLASLB = 'SQHDO  '
C
C     Diamagnetic susceptibility with common gauge origin
C
      ELSE IF(LABINT(3:8).EQ.'SUSCGO') THEN
        CLASLB = 'DSUSCGO'
C
C     Nuclear shielding integrals with common gauge origin
C
      ELSE IF(LABINT(4:7).EQ.'NSCO') THEN
        CLASLB='NSTCGO '
C
C     Cosine and sine integrals
C
      ELSE IF(LABINT(2:7).EQ.'COSINE'.OR.LABINT(2:5).EQ.'SINE') THEN
        CLASLB = 'EXPIKR '
C
C     Mass velocity integrals
C
      ELSE IF(LABINT.EQ.'MASSVELO') THEN
        CLASLB = 'MASSVEL'
C
C     Darwin type integrals
C
      ELSE IF(LABINT(1:6).EQ.'DARWIN') THEN
        CLASLB = 'DARWIN '
C
C     First order magnetic field derivatives of electric field
C
      ELSE IF(LABINT(3:5).EQ.'CM1') THEN
        CLASLB = 'CM1    '
C
C     Second order magnetic field derivatives of electric field
C
      ELSE IF(LABINT(3:5).EQ.'CM2') THEN
        CLASLB = 'CM2    '
C
C     Half-derivative overlap integrals not to be anti-symmetrized
C
      ELSE IF(LABINT(1:6).EQ.'SQHDOR') THEN
        CLASLB = 'SQHDOR '
C
C
C
      ELSE IF(LABINT.EQ.'SQOVLAP') THEN
        CLASLB = 'SQOVLAP'
C
C     non-traceless electric field third derivatives
C
      ELSE IF(LABINT(5:5).EQ.'v') THEN
        CLASLB = 'EFTNTL '
C
C     Electric field gradient at the individual nuclei, cartesian
C
      ELSE IF(LABINT(3:5).EQ.'EFG') THEN
        CLASLB = 'ELFGRDC'
CMI
      ELSE IF(LABINT(2:6).EQ.'RM1H2') THEN
        CLASLB = 'RM1H2'
CMI
      ELSE IF(LABINT(3:7).EQ.'RM2H2') THEN
        CLASLB = 'RM2H2'
CMI
      ELSE IF(LABINT(2:6).EQ.'RM1H3') THEN
        CLASLB = 'RM1H3'
CMI
      ELSE IF(LABINT(3:7).EQ.'RM2H3') THEN
        CLASLB = 'RM2H3'

CMI ...contributions from RDSUSLL term
      ELSE IF(LABINT(3:7).EQ.'RDSUL') THEN
        CLASLB = 'RDSUSLL'
CMI
      ELSE IF(LABINT(2:6).EQ.'RM1RN') THEN
        CLASLB = 'RM1RN'

CMI  London orbital contribution to exp.value of nuclear shielding tensor integrals
      ELSE IF(LABINT(4:7).EQ.'RNST') THEN
        CLASLB = 'RM1N1H '

CMI  CAP integrals/labels to be used in the MOLTRA/PRPTRA input
      ELSE IF(LABINT(1:6).EQ.'CAP_RE') THEN
        CLASLB = 'CAP_RE  '
      ELSE IF(LABINT(1:6).EQ.'CAP_IM') THEN
        CLASLB = 'CAP_IM  '
      ELSE IF(LABINT(1:6).EQ.'CAP_OV') THEN
        CLASLB = 'CAP_OVL '
      ELSE IF(LABINT(2:7).EQ.'CAPD1R') THEN
        CLASLB = 'CAPD1R' 
      ELSE IF(LABINT(2:7).EQ.'CAPD1I') THEN
        CLASLB = 'CAPD1I' 
      ELSE IF(LABINT(2:7).EQ.'CAPDVE') THEN
        CLASLB = 'CAPDVE' 
      ELSE IF(LABINT(1:6).EQ.'CAPD_V') THEN
        CLASLB = 'CAPD_VD' 
CMI  EDM operator
      ELSE IF(LABINT(1:3).EQ.'EDM') THEN
        CLASLB = 'KINENER'
C
C     Unknown LABINT
C   
      ELSE
       CLASLB = 'Unknown'
      ENDIF

C  Control print out
      IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/,2A)')
     &  'GETCLS called for LABINT=',LABINT
        WRITE(LUPRI,'(7X,2A)')
     &  ' assigned  CLASLB=',CLASLB
      ENDIF
C
      CALL QEXIT('GETCLS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xpldef */
      SUBROUTINE XPLDEF(PLABEL,CLASLB,IORDER,PCOMB,INDX,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define property label: symmetries etc..
C
C     Called from: integ.function INDXPL
C
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
      CHARACTER PLABEL*8,CLASLB*7,PCOMB*4,LABINT(15*MXCENT)*8
      DIMENSION WORK(LWORK)
C
      CALL QENTER('XPLDEF')
#include "memint.h"
C
      MAXTYP = 15*MXCENT
      IF ((CLASLB.EQ.'ELFGRDC').OR.(CLASLB.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE
         NINTAD = MAXTYP
      END IF
C     Memory allocation
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KATOM ,NUCIND,WORK,KFREE,LFREE)
C
      CALL XPLDE1(PLABEL,CLASLB,PCOMB,IORDER,WORK(KINTRP),
     &            WORK(KINTAD),LABINT,INDX,WORK(KATOM),
     &            WORK(KFREE),LFREE,IPRINT)
C
      CALL QEXIT('XPLDEF')      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xplde1 */
      SUBROUTINE XPLDE1(PLABEL,CLASLB,PCOMB,IORDER,INTREP,INTADR,
     &           LABINT,INDX,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define property label: symmetries etc..
C
C    Entering: PLABEL
C
C
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "dcbprl.h"
#include "cbiher.h"
#include "nuclei.h"
C#include "maxorb.h"
C#include "symmet.h"

      LOGICAL ANTI, TRASPH, SOLVNT, SQUARE,LABERR,
     &        DOATOM(*)
      CHARACTER PLABEL*8,CLASLB*7,LABINT*8,PCOMB*4
      DIMENSION INTREP(*),INTADR(*),LABINT(*),WORK(LWORK)
      CALL QENTER('XPLDE1')
C
C
C     Get information on this class of integrals
C
      IF    (CLASLB.EQ.'MOLFLD ') THEN
        NOPTYP    = 1
        LABINT(1) = 'MOLFIELD'
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASLB.EQ.'FDEVEMB') THEN
        NOPTYP    = 1
        LABINT(1) = 'FDEVEMB '
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASLB.EQ.'BETAMAT') THEN
        NOPTYP    = 1
        LABINT(1) = 'BETAMAT '
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASLB .EQ. 'SPHMOM ') THEN
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        CALL SMTYP(IORDER,NOPTYP,INTREP,LABINT)
      ELSEIF(CLASLB .EQ. 'EFFDEN ') THEN
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        CALL SETATM(DOATOM,NATOM,-99)
        CALL EFFDENTYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM) 
      ELSE
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO

        IF (PLABEL.EQ.'EDM') THEN
          PLABEL='KINENERG'
          WRITE(LUPRI,*)
     &     'EDM: PLABEL SET to KINENERG, factor -4*c*field_strenght'
        ENDIF

        IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/,2X,A)') '** XPLDE1 before PR1DIR calling'
          WRITE(LUPRI,'(2X,A,I3)') 'INTTYP=',INTTYP
          WRITE(LUPRI,*) 'CLASLB=',CLASLB
          WRITE(LUPRI,*) 'PLABEL=',PLABEL
          WRITE(LUPRI,*) 'ANTI=',ANTI
          WRITE(LUPRI,*) 'SQUARE=',SQUARE
          WRITE(LUPRI,*) 'TRIANG=',TRIANG
        ENDIF

       ! ... call abacus routine...
        CALL PR1DIR(WORK,LWORK,CLASLB,INTTYP,NOPTYP,INTREP,
     &            ANTI,SQUARE,INTADR,LABINT,TRIANG,TRASPH,
     &            SOLVNT,IORDER,DOATOM,NATOM,NBAST,NELMNT,IPRINT)

        IF (IPRINT.GE.5) THEN
           WRITE(LUPRI,'(/,2X,A)')
     &      '*** XPLDE1 after PR1DIR calling ***'
           WRITE(LUPRI,'(2X,A,I3)') 'INTTYP=',INTTYP
           WRITE(LUPRI,*) 'CLASLB=',CLASLB
           WRITE(LUPRI,*) 'PLABEL=',PLABEL
           WRITE(LUPRI,*) 'ANTI=',ANTI
           WRITE(LUPRI,*) 'SQUARE=',SQUARE
           WRITE(LUPRI,*) 'TRIANG=',TRIANG
           WRITE(LUPRI,*) 'NOPTYP=',NOPTYP
           DO ICMP = 1,NOPTYP
             WRITE(LUPRI,*) 'LABINT(1-NOPTYP):',ICMP,LABINT(ICMP)
           ENDDO
        ENDIF
      ENDIF
C
C     Loop over components of operator class
C     ======================================   
C
      LABERR = .TRUE.
      DO ICMP = 1, NOPTYP
        NPRPLBL = NPRPLBL + 1
        IF(NPRPLBL.GT.MAXPRPLBL) THEN
          CALL QUIT('XPLDEF:Redimension MAXPRPLBL !')
        ENDIF
        IF(PLABEL.EQ.LABINT(ICMP)) THEN
        !mi: the only place to reset LABERR, more printout needed for "else" branch
          LABERR = .FALSE.
          INDX   = NPRPLBL
        ELSE
          IF (IPRINT.GE.7) THEN
            WRITE(LUPRI,*) 'PLABEL.NE.LABINT(ICMP), ICMP=',ICMP
            WRITE(LUPRI,*)
     &      'PLABEL=',PLABEL,' LABINT(ICMP)=',LABINT(ICMP)
          ENDIF
        ENDIF
        PRPLBL (NPRPLBL) = LABINT(ICMP)
        PDOINT (NPRPLBL) = PCOMB
        IPRLREP(NPRPLBL) = INTREP(ICMP)
!       .....  
        IF(SQUARE) THEN
C
C       Non - symmetric operator
C
          IPRLTYP(NPRPLBL) = 0
        ELSEIF(ANTI) THEN
C
C       Anti-symmetric operator
C
          IPRLTYP(NPRPLBL) = -1
        ELSE
C
C       Symmetric operator
C
          IPRLTYP(NPRPLBL) =  1
        ENDIF
      ENDDO
C
C     Property label is not among the listed labels
C
      IF(LABERR) THEN
        WRITE(LUPRI,'(/,A,A8,A)') '*** XPLDE1 ERROR: Property label ',
     &        PLABEL,' not found ! *****'
        WRITE(LUPRI,'(A,A7)') 'Class label      : ',CLASLB
        WRITE(LUPRI,'(A,A8)') '  Property labels: ',LABINT(1)
        WRITE(LUPRI,'(19X,A8)') (LABINT(I),I=2,NOPTYP)
        CALL QUIT('XPLDE1:  error - operator label not found')
      ENDIF

C
      CALL QEXIT('XPLDE1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpsym */
      SUBROUTINE PRPSYM(INDXPR,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     This subroutine finds the symmetry of the unique operator no.
C     INDXPR in COMMON /DCBPRP/.
C
C     The routine assumes Hermitian operators and will therefore
C     consider symmetric operators purely real and antisymmetric operators
C     purely imaginary. Non-symmetric operators results in error and
C     termination.
C   
C     If an operators involves more than one label, the logical variable
C     SAMSYM is used to check that the components have identical symmetry
C     under space and time.
C
C     Written by T.Saue, January 1995
C     Last revision: tsaue - May 14 1996
C
C    MI,HJAaJ(TEC),2003 - added changes for LAO-based  integrals
C     - routine assumed general operators which can be symmetric, antisymm. 
C       and square
C     - removing the program stop when indeterminate symmetry of operator
C     - print out
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "dcbprl.h"
#include "dcbxpr.h"
#include "symmet.h"
#include "dgroup.h"
#include "pgroup.h"
C
      LOGICAL DIFSYM,DIFTIM, SKIP_JM4
      DIMENSION ICREP(3),ICTIM(3)
#include "ibtfun.h"

      CALL QENTER('PRPSYM')
C
C     Investigate property label no. INDXPR
C
      ITYP     = IPRPTYP(INDXPR)
      NLBL     = MCMP(ITYP)
      DIFSYM   = .FALSE.
      DIFTIM   = .FALSE.
      IAIND    = JM4(1,ITYP)
      IBIND    = IPRPLBL(1,INDXPR)
C
C     Point group symmetry
C     ====================
C
C     Matrix operator
      IAREP    = JM4REP(IAIND)
C     Scalar operator
      IBREP    = IPRLREP(IBIND)
C     Total point gropu symmetry of the operator, first component..
      ICREP(1) = IBTXOR(IAREP,IBREP)
C
C     Time reversal symmetry
C     ======================
C
C     Matrix operator :
C - assign IATIM according to the 
C  IAIND = JM4(1,ITYP) = JM4(1,IPRPTYP(INDXPR)) = JM4(1,1-21)
C  JM4 is assigned in abacus/hergroup.F
C
C JM4(1,ITYP) = 0 ... scalar operator
C               4 ... gamma5 P
C               16..18 ????
C
C Miro: for determination of IATIm why not not use rather ITYP ?
C
      SKIP_JM4 = .TRUE.
!     SKIP_JM4 = .FALSE.

CMI ... catch possible exotic IAND...
      IF (IAIND.GE.16 .AND.IAIND.LE.18) THEN 
        CALL QUIT('quit here due to IAIND !')
      ENDIF

CMI  remove this part (THEN block) if everythink goes well and keep only ELSE branch
      IF (.NOT.SKIP_JM4) THEN

      IF (IAIND.EQ. 0 .OR. IAIND.EQ.4 .OR.
     &   (IAIND.GE.16 .AND.IAIND.LE.18)) THEN
C       DIAGONAL, GAMMA5, i*BETA*[XYZ]ALPHA
        IATIM = 1
      ELSE
C       [XYZ]SIGMA, [XYZ]ALPHA, but i*BETA*GAMMA5 has to be treated extra
        IATIM = -1
      ENDIF

CMI    exception for the i*BETA*GAMMA5 operator
      IF (ITYP.EQ.21) THEN
        IATIM=-1 !
         WRITE(LUPRI,'(2X,A)')
     &    'i*BETA*GAMMA5 operator, setting IATIM=-1' 
!    &    'FYI: i*BETA*GAMMA5 operator' 
      ENDIF

      ELSE
!  -------------------------------------------------------------------------------------------------
!   Simplification branch:
!       ...only  scalar , gamma5, and beta matrix operators have IATIM=+1, the rest has IATIM=-1
! FYI: only this part (ELSE block) is to be kept for the future
!  --------------------------------------------------------------------------------------------------
       IF (ITYP.EQ.1.OR.ITYP.EQ.9.OR.ITYP.EQ.19) THEN
        IATIM=+1
       ELSE ! the rest ....
        IATIM=-1
       ENDIF

      ENDIF

CMI - TODO force change of the time rev.symm. if wanted
!     IF (IPRPTYP(INDXPR).EQ.2.OR.
!    &   IPRPTYP(INDXPR).EQ.3.OR.IPRPTYP(INDXPR).EQ.4) THEN
!       IATIM=+1
!       WRITE(LUPRI,*) 'forcing IATIM=1 for IPRPTYP(INDXPR)!'
!     ENDIF


C     Scalar operator
      IBTIM    = IPRLTYP(IBIND)
C     Total time reversal symmetry of the 4-comp.operator
      ICTIM(1) = IATIM*IBTIM

      IF (IPRINT.GE.15) THEN
        CALL HEADER('Output from PRPSYM',-1)
        WRITE(LUPRI,'(2X,4(A,I2))')
     &  'IAIND=',IAIND,' IBIND=',IBIND,
     &  ' IATIM=',IATIM,' IBTIM=',IBTIM
        WRITE(LUPRI,'(2X,1(A,I2))')
     &  'ICTIM(1)=',ICTIM(1)
      ENDIF

C
C     Check against other components of operator
C     ==========================================
      DO 20 ILBL = 2,NLBL
        IAIND       = JM4(ILBL,ITYP)
        IBIND       = IPRPLBL(ILBL,INDXPR)
        IAREP       = JM4REP(IAIND)
        IBREP       = IPRLREP(IBIND)
        ICREP(ILBL) = IBTXOR(IAREP,IBREP)
        DIFSYM      = DIFSYM.AND.(ICREP(ILBL).NE.ICREP(1))
        IATIM       = -1
        IF(IAIND.EQ.0) IATIM = 1
        IBTIM       = IPRLTYP(IBIND)
        ICTIM(ILBL) = IATIM*IBTIM
        DIFTIM      = DIFTIM.AND.(ICTIM(ILBL).NE.ICTIM(1))
   20 CONTINUE
!     error message
      IF(DIFSYM.OR.DIFTIM) THEN
C
C         The components of the operator have not the same symmetry; abort
C
        CALL HEADER('PRPSYM: Abort !!!',-1)
        WRITE(LUPRI,'(A,I5,A)') 'Components of operator no.',INDXPR,
     &    ' have not the same symmetry :'
        WRITE(LUPRI,'(A,I5)') 'Operator type:',ITYP
        WRITE(LUPRI,'(32X,A12,23X,A4)')
     &     'boson irrep','trev'
        DO 30 ILBL = 1,NLBL
          IAREP       = JM4REP(JM4(ILBL,ITYP))
          IBREP       = IPRLREP(IPRPLBL(ILBL,INDXPR))
          WRITE(LUPRI,'(3X,A8,3(3X,A3,A1,I1,A1,3X,A4),I5)')
     &      PRPLBL(IPRPLBL(ILBL,INDXPR)),
     &      REP(IAREP),'(',(IAREP+1),')','X',
     &      REP(IBREP),'(',(IBREP+1),')','=',
     &      REP(ICREP(ILBL)),'(',(ICREP(ILBL)+1),')','  |  ',
     &      ICTIM(ILBL)
   30   CONTINUE
        WRITE(LUPRI,*) 'FYI: matrix oper, IAREP=JM4REP(:)=',JM4REP
        WRITE(LUPRI,*) '     scalar oper,IBREP=IPRLREP(:)=',IPRLREP
        CALL QUIT('PRPSYM: Components of different symmetry !!!')
      ENDIF

CMI   ... assign symmetries
      IPRPSYM(INDXPR) = ICREP(1)+1
      IPRPTIM(INDXPR) = ICTIM(1)

CMI ... control printout is always worth 
      IF (IPRINT.GE.15) THEN
        !CALL HEADER('Output from PRPSYM',-1)

        WRITE(LUPRI,'(2X,A,I3,A,I2,A,I2)')
     &  'Operator index, INDXPR =',INDXPR,
     &  ' and operator type (1-',NOPTYP,
     &  '), IPRPTYP(INDXPR) =',IPRPTYP(INDXPR)

        WRITE(LUPRI,'(2X,A,I2)')
     & 'symmetry of scalar operator under matrix transpose(-1,0,1),'//
     & ' IPRLTYP() =',IPRLTYP(INDXPR)
        WRITE(LUPRI,'(2X,A,I1)')
     &   'Number of operator components, NLBL (1-3) =',NLBL   
        WRITE(LUPRI,'(2X,3(A,I2))')
     &  'IAIND=',IAIND,' IATIM=',IATIM,' IBTIM=',IBTIM
        WRITE(LUPRI,'(20X,A,14X,A)')
     &     '--- boson irrep ---','--- time rever symm ---'
        DO ILBL = 1, NLBL
          IAREP       = JM4REP(JM4(ILBL,ITYP))
          IBREP       = IPRLREP(IPRPLBL(ILBL,INDXPR))
          WRITE(LUPRI,'(3X,A8,1X,A4,1X,3(3X,A3,A1,I1,A1,3X,A4),I5)')
     &      PRPLBL(IPRPLBL(ILBL,INDXPR)),PDOINT(IPRPLBL(ILBL,INDXPR)),
     &      REP(IAREP),'(',(IAREP+1),')','X',
     &      REP(IBREP),'(',(IBREP+1),')','=',
     &      REP(ICREP(ILBL)),'(',(ICREP(ILBL)+1),')','  |  ',
     &      ICTIM(ILBL)
        ENDDO
        WRITE(LUPRI,'(2X,A,I2,A,I2/)')
     & '* assigned IPRPSYM() =',IPRPSYM(INDXPR),
     & ' and time rev symm, IPRPTIM() =',IPRPTIM(INDXPR)
      ENDIF
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Output from PRPSYM',-1)
        CALL WRIXPR(0,INDXPR)
      ENDIF
C
      CALL QEXIT('PRPSYM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrixpr */
      SUBROUTINE WRIXPR(IND,INDXPR)
C***********************************************************************
C
C     Write out data on a operator specified by pointer INDXPR
C
C     Written by T.Saue - May 1996
C     Last revision: May 18 1996 - tsaue
C                    M.Ilias - Febr.2007 - new operator type 20
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0)
C
#include "dcbxpr.h"
#include "dcbprl.h"
#include "pgroup.h"
#include "dgroup.h"
      CHARACTER MXFORM*6,FMT*6,OPTYPE*12,TTYP(-1:1)*2
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
C
C     Operator name and symmetry
C     ==========================
C
      ISYM = IPRPSYM(INDXPR)
      ITIM = IPRPTIM(INDXPR)
      IF(IND.EQ.0) THEN
        WRITE(LUPRI,'(6X,A16,5X,A3,5X,A2)') 
     &    PRPNAM(INDXPR),REP(ISYM-1),TTYP(ITIM)
      ELSE
        WRITE(LUPRI,'(2X,I3,2X,A16,5X,A3,5X,A2)') 
     &   IND,PRPNAM(INDXPR),REP(ISYM-1),TTYP(ITIM)
      ENDIF
      CALL PRSYMB(LUPRI,'.',75,0)
C
C     Operator type
C     =============
C
      ITYP = IPRPTYP(INDXPR)
      IF(ITYP.EQ.0) THEN
        WRITE(LUPRI,'(6X,A,A8)') 'Operator read from matrix ',
     &        PRPLBL(IPRPLBL(1,INDXPR))
        FMT = MXFORM(ABS(FACPRP(1,INDXPR)),16)
        WRITE(LUPRI,'(6X,A,4X,'//FMT//')')
     &      'Factor     : ',FACPRP(1,INDXPR)
      GOTO 30
      ENDIF
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,
     &       13,14,15,16,17,18,19,20,21), ITYP
        WRITE(LUPRI,'(6X,A,I5)')
     &  'Program error: unknown operator type:',ITYP
        CALL QUIT('WRIXPR: unknown operator type')
    1 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'DIAGONAL : scalar operator'
      GOTO 200
    2 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XALPHA   : [alpha_x]P'
      GOTO 200
    3 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YALPHA   : [alpha_y]P'
      GOTO 200
    4 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZALPHA   : [alpha_z]P'
      GOTO 200
    5 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XAVECTOR : [alpha x P]_x'
      GOTO 200
    6 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YAVECTOR : [alpha x P]_y'
      GOTO 200
    7 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZAVECTOR : [alpha x P]_z'
      GOTO 200
    8 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ALPHADOT : [alpha.P]'
      GOTO 200
    9 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'GAMMA5   : [gamma_5]P'
      GOTO 200
   10 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XSIGMA   : [Sigma_x]P'
      GOTO 200
   11 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YSIGMA   : [Sigma_y]P'
      GOTO 200
   12 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZSIGMA   : [Sigma_z]P'
      GOTO 200
   13 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XBETASIG : [beta Sigma_x]P'
      GOTO 200
   14 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YBETASIG : [beta Sigma_y]P'
      GOTO 200
   15 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZBETASIG : [beta Sigma_z]P'
      GOTO 200
   16 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XiBETAAL : i[beta alpha_x]P'
      GOTO 200
   17 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YiBETAAL : i[beta alpha_y]P'
      GOTO 200
   18 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZiBETAAL : i[beta alpha_z]P'
      GOTO 200
   19 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'BETA     : scalar operator'
      GOTO 200
   20 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'SIGMADOT : [sigma.P]'
      GOTO 200
   21 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'BETAGAMMA5 : i[beta gamma5]P'
      GOTO 200
C
  200 CONTINUE
C
C     Labels and factors
C     ==================
C
      NLBL   = MCMP(ITYP)
      DO I = 1,NLBL
        FMT = MXFORM(ABS(FACPRP(I,INDXPR)),16)
        KLBL = IPRPLBL(I,INDXPR)
        IF    (IPRLTYP(KLBL).EQ. 1) THEN
          OPTYPE = ' (real)     '
        ELSEIF(IPRLTYP(KLBL).EQ.-1) THEN
          OPTYPE = ' (imaginary)'
        ELSE
          OPTYPE = ' (complex)  '
        ENDIF
        IF(I.EQ.1) THEN
          WRITE(LUPRI,'(6X,A,A8,1X,A4,4X,'//FMT//',A12)')
     &      'Labels and factors     : ',
     &      PRPLBL(IPRPLBL(I,INDXPR)),PDOINT(IPRPLBL(I,INDXPR)),
     &      FACPRP(I,INDXPR),OPTYPE
        ELSE
          WRITE(LUPRI,'(31X,A8,1X,A4,4X,'//FMT//',A12)')
     &      PRPLBL(IPRPLBL(I,INDXPR)),PDOINT(IPRPLBL(I,INDXPR)),
     &      FACPRP(I,INDXPR),OPTYPE
        ENDIF
      ENDDO
 30   CONTINUE
      CALL PRSYMB(LUPRI,'.',75,0)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1ind */
      SUBROUTINE OP1IND(OP1KEY,IND1OP,IPR1OP,N1OPER,INDXPR,MAX1OP)
C***********************************************************************
C
C     INDXPR points to an operator in XCBPRP. Search through
C     IPR1OP for elements pointing to this operator. If not found,
C     extend the list by one element.
C
C     Written by T.Saue June 28 1996
C     Last revision: July 8 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER OP1KEY*6
      DIMENSION IPR1OP(MAX1OP)
C
C     Search for INDXPR
C     =================
C
      DO J = 1,N1OPER
        IF(INDXPR.EQ.IPR1OP(J)) THEN
          IND1OP = J
          RETURN
        ENDIF
      ENDDO
C
C     Not found --> Add to list
C     =========================
C
      N1OPER = N1OPER + 1
      IND1OP = N1OPER
      IF(N1OPER.GT.MAX1OP) THEN
        WRITE(LUPRI,'(A6,A)') OP1KEY,' out of bounds in OP1IND'
        WRITE(LUPRI,'(A,I5)') 'Current value of MAX1OP: ',MAX1OP
        CALL QUIT('OP1IND: Out of bounds !')
      ELSE
        IPR1OP(IND1OP) = INDXPR
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1srt */
      SUBROUTINE OP1SRT(OP1KEY,IPOP1U,IPOP1S,NOP1T,NOP1,JOP1,
     &           WORK,LWORK)
C***********************************************************************
C
C     Sort a list of operators in symmetries
C
C     Written by T.Saue June 28 1996
C     Last revision - June 26 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER OP1KEY*6
      DIMENSION IPOP1U(*),IPOP1S(*),NOP1(8,2),JOP1(8,2),WORK(LWORK)
C
      CALL QENTER(OP1KEY)
#include "memint.h"
      NBUF = NOP1T*8
      CALL MEMGET('INTE',KIBUF,NBUF,WORK,KFREE,LFREE)
      CALL OP1SR1(IPOP1U,IPOP1S,NOP1T,NOP1,JOP1,WORK(KIBUF))
      CALL QEXIT (OP1KEY)      
C   
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1sr1 */
      SUBROUTINE OP1SR1(IPOP1U,IPOP1S,NOP1T,NOP1,JOP1,IBUF)
C***********************************************************************
C
C     Sort a list of operators in symmetries
C
C     Written by T.Saue June 28 1996
C     Last revision - June 26 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxpr.h"
#include "dgroup.h"
      DIMENSION IPOP1U(NOP1T),IPOP1S(NOP1T),NOP1(8),
     &          JOP1(8),IBUF(NOP1T,8)
C
C     Sort in buffer
C
      DO J = 1,NOP1T
        INDXPR = IPOP1U(J)
        ISYM   = IPRPSYM(INDXPR)
        NOP1(ISYM) = NOP1(ISYM) + 1
        IBUF(NOP1(ISYM),ISYM) = J
      ENDDO
C   
C     Pack results in IPOP1U
C
      NFC  = NBSYM/NFSYM
      IOFF = 0
      DO IFRP = 1,NFSYM
        DO IS = 1,NFC
          ISYM = JFSYM(IS,IFRP)
          NBUF = NOP1(ISYM)
          JOP1(ISYM) = IOFF
          IF(NBUF.GT.0) THEN
            DO J = 1,NBUF
              IPOP1S(IOFF+J) = IBUF(J,ISYM)
            ENDDO            
            IOFF = IOFF + NBUF
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck INIGRD */
      SUBROUTINE INIGRD
C***********************************************************************
C
C     Initialization for molecular gradient
C
C     Written by J. Thyssen - 1997/09/01
C     Last revision: jth - 1997/09/01
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
C
#include "maxorb.h"
#include "mxcent.h"
C
#include "dcbgrd.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbprp.h"
C
#include "ibtfun.h"
C
C     Initialize /CBIGRD/
C     ===================
C
      IPRGRD = IPRPRP
C
C     RELGRD is used in DHFGRD to control integrals
C
      RELGRD = .FALSE.
      DONGRD = .FALSE.
      IGRD_INTFLG = INTGEN
      IF (LEVYLE) IGRD_INTFLG = IBTAND ( 1, IGRD_INTFLG )
      DOTRCK    = .FALSE.
      SCRGRD    = -1.0D0      
C
      RETURN
      END
C***********************************************************************
C  /* Deck grdinp */
      SUBROUTINE GRDINP(WORD,RESET)
C***********************************************************************
C
C     Input section for molecular gradient
C
C     Written by J. Thyssen - 1997/07/09
C     Last revision: jth - 1997/07/09
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0=0.0D0)
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 5)
C
#include "dcbham.h"
#include "dcbgrd.h"
#include "dcbprp.h"
C
      LOGICAL SET, NEWDEF, RESET, LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &          PNAME*16, PLABEL(3)*8, STR*3
      DIMENSION PFAC(3),II(3)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.INTFLG','.TRICK ','.SCREEN','.NUMGRA'/
      DATA SET/.FALSE./
C
      NEWDEF = (WORD .EQ. '*MOLGRD')
      IF (SET) THEN
         IF (NEWDEF)
     &      CALL QUIT('Only one "*MOLGRD" input section allowed')
C        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
      CALL INIGRD
C
C
C     Process input for DCBEXP
C     ========================
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
                     GO TO (1,2,3,4,5), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in GRDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in GRDINP.')
    1          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRGRD
               GO TO 100
    2          CONTINUE
C&&&& INTFLG: Which two-electron gradient to calculate
                  IF (GAUNT) THEN 
                     READ(LUCMD,*) ILL,ISL,ISS,IGT
                  ELSE
                     READ(LUCMD,*) ILL,ISL,ISS
                  END IF
                  READ(LUCMD,*) ILL, ILS, ISS
                  IGRD_INTFLG = ILL + 2 * ILS + 4 * ISS + 8 * IGT
               GO TO 100
    3          CONTINUE
C&&&& TRICK : Skip LS- and SS- integrals if contribution negligible
                  DOTRCK = .TRUE.
               GO TO 100
    4          CONTINUE
C&&&& SCREEN: Threshold for screening on the gradient
                  READ (LUCMD,*) SCRGRD
               GO TO 100
    5          CONTINUE
C&&&& NUMGRA: Use numerical gradient
                  DONGRD = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in GRDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in GRDINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
      IF (MGRAD) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(A)')
     &        ' GRDINP: Calculation of molecular gradient.'
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A,I3)')
     &        '* Print level                       : ', IPRGRD
         IF ( DONGRD ) THEN
            WRITE(LUPRI,'(1X,A)')
     &         '* Gradient evaluated with numerical differentiation.'
         ELSE
            if (bss) then
              call quit('analytical molecular '//
     &                  'gradient not implemented for BSS/X2C/IOTC')
            end if
            WRITE(LUPRI,'(1X,A)')
     &         '* Gradient evaluated analytically.'
         END IF
         WRITE(LUPRI,'(1X,A)')
     &        '* Contributions to the molecular gradient from'
         IF ( LBIT( IGRD_INTFLG, 1 ) )
     &         WRITE(LUPRI,'(3X,A)') '- LL-integrals.'
         IF ( LBIT( IGRD_INTFLG, 2 ) ) THEN
            IF ( DOTRCK ) THEN
               WRITE(LUPRI,'(3X,A)')
     &            '- LS-integrals (skipped if estimated to be small).'
            ELSE
               WRITE(LUPRI,'(3X,A)') '- LS-integrals.'
            END IF
         END IF
         IF ( LBIT( IGRD_INTFLG, 3 ) ) THEN
            IF ( DOTRCK ) THEN
               WRITE(LUPRI,'(3X,A)')
     &            '- SS-integrals (skipped if estimated to be small).'
            ELSE
               WRITE(LUPRI,'(3X,A)') '- SS-integrals.'
            END IF
         END IF
         IF ( LBIT( IGRD_INTFLG, 4 ) ) THEN
            IF ( DOTRCK ) THEN
               WRITE(LUPRI,'(3X,A)')
     &            '- GT-integrals (skipped if estimated to be small).'
            ELSE
               WRITE(LUPRI,'(3X,A)') 
     &         '- GT-integrals (WARNING: gradient not tested yet)'
            END IF
         END IF
         IF(SCRGRD.GT.D0) THEN
            WRITE(LUPRI,'(A,1P,E8.2)')
     &         ' * Screening threshold in gradient calculation ',SCRGRD
         ELSE
            WRITE(LUPRI,'(A)')
     &         ' * No screening in the calculation of '//
     &         'two-electron integrals.'
         END IF
         CALL PRSYMB(LUPRI,'-',75,0)
      END IF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck minprp */
      SUBROUTINE MINPRP(WORK,LWORK)
C***********************************************************************
C     Calculate specified properties at final, converged geometry
C     (called when .OPTIMIZE specified if geometry converges)
C
#include "implicit.h"
#include "dcbxpr.h"
#include "priunit.h"
#include "dcbgen.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL FOPEN
C
      CALL QENTER('MINPRP')
C
      NPRPS = 0
C
C
C     Open LU1INT if not opened
C
      IF (DOPRP) THEN
         INQUIRE(FILE='AOPROPER',OPENED=FOPEN)
         IF (.NOT. FOPEN) THEN
#if defined (SYS_IBM)
            OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
            OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
         END IF
C
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PRPINP('**PRP F',WORK,LWORK)
         CLOSE(LUCMD)
C
         CALL ONEGEN(WORK,LWORK)
C
C        Close AOPROPER if it was closed on entry to
C        MINPRP.
C
         IF (.NOT. FOPEN) CLOSE(LU1INT,STATUS = 'KEEP')
C
C
C***********************************************************************
C*****  R E S P O N S E    M O D U L E  ********************************
C***********************************************************************
C
         CALL PAMPRP(WORK,LWORK)
C
      END IF
C
      CALL QEXIT('MINPRP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
