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
C VERSION : $Revision: 1.44.2.3 $
C DATE    : $Date: 2004/09/01 08:59:43 $
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
C       Written by Trond Saue and Hans J.Aa.Jensen
C       Last revision: June 14 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 8, NTABLE = 25)
#include "maxorb.h"
#include "mxcent.h"
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7, WRDSRC*7
      DIMENSION WORK(LWORK)
C
#include "orgcom.h"
#include "dcbprp.h"
C
      DATA TABDIR /'*END OF','*EXPECT','*LINEAR','*EXCITA',
     &             '*MOLGRD','*NMR   ','*ESR   ','*QUADRA'/
      DATA TABLE  /'.DIPOLE','.QUADRU','.POLARI','.MAGNET',
     &             '.NMR   ','.SPIN-S','.SHIELD','.PRINT ',
     &             '.ABUNDA','.EFG   ','.NQCC  ','.DSO   ',
     &             '.MOLGRD','.NSTDIA','.ESR   ','.PVC   ',
     &             '.EFT   ','.EFTNTL','.PVCSHI','.EFN   ',
     &             '.WAVE F','.RHONUC','.FIRST ','.TWO-PH',
     &             '.VERDET'/
C
      CALL QENTER('PRPINP')
#include "memint.h"
C
      NJOBS = 0
C
C     Initialize /CBIPRP/
C     ===================
C
      IPRPRP   = 0
      DOEXP    = .FALSE.
      DOESR    = .FALSE.
      DOXLR    = .FALSE.
      DOXPP    = .FALSE.
      DOXQR    = .FALSE.
      DOTPA    = .FALSE.
      DOVER    = .FALSE.
      DIPOLE   = .FALSE.
      QUADRU   = .FALSE.
      EPOLAR   = .FALSE.
      MSUSCP   = .FALSE.
      SHIELD   = .FALSE.
      PVCSHI   = .FALSE.
      SPNSPN   = .FALSE.
      EFG      = .FALSE.
      EFN      = .FALSE.
      EFT      = .FALSE.
      EFTNTL   = .FALSE.
      NQCC     = .FALSE.
      DSO      = .FALSE.
      NSTDIA   = .FALSE.
      PVC      = .FALSE.
      RHONUC   = .FALSE.
      ABUND    = D1
      NDSOINT  = 0
      MGRAD    = .FALSE.
      ESRGTENS = .FALSE.
      ESRANISO = .FALSE.
      NPRPWF   = 1
      PRPWF(1) = 'DHF'
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
     &                117,118,119,120,121,122,123,124,125),I
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
            DOXLR  = .TRUE.
            MSUSCP = .TRUE.
C            CALL DEF_MSUSCP(WORK,LWORK)   
            GO TO 100
  105    CONTINUE
C&&&& NMR parameters:
            DOXLR  = .TRUE.
            SHIELD = .TRUE.
C            CALL DEF_SHIELD(WORK,LWORK)
            SPNSPN = .TRUE.
C            CALL DEF_SPNSPN(WORK,LWORK)                  
            GO TO 100
  106    CONTINUE
C&&&& SPIN-S: NMR indirect spin-spin coupling constants:
            DOXLR  = .TRUE.
            SPNSPN = .TRUE.
C            CALL DEF_SPNSPN(WORK,LWORK)                  
            GO TO 100
  107    CONTINUE
C&&&& SHIELD: NMR shielding constants:
            DOXLR = .TRUE.
            SHIELD = .TRUE.
C            CALL DEF_SHIELD(WORK,LWORK)
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
C&&&& SHIELD: NMR shielding constants with PVC contributions
            DOXLR = .TRUE.
            SHIELD = .TRUE.
            PVCSHI = .TRUE.
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
#if defined (MOD_XPP)
            CALL DEF_TPA(WORK,LWORK,DOTPA)
#else
            CALL QUIT('TPA not included in this version')
#endif
            GO TO 100
C&&&& .VERDET - Calculates the Verdet constant with the QR module
 125     CONTINUE
            DOVER = .TRUE.
            DOXQR = .TRUE.
            CALL DEF_VERDET(WORK,LWORK)
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
      WRITE (LUPRI,'(1X,A,3F18.10)')
     &     '* Operator center:', (ORIGIN(I),I=1,3),
     &     '* Gauge origin   :', (GAGORG(I),I=1,3),
     &     '* Dipole origin  :', (DIPORG(I),I=1,3)
C
C     Process input for various program sections
C     ==========================================
C
C     hjaaj Oct 2003:
C        first initialize common /dcbexp/
C        because it may be modified both by EXPINP
C        and elsewhere for diamag. terms for magnetic properties
C
      CALL EXPINP('INITIALIZE',.TRUE.,WORK,LWORK)
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
               GO TO (1,2,3,4,5,6,7,8), I
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
#if defined (MOD_XPP)
        CALL XPPINP(WORD,.FALSE.,WORK,LWORK)
#else
        CALL QUIT('*EXCITA not included in this version')
#endif
        GO TO 200
    5 CONTINUE
        CALL GRDINP(WORD,.FALSE.)
        GO TO 200
    6 CONTINUE
        CALL NMRINP(WORD,.FALSE.)
        GO TO 200
    7 CONTINUE
#if defined (MOD_ESR)
        CALL ESRINP(WORD,.FALSE.,WORK,LWORK)
#else
        CALL QUIT('*ESR    not included in this version')
#endif
        GO TO 200
    8 CONTINUE
         CALL XQRINP(WORD,.FALSE.,WORK,LWORK)
         GO TO 200
C
    1 CONTINUE
      WORD = 'INITLZE'
      IF(DOEXP) CALL EXPINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOXLR) CALL XLRINP(WORD,.TRUE.,WORK,LWORK)
      CALL NMRINP(WORD,.TRUE.)
      CALL GRDINP(WORD,.TRUE.)
#if defined (MOD_XPP)
      IF(DOXPP) CALL XPPINP(WORD,.TRUE.,WORK,LWORK)
#endif
      IF(DOXQR) CALL XQRINP(WORD,.TRUE.,WORK,LWORK)
#if defined (MOD_ESR)
      IF(DOESR) CALL ESRINP(WORD,.TRUE.,WORK,LWORK)
#endif
C
C     Define special operators
C
      CALL PRPDEF(WORK,LWORK)
C
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
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
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
      IF(EFG) THEN    
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
C    2. XALPHA
C       [alpha_x]P    * x-component of alpha times scalar operator
C    3. YALPHA
C       [alpha_y]P    * y-component of alpha times scalar operator
C    4. ZALPHA
C       [alpha_z]P    * z-component of alpha times scalar operator
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
C
#include "implicit.h"
#include "priunit.h"
      CHARACTER*7 WORD
C
      PARAMETER (NTABKW = 23)
      CHARACTER*8 TABKW(NTABKW)
      CHARACTER   LINE*80, KEYWD*8, TOKEN16*16, BLANK8*8, BLANK16*16
      DIMENSION PFAC(3),WORK(LWORK)
      CHARACTER   PLABEL(3)*8, PNAME*16
      LOGICAL     CLABEL
C
      PARAMETER  (D1 = 1.0D0)
      PARAMETER (BLANK8='        ', BLANK16='                ' )
      SAVE TABKW
      DATA TABKW/'CLASS   ', 'FACTORS ', 'COMFACTO', 'DIAGONAL',
     &           'XALPHA  ', 'YALPHA  ', 'ZALPHA  ', 'XAVECTOR',
     &           'YAVECTOR', 'ZAVECTOR', 'ALPHADOT', 'GAMMA5  ',
     &           'XSIGMA  ', 'YSIGMA  ', 'ZSIGMA  ', 'XBETASIG',
     &           'YBETASIG', 'ZBETASIG', 'XiBETAAL', 'YiBETAAL',
     &           'ZiBETAAL', 'CFACTORS', 'COMCFACT'/
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
C
      IKW = 0
C
C     *****************
C     *** Read line ***
C     *****************
C
  100 READ (LUCMDX,'(A)',END=8000) LINE
      IF (IPRINT .GT. 100) THEN
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
           IF (I .GT. 3) THEN
              IPTYP = I - 3
              NLABEL= 0
              NFAC  = 0
              IKW = 2
           END IF
           GO TO (4,5,6,1,1,1,1,2,2,2,3,1,1,1,1,1,1,1,1,1,1), I
         END IF
         END DO
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')
     &     ' XPRINP: INPUT ERROR following '//WORD,
     &     ' Keyword "'//KEYWD(1:LKW)//'" not recognized.'
         GO TO 300
C&&&&    One component: DIAGONAL, [XYZ]ALPHA, GAMMA5, [XYZ]SIGMA,
C&&&&                   [XYZ]BETASIG, [XYZ]iBETAAL
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
C&&&&    Class keyword
    4    CONTINUE
           CLABEL = .TRUE.
           IKW = 3
         GO TO 300
C&&&&    Read a multiplicative factor for each component
    5    CONTINUE
           IKW  = 4
         GO TO 300
C&&&&    Read a common factor for each component
    6    CONTINUE
           IKW  = 5
         GO TO 300
C
C     Read property labels for each component
C     =======================================
C
      ELSE IF (IKW .EQ. 2) THEN
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
CHJMAERKE: implement interpretation of class keyword
C
C     Read factor for each component
C     ==============================
C
      ELSE IF (IKW .EQ. 4) THEN
         NFAC = NFAC + 1
         JTKN = ITKN - 1 + LTKN
         I    = 16 - LTKN + 1
         TOKEN16 = BLANK16
         TOKEN16(I:16) = LINE(ITKN:JTKN)
         READ(TOKEN16,'(F16.0)') PFAC(NFAC)
         IF (NFAC .GE. NPCOMP) IKW = 1
C
C     Read a common factor for all components
C     =======================================
      ELSE IF (IKW .EQ. 5) THEN
         JTKN = ITKN - 1 + LTKN
         I    = 16 - LTKN + 1
         TOKEN16 = BLANK16
         TOKEN16(I:16) = LINE(ITKN:JTKN)
         READ(TOKEN16,'(F16.0)') COMFAC
         PFAC(1) = COMFAC*PFAC(1)
         PFAC(2) = COMFAC*PFAC(2)
         PFAC(3) = COMFAC*PFAC(3)
         IKW = 1
C
C     Unknown IKW
C     ===========
C
      ELSE
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
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xprind */
      SUBROUTINE XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
C*********************************************************************** 
C  Given a 4-component operator this routine returns an index
C  INDXPR to the list of operators in COMMON block XCBPRP and
C  returns its symmetry under time reversal(ITRXPR) and its
C  spatial symmetry(ISYXPR)
C
C  Copyright (c) 27-Jan-1995 Hans Joergen Aa. Jensen
C  Revised May 14 1996 - tsaue
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
      CALL QENTER('XPRIND')
      INDXPR = -1
C
C     Search for this property operator in existing list
C
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
C
C     This property operator has not been defined previously.
C     Append the new property operator to /XCBPRP/ list.
C
      CALL XPRADD(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,INDXPR,
     &            WORK,LWORK,IPRINT)
C
 1000 CONTINUE
      ISYXPR = IPRPSYM(INDXPR)
      ITRXPR = IPRPTIM(INDXPR)
      CALL QEXIT('XPRIND')
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER PNAME*16, PLABEL(NPCOMP)*8
      DIMENSION PFAC(NPCOMP),WORK(LWORK)
      LOGICAL DELSS
C
C Used from common blocks:
C  XCBPRP: everything
C  XCBPRI: LUPRI,NINFO
C
#include "dcbxpr.h"
#include "dcbham.h"
C
      CHARACTER*4 PCOMB
C
      CALL QENTER('XPRADD')
C
C     This property operator has not been defined previously.
C     Append the new property operator to /XCBPRP/ list.
C
      IF (PLABEL(1)(1:7).EQ.'BETAMAT') THEN
Chj      PCOMB = 'FFFT'
         PCOMB = '000+'
      ELSE
         IF (IPTYP.LT.1 .OR. IPTYP.GT.18)
     &      CALL QUIT('Illegal IPTYP')
         GO TO (1,2,2,2,2,2,2,2,2,
     &          1,1,1,3,3,3,4,4,4),IPTYP
C
C        DIAGONAL and SIGMA:
C
    1    PCOMB = '+00+'
Chj      PCOMB = 'TFFT'
Chj      IF (LEVYLE) PCOMB = 'TFFF'
CJTH DEBUG!!
#ifdef UNDEF
         write(6,*) 'debug debug: operator ',PNAME
         IF (PNAME(1:2) .EQ. 'EF') THEN
Chj         PCOMB = 'TFFF'
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
C        BETA and BETA SIGMA:
C
    3    PCOMB = '+00-'
         GO TO 100
C
C        i BETA ALPHA and i BETA GAMMA5:
C
    4    PCOMB = '0+-0'
         GO TO 100
C
  100    CONTINUE
         IF (LEVYLE) THEN
           DELSS=.FALSE.
           IF (PCOMB(4:4).NE.'0') DELSS=.TRUE.
           PCOMB(4:4) = '0'
         ENDIF
      END IF

      IF (IPRINT .GE. 7) THEN
         WRITE (LUPRI,'(/A/2A,/A,A4,/A,I5/A,I5/A,1P,3(D20.10))')
     &   'XPRADD called, property operator parameters:',
     &   '  PNAME   =',PNAME,
     &   '  PCOMB   =',PCOMB,
     &   '  IPTYP   =',IPTYP,
     &   '  NPCOMP  =',NPCOMP,
     &   '  PFAC    =',(PFAC(I),I=1,NPCOMP)
         WRITE (LUPRI,'(A,3(A8,2X))')
     &   '  PLABEL  =',(PLABEL(I),I=1,NPCOMP)
         IF (DELSS) WRITE(LUPRI,'(2X,A)')
     &    'Levy-Leblond nonrelativistic mode, SS block zeroed!'
      END IF

C
C     Define property
C
      NPRPS  = NPRPS + 1
      IF(NPRPS.GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &   ' ** ERROR in XPRADD **  Redimension MAXPRPS to ',NPRPS
        CALL QUIT('XPRADD: Too may properties !')
      ENDIF
      INDXPR = NPRPS
C     Property name
      PRPNAM(NPRPS)  = PNAME
C     Property type
      IPRPTYP(NPRPS) = IPTYP
C     Pointers to property labels
      DO 600 J = 1,NPCOMP
         FACPRP(J,NPRPS)  = PFAC(J)
         INDX = INDXPL(PLABEL(J),PCOMB,WORK,LWORK,IPRINT)
         IF (INDX .LE. 0) THEN
            INDXPR = -1
            WRITE (LUPRI,'(3A)')
     &      'PLABEL, PCOMB not found in /XCBPRL/ : ', PLABEL(J), PCOMB
         END IF
         IPRPLBL(J,NPRPS) = INDX
  600 CONTINUE
      IF (INDXPR .EQ. -1) CALL QUIT('XPRADD: PLABEL, PCOMB not found')
C     Symmetry under time reversal and point group symmetry
      CALL PRPSYM(INDXPR,IPRINT)
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
      CHARACTER PLABEL*8, PCOMB*4,CLASSLB*7
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
        CLASSLB = 'MOLFLD '
        IORDER  = 0
      ELSEIF(PLABEL(1:7).EQ.'BETAMAT') THEN
        CLASSLB = 'BETAMAT'
        IORDER   = 0
      ELSE
        CALL GETCLS(PLABEL,CLASSLB,IORDER)
        IF(CLASSLB.EQ.'Unknown') THEN
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
      CLSINT(NPRPCLS) = CLASSLB
      CLSCMB(NPRPCLS) = PCOMB
      CLSCAL(NPRPCLS) = .TRUE.
      IORDCL(NPRPCLS) = IORDER
      CALL XPLDEF(PLABEL,CLASSLB,IORDER,PCOMB,INDX,
     &                WORK,LWORK,IPRINT)
      INDXPL = INDX
 9999 CONTINUE
      CALL QEXIT('INDXPL')
      RETURN
C
 1000 CONTINUE
 1010 CONTINUE
      WRITE (LUPRI,'(/A/2A/2A)')
     &'INDXPL ERROR: illegal PCOMB or PDOINT',
     &'    PCOMB    = ',PCOMB,
     &'    PDOINT   = ',PDOINT(I)
      STOP 'INDXPL ERROR: illegal PCOMB or PDOINT(I)'
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getcls */
      SUBROUTINE GETCLS(LABINT,CLASSLB,IORDER)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Given a file label for a given integral, this subroutine
C     returns the label for the general class of integrals
C     the integral belongs to. For moment integrals the order is returned
C     as well.
C
C     Jan-18-1995 - tsaue
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      CHARACTER LABINT*8,CLASSLB*7
C
C     Overlap integral
C
      IF(LABINT.EQ.'OVERLAP') THEN
        CLASSLB = 'OVERLAP'
        RETURN
      ENDIF
C
C     Dipole length/velocity integrals
C
      IF(LABINT(2:4).EQ.'DIP') THEN
        IF(LABINT(5:7).EQ.'VEL') THEN
          CLASSLB = 'DIPVEL '
        ELSE
          CLASSLB = 'DIPLEN '
        ENDIF
        RETURN
      ENDIF
C
C     Quadrupole moments integrals
C
      IF(LABINT(3:8).EQ.'QUADRU') THEN
        CLASSLB = 'QUADRUP'
        RETURN
      ENDIF
C
C     Spin-orbit integrals
C
      IF(LABINT(3:8).EQ.'SPNORB') THEN
        CLASSLB = 'SPNORB '
        RETURN
      ENDIF
C
C     Second moments integrals
C
      IF(LABINT(3:8).EQ.'SECMOM') THEN
        CLASSLB = 'SECMOM '
        RETURN
      ENDIF
C
C     Traceless theta quadrupole integrals
C
      IF(LABINT(3:7).EQ.'THETA') THEN
        CLASSLB = 'THETA  '
        RETURN
      ENDIF
C
C     Cartesian/spherical/electronic solvent moments integrals
C
      IF(LABINT(2:2).EQ.'M') THEN
        IF(LABINT(1:1).EQ.'C') THEN
          CLASSLB = 'CARMOM '    
          READ(LABINT(3:4),'(I2)') NX
          READ(LABINT(5:6),'(I2)') NY
          READ(LABINT(7:8),'(I2)') NZ
          IORDER = NX + NY + NZ
        ELSE
          CLASSLB = 'SPHMOM '
          READ(LABINT(3:4),'(I2)') IORDER
        ENDIF
        RETURN
      ENDIF
C
C     One-electron Fermi contact integrals
C
      IF(LABINT(1:2).EQ.'FC') THEN
        CLASSLB = 'FERMI C'
        RETURN
      ENDIF
C
C     Paramagnetic spin-orbit integrals
C
      IF(LABINT(1:3).EQ.'PSO') THEN
        CLASSLB = 'PSO    '
        RETURN
      ENDIF
C     
C     Spin-dipole integrals
C
      IF(LABINT(1:2).EQ.'SD') THEN
C
C     Spin-dipole + Fermi contact integrals
C
        IF(LABINT(3:3).EQ.'C') THEN
          CLASSLB = 'SDFC   '
          RETURN
        ELSE
C     
C     Spin-dipole integrals
C
          CLASSLB = 'SPIN-DI'
        ENDIF
        RETURN
      ENDIF
C
C     Diamagnetic spin-orbit integrals
C
      IF(LABINT(1:3).EQ.'DSO') THEN
        CLASSLB = 'DSO    '
        RETURN
      ENDIF
C
C     Half-derivative overlap integrals
C
      IF(LABINT(1:3).EQ.'HDO') THEN
        CLASSLB = 'HDO   '
        RETURN
      ENDIF
C
C     Contribution from overlap matrix to magnetic properties
C
      IF(LABINT(1:5).EQ.'dS/dB') THEN
        IF(LABINT(6:6).EQ.'2') THEN
C
C       Contribution from overlap matrix to magnetic properties
C
          CLASSLB = 'S2MAG  '
        ELSE
C
C       Second order contribution from overlap matrix to 
C       magnetic properties
C
          CLASSLB = 'S1MAG  '
        ENDIF
        RETURN
      ENDIF
C
C     Angular momentum around the nuclei
C
      IF(LABINT(2:7).EQ.'ANGLON') THEN
        CLASSLB = 'ANGLON '
        RETURN
      ENDIF
C
C     Electronic angular momentum around origin
C
      IF(LABINT(2:7).EQ.'ANGMOM') THEN
        CLASSLB = 'ANGMOM '
        RETURN
      ENDIF
C
C     London orbital contribution to angular momentum
C
      IF(LABINT(2:7).EQ.'LONMOM') THEN
        CLASSLB = 'LONMOM '
        RETURN
      ENDIF
C
C     One-electron contributions to magnetic moment
C
      IF(LABINT(1:5).EQ.'dh/dB') THEN
        CLASSLB = 'MAGMOM '
        RETURN
      ENDIF
C
C     Electronic kinetic energy
C
      IF(LABINT(1:7).EQ.'KINENER') THEN   
        CLASSLB = 'KINENER'
        RETURN
      ENDIF
C
C     Diamagnetic susceptibility
C
      IF(LABINT(4:6).EQ.'SUS') THEN
        IF    (LABINT(7:8).EQ.'NL') THEN
C
C         Diamagnetic susceptibility without London contribution
C
          CLASSLB = 'DSUSNOL'
        ELSEIF(LABINT(7:8).EQ.'LL') THEN
C
C         Angular London orbital contribution to diamagnetic
C         susceptibility
C
          CLASSLB = 'DSUSLAN'
        ELSE
C
C       Angular London orbital contribution to diamagnetic 
C       susceptibility
C
          CLASSLB = 'DSUSLH'
        ENDIF
        RETURN
      ENDIF
C
C     Angular London orbital contribution to diamagnetic susceptibility
C
      IF(LABINT(3:8).EQ.'dh/dB2') THEN
        CLASSLB = 'DIASUS '
        RETURN
      ENDIF
C
C     Nuclear shielding integrals without London orbital contribution
C
      IF(LABINT(3:6).EQ.'NSNL') THEN
        CLASSLB = 'NUCSNLO'
        RETURN
      ENDIF
C
C     London orbital contribution to nuclear shielding tensor integrals
C
      IF(LABINT(3:6).EQ.'NSLO') THEN
        CLASSLB = 'NUCSLO '
        RETURN
      ENDIF
C
C     Nuclear shielding tensor integrals
C
      IF(LABINT(4:6).EQ.'NST') THEN
        CLASSLB = 'NUCSHI '
        RETURN
      ENDIF
C
C     Electric field at the individual nuclei
C
      IF(LABINT(1:3).EQ.'NEF') THEN
        CLASSLB = 'NEFIELD'
        RETURN
      ENDIF
C
C     traceless electric field third derivatives
C
      IF(LABINT(5:5).EQ.'V') THEN
        CLASSLB = 'EFT    '
        RETURN
      ENDIF
C
C     Parity violation - chirality operator
C
      IF(LABINT(1:3).EQ.'PVC') THEN
        CLASSLB = 'PVC    '
        RETURN
      ENDIF
C
C     Electric field gradient at the individual nuclei, spherical
C
      IF(LABINT(3:5).EQ.'FGS') THEN
        CLASSLB = 'ELFGRDS'
        RETURN
      ENDIF
C
C     Bra-differentiation of overlap matrix with respect to magnetic 
C     field
C
      IF(LABINT(1:7).EQ.'d<S|/dB') THEN
        CLASSLB = 'S1MAGL'
        RETURN
      ENDIF
C
C     Ket-differentiation of overlap matrix with respect to 
C     magnetic field
C
      IF(LABINT(1:7).EQ.'d|S>/dB') THEN
        CLASSLB = 'S1MAGR'
        RETURN
      ENDIF
C
C     Ket-differentation of HDO-integrals with respect to magnetic field
C
      IF(LABINT(4:6).EQ.'HDB') THEN
        CLASSLB = 'HDOBR  '
        RETURN
      ENDIF
C
C     Potential energy at the nuclei
C
      IF(LABINT(1:5).EQ.'POT.E') THEN
        CLASSLB = 'NUCPOT '
        RETURN
      ENDIF
C
C     Half B-differentiated overlap matrix
C
      IF(LABINT(2:5).EQ.'HBDO') THEN
        CLASSLB = 'HBDO   '
        RETURN
      ENDIF
C
C     Half-derivative overlap integrals not to be antisymmetrized
C
      IF(LABINT(1:5).EQ.'SQHDO') THEN
        CLASSLB = 'SQHDO  '
        RETURN
      ENDIF
C
C     Diamagnetic susceptibility with common gauge origin
C
      IF(LABINT(3:8).EQ.'SUSCGO') THEN
        CLASSLB = 'DSUSCGO'
        RETURN
      ENDIF
C
C     Nuclear shielding integrals with common gauge origin
C
      IF(LABINT(3:6).EQ.'NSCO') THEN
        CLASSLB='NSTCGO '
        RETURN
      ENDIF
C
C     Cosine and sine integrals
C
      IF(LABINT(2:7).EQ.'COSINE'.OR.LABINT(2:5).EQ.'SINE') THEN
        CLASSLB = 'EXPIKR '
        RETURN
      ENDIF
C
C     Mass velocity integrals
C
      IF(LABINT.EQ.'MASSVELO') THEN
        CLASSLB = 'MASSVEL'
        RETURN
      ENDIF
C
C     Darwin type integrals
C
      IF(LABINT(1:6).EQ.'DARWIN') THEN
        CLASSLB = 'DARWIN '
        RETURN
      ENDIF
C
C     First order magnetic field derivatives of electric field
C
      IF(LABINT(3:5).EQ.'CM1') THEN
        CLASSLB = 'CM1    '
        RETURN
      ENDIF
C
C     Second order magnetic field derivatives of electric field
C
      IF(LABINT(3:5).EQ.'CM2') THEN
        CLASSLB = 'CM2    '
        RETURN
      ENDIF
C
C     Half-derivative overlap integrals not to be anti-symmetrized
C
      IF(LABINT(1:6).EQ.'SQHDOR') THEN
        CLASSLB = 'SQHDOR '
        RETURN
      ENDIF
C
C
C
      IF(LABINT.EQ.'SQOVLAP') THEN
        CLASSLB = 'SQOVLAP'
        RETURN
      ENDIF
C
C     non-traceless electric field third derivatives
C
      IF(LABINT(5:5).EQ.'v') THEN
        CLASSLB = 'EFTNTL '
        RETURN
      ENDIF
C
C     Electric field gradient at the individual nuclei, cartesian
C
      IF(LABINT(3:5).EQ.'EFG') THEN
        CLASSLB = 'ELFGRDC'
        RETURN
      ENDIF
C
C     Unknown LABINT
C
      CLASSLB = 'Unknown'
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xpldef */
      SUBROUTINE XPLDEF(PLABEL,CLASSLB,IORDER,PCOMB,INDX,
     &                  WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define property label: symmetries etc..
C
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
      CHARACTER PLABEL*8,CLASSLB*7,PCOMB*4,LABINT(15*MXCENT)*8
      DIMENSION WORK(LWORK)
C
      CALL QENTER('XPLDEF')
#include "memint.h"
C
      MAXTYP = 15*MXCENT
      IF ((CLASSLB.EQ.'ELFGRDC').OR.(CLASSLB.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE
         NINTAD = MAXTYP
      END IF
C     Memory allocation
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KATOM ,NUCIND,WORK,KFREE,LFREE)
C
      CALL XPLDE1(PLABEL,CLASSLB,PCOMB,IORDER,WORK(KINTRP),
     &            WORK(KINTAD),LABINT,INDX,WORK(KATOM),
     &            WORK(KFREE),LFREE,IPRINT)
C     Memory deallocation
      CALL MEMREL('XPLDEF.XPLDE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('XPLDEF')      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xplde1 */
      SUBROUTINE XPLDE1(PLABEL,CLASSLB,PCOMB,IORDER,INTREP,INTADR,
     &           LABINT,INDX,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define property label: symmetries etc..
C
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "dcbprl.h"
#include "cbiher.h"
#include "nuclei.h"
      LOGICAL ANTI, TRASPH, SOLVNT, SQUARE,LABERR,
     &        DOATOM(*)
      CHARACTER PLABEL*8,CLASSLB*7,LABINT*8,PCOMB*4
      DIMENSION INTREP(*),INTADR(*),LABINT(*),WORK(LWORK)
C
C
C     Get information on this class of integrals
C
      IF    (CLASSLB.EQ.'MOLFLD ') THEN
        NOPTYP    = 1
        LABINT(1) = 'MOLFIELD'
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSEIF(CLASSLB.EQ.'BETAMAT') THEN
        NOPTYP    = 1
        LABINT(1) = 'BETAMAT '
        INTREP(1) = 0
        SQUARE    = .FALSE.
        ANTI      = .FALSE.
        IORDER    = 0
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
      ELSE
        DO I = 1,NUCIND
          DOATOM(I) = .TRUE.
        ENDDO
        CALL PR1DIR(WORK,LWORK,CLASSLB,INTTYP,NOPTYP,INTREP,
     &            ANTI,SQUARE,INTADR,LABINT,TRIANG,TRASPH,
     &            SOLVNT,IORDER,DOATOM,NATOM,NBAST,NELMNT,IPRINT)
      ENDIF
C
C     Loop over components of operator class
C     ======================================   
C
      LABERR = .TRUE.
      DO ICMP = 1,NOPTYP
        NPRPLBL = NPRPLBL + 1
        IF(NPRPLBL.GT.MAXPRPLBL) THEN
          CALL QUIT('XPLDEF:Redimension MAXPRPLBL !!')
        ENDIF
        IF(PLABEL.EQ.LABINT(ICMP)) THEN
          LABERR = .FALSE.
          INDX   = NPRPLBL
        ENDIF
        PRPLBL (NPRPLBL) = LABINT(ICMP)
        PDOINT (NPRPLBL) = PCOMB
        IPRLREP(NPRPLBL) = INTREP(ICMP)
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
        WRITE(LUPRI,'(A,A8,A)') '*** XPLDEF ERROR: Property label ',
     &        PLABEL,' not found ! *****'
        WRITE(LUPRI,'(A,A7)') 'Class label      : ',CLASSLB
        WRITE(LUPRI,'(A,A8)') '  Property labels: ',LABINT(1)
        WRITE(LUPRI,'(19X,A8)') (LABINT(I),I=2,NOPTYP)
        CALL QUIT('XPLDEF: Label error')
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpsym */
      SUBROUTINE PRPSYM(INDXPR,IPRINT)
C***********************************************************************
C
C     This subroutine finds the symmetry of the unique operator no.
C     INDXPR in COMMON /SCBPRP/.
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
C***********************************************************************
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
      LOGICAL DIFSYM,DIFTIM
      DIMENSION ICREP(3),ICTIM(3)
#include "ibtfun.h"
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
C     Total operator
      ICREP(1) = IBTXOR(IAREP,IBREP)
C
C     Time reversal symmetry
C     ======================
C
C     Matrix operator
      IF(IAIND.EQ.0 .OR. IAIND.EQ.4 .OR.
     &   (IAIND.GE.16 .AND.IAIND.LE.18) ) THEN
C       DIAGONAL, GAMMA5, i*BETA*[XYZ]ALPHA:
        IATIM = 1
      ELSE
C       [XYZ]SIGMA, [XYZ]ALPHA:
        IATIM = -1
      ENDIF
C     Scalar operator
      IBTIM    = IPRLTYP(IBIND)
C     Total operator
      ICTIM(1) = IATIM*IBTIM
      IF(ICTIM(1).EQ.0) THEN
        CALL HEADER('PRPSYM: Abort !!!',-1)
        WRITE(LUPRI,'(A,I5,A)') 
     &    'Operator no.',INDXPR,
     &    ' has indeterminate time reversal symmetry !'
        CALL QUIT('PRPSYM: Indeterminate time reversal symmetry !')
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
        CALL QUIT('PRPSYM: Components of different symmetry !!!')
      ENDIF
      IPRPSYM(INDXPR) = ICREP(1)+1
      IPRPTIM(INDXPR) = ICTIM(1)
C
      IF(IPRINT.GE.10) THEN
        CALL HEADER('Output from PRPSYM',-1)
        CALL WRIXPR(0,INDXPR)
      ENDIF
C
      RETURN
C
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
        WRITE(LUPRI,'(2X,I2,2X,A16,5X,A3,5X,A2)') 
     &   IND,PRPNAM(INDXPR),REP(ISYM-1),TTYP(ITIM)
      ENDIF
      CALL PRSYMB(LUPRI,'.',75,0)
C
C     Operator type
C     =============
C
      ITYP = IPRPTYP(INDXPR)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 ), ITYP
        WRITE(LUPRI,'(6X,A,I5)')
     &  'Program error: unknown operator type:',ITYP
        CALL QUIT('WRIXPR: unknown operator type')
    1 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'DIAGONAL : scalar operator'
      GOTO 20
    2 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XALPHA   : [alpha_x]P'
      GOTO 20
    3 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YALPHA   : [alpha_y]P'
      GOTO 20
    4 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZALPHA   : [alpha_z]P'
      GOTO 20
    5 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XAVECTOR : [alpha x P]_x'
      GOTO 20
    6 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YAVECTOR : [alpha x P]_y'
      GOTO 20
    7 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZAVECTOR : [alpha x P]_z'
      GOTO 20
    8 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ALPHADOT : [alpha.P]'
      GOTO 20
    9 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'GAMMA5   : [gamma_5]P'
      GOTO 20
   10 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XSIGMA   : [Sigma_x]P'
      GOTO 20
   11 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YSIGMA   : [Sigma_y]P'
      GOTO 20
   12 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZSIGMA   : [Sigma_z]P'
      GOTO 20
   13 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XBETASIG : [beta Sigma_x]P'
      GOTO 20
   14 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YBETASIG : [beta Sigma_y]P'
      GOTO 20
   15 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZBETASIG : [beta Sigma_z]P'
      GOTO 20
   16 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'XiBETAAL : i[beta alpha_x]P'
      GOTO 20
   17 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'YiBETAAL : i[beta alpha_y]P'
      GOTO 20
   18 CONTINUE
        WRITE(LUPRI,'(6X,A,A)') 'Operator type ',
     &   'ZiBETAAL : i[beta alpha_z]P'
      GOTO 20
C
   20 CONTINUE
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
          WRITE(LUPRI,'(6X,A,A8,4X,'//FMT//',A12)')
     &      'Labels and factors     : ',
     &      PRPLBL(IPRPLBL(I,INDXPR)),
     &      FACPRP(I,INDXPR),OPTYPE
        ELSE
          WRITE(LUPRI,'(31X,A8,4X,'//FMT//',A12)')
     &      PRPLBL(IPRPLBL(I,INDXPR)),
     &      FACPRP(I,INDXPR),OPTYPE
        ENDIF
      ENDDO
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
                  READ(LUCMD,*) ILL, ILS, ISS
                  IGRD_INTFLG = ILL + 2 * ILS + 4 * ISS
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
