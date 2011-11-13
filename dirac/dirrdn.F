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
C VERSION : $Revision: 1.151.2.10 $
C DATE    : $Date: 2004/09/01 08:08:51 $
C FILE    : dirrdn.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck paminp */
      SUBROUTINE PAMINP(WORK,LWORK)
C***********************************************************************
C     
C     Driver for input module of DIRAC
C     Written by T.Saue April 26 1996
C     Last revision: 260496 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "cbirea.h"
      DIMENSION WORK(*)
C
C     Initializations
C
      CALL PAMINI(WORK,LWORK)
C
C     Specify job
C
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      CALL JOBINP
      CLOSE(LUCMD)
C
C     General input 
C
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      CALL GENINP(WORK,LWORK)
      CLOSE(LUCMD)
C
C     Integral and basis input
C
      OPEN(LUCMD,FILE = 'DIRAC.INP')
      CALL HERINP(WORK,LWORK)
      CLOSE(LUCMD)
C
      DOJACO = DOJACO.AND.(NZ.EQ.1)
C
C     Set DCB common blocks - first time
C     ==================================
C     (will be done second time in PAMSET)
C
      CALL SETDC1(IPREAD)
C
C     Make primitive labels and use label
C     information to list SO-basis
C     ===================================
C
      CALL GETLAB(WORK,LWORK,IPREAD)
C
C     Define Hamiltonian
C     Redefine INTGEN if e.g. Levy-Leblond
C
      OPEN(LUCMD,FILE='DIRAC.INP')
      CALL HAMINP(WORK,LWORK)
      CLOSE(LUCMD)
C
C     Give estimate of memory needed for two-electron integrals
C
      CALL TWOMEM(WORK,LWORK)
C
C
C     Wave function input
C
C     Call PSIINP always to get wave function information for
C     other modules. This could be changed if wave function information
C     is saved on file.  TODO ?
C
      IF ( DOPSI .OR. .NOT. INPTES ) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PSIINP(WORK,LWORK)
         CLOSE(LUCMD)
      END IF
C
C     Input of wave function analysis module
C HJAaJ May 2002:
C   more logical to move CALL PAMANA to PAMPSI, so it follows directly DHF/DFT in output
C   even if e.g. MP2 or RELCCSD is called, because .ANALYZE is currently an analysis
C   of the DHF/DFT wave function, not of the following output.
C   - in a way *ANALYZE ought to be changed to *DHF ANALYZE and
C     .ANALYZE moved to *DHF INPUT ??
C   - one could use the analysis for other wave functions with some changes ?
C
      IF(DOANA) THEN
        OPEN(LUCMD,FILE = 'DIRAC.INP')
        CALL ANWINP('**ANALY',WORK,LWORK)
        CLOSE(LUCMD)
      ENDIF
C
C     Input of property module
C
      IF(DOPRP) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PRPINP('**PROPE',WORK,LWORK)
         CLOSE(LUCMD)
      ENDIF
C
C     Input of transformation module
C
      IF(DOTRA) THEN
        OPEN(LUCMD,FILE = 'DIRAC.INP')
        CALL TRAINP(WORK,LWORK)
        CLOSE(LUCMD)
      ENDIF
C
C     Consistency checks
C
      CALL CONINP
C
      CALL TITLER('End of input processing','*',124)
C
C     If input test , exit here
C     =========================
C
CLV   Replaced stop by (not so elegant) QUIT because we otherwise 
CLV   fail miserably in MPI...
      IF(INPTES) CALL QUIT('Normal exit, only input check was asked.')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck geninp */
      SUBROUTINE GENINP(WORK,LWORK)
C***********************************************************************
C
C <<< General Input for DIRAC >>>
C
C       This routine is based on the corresponding HERINP routine
C       in HERMIT
C
C       Written by Trond Saue, Tromsoe October 1993
C       Last revision: Nov 14 1994 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 8, NTABLE = 24)
C
      LOGICAL LBIT
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7, MXFORM*6, FMT*6
      DIMENSION WORK(LWORK)
      LOGICAL   NOGENINP
C
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "codata.h"
C
      DATA TABDIR /'*END OF','*PARALL','*XXXXXX','*XXXXXX',
     &             '*XXXXXX','*XXXXXX','*XXXXXX','*XXXXXX'/
      DATA TABLE  /'.DIRECT','.CVALUE','.XXXXXX','.DOJACO',
     &             '.PRINT ','.PCMOIN','.PCMOUT','.LINDEP',
     &             '.SPHTRA','.XXXXXX','.XXXXXX','.XXXXXX',
     &             '.ACMOUT','.XXXXXX','.XXXXXX','.FAMILY',
     &             '.DENSFI','.XXXXXX','.XXXXXX','.XXXXXX',
     &             '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
C     Initialise NOGENINP
C     NOGENINP is set, if there is no **GENERAL
C
      NOGENINP = .FALSE.
C
C     Read menu file
C     ==============
C     **** Find General input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 CONTINUE
      READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      IF (WORD .EQ. '**GENER') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         NOGENINP = .TRUE.
         GOTO 180
C         CALL QUIT('End of file on LUCMD, no DIRAC input found')
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no Dirac input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /DCBGEN/
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
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in GENINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in GENINP.')
  101    CONTINUE
C&&&& DIRECT: Direct calculation
                 READ(LUCMD,*) ILLDIR,ISLDIR,ISSDIR
                 IDFLAG = ILLDIR+2*ISLDIR+4*ISSDIR
                 DIRSET = .TRUE.
                 IF(IDFLAG.GT.0) THEN
                   DIRCAL = .TRUE.
                 ELSE
                   DIRCAL = .FALSE.
                 ENDIF
            GO TO 100
  102    CONTINUE
C&&&& CVALUE: Set a different value of light
            READ(LUCMD,*) CVAL
            GO TO 100
  103    CONTINUE
            GO TO 100
  104    CONTINUE
C&&& DOJACO: Do Jacobi diagonalization if real group
C            (to avoid mixing of symmetries)
               DOJACO = .TRUE.
            GO TO 100
  105    CONTINUE
C&&&& PRINT : General print level
            READ(LUCMD,*) IPRGEN
            GO TO 100
  106    CONTINUE
C&&&& PCMOIN: Read coefficients from formatted DFPCMO and
C             write to unformatted DFCOEF
            DOPIN = .TRUE.
            GO TO 100
  107    CONTINUE
C&&&& PCMOUT: Write coefficients to formatted DFPCMO
            DOPUT = .TRUE.
            GO TO 100
  108    CONTINUE
C&&& LINDEP: Thresholds for linear dependence (eigenvalues of overlap matrix)
            READ(LUCMD,*) STOL(1),STOL(2)
            GO TO 100
  109    CONTINUE
C&&&& SPHTRA: Embed transformation to spherical harmonics in
C&&&&         MO-transformation
C     This is transformation to
C      
C     1. Spherical harmonics for Large components
C     2. Scalar restricted kinetic balance type for Small components
C     
            READ(LUCMD,*) ILSPH,ISSPH
            ISPHTR = ILSPH+2*ISSPH
            GO TO 100
  110     CONTINUE
            GO TO 100
  111     CONTINUE
            GO TO 100
  112     CONTINUE
            GO TO 100
  113    CONTINUE
C&&&& ACMOUT: Print coefficients without symmetry (the C1 group).
            DOACUT = .TRUE.    
            GO TO 100
  114    CONTINUE
            GO TO 100
  115    CONTINUE
            GO TO 100
  116    CONTINUE
C&&&& FAMILY: Family basis (Not yet implemented)
            FAMILY = .TRUE.
            GO TO 100
  117    CONTINUE
C&&&& DENSFIT: Density fitting (Not yet implemented)
            DENSFIT = .TRUE.
            GO TO 100
  118    CONTINUE
            GO TO 100
  119    CONTINUE
            GO TO 100
  120    CONTINUE
            GO TO 100
  121    CONTINUE
            GO TO 100
  122    CONTINUE
            GO TO 100
  123    CONTINUE
            GO TO 100
  124    CONTINUE
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in GENINP.')
      END IF
  180 CONTINUE
C
      IF (PARCAL) THEN
ctec we require that it is a direct calculation
         DIRCAL = .TRUE.
         ILLDIR = 1
         ISLDIR = 1
         ISSDIR = 1
         IGTDIR = 0
         IDFLAG = 7
      ENDIF
C
C
C
C     Print section
C     =============
C
      CALL TITLER('General DIRAC set-up','*',126)
C
      WRITE(LUPRI,'(1X,A,F18.7)')
     +   '* The speed of light : ',CVEL
      IF(CVAL.NE.CVEL) THEN
        FMT = MXFORM(CVAL,14)
        WRITE(LUPRI,'(3X,A,1P,D9.2,0P,A,'//FMT//')')
     +   '* The speed of light scaled by ',(CVAL/CVEL),
     +   ' to : ',CVAL
      ENDIF

      IF(PARCAL) WRITE(LUPRI,'(1X,A)') '* Parallel run'
      IF(IDFLAG.GT.0) THEN
        WRITE(LUPRI,'(1X,A)') '* Direct evaluation of '//
     +    'the following two-electron integrals:'
        IF(ILLDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- LL-integrals'
        IF(ISLDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- SL-integrals'
        IF(ISSDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- SS-integrals'
        IF(IGTDIR.EQ.1) WRITE(LUPRI,'(3X,A)') '- GT-integrals'
      ENDIF
C
      IF(DENSFIT) 
     &  WRITE(LUPRI,'(1X,A)') '* Employ density fitting techniques'
C
      IF(LBIT(ISPHTR,1)) THEN
        WRITE(LUPRI,'(1X,A)')
     &   '* Spherical transformation embedded in MO-transformation',
     &   '  for large components'
        IF(ISPHTR.EQ.3) THEN
           WRITE(LUPRI,'(1X,A)')
     &          '* Transformation to scalar RKB basis embedded in',
     &          '  MO-transformation for small components'
        ENDIF
      ELSEIF(ISPHTR.EQ.2) THEN
         WRITE(LUPRI,'(/1X,A)')
     &    'RKB with pre-projection not meaningful with non-spherical'//
     &    ' L component basis !'
         CALL QUIT('GENINP: Pre-RKB with non-spherical L-comp!')
      ENDIF
      WRITE(LUPRI,'(1X,A,2(/3X,A,1P,D9.2))')
     +     '* Thresholds for linear dependence:',
     +     'Large components:  ',STOL(1),
     +     'Small components:  ',STOL(2)
      IF(DOJACO) WRITE(LUPRI,'(1X,A)')
     &    '* Matrix diagonalization by the Jacobi method '//
     &    'for real groups.'
      IF(DOPIN) WRITE(LUPRI,'(1X,A,A)')
     &     '* MO-coefficients read from formatted file ',
     &     ' DFPCMO and written to DFCOEF'
      IF(DOPUT) WRITE(LUPRI,'(1X,A,A)')
     &     '* MO-coefficients written to formatted file DFPCMO'
      WRITE(LUPRI,'(1X,A,I3)')
     &     '* General print level   : ',IPRGEN
C
C     Process input for various program sections
C     ==========================================
C
      IF (NOGENINP) GOTO 300
  200 PROMPT = WORD(1:1)
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
         CALL QUIT('Illegal directory in GENINP.')
      ELSE
         WRITE (LUPRI,'(/,4A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in GENINP, error in prompt.')
      END IF
    2 CONTINUE
#if defined (VAR_MPI)
      CALL PRLINP(WORD)
#endif
        GO TO 200
    3 CONTINUE
        GO TO 200
    4 CONTINUE
        GO TO 200
    5 CONTINUE
        GO TO 200
    6 CONTINUE
        GO TO 200
    7 CONTINUE
        GO TO 200
    8 CONTINUE
        GO TO 200
C
    1 CONTINUE
  300 CONTINUE
#if defined (VAR_MPI)
      CALL PRLINP(WORD)
#endif
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck jobinp */
      SUBROUTINE JOBINP
C***********************************************************************
C
C <<< Job input for DIRAC >>>
C
C     Written by Trond Saue May 1996
C     Last revision: 1997/09/15 - jth (OPTIMIZE added)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (NTABLE = 12, NDIR = 4)
C
      CHARACTER WORD*7, PROMPT*1, TABLE(NTABLE)*7, TABDIR(NDIR)*7,
     &          WORD1*7
C
#include "dcbgen.h"
C
      DATA TABDIR /'*END OF','*OPTIMI','*MINIMI','xXXXXXX'/
C
      DATA TABLE  /'.TITLE ','.INPTES','.4INDEX','.ONLY I',
     &             '.WAVE F','.ANALYZ','.PROPER','.XXXXXX',
     &             '.NO4IND','xXXXXXX','.OPTIMI','.MINIMI'/
C
      NJOBS = 0
C
C     Read menu file
C     ==============
C     **** Find Dirac input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
C
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      IF (WORD .EQ. '**DIRAC') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         CALL QUIT('End of file on LUCMD, no Dirac input found')
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no Dirac input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process job input for COMMON  /CBIPAM/
C     ======================================
C
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in JOBINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in JOBINP.')
  101    CONTINUE
C&&&& TITLE: Title of DHF run
            READ(LUCMD,'(A50)') TITLE
            GO TO 100
  102    CONTINUE
C&&&& INPTES: Input test
            INPTES = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  103    CONTINUE
C&&&& activate 4-INDEX transformation (on by default for CC and CI)
            DOTRA = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  104    CONTINUE
C&&&& ONLY INTEGRALS: Only run integrals
            DOHRM = .TRUE.
            GO TO 100
  105    CONTINUE
C&&&& WAVE FUNCTION: call wavefunction module
            DOPSI = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  106    CONTINUE
C&&&& ANALYZE : call analysis module
            DOANA = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  107    CONTINUE
C&&&& PROPERTIES:  Call property module
            DOPRP = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  108    CONTINUE
            GO TO 100
  109    CONTINUE
C&&&& deactivate 4-INDEX transformation (on by default for CC and CI)
            NOTRA = .TRUE.
            GO TO 100
  110    CONTINUE
            GO TO 100
 111     CONTINUE
C&&&& OPTIMIZE: Geometry optimization
            OPTIMI = .TRUE.
            GOTO 100
 112     CONTINUE
C&&&& MINIMIZE: Geometry optimization 
Ctsaue: retained for backwards compatibility
            OPTIMI = .TRUE.
            GOTO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in JOBINP.')
      END IF
C
C     Print section
C     ==============
C     
  180 CONTINUE
      DO I = 50,1,-1
      IF(TITLE(I:I).NE.' ') GOTO 181
      ENDDO
  181 CONTINUE
      IN = (272-I)/2
      CALL TITLER(TITLE(1:I),'*',IN)
      IF(NJOBS.GT.0) THEN
        WRITE(LUPRI,'(1X,A)') 'Jobs in this run:'
      ENDIF
      IF(DOHRM) THEN
        WRITE(LUPRI,'(3X,A)') '* Only run one-electron integrals'
        GOTO 666
      ENDIF
      IF(DOPSI)   WRITE(LUPRI,'(3X,A)') '* Wave function'
      IF(DOANA)   WRITE(LUPRI,'(3X,A)') '* Analysis'
      IF(DOPRP)   WRITE(LUPRI,'(3X,A)') '* Properties'
      IF(DOTRA) WRITE(LUPRI,'(3X,A)') '* Transformation to Molecu'//
     &    'lar Spinor basis'
      IF(INPTES)  WRITE(LUPRI,'(3X,A)') '* Input test run only '
      IF(OPTIMI)  WRITE(LUPRI,'(3X,A/3X,A)') '* Geometry optimi'//
     &    'zation of DHF-wavefunction','  IMPORTANT'//
     &    ': *WAVE F must be present for readin of wavefunction.'
C
 666  CONTINUE
C
C     Process input for various program sections
C     ==========================================
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/3A/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in JOBINP.')
      ELSE
         WRITE (LUPRI,'(/4A/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in JOBINP, error in prompt.')
      END IF
    2 CONTINUE
C       *OPTIMI
        CALL OPINPU(WORD)
        GO TO 200
    3 CONTINUE
Ctsaue: *MINIMI for backwards compatibility
        CALL OPINPU(WORD)
        GO TO 200
    4 CONTINUE
        GO TO 200
    1 CONTINUE
C     **<something> or *END OF
      CALL OPINPU(WORD)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfinp */
      SUBROUTINE DHFINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for DHF-module
C
C     Written by T.Saue - November 1992
C     Last revision: 1997/09/15 - jth
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 40)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF, LBIT, RESET, NOAUTOCC
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER CTEMP*72, LINE*80, TEXT*20
      DIMENSION WORK(LWORK)
C
#include "nuclei.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbreo.h"
C
      SAVE SET
      DATA TABLE /'.PROJEC','.PRINT ','.MAXITR','.ERGCNV','.EVCCNV',
     &            '.FCKCNV','.NODSCF','.PHCOEF','.FIXDIF','.NODAMP',
     &            '.NODIIS','.MXDIIS','.CNVINT','.ITRINT','.INTFLG',
     &            '.TRIVEC','.TRIFCK','.CLOSED','.EIGPRI','.OPEN S',
     &            '.AUTOCC','.OWNBAS','.NOBNCR','.DIISTH','.DAMPFC',
     &            '.OVLSEL','.NODYNS','.REORDE','.DIISMO','.DIISAO',
     &            '.POST D','.MAX MA','.MAX MI','.NOQCDH','.SKIPEE',
     &            '.SKIPEP','.MAX BA','.PRJTHR','.FROZEN','.XXXXXX'/
      DATA SET/.FALSE./
C
#include "sdpre.h"
#include "ibtfun.h"
#include "memint.h"
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
C *** Initialize in dcbdhf.h ***
C
C 1.  Initialize /DCCDHF/
C     ===================
C
      DHF_INTTYP= '         '
      CACC      = '        '
C
C 2.  Initialize /DCLDHF/
C     ===================
C
      ERGCNV    = .FALSE.
      EVCCNV    = .TRUE.
      FCKCNV    = .FALSE.
      DHFCONV(1)= .FALSE.
      DHFCONV(2)= .FALSE.
      DHFEXIT   = .FALSE.
      TRIVEC    = .TRUE.
      TRIFCK    = .FALSE.
      DODSCF    = .FALSE.
      PHCOEF    = .FALSE.
      FIXDIF    = .FALSE.
      IF(IDFLAG.GT.0) DODSCF    = .TRUE.
      DODAMP    = .TRUE.
      DODIIS    = .TRUE.
      DIISON    = .FALSE.
      OVLSEL    = .FALSE.
      DYNSEL    = .FALSE.
      DIISAO    = .FALSE.
      DIISMO    = .TRUE.
      AUTOCC    = .FALSE.
      BNCRON    = .TRUE.
      INIOCC    = .TRUE.
      PROJEC  = .FALSE.
      PROOWN  = .FALSE.
      NOQCDHF   = .FALSE.
C     ... default do QCDHF if DIIS does not converge
      DHF_SKIPEE = .FALSE.
      DHF_SKIPEP = .FALSE.
C
C 3.  Initialize /DCIDHF/
C     ===================
C
      MAXITR    = 50
      IPRSCF    = IPRGEN
      MXDIIS    = 10
      NBELM     = 1
      ITDIIS    = 0
      ILLINT    = IBTAND(INTGEN,1)
      ISLINT    = IBTAND(INTGEN/2,1)
      ISSINT    = IBTAND(INTGEN/4,1)
      IGTINT    = IBTAND(INTGEN/8,1)
      INTBUF    = 0
      INTFLG    = 0
      ITRINT(1) = 1
      ITRINT(2) = 1
      IPREIG    = 1
      DHF_INTTYP= 'Bare nuc.'
      NEWOCC  = 0
      MXMACRO = 25
      MXMICRO = 25
      MAXBCK  = 5
C
C 4.  Initialize /DCRDHF/
C     ===================
C
      SCFCNV(1) = 1.0D-7
      SCFCNV(2) = 1.0D-6
      ERGVAL    = DUMMY
      EVCVAL    = D0
      FCKVAL    = D0
      DHFERG    = D0
      ELERGY    = D0
      E1PART    = D0
      E2PART    = D0
      ERGBUF    = D0
      CONVRG    = D0
      TDF2      = D0
      TDDG      = D0
      BMCOND    = D0
      DAMPFC    = 0.25D0
      DIISTH    = DUMMY
      CNVINT(1) = DUMMY
      CNVINT(2) = DUMMY
      PRJTHR  = 1.0D-10
      CALL DZERO(SCFTID,8)
      CALL IZERO(ITRSCF,8)
C
C 5.  Initialize /DCODHF/
C     ===================
C
      NELEC_DHF(1) = 0
      NELEC_DHF(2) = 0
      NASH_DHF(1)  = 0
      NASH_DHF(2)  = 0
      CALL IZERO(NACSH,2*MXOPEN)
      NOPEN    = 0
      NFRO_DHF(1)  = 0
      NFRO_DHF(2)  = 0
      CALL DZERO(DA,1+MXOPEN)
      CALL DZERO(DF,1+MXOPEN)
      CALL DZERO(DALPHA,1+MXOPEN)
C
C *** Initialize /CBIREO/ in dcbreo.h
C     ===================
C
      LMOORD    = .FALSE.
      LMOORF    = .FALSE.
C
C     Process input for DHFCAL
C     ========================
C
      NEWDEF = (WORD .EQ. '*DHFCAL')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,
     &                     31,32,33,34,35,36,37,38,39,40), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DHFINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DHFINP.')
    1          CONTINUE
C&&&& PROJECtion: Reduce transformation matrix by projection
               PROJEC = .TRUE.
               READ(LUCMD,*) NPRJREF
               IF(NPRJREF.GT.MPRJREF) THEN
                 WRITE(LUPRI,'(A,I5)')
     &           '* Number of reference fragments specified:',NPRJREF,
     &           '* Current maximum                        :',MPRJREF
                 CALL QUIT('DHFINP:Too many fragments !')
               ENDIF
               DO J = 1,NPRJREF
                 READ(LUCMD,'(A6)') PRJFIL(J)
                 READ(LUCMD,'(A72)') (VCPROJ(I,J),I=1,NFSYM)
               ENDDO
               GO TO 100
    2          CONTINUE
C&&& PRINT : General print level in DHF-module
                  READ(LUCMD,*) IPRSCF
               GO TO 100
    3          CONTINUE
C&&& MAXITR: Maximum number of SCF-iterations
                  READ(LUCMD,*) MAXITR
               GO TO 100
    4          CONTINUE
C&&& ERGCNV: Convergence on energy difference
                  READ(LUCMD,'(A80)') LINE
                  READ(LINE,*,END=41,ERR=41) SCFCNV(1),SCFCNV(2)
                  GOTO 42
   41             READ(LINE,*) SCFCNV(1)
                  SCFCNV(2) = SCFCNV(1)
   42             CONTINUE
                  EVCCNV = .FALSE.
                  ERGCNV = .TRUE.
                  FCKCNV = .FALSE.
               GO TO 100
    5          CONTINUE
C&&& EVCCNV: Convergence on approximate gradient (error vector)
                  READ(LUCMD,'(A80)') LINE
                  READ(LINE,*,END=51,ERR=51) SCFCNV(1),SCFCNV(2)
                  GOTO 52
   51             READ(LINE,*) SCFCNV(1)
                  SCFCNV(2) = SCFCNV(1)
   52             CONTINUE
                  EVCCNV = .TRUE.
                  ERGCNV = .FALSE.
                  FCKCNV = .FALSE.
               GO TO 100
    6          CONTINUE
C&&& FCKCNV: Convergence on total Fock matrix
                  READ(LUCMD,'(A80)') LINE
                  READ(LINE,*,END=61,ERR=61) SCFCNV(1),SCFCNV(2)
                  GOTO 62
   61             READ(LINE,*) SCFCNV(1)
                  SCFCNV(2) = SCFCNV(1)
   62             CONTINUE
                  EVCCNV = .FALSE.
                  ERGCNV = .FALSE.
                  FCKCNV = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& NODSCF: Do not employ differential matrices in SCF-iterations
                  DODSCF = .FALSE.
               GO TO 100
    8          CONTINUE
C&&&& PHCOEF: Phase adjustment of coefficients
                  PHCOEF = .TRUE.
               GO TO 100
    9          CONTINUE
C&&&& FIXDIF: Use fixed value of DIFDEN
                  FIXDIF = .TRUE.
               GO TO 100
   10          CONTINUE
C&&& NODAMP: No damping of Fock matrix
                  DODAMP = .FALSE.
               GO TO 100
   11          CONTINUE
C&&& NODIIS: Do not perform DIIS iterations
                  DODIIS = .FALSE.
               GO TO 100
   12          CONTINUE
C&&& MXDIIS: Maximum size of B-matrix in DIIS
                  READ(LUCMD,*) MXDIIS
                  DODIIS = .TRUE.
               GO TO 100
   13          CONTINUE
C&&& CNVINT: Convergence thresholds for adding SL- and SS-integrals
                  READ(LUCMD,*) CNVINT(1),CNVINT(2)
               GO TO 100
   14          CONTINUE
C&&& ITRINT: Number of iterations before adding SL- and SS-integrals
                  READ(LUCMD,*) ITRINT(1),ITRINT(2)
               GO TO 100
   15          CONTINUE
C&&& INTFLG: Specify what two-itegrals should be included in this run
                  IF (IGTINT.EQ.1) THEN
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT,IGTINT
                  ELSE
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  END IF
               GO TO 100
   16          CONTINUE
C&&& TRIVEC: Read trial vectors from file
                  TRIVEC = .TRUE.
                  TRIFCK = .FALSE.
               GO TO 100
   17          CONTINUE
C&&& TRIFCK: Start on two-electron Fock matrix read from file
                  TRIFCK = .TRUE.
                  TRIVEC = .FALSE.
               GO TO 100
   18          CONTINUE
C&&&& number of closed shell electrons
                  READ (LUCMD,*) (NELEC_DHF(I),I=1,NFSYM)
                  INIOCC = .FALSE.
               GO TO 100
   19          CONTINUE
C&&&& EIGPRI: Print control for electron and positron solutions
                  READ(LUCMD,*) IL,IS
                  IPREIG = IL+2*IS
               GO TO 100
   20          CONTINUE
C&&& ACTIVE SHELLS: read in number of active shells
C                   and their occupation
C             Read in number of open shells
C
               READ(LUCMD,*) NOPEN
               IF (NOPEN.LE.0.OR.NOPEN.GT.MXOPEN) THEN
                  WRITE(LUPRI,'(1X,A,I2,A,I2)') 
     &               'DHFINP: NOPEN must be between ',
     &               1,' and ',MXOPEN
                  CALL QUIT('DHFINP ERROR: Open shell error!')
               END IF
C
C              Read in number of active shells and their occ.
C              Format: 
C              If inversion symmetry: A/K,L
C              If not               : A/K
C              If inversion symmetry and the format is specified
C              as A/K then L is assumed to be zero.
C
               DO IOPEN = 1,NOPEN
                  READ(LUCMD,'(A80)') LINE
C
                  ISLASH = INDEX(LINE,'/')
                  IF (ISLASH .LE. 0) THEN
                     WRITE (LUPRI,'(/A,I4/A/A)')
     &               'ERROR for *DHFCAL .OPEN SHELL shell no.',IOPEN,
     &               'the defining input line does not contain a "/":',
     &               LINE
                     CALL QUIT(
     &               'Input error for .OPEN SHELL under *DHFCAL')
                  END IF
                  TEXT = LINE(1:ISLASH-1)
                  READ (TEXT,'(I20)') IELEC
C
                  IPERIOD= INDEX(LINE,',')
                  IF (IPERIOD .GT. 0) THEN
                     IF (IPERIOD .LE. ISLASH+1) THEN
                        CALL QUIT('ERROR: Misplaced /')
                     END IF
                     TEXT = LINE(ISLASH+1:IPERIOD-1)
                     READ (TEXT,'(I20)') IAC1
                     TEXT = LINE(IPERIOD+1:)
                     READ (TEXT,'(I20)') IAC2
                  ELSE
                     TEXT = LINE(ISLASH+1:)
                     READ (TEXT,'(I20)') IAC1
                     IAC2 = 0
                  END IF
                  IF (NFSYM .EQ. 1) IAC2 = 0
                  IF ( (MOD(IAC1,2) .NE. 0) .OR.
     &                 (MOD(IAC2,2) .NE. 0) .OR.
     &                 (IAC1 .LT. 0) .OR. (IAC2 .LT. 0) .OR.
     &                 (IAC1+IAC2 .EQ. 0)) THEN
                     WRITE(LUPRI,'(A,2I3)')
     &                  'ERROR: Number of spinors must be a '//
     &                  'positive even number:',IAC1,IAC2
                     CALL QUIT('ERROR: Number of spinors must be '//
     &                         'a positive even number')
                  END IF
c                 IF (IELEC .LE. 0 .OR. IELEC .GT. IAC1+IAC2) THEN
                  IF (IELEC .LT. 0 .OR. IELEC .GT. IAC1+IAC2) THEN
                     WRITE(LUPRI,'(A,I3,A,I3)')
     &                  'ERROR: Number of electrons (',IELEC,
     &                  ') must be a positive even number less than',
     &                  IAC1+IAC2
                     CALL QUIT('ERROR: in number of electrons')
                  END IF
                  NACSH(1,IOPEN) = IAC1/2
                  NACSH(2,IOPEN) = IAC2/2
                  DF(IOPEN) = SDPRE(IELEC)/SDPRE(IAC1+IAC2)
               END DO
               INIOCC = .FALSE.
               GO TO 100
   21          CONTINUE
C&&& AUTOCC: Automatic DHF occupation
                  AUTOCC = .TRUE.
               GO TO 100
   22          CONTINUE
C&&&& OWNBASis: Use local basis for projection
               PROOWN = .TRUE.
               GO TO 100
   23          CONTINUE
C&&&& BNCRON : Improved screening for use with bare nucleus start guess
                  BNCRON = .FALSE.
               GO TO 100
   24          CONTINUE
C&&& DIISTH: Change default threshold for DIIS
                    READ(LUCMD,*) DIISTH
               GO TO 100
   25          CONTINUE
C&&& DAMPFC: Change default damping factor
                    READ(LUCMD,*) DAMPFC
               GO TO 100
   26          CONTINUE
C&&& OVLSEL: Activate overlap selection 
C            (default is dynamic overlap selection, see below)
                    OVLSEL = .TRUE.
                    DYNSEL = .TRUE.
               GO TO 100
   27          CONTINUE
C&&& NODYNSEL: No Dynamic overlap selection 
C            (new overlap selection matrix formed in each iteration)
                    DYNSEL = .FALSE.
               GO TO 100
   28          CONTINUE
C&&& LMOORD: Reorder MO orbitals (only on restart on coefficients)
C                 First read the number of orbital reorderings
C                 in each fermion corep.
C
                  LMOORD = .TRUE.
                  IREORD(1) = 0
                  IREORD(2) = 0
                  DO I = 1,NFSYM
                     READ(LUCMD,'(A72)') CTEMP
                     IREORD(I) = -1
                     CALL NUMLST(CTEMP,IMOORD(1,I),MXREORD,1,MXREORD,
     &                           I,IREORD(I),WORK,LWORK)
                  END DO
               GO TO 100
   29          CONTINUE
C&&& DIISMO: Do tsaue/jon DIIS (DIIS in MO-basis, only works for
C                               closed shell)
                  DIISAO = .FALSE.
                  DIISMO = .TRUE.
               GO TO 100
   30          CONTINUE
C&&& DIISAO: Do jth Dalton like DIIS 
                  DIISAO = .TRUE.
                  DIISMO = .FALSE.
               GO TO 100
   31          CONTINUE
C&&& POST DHF REORDER: 
C                 Reorder finished MO orbitals
C                 First read the number of orbital reorderings
C                 in each fermion corep.
C
                  LMOORF = .TRUE.
                  IREORF(1) = 0
                  IREORF(2) = 0
                  DO I = 1,NFSYM
                     READ(LUCMD,'(A72)') CTEMP
                     IREORF(I) = -1
                     CALL NUMLST(CTEMP,IMOORF(1,I),MXREORD,1,MXREORD,
     &                           I,IREORF(I),WORK,LWORK)
                  END DO
               GO TO 100
   32          CONTINUE
C&&& MAX MACRO ITERATIONS: Maximum number of macro iterations for QC-DHF
                  READ(LUCMD,*) MXMACRO
               GO TO 100
   33          CONTINUE
C&&& MAX MICRO ITERATIONS: Maximum number of micro iterations for QC-DHF
                  READ(LUCMD,*) MXMICRO
               GO TO 100
   34          CONTINUE
C&&& NOQCDHF: No QC DHF
                  NOQCDHF = .TRUE.
               GO TO 100
   35          CONTINUE
C&&& SKIPEE: no e-e rotations in 2nd order optimization
                  DHF_SKIPEE = .TRUE.
               GO TO 100
   36          CONTINUE
C&&& SKIPEP: no e-p rotations in 2nd order optimization
                  DHF_SKIPEP = .TRUE.
               GO TO 100
 37            CONTINUE
C&&& MAX BACKSTEPS: max. number of backsteps
                  READ(LUCMD,*) MAXBCK
               GOTO 100
 38            CONTINUE
C&&& PRJTHR: Smallest norm accepted in projection step
                  READ(LUCMD,*) PRJTHR
               GOTO 100
 39               CONTINUE
C&&& FROZEN: Frozen orbitals
               DO I = 1,NFSYM
                 READ(LUCMD,'(A72)') VCFROZ(I)
                 CALL NUMLS3(VCFROZ(I),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),NFRO_DHF(I),NZERO)
                 NFRO_DHF(I)=NFRO_DHF(I)-NZERO
               ENDDO
               GOTO 100
 40               CONTINUE
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in DHFINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in DHFINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Determine occupation
C     ====================
C
C     If there only is one symmetry there is no need to change 
C     occupation during SCF cycles.
C
      IF (NFSYM .EQ. 1) AUTOCC = .FALSE.
C
      IF (INIOCC) THEN
C
         CALL RMOLCHR(ICHRG)
         WRITE(LUPRI,'(/,1X,A,3(/A,I3))')
     &      '* Initial Automatic occupation based on:',
     &      '  Total charge of atoms = ',ICHRG,
     &      '  Charge of molecule    = ',KCHARG,
     &      '  ie. no. of electrons  = ',ICHRG-KCHARG
C
         IF (NFSYM .EQ. 1) THEN
C
C           The easy case: 
C              all electrons goes into the only symmetry!
C
            INIOCC = .FALSE.
            NELEC_DHF(1) = ICHRG - KCHARG
            IF (MOD(NELEC_DHF(1),2).EQ.1) THEN
               NOPEN = 1
               NACSH(1,1) = 1
               DF(1) = DP5
               NELEC_DHF(1) = NELEC_DHF(1) - 1
               WRITE(LUPRI,'(/,6X,A)')
     &            'WARNING! INIOCC assumes 1 electron in 1 open shell!'
            END IF
C
         ELSE
C
C           The difficult case:
C              The occupation is determined in the DHF routines.
C
            DF(0) = D1
            DA(0) = D1
            DALPHA(0) = D0
            NELECT_DHF = ICHRG - KCHARG
            IF (MOD(NELECT_DHF,2).EQ.1)  THEN
               NOPEN = 1
C
               DF(1) = DP5
               DA(1) = D0
               DALPHA(1) = D2
C
               NELECT_DHF = NELECT_DHF - 1
               NAELEC_DHF = 1
               WRITE(LUPRI,'(/,6X,A)')
     &            'WARNING! INIOCC assumes 1 electron in 1 open shell!'
            END IF
C
C           Since guessing the occupation is difficult, we should allow the programme to change
C           occupation during the SCF.
C
            AUTOCC = .TRUE.
         END IF
      END IF
      IF (.NOT. INIOCC) THEN
         DF(0) = D1
         DA(0) = D1
         DALPHA(0) = D0
         DTEMP = D0
         NASHMFT = 0
         DO IFRP = 1,NFSYM
            DO IOPEN = 1,NOPEN
               NASH_DHF(IFRP) = NASH_DHF(IFRP) + NACSH(IFRP,IOPEN)
               NACSHMF(IFRP,IOPEN) = NACSH(IFRP,IOPEN)
               NASHMFT = NASHMFT + NACSHMF(IFRP,IOPEN)
               DTEMP          = DTEMP +
     &                          NACSH(IFRP,IOPEN)*D2*DF(IOPEN)
            END DO
            IF (NELEC_DHF(IFRP).LT.0) 
     &         CALL QUIT('DHFINP: negative number of '//
     &                   'closed shell electrons!')
            IF (MOD(NELEC_DHF(IFRP),2).NE.0) 
     &         CALL QUIT('DHFINP: odd number of '//
     &                   'closed shell electrons!')
            NISH_DHF(IFRP)  = NELEC_DHF(IFRP)/2
         END DO
         NAELEC_DHF = NINT(DTEMP)
         DO IOPEN = 1,NOPEN
C   
C           We do allow fractional occupation of 1
C           This corresponds to a inactive electron
C           and gives DALPHA = 0
C   
            IF (ABS(DF(IOPEN)-D1).LE.DTHRS) THEN
               WRITE(LUPRI,'(1X,A,I2,A,F6.4)')  
     &            'WARNING: Open shell no. ',IOPEN,
     &            ' has frac. occ. of ',DF(IOPEN)
               WRITE(LUPRI,'(1X,2A)')
     &            'This corresponds to an inactive shell. So ALPHA and',
     &            ' A is reset.'
               DA(IOPEN)     = D1
               DF(IOPEN)     = D1
               DALPHA(IOPEN) = D0
            ELSE IF (ABS(DF(IOPEN)) .LE. DTHRS) THEN
               WRITE(LUPRI,'(1X,A,I2,A,F6.4)')  
     &            'WARNING: Open shell no. ',IOPEN,
     &            ' has frac. occ. of ',DF(IOPEN)
               WRITE(LUPRI,'(1X,2A)')
     &            'This corresponds to an secondary shell. So ALPHA',
     &            ' and A is reset.'
               DA(IOPEN)     = D0
               DF(IOPEN)     = D0
               DALPHA(IOPEN) = D0
            ELSE
               DM        = D2*(NACSH(1,IOPEN)+NACSH(2,IOPEN))
               DN        = DM*DF(IOPEN)
               DA(IOPEN) = (DM*(DN-D1))/(DN*(DM-D1)) 
               DALPHA(IOPEN) = (D1-DA(IOPEN))/(D1-DF(IOPEN))
            END IF
         ENDDO
         NASHT_DHF  = NASH_DHF(1)  + NASH_DHF(2)
         NFROT_DHF  = NFRO_DHF(1)  + NFRO_DHF(2)
         NELECT_DHF = NELEC_DHF(1) + NELEC_DHF(2)
      END IF
      IF(PROJEC) THEN
        WRITE(LUPRI,'(A/A,1P,E13.4)')
     &      ' * Variational space reduced by projection'//
     &          ' using a set of reference orbitals.',
     &      '   - Smallest norm accepted:',PRJTHR
        IF(PROOWN) THEN
          WRITE(LUPRI,'(3X,A)') 
     &     '- Reference orbitals are obtained from local bases !',
     &     '  The indexing of reference spinor sets is assumed ',
     &     '  to follow list of symmetry independent nuclei.'
          DO J = 1,NPRJREF
            WRITE(LUPRI,'(3X,A,I3,A,A4)') '- Reference spinor set ',J,
     &             ' --> ',NAMN(J)
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)') 
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)') 
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO J = 1,NPRJREF
            WRITE(LUPRI,'(3X,A,I3)') '- Reference spinor set ',J
            WRITE(LUPRI,'(3X,A,A6)') '  read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(5X,A,A3)') 
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(5X,A,A3,A,A72)') 
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C     If no DIIS set MXDIIS to zero
C     =============================
C
      IF(.NOT.DODIIS) THEN
        MXDIIS = 0
      ENDIF
C
      IF(DOLVC .AND.
     &   .NOT.(ONECAP.AND.(INTV1C.EQ.2.OR.INTV1C.EQ.3)) ) ISSINT = 0
      IF (LEVYLE) THEN
        ISLINT = 0
        ISSINT = 0
        IGTINT = 0
      ENDIF
      IF (ZORA.AND..NOT.ZORA4) THEN
        ISSINT = 0
      ENDIF
      INTDEF = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
      IF(CNVINT(1).LT.DUMMY) ITRINT(1) = 1
      IF(CNVINT(2).LT.DUMMY) ITRINT(2) = 1
C
C     tsaue - 
C     Temporary fix: load HF occupation numbers into dcborb.h
C     
      CALL TDHFORB

C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(DFTHAM) THEN
        WRITE(LUPRI,'(A)') 
     &   ' DHFINP: Set-up for Kohn-Sham calculation:'
      ELSE
        WRITE(LUPRI,'(A)') 
     &   ' DHFINP: Set-up for Hartree-Fock calculation:'
      ENDIF 
      CALL PRSYMB(LUPRI,'=',75,0)
c
c     We disable evaluation of LS, SS, and GT integrals in 
c     One-center app. type 4 /jkp
c
      IF (ONECAP.AND.(INTV1C.EQ.4).AND.INTDEF.NE.1) THEN
         INTDEF = 1
         WRITE(LUPRI,'(A)') 
     &     ' * INTFLG reset to `1 0 0 0` since ONECAP type 4 has '//
     &     'been specified'
      END IF
C
C     Print info about #electrons etc.
C
      WRITE(LUPRI,'(A,I1)') ' * Number of fermion irreps: ',NFSYM
C
      IF (.NOT. INIOCC) CALL PROCC()
      IF(NFROT_DHF.GT.0) THEN
        WRITE(LUPRI,'(A)') ' * Number of frozen shells:'
        WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion irrep ',I,':',NFRO_DHF(I),I=1,NFSYM)
      ENDIF
      IF (NELECT_DHF+NAELEC_DHF .EQ. 0) THEN
         CALL QUIT('DHFINP: Zero electrons not allowed!')
      END IF
C
C
      IF(.NOT.DOSCF) GOTO 999
C
      IF (NOPEN.GT.0.AND.ZORA) THEN
          WRITE(LUPRI,'(A)') 
     &       'DHFINP: Open shell ZORA is not implemented yet, sorry.. '
          CALL QUIT('DHFINP ERROR: asked for ZORA Open shell calc.')
      ENDIF
C
      WRITE(LUPRI,'(A,I3)') ' * General print level   : ',IPRSCF
      WRITE(LUPRI,'(A)') ' ***** TRIAL FUNCTION *****'
      IF   (TRIFCK) THEN
        WRITE(LUPRI,'(A)')
     +     ' * Trial two-electron Fock matrix read from file DF2FCK'
      ELSEIF(TRIVEC) THEN
        WRITE(LUPRI,'(A)')
     +   ' * Trial vectors read from file DFCOEF'
      ELSE
        WRITE(LUPRI,'(A)')
     +   ' * Trial vectors generated by the bare '//
     +   'nucleus approximation'
      ENDIF
      WRITE(LUPRI,'(A)') ' ***** CONVERGENCE CRITERIA *****'
      IF(ERGCNV) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on total energy.',
     &    '   Desired convergence:',SCFCNV(1),
     &    '   Allowed convergence:',SCFCNV(2)
      ELSEIF(FCKCNV) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on total Fock matrix.',
     &    '   Desired convergence:',SCFCNV(1),
     &    '   Allowed convergence:',SCFCNV(2)
      ELSEIF(EVCCNV) THEN
        WRITE(LUPRI,'(A/2(A,1P,D9.3/))')
     &    ' * Convergence on norm of error vector (gradient).',
     &    '   Desired convergence:',SCFCNV(1),
     &    '   Allowed convergence:',SCFCNV(2)
      ENDIF
      WRITE(LUPRI,'(A)') ' ***** CONVERGENCE CONTROL *****'
      IF(DODSCF) THEN
        WRITE(LUPRI,'(A)') ' * Fock matrix constructed using'//
     &       ' differential density matrix'
        IF(FIXDIF) THEN
          WRITE(LUPRI,'(4X,A)')
     &       'with fixed parameter 1.0'
        ELSE
          WRITE(LUPRI,'(4X,A)')
     &       'with optimal parameter.'
        ENDIF
      ENDIF
      IF(DODIIS) THEN
         IF (DIISAO .OR. NASHT_DHF.GT.0) THEN
           WRITE(LUPRI,'(A)') ' * DIIS (in AO basis)'
         ELSEIF(DIISMO) THEN
           WRITE(LUPRI,'(A)') ' * DIIS (in MO basis)'
         ENDIF
         WRITE(LUPRI,'(A,1P,D8.2)')
     &      ' * DIIS will be activated when convergence reaches : ',
     &      DIISTH
         WRITE(LUPRI,'(3X,A,I5)') '- Maximum size of B-matrix:',
     &      MXDIIS
         IF(DODAMP) THEN
           WRITE(LUPRI,'(A/A,F5.3)')
     &     ' * Damping of Fock matrix when DIIS is not activated. ',
     &     '   Weight of old matrix    : ',DAMPFC
         ENDIF
      ELSEIF(DODAMP) THEN
         WRITE(LUPRI,'(A,F5.3)')
     &     ' * Damping of Fock matrix. Weight of old matrix    : ',
     &     DAMPFC
      E NDIF
      WRITE(LUPRI,'(A,I4)')
     +   ' * Maximum number of SCF iterations  : ',MAXITR
C
C     Quadratic convergent DHF?
C
#ifdef MOD_KRMC
      IF (NOQCDHF) THEN
         WRITE(LUPRI,'(A)')
     &   ' * No quadratic convergent Hartree-Fock'
      ELSE
         WRITE(LUPRI,'(A,2(/A,I3))')
     &   ' * Quadratic convergent Hartree-Fock',
     &   '   - Maximum number of macro iterations:',MXMACRO,
     &   '   - Maximum number of micro iterations:',MXMICRO
      END IF
#endif
C
C     Reordering
c
      IF (LMOORD) THEN
         WRITE(LUPRI,'(A)') ' * Reorder start orbitals rbitals...'
         WRITE(LUPRI,'(3X,A)') '- new ordering of start orbitals...'
         DO I = 1,NFSYM
            WRITE(LUPRI,'(4X,A,I1)') 'Fermion corep. ',I
            IF (IREORD(I) .GT. 0) THEN
               WRITE(LUPRI,'(5X,50(I3,A))') 
     &            (IMOORD(J,I),',',J=1,IREORD(I))
            ELSE
               WRITE(LUPRI,'(5X,A)') 'No reordering'
            END IF
         END DO
      END IF
      IF (LMOORF) THEN
         WRITE(LUPRI,'(A)') ' * Reorder converged orbitals...'
         WRITE(LUPRI,'(3X,A)') '- new ordering of converged orbitals...'
         DO I = 1,NFSYM
            WRITE(LUPRI,'(4X,A,I1)') 'Fermion corep. ',I
            IF (IREORF(I) .GT. 0) THEN
               WRITE(LUPRI,'(5X,50(I3,A))') 
     &            (IMOORF(J,I),',',J=1,IREORF(I))
            ELSE
               WRITE(LUPRI,'(5X,A)') 'No reordering'
            END IF
         END DO
      END IF
      IF (DYNSEL .AND. .NOT. OVLSEL) OVLSEL = .TRUE.
      IF (OVLSEL) THEN
         IF (AUTOCC) 
     &      CALL QUIT('WARNING! .AUTOCC combined with overlap '//
     &           'selection does not work')
         IF(.NOT. DYNSEL) THEN
            WRITE(LUPRI,'(A)') ' * Vector selection based on overlap'
         ELSE
            WRITE(LUPRI,'(A/3X,A)') 
     &         ' * Vector selection based on overlap ',
     &            '(dynamic update of overlap-selection matrix)'
         END IF
      END IF
      IF (AUTOCC .AND. NASHT_DHF .GT. 1) THEN
         WRITE(LUPRI,'(A)')
     &      ' WARNING! .AUTOCC and NASHT > 1 does not work',
     &      '  ==> turning .AUTOCC off'
         AUTOCC = .FALSE.
      END IF
      IF (AUTOCC .AND. NFSYM .GT. 1) THEN
         WRITE(LUPRI,'(A)') 
     &      ' * DHF occupation is allowed to change during SCF cycles.'
      END IF
      WRITE(LUPRI,'(A)')
     +    ' * Contributions from 2-electron integrals to Fock matrix:'
      IF(LBIT(INTDEF,1)) THEN
        WRITE(LUPRI,'(3X,A)') 'LL-integrals.'
      ENDIF
      IF(LBIT(INTDEF,2)) THEN
        IF(CNVINT(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    'SL-integrals  below SCF convergence ',CNVINT(1)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'SL-integrals from iteration ',ITRINT(1)
        ENDIF
      ENDIF
      IF(LBIT(INTDEF,3)) THEN
        IF(CNVINT(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,1P,D8.1)')
     +    'SS-integrals  below SCF convergence ',CNVINT(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'SS-integrals from iteration ',ITRINT(2)
        ENDIF
      ENDIF
C
C     Rotations in QC DHF
C
      IF(DHF_SKIPEE) THEN
         WRITE(LUPRI,'(A)')
     &      ' * NB!!! No e-e rotations in 2nd order optimization'
      ENDIF
      IF(DHF_SKIPEP) THEN
         WRITE(LUPRI,'(A)')
     &      ' * NB!!! No e-p rotations in 2nd order optimization.'
      ENDIF
C
      IF ((ONECNV .NE. D0) .AND. (ONECNV .LT. SCFCNV(1))) THEN
         WRITE(LUPRI,'(2(/A,D12.6),/A)')
     &   ' ** ONECAP-WARNING : ONECNV = ',ONECNV,
     &   '                     SCFCNV(1) = ',SCFCNV(1),
     &   '                     ONECNV re-initialized to SCFCNV(1)'
         ONECNV = SCFCNV(1)
      END IF
C
      WRITE(LUPRI,'(A)') ' ***** OUTPUT CONTROL *****'
      IF(IPREIG.EQ.0) WRITE(LUPRI,'(A)')
     +  ' * No eigenvalues written out.'
      IF(IPREIG.EQ.1) WRITE(LUPRI,'(A)')
     +  ' * Only electron eigenvalues written out.'
      IF(IPREIG.EQ.2) WRITE(LUPRI,'(A)')
     +  ' * Only positron eigenvalues written out.'
      IF(IPREIG.EQ.3) WRITE(LUPRI,'(A)')
     +  ' * Both electron and positron component eigenvalues '
     +  //'written out.'

C
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck idnnuc */
      SUBROUTINE IDNNUC(IPRINT)
C***********************************************************************
C
C       Count number of nuclei with identical coordinates
C
C       T.Saue, March 24 1994, University of Odense
C       Last revision : Nov 14 1994 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"
#include "shells.h"
C
      CALL QENTER('IDNNUC')
C
C     Identify nuclei with identical coordinates and give
C     pointer to the first in the series
C     ===================================================
C
      CALL IZERO(INCENT,NUCIND)
      CALL IZERO(INUNIQ,NUCIND)
      NCLOSE = 0
      NUCNET = 0
      DO 20 I = 1,NUCIND
      IF(INCENT(I).NE.0) GOTO 20
        INCENT(I) = I
        NUCNET = NUCNET + 1
        INUNIQ(NUCNET) = I
        DO 30 J = (I+1),NUCIND
        IF(INCENT(J).EQ.J) GOTO 30
          DIST = (CORD(1,I)-CORD(1,J))**2
     +          +(CORD(2,I)-CORD(2,J))**2
     +          +(CORD(3,I)-CORD(3,J))**2
          IF(DIST.LT.0.1) THEN
            INCENT(J) = I
            NCLOSE = NCLOSE + 1
          ENDIF
   30   CONTINUE
   20 CONTINUE
      IF(NUCNET.LT.NUCIND) THEN
        CALL HEADER('Output from IDNNUC',-1)
        WRITE(LUPRI,'(3X,A,I5/9X,A,I5,A)')
     +    '* The number of centers with unique coordinates: ',
     +    NUCNET,
     +    '(out of ',NUCIND,' centers in total)'
        WRITE(LUPRI,'(3X,A)') 'Pointer INCENT:'
        WRITE(LUPRI,'(6X,2I5)') (I,INCENT(I),I=1,NUCIND)
        WRITE(LUPRI,'(3X,A)') 'Pointer INUNIQ:'
        WRITE(LUPRI,'(6X,2I5)') (I,INUNIQ(I),I=1,NUCNET)
      ENDIF
C
      CALL QEXIT('IDNNUC')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fldinp */
      SUBROUTINE FLDINP(WORD,INPERR)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER WORD*7
C
      WRITE(LUPRI,'(A)') 'FLDINP: Sorry, not yet written !'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiinp */
      SUBROUTINE PSIINP(WORK,LWORK)
C***********************************************************************
C
C <<< Wave function Input for DIRAC >>>
C
C
C       Written by Trond Saue, Tromsoe April 1996
C       Last revision: April 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 12, NTABLE = 12)
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7
      DIMENSION WORK(LWORK)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
C
      DATA TABDIR /'*END OF','*DHFCAL','*MP2CAL','*LUCITA',
     &             '*MVOCAL','*RESOLV','*KRMCSC','*STEP C',
     &             '*OPTIMI','XXXXXXX','XXXXXXX','XXXXXXX'/
      DATA TABLE  /'.DHF   ','.MP2   ','.LUCITA','.XXXXXX',
     &             '.XXXXXX','.RELCCS','.DIRRCI','.MVO   ',
     &             '.RESOLV','.XXXXXX','.LUCIAR','.KRMCSC'/
C
      NJOBS = 0
C
C     Initialize /CBLPSI/
C     ===================
C
      DOSCF   = .FALSE.
      DODHF   = .FALSE.
      DODFT   = .FALSE.
      DOMP2   = .FALSE.
      DORES   = .FALSE.
      DOCCM   = .FALSE.
      DOCIM   = .FALSE.
      DOLUCT  = .FALSE.
      DOKRMC  = .FALSE.
      DOLUCIAR= .FALSE.
C
C     Read menu file
C     ==============
C     **** Find Wave function input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      IF (WORD .EQ. '**WAVE ') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         CALL QUIT('End of file on LUCMD, no wave function input found')
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no wave function input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /CBIPSI/
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
     &                109,110,111,112),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in PSIINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in PSIINP.')
C&&&& DHF: Perform Hartree-Fock or Kohn-Sham calculation
  101    CONTINUE
            DOSCF = .TRUE.
            IF(DFTHAM) THEN
              DODFT = .TRUE.
            ELSE
              DODHF = .TRUE.
            ENDIF
            NJOBS = NJOBS +1
            GO TO 100
C&&&& MP2: Perform second-order Moller-Plesset calculation
  102    CONTINUE
            DOMP2 = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
C&&&& LUCITA: Perform spinfree CI calculation
  103    CONTINUE
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOLUCT = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  104    CONTINUE
C&&&& XXXXXX: Empty space for new option (removed LVCOR here)
            GO TO 100
  105    CONTINUE
C&&&& XXXXXX: Empty space for new option (removed TSCOR here)
            GO TO 100
  106    CONTINUE
C&&&& RELCCSD: Perfrom CCSD, CCSD(T) or MP2 with relccsd
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOCCM = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  107    CONTINUE
C&&&& DIRRCI: Perform CISD or COSCI
            IF (.NOT.NOTRA) DOTRA = .TRUE.
            DOCIM = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  108    CONTINUE
C&&&& MVO: Calculate modified virtual orbitals
            DOMVO = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  109    CONTINUE
C&&&& RESOLVE: Resolve open-shell states
            DORES = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  110    CONTINUE
C&&&& XXXXXX: Empty space for new option (removed surplus LUCITA entry!)
            GO TO 100
  111    CONTINUE
C&&&& LUCIAREL: Call LUCIAREL
            DOLUCIAR = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  112    CONTINUE
C&&&& KRMCSCF: KR-MCSCF
            DOKRMC = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in PSIINP.')
      END IF
  180 CONTINUE
C
C     Print section
C     =============
C
      CALL TITLER('Wave function module','*',126)
      IF(NJOBS.GT.0) THEN
        WRITE(LUPRI,'(A)') ' Jobs in this run (in execution order):'
      ENDIF
      IF(DODHF) WRITE(LUPRI,'(A)') ' * Hartree-Fock calculation'
      IF(DODFT) WRITE(LUPRI,'(A)') ' * Density functional '//
     &    'calculation (Kohn-Sham method)'
#ifndef MOD_DFT
      IF(DODFT) THEN
        CALL QUIT('*** DFT is not included in this version')        
      ENDIF
#endif
      IF(DORES) WRITE(LUPRI,'(A)') 
     &   ' * Followed by resolution of open-shell states'
      IF(DOMP2) WRITE(LUPRI,'(A)') ' * MP2 calculation'
      IF(DOMVO) WRITE(LUPRI,'(A)') 
     &   ' * Generate modified virtual orbitals'
      IF(DOKRMC) WRITE(LUPRI,'(A)')
     &   ' * Kramers restricted MCSCF calculation'
#ifndef MOD_KRMC
      IF(DOKRMC .OR. DOLUCIAR)
     &  CALL QUIT('*** KRMCSCF is not included in this version')        
#endif
      IF(DOLUCIAR) THEN
C hj aug03: TODO: idea was/is to maybe follow MCSCF by big CI 
         WRITE(LUPRI,'(A)') ' * Run LUCIAREL CI code'
         WRITE(LUPRI,'(A)') '   ERROR: not implemented yet'
         CALL QUIT('.LUCIAREL as stand-alone not implemented yet')
      END IF
      IF(DOCCM) WRITE(LUPRI,'(A)') ' * Run RELCCSD code'
      IF(DOCIM) WRITE(LUPRI,'(A)') ' * Run DIRRCI CI code'
      IF(DOLUCT)WRITE(LUPRI,'(A)') ' * Run LUCITA CI code'
C
C     Process input for various program sections
C     ==========================================
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in PSIINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in PSIINP, error in prompt.')
      END IF
    2 CONTINUE
        CALL DHFINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    3 CONTINUE
        CALL MP2INP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    4 CONTINUE
        CALL LUCTINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    5 CONTINUE
        CALL MVOINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    6 CONTINUE
        CALL RESINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    7 CONTINUE
#ifdef MOD_KRMC
        CALL KRMCINP(WORD,.FALSE.,WORK,LWORK)
#endif
        GO TO 200
    8 CONTINUE
#ifdef MOD_KRMC
        CALL RSTPINP(WORD,.FALSE.,WORK,LWORK)
#endif
        GO TO 200
    9 CONTINUE
#ifdef MOD_KRMC
        CALL ROPTINP(WORD,.FALSE.,WORK,LWORK)
#endif
        GO TO 200
C
    1 CONTINUE
C
C     For geometry optimizations to work the INPUT routines
C     must always be called.
C
C     IF (DODHF) CALL DHFINP(WORD,.TRUE.,WORK,LWORK)
C     IF (DOMVO) CALL MVOINP(WORD,.TRUE.,WORK,LWORK)
C     IF (DORES) CALL RESINP(WORD,.TRUE.,WORK,LWORK)
C     IF (DOMP2) CALL MP2INP(WORD,.TRUE.,WORK,LWORK)
      CALL DHFINP(WORD,.TRUE.,WORK,LWORK)
      CALL MVOINP(WORD,.TRUE.,WORK,LWORK)
      CALL RESINP(WORD,.TRUE.,WORK,LWORK)
      CALL MP2INP(WORD,.TRUE.,WORK,LWORK)
      CALL LUCTINP(WORD,.TRUE.,WORK,LWORK)
#ifdef MOD_KRMC
      CALL KRMCINP(WORD,.TRUE.,WORK,LWORK)
      CALL RSTPINP(WORD,.TRUE.,WORK,LWORK)
      CALL ROPTINP(WORD,.TRUE.,WORK,LWORK)
#endif
C
C     Input of relccsd module
C
      IF(DOCCM) THEN
C        DIRAC.INP is already open.
C        OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL CCINPT(LUCMD,LUPRI)
         CLOSE(LUCMD)
      ENDIF
C
C     Input of dirrci module
C
      IF(DOCIM) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL CIMINP(WORK,LWORK)
         CLOSE(LUCMD)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCTINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C       Input section for LUCITA module
C
C       Written by Timo Fleig, December 2001
C       (modified DHFINP routine)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "luctinp.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
      LOGICAL   SET,RESET,NEWDEF
      CHARACTER WORD*7,WORD1*7,PROMPT*1,TABLE(NTABLE)*7
C
      DIMENSION WORK(LWORK)
      DIMENSION IMOKW(NTABLE)
C
C
      DATA TABLE /'.TITLE ','.INIWFC','.CITYPE','.NROOTS','.SYMMET',
     &            '.NACTEL','.MULTIP','.PRINTG','.PRINTL','.SZCALC',
     &            '.INACTI','.GASSHE','.GASSPC','.FROORB','.RAS1  ',
     &            '.RAS2  ','.RAS3  ','.DENSI ','.RSTRCI','.XXXXXX'/
      DATA SET/.FALSE./
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
C     Initialize keyword status flags driving interface.
C     (Therefore local keywords do not have to be initialized!)
C
      DO IKW = 1,NTABLE
         IMOKW(IKW) = 0
      END DO
C
C
C     Process input for LUCITA
C     ========================
C
      NEWDEF = (WORD .EQ. '*LUCITA')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in LUCITA.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LUCTINP.')
C
C&&& TITLE : One title line for CI specification
    1          CONTINUE
                  READ(LUCMD,'(A)') TITLUC
                  IMOKW(1) = 1
               GO TO 100
    2          CONTINUE
C&&& INIWFC: Initial DHF wave function (closed- or open-shell)
                  READ(LUCMD,'(A)') WAFFCD
                  IMOKW(2) = 1
               GO TO 100
    3          CONTINUE
C&&& CITYPE: Type of CI calculation (FCI, SD, SDTQ, RAS, GAS)
                  READ(LUCMD,'(A)') CALCTP
                  IMOKW(3) = 1
               GO TO 100
    4          CONTINUE
C&&& NROOTS: Number of eigenvalues
                  READ(LUCMD,*) NROOTD
                  IMOKW(4) = 1
               GO TO 100
    5          CONTINUE
C&&& SYMMET: State symmetry
                  READ(LUCMD,*) ISSYMD
                  IMOKW(5) = 1
               GO TO 100
    6          CONTINUE
C&&& NACTEL: Number of active electrons
                  READ(LUCMD,*) NACTED
                  IMOKW(6) = 1
               GO TO 100
    7          CONTINUE
C&&& MULTIP: State spin multiplicity
                  READ(LUCMD,*) IMULTD
                  IMOKW(7) = 1
               GO TO 100
    8          CONTINUE
C&&& PRINTG: Global LUCIA print flag
                  READ(LUCMD,*) IPRNGD
                  IMOKW(8) = 1
               GO TO 100
    9          CONTINUE
C&&& PRINTL: Local print flag DIRLUC
                  READ(LUCMD,*) IPRNLD
                  IMOKW(9) = 1
               GO TO 100
   10          CONTINUE
C&&& SZCALC: Rough size of calculation
                  READ(LUCMD,'(A)') SZCALD
                  IMOKW(10) = 1
               GO TO 100
   11          CONTINUE
C&&& INACTI: Inactive function space
                  READ(LUCMD,'(A)') CRDINA
                  IMOKW(11) = 1
               GO TO 100
   12          CONTINUE
C&&& GASSHE: Number of GAS and distribution
                  READ(LUCMD,*) INGASD
                  DO IR = 1,INGASD,1
                    READ(LUCMD,'(A)') CRDGAS(IR)
                  END DO
                  IMOKW(12) = 1
               GO TO 100
   13          CONTINUE
C&&& GASSPC: Cumulated max. and min. GAS occupations
                  READ(LUCMD,*) NSEQCD
                  DO IR = 1,INGASD,1
                    READ(LUCMD,'(A)') CRDGOC(IR)
                  END DO
                  IMOKW(13) = 1
               GO TO 100
   14          CONTINUE
C&&& FROORB: Frozen orbital function space
                  READ(LUCMD,'(A)') CRDFRO
                  write(6,'(6A,1A72)') 'CRDFRO',CRDFRO
                  IMOKW(14) = 1
               GO TO 100
   15          CONTINUE
C&&& RAS1  : RAS1 function space
                  READ(LUCMD,'(A)') CRDRS1
                  READ(LUCMD,*) MXHL1D
                  IMOKW(15) = 1
               GO TO 100
   16          CONTINUE
C&&& RAS2  : RAS2 function space
                  READ(LUCMD,'(A)') CRDRS2
                  IMOKW(16) = 1
               GO TO 100
   17          CONTINUE
C&&& RAS3  : RAS3 function space
                  READ(LUCMD,'(A)') CRDRS3
                  READ(LUCMD,*) MXEL3D
                  IMOKW(17) = 1
               GO TO 100
   18          CONTINUE
C&&& DENSI : Level of density matrix calculation (1-,1- and 2-particle)
                  READ(LUCMD,*) IDENSD
                  IMOKW(18) = 1
               GO TO 100
   19          CONTINUE
C&&& RSTRCI: Restart option from CI vector
                  READ(LUCMD,*) IRSTLT
                  IMOKW(19) = 1
               GO TO 100
   20          CONTINUE
C&&& XXXXXX: Free space for additional flag
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in DHFINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in DHFINP.')
            END IF
      END IF
  300 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ANWINP */
      SUBROUTINE ANWINP(WRDSRC,WORK,LWORK)
C***********************************************************************
C
C <<< Analysis Input for DIRAC >>>
C
C
C       Written by Trond Saue, April 1996
C       Last revision: April 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 8, NTABLE = 8)
#include "maxorb.h"
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7,WRDSRC*7
      DIMENSION WORK(LWORK)
C
#include "dcblab.h"
#include "dcbana.h"
C
      DATA TABDIR /'*END OF','*MULPOP','*PROJEC','*PRIVEC',
     &             '*RHO1  ','*DENSIT','*WT1   ','*XXXXXX'/
      DATA TABLE  /'.MULPOP','.PROJEC','.PRIVEC','.RHO1  ',
     &             '.DENSIT','.WT1   ','.XXXXXX','.XXXXXX'/
C
      CALL QENTER('ANWINP')
#include "memint.h"
C
      NJOBS = 0
C
C     Initialize /CBLPRP/
C     ===================
C
      DOVEC   = .FALSE.
      DOPOP   = .FALSE.
      DOPRJ   = .FALSE.
      DO1RHO  = .FALSE.
      DO3RHO  = .FALSE.
      DO1WT   = .FALSE.
C
C     Read menu file
C     ==============
C
C
C     **** Find Analysis input *****
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
         WRITE(LUPRI,'(A)') 'End of file on LUCMD,'//
     &                      ' no analysis input found'
         CALL QEXIT('ANWINP')
         RETURN
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no analysis input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /CBIANA/
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
               GO TO (101,102,103,104,105,106,107,108),I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in ANWINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in ANWINP.')
C&&&& MULPOP: Perform Mulliken population analysis
  101    CONTINUE
            DOPOP = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
C&&&& PROJECTION: Perform projection of current solution onto another solution
  102    CONTINUE
            DOPRJ = .TRUE.
            NJOBS = NJOBS +1
            GO TO 100
  103    CONTINUE
C&&&& PRIVEC: Print vectors
            DOVEC = .TRUE.
            NJOBS = NJOBS +1
            GO TO 100
  104    CONTINUE
C&&&& RHO1: Punch density along the lines between centers
            DO1RHO = .TRUE.
            NJOBS = NJOBS +1
            GO TO 100
  105    CONTINUE
C&&&& DENSITY: Write density to unformatted file (Gaussian cube format)
            DO3RHO = .TRUE.
            GO TO 100
  106    CONTINUE
C&&&& RHO1: Punch Becke integration weights along the lines between atomic centers
            DO1WT = .TRUE.
            NJOBS = NJOBS +1
            GO TO 100
  107    CONTINUE
            GO TO 100
  108    CONTINUE
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 180
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in ANWINP.')
      END IF
  180 CONTINUE
C
C     Print section
C     =============
C
      CALL TITLER('Analysis module','*',129)
      IF(NJOBS.GT.0) THEN
        WRITE(LUPRI,'(1X,A)') 'Jobs in this run:'
      ENDIF
      IF(DOVEC) WRITE(LUPRI,'(1X,A)') '* Write vectors'
      IF(DOPOP) WRITE(LUPRI,'(1X,A)') 
     &     '* Mulliken population analysis'
      IF(DOPRJ) WRITE(LUPRI,'(1X,A)') 
     &     '* Projection onto another solution '
      IF(DO3RHO) WRITE(LUPRI,'(1X,A)') 
     &    '* Write out density for visualization (Gaussian cube format)'
      IF(DO1RHO) WRITE(LUPRI,'(1X,A)') 
     &     '* Punch density along lines between nuclear centers'
      IF(DO1WT) WRITE(LUPRI,'(1X,A)')
     &     '* Punch Becke weights along lines between nuclear centers'
C
C
C     Process input for various program sections
C     ==========================================
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
         CALL QUIT('Illegal directory in ANWINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in ANWINP, error in prompt.')
      END IF
    2 CONTINUE
        CALL POPINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    3 CONTINUE
        CALL PRJINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    4 CONTINUE
        CALL VECINP(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    5 CONTINUE
        CALL RHO1IN(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    6 CONTINUE
        CALL DENSIN(WORD,.FALSE.,WORK,LWORK)
        GO TO 200
    7 CONTINUE
#if defined (MOD_DFT)
        CALL WT1IN(WORD,.FALSE.,WORK,LWORK)
#endif
        GO TO 200
    8 CONTINUE
        GO TO 200
C
    1 CONTINUE
      IF(DOPOP) CALL POPINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOPRJ) CALL PRJINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOVEC) CALL VECINP(WORD,.TRUE.,WORK,LWORK)
      IF(DO1RHO) CALL RHO1IN(WORD,.TRUE.,WORK,LWORK)
      IF(DO3RHO) CALL DENSIN(WORD,.TRUE.,WORK,LWORK)
#if defined (MOD_DFT)
      IF(DO1WT) CALL WT1IN(WORD,.TRUE.,WORK,LWORK)
#endif
C
      CALL QEXIT('ANWINP')
      RETURN
C
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popinp */
      SUBROUTINE POPINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for Mulliken population analysis
C
C     Written by T.Saue - April 1996
C     Last revision: April 27 - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbpop.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      LOGICAL SET, NEWDEF,ORBDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,LINE*72
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.AOLAB ','.LABDEF','.VECPOP','.XXXXXX',
     &            '.PRINT ','.INDSML','.NETPOP','.XXXXXX'/
      DATA SET/.FALSE./
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
#include "memint.h"
C
C
C     Initialize /CBIPOP/
C     ===================
C
      IPRPOP    = 0
      DONETP  = .FALSE.
      ADDSML  = .TRUE.
      ORBDEF  = .TRUE.
      ILABDF  = 2
C
C     Process input from CBIPOP
C     =========================
C
      NEWDEF = (WORD .EQ. '*MULPOP')
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
     *            '" not recognized in ANWINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in POPINP.')
    1          CONTINUE
C&&&& AOLAB:  Base definition of labels on AO-labels:
                 ILABDF = 1
               GO TO 100
    2          CONTINUE
C&&&& LABDEF: Define groups of primitive labels
                 LABDEF = .TRUE.
                 READ(LUCMD,*) NPOPLAB
                 CALL MEMGET('INTE',KIBUF,MAXLAB,WORK,KFREE,LFREE)
                 CALL PROLAB(IPOPLAB,POPLAB,MAXLAB,NPOPLAB,LINE,LUCMD,
     &                       WORK(KIBUF),WORK(KFREE),LFREE)
                 CALL MEMREL('LABDEF',WORK,KWORK,KWORK,KFREE,LFREE)
               GO TO 100
    3          CONTINUE
C&&&& VECPOP: Number of spinors to analyze
                  READ(LUCMD,'(A72)') (VECPOP(I),I=1,NFSYM)
                  ORBDEF = .FALSE.
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE
C&&&& IPRPOP: Print level
                  READ(LUCMD,*) IPRPOP
               GO TO 100
    6          CONTINUE
C&&&& INDSML: Do not collect labels of small components of given center/symmetry
                  ADDSML = .FALSE.
               GO TO 100
    7          CONTINUE
C&&&& NETPOP: Mulliken Net populations
                 DONETP = .TRUE.
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ANWINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ANWINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      IF(ORBDEF) THEN
        DO I = 1,NFSYM
          VECPOP(I) = ' '
          IF(NOCC(I).GT.0) THEN
            IF(NOCC(I).EQ.1) THEN
              WRITE(VECPOP(I),'(I4)') 1
            ELSE
              WRITE(VECPOP(I),'(I4,A2,I4)') 1,'..',NOCC(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(.NOT.DOPOP) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'POPINP: Mulliken population analysis'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') '* Gross populations'
      IF(DONETP) WRITE(LUPRI,'(1X,A)') 
     &     '* Net populations and overlap'
      IF(ILABDF.EQ.1) THEN
        WRITE(LUPRI,'(1X,A)') '* Label definitions based on AO-labels'
      ELSE
        WRITE(LUPRI,'(1X,A)') '* Label definitions based on SO-labels'
      ENDIF
      IF(.NOT.ADDSML) WRITE(LUPRI,'(3X,A)')
     &      '(individual small component labels will be used)'
      WRITE(LUPRI,'(1X,A)') '* Number of spinors analyzed:'
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPOP(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)') 
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPOP(I)
        ENDIF
      ENDDO
      WRITE(LUPRI,'(1X,A,I5)') '* Print level:',IPRPOP
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjinp */
      SUBROUTINE PRJINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for module for wavefunction projection
C
C     Written by T.Saue - April 1996
C     Last revision: Apr 26 - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbprj.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "nuclei.h"
C
#include "dcbxpr.h"
      LOGICAL SET, NEWDEF,ORBDEF,RESET,SYMSTP
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.PRINT ','.PROTHR','.VECPRJ','.VECREF',
     &            '.WGPOP','.OWNBAS','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
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
C     Initialize /CBIPRJ/
C     ===================
C
      IPRPRJ   = 0
C     Default is to project all occupied vectors onto all reference vectors
      ORBDEF = .TRUE.
      NREFS    = 0
      PROTHR   = 1.0D-3
      LWGPOP = .FALSE.
      OWNBAS = .FALSE.
C
C     Process input from CBIPRJ
C     =========================
C
      NEWDEF = (WORD .EQ. '*PROJEC')
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
     *            '" not recognized in PRJINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in PRJINP.')
    1          CONTINUE
C&&&& PRINT: Print level for projections
                 READ(LUCMD,*) IPRPRJ
               GO TO 100
    2          CONTINUE
C&&&& PROTHR: Threshold for printout of projection coefficients
                 READ(LUCMD,*) PROTHR
               GO TO 100
    3          CONTINUE
C&&&& VECPRJ: Spinors to project onto reference vectors
                  READ(LUCMD,'(A72)') (VECPRJ(I),I=1,NFSYM)
                  ORBDEF = .FALSE.
               GO TO 100
    4          CONTINUE
C&&&& VECREF: Reference vectors
                  READ(LUCMD,*) NREFS
                  IF((NREFS+1).GT.MAXREF) THEN
                    WRITE(LUPRI,'(A,I5)')
     &         '* Number of reference fragments specified:',NREFS,
     &         '* Current maximum                        :',MAXREF
                    CALL QUIT('PRJINP:Too many fragments !')
                  ENDIF
                  DO J = 1,NREFS
                    READ(LUCMD,'(A6)') REFFIL(J)
                    READ(LUCMD,'(A72)') (VECREF(I,J),I=1,NFSYM)
                  ENDDO
                  REFFIL(NREFS+1) = 'Polari'
               GO TO 100
    5          CONTINUE
C&&& WGPOP: Split overlap densities according to weight of contributions
                 LWGPOP = .TRUE.
               GO TO 100
    6          CONTINUE
C&&& OWNBAS: Fragments (defined from NUCIND)  are calculated in their own basis sets
                 OWNBAS = .TRUE.
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in PRJINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in PRJINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      IF(ORBDEF) THEN
        DO I = 1,NFSYM
          VECPRJ(I) = ' '
          IF(NOCC(I).GT.0) THEN
            IF(NOCC(I).EQ.1) THEN
              WRITE(VECPRJ(I),'(I4)') 1
            ELSE
              WRITE(VECPRJ(I),'(I4,A2,I4)') 1,'..',NOCC(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(.NOT.DOPRJ) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'PRJINP: Projection onto another solution'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,1P,D9.2)')
     &       '* Threshold for printing projection coefficients:',PROTHR
      WRITE(LUPRI,'(1X,A)') '* Number of spinors projected:'
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPRJ(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)') 
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPRJ(I)
        ENDIF
      ENDDO
      IF(OWNBAS) THEN
        WRITE(LUPRI,'(1X,A)') 
     &   '* Fragments are calculated in their individual bases !',
     &   '    The indexing of reference spinor sets is assumed ',
     &   '    to follow list of symmetry independent nuclei.'
        DO J = 1,NREFS
          WRITE(LUPRI,'(1X,A,I3,A,A6)') '* Reference spinor set ',J,
     &           ' --> ',NAMDEP(J)
          WRITE(LUPRI,'(1X,A,A6)') ' - read from file ',REFFIL(J)
          DO I = 1,NFSYM
            NVEC = 0
            CALL  NUMLST(VECREF(I,J),IDUMMY,NFBAS(I,0),
     &                -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
            IF(NVEC.EQ.0) THEN
              WRITE(LUPRI,'(4X,A,A3)') 
     &           '- No orbitals in fermion ircop ',FREP(I)
            ELSE
              WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &        '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I,J)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO J = 1,NREFS
          WRITE(LUPRI,'(1X,A,I3)') '* Reference spinor set ',J
          WRITE(LUPRI,'(1X,A,A6)') ' - read from file ',REFFIL(J)
          DO I = 1,NFSYM
            NVEC = 0
            CALL  NUMLST(VECREF(I,J),IDUMMY,NFBAS(I,0),
     &                -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
            IF(NVEC.EQ.0) THEN
              WRITE(LUPRI,'(4X,A,A3)') 
     &           '- No orbitals in fermion ircop ',FREP(I)
            ELSE
              WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &        '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I,J)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF(LWGPOP) WRITE(LUPRI,'(A,A)') '* Overlap densities split ',
     &  'between contributions according to their weight.'
      WRITE(LUPRI,'(1X,A,I5)') '* Print level:',IPRPRJ
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vecinp */
      SUBROUTINE VECINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for vector print
C
C     Written by T.Saue - May 1996
C     Last revision: May 2 - 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbvec.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRICMP','.AOLAB ','.VECPRI','.PRINT ',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
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
C     Initialize /CBIVEC/
C     ===================
C
      IPRCMP    = 1
      DO I = 1,NFSYM
        VECPRI(I) = ' '
        IF(NOCC(I).EQ.0) THEN
          WRITE(VECPRI(I),'(I1)') 1
        ELSEIF(NOCC(I).GT.0) THEN
          WRITE(VECPRI(I),'(I1,A2,I4)') 1,'..',NOCC(I)
        ENDIF
      ENDDO
      ILABDF  = 2
      IPRVEC  = 0
C
C     Process input from CBIPOP
C     =========================
C
      NEWDEF = (WORD .EQ. '*PRIVEC')
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
     *            '" not recognized in VECINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in VECINP.')
    1          CONTINUE
C&&&& PRICMP: Print control for large and small component vectors
                  READ(LUCMD,*) IL,IS
                  IPRCMP = IL+2*IS
               GO TO 100
    2          CONTINUE
C&&&& AOLAB : Print coefficients in AO-basis
                  ILABDF = 1
               GO TO 100
    3          CONTINUE
C&&&& VECPRI: Number of spinors to print
                  READ(LUCMD,'(A72)') (VECPRI(I),I=1,NFSYM)
               GO TO 100
    4          CONTINUE
C&&&& PRINT: General print level
                  READ(LUCMD,*) IPRVEC
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
     *            '" not recognized in VECINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in VECINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DOVEC) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'VECINP: Vector print'
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(ILABDF.EQ.1) THEN
        WRITE(LUPRI,'(1X,A)') '* Coefficients written in AO-basis'
      ELSE
        WRITE(LUPRI,'(1X,A)') '* Coefficients written in SO-basis'
      ENDIF
      WRITE(LUPRI,'(1X,A)') '* Vector print:'
      NTOT = 0
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPRI(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)') 
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          NTOT = NTOT + NVEC
          WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPRI(I)
        ENDIF
      ENDDO
      IF(NTOT.EQ.0) IPRCMP = 0
      IF(IPRCMP.EQ.1) WRITE(LUPRI,'(1X,A)')
     +  '* Only large component coefficients written out.'
      IF(IPRCMP.EQ.2) WRITE(LUPRI,'(1X,A)')
     +  '* Only small component coefficients written out.'
      IF(IPRCMP.EQ.3) WRITE(LUPRI,'(1X,A)')
     +  '* Both large and small component coefficients written out.'
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamini */
      SUBROUTINE PAMINI(WORK,LWORK)
C***********************************************************************
C     
C     Initialize data in DIRAC
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
#include "maxaqn.h"
#include "codata.h"
#include "maxorb.h"
C
#include "siripc.h"
#include "huckel.h"
#include "dcbgen.h"
#include "ccom.h"
#include "cbirea.h"
#include "dcbxpr.h"
#include "dcbcls.h"
#include "dcbprl.h"
      DIMENSION WORK(LWORK)
C
C     Initialize temporary variables
C     ==============================
C
      ILSPH  = 1
      ISSPH  = 1
      ILLDIR = 1
      ISLDIR = 1
      ISSDIR = 1
      IGTDIR = 0
      ILLINT = 1
      ISLINT = 1
      ISSINT = 1
      IGTINT = 0
C
C     Initialize /DCBGEN/
C     ===================
C
      TITLE    = 'DIRAC: No title specified !!!'
      FAMILY   = .FALSE.
      INPTES   = .FALSE.
      VACUUM   = .FALSE.
      NUCVAC   = .FALSE.
      DOPSI    = .FALSE.
      DOANA    = .FALSE.
      DOPRP    = .FALSE.
      DOHRM    = .FALSE.
      DOJACO   = .FALSE.
      DOPIN    = .FALSE.
      DOPUT    = .FALSE.
      DOACUT   = .FALSE.
      DOTRA    = .FALSE.
      NOTRA    = .FALSE.
      OPTIMI   = .FALSE.
      DENSFIT  = .FALSE.
      IPRGEN   = 0
      ISPHTR   = ILSPH+2*ISSPH
      IDFLAG   = ILLDIR+2*ISLDIR+4*ISSDIR+4*IGTDIR
      INTGEN   = ILLINT+2*ISLINT+4*ISSINT+4*IGTINT
      CVAL     = CVEL
      STOL(1)  = 1.0D-6
      STOL(2)  = 1.0D-8
      PANAS  = D0
      HRINPC = .FALSE.
      RDINPC = .FALSE.
      IF(IDFLAG.GT.0) THEN
        DIRCAL = .TRUE.
      ELSE
        DIRCAL = .FALSE.
      ENDIF
C --->
C
C     Test : disable automatic conversional evaluation of
C            one-center LS and SS integrals in ONECAP models. 
C            These can also take up quite some space !!!   /jkp
C
      DIRSET = .TRUE.
C <---
      IPRUSR = 0
      USRIPR = .FALSE.
      SEGBAS = .TRUE.
      NEWBAS = .TRUE.
      NEWPRP = .TRUE.
      RELCAL = .TRUE.
C
C     Initialize SIRIPC from DALTON
C
      NEWGEO = .TRUE.
      INPPRC = .FALSE.
C
C     Initialize /HUCKEL/from DALTON
C
      ADDSTO = .FALSE.
      HUCCNT = 2.0D0
      CALL IZERO(NHUCAO,8)
      CALL IZERO(IHUCPT,MXSHEL)
C
C     Initialize XCBCLS - list of integrals to calculate
C     ==================================================
      NPRPCLS   = 0
C
C     Initialize XCBPRL - list of integral labels
C     ===========================================
C
      NPRPLBL    = 0
C
C     Initialize XCBXPR - list of properties
C     ======================================
      NPRPS        = 0
C
#ifdef MOD_KRMC
      CALL OPTINI()
#endif
#ifdef MOD_ERI
      CALL FCK3INI
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck haminp */
      SUBROUTINE HAMINP(WORK,LWORK)
C***********************************************************************
C
C     Input section for modification of Dirac-Coulomb Hamiltonian
C
C     Written by T.Saue - May 1996
C     Last revision: May 27 1998 - tsaue and kmi
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 28,NDIR = 8)
      PARAMETER (D0 = 0.00D00, D1 = 1.0D0,DM2 = -2.0D0,D2 = 2.00D00)
C
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "gencon.h"
#include "nuclei.h"
#if defined (MOD_DFT)
#include "dftcom.h"
Cdirac#include "functionals.h"
#endif
C
      LOGICAL SYMSTP,NOINPUT,SETFUN,SUCCESS
      CHARACTER PROMPT*1, WORD*7, TABDIR(NDIR)*7, TABLE(NTABLE)*7, 
     &          WORD1*7, LINE*80
      DIMENSION KSOP(0:7),WORK(LWORK)
C
      SAVE SETFUN
      DATA TABDIR /'*END OF','*SOLVEN','*DFT   ','XXXXXXX',
     &             'XXXXXXX','XXXXXXX','XXXXXXX','XXXXXXX'/
      DATA TABLE /'.PRINT ','.ONESYS','.OPERAT','.GAUNT ',
     &            '.LVCORR','.LVNEW ','.SOLVEN','.FREEPJ',
     &            '.VEXTPJ','.URKBAL','.NOSMLV','.LEVY-L',
     &            '.SPINFR','.ZORA  ','.SMLV1C','.ALDA  ',
     &            '.DFT   ','.HFXFAC','.SCQSET','.ONECNV',
     &            '.ONECAP','.SPINF2','.NOSPIN','.XXXXXX',
     &            '.INTFLG','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SETFUN/.FALSE./
C
#include "ibtfun.h"
C
      SYMSTP = .FALSE.
      INPERR = 0
      CALL TITLER('Hamiltonian defined','*',127)      
C
C
C     Initialize /CBIHAM/
C     ===================
C
      IPRHAM  = 0
      QED     = .FALSE.
      ONESYS  = .FALSE.
      FREEPJ  = .FALSE.
      VEXTPJ  = .FALSE.
      URKBAL  = .FALSE.
      NOSMLV  = .FALSE.
      LEVYLE  = .FALSE.
      SPINFR  = .FALSE.
      NOSPIN  = .FALSE.
      SPINFR2 = .FALSE.
      DOLVC   = .FALSE.
      LVNEW   = .FALSE.
      DFTHAM  = .FALSE.
      SOLVEN  = .FALSE.
      ZORA    = .FALSE.
      ZORA4   = .FALSE.
      ZORASC  = .FALSE.
      ONECAP  = .FALSE.
      SMLV1C  = .FALSE.
      ONECOFF = .FALSE.
      INTV1C  = 0
      ONECNV  = D0
      N1OPER  = 0
      SSMTRC  = D1
Ctsaue: Exchange factor must not be reset during geometry optimizations.
      IF (.NOT. SETFUN) HFXFAC  = D1
C ... get integral flags in advance
      ILLINT    = IBTAND(INTGEN,1)
      ISLINT    = IBTAND(INTGEN/2,1)
      ISSINT    = IBTAND(INTGEN/4,1)
      IGTINT    = IBTAND(INTGEN/8,1)

      CALL OP1INI(WORK,LWORK)
      NOINPUT = .FALSE.
#if defined (MOD_DFT)
C
C     Initialize from DFTCOM
C     ----------------------
C
      DFTRUN = .FALSE.
      DFTADD = .FALSE.
#endif
C
C     Read menu file
C     ==============
C     **** Find Hamiltonian input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      IF (WORD .EQ. '**HAMIL') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         NOINPUT = .TRUE.
         GOTO 300
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no HAMIL input found')
  930 CONTINUE
      WORD1 = WORD(1:1)
C
C     Process input for COMMON  /CBIHAM/
C     ==================================
C
C     Keep track of number of changes and added finite field operators
      I_FINFIELD = 0
      ICHANG  = 0
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124,125,126,127,128),I
            END IF
  99    CONTINUE
        IF (WORD .EQ. '.OPTION') THEN
           CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
           CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
           GO TO 100
         END IF
         WRITE (LUPRI,'(/,3A,/)')
     *       ' Keyword ',WORD,' not recognized in HAMINP.'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword in HAMINP.')
  101    CONTINUE
C&&&& PRINT: Print level for Hamiltonian
           READ(LUCMD,*) IPRHAM
           GO TO 100
  102    CONTINUE
C&&&& ONESYS: Ignore two-electron operator
           ONESYS = .TRUE.
           ICHANG = ICHANG + 1
           GO TO 100
  103    CONTINUE
C&&&& OPERATOR: Read in operator
           CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                 WORK,LWORK,IPRGEN)
           N1OPER = N1OPER + 1
           IF(N1OPER.GT.MAX1OPER) THEN
             WRITE(LUPRI,'(A)') 
     &       '***** HAMINP ERROR *****',
     &       'Redimension MAX1OPER of /CBIHAM/ !'
             CALL QUIT('HAMINP: beyond MAX1OPER !')
           ENDIF
           IPR1OP(N1OPER) = INDXPR
           ICHANG = ICHANG + 1
           I_FINFIELD = I_FINFIELD + 1
           GO TO 100
  104    CONTINUE
C&&&& GAUNT: Gaunt operator
CTROND          Must be implemented......                  
            CALL QUIT('.GAUNT not implemented')
            GO TO 100
  105     CONTINUE
C&&&& LVCORR: Model interatomic SS-integral contribution by classical
C             repulsion of small component atomic charges
            DOLVC = .TRUE.
            GO TO 100
  106     CONTINUE
C&&&& LVNEW: Get Small component charge via Mulliken analysis
            DOLVC = .TRUE.
            LVNEW = .TRUE.
            GO TO 100
  107     CONTINUE
C&&& SOLVENT: Solvent effects modelled by spherical cavity in dielectric medium
            SOLVEN = .TRUE.
            GO TO 100
  108     CONTINUE
C&&& FREEPJ: Reduce variational space by projecting out free positronic solutions
             FREEPJ = .TRUE.
             ICHANG = ICHANG + 1
            GO TO 100
  109     CONTINUE
C&&& VEXTPJ: Reduce variational space by projecting out positronic solutions of
C&&&         external field
             VEXTPJ = .TRUE.
             ICHANG = ICHANG + 1
             GO TO 100
  110      CONTINUE
C&&&& URKBAL : No restriction on kinetic balance; Unphysical solutions not removed
             URKBAL = .TRUE.
             GO TO 100
  111      CONTINUE
C&&& NOSMLV: Turn off small nuclear attraction integrals
             NOSMLV = .TRUE.
             ICHANG = ICHANG + 1
             GO TO 100
  112      CONTINUE
C&&& LEVYLE : Use Levy-Leblond Hamiltonian
             LEVYLE = .TRUE.
C            Levy-Leblond implies also spinfree, specify that
             SPINFR = .TRUE.
             NOSMLV = .TRUE.
             SSMTRC = D0
C            Check that we have not defined operators yet, as they need to be modified in Levy-Leblond calcs.
             IF (I_FINFIELD.GT.0) 
     &          CALL QUIT ('give .LEVY-LEBLOND option before .OPERAT')
C            No modifications of Hamiltonian, only new definition
C            ICHANG = ICHANG + 1
             GO TO 100
  113      CONTINUE
C&&& SPINFREE : Zero out spin-orbit terms
             SPINFR = .TRUE.
             ICHANG = ICHANG + 1
             GOTO 100
  114      CONTINUE
C&&& ZORA   : Use ZORA approximation
             READ(LUCMD,*) IZORA4,IZORAS
             ZORA   = .TRUE.
             ZORA4  = (IZORA4.NE.0)
             ZORASC = (IZORAS.NE.0)
             SSMTRC = D0
             ICHANG = ICHANG + 1
             GO TO 100
  115      CONTINUE
C&&& SMLV1C : Neglect potential for off-diagonal SS-blocks
C    jkp    : This keyword enables Luuk's one-center
C             projection scheme
             SMLV1C = .TRUE.
             ONECAP = .TRUE.
             INTV1C = 1
             ICTLV1C(2) = 1
             ICHANG = ICHANG + 1
             GO TO 100
  116      CONTINUE
C&&& ALDA   : Kohn-Sham (DFT) Hamiltonian w/ LDA functional for LR.
             DFTHAM = .TRUE.
#if defined (MOD_DFT)
             CALL DOALDA(.TRUE.)
#endif
             GOTO 100
  117      CONTINUE
C&&& DFT    : Kohn-Sham (DFT) Hamiltonian
#if defined (MOD_DFT)
               READ(LUCMD,'(A80)') LINE
               IF (.NOT. SETFUN) THEN
C              ... hjaaj: DFTINPUT accumulates def. of functional
C                  so if we call more than once (as in geom.opt.)
C                  we get n*functional for n'th call !!!! /3-Jun-2003 hjaaj
                  CALL DFTINPUT(LINE, SUCCESS)
                  IF(.NOT.SUCCESS) THEN
                     CALL DFTLISTFUNCS
                     CALL QUIT("DFT FUNCTIONAL Read failed.")
                  ENDIF
               END IF
               DFTHAM = .TRUE.
               DFTADD = .TRUE.
               DFTRUN = .TRUE.
               SETFUN = .TRUE.
#endif
            GO TO 100
  118     CONTINUE
C&&& HFXFAC: weight of pure exchange
             READ(LUCMD,*) HFXFAC
            GO TO 100
  119     CONTINUE
C&&& SCQSET: Reset small component charge for nuclei IZ (for onecap)
             READ(LUCMD,*) IZ, SCQVAL
             CALL SCQSET(IZ,SCQVAL)
            GO TO 100
  120     CONTINUE
C&&& ONECNV: Set threshold for turning off ONECAP model
            READ(LUCMD,*) ONECNV
            GO TO 100
  121     CONTINUE
C&&& ONECAP : Approximate potential for off-diagonal SS-blocks
C             if no model specified, goto 199, choose onecap 2
            READ(LUCMD,'(I10)',END=199,ERR=199) JNTV1C
            INTV1C = MOD(JNTV1C,10)
            ICTLV1C(2) = JNTV1C / 10
            GOTO 198
  199       INTV1C = 2
            ICTLV1C(2) = 0
            BACKSPACE LUCMD
  198       CONTINUE
            ONECAP  = .TRUE.
            IF(INTV1C.EQ.2) DOLVC = .TRUE.
C
C   ... Whether to evaluate integrals in a direct or concentional
C       matter has not been specified in **GENEREL
C       We write LS and SS integrals to disk by default in the
C       ONECAP type 2 & 3 model.            /jkp
C
C       Switched off in this version [DIRSET = .TRUE.]
C       - even one-center integrals become quite numerous in some cases!!!
            IF (.NOT.DIRSET .AND. INTV1C.NE.2) THEN
               ILLDIR = 1
               ISLDIR = 0
               ISSDIR = 0
               IDFLAG = ILLDIR+2*ISLDIR+4*ISSDIR
            END IF
            ICHANG = ICHANG + 1
            GO TO 100
  122     CONTINUE
C&&& SPINF2 : debug: Zero out 2-electron spin-orbit terms
            SPINFR2 = .TRUE.
            ICHANG = ICHANG + 1
            GO TO 100
  123     CONTINUE
C&&& NOSPIN: Take out all spin-interactions
             SPINFR = .TRUE.
             NOSPIN = .TRUE.
             ICHANG = ICHANG + 1
            GO TO 100
  124     CONTINUE
            GO TO 100
  125     CONTINUE
C&&& INTFLG: Specify what two-integrals should be included in this run
            IF (IGTINT .EQ. 1) THEN
               READ(LUCMD,*) ILLINT,ISLINT,ISSINT,IGTINT
            ELSE
               READ(LUCMD,*) ILLINT,ISLINT,ISSINT
            END IF
            GO TO 100
  126     CONTINUE
            GO TO 100
  127     CONTINUE
            GO TO 100
  128     CONTINUE
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
        GO TO 300
      ELSE
        WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *         '" not recognized in HAMINP.'
        CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
        CALL QUIT('Illegal prompt in HAMINP.')
      END IF
  300 CONTINUE
C
C     Remove inconsistencies
C     =====================
C
      IF (NOSMLV .AND. ONECAP) THEN
        IF(INTV1C.EQ.1.AND.ICTLV1C(2).EQ.1) THEN
          WRITE(LUPRI,'(A)')
     &  ' INFO: .SMLV1C ignored because .LEVY-L or .NOSMLV specified'
        ELSE
          WRITE(LUPRI,'(A)')
     &  ' INFO: .ONECAP ignored because .LEVY-L or .NOSMLV specified'
        ENDIF
        ONECAP = .FALSE.
      END IF
      IF (SMLV1C.AND.ONECAP.AND.(INTV1C.NE.1).AND.(ICTLV1C(2).NE.1))THEN
C     ... User asked for both SMLV1C and a non-equivalent ONECAP model
C         The ONECAP model overrides SMLV1C !!!
          WRITE(LUPRI,'(A)') ' INFO: .SMLV1C ignored - '//
     &   'another ONECAP model was specified.'
      ENDIF
      IF (LEVYLE) THEN
         IF(FREEPJ) CALL QUIT('Conflicting options : LEVYLE & FREEPJ !')
         IF(VEXTPJ) CALL QUIT('Conflicting options : LEVYLE & VEXTPJ !')
      ENDIF
C
      IF( (FREEPJ.OR.VEXTPJ.OR.LEVYLE) .AND. URKBAL )
     $   CALL QUIT('FREEPJ, VEXTPJ, or LEVYLE require URKBAL = .FALSE.')
C
      IF(FREEPJ.AND.VEXTPJ) THEN
        WRITE(LUPRI,'(A)') '***** HAMINP ERROR *****'
        WRITE(LUPRI,'(A)') 'You cannot specify both VEXTPJ and FREEPJ !'
        CALL QUIT('Conflicting options : VEXTPJ and FREEPJ !')
      ENDIF
C
C     Print section
C     =============
C
 400  CONTINUE
      WRITE(LUPRI,'(A,I5)') ' * Print level :',IPRHAM
C
      IF(LEVYLE) THEN
        WRITE(LUPRI,'(A)') ' * Levy-Leblond Hamiltonian'
      ELSE IF (IGTINT .NE. 0) THEN
        WRITE(LUPRI,'(A)') ' * Dirac-Coulomb-Gaunt Hamiltonian'
      ELSE
        WRITE(LUPRI,'(A)') ' * Dirac-Coulomb Hamiltonian'
      ENDIF
      IF(ICHANG.GT.0) WRITE(LUPRI,'(A)') 
     &    '   with the following modifications:'
      IF(SPINFR) THEN
         WRITE(LUPRI,'(3X,A)') '- Spin-orbit interactions neglected'
      ELSE IF(SPINFR2) THEN
         WRITE(LUPRI,'(3X,A)') '- 2-el. spin-orbit terms zeroed out '//
     &      '(NB : ONLY in the DHF calculation !!)'
      END IF
#if defined (MOD_DFT)
      IF(DFTHAM) CALL DFTREPORT()
#endif
      IF(DOLVC) ISSINT = 0
C
C     Set integral flags
C
      INTGEN = 1*ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
C
C     Write info about the modelling of SS interactions
C
      IF(DOLVC) THEN
C
C        Calculate LV-correction
C        (SS integrals modelled by Coulombic correction)
C        ===============================================
C
         WRITE(LUPRI,'(A)')
     &   ' * SS integrals neglected:',
     &   '   Interatomic Coulombic SS-contributions modelled by',
     &   '   classical repulsion of small component atomic charges'
         IF (LVNEW) THEN
            WRITE(LUPRI,'(A)')
     &   '   found from Mulliken population analysis.'
         ELSE
            WRITE(LUPRI,'(A)')
     &   '   using tabulated charges.'
         END IF
      END IF
C
C     Write information about default integrals flags
C
      WRITE(LUPRI,9001) '* Default integral flags passed to all modules'
      WRITE(LUPRI,9002) '- LL-integrals: ',ILLINT
      WRITE(LUPRI,9002) '- LS-integrals: ',ISLINT
      WRITE(LUPRI,9002) '- SS-integrals: ',ISSINT
      WRITE(LUPRI,9002) '- GT-integrals: ',IGTINT
 9001 FORMAT(1X,A)
 9002 FORMAT(3X,A,L1)
C
      WRITE(LUPRI,'(A)') ' * Basis set:'
      IF ( DIRCON(1) ) THEN
         WRITE(LUPRI,'(A)')'   - contracted large component basis set'
      ELSE
         WRITE(LUPRI,'(A)')'   - uncontracted large component basis set'
      END IF
      IF ( DIRCON(2) ) THEN
         WRITE(LUPRI,'(A)')'   - contracted small component basis set'
      ELSE
         WRITE(LUPRI,'(A)')'   - uncontracted small component basis set'
      END IF
      IF ( DIRCON(2) .AND. .NOT. URKBAL ) THEN
         WRITE(LUPRI,'(A)') '**** ERROR in HAMINP ****'
         WRITE(LUPRI,'(A)') 'Contracted small component basis set '//
     &      'and restricted kinetic balance is not implemented.'
         CALL QUIT('*** ERROR in HAMINP ****')
      END IF
C
C     Remove integral classes not needed for this Hamiltonian
C
      IF (LEVYLE) THEN
        ISLINT = 0
        ISSINT = 0
        IGTINT = 0
      ENDIF
      IF (ZORA.AND..NOT.ZORA4) THEN
        ISSINT = 0
      ENDIF
      INTGEN = 1*ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
C
      IF(URKBAL) THEN
        WRITE(LUPRI,1000) 'unrestricted'
      ELSE
        WRITE(LUPRI,1000) 'restricted'
      ENDIF
C
C     Projections on basis set
C
      IF(ISPHTR.EQ.3) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(A)')
     &      ' * Default RKB projection:',
     &      '   1: Pre-projection in scalar basis',
     &      '   2: Removal of unphysical solutions'//
     &      ' (via diagonalization of free particle Hamiltonian)'
         ELSE
            WRITE(LUPRI,'(A)')
     &      ' * RKB pre-projection in scalar basis',
     &      '  without removing unphysical'//
     &      ' solutions for no external field',
     &      '  This gives extended kinetic balance and is'//
     &      '  only recommended for CONTRACTED LC basis sets'
         ENDIF
      ELSEIF(ISPHTR.EQ.1) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(A)')
     &      ' * RKB projection:',
     &      '   1: No Pre-projection in scalar basis',
     &      '   2: Removal of unphysical solutions'//
     &      ' (via diagonalization of free particle Hamiltonian)',
     &      ' *** WARNING **** This procedure may give severe linear' 
     &    //'  dependence problems for large basis sets'
         ELSE
            WRITE(LUPRI,'(A)')
     &      '  Spherical large component basis and Cartesian'//
     &           ' small components ',
     &      ' *** WARNING **** This procedure may give severe linear' 
     &    //'  dependence problems for large basis sets'
         ENDIF
      ELSEIF(ISPHTR.EQ.0) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(A)')
     &      ' * Non-standard RKB projection:',
     &      '   1: Cartesian large components',
     &      '   2: Removal of unphysical solutions'//
     &      ' (via diagonalization of free particle Hamiltonian)'
         ELSE
            WRITE(LUPRI,'(A)') ' * Unrestricted Cartesian basis'
         ENDIF
      ELSE
         WRITE(LUPRI,'(/A,I5)')
     &        ' ERROR: ISPHTR not legal:',ISPHTR
         CALL QUIT('HAMINP: Invalid ISPHTR!')
      ENDIF
C
C
C   ... Check that the selected One-center model is reasonable,
C       and write appropriate output.
C
      IF(.NOT.LEVYLE .AND. ONECAP) THEN
         IF (INTV1C .LE. 0 .OR. INTV1C .GT. 4) THEN
            CALL QUIT('  Invalid one-center model specified.')
         END IF
         IF(INTV1C.EQ.1) THEN
            WRITE(LUPRI,'(1X,A)')
     &      '* Only one-center SS parts of the e-e repulsion and n-e',
     &      '  attraction are calculated explicitly.'
            IF (ICTLV1C(2) .EQ. 1)
     &          WRITE(LUPRI,'(1X,A)')
     &       '  - Using a projection scheme for the multicenter parts.'
            IF (ICTLV1C(2) .EQ. 2)
     &          WRITE(LUPRI,'(1X,A)')
     &       '  - Testing the simulation of terms neglected in this',
     &       '    way (cheap alternative to Luuk`s projection scheme)'
         ELSE IF (INTV1C .EQ. 2) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '* One-center approximation (model 2):',
     &          '  - Multicenter SS and LS integrals are modelled by',
     &          '    classical interaction of atomic charges'
             IF(ICTLV1C(2) .EQ. 0) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using Mulliken small component charges.'
             ELSE IF (ICTLV1C(2) .EQ. 1) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using tabulated small component charges.'
             ELSE IF(ICTLV1C(2) .EQ. 2) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using Mulliken small component charges.',
     &          '  - Estimating error associated with this model.'
             ELSE IF (ICTLV1C(2) .EQ. 3) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Testing HJJs penalty suggestion.'
             ELSE
                WRITE(LUPRI,'(1X,/A/)') 
     &          '    Wrong one-center model specification !'
                CALL QUIT('    HERRDN : ERROR in readin. ')
             END IF
         ELSE IF(INTV1C.EQ.3) THEN
             WRITE(LUPRI,'(1X,A)')
     &       '* One-center approximation (model 3):',
     &    '  - Multicenter SS and LS modelled by Coulomb attraction',
     &    '    integrals using large and small component charges',
     &    '    based on the density matrix.'
             IF(ICTLV1C(2) .EQ. 0) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Charges are distributed on the nuclei and on points',
     &    '    between the nuclei where the contribution to the',
     &    '    interaction from the dipole moment is zero.'
             ELSE IF(ICTLV1C(2) .EQ. 1) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Option not available (no tabulated large component '//
     &         'charges).'
               CALL QUIT('    HERRDN : ERROR in readin. ')
             ELSE IF(ICTLV1C(2) .EQ. 2) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Using atomic centered Mulliken large and small'//
     &        ' component charges.'
             ELSE IF(ICTLV1C(2) .EQ. 3) THEN
                WRITE(LUPRI,'(1X,A)')
     &    '  - Charges are distributed on the nuclei and on the ',
     &    '    midpoints between the nuclei.'
             ELSE
                WRITE(LUPRI,'(1X,/A/)') 
     &          '    Invalid one-center model specified !'
                CALL QUIT('    ONECAP : ERROR in readin. ')
             END IF
         ELSE IF(INTV1C.EQ.4) THEN
             WRITE(LUPRI,'(1X,A)')
     &       '* One-center approximation (model 4):',
     &    '  - All SS and LS integrals modelled by Coulomb attraction',
     &    '    integrals using small component atomic charges'
             IF(ICTLV1C(2) .EQ. 0) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using Mulliken small component charges.'
             ELSE IF (ICTLV1C(2) .EQ. 1) THEN
                WRITE(LUPRI,'(1X,A)')
     &          '  - Using tabulated small component charges.'
             ELSE
                WRITE(LUPRI,'(1X,/A/)') 
     &          '    Wrong ONECAP type 4 specification !'
                CALL QUIT('    HERRDN : ERROR in readin. ')
             END IF
         END IF
         IF (.NOT.DIRSET .AND. INTV1C.NE.2) THEN
             WRITE(LUPRI,'(1X,A)')
     &    '  - Default direct evaluation of all integrals changed',
     &    '    by ONECAP since nothing was specified in **GENEREL.',
     &    '    -->  LS and SS integrals are written to disk.'
         END IF
         IF(ONECNV .GT. D0) THEN
           WRITE(LUPRI,'(1X,A,/1X,A,1P,D15.6)')
     &  '  - One-center app. will be turned off when wavefunction is',
     &  '    converged beyond the specified threshold of : ONECNV = ',
     &       ONECNV
         END IF
C
      END IF
C
      IF(N1OPER.GT.0) THEN
        WRITE(LUPRI,'(3X,A)') '- The following one-electron '//
     &    'operators have been added:'
        DO I = 1,N1OPER
          INDXPR = IPR1OP(I)
          CALL PRSYMB(LUPRI,'-',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &     'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'-',75,0)
          CALL WRIXPR(I,INDXPR)
C         Operator must be totally symmetric
          ISYXPR = IPRPSYM(INDXPR)-1
          IF(ISYXPR.NE.0) THEN
             SYMSTP = .TRUE.
             WRITE(LUPRI,'(4X,A)') '***** HAMINP ERROR *****'
             WRITE(LUPRI,'(4X,A,A)') 
     &         'The operator is not totally symmetric under ',
     &         'the molecular point group'
               CALL REDSYM(ISYXPR,IGROUP,KSOP,NELM)
             WRITE(LUPRI,'(4X,A)') '* Possible actions:'
             WRITE(LUPRI,'(6X,A)') 
     &       '1. Include operator as perturbation.'
             WRITE(LUPRI,'(6X,A,A3,A,A3,10(1X,A3))')
     &       '2. Reduce symmetry: ',GROUP,' --> ',
     &       GROUPS(IGROUP),' :{',
     &       (SYMOP(KSOP(J-1)),J=1,NELM),'} '
          ENDIF
C         Operator must be symmetric under time reversal
          ITRXPR = IPRPTIM(INDXPR)
          IF(ITRXPR.NE.1) THEN
            SYMSTP = .TRUE.
            WRITE(LUPRI,'(4X,A)') '***** HAMINP ERROR *****'
            WRITE(LUPRI,'(4X,A,A)') 
     &         'The operator is not totally symmetric under ',
     &         'time reversal'
             WRITE(LUPRI,'(4X,A)') '* Possible actions:'
             WRITE(LUPRI,'(6X,A)') 
     &       '1. Include operator as perturbation.'
          ENDIF
        ENDDO
      ENDIF      
      IF(ONESYS) THEN
        WRITE(LUPRI,'(3X,A)') '- Two-electron part ignored '//
     &   '(one-electron system)'
      ENDIF
      IF (.NOT.LEVYLE) THEN
         IF(FREEPJ) THEN
           WRITE(LUPRI,'(3X,A)')
     &        '- Free positronic solutions projected out'
         ENDIF
         IF(VEXTPJ) THEN
           WRITE(LUPRI,'(3X,A)') '- External potential positronic '//
     &          'solutions projected out'
         ENDIF
         IF(ZORA.AND..NOT.ZORA4) WRITE(LUPRI,'(3X,A,1P,D8.2)')
     &      '- ZORA approximation using approximate density'
         IF(ZORA4) WRITE(LUPRI,'(3X,A,1P,D8.2)')
     &      '- ZORA approximation using the full density'
         IF(ZORASC) WRITE(LUPRI,'(3X,A,1P,D8.2)')
     &      '- Calculating and printing scaled ZORA orbital energies'
         IF(NOSMLV) WRITE(LUPRI,'(3X,A)')
     &      '- No nuclear attraction integrals for small components !'
      ENDIF
C
      IF(PROJEC) THEN
        WRITE(LUPRI,'(A)')
     &      ' * Variational space reduced by projection'//
     &      '  using a set of reference orbitals.'
        WRITE(LUPRI,'(3X,A,D12.4)') '- Smallest norm accepted: ',PRJTHR
        IF(PROOWN) THEN
          WRITE(LUPRI,'(3X,A)') 
     &     '- Reference orbitals are obtained from local bases !',
     &     '  The indexing of reference spinor sets is assumed ',
     &     '  to follow list of symmetry independent nuclei.'
          DO J = 1,NPRJREF
            WRITE(LUPRI,'(3X,A,I3,A,A6)') '- Reference spinor set ',J,
     &             ' --> ',NAMDEP(J)
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)') 
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)') 
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO J = 1,NPRJREF
            WRITE(LUPRI,'(3X,A,I3)') '- Reference spinor set ',J
            WRITE(LUPRI,'(3X,A,A6)') '  read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(5X,A,A3)') 
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(5X,A,A3,A,A72)') 
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF(PROJEC) THEN
        WRITE(LUPRI,'(1X,A,A)')
     &      '* Variational space reduced by projection',
     &      '  using a set of reference orbitals.'
        WRITE(LUPRI,'(3X,A,D12.4)') '- Smallest norm accepted: ',PRJTHR
        IF(PROOWN) THEN
          WRITE(LUPRI,'(3X,A)') 
     &     '- Reference orbitals are obtained from local bases !',
     &     '  The indexing of reference spinor sets is assumed ',
     &     '  to follow list of symmetry independent nuclei.'
          DO J = 1,NPRJREF
            WRITE(LUPRI,'(3X,A,I3,A,A6)') '- Reference spinor set ',J,
     &             ' --> ',NAMDEP(J)
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)') 
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)') 
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO J = 1,NPRJREF
            WRITE(LUPRI,'(3X,A,I3)') '- Reference spinor set ',J
            WRITE(LUPRI,'(5X,A,A6)') ' - read from file ',PRJFIL(J)
            DO I = 1,NFSYM
              NVEC = 0
              CALL  NUMLST(VCPROJ(I,J),IDUMMY,NFBAS(I,0),
     &                  -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
              IF(NVEC.EQ.0) THEN
                WRITE(LUPRI,'(7X,A,A3)') 
     &             '- No orbitals in fermion ircop ',FREP(I)
              ELSE
                WRITE(LUPRI,'(7X,A,A3,A,A72)') 
     &          '- Orbitals in fermion ircop ',FREP(I),' :',VCPROJ(I,J)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF (NOINPUT) GOTO 500
C
C     Process input for various program sections
C     ==========================================
C
  200 PROMPT = WORD(1:1)
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
         CALL QUIT('Illegal directory in HAMINP.')
      ELSE
         WRITE (LUPRI,'(/,4A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in HAMINP, error in prompt.')
      END IF
    2 CONTINUE
        CALL SOLVIN(WORD,WORK,LWORK)
        GO TO 200
    3 CONTINUE
#if defined (MOD_DFT)
        CALL DFTINP(WORD,.FALSE.,WORK,LWORK)
#else
C       hjaaj 3jun2003: dirac.x went into an infinite loop here
C       if *DFT in input but MOD_DFT not defined !!!
        CALL QUIT(
     &  '*DFT specified, but DFT is not included in this version')
#endif
        GO TO 200
    4 CONTINUE
        GO TO 200
    5 CONTINUE
        GO TO 200
    6 CONTINUE
        GO TO 200
    7 CONTINUE
        GO TO 200
    8 CONTINUE
        GO TO 200
C
    1 CONTINUE
C
#if defined (MOD_DFT)
Chjaaj: we need to be sure DFT is initialized (and SET reset for geom.opt.)
CJuly2003
      CALL DFTINP(WORD,.TRUE.,WORK,LWORK)
#endif
      IF(SOLVEN) CALL SOLVIN(WORD,WORK,LWORK)
C
 500  CONTINUE
C
      IF(SYMSTP) CALL QUIT('HAMINP: Non-symmetric operator !')
C
      RETURN
C
 1000 FORMAT (//' Information about the ',A,' kinetic balance scheme:')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1ini */
      SUBROUTINE OP1INI(WORK,LWORK)
C***********************************************************************
C
C     Initialize one-electron integrals for standard 
C     Dirac-Coulomb Hamiltonian
C
C     Written by T.Saue May 21 1996
C     Last reviosn : may 21 1996
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,DM2 = -2.0D0,DHALF=0.5D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3),WORK(LWORK)
C
C    Overlap integrals
C    =================
C
      PNAME            = 'Overlap matrix  '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'OVERLAP '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRHAM)
      IPOVRLAP        = INDXPR
C
C     Nuclear attraction integrals
C     ============================
C
      PNAME            = 'Nuc. attraction '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'MOLFIELD'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRHAM)
      IPMOLFLD         = INDXPR
C
C     Beta matrix integrals
C     =====================
C
      PNAME            = 'Beta matrix    '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = DM2*CVAL*CVAL
      PLABEL(1)        = 'BETAMAT '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRHAM)
      IPBETAMT         = INDXPR
C
C     Kinetic energy integrals
C     ========================
C
      PNAME            = 'Kinetic energy  '
      IPTYP            = 8
      NPCOMP           = 3
      PFAC(1)          = -CVAL
      PFAC(2)          = -CVAL
      PFAC(3)          = -CVAL
      PLABEL(1)        = 'XDIPVEL '
      PLABEL(2)        = 'YDIPVEL '
      PLABEL(3)        = 'ZDIPVEL '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRHAM)
      IPKINERG         = INDXPR
      IF (LINEAR) THEN
C
C        Angular momentum integrals (used in case of atomic or linear symmetry !)
C        ========================================================================
C
         PNAME            = 'Orbital z-momentum'
         IPTYP            = 1
         NPCOMP           = 1
         PFAC(1)          = -D1
         PLABEL(1)        = 'ZANGMOM '
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                  INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRHAM)
         IPANGMOM         = INDXPR
C
         PNAME            = 'Spin z-momentum'
         IPTYP            = 12     
         NPCOMP           = 1
         PFAC(1)          = DHALF
         PLABEL(1)        = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                  INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRHAM)
         IPSPNMOM         = INDXPR
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck CIMINP */
      SUBROUTINE CIMINP(WORK,LWORK)
C***********************************************************************
C
C     Dummy routine while the interface is being implemented
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TWOMEM */
      SUBROUTINE TWOMEM(WORK,LWORK)
C***********************************************************************
C
C     Check for memory for processing of two-electron integrals
C     Written by T.Saue June 12 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
C
#include "dcbgen.h"
#include "cbihr2.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL LBIT,DIRFCK
C
      CALL TITLER('Memory required by HERMIT','*',127)      
      CALL PAOVEC(WORK,LWORK,0,IPRTWO)
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') 
     &   ' Maximum memory load for two-electron integral processing:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     always DIRFCK false to calculate memory for SO-integrals
C     which is what we normally use in Dirac (both SOFOCK for
C     Fock matrices and for MOLTRA). DIRFCK true will calculate
C     the correct for AOFOCK only. /hjaaj May 2001
C
      DIRFCK = .FALSE.
      WRITE(LUPRI,'(A)') ' * LL-integrals:'
      CALL MEMLOP(1,DIRFCK)  
      WRITE(LUPRI,'(A)') ' * SL-integrals:'
      CALL MEMLOP(2,DIRFCK)  
      WRITE(LUPRI,'(A)') ' * SS-integrals:'
      CALL MEMLOP(3,DIRFCK)  
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MINANA */
      SUBROUTINE MINANA(WORK,LWORK)
C***********************************************************************
C
C     Readin analysis input for the final geometry
C     Call analysis module
C
C     Written by J. Thyssen, Nov 13 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('MINANA')
C
      IF(DOANA) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL ANWINP('**ANA F',WORK,LWORK)
         CLOSE(LUCMD)
C
C***********************************************************************
C*****  A N A L Y S I S    M O D U L E  ********************************
C***********************************************************************
C
         IF(DOANA) CALL PAMANA(WORK,LWORK)
C
      ENDIF
C
      CALL QEXIT('MINANA')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prolab */
      SUBROUTINE PROLAB(ILAB,LAB,NDIM,NLAB,LINE,IUNIT,IBUF,WORK,LWORK)
C***********************************************************************
C
C     Process labels:
C       Read a string of new label name and list of primitive
C       labels. 
C         From pointer ILAB relating proimitve labels with new
C       labels.
C
C     Written by T.Saue - Feb 3 1998
C
C***********************************************************************
#include"implicit.h"
#include"priunit.h"
C
      CHARACTER LINE*72,LAB(*)*12
      DIMENSION ILAB(*),IBUF(*),WORK(*)
C
      CALL IZERO(ILAB,NDIM)
      DO I = 1,NLAB
        NGRP = 1
        LINE = ' '
        READ(IUNIT,'(A72)') LINE
        LAB(I) = LINE(1:12)
        CALL NUMLST(LINE(12:72),IBUF,NDIM,1,NDIM,1,NGRP,WORK,LWORK)
        DO J = 1,NGRP
          ILAB(IBUF(J)) = I
        ENDDO
      ENDDO
C
      RETURN
      END
#ifdef NOT_USED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck READREORDER*/
      SUBROUTINE READREORDER(STRING,IREORD,IMOORD,MXREORD)
C***********************************************************************
C
C     Read in new order of orbitals
C
C     Written by J. Thyssen, Sep 1 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION IMOORD(MXREORD)
      CHARACTER STRING*72

      string(72:72) = char(0)
      call creadreorder(string,ireord,imoord,mxreord)
      RETURN
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mvoinp */
      SUBROUTINE MVOINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for generation of modified virtual orbitals
C
C     Written by T.Saue - June 1998
C     Last revision: June 16 1998
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
#include "dcbmvo.h"
C
      LOGICAL SET, NEWDEF, RESET,LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.VECMVO','.PRINT ','.INTFLG','.F2CONT',
     &            '.ADDREP','.IONMVO','.POTMVO','.XXXXXX'/
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
C     Initialize /CBIMVO/
C     ===================
C
      IPRMVO = 0
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      IFCMVO = 1
      IFEMVO = 1
      ADDREP = .FALSE.
      IONMVO = .TRUE.
      DO I = 1,NFSYM
        VECMVO(I) = ' '
        IF(NOCC(I).EQ.0) THEN
          WRITE(VECMVO(I),'(I1)') 1
        ELSEIF(NOCC(I).GT.0) THEN
          WRITE(VECMVO(I),'(I1,A2,I4)') 1,'..',NOCC(I)
        ENDIF
      ENDDO
C
C     Process input from CBIMVO
C     =========================
C
      NEWDEF = (WORD .EQ. '*MVOCAL')
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
     *            '" not recognized in MVOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in MVOINP.')
    1          CONTINUE
C&&&& VECMVO: Spinors defining potential
                  READ(LUCMD,'(A72)') (VECMVO(I),I=1,NFSYM)
               GO TO 100
    2          CONTINUE
C&&&& PRINT: General print level
                  READ(LUCMD,*) IPRMVO
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
                  READ(LUCMD,*) IFCMVO,IFEMVO
               GO TO 100
    5          CONTINUE
C&&& ADDREP: Add remaining occupied orbitals as repulsion
                  ADDREP = .TRUE.
               GO TO 100
    6          CONTINUE
C&&& IONMVO: Generate MVOS by ionic Fock operator
                  IONMVO = .TRUE.
               GO TO 100
    7          CONTINUE
C&&& POTMVO: Generate MVOS from Hartree-Fock potential
                  IONMVO = .FALSE.
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in MVOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in MVOINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      INTMVO = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
      IF2MVO = IFCMVO + 2*IFEMVO
C
C     Print section
C     =============
C
      IF(.NOT.DOMVO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' MVOGEN: Modified virtual orbitals'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I3)') ' * General print level:',IPRMVO
      IF(IONMVO) THEN
        WRITE(LUPRI,'(A)') ' * MVOs generated by ionic Fock operator'
      ELSE
        WRITE(LUPRI,'(A)') ' * MVOs generated by from DHF potential:'
      ENDIF
      WRITE(LUPRI,'(A)') ' * Vectors defining the DHF potential:'
      NTOT = 0
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECMVO(I),IDUMMY,NFBAS(I,0),
     &            -NFBAS(I,2),NFBAS(I,1),I,NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)') 
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          NTOT = NTOT + NVEC
          WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &    '- Orbitals in fermion ircop ',FREP(I),' :',VECMVO(I)
        ENDIF
      ENDDO
      IF(ADDREP.AND..NOT.IONMVO) THEN
        WRITE(LUPRI,'(A)')
     &     ' * Remaining occupied orbitals added as repulsion.'
      ENDIF
      IF(ILLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No LL contributions'
      IF(ISLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SL contributions'
      IF(ISSINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SS contributions'
      IF(IGTINT.EQ.0) WRITE(LUPRI,'(A)') ' * No GT contributions'
      IF(IF2MVO.NE.3) THEN
        WRITE(LUPRI,'(A)')' * Two-electron contributions restricted to:'
        IF(LBIT(IF2MVO,1)) WRITE(LUPRI,'(3X,A)')
     &     '+ direct contribution'
        IF(LBIT(IF2MVO,2)) WRITE(LUPRI,'(3X,A)')
     &     '+ exchange contribution'
      ENDIF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck coninp */
      SUBROUTINE CONINP
C***********************************************************************
C
C <<< Check concistency of the specified input >>>
C
C     Written by Luuk Visscher may 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "mxmcgas.h"
#include "dcbkrmc.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbana.h"
#include "dcbprp.h"
#include "dgroup.h"
      LOGICAL SS_COR, SS_COR_NEEDED
#include "ibtfun.h"
C
      CALL TITLER('Input consistency checks','*',127)
      IF (SPINFR.AND..NOT.LEVYLE.AND.DOMP2)
     &   CALL QUIT ('Spinfree MP2 not implemented')
      IF (SPINFR.AND.DOPRP)WRITE (LUPRI,1000)'Spinfree properties'
      IF (DORES .AND.DODFT)
     &   CALL QUIT ('.RESOLVE does not work with DFT')
      IF (DOXLR .AND.DODFT)
     &   CALL QUIT ('Linear response does not work with DFT')

C hjaaj April 2002: why not DORES .and. DOTRA ??? TODO
C     the input routines use different common blocks,
C     so I think running both should be OK. I have modified the test:
C     if DORES or DOMP2 integral transformation, then DOTRA necessary
C     for DOCCM or DOCIM as integrals on disk now cannot be correct.
Chj   IF (DORES.AND.DOTRA)
      IF ( (DORES.OR.DOMP2) .AND. ((DOCCM.OR.DOCIM) .AND. .NOT.DOTRA) )
     &CALL QUIT('4INDEX for CI/CC must be specified when RESOLV or MP2')
C
C     Correct input from the naive user who neglects the SS without applying a correction.
C     The experienced user may subtract the LVCOR later on.
C
C     Check for Hamiltonians were SS do not appear anyway..
      SS_COR_NEEDED = .NOT.(LEVYLE.OR.(ZORA.AND..NOT.ZORA4))
C     .. or are computed exactly
      SS_COR_NEEDED = SS_COR_NEEDED.AND.(IBTAND(INTGEN/4,1).EQ.0)
C     Check for possible ways of correcting..
      SS_COR = DOLVC.OR.DOTSC
      IF (SS_COR_NEEDED.AND..NOT.SS_COR) THEN
         DOLVC = .TRUE.
         WRITE (LUPRI,'(A/A)') 
     &     " * WARNING: User asked for neglect of SS-integrals "//
     &     "without invoking a correction.",
     &     "   Activated SCC method (.LVCORR), "//
     &     "see L. Visscher TCA 98 (1997) 68."
      ENDIF
C
C     Check for LUCIA and parallel
C
      IF (PARCAL .AND. DOKRMC) THEN
         IF (KRMC_CIPROGRAM .EQ. 'LUCIAREL')
     &   CALL QUIT ('LUCIAREL cannot run in parallel mode')
      ENDIF
C
      RETURN
 1000 FORMAT (/' **** WARNING **** ',A,' not tested.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck procc */
      SUBROUTINE PROCC
C***********************************************************************
C
C     Print DHF occupation information
C
C     Written by J. Thyssen Aug 12 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D2=2.0D00)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      IF (NAELEC_DHF .EQ. 0) THEN
         WRITE(LUPRI,'(1X,A,I4,A)') 
     &      '* Closed shell DHF calculation with ',NELECT_DHF,
     &      ' electrons in'
         DO I = 1,NFSYM
            WRITE(LUPRI,'(3X,I5,A,I1)')
     &          NISH_DHF(I),' orbitals in Fermion irrep ',I
         END DO
      ELSE
         WRITE(LUPRI,'(A/A//40X,A,2(/24X,A))') 
     &      ' * Open shell DHF calculation.',
     &      ' * Shell specifications:',
     &      'Orbitals',
     &      '#electrons  irrep 1  irrep 2     f        a      alpha',
     &      '----------  -------  -------  -------  -------  -------'
         IF (NFSYM .EQ. 1) THEN
            WRITE(LUPRI,1001) 
     &         NELECT_DHF,NISH_DHF(1),DF(0),DA(0),DALPHA(0)
         ELSE
            WRITE(LUPRI,1002) 
     &         NELECT_DHF,NISH_DHF(1),NISH_DHF(2),DF(0),DA(0),DALPHA(0)
         END IF
         DO IOPEN = 1,NOPEN
            IF (NFSYM .EQ. 1) THEN
               WRITE(LUPRI,1011) 
     &            IOPEN,
     &            NINT(DF(IOPEN)*D2*(NACSH(1,IOPEN)+NACSH(2,IOPEN))),
     &            NACSH(1,IOPEN),DF(IOPEN),DA(IOPEN),DALPHA(IOPEN)
            ELSE
               WRITE(LUPRI,1012) 
     &            IOPEN,
     &            NINT(DF(IOPEN)*D2*(NACSH(1,IOPEN)+NACSH(2,IOPEN))),
     &            NACSH(1,IOPEN),NACSH(2,IOPEN),
     &            DF(IOPEN),DA(IOPEN),DALPHA(IOPEN)
            END IF
         END DO
         WRITE(LUPRI,1020)
         IF (NFSYM .EQ. 1) THEN
            WRITE(LUPRI,1031)
     &         NELECT_DHF+NAELEC_DHF,NISH_DHF(1)+NASH_DHF(1)
         ELSE
            WRITE(LUPRI,1032)
     &         NELECT_DHF+NAELEC_DHF,NISH_DHF(1)+NASH_DHF(1),
     &         NISH_DHF(2)+NASH_DHF(2)
         END IF
         WRITE(LUPRI,'(/,3X,A,/)')
     &      'f is the fraction occupation; a and alpha open shell '//
     &      'coupling coefficients.'
      END IF
      RETURN
 1001 FORMAT(3X,'Closed shell',12X,I4,7X,I4,6X,'N/A ',3(3X,F6.4))
 1002 FORMAT(3X,'Closed shell',12X,I4,7X,I4,5X,I4,1X,3(3X,F6.4))
 1011 FORMAT(3X,'Open shell no.',I2,8X,I4,7X,I4,6X,'N/A ',3(3X,F6.4))
 1012 FORMAT(3X,'Open shell no.',I2,8X,I4,7X,I4,5X,I4,1X,3(3X,F6.4))
 1020 FORMAT(3X,76('-'))
 1031 FORMAT(3X,'Total',19X,I4,7X,I4)
 1032 FORMAT(3X,'Total',19X,I4,7X,I4,5X,I4)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rho1in */
      SUBROUTINE RHO1IN(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for 1-dimensional density plot
C
C     Written by T.Saue - May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbrho1.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.MESH  ','.PRINT ','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
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
C     Initialize /DCBRHO1/
C     ===================
C
      DSTEP  = 0.01D0
      IPRHO1 = 0
C
C     Process input from DCBRHO1
C     ==========================
C
      NEWDEF = (WORD .EQ. '*RHO1  ')
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
     *            '" not recognized in RHO1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DHFINP.')
    1          CONTINUE
C&&&& MESH; Step length for grid
                  READ(LUCMD,*) DSTEP
               GO TO 100
    2          CONTINUE
C&&&& PRINT: Print level
                  READ(LUCMD,*) IPRHO1
               GO TO 100
    3          CONTINUE
               GO TO 100
    4          CONTINUE
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
     *            '" not recognized in RHO1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in RHO1IN.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DO1RHO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' RHO1IN: 1-dimensional density plot:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I4)') ' * Print level: ',IPRHO1
      WRITE(LUPRI,'(A,F7.4)') ' * Grid mesh (in Angstroms): ',DSTEP
      WRITE(LUPRI,'(A)') 
     &    ' (Remember that lines can be defined by ghost centers)'
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densin */
      SUBROUTINE DENSIN(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for 3-dimensional density plot
C
C     Written by T.Saue - Nov 5 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbrho.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.NCUBE ','.PRINT ','.CUBADJ','.DOCUBE',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
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
C     Initialize /DCBRHO1/
C     ===================
C
      ILSCUB = 3
      IPRRHO = 0
      NCUBE(1) = 80
      NCUBE(2) = 80
      NCUBE(3) = 80
      CUBADJ(1) = 4.0D0
      CUBADJ(2) = 8.0D0
C
C     Process input from DCBRHO
C     =========================
C
      NEWDEF = (WORD .EQ. '*DENSIT')
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
     *            '" not recognized in DENSIN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DENSIN.')
    1          CONTINUE
C&&&& NCUBE: Number of points along the sides of the cube
               READ(LUCMD,*) NCUBE(1),NCUBE(2),NCUBE(3)
               GO TO 100
    2          CONTINUE
C&&&& PRINT: Print level
                  READ(LUCMD,*) IPRRHO
               GO TO 100
    3          CONTINUE
C&&&& CUBADJ: Adjust size of cube to hopefully include 'most' density.
C&&&&         CUBADJ(1) : adjust startpoint of cube
C&&&&         CUBADJ(2) : adjust stepvectors of cube
               READ(LUCMD,*)  CUBADJ(1), CUBADJ(2)
               GO TO 100
    4          CONTINUE
C&&&  DOCUBE: Specify generation of cube-file for L- and S-comp.
               READ(LUCMD,*) ILC,ISC
               ILSCUB = ILC + 2*ISC
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
     *            '" not recognized in DENSIN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in DENSIN.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DO3RHO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'DENSITY: 3-dimensional density plot:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I4)') '* Print level: ',IPRRHO
  999 CONTINUE
      RETURN
C
      END
#if defined (MOD_DFT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wt1in */
      SUBROUTINE WT1IN(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for 1-dimensional plot of Becke weights
C
C     Written by T.Saue - May 9 2000
C
C     Modified for WT1IN by O. Fossgaard - May 2003
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbana.h"
#include "dcbwt1.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.MESH  ','.PRINT ','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
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
C     Initialize /DCBWT1/
C     ===================
C
      DSTEP  = 0.01D0
      IPWT1 = 0
C
C     Process input from DCBBWT1
C     ==========================
C
      NEWDEF = (WORD .EQ. '*WT1   ')
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
     *            '" not recognized in WT1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DHFINP.')
    1          CONTINUE
C&&&& MESH; Step length for grid
                  READ(LUCMD,*) DSTEP
               GO TO 100
    2          CONTINUE
C&&&& PRINT: Print level
                  READ(LUCMD,*) IPWT1
               GO TO 100
    3          CONTINUE
               GO TO 100
    4          CONTINUE
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
     *            '" not recognized in WT1IN.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in WT1IN.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.DO1WT) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' WT1IN: 1-dimensional plot of Becke weights:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I4)') ' * Print level: ',IPWT1
      WRITE(LUPRI,'(A,F7.4)') ' * Grid mesh (in Angstroms): ',DSTEP
      WRITE(LUPRI,'(A)') 
     &    ' (Remember that lines can be defined by ghost centers)'
  999 CONTINUE
      RETURN
C
      END
#endif
