C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirrdn.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck paminp */
      SUBROUTINE PAMINP(WORK,LWORK)
C*****************************************************************************
C     
C     Driver for input module of DIRAC
C     Written by T.Saue April 26 1996
C     Last revision: 260496 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dgroup.h>
#include <cbirea.h>
      DIMENSION WORK(*)
C    
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
      CALL SETDC1(IPREAD)
C
C     Define Hamiltonian
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
      IF(DOPSI.OR.MINIMI) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL PSIINP(WORK,LWORK)
         CLOSE(LUCMD)
      ENDIF
C
C     Input of wave function analysis module
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
#if defined (MOD_TRA)
C
C     Input of transformation module
C
      IF(DOTRA) THEN
        OPEN(LUCMD,FILE = 'DIRAC.INP')
        CALL TRAINP(WORK,LWORK)
        CLOSE(LUCMD)
      ENDIF
#endif
#if defined (MOD_CC)
C
C     Input of relccsd module
C
      IF(DOCCM) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL CCINPT(LUCMD,LUPRI)
         CLOSE(LUCMD)
      ENDIF
#endif
#if defined (MOD_CI)
C
C     Input of dirrci module
C
      IF(DOCIM) THEN
         OPEN(LUCMD,FILE = 'DIRAC.INP')
         CALL CIMINP(WORK,LWORK)
         CLOSE(LUCMD)
      ENDIF
#endif
C
      CALL TITLER('End of input processing','*',124)
C
C     If input test , exit here
C     =========================
C
      IF(INPTES) STOP '***End of DIRAC input test ***'
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck geninp */
      SUBROUTINE GENINP(WORK,LWORK)
C*****************************************************************************
C
C <<< General Input for DIRAC >>>
C
C       This routine is based on the corresponding HERINP routine
C       in HERMIT
C
C       Written by Trond Saue, Tromsoe October 1993
C       Last revision: Nov 14 1994 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <alphac.h>
#include <mxcent.h>
#include <mxorb.h>
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 8, NTABLE = 24)
C
      LOGICAL LBIT
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7,FMT*6,MXFORM*6
      DIMENSION WORK(LWORK)
C
#include <dcbgen.h>
#include <dcbxpr.h>
#include <dcbcls.h>
#include <dcbprl.h>
C
      DATA TABDIR /'*END OF','*PARALL','*XXXXXX','*XXXXXX',
     &             '*XXXXXX','*XXXXXX','*XXXXXX','*XXXXXX'/
      DATA TABLE  /'.DIRECT','.CVALUE','.PARALL','.DOJACO',
     &             '.PRINT ','.PCMOIN','.PCMOUT','.LINDEP',
     &             '.SPHTRA','.XXXXXX','.XXXXXX','.XXXXXX',
     &             '.XXXXXX','.VACUUM','.NUCVAC','.FAMILY',
     &             '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &             '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
C     Read menu file
C     ==============
C     **** Find General input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 CONTINUE
      READ (LUCMD,'(A7)',END=180,ERR=920) WORD
      IF (WORD .EQ. '**GENER') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
C  910 CONTINUE
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
C&&&& PARALL: Parallel Direct calculation
            PARCAL = .TRUE.
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
  110    CONTINUE
            GO TO 100
  111    CONTINUE
            GO TO 100
  112    CONTINUE
            GO TO 100
  113    CONTINUE
            GO TO 100
  114    CONTINUE
C&&&& VACUUM:QED-Fock (vacuum polarization included)
            VACUUM = .TRUE.
            GO TO 100
  115    CONTINUE
C&&&& NUCVAC: Uehling - potential
            NUCVAC = .TRUE.
            GO TO 100
  116    CONTINUE
C&&&& FAMILY: Family basis (Not yet implemented)
            FAMILY = .TRUE.
            GO TO 100
  117    CONTINUE
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
         IDFLAG = 7
      ENDIF
C
C     Print section
C     =============
C
      CALL TITLER('General DIRAC set-up','*',126)
      WRITE(LUPRI,'(1X,A,F18.7)')
     +   '* The speed of light : ',CVEL
      IF(CVAL.NE.CVEL) THEN
        FMT = MXFORM(CVAL,14)
        WRITE(LUPRI,'(3X,A,E9.3,A,'//FMT//')')
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
      ENDIF
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
      WRITE(LUPRI,'(1X,A,2(/3X,A,E9.3))')
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
#if defined (VAR_MPI) || defined (VAR_PVM)
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
#if defined (VAR_MPI) || defined (VAR_PVM)
      CALL PRLINP(WORD)
#endif
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck jobinp */
      SUBROUTINE JOBINP
C*****************************************************************************
C
C <<< Job input for DIRAC >>>
C
C     Written by Trond Saue May 1996
C     Last revision: 1997/09/15 - jth (MINIMIZE added)
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (NTABLE = 12, NDIR = 4)
C
      CHARACTER WORD*7, PROMPT*1, TABLE(NTABLE)*7, TABDIR(NDIR)*7,
     &          WORD1*7
C
#include <dcbgen.h>
C
      DATA TABDIR /'*END OF','*MINIMI','*XXXXXX','*XXXXXX'/
C
      DATA TABLE  /'.TITLE ','.INPTES','.4INDEX','.INTGRL',
     &             '.WAVE F','.ANALYZ','.PROPER','.HEINZ ',
     &             '.RELCCS','.DIRRCI','.MINIMI','.XXXXXX'/
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
         CALL QUIT('End of file on LUCMD, no DIRAC input found')
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
C&&&& $INDEX transformation
            DOTRA = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  104    CONTINUE
C&&&& INTGRL: Only run integrals
            DOHRM = .TRUE.
            GO TO 100
  105    CONTINUE
C&&&& WAVE FUNCTION: call wavefunction module
            DOPSI = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  106    CONTINUE
C&&&& ANALYSIS : call anaysis module
            DOANA = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  107    CONTINUE
C&&&& PROPERTIES:  Call property module
            DOPRP = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  108    CONTINUE
C&&& HEINZ: Solve one-electronsystems using free-particle projection operators
            DOHEINZ = .TRUE.
            GO TO 100
  109    CONTINUE
C&&&& RELCCSD:  Call relccsd module
            DOCCM = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  110    CONTINUE
C&&&& DIRRCI:  Call dirrci module
            DOCIM = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
 111     CONTINUE
C&&&& MINIMIZE: Geometry optimization
            MINIMI = .TRUE.
            GOTO 100
 112     CONTINUE
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
#if defined (MOD_TRA)
      IF(DOTRA) WRITE(LUPRI,'(3X,A)') '* Transformation to Molecu'//
     &    'lar Spinor basis'
#endif
#if defined (MOD_CC)
      IF(DOCCM) WRITE(LUPRI,'(3X,A)') '* Run RELCCSD code'
#endif
#if defined (MOD_CI)
      IF(DOCIM) WRITE(LUPRI,'(3X,A)') '* Run DIRCCI code'
#endif
      IF(DOHEINZ) WRITE(LUPRI,'(3X,A)') '* One-electron systems with'//
     &    ' free particle projection and point nucleus'
      IF(INPTES)  WRITE(LUPRI,'(3X,A)') '* Input test run only '
      IF(MINIMI)  WRITE(LUPRI,'(3X,A,/,3X,A)') '* Geometry optimi'//
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
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in JOBINP.')
      ELSE
         WRITE (LUPRI,'(/,4A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in JOBINP, error in prompt.')
      END IF
    2 CONTINUE
        CALL MININP(WORD)
        GO TO 200
    3 CONTINUE
        GO TO 200
    4 CONTINUE
        GO TO 200
    1 CONTINUE
      CALL MININP(WORD)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfinp */
      SUBROUTINE DHFINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for DHF-module
C
C     Written by T.Saue - November 1992
C     Last revision: 1997/09/15 - jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (NTABLE = 30)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0)
C
      LOGICAL SET, NEWDEF, LBIT, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbdhf.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbreo.h>
C
      SAVE SET
      DATA TABLE /'.XXXXXX','.PRINT ','.MAXITR','.ERGCNV','.EVCCNV',
     &            '.FCKCNV','.NODSCF','.PHCOEF','.FIXDIF','.NODAMP',
     &            '.NODIIS','.MXDIIS','.CNVINT','.ITRINT','.INTFLG',
     &            '.TRIVEC','.TRIFCK','.XXXXXX','.EIGPRI','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.DIISTH','.DAMPFC',
     &            '.OVLSEL','.XXXXXX','.REORDE','.DIISMO','.DIISAO'/
      DATA SET/.FALSE./
C
#include <memint.h>
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
C     Initialize /CBCDHF/
C     ===================
C
      INTTYP    = '         '
      CACC      = '        '
C
C     Initialize /CBLDHF/
C     ===================
C
      ERGCNV    = .FALSE.
      EVCCNV    = .TRUE.
      FCKCNV    = .FALSE.
      SUCCES    = .FALSE.
      EXIT      = .FALSE.
      TRIVEC    = .TRUE.
      TRIFCK    = .FALSE.
      DODSCF    = .FALSE.
      PHCOEF    = .FALSE.
      FIXDIF    = .FALSE.
      IF(IDFLAG.GT.0) DODSCF    = .TRUE.
      DODAMP    = .TRUE.
      DAMPON    = .FALSE.
      DODIIS    = .TRUE.
      DIISON    = .FALSE.
      OVLSEL    = .FALSE.
      DIISAO    = .FALSE.
      DIISMO    = .TRUE.
C
C     Initialize /CBIDHF/
C     ===================
C
      MAXITR    = 50
      IPRSCF    = IPRGEN
      MXDIIS    = 10
      NBELM     = 1
      ITDIIS    = 0
      ILLINT    = 1
      ISLINT    = 1
      ISSINT    = 1
      INTBUF    = 0
      INTFLG    = 0
      ITRINT(1) = 1
      ITRINT(2) = 1
      IPREIG    = 1
      INTTYP  = 'Bare nuc.'
C
C     Initialize /CBRDHF/
C     ===================
C
      SCFCNV    = 1.0D-6
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
      CALL DZERO(SCFTID,8)
      CALL IZERO(ITRSCF,8)
C
C     Initialize /CBIREO/
C     ===================
C
      LMOORD    = .FALSE.
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
     &                     21,22,23,24,25,26,27,28,29,30), I
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
                  READ(LUCMD,*) SCFCNV
                  EVCCNV = .FALSE.
                  ERGCNV = .TRUE.
                  FCKCNV = .FALSE.
               GO TO 100
    5          CONTINUE
C&&& EVCCNV: Convergence on approximate gradient (error vector)
                  READ(LUCMD,*) SCFCNV
                  EVCCNV = .TRUE.
                  ERGCNV = .FALSE.
                  FCKCNV = .FALSE.
               GO TO 100
    6          CONTINUE
C&&& FCKCNV: Convergence on total Fock matrix
                  READ(LUCMD,*) SCFCNV
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
                  READ(LUCMD,*) ILLINT,ISLINT,ISSINT
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
               GO TO 100
   19          CONTINUE
C&&&& EIGPRI: Print control for electron and positron solutions
                  READ(LUCMD,*) IL,IS
                  IPREIG = IL+2*IS
               GO TO 100
   20          CONTINUE
               GO TO 100
   21          CONTINUE
               GO TO 100
   22          CONTINUE
               GO TO 100
   23          CONTINUE
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
                    OVLSEL = .TRUE.
               GO TO 100
   27          CONTINUE
               GO TO 100
   28          CONTINUE
C&&& LMOORD: Reorder MO orbitals (only on restart on coefficients)
C                 First read the number of orbital reorderings
C                    in each fermion corep.
C
                  LMOORD = .TRUE.
                  READ(LUCMD,*) (IREORD(I),I=1,NFSYM)
C
C                 Then read in the reordering. 
C                    Ie. >.REORDER MO's
C                        > 3 0
C                        > 1 2 2 3 3 1
C                    Means: in fermion symmetry substitue orbital
C                           1 with 2, 2 with 3 and 3 with 1
C                           Before: 1 2 3
C                           Now   : 2 3 1
C
                  DO I = 1,NFSYM
                     IF (IREORD(I).GT.MXREORD) 
     &                  CALL QUIT('You cannot reorder that many'//
     &                            ' orbitals...')
                     IF (IREORD(I).GT.0) 
     &                  READ(LUCMD,*) (IMOORD(J,I),J=1,2*IREORD(I))
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
#ifdef OPEN_MOD
                  DIISAO = .TRUE.
                  DIISMO = .FALSE.
#else
                  CALL QUIT('DIISAO NOT IN THIS DISTRIBUTION')
#endif
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
C
C     If no DIIS set MXDIIS to zero
C     =============================
C
      IF(.NOT.DODIIS) THEN
        MXDIIS = 0
      ENDIF
C
C
      INTDEF = ILLINT + 2*ISLINT + 4*ISSINT
      IF (LEVYLE) INTDEF = ILLINT
      IF(CNVINT(1).LT.DUMMY) ITRINT(1) = 1
      IF(CNVINT(2).LT.DUMMY) ITRINT(2) = 1
C
C     Print section
C     =============
C
      IF(.NOT.DODHF) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &   'DHFINP: Set-up for Dirac-Fock calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I3)') '* General print level   : ',IPRSCF
      WRITE(LUPRI,'(1X,A)') '***** TRIAL FUNCTION *****'
      IF   (TRIFCK) THEN
        WRITE(LUPRI,'(1X,A)')
     +     '* Trial two-electron Fock matrix read from file DF2FCK'
      ELSEIF(TRIVEC) THEN
        WRITE(LUPRI,'(1X,A)')
     +   '* Trial vectors read from file DFCOEF'
      ELSE
        WRITE(LUPRI,'(1X,A)')
     +   '* Trial vectors generated by the bare '//
     +   'nucleus approximation'
      ENDIF
      WRITE(LUPRI,'(1X,A)') '***** CONVERGENCE CRITERIA *****'
      WRITE(LUPRI,'(1X,A,I4)')
     +   '* Maximum number of iterations: ',MAXITR
      IF(ERGCNV) THEN
        WRITE(LUPRI,'(1X,A,A,E9.3)')
     &    '* Convergence on total energy.',
     &    ' Convergence threshold:',SCFCNV
      ELSEIF(FCKCNV) THEN
        WRITE(LUPRI,'(1X,A,A,E9.3)')
     &    '* Convergence on total Fock matrix.',
     &    ' Convergence threshold:',SCFCNV
      ELSEIF(EVCCNV) THEN
      WRITE(LUPRI,'(1X,A,A,E9.3)')
     &    '* Convergence on error vector(approximate gradient).',
     &    ' Convergence threshold:',SCFCNV
      ENDIF
      WRITE(LUPRI,'(1X,A)') '***** CONVERGENCE CONTROL *****'
      IF(DODSCF) THEN
        WRITE(LUPRI,'(1X,A)') '* Fock matrix constructed using'//
     &       ' differential density matrix'
        IF(FIXDIF) THEN
          WRITE(LUPRI,'(4X,A,F4.1)')
     &       'with fixed parameter ',D1
        ELSE
          WRITE(LUPRI,'(4X,A)')
     &       'with optimal parameter .'
        ENDIF
      ENDIF
      IF(DODIIS) THEN
         IF (NASHT .GT. 0) THEN
            DIISAO = .TRUE.
            DIISMO = .FALSE.
         ENDIF
         IF (DIISMO) WRITE(LUPRI,'(1X,A)') '* DIIS (dirac-like)'
         IF (DIISAO) WRITE(LUPRI,'(1X,A)') '* DIIS (dalton-like)'
        WRITE(LUPRI,'(1X,A,D8.2)')
     &      '* DIIS will be activated when convergence reaches : ',
     &      DIISTH
        WRITE(LUPRI,'(3X,A,I5)') '- Maximum size of B-matrix:',
     &      MXDIIS
        IF(DODAMP) THEN
          WRITE(LUPRI,'(1X,A,F5.3)')
     &     '* Damping of Fock matrix when DIIS is not activated. '//
     &     'Weight of old matrix    : ',DAMPFC
        ENDIF
      ELSEIF(DODAMP) THEN
          WRITE(LUPRI,'(1X,A,F5.3)')
     &     '* Damping of Fock matrix. Weight of old matrix    : ',
     &     DAMPFC
      ENDIF
      IF (LMOORD) THEN
          WRITE(LUPRI,'(1X,A)') '* Reorder orbitals...'
          DO I = 1,NFSYM
             WRITE(LUPRI,'(4X,A,I1)') 'Fermion corep. ',I
             DO J = 1,IREORD(I)*2,2
                WRITE(LUPRI,'(5X,A,I3,A,I3)') 'Substitute orbital ',
     &            IMOORD(J,I),' with orbital ',IMOORD(J+1,I)
             END DO
             IF (IREORD(I).EQ.0) WRITE(LUPRI,'(5X,A)') 'No reordering'
          END DO
      END IF
      IF(OVLSEL) WRITE(LUPRI,'(1X,A)')
     &    '* Vector selection based on overlap'
      WRITE(LUPRI,'(1X,A)')
     +    '* Contributions from 2-electron integrals to Fock matrix:'
      IF(LBIT(INTDEF,1)) THEN
        WRITE(LUPRI,'(3X,A)') 'LL-integrals.'
      ENDIF
      IF(LBIT(INTDEF,2)) THEN
        IF(CNVINT(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,E8.1)')
     +    'SL-integrals  below SCF convergence ',CNVINT(1)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'SL-integrals from iteration ',ITRINT(1)
        ENDIF
      ENDIF
      IF(LBIT(INTDEF,3)) THEN
        IF(CNVINT(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,E8.1)')
     +    'SS-integrals  below SCF convergence ',CNVINT(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    'SS-integrals from iteration ',ITRINT(2)
        ENDIF
      ENDIF
C
      WRITE(LUPRI,'(1X,A)') '***** OUTPUT CONTROL *****'
      IF(IPREIG.EQ.0) WRITE(LUPRI,'(1X,A)')
     +  '* No eigenvalues written out.'
      IF(IPREIG.EQ.1) WRITE(LUPRI,'(1X,A)')
     +  '* Only electron eigenvalues written out.'
      IF(IPREIG.EQ.2) WRITE(LUPRI,'(1X,A)')
     +  '* Only positron eigenvalues written out.'
      IF(IPREIG.EQ.3) WRITE(LUPRI,'(1X,A)')
     +  '* Both electron and positron component eigenvalues '
     &  //'written out.'
C
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck idnnuc */
      SUBROUTINE IDNNUC(IPRINT)
C*****************************************************************************
C
C	Count number of nuclei with identical coordinates
C
C       T.Saue, March 24 1994, University of Odense
C       Last revision : Nov 14 1994 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
#include <nuclei.h>
#include <symmet.h>
#include <ccom.h>
#include <shells.h>
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fldinp */
      SUBROUTINE FLDINP(WORD,INPERR)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER WORD*7
C
      WRITE(LUPRI,'(A)') 'FLDINP: Sorry, not yet written !'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiinp */
      SUBROUTINE PSIINP(WORK,LWORK)
C*****************************************************************************
C
C <<< Wave function Input for DIRAC >>>
C
C
C       Written by Trond Saue, Tromsoe April 1996
C       Last revision: April 27 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 8, NTABLE = 8)
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7
      DIMENSION WORK(LWORK)
C
#include <dcbgen.h>
#include <dcbpsi.h>
C
      DATA TABDIR /'*END OF','*DHFCAL','*MP2CAL','*DFTCAL',
     &             '*XXXXXX','*XXXXXX','*XXXXXX','*XXXXXX'/
      DATA TABLE  /'.DHF   ','.MP2   ','.DFT   ','.LVCORR',
     &             '.TSCORR','.XXXXXX','.XXXXXX','.XXXXXX'/
C
      NJOBS = 0
C
C     Initialize /CBLPSI/
C     ===================
C
      DODHF   = .FALSE.
      DOMP2   = .FALSE.
      DODFT   = .FALSE.
      DOLVC   = .FALSE.
      DOTSC   = .FALSE.
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
               GO TO (101,102,103,104,105,106,107,108),I
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
C&&&& DHF: Perform Dirac-Fock calculation
  101    CONTINUE
            DODHF = .TRUE.
            NJOBS = NJOBS +1
            GO TO 100
C&&&& MP2: Perform second-order Moller-Plesset calculation
  102    CONTINUE
            DOMP2 = .TRUE.
            NJOBS = NJOBS + 1
            GO TO 100
  103    CONTINUE
C&&&& DFT: Perform density functional calculation by the Panas method
            DODFT = .TRUE.
            NJOBS = NJOBS +1
            GO TO 100
  104    CONTINUE
C&&&& LVCORR: Model interatomic SS-integral contribution by classical
C             repulsion of small component atomic charges
            DOLVC = .TRUE.
            GO TO 100
  105    CONTINUE
C&&&& TSCORR: Model SS-integrals by SS nuclear attraction integrals
C             using small component atomic charges
            DOTSC = .TRUE.
            GO TO 100
  106    CONTINUE
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
         CALL QUIT('Illegal prompt in PSIINP.')
      END IF
  180 CONTINUE
C
C     Print section
C     =============
C
      CALL TITLER('Wave function module','*',126)
      IF(NJOBS.GT.0) THEN
        WRITE(LUPRI,'(1X,A)') 'Jobs in this run:'
      ENDIF
      IF(DODHF) WRITE(LUPRI,'(1X,A)') '* Dirac-Fock calculation'
      IF(DODFT) WRITE(LUPRI,'(1X,A)') '* Density functional '//
     &    'calculation (Panas method)'
#if defined (MOD_TRA)
      IF(DOMP2) WRITE(LUPRI,'(1X,A)') '* MP2 calculation'
#endif
      IF(DOLVC) WRITE(LUPRI,'(1X,A)')
     &   '* SS integrals neglected: ',
     &   '  Interatomic Coulombic SS-contributions modelled ',
     &   '  by  classical repulsion of small component atomic charges.'
      IF(DOTSC) WRITE(LUPRI,'(1X,A)')
     &   '* SS integrals neglected: ',
     &   '  Interatomic Coulombic SS-contributions modelled ',
     &   '  by  nuclear attraction integrals using ',
     &   '  small component atomic charges.'
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
#if defined (MOD_TRA)
        CALL MP2INP(WORD,WORK,LWORK)
#endif
        GO TO 200
    4 CONTINUE
C        CALL DFTINP(WORD)
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
      IF (DODHF) CALL DHFINP(WORD,.TRUE.,WORK,LWORK)
#if defined (MOD_TRA)
      IF (DOMP2) CALL MP2INP(WORD,WORK,LWORK)
#endif
C      CALL DFTINP(WORD)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ANWINP */
      SUBROUTINE ANWINP(WRDSRC,WORK,LWORK)
C*****************************************************************************
C
C <<< Analysis Input for DIRAC >>>
C
C
C       Written by Trond Saue, April 1996
C       Last revision: April 27 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0,D1=1.0D0)
      PARAMETER (NDIR = 8, NTABLE = 8)
#include <mxorb.h>
C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7,WRDSRC*7
      DIMENSION WORK(LWORK)
C
#include <dcblab.h>
#include <dcbana.h>
C
      DATA TABDIR /'*END OF','*MULPOP','*PROJEC','*PRIVEC',
     &             '*XXXXXX','*XXXXXX','*XXXXXX','*XXXXXX'/
      DATA TABLE  /'.MULPOP','.PROJEC','.PRIVEC','.XXXXXX',
     &             '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
      CALL QENTER('ANWINP')
#include <memint.h>
C
      NJOBS = 0
C
C     Initialize /CBLPRP/
C     ===================
C
      DOVEC   = .FALSE.
      DOPOP   = .FALSE.
      DOPRJ   = .FALSE.
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
            GO TO 100
  105    CONTINUE
            GO TO 100
  106    CONTINUE
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
        GO TO 200
    6 CONTINUE
        GO TO 200
    7 CONTINUE
        GO TO 200
    8 CONTINUE
        GO TO 200
C
    1 CONTINUE
      IF(DOPOP) CALL POPINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOPRJ) CALL PRJINP(WORD,.TRUE.,WORK,LWORK)
      IF(DOVEC) CALL VECINP(WORD,.TRUE.,WORK,LWORK)
C
      CALL QEXIT('ANWINP')
      RETURN
C
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck popinp */
      SUBROUTINE POPINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for Mulliken population analysis
C
C     Written by T.Saue - April 1996
C     Last revision: April 27 - 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxorb.h>
      PARAMETER (NTABLE = 8)
C
#include <dcbgen.h>
#include <dcbana.h>
#include <dcbpop.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
C
      LOGICAL SET, NEWDEF,ORBDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,LINE*72
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.AOLAB ','.LABDEF','.VECPOP','.XXXXXX',
     &            '.PRINT ','.ADDSML','.NETPOP','.XXXXXX'/
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
C     Initialize /CBIPOP/
C     ===================
C
      IPRPOP    = 0
      DONETP  = .FALSE.
      ADDSML  = .FALSE.
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
                 CALL IZERO(IPOPLAB,MAXLAB)
                 READ(LUCMD,*) NPOPLAB
                 IOFF = 1
                 DO I = 1,NPOPLAB
                   NGRP = 1
                   NOFF = MAXLAB - IOFF + 1
                   LINE = ' '
                   READ(LUCMD,'(A72)') LINE
                   POPLAB(I) = LINE(1:12)
                   CALL NUMLST(LINE(12:72),IPOPLAB(IOFF),NOFF,1,
     &                         MAXLAB,NGRP,WORK,LWORK)
                 ENDDO
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
C&&&& ADDSML: Add default group label for small components
                  ADDSML = .TRUE.
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
      IF(ADDSML) WRITE(LUPRI,'(3X,A)')
     &      '(default small component groups labels will be added)'
      WRITE(LUPRI,'(1X,A)') '* Number of spinors analyzed:'
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPOP(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),NVEC,WORK,LWORK)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjinp */
      SUBROUTINE PRJINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for module for wavefunction projection
C
C     Written by T.Saue - April 1996
C     Last revision: Apr 26 - 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxorb.h>
      PARAMETER (NTABLE = 8)
C
#include <dcbgen.h>
#include <dcbana.h>
#include <dcbprj.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
C
      LOGICAL SET, NEWDEF,ORBDEF(2), RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.PRINT ','.PROTHR','.VECPRJ','.VECREF',
     &            '.REFFIL','.XXXXXX','.XXXXXX','.XXXXXX'/
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
      ORBDEF(1) = .TRUE.
      ORBDEF(2) = .TRUE.
      REFFIL   = 'DFCOEF'
      PROTHR   = 1.0D-3
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
                  ORBDEF(1) = .FALSE.
               GO TO 100
    4          CONTINUE
C&&&& VECREF: Reference vectors
                  READ(LUCMD,'(A72)') (VECREF(I),I=1,NFSYM)
                  ORBDEF(2) = .FALSE.
               GO TO 100
    5          CONTINUE
C&&&& REFFIL: Name file of reference vectors
                 READ(LUCMD,'(A6)') REFFIL
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
     *            '" not recognized in PRJINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in PRJINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      IF(ORBDEF(1)) THEN
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
      IF(ORBDEF(2)) THEN
        DO I = 1,NFSYM
          VECREF(I) = 'all spinors'
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
      WRITE(LUPRI,'(1X,A,A6)') '* Current solution projected '//
     &       'down onto coefficients read from file ',REFFIL
      WRITE(LUPRI,'(1X,A,1P,E9.2)')
     &       '* Threshold for printing projection coefficients:',PROTHR
      WRITE(LUPRI,'(1X,A)') '* Number of spinors projected:'
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECPRJ(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)') 
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECPRJ(I)
        ENDIF
      ENDDO
      WRITE(LUPRI,'(1X,A)') '* Reference spinors :'
      DO I = 1,NFSYM
        NVEC = 0
        CALL  NUMLST(VECREF(I),IDUMMY,NFBAS(I,0),
     &              -NFBAS(I,2),NFBAS(I,1),NVEC,WORK,LWORK)
        IF(NVEC.EQ.0) THEN
          WRITE(LUPRI,'(4X,A,A3)') 
     &       '- No orbitals in fermion ircop ',FREP(I)
        ELSE
          WRITE(LUPRI,'(4X,A,A3,A,A72)') 
     &      '- Orbitals in fermion ircop ',FREP(I),' :',VECREF(I)
        ENDIF
      ENDDO
      WRITE(LUPRI,'(1X,A,I5)') '* Print level:',IPRPRJ
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vecinp */
      SUBROUTINE VECINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for vector print
C
C     Written by T.Saue - May 1996
C     Last revision: May 2 - 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (NTABLE = 8)
C
#include <dcbgen.h>
#include <dcbana.h>
#include <dcbvec.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRICMP','.AOLAB ','.VECPRI','.IPRINT',
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
               CALL QUIT('Illegal keyword in DHFINP.')
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
     &              -NFBAS(I,2),NFBAS(I,1),NVEC,WORK,LWORK)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamini */
      SUBROUTINE PAMINI(WORK,LWORK)
C*****************************************************************************
C     
C     Initialize data in DIRAC
C     Written by T.Saue May 14 1996
C     Last revision: May 14 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
#include <maxaqn.h>
#include <alphac.h>
#include <mxorb.h>
C
#include <siripc.h>
#include <huckel.h>
#include <dcbgen.h>
#include <ccom.h>
#include <cbirea.h>
#include <dcbxpr.h>
#include <dcbcls.h>
#include <dcbprl.h>
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
      DOHEINZ  = .FALSE.
      DOHRM    = .FALSE.
      DOJACO   = .FALSE.
      DOPIN    = .FALSE.
      DOPUT    = .FALSE.
      DOTRA    = .FALSE.
      DOCCM    = .FALSE.
      DOCIM    = .FALSE.
      MINIMI   = .FALSE.
      IPRGEN   = 0
      ISPHTR   = ILSPH+2*ISSPH
      IDFLAG   = ILLDIR+2*ISLDIR+4*ISSDIR
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
      IPRUSR = 0
      SEGBAS = .TRUE.
      LENBAS = 0
      NEWBAS = .TRUE.
      NEWPRP = .TRUE.
      BASDIR = 'Not implemented !!!!'            
      RELCAL = .TRUE.
      PARCAL = .FALSE.
C
C     Initialize SIPIPC from DALTON
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
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck haminp */
      SUBROUTINE HAMINP(WORK,LWORK)
C*****************************************************************************
C
C     Input section for modification of Dirac-Coulomb Hamiltonian
C
C     Written by T.Saue - May 1996
C     Last revision: May 14 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (NTABLE = 16)
      PARAMETER (D0 = 0.00D00, D1 = 1.0D0,DM2 = -2.0D0,D2 = 2.00D00,
     &           DTHRS = 1.0D-13)
C
#include <dcbgen.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
#include <pgroup.h>
#include <dcbcls.h>
#include <dcbxpr.h>
#include <dcbprl.h>
C
      LOGICAL SYMSTP,LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION KSOP(0:7),WORK(LWORK)
      CHARACTER A*80
C
      DATA TABLE /'.PRINT ','.ONESYS','.OPERAT','.GAUNT ',
     &            '.NELECT','.OPEN S','.XXXXXX','.FREEPJ',
     &            '.VEXTPJ','.URKBAL','.NOSMLV','.LEVY-L',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
      SYMSTP = .FALSE.
      CALL TITLER('Hamiltonian defined','*',127)      
C
C
C     Initialize /CBIHAM/
C     ===================
C
      IPRHAM  = 0
      ONESYS  = .FALSE.
      FREEPJ  = .FALSE.
      VEXTPJ  = .FALSE.
      URKBAL  = .FALSE.
      NOSMLV  = .FALSE.
      LEVYLE  = .FALSE.
      N1OPER  = 0
      SSMTRC  = D1
      CALL OP1INI(WORK,LWORK)
C
C     Initialize /DCBORB/
C     ===================
C
      NELEC(1) = 0
      NELEC(2) = 0
      NASH(1)  = 0
      NASH(2)  = 0
      NAELEC(1)= 0
      NAELEC(2)= 0
      CALL IZERO(NACSH,2*(1+MXOPEN))
      NOPEN    = 0
      NFRO(1)  = 0
      NFRO(2)  = 0
      CALL DZERO(DA,1+MXOPEN)
      CALL DZERO(DF,1+MXOPEN)
      CALL DZERO(DALPHA,1+MXOPEN)
C
C     Read menu file
C     ==============
C     **** Find Hamiltonian input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=300,ERR=920) WORD
      IF (WORD .EQ. '**HAMIL') THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no HAMIL input found')
  930 CONTINUE
      WORD1 = WORD
C
C     Process input for COMMON  /CBIHAM/
C     ==================================
C
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
     &                111,112,113,114,115,116),I
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
           GO TO 100
  104    CONTINUE
C&&&& GAUNT: Gaunt operator
CTROND          Must be implemented......                  
            GO TO 100
  105     CONTINUE
C&&&& Number of electrons
            READ (LUCMD,*) (NELEC(I),I=1,NFSYM)
            GO TO 100
  106     CONTINUE
C&&& ACTIVE SHELLS: read in number of active shells
C                   and their occupation
C             Read in number of open shells
C
#ifndef MOD_OPEN
               CALL QUIT('OPEN SHELL NOT AVAILABLE IN THIS VERSION')
#else
               READ(LUCMD,*) NOPEN
               IF (NOPEN.LE.0.OR.NOPEN.GT.MXOPEN) THEN
                  WRITE(LUPRI,'(1X,A,I2,A,I2)') 
     &               'HAMINP: NOPEN must be between ',
     &               1,' and ',MXOPEN
                  CALL QUIT('HAMINP ERROR: Open shell error!')
               END IF
C
C              Read in number of active shells and their occ.
C
               DO IOPEN = 1,NOPEN
                  READ(LUCMD,*) (NACSH(J,IOPEN),J=1,NFSYM)
                  READ(LUCMD,'(A)') A
                  CALL FRAC2I(DF(IOPEN),A)
                  IF (DF(IOPEN).GT.(D1+DTHR)) THEN
                     WRITE(LUPRI,'(A,I2,A,F6.4,A)')
     &                  'HAMINP: Frac. occ. for shell ',IOPEN,
     &                  ' greater than 1 (',DF(IOPEN),')'
                     CALL QUIT('HAMINP ERROR: Open shell error!')
                  END IF
               END DO
#endif
            GO TO 100
  107     CONTINUE
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
             NOSMLV = .TRUE.
             SSMTRC = D0
             IF (NPRPS.GT.4) 
     &          CALL QUIT ('give .LEVY-LE option before .OPERAT')
             ICHANG = ICHANG + 1
             GO TO 100
  113      CONTINUE
               GOTO 100
  114      CONTINUE
              GOTO 100
  115      CONTINUE
              GOTO 100
  116      CONTINUE
              GOTO 100
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
C     Determine occupation
C     ====================
C
      DF(0) = D1
      DA(0) = D1
      DALPHA(0) = D0
      DO IFRP = 1,NFSYM
         NISH(IFRP)  = NELEC(IFRP)
         DTEMP = D0
         DO IOPEN = 1,NOPEN
            NASH(IFRP)   = NASH(IFRP) + NACSH(IFRP,IOPEN)
            DTEMP        = DTEMP +
     &                     NACSH(IFRP,IOPEN)*D2*DF(IOPEN)
         END DO
C         IF (ABS(MOD(DTEMP,D1)).GE.1.0D-10) THEN
C            WRITE(LUPRI,'(A,F,A)') 
C     &         'ERROR! Number of active electrons not an integer (',
C     &         DTEMP,')'
C            CALL QUIT('HAMINP ERROR: Open shell error!')
C         END IF
         NAELEC(IFRP) = NINT(DTEMP)
         IF (NELEC(IFRP).LT.0) 
     &      CALL QUIT('HAMINP: negative number of '//
     &                'closed shell electrons!')
         IF(MOD(NELEC(IFRP),2).EQ.1) GOTO 1000
         NISH(IFRP)  = NELEC(IFRP)/2
         NOCC(IFRP)  = NISH(IFRP) + NASH(IFRP)
      END DO
      DO IOPEN = 1,NOPEN
C
C        We do     allow fractional occupation of 1
C           This corresponds to a inactive electron
C           and gives DALPHA = 0
C
         IF (ABS(DF(IOPEN)-D1).LE.DTHR) THEN
            WRITE(LUPRI,'(X,A,I2,A,F6.4)')  
     &         'WARNING: Open shell no. ',IOPEN,
     &         ' has frac. occ. of ',DF(IOPEN)
            WRITE(LUPRI,'(X,A)')
     &         'This corresponds to an inactive shell. So ALPHA and '//
     &         'A is reset.'
            DA(IOPEN)     = D1
            DF(IOPEN)     = D1
            DALPHA(IOPEN) = D0
         ELSE
            DM        = D2*(NACSH(1,IOPEN)+NACSH(2,IOPEN))
            DN        = DM*DF(IOPEN)
            DA(IOPEN) = (DM*(DN-D1))/(DN*(DM-D1)) 
            DALPHA(IOPEN) = (D1-DA(IOPEN))/(D1-DF(IOPEN))
         END IF
      ENDDO
      NASHT = NASH(1) + NASH(2)
      NFROT = NFRO(1) + NFRO(2)
      NELECT  = NELEC(1) + NELEC(2)
C
C     Remove inconsitencies
C     =====================
C
      IF (LEVYLE) THEN
         IF(FREEPJ) CALL QUIT('Conflicting options : LEVYLE & FREEPJ !')
         IF(VEXTPJ) CALL QUIT('Conflicting options : LEVYLE & VEXTPJ !')
      ENDIF
C
      IF(FREEPJ.OR.VEXTPJ.OR.LEVYLE) URKBAL = .FALSE.
      IF(FREEPJ.AND.VEXTPJ) THEN
        WRITE(LUPRI,'(A)') '***** HAMINP ERROR *****'
        WRITE(LUPRI,'(A)') 'You cannot specify both VEXTPJ and FREEPJ !'
        CALL QUIT('HAMINP: Both VEXTPJ and FREEPJ !')
      ENDIF
C
C     Print section
C     =============
C
 400  CONTINUE
      WRITE(LUPRI,'(1X,A)') '* Number of inactive electrons:'
      WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion corep ',I,':',NELEC(I),I=1,NFSYM)
      WRITE(LUPRI,'(1X,A)') '* Number of active electrons:'
      WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion corep ',I,':',NAELEC(I),I=1,NFSYM)
      WRITE(LUPRI,'(1X,A)') '* Total number of electrons:'
      WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion corep ',I,':',NAELEC(I)+NELEC(I),I=1,NFSYM)
      IF (NFSYM.GT.1) THEN
         ISUM = 0
         DO I = 1,NFSYM
            ISUM = ISUM + NAELEC(I)+NELEC(I)
         END DO
         WRITE(LUPRI,'(3X,A,3X,I3)')
     &      'Sum of coreps   :',ISUM
      END IF
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(1X,A)') '* Number of inactive shells:'
      WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion irrep ',I,':',NISH(I),I=1,NFSYM)
      WRITE(LUPRI,'(1X,A,I2)')
     &   '* Number of open shells : ',NOPEN
      DO IOPEN = 1,NOPEN
         WRITE(LUPRI,'(3X,A,I2)') 'Open shell number ',IOPEN
         WRITE(LUPRI,'(5X,A)') 'Number of active shells:'
         WRITE(LUPRI,'(5X,A,I2,A,3X,I3)')
     &      ('Fermion irrep ',I,':',NACSH(I,IOPEN),I=1,NFSYM)
         WRITE(LUPRI,'(7X,A,F6.4,A/7X,A,F6.4,A/7X,A,F6.4,A)')
     &       '(Fractional occupation :      ',DF(IOPEN),',',
     &       ' Coupling coefficient ALPHA : ',DALPHA(IOPEN),',',
     &       ' Coupling coefficient A     : ',DA(IOPEN),')'
      END DO
      IF(NFROT.GT.0) THEN
        WRITE(LUPRI,'(1X,A)') '* Number of frozen shells:'
        WRITE(LUPRI,'(3X,A,I2,A,3X,I3)')
     &   ('Fermion irrep ',I,':',NFRO(I),I=1,NFSYM)
      ENDIF
      WRITE(LUPRI,'(1X,A,I5)') '* Print level :',IPRHAM
C
      IF(URKBAL) THEN
        WRITE(LUPRI,'(1X,A)') '* Unrestricted kinetic balance'
      ELSE
        WRITE(LUPRI,'(1X,A)') '* Restricted kinetic balance'
      ENDIF
C
C     Projections on basis set
C
      IF(ISPHTR.EQ.3) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(1X,A)')
     &      '* Default RKB projection:',
     &      '  1: Pre-projection on scalar basis',
     &      '  2: Removing unphysical solutions for no external field'
         ELSE
            WRITE(LUPRI,'(1X,A)')
     &           '*** WARNING ****',
     &           ' RKB pre-projection without removing unphysical '//
     &           'solutions for no external field',
     &           ' I HOPE YOU NOW WHAT YOU ARE DOING!'
         ENDIF
      ELSEIF(ISPHTR.EQ.1) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(1X,A)')
     &      '* RKB projection:',
     &      '  1: No Pre-projection on scalar basis',
     &      '  2: Removing unphysical solutions for no external field',
     &      '*** WARNING ****, Ok, but may give severe linear',
     &           '   dependence problems for large basis sets'
         ELSE
            WRITE(LUPRI,'(1X,A)')
     &           '*** WARNING ****',
     &           ' Spherical large component basis and Cartesian'//
     &           ' small components',
     &           ' I HOPE YOU NOW WHAT YOU ARE DOING!'
         ENDIF
      ELSEIF(ISPHTR.EQ.0) THEN
         IF(.NOT.URKBAL) THEN
            WRITE(LUPRI,'(1X,A)')
     &      '* RKB projection:',
     &      '  1: Cartesian large components',
     &      '  2: Removing unphysical solutions for no external field'
         ELSE
            WRITE(LUPRI,'(1X,A)')
     &           '* Unresticted Cartesian basis'
         ENDIF
      ELSE
         WRITE(LUPRI,'(1X,A,I5)')
     &        'ISPHTR not legal:',ISPHTR
         CALL QUIT('HAMINP: Wrong ISPHTR!')
      ENDIF
C
      IF(.NOT.LEVYLE) THEN
        WRITE(LUPRI,'(1X,A)') '* Dirac-Coulomb Hamiltonian'
      ELSE
        WRITE(LUPRI,'(1X,A)') '* Levy-Leblond Hamiltonian'
      ENDIF
      IF(ICHANG.GT.0) WRITE(LUPRI,'(1X,A)') 
     &    'with the following modifications:'
      IF(N1OPER.GT.0) THEN
        WRITE(LUPRI,'(1X,A,A)') '* The following one-electron ',
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
        WRITE(LUPRI,'(1X,A,A)') '* Two-electron part ignored ',
     &   '(one-electron system)'
      ENDIF
      IF (.NOT.LEVYLE) THEN
         IF(FREEPJ) THEN
           WRITE(LUPRI,'(1X,A)')
     &        '* Free positronic solutions projected out'
         ENDIF
         IF(VEXTPJ) THEN
           WRITE(LUPRI,'(1X,A)') '* External potential positronic '//
     &          'solutions projected out'
         ENDIF
         IF(NOSMLV) WRITE(LUPRI,'(1X,A)')
     &      '* No nuclear attraction integrals for small components !'
      ENDIF
C
C
C     Calculate LV-correction
C     (SS integrals modelled by Coulombic correction)
C     ===============================================
C
      CALL LVCORR(WORK,LWORK,IPRHAM)
C
      IF(SYMSTP) CALL QUIT('HAMINP: Non-symmetric operator !')
C
C     Generating Transformation matrices for 
C     spherical basis and scalar RKB type basis
C     
      IF(LBIT(ISPHTR,1)) CALL SPHGEN(WORK,LWORK,IPRHAM)
      IF(LBIT(ISPHTR,2)) CALL RKBGEN(WORK,LWORK,IPRHAM)
C
      RETURN
C
 1000 CONTINUE
      CALL QUIT('OCCINP:Inconsistency between electrons '//
     &          'and active shells!')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1ini */
      SUBROUTINE OP1INI(WORK,LWORK)
C*****************************************************************************
C
C     Initialize one-electron integrals for standard Dirac-Coulomb Hamiltonian
C
C     Written by T.Saue May 21 1996
C     Last reviosn : may 21 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0,DM2 = -2.0D0)
C
#include <dcbgen.h>
#include <dcbham.h>
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <cbihr2.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL LBIT,TEST
C
      CALL PAOVEC(WORK,LWORK,0,IPRTWO)
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &   'Maximum memory load for two-integral processing in HERMIT:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') '* LL-integrals:'
      TEST = LBIT(IDFLAG,1)
      CALL MEMLOP(1,TEST)  
      WRITE(LUPRI,'(1X,A)') '* SL-integrals:'
      TEST = LBIT(IDFLAG,2)
      CALL MEMLOP(2,TEST)  
      WRITE(LUPRI,'(1X,A)') '* SS-integrals:'
      TEST = LBIT(IDFLAG,3)
      CALL MEMLOP(3,TEST)  
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MINANA */
      SUBROUTINE MINANA(WORK,LWORK)
C***********************************************************************
C
C     Readin analasis input for the final geometry
C
C     Written by J. Thyssen, Nov 13 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('MINANA')
C
      IF(DOANA) THEN
        OPEN(LUCMD,FILE = 'DIRAC.INP')
        CALL ANWINP('**ANA F',WORK,LWORK)
        CLOSE(LUCMD)
      ENDIF
C
      CALL QEXIT('MINANA')
C
      RETURN
      END
