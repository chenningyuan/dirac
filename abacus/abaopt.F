C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id$
C
#ifdef UNDEF
/* Comdeck log */
950915-vebjornb: New module taking care of optimizations.
                 It uses both the old cbiwlk.h and the new optinf.h
#endif
C  /* Deck mininp */
      SUBROUTINE MININP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <optinf.h>
#include <cbiwlk.h>
Cjth
#ifdef PRG_DIRAC
#include <dcbgen.h>
#include <dcbgrd.h>
#else
#include <gnrinf.h>
#include <infinp.h>
#endif
      PARAMETER (NTABLE = 48, ITRMX = 25, MXREJ = 5, D0 = 0.0D0,
     &                              DEFTHR = 1.0D-6, DEFDSP = 1.0D-3)
      LOGICAL NEWDEF, TRSTCH, FIRST
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, STR*3
      DIMENSION II(3)
C
      SAVE FIRST
      DATA TABLE /'.PRINT ', '.MAX IT', '.TRUSTR', '.TR FAC',
     &            '.TR LIM', '.MAX RE', '.NOTRUS', '.ENERGY',
     &            '.GRADIE', '.STEP T', '.CONDIT', '.NOBREA',
     &            '.SP BAS', '.PREOPT', '.VISUAL', '.VRML  ',
     &            '.SYMTHR', '.INTERA', '.VR-BON', '.VR-EIG',
     &            '.INITHE', '.INITEV', '.HESFIL', '.REJINI', 
     &            '.STEEPD', '.RANKON', '.PSB   ', '.DFP   ',
     &            '.BFGS  ', '.NEWTON', '.QUADSD', '.SCHLEG',
     &            '.NOREIN', '.BAKER ', '.M-BFGS', '.CARTES',
     &            '.REDINT', '.INIRED', '.1STORD', '.2NDORD',
     &            '.GRDINI', '.DISPLA', '.DOTRCK', '.TWOGRD',
     &            '.IGTHRS', '.XXXXXX', '.XXXXXX', '.XXXXXX'/
C
      DATA FIRST /.TRUE./
C
      IF (.NOT. FIRST) THEN
         IF ((WORD .NE. '*END OF') .AND. (WORD(1:2) .NE. '**')) THEN
 969        READ(LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GOTO 969
         END IF
         RETURN
      END IF
      FIRST = .FALSE.
C
      NEWDEF = (WORD .EQ. '*MINIMI')
C
      ICHANG = 0
      IPDEF  = 0
      NOTRST = .FALSE.
      NOBRKS = .FALSE.
      BRKSYM = .FALSE.
      ITRBRK = -1
      NWSYMM = .FALSE.
      NEWSYM = .FALSE.
      DOSPE  = .FALSE.
      DOPRE  = .FALSE.
      FINPRE = .FALSE.
      KEEPHE = .FALSE.
      REINIH = .TRUE.
      REJINI = .FALSE.
      GRDINI = .FALSE.
      CHGRDT = .FALSE.
      NUMPRE = 0
      IPRE   = 0
      STEEPD = .FALSE.
      RANKON = .FALSE.
      PSB    = .FALSE.
      DFP    = .FALSE.
      BFGS   = .FALSE.
      MBFGS  = .FALSE.
      SCHLEG = .FALSE.
      NEWTON = .FALSE.
      QUADSD = .FALSE.
#if defined(PRG_DIRAC)
      FSTORD = .TRUE.
      SNDORD = .FALSE.
C
C     Initialize dcbgrd.h
C
      CALL INIGRD
      DOTRCK = .TRUE.
#else
      FSTORD = .FALSE.
      SNDORD = .TRUE.
#endif
      REDINT = .FALSE.
      CARTCO = .FALSE.
      INRDHS = .FALSE.
      INITHS = .FALSE.
      HSFILE = .FALSE.
      EVLINI = -1.0D0
      VISUAL = .FALSE.
      VRML   = .FALSE.
      VRBOND = .FALSE.
      VREIGV = .FALSE.
      INTRAC = .FALSE.
      TRSTCH = .FALSE.
      IPRINT = 0
      TRSTRA = 0.5D0
      TRSTIN = 1.2D0
      TRSTDE = 0.7D0
      RTENBD = 0.4D0
      RTENGD = 0.8D0
      RTRJMN = 0.1D0
      RTRJMX = 3.0D0
      GECONV = .FALSE.
      GETOLL = DEFTHR
      BAKER  = .FALSE.
      ITRNMR = 0
      ITRMAX = ITRMX
      MAXREJ = MXREJ
      DISPLA = DEFDSP
      GRDTHR = DEFTHR
      THRSTP = DEFTHR
      THRERG = DEFTHR
      THRSYM = SQRT(DEFTHR)
      ICONDI = 2
      NCRTOT = 0
      NCART  = 0
      ENERGY = D0
      ERGOLD = D0
      ERGPRD = D0
      ERGPRO = D0
      STPNRM = D0
      STPNRO = D0
      GRADNM = D0
      CALL DZERO(STPDIA,2*MXCOOR)
      CALL DZERO(STPSYM,2*MXCOOR)
      CALL DZERO(STPINT,2*MXCOOR)
      CALL DZERO(GRDDIA,2*MXCOOR)
      CALL DZERO(EVAL  ,2*MXCOOR)
      CALL DZERO(EVALOL,2*MXCOOR)
      DO 123 I = 0, 7
         CNDHES(I) = D0
         INDHES(I) = 0
 123  CONTINUE
C
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
         READ (LUCMD, '(A7)') WORD
         PROMPT = WORD(1:1)
         IF ((PROMPT .EQ. '!') .OR. (PROMPT .EQ. '#')) THEN
            GOTO 100
         ELSE IF (PROMPT .EQ. '.') THEN
            ICHANG = ICHANG + 1
            DO 200 I = 1, NTABLE
               IF (TABLE(I) .EQ. WORD) THEN
                  GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                  11,12,13,14,15,16,17,18,19,20,
     &                  21,22,23,24,25,26,27,28,29,30,
     &                  31,32,33,34,35,36,37,38,39,40,
     &                  41,42,43,44,45,46,47,48), I
               END IF
 200        CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GOTO 100
            END IF
            WRITE(LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &                         '" not recognized in MININP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' INPUT KEYWORDS',LUPRI)
            CALL QUIT('Illegal keyword in MININP.')
 1          READ(LUCMD,*) IPRINT
            IPRGRD = IPRINT
            GOTO 100
 2          READ(LUCMD,*) ITRMAX
            GOTO 100
 3          READ(LUCMD,*) TRSTRA
            TRSTCH = .TRUE.
            GOTO 100
 4          READ(LUCMD,*) TRSTIN, TRSTDE
            TRSTCH = .TRUE.
            GOTO 100
 5          READ(LUCMD,*) RTENBD, RTENGD, RTRJMN, RTRJMX
            TRSTCH = .TRUE.
            GOTO 100
 6          READ(LUCMD,*) MAXREJ
            GOTO 100
 7          NOTRST = .TRUE.
            GOTO 100
 8          READ(LUCMD,*) THRERG
            CHGRDT = .TRUE.
            GOTO 100
 9          READ(LUCMD,*) GRDTHR
            CHGRDT = .TRUE.
            GOTO 100
 10         READ(LUCMD,*) THRSTP
            CHGRDT = .TRUE.
            GOTO 100
 11         READ(LUCMD,*) ICONDI
            GOTO 100
 12         NOBRKS = .TRUE.
            GOTO 100
 13         READ(LUCMD,*) SPBSTX
            DOSPE = .TRUE.
            GOTO 100
 14         READ(LUCMD,*) NUMPRE
            IF ((NUMPRE .LT. 1) .OR. (NUMPRE .GT. MAXPRE)) THEN
               WRITE(LUPRI,'(/,A,I2,A/)')
     &        ' Number of preoptimization sets must be between 1 and',
     &            MAXPRE,'!'
               CALL QUIT('Illegal number of preoptimization sets.')
            ELSE
               DO 144 I = 1, NUMPRE 
                  READ(LUCMD,*) PREBTX(I)
 144           CONTINUE
               DOPRE = .TRUE.
            END IF
            GOTO 100
 15         VISUAL = .TRUE.
            GOTO 100
 16         VRML = .TRUE.
            GOTO 100
 17         READ(LUCMD,*) THRSYM
            GOTO 100
 18         INTRAC = .TRUE.
            GOTO 100
 19         VRBOND = .TRUE.
            GOTO 100
 20         VREIGV = .TRUE.
            GOTO 100
 21         INITHS = .TRUE.
            GOTO 100
 22         READ(LUCMD,*) EVLINI
            GOTO 100
 23         HSFILE = .TRUE.
            GOTO 100
 24         REJINI = .TRUE.
            GOTO 100
 25         STEEPD = .TRUE.
            GOTO 100
 26         RANKON = .TRUE.
            GOTO 100
 27         PSB = .TRUE.
            GOTO 100
 28         DFP    = .TRUE.
            GOTO 100
 29         BFGS   = .TRUE.
            GOTO 100
 30         NEWTON = .TRUE.
            GOTO 100
 31         QUADSD = .TRUE.
            GOTO 100
 32         SCHLEG = .TRUE.
            GOTO 100
 33         REINIH = .FALSE.
            GOTO 100
 34         BAKER  = .TRUE.
            GOTO 100
 35         MBFGS  = .TRUE.
            GOTO 100
 36         CARTCO = .TRUE.
            GOTO 100
 37         REDINT = .TRUE.
            GOTO 100
 38         INRDHS = .TRUE.
            GOTO 100
 39         FSTORD = .TRUE.
            SNDORD = .FALSE.
            GOTO 100
 40         SNDORD = .TRUE.
            GOTO 100
 41         GRDINI = .TRUE.
            GOTO 100
 42         READ (LUCMD,*) DISPLA
Cjth - removed ??            GOTO 100
            GOTO 100
 43         READ (LUCMD,*) ITEMP
            IF ( ITEMP .GT. 0) THEN
               DOTRCK = .TRUE.
            ELSE
               DOTRCK = .FALSE.
            END IF
            GOTO 100
 44         CONTINUE
            READ(LUCMD,*) II(1),II(2),II(3)
            DO III=1,3
               IF ( II(III) .GT. 0 ) THEN
                  DO2GRD(III) = .TRUE.
               ELSE
                  DO2GRD(III) = .FALSE.
               END IF
             END DO
             GOTO 100
 45         CONTINUE
            READ(LUCMD,*) IGTHRS
            GOTO 100
 46         CONTINUE
            GOTO 100
 47         CONTINUE
            GOTO 100
 48         CONTINUE
            GOTO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GOTO 300
         ELSE
            WRITE(LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &           '" not recognized in MININP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in MININP.')
         END IF
      END IF
 300  CONTINUE
C
      IF ((.NOT. MINIMI) .AND. (ICHANG .GT. 0)) THEN
         WRITE(LUPRI,'(//A/)')
     &          ' Minimization not requested, *MINIMI input ignored.'
         RETURN
      END IF
C
      IF (.NOT. MINIMI) RETURN
      CALL HEADER('Chosen parameters for MINIMI:',0)
C
C     Check for visualization
C
      IF (VISUAL) THEN
         WRITE(LUPRI,'(A)') ' Visualization has been ' //
     &        'requested. No minimization will be done.'
         WRITE(LUPRI,'(A)')
     &        ' VRML-file of geometry will be created.'
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' Bonds will be drawn between nearby atoms.'
         IF (VREIGV) THEN
            WRITE(LUPRI,'(A)') ' Eigenvectors can ' //
     &           'only be visualized during a minimization.'
            WRITE(LUPRI,'(A)') ' Keyword will be ignored!'
            VREIGV = .FALSE.
         END IF
         WRITE(LUPRI,'(A)') ' Other keywords will be ignored!'
         RETURN
      END IF
C     
C     The type of minimization is determined, BFGS is default.
C
C     For the relativistic case,  2nd order methods are .NOT.
C     implemented, so if 2nd order is requested then we quit
C
C     
      ITYP = 0
      IF (STEEPD) THEN
         WRITE(LUPRI,'(A)')
     &        ' 1st order steepest descent method will be used.'
         ITYP = ITYP + 1
      END IF
      IF (RANKON) THEN
         WRITE(LUPRI,'(A)')
     &        ' 1st order method with rank one update will be used.'
         ITYP = ITYP + 1
      END IF
      IF (PSB) THEN
         WRITE(LUPRI,'(A)')
     &        ' 1st order method with PSB update will be used.'
         ITYP = ITYP + 1
      END IF
      IF (DFP) THEN
         WRITE(LUPRI,'(A)')
     &        ' 1st order method with DFP update will be used.'
         ITYP = ITYP + 1
      END IF
      IF (BFGS) THEN
         WRITE(LUPRI,'(A)')
     &        ' 1st order method with BFGS update will be used.'
         ITYP = ITYP + 1
      END IF
      IF (MBFGS) THEN
         WRITE(LUPRI,'(A)') ' 1st order method with ' //
     &        '"multiple BFGS" updating scheme will be used.'
         ITYP = ITYP + 1
      END IF
      IF (SCHLEG) THEN
         WRITE(LUPRI,'(A)')
     &  ' 1st order method with Schlegels updating scheme will be used.'
         ITYP = ITYP + 1
      END IF
      IF (NEWTON) THEN
         WRITE(LUPRI,'(A)')
     &        ' 2nd order level-shifted Newton method will be used.'
         ITYP = ITYP + 1
         CALL QUIT('SECOND ORDER METHODS NOT IMPLEMENTED IN DIRAC')
      END IF
      IF (QUADSD) THEN
         WRITE(LUPRI,'(A)')
     &        ' 2nd order quadratic steepest ' //
     &        'descent method will be used.'
         ITYP = ITYP + 1
         CALL QUIT('SECOND ORDER METHODS NOT IMPLEMENTED IN DIRAC')
      END IF
      IF (FSTORD) THEN
         WRITE(LUPRI,'(A)')
     &        ' Default 1st order method will be used:   BFGS update.'
         BFGS = .TRUE.
CKR         REJINI = .TRUE.
         ITYP = ITYP + 1
      END IF
      IF (SNDORD) THEN
         WRITE(LUPRI,'(A)')
     &        ' Default 2nd order method will be used:   Newton method.'
         NEWTON = .TRUE.
         ITYP = ITYP + 1
         CALL QUIT('SECOND ORDER METHODS NOT IMPLEMENTED IN DIRAC')
      END IF
      IF (ITYP .EQ. 0) THEN
         WRITE(LUPRI,'(A)')
     &        ' Using default optimization method:   BFGS update.'
         FSTORD = .TRUE.
         BFGS = .TRUE.
      ELSE IF (ITYP .GT. 1) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)') ' WARNING! More than one ' // 
     &        'minimization method has been selected.'
         WRITE(LUPRI,'(A)') ' Aborting.'            
         WRITE(LUPRI,*)
         CALL QUIT
     &        ('More than one minimization method chosen in MININP.')
      END IF
C     
      IF (REDINT .AND. CARTCO) THEN
         WRITE(LUPRI,'(A)') ' Both Cartesian coordinates and ' //
     &        'redundant internal coordinates has been requested.'
         WRITE(LUPRI,'(A)') ' Keyword .REDINT is ignored.'
         WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &        'in Cartesian coordinates.'
         REDINT = .FALSE.
      ELSE IF (REDINT) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &        'in redundant internal coordinates.'
      ELSE IF (CARTCO) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &        'in Cartesian coordinates.'
      ELSE
         IF (FSTORD) THEN
            REDINT = .TRUE.
            WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &           'in redundant internal coordinates.'
         ELSE
            CARTCO = .TRUE.
            WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &           'in Cartesian coordinates.'
         END IF
      END IF
C     
      IF (HSFILE) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .HESFIL only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Initial Hessian will be read from file.'
         END IF
      END IF
      IF (INITHS) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITHE only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (HSFILE) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITHE has no effect when .HESFIL ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            INITHS = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Initial Hessian will be calculated.'
         END IF
      END IF
      IF (INRDHS) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (HSFILE) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED has no effect when .HESFIL ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (INITHS) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED has no effect when .INITHE ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            INRDHS = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &      ' Initial Hessian will be diagonal in internal coordinates.'
         END IF
      END IF
      IF (EVLINI .GT. 0.0D0) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (HSFILE) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV has no effect when .HESFIL ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (INITHS) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV has no effect when initial ' //
     &           'Hessian is calculated.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE
            WRITE(LUPRI,'(A,F10.6)') ' Initial diagonal ' //
     &           'Hessian will have elements equal to: ', EVLINI
         END IF
      ELSE IF (.NOT. (REDINT .OR. INRDHS)) THEN
         EVLINI = 1.0D0
      END IF
      WRITE(LUPRI,*)
      IF (.NOT. REINIH) THEN
         WRITE(LUPRI,'(A)') ' Hessian will NOT be reinitialized' //
     &        ' when a non-zero index is encountered.'
         WRITE(LUPRI,*)
      END IF
      IF (REJINI) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .REJINI only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(A)') ' Hessian will be reinitialized' //
     &           ' after rejected steps.'
            WRITE(LUPRI,*)
         END IF
      END IF
      IF (GRDINI) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .GRDINI only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(A)') ' Hessian will be reinitialized' //
     &           ' when the norm of the gradient increases.'
            WRITE(LUPRI,*)
         END IF
      END IF
C     
      IF (IPRINT .NE. IPDEF) THEN
         WRITE(LUPRI,'(A,I10)') ' Print level in MINIMI  :',IPRINT
      END IF
      IF (ITRMAX .NE. ITRMX) THEN
         WRITE(LUPRI,'(A,I10)') ' Maximum # of iterations:',ITRMAX
      END IF
      IF (MAXREJ .NE. MXREJ) THEN
         WRITE(LUPRI,'(A,I10)') ' Max # of rejected steps:',MAXREJ
      END IF
      IF (VRML) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)') ' VRML-files of initial and final ' //
     &        'geometry will be created.'
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' Bonds will be drawn between nearby atoms.'
         IF (VREIGV) WRITE(LUPRI,'(A)')
     &        ' Eigenvectors will be visualized.'
      ELSE 
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' .VR-BOND only has effect when .VRML is specified.' // 
     &        ' Keyword ignored.'
         IF (VREIGV) WRITE(LUPRI,'(A)')
     &        ' .VR-EIG only has effect when .VRML is specified.' // 
     &        ' Keyword ignored.'
      END IF
      IF (INTRAC) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)')
     &        ' Interactive run: Information for each iteration is',
     &        '                  written to unit LUERR.'
      END IF
      IF (DOPRE) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)') ' Preoptimization will be' //
     &        ' performed with the basis set(s):'
         DO 244 I = 1, NUMPRE
            WRITE(LUPRI,'(A,A60)') '      ',PREBTX(I)
 244     CONTINUE
      END IF
      IF (DOSPE) THEN
         WRITE(LUPRI,'(A)') ' Single point energy will be' //
     &        ' calculated using the basis:'
         WRITE(LUPRI,'(A,A60)') '      ',SPBSTX
      END IF
      IF (TRSTCH) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(/A/A/,(A,F10.4))')
     &        ' Restricted step control parameters',
     &        ' ----------------------------------',
     &        ' Initial trust radius   :',TRSTRA,
     &        ' Trust radius increment :',TRSTIN,
     &        ' Trust radius decrement :',TRSTDE,
     &        ' Bad prediction ratio   :',RTENBD,
     &        ' Good prediction ratio  :',RTENGD,
     &        ' Rejection ratio, low   :',RTRJMN,
     &        ' Rejection ratio, high  :',RTRJMX
      END IF
      WRITE(LUPRI,*)
      IF (NOTRST) THEN
         WRITE(LUPRI,*) 'No trust region will be used for steps.'
         WRITE(LUPRI,*)
      END IF
      IF (NOBRKS) THEN
         WRITE(LUPRI,*)
     &        'Symmetry will not be broken during minimization.'
         WRITE(LUPRI,*)
      END IF
      IF (THRSYM .NE. SQRT(DEFTHR)) THEN
         WRITE(LUPRI,'(A,D13.2)')
     &        ' Threshold for symmetry-break set to       : ', THRSYM
         IF ((THRSYM .LT. D0) .OR. (THRSYM .GT. 0.1D0)) THEN
            WRITE(LUPRI,*) 'Threshold negative or larger than 0.1'
            WRITE(LUPRI,'(A,D13.2)')
     &           ' Threshold reset to:', SQRT(DEFTHR)
            THRSYM = SQRT(DEFTHR)
         END IF
      END IF
      IF (GRDTHR .NE. DEFTHR) THEN
         WRITE(LUPRI,'(A,D13.2)')
     &        ' Convergence threshold for gradient set to : ', GRDTHR
         IF ((GRDTHR .LT. D0) .OR. (GRDTHR .GT. 0.1D0)) THEN
            WRITE(LUPRI,*) 'Threshold negative or larger than 0.1'
            WRITE(LUPRI,'(A,D13.2)')
     &           ' Threshold reset to:', DEFTHR
            GRDTHR = DEFTHR
         END IF
      END IF
      IF (THRERG .NE. DEFTHR) THEN
         WRITE(LUPRI,'(A,D13.2)')
     &        ' Convergence threshold for energy set to   : ', THRERG
         IF ((THRERG .LT. D0) .OR. (THRERG .GT. 0.1D0)) THEN
            WRITE(LUPRI,*) 'Threshold negative or larger than 0.1'
            WRITE(LUPRI,'(A,D13.2)')
     &           ' Threshold reset to:', DEFTHR
            THRERG = DEFTHR
         END IF
      END IF
      IF (THRSTP .NE. DEFTHR) THEN
         WRITE(LUPRI,'(A,D13.2)')
     &        ' Convergence threshold for step set to     : ', THRSTP
         IF ((THRSTP .LT. D0) .OR. (THRSTP .GT. 0.1D0)) THEN
            WRITE(LUPRI,*) 'Threshold negative or larger than 0.1'
            WRITE(LUPRI,'(A,D13.2)')
     &           ' Threshold reset to:', DEFTHR
            THRSTP = DEFTHR
         END IF
      END IF
      IF (DISPLA .NE. DEFDSP) THEN
         WRITE (LUPRI,'(A,D13.2)')
     &        ' In case of numerical gradients, displacements are : ',
     &        DISPLA
      END IF
      IF (ICONDI .NE. 2) THEN
         WRITE(LUPRI,'(A,I3)')
     &        ' Number of convergence criteria set to     : ', ICONDI
         IF ((ICONDI .LT. 1) .OR. (ICONDI .GT. 3)) THEN
            WRITE(LUPRI,*) 'Acceptable values are 1, 2 and 3.'
            ICONDI = 2
            WRITE(LUPRI,'(A,I3)')
     &           ' Number reset to:', ICONDI
         END IF
      END IF
      IF (BAKER) THEN
         WRITE(LUPRI,'(A)')
     &        ' Baker''s convergence criteria will be used.'
         WRITE(LUPRI,*)
      END IF
#ifdef PRG_DIRAC
      CALL CYESNO(STR,DO2GRD(1))
      WRITE(LUPRI,'(1X,A,A)') 'Calculate LL-two-electron gradient: ',
     &                           STR
      CALL CYESNO(STR,DO2GRD(2))
      WRITE(LUPRI,'(1X,A,A)') 'Calculate LS-two-electron gradient: ',
     &                           STR
      CALL CYESNO(STR,DO2GRD(3))
      WRITE(LUPRI,'(1X,A,A)') 'Calculate SS-two-electron gradient: ',
     &                           STR
      CALL CYESNO(STR,DOTRCK)
      WRITE(LUPRI,'(1X,A,A)') 'Skip LS/SS-two-el. grd. if small  : ',
     &                           STR
      DTMP = 10.0D00**(-IGTHRS)
      WRITE(LUPRI,'(1X,A,1P,E8.2)')
     &   'Screening threshold in gradient calculation ',DTMP
      IF (IGTHRS .GE. 16) THEN
         WRITE(LUPRI,'(1X,A,1P,E8.2)')
     &      'WARNING: Screening turned off...'
         IGTHRS    = 0
      END IF
#endif
C
      RETURN
      END

C  /* Deck optmin */
      SUBROUTINE OPTMIN(WORK,LWORK,IWOFF,WRKDLM)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
      DIMENSION WORK(LWORK)
      KGRAD  = 1
      KSTEP  = KGRAD  + MXCOOR
      KHESS  = KSTEP  + MXCOOR
      KAHESS = KHESS  + MXCOOR*MXCOOR
      KTMAT  = KAHESS + MXCOOR*MXCOOR
      KPROJ  = KTMAT  + MXCOOR*MXCOOR
      KTEMPM = KPROJ  + MXCOOR*MXCOOR
      KTMPM2 = KTEMPM + 4*MXCOOR*MXCOOR
      KTMPM3 = KTMPM2 + 4*MXCOOR*MXCOOR
      KTMPM4 = KTMPM3 + 4*MXCOOR*MXCOOR
      KTMPM5 = KTMPM4 + 4*MXCOOR*MXCOOR
      KEVEC  = KTMPM5 + 4*MXCOOR*MXCOOR
      KSCLVC = KEVEC  + 4*MXCOOR*MXCOOR
      KGAMMA = KSCLVC + MXCOOR
      KGRDOL = KGAMMA + 2*MXCOOR
      KGRDMT = KGRDOL + 2*MXCOOR
      KSTPMT = KGRDMT + 2*MXCOOR*25
      KHESOL = KSTPMT + 2*MXCOOR*25
      KINFO  = KHESOL + 4*MXCOOR*MXCOOR
      KWILBM = KINFO  + (ITRMAX + 1)*6
      KBMTIN = KWILBM + 2*MXCOOR*MXCOOR
      KPJINM = KBMTIN + 2*MXCOOR*MXCOOR
      KHESIN = KPJINM + 4*MXCOOR*MXCOOR
      KWRK   = KHESIN + 4*MXCOOR*MXCOOR
      LWRK   = LWORK  - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('OPTMIN',' ',KWRK,LWORK)
#ifndef PRG_DIRAC
      IF (NEWTON .OR. QUADSD .OR. VISUAL) THEN
         CALL RN2ORD(IWOFF,WRKDLM,WORK(KGRAD),WORK(KSTEP),WORK(KHESS),
     &        WORK(KAHESS),WORK(KTMAT),WORK(KPROJ),WORK(KTEMPM),
     &        WORK(KTMPM2),WORK(KTMPM3),WORK(KTMPM4),WORK(KTMPM5),
     &        WORK(KEVEC), WORK(KSCLVC),WORK(KGAMMA),WORK(KGRDOL),
     &        WORK(KGRDMT),WORK(KSTPMT),WORK(KHESOL),WORK(KINFO),
     &        WORK(KWILBM),WORK(KBMTIN),WORK(KPJINM),WORK(KHESIN),
     &        WORK(KWRK),LWRK)
      ELSE
#endif
         CALL RN1ORD(IWOFF,WRKDLM,WORK(KGRAD),WORK(KSTEP),WORK(KHESS),
     &        WORK(KAHESS),WORK(KTMAT),WORK(KPROJ),WORK(KTEMPM),
     &        WORK(KTMPM2),WORK(KTMPM3),WORK(KTMPM4),WORK(KTMPM5),
     &        WORK(KEVEC), WORK(KSCLVC),WORK(KGAMMA),WORK(KGRDOL),
     &        WORK(KGRDMT),WORK(KSTPMT),WORK(KHESOL),WORK(KINFO),
     &        WORK(KWILBM),WORK(KBMTIN),WORK(KPJINM),WORK(KHESIN),
     &        WORK(KWRK),LWRK)
#ifndef PRG_DIRAC
      END IF
#endif
      RETURN
      END

#ifndef PRG_DIRAC
C  /* Deck rn2ord */
      SUBROUTINE RN2ORD(IWOFF,WRKDLM,EGRAD,CSTEP,EHESS,ALLHES,TMAT,
     &     PROJOP,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,EVEC,SCLVEC,
     &     GAMMA,GRDOLD,GRDMAT,STPMAT,HESOLD,GEINFO,WILBMT,BMTINV,
     &     PJINMT,HESINT,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <optinf.h>
#include <cbirea.h>
#include <molinp.h>
#include <abainf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <prkoor.h>
#include <trkoor.h>
#include <cbiwlk.h>
#include <hrunit.h>
#include <units.h>
      PARAMETER (LUERR = 0)
      LOGICAL EXHER, EXSIR, EXABA, MINEND, INDXOK, STATPO
      LOGICAL REJGEO, TRU, FAL, TMPLOG
      CHARACTER TMPLIN*80, WRDRSP*7
      DIMENSION EGRAD(MXCOOR), CSTEP(MXCOOR)
      DIMENSION EHESS(MXCOOR,MXCOOR), ALLHES(MXCOOR,MXCOOR)
      DIMENSION TMAT(MXCOOR,MXCOOR), PROJOP(MXCOOR,MXCOOR)
      DIMENSION TMPMAT(4*MXCOOR*MXCOOR), TMPMT2(4*MXCOOR*MXCOOR)
      DIMENSION TMPMT3(4*MXCOOR*MXCOOR), TMPMT4(4*MXCOOR*MXCOOR)
      DIMENSION TMPMT5(4*MXCOOR*MXCOOR)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR), SCLVEC(MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION GRDOLD(2*MXCOOR), GRDMAT(25,2*MXCOOR)
      DIMENSION STPMAT(25,2*MXCOOR), HESOLD(2*MXCOOR,2*MXCOOR)
C
C     The array geinfo contains optimization information for each
C     iteration. The first index is the iteration, the second gives   
C     the property:   1  -  Energy
C                     2  -  Gradient norm
C                     3  -  Index of Hessian
C                           (a negative index indicates symmetry break)
C                     4  -  Step length
C                     5  -  Trust radius
C                     6  -  # rejected steps
C
      DIMENSION GEINFO(0:ITRMAX,6)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR), BMTINV(2*MXCOOR,MXCOOR)
      DIMENSION PJINMT(2*MXCOOR,2*MXCOOR), HESINT(2*MXCOOR,2*MXCOOR)
      DIMENSION WORK(LWORK)
C
C     Initialization of variables.
C
      THRLDP = 1.0D-4
      THRIND = 5.0D-4
      TOLST  = 1.0D-5
      IPRWLK = IPRINT
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      CALL DZERO(GEINFO,(ITRMAX+1)*6)
      GEINFO(0,5) = TRSTRA
C
C     Initialization of new memory trap for programs.
C
      WORK(IWOFF+1) = WRKDLM
C
C     If visualization is requested, the program stops afterwards.
C
      IF (VISUAL) THEN
         CALL VISULZ(WORK,LWORK,IWOFF,WRKDLM,TMPMAT,EVEC)
         CALL QUIT('*** End of visualization ***')
         STOP '*** End of visualization ***'
      END IF
C
      CALL TITLER('2nd Order Geometry Minimization','>',110)
C
C     Perform preoptimization if requested.
C
      IF (DOPRE) CALL INIPRE(WORK,LWORK,IWOFF,WRKDLM,TMPMAT,EVEC)
C     
C     Calculate energy, gradient and Hessian.
C
      CALL GTHESS(EGRAD,EHESS,ALLHES,EXHER,EXSIR,EXABA,
     &                                        WORK,LWORK,IWOFF,WRKDLM)
C
C     Write Hessian to file (for 1st order restarts)
C
      CALL PNCHES(HESINT)
C
C     Make VRML-file of initial geometry if requested
C
      IF (VRML) CALL MKVRML(.FALSE.,TMPMAT,EVEC)
C
C     Count coordinates.
C
      CALL CNTCRD
C
C     Check if redundant internal coordinates should be used.
C
      IF (REDINT) THEN
         CALL INIRED(WILBMT,BMTINV,PJINMT,TMPMAT,TMPMT2,TMPMT3,TMPMT4)
      END IF
C
      CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
C
C     Construct projection operator and use it.
C     Then diagonalize Hessian.
C
      IF (REDINT) THEN
         CALL CGHINT(SCLVEC,HESOLD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,
     &               BMTINV,HESINT)
         CALL PRJINT(TMPMAT,PJINMT,HESINT)
         CALL DIAINT(EVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,THRIND,HESINT)
      ELSE
         CALL MKSCVC(SCLVEC)
         CALL PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
         CALL DIAHES(EGRAD,EHESS,ALLHES,TMAT,THRIND,
     &                                       EVEC,TMPMT2,WORK,LWORK)
      END IF
      GEINFO(0,1) = ENERGY
      GEINFO(0,3) = INDTOT*1.0D0
C
C     Determine step, check for convergence, print output and
C     and update geometry.
C
      IF (REDINT) THEN
         CALL FNSTIN(SCLVEC,EVEC,TMPMAT,TMPMT2,TMPMT3,CSTEP,BMTINV)
      ELSE
         CALL FNDSTP(EVEC,TMPMAT,TMPMT2,CSTEP)
      END IF
      GECONV = MINEND(SCLVEC)
      CALL PRIALL
      IF (.NOT. GECONV) CALL FNDGEO(CSTEP,EGRAD,TMPMAT,TMPMT2,
     &     EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM,IREJ,GEINFO)
      GEINFO(0,2) = GRADNM
      GEINFO(0,4) = STPNRM
      GEINFO(1,5) = TRSTRA
      GEINFO(0,6) = IREJ*1.0D0
C
C     Interactive runs print out useful information to LUERR
C
      IF (INTRAC) WRITE(LUERR,'(A,F14.6,A,F10.6,A,F10.6,A,I2)')
     &     ' Energy:', ENERGY, '     Grad.:', GRADNM,
     &     '     Step:', STPNRM,'     Index:', INDTOT
C
C     Make VRML-file of next geometry if requested
C
      IF (VRML) CALL MKVRML(.TRUE.,TMPMAT,EVEC)
C
C     Check if symmetry should be broken.
C
      IF (BRKSYM .AND. (.NOT. NOBRKS)) THEN
         CALL DOBRK(EXHER,EXSIR,EXABA,GEINFO,EVEC,CSTEP,TMPMAT,TMPMT2)
C
C     Check if preoptimization is finished
C
      ELSE IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
         CALL ENDPRE(EXHER,EXSIR,EXABA)
         KEEPHE = .FALSE.
      END IF
C
C     DO-WHILE loop that runs until geometry has converged or
C     maximum number of iterations is reached.
C
 10   CONTINUE
      IF ((ITRNMR .LT. ITRMAX) .AND. (.NOT. GECONV)) THEN
         ITRNMR = ITRNMR + 1
         ERGPRO = ERGPRD
         DO 20 I = 1, NCRTOT
            EVALOL(I) = EVAL(I)
 20      CONTINUE
C
C     Memory trap has to be reinitialized, then we go through the same
C     procedure as for the first iteration.
C
         WORK(IWOFF+1) = WRKDLM
         CALL GTHESS(EGRAD,EHESS,ALLHES,EXHER,EXSIR,EXABA,
     &                                    WORK,LWORK,IWOFF,WRKDLM)
         CALL PNCHES(HESINT)
C
C     If redundant internal coordinates are used, Wilson's B matrix,
C     its derivative and its inverse must be updated.
C
         IF (REDINT) THEN
            CALL GETWIL(TMPMAT,WILBMT)
            CALL GETDWL(TMPMAT,TMPMT2,TMPMT3,WILBMT)
            CALL GTBINV(TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTINV,
     &                  PJINMT)
         END IF
C
C     If new symmetry has been applied, coordinates has to counted again.
C
         IF (NWSYMM) THEN
            CALL CNTCRD
            CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
            NWSYMM = .FALSE.
         END IF
         IF (REDINT) THEN
            CALL CGHINT(SCLVEC,HESOLD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,
     &                  WILBMT,BMTINV,HESINT)
            CALL PRJINT(TMPMAT,PJINMT,HESINT)
            CALL DIAINT(EVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,THRIND,HESINT)
         ELSE
            CALL PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
            CALL DIAHES(EGRAD,EHESS,ALLHES,TMAT,THRIND,
     &           EVEC,TMPMT2,WORK,LWORK)
         END IF
         GEINFO(ITRNMR,1) = ENERGY
         GEINFO(ITRNMR,3) = INDTOT*1.0D0
         IF (REDINT) THEN
            CALL FNSTIN(SCLVEC,EVEC,TMPMAT,TMPMT2,TMPMT3,CSTEP,BMTINV)
         ELSE
            CALL FNDSTP(EVEC,TMPMAT,TMPMT2,CSTEP)
         END IF
         GECONV = MINEND(SCLVEC)
         CALL PRIALL
         IF ((.NOT. GECONV) .AND. ((.NOT. BRKSYM) .OR. NOBRKS))
     &        CALL FNDGEO(CSTEP,EGRAD,TMPMAT,TMPMT2,EXHER,
     &                EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM,IREJ,GEINFO)
         GEINFO(ITRNMR,2) = GRADNM
         GEINFO(ITRNMR,4) = STPNRM
         IF (ITRNMR .LT. ITRMAX) GEINFO(ITRNMR+1,5) = TRSTRA
         GEINFO(ITRNMR,6) = IREJ*1.0D0
         IF (INTRAC) WRITE(LUERR,'(A,F14.6,A,F10.6,A,F10.6,A,I2)')
     &        ' Energy:', ENERGY, '     Grad.:', GRADNM,
     &        '     Step:', STPNRM,'     Index:', INDTOT
C     
C     Check if symmetry should be broken
C
         IF (BRKSYM .AND. (.NOT. NOBRKS)) THEN
            CALL DOBRK(EXHER,EXSIR,EXABA,GEINFO,EVEC,CSTEP,
     &                                                 TMPMAT,TMPMT2)
C
C     Check if preoptimization is finished
C
         ELSE IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
            CALL ENDPRE(EXHER,EXSIR,EXABA)
            KEEPHE = .FALSE.
         END IF
C
C     Make VRML-file of next geometry if requested
C
         IF (VRML) CALL MKVRML(.TRUE.,TMPMAT,EVEC)
         GOTO 10
C
C     Finished case 1: Geometry has converged.
C
      ELSE IF (GECONV) THEN
C
C     We check to see if there is any properties to be calculated
C
Cjth         WORK(IWOFF+1) = WRKDLM
Cjth         WRINDX = .TRUE.
Cjth         CALL ABAINP('**FINAL',WORK(IWOFF + 2),LWORK)
Cjth         CALL EXEABA(WORK,LWORK,IWOFF,WRKDLM)
Cjth
CjthC     We also check to see if we have requested a RESPONSE calculations
CjthC
Cjth            REWIND (LUCMD)
CjthC
CjthC     Let us first see if we can find any **RESPONS input
CjthC
Cjth 133        CONTINUE
Cjth            READ (LUCMD,'(A7)',END=134,ERR=134) WRDRSP
Cjth            IF (WRDRSP .EQ. '**RESPO') THEN
Cjth               REWIND (LUCMD)
Cjth               RNRESP = .TRUE.
Cjth               CALL EXESIR(WORK,LWORK,IWOFF,WRKDLM)
Cjth               CALL RSPDRV(WORK(IWOFF + 2),LWORK)
Cjth            ELSE
Cjth               GOTO 133
Cjth            END IF
Cjth 134        CONTINUE 
C
C     Single point energy is calculated if requested.
C
         IF (DOSPE) CALL SPNRGY(GEINFO,EXHER,EXSIR,EXABA,WORK,LWORK,
     &                                                   IWOFF,WRKDLM)
C
C     Final results are printed, partially through PRIINF.
C
         CALL TITLER(' End of Minimization ','<',120)
         CALL PRIINF(GEINFO)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Geometry converged in ', ITRNMR+1,
     &        ' iterations!'
         IF (ITRBRK .GE. 0) THEN
            WRITE(LUPRI,*)
     &    'Please note that symmetry was broken during minimization.'
         ELSE IF (NOBRKS .AND. BRKSYM) THEN
            WRITE(LUPRI,*)
     &    'As specified, symmetry was not broken. Minimum was reached'
            WRITE(LUPRI,*)
     &    'within the given symmetry. Please note that breaking this'
            WRITE(LUPRI,*)
     &    'symmetry will decrease the energy further.'
         END IF
         ENERGY = GEINFO(ITRNMR,1)
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,F14.6,A)')
     &          ' Energy at final geometry is       : ',ENERGY,' a.u.'
         ERGDIF = ENERGY - GEINFO (0,1)
         WRITE(LUPRI,'(A,F14.6,A)')
     &          ' Energy change during minimization : ',ERGDIF,' a.u.'
         ERGDIF = ERGDIF * XKJMOL
         WRITE(LUPRI,'(A,F14.3,A)')
     &        '                                     ',ERGDIF,' kJ/mol'
         WRITE(LUPRI,*)
         IF (DOPRE) THEN
            WRITE(LUPRI,'(A)') ' Preoptimization was performed using'//
     &           ' the basis set(s):'
            DO 111 I = 1, IPRE-1
               WRITE(LUPRI,'(A,A60)') '       ',PREBTX(I)
 111        CONTINUE
            WRITE(LUPRI,*)
         END IF
         IF (DOSPE) THEN
            ENERGY = GEINFO(ITRNMR+1,1)
            WRITE(LUPRI,'(A,A60)') ' Using the basis ',SPBSTX
            WRITE(LUPRI,'(A,F14.6,A)')
     &          ' single point energy was calculated: ',ENERGY,' a.u.'
            WRITE(LUPRI,*)
         END IF
C
C     Finished case 2: Exceeded maximum number of iterations.
C
      ELSE
C     No single point energy has been calculated.
         TMPLOG = DOSPE
         DOSPE = .FALSE.
         CALL TITLER('Minimization Control Center','<',120)
         CALL PRIINF(GEINFO)
         DOSPE = TMPLOG
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Geometry has NOT converged!'
         WRITE(LUPRI,*) 'Maximum number of iterations (', ITRMAX,
     &                                        ') has been reached and'
         WRITE(LUPRI,*) 'minimization halted. Increase number or ',
     &                                   'restart from last geometry.'
         IF (DOSPE) WRITE(LUPRI,*) 'No single point energy has been '
     &                          // 'calculated.'
         WRITE(LUPRI,*)
      END IF
      RETURN
      END
#endif

C  /* Deck rn1ord */
      SUBROUTINE RN1ORD(IWOFF,WRKDLM,EGRAD,CSTEP,EHESS,ALLHES,TMAT,
     &     PROJOP,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,EVEC,SCLVEC,
     &     GAMMA,GRDOLD,GRDMAT,STPMAT,HESOLD,GEINFO,WILBMT,BMTINV,
     &     PJINMT,HESINT,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <optinf.h>
#include <cbirea.h>
#include <molinp.h>
#include <abainf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <prkoor.h>
#include <trkoor.h>
#include <cbiwlk.h>
#include <hrunit.h>
#include <units.h>
Cjth
#include <nuclei.h>
      PARAMETER (LUERR = 0)
      LOGICAL EXHER, EXSIR, EXABA, MINEND, INDXOK, STATPO
      LOGICAL REJGEO, TRU, FAL, TMPLOG
      CHARACTER TMPLIN*80, WRDRSP*7
      DIMENSION EGRAD(MXCOOR), CSTEP(MXCOOR)
      DIMENSION EHESS(MXCOOR,MXCOOR), ALLHES(MXCOOR,MXCOOR)
      DIMENSION TMAT(MXCOOR,MXCOOR), PROJOP(MXCOOR,MXCOOR)
      DIMENSION TMPMAT(4*MXCOOR*MXCOOR), TMPMT2(4*MXCOOR*MXCOOR)
      DIMENSION TMPMT3(4*MXCOOR*MXCOOR), TMPMT4(4*MXCOOR*MXCOOR)
      DIMENSION TMPMT5(4*MXCOOR*MXCOOR)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR), SCLVEC(MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION GRDOLD(2*MXCOOR), GRDMAT(25,2*MXCOOR)
      DIMENSION STPMAT(25,2*MXCOOR), HESOLD(2*MXCOOR,2*MXCOOR)
C
C     The array geinfo contains optimization information for each
C     iteration. The first index is the iteration, the second gives   
C     the property:   1  -  Energy
C                     2  -  Gradient norm
C                     3  -  Index of Hessian
C                           (a negative index indicates symmetry break)
C                     4  -  Step length
C                     5  -  Trust radius
C                     6  -  # rejected steps
C
      DIMENSION GEINFO(0:ITRMAX,6)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR), BMTINV(2*MXCOOR,MXCOOR)
      DIMENSION PJINMT(2*MXCOOR,2*MXCOOR), HESINT(2*MXCOOR,2*MXCOOR)
      DIMENSION WORK(LWORK)
C
C     Initialization of variables.
C
      THRLDP = 1.0D-4
      THRIND = 5.0D-4
      TOLST  = 1.0D-5
      IPRWLK = IPRINT
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      CALL DZERO(GEINFO,(ITRMAX+1)*6)
      GEINFO(0,5) = TRSTRA
C
C     Initialization of new memory trap for programs.
C
      WORK(IWOFF+1) = WRKDLM
      CALL TITLER('1st Order Geometry Minimization','>',110)
C
C     Perform preoptimization if requested.
C
      IF (DOPRE) CALL INIPRE(WORK,LWORK,IWOFF,WRKDLM,TMPMAT,EVEC)
C     
C     Calculate energy and gradient, initial Hessian is calculated
C     if requested.
C
      IF (INITHS) THEN
         CALL GTHESS(EGRAD,EHESS,ALLHES,EXHER,EXSIR,EXABA,
     &                                       WORK,LWORK,IWOFF,WRKDLM)
         
      ELSE
         CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
      END IF
C
C     Make VRML-file of initial geometry if requested
C
      IF (VRML) CALL MKVRML(.FALSE.,TMPMAT,EVEC)
C
C     Count coordinates.
C
      CALL CNTCRD
C
C     Check if redundant internal coordinates should be used.
C
      IF (REDINT .OR. INRDHS) THEN
         CALL INIRED(WILBMT,BMTINV,PJINMT,TMPMAT,TMPMT2,TMPMT3,TMPMT4)
      END IF
C
C     Initialize Hessian
C
      CALL INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,
     &            WILBMT,BMTINV,HESINT)
      CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
C
C     Construct projection operator and use it.
C     Then diagonalize Hessian.
C
      IF (REDINT) THEN
         CALL PRJINT(TMPMAT,PJINMT,HESINT)
         CALL DIAINT(EVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,THRIND,HESINT)
      ELSE
         CALL PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
         CALL DIAHES(EGRAD,EHESS,ALLHES,TMAT,THRIND,
     &                                     EVEC,TMPMT2,WORK,LWORK)
      END IF
      GEINFO(0,1) = ENERGY
      GEINFO(0,3) = INDTOT*1.0D0
C
C     Determine step, check for convergence, print output and
C     and update geometry.
C
      IF (REDINT) THEN
         CALL FNSTIN(SCLVEC,EVEC,TMPMAT,TMPMT2,TMPMT3,CSTEP,BMTINV)
      ELSE
         CALL FNDSTP(EVEC,TMPMAT,TMPMT2,CSTEP)
      END IF
      GECONV = MINEND(SCLVEC)
      CALL PRIALL
      IF (.NOT. GECONV) CALL FNDGEO(CSTEP,EGRAD,TMPMAT,TMPMT2,
     &     EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM,IREJ,GEINFO)
      GEINFO(0,2) = GRADNM
      GEINFO(0,4) = STPNRM
      GEINFO(1,5) = TRSTRA
      GEINFO(0,6) = IREJ*1.0D0
C
C     Interactive runs print out useful information to LUERR
C
      IF (INTRAC) WRITE(LUERR,'(A,F14.6,A,F10.6,A,F10.6,A,I2)')
     &     ' Energy:', ENERGY, '     Grad.:', GRADNM,
     &     '     Step:', STPNRM,'     Index:', INDTOT
C
C     Make VRML-file of next geometry if requested
C
      IF (VRML) CALL MKVRML(.TRUE.,TMPMAT,EVEC)
C
C     Check if preoptimization is finished
C
      IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
         CALL ENDPRE(EXHER,EXSIR,EXABA)
      END IF
C
C
C     DO-WHILE loop that runs until geometry has converged or
C     maximum number of iterations is reached.
C
 10   CONTINUE
      IF ((ITRNMR .LT. ITRMAX) .AND. (.NOT. GECONV)) THEN
         ITRNMR = ITRNMR + 1
         ERGPRO = ERGPRD
         DO 20 I = 1, NCRTOT
            EVALOL(I) = EVAL(I)
 20      CONTINUE
C
C     Memory trap has to be reinitialized, then we go through the same
C     procedure as for the first iteration.
C
         WORK(IWOFF+1) = WRKDLM
         CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
C
C     If redundant internal coordinates are used, Wilson's B matrix,
C     its derivative and its inverse must be updated.
C
         IF (REDINT) THEN
            CALL GETWIL(TMPMAT,WILBMT)
            CALL GETDWL(TMPMAT,TMPMT2,TMPMT3,WILBMT)
            CALL GTBINV(TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTINV,
     &                  PJINMT)
         END IF
C
C     If new symmetry has been applied, coordinates has to be counted
C     again and Hessian reinitialized.
C
         IF (NWSYMM) THEN
            CALL CNTCRD
            CALL INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,
     &                  TMPMT3,TMPMT4,WILBMT,BMTINV,HESINT)
            CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
            NWSYMM = .FALSE.
         ELSE
C
C     Otherwise the Hessian is updated.
C
            CALL UPDHES(SCLVEC,GRDOLD,GRDMAT,STPMAT,HESOLD,GAMMA,TMPMAT,
     &           TMPMT2,TMPMT3,TMPMT4,TMPMT5,WILBMT,BMTINV,HESINT)
            CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
         END IF
C
         IF (REDINT) THEN
            CALL PRJINT(TMPMAT,PJINMT,HESINT)
            CALL DIAINT(EVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,THRIND,HESINT)
         ELSE
            CALL PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
            CALL DIAHES(EGRAD,EHESS,ALLHES,TMAT,THRIND,
     &                                        EVEC,TMPMT2,WORK,LWORK)
         END IF
C
C     If the Hessian has negative eigenvalues, the Hessian should be
C     reinitialized, because minimization requires a positive
C     definite Hessian.
C
         IF (REINIH .AND. (INDTOT .GT. 0)) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)'***** NOTE! *****'
            WRITE(LUPRI,*)
     &      'Due to non-zero index, Hessian is reinitialized.'
            WRITE(LUPRI,*)
            CALL INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,
     &                  TMPMT3,TMPMT4,WILBMT,BMTINV,HESINT)
            CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
            IF (REDINT) THEN
               CALL PRJINT(TMPMAT,PJINMT,HESINT)
               CALL DIAINT(EVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,THRIND,
     &                     HESINT)
            ELSE
               CALL PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
               CALL DIAHES(EGRAD,EHESS,ALLHES,TMAT,THRIND,
     &                                           EVEC,TMPMT2,WORK,LWORK)
            END IF
         END IF
C
         GEINFO(ITRNMR,1) = ENERGY
         GEINFO(ITRNMR,3) = INDTOT*1.0D0
         IF (REDINT) THEN
            CALL FNSTIN(SCLVEC,EVEC,TMPMAT,TMPMT2,TMPMT3,CSTEP,BMTINV)
         ELSE
            CALL FNDSTP(EVEC,TMPMAT,TMPMT2,CSTEP)
         END IF
         GECONV = MINEND(SCLVEC)
         CALL PRIALL
         IF ((.NOT. GECONV) .AND. ((.NOT. BRKSYM) .OR. NOBRKS))
     &        CALL FNDGEO(CSTEP,EGRAD,TMPMAT,TMPMT2,EXHER,
     &                EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM,IREJ,GEINFO)
         GEINFO(ITRNMR,2) = GRADNM
         GEINFO(ITRNMR,4) = STPNRM
         IF (ITRNMR .LT. ITRMAX) GEINFO(ITRNMR+1,5) = TRSTRA
         GEINFO(ITRNMR,6) = IREJ*1.0D0
C
C     Check if rejected steps should cause reinitialization of Hessian.
C
         IF (REJINI .AND. (IREJ .GE. 1)) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)'***** NOTE! *****'
            WRITE(LUPRI,*)
     &           'Due to rejected step, Hessian is reinitialized.'
            WRITE(LUPRI,*)
            CALL INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,
     &           TMPMT3,TMPMT4,WILBMT,BMTINV,HESINT)
            TRSTRA = GEINFO(0,5)
            GEINFO(ITRNMR+1,5) = TRSTRA
         END IF
C
C     Check if increase of gradient norm should cause reinitialization of Hessian.
C     Reinitialization occurs when the norm of the gradient is larger than the
C     norm of the gradient two iterations earlier.
C
         IF (GRDINI .AND. (ITRNMR .GE. 2) .AND.
     &        (GEINFO(ITRNMR,2) .GE. GEINFO(ITRNMR-2,2))) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)'***** NOTE! *****'
            WRITE(LUPRI,*)'Due to increasing gradient norm, ' //
     &           'Hessian is reinitialized.'
            WRITE(LUPRI,*)
            CALL INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,
     &           TMPMT3,TMPMT4,WILBMT,BMTINV,HESINT)
            TRSTRA = GEINFO(0,5)
            GEINFO(ITRNMR+1,5) = TRSTRA
         END IF
C
         IF (INTRAC) WRITE(LUERR,'(A,F14.6,A,F10.6,A,F10.6,A,I2)')
     &        ' Energy:', ENERGY, '     Grad.:', GRADNM,
     &        '     Step:', STPNRM,'     Index:', INDTOT
C
C     Check if preoptimization is finished
C
         IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
            CALL ENDPRE(EXHER,EXSIR,EXABA)
         END IF
C
C     Make VRML-file of next geometry if requested
C
         IF (VRML) CALL MKVRML(.TRUE.,TMPMAT,EVEC)
         GOTO 10
C
C     Finished case 1: Geometry has converged.
C
      ELSE IF (GECONV) THEN
C
C     We check to see if there is any properties to be calculated
C
Cjth - no we don't!!
Cjth         WRINDX = .TRUE.
Cjth         CALL ABAINP('**FINAL',WORK(IWOFF + 2),LWORK)
Cjth         CALL EXEABA(WORK,LWORK,IWOFF,WRKDLM)
CjthC
CjthC     We also check to see if we have requested a RESPONSE calculations
CjthC
Cjth         REWIND (LUCMD)
CjthC
CjthC     Let us first see if we can find any **RESPONS input
CjthC
Cjth 133     CONTINUE
Cjth         READ (LUCMD,'(A7)',END=134,ERR=134) WRDRSP
Cjth         IF (WRDRSP .EQ. '**RESPO') THEN
Cjth            REWIND (LUCMD)
Cjth            RNRESP = .TRUE.
Cjth            CALL EXESIR(WORK,LWORK,IWOFF,WRKDLM)
Cjth            CALL RSPDRV(WORK(IWOFF + 2),LWORK)
Cjth         ELSE
Cjth            GOTO 133
Cjth         END IF
Cjth 134     CONTINUE
C
C     Single point energy is calculated if requested.
C
         IF (DOSPE) CALL SPNRGY(GEINFO,EXHER,EXSIR,EXABA,WORK,LWORK,
     &                                                   IWOFF,WRKDLM)
C
C     Final results are printed, partially through PRIINF.
C
         CALL TITLER(' End of Minimization ','<',120)
         CALL PRIINF(GEINFO)
C        Print angles, bondlength etc.
         CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,.FALSE.,
     &               WORK(IWOFF+2),LWORK)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Geometry converged in ', ITRNMR+1,
     &        ' iterations!'
         IF (ITRBRK .GE. 0) THEN
            WRITE(LUPRI,*)
     &    'Please note that symmetry was broken during minimization.'
         ELSE IF (NOBRKS .AND. BRKSYM) THEN
            WRITE(LUPRI,*)
     &    'As specified, symmetry was not broken. Minimum was reached'
            WRITE(LUPRI,*)
     &    'within the given symmetry. Please note that breaking this'
            WRITE(LUPRI,*)
     &    'symmetry will decrease the energy further.'
         END IF
         ENERGY = GEINFO(ITRNMR,1)
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,F14.6,A)')
     &          ' Energy at final geometry is       : ',ENERGY,' a.u.'
         ERGDIF = ENERGY - GEINFO (0,1)
         WRITE(LUPRI,'(A,F14.6,A)')
     &          ' Energy change during minimization : ',ERGDIF,' a.u.'
         ERGDIF = ERGDIF * XKJMOL
         WRITE(LUPRI,'(A,F14.3,A)')
     &        '                                     ',ERGDIF,' kJ/mol'
         WRITE(LUPRI,*)
         IF (DOPRE) THEN
            WRITE(LUPRI,'(A)') ' Preoptimization was performed using'//
     &           ' the basis set(s):'
            DO 111 I = 1, IPRE-1
               WRITE(LUPRI,'(A,A60)') '       ',PREBTX(I)
 111        CONTINUE
            WRITE(LUPRI,*)
         END IF
         IF (DOSPE) THEN
            ENERGY = GEINFO(ITRNMR+1,1)
            WRITE(LUPRI,'(A,A60)') ' Using the basis ',SPBSTX
            WRITE(LUPRI,'(A,F14.6,A)')
     &          ' single point energy was calculated: ',ENERGY,' a.u.'
            WRITE(LUPRI,*)
         END IF
C
C     Finished case 2: Exceeded maximum number of iterations.
C
      ELSE
C     No single point energy has been calculated.
         TMPLOG = DOSPE
         DOSPE = .FALSE.
         CALL TITLER('Minimization Control Center','<',120)
         CALL PRIINF(GEINFO)
C        Print angles, bondlength etc.
         CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,.FALSE.,
     &               WORK(IWOFF+1),LWORK)
         DOSPE = TMPLOG
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Geometry has NOT converged!'
         WRITE(LUPRI,*) 'Maximum number of iterations (', ITRMAX,
     &                                        ') has been reached and'
         WRITE(LUPRI,*) 'minimization halted. Increase number or ',
     &                                   'restart from last geometry.'
         IF (DOSPE) WRITE(LUPRI,*) 'No single point energy has been '
     &                          // 'calculated.'
         WRITE(LUPRI,*)
      END IF
      RETURN
      END

C  /* Deck priall */
      SUBROUTINE PRIALL
C
C     Prints important information for the current iteration.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
#include <trkoor.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <symmet.h>
#include <nuclei.h>
      CHARACTER*1 MARK
      CALL TITLER('Minimization Control Center',':',115)
      CALL HEADER('Optimization information',-1)
      WRITE(LUPRI,'(A,I7)')
     &     ' Iteration number               : ',ITRNMR
      MARK = ' '
      IF (BRKSYM .AND. GECONV) MARK = '*'
      IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) MARK = '*'
      WRITE(LUPRI,'(A,L1,A1)')
     &      ' End of minimization            :       ', GECONV, MARK
      WRITE(LUPRI,678) ' Energy at this geometry is     : ', ENERGY
      IF ((ITRNMR .GT. 0) .AND. (ITRBRK .LT. (ITRNMR - 1))) THEN
         ERGDIF = ENERGY - ERGOLD
         IF (ABS(ERGPRO) .GT. 1.0D-10) THEN
            RATIO = ERGDIF / ERGPRO 
         ELSE
            RATIO = 1.0D0
         END IF
 678     FORMAT(A,F14.6)
         WRITE(LUPRI,678)' Energy change from last geom.  : ', ERGDIF
         WRITE(LUPRI,678)' Predicted change               : ', ERGPRO
         WRITE(LUPRI,678)' Ratio, actual/predicted change : ', RATIO
      END IF
      WRITE(LUPRI,678) ' Norm of gradient               : ', GRADNM
      DO 923 I = 0, MAXREP
         IF (INDHES(I) .GT. 0) WRITE(LUPRI,'(A,I1,A,I7)')
     &        ' Hessian index (',I,')              : ', INDHES(I)
 923  CONTINUE
      WRITE(LUPRI,'(A,I7)') ' Total Hessian index            : ',INDTOT
      IF (BRKSYM .AND. GECONV) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '*) Within given symmetry.'
      ELSE IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '*) End of preoptimization.'
      END IF
      CALL HEADER('Change in geometry',-1)
      WRITE(LUPRI,678) ' Updated trust radius           : ', TRSTRA
      WRITE(LUPRI,678) ' Norm of step                   : ', STPNRM
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Eigenvalues',-1)
         WRITE(LUPRI,*)
     &        '  #      Current value  Previous value      Change   '
         WRITE(LUPRI,*)
     &        '-----------------------------------------------------'
         DO 10 I = 1, NCRTOT
            WRITE(LUPRI,'(I4,3F16.6)') I, EVAL(I), EVALOL(I),
     &                                       EVAL(I)-EVALOL(I)
 10      CONTINUE
      END IF
      RETURN
      END

C  /* Deck priinf */
      SUBROUTINE PRIINF(GEINFO)
C
C     Prints important information at the end of an optimization.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <optinf.h>
      DIMENSION GEINFO(0:ITRMAX,6)
      CHARACTER NWSMRK*1,BRKMRK*1
      TMPNRG = 0.0D0
#ifdef PRG_DIRAC
      CALL HEADE2('Final geometry','@',-1)
      CALL RPRIGEO(CORD,'@')
#else
      CALL HEADER('Final geometry','@',-1)
      CALL PRIGEO(CORD)
#endif
      WRITE(LUPRI,*)
     & 'Iter     Energy     Change    GradNorm  Index   StepLen'
     & // '    TrustRad #Rej'
      WRITE(LUPRI,*)
     & '------------------------------------------------------'
     & // '-----------------'
      DO 10 I = 0, ITRNMR
         IF (I .GT. 0) TMPNRG = GEINFO(I,1) - GEINFO(I-1,1)
C
C     There's two special marks for each iteration, one (*) for a Newton
C     Newton step (that is a step smaller than the trust radius) and
C     one (x) for the breaking of symmetry.
C
         NWSMRK = ' '
         IF ((GEINFO(I,4)+1.0D-6) .LT. GEINFO(I,5)) NWSMRK = '*'
         BRKMRK = ' '
         IF (NINT(GEINFO(I,3)) .LT. 0) BRKMRK = 'x'
         WRITE(LUPRI,'(I4,A1,F12.6,2F11.6,I5,F12.6,A1,F11.6,I4)') I,
     &        BRKMRK,GEINFO(I,1),TMPNRG,GEINFO(I,2),
     &        NINT(ABS(GEINFO(I,3))),GEINFO(I,4),NWSMRK,GEINFO(I,5),
     &        NINT(GEINFO(I,6))
 10   CONTINUE
C
C     We also write the single point energy, if calculated.
C
      IF (DOSPE) WRITE(LUPRI,'(A5,F12.6,F11.6)') '  SP ',
     &          GEINFO(ITRNMR+1,1),GEINFO(ITRNMR+1,1)-GEINFO(ITRNMR,1)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '*) Newton step taken.'
      IF (ITRBRK .GE. 0) WRITE(LUPRI,*)
     &                 'x) Symmetry was broken after this iteration.'
      RETURN
      END

C  /* Deck gtnrgy */
      SUBROUTINE GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
C
C     Retreives the energy of the current geometry
C     (by running HERMIT and SIRIUS).
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <infopt.h>
#include <optinf.h>
#include <siripc.h>
#ifdef PRG_DIRAC
#include <dcbdhf.h>
#include <dcbgen.h>
#else
#include <gnrinf.h>
#include <infinp.h>
#endif
      LOGICAL EXHER, EXSIR, EXABA
      DIMENSION WORK(LWORK)
C
C     jth - relativistic
C
C
C     Integral and basis input
C
      IF ( .NOT. EXSIR ) THEN
         CALL EXEDIR(WORK(IWOFF+1),LWORK)
         ENERGY = DHFERG
         EXHER = .TRUE.
         EXSIR = .TRUE.
      END IF
      NEWGEO = .TRUE.
C
C      DIMENSION WORK(LWORK)
C      IF (.NOT. EXHER) THEN
C         CALL EXEHER(WORK,LWORK,IWOFF,WRKDLM)
C         EXHER = .TRUE.
C      END IF
C      NEWGEO = .TRUE.
C      IF (.NOT. EXSIR) THEN
C         CALL EXESIR(WORK,LWORK,IWOFF,WRKDLM)
C         IF (DOCCSD) THEN
C            CALL EXECC(WORK,LWORK,IWOFF,WRKDLM)
C         ENDIF
C         EXSIR = .TRUE.
C      END IF
C      IF (FLAG(16)) THEN
C         WRITE (LUPRI,'(/A)') ' Solvent model not allowed with new'//
C     &        'optimization module. Use .WALK'
C         CALL QUIT('Solvent incompatible with .MINIMIZE')
C      END IF
C     
C
C     In Dirac THRGRD is called SCFCNV if EVCCNV is .true.
C
      IF ( EVCCNV .OR. FCKCNV ) THEN
         THRGRD = SCFCNV
      ELSE 
         THRGRD = SQRT(SCFCNV)
      END IF
C 
      IF ((2.0D0*THRGRD) .GT. GRDTHR .AND. .NOT. CHGRDT) THEN
         GRDTHR = 2.0D0*THRGRD
         THRSTP = GRDTHR
         THRERG = GRDTHR
         THRSYM = SQRT(THRERG)
         WRITE (LUPRI,'(/A)') ' WARNING>>>> Due to convergence '//
     &        'thresholds for the wave function'
         WRITE (LUPRI,'(A)') ' WARNING>>>> thresholds for convergence'
     &        //' of geometry optimization has been reset'
         WRITE (LUPRI,'(/A,/,3(/,20X,A,F10.8))') ' New thresholds:',
     &        'Gradient norm  ',GRDTHR,'Step norm      ',THRSTP,
     &        'Energy change  ',THRERG
      END IF
Cjth - nonrel
C      IF (ISTATE .GT. 1 .AND. ITRNMR .EQ. 0) THEN
C         WRITE (LUPRI,'(/A)') ' WARNING>>>> As this is an excited '//
C     &        'optimization I have updated some thresholds'
C         IF (TRSTRA .EQ. 0.50D0) TRSTRA = 0.25D0
C         IF (RTENBD .EQ. 0.4D0) RTENBD = 0.85D0
C         IF (RTRJMN .EQ. 0.1D0) RTRJMN = 0.80D0
C         IF (RTENGD .EQ. 0.80D0) RTENGD = 0.90D0
C         IF (RTRJMX .EQ. 3.0D0) RTRJMX = 1.20D0
C         WRITE(LUPRI,'(/A/A/,(A,F10.4))')
C     &        ' Restricted step control parameters',
C     &        ' ----------------------------------',
C     &        ' Initial trust radius   :',TRSTRA,
C     &        ' Trust radius increment :',TRSTIN,
C     &        ' Trust radius decrement :',TRSTDE,
C     &        ' Bad prediction ratio   :',RTENBD,
C     &        ' Good prediction ratio  :',RTENGD,
C     &        ' Rejection ratio, low   :',RTRJMN,
C     &        ' Rejection ratio, high  :',RTRJMX
C      END IF
C      IF ((DOMP2 .AND. .NOT. DOMC) .OR. DOCCSD .OR. ((DOCI .OR. DOCINO)
C     &    .AND. .NOT. DOMC)) THEN
C         ENERGY = ECORR
C      ELSE
C         ENERGY = EMCSCF
C      END IF
      RETURN
      END

C  /* Deck gtgrad */
      SUBROUTINE GTGRAD(EGRAD,EXHER,EXSIR,EXABA,
     &                             WORK,LWORK,IWOFF,WRKDLM)
C
C     Retreives the energy and gradient of the current geometry
C     (by running HERMIT, SIRIUS and ABACUS).
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxorb.h>
#include <taymol.h>
#include <maxaqn.h>
#include <symmet.h>
#include <nuclei.h>
#include <abainf.h>
#include <optinf.h>
Cjth
#ifdef PRG_DIRAC
#include <dcbgrd.h>
#include <dcbdhf.h>
#include <dcbgen.h>
#else
#include <gnrinf.h>
#include <infinp.h>
#endif
      LOGICAL EXHER, EXSIR, EXABA, LSAVE1, LSAVE2(3)
      DIMENSION WORK(LWORK)
      DIMENSION EGRAD(MXCOOR)
C
C     
C     Get energy if requested
C
      CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
C
C
C     Get gradient
C
      CALL DZERO(EGRAD,MXCOOR)
      CALL DZERO(EGRAD,MXCOOR)
C
C     Calculate gradient
C
      CALL PSIGRD(WORK(IWOFF+1),LWORK)
C
      EXABA = .TRUE.
      ENERGY = DHFERG
C     Non-relativistic:...
C      CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
CC
CC     Let's see if this is a wave function we know, if not: Do it numerically
CC
C      IF ((DOMP2 .AND. .NOT. DOMC) .OR. DOCCSD .OR. ((DOCI .OR. DOCINO)
C     &    .AND. .NOT. DOMC)) THEN
C         CALL NUMGRD(WORK,LWORK,IWOFF,WRKDLM)
C      ELSE
C         WRINDX = .TRUE.
C         MOLGRD = .TRUE.
C         IF (ITRNMR .EQ. 0) THEN
C            CALL ABAINP('**START',WORK(IWOFF + 2),LWORK)
C         ELSE
C            CALL ABAINP('**PROPE',WORK(IWOFF + 2),LWORK)
C         END IF
C         CALL EXEABA(WORK,LWORK,IWOFF,WRKDLM)
C         EXABA = .TRUE.
C         MOLGRD = .FALSE.
C      END IF
      CALL TRAGRD(GRDMOL,EGRAD,WORK(IWOFF + 2),
     &            WORK(IWOFF + MXCOOR*MXCOOR + 2),MXCOOR,3*NUCDEP)
      RETURN
      END
C  /* Deck gthess */
      SUBROUTINE GTHESS(EGRAD,EHESS,ALLHES,EXHER,EXSIR,EXABA,
     &                                 WORK,LWORK,IWOFF,WRKDLM)
C
C     Retreives the energy, gradient and hessian of the current geometry
C     (by running HERMIT, SIRIUS and ABACUS).
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <taymol.h>
#include <maxaqn.h>
#include <symmet.h>
#include <nuclei.h>
#include <abainf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
#include <trkoor.h>
#include <cbiwlk.h>
Cjth - this routine is .NOT. working with DIRAC!!!!
#if !defined (PRG_DIRAC)
      LOGICAL EXHER, EXSIR, EXABA
      DIMENSION WORK(LWORK)
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION ALLHES(NCRTOT*NCRTOT)
      CALL DZERO(EGRAD,MXCOOR)
      CALL DZERO(EHESS,MXCOOR*MXCOOR)
      CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
      WRINDX = .TRUE.
      MOLGRD = .TRUE.
      MOLHES = .TRUE.
      IF (ITRNMR .EQ. 0) THEN
         CALL ABAINP('**START',WORK(IWOFF + 2),LWORK)
      ELSE
         CALL ABAINP('**PROPE',WORK(IWOFF + 2),LWORK)
      END IF
      CALL EXEABA(WORK,LWORK,IWOFF,WRKDLM)
      EXABA = .TRUE.
      START = .FALSE.
      MOLGRD = .FALSE.
      MOLHES = .FALSE.
C
C     For all iterations except the first, MINCGH is run here.
C
      IF ((ITRNMR .GT. 0) .AND. (.NOT. NWSYMM)) THEN
         CALL MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
      END IF
#endif
      RETURN
      END
C  /* Deck inihes */
      SUBROUTINE INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,
     &                  TMPMT3,TMPMT4,WILBMT,BMTINV,HESINT)
C
C     Initializes Hessian
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <optinf.h>
#include <taymol.h>
#include <nuclei.h>
#include <symmet.h>
      DIMENSION SCLVEC(MXCOOR), GRDOLD(2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), TMPMAT(2*MXCOOR,2*MXCOOR)
      DIMENSION TMPMT2(2*MXCOOR,2*MXCOOR), TMPMT3(2*MXCOOR,2*MXCOOR)
      DIMENSION TMPMT4(2*MXCOOR,2*MXCOOR), WILBMT(2*MXCOOR,MXCOOR)
      DIMENSION BMTINV(2*MXCOOR,2*MXCOOR), HESINT(2*MXCOOR,2*MXCOOR)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL DZERO(HESOLD,4*MXCOOR*MXCOOR)
      CALL DZERO(GRDOLD,2*MXCOOR)
      CALL DZERO(SCLVEC,MXCOOR)
C
C     We initialize the scaling vector (the inverse of the
C     normalization vector in WLKCGH).
C
      DO 1 IREP = 0, MAXREP
         DO 2 ICENT = 1, NUCIND
#include <scalar.h>
            DO 3 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  SCLVEC(ISCOOR) = SQRT(FMULT(ISTBNU(ICENT)))
               END IF
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
C
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from INIHES','*',103)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Scaling vector',-1)
            CALL OUTPUT(SCLVEC,1,1,1,NCRTOT,1,MXCOOR,1,LUPRI)
         END IF
      END IF
C
C     Check if Hessian should be read from file
C
      IF (HSFILE) THEN
         CALL REAHES(HESINT,IERR)
         IF (IERR .EQ. -1) THEN
            CALL QUIT('Unable to open the file DALTON.HES.')
         ELSE IF (IERR .EQ. -2) THEN
            CALL QUIT('The Hessian in DALTON.HES has ' //
     &                                     'wrong dimensions.')
         END IF
         IF (REDINT) THEN
            IF (MAXREP .GT. 0) THEN
               CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
               CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                     NCRREP(0,1),3*NUCDEP)
               CALL GX2GQ(TMPMAT,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
            ELSE
               CALL GX2GQ(GRDMOL,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
            END IF
            DO 5 J = 1, IINTCR
               GRDOLD(J) = GRDINT(J)
               DO 7 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 7             CONTINUE
 5          CONTINUE
         ELSE
            DO 10 I = NCART+1, NCRTOT
               HESMOL(I,I) = SCLVEC(I)*SCLVEC(I)
 10         CONTINUE
            DO 12 J = 1, NCRTOT
               GRDOLD(J) = GRDMOL(J)
               DO 14 I = 1, NCRTOT
                  HESOLD(I,J) = HESMOL(I,J)/(SCLVEC(I)*SCLVEC(J))
 14            CONTINUE
 12         CONTINUE
         END IF
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian has been read from file.'
         WRITE(LUPRI,*)
         HSFILE = .FALSE.
C
C     Check if Hessian has been calculated
C
      ELSE IF (INITHS) THEN
         IF (REDINT) THEN
            IF (MAXREP .GT. 0) THEN
               CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
               CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                     NCRREP(0,1),3*NUCDEP)
               CALL GX2GQ(TMPMAT,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian gradient',-1)
                  CALL OUTPUT(TMPMAT,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
                  CALL HEADER('Internal gradient',-1)
                  CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
               END IF
               CALL DZERO(HESOLD,4*MXCOOR*MXCOOR)
               CALL TRAHES(HESMOL,MXCOOR,HESOLD,TMPMT2,TMPMT3,
     &                     MXCOOR,3*NUCDEP,1)
               CALL HX2HQ(TMPMAT,TMPMT2,TMPMT3,TMPMT4,HESOLD,
     &                    GRDINT,HESINT,WILBMT,BMTINV)
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian Hessian',-1)
                  CALL OUTPUT(HESOLD,1,ICRTCR,1,ICRTCR,
     &                                      MXCOOR,MXCOOR,1,LUPRI)
               END IF
               CALL DZERO(HESOLD,4*MXCOOR*MXCOOR)
            ELSE
               CALL GX2GQ(GRDMOL,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian gradient',-1)
                  CALL OUTPUT(GRDMOL,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
                  CALL HEADER('Internal gradient',-1)
                  CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
               END IF
               CALL HX2HQ(TMPMAT,TMPMT2,TMPMT3,TMPMT4,HESMOL,
     &                    GRDINT,HESINT,WILBMT,BMTINV)
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian Hessian',-1)
                  CALL OUTPUT(HESMOL,1,ICRTCR,1,ICRTCR,
     &                                      MXCOOR,MXCOOR,1,LUPRI)
               END IF
            END IF
            DO 15 J = 1, IINTCR
               DO 16 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 16            CONTINUE
               GRDOLD(J) = GRDINT(J)
 15         CONTINUE
         ELSE
            DO 18 J = 1, NCRTOT
               DO 19 I = 1, NCRTOT
                  HESOLD(I,J) = HESMOL(I,J)/(SCLVEC(I)*SCLVEC(J))
 19            CONTINUE
               GRDOLD(J) = GRDMOL(J)
 18         CONTINUE
         END IF
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian has been calculated.'
         WRITE(LUPRI,*)
C
C     Next line makes reinitialization possible:
C
         INITHS = .FALSE.
C
C     Check if initial Hessian should be diagonal in redundant internal
C     coordinates.
C
      ELSE IF (INRDHS .OR. REDINT) THEN
         IF (MAXREP .GT. 0) THEN
            CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
            CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                  NCRREP(0,1),3*NUCDEP)
            CALL GX2GQ(TMPMAT,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
         ELSE
            CALL GX2GQ(GRDMOL,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
         END IF
C
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Cartesian gradient',-1)
            CALL OUTPUT(GRDMOL,1,1,1,NCART,1,MXCOOR,1,LUPRI)
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Internal gradient',-1)
            CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         END IF
C
C     The default is different values for bonds and other internal
C     coordinates.
C
         CALL DZERO(HESINT,4*MXCOOR*MXCOOR)
         IF (EVLINI .LE. 0.0D0) THEN
            DO 20 I = 1, IINTCR
               IF (INTCRD(I,1) .LT. 10) THEN
                  HESINT(I,I) = 0.5D0
               ELSE
                  HESINT(I,I) = 0.2D0
               END IF
 20         CONTINUE
         ELSE
            DO 22 I = 1, IINTCR
               HESINT(I,I) = EVLINI
 22         CONTINUE
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Internal Hessian',-1)
            CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,
     &                                       2*MXCOOR,2*MXCOOR,1,LUPRI)
         END IF
         IF (INRDHS) THEN
            CALL DZERO(HESMOL,MXCOOR*MXCOOR)
            CALL HQ2HX(TMPMAT,TMPMT2,TMPMT3,HESINT,GRDINT,HESMOL,WILBMT)
            IF (MAXOPR .GT. 1) THEN
               CALL DZERO(TMPMT4,4*MXCOOR*MXCOOR)
               DO 100 J = 1, ICRTCR
                  DO 102 I = 1, ICRTCR
                     TMPMT4(I,J) = HESMOL(I,J)
 102              CONTINUE
 100           CONTINUE
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Expanded Cartesian Hessian',-1)
                  CALL OUTPUT(HESMOL,1,ICRTCR,1,ICRTCR,MXCOOR,MXCOOR,
     &                                                       1,LUPRI)
               END IF
               DO 715 J = 1, ICRTCR
                  DO 717 I = 1, J-1
                     HESMOL(J,I) = HESMOL(I,J)
 717              CONTINUE
 715           CONTINUE
C
C     Collect the symmetry-to-cartesian transformation matrix in TMPMT3
C
               CALL TRACOR(TMPMT2,TMPMT3,1,ICRTCR,0)
C
               CALL MPATB(TMPMT3,ICRTCR,ICRTCR,ICRTCR,ICRTCR,
     &                    HESMOL,ICRTCR,ICRTCR,MXCOOR,MXCOOR,
     &                    TMPMAT,              MXCOOR,MXCOOR)
               CALL DZERO(HESMOL,MXCOOR*MXCOOR)
               CALL  MPAB(TMPMAT,ICRTCR,ICRTCR,MXCOOR,MXCOOR,
     &                    TMPMT3,ICRTCR,ICRTCR,ICRTCR,ICRTCR,
     &                    HESMOL,              MXCOOR,MXCOOR)
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Totally Symmetric Cartesian Hessian',-1)
                  CALL OUTPUT(HESMOL,1,NCART,1,NCART,MXCOOR,MXCOOR,
     &                                                         1,LUPRI)
               END IF
            END IF
            DO 25 J = 1, ICRTCR
               GRDOLD(J) = GRDMOL(J)
               DO 27 I = 1, ICRTCR
                  HESOLD(I,J) = HESMOL(I,J)*SCLVEC(I)*SCLVEC(J)
                  HESMOL(I,J) = HESMOL(I,J)*SCLVEC(I)*SCLVEC(I)
     &                                             *SCLVEC(J)*SCLVEC(J)
 27            CONTINUE
 25         CONTINUE
         ELSE
            DO 30 J = 1, IINTCR
               GRDOLD(J) = GRDINT(J)
               DO 32 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 32            CONTINUE
 30         CONTINUE
         END IF
C
C     Otherwise the Hessian is set equal to a diagonal matrix. When
C     symmetry is used, the Hessian must be scaled so that it will be
C     correct after normalization.
C
      ELSE
         CALL DZERO(HESMOL,MXCOOR*MXCOOR)
         DO 40 I = 1, NCRTOT
            HESMOL(I,I) = EVLINI*SCLVEC(I)*SCLVEC(I)
            HESOLD(I,I) = EVLINI
            GRDOLD(I) = GRDMOL(I)
 40      CONTINUE
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian is a diagonal matrix.'
         WRITE(LUPRI,*)
      END IF
C
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Initial Hessian',-1)
         IF (REDINT) THEN
            CALL OUTPUT(HESOLD,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,
     &                                                       1,LUPRI)
         ELSE
            CALL OUTPUT(HESOLD,1,NCART,1,NCART,2*MXCOOR,2*MXCOOR,
     &                                                       1,LUPRI)
         END IF
      END IF
C
C     We punch the Hessian to file
C
      CALL PNCHES(HESINT)
      RETURN
      END

C  /* Deck updhes */
      SUBROUTINE UPDHES(SCLVEC,GRDOLD,GRDMAT,STPMAT,HESOLD,GAMMA,TMPMAT,
     &     TMPMT2,TMPMT3,TMPMT4,TMPMT5,WILBMT,BMTINV,HESINT)
C
C     Controls Hessian updates
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <optinf.h>
#include <symmet.h>
#include <nuclei.h>
#include <taymol.h>
      DIMENSION SCLVEC(MXCOOR), GRDOLD(2*MXCOOR)
      DIMENSION GRDMAT(25,2*MXCOOR), STPMAT(25,2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION TMPMAT(4*MXCOOR*MXCOOR), TMPMT2(4*MXCOOR*MXCOOR)
      DIMENSION TMPMT3(4*MXCOOR*MXCOOR), TMPMT4(4*MXCOOR*MXCOOR)
      DIMENSION TMPMT5(4*MXCOOR*MXCOOR), WILBMT(2*MXCOOR,MXCOOR)
      DIMENSION BMTINV(2*MXCOOR,MXCOOR), HESINT(2*MXCOOR,2*MXCOOR)
      LOGICAL RESET
C
C     If redundant internal coordinates are used, we have to transform
C     the calculated Cartesian gradient
C
      IF (REDINT) THEN
         IF (MAXREP .GT. 0) THEN
            CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
            CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                  NCRREP(0,1),3*NUCDEP)
            CALL GX2GQ(TMPMAT,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
         ELSE
            CALL GX2GQ(GRDMOL,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
         END IF
      END IF
C
C     When a new basis set is used in preoptimization, the Hessian is
C     kept as it is, because no step is taken.
C
      IF (KEEPHE) THEN
         CALL UPDOLD(SCLVEC,GRDOLD,HESOLD,HESINT)
C     
C     The requested Hessian update method is used
C
      ELSE 
         IF (BFGS) THEN 
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPBFGS')
            IF (REDINT) THEN
               CALL UPBFGS(STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,2*MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPBFGS(STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (DFP) THEN 
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPDFP ')
            IF (REDINT) THEN
               CALL UPDFP(STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,2*MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPDFP(STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (PSB) THEN 
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPPSB ')
            IF (REDINT) THEN
               CALL UPPSB(STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,2*MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPPSB(STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (RANKON) THEN 
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPRNKO')
            IF (REDINT) THEN
               CALL UPRNKO(STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,2*MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPRNKO(STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (SCHLEG) THEN
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPSCHL')
            RESET = .FALSE.
            IF ((ITRNMR .LT. 2) .OR. (ITRBRK .EQ. (ITRNMR-2)))
     &           RESET = .TRUE.
            IF (REDINT) THEN
               CALL UPSCHL(STPINT,STPMAT,GAMMA,GRDMAT,HESOLD,HESINT,
     &              IINTCR,2*MXCOOR,TMPMAT(1),TMPMAT(2*MXCOOR+1),
     &              TMPMT2,TMPMT3,NCART-NPROJ,RESET,GRADNM,IPRINT)
            ELSE
               CALL UPSCHL(STPSYM,STPMAT,GAMMA,GRDMAT,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT(1),TMPMAT(2*MXCOOR+1),
     &              TMPMT2,TMPMT3,NCART-NPROJ,RESET,GRADNM,IPRINT)
            END IF
         ELSE IF (MBFGS) THEN
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPMULT')
            RESET = .FALSE.
            IF ((ITRNMR .LT. 2) .OR. (ITRBRK .EQ. (ITRNMR-2)))
     &           RESET = .TRUE.
            IF (REDINT) THEN
               CALL UPMULT(STPINT,STPMAT,GAMMA,GRDMAT,HESOLD,HESINT,
     &              IINTCR,2*MXCOOR,TMPMAT(1),TMPMAT(2*MXCOOR+1),TMPMT2,
     &              TMPMT3,TMPMT4,TMPMT5,NCART-NPROJ,RESET,
     &              GRADNM,IPRINT)
            ELSE
               CALL UPMULT(STPSYM,STPMAT,GAMMA,GRDMAT,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT(1),TMPMAT(2*MXCOOR+1),TMPMT2,
     &              TMPMT3,TMPMT4,TMPMT5,NCART-NPROJ,RESET,
     &              GRADNM,IPRINT)
            END IF
         ELSE IF (STEEPD) THEN
            CALL UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,'UPSTPD')
            CALL INIHES(SCLVEC,GRDOLD,HESOLD,TMPMAT,TMPMT2,
     &                  TMPMT3,TMPMT4,WILBMT,BMTINV,HESINT)

            IF (.NOT. REDINT) THEN
               DO 10 J = 1, NCART
                  DO 12 I = 1, NCART
                     HESMOL(I,J) = HESMOL(I,J)/(SCLVEC(I)*SCLVEC(J))
 12               CONTINUE
 10            CONTINUE
            END IF
         END IF
C
C     In all higher symmetries, we set all diagonal elements to 1
C
         IF (.NOT. REDINT) THEN
            DO 15 I = (NCART + 1), NCRTOT
               HESMOL(I,I) = SCLVEC(I)*SCLVEC(I)
 15         CONTINUE
C
C     We copy the updated Hessian to HESOLD and the gradient to GRDOLD,
C     making everything ready for the next iteration. The updated Hessian
C     then has to be "off-scaled" before it's sent to MINCGH.
C
            CALL DZERO (HESOLD,4*MXCOOR*MXCOOR)
            CALL DZERO (GRDOLD,MXCOOR)
            DO 20 J = 1, NCART
               DO 22 I = 1, NCART
                  HESOLD(I,J) = HESMOL(I,J)
                  HESMOL(I,J) = HESMOL(I,J)*SCLVEC(I)*SCLVEC(J)
 22            CONTINUE
               GRDOLD(J) = GRDMOL(J)
 20         CONTINUE
         ELSE
            CALL DZERO (HESOLD,4*MXCOOR*MXCOOR)
            CALL DZERO (GRDOLD,MXCOOR)
            DO 30 J = 1, IINTCR
               DO 32 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 32            CONTINUE
               GRDOLD(J) = GRDINT(J)
 30         CONTINUE
         END IF
      END IF
C     
C     We punch the Hessian to file
C     
      CALL PNCHES(HESINT)
      RETURN
      END

C  /* Deck upoutp */
      SUBROUTINE UPOUTP(SCLVEC,GRDOLD,GAMMA,HESOLD,UPDTXT)
C
C     Some common output for the updating subroutines.
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <optinf.h>
#include <taymol.h>
#include <nuclei.h>
#include <symmet.h>
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      DIMENSION SCLVEC(MXCOOR), GRDOLD(2*MXCOOR) 
      DIMENSION GAMMA(2*MXCOOR), HESOLD(2*MXCOOR,2*MXCOOR) 
      CHARACTER UPDTXT*6, OUTTXT*18
C
      NCOOR = NCART
      IF (REDINT) NCOOR = IINTCR
      OUTTXT = 'Output from ' // UPDTXT
      IF (IPRINT .GT. IPRMIN) THEN
         CALL TITLER(OUTTXT,'*',103)
         IF ((IPRINT .GE. IPRDBG) .AND. (.NOT. REDINT)) THEN
            CALL HEADER('Scaling vector in UPOUTP',-1)
            CALL OUTPUT(SCLVEC,1,1,1,NCART,1,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Step from last geometry',-1)
         IF (REDINT) THEN
            CALL OUTPUT(STPINT,1,1,1,NCOOR,1,2*MXCOOR,1,LUPRI)
         ELSE
            CALL OUTPUT(STPSYM,1,1,1,NCOOR,1,MXCOOR,1,LUPRI)
         END IF
         CALL HEADER('Previous Hessian',-1)
         CALL OUTPUT(HESOLD,1,NCOOR,1,NCOOR,2*MXCOOR,2*MXCOOR,1,LUPRI)
         CALL HEADER('Gradient at last geometry',-1)
         CALL OUTPUT(GRDOLD,1,1,1,NCOOR,1,2*MXCOOR,1,LUPRI)
         CALL HEADER('Gradient at current geometry',-1)
         IF (REDINT) THEN
            CALL OUTPUT(GRDINT,1,1,1,NCOOR,1,2*MXCOOR,1,LUPRI)
         ELSE
            CALL OUTPUT(GRDMOL,1,1,1,NCOOR,1,MXCOOR,1,LUPRI)
         END IF
      END IF
C     
C     The gradient difference (gamma) is calculated. The vector
C     has to be scaled, to make it "compatible" with the step
C     vector and the Hessian.
C     
      CALL DZERO(GAMMA,2*MXCOOR)
      DO 20 I = 1, NCOOR
         IF (REDINT) THEN
            GAMMA(I) = GRDINT(I) - GRDOLD(I)
         ELSE
            GAMMA(I) = (GRDMOL(I) - GRDOLD(I))/SCLVEC(I)
         END IF
 20   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Scaled gradient difference',-1)
         CALL OUTPUT(GAMMA,1,1,1,NCOOR,1,2*MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck updold */
      SUBROUTINE UPDOLD(SCLVEC,GRDOLD,HESOLD,HESINT)
C
C     The previous Hessian is copied to MOLHES and GRDOLD is updated.
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <optinf.h>
#include <taymol.h>
#include <nuclei.h>
#include <symmet.h>
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      DIMENSION SCLVEC(MXCOOR), GRDOLD(2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESINT(2*MXCOOR,2*MXCOOR)
C
      IF (REDINT) THEN
         CALL DZERO(HESINT,4*MXCOOR*MXCOOR)
         CALL DZERO(GRDOLD,2*MXCOOR)
         DO 10 J = 1, IINTCR
            DO 15 I = 1, IINTCR
               HESINT(I,J) = HESOLD(I,J)
 15         CONTINUE
            GRDOLD(J) = GRDINT(J)
 10      CONTINUE
      ELSE
         CALL DZERO(HESMOL,MXCOOR*MXCOOR)
         CALL DZERO(GRDOLD,2*MXCOOR)
         DO 20 J = 1, NCART
            DO 25 I = 1, NCART
               HESMOL(I,J) = HESOLD(I,J)*SCLVEC(I)*SCLVEC(J)
 25         CONTINUE
            GRDOLD(J) = GRDMOL(J)
 20      CONTINUE
         DO 30 I = NCART+1, NCRTOT
            HESMOL(I,I) = SCLVEC(I)*SCLVEC(I)
 30      CONTINUE
      END IF
C
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from UPDOLD','*',103)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Unchanged Hessian',-1)
            CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
            CALL HEADER('Gradient at current geometry',-1)
            CALL OUTPUT(GRDOLD,1,1,1,ICRD,1,MCRD,1,LUPRI)
         ELSE
            WRITE(LUPRI,*)
     &           'Hessian will not be updated in this iteration.'
         END IF
      END IF
      KEEPHE = .FALSE.
      RETURN
      END

C  /* Deck upbfgs */
      SUBROUTINE UPBFGS(DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using BFGS method
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(2*MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     ******************************************************************
C     The BFGS formula is:
C
C                                      T                         T
C                     gamma(n)*gamma(n)    B(n)*delta(n)*delta(n) *B(n)
C     B(n+1) = B(n) + ------------------ - ----------------------------
C                             T                       T
C                     gamma(n) *delta(n)      delta(n) *B(n)*delta(n)
C
C     where
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     The terms in the formula is evaluated one by one below.
C     ******************************************************************
C
C     First we calculate (gamma^T*delta)
C
      FAC = D0
      DO 10 I = 1, ICRD
         FAC = FAC + GAMMA(I)*DELTA(I)
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC 
      END IF
C
C     (gamma*gamma^T)/(gamma^T*delta) is calculated
C     
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 20 J = 1, ICRD
         DO 22 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*GAMMA(J)/FAC
 22      CONTINUE
 20   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(gamma*gamma^T)/(gamma^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The first two terms of the BFGS formula is placed in TMPMT2.
C
      CALL DZERO(TMPMT2,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMT2(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Sum of first two terms',-1)
         CALL OUTPUT(TMPMT2,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     We place (B*delta) in TMPVEC and calculate (delta^T*B*delta)
C     
      CALL DZERO(TMPVEC,MCRD)
      FAC = D0
      DO 40 I = 1, ICRD
         TMP = D0
         DO 42 J = 1, ICRD
            TMP = TMP + HESOLD(I,J)*DELTA(J)
 42      CONTINUE
         TMPVEC(I) = TMP
         FAC = FAC + DELTA(I)*TMP
 40   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(H*delta)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         CALL HEADER('(delta^T*B*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC 
      END IF
C
C     Since the Hessian is symmetric, the elements in the row vector
C     (delta^T*B) must be equal to the elements in the column vector
C     (B*delta). Calculation of (B*delta*delta^T*B)/(delta^T*B*delta)
C     is therefore quite simple.
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 50 J = 1, ICRD
         DO 52 I = 1, ICRD
            TMPMAT(I,J) = TMPVEC(I)*TMPVEC(J)/FAC
 52      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(B*delta*delta^T*B)/(delta^T*B*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      DO 60 J = 1, ICRD
         DO 62 I = 1, ICRD
            HESNEW(I,J) = TMPMT2(I,J) - TMPMAT(I,J)
 62      CONTINUE
 60   CONTINUE
      IF (IPRINT .GT. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck updfp */
      SUBROUTINE UPDFP(DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using DFP method
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(2*MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     ******************************************************************
C     The DFP formula is:
C
C                                  T                                   T
C                          delta(n) *B(n)*delta(n)    gamma(n)*gamma(n)
C     B(n+1) = B(n) + (1 + -----------------------) * ------------------
C                                    T                         T
C                            gamma(n) *delta(n)       gamma(n) *delta(n)
C
C                               T                              T
C              gamma(n)*delta(n) *B(n) + B(n)*delta(n)*gamma(n)
C            - -------------------------------------------------
C                                    T
C                            gamma(n) *delta(n)
C
C     where
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     The terms in the formula is evaluated one by one below.
C     ******************************************************************
C
C     First we calculate (gamma^T*delta)
C
      FAC1 = D0
      DO 10 I = 1, ICRD
         FAC1 = FAC1 + GAMMA(I)*DELTA(I)
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC1
      END IF
C
C     (delta^T*B*delta) is calculated
C      
      FAC2 = D0
      DO 20 I = 1, ICRD
         DUMMY = D0
         DO 22 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 22      CONTINUE
         FAC2 = FAC2 + DUMMY*DELTA(I)
 20   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*B*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
      FAC2 = 1.0D0 + FAC2/FAC1
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(1+(delta^T*B*delta)/(gamma^T*delta))',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
C
C     We calculate (gamma*gamma^T)/(gamma^T*delta)
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 25 J = 1, ICRD
         DO 27 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*GAMMA(J)/FAC1
 27      CONTINUE
 25   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(gamma*gamma^T)/(gamma^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The first two terms of the DFP formula is placed in TMPMT2.
C
      CALL DZERO(TMPMT2,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMT2(I,J) = HESOLD(I,J) + FAC2*TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Sum of first two terms',-1)
         CALL OUTPUT(TMPMT2,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Next is (delta^T*B)
C
      CALL DZERO(TMPVEC,MCRD)
      DO 40 J = 1, ICRD
         DO 42 I = 1, ICRD
            TMPVEC(J) = TMPVEC(J) + DELTA(I)*HESOLD(I,J)
 42      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*B)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
C
C     Then (gamma*delta^T*B + B*delta*gamma^T)
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 45 J = 1, ICRD
         DO 47 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*TMPVEC(J) + TMPVEC(I)*GAMMA(J)
 47      CONTINUE
 45   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma*delta^T*B + B*delta*gamma^T)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     And then (gamma*delta^T*B + B*delta*gamma^T)/(gamma^T*delta)
C
      DO 50 J = 1, ICRD
         DO 52 I = 1, ICRD
            TMPMAT(I,J) = TMPMAT(I,J)/FAC1
 52      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(gamma*delta^T*B + B*delta*gamma^T)/' //
     &                                           '(gamma^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      DO 60 J = 1, ICRD
         DO 62 I = 1, ICRD
            HESNEW(I,J) = TMPMT2(I,J) - TMPMAT(I,J)
 62      CONTINUE
 60   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck uppsb */
      SUBROUTINE UPPSB(DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using Powell-symmetric-Broyden
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(2*MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     **********************************************************************
C     The Powell-symmetric-Broyden formula is:
C
C                               T                        T
C                      (delta(n) *delta(n))*T(n)*delta(n)
C     B(n+1) = B(n) +  -----------------------------------
C                                     T          2
C                            (delta(n) *delta(n))
C
C                T                        T      T                            T
C       (delta(n) *delta(n))*delta(n)*T(n) -(T(n) *delta(n))*delta(n)*delta(n)
C     + -----------------------------------------------------------------------
C                                        T          2
C                               (delta(n) *delta(n))
C
C     where
C               T(n)     = gamma(n) - B(n)*delta(n)
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     **********************************************************************
C
C     First we calculate (T = gamma-B*delta)
C
      CALL DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DUMMY = D0
         DO 12 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - DUMMY
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(T = gamma-B*delta)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
C
C     We calculate (delta^T*delta)
C
      FAC1 = D0
      DO 20 I = 1, ICRD
         FAC1 = FAC1 + DELTA(I)*DELTA(I)
 20   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC1
      END IF
C
C     Then (T^T*delta)
C
      FAC2 = D0
      DO 22 I = 1, ICRD
         FAC2 = FAC2 + TMPVEC(I)*DELTA(I)
 22   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(T^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
C
C     and (delta^T*delta)^2
C
      FAC3 = FAC1*FAC1
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*delta)^2',-1)
         WRITE(LUPRI,'(A,G16.6)') 'Value :    ', FAC3
      END IF
C
C     Then we construct the complete term
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMAT(I,J) = (  FAC1*TMPVEC(I)*DELTA(J)
     &                     + FAC1*DELTA(I)*TMPVEC(J)
     &                     - FAC2*DELTA(I)*DELTA(J) )/FAC3
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Second term of formula',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      DO 40 J = 1, ICRD
         DO 45 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 45      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck uprnko */
      SUBROUTINE UPRNKO(DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using rank one method
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(2*MXCOOR), GAMMA(2*MXCOOR)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     ******************************************************************
C     The rank one formula is:
C
C                                                                     T
C                     (gamma(n)-B(n)*delta(n))(gamma(n)-B(n)*delta(n))
C     B(n+1) = B(n) + -------------------------------------------------
C                                                   T
C                           (gamma(n)-B(n)*delta(n)) *delta(n)
C
C     where
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     The terms in the formula is evaluated one by one below.
C     ******************************************************************
C
C     First we calculate (gamma-B*delta)
C
      CALL DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DUMMY = D0
         DO 12 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - DUMMY
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma-B*delta)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
C
C     (gamma-B*delta)(gamma-B*delta)^T is calculated
C     
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 20 J = 1, ICRD
         DO 22 I = 1, ICRD
            TMPMAT(I,J) = TMPVEC(I)*TMPVEC(J)
 22      CONTINUE
 20   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma-B*delta)(gamma-B*delta)^T',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     ((gamma-B*delta)^T*delta) is calculated
C     
      FAC = D0
      DO 25 I = 1, ICRD
         FAC = FAC + TMPVEC(I)*DELTA(I)
 25   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('((gamma-B*delta)^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC 
      END IF
C
C     And we complete the term
C     ((gamma-B*delta)(gamma-B*delta)^T)/((gamma-B*delta)^T*delta)
C
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMAT(I,J) = TMPMAT(I,J)/FAC
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('((gamma-B*delta)(gamma-B*delta)^T)/' //
     &                                 '((gamma-B*delta)^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      CALL DZERO(HESNEW,MCRD*MCRD)
      DO 40 J = 1, ICRD
         DO 42 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 42      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upschl */
      SUBROUTINE UPSCHL(DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPVEC,TMPVC2,RMAT,BMAT,MAXDIM,RESET,GNRM,IPRINT)
C
C     Updates the Hessian using Schlegel's updating scheme
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(2*MXCOOR), STPMAT(25,MCRD)
      DIMENSION GAMMA(2*MXCOOR), GRDMAT(25,MCRD)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION RMAT(MCRD,MCRD), BMAT(MCRD,MCRD)
      LOGICAL RESET
C
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(TMPVC2,MCRD)
      CALL DZERO(RMAT,MCRD*MCRD)
      CALL DZERO(BMAT,MCRD*MCRD)
C
C     First we have to transfer the last step and gradient difference to
C     STPMAT and GRDMAT respectively.
C
      IF (RESET) THEN
         CALL DZERO(GRDMAT,25*MCRD)
         CALL DZERO(STPMAT,25*MCRD)
         INUM = 1
C
C     After 25 iterations we have to discard the first entries.
C
      ELSE IF (STPMAT(25,1) .GT. 1.0D10) THEN
         DO 10 I = 1, 24
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(25,1) = 0.0D0
         INUM = 24
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0D10) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
C
C     We also have to check the dimension of our variational space,
C     if the number of displacement vectors exceeds this number, we
C     have to remove some.
C
         IF (INUM .GT. MAXDIM) THEN
            DO 15 I = 1, INUM -1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
         END IF
      END IF
C
C     Then we update all vectors
C
      SNMLST = 0.0D0
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)

         write(lupri,*) 'SNMLST:', snmlst

         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1D10
C
C     **********************************************************************
C     Schlegel's scheme goes as follows:
C
C                    i-1            t
C     r' = (x -x )-  SUM  r ((x -x ) r ) ;   r = r'/|r'| ;   j = i-1, i-2, ...
C      j     j  i   m=j+1  m   j  i   m       j   j   j
C
C                    t      i-1             t               t
C     b   = [ (g -g ) r  -  SUM  b  ((x -x ) r ) ] / (x -x ) r  ;
C      jk       j  i   k   m=j+1  mk   j  i   m        j  i   j
C
C     b   = b   ;   j <= k = i-1, i-2, ...
C      kj    jk
C
C                                       t
C     B  = B   + SUM (b  - r B   r ) r r
C      i    i-1  j k   jk   j i-1 k   j k
C
C     **********************************************************************
C
C     We start by calculating the orthonormal displacement vectors r(j)
C
      DO 20 II = INUM, 1, -1
         TMPNRM = 0.0D0
         DO 22 I = 1, ICRD
            RMAT(II,I) = STPMAT(II,I)
            TMPNRM = TMPNRM + STPMAT(II,I)*STPMAT(II,I)
            TMPVC2(I) = RMAT(II,I)
 22      CONTINUE
         TMPNRM = SQRT(TMPNRM)
         DO 26 I = II+1, INUM
            TMP = 0.0D0
            DO 28 J = 1, ICRD
               TMP = TMP + STPMAT(II,J)*RMAT(I,J)
 28         CONTINUE
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('(x(j)-x(i))^T * r(m)',-1)
               WRITE(LUPRI,*) TMP
            END IF
            DO 30 J = 1, ICRD
               RMAT(II,J) = RMAT(II,J) - RMAT(I,J)*TMP
               TMPVC2(J) = RMAT(II,J)
 30         CONTINUE
 26      CONTINUE
         IF(IPRINT .GE. IPRDBG) THEN
            CALL HEADER('r(j)''',-1)
            CALL OUTPUT(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         TMP = 0.0D0
         DO 32 J = 1, ICRD
            TMP = TMP + RMAT(II,J)*RMAT(II,J)
 32      CONTINUE
         TMP = SQRT(TMP)
C
C     If a vector is in a space already spanned by other vectors
C     (that is less than 15% of the vector is outside this space),
C     all elements are set equal to a very small number.
C
         IF ((TMP .LT. 0.05D0*TMPNRM) .OR. 
     &        (TMPNRM .GE. 1.0D2*SNMLST) .OR.
     &        (GNRM .GE. 0.1D0)) THEN
            GRDMAT(II,1) = 1.1D10
            DO 35 J = 1, ICRD
               RMAT(II,J) = 1.0D-25
 35         CONTINUE
            TMP = 1.0D0
         END IF
         DO 40 J = 1, ICRD
            RMAT(II,J) = RMAT(II,J)/TMP
            TMPVC2(J)  = RMAT(II,J)
 40      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Original displacement vector',-1)
            CALL OUTPUT(STPMAT,II,II,1,ICRD,25,2*MXCOOR,1,LUPRI)
            CALL HEADER('Orthonormalized vector',-1)
            CALL OUTPUT(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
 20   CONTINUE
C
C     All redundant vectors are removed
C
      DO 45 I = INUM, 1, -1
         IF (GRDMAT(I,1) .GT. 1.0D10) THEN
            DO 47 II = I, INUM
               DO 48 III = 1, ICRD
                  RMAT(II,III)   = RMAT(II+1,III)
                  GRDMAT(II,III) = GRDMAT(II+1,III)
                  STPMAT(II,III) = STPMAT(II+1,III)
 48            CONTINUE
 47         CONTINUE
            INUM = INUM -1
         END IF
 45   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Matrix of orthonormalized vectors',-1)
         CALL OUTPUT(RMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Next task is to calculate the matrix elements b(jk)
C
      DO 50 K = INUM, 1, -1
         DO 52 J = K, 1, -1
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('g(j)-g(i)',-1)
               CALL OUTPUT(GRDMAT,J,J,1,ICRD,25,MCRD,1,LUPRI)
            END IF
            BMAT(J,K) = 0.0D0
            DO 58 I = 1, ICRD
               BMAT(J,K) = BMAT(J,K) + GRDMAT(J,I)*RMAT(K,I)
 58         CONTINUE
C
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('x(j)-x(i)',-1)
               CALL OUTPUT(STPMAT,J,J,1,ICRD,25,MCRD,1,LUPRI)
            END IF
C
            DO 70 M = J+1, INUM
               TMP = 0.0D0
               DO 72 I = 1, ICRD
                  TMP = TMP + STPMAT(J,I)*RMAT(M,I)
 72            CONTINUE
               BMAT(J,K) = BMAT(J,K) - BMAT(M,K)*TMP
 70         CONTINUE
            TMP = 0.0D0
            DO 75 I = 1, ICRD
               TMP = TMP + STPMAT(J,I)*RMAT(J,I)
 75         CONTINUE
            BMAT(J,K) = BMAT(J,K)/TMP
            BMAT(K,J) = BMAT(J,K)
 52      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('B-matrix of coefficients',-1)
         CALL OUTPUT(BMAT,1,INUM,1,INUM,MCRD,MCRD,1,LUPRI)
      END IF
C
C     And finally we are ready to update the Hessian
C
C     The elements b(jk) are replaced by b(jk)-r(j)^T*B*r(k)
C
      DO 80 K = 1, INUM
         DO 82 I = 1, ICRD
            TMPVEC(I) = 0.0D0
            DO 83 II = 1, ICRD
               TMPVEC(I) = TMPVEC(I) + HESOLD(I,II)*RMAT(K,II)
 83         CONTINUE
 82      CONTINUE
         IF(IPRINT .GE. IPRDBG) THEN
            CALL HEADER('B*r(k)^T',-1)
            CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         DO 85 J = 1, INUM
            TMP = 0.0D0
            DO 87 I = 1, ICRD
               TMP = TMP + RMAT(J,I)*TMPVEC(I)
 87         CONTINUE
            BMAT(J,K) = BMAT(J,K) - TMP
 85      CONTINUE
 80   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('b(jk)-r(j)^T*B*r(k)',-1)
         CALL OUTPUT(BMAT,1,INUM,1,INUM,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The old Hessian is moved to HESNEW, so that HESOLD can be used for
C     temporary storage.
C
      CALL DZERO(HESNEW,MCRD*MCRD)
      DO 90 J = 1, ICRD
         DO 92 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J)
 92      CONTINUE
 90   CONTINUE
      CALL DZERO(TMPVEC,MCRD)
C
C     Then we construct the matrices r(j)*r(k)^T
C
      DO 100 J = 1, INUM
         DO 110 K = 1, INUM
            CALL DZERO(HESOLD,4*MXCOOR*MXCOOR)
            DO 114 IJ = 1, ICRD
               DO 116 IK = 1, ICRD
                  HESOLD(IJ,IK) = RMAT(J,IJ)*RMAT(K,IK)
 116           CONTINUE
 114        CONTINUE
            IF (IPRINT .GE. IPRDBG) THEN
               CALL HEADER('r(j)*r(k)^T',-1)
               CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,2*MXCOOR,2*MXCOOR,
     &              1,LUPRI)
            END IF
            DO 120 II = 1, ICRD
               DO 122 JJ = 1, ICRD
                  HESNEW(II,JJ) = HESNEW(II,JJ)+BMAT(J,K)*HESOLD(II,JJ)
 122           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upmult */
      SUBROUTINE UPMULT(DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,ICRD,
     &     MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,MAXDIM,RESET,
     &     GNRM,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(2*MXCOOR), STPMAT(25,MCRD)
      DIMENSION GAMMA(2*MXCOOR), GRDMAT(25,MCRD)
      DIMENSION HESOLD(2*MXCOOR,2*MXCOOR), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION TMPMAT(4*MXCOOR*MXCOOR), TMPMT2(4*MXCOOR*MXCOOR)
      DIMENSION RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET
C
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(TMPVC2,MCRD)
      CALL DZERO(RMAT,MCRD*MCRD)
      CALL DZERO(GMAT,MCRD*MCRD)
C
C     First we have to transfer the last step and gradient difference to
C     STPMAT and GRDMAT respectively.
C
      IF (RESET) THEN
         CALL DZERO(GRDMAT,25*MCRD)
         CALL DZERO(STPMAT,25*MCRD)
         INUM = 1
C
C     After 25 iterations we have to discard the first entries.
C
      ELSE IF (STPMAT(25,1) .GT. 1.0D10) THEN
         DO 10 I = 1, 24
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(25,1) = 0.0D0
         INUM = 24
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0D10) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
C
C     We also have to check the dimension of our variational space,
C     if the number of displacement vectors exceeds this number, we
C     have to remove some.
C
         IF (INUM .GT. MAXDIM) THEN
            DO 15 I = 1, INUM-1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
         END IF
      END IF
C
C     Then we update all vectors
C
      SNMLST = 0.0D0
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1D10
C
C     We start by calculating the orthonormal displacement vectors r(j)
C
      DO 20 II = INUM, 1, -1
         TMPNRM = 0.0D0
         DO 22 I = 1, ICRD
            RMAT(II,I) = STPMAT(II,I)
            TMPNRM = TMPNRM + STPMAT(II,I)*STPMAT(II,I)
            TMPVC2(I) = RMAT(II,I)
 22      CONTINUE
         TMPNRM = SQRT(TMPNRM)
         DO 26 I = II+1, INUM
            TMP = 0.0D0
            DO 28 J = 1, ICRD
               TMP = TMP + STPMAT(II,J)*RMAT(I,J)
 28         CONTINUE
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('(x(j)-x(i))^T * r(m)',-1)
               WRITE(LUPRI,*) TMP
            END IF
            DO 30 J = 1, ICRD
               RMAT(II,J) = RMAT(II,J) - RMAT(I,J)*TMP
               TMPVC2(J) = RMAT(II,J)
 30         CONTINUE
 26      CONTINUE
         IF(IPRINT .GE. IPRDBG) THEN
            CALL HEADER('r(j)''',-1)
            CALL OUTPUT(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         TMP = 0.0D0
         DO 32 J = 1, ICRD
            TMP = TMP + RMAT(II,J)*RMAT(II,J)
 32      CONTINUE
         TMP = SQRT(TMP)
C
C     If a vector is in a space already spanned by other vectors
C     (that is less than 15% of the vector is outside this space),
C     all elements are set equal to a very small number.
C
         IF ((TMP .LT. 0.05D0*TMPNRM) .OR. 
     &        (TMPNRM .GE. 1.0D2*SNMLST) .OR.
     &        (GNRM .GE. 0.1D0)) THEN
            GRDMAT(II,1) = 1.1D10
            DO 35 J = 1, ICRD
               RMAT(II,J) = 1.0D-25
 35         CONTINUE
            TMP = 1.0D0
         END IF
         DO 40 J = 1, ICRD
            RMAT(II,J) = RMAT(II,J)/TMP
            TMPVC2(J)  = RMAT(II,J)
 40      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Original displacement vector',-1)
            CALL OUTPUT(STPMAT,II,II,1,ICRD,25,2*MXCOOR,1,LUPRI)
            CALL HEADER('Orthonormalized vector',-1)
            CALL OUTPUT(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
 20   CONTINUE
C
C     All redundant vectors are removed
C
      DO 45 I = INUM, 1, -1
         IF (GRDMAT(I,1) .GT. 1.0D10) THEN
            DO 47 II = I, INUM
               DO 48 III = 1, ICRD
                  RMAT(II,III)   = RMAT(II+1,III)
                  GRDMAT(II,III) = GRDMAT(II+1,III)
                  STPMAT(II,III) = STPMAT(II+1,III)
 48            CONTINUE
 47         CONTINUE
            INUM = INUM -1
         END IF
 45   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Matrix of orthonormalized vectors',-1)
         CALL OUTPUT(RMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The next step is to transform all the gradients
C
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Original gradient vectors',-1)
         CALL OUTPUT(GRDMAT,1,INUM,1,ICRD,25,2*MXCOOR,1,LUPRI)
      END IF
      DO 50 J = 1, INUM
         DLEN2 = D0
         DO 52 JJ = 1, ICRD
            DLEN2 = DLEN2 + STPMAT(J,JJ)*STPMAT(J,JJ)
 52      CONTINUE
         DO 54 I = 1, INUM
            DPROD = D0
            DO 56 II = 1, ICRD
               DPROD = DPROD + RMAT(I,II)*STPMAT(J,II)
 56         CONTINUE
            FAC = DPROD/DLEN2
            DO 58 K = 1, ICRD
               GMAT(I,K) = GMAT(I,K) + FAC*GRDMAT(J,K)
 58         CONTINUE
 54      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Transformed gradient vectors',-1)
         CALL OUTPUT(GMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we use the UPBFGS subroutine to do the updating,
C     suppressing output.
C
      DO 60 K = 1, INUM
         CALL DZERO(DELTA,2*MXCOOR)
         CALL DZERO(GAMMA,2*MXCOOR)
         DO 62 I = 1, ICRD
             DELTA(I) = -RMAT(K,I)
             GAMMA(I) = -GMAT(K,I)
 62      CONTINUE
         CALL UPBFGS(DELTA,GAMMA,HESOLD,HESNEW,
     &        ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,0)
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Updating the Hessian',-1)
            CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         END IF
         CALL DZERO(HESOLD,4*MXCOOR*MXCOOR)
         DO 65 J = 1, ICRD
            DO 67 I = 1, ICRD
               HESOLD(I,J) = HESNEW(I,J)
 67         CONTINUE
 65      CONTINUE
 60   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck mincgh */
C
C     Calls WLKCGH to do some transformations.
C
      SUBROUTINE MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <nuclei.h>
#include <abainf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
#include <cbiwlk.h>
      LOGICAL TRU, FAL
      DIMENSION WORK(LWORK)
      DIMENSION EGRAD(MXCOOR)
      DIMENSION EHESS(MXCOOR,MXCOOR), ALLHES(NCRTOT*NCRTOT)
      IF (REDINT) RETURN
      CALL DZERO(EGRAD,MXCOOR)
      CALL DZERO(EHESS,MXCOOR*MXCOOR)
      CALL DZERO(ALLHES,NCRTOT*NCRTOT)
      TRU = .TRUE.
      FAL = .FALSE.
      CALL WLKCGH(EGRAD,ALLHES,WORK(1),WORK(4),WORK(4+3*NCART),
     &              LWORK-(4+3*NCART),DOREPW,TRU,FAL,NCRTOT,
     &              NCRTOT*NCRTOT,0,NCART,FAL,IPRINT)
      JI = 1
      DO 10 I = 1, NCART
         DO 20 J = 1, NCART
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 20      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck cntcrd */
C
C     Counts the coordinates of all symmetries.
C
      SUBROUTINE CNTCRD
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <prkoor.h>
#include <trkoor.h>
#include <cbiwlk.h>
      NCRTOT = 0
      NTMAT = 0
      NCART = NCRREP(0,1)
      NPROJ = NTRREP(0)
      DO 5 I = 0, 7
         DOREPW(I) = .FALSE.
 5    CONTINUE
      IF (IPRINT .GT. 5) CALL HEADER('Coordinates in CNTCRD',-1)
      DO 10 I  = 0, MAXREP
         NPRREP(I) = NTRREP(I)
         NCR = NCRREP(I,1)
         NPR = NTRREP(I)
         DOREPW(I) = (NCR .GT. NPR)
         IF (DOREPW(I)) THEN
            NCRTOT = NCRTOT + NCR
            NTMAT = NTMAT + NCR*NPR
         END IF
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,'(/A,I2)')' Symmetry', I+1
            WRITE(LUPRI,'(A,I5)')' Cartesian coordinates:', NCR
            WRITE(LUPRI,'(A,I5)')' External coordinates :', NPR
            WRITE(LUPRI,'(A,I5)')' Internal coordinates :', NCR-NPR
            WRITE(LUPRI,'(A,L5)')' Dorepw               :', DOREPW(I)
         END IF
 10   CONTINUE
      IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A,I5)')
     &  ' Total number of Cartesian coordinates:', NCRTOT
      RETURN
      END

C  /* Deck gttmat */
C
C     Construct, scale and orthogonalize T matrix.
C
      SUBROUTINE GTTMAT(TMAT,TMPMAT,NCR,NPR,IREP,THRLDP)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
      DIMENSION TMAT(NCR,NPR), TMPMAT(MXCOOR)
      CALL DZERO(TMAT,NCR*NPR)
      CALL GETTRO(TMAT,TMPMAT,NCR,NPR,'BOTH','TORTHO',
     &                                         'CT',IREP,IPRINT)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('T matrix in GTTMAT',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
C
C     No scaling is done yet...
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Scaled T matrix in GTTMAT',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      NPR1 = NPR
      CALL ORTVEC(0,NPR1,NCR,THRLDP,TMAT,TMPMAT)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Orthogonalized T matrix in GTTMAT',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      IF (NPR1 .NE. NPR) THEN
         WRITE(LUPRI,'(//,2(A,I1),A,/A)')
     *      ' Number of trarot vectors reduced from ', NPR,
     *      ' to ', NPR1, ' in ORTVEC called from WLKPRJ.',
     *      ' Program cannot proceed .'
         CALL QUIT('Insufficient number of trarot vectors in WLKPRJ.')
      END IF
      RETURN
      END

C  /* Deck projgh */
      SUBROUTINE PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
#include <cbiwlk.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <trkoor.h>
#include <prkoor.h>
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION ALLHES(NCRTOT*NCRTOT), TMAT(NTMAT)
      DIMENSION TMPMAT(MXCOOR,MXCOOR), TMPMT2(MXCOOR)
      DIMENSION PROJOP(MXCOOR,MXCOOR)
      IF (IPRINT .GT. 5) CALL TITLER('Output from PROJGH','*',103)
      IHESS = 1
      ITMAT = 1
      DO 10 IREP = 0, MAXREP
         IF (DOREPW(IREP)) THEN
            NCR = NCRREP(IREP,1)
            NPR = NPRREP(IREP)
            IF (IPRINT .GT. 5) THEN
               WRITE(LUPRI,'(1X,A,I2)') 'Symmetry ', IREP+1
               WRITE(LUPRI,'(1X,A,I2)') 'NCR      ', NCR
               WRITE(LUPRI,'(1X,A,I2)') 'NPR      ', NPR
            END IF
            CALL PRJGH1(EGRAD,ALLHES(IHESS),TMAT(ITMAT),PROJOP,
     &                           TMPMAT,TMPMT2,NCR,NPR,IREP,THRLDP)
            IHESS = IHESS + NCR*NCR
            ITMAT = ITMAT + NPR*NCR
         END IF
 10   CONTINUE
      JI = 1
      DO 20 I = 1, NCART
         DO 30 J = 1, NCART
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 30      CONTINUE
 20   CONTINUE            
      RETURN
      END

C  /* Deck prjgh1 */
C
C     Construct projection operator and use it on gradient and
C     Hessian to remove both rotation and translation.
C
      SUBROUTINE PRJGH1(EGRAD,ALLHES,TMAT,PROJOP,TMPMAT,TMPMT2,
     &                                        NCR,NPR,IREP,THRLDP)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
      DIMENSION EGRAD(NCR), ALLHES(NCR,NCR)
      DIMENSION TMAT(NCR,NPR), PROJOP(NCR,NCR)
      DIMENSION TMPMAT(MXCOOR,MXCOOR), TMPMT2(MXCOOR)
      CALL DZERO(PROJOP,NCR*NCR)
      CALL DZERO(TMPMAT,MXCOOR*MXCOOR)
C
C     Get translation and rotation matrix.
C
      CALL GTTMAT(TMAT,TMPMT2,NCR,NPR,IREP,THRLDP)
C
C     Construct operator.
C
      CALL DUNIT(PROJOP,NCR)
      CALL SMPABT(TMAT,NCR,NPR,NCR,NPR,TMAT,NCR,NPR,NCR,NPR,
     &                                           PROJOP,NCR,NCR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Unprojected gradient in PROJGH',-1)
         CALL OUTPUT(EGRAD,1,1,1,NCR,1,NCR,1,LUPRI)
         CALL HEADER('Unprojected Hessian in PROJGH',-1)
         CALL OUTPUT(ALLHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         CALL HEADER('Projection operator in PROJGH',-1)
         CALL OUTPUT(PROJOP,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
C
C     Do projection.
C
      IF (IREP .EQ. 0) THEN
         CALL MPAB(PROJOP,NCR,NCR,NCR,NCR,EGRAD,NCR,1,NCR,1,
     &                                                   TMPMAT,NCR,1)
         CALL DCOPY(NCR,TMPMAT,1,EGRAD,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('Projected gradient in PROJGH',-1)
            CALL OUTPUT(EGRAD,1,1,1,NCR,1,NCR,1,LUPRI)
         END IF
      END IF
      CALL MPAB(PROJOP,NCR,NCR,NCR,NCR,ALLHES,NCR,NCR,NCR,NCR,
     &                                                 TMPMAT,NCR,NCR)
      CALL MPAB(TMPMAT,NCR,NCR,NCR,NCR,PROJOP,NCR,NCR,NCR,NCR,
     &                                                 ALLHES,NCR,NCR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Projected Hessian in PROJGH',-1)
         CALL OUTPUT(ALLHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck diahes */
C
C     Diagonalize Hessian
C
      SUBROUTINE DIAHES(EGRAD,EHESS,ALLHES,TMAT,THRIND,
     &                                         EVEC,EVCTMP,WORK,LWORK)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <cbiwlk.h>
#include <optinf.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION ALLHES(NCRTOT*NCRTOT), TMAT(NTMAT)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR), EVCTMP(NCRTOT*NCRTOT)
      DIMENSION WORK(LWORK)
      LOGICAL PRJTRO
      PRJTRO = .TRUE.
      CALL DZERO(GRDDIA,2*MXCOOR)
      CALL WLKEIG(EGRAD,ALLHES,EVAL,EVCTMP,GRDDIA,TMAT,THRLDP,
     &            THRIND,WORK,CNDHES,INDHES,LWORK,IZEROG,NZEROG,
     &            DOREPW,NCRTOT,NCRTOT*NCRTOT,NTMAT,PRJTRO,IPRINT)
      INDTOT = ISUM(MAXREP+1,INDHES(0),1)
      JI = 1
C
C     Diagonal hessian and eigenvalues are copied.
C
      DO 10 I = 1, NCART
         DO 20 J = 1, NCART
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 20      CONTINUE
 10   CONTINUE
      JI = 1
      II = 0
      DO 30 ISYM = 0, MAXREP
         DO 40 I = 1, NCRREP(ISYM,1)
            DO 50 J = 1, NCRREP(ISYM,1)
               EVEC(II+J,II+I) = EVCTMP(JI)
               JI = JI + 1
 50         CONTINUE
 40      CONTINUE
         II = II + NCRREP(ISYM,1)
 30   CONTINUE
      RETURN
      END

c  /* Deck uptrad */
      SUBROUTINE UPTRAD(REJGEO)
C
C     Updates the trust radius and checks if step should be rejected.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <optinf.h>
      PARAMETER (DP25=0.25D0, THDE=1.0D-10)
      LOGICAL REJGEO
      ERGDIF = ENERGY - ERGOLD
      REJGEO = .FALSE.
C
C     The ratio between actual and predicted energy change is calculated,
C     and this ratio is then used to update the trust radius.
C
      IF (ABS(ERGPRD) .GT. THDE) THEN
         RATIO = ERGDIF/ERGPRD
         IF (IPRINT .GT. 2) THEN
            CALL HEADER('Energy difference to previous geometry:',-1)
            WRITE(LUPRI,'(/A),/)')
     &      '      Actual           /  Predicted       =    Ratio '
            WRITE(LUPRI,'(5X,F15.10,2(A,F15.10))')
     &            ERGDIF,'  /  ',ERGPRD,'  =  ',RATIO
         END IF
      ELSE
         RATIO = 1.0D0
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(3(/A),/5X,1P,2D16.6)')
     &          ' Close to convergence, ratio set to one.',
     &          ' Energy difference to previous geometry:',
     &          ' actual and predicted:', ERGDIF, ERGPRD
         END IF
      END IF
      IF ((RATIO .LT. RTRJMN) .OR. (RATIO .GT. RTRJMX)) THEN
         WRITE(LUPRI,'(/A/A,F12.6)')
     &           ' Step rejected because ratio between predicted',
     &           ' and actual energy change is : ', RATIO
         REJGEO = .TRUE.
      ELSE IF (RATIO .LT. RTENBD) THEN
         WRITE(LUPRI,*) 'Trust radius decreased due to bad ratio.'
         TRSTRA = TRSTDE*STPNRM
      ELSE IF (ABS(RATIO-1.0D0) .LT. DP25*(1.0D0-RTENGD)) THEN
         WRITE(LUPRI,*)
     &       'Trust radius squarely increased due to very good ratio.'
         TRSTRA = MAX(TRSTIN*TRSTIN*STPNRM,TRSTRA)
         IF (ISTATE .GT. 1 .AND. TRSTRA .GT. 0.30D0) TRSTRA = 0.30D0
      ELSE IF (RATIO .GT. RTENGD) THEN
         WRITE(LUPRI,*) 'Trust radius increased due to good ratio.'
         TRSTRA = MAX(TRSTIN*STPNRM,TRSTRA)
         IF (ISTATE .GT. 1 .AND. TRSTRA .GT. 0.30D0) TRSTRA = 0.30D0
      ELSE
         WRITE(LUPRI,*) 'Trust radius set equal to norm of step.'
         TRSTRA = STPNRM
         IF (ISTATE .GT. 1 .AND. TRSTRA .GT. 0.30D0) TRSTRA = 0.30D0
      END IF
      IF (.NOT. REJGEO)
     &     WRITE(LUPRI,'(A,F10.5)') ' Updated trust radius', TRSTRA
      RETURN
      END

C  /* Deck fndstp */
      SUBROUTINE FNDSTP(EVEC,TMPMAT,EVCTMP,CSTEP)
C
C     This routine calculates the step that should be taken to obtain
C     the next geometry.
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
      PARAMETER (D0 = 0.0D0 , DP5 = 0.5D0)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR),EVCTMP(NCART,NCART)
      DIMENSION CSTEP(MXCOOR)
      EXTERNAL WLKSTL
      CALL DZERO(STPDIA,NCART)
      CALL DZERO(STPSYM,NCART)
      NVEC = NCART-NPROJ
      IF (IPRINT .GT. 2) CALL TITLER('Output from FNDSTP','*',103)
C
C     We take a copy of the eigenvectors, to get correct
C     matrix dimensions.
C
      DO 10 I = 1, NCART
         DO 20 J = 1, NCART
            EVCTMP(J,I) = EVEC(J,I)
 20      CONTINUE
 10   CONTINUE
C
C     Newton step is calculated.
C
      DO 30 I = 1, NVEC
         STPDIA(I) = -GRDDIA(I)/EVAL(I)
 30   CONTINUE
      STPNRM = SQRT(DDOT(NVEC,STPDIA,1,STPDIA,1))
      GRADNM = SQRT(DDOT(NVEC,GRDDIA,1,GRDDIA,1))
C
C     If Newton step is larger that trust radius, we take a step
C     to the boundary.
C
      IF ((STPNRM .GT. TRSTRA) .AND. (.NOT. NOTRST)) THEN
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,'(/A,F15.10)')' Norm of Newton step:', STPNRM
            WRITE(LUPRI,'(A,F15.10/)')' Trust radius       :', TRSTRA
         END IF
C
C     Simple rip-off of level-shift from ABAWLK.
C
         GRDMIN = ABS(GRDDIA(1))
         HESMIN = EVAL(1)
         XMIN = GRADNM/TRSTRA - MIN(HESMIN,D0)
         XMAX = MAX(D0, -HESMIN + DP5*GRDMIN/TRSTRA)
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,'(/A,2(1P,D12.5,2X))')
     &          ' XMIN AND XMAX BEFORE WLKBIS (WLKSTL):  ',XMIN,XMAX
         END IF
         CALL FLSHFO(LUPRI)
         CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA,EVAL,TRSTRA,NVEC,
     &                                                  WLKSTL,IFAIL)
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,'(A,2(1P,D12.5,2X))')
     &          ' XMIN and XMAX after WLKBIS (WLKSTL):   ',XMIN,XMAX
         END IF
         IF (IPRINT .GT. 3) THEN
            WRITE(LUPRI,'(A,1P,D12.5)')' Level shift parameter: ',RNU
         END IF
         IF (IFAIL .EQ. 0) THEN
            WRITE(LUPRI,'(/A/,2(A,F10.6))')
     &        '  *** Wrong interval specified in WLKBIS (WLKSTL) ***',
     &        ' XMAX = ', XMAX, ' XMIN = ', XMIN
         ELSE IF (IFAIL .EQ. 1) THEN
            WRITE(LUPRI,'(3(/A))')
     &      '  *** WARNING WLKBIS (WLKSTL) ***',
     &      ' Desired accuracy not obtained in the specified maximum',
     &      ' number of iterations.'
         END IF
         DO 37 I = 1, NVEC
            STPDIA(I) = -GRDDIA(I)/(EVAL(I) + RNU)
 37      CONTINUE
         STPNRM = SQRT(DDOT(NVEC,STPDIA,1,STPDIA,1))
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,'(/A,F15.10)')' Norm, boundary step:', STPNRM
         END IF
      END IF
C
C     Energy is predicted, will be used later to update trust radius.
C
      ERGPRD = DDOT(NVEC,GRDDIA,1,STPDIA,1)
     &         + 0.5D0*DV3DOT(NVEC,STPDIA,EVAL,STPDIA)
      WRITE(LUPRI,'(/A,F25.15)') ' Predicted energy change',ERGPRD
C
C     If the predicted energy is positive, it means the Newton step is
C     towards a maximum/saddle point. We then simply reverse the
C     step direction (a bit dirty, but seems to work).
C
      IF (ERGPRD .GT. 0.0D0) THEN
         WRITE(LUPRI,*) 'Reversing step!'
         DO 40 I = 1, NVEC
            STPDIA(I) = -STPDIA(I)
 40      CONTINUE
         ERGPRD = DDOT(NVEC,GRDDIA,1,STPDIA,1)
     &        + 0.5D0*DV3DOT(NVEC,STPDIA,EVAL,STPDIA)
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(A,F25.15)')
     &           ' New pred. energy change',ERGPRD
         END IF
      END IF
      CALL HEADER('Step in diagonal representation',1)
      CALL OUTPUT(STPDIA,1,1,1,NVEC,1,2*MXCOOR,1,LUPRI)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvector basis',1)
         CALL OUTPUT(EVCTMP,1,NCART,1,NVEC,NCART,NCART,1,LUPRI)
      END IF
      DO 50 I = 1, NVEC
         CALL DAXPY(NCART,STPDIA(I),EVEC(1,I),1,STPSYM,1)
 50   CONTINUE
      STPNRM = DDOT(NCART,STPSYM,1,STPSYM,1)
      STPNRM = SQRT(STPNRM)
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Cartesian symmetry step vector',-1)
         CALL OUTPUT(STPSYM,1,1,1,NCART,1,NCART,1,LUPRI)
         WRITE(LUPRI,'(/A,F15.10/)') ' Norm of step:', STPNRM
      END IF
      CALL WLKCOR(STPSYM,CSTEP,NCART,MXCOOR,IPRINT)
      RETURN
      END

C  /* Deck fndgeo */
      SUBROUTINE FNDGEO(CSTEP,EGRAD,COONEW,COOOLD,EXHER,EXSIR,EXABA,
     &                       WORK,LWORK,IWOFF,WRKDLM,IREJ,GEINFO)
C
C     If the step is acceptable, the geometry is updated and written to file.
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
      PARAMETER (LUMOL = 9)
      DIMENSION CSTEP(MXCOOR), EGRAD(MXCOOR)
      DIMENSION COONEW(3,MXCENT), COOOLD(3,MXCENT)
      DIMENSION WORK(LWORK)
      DIMENSION GEINFO(0:ITRMAX,6)
      CHARACTER*10 FILENM
      LOGICAL EXHER,EXSIR,EXABA,REJGEO
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      REJGEO = .TRUE.
      ERGOLD = ENERGY
      IREJ   = 0
      IJ = 1
      DO 10 J = 1, NUCIND
         DO 20 I = 1, 3
            COOOLD(I,J) = CORD(I,J)
            COONEW(I,J) = CORD(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 20      CONTINUE
 10   CONTINUE
C
C     Here we start a loop to obtain acceptable step, note that REJGEO
C     is initially set TRUE to enter the loop.
C
 50   CONTINUE
      IF ((IREJ .LE. MAXREJ) .AND. REJGEO) THEN
         IF (IPRINT .GT. 2) THEN
            CALL HEADER('New geometry',-1)
            CALL PRIGEO(COONEW)
         END IF
         CALL WLKMOL(COONEW)
C     
C     Write updated files.
C
         CALL PNCMOL
C
C     Memory trap has to be reinitialized.
C
         WORK(IWOFF+1) = WRKDLM
C
C     Calculate energy at new geometry, which is compared to predicted
C     energy(change) in UPTRAD.
C
         CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
         IF (IREJ .EQ. 0) THEN
            CALL UPTRAD(REJGEO)
C
C     After the first failure, we are satisfied if the new energy is below
C     the last. No comparison with predicted energy is done.
C
         ELSE IF (ENERGY .GT. ERGOLD) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
     &           'Step rejected because energy is increasing.'
            WRITE(LUPRI,'(A,F10.5)')' Updated trust radius', TRSTRA
            REJGEO = .TRUE.
         ELSE
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
     &           'Acceptable step has been found.'
            REJGEO = .FALSE.
         END IF
         IF (REJGEO) THEN
            IREJ = IREJ + 1
            EXHER  = .FALSE.
            EXSIR  = .FALSE.
            RDINPC = .FALSE.
            RDMLIN = .FALSE.
            HRINPC = .FALSE.
C     
C     Line search based on quadratic model
C
            GRADDI = 0.0D0
            DO 60 I = 1, NCART
               GRADDI = GRADDI + GRDDIA(I)*STPDIA(I)
 60         CONTINUE
            GRADDI = GRADDI/STPNRM
C
            IF (IPRINT .GT. 12) THEN
               CALL HEADER('Line search based on quadratic model',-1)
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Energy at last geometry     : ', ERGOLD
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Energy at rejected geometry : ', ENERGY
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Norm of rejected step       : ', STPNRM
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Norm of gradient            : ', GRADNM
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Gradient along step         : ', GRADDI
            END IF
C
C     The minimum for a quadratic model is calculated with the formula
C                    -f'(0)
C     x     =  -------------------
C      min     2*(f(1)-f(0)-f'(0))
C
            FAC = -0.5D0*GRADDI/(ENERGY-ERGOLD-GRADDI)
C
C     If the factor found is very small or very large, we don't trust
C     it. The factor is replaced by "safer" (but rather atbitrary) numbers.
C
            IF (FAC .LT. 0.1D0) FAC = 0.25D0
            IF (FAC .GT. 0.9D0) FAC = 0.75D0
C
C     We have to update both steps and their norm.
C
            DO 70 I = 1, IINTCR
               STPINT(I) = STPINT(I)*FAC
 70         CONTINUE
            DO 75 I = 1, NCART
               STPDIA(I) = STPDIA(I)*FAC
               STPSYM(I) = STPSYM(I)*FAC
 75         CONTINUE
            STPNRM = STPNRM*FAC
C
C     We also set the trust radius equal to the new norm
C
            TRSTRA = STPNRM
C
            WRITE(LUPRI,'(A,F12.6)')
     &           ' Minimum for quadratic model : ', FAC
            WRITE(LUPRI,'(A,F12.6)')
     &           ' Norm of new step            : ', STPNRM
C
C     Finally we construct a new geometry based on the factor found
C
            DO 80 J = 1, NUCIND
               DO 85 I = 1, 3
                  COONEW(I,J)=FAC*COONEW(I,J)+(1.0D0-FAC)*COOOLD(I,J)
 85            CONTINUE
 80         CONTINUE
         END IF
         GOTO 50
      ELSE IF (REJGEO) THEN
C
C     Maximum number of allowed rejections reached
C
         GEINFO(ITRNMR,4) = STPNRM
         IF (ITRNMR .LT. ITRMAX) GEINFO(ITRNMR+1,5) = TRSTRA
         GEINFO(ITRNMR,6) = IREJ*1.0D0
         CALL PRIINF(GEINFO)
C        Print angles, bondlength etc.
         CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,.FALSE.,
     &               WORK(IWOFF+1),LWORK)
         WRITE(LUPRI,*) 'Maximum number of rejected steps (',MAXREJ,
     &                    ') reached.'
         WRITE(LUPRI,*) 'No accepatble step found. Aborting.'
         CALL QUIT('*** FNDGEO *** No acceptable step found.')
      END IF
      RETURN
      END

C  /* Deck dobrk */
      SUBROUTINE DOBRK(EXHER,EXSIR,EXABA,GEINFO,EVEC,
     &                                        CSTEP,TMPMAT,NONT)
C
C     This procedure does the actual breaking of symmetry.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <molinp.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
#include <cbiwlk.h>
#include <cbirea.h>
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
#include <huckel.h>
      PARAMETER (LUMOL = 9)
      LOGICAL EXHER, EXSIR, EXABA, BRKALS
      CHARACTER FILENM*10, TMPLN*80
      DIMENSION GEINFO(0:ITRMAX,6), NONT(MXCENT)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR), CSTEP(MXCOOR)
      DIMENSION TMPMAT(NCRTOT*NCRTOT)
#include <ibtfun.h>
      CALL DZERO(STPDIA,2*MXCOOR)
      CALL DZERO(CSTEP,MXCOOR)
      ISYMBR = MAXREP
      START = .TRUE.
      BRKALS = .FALSE.
      CALL IZERO(NUCNUM, MXCENT*8)
C
C     We find the highest symmetry number with a non-zero index.
C     Some test runs indicated this as the most effective way to
C     decrease the energy.
C
 10   CONTINUE
      IF ((INDHES(ISYMBR) .LT. 1) .OR. (.NOT. DOREPW(ISYMBR))) THEN
         ISYMBR = ISYMBR - 1
         IF (ISYMBR .LT. 0) THEN
C
C     If breaking the symmetry causes problems, we remove all symmetry
C
            IF (MAXOPR .GT. 0) THEN
               BRKALS = .TRUE.
               GOTO 77
C
C     If all symmetry has been removed and we still have trouble, it's
C     time to give up!
C
            ELSE
               CALL QUIT
     &         ('*** DOBRK *** Breaking of symmetry was unsuccesful.')
            END IF
         END IF
         GOTO 10
      END IF
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)'***** NOTE! *****'
      WRITE(LUPRI,*)
     &     'Due to non-zero index of total Hessian, symmetry #',ISYMBR
      WRITE(LUPRI,*)'has to be broken to minimize energy!'
      WRITE(LUPRI,*)
C
C     The eigenvectors of the symmetry to be broken, are copied to TMPMAT.
C
      II = 0
      DO 15 I = 0, ISYMBR - 1
         II = II + NCRREP(I,1)
 15   CONTINUE
      JI = 1
      DO 20 I = 1, NCRREP(ISYMBR,1)
         DO 25 J = 1, NCRREP(ISYMBR,1)
            TMPMAT(JI) = EVEC(II+J,II+I)
            JI = JI + 1
 25      CONTINUE
 20   CONTINUE
      NCR = NCRREP(ISYMBR,1)
      NVC = NCR - NTRREP(ISYMBR)
      DO 27 I = 1, NVC
         STPDIA(I) = 1.0D0/SQRT(NVC*1.0D0)
 27   CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Geometry before break of symmetry',1)
         CALL PRIGEO(CORD)
         CALL HEADER('Step in diagonal representation',1)
         CALL OUTPUT(STPDIA,1,1,1,NVC,1,2*MXCOOR,1,LUPRI)
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvector basis',1)
         CALL OUTPUT(TMPMAT,1,NCR,1,NVC,NCRTOT,NCRTOT,1,LUPRI)
      END IF
C
C     We determine the cartesian symmetry-breaking step.
C
      DO 30 I = 1, NVC
         CALL DAXPY(NCR,STPDIA(I),TMPMAT((I-1)*NCR+1),1,CSTEP,1)
 30   CONTINUE
      CALL DZERO(TMPMAT, NCRTOT*NCRTOT)
      DO 32 I = 1, NCR
         TMPMAT(I) = CSTEP(I)*1.00D0
 32   CONTINUE
      STPNRM = SQRT(DDOT(NCR,TMPMAT,1,TMPMAT,1))
      CALL DZERO(CSTEP,MXCOOR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Cartesian symmetry step vector',1)
         CALL OUTPUT(TMPMAT,1,1,1,NCR,1,NCR,1,LUPRI)
         WRITE(LUPRI,'(/A,F15.10/)') ' Norm of step:', STPNRM
      END IF
      DO 40 IATOM = 1, NUCIND
         DO 42 ICOOR = 1, 3
            ICCOOR = 3*(IATOM - 1) + ICOOR
            ISCOOR = IPTCNT(ICCOOR,ISYMBR,1)
            IF (ISCOOR .GT. 0) THEN
               CSTEP(ICCOOR)=TMPMAT(ICCOOR)/SQRT(FMULT(ISTBNU(IATOM)))
            END IF
 42      CONTINUE
 40   CONTINUE
C
C     Occasionally the symmetry breaking step is zero, and another
C     symmetry must be chosen for breaking.
C
      TMP = 0.0D0
      DO 31 I = 1, NCR
         TMP = TMP + CSTEP(I)*CSTEP(I)
 31   CONTINUE
      IF (TMP .LE. ZERGRD) THEN
         WRITE(LUPRI,*)
     &          'Zero step vector found, breaking another symmetry'
         ISYMBR = ISYMBR - 1
         GOTO 10
      END IF
C
C     The new geometry is calculated.
C
 77   CONTINUE
      IJ = 1
      DO 45 J = 1, NUCIND
         DO 47 I = 1, 3
            CORD(I,J) = CORD(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 47      CONTINUE
 45   CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Cartesian step vector in non-symmetry basis',1)
         CALL PRIGEO(CSTEP)
         CALL HEADER('Geometry after break of symmetry',1)
         CALL PRIGEO(CORD)
      END IF
      CALL WLKMOL(CORD)
C
C     We remove the number of symmetry operations, so that new symmetry
C     will be added in the next iteration. If we have broken the symmetry
C     earlier but ends up with the same energy, we're having trouble.
C     This is brutally resolved by removing all symmetry, that is the
C     molecule has to be minimized withim the C1 point group.
C
      TMPLN = MLINE(1)
      BASIS = (TMPLN(1:5) .EQ. 'BASIS')
      I = 4
      IF (BASIS) I = 5
      TMPLN = MLINE(I)
      READ(TMPLN(1:5),'(I5)') NONTYP 
      WRITE(TMPLN(10:20), '(A11)') '           '
      II = 0
      DO 49 J = 0, ITRNMR-1
         IF ((GEINFO(J,3) .LT. -0.5D0) .AND. 
     &        (ABS(ENERGY-GEINFO(J,1)) .LT. THRERG*10)) II=II+1
 49   CONTINUE
      IF (II .GT. 0) BRKALS = .TRUE.
      IF (BRKALS) WRITE(TMPLN(10:20), '(A11)') '0          '
      MLINE(I) = TMPLN
C
      IATOM = 1
      I = NCLINE(IATOM)-1
      DO 50 J = 1, NONTYP
         READ(MLINE(I),'(BN,6X,F4.1,I5)') Q, NONT(J)
         IATOM = IATOM + NONT(J)
         I = NCLINE(IATOM)-1
 50   CONTINUE
C
C     We run over all atom types and all symmetry independant centres and
C     and expand them to all atoms. The molecule input is modified
C     according to this.
C
      IATOM = 1
      DO 60 ITYP = 1, NONTYP
         IXTRA = 0
         DO 70 I = 1, NONT(ITYP)
            MULCNT = ISTBNU(IATOM)
            DO 80 ISYMOP = 1, MAXOPR
               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  COOX=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  COOY=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  COOZ=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
                  INSLIN = NCLINE(IATOM)
                  DO 90 J = NMLINE, INSLIN, -1
                     MLINE(J+1) = MLINE(J)
 90               CONTINUE
                  DO 95 J = (IATOM+1), NUCIND
                     NCLINE(J) = NCLINE(J) + 1
 95               CONTINUE
                  TMPLN = MLINE(INSLIN + 1)
                  WRITE(TMPLN(5:64),'(3F20.16)') COOX,COOY,COOZ
                  MLINE(INSLIN + 1) = TMPLN
                  NMLINE = NMLINE + 1
                  IXTRA = IXTRA + 1
               END IF
 80         CONTINUE
            IATOM = IATOM + 1
 70      CONTINUE
         TMPLN = MLINE(NCLINE(IATOM - NONT(ITYP)) - 1)
         WRITE(TMPLN(11:15),'(I5)') NONT(ITYP) + IXTRA
         MLINE(NCLINE(IATOM - NONT(ITYP)) - 1) = TMPLN
 60   CONTINUE
C     
C     Write updated geometry to files.
C
      CALL PNCMOL
C
C     To mark that symmetry was broken, the index of the hessian is
C     given a negative sign. This is interpreted in PRIINF.
C
      GEINFO(ITRNMR,3) = -ABS(GEINFO(ITRNMR,3))
C
C     Several variables and arrays has to be modified/reset to be able
C     to continue calculation with new symmetry.
C
      BRKSYM = .FALSE.
      ITRBRK = ITRNMR
      GECONV = .FALSE.
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      NEWSYM = .TRUE.
      ADDSTO = .TRUE.
      KEEPHE = .FALSE.
      ERGOLD =  ENERGY
      CALL IZERO(NUCNUM, MXCENT*8)
      CALL IZERO(NCRREP, 16)
      CALL IZERO(IPTCNT, MXCENT*48)
      CALL IZERO(NAXREP, 16)
      CALL IZERO(INDHES, 8)
      RETURN
      END

C  /* Deck visulz */
      SUBROUTINE VISULZ(WORK,LWORK,IWOFF,WRKDLM,TMPMAT,EVEC)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <optinf.h>
#include <cbirea.h>
#include <molinp.h>
#include <abainf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <hrunit.h>
      CHARACTER TMPLIN*80
      DIMENSION WORK(LWORK), TMPMAT(MXCOOR,MXCOOR)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR)
C
C     HERINP has to be run to process geometry input.
C
      LUONEL = 16
      CALL HERINP(WORK,LWORK)
      CLOSE(24)
      WORK(IWOFF+1) = WRKDLM
C     
C     Make VRML-file of geometry
C     
      CALL MKVRML(.FALSE.,TMPMAT,EVEC)
      RETURN
      END

C  /* Deck inipre */
      SUBROUTINE INIPRE(WORK,LWORK,IWOFF,WRKDLM,TMPMAT,EVEC)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <optinf.h>
#include <cbirea.h>
#include <molinp.h>
#include <abainf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <hrunit.h>
      CHARACTER TMPLIN*80
      DIMENSION WORK(LWORK), TMPMAT(MXCOOR,MXCOOR)
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR)
C
C     When preoptimization is scheduled, we have to change the basis set.
C     HERINP has to be run before we can do the first change.
C
      LUONEL = 16
      CALL HERINP(WORK,LWORK)
      CLOSE(LUONEL)
      CLOSE(24)
      WORK(IWOFF+1) = WRKDLM
C
C     Basis set library has to be used for preoptimization to work
C
      IF (BASIS) THEN
         TMPLIN = MLINE(2)
         MLINE(2) = PREBTX(1)
         PREBTX(1) = TMPLIN
         ITRNMR = ITRNMR - 1
         CALL PNCMOL
         ITRNMR = ITRNMR + 1
      ELSE
         DOPRE = .FALSE.
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '*** WARNING! ***'
         WRITE(LUPRI,*) 'Preoptimization can only be done when ' //
     &        'basis set library is used.'
         WRITE(LUPRI,*) 'Preoptimization has been turned off!'
         WRITE(LUPRI,*)
      END IF
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      IPRE = IPRE + 1
C
C     Make VRML-file of initial geometry if requested
C
      IF (VRML) CALL MKVRML(.FALSE.,TMPMAT,EVEC)
      RETURN
      END

C  /* Deck endpre */
      SUBROUTINE ENDPRE(EXHER,EXSIR,EXABA)
C
C     This procedure ends preoptimization and starts optimization with
C     the "main" basis set.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <molinp.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
#include <cbiwlk.h>
#include <cbirea.h>
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
      CHARACTER TMPLIN*80
      LOGICAL EXHER, EXSIR, EXABA
#include <ibtfun.h>
      TMPLIN = MLINE(2)
      IF (IPRE .LT. NUMPRE) THEN
         MLINE(2) = PREBTX(IPRE+1)
         PREBTX(IPRE+1) = PREBTX(IPRE)
         PREBTX(IPRE) = TMPLIN
      ELSE
         MLINE(2) = PREBTX(IPRE)
         PREBTX(IPRE) = TMPLIN
         FINPRE = .TRUE.
      END IF
C
C     If all symmetry has been removed during preoptimization, we try
C     to detect it again after each preoptimization phase.
C
      TMPLIN = MLINE(5)
      IF ((MAXOPR .EQ. 0) .AND. (ITRBRK .GE. 0)) THEN
         WRITE(TMPLIN(10:20), '(A11)') '           '
         MLINE(5) = TMPLIN
         NWSYMM = .TRUE.
      END IF
      CALL PNCMOL
      IPRE = IPRE + 1
      GECONV = .FALSE.
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      NEWSYM = .TRUE.
      KEEPHE = .TRUE.
      ERGOLD =  ENERGY
      CALL IZERO(NUCNUM, MXCENT*8)
      CALL IZERO(NCRREP, 16)
      CALL IZERO(IPTCNT, MXCENT*48)
      CALL IZERO(NAXREP, 16)
      CALL IZERO(INDHES, 8)
      RETURN
      END

C  /* Deck mkscvc */
      SUBROUTINE MKSCVC(SCLVEC)
C
C     Makes scaling vector that's necessary for manipulating
C     gradient and step vectors.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <optinf.h>
#include <nuclei.h>
#include <symmet.h>
      DIMENSION SCLVEC(MXCOOR)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We initialize the scaling vector (the inverse of the
C     normalization vector in WLKCGH).
C
      CALL DZERO(SCLVEC,MXCOOR)
      DO 1 IREP = 0, MAXREP
         DO 2 ICENT = 1, NUCIND
#include <scalar.h>
            DO 3 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  SCLVEC(ISCOOR) = SQRT(FMULT(ISTBNU(ICENT)))
               END IF
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from MKSCVC','*',103)
         CALL HEADER('Scaling vector',-1)
         CALL OUTPUT(SCLVEC,1,1,1,NCRTOT,1,MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck minend */
      LOGICAL FUNCTION MINEND(SCLVEC)
C
C     Determines if the end of the minimization has been reached.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <taymol.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
      DIMENSION SCLVEC(MXCOOR)
      LOGICAL CNVERG, CNVGRD, CNVSTP, INDXOK
      CHARACTER*4 LOGTXT
      CNVERG = .FALSE.
      CNVGRD = .FALSE.
      CNVSTP = .FALSE.
      ICOOR = NCART
      IF (REDINT) ICOOR = IINTCR
C
C     First we check the different convergence criterias.
C
      IF (.NOT. BAKER) THEN
         ICONV = 0
         EDIFF = ABS(ERGOLD-ENERGY)
         IF (EDIFF .LE. THRERG) THEN
            CNVERG = .TRUE.
            ICONV = ICONV + 1
         END IF
         IF (GRADNM .LE. GRDTHR) THEN
            CNVGRD = .TRUE.
            ICONV = ICONV + 1
         END IF
         IF (STPNRM .LE. THRSTP) THEN
            CNVSTP = .TRUE.
            ICONV = ICONV + 1
         END IF
         INDXOK = (INDTOT .EQ. 0)
         MINEND = INDXOK .AND. (ICONV .GE. ICONDI)
         IF (NOBRKS) MINEND = (ICONV .GE. ICONDI)
C
C     Alternatively we use the convergence criteria of Baker
C
      ELSE
         IF (REDINT) THEN
            CALL MAXELM(GRDINT,IINTCR,SCLVEC,0,GRDMAX)
            CALL MAXELM(STPINT,IINTCR,SCLVEC,0,STPMAX)
         ELSE
            CALL MAXELM(GRDMOL,ICOOR,SCLVEC,2,GRDMAX)
            CALL MAXELM(STPSYM,ICOOR,SCLVEC,1,STPMAX)
         END IF
         EDIFF  = ABS(ERGOLD-ENERGY)
         INDXOK = (INDTOT .EQ. 0)
C
         CNVGRD = (GRDMAX .LE. 3.0D-4)
         CNVERG = (EDIFF  .LE. 1.0D-6)
         IF ((ITRNMR .LT. 1) .OR. (ITRBRK .EQ. (ITRNMR-1)))
     &          CNVERG = .FALSE.
         CNVSTP = (STPMAX .LE. 3.0D-4)
         MINEND = (INDXOK .AND. CNVGRD) .AND. (CNVERG .OR. CNVSTP)
         IF (NOBRKS) MINEND =  (CNVGRD .AND. (CNVERG .OR. CNVSTP))
      END IF
C
C     If the energy difference is below THRSYM and the total Hessian is
C     non-zero, we can be pretty sure that we're approaching a saddle
C     point, and the symmetry should be broken. Only applies to
C     second order methods.
C
      IF ((NEWTON .OR. QUADSD) .AND. (INDTOT .GT. 0) .AND.
     &     (EDIFF .LE. MAX(THRERG,THRSYM)) .AND.
     &     ((ITRNMR-1) .NE. ITRBRK)) THEN
         BRKSYM = .TRUE.
         NWSYMM = .TRUE.
      END IF
C
C     Output from the testing is written.
C
      IF (IPRINT .GT. 2) THEN
         CALL TITLER('Output from MINEND','*',103)
         LOGTXT = 'no  '
         IF (CNVERG) LOGTXT = 'yes '
         IF ((ITRNMR .LT. 1) .OR. (ITRBRK .EQ. (ITRNMR-1)))
     &        LOGTXT = 'N/A '
         WRITE(LUPRI,'(/A,1P,A5)') ' Energy converged      ',LOGTXT
         LOGTXT = 'no  '
         IF (CNVGRD) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' Gradient converged    ',LOGTXT
         LOGTXT = 'no  '
         IF (CNVSTP) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' Step converged        ',LOGTXT
         IF (.NOT. BAKER) THEN
            WRITE(LUPRI,'(A,1P,I3)') ' Conditions fullfilled ',ICONV
            WRITE(LUPRI,'(A,1P,I3)') ' Required conditions   ',ICONDI
         END IF
         WRITE(LUPRI,'(A,1P,I3)') ' Totally sym. index    ',INDHES(0)
         WRITE(LUPRI,'(A,1P,I3)') ' Hessian index         ',INDTOT
         LOGTXT = 'no  '
         IF (MINEND) LOGTXT = 'yes '
         IF (BRKSYM) LOGTXT = 'yes*'
         WRITE(LUPRI,'(A,1P,A5)') ' End of minimization   ',LOGTXT
         IF (BRKSYM) WRITE(LUPRI,'(/A)') ' *) Within given symmetry.'
      END IF
      RETURN
      END

C  /* Deck spnrgy */
      SUBROUTINE SPNRGY(GEINFO,EXHER,EXSIR,EXABA,WORK,LWORK,
     &                                                IWOFF,WRKDLM)
C
C     Calculates single point energy with another basis if requested.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
      DIMENSION GEINFO(0:ITRMAX,6)
      LOGICAL EXHER,EXSIR,EXABA,REJGEO
      DIMENSION WORK(LWORK)
      EXHER = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      MLINE(2) = SPBSTX
C
C     New input files are written with desired basis. Iteration number
C     has to be decreased first in order to write to the correct filename.
C
      CALL PNCMOL
C
C     Initialization of new memory trap for programs.
C
      WORK(IWOFF+1) = WRKDLM
C
      NEWSYM = .TRUE.
      CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
      GEINFO(ITRNMR+1,1) = ENERGY
      RETURN
      END
 
C  /* Deck pncmol */
      SUBROUTINE PNCMOL(TMPMAT)
C     
C     Punch MOLECULE input with updated coordinates to MOLECULE.INP
C     and XXXX_mol.inp + optional VRML-file
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
      PARAMETER (LUMOL = 9)
      DIMENSION TMPMAT(MXCOOR)
      CHARACTER*12 FILENM
      CLOSE(LUMOL)
#if    defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) \
    || defined (SYS_PARAGON) || defined (SYS_CRAY)   || defined (SYS_IRIX) \
    || defined (SYS_SUN)     || defined (SYS_LINUX)  || defined (SYS_HPUX) 
#define SYS_MANY
#endif
#if defined (SYS_VAX)
      OPEN(LUMOL,STATUS='NEW',FORM='FORMATTED',
     &                                 CARRIAGECONTROL='LIST')
#endif
#ifdef SYS_MANY
      OPEN(LUMOL,FILE='MOLECULE.INP',STATUS='UNKNOWN',
     &                                          FORM='FORMATTED')
#endif
#if defined (SYS_IBM)
      OPEN(LUMOL,FILE='MOLINP',STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#ifndef SYS_MANY
      OPEN(LUMOL,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
      DO 200 IMLINE = 1,NMLINE
         WRITE(LUMOL,'(A)') MLINE(IMLINE)
 200  CONTINUE
      CLOSE(LUMOL)
      FILENM = 'XXXX_mol.inp'
      WRITE(FILENM(1:4),'(I4)') (ITRNMR + 1)
      IF (ITRNMR .LT. 9)   WRITE(FILENM(3:3),'(A1)') '0'
      IF (ITRNMR .LT. 99)  WRITE(FILENM(2:2),'(A1)') '0'
      IF (ITRNMR .LT. 999) WRITE(FILENM(1:1),'(A1)') '0'
#ifdef SYS_MANY
      OPEN(LUMOL,FILE=FILENM, STATUS='UNKNOWN',FORM='FORMATTED')
      DO 300 IMLINE = 1,NMLINE
         WRITE(LUMOL,'(A)') MLINE(IMLINE)
 300  CONTINUE
      CLOSE(LUMOL)
#endif
#ifdef SYS_MANY
#undef SYS_MANY
#endif
      RETURN
      END

C  /* Deck pnches */
      SUBROUTINE PNCHES(HESINT)
C     
C     Punch molecular Hessian to the file DALTON.HES, this file can be
C     used to obtain initial Hessian (1st order methods) for restarts
C     or other runs.
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <taymol.h>
      PARAMETER (LUHES = 58)
      DIMENSION HESINT(2*MXCOOR,2*MXCOOR)
      CLOSE(LUHES)
#if defined (SYS_VAX)
      OPEN(LUHES,STATUS='NEW',FORM='FORMATTED',
     &                                 CARRIAGECONTROL='LIST')
#endif
#if defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_SUN)
      OPEN(LUHES,FILE='DALTON.HES',STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#if defined (SYS_IBM)
      OPEN(LUHES,FILE='MOLHES',STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#if !defined (SYS_CRAY) && !defined (SYS_VAX) && !defined (SYS_ALLIANT) && !defined (SYS_CONVEX) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_IBM) && !defined (SYS_IRIX) && !defined (SYS_SUN)
      OPEN(LUHES,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
      IF (REDINT) THEN
         WRITE(LUHES,*) IINTCR
         WRITE(LUHES,*)
         DO 10 J = 1, IINTCR
            DO 12 I = 1, IINTCR
               WRITE(LUHES,*) HESINT(I,J)
 12         CONTINUE
            WRITE(LUHES,*)
 10      CONTINUE
      ELSE
         WRITE(LUHES,*) NCART
         WRITE(LUHES,*)
         DO 15 J = 1, NCART
            DO 17 I = 1, NCART
               WRITE(LUHES,*) HESMOL(I,J)
 17         CONTINUE
            WRITE(LUHES,*)
 15      CONTINUE
      END IF
      CLOSE(LUHES)
      RETURN
      END

C  /* Deck reahes */
      SUBROUTINE REAHES(HESINT,IERR)
C     
C     Read molecular Hessian from the file DALTON.HES, which is then
C     used as initial Hessian in 1st order methods/restarts. IERR is
C     returned with the value 0 if everything is OK. -1 indicates that
C     the file cannot be opened, -2 that the Hessian in the file has wrong
C     wrong dimensions.
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <taymol.h>
      LOGICAL HESEXS
      PARAMETER (LUHES = 58)
      DIMENSION HESINT(2*MXCOOR,2*MXCOOR)
      IERR = 0
      CALL DZERO(HESMOL,MXCOOR*MXCOOR)
      CLOSE(LUHES)
      INQUIRE(FILE='DALTON.HES',EXIST=HESEXS)
      IF (.NOT. HESEXS) THEN
         IERR = -1
         RETURN
      END IF
#if defined (SYS_VAX)
      OPEN(LUHES,STATUS='NEW',FORM='FORMATTED',
     &                                 CARRIAGECONTROL='LIST')
#endif
#if defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_SUN)
      OPEN(LUHES,FILE='DALTON.HES',STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#if defined (SYS_IBM)
      OPEN(LUHES,FILE='DALTON.HES',STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#if !defined (SYS_CRAY) && !defined (SYS_VAX) && !defined (SYS_ALLIANT) && !defined (SYS_CONVEX) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_IBM) && !defined (SYS_IRIX) && !defined (SYS_SUN)
      OPEN(LUHES,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
      READ(LUHES,*) IDIM
      READ(LUHES,*)
      IF (REDINT) THEN
         IF (IDIM .NE. IINTCR) THEN
            IERR = -2
            RETURN
         END IF
         DO 10 J = 1, IINTCR
            DO 12 I = 1, IINTCR
               READ(LUHES,*) HESINT(I,J)
 12         CONTINUE
            READ(LUHES,*)
 10      CONTINUE
      ELSE
         IF (IDIM .NE. NCART) THEN
            IERR = -2
            RETURN
         END IF
         DO 15 J = 1, NCART
            DO 17 I = 1, NCART
               READ(LUHES,*) HESMOL(I,J)
 17         CONTINUE
            READ(LUHES,*)
 15      CONTINUE
      END IF
      CLOSE(LUHES)
      RETURN
      END

C  /* Deck maxelm */
      SUBROUTINE MAXELM(VEC,IDIM,SCLVEC,ISCL,ELMMX)
C     
C     Finds the largest elemement (absolute value) of the vector VEC
C     of dimension IDIM. The value is returned through the
C     variable ELMMX.
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      DIMENSION VEC(IDIM),SCLVEC(MXCOOR)
      LOGICAL SCALE
      ELMMX = ABS(VEC(1))
      IF (ISCL .GE. 1) ELMMX = ELMMX/SCLVEC(1)
      IF (ISCL .GE. 2) ELMMX = ELMMX/SCLVEC(1)
      DO 10 I = 2, IDIM
         ELM = ABS(VEC(I))
         IF (ISCL .GE. 1) ELM = ELM/SCLVEC(I)
         IF (ISCL .GE. 2) ELM = ELM/SCLVEC(I)
         IF (ELM .GT. ELMMX) ELMMX = ELM
 10   CONTINUE
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck numgrd */
      SUBROUTINE NUMGRD(WORK,LWORK,IWOFF,WRKDLM)
C
C     Performs a numerical differentiation in order to get an 
C     finite-difference molecular gradient for use in geometry optimization
C     of MP2, CI or CC wave functions for which there does not exist 
C     analytical implementations
C
C     If I have understood WLKINC correctly, it will do a totally symmetric
C     gradient estimation (thus all we need for a geometry optimization, and
C     no need to run without symmetry).
C
C     K.Ruud, Feb.-6 1997
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
C
      PARAMETER (LUMOL = 9, DP5 = 0.50D0, D100=100.0D0)
      LOGICAL EXHER, EXSIR, EXABA, BIG
      DIMENSION WORK(LWORK), COOR(3,MXCENT)
C
#include <molinp.h>
#include <taymol.h>
#include <symmet.h>
#include <nuclei.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <optinf.h>
#include <infopt.h>
C
      ESTART = ECORR
      NCART  = NCRREP(0,1)
      NCRIND = 3*NUCIND
      KSTPSY = IWOFF + 2
      KSTPCA = KSTPSY + NCART
      KLAST  = KSTPCA + NCRIND
      IPRINT = 0
      IF (.NOT. CHGRDT) THEN
         GRDTHR = 1.0D-4
         THRSTP = GRDTHR
         THRERG = GRDTHR
         THRSYM = SQRT(THRERG)
         WRITE (LUPRI,'(/A)') ' WARNING>>>> Due to limitations '//
     &        'of the accuracy of the numerical gradients'
         WRITE (LUPRI,'(A)') ' WARNING>>>> thresholds for convergence'
     &        //' of geometry optimization has been reset'
         WRITE (LUPRI,'(/A,/,3(/,20X,A,F10.8))') ' New thresholds:',
     &        'Gradient norm  ',GRDTHR,'Step norm      ',THRSTP,
     &        'Energy change  ',THRERG
      END IF
C
C     Start loop over displaced geometries
C
      DO IPOINT = 0, 2*NCART
         IF (IPOINT .NE. 0) THEN
C     
C     Do numerical differentiation
C
            SCLFCK = 1.0D0
            I = (IPOINT+1)/2
            DO 110 IATOM = 1, NUCIND
               DO 111 ICOOR = 1, 3
                  ICCOOR = 3*(IATOM - 1) + ICOOR
                  ISCOOR = IPTCNT(ICCOOR,0,1)
                  IF (ISCOOR .EQ. I) THEN
                     SCLFCK = SQRT(FMULT(ISTBNU(IATOM)))
                  END IF
 111           CONTINUE
 110        CONTINUE
C
            IF (MOD(IPOINT,2) .EQ. 1) THEN
               GRDMOL(I) = (ECORR - ESTART)*SCLFCK/DISPLA
            ELSE
               GRDMOL(I) = DP5*(GRDMOL(I) 
     &                   - SCLFCK*(ECORR - ESTART)/DISPLA)
            END IF
         END IF
         CALL WLKINC(WORK(KSTPSY),WORK(KSTPCA),COOR,DISPLA,NCART,
     &               NCRIND,IPOINT,IPRINT)
C
C     Update geometry in common block
C
         MLINE(NMLAU)(20:20) = ' '
         DO 100 N = 1, NUCIND
            CRX = COOR(1,N)
            CRY = COOR(2,N)
            CRZ = COOR(3,N)
            BIG = (ABS(CRX) .GT. D100 .OR.
     *           ABS(CRY) .GT. D100 .OR.
     *           ABS(CRZ) .GT. D100)
            NC = NCLINE(N)
            IF (NC .NE. 0) THEN
               READ  (MLINE(NC),9100) NAME
               IF (BIG) THEN
                  WRITE (MLINE(NC),9200) NAME,CRX,CRY,CRZ,'*'
               ELSE
                  WRITE (MLINE(NC),9300) NAME,CRX,CRY,CRZ,'*'
               END IF
            END IF
 100     CONTINUE
C
C Punch MOLECULE input with updated coordinates to LUMOL
C
#if    defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) \
    || defined (SYS_PARAGON) || defined (SYS_CRAY)   || defined (SYS_IRIX) \
    || defined (SYS_SUN)     || definbed(SYS_LINUX)  || defined (SYS_HPUX) 
#define SYS_MANY
#endif
#if defined (SYS_VAX)
         OPEN(LUMOL,STATUS='NEW',FORM='FORMATTED',
     *        CARRIAGECONTROL='LIST')
#endif
#ifdef SYS_MANY
         OPEN(LUMOL,FILE='MOLECULE.INP',STATUS='UNKNOWN',
     &        FORM='FORMATTED')
#endif
#if defined (SYS_IBM)
         OPEN(LUMOL,FILE='ABAMOL',STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#ifndef SYS_MANY
         OPEN(LUMOL,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
         REWIND (LUMOL)
         DO 200 IMLINE = 1,NMLINE
            WRITE (LUMOL,'(A)') MLINE(IMLINE)
 200     CONTINUE
         CLOSE (UNIT=LUMOL)
C
         EXHER = .FALSE.
         EXSIR = .FALSE.
         EXABA = .FALSE.
         RDINPC = .FALSE.
         CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
      END DO
C
 9100 FORMAT (A4)
 9200 FORMAT (A4,3F20.10,7X,A1)
 9300 FORMAT (A4,3F20.16,7X,A1)
C
#ifdef SYS_MANY
#undef SYS_MANY
#endif
      RETURN
      END   
#endif
