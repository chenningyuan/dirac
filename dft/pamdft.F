C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dftinp */
      SUBROUTINE DFTINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for DFT-module
C
C     Written by T.Saue and T.Helgaker - feb 2000
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 9)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF, LBIT, RESET, NOAUTOCC,SUCCESS
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER CTEMP*72, LINE*80, TEXT*20
      DIMENSION WORK(LWORK)
C
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dfterg.h"
#include "dftcom.h"
#include "functionals.h"
C
      SAVE SET
      DATA TABLE /'.DFTFUN','.FROMVX','.DFTAC ','.RADINT',
     &            '.ANGINT','.PRINT ','.INTCHK','.NOPRUN',
     &            '.ANGMIN'/
      DATA SET/.FALSE./
C
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
C     Initialize /DFTCOM/ and integration parameters
C     ==============================================
C
      IPRDFT = 2
      DFTRUN = .FALSE.
      DFTADD = .FALSE.
      GRDONE = .FALSE.
      DFTPOT = .FALSE.
      DFTORD = .FALSE.
      DFTASC = .FALSE.
      RADINT = 1.0D-8
      ANGINT = 41
      INTCHK = 0
      NOPRUN = .FALSE.
      ANGMIN = 15  
C
C     Process input for DFTCAL
C     ========================
C
      NEWDEF = (WORD .EQ. '*DFT   ')
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DFTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DFTINP.')
 1             CONTINUE
C&&&: DFT functional selection.
               READ(LUCMD,'(A80)') LINE
               CALL DFTINPUT(LINE, SUCCESS)
               IF(.NOT.SUCCESS) THEN
                  CALL DFTLISTFUNCS
                  CALL QUIT("DFT FUNCTIONAL Read failed.")
               ENDIF
               DFTADD = .TRUE.
               DFTRUN = .TRUE.
               GO TO 100
 2             CONTINUE
C&&& DFTPOT: calculate Kohn-Sham matrix directly using v_xc
                 DFTPOT = .TRUE.
                 GO TO 100
 3               CONTINUE
C&&& DFTAC: asymptotic correction
                 DFTASC = .TRUE.
                 DFTPOT = .TRUE.
                 READ (LUCMD,*) DFTIPT, DFTBR1, DFTBR2
                 GO TO 100
 4               CONTINUE
C&&& RADINT: Give relative precision of the radial integration
                 READ(LUCMD,*) RADINT
                 GO TO 100
 5               CONTINUE
C&&& ANGINT: Give precision of the angular integration.
C            The integration is exact for spherical harmonics up to
C            order ANGINT
            READ(LUCMD,*) ANGINT
            GO TO 100
 6          CONTINUE
C&&& PRINT: Print level for DFT-module
                 READ(LUCMD,*) IPRDFT
                 GO TO 100
 7               CONTINUE
C&&& INTCHK: Test grid by analyzing numerical SMAT
              READ(LUCMD,*) INTCHK
                 GO TO 100
 8               CONTINUE
C&&& NOPRUN: Turn off pruning of angular grid.
              NOPRUN = .TRUE.
                 GO TO 100
 9               CONTINUE
C&&& ANGMIN: Set minimum order of Lebedev scheme after pruning
                 READ(LUCMD,*) ANGMIN
                 GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in DFTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in DFTINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &   'DFTINP: Set-up for Kohn-Sham calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I5)') '* Print level: ',IPRDFT
C
      CALL DFTREPORT()
 90   CONTINUE
      WRITE(LUPRI,'(A)')
     &     '* Radial quadrature according to ',
     &     '     R.Lindh, P.-Aa. Malmqvist and L.Gagliardi'
      WRITE(LUPRI,'(A,E12.5)') 
     &     '  Precision of radial quadrature set to: ', RADINT
      WRITE(LUPRI,'(A)')
     &     '* Angular quadrature according to Lebedev scheme.'
      WRITE(LUPRI,'(A,I2)')
     &     '  Angular integration exact up to order L = ', ANGINT
      IF(NOPRUN) THEN
         WRITE(LUPRI,'(A)')
     &        '* Pruning of angular grid turned off.'
      END IF
      IF (INTCHK.GT.0) THEN
         WRITE(LUPRI,'(A)')
     &        '* Numerical integration will be tested.'
      END IF
C
  999 CONTINUE
      RETURN
C
      END
