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
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (NTABLE = 15)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF, LBIT, RESET, NOAUTOCC
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER CTEMP*72, LINE*80, TEXT*20
      DIMENSION WORK(LWORK)
C
#include <mxcent.h>
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbham.h>
#include <dfterg.h>
#include <dftcom.h>
C
      SAVE SET
      DATA TABLE /'.XALPHA','.HFXFAC','.LDA   ','.BLYP  ','.B3LYP ',
     &            '.FROMVX','.WAH   ','.DFTWEI','.QUADRA','.DFTAC ',
     &            '.RADINT','.PRINT ','.XXXXXX','.XXXXXX','.XXXXXX'/
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
C     Initialize /DFTCOM/
C     ===================
C
      IQDTYP = 1
      IPRDFT = 0
      WDFTX  = D0
      WDFTC  = D0
      WDFTL  = D0
      WDFTB  = D0
      DFTRUN = .FALSE.
      DFTADD = .FALSE.
      GRDONE = .FALSE.
      DFTPOT = .FALSE.
      DFTORD = .FALSE.
      DFTASC = .FALSE.
      RADINT = 1.0D-8
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15),I
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
    1          CONTINUE
C&&&: X alpha calculation
                 READ(LUCMD,*) ALPHA
                 ITYP = 1
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 WDFTX  = 1.5D0*ALPHA
                 WDFTB  = D0
                 WDFTC  = D0
                 WDFTL  = D0
                 HFXFAC = D0
                 GO TO 100
    2          CONTINUE
C&&& HFXFAC: weight of pure exchange
                 READ(LUCMD,*) HFXFAC
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 GO TO 100
    3          CONTINUE
C&&&: LDA local density approximation
                 ITYP = 1
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 WDFTX  = D1
                 WDFTB  = D0
                 WDFTC  = D1
                 WDFTL  = D0
                 HFXFAC = D0
                 GO TO 100
    4          CONTINUE
C&&&: BLYP: B88 exchange + LYP correlation
                 ITYP = 2
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 WDFTX  = D1
                 WDFTB  = D1
                 WDFTC  = D0
                 WDFTL  = D1
                 HFXFAC = D0
                 GO TO 100
    5          CONTINUE
C&&&: B3LYP: hybrid functional
                 ITYP = 3
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 WDFTX  = 0.80D0
                 WDFTB  = 0.72D0
                 WDFTL  = 0.81D0
                 WDFTC  = D1 - WDFTL
                 HFXFAC = D1 - WDFTX
                 GO TO 100
    6          CONTINUE
C&&& DFTPOT: calculate Kohn-Sham matrix directly using v_xc
                 DFTPOT = .TRUE.
                 GO TO 100
    7          CONTINUE
C&&& WAH: Watson-Amis-Handy functional
                 ITYP = 4
                 DFTORD = .TRUE.
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 WDFTX  = 0.95D0
                 WDFTB  = 0.72D0
                 WDFTL  = 0.81D0
                 WDFTC  = D1 - WDFTL
                 HFXFAC = D1 - WDFTX
                 GO TO 100
    8          CONTINUE
C&&& OWN: your own hybrid functional
                 ITYP = 5
                 DFTADD = .TRUE.
                 DFTRUN = .TRUE.
                 WDFTX  = 0.95D0
                 WDFTB  = 0.72D0
                 WDFTL  = 0.81D0
                 READ(LUCMD,*) WDFTX
                 WDFTC  = D1 - WDFTL
                 HFXFAC = D1 - WDFTX
                 GO TO 100
    9          CONTINUE
C&&& QUADRA: number of quadrature points
C            IQUAD = 2    very coarse
C            IQUAD = 4    coarse
C            IQUAD = 5    medium
C            IQUAD = 1    dense
C            IQUAD = 7    very dense
                 READ (LUCMD,*) IQDTYP
                 GO TO 100
   10          CONTINUE
C&&& DFTAC: asymptotic correction
                 DFTASC = .TRUE.
                 DFTPOT = .TRUE.
                 READ (LUCMD,*) DFTIPT, DFTBR1, DFTBR2
                 GO TO 100
 11              CONTINUE
C&&& RADINT: Give relative precision of the radial integration
                 READ(LUCMD,*) RADINT
                 GO TO 100
 12              CONTINUE
C&&& PRINT: Print level for DFT-module
                 READ(LUCMD,*) IPRDFT
                 GO TO 100
 13              CONTINUE
                 GO TO 100
 14              CONTINUE
                 GO TO 100
 15              CONTINUE
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
      WRITE(LUPRI,'(1X,A)') '* Print level: ',IPRDFT
C
      GO TO (91,92,93,94,95),ITYP
 91   CONTINUE
        DFTTYP = 'LDA   '
        WRITE(LUPRI,'(1X,A)') '* Local density approximation'
        WRITE(LUPRI,'(4X,A)') '- Slater exchange'
        WRITE(LUPRI,'(4X,A)') '- VWN correlation functional'
      GOTO 90
 92   CONTINUE
        DFTTYP = 'BLYP  '
        WRITE(LUPRI,'(1X,A)') '* BLYP functional'
      GOTO 90
 93   CONTINUE
        DFTTYP = 'B3LYP '
        WRITE(LUPRI,'(1X,A)') '* B3LYP functional'
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- Slater exchange                :',WDFTX
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- Becke(B88) exchange correction :',WDFTB
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- exact exchange                 :',HFXFAC
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- WVN correlation                :',(D1-WDFTL)
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- LYP correlation correction     :',WDFTL
      GOTO 90
 94   CONTINUE
        DFTTYP = 'WAH   '
        WRITE(LUPRI,'(1X,A)') '* Watson-Amis-Handy functional'
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- Slater exchange                :',WDFTX
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- Becke(B88) exchange correction :',WDFTB
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- exact exchange                 :',HFXFAC
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- WVN correlation                :',(D1-WDFTL)
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- LYP correlation correction     :',WDFTL
      GOTO 90
 95   CONTINUE
        DFTTYP = 'OWN   '
        WRITE(LUPRI,'(1X,A)') '* User defined functional'
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- Slater exchange               :',WDFTX
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- Becke(B88) exchange correction:',WDFTB
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- exact exchange                :',HFXFAC
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- WVN correlation               :',(D1-WDFTL)
        WRITE(LUPRI,'(4X,A,F6.2)') 
     &     '- LYP correlation correction    :',WDFTL
      GOTO 90
 90   CONTINUE
      WRITE(LUPRI,'(A)')
     &  '* Radial quadrature according to ',
     &  '     R.Lindh, P.-Aa. Malmqvist and L.Gagliardi'
      WRITE(LUPRI,'(1X,A,E12.5)') 
     &  '  Precision of radial quadrature set to: ', RADINT
      WRITE(LUPRI,'(A)')
     &  '* Angular quadrature accroding to Lebedev scheme.'
C
  999 CONTINUE
      RETURN
C
      END
