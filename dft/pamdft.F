C
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
      PARAMETER (NTABLE = 20)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF, LBIT, RESET, NOAUTOCC
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, AC_WORD*4
      DIMENSION WORK(LWORK)
C
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbham.h"
#include "dfterg.h"
#include "dftcom.h"
#include "dftac.h"
#include "cdft.h"
#include "dftbuf.h"
#include "densfit.h"
#include "functionals.h"
C
      SAVE SET
      DATA TABLE /'.EXPORT','.FROMVX','.IMPORT','.RADINT',
     &            '.ANGINT','.PRINT ','.INTCHK','.NOPRUN',
     &            '.ANGMIN','.ATSIZE','.ERRELS','.GRAC  ',
     &            '.SAOP  ','.SAOP! ','.SDFT  ','.VRG   ',
     &            '.COLLIN','.GAUNTS','.ZIPGRI','.XXXXXX'/
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
C
C     Initialize /DFTCOM/ and integration parameters
C     ==============================================
C
      IPRDFT = 2
      GRDONE = .FALSE.
C     ... GRDONE false: calculate grid
C         (this also ensures new grid is calculated for each
C          geometry in a geometry optimization) /hjaaj
      DFTPOT = .FALSE.
      DFTORD = .FALSE.
      DFTASC = .FALSE.
      RADINT = 1.0D-13
      ANGINT = 41
      INTCHK = 0
      DFTHR0 = 1.0D-9
      DFTHRL = 1.0D-10
      DFTELS = 1.0D-3
      NOPRUN = .FALSE.
      ANGMIN = 15
      ATSIZE = .FALSE.
      import_grid = .false.
      export_grid = .false.
      ZIPGRID = .FALSE.
C
C
C     initialize /DFTAC/
C     ==================
C
      LBIT_DFTAC           = 0
      ICOUNTER_DFTAC_DEBUG = 0
C
C

!     initialize cdft.h
      LBIT_CDFT = 0

C
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
                     GOTO ( 1, 2, 3, 4,
     &                      5, 6, 7, 8,
     &                      9,10,11,12,
     &                     13,14,15,16,
     &                     17,18,19,20), I
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

!              export grid
 1             continue
               export_grid = .true.
               read(lucmd,'(A30)') gridfile
               go to 100

 2             CONTINUE
C&&& DFTPOT: calculate Kohn-Sham matrix directly using v_xc
                 DFTPOT = .TRUE.
                 GO TO 100

!              import grid
 3             continue
               import_grid = .true.
               read(lucmd,'(A30)') gridfile
               go to 100

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
 10              CONTINUE
C&&& ATSIZE: Estimate relative atomic sizes for use in the Becke
C            partitioning scheme from atomic contributions
C            to the density.
                 ATSIZE = .TRUE.
                 GO TO 100
 11              CONTINUE
C&&& ERRELS: Threshold for accepted error in integrated number of electrons
                 READ(LUCMD,*) DFTELS
                 GO TO 100
 12              CONTINUE
C&&& GRAC: gradient regulated asymptotic correction
                 LBIT_DFTAC = LBIT_DFTAC + 1
                 READ(LUCMD,'(A4)') AC_WORD
                 IF(AC_WORD .EQ. 'LB94') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 32
                 ELSEIF(AC_WORD .EQ. 'LBAL') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 64
                 ELSEIF(AC_WORD .EQ. 'LBAl') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 64
                 ELSEIF(AC_WORD .EQ. 'LBal') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 64
                 ELSE
                   CALL QUIT('keyword following .GRAC not recognized')
                 ENDIF
                 READ(LUCMD,*) GRAC_ALPHA, GRAC_BETA,
     &                         DFTAC_IP, DFTAC_THRESHOLD
                 GOTO 100
 13              CONTINUE
C&&& SAOP: statistical averaging of (model) orbital potentials
                 LBIT_DFTAC = LBIT_DFTAC + 4
                 READ(LUCMD,'(A4)') AC_WORD
                 IF(AC_WORD .EQ. 'LB94') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 32
                 ELSEIF(AC_WORD .EQ. 'LBAL') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 64
                 ELSEIF(AC_WORD .EQ. 'LBAl') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 64
                 ELSEIF(AC_WORD .EQ. 'LBal') THEN
                   LBIT_DFTAC = LBIT_DFTAC + 64
                 ELSE
                   CALL QUIT('keyword following .SAOP not recognized')
                 ENDIF
                 GOTO 100
 14              CONTINUE
C&&& SAOP!: "original" SAOP as in JCP 112, 1344 (2000).
C                             SAOP
                 LBIT_DFTAC = 4
C                                 SAOP with response part
     &                          + 16
C                                      asymptotic part is LBalpha
     &                               + 64
                 DFTHAM     = .TRUE.
                 ALDAON     = .TRUE.
#ifdef MOD_DFT
                 CALL DOALDA(.TRUE.)
#endif
                 GOTO 100

 15              CONTINUE
                 CALL QUIT('missing CDFT routines')
                 GOTO 100

 16              CONTINUE
                 CALL QUIT('missing CDFT routines')
                 GOTO 100

 17              CONTINUE
                 CALL QUIT('missing CDFT routines')
                 GOTO 100
C     .GAUNTSCALE
 18              CONTINUE
                 LSCALE_DFT_GAUNT=.TRUE.
                 GOTO 100

 19              CONTINUE
                 ZIPGRID = .TRUE.
                 READ(LUCMD,'(A30)') ZIPGRID_TEXT
                 GOTO 100

 20              CONTINUE
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
      IF(.NOT.DFTHAM) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   'DFTINP: Set-up for Kohn-Sham calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I5)') ' * Print level: ',IPRDFT
C
 90   CONTINUE
      WRITE(LUPRI,'(A)')
     &     ' * Radial quadrature according to',
     &     '          R. Lindh, P.-Aa. Malmqvist and L. Gagliardi'
      WRITE(LUPRI,'(A,E12.5)')
     &     '   Precision of radial quadrature set to: ', RADINT
      WRITE(LUPRI,'(A,E12.5)')
     &    '   Accepted error in integrated number of electrons: ',DFTELS
      WRITE(LUPRI,'(A)')
     &     ' * Angular quadrature according to Lebedev scheme.'
      WRITE(LUPRI,'(A,I2)')
     &     '   Angular integration exact up to order L = ', ANGINT
      IF(NOPRUN) THEN
         WRITE(LUPRI,'(A)')
     &        ' * Pruning of angular grid turned off.'
      END IF
      IF (INTCHK.GT.0) THEN
         WRITE(LUPRI,'(A)')
     &        ' * Numerical integration will be tested.'
      END IF
C
      IF(DFTPOT) THEN
        WRITE(LUPRI,'(A)')
     &        ' * Explicit exchange-correlation potential (.FROMVX):'
        WRITE(LUPRI,'(A)')
     &        '   WARNING: this option is not available in this version'
        WRITE(LUPRI,'(A)')
     &        '   and will be ignored (DFTPOT = .FALSE.)'
        DFTPOT = .FALSE.
      ENDIF
C
C             GRAC
      IF(LBIT(LBIT_DFTAC,1)) THEN
        WRITE(LUPRI,'(A)')
     &        ' * Gradient regulated asymptotic correction:'
        WRITE(LUPRI,'(A,F12.6)') '   alpha = ', GRAC_ALPHA
        WRITE(LUPRI,'(A,F12.6)') '   beta  = ', GRAC_BETA
        WRITE(LUPRI,'(A)')
     &        '   threshold for difference in HOMO eigenvalue'
        WRITE(LUPRI,'(A,E12.6)')
     &        '   below which asymptotic correction is switched on: ',
     &        DFTAC_THRESHOLD
        WRITE(LUPRI,'(A,E12.6)')
     &        '   bulk potential will be shifted by (IP): ', DFTAC_IP
      ENDIF
C             SAOP
      IF(LBIT(LBIT_DFTAC,3)) THEN
        WRITE(LUPRI,'(A)')
     &        ' * Statistical averaging of (model) orbital potentials:'
      ENDIF
C             SAOP with response part
      IF(LBIT(LBIT_DFTAC,5)) THEN
        WRITE(LUPRI,'(A)')
     &        '   original SAOP as in JCP 112, 1344 (2000) chosen'
        WRITE(LUPRI,'(A)')
     &        '   this activates ALDA kernel'
        WRITE(LUPRI,'(A)')
     &        '   the functional under .DFT is expected to be GLLBhole'
      ENDIF
C
      IF(LBIT(LBIT_DFTAC,6)) THEN
        WRITE(LUPRI,'(A)') '   the asymptotic potential is LB94'
      ENDIF
      IF(LBIT(LBIT_DFTAC,7)) THEN
        WRITE(LUPRI,'(A)') '   the asymptotic potential is LBalpha'
      ENDIF
C
  999 CONTINUE
C
C
      RETURN
      END
