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
      PARAMETER (NTABLE = 12)
      PARAMETER (RTOL = 1.0D-15,D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF, LBIT, RESET, NOAUTOCC
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
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
      DATA TABLE /'.XXXXXX','.FROMVX','.DFTAC ','.RADINT',
     &            '.ANGINT','.PRINT ','.INTCHK','.NOPRUN',
     &            '.ANGMIN','.ATSIZE','.ERRELS','.XXXXXX'/
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
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12),I
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
               CALL QUIT('Functional selection moved to **HAMILT')
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
      IF(.NOT.DODFT) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &   'DFTINP: Set-up for Kohn-Sham calculation:'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I5)') '* Print level: ',IPRDFT
C
 90   CONTINUE
      WRITE(LUPRI,'(A)')
     &     '* Radial quadrature according to ',
     &     '     R.Lindh, P.-Aa. Malmqvist and L.Gagliardi'
      WRITE(LUPRI,'(A,E12.5)') 
     &     '  Precision of radial quadrature set to: ', RADINT
      WRITE(LUPRI,'(A,E12.5)') 
     &     '  Accepted error in integrated number of electrons: ',DFTELS
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FORTWRT(str, slength)
#include "priunit.h"
      INTEGER unit, slength
      CHARACTER str(slength)
      WRITE(LUPRI, *) str
      END
