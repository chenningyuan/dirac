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
      PARAMETER (NTABLE = 28)
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
     &            '.COLLIN','.GAUNTS','.ZIPGRI','.ALDA  ',
     &            '.ALDA+ ','.ALDA- ','.XALDA ','.XALDA+',
     &            '.XALDA-','.XXXXXX','.XXXXXX','.XXXXXX'/
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

      alda_hs = .false.
      alda_ha = .false.
      xalda   = .false.


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
     &                     17,18,19,20,
     &                     21,22,23,24,
     &                     25,26,27,28), I
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
                 alda_hs = .true.
                 alda_ha = .true.
                 call setaldahs(.true.) !activate ALDA for the H+ part
                 call setaldaha(.true.) !activate ALDA for the H- part
                 GOTO 100

 15              CONTINUE
#ifdef MOD_CDFT
                 CALL SET_LBIT(LBIT_CDFT,1)
#else
                 CALL QUIT('missing CDFT routines')
#endif /* MOD_CDFT */
                 GOTO 100

 16              CONTINUE
#ifdef MOD_CDFT
                 CALL SET_LBIT(LBIT_CDFT,2)
#else
                 CALL QUIT('missing CDFT routines')
#endif /* MOD_CDFT */
                 GOTO 100

 17              CONTINUE
#ifdef MOD_CDFT
                 CALL SET_LBIT(LBIT_CDFT,3)
#else
                 CALL QUIT('missing CDFT routines')
#endif /* MOD_CDFT */
                 GOTO 100
C     .GAUNTSCALE
 18              CONTINUE
                 LSCALE_DFT_GAUNT=.TRUE.
                 GOTO 100

 19              CONTINUE
                 ZIPGRID = .TRUE.
                 READ(LUCMD,'(A30)') ZIPGRID_TEXT
                 GOTO 100

!                .ALDA
 20              continue
                 alda_hs = .true.
                 alda_ha = .true.
                 call setaldahs(.true.) !activate ALDA for the H+ part
                 call setaldaha(.true.) !activate ALDA for the H- part
                 go to 100

!                .ALDA+
!                use ALDA only for the      hermitian part
 21              continue
                 alda_hs = .true.
                 call setaldahs(.true.) !activate ALDA for the H+ part
                 go to 100

!                .ALDA-
!                use ALDA only for the anti-hermitian part
 22              continue
                 alda_ha = .true.
                 call setaldaha(.true.) !activate ALDA for the H- part
                 go to 100

!                .XALDA
!                use        (1-f)*S + VWN + f*HFx
!                instead of       S + VWN
 23              continue
                 xalda   = .true.
                 alda_hs = .true.
                 alda_ha = .true.
                 call setxalda(.true.)
                 call setaldahs(.true.) !activate ALDA for the H+ part
                 call setaldaha(.true.) !activate ALDA for the H- part
                 go to 100

!                .XALDA+
!                use XALDA only for the      hermitian part
 24              continue
                 xalda   = .true.
                 alda_hs = .true.
                 call setxalda(.true.)
                 call setaldahs(.true.) !activate ALDA for the H+ part
                 go to 100

!                .XALDA-
!                use XALDA only for the anti-hermitian part
 25              continue
                 xalda   = .true.
                 alda_ha = .true.
                 call setxalda(.true.)
                 call setaldaha(.true.) !activate ALDA for the H- part
                 go to 100

 26              continue
                 go to 100

 27              continue
                 go to 100

 28              continue
                 go to 100

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


!     consistency
!     ===========

      if (dftpot) then
        call quit('explicit xc potential (.FROMVX) not available')
      end if

      if (lbit(lbit_cdft, 3) .and. .not. lbit(lbit_cdft, 1)) then
        call quit('please use .COLLINEAR together with .SDFT')
      end if


!     print section
!     =============

      if (.not. dftham) go to 999

      call header('Set-up for Kohn-Sham calculation:', 2)


!     print level
!     -----------

      write(lupri, '(a, i5)') ' * Print level: ', iprdft


!     integration grid
!     ----------------

      write(lupri, '(a)')
     &     ' * Numerical integration grid:'

      if (import_grid) then
      write(lupri, '(a)')
     &     '   - integration grid imported from file'
      else

      write(lupri, '(a)')
     &     '   - radial quadrature according to '//
     &     'R. Lindh, P.-Aa. Malmqvist, and L. Gagliardi'
      write(lupri, '(a, e12.5)')
     &     '     precision of radial quadrature set to: ', radint
      write(lupri, '(a, i2)')
     &     '   - angular quadrature using the Lebedev scheme, '//
     &     'exact up to order L = ', angint

      if (noprun) then
      write(lupri, '(a)')
     &     '   - pruning of angular grid turned off'
      end if

      if (intchk .gt. 0) then
      write(lupri, '(a)')
     &     '   - numerical integration will be tested'
      end if

      write(lupri, '(a, e12.5)')
     &     '   - accepted error in integrated number '//
     &     'of electrons: ', dftels

      if (atsize) then
      write(lupri, '(a)')
     &     '   - estimate relative atomic sizes for use in the Becke'
      write(lupri, '(a)')
     &     '     partitioning scheme from atomic contributions '//
     &     'to the density'
      end if

      end if !if (import_grid)

      if (export_grid) then
      write(lupri, '(a)')
     &     '   - integration grid exported to file'
      end if


!     asymptotic correction
!     ---------------------

!     GRAC
      if (lbit(lbit_dftac, 1)) then
      write(lupri, '(a)')
     &     ' * Gradient regulated asymptotic correction:'
      write(lupri, '(a, e12.5)')
     &     '   - alpha = ', grac_alpha
      write(lupri, '(a, e12.5)')
     &     '   - beta  = ', grac_beta
      write(lupri, '(a)')
     &     '   - threshold for difference in HOMO eigenvalue'
      write(lupri, '(a, e12.5)')
     &     '     below which asymptotic correction is switched on: ', 
     &     dftac_threshold
      write(lupri, '(a, e12.5)')
     &     '   - bulk potential will be shifted by (IP): ', dftac_ip
      end if

!     SAOP
      if (lbit(lbit_dftac, 3)) then
      write(lupri, '(a)')
     &     ' * Statistical averaging of (model) orbital potentials:'
      end if

!     SAOP with "response part"
      if (lbit(lbit_dftac, 5)) then
      write(lupri, '(a)')
     &     '   - original SAOP as in JCP 112, 1344 (2000) chosen'
      write(lupri, '(a)')
     &     '   - this activates ALDA xc kernel'
      write(lupri, '(a)')
     &     '   - the functional under .DFT is expected to be GLLBhole'
      end if

!     asymptotic potential
      if (lbit(lbit_dftac, 6)) then
      write(lupri, '(a)')
     &     '   - the asymptotic potential is LB94'
      end if
      if (lbit(lbit_dftac, 7)) then
      write(lupri, '(a)')
     &     '   - the asymptotic potential is LBalpha'
      end if


!     SDFT
!     ----

      if (lbit(lbit_cdft, 1)) then

      write(lupri, '(a)')
     &     ' * Spin density contribution in response calculations:'
      if (lbit(lbit_cdft, 3)) then
      write(lupri, '(a)')
     &     '   - use the collinear approximation as '//
     &     'a definition of the spin density'
      else
      write(lupri, '(a)')
     &     '   - use the norm of the spin magnetization vector as '//
     &     'a definition of the spin density (noncollinear definition)'
      end if

      else
      write(lupri, '(a)')
     &     ' * No spin density contribution in response calculations.'
      write(lupri, '(a)')
     &     '     (if you want the spin density contribution, use .SDFT)'
      end if


!     ALDA
!     ----

      if (alda_hs .or. alda_ha) then

      write(lupri, '(a)')
     &     ' * Adiabatic local density approximation (ALDA):'
      write(lupri, '(a)')
     &     '   - approximate all functional derivatives '//
     &     'beyond the xc potential by SVWN derivatives'

      if (alda_hs .and. .not. alda_ha) then
      write(lupri, '(a)')
     &     '   - use ALDA only for the Hermitian part '//
     &     '(density contribution)'
      write(lupri, '(a)')
     &     '   - use the proper xc kernel for the anti-Hermitian '//
     &     'part (spin density contribution)'
      end if

      if (alda_ha .and. .not. alda_hs) then
      write(lupri, '(a)')
     &     '   - use ALDA only for the anti-Hermitian part '//
     &     '(spin density contribution)'
      write(lupri, '(a)')
     &     '   - use the proper xc kernel for the Hermitian '//
     &     'part (density contribution)'
      end if

      if (xalda) then
      write(lupri, '(a)')
     &     '   - keep the fraction of exact exchange of the '//
     &     'xc functional in the solution of the response equation'
      else
      write(lupri, '(a)')
     &     '   - for hybrid functionals exact exchange is '//
     &     'switched off in the solution of the response equation'
      write(lupri, '(a)')
     &     '     (if you want to keep it, use .XALDA)'
      end if

      end if !if (alda_hs .or. alda_ha)


!     thresholds
!     ----------

      write(lupri, '(a)')
     &     ' * Thresholds:'
      write(lupri, '(a, e12.5)')
     &     '   - DFTHR0 = ', dfthr0
      write(lupri, '(a, e12.5)')
     &     '   - DFTHRL = ', dfthrl

  999 continue
      return
      end
