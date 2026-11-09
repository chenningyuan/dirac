!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

C FILE    : pamprp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamprp */
      SUBROUTINE PAMPRP()
!      use dcbpsi_mod, only: DOEXACC
!      use dcbprp_mod, only: DOESR, NPRP_WF, PRP_WF
      use eom_driver , only: number_of_states
C***********************************************************************
C
C     Driver for PAMPRP_1.
C
C     Loop over wave functions and call PAMPRP_1 for each of them.
C
C
C***********************************************************************
#ifdef MOD_EXACORR
      use exacorr_cfg, only: exa_input
#endif
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
#include "implicit.h"
      INTEGER :: I
C
C Used in common blocks:
C dcbprp.h: NPRP_WF, PRP_WF(*)
C
#include "dcbpsi.h"
#include "dcbprp.h"
C
      CHARACTER WF*4, WFlabel*6
      INTEGER state
      INTEGER J
C
      CALL QENTER('PAMPRP')
C
      CALL TITLER('Property module','*',125)
      WRITE(LUPRI,9000)
 9000 FORMAT(//,
     &     3X,'This is output from the Dirac property module:',/,
     &    /3X,'* HF & DFT first order properties',/,
     &     6X,'Trond Saue',//,
     &     3X,'* First-order ESR properties',/,
     &     6X,'Hans Joergen Aa. Jensen et al.',//,
     &     3X,'* MP2 first order properties:',/,
     &     6X,'J. N. P. van Stralen, L. Visscher, C. V. Larsen'
     &       ,' and H. J. Aa Jensen, Chem. Phys. 311 (2005) 81.',//
     &     3X,'* KR-RPA second-order properties',/,
     &     6X,'Hans Joergen Aa. Jensen and Trond Saue',//,
     &     3X,'* KR-QR third order properties',/,
     &     6X,'Patrick Norman and Hans Joergen Aa. Jensen',//,
     &     3X,'* Oscillator strengths beyond the electric '
     &        'dipole approximation',/,
     &     6X,'N.H. List, T.R.L Melin, M. van Horn, T. Saue, '
     &        'J. Chem. Phys. 152 (2020) 184110'//
     &     3X,'* Molecular gradient',/,
     &     6X,'Joern Thyssen',//,
     &     3X,'* Additional contributions from:',/,
     &     6X,'Thomas Enevoldsen, Miroslav Ilias (London orbitals)'//
     &     3X,'* First-order exp. val. contribution from:',/,
     &     6X,'Chima Chibueze (KU DHF and KU CCSD)'//)
C
      DO I = 1, NPRP_WF
         WF = PRP_WF(I)
C
         if (.not.doexacc) then
           IF (WF .EQ. 'DHF' .OR. WF .EQ. 'HF' .OR. WF .EQ. 'DFT') THEN
              CALL TDHFORB('CASCI')
              ! for doing a "RESOLVE" with GASCIP code,
              ! CI information is only used if requested.
              CALL SETDC2(0)
           ELSE IF (WF .EQ. 'KRMC') THEN
              CALL TKRMCORB
              CALL SETDC2(0)
           ELSE IF (WF .EQ. 'MP2 ') THEN
              IF (DOESR) cycle
              CALL SETDC2(0)
           ELSE IF (WF .EQ. 'CCSD') THEN
              CALL SETDC2(0)
           ELSE IF (WF .EQ. 'EOM ') THEN
              CALL SETDC2(0)
           ELSE
              CALL QUIT('*** ERROR in PAMPRP *** ' //
     &           'Unknown wave function: ' // WF)
           END IF
         endif !(if .not.doexacc)

         WFlabel = WF
#ifdef MOD_EXACORR
         if ((WF.EQ.'CCSD').AND.(exa_input%lambda)) then
           select case (exa_input%iterative_level)
             case (2)
             ! Distinguish between CCSD and CCD.
             ! CC2 not implemented for first order properties.
             if (exa_input%iterative_scheme == 20) WFlabel = 'CCD'
             case (3)
             WFlabel = 'CCSDT'
             case (4)
             WFlabel = 'CCSDTQ'
             case default
             WFlabel = WF
           end select
         end  if
#endif
         CALL TITLER('Properties for ' // trim(WFlabel) //
     &               ' wave function','*',110)

         IF (WF .EQ. 'EOM ') THEN
           CALL number_of_states(state)
             write(*,*)'# states', state
            do j = 1, state

             CALL PAMPRP_1(WF, j)

            enddo

         ELSE
           CALL PAMPRP_1(WF,1)
         ENDIF
      END DO
C
      CALL QEXIT('PAMPRP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampr1 */
      SUBROUTINE PAMPRP_1(WF,STATE)
C***********************************************************************
C
C     NOTE
C
C     Copyright (c) 1995
C     Authors: Hans Joergen Aa. Jensen and Trond Saue
C
C     This is an experimental code for the evaluation of relativistic
C     SCF molecular linear response properties.  The authors accept
C     no responsibility for the performance of the code or for the
C     correctness of the results.
C
C     If results obtained with this code are published, an appropriate
C     citation would be:  H. J. Aa. Jensen and T. Saue,
C     "XXX version 1.0 1995,
C      a Dirac-Fock molecular linear response properties program".
C
C     Written by T.Saue June 24 1996
C     Revision: 1997/09/15 - jth (addition of gradient code)
C
C***********************************************************************

      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: DOESR, DOEXCPRP, DOEXP, DOSTEX, DOTPA, DOXLR, DOXPP, DOXQR, HDF5IO, MGRAD
!      use dcbexp_mod, only: NEXPP
!      use dcbesr_mod, only: IPRESR, MULTIPESR, NESRP
!      use dcbpsi_mod, only: DOCCM, DOEXACC
      use memory_allocator
      use dirac_cfg
      use fde_mod

#include "implicit.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbesr.h"
#include "dcbpsi.h"
C
      real(8), allocatable :: WORK(:)
      CHARACTER WF*4
      CHARACTER WFEXPORT*4
C     INTEGER, OPTIONAL :: STATE
      INTEGER :: STATE
C
      CALL QENTER('PAMPRP_1')
C
      call legacy_lwork_get(LWORK)
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
      call alloc(WORK,LWORK,id='WORK in PAMPRP_1')
C
C     Make sure common block dcbidx is set
C
      CALL SETDCBIDX()
C
C***********************************************************************
C*****  M O L E C U L A R    G R A D I E N T ***************************
C***********************************************************************
C
      IF (MGRAD) THEN
         IF (WF .NE. 'DHF' .AND. WF .NE. 'HF') THEN
            CALL QUIT('*** ERROR in PAMPRP_1 *** ' //
     &           'Molecular gradient not implemented for ' // WF )
         END IF
         CALL RMOLGRD(WORK(KFREE),LFREE)
      END IF
C
C***********************************************************************
C*****  E X P E C T A T I O N  V A L U E S *****************************
C***********************************************************************
C
      IF (DOEXP) THEN
         CALL MEMGET('REAL',KEXP,NEXPP,WORK,KFREE,LFREE)
C     allow choice to use hdf5 file(s) to be used for expectation values
C     that is supported by the exacorr module. not yet the default as there 
C     are some issues to resolve for construction density matrices,
C     notably in property calculations using mp2 frozen natural orbitals
         if (doexacc.and.hdf5io) then
            call kuprpexp(wf, work(kexp), work(kfree), lfree, state)
         else
             CALL PRPEXP(WF,WORK(KEXP),WORK(KFREE),LFREE,STATE)
         endif
      ELSE
         KEXP = KFREE
      END IF
C
C***********************************************************************
C*****  E S R   P R O P E R T I E S      *******************************
C***********************************************************************
C
      IF (DOESR) THEN
         KFRSAV = KFREE
         CALL ESR_INIT(WF,IPRESR,WORK,KFREE,LFREE)
C
         LENESR = 2*MULTIPESR*MULTIPESR*MAX(4,NESRP)
C        max(4,nesrp): make sure sufficient space for ESRPHASE
         CALL MEMGET2('REAL','ESR',KESR,LENESR,WORK,KFREE,LFREE)
         CALL PRPESR(WF,WORK(KESR),WORK,KFREE,LFREE)
         CALL MEMREL('after PRPESR',WORK,KWORK,KFRSAV,KFREE,LFREE)
C        ... remove this MEMREL if WORK(KESR) to be used in PRPOUT below
      END IF
C
C***********************************************************************
C*****  P P   E X C I T A T I O N   E N E R G I E S  *******************
C***********************************************************************
C
      IF (DOXPP) THEN
         IF (WF .EQ. 'DHF' .OR. WF .EQ. 'HF' .OR. WF .EQ. 'DFT') THEN
            CALL PRPXPP(WORK,KPPF,KOMEGA,KPPCNV,KPPO,KSSYM,KFREE,LFREE)
         ELSE
            CALL QUIT('*** ERROR in PAMPRP_1 *** ' //
     &           'PP excitation energies not implemented for ' // WF )
         END IF
      END IF
C
C***********************************************************************
C*****  L I N E A R   R E S P O N S E  *********************************
C***********************************************************************
C
      IF(DOXLR) THEN
         IF (WF .EQ. 'DHF' .OR. WF .EQ. 'HF' .OR. WF .EQ. 'DFT') THEN
            CALL PRPXLR(WORK,KLRF,KLRCNV,KLRO,KFREE,LFREE)
         ELSE
            CALL QUIT('*** ERROR in PAMPRP_1 *** ' //
     &           'Linear response not implemented for ' // WF )
         END IF
      ELSE
        KLRO   = KFREE
        KLRF   = KFREE
        KLRCNV = KFREE
      ENDIF
C
C***********************************************************************
C*****  Q U A D R A T I C   R E S P O N S E  ***************************
C***********************************************************************
C
      IF (DOXQR .OR. DOTPA .OR. DOEXCPRP) THEN
         IF (WF .NE. 'DHF' .AND. WF .NE. 'HF' .AND. WF .NE. 'DFT')
     &      CALL QUIT('*** ERROR in PAMPRP_1 *** ' //
     &         'Quadratic response not implemented for ' // WF )
         IF (DOXQR) CALL PRPXQR(WORK(KFREE),LFREE)
         IF (DOTPA) CALL PRPTPA(WORK(KFREE),LFREE)
         IF (DOEXCPRP) CALL EXCPRP(WORK(KFREE),LFREE)
      ENDIF
C
C***********************************************************************
C*****  S T E X  *******************************************************
C***********************************************************************
C
      IF (DOSTEX) THEN
         CALL STEX(WORK(KFREE),LFREE)
      ENDIF

C
C***********************************************************************
C*****  standardized P R O P E R T Y   O U T P U T  ********************
C***********************************************************************
C
C hjaaj sep2000:
C     in this version only for DOEXP and for DOXLR;
C     maybe a good idea to assemble final output from various
C     sections in a file which is then copied to output here ?
C
      IF (DOEXP .OR. DOXLR .OR. DOCCM) THEN
         CALL PRPOUT(WORK(KEXP),WORK(KLRF),WORK(KLRCNV),
     &            WORK,KFREE,LFREE)
      END IF
C
C**************************************************************************************
C*****  OUTPUT of density, derivatives of density, and coulomb potential on a grid ****
C**************************************************************************************
C
      if (dirac_cfg_fde_export) then
         CALL TITLER('Output FDE data to file','*',118)
         call fde_export_to_file(wf)
      END IF


C
      call dealloc(WORK)
      CALL QEXIT('PAMPRP_1')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_dipol */
      SUBROUTINE DEF_DIPOLE(INDDIP,IPRINT)
C***********************************************************************
C
C     Define dipole moment operator
C     Written by T.Saue - May 1996
C     Last revision: May 27 1996 - tsaue
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcbxpr_mod, only: PRPNAM
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: DM1 = -1.0D0
#include "dcbxpr.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION INDDIP(3),PFAC(3)
C
C
      PNAME     = 'Dipole length: X'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDDIP(1),ISYXPR,ITRXPR,IPRINT)
C
      PNAME      = 'Dipole length: Y'
      IPTYP      = 1
      NPCOMP     = 1
      PFAC(1)    = DM1
      PLABEL(1)  = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDDIP(2),ISYXPR,ITRXPR,IPRINT)
C
      PNAME      = 'Dipole length: Z'
      IPTYP      = 1
      NPCOMP     = 1
      PFAC(1)    = DM1
      PLABEL(1)  = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDDIP(3),ISYXPR,ITRXPR,IPRINT)
C
      IF(IPRINT.GE.1) THEN
        DO I = 1,3
          INDXPR = INDDIP(I)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Dipole component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_quadru */
      SUBROUTINE DEF_QUADRU(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define quadrupole moment operator
C     Written by T.Saue - May 1996
C     Last revision: May 27 1996 - tsaue
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT
!      use dcbxpr_mod, only: PRPNAM
!      use dcbexp_mod, only: IPQUADRU, LEXPP, MAXEXP, NEXPP
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: D1 = 1.0D0
#include "dcbxpr.h"
#include "dcbexp.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Quadrupole moment'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      PNAME     = 'Quadrupol m.: XX'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XXTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XXQUAD',IPQUADRU(1),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: XY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XYTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XYQUAD',IPQUADRU(2),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: XZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XZTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XZQUAD',IPQUADRU(3),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: YY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YYTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YYQUAD',IPQUADRU(4),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: YZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YZTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YZQUAD',IPQUADRU(5),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: ZZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'ZZTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('ZZQUAD',IPQUADRU(6),LEXPP,NEXPP,INDXPR,MAXEXP)
C
C     Print section
C     =============
C
      IF(IPRINT.GE.1) THEN
        DO I = 1,6
          INDXPR = LEXPP(IPQUADRU(I))
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_epolar */
      SUBROUTINE DEF_EPOLAR(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operators for electric dipole polarizabilities
C
C     Written by T.Saue - June 28 1996
C     Last revision: June 28 - tsaue
C                    July 2006 - M.Ilias, Prievidza
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: IPEPOLAR, IPRPRP
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0
#include "dcbprp.h"
#include "dcbxlr.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
C
      CALL QENTER('DEF_EPOLAR')

      PNAME     = 'X-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(1,1) = INDOP1
      CALL OP1IND('XEPOLB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IPEPOLAR(1,2) = INDOP1
C
      PNAME     = 'Y-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(2,1) = INDOP1
      CALL OP1IND('YEPOLB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IPEPOLAR(2,2) = INDOP1
C
      PNAME     = 'Z-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('ZEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(3,1) = INDOP1
      CALL OP1IND('ZEPOLB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IPEPOLAR(3,2) = INDOP1

C  Print section

      IF (IPRPRP.GE.1) THEN

       CALL PRSYMB(LUPRI,'=',75,0)
       WRITE(LUPRI,'(1X,A)')
     &  'DEF_EPOLAR: Electric dipole polarizability'

       CALL PRSYMB(LUPRI,'-',75,0)
       WRITE(LUPRI,'(1X,A)')
     &'** A operators for the electric dipole polarizability **'
       CALL PRSYMB(LUPRI,'-',75,0)
       DO I = 1, 3
        INDOP1 = IPEPOLAR(I,1)
        INDXPR = LLRAPU(INDOP1)
        CALL WRIXPR(INDOP1,INDXPR)
       ENDDO

       CALL PRSYMB(LUPRI,'-',75,0)
       WRITE(LUPRI,'(1X,A)')
     &'** B operators for the electric dipole polarizability **'
       CALL PRSYMB(LUPRI,'-',75,0)
       DO I = 1, 3
        INDOP1 = IPEPOLAR(I,2)
        INDXPR = LLRBPU(INDOP1)
        CALL WRIXPR(INDOP1,INDXPR)
       ENDDO

      ENDIF

      CALL QEXIT('DEF_EPOLAR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_EFN */
      SUBROUTINE DEF_EFN(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operator for electric field
C
C     Written by T.Saue Feb 26 2003
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
      CALL DEF_EFN1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_efn1 */
      SUBROUTINE DEF_EFN1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C*****************************************************************************
C
C     Define electric field at individual nuclei
C     Based on EFNTYP
C     Written by T. Saue Feb 26 2003
C
C*****************************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: IPTCNT, MAXREP
!      use pgroup_mod, only: REP
      use chrnos_mod, only: CHRNOS
      use chrxyz_mod, only: CHRXYZ
      use nuclei_mod, only: NAMN, NUCIND
!      use dcbxpr_mod, only: PRPNAM
!      use dcbexp_mod, only: IPEFN, LEXPP, MAXEXP, NEXPP
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbexp.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND)
      INTEGER :: IFIRST, ISECND, ITHIRD
C
      ITYP = 0
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      DO 100 IATOM = 1, NUCIND
      IF (DOATOM(IATOM)) THEN
        DO 50 IREP = 0, MAXREP
          DO 200 ICOOR = 1, 3
            ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
            IF(ISCOOR .GT. 0) THEN
              ITYP = ITYP + 1
              IFIRST=ISCOOR/100
              ISECND=MOD(ISCOOR,100)/10
              ITHIRD = MOD(MOD(ISCOOR,100),10)
              PLABEL(1) = 'NEF '//CHRNOS(IFIRST)
     &                    //CHRNOS(ISECND)//CHRNOS(ITHIRD)//' '
              PNAME    = 'NEF '//CHRNOS(IFIRST)
     &                    //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                //NAMN(IATOM)//REP(IREP)//CHRXYZ(ICOOR)
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND(CHRXYZ(ICOOR)//'EFN',
     &                 IPEFN(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
            END IF
 200      CONTINUE
 50     CONTINUE
      END IF
 100  CONTINUE
C
      IF(IPRINT.GE.1) THEN
        DO I = 1,ITYP
          INDXPR = LEXPP(IPEFN(I))
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_EFG */
      SUBROUTINE DEF_EFG(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operator for electric field gradient
C
C     Written by T.Saue Feb 19 1997
C     Last revision: J. Thyssen - Oct 13 1998
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     EFGINT: INTTYP = 30 (cart. el. field grad. at indiv. nuclei)
C     (-30: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-30)
      CALL DEF_EFG1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_efg */
      SUBROUTINE DEF_EFG1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define electric field gradient
C     Written by T.Enevoldsen - Sep 1996
C     Last revision: J. Thyssen - Oct 13 1998

C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbxpr_mod, only: PRPNAM
!      use dcbexp_mod, only: IPEFG, LEXPP, MAXEXP, NEXPP
      use nuclei_mod, only: ISTBNU, NUCIND
      use symmet_mod, only: ISYMAX, MAXREP
      use chrnos_mod, only: CHRNOS
      use chrxyz_mod, only: CHRXYZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: D1 = 1.0D0
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND)
#include "dcbxpr.h"
#include "dcbexp.h"
C
C
      ITYP = 0
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
         DO 200 ICOOR1 = 1, 3
         DO 200 ICOOR2 = ICOOR1, 3
            ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IAND(ISTBNU(IATOM),IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                  IOFF = IOFF + 1
                  ITYP = ITYP + 1
                  PNAME        = 'EFG: '//CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)
     *                 //CHRNOS(IATOM/10)
     &                 //CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  PLABEL(1) = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
                  CALL OP1IND(CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'
     *                 ,IPEFG(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      IF(IPRINT.GE.1) THEN
        DO I = 1,ITYP
          INDXPR = LEXPP(IPEFG(I))
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_PVC */
      SUBROUTINE DEF_PVC(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operator for parity violation - chirality
C
C     Written by J. Thyssen - Oct 7 1999
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     PVC: INTTYP = 62
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)
      CALL DEF_PVC1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_pvc */
      SUBROUTINE DEF_PVC1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define operator for parity violation - chirality
C
C     Written by J. Thyssen - Oct 7 1999
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use codata_mod, only: GFAC
      use pi_mod, only: PI
!      use dcbxpr_mod, only: PRPNAM
!      use dcbexp_mod, only: IPPVC, LEXPP, MAXEXP, NEXPP
      use nuclei_mod, only: GAUNUC, ISTBNU, NAMDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTNUC, MAXREP
      use chrnos_mod, only: CHRNOS
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
#include "codata.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0, D3 = 3.0D0, D4 = 4.0D0
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND)
#include "dcbxpr.h"
#include "dcbexp.h"
C
C
      ITYP = 0
      IPTYP     = 9
      NPCOMP    = 1
      IF(GAUNUC) THEN
        PFAC(1) = D1
      ELSE
        PFAC(1) = D3/(D4*GFAC*PI)
      ENDIF
      DO 100 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 200 IREP = 0, MAXREP
               IF (IAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
                  ITYP = ITYP + 1
                  PNAME = 'PVC: '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                    //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                    CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
C....             Gaussian nucleus
                  IF(GAUNUC) THEN
                    PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
C....             Point charge
                  ELSE
                    PLABEL(1) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  ENDIF
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
                  CALL OP1IND('PVC   '
     *                 ,IPPVC(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
C
      IF(IPRINT.GE.1) THEN
        DO I = 1,ITYP
          INDXPR = LEXPP(IPPVC(I))
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_EFT */
      SUBROUTINE DEF_EFT(PNAME,PLABEL,PFAC,TRLESS,IPRINT)
C***********************************************************************
C
C     Define operator for electric field third derivatives
C
C     TRLESS: traceless or non-traceless
C
C     Written by: J. Thyssen - Sept 30 1999
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   TRLESS
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     EFTINT: INTTYP = 52 or INTTYP = 53
      IF (TRLESS) THEN
         CALL SETATM(DOATOM,NATOM,52)
      ELSE
         CALL SETATM(DOATOM,NATOM,53)
      END IF
      CALL DEF_EFT1(PNAME,PLABEL,PFAC,DOATOM,TRLESS,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_efg */
      SUBROUTINE DEF_EFT1(PNAME,PLABEL,PFAC,DOATOM,TRLESS,IPRINT)
C***********************************************************************
C
C     Define operator for fourth derivative nuclear attraction
C     integrals.
C
C     Written by: J. Thyssen - Sept 30 1999
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbxpr_mod, only: PRPNAM
!      use dcbexp_mod, only: IPEFT, IPEFTNTL, LEXPP, MAXEXP, NEXPP
      use nuclei_mod, only: ISTBNU, NUCIND
      use symmet_mod, only: ISYMAX, MAXREP
      use chrnos_mod, only: CHRNOS
      use chrxyz_mod, only: CHRXYZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: D1 = 1.0D0
      CHARACTER PNAME*16, PLABEL(3)*8,
     &          PNPREFIX*6, PLBNAME*1
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND),TRLESS
#include "dcbxpr.h"
#include "dcbexp.h"
Ctec
Ctec
      ITYP = 0
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      IF (TRLESS) THEN
         PNPREFIX = 'EFT'
         PLBNAME  = 'V'
      ELSE
         PNPREFIX = 'EFTNTL'
         PLBNAME  = 'v'
      END IF
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
         DO 200 ICOOR1 = 1, 3
         DO 201 ICOOR2 = ICOOR1, 3
         DO 202 ICOOR3 = ICOOR2, 3
         DO 203 ICOOR4 = ICOOR3, 3
            ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            ISYMIJ = IEOR(ISYMIJ          ,ISYMAX(ICOOR3,1))
            ISYMIJ = IEOR(ISYMIJ          ,ISYMAX(ICOOR4,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IAND(ISTBNU(IATOM),IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                  IOFF = IOFF + 1
                  ITYP = ITYP + 1
                  PNAME        = PNPREFIX//
     &               CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//
     &               CHRXYZ(ICOOR3)//CHRXYZ(ICOOR4)//
     &               CHRNOS(IATOM/10)//
     &               CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  PLABEL(1) =
     &                 CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//
     &                 CHRXYZ(ICOOR3)//CHRXYZ(ICOOR4)//
     &                 PLBNAME//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
                  IF (TRLESS) THEN
                     CALL OP1IND(
     &                  CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//
     &                  CHRXYZ(ICOOR3)//CHRXYZ(ICOOR4)//'V',
     &                  IPEFT(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
                  ELSE
                     CALL OP1IND(
     &                  CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//
     &                  CHRXYZ(ICOOR3)//CHRXYZ(ICOOR4)//'v',
     &                  IPEFTNTL(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
                  END IF
               END IF
 300        CONTINUE
 203     CONTINUE
 202     CONTINUE
 201     CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      IF(IPRINT.GE.1) THEN
         DO I = 1,ITYP
            IF (TRLESS) THEN
               INDXPR = LEXPP(IPEFT(I))
            ELSE
               INDXPR = LEXPP(IPEFTNTL(I))
            END IF
            CALL PRSYMB(LUPRI,'.',75,0)
            WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &        'Component no.',I,':',PRPNAM(INDXPR)
            CALL PRSYMB(LUPRI,'.',75,0)
            CALL WRIXPR(I,INDXPR)
         ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_SPNSPN */
      SUBROUTINE DEF_SPNSPN(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operator for indirect spin-spin coupling
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
      CALL DEF_SPNSP1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_SPNSP1 */
      SUBROUTINE DEF_SPNSP1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define operator for indirect spin-spin coupling
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NAMN, NUCDEG, NUCDEP, NUCIND
!      use pgroup_mod, only: REP
      use symmet_mod, only: IPTCNT, ISYMAX, MAXREP
!      use dcbgen_mod, only: CVAL
!      use dcbprp_mod, only: ABUND, IPSPNSPN
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT, TRIAB
      use chrxyz_mod, only: CHRXYZ
      use chrnos_mod, only: CHRNOS
#include "implicit.h"
      INTEGER :: I
C
!#include "dgroup.h"
#include "pgroup.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
      LOGICAL DOATOM(NUCIND),SAME
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3),IOP(3)
C
C
      NSCOOR = 3*NUCDEP
      CALL IZERO(IPSPNSPN,2*NSCOOR)
C
C     Here we set up the operators needed for the linear reponse
C
C     Define left-hand operators: (-c) i.alpha x r/r^3
C     ================================================
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
      PFAC(1)  = -CVAL
      PFAC(2)  = -CVAL
      NPCOMP   = 2
C
C     Irreps
C
      DO 100 IREPO = 0, MAXREP
C
C       Atoms
C
        DO 200 IATOM1 = 1, NUCIND
        IF(DOATOM(IATOM1)) THEN
          MXATM2 = NUCIND
          IF (TRIAB) MXATM2 = IATOM1
C
C         Cartesian directions
C
          DO 300 ICOOR1 = 1, 3
            ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
            IF (ISCOR1 .GE. 1000) THEN
              CALL QUIT(
     &'DEF_SPNSP1: Max 1000 coords in this implementation !')
            ENDIF
            IF (ISCOR1 .GT. 0) THEN
              IFIRST = ISCOR1/100
              ISECND = MOD(ISCOR1,100)/10
              ITHIRD = MOD(MOD(ISCOR1,100),10)
              PNAME    = 'NUCSP'//CHRNOS(IFIRST)
     &                   //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                   //NAMN(IATOM1)//REP(IREPO)//CHRXYZ(ICOOR1)
              IPTYP    = IOP(ICOOR1)
              IC1      = MOD(ICOOR1+1,3) + 1
              IC2      = MOD(ICOOR1  ,3) + 1
              IREP1    = IEOR(IREPO,ISYMAX(IC2,1))
              IREP2    = IEOR(IREPO,ISYMAX(IC1,1))
              IS1      = IPTCNT(3*(IATOM1 - 1) + IC1,IREP1,1)
              IS2      = IPTCNT(3*(IATOM1 - 1) + IC2,IREP2,1)
              IF((IS1.LE.0).OR.(IS2.LE.0))
     &          WRITE(LUPRI,'(A,5I5)')
     &          'DEF_SPNSPN: sym. error !',IATOM1,ICOOR1,IC1,IC2
              IFIRST = IS1/100
              ISECND = MOD(IS1,100)/10
              ITHIRD = MOD(MOD(IS1,100),10)
              PLABEL(1)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              IFIRST = IS2/100
              ISECND = MOD(IS2,100)/10
              ITHIRD = MOD(MOD(IS2,100),10)
              PLABEL(2)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND('XSPNA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
              IPSPNSPN(1,ISCOR1) = INDOP1
C
C             Define right-hand operators: (-c) i.alpha x r/r^3
C             =================================================
C
              DO 400 IATOM2 = 1, MXATM2
              IF (DOATOM(IATOM2)) THEN
                SAME = IATOM1.EQ.IATOM2
                IF(SAME.AND.NUCDEG(IATOM1).EQ.1) GOTO 400
                SAME = TRIAB .AND. SAME
                MXCR2 = 3
                IF (SAME) MXCR2 = ICOOR1
                DO 500 ICOOR2 = 1, MXCR2
                  ISCOR2 = IPTCNT(3*(IATOM2-1)+ICOOR2,IREPO,2)
                  IF (ISCOR2 .GT. 0) IPSPNSPN(2,ISCOR2) = 1
  500           CONTINUE
              END IF
  400       CONTINUE
            ENDIF
  300     CONTINUE
        END IF
  200   CONTINUE
  100 CONTINUE
C
C
C     Define right-hand operators: (-c) i.alpha x r/r^3
C     =================================================
C
      DO I = 1,NSCOOR
        IF(IPSPNSPN(2,I).GT.0) THEN
#if defined (XH4)
          IF((I.EQ.1).OR.(I.EQ.5).OR.(I.EQ.6)) THEN
#endif
          INDOP1 = IPSPNSPN(1,I)
          INDXPR = LLRAPU(INDOP1)
          CALL OP1IND('XSPNB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
          IPSPNSPN(2,I) = INDOP1
#if defined (XH4)
          END IF
#endif
        ENDIF
      ENDDO
C
C     Print section
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' NMR indirect spin-spin coupling:'
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(TRIAB) THEN
        WRITE(LUPRI,'(A)')
     &    '*** WARNING ! Assuming symmetry J(AB) = J(BA)'
      ENDIF
      WRITE (LUPRI,'(A,F9.4,A)')
     &   ' Spin-spin couplings printed for atoms with '//
     &   'abundance greater than', ABUND,' %'
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(A)') ' ** A operators **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOOR
        INDOP1 = IPSPNSPN(1,I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRAPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(A)') ' ** B operators **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOOR
#if defined (XH4)
      IF((I.EQ.1).OR.(I.EQ.5).OR.(I.EQ.6)) THEN
#endif
        INDOP1 = IPSPNSPN(2,I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
#if defined (XH4)
      END IF
#endif
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prpdef */
      SUBROUTINE PRPDEF(WORK,LWORK)
C***********************************************************************
C
C     Set various properties
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
!      use dcbprp_mod, only: BXLAO, BYLAO, BZLAO, DIPOLE, DOEXP, DSO, EFFDE2, EFFDEN, EFG, EFN, EFT, EFTNTL, EPOLAR, ESRGTENS, ESR_HFCC, IPRPRP, LONDON, MSUSCDIA, MSUSCP, NQCC, NSTDIA, OPTROT, PVC, PVC_EFG, PVC_SHIELD, PVC_SPINRO, PVC_SPINSPIN, QUADRU, RHONUC, ROTG, SHIELD, SPINRO, SPNSPN
!      use dcbxlr_mod, only: XLR_SKIPEP
!      use dcbexp_mod, only: IPDIPOLE, IPREXP, LEXPP, MAXEXP, NEXPP
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
      use cbiher_mod, only: DSO, QUADRU
#include "implicit.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3), INDDIP(3)
#include "dcbprp.h"
#include "dcbxlr.h"
#include "dcbexp.h"
C
C Make sure that the diamagnetic parts are calculated if SKIPEP is true
C
      IF (XLR_SKIPEP) THEN
C      .... NMR shieldings
#ifdef MOD_LAO_REARRANGED
         IF (SHIELD.OR.shielding_rearrange) THEN
#else
         IF (SHIELD) THEN
#endif
            DOEXP = .TRUE.
            NSTDIA = .TRUE.
         END IF
C      .... spin-spin coupling
         IF (SPNSPN) THEN
            DOEXP = .TRUE.
            DSO = .TRUE.
         END IF
C      .... magnetic susceptibilities
         IF (MSUSCP) THEN
            DOEXP = .TRUE.
            MSUSCDIA = .TRUE.
         END IF
      END IF
C
C     Dipole moment
C
      IF(DIPOLE) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Dipole moment'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_DIPOLE(INDDIP,IPRPRP)
        CALL OP1IND('XDIPOL',IPDIPOLE(1),LEXPP,NEXPP,INDDIP(1),MAXEXP)
        CALL OP1IND('YDIPOL',IPDIPOLE(2),LEXPP,NEXPP,INDDIP(2),MAXEXP)
        CALL OP1IND('ZDIPOL',IPDIPOLE(3),LEXPP,NEXPP,INDDIP(3),MAXEXP)
      ENDIF
C
C     Quadrupole moment
C
      IF(QUADRU) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Quadrupole moments'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_QUADRU(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Electric fields at nuclei
C
      IF(EFN) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electric field at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFN(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Electric field gradients at nuclei
C
      IF(EFG) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electric field gradients'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFG(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Electronic density at nuclei
C
      IF(RHONUC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electronic density at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_RHONUC(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Parity violation - Chirality
C
      IF(PVC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Parity violation - Chirality'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_PVC(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Non-traceless electric field third derivatives
C
      IF(EFTNTL) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &     'Non-traceless electric field third derivatives'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFT(PNAME,PLABEL,PFAC,.FALSE.,IPRPRP)
      ENDIF
C
C     Traceless electric field third derivatives
C
      IF(EFT) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &     'Traceless electric field third derivatives'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFT(PNAME,PLABEL,PFAC,.TRUE.,IPRPRP)
      ENDIF
C
C     Nuclear quadrupole coupling
C
      IF(NQCC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Nuclear quadrupole moments'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFG(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Electric polarizabilities
C
      IF(EPOLAR) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electric dipole polarizabilities'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EPOLAR(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Diamagnetic spin-orbit
C
      IF(DSO) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &    'NMR diamagnetic spin-orbit contribution:'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_DSO(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Indirect spin-spin coupling
C
      IF(SPNSPN) THEN
        CALL DEF_SPNSPN(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Magnetizabilities (magnetic susceptibilities), relaxation,expectation value,
C     and if wished, diamagnetic terms
C
      IF(MSUSCP) THEN
        CALL DEF_MSUSCP()
      ENDIF
C
C     Nuclear shielding or ESR g-tensor or ESR anisotropy
C
#ifdef MOD_LAO_REARRANGED
      IF (SHIELD .OR. shielding_rearrange .OR.
     &    ESRGTENS .OR. ESR_HFCC) THEN
#else
      IF (SHIELD .OR. ESRGTENS .OR. ESR_HFCC) THEN
#endif
        CALL DEF_SHIELD(IPRPRP)
      ENDIF
C
C     Diamagnetic contribution to shieldings
C
      IF(NSTDIA) THEN
        CALL DEF_NSTDIA(PNAME,PLABEL,PFAC,IPREXP)
#ifdef MOD_LAO_REARRANGED
      ELSE IF (SHIELD .OR. shielding_rearrange) THEN
#else
      ELSE IF (SHIELD) THEN
#endif
      ! message has relevance if wanted NMR shieldings
        WRITE(LUPRI,'(2X,A)')
     &  '...Aucars diamagnetic'//
     &  ' contribution to NMR shieldings not calculated.'
      ENDIF
C
C     Optical rotation
C
      IF (OPTROT) THEN
        CALL DEF_OPTROT(IPRPRP)
      ENDIF
C
C     Molecular rotational g-tensor
C
      IF(ROTG) THEN
         IF(LONDON) CALL QUIT('.ROTG only implemented for CGO')
         CALL DEF_ROTG(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Nuclear spin-rotation constant
C
      IF(SPINRO) THEN
         CALL DEF_SPINRO(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     PV contribution to NMR shielding
C
      IF(PVC_SHIELD) THEN
         IF(LONDON) CALL QUIT('.PVCSHI only implemented for CGO')
         CALL DEF_PVC_SHIELD(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     PV contribution to spin-spin coupling
C
      IF(PVC_SPINSPIN) THEN
         CALL DEFINE_PVC_SPINSPIN()
      ENDIF
C
C     PV contribution to nuclear spin-rotation constant
C
      IF(PVC_SPINRO) THEN
         CALL DEF_PVC_SPINRO(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     PV contribution to EFG
C
      IF(PVC_EFG) THEN
        CALL DEF_PVC_EFG(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
! gosia:
! calculate << B_i; B_i >>_0 with LAO
! only for one component of magnetic field, i = x, y or z
! i need it temporarily in visual
      IF (BXLAO .OR. BYLAO .OR. BZLAO) THEN
        CALL DEF_BLAO(WORK,LWORK,IPRPRP)
      ENDIF
C
C     Effective electronic densities at the nuclei - new code
C
      IF(EFFDEN) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Effective electronic density at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFFDEN(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF
C
C     Effective electronic densities at the nuclei - old code
C
      IF(EFFDE2) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Effective electronic density at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFF_DEN2(PNAME,PLABEL,PFAC,IPRPRP)
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prpout */
      SUBROUTINE PRPOUT(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output of standardized properties
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
!      use dcbexp_mod, only: NEXPP
!      use dcbprp_mod, only: DIPOLE, DSO, EFG, EFN, EFT, EPOLAR, MSUSCP, NQCC, OPTROT, PVC, PVC_EFG, PVC_SHIELD, PVC_SPINRO, PVC_SPINSPIN, QUADRU, ROTG, SHIELD, SPINRO, SPNSPN
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
#ifdef HAS_PELIB
        use pelib_dirac_interface, only: use_pelib, pelib_ifc_gspol,
     &                                   pelib_ifc_do_lf
#endif
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
      use cbiher_mod, only: DSO, QUADRU
#include "implicit.h"
      DIMENSION EXPVAL(*),ABLRF(*),ABCNV(*),WORK(*)
#include "dcbexp.h"
#include "dcbprp.h"
      CALL PRSYMB(LUPRI,'-',76,0)
C
      KFRSAV=KFREE
C
C     Memory Allocation
C     -----------------
C
      CALL MEMGET('REAL',KEXP_CORR,NEXPP,WORK,KFREE,LFREE)
C
      CALL DZERO(WORK(KEXP_CORR),NEXPP)
C
C
C     Dipole moment
C
      IF(DIPOLE) THEN
        CALL PRP_DIPOLE(EXPVAL,WORK(KEXP_CORR),.FALSE.,WORK,KFREE,LFREE)
      ENDIF
C
C     Quadrupole moment
C
      IF(QUADRU) THEN
        CALL PRSYMB(LUPRI,'-',76,0)
C        CALL PRP_QUADRU(EXPVAL,WORK,LWORK)
      ENDIF
C
C     Parity violation: Chirality
C
      IF(PVC) THEN
        CALL PRP_PVC(EXPVAL,WORK,KFREE,LFREE)
      ENDIF
C
C     Electric fields at individual nuclei
C
      IF(EFN) THEN
        CALL PRP_EFN(EXPVAL,WORK(KEXP_CORR),.FALSE.,WORK,KFREE,LFREE)
      ENDIF
C
C     Electric field gradients
C
      IF(EFG) THEN
        CALL PRP_EFG(EXPVAL,WORK(KEXP_CORR),.FALSE.,WORK,KFREE,LFREE)
      ENDIF
C
C     non-traceless EFT's (electric field third derivatives)
C
      IF(EFT) THEN
        CALL PRP_EFT(EXPVAL,WORK(KEXP_CORR),.FALSE.,.FALSE.,WORK(KFREE),
     &               LFREE)
      ENDIF
C
C     Traceless EFT's (electric field third derivatives)
C
      IF(EFT) THEN
        CALL PRP_EFT(EXPVAL,WORK(KEXP_CORR),.TRUE.,.FALSE.,WORK(KFREE),
     &               LFREE)
      ENDIF
C
C     Nuclear quadrupole coupling
C
      IF(NQCC) THEN
        CALL PRP_NQCC(EXPVAL,WORK(KEXP_CORR),.FALSE.,WORK,KFREE,LFREE)
      ENDIF
C
C
      IF(EPOLAR) THEN
        CALL PRP_EPOLAR(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
C
C     Indirect spin-spin coupling
C
      IF(SPNSPN) THEN
        CALL PRP_SPNSPN(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
C
C     NMR diamagnetic spin-orbit
C
      IF(DSO) THEN
        CALL PRP_DSO(EXPVAL,WORK,KFREE,LFREE)
      ENDIF
C
C     NMR shielding
C
#ifdef MOD_LAO_REARRANGED
      IF(SHIELD.or.shielding_rearrange) THEN
#else
      IF(SHIELD) THEN
#endif
        CALL PRP_SHIELD(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
C
C     Magnetic susceptibilities
C
      IF(MSUSCP) THEN
        CALL PRP_MSUSCP(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
C
C     Nuclear spin-rotation
C
      IF(SPINRO) THEN
        CALL PRP_SPINRO(EXPVAL,ABLRF,ABCNV)
      ENDIF
C
C     Molecular rotational g-tensor
C
      IF(ROTG) THEN
        CALL PRP_ROTG(ABLRF,ABCNV)
      ENDIF
C
C     Optical rotation
C
      IF(OPTROT) THEN
        CALL PRP_OPTROT(ABLRF,ABCNV)
      ENDIF
C
C     PV contribution to NMR shielding
C
      IF(PVC_SHIELD)
     &  CALL PRP_PVC_SHIELD(ABLRF,ABCNV)
C
C     PV contribution to spin-spin coupling
C
      IF(PVC_SPINSPIN)
     &  CALL PRINTOUT_PVC_SPINSPIN(ABLRF,ABCNV,WORK,KFREE,LFREE)
C
C     PV contribution to nuclear spin-rotation constant
C
      IF(PVC_SPINRO)
     &  CALL PRP_PVC_SPINRO(ABLRF,ABCNV)
C
C     PV contribution to EFG
C
      IF(PVC_EFG)
     &  CALL PRP_PVC_EFG(ABLRF,ABCNV)
C
#ifdef HAS_PELIB
      if (use_pelib()) then
         write(lupri,'(/A)')
     &      'Polarizable environment included using PElib'//
     &      ' (https://gitlab.com/pe-software/pelib-public).'
         if (pelib_ifc_gspol()) then
            write(lupri,'(/A)')
     &      'Response properties calculated with dynamical'//
     &      ' polarizable embedding disabled (.GSPOL)'
         else
            write(lupri,'(/A)')
     &      'Response properties calculated with full inclusion of'//
     &      ' polarizable environment.'
         endif
         if (pelib_ifc_do_lf()) write(lupri,'(/A)')
     &   'Effective external field (EEF) correction applied.'
      endif
#endif
C
C     Memory deallocation
C     -------------------
C
      CALL MEMREL('PRPOUT',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prp_dipol */
      SUBROUTINE PRP_DIPOLE(EXPVAL,EXP_CORR,CORR,WORK,KFREE,LFREE)
C*****************************************************************************
C
C      Calculate dipole moment
C
C      Written by T.Saue - May 30 1996
C      Revision July 8 1996 - tsaue
C      Last revision Joost van Stralen - June 11 2002
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use codata_mod, only: DEBYE
      use mxcent_mod, only: MXCENT, MXCOOR
      use dipole_mod, only: DIPME, DIPMN
!      use dcbexp_mod, only: IPDIPOLE
#include "implicit.h"
      INTEGER :: I
#include "codata.h"
#include "dcbexp.h"
C
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL CORR
C
C     Electronic contribution to dipole moment
C     ========================================
C
      DO I = 1,3
        DIPME(I) = EXPVAL(IPDIPOLE(I))
      ENDDO
C
C     Nuclear contribution to dipole moment
C     ========================================
C
      CALL DIPNUC(0,.FALSE.)
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(A/)')
     &     '* Dipole moment:'
      WRITE(LUPRI,'(10X,3(A12,13X))')
     &      'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(10X,3(A12,13X))')
     &      'contribution','contribution','contribution'
      CALL PRSYMB(LUPRI,'-',76,0)
      IOFF = ICHAR('w')
      ! print in debye
      DO I = 1,3
         WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &        CHAR(IOFF+I),(DEBYE*DIPME(I)),' Debye',
     &                   (DEBYE*DIPMN(I)),' Debye',
     &        (DEBYE*(DIPME(I)+DIPMN(I))),' Debye'
      ENDDO
      CALL PRSYMB(LUPRI,'-',76,0)
      ! print in a.u.
      DO I = 1,3
         WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &        CHAR(IOFF+I),(DIPME(I)),' a.u. ',
     &                   (DIPMN(I)),' a.u. ',
     &        (DIPME(I)+DIPMN(I)),' a.u.'
      ENDDO
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(/A,F12.8,A)') ' 1 a.u. = ',DEBYE,' Debye'
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRP_NQCC */
      SUBROUTINE PRP_NQCC(EXPVAL,EXP_CORR,CORR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. Thyssen - Oct 13 1998
C     Last revision by Joost van Stralen - June 11 2002
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
!#include "dcbprp.h"
!#include "dcbexp.h"
!#include "dcbxlr.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      logical, allocatable :: DOATOM(:)
      LOGICAL CORR
C
      KFRSAV = KFREE
      allocate(DOATOM(NUCIND))
      CALL MEMGET('REAL',KAXIS,9*MXCENT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,9*MXCENT,WORK,KFREE,LFREE)
C
      CALL SETATM(DOATOM,NATOM,-30)
      CALL PRP_NQCC1(EXPVAL,EXP_CORR,CORR,DOATOM,WORK(KAXIS),
     &               WORK(KTEMP))
      deallocate(DOATOM)
      CALL MEMREL('PRP_NQCC',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prp_nqcc */
      SUBROUTINE PRP_NQCC1(EXPVAL,EXP_CORR,CORR,DOATOM,AXIS,TEMP)
C***********************************************************************
C
C      Calculate nuclear quadrupole coupling constants
C
C      Written by T.Enevoldsen  - Sep 1996
C      Last revision by Joost van Stralen - June 11 2002
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: EFG, NQCCCART
      use nuclei_mod, only: IZATOM, NAMDEP, NUCDEG, NUCIND
!      use nqcc_mod, only: EFGEL, EFGNUC, EFGTOT
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, J
C
#include "dcbprp.h"
!#include "dcbexp.h"
      REAL(dp), PARAMETER :: CONV = 234.9647417D0
      DIMENSION EXPVAL(*),AXIS(*),TEMP(*),PVAL(9),TMP(3),ITMP(3)
      DIMENSION EXP_CORR(*)
      LOGICAL   DOATOM(*), CORR
      CHARACTER*4 HZUNIT
#include "nqcc.h"
C
      IF (.NOT.EFG) THEN
C
C     Electronic contribution to electric field gradient
C     ==================================================
C
         CALL EL_EFG(0,EXPVAL,DOATOM,TEMP)
C
C     Nuclear contribution to electric field gradient
C     ===============================================
C
         CALL NUC_EFG(0,DOATOM)
C
      ENDIF
C
C print section
C
      IOFF = ICHAR('w')
      CALL PRSYMB(LUPRI,'-',76,0)
      IF (.NOT.NQCCcart) THEN
         WRITE(LUPRI,'(A/)')
     &'* Nuclear quadrupole coupling constants (principal components):'
      ELSE IF (NQCCcart) THEN
         WRITE(LUPRI,'(A/)')
     &   '* Nuclear quadrupole coupling constants'//
     &   '(cartesian axes components):'
      END IF
      WRITE(LUPRI,'(14X,A18)') 'Total contribution'
      CALL PRSYMB(LUPRI,'-',72,2)
      IATOMD = 0
      NATOM  = 0
      IATOMS = 0
      DO 100 IATOM = 1, NUCIND
      IF (DOATOM(IATOM)) THEN
        DO 110 IDEG = 1,NUCDEG(IATOM)
          IATOMD = IATOMD + 1
          NATOM = NATOM + 1
          DO I=1,3
            DO J=1,3
              EFGTOT(I,J,NATOM)=EFGNUC(I,J,NATOM)+EFGEL(I,J,NATOM)
            ENDDO
          ENDDO
C
C         Diagonalize to get principal components
C
          IJ = 1
          DO I = 1, 3
            DO J = 1, I
              PVAL(IJ) = EFGTOT(I,J,IATOM)
              IJ = IJ + 1
            ENDDO
          ENDDO

      IF (.NOT.NQCCcart) THEN
         CALL JACO(PVAL,AXIS,3,3,3,TMP,ITMP)
         EFGTOT(1,1,IATOM) = PVAL(1)
         EFGTOT(2,2,IATOM) = PVAL(3)
         EFGTOT(3,3,IATOM) = PVAL(6)
C
C        Reorder to calculate eta:
C           |xx|<|yy|<|zz|
C
         PVAL(2) = PVAL(3)
         PVAL(3) = PVAL(6)
         IZ = IDAMAX(3,PVAL(1),1)
         CALL DSWAP(1,PVAL(3),1,PVAL(IZ),1)
         IX = IDAMIN(3,PVAL(1),1)
         CALL DSWAP(1,PVAL(1),1,PVAL(IX),1)
         ETA = (PVAL(1) - PVAL(2))/ PVAL(3)
      ELSE IF (NQCCcart) THEN
         ETA = (PVAL(1) - PVAL(3))/ PVAL(6) !(xx-yy)/zz
      END IF

          NUCCHA = IZATOM(IATOM)
          DO ISO = 1, 5
            QMOM = DISOTP(NUCCHA,ISO,'QMOM')
            IF (QMOM .NE. 0.D0) THEN
              IATOMS = IATOMS + 1
              FACT=CONV*QMOM
C
C             Change to GHz when abs values bigger then 1.0D4
C
              IF ((ABS(EFGTOT(1,1,IATOM))*FACT.GT.1.0D4).OR.
     &            (ABS(EFGTOT(2,2,IATOM))*FACT.GT.1.0D4).OR.
     &            (ABS(EFGTOT(3,3,IATOM))*FACT.GT.1.0D4)) THEN
                FACT   = FACT * 1.0D-3
                HZUNIT = ' GHz'
              ELSE
                HZUNIT = ' MHz'
              ENDIF
              SPIN = DISOTP(NUCCHA,ISO,'SPIN')
              MASS = NINT(DISOTP(NUCCHA,ISO,'A'))
              WRITE (LUPRI,9000)
     &          'Nucleus: ',NAMDEP(IATOMD),
     &          'Proton charge:',NUCCHA,
     &          'Mass: ',MASS,'Spin: ',
     &           SPIN,'Quadrupole moment (b): ',QMOM
              DO I=1,3
                WRITE(LUPRI,'(2X,A3,6X,3(3X,F16.10,A))')
     &                'X'//CHAR(IOFF+I)//CHAR(IOFF+I),
     &                EFGTOT(I,I,IATOM)*FACT,HZUNIT
              ENDDO
              WRITE (LUPRI,'(2X,A3,9X,F16.8)') 'Eta', ETA
            ENDIF
          ENDDO
 110    CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
 100  CONTINUE
      IF(IATOMS.EQ.0) THEN
        WRITE(LUPRI,'(2X,A)')
     & 'No active atoms with quadrupole moment.'
      ENDIF
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(/A,F12.8,A)') '  1 a.u. = ',CONV,' MHz'
      RETURN
 9000 FORMAT (/2X,A9,A4,1X,A,I3,3X,A,I3,3X,A6,F4.1,3X,A23,F12.7,/)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck el_efg */
      SUBROUTINE EL_EFG(IPRINT,EXPVAL,DOATOM,TEMP)
C***********************************************************************
C
C     Calculates electronic contribution to electric field gradients
C
C     Written by T.Enevoldsen  - Sep 1996
C     Last revision
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbexp_mod, only: IPEFG
      use nuclei_mod, only: ISTBNU, NAMDEP, NUCDEG, NUCIND
      use symmet_mod, only: ISYMAX, MAXREP, PT
!      use nqcc_mod, only: EFGEL
#include "implicit.h"
      INTEGER :: I, J
C
#include "dcbexp.h"
      DIMENSION EXPVAL(*), TEMP(3,3,MXCENT)
      LOGICAL DOATOM(*)
#include "nqcc.h"
      CALL DZERO(EFGEL,9*MXCENT)
      CALL DZERO(TEMP,9*MXCENT)
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
         DO 200 ICOOR1 = 1, 3
            DO 210 ICOOR2 = ICOOR1, 3
               ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
               DO 300 IREPC = 0, MAXREP
                  IF (IAND(ISTBNU(IATOM),
     &                IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                     ITYP = ITYP + 1
                     TEMP(ICOOR1,ICOOR2,IATOM) =
     &                  TEMP(ICOOR1,ICOOR2,IATOM) + EXPVAL(IPEFG(ITYP))
     &                                             /NUCDEG(IATOM)
                  END IF
 300           CONTINUE
 210        CONTINUE
 200     CONTINUE
 100  CONTINUE
      NATOM =  0
      IATOMD = 0
      DO 110 IATOM = 1, NUCIND
      IF(DOATOM(IATOM)) THEN
        DO 310 IREPC = 0, MAXREP
        IF (IAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          DO 220 I = 1, 3
            DO 230 J = I, 3
              EFGEL(I,J,NATOM) = TEMP(I,J,IATOM)*
     &                           PT(IAND(IEOR(ISYMAX(I,1),
     &                           ISYMAX(J,1)),IREPC))
              IF (I .NE. J) EFGEL(J,I,NATOM) = EFGEL(I,J,NATOM)
 230        CONTINUE
 220      CONTINUE
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,A,A/)')
     &        'Electronic EFG Tensor of nucleus ',NAMDEP(IATOMD)
              CALL OUTPUT(EFGEL(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
            END IF
        END IF
 310    CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
 110  CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRP_EFG */
      SUBROUTINE PRP_EFG(EXPVAL,EXP_CORR,CORR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. Thyssen - Oct 13 1998
C     Last revision by Joost van Stralen - June 11 2002
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: IPRPRP
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
#include "dcbprp.h"
!#include "dcbexp.h"
!#include "dcbxlr.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL CORR
      logical, allocatable :: DOATOM(:)
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      allocate(DOATOM(NUCIND))
      CALL MEMGET('REAL',KAXIS,9*MXCENT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,9*MXCENT,WORK,KFREE,LFREE)
C
C     (-30: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-30)
      CALL PRP_EFG1(EXPVAL,EXP_CORR,DOATOM,WORK(KAXIS),WORK(KTEMP),
     &              CORR,IPRPRP)
      deallocate(DOATOM)
      CALL MEMREL('PRP_EFG',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRP_EFG1(EXPVAL,EXP_CORR,DOATOM,AXIS,TEMP,CORR,IPRINT)
C***********************************************************************
C
C     Calculate nuclear quadrupole coupling constants
C
C     Written by T.Enevoldsen  - Sep 1996
C     Revision by J. Thyssen - Oct 13 1998
C     Last revision by Joost van Stralen - June 11 2002
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbexp_mod, only: IPREXP
      use nuclei_mod, only: NAMDEP, NUCDEG, NUCIND
!      use nqcc_mod, only: EFGEL, EFGNUC, EFGTOT, EFG_METHOD
#include "implicit.h"
      INTEGER :: I, J
C
#include "dcbexp.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),AXIS(3,3,MXCENT),TEMP(*)
      DIMENSION TMP(3),ITMP(3),PVAL(9)
      LOGICAL   DOATOM(*), CORR
      CHARACTER FMTEFG*80, FMTEFGA*6, FMTEFGB*6, FMTEFGC*6
#include "nqcc.h"
C
C
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(A/)') '* Electric field gradients :'
C
      CALL EL_EFG(IPREXP,EXPVAL,DOATOM,TEMP)
      CALL NUC_EFG(IPREXP,DOATOM)
C
C     print section
      IOFF = ICHAR('w')
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(12X,A/)')
     &     'Individual (non-zero) components'
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      CALL PRSYMB(LUPRI,'-',72,2)
      NATOM  = 0
      IATOMD = 0
      DO IATOM = 1, NUCIND
      IF(DOATOM(IATOM)) THEN
        DO IDEG = 1,NUCDEG(IATOM)
          IATOMD = IATOMD + 1
          NATOM  = NATOM + 1
          WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(IATOMD)
          DO I = 1,3
            DO J = 1,3
              EFGTOT(I,J,NATOM)=EFGNUC(I,J,NATOM)+EFGEL(I,J,NATOM)
              EFG_METHOD(I,J,NATOM) = EFGTOT(I,J,NATOM)
              EFGPRINT1 = EFGEL(I,J,NATOM)
              EFGPRINT2 = EFGNUC(I,J,NATOM)
              EFGPRINT3 = EFGTOT(I,J,NATOM)
              IF (ABS(EFG_METHOD(I,J,NATOM)).GT.1.0D-12) THEN
                IF (ABS(EFG_METHOD(I,I,NATOM)).GT.1.0D+04) THEN
                   FMTEFG = '(2X,A3,3(3X,E16.10,A))'
                ELSE
                   FMTEFG = '(2X,A3,3(3X,F16.10,A))'
                ENDIF
                WRITE(LUPRI,FMTEFG)
     &            'q'//CHAR(IOFF+I)//CHAR(IOFF+J),EFGPRINT1,
     &            ' au ',EFGPRINT2,' au ',EFGPRINT3,' au '
              ENDIF
            END DO
          END DO
        END DO
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
      END DO
C
C           diagonalize total efg tensor
C
      WRITE(LUPRI,'(/,A)') '  '
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(12X,A)')
     &     'Total contribution to principal components '
      CALL PRSYMB(LUPRI,'-',72,2)
      IATOMD = 0
      NATOM  = 0
      DO IATOM = 1, NUCIND
      IF(DOATOM(IATOM)) THEN
        DO IDEG = 1,NUCDEG(IATOM)
          IATOMD = IATOMD + 1
          NATOM  = NATOM + 1
          CALL DUNIT(AXIS(1,1,NATOM),3)
          IJ = 1
          DO I = 1, 3
            DO J = 1, I
              PVAL(IJ) = EFGTOT(I,J,NATOM)
              IJ = IJ + 1
            ENDDO
          ENDDO
          CALL JACO(PVAL,AXIS(1,1,NATOM),3,3,3,TMP,ITMP)
          EFGTOT(1,1,NATOM) = PVAL(1)
          EFGTOT(2,2,NATOM) = PVAL(3)
          EFGTOT(3,3,NATOM) = PVAL(6)
          WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(IATOMD)
          DO I=1,3
            FMTEFG = '(2X,A3,3(3X,E16.10,A))'
            WRITE(LUPRI,FMTEFG)
     &            'q'//CHAR(IOFF+I)//CHAR(IOFF+I)
     &            ,EFGTOT(I,I,NATOM),' au '
          ENDDO
C
C         print section principle axis
C
          IF (IPRINT .GE. 4) THEN
            WRITE(LUPRI,'(A/)')
     &         '  Principal axis to total contribution'
            WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(IATOMD)
            CALL OUTPUT(AXIS(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
          END IF
        END DO
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(/A,1PE12.5,A,/)')
     &      '  1 a.u. = ',9.71736E21,' V*m**-2'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck prp_epolar */
      SUBROUTINE PRP_EPOLAR(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
C  *********************************************************************
C  *                                                                   *
C  *  electric dipole polarizability output                            *
C  *                                                                   *
C  *********************************************************************
C  *                            written                   revision     *
C  *  radovan bast              december 2004             07/02/2005   *
C  *                                                                   *
C  *********************************************************************
C
!      use codata_mod, only: XTANG
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
      use consts_mod, only: D0
!      use dcbprp_mod, only: DOVDW, IPEPOLAR
!      use dcbxlr_mod, only: BFREQ, MAXFLR, NBFREQ, NLRAPT, NLRBPT, THCLR, XLR_SKIPEP, XLR_UNCOUP
#include "implicit.h"
      INTEGER :: I, J
#include "codata.h"
C
      DOUBLE PRECISION DFAC
C
#include "dcbprp.h"
#include "dcbxlr.h"
C
      CHARACTER CNVTEXT(3)*17
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),TMPMAT(3,3),WORK(*)
C
      dimension average(maxflr)
C
C     conversion factor alpha/a.u. > 1.0D30*alpha/m**3
      DFAC = XTANG**3
C
      CALL AROUND('Electric dipole polarizability')
      WRITE(LUPRI,'(/A,F12.8,A/)') ' 1 a.u = ',DFAC,' angstrom**3'
C
C
C     print header
C     ============
C
      WRITE(LUPRI,'(/A/)')
     &   '@   Elements of the electric dipole polarizability tensor'
C
      DO IFREQ = 1,NBFREQ
        CALL DZERO(TMPMAT,9)
C
C
C       print frequency
C       ===============
C
        IF(.NOT.(NBFREQ.EQ.1)) THEN
          WRITE(LUPRI,'(/A,F15.8,A/)')
     &       '@  frequency', BFREQ(IFREQ), ' a.u.'
        ENDIF
C
C
C       fill temporary matrix
C       =====================
C
        ! edh: ABLRF appears to have different dim. for damped response?
        ! Check if this output still works...
        DO I = 1,3
          DO J = 1,3
            KI = IPEPOLAR(I,1)
            KJ = IPEPOLAR(J,1)
C
C           e-e
            TMPMAT(I,J) = -ABLRF(IFREQ,KI,KJ,2)
C
C           don't look for e-p in case of XLR_SKIPEP
            IF(.NOT.XLR_SKIPEP) THEN
C
C             e-p
              TMPMAT(I,J) = TMPMAT(I,J)
     &                    - ABLRF(IFREQ,KI,KJ,3)
            ENDIF
          ENDDO
        ENDDO
C
C
C       check if converged
C       ==================
C
        DO J = 1,3
          IF(XLR_UNCOUP) THEN
            CNVTEXT(J) = '(uncoupled calc.)'
          ELSEIF(ABCNV(IFREQ,J).LE.THCLR) THEN
            CNVTEXT(J) = '(converged)'
          ELSE
            CNVTEXT(J) = '(NOT converged!)'
          ENDIF
        ENDDO
C
C       calculate tensor trace
C       ==========================
C
        average(ifreq)=(tmpmat(1,1)+tmpmat(2,2)+tmpmat(3,3))/3.0d0
C
C
C       printout a.u.
C       =============
C
        WRITE(LUPRI,'(A3,A3,F22.8,A5,A20)') '@  ',
     &              'xx', TMPMAT(1,1), 'a.u.', CNVTEXT(1)
        WRITE(LUPRI,'(A3,A3,F22.8,A5,A20)') '@  ',
     &              'yy', TMPMAT(2,2), 'a.u.', CNVTEXT(2)
        WRITE(LUPRI,'(A3,A3,F22.8,A5,A20)') '@  ',
     &              'zz', TMPMAT(3,3), 'a.u.', CNVTEXT(3)
        WRITE(LUPRI,*)
        WRITE(LUPRI,'(A3,A8,F17.8,A5)') '@  ',
     &              'average', average(ifreq), 'a.u.'


!       calculate anisotropy
!       ====================

        anisotropy_squared = 0.0d0
        do i = 1, 3
          do j = 1, 3
              anisotropy_squared = anisotropy_squared
     &                         + (3.0d0/2.0d0)*tmpmat(i, j)*tmpmat(i, j)
     &                         - (1.0d0/2.0d0)*tmpmat(i, i)*tmpmat(j, j)
          end do
        end do
C
        if (abs(anisotropy_squared) .lt. 1.0D-12) then
           anisotropy = 0.0d0
        else
           anisotropy = anisotropy_squared**0.5d0
        end if
C
C       print anisotropy
C       ================
C
        WRITE(LUPRI,'(A3,A11,F9.3,A10)') '@  ',
     &              'anisotropy', anisotropy, 'a.u.'
C
C       printout SI
C       ===========
C
        WRITE(LUPRI,'(/A3,A3,F22.8,A12)') '@  ',
     &              'xx', TMPMAT(1,1)*DFAC, 'angstrom**3'
        WRITE(LUPRI,'(A3,A3,F22.8,A12)') '@  ',
     &              'yy', TMPMAT(2,2)*DFAC, 'angstrom**3'
        WRITE(LUPRI,'(A3,A3,F22.8,A12)') '@  ',
     &              'zz', TMPMAT(3,3)*DFAC, 'angstrom**3'
        WRITE(LUPRI,*)
        WRITE(LUPRI,'(A3,A8,F17.8,A12)') '@  ',
     &              'average', average(ifreq)*DFAC, 'angstrom**3'
        WRITE(LUPRI,'(A3,A11,F9.3,A17)') '@  ',
     &              'anisotropy', anisotropy*DFAC, 'angstrom**3'
        WRITE(LUPRI,*)
C
        IF(.NOT.(IFREQ.EQ.NBFREQ)) THEN
          WRITE(LUPRI,'(A)') '@   ---------'
        ENDIF
      ENDDO
C
      CALL PRSYMB(LUPRI,'-',46,0)
C
      if (dovdw) then
         call van_der_waals(nbfreq-1,average(2),work(kfree),lfree)
      end if
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prp_eft */
      SUBROUTINE PRP_EFT(EXPVAL,EXP_CORR,TRLESS,CORR,WORK,LWORK)
C***********************************************************************
C
C     Written by J. Thyssen - Oct 13 1998
C     Last revision: Joost van Stralen - June 7 2002
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCDEP, NUCIND
#include "implicit.h"
C
!#include "dcbprp.h"
!#include "dcbexp.h"
!#include "dcbxlr.h"
C
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL   TRLESS
      logical, allocatable :: DOATOM(:)
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
C     Allocate memory
C
      allocate(DOATOM(NUCIND))
      CALL MEMGET('REAL',KTEMP    ,15*NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP_CORR,15*NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KNUC     ,15*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEL      ,15*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCORR    ,15*NUCDEP,WORK,KFREE,LFREE)
C
      IF (TRLESS) THEN
         CALL SETATM(DOATOM,NATOM,52)
      ELSE
         CALL SETATM(DOATOM,NATOM,53)
      END IF
      CALL PRP_EFT1(EXPVAL,EXP_CORR,DOATOM,WORK(KTEMP),
     &              WORK(KTMP_CORR),WORK(KEL),WORK(KNUC),WORK(KCORR),
     &              TRLESS,CORR)
      deallocate(DOATOM)
      CALL MEMREL('PRP_EFT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prp_EFT */
      SUBROUTINE PRP_EFT1(EXPVAL,EXP_CORR,DOATOM,TEMPEFT,TMPEFT_CORR,
     &                    EFTEL,EFTNUC,EFTCORR,TRLESS,CORR)
C***********************************************************************
C
C     Written by J. Thyssen - Oct 19 1999
C     Last revision: Joost van Stralen - June 7 2002
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbexp_mod, only: IPEFT, IPEFTNTL
      use nuclei_mod, only: ISTBNU, NAMDEP, NUCDEG, NUCDEP, NUCIND
      use chrxyz_mod, only: CHRXYZ
      use symmet_mod, only: ISYMAX, MAXREP, PT
#include "implicit.h"
      INTEGER :: I
C
      PARAMETER (D0  =  0.0D00, D1 = 1.0D00, D3 = 3.0D00,
     &           D4  =  4.0D00, D7 = 7.0D00,D27=27.0D00,
     &           D35 = 35.0D00, DTHR = 1.0D-10)
C
!#include "dcbprp.h"
#include "dcbexp.h"
!#include "dcbxlr.h"
      LOGICAL   DOATOM(*), TRLESS, CORR
      DIMENSION EXPVAL(*),EXP_CORR(*),TEMPEFT(15,*),TMPEFT_CORR(15,*),
     &          EFTEL(15,*),EFTNUC(15,*),EFTCORR(15,*)
      DIMENSION ITR(3,6),DTRNUC(6),DTREL(6),DTRCORR(6)
      DATA ITR  / 1, 4, 6,
     &            2, 7, 9,
     &            3, 8,10,
     &            4,11,13,
     &            5,12,14,
     &            6,13,15/
C
C
      CALL QENTER('PRP_EFT')
C
      CALL PRSYMB(LUPRI,'-',76,0)
      IF (TRLESS) THEN
         WRITE(LUPRI,'(/A/)')
     &        '* Traceless electric field third derivative :'
      ELSE
         WRITE(LUPRI,'(/A/)')
     &        '* Non-traceless electric field third derivative :'
      END IF
C
C     **************************************
C     ***  Nuclear contribution to EFT:  ***
C     **************************************
C
      CALL NUCEFT(EFTNUC)
C
C     *****************************************
C     ***  Electronic contribution to EFT:  ***
C     *****************************************
C
C     Gather expectation values and take
C     symmetry combinations.
C
      CALL DZERO(TEMPEFT    ,15*NUCIND)
      CALL DZERO(EFTEL      ,15*NUCDEP)
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
c        IF (.NOT. DOATOM(IATOM)) GOTO 100
         IELEM = 0
         DO ICOOR1 = 1, 3
         DO ICOOR2 = ICOOR1, 3
         DO ICOOR3 = ICOOR2, 3
         DO ICOOR4 = ICOOR3, 3
            IELEM = IELEM + 1
            ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            ISYMIJ = IEOR(ISYMIJ          ,ISYMAX(ICOOR3,1))
            ISYMIJ = IEOR(ISYMIJ          ,ISYMAX(ICOOR4,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IAND(ISTBNU(IATOM),IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                  ITYP = ITYP + 1
                  IF (TRLESS) THEN
                     IPE = IPEFT(ITYP)
                  ELSE
                     IPE = IPEFTNTL(ITYP)
                  END IF
                  TEMPEFT(IELEM,IATOM) =
     &               TEMPEFT(IELEM,IATOM) +
     &               EXPVAL(IPE) / NUCDEG(IATOM)
               END IF
 300        CONTINUE
         END DO
         END DO
         END DO
         END DO
 100  CONTINUE
C
C     Assign value for individual centers
C
      NATOM = 0
      DO 110 IATOM = 1,NUCIND
c        IF (.NOT. DOATOM(IATOM)) GOTO 110
         DO IREPC = 0,MAXREP
            IF (IAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
               NATOM = NATOM + 1
               IELEM = 0
               DO ICOOR1 = 1, 3
               DO ICOOR2 = ICOOR1, 3
               DO ICOOR3 = ICOOR2, 3
               DO ICOOR4 = ICOOR3, 3
                  IELEM = IELEM + 1
                  ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
                  ISYMIJ = IEOR(ISYMIJ          ,ISYMAX(ICOOR3,1))
                  ISYMIJ = IEOR(ISYMIJ          ,ISYMAX(ICOOR4,1))
                  EFTEL(IELEM,NATOM) =
     &               TEMPEFT(IELEM,IATOM) * PT(IAND(ISYMIJ,IREPC))
               END DO
               END DO
               END DO
               END DO
            END IF
         END DO
 110  CONTINUE
C
C     *****************************************
C     ***  Output section 1: non-traceless  ***
C     *****************************************
C
      IF (TRLESS) GOTO 200
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(14X,A/)')
     &     'Non-traceless electric field third derivatives'
      WRITE(LUPRI,'(14X,A/)')
     &     'Individual (non-zero) components'
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &  'contribution','contribution','contribution'
      CALL PRSYMB(LUPRI,'-',72,2)
      DO IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(IATOM)
            IELEM = 0
            DO ICOOR1=1,3
            DO ICOOR2=ICOOR1,3
            DO ICOOR3=ICOOR2,3
            DO ICOOR4=ICOOR3,3
               IELEM = IELEM + 1
               EFTTOT = EFTNUC(IELEM,IATOM)+EFTEL(IELEM,IATOM)
               IF (ABS(EFTTOT).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A5,3(3X,1P,D16.9,A))')
     &                 'h'//CHRXYZ(-ICOOR1)//CHRXYZ(-ICOOR2)
     &                    //CHRXYZ(-ICOOR3)//CHRXYZ(-ICOOR4),
     &                 EFTEL (IELEM,IATOM),' au ',
     &                 EFTNUC(IELEM,IATOM),' au ',
     &                 EFTTOT             ,' au '
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDDO
C
      DO IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
C
C          Trace 1 = V_xxaa = V_xxxx + V_xxyy + V_xxzz
C          Trace 2 = V_xyaa = V_xxxy + V_xyyy + V_xyzz
C          Trace 3 = V_xzaa = V_xxxz + V_xyyz + V_xzzz
C          Trace 4 = V_yyaa = V_xxyy + V_yyyy + V_yyzz
C          Trace 5 = V_yzaa = V_xxyz + V_yyyz + V_yzzz
C          Trace 6 = V_zzaa = V_xxzz + V_yyzz + V_zzzz
C          (aa is short for sum over xx,yy, and zz)
C
           DO I = 1,6
               DTRNUC(I) = EFTNUC(ITR(1,I),IATOM) +
     &                     EFTNUC(ITR(2,I),IATOM) +
     &                     EFTNUC(ITR(3,I),IATOM)
               DTREL(I)  = EFTEL (ITR(1,I),IATOM) +
     &                     EFTEL (ITR(2,I),IATOM) +
     &                     EFTEL (ITR(3,I),IATOM)
            END DO
C
C           Subtract linear kombinations of traces...
C           Documented in J. Thyssen et al. (unpublished)
C
C           V_xxxx = V_xxxx - 27/35 V_xxaa +  3/35 V_yyaa +  3/35 V_zzaa
C           V_xxxy = V_xxxy -  3/ 7 V_xyaa
C           V_xxxz = V_xxxz -  3/ 7 V_xzaa
C           V_xxyy = V_xxyy -  4/35 V_xxaa -  4/35 V_yyaa +  1/35 V_zzaa
C           V_xxyz = V_xxyz -  1/ 7 V_yzaa
C           V_xxzz = V_xxzz -  4/35 V_xxaa +  1/35 V_yyaa -  4/35 V_zzaa
C           V_xyyy = V_xyyy -  3/ 7 V_xyaa
C           V_xyyz = V_xyyz -  1/ 7 V_xzaa
C           V_xyzz = V_xyzz -  1/ 7 V_xyaa
C           V_xzzz = V_xzzz -  3/ 7 V_xzaa
C           V_yyyy = V_yyyy +  3/35 V_xxaa - 27/35 V_yyaa +  3/35 V_zzaa
C           V_yyyz = V_yyyz -  3/ 7 V_yzaa
C           V_yyzz = V_yyzz +  1/35 V_xxaa -  4/35 V_yyaa -  4/35 V_zzaa
C           V_yzzz = V_yzzz -  3/ 7 V_yzaa
C           V_zzzz = V_yyyy +  3/35 V_xxaa +  3/35 V_yyaa - 27/35 V_zzaa
C
            EFTEL( 1,IATOM) = EFTEL( 1,IATOM) -
     &         D27/D35*DTREL(1) +  D3/D35*DTREL(4) +  D3/D35*DTREL(6)
            EFTEL( 2,IATOM) = EFTEL( 2,IATOM) - D3/D7*DTREL(2)
            EFTEL( 3,IATOM) = EFTEL( 3,IATOM) - D3/D7*DTREL(3)
            EFTEL( 4,IATOM) = EFTEL( 4,IATOM) +
     &          (-D4/D35)*DTREL(1) -  D4/D35*DTREL(4) +  D1/D35*DTREL(6)
            EFTEL( 5,IATOM) = EFTEL( 5,IATOM) - D1/D7*DTREL(5)
            EFTEL( 6,IATOM) = EFTEL( 6,IATOM) +
     &          (-D4/D35)*DTREL(1) +  D1/D35*DTREL(4) -  D4/D35*DTREL(6)
            EFTEL( 7,IATOM) = EFTEL( 7,IATOM) - D3/D7*DTREL(2)
            EFTEL( 8,IATOM) = EFTEL( 8,IATOM) - D1/D7*DTREL(3)
            EFTEL( 9,IATOM) = EFTEL( 9,IATOM) - D1/D7*DTREL(2)
            EFTEL(10,IATOM) = EFTEL(10,IATOM) - D3/D7*DTREL(3)
            EFTEL(11,IATOM) = EFTEL(11,IATOM) +
     &           D3/D35*DTREL(1) - D27/D35*DTREL(4) +  D3/D35*DTREL(6)
            EFTEL(12,IATOM) = EFTEL(12,IATOM) - D3/D7*DTREL(5)
            EFTEL(13,IATOM) = EFTEL(13,IATOM) +
     &           D1/D35*DTREL(1) -  D4/D35*DTREL(4) -  D4/D35*DTREL(6)
            EFTEL(14,IATOM) = EFTEL(14,IATOM) - D3/D7*DTREL(5)
            EFTEL(15,IATOM) = EFTEL(15,IATOM) +
     &           D3/D35*DTREL(1) +  D3/D35*DTREL(4) - D27/D35*DTREL(6)
         END IF
C
      END DO
C
C     *************************************
C     ***  Output section 2: traceless  ***
C     *************************************
C
 200  CONTINUE
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(14X,A)')
     &     '"traceless" electric field third derivatives',
     &     'after applying the Laplace equation: habcc = 0',
     &     '(einstein summation convention).',
     &     'See, for example:',
     &     ' A. D. McLean & M. Yoshimine, JCP 47, 1927 (1967)',
     &     ' J. Thyssen, P. Schwerdfeger, M. Bender, W. Nazarewicz,',
     &     '    and P. B. Semmes, PRA 63, 022505 (2001).'
      WRITE(LUPRI,'(14X,A/)')
     &     'Individual (non-zero) components'
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      CALL PRSYMB(LUPRI,'-',72,2)
      DO IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(IATOM)
            IELEM = 0
            DO ICOOR1=1,3
            DO ICOOR2=ICOOR1,3
            DO ICOOR3=ICOOR2,3
            DO ICOOR4=ICOOR3,3
               IELEM = IELEM + 1
               EFTTOT = EFTNUC(IELEM,IATOM)+EFTEL(IELEM,IATOM)
               IF (ABS(EFTTOT).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A5,3(3X,F16.10,A))')
     &                 'h'//CHRXYZ(-ICOOR1)//CHRXYZ(-ICOOR2)
     &                    //CHRXYZ(-ICOOR3)//CHRXYZ(-ICOOR4),
     &                 EFTEL (IELEM,IATOM),' au ',
     &                 EFTNUC(IELEM,IATOM),' au ',
     &                 EFTTOT             ,' au '
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDDO
C
C     Double checkt to see if Laplace equation is fullfilled.
C     It should be, or my formulas are wrong!
C
      DO IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
C
            DO I = 1,6
               DTRNUC(I) = EFTNUC(ITR(1,I),IATOM) +
     &                     EFTNUC(ITR(2,I),IATOM) +
     &                     EFTNUC(ITR(3,I),IATOM)
               DTREL(I)  = EFTEL (ITR(1,I),IATOM) +
     &                     EFTEL (ITR(2,I),IATOM) +
     &                     EFTEL (ITR(3,I),IATOM)
               IF (DTREL (I) .GE. DTHR .OR.
     &             DTRNUC(I) .GE. DTHR) THEN
                  WRITE (LUPRI ,'(/2X,A9,A4,/)')
     &               'Nucleus: ',NAMDEP(IATOM)
                  WRITE(LUPRI,'(2X,A,I1,A,1P,D20.10)')
     &               'Laplace equation ',I,' (elec. part): ',DTREL(I)
                  WRITE(LUPRI,'(2X,A,I1,A,1P,D20.10)')
     &               'Laplace equation ',I,' (nuc. part) : ',DTRNUC(I)
               END IF
            END DO
         END IF
C
C
      END DO
C
C
      CALL QEXIT('PRP_EFT')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck nuceft */
      SUBROUTINE NUCEFT(EFTNUC)
C***********************************************************************
C
C     Written by J. Thyssen - Nov 15 1999
C     Nuclear contribution to electric field third-derivatives.
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CHARGE, CORD, ISTBNU, NUCDEP, NUCIND
      use symmet_mod, only: ISYMAX, MAXOPR, PT
#include "implicit.h"
      INTEGER :: I, KK
C
      PARAMETER (D3 = 3.0D0, D9 = 9.0D00, D15 = 15.0D0,
     &           D45 = 45.0D00, D90 = 90.0D00, D105 = 105.0D00)
      DIMENSION GEOM(3*8*MXCENT), EFTNUC(15,*)
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
C
      CALL DZERO(EFTNUC,15*NUCDEP)
      KK = 1
      NATOM = 1
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               LL = 1
               DO 100 I2 = 1, NUCIND
                  IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
                  DO 110 ISYM = 0, MAXOPR
                     IF (IAND(ISYM,ISTBNU(I2)) .EQ. 0) THEN
                        IF (LL .NE. KK) THEN
                           XCOOR = GEOM(LL) - GEOM(KK)
                           YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                           ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                           R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                           R5 = R2*R2*SQRT(R2)
                           R7 = R5*R2
                           R9 = R7*R2
C
C                          V_xxxx
C
                           EFTNUC( 1,NATOM) = EFTNUC( 1,NATOM)
     &                        + CHARGE(I2) * (
     &                            D9 / R5
     &                        -  D90 * (XCOOR**2) / R7
     &                        + D105 * (XCOOR**4) / R9 )
C
C                          V_xxxy
C
                           EFTNUC( 2,NATOM) = EFTNUC( 2,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D45 * XCOOR*YCOOR / R7
     &                        + D105 * (XCOOR**3) * YCOOR / R9 )
C
C                          V_xxxz
C
                           EFTNUC( 3,NATOM) = EFTNUC( 3,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D45 * XCOOR*ZCOOR / R7
     &                        + D105 * (XCOOR**3) * ZCOOR / R9 )
C
C                          V_xxyy
C
                           EFTNUC( 4,NATOM) = EFTNUC( 4,NATOM)
     &                        + CHARGE(I2) * (
     &                            D3 / R5
     &                        -  D15 * XCOOR**2 / R7
     &                        -  D15 * YCOOR**2 / R7
     &                        + D105 * (XCOOR**2) * (YCOOR**2) / R9 )
C
C                          V_xxyz
C
                           EFTNUC( 5,NATOM) = EFTNUC( 5,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D15 * YCOOR * ZCOOR / R7
     &                        + D105 * (XCOOR**2) * YCOOR * ZCOOR / R9 )
C
C                          V_xxzz
C
                           EFTNUC( 6,NATOM) = EFTNUC( 6,NATOM)
     &                        + CHARGE(I2) * (
     &                            D3 / R5
     &                        -  D15 * XCOOR**2 / R7
     &                        -  D15 * ZCOOR**2 / R7
     &                        + D105 * (XCOOR**2) * (ZCOOR**2) / R9 )
C
C                          V_xyyy
C
                           EFTNUC( 7,NATOM) = EFTNUC( 7,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D45 * XCOOR*YCOOR / R7
     &                        + D105 * (YCOOR**3) * XCOOR / R9 )
C
C                          V_xyyz
C
                           EFTNUC( 8,NATOM) = EFTNUC( 8,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D15 * XCOOR * ZCOOR / R7
     &                        + D105 * (YCOOR**2) * XCOOR * ZCOOR / R9 )
C
C                          V_xyzz
C
                           EFTNUC( 9,NATOM) = EFTNUC( 9,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D15 * XCOOR * YCOOR / R7
     &                        + D105 * (ZCOOR**2) * XCOOR * YCOOR / R9 )
C
C                          V_xzzz
C
                           EFTNUC(10,NATOM) = EFTNUC(10,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D45 * XCOOR*ZCOOR / R7
     &                        + D105 * (ZCOOR**3) * XCOOR / R9 )
C
C                          V_yyyy
C
                           EFTNUC(11,NATOM) = EFTNUC(11,NATOM)
     &                        + CHARGE(I2) * (
     &                            D9 / R5
     &                        -  D90 * (YCOOR**2) / R7
     &                        + D105 * (YCOOR**4) / R9 )
C
C                          V_yyyz
C
                           EFTNUC(12,NATOM) = EFTNUC(12,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D45 * YCOOR*ZCOOR / R7
     &                        + D105 * (YCOOR**3) * ZCOOR / R9 )
C
C                          V_yyzz
C
                           EFTNUC(13,NATOM) = EFTNUC(13,NATOM)
     &                        + CHARGE(I2) * (
     &                            D3 / R5
     &                        -  D15 * YCOOR**2 / R7
     &                        -  D15 * ZCOOR**2 / R7
     &                        + D105 * (YCOOR**2) * (ZCOOR**2) / R9 )
C
C                          V_yzzz
C
                           EFTNUC(14,NATOM) = EFTNUC(14,NATOM)
     &                        + CHARGE(I2) * (
     &                        -  D45 * YCOOR*ZCOOR / R7
     &                        + D105 * (ZCOOR**3) * YCOOR / R9 )
C
C                          V_zzzz
C
                           EFTNUC(15,NATOM) = EFTNUC(15,NATOM)
     &                        + CHARGE(I2) * (
     &                            D9 / R5
     &                        -  D90 * (ZCOOR**2) / R7
     &                        + D105 * (ZCOOR**4) / R9 )
                        END IF
                        LL = LL + 3
                     END IF
 110              CONTINUE
 100           CONTINUE
               NATOM = NATOM + 1
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRP_PVC */
      SUBROUTINE PRP_PVC(EXPVAL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. Thyssen - Oct 12 1999
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: IPRPRP
      use nuclei_mod, only: NUCDEP, NUCIND
#include "implicit.h"
C
#include "dcbprp.h"
!#include "dcbexp.h"
!#include "dcbxlr.h"
      DIMENSION EXPVAL(*),WORK(*)
      logical, allocatable :: DOATOM(:)
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      allocate(DOATOM(NUCIND))
      CALL MEMGET('REAL',KTEMP,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KMPV,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQV,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIS,NUCDEP,WORK,KFREE,LFREE)
C
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)
      CALL PRP_PVC1(EXPVAL,DOATOM,WORK(KMPV),WORK(KQV),
     &              WORK(KIS),WORK(KTEMP),IPRPRP)
      deallocate(DOATOM)
      CALL MEMREL('PRP_PVC',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRP_PVC1(EXPVAL,DOATOM,DMPV,QW,IS,TEMP,IPRINT)
C***********************************************************************
C
C     Calculate parity violation - chirality
C
C     Written by J. Thyssen - Oct 12 1999
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use consts_mod, only: D0, D1, D2, D4
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbexp_mod, only: IPPVC
      use nuclei_mod, only: ISOTOP, IZATOM, NAMDEP, NUCDEG, NUCDEP,
     &  NUCIND
!      use codata_mod, only: GFERMI, S2THETAW
#include "implicit.h"
#include "dcbexp.h"
      DIMENSION EXPVAL(*),TEMP(*),DMPV(*),QW(*),
     &          IS(*)
      LOGICAL   DOATOM(*)
!#include "nqcc.h"
#include "codata.h"
C
      CALL QENTER('PRP_PVC')
C
C     Gather expectation values into TEMP
C     (only over symmetry distinct centers)
C
      DFAC = GFERMI/(D2*SQRT(D2))
C
      CALL DZERO(TEMP,NUCDEP)
      ITYP = 0
      DO 100 IATOM = 1,NUCIND
      IF(DOATOM(IATOM)) THEN
        DO IDEG = 1,NUCDEG(IATOM)
          ITYP = ITYP + 1
          TEMP(IATOM) = TEMP(IATOM) +
     &         EXPVAL(IPPVC(ITYP)) / NUCDEG(IATOM)
        END DO
      ENDIF
 100  CONTINUE
C
C     For each atom calculate exp. value
C
      DTOT = D0
      NATOM = 0
      DO 110 IATOM = 1, NUCIND
      IF(DOATOM(IATOM)) THEN
        DO IDEG = 1,NUCDEG(IATOM)
          NATOM =  NATOM + 1
          DMPV(NATOM) = TEMP(IATOM)
          ICHARGE = IZATOM(IATOM)
          IMASS = NINT(DISOTP(ICHARGE,ISOTOP(IATOM),'MASS'))
          NPROTON  = ICHARGE
          NNEUTRON = IMASS - ICHARGE
          IS(NATOM) = IMASS
          QW(NATOM) = (D1 - D4*S2THETAW) * NPROTON -
     &            D1 * NNEUTRON
           DTOT = DTOT + DMPV(NATOM)*QW(NATOM)*DFAC
        END DO
      ENDIF
 110  CONTINUE
C
C     Output section
C
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(A/)')
     &     '* Parity violation - chirality E(PV):'
      CALL PRSYMB(LUPRI,'-',76,0)
C
      WRITE(LUPRI,9000)
     &   'Nucleus','Isotope','integral    ','factor  ','contrib.    '
      WRITE(LUPRI,9000) ' ',' ','MPV(n)/au    ','QW(n)    '
      NATOM  = 0
      IATOMD = 0
      DO 120 IATOM = 1, NUCIND
      IF(DOATOM(IATOM)) THEN
        DO IDEG = 1,NUCDEG(IATOM)
          IATOMD = IATOMD + 1
          NATOM =  NATOM + 1
          WRITE(LUPRI,9010)
     &            NAMDEP(IATOMD),IS(NATOM),DMPV(NATOM),QW(NATOM),
     &            DMPV(NATOM)*QW(NATOM)*DFAC
        END DO
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
 120  CONTINUE
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,9020) 'Total EPV = ',DTOT
 9000 FORMAT(4X,A8,'  ',A7,' ',A20,'  ',A10,'  ',A20)
 9010 FORMAT(4X,A8,'   ',I4,'   ',1P,E20.12,'  ',0P,F11.5,'  ',
     &       1P,E20.12)
 9020 FORMAT(4X,A12,41X,1P,E20.12)
C
C
C
      CALL QEXIT('PRP_PVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_SPNSPN */
      SUBROUTINE PRP_SPNSPN(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR indirect spin-spin coupling
C     This routine is modelled after corresponding
C     non-relativistic routine , so thanks to Kenneth Ruud !
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbxlr_mod, only: NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEP
      use nuclei_mod, only: NUCDEP, NUCIND
#include "implicit.h"
C
!#include "dcbprp.h"
!#include "dcbexp.h"
#include "dcbxlr.h"
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &          ABCNV(NBFREQ,NLRBPT),WORK(*)
      logical, allocatable :: DOATOM(:)
      KFRSAV = KFREE
      IF(XLR_SKIPEP) THEN
        IC = 2
      ELSE
        IC = 1
      ENDIF
C
C     Allocate memory
C
      NSCOOR  = NUCDEP*3
      N2SCOR  = NSCOOR*NSCOOR
      NNDEP   = NUCDEP*(NUCDEP+1)/2
      allocate(DOATOM(NUCIND))
      CALL MEMGET('INTE',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,IC*N2SCOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHES,IC*N2SCOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSPAR,5*NNDEP  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NSCOOR*2 ,WORK,KFREE,LFREE)
C
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
      CALL PRP_SPNSP1(EXPVAL,ABLRF,ABCNV,
     &                DOATOM,WORK(KIATM),WORK(KRATM),
     &                WORK(KSMAT),WORK(KCHES),WORK(KSPAR),
     &                WORK(KIBUF),WORK,KFREE,LFREE)
      deallocate(DOATOM)
      CALL MEMREL('PRP_SPNSPN',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_SPNSP1 */
      SUBROUTINE PRP_SPNSP1(EXPVAL,ABLRF,ABCNV,DOATOM,IATINF,RATINF,
     &                      SPNMAT,CHESS,SPNPAR,IBUF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR indirect spin-spin coupling
C     This routine is modelled after corresponding
C     non-relativistic routine , so thanks to Kenneth Ruud !
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: ABUND, IPRPRP, IPSPNSPN
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, TRIAB
      use nuclei_mod, only: IZATOM, NAMDEP, NAMN, NUCDEG, NUCDEP, NUCIND
!      use dgroup_mod, only: IPQTOQ
!      use codata_mod, only: ALPHA2, NMAGN, NMAGNAU, PLANCK, XTHZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, J, K
      REAL(dp), PARAMETER :: D3 = 3.0D0, D0 = 0.0D0
      DOUBLE PRECISION AUTOHZ,REDHZ
C
#include "dcbprp.h"
!#include "dcbxpr.h"
!#include "dcbprl.h"
!#include "dcbexp.h"
#include "dcbxlr.h"
#include "dgroup.h"
#include "codata.h"
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &          ABCNV(NBFREQ,NLRBPT),
     &          IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          SPNMAT(3*NUCDEP,3*NUCDEP),
     &          CHESS(3*NUCDEP,3*NUCDEP),
     &          SPNPAR(5,*),IBUF(3*NUCDEP,2),WORK(*)
C
      CALL QENTER('PRP_SPNSP1')
      CALL TITLER
     &    ('INDIRECT NUCLEAR SPIN-SPIN-COUPLINGS','*',116)
C
      AUTOHZ = ALPHA2*ALPHA2*nmagnau*nmagnau*xthz
      REDHZ=1E19*nmagn*nmagn/planck
C
      KFRSAV = KFREE
      NSCOOR = 3*NUCDEP
      NCHES  = NSCOOR*NSCOOR
C
C     ***************************************
C     ***** Extract nuclear information *****
C    ****************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &  '@ name  charge  isotope   mass         abundance     g factor'
      CALL PRSYMB(LUPRI,'-',61,0)
      NSPEC = 0
      IDEP  = 0
      DO INUC = 1,NUCIND
      IDEG = NUCDEG(INUC)
      IF(DOATOM(INUC)) THEN
        TEST = .FALSE.
        IZ   = IZATOM(INUC)
        NISO = 0
        DO ISO = 1,5
          GVAL = DISOTP(IZ,ISO,'GVAL')
          IF(GVAL.NE.0) THEN
            ABND = DISOTP(IZ,ISO,'ABUNDANCE')
            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = GVAL
            ENDIF
          ENDIF
        ENDDO
        IF(NISO.GT.0) THEN
          IATINF(1,NSPEC+1) = IDEP + 1
          IATINF(2,NSPEC+1) = IZ
          IATINF(3,NSPEC+1) = NISO
          WRITE(LUPRI,
     &      '(A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &       '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
          DO J = 2,NISO
            WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &         '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
          ENDDO
          DO K = 2,IDEG
            IATINF(1,NSPEC+K) = IDEP + K
            IATINF(2,NSPEC+K) = IZ
            IATINF(3,NSPEC+K) = NISO
            DO L = 1,NISO
              RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
              RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
              RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
            ENDDO
          ENDDO
          NSPEC = NSPEC + IDEG
        ENDIF
      ENDIF
      IDEP = IDEP + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        CALL QEXIT('PRP_SPNSP1')
        RETURN
      ENDIF
C
C     *** Extract spin-spin couplings ***
C
      NA = 0
      NB = 0
      DO I = 1,NSCOOR
        IF(IPSPNSPN(1,I).GT.0) THEN
          NA = NA + 1
          IBUF(NA,1) = I
        ENDIF
        IF(IPSPNSPN(2,I).GT.0) THEN
          NB = NB + 1
          IBUF(NB,2) = I
        ENDIF
      ENDDO
C
C     Fill spin-spin matrices
C
#if defined (XH4)
      WRITE(LUPRI,'(A)') '*** WARNING *** Using code for XH4 !'
#endif
      DO IFREQ = 1,NBFREQ
C
        IF(BFREQ(IFREQ).EQ.D0) THEN
C          WRITE(LUPRI,'(A)') '*** Static values:'
        ELSE
          WRITE(LUPRI,'(A,F15.8,2X,A)')
     &    '*** Frequency    :',BFREQ(IFREQ),'a.u.'
        ENDIF
C
C       Fill spin-matrix
C
        CALL DZERO(SPNMAT,NCHES)
        DO IB = 1,NB
          JB = IBUF(IB,2)
          KB = IPSPNSPN(2,JB)
          DO IA = 1,NA
            JA = IBUF(IA,1)
            KA = IPSPNSPN(1,JA)
            SPNMAT(JA,JB) = ABLRF(IFREQ,KA,KB)
            IF(TRIAB) SPNMAT(JB,JA) = ABLRF(IFREQ,KA,KB)
          ENDDO
        ENDDO
#if defined (XH4)
      SPNMAT( 6, 4)   = SPNMAT(4,6)
      SPNMAT( 5, 4)   = SPNMAT(4,5)
      SPNMAT( 2, 2)   = SPNMAT(1,1)
      SPNMAT( 3, 3)   = SPNMAT(1,1)
      SPNMAT(14,12)   = SPNMAT(4,5)
      SPNMAT(12,14)   = SPNMAT(4,5)
      SPNMAT(11, 8)   = SPNMAT(4,5)
      SPNMAT( 8,11)   = SPNMAT(4,5)
      SPNMAT(14,13)   = SPNMAT(6,5)
      SPNMAT(13,14)   = SPNMAT(6,5)
      SPNMAT(11,10)   = SPNMAT(6,5)
      SPNMAT(10,11)   = SPNMAT(6,5)
      SPNMAT(11, 9)   = SPNMAT(6,5)
      SPNMAT( 9,11)   = SPNMAT(6,5)
      SPNMAT(15,14)   = SPNMAT(6,5)
      SPNMAT(14,15)   = SPNMAT(6,5)
      SPNMAT( 7, 4)   = SPNMAT(4,6)
      SPNMAT( 4, 7)   = SPNMAT(4,6)
      SPNMAT(13,12)   = SPNMAT(4,6)
      SPNMAT(12,13)   = SPNMAT(4,6)
      SPNMAT( 9, 8)   = SPNMAT(4,6)
      SPNMAT( 8, 9)   = SPNMAT(4,6)
      SPNMAT(15,12)   = SPNMAT(4,6)
      SPNMAT(12,15)   = SPNMAT(4,6)
      SPNMAT(10, 8)   = SPNMAT(4,6)
      SPNMAT( 8,10)   = SPNMAT(4,6)
      SPNMAT( 7, 7)   = SPNMAT(6,6)
      SPNMAT( 9, 9)   = SPNMAT(6,6)
      SPNMAT(10,10)   = SPNMAT(6,6)
      SPNMAT(15,15)   = SPNMAT(6,6)
      SPNMAT(13,13)   = SPNMAT(6,6)
      SPNMAT(15,13)   = SPNMAT(7,6)
      SPNMAT(13,15)   = SPNMAT(7,6)
      SPNMAT(10, 9)   = SPNMAT(7,6)
      SPNMAT( 9,10)   = SPNMAT(7,6)
      SPNMAT(14,14)   = SPNMAT(5,5)
      SPNMAT(11,11)   = SPNMAT(5,5)
#endif
        IF (IPRPRP.GE.4) THEN
          CALL HEADER('Raw coupling constants',1)
          CALL PRIH2S (SPNMAT,'SPNSPN',WORK,KFREE,LFREE)
          CALL PRQMAT(SPNMAT,NSCOOR,NSCOOR,
     &                NSCOOR,NSCOOR,1,IPQTOQ(1,0),LUPRI)
        ENDIF
C
C       Transform to non-symmetry basis
C
        CALL DZERO(CHESS,NCHES)
Cjth - changed because of the new version of TRACOR ( => new version of TRAHES)
        CALL TRAHES2(SPNMAT,NSCOOR,CHESS,NSCOOR,NSCOOR,2,
     &               WORK(KFREE),LFREE)
        IF (IPRPRP.GE.4) THEN
          CALL HEADER('Nsym coupling constants',1)
          CALL PRIH2S (CHESS,'SPNSPN',WORK,KFREE,LFREE)
          CALL PRQMAT(CHESS,NSCOOR,NSCOOR,
     &                NSCOOR,NSCOOR,1,IPQTOQ(1,0),LUPRI)
        ENDIF
Cold code        CALL TRAHES(SPNMAT,NSCOOR,CHESS,NSCOOR,NSCOOR,2)
C

C       Print indirect spin-spin coupling tensors (Hz)
C       for each pair of nuclei and isotope combination.
        IF (IPRPRP .GE. 2) THEN
          CALL HEADER
     &    ('Indirect Spin-spin tensors (Hz):',1)
C         Loop over all unique pairs of nuclei
          DO I = 1,NSPEC
            IISO   = IATINF(3,I)
            IATOM1 = IATINF(1,I)
            DO J = 1,(I-1)
              JISO   = IATINF(3,J)
              IATOM2 = IATINF(1,J)
C
C             Packed index for storing pairs
              IATIJ = IATOM1*(IATOM1 -1)/2 + IATOM2
C
C             Locate the 3x3 spin-spin tensor block corresponding
C             to the atom pair (IATOM1,IATOM2) for the global (NCHES x NCHES) matrix CHESS.
              IA = 3*(IATOM1-1)
              JA = 3*(IATOM2-1)
C
C             Loop over all isotope combinations for the curent atom pair.
              DO I3 = 1,IISO
                IF(IATOM1.EQ.IATOM2) THEN
                  J3END = I3 - 1
                ELSE
                  J3END = JISO
                ENDIF
C
                DO J3 = 1,J3END
C                 Atomic units to Hz using the gyromagnetic ratios
C                 of the selected isotopes.
                  FACTOR = AUTOHZ * RATINF(3,I3,I) * RATINF(3,J3,J)
C
C                 Print atom/isotope information.
                  WRITE(LUPRI,'(/,3X,A5,3X,A3,2X,A5,3X,A3)')
     &            'Atom1','iso','Atom2','iso'
C
                  WRITE(LUPRI,'(3X,A6,2X,A3,2X,A6,2X,A3)')
     &            '------','---','------','---'
C
                  WRITE(LUPRI,'(3X,A6,2X,I3,2X,A6,2X,I3)')
     &            NAMDEP(IATOM1),I3,NAMDEP(IATOM2),J3
C
C                 Print the 3x3 indirect spin-spin coupling tensor
C                 in Cartesian coordinates (Hz).
                  WRITE (LUPRI,'(/,18X,3(A,11X))') 'x', 'y', 'z'
                  WRITE(LUPRI,'(3(3X,A,6X,3F12.4,/))')
     &            'x',(FACTOR*CHESS(IA+1,JA+JC),JC=1,3),
     &            'y',(FACTOR*CHESS(IA+2,JA+JC),JC=1,3),
     &            'z',(FACTOR*CHESS(IA+3,JA+JC),JC=1,3)
C
                END DO
              END DO
C
C             Separator between different atom pairs.
              IF (I .GT. 1) CALL PRSYMB(LUPRI,'=',46,0)
C
            END DO
          END DO
        END IF
C
C
C       Analyze results
C
        CALL AROUND('Final spin-spin-couplings J (Hz)')
CTROND  Insert dimensions !!
        WRITE(LUPRI,'(2A15,5(2X,A12))')
     &  'At1  #  iso  : ','At2  #  iso  : ',
     &  'Isotropic   ','Anisotropic ','Asymmetry   ',
     &  'S parameter ','A parameter '
        CALL PRSYMB(LUPRI,'-',100,0)
        DO I = 1,NSPEC
          IISO   = IATINF(3,I)
          IATOM1 = IATINF(1,I)
          DO J = 1,(I-1)
            JISO   = IATINF(3,J)
            IATOM2 = IATINF(1,J)
            IATIJ = IATOM1*(IATOM1 -1)/2 + IATOM2
            CALL SPIANA(CHESS,SPNPAR(1,IATIJ),IATOM1,IATOM2,
     &                      WORK,KFREE,LFREE)
            DO I3 = 1,IISO
              IF(IATOM1.EQ.IATOM2) THEN
                J3END = I3 - 1
              ELSE
                J3END = JISO
              ENDIF
              DO J3 = 1,J3END
                FACTOR = AUTOHZ*RATINF(3,I3,I)*RATINF(3,J3,J)
                WRITE(LUPRI,'(2(A6,I5,2X,A1,1X),5(2X,F12.4))')
     &                NAMDEP(IATOM1),I3,':',NAMDEP(IATOM2),J3,':',
     &                FACTOR*SPNPAR(1,IATIJ),
     &                FACTOR*SPNPAR(2,IATIJ),
     &                       SPNPAR(3,IATIJ),
     &                FACTOR*SPNPAR(4,IATIJ),
     &                FACTOR*SPNPAR(5,IATIJ)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL PRSYMB(LUPRI,'-',100,0)
!radovan: print out reduced SSCC in SI units
        FACTOR = AUTOHZ/REDHZ
        CALL AROUND('Final reduced spin-spin-couplings K '//
     &              '(SI:10**19 m-2 kg s-2 A-2)')
        WRITE(LUPRI,'(2A15,5(2X,A12))')
     &  'At1  #  iso  : ','At2  #  iso  : ',
     &  'Isotropic   ','Anisotropic ','Asymmetry   ',
     &  'S parameter ','A parameter '
        CALL PRSYMB(LUPRI,'-',100,0)
        DO I = 1,NSPEC
          IISO   = IATINF(3,I)
          IATOM1 = IATINF(1,I)
          DO J = 1,(I-1)
            JISO   = IATINF(3,J)
            IATOM2 = IATINF(1,J)
            IATIJ = IATOM1*(IATOM1 -1)/2 + IATOM2
            CALL SPIANA(CHESS,SPNPAR(1,IATIJ),IATOM1,IATOM2,
     &                      WORK,KFREE,LFREE)
            DO I3 = 1,IISO
              IF(IATOM1.EQ.IATOM2) THEN
                J3END = I3 - 1
              ELSE
                J3END = JISO
              ENDIF
              DO J3 = 1,J3END
                WRITE(LUPRI,'(2(A6,I5,2X,A1,1X),5(2X,F12.4))')
     &                NAMDEP(IATOM1),I3,':',NAMDEP(IATOM2),J3,':',
     &                FACTOR*SPNPAR(1,IATIJ),
     &                FACTOR*SPNPAR(2,IATIJ),
     &                       SPNPAR(3,IATIJ),
     &                FACTOR*SPNPAR(4,IATIJ),
     &                FACTOR*SPNPAR(5,IATIJ)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL PRSYMB(LUPRI,'-',100,0)

      ENDDO
      CALL QEXIT('PRP_SPNSP1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spiana */
      SUBROUTINE SPIANA(PRIBAS,SPNPAR,IATOM1,IATOM2,
     &                  WORK,KFREE,LFREE)
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCDEP
#include "implicit.h"
      INTEGER :: I, J
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D3 = 3.0D0, D1 = 1.0D0)
      DIMENSION PRIBAS(3*NUCDEP,3*NUCDEP),SPNPAR(5),COUPLM(3,3),
     &          PVAL(3), AXES(3,3), SYM(3,3), ANTI(3,3),WORK(*)
C
      CALL DZERO(SPNPAR,5)
C
C     Extract coupling information for the atomic centers
C
      IOFF1 = 3*(IATOM1 - 1)
      IOFF2 = 3*(IATOM2 - 1)
      DO ICOOR2 = 1, 3
        DO ICOOR1 = 1, 3
          COUPLM(ICOOR1,ICOOR2) = PRIBAS(IOFF1+ICOOR1,IOFF2+ICOOR2)
        ENDDO
      ENDDO
C
C     Extract symmetric and anti-symmetric parts
C
      SPNPAR(1) = (COUPLM(1,1) + COUPLM(2,2) + COUPLM(3,3))/D3
      DO I = 1, 3
        DO J = 1, 3
          ANTI(I,J) = (COUPLM(I,J) - COUPLM(J,I))/D2
          SYM (I,J) = (COUPLM(I,J) + COUPLM(J,I))/D2
        ENDDO
        SYM(I,I) = SYM(I,I) - SPNPAR(1)
      ENDDO
C
C     Diagonalize symmetric part to determine principal axis system
C     IZ(1) is largest absolute eigenvalue, IY(3) smallest, IX
C
      CALL RSJACO(3,3,3,SYM,PVAL,1,-2,0,AXES)
C
      SPNPAR(2) = PVAL(1) - (PVAL(2) + PVAL(3))/D2
      IF (PVAL(1) .EQ. 0) THEN
        SPNPAR(3) = D0
      ELSE
        SPNPAR(3) = (PVAL(2) - PVAL(3))/PVAL(1)
      END IF
      SPNPAR(5) = SQRT(DDOT(9,ANTI,1,ANTI,1)/D2)
      SPNPAR(4) = ABS(SPNPAR(2))*SQRT(D1+(SPNPAR(3)*SPNPAR(3))/D3)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prih2s */
      SUBROUTINE PRIH2S(SHESS,KEY,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Adapted from Nelson's OUTPAK for Hessians March 1985 tuh
C     Revised 16-Dec-1983 by Hans Jorgen Aa. Jensen.
C     Jun 28 1988 tuh - modified for symmetry
C     Less memory intensive version Feb 1997 - tsaue
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NATOMS, NUCDEP
#include "implicit.h"
C
      CHARACTER*(*) KEY
      DIMENSION SHESS(*),WORK(*)
C
      KFRSAV = KFREE
      IF (KEY .EQ. 'CENTERS') THEN
         NROW = 3*NUCDEP
         ITYPE = 1
      ELSE IF (KEY .EQ. 'ATOMS') THEN
         NROW = 3*NATOMS
         ITYPE = 1
      ELSE IF (KEY .EQ. 'SPNSPN') THEN
         NROW = 3*NUCDEP
         ITYPE = 2
      ELSE
         WRITE (LUPRI,'(//A/2A//)') ' >>> ERROR in PRIH2S ',
     &      '     - incorrect specification of keyword:', KEY
         CALL QUIT('Incorrect keyword in PRIH2S.')
      END IF
      NDIM = NROW*NROW
      CALL MEMGET('REAL',KBUF,NDIM,WORK,KFREE,LFREE)
      CALL PRIH21(SHESS,WORK(KBUF),NROW,ITYPE,WORK(KFREE),LFREE)
      CALL MEMREL('PRIH2S',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prih21 */
      SUBROUTINE PRIH21(SHESS,CHESS,NROW,ITYPE,WORK,LWORK)
C***********************************************************************
C
C     Adapted from Nelson's OUTPAK for Hessians March 1985 tuh
C     Revised 16-Dec-1983 by Hans Jorgen Aa. Jensen.
C     Jun 28 1988 tuh - modified for symmetry
C     Less memory intensive version Feb 1997 - tsaue
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: NAMEX
      use symmet_mod, only: IPTCOR, MAXREP, NCRREP
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, J, K
      INTEGER, PARAMETER :: KCOL=6
      REAL(dp), PARAMETER :: D0=0.D00
      DIMENSION SHESS(NROW,NROW), CHESS(NROW,NROW)
      DIMENSION WORK(LWORK)
C
      INTEGER BEGIN
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
      CALL QENTER('PRIH21')
C
      IF (MAXREP .EQ. 0) THEN
         LAST = MIN(NROW,KCOL)
         BEGIN= 1
   50    CONTINUE
         WRITE (LUPRI,1000) (NAMEX(I),I = BEGIN,LAST)
         WRITE (LUPRI,'()')
         NCOL = 1
         DO 100 K = BEGIN,NROW
            DO 200 I = 1,NCOL
               IF (SHESS(K,(BEGIN-1)+I) .NE. D0) GO TO 400
  200       CONTINUE
            GO TO 300
  400       WRITE (LUPRI,2000) NAMEX(K),
     &            (SHESS(K,(BEGIN-1)+J),J=1,NCOL)
            IF (MOD(K,3) .EQ. 0) WRITE (LUPRI,'()')
  300       IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
  100    CONTINUE
         WRITE (LUPRI,'()')
         LAST = MIN(LAST+KCOL,NROW)
         BEGIN= BEGIN+NCOL
         IF (BEGIN.LE.NROW) GO TO 50
      ELSE
         NOFF = 0
         DO 500 IREP = 0, MAXREP
C            IF (DOSYM(IREP+1) .AND. (NCRREP(IREP,ITYPE) .GT. 0)) THEN
            IF ((NCRREP(IREP,ITYPE) .GT. 0)) THEN
               WRITE (LUPRI,'(/11X,A,I1/)') 'Symmetry ',IREP + 1
               NROW  = NOFF + NCRREP(IREP,ITYPE)
               LAST  = MIN(NROW,NOFF + KCOL)
               BEGIN = NOFF + 1
C
  550          CONTINUE
                  WRITE(LUPRI,1000)(NAMEX(IPTCOR(I,ITYPE)),I=BEGIN,LAST)
                  WRITE(LUPRI,'()')
                  NCOL = 1
                  DO 40 K = BEGIN, NROW
                     WRITE (LUPRI,2000) NAMEX(IPTCOR(K,ITYPE)),
     &                     (SHESS(K,(BEGIN-1)+J),J=1,NCOL)
                     IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40             CONTINUE
                  WRITE (LUPRI,'()')
                  LAST  = MIN(LAST + KCOL,NROW)
                  BEGIN = BEGIN + NCOL
               IF (BEGIN.LE.NROW) GO TO 550
            END IF
            NOFF = NOFF + NCRREP(IREP,ITYPE)
  500    CONTINUE
C
C        Print in non-symmetry basis
C
         WRITE (LUPRI,'(//)')
         CALL DZERO(CHESS,NROW*NROW)
         CALL TRAHES2(SHESS,NROW,CHESS,NROW,NROW,ITYPE,
     &                WORK(KFREE),LFREE)
Cold code         CALL TRAHES(SHESS,NROW,CHESS,NROW,NROW,ITYPE)
         CALL PR2DER(CHESS,NROW,NROW,LUPRI)
      END IF
      WRITE (LUPRI,'()')
C
      CALL QEXIT('PRIH21')
C
      CALL FLSHFO(LUPRI)
      RETURN
 1000 FORMAT (8X,6(3X,A6,3X),(3X,A6,3X))
 2000 FORMAT (1X,A6,6F12.6)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_dso */
      SUBROUTINE DEF_DSO(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define diamagnetic spin-orbit operator
C     Written by T.Saue - Feb 21 1997
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     DSOINT: INTTYP = 12
C     (-12: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-12)
      CALL DEF_DSO1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_DSO1 */
      SUBROUTINE DEF_DSO1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define diamagnetic spin-orbit operator
C     Written by T.Saue - Feb 21 1997
C
C***********************************************************************
      use cbiher_mod, only: NPQUAD, TRIANG
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NAMN, NUCIND
      use symmet_mod, only: IPTCNT, MAXREP
!      use dcbprp_mod, only: IPDSO, NDSOIN, NPQUAD, TRIANG
!      use dcbexp_mod, only: LEXPP, MAXEXP, NEXPP
      use chrnos_mod, only: CHRNOS
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: DM1 = -1.0D0
C
!#include "dgroup.h"
!#include "pgroup.h"
!#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbexp.h"
      LOGICAL DOATOM(NUCIND),SAME
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3),IOP(3)
C
      WRITE(LUPRI,'(3X,A,A/A,I5)')
     &  'Diamagnetic spin-orbit integrals are generated by ',
     &  'Gaussian quadrature.',
     &  'Number of points: ',NPQUAD
      NDSOIN  = 0
      IPTYP   = 1
      NPCOMP  = 1
      PFAC(1) = DM1
C
C     Irreps
C
      DO 100 IREPO = 0, MAXREP
C
C        Atoms
C
         DO 200 IATOM1 = 1, NUCIND
         IF (DOATOM(IATOM1)) THEN
            MXATM2 = NUCIND
            IF (TRIANG) MXATM2 = IATOM1
            DO 400 IATOM2 = 1, MXATM2
            IF (DOATOM(IATOM2)) THEN
               SAME = TRIANG .AND. IATOM1.EQ.IATOM2
C
C              Cartesian directions
C
               DO 500 ICOOR1 = 1, 3
                  ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
                  IF (ISCOR1 .GT. 0) THEN
                     MXCR2 = 3
                     IF (SAME) MXCR2 = ICOOR1
                     DO 600 ICOOR2 = 1, MXCR2
                        ISCOR2 = IPTCNT(3*(IATOM2-1)+ICOOR2,IREPO,2)
                        IF (ISCOR2 .GT. 0) THEN
                           NDSOIN = NDSOIN + 1
                           PNAME = 'DS'
     &                       //NAMN(IATOM1)
     &                       //CHRNOS(ISCOR1/100)
     &                       //CHRNOS(MOD(ISCOR1,100)/10)
     &                       //CHRNOS(MOD(MOD(ISCOR1,100),10))
     &                       //NAMN(IATOM2)
     &                       //CHRNOS(ISCOR2/100)
     &                       //CHRNOS(MOD(ISCOR2,100)/10)
     &                       //CHRNOS(MOD(MOD(ISCOR2,100),10))
                           PLABEL(1) = 'DS'
     &                       //CHRNOS(ISCOR1/100)
     &                       //CHRNOS(MOD(ISCOR1,100)/10)
     &                       //CHRNOS(MOD(MOD(ISCOR1,100),10))
     &                       //CHRNOS(ISCOR2/100)
     &                       //CHRNOS(MOD(ISCOR2,100)/10)
     &                       //CHRNOS(MOD(MOD(ISCOR2,100),10))
                           CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                         INDXPR,ISYXPR,ITRXPR,IPRINT)
                           IF(IPRINT.GE.1) THEN
                             CALL WRIXPR(NDSOIN,INDXPR)
                           ENDIF
                           CALL OP1IND('XDSO',IPIND,LEXPP,NEXPP,
     &                                 INDXPR,MAXEXP)
                           IPDSO(NDSOIN) = IPIND
                        END IF
  600                CONTINUE
                  END IF
  500          CONTINUE
            END IF
  400       CONTINUE
         END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_DSO */
      SUBROUTINE PRP_DSO(EXPVAL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR indirect spin-spin coupling
C     This routine is modelled after corresponding
C     non-relativistic routine , so thanks to Kenneth Ruud !
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCDEP, NUCIND
#include "implicit.h"
C
!#include "dcbprp.h"
!#include "dcbexp.h"
!#include "dcbxlr.h"
      DIMENSION EXPVAL(*),WORK(*)
      logical, allocatable :: DOATOM(:)
      KFRSAV = KFREE
      NSCOOR  = NUCDEP*3
      N2SCOR  = NSCOOR*NSCOOR
      NNDEP   = NUCDEP*(NUCDEP+1)/2
      allocate(DOATOM(NUCIND))
      CALL MEMGET('INTE',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,N2SCOR   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHES,N2SCOR   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSPAR,5*NNDEP  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NSCOOR*2 ,WORK,KFREE,LFREE)
C
C     (-12: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-12)
      CALL PRP_DSO1(EXPVAL,
     &              DOATOM,WORK(KIATM),WORK(KRATM),
     &              WORK(KSMAT),WORK(KCHES),WORK(KSPAR),
     &              WORK(KIBUF),WORK,KFREE,LFREE)
      deallocate(DOATOM)
      CALL MEMREL('PRP_DSO',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_DSO1 */
      SUBROUTINE PRP_DSO1(EXPVAL,DOATOM,IATINF,RATINF,
     &              SPNMAT,CHESS,SPNPAR,IBUF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR indirect spin-spin coupling
C     This routine is modelled after corresponding
C     non-relativistic routine , so thanks to Kenneth Ruud !
C
C     Written by T.Saue Feb 19 1997
C     Last revised Jan 19 1998 - tsaue
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: ABUND, IPDSO, IPRPRP, NDSOIN
!      use dcbxpr_mod, only: IPRPLBL
!      use dcbprl_mod, only: PRPLBL
!      use dcbexp_mod, only: LEXPP
      use nuclei_mod, only: IZATOM, NAMDEP, NAMN, NUCDEG, NUCDEP, NUCIND
!      use dgroup_mod, only: IPQTOQ
!      use codata_mod, only: ALPHA2, NMAGN, NMAGNAU, PLANCK, XTHZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, J, K
      REAL(dp), PARAMETER :: D3 = 3.0D0, D0 = 0.0D0
      DOUBLE PRECISION AUTOHZ,REDHZ
C
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "dgroup.h"
#include "codata.h"
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION EXPVAL(*),
     &          IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          SPNMAT(3*NUCDEP,3*NUCDEP),
     &          CHESS(3*NUCDEP,3*NUCDEP),
     &          SPNPAR(5,*),IBUF(3*NUCDEP,2),WORK(*)
C
C
      AUTOHZ = ALPHA2*ALPHA2*nmagnau*nmagnau*xthz
      REDHZ=1E19*nmagn*nmagn/planck
C
      KFRSAV = KFREE
      CALL QENTER('PRP_DSO1')
      CALL TITLER
     &    ('DSO CONTRIBUTIONS','*',116)
      NSCOOR = 3*NUCDEP
      NCHES  = NSCOOR*NSCOOR
C
C
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &  '@ name  charge  isotope   mass         abundance     g factor'
      CALL PRSYMB(LUPRI,'-',61,0)
      NSPEC = 0
      IDEP  = 0
      DO INUC = 1,NUCIND
      IDEG = NUCDEG(INUC)
      IF(DOATOM(INUC)) THEN
        TEST = .FALSE.
        IZ   = IZATOM(INUC)
        NISO = 0
        DO ISO = 1,5
          GVAL = DISOTP(IZ,ISO,'GVAL')
          IF(GVAL.NE.0) THEN
            ABND = DISOTP(IZ,ISO,'ABUNDANCE')
            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = GVAL
            ENDIF
          ENDIF
        ENDDO
        IF(NISO.GT.0) THEN
          IATINF(1,NSPEC+1) = IDEP + 1
          IATINF(2,NSPEC+1) = IZ
          IATINF(3,NSPEC+1) = NISO
          WRITE(LUPRI,
     &      '(A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &       '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
          DO J = 2,NISO
            WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &         '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
          ENDDO
          DO K = 2,IDEG
            IATINF(1,NSPEC+K) = IDEP + K
            IATINF(2,NSPEC+K) = IZ
            IATINF(3,NSPEC+K) = NISO
            DO L = 1,NISO
              RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
              RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
              RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
            ENDDO
          ENDDO
          NSPEC = NSPEC + IDEG
        ENDIF
      ENDIF
      IDEP = IDEP + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        CALL QEXIT('PRP_DSO1')
        RETURN
      ENDIF
C
C
C     Fill spin-spin matrices
C
      CALL DZERO(SPNMAT(1,1),NCHES)
      DO I = 1,NDSOIN
        IOP = IPDSO(I)
        IPL = IPRPLBL(1,LEXPP(IOP))
        READ(PRPLBL(IPL)(3:5),'(I3)') IA
        READ(PRPLBL(IPL)(6:8),'(I3)') IB
        SPNMAT(IA,IB) = EXPVAL(IOP)
        SPNMAT(IB,IA) = EXPVAL(IOP)
      ENDDO
      IF (IPRPRP.GE.4) THEN
        CALL HEADER('Raw DSO coupling constants',1)
        CALL PRIH2S (SPNMAT,'SPNSPN',WORK,KFREE,LFREE)
        CALL PRQMAT(SPNMAT,NSCOOR,NSCOOR,
     &              NSCOOR,NSCOOR,1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Transform to non-symmetry basis
C
      CALL DZERO(CHESS,NCHES)
Cjth - changed because of the new version of TRACOR
C      ( => new version of TRAHES)
      CALL TRAHES2(SPNMAT,NSCOOR,CHESS,NSCOOR,NSCOOR,2,
     &               WORK(KFREE),LFREE)
      IF (IPRPRP.GE.4) THEN
        CALL HEADER('Nsym coupling constants',1)
        CALL PRIH2S (CHESS,'SPNSPN',WORK,KFREE,LFREE)
        CALL PRQMAT(CHESS,NSCOOR,NSCOOR,
     &              NSCOOR,NSCOOR,1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C       Analyze results
C
      CALL AROUND('DSO contribution to SSCC J (Hz)')
CTROND  Insert dimensions !!
      WRITE(LUPRI,'(2A15,5(2X,A12))')
     &  'At1  #  iso  : ','At2  #  iso  : ',
     &  'Isotropic   ','Anisotropic ','Asymmetry   ',
     &  'S parameter ','A parameter '
      CALL PRSYMB(LUPRI,'-',100,0)
      DO I = 1,NSPEC
        IISO   = IATINF(3,I)
        IATOM1 = IATINF(1,I)
        DO J = 1,(I-1)
          JISO   = IATINF(3,J)
          IATOM2 = IATINF(1,J)
          IATIJ = IATOM1*(IATOM1 -1)/2 + IATOM2
          CALL SPIANA(CHESS,SPNPAR(1,IATIJ),IATOM1,IATOM2,
     &                    WORK,KFREE,LFREE)
          DO I3 = 1,IISO
            IF(IATOM1.EQ.IATOM2) THEN
              J3END = I3 - 1
            ELSE
              J3END = JISO
            ENDIF
            DO J3 = 1,J3END
              FACTOR = AUTOHZ*RATINF(3,I3,I)*RATINF(3,J3,J)
              WRITE(LUPRI,'(2(A6,I5,2X,A1,1X),5(2X,F12.4))')
     &              NAMDEP(IATOM1),I3,':',NAMDEP(IATOM2),J3,':',
     &              FACTOR*SPNPAR(1,IATIJ),
     &              FACTOR*SPNPAR(2,IATIJ),
     &                     SPNPAR(3,IATIJ),
     &              FACTOR*SPNPAR(4,IATIJ),
     &              FACTOR*SPNPAR(5,IATIJ)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL PRSYMB(LUPRI,'-',100,0)
!radovan: print out reduced SSCC in SI units
      CALL AROUND('DSO contribution to reduced SSCC K '//
     &     '(SI:10**19 m-2 kg s-2 A-2)')
      WRITE(LUPRI,'(2A15,5(2X,A12))')
     &  'At1  #  iso  : ','At2  #  iso  : ',
     &  'Isotropic   ','Anisotropic ','Asymmetry   ',
     &  'S parameter ','A parameter '
      CALL PRSYMB(LUPRI,'-',100,0)
      DO I = 1,NSPEC
        IISO   = IATINF(3,I)
        IATOM1 = IATINF(1,I)
        DO J = 1,(I-1)
          JISO   = IATINF(3,J)
          IATOM2 = IATINF(1,J)
          IATIJ = IATOM1*(IATOM1 -1)/2 + IATOM2
          CALL SPIANA(CHESS,SPNPAR(1,IATIJ),IATOM1,IATOM2,
     &                    WORK,KFREE,LFREE)
          DO I3 = 1,IISO
            IF(IATOM1.EQ.IATOM2) THEN
              J3END = I3 - 1
            ELSE
              J3END = JISO
            ENDIF
            DO J3 = 1,J3END
              FACTOR = AUTOHZ/REDHZ
              WRITE(LUPRI,'(2(A6,I5,2X,A1,1X),5(2X,F12.4))')
     &              NAMDEP(IATOM1),I3,':',NAMDEP(IATOM2),J3,':',
     &              FACTOR*SPNPAR(1,IATIJ),
     &              FACTOR*SPNPAR(2,IATIJ),
     &                     SPNPAR(3,IATIJ),
     &              FACTOR*SPNPAR(4,IATIJ),
     &              FACTOR*SPNPAR(5,IATIJ)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL PRSYMB(LUPRI,'-',100,0)

      CALL QEXIT('PRP_DSO1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck trahes2 */
      SUBROUTINE TRAHES2(SHESS,NDIMS,CHESS,NDIMC,NCOOR,ITYPE,
     &                   WORK,LWORK)
C***********************************************************************
C
C     Interface between the old version of TRAHES and the new version.
C     Written by J. Thyssen - 1997/09/18
C
C***********************************************************************
#include "implicit.h"
C
C
      DIMENSION WORK(LWORK)
      DIMENSION SHESS(NDIMS,NDIMS), CHESS(NDIMC,NDIMC)
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
      CALL QENTER('TRAHES2')
C
      CALL MEMGET('REAL',KCSTRA,NCOOR*NCOOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,NCOOR*NCOOR,WORK,KFREE,LFREE)
      CALL TRAHES(SHESS,NDIMS,CHESS,WORK(KCSTRA),WORK(KSCTRA),
     &            NDIMC,NCOOR,ITYPE)
      CALL MEMREL('TRAHES2',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('TRAHES2')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_SHIELD */
      SUBROUTINE DEF_SHIELD(IPRINT)
C***********************************************************************
C
C     Define operators for NMR shieldings or ESR
C     DEF_SPNSPN by T. Saue used as skeleton
C     Written by T. Enevoldsen - Oct. 1997
C     Last revision: 7. Nov. 1997
C                    22. Mar. 2000 (added ESR options, Anette Noerager + HJAaJ)
C
C   2002,2003/MI&HJAaJ - added TEC's changes for London atomic orbitals
C
C***********************************************************************
!      use codata_mod, only: G1FAC
!      use pgroup_mod, only: REP
!      use dcbgen_mod, only: CVAL
!      use dcbprp_mod, only: DOESR, DOEXP, ESRGTENS, ESR_HFCC, IPESR, IPRNST, IPSHIELD, LONDON, NNSTIN, SHIELD
!      use dcbesr_mod, only: IPSIGMA, LESRP, MAXESR, NESRP
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
!      use dcbexp_mod, only: LEXPP, MAXEXP, NEXPP
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use symmet_mod, only: IPTCNT, ISYMAX, MAXREP
      use consts_mod, only: D0
      use chrxyz_mod, only: CHRXYZ
      use chrnos_mod, only: CHRNOS
      use nuclei_mod, only: NAMN, NUCDEP, NUCIND
#include "implicit.h"
      INTEGER :: I
#include "codata.h"
!#include "dgroup.h"
#include "pgroup.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "dcbexp.h"
!#include "dcbnmr.h"

      logical, allocatable :: DOATOM(:)
      CHARACTER PNAME*16, PLABEL(3)*8,LABEL*4
      REAL*8    PFAC(3),IOP(3)
C
      CALL QENTER('DEF_SHIELD')
C
      allocate(DOATOM(NUCIND))
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)

      NNSTIN = 0
      NSCOOR = 3*NUCDEP
      CALL IZERO(IPSHIELD,3+MXCOOR)
C
C     Here we set up the operators needed for the linear reponse
C
C     Define left-hand operators: (-c) i.alpha x r/r^3
C     ================================================
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
C
CMI previously were here + signs (old TEC's convention)
CMI but the - sign is correct, so everything is rearranged
C
      PFAC(1)  = -CVAL
      PFAC(2)  = -CVAL
      NPCOMP   = 2
C
C     Irreps
C
      DO 100 IREPO = 0, MAXREP
C
C       Atoms
C
        DO 200 IATOM1 = 1, NUCIND
        IF(DOATOM(IATOM1)) THEN
          MXATM2 = NUCIND
C          IF (TRIAB) MXATM2 = IATOM1
C
C         Cartesian directions
C
          DO 300 ICOOR1 = 1, 3
            ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
            IF (ISCOR1 .GE. 1000) THEN
              CALL QUIT(
     &'DEF_SHIELD: Max 1000 coords in this implementation !')
            ENDIF
            IF (ISCOR1 .GT. 0) THEN
              IFIRST = ISCOR1/100
              ISECND = MOD(ISCOR1,100)/10
              ITHIRD = MOD(MOD(ISCOR1,100),10)
              PNAME    = 'NUCSP'//CHRNOS(IFIRST)
     &                   //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                   //NAMN(IATOM1)//REP(IREPO)//CHRXYZ(ICOOR1)
              IPTYP    = IOP(ICOOR1)
              IC1      = MOD(ICOOR1+1,3) + 1
              IC2      = MOD(ICOOR1  ,3) + 1
              IREP1    = IEOR(IREPO,ISYMAX(IC2,1))
              IREP2    = IEOR(IREPO,ISYMAX(IC1,1))
              IS1      = IPTCNT(3*(IATOM1 - 1) + IC1,IREP1,1)
              IS2      = IPTCNT(3*(IATOM1 - 1) + IC2,IREP2,1)
              IF((IS1.LE.0).OR.(IS2.LE.0))
     &          WRITE(LUPRI,'(A,5I5)')
     &          'DEF_SHIELD: sym. error !',IATOM1,ICOOR1,IC1,IC2
              IFIRST = IS1/100
              ISECND = MOD(IS1,100)/10
              ITHIRD = MOD(MOD(IS1,100),10)
              PLABEL(1)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              IFIRST = IS2/100
              ISECND = MOD(IS2,100)/10
              ITHIRD = MOD(MOD(IS2,100),10)
              PLABEL(2)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
#ifdef MOD_LAO_REARRANGED
              IF (SHIELD .OR. shielding_rearrange) THEN
#else
              IF (SHIELD) THEN
#endif
                CALL OP1IND('XSHIELDA',INDOP1,LLRAPU,NLRAPT,
     &                      INDXPR,MAXLLR)
                IPSHIELD(ISCOR1) = INDOP1
              END IF
              IF (ESR_HFCC) THEN
                CALL OP1IND('ESR_HFCC',INDOP1,LESRP,NESRP,INDXPR,MAXESR)
                IPESR(ISCOR1,1) = INDOP1
              END IF

C Expectation value contribution is always from the London phase factor
C (Note! it is not the diamagnetic term!)
              LABEL = 'RNST'
              IF (LONDON) THEN
C                ... make sure the PRPEXP is always called !!! /hjaaj nov2002
                 DOEXP = .TRUE.
                 DO ICOOR2 = 1, 3
                    IPTYP = IOP(ICOOR2)
                    NNSTIN = NNSTIN + 1
                    ICOORP1 = MOD(ICOOR2+1,3) + 1
                    ICOORP2 = MOD(ICOOR2,3) + 1
                    IFIRST = ISCOR1/100
                    ISECND = MOD(ISCOR1,100)/10
                    ITHIRD = MOD(MOD(ISCOR1,100),10)
                    PNAME = LABEL//NAMN(IATOM1)//CHRNOS(IFIRST)//
     &                    CHRNOS(ISECND)//CHRNOS(ITHIRD)//CHRXYZ(ICOOR2)
                    PLABEL(1) = CHRNOS(IFIRST)//CHRNOS(ISECND)//
     &                          CHRNOS(ITHIRD)//LABEL//
     &                          CHRXYZ(ICOORP1)
                    PLABEL(2) = CHRNOS(IFIRST)//CHRNOS(ISECND)//
     &                          CHRNOS(ITHIRD)//LABEL//
     &                          CHRXYZ(ICOORP2)
                    CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                          INDXPR,ISYXPR,ITRXPR,IPRINT)
                    IF (IPRINT.GE.5) THEN
                       !  print out assigned operators of wanted property
                       write(LUPRI,'(4X,A)')
     &                 "London atomic orbitals expectation value term:"
                       CALL WRIXPR(NNSTIN,INDXPR)
                    ENDIF
                    CALL OP1IND(LABEL,IPIND,LEXPP,NEXPP,INDXPR,MAXEXP)
                    IPRNST(NNSTIN) = IPIND
                 END DO
              END IF
            ENDIF
  300     CONTINUE
        END IF
  200   CONTINUE
  100 CONTINUE
C
C
C     Define right-hand operators: -(c/2) i.alpha x r
C     ===============================================
C     (PFAC(1:2) for Dirac mag. field, PFAC(3) for g_1 QED correction)
C
      IF (LONDON) THEN
        CALL DEF_D1HBLOND(1,IPRINT)
      ELSE
C
C  Conventional atomic orbitals (common gauge origin approach)
C
      PFAC(1)  = -0.5D0*CVAL
      PFAC(2)  = -0.5D0*CVAL
      PFAC(3)  =  0.5D0*G1FAC
      NCOMP = 2
C X component of B-field
      PNAME     = 'X magnetic field'
C     PNAME     = 'MAGFLDX'
      IPTYP     = IOP(1)
      PLABEL(1) = 'ZDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      IF (SHIELD) CALL OP1IND('XSHIELB',IPSHIELD(1+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IF (ESRGTENS) THEN
         CALL OP1IND('X magfield',IPESR(1+MXCOOR,1),
     &      LESRP,NESRP,INDXPR,MAXESR)
         PNAME    =  'g_1 MAGFLDX'
         IPTYP    =  13
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,PFAC(3),PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('X g_1 magfield',IPESR(1+MXCOOR,2),
     &      LESRP,NESRP,INDXPR,MAXESR)
      END IF
      IF (DOESR) THEN
C        Sigma_X:
         PNAME    =  'x-spin'
         IPTYP    =  10
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,0.5D0,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('x spin',IPSIGMA(1),LESRP,NESRP,INDXPR,MAXESR)
      END IF
#if !defined (XH4)
C Y component of B-field
      PNAME     = 'Y magnetic field'
C     PNAME     = 'MAGFLDY'
      IPTYP     = IOP(2)
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      IF (SHIELD) CALL OP1IND('YSHIELB',IPSHIELD(2+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IF (ESRGTENS) THEN
         CALL OP1IND('Y magfield',IPESR(2+MXCOOR,1),
     &      LESRP,NESRP,INDXPR,MAXESR)
         PNAME    =  'g_1 MAGFLDY'
         IPTYP    =  14
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,PFAC(3),PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('Y g_1 magfield',IPESR(2+MXCOOR,2),
     &      LESRP,NESRP,INDXPR,MAXESR)
      END IF
      IF (DOESR) THEN
C        Sigma_Y:
         PNAME    =  'y-spin'
         IPTYP    =  11
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,0.5D0,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('y spin',IPSIGMA(2),LESRP,NESRP,INDXPR,MAXESR)
      END IF
C Z component of B-field
      PNAME     = 'Z magnetic field'
C     PNAME     = 'MAGFLDZ'
      IPTYP     = IOP(3)
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      IF (SHIELD) CALL OP1IND('ZSHIELB',IPSHIELD(3+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IF (ESRGTENS) THEN
         CALL OP1IND('Z magfield',IPESR(3+MXCOOR,1),
     &      LESRP,NESRP,INDXPR,MAXESR)
         PNAME    =  'g_1 MAGFLDZ'
         IPTYP    =  15
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,PFAC(3),PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('Z g_1 magfield',IPESR(3+MXCOOR,2),
     &      LESRP,NESRP,INDXPR,MAXESR)
      END IF
      IF (DOESR) THEN
C        Sigma_Z:
         PNAME    =  'z-spin'
         IPTYP    =  12
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,0.5D0,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('z spin',IPSIGMA(3),LESRP,NESRP,INDXPR,MAXESR)
      END IF
#endif
      ENDIF
C
C     Print section for NMR shielding:
C
#ifdef MOD_LAO_REARRANGED
      IF (SHIELD .OR. shielding_rearrange) THEN
#else
      IF (SHIELD) THEN
#endif
C
C     Print section for NMR shieldings:
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'NMR shielding:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     Print reference
C
C     A operators
C     ===========
C
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(1X,A)')
     &   '** A operators for NMR shielding **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NSCOOR
          INDOP1 = IPSHIELD(I)
          IF(INDOP1.GT.0) THEN
            INDXPR = LLRAPU(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO
C
C     B operators
C     ===========
C
        CALL PRSYMB(LUPRI,'-',75,0)
        !WRITE(LUPRI,'(A)') ' ** B operators for NMR shielding or ESR **'
        IF (LONDON) THEN
         WRITE(LUPRI,'(1X,A)')
     &   '** B operators for NMR shielding, LAO approach **'
        ELSE
         WRITE(LUPRI,'(1X,A)')
     &   '** B operators for NMR shielding, CGO approach **'
        ENDIF
        CALL PRSYMB(LUPRI,'-',75,0)
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,3
          INDOP1 = IPSHIELD(I+MXCOOR)
          IF(INDOP1.GT.0) THEN
            INDXPR = LLRBPU(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO

      ENDIF
C
C     Print section for ESR parameters
C
      IF (ESRGTENS) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(A)') ' ** ESR g tensor operators **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,3
          INDOP1 = IPESR(I+MXCOOR,1)
          IF(INDOP1.GT.0) THEN
            INDXPR = LESRP(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO
      END IF
C
      IF (ESR_HFCC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(A)')
     &     ' ** ESR nuclear hyperfine coupling operators **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NSCOOR
          INDOP1 = IPESR(I,1)
          IF(INDOP1.GT.0) THEN
            INDXPR = LESRP(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO
      END IF
C
C     Expectation values from London orbitals (not diamagnetic terms!)
C     ===================================================================
C
      IF (LONDON) THEN
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(A)')
     & ' ** NMR shieldings expectation values, LAO **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NNSTIN
          INDOP1 = IPRNST(I)
          IF(INDOP1.GT.0) THEN
            INDXPR = LEXPP(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO
      ENDIF

      deallocate(DOATOM)
      CALL QEXIT('DEF_SHIELD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_OSCILL */
      SUBROUTINE DEF_OSCILL(IPRINT)
C***********************************************************************
C
C     Define operator for oscillator strength calculations
C
C     Written by S. Knecht Sep 2008
C
C***********************************************************************
      use consts_mod, only: D0, D1
!      use krciprop_mod, only: LPROP_KRCI, MXPROP_KRCI, NPROP_KRCI
#include "implicit.h"
      INTEGER :: II
C
!#include "maxash.h"
!#include "dcbgen.h"
!#include "dcbprp.h"
!#include "dcbxpr.h"
#include "krciprop.h"
!#include "dcborb.h"
!#include "dcbbas.h"
!#include "dgroup.h"
!#include "pgroup.h"
C
#ifdef CODE_NOT_WORKING
      CHARACTER PNAME(6)*16, PLABEL(6)*8
      DIMENSION PFAC(6),IOP(6)
C
C     define operator type, factor, name and label
C
      NPCOMP = 1
      PFAC(1)  = D1
      PFAC(2)  = D1
      PFAC(3)  = D1
      PFAC(4)  = D1
      PFAC(5)  = D1
      PFAC(6)  = D1
C     three components of dipole length integrals
      PNAME(1)  = 'X dipole length'
      PNAME(2)  = 'Y dipole length'
      PNAME(3)  = 'Z dipole length'
      PNAME(4)  = 'X dip vel'
      PNAME(5)  = 'Y dip vel'
      PNAME(6)  = 'Z dip vel'
      DO II  = 1, 6
        IOP(II) = 1
      END DO
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      PLABEL(3) = 'ZDIPLEN'
      PLABEL(4) = 'XDIPVEL'
      PLABEL(5) = 'YDIPVEL'
      PLABEL(6) = 'ZDIPVEL'
      DO II = 1, 6
         CALL XPRIND(PNAME(II),IOP(II),NPCOMP,PFAC(II),PLABEL(II),
     &              INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('KRCI PROPDIP',LUCIINDTRP,LPROP_KRCI,
     &               NPROP_KRCI,INDXPR,MXPROP_KRCI)
      END DO
#endif
      CHARACTER PNAME(3)*16, PLABEL(3)*8
      DIMENSION PFAC(3),IOP(3)
C
C     define operator type, factor, name and label
C
      NPCOMP = 1
      PFAC(1)  = 1.0d0 ! we introduce the factor -1.0d0 later in src/krmc/krmccan.F to enforce our "local" dipole length operator
      PFAC(2)  = 1.0d0 ! we introduce the factor -1.0d0 later in src/krmc/krmccan.F to enforce our "local" dipole length operator
      PFAC(3)  = 1.0d0 ! we introduce the factor -1.0d0 later in src/krmc/krmccan.F to enforce our "local" dipole length operator
C     three components of dipole length integrals
      PNAME(1)  = 'X dipole length'
      PNAME(2)  = 'Y dipole length'
      PNAME(3)  = 'Z dipole length'
      DO II  = 1, 3
        IOP(II) = 1
      END DO
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      PLABEL(3) = 'ZDIPLEN'
      DO II = 1, 3
         CALL XPRIND(PNAME(II),IOP(II),NPCOMP,PFAC(II),PLABEL(II),
     &              INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('KRCI PROPDIP',LUCIINDTRP,LPROP_KRCI,
     &               NPROP_KRCI,INDXPR,MXPROP_KRCI)
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_OMEGAQUANT */
      SUBROUTINE DEF_OMEGAQUANT(IPRINT)
C***********************************************************************
C
C     Define operator for determination of omega quantum numbers
C
C     Written by S. Knecht Sep 2008
C
C***********************************************************************
      use consts_mod, only: D1, DP5
!      use krciprop_mod, only: LPROP_KRCI, MXPROP_KRCI, NPROP_KRCI
#include "implicit.h"
      INTEGER :: II
C
!#include "maxash.h"
!#include "dcbgen.h"
!#include "dcbprp.h"
!#include "dcbxpr.h"
#include "krciprop.h"
!#include "dcborb.h"
!#include "dcbbas.h"
!#include "dgroup.h"
!#include "pgroup.h"
C
      CHARACTER PNAME(2)*16, PLABEL(2)*8
      DIMENSION PFAC(2),IOP(2)
C
C     define operator type, factor, name and label
C
      NPCOMP = 1
      PFAC(1)  = DP5
      PFAC(2)  = D1
C     define s_z and l_z as j_z = s_z + l_z
      PNAME(1) = 'z-spin'
      PNAME(2) = 'z-angular mom'
      IOP(1)   = 12
      IOP(2)   = 1
      PLABEL(1) = 'OVERLAP'
      PLABEL(2) = 'ZANGMOM'
      DO II = 1, 2
         CALL XPRIND(PNAME(II),IOP(II),NPCOMP,PFAC(II),PLABEL(II),
     &              INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('KRCI PROPOME',LUCIINDTRP,LPROP_KRCI,
     &               NPROP_KRCI,INDXPR,MXPROP_KRCI)
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_EEDM */
      SUBROUTINE DEF_EEDM(IPRINT)
C***********************************************************************
C
C     Define operator for electron EDM calculation
C
C     Written by T. Fleig, May 2012
C     and Malaya K. Nayak, Dec 2019
C
C***********************************************************************
      use consts_mod, only: D0
!      use dcbgen_mod, only: CVAL
!      use dcbxpr_mod, only: PRPNAM
!      use krciprop_mod, only: LPROP_KRCI, MXPROP_KRCI, NPROP_KRCI
#include "implicit.h"
C
!#include "maxash.h"
#include "dcbgen.h"
!#include "dcbprp.h"
#include "dcbxpr.h"
#include "krciprop.h"
!#include "dcborb.h"
!#include "dcbbas.h"
!#include "dgroup.h"
!#include "pgroup.h"
C
C     define operator type, factor, name and label
C     4*ic gamma^0 gamma^5 p^2/2m
C
      CHARACTER, PARAMETER :: PNAME*16 = 'P2-EDM  '
      CHARACTER, PARAMETER :: PLABEL*8 = 'KINENERG'
      INTEGER,   PARAMETER :: NPCOMP = 1, IOP = 23
      REAL*8               :: PFAC
C
CTF   Since we are working with the p^2 form of the EDM operator, there
CTF   is no explicit reference to Atom1 or Atom2 (etc.). We evaluate
CTF   over EDM integrals in the molecular spinor basis.
C
      PFAC = 4.0D0*CVAL
      CALL XPRIND(PNAME,IOP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('KRCI PROPEDM',LUCIINDTRP,LPROP_KRCI,
     &            NPROP_KRCI,INDXPR,MXPROP_KRCI)
      if (IPRINT.GE.1) then
        write(6,*)' (DEF_EEDM) testing ',LPROP_KRCI(NPROP_KRCI),
     &             NPROP_KRCI,INDXPR,PRPNAM(INDXPR)
        print*,'ISYXPR,ITRXPR,IPTYP,NPCOMP',ISYXPR,ITRXPR,IPTYP,NPCOMP
      end if
C The output should be
C   ISYXPR= 4  ITRXPR=-1  IPTYP=21  NPCOMP= 1
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_MHYP */
      SUBROUTINE DEF_MHYP(IPRINT)
C***********************************************************************
C
C     Define operator for magnetic hyperfine constants
C
C     Written by Malaya K. Nayak August 30 2013
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use consts_mod, only: D0
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbxpr_mod, only: PRPNAM
!      use krciprop_mod, only: LPROP_KRCI, MXPROP_KRCI, NPROP_KRCI
      use maxaqn_mod, only: MXAQN, MXQN
      use nuclei_mod, only: NUCIND
      use symmet_mod, only: IPTCNT, ISYMAX, MAXREP
      use chrxyz_mod, only: CHRXYZ
      use chrnos_mod, only: CHRNOS
#include "implicit.h"
C
!#include "maxash.h"
!#include "dcbgen.h"
!#include "dcbprp.h"
#include "dcbxpr.h"
#include "krciprop.h"
!#include "dgroup.h"
!#include "pgroup.h"
!#include "dcbxlr.h"
!#include "dcborb.h"
!#include "dcbbas.h"
      LOGICAL   SAME
      CHARACTER PNAME*16, PLABEL(2)*8
      REAL*8    PFAC(2)
      INTEGER   IOP(3)
C
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     MHYP: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)

C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
      PFAC(1)  = 1.0D0
      PFAC(2)  = 1.0D0
      NPCOMP   = 2
C
C     Irreps
C
      DO 100 IREPO = 0, MAXREP
C
C       Atoms
C
C       Define the Hyperfine operators
C
        DO 200 IATOM = 1, NUCIND
        IF(DOATOM(IATOM)) THEN
C
C         Cartesian directions
C
          DO 300 ICOOR = 1, 3
            ISCOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREPO,2)
            IF (ISCOR .GE. 1000) THEN
              CALL QUIT(
     &'DEF_MHYP1: Max 1000 coords in this implementation !')
            ENDIF
            IF (ISCOR .GT. 0) THEN
              IFIRST = ISCOR/100
              ISECND = MOD(ISCOR,100)/10
              ITHIRD = MOD(MOD(ISCOR,100),10)
C              PNAME    = 'HYP'//CHRNOS(IFIRST)
C     &                   //CHRNOS(ISECND)//CHRNOS(ITHIRD)
C     &                   //NAMN(IATOM)//REP(IREPO)//CHRXYZ(ICOOR)
              PNAME    = CHRXYZ(ICOOR)//CHRNOS(IATOM)//'-HYP'
              IPTYP    = IOP(ICOOR)
              IC1      = MOD(ICOOR+1,3) + 1
              IC2      = MOD(ICOOR  ,3) + 1
              IREP1    = IEOR(IREPO,ISYMAX(IC2,1))
              IREP2    = IEOR(IREPO,ISYMAX(IC1,1))
              IS1      = IPTCNT(3*(IATOM - 1) + IC1,IREP1,1)
              IS2      = IPTCNT(3*(IATOM - 1) + IC2,IREP2,1)
              IF((IS1.LE.0).OR.(IS2.LE.0))
     &          WRITE(LUPRI,'(A,5I5)')
     &          'DEF_MHYP1: sym. error !',IATOM,ICOOR,IC1,IC2
              IFIRST = IS1/100
              ISECND = MOD(IS1,100)/10
              ITHIRD = MOD(MOD(IS1,100),10)
              PLABEL(1)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              IFIRST = IS2/100
              ISECND = MOD(IS2,100)/10
              ITHIRD = MOD(MOD(IS2,100),10)
              PLABEL(2)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND('KRCI PROPHYP',LUCIINDTRP,LPROP_KRCI,
     &           NPROP_KRCI,INDXPR,MXPROP_KRCI)
              IF (IPRINT.GE.2) THEN
               write(6,*)'Hyperfine DEF Testing',LPROP_KRCI(NPROP_KRCI),
     &             NPROP_KRCI,INDXPR,PRPNAM(INDXPR)
              END IF
            ENDIF
  300     CONTINUE
        END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      deallocate(DOATOM)
      END ! DEF_MHYP
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_ENSPS */
      SUBROUTINE DEF_ENSPS(IPRINT)
C***********************************************************************
C
C     Define operator for electron-Nucleus Scalar-Pseudoscalar
C     P,T-odd interaction constant determining parameter
C
C     Written by Malaya K. Nayak Feb 26 2014
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use codata_mod, only: GFAC
      use pi_mod, only: PI
      use consts_mod, only: D1, D3, D4
!      use krciprop_mod, only: LPROP_KRCI, MXPROP_KRCI, NPROP_KRCI
!      use dcbxpr_mod, only: PRPNAM
      use nuclei_mod, only: GAUNUC, ISTBNU, NAMDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTNUC, MAXREP
      use chrnos_mod, only: CHRNOS
#include "implicit.h"
#include "codata.h"
!#include "maxash.h"
!#include "dcbgen.h"
!#include "dcbprp.h"
#include "krciprop.h"
!#include "dgroup.h"
!#include "pgroup.h"
#include "dcbxpr.h"
!#include "dcbexp.h"
C
      CHARACTER PNAME*16, PLABEL*8
      REAL*8    PFAC(2)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     ENSPS: INTTYP = 62 (nuclear charge density normalised)
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)

      IPTYP     = 23
      NPCOMP    = 1
      IF(GAUNUC) THEN
        PFAC(1) = D1
      ELSE
        PFAC(1) = D3/(D4*GFAC*PI)
      ENDIF
      DO 100 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 200 IREP = 0, MAXREP
               IF (IAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
C                  PNAME = 'SPS: '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
C     &                    //CHRNOS(IPTNUC(IATOM,IREP)/10)//
C     &                    CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  PNAME = 'A'//CHRNOS(IATOM)//'-SPS'
C....             Gaussian nucleus
                  IF(GAUNUC) THEN
                    PLABEL = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
C....             Point charge
                  ELSE
                    PLABEL = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  ENDIF
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
                  CALL OP1IND('KRCI PROPENSPS',LUCIINDTRP,LPROP_KRCI,
     &           NPROP_KRCI,INDXPR,MXPROP_KRCI)
       IF (IPRINT.GE.1) THEN
        write(6,*)'e-N P,T-odd DEF Testing',LPROP_KRCI(NPROP_KRCI),
     &             NPROP_KRCI,INDXPR,PRPNAM(INDXPR)
       END IF
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
C
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_NMQM */
      SUBROUTINE DEF_NMQM(IPRINT)
C***********************************************************************
C
C     Define nuclear quadruple moment interaction constant
C
C     Written by Malaya K. Nayak - March 20 2015
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbxpr_mod, only: PRPNAM
      use nuclei_mod, only: ISTBNU, NUCIND
      use symmet_mod, only: ISYMAX, MAXREP
      use chrnos_mod, only: CHRNOS
      use chrxyz_mod, only: CHRXYZ
!      use krciprop_mod, only: LPROP_KRCI, MXPROP_KRCI, NPROP_KRCI
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0, D3 = 3.0D0, D3INV = D1/D3
      CHARACTER PNAME*16, PLABEL(2)*8
      REAL*8    PFAC(2)
#include "dcbxpr.h"
!#include "dcbexp.h"
C
#include "krciprop.h"
C
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     NMQM: INTTYP = 30 (cart. el. field grad. at indiv. nuclei)
C     (-30: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-30)

      ITYP = 0
C     ZAVECTOR
      IPTYP     = 7
      NPCOMP    = 2
      PFAC(1)   = -D3INV
      PFAC(2)   = -D3INV
CMKN  Factor and Sign included, see Eq. (19) of reference on TaN P,T-odd
Cproperties.
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
                ICOOR1 = 2
                ICOOR2 = 3
            ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IAND(ISTBNU(IATOM),IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                  IOFF = IOFF + 1
                  ITYP = ITYP + 1
C                  PNAME        = 'EFG:
C                  '//CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)
C     *                 //CHRNOS(IATOM/10)
C     &                 //CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  PNAME    = CHRXYZ(ICOOR2)//CHRNOS(IATOM)//'-MQM'
                  PLABEL(1) = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  PLABEL(2) = CHRXYZ(ICOOR1-1)//CHRXYZ(ICOOR2)//'EFG'//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                 INDXPR,ISYXPR,ITRXPR,IPRINT)
                  CALL OP1IND('KRCI PROPNMQM',LUCIINDTRP,LPROP_KRCI,
     &                 NPROP_KRCI,INDXPR,MXPROP_KRCI)
C                  CALL OP1IND(CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'
C     *                 ,IPEFG(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
      IF(IPRINT.GE.2) THEN
        write(6,*)'NMQM DEF Testing',LPROP_KRCI(NPROP_KRCI),
     &             NPROP_KRCI,INDXPR,PRPNAM(INDXPR)
      ENDIF
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_PRP_NMQM */
      SUBROUTINE DEF_PRP_NMQM(LTRPP,NTRPP,MAXTRP,IPRINT)
C***********************************************************************
C
C     Define operator for nuclear magnetic quadruple moment constant
C
C     Written by Malaya K. Nayak July 19 2018
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbxpr_mod, only: PRPNAM
      use nuclei_mod, only: ISTBNU, NUCIND
      use symmet_mod, only: ISYMAX, MAXREP
      use chrnos_mod, only: CHRNOS
      use chrxyz_mod, only: CHRXYZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER   LTRPP,NTRPP,MAXTRP,IPRINT
#include "dcbxpr.h"
!#include "dcbexp.h"
C
C
      REAL(dp), PARAMETER :: D3INV = 1.0D0/3.0D0
      CHARACTER PNAME*16, PLABEL(2)*8
      REAL*8    PFAC(2)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     NMQM: INTTYP = 30 (cart. el. field grad. at indiv. nuclei)
C     (-30: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-30)
C
      ITYP = 0
CMKN  ZAVECTOR
      IPTYP     = 7
      NPCOMP    = 2
      PFAC(1)   = -D3INV
      PFAC(2)   = -D3INV
CMKN  Factor and Sign included, see Eq. (19) of reference on TaN P,T-odd properties.
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
            ICOOR1 = 2
            ICOOR2 = 3
            ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IAND(ISTBNU(IATOM),IEOR(IREPC,ISYMIJ)).EQ.0) THEN
                  IOFF = IOFF + 1
                  ITYP = ITYP + 1
                  PNAME    = CHRXYZ(ICOOR2)//CHRNOS(IATOM)//'-MQM'
                  PLABEL(1) = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  PLABEL(2) = CHRXYZ(ICOOR1-1)//CHRXYZ(ICOOR2)//'EFG'//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                        INDXPR,ISYXPR,ITRXPR,IPRINT)
                  CALL OP1IND('TRPINP NMQM',INDTRP,LTRPP,NTRPP,
     &                        INDXPR,MAXTRP)
      IF(IPRINT.GE.1) THEN
      write(6,*)'NMQM DEF Testing',LTRPP,NTRPP,INDXPR,PRPNAM(INDXPR)
      ENDIF
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_PVC_SHIELD */
      SUBROUTINE DEF_PVC_SHIELD(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operators for PV contributions to NMR shieldings
C-----------------------------------------------------------------------
C     Radovan Bast                         Revision:      september 2006
C     Agustin Aucar                        Last revision: november  2022
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use priunit_mod, only: LUPRI
      use consts_mod, only: D0, D3, D4
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use codata_mod, only: GFAC
      use pi_mod, only: PI
      use nuclei_mod, only: GAUNUC, ISTBNU, NAMDEP, NUCDEP, NUCIND,
     &  NUCPRE
      use symmet_mod, only: IPTCNT, IPTNUC, MAXREP
!      use dcbgen_mod, only: CVAL
!      use dcbprp_mod, only: IP_PVC_SHIELD
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
      use chrnos_mod, only: CHRNOS
      use memory_allocator
C
      implicit none
C
      integer I,J
C
#include "codata.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
C

      INTEGER IPRINT,IATOM,IXYZ,NATOM,ISCOOR,NSCOOR,NPCOMP,
     &        INDOP1,INDXPR,IOP,IPTYP,IREP,ISYXPR,ITRXPR
C
      DOUBLE PRECISION PFAC
C
      logical, allocatable :: DOATOM(:)
C
      CHARACTER PNAME*16, PLABEL(2)*8
      DIMENSION PFAC(2),IOP(3)
C
      CALL QENTER('DEF_PVC_SHIELD')
C
      NSCOOR = 3*NUCDEP
C
      CALL IZERO(IP_PVC_SHIELD,3+MXCOOR)
C
      allocate(DOATOM(NUCIND))
C
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
C
C     Here we set up the operators needed for the linear reponse calculation
C
C     Define left-hand operators: nuclear-spin-dependent PV: c i.alpha rho(r)
C     =======================================================================
C
C     XALPHA
      IOP(1) = 2
C     YALPHA
      IOP(2) = 3
C     ZALPHA
      IOP(3) = 4
C
      NPCOMP  = 1    ! one component
      PFAC(1) = CVAL ! factor
C
C     this prefactor (CVAL) guarantees the NRL employes in the literature
C     in the printout section we will divide by CVEL
C
C     Point nucleus
      IF(.NOT.GAUNUC) PFAC(1) = PFAC(1)*D3/(D4*GFAC*PI)
C
      DO IATOM = 1,NUCDEP
        IF(DOATOM(IATOM)) THEN
          DO IREP = 0,MAXREP
            IF(IAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
C
              IF(GAUNUC) THEN ! Gaussian nuclear model
                PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ELSE ! Point nucleus
                PLABEL(1) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ENDIF
              PNAME = 'PVCa '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                       //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                       //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
C
              DO IXYZ = 1,3
                ISCOOR = IPTCNT(3*(IATOM-1)+IXYZ,IREP,2)
C
                IF(ISCOOR.GT.0) THEN
                IPTYP = IOP(IXYZ)
                PNAME(5:5) = CHRNOS(IXYZ)
                CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                      INDXPR,ISYXPR,ITRXPR,IPRINT)
                CALL OP1IND(PNAME(1:6),IP_PVC_SHIELD(ISCOOR),
     &                      LLRAPU,NLRAPT,INDXPR,MAXLLR)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C
C     Define right-hand operators: -(c/2) i.alpha x r
C     ===============================================
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
C
      NPCOMP  = 2              ! two components
      PFAC(1) = -0.5D0*CVAL    ! factor
      PFAC(2) = -0.5D0*CVAL    ! factor
C
C     x component of B-field
      PNAME     = 'X magnetic field'
C     PNAME     = 'MAGFLDX'
      IPTYP     = IOP(1)
      PLABEL(1) = 'ZDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XSHIELB',IP_PVC_SHIELD(1+MXCOOR),
     &            LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
C     y component of B-field
      PNAME     = 'Y magnetic field'
C     PNAME     = 'MAGFLDY'
      IPTYP     = IOP(2)
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YSHIELB',IP_PVC_SHIELD(2+MXCOOR),
     &            LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
C     z component of B-field
      PNAME     = 'Z magnetic field'
C     PNAME     = 'MAGFLDZ'
      IPTYP     = IOP(3)
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('ZSHIELB',IP_PVC_SHIELD(3+MXCOOR),
     &            LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
C     Print section for parity-violation contribution to NMR shielding:
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Parity-violation contribution to '//
     &  'NMR shielding:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     Print reference
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   '** A operators for PVC to NMR shielding **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOOR
        INDOP1 = IP_PVC_SHIELD(I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRAPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   '** B operators for PVC to NMR shielding **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3
        INDOP1 = IP_PVC_SHIELD(I+MXCOOR)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C
      deallocate(DOATOM)
C
      CALL QEXIT('DEF_PVC_SHIELD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DEFINE_PVC_SPINSPIN */
      SUBROUTINE DEFINE_PVC_SPINSPIN()
C=======================================================================
C     define operators for
C     PV contribution to spin-spin coupling
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
C
      logical, allocatable :: DOATOM(:)
C
C
      allocate(DOATOM(NUCIND))
      CALL SETATM(DOATOM,NATOM,-29)
      CALL DEFINE_PVC_SPINSPIN_SUB(DOATOM)
      deallocate(DOATOM)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DEFINE_PVC_SPINSPIN_SUB */
      SUBROUTINE DEFINE_PVC_SPINSPIN_SUB(DOATOM)
C=======================================================================
C     define operators for
C     PV contribution to spin-spin coupling
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use codata_mod, only: GFAC
      use pi_mod, only: PI
      use nuclei_mod, only: GAUNUC, ISTBNU, NAMDEP, NAMN, NUCDEP,
     &  NUCIND, NUCPRE
!      use pgroup_mod, only: REP
      use symmet_mod, only: IPTCNT, IPTNUC, ISYMAX, MAXREP
!      use dcbgen_mod, only: CVAL
!      use dcbprp_mod, only: IP_PVC_SPINSPIN
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
      use chrxyz_mod, only: CHRXYZ
      use chrnos_mod, only: CHRNOS
#include "implicit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "codata.h"
!#include "dgroup.h"
#include "pgroup.h"
#include "dcbgen.h"
#include "dcbprp.h"
!#include "dcbesr.h"
#include "dcbxlr.h"
!#include "dcbexp.h"
!#include "dcbnmr.h"
C
      LOGICAL   DOATOM(NUCIND)
      CHARACTER PNAME*16,
     &          PLABEL(3)*8
      DIMENSION PFAC(3),IOP(3)
C
C
C
C     right hand
C     ==========
C
      CALL IZERO(IP_PVC_SPINSPIN,2*MXCOOR)
C
C     xavector
      IOP(1) = 5
C     yavector
      IOP(2) = 6
C     zavector
      IOP(3) = 7
C
      PFAC(1)  = CVAL
      PFAC(2)  = CVAL
C
      IREP = 0
      DO IATOM = 1,NUCDEP
        IF(DOATOM(IATOM)) THEN
          DO IXYZ = 1,3
            ISCOOR = IPTCNT(3*(IATOM-1)+IXYZ,IREP,2)
            IF(ISCOOR.GT.0) THEN
C
              IFIRST = ISCOOR/100
              ISECND = MOD(ISCOOR,100)/10
              ITHIRD = MOD(MOD(ISCOOR,100),10)
              PNAME    = 'NUCSP'//CHRNOS(IFIRST)
     &                   //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                   //NAMN(IATOM)//REP(IREP)//CHRXYZ(IXYZ)
C
              IPTYP = IOP(IXYZ)
              IC1   = MOD(IXYZ+1,3)+1
              IC2   = MOD(IXYZ  ,3)+1
              IREP1 = IEOR(IREP,ISYMAX(IC2,1))
              IREP2 = IEOR(IREP,ISYMAX(IC1,1))
              IS1   = IPTCNT(3*(IATOM-1)+IC1,IREP1,1)
              IS2   = IPTCNT(3*(IATOM-1)+IC2,IREP2,1)
C
              IFIRST = IS1/100
              ISECND = MOD(IS1,100)/10
              ITHIRD = MOD(MOD(IS1,100),10)
              PLABEL(1)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              IFIRST = IS2/100
              ISECND = MOD(IS2,100)/10
              ITHIRD = MOD(MOD(IS2,100),10)
              PLABEL(2)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
C
              CALL XPRIND(PNAME,IPTYP,2,PFAC,PLABEL,
     &                    INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND('XSPNA',INDEXB,LLRBPU,NLRBPT,INDXPR,MAXLLR)
              IP_PVC_SPINSPIN(1,ISCOOR) = INDEXB
              CALL WRIXPR(INDEXB,INDXPR)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C
C     left hand
C     =========
C
      WRITE(LUPRI,'(A)') 'Parity-violating contribution'
      WRITE(LUPRI,'(A)') 'to spin-spin coupling'
C
C     xalpha
      IOP(1) = 2
      PFAC(1) = CVAL
C     yalpha
      IOP(2) = 3
      PFAC(2) = CVAL
C     zalpha
      IOP(3) = 4
      PFAC(3) = CVAL
C
C     this prefactor (CVAL) guarantees the NRL employes in the literature
C     in the printout section we will divide by CVEL
C
      IF(.NOT.GAUNUC) THEN ! Point nucleus
        DO IXYZ = 1,3
          PFAC(IXYZ) = PFAC(IXYZ)*D3/(D4*GFAC*PI)
        ENDDO
      ENDIF
C
      DO IATOM = 1,NUCDEP
        IF(DOATOM(IATOM)) THEN
          DO IREP = 0,MAXREP
            IF(IAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
C
              IF(GAUNUC) THEN ! Gaussian nuclear model
                PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ELSE ! Point nucleus
                PLABEL(1) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ENDIF
              PNAME = 'PVCa '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                       //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                       //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
C
              DO IXYZ = 1,3
C
                IPTYP = IOP(IXYZ)
                PNAME(5:5) = CHRNOS(IXYZ)
                CALL XPRIND(PNAME,IPTYP,1,PFAC,PLABEL,
     &                      INDXPR,ISYXPR,ITRXPR,IPRINT)
                CALL OP1IND(PNAME(1:6),INDEXA,LLRAPU,NLRAPT,
     &                      INDXPR,MAXLLR)
                ISCOOR = IPTCNT(3*(IATOM-1)+IXYZ,IREP,2)
                IF(ISCOOR.GT.0) THEN
                  IP_PVC_SPINSPIN(2,ISCOOR) = INDEXA
                ENDIF
                CALL WRIXPR(INDEXA,INDXPR)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_PVC_EFG */
      SUBROUTINE DEF_PVC_EFG(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operators for PV contributions to EFG
C
C     Written by Juan J. Aucar - March 2025
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: GAUNUC, ISTBNU, NAMDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTCNT, IPTNUC, ISYMAX, MAXREP
      use chrxyz_mod, only: CHRXYZ
!      use dcbprp_mod, only: IP_PVC_EFG
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
      use consts_mod, only: D1, D3, D4
!      use codata_mod, only: GFAC
      use pi_mod, only: PI
      use chrnos_mod, only: CHRNOS
      use memory_allocator
C
      implicit none
C
      integer I,J
C
!#include "pgroup.h"
#include "dcbprp.h"   ! common block for expectation values
#include "dcbxlr.h"   ! common block for linear reponse
!#include "dcbgen.h"
#include "codata.h"
!#include "dcbexp.h"
C
      INTEGER IPRINT,IATOM,IXYZ,NATOM,
     &        IC1,IC2,ICOOR,IFIRST,INDOP1,
     &        INDXPR,IOP,IPTYP,IREP,IREP1,IREP2,IS1,
     &        IS2,ISCOOR,ISECND,NSCOORA,NSCOORB, INDEXA,
     &        INDEXB, ISYXPR,ITRXPR,NPCOMP,
     &        ITYP, ICOOR1, ICOOR2, ISYMIJ, IOFF, IREPC
C
      DOUBLE PRECISION PFAC(NUCIND)
C
      logical, allocatable :: DOATOM(:)
C
      CHARACTER PNAME*16, PLABEL(2)*8
      DIMENSION IOP(1)
C
      CALL QENTER('DEF_PVC_EFG')
C
      NSCOORA = 9*NUCIND
      NSCOORB = NUCIND
C
      CALL IZERO(IP_PVC_EFG,MXCOOR)
      allocate(DOATOM(NUCIND))
C
C     (accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)

C     Open print section for PVC to EFG
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Parity-violation contribution to '//
     &'electric field gradient:'
      CALL PRSYMB(LUPRI,'=',75,0)


C
C     Define left-hand operators: EFG components
C     ===========================================
C
        PFAC = 1
        ITYP = 0
        IPTYP     = 1
        NPCOMP    = 1
        PFAC   = D1
        DO 100 IATOM = 1, NUCIND
        IF (.NOT. DOATOM(IATOM)) GOTO 100
        DO 200 ICOOR1 = 1, 3
        DO 200 ICOOR2 = ICOOR1, 3
          ISYMIJ = IEOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
          IOFF = 0
          DO 300 IREPC = 0, MAXREP
            IF (IAND(ISTBNU(IATOM),IEOR(IREPC,ISYMIJ)).EQ.0) THEN
              IOFF = IOFF + 1
              ITYP = ITYP + 1
              PNAME = 'EFG: '//CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)
     &                       //CHRNOS(IATOM/10)
     &                       //CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
              PLABEL(1) = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     &        CHRNOS(IATOM/10)//
     &        CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)

!             Set A operators
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND(CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG',
     &        INDEXA, LLRAPU,NLRAPT,INDXPR,MAXLLR)
              ISCOOR = IPTCNT((IATOM),IREPC,2)
              IF(ISCOOR.GT.0) THEN
                IP_PVC_EFG(1,ITYP) = INDEXA
              ENDIF
C              CALL WRIXPR(INDEXB,INDXPR)
            END IF
 300      CONTINUE
 200    CONTINUE
 100  CONTINUE
C
C
C     Define right-hand operators: nuclear-spin-independent PV
C     =======================================================
C
      ITYP = 0
      NPCOMP  = 1
      PFAC = 1
C
      IF(.NOT.GAUNUC) THEN ! Point nucleus
          PFAC = PFAC *D3/(D4*GFAC*PI)
      ENDIF
C
      DO IATOM = 1,NUCIND
        IF(DOATOM(IATOM)) THEN
          DO IREP = 0,MAXREP
            IF(IAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
C
              IF(GAUNUC) THEN ! Gaussian nuclear model
                PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ELSE ! Point nucleus
                PLABEL(1) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ENDIF
              PNAME = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                       //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                       //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
C             Set B operators
              IPTYP = 9
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                    INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND(PNAME(1:6),INDEXB,LLRBPU,NLRBPT,
     &                    INDXPR,MAXLLR)
              ISCOOR = IPTCNT((IATOM),IREP,2)
              IF(ISCOOR.GT.0) THEN
                IP_PVC_EFG(2,ISCOOR) = INDEXB
              ENDIF
C              CALL WRIXPR(INDEXA,INDXPR)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C     Print section for A and B operators
C     ===================================
C
C     A operators
C     -----------
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     & '** A operators for PV electric field gradient **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOORA
        INDOP1 = IP_PVC_EFG(1,I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRAPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     B operators
C     -----------
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &'** B operators for PV electric field gradient **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOORB
        INDOP1 = IP_PVC_EFG(2,I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C
      deallocate(DOATOM)
C
      CALL QEXIT('DEF_PVC_EFG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_NSTDIA */
      SUBROUTINE DEF_NSTDIA(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define diamagnetic contribution to nuclear shielding tensor (NST)
C     (based on the Sternheimer approximation, see
C      Aucar et al., JCP 110 (1999) 6208) )
!
!     Written by T. Enevoldsen Oct. 1997,
!     updates by MI&HJJ 2003,2004
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     NSTDIACGO: INTTYP = 38
C     (-38: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-38)
      CALL DEF_NSTDIA1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_NSTDIA1 */
      SUBROUTINE DEF_NSTDIA1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define diamagnetic contribution to shieldings
C     Written by T. Enevoldsen Oct. 1997
C
C  MI&HJJ/2002,2003 - extended for LAO
C
C***********************************************************************
      use cbiher_mod, only: TRIANG
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NAMN, NUCIND
      use symmet_mod, only: IPTCNT, MAXREP
!      use dcbprp_mod, only: IPNSTDIA, LONDON, NNSTIN, TRIANG
!      use dcbexp_mod, only: LEXPP, MAXEXP, NEXPP
      use chrxyz_mod, only: CHRXYZ
      use chrnos_mod, only: CHRNOS
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0
C
!#include "dgroup.h"
!#include "pgroup.h"
!#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbexp.h"
      LOGICAL DOATOM(NUCIND),SAME
      CHARACTER PNAME*16, PLABEL(3)*8, LABEL*4
      DIMENSION PFAC(3),IOP(3)

      CALL QENTER('DEF_NSTDIA1')

      IF (LONDON) THEN
         LABEL = 'NSNL'
      ELSE
         LABEL = 'NSCO'
      END IF

C  always do control print out
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)')
     &' NMR diamagnetic shielding contributions '//
     &' (Aucar et al., JCP 110 (1999) 6208)'
       IF (LONDON) THEN
        WRITE(LUPRI,'(2X,A)')
     &'- for London atomic orbitals (invariant gauge origin)'
       ELSE
        WRITE(LUPRI,'(2X,A)')
     &'- for conventional atomic orbitals (common gauge origin)'
       ENDIF
      CALL PRSYMB(LUPRI,'=',75,0)
C
      NNSTIN  = 0
      IPTYP   = 1
      NPCOMP  = 1
      PFAC(1) = D1
C
C     Irreps
C
      DO 100 IREPO = 0, MAXREP
C
C       Atoms
C
        DO 200 IATOM1 = 1, NUCIND
        IF (DOATOM(IATOM1)) THEN
          MXATM2 = NUCIND
          IF (TRIANG) MXATM2 = IATOM1
C
C         Cartesian directions
C
          DO 500 ICOOR1 = 1, 3
            ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
            IF (ISCOR1 .GT. 0) THEN
              DO 600 ICOOR2 = 1, 3
                NNSTIN = NNSTIN + 1
                IFIRST = ISCOR1/100
                ISECND = MOD(ISCOR1,100)/10
                ITHIRD = MOD(MOD(ISCOR1,100),10)
                PNAME = LABEL//NAMN(IATOM1)//CHRNOS(IFIRST)
     &                  //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                  //CHRXYZ(ICOOR2)
                PLABEL(1) = CHRNOS(IFIRST)//CHRNOS(ISECND)//
     &                      CHRNOS(ITHIRD)//LABEL//CHRXYZ(ICOOR2)
                CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                      INDXPR,ISYXPR,ITRXPR,IPRINT)
                IF(IPRINT.GE.0) CALL WRIXPR(NNSTIN,INDXPR)
                CALL OP1IND(LABEL,IPIND,LEXPP,NEXPP,
     &                      INDXPR,MAXEXP)
                IPNSTDIA(NNSTIN) = IPIND
 600          CONTINUE
            END IF
 500      CONTINUE
        END IF
 200    CONTINUE
 100  CONTINUE

      CALL QEXIT('DEF_NSTDIA1')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_msuscp */
      SUBROUTINE DEF_MSUSCP()
C****************************************************************************
C
C    PURPOSE:
C  ===========
C   Define operators for calculating the magnetic susceptibility
C   both for response part and the expectation value (& diamagnetic) terms.
C
C   Written by MI, jan.2003 (based on old TEC's routines)
C   Last modifications: MI, April 2006 (for IMAGOUT)
C
C****************************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbgen_mod, only: CVAL
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
!      use dcbxpr_mod, only: PRPNAM
!      use dcbprp_mod, only: DOEXP, IMAGOUT, IPLONDON, IPMSUSC, IPRPRP, LONDON, MSUSCDIA
!      use dcbexp_mod, only: IPDSUSC, LEXPP, MAXEXP, NEXPP
!      use dcbham_mod, only: RDSUSLLMOD
#include "implicit.h"
      INTEGER :: I, J
C
      PARAMETER(DI2 = 0.5D0, DI4 = 0.25D0, DM1 = -1.0D0,
     &           D2 = 2.0D0, DM2 = -2.0D0,  D1 = 1.0D0)
#include "dcbgen.h"

#include "dcbxlr.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbexp.h"
!#include "dcbsusc.h"
#include "dcbham.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3),IOP(3)

      CALL QENTER('DEF_MSUSCP')
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     &  'DEF_MSUSCP: Magnetic susceptibility'
       IF (LONDON) THEN
          WRITE(LUPRI,'(1X,A)')
     &  '- with London atomic orbitals (LAO).'
       ELSE
          WRITE(LUPRI,'(1X,A)')
     &  '- with conventional atomic orbitals (CGO)'
       ENDIF
       CALL PRSYMB(LUPRI,'=',75,0)

       IF (IMAGOUT.NE.0) THEN
         IF (IMAGOUT.NE.1.AND.IMAGOUT.NE.2.AND.IMAGOUT.NE.3) THEN
          WRITE(LUPRI,*)
     &    '...IMAGOUT=',IMAGOUT,' must be 1/2/3, set to 0'
          IMAGOUT=0
         ELSE
          WRITE(LUPRI,*)
     &    '...hard setting of the output for magnetizabilities'//
     &    '...IMAGOUT=',IMAGOUT
         ENDIF
       ENDIF
C
C  ** Split for conventional and London orbitals
C
      IPRINT = IPRPRP
      CALL IZERO(IPMSUSC,6)
C
      IF (LONDON) THEN ! use London AO's
C ... always calculate expectation value terms
       DOEXP = .TRUE.
C ... define operators for response and exp.value parts
       CALL DEF_D1HBLOND(2,LPRINT)
      ELSE ! Conventional AO
C      ... XAVECTOR
       IOP(1) = 5
C      ... YAVECTOR
       IOP(2) = 6
C      ... ZAVECTOR
       IOP(3) = 7
C
C     Define operators for conventional orbitals
       NPCOMP   = 2
       PFAC(1)  = -DI2*CVAL
       PFAC(2)  = -DI2*CVAL
C   ....X componete of B-field
C      PNAME     = 'MAGFLDX'
       PNAME     = 'X magnetic field'
       IPTYP     = IOP(1)
       PLABEL(1) = 'ZDIPLEN'
       PLABEL(2) = 'YDIPLEN'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XFIELDA',IPMSUSC(1,1),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
       CALL OP1IND('XFIELDB',IPMSUSC(2,1),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
#if !defined (XH4)
C Y componete of B-field
C      PNAME     = 'MAGFLDY'
       PNAME     = 'Y magnetic field'
       IPTYP     = IOP(2)
       PLABEL(1) = 'XDIPLEN'
       PLABEL(2) = 'ZDIPLEN'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XFIELDA',IPMSUSC(1,2),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
       CALL OP1IND('XFIELDB',IPMSUSC(2,2),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
C Z componete of B-field
C      PNAME     = 'MAGFLDZ'
       PNAME     = 'Z magnetic field'
       IPTYP     = IOP(3)
       PLABEL(1) = 'YDIPLEN'
       PLABEL(2) = 'XDIPLEN'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XFIELDA',IPMSUSC(1,3),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
       CALL OP1IND('XFIELDB',IPMSUSC(2,3),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
#endif
       ENDIF

C----------------------------------------------------------
C   Print section of operators for both LAO and conventional AO
C----------------------------------------------------------
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     & '** A operators for magnetic susceptibility **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3
       IF (LONDON) THEN
         DO J = 1, 4
           CALL WRIXPR(0,IPLONDON(J,I))
         ENDDO
         WRITE(LUPRI,'(2X,A)')
     &'Including g(B) and {T(1),H(0)} terms. '
       ELSE
       INDOP1 = IPMSUSC(1,I)
       IF(INDOP1.GT.0) THEN
         INDXPR = LLRAPU(INDOP1)
         CALL WRIXPR(INDOP1,INDXPR)
       ENDIF
       ENDIF
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     & '** B operators for magnetic susceptibility **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1, 3
       IF (LONDON) THEN
         DO J = 1, 4
           CALL WRIXPR(0,IPLONDON(J,I))
         ENDDO
         WRITE(LUPRI,'(2X,A)')
     &'Including g(B) and {T(1),H(0)} terms. '
       ELSE
        INDOP1 = IPMSUSC(2,I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
       ENDIF
      ENDDO

C ----------------------------------------------
C Define operators for the expectation value
C ----------------------------------------------
      IF (MSUSCDIA.AND.(.NOT.LONDON)) THEN
       DOEXP = .TRUE. ! Calculate the expectation value
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     & 'Diamagnetic susceptibility terms (CGO):'
      CALL PRSYMB(LUPRI,'=',75,0)
C
      PNAME     = 'DSUSC_CGO-XX'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XXSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XXDSUSC',IPDSUSC(1),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-XY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XYSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XYDSUSC',IPDSUSC(2),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-XZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XZDSUSC',IPDSUSC(3),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-YY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YYSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YYDSUSC',IPDSUSC(4),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-YZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YZDSUSC',IPDSUSC(5),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-ZZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'ZZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('ZZDSUSC',IPDSUSC(6),LEXPP,NEXPP,INDXPR,MAXEXP)
C
C      Print section for CGO diamag.susceptibilities
C    ==================================================
C
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,'(1X,A)')
     & 'Diamagnetic susceptibility expectation value terms (CGO):'
        DO I = 1, 6
          INDXPR = LEXPP(IPDSUSC(I))
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF

      ELSE IF ((.NOT.MSUSCDIA).AND.(.NOT.LONDON)) THEN
C  ... Do not calculate exp. value - diamag.term
C  which is fully  involved in paramag. term when rotating electrons and positrons

        DOEXP = .FALSE.
        WRITE(LUPRI,'(1X,A)')
     & 'Diamagnetic CGO susceptibility expectation value term'//
     & ' is  not calculated (MSUSCDIA=.false.)!'

! ... London atomic orbitals - always calculate the expectation value
      ELSE IF (LONDON) THEN
C  ... Here take care of one-electron terms
C ... first 3 terms from the <r r~ hd> - RM2H1+RM2H2+RM2H3
C  ... 6 terms  Qmn.<r r~ beta' c^2>.Qmn - 6 diagonal operators
C     RM2H1 terms, IPDSUSC(1-6)

!MI note, Oct.2012: factor -1/2 = (1/4)*(-2),where 1/4 is from integral
! definition (1/4Q_MN<rr^T h_D>Q_MN),-2 is in beta' prescription

        PNAME     = 'LAO-XXRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS_dB2XX'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XXRM2H1',IPDSUSC(1),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS_dB2XY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XYRM2H1',IPDSUSC(2),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS_dB2XZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XZRM2H1',IPDSUSC(3),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS_dB2YY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YYRM2H1',IPDSUSC(4),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS_dB2YZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YZRM2H1',IPDSUSC(5),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS_dB2ZZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZZRM2H1',IPDSUSC(6),LEXPP,NEXPP,INDXPR,MAXEXP)

C ... RM2H2 terms, IPDSUSC(7-12)
!MI note,Oct.2012: in abacus(routine RM2H2) they already have factor 1/4
        PNAME     = 'LAO-XXRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XXRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XXRM2H2',IPDSUSC(7),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XYRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XYRM2H2',IPDSUSC(8),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XZRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XZRM2H2',IPDSUSC(9),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YYRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YYRM2H2',IPDSUSC(10),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YZRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YZRM2H2',IPDSUSC(11),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'ZZRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZZRM2H2',IPDSUSC(12),LEXPP,NEXPP,INDXPR,MAXEXP)

C -----------------------------------
C ... RM2H3 terms, IPDSUSC(13-18)
! MI note, Oct.2012: factor of 1/4 is in corresponding abacus routine RM2H3
C -----------------------------------
        PNAME     = 'LAO-XXRM2H3'
        IPTYP     = 8
        NPCOMP    = 3
        PFAC(1)  = CVAL
        PFAC(2)  = CVAL
        PFAC(3)  = CVAL
        PLABEL(1) = 'XXRM2H3X'
        PLABEL(2) = 'XXRM2H3Y'
        PLABEL(3) = 'XXRM2H3Z'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XXRM2H3',IPDSUSC(13),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRM2H3'
        IPTYP     = 8
        NPCOMP    = 3
        PFAC(1)  = CVAL
        PFAC(2)  = CVAL
        PFAC(3)  = CVAL
        PLABEL(1) = 'XYRM2H3X'
        PLABEL(2) = 'XYRM2H3Y'
        PLABEL(3) = 'XYRM2H3Z'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XYRM2H3',IPDSUSC(14),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRM2H3'
        IPTYP     = 8
        NPCOMP    = 3
        PFAC(1)  = CVAL
        PFAC(2)  = CVAL
        PFAC(3)  = CVAL
        PLABEL(1) = 'XZRM2H3X'
        PLABEL(2) = 'XZRM2H3Y'
        PLABEL(3) = 'XZRM2H3Z'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XZRM2H3',IPDSUSC(15),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRM2H3'
        IPTYP     = 8
        NPCOMP    = 3
        PFAC(1)  = CVAL
        PFAC(2)  = CVAL
        PFAC(3)  = CVAL
        PLABEL(1) = 'YYRM2H3X'
        PLABEL(2) = 'YYRM2H3Y'
        PLABEL(3) = 'YYRM2H3Z'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YYRM2H3',IPDSUSC(16),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRM2H3'
        IPTYP     = 8
        NPCOMP    = 3
        PFAC(1)  = CVAL
        PFAC(2)  = CVAL
        PFAC(3)  = CVAL
        PLABEL(1) = 'YZRM2H3X'
        PLABEL(2) = 'YZRM2H3Y'
        PLABEL(3) = 'YZRM2H3Z'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YZRM2H3',IPDSUSC(17),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRM2H3'
        IPTYP     = 8
        NPCOMP    = 3
        PFAC(1)  = CVAL
        PFAC(2)  = CVAL
        PFAC(3)  = CVAL
        PLABEL(1) = 'ZZRM2H3X'
        PLABEL(2) = 'ZZRM2H3Y'
        PLABEL(3) = 'ZZRM2H3Z'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZZRM2H3',IPDSUSC(18),LEXPP,NEXPP,INDXPR,MAXEXP)

C ... RDSUSLL terms,  IPDSUSC(19-27)
C ... off-diagonal tensor elements have to be symmetrized !
! see abacus/RDSUSAN, where they are multiplied with 1/2 factor
        IF (RDSUSLLMOD) THEN
           ! factor for Dalton-nonrelatic form
           FACLL = -2.0D0
        ELSE
!radovan:  why is this factor different?
!          in which situation is this needed?
!miro: this is to be carefully examined, at least this is
! factor of the full relativistic form of the RDSUSLL operator
           FACLL = -1.0D0
        ENDIF

        PNAME     = 'LAO-XXRDSUSLL'
        IPTYP     = 5
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'XXRDSULZ'
        PLABEL(2) = 'XXRDSULY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XXRDSUL',IPDSUSC(19),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRDSUSLL'
        IPTYP     = 6
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'XYRDSULX'
        PLABEL(2) = 'XYRDSULZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XYRDSUL',IPDSUSC(20),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRDSUSLL'
        IPTYP     = 7
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'XZRDSULY'
        PLABEL(2) = 'XZRDSULX'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XZRDSUL',IPDSUSC(21),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YXRDSUSLL'
        IPTYP     = 5
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'YXRDSULZ'
        PLABEL(2) = 'YXRDSULY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YXRDSUL',IPDSUSC(22),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRDSUSLL'
        IPTYP     = 6
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'YYRDSULX'
        PLABEL(2) = 'YYRDSULZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YYRDSUL',IPDSUSC(23),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRDSUSLL'
        IPTYP     = 7
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'YZRDSULY'
        PLABEL(2) = 'YZRDSULX'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YZRDSUL',IPDSUSC(24),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZXRDSUSLL'
        IPTYP     = 5
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'ZXRDSULZ'
        PLABEL(2) = 'ZXRDSULY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZXRDSUL',IPDSUSC(25),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZYRDSUSLL'
        IPTYP     = 6
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'ZYRDSULX'
        PLABEL(2) = 'ZYRDSULZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZYRDSUL',IPDSUSC(26),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRDSUSLL'
        IPTYP     = 7
        NPCOMP    = 2
        PFAC(1)   = FACLL*CVAL
        PFAC(2)   = FACLL*CVAL
        PLABEL(1) = 'ZZRDSULY'
        PLABEL(2) = 'ZZRDSULX'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZZRDSUL',IPDSUSC(27),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)
!-----------------------------------------------------------------
! RDSUSNL terms,  pointers  IPDSUSC(28-33),
! represent gauge-origin independent Sternheim diamagnetic term
! for LAO-magnetizabilities
!
! include them 1.RKB-by default,  2.UKB - when specified (.MSUSCD)
!-----------------------------------------------------------------
      IF (MSUSCDIA) THEN
        PNAME     = 'LAO-XXRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XXDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XXDSUSN',IPDSUSC(28),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XYDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XYDSUSN',IPDSUSC(29),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XZDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XZDSUSN',IPDSUSC(30),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YYDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YYDSUSN',IPDSUSC(31),
     &                LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YZDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YZDSUSN',IPDSUSC(32),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'ZZDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZZDSUSN',IPDSUSC(33),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

C ... for diamagnetic terms, we also need DSUSCGM DALTON analogues, pointers IPDSUSC(34-39)

        PNAME     = 'DSUSC_CGO-XX'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM1
        PLABEL(1) = 'XXSUSCGO'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XXDSUSC',IPDSUSC(34),LEXPP,NEXPP,INDXPR,MAXEXP)
C
        PNAME     = 'DSUSC_CGO-XY'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM1
        PLABEL(1) = 'XYSUSCGO'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XYDSUSC',IPDSUSC(35),LEXPP,NEXPP,INDXPR,MAXEXP)
C
        PNAME     = 'DSUSC_CGO-XZ'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM1
        PLABEL(1) = 'XZSUSCGO'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('XZDSUSC',IPDSUSC(36),LEXPP,NEXPP,INDXPR,MAXEXP)
C
        PNAME     = 'DSUSC_CGO-YY'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM1
        PLABEL(1) = 'YYSUSCGO'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YYDSUSC',IPDSUSC(37),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'DSUSC_CGO-YZ'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM1
        PLABEL(1) = 'YZSUSCGO'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('YZDSUSC',IPDSUSC(38),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'DSUSC_CGO-ZZ'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM1
        PLABEL(1) = 'ZZSUSCGO'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
        CALL OP1IND('ZZDSUSC',IPDSUSC(39),LEXPP,NEXPP,INDXPR,MAXEXP)

      ENDIF ! of MSUSCDIA
C
C       Print section for LAO exp.values magnetizabilities
C    =======================================================
C
       IF(IPRINT.GE.1) THEN
         WRITE(LUPRI,'(2X,A)') '*** Expectation value terms '//
     &'for magnetizabilities with LAO ***'
         NTERMS_EXP=27
         DO I = 1, NTERMS_EXP
           INDXPR = LEXPP(IPDSUSC(I))
           CALL PRSYMB(LUPRI,'.',75,0)
           WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
           CALL PRSYMB(LUPRI,'.',75,0)
           CALL WRIXPR(I,INDXPR)
         ENDDO
         IF (MSUSCDIA) THEN
           WRITE(LUPRI,'(2X,A)') '*** Aucars diamagnetic terms '//
     &'for magnetizabilities with LAO ***'
           NTERMS_EXP_DIAMAG=39
           DO I = NTERMS_EXP+1,NTERMS_EXP_DIAMAG
             INDXPR = LEXPP(IPDSUSC(I))
             CALL PRSYMB(LUPRI,'.',75,0)
             WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
             CALL PRSYMB(LUPRI,'.',75,0)
             CALL WRIXPR(I,INDXPR)
           ENDDO
         ELSE
           WRITE(LUPRI,'(2X,A)') 'No Aucars diamagnetic terms '//
     &     'for magnetizabilities with LAO. Good so.'
        ENDIF

       ENDIF
      ENDIF

      CALL QEXIT('DEF_MSUSCP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck PRP_MSUSCP */
      SUBROUTINE PRP_MSUSCP(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C    PURPOSE: Output routine for magnetic susceptibilities
C
C    It is modelled after TEC's routines
C    and after the dalton routine 'abacus/abadrv.F/SUSRES'
C
C    Written by MI, jan. 2003
C    Last modifications: MI, April 2006; MI, Nov. 2007, Tel Aviv
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbprp_mod, only: IMAGOUT, IPLONDON, IPMSUSC, IPRPRP, LONDON, MSUSCDIA
!      use dcbexp_mod, only: IPDSUSC, LEXPP, NEXPP
!      use dcbxlr_mod, only: BFREQ, LLRAPU, LLRBPU, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEP
!      use dcbxpr_mod, only: PRPNAM
      use symmet_mod, only: IPTAX
!      use dcbsusc_mod, only: SUS2EL, SUSCOM, SUSDIA, SUSFS, SUSFSY, SUSREL, SUSTOT
!      use codata_mod, only: ECHARGE, EMASS, XFMOL, XTANG
#include "implicit.h"
      INTEGER :: I, II, J, K
!#include "dgroup.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "dcbxpr.h"
#include "dcbsusc.h"
#include "codata.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D0=0.0D0,
     &           THRESHM = 1.0D-5)

      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &     ABCNV(NBFREQ,NLRBPT),WORK(*)
      LOGICAL D12, D23, MAGOUT1,MAGOUT2,MAGOUT3
      DIMENSION WRK(3), IWRK(3), AXES(3,3), PVAL(6),PVAL1(6),
     &         PVAL2(6),DIASUSC(3,3),PARSUS(3,3)
      DIMENSION SUSREL_EP(3,3,2),SUSDIA_EX(3,3), EXPVAL_NS(9)
      CHARACTER*6 ABC(3)
      CHARACTER*1 JCA,JCB
      INTEGER ICONVXYZ
      DOUBLE PRECISION SI,CGS
      DATA ABC /'khi_aa','khi_bb','khi_cc'/

      CALL QENTER('PRP_MSUSCP')

      SI = (XTANG*XTANG*ECHARGE*ECHARGE/EMASS)*1.0D-21
      CGS = SI*XFMOL*1.0D6

C=============================================================================
C  Call routine for calculation of all remaining exp.value terms
C  when using London atomic orbitals.
C=============================================================================
      CALL DZERO(SUS2EL,9)
      CALL DZERO(SUSFSY,9)
      CALL DZERO(SUSFS,9)
C     ... evaluate expectation value terms first
      IF (LONDON) THEN
        CALL EXPVSUSLO(WORK(KFREE),LFREE)
      ENDIF

      KFRSAV = KFREE
      CALL TITLER('MAGNETIZABILITIES','*',124)

      IF (.NOT.XLR_SKIPEP.AND.MSUSCDIA) THEN
         IF (LONDON) THEN
          WRITE(LUPRI,'(2X,A)')
     &  'For LAO & e-p rotations the diamag.term is not (yet) defined !'
          WRITE(LUPRI,'(2X,A)')
     &  'MSUSCDIA for magnetizabilities is set to .false. ! '
          MSUSCDIA=.FALSE.
         ENDIF
      ENDIF

      CALL DZERO(SUSDIA,9)
      CALL DZERO(SUSDIA_EX,9)

      IF (.NOT.LONDON) THEN
        IF (MSUSCDIA) THEN
          IJ = 0
          DO I = 1, 3
          DO J = I, 3
            IJ = IJ + 1
            IJX = IPDSUSC(IJ)
            IF (IJ.GT.NEXPP) WRITE(LUPRI,*)
     &   "     PRP_MSUSCP: NEXPP is exceeded!"
           IF (XLR_SKIPEP) THEN
C  -1 is for the compatibility with Dalton
             SUSDIA(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJX)
             SUSDIA(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJX)
           ELSE
             SUSDIA_EX(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJX)
             SUSDIA_EX(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJX)
           ENDIF
          ENDDO
          ENDDO

          IF (XLR_SKIPEP.AND.IPRPRP.GE.1) THEN
           CALL HEADER('Expectation values part of susceptibilities',-1)
           CALL OUTPUT(SUSDIA,1,3,1,3,3,3,1,LUPRI)
           CALL FLSHFO(LUPRI)
          ENDIF
        ENDIF ! IF (MSUSCDIA) THEN
      ELSE IF (LONDON) THEN
C ... accumulate one-electron expectation values into SUSDIA
C.......................................................................................
C   IPDSUSC(1-6)...pqRM2H1, IPDSUSC(7-12)...pqRM2H2, IPDSUSC(13-18)...pqRM3H3,
C   IPDSUSC(19-27)...pqRDUSLL, IPDSUSC(28-33)...pqRDSUSNL, IPDSUSC(34-39)...pqSUSCGO
C
        CALL DZERO(SUSDIA,9)
!gosia: it seems that here it is for ipdsusc(1) - ipdsusc(18):
        IJ = 0
        NEXPVALTERMS=3 ! RM2H1 + RM2H2 + RM2H3: second magnetic derivatives of the Dirac Hamiltonian
        DO K = 1, NEXPVALTERMS
        DO I = 1, 3
        DO J = I, 3
          IJ = IJ + 1
          IJX = IPDSUSC(IJ)
          SUSDIA(IPTAX(J,2),IPTAX(I,2)) = EXPVAL(IJX) +
     &    SUSDIA(IPTAX(J,2),IPTAX(I,2))
          SUSDIA(IPTAX(I,2),IPTAX(J,2)) = SUSDIA(IPTAX(J,2),IPTAX(I,2))
        ENDDO
        ENDDO
        ENDDO

!gosia: it seems that here it is for ipdsusc(19) - ipdsusc(27):
C -------------------------------------------------------
C     RDSUSLL - first magnetic derivative of London magnetic dipole operator
C     These are 9 tensor terms, of which nondiagonal
C     have to be symmetrized according to derived formulas
C -------------------------------------------------------
C
C      ...  EXPVAL(IPDSUSC(19..27))
C
C              19  20  21
C              22  23  24
C              25  26  27
C

C  ...   so, symmetrize resulting off-diagonal RDSUSLL terms ...
        PLLXY = (EXPVAL(IPDSUSC(20))+EXPVAL(IPDSUSC(22)))/D2
        PLLXZ = (EXPVAL(IPDSUSC(21))+EXPVAL(IPDSUSC(25)))/D2
        PLLYZ = (EXPVAL(IPDSUSC(24))+EXPVAL(IPDSUSC(26)))/D2

C       ... store non-symmetrized data ...
        IIC=1
        DO II = 19, 27
         EXPVAL_NS(IIC) = EXPVAL(IPDSUSC(II))
         IIC = IIC + 1
        ENDDO

C ......... store RDSUSLL symmetrized off-diagonal terms ...
        EXPVAL(IPDSUSC(20))=PLLXY
        EXPVAL(IPDSUSC(22))=PLLXY
        EXPVAL(IPDSUSC(21))=PLLXZ
        EXPVAL(IPDSUSC(25))=PLLXZ
        EXPVAL(IPDSUSC(24))=PLLYZ
        EXPVAL(IPDSUSC(26))=PLLYZ

C    ... control print-out of both symmetrized and nosymmetrized RDSUSLL terms
        IF (IPRPRP.GE.3) THEN
         WRITE(LUPRI,'(/,10X,A)')
     &   'RDSUSLL contributions to the one-electron expectation values:'
         WRITE(LUPRI,'(A)')
     &   '-----------------------------------------------------------'//
     &   '------------'
         WRITE(LUPRI,'(4X,A)')
     &   ' Operator          nonsymmetrized (original)'//
     &   '                   symmetrized (for results)   '
         IIC = 1
         DO II = 19, 27
           INDXPR = LEXPP(IPDSUSC(II))
           WRITE(LUPRI,'(2X,A16,1X,D15.8,4X,D15.8)')
     &     PRPNAM(INDXPR),EXPVAL_NS(IIC),EXPVAL(IPDSUSC(II))
           IIC = IIC + 1
         ENDDO
        ENDIF

        ! store the RDSUSLL expectation value contributions
        DO I = 1, 3
        DO J = 1, 3
         IJ = IJ + 1
         IJX = IPDSUSC(IJ)
         SUSDIA(IPTAX(I,2),IPTAX(J,2)) = EXPVAL(IJX) +
     &   SUSDIA(IPTAX(I,2),IPTAX(J,2))
        ENDDO
        ENDDO

C ------------------------------------------------------------------------------
C  LAST 6 Aucars (Sternheim) diamagnetic terms, RDSUSNL, extracted if wanted
C ------------------------------------------------------------------------------
        IF (MSUSCDIA) THEN
          DO I = 1, 3
          DO J = I, 3
           IJ = IJ + 1
           IJX = IPDSUSC(IJ)
           SUSDIA(IPTAX(J,2),IPTAX(I,2)) =  EXPVAL(IJX) +
     &     SUSDIA(IPTAX(J,2),IPTAX(I,2))
           SUSDIA(IPTAX(I,2),IPTAX(J,2)) =
     &     SUSDIA(IPTAX(J,2),IPTAX(I,2))
C          ... save them only when e-p rotations are active!
           IF (.NOT.XLR_SKIPEP) THEN
             SUSDIA_EX(IPTAX(J,2),IPTAX(I,2)) = EXPVAL(IJX)
             SUSDIA_EX(IPTAX(I,2),IPTAX(J,2)) = EXPVAL(IJX)
           ENDIF
         ENDDO
         ENDDO

C ... This is needed for compatibility with Dalton
!gosia: but suscom is not used in following print of one-electron exp.value!
         DO I = 1, 3
         DO J = I, 3
          IJ = IJ + 1
          IJX = IPDSUSC(IJ)
          SUSCOM(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJX)
          SUSCOM(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJX)
         ENDDO
         ENDDO
       ENDIF

! ...  Write output - always print out partial results for compatibility with Dalton
       IF (IPRPRP.GE.0) THEN
        CALL HEADER('One-electron expectation values',-1)
        CALL POLPRI(SUSDIA,'   ',-2)
        AVERAG = (SUSDIA(1,1) + SUSDIA(2,2) + SUSDIA(3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F20.6)') ' Average value:',AVERAG

        CALL HEADER('Two-electron expectation values',-1)
        CALL POLPRI(SUS2EL,'   ',-2)
        AVERAG = (SUS2EL(1,1) + SUS2EL(2,2) + SUS2EL(3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F20.6)') ' Average value:',AVERAG

        CALL HEADER('Highest-order reorthonormalization',-1)
        CALL POLPRI(SUSFS,'   ',-2)
        AVERAG = (SUSFS (1,1) + SUSFS (2,2) + SUSFS (3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F20.6)') ' Average value:',AVERAG
C
        CALL HEADER('Lowest-order reorthonormalization',-1)
        CALL POLPRI(SUSFSY,'   ',-2)
        AVERAG = (SUSFSY(1,1) + SUSFSY(2,2) + SUSFSY(3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F20.6)') ' Average value:',AVERAG
        CALL FLSHFO(LUPRI)
       ENDIF
      ENDIF
C------------------------------------------------------------
C  ... get relaxation parts of susceptibilities
C     for both conventional and London atomic orbitals
C------------------------------------------------------------
      DO IFREQ = 1, NBFREQ
        CALL PRSYMB(LUPRI,'-',72,2)
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
           WRITE(LUPRI,'(A,F15.8,2X,A)')
     &           '*** Frequency    :',BFREQ(IFREQ),'a.u.'
        ENDIF
        CALL DZERO(SUSREL,9)
        CALL DZERO(SUSREL_EP,18)
C ===========================================================================
C ... extract the relaxation (response) part of magnetic susceptibilities...
C ===========================================================================
       DO IB = 1,3
         IF (LONDON) THEN
          IBX = IPLONDON(1,IB)
          READ (PRPNAM((IBX)),'(4X,A1)') JCB
          JB = ICONVXYZ(JCB)
          IF (JB.NE.1.AND.JB.NE.2.AND.JB.NE.3)  THEN
           write(lupri,*)
     &    'PRP_MSUSCP: before ICONVXYZ the name...',PRPNAM(LLRBPU(IBX))
           write(lupri,*) 'got JCB=',JCB
           CALL QUIT('PRP_MSUSCP: wrong extracted number in LAO !')
          ENDIF
         ELSE
          IBX = IPMSUSC(2,IB)
CMI  ...    extract first character, like 'X magnetic field'
          READ (PRPNAM(LLRBPU(IBX)),'(A1)') JCB
          JB = ICONVXYZ(JCB)
          IF (JB.NE.1.AND.JB.NE.2.AND.JB.NE.3)  THEN
           write(lupri,*)
     &     'PRP_MSUSCP: operator name...',PRPNAM(LLRBPU(IB))
           write(lupri,*) 'extracted JCB,JB=',JCB,JB
           CALL QUIT('PRP_MSUSCP: wrong extracted number !')
          ENDIF
         ENDIF

         KB = IPMSUSC(2,JB)
         DO IA = 1, 3
          IF (LONDON) THEN
            IAX = IPLONDON(1,IA)
            READ (PRPNAM((IAX)),'(4X,A1)') JCA
            JA = ICONVXYZ(JCA)
          ELSE
            IAX = IPMSUSC(1,IA)
            READ (PRPNAM(LLRAPU(IAX)),'(A1)') JCA
            JA = ICONVXYZ(JCA)
            IF (JA.NE.1.AND.JA.NE.2.AND.JA.NE.3)  THEN
              write(lupri,*)
     &        'PRP_MSUSCP: operator..',PRPNAM(LLRAPU(IA))
              write(lupri,*) 'extracted JCA,JA=',JCA,JA
              CALL QUIT('PRP_MSUSCP: wrong extracted number !')
            ENDIF
          ENDIF
          ! store the relaxation part
          KA = IPMSUSC(1,IA)
          SUSREL(IPTAX(IB,2),IPTAX(IA,2)) = ABLRF(IFREQ,KA,KB,1)
          IF (.NOT.XLR_SKIPEP) THEN
C          save e-e "paramagnetic" and e-p "diamagnetic"
C          contributions from the linear response function
           SUSREL_EP(IPTAX(IB,2),IPTAX(IA,2),1) = ABLRF(IFREQ,KA,KB,2)
           SUSREL_EP(IPTAX(IB,2),IPTAX(IA,2),2) = ABLRF(IFREQ,KA,KB,3)
          ENDIF
         ENDDO
       ENDDO

         CALL HEADER('Relaxation part of susceptibilities',-1)
         IF (LONDON) THEN
C    .... DALTON type of output
          CALL POLPRI(SUSREL,'   ',-2)
         ELSE
          CALL OUTPUT(SUSREL,1,3,1,3,3,3,1,LUPRI)
         ENDIF
         AVERAG = (SUSREL(1,1) + SUSREL(2,2) + SUSREL(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F20.6)') ' Average value:',AVERAG
         IF (.NOT.XLR_SKIPEP) THEN
C ....write out e-e and e-p contributions...
          CALL HEADER('e-p "diamagnetic" and e-e "paramagnetic"'//
     &      ' contributions (au):',1)
          DO IROW = 1, 3
           WRITE(LUPRI,'(2X,I1,1X,3F20.12,3X,3F20.12)')
     &       IROW,(SUSREL_EP(IROW,ICOL,2),ICOL=1,3),
     &            (SUSREL_EP(IROW,ICOL,1),ICOL=1,3)
          ENDDO
          AVERAG_EE =
     &      (SUSREL_EP(1,1,1)+SUSREL_EP(2,2,1)+SUSREL_EP(3,3,1))/D3
          AVERAG_EP =
     &      (SUSREL_EP(1,1,2)+SUSREL_EP(2,2,2)+SUSREL_EP(3,3,2))/D3
          WRITE(LUPRI,'(/,1X,A,F20.6,4X,A,F20.6)')
     &    'e-p "diamagnetic" contribution:',AVERAG_EP,
     &    'e-e "paramagnetic":',AVERAG_EE
!        ... write out MSUSCDIA diamagnetic terms for comparison
          IF (MSUSCDIA) THEN
           CALL HEADER('.MSUSCD diamagnetic term for comparison '//
     &                 'with e-p "diamagnetic" term (au):',1)
           CALL OUTPUT(SUSDIA_EX,1,3,1,3,3,3,1,LUPRI)
           AVERAG_EX =
     &      (SUSDIA_EX(1,1) + SUSDIA_EX(2,2) + SUSDIA_EX(3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F20.6)')
     &      '.MSUSCDIA average value:',AVERAG_EX
          ENDIF
C
          IF (.NOT.LONDON) THEN
C        ... copy e-e "paramagnetic" into SUSREL and e-p "diamagnetic" into SUSDIA
            CALL DCOPY(9,SUSREL_EP(1,1,1),1,SUSREL,1)
            CALL DCOPY(9,SUSREL_EP(1,1,2),1,SUSDIA,1)
          ELSE
C          ... case with London atomic orbitals...
C          ... copy the e-e "diamagnetic" contribution into SUSCOM
!gosia: line before was commented out, but then the diamagnetic contribution with LAO
!is always printed as 0
           CALL DCOPY(9,SUSREL_EP(1,1,2),1,SUSCOM,1)
          ENDIF
         ENDIF
      ENDDO

      CALL DZERO(SUSTOT, 9)
      CALL DZERO(DIASUSC,9)
      CALL DZERO(PARSUS, 9)

      IF (.NOT.LONDON) THEN
         DO J = 1, 3
         DO I = 1, 3
            DIASUSC(I,J) = SUSDIA(I,J)
            PARSUS(I,J)  = SUSREL(I,J)
            SUSTOT(I,J)  = SUSDIA(I,J) + SUSREL(I,J)
         ENDDO
         ENDDO
      ELSE
C        ... the case with London atomic orbitals ...
         DO I = 1, 3
         DO J = 1, 3
           DIASUSC(I,J) = DIASUSC(I,J) + SUSCOM(I,J)
           SUSTOT(I,J)  = SUSTOT(I,J) + SUSREL(I,J) ! response part
     &                                + SUSDIA(I,J) ! one-electron expectation values
     &                                + SUS2EL(I,J) ! two-electron expectation values
     &                                + SUSFS (I,J) ! high-order renormalization terms
     &                                + SUSFSY(I,J) ! low-order renormalization terms
           PARSUS(I,J) = SUSTOT(I,J) - DIASUSC(I,J)
         ENDDO
         ENDDO
      ENDIF

C
C     Principal values
C     ================
C
      CALL DUNIT(AXES,3)
      IJ = 1
      DO 300 I = 1, 3
         DO 310 J = 1, I
!gosia: todo check if these "magnetizability outputs 1/2/3" have sense with LAOs
!            (diasusc=0 for some keyword combination)
            PVAL1(IJ) = DIASUSC(I,J)
            PVAL2(IJ) = PARSUS(I,J)
            PVAL (IJ) = SUSTOT(I,J)
            IJ = IJ + 1
 310     CONTINUE
 300  CONTINUE
      CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
      PVAL(1) = -PVAL(1)
      PVAL(2) = -PVAL(3)
      PVAL(3) = -PVAL(6)
      CALL ORDER2(AXES,PVAL,3,3)

      D12DIFF = ABS(PVAL(1)-PVAL(2))
      D23DIFF = ABS(PVAL(2)-PVAL(3))

      D12 = (D12DIFF .LT. THRESHM)
      D23 = (D23DIFF .LT. THRESHM)

      MAGOUT1 = .FALSE.
      MAGOUT2 = .FALSE.
      MAGOUT3 = .FALSE.
      IF (IMAGOUT.EQ.0) THEN
       MAGOUT1 = D12. AND. D23
       MAGOUT2 = D12. OR.  D23
       MAGOUT3 = .TRUE.
      ELSE
       write(lupri,*)
       IF (D12.AND.D23) THEN
         write(lupri,'(3X,A)')
     &   'forced output, otherwise would be magnetizability output 1'
       ELSE IF (D12. OR.  D23) THEN
         write(lupri,'(3X,A)')
     &   'forced output, otherwise would be magnetizability output 2'
       ELSE
         write(lupri,'(3X,A)')
     &   'forced output, otherwise would be magnetizability output 3'
       ENDIF

       IF (IMAGOUT.EQ.1) MAGOUT1 = .TRUE.
       IF (IMAGOUT.EQ.2) MAGOUT2 = .TRUE.
       IF (IMAGOUT.EQ.3) MAGOUT3 = .TRUE.
      ENDIF
C
      CALL JACO(PVAL1,DUMMY,3,3,0,WRK,IWRK)
      CALL JACO(PVAL2,DUMMY,3,3,0,WRK,IWRK)
      DIAMAG = -(PVAL1(1) + PVAL1(3) + PVAL1(6))/D3
      PARAMA = -(PVAL2(1) + PVAL2(3) + PVAL2(6))/D3

      WRITE (LUPRI,'(2(/1X,A))')
     &      ' Units:  JT**(-2)/(10**(-30))   (SI) ',
     &      '         ppm cm**(3) mol**(-1)  (cgs)'

      write(lupri,'(/4X,A,F18.14,A)') '1 au=',CGS,' (cgs).'
      write(lupri,'(4X,A,F18.14,A)') '1 au=',SI*1.0D31,' (SI)'

      IF (.NOT.LONDON) THEN
         WRITE (LUPRI,'(/,1X,A)')
     &  ' Calculations without London atomic orbitals.'
      ELSE
         WRITE (LUPRI,'(/,1X,A)')
     &  ' London atomic orbitals used in calculations.'
      END IF

CMI   ...there are three defferent types of the magnetizability output
CMI  ====================================================================
      IF (MAGOUT1) THEN
         AVE = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         CALL HEADER('Magnetizabilities        '
     &             //'             au          '
     &             //'             SI          '
     &             //'            cgs',1)
         WRITE (LUPRI,'(3(1X,A,3F20.4,/))')
     &      ' Magnetizability:           ', AVE,SI*AVE*1.0D31,CGS*AVE,
     &      ' Diamagnetic contribution:  ',DIAMAG,SI*DIAMAG*1.0D31,
     &                                     CGS*DIAMAG,
     &      ' Paramagnetic contribution: ',PARAMA,SI*PARAMA*1.0D31,
     &                                     CGS*PARAMA
         WRITE (LUPRI,'(/,1X,A)') ' Magnetizability is spherical.'
      ELSE IF (MAGOUT2) THEN
         AVE = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         IF (D12) THEN
            PAR = PVAL(3)
            PER = (PVAL(1) + PVAL(2))/D2
         ELSE
            PAR = PVAL(1)
            PER = (PVAL(2) + PVAL(3))/D2
         END IF
         ANI = PAR - PER
         CALL HEADER('Magnetizabilities        '
     &             //'             au          '
     &             //'             SI          '
     &             //'            cgs',1)
         WRITE (LUPRI,'(6(1X,A,3F20.4,/))')
     &      ' Isotropic magnetizability: ',AVE,SI*AVE*1.0D31,CGS*AVE,
     &      ' Diamagnetic contribution:  ',DIAMAG,SI*DIAMAG*1.0D31,
     &                                     CGS*DIAMAG,
     &      ' Paramagnetic contribution: ',PARAMA,SI*PARAMA*1.0D31,
     &                                     CGS*PARAMA,
     &      ' Parallel component:        ',PAR,SI*PAR*1.0D31,CGS*PAR,
     &      ' Perpendicular component:   ',PER,SI*PER*1.0D31,CGS*PER,
     &      ' Anisotropy:                ',ANI,SI*ANI*1.0D31,CGS*ANI
         WRITE (LUPRI,'(1X,A)') ' Magnetizability is cylindrical.'
         CALL HEADER('Principal values (au, SI, and cgs) and axes:',1)
         DO 700 I = 1, 3
            WRITE (LUPRI,'(2X,A,1X,3F12.4,2X,3F10.4)')
     &         ABC(I),PVAL(I),SI*PVAL(I)*1.0D31,CGS*PVAL(I),
     &         (AXES(IPTAX(J,2),I),J=1,3)
 700     CONTINUE
      ELSE IF (MAGOUT3) THEN
         CALL DSWAP(1,PVAL(2),1,PVAL(1),1)
         CALL DSWAP(3,AXES(1,2),1,AXES(1,1),1)
         AVE  = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         ANI1 = PVAL(3) - (PVAL(1) + PVAL(2))/D2
         ANI2 = PVAL(2) - (PVAL(1) + PVAL(3))/D2
         CALL HEADER('Magnetizabilities        '
     &             //'             au          '
     &             //'             SI          '
     &             //'            cgs',1)
         WRITE (LUPRI,'(5(1X,A,3F20.4,/))')
     &      ' Isotropic magnetizability: ',AVE, SI*AVE*1.0D31, CGS*AVE,
     &      ' Diamagnetic contribution:  ',DIAMAG,SI*DIAMAG*1.0D31,
     &                                     CGS*DIAMAG,
     &      ' Paramagnetic contribution: ',PARAMA,SI*PARAMA*1.0D31,
     &                                     CGS*PARAMA,
     &      ' 1st anisotropy:            ',ANI1,SI*ANI1*1.0D31,CGS*ANI1,
     &      ' 2nd anisotropy:            ',ANI2,SI*ANI2*1.0D31,CGS*ANI2
         CALL HEADER('Principal values (au, SI, and cgs) and axes:',1)
         DO 800 I = 1, 3
            WRITE (LUPRI,'(2X,A,1X,3F12.4,2X,3F10.4)')
     &         ABC(I),PVAL(I),SI*PVAL(I)*1.0D31,CGS*PVAL(I),
     &         (AXES(IPTAX(J,2),I),J=1,3)
 800     CONTINUE
      ELSE
       write(lupri,*) 'PRP_MSUSCP: blind branch of the output !!!'
      END IF
C
      CALL HEADER('Total magnetizability tensor (au)',-1)
csonia 04/10/1995
csonia 04/10/1995
      CALL DSCAL(9,-D1,SUSTOT,1)
      CALL POLPRI(SUSTOT,'   ',-2)
CMI   CALL OUTPUT(SUSTOT,1,3,1,3,3,3,1,LUPRI)
C     CALL DSCAL(9,-D1,SUSTOT,1)
C
Cspas: also print diamagnetic and paramagnetic tensors
C
      CALL HEADER('Diamagnetic magnetizability tensor (au)',-1)
      CALL DSCAL(9,-D1,DIASUSC,1)
CMI   CALL OUTPUT(DIASUSC,1,3,1,3,3,3,1,LUPRI)
      CALL POLPRI(DIASUSC,'   ',-2)
C     CALL DSCAL(9,-D1,DIASUSC,1)

      CALL HEADER('Paramagnetic magnetizability tensor (au)',-1)
      CALL DSCAL(9,-D1,PARSUS,1)
CMI   CALL OUTPUT(PARSUS,1,3,1,3,3,3,1,LUPRI)
      CALL POLPRI(PARSUS,'   ',-2)
C     CALL POLPRI(PARSUS,'PRI',-2)
C     CALL DSCAL(9,-D1,PARSUS,1)

      CALL FLSHFO(LUPRI)
      CALL QEXIT('PRP_MSUSCP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_SHIELD */
      SUBROUTINE PRP_SHIELD(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR shieldings
C     This routine is modelled after PRP_SPNSPN
C
C     Written by T.Enevoldsen Nov. 7 1997
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbxlr_mod, only: NBFREQ, NLRAPT, NLRBPT
      use nuclei_mod, only: NUCDEP, NUCIND
#include "implicit.h"
C
!#include "dcbprp.h"
!#include "dcbexp.h"
#include "dcbxlr.h"
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &     ABCNV(NBFREQ,NLRBPT),WORK(*)
      logical, allocatable :: DOATOM(:)
      KFRSAV = KFREE
      NCMAT =  NUCDEP*9
      CALL MEMGET('REAL',KTMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATLO,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATRE,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATREE,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATREP,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATREO,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT,2*NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KANIS,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KASYM,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSKEW,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSPAN,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMPMAT,2*NCMAT,WORK,KFREE,LFREE)
C
C     (SETATM code -29: accept .SELECT specification)
      allocate(DOATOM(NUCIND))
      CALL SETATM(DOATOM,NATOM,-29)
      CALL PRP_SHIELD1(EXPVAL,ABLRF,WORK(KTMPMAT),WORK(KTMAT),
     &     WORK(KDMATLO),WORK(KDMATRE),WORK(KDMATREE),WORK(KDMATREP),
     &     WORK(KDMATREO),
     &     WORK(KDMAT),WORK(KPMAT),WORK(KTISO),WORK(KDISO),WORK(KPISO),
     &     WORK(KANIS),WORK(KASYM),WORK(KSKEW),WORK(KSPAN),WORK(KIBUF),
     &     DOATOM,WORK(KCSTRA),WORK(KSCTRA),WORK,KFREE,LFREE)
      deallocate(DOATOM)
      CALL MEMREL('PRP_SHIELD',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_SHIELD1 */
      SUBROUTINE PRP_SHIELD1(EXPVAL,ABLRF,TMPMAT,TMAT,
     &      DMATLO,DMATRE,DMATREE,DMATREP,DMATREO,DMAT,PMAT,TISO,
     &      DISO,PISO,ANIS,ASYM,SKEW,SPAN,IBUFA,
     &      DOATOM,CSTRA,SCTRA,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR shieldings
C
C     Written by T.Enevoldsen Nov. 7 1997
C
C     Updated by MI,2003,2004 for LAO based calculation
C     of NMR shieldings.
C
C***********************************************************************
!      use dcbprp_mod, only: IPNSTDIA, IPRNST, IPRPRP, IPSHIELD, LONDON, NNSTIN, NOORTH, NSTDIA, SEPTEP
!      use dcborb_mod, only: N2ORBXQ
!      use dcbxpr_mod, only: IPRPLBL, PRPNAM
!      use dcbprl_mod, only: PRPLBL
!      use dcbexp_mod, only: LEXPP
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEE, XLR_SKIPEP
!      use dcbnmr_mod, only: EPREORTH
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
      use memory_allocator
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use symmet_mod, only: IPTAX, MAXOPR
      use nuclei_mod, only: ISTBNU, NAMDEP, NAMDPX, NUCDEP, NUCIND
#include "implicit.h"
      INTEGER :: I
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "dcbnmr.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, D2I = 0.50D0)
      LOGICAL TEST,DOATOM(NUCIND), ee, ep,occ
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &     TMPMAT(3,3*NUCDEP,2),TMAT(3,3,NUCDEP), DMAT(3,3,NUCDEP),
     &     DMATLO(3,3,NUCDEP),DMATRE(3,3,NUCDEP),
     &     DMATREE(3,3,NUCDEP),DMATREP(3,3,NUCDEP),DMATREO(3,3,NUCDEP),
     &     PMAT(3,3,NUCDEP,2),TISO(NUCDEP),DISO(NUCDEP),PISO(NUCDEP),
     &     ANIS(NUCDEP),ASYM(NUCDEP),SKEW(NUCDEP),SPAN(NUCDEP),
     &     IBUFA(3*NUCDEP),IBUFB(3),CSTRA(NUCDEP*NUCDEP*9),
     &     SCTRA(NUCDEP*NUCDEP*9),WORK(*)
      INTEGER ICONVXYZ
      CHARACTER*1 CJB
C
      CALL QENTER('PRP_SHIELD1')
C
      CALL TITLER('NMR SHIELDING TENSOR','*',116)

#ifdef MOD_LAO_REARRANGED
! gosia: if .SHIEL2 then calculate NMR shielding + LAO from rearranged
! equation (Eq. 40-42 in JCP 131, 124119 (2009))
      if (shielding_rearrange) then
        write(*, *) 'NMR shielding tensor calculated ',
     &               'from rearranged equations ',
     &               '(Eq.40-42 in JCP 131, 124119 (2009))'
      end if
! if .LAOMOD, then NMR shielding is calculated from rearranged equations
! but also LR equations are modified
      if (lao_lr_rearrange) then
        write(*, *) 'NMR shielding rearranged, LR equations in LAO ',
     &               'basis are rearranged too (fully connection-',
     &               'independent approach)'
      end if
#endif

C
      KFRSAV = KFREE
      NSCOOR = 3*NUCDEP
      NCMAT  = 9*NUCDEP
      CALL DZERO(TMPMAT,NCMAT)
      CALL DZERO(DMAT,NCMAT)
      CALL DZERO(DMATLO,NCMAT)
      CALL DZERO(DMATRE,NCMAT)
      CALL DZERO(DMATREE,NCMAT)
      CALL DZERO(DMATREP,NCMAT)
      CALL DZERO(TMAT,NCMAT)
      CALL DZERO(PMAT,2*NCMAT)
C--------------------------------------------------------------------------
C     NSTDIA: Nuclear Shielding Tensor DIAmagnetic contribution
C     Extract Diamagnetic approximation if calculated
C     (used for total shielding if e-p rotations are not included)
C       -- this becomes the
C     A**2 diamagnetic operator in the non-relativistic limit
C     (ref: Aucar et al., JCP 110 (1999) 6208) /hjaaj Mar 2003
C
      IF(NSTDIA) THEN
        DO I = 1,NNSTIN
          IOP = IPNSTDIA(I)
          IPL = IPRPLBL(1,LEXPP(IOP))
CMI       ... take care of the correct form of PRPLBL(IPL) name
          READ(PRPLBL(IPL)(1:3),'(I3)') IA
          IB = MOD(I+2,3) + 1
          TMPMAT(IPTAX(IB,2),IA,1) = EXPVAL(IOP)
        ENDDO
        IF(IPRPRP.GE.2) THEN
C        ...for what do we have calculated diamagnetic terms
          IF (XLR_SKIPEP) THEN
            CALL HEADER('Diamagnetic part of shielding'//
     &      ' in symmetry coordinates (ppm)',1)
          ELSE
            CALL HEADER('Diamagnetic terms for comparison'//
     &      ' with e-p "diamagnetic" in symmetry coordinates (ppm)',1)
          ENDIF
          CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
        ENDIF
C       ... extract the diamagnetic part into DMAT
        CALL TRADIP(TMPMAT,DMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
      END IF ! IF(NSTDIA) THEN
C ------------------------------------------------------
C  Get the London orbital contributions to the exp.value
C  ... these are RNST terms
C ------------------------------------------------------
      IF(LONDON) THEN
        CALL DZERO(TMPMAT,NCMAT)
        CALL DZERO(PMAT  ,NCMAT)
        DO I = 1,NNSTIN
          IOP = IPRNST(I)
          IPL = LEXPP(IOP)
CMI ... take care of the correct form of PRPLBL(IPL) name
          READ(PRPNAM(IPL)(10:11),'(I2)') IA
          IB = MOD(I+2,3) + 1
          TMPMAT(IPTAX(IB,2),IA,1) = -EXPVAL(IOP)
        ENDDO
        IF (IPRPRP.GE.2) THEN
          CALL HEADER('London part of shielding (exp.value)'//
     &               ' in symmetry coordinates (ppm)',1)
          CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
        ENDIF
C       ... extract London part of shielding into the PMAT
        CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C       ... add the London exp.value parts of shield. to DMATLO
CTROND: DCOPY is not needed, use DMATLO directly
        CALL DAXPY(NCMAT,D1,PMAT,1,DMATLO,1)
      END IF
C--------------------------------------------------------------------
C  Get the (expectation value) reorthonormalization terms, {T,h(mK)}
C--------------------------------------------------------------------
      IF(LONDON.AND.(.NOT.NOORTH)) THEN

#ifdef MOD_LAO_REARRANGED
        if (SEPTEP) then

         if (.not. XLR_SKIPEE) then
! occ-occ:
           ee = .false.
           ep = .false.
           occ = .true.
           CALL DZERO(TMPMAT,NCMAT)
           CALL DZERO(PMAT,NCMAT)
           CALL GETCON_septep(ee,ep,occ,TMPMAT,WORK,KFREE,LFREE,IPRPRP)
           CALL DSCAL(NCMAT,DM1,TMPMAT,1)
           IF (IPRPRP.GE.2) THEN
             CALL HEADER('reorthonormalization contribution '//
     &            'from occ-occ ("ii") part of T matrix',1)
             CALL FCPRI (TMPMAT,'SIGMAS',CSTRA,SCTRA)
           ENDIF
           CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
           CALL DAXPY(NCMAT,D1,PMAT,1,DMATREO,1)
! e-occ:
           ee = .true.
           ep = .false.
           occ = .false.
           CALL DZERO(TMPMAT,NCMAT)
           CALL DZERO(PMAT,NCMAT)
           CALL GETCON_septep(ee,ep,occ,TMPMAT,WORK,KFREE,LFREE,IPRPRP)
           CALL DSCAL(NCMAT,DM1,TMPMAT,1)
           IF (IPRPRP.GE.2) THEN
             CALL HEADER('reorthonormalization contribution from '//
     &      'occ-"secondary electronic" ("ia+") part of T matrix',1)
             CALL FCPRI (TMPMAT,'SIGMAS',CSTRA,SCTRA)
           ENDIF
           CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
           CALL DAXPY(NCMAT,D1,PMAT,1,DMATREE,1)

         end if
         if (EPREORTH) then
! p-occ:
           ee = .false.
           ep = .true.
           occ = .false.
           CALL DZERO(TMPMAT,NCMAT)
           CALL DZERO(PMAT,NCMAT)
           CALL GETCON_septep(ee,ep,occ,TMPMAT,WORK,KFREE,LFREE,IPRPRP)
           CALL DSCAL(NCMAT,DM1,TMPMAT,1)
           IF (IPRPRP.GE.2) THEN
             CALL HEADER('reorthonormalization contribution from '//
     &      'occ-"secondary positronic" ("ia-") part of T matrix',1)
             CALL FCPRI (TMPMAT,'SIGMAS',CSTRA,SCTRA)
           ENDIF
           CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
           CALL DAXPY(NCMAT,D1,PMAT,1,DMATREP,1)

         end if

        end if
#endif

        CALL DZERO(TMPMAT,NCMAT)
        CALL DZERO(PMAT,NCMAT)
C  ...  get the expect. value {T,h(mK)} reorthonormalization term
!  ...  gosia: if connection-independent formulation of shieldings,
!  ...  then the -{S^B, h(mK)} expectation value is calculated here
!todo from Trond's email:
!check the use of GETCON in  PRP_SHIELD1(pamprp.F). Can it be completely removed ?
! gosia: no, it's one of contributions to NMR shielding
        CALL MEMGET('REAL',KW   ,3*N2ORBXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KCON ,N2ORBXQ  ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KDMAT,N2ORBXQ  ,WORK,KFREE,LFREE)
        CALL GETCON(WORK(KW),WORK(KCON),WORK(KDMAT),TMPMAT,
     &              WORK,KFREE,LFREE,IPRPRP)
        CALL MEMREL('PRP_SHIELD1_GETCON',WORK,1,KW,KFREE,LFREE)
CMI ... sign convention !!!
        CALL DSCAL(NCMAT,DM1,TMPMAT,1)
        IF (IPRPRP.GE.2) THEN
#ifdef MOD_LAO_REARRANGED
          if (lao_lr_rearrange .or. shielding_rearrange) then
            CALL HEADER('Contribution to the shielding from overlap '//
     &      'magnetic derivative (expect.value) in symmetry '//
     &      'coordinates (ppm)',1)
            CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
          else
#endif
            CALL HEADER('Reorthonormalization part of the '//
     &      'shielding (expect.value) in symmetry coordinates (ppm)',1)
            CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
#ifdef MOD_LAO_REARRANGED
          end if
#endif
        ENDIF
C       ... extract the reorthonormalization part into the PMAT
        CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C       ... add the reorthonormalization terms to DMATRE
CTROND...DAXPY not needed !!!! use dmatre directly
        CALL DAXPY(NCMAT,D1,PMAT,1,DMATRE,1)
      END IF ! (LONDON.AND.(.NOT.NOORTH))
C
C     *** Extract shieldings ***
C
C     The NSCOOR nucl. magn. dipole operators as A operators in <<A;B>> :
      NA = 0
      DO I = 1,NSCOOR
        IF(IPSHIELD(I).GT.0) THEN
          NA = NA + 1
          IBUFA(NA) = I
        ENDIF
      END DO
CMI/HJ internal consistency check
      IF (NA.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_SHIELD1: NA=',NA
        CALL QUIT('PRP_SHIELD1:Wrong value of NA!')
      ENDIF
C     The 3 B field directions as B operators in <<A;B>> :
      NB = 0
      DO I = 1, 3
        IF(IPSHIELD(I+MXCOOR).GT.0) THEN
          NB = NB + 1
          IBUFB(NB) = I
        ENDIF
      END DO
      IF (NB.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_SHIELD1: NB=',NB
        CALL QUIT('PRP_SHIELD1:Wrong value of NB!')
      ENDIF
C
C     Fill shielding matrices
C
      IF (XLR_SKIPEE) THEN
        IABLRF_EE = -999999999
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 2
        END IF
      ELSE
        IABLRF_EE = 2
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 3
        END IF
      END IF
#if defined (XH4)
      WRITE(LUPRI,'(A)') '*** WARNING *** Using code for XH4 !'
#endif
      DO IFREQ = 1,NBFREQ
        CALL PRSYMB(LUPRI,'-',72,2)
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
           WRITE(LUPRI,'(A,F15.8,A)')
     &          '*** Frequency    :',BFREQ(IFREQ),' a.u.'
        ENDIF
        CALL DZERO(TMPMAT,2*NCMAT)
        DO IB = 1,NB
          JB = IBUFB(IB)
          KB = IPSHIELD(MXCOOR+JB)
          JB = IPTAX(JB,2)
          DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPSHIELD(JA)
C...        First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         TMPMAT(JB,JA,1) = ABLRF(IFREQ,KA,KB,IABLRF_EE)
C...        then extract the e-p "diamagnetic" contribution
            IF (.NOT. XLR_SKIPEP)
     &         TMPMAT(JB,JA,2) = ABLRF(IFREQ,KA,KB,IABLRF_EP)
          ENDDO
        ENDDO
        IF (LONDON) THEN
        ! ... different sign convention for London operators
           CALL DSCAL(2*NCMAT,DM1,TMPMAT,1)
           IF (IPRPRP .GT. 5) THEN
              write(lupri,*) 'Changing sign on London contributions e-e'
              call output(TMPMAT,1,NB,1,NA,NB,NA,-1,LUPRI)
              write(lupri,*) 'Changing sign on London contributions e-p'
              call output(TMPMAT(1,1,2),1,NB,1,NA,NB,NA,-1,LUPRI)
           END IF
        END IF
#if defined (XH4)
        TMPMAT(2,8,1)  = TMPMAT(1,4,1)
        TMPMAT(2,9,1)  = TMPMAT(1,6,1)
        TMPMAT(2,10,1) = TMPMAT(1,7,1)
        TMPMAT(2,11,1) = TMPMAT(1,5,1)
        TMPMAT(3,12,1) = TMPMAT(1,4,1)
        TMPMAT(3,13,1) = TMPMAT(1,6,1)
        TMPMAT(3,15,1) = TMPMAT(1,7,1)
        TMPMAT(3,14,1) = TMPMAT(1,5,1)
        TMPMAT(2,8,2)  = TMPMAT(1,4,2)
        TMPMAT(2,9,2)  = TMPMAT(1,6,2)
        TMPMAT(2,10,2) = TMPMAT(1,7,2)
        TMPMAT(2,11,2) = TMPMAT(1,5,2)
        TMPMAT(3,12,2) = TMPMAT(1,4,2)
        TMPMAT(3,13,2) = TMPMAT(1,6,2)
        TMPMAT(3,15,2) = TMPMAT(1,7,2)
        TMPMAT(3,14,2) = TMPMAT(1,5,2)
#endif
C
C     Transform to C1 basis
C
        CALL DZERO(PMAT,2*NCMAT)
        CALL DZERO(TMAT,NCMAT)
        IF (.NOT. XLR_SKIPEE) THEN
          IF (IPRPRP.GE.2) THEN
            CALL HEADER('e-e response part of shielding'//
     &           ' in symmetry coordinates (ppm)',1)
            CALL FCPRI (TMPMAT,'SIGMAS',CSTRA,SCTRA)
          ENDIF
C         .... extract the e-e response part into PMAT(...,1)
          CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C         .... save e-e response part into TMAT
          CALL DAXPY(9*NUCDEP,D1,PMAT(1,1,1,1),1,TMAT,1)
        END IF
C       ... include London + reorthonormalization to total shielding
        IF (LONDON) THEN
          CALL DAXPY(9*NUCDEP,D1,DMATLO,1,TMAT,1)
          CALL DAXPY(9*NUCDEP,D1,DMATRE,1,TMAT,1)
        ENDIF
        IF (XLR_SKIPEP) THEN
C       ... add the diamagnetic term(s),DMAT, to TMAT(contains paramagnetic terms)
            CALL DAXPY(9*NUCDEP,D1,DMAT,1,TMAT,1)
        ELSE
C       ... add (diamagnetic-NSNL)+London+reortho terms to the TMAT
          IF (IPRPRP.GE.2) THEN
             CALL HEADER('e-p response part of shielding'//
     &          ' in symmetry coordinates (ppm)',1)
             CALL FCPRI (TMPMAT(1,1,2),'SIGMAS',CSTRA,SCTRA)
          ENDIF
C         ... extract the e-p response part into PMAT(...,2)
          CALL TRADIP(TMPMAT(1,1,2),PMAT(1,1,1,2),
     &                CSTRA,SCTRA,3*NUCDEP,2,2)
C         ... add the e-p response PMAT(...,2) into TMAT
C         so the TMAT now contains total response
          CALL DAXPY(9*NUCDEP,D1,PMAT(1,1,1,2),1,TMAT,1)
        END IF
C
C     Analyze results
C
        IATOM = 0
        DO 100 I = 1, NUCIND
          DO 200 ISYMOP = 0, MAXOPR
          IF (IAND(ISTBNU(I),ISYMOP).EQ.0) THEN
            IATOM = IATOM + 1
            IF (DOATOM(I)) THEN
              CALL SHIANA(TMAT(1,1,IATOM),DMAT(1,1,IATOM),
     &             DMATLO(1,1,IATOM),DMATRE(1,1,IATOM),
     &             DMATREE(1,1,IATOM), DMATREP(1,1,IATOM),
     &             DMATREO(1,1,IATOM),
     &             PMAT(1,1,IATOM,1),PMAT(1,1,IATOM,2),
     &             TISO(IATOM),DISO(IATOM),PISO(IATOM),
     &             ANIS(IATOM),ASYM(IATOM),SKEW(IATOM),
     &             SPAN(IATOM),NAMDEP(IATOM),NAMDPX(3*(IATOM-1)+1),
     &             IPRPRP)
            END IF
          END IF
 200      CONTINUE
 100    CONTINUE
C
        CALL AROUND('Summary of chemical shieldings')
        WRITE (LUPRI,'(A)') ' Definitions from J.Mason, Solid '//
     &      'state Nuc.Magn.Res. 2 (1993), 285'
        WRITE (LUPRI,'()')
        WRITE (LUPRI,'(A/A/A/A)')
     &      '@            isotropic shielding',
     &      '@        ----------------------------',
     &      '@atom          total          dia         para'//
     &      '            skew'//
     &      '            span            anis            asym',
     &      '@-----------------------------------------------------'//
     &      '-----------------------'

        IATOM = 0
        DO 300 I = 1, NUCIND
          DO 400 ISYMOP = 0, MAXOPR
          IF (IAND(ISTBNU(I),ISYMOP).EQ.0) THEN
            IATOM = IATOM + 1
            IF (DOATOM(I)) THEN
              WRITE (LUPRI,'(A,A,7F14.4)') '@',
     &        NAMDEP(IATOM),TISO(IATOM),DISO(IATOM),PISO(IATOM),
     &        SKEW(IATOM),SPAN(IATOM),ANIS(IATOM),ASYM(IATOM)
            END IF
          END IF
 400      CONTINUE
 300    CONTINUE
        WRITE(LUPRI,'(A,/)')
     &        '@-----------------------------------------------------'//
     &        '-----------------------'
C
      END DO
      CALL QEXIT('PRP_SHIELD1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck SHIANA */
      SUBROUTINE SHIANA(TMAT,DMAT,DMATLO,DMATRE,DMATREE,DMATREP,
     &                  DMATREO,PMATEE,PMATEP,
     &                  AVETOT,AVEDIA,AVEPAR,ANIS,ASYM,
     &                  SKEW,SPAN,NAME,NAMEX,IPRINT)
C***********************************************************************
C
C     Analyze routine for NMR shieldings
C     This routine is modelled after SHIANA in DALTON
C
C     INPUT:
C        TMAT   - total shielding
C        DMAT   - diamagnetic term if calculated as expectation value
C        DMATLO - London term
C        DMATRE - reorthonormalization term
C        PMATEE - e-e (++) - contribution to response function
C        PMATEP - e-p (+-) - contribution to response function
C     Written by T.Enevoldsen Nov. 14 1997
C
CMI/Dec '04  DMATLO is reserved for future use
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: LONDON, NSTDIA, SEPTEP
!      use dcbxlr_mod, only: THCLR, XLR_SKIPEP
!      use codata_mod, only: ALPHA2
#include "implicit.h"
      INTEGER :: I, J, K
#include "dcbprp.h"
#include "dcbxlr.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      CHARACTER NAME*6, NAMEX(3)*6, SGN*3
      INTEGER Z
      DOUBLE PRECISION FACTOR
      DIMENSION TMAT(3,3), DMAT(3,3), PMATEE(3,3), PMATEP(3,3),
     &          SYM(3,3),  ANTI(3,3), AXES(3,3), DMATLO(3,3),
     &          DMATRE(3,3),DMATREE(3,3),DMATREP(3,3),DMATREO(3,3),
     &          PVAL(6), WRK(3), IWRK(3), ABSP(3)
C
      FACTOR = 1.D6*ALPHA2
C
      CALL DSCAL(9,FACTOR,TMAT,1)
      CALL DSCAL(9,FACTOR,DMAT,1)
      CALL DSCAL(9,FACTOR,DMATLO,1)
      CALL DSCAL(9,FACTOR,DMATRE,1)
      CALL DSCAL(9,FACTOR,DMATREE,1)
      CALL DSCAL(9,FACTOR,DMATREP,1)
      CALL DSCAL(9,FACTOR,DMATREO,1)
      CALL DSCAL(9,FACTOR,PMATEE,1)
      CALL DSCAL(9,FACTOR,PMATEP,1)
C     Total shielding
      AVETOT = (TMAT(1,1) + TMAT(2,2) + TMAT(3,3))/D3
C     Diamagnetic part
      IF (XLR_SKIPEP) THEN
C     ... use NSTDIA expectation value approximation
         AVEDIA = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
      ELSE
C     ... use result of e-p response equations
        AVEDIA  = (PMATEP(1,1) + PMATEP(2,2) + PMATEP(3,3))/D3
        IF(NSTDIA) THEN
          AVEDIA2 = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
        ENDIF
      END IF
      IF(LONDON) THEN
C       ...add London + reorthonormalization contribution
          AVEDIA1 = (DMATLO(1,1) + DMATLO(2,2) + DMATLO(3,3))/D3
          AVEDIA3 = (DMATRE(1,1) + DMATRE(2,2) + DMATRE(3,3))/D3
          AVEDIA3E = (DMATREE(1,1) + DMATREE(2,2) + DMATREE(3,3))/D3
          AVEDIA3P = (DMATREP(1,1) + DMATREP(2,2) + DMATREP(3,3))/D3
          AVEDIA3O = (DMATREO(1,1) + DMATREO(2,2) + DMATREO(3,3))/D3
          AVEDIA = AVEDIA + AVEDIA1 + AVEDIA3
      ENDIF
C     Paramagnetic part
      AVEPAR = (PMATEE(1,1) + PMATEE(2,2) + PMATEE(3,3))/D3
C
C     Parameters calculated in accordance with J.Mason, Solid State
C     Nuc.Magn.Resonance, 2 (1993), 285
C
      DO 100 I = 1, 3
         DO 200 J = 1, 3
            ANTI(I,J) = (TMAT(I,J) - TMAT(J,I))/D2
            SYM (I,J) = (TMAT(I,J) + TMAT(J,I))/D2
  200    CONTINUE
  100 CONTINUE
      CALL DUNIT(AXES,3)
      IJ = 1
      DO 300 I = 1, 3
         DO 310 J = 1, I
            PVAL(IJ) = SYM(I,J)
            IJ = IJ + 1
  310    CONTINUE
  300 CONTINUE
      CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
      PVAL(1) = PVAL(1)
      PVAL(2) = PVAL(3)
      PVAL(3) = PVAL(6)
C
      IZ = IDMAX(3,PVAL,1)
      IX = IDMIN(3,PVAL,1)
      DO 320 I = 1, 3
        IF (I .NE. IZ .AND. I .NE. IX) IY = I
 320  CONTINUE
C
      SPAN  = PVAL(IZ) - PVAL(IX)
      SKEW  = D3*(AVETOT - PVAL(IY))/SPAN
      ANIS = PVAL(IZ) - (PVAL(IX) + PVAL(IY))/D2
      ASYM = (PVAL(IY) - PVAL(IX))/(PVAL(IZ) - AVETOT)
C This is done to avoid numerical instabilities
      IF (ABS(SPAN).LT.THCLR) THEN
        SPAN = D0
        SKEW = D0
        ANIS = D0
        ASYM = D1
      ELSE
        IF (ABS(PVAL(IY)-PVAL(IX)).LT.THCLR) THEN
          SKEW = D1
        END IF
        IF (ABS(PVAL(IZ)-PVAL(IY)).LT.THCLR) THEN
          SKEW = -D1
        END IF
      END IF
      IF (ABS((PVAL(IZ) - AVETOT)).LT.THCLR) ASYM = D1
C
      IF (IPRINT.GE.1) THEN
        CALL TITLER('Chemical shielding for '//NAME//':','=',1)
        WRITE (LUPRI,'(5(/,1X,A,F12.4,A))')
     &   ' Shielding constant:',AVETOT,' ppm',
     &   ' Anisotropy:        ',ANIS,  ' ppm',
     &   ' Span:              ',SPAN,  ' ppm',
     &   ' Asymmetry:         ',ASYM,  '    ',
     &   ' Skew:              ',SKEW,  '    '
        CALL HEADER('Total shielding tensor (ppm):',1)
        WRITE (LUPRI,'(18X,3(A,13X),/)') 'Bx', 'By', 'Bz'
        DO ICOOR = 1,3
          WRITE (LUPRI,'(2X,A6,3F15.8)')
     &          NAMEX(ICOOR), (TMAT(K,ICOOR),K=1,3)
        ENDDO
C
        IF(XLR_SKIPEP) THEN
C   XLR_SKIPEP: Diamagnetic contribution (Sternheimer approx.)
          CALL HEADER
     &    ('Diamagnetic and paramagnetic contributions (ppm):',1)
          WRITE (LUPRI,'(15X,3(A,9X),4X,3(A,9X),/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
          DO ICOOR = 1, 3
            WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &         (DMAT(K,ICOOR),K=1,3), (PMATEE(K,ICOOR),K=1,3)
          END DO
          WRITE (LUPRI,'(/A,F14.6,8X,A,F14.6)')
     &          ' Diamagnetic contribution:',AVEDIA,
     &          ' Paramagnetic: ',AVEPAR
        ELSE
C   ... make separate print out of e-e and e-p response contributions
          CALL HEADER('e-p "diamagnetic" and e-e "paramagnetic"'//
     &      ' response contributions (ppm):',1)
          WRITE (LUPRI,'(15X,3(A,9X),4X,3(A,9X)/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
          DO ICOOR = 1, 3
            WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &            (PMATEP(K,ICOOR),K=1,3), (PMATEE(K,ICOOR),K=1,3)
          END DO
!gosia: below we want to print isotropic value calculated for PMATEP tensor printed above
!so if we use london, we cannot use AVEDIA, as this already contains other terms than response,
!and those are printed separately
          WRITE (LUPRI,'(/A,F14.6,8X,A,F14.6)')
     &          ' e-p "diamagnetic" contribution:',
     &          (PMATEP(1,1) + PMATEP(2,2) + PMATEP(3,3))/D3,
     &          ' e-e "paramagnetic": ',AVEPAR
          IF (NSTDIA) THEN
            IF (LONDON) THEN
              CALL HEADER('.NSTDIA London diamagnetic term '//
     &       'for comparison  with e-p "diamagnetic" term (ppm):',1)
            ELSE
              CALL HEADER('.NSTDIA CGO diamagnetic term '//
     &        ' for comparison with e-p "diamagnetic" term (ppm):',1)
            END IF
            WRITE (LUPRI,'(15X,3(A,9X)/)') 'Bx', 'By', 'Bz'
            DO ICOOR = 1, 3
              WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &             (DMAT(K,ICOOR),K=1,3)
            END DO
            WRITE (LUPRI,'(/A,F14.6)')
     &          ' .NSTDIA diamagnetic contribution, for comparison:',
     &          AVEDIA2
          ENDIF
        ENDIF
        IF(LONDON) THEN
          CALL HEADER('London contributions (ppm):',1)
          DO ICOOR = 1, 3
            WRITE (LUPRI,'(2X,A6,3F11.4)') NAMEX(ICOOR),
     &            (DMATLO(K,ICOOR),K=1,3)
          END DO
          WRITE (LUPRI,'(/A,F14.6)') ' London contribution:',AVEDIA1
          CALL HEADER('Reorthonormalization contributions (ppm):',1)
          DO ICOOR = 1, 3
            WRITE (LUPRI,'(2X,A6,3F11.4)') NAMEX(ICOOR),
     &            (DMATRE(K,ICOOR),K=1,3)
          END DO
          WRITE (LUPRI,'(/A,F14.6)')
     &    ' Reorthonormalization contribution:',AVEDIA3
          if (SEPTEP) then
           write(LUPRI,'(/A60,F14.6)') '...from "occ-secondary '//
     &            'electronic (ia+)" blocks of T = ',
     &                                AVEDIA3E
           write(LUPRI,'(/A60,F14.6)') '...from "occ-secondary '//
     &            'positronic (ia-)" blocks of T = ',
     &                                AVEDIA3P
           write(LUPRI,'(/A60,F14.6)') '...from "occ-occ (ij)" '//
     &            'blocks of T = ',
     &                                AVEDIA3O
          end if
        ENDIF
C
        CALL HEADER
     &   ('Antisymmetric and traceless symmetric parts (ppm):',1)
        WRITE (LUPRI,'(15X,3(A,9X),4X,3(A,9X),/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
        DO 600 ICOOR = 1, 3
          WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &         (ANTI(K,ICOOR),K=1,3), (SYM(K,ICOOR),K=1,3)
  600   CONTINUE
        CALL HEADER('Principal values and axes:',1)
        DO 700 I = 1, 3
          IF (PVAL(I) .GE. D0) THEN
            SGN = '  +'
          ELSE
            SGN = '  -'
          END IF
          WRITE (LUPRI,'(2X,A,I1,2X,F12.6,A,F8.2,A,F7.2,A,3X,3F10.6)')
     &      NAME,I,PVAL(I) + AVETOT,'  =',AVETOT,SGN,
     &      ABS(PVAL(I)),':', (AXES(J,I),J=1,3)
  700   CONTINUE
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PRP_PVC_SHIELD */
      SUBROUTINE PRP_PVC_SHIELD(ABLRF,ABCNV)
C***********************************************************************
C
C     Output routine for Parity-violation contribution
C     to NMR shielding constants
C
C     Radovan Bast                         Revision:      september 2006
C     Agustin Aucar                        Last revision: november  2022
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: IPRPRP, IP_PVC_SHIELD
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEE, XLR_SKIPEP
      use nuclei_mod, only: NATOMS, NFLOAT, NUCDEP, NUCIND
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: IPTAX
      use consts_mod, only: D0
      use memory_allocator
C
      implicit none
C
      integer I,J
C
#include "dcbprp.h"
#include "dcbxlr.h"
!#include "codata.h"
C
      logical, allocatable :: DOATOM(:)

      double precision, allocatable :: PMAT(:,:,:,:),TMPMAT(:,:,:),
     &                                 CSTRA(:),SCTRA(:),GVAL(:)
C
      INTEGER IBUFA(3*NUCDEP),IBUFB(3),IFREQ,NCMAT,NATOM,JATOM,
     &        NSCOOR,IA,IB,JA,JB,KA,KB,NA,NB,IABLRF_EE,IABLRF_EP
      DOUBLE PRECISION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &                 ABCNV(NBFREQ,NLRBPT)
C
      JATOM  = NATOMS + NFLOAT
      NCMAT  = 9*NUCDEP
      NSCOOR = 3*NUCDEP
      allocate(DOATOM(NUCIND))
      allocate(PMAT(3,3,NUCDEP,2))
      allocate(CSTRA(NSCOOR*NSCOOR))
      allocate(SCTRA(NSCOOR*NSCOOR))
      allocate(TMPMAT(3,NSCOOR,2))
      allocate(GVAL(NUCDEP))
C
      CALL QENTER('PRP_PVC_SHIELD')
C
      CALL TITLER('PARITY-VIOLATION CONTRIBUTION TO '//
     &  'NMR SHIELDING TENSOR','*',116)
C
      WRITE(LUPRI,'(10X,A)')
     & 'According to A. L. Barra et al., '//
     & 'Europhys. Lett. 5 (1988) p.217'
      WRITE(LUPRI,'(19X,A)')
     & 'and R. Bast et al., JCP 125 (2006) p.064504'
      CALL PRSYMB(LUPRI,'-',81,0)
      WRITE(LUPRI,*) ''
C
C     INTTYP = 62 (Integrals for Parity Violation - chirality)
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C     **********************************************************************
C     ********** Extract PV contribution to NMR shielding tensors **********
C     **********************************************************************
C
C     The NUCDEP PV operators as A operators in <<A;B>> :
      NA = 0
      DO I = 1,NSCOOR
        IF(IP_PVC_SHIELD(I).GT.0) THEN
          NA = NA + 1
          IBUFA(NA) = I
        ENDIF
      END DO
      IF (NA.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_SHIELD: NA=',NA
        CALL QUIT('PRP_PVC_SHIELD: Wrong value of NA!')
      ENDIF
C
C     The 3 magnetic moment directions as B operators in <<A;B>> :
      NB = 0
      DO I = 1,3
        IF(IP_PVC_SHIELD(I+MXCOOR).GT.0) THEN
          NB = NB + 1
          IBUFB(NB) = I
        ENDIF
      END DO
      IF (NB.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_SHIELD: NB=',NB
        CALL QUIT('PRP_PVC_SHIELD: Wrong value of NB!')
      ENDIF
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C--------------------------------------------------------------------
C
C     LINEAR RESPONSE PART OF PVC-shielding TENSOR
C
C     Fill matrices of parity-violation contribution to NMR shieldings
C
      IF (XLR_SKIPEE) THEN
        IABLRF_EE = -999999999
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 2
        END IF
      ELSE
        IABLRF_EE = 2
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 3
        END IF
      END IF
C
C
      DO IFREQ = 1,NBFREQ
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
           CALL PRSYMB(LUPRI,'-',72,2)
           WRITE(LUPRI,'(A,F15.8,A)')
     &          '*** Frequency    :',BFREQ(IFREQ),' a.u.'
        ENDIF
        CALL DZERO(TMPMAT,2*NCMAT)
C
        DO IB = 1,NB
          JB = IBUFB(IB)
          KB = IP_PVC_SHIELD(MXCOOR+JB)
          JB = IPTAX(JB,2)
          DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IP_PVC_SHIELD(JA)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         TMPMAT(JB,JA,1) = ABLRF(IFREQ,KA,KB,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         TMPMAT(JB,JA,2) = ABLRF(IFREQ,KA,KB,IABLRF_EP)
          ENDDO
        ENDDO
C
C     Transform to C1 basis
C
        CALL DZERO(PMAT,2*NCMAT)
C
        IF (.NOT. XLR_SKIPEE) THEN
C           .... extract the e-e response part into PMAT(...,1)
            CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,NSCOOR,2,2)
        END IF
        IF (.NOT.XLR_SKIPEP) THEN
C           ... extract the e-p response part into PMAT(...,2)
            CALL TRADIP(TMPMAT(1,1,2),PMAT(1,1,1,2),
     &                CSTRA,SCTRA,NSCOOR,2,2)
        END IF
      END DO
C
C--------------------------------------------------------------------
      CALL PVCSHI(PMAT,JATOM,IPRPRP,DOATOM)
C--------------------------------------------------------------------
C
      CALL QEXIT('PRP_PVC_SHIELD')
C
      deallocate(DOATOM)
      deallocate(GVAL)
      deallocate(PMAT)
      deallocate(CSTRA)
      deallocate(SCTRA)
      deallocate(TMPMAT)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pvcshi */
      SUBROUTINE PVCSHI(PMAT,KATOM,IPRINT,DOATOM)
C***********************************************************************
C
C     Analyze routine for the parity-violation contribution
C     to NMR shielding constants
C
C     This routine is modelled after SPRRES subroutine
C
C     Agustin Aucar - November 2022
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CHARGE, ISOTOP, ISTBNU, IZATOM, NAMDEP,
     &  NUCDEP, NUCIND
      use symmet_mod, only: MAXOPR
!      use spinro_mod, only: GTRANT
!      use codata_mod, only: CVEL, GFERMI, S2THETAW, XFMP
      use consts_mod, only: D0, D1, D2, D4, D6
      implicit none
C
      integer I,J
C
#include "spinro.h"
!#include "dgroup.h"
!#include "dcbham.h"
#include "codata.h"
C
      integer K,KATOM,IPRINT,IATOM,ISYMOP,JATOM,NATTYP,NUMIS,IJ
C
      double precision FACTOR,FAC,DISOTP,PMAT(3,3,KATOM,2),
     &                 GTRANPEE(3,3,KATOM),GTRANPEP(3,3,KATOM),
     &                 GVAL(KATOM),SYM_TENSOR(3,3,KATOM),
     &                 PRINCIPAL_VAL(9),PRINCIPAL_AX(3,3),
     &                 SOME(3),ISOME(3),AMASS(KATOM)
C
      LOGICAL DOATOM(NUCIND)
C
      FACTOR = -GFERMI*(D1-D4*S2THETAW)*D2*XFMP*1.0D6/(CVEL*SQRT(D2))
C
C     Set up g-values of all nuclei
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         NATTYP = IZATOM(IATOM)
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
         IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
           JATOM = JATOM + 1
           AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
           GVAL(JATOM)  = DISOTP(NATTYP,NUMIS,'GVAL')
         END IF
 110     CONTINUE
 100  CONTINUE
C
C    Print PVC to NMR shielding tensors
C
      IATOM = 0
      DO 98 K = 1, NUCIND
       DO 200 ISYMOP = 0, MAXOPR
        IF (IAND(ISTBNU(K),ISYMOP).EQ.0) THEN
         IATOM = IATOM + 1
         IF (DOATOM(K)) THEN
          IF (CHARGE(IATOM) .GT. D0) THEN
C
C    Print g-values of selected nuclei
C
           IF (GVAL(IATOM) .NE. D0) THEN
            FAC=FACTOR/GVAL(IATOM)
            IF (IPRINT .GE. 0) THEN
             CALL HEADER('Parity-violation contribution to '//
     &             'NMR shielding tensor (in ppm) for '     //
     &                     NAMDEP(IATOM),-1)
             WRITE (LUPRI,'(2X,A17,F12.6,/)')
     &                 'Nuclear g-value: ',GVAL(IATOM)
            END IF
           ELSE   !   if GVAL = 0
            FAC=FACTOR
            IF (IPRINT .GE. 0) THEN
             CALL HEADER('Artificial parity-violation contribution '//
     &             'to NMR shielding tensor (in ppm) for '          //
     &                     NAMDEP(IATOM),-1)
             WRITE (LUPRI,'(2X,A,F8.3)')
     &  'WARNING: Nuclear g-value NOT available for isotope of mass ',
     &  AMASS(IATOM)
             WRITE (LUPRI,'(2X,A,/)')
     &  'Nuclear g-value is artificially set equal to 1'
            END IF
           END IF
C
           CALL DZERO(GTRANT,9*NUCDEP)
           CALL DZERO(GTRANPEE,9*NUCDEP)
           CALL DZERO(GTRANPEP,9*NUCDEP)
C
           DO 50 I = 1, 3
            DO 50 J = 1, 3
C           ee and pp parts of the PVC-shieldings
                GTRANPEE(I,J,IATOM) = FAC*PMAT(I,J,IATOM,1)
                GTRANPEP(I,J,IATOM) = FAC*PMAT(I,J,IATOM,2)
C           Total PV contribution to shieldings
                GTRANT(I,J,IATOM) = GTRANPEE(I,J,IATOM)
     &                         + GTRANPEP(I,J,IATOM)
 50        CONTINUE
C
           IF (IPRINT .GE. 0) THEN
C
C           Total PVC-shieldings
C
             CALL HEADER('Total PVC to NMR shielding '//
     &              'tensor',-1)
             WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'x', (GTRANT(1,I,IATOM),I = 1,3),
     &                 'y', (GTRANT(2,I,IATOM),I = 1,3),
     &                 'z', (GTRANT(3,I,IATOM),I = 1,3)
             WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANT(1,1,IATOM)
     &              +GTRANT(2,2,IATOM)+GTRANT(3,3,IATOM))/3
           END IF
C
C           Linear response contributions to PVC-shieldings
C
           IF (IPRINT .GE. 1) THEN
             CALL HEADER('sigma^PV(e-e)',-1)
             WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'x', (GTRANPEE(1,I,IATOM),I = 1,3),
     &                 'y', (GTRANPEE(2,I,IATOM),I = 1,3),
     &                 'z', (GTRANPEE(3,I,IATOM),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANPEE(1,1,IATOM)
     &              +GTRANPEE(2,2,IATOM)+GTRANPEE(3,3,IATOM))/3
             CALL HEADER('sigma^PV(e-p)',-1)
             WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'x', (GTRANPEP(1,I,IATOM),I = 1,3),
     &                 'y', (GTRANPEP(2,I,IATOM),I = 1,3),
     &                 'z', (GTRANPEP(3,I,IATOM),I = 1,3)
             WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANPEP(1,1,IATOM)
     &              +GTRANPEP(2,2,IATOM)+GTRANPEP(3,3,IATOM))/3
           END IF
C
C           PVC-shieldings in principal axes
C
           CALL DUNIT(PRINCIPAL_AX,3)
           CALL DZERO(PRINCIPAL_VAL,9)
           IJ = 1
           DO I = 1,3
            DO J = 1,I
              SYM_TENSOR(I,J,IATOM) =
     &                  (GTRANT(I,J,IATOM) + GTRANT(J,I,IATOM))/D2
C
              PRINCIPAL_VAL(IJ) = SYM_TENSOR(I,J,IATOM)
              IJ = IJ + 1
            ENDDO
           ENDDO
           CALL JACO(PRINCIPAL_VAL,PRINCIPAL_AX,3,3,3,SOME,ISOME)
C
           IF (IPRINT .GE. 0) THEN
             CALL HEADER('Principal components',-1)
             WRITE (LUPRI,'(10X,E16.8)')
     &                 PRINCIPAL_VAL(1),
     &                 PRINCIPAL_VAL(3),
     &                 PRINCIPAL_VAL(6)
             WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (PRINCIPAL_VAL(1)
     &              +PRINCIPAL_VAL(3)+PRINCIPAL_VAL(6))/3
             CALL HEADER('Principal axes',-1)
             WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'x', (PRINCIPAL_AX(1,I),I=1,3),
     &                 'y', (PRINCIPAL_AX(2,I),I=1,3),
     &                 'z', (PRINCIPAL_AX(3,I),I=1,3)
           END IF
          END IF    ! (CHARGE(IATOM).GT.D0)
         END IF     ! (DOATOM(K))
        END IF
 200   CONTINUE
 98   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PRINTOUT_PVC_SPINSPIN */
      SUBROUTINE PRINTOUT_PVC_SPINSPIN(ABLRF,ABCNV,WORK,KFREE,LFREE)
C=======================================================================
C     PV contribution to spin-spin coupling
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
C
      DIMENSION WORK(*)
C
C
      CALL MEMGET('LOGI',KDOATOM,NUCIND,WORK,KFREE,LFREE)
      CALL SETATM(WORK(KDOATOM),NATOM,-29)
      CALL PRINTOUT_PVC_SPINSPIN_SUB(ABLRF,ABCNV,
     &                               WORK(KDOATOM),
     &                               WORK,KFREE,LFREE)
      CALL MEMREL('PRINTOUT_PVC_SPINSPIN',WORK,1,KFREE,KFREE,LFREE)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PRINTOUT_PVC_SPINSPIN_SUB */
      SUBROUTINE PRINTOUT_PVC_SPINSPIN_SUB(ABLRF,ABCNV,DOATOM,
     &                                     WORK,KFREE,LFREE)
C=======================================================================
C     PV contribution to spin-spin coupling
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: IP_PVC_SPINSPIN
!      use dcbxlr_mod, only: NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEP
      use nuclei_mod, only: IZATOM, NAMDEP, NUCDEP
!      use codata_mod, only: CVEL, GFERMI, S2THETAW, XFMP, XTHZ
      use pi_mod, only: PI
#include "implicit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "dcbprp.h"
!#include "dcbexp.h"
#include "dcbxlr.h"
!#include "dcbxpr.h"
!#include "dcbprl.h"
#include "codata.h"
C
      LOGICAL   DOATOM(NUCDEP)
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),
     &          WORK(*)
      DOUBLE PRECISION FACTOR
C
!#include "dgroup.h"
C
      FACTOR = (GFERMI*(D1-D4*S2THETAW))/(CVEL*SQRT(D2))
      FACTOR = (FACTOR*XTHZ)/(D2*PI*CVEL*CVEL)
      FACTOR = FACTOR*1.0D9
C
      CALL TITLER('PARITY-VIOLATION CONTRIBUTION TO '//
     &  'SPIN-SPIN COUPLING','*',116)
C
      WRITE(LUPRI,'(10X,A)')
     & 'According to A. L. Barra et al., '//
     & 'Europhys. Lett. 5 (1988) p.217'
      CALL PRSYMB(LUPRI,'-',81,0)
      WRITE(LUPRI,*) ''
C
      CALL HEADER('Total isotropic PVC to spin-spin J-couplings '//
     &              '(in nHz)',-1)
C
C     NOTE that we here exploit that we are in C1 symmetry,
C     so that all nuclei are symmetry-independent ! (tsaue)
      DO IATOM = 1,NUCDEP-1
      IF(DOATOM(IATOM)) THEN
        DO JATOM = IATOM+1,NUCDEP
        IF(DOATOM(JATOM)) THEN
C
C
C         extract nuclear g factors
C         =========================
C
          IZ = IZATOM(IATOM)
          ABUNDANCE = D0
          DO ISO = 1,5
            GVAL_IATOM = DISOTP(IZ,ISO,'GVAL')
            IF(GVAL_IATOM .NE. D0) GOTO 1
          ENDDO
 1        CONTINUE
          IZ = IZATOM(JATOM)
          ABUNDANCE = D0
          DO ISO = 1,5
            GVAL_JATOM = DISOTP(IZ,ISO,'GVAL')
            IF(GVAL_JATOM .NE. D0) GOTO 2
          ENDDO
 2        CONTINUE
C
          RATIO_IATOM = GVAL_IATOM/(D2*XFMP)
          RATIO_JATOM = GVAL_JATOM/(D2*XFMP)
C
C
C         extract response
C         ================
C
          RESPONSE = D0
          DO IXYZ = 1,3
            KA = IP_PVC_SPINSPIN(1,(IATOM*3)-3+IXYZ)
            KB = IP_PVC_SPINSPIN(2,(JATOM*3)-3+IXYZ)
            RESPONSE = RESPONSE
     &               + RATIO_JATOM*ABLRF(1,KA,KB,2)
     &               + RATIO_IATOM*ABLRF(1,KB,KA,2)
            IF(.NOT. XLR_SKIPEP) THEN
              RESPONSE = RESPONSE
     &                 + RATIO_JATOM*ABLRF(1,KA,KB,3)
     &                 + RATIO_IATOM*ABLRF(1,KB,KA,3)
            ENDIF
          ENDDO
          RESPONSE = RESPONSE*FACTOR/D3
C
          WRITE(LUPRI,'(2A10,E20.10)')
     &          NAMDEP(IATOM), NAMDEP(JATOM), RESPONSE
        ENDIF
        ENDDO
      ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'-',46,0)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PRP_PVC_EFG */
      SUBROUTINE PRP_PVC_EFG(ABLRF,ABCNV)
C=======================================================================
C     PV contributions to EFG
C-----------------------------------------------------------------------
C     Juan J. Aucar                                    March 2025
C     doi Ref.: https://doi.org/10.1039/D4CP04840G
C=======================================================================
      use priunit_mod, only: LUPRI
      use consts_mod, only: D0, D1, D2, D4
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: ABUND, IPRPRP, IP_PVC_EFG
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEE, XLR_SKIPEP
      use nuclei_mod, only: GAUNUC, ISOTOP, IZATOM, NAMDEP, NAMN,
     &  NUCDEG, NUCIND
!      use codata_mod, only: GFERMI, S2THETAW
      use memory_allocator
C
      implicit none
C
      integer I,J
C
#include "dcbprp.h"
#include "dcbxlr.h"
#include "codata.h"
C
      logical, allocatable :: DOATOM(:)
      integer, allocatable :: ISOT(:)
      double precision, allocatable :: PVEFGMAT(:,:,:), FACTORS(:)
C
      LOGICAL TEST
      INTEGER IBUF(6*NUCIND,2),IXYZ,JXYZ,IATOM,JATOM,L,NATOM,
     &        IATINF(3,NUCIND),IDEG,IDEP,IFREQ,INUC,ISO,IZ,NISO,
     &        IA,IB,JA,JB,K,KA,KB,NA,NB,NSCOOR,NCHES,NSPEC,
     &        IABLRF_EE,IABLRF_EP,MASSNUM(NUCIND),ISOFLAG,
     &        NSCOORA, NSCOORB, ITENSOR, NPROTON,NNEUTRON,ICHARGE,
     &        IMASS
      DOUBLE PRECISION FACTOR,DISOTP,ABND,SELEC_ABND,GVALUE,
     &                 ABCNV(NBFREQ,NLRBPT),RATINF(3,5,NUCIND),
     &                 ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &                 RESPEE(3,3,NUCIND),RESPEP(3,3,NUCIND),
     &                 RESP(3,3,NUCIND),QW, QW_SELECTED(NUCIND),
     &                 IMASS_SELECTED(NUCIND)
C
      NSCOOR = 3*NUCIND
      NSCOORA = 6*NUCIND
      NSCOORB = NUCIND
      NCHES  = NSCOOR*NSCOOR
      allocate(DOATOM(NUCIND))
      allocate(PVEFGMAT(6*NUCIND,NUCIND,2))
      allocate(FACTORS(NUCIND))
C
      CALL QENTER('PRP_PVC_EFG')
C
      FACTOR = 1
      FACTORS = 1
      CALL DZERO(PVEFGMAT,6*NUCIND*NUCIND*2)
C
      CALL TITLER('PARITY-VIOLATION CONTRIBUTION TO '//
     &  'EFG','*',116)
!      WRITE(LUPRI,*) ''
      WRITE(LUPRI,'(6X,A,/,6X,A,/,6X,A)')
     & 'According to J. J. Aucar and A. F. Maldonado,',
     & 'Parity Violation Effects on the Electric Field Gradient,',
     & 'Phys. Chem. Chem. Phys., 27, 7594 (2025)'
      CALL PRSYMB(LUPRI,'-',56,6)
      WRITE(LUPRI,*) ''

C
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
C
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      ISOFLAG = 0
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',63,0)
      WRITE(LUPRI,'(A,9X,A,11X,A)')
     &  '@ name  charge  isotope   mass','abundance','Qw'
      CALL PRSYMB(LUPRI,'-',63,0)
      NSPEC = 0
      IDEP  = 0
      DO INUC = 1,NUCIND
        IDEG = NUCDEG(INUC)
        IF(DOATOM(INUC)) THEN
          TEST = .FALSE.
          IZ   = IZATOM(INUC)
          NISO = 0
!         Abundance of the selected isotope
          SELEC_ABND = DISOTP(IZ,ISOTOP(INUC),'ABUNDANCE')
          IF (SELEC_ABND.LT.ABUND) ISOFLAG = ISOFLAG + 1
          DO ISO = 1,5
            ABND   = DISOTP(IZ,ISO,'ABUNDANCE')
            IMASS = DISOTP(IZ,ISO,'A')
            NPROTON  = IZ
            NNEUTRON = IMASS - IZ
            QW = (D1 - D4*S2THETAW) * NPROTON -
     &      D1 * NNEUTRON

            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = QW
            ENDIF
          ENDDO

!         Store values for selected isotopes
          IMASS = DISOTP(IZ,ISOTOP(INUC),'A')
          NNEUTRON = IMASS - IZ
          QW_SELECTED(INUC) = (D1 - D4*S2THETAW) * NPROTON -
     &      D1 * NNEUTRON
          IMASS_SELECTED(INUC) = IMASS

!         Print section for isotopes above abundance threshold
          IF(NISO.GT.0) THEN
            IATINF(1,NSPEC+1) = IDEP + 1
            IATINF(2,NSPEC+1) = IZ
            IATINF(3,NSPEC+1) = NISO
            WRITE(LUPRI,
     &      '(A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F12.4)')
     &      '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
            DO J = 2,NISO
              WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F12.4)')
     &        '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
            ENDDO
            DO K = 2,IDEG
              IATINF(1,NSPEC+K) = IDEP + K
              IATINF(2,NSPEC+K) = IZ
              IATINF(3,NSPEC+K) = NISO
              DO L = 1,NISO
                RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
                RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
                RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
              ENDDO
            ENDDO
            NSPEC = NSPEC + IDEG
          ENDIF
        ENDIF
        IDEP = IDEP + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
C
      IF (ISOFLAG.GT.0) THEN
      write(lupri,'(/,A)')
     &  'Some selected isotopes are below the abundances threshold.'
      write(lupri,'(A)')
     &  'Dimish the abundances keyword under PROPERTIES to get proper'
     &  //' nuclear information.'
      END IF
C
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        CALL QEXIT('PRP_PVC_EFG')
        RETURN
      ENDIF


!     Print section for selected isotopes
      WRITE(LUPRI,'(/,A)') '@ Selected nuclear species:'
      CALL PRSYMB(LUPRI,'=',63,0)
      WRITE(LUPRI,'(A,9X,A,11X,A)')
     &  '@ name  charge  isotope   mass','abundance','Qw'
      CALL PRSYMB(LUPRI,'-',63,0)

      DO INUC = 1, NUCIND
        IF(DOATOM(INUC)) THEN
        IZ   = IZATOM(INUC)
        ABND   = DISOTP(IZ,ISOTOP(INUC),'ABUNDANCE')
        write(lupri,'(A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F12.4)')
     &    '@',NAMN(INUC),IZ,ISOTOP(INUC),IMASS_SELECTED(INUC),
     &    ABND,QW_SELECTED(INUC)

!       FACTOR: Qw * Gf [au] / (2 * raiz(2))
        FACTORS(INUC) = (QW_SELECTED(INUC)*GFERMI)/(2*SQRT(D2))
        END IF
      END DO
      CALL PRSYMB(LUPRI,'-',61,0)
      WRITE(LUPRI,'(3/)')

C
C     **********************************************************************
C     ******************* Extract PV contributions to EFG ******************
C     **********************************************************************
C
      NA = 0
      NB = 0

C     The NSCOORA PV operators as A operators in <<A;B>> :
      DO I = 1,NSCOORA
        IF(IP_PVC_EFG(1,I).GT.0) THEN
          NA = NA + 1
          IBUF(NA,1) = I
        ENDIF
      END DO
C     The NSCOORB EFG operators as B operators in <<A;B>> :
      DO I = 1,NSCOORB
        IF(IP_PVC_EFG(2,I).GT.0) THEN
          NB = NB + 1
          IBUF(NB,2) = I
        ENDIF
      END DO
C
      IF (NA.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_EFG: NA=',NA
        CALL QUIT('PRP_PVC_EFG: Wrong value of NA!')
      ENDIF
      IF (NB.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_EFG: NB=',NB
        CALL QUIT('PRP_PVC_EFG: Wrong value of NB!')
      ENDIF
C
C     LINEAR RESPONSE PART OF PVC EFG
C     Fill matrices of parity-violation contribution to EFG
C
      IF (XLR_SKIPEE) THEN
        IABLRF_EE = -999999999
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 2
        END IF
      ELSE
        IABLRF_EE = 2
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 3
        END IF
      END IF
C
C
      DO IFREQ = 1,NBFREQ
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
          CALL PRSYMB(LUPRI,'-',72,2)
          WRITE(LUPRI,'(A,F15.8,A)')
     &         '*** Frequency    :',BFREQ(IFREQ),' a.u.'
        ENDIF
C
        DO IA = 1,NA
          JA = IBUF(IA,1)
          KA = IP_PVC_EFG(1,JA)
          DO IB = 1,NB
            JB = IBUF(IB,2)
            KB = IP_PVC_EFG(2,JB)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &      PVEFGMAT(JA,JB,1) = ABLRF(IFREQ,KA,KB,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &      PVEFGMAT(JA,JB,2) = ABLRF(IFREQ,KA,KB,IABLRF_EP)
          ENDDO
        ENDDO
      ENDDO
C
C     Analyze results
C
C     Here we exploit that we are in C1 symmetry,
C     so that all nuclei are symmetry-independent!
C
C     Initialize response matrices
      CALL DZERO(RESPEE,3*3*NUCIND)
      CALL DZERO(RESPEP,3*3*NUCIND)
      CALL DZERO(RESP,3*3*NUCIND)
C
C     Print some notes
      IF(.NOT.GAUNUC) THEN ! Point nucleus
        write(lupri,'(A,/,A)')
     &  'The PVC to EFG tensor strongly depends on the nuclear model.',
     &  'Gaussian model is suggested over point-like model.'
      ENDIF


C     Extract response
C     ================
C     First index in PVEFGMAT is related to the EFG operators
C     second index in PVEFGMAT is related to the PV spin-independent operators
      ITENSOR = 1
      DO IATOM = 1,NUCIND
        IF(DOATOM(IATOM)) THEN
          DO IXYZ = 1,3
          DO JXYZ = IXYZ,3
C           Get the ee and pp parts of the PVC to EFG
            RESPEE(JXYZ,IXYZ,IATOM) = FACTORS(IATOM)*
     &      PVEFGMAT(ITENSOR,IATOM,1)
            RESPEE(IXYZ,JXYZ,IATOM) = RESPEE(JXYZ,IXYZ,IATOM)
            RESPEP(JXYZ,IXYZ,IATOM) = FACTORS(IATOM)*
     &      PVEFGMAT(ITENSOR,IATOM,2)
            RESPEP(IXYZ,JXYZ,IATOM) = RESPEP(JXYZ,IXYZ,IATOM)
            ITENSOR = ITENSOR+1
          ENDDO
          ENDDO
!         Total PV contribution to EFG
          RESP = RESPEE + RESPEP
        ENDIF
      ENDDO

      WRITE(lupri,'(/)')
      IF (IPRPRP .LT. 1) CALL HEADER
     & ('Total PV contributions to EFG tensor (a.u.):',1)

C     Print the response tensor.
      DO IATOM = 1,NUCIND
        IF(DOATOM(IATOM)) THEN
          IF (IPRPRP .GE. 1) CALL PRSYMB(LUPRI,'=',91,0)
          WRITE (LUPRI ,'(2X,A9,A4)')'Nucleus: ',NAMDEP(IATOM)
          IF (IPRPRP .GE. 1) CALL HEADER
     &    ('Total PV contributions to EFG tensor (a.u.):',1)
          WRITE (LUPRI,'(19X,3(A,13X))') 'x', 'y', 'z'
          WRITE(LUPRI,'(3(3X,A,6X,3E14.6,/))')
     &    'x', (RESP(1,I,IATOM), I = 1,3),
     &    'y', (RESP(2,I,IATOM), I = 1,3),
     &    'z', (RESP(3,I,IATOM), I = 1,3)
C
C       Print out of e-e and e-p response contributions
          IF (IPRPRP .GE. 1) THEN
            CALL HEADER('e-p and e-e '//
     &            'response contributions (a.u.):',1)
            WRITE (LUPRI,'(12X,3(A,13X),4X,3(A,13X))')
     &             'x', 'y', 'z', 'x','y','z'
            WRITE (LUPRI,'(3(2X,A,3E14.6,4X,3E14.6,/))'),
     &  'x', (RESPEP(1,I,IATOM), I = 1,3),(RESPEE(1,I,IATOM), I = 1,3),
     &  'y', (RESPEP(2,I,IATOM), I = 1,3),(RESPEE(2,I,IATOM), I = 1,3),
     &  'z', (RESPEP(3,I,IATOM), I = 1,3),(RESPEE(3,I,IATOM), I = 1,3)
          END IF
C
        END IF
      END DO
      IF (IPRPRP .GE. 1) CALL PRSYMB(LUPRI,'=',91,0)
C
      CALL QEXIT('PRP_PVC_EFG')
C
      deallocate(FACTORS)
      deallocate(PVEFGMAT)
      deallocate(DOATOM)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRP_EFN */
      SUBROUTINE PRP_EFN(EXPVAL,EXP_CORR,CORR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by T. Saue Feb 26 2003
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: IPRPRP
      use nuclei_mod, only: NUCDEP, NUCIND
#include "implicit.h"
C
#include "dcbprp.h"
!#include "dcbexp.h"
!#include "dcbxlr.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL CORR
      logical, allocatable :: DOATOM(:)
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      allocate(DOATOM(NUCIND))
      CALL MEMGET('REAL',KEFN ,9*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,3*NUCDEP,WORK,KFREE,LFREE)
C
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
      CALL PRP_EFN1(EXPVAL,EXP_CORR,WORK(KEFN),WORK(KTEMP),DOATOM,
     &              CORR,IPRPRP)
      deallocate(DOATOM)
      CALL MEMREL('PRP_EFN',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRP_EFN1(EXPVAL,EXP_CORR,EFNVEC,TEMP,DOATOM,CORR,
     &                    IPRINT)
C*****************************************************************************
C
C     Calculate electric field at individual nuclei
C
C     Written by T.Saue - Feb 26 2003
C
C*****************************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbexp_mod, only: IPREXP
      use nuclei_mod, only: NAMDEP, NUCDEG, NUCDEP, NUCIND
#include "implicit.h"
      INTEGER :: I
C
#include "dcbexp.h"
      LOGICAL   DOATOM(*), CORR
      DIMENSION EXPVAL(*),EXP_CORR(*),EFNVEC(3,NUCDEP,3),TEMP(3,NUCIND)
C
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(A/)') '* Electric field at nuclei:'
C
      CALL EL_EFN(IPREXP,EXPVAL,DOATOM,EFNVEC(1,1,1),TEMP)
      CALL NUC_EFN(IPREXP,DOATOM,EFNVEC(1,1,2),TEMP)
C
C     print section
      IOFF = ICHAR('w')
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(12X,A/)')
     &     'Individual (non-zero) components'
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      CALL PRSYMB(LUPRI,'-',72,2)
      NATOM  = 0
      IATOMD = 0
      DO IATOM = 1, NUCIND
      IF(DOATOM(IATOM)) THEN
        DO IDEG = 1,NUCDEG(IATOM)
          IATOMD = IATOMD + 1
          NATOM  = NATOM + 1
          WRITE (LUPRI ,'(/2X,A9,A4,/)')'Nucleus: ',NAMDEP(IATOMD)
          DO I = 1,3
            EFNTOT=EFNVEC(I,NATOM,1)+EFNVEC(I,NATOM,2)
            WRITE(LUPRI,'(2X,A1,3(3X,F16.10,A))') CHAR(IOFF+I),
     &           EFNVEC(I,NATOM,1),' au ',EFNVEC(I,NATOM,2),' au ',
     &           EFNTOT,' au '
          END DO
        ENDDO
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
      END DO
      CALL PRSYMB(LUPRI,'-',72,2)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck el_efn */
      SUBROUTINE EL_EFN(IPRINT,EXPVAL,DOATOM,EFNEL,TEMP)
C*****************************************************************************
C
C     Calculates electronic contribution to electric fields
C     Based on EL_EFG
C
C     Written by T.Saue Feb 26 3004
C
C*****************************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbexp_mod, only: IPEFN
      use nuclei_mod, only: ISTBNU, NAMDEP, NUCDEG, NUCDEP, NUCIND
      use symmet_mod, only: ISYMAX, MAXREP, PT
#include "implicit.h"
      INTEGER :: I
C
#include "dcbexp.h"
      DIMENSION EXPVAL(*), TEMP(3,NUCIND),EFNEL(3,NUCDEP)
      LOGICAL DOATOM(*)
!#include "nqcc.h"
      CALL DZERO(EFNEL,3*NUCDEP)
      CALL DZERO(TEMP,3*NUCIND)
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
      IF (.NOT. DOATOM(IATOM)) GOTO 100
        DO 50 IREPC = 0, MAXREP
          DO 200 ICOOR = 1, 3
          IF (IAND(ISTBNU(IATOM),
     &        IEOR(IREPC,ISYMAX(ICOOR,1))).EQ.0) THEN
            ITYP = ITYP + 1
            TEMP(ICOOR,IATOM) = TEMP(ICOOR,IATOM)
     &                        + EXPVAL(IPEFN(ITYP))/NUCDEG(IATOM)
          END IF
 200      CONTINUE
  50    CONTINUE
 100  CONTINUE
      NATOM  = 0
      IATOMD = 0
      DO 110 IATOM = 1, NUCIND
      IF (DOATOM(IATOM)) THEN
        DO 310 IREPC = 0, MAXREP
        IF (IAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          DO 220 I = 1, 3
            EFNEL(I,NATOM) = TEMP(I,IATOM)*
     &                       PT(IAND(ISYMAX(I,1),IREPC))
 220      CONTINUE
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,A,A/)')
     &       'Electronic electric field at nucleus ',NAMDEP(IATOMD)
            IOFF = ICHAR('w')
            DO I = 1,3
              WRITE(LUPRI,'(4X,A1,A,E20.10)') CHAR(IOFF+I),
     &            '-component: ',EFNEL(I,NATOM)
            ENDDO
          END IF
        ENDIF
 310    CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(IATOM)
      ENDIF
 110  CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_RHONUC */
      SUBROUTINE DEF_RHONUC(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operator for electronic density at nuclei
C
C     Written by T.Saue May 8 2004
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     INTTYP = 9: One-electron Fermi contact integrals
      CALL SETATM(DOATOM,NATOM,-9)
      CALL DEF_RHONU1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_RHONU1 */
      SUBROUTINE DEF_RHONU1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define electronic density at nuclei
C     Written by T.Saue May 8 2004
C
C***********************************************************************
      use pi_mod, only: PI
!      use codata_mod, only: GFAC
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: ISTBNU, NAMDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTNUC, MAXREP
      use chrnos_mod, only: CHRNOS
!      use dcbexp_mod, only: LEXPP, MAXEXP, NEXPP
      use kinds_mod, only: dp
#include "implicit.h"
#include "codata.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0, D4 = 4.0D0, D3 = 3.0D0
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND)
C
!#include "dcbxpr.h"
#include "dcbexp.h"
C
C
      DFAC = D4*GFAC*PI/D3
C
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1/DFAC
      DO 100 IREP = 0, MAXREP
        DO 200 IATOM = 1, NUCIND
        IF (DOATOM(IATOM)) THEN
          IF (IAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
            PLABEL(1) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                        //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
            PNAME = 'Rho at nuc '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                            //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                            //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
            CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                   INDXPR,ISYXPR,ITRXPR,IPRINT)
            CALL OP1IND('RHONUC',IDUM,LEXPP,NEXPP,INDXPR,MAXEXP)
          END IF
        END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_rotg */
      SUBROUTINE DEF_ROTG(PNAME,PLABEL,PFAC,IPRINT)
C*****************************************************************************
C
C     Define operators for molecular rotational g-tensor
C     DEF_SHIELD as skeleton
C     and based on DEF_SPINRO as well (Agustin Aucar)
C
C     Written by T. Enevoldsen - Jul. 1998
C     Last revision: MI, jan.2003 - added from TEC's code
C                    MI, Sept.2006 - making functional
C                    Agustin Aucar, November 2019
C
C*****************************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbgen_mod, only: CVAL
!      use dcbprp_mod, only: IPROTG
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
!      use rotg_mod, only: LINDET
!      use dgroup_mod, only: LINEAR
      use memory_allocator
C
      implicit none
C
      integer I
C
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
#include "rotg.h"
#include "dgroup.h"
C
      INTEGER IPRINT,INDOP1,INDXPR,IOP,IPTYP,ISYXPR,ITRXPR,NPCOMP
      DOUBLE PRECISION PFAC,DL,DS
C
      CHARACTER PNAME*16, PLABEL(2)*8
      DIMENSION PFAC(2),IOP(3)
C     ANGMOM integrals correspond to -iL, so we set DL = -1
C     XSIGMA,YSIGMA and ZSIGMA come with imaginary i, so we set DS positive (DS = 0.5)
      PARAMETER(DL = -1.0D0,DS = 0.5D0)
C
      CALL QENTER('DEF_ROTG')
C
      LINDET=.FALSE.    ! to detect if linear symmetry is used
C
      CALL IZERO(IPROTG,9)
C
C     Here we set up the operators needed for the linear reponse part
C
C     Define left-hand operators: i.L and i.S
C     =======================================
C
C     Orbital angular momentum
C
       PNAME     = 'X Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'XANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XANGMOM ',IPROTG(1),LLRAPU,NLRAPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Y Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'YANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('YANGMOM ',IPROTG(2),LLRAPU,NLRAPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Z Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'ZANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('ZANGMOM ',IPROTG(3),LLRAPU,NLRAPT,
     &                       INDXPR,MAXLLR)
C
C     Spin angular momentum
C
       PNAME     = 'X Spin ang.mom. '
       IPTYP     = 10   ! Sigma_x
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XSPINMOM',IPROTG(4),LLRAPU,NLRAPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Y Spin ang.mom. '
       IPTYP     = 11   ! Sigma_y
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('YSPINMOM',IPROTG(5),LLRAPU,NLRAPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Z Spin ang.mom. '
       IPTYP     = 12   ! Sigma_z
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('ZSPINMOM',IPROTG(6),LLRAPU,NLRAPT,
     &                       INDXPR,MAXLLR)
C
C
C     Define right-hand operators: -(c/2) i.alpha x r
C     ===============================================
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
C
      PFAC(1)  = -0.5D0*CVAL
      PFAC(2)  = -0.5D0*CVAL
      NPCOMP   = 2
C
      IF (LINEAR) THEN  ! linear molecule detected using symmetry
C
Ciaa   LINDET is used to record that LINEAR symmetry is used,
C      even when in further steps it could be changed, as occur
C      for UKB, spin-free and ZORA calculations.
       LINDET=.TRUE.    ! save that linear symmetry is used
C
C        X component of B-field
C
         PNAME     = 'X magnetic field'
         IPTYP     = IOP(1)
         PLABEL(1) = 'ZDIPLEN'
         PLABEL(2) = 'YDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('XFIELD ',IPROTG(7),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
C
      ELSE ! if LINEAR is not TRUE
C
C        X component of B-field
C
         PNAME     = 'X magnetic field'
         IPTYP     = IOP(1)
         PLABEL(1) = 'ZDIPLEN'
         PLABEL(2) = 'YDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('XFIELD ',IPROTG(7),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
C
C        Y component of B-field
C
         PNAME     = 'Y magnetic field'
         IPTYP     = IOP(2)
         PLABEL(1) = 'XDIPLEN'
         PLABEL(2) = 'ZDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('YFIELD ',IPROTG(8),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
C
C        Z component of B-field
C
         PNAME     = 'Z magnetic field'
         IPTYP     = IOP(3)
         PLABEL(1) = 'YDIPLEN'
         PLABEL(2) = 'XDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
         CALL OP1IND('ZFIELD ',IPROTG(9),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
C
      END IF !LINEAR
C
C
C     Print section for molecular rotational g-tensor
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Molecular rotational g-tensor:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     Print reference
C
      WRITE(LUPRI,'(4X,A)')
     & 'Reference: I. A. Aucar et al., JCP 141, 194103 (2014)'
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     & '** A operators for molecular rotational g-tensor **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,6
        INDOP1 = IPROTG(I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRAPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     & '** B operators for molecular rotational g-tensor **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3
        INDOP1 = IPROTG(I+6)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
      CALL QEXIT('DEF_ROTG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ROTG */
      SUBROUTINE PRP_ROTG(ABLRF,ABCNV)
C***********************************************************************
C
C    Output routine for molecular rotational g-tensor
C
C   Written by T.Enevoldsen Jul. 1998
C
C   Revisions: MI, jan.2003 - added from TEC's 'left' code
C              MI, Prievidza/Odense, August 2006
C              Agustin Aucar, November 2019
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use orgcom_mod, only: CMXYZ, GAGORG
!      use dcbprp_mod, only: IPROTG, IPRPRP
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEE, XLR_SKIPEP
      use nuclei_mod, only: NATOMS, NFLOAT, NUCDEP
!      use rotg_mod, only: GFANUC, LINDET
!      use dgroup_mod, only: LINEAR
!      use dcbham_mod, only: SPINFR, URKBAL, ZORA
      use consts_mod, only: D0
      use memory_allocator
C
      implicit none
C
      integer I,J
C
!#include "dcborb.h"
#include "dcbprp.h"
!#include "dcbxpr.h"
!#include "dcbprl.h"
#include "dcbxlr.h"
#include "rotg.h"
#include "dgroup.h"
#include "dcbham.h"
C
      integer, allocatable :: NAT(:),ISOT(:)
      double precision, allocatable :: PMATL(:,:,:),PMATS(:,:,:),
     &                                 GEOM(:,:),MASS(:),AMASS(:)
C
      integer IBUFA1(3),IBUFA2(3),IBUFB(3),IA1,IA2,IB,JA1,JA2,JB,
     &        KA1,KA2,KB,NA1,NA2,NB,JATOM,IABLRF_EE,IABLRF_EP,IFREQ
      double precision ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &                 ABCNV(NBFREQ,NLRBPT)
C
      JATOM = NATOMS + NFLOAT
      allocate(PMATL(3,3,2)) ! L part of g-tensor for ee and ep contributions.
      allocate(PMATS(3,3,2)) ! S part of g-tensor for ee and ep contributions.
      allocate(GEOM(NUCDEP,3))
      allocate(MASS(NATOMS))
      allocate(NAT(NATOMS))
      allocate(ISOT(NATOMS))
      allocate(AMASS(NUCDEP))
C
      CALL QENTER('PRP_ROTG')
C
      CALL TITLER('MOLECULAR ROTATIONAL G-TENSOR','*',116)
      WRITE(LUPRI,'(13X,A)')
     & 'According to I. A. Aucar et al., JCP 141 (2014) p.194103'
      CALL PRSYMB(LUPRI,'-',81,0)
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
C     *** Extract molecular rotational g-tensor ***
C
C     The 3 L orbital angular momentum directions as A operators in <<A;B>> :
      NA1 = 0
      DO I = 1,3
        IF(IPROTG(I).GT.0) THEN
          NA1 = NA1 + 1
          IBUFA1(NA1) = I
        ENDIF
      END DO
      IF (NA1.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_ROTG: NA1=',NA1
        CALL QUIT('PRP_ROTG: Wrong value of NA1!')
      ENDIF
C
C     The 3 S spin angular momentum directions as A operators in <<A;B>> :
      NA2 = 0
      DO I = 1,3
        IF(IPROTG(I+3).GT.0) THEN
          NA2 = NA2 + 1
          IBUFA2(NA2) = I
        ENDIF
      END DO
      IF (NA2.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_ROTG: NA2=',NA2
        CALL QUIT('PRP_ROTG: Wrong value of NA2!')
      ENDIF
C
Ciaa   As in the current implementation UKB, spin-free and ZORA
C      calculations of energy cannot be performed employing symmetry
C      (they set LINEAR to false in such cases, see subroutine GLINSM)
C      LINDET allows to recover LINEAR symmetry defined in DEF_ROTG
C      if it is lost due to the use of UKB, SF or ZORA functionalities
      IF ((URKBAL.or.SPINFR.or.ZORA).and.LINDET) LINEAR=.TRUE.
C
C     The 3 (c/2) alpha x r operators as B operators in <<A;B>> :
      NB = 0
      DO I = 1,3
        IF(IPROTG(I+6).GT.0) THEN
          NB = NB + 1
          IBUFB(NB) = I
        ENDIF
      END DO
      IF (NB.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_ROTG: NB=',NB
        CALL QUIT('PRP_ROTG: Wrong value of NB!')
      ENDIF
C
C--------------------------------------------------------------------
C
C     NUCLEAR PART OF MOLECULAR ROTATIONAL G-TENSOR
C
      CALL DZERO(GFANUC,9)
C
      CALL CMMASS(GEOM,MASS,NAT,ISOT,IPRPRP)
      CALL NUCMGF(GEOM,CMXYZ,GAGORG,IPRPRP)
C
C--------------------------------------------------------------------
C
C     LINEAR RESPONSE PART OF MOLECULAR ROTATIONAL G-TENSOR
C
C     Fill molecular rotational g-tensor matrices
C
      IF (XLR_SKIPEE) THEN
        IABLRF_EE = -999999999
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 2
        END IF
      ELSE
        IABLRF_EE = 2
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 3
        END IF
      END IF
C
C
      DO IFREQ = 1,NBFREQ
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
           CALL PRSYMB(LUPRI,'-',72,2)
           WRITE(LUPRI,'(A,F15.8,A)')
     &          '*** Frequency    :',BFREQ(IFREQ),'a.u.'
        ENDIF
        CALL DZERO(PMATL,2*9)
        CALL DZERO(PMATS,2*9)
C
        DO IB = 1,NB
          JB = IBUFB(IB)
          KB = IPROTG(6+JB)
          DO IA1 = 1,NA1
            JA1 = IBUFA1(IA1)
            KA1 = IPROTG(JA1)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         PMATL(JB,JA1,1) = ABLRF(IFREQ,KA1,KB,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         PMATL(JB,JA1,2) = ABLRF(IFREQ,KA1,KB,IABLRF_EP)
          ENDDO
C
          DO IA2 = 1,NA2
            JA2 = IBUFA2(IA2)
            KA2 = IPROTG(3+JA2)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         PMATS(JB,JA2,1) = ABLRF(IFREQ,KA2,KB,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         PMATS(JB,JA2,2) = ABLRF(IFREQ,KA2,KB,IABLRF_EP)
          ENDDO
        ENDDO
      END DO
C
C
      WRITE (LUPRI,'(/,A,18X,3F20.15)')
     &    ' Gauge origin (a.u.):', (GAGORG(I),I=1,3)
C
C--------------------------------------------------------------------
      CALL MGFRES(GEOM,AMASS,JATOM,PMATL,PMATS,IPRPRP)
C--------------------------------------------------------------------
      CALL QEXIT('PRP_ROTG')
C
      deallocate(AMASS)
      deallocate(ISOT)
      deallocate(NAT)
      deallocate(MASS)
      deallocate(GEOM)
      deallocate(PMATS)
      deallocate(PMATL)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucmgf */
      SUBROUTINE NUCMGF(GEOM,CMXYZ,ORIGIN,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Calculate nuclear contribution to molecular g-factors
C     K.Ruud, June-94
C
C***********************************************************************
C
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CHARGE, ISTBNU, NATOMS, NFLOAT, NUCIND
      use symmet_mod, only: MAXOPR
!      use rotg_mod, only: GFANUC
      use priunit_mod, only: LUPRI
      use consts_mod, only: D0
      implicit none
C
      integer I,J
C
#include "rotg.h"
C
      double precision GEOM,CMXYZ,ORIGIN
C
      DIMENSION GEOM(3*(NATOMS+NFLOAT)),CMXYZ(3),ORIGIN(3)
C
      integer IX,IY,IZ,IPRINT,ISYMOP,KK
      double precision XCOOR1,YCOOR1,ZCOOR1,XCOOR2,YCOOR2,ZCOOR2
C
      IF (IPRINT.GE.5) THEN
       CALL HEADER(
     & 'NUCMGF: Entering parameters to calculate nuclear'//
     & ' contribution to molecular g-factor',-1)
       write(lupri,*) 'Gauge origin=',ORIGIN(1),ORIGIN(2),ORIGIN(3)
       do i=1, (3*(NATOMS+NFLOAT))
        write(lupri,*) 'Entering geom(',i,')=',GEOM(i)
       enddo
      ENDIF
C
      CALL DZERO(GFANUC,9)
C
      KK = 1
      DO 10 I = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
              IF (CHARGE(I) .GT. D0) THEN
               XCOOR1 = GEOM(KK) - ORIGIN(1)
               YCOOR1 = GEOM(KK + 1) - ORIGIN(2)
               ZCOOR1 = GEOM(KK + 2) - ORIGIN(3)
               XCOOR2 = GEOM(KK) - CMXYZ(1)
               YCOOR2 = GEOM(KK + 1) - CMXYZ(2)
               ZCOOR2 = GEOM(KK + 2) - CMXYZ(3)
               GFANUC(1,1) = GFANUC(1,1) + CHARGE(I)*(YCOOR1*YCOOR2
     &                                       + ZCOOR1*ZCOOR2)
               GFANUC(2,2) = GFANUC(2,2) + CHARGE(I)*(XCOOR1*XCOOR2
     &                                       + ZCOOR1*ZCOOR2)
               GFANUC(3,3) = GFANUC(3,3) + CHARGE(I)*(XCOOR1*XCOOR2
     &                                       + YCOOR1*YCOOR2)
               GFANUC(1,2) = GFANUC(1,2) - CHARGE(I)*XCOOR1*YCOOR2
               GFANUC(1,3) = GFANUC(1,3) - CHARGE(I)*XCOOR1*ZCOOR2
               GFANUC(2,3) = GFANUC(2,3) - CHARGE(I)*YCOOR1*ZCOOR2
               KK = KK + 3
              END IF
            END IF
 20      CONTINUE
 10   CONTINUE
      GFANUC(2,1) = GFANUC(1,2)
      GFANUC(3,1) = GFANUC(1,3)
      GFANUC(3,2) = GFANUC(2,3)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mgfres */
      SUBROUTINE MGFRES(GEOM,AMASS,KATOM,PMATL,PMATS,IPRINT)
C***********************************************************************
C
C     Output routine for molecular rotational g-tensor
C     This routine is modelled after MGFRES (Prof. K Ruud) in DALTON,
C     and based on SPRRES routine (DIRAC).
C
C   INPUT:
C     PMATL - contribution from orbital angular moment to linear response term
C     PMATS - contribution from spin angular moment to linear response term
C
C     Agustin Aucar - November 2019
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CORD, ISOTOP, ISTBNU, IZATOM, NUCIND
      use symmet_mod, only: IPTAX, ISYMAX, MAXOPR, PT
!      use rotg_mod, only: GFANUC
      use orgcom_mod, only: CMXYZ
!      use dgroup_mod, only: LINEAR
!      use dcbprp_mod, only: LONDON
!      use codata_mod, only: PMASS, XFAMU
      use consts_mod, only: D1, D2
      implicit none
C
      integer I,J
C
#include "rotg.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "codata.h"
C
      integer IATOM,IX,IY,IPRINT,ISYMOP,JATOM,JX,JY,
     &        KATOM,NUMIS,ICOOR,NATTYP
C
      double precision GEOM(KATOM,3),AMASS(KATOM),ANGUMOM(3),
     &                 TINERT(9),OMEGA(3),EIGVAL(3),EIGVEC(3,3),
     &                 PMATL(3,3,2),PMATS(3,3,2),GFATOT,
     &                 GFAN,GFALR,GFALRLEE,GFALRLEP,GFALRSEE,
     &                 GFALRSEP,GTRANT(3,3),GTRANN(3,3),
     &                 GTRANP(3,3),GTRANPLEE(3,3),GTRANPLEP(3,3),
     &                 GTRANPSEE(3,3),GTRANPSEP(3,3),CONST,DISOTP
C
      LOGICAL PLANAR,LINMOL
C
Ciaa  Nuclear rotational angular momentum could be defined here to obtain Omega afterwards
      DO 300 I = 1,3
        ANGUMOM(I) = D1
 300  CONTINUE
C
C
C     Determine principal moments of inertia
C     ======================================
C
C     Set up coordinates (relative to center of mass) and masses of all nuclei
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         NATTYP = IZATOM(IATOM)
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
         IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
           JATOM = JATOM + 1
           AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
           DO 120 ICOOR = 1, 3
             GEOM(JATOM,ICOOR) = PT(IAND(ISYMAX(ICOOR,1),ISYMOP))
     &                          *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
 120       CONTINUE
         END IF
 110     CONTINUE
 100  CONTINUE
C
      CALL WLKDIN(GEOM,AMASS,KATOM,ANGUMOM,TINERT,OMEGA,EIGVAL,EIGVEC,
     &            .TRUE.,PLANAR,LINMOL,IPRINT)
C     On output
C        TINERT - inertia tensor
C        EIGVAL - inverse eigenvalues (for linear molecules, in descending order)
C        EIGVEC - corresponding eigenvectors
C        PLANAR=T : planar molecule
C        LINMOL=T : linear molecule (EIGVAL(1) = 0.0)
C
C
Ciaa  To obtain Omega, the nuclear rot. angular momentum (ANGUMOM) should be defined before
      IF (IPRINT .GE. 8) THEN
       CALL HEADER
     &    ('Nuclear angular velocity (a.u.)',0)
       WRITE (LUPRI,'(3X,A,2X,F16.8)')
     &       'Omega_A',OMEGA(1),
     &       'Omega_B',OMEGA(2),
     &       'Omega_C',OMEGA(3)
      END IF
C
C    Print inertia tensor or moment
C
      IF (IPRINT .GE. 0) THEN
       IF (LINEAR) THEN      ! Linear molecule detected using symmetry
              WRITE (LUPRI,'(/,2X,A,F16.8,/)')
     &         'Moment of inertia (a.u.): ', XFAMU/EIGVAL(3)
       ELSE
        IF(LINMOL) THEN      ! Linear molecule (without symmetry detection)
         WRITE (LUPRI,'(/,2X,A,F16.8,/)')
     &         'Moment of inertia (a.u.): ', XFAMU/EIGVAL(3)
        ELSE                 ! Non-linear molecule
         CALL HEADER
     &    ('Principal moments of inertia (a.u.) and principal axes',0)
         WRITE (LUPRI,'(3X,A,2X,F16.8,6X,3F16.8)')
     &       'IA',XFAMU/EIGVAL(1), (EIGVEC(1,I),I = 1,3),
     &       'IB',XFAMU/EIGVAL(2), (EIGVEC(2,I),I = 1,3),
     &       'IC',XFAMU/EIGVAL(3), (EIGVEC(3,I),I = 1,3)
         WRITE (LUPRI,'(/,2X,A)') 'Molecular rotational g-tensor '    //
     &                'elements are given in the principal axis system.'
        END IF
       END IF
      END IF
C
C    Print molecular rotational g-tensors
C
      IF (IPRINT .GE. 0) THEN
         IF (.NOT.LONDON) THEN
            WRITE (LUPRI,'(/,1X,A)')' Calculation without '//
     &                               'London orbitals.'
         ELSE
            WRITE (LUPRI,'(/,1X,A)')' London orbitals used.'
         END IF
      END IF
C
      IF (IPRINT .GE. 0) THEN
        CALL HEADER('Molecular rotational g-tensor',-1)
      END IF
C
C     Transform molecular rotational g-tensor to principal axis system,
C     and multiply by inverse moment of inertia and proton mass
C
C     Linear molecule detected using symmetry
C
      IF (LINEAR) THEN
        GFAN     = GFANUC(1,1)*PMASS*EIGVAL(3)
        GFALRLEE = PMATL(1,1,1)*PMASS*EIGVAL(3)*D2
        GFALRLEP = PMATL(1,1,2)*PMASS*EIGVAL(3)*D2
        GFALRSEE = PMATS(1,1,1)*PMASS*EIGVAL(3)*D2
        GFALRSEP = PMATS(1,1,2)*PMASS*EIGVAL(3)*D2
C       Total linear response contribution to GFA
        GFALR=GFALRLEE+GFALRLEP+GFALRSEE+GFALRSEP
C       Total molecular rotational g-factor
        GFATOT=GFAN+GFALR
C
        IF (IPRINT .GE. 0)
     &  WRITE(LUPRI,'((1X,A,F16.8,/),/,2(1X,A,F16.8,/),(1X,A,F16.8))')
     & ' Total molecular rotational g-factor          :',GFATOT,
     & '   Nuclear contribution to g-factor (g^nuc)   :',GFAN,
     & '   Electronic contribution to g-factor (g^LR) :',GFALR
C
        IF (IPRINT .GE. 1)
     &  WRITE(LUPRI,'(3(1X,A,F16.8,/),(1X,A,F16.8))')
     & '           g^LR-L(e-e) :',GFALRLEE,
     & '           g^LR-S(e-e) :',GFALRSEE,
     & '           g^LR-L(e-p) :',GFALRLEP,
     & '           g^LR-S(e-p) :',GFALRSEP
C
C     Linear symmetry is not detected
C
      ELSE
C
        CALL DZERO(GTRANT,9)
        CALL DZERO(GTRANN,9)
        CALL DZERO(GTRANP,9)
        CALL DZERO(GTRANPLEE,9)
        CALL DZERO(GTRANPLEP,9)
        CALL DZERO(GTRANPSEE,9)
        CALL DZERO(GTRANPSEP,9)
        DO 50 I = 1, 3
         DO 50 J = 1, 3
          DO 60 IX = 1, 3
           DO 60 IY = 1, 3
            JX = IPTAX(IX,2)
            JY = IPTAX(IY,2)
            CONST = EIGVAL(I)*PMASS
            GTRANN(I,J)    = GTRANN(I,J)    + EIGVEC(IX,I)
     &                   *GFANUC(JX,JY)*EIGVEC(IY,J)*CONST
            GTRANPLEE(I,J) = GTRANPLEE(I,J) + EIGVEC(IX,I)
     &                   *PMATL(JX,JY,1)*EIGVEC(IY,J)*CONST*D2
            GTRANPLEP(I,J) = GTRANPLEP(I,J) + EIGVEC(IX,I)
     &                   *PMATL(JX,JY,2)*EIGVEC(IY,J)*CONST*D2
            GTRANPSEE(I,J) = GTRANPSEE(I,J) + EIGVEC(IX,I)
     &                   *PMATS(JX,JY,1)*EIGVEC(IY,J)*CONST*D2
            GTRANPSEP(I,J) = GTRANPSEP(I,J) + EIGVEC(IX,I)
     &                   *PMATS(JX,JY,2)*EIGVEC(IY,J)*CONST*D2
 60       CONTINUE
C         Total linear response contribution to molecular rotational g-tensor
          GTRANP(I,J) = GTRANPLEE(I,J) + GTRANPLEP(I,J)
     &                   + GTRANPSEE(I,J) + GTRANPSEP(I,J)
C         Total molecular rotational g-tensor
          GTRANT(I,J) = GTRANN(I,J) + GTRANP(I,J)
 50     CONTINUE
C
C    At least one inertia tensor eigenvalue equal to zero (linear molecule):
        IF(LINMOL) THEN
          IF (IPRINT .GE. 0)
     &  WRITE(LUPRI,'((1X,A,F16.8,/),/,2(1X,A,F16.8,/),(1X,A,F16.8))')
     &  ' Total molecular rotational g-factor          :',GTRANT(3,3),
     &  '   Nuclear contribution to g-factor (g^nuc)   :',GTRANN(3,3),
     &  '   Electronic contribution to g-factor (g^LR) :',GTRANP(3,3)
C
          IF (IPRINT .GE. 1)
     &  WRITE(LUPRI,'(3(1X,A,F16.8,/),(1X,A,F16.8))')
     &  '           g^LR-L(e-e) :',GTRANPLEE(3,3),
     &  '           g^LR-S(e-e) :',GTRANPSEE(3,3),
     &  '           g^LR-L(e-p) :',GTRANPLEP(3,3),
     &  '           g^LR-S(e-p) :',GTRANPSEP(3,3)
C
C    Inertia tensor of a non-linear molecule (with and without symmetry):
        ELSE
          IF (IPRINT .GE. 0) THEN
            CALL HEADER('Total molecular rotational g-tensor',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANT(1,I),I = 1,3),
     &               'B', (GTRANT(2,I),I = 1,3),
     &               'C', (GTRANT(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANT(1,1)
     &            +GTRANT(2,2)+GTRANT(3,3))/3
            CALL HEADER('Nuclear contribution (g^nuc)',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANN(1,I),I = 1,3),
     &               'B', (GTRANN(2,I),I = 1,3),
     &               'C', (GTRANN(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANN(1,1)
     &            +GTRANN(2,2)+GTRANN(3,3))/3
            CALL HEADER('Electronic contribution (g^LR)',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANP(1,I),I = 1,3),
     &               'B', (GTRANP(2,I),I = 1,3),
     &               'C', (GTRANP(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANP(1,1)
     &            +GTRANP(2,2)+GTRANP(3,3))/3
          END IF
C
C         Print detailed linear response contributions to molecular rotational g-tensor
C
          IF (IPRINT .GE. 1) THEN
            CALL HEADER('g^LR-L(e-e)',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANPLEE(1,I),I = 1,3),
     &               'B', (GTRANPLEE(2,I),I = 1,3),
     &               'C', (GTRANPLEE(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANPLEE(1,1)
     &            +GTRANPLEE(2,2)+GTRANPLEE(3,3))/3
            CALL HEADER('g^LR-S(e-e)',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANPSEE(1,I),I = 1,3),
     &               'B', (GTRANPSEE(2,I),I = 1,3),
     &               'C', (GTRANPSEE(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANPSEE(1,1)
     &            +GTRANPSEE(2,2)+GTRANPSEE(3,3))/3
            CALL HEADER('g^LR-L(e-p)',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANPLEP(1,I),I = 1,3),
     &               'B', (GTRANPLEP(2,I),I = 1,3),
     &               'C', (GTRANPLEP(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANPLEP(1,1)
     &              +GTRANPLEP(2,2)+GTRANPLEP(3,3))/3
            CALL HEADER('g^LR-S(e-p)',-1)
            WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &               'A', (GTRANPSEP(1,I),I = 1,3),
     &               'B', (GTRANPSEP(2,I),I = 1,3),
     &               'C', (GTRANPSEP(3,I),I = 1,3)
            WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &               'iso', (GTRANPSEP(1,1)
     &            +GTRANPSEP(2,2)+GTRANPSEP(3,3))/3
          END IF
        END IF  ! (LINMOL)
      END IF    ! (LINEAR)
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prop2bss */
      SUBROUTINE PROP2BSS(VMAT,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Run the list of four-component property operators
C             in the DCBPRP block and do the picture change
C             transformation of selected. Do not transform
C             certain operators (like overlap, connection matrixes...)
C
C             However, all operators of the list are written into BSSMAT in LL-blocks,
C             either picture change transformed or not. Some of them in LL form
C            (like kinetic energy, beta matrix) are useless (MI/Nov 07).
C             TODO: select operators which do not have to be picked up in LL block.
C
C    Called from: MAKE_H2C
C      (if TWOCOMPBSS is true, ie in true two-comp. mode when IOTC Hamiltonian is on)
C
C    On input:
C              VMAT - empty field (N2BBASXQ) to be utilized inside
C              IPRINT - print level (=IPRHAM)
C
C    On output:
C              BSS/DKn integrals in LL_only matrix are written to the BSSMAT file
C              under the label P2C_xxxx.
C
C    Written by Miro Ilias, Strasbourg, March 2006
C    Modifications: Miro Ilias, Tel Aviv, Nov. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXQN
!      use dgroup_mod, only: NZ
!      use dcbbas_mod, only: NTBAS
!      use dcbgen_mod, only: LUBSS
!      use dcbxpr_mod, only: IPRPSYM, NPRPS, PRPNAM
!      use ccom_mod, only: DOSPHE
#include "implicit.h"
      INTEGER :: I
#include "dgroup.h"
!#include "dcborb.h"
#include "dcbbas.h"
!#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
!#include "dcbprp.h"
!#include "dcbprl.h"
#include "ccom.h"
C
      DIMENSION WORK(LWORK),VMAT(*)
      LOGICAL IS_PCTRA
      CHARACTER*8 PRPLABEL
      CHARACTER*4 NUMSTR
C
      CALL QENTER('PROP2BSS')
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK

C     ... allocate for LL-only operator
      N2BBASXQ_L = NTBAS(1)*NTBAS(1)*NZ
      CALL MEMGET('REAL',KPRLL,N2BBASXQ_L,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ       ,WORK,KFREE,LFREE)
      IF (IPRINT.GE.3) THEN
        CALL HEADER('*** Output from PROP2BSS ***',-1)
      ENDIF
C     ... run over all operators in the list
      DO I = 1, NPRPS
C       ... get info if this operator has to be picture change transformed
        CALL DECIDEPCTR(IS_PCTRA,PRPNAM(I),IPRINT)
        IREP = IPRPSYM(I) - 1
        IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,*)
          CALL WRIXPR(I,I)
          IF (IS_PCTRA) THEN
            WRITE(LUPRI,'(2X,A)')
     &      'Operator has to be picture change transformed.'
          ELSE
            WRITE(LUPRI,'(2X,A)')
     &      'Operator need NOT to be picture change transformed.'
          ENDIF
        ENDIF
C       ... get the operator either in the full four-component form,
C       .... or already in  the picture change transformed shape
C       ...   (occupying in fact only  LL block of the four-component matrix)
        INDXPR = I
        CALL PRPMSAO(INDXPR,VMAT,.TRUE.,WORK,WORK(KFIRST),
     &               WORK(KFREE),LFREE,IPRINT)
C       .. extract LL block from VMAT into WORK(KPRLL)
        CALL EXTRLL(IREP,VMAT,WORK(KPRLL),IPRINT)
C      ... save WORK(KPRLL) into BSSMAT under label P2C_xxxx, with xxxx=INDXPR
        CALL NUM2STR(INDXPR,NUMSTR)
        PRPLABEL = 'P2C_'//NUMSTR
caspg
        if (dosphe) then
           CALL SPH_IOTC(WORK(KPRLL),VMAT,NSPH,IPRINT,WORK(KFREE),LFREE)
           NDIM=NSPH*NSPH*NZ
           CALL W2BSSMAT(1,.TRUE.,PRPLABEL,'PROP2BSS',
     &              VMAT,NDIM,LUBSS,IPRINT)
        else
           CALL W2BSSMAT(1,.TRUE.,PRPLABEL,'PROP2BSS',
     &              WORK(KPRLL),N2BBASXQ_L,LUBSS,IPRINT)
        endif

        IF (IPRINT.GE.2) THEN
          WRITE(LUPRI,'(/2X,A,I3,A)')
     &    'PROP2BSS: 2comp. operator (indxpr=',INDXPR,') '//
     &    PRPLABEL//'->'//PRPNAM(I)//
     &    ' was written to the BSSMAT file in LL block.'
          IF (IS_PCTRA) THEN
            WRITE(LUPRI,'(2X,A)')
     &      '...operator was picture change transformed'
          ELSE
C            ...
            WRITE(LUPRI,'(2X,A)')
     &      '...operator was NOT picture change transformed'
          ENDIF
          if (dosphe) then
           WRITE(LUPRI,'(2X,A)')' A SPHERICAL representation is used.'
          endif
        ENDIF
caspg
      ENDDO
C
      CALL MEMREL('KPRLL',WORK,KWORK,KPRLL,KFREE,LFREE)
C
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'PROP2BSS: EOFLABEL was renewed after saving all 2c operators '//
     & 'into BSSMAT.'
      ENDIF

      CALL QEXIT('PROP2BSS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_OPTROT */
      SUBROUTINE DEF_OPTROT(IPRINT)
C***********************************************************************
C
C     Define operator for optical rotation; based on DEF_SHIELD
C     and DEF_EPOLAR
C     Written by T. Saue and K. Ruud Jan 8 2008
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use codata_mod, only: G1FAC
!      use dcbgen_mod, only: CVAL
!      use dcbprp_mod, only: IPEPOLAR, IPSHIELD, LONDON
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: DM1 = -1.0D0, DMI2 = -0.5D0
C
#include "codata.h"
!#include "dgroup.h"
!#include "pgroup.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
!#include "dcbexp.h"
!#include "dcbnmr.h"
      CHARACTER PNAME*16, PLABEL(3)*8,LABEL*4
      REAL*8    PFAC(3)
      INTEGER   IOP(3)
C
      CALL QENTER('DEF_OPTROT')
C
C     First define left-hand operators
C
      PNAME     = 'X-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(1,1) = INDOP1
C
      PNAME     = 'Y-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(2,1) = INDOP1
C
      PNAME     = 'Z-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('ZEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(3,1) = INDOP1
C
C     Define right hand operators
C     (PFAC(1:2) for Dirac mag. field, PFAC(3) for g_1 QED correction)
C

      IF (LONDON) THEN
        CALL DEF_D1HBLOND(1,IPRINT)
      ELSE
C
C     Conventional atomic orbitals
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
      NPCOMP = 2
CTROND: Using DALTON sign convention for now. To fix later.
      PFAC(1)  = -DMI2*CVAL
      PFAC(2)  = -DMI2*CVAL
      PFAC(3)  = -DMI2*G1FAC
C X component of B-field
      PNAME     = 'X magnetic field'
      IPTYP     = IOP(1)
      PLABEL(1) = 'ZDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('XSHIELB',IPSHIELD(1+MXCOOR),
     &             LLRBPU,NLRBPT,INDXPR,MAXLLR)
#if !defined (XH4)
C Y component of B-field
      PNAME     = 'Y magnetic field'
      IPTYP     = IOP(2)
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('YSHIELB',IPSHIELD(2+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
C Z component of B-field
      PNAME     = 'Z magnetic field'
      IPTYP     = IOP(3)
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
      CALL OP1IND('ZSHIELB',IPSHIELD(3+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
#endif
      ENDIF
C
C     Print section for optical rotation:
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' Optical rotation :'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     A operators
C     ===========
C
       CALL PRSYMB(LUPRI,'-',75,0)
       WRITE(LUPRI,'(1X,A)')
     &'** A operators for optical rotation **'
       CALL PRSYMB(LUPRI,'-',75,0)
       DO I = 1, 3
        INDOP1 = IPEPOLAR(I,1)
        INDXPR = LLRAPU(INDOP1)
        CALL WRIXPR(INDOP1,INDXPR)
       ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(A)') ' ** B operators optical rotation **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3
        INDOP1 = IPSHIELD(I+MXCOOR)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO

      CALL QEXIT('DEF_OPTROT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_optrot */
      SUBROUTINE PRP_OPTROT(ABLRF,ABCNV)
C***********************************************************************
C
C     Output routine for optical rotation, based on LNROUT from DALTON
C     Written by T. Saue in collaboration with K. Ruud Feb 2008
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: IPEPOLAR, IPSHIELD, LONDON
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT
!      use codata_mod, only: XFMOL, XTANG, XTKAYS, XTNM
      use pi_mod, only: PI
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: D0=0.0D0, D1=1.0D0, D3=3.0D0
      DOUBLE PRECISION FACTOR
C
#include "dcbprp.h"
#include "dcbxlr.h"
#include "codata.h"
      REAL*8    ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),TMPMAT(3,3)
      EXTERNAL TOTMASS
C
      FACTOR = (288.0D-30)*(PI**2)*XFMOL*(XTANG**4)
C
      IF(LONDON) THEN
        CALL AROUND('Optical rotation (London)')
      ELSE
        CALL AROUND('Optical rotation (no London)')
      ENDIF
      TMASS = TOTMASS()
      WRITE(6,*) '* Molecular mass:',TMASS,FACTOR
      FAC = FACTOR/TMASS
      WRITE(LUPRI,*)
      DO IFREQ = 1,NBFREQ
        FREQ = BFREQ(IFREQ)*XTKAYS
        WRITE(LUPRI,'(/A,3(F12.6,A))')
     &        '* Frequency: ',BFREQ(IFREQ),' au  = ',
     &         XTNM/BFREQ(IFREQ), ' nm = ',FREQ,' cm-1'
        WRITE(LUPRI,*)
C
C       Extract mixed electric dipole-magnetic dipole polarizability
C
        DO I2 = 1,3
          IB = IPSHIELD(I2+MXCOOR)
          DO I1 = 1,3
            IE = IPEPOLAR(I1,1)
            TMPMAT(IE,IB) = -ABLRF(IFREQ,IE,IB,1)
          ENDDO
        ENDDO
        WRITE(LUPRI,'(/10X,A/)')
     &    'G tensor:'
        WRITE(LUPRI,9009) 'Bx','By','Bz'
        WRITE(LUPRI,9008) 'Ex',(TMPMAT(1,I), I=1,3)
        WRITE(LUPRI,9008) 'Ey',(TMPMAT(2,I), I=1,3)
        WRITE(LUPRI,9008) 'Ez',(TMPMAT(3,I), I=1,3)
        OR_BETA = (TMPMAT(1,1)+TMPMAT(2,2)+TMPMAT(3,3))/D3/BFREQ(IFREQ)
        OR_ALPHA = FAC*FREQ*FREQ*OR_BETA
        WRITE(LUPRI,'(4(/10X,A,F12.6))')
     &    'Beta             :',OR_BETA,
     &    'Optical rotation :',OR_ALPHA
      ENDDO
 9008 FORMAT (16X,A,2X,3F12.7)
 9009 FORMAT (27X,3(A,10X)/)

      RETURN
C
      end

      subroutine def_blao(work, lwork, iprint)
!     ------------------------------------------------------------------
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: BXLAO, BYLAO, BZLAO, DOEXP, IPLONDON, IP_BLAO
#include "implicit.h"
      INTEGER :: J
      PARAMETER(D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0,
     &          DI2 = +0.5D00,DM1 = - 1.0D0)
!#include "codata.h"
!#include "dgroup.h"
!#include "pgroup.h"
!#include "dcbgen.h"
#include "dcbprp.h"
!#include "dcbesr.h"
!#include "dcbxlr.h"
!#include "dcbexp.h"
!#include "dcbnmr.h"
      CHARACTER PNAME*16, PLABEL(3)*8,LABEL*4,chcomp*1
      DIMENSION WORK(LWORK)
!#include "dcbxpr.h"
!#include "dcbsusc.h"
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK

!     ------------------------------------------------------------------

      CALL HEADER(
     &     'DEF_BLAO: Operators in case the perturbation is B_i' //
     &     ' with LAO ',-1)

      CALL IZERO(IP_BLAO, 6)
      CALL IZERO(IPLONDON, 12)

      doexp = .true.
      call def_d1hblond(3, work, lwork, iprint)


      if(bxlao) then
        icomp  = 1
        chcomp = 'X'
      else if(bylao) then
        icomp  = 2
        chcomp = 'Y'
      else if(bzlao) then
        icomp  = 3
        chcomp = 'Z'
      else
        call quit('programming error in def_blao')
      endif

!     print section:
!     ==============
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,1X,A)') ' BLAO, component :', chcomp
      CALL PRSYMB(LUPRI,'=',75,0)
!
!     A operators
!     ===========
!
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A,1X,A)')
     &      '** A operators for blao, component **', chcomp
      CALL PRSYMB(LUPRI,'-',75,0)
      DO J = 1, 4
        CALL WRIXPR(0,IPLONDON(J,ICOMP))
      ENDDO
!
!     B operators
!     ===========
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A,1X,A)')
     &      '** B operators for blao, component **', chcomp
      CALL PRSYMB(LUPRI,'-',75,0)
      DO J = 1, 4
        CALL WRIXPR(0,IPLONDON(J,icomp))
      ENDDO


      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_EFF_DEN2(PNAME,PLABEL,PFAC,IPRINT)
C**********************************************************
C     Define effective density operator
C     Written by R. van Meer September 2010
C**********************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"

      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     (-99: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-99)
      CALL DEF_EFF_DEN21(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_EFF_DEN21(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C**********************************************************
C     Define effective density operator
C     Written by R. van Meer September 2010
C**********************************************************
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: ISTBNU, NAMDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTNUC, MAXREP
      use chrnos_mod, only: CHRNOS
!      use dcbexp_mod, only: LEXPP, MAXEXP, NEXPP
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND)
C
!#include "dcbxpr.h"
#include "dcbexp.h"
C
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      DO 100 IREP = 0, MAXREP
        DO 200 IATOM = 1, NUCIND
          IF (DOATOM(IATOM)) THEN
            IF (IAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
              PLABEL(1) = 'ED '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                          //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                          //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              PNAME = 'Eff._Dens '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                              //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                              //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                     INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND('EFFDEN',IDUM,LEXPP,NEXPP,INDXPR,MAXEXP)
            END IF
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_EFFDEN */
      SUBROUTINE DEF_EFFDEN(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operator for effective density at nuclei
C     see S. Knecht, S. Fux, R. van Meer, L. Visscher, M. Reiher
C     and T. Saue, Theor. Chem. Acc. 129 (2011) 631
C     Based on DEF_PVC
C
C     Written by T. Saue July 25 2012
C
C***********************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCIND
#include "implicit.h"
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      logical, allocatable :: DOATOM(:)
C
      allocate(DOATOM(NUCIND))
C     PVC: INTTYP = 62
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)
      CALL DEF_EFFDEN1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
      deallocate(DOATOM)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_effden1 */
      SUBROUTINE DEF_EFFDEN1(PNAME,PLABEL,PFAC,DOATOM,IPRINT)
C***********************************************************************
C
C     Define operator for parity violation - chirality
C
C     Written by J. Thyssen - Oct 7 1999
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbxpr_mod, only: PRPNAM
!      use dcbexp_mod, only: IPEFFDEN, LEXPP, MAXEXP, NEXPP
      use nuclei_mod, only: GAUNUC, ISTBNU, NAMDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTNUC, MAXREP
      use chrnos_mod, only: CHRNOS
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
!#include "codata.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0, D3 = 3.0D0, D4 = 4.0D0
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3)
      LOGICAL   DOATOM(NUCIND)
#include "dcbxpr.h"
#include "dcbexp.h"
C
C
      ITYP = 0
      IPTYP     = 1
      NPCOMP    = 1
      IF(GAUNUC) THEN
        PFAC(1) = D1
      ELSE
        WRITE(6,*) 'WARNING: Effective density at nucleus is not'//
     &       'well defined for point nucleus,',
     &       'due to singularity in the relativistic case. Skipping.'
        RETURN
      ENDIF
      DO 100 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 200 IREP = 0, MAXREP
               IF (IAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
                  ITYP = ITYP + 1
                  PNAME = 'EFFD:'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                    //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                    CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,IPRINT)
                  CALL OP1IND('PVC   '
     &                 ,IPEFFDEN(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
C
      IF(IPRINT.GE.1) THEN
        DO I = 1,ITYP
          INDXPR = LEXPP(IPEFFDEN(I))
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_EPOLE(L,IPEPOL,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Define electric Cartesian multipole operators of order L
C     For given L there are (L+2)(L+1)/2 linearly
C     independent Cartesian operators specified by powers (NX,NY,NZ)
C     of the Cartesian coordinates. They can therefore be generated by
C     looping over the lower triangle of a square matrix of dimension L+1.
C     The loop below is inspired by the routine CARPOW, but simplified
C     by starting rows and columns by index zero. From an element (a,b) of
C     such a shifted matrix the values of the power triplet is given by
C
C     NX = L-a, NY =  a - b , NZ = b
C
C     Note that this routine does not check if IPEPOL goes out of
C     bounds. This has to be done outside the routine.
C
C     Written by Trond Saue March 29 2013
C***********************************************************************
      use chrnos_mod, only: CHRNOS
#include "implicit.h"
      INTEGER :: II
      DIMENSION IPEPOL(*)
      CHARACTER PNAME*16, PLABEL*8
C
      IPTYP       = 1
      NPCOMP      = 1
      PFAC        = -1.0D0
      PLABEL(1:2) = 'CM'
      PNAME       = 'EpolL'//CHRNOS(L/10)//CHRNOS(MOD(L,10))
     &              //'X00Y00Z00'
      II = 0
      DO IA = 0,L
        NX = L - IA
        DO IB = 0,IA
          II = II + 1
          NY = IA - IB
          NZ = IB
          PLABEL(3:8) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
          PNAME( 9:10) = PLABEL(3:4)
          PNAME(12:13) = PLABEL(5:6)
          PNAME(15:16) = PLABEL(7:8)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPEPOL(II),ISYXPR,ITRXPR,IPRINT)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_VPOLE(L,IPVPOL,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Define electric multipole moment in the generalized velocity representation.
C     Multipole moments of this type are given by
C     
C     Q^[n+1]_(j1...jn;p) = -e/omega*r_j1...r_(jn-1)(r_jn*(ca_p) + n*(ca_jn)r_p)
C
C     where r_i represent the ith component of a Cartesian coordinate and ca_i
C     the ith alpha matrix multiplied with the speed of light.
C     For a given combination of (j1...j_(n-1)), we can generate nine unique
C     multipole moments
C
C     1: (n+1)x*ca_x               2: (n+1)y*ca_y            3: (n+1)z*ca_z 
C     4: (z*(ca_y) + n*(ca_z)y)    5: (y*(ca_z) + n*(ca_y)z) 6: (x*(ca_z) + n*(ca_x)z) 
C     7: (z*(ca_x) + n*(ca_z)x)    8: (y*(ca_x) + n*(ca_y)x) 9: (x*(ca_y) + n*(ca_x)y)
C
C     Note that for a given set of (j1...j_(n-1)), the multipole moments are stored in the
C     order above. To each of these terms, we associate the following label
C
C     1: xx    2: yy    3: zz
C     4: zy    5: yz    6: xz 
C     7: zx    8: yx    9: xy
C
C     For example, a multipole moment of the form 
C
C     Q^[5]_(xyyz;x) = -e/omega*xy^2(z*(ca_x) + n*(ca_z)x)
C
C     gets the following name 
C
C     VpolzxX01Y02Z00
C
C     This operator is symmetric with respect to the indices (j1...j_(n-1)),
C     thus the possible combinations are stored in lower triangular form.
C     Therefore, at order L, we have 9*(L-1)L/2 unique multipoles.
C     
C     However, L=1 is an exception to this rule.
C     At this order, the multipole operator reduces to the set (ca_x,ca_y,ca_z),
C
C     Note that this routine does not check if IPVPOL goes out of 
C     bounds. This has to be done outside the routine.
C
C     Written by Martin van Horn July 16 2021
C***********************************************************************
      use chrnos_mod, only: CHRNOS
!      use dcbgen_mod, only: CVAL
#include "implicit.h"
      INTEGER :: I, II
#include "dcbgen.h"
      DIMENSION IPVPOL(*), PFAC(2)
      REAL CARPOW(3,3)
      CHARACTER PNAME*16, PLABEL(2)*8,COMPLAB(3)*2  
C
      II = 0
      PFAC  = CVAL
      PNAME       = 'VpolxxX00Y00Z00'
      PLABEL = 'CM000000'
      IF (L.EQ.1) THEN ! for L=1 the set of multipoles reduces to (ca_x,ca_y,ca_z)
            COMPLAB = ['x1','y1','z1']
            NPCOMP = 1
            IPTYP = 2
            DO IP = 1, 3
                  PNAME(6:7) = COMPLAB(IP)
                  II = II + 1 
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC(1),PLABEL(1),
     &            IPVPOL(II),ISYXPR,ITRXPR,IPRINT)
                  IPTYP = IPTYP + 1
            ENDDO
      ELSE 
            IPTYP = 2
            COMPLAB = ['xx','yy','zz']
            CARPOW = 0.0D0    
            DO I = 1,3
                  CARPOW(I,I) = 1
            ENDDO
            DO IA = 0,L-2   
                  IX = L - 2 - IA
                  DO IB = 0,IA     !Loops over lower triangle
                        IY = IA - IB
                        IZ = IB
                        IPTYP = 2
                        DO IP = 1,3  !define multipole moments xx yy and zz
                              PNAME(5:6) = COMPLAB(IP)                                           
                              PNAME( 8:9) = CHRNOS(IX/10)
     &                                     //CHRNOS(MOD(IX,10))
                              PNAME(11:12) = CHRNOS(IY/10)
     &                                     //CHRNOS(MOD(IY,10))
                              PNAME(14:15) = CHRNOS(IZ/10)
     &                                     //CHRNOS(MOD(IZ,10))     
                              NPCOMP = 1
                              PFAC(1) = L*CVAL
                              NX = IX + CARPOW(1,IP)
                              NY = IY + CARPOW(2,IP)
                              NZ = IZ + CARPOW(3,IP)
                              PLABEL(1)(3:8) = CHRNOS(NX/10)
     &                        //CHRNOS(MOD(NX,10))//CHRNOS(NY/10)
     &                        //CHRNOS(MOD(NY,10))//CHRNOS(NZ/10)
     &                        //CHRNOS(MOD(NZ,10))
                              II = II + 1
                              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC(1)
     &                       ,PLABEL(1),IPVPOL(II),ISYXPR,ITRXPR,IPRINT)          
                              IPTYP = IPTYP + 1
                        ENDDO
                        COMPLAB = ['zy','xz','yx']
                        DO IP = 1,3    !define multipole moments zy yz xz zx yx xy
                              NPCOMP = 2
                              PFAC(1) = CVAL 
                              PFAC(2) = -(L-1)*CVAL
                              PNAME(5:6) = COMPLAB(IP)
                              PNAME( 8:9) = CHRNOS(IX/10)
     &                                     //CHRNOS(MOD(IX,10))
                              PNAME(11:12) = CHRNOS(IY/10)
     &                                     //CHRNOS(MOD(IY,10))
                              PNAME(14:15) = CHRNOS(IZ/10)
     &                                     //CHRNOS(MOD(IZ,10))    
                              NX = IX + CARPOW(IP,2)
                              NY = IY + CARPOW(IP,3)
                              NZ = IZ + CARPOW(IP,1)
                              PLABEL(1)(3:8) = CHRNOS(NX/10)
     &                               //CHRNOS(MOD(NX,10))//CHRNOS(NY/10)
     &                               //CHRNOS(MOD(NY,10))//CHRNOS(NZ/10)
     &                               //CHRNOS(MOD(NZ,10))                                     
                              NX = IX + CARPOW(3,IP)
                              NY = IY + CARPOW(1,IP)
                              NZ = IZ + CARPOW(2,IP)
                              PLABEL(2)(3:8) = CHRNOS(NX/10)
     &                        //CHRNOS(MOD(NX,10))//CHRNOS(NY/10)
     &                        //CHRNOS(MOD(NY,10))//CHRNOS(NZ/10)
     &                        //CHRNOS(MOD(NZ,10))
                              II = II + 1
                              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL
     &                           ,IPVPOL(II),ISYXPR,ITRXPR,IPRINT)    
                              PFAC(1) = (L-1)*CVAL 
                              PFAC(2) = -CVAL
                              PNAME(5:6) = COMPLAB(IP)(2:2)
     &                                   //COMPLAB(IP)(1:1)
                              II = II + 1
                              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC
     &                          ,PLABEL,IPVPOL(II),ISYXPR,ITRXPR,IPRINT)
                              IPTYP = IPTYP + 1
                        ENDDO
                  ENDDO
            ENDDO
      ENDIF
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_MPOLE(L,IPMPOL,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Define magnetic Cartesian multipole operators of order L
C
C     The general formula is
C
C       L/(L+1)*x^(NX)*y^(NY)*z^(NZ)*(RVEC times JVEC)_q
C
C     where the Cartesian powers (NX,NY,NZ) sum to L-1,
C     RVEC is the position vector and JVEC the operator for the
C     current density vector. In the relativistic domain the current
C     density vector is -c*ALPHAVEC, where the sign comes from
C     the electron charge, c is the speed of light and ALPHAVEC is the
C     vector of Dirac alpha matrices.
C
C     A magnetic multipole of order L is symmetric in all indices
C     except q, so there are 3(L+1)L/2 linearly independent Cartesian operators
C     specified by the powers (NX,NY,NZ) of the Cartesian coordinates.
C     They can therefore be generated by looping over the lower triangle
C     of a square matrix of dimension L.
C     The loop below is inspired by the routine CARPOW, but simplified
C     by starting rows and columns by index zero. Let LM = L - 1.
C     From an element (a,b) of such a shifted matrix the values of the
C     power triplet is given by
C
C     NX = LM-a, NY =  a - b , NZ = b
C
C     Note that this routine does not check if IPMPOL goes out of
C     bounds. This has to be done outside the routine.
C
C     Written by Trond Saue March 29 2013
C***********************************************************************
      use chrnos_mod, only: CHRNOS
!      use dcbgen_mod, only: CVAL
#include "implicit.h"
      INTEGER :: II
#include "dcbgen.h"
      DIMENSION IPMPOL(*),PFAC(2)
      CHARACTER PNAME*16, PLABEL(2)*8
C
      NPCOMP      = 2
      PFAC(1)     = CVAL*L/(L+1)
      PFAC(2)     = PFAC(1)
      PLABEL(1)(1:2) = 'CM'
      PLABEL(2)(1:2) = 'CM'
      PNAME       = 'MpolL'//CHRNOS(L/10)//CHRNOS(MOD(L,10))
     &              //'X00Y00Z00'
      LM = L-1
      II = 0
      DO IA = 0,LM
        NX = LM - IA
        DO IB = 0,IA
          NY = IA - IB
          NZ = IB
          PNAME( 9:10) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
          PNAME(12:13) = CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
          PNAME(15:16) = CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
C.........X-component
          II = II + 1
          IPTYP = 5 ! XAVECTOR
          PNAME(5:5) = 'X'
          NQ = NZ+1
          PLABEL(1)(3:8) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                   //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                   //CHRNOS(NQ/10)//CHRNOS(MOD(NQ,10))
          NQ = NY+1
          PLABEL(2)(3:8) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                   //CHRNOS(NQ/10)//CHRNOS(MOD(NQ,10))
     &                   //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPMPOL(II),ISYXPR,ITRXPR,IPRINT)
C.........Y-component
          II = II + 1
          IPTYP = 6 ! YAVECTOR
          PNAME(5:5) = 'Y'
          NQ = NX+1
          PLABEL(1)(3:8) = CHRNOS(NQ/10)//CHRNOS(MOD(NQ,10))
     &                   //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                   //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
          NQ = NZ+1
          PLABEL(2)(3:8) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                   //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                   //CHRNOS(NQ/10)//CHRNOS(MOD(NQ,10))
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPMPOL(II),ISYXPR,ITRXPR,IPRINT)
C.........Z-component
          II = II + 1
          IPTYP = 7 ! ZAVECTOR
          PNAME(5:5) = 'Z'
          NQ = NY+1
          PLABEL(1)(3:8) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                   //CHRNOS(NQ/10)//CHRNOS(MOD(NQ,10))
     &                   //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
          NQ = NX+1
          PLABEL(2)(3:8) = CHRNOS(NQ/10)//CHRNOS(MOD(NQ,10))
     &                   //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                   //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPMPOL(II),ISYXPR,ITRXPR,IPRINT)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_BEDCON(L,IPBEDCON,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Define contribution to order L in the wave vector from the
C     full light-matter interaction
C
C     The general formula is
C
C       x^(NX)*y^(NY)*z^(NZ)*(c*ALPHAVEC)_p/(L!)
C
C     where the Cartesian powers (NX,NY,NZ) sum to L,
C     c is the speed of light and ALPHAVEC is the
C     vector of Dirac alpha matrices.
C
C     A contribution of order L is symmetric in all indices
C     except p, so there are 3(L+2)(L+1)/2 linearly independent Cartesian operators
C     specified by the powers (NX,NY,NZ) of the Cartesian coordinates.
C     They can therefore be generated by looping over the lower triangle
C     of a square matrix of dimension L+1.
C     The loop below is inspired by the routine CARPOW, but simplified
C     by starting rows and columns by index zero. From an element (a,b) of
C     such a shifted matrix the values of the power triplet is given by
C
C     NX = L-a, NY =  a - b , NZ = b
C
C     Note that this routine does not check if IPBEDCON goes out of
C     bounds. This has to be done outside the routine.
C
C     Written by Trond Saue October 5 2018
C***********************************************************************
      use chrnos_mod, only: CHRNOS
!      use dcbgen_mod, only: CVAL
#include "implicit.h"
      INTEGER :: I, II
#include "dcbgen.h"
      DIMENSION IPBEDCON(*)
      CHARACTER PNAME*16, PLABEL*8
C
C      IF(MOD(L,2).EQ.0) THEN
C         K = L/2
C      ELSE
C         K = (L-1)/2
C      ENDIF
C      PFAC = CVAL*(-1.0D0)**K
      PFAC = CVAL
      DO I = 2,L
        PFAC = PFAC/I
      ENDDO
      NPCOMP         = 1
      PLABEL(1:2) = 'CM'
      PNAME       = 'BEDcL'//CHRNOS(L/10)//CHRNOS(MOD(L,10))
     &              //'X00Y00Z00'
      II = 0
      DO IA = 0,L
        NX = L - IA
        DO IB = 0,IA
          NY = IA - IB
          NZ = IB
          PLABEL(3:8) = CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
          PNAME( 9:10) = PLABEL(3:4)
          PNAME(12:13) = PLABEL(5:6)
          PNAME(15:16) = PLABEL(7:8)
C.........X-component
          II = II + 1
          IPTYP = 2 ! XALPHA
          PNAME(5:5) = 'X'
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPBEDCON(II),ISYXPR,ITRXPR,IPRINT)
C.........Y-component
          II = II + 1
          IPTYP = 3 ! YALPHA
          PNAME(5:5) = 'Y'
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPBEDCON(II),ISYXPR,ITRXPR,IPRINT)
C.........Z-component
          II = II + 1
          IPTYP = 4 ! ZALPHA
          PNAME(5:5) = 'Z'
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            IPBEDCON(II),ISYXPR,ITRXPR,IPRINT)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_efg */
      SUBROUTINE NUC_EFG(IPRINT,DOATOM)
C***********************************************************************
C
C     Calculates nuclear contribution to electric field gradients
C
C     Written by T.Enevoldsen  - Sep 1996
C     Last revision by W.A. de Jong - April 1999
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CHARGE, CORD, ISTBNU, NAMDEP, NUCDEG, NUCIND
      use symmet_mod, only: ISYMAX, MAXOPR, PT
!      use nqcc_mod, only: EFGNUC
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: KK
      REAL(dp), PARAMETER :: D0=0.0D0, D3 = 3.0D0
      DIMENSION  GEOM(3*8*MXCENT)
      LOGICAL    DOATOM(*)
#include "nqcc.h"
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
C
      CALL DZERO(EFGNUC,9*MXCENT)
      KK = 1
      NATOM = 0
      IATOMD = 0
      DO 10 I1 = 1, NUCIND
      IF(DOATOM(I1)) THEN
        DO 20 IDEG1 = 1,NUCDEG(I1)
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          LL = 1
          DO 100 I2 = 1, NUCIND
          IF (CHARGE(I2).EQ.D0) GOTO 100
            DO 110 IDEG2=1,NUCDEG(I2)
              IF (LL .NE. KK) THEN
                XCOOR = GEOM(LL) - GEOM(KK)
                YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                R5 = R2*R2*SQRT(R2)
                EFGNUC(1,1,NATOM) = EFGNUC(1,1,NATOM)
     &                            + CHARGE(I2)*(D3*XCOOR*XCOOR - R2)/R5
                EFGNUC(2,2,NATOM) = EFGNUC(2,2,NATOM)
     &                            + CHARGE(I2)*(D3*YCOOR*YCOOR - R2)/R5
                EFGNUC(3,3,NATOM) = EFGNUC(3,3,NATOM)
     &                            + CHARGE(I2)*(D3*ZCOOR*ZCOOR - R2)/R5
                EFGNUC(1,2,NATOM) = EFGNUC(1,2,NATOM)
     &                            + CHARGE(I2)*D3*XCOOR*YCOOR/R5
                EFGNUC(1,3,NATOM) = EFGNUC(1,3,NATOM)
     &                            + CHARGE(I2)*D3*XCOOR*ZCOOR/R5
                EFGNUC(2,3,NATOM) = EFGNUC(2,3,NATOM)
     &                            + CHARGE(I2)*D3*YCOOR*ZCOOR/R5
              END IF
              LL = LL + 3
 110        CONTINUE
 100      CONTINUE
          EFGNUC(2,1,NATOM) = EFGNUC(1,2,NATOM)
          EFGNUC(3,1,NATOM) = EFGNUC(1,3,NATOM)
          EFGNUC(3,2,NATOM) = EFGNUC(2,3,NATOM)
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,2A,/)')
     &             'Nuclear EFG Tensor of nucleus: ',NAMDEP(IATOMD)
            CALL OUTPUT(EFGNUC(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
          END IF
          KK = KK + 3
 20     CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(I1)
        KK = KK + 3*NUCDEG(I1)
      ENDIF
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_nef */
      SUBROUTINE NUC_EFN(IPRINT,DOATOM,EFNNUC,GEOM)
C***********************************************************************
C
C     Calculates nuclear contribution to the electric field at
C     individual nuclei
C
C     Based on NUC_EFG
C     Written by T.Saue  - Feb 2003
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use consts_mod, only: D0
      use nuclei_mod, only: CHARGE, CORD, ISTBNU, NAMDEP, NUCDEG,
     &  NUCDEP, NUCIND
      use symmet_mod, only: ISYMAX, MAXOPR, PT
#include "implicit.h"
      INTEGER :: K, KK
      LOGICAL    DOATOM(*)
      DIMENSION  GEOM(3*NUCDEP),EFNNUC(3,NUCDEP)
C
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
      CALL DZERO(EFNNUC,3*NUCDEP)
      KK = 1
      NATOM  = 0
      IATOMD = 0
      DO 10 I1 = 1, NUCIND
      IF(DOATOM(I1)) THEN
        DO 20 IDEG1 = 1,NUCDEG(I1)
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          LL = 1
          DO 100 I2 = 1, NUCIND
            IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
            DO 110 IDEG2 = 1,NUCDEG(I2)
              IF (LL .NE. KK) THEN
                XCOOR = GEOM(LL) - GEOM(KK)
                YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                R3 = R2*SQRT(R2)
                EFNNUC(1,NATOM) = EFNNUC(1,NATOM)- CHARGE(I2)*XCOOR/R3
                EFNNUC(2,NATOM) = EFNNUC(2,NATOM)- CHARGE(I2)*YCOOR/R3
                EFNNUC(3,NATOM) = EFNNUC(3,NATOM)- CHARGE(I2)*ZCOOR/R3
              END IF
              LL = LL + 3
 110        CONTINUE
 100      CONTINUE
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,2A,/)')
     &        'Nuclear electric field at nucleus: ',NAMDEP(IATOMD)
            IOFF = ICHAR('w')
            DO K = 1,3
               WRITE(LUPRI,'(4X,A1,A,E20.10)') CHAR(IOFF+K),
     &              '-component: ',EFNNUC(K,NATOM)
            ENDDO
          END IF
          KK = KK + 3
 20     CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(I1)
        KK = KK + 3*NUCDEG(I1)
      ENDIF
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_BEDCOMP(IBED,NCOMP,IPRINT)
C***********************************************************************
C
C     Define integrals for full light-matter interaction
C
C***********************************************************************
!      use dgroup_mod, only: MCMP
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: II
C
      REAL(dp), PARAMETER :: DP5=0.5D0
!#include "dcbgen.h"
#include "dgroup.h"
      CHARACTER PNAME*16
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      DIMENSION IBED(*)
C
C
C     Allocate placeholders for integral labels and irreps
C
      allocate(LABINT(8))
      allocate(INTREP(8))
C
C     Get number of scalar integrals (NCOMP) and integral labels (and symmetries)
C     NCOMP is an even number; first half is real part, then imaginary part
C
      CALL CEXPINT_INFO(NCOMP,INTREP,LABINT)
      II = 0
C
C     x-components
C
      IPTYP  = 2
      NPCOMP = MCMP(IPTYP)
      PFAC   = DP5
      DO ICOMP = 1,NCOMP
          II = II + 1
          PNAME='BEDintX_'//LABINT(ICOMP)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,LABINT(ICOMP),
     &            IBED(II),ISYXPR,ITRXPR,IPRINT)
      ENDDO
C
C     y-components
C
      IPTYP  = 3
      NPCOMP = MCMP(IPTYP)
      PFAC   = DP5
      DO ICOMP = 1,NCOMP
          II = II + 1
          PNAME='BEDintY_'//LABINT(ICOMP)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,LABINT(ICOMP),
     &            IBED(II),ISYXPR,ITRXPR,IPRINT)
      ENDDO
C
C     z-components
C
      IPTYP  = 4
      NPCOMP = MCMP(IPTYP)
      PFAC   = DP5
      DO ICOMP = 1,NCOMP
          II = II + 1
          PNAME='BEDintZ_'//LABINT(ICOMP)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,LABINT(ICOMP),
     &            IBED(II),ISYXPR,ITRXPR,IPRINT)
      ENDDO
C     Memory deallocation
      deallocate(LABINT)
      deallocate(INTREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_BED(IBED,NCOMP,UPOL,IPRINT)
C***********************************************************************
C
C     Define integrals for full light-matter interaction
C
C***********************************************************************
!      use dgroup_mod, only: MCMP
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: II
C
      REAL(dp), PARAMETER :: DP5=0.5D0
!#include "dcbgen.h"
#include "dgroup.h"
      CHARACTER PNAME*16
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      DIMENSION IBED(*)
C
C
C     Allocate placeholders for integral labels and irreps
C
      allocate(LABINT(8))
      allocate(INTREP(8))
C
C     Get number of scalar integrals (NCOMP) and integral labels (and symmetries)
C     NCOMP is an even number; first half is real part, then imaginary part
C
      CALL CEXPINT_INFO(NCOMP,INTREP,LABINT)
      II = 0
C
C     x-components
C
      IPTYP  = 2
      NPCOMP = MCMP(IPTYP)
      PFAC   = DP5
      DO ICOMP = 1,NCOMP
          II = II + 1
          PNAME='BEDintX_'//LABINT(ICOMP)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,LABINT(ICOMP),
     &            IBED(II),ISYXPR,ITRXPR,IPRINT)
      ENDDO
C
C     y-components
C
      IPTYP  = 3
      NPCOMP = MCMP(IPTYP)
      PFAC   = DP5
      DO ICOMP = 1,NCOMP
          II = II + 1
          PNAME='BEDintY_'//LABINT(ICOMP)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,LABINT(ICOMP),
     &            IBED(II),ISYXPR,ITRXPR,IPRINT)
      ENDDO
C
C     z-components
C
      IPTYP  = 4
      NPCOMP = MCMP(IPTYP)
      PFAC   = DP5
      DO ICOMP = 1,NCOMP
          II = II + 1
          PNAME='BEDintZ_'//LABINT(ICOMP)
          CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,LABINT(ICOMP),
     &            IBED(II),ISYXPR,ITRXPR,IPRINT)
      ENDDO
C     Memory deallocation
      deallocate(LABINT)
      deallocate(INTREP)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_spinrot */
      SUBROUTINE DEF_SPINRO(PNAME,PLABEL,PFAC,IPRINT)
C***********************************************************************
C
C     Define operators for nuclear spin-rotation constants
C
C     Agustin Aucar and Trond Saue - September 2013
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: NAMN, NUCDEP, NUCIND
      use symmet_mod, only: IPTCNT, ISYMAX, MAXREP
!      use dgroup_mod, only: LINEAR
!      use pgroup_mod, only: REP
      use chrxyz_mod, only: CHRXYZ
!      use dcbexp_mod, only: IPEFN, LEXPP, MAXEXP, NEXPP
!      use dcbprp_mod, only: IPSPINRO
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
!      use dcbgen_mod, only: CVAL
!      use spinro_mod, only: LINDET
      use consts_mod, only: D0, D1
      use chrnos_mod, only: CHRNOS
      use memory_allocator
C
      implicit none
C
      integer I,J
C
#include "dgroup.h"
#include "pgroup.h"
#include "dcbexp.h"
#include "dcbprp.h"   ! common block for expectation values
#include "dcbxlr.h"   ! common block for linear reponse
#include "dcbgen.h"
#include "spinro.h"
C
      INTEGER IPRINT,NATOM,IATOM,IATOM1,
     &        IC1,IC2,ICOOR,ICOOR1,IFIRST,ILOW,INDOP1,
     &        INDXPR,IOP,IPTYP,IREP,IREP1,IREP2,IREPO,IS1,
     &        IS2,ISCOOR,ISECND,NSCOOR,ISCOR1,
     &        ISYXPR,ITHIRD,ITRXPR,ITYP,NPCOMP
C
      DOUBLE PRECISION PFAC,DL,DS
C
      logical, allocatable :: DOATOM(:)
C
      CHARACTER PNAME*16, PLABEL(2)*8
      DIMENSION PFAC(2),IOP(3)
C     ANGMOM integrals correspond to -iL, so we set DL = -1
C     XSIGMA,YSIGMA and ZSIGMA come with imaginary i, so we set DS positive (DS = 0.5)
      PARAMETER(DL = -1.0D0,DS = 0.5D0)
C
      CALL QENTER('DEF_SPINRO')
C
      NSCOOR = 3*NUCDEP
C
      CALL IZERO(IPEFN,MXCOOR)
      CALL IZERO(IPSPINRO,6+MXCOOR)
C
      allocate(DOATOM(NUCIND))
C
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
C
C     Here we set up the operators needed for the expectation value
C
      ITYP    = 0
      IPTYP   = 1  ! diagonal operator
      NPCOMP  = 1  ! number of components
      PFAC(1) = D1 ! factor
      ILOW    = 1
      IF (LINEAR) ILOW = 3
C
C     Irreps
C
      DO 100 IREP = 0, MAXREP     ! loop over irreps to generate full list of nuclei
C
C       Atoms
C
        DO 200 IATOM = 1, NUCIND  ! loop over symmetry-independent nuclei
C
C         Cartesian directions
C
          DO 500 ICOOR = ILOW, 3  ! loop over coordinates
            ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
            IF(ISCOOR .GT. 0) THEN
              ITYP   = ITYP + 1
              IFIRST = ISCOOR/100
              ISECND = MOD(ISCOOR,100)/10
              ITHIRD = MOD(MOD(ISCOOR,100),10)
              PLABEL(1) = 'NEF '//CHRNOS(IFIRST)
     &                    //CHRNOS(ISECND)//CHRNOS(ITHIRD)//' '
              PNAME     = 'NEF '//CHRNOS(IFIRST)
     &                    //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                    //NAMN(IATOM)//REP(IREP)//CHRXYZ(ICOOR)
        IF (DOATOM(IATOM)) THEN   ! selection of atoms
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                    INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND(CHRXYZ(ICOOR)//'EFN',
     &                    IPEFN(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
        ELSE
              IPEFN(ITYP) = -1
        END IF
            END IF
 500      CONTINUE
 200    CONTINUE
 100  CONTINUE
C
C     Here we set up the operators needed for the linear reponse part
C
C     Define left-hand operators: (-c) i.alpha x r/r^3
C     ================================================
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
C
      PFAC(1)  = -CVAL
      PFAC(2)  = -CVAL
      NPCOMP   = 2
C
C     Irreps
C
      DO IREPO = 0, MAXREP
C
C       Atoms
C
        DO IATOM1 = 1, NUCIND
        IF(DOATOM(IATOM1)) THEN
C
C         Cartesian directions
C
          DO ICOOR1 = 1, 3
            ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
            IF (ISCOR1 .GE. 1000) THEN
              CALL QUIT(
     &'DEF_SPINRO: Max 1000 coords in this implementation !')
            ENDIF
            IF (ISCOR1 .GT. 0) THEN
              IFIRST = ISCOR1/100
              ISECND = MOD(ISCOR1,100)/10
              ITHIRD = MOD(MOD(ISCOR1,100),10)
              PNAME    = 'NUCSP'//CHRNOS(IFIRST)
     &                   //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                   //NAMN(IATOM1)//REP(IREPO)//CHRXYZ(ICOOR1)
              IPTYP    = IOP(ICOOR1)
              IC1      = MOD(ICOOR1+1,3) + 1
              IC2      = MOD(ICOOR1  ,3) + 1
              IREP1    = IEOR(IREPO,ISYMAX(IC2,1))
              IREP2    = IEOR(IREPO,ISYMAX(IC1,1))
              IS1      = IPTCNT(3*(IATOM1 - 1) + IC1,IREP1,1)
              IS2      = IPTCNT(3*(IATOM1 - 1) + IC2,IREP2,1)
              IF((IS1.LE.0).OR.(IS2.LE.0))
     &          WRITE(LUPRI,'(A,5I5)')
     &          'DEF_SPINRO: sym. error !',IATOM1,ICOOR1,IC1,IC2
              IFIRST = IS1/100
              ISECND = MOD(IS1,100)/10
              ITHIRD = MOD(MOD(IS1,100),10)
              PLABEL(1)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              IFIRST = IS2/100
              ISECND = MOD(IS2,100)/10
              ITHIRD = MOD(MOD(IS2,100),10)
              PLABEL(2)= 'NEF '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                     //CHRNOS(ITHIRD)//' '
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
              CALL OP1IND('XSPINRA',IPSPINRO(ISCOR1),LLRAPU,NLRAPT,
     &                      INDXPR,MAXLLR)
            ENDIF
          ENDDO

        END IF
        ENDDO
      ENDDO
C
C
C     Define right-hand operators: i.L and i.S
C     ========================================
C
      LINDET=.FALSE.    ! to detect if linear symmetry is used
C
      IF (LINEAR) THEN  ! linear molecule detected using symmetry
C
Ciaa   LINDET is used to record that LINEAR symmetry is used,
C      even when in further steps it could be changed, as occur
C      for UKB, spin-free and ZORA calculations.
       LINDET=.TRUE.    ! save that linear symmetry is used
C
C     Orbital angular momentum
C
       PNAME     = 'X Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'XANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XANGMOM ',IPSPINRO(1+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
C     Spin angular momentum
C
       PNAME     = 'X Spin ang.mom. '
       IPTYP     = 10   ! Sigma_x
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XSPINMOM',IPSPINRO(4+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
      ELSE ! if LINEAR is not TRUE
C
C     Orbital angular momentum
C
       PNAME     = 'X Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'XANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XANGMOM ',IPSPINRO(1+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Y Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'YANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('YANGMOM ',IPSPINRO(2+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Z Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC(1)   = DL   ! factor
       PLABEL(1) = 'ZANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('ZANGMOM ',IPSPINRO(3+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
C     Spin angular momentum
C
       PNAME     = 'X Spin ang.mom. '
       IPTYP     = 10   ! Sigma_x
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XSPINMOM',IPSPINRO(4+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Y Spin ang.mom. '
       IPTYP     = 11   ! Sigma_y
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('YSPINMOM',IPSPINRO(5+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Z Spin ang.mom. '
       IPTYP     = 12   ! Sigma_z
       NPCOMP    = 1    ! one component
       PFAC(1)   = DS   ! factor
       PLABEL(1) = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('ZSPINMOM',IPSPINRO(6+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
      END IF !LINEAR
C
C
C     Print section for spin-rotation:
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Nuclear spin-rotation:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     Print reference
C
C     Expectation values operators
C     ============================
C
      CALL PRSYMB(LUPRI,'-',75,0)
         WRITE(LUPRI,'(1X,A)')
     &   '** Expectation value operators for nuclear spin-rotation **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOOR
        INDOP1 = IPEFN(I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LEXPP(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   '** A operators for nuclear spin-rotation **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NSCOOR
        INDOP1 = IPSPINRO(I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRAPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   '** B operators for nuclear spin-rotation **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,6
        INDOP1 = IPSPINRO(I+MXCOOR)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
      deallocate(DOATOM)
C
      CALL QEXIT('DEF_SPINRO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_SPINRO */
      SUBROUTINE PRP_SPINRO(EXPVAL,ABLRF,ABCNV)
C***********************************************************************
C
C     Output routine for nuclear spin-rotation constants
C
C     Agustin Aucar and Trond Saue - September 2013
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: ABUND, IPRPRP, IPSPINRO
!      use dcbexp_mod, only: IPEFN
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEE, XLR_SKIPEP
      use nuclei_mod, only: IZATOM, NAMN, NATOMS, NFLOAT, NUCDEG,
     &  NUCDEP, NUCIND
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: IPTAX
      use orgcom_mod, only: CMXYZ
!      use spinro_mod, only: LINDET, SPRDNL
!      use dgroup_mod, only: LINEAR
!      use dcbham_mod, only: SPINFR, URKBAL, ZORA
      use consts_mod, only: D0
      use memory_allocator
C
      implicit none
C
      integer I,J
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "spinro.h"
#include "dgroup.h"
#include "dcbham.h"
C
      logical, allocatable :: DOATOM(:)
      integer, allocatable :: NAT(:),ISOT(:)
      double precision, allocatable :: PMATL(:,:,:,:),PMATS(:,:,:,:),
     &                                 TMPMATL(:,:,:),TMPMATS(:,:,:),
     &                                 CSTRA(:),SCTRA(:),GEOM(:,:),
     &                                 MASS(:),AMASS(:),GVAL(:),
     &                                 STAVEC(:)
C
      LOGICAL TEST
      INTEGER IBUFA(3*NUCDEP),IBUFB1(3),IBUFB2(3),IBUFEXP(3*NUCDEP),
     &        IATINF(3,NUCDEP),IA,IB1,IB2,IATOM,IDEG,IDEP,IFREQ,
     &        INUC,ISO,ISYMOP,ITYP,IZ,JA,JATOM,JB1,JB2,K,KA,KB1,KB2,
     &        L,NA,NATOM,NB1,NB2,NCMAT,NEXP,NISO,NSCOOR,NSPEC,IEXP,
     &        JEXP,KEXP,IABLRF_EE,IABLRF_EP
      DOUBLE PRECISION RATINF(3,5,NUCDEP),DISOTP,ABCNV(NBFREQ,NLRBPT),
     &                 EXPVAL(3*NUCDEP),ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &                 ABND,GVALUE
C
      JATOM  = NATOMS + NFLOAT
      NSCOOR = 3*NUCDEP
      NCMAT  = 9*NUCDEP
      allocate(DOATOM(NUCIND))
      allocate(STAVEC(NSCOOR))
      allocate(PMATL(3,3,NUCDEP,2)) ! L part of NSR for ee and ep contributions.
      allocate(PMATS(3,3,NUCDEP,2)) ! S part of NSR for ee and ep contributions.
      allocate(CSTRA(NUCDEP*NUCDEP*9))
      allocate(SCTRA(NUCDEP*NUCDEP*9))
      allocate(TMPMATL(3,3*NUCDEP,2))
      allocate(TMPMATS(3,3*NUCDEP,2))
      allocate(GEOM(NUCDEP,3))
      allocate(MASS(NATOMS))
      allocate(NAT(NATOMS))
      allocate(AMASS(NUCDEP))
      allocate(ISOT(NUCDEP))
      allocate(GVAL(NUCDEP))
C
      CALL QENTER('PRP_SPINRO')
C
      CALL TITLER('NUCLEAR SPIN-ROTATION TENSOR','*',116)
      WRITE(LUPRI,'(13X,A)')
     & 'According to I. A. Aucar et al., JCP 136 (2012) p.204119'
      CALL PRSYMB(LUPRI,'-',81,0)
      WRITE(LUPRI,*) ''
C
C
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-29)
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &  '@ name  charge  isotope   mass         abundance     g factor'
      CALL PRSYMB(LUPRI,'-',61,0)
      NSPEC = 0
      IDEP  = 0
      DO INUC = 1,NUCIND
      IDEG = NUCDEG(INUC)
      IF(DOATOM(INUC)) THEN
        TEST = .FALSE.
        IZ   = IZATOM(INUC)
        NISO = 0
        DO ISO = 1,5
            GVALUE = DISOTP(IZ,ISO,'GVAL')
            ABND = DISOTP(IZ,ISO,'ABUNDANCE')
            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = GVALUE
            ENDIF
        ENDDO
        IF(NISO.GT.0) THEN
          IATINF(1,NSPEC+1) = IDEP + 1
          IATINF(2,NSPEC+1) = IZ
          IATINF(3,NSPEC+1) = NISO
          WRITE(LUPRI,
     &      '(A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &       '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
          DO J = 2,NISO
            WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &         '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
          ENDDO
          DO K = 2,IDEG
            IATINF(1,NSPEC+K) = IDEP + K
            IATINF(2,NSPEC+K) = IZ
            IATINF(3,NSPEC+K) = NISO
            DO L = 1,NISO
              RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
              RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
              RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
            ENDDO
          ENDDO
          NSPEC = NSPEC + IDEG
        ENDIF
      ENDIF
      IDEP = IDEP + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        CALL QEXIT('PRP_SPINRO')
        RETURN
      ENDIF
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
C     *** Extract nuclear spin-rotation constants ***
C
C     The NSCOOR nucl. magn. dipole operators as A operators in <<A;B>> :
      NA = 0
      DO I = 1,NSCOOR
        IF(IPSPINRO(I).GT.0) THEN
          NA = NA + 1
          IBUFA(NA) = I
        ENDIF
      END DO
      IF (NA.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_SPINRO: NA=',NA
        CALL QUIT('PRP_SPINRO: Wrong value of NA!')
      ENDIF
C
Ciaa   As in the current implementation UKB, spin-free and ZORA
C      calculations of energy cannot be performed employing symmetry
C      (they set LINEAR to false in such cases, see subroutine GLINSM)
C      LINDET allows to recover LINEAR symmetry defined in DEF_SPINRO
C      if it is lost due to the use of UKB, SF or ZORA functionalities
      IF ((URKBAL.or.SPINFR.or.ZORA).and.LINDET) LINEAR=.TRUE.
C
C     The 3 L orbital angular momentum directions as B operators in <<A;B>> :
      NB1 = 0
      DO I = 1,3
        IF(IPSPINRO(I+MXCOOR).GT.0) THEN
          NB1 = NB1 + 1
          IBUFB1(NB1) = I
        ENDIF
      END DO
      IF (NB1.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_SPINRO: NB1=',NB1
        CALL QUIT('PRP_SPINRO: Wrong value of NB1!')
      ENDIF
C
C     The 3 S spin angular momentum directions as B operators in <<A;B>> :
      NB2 = 0
      DO I = 1,3
        IF(IPSPINRO(I+3+MXCOOR).GT.0) THEN
          NB2 = NB2 + 1
          IBUFB2(NB2) = I
        ENDIF
      END DO
      IF (NB2.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_SPINRO: NB2=',NB2
        CALL QUIT('PRP_SPINRO: Wrong value of NB2!')
      ENDIF
C
C     The NSCOOR nucl. electric field operators as NEF operators in < 0 | NEF | 0 > :
      NEXP = 0
      DO ITYP = 1,NSCOOR
        IF(IPEFN(ITYP).GT.0 .OR. IPEFN(ITYP).EQ.-1) THEN
          NEXP = NEXP + 1
          IBUFEXP(NEXP) = ITYP
        ENDIF
      END DO
      IF (NEXP.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_SPINRO: NEXP=',NEXP
        CALL QUIT('PRP_SPINRO: Wrong value of NEXP!')
      ENDIF
C
C--------------------------------------------------------------------
C
C     NUCLEAR PART OF NSR TENSOR
C
      CALL CMMASS(GEOM,MASS,NAT,ISOT,IPRPRP)
      CALL NUCSPR(GEOM,CMXYZ,IPRPRP)
C
C-------------------------------------------------------------------
C
C     MEAN VALUE PART OF NSR TENSOR
C
      CALL DZERO(STAVEC,NSCOOR)
C
      DO IEXP = 1,NEXP
        JEXP = IBUFEXP(IEXP)
        KEXP = IPEFN(JEXP)
        IF(KEXP.NE.-1) STAVEC(IEXP) = EXPVAL(KEXP)
      END DO
C
      CALL DIASPR(STAVEC,EXPVAL,CSTRA,SCTRA,NCMAT,IPRPRP)
C
C--------------------------------------------------------------------
C
C     LINEAR RESPONSE PART OF NSR TENSOR
C
C     Fill nuclear spin-rotation matrices
C
      IF (XLR_SKIPEE) THEN
        IABLRF_EE = -999999999
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 2
        END IF
      ELSE
        IABLRF_EE = 2
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 3
        END IF
      END IF
C
C
      DO IFREQ = 1,NBFREQ
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
           CALL PRSYMB(LUPRI,'-',72,2)
           WRITE(LUPRI,'(A,F15.8,A)')
     &          '*** Frequency    :',BFREQ(IFREQ),' a.u.'
        ENDIF
        CALL DZERO(TMPMATL,2*NCMAT)
        CALL DZERO(TMPMATS,2*NCMAT)
C
        DO IB1 = 1,NB1
          JB1 = IBUFB1(IB1)
          KB1 = IPSPINRO(MXCOOR+JB1)
          JB1 = IPTAX(JB1,2)
          DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPSPINRO(JA)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         TMPMATL(JB1,JA,1) = ABLRF(IFREQ,KA,KB1,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         TMPMATL(JB1,JA,2) = ABLRF(IFREQ,KA,KB1,IABLRF_EP)
          ENDDO
        ENDDO
C
        DO IB2 = 1,NB2
          JB2 = IBUFB2(IB2)
          KB2 = IPSPINRO(MXCOOR+3+JB2)
          JB2 = IPTAX(JB2,2)
          DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPSPINRO(JA)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         TMPMATS(JB2,JA,1) = ABLRF(IFREQ,KA,KB2,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         TMPMATS(JB2,JA,2) = ABLRF(IFREQ,KA,KB2,IABLRF_EP)
          ENDDO
        ENDDO
C
C     Transform to C1 basis
C
        CALL DZERO(PMATL,2*NCMAT)
        CALL DZERO(PMATS,2*NCMAT)
C
        IF (.NOT. XLR_SKIPEE) THEN
C           .... extract the e-e response part into PMAT(...,1)
            CALL TRADIP(TMPMATL,PMATL,CSTRA,SCTRA,3*NUCDEP,2,2)
            CALL TRADIP(TMPMATS,PMATS,CSTRA,SCTRA,3*NUCDEP,2,2)
        END IF
        IF (.NOT.XLR_SKIPEP) THEN
C           ... extract the e-p response part into PMAT(...,2)
            CALL TRADIP(TMPMATL(1,1,2),PMATL(1,1,1,2),
     &                CSTRA,SCTRA,3*NUCDEP,2,2)
            CALL TRADIP(TMPMATS(1,1,2),PMATS(1,1,1,2),
     &                CSTRA,SCTRA,3*NUCDEP,2,2)
        END IF
      END DO
C
C--------------------------------------------------------------------
      CALL SPRRES(GEOM,AMASS,GVAL,SPRDNL,PMATL,PMATS,JATOM,
     &            IPRPRP,DOATOM)
C--------------------------------------------------------------------
      CALL QEXIT('PRP_SPINRO')
C
      deallocate(AMASS)
      deallocate(ISOT)
      deallocate(GVAL)
      deallocate(GEOM)
      deallocate(MASS)
      deallocate(NAT)
      deallocate(DOATOM)
      deallocate(PMATL)
      deallocate(PMATS)
      deallocate(CSTRA)
      deallocate(SCTRA)
      deallocate(STAVEC)
      deallocate(TMPMATL)
      deallocate(TMPMATS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ANGLES(THETA,PHI,CHI,WAVVEC,POLVEC)
C***********************************************************************
C     Define orientation of wave and polarization vector,
C     based on angles read in degrees.
C     For this we use a right-handed coordinate system defined by
C        e_r     = ( sin[Theta}cos[Phi], sin[Theta]sin[Phi], cos[Theta])
C        e_Theta = ( cos[Theta]cos[Phi], cos[Theta]sin[Phi],-sin[Theta])
C        e_Phi   = (-sin[Phi]          , cos[Phi]          ,0          )
C     When Theta and Phi are set to zero we have
C        (e_r,e_Theta,e_phi) --> (e_z, e_x, e_y)
C     We align the wave vector with e_r. The polarization vector epsilon
C     is then in the perpendicular plane spanned by e_Theta and e_Phi,
C     so we write
C        epsilon = cos[Chi]e_Theta + sin[Chi]e_Phi
C     where we have introduced a third angle Theta
C
C     Written by Trond Saue, July 18 2017
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcbxpp_mod, only: BEDECD
#include "implicit.h"
      INTEGER :: I
#include "dcbxpp.h"
      DIMENSION WAVVEC(3),POLVEC(3,2)
C
      FAC      = ACOS(-1.0D0)/180.0D0
      THETA2   = THETA*FAC
      PHI2     = PHI*FAC
      CHI2     = CHI*FAC
      WAVVEC(1) = SIN(THETA2)*COS(PHI2)
      WAVVEC(2) = SIN(THETA2)*SIN(PHI2)
      WAVVEC(3) = COS(THETA2)
      W1       = COS(CHI2)
      W2       = SIN(CHI2)
      POLVEC(1,1)  = W1*COS(THETA2)*COS(PHI2) - W2*SIN(PHI2)
      POLVEC(2,1)  = W1*COS(THETA2)*SIN(PHI2) + W2*COS(PHI2)
      POLVEC(3,1)  =-W1*SIN(THETA2)
      IF (BEDECD) THEN
        POLVEC(1,2)  = -W2*COS(THETA2)*COS(PHI2) - W1*SIN(PHI2)
        POLVEC(2,2)  = -W2*COS(THETA2)*SIN(PHI2) + W1*COS(PHI2)
        POLVEC(3,2)  = W2*SIN(THETA2)
      ENDIF
C.....Output section
      IF (BEDECD) THEN
        WRITE(LUPRI,'(1X,A,3F18.10)')
     &   '* Unit wave vector   :', (WAVVEC(I),I=1,3),
     &   '* Polarization Vector 1:', (POLVEC(I,1),I=1,3),
     &   '* Polarization vector 2:', (POLVEC(I,2),I=1,3)
      ELSE

      WRITE (LUPRI,'(1X,A,3F18.10)')
     &     '* Unit wave vector   :', (WAVVEC(I),I=1,3),
     &     '* Polarization vector:', (POLVEC(I,1),I=1,3)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck  nucspr */
      SUBROUTINE NUCSPR(GEOM,ORIGIN,IPRINT)
C***********************************************************************
C
C     Calculate nuclear contribution to nuclear spin-rotation constants
C     K.Ruud, October-94
C
C***********************************************************************
C
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use nuclei_mod, only: CHARGE, ISTBNU, NATOMS, NFLOAT, NUCDEP,
     &  NUCIND
      use symmet_mod, only: MAXOPR
!      use spinro_mod, only: SPRNUC, SPRNUC2, SPRNUC3
      use consts_mod, only: D0
      implicit none
C
      integer I,J
C
#include "spinro.h"
C
      double precision GEOM,ORIGIN
C
      DIMENSION GEOM(3*(NATOMS+NFLOAT)),ORIGIN(3)
C
      integer IPRINT,IATOM1,IATOM2,ISYM1,ISYM2,
     &        KK,LL,NATOM1
C
      double precision XDIF,YDIF,ZDIF,XDIFCM,YDIFCM,ZDIFCM,
     &                 XDIFCM2,YDIFCM2,ZDIFCM2,ZVAL2,PFAC
C
      IF (IPRINT.GE.5) THEN
       CALL HEADER(
     & 'NUCSPR: Entering parameters to calculate nuclear'//
     & ' contribution to nuclear spin-rotation',-1)
       write(lupri,*) 'ORIGIN=',ORIGIN(1),ORIGIN(2),ORIGIN(3)
       do i=1, (3*(NATOMS+NFLOAT))
        write(lupri,*) 'Entering geom(',i,')=',GEOM(i)
       enddo
      ENDIF
C
      CALL DZERO(SPRNUC,9*NUCDEP)
      CALL DZERO(SPRNUC2,9*NUCDEP)
      CALL DZERO(SPRNUC3,9*NUCDEP)
C
      KK = 1
      NATOM1 = 1
      DO 10 IATOM1 = 1, NUCIND
         DO 20 ISYM1 = 0, MAXOPR
            IF (IAND(ISYM1,ISTBNU(IATOM1)) .EQ. 0) THEN
               IF (CHARGE(IATOM1) .GT. D0) THEN
                  LL = 1
                  DO 30 IATOM2 = 1, NUCIND
                     DO 40 ISYM2 = 0, MAXOPR
                        IF (IAND(ISYM2,ISTBNU(IATOM2)) .EQ. 0) THEN
                           ZVAL2 = CHARGE(IATOM2)
                           IF (ZVAL2 .GT. D0) THEN
                              IF (KK .NE. LL) THEN
C
                                 XDIF = GEOM(KK) - GEOM(LL)
                                 YDIF = GEOM(KK + 1) - GEOM(LL + 1)
                                 ZDIF = GEOM(KK + 2) - GEOM(LL + 2)
C
                                 XDIFCM = GEOM(KK) - ORIGIN(1)
                                 YDIFCM = GEOM(KK + 1) - ORIGIN(2)
                                 ZDIFCM = GEOM(KK + 2) - ORIGIN(3)
C
                                 XDIFCM2 = GEOM(LL) - ORIGIN(1)
                                 YDIFCM2 = GEOM(LL + 1) - ORIGIN(2)
                                 ZDIFCM2 = GEOM(LL + 2) - ORIGIN(3)
C
                                 PFAC  = ZVAL2/(SQRT(XDIF*XDIF
     &                                + YDIF*YDIF + ZDIF*ZDIF))**3
C
                                SPRNUC(1,1,NATOM1) = SPRNUC(1,1,NATOM1)
     &                                + (YDIF*YDIF + ZDIF*ZDIF)*PFAC
                                SPRNUC(2,2,NATOM1) = SPRNUC(2,2,NATOM1)
     &                                + (XDIF*XDIF + ZDIF*ZDIF)*PFAC
                                SPRNUC(3,3,NATOM1) = SPRNUC(3,3,NATOM1)
     &                                + (XDIF*XDIF + YDIF*YDIF)*PFAC
                                SPRNUC(2,1,NATOM1) = SPRNUC(2,1,NATOM1)
     &                                - YDIF*XDIF*PFAC
                                SPRNUC(3,1,NATOM1) = SPRNUC(3,1,NATOM1)
     &                                - ZDIF*XDIF*PFAC
                                SPRNUC(1,2,NATOM1) = SPRNUC(1,2,NATOM1)
     &                                - XDIF*YDIF*PFAC
                                SPRNUC(3,2,NATOM1) = SPRNUC(3,2,NATOM1)
     &                                - ZDIF*YDIF*PFAC
                                SPRNUC(1,3,NATOM1) = SPRNUC(1,3,NATOM1)
     &                                - XDIF*ZDIF*PFAC
                                SPRNUC(2,3,NATOM1) = SPRNUC(2,3,NATOM1)
     &                                - YDIF*ZDIF*PFAC
C
C     Thomas precesion contribution to nuclear spin-rotation tensor
C
                              SPRNUC2(1,1,NATOM1) = SPRNUC2(1,1,NATOM1)
     &                               - (YDIFCM*YDIF + ZDIFCM*ZDIF)*PFAC
                              SPRNUC2(2,2,NATOM1) = SPRNUC2(2,2,NATOM1)
     &                               - (XDIFCM*XDIF + ZDIFCM*ZDIF)*PFAC
                              SPRNUC2(3,3,NATOM1) = SPRNUC2(3,3,NATOM1)
     &                               - (XDIFCM*XDIF + YDIFCM*YDIF)*PFAC
                              SPRNUC2(2,1,NATOM1) = SPRNUC2(2,1,NATOM1)
     &                               + YDIFCM*XDIF*PFAC
                              SPRNUC2(3,1,NATOM1) = SPRNUC2(3,1,NATOM1)
     &                               + ZDIFCM*XDIF*PFAC
                              SPRNUC2(1,2,NATOM1) = SPRNUC2(1,2,NATOM1)
     &                               + XDIFCM*YDIF*PFAC
                              SPRNUC2(3,2,NATOM1) = SPRNUC2(3,2,NATOM1)
     &                               + ZDIFCM*YDIF*PFAC
                              SPRNUC2(1,3,NATOM1) = SPRNUC2(1,3,NATOM1)
     &                               + XDIFCM*ZDIF*PFAC
                              SPRNUC2(2,3,NATOM1) = SPRNUC2(2,3,NATOM1)
     &                               + YDIFCM*ZDIF*PFAC
C
C     Nuclear contribution to NSR for molecules in equilibrium
C
                            SPRNUC3(1,1,NATOM1) = SPRNUC3(1,1,NATOM1)
     &                             - (YDIFCM2*YDIF + ZDIFCM2*ZDIF)*PFAC
                            SPRNUC3(2,2,NATOM1) = SPRNUC3(2,2,NATOM1)
     &                             - (XDIFCM2*XDIF + ZDIFCM2*ZDIF)*PFAC
                            SPRNUC3(3,3,NATOM1) = SPRNUC3(3,3,NATOM1)
     &                             - (XDIFCM2*XDIF + YDIFCM2*YDIF)*PFAC
                            SPRNUC3(2,1,NATOM1) = SPRNUC3(2,1,NATOM1)
     &                             + YDIFCM2*XDIF*PFAC
                            SPRNUC3(3,1,NATOM1) = SPRNUC3(3,1,NATOM1)
     &                             + ZDIFCM2*XDIF*PFAC
                            SPRNUC3(1,2,NATOM1) = SPRNUC3(1,2,NATOM1)
     &                             + XDIFCM2*YDIF*PFAC
                            SPRNUC3(3,2,NATOM1) = SPRNUC3(3,2,NATOM1)
     &                             + ZDIFCM2*YDIF*PFAC
                            SPRNUC3(1,3,NATOM1) = SPRNUC3(1,3,NATOM1)
     &                             + XDIFCM2*ZDIF*PFAC
                            SPRNUC3(2,3,NATOM1) = SPRNUC3(2,3,NATOM1)
     &                             + YDIFCM2*ZDIF*PFAC
                              END IF
                              LL = LL + 3
                           END IF
                        END IF
 40                  CONTINUE
 30               CONTINUE
                  KK = KK + 3
               END IF
               NATOM1 = NATOM1 + 1
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diaspr */
      SUBROUTINE DIASPR(EXPVAL,ELFLD,CSTRA,SCTRA,NCMAT,IPRINT)
C***********************************************************************
C
C     Routine for calculating the expectation value contribution to the
C     nuclear spin-rotation constant, K.Ruud - Oct.-95
C
C***********************************************************************
C
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CORD, ISTBNU, NUCDEP, NUCIND
      use symmet_mod, only: ISYMAX, MAXOPR, PT
!      use spinro_mod, only: SPRDNL
      use orgcom_mod, only: CMXYZ
      implicit none
C
      integer I,J
C
#include "spinro.h"
C
      integer IPRINT,IATOM,ISYMOP,NATOM,NCMAT
C
      double precision EXPVAL(3*NUCDEP), ELFLD(3,NUCDEP),
     &                 CSTRA(3*NUCDEP,3*NUCDEP),
     &                 SCTRA(3*NUCDEP,3*NUCDEP),
     &                 COORX,COORY,COORZ
C
      CALL DZERO(SPRDNL,NCMAT)
C
C     We start by sorting the electric field results
C
      CALL TRACOR(CSTRA,SCTRA,1,3*NUCDEP,0)
      CALL DGEMM('N','T',1,3*NUCDEP,3*NUCDEP,1.D0,
     &           EXPVAL,1,
     &           SCTRA(1,1),3*NUCDEP,0.D0,
     &           ELFLD,1)
C
C     Then we construct the expectation value contribution
C
      NATOM = 1
      DO 10 IATOM = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
         IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
            COORX = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
     &            - CMXYZ(1)
            COORY = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
     &            - CMXYZ(2)
            COORZ = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
     &            - CMXYZ(3)
            SPRDNL(1,1,NATOM) = (ELFLD(2,NATOM)*COORY
     &                        +  ELFLD(3,NATOM)*COORZ)
            SPRDNL(2,2,NATOM) = (ELFLD(1,NATOM)*COORX
     &                        +  ELFLD(3,NATOM)*COORZ)
            SPRDNL(3,3,NATOM) = (ELFLD(1,NATOM)*COORX
     &                        +  ELFLD(2,NATOM)*COORY)
            SPRDNL(1,2,NATOM) = -ELFLD(2,NATOM)*COORX
            SPRDNL(1,3,NATOM) = -ELFLD(3,NATOM)*COORX
            SPRDNL(2,1,NATOM) = -ELFLD(1,NATOM)*COORY
            SPRDNL(2,3,NATOM) = -ELFLD(3,NATOM)*COORY
            SPRDNL(3,1,NATOM) = -ELFLD(1,NATOM)*COORZ
            SPRDNL(3,2,NATOM) = -ELFLD(2,NATOM)*COORZ
            NATOM = NATOM + 1
         END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sprres */
      SUBROUTINE SPRRES(GEOM,AMASS,GVAL,DMAT,PMATL,PMATS,KATOM,
     &                  IPRINT,DOATOM)
C***********************************************************************
C
C     Analyze routine for nuclear spin-rotation constants
C     This routine is modelled after SPRRES (Prof K. Ruud) in DALTON
C
C   INPUT:
C     DMAT  - expectation value term of NSR
C     PMATL - contribution from orbital angular moment to linear response term
C     PMATS - contribution from spin angular moment to linear response term
C
C     Agustin Aucar and Trond Saue - June  2014
C     Last revision: Agustin Aucar - March 2024
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CHARGE, CORD, ISOTOP, ISTBNU, IZATOM,
     &  NAMDEP, NUCDEP, NUCIND
      use symmet_mod, only: IPTAX, ISYMAX, MAXOPR, PT
!      use spinro_mod, only: GTRANT, SPRNUC, SPRNUC2, SPRNUC3
      use orgcom_mod, only: CMXYZ
!      use dgroup_mod, only: LINEAR
!      use codata_mod, only: ALPHA2, PMASS, XFAMU, XTHZ
      use consts_mod, only: D0, D1, D2, D6
      implicit none
C
      integer I,J
C
#include "spinro.h"
#include "dgroup.h"
#include "codata.h"
C
      integer K,KATOM,IPRINT,IATOM,ICOOR,ISYMOP,IX,IY,
     &        JATOM,JX,JY,NATTYP,NUMIS,KHZPPM
C
      double precision FAC,FACPPM,FACTOT,FACTHOMAS,CONST,
     &                 GEOM(KATOM,3),AMASS(KATOM),GVAL(KATOM),
     &                 EIGVAL(3),EIGVEC(3,3),OMEGA(3),TINERT(9),
     &                 DISOTP,ANGUMOM(3),DMAT(3,3,KATOM),
     &                 PMATL(3,3,KATOM,2),PMATS(3,3,KATOM,2),
     &                 NSRTOT,NSRNUC,NSRNUCTHO,NSREXV,NSREXVTHO,NSRLR,
     &                 NSRLRLEE,NSRLRLEP,NSRLRSEE,NSRLRSEP,NSRNUCEQ,
     &                 GTRANN(3,3,KATOM),GTRANNTHO(3,3,KATOM),
     &                 GTRANNUC(3,3,KATOM),GTRAND(3,3,KATOM),
     &                 GTRANEXP(3,3,KATOM),GTRANDTHO(3,3,KATOM),
     &                 GTRANPLEE(3,3,KATOM),GTRANPLEP(3,3,KATOM),
     &                 GTRANPSEE(3,3,KATOM),GTRANPSEP(3,3,KATOM),
     &                 GTRANP(3,3,KATOM),GTRANNEQ(3,3,KATOM),
     &                 GTRANTNEQ(3,3,KATOM)
C
      LOGICAL PLANAR,LINMOL,DOATOM(NUCIND)
C
      FACPPM = ALPHA2/D2*1.0D6
      FAC = ALPHA2/(D2*PMASS*XFAMU*XFAMU)*XTHZ*1.0D-3
C
C
Ciaa  Nuclear rotational angular momentum could be defined here to obtain Omega afterwards
      DO 300 I = 1,3
        ANGUMOM(I) = D1
 300  CONTINUE
C
C
C     Determine principal moments of inertia
C     ======================================
C
C     Set up coordinates (relative to center of mass), g-values and masses of all nuclei
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         NATTYP = IZATOM(IATOM)
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
         IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
           JATOM = JATOM + 1
           AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
           GVAL(JATOM)  = DISOTP(NATTYP,NUMIS,'GVAL')
           DO 120 ICOOR = 1, 3
             GEOM(JATOM,ICOOR) = PT(IAND(ISYMAX(ICOOR,1),ISYMOP))
     &                          *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
 120       CONTINUE
         END IF
 110     CONTINUE
 100  CONTINUE
C
      CALL WLKDIN(GEOM,AMASS,KATOM,ANGUMOM,TINERT,OMEGA,EIGVAL,EIGVEC,
     &            .TRUE.,PLANAR,LINMOL,IPRINT)
C     On output
C        TINERT - inertia tensor
C        EIGVAL - inverse eigenvalues (for linear molecules, in descending order)
C        EIGVEC - corresponding eigenvectors
C        PLANAR=T : planar molecule
C        LINMOL=T : linear molecule (EIGVAL(1) = 0.0)
C
C     Transform each nuclear spin-rotation constant to principal axis system,
C     and multiply by inverse moment of inertia and various constants
C
Ciaa  To obtain Omega, the nuclear rot. angular momentum (ANGUMOM) should be defined before
      IF (IPRINT .GE. 8) THEN
       CALL HEADER
     &    ('Nuclear angular velocity (a.u.)',0)
       WRITE (LUPRI,'(3X,A,2X,F16.8)')
     &       'Omega_A',OMEGA(1),
     &       'Omega_B',OMEGA(2),
     &       'Omega_C',OMEGA(3)
      END IF
C
      WRITE (LUPRI,'(/,2X,A,/,11X,A,//,11X,A,/,11X,A,/,11X,A)')
     &      'WARNING: Nuclear contributions to the NSR'           //
     &              ' tensors/constants are','only valid for'     //
     &              ' molecules in their equilibrium geometry.',
     &              'To get valid results out of the equilibrium,'//
     &              ' please use','the .PRINT keyword of the'     //
     &              ' **PROPERTIES module with',
     &              'at least a value of 4.'
C
C    Print inertia tensor or moment
C
      IF (IPRINT .GE. 0) THEN
       IF (LINEAR) THEN      ! Linear molecule detected using symmetry
              WRITE (LUPRI,'(/,2X,A,F16.8,/)')
     &         'Moment of inertia (a.u.): ', XFAMU/EIGVAL(3)
       ELSE
        IF(LINMOL) THEN      ! Linear molecule (without symmetry detection)
         WRITE (LUPRI,'(/,2X,A,F16.8,/)')
     &         'Moment of inertia (a.u.): ', XFAMU/EIGVAL(3)
        ELSE                 ! Non-linear molecule
         CALL HEADER
     &    ('Principal moments of inertia (a.u.) and principal axes',0)
         WRITE (LUPRI,'(3X,A,2X,F16.8,6X,3F16.8)')
     &       'IA',XFAMU/EIGVAL(1), (EIGVEC(1,I),I = 1,3),
     &       'IB',XFAMU/EIGVAL(2), (EIGVEC(2,I),I = 1,3),
     &       'IC',XFAMU/EIGVAL(3), (EIGVEC(3,I),I = 1,3)
         WRITE (LUPRI,'(/,2X,A)') 'Nuclear spin-rotation tensor'      //
     &               ' elements are given in the principal axis system.'
        END IF
       END IF
      END IF
C
C    Print nuclear spin-rotation constants or tensors
C
      IATOM = 0
      DO 98 K = 1, NUCIND
       DO 200 ISYMOP = 0, MAXOPR
        IF (IAND(ISTBNU(K),ISYMOP).EQ.0) THEN
         IATOM = IATOM + 1
         IF (DOATOM(K)) THEN
          IF (CHARGE(IATOM).GT. D0) THEN
C
C    Print g-values of selected nuclei
C
          IF (GVAL(IATOM) .NE. D0) THEN
           FACTHOMAS=(CHARGE(IATOM)*PMASS)/(AMASS(IATOM)*GVAL(IATOM))
           IF (IPRINT .GE. 0) THEN
            CALL HEADER('Nuclear spin-rotation constant (in kHz) for '//
     &                   NAMDEP(IATOM),-1)
            WRITE (LUPRI,'(2X,A17,F12.6,/)')
     &                 'Nuclear g-value: ',GVAL(IATOM)
           END IF
          ELSE   !   if GVAL = 0
           FACTHOMAS=(CHARGE(IATOM)*PMASS)/(AMASS(IATOM))
           IF (IPRINT .GE. 0) THEN
          CALL HEADER('Artificial nuclear spin-rotation constant'//
     &                ' (in kHz) for '//
     &                 NAMDEP(IATOM),-1)
            WRITE (LUPRI,'(2X,A,F8.3)')
     &  'WARNING: Nuclear g-value NOT available for isotope of mass ',
     &  AMASS(IATOM)
            WRITE (LUPRI,'(2X,A,/)')
     &  'Nuclear g-value is artificially set equal to 1'
           END IF
          END IF
C
C
C         Linear molecule detected using symmetry
C
          IF (LINEAR) THEN
C
           IF (GVAL(IATOM) .NE. D0) THEN
            FACTOT = EIGVAL(3)*GVAL(IATOM)*FAC
           ELSE
            FACTOT = EIGVAL(3)*FAC  ! GVAL(IATOM) = 1
           END IF
C
           DO 46 KHZPPM= 1,2          ! Results in 1: kHz ; 2: ppm
            IF (KHZPPM.EQ.2) THEN
            FACTOT = FACPPM
             IF (IPRINT .GE. 2) THEN
            CALL HEADER('Nuclear spin-rotation constant (in ppm) for '//
     &                    NAMDEP(IATOM),-1)
             END IF
            END IF
C           Nuclear contribution to NSR for molecules in equilibrium
            NSRNUCEQ = SPRNUC3(1,1,IATOM)*FACTOT
C           Linear response contributions to NSR
            NSRLRLEE = PMATL(1,1,IATOM,1)*FACTOT
            NSRLRLEP = PMATL(1,1,IATOM,2)*FACTOT
            NSRLRSEE = PMATS(1,1,IATOM,1)*FACTOT
            NSRLRSEP = PMATS(1,1,IATOM,2)*FACTOT
C           Total linear response contribution to NSR
            NSRLR=NSRLRLEE+NSRLRLEP+NSRLRSEE+NSRLRSEP
C
C           Total NSR for molecules in equilibrium
C
            NSRTOT=NSRNUCEQ+NSRLR
C
C           For molecules in non-equilibrium (general case)
C
C           Nuclear contributions to NSR
            NSRNUC   = SPRNUC(1,1,IATOM)*FACTOT
            NSRNUCTHO= SPRNUC2(1,1,IATOM)*FACTOT*FACTHOMAS
            NSRNUC   = NSRNUC + NSRNUCTHO
C           Expectation value contributions to NSR
            NSREXV   = DMAT(1,1,IATOM)*FACTOT
            NSREXVTHO= DMAT(1,1,IATOM)*FACTOT*FACTHOMAS
            NSREXV   = NSREXV - NSREXVTHO
C
C
             IF (IPRINT .GE. 2*KHZPPM-2)
     &  WRITE(LUPRI,'((1X,A,F16.8,/),/,1X,A,F16.8,/,1X,A,F16.8)')
     & ' Total spin-rotation constant (SRC)      :',NSRTOT,
     & '   Nuclear contribution to SRC (M^nuc)   :',NSRNUCEQ,
     & '   Electronic contribution to SRC (M^LR) :',NSRLR
C
            IF (IPRINT .GE. 2*KHZPPM-1)
     &  WRITE(LUPRI,'(3(1X,A,F16.8,/),(1X,A,F16.8))')
     & '           M^LR-L(e-e) :',NSRLRLEE,
     & '           M^LR-S(e-e) :',NSRLRSEE,
     & '           M^LR-L(e-p) :',NSRLRLEP,
     & '           M^LR-S(e-p) :',NSRLRSEP
C
            IF (IPRINT .GE. 4)
     &  WRITE(LUPRI,'(/,A,3(/,1X,A,F16.8),/,1X,A,F16.8,/,A)')
     & '******************************************************',
     & '********      M^NU         ********:',NSRNUC,
     & '********      M^EV         ********:',NSREXV,
     & '********      M^LR         ********:',NSRLR,
     & '********      M^total      ********:',NSRNUC+NSREXV+NSRLR,
     & '******************************************************'
C
 46        CONTINUE
C
C         Linear symmetry is not detected
C
          ELSE
C
           IF (GVAL(IATOM) .NE. D0) THEN
            FACTOT = GVAL(IATOM)*FAC
           ELSE
            FACTOT = FAC  ! GVAL(IATOM) = 1
           END IF
C
           DO 47 KHZPPM= 1,2          ! Results in 1: kHz ; 2: ppm
            IF (KHZPPM.EQ.2.AND.IPRINT.GE.2) THEN
            CALL HEADER('Nuclear spin-rotation constant (in ppm) for '//
     &                    NAMDEP(IATOM),-1)
            END IF
C
            CALL DZERO(GTRANT,9*NUCDEP)
            CALL DZERO(GTRANNEQ,9*NUCDEP)
            CALL DZERO(GTRANN,9*NUCDEP)
            CALL DZERO(GTRANNTHO,9*NUCDEP)
            CALL DZERO(GTRANNUC,9*NUCDEP)
            CALL DZERO(GTRAND,9*NUCDEP)
            CALL DZERO(GTRANDTHO,9*NUCDEP)
            CALL DZERO(GTRANEXP,9*NUCDEP)
            CALL DZERO(GTRANP,9*NUCDEP)
            CALL DZERO(GTRANPLEE,9*NUCDEP)
            CALL DZERO(GTRANPLEP,9*NUCDEP)
            CALL DZERO(GTRANPSEE,9*NUCDEP)
            CALL DZERO(GTRANPSEP,9*NUCDEP)
C
            DO 50 I = 1, 3
             DO 50 J = 1, 3
              DO 60 IX = 1, 3
               DO 60 IY = 1, 3
                JX = IPTAX(IX,2)
                JY = IPTAX(IY,2)
                IF (KHZPPM.EQ.1) THEN
                 CONST = EIGVAL(I)*FACTOT
                ELSE
                 CONST = FACPPM
                END IF
C           Nuclear contribution to NSR for molecules in equilibrium
                GTRANNEQ(I,J,IATOM)=GTRANNEQ(I,J,IATOM)+EIGVEC(IX,I)
     &                       *SPRNUC3(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
C           Linear response contributions to NSR
                GTRANPLEE(I,J,IATOM) = GTRANPLEE(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATL(JX,JY,IATOM,1)*EIGVEC(IY,J)*CONST
                GTRANPLEP(I,J,IATOM) = GTRANPLEP(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATL(JX,JY,IATOM,2)*EIGVEC(IY,J)*CONST
                GTRANPSEE(I,J,IATOM) = GTRANPSEE(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATS(JX,JY,IATOM,1)*EIGVEC(IY,J)*CONST
                GTRANPSEP(I,J,IATOM) = GTRANPSEP(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATS(JX,JY,IATOM,2)*EIGVEC(IY,J)*CONST
C           Total linear response contribution to NSR
                GTRANP(I,J,IATOM) = GTRANPLEE(I,J,IATOM)
     &                       + GTRANPLEP(I,J,IATOM)
     &                       + GTRANPSEE(I,J,IATOM)
     &                       + GTRANPSEP(I,J,IATOM)
C
C           Total NSR for molecules in equilibrium
C
            GTRANT(I,J,IATOM) = GTRANNEQ(I,J,IATOM) + GTRANP(I,J,IATOM)
C
C           For molecules in non-equilibrium (general case)
C
C           Nuclear contributions to NSR
C
                GTRANN(I,J,IATOM)=GTRANN(I,J,IATOM)+EIGVEC(IX,I)
     &                       *SPRNUC(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
                GTRANNTHO(I,J,IATOM)=GTRANNTHO(I,J,IATOM)+EIGVEC(IX,I)
     &                       *SPRNUC2(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
     &                       *FACTHOMAS
C           Total nuclear contribution to NSR
                GTRANNUC(I,J,IATOM) = GTRANN(I,J,IATOM)
     &                               +GTRANNTHO(I,J,IATOM)
C           Expectation value contributions to NSR
                GTRAND(I,J,IATOM)=GTRAND(I,J,IATOM)+EIGVEC(IX,I)
     &                       *DMAT(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
                GTRANDTHO(I,J,IATOM)=GTRANDTHO(I,J,IATOM)+EIGVEC(IX,I)
     &                       *DMAT(JX,JY,IATOM)*EIGVEC(IY,J)*CONST
     &                       *(-FACTHOMAS)
C           Total expectation value contribution to NSR
                GTRANEXP(I,J,IATOM) = GTRAND(I,J,IATOM)
     &                               +GTRANDTHO(I,J,IATOM)
C           TOTAL NSR for general case
                GTRANTNEQ(I,J,IATOM) = GTRANNUC(I,J,IATOM)
     &                       + GTRANEXP(I,J,IATOM) + GTRANP(I,J,IATOM)
 60           CONTINUE
 50         CONTINUE
C
C    At least one inertia tensor eigenvalue equal to zero (linear molecule):
            IF(LINMOL) THEN
C
             IF (IPRINT .GE. 2*KHZPPM-2)
     &  WRITE(LUPRI,'((1X,A,F16.8,/),/,1X,A,F16.8,/,1X,A,F16.8)')
     & ' Total spin-rotation constant (SRC)      :',GTRANT(3,3,IATOM),
     & '   Nuclear contribution to SRC (M^nuc)   :',GTRANNEQ(3,3,IATOM),
     & '   Electronic contribution to SRC (M^LR) :',GTRANP(3,3,IATOM)
C
            IF (IPRINT .GE. 2*KHZPPM-1)
     &  WRITE(LUPRI,'(3(1X,A,F16.8,/),(1X,A,F16.8))')
     & '           M^LR-L(e-e) :',GTRANPLEE(3,3,IATOM),
     & '           M^LR-S(e-e) :',GTRANPSEE(3,3,IATOM),
     & '           M^LR-L(e-p) :',GTRANPLEP(3,3,IATOM),
     & '           M^LR-S(e-p) :',GTRANPSEP(3,3,IATOM)
C
            IF (IPRINT .GE. 4)
     &  WRITE(LUPRI,'(/,A,4(/,1X,A,F16.8),/,A)')
     & '******************************************************',
     & '********      M^NU         ********:',GTRANNUC(3,3,IATOM),
     & '********      M^EV         ********:',GTRANEXP(3,3,IATOM),
     & '********      M^LR         ********:',GTRANP(3,3,IATOM),
     & '********      M^total      ********:',GTRANTNEQ(3,3,IATOM),
     & '******************************************************'
C
C
C    Inertia tensor of a non-linear molecule (with and without symmetry):
            ELSE
             IF (IPRINT .GE. 2*KHZPPM-2) THEN
C
C           Total NSR for molecules in equilibrium
C
             CALL HEADER('Total nuclear spin-rotation tensor (NSRT)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANT(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANT(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANT(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANT(1,1,IATOM)
     &              +GTRANT(2,2,IATOM)+GTRANT(3,3,IATOM))/3
C           Nuclear contribution to NSR for molecules in equilibrium
              CALL HEADER('Nuclear contribution to NSRT (M^nuc)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANNEQ(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANNEQ(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANNEQ(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANNEQ(1,1,IATOM)
     &              +GTRANNEQ(2,2,IATOM)+GTRANNEQ(3,3,IATOM))/3
C           Total linear response contribution to NSR
              CALL HEADER('Electronic contribution to NSRT (M^LR)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANP(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANP(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANP(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANP(1,1,IATOM)
     &              +GTRANP(2,2,IATOM)+GTRANP(3,3,IATOM))/3
             END IF
C
C           Linear response contributions to NSR
C
             IF (IPRINT .GE. 2*KHZPPM-1) THEN
              CALL HEADER('M^LR-L(e-e)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANPLEE(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPLEE(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPLEE(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANPLEE(1,1,IATOM)
     &              +GTRANPLEE(2,2,IATOM)+GTRANPLEE(3,3,IATOM))/3
              CALL HEADER('M^LR-S(e-e)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANPSEE(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPSEE(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPSEE(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANPSEE(1,1,IATOM)
     &              +GTRANPSEE(2,2,IATOM)+GTRANPSEE(3,3,IATOM))/3
              CALL HEADER('M^LR-L(e-p)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANPLEP(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPLEP(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPLEP(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANPLEP(1,1,IATOM)
     &              +GTRANPLEP(2,2,IATOM)+GTRANPLEP(3,3,IATOM))/3
              CALL HEADER('M^LR-S(e-p)',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANPSEP(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPSEP(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPSEP(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANPSEP(1,1,IATOM)
     &              +GTRANPSEP(2,2,IATOM)+GTRANPSEP(3,3,IATOM))/3
             END IF
C
C           For molecules in non-equilibrium (general case)
C
C           Total nuclear contribution to NSR
             IF (IPRINT .GE. 4) THEN
              WRITE(LUPRI,'(/,/,A)')
     & '*************************************************************'
              CALL HEADER('M^NU',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANNUC(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANNUC(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANNUC(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANNUC(1,1,IATOM)
     &              +GTRANNUC(2,2,IATOM)+GTRANNUC(3,3,IATOM))/3
C           Total expectation value contribution to NSR
              CALL HEADER('M^EV',-1)
              WRITE (LUPRI,'(3X,A,6X,3F16.8)')
     &                 'A', (GTRANEXP(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANEXP(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANEXP(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,F16.8)')
     &                 'iso', (GTRANEXP(1,1,IATOM)
     &              +GTRANEXP(2,2,IATOM)+GTRANEXP(3,3,IATOM))/3
              WRITE(LUPRI,'(A)')
     & '*************************************************************'
             END IF
            END IF  ! (LINMOL)
 47        CONTINUE
          END IF    ! (LINEAR)
          END IF    ! (CHARGE(IATOM).GT.D0)
         END IF     ! (DOATOM(K))
        END IF
 200   CONTINUE
 98   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_pvc_spinro */
      SUBROUTINE DEF_PVC_SPINRO(PNAME,PLABEL,PFAC,IPRINT)
C=======================================================================
C
C     Define operators for parity-violation contribution
C     to nuclear spin-rotation constants
C
C     Agustin Aucar - November 2022
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: GAUNUC, NAMDEP, NUCDEP, NUCIND, NUCPRE
      use symmet_mod, only: IPTCNT, IPTNUC, MAXREP
!      use dgroup_mod, only: LINEAR
!      use dcbprp_mod, only: IP_PVC_SPINRO
!      use dcbxlr_mod, only: LLRAPU, LLRBPU, MAXLLR, NLRAPT, NLRBPT
!      use dcbgen_mod, only: CVAL
!      use spinro_mod, only: LINDET
      use consts_mod, only: D0, D3, D4
!      use codata_mod, only: GFAC
      use pi_mod, only: PI
      use chrnos_mod, only: CHRNOS
      use memory_allocator
C
      implicit none
C
      integer I,J
C
#include "dgroup.h"
#include "dcbprp.h"
#include "dcbxlr.h"     ! Common block for linear reponse
#include "dcbgen.h"
#include "spinro.h"
#include "codata.h"
C
      INTEGER IPRINT,NATOM,IATOM,ICOOR,INDXPR,
     &        IOP,IPTYP,IREPO,ISCOR,ISYXPR,
     &        INDOP1,ITRXPR,NPCOMP
      DOUBLE PRECISION PFAC,DL,DS
      logical, allocatable :: DOATOM(:)
      CHARACTER PNAME*16, PLABEL*8
      DIMENSION IOP(3)
C
C     ANGMOM integrals correspond to -iL, so we set DL = -1
C     XSIGMA,YSIGMA and ZSIGMA come with imaginary i, so we set DS positive (DS = 0.5)
      PARAMETER(DL = -1.0D0,DS = 0.5D0)
C
      CALL QENTER('DEF_PVC_SPINRO')
C
      CALL IZERO(IP_PVC_SPINRO,6+MXCOOR)
C
      allocate(DOATOM(NUCIND))
C
C     INTTYP = 62 (Integrals for Parity Violation - chirality)
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)
C
C
C     Here we set up the operators needed for the linear reponse part
C
C     Define left-hand operators: nuclear-spin-dependent PV: c i.alpha rho(r)
C     =======================================================================
C
C     XALPHA
      IOP(1) = 2
C     YALPHA
      IOP(2) = 3
C     ZALPHA
      IOP(3) = 4
C
      NPCOMP = 1    ! one component
      PFAC   = CVAL ! factor
C
C     Point nucleus
      IF(.NOT.GAUNUC) PFAC = PFAC*D3/(D4*GFAC*PI)
C
C     Irreps
      DO IREPO = 0, MAXREP
C       Atoms
        DO IATOM = 1, NUCIND
          IF(DOATOM(IATOM)) THEN
C           Cartesian directions
            DO ICOOR = 1, 3
              ISCOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREPO,2)
              IF (ISCOR .GE. 1000) THEN
                CALL QUIT(
     &      'DEF_PVC_SPINRO: Max 1000 coords in this implementation !')
              ENDIF
              IF (ISCOR .GT. 0) THEN
                PNAME    = 'PVCa '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                            //CHRNOS(IPTNUC(IATOM,IREPO)/10)
     &                            //CHRNOS(MOD(IPTNUC(IATOM,IREPO),10))
                PNAME(5:5) = CHRNOS(ICOOR)
                IPTYP  = IOP(ICOOR)
                IF(GAUNUC) THEN ! Gaussian nuclear model
                  PLABEL= 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                         //CHRNOS(IPTNUC(IATOM,IREPO)/10)
     &                         //CHRNOS(MOD(IPTNUC(IATOM,IREPO),10))
                ELSE            ! Point nuclear model
                  PLABEL= 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                         //CHRNOS(IPTNUC(IATOM,IREPO)/10)
     &                         //CHRNOS(MOD(IPTNUC(IATOM,IREPO),10))
                ENDIF
                CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                      INDXPR,ISYXPR,ITRXPR,IPRINT)
                CALL OP1IND(PNAME(1:6),IP_PVC_SPINRO(ISCOR),
     &                      LLRAPU,NLRAPT,INDXPR,MAXLLR)
              ENDIF
            ENDDO
          END IF
        ENDDO
      ENDDO
C
C
C     Define right-hand operators: i.L and i.S
C     ========================================
C
      LINDET=.FALSE.    ! to detect if linear symmetry is used
C
      IF (LINEAR) THEN  ! linear molecule detected using symmetry
C
Ciaa   LINDET is used to record that LINEAR symmetry is used,
C      even when in further steps it could be changed, as occur
C      for UKB, spin-free and ZORA calculations.
       LINDET=.TRUE.    ! save that linear symmetry is used
C
C     Orbital angular momentum
C
       PNAME     = 'X Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC      = DL   ! factor
       PLABEL    = 'XANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XANGMOM ',IP_PVC_SPINRO(1+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
C     Spin angular momentum
C
       PNAME     = 'X Spin ang.mom. '
       IPTYP     = 10   ! Sigma_x
       NPCOMP    = 1    ! one component
       PFAC      = DS   ! factor
       PLABEL    = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XSPINMOM',IP_PVC_SPINRO(4+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
      ELSE ! if LINEAR is not TRUE
C
C     Orbital angular momentum
C
       PNAME     = 'X Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC      = DL   ! factor
       PLABEL    = 'XANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XANGMOM ',IP_PVC_SPINRO(1+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Y Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC      = DL   ! factor
       PLABEL    = 'YANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('YANGMOM ',IP_PVC_SPINRO(2+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Z Orb. ang. mom.'
       IPTYP     = 1    ! diagonal operator
       NPCOMP    = 1    ! one component
       PFAC      = DL   ! factor
       PLABEL    = 'ZANGMOM'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('ZANGMOM ',IP_PVC_SPINRO(3+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
C     Spin angular momentum
C
       PNAME     = 'X Spin ang.mom. '
       IPTYP     = 10   ! Sigma_x
       NPCOMP    = 1    ! one component
       PFAC      = DS   ! factor
       PLABEL    = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('XSPINMOM',IP_PVC_SPINRO(4+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Y Spin ang.mom. '
       IPTYP     = 11   ! Sigma_y
       NPCOMP    = 1    ! one component
       PFAC      = DS   ! factor
       PLABEL    = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('YSPINMOM',IP_PVC_SPINRO(5+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
C
       PNAME     = 'Z Spin ang.mom. '
       IPTYP     = 12   ! Sigma_z
       NPCOMP    = 1    ! one component
       PFAC      = DS   ! factor
       PLABEL    = 'OVERLAP'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,IPRINT)
       CALL OP1IND('ZSPINMOM',IP_PVC_SPINRO(6+MXCOOR),LLRBPU,NLRBPT,
     &                       INDXPR,MAXLLR)
      END IF !LINEAR
C
C
C     Print section for parity-violation contribution to nuclear spin-rotation:
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Parity-violation contribution to '//
     &  'nuclear spin-rotation constant:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     Print reference
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   '** A operators for PVC to nuclear spin-rotation **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3*NUCDEP
        INDOP1 = IP_PVC_SPINRO(I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRAPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   '** B operators for PVC to nuclear spin-rotation **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,6
        INDOP1 = IP_PVC_SPINRO(I+MXCOOR)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
C
      deallocate(DOATOM)
C
      CALL QEXIT('DEF_PVC_SPINRO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* deck prp_pvc_spinro */
      SUBROUTINE PRP_PVC_SPINRO(ABLRF,ABCNV)
C***********************************************************************
C
C     Output routine for Parity-violation contribution
C     to nuclear spin-rotation constants
C
C     Agustin Aucar - November 2022
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use dcbprp_mod, only: ABUND, IPRPRP, IP_PVC_SPINRO
!      use dcbxlr_mod, only: BFREQ, NBFREQ, NLRAPT, NLRBPT, XLR_SKIPEE, XLR_SKIPEP
      use nuclei_mod, only: IZATOM, NAMN, NATOMS, NFLOAT, NUCDEG,
     &  NUCDEP, NUCIND
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: IPTAX
!      use spinro_mod, only: LINDET
!      use dgroup_mod, only: LINEAR
!      use dcbham_mod, only: SPINFR, URKBAL, ZORA
      use consts_mod, only: D0
      use memory_allocator
C
      implicit none
C
#include "dcbprp.h"
#include "dcbxlr.h"
#include "spinro.h"
#include "dgroup.h"
#include "dcbham.h"
C
      logical, allocatable :: DOATOM(:)
      integer, allocatable :: NAT(:),ISOT(:)
      double precision, allocatable :: PMATL(:,:,:,:),PMATS(:,:,:,:),
     &                                 TMPMATL(:,:,:),TMPMATS(:,:,:),
     &                                 CSTRA(:),SCTRA(:),GEOM(:,:),
     &                                 MASS(:),AMASS(:)
C
      LOGICAL TEST
      INTEGER IBUFA(3*NUCDEP),IBUFB1(3),IBUFB2(3),
     &        IATINF(3,NUCDEP),IDEG,IDEP,IFREQ,INUC,ISO,IZ,JATOM,NATOM,
     &        I,IA,IB1,IB2,J,JA,JB1,JB2,K,KA,KB1,KB2,NA,NB1,NB2,
     &        L,NCMAT,NSCOOR,NISO,NSPEC,IABLRF_EE,IABLRF_EP
      DOUBLE PRECISION RATINF(3,5,NUCDEP),DISOTP,ABND,GVALUE,
     &                 ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &                 ABCNV(NBFREQ,NLRBPT)
C
      JATOM  = NATOMS + NFLOAT
      NCMAT  = 9*NUCDEP
      NSCOOR = 3*NUCDEP
      allocate(DOATOM(NUCIND))
      allocate(PMATL(3,3,NUCDEP,2)) ! L part of NSR for ee and ep contributions.
      allocate(PMATS(3,3,NUCDEP,2)) ! S part of NSR for ee and ep contributions.
      allocate(CSTRA(NSCOOR*NSCOOR))
      allocate(SCTRA(NSCOOR*NSCOOR))
      allocate(TMPMATL(3,NSCOOR,2))
      allocate(TMPMATS(3,NSCOOR,2))
      allocate(GEOM(NUCDEP,3))
      allocate(MASS(NATOMS))
      allocate(NAT(NATOMS))
      allocate(AMASS(NUCDEP))
      allocate(ISOT(NUCDEP))
C
      CALL QENTER('PRP_PVC_SPINRO')
C
      CALL TITLER('PARITY-VIOLATION CONTRIBUTION TO '//
     &  'NUCLEAR SPIN-ROTATION TENSOR','*',116)
C
      WRITE(LUPRI,'(13X,A)')
     & 'According to I. A. Aucar et al., JCP 155 (2021) p.134307'
      CALL PRSYMB(LUPRI,'-',81,0)
      WRITE(LUPRI,*) ''
C
C     INTTYP = 62 (Integrals for Parity Violation - chirality)
C     (-62: accept .SELECT specification)
      CALL SETATM(DOATOM,NATOM,-62)
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &  '@ name  charge  isotope   mass         abundance     g factor'
      CALL PRSYMB(LUPRI,'-',61,0)
      NSPEC = 0
      IDEP  = 0
      DO INUC = 1,NUCIND
      IDEG = NUCDEG(INUC)
      IF(DOATOM(INUC)) THEN
        TEST = .FALSE.
        IZ   = IZATOM(INUC)
        NISO = 0
        DO ISO = 1,5
            GVALUE = DISOTP(IZ,ISO,'GVAL')
            ABND = DISOTP(IZ,ISO,'ABUNDANCE')
            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = GVALUE
            ENDIF
        ENDDO
        IF(NISO.GT.0) THEN
          IATINF(1,NSPEC+1) = IDEP + 1
          IATINF(2,NSPEC+1) = IZ
          IATINF(3,NSPEC+1) = NISO
          WRITE(LUPRI,
     &      '(A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &       '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
          DO J = 2,NISO
            WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F10.6)')
     &         '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
          ENDDO
          DO K = 2,IDEG
            IATINF(1,NSPEC+K) = IDEP + K
            IATINF(2,NSPEC+K) = IZ
            IATINF(3,NSPEC+K) = NISO
            DO L = 1,NISO
              RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
              RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
              RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
            ENDDO
          ENDDO
          NSPEC = NSPEC + IDEG
        ENDIF
      ENDIF
      IDEP = IDEP + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        CALL QEXIT('PRP_PVC_SPINRO')
        RETURN
      ENDIF
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C     **********************************************************************
C     ***** Extract PV contribution to nuclear spin-rotation constants *****
C     **********************************************************************
C
C     The NUCDEP PV operators as A operators in <<A;B>> :
      NA = 0
      DO I = 1,NSCOOR
        IF(IP_PVC_SPINRO(I).GT.0) THEN
          NA = NA + 1
          IBUFA(NA) = I
        ENDIF
      END DO
      IF (NA.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_SPINRO: NA=',NA
        CALL QUIT('PRP_PVC_SPINRO: Wrong value of NA!')
      ENDIF
C
Ciaa   As in the current implementation UKB, spin-free and ZORA
C      calculations of energy cannot be performed employing symmetry
C      (they set LINEAR to false in such cases, see subroutine GLINSM)
C      LINDET allows to recover LINEAR symmetry defined in DEF_PVC_SPINRO
C      if it is lost due to the use of UKB, SF or ZORA functionalities
      IF ((URKBAL.or.SPINFR.or.ZORA).and.LINDET) LINEAR=.TRUE.
C
C     The 3 L orbital angular momentum directions as B operators in <<A;B>> :
      NB1 = 0
      DO I = 1,3
        IF(IP_PVC_SPINRO(I+MXCOOR).GT.0) THEN
          NB1 = NB1 + 1
          IBUFB1(NB1) = I
        ENDIF
      END DO
      IF (NB1.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_SPINRO: NB1=',NB1
        CALL QUIT('PRP_PVC_SPINRO: Wrong value of NB1!')
      ENDIF
C
C     The 3 S spin angular momentum directions as B operators in <<A;B>> :
      NB2 = 0
      DO I = 1,3
        IF(IP_PVC_SPINRO(I+3+MXCOOR).GT.0) THEN
          NB2 = NB2 + 1
          IBUFB2(NB2) = I
        ENDIF
      END DO
      IF (NB2.LE.0) THEN
        WRITE(LUPRI,'(2X,A,I2)') 'PRP_PVC_SPINRO: NB2=',NB2
        CALL QUIT('PRP_PVC_SPINRO: Wrong value of NB2!')
      ENDIF
C
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C--------------------------------------------------------------------
C
      CALL CMMASS(GEOM,MASS,NAT,ISOT,IPRPRP)
C
C     LINEAR RESPONSE PART OF PVC-NSR TENSOR
C
C     Fill matrices of parity-violation contribution to nuclear spin-rotation
C
      IF (XLR_SKIPEE) THEN
        IABLRF_EE = -999999999
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 2
        END IF
      ELSE
        IABLRF_EE = 2
        IF (XLR_SKIPEP) THEN
          IABLRF_EP = -999999999
        ELSE
          IABLRF_EP = 3
        END IF
      END IF
C
C
      DO IFREQ = 1,NBFREQ
        IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
           CALL PRSYMB(LUPRI,'-',72,2)
           WRITE(LUPRI,'(A,F15.8,A)')
     &          '*** Frequency    :',BFREQ(IFREQ),' a.u.'
        ENDIF
        CALL DZERO(TMPMATL,2*NCMAT)
        CALL DZERO(TMPMATS,2*NCMAT)
C
        DO IB1 = 1,NB1
          JB1 = IBUFB1(IB1)
          KB1 = IP_PVC_SPINRO(MXCOOR+JB1)
          JB1 = IPTAX(JB1,2)
          DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IP_PVC_SPINRO(JA)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         TMPMATL(JB1,JA,1) = ABLRF(IFREQ,KA,KB1,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         TMPMATL(JB1,JA,2) = ABLRF(IFREQ,KA,KB1,IABLRF_EP)
          ENDDO
        ENDDO
C
        DO IB2 = 1,NB2
          JB2 = IBUFB2(IB2)
          KB2 = IP_PVC_SPINRO(MXCOOR+3+JB2)
          JB2 = IPTAX(JB2,2)
          DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IP_PVC_SPINRO(JA)
C           First extract the e-e contribution
            IF (.NOT. XLR_SKIPEE)
     &         TMPMATS(JB2,JA,1) = ABLRF(IFREQ,KA,KB2,IABLRF_EE)
C           then extract the e-p contribution
            IF (.NOT. XLR_SKIPEP)
     &         TMPMATS(JB2,JA,2) = ABLRF(IFREQ,KA,KB2,IABLRF_EP)
          ENDDO
        ENDDO
C
C     Transform to C1 basis
C
        CALL DZERO(PMATL,2*NCMAT)
        CALL DZERO(PMATS,2*NCMAT)
C
        IF (.NOT. XLR_SKIPEE) THEN
C           .... extract the e-e response part into PMAT(...,1)
            CALL TRADIP(TMPMATL,PMATL,CSTRA,SCTRA,NSCOOR,2,2)
            CALL TRADIP(TMPMATS,PMATS,CSTRA,SCTRA,NSCOOR,2,2)
        END IF
        IF (.NOT.XLR_SKIPEP) THEN
C           ... extract the e-p response part into PMAT(...,2)
            CALL TRADIP(TMPMATL(1,1,2),PMATL(1,1,1,2),
     &                CSTRA,SCTRA,NSCOOR,2,2)
            CALL TRADIP(TMPMATS(1,1,2),PMATS(1,1,1,2),
     &                CSTRA,SCTRA,NSCOOR,2,2)
        END IF
      END DO
C
C--------------------------------------------------------------------
      CALL PVCNSR(GEOM,AMASS,PMATL,PMATS,JATOM,IPRPRP,DOATOM)
C--------------------------------------------------------------------
C
      CALL QEXIT('PRP_PVC_SPINRO')
C
      deallocate(AMASS)
      deallocate(ISOT)
      deallocate(GEOM)
      deallocate(MASS)
      deallocate(NAT)
      deallocate(DOATOM)
      deallocate(PMATL)
      deallocate(PMATS)
      deallocate(CSTRA)
      deallocate(SCTRA)
      deallocate(TMPMATL)
      deallocate(TMPMATS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pvcnsr */
      SUBROUTINE PVCNSR(GEOM,AMASS,PMATL,PMATS,KATOM,IPRINT,DOATOM)
C***********************************************************************
C
C     Analyze routine for the parity-violation contribution
C     to nuclear spin-rotation constants
C
C     This routine is modelled after SPRRES subroutine
C
C   INPUT:
C     PMATL - contribution from orbital angular moment to linear response term
C     PMATS - contribution from spin angular moment to linear response term
C
C     Agustin Aucar                - November 2022
C     Last revision: Agustin Aucar - March    2024
C
C***********************************************************************
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use nuclei_mod, only: CHARGE, CORD, ISOTOP, ISTBNU, IZATOM,
     &  NAMDEP, NUCDEP, NUCIND
      use symmet_mod, only: IPTAX, ISYMAX, MAXOPR, PT
!      use spinro_mod, only: GTRANT
      use orgcom_mod, only: CMXYZ
!      use codata_mod, only: CVEL, GFERMI, S2THETAW, XFAMU, XTHZ
      use consts_mod, only: D0, D1, D2, D4
      implicit none
C
      integer I,J
C
#include "spinro.h"
#include "codata.h"
C
      integer K,KATOM,IPRINT,IATOM,ICOOR,ISYMOP,IX,IY,
     &        JATOM,JX,JY,NATTYP,NUMIS
C
      double precision FACTOR,CONST,GEOM(KATOM,3),AMASS(KATOM),
     &                 EIGVAL(3),EIGVEC(3,3),LEIGVAL(3),LEIGVEC(3,3),
     &                 OMEGA(3),TINERT(9),DISOTP,ANGUMOM(3),
     &                 PMATL(3,3,KATOM,2),PMATS(3,3,KATOM,2),
     &                 NSRTOT,NSRLRLEE,NSRLRLEP,NSRLRSEE,NSRLRSEP,
     &                 GTRANPLEE(3,3,KATOM),GTRANPLEP(3,3,KATOM),
     &                 GTRANPSEE(3,3,KATOM),GTRANPSEP(3,3,KATOM)
C
      LOGICAL PLANAR,LINMOL,DOATOM(NUCIND),LINNOZ
C
      FACTOR = (-GFERMI*(D1-D4*S2THETAW)/(CVEL*SQRT(D2)))*(XTHZ/XFAMU)
C
C
      LINNOZ = .FALSE. ! Only true for linear molecules, and when they are along the z axis
C
Ciaa  Nuclear rotational angular momentum could be defined here to obtain Omega afterwards
      DO 300 I = 1,3
        ANGUMOM(I) = D1
 300  CONTINUE
C
C
C     Determine principal moments of inertia
C     ======================================
C
C     Set up coordinates (relative to center of mass), g-values and masses of all nuclei
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         NATTYP = IZATOM(IATOM)
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
         IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
           JATOM = JATOM + 1
           AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')
           DO 120 ICOOR = 1, 3
             GEOM(JATOM,ICOOR) = PT(IAND(ISYMAX(ICOOR,1),ISYMOP))
     &                          *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
 120       CONTINUE
         END IF
 110     CONTINUE
 100  CONTINUE
C
      CALL WLKDIN(GEOM,AMASS,KATOM,ANGUMOM,TINERT,OMEGA,EIGVAL,EIGVEC,
     &            .TRUE.,PLANAR,LINMOL,IPRINT)
C     On output
C        TINERT - inertia tensor
C        EIGVAL - inverse eigenvalues (for linear molecules, in descending order)
C        EIGVEC - corresponding eigenvectors
C        PLANAR=T : planar molecule
C        LINMOL=T : linear molecule (EIGVAL(1) = 0.0)
C
C     Transform each PVC to NSR constant to principal axis system,
C     and multiply by inverse moment of inertia and various constants
C
      IF(LINMOL) THEN
        LEIGVAL(1)=EIGVAL(3)
        LEIGVAL(3)=EIGVAL(1)
        DO I=1,3
          LEIGVEC(I,1)=EIGVEC(I,3)
          LEIGVEC(I,3)=EIGVEC(I,1)
        END DO
        EIGVAL(1)=LEIGVAL(1)
        EIGVAL(3)=LEIGVAL(3)
        DO I=1,3
          EIGVEC(I,1)=LEIGVEC(I,1)
          EIGVEC(I,3)=LEIGVEC(I,3)
        END DO
        IF((EIGVEC(1,1).NE.D1).OR.
     &     (EIGVEC(2,2).NE.D1).OR.
     &     (EIGVEC(3,3).NE.D1)) THEN
          LINNOZ=.TRUE.
        END IF
      END IF
C
      IF(LINNOZ) THEN
        WRITE (LUPRI,'(/,2X,A)') '*** WARNING: In the current'//
     & ' implementation of PVCSR, the linear molecule should'//
     & ' be placed along the z axis. ***'
      ELSE
Ciaa  To obtain Omega, the nuclear rot. angular momentum (ANGUMOM) should be defined before
      IF (IPRINT .GE. 8) THEN
       CALL HEADER
     &    ('Nuclear angular velocity (a.u.)',0)
       WRITE (LUPRI,'(3X,A,2X,F16.8)')
     &       'Omega_A',OMEGA(1),
     &       'Omega_B',OMEGA(2),
     &       'Omega_C',OMEGA(3)
      END IF
C
C    Print inertia tensor or moment
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER
     &    ('Principal moments of inertia (a.u.) and principal axes',0)
         WRITE (LUPRI,'(3X,A,2X,F16.8,6X,3F16.8)')
     &       'IA',XFAMU/EIGVAL(1), (EIGVEC(1,I),I = 1,3),
     &       'IB',XFAMU/EIGVAL(2), (EIGVEC(2,I),I = 1,3),
     &       'IC',XFAMU/EIGVAL(3), (EIGVEC(3,I),I = 1,3)
         WRITE (LUPRI,'(/,2X,A)') 'PVC to nuclear spin-rotation '     //
     &         'tensor elements are given in the principal axis system.'
      END IF
C
C    Print PVC to nuclear spin-rotation constants or tensors
C
      IATOM = 0
      DO 98 K = 1, NUCIND
       DO 200 ISYMOP = 0, MAXOPR
        IF (IAND(ISTBNU(K),ISYMOP).EQ.0) THEN
         IATOM = IATOM + 1
         IF (DOATOM(K)) THEN
          IF (CHARGE(IATOM) .GT. D0) THEN
           IF (IPRINT .GE. 0) THEN
            CALL HEADER('Parity-violation contribution to '    //
     &            'nuclear spin-rotation constant (in Hz) for '//
     &                     NAMDEP(IATOM),-1)
            WRITE (LUPRI,'(2X,A14,F12.6,/)')
     &                 'Nuclear mass: ', AMASS(IATOM)
           END IF
C
            CALL DZERO(GTRANT,9*NUCDEP)
            CALL DZERO(GTRANPLEE,9*NUCDEP)
            CALL DZERO(GTRANPLEP,9*NUCDEP)
            CALL DZERO(GTRANPSEE,9*NUCDEP)
            CALL DZERO(GTRANPSEP,9*NUCDEP)
C
            DO 50 I = 1, 3
             DO 50 J = 1, 3
              DO 60 IX = 1, 3
               DO 60 IY = 1, 3
                JX = IPTAX(IX,2)
                JY = IPTAX(IY,2)
                CONST = FACTOR*EIGVAL(I)
C           ee and pp parts of the PVC-NSR
                GTRANPLEE(I,J,IATOM) = GTRANPLEE(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATL(JX,JY,IATOM,1)*EIGVEC(IY,J)*CONST
                GTRANPLEP(I,J,IATOM) = GTRANPLEP(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATL(JX,JY,IATOM,2)*EIGVEC(IY,J)*CONST
                GTRANPSEE(I,J,IATOM) = GTRANPSEE(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATS(JX,JY,IATOM,1)*EIGVEC(IY,J)*CONST
                GTRANPSEP(I,J,IATOM) = GTRANPSEP(I,J,IATOM)
     &                       + EIGVEC(IX,I)
     &                        *PMATS(JX,JY,IATOM,2)*EIGVEC(IY,J)*CONST
 60           CONTINUE
C           Total PV contribution to NSR
              GTRANT(I,J,IATOM) = GTRANPLEE(I,J,IATOM)
     &                       + GTRANPLEP(I,J,IATOM)
     &                       + GTRANPSEE(I,J,IATOM)
     &                       + GTRANPSEP(I,J,IATOM)
 50         CONTINUE
C
             IF (IPRINT .GE. 0) THEN
C
C           Total PVC-NSR for molecules in equilibrium
C
              CALL HEADER('Total PVC to nuclear spin-rotation '//
     &              'constant',-1)
              WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'A', (GTRANT(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANT(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANT(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANT(1,1,IATOM)
     &              +GTRANT(2,2,IATOM)+GTRANT(3,3,IATOM))/3
             END IF
C
C           Linear response contributions to PVC-NSR
C
             IF (IPRINT .GE. 1) THEN
              CALL HEADER('M^LR-L(e-e)',-1)
              WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'A', (GTRANPLEE(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPLEE(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPLEE(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANPLEE(1,1,IATOM)
     &              +GTRANPLEE(2,2,IATOM)+GTRANPLEE(3,3,IATOM))/3
              CALL HEADER('M^LR-S(e-e)',-1)
              WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'A', (GTRANPSEE(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPSEE(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPSEE(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANPSEE(1,1,IATOM)
     &              +GTRANPSEE(2,2,IATOM)+GTRANPSEE(3,3,IATOM))/3
              CALL HEADER('M^LR-L(e-p)',-1)
              WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'A', (GTRANPLEP(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPLEP(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPLEP(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANPLEP(1,1,IATOM)
     &              +GTRANPLEP(2,2,IATOM)+GTRANPLEP(3,3,IATOM))/3
              CALL HEADER('M^LR-S(e-p)',-1)
              WRITE (LUPRI,'(3X,A,6X,3E16.8)')
     &                 'A', (GTRANPSEP(1,I,IATOM),I = 1,3),
     &                 'B', (GTRANPSEP(2,I,IATOM),I = 1,3),
     &                 'C', (GTRANPSEP(3,I,IATOM),I = 1,3)
              WRITE (LUPRI,'(/,3X,A,4X,E16.8)')
     &                 'iso', (GTRANPSEP(1,1,IATOM)
     &              +GTRANPSEP(2,2,IATOM)+GTRANPSEP(3,3,IATOM))/3
             END IF
C
          END IF    ! (CHARGE(IATOM).GT.D0)
         END IF     ! (DOATOM(K))
        END IF
 200   CONTINUE
 98   CONTINUE
      END IF        ! LINNOZ
      RETURN
      END
