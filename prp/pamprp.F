!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C FILE    : pamprp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamprp */
      SUBROUTINE PAMPRP(WORK,LWORK)
C***********************************************************************
C
C     Driver for PAMPR1.
C
C     Loop over wave functions and call PAMPR1 for each of them.
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
C Used in common blocks:
C dcbprp.h: NPRPWF, PRPWF(*)
C
#include "dcbprp.h"
C
      DIMENSION WORK(LWORK)
C
      CHARACTER WF*4
C
      CALL QENTER('PAMPRP')
#include "memint.h"
C
      CALL TITLER('Property module','*',125)
      WRITE(LUPRI,9000)
 9000 FORMAT(//,
     &     3X,'This is output from the Dirac property module:',/,
     &    /3X,'HF & DFT first order properties',/,
     &     3X,'Trond Saue',//,
     &     3X,'First-order ESR properties',/,
     &     3X,'Hans Joergen Aa. Jensen and Anette Noerager',//,
     &     3X,'MP2 first order properties:',/,
     &     3X,'J. N. P. van Stralen, L. Visscher, C. V. Larsen' 
     &       ,' and H. J. Aa Jensen, Chem. Phys. 311 (2005) 81.',//
     &     3X,'KR-RPA second-order properties',/,
     &     3X,'Hans Joergen Aa. Jensen and Trond Saue',//,
     &     3X,'KR-QR third order properties',/,
     &     3X,'Patrick Norman and Hans Joergen Aa. Jensen',//,
     &     3X,'Molecular gradient',/,
     &     3X,'Joern Thyssen',//,
     &     3X,'Additional contributions from:',/,
     &     3X,'Thomas Enevoldsen, Miroslav Ilias (London orbitals)'//)
C
      DO I = 1, NPRPWF
         WF = PRPWF(I)
C
         IF (WF .EQ. 'DHF' .OR. WF .EQ. 'DFT') THEN
            CALL TDHFORB('FULLCI')
            CALL SETDC2(0)
         ELSE IF (WF .EQ. 'KRMC') THEN
            CALL TKRMCORB
            CALL SETDC2(0)
         ELSE IF (WF .EQ. 'MP2 ') THEN
            CALL SETDC2(0)
         ELSE
            CALL QUIT('*** ERROR in PAMPRP *** ' //
     &           'Unknown wave function: ' // WF)
         END IF
         CALL TITLER('Properties for ' // WF // ' wave function',
     &        '*',110)
         CALL PAMPR1(WF,WORK,LWORK)
      END DO
C
      CALL QEXIT('PAMPRP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampr1 */
      SUBROUTINE PAMPR1(WF,WORK,LWORK)
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
C     The code (in whole or part) is not to be reproduced for
C     further distribution without the written permission of the authors.
C
C     If results obtained with this code are published, an appropriate
C     citation would be:  H. J. Aa. Jensen and T. Saue,
C     "XXX version 1.0 1995,
C      a Dirac-Fock molecular linear response properties program".
C
C     This note is not to be removed or altered.
C
C     Written by T.Saue June 24 1996
C     The last revision: 1997/09/15 - jth (addition of gradient code)
C
C***********************************************************************

      use dirac_cfg
#ifdef MOD_OPENRSP
      use openrsp_main
#endif

!     USE FDE_MOD, ONLY: DO_GRID_OUT,FDE_GRIDOUT
      USE FDE_MOD ! MI: gfortran-4.2 does not accept "ONLY" (yet)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "dcborb.h"
#include "dcbpsi.h"
C
      DIMENSION WORK(LWORK)
      CHARACTER WF*4
      CHARACTER WFEXPORT*4
C
      CALL QENTER('PAMPR1')
C
#include "memint.h"
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
         IF (WF .NE. 'DHF') THEN
            CALL QUIT('*** ERROR in PAMPR1 *** ' //
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
         CALL PRPEXP(WF,WORK(KEXP),WORK(KFREE),LFREE)
      ELSE
         KEXP = KFREE
      END IF
C
C***********************************************************************
C*****  E S R   P R O P E R T I E S      *******************************
C***********************************************************************
C
#ifdef MOD_ESR
      IF (DOESR) THEN
         KFRSAV = KFREE
         CALL ESR_INIT(WF,IPRESR,N_CIESR,SPINPOLESR,WORK,KFREE,LFREE)
C
         LENESR = 2*MULTIPESR*MULTIPESR*MAX(4,NESRP)
C        max(4,nesrp): make sure sufficient space for ESRPHASE
         CALL MEMGET('REAL',KESR,LENESR,WORK,KFREE,LFREE)
         CALL PRPESR(WF,WORK(KESR),WORK,KFREE,LFREE)
         CALL MEMREL('after PRPESR',WORK,KWORK,KFRSAV,KFREE,LFREE)
C        ... remove this MEMREL if WORK(KESR) to be used in PRPOUT below
      END IF
#endif
C
C***********************************************************************
C*****  P P   E X C I T A T I O N   E N E R G I E S  *******************
C***********************************************************************
C
      IF (DOXPP) THEN
         IF (WF. EQ. 'DHF') THEN
            CALL PRPXPP(WORK,KPPF,KOMEGA,KPPCNV,KPPO,KSSYM,KFREE,LFREE)
         ELSE
            CALL QUIT('*** ERROR in PAMPR1 *** ' //
     &           'PP excitation energies not implemented for ' // WF )
         END IF
      END IF
C
C***********************************************************************
C*****  L I N E A R   R E S P O N S E  *********************************
C***********************************************************************
C
      IF(DOXLR) THEN
         IF (WF. EQ. 'DHF') THEN
            CALL PRPXLR(WORK,KLRF,KLRCNV,KLRO,KFREE,LFREE)
         ELSE
            CALL QUIT('*** ERROR in PAMPR1 *** ' //
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
         IF (WF. NE. 'DHF') 
     &        CALL QUIT('*** ERROR in PAMPR1 *** ' //
     &        'Quadratic response not implemented for ' // WF )
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
C JHS 17 apr 2008
C
      call get_xml_export_level(wfexport)
      IF (DO_GRID_OUT.and.(WF.eq.wfexport)) THEN
#ifdef MOD_XML
         CALL FDE_GRIDOUT(WORK,LFREE)
#else
         call quit("FDE density export unsupported in this version.")
#endif
      END IF

#ifdef MOD_OPENRSP
      if (dirac_cfg_openrsp) then
!       main openrsp driver
        call openrsp_driver(work(kfree), lfree)
!       initialize dirac interface
        call init_dirac_interface()
      end if
#endif

C
CMI   CALL MEMREL('PAMPRP',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PAMPR1')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_dipol */
      SUBROUTINE DEF_DIPOLE(INDDIP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define dipole moment operator
C     Written by T.Saue - May 1996
C     Last revision: May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(DM1 = -1.0D0)
#include "dcbexp.h"
#include "dcbxpr.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION INDDIP(3),WORK(LWORK),PFAC(3)
C
C
      PNAME     = 'Dipole length: X'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDDIP(1),ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
C
      PNAME      = 'Dipole length: Y'
      IPTYP      = 1
      NPCOMP     = 1
      PFAC(1)    = DM1
      PLABEL(1)  = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDDIP(2),ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
C
      PNAME      = 'Dipole length: Z'
      IPTYP      = 1
      NPCOMP     = 1
      PFAC(1)    = DM1
      PLABEL(1)  = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &            INDDIP(3),ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_QUADRU(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define quadrupole moment operator
C     Written by T.Saue - May 1996
C     Last revision: May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D1 = 1.0D0)
#include "dcbxpr.h"
#include "dcbexp.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XXQUAD',IPQUADRU(1),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: XY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XYTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XYQUAD',IPQUADRU(2),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: XZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XZTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XZQUAD',IPQUADRU(3),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: YY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YYTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YYQUAD',IPQUADRU(4),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: YZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YZTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YZQUAD',IPQUADRU(5),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'Quadrupol m.: ZZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'ZZTHETA '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_EPOLAR(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operators for electric dipole polarizabilities
C    
C     This is only linear response calculation.
C
C     Written by T.Saue - June 28 1996
C     Last revision: June 28 - tsaue
C                    July 2006 - M.Ilias, Prievidza
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D1 = 1.0D0)
#include "dcbprp.h"
#include "dcbxlr.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
#include "dcbibt.h"
C
      CALL QENTER('DEF_EPOLAR')

      PNAME     = 'X-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_EFN(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for electric field 
C
C     Written by T.Saue Feb 26 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     EFNINT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-29)
      CALL DEF_EFN1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_EFN',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_efn1 */
      SUBROUTINE DEF_EFN1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Define electric field at individual nuclei
C     Based on EFNTYP
C     Written by T. Saue Feb 26 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
#include "symmet.h"
#include "pgroup.h"
#include "chrnos.h"
#include "chrxyz.h"
#include "nuclei.h"
#include "dcbxpr.h"
#include "dcbexp.h"
      PARAMETER(D1 = 1.0D0)
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   DOATOM(NUCIND)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_EFG(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for electric field gradient
C
C     Written by T.Saue Feb 19 1997
C     Last revision: J. Thyssen - Oct 13 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     EFGINT: INTTYP = 30 (cart. el. field grad. at indiv. nuclei)
C     (-30: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-30)
      CALL DEF_EFG1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_EFG',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_efg */
      SUBROUTINE DEF_EFG1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define electric field gradient
C     Written by T.Enevoldsen - Sep 1996
C     Last revision: J. Thyssen - Oct 13 1998

C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0)
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   DOATOM(NUCIND)
#include "dcbxpr.h"
#include "dcbexp.h"
C
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "chrxyz.h"
#include "ibtfun.h"
C
      ITYP = 0
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
         DO 200 ICOOR1 = 1, 3
         DO 200 ICOOR2 = ICOOR1, 3
            ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IBTAND(ISTBNU(IATOM),IBTXOR(IREPC,ISYMIJ)).EQ.0) THEN
                  IOFF = IOFF + 1
                  ITYP = ITYP + 1
                  PNAME        = 'EFG: '//CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)
     *                 //CHRNOS(IATOM/10)
     &                 //CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  PLABEL(1) = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     &                 CHRNOS(IATOM/10)//
     &                 CHRNOS(MOD(IATOM,10))//CHRNOS(IOFF)
                  CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_PVC(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for parity violation - chirality
C
C     Written by J. Thyssen - Oct 7 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     PVC: INTTYP = 62
C     (-62: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-62)
      CALL DEF_PVC1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_PVC',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_pvc */
      SUBROUTINE DEF_PVC1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for parity violation - chirality
C
C     Written by J. Thyssen - Oct 7 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
#include "gfac.h"
#include "pi.h"
      PARAMETER(D1 = 1.0D0,D3 = 3.0D0,D4 = 4.0D0)
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   DOATOM(NUCIND)
#include "dcbxpr.h"
#include "dcbexp.h"
C
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "chrxyz.h"
#include "ibtfun.h"
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
               IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_EFT(PNAME,PLABEL,PFAC,TRLESS,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for electric field third derivatives
C
C     TRLESS: traceless or non-traceless
C
C     Written by: J. Thyssen - Sept 30 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   TRLESS
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     EFTINT: INTTYP = 52 or INTTYP = 53
      IF (TRLESS) THEN
         CALL SETATM(WORK(KDATOM),NATOM,52)
      ELSE
         CALL SETATM(WORK(KDATOM),NATOM,53)
      END IF
      CALL DEF_EFT1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                TRLESS,WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_EFT',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_efg */
      SUBROUTINE DEF_EFT1(PNAME,PLABEL,PFAC,DOATOM,TRLESS,
     &                    WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for fourth derivative nuclear attraction
C     integrals.
C
C     Written by: J. Thyssen - Sept 30 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0)
      CHARACTER PNAME*16, PLABEL(3)*8,
     &          PNPREFIX*6, PLBNAME*1
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   DOATOM(NUCIND),TRLESS
#include "dcbxpr.h"
#include "dcbexp.h"
Ctec
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "chrxyz.h"
#include "ibtfun.h"
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
            ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            ISYMIJ = IBTXOR(ISYMIJ          ,ISYMAX(ICOOR3,1))
            ISYMIJ = IBTXOR(ISYMIJ          ,ISYMAX(ICOOR4,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IBTAND(ISTBNU(IATOM),IBTXOR(IREPC,ISYMIJ)).EQ.0) THEN
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      SUBROUTINE DEF_SPNSPN(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for indirect spin-spin coupling
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-29)
      CALL DEF_SPNSP1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_SPNSPN',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_SPNSP1 */
      SUBROUTINE DEF_SPNSP1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for indirect spin-spin coupling
C
C     Written by T.Saue Feb 19 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
#include "chrxyz.h"
      LOGICAL DOATOM(NUCIND),SAME
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"
C
C
      NSCOOR = 3*NUCDEP
      CALL IZERO(IPSPNSPN,2*NSCOOR)
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
C
C       First define left-hand operators
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
              IREP1    = IBTXOR(IREPO,ISYMAX(IC2,1))
              IREP2    = IBTXOR(IREPO,ISYMAX(IC1,1))
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
     &           INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
              CALL OP1IND('XSPNA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
              IPSPNSPN(1,ISCOR1) = INDOP1
C
C             Right hand operators
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
C     Define right hand operators
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
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
         IF (SHIELD) THEN
            DOEXP = .TRUE.
            NSTDIA = .TRUE.
         END IF
C      .... spin-spin coupling
         IF (SPNSPN) THEN
            DOEXP = .TRUE.
            DSO = .TRUE.
         END IF
C      .... magnetic susceptibilities
#ifdef MOD_MAGNETIZ
         IF (MSUSCP) THEN
            DOEXP = .TRUE.
            MSUSCDIA = .TRUE.
         END IF
#endif
      END IF
C
C     Dipole moment
C
      IF(DIPOLE) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Dipole moment'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_DIPOLE(INDDIP,WORK,LWORK,IPRPRP)
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
        CALL DEF_QUADRU(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Electric fields at nuclei
C
      IF(EFN) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electric field at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFN(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Electric field gradients at nuclei
C
      IF(EFG) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electric field gradients'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFG(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Electronic density at nuclei
C
      IF(RHONUC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electronic density at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_RHONUC(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Parity violation - Chirality
C
      IF(PVC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Parity violation - Chirality'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_PVC(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Non-traceless electric field third derivatives
C
      IF(EFTNTL) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &     'Non-traceless electric field third derivatives'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFT(PNAME,PLABEL,PFAC,.FALSE.,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Traceless electric field third derivatives
C
      IF(EFT) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &     'Traceless electric field third derivatives'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFT(PNAME,PLABEL,PFAC,.TRUE.,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Nuclear quadrupole coupling
C
      IF(NQCC) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Nuclear quadrupole moments'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFG(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Electric polarizabilities
C
      IF(EPOLAR) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Electric dipole polarizabilities'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EPOLAR(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Diamagnetic spin-orbit
C
      IF(DSO) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &    'NMR diamagnetic spin-orbit contribution:'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_DSO(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Indirect spin-spin coupling
C
      IF(SPNSPN) THEN
        CALL DEF_SPNSPN(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF
C
C     Magnetizabilities
C
#ifdef MOD_MAGNETIZ
      IF(MSUSCP) THEN
        CALL DEF_MSUSCP(WORK,LWORK)
      ENDIF
#endif
C
C     Nuclear shielding or ESR g-tensor or ESR anisotropy
C
      IF (SHIELD .OR. ESRGTENS .OR. ESRANISO) THEN
        CALL DEF_SHIELD(WORK,LWORK,IPRPRP)
      ENDIF
C
C     Optical rotation
C
      IF (OPTROT) THEN
        CALL DEF_OPTROT(WORK,LWORK,IPRPRP)
      ENDIF
C
C
C     PV contribution to NMR shielding
      IF(PVC_SHIELD) CALL DEFINE_PVC_SHIELD(WORK,LWORK)
C     PV contribution to spin-spin coupling
      IF(PVC_SPINSPIN) CALL DEFINE_PVC_SPINSPIN(WORK,LWORK)
C
C
C     Diamagnetic contribution to shieldings
C
      IF(NSTDIA) THEN
        CALL DEF_NSTDIA(PNAME,PLABEL,PFAC,WORK,LWORK,IPREXP)
      ENDIF
C
      IF (ROTG) THEN
       CALL DEF_ROTG(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
      ENDIF

! gosia:
! calculate << B_i; B_i >>_0 with LAO
! only for one component of magnetic field, i = x, y or z 
! i need it temporarily in visual
      IF (BXLAO .OR. BYLAO .OR. BZLAO) THEN
        CALL DEF_BLAO(WORK, LWORK, IPRPRP)
      ENDIF
C
C     Effective electronic densities at the nuclei
C
      IF(EFFDEN) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)') 'Effective electronic density at nuclei'
        CALL PRSYMB(LUPRI,'=',75,0)
        CALL DEF_EFF_DENS(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
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
      IF(SHIELD) THEN
        CALL PRP_SHIELD(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
C
C
C     PV contribution to NMR shielding
      IF(PVC_SHIELD) 
     &  CALL PRINTOUT_PVC_SHIELD(ABLRF,ABCNV,WORK,KFREE,LFREE)
C     PV contribution to spin-spin coupling
      IF(PVC_SPINSPIN)
     &  CALL PRINTOUT_PVC_SPINSPIN(ABLRF,ABCNV,WORK,KFREE,LFREE)
C
C
C     Magnetic susceptibilities
C
#ifdef MOD_MAGNETIZ
      IF (MSUSCP) THEN
        CALL PRP_MSUSCP(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
#endif
C
C     Molecular rotational g-tensor
C
      IF (ROTG) THEN
C       CALL PRP_ROTG(EXPVAL,ABLRF,ABCNV,WORK,KFREE,LFREE)
        CALL PRP_ROTG(ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
C
C     Optical rotation
C
      IF (OPTROT) THEN
        CALL PRP_OPTROT(ABLRF,ABCNV,WORK,KFREE,LFREE)
      ENDIF
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
#include "implicit.h"
#include "priunit.h"
#include "debye.h"
C
#include "mxcent.h"
#include "dipole.h"
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
      DO I = 1,3
         WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &        CHAR(IOFF+I),(DEBYE*DIPME(I)),' Debye',
     &                   (DEBYE*DIPMN(I)),' Debye',
     &        (DEBYE*(DIPME(I)+DIPMN(I))),' Debye'
      ENDDO
      CALL PRSYMB(LUPRI,'-',76,0)
      WRITE(LUPRI,'(/A,F12.8,A)') ' 1 a.u = ',DEBYE,' Debye'
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL CORR
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      CALL MEMGET('LOGI',KATOM,NUCIND  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAXIS,9*MXCENT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,9*MXCENT,WORK,KFREE,LFREE)
C
      CALL SETATM(WORK(KATOM),NATOM,-30)
      CALL PRP_NQCC1(EXPVAL,EXP_CORR,CORR,WORK(KATOM),WORK(KAXIS),
     &               WORK(KTEMP))
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
      PARAMETER(CONV = 234.9647417D0)
      DIMENSION EXPVAL(*),AXIS(*),TEMP(*),PVAL(9),TMP(3),ITMP(3)
      DIMENSION EXP_CORR(*)
      LOGICAL   DOATOM(*), CORR
      CHARACTER*4 HZUNIT
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
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
      WRITE(LUPRI,'(A/)')
     &'* Nuclear quadrupole coupling constants (principal components):'
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
          CALL JACO(PVAL,AXIS,3,3,3,TMP,ITMP)
          EFGTOT(1,1,IATOM) = PVAL(1)
          EFGTOT(2,2,IATOM) = PVAL(3)
          EFGTOT(3,3,IATOM) = PVAL(6)
C
C         Reorder to calculate eta:
C            |xx|<|yy|<|zz|
C
          PVAL(2) = PVAL(3)
          PVAL(3) = PVAL(6)
          IZ = IDAMAX(3,PVAL(1),1)
          CALL DSWAP(1,PVAL(3),1,PVAL(IZ),1)
          IX = IDAMIN(3,PVAL(1),1)
          CALL DSWAP(1,PVAL(1),1,PVAL(IX),1)
          ETA = (PVAL(1) - PVAL(2))/ PVAL(3)
          NUCCHA = NINT(CHARGE(IATOM))
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
     &           SPIN,'Quadrupole moment: ',QMOM
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
        WRITE(6,'(2X,A)') 
     & 'No active atoms with quadrupole moment.'
      ENDIF
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(/A,F12.8,A)') '  1 a.u. = ',CONV,' MHz'
      RETURN
 9000 FORMAT (/2X,A9,A4,1X,A,I3,3X,A,I3,3X,A6,F4.1,3X,A19,F12.7,/)
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
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      DIMENSION EXPVAL(*), TEMP(3,3,MXCENT)
      LOGICAL DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
      CALL DZERO(EFGEL,9*MXCENT)
      CALL DZERO(TEMP,9*MXCENT)
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 100
         DO 200 ICOOR1 = 1, 3
            DO 210 ICOOR2 = ICOOR1, 3
               ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
               DO 300 IREPC = 0, MAXREP
                  IF (IBTAND(ISTBNU(IATOM),
     &                IBTXOR(IREPC,ISYMIJ)).EQ.0) THEN
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
        IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          DO 220 I = 1, 3
            DO 230 J = I, 3
              EFGEL(I,J,NATOM) = TEMP(I,J,IATOM)*
     &                           PT(IBTAND(IBTXOR(ISYMAX(I,1),
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL CORR
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      CALL MEMGET('LOGI',KATOM,NUCIND  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAXIS,9*MXCENT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,9*MXCENT,WORK,KFREE,LFREE)
C
C     (-30: accept .SELECT specification)
      CALL SETATM(WORK(KATOM),NATOM,-30)
      CALL PRP_EFG1(EXPVAL,EXP_CORR,WORK(KATOM),WORK(KAXIS),WORK(KTEMP),
     &              CORR,IPRPRP)
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
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),AXIS(3,3,MXCENT),TEMP(*)
      DIMENSION TMP(3),ITMP(3),PVAL(9)
      LOGICAL   DOATOM(*), CORR
      CHARACTER FMTEFG*80, FMTEFGA*6, FMTEFGB*6, FMTEFGC*6
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
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
      use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
C     conversion factor alpha/a.u. > 1.0D30*alpha/m**3
      PARAMETER (DFAC = XTANG**3)
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
      CALL AROUND('Electric dipole polarizability')
      WRITE(LUPRI,*)
C
C
C     print header
C     ============
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A)') '@'
      WRITE(LUPRI,'(A)') '@   Elements of the electric dipole'
      WRITE(LUPRI,'(A)') '@   polarizability tensor'
      WRITE(LUPRI,'(A)') '@'
      WRITE(LUPRI,*)
C
      DO IFREQ = 1,NBFREQ
        CALL DZERO(TMPMAT,9)
C
C
C       print frequency
C       ===============
C
        IF(.NOT.(NBFREQ.EQ.1)) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,'(A3,A10,F15.8,A5)') '@  ',
     &                         'frequency', BFREQ(IFREQ), 'a.u.'
          WRITE(LUPRI,'(A)') '@'
          WRITE(LUPRI,*)
        ENDIF
C
C
C       fill temporary matrix
C       =====================
C
        DO I = 1,3
          DO J = 1,3
            KI = IPEPOLAR(I,1)
            KJ = IPEPOLAR(J,1)
C
C           e-e
            TMPMAT(I,J) = -D1*ABLRF(IFREQ,KI,KJ,2)
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
        WRITE(LUPRI,'(A)') '@'
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
C
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL   TRLESS
C
#include "memint.h"
C
C     Allocate memory
C
      CALL MEMGET('LOGI',KATOM    ,NUCIND  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP    ,15*NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP_CORR,15*NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KNUC     ,15*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEL      ,15*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCORR    ,15*NUCDEP,WORK,KFREE,LFREE)
C
      IF (TRLESS) THEN
         CALL SETATM(WORK(KATOM),NATOM,52)
      ELSE
         CALL SETATM(WORK(KATOM),NATOM,53)
      END IF
      CALL PRP_EFT1(EXPVAL,EXP_CORR,WORK(KATOM),WORK(KTEMP),
     &              WORK(KTMP_CORR),WORK(KEL),WORK(KNUC),WORK(KCORR),
     &              TRLESS,CORR)
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
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
      PARAMETER (D0  =  0.0D00, D1 = 1.0D00, D3 = 3.0D00,
     &           D4  =  4.0D00, D7 = 7.0D00,D27=27.0D00,
     &           D35 = 35.0D00, DTHR = 1.0D-10)
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
#include "chrxyz.h"
#include "symmet.h"
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
#include "ibtfun.h"
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
            ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
            ISYMIJ = IBTXOR(ISYMIJ          ,ISYMAX(ICOOR3,1))
            ISYMIJ = IBTXOR(ISYMIJ          ,ISYMAX(ICOOR4,1))
            IOFF = 0
            DO 300 IREPC = 0, MAXREP
               IF (IBTAND(ISTBNU(IATOM),IBTXOR(IREPC,ISYMIJ)).EQ.0) THEN
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
            IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
               NATOM = NATOM + 1
               IELEM = 0
               DO ICOOR1 = 1, 3
               DO ICOOR2 = ICOOR1, 3
               DO ICOOR3 = ICOOR2, 3
               DO ICOOR4 = ICOOR3, 3
                  IELEM = IELEM + 1
                  ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
                  ISYMIJ = IBTXOR(ISYMIJ          ,ISYMAX(ICOOR3,1))
                  ISYMIJ = IBTXOR(ISYMIJ          ,ISYMAX(ICOOR4,1))
                  EFTEL(IELEM,NATOM) =
     &               TEMPEFT(IELEM,IATOM) * PT(IBTAND(ISYMIJ,IREPC))
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
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D3 = 3.0D0, D9 = 9.0D00, D15 = 15.0D0,
     &           D45 = 45.0D00, D90 = 90.0D00, D105 = 105.0D00)
      DIMENSION GEOM(3*8*MXCENT), EFTNUC(15,*)
#include "nuclei.h"
#include "symmet.h"
C
#include "ibtfun.h"
C
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
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
            IF (IBTAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               LL = 1
               DO 100 I2 = 1, NUCIND
                  IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
                  DO 110 ISYM = 0, MAXOPR
                     IF (IBTAND(ISYM,ISTBNU(I2)) .EQ. 0) THEN
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),WORK(*)
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      CALL MEMGET('LOGI',KATOM,NUCIND  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KMPV,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQV,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIS,NUCDEP,WORK,KFREE,LFREE)
C
C     (-62: accept .SELECT specification)
      CALL SETATM(WORK(KATOM),NATOM,-62)
      CALL PRP_PVC1(EXPVAL,WORK(KATOM),WORK(KMPV),WORK(KQV),
     &              WORK(KIS),WORK(KTEMP),IPRPRP)
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
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D00, D1 = 1.0D00, D2 = 2.0D00, D4 = 4.0D00)
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      DIMENSION EXPVAL(*),TEMP(*),DMPV(*),QW(*),
     &          IS(*)
      LOGICAL   DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
C
      CALL QENTER('PRP_PVC')
C
C     Gather expectation values into TEMP
C     (only over symmetry distinct centers)
C
      STHETAW2 = 0.2319D0
      DFAC = 2.22255D-14 / (D2 * SQRT(D2))
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
          ICHARGE = NINT(CHARGE(IATOM))
          IMASS = NINT(DISOTP(ICHARGE,1,'MASS'))
          NPROTON  = ICHARGE
          NNEUTRON = IMASS - ICHARGE
          IS(NATOM) = IMASS
          QW(NATOM) = (D1 - D4*STHETAW2) * NPROTON -
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
     &   'Nucleus','Isotop','integral    ','factor  ','contrib.    '
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
 9000 FORMAT(4X,A8,'  ',A6,'  ',A20,'  ',A10,'  ',A20)
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &          ABCNV(NBFREQ,NLRBPT),WORK(*)
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
      CALL MEMGET('LOGI',KATOM,NUCIND   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,IC*N2SCOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHES,IC*N2SCOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSPAR,5*NNDEP  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NSCOOR*2 ,WORK,KFREE,LFREE)
C
C     (-29: accept .SELECT specification)
      CALL SETATM(WORK(KATOM),NATOM,-29)
      CALL PRP_SPNSP1(EXPVAL,ABLRF,ABCNV,
     &                WORK(KATOM),WORK(KIATM),WORK(KRATM),
     &                WORK(KSMAT),WORK(KCHES),WORK(KSPAR),
     &                WORK(KIBUF),WORK,KFREE,LFREE)
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
      use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dummy.h"
      PARAMETER (D3 = 3.0D0, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *6.5796838999D15)
C
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
#include "symmet.h"
#include "dgroup.h"
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &          ABCNV(NBFREQ,NLRBPT),
     &          IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          SPNMAT(3*NUCDEP,3*NUCDEP),
     &          CHESS(3*NUCDEP,3*NUCDEP),
     &          SPNPAR(5,*),IBUF(3*NUCDEP,2),WORK(*)
C
#include "ibtfun.h"
      CALL QENTER('PRP_SPNSP1')
      CALL TITLER
     &    ('DIRAC  - INDIRECT NUCLEAR SPIN-SPIN-COUPLINGS','*',116)
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
        IZ   = NINT(CHARGE(INUC))
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
        CALL AROUND('Final reduced spin-spin-couplings K (SI)')
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
                FACTOR = AUTOHZ/3.85000804D-1
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
      use codata
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D3 = 3.0D0, D1 = 1.0D0)
#include "mxcent.h"
#include "nuclei.h"
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
      CALL RSJACO(3,3,3,SYM,PVAL,1,-2,0,AXES,WORK(KFREE),LFREE)
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
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "abainf.h"
#include "nuclei.h"
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
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0=0.D00, KCOL=6)
      DIMENSION SHESS(NROW,NROW), CHESS(NROW,NROW)
      DIMENSION WORK(LWORK)
#include "abainf.h"
#include "nuclei.h"
#include "symmet.h"
C
      INTEGER BEGIN
C
#include "memint.h"
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
      SUBROUTINE DEF_DSO(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define diamagnetic spin-orbit operator
C     Written by T.Saue - Feb 21 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     DSOINT: INTTYP = 12
C     (-12: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-12)
      CALL DEF_DSO1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_DSO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_DSO1 */
      SUBROUTINE DEF_DSO1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define diamagnetic spin-orbit operator
C     Written by T.Saue - Feb 21 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "chrxyz.h"
      LOGICAL DOATOM(NUCIND),SAME
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"
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
                           PNAME = 'DSO '//NAMN(IATOM1)
     &                      //CHRNOS(ISCOR1/10)//CHRNOS(MOD(ISCOR1,10))
     &                      //NAMN(IATOM2)
     &                      //CHRNOS(ISCOR2/10)//CHRNOS(MOD(ISCOR2,10))
                           PLABEL(1) = 'DSO '
     &                       //CHRNOS(ISCOR1/10)//CHRNOS(MOD(ISCOR1,10))
     &                       //CHRNOS(ISCOR2/10)//CHRNOS(MOD(ISCOR2,10))
                           CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                         INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),WORK(*)
      KFRSAV = KFREE
      NSCOOR  = NUCDEP*3
      N2SCOR  = NSCOOR*NSCOOR
      NNDEP   = NUCDEP*(NUCDEP+1)/2
      CALL MEMGET('LOGI',KATOM,NUCIND   ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,N2SCOR   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHES,N2SCOR   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSPAR,5*NNDEP  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,NSCOOR*2 ,WORK,KFREE,LFREE)
C
C     (-12: accept .SELECT specification)
      CALL SETATM(WORK(KATOM),NATOM,-12)
      CALL PRP_DSO1(EXPVAL,
     &              WORK(KATOM),WORK(KIATM),WORK(KRATM),
     &              WORK(KSMAT),WORK(KCHES),WORK(KSPAR),
     &              WORK(KIBUF),WORK,KFREE,LFREE)
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
      use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dummy.h"
      PARAMETER (PRMASS = 1.007276470D0, D3 = 3.0D0, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PRMASS*PRMASS)
     &                   *6.5796838999D15)
C
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "nuclei.h"
#include "symmet.h"
#include "dgroup.h"
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION EXPVAL(*),
     &          IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          SPNMAT(3*NUCDEP,3*NUCDEP),
     &          CHESS(3*NUCDEP,3*NUCDEP),
     &          SPNPAR(5,*),IBUF(3*NUCDEP,2),WORK(*)
C
#include "ibtfun.h"
C
      KFRSAV = KFREE
      CALL QENTER('PRP_DSO1')
      CALL TITLER
     &    ('DIRAC  - DSO CONTRIBUTIONS','*',116)
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
        IZ   = NINT(CHARGE(INUC))
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
        READ(PRPLBL(IPL)(5:6),'(I2)') IA
        READ(PRPLBL(IPL)(7:8),'(I2)') IB
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
      CALL AROUND('DSO contribution to reduced SSCC K (SI)')
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
              FACTOR = AUTOHZ/3.85000804D-1
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
#include "mxcent.h"
C
      DIMENSION WORK(LWORK)
      DIMENSION SHESS(NDIMS,NDIMS), CHESS(NDIMC,NDIMC)
C
#include "memint.h"
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
      SUBROUTINE DEF_SHIELD(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for NMR shieldings
C
C     Written by T. Enevoldsen Nov 7 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION WORK(LWORK)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     EF1INT: INTTYP = 29 (electric field at individual nuclei)
C     (-29: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-29)
      CALL DEF_SHIELD1(WORK(KDATOM),WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_SHIELD',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_shield1 */
      SUBROUTINE DEF_SHIELD1(DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operators for shieldings
C     DEF_SPNSP1 by T. Saue used as skeleton
C     Written by T. Enevoldsen - Oct. 1997
C     Last revision: 7. Nov. 1997
C                    22. Mar. 2000 (added ESR options, Anette Noerager + HJAaJ)
C
C   2002,2003/MI&HJAaJ - added TEC's changes for London atomic orbitals
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0,
     &          DMI2 = -0.5D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "gfac.h"
#include "pi.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "dcbexp.h"
#include "dcbnmr.h"
#include "chrxyz.h"
      LOGICAL DOATOM(NUCIND)
      CHARACTER PNAME*16, PLABEL(3)*8,LABEL*4,LABELC*7
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"
C
      CALL QENTER('DEF_SHIELD1')
C
      NNSTIN = 0
      NSCOOR = 3*NUCDEP
      CALL IZERO(IPSHIELD,3+MXCOOR)
C
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
CMI previously were here + signs (old TEC's convention)
CMI but the - sign is correct, so everything is rearranged
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
C       First define left-hand operators
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
     &'DEF_SHIELD1: Max 1000 coords in this implementation !')
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
              IREP1    = IBTXOR(IREPO,ISYMAX(IC2,1))
              IREP2    = IBTXOR(IREPO,ISYMAX(IC1,1))
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
     &           INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
              IF (SHIELD) THEN
                CALL OP1IND('XSHIELDA',INDOP1,LLRAPU,NLRAPT,
     &                      INDXPR,MAXLLR)
                IPSHIELD(ISCOR1) = INDOP1
              END IF
              IF (ESRANISO) THEN
                CALL OP1IND('ESRANISO',INDOP1,LESRP,NESRP,INDXPR,MAXESR)
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
     &                          INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
                    IF(IPRINT.GE.1) THEN
                       CALL WRIXPR(NNSTIN,INDXPR)
                    ENDIF
                    CALL OP1IND(LABEL,IPIND,LEXPP,NEXPP,
     &                          INDXPR,MAXEXP)
                    IPRNST(NNSTIN) = IPIND
                 END DO
              END IF
            ENDIF
  300     CONTINUE
        END IF
  200   CONTINUE
  100 CONTINUE
C
C     Define right hand operators
C     (PFAC(1:2) for Dirac mag. field, PFAC(3) for g_1 QED correction)
C

      IF (LONDON) THEN
        CALL DEF_D1HBLOND(1,WORK,LWORK,IPRINT)
      ELSE
C
C         Conventional atomic orbitals
C
      NCOMP = 2
      PFAC(1)  = -DMI2*CVAL
      PFAC(2)  = -DMI2*CVAL
      PFAC(3)  = -DMI2*G1FAC
C X component of B-field
      PNAME     = 'X magnetic field'
C     PNAME     = 'MAGFLDX'
      IPTYP     = IOP(1)
      PLABEL(1) = 'ZDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IF (SHIELD) CALL OP1IND('XSHIELB',IPSHIELD(1+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IF (ESRGTENS) THEN
         CALL OP1IND('X magfield',IPESR(1+MXCOOR,1),
     &      LESRP,NESRP,INDXPR,MAXESR)
         PNAME    =  'g_1 MAGFLDX'
         IPTYP    =  13
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,PFAC(3),PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('X g_1 magfield',IPESR(1+MXCOOR,2),
     &      LESRP,NESRP,INDXPR,MAXESR)
      END IF
      IF (DOESR) THEN
C        Sigma_X:
         PNAME    =  'x-spin'
         IPTYP    =  10
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,0.5D0,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IF (SHIELD) CALL OP1IND('YSHIELB',IPSHIELD(2+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IF (ESRGTENS) THEN
         CALL OP1IND('Y magfield',IPESR(2+MXCOOR,1),
     &      LESRP,NESRP,INDXPR,MAXESR)
         PNAME    =  'g_1 MAGFLDY'
         IPTYP    =  14
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,PFAC(3),PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('Y g_1 magfield',IPESR(2+MXCOOR,2),
     &      LESRP,NESRP,INDXPR,MAXESR)
      END IF
      IF (DOESR) THEN
C        Sigma_Y:
         PNAME    =  'y-spin'
         IPTYP    =  11
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,0.5D0,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('y spin',IPSIGMA(2),LESRP,NESRP,INDXPR,MAXESR)
      END IF
C Z component of B-field
      PNAME     = 'Z magnetic field'
C     PNAME     = 'MAGFLDZ'
      IPTYP     = IOP(3)
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IF (SHIELD) CALL OP1IND('ZSHIELB',IPSHIELD(3+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
      IF (ESRGTENS) THEN
         CALL OP1IND('Z magfield',IPESR(3+MXCOOR,1),
     &      LESRP,NESRP,INDXPR,MAXESR)
         PNAME    =  'g_1 MAGFLDZ'
         IPTYP    =  15
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,PFAC(3),PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('Z g_1 magfield',IPESR(3+MXCOOR,2),
     &      LESRP,NESRP,INDXPR,MAXESR)
      END IF
      IF (DOESR) THEN
C        Sigma_Z:
         PNAME    =  'z-spin'
         IPTYP    =  12
         PLABEL(1) = 'OVERLAP '
         CALL XPRIND(PNAME,IPTYP,1,0.5D0,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('z spin',IPSIGMA(3),LESRP,NESRP,INDXPR,MAXESR)
      END IF
#endif
      ENDIF
C
C     Print section for NMR shielding:
C
      IF (SHIELD) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(A)') ' NMR shieldings:'
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE (LUPRI,'(A,F9.4,A)')
     &   ' Shieldings printed for atoms with '//
     &   'abundance greater than', ABUND,' %'
C
C     A operators
C     ===========
C
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(A)') ' ** A operators for shielding **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NSCOOR
          INDOP1 = IPSHIELD(I)
          IF(INDOP1.GT.0) THEN
            INDXPR = LLRAPU(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO
C
C
C     B operators
C     ===========
C
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(A)') ' ** B operators for shielding or ESR **'
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
      IF (ESRANISO) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(A)') ' ** ESR anisotropy operators **'
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
     & ' ** NMR shieldings expectation values (London) **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NNSTIN
          INDOP1 = IPRNST(I)
          IF(INDOP1.GT.0) THEN
            INDXPR = LEXPP(INDOP1)
            CALL WRIXPR(INDOP1,INDXPR)
          ENDIF
        ENDDO
      ENDIF

      CALL QEXIT('DEF_SHIELD1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_OSCILL */
      SUBROUTINE DEF_OSCILL(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for oscillator strength calculations
C
C     Written by S. Knecht Sep 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "krciprop.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
C
#ifdef CODE_NOT_WORKING
      CHARACTER PNAME(6)*16, PLABEL(6)*8
      DIMENSION WORK(LWORK),PFAC(6),IOP(6)
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
     &              INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('KRCI PROPDIP',LUCIINDTRP,LPROP_KRCI,
     &               NPROP_KRCI,INDXPR,MXPROP_KRCI)
      END DO
#endif
      CHARACTER PNAME(3)*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
C
C     define operator type, factor, name and label
C
      NPCOMP = 1
      PFAC(1)  = D1
      PFAC(2)  = D1
      PFAC(3)  = D1
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
     &              INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('KRCI PROPDIP',LUCIINDTRP,LPROP_KRCI,
     &               NPROP_KRCI,INDXPR,MXPROP_KRCI)
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_OMEGAQUANT */
      SUBROUTINE DEF_OMEGAQUANT(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for determination of omega quantum numbers
C
C     Written by S. Knecht Sep 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "krciprop.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
C
      CHARACTER PNAME(2)*16, PLABEL(2)*8
      DIMENSION WORK(LWORK),PFAC(2),IOP(2)
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
     &              INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('KRCI PROPOME',LUCIINDTRP,LPROP_KRCI,
     &               NPROP_KRCI,INDXPR,MXPROP_KRCI)
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DEFINE_PVC_SHIELD */
      SUBROUTINE DEFINE_PVC_SHIELD(WORK,LWORK)
C=======================================================================
C     define operators for 
C     PV contribution to NMR shielding
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION WORK(LWORK)
#include "memint.h"
C
C
      CALL MEMGET('LOGI',KDOATOM,NUCIND,WORK,KFREE,LFREE)
      CALL SETATM(WORK(KDOATOM),NATOM,-62)
      CALL DEFINE_PVC_SHIELD_SUB(WORK(KDOATOM),WORK(KFREE),LFREE)
      CALL MEMREL('DEFINE_PVC_SHIELD',WORK,1,KFREE,KFREE,LFREE)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DEFINE_PVC_SPINSPIN */
      SUBROUTINE DEFINE_PVC_SPINSPIN(WORK,LWORK)
C=======================================================================
C     define operators for 
C     PV contribution to spin-spin coupling
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION WORK(LWORK)
#include "memint.h"
C
C
      CALL MEMGET('LOGI',KDOATOM,NUCIND,WORK,KFREE,LFREE)
      CALL SETATM(WORK(KDOATOM),NATOM,-29)
      CALL DEFINE_PVC_SPINSPIN_SUB(WORK(KDOATOM),WORK(KFREE),LFREE)
      CALL MEMREL('DEFINE_PVC_SPINSPIN',WORK,1,KFREE,KFREE,LFREE)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DEFINE_PVC_SHIELD_SUB */
      SUBROUTINE DEFINE_PVC_SHIELD_SUB(DOATOM,WORK,LWORK)
C=======================================================================
C     define operators for
C     PV contribution to NMR shielding
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "gfac.h"
#include "pi.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "dcbexp.h"
#include "dcbnmr.h"
#include "chrxyz.h"
C
      LOGICAL   DOATOM(NUCIND)
      CHARACTER PNAME*16,
     &          PLABEL(3)*8
      DIMENSION WORK(LWORK),
     &          PFAC(3),
     &          IOP(3)
C
#include "chrnos.h"
#include "ibtfun.h"
C
C
C     right hand
C     ==========
C
      CALL IZERO(IP_PVC_SHIELD,3+MXCOOR)
C
C     xavector
      IOP(1) = 5
C     yavector
      IOP(2) = 6
C     zavector
      IOP(3) = 7
C
      PFAC(1)   = 0.5D0*CVAL
      PFAC(2)   = 0.5D0*CVAL
      PFAC(3)   = 0.5D0*G1FAC
C
C     x - component
      PNAME     = 'MAGFLDX'
      IPTYP     = IOP(1)
      PLABEL(1) = 'ZDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,2,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XSHIELB',IP_PVC_SHIELD(1+MXCOOR),
     &            LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
C     y - component
      PNAME     = 'MAGFLDY'
      IPTYP     = IOP(2)
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,2,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YSHIELB',IP_PVC_SHIELD(2+MXCOOR),
     &            LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
C     z - component
      PNAME     = 'MAGFLDZ'
      IPTYP     = IOP(3)
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,2,PFAC,PLABEL,
     &            INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('ZSHIELB',IP_PVC_SHIELD(3+MXCOOR),
     &            LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
C
C     left hand
C     =========
C
      CALL IZERO(IP_PVC_SHIELD,MXCOOR)
      WRITE(LUPRI,'(A)') 'Parity-violating contribution'
      WRITE(LUPRI,'(A)') 'to NMR shielding'
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
      IF(.NOT.GAUNUC) THEN !point nucleus
        DO IXYZ = 1,3
          PFAC(IXYZ) = PFAC(IXYZ)*D3/(D4*GFAC*PI)
        ENDDO
      ENDIF
C
      DO IATOM = 1,NUCDEP
        IF(DOATOM(IATOM)) THEN
          DO IREP = 0,MAXREP
            IF(IBTAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
C
              IF(GAUNUC) THEN !gaussian model
                PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ELSE !point nucleus
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
     &                      INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
                CALL OP1IND(PNAME(1:6),INDEXA,LLRAPU,NLRAPT,
     &                      INDXPR,MAXLLR)
                ISCOOR = IPTCNT(3*(IATOM-1)+IXYZ,IREP,2)
                IF(ISCOOR.GT.0) THEN
                  IP_PVC_SHIELD(ISCOOR) = INDEXA
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
C     /* deck DEFINE_PVC_SPINSPIN_SUB */
      SUBROUTINE DEFINE_PVC_SPINSPIN_SUB(DOATOM,WORK,LWORK)
C=======================================================================
C     define operators for
C     PV contribution to spin-spin coupling
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "gfac.h"
#include "pi.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "dcbexp.h"
#include "dcbnmr.h"
#include "chrxyz.h"
C
      LOGICAL   DOATOM(NUCIND)
      CHARACTER PNAME*16,
     &          PLABEL(3)*8
      DIMENSION WORK(LWORK),
     &          PFAC(3),
     &          IOP(3)
C
#include "chrnos.h"
#include "ibtfun.h"
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
              IREP1 = IBTXOR(IREP,ISYMAX(IC2,1))
              IREP2 = IBTXOR(IREP,ISYMAX(IC1,1))
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
     &                    INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
      IF(.NOT.GAUNUC) THEN !point nucleus
        DO IXYZ = 1,3
          PFAC(IXYZ) = PFAC(IXYZ)*D3/(D4*GFAC*PI)
        ENDDO
      ENDIF
C
      DO IATOM = 1,NUCDEP
        IF(DOATOM(IATOM)) THEN
          DO IREP = 0,MAXREP
            IF(IBTAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
C
              IF(GAUNUC) THEN !gaussian model
                PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                           //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                           //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              ELSE !point nucleus
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
     &                      INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
C  /* Deck DEF_NSTDIA */
      SUBROUTINE DEF_NSTDIA(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define diamagnetic contribution to shieldings
C     (based on the Sternheimer approximation, see
C      Aucar et al., JCP 110 (1999) 6208) )
C     Written by T. Enevoldsen Oct. 1997
C    MI&HJJ 2003,2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     NSTDIACGO: INTTYP = 38
C     (-38: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-38)
      CALL DEF_NSTDIA1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_NSTDIA',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_NSTDIA1 */
      SUBROUTINE DEF_NSTDIA1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define diamagnetic contribution to shieldings
C     Written by T. Enevoldsen Oct. 1997
C
C  MI&HJJ/2002,2003 - extended for LAO 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "chrxyz.h"
      LOGICAL DOATOM(NUCIND),SAME
      CHARACTER PNAME*16, PLABEL(3)*8, LABEL*4
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"

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
     &' (Aucar et al, JCP 110 (1999) 6208)'
       IF (LONDON) THEN
        WRITE(LUPRI,'(2X,A)')
     &'- for London atomic orbitals (invariant gauge origin)'
       ELSE 
        WRITE(LUPRI,'(2X,A)')
     &'- for conventional atomic orbitals (common gauge origin)'
       ENDIF
      CALL PRSYMB(LUPRI,'=',75,0)
C
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
C        Atoms
C
         DO 200 IATOM1 = 1, NUCIND
            IF (DOATOM(IATOM1)) THEN
               MXATM2 = NUCIND
               IF (TRIANG) MXATM2 = IATOM1
C
C              Cartesian directions
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
     &                       //CHRNOS(ISECND)//CHRNOS(ITHIRD)
     &                       //CHRXYZ(ICOOR2)
                        PLABEL(1) = CHRNOS(IFIRST)//CHRNOS(ISECND)//
     &                       CHRNOS(ITHIRD)//LABEL//CHRXYZ(ICOOR2)
                        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                       INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
                        IF(IPRINT.GE.0) THEN
                           CALL WRIXPR(NNSTIN,INDXPR)
                        ENDIF
                        CALL OP1IND(LABEL,IPIND,LEXPP,NEXPP,
     &                       INDXPR,MAXEXP)
                        IPNSTDIA(NNSTIN) = IPIND
 600                 CONTINUE
                  END IF
 500           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE

      CALL QEXIT('DEF_NSTDIA1')   
      RETURN
      END

#ifdef MOD_MAGNETIZ
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_msuscp */
      SUBROUTINE DEF_MSUSCP(WORK,LWORK)
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
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER(DI2 = 0.5D0, DI4=0.25D0, DM1 = -1.0D0,
     &           D2 = 2.0D0, DM2 = -2.0D0, D1 = 1.0D0)
#include "dcbgen.h"

#include "dcbxlr.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbsusc.h"
#include "dcbham.h"
      CHARACTER PNAME*16, PLABEL(3)*8, LABELC*7
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)

      CALL QENTER('DEF_MSUSCP')
#include "memint.h"

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
    
      IF (LONDON) THEN ! use London AO's
C ... always calculate expectation value terms 
       DOEXP = .TRUE. 
C ... define operators for response and exp.value parts
       CALL DEF_D1HBLOND(2,WORK,LWORK,LPRINT)  
      ELSE ! Conventional AO

C      ... XAVECTOR
       IOP(1) = 5
C      ... YAVECTOR
       IOP(2) = 6
C      ... ZAVECTOR
       IOP(3) = 7
C
C     Define  operators for conventional orbitals
       NPCOMP   = 2
       PFAC(1)  = DI2*CVAL
       PFAC(2)  = DI2*CVAL
C   ....X componete of B-field
C      PNAME     = 'MAGFLDX'
       PNAME     = 'X magnetic field'
       IPTYP     = IOP(1)
       PLABEL(1) = 'ZDIPLEN'
       PLABEL(2) = 'YDIPLEN'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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

C     IF (XLR_SKIPEP.AND.(.NOT.LONDON)) THEN
      IF (MSUSCDIA.AND.(.NOT.LONDON)) THEN

       DOEXP = .TRUE. ! Calculate the expectation value !
C
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XXDSUSC',IPDSUSC(1),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-XY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XYSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XYDSUSC',IPDSUSC(2),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-XZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XZDSUSC',IPDSUSC(3),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-YY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YYSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YYDSUSC',IPDSUSC(4),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-YZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YZDSUSC',IPDSUSC(5),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-ZZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'ZZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('ZZDSUSC',IPDSUSC(6),LEXPP,NEXPP,INDXPR,MAXEXP)
C
C      Print section for CGO diamag.susceptibilities
C    ==================================================
C
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,'(1X,A)')
     & 'Diamagnetic susceptibility expectation value terms (no LAO):'
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
C    ... DO not calculate exp. value - diamag.term 
C        which is fully  involved in paramag. term when using URKBAL !

        DOEXP = .FALSE. 
        WRITE(LUPRI,'(1X,A)')
     & 'Diamagnetic susceptibility expectation value term'//
     & ' is  not calculated (MSUSCDIA=.false.)!'

C ... the case of London orbitals - do always calculate the exp.value !
      ELSE IF (LONDON) THEN
C  ... Here take care of one-electron terms 
C ... first 3 terms from the <r r~ hd> - RM2H1+RM2H2+RM2H3
C  ... 6 terms  Qmn.<r r~ beta' c^2>.Qmn - 6 diagonal operators
C     RM2H1 terms, IPDSUSC(1-6)

        PNAME     = 'LAO-XXRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS/dB2XX'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XXRM2H1',IPDSUSC(1),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS/dB2XY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XYRM2H1',IPDSUSC(2),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS/dB2XZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XZRM2H1',IPDSUSC(3),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS/dB2YY'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('YYRM2H1',IPDSUSC(4),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS/dB2YZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('YZRM2H1',IPDSUSC(5),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRM2H1'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = DM2*CVAL*CVAL
        PLABEL(1) = 'dS/dB2ZZ'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('ZZRM2H1',IPDSUSC(6),LEXPP,NEXPP,INDXPR,MAXEXP)

C ... RM2H2 terms, IPDSUSC(7-12)
        PNAME     = 'LAO-XXRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XXRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XXRM2H2',IPDSUSC(7),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XYRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XYRM2H2',IPDSUSC(8),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XZRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XZRM2H2',IPDSUSC(9),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YYRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('YYRM2H2',IPDSUSC(10),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YZRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('YZRM2H2',IPDSUSC(11),LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRM2H2'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'ZZRM2H2'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('ZZRM2H2',IPDSUSC(12),LEXPP,NEXPP,INDXPR,MAXEXP)

C -----------------------------------
C ... RM2H3 terms, IPDSUSC(13-18)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('ZZRM2H3',IPDSUSC(18),LEXPP,NEXPP,INDXPR,MAXEXP)

C ... RDSUSLL terms,  IPDSUSC(19-27)
C ... off-diagonal terms need to be symmetrized !
        IF (RDSUSLLMOD) THEN
           FACLL = -2.0D0
        ELSE
!radovan:  why is this factor different?
!          in which situation is this needed?
!miro: this is to be carefully examined
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('ZZRDSUL',IPDSUSC(27),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

C-----------------------------------------------------------------
C ... RDSUSNL terms,  pointers in IPDSUSC(28-33)
C   are TRUE diamagnetic terms according to the Sternheim  
C-----------------------------------------------------------------
C      IF (XLR_SKIPEP) THEN
C ... do include them 1.RKB-by default,  2.UKB - when specified (.MSUSCD)
       IF (MSUSCDIA) THEN
        PNAME     = 'LAO-XXRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XXDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XXDSUSN',IPDSUSC(28),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XYRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XYDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XYDSUSN',IPDSUSC(29),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-XZRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'XZDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('XZDSUSN',IPDSUSC(30),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YYRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YYDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('YYDSUSN',IPDSUSC(31),
     &                LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-YZRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'YZDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('YZDSUSN',IPDSUSC(32),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

        PNAME     = 'LAO-ZZRDSUSNL'
        IPTYP     = 1
        NPCOMP    = 1
        PFAC(1)   = D1
        PLABEL(1) = 'ZZDSUSNL'
        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
        CALL OP1IND('ZZDSUSN',IPDSUSC(33),
     &               LEXPP,NEXPP,INDXPR,MAXEXP)

C ... we need also DSUSCGM DALTON analogues, pointers IPDSUSC(34-39) 

      PNAME     = 'DSUSC_CGO-XX'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XXSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XXDSUSC',IPDSUSC(34),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-XY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XYSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XYDSUSC',IPDSUSC(35),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-XZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XZDSUSC',IPDSUSC(36),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-YY'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YYSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YYDSUSC',IPDSUSC(37),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-YZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YZDSUSC',IPDSUSC(38),LEXPP,NEXPP,INDXPR,MAXEXP)
C
      PNAME     = 'DSUSC_CGO-ZZ'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'ZZSUSCGO'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('ZZDSUSC',IPDSUSC(39),LEXPP,NEXPP,INDXPR,MAXEXP)

       ENDIF
C
C      Print section for LAO exp.values susceptibilities
C    =======================================================
C
       IF(IPRINT.GE.1) THEN
         WRITE(LUPRI,'(2X,A)') 'DEF_MSUSCP'//
     &': Expec.value terms term for the magnet. susc with LAO:'
C        IF (XLR_SKIPEP) THEN 
         IF (MSUSCDIA) THEN 
           NOP_DSUSC=39
         ELSE 
           NOP_DSUSC=27
         ENDIF
         DO I = 1, NOP_DSUSC
           INDXPR = LEXPP(IPDSUSC(I))
           CALL PRSYMB(LUPRI,'.',75,0)
           WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Component no.',I,':',PRPNAM(INDXPR)
           CALL PRSYMB(LUPRI,'.',75,0)
           CALL WRIXPR(I,INDXPR)
         ENDDO
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
      use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "dcbxpr.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbsusc.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D0=0.0D0,
     &           SI = (XTANG*XTANG*ECHARGE*ECHARGE/EMASS)*1.0D-21,
     &           CGS = SI*XFMOL*1.0D6, THRESHM = 1.0D-5)

      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &     ABCNV(NBFREQ,NLRBPT),WORK(*)
      LOGICAL D12, D23, MAGOUT1,MAGOUT2,MAGOUT3
      DIMENSION WRK(3), IWRK(3), AXES(3,3), PVAL(6),PVAL1(6),
     &         PVAL2(6),DIASUSC(3,3),PARSUS(3,3)
      DIMENSION SUSREL_EP(3,3,2),SUSDIA_EX(3,3), EXPVAL_NS(9)
      CHARACTER*6 ABC(3)
      CHARACTER*1 JCA,JCB
      INTEGER ICONVXYZ
      DATA ABC /'khi_aa','khi_bb','khi_cc'/

      CALL QENTER('PRP_MSUSCP')

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
     &   "PRP_MSUSCP: NEXPP is exceeded!"
         IF (XLR_SKIPEP) THEN
C  -1 is for the compatibility with Dalton
          SUSDIA(IPTAX(J,2),IPTAX(I,2)) = -D1*EXPVAL(IJX)
          SUSDIA(IPTAX(I,2),IPTAX(J,2)) = -D1*EXPVAL(IJX)
         ELSE                
          SUSDIA_EX(IPTAX(J,2),IPTAX(I,2)) = -D1*EXPVAL(IJX)
          SUSDIA_EX(IPTAX(I,2),IPTAX(J,2)) = -D1*EXPVAL(IJX)
         ENDIF
        ENDDO
        ENDDO

        IF (XLR_SKIPEP.AND.IPRPRP.GE.1) THEN
         CALL HEADER('Expectation values part of susceptibilities',-1)
         CALL OUTPUT(SUSDIA,1,3,1,3,3,3,1,LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF 

       ENDIF 

      ELSE IF (LONDON) THEN 
C ... accumulate one-electron expectation values into SUSDIA
C.......................................................................................
C   IPDSUSC(1-6)...pqRM2H1, IPDSUSC(7-12)...pqRM2H2, IPDSUSC(13-18)...pqRM3H3, 
C   IPDSUSC(19-27)...pqRDUSLL, IPDSUSC(28-33)...pqRDSUSNL, IPDSUSC(34-39)...pqSUSCGO
C
        CALL DZERO(SUSDIA,9)
        IJ = 0
        NEXPVALTERMS=3 ! RM2H1 + RM2H2 + RM2H3
        DO K = 1, NEXPVALTERMS
        DO I = 1, 3
        DO J = I, 3 
          IJ = IJ + 1
          IJX = IPDSUSC(IJ) 
          SUSDIA(IPTAX(J,2),IPTAX(I,2)) = EXPVAL(IJX) +
     &   SUSDIA(IPTAX(J,2),IPTAX(I,2))  
          SUSDIA(IPTAX(I,2),IPTAX(J,2)) = 
     &   SUSDIA(IPTAX(J,2),IPTAX(I,2))
        ENDDO
        ENDDO
        ENDDO

C -------------------------------------------------------
C  ...  RDSUSLL - 9 terms which have to be symmetrized 
C -------------------------------------------------------
C
C      ...  EXPVAL(IPDSUSC(19..27))
C
C              19  20  21
C              22  23  24
C              25  26  27
C

C       ...     symmetrize off-diagonal terms ...
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

C       ... control print-out of symmetrized and nosymmetrized RDSUSLL data as well
        IF (IPRPRP.GE.2) THEN
         WRITE(LUPRI,'(/,10X,A)')
     &   'RDSUSLL contributions to the expect.values:'
         WRITE(LUPRI,'(A)')
     &   '-----------------------------------------------------------'//
     &   '------------'
         WRITE(LUPRI,'(4X,A)')
     &   ' Operator             '//
     &   '  symmetrized       diff. from unsymmetrized'
         IIC = 1
         DO II = 19, 27
           INDXPR = LEXPP(IPDSUSC(II))
           WRITE(LUPRI,'(2X,A16,1X,F24.12,4X,D12.5)')
     &     PRPNAM(INDXPR),EXPVAL(IPDSUSC(II)),
     &     ( EXPVAL(IPDSUSC(II))-EXPVAL_NS(IIC) )
           IIC = IIC + 1
         ENDDO
        ENDIF

        DO I = 1, 3
        DO J = 1, 3 
         IJ = IJ + 1
         IJX = IPDSUSC(IJ) 
         SUSDIA(IPTAX(I,2),IPTAX(J,2)) = EXPVAL(IJX) +
     &   SUSDIA(IPTAX(I,2),IPTAX(J,2)) 
        ENDDO
        ENDDO

C ----------------------------------------------------------------------
C         LAST 6 true diamagnetic terms (RDSUSNL) 
C ----------------------------------------------------------------------
        IF (MSUSCDIA) THEN  
         DO I = 1, 3
         DO J = I, 3 
          IJ = IJ + 1
         IJX = IPDSUSC(IJ) 
          SUSDIA(IPTAX(J,2),IPTAX(I,2)) = EXPVAL(IJX) +
     &     SUSDIA(IPTAX(J,2),IPTAX(I,2)) 
          SUSDIA(IPTAX(I,2),IPTAX(J,2)) =  
     &     SUSDIA(IPTAX(J,2),IPTAX(I,2))
C         ... save them only when e-p rotations are active!
          IF (.NOT.XLR_SKIPEP) THEN
            SUSDIA_EX(IPTAX(J,2),IPTAX(I,2)) = EXPVAL(IJX)
            SUSDIA_EX(IPTAX(I,2),IPTAX(J,2)) = EXPVAL(IJX)
          ENDIF
         ENDDO
         ENDDO

C ... This is needed for compatibility with Dalton 
         DO I = 1, 3
         DO J = I, 3 
          IJ = IJ + 1
         IJX = IPDSUSC(IJ) 
          SUSCOM(IPTAX(J,2),IPTAX(I,2)) = -EXPVAL(IJX)
          SUSCOM(IPTAX(I,2),IPTAX(J,2)) = -EXPVAL(IJX)
         ENDDO
         ENDDO

        ENDIF

C ... Write output - always print out partial results
       IF (IPRPRP.GE.0) THEN
        CALL HEADER('One-electron expectation values',-1)
        CALL POLPRI(SUSDIA,'   ',-2)
        AVERAG = (SUSDIA(1,1) + SUSDIA(2,2) + SUSDIA(3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG

        CALL HEADER('Two-electron expectation values',-1)
        CALL POLPRI(SUS2EL,'   ',-2)
        AVERAG = (SUS2EL(1,1) + SUS2EL(2,2) + SUS2EL(3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG

        CALL HEADER('Highest-order reorthonormalization',-1)
        CALL POLPRI(SUSFS,'   ',-2)
        AVERAG = (SUSFS (1,1) + SUSFS (2,2) + SUSFS (3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
C
        CALL HEADER('Lowest-order reorthonormalization',-1)
        CALL POLPRI(SUSFSY,'   ',-2)
        AVERAG = (SUSFSY(1,1) + SUSFSY(2,2) + SUSFSY(3,3))/D3
        WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
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
     &    'PRP_MSUSCP: operator..',PRPNAM(LLRAPU(IA))
           write(lupri,*) 'extracted JCA,JA=',JCA,JA
           CALL QUIT('PRP_MSUSCP: wrong extracted number !')
           ENDIF
          ENDIF

          KA = IPMSUSC(1,IA)
          SUSREL(IPTAX(IB,2),IPTAX(IA,2)) = ABLRF(IFREQ,KA,KB,1)

          IF (.NOT.XLR_SKIPEP) THEN
C   .... save e-e "paramagnetic" and e-p "diamagnetic"
C        contributions from the linear response function
           SUSREL_EP(IPTAX(IB,2),IPTAX(IA,2),1) = ABLRF(IFREQ,KA,KB,2)
           SUSREL_EP(IPTAX(IB,2),IPTAX(IA,2),2) = ABLRF(IFREQ,KA,KB,3)
          ENDIF
         ENDDO
       ENDDO

C    .... DALTON type of output
         IF (LONDON) THEN  
          CALL HEADER('Relaxation',-1)
          CALL POLPRI(SUSREL,'   ',-2)
         ELSE
          CALL HEADER('Relaxation part of susceptibilities',-1)
          CALL OUTPUT(SUSREL,1,3,1,3,3,3,1,LUPRI)
         ENDIF
         AVERAG = (SUSREL(1,1) + SUSREL(2,2) + SUSREL(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
         IF (.NOT.XLR_SKIPEP) THEN
C ....write out e-e and e-p contributions...
          CALL HEADER('e-p "diamagnetic" and e-e "paramagnetic"'//
     &      ' contributions (au):',1)
          DO IROW = 1, 3
           WRITE(LUPRI,'(2X,I1,1X,3F12.8,3X,3F12.8)')
     &       IROW,(SUSREL_EP(IROW,ICOL,2),ICOL=1,3),
     &            (SUSREL_EP(IROW,ICOL,1),ICOL=1,3)
          ENDDO
          AVERAG_EE =
     &      (SUSREL_EP(1,1,1)+SUSREL_EP(2,2,1)+SUSREL_EP(3,3,1))/D3
          AVERAG_EP =
     &      (SUSREL_EP(1,1,2)+SUSREL_EP(2,2,2)+SUSREL_EP(3,3,2))/D3
          WRITE(LUPRI,'(/,1X,A,F12.6,4X,A,F12.6)')
     &    'e-p "diamagnetic" contribution:',AVERAG_EP,
     &    'e-e "paramagnetic":',AVERAG_EE
C        ... write out MSUSCDIA diamagnetic terms for comparison
          IF (MSUSCDIA) THEN
           CALL HEADER('.MSUSCD diamagnetic term for comparison '//
     &                 'with e-p "diamagnetic" term (au):',1)
           CALL OUTPUT(SUSDIA_EX,1,3,1,3,3,3,1,LUPRI)
           AVERAG_EX = 
     &      (SUSDIA_EX(1,1) + SUSDIA_EX(2,2) + SUSDIA_EX(3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F12.6)')
     &      '.MSUSCDIA average value:',AVERAG_EX
          ENDIF
 
          IF (.NOT.LONDON) THEN
C        ... copy e-e "paramagnetic" into SUSREL and e-p "diamagnetic" into SUSDIA
            CALL DCOPY(9,SUSREL_EP(1,1,1),1,SUSREL,1)
            CALL DCOPY(9,SUSREL_EP(1,1,2),1,SUSDIA,1)
          ELSE
C          ... case with London atomic orbitals...
C          ... copy the e-e "diamagnetic" contribution into SUSCOM
C           CALL DCOPY(9,SUSREL_EP(1,1,2),1,SUSCOM,1)
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
           SUSTOT(I,J)  = SUSTOT(I,J) + SUSREL(I,J)
     &                                + SUSDIA(I,J)
     &                                + SUS2EL(I,J)
     &                                + SUSFS (I,J)
     &                                + SUSFSY(I,J)
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
CMI   IF (D12 .AND. D23) THEN
      IF (MAGOUT1) THEN
         AVE = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         write(lupri,'(3X,A)') '...magnetizability output 1...'
         WRITE (LUPRI,'(3(/,1X,A,3(F12.4,A)))')
     &      ' Magnetizability:           ',    AVE,' (au), ',
     &        SI*AVE*1.0D31,' (SI), ', CGS*AVE,' (cgs)',
     &      ' Diamagnetic contribution:  ',DIAMAG,' (au), ',
     &        SI*DIAMAG*1.0D31,' (SI), ',CGS*DIAMAG,' (cgs)',
     &      ' Paramagnetic contribution: ',PARAMA,' (au), ',
     &        SI*PARAMA*1.0D31,' (SI), ',CGS*PARAMA,' (cgs)'
         WRITE (LUPRI,'(/,1X,A)') ' Magnetizability is spherical.'
CMI   ELSE IF (D12 .OR. D23) THEN
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
         CALL HEADER('Magnetizabilities           '
     &             //'         au          SI         cgs',1)
         write(lupri,'(3X,A)') '...magnetizability output 2...'
         WRITE (LUPRI,'(6(1X,A,3F12.4,/))')
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
CMI   ELSE
      ELSE IF (MAGOUT3) THEN
         CALL DSWAP(1,PVAL(2),1,PVAL(1),1)
         CALL DSWAP(3,AXES(1,2),1,AXES(1,1),1)
         AVE  = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         ANI1 = PVAL(3) - (PVAL(1) + PVAL(2))/D2
         ANI2 = PVAL(2) - (PVAL(1) + PVAL(3))/D2
         CALL HEADER('Magnetizabilities           '
     &             //'         au          SI         cgs',1)
         write(lupri,'(3X,A)') '...magnetizability output 3...'
         WRITE (LUPRI,'(5(1X,A,3F12.4,/))')
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
C      call quit('PRP_MSUSCP: blind branch of the output !!!')
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
#endif

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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &     ABCNV(NBFREQ,NLRBPT),WORK(*)
      KFRSAV = KFREE
      NCMAT =  NUCDEP*9
      CALL MEMGET('REAL',KTMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATLO,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT,2*NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KANIS,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KASYM,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSKEW,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSPAN,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KATOM,NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMPMAT,2*NCMAT,WORK,KFREE,LFREE)
C
C     (-29: accept .SELECT specification)
      CALL SETATM(WORK(KATOM),NATOM,-29)
      CALL PRP_SHIELD1(EXPVAL,ABLRF,WORK(KTMPMAT),WORK(KTMAT),
     &     WORK(KDMATLO),
     &     WORK(KDMAT),WORK(KPMAT),WORK(KTISO),WORK(KDISO),WORK(KPISO),
     &     WORK(KANIS),WORK(KASYM),WORK(KSKEW),WORK(KSPAN),WORK(KIBUF),
     &     WORK(KATOM),WORK(KCSTRA),WORK(KSCTRA),WORK,KFREE,LFREE)
      CALL MEMREL('PRP_SHIELD',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_SHIELD1 */
      SUBROUTINE PRP_SHIELD1(EXPVAL,ABLRF,TMPMAT,TMAT,
     &      DMATLO,DMAT,PMAT,TISO
     &     ,DISO,PISO,ANIS,ASYM,SKEW,SPAN,IBUFA
     &     ,DOATOM,CSTRA,SCTRA,WORK,KFREE,LFREE)
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
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "dcbxlr.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, D2I = 0.50D0)
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &     TMPMAT(3,*),TMAT(3,3,NUCDEP), DMAT(3,3,NUCDEP),
     &     DMATLO(3,3,NUCDEP),
     &     PMAT(3,3,NUCDEP,2),TISO(NUCDEP),DISO(NUCDEP),PISO(NUCDEP),
     &     ANIS(NUCDEP),ASYM(NUCDEP),SKEW(NUCDEP),SPAN(NUCDEP),
     &     IBUFA(3*NUCDEP),IBUFB(3),CSTRA(NUCDEP*NUCDEP*9),
     &     SCTRA(NUCDEP*NUCDEP*9),WORK(*)
      INTEGER ICONVXYZ
      CHARACTER*1 CJB
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
      CALL QENTER('PRP_SHIELD1')
C 
      CALL TITLER('DIRAC  - NMR SHIELDINGS','*',116)
C
      KFRSAV = KFREE
      NSCOOR = 3*NUCDEP
      NCMAT = 9*NUCDEP
      CALL DZERO(TMPMAT,NCMAT)
      CALL DZERO(DMAT,NCMAT)
      CALL DZERO(DMATLO,NCMAT)
      CALL DZERO(TMAT,NCMAT)
      CALL DZERO(PMAT,2*NCMAT)

      IF (.NOT.XLR_SKIPEP.AND.NSTDIA) THEN
         IF (LONDON) THEN
          WRITE(LUPRI,'(2X,A)')
     &  'For LAO & e-p rotations the diamag.term is not (yet) defined !'
          WRITE(LUPRI,'(2X,A)')
     &  'NSTDIA set to .false. ! '
          NSTDIA=.FALSE.
         ENDIF
      ENDIF
C--------------------------------------------------------------------------
C     Extract Diamagnetic approximation if calculated
C     (used for total shielding if e-p rotations are not included)
C       -- this becomes the
C     A**2 diamagnetic operator in the non-relativistic limit
C     (ref: Aucar et al., JCP 110 (1999) 6208) /hjaaj Mar 2003
C
C ... extract diamag. term only for e-e rotations in lin.response function
C     IF(NSTDIA.AND.XLR_SKIPEP) THEN
      IF(NSTDIA) THEN
         DO I = 1,NNSTIN
            IOP = IPNSTDIA(I)
            IPL = IPRPLBL(1,LEXPP(IOP))
CMI       ... take care of the correct form of PRPLBL(IPL) name
            READ(PRPLBL(IPL)(1:3),'(I3)') IA
            IB = MOD(I+2,3) + 1
            TMPMAT(IPTAX(IB,2),IA) = EXPVAL(IOP)
         ENDDO
         IF (IPRPRP.GE.2) THEN
C         ...for what do we have calculated diamagnetic terms
          IF (XLR_SKIPEP) THEN
            CALL HEADER('Diamagnetic part of shielding'//
     &    ' in symmetry coordinates (ppm)',1)
          ELSE
            CALL HEADER('Diamagnetic terms for comparison'//
     &    ' with e-p "diamagnetic" in symmetry coordinates (ppm)',1)
          ENDIF
            CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
         ENDIF
C        ... DMAT contains extracted diamagnetic part (NSCO or NSNL)
         IF ( (LONDON.AND.XLR_SKIPEP).OR.(.NOT.LONDON) ) THEN
            CALL TRADIP(TMPMAT,DMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C        ENDIF
         ELSE
C        ... DMATLO contains extracted 'diamagnetic' term (NSNL)
          CALL TRADIP(TMPMAT,DMATLO,CSTRA,SCTRA,3*NUCDEP,2,2)
         ENDIF
      END IF

C ------------------------------------------------------
C  Get the London orbital contributions to the exp.value
C  ... these are RNST terms
C ------------------------------------------------------
      IF(LONDON) THEN
         CALL DZERO(TMPMAT,NCMAT)
         CALL DZERO(PMAT,NCMAT)
         DO I = 1,NNSTIN
            IOP = IPRNST(I)
            IPL = LEXPP(IOP)
CMI ... take care of the correct form of PRPLBL(IPL) name
            READ(PRPNAM(IPL)(10:11),'(I2)') IA
            IB = MOD(I+2,3) + 1
            TMPMAT(IPTAX(IB,2),IA) = -EXPVAL(IOP)
         ENDDO
         IF (IPRPRP.GE.2) THEN
            CALL HEADER('London part of shielding (exp.value)'//
     &           ' in symmetry coordinates (ppm)',1)
            CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
         ENDIF
C         ... extract London part of shielding into the PMAT
         CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C         ... add the London exp.value parts of shield. to to the DMAT
         CALL DAXPY(NCMAT,D1,PMAT,1,DMAT,1)
      END IF

C--------------------------------------------------------------------
C  Get the (expectation value) reorthonormalization terms, {T,h(mK)}
C--------------------------------------------------------------------
      IF(LONDON.AND.(.NOT.NOORTH)) THEN
         CALL MEMGET('REAL',KW,3*N2ORBXQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KCON,N2ORBXQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KDMAT,N2ORBXQ,WORK,KFREE,LFREE)
         CALL DZERO(TMPMAT,NCMAT)
         CALL DZERO(PMAT,NCMAT)
C  ...  get the expect. value {T,h(mK)} reorthonormalization term
!todo from Trond's email:
!check the use of GETCON in  PRP_SHIELD1(pamprp.F). Can it be completely removed ?
         CALL GETCON(WORK(KW),WORK(KCON),WORK(KDMAT),TMPMAT,
     &               WORK,KFREE,LFREE,IPRPRP)
CMI   ...  sign convention !!!
         CALL DSCAL(NCMAT,DM1,TMPMAT,1)
         IF (IPRPRP.GE.2) THEN
          CALL HEADER('Reorthonormalization part of the '// 
     &  'shielding (expect.value) in symmetry coordinates (ppm)',1)
          CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
         ENDIF
C         ... extract the reorthonormalization part into the PMAT
         CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C         ... add the reorthonormalization terms to the DMAT
         CALL DAXPY(NCMAT,D1,PMAT,1,DMAT,1)
      END IF
C
C     *** Extract shieldings ***
C
C  The NSCOOR nucl. magn. dipole operators as A operators in <<A;B>> :
C
      NA = 0
      NB = 0
      DO I = 1,NSCOOR
         IF(IPSHIELD(I).GT.0) THEN
            NA = NA + 1
            IBUFA(NA) = I
         ENDIF
      END DO
CMI/HJ internal consistency check
      IF (NA.LE.0) THEN
       WRITE(LUPRI,'(2X,A,I2)')
     & 'PRP_SHIELD1: NA=',NA
       CALL QUIT('PRP_SHIELD1:Wrong value of NA!')
      ENDIF

C    The 3 B field directions as B operators in <<A;B>> :
      DO I = 1, 3
         IF(IPSHIELD(I+MXCOOR).GT.0) THEN
            NB = NB + 1
            IBUFB(NB) = I
         ENDIF
      END DO
      IF (NB.LE.0) THEN
       WRITE(LUPRI,'(2X,A,I2)')
     & 'PRP_SHIELD1: NB=',NB
       CALL QUIT('PRP_SHIELD1:Wrong value of NB!')
      ENDIF
C
C     Fill shielding matrices
C
      IF (XLR_SKIPEE) THEN
         IF (XLR_SKIPEP) THEN
           IABLRF_EE = -999999999
           IABLRF_EP = -999999999
         ELSE
           IABLRF_EE = -999999999
           IABLRF_EP = 2
         END IF
      ELSE
         IF (XLR_SKIPEP) THEN
           IABLRF_EE = 2
           IABLRF_EP = -999999999
         ELSE
           IABLRF_EE = 2
           IABLRF_EP = 3
         END IF
      END IF
#if defined (XH4)
      WRITE(LUPRI,'(A)') '*** WARNING *** Using code for XH4 !'
#endif
      DO IFREQ = 1,NBFREQ
C
         CALL PRSYMB(LUPRI,'-',72,2)
         IF(.NOT.BFREQ(IFREQ).EQ.D0) THEN
            WRITE(LUPRI,'(A,F15.8,2X,A)')
     &           '*** Frequency    :',BFREQ(IFREQ),'a.u.'
         ENDIF

         CALL DZERO(TMPMAT,2*NCMAT)
         DO IB = 1,NB
CTROND           IBINDEX = LLRBPU(IPSHIELD(MXCOOR+IBUFB(IB)))
CTROND           IF (LONDON) THEN
CTROND             READ (PRPNAM(IBINDEX),'(4X,A1)') CJB
CTROND             JB = ICONVXYZ(CJB)
CTROND             IF (JB.NE.1.AND.JB.NE.2.AND.JB.NE.3) THEN
CTROND               write(lupri,*) '...operator name:',PRPNAM(IBINDEX)
CTROND               write(lupri,*) 'extracted CJB,JB:',CJB,JB
CTROND               CALL QUIT('PRP_SHIELD1: wrong index !')
CTROND             ENDIF
CTROND           ELSE
CMI  ...     first character, like 'X magnetic field'
CTROND
CTROND             WRITE(6,*) 'test..',PRPNAM(IBINDEX)
CTROND             READ (PRPNAM(IBINDEX),'(A1)') CJB 
CTROND            !write(lupri,*) 'PRPNAM(IBINDEX)=',PRPNAM(IBINDEX)
CTROND             JB = ICONVXYZ(CJB)
CTROND             IF (JB.NE.1.AND.JB.NE.2.AND.JB.NE.3) THEN
CTROND               write(lupri,*) '...operator name:',PRPNAM(IBINDEX)
CTROND               write(lupri,*) 'extracted CJB,JB:',CJB,JB
CTROND               CALL QUIT('PRP_SHIELD1: wrong index !')
CTROND             ENDIF
CTROND           ENDIF
            JB = IBUFB(IB)
            KB = IPSHIELD(JB+MXCOOR)
            DO IA = 1,NA
               JA = IBUFA(IA)
               KA = IPSHIELD(JA)
C...           First extract the e-e contribution
               IF (.NOT. XLR_SKIPEE) TMPMAT(IPTAX(JB,2),JA) =
     &            -ABLRF(IFREQ,KA,KB,IABLRF_EE)
C...           then extract the e-p "diamagnetic" contribution
               IF (.NOT. XLR_SKIPEP) TMPMAT(IPTAX(JB,2),JA+NCMAT) =
     &            -ABLRF(IFREQ,KA,KB,IABLRF_EP)
            ENDDO
         ENDDO
#if defined (XH4)
         TMPMAT(2,8)  = TMPMAT(1,4)
         TMPMAT(2,9)  = TMPMAT(1,6)
         TMPMAT(2,10) = TMPMAT(1,7)
         TMPMAT(2,11) = TMPMAT(1,5)
         TMPMAT(3,12) = TMPMAT(1,4)
         TMPMAT(3,13) = TMPMAT(1,6)
         TMPMAT(3,15) = TMPMAT(1,7)
         TMPMAT(3,14) = TMPMAT(1,5)
         TMPMAT(2,8+NCMAT)  = TMPMAT(1,4+NCMAT)
         TMPMAT(2,9+NCMAT)  = TMPMAT(1,6+NCMAT)
         TMPMAT(2,10+NCMAT) = TMPMAT(1,7+NCMAT)
         TMPMAT(2,11+NCMAT) = TMPMAT(1,5+NCMAT)
         TMPMAT(3,12+NCMAT) = TMPMAT(1,4+NCMAT)
         TMPMAT(3,13+NCMAT) = TMPMAT(1,6+NCMAT)
         TMPMAT(3,15+NCMAT) = TMPMAT(1,7+NCMAT)
         TMPMAT(3,14+NCMAT) = TMPMAT(1,5+NCMAT)
#endif
C
C     Transform to non-symmetry basis
C
         CALL DZERO(PMAT,2*NCMAT)
         CALL DZERO(TMAT,NCMAT)
         IF (.NOT. XLR_SKIPEE) THEN
            IF (IPRPRP.GE.2) THEN
               CALL HEADER('e-e response part of shielding'//
     &            ' in symmetry coordinates (ppm)',1)
               CALL FCPRI (TMPMAT,'SIGMAS',CSTRA,SCTRA)
            ENDIF
C           .... extract the e-e response part into PMAT(...,1)
            CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C           .... save e-e response part into TMAT
            CALL DAXPY(9*NUCDEP,D1,PMAT(1,1,1,1),1,TMAT,1)
         END IF
CMI/HJ  ... include diamagnetic(=Sternheimer) term if XLR_SKIPEP
C       ... include London exp.value correction to paramagnetic term if LONDON
         IF (XLR_SKIPEP) THEN
C       ... add the diamagnetic term(s),DMAT, to TMAT(contains paramagnetic terms)
            CALL DAXPY(9*NUCDEP,D1,DMAT,1,TMAT,1)
         ELSE
C       ... add (diamagnetic-NSNL)+London+reortho terms to the TMAT
            IF (LONDON) CALL DAXPY(9*NUCDEP,D1,DMAT,1,TMAT,1)
            IF (IPRPRP.GE.2) THEN
               CALL HEADER('e-p response part of shielding'//
     &            ' in symmetry coordinates (ppm)',1)
               CALL FCPRI (TMPMAT(1,1+NCMAT),'SIGMAS',CSTRA,SCTRA)
            ENDIF
C          ... extract the e-p response part into PMAT(...,2)
            CALL TRADIP(TMPMAT(1,1+NCMAT),PMAT(1,1,1,2),
     &                  CSTRA,SCTRA,3*NUCDEP,2,2)
C          ... add the e-p response PMAT(...,2) into TMAT
C           so the TMAT now contains total response
            CALL DAXPY(9*NUCDEP,D1,PMAT(1,1,1,2),1,TMAT,1)
         END IF
C
C     Analyze results
C
C ... for analyzing e-p "diamagnetic" contribution 
C  NoLONDON - DMAT contains pure diamagnetic term,NSCO
C  LONDON - DMAT contains London+reortho terms
         IF (LONDON.AND..NOT.XLR_SKIPEP) 
C    &    CALL DCOPY(NCMAT,DMATLO,1,DMAT,1)
     &    CALL DAXPY(NCMAT,D1,DMAT,1,DMATLO,1)

         IATOM = 0
         DO 100 I = 1, NUCIND
            DO 200 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(I)) THEN
                     CALL SHIANA(TMAT(1,1,IATOM),DMAT(1,1,IATOM),
     &                 DMATLO(1,1,IATOM),
     &                 PMAT(1,1,IATOM,1),PMAT(1,1,IATOM,2),
     &                 TISO(IATOM),DISO(IATOM),
     &                 PISO(IATOM),ANIS(IATOM),ASYM(IATOM),SKEW(IATOM),
     &                 SPAN(IATOM),NAMDEP(IATOM),NAMDPX(3*(IATOM-1)+1),
     &                 IPRPRP)
                  END IF
               END IF
 200        CONTINUE
 100     CONTINUE
C
         CALL AROUND('Summary of chemical shieldings')
         WRITE (LUPRI,'(A)') ' Definitions from J.Mason, Solid '//
     &        'state Nuc.Magn.Res. 2 (1993), 285'
         WRITE (LUPRI,'()')
         WRITE (LUPRI,'(A/A/A/A)')
     &        '@            isotropic shielding',
     &        '@        ----------------------------',
     &        '@atom          total          dia         para'//     
     &        '            skew'//
     &        '            span            anis            asym',
     &        '@-----------------------------------------------------'//
     &        '-----------------------'

         IATOM = 0
         DO 300 I = 1, NUCIND
            DO 400 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(I)) THEN
                      WRITE (LUPRI,'(A,A,7F14.4)') '@',
     &                 NAMDEP(IATOM),TISO(IATOM),DISO(IATOM),PISO(IATOM)
     &                 ,SKEW(IATOM),SPAN(IATOM),ANIS(IATOM),ASYM(IATOM)
                   END IF
               END IF
 400        CONTINUE
 300     CONTINUE
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
      SUBROUTINE SHIANA(TMAT,DMAT,DMATLO,PMATEE,PMATEP,
     &                  AVETOT,AVEDIA,AVEPAR,ANIS,ASYM,
     &                  SKEW,SPAN,NAME,NAMEX,IPRINT)
C***********************************************************************
C
C     Analyze routine for NMR shieldings
C     This routine is modelled after SHIANA in DALTON
C
C     Written by T.Enevoldsen Nov. 14 1997
C
CMI/Dec '04  DMATLO is reserved for future use
C***********************************************************************
      use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbxlr.h"
      PARAMETER (FACTOR = 1.D6*ALPHA2)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      CHARACTER NAME*6, NAMEX(3)*6, SGN*3
      INTEGER Z
      DIMENSION TMAT(3,3), DMAT(3,3), PMATEE(3,3), PMATEP(3,3),
     &          SYM(3,3),  ANTI(3,3), AXES(3,3), DMATLO(3,3),
     &          PVAL(6), WRK(3), IWRK(3), ABSP(3)
#include "chrxyz.h"
C
      CALL DSCAL(9,FACTOR,TMAT,1)
      CALL DSCAL(9,FACTOR,DMAT,1)
      CALL DSCAL(9,FACTOR,DMATLO,1)
      CALL DSCAL(9,FACTOR,PMATEE,1)
      CALL DSCAL(9,FACTOR,PMATEP,1)

C     .... total shielding
      AVETOT = (TMAT(1,1) + TMAT(2,2) + TMAT(3,3))/D3
      IF (XLR_SKIPEP) THEN
C        ... use NSTDIA expectation value approximation
         AVEDIA = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
      ELSE
C        ... use result of e-p response equations
C        AVEDIA = (PMATEP(1,1) + PMATEP(2,2) + PMATEP(3,3))/D3
        IF (.NOT.LONDON) THEN
         AVEDIA  = (PMATEP(1,1) + PMATEP(2,2) + PMATEP(3,3))/D3
         AVEDIA2 = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
        ELSE
         AVEDIA  = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
         AVEDIA1 = (PMATEP(1,1) + PMATEP(2,2) + PMATEP(3,3))/D3
         AVEDIA2 = (DMATLO(1,1) + DMATLO(2,2) + DMATLO(3,3))/D3
        ENDIF
      END IF
      IF (.NOT.LONDON) THEN
        AVEPAR = (PMATEE(1,1) + PMATEE(2,2) + PMATEE(3,3))/D3
      ELSE
        AVEPAR = (PMATEE(1,1) + PMATEE(2,2) + PMATEE(3,3))/D3
     &         + (PMATEP(1,1) + PMATEP(2,2) + PMATEP(3,3))/D3
      ENDIF
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
C This is done to avoid numerical insatbilities
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
      DO 400 ICOOR = 1, 3
         WRITE (LUPRI,'(2X,A6,3F15.8)')
     &          NAMEX(ICOOR), (TMAT(K,ICOOR),K=1,3)
  400 CONTINUE
C
C     IF (XLR_SKIPEP.OR.LONDON) THEN
      IF (XLR_SKIPEP) THEN
C   XLR_SKIPEP: Diamagnetic contribution (Sternheimer approx.)
         CALL HEADER('Diamagnetic and paramagnetic contributions (ppm):'
     &      ,1)
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
     &      ' contributions (ppm):',1)
         WRITE (LUPRI,'(15X,3(A,9X),4X,3(A,9X)/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
         DO ICOOR = 1, 3
            WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &         (PMATEP(K,ICOOR),K=1,3), (PMATEE(K,ICOOR),K=1,3)
         END DO
         IF (.NOT.LONDON) THEN
           AVEDIAX=AVEDIA
           AVEPARX=AVEPAR
         ELSE
           AVEDIAX=AVEDIA1
           AVEPARX=AVEPAR-AVEDIA1
         ENDIF
         WRITE (LUPRI,'(/A,F14.6,8X,A,F14.6)')
     &          ' e-p "diamagnetic" contribution:',AVEDIAX,
     &          ' e-e "paramagnetic": ',AVEPARX
         IF (NSTDIA) THEN
            IF (.NOT.LONDON) THEN
              CALL HEADER('.NSTDIA CGO diamagnetic term '//
     &        ' for comparison with e-p "diamagnetic" term (ppm):',1)
            ELSE
              CALL HEADER('.NSTDIA London diamagnetic term '//
     &        'for comparison  with e-p "diamagnetic" term (ppm):',1)
            ENDIF
            WRITE (LUPRI,'(15X,3(A,9X)/)')
     &       'Bx', 'By', 'Bz'
            DO ICOOR = 1, 3
             IF (.NOT.LONDON) THEN
               WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &           (DMAT(K,ICOOR),K=1,3)
             ELSE 
               WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &           (DMATLO(K,ICOOR),K=1,3)
             ENDIF
            END DO
            WRITE (LUPRI,'(/A,F14.6,8X,A,F14.6)')
     &          ' .NSTDIA diamagnetic contribution, for comparison:',
     &          AVEDIA2
         END IF
      END IF
C
      CALL HEADER
     &   ('Antisymmetric and traceless symmetric parts (ppm):',1)
      WRITE (LUPRI,'(15X,3(A,9X),4X,3(A,9X),/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
      DO 600 ICOOR = 1, 3
         WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &         (ANTI(K,ICOOR),K=1,3), (SYM(K,ICOOR),K=1,3)
  600 CONTINUE
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
  700 CONTINUE
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PRINTOUT_PVC_SHIELD */
      SUBROUTINE PRINTOUT_PVC_SHIELD(ABLRF,ABCNV,WORK,KFREE,LFREE)
C=======================================================================
C     PV contribution to NMR shielding
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION WORK(*)
C
C
      CALL MEMGET('LOGI',KDOATOM,NUCIND,WORK,KFREE,LFREE)
      CALL SETATM(WORK(KDOATOM),NATOM,-62)
      CALL PRINTOUT_PVC_SHIELD_SUB(ABLRF,ABCNV,WORK(KDOATOM),
     &                             WORK,KFREE,LFREE)
      CALL MEMREL('PRINTOUT_PVC_SHIELD',WORK,1,KFREE,KFREE,LFREE)
C
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "nuclei.h"
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
C     /* deck PRINTOUT_PVC_SHIELD_SUB */
      SUBROUTINE PRINTOUT_PVC_SHIELD_SUB(ABLRF,ABCNV,DOATOM,
     &                                   WORK,KFREE,LFREE)
C=======================================================================
C     PV contribution to NMR shielding
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
      use codata
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      PARAMETER ( GFERMI = 2.22254D-14 )
      PARAMETER ( STHETAW2 = 0.2319D0 )
C
      LOGICAL   DOATOM(NUCDEP)
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),
     &          RESPONSE_MATRIX(3,3),
     &          WORK(*)
C
#include "symmet.h"
#include "dgroup.h"
#include "ibtfun.h"
C
C
      FACTOR = (GFERMI*(D1-D4*STHETAW2)*1.0D6)/(CVEL*SQRT(D2))
C
      CALL AROUND('PV NMR shielding tensor (in ppm)')
      WRITE(LUPRI,*)
C
      DO IATOM = 1,NUCDEP
        IF(DOATOM(IATOM)) THEN
C
C
C         extract response
C         ================
C
          DO IXYZ = 1,3
            DO JXYZ = 1,3
              KA = IP_PVC_SHIELD((IATOM*3)-3+IXYZ)
              KB = IP_PVC_SHIELD(JXYZ+MXCOOR)
              RESPONSE_MATRIX(IXYZ,JXYZ) = ABLRF(1,KA,KB,2)
              IF(.NOT. XLR_SKIPEP) THEN
                RESPONSE_MATRIX(IXYZ,JXYZ) = RESPONSE_MATRIX(IXYZ,JXYZ)
     &                                     + ABLRF(1,KA,KB,3)
              ENDIF
            ENDDO
          ENDDO
C
C
C         extract nuclear g factor
C         ========================
C
          IZ = NINT(CHARGE(IATOM))
          ABUNDANCE = D0
          DO ISO = 1,5
            GVAL = DISOTP(IZ,ISO,'GVAL')
            IF(GVAL .NE. D0) GOTO 1
          ENDDO
 1        CONTINUE
C
          CALL DSCAL(9,FACTOR,RESPONSE_MATRIX,1) !scale by FACTOR
C
          WRITE(LUPRI,'(2A10)') 'nucleus', NAMDEP(IATOM)
          WRITE(LUPRI,'(A10,E20.10)') 'g-factor', GVAL
C
          RATIO = GVAL/(D2*XFMP)
          IF(RATIO .NE. D0) THEN
            BY_RATIO = D1/RATIO
            CALL DSCAL(9,BY_RATIO,RESPONSE_MATRIX,1) !divide by magnetogyric ratio 
          ELSE
            WRITE(LUPRI,'(A)') '  g-factor is zero'
            WRITE(LUPRI,'(A)') '  therefore no division'
            WRITE(LUPRI,'(A)') '  by magnetogyric ratio'
          ENDIF
C
          CALL PRINTOUT_TENSOR(RESPONSE_MATRIX)
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'-',46,0)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PRINTOUT_TENSOR */
      SUBROUTINE PRINTOUT_TENSOR(TENSOR)
C=======================================================================
C     print tensor
C-----------------------------------------------------------------------
C     radovan bast                           last revision: october 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
      DIMENSION TENSOR(3,3),
     &          SYM_TENSOR(3,3),
     &          PRINCIPAL_VAL(9),
     &          PRINCIPAL_AX(3,3),
     &          SOME(3),
     &          ISOME(3)

      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A22,2A20)') 'x', 'y', 'z'
      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A10,3E20.10)') 'x', (TENSOR(1,I),I=1,3)
      WRITE(LUPRI,'(A10,3E20.10)') 'y', (TENSOR(2,I),I=1,3)
      WRITE(LUPRI,'(A10,3E20.10)') 'z', (TENSOR(3,I),I=1,3)

      AVG_TRACE = D0
      DO I = 1,3
        AVG_TRACE = AVG_TRACE + TENSOR(I,I)
      ENDDO
      AVG_TRACE = AVG_TRACE/D3

      WRITE(LUPRI,'(A10,3A20,E20.10)')
     &                               'trace/3', ' ', ' ', ' ', AVG_TRACE

      DO I = 1,3
        DO J = 1,3
          SYM_TENSOR(I,J) = (TENSOR(I,J) + TENSOR(J,I))/D2
        ENDDO
      ENDDO
      CALL DUNIT(PRINCIPAL_AX,3)
      IJ = 1
      DO I = 1,3
        DO J = 1,I
          PRINCIPAL_VAL(IJ) = SYM_TENSOR(I,J)
          IJ = IJ + 1
        ENDDO
      ENDDO
      CALL JACO(PRINCIPAL_VAL,PRINCIPAL_AX,3,3,3,SOME,ISOME)

      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A30)') 'principal components'
      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(E30.10)') PRINCIPAL_VAL(1)
      WRITE(LUPRI,'(E30.10)') PRINCIPAL_VAL(3)
      WRITE(LUPRI,'(E30.10)') PRINCIPAL_VAL(6)
      WRITE(LUPRI,'(A)')

      AVG = PRINCIPAL_VAL(1) + PRINCIPAL_VAL(3) + PRINCIPAL_VAL(6)
      AVG = AVG/D3

      WRITE(LUPRI,'(A10,E20.10)') 'sum/3', AVG
      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A30)') 'principal axes'
      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A10,3E20.10)') 'x', (PRINCIPAL_AX(1,I),I=1,3)
      WRITE(LUPRI,'(A10,3E20.10)') 'y', (PRINCIPAL_AX(2,I),I=1,3)
      WRITE(LUPRI,'(A10,3E20.10)') 'z', (PRINCIPAL_AX(3,I),I=1,3)
      WRITE(LUPRI,'(A)')

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
      use codata
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      PARAMETER ( GFERMI = 2.22254D-14 )
      PARAMETER ( STHETAW2 = 0.2319D0 )
      PARAMETER ( AU2HZ = 6.579683920721D15 )
C
      LOGICAL   DOATOM(NUCDEP)
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),
     &          WORK(*)
C
#include "symmet.h"
#include "dgroup.h"
#include "ibtfun.h"
C
C
      FACTOR = (GFERMI*(D1-D4*STHETAW2))/(CVEL*SQRT(D2))
      FACTOR = (FACTOR*AU2HZ)/(D2*PI*CVEL*CVEL)
      FACTOR = FACTOR*1.0D9
C
      CALL AROUND('PV contribution to spin-spin coupling (in nHz)')
      WRITE(LUPRI,*)
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
          IZ = NINT(CHARGE(IATOM))
          ABUNDANCE = D0
          DO ISO = 1,5
            GVAL_IATOM = DISOTP(IZ,ISO,'GVAL')
            IF(GVAL_IATOM .NE. D0) GOTO 1
          ENDDO
 1        CONTINUE
          IZ = NINT(CHARGE(JATOM))
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
C/* Deck PRP_EFN */
      SUBROUTINE PRP_EFN(EXPVAL,EXP_CORR,CORR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by T. Saue Feb 26 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION EXPVAL(*),EXP_CORR(*),WORK(*)
      LOGICAL CORR
C
      KFRSAV = KFREE
C
C     Allocate memory
C
      CALL MEMGET('LOGI',KATOM,NUCIND  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEFN ,9*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTEMP,3*NUCDEP,WORK,KFREE,LFREE)
C
C     (-29: accept .SELECT specification)
      CALL SETATM(WORK(KATOM),NATOM,-29)
      CALL PRP_EFN1(EXPVAL,EXP_CORR,WORK(KEFN),WORK(KTEMP),WORK(KATOM),
     &              CORR,IPRPRP)
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
#include "implicit.h"
#include "priunit.h"
C     
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      LOGICAL   DOATOM(*), CORR
      DIMENSION EXPVAL(*),EXP_CORR(*),EFNVEC(3,NUCDEP,3),TEMP(3,NUCIND)
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
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
#include "implicit.h"
#include "priunit.h"
C     
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      DIMENSION EXPVAL(*), TEMP(3,NUCIND),EFNEL(3,NUCDEP)
      LOGICAL DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
      CALL DZERO(EFNEL,3*NUCDEP)
      CALL DZERO(TEMP,3*NUCIND)
      ITYP = 0
      DO 100 IATOM = 1, NUCIND
      IF (.NOT. DOATOM(IATOM)) GOTO 100
        DO 50 IREPC = 0, MAXREP
          DO 200 ICOOR = 1, 3
          IF (IBTAND(ISTBNU(IATOM),
     &        IBTXOR(IREPC,ISYMAX(ICOOR,1))).EQ.0) THEN
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
        IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          DO 220 I = 1, 3
            EFNEL(I,NATOM) = TEMP(I,IATOM)*
     &                       PT(IBTAND(ISYMAX(I,1),IREPC))
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
      SUBROUTINE DEF_RHONUC(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for electronic density at nuclei
C
C     Written by T.Saue May 8 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     INTTYP = 9: One-electron Fermi contact integrals
      CALL SETATM(WORK(KDATOM),NATOM,-9)
      CALL DEF_RHONU1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_RHONUC',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_RHONU1 */
      SUBROUTINE DEF_RHONU1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define electronic density at nuclei
C     Written by T.Saue May 8 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "gfac.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D1 = 1.0D0, D4 = 4.0D0, D3 = 3.0D0)
      PARAMETER (DFAC = D4*GFAC*PI/D3)
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   DOATOM(NUCIND)
C
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "dcbxpr.h"
#include "dcbexp.h"
C
#include "ibtfun.h"
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1/DFAC
      DO 100 IREP = 0, MAXREP
        DO 200 IATOM = 1, NUCIND
        IF (DOATOM(IATOM)) THEN
          IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
            PLABEL(1) = 'FC '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                        //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
            PNAME = 'Rho at nuc '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                            //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                            //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
            CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                   INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
            CALL OP1IND('RHONUC',IDUM,LEXPP,NEXPP,INDXPR,MAXEXP)
          END IF
        END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_rotg */
      SUBROUTINE DEF_ROTG(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Define operators for rotational g-tensor
C     DEF_SHIELD1 as skeleton
C
C     Written by T. Enevoldsen - Jul. 1998
C     Last revision: MI, jan.2003 - added from TEC's code
C                    MI, Sept.2006 - making functional
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
      CHARACTER PNAME*16, PLABEL(2)*8
      DIMENSION WORK(LWORK),PFAC(2),IOP(3)
      PARAMETER(DI2 = 0.5D0)
C
C
      CALL IZERO(IPROTG,6)
C     XAVECTOR
      IOP(1) = 5
C     YAVECTOR
      IOP(2) = 6
C     ZAVECTOR
      IOP(3) = 7
C
C     Define  operators
C     Conventional orbitals
         NPCOMP   = 2
         PFAC(1)  = DI2*CVAL
         PFAC(2)  = DI2*CVAL
C X componete of B-field
         PNAME     = 'X magnetic field'
         IPTYP     = IOP(1)
         PLABEL(1) = 'ZDIPLEN'
         PLABEL(2) = 'YDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('XFIELDA',IPROTG(1,1),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
         CALL OP1IND('XFIELDB',IPROTG(2,1),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
#if !defined (XH4)
C Y componete of B-field
         PNAME     = 'Y magnetic field'
         IPTYP     = IOP(2)
         PLABEL(1) = 'XDIPLEN'
         PLABEL(2) = 'ZDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('XFIELDA',IPROTG(1,2),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
         CALL OP1IND('XFIELDB',IPROTG(2,2),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
C Z componete of B-field
         PNAME     = 'Z magnetic field'
         IPTYP     = IOP(3)
         PLABEL(1) = 'YDIPLEN'
         PLABEL(2) = 'XDIPLEN'
         CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
         CALL OP1IND('XFIELDA',IPROTG(1,3),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
         CALL OP1IND('XFIELDB',IPROTG(2,3),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
#endif
C
C     Print section
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 'Rotational g-tensor:'
      CALL PRSYMB(LUPRI,'=',75,0)
C
C      Print reference
C
      WRITE(LUPRI,'(4X,A)')
     & 'Reference: T.Enevoldsen et al, JCP 114 (2001) p.84'
     
      WRITE(LUPRI,'(1X,A)')
     & 'NOTE: You must set .USECM in the *NMR section '//
     & ' to get correct values of rot_g_tensor.'
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)') '** A operators **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3
        INDOP1 = IPROTG(1,I)
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
      WRITE(LUPRI,'(1X,A)') '** B operators **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,3
        INDOP1 = IPROTG(2,I)
        IF(INDOP1.GT.0) THEN
          INDXPR = LLRBPU(INDOP1)
          CALL WRIXPR(INDOP1,INDXPR)
        ENDIF
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ROTG */
      SUBROUTINE PRP_ROTG(ABLRF,ABCNV,WORK,KFREE,LFREE)
C***********************************************************************
C     
C    Output routine for rotational g-tensor
C
C   Written by T.Enevoldsen Jul. 1998
C     
C  Revisions: MI, jan.2003 - added from TEC's 'left' code
C             MI, Prievidza/Odense, August 2006
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "orgcom.h"
#include "dcborb.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbxlr.h"
CMI ... added
#include "nuclei.h"
#include "symmet.h"
#include "molgfa.h"
#include "iratdef.h"

      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT),ABCNV(NBFREQ,NLRBPT),
     &     WORK(*),IBUFA(3),IBUFB(3)
      INTEGER ICONVXYZ
      CHARACTER AX*1,BX*1

      CALL QENTER('PRP_ROTG')
C     
      CALL TITLER('DIRAC  - ROTATIONAL g-tensor','*',116)
C     
      KFRSAV = KFREE
C     
#if defined (XH4)
      WRITE(LUPRI,'(A)') '*** WARNING *** Using code for XH4 !'
#endif
C     
C     *** Extract g-tensor ***
C     
      NA = 0
      NB = 0
      DO I = 1, 3
         IF(IPROTG(1,I).GT.0) THEN
            NA = NA + 1
            IBUFA(NA) = IPROTG(1,I)
         ENDIF
         IF(IPROTG(2,I).GT.0) THEN
            NB = NB + 1
            IBUFB(NB) = IPROTG(2,I)
         ENDIF
      END DO

      DO IFREQ = 1,NBFREQ

         CALL PRSYMB(LUPRI,'-',72,2)

         IF(BFREQ(IFREQ).NE.0.0D0) THEN
            WRITE(LUPRI,'(A,F15.8,2X,A)') 
     &      '*** Frequency    :',BFREQ(IFREQ),'a.u.'
         ENDIF

         DO IB = 1,NB
            READ (PRPNAM(LLRBPU(IBUFB(IB))),'(A1)') BX

            KB = IBUFB(ICONVXYZ(BX))
            KB = IPTAX(KB,2)

            DO IA = 1,NA

               READ (PRPNAM(LLRAPU(IBUFA(IA))),'(A1)') AX

               KA = IBUFA(ICONVXYZ(AX))
               KA = IPTAX(KA,2)

               VALUE = -4*ABLRF(IFREQ,KA,KB)

               ELMGF (KA,KB) = ABLRF(IFREQ,IBUFA(IA),IBUFB(IB))
               ELMGF (KB,KA) = ELMGF (KA,KB)

       WRITE(LUPRI,'(5A,F15.8,2X,A)') '@ g(',AX,',',BX,') =',
     &             VALUE,'a.u.'
            ENDDO
         ENDDO
C         
         IF (IPRPRP.GE.7) THEN
         CALL HEADER('Electronic contribution to molecular'//
     &                  ' g-factor',-1)
         CALL OUTPUT(ELMGF,1,3,1,3,3,3,1,LUPRI)
         write(lupri,*) 'POLPRI output:'
         CALL POLPRI(ELMGF,'   ',-2)
         ENDIF

      END DO

      CALL MEMCHK('PRP_ROTG 1',WORK,1)

      IGEOMSIZE = 3*(NATOMS+NFLOAT)
      CALL MEMGET('REAL',KGEOM  ,IGEOMSIZE, WORK,KFREE,LFREE)

      CALL MEMCHK('PRP_ROTG 2',WORK,1)

      CALL MEMGET('REAL',KMASS  ,  NATOMS ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNAT   ,  NATOMS ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNUMIS ,  NATOMS ,WORK,KFREE,LFREE)


      DO I = 1, MXCENT
            ISOTOP(I) = 1
      END DO
      CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),ISOTOP,
     &            IPRPRP)

       CALL DCOPY(3,CMXYZ,1,GAGORG,1)
       CALL NUCMGF(WORK(KGEOM),GAGORG,LUPRI,IPRPRP)

CMI ... from DALTON
         JATOM  = NATOMS + NFLOAT
         CALL MEMGET('REAL',KGEOM  ,3*JATOM ,WORK,KFREE,LFREE)
C        KGEOM  = 1
C        KAMASS = KGEOM + 3*JATOM
         CALL MEMGET('REAL',KAMASS ,JATOM ,WORK,KFREE,LFREE)

         KLAST  = KAMASS + JATOM
C        LWRK   = LWORK - KLAST + 1
         LWRK   = LFREE - KLAST + 1
CMI      IF (KLAST .GT. LWORK) CALL
         IF (KLAST .GT. LFREE) CALL
     &        QUIT('PRP_ROTG...error allocations !!! ')

        CALL MEMCHK('prp_rotg 0.1',WORK,1)

      DO I = 1, MXCENT
            ISOTOP(I) = 1
      END DO
      CALL MGFRES(WORK(KGEOM),WORK(KAMASS),JATOM,IPRPRP)


C ... release the memmory...
      CALL MEMREL('GETCON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('PRP_ROTG')
      RETURN
      END

C  /* Deck nucmgf */
      SUBROUTINE NUCMGF(GEOM,ORIGIN,LUPRI,IPRINT)
C
C     Calculate nuclear contribution to molecular g-factors
C     K.Ruud, June-94
C
CMI Added from dalton-2.0, May 2006, Strasbourg
C
C calling in DALTON:  IF (MOLGFA) CALL NUCMGF(WORK(KGEOM),GAGORG,LUPRI,IPRDEF)
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <nuclei.h>
      DIMENSION ORIGIN(3), GEOM(3*(NATOMS+NFLOAT))
#include <symmet.h>
#include <molgfa.h>
#include <ibtfun.h>
C

      IF (IPRINT.GE.5) THEN
       CALL HEADER(
     & 'NUCMGF: Entering parameters to calculate nuclear'//
     & ' contribution to molecular g-factor',-1)
       write(lupri,*) 'ORIGIN=',ORIGIN(1),ORIGIN(2),ORIGIN(3)
       do i=1, (3*(NATOMS+NFLOAT))
        write(lupri,*) 'Entering geom(',i,')=',GEOM(i)
       enddo
      ENDIF

      CALL DZERO(GFANUC,9)
      KK = 1
      NATOM = 1
      IX = IPTAX(1,2)
      IY = IPTAX(2,2)
      IZ = IPTAX(3,2)
      DO 10 I = 1, NUCIND
         IF (IPRINT.GE.6) THEN
          write(lupri,*) I,'. nuclei with charge=',CHARGE(I)
         ENDIF
         DO 20 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(I)) .EQ. 0) THEN
               XCOOR = GEOM(KK) - ORIGIN(1)
               YCOOR = GEOM(KK + 1) - ORIGIN(2)
               ZCOOR = GEOM(KK + 2) - ORIGIN(3)
               R2    = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
               GFANUC(IX,IX) = GFANUC(IX,IX) + CHARGE(I)*(YCOOR*YCOOR
     &                                       + ZCOOR*ZCOOR)
               GFANUC(IY,IY) = GFANUC(IY,IY) + CHARGE(I)*(XCOOR*XCOOR
     &                                       + ZCOOR*ZCOOR)
               GFANUC(IZ,IZ) = GFANUC(IZ,IZ) + CHARGE(I)*(XCOOR*XCOOR
     &                                       + YCOOR*YCOOR)
               GFANUC(IX,IY) = GFANUC(IX,IY) - CHARGE(I)*XCOOR*YCOOR
               GFANUC(IX,IZ) = GFANUC(IX,IZ) - CHARGE(I)*XCOOR*ZCOOR
               GFANUC(IY,IZ) = GFANUC(IY,IZ) - CHARGE(I)*YCOOR*ZCOOR
               KK = KK + 3
            END IF
 20      CONTINUE
 10   CONTINUE
      GFANUC(IY,IX) = GFANUC(IX,IY)
      GFANUC(IZ,IX) = GFANUC(IX,IZ)
      GFANUC(IZ,IY) = GFANUC(IY,IZ)
CMI   ... final output ...
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Nuclear contribution to molecular g-factor',-1)
         CALL OUTPUT(GFANUC,1,3,1,3,3,3,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck mgfres */
C
C
C    Output routine for molecular g-factor
C
C
      SUBROUTINE MGFRES(GEOM,AMASS,KATOM,IPRINT)
      use codata
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
C#include <mxorb.h>
#include <maxorb.h>
C#include <pmass.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0,
     &           THRESH = 1.0D-5)
      LOGICAL PLANAR, LINEAR, D12, D23
CC    DIMENSION ANGMOM(3), TINERT(9), OMEGA(3,3), WORK(LWORK),
      DIMENSION ANGMOM(3), TINERT(9), OMEGA(3,3), 
     &          EIGVAL(3), EIGVEC(3,3),
     &          GTRAN1(3,3), GTRAN2(3,3), GTRAN3(3,3), GEOM(KATOM,3),
     &          AMASS(KATOM)
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
C#include <suscpt.h>
#include "dcbsusc.h"
C#include "cbiher.h"
#include <molgfa.h>
#include <orgcom.h>
#include <ibtfun.h>
C

      CALL QENTER('MGFRES')

C     CALL MEMCHK('MGFRES 0.9',WORK,1)

      DO I = 1, MXCENT
          ISOTOP(I) = 1
      END DO

      NOLOND = .TRUE.

      D12 = .FALSE.
      D23 = .FALSE.
C
      IF (IPRINT .GE. 0) THEN
        CALL TITLER('ABACUS - MOLECULAR g-TENSOR','*',124)
      ENDIF

      CALL DZERO(TOTMGF,9)
C
C     Add terms to diamagnetic molecular g factor
C
C     write(lupri,*) 'NOLOND=',NOLOND,' IPRINT=',IPRINT

      IF (.NOT. NOLOND) THEN
         DO 10 I = 1, 3
         DO 10 J = 1, 3
            GFACDIA(I,J) = GFACDIA(I,J) + SUS2EL(I,J) + SUSFS(I,J)
     &                  + SUSFSY(I,J)
 10      CONTINUE 
      END IF
      DO 200 I = 1, 3
         ANGMOM(I) = D1
         DO 200 J = 1, 3
            ELMGF(I,J)  = D4*ELMGF(I,J)
            GFACDIA(I,J) = D4*GFACDIA(I,J)
            TOTMGF(I,J) = GFANUC(I,J) + ELMGF(I,J) + GFACDIA(I,J)
 200  CONTINUE
C
      IF (IPRINT .GE. 2) THEN
         CALL HEADER(
     &   'Nuclear contribution * proton mass * moment of inertia',-1)
         CALL POLPRI(GFANUC,'   ',-2)
         AVENUC = (GFANUC(1,1) + GFANUC(2,2) + GFANUC(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVENUC
C
         CALL HEADER(
     &   'Diamagnetic contribution * proton mass * moment of inertia',
     &   -1)
C        CALL POLPRI(GFACDI,'   ',-2)
         CALL POLPRI(GFACDIA,'   ',-2)
CMI      AVEDIA = (GFACDI(1,1) + GFACDI(2,2) + GFACDI(3,3))/D3
         AVEDIA = (GFACDIA(1,1) + GFACDIA(2,2) + GFACDIA(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVEDIA
C
         CALL HEADER(
     &   'Electronic contribution * proton mass * moment of inertia',-1)
         CALL POLPRI(ELMGF,'   ',-2)
         AVEEL = (ELMGF(1,1) + ELMGF(2,2) + ELMGF(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVEEL
C
         CALL HEADER(
     &   'Total molecular g-factor * proton mass * moment of inertia',
     &   -1)
         CALL POLPRI(TOTMGF,'   ',2)
         AVETOT = (TOTMGF(1,1) + TOTMGF(2,2) + TOTMGF(3,3))/D3
         WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVETOT
      END IF
C
C     Determine principal moments of inertia
C
      JATOM = 1
      DO 100 IATOM = 1, NUCIND
         NATTYP = NINT(CHARGE(IATOM))
         NUMIS  = ISOTOP(IATOM)
         DO 110 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               AMASS(JATOM) = DISOTP(NATTYP,NUMIS,'MASS')

          IF (IPRINT.GE.9) THEN
           write(lupri,*) 'NATTYP,NUMIS,DISOTP=',NATTYP,NUMIS,
     &     DISOTP(NATTYP,NUMIS,'MASS')
           write(lupri,*)
     &    'AMASS(',JATOM,')=',AMASS(JATOM) 
          ENDIF

               DO 120 ICOOR = 1, 3
                  GEOM(JATOM,ICOOR) = PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
     &                               *CORD(ICOOR,IATOM) - CMXYZ(ICOOR)
             if (IPRINT.GE.8) THEN
               write(lupri,*)
     &         'geom(',JATOM,';',ICOOR,')=',
     &         GEOM(JATOM,ICOOR)
             endif

 120           CONTINUE
               JATOM = JATOM + 1
            END IF
 110     CONTINUE
 100  CONTINUE

      CALL WLKDIN(GEOM,AMASS,KATOM,ANGMOM,TINERT,OMEGA,EIGVAL,EIGVEC,
     &            .TRUE.,PLANAR,LINEAR,IPRINT)

C
C     Transform molecular g-factor to principal axis system, and multiply by
C     inverse moment of inertia and proton mass
C
      IF (IPRINT .GE. 0) THEN
         IF (NOLOND) THEN
            WRITE (LUPRI,'(/,1X,A)')' Calculation without '//
     &                               'London orbitals.'
         ELSE
            WRITE (LUPRI,'(/,1X,A)')' London orbitals used.'
         END IF
      END IF
C
C     All components of the moment of inertia not defined for linear molecules
C
      IF (LINEAR) THEN
         IF (ABS(TOTMGF(1,1)) .GT. THRESH) THEN
            AVETOT = TOTMGF(1,1)*PMASS*EIGVAL(3)
CMI         AVEDIA = GFACDI(1,1)*PMASS*EIGVAL(3)
            AVEDIA = GFACDIA(1,1)*PMASS*EIGVAL(3)
            AVEEL  = ELMGF (1,1)*PMASS*EIGVAL(3)
            AVENUC = GFANUC(1,1)*PMASS*EIGVAL(3)
            CALL DZERO(GTRAN,9)
            GTRAN(1,1) = AVETOT
         ELSE
            AVETOT = TOTMGF(2,2)*PMASS*EIGVAL(3)
CMI         AVEDIA = GFACDI(2,2)*PMASS*EIGVAL(3)
            AVEDIA = GFACDIA(2,2)*PMASS*EIGVAL(3)
            AVEEL  = ELMGF (2,2)*PMASS*EIGVAL(3)
            AVENUC = GFANUC(2,2)*PMASS*EIGVAL(3)
            CALL DZERO(GTRAN,9)
            GTRAN(1,1) = AVETOT
         END IF
CMI       IF (IPRDEF .GE. 0) THEN
        IF (IPRINT .GE. 0) THEN
            WRITE (LUPRI,'(2X,A,F12.6)') 'Moment of inertia :',
     &           D1/EIGVAL(3)
            CALL HEADER('Molecular g-factor',1)
            WRITE(LUPRI,'(4(1X,A,F12.8,/))')
     &           ' Isotropic g-factor          :',AVETOT,
     &           ' Nuclear part of g-factor    :',AVENUC,
     &           ' Diamagnetic part of g-factor:',AVEDIA,
     &           ' Electronic part of g-factor :',AVEEL
            WRITE (LUPRI,'(1X,A)') 'The molecule is linear'
         END IF
      ELSE
         IF (IPRINT .GE. 0) THEN
            CALL HEADER
     &    ('Principal moments of inertia (a.u.) and principle axis',0)
            WRITE (LUPRI,'(3X,A,F12.6,6X,3F12.6)')
     &           'IA',D1/EIGVAL(1), (EIGVEC(I,1),I = 1,3),
     &           'IB',D1/EIGVAL(2), (EIGVEC(I,2),I = 1,3),
     &           'IC',D1/EIGVAL(3), (EIGVEC(I,3),I = 1,3)
         END IF
         CALL DZERO(GTRAN,9)
         CALL DZERO(GTRAN1,9)
         CALL DZERO(GTRAN2,9)
         CALL DZERO(GTRAN3,9)
         DO 50 I = 1, 3
            DO 50 J = 1, 3
               DO 60 IX = 1, 3
                  DO 60 IY = 1, 3
                     JX = IPTAX(IX,2)
                     JY = IPTAX(IY,2)
                     CONST = EIGVAL(I)*PMASS
                     GTRAN(I,J)  = GTRAN(I,J) + EIGVEC(IX,I)
     &                            *TOTMGF(JX,JY)*EIGVEC(IY,J)*CONST
                     GTRAN1(I,J) = GTRAN1(I,J) + EIGVEC(IX,I)
     &                            *GFANUC(JX,JY)*EIGVEC(IY,J)*CONST
                     GTRAN2(I,J) = GTRAN2(I,J) + EIGVEC(IX,I)
     &                            *ELMGF(JX,JY)*EIGVEC(IY,J)*CONST
                     GTRAN3(I,J) = GTRAN3(I,J) + EIGVEC(IX,I)
     &                            *GFACDIA(JX,JY)*EIGVEC(IY,J)*CONST
 60            CONTINUE
 50      CONTINUE
         DIFF1 = ABS(GTRAN(1,1) - GTRAN(2,2))
         DIFF2 = ABS(GTRAN(2,2) - GTRAN(3,3))
         IF (DIFF1 .LT. THRESH) D12 = .TRUE.
         IF (DIFF2 .LT. THRESH) D23 = .TRUE.
         AVETOT = (GTRAN(1,1)  + GTRAN(2,2)  + GTRAN(3,3) )/D3
         AVENUC = (GTRAN1(1,1) + GTRAN1(2,2) + GTRAN1(3,3))/D3
         AVEEL  = (GTRAN2(1,1) + GTRAN2(2,2) + GTRAN2(3,3))/D3
         AVEDIA = (GTRAN3(1,1) + GTRAN3(2,2) + GTRAN3(3,3))/D3
         IF (IPRINT .GE. 0) THEN
         IF (D12 .AND. D23) THEN
            WRITE (LUPRI,'(4(/,1X,A,F12.8))')
     &           ' Molecular g-factor      :  ',AVETOT,
     &           ' Nuclear contribution    :  ',AVENUC,
     &           ' Diamagnetic contribution:  ',AVEDIA,
     &           ' Electronic contribution :  ',AVEEL
            WRITE (LUPRI,'(/,1X,A)') ' Molecular g-factor is spherical.'
         ELSE IF (D12 .OR. D23) THEN
            IF (D12) THEN
               PAR = GTRAN(3,3)
               PER = (GTRAN(1,1) + GTRAN(2,2))/D2
            ELSE
               PAR = GTRAN(1,1)
               PER = (GTRAN(2,2) + GTRAN(3,3))/D2
            END IF
            ANI = PAR - PER
            CALL HEADER('Molecular g-factor',1)
            WRITE (LUPRI,'(6(1X,A,F12.8,/))')
     &           ' Isotropic g-factor       : ',AVETOT,
     &           ' Nuclear contribution     : ',AVENUC,
     &           ' Electronic contribution  : ',AVEEL,
     &           ' Parallel component:        ',PAR,
     &           ' Perpendicular component:   ',PER,
     &           ' Anisotropy:                ',ANI
            WRITE (LUPRI,'(1X,A)') ' Molecular g-factor is cylindrical.'
            CALL HEADER('Molecular g-factor in principal axis system '//
     &                  '(a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN(I,1),I = 1,3),
     &         'B', (GTRAN(I,2),I = 1,3),
     &         'C', (GTRAN(I,3),I = 1,3)
            CALL HEADER('Nuclear contribution '//
     &                  'in principle axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN1(I,1),I = 1,3),
     &         'B', (GTRAN1(I,2),I = 1,3),
     &         'C', (GTRAN1(I,3),I = 1,3)
            CALL HEADER('Electronic contribution '//
     &                  'in principle axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN2(I,1),I = 1,3),
     &         'B', (GTRAN2(I,2),I = 1,3),
     &         'C', (GTRAN2(I,3),I = 1,3)
            CALL HEADER('Diamagnetic contribution '//
     &                  'in principle axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN3(I,1),I = 1,3),
     &         'B', (GTRAN3(I,2),I = 1,3),
     &         'C', (GTRAN3(I,3),I = 1,3)
         ELSE
            ANI1 = GTRAN(1,1) - (GTRAN(2,2) + GTRAN(3,3))/D2
            ANI2 = GTRAN(2,2) - (GTRAN(1,1) + GTRAN(3,3))/D2
            CALL HEADER('Molecular g-tensor',1)
            WRITE (LUPRI,'(6(1X,A,F12.8,/))')
     &           ' Isotropic g-tensor       : ',AVETOT,
     &           ' Nuclear contribution     : ',AVENUC,
     &           ' Diamagnetic contribution : ',AVEDIA,
     &           ' Electronic contribution  : ',AVEEL,
     &           ' 1st anisotropy:            ',ANI1,
     &           ' 2nd anisotropy:            ',ANI2
            CALL HEADER('Molecular g-tensor in principal axis system '//
     &                  '(a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN(I,1),I = 1,3),
     &         'B', (GTRAN(I,2),I = 1,3),
     &         'C', (GTRAN(I,3),I = 1,3)
            CALL HEADER('Nuclear contribution '//
     &                  'in principle axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN1(I,1),I = 1,3),
     &         'B', (GTRAN1(I,2),I = 1,3),
     &         'C', (GTRAN1(I,3),I = 1,3)
            CALL HEADER('Electronic contribution '//
     &                  'in principle axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN2(I,1),I = 1,3),
     &         'B', (GTRAN2(I,2),I = 1,3),
     &         'C', (GTRAN2(I,3),I = 1,3)
            CALL HEADER('Diamagnetic contribution '//
     &                  'in principle axis system (a.u.):',1)
            WRITE (LUPRI,'(3X,A,6X,3F12.8)')
     &         'A', (GTRAN3(I,1),I = 1,3),
     &         'B', (GTRAN3(I,2),I = 1,3),
     &         'C', (GTRAN3(I,3),I = 1,3)
         END IF
         END IF
      END IF

      CALL QEXIT('MGFRES')
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
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "cbihr1.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbprl.h"
#include "dummy.h"
#include "ccom.h"
C
      DIMENSION WORK(LWORK),VMAT(*)
      LOGICAL IS_SQUARE, IS_PCTRA
      CHARACTER*8 PRPLABEL
      CHARACTER*4 NUMSTR
C
      CALL QENTER('PROP2BSS')
#include "memint.h"

C     ... allocate for LL-only operator 
      N2BBASXQ_L = NTBAS(1)*NTBAS(1)*NZ
      CALL MEMGET('REAL',KPRLL,N2BBASXQ_L,WORK,KFREE,LFREE)

      IF (IPRINT.GE.3) THEN
        CALL HEADER('*** Output from PROP2BSS ***',-1)
      ENDIF

      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',FILE='AOPROPER')
      REWIND (LU1INT)

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

        IF (IPRPTIM(I).EQ.0) THEN
          IDIMENS = N2BBASX
          IS_SQUARE = .TRUE.
        ELSE IF (IPRPTIM(I).EQ.1.OR.IPRPTIM(I).EQ.-1) THEN
          IDIMENS = NNBBASX
          IS_SQUARE = .FALSE.
        ELSE
          CALL QUIT('PROP2BSS: wrong parameter ! ')
        ENDIF
C
C       Memory allocation
C
        CALL MEMGET('REAL',KPTRI  ,NZ*IDIMENS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KONEIN ,   IDIMENS,WORK,KFREE,LFREE)
        CALL MEMGET('LOGI',KFIRST ,NZ        ,WORK,KFREE,LFREE)

C       ... get the operator either in the full four-component form, 
C       .... or already in  the picture change transformed shape
C       ...   (occupying in fact only  LL block of the four-component matrix)
        INDXPR = I
        IF (IS_SQUARE) THEN
          CALL PRPMSAO(INDXPR,IDIMENS,DUMMY,VMAT,
     &         WORK(KONEIN),WORK(KFIRST),WORK(KFREE),
     &         LFREE,IPRINT)   
        ELSE
          CALL PRPMSAO(INDXPR,IDIMENS,WORK(KPTRI),VMAT,
     &         WORK(KONEIN),WORK(KFIRST),WORK(KFREE),
     &         LFREE,IPRINT)   
        ENDIF
C     ... deallocate the memory...
        CALL MEMREL('KPTRI',WORK,KWORK,KPTRI,KFREE,LFREE)

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
C     Close AOPROPER
C
      CLOSE(LU1INT,STATUS='KEEP')
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
      SUBROUTINE DEF_OPTROT(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for optical rotation; based on DEF_SHIELD
C     and DEF_EPOLAR
C     Written by T. Saue and K. Ruud Jan 8 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION WORK(LWORK)
C
#include "memint.h"
      CALL DEF_OPTROT1(WORK(KFREE),LFREE,IPRINT)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_optrot1 */
      SUBROUTINE DEF_OPTROT1(WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define operator for optical rotation; based on DEF_SHIELD
C     and DEF_EPOLAR
C     Written by T. Saue and K. Ruud Jan 8 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0,DMI2 = -0.5D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "gfac.h"
#include "pi.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
#include "dcbexp.h"
#include "dcbnmr.h"
#include "chrxyz.h"
      CHARACTER PNAME*16, PLABEL(3)*8,LABEL*4,LABELC*7
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"
C
      CALL QENTER('DEF_OPTROT1')
C
C     First define left-hand operators
C
      PNAME     = 'X-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(1,1) = INDOP1
C
      PNAME     = 'Y-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(2,1) = INDOP1
C
      PNAME     = 'Z-Dipole length '
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = DM1
      PLABEL(1) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('ZEPOLA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      IPEPOLAR(3,1) = INDOP1
C
C     Define right hand operators
C     (PFAC(1:2) for Dirac mag. field, PFAC(3) for g_1 QED correction)
C

      IF (LONDON) THEN
        CALL DEF_D1HBLOND(1,WORK,LWORK,IPRINT)
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
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XSHIELB',IPSHIELD(1+MXCOOR),
     &             LLRBPU,NLRBPT,INDXPR,MAXLLR)
#if !defined (XH4)
C Y component of B-field
      PNAME     = 'Y magnetic field'
      IPTYP     = IOP(2)
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YSHIELB',IPSHIELD(2+MXCOOR),
     &      LLRBPU,NLRBPT,INDXPR,MAXLLR)
C Z component of B-field
      PNAME     = 'Z magnetic field'
      IPTYP     = IOP(3)
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
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

      CALL QEXIT('DEF_OPTROT1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_optrot */
      SUBROUTINE PRP_OPTROT(ABLRF,ABCNV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for optical rotation, based on LNROUT from DALTON
C     Written by T. Saue in collaboration with K. Ruud Feb 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxlr.h"
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),WORK(*)
      KFRSAV = KFREE
      CALL MEMGET('REAL',KTEMP,9,WORK,KFREE,LFREE)
      CALL PRP_OPTRO1(ABLRF,ABCNV,WORK(KTEMP),
     &            WORK,KFREE,LFREE)
      CALL MEMREL('PRP_OPTROT',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prp_optrot */
      SUBROUTINE PRP_OPTRO1(ABLRF,ABCNV,TMPMAT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for optical rotation, based on LNROUT from DALTON
C     Written by T. Saue in collaboration with K. Ruud Feb 2008
C
C***********************************************************************
      use codata
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D3=3.0D0,
     &          FACTOR = (288.0D-30)*(PI**2)*XFMOL*(XTANG**4))
#include "mxcent.h"
C
#include "dcbprp.h"
#include "dcbxlr.h"
#include "nuclei.h"
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),
     &          ABCNV(NBFREQ,NLRBPT),TMPMAT(3,3),WORK(*)
      EXTERNAL TOTMASS
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0,
     &          DMI2 = -0.5D0,DI2 = +0.5D00,DM1 = - 1.0D0)
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "gfac.h"
#include "pi.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "dcbexp.h"
#include "dcbnmr.h"
#include "chrxyz.h"
      CHARACTER PNAME*16, PLABEL(3)*8,LABEL*4,LABELC*7,chcomp*1
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
      integer a, b, c
#include "chrnos.h"
#include "ibtfun.h"
#include "dcbxpr.h"
#include "dcbsusc.h"
#include "memint.h"

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
        CALL WRIXPR(0,IPLONDON(J,icomp))
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
      SUBROUTINE DEF_EFF_DENS(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C**********************************************************
C     Define effective density operator
C     Written by R. van Meer September 2010
C**********************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"

      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
C
#include "memint.h"
      CALL MEMGET('LOGI',KDATOM,NUCIND,WORK,KFREE,LFREE)
C     (-99: accept .SELECT specification)
      CALL SETATM(WORK(KDATOM),NATOM,-99)
      CALL DEF_EFF_DENS1(PNAME,PLABEL,PFAC,WORK(KDATOM),
     &                WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('DEF_EFF_DENS',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN  
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DEF_EFF_DENS1(PNAME,PLABEL,PFAC,DOATOM,WORK,LWORK,
     &               IPRINT)
C**********************************************************
C     Define effective density operator
C     Written by R. van Meer September 2010
C**********************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D1 = 1.0D0)
C
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3)
      LOGICAL   DOATOM(NUCIND)
C
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "dcbxpr.h"
#include "dcbexp.h"
C
#include "ibtfun.h"
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      DO 100 IREP = 0, MAXREP
        DO 200 IATOM = 1, NUCIND
          IF (DOATOM(IATOM)) THEN
            IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
              PLABEL(1) = 'ED '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                          //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                          //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              PNAME = 'Eff._Dens '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                              //CHRNOS(IPTNUC(IATOM,IREP)/10)
     &                              //CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                     INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
              CALL OP1IND('EFFDEN',IDUM,LEXPP,NEXPP,INDXPR,MAXEXP)
            END IF
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
