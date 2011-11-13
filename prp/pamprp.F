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
C VERSION : $Revision: 1.60.2.4 $
C DATE    : $Date: 2004/09/01 13:26:01 $
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
     &     3X,'First order properties',/,
     &     3X,'Trond Saue',//,
     &     3X,'KR-RPA second-order properties',/,
     &     3X,'Hans Joergen Aa. Jensen and Trond Saue',//,
     &     3X,'KR-QR third order properties',/,
     &     3X,'Patrick Norman and Hans Joergen Aa. Jensen',//,
     &     3X,'Molecular gradient',/,
     &     3X,'Joern Thyssen',//,
     &     3X,'First-order ESR properties',/,
     &     3X,'Hans Joergen Aa. Jensen and Anette Noerager',//,
     &     3X,'Additional contributions from:',/,
     &     3X,'Thomas Enevoldsen and Luuk Visscher'//)
C
C
      DO I = 1, NPRPWF
         WF = PRPWF(I)
C
         IF (WF .EQ. 'DHF' .OR. WF .EQ. 'DFT') THEN
            CALL TDHFORB
            CALL SETDC2(0)
         ELSE IF (WF .EQ. 'KRMC') THEN
            CALL TKRMCORB
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
#if defined (MOD_ESR)
      IF (DOESR) THEN
         KFRSAV = KFREE
         CALL ESR_INIT(WF,IPRESR,NCIESR,WORK,KFREE,LFREE)
C
         CALL MEMGET('REAL',KESR,2*NCIESR*NCIESR*NESRP,WORK,KFREE,LFREE)
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
#if defined (MOD_XPP)
      IF (DOXPP) THEN
         IF (WF. EQ. 'DHF') THEN
            CALL PRPXPP(WORK,KPPF,KOMEGA,KPPCNV,KPPO,KFREE,LFREE)
         ELSE
            CALL QUIT('*** ERROR in PAMPR1 *** ' //
     &           'PP excitation energies not implemented for ' // WF )
         END IF
      END IF
#endif
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
      IF (DOXQR .OR. DOTPA) THEN
         IF (WF. NE. 'DHF') 
     &        CALL QUIT('*** ERROR in PAMPR1 *** ' //
     &        'Quadratic response not implemented for ' // WF )
         IF (DOXQR) CALL PRPXQR(WORK(KFREE),LFREE)
#if defined (MOD_XPP)
         IF (DOTPA) CALL PRPTPA(WORK(KFREE),LFREE)
#endif
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
      IF (DOEXP .OR. DOXLR) THEN
         CALL PRPOUT(WORK(KEXP),WORK(KLRF),WORK(KLRCNV),
     &            WORK,KFREE,LFREE)
         IF (DOCCM) THEN
            CALL PRPOUT_CORR(WORK(KEXP),WORK,KFREE,LFREE)
         END IF
      END IF
C
      CALL MEMREL('PAMPRP',WORK,KWORK,KWORK,KFREE,LFREE)
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
C     Written by T.Saue - June 28 1996
C     Last revision: June 28 - tsaue
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
      IPEPOLAR(1,2) = INDOP1
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
C
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
              PLABEL(1) = 'NEF '//CHRNOS(ISCOOR/10)
     &                          //CHRNOS(MOD(ISCOOR,10))//'  '
              PNAME    = 'NEF   '//CHRNOS(ISCOOR/10)
     &                           //CHRNOS(MOD(ISCOOR,10))
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
      IPTYP     = 9
      NPCOMP    = 1
      PFAC(1)   = D1
      DO 100 IATOM = 1, NUCIND
         IF (DOATOM(IATOM)) THEN
            DO 200 IREP = 0, MAXREP
               IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
                  ITYP = ITYP + 1
                  PNAME = 'PVC: '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                    //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                    CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
                  PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
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
      IF (TRLESS) THEN
         CALL SETATM(WORK(KDATOM),NATOM,IDUM)
      ELSE
         CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
            IF (ISCOR1 .GT. 0) THEN
C Is ITYP used anywhere? /  Feb 25 1998 - jth
C             ITYP = ITYP + 1
              PNAME    = 'NUCSPI'//CHRNOS(ISCOR1/10)
     &                         //CHRNOS(MOD(ISCOR1,10))
     &                //NAMN(IATOM1)//REP(IREPO)//CHRXYZ(ICOOR1)
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
              PLABEL(1)= 'NEF '//CHRNOS(IS1/10)
     &                     //CHRNOS(MOD(IS1,10))//'  '
              PLABEL(2)= 'NEF '//CHRNOS(IS2/10)
     &                     //CHRNOS(MOD(IS2,10))//'  '
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
      WRITE(LUPRI,'(1X,A)') 'NMR indirect spin-spin coupling:'
      CALL PRSYMB(LUPRI,'=',75,0)
      IF(TRIAB) THEN
        WRITE(LUPRI,'(A)')
     &    '*** WARNING ! Assuming symmetry J(AB) = J(BA)'
      ENDIF
      WRITE (LUPRI,'(A,E9.2)')
     &   ' Spin-spin couplings printed for atoms with '//
     &   'abundance greater than :', ABUND
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)') '** A operators **'
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
      WRITE(LUPRI,'(1X,A)') '** B operators **'
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
         IF (SHIELD) THEN
            DOEXP = .TRUE.
            NSTDIA = .TRUE.
         END IF
         IF (SPNSPN) THEN
            DOEXP = .TRUE.
            DSO = .TRUE.
         END IF
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
      IF(MSUSCP) THEN
        CALL DEF_MSUSCP(WORK,LWORK)   
      ENDIF
C
C     Nuclear shielding or ESR g-tensor or ESR anisotropy
C
      IF (SHIELD .OR. ESRGTENS .OR. ESRANISO) THEN
        CALL DEF_SHIELD(WORK,LWORK,IPRPRP)
      ENDIF
C
C     Diamagnetic contribution to shieldings
C
      IF(NSTDIA) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(1X,A)')
     &    'NMR diamagnetic shielding contribution:'
        CALL PRSYMB(LUPRI,'=',75,0)
CMI     CALL DEF_NSTDIA(PNAME,PLABEL,PFAC,WORK,LWORK,IPRPRP)
        CALL DEF_NSTDIA(PNAME,PLABEL,PFAC,WORK,LWORK,IPREXP)
      ENDIF
C
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
C     Electric polarizabilities
C
      IF(EPOLAR) THEN
C        CALL PRP_EPOLAR(ABLRF,ABCNV,WORK,KFREE,LFREE)
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
C     Memory deallocation
C     -------------------
C
      CALL MEMREL('PRPOUT',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
ctec      CALL PRSYMB(LUPRI,'-',76,0)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prpout_corr */
      SUBROUTINE PRPOUT_CORR(EXPVAL,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output of static 1st order properties at correlated level
C
C     Needs to be generalized
C
C     Written by Joost van Stralen June 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      DIMENSION EXPVAL(*),WORK(*)
      LOGICAL TOBE
#include "dcbprp.h"
#include "dcbexp.h"
C
      CALL PRSYMB(LUPRI,'-',76,0)
C
      CALL TITLER('MP2 expectation values','*',124)
      WRITE(LUPRI,'(14X,A)')
     &     'MP2 expectation values are documented in:',
     &     'J. N. P. van Stralen, L. Visscher, C. V. Larsen',
     &     'and H. J. Aa Jensen, Chem. Phys. in press.',
     &     ' '
C
      KFRSAV=KFREE
C
C     Memory Allocation
C     -----------------
C
      CALL MEMGET('REAL',KEXP_CORR,NEXPP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEXP_TOT,NEXPP,WORK,KFREE,LFREE)
C
C     Read correlation contribution to expectation value
C     ==================================================
C
      CALL READEXPV(WORK(KEXP_CORR),TOBE)
      IF (.NOT.TOBE) THEN
         CALL MEMREL('PRPOUT_CORR',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
         RETURN
      ENDIF
C
C     Calculate total expectation value
C     (HF + correlation)
C     ==================================================
C
      CALL EXPTOT(EXPVAL,WORK(KEXP_CORR),WORK(KEXP_TOT))
C
C     Dipole moment
C
      IF(DIPOLE) THEN
        CALL PRP_DIPOLE(EXPVAL,WORK(KEXP_CORR),.TRUE.,WORK,KFREE,LFREE)
      ENDIF
C
C     Parity violation: Chirality
C
      IF(PVC) THEN
        CALL PRP_PVC(WORK(KEXP_TOT),WORK,KFREE,LFREE)
      ENDIF
C
C     Electric fields at individual nuclei
C
      IF(EFN) THEN
        CALL PRP_EFN(EXPVAL,WORK(KEXP_CORR),.TRUE.,WORK,KFREE,LFREE)
      ENDIF
C
C     Electric field gradients
C
      IF(EFG) THEN
        CALL PRP_EFG(EXPVAL,WORK(KEXP_CORR),.TRUE.,WORK,KFREE,LFREE)
      ENDIF
C
C     non-traceless EFT's (electric field third derivatives)
C
      IF(EFT) THEN
        CALL PRP_EFT(EXPVAL,WORK(KEXP_CORR),.FALSE.,.TRUE.,WORK(KFREE),
     &               LFREE)
      ENDIF
C
C     Traceless EFT's (electric field third derivatives)
C
      IF(EFT) THEN
        CALL PRP_EFT(EXPVAL,WORK(KEXP_CORR),.TRUE.,.TRUE.,WORK(KFREE),
     &               LFREE)
      ENDIF
C
C     Nuclear quadrupole coupling
C
      IF(NQCC) THEN
        CALL PRP_NQCC(EXPVAL,WORK(KEXP_CORR),.TRUE.,WORK,KFREE,LFREE)
      ENDIF
C
C     NMR diamagnetic spin-orbit
C
      IF(DSO) THEN
        CALL PRP_DSO(WORK(KEXP_TOT),WORK,KFREE,LFREE)
      ENDIF
C
C     Memory deallocation
C     -------------------
C
      CALL MEMREL('PRPOUT_CORR',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
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
C     Correlation contribution to dipole moment
C     =========================================
C
      IF (CORR) THEN
         DO I = 1,3
           DIPMCORR(I) = EXP_CORR(IPDIPOLE(I))
         ENDDO
      ENDIF
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(A/)')
     &     '* Dipole moment:'
      IF (CORR) THEN
         WRITE(LUPRI,'(12X,3(A12,11X))')
     &        'Hartree-Fock','correlation ','Total       '
         WRITE(LUPRI,'(12X,3(A12,11X))')
     &        'contribution','contribution','MP2 value   '
      ELSE
         WRITE(LUPRI,'(10X,3(A12,13X))')
     &      'Electronic  ','Nuclear     ','Total       '
         WRITE(LUPRI,'(10X,3(A12,13X))')
     &      'contribution','contribution','contribution'
      ENDIF
      CALL PRSYMB(LUPRI,'-',76,0)
      IOFF = ICHAR('w')
      IF (CORR) THEN
         DO I = 1,3
           WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &          CHAR(IOFF+I),
     &          (DEBYE*(DIPME(I)+DIPMN(I))),' Debye',
     &          (DEBYE*DIPMCORR(I)),' Debye',
     &          (DEBYE*(DIPME(I)+DIPMN(I)+DIPMCORR(I))),' Debye'
         ENDDO
      ELSE
        DO I = 1,3
           WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &          CHAR(IOFF+I),(DEBYE*DIPME(I)),' Debye',
     &                     (DEBYE*DIPMN(I)),' Debye',
     &          (DEBYE*(DIPME(I)+DIPMN(I))),' Debye'
         ENDDO
      ENDIF
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
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
C     correlation contribution to electric field gradient
C     ===================================================
C
         IF (CORR) THEN
            CALL CORR_EFG(0,EXP_CORR,DOATOM,TEMP)
         ENDIF
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
      NATOM = 0
      DO 100 IATOM = 1, NUCIND
         DO 110 ISYM = 0,MAXOPR
            IF (IBTAND(ISTBNU(IATOM),ISYM) .EQ. 0) THEN
               NATOM = NATOM + 1
               DO I=1,3
                  DO J=1,3
                     EFGTOT(I,J,NATOM)=
     &                  EFGNUC(I,J,NATOM)+EFGEL(I,J,NATOM)
                     IF (CORR) THEN
                        EFG_METHOD(I,J,NATOM)=
     &                     EFGTOT(I,J,NATOM)+EFGCORR(I,J,NATOM)
                     ENDIF
                  ENDDO
               ENDDO
C
C              Diagonalize to get principal components
C
               IJ = 1
               DO I = 1, 3
                  DO J = 1, I
                     IF (CORR) THEN
                        PVAL(IJ) = EFG_METHOD(I,J,IATOM)
                     ELSE
                        PVAL(IJ) = EFGTOT(I,J,IATOM)
                     ENDIF
                     IJ = IJ + 1
                  ENDDO
               ENDDO
               CALL JACO(PVAL,AXIS,3,3,3,TMP,ITMP)
               IF (CORR) THEN
                  EFG_METHOD(1,1,IATOM) = PVAL(1)
                  EFG_METHOD(2,2,IATOM) = PVAL(3)
                  EFG_METHOD(3,3,IATOM) = PVAL(6)
               ELSE
                  EFGTOT(1,1,IATOM) = PVAL(1)
                  EFGTOT(2,2,IATOM) = PVAL(3)
                  EFGTOT(3,3,IATOM) = PVAL(6)
               ENDIF
C
C              Reorder to calculate eta:
C              |xx|<|yy|<|zz|
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
                     FACT=CONV*QMOM
C
C                    Change to GHz when abs values bigger then 1.0D4
C
                     IF ((ABS(EFGTOT(1,1,IATOM))*FACT.GT.1.0D4).OR.
     &                   (ABS(EFGTOT(2,2,IATOM))*FACT.GT.1.0D4).OR.
     &                   (ABS(EFGTOT(3,3,IATOM))*FACT.GT.1.0D4)) THEN
                        FACT   = FACT * 1.0D-3
                        HZUNIT = ' GHz'
                     ELSE
                        HZUNIT = ' MHz'
                     ENDIF
                     SPIN = DISOTP(NUCCHA,ISO,'SPIN')
                     MASS = NINT(DISOTP(NUCCHA,ISO,'A'))
                     WRITE (LUPRI,9000)
     &                  'Nucleus: ',NAMDEP(IATOM),
     &                  'Mass: ',MASS,'Spin: ',
     &                  SPIN,'Quadrupole moment: ',QMOM
                     DO I=1,3
                        IF (CORR) THEN
                           WRITE(LUPRI,'(2X,A3,6X,3(3X,F16.10,A))')
     &                        'X'//CHAR(IOFF+I)//CHAR(IOFF+I),
     &                        EFG_METHOD(I,I,IATOM)*FACT,HZUNIT
                        ELSE
                           WRITE(LUPRI,'(2X,A3,6X,3(3X,F16.10,A))')
     &                        'X'//CHAR(IOFF+I)//CHAR(IOFF+I),
     &                        EFGTOT(I,I,IATOM)*FACT,HZUNIT
                        ENDIF
                     ENDDO
                     WRITE (LUPRI,'(2X,A3,9X,F16.8)')
     &                  'Eta', ETA
                  ENDIF
               ENDDO
            END IF
 110     CONTINUE
 100  CONTINUE
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(/A,F12.8,A)') '  1 a.u. = ',CONV,' MHz'
      RETURN
 9000 FORMAT (/2X,A9,A4,3X,A6,I3,3X,A6,F4.1,3X,A19,F12.7,/)
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
      NATOM = 1
      DO 110 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 110
         DO 310 IREPC = 0, MAXREP
            IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
               DO 220 I = 1, 3
                  DO 230 J = I, 3
                     EFGEL(I,J,NATOM) = TEMP(I,J,IATOM)*
     &                         PT(IBTAND(IBTXOR(ISYMAX(I,1),
     &                            ISYMAX(J,1)),IREPC))
                     IF (I .NE. J) EFGEL(J,I,NATOM) = EFGEL(I,J,NATOM)
 230              CONTINUE
 220           CONTINUE
               IF (IPRINT .GE. 2) THEN
                  WRITE (LUPRI ,'(/2X,A,A/)')
     &               'Electronic EFG Tensor of nucleus ',NAMDEP(NATOM)
                  CALL OUTPUT(EFGEL(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
               END IF
               NATOM = NATOM + 1
            END IF
 310     CONTINUE
 110  CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck corr_efg */
      SUBROUTINE CORR_EFG(IPRINT,EXPVAL_CORR,DOATOM,TEMP)
C*****************************************************************************
C
C     Calculates correlation contribution to electric field gradients
C     Based on EL_EFG
C
C     Written by Joost van Stralen  - June 11 2002 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      DIMENSION EXPVAL_CORR(*), TEMP(3,3,MXCENT)
      LOGICAL DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
      CALL DZERO(EFGCORR,9*MXCENT)
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
     &               TEMP(ICOOR1,ICOOR2,IATOM) +EXPVAL_CORR(IPEFG(ITYP))
     &                                             /NUCDEG(IATOM)
                  END IF
 300           CONTINUE
 210        CONTINUE
 200     CONTINUE
 100  CONTINUE
      NATOM = 1
      DO 110 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 110
         DO 310 IREPC = 0, MAXREP
            IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
               DO 220 I = 1, 3
                  DO 230 J = I, 3
                     EFGCORR(I,J,NATOM) = TEMP(I,J,IATOM)*
     &                         PT(IBTAND(IBTXOR(ISYMAX(I,1),
     &                            ISYMAX(J,1)),IREPC))
                     IF (I .NE. J) EFGCORR(J,I,NATOM)=EFGCORR(I,J,NATOM)
 230              CONTINUE
 220           CONTINUE
               IF (IPRINT .GE. 2) THEN
                  WRITE (LUPRI ,'(/2X,A,A/)')
     &               'MP2 corr. EFG Tensor of nucleus ',NAMDEP(NATOM)
                  CALL OUTPUT(EFGCORR(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
               END IF
               NATOM = NATOM + 1
            END IF
 310     CONTINUE
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
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
      IF (CORR) THEN
         CALL CORR_EFG(IPREXP,EXP_CORR,DOATOM,TEMP)
      END IF
C
C     print section
      IOFF = ICHAR('w')
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(12X,A/)')
     &     'Individual (non-zero) components'
      IF (CORR) THEN
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'Hartree-Fock','correlation ','Total       '
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'contribution','contribution','MP2 value   '
      ELSE
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      END IF
      CALL PRSYMB(LUPRI,'-',72,2)
      KATOM = 0
      DO JATOM = 1, NUCIND
      NDEG = NUCDEG(JATOM)
      IF (DOATOM(JATOM)) THEN
        DO IDEG = 1,NDEG
          IATOM = KATOM + IDEG
          WRITE (LUPRI ,'(/2X,A9,A6,/)')'Nucleus: ',NAMDEP(IATOM)
          DO I = 1,3
            DO J = 1,3
              IF (CORR) THEN
              EFGTOT(I,J,IATOM)=EFGNUC(I,J,IATOM)+EFGEL(I,J,IATOM)
              EFG_METHOD(I,J,IATOM)=EFGTOT(I,J,IATOM)+EFGCORR(I,J,IATOM)
              IF (ABS(EFG_METHOD(I,J,IATOM)).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A3,3(3X,F16.10,A))')
     &              'q'//CHAR(IOFF+I)//CHAR(IOFF+J),EFGTOT(I,J
     &              ,IATOM),' au ',EFGCORR(I,J,IATOM),' au '
     &              ,EFG_METHOD(I,J,IATOM),' au '
              ELSE
                EFGTOT(I,J,IATOM)=EFGNUC(I,J,IATOM)+EFGEL(I,J,IATOM)
                IF (ABS(EFGTOT(I,J,IATOM)).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A3,3(3X,F16.10,A))')
     &              'q'//CHAR(IOFF+I)//CHAR(IOFF+J),EFGEL(I,J
     &              ,IATOM),' au ',EFGNUC(I,J,IATOM),' au '
     &              ,EFGTOT(I,J,IATOM),' au '
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      KATOM = KATOM + NDEG
      ENDDO
C
C           diagonalize total efg tensor
C
      WRITE(LUPRI,'(/,A)') '  '
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(12X,A)')
     &     'Total contribution to principal components '
      CALL PRSYMB(LUPRI,'-',72,2)
      KATOM = 0
      DO JATOM = 1, NUCIND
      NDEG = NUCDEG(JATOM)
      IF (DOATOM(JATOM)) THEN
        DO IDEG = 1,NDEG
          IATOM = KATOM + IDEG
          CALL DUNIT(AXIS(1,1,IATOM),3)
          IJ = 1
          DO I = 1, 3
            DO J = 1, I
            IF (CORR) THEN
              PVAL(IJ) = EFG_METHOD(I,J,IATOM)
            ELSE
              PVAL(IJ) = EFGTOT(I,J,IATOM)
            END IF
              IJ = IJ + 1
            ENDDO
          ENDDO
          CALL JACO(PVAL,AXIS(1,1,IATOM),3,3,3,TMP,ITMP)
          IF (CORR) THEN
            EFG_METHOD(1,1,IATOM) = PVAL(1)
            EFG_METHOD(2,2,IATOM) = PVAL(3)
            EFG_METHOD(3,3,IATOM) = PVAL(6)
          ELSE
            EFGTOT(1,1,IATOM) = PVAL(1)
            EFGTOT(2,2,IATOM) = PVAL(3)
            EFGTOT(3,3,IATOM) = PVAL(6)
          ENDIF
          WRITE (LUPRI ,'(/2X,A9,A6,/)')'Nucleus: ',NAMDEP(IATOM)
          DO I=1,3
            IF (CORR) THEN
              WRITE(LUPRI,'(2X,A3,3(3X,F16.10,A))')
     &                 'q'//CHAR(IOFF+I)//CHAR(IOFF+I)
     &                 ,EFG_METHOD(I,I,IATOM),' au '
            ELSE
              WRITE(LUPRI,'(2X,A3,3(3X,F16.10,A))')
     &                 'q'//CHAR(IOFF+I)//CHAR(IOFF+I)
     &                 ,EFGTOT(I,I,IATOM),' au '
            END IF
          ENDDO
C
C           print section principle axis
C
          IF (IPRINT .GE. 4) THEN
            WRITE(LUPRI,'(A/)')
     &         '  Principal axis to total contribution'
            WRITE (LUPRI ,'(/2X,A9,A6,/)')'Nucleus: ',NAMDEP(IATOM)
            CALL OUTPUT(AXIS(1,1,IATOM),1,3,1,3,3,3,1,LUPRI)
          END IF
        ENDDO
      ENDIF
      KATOM = KATOM + NDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(/A,1PE12.5,A,/)')
     &      '  1 a.u. = ',9.71736E21,' V*m**-2'
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
      LOGICAL   CORR 		!add correlation contribution?
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
         CALL SETATM(WORK(KATOM),NATOM,IDUM)
      ELSE
         CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
      CALL DZERO(TMPEFT_CORR,15*NUCIND)
      CALL DZERO(EFTCORR    ,15*NUCDEP)
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
                  IF (CORR) THEN
                  TMPEFT_CORR(IELEM,IATOM) =
     &               TMPEFT_CORR(IELEM,IATOM) +
     &               EXP_CORR(IPE) / NUCDEG(IATOM)
                  ENDIF
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
                  IF (CORR) THEN
                  EFTCORR(IELEM,NATOM) =
     &               TMPEFT_CORR(IELEM,IATOM) * PT(IBTAND(ISYMIJ,IREPC))
                  END IF
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
      IF (CORR) THEN
         WRITE(LUPRI,'(14X,3(A12,11X))')
     &        'Hartee-Fock ','correlation  ','Total       '
         WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'contribution','contribution','MP2 value   '
      ELSE
         WRITE(LUPRI,'(14X,3(A12,11X))')
     &        'Electronic  ','Nuclear     ','Total       '
         WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      ENDIF
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
C     EFT correlated method
               IF (CORR) THEN
               EFTHF = EFTNUC(IELEM,IATOM)+EFTEL(IELEM,IATOM)
               EFT_METHOD = EFTHF + EFTCORR(IELEM,IATOM)
               IF (ABS(EFT_METHOD).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A5,3(3X,1P,D16.9,A))')
     &                 'h'//CHRXYZ(-ICOOR1)//CHRXYZ(-ICOOR2)
     &                    //CHRXYZ(-ICOOR3)//CHRXYZ(-ICOOR4),
     &                 EFTHF                  ,' au ',
     &                 EFTCORR(IELEM,IATOM)   ,' au ',
     &                 EFT_METHOD             ,' au '
C     EFT Hartee-Fock
               ELSE
               EFTTOT = EFTNUC(IELEM,IATOM)+EFTEL(IELEM,IATOM)
               IF (ABS(EFTTOT).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A5,3(3X,1P,D16.9,A))')
     &                 'h'//CHRXYZ(-ICOOR1)//CHRXYZ(-ICOOR2)
     &                    //CHRXYZ(-ICOOR3)//CHRXYZ(-ICOOR4),
     &                 EFTEL (IELEM,IATOM),' au ',
     &                 EFTNUC(IELEM,IATOM),' au ',
     &                 EFTTOT             ,' au '
               END IF
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
               IF (CORR) THEN
               DTRCORR(I)  = EFTCORR (ITR(1,I),IATOM) +
     &                       EFTCORR (ITR(2,I),IATOM) +
     &                       EFTCORR (ITR(3,I),IATOM)
               END IF
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
         IF (CORR) THEN
            EFTCORR( 1,IATOM) = EFTCORR( 1,IATOM) -
     &         D27/D35*DTRCORR(1) +  D3/D35*DTRCORR(4) +
     &         D3/D35*DTRCORR(6)
            EFTCORR( 2,IATOM) = EFTCORR( 2,IATOM) - D3/D7*DTRCORR(2)
            EFTCORR( 3,IATOM) = EFTCORR( 3,IATOM) - D3/D7*DTRCORR(3)
            EFTCORR( 4,IATOM) = EFTCORR( 4,IATOM) +
     &          (-D4/D35)*DTRCORR(1) -  D4/D35*DTRCORR(4) +
     &           D1/D35*DTRCORR(6)
            EFTCORR( 5,IATOM) = EFTCORR( 5,IATOM) - D1/D7*DTRCORR(5)
            EFTCORR( 6,IATOM) = EFTCORR( 6,IATOM) +
     &          (-D4/D35)*DTRCORR(1) +  D1/D35*DTRCORR(4) -
     &           D4/D35*DTRCORR(6)
            EFTCORR( 7,IATOM) = EFTCORR( 7,IATOM) - D3/D7*DTRCORR(2)
            EFTCORR( 8,IATOM) = EFTCORR( 8,IATOM) - D1/D7*DTRCORR(3)
            EFTCORR( 9,IATOM) = EFTCORR( 9,IATOM) - D1/D7*DTRCORR(2)
            EFTCORR(10,IATOM) = EFTCORR(10,IATOM) - D3/D7*DTRCORR(3)
            EFTCORR(11,IATOM) = EFTCORR(11,IATOM) +
     &           D3/D35*DTRCORR(1) - D27/D35*DTRCORR(4) +
     &           D3/D35*DTRCORR(6)
            EFTCORR(12,IATOM) = EFTCORR(12,IATOM) - D3/D7*DTRCORR(5)
            EFTCORR(13,IATOM) = EFTCORR(13,IATOM) +
     &           D1/D35*DTRCORR(1) -  D4/D35*DTRCORR(4) -
     &           D4/D35*DTRCORR(6)
            EFTCORR(14,IATOM) = EFTCORR(14,IATOM) - D3/D7*DTRCORR(5)
            EFTCORR(15,IATOM) = EFTCORR(15,IATOM) +
     &           D3/D35*DTRCORR(1) +  D3/D35*DTRCORR(4) -
     &           D27/D35*DTRCORR(6)
         END IF
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
      IF (CORR) THEN
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'Hartree-Fock','correlation ','Total       '
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'contribution','contribution','MP2 value   '
      ELSE
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(14X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      END IF
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
C     EFT correlated method
               IF (CORR) THEN
               EFTHF = EFTNUC(IELEM,IATOM)+EFTEL(IELEM,IATOM)
               EFT_METHOD = EFTHF + EFTCORR(IELEM,IATOM)
               IF (ABS(EFT_METHOD).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A5,3(3X,F16.10,A))')
     &                 'h'//CHRXYZ(-ICOOR1)//CHRXYZ(-ICOOR2)
     &                    //CHRXYZ(-ICOOR3)//CHRXYZ(-ICOOR4),
     &                 EFTHF               ,' au ',
     &                 EFTCORR(IELEM,IATOM),' au ',
     &                 EFT_METHOD          ,' au '
               ELSE
C     EFT Hartree-Fock
               EFTTOT = EFTNUC(IELEM,IATOM)+EFTEL(IELEM,IATOM)
               IF (ABS(EFTTOT).GT.1.0D-12)
     &            WRITE(LUPRI,'(2X,A5,3(3X,F16.10,A))')
     &                 'h'//CHRXYZ(-ICOOR1)//CHRXYZ(-ICOOR2)
     &                    //CHRXYZ(-ICOOR3)//CHRXYZ(-ICOOR4),
     &                 EFTEL (IELEM,IATOM),' au ',
     &                 EFTNUC(IELEM,IATOM),' au ',
     &                 EFTTOT             ,' au '
               END IF
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
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
         IF (.NOT. DOATOM(IATOM)) GOTO 100
         DO IREP = 0,MAXREP
            IF (IBTAND(ISTBNU(IATOM),IREP) .EQ. 0) THEN
               ITYP = ITYP + 1
c              write(6,*) 'jth iatom',
c    &            iatom,ippvc(ityp),ityp,
c    &            expval(ippvc(ityp)),nucdeg(iatom)
               TEMP(IATOM) = TEMP(IATOM) +
     &            EXPVAL(IPPVC(ITYP)) / NUCDEG(IATOM)
            END IF
         END DO
 100  CONTINUE
C
C     For each atom calculate exp. value
C
      DTOT = D0
      NATOM = 0
      DO 110 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 110
         DO IREP = 0 ,MAXREP
            IF (IBTAND(ISTBNU(IATOM),IREP).EQ.0) THEN
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
            END IF
         END DO
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
      NATOM = 0
      DO 120 IATOM = 1, NUCIND
         IF (.NOT. DOATOM(IATOM)) GOTO 120
         DO IREP = 0 ,MAXREP
            IF (IBTAND(ISTBNU(IATOM),IREP).EQ.0) THEN
               NATOM =  NATOM + 1
               WRITE(LUPRI,9010)
     &            NAMDEP(NATOM),IS(NATOM),DMPV(NATOM),QW(NATOM),
     &            DMPV(NATOM)*QW(NATOM)*DFAC
            END IF
         END DO
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
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dummy.h"
#include "codata.h"
      PARAMETER (PRMASS = 1.007276470D0, D3 = 3.0D0, D0 = 0.0D0)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PRMASS*PRMASS)
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
C          WRITE(LUPRI,'
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
        CALL AROUND('DIRAC - Final spin-spin-couplings')
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
C
      ENDDO
      CALL QEXIT('PRP_SPNSP1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck spiana */
      SUBROUTINE SPIANA(PRIBAS,SPNPAR,IATOM1,IATOM2,
     &                  WORK,KFREE,LFREE)
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dummy.h"
#include "codata.h"
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
      CALL AROUND('DIRAC - DSO contribution to spin-spin-couplings')
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
C
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "gfac.h"
#include "nuclei.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbesr.h"
#include "dcbxlr.h"
#include "chrxyz.h"
      LOGICAL DOATOM(NUCIND)
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"
C
C
      NSCOOR = 3*NUCDEP
      CALL IZERO(IPSHIELD,3+MXCOOR)
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
C
C         Cartesian directions
C
          DO 300 ICOOR1 = 1, 3
            ISCOR1 = IPTCNT(3*(IATOM1 - 1) + ICOOR1,IREPO,2)
            IF (ISCOR1 .GT. 0) THEN
              PNAME    = 'NUCSPI'//CHRNOS(ISCOR1/10)
     &                         //CHRNOS(MOD(ISCOR1,10))
     &                //NAMN(IATOM1)//REP(IREPO)//CHRXYZ(ICOOR1)
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
              PLABEL(1)= 'NEF '//CHRNOS(IS1/10)
     &                     //CHRNOS(MOD(IS1,10))//'  '
              PLABEL(2)= 'NEF '//CHRNOS(IS2/10)
     &                     //CHRNOS(MOD(IS2,10))//'  '
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
            ENDIF
  300     CONTINUE
        END IF
  200   CONTINUE
  100 CONTINUE
C
C     Define right hand operators
C     (PFAC(1:2) for Dirac mag. field, PFAC(3) for g_1 QED correction)
C
      PFAC(1)  = -0.5D0*CVAL
      PFAC(2)  = -0.5D0*CVAL
      PFAC(3)  =  0.5D0*G1FAC
C X component of B-field
      PNAME     = 'MAGFLDX1'
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
      PNAME     = 'MAGFLDY2'
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
      PNAME     = 'MAGFLDZ3'
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
C
C     Print section for NMR shielding:
C
      IF (SHIELD) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(A)') ' NMR shieldings:'
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE (LUPRI,'(A,E9.2)')
     &   ' Shieldings printed for atoms with '//
     &   'abundance greater than :', ABUND
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
C     Parity violating contributions to shielding
C
      IF(PVCSHI) THEN
        WRITE(LUPRI,'(A)') 'Parity violating contributions..'
        ITYP = 0
C.......Xalpha
        IOP(1)  = 2
        PFAC(1) = CVAL
C.......Yalpha
        IOP(2)  = 3
        PFAC(2) = CVAL
C.......Zalpha
        IOP(3) = 4
        PFAC(2) = CVAL
C
        NPCOMP    = 1
        DO 110 IATOM = 1, NUCIND
        IF (DOATOM(IATOM)) THEN
          DO 120 IREP = 0, MAXREP
          IF (IBTAND(IREP,ISTBNU(IATOM)).EQ.0) THEN
            ITYP = ITYP + 1
            PLABEL(1) = 'PVC'//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                        //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                        CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
            PNAME = 'PVCa '//NAMDEP(NUCPRE(IATOM)+1)(1:3)
     &                    //CHRNOS(IPTNUC(IATOM,IREP)/10)//
     &                    CHRNOS(MOD(IPTNUC(IATOM,IREP),10))
            DO I = 1,3
              IPTYP=IOP(I)
              PNAME(5:5) = CHRNOS(I)
              CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &             INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
              CALL OP1IND(PNAME(1:6),INDOP1,LLRAPU,NLRAPT,
     &                      INDXPR,MAXLLR)
              CALL WRIXPR(INDOP1,INDXPR)
C              CALL OP1IND('PVC   '
C     &                 ,IPPVC(ITYP),LEXPP,NEXPP,INDXPR,MAXEXP)
            ENDDO
          END IF
  120     CONTINUE
        END IF
  110   CONTINUE
      ENDIF
C
C     B operators
C     ===========
C
        CALL PRSYMB(LUPRI,'-',75,0)
        WRITE(LUPRI,'(A)') ' ** B operators for shielding **'
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
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck DEF_MSUSCP */
      SUBROUTINE DEF_MSUSCP(WORK,LWORK)   
C***********************************************************************
C
C     Define operators for magnetizabilities
C
C     T. Saue Jan 20 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbxlr.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(2)
C
      PFAC(1)  = -0.5D0*CVAL
      PFAC(2)  = -0.5D0*CVAL
      NPCOMP = 2
C X component of B-field
      PNAME     = 'MAGFLDX1'
      IPTYP     = 5
      PLABEL(1) = 'ZDIPLEN'
      PLABEL(2) = 'YDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('XMAGNA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      CALL OP1IND('XMAGNB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
C Y component of B-field
      PNAME     = 'MAGFLDY2'
      IPTYP     = 6
      PLABEL(1) = 'XDIPLEN'
      PLABEL(2) = 'ZDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('YMAGNA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      CALL OP1IND('YMAGNB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
C Z component of B-field
      PNAME     = 'MAGFLDZ3'
      IPTYP     = 7
      PLABEL(1) = 'YDIPLEN'
      PLABEL(2) = 'XDIPLEN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      CALL OP1IND('ZMAGNA',INDOP1,LLRAPU,NLRAPT,INDXPR,MAXLLR)
      CALL OP1IND('ZMAGNB',INDOP1,LLRBPU,NLRBPT,INDXPR,MAXLLR)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEF_NSTDIA */
      SUBROUTINE DEF_NSTDIA(PNAME,PLABEL,PFAC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Define diamagnetic contribution to shieldings
C     Written by T. Enevoldsen Oct. 1997
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK),PFAC(3),IOP(3)
#include "chrnos.h"
#include "ibtfun.h"
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
                        PNAME = 'NSCO '//NAMN(IATOM1)//' '
     &                       //CHRNOS(ISCOR1/10)//CHRNOS(MOD(ISCOR1,10))
     &                       //CHRXYZ(ICOOR2)
                        PLABEL(1) = CHRNOS(ISCOR1/10)//
     &                       CHRNOS(MOD(ISCOR1,10))//'NSCO'//
     &                       ' '//CHRXYZ(ICOOR2)
                        CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &                       INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
                        IF(IPRINT.GE.1) THEN
                           CALL WRIXPR(NNSTIN,INDXPR)
                        ENDIF
                        CALL OP1IND('XNSGO',IPIND,LEXPP,NEXPP,
     &                       INDXPR,MAXEXP)
                        IPNSTDIA(NNSTIN) = IPIND
 600                 CONTINUE
                  END IF
 500           CONTINUE
            END IF
 200     CONTINUE
 100  CONTINUE
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
      CALL MEMGET('REAL',KPMAT,NCMAT,WORK,KFREE,LFREE)
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
      CALL MEMGET('REAL',KTMPMAT,NCMAT,WORK,KFREE,LFREE)
C
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
      CALL PRP_SHIELD1(EXPVAL,ABLRF,WORK(KTMPMAT),WORK(KTMAT),
     &     WORK(KDMAT),WORK(KPMAT),WORK(KTISO),WORK(KDISO),WORK(KPISO),
     &     WORK(KANIS),WORK(KASYM),WORK(KSKEW),WORK(KSPAN),WORK(KIBUF),
     &     WORK(KATOM),WORK(KCSTRA),WORK(KSCTRA),WORK,KFREE,LFREE)
      CALL MEMREL('PRP_SHIELD',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_SHIELD1 */
      SUBROUTINE PRP_SHIELD1(EXPVAL,ABLRF,TMPMAT,TMAT,DMAT,PMAT,TISO
     &     ,DISO,PISO,ANIS,ASYM,SKEW,SPAN,IBUFA
     &     ,DOATOM,CSTRA,SCTRA,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Output routine for NMR shieldings
C     C
C     Written by T.Enevoldsen Nov. 7 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "dcbxlr.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION EXPVAL(*),ABLRF(NBFREQ,NLRAPT,NLRBPT),
     &     TMPMAT(3,*),TMAT(3,3,NUCDEP), DMAT(3,3,NUCDEP),
     &     PMAT(3,3,NUCDEP),TISO(NUCDEP),DISO(NUCDEP),PISO(NUCDEP),
     &     ANIS(NUCDEP),ASYM(NUCDEP),SKEW(NUCDEP),SPAN(NUCDEP),
     &     IBUFA(3*NUCDEP),IBUFB(3),CSTRA(NUCDEP*NUCDEP*9),
     &     SCTRA(NUCDEP*NUCDEP*9),WORK(*)
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
      CALL TITLER('DIRAC  - NMR SHIELDINGS','*',116)
C
      KFRSAV = KFREE
      NSCOOR = 3*NUCDEP
      NCMAT = 9*NUCDEP
      CALL DZERO(TMPMAT,NCMAT)
      CALL DZERO(DMAT,NCMAT)
C
C     Extract Diamagnetic part if e-p rotations are not included
C
      IF(XLR_SKIPEP) THEN
         CALL DZERO(DMAT(1,1,1),NCMAT)
         DO I = 1,NNSTIN
            IOP = IPNSTDIA(I)
            IPL = IPRPLBL(1,LEXPP(IOP))
            READ(PRPLBL(IPL)(:2),'(I2)') IA
            IB = MOD(I+2,3) + 1
            TMPMAT(IPTAX(IB,2),IA) = EXPVAL(IOP)
         ENDDO
         IF (IPRPRP.GE.4) THEN
            CALL HEADER('Diamagnetic part of shielding'//
     &           ' in symmetry coordinates (ppm)',1)
            CALL FCPRI(TMPMAT,'SIGMAS',CSTRA,SCTRA)
         ENDIF
         CALL TRADIP(TMPMAT,DMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
      END IF
C
C     *** Extract shieldings ***
C
      NA = 0
      NB = 0
      DO I = 1,NSCOOR
         IF(IPSHIELD(I).GT.0) THEN
            NA = NA + 1
            IBUFA(NA) = I
         ENDIF
      END DO
C
      DO I = 1, 3
         IF(IPSHIELD(I+MXCOOR).GT.0) THEN
            NB = NB + 1
            IBUFB(NB) = I
         ENDIF
      END DO

C
C     Fill shielding matrices
C
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

         CALL DZERO(TMPMAT,NCMAT)
         DO IB = 1,NB
            READ (PRPNAM(LLRBPU(IBUFB(IB))),'(7X,I1)') JB
            KB = IPSHIELD(JB+MXCOOR)
            DO IA = 1,NA
               JA = IBUFA(IA)
               KA = IPSHIELD(JA)
               TMPMAT(IPTAX(JB,2),JA) = ABLRF(IFREQ,KA,KB)
            ENDDO
         ENDDO
#if defined (XH4)
         TMPMAT(2,8) = TMPMAT(1,4)
         TMPMAT(2,9) = TMPMAT(1,6)
         TMPMAT(2,10) = TMPMAT(1,7)
         TMPMAT(2,11) = TMPMAT(1,5)
         TMPMAT(3,12) = TMPMAT(1,4)
         TMPMAT(3,13) = TMPMAT(1,6)
         TMPMAT(3,15) = TMPMAT(1,7)
         TMPMAT(3,14) = TMPMAT(1,5)
#endif
         IF (IPRPRP.GE.4) THEN
            CALL HEADER('Response part of shielding'//
     &           ' in symmetry coordinates (ppm)',1)
            CALL FCPRI (TMPMAT,'SIGMAS',CSTRA,SCTRA)
         ENDIF
C
C     Transform to non-symmetry basis
C
         CALL DZERO(PMAT(1,1,1),NCMAT)
         CALL TRADIP(TMPMAT,PMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C
C     Analyze results
C
         CALL DCOPY(9*NUCDEP,PMAT,1,TMAT,1)
         IF (XLR_SKIPEP) CALL DAXPY(9*NUCDEP,D1,DMAT,1,TMAT,1)
         IATOM = 0
         DO 100 I = 1, NUCIND
            DO 200 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(I)) THEN
                     CALL SHIANA(TMAT(1,1,IATOM),DMAT(1,1,IATOM),
     &                 PMAT(1,1,IATOM),TISO(IATOM),DISO(IATOM),
     &                 PISO(IATOM),ANIS(IATOM),ASYM(IATOM),SKEW(IATOM),
     &                 SPAN(IATOM),NAMDEP(IATOM),NAMDPX(3*(IATOM-1)+1),
     &                 IPRPRP)
                  END IF
               END IF
 200        CONTINUE
 100     CONTINUE
C
         CALL AROUND('Summary of chemical shieldings')
         WRITE (LUPRI,'(1X,A)') 'Definitions from J.Mason, Solid '//
     &        'state Nuc.Magn.Res. 2 (1993), 285'
         WRITE (LUPRI,'()')
         WRITE (LUPRI,'(A,/,A,/,A,/,A)')
     &        '@            isotropic shielding',
     &        '@        ----------------------------',
     &        '@atom       total       dia      para      skew'//
     &        '      span      anis      asym',
     &        '@-----------------------------------------------------'//
     &        '-----------------------'

         IATOM = 0
         DO 300 I = 1, NUCIND
            DO 400 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(I)) THEN
                      WRITE (LUPRI,'(A,A,7F10.4)') '@',
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
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck SHIANA */
      SUBROUTINE SHIANA(TMAT,DMAT,PMAT,AVETOT,AVEDIA,AVEPAR,ANIS,ASYM,
     &                  SKEW,SPAN,NAME,NAMEX,IPRINT)
C***********************************************************************
C
C     Analyze routine for NMR shieldings
C     This routine is modelled after SHIANA in DALTON
C
C     Written by T.Enevoldsen Nov. 14 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
#include "dcbxlr.h"
      PARAMETER (FACTOR = 1.D6*ALPHA2)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      CHARACTER NAME*6, NAMEX(3)*6, SGN*3
      INTEGER Z
      DIMENSION TMAT(3,3), DMAT(3,3), PMAT(3,3), SYM(3,3), ANTI(3,3),
     &          AXES(3,3), PVAL(6), WRK(3), IWRK(3), ABSP(3)
#include "chrxyz.h"
C
      CALL DSCAL(9,FACTOR,TMAT,1)
      CALL DSCAL(9,FACTOR,DMAT,1)
      CALL DSCAL(9,FACTOR,PMAT,1)
C
      AVETOT = (TMAT(1,1) + TMAT(2,2) + TMAT(3,3))/D3
      AVEDIA = (DMAT(1,1) + DMAT(2,2) + DMAT(3,3))/D3
      AVEPAR = (PMAT(1,1) + PMAT(2,2) + PMAT(3,3))/D3
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
      IF (XLR_SKIPEP) THEN
      CALL HEADER('Diamagnetic and paramagnetic contributions (ppm):',1)
      WRITE (LUPRI,'(15X,3(A,9X),4X,3(A,9X),/)')
     &    'Bx', 'By', 'Bz', 'Bx','By','Bz'
      DO 500 ICOOR = 1, 3
         WRITE (LUPRI,'(2X,A6,3F11.4,4X,3F11.4)') NAMEX(ICOOR),
     &         (DMAT(K,ICOOR),K=1,3), (PMAT(K,ICOOR),K=1,3)
  500 CONTINUE
      WRITE (LUPRI,'(/,1X,A,F14.6,8X,A,F14.6)')
     &          ' Diamagnetic contribution:',AVEDIA,
     &          ' Paramagnetic: ',AVEPAR
      END IF
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
      CALL SETATM(WORK(KATOM),NATOM,IDUM)
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
      IF (CORR) THEN
         CALL CORR_EFN(IPREXP,EXP_CORR,DOATOM,EFNVEC(1,1,3),TEMP)
      END IF
C     
C     print section
      IOFF = ICHAR('w')
      CALL PRSYMB(LUPRI,'-',72,2)
      WRITE(LUPRI,'(12X,A/)')
     &     'Individual (non-zero) components'
      IF (CORR) THEN
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'Hartree-Fock','correlation ','Total       '
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'contribution','contribution','MP2 value   '
      ELSE
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(12X,3(A12,11X))')
     &     'contribution','contribution','contribution'
      END IF
      CALL PRSYMB(LUPRI,'-',72,2)
      KATOM = 0
      DO JATOM = 1, NUCIND
      NDEG = NUCDEG(JATOM)
      IF (DOATOM(JATOM)) THEN
        DO IDEG = 1,NDEG
          IATOM = KATOM + IDEG
          WRITE (LUPRI ,'(/2X,A9,A6,/)')'Nucleus: ',NAMDEP(IATOM)
          DO I = 1,3
            IF (CORR) THEN
              EFNTOT1=EFNVEC(I,IATOM,1)+EFNVEC(I,IATOM,2)
              EFNTOT2 = EFNTOT1 + EFNVEC(I,IATOM,3)
              WRITE(LUPRI,'(2X,A1,3(3X,F16.10,A))') CHAR(IOFF+I),
     &             EFNTOT1,' au ',EFNVEC(I,IATOM,3),' au ',
     &             EFNTOT2,' au '
            ELSE
              EFNTOT=EFNVEC(I,IATOM,1)+EFNVEC(I,IATOM,2)
CTROND                 IF (ABS(EFGTOT(I,J,IATOM)).GT.1.0D-12)
              WRITE(LUPRI,'(2X,A1,3(3X,F16.10,A))') CHAR(IOFF+I),
     &             EFNVEC(I,IATOM,1),' au ',EFNVEC(I,IATOM,2),' au ',
     &             EFNTOT,' au '
CTROND                 ENDIF
            ENDIF
          ENDDO
        END DO
      ENDIF
      KATOM = KATOM + NDEG
      END DO
      CALL PRSYMB(LUPRI,'-',72,2)
CTROND      WRITE(LUPRI,'(/A,1PE12.5,A,/)') 
CTROND     &      '  1 a.u. = ',9.71736E21,' V*m**-2' 
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
      NATOM = 1
      DO 110 IATOM = 1, NUCIND
        IF (.NOT. DOATOM(IATOM)) GOTO 110
        DO 310 IREPC = 0, MAXREP
          IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
            DO 220 I = 1, 3
              EFNEL(I,NATOM) = TEMP(I,IATOM)*
     &                         PT(IBTAND(ISYMAX(I,1),IREPC))
 220        CONTINUE
            IF (IPRINT .GE. 2) THEN
              WRITE (LUPRI ,'(/2X,A,A/)')
     &         'Electronic electronic field at nucleus ',NAMDEP(NATOM)
              IOFF = ICHAR('w')
              DO I = 1,3
                WRITE(LUPRI,'(4X,A1,A,E20.10)') CHAR(IOFF+I),
     &              '-component: ',EFNEL(I,NATOM)
              ENDDO
            END IF
            NATOM = NATOM + 1
          END IF
 310    CONTINUE
 110  CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck corr_efn */
      SUBROUTINE CORR_EFN(IPRINT,EXPVAL_CORR,DOATOM,EFNCORR,TEMP)
C*****************************************************************************
C           
C     Calculates correlation contribution to electric fields
C     Based on EL_EFN
C             
C     Written by Joost van Stralen - June 9 2004
C               
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C           
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbexp.h"
      DIMENSION EXPVAL_CORR(*), TEMP(3,NUCIND),EFNCORR(3,NUCDEP)
      LOGICAL DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
      CALL DZERO(EFNCORR,3*NUCDEP)
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
     &                        + EXPVAL_CORR(IPEFN(ITYP))/NUCDEG(IATOM)
          END IF
 200      CONTINUE
  50    CONTINUE
 100  CONTINUE
      NATOM = 1
      DO 110 IATOM = 1, NUCIND
        IF (.NOT. DOATOM(IATOM)) GOTO 110
        DO 310 IREPC = 0, MAXREP
          IF (IBTAND(ISTBNU(IATOM),IREPC).EQ.0) THEN
            DO 220 I = 1, 3
              EFNCORR(I,NATOM) = TEMP(I,IATOM)*
     &                         PT(IBTAND(ISYMAX(I,1),IREPC))
 220        CONTINUE
            IF (IPRINT .GE. 2) THEN
              WRITE (LUPRI ,'(/2X,A,A/)')
     &         'MP2 corr. electronic field at nucleus ',NAMDEP(NATOM)
              IOFF = ICHAR('w')
              DO I = 1,3
                WRITE(LUPRI,'(4X,A1,A,E20.10)') CHAR(IOFF+I),
     &              '-component: ',EFNCORR(I,NATOM)
              ENDDO
            END IF
            NATOM = NATOM + 1
          END IF
 310    CONTINUE
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
      CALL SETATM(WORK(KDATOM),NATOM,IDUM)
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

