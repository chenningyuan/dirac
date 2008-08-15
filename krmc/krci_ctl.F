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
C VERSION : $Revision: 7316 $
C DATE    : $Date: 2008-08-07 16:46:22 +0200 (Do, 07 Aug 2008) $
C FILE    : krmcopt.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCI_CALC(WFTYP,WFCONV,WORK,LWORK)
C***********************************************************************
C
C     Driver routine to KRCI_CALC_1
C     Get CI wave function.
C     Either LUCIAREL or GOSCI or ...
C
C     Based on PSIOPT.
C
C     Written by S. Knecht and H. J. Aa. Jensen - Aug 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
      PARAMETER(D0=0.0D0)
C
      CHARACTER WFTYP*(*), KRCLCTYP*8
      LOGICAL   WFCONV
      DIMENSION WORK(LWORK)
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "dcbgen.h"
C
      LOGICAL EX
      INTEGER RENAME
C
#include "memint.h"
C
      CALL QENTER('KRCI_CALC')
C
C
C     Determine the wave function type
C     --------------------------------
C
      IF ( WFTYP .EQ. 'SCF  ' .or. WFTYP .EQ. 'MCSCF') THEN
C
C        that should not happen...
C
         WRITE(LUPRI,'(/A)')
     $   '*** ERROR in KRCI_CALC *** wrong wave function type '
         CALL QUIT(
     &   '*** ERROR in KRCI_CALC: wave function type not allowed ***')
C
      ELSE IF ( WFTYP .EQ. 'KR-CI ') THEN
        MCTYPE = JGAS
      END IF
C
C     Say Hi!
C     -------
C
      CALL KRCIHI
C
C     ***************************************************
C     *** Open files needed *****************************
C     ***************************************************
C     *** Set logical units ***
C
      CALL SETKRCI(LUPRI)

C
C     ***************************************************
C     *** Initialize configurational and orbital      ***
C     *** parameters.                                 ***
C     ***************************************************
C
C
      CALL SETDCBIDX()
      CALL RSETWOP(WORK,KFREE,LFREE) ! needed ? no orb-rot in CI !!!
      CALL FLSHFO(LUPRI)
C
#if defined (VAR_MPI2)
C
C     ***************************************************
C     *** Initialize common block information         ***
C     ***************************************************
C
      CALL MPIXCALC_TYPESZ()
C
#endif
C
C     ***************************************************
C     *** Memory allocation *****************************
C     ***************************************************
C
      CALL MEMGET('REAL',   KCMO  ,NCMOTQ  ,WORK,KFREE,LFREE)
      IF (OPT_CIPROGRAM .EQ. 'LUCIAREL') THEN
         CALL MEMGET('REAL',KCREF ,0       ,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KCREF ,NZCONFQ ,WORK,KFREE,LFREE)
      END IF
C
C     ***************************************************
C     *** Call DRIVER for CI calculation ****************
C     ***************************************************
C
      CALL KRCI_CALC_1(WORK(KCMO),WORK(KCREF),
     &                 WFTYP,WFCONV,WORK,KFREE,LFREE)
C
C     ***************************************************
C     *** Close files ***********************************
C     ***************************************************
C
      CLOSE(LUKRM3,STATUS='KEEP')
      IF (LUKRM5.GT.0) CLOSE(LUKRM5,STATUS='KEEP')
C
C     Save orbitals on file KRMCSCF if it does not exist.
C     Other DIRAC modules may expect to find it, e.g., 
C     MULLPOP.
C
      INQUIRE ( FILE = 'KRMCSCF', EXIST = EX )
      IF( .NOT. EX ) THEN
        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','KRCI  ')
        CALL NEWLAB('*KRCI   ',LUKRMC,LUPRI)
        CALL WRTKRMC(LUKRMC,'NEWORB  ',WORK(KCMO),NCMOTQ)
        CLOSE(LUKRMC,STATUS='KEEP')
      END IF
C
C
C     ***************************************************
C     *** Memory deallocation ***************************
C     ***************************************************
C
      CALL MEMREL('KRCI_CALC',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('KRCI_CALC')
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRCI_CALC_1(CMO,CREF,WFTYP,WFCONV,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get CI wave function.
C
C     Based on PSIOP1.
C
C     Written by S. Knecht and H. J. Aa. Jensen - Aug 2008
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "mx2orb.h"
#include "dcbnrt.h"
#include "dcbpsi.h"
#include "dcbkrmc_itinfo.h"
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "frame.h"
#include "krciprop.h"
C
      LOGICAL   WFCONV, DO4INDEX, TOBE
      DIMENSION WORK(*), CMO(*), CREF(*)
      CHARACTER WFTYP*(*), SECTID*12, CPUTID*12, WALLTID*12
C
      CALL QENTER('KRCI_CALC_1')
C
      CALL GETTIM(CPUTOT1,WALLTOT1)
      KFRSAV = KFREE
C
C
C     ***************************************************
C     *** Initialize some parameters ********************
C     ***************************************************
C
      DO I = 1,2
         IOPT_ISYMOP(I) = 1
         IOPT_IFCKOP(I) = 1
         IOPT_IHRMOP(I) = 1
      END DO
      JTRLVL = 0
      ITMAC  = 0
      DO4INDEX = ( MCTYPE .GE. JMCMIN )
C     ... MCTYPE .lt. JMCMIN is SCF
      EMCOLD = D0
      WFCONV = .FALSE.
C
      CALL DZERO(GNORM,5)
C
C     Integrals for KR-CI calculation.
C
      CALL INTCON(IOPT_INTFLG,IOPT_INTBUF,IOPT_INTDEF,GNORM(5),
     &            OPT_CNVINT,ITMAC,IOPT_ITRINT,OPT_INTTYP)
C
C
C     ***************************************************
C     *** Get CI vectors ********************************
C     ***************************************************
C
C
#ifdef UNDEF
      write(LUPRI,*) '(1) iopt_intflg=',iopt_intflg
      write(LUPRI,*) '(1) iopt_intbuf=',iopt_intflg
      write(LUPRI,*) '(1) iopt_intdef=',iopt_intflg
      write(LUPRI,*) '(1) gnorm(5)=',gnorm(5)
      write(LUPRI,*) '(1) opt_cnvint=',opt_cnvint
      write(LUPRI,*) '(1) itmac=',itmac
      write(LUPRI,*) '(1) iopt_itrint=',iopt_itrint
      write(LUPRI,*) '(1) opt_inttyp=',opt_inttyp
#endif
C     Do not write orbitals to file KRMCSCF in ROPTST.
C     KRMCSCF may not be opened.
      IWRT = 0
C
C     Read start orbitals and get CI vectors:
C
      CALL ROPTST(IOPTST,IWRT,CMO,CREF,IPROPT,WORK,KFREE,LFREE)
C     CALL ROPTST(ISTART,IWRT,CMO,CREF,IPRINT,WORK,KFREE,LFREE)
C
C     only 4-index integral transformation done
C
      IF (KTRLVL.eq.5) GOTO 3000
C
C     ***************************************************
C     *** KR-CI property calculation ********************
C     ***************************************************
C
      IF( DOPROP_KRCI ) CALL RPROP_KRCI(WORK(KFREE),LFREE)
C
C
C     ***************************************************
C     *** KR-CI finished - no task left *****************
C     ***************************************************
C
C
 3000 CONTINUE
C
      WFCONV = .TRUE.
C
      CALL GETTIM(CPUTOT2,WALLTOT2)
      CPUTID = SECTID(CPUTOT2-CPUTOT1)
      WALLTID = SECTID(WALLTOT2-WALLTOT1)
      WRITE(LUPRI,9300) CPUTID,WALLTID
C
 9300 FORMAT(/'>>> TOTAL CPU (WALL) TIME IN KR-CI: ',A,'(',A,')')
C
      CALL MEMREL('KRCI_CALC_1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('KRCI_CALC_1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krcihi */
      SUBROUTINE KRCIHI
C***********************************************************************
C
C     Write welcome message.
C
C     Written by J. Thyssen - Jan 18 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "mx2orb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      CALL TITLER('KR-CI calculation ','*',199)
C
      WRITE(LUPRI,9000)
 9000 FORMAT(
     &     3X,'This is output from DIRAC KR-CI',/,
     &     3X,'- a relativistic four-component CI',
     &     1X,'wave function program.',
     &     ///,
     &     3X,'General structure: ',/,
     &     3X,'  Stefan Knecht and Hans Joergen Aa. Jensen ',//,
     &     3X,'Integral transformation: ',/,
     &     3X,'  Luuk Visscher, Jon K. Laerdahl, and Trond Saue',//,
     &     3X,'GOSCI CI code: ',/,
     &     3X,'  Joern Thyssen and Hans Joergen Aa. Jensen',//,
     &     3X,'LUCIAREL CI code: ',/,
     &     3X,'  Timo Fleig and Jeppe Olsen',//,
     &     3X,'Parallel LUCIAREL CI code: ',/,
     &     3X,'  Stefan Knecht, Hans Joergen Aa. Jensen and Timo Fleig',
     &     //,
     &     1X,79('*'))
!
      WRITE(LUPRI,9110)
 9110 FORMAT(/' This module is published in:',//,
     &     4X,'GOSCI: J Thyssen, T Fleig, and H J Aa Jensen',/,
     &    20X,'   J Chem Phys 129, 034109 (2008), suppl. material.',/,
     &     4X,'DIRAC-LUCIAREL: T Fleig, J Olsen, and L Visscher',/,
     &    20X,'   J Chem Phys, 119,6 (2003) 2963',/,
     &     4X,'PARALLEL LUCIAREL:',/,
     &     4X,'                S. Knecht, H J Aa Jensen, and T Fleig',/,
     &    20X,'   manuscript in preparation.', //,
     &     1X,79('*'),/)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setkrci */
      SUBROUTINE SETKRCI(LUPRI)
C***********************************************************************
C
C     Assign logical units and open files.
C
C     Based on SETKRMC.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcbpsi.h"
C
      LOGICAL EX
C
C     LUFCK1 definition needs to be done here
C     as SETDHF is not called in the KRCI run:
C
      LUFCK1 = 10
C
C     Units 50-59 are reserved for KRMC and KRCI.
C     Units 50-56 are opened here, 57-58 in krmctra.F.
C     (hjaaj March 2002).
C     Unit 59 is currently used for KRMCOLD, restarting option.
C     (tf August 2003)
C
C     LUKRMC ('KRMCSCF') not used in KRCI, but not touched here
C            because LUKRMC = 50 is set globally in dirac/dirset.F
C     LUKRM1 ('KRMC_KAPPA') not used in KRCI
C     LUKRM2 ('KRMC_CMO') not used in KRCI
C     LUKRM3 ('KRMC_FOCK') used in KRCI for FOCK matrix and propertz
C            matrices
C     LUKRM4 ('KRMC_LRED') not used in KRCI
C     LUITFO ('KRMC_ITINFO') not used in KRCI
C
      LUKRM1 = -51
      LUKRM2 = -52
      LUKRM3 =  53 ! note the invisible plus
      LUKRM4 = -54
      LUITFO = -56
C
C     LUKRM5 is used for reading KR-MCSCF orbitals in ROPTST.
C     First check if KRMCSCF available, meaning this calculation follows
C     a KR-MCSCF calculation. Then we will use these orbitals. If
C     KRMCSCF is not available, we check for 'KRMCOLD' with orbitals
C     from a previous calculation.
C
      LUKRM5 = 59
      INQUIRE ( FILE = 'KRMCSCF', EXIST = EX )
      IF( EX ) THEN
         CALL OPNFIL(LUKRM5,'KRMCSCF','OLD','STKRCI')
         REWIND LUKRM5
      ELSE
         INQUIRE ( FILE = 'KRMCOLD', EXIST = EX )
         IF (EX) THEN
            CALL OPNFIL(LUKRM5,'KRMCOLD','OLD','STKRCI')
            REWIND LUKRM5
         ELSE
            LUKRM5 = -59
         END IF
      END IF
C
C     Delete 2-el integral files from KRMCSCF or from MOLTRA, thez
C     are probably not correct for KRCI. However, if neither KRMCSCF or
C     MOLTRA have been active in this Dirac run, the we assume we are
C     restarting KRCI and that any integrals are correct. (aug08, HJAaJ)
C
      IF (DOKRMC .OR. DOTRA) THEN
         CALL OPNFIL(LUKRM3,'4INDINFO','UNKNOWN','STKRCI')
         CLOSE(LUKRM3,STATUS='DELETE')
         CALL OPNFIL(LUKRM3,'KRMC_FOCK','UNKNOWN','STKRCI')
         CLOSE(LUKRM3,STATUS='DELETE')
      END IF
C
C     If old 'KRMC_FOCK' still exists, then we are restarting a KRCI.
C
      INQUIRE ( FILE = 'KRMC_FOCK', EXIST = EX )
      IF( EX ) THEN
        CALL OPNFIL(LUKRM3,'KRMC_FOCK','OLD','STKRCI')
      ELSE
        CALL OPNFIL(LUKRM3,'KRMC_FOCK','UNKNOWN','STKRCI')
        CALL NEWLAB('SODLABEL',LUKRM3,LUPRI)
      END IF
      REWIND LUKRM3
C
      RETURN
      END
C=======================================================================
C     end of file krci_ctl.F
C=======================================================================
