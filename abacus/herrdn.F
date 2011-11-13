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
C VERSION : $Revision: 1.89.2.2 $
C DATE    : $Date: 2004/09/01 08:04:10 $
C FILE    : herrdn.F
C
C
#ifdef UNDEF
/* Comdeck log */
960705-kr:
-- We use the small component basis for keeping a possible extra STO-3G
   basis needed when using Hueckel initial guess.
951115-hjaaj:
-- no warning for THRS .le. 0 if input field blank
-- implemented free format reads for IRIX (has removed character reads
   which does not work on IRIX and gives nulls on AIX and ?).
CNTINP: removed unused 'NAVN'; always read NAMN with '(A4)' and never
   with free format; if free format doesn't work then try formatted;
ZMAT: Read NAMN with '(A4)', allow free input for IRIX;
   for CRAY formatted input made last field for CHARGE very long;
   removed old END=204 branches which were reminiscences from when
   the input was read directly from file.
BUILDZ: changed format and output of internal coordinates;
   removed printing of zero cartesian coords. before printing of
   the correct cartesian coords.
950428-kr/hjaaj:
BASOUT,ORBOUT: removed NPRIM,NPRIMD,NORB,NORBD in parameter list to
 ORBOUT; we think they are related to an old organization of
 small component basis (before IQM(*,2) etc.)
941219-hjaaj:
BASINP: s/IQMAX/NHTYP/
950208-vebjorn
READI1: NSYMOP is first read as a text to determine if symmetry is to
        be added (indicated by blanks).
SYMINP: Has been moved after BASINP, because symmetry has to be added before
        SYMINP is executed. The full group string CLASS has been added
        as a parameter to the subroutine.
SYMADD: Completely new module that adds symmetry on request from user. The
        full group of the molecule is returned as well.
BASINP: Has been moved above SYMADD and SYMINP. As a result the updating
        of the array ISTBNU had to be moved to BASPRO. As this is processing
        rather than input, the move makes good sense.
BASPRO: Now takes care of the generation of the symmetry-dependant atoms
        through the array ISTBNU.
941119-tsaue:
READIN: added extra index to arrays IQM,NBLCK,JCO,NUC,NRC,SEG,NOINT,ALPHA,
        CPRIMU ta accomodat both large and small components
ACPORB: implemented readin of exponents only
SYMGRP: New routine determining irreps, characters etc.
DBLGRP: new routine for handling double group symmetry
SYMOUT: moved test for IPREAD out of routine,
   for IPREAD=0 the routine was pretty meaningless
BASINP: major modifications - check it carefully out !!!
        may consider replacing NOORBT with MAXIQM wich has
        more information
CNTINP: info on angular momentum moved to BASINP
941114-hjaaj:
ORBINX: modified CALL SEGORB so no line is longer than 72 char
        (ORBINX not used)
940207-hjaaj
READI1: 'ANG' was not initialized; errors in 'ANG' message format
931220-hjaaj
WRONEL: corrected calc. of INAMN and IGTO for LUONEL
931022-hjaaj
REAINP: INPTST = IPREAD .gt. 3 (INPTST was not defined)
expanded tabs with 'expand'
READI1: If ID3 then reset to blank (for abawalk)
BASINP: changed format for nuclear charge from F4.1 to F4.0
 (such that an entry '6' is read as 6. and not as 0.6 )
CNTINP: read in old 3F10 format if error reading coordinates
        changed format to 3F20.0 instead of 3F20.15
#endif
C  /* Deck reainp */
      SUBROUTINE REAINP(WORD,IPRDEF,RELCAL,TSTINP)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxsymm.h"
C
C     Note that the parameter MAXPRD also must be set in subroutine REAINI
C
#if defined (SYS_T3D)
      PARAMETER ( MAXPRD = 22 )
#else
      PARAMETER ( MAXPRD = 30 )
#endif
      PARAMETER ( NTABLE = 10, D1 = 1.0D0)
      LOGICAL NEWDEF, RELCAL, TSTINP
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "cbirea.h"
#include "ccom.h"
C
      DATA TABLE /'.PRINT ', '.OLDNOR', '.MAXPRI', '.BIGVEC',
     *            '.MOLINP', '.SYMTHR', '.CM FUN', '.ZCMVAL',
     *            '.MOLPRI', '.UNCONT'/
C
      NEWDEF = (WORD .EQ. '*READIN')
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
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in REAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in REAINP.')
    1          CONTINUE
                  READ (LUCMD, *) IPREAD
                  IF (IPREAD .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    2          CONTINUE
                  OLDNRM = .TRUE.
               GO TO 100
    3          CONTINUE
                  READ (LUCMD, *) MAXPRI
               GO TO 100
    4          CONTINUE
                  BIGVC = .TRUE.
               GO TO 100
    5          CONTINUE
                  LUMLCL = LUCMD
               GO TO 100
 6             CONTINUE
                  READ (LUCMD, *) TOLLRN
               GO TO 100
 7             CONTINUE
                  READ (LUCMD,*) LCMMAX, CMSTR, CMEND
                  NCMSTR = NINT(2*CMSTR)
                  NCMEND = NINT(2*CMEND)
               GO TO 100
 8             CONTINUE
                  READ (LUCMD,*) ZCMVAL
               GO TO 100
 9             CONTINUE
                  PRIBAS = .TRUE.
               GO TO 100
 10            CONTINUE
                  UNCONT = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in REAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in REAINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for READIN:',0)
         IF (IPREAD .NE. IPRDEF) THEN
            WRITE (LUPRI,'(A,I5)') ' Print level in READIN:',IPREAD
         END IF
         IF (OLDNRM) THEN
            WRITE (LUPRI,'(A)')
     &         ' Dalton 1.0 normalization of spherical GTOs used.'
         END IF
         IF (TOLLRN .NE. 0.0D0) THEN
            WRITE (LUPRI,'(A,G9.2)') ' Symmetry threshold: ',TOLLRN
         END IF
         IF (LUMLCL .EQ. LUCMD) THEN
            WRITE (LUPRI,'(/A)') ' MOLECULE read from input file'
         END IF
         IF (UNCONT) THEN
            WRITE (LUPRI,'(/A)')
     &         ' Uncontracted basis forced, irrespective of '//
     &         'basis input file.'
         END IF
         IF (MAXPRI .NE. MAXPRD) THEN
            WRITE (LUPRI,'(/A,I5)')
     *         ' Maximum number of primitives per integral block :',
     *         MAXPRI
         END IF
         IF (BIGVC) THEN
            WRITE (LUPRI,'(/A)')
     *         ' Primitives from different centers treated '//
     *         'simultaneously.'
            WRITE (LUPRI,'(/A)')
     &       ' Option .BIGVEC not implemented in present version.'
             CALL QUIT('Option .BIGVEC not allowed.')
         END IF
         IF (LCMMAX .GT. -1) THEN
            WRITE (LUPRI,'(/A,I5)')
     &           ' Center of mass basis functions used with l(max) =',
     &           LCMMAX
            WRITE (LUPRI,'(2(A,F4.1))')
     &           ' Sequence starts at ',CMSTR,' and ends at ', CMEND
            IF (ZCMVAL .NE. D1)
     &           WRITE (LUPRI,'(A,F3.1)') ' Net charge of ionic core  ',
     &           ZCMVAL
            IF (LCMMAX .GT. 4) THEN
               WRITE (LUPRI,'(/A)')
     &              ' Center of mass bf. not defined for l .gt. 4'
               CALL QUIT('Error in *READIN input processing')
            END IF
         END IF
         IF (PRIBAS) THEN
            WRITE (LUPRI,'(/A)')
     &           ' The molecular input will be dumped on DALTON.BAS'
         END IF
C        IF (GENCON) THEN
C           WRITE (LUPRI,'(/A)')
C    *         ' Routines for general contraction used.'
C        END IF
         WRITE (LUPRI,'()')
      END IF
      RETURN
      END
C  /* Deck reaini */
      SUBROUTINE REAINI(IPRDEF,RELCAL,TSTINP)
C
C     Define THRS in CCOM
C
C
C     Initialize /CBIREA/
C
#include "implicit.h"
C
C     Note that MAXPRD also must me set in subroutine REAINP
C
#if defined (SYS_T3D)
      PARAMETER (MAXPRD = 22)
#else
      PARAMETER (MAXPRD = 30)
#endif
      LOGICAL RELCAL, TSTINP
#if defined (VAR_MPI)
#include "infpar.h"
      INCLUDE 'mpif.h'
#endif
#include "maxaqn.h"
#include "ccom.h"
#include "mxsymm.h"
#include "cbirea.h"
      THRS = 1.0D-15
C
      IPREAD = IPRDEF
      LUMLCL = 9
      BIGVC  = .FALSE.
      DIRAC  = RELCAL
      GENCON = .TRUE.
      INPTST = TSTINP
      BASIS  = .FALSE.
      ATOMBA = .FALSE.
      PRIBAS = .FALSE.
      OLDNRM = .FALSE.
      TOLLRN = 0.0D0
      ZCMVAL = 1.0D0
      LCMMAX = -1
C
C     In MPI runs this initialization needs to be done by the master
C     because this is the only provess to read the DIRAC.INP file.
C
C     TODO : Check that the other variables are indeed not modified
C     in DIRAC because the slaves will always get the default, unless
C     the variable is explictly communicated after reading !!
C
#if defined (VAR_MPI)
      IF (MYTID.EQ.MPARID) THEN
C        This is the master, initialize as usual
C        Slaves are already initialized before coming here.
         UNCONT = .FALSE.
         MAXPRI = MAXPRD
      ENDIF
#else
      UNCONT = .FALSE.
      MAXPRI = MAXPRD
#endif
      RETURN
      END
C  /* Deck readin */
      SUBROUTINE READIN(WORK,LWORK,HERMIT)
C     Based on
C       READIN:Input processing routine for  M O L E C U L E.
C         Jan Almlof, Stockholm, Dec. 1971.
C         - adaption for ABACUS (T.Helgaker, University of Oslo)
C         - Occams razor by T.Saue, University of Oslo, March 10 - 1993
C           Major surgery includes:
C               - common block CCOM :
C                       - added KHK(MXQN),KCK(MXQN),NHKOFF(MXQN)
C                       - added GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6)
C               - common block NUCLEI:
C                       - added GAUNUC, GNUEXP(MXCENT)
C                       - added NAMN(MXCENT)
C               - NEW common block FRAME
C                       - added POTNUC, DPNUC(3)
C               - common block CBIREA:
C                       - IPRINT changed to IPREAD
C                       - added INPTST,FAMILY
C               - common block SHELLS:
C                       - added NLRGSH,NSMLSH,NLARGE,NSMALL,NORBS
C                       - added LCLASS(MXSHEL),NOTWOC(MXSHEL)
C               - common block MOLINP: deleted !!!
C               - common block MOLINC: deleted !!!
C               - common block HRUNIT: deleted !!!
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "abainf.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbirea.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL HERMIT,NEWFILE,SET
      CHARACTER KEYWRD*6
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
#include "infpar.h"
C
C     A SET variable has been replaced by RDINPC. This allows external control
C     of READIN processing. The processing should be performed once every
C     iteration.
C
      IF (RDINPC) RETURN
      RDINPC = .TRUE.
C
      IF (IPREAD .GE. 11) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPREAD .GT.  0) THEN
         CALL TITLER('Output from READIN input processing','*',118)
      END IF
      LUINFO = LUMLCL
C
C     PARHER is true for slaves, false for master.
C     If (PARHER) then RDLINE reads from common block transferred
C     from master, not from file LUINFO.
C
      NEWFILE = LUINFO .NE. LUCMD .AND. .NOT. PARHER
      IF (NEWFILE) THEN
#if !defined (SYS_IBM)
         OPEN (LUINFO,FILE='MOLECULE.INP',STATUS='OLD',
     *         FORM='FORMATTED',ERR=5000)
#endif
#if defined (SYS_IBM)
         OPEN (LUINFO,FILE='MOLINP',STATUS='OLD',
     *         FORM='FORMATTED',ERR=5000)
#endif
      END IF
C
      IF (.NOT.PARHER) REWIND(LUINFO)
      NMLINE = 0
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A6)') KEYWRD
#if !defined (PRG_DIRAC)
      IF (KEYWRD .EQ. 'BASIS ') THEN
         BASIS = .TRUE.
      ELSE IF (KEYWRD .EQ. 'ATOMBA') THEN
         ATOMBA = .TRUE.
      ELSE IF (KEYWRD .NE. 'INTGRL') THEN
         GOTO 5010
      END IF
#else
C     In Dirac we do the readin the same way whatever KEYWRD is.
#endif
C
C     ********************************************
C     ***** Allocate memory **********************
C     ********************************************
C
      IF (DIRAC) THEN
C        Number of basis sets : large, small, fit LL dens, fit SS dens
         KSETS = 4
      ELSE
C        Number of basis sets : regular, Huckel, fit
         KSETS = 3
      ENDIF
      KATOM  = MXCENT
      KANG   = 2*MXQN+1
      KPRIM  = MAXPRI
      KCMAT  = KPRIM*KPRIM
      KNONT  = 1
      KIQM   = KNONT  +  KATOM
      KBLCK  = KIQM   +  KATOM*KSETS
      KJCO   = KBLCK  +  KATOM*KSETS
      KJCO2  = KJCO   +  KATOM*KANG*KSETS
      KLAST  = KJCO2  +  KATOM*KANG
C hjj: take as many blocks as fit in memory
C      (to be able to run tests with memory < 1 mio DP words)
      LWRK   = LWORK  -  KLAST + 1
      LBLOCK = 6 + 2*KPRIM + 4*KCMAT
C      (we reserve size of one block for later allocations)
      LRESRV = MAX(LBLOCK, 3*MXAQN+MXCORB, 2*MXCOOR*MXCOOR) + 10000
C                          ^- for sympro   ^-- for nucpro, pricar
      KBLOCK = (LWRK-LRESRV)/LBLOCK
      KBLOCK = MAX(KBLOCK,10)
C              ... minimum 10 blocks
      KBLOCK = MIN(KBLOCK,MXSHEL)
C              ... maximum MXSHEL blocks
C              (otherwise the program will stop later;
C              and no reason to allocate all of memory here).
C
      KNUC   = KLAST
      KNRC   = KNUC   +  KBLOCK*KSETS
      KSEG   = KNRC   +  KBLOCK*KSETS
      KA     = KSEG   +  KBLOCK*KSETS
      KC     = KA     +  KBLOCK*KPRIM*KSETS
      KCC    = KC     +  KBLOCK*KCMAT*KSETS
      KLAST  = KCC    +  KBLOCK*KCMAT*KSETS
      IF(KLAST.GT.LWORK) GOTO 5020
      LWRK   = LWORK  -  KLAST + 1
C
C     ********************************************
C     ***** Read input file and process data *****
C     ********************************************
C
      CALL READI1(LUINFO,WORK(KLAST),LWRK,WORK(KNONT),WORK(KIQM),
     &            WORK(KBLCK),WORK(KJCO),WORK(KJCO2),
     &            WORK(KNUC),WORK(KNRC),WORK(KSEG),WORK(KA),
     &            WORK(KC),WORK(KCC),KATOM,KANG,KSETS,KBLOCK,KPRIM, 
     &            HERMIT)
C     CALL RDLINE(LUINFO)
C     READ (MLINE(NMLINE),'(A6)') KEYWRD
C     IF (KEYWRD .NE. 'FINISH') GOTO 5030
      IF (IPREAD .GT. 20) CALL REAPRI
      RETURN
C
C Error messages:
C
 5000 CONTINUE
      WRITE (LUPRI,'(//A)') ' *** ERROR *** unable to open MOLECULE.INP'
      CALL QUIT(' *** ERROR (READIN) unable to open input file LUINFO ')
 5010 CONTINUE
      CALL QUIT('Wrong keyword. Expecting INTGRL or BASIS')
 5020 CONTINUE
      CALL STOPIT('READIN','READI1',KLAST,LWORK)
C5030 CONTINUE
      CALL QUIT('Wrong keyword. Expecting FINISH')
      END
C  /* Deck readi1 */
      SUBROUTINE READI1(LUINFO,WORK,LWORK,NONT,IQM,NBLCK,JCO,JCO2,
     &                  NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                  KATOM,KANG,KSETS,KBLOCK,KPRIM,HERMIT)
C******************************************************************************
C
C     Reads and processes input file of molecular data
C
C       **** Temporary variables used in READIN ****
C
C       KATOM   - atomic types
C       KANG    - angular momentum
C       KBLOCK  - AO-block
C       KPRIM   - primitives
C
C      (Large/small)
C       NONT            (KATOM)         - number of symmetry independent
C                                         centres for this atomic type
C       IQM             (KATOM)         - highest L-value
C       NBLCK           (KATOM)         - number of AO-blocks for a given
C                                         atomic type
C       JCO             (KANG,KATOM)    - number of AO-blocks for a given
C                                         atomic type and L-value
C       NUC             (KBLOCK)        - number of uncontracted functions
C                                         in a given AO-block
C       NRC             (KBLOCK)        - number of contracted functions
C                                         in a given AO-block
C       SEG             (KBLOCK)        - TRUE: segmented contraction
C       ALPHA           (KPRIM,KBLOCK)  - exponents
C       CPRIM           (KPRIM,KPRIM,KBLOCK)  - normalized contraction
C                                               coefficients
C       CPRIMU          (KPRIM,KPRIM,KBLOCK)  - contraction coefficients
C
C     LOCAL VARIABLES
C       IPRIMA - runs over primitives for a given atomic type
C       ISHELL - runs over all shells
C       IPRIM - runs over all primitives
C       NUCIND2 - runs over all symmetry independent centers
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0, THRMIN = 1.D-15)
#include "clsfmm.h"
#include "cbirea.h"
      LOGICAL SEG,HERMIT,DOOWN,ZMAT
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),NBLCK(KATOM,KSETS),
     &          JCO(KANG,KATOM,KSETS),JCO2(KANG,KATOM),
     &          NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),SEG(KBLOCK,KSETS),
     &          ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS),
     &          WORK(LWORK)
C
C       Local variables
C
      CHARACTER*80 BASREF(10,MXCENT,MXBSETS)
C
      CHARACTER*1 KASYM(3,3),ID3,SPDCAR,CRT
      CHARACTER*72 TTITLE(2)
      CHARACTER*2 SYMTXT
      CHARACTER*15 CLASS
      DIMENSION IFXYZ(3)
      LOGICAL ANG, ADDSYM
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "infpar.h"
#include "dcbrkb.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
      LOGICAL LINEAR
#endif
#include "ibtfun.h"
      ZMAT = .FALSE.
      CALL IZERO(IRKB,MXPRIM)
C
C*****************************************************************************
C     Card 1: Keyword (INTGRL/BASIS) - read in READIN
C*****************************************************************************
C*****************************************************************************
C     Card 2: If basis set library, read additional card containing basis set
C*****************************************************************************
      IF (BASIS) THEN
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(A)') BASNAM
         IF (IPREAD .GT. 1) WRITE (LUPRI,'(/A/5X,A/A)') ' Basis set',
     &      BASNAM,' from the basis set library will be used.'
      END IF
C*****************************************************************************
C     Card 2-3: Title cards
C*****************************************************************************
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A72)') TTITLE(1)
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A72)') TTITLE(2)
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Title Cards',1)
         WRITE (LUPRI,'(2X,A,/,2X,A)') TTITLE(1),TTITLE(2)
      END IF
C*****************************************************************************
C     Card 4:
C       CRT     - flag for spherical harmonics or "your own scheme"
C       NONTYP  - number of atomic types
C       SYMTXT  - Symmetry operations read as a text to test if
C                 symmetry is to be added automatically
C       KCHARG  - Molecular charge
C       NSYMOP  - number of independent twofold symmetry operations of
C                 point group to be used
C       KASYM   - 3x3 character array specifying the basic symmetry operations;
C                 specify axis that are reversed under the operation
C       THRS    - integral threshold
C******************************************************************************
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)')
     &     CRT,NONTYP,KCHARG,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),ID3,THRS
      IF(NONTYP.GT.KATOM) GOTO 5000
C*****************************************************************************
C Check integral threshold: THRS must be .gt. 0
C       ( log(thrs*thrs) taken in oneint) (851005-hjaaj)
C*****************************************************************************
C
      IF (MLINE(NMLINE)(21:30) .EQ. '          ') THEN
C        nothing specified, default minimum value used
         THRS = THRMIN
      ELSE IF (THRS .LT. THRMIN) THEN
         IF (.NOT.PARHER) THEN
            WRITE (LUPRI,'(/2X,A)') '*** WARNING from READIN ***'
            WRITE (LUPRI,'(2(/2X,A,1P,D12.2))')
     &      'Threshold for discarding integrals was',THRS,
     &      'Threshold is reset to minimum value   ',THRMIN
         END IF
         THRS = THRMIN
      ENDIF
C*****************************************************************************
C Check if symmetry is to be added automatically
C       ADDSYM  - flag for automatic determination of symmetry
C*****************************************************************************
      ADDSYM = (SYMTXT .EQ. '  ')
      IF (.NOT. ADDSYM) THEN
         READ(SYMTXT,'(I2)') NSYMOP
      ELSE
         NSYMOP = 0
      END IF
      NMLAU = NMLINE
      IF (ID3 .NE. ' ') THEN
        ANG = .TRUE.
        IF (.NOT. PARHER)
     &       WRITE (LUPRI,'(/2X,A/10X,A,F11.8,A2)')
     &          'Coordinates are entered in Angstroms'//
     &          ' and converted to atomic units.',
     &          '- Conversion factor : 1 bohr =',XTANG,' A'
      ELSE
         ANG = .FALSE.
      ENDIF
C*****************************************************************************
C Process Cartesian to spherical or your own transformation
C     note: SPHINP will read additional input lines if DOOWN true
C*****************************************************************************
      DOCART = CRT .EQ. 'C' .OR. CRT .EQ. 'c'
      DOOWN  = CRT. EQ. 'X'
      IF (DIRAC .AND. .NOT. (DOCART .OR. DOOWN)) THEN
         IF (IPREAD .GT. 0) THEN
            WRITE(LUPRI,'(/3(A/))')
     &      'WARNING!!! Cartesian Gaussians must be specified',
     &      '           in molecule input file.',
     &      '           DIRAC proceeds with readin assuming '//
     &      '           Cartesian Gaussians.'
         ENDIF
         DOCART = .TRUE.
      END IF
      CALL SPHINP(LUINFO,WORK,LWORK,DOOWN)
C
C     ************************************************
C     ***** Read orbital and geometry input data *****
C     ************************************************
C
      KISGEN = 1
      KLAST  = KISGEN + KBLOCK*KSETS
      LWRK   = LWORK  -  KLAST + 1
#if !defined (MOD_DNF)
      CALL BASINP(LUINFO,WORK(KLAST),LWRK,NONTYP,NONT,
     &            IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,WORK(KISGEN),
     &            KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG,ZMAT)
#else
      CALL CENINP(LUINFO,WORK(KLAST),LWRK,NONTYP,NONT,
     &            IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,WORK(KISGEN),
     &            KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG)
#endif
C
C     *************************************************
C     ***** Process ZMAT coordinates if requested *****
C     *************************************************
C
      IF (ZMAT) CALL BUILDZ(IPREAD,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,
     &                  BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,
     &                  KANG,KBLOCK,KPRIM,DOOWN)

C
C     *************************************************
C     ***** Add symmetry to molecule if requested *****
C     *****   Returns CLASS, the full symmetry    *****
C     *************************************************
C
      LINEAR = .FALSE.
      CLASS = 'N/A'
      IF (ADDSYM) THEN
         CALL SYMADD(WORK,LWORK,NONTYP,NSYMOP,NONT,KATOM,KASYM,CLASS)
         IF (CLASS(3:4).EQ.'oo') LINEAR = .TRUE.
      END IF
C
C     ****************************************************
C     ***** Adds an STO-3G basis if Extended Hueckel *****
C     *****        wanted as starting orbitals       *****
C     ****************************************************
C
C      IF (ADDSTO) CALL STOADD()
C
C     ***************************************
C     ***** Process symmetry input data *****
C     ***************************************
C
      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)
C
C     ***************************************************
C     ***** Process orbital and geometry input data *****
C     ****************************************************
C
      CALL BASPRO(WORK,LWORK,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,NUC,NRC,
     &            SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,KBLOCK,KPRIM,
     &            DOOWN)
C
C     *************************************************
C     ***** Print orbital and geometry input data *****
C     *************************************************
C
      CALL BASOUT(WORK,LWORK,NONTYP,NONT,IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &            KBLOCK,KPRIM,HERMIT)
C
      IF (IPREAD .GE. 2) THEN
         CALL HEADER('Copy of input to READIN',0)
         WRITE (LUPRI,'(A)') (MLINE(I), I = 1, NMLINE)
         WRITE (LUPRI,'(/)')
      END IF

C     Default threshold for distinguishing classical and
C     non-classical contributions to Coulomb integrals

      THRCLS = THRS
C
C     ****************************************
C     ***** Output on LUONEL *****************
C     ****************************************
C
      IF (HERMIT) CALL WRONEL(TTITLE,NONTYP,NONT,IQM,IFXYZ,KATOM,JCO2,
     &                        KANG)
C
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> READI1 error, no. of atomic types    ',NONTYP,
     *   '                   current maximum number ',KATOM
      CALL QUIT('Too many atomic types')
      END
C  /* Deck baspar */
      SUBROUTINE BASPAR(NSYMOP)
C*****************************************************************************
C
C       Set various basis parameters:
C               KHK(J) - number of spherical (cartesian) components for given J
C               KCK(J) - number of Cartesian components for given J
C                      - tabulate incomplete Gamma function
C                      - determine Cartesian powers
C               NHKOFF(J) - offset for components in list of l-functions
C               MAXOPR - maximum number of operations to loop over
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
#include "cbirea.h"
C
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "aosotr.h"
#include "ibtfun.h"
C
C     ***** Tabulate Incomplete Gamma Function *****
C
      JMAX = MAX(4,4*(NHTYP - 1) + 2)
      CALL GAMTAB(JMAX)
C
C     ***** Cartesian powers *****
C
      CALL CARPOW
C
C     ***** NHKOFF() *****
C     - offset for components in list of l-functions
C
      IOFF = 0
      DO 100 I = 1, NHTYP
         NHKOFF(I) = IOFF
         IOFF = IOFF + KHK(I)
  100 CONTINUE
C
C     ***** MAXOPR *****
C     - maximum number of operations necessary to loop over
C
      MAXOPR = 1
      IF (NSYMOP .GT. 0) THEN
         MAXLO = 7
         DO 200 I = 1, NUCIND
            MAXLO = IBTAND(ISTBNU(I),MAXLO)
  200    CONTINUE
         II = 1
         DO 300 I = 1, NSYMOP
            IF(IBTAND(II,MAXLO).EQ.0) MAXOPR = 2*II
            II = 2*II
  300    CONTINUE
      END IF
      MAXOPR = MAXOPR - 1
C
      END
C  /* Deck basinp */
      SUBROUTINE BASINP(LUINFO,WORK,LWORK,NONTYP,NONT,IQM,NBLCK,
     &                  JCO,BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &                  KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG,ZMTTST)
C*****************************************************************************
C
C       Read orbital and geometry input data
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (DSM = 1.0D-30)
#include "codata.h"
#include "cbirea.h"
      DIMENSION WORK(LWORK)
      LOGICAL SEG,ANG,NOORBTS,NOFITFS,ZMTTST,CNTBAS
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS),
     &          ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS),ISGEN(KBLOCK)
      CHARACTER SPDCAR*1,BSET*5,ZMATL*1,BASFIL*20,BSKEYWORD*14
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      DIMENSION IBLOCK(MXBSETS)
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbisol.h"
#include "huckel.h"
#include "aosotr.h"
#include "gencon.h"
#include "dcbgen.h"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Additional scratch arrays needed for MULTIBASIS option
C
      PARAMETER (MAXBSFIL=5,KNG=MXAQN,KBLK=500,KPRM=30)
C
      INTEGER IQMT
      DIMENSION JCOT(KNG,MAXBSFIL)
      DIMENSION IBLKAGNMOM(KBLK)
      DIMENSION NUCT(KBLK)
      DIMENSION NRCT(KBLK)
      DIMENSION SEGT(KBLK)
      DIMENSION ALPHAT(KPRM,KBLK)
      DIMENSION CPRIMT(KPRM,KPRM,KBLK)
      DIMENSION CPRIMUT(KPRM,KPRM,KBLK)
      DIMENSION ISGENT(KBLK)
      CHARACTER*80 BASREFT(10,KATOM,MAXBSFIL)
      CHARACTER*20 BSFIL(MAXBSFIL)
      LOGICAL SEGT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
#include "ibtfun.h"
      CNTBAS = .FALSE.
      NUCIND = 0
      NHTYP  = 0
      IBLOCK(1) = 1
      IBLOCK(2) = 1
      IBLOCK(3) = 1
      IBLOCK(4) = 1
      CALL IZERO(IQM,KATOM*4)
      DIRCON(1) = .FALSE.
      DIRCON(2) = .FALSE.
C
C     Run over atomic centers
C     =======================
C
      DO 100 I = 1, NONTYP
         NOORBTS = .TRUE.
         NOFITFS = .TRUE.
         INUC = NUCIND
C
C*****************************************************************************
C        Card 6: ATOMIC DATA
C         Q     - nuclear charge
C         NONT  - number of symmetry independent atoms
C         IQM   - maximum angular quantum number (s=1,p=2,d=3 etc.)
C         NHTYP - maximum angular quantum number for ALL orbitals
C         JCO   - number of AO-blocks for each l-value
C*****************************************************************************
C
         CALL RDLINE(LUINFO)
C
C        Test if Z-matrix for first line (I.eq.1):
C

         IF (I .EQ. 1 .AND. .NOT. DIRAC) THEN
           READ (MLINE(NMLINE),'(A1)') ZMATL
           IF (ZMATL .NE. ' ') THEN
              ZMTTST = .TRUE.
              WRITE (LUPRI,'(/A)')
     &        'Atoms and their geometry is read in Z-matrix format.'
              IF (.NOT. BASIS) THEN
                 WRITE (LUPRI,'(/A)') '*** READI1 ERROR:'
                 CALL QUIT('Z-matrix only implemented with BASIS')
              END IF
           END IF
         END IF
C
C       Read/process data on centers/blocks
C       ===================================
C
        QEXP = D0
        IF(.NOT.ZMTTST) THEN
          IF (DIRAC) THEN
            ADDSTO = .FALSE.
            READ (MLINE(NMLINE),'(BN,1X,F9.0,I5,F20.5)')
     &            Q,NONT(I),QEXP
          ELSE IF (.NOT.ZMTTST) THEN
            NSETS = 1
            IF (ADDSTO) NSETS = 2
            BSET  = 'Basis'
            READ (MLINE(NMLINE),'(BN,1X,F9.0,25I5)') Q,NONT(I),
     &                           IQM(I,1), (JCO(J,I,1), J=1, IQM(I,1))
            IF (NINT(Q) .GT. 18) THEN
              NSETS = 1
              ADDSTO = .FALSE.
              WRITE (LUPRI,'(/A)') ' No Huckel start guess due to '//
     &             'nuclear charge larger than 18'
            END IF
          ENDIF
        ENDIF
C
C       Read/process data on symmetry independent centers
C       =================================================
C
        IF (ZMTTST) THEN
           CALL ZMAT(NONT,KATOM)
           Q = CHARGE(NUCIND)
        ELSE
           IF (INPTST) THEN
              WRITE(LUPRI,'(/A,I5)') 'Atomic type: ',I
              WRITE(LUPRI,'(3X,A,F4.1)') 'Nuclear charge:  ',Q
              WRITE(LUPRI,'(3X,A,I5)') 'Number of symmetry '//
     &             'independent centres: ',NONT(I)
           ENDIF
           CALL CNTINP(LUINFO,NONT(I),Q,QEXP,ANG)
C
C
C          Initialize reference strings for this atom type.
C
           BASREF(1,I,1) = ' Reference not found in input'
           DO ITYP = 1, 3
               DO IREF = 2, 10
                  BASREF(IREF,I,ITYP) = ' Not initialized'
               ENDDO
           ENDDO
C
C          Read the basis set. The input differs in DALTON/DIRAC due to
C          the optional reading of small component basis sets. The code 
C          switch starts here with first the processing of the DIRAC format.
C
           IF(DIRAC) THEN
C
C            Construct basis set for dummy centers if requested
C            ==================================================
C
             NQ = NINT(Q + 0.01D0)
             IF ((LCMMAX .GT. -1) .AND. (NQ .EQ. 0)) THEN
               WRITE(LUPRI,'(/A)') 'BASINP: CM basis not implemented !'
               CALL QUIT('BASINP: CM basis not implemented in Dirac !')
C
C            Get basis set from basis set library
C            ===================================
C
             ELSE IF (BASIS) THEN
               WRITE(LUPRI,'(/A/A)')
     &           'BASINP: "BASIS" in first line not allowed in Dirac.',
     &           'Use "LARGE BASIS <basis>" for each atom type instead.'
               CALL QUIT(
     &           'BASINP: "BASIS" in first line not allowed in Dirac!')
C
             ELSE
C
C            Read basis set info from file MOLECULE.INP
C            ==========================================
C
C              Start with LARGE Component
C              ==========================
C
C              The following seven possibilities exist:
C
C              LARGE INTGRL    2    1    1  (IBSFLAG = 1)
C              LARGE    2    1    1 (IBSFLAG = 1)
C              LARGE EXPLICIT    2    1    1 (IBSFLAG = 1)
C              LARGE MOLFBAS H.bas (IBSFLAG = 2)
C              LARGE BASIS basisset (IBSFLAG = 3)
C              LARGE GEOM alpha beta N 2 1 1 (IBSFLAG = 4)
C              LARGE EVENTEMP alpha beta N 2 1 1 (IBSFLAG = 4)
c               (alias for LARGE GEOM)
C              LARGE WELLTEMP alpha beta gamma delta N 2 1 1 (IBSFLAG = 4)
C              LARGE FAMILY    N    2    1    1 (IBSFLAG = 4)
C              LARGE DUALFAMILY    N    2    1    1 (IBSFLAG = 4)
C
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),'(A5,X,A14)') BSET,BSKEYWORD
C
               IF ( BSET .NE. 'LARGE') THEN
                  WRITE(LUPRI,'(A,I5,A,/,A)')
     $                 'Expected keyword LARGE, read line no. ',
     $                 NMLINE,':',MLINE(NMLINE)
                  GOTO 5100
               END IF
C
C              Initialize basis set flag
C
               IBSFLAG = 0
C
C              Check what kind of basis set we have
C
               IF ( BSKEYWORD(1:6) .EQ. 'INTGRL' ) THEN
                  IBSFLAG = 1
                  READ (MLINE(NMLINE),'(12X,12I5)')
     &                 IQM(I,1),(JCO(K,I,1),K=1,IQM(I,1))
               ELSE IF ( BSKEYWORD(1:8) .EQ. 'EXPLICIT') THEN
                  IBSFLAG = 1
                  READ (MLINE(NMLINE),'(14X,12I5)')
     &                 IQM(I,1),(JCO(K,I,1),K=1,IQM(I,1))
               ELSE IF ( BSKEYWORD(1:7) .EQ. 'MOLFBAS') THEN
                  IBSFLAG = 2
                  READ (MLINE(NMLINE),'(14X,A20)') BASFIL
               ELSE IF ( BSKEYWORD(1:5) .EQ. 'BASIS') THEN
                  IBSFLAG = 3
                  READ (MLINE(NMLINE),'(12X,A20)') BASFIL
               ELSE IF ( BSKEYWORD(1:8) .EQ. 'WELLTEMP' ) THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $                 (FAMPAR(K),K=1,4),NFAMEXP(1),
     $                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,0)
               ELSE IF ( BSKEYWORD(1:8) .EQ. 'EVENTEMP' ) THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $                 (FAMPAR(K),K=1,2),NFAMEXP(1),
     $                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,0)
               ELSE IF ( BSKEYWORD(1:4) .EQ. 'GEOM' ) THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE)(11:),*,ERR=9000)
     $                 (FAMPAR(K),K=1,2),NFAMEXP(1),
     $                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  FAMPAR(3) = D0
                  FAMPAR(4) = D0
                  CALL GENFAMEXP(LUINFO,0)
               ELSE IF ( BSKEYWORD(1:6) .EQ. 'FAMILY') THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE),'(12X,12I5)',ERR=9000)
     $                 NFAMEXP(1), IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,1)
               ELSE IF ( BSKEYWORD(1:10) .EQ. 'DUALFAMILY') THEN
                  IBSFLAG = 4
                  READ (MLINE(NMLINE),'(16X,12I5)',ERR=9000)
     $                 NFAMEXP(1), NFAMEXP(2),
     &                 IQM(I,1), (JCO(K,I,1),K=1,IQM(I,1))
                  CALL GENFAMEXP(LUINFO,2)
C
               ELSE IF ( BSKEYWORD(1:10) .EQ. 'MULTIBASIS') THEN
C
                  IBSFLAG = 5
                  READ (MLINE(NMLINE)(18:),*,ERR=9000)
     &                   NFILES,(BSFIL(K),K=1,NFILES)
                  IF ( NFILES .EQ. 0 ) THEN
                     WRITE(*,*) 
     &                 'PLEASE SPECIFY THE NUMBER OF BASIS SET FILES'
                     GOTO 9000
                  END IF
C
               ELSE
C                 Assume it is 'LARGE   2    1....'
                  IBSFLAG = 1
                  READ (MLINE(NMLINE),'(5X,12I5)')
     &                 IQM(I,1),(JCO(K,I,1),K=1,IQM(I,1))
               END IF
C
               IF ( IBSFLAG .EQ. 1 ) THEN
C
C                 Read explicitly given basis set.
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(1)
                  CALL GTOINP(LUINFO,IQM(I,1),JCO(1,I,1),
     &                 NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &                 SEG(IBLOCK(1),1),
     &                 ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &                 CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &                 NBLOCK,KAOVEC,KPRIM)
                  NBLCK(I,1) = NBLOCK
                  BASREF(1,I,1) = 
     &            "Basis set typed explicitly in input file "
C
               ELSE IF ( IBSFLAG .EQ. 2 ) THEN
C
C                 Molfdir basis set
C
                  KBFPR = 1
                  KBFCO = KBFPR + KPRIM
                  KLAST = KBFCO + KPRIM*KPRIM
                  IF (KLAST.GT.LWORK)
     &               CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
                  DIRCON(1) = .TRUE.
                  CALL BFGINP(I,IBLOCK(1),1,WORK(KBFPR),WORK(KBFCO),
     &                IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,BASFIL,
     &                BASREF(1,1,1))
               ELSE IF ( IBSFLAG .EQ. 3 ) THEN
C
                  BASNAM = BASFIL
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(1)

                  CALL BASLIB(IQM(I,1),JCO(1,I,1),NUC(IBLOCK(1),1),
     &                 NRC(IBLOCK(1),1),SEG(IBLOCK(1),1),
     &                 ALPHA(1,IBLOCK(1),1),
     &                 CPRIM(1,1,IBLOCK(1),1),CPRIMU(1,1,IBLOCK(1),1),
     &                 NBLOCK,ISGEN(IBLOCK(1)),KAOVEC,KPRIM,Q,1,DSM,
     &                 BASREF(1,I,1))
                  NBLCK(I,1) = NBLOCK
C
               ELSE IF ( IBSFLAG .EQ. 4 ) THEN
C
C                 Generate well-temperered/family basis set
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(1)
                  CALL FAMBAS(LUINFO,IQM(I,1),JCO(1,I,1),
     &                 NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &                 SEG(IBLOCK(1),1),
     &                 ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &                 CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &                 NBLOCK,KAOVEC,KPRIM)
                  NBLCK(I,1) = NBLOCK
                  BASREF(1,I,1) = 
     &   "Well-tempered basis set typed explicitly in input file "
C
               ELSE IF ( IBSFLAG .EQ. 5 ) THEN
C
C reading of basis sets from a number of files, including when the files
C do not necessarily have the same maximum and/or minimum orbital 
C angular momenta (last modif. 11/09/2003 by gomes) 
C
C we will move through the basis set files using the dummy variable iblock1
C 
                  IBLOCK1 = 1
                  MAXIQM  = 0
C
C we make sure we don't have any junk on the temporary variables 
C
                  CALL IZERO(JCOT,KNG*NFILES)
                  CALL IZERO(NUCT,KBLK)
                  CALL IZERO(NRCT,KBLK)
C
                  BASREF(1,I,1) =
     &            "Basis set read from multiple files"
C
                  DO 463 M=1,NFILES
                     BASNAM = BSFIL(M)
                     KAOVEC  = KBLOCK + 1 - IBLOCK1 
C
                     BASREF(3*(M-1)+2,I,1) = " " 
                     BASREF(3*(M-1)+3,I,1) = " ref. in file "//BSFIL(M) 

                     CALL BASLIB(IQMT,JCOT(1,M),
     &                   NUCT(IBLOCK1),
     &                   NRCT(IBLOCK1),SEGT(IBLOCK1),
     &                   ALPHAT(1,IBLOCK1),
     &                   CPRIMT(1,1,IBLOCK1),
     &                   CPRIMUT(1,1,IBLOCK1),
     &                   NBLOCK,ISGENT(IBLOCK1),KAOVEC,KPRIM,Q,1,DSM,
     &                   BASREFT(1,I,M))

                     BASREF(3*(M-1)+4,I,1) = BASREFT(1,I,M)
C
C to help in the reordering, we create a table mapping the orbital
C angular momentum of a block to the number of the block in the
C temp storage 
C 
                     IDBLK = IBLOCK1
                     IDANG = 1
                     DO 464 IRR=1,IQMT
                        DO 465 IRS=1,JCOT(IRR,M)
                           IBLKAGNMOM(IDBLK) = IDANG
                           IDBLK = IDBLK + 1
  465                   CONTINUE 
                        IDANG = IDANG + 1
  464                CONTINUE 
C
                     IF ( IQMT .GT. MAXIQM ) THEN
                        MAXIQM = IQMT
                     ENDIF
                     IBLOCK1 = IBLOCK1 + NBLOCK
  463             CONTINUE
C
C after reading all files, we reorder the data wrt orbital angular 
C momentum. in this reordering variable kt is the "final" index of a given block
C
                  IQM(I,1) = MAXIQM
                  KT = IBLOCK(1)
                  DO 466 IQ=1,MAXIQM
                  DO 467 IR=1,IBLOCK1 - 1
C
C in filling up the permanent storage with data, care should be taken to skip 
C zero-sized blocks that arise when we have only parts of a basis
C set in one file, as when the polarization functions are in one file and
C the scf set is in another...
C
                       IF ( IBLKAGNMOM(IR) .EQ. IQ ) THEN
                         IF ( NUCT(IR) .GT. 0 ) THEN

                              JCO(IQ,I,1)= JCO(IQ,I,1)+ 1
                              NUC(KT,1)  = NUC(KT,1) + NUCT(IR)
                              NRC(KT,1)  = NRC(KT,1) + NRCT(IR)
                              SEG(KT,1)  = SEGT(IR)
                              ISGEN(KT)  = ISGENT(IR)

                              DO 468 IX=1, NUCT(IR)
                                 ALPHA(IX,KT,1)      =  ALPHAT(IX,IR)
                                 DO 468 JX=1, NRCT(IR) 
                                  CPRIM(IX,JX,KT,1) =  CPRIMT(IX,JX,IR)
                                  CPRIMU(IX,JX,KT,1)= CPRIMUT(IX,JX,IR)
  468                         CONTINUE
                              KT = KT + 1
                           ENDIF 
                        ENDIF
  467                CONTINUE
  466             CONTINUE
C
C finally, we update nblck(i,1) with the total number of non-zero sized
C blocks read, subtracting one only because on the loops above kt was 
C incremented after the last block was put in place ...
C
                 NBLCK(I,1) = KT - 1 
C
C ... and we are done
C
               ENDIF
C
C              Print out basis set
C
               NOORBTS    = NOORBTS.AND.IQM(I,1).EQ.0
               IF (IQM(I,1).GT.0) THEN
                  IF(INPTST) THEN
                     WRITE(LUPRI,'(3X,A)') 'Large component set:'
                     WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                    'Max.ang.quantum no.:',(IQM(I,1)-1),
     &                    '  Blocks:',(JCO(K,I,1), K=1,IQM(I,1))
                  ENDIF
                  NHTYP     = MAX(NHTYP,IQM(I,1))
                  IF(NHTYP.GT.MXQN) GOTO 5000
               ENDIF
C
C              Check if large component basis is contracted
C              ============================================
C
               DO 210 JBLCK = 1, NBLCK(I,1)
                  IF ((NUC(JBLCK,1).NE.NRC(JBLCK,1)).AND.
     &                (NRC(JBLCK,1).GT.0)) CNTBAS = .TRUE.
  210          CONTINUE
               IF ( CNTBAS .AND.
     $              (IBSFLAG .EQ. 3 .OR. IBSFLAG .EQ. 1 )) THEN
                  DIRCON(1) = .TRUE.
               END IF
C
C
C
C              Next is the SMALL Component
C              ==========================
C
C
C              The following four possibilities exist:
C
C              SMALL INTGRL     2    1    1  (IBSFLAG = 1)
C              SMALL    2    1    1 (IBSFLAG = 1)
C              SMALL EXPLICIT     2    1    1  (IBSFLAG = 1)
C              SMALL MOLFBAS H.bas (IBSFLAG = 2)
C              SMALL KINBAL (IBSFLAG = 4)
C
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),'(A5,X,A7)') BSET,BSKEYWORD
C
C              Initialize basis set flag
C
               IBSFLAG = 0
C
C              Check whether the small component is specified in input
C
               IF ( BSET .NE. 'SMALL' ) THEN
C
C                    No 'SMALL' keyword given.
C                    Assume kinetic balance
C
C                    Discard the line we just read.
C
                     NMLINE = NMLINE - 1
                     IBSFLAG = 4
C
               ELSE
C
C
C                 Check what kind of basis set we have
C
                  IF ( BSKEYWORD(1:6) .EQ. 'INTGRL' ) THEN
                     IBSFLAG = 1
                     READ (MLINE(NMLINE),'(12X,12I5)')
     &                     IQM(I,2),(JCO(K,I,2),K=1,IQM(I,2))
                  ELSE IF ( BSKEYWORD(1:8) .EQ. 'EXPLICIT') THEN
                     IBSFLAG = 1
                     READ (MLINE(NMLINE),'(14X,12I5)')
     &                     IQM(I,2),(JCO(K,I,2),K=1,IQM(I,2))
                  ELSE IF ( BSKEYWORD .EQ. 'MOLFBAS') THEN
                     IBSFLAG = 2
                     READ (MLINE(NMLINE),'(14X,A20)') BASFIL
                  ELSE IF ( BSKEYWORD(1:6) .EQ. 'KINBAL' ) THEN
                     IBSFLAG = 4
                  ELSE
C                    Assume it is 'SMALL   2    1....'
                     IBSFLAG = 1
                     READ (MLINE(NMLINE),'(BN,A5,12I5)')
     &                    BSET,IQM(I,2),(JCO(K,I,2),K=1,IQM(I,2))
                  END IF
               END IF
C
C              We know the type, call the appropriate reader
C
               IF ( IBSFLAG .EQ. 1 ) THEN
C
C                 Read explicitly given small component basis set.
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(2)
                  CALL GTOINP(LUINFO,IQM(I,2),JCO(1,I,2),
     &                 NUC(IBLOCK(2),2),NRC(IBLOCK(2),2),
     &                 SEG(IBLOCK(2),2),
     &                 ALPHA(1,IBLOCK(2),2),CPRIM(1,1,IBLOCK(2),2),
     &                 CPRIMU(1,1,IBLOCK(2),2),ISGEN(IBLOCK(2)),
     &                 NBLOCK,KAOVEC,KPRIM)
                  NBLCK(I,2) = NBLOCK
                  BASREF(1,I,2) = 
     &            "Basis set typed explicitly in input file "
C
               ELSE IF ( IBSFLAG .EQ. 2 ) THEN
C
C                 Read molfdir small component basis set.
C
                  KBFPR = 1
                  KBFCO = KBFPR + KPRIM
                  KLAST = KBFCO + KPRIM*KPRIM
                  IF (KLAST.GT.LWORK)
     &                 CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
                  CALL BFGINP(I,IBLOCK(2),2,WORK(KBFPR),WORK(KBFCO),
     &                IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,BASFIL,
     &                BASREF(1,1,2))
C
               ELSE IF ( IBSFLAG .EQ. 4 ) THEN
C
C                 When uncontracted we generate the small component basis
C                 using the kinetic balance relation
C                 ...or...
C                 Make uncontracted small component basis
C                 set from kinetic balance.
C                 =======================================================
C
                  CALL KINBAL(I,IBLOCK(1),IBLOCK(2),
     &                 IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &                 NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,INPTST,
     &                 CNTBAS)
                  BASREF(1,I,2) = ' Derived from large component'
C
               ENDIF
               IBLOCK(1)  = IBLOCK(1) + NBLCK(I,1)
               IBLOCK(2)  = IBLOCK(2) + NBLCK(I,2)
               NOORBTS    = NOORBTS.AND.IQM(I,2).EQ.0
C
C              Print small comp. basis set
C
               IF (IQM(I,2).GT.0) THEN
                 IF(INPTST) THEN
                   WRITE(LUPRI,'(3X,A)') 'Small component set:'
                   WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                 'Max.ang.quantum no.:',(IQM(I,2)-1),
     &                 '  Blocks:',(JCO(K,I,2), K=1,IQM(I,2))
                 ENDIF
                 NHTYP     = MAX(NHTYP,IQM(I,2))
                 IF(NHTYP.GT.MXQN) GOTO 5000
               ENDIF
C
C              Check if small component basis is contracted
C              ============================================
C
               DO JBLCK = 1, NBLCK(I,2)
                  IF ((NUC(JBLCK,2).NE.NRC(JBLCK,2)).AND.
     &                (NRC(JBLCK,2).GT.0)) CNTBAS = .TRUE.
               END DO
               IF ( CNTBAS .AND. IBSFLAG .NE. 4 ) DIRCON(2) = .TRUE.
C
             ENDIF
           ELSE
C
C            Dalton (NOT DIRAC) case
C
             DO 300 J = 1,NSETS
C
C              Construct basis set for dummy centers if requested
C              ==================================================
C
               NQ = NINT(Q + 0.01D0)
               IF ((LCMMAX .GT. -1) .AND. (NQ .EQ. 0)) THEN
                  IQM(I,J) = LCMMAX + 1
                  KAOVEC = KBLOCK + 1 - IBLOCK(J)
                  CALL CMBAS(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                 NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                 ALPHA(1,IBLOCK(J),J),
     &                 CPRIM(1,1,IBLOCK(J),J),CPRIMU(1,1,IBLOCK(J),J),
     &                 NBLOCK,KAOVEC,KPRIM)
C
C              Get basis set from basis set library
C              ===================================
C
               ELSE IF ((BASIS .OR. (ADDSTO .AND. J .EQ. 2))) THEN
                  KAOVEC = KBLOCK + 1 - IBLOCK(2)
                  CALL BASLIB(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                  NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                  ALPHA(1,IBLOCK(J),J),
     &                  CPRIM(1,1,IBLOCK(J),J),CPRIMU(1,1,IBLOCK(J),J),
     &                  NBLOCK,ISGEN(IBLOCK(J)),KAOVEC,KPRIM,Q,J,DSM,
     &                  BASREF(1,I,1))
C
               ELSE
C
C              Read basis set from file
C              ========================
C
                  KAOVEC = KBLOCK + 1 - IBLOCK(J)
                  CALL GTOINP(LUINFO,IQM(I,J),JCO(1,I,J),
     &                NUC(IBLOCK(J),J),NRC(IBLOCK(J),J),
     &                SEG(IBLOCK(J),J),
     &                ALPHA(1,IBLOCK(J),J),CPRIM(1,1,IBLOCK(J),J),
     &                CPRIMU(1,1,IBLOCK(J),J),IDUMMY,
     &                NBLOCK,KAOVEC,KPRIM)
               ENDIF
               NBLCK(I,J) = NBLOCK
               IBLOCK(J)  = IBLOCK(J) + NBLCK(I,J)
               NOORBTS    = NOORBTS.AND.IQM(I,J).EQ.0
               IF (IQM(I,J).GT.0) THEN
                 IF(INPTST) THEN
                    WRITE(LUPRI,'(3X,A5,A)') BSET,' set:'
                    WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                   'Max.ang.quantum no.:',(IQM(I,J)-1),
     &                   '  Blocks:',(JCO(K,I,J), K=1,IQM(I,J))
                 ENDIF
C
C     Check angular momentum quantum number
C
                 NHTYP     = MAX(NHTYP,IQM(I,J))
                 IF(NHTYP.GT.MXQN) GOTO 5000
               ENDIF
 300         CONTINUE
           ENDIF
C
           DO 400 K = 1,NONT(I)
              NOORBT(K+INUC) = NOORBTS
 400       CONTINUE
C       End of switch DALTON/DIRAC
        END IF 
C
        DO ISET = 3, 4
C
C       Read basis sets for density fitting
C       ===================================
C
C       The following possibilities exist:
C
C       FTSET INTGRL     2    1    1  (IBSFLAG = 1)
C       FTSET     2    1    1 (IBSFLAG = 1)
C       FTSET DEFAULT (IBSFLAG = 2)
C       FTSET basisset (IBSFLAG = 3)
C       FTSET GEOM alpha beta N 2 1 1 (IBSFLAG = 4)
C       FTSET EVENTEMP alpha beta N 2 1 1 (IBSFLAG = 4)
c               (alias for FTSET GEOM)
C       FTSET WELLTEMP alpha beta gamma delta N 2 1 1 (IBSFLAG = 4)
C       FTSET FAMILY    N    2    1    1 (IBSFLAG = 4)
C       FTSET DUALFAMILY    N    2    1    1 (IBSFLAG = 4)
C 
C       We also distinguish between the keywords :
C       FTSET : fitset to be used for both large and small component densities
C       FTLRG : fitset to be used for only the large component density
C       FTSML : fitset to be used for only the small component density
C
        CALL RDLINE(LUINFO)
        READ (MLINE(NMLINE),'(A5,X,A14)') BSET,BSKEYWORD
C
C       Initialize basis set flag
C
        IBSFLAG = 0
C
C       Check whether the fit function set is specified in input
C                 
        IF ( BSET .NE. 'FTSET' ) THEN
C
C          No 'FTSET' keyword given.
C          Take standard set
C
C          Discard the line we just read.
C
           NMLINE = NMLINE - 1   
           IBSFLAG = 2     
C
        ELSE 
C
C          Check what kind of basis set we have
C
           IF ( BSKEYWORD(1:6) .EQ. 'INTGRL' ) THEN
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(12X,12I5)')
     &             IQM(I,ISET),(JCO(K,I,ISET),K=1,IQM(I,ISET))
           ELSE IF ( BSKEYWORD(1:8) .EQ. 'EXPLICIT') THEN
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(14X,12I5)')
     &             IQM(I,ISET),(JCO(K,I,ISET),K=1,IQM(I,ISET))
           ELSE IF ( BSKEYWORD(1:7) .EQ. 'DEFAULT') THEN
              IBSFLAG = 2
           ELSE IF ( BSKEYWORD(1:5) .EQ. 'BASIS') THEN
              IBSFLAG = 3
              READ (MLINE(NMLINE),'(12X,A20)') BASFIL
           ELSE IF ( BSKEYWORD(1:8) .EQ. 'WELLTEMP' ) THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $             (FAMPAR(K),K=1,4),NFAMEXP(1),
     $             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,0)
           ELSE IF ( BSKEYWORD(1:8) .EQ. 'EVENTEMP' ) THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE)(16:),*,ERR=9000)
     $             (FAMPAR(K),K=1,2),NFAMEXP(1),
     $             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,0)
           ELSE IF ( BSKEYWORD(1:4) .EQ. 'GEOM' ) THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE)(11:),*,ERR=9000)
     $             (FAMPAR(K),K=1,2),NFAMEXP(1),
     $             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              FAMPAR(3) = D0
              FAMPAR(4) = D0
              CALL GENFAMEXP(LUINFO,0)
           ELSE IF ( BSKEYWORD(1:6) .EQ. 'FAMILY') THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE),'(11X,12I5)',ERR=9000)
     $        NFAMEXP(1), IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,1)
           ELSE IF ( BSKEYWORD(1:10) .EQ. 'DUALFAMILY') THEN
              IBSFLAG = 4
              READ (MLINE(NMLINE),'(15X,12I5)',ERR=9000)
     $             NFAMEXP(1), NFAMEXP(2),
     &             IQM(I,ISET), (JCO(K,I,ISET),K=1,IQM(I,ISET))
              CALL GENFAMEXP(LUINFO,2)
C
           ELSE
C             Assume it is 'FTSET   2    1....'
              IBSFLAG = 1
              READ (MLINE(NMLINE),'(5X,12I5)')
     &             IQM(I,ISET),(JCO(K,I,ISET),K=1,IQM(I,ISET))
           END IF
C
           IF ( IBSFLAG .EQ. 1 ) THEN
C
C             Read explicitly given basis set.
C
              KAOVEC = KBLOCK + 1 - IBLOCK(ISET)
              CALL GTOINP(LUINFO,IQM(I,ISET),JCO(1,I,ISET),
     &          NUC(IBLOCK(ISET),ISET),NRC(IBLOCK(ISET),ISET),
     &          SEG(IBLOCK(ISET),ISET),
     &          ALPHA(1,IBLOCK(ISET),ISET),CPRIM(1,1,IBLOCK(ISET),ISET),
     &          CPRIMU(1,1,IBLOCK(ISET),ISET),ISGEN(IBLOCK(ISET)),
     &          NBLOCK,KAOVEC,KPRIM)
              NBLCK(I,ISET) = NBLOCK
              BASREF(1,I,ISET) = 
     &        "Fit set typed explicitly in input file   "
C
           ELSE IF ( IBSFLAG .EQ. 2 ) THEN
              BASREF(1,I,ISET) = ' Default fit function set    '
C
           ELSE IF ( IBSFLAG .EQ. 3 ) THEN
C
              BASNAM = BASFIL
C
              KAOVEC = KBLOCK + 1 - IBLOCK(ISET)
              CALL BASLIB(IQM(I,ISET),JCO(1,I,ISET),
     &             NUC(IBLOCK(ISET),ISET),
     &             NRC(IBLOCK(1),ISET),SEG(IBLOCK(ISET),ISET),
     &             ALPHA(1,IBLOCK(ISET),ISET),
     &             CPRIM(1,1,IBLOCK(ISET),ISET),
     &             CPRIMU(1,1,IBLOCK(ISET),ISET),
     &             NBLOCK,ISGEN(IBLOCK(ISET)),KAOVEC,KPRIM,Q,1,DSM,
     &             BASREF(1,I,ISET))
              NBLCK(I,ISET) = NBLOCK
C
           ELSE IF ( IBSFLAG .EQ. 4 ) THEN
C
C             Generate well-temperered/family basis set
C
              KAOVEC = KBLOCK + 1 - IBLOCK(ISET)
              CALL FAMBAS(LUINFO,IQM(I,ISET),JCO(1,I,ISET),
     &             NUC(IBLOCK(ISET),ISET),
     &             NRC(IBLOCK(ISET),ISET),SEG(IBLOCK(ISET),ISET),
     &             ALPHA(1,IBLOCK(ISET),ISET),
     &             CPRIM(1,1,IBLOCK(ISET),ISET),
     &             CPRIMU(1,1,IBLOCK(ISET),ISET),
     &             ISGEN(IBLOCK(ISET)),
     &             NBLOCK,KAOVEC,KPRIM)
              NBLCK(I,ISET) = NBLOCK
              BASREF(1,I,ISET) = 
     &   "Well-tempered fit set typed explicitly in input file "
C
           ENDIF
C
C          Print out fit set
C
           NOFITFS  = NOFITFS.AND.IQM(I,ISET).EQ.0
           IF (IQM(I,ISET).GT.0) THEN
              IF(INPTST) THEN
                 WRITE(LUPRI,'(3X,A)') 'Density fit set:'
                 WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                'Max.ang.quantum no.:',(IQM(I,ISET)-1),
     &                '  Blocks:',(JCO(K,I,ISET), K=1,IQM(I,ISET))
              ENDIF
              NHTYPF    = MAX(NHTYPF,IQM(I,ISET))
C             Since it's the density the angular momentum may be twice
C             the maximum value for the orbitals : 2 * MXQN - 2
              IF(NHTYPF-1.GT.2*MXQN-2) GOTO 5001
           ENDIF
C
C       (no need to check whether the set is contracted, so code differs
C        from analogous basis set read : this is not an error !)
C
        ENDIF
        IBLOCK(ISET)  = IBLOCK(ISET) + NBLCK(I,ISET)
C
C       End of loop over small/large fitsets.
        ENDDO
C
C     End of loop over unique atom types
  100 CONTINUE
C
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &     '*  Input specifies highest orbital of atomic type ',
     &     I,' as "',SPDCAR(NHTYP-1),'".',
     &     ' Highest allowed orbital in this version: ',
     &     SPDCAR(MXQN-1),
     &     ' Increase MXQN from',MXQN,' to',NHTYP/2,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 5001 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &     '*  Input specifies highest fit function of atomic type ',
     &     I,' as "',SPDCAR(NHTYPF-1),'".',
     &     ' Highest allowed fit function in this version: ',
     &     SPDCAR(2*MXQN-2),
     &     ' Increase MXQN from',MXQN,' to',NHTYPF,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 5100 CONTINUE
        CALL QUIT('Error reading keyword LARGE.')
 5200 CONTINUE
        CALL QUIT('Error reading keyword SMALL or SBFIL.')
 9000 CONTINUE
      WRITE(LUPRI,'(A,/,A,I3,/,A)')
     $     '*** ERROR in READIN ***',
     &     'Error reading line ',NMLINE,
     &     MLINE(NMLINE)
         CALL QUIT('*** ERROR in READIN ***')
      END
C  /* Deck cntinp */
      SUBROUTINE CNTINP(LUINFO,NONT,Q,QEXP,ANG)
C*****************************************************************************
C
C        Read and process data about symmetry independent centers
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "codata.h"
      PARAMETER(D0 = 0.0D00, CORMAX = 1.D5)
#if defined (SYS_HPUX)
      PARAMETER (LUERR = 7)
#else
      PARAMETER (LUERR = 0)
#endif
#include "cbirea.h"
C
      LOGICAL ANG
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IF (GAUNUC.AND.NINT(Q).NE.0) THEN
        IF(QEXP.EQ.D0) THEN
           CALL NUCSIZ(Q,GEXP)
        ELSE
           GEXP = QEXP
        ENDIF
C
C       Print out warning if GEXP < 1.0D8
C
        IF (GEXP .LE. 1.0D8) THEN
           WRITE(LUPRI,*) 'WARNING!!!'
           WRITE(LUPRI,*) 'WARNING!!!'
           WRITE(LUPRI,*) 'WARNING!!!'
           WRITE(LUPRI,*) 'Nuclear gaussian exponent small: ',GEXP
        END IF
        WRITE(LUPRI,'(/A,F8.3,A,1P,D20.10)')
     &       '  Nuclear Gaussian exponent for atom of charge',
     &       Q,' :',GEXP
      ELSE
         GEXP = D0
      END IF
      NUCSTR = NUCIND + 1
      IF (INPTST) WRITE(LUPRI,'(3X,A)') 'Symmetry independent centers:'
      DO 100 N = 1,NONT
         NUCIND = NUCIND + 1
         IF (NUCIND.GT.MXCENT) GOTO 5000
C
C*****************************************************************************
C        Card 7:
C           NAMN - name of symmetry independent center
C           CORR - coordinates of symmetry independent center
C*****************************************************************************
C
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(A4)') NAMN(NUCIND)
C        951115-hjaaj: NAMN() must never be read with free format
C          (otherwise it will on e.g. AIX be filled with nulls:
C           e.g. ' C  ' would become 'C\NULL\NULL\NULL'; and on
C           e.g. IRIX free format read of char. var. is not defined).
#if !defined (SYS_CRAY) && !defined (VAR_NOFREE) && !defined (SYS_T3D)
         READ (MLINE(NMLINE)(5:80),*,ERR=101) (CORD(J,NUCIND),J = 1, 3)
         GO TO 104
#endif
  101    READ (MLINE(NMLINE)(1:80),'(BN,4X,3F20.0)',ERR=102)
     &        (CORD(J,NUCIND), J = 1,3)
         GO TO 104
  102    READ (MLINE(NMLINE)(1:80),'(BN,4X,3F10.0)',ERR=103)
     &           (CORD(J,NUCIND), J = 1,3)
         GO TO 104
  103    CONTINUE
            WRITE(LUPRI,'(/A,I5/A,I5,A)')
     &      ' ERROR: Unable to read Cartesian coordinates of atom no.',
     &      NUCIND,' from line',NMLINE,' in the MOLECULE input file:'
            WRITE(LUPRI,'(A)') MLINE(NMLINE)
         CALL QUIT('ERROR reading atomic coordinates in MOLECULE input')
  104    CONTINUE
         IF(INPTST) THEN
            WRITE(LUPRI,'(6X,A4,3F20.15)') NAMN(NUCIND),
     &          (CORD(J,NUCIND), J = 1,3)
         ENDIF
         NCLINE(NUCIND) = NMLINE
         NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
         NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
         NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
         DO 200 J = 1,3
            IF(ANG) CORD(J,NUCIND) = CORD(J,NUCIND)/XTANG
            IF(ABS(CORD(J,NUCIND)).GT.CORMAX) GOTO 5010
  200   CONTINUE
C
C*****************************************************************************
C       CHARGE  - charge of center
C       NOORBT  - TRUE: no orbitals on this center
C       GNUEXP  - exponent of Gaussian nuclear charge distribution
C*****************************************************************************
C
         CHARGE(NUCIND) = Q
         GNUEXP(NUCIND) = GEXP
  100 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A/A,I5)')
     &    ' Too many atomic centers: MXCENT exceed in CNTINP,',
     &    ' Current limit:',MXCENT
        CALL QUIT('*** ERROR *** MXCENT exceeded in CNTINP')
 5010 CONTINUE
        WRITE (LUPRI,'(A,E12.6,A/A/A,E12.6)')
     &    ' Atomic coordinate ',CORD(J,NUCIND),
     &    ' too large in CNTINP.',
     &    ' Note: Program is unstable for large coordinates.',
     &    ' Maximum coordinate value:',CORMAX
        CALL QUIT('*** ERROR: Atomic coordinate too large in CNTINP')
      END
C  /* Deck baspro */
      SUBROUTINE BASPRO(WORK,LWORK,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,NUC,
     &                  NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &                  KBLOCK,KPRIM,DOOWN)
C
C*****************************************************************************
C
C       Process orbital and geometry input data
C
C*****************************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
#include "cbirea.h"
#include "dcbgen.h"
      DIMENSION WORK(LWORK)
      LOGICAL SEG,DOOWN
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &     NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &     SEG(KBLOCK,KSETS),
     &     ALPHA(KPRIM,KBLOCK,KSETS),CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &     CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS)
      CHARACTER*1 SPDCAR
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "primit.h"
#include "huckel.h"
#include "aosotr.h"
#include "gencon.h"
#include "ibtfun.h"
C
C*****************************************************************************
C       MULK    - bitstring of basic operations that stabilise a center
C*****************************************************************************
C
      II = 1
      DO 100 I = 1, NONTYP
         DO 110 J = 1, NONT(I)
            MULK = 0
            LL = 1
            DO 140 L = 1, NSYMOP
               DO 150 M = 1,3
                  IF(IBTAND(LL,ISYMAX(M,1)) .NE. 0) THEN
C                 Basic operation LL moves coordinate M
                     IF(ABS(CORD(M,II)).GE.1.0D-06) GOTO 140
                  END IF
 150           CONTINUE
               MULK = MULK + LL
 140        LL = 2*LL
            ISTBNU(II) = MULK
            II = II + 1
 110     CONTINUE
 100  CONTINUE
C*****************************************************************************
C       Set various basis parameters
C*****************************************************************************
      CALL BASPAR(NSYMOP)
C*****************************************************************************
C       Process nuclear data
C*****************************************************************************
      CALL NUCPRO(WORK,LWORK)
C*****************************************************************************
C       Process orbital data for large components
C*****************************************************************************
      LCOMP  = 1
      ISHELL = 0
      IPRIM  = 0
      IPRIMD = 0
      IORB   = 0
      IORBD  = 0
      CALL ORBPRO(NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),NUC(1,1),
     &            NRC(1,1),SEG(1,1),
     &            ALPHA(1,1,1),CPRIM(1,1,1,1),CPRIMU(1,1,1,1),
     &            KATOM,KANG,KBLOCK,KPRIM,ISHELL,
     &            IPRIM,IPRIMD,IORB,IORBD,LCOMP)
      NLRGSH = ISHELL
      NLARGE = IORBD
      NORBL  = IORB
      NPLRG  = IPRIMD
      NPLSH  = IPRIM
CTROND      NPLRG  = IPRIM
C*****************************************************************************
C       Process orbital data for small components
C*****************************************************************************
      LCOMP = 2
      CALL ORBPRO(NONTYP,NONT,IQM(1,2),NBLCK(1,2),JCO(1,1,2),NUC(1,2),
     &            NRC(1,2),SEG(1,2),
     &            ALPHA(1,1,2),CPRIM(1,1,1,2),CPRIMU(1,1,1,2),
     &            KATOM,KANG,KBLOCK,KPRIM,
     &            ISHELL,IPRIM,IPRIMD,IORB,IORBD,LCOMP)
      IF (ISHELL.GT. MXSHEL) GOTO 5000
      IF (IORBD .GT. MXCORB) GOTO 5010
      IF (IPRIM .GT. MXPRIM) GOTO 5020
      IF (ADDSTO) THEN
         NSMLSH = ISHELL - NLRGSH
         NORBS  = IORB
         NSMALL = IORBD - NLARGE
         NPSML  = IPRIMD - NPLRG
         NPSHEL = IPRIM - NPLRG
C
         KMAX   = NLRGSH
         NBASIS = NLARGE
         NPBAS  = NPLRG
C
         NHUCSH = NSMLSH
         NHUCBA = NSMALL
         NHUCNP = NPSML
C
      ELSE
         KMAX   = ISHELL
         NSMLSH = KMAX - NLRGSH
         NORBS  = IORB
         NBASIS = IORBD
         IF (IPREAD .GT. 20) THEN
            WRITE(LUPRI,'(A)') 'PRIEXP:'
            WRITE(LUPRI,'(I5,F15.5)') (KK,PRIEXP(KK),KK=1,IPRIM)
         END IF
         NSMALL = NBASIS - NLARGE
         NPBAS  = IPRIMD
         NPSML   = IPRIMD - NPLRG
         NPSHEL  = IPRIM
         NPSSH   = IPRIM - NPLSH
      END IF
C*****************************************************************************
C       Process fit set data for first (LL or LL+SS) density fitset
C*****************************************************************************
      LCOMP_FIT  = 1
      ISHELL_FIT = 0
      IPRIM_FIT  = 0
      IPRIMD_FIT = 0
      IORB_FIT   = 0
      IORBD_FIT  = 0
      CALL FITPRO(NONTYP,NONT,IQM(1,3),NBLCK(1,3),JCO(1,1,3),NUC(1,3),
     &            NRC(1,3),SEG(1,3),
     &            ALPHA(1,1,3),CPRIM(1,1,1,3),CPRIMU(1,1,1,3),
     &            KATOM,KANG,KBLOCK,KPRIM,
     &            ISHELL_FIT,IPRIM_FIT,IPRIMD_FIT,
     &            IORB_FIT,IORBD_FIT,LCOMP_FIT)
      NLRGSH_FIT = ISHELL_FIT
      NLARGE_FIT = IORBD_FIT
      NORBL_FIT  = IORB_FIT
      NPLRG_FIT  = IPRIMD_FIT
      NPLSH_FIT  = IPRIM_FIT
C*****************************************************************************
C       Process fit set data for secundary (SS) density fitset
C*****************************************************************************
      LCOMP_FIT = 2
      CALL FITPRO(NONTYP,NONT,IQM(1,4),NBLCK(1,4),JCO(1,1,4),NUC(1,4),
     &            NRC(1,4),SEG(1,4),
     &            ALPHA(1,1,4),CPRIM(1,1,1,4),CPRIMU(1,1,1,4),
     &            KATOM,KANG,KBLOCK,KPRIM,
     &            ISHELL_FIT,IPRIM_FIT,IPRIMD_FIT,
     &            IORB_FIT,IORBD_FIT,LCOMP_FIT)
      IF (ISHELL_FIT.GT. MXSHEL) GOTO 5000
      IF (IORBD_FIT .GT. MXCORB) GOTO 5010
      IF (IPRIM_FIT .GT. MXPRIM) GOTO 5020
      KMAX_FIT   = ISHELL_FIT
      NSMLSH_FIT = KMAX_FIT - NLRGSH_FIT
      NBASIS_FIT = IORBD_FIT
      NSMALL_FIT = NBASIS_FIT - NLARGE_FIT
      NPBAS_FIT  = IPRIMD_FIT
      NPSML_FIT   = IPRIMD_FIT - NPLRG_FIT
      NPSHEL_FIT  = IPRIM_FIT
      NPSSH_FIT   = IPRIM_FIT - NPLSH_FIT
C
C     ***************************************
C     ***** Orbital Symmetry Processing *****
C     ***************************************
C
      KKVAL = 1
      KMVAL = KKVAL + MXAQN
      KNVAL = KMVAL + MXAQN
      KIREP = KNVAL + MXAQN
      KLAST = KIREP + MXCORB
      IF (KLAST .GT. LWORK) GOTO 5030
      CALL SYMPRO(WORK(KKVAL),WORK(KMVAL),WORK(KNVAL),WORK(KIREP),
     &            DOOWN)
C
C FIXME: In principle kinetic balance can be made for contracted
C FIXME: basis sets iff the small component basis set exponents
C FIXME: are kinetically balanced. Currently, the code resets
C FIXME: ISPHTR if the small component basis set is explicitly given
C FIXME: (IBSFLAG .NE. 4, see INFO: above)./hjaaj+jth sep 2000
C
      IF (DIRAC.AND.ISPHTR.EQ.3)
     $     CALL RKBREL(IPREAD)
C
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, number of contracted shells      ',ISHELL,
     *  '                   current maximum number (MXSHEL)  ',MXSHEL
        CALL QUIT('BASPRO: Too many contracted shells')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, number of contracted basis functions',IORBD,
     *  '                   current maximum number (MXCORB)     ',MXCORB
        CALL QUIT('BASPRO: Too many contracted basis functions')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, number of primitive basis shells   ',IPRIM,
     *  '                   current maximum number (MXPRIM)    ',MXPRIM
        CALL QUIT('BASPRO: Too many primitive shells')
 5030 CONTINUE
        CALL STOPIT('BASPRO','SYMPRO',KLAST,LWORK)
      END
C  /* Deck basout */
      SUBROUTINE BASOUT(WORK,LWORK,NONTYP,NONT,IQM,NBLCK,JCO,BASREF,
     &                  NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &                  KBLOCK,KPRIM,HERMIT)
C*****************************************************************************
C
C     Print orbital and geometry input data
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
#include "cbirea.h"
#include "dcbgen.h"
      DIMENSION WORK(LWORK)
      LOGICAL SEG,HERMIT,NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &     NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &     SEG(KBLOCK,KSETS),
     &     ALPHA(KPRIM,KBLOCK,KSETS),CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &     CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS)
      CHARACTER SPDCAR
      CHARACTER*(6*MXQN+3) TSTRING(MXBSETS)
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      CHARACTER*31 SET_TYPE(MXBSETS)
      DIMENSION NPRIM(MXBSETS),NBAS(MXBSETS)
      DIMENSION IBLOCK(MXBSETS)
#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "primit.h"
#include "symmet.h"
#include "aosotr.h"
#include "ibtfun.h"
C
      IF(DIRAC) THEN
         SET_TYPE(1) = 'large component basis functions'
         SET_TYPE(2) = 'small component basis functions'
      ELSE
         SET_TYPE(1) = 'basis functions                '
         SET_TYPE(2) = 'auxiliary basis functions      '
      ENDIF
      SET_TYPE(3) = 'primary density fit functions  '
      SET_TYPE(4) = 'secundary density fit functions'
C
C     ******************************************
C     ***** Printing of basis information  *****
C     ******************************************
C
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Atoms and basis sets',1)
         WRITE (LUPRI,'(1X,A,I5  )') ' Number of atom types: ',NONTYP
         WRITE (LUPRI,'(1X,A,I5,/)') ' Total number of atoms:',NUCDEP
      END IF
C
      DO ISET = 1, KSETS
         IBLOCK(ISET) = 1
      ENDDO
      ICENT  = 0
      IF (IPREAD .GT. 0) THEN
         WRITE (LUPRI,'(2X,A)')
     &        'label    atoms   charge   prim    cont     basis   '
         WRITE (LUPRI,'(2X,70A1)') ('-',I=1,70)
      END IF
      NCHTOT = 0
      DO 100 I = 1,NONTYP
C
C        Convert information to a string that can be understood by the innocent user
C
         DO ISET = 1, KSETS
            TSTRING(ISET) = ' '
            IF (IQM(I,ISET).GT.0) THEN
               CALL BASTYP(IQM(I,ISET),JCO(1,I,ISET),
     &                     NRC(IBLOCK(ISET),ISET),
     &                     NUC(IBLOCK(ISET),ISET),
     &                     NPRIM(ISET),NBAS(ISET),TSTRING(ISET))
            ELSE
               TSTRING(ISET) = 'No '//SET_TYPE(ISET)//' attached'
               NPRIM(ISET) = 0
               NBAS(ISET) = 0
            END IF
            IBLOCK(ISET) = IBLOCK(ISET) + NBLCK(I,ISET)
         END DO
C
         DO 20 N = 1,NONT(I)
            ICENT  = ICENT + 1
            NCHARG = NINT(CHARGE(ICENT))
            NCHTOT = NCHTOT + NUCDEG(ICENT)*NCHARG
            IF (IPREAD .GT. 0) THEN
               IF(DIRAC) THEN
                  WRITE (LUPRI,'(2X,A4,1X,4I8,6X,2A)') NAMN(ICENT),
     &                   NUCDEG(ICENT),NCHARG,NPRIM(1),NBAS(1),
     &                   'L  - ',TSTRING(1)
                  WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &                  NPRIM(2),NBAS(2),'S  - ',TSTRING(2)
                  IF (NPRIM(3).GT.0)
     &               WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &                  NPRIM(3),NBAS(3),'FL - ',TSTRING(3)
                  IF (NPRIM(4).GT.0)
     &               WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &                  NPRIM(4),NBAS(4),'FS - ',TSTRING(4)
               ELSE
                  WRITE (LUPRI,'(2X,A4,1X,4I8,6X,A)') NAMN(ICENT),
     &                 NUCDEG(ICENT),NCHARG,NPRIM(1),NBAS(1),TSTRING(1)
                  WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &                  NPRIM(3),NBAS(2),'Aux - ',TSTRING(2)
                  WRITE(LUPRI,'(23X,2I8,6X,2A)')
     &                  NPRIM(3),NBAS(3),'Fit - ',TSTRING(3)
               ENDIF
            END IF
   20    CONTINUE
  100 CONTINUE
      IF (IPREAD .GT. 0) THEN
         IF(DIRAC) THEN
            CALL PRSYMB(LUPRI,'-',70,2)
            WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPLRG,NLARGE,'L - large components'
            WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPSML,NSMALL,'S - small components'
         ENDIF
         CALL PRSYMB(LUPRI,'-',70,2)
         WRITE (LUPRI,'(2X,A,I7,3I8,6X,/)')
     &         'total:',NUCDEP,NCHTOT,NPBAS,NBASIS
         IF (NHTYP .GT. 2) THEN
            IF(DOCART) THEN
               WRITE (LUPRI,'(2X,A)') 'Cartesian basis used.'
            ELSE
               WRITE (LUPRI,'(2X,A)') 'Spherical harmonic basis used.'
            ENDIF
         END IF
         WRITE (LUPRI,'(1X,A,1P,D10.2)') ' Threshold for integrals:',
     &                                   THRS
      END IF
C
C     Print the references to the basis sets

      IF (IPREAD .GT. 0) CALL PRBASREF(LUPRI,NONTYP,BASREF)
C
C     ******************************************
C     ***** Printing of atomic coordinates *****
C     ******************************************
C
      CALL PRICAR(WORK,LWORK)
C
 1010 FORMAT(11X,I3,' large AO-blocks of ',A1,' type')
 1020 FORMAT(11X,I3,' small AO-blocks of ',A1,' type')
      IF (IPREAD .GE. 2) THEN
         KIPCON = 1
         KLAST  = KIPCON + KPRIM
         IF (KLAST.GT.LWORK) CALL STOPIT('BASOUT','ORBOUT',KLAST,LWORK)
C*****************************************************************************
C        Print information about large component basis
C*****************************************************************************
         IF(NLARGE.GT.0) THEN
            IF (DIRAC) CALL HEADER('Large Components Basis',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(1,NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &                  NUC(1,1),NRC(1,1),SEG(1,1),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  CPRIMU,NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
               CALL ORBOUT(1,NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &                     NUC(1,1),NRC(1,1),SEG(1,1),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     CPRIMU,NRMPRI)
            END IF
         END IF
C*****************************************************************************
C        Print information about small component basis
C*****************************************************************************
         IF (NBASIS.GT.NLARGE) THEN
            CALL HEADER('Small Components Basis',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(2,NONTYP,NONT,IQM(1,2),NBLCK(1,2),JCO(1,1,2),
     &                  NUC(1,2),NRC(1,2),SEG(1,2),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  CPRIMU(1,1,1,2),NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
               CALL ORBOUT(2,NONTYP,NONT,IQM(1,2),NBLCK(1,2),JCO(1,1,2),
     &                     NUC(1,2),NRC(1,2),SEG(1,2),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     CPRIMU(1,1,1,2),NRMPRI)
            END IF
         END IF
      END IF
C
C     Symmetry basis:
C
      IF (MAXREP.GT.0.AND.IPREAD.GT.2) CALL SYMOUT
C
C     Symmetrized magnetic moments
C
CTROND      IF (HERMIT) CALL MAGCOR(IPREAD)
      CALL MAGCOR(IPREAD)
C
C     ***************************
C     ***** Determine IAOAO *****
C     ***************************
C
      CALL AOTOAO(IAOAO,IPREAD)
C
C     **************************
C     **** Determine ICNTAO ****
C     **************************
C
      CALL CNTAO(IPREAD)
      END
C  /* Deck bastyp */
      SUBROUTINE BASTYP(IQM,JCO,NRC,NUC,NPRIM,NBAS,TSTRING)
C*****************************************************************************
C
C     Generate string : Uncontracted/contracted basis functions
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
      CHARACTER TSTRING*(6*MXAQN+3),SPDCAR*1
      DIMENSION JCO(IQM),NRC(*),NUC(*)
#include "ccom.h"
#include "gencon.h"
C*****************************************************************************
      IND = 1
      TSTRING(IND:IND) = '['
C*****************************************************************************
C * P R I M I T I V E S:
C*****************************************************************************
      IBLOCK = 0
      NPRIM = 0
      DO 100 I = 1,IQM
         IPRIM  = 0
         ISHELL = 0
         NCOMP = KHK(I)
C
         DO 200 J = 1,JCO(I)
            IBLOCK = IBLOCK + 1
            IPRIM  = IPRIM  + NUC(IBLOCK)
            NPRIM = NPRIM + NCOMP*NUC(IBLOCK)
  200    CONTINUE
C
         IF (IPRIM.GT.0) THEN
            IF (IPRIM.GE.10) THEN
               NDIG = INT(LOG10(FLOAT(IPRIM))+1)
               ITEN = 10**NDIG
               DO 300 K = NDIG,1,-1
                  ITEN = ITEN/10
                  IND = IND + 1
                  IDIG = ICHAR('0') + INT(IPRIM/ITEN)
                  TSTRING(IND:IND) = CHAR(IDIG)
                  IPRIM=MOD(IPRIM,ITEN)
  300          CONTINUE
            ELSE
               IND = IND + 1
               IDIG = ICHAR('0') + IPRIM
               TSTRING(IND:IND) = CHAR(IDIG)
            ENDIF
            IND = IND + 1
            TSTRING(IND:IND) = SPDCAR(I-1)
         END IF
  100 CONTINUE
C*****************************************************************************
      IND = IND + 1
      TSTRING(IND:IND) = '|'
C*****************************************************************************
C * S H E L L S:
C*****************************************************************************
      IBLOCK = 0
      NBAS = 0
      DO 400 I = 1,IQM
         ISHELL = 0
         NCOMP = KHK(I)
         DO 500 J = 1,JCO(I)
            IBLOCK = IBLOCK + 1
            ISHELL  = ISHELL  + NRC(IBLOCK)
            NBAS = NBAS + NCOMP*NRC(IBLOCK)
  500    CONTINUE
         IF(ISHELL.GT.0) THEN
           IF (ISHELL.GE.10) THEN
              NDIG = INT(LOG10(FLOAT(ISHELL))+1)
              ITEN = 10**NDIG
              DO 600 K = NDIG,1,-1
                 ITEN = ITEN/10
                 IND = IND + 1
                 IDIG = ICHAR('0') + INT(ISHELL/ITEN)
                 TSTRING(IND:IND) = CHAR(IDIG)
                 ISHELL=MOD(ISHELL,ITEN)
  600         CONTINUE
           ELSE
              IND = IND + 1
              IDIG = ICHAR('0') + ISHELL
              TSTRING(IND:IND) = CHAR(IDIG)
           ENDIF
           IND = IND + 1
           TSTRING(IND:IND) = SPDCAR(I-1)
         ENDIF
  400 CONTINUE
C*****************************************************************************
      IND = IND + 1
      TSTRING(IND:IND) = ']'
      RETURN
      END
C  /* Deck gtoinp */
      SUBROUTINE GTOINP(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,ISGEN,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 SPDCAR, FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
C
C
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),ISGEN(KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
#include "dcbgen.h"
#include "gencon.h"
C
      IF(DIRAC) THEN
C
C     DIRAC allows user to split blocks !
C
        NBLOCK = 0
        DO 100 I = 1, IQM
          DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
            CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
            CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C           ISGEN - bit packed information on generation of small
C                   components using the kinetic balance relation
C             1 - downwards: small L+1 from large L
C             2 - upwards  : small L-1 from large L
C             0 - if not specified same as 3 (both down and up)
C
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,A,I4,2I5)')
     &              FRMT,NUCIJ,NRCIJ,ISGEN(NBLOCK)
            IF (INPTST) THEN
              WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
              WRITE(LUPRI,'(6X,A,I5)')
     $             'Primitives:  ',NUCIJ,
     &             'Shells    :  ',NRCIJ,
     $             'KinBal    :  ',ISGEN(NBLOCK)
            END IF
            SEGIJ = NRCIJ .GE. 0
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C     Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &           CPRIMU(1,1,NBLOCK),KPRIM)
C
C           Identify segmented contractions
C
            CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
            SEG(NBLOCK) = SEGIJ
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &             SEG(NBLOCK),KPRIM,DSM)
C
C            Normalize orbitals
C
             CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),
     &                   CPRIM(1,1,NBLOCK),
     &                   CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
            NUC(NBLOCK) = NUCIJ
            NRC(NBLOCK) = NRCIJ
  200      CONTINUE
  100   CONTINUE
      ELSE
      NBLOCK = 0
      DO 300 I = 1, IQM
         NBLOCK = NBLOCK + 1
         IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
         CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
         NUCI = 0
         NRCI = 0
         SEGIJ  = .TRUE.
         DO 400 J = 1, JCO(I)
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C
            CALL RDLINE(LUINFO)
            IF (J. EQ. 1) THEN
               READ (MLINE(NMLINE),'(BN,A,I4,1I5)')
     &              FRMT,NUCIJ,NRCIJ
            ELSE
               READ (MLINE(NMLINE),'(BN,A,I4,1I5)') FRMT,NUCIJ,NRCIJ
            END IF
            IF (INPTST) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
               WRITE(LUPRI,'(6X,A,I5)')
     $              'Primitives:  ',NUCIJ,
     &              'Shells    :  ',NRCIJ
            END IF
            SEGIJ = SEGIJ .AND. (NRCIJ .GE. 0)
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCI + NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCI + NRCIJ.GT.KPRIM) GOTO 5020
C
C           Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(NUCI+1,NBLOCK),
     &                  CPRIMU(NUCI+1,NRCI+1,NBLOCK),KPRIM)
C
            NUCI = NUCI + NUCIJ
            NRCI = NRCI + NRCIJ
  400    CONTINUE
         NUC(NBLOCK) = NUCI
         NRC(NBLOCK) = NRCI
C
         JCO(I) = 1
C
C        Identify segmented contractions
C
         CALL SEGORB(SEGIJ,NUCI,NRCI,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
         SEG(NBLOCK) = SEGIJ
C
C        Reorder primitive orbitals
C
         CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCI,NRCI,
     &               SEG(NBLOCK),KPRIM,DSM)
C
C        Normalize orbitals
C
         CALL NRMORB(I,NRCI,NUCI,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &               CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  300 CONTINUE
      ENDIF
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> GTOINP error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error,'//
     *  'number of primitives per block      ',NUCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',
     *  NRCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
C  /* Deck gtoinp */
      SUBROUTINE FAMBAS(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,ISGEN,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 SPDCAR, FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
C
C
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),ISGEN(KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
#include "dcbgen.h"
#include "gencon.h"
C
      NBLOCK = 0
      DO 100 I = 1, IQM
         DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
            CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
            CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
            ISGEN (NBLOCK) = 0
C
C           Get exponents
C
            IC = MOD(I + 1, 2) + 1
            CALL FAMCOPY(LUINFO,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &           CPRIMU(1,1,NBLOCK),KPRIM,IC)
C
            IF (INPTST) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
               WRITE(LUPRI,'(6X,A,I5)')
     $              'Primitives:  ',NUCIJ,
     &              'Shells    :  ',NRCIJ,
     $              'KinBal    :  ',ISGEN(NBLOCK)
            END IF
            SEGIJ = NRCIJ .GE. 0
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C
C           Identify segmented contractions
C
            CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
            SEG(NBLOCK) = SEGIJ
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &           SEG(NBLOCK),KPRIM,DSM)
C
C           Normalize orbitals
C
            CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),
     &           CPRIM(1,1,NBLOCK),
     &           CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
            NUC(NBLOCK) = NUCIJ
            NRC(NBLOCK) = NRCIJ
  200      CONTINUE
  100   CONTINUE
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6,A)')
     *  ' >>> FAMBAS error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> FAMBAS error,'//
     *  'number of primitives per block      ',NUCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> FAMBAS error, number of contracted functions      ',
     *  NRCIJ,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
C  /* Deck famcopy */
      SUBROUTINE FAMCOPY(LUINFO,NUCIJ,NRCIJ,ALPHA,CPRIMU,KPRIM,IC)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D1 = 1.0D0)
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
#include "cbirea.h"
      DIMENSION ILIST(MXFAMEXP),TMP(MXFAMEXP+2)
#include "molinp.h"
C
C     Card 9:
C
C     ALPHA  - exponent of primitive
C     CPRIMU - unnormalized contraction coefficients
C
C
C     NUMLST needs an work array for dimension NFALPHA.
C
C     Read string in format as described in SUBR NUMLST
C     Eg.: 1..20
C
      CALL RDLINE(LUINFO)
      NUCIJ = 1
      CALL NUMLST(MLINE(NMLINE)(1:72),ILIST,NFAMEXP(IC),1,NFAMEXP(IC),
     &     1,NUCIJ,TMP,MXFAMEXP)
      NRCIJ = NUCIJ
      DO 300 L = 1, NUCIJ
         ALPHA(L) = FAMEXP(ILIST(L), IC)
         CPRIMU(L,L) = D1
  300 CONTINUE
C
C     Test print
C
      IF (INPTST) THEN
         WRITE(LUPRI,'(A)') 'Contraction matrix:'
         DO 400 L = 1, NUCIJ
            WRITE (LUPRI,'(E16.10)') ALPHA(L)
            WRITE (LUPRI,'(6(2X,F10.8))') (CPRIMU(L,M),M=1,NRCIJ)
  400    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck cmbas */
      SUBROUTINE CMBAS(IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC),
     &          BASA(0:4), BASB(0:4)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
      DATA (BASA(I),I=0,4) /0.584342, 0.452615, 0.382362, 0.337027,
     &                      0.304679/
      DATA (BASB(I),I=0,4) /0.424483, 0.309805, 0.251333, 0.215013,
     &                      0.189944/
C
      NBLOCK = 0
      DO 100 I = 1, IQM
         NBLOCK = NBLOCK + 1
         IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
         CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
         NUCI = 0
         NRCI = 0
         SEGIJ  = .TRUE.
         DO 110 J = NCMSTR, NCMEND
            NUCI = NUCI + 1
            NRCI = NRCI + 1
            IF (NUCI.GT.KPRIM) GOTO 5010
            IF (NRCI.GT.KPRIM) GOTO 5020
            ALPHA(NUCI,NBLOCK) = ((ZCMVAL/J)/
     &                           ((BASA(I-1)*J)/2+BASB(I-1)))**2
            CPRIMU(NUCI,NRCI,NBLOCK) = 1.0D0
  110    CONTINUE
         NUC(NBLOCK) = NUCI
         NRC(NBLOCK) = NRCI
C
         JCO(I) = 1
         SEG(NBLOCK) = SEGIJ
C
C        Normalize orbitals
C
         CALL NRMORB(I,NRCI,NUCI,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &               CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  100 CONTINUE
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> CMBAS  error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error,'//
     *  'number of primitives per block      ',NUCI,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error, number of contracted functions      ',
     *  NRCI,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
#ifdef NOT_USED
C old version of GTOINP ??? -hjaaj July 2000
C  /* Deck orbinx */
      SUBROUTINE ORBINX(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 SPDCAR, FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
      LOGICAL SEG, FREE
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
C
      NBLOCK = 0
      DO 100 I=1,IQM
         DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,A,I4,I5)')
     &           FRMT,NUC(NBLOCK),NRC(NBLOCK)
            IF (INPTST) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
               WRITE(LUPRI,'(6X,A,I5)') 'Primitives:  ',NUC(NBLOCK),
     &                                  'Shells    :  ',NRC(NBLOCK)
            END IF
            SEG(NBLOCK) = NRC(NBLOCK).GE.0
            NUC(NBLOCK) = ABS(NUC(NBLOCK))
            NRC(NBLOCK) = ABS(NRC(NBLOCK))
            NUCIJ = NUC(NBLOCK)
            NRCIJ = NRC(NBLOCK)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C           Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),KPRIM)
C
C           Identify segmented contractions
C
            CALL SEGORB(SEG(NBLOCK),NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),
     &                  KPRIM,DSM)
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &                  SEG(NBLOCK),KPRIM,DSM)
C
C           Normalize orbitals
C
            CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  200    CONTINUE
  100 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> GTOINP error, number of AO-blocks            ',NBLOCK,
     *  '                   number of blocks not used yet  ',KAOVEC,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error,'//
     *  'number of primitives per block      ',NUC(NBLOCK),
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',
     *  NRC(NBLOCK),
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
#endif /* NOT_USED */
C  /* Deck acporb */
      SUBROUTINE ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA,CPRIMU,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0)
      CHARACTER*1 FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
      CHARACTER*7 FMT3
      LOGICAL FREE
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
#include "cbirea.h"
#include "molinp.h"
C
C     Card 9:
C
C     ALPHA  - exponent of primitive
C     CPRIMU - unnormalized contraction coefficients
C
      IF (FRMT .EQ. 'F' .OR. FRMT .EQ. 'f') THEN
         FREE = .TRUE.
      ELSE IF (FRMT .EQ. 'H' .OR. FRMT .EQ. 'h') THEN
         FREE = .FALSE.
         NCOL = 4
         FMT1 = '(4F20.8)'
         FMT2 = '(20X,3F20.8)'
         FMT3 = '(F20.8)'
      ELSE
         FREE = .FALSE.
         NCOL = 8
         FMT1 = '(8F10.4)'
         FMT2 = '(10X,7F10.4)'
         FMT3 = '(F10.4)'
      END IF
      IF(NRCIJ.NE.0) THEN
        DO 100 L = 1, NUCIJ
           CALL RDLINE(LUINFO)
           IF (FREE) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
              ISTART = 1
              CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALPHA(L),'REA')
              DO M = 1, NRCIJ
                 CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,
     &                       CPRIMU(L,M),'REA')
              END DO
#else
              READ (MLINE(NMLINE),*) ALPHA(L), (CPRIMU(L,M),M=1,NRCIJ)
#endif
           ELSE
              READ (MLINE(NMLINE),FMT1) ALPHA(L),
     &           (CPRIMU(L,M), M = 1, MIN(NRCIJ,NCOL - 1))
              DO 200 K = 2, (NRCIJ - 1)/(NCOL - 1) + 1
                 CALL RDLINE(LUINFO)
                 READ (MLINE(NMLINE),FMT2) (CPRIMU(L,M),
     &              M = 1 + (NCOL-1)*(K-1), MIN(NRCIJ,(NCOL-1)*K))
  200         CONTINUE
           END IF
  100   CONTINUE
      ELSE
C
C     Read only exponents
C
        NRCIJ = NUCIJ
        DO 300 L = 1, NUCIJ
           CALL RDLINE(LUINFO)
           IF (FREE) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
              ISTART = 1
              CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALPHA(L),'REA')
#else
              READ (MLINE(NMLINE),*) ALPHA(L)
#endif
           ELSE
              READ (MLINE(NMLINE),FMT3) ALPHA(L)
           END IF
        CPRIMU(L,L) = D1
  300   CONTINUE
      ENDIF
C
C     Test print
C
      IF (INPTST) THEN
         WRITE(LUPRI,'(A)') 'Contraction matrix:'
         DO 400 L = 1, NUCIJ
            WRITE (LUPRI,'(E16.10)') ALPHA(L)
            WRITE (LUPRI,'(6(2X,F10.8))') (CPRIMU(L,M),M=1,NRCIJ)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck segorb */
      SUBROUTINE SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU,KPRIM,DSM)
#include "implicit.h"
#include "priunit.h"
      LOGICAL SEGIJ
      DIMENSION CPRIMU(KPRIM,KPRIM)
C
      DO 100 L = 1, NUCIJ
         NONZER = 0
         DO 200 M = 1, NRCIJ
            IF (ABS(CPRIMU(L,M)).GT.DSM) NONZER = NONZER + 1
  200    CONTINUE
         SEGIJ = SEGIJ .AND. (NONZER .LE. 1)
  100 CONTINUE
      RETURN
      END
C  /* Deck nrmorb */
      SUBROUTINE NRMORB(I,NRCIJ,NUCIJ,ALPHA,CPRIM,CPRIMU,KPRIM,NBLOCK)
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &           DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0,
     &           THRMIN = 1.D-17)
      CHARACTER*1 SPDCAR
      DIMENSION ALPHA(KPRIM), CPRIM(KPRIM,KPRIM), CPRIMU(KPRIM,KPRIM)
C
      PIPPI = (DP5/PI)**DP75
      DO 100 N = 1, NRCIJ
         SUM = D0
         DO 200 L = 1, NUCIJ
         DO 200 M = 1, NUCIJ
            T = D2*SQRT(ALPHA(L)*ALPHA(M))/(ALPHA(L)+ALPHA(M))
            SUM = SUM + CPRIMU(L,N)*CPRIMU(M,N)*(T**(I + DP5))
  200    CONTINUE
         IF (SQRT(SUM) .LT. THRMIN) GOTO 1000
         SUM=D1/SQRT(SUM)
         DO 300 L=1, NUCIJ
            CPRIM(L,N)=CPRIMU(L,N)*SUM*(D4*ALPHA(L))**(DP5*I+DP25)*PIPPI
  300    CONTINUE
  100 CONTINUE
      RETURN
 1000 CONTINUE
        WRITE (LUPRI,'(/1X,A,2(I3,A),A1,A/A,2I4)')
     &    'INPUT ERROR: CGTO no.',N,' for block',NBLOCK,
     &    ' of ',SPDCAR(I-1),' type has zero norm.',
     &    ' Contraction matrix: ',NUCIJ,NRCIJ
         CALL QUIT('CGTO with zero norm.')
      END
C  /* Deck orbpro */
      SUBROUTINE ORBPRO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM,
     &           ISHELL,IPRIM,IPRIMD,IORB,IORBD,LCOMP)
C******************************************************************************
C
C     Process basis data
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0)
C
#include "cbirea.h"
      LOGICAL SEG,SPHER
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
#include "ibtfun.h"
      ICENT  = 0
      IBLOCK = 0
      DO 10 I = 1,NONTYP
         DO 20 N = 1,NONT(I)
            ICENT = ICENT + 1
            KBCH = IBLOCK
            NDEG = NUCDEG(ICENT)
            DO 30 J = 1,IQM(I)
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
               DO 40 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NCONT = NRC(KBCH)
                  IF(NCONT.GT.MXCONT) GOTO 5000
C                 ...  dimension PRICCF(MXPRIM,MXCONT)
                  DO 50 KK = 1, NCONT
                     KKK = KKK + 1
                     ISHELL = ISHELL + 1
C*****************************************************************************
C     Data on the AO-block associated with a given shell(ISHELL):
C       NUCO    - number of uncontracted functions
C       NRCO    - number of contracted functions
C     Data on a given shell in an AO-block:
C       NCENT   - index of symmetry independent center
C       NUMCF   - index of shell in AO-block
C       NBCH    - index of block in AO-vector
C       ISTBAO  - stabiliser: basic sym. op. that do not move center
C       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
C       KHKT    - number of spherical (Cartesian) components
C       KCKT    - number of Cartesian components
C       SEGM    - segmented contraction
C       LCLASS  - component: large(1) or small(0)
C       CENT    - coordinates of center
C*****************************************************************************
                     NUCO  (ISHELL) = NUC(KBCH)
                     NRCO  (ISHELL) = NCONT
                     JSTRT (ISHELL) = IPRIM
                     NSTRT (ISHELL) = IORB
                     KSTRT (ISHELL) = IORBD
                     NCENT (ISHELL) = ICENT
                     NUMCF (ISHELL) = KK
                     NUMCFT(ISHELL) = KKK
                     NBCH  (ISHELL) = KBCH
                     SHARE (ISHELL) = .FALSE.
                     ISTBAO(ISHELL) = ISTBNU(ICENT)
                     NHKT  (ISHELL) = J
                     KHKT  (ISHELL) = NCOMP
                     SPHR  (ISHELL) = SPHER
                     KCKT  (ISHELL) = NCCMP
                     SEGM  (ISHELL) = SEG(KBCH)
                     LCLASS(ISHELL) = LCOMP
                     CENT(ISHELL,1,1) = CORD(1,ICENT)
                     CENT(ISHELL,2,1) = CORD(2,ICENT)
                     CENT(ISHELL,3,1) = CORD(3,ICENT)
                     IORB  = IORB  + NCOMP
                     IORBD = IORBD + NCOMP*NDEG
                     IF (.NOT. SEGM(ISHELL)) SEGBAS = .FALSE.
   50             CONTINUE
                  DO 60 L = 1, NUC(KBCH)
                     IPRIM = IPRIM + 1
C*****************************************************************************
C       PRIEXP  - exponent of primitive shell
C       PRICCF  - normalized contraction coefficient
C       PRICRX  - x-coordinate of center
C       PRICRY  - y-coordiante of center
C       PRICRZ  - z-coordinate of center
C*****************************************************************************
                     PRIEXP(IPRIM) = ALPHA(L,KBCH)
                     DO 70  M = 1, NCONT
                        PRICCF(IPRIM,M) = CPRIM(L,M,KBCH)
   70                CONTINUE
                     PRICRX(IPRIM) = CORD(1,ICENT)
                     PRICRY(IPRIM) = CORD(2,ICENT)
                     PRICRZ(IPRIM) = CORD(3,ICENT)
                     IPRIMD = IPRIMD + NCOMP*NDEG
   60             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         IBLOCK = IBLOCK + NBLCK(I)
   10 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> ORBPRO error, '//
     *   'no. of contracted functions per AO-block',NCONT,
     *   '                   '//
     *   'current maximum number                  ',MXCONT
      CALL QUIT('Too many contracted functions per AO-block')
      END
C  /* Deck fitpro */
      SUBROUTINE FITPRO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM,
     &           ISHELL_FIT,IPRIM_FIT,IPRIMD_FIT,
     &           IORB_FIT,IORBD_FIT,LCOMP_FIT)
C******************************************************************************
C
C     Process basis data
C     written by Luuk Visscher, trivial modification of orbpro
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0)
C
#include "cbirea.h"
      LOGICAL SEG,SPHER
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
#include "ccom.h"
#include "nuclei.h"
#include "fitprimit.h"
#include "fitshells.h"
#include "symmet.h"
#include "aosotr.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
#include "ibtfun.h"
      ICENT  = 0
      IBLOCK = 0
      DO 10 I = 1,NONTYP
         DO 20 N = 1,NONT(I)
            ICENT = ICENT + 1
            KBCH = IBLOCK
            NDEG = NUCDEG(ICENT)
            DO 30 J = 1,IQM(I)
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
               DO 40 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NCONT = NRC(KBCH)
                  IF(NCONT.GT.MXCONT) GOTO 5000
C                 ...  dimension PRICCF(MXPRIM,MXCONT)
                  DO 50 KK = 1, NCONT
                     KKK = KKK + 1
                     ISHELL_FIT = ISHELL_FIT + 1
C*****************************************************************************
C     Data on the AO-block associated with a given shell(ISHELL):
C       NUCO    - number of uncontracted functions
C       NRCO    - number of contracted functions
C     Data on a given shell in an AO-block:
C       NCENT   - index of symmetry independent center
C       NUMCF   - index of shell in AO-block
C       NBCH    - index of block in AO-vector
C       ISTBAO  - stabiliser: basic sym. op. that do not move center
C       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
C       KHKT    - number of spherical (Cartesian) components
C       KCKT    - number of Cartesian components
C       SEGM    - segmented contraction
C       LCLASS  - component: large(1) or small(0)
C       CENT    - coordinates of center
C*****************************************************************************
                     NUCO_FIT  (ISHELL_FIT) = NUC(KBCH)
                     NRCO_FIT  (ISHELL_FIT) = NCONT
                     JSTRT_FIT (ISHELL_FIT) = IPRIM_FIT
                     NSTRT_FIT (ISHELL_FIT) = IORB_FIT
                     KSTRT_FIT (ISHELL_FIT) = IORBD_FIT
                     NCENT_FIT (ISHELL_FIT) = ICENT
                     NUMCF_FIT (ISHELL_FIT) = KK
                     NUMCFT_FIT(ISHELL_FIT) = KKK
                     NBCH_FIT  (ISHELL_FIT) = KBCH
                     SHARE_FIT (ISHELL_FIT) = .FALSE.
                     ISTBAO_FIT(ISHELL_FIT) = ISTBNU(ICENT)
                     NHKT_FIT  (ISHELL_FIT) = J
                     KHKT_FIT  (ISHELL_FIT) = NCOMP
                     SPHR_FIT  (ISHELL_FIT) = SPHER
                     KCKT_FIT  (ISHELL_FIT) = NCCMP
                     SEGM_FIT  (ISHELL_FIT) = SEG(KBCH)
                     LCLASS_FIT(ISHELL_FIT) = LCOMP_FIT
                     CENT_FIT(ISHELL_FIT,1,1) = CORD(1,ICENT)
                     CENT_FIT(ISHELL_FIT,2,1) = CORD(2,ICENT)
                     CENT_FIT(ISHELL_FIT,3,1) = CORD(3,ICENT)
                     IORB_FIT  = IORB_FIT  + NCOMP
                     IORBD_FIT = IORBD_FIT + NCOMP*NDEG
   50             CONTINUE
                  DO 60 L = 1, NUC(KBCH)
                     IPRIM_FIT = IPRIM_FIT + 1
C*****************************************************************************
C       PRIEXP  - exponent of primitive shell
C       PRICCF  - normalized contraction coefficient
C       PRICRX  - x-coordinate of center
C       PRICRY  - y-coordiante of center
C       PRICRZ  - z-coordinate of center
C*****************************************************************************
                     PRIEXP_FIT(IPRIM_FIT) = ALPHA(L,KBCH)
                     DO 70  M = 1, NCONT
                        PRICCF_FIT(IPRIM_FIT,M) = CPRIM(L,M,KBCH)
   70                CONTINUE
                     PRICRX_FIT(IPRIM_FIT) = CORD(1,ICENT)
                     PRICRY_FIT(IPRIM_FIT) = CORD(2,ICENT)
                     PRICRZ_FIT(IPRIM_FIT) = CORD(3,ICENT)
                     IPRIMD_FIT = IPRIMD_FIT + NCOMP*NDEG
   60             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         IBLOCK = IBLOCK + NBLCK(I)
   10 CONTINUE
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> FITPRO error, '//
     *   'no. of contracted functions per AO-block',NCONT,
     *   '                   '//
     *   'current maximum number                  ',MXCONT
      CALL QUIT('Too many contracted functions per AO-block')
      END
C  /* Deck orbout */
      SUBROUTINE ORBOUT(ICLS,NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           IPCON,KATOM,KANG,KBLOCK,KPRIM,CPRIMU,NRMPRI)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0)
C
#include "ccom.h"
#include "cbirea.h"
#include "nuclei.h"
#include "primit.h"
      CHARACTER CHRSEG*10,DEG*2
      LOGICAL SEG, NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK),
     &          SEG(KBLOCK)
      DIMENSION IPCON(KPRIM)
C
C     P R I M I T I V E S
C     ===================
C
      IF (IPREAD .GT. 0) THEN
       IF (.NOT.NRMPRI) THEN
        CALL HEADER('Orbital exponents and contraction coefficients',1)
       ELSE
        CALL HEADER
     &  ('Orbital exponents and normalized contraction coefficients',1)
       END IF
      END IF
      IF(ICLS.EQ.1) THEN
        IPRIMD = 0
        IPRIM =  0
      ELSE
        IPRIMD = NPLRG
        IPRIM  = NPLSH
      ENDIF
      ICENT  = 0
      IBLOCK = 0
      DO 100 I = 1, NONTYP
         DO 110 N = 1, NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = IBLOCK
            DO 200 J = 1, IQM(I)
            DO 200 K = 1, JCO(J,I)
               KBCH = KBCH + 1
               NNUC  = NUC(KBCH)
               NNRC  = NRC(KBCH)
            IF (NNUC .EQ. 0) GO TO 200
               ITYP = NHKOFF(J)
               IPSTRT = IPRIM + 1
               IPRIM =  IPRIM + NNUC
               IF (SEG(KBCH)) THEN
                  CHRSEG = 'seg. cont.'
               ELSE
                  CHRSEG = 'gen. cont.'
               END IF
               IF (NNRC .LE. 5) THEN
                  ASSIGN 1000 TO IFMTA
                  ASSIGN 1010 TO IFMTB
                  ASSIGN 1020 TO IFMTC
                  ASSIGN 1030 TO IFMTD
               ELSE
                  ASSIGN 1005 TO IFMTA
                  ASSIGN 1015 TO IFMTB
                  ASSIGN 1025 TO IFMTC
                  ASSIGN 1035 TO IFMTD
               END IF
               DEG = '  '
               DO 400 L = 1, NDEG
                 IF (NDEG .GT. 1) WRITE(DEG,'(A1,I1)') '#',L
                 ITYP = NHKOFF(J)
                 DO 300 ICOMP = 1, KHK(J)
                     ITYP = ITYP + 1
                     IPRIMD = IPRIMD + 1
                     IF (IPREAD .GT. 0) THEN
                       WRITE (LUPRI,IFMTA) NAMN(ICENT),DEG,
     &                                     GTOTYP(ITYP),
     &                                     IPRIMD,PRIEXP(IPSTRT),
     &                                     (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                     ENDIF
                     IF (NNUC .GT. 1) THEN
                        IPRIMD = IPRIMD + 1
                        IF (IPREAD .GT. 0) THEN
                           WRITE (LUPRI,IFMTC) CHRSEG,
     &                                     IPRIMD,PRIEXP(IPSTRT+1),
     &                                     (CPRIMU(2,MM,KBCH),MM=1,NNRC)
                        END IF
                     END IF
                     DO 410 M = 3, NNUC
                        IPRIMD = IPRIMD + 1
                        IF (IPREAD .GT. 0) THEN
                          WRITE (LUPRI,IFMTD) IPRIMD,PRIEXP(IPSTRT-1+M),
     &                                     (CPRIMU(M,MM,KBCH),MM=1,NNRC)
                       END IF
  410              CONTINUE
  300            CONTINUE
  400          CONTINUE
  200       CONTINUE
  110    CONTINUE
         IBLOCK = IBLOCK + NBLCK(I)
  100 CONTINUE
 1000 FORMAT(/2X,A4,A2,1X,A4,I4,1P,D18.10,0P,2X,5F8.4)
 1010 FORMAT(/2X,A4,3X,   A4,I4,1P,D18.10,0P,2X,5F8.4)
 1020 FORMAT( 3X,A10,        I4,1P,D18.10,0P,2X,5F8.4)
 1030 FORMAT(13X,            I4,1P,D18.10,0P,2X,5F8.4)
 1005 FORMAT(/2X,A4,A2,1X,A4,I4,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))
 1015 FORMAT(/2X,A4,3X,   A4,I4,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))
 1025 FORMAT( 3X,A10,        I4,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))
 1035 FORMAT(13X,            I4,1P,D18.10,0P,2X,5F8.4/,(37X,5F8.4))
C
      IF (NRMPRI) RETURN
C
C       C O N T R A C T I O N
C       =====================
C
C       Looping is over
C                       - atomic type
C                         - symmetry independent center
C                           - shell
C                             - component
C                               - degeneracy of center
C
      CALL HEADER('Contracted Orbitals',1)
      IF(ICLS.EQ.1) THEN
        IPRIMD  = 0
        IORBD   = 0
        IPRIM   = 0
      ELSE
        IPRIMD  = NPLRG
        IORBD   = NLARGE
        IPRIM   = NPLSH
      ENDIF
      ICENT   = 0
      IBLOCK  = 0
      DO 500 I = 1,NONTYP
         DO 510 N = 1,NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = IBLOCK
            DO 600 J = 1,IQM(I)
               NCOMP = KHK(J)
               DO 610 K = 1,JCO(J,I)
                  KBCH = KBCH + 1
                  NNUC = NUC(KBCH)
                  NNRC  = NRC(KBCH)
                  IPSTRT = IPRIM + 1
                  IPRIM  = IPRIM + NNUC
                  DO 700 L = 1,NNRC
                     JPRIM = 0
                     JPRIMD = IPRIMD
                     DO 750 M = IPSTRT,IPRIM
                        JPRIMD = JPRIMD + 1
                        IF(PRICCF(M,L).NE.0.00D0) THEN
                           JPRIM = JPRIM + 1
                           IPCON(JPRIM) = JPRIMD
                        ENDIF
  750                CONTINUE
                     IOFF = 0
                     DEG = '  '
                     DO 810 LL = 1,NDEG
                       IF(NDEG.GT.1) WRITE(DEG,'(A1,I1)') '#',LL
                       ITYP = NHKOFF(J)
                       DO 800 ICOMP = 1, NCOMP
                         ITYP = ITYP + 1
                         IORBD = IORBD + 1
                         IF (IPREAD .GT. 0) THEN
                            WRITE(LUPRI,1040)
     &                      IORBD,NAMN(ICENT),DEG,GTOTYP(ITYP),
     &                      ((IPCON(M)+IOFF),M=1,JPRIM)
                         END IF
                        IOFF = IOFF + NNUC
  800                  CONTINUE
  810                CONTINUE
  700             CONTINUE
                  IPRIMD = IPRIMD + NNUC*NDEG*NCOMP
  610          CONTINUE
  600       CONTINUE
  510    CONTINUE
         IBLOCK = IBLOCK + NBLCK(I)
  500 CONTINUE
      IF (IPREAD .GT. 0) THEN
         WRITE (LUPRI,'(/)')
      END IF
 1040 FORMAT(I4,2X,A4,A2,2X,A4,24I4)
 1050 FORMAT(I4,2X,A4,4X,A4,12I6)
C
      RETURN
      END
C  /* Deck nucpro */
      SUBROUTINE NUCPRO(WORK,LWORK)
C******************************************************************************
C
C     Process molecular data
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
#include "cbirea.h"
      DIMENSION WORK(LWORK)
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
C
C
#include "ibtfun.h"
C
C     Statement function to "convert" 0 - 7 to a printable
C     "binary" form
C
      IFAKBN(I) = 25*IBTAND(I,4) + 5*IBTAND(I,2) + IBTAND(I,1)
C
C     *************************************************************
C     ***** Calculate Distances ,Nuclear Potential Energy and *****
C     ***** Nuclear Contributions to Dipole Moment            *****
C     * LOOPS:
C     * 10      - loop over first symmetry independent centre
C     * 20+IF   - loop over dependent centers
C     * 30      - loop over second symmetry independent centre
C     * 40+IF   - loop over dependent centers
C     ************************************************************
C
      DPNUC(1) = D0
      DPNUC(2) = D0
      DPNUC(3) = D0
      NUCDEP    = 0
      NUCPRE(1) = 0
      NFLOAT    = 0
      DO 100 N = 1, NUCIND
         IVARB = ISTBNU(N)
         NDEG  = MULT(IVARB)
         IF (IPREAD .GE. 5)
     &      WRITE (LUPRI,1000) NAMN(N),(CORD(I,N),I=1,3),
     &                        IFAKBN(ISTBNU(N)), NDEG - 1
         NUCDEG(N) = NDEG
         NUCPRE(N + 1) = NUCPRE(N) + NDEG
         IF (NINT(CHARGE(N)) .EQ. 0 .AND. NAMN(N) .NE. 'cav')
     &        NFLOAT = NFLOAT + 1
         DPNUC(1) = DPNUC(1) + CORD(1,N)*CHARGE(N)*FMULT(IVARB)
         DPNUC(2) = DPNUC(2) + CORD(2,N)*CHARGE(N)*FMULT(IVARB)
         DPNUC(3) = DPNUC(3) + CORD(3,N)*CHARGE(N)*FMULT(IVARB)
         II = 0
         DO 200 LA = 0,MAXREP
            IF (IBTAND(IVARB,LA) .EQ. 0) THEN
               II = II + 1
               NUCDEP = NUCDEP + 1
               NUCNUM(N,LA+1) = NUCDEP
               IF (NDEG .EQ. 1 ) THEN
                  NAMDEP(NUCDEP)     = NAMEX(3*N)(1:4)//'  '
                  NAMDPX(3*NUCDEP-2) = NAMEX(3*N-2)
                  NAMDPX(3*NUCDEP-1) = NAMEX(3*N-1)
                  NAMDPX(3*NUCDEP  ) = NAMEX(3*N  )
                  NDEGNM(NUCDEP)     = 1
               ELSE
                  ICHARD = ICHAR('0') + II
                  NAMDEP(NUCDEP)     = NAMEX(3*N)(1:4)//' '
     &                //CHAR(ICHARD)
                  NAMDPX(3*NUCDEP-2) = NAMEX(3*N-2)(1:3)
     &                //CHAR(ICHARD)//NAMEX(3*N-2)(5:6)
                  NAMDPX(3*NUCDEP-1) = NAMEX(3*N-1)(1:3)
     &                //CHAR(ICHARD)//NAMEX(3*N-1)(5:6)
                  NAMDPX(3*NUCDEP  ) = NAMEX(3*N  )(1:3)
     &                //CHAR(ICHARD)//NAMEX(3*N  )(5:6)
                  NDEGNM(NUCDEP)     = II
               END IF
            ELSE
               NUCNUM(N,LA+1) = 0
            END IF
  200    CONTINUE
  100 CONTINUE
      NATOMS = NUCDEP - NFLOAT
      IF (NUCDEP .GT. MXCENT) GOTO 5010
C
C     Nuclear repulsion energy
C     ========================
C
      POTNUC =  GETPOT(CHARGE)
C
C     **********************************************************
C     ***** Calculate symmetry-adapted nuclear coordinates *****
C     **********************************************************
C
      CALL SYMNCO(LWORK,WORK)
C
      RETURN
C
C     Error messages:
C
 5010 CONTINUE
      WRITE (LUPRI,'(///1X,A,/,2(1X,A,I3/),/,1X,A,I3,A)')
     &      ' Allowed number of atoms exceeded.',
     &      ' Number of atoms in input:',NUCDEP,
     &      ' Number of atoms allowed: ',MXCENT,
     &      ' Increase MXCENT to ',NUCDEP,'.'
      CALL QUIT('Too many atoms in READIN.')
C
 1000 FORMAT(2X,A4,5X,3F10.5,6X,/,
     &       16X, 'Stabilizer ',I3,',   with ',I2,
     &       ' symmetry equivalent atoms',/)
      END
C  /* Deck syminp */
      SUBROUTINE SYMINP(NSYMOP,KASYM,IFXYZ,CLASS)
C*****************************************************************************
C
C     This subroutine sets up point group symmetry and looks at
C     the behaviour of principal axes and rotations under point group
C     symmetry
C
C     tsaue - 940825 - major revision: tagged on SYMGRP + polish
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "ccom.h"
#include "cbirea.h"
#include "chrxyz.h"
      CHARACTER*1 KASYM(3,3)
      CHARACTER*(*) CLASS
      DIMENSION IFXYZ(3),IS(0:7),IGEN(3)
      DATA IS/0,1,1,2,1,2,2,3/
#include "sdpre.h"
#include "ibtfun.h"
C
C     Initialization
C     ==============
C
C     PT is parity of a bitstring:
C       1 for an even number of ones: 000,011,110,101
C      -1 for an odd  number of ones: 001,010,100,111
C
      PT(0) =  1.0D0
      PT(1) = -1.0D0
      PT(2) = -1.0D0
      PT(3) =  1.0D0
      PT(4) = -1.0D0
      PT(5) =  1.0D0
      PT(6) =  1.0D0
      PT(7) = -1.0D0
      DO 5 I = 1,3
         IFXYZ(I)    = 0
         ISYMAX(I,1) = 0
         IGEN(I)     = 0
    5 CONTINUE
C
C     Determine:
C     IGEN(I)   - basic operations
C     ISYMAX(I) - behavior of principal axes under basic operations
C     ===============================================================
C
      IAXIS = 0
      MAXREP = 2**NSYMOP - 1
      IF (NSYMOP.GT.0.AND.IPREAD .GT. 0) THEN
         CALL HEADER('Symmetry Operations',1)
         WRITE (LUPRI,'(1X,A,I2,/)') ' Symmetry operations:',NSYMOP
      END IF
      DO 100 J = 1,NSYMOP
        DO 110 I = 1,3
        IF (KASYM(I,J).NE.' ') THEN
          K = ICHAR(KASYM(I,J)) - ICHAR('W')
          IGEN(J)     = IGEN(J)     + 2**(K-1)
          ISYMAX(K,1) = ISYMAX(K,1) + 2**(J-1)
        END IF
  110   CONTINUE
        IAXIS = IBTOR(IAXIS,IGEN(J))
  100 CONTINUE
C
C     Determine IFXYZ
C     ===============
C     Do we really need it ????
C
      IND = 0
      DO 140 I = 1,NSYMOP
        IND = IBTOR(IND,IGEN(I))
  140 CONTINUE
      DO 141 I = 1,3
        IFXYZ(I) = IBTAND(IBTSHR(IND,(I-1)),1)
  141 CONTINUE
C
C     Determine:
C     ISYMAX(I,2) - behaviour of principal rotations under basic operations
C     =====================================================================
C
      ISYMAX(1,2) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
      ISYMAX(2,2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
      ISYMAX(3,2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
C
C     Determine:
C     IPTAX   - coordinate axis: pointer analogous to IPTSYM
C     NAXREP  - number of coordinate axis in each symmetry
C     ========================================================
C
      CALL IZERO(IPTAX,6)
      CALL IZERO(IPTXYZ(1,0,1),48)
      DO 200 ITYPE = 1,2
        IPTAXI = 0
        DO 205 IREP = 0, MAXREP
           NAXIS = 0
           DO 210 ICOOR = 1, 3
              IF (IBTXOR(IREP,ISYMAX(ICOOR,ITYPE)) .EQ. 0) THEN
                 NAXIS  = NAXIS + 1
                 IPTAXI = IPTAXI + 1
                 IPTAX(ICOOR,ITYPE) = IPTAXI
                 IPTXYZ(NAXIS,IREP,ITYPE) = ICOOR
              END IF
  210      CONTINUE
           NAXREP(IREP,ITYPE) = NAXIS
  205   CONTINUE
  200 CONTINUE
      IF (IPREAD.GT. 20) THEN
         WRITE (LUPRI,'(1X,A,3I5)') ' IPTAX(*,1)  ',
     &                               (IPTAX(I,1),I=1,3)
         WRITE (LUPRI,'(1X,A,3I5)') ' IPTAX(*,2)  ',
     &                               (IPTAX(I,2),I=1,3)
         WRITE (LUPRI,'(1X,A,8I5)') ' NAXREP(*,1) ',
     &                               (NAXREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(1X,A,8I5)') ' NAXREP(*,2) ',
     &                               (NAXREP(I,2),I=0,MAXREP)
      END IF
C
C     Determine group & properties
C
      CALL SYMGRP(IGEN,NSYMOP,CLASS)
C
C     For DIRAC: double group symmetry
C
      IF(DIRAC) CALL DBLGRP
C
C     Determine MULT(I) and FMULT(I) - multiplicity of center
C     =======================================================
C
      DO 101 I = 0,7
         MULT(I)  = 2**MAX(0,NSYMOP-IS(I))
         FMULT(I) = SDPRE(MULT(I))
  101 CONTINUE
C
      RETURN
C
      END
C  /* Deck sympro */
      SUBROUTINE SYMPRO(LVALUE,MVALUE,NVALUE,IRREP,DOOWN)
C*****************************************************************************
C
C     SYMPRO generates the symmetry orbitals (SOs) from the list of
C     symmetry independent atomic orbitals(AOs).
C     A given SO of irrep L is given as:
C
C         SO(L,a) = SUM{i} [CTRAN(NSORB,i)*G(i)*AO(NSORB)]
C
C     Here {G(i)} is a set of symmetry operators transforming AO(NSORB)
C     into all corresponding orbitals centered on symmetry related sites
C     in the molecule. The loop over i is only allowed to encounter
C     elements that do not contain any basic operations belonging to the
C     stabiliser of the centre; thus any G uniquely defines transformation
C     to a given centre.
C
C     By allowing all symmetry operations to work on the list of
C     symmetry independent AOs a list of symmetry dependent AOs
C     is generated (counted by NAORB). Redundancies are removed by
C     limiting symmetry operations to the set G(i) and the list thus
C     reduced to the list of non-trivial symmetry dependent AOs
C     (counted by NSORB).
C
C******************************************************************************
C
C     This subroutine determines the following parameters:
C
C     ISYMAO  -  indicates how an atomic orbital behaves under the basic
C                symmetry operations. When a basic operation changes the
C                sign of the AO (centered at origo) the corresponding bit
C                gets the value 1.
C
C     CTRAN   -  transformation coefficients (+1 or -1)
C     IPIND   -  packed integer
C     IPTSYM  -  pointer from redundant list of symmetry dependent AOs
C                to list of non-trivial symmetry dependent AOs. The
C                pointer is zero if the AO does not contribute to the
C                irrep.
C     ITRAN   -  given an operation G(i): pointer from list of symmetry
C                dependent AOs to list of symmetry independent AOs
C     JPRX    -  AO label
C     JTRAN   -  index of operation G moving a given AO to another centre
C     MAMN    -  name of AO centre
C     NPARNU  -  offset pointer from non-symmetric operators for given
C                irrep
C     NPARSU  -  offset pointer for symmetry dependent AOs for given irrep
C     NSORB    -  number of symmetry dependent AOs
C     NSYM    -  number of irreps that the AOs contribute to
C*************************************************************************
C
C     The code originally resided inside READIN
C     tuh 120988
C
C     Polished and annotated - tsaue March 10 1993
C*************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0)
C
      LOGICAL DOOWN
      DIMENSION LVALUE(MXAQN), MVALUE(MXAQN), NVALUE(MXAQN),
     &          IRREP(MXCORB)
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pincom.h"
#include "aosotr.h"
#include "sphtrm.h"
#include "huckel.h"
#include "cbirea.h"
C
#include "ibtfun.h"
C*****************************************************************************
C Determine ISYMAO:
C       ISYMAO is a bitstring indicating how atomic orbitals centered in
C       origo behaves under the basic symmetry operations
C*****************************************************************************
C
C     I. Cartesian functions
C     ----------------------
C
      IF (DOCART) THEN
         DO 100 I = 1, NHTYP
            CALL LMNVAL(I,I*(I + 1)/2,LVALUE,MVALUE,NVALUE)
            DO 100 J = 1, I*(I + 1)/2
                LVAR = MOD(LVALUE(J),2)*ISYMAX(1,1)
                MVAR = MOD(MVALUE(J),2)*ISYMAX(2,1)
                NVAR = MOD(NVALUE(J),2)*ISYMAX(3,1)
                ISYMAO(I,J) = IBTXOR(LVAR,IBTXOR(MVAR,NVAR))
  100    CONTINUE
C
C     II. Your own scheme
C     -------------------
C
      ELSE IF (DOOWN) THEN
         IJK = 0
         DO 200 I = 1, NHTYP
            DO 210 K = 1, KHK(I)
               ISYMAO(I,K) = -1
  210       CONTINUE
            CALL LMNVAL(I,I*(I + 1)/2,LVALUE,MVALUE,NVALUE)
            DO 220 J = 1, I*(I + 1)/2
                LVAR = MOD(LVALUE(J),2)*ISYMAX(1,1)
                MVAR = MOD(MVALUE(J),2)*ISYMAX(2,1)
                NVAR = MOD(NVALUE(J),2)*ISYMAX(3,1)
                DO 220 K = 1, KHK(I)
                   IJK = IJK + 1
                   IF (ABS(CSP(IJK)).GT.D0) THEN
                      IF (ISYMAO(I,K) .EQ. -1) THEN
                         ISYMAO(I,K) = IBTXOR(LVAR,IBTXOR(MVAR,NVAR))
                      ELSE
                         WRITE (LUPRI,'(1X,A,/,A,/A,2I5)')
     &                   ' Incorrect Cartesian transformation input:',
     &                   ' Components belonging to different irreps '
     &                   //'have been combined.',
     &                   ' Angular quantum number and component:',I-1,K
                         CALL QUIT('Error in input')
                      END IF
                   END IF
  220       CONTINUE
  200    CONTINUE
      ELSE
C
C     III. Spherical harmonics
C     ------------------------
C
         DO 300 I = 0, NHTYP - 1
            II = I + 1
            IF (I .EQ. 0) THEN
               ISYMAO(II,1) = IREPLM(0,0)
            ELSE IF (I .EQ. 1) THEN
               ISYMAO(II,1) = IREPLM(1, 1)
               ISYMAO(II,2) = IREPLM(1,-1)
               ISYMAO(II,3) = IREPLM(1, 0)
            ELSE
               DO 310 J = -I,I
                  JJ = J + I + 1
                  ISYMAO(II,JJ) = IREPLM(I,J)
  310          CONTINUE
            END IF
  300    CONTINUE
      END IF
C*****************************************************************************
C Determine transformation from AOs to SOs :
C   * DO-loops:
C       400     - run over irreducible representations
C       410     - run over shells
C       420     - run over components, thus (60+70):run over orbitals
C       430     - run over the symmetryoperations that have no part in the
C               stabilizer of the center, that is run over a unique set
C               of symmetry operations transforming between symmetry
C               dependent centers
C   * Counting variables:
C       NSORB   -
C       NAORBD  -
C       NAORB   -
C       NAA     -
C       NSORBD  -
C       JKB     -
C       IRREPN  -
C*****************************************************************************
      IHNMAO    = 0
      NSORB     = 0
      NPARSU(1) = 0
      NSYM      = 0
      NPARLA    = 0
      CALL IZERO(NAOS,8)
      CALL IZERO(NCOS,16)
      CALL IZERO(IRREP,MXCORB)
C
C     Run over Irreducible Representations LAMBDA
C
      DO 400 LA = 0, MAXREP
         NAORBD = 0
         NAORB  = 0
         IBLOCK = 0
         DO 410 IA = 1, KMAX
            MULA   = ISTBAO(IA)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            NUMCFA = NUMCF(IA)
            IC     = LCLASS(IA)
            IF (NUMCFA .EQ. 1) IBLOCK = IBLOCK + 1
            IOFF = NAORBD
            DO 420 NA = 1, KHKTA
               NAORB = NAORB+1
               IVARB = IBTXOR(LA,ISYMAO(NHKTA,NA))
C
C              If orbital contributes to this representation:
C
c                (IVARB is a bitstring of basic operations with
C                1 in the positions where the basic operations
C                has a different parity for the irrep and the
C                AO(when centered in origo). These basic
C                operations with different parities must then
C                not be part of the stabilizer of the center
C                for the AO).
C
               IF (IBTAND(MULA,IVARB) .EQ. 0) THEN
                  NSORB = NSORB + 1
                  JKB = 0
C
C                 loop over symmetry dependent centers
C
                  DO 430 KB = 0, MAXREP
                  IF (IBTAND(KB,MULA) .EQ. 0) THEN
                     JKB              = JKB + 1
                     NAORBD           = NAORBD + 1
                     CTRAN(NSORB,JKB) = PT(IBTAND(KB,IVARB))
                     ITRAN(NSORB,JKB) = IOFF + (JKB-1)*KHKTA + NA
                  END IF
  430             CONTINUE
                  NAOS(LA+1)    = NAOS(LA+1) + 1
                  NCOS(LA+1,IC) = NCOS(LA+1,IC) + 1
                  IPTYP(NSORB)  = NHKOFF(NHKTA) + NA
                  IPCEN(NSORB)  = NCENT(IA)
                  ICLASS(NSORB) = IC
#ifndef PRG_DIRAC
                  IF (IC .EQ. 2 .AND. ADDSTO) THEN
                     IHNMAO = IHNMAO + 1
                     IHUCPT(IA) = IHNMAO
                     HUCCHR = CHARGE(NCENT(IA))
                     NHUCCO = NUMCF(IA)
                     CALL HUCFUN(HUCCHR,NHKTA,NHUCCO,HUCEXC(IHNMAO))
                  END IF
#endif
                  IRREPN        = IRREP(NAORB)
                  IRREP(NAORB)  = IRREPN + 1
C                  IPIND(NSORB)  = IA*2**8 + NA*2**4 + IRREPN
                  IPIND(NSORB)  = IA*2**16 + NA*2**8 + IRREPN
                  JTRAN(NSORB)  = JKB
                  INDGEN(NSORB) = IRREPN*2**24 + NA*2**16 + NUMCFA*2**8
     &                                                    + IBLOCK
                  IPTSYM(NAORB,LA) = NSORB
               ELSE
                  IPTSYM(NAORB,LA) = 0
                  NAORBD = NAORBD + MULT(MULA)
               END IF
  420       CONTINUE
  410    CONTINUE
         IF (LA .GT. 0) NPARSU(LA+1) = NPARSU(LA) + NPARLA
         NBASI  = NAOS(LA+1)
         NPARLA = NBASI*(NBASI + 1)/2
         IF (NBASI .GT. 0) NSYM = NSYM + 1
  400 CONTINUE
C
C     Compute symmetry offsets for non-symmetric operators
C     (Note case IREPO = 0 is thereby excluded)
C
      DO 500 IREPO = 1,MAXREP
         IBLK = 0
         DO 500 IREPA = 0,MAXREP
            IREPB = IBTXOR(IREPO,IREPA)
            IF (IREPA .GT. IREPB) THEN
               NPARNU(IREPO+1,IREPA+1) = IBLK
               IBLK = IBLK + NAOS(IREPA+1)*NAOS(IREPB+1)
            ENDIF
 500  CONTINUE
C
C     Compute symmetry offsets for basis functions
C
      IBAS = 0
      DO ISYM = 1,MAXREP+1
        DO IC = 1,2
          ICOS(ISYM,IC) = IBAS
          IBAS = IBAS + NCOS(ISYM,IC)
        ENDDO
      ENDDO
      DO JC = 1,2
        DO IC = 1,2
          DO JSYM = 1,MAXREP+1
            DO ISYM = 1,MAXREP+1
              I2COSX(ISYM,JSYM,IC,JC) =
     &               NSORB*ICOS(JSYM,JC)+ICOS(ISYM,IC)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C  /* Deck symout */
      SUBROUTINE SYMOUT
C*****************************************************************************
C
C     Print routine for symmetry
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pgroup.h"
#include "pincom.h"
#include "cbirea.h"
#include "aosotr.h"
#include "chrsgn.h"
      CHARACTER COMP(2)*1
C
#include "ibtfun.h"
C*****************************************************************************
C       S Y M M E T R Y      O R B I T A L S
C*****************************************************************************
      CALL HEADER('Symmetry Orbitals',1)
      WRITE (LUPRI,'(A,6X,8I5)')
     &  '  Number of orbitals in each symmetry: ',(NAOS(I),I=1,MAXREP+1)
      IF(DIRAC) THEN
        COMP(1) = 'L'
        COMP(2) = 'S'
        WRITE (LUPRI,'(A,8I5)')
     &   '  Number of large orbitals in each symmetry: ',
     &    (NCOS(I,1),I=1,MAXREP+1)
        WRITE (LUPRI,'(A,8I5)')
     &   '  Number of small orbitals in each symmetry: ',
     &    (NCOS(I,2),I=1,MAXREP+1)
      ELSE
        COMP(1) = ' '
        COMP(2) = ' '
      ENDIF
      I=0
      DO 10 LA=1,MAXREP+1
        NBI=NAOS(LA)
        IF (NBI.GT.0) THEN
          WRITE (LUPRI,'(//A,2X,A3,A1,I2,A1/)')
     &          '  Symmetry',REP(LA-1),'(',LA,')'
          DO 20 L=1,NBI
            I=I+1
            ICENT = IPCEN(I)
            J=NUCDEG(ICENT)
            WRITE (LUPRI,1060) I,NAMN(ICENT),COMP(ICLASS(I)),
     +         GTOTYP(IPTYP(I)), ITRAN(I,1),
     *         (CHRSGN(NINT(CTRAN(I,K))),ITRAN(I,K),K=2,J)
   20     CONTINUE
        ELSE
          WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1)')
     *      'No orbitals in symmetry',REP(LA-1),'(',LA,')'
        END IF
  10  CONTINUE
      IF (IPREAD .GT. 10) THEN
         WRITE(LUPRI,'(10X,A,/)') 'Symmetry pointer indices'
         WRITE(LUPRI,'(4X,8I4)')
     &      ((IPTSYM(I,J),J = 0,MAXREP),I=1,NORBS)
         IF (NFLOAT .GT. 0) WRITE (LUPRI,'(/10X,A,I4)')
     *       'Number of floating orbitals :',NFLOAT
      END IF
      WRITE (LUPRI,'(/1X,A,3(2X,A3,A1,I1,A1))')
     &        ' Symmetries of electric field:',
     &        (REP(ISYMAX(I,1)),'(',ISYMAX(I,1)+1,')',I=1,3)
      WRITE (LUPRI,'(/1X,A,3(2X,A3,A1,I1,A1))')
     &        ' Symmetries of magnetic field:',
     &        (REP(ISYMAX(I,2)),'(',ISYMAX(I,2)+1,')',I=1,3)
 1060 FORMAT(I5,3X,A6,1X,A1,3X,A4,5X,I5,7(2X,A,1X,I5))
      RETURN
      END
C  /* Deck symnco */
      SUBROUTINE SYMNCO(LWORK,WORK)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
C
#include "chrxyz.h"
#include "cbirea.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION WORK(LWORK)
#include "ibtfun.h"
      LENGTH = 48*NUCIND
      CALL IZERO(IPTCNT,LENGTH)
      IOFFT = 0
      IOFFR = 0
      IOFFTN = 0
      DO 100 IREP = 0,MAXREP
         ICENTA = 0
         ICENTB = 0
         DO 110 ICENT = 1,NUCIND
            MULC = ISTBNU(ICENT)
            DO 200 IDIRN = 1,3
               ICENTA = ICENTA + 1
               IF (IBTAND(MULC,IBTXOR(IREP,ISYMAX(IDIRN,1))).EQ.0)THEN
                  IOFFT = IOFFT + 1
                  IPTCNT(ICENTA,IREP,1) = IOFFT
               ENDIF
               IF (IBTAND(MULC,IBTXOR(IREP,ISYMAX(IDIRN,2))).EQ.0)THEN
                  IOFFR = IOFFR + 1
                  IPTCNT(ICENTA,IREP,2) = IOFFR
               ENDIF
  200       CONTINUE
            DO 300 IDIR1 = 1, 3
            DO 300 IDIR2 = IDIR1, 3
               ICENTB = ICENTB + 1
               ISYMIJ = IBTXOR(ISYMAX(IDIR1,1),ISYMAX(IDIR2,1))
               IF (IBTAND(MULC,IBTXOR(IREP,ISYMIJ)) .EQ. 0) THEN
                  IOFFTN = IOFFTN + 1
               ENDIF
  300       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IOFF = 0
      DO 400 IREP = 0, MAXREP
         DO 410 ICENT = 1, NUCIND
            IPTNUC(ICENT,IREP) = 0
  410    CONTINUE
         DO 420 ICENT = 1, NUCIND
            MULC = ISTBNU(ICENT)
            IF (IBTAND(MULC,IREP) .EQ. 0) THEN
               IOFF = IOFF + 1
               IPTNUC(ICENT,IREP) = IOFF
            END IF
  420    CONTINUE
  400 CONTINUE
C
      IF (IPREAD .GE. 20) THEN
         CALL HEADER('Symmetry-adapted nuclear coordinates',-1)
         WRITE(LUPRI,'(11X,8I4)') (IREP, IREP = 0,MAXREP)
         ICENTA = 0
         DO 500 ICENT = 1,NUCIND
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     *      NAMN(ICENT),'X', (IPTCNT(ICENTA+1,IREP,1),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     *      NAMN(ICENT),'Y', (IPTCNT(ICENTA+2,IREP,1),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     *      NAMN(ICENT),'Z', (IPTCNT(ICENTA+3,IREP,1),IREP = 0,MAXREP)
            ICENTA = ICENTA + 3
  500    CONTINUE
         CALL HEADER('Symmetry-adapted nuclear magnetic moments',-1)
         WRITE(LUPRI,'(11X,8I4)') (IREP, IREP = 0,MAXREP)
         ICENTA = 0
         DO 510 ICENT = 1,NUCIND
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     *      NAMN(ICENT),'X', (IPTCNT(ICENTA+1,IREP,2),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     *      NAMN(ICENT),'Y', (IPTCNT(ICENTA+2,IREP,2),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     *      NAMN(ICENT),'Z', (IPTCNT(ICENTA+3,IREP,2),IREP = 0,MAXREP)
            ICENTA = ICENTA + 3
  510    CONTINUE
         CALL HEADER('Symmetry-adapted nuclei',-1)
         WRITE (LUPRI,'(11X,8I4)') (IREP, IREP = 0, MAXREP)
         DO 520 ICENT = 1, NUCIND
            WRITE(LUPRI,'(4X,A1,2X,A1,8I4)')
     &      NAMN(ICENT),CHRXYZ(1),(IPTNUC(ICENT,IREP),IREP = 0,MAXREP)
  520    CONTINUE
      END IF
C
C     *****************************************************************
C     **** Transformation matrix Cartesian to symmetry coordinates ****
C     *****************************************************************
C
      NCOOR  = 3*NUCDEP
      KCSTRA = 1
      KSCTRA = KCSTRA + NCOOR*NCOOR
      KLAST  = KSCTRA + NCOOR*NCOOR
      IF (KLAST.GT.LWORK) CALL STOPIT('NUCPRO','TRACOR',KLAST,LWORK)
      CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,NCOOR,IPREAD)
      CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),2,NCOOR,IPREAD)
      CALL TRACR(NCOOR)
      END
C  /* Deck pricar */
      SUBROUTINE PRICAR(WORK,LWORK)
C
C     This subroutine prints information about atomic coordinates
C
C     tuh 081188 - Bush elected
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "chrxyz.h"
#include "chrsgn.h"
C
      CHARACTER*4 NAME
      INTEGER NBASE(8), NSIGN(8)
      DIMENSION WORK(LWORK)
C
#include "nuclei.h"
#include "frame.h"
#include "symmet.h"
#include "pgroup.h"
#include "cbirea.h"
C
#include "chrnos.h"
#include "ibtfun.h"
C
#include "memint.h"
      NCOOR = 3*NUCDEP
C
C     *********************************
C     ***** Cartesian Coordinates *****
C     *********************************
C
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Cartesian Coordinates',1)
         WRITE (LUPRI,'(A,I3//)')'  Total number of coordinates:',NCOOR
         ICRX = 1
         ICRY = 2
         ICRZ = 3
         DO 100 ICENT = 1, NUCIND
            MULCNT = ISTBNU(ICENT)
            NAME   = NAMEX(3*ICENT)(1:4)
            IF (MULT(MULCNT) .EQ. 1) THEN
               WRITE (LUPRI,'(I4,3X,A,5X,A,3X,F15.10)')
     *               ICRX, NAME,CHRXYZ(-1), CORD(1,ICENT)
               WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')
     *               ICRY, CHRXYZ(-2), CORD(2,ICENT)
               WRITE (LUPRI,'(I4,12X,A,3X,F15.10,/)')
     *               ICRZ, CHRXYZ(-3), CORD(3,ICENT)
               ICRX = ICRX + 3
               ICRY = ICRY + 3
               ICRZ = ICRZ + 3
            ELSE
               JATOM = 0
               DO 200 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                     JATOM = JATOM + 1
                     WRITE (LUPRI,'(I4,3X,A,I2,3X,A,3X,F15.10)')
     *                     ICRX, NAME,JATOM,CHRXYZ(-1),
     *                     PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                     WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')
     *                     ICRY,CHRXYZ(-2),
     *                     PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                     WRITE (LUPRI,'(I4,12X,A,3X,F15.10,/)')
     *                     ICRZ,CHRXYZ(-3),
     *                     PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                     ICRX = ICRX + 3
                     ICRY = ICRY + 3
                     ICRZ = ICRZ + 3
                  END IF
  200          CONTINUE
            END IF
  100    CONTINUE
      END IF
C
C     ********************************
C     ***** Symmetry Coordinates *****
C     ********************************
C
      IF (LWORK .LT. 2*NCOOR*NCOOR) CALL STOPIT('PRICAR','TRACOR',
     &                                          2*NCOOR*NCOOR,LWORK)
      CALL MEMGET('REAL',KCS ,NCOOR*NCOOR,   WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSC ,NCOOR*NCOOR,   WORK,KFREE,LFREE)
      CALL TRACOR(WORK(KCS),WORK(KSC),1,NCOOR,IPREAD)
      IF (MAXREP .GT. 0) THEN
         IF (IPREAD .GT. 0) THEN
            CALL HEADER('Symmetry Coordinates',1)
            WRITE (LUPRI,'(1X,A,8I3)')
     *            ' Number of coordinates in each symmetry: ',
     *            (NCRREP(I,1),I=0,MAXREP)
         END IF
         DO 300 ISYM = 0, MAXREP
         IF (NCRREP(ISYM,1) .GT. 0) THEN
            IF (IPREAD .GT. 0) THEN
               WRITE (LUPRI,'(//1X,A,I2/)') ' Symmetry', ISYM + 1
            END IF
            DO 400 IATOM = 1, NUCIND
               DO 500 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM,1)
                  IF (ISCOOR .GT. 0) THEN
                     NB = 0
                     DO 600 I = 1, NCOOR
                        IADR = KCS + (I - 1)*NCOOR + ISCOOR - 1
                        NSGN = NINT(WORK(IADR))
                        IF (NSGN .NE. 0) THEN
                           NB = NB + 1
                           NBASE(NB) = I
                           NSIGN(NB) = NSGN
                        END IF
  600                CONTINUE
                     IF (IPREAD .GT. 0) THEN
                     IF (NB .EQ. 1) THEN
                        WRITE (LUPRI,'(1X,I3,3X,A,2X,A,3X,I2)') ISCOOR,
     &                     NAMEX(3*IATOM)(1:4), CHRXYZ(-ICOOR), NBASE(1)
                     ELSE
                        WRITE (LUPRI,'(1X,I3,3X,A,2X,A,3X,A,I2,'
     &                   //CHRNOS(NB-1)//'(2X,A,1X,I2),A,I1)')
     &                   ISCOOR,NAMEX(3*IATOM)(1:4),CHRXYZ(-ICOOR),' [',
     &                   NBASE(1),(CHRSGN(NSIGN(I)),NBASE(I),I=2,NB),
     &                   ' ]/',NB
                     END IF
                     END IF
                  END IF
  500          CONTINUE
  400       CONTINUE
         END IF
  300    CONTINUE
      END IF
      CALL MEMREL('PRICAR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Print internuclear distances (if not an atom)
C
      IF (IPREAD .GT. 0 .AND. NUCDEP .GT. 1)
     &   CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,.FALSE.,WORK,LWORK)
C
      IF (IPREAD .GT. 0) WRITE(LUPRI,'(A,F18.12)')
     &   '  Nuclear repulsion energy :',POTNUC
      RETURN
      END
C  /* Deck priord */
      SUBROUTINE PRIORD(ALPHA,CPRIMU,NPRI,NRCI,SEG,KPRIM,DSM)
#include "implicit.h"
#include "maxorb.h"
      LOGICAL SEG
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
C
C     *****************************
C     *** Segmented contraction ***
C     *****************************
C
      IDONE = 0
      IF (SEG) THEN
         DO 100 ICONTR = 1, NRCI
            ISTART = IDONE + 1
            NLEFT  = NPRI - IDONE
C
C           Find first primitive
C           ====================
C
            IMXA = ISTART + IDAMAX(NLEFT,ALPHA(ISTART),1) - 1
            CALL DSWAP(1,ALPHA(ISTART),1,ALPHA(IMXA),1)
            CALL DSWAP(NRCI,CPRIMU(ISTART,1),KPRIM,CPRIMU(IMXA,1),
     &                 KPRIM)
C
C           Find corresponding contracted function
C           ======================================
C
            IMXC = IDAMAX(NRCI,CPRIMU(ISTART,1),KPRIM)
            CALL DSWAP(NPRI,CPRIMU(1,ICONTR),1,CPRIMU(1,IMXC),1)
C
C           Collect other primitives contributing to same contracted
C           ========================================================
C
            IPRI = 1
            DO 200 I = ISTART + 1, NPRI
               IF (ABS(CPRIMU(I,ICONTR)) .GT. DSM) THEN
                  CALL DSWAP(1,ALPHA(I),1,ALPHA(ISTART+IPRI),1)
                  CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                            CPRIMU(ISTART+IPRI,1),KPRIM)
                  IPRI = IPRI + 1
               END IF
  200       CONTINUE
C
C           Sort primitives
C           ===============
C
            IF (IPRI .GT. 2) THEN
               DO 300 I = ISTART + 1, ISTART + IPRI - 2
                  DO 400 J = I + 1, ISTART + IPRI - 1
                     IF (ALPHA(J) .GT. ALPHA(I)) THEN
                        CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
                        CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                                  CPRIMU(J,1),KPRIM)
                     END IF
  400             CONTINUE
  300          CONTINUE
            END IF
C
            IDONE = IDONE + IPRI
  100    CONTINUE
C
C     ***************************
C     *** General contraction ***
C     ***************************
C
      ELSE
         DO 500 I = 1, NPRI - 1
            DO 600 J = I + 1, NPRI
            IF(ALPHA(J) .GT. ALPHA(I)) THEN
               CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
               CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,CPRIMU(J,1),KPRIM)
           END IF
  600      CONTINUE
  500   CONTINUE
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucsiz */
      SUBROUTINE NUCSIZ(CHARGE,DELTA)
C.......................................................................
C
C     For a Nucleus of Atomic no. "CHARGE", Calculate the Exponent DELTA
C     for a Gaussian Charge Distribution Approximation to the Finite
C     Nucleus. Uses the Formula and Nuclear Masses from
C     REHE Newsletter No. 13 (14. June 1995), L. Visscher and K. Dyall.
C     For Charge.GT.109 Nuclear Mass = CHARGE*2.556.
C     Conversion Factor Bohr/fm 52917.7249.
C                                      T. Saue and J.K. Laerdahl 12.02.1996
C.......................................................................
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP3=1.0D0/3.0D0)
      PARAMETER (FAC=1.5D10,AF=0.529177249D0,BF=0.836D0,CF=0.570D0)
      PARAMETER (FAC2=2.556D0)
      DIMENSION ASYMB(109), AMASS(109)
      CHARACTER*3 ASYMB
      DATA (ASYMB(I),I = 1,109)
     1/'H  ', 'He ', 'Li ', 'Be ', 'B  ',
     2 'C  ', 'N  ', 'O  ', 'F  ', 'Ne ',
     3 'Na ', 'Mg ', 'Al ', 'Si ', 'P  ',
     4 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ',
     5 'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ',
     6 'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ',
     7 'Ga ', 'Ge ', 'As ', 'Se ', 'Br ',
     8 'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ',
     9 'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ',
     O 'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ',
     1 'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ',
     2 'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ',
     3 'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ',
     4 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ',
     5 'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ',
     6 'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ',
     7 'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ',
     8 'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ',
     9 'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ',
     O 'Cm ', 'Bk ', 'Cf ', 'Es ', 'Fm ',
     1 'Md ', 'No ', 'Lw ', 'Uuq', 'Uup',
     2 'Uuh', 'Uus', 'Uuo', 'Uue'/
      DATA (AMASS(I),I = 1,109)
     1/ 1.D0, 4.D0, 7.D0, 9.D0, 11.D0,
     2  12.D0, 14.D0, 16.D0, 19.D0, 20.D0,
     3  23.D0, 24.D0, 27.D0, 28.D0, 31.D0,
     4  32.D0, 35.D0, 40.D0, 39.D0, 40.D0,
     5  45.D0, 48.D0, 51.D0, 52.D0, 55.D0,
     6  56.D0, 59.D0, 58.D0, 63.D0, 64.D0,
     7  69.D0, 74.D0, 75.D0, 80.D0, 79.D0,
     8  84.D0, 85.D0, 88.D0, 89.D0, 90.D0,
     9  93.D0, 98.D0, 98.D0, 102.D0, 103.D0,
     O  106.D0, 107.D0, 114.D0, 115.D0, 120.D0,
     1  121.D0, 130.D0, 127.D0, 132.D0, 133.D0,
     2  138.D0, 139.D0, 140.D0, 141.D0, 144.D0,
     3  145.D0, 152.D0, 153.D0, 158.D0, 159.D0,
     4  162.D0, 162.D0, 168.D0, 169.D0, 174.D0,
     5  175.D0, 180.D0, 181.D0, 184.D0, 187.D0,
     6  192.D0, 193.D0, 195.D0, 197.D0, 202.D0,
     7  205.D0, 208.D0, 209.D0, 209.D0, 210.D0,
     8  222.D0, 223.D0, 226.D0, 227.D0, 232.D0,
     9  231.D0, 238.D0, 237.D0, 244.D0, 243.D0,
     O  247.D0, 247.D0, 251.D0, 252.D0, 257.D0,
     1  258.D0, 259.D0, 262.D0, 261.D0, 262.D0,
     2  263.D0, 262.D0, 265.D0, 266.D0/
C
      INUC = CHARGE+0.2D0
      A = CHARGE*FAC2
      IF (INUC.LE.109) A = AMASS(INUC)
      A = A**DP3
      DELTA = (AF/(BF*A+CF))
      DELTA = FAC*DELTA*DELTA
C      IF (INUC.LE.109) THEN
C         WRITE(LUPRI,'(2X,A,I4,A,A3,A)')
C     +        'Element', INUC, ', ',ASYMB(INUC), ':'
C         WRITE(LUPRI,'(10X,A,F6.2)')
C     +        'Nuclear Mass:', AMASS(INUC)
C      ELSE
C         WRITE(LUPRI,'(2X,A,I4,A)')
C     +        'Element', INUC, ':'
C         WRITE(LUPRI,'(10X,A,F6.2)')
C     +        'Nuclear Mass:', CHARGE*FAC2
C      ENDIF
C      WRITE(LUPRI,'(10X,A,D20.11,/)')
C     +     'Exponent for Gaussian Charge Distribution : ',
C     +     DELTA
      RETURN
      END
C  /* Deck sphlab */
      SUBROUTINE SPHLAB(IORDER,LABINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      CHARACTER LABINT(MXQN**2)*4
      CHARACTER SPDCAR*1, SIGNJ*1
#include "chrnos.h"
C
      IOFF = 0
      DO 100 I = 0, IORDER
         IF (I .EQ. 0) THEN
            LABINT(1) = '1s  '
            IOFF = IOFF + 1
         ELSE IF (I .EQ. 1) THEN
            LABINT(2) = '2px '
            LABINT(3) = '2py '
            LABINT(4) = '2pz '
            IOFF = IOFF + 3
         ELSE
            DO 200 J = -I, I
               IADR = IOFF + J + I + 1
               IF (J .LT. 0) THEN
                  SIGNJ = '-'
               ELSE IF (J .EQ. 0) THEN
                  SIGNJ = ' '
               ELSE
                  SIGNJ = '+'
               END IF
               LABINT(IADR) = CHRNOS(I+1)//SPDCAR(I)//
     &                        CHRNOS(ABS(J))//SIGNJ
  200       CONTINUE
            IOFF = IOFF + 2*I + 1
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck carlab */
      SUBROUTINE CARLAB(IORDER,LABINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION IX(MXAQN), IY(MXAQN), IZ(MXAQN)
      CHARACTER LABINT(MXQN*(MXQN+1)*(MXQN+2)/6)*4
      CHARACTER SPDCAR*1
#include "chrnos.h"
#include "ibtfun.h"
      J = 0
      DO 100 I = 0, IORDER
         IF (I .LE. 3) THEN
            J = J + (I + 1)*(I + 2)/2
            CALL SETCLB(I,LABINT,MXQN*(MXQN+1)*(MXQN+2)/6)
         ELSE
            CALL LMNVAL(I+1,(I + 1)*(I + 2)/2,IX,IY,IZ)
            DO 200 K = 1, (I + 1)*(I + 2)/2
               J  = J + 1
               NX = IX(K)
               NY = IY(K)
               NZ = IZ(K)
               LABINT(J) = SPDCAR(I)//CHRNOS(NX)//CHRNOS(NY)//CHRNOS(NZ)
  200       CONTINUE
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck aotoao */
      SUBROUTINE AOTOAO(IAOAO,IPREAD)
C
C     This subroutine sets up pointer IAOAO which convertes between
C     two different orderings of AO's: from the ordering
C     of MOLECULE to an ordering in which the outer loop is over
C     atoms. (This orderings are identical when no symmetry is used.)
C
C     The purpose of this ordering is to make the AO's appear in the
C     same order as in the corresponding run with no symmetry.
C
C     tuh 120988
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      DIMENSION IAOAO(MXCORB)
C
#include "symmet.h"
#include "nuclei.h"
#include "shells.h"
C
#include "ibtfun.h"
C
      IF (IPREAD .GE. 10) CALL TITLER('Output from AOTOAO','*',103)
C
C     Loop over all atoms
C
      IAO   = 0
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         ISTABA = ISTBNU(IATOM)
         DO 200 ISYMOP = 0, MAXREP
         IF (IBTAND(ISYMOP,ISTABA) .EQ. 0) THEN
            JATOM = JATOM + 1
C
C           Loop over all orbitals
C
            JAO = 0
            DO 300 ISHELL = 1, KMAX
               ICENT  = NCENT(ISHELL)
               ISTABO = ISTBAO(ISHELL)
               DO 400 ICMP = 1, KHKT(ISHELL)
                  DO 500 JSYMOP = 0, MAXREP
                  IF (IBTAND(JSYMOP,ISTABO) .EQ. 0) THEN
                     JAO = JAO + 1
                     JCENT = NUCNUM(ICENT,JSYMOP+1)
                     IF (JATOM .EQ. JCENT) THEN
                        IAO = IAO + 1
                        IAOAO(JAO) = IAO
                     END IF
                  END IF
  500             CONTINUE
  400          CONTINUE
  300       CONTINUE
C
C           End loop over orbitals
C
         END IF
  200    CONTINUE
  100 CONTINUE
C
C     End loop over atoms
C
      IF (IPREAD .GE. 10) THEN
         CALL HEADER('I - IAOAO',6)
         DO 600 I = 1, IAO
            WRITE (LUPRI,'(4X,I5,1X,I5)') I, IAOAO(I)
  600    CONTINUE
      END IF
      RETURN
      END
C/* Deck cntao */
      SUBROUTINE CNTAO(IPRINT)
C 
C     tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
C
#include <symmet.h>
#include <nuclei.h>
#include <shells.h>
C
#include <ibtfun.h>
C
      IORB = 0
      DO 100 LA = 0, MAXREP
      DO 200 IA = 1, KMAX
         DO 300 NA = 1, KHKT(IA)
            IVARB = IBTXOR(LA,ISYMAO(NHKT(IA),NA))
            IF (IBTAND(ISTBAO(IA),IVARB) .EQ. 0) THEN
               IORB = IORB + 1
               ICNTAO(IORB)  = NCENT(IA)
            END IF
  300    CONTINUE
  200 CONTINUE
  100 CONTINUE
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ICNTAO',-1)
         DO 400 I = 1, IORB
            WRITE (LUPRI,'(4X,I5,1X,I5)') I, ICNTAO(I)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck magcor */
      SUBROUTINE MAGCOR(IPREAD)
C
C     Printing of symmetrized nuclear magnetic moments
C
C     tuh March 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "chrxyz.h"
#include "chrsgn.h"
C
      INTEGER NBASE(8), NSIGN(8), NRTREP(0:7)
C
#include "cbiher.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
C
#include "ibtfun.h"
C
      IF (NMRISS .AND. MAXREP .GT. 0) THEN
CTROND         IF (IPREAD .GT. 0) THEN
         IF (IPREAD .GE. 0) THEN
            WRITE (LUPRI,'(/)')
            CALL HEADER('Symmetrized nuclear magnetic moments',1)
         END IF
         DO 100 IREP = 0, MAXREP
            NRTREP(IREP) = 0
            DO 110 ICENT = 1, NUCIND
               DO 120 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR,IREP,2)
                  IF (ISCOOR .GT. 0) NRTREP(IREP) = NRTREP(IREP) + 1
 120           CONTINUE
 110        CONTINUE
 100     CONTINUE
CTROND         IF (IPREAD .GT. 0) THEN
         IF (IPREAD .GE. 0) THEN
            WRITE (LUPRI,'(1X,A,8I3)')
     &           ' Number of components in each symmetry: ',
     &           (NRTREP(I),I=0,MAXREP)
         END IF
C
         DO 200 IREP = 0, MAXREP
         IF (NRTREP(IREP) .GT. 0) THEN
CTROND            IF (IPREAD .GT. 0) THEN
            IF (IPREAD .GE. 0) THEN
              WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1/)')
     &           'Symmetry',REP(IREP),'(',IREP,')'
            END IF
            DO 210 ICENT = 1, NUCIND
               DO 220 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR,IREP,2)
                  IF (ISCOOR .GT. 0) THEN
                     IVAR = IBTXOR(IREP,ISYMAX(ICOOR,2))
                     NB = 0
                     DO 230 ISYMOP = 0, MAXOPR
                     IF (IBTAND(ISYMOP,ISTBNU(ICENT)) .EQ. 0) THEN
                        NB = NB + 1
                        NSIGN(NB) = NINT(PT(IBTAND(ISYMOP,IVAR)))
                        NBASE(NB) = 3*(NUCNUM(ICENT,ISYMOP+1)-1) + ICOOR
                     END IF
 230                 CONTINUE
CTROND                     IF (IPREAD .GT. 0) THEN
                     IF (IPREAD .GE. 0) THEN
                        WRITE (LUPRI,
     &                  '(2X,A,I3,3X,A,2X,A,3X,I2,7(2X,A,1X,I2))')
     &                  'I',ISCOOR, NAMEX(3*ICENT)(1:4), CHRXYZ(-ICOOR),
     &                  NBASE(1), (CHRSGN(NSIGN(I)),NBASE(I),I=2,NB)
                     END IF
                  END IF
  220          CONTINUE
  210       CONTINUE
         END IF
  200    CONTINUE
      END IF
      RETURN
      END
C  /* Deck tracr */
      SUBROUTINE TRACR(NCOOR)
C
C     Set up IPTCOR - points from symmetry coordinate to generating
C     Cartesian coordinate. Extracted from old TRACOR routine by
C     K.Ruud, Dec-96
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
C
      DO 400 ITYPE = 1, 2
         IJ = 0
         DO 410 J = 0, MAXREP
            DO 420 I = 1, NCOOR
               IF (IPTCNT(I,J,ITYPE) .NE. 0) THEN
                  IJ = IJ + 1
                  IPTCOR(IJ,ITYPE) = I
               END IF
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
      RETURN
      END
C  /* Deck setclb */
      SUBROUTINE SETCLB(I,LABINT,NDIM)
#include "implicit.h"
      CHARACTER LABINT(NDIM)*4
      IF (I .EQ. 0) THEN
         LABINT(1) = 's   '
      ELSE IF (I .EQ. 1) THEN
         LABINT(2) = 'px  '
         LABINT(3) = 'py  '
         LABINT(4) = 'pz  '
      ELSE IF (I .EQ. 2) THEN
         LABINT(5)  = 'dxx '
         LABINT(6)  = 'dxy '
         LABINT(7)  = 'dxz '
         LABINT(8)  = 'dyy '
         LABINT(9)  = 'dyz '
         LABINT(10) = 'dzz '
      ELSE IF (I .EQ. 3) THEN
         LABINT(11)  = 'fxxx'
         LABINT(12)  = 'fxxy'
         LABINT(13)  = 'fxxz'
         LABINT(14)  = 'fxyy'
         LABINT(15)  = 'fxyz'
         LABINT(16)  = 'fxzz'
         LABINT(17)  = 'fyyy'
         LABINT(18)  = 'fyyz'
         LABINT(19)  = 'fyzz'
         LABINT(20)  = 'fzzz'
      END IF
      RETURN
      END
C  /* Deck spdcar */
      FUNCTION SPDCAR(IORDER)
#include "implicit.h"
      CHARACTER*1 SPDCAR
      IF (IORDER .EQ. 0) THEN
         SPDCAR = 's'
      ELSE IF (IORDER .EQ. 1) THEN
         SPDCAR = 'p'
      ELSE IF (IORDER .EQ. 2) THEN
         SPDCAR = 'd'
      ELSE IF (IORDER .LE. 6) THEN
         SPDCAR = CHAR(ICHAR('f') + IORDER - 3)
      ELSE
C
C        Note: j is not used - therefore this special case
C
         SPDCAR = CHAR(ICHAR('k') + IORDER - 7)
      END IF
      RETURN
      END
C  /* Deck sphinp */
      SUBROUTINE SPHINP(LUINFO,WORK,LWORK,DOOWN)
C
C     Calculate spherical harmonic transformation coefficients
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1 = 1.0D0)
#include "mxcent.h"
#include "maxaqn.h"
#include "sphtrm.h"
#include "molinp.h"
#include "ccom.h"
#include "cbirea.h"
      CHARACTER*1 KEY
      LOGICAL DOOWN
      DIMENSION WORK(LWORK)
C
      IF (IPREAD .GT. 4) THEN
         WRITE(LUPRI,'(1X,A)') ' Calculation of transformation '
     &                         //'matrix for spherical integrals.'
      END IF
C
      IOFF = 1
      DO 100 I = 1, MXQN
         KCK(I)    = I*(I+1)/2
         ISPADR(I) = IOFF
         IOFF      = IOFF + KCK(I)**2
         IF (DOOWN) THEN
            SPH(I) = .TRUE.
         ELSE
            IF (DOCART) THEN
               KHK(I) = I*(I+1)/2
               SPH(I) = .FALSE.
            ELSE
               KHK(I) = 2*I - 1
               SPH(I) = I .GT. 2
            END IF
         END IF
  100 CONTINUE
C
      SPHNRM = .TRUE.
C
      IF (DOCART) THEN
         CALL CARPOW
         CALL CARLAB(MXQN-1,GTOTYP)
      ELSE
         CALL DZERO(CSP,NCSP)
         IF (.NOT.DOOWN) THEN
            CALL SPHLAB(MXQN-1,GTOTYP)
            KCSPT = 1
            KWRK  = KCSPT + (MXQN*(MXQN + 1)/2)*(2*MXQN - 1)
            IF (KWRK .GT. LWORK) CALL STOPIT('SPHINP',' ',KWRK,LWORK)
            LWRK  = LWORK - KWRK + 1
            DO 200 LVAL = 0, MXQN - 1
               NLM = 2*LVAL + 1
               NXYZ = (LVAL+1)*(LVAL+2)/2
               IF (LVAL .EQ. 0) THEN
                  CSP(1)  = D1
               ELSE IF (LVAL .EQ. 1) THEN
                  CSP (2) = D1
                  CSP (6) = D1
                  CSP(10) = D1
               ELSE
                  MINTEG = 1
                  SPHNRM = .FALSE.
                  CALL SPHCOM(LVAL,WORK(KCSPT),NLM,NXYZ,1,MINTEG,
     &                        WORK(KWRK),LWRK,IPREAD)
                  CALL MTRSP(NXYZ,NLM,WORK(KCSPT),NXYZ,
     &                       CSP(ISPADR(LVAL+1)),NLM)
               END IF
  200       CONTINUE
         ELSE
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(I5)') MAXSPD
            IOFF = 1
            DO 300 I = 1, MAXSPD
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),'(4X,A1,I5)')  KEY, KHK(I)
               CALL RDCSP(LUINFO,CSP(ISPADR(I)),GTOTYP(IOFF),
     &           KHK(I),KCK(I))
               IOFF = IOFF + KHK(I)
  300       CONTINUE
         END IF
         IF (DOOWN .OR. IPREAD .GT. 4) THEN
            CALL HEADER('Cartesian transformation matrix',-1)
            DO 400 I = 1, MXQN
               KHKI = KHK(I)
               KCKI = KCK(I)
               WRITE (LUPRI,'(1X,A,I5)')
     &               ' Coefficients for angular quantum number ',I-1
               CALL OUTPUT(CSP(ISPADR(I)),1,KHKI,1,KCKI,KHKI,KCKI,1,
     &                     LUPRI)
 400        CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck rdcsp */
      SUBROUTINE RDCSP(LUINFO,CSP,GTOTYP,KHKI,KCKI)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
      CHARACTER*4 GTOTYP(*)
      DIMENSION CSP(KHKI,KCKI)
      DO 100 I = 1, KHKI
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(1X,A,15F5.3)')
     &          GTOTYP(I), (CSP(I,J),J=1,KCKI)
  100 CONTINUE
      RETURN
      END
C  /* Deck rdline */
      SUBROUTINE RDLINE(LUINFO)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "molinp.h"
#include "infpar.h"
C
C     Reads a line from MOLECULE.INP
C
C     We keep track of which line we've read last time (NMLN_LASTREAD).
C     This allows some routines to read a line, check if it want to use
C     that line, and then discard it by decrementing NMLINE.
C     Note that we only allow decrements of one /jth-20000621.
C     (This is used to check if SMALL basis keyword has been specified,
C     or if default - kinetic balance - should be used.)
C
      DATA NMLN_LASTREAD /0/
      SAVE NMLN_LASTREAD
C
C     If NMLINE is zero, reset NMLN_LASTREAD since this is a new readin,
C     e.g. during a geometry optimization.
C
      IF ( NMLINE .EQ. 0) NMLN_LASTREAD = 0
C
   10 CONTINUE
      NMLINE = NMLINE + 1
      IF ( NMLINE .LT. NMLN_LASTREAD ) THEN
         WRITE(LUPRI,'(//A,2(/A,I5))')
     $        '*** ERROR ***: NMLINE < NMLN_LASTREAD in RDLINE',
     $        'NMLINE = ',NMLINE,
     $        'NMLN_LASTREAD = ',NMLN_LASTREAD
         CALL QUIT('*** ERROR reading basis')
      ENDIF
      IF (.NOT. PARHER .AND. NMLINE .GT. NMLN_LASTREAD)
     $     READ (LUINFO,'(A)',END=11) MLINE(NMLINE)
      IF (NMLINE.GT.KMLINE) THEN
         WRITE (LUPRI,'(//A,/A,I5)')
     &        ' READIN ERROR: Too many lines in input -,',
     &        ' Increase dimension in /MOLINP/ to at least',NMLINE
         CALL QUIT('*** ERROR *** Too many input lines for READIN')
      END IF
      IF((MLINE(NMLINE)(1:1).EQ.'#') .OR.
     $     (MLINE(NMLINE)(1:1).EQ.'!') .OR.
     $     (MLINE(NMLINE)(1:1) .EQ. '$' ))
     &     GOTO 10
      NMLN_LASTREAD = NMLINE
      RETURN
   11 MLINE(NMLINE)='        '
      RETURN
      END
C  /* Deck wronel */
      SUBROUTINE WRONEL(TITLE,NONTYP,NONT,IQM,IFXYZ,KATOM,JCO2,KANG)
C
C     Write interface records
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      CHARACTER*72 TITLE(2)
      CHARACTER*8 CDATE,CTIME
      CHARACTER*8 TABLE1(4), TABLE2(4), TABLE3(4), TABLE4(4), TABLE5(4)
      DIMENSION QPOL(6), QQ(3), IFXYZ(3), JFXYZ(3)
      DIMENSION CHRN(MXCENT), COOO(MXCENT,3), IPRXYZ(MXCORB),
     &          INAMN(MXCORB), IGTO(MXCORB), RTITLE(24)
      DIMENSION NONT(KATOM), IQM(KATOM), JCO2(KANG,KATOM)
      DATA TABLE1 /'********','        ','        ','ISORDK  '/,
     &     TABLE2 /'********','        ','        ','SCFINP  '/,
     &     TABLE3 /'********','        ','        ','SYMINPUT'/,
     &     TABLE4 /'********','        ','        ','TESTDATA'/,
     &     TABLE5 /'********','        ','        ','SPHERICA'/
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "hrunit.h"
#include "nuclei.h"
#include "frame.h"
#include "sphtrm.h"
#include "ibtfun.h"
C
      CALL GETDAT(CDATE,CTIME)
      TABLE1(2) = CDATE
      TABLE2(2) = CDATE
      TABLE3(2) = CDATE
      TABLE4(2) = CDATE
      TABLE5(2) = CDATE
      TABLE1(3) = CTIME
      TABLE2(3) = CTIME
      TABLE3(3) = CTIME
      TABLE4(3) = CTIME
      TABLE5(3) = CTIME
C
      READ (TITLE(1),'(12A6)') (RTITLE(I),I= 1,12)
      READ (TITLE(2),'(12A6)') (RTITLE(I),I=13,24)
C
      ICENT = 0
      DO 100 N = 1, NUCIND
      DO 100 IREP = 0,MAXREP
         IF (IBTAND(ISTBNU(N),IREP) .EQ. 0) THEN
            ICENT = ICENT + 1
            CHRN(ICENT)   = CHARGE(N)
            COOO(ICENT,1) = PT(IBTAND(ISYMAX(1,1),IREP))*CORD(1,N)
            COOO(ICENT,2) = PT(IBTAND(ISYMAX(2,1),IREP))*CORD(2,N)
            COOO(ICENT,3) = PT(IBTAND(ISYMAX(3,1),IREP))*CORD(3,N)
      END IF
  100 CONTINUE
C
      I = 0
      DO 200 ISYMOP = 0, MAXREP
         NBI = NAOS(ISYMOP+1)
         DO 220 L = 1,NBI
            I = I + 1
            READ(  NAMN(IPCEN(I)),'(A4)') INAMN(I)
            READ(GTOTYP(IPTYP(I)),'(A4)') IGTO(I)
  220    CONTINUE
  200 CONTINUE
C
C     0 - unlabeled first record
C     ==========================
C
      REWIND LUONEL
      WRITE (LUONEL) RTITLE,MAXREP+1,(NAOS(I),I=1,MAXREP+1), POTNUC
C
C     1 - ISORDK
C     ==========
C
      WRITE (LUONEL) TABLE1
      WRITE (LUONEL) DUMMY, DUMMY, DUMMY, DUMMY
      WRITE (LUONEL) CHRN, NUCDEP, COOO
C
C     2 - SCFINP
C     ==========
C
      WRITE (LUONEL) TABLE2
      WRITE (LUONEL) RTITLE,
     &               MAXREP+1, (NAOS(I), I = 1,MAXREP+1),
     &               POTNUC,
     &               KMAX, (NCENT(I), I = 1,KMAX),
     &               NBASIS, ( JTRAN(I),             I = 1,NBASIS),
     &                       ((ITRAN(I,J), J = 1,8), I = 1,NBASIS),
     &                       ((CTRAN(I,J), J = 1,8), I = 1,NBASIS),
     &               NBASIS, ( INAMN(I),             I = 1,NBASIS),
     &                       (IPRXYZ(I),             I = 1,NBASIS),
     &               DPNUC,
     &               NUCDEP, ((COOO(I,J),J=1,3),I=1,NUCDEP),
CTROND: JFXYZ, QPOL and QQ is not defined !!!
     &               IFXYZ, DUMMY, QPOL, QQ, JFXYZ
C
C     3 - SYMINPUT
C     ============
C
      WRITE (LUONEL) TABLE3
      WRITE (LUONEL) NBASIS,(INAMN(I),I=1,NBASIS),(IGTO(I),I=1,NBASIS),
     &               POTNUC,DUMMY,DUMMY,DUMMY
C
C     4 - TESTDATA
C     ============
C
C     Information for testing magnetic integrals
C
      WRITE (LUONEL) TABLE4
      WRITE (LUONEL) NONTYP, (NONT(I), I = 1, NONTYP),
     &               (IQM(I), (JCO2(I,J),J=1,IQM(I)), I = 1, NONTYP),
     &               (IDUMMY,I=1,5)
C
C     5 - SPHERICA
C     ============
C
C     Information for testing spherical transformation
C
      IDOCRT = 1
      IF (.NOT.DOCART) IDOCRT = 0
      WRITE (LUONEL) TABLE5
      WRITE (LUONEL) KMAX, IDOCRT, (NHKT(I),I = 1,KMAX), (IDUMMY,I=1,5)
      WRITE (LUONEL) MAXREP, (ISTBAO(I),I=1,KMAX), (IDUMMY,I=1,6)
      WRITE (LUONEL) MXQN, MXAQN, ((ISYMAO(I,J),I=1,MXQN),J=1,MXAQN)
      WRITE (LUONEL) (CSP(I),I=1,NCSP),(ISPADR(I),I=1,MXQN)
C
      RETURN
      END
C  /* Deck reapri */
      SUBROUTINE REAPRI
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"
C
#include "abainf.h"
#include "cbirea.h"
C
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
C
      CALL HEADER('Test output from REAPRI',-1)
      WRITE (LUPRI,'(A,1P,D15.6)') ' THRS  ', THRS
      WRITE (LUPRI,'(A,I5)') ' NHTYP ', NHTYP
      WRITE (LUPRI,'(A,L5)') ' DOCART', DOCART
      print *, 'nucind,nucdep',NUCIND, NUCDEP
      print *, 'charge ',(CHARGE(i),i=1,nucind)
      print *, 'cord1',(cord(1,i),i=1,nucind)
      print *, 'cord2',(cord(2,i),i=1,nucind)
      print *, 'cord2',(cord(3,i),i=1,nucind)
      print *, 'noorbt ',(noorbt(i),i=1,nucind)
      print *, 'nucpre ',(nucpre(i),i=1,nucind)
      print *, 'nucnum ',((NUCNUM(i,j),i=1,nucind),j=1,8)
      print *, 'nucdeg ',(nucdeg(i),i=1,nucind)
      print *, 'istbnu ',(istbnu(i),i=1,nucind)
      print *, 'ntraco,itraco',NTRACO, (ITRACO(i),i=1,3)
      print *, 'NATOMS, NFLOAT, NBASIS, NPBAS',
     &          NATOMS, NFLOAT, NBASIS, NPBAS
      print *, 'namex ',(namex(i),i=1,3*nucind)
      print *, 'namdep',(namdep(i),i=1,nucdep)
      print *, 'namdpx',(namdpx(i),i=1,nucdep)
      call header('priexp',-1)
      call output(priexp,1,1,1,mxprim,1,mxprim,1,lupri)
      call header('priccf',-1)
      call output(priccf,1,mxprim,1,mxcont,mxprim,mxcont,1,lupri)
      call header('pricrx',-1)
      call output(pricrx,1,mxprim,1,3,mxprim,3,1,lupri)
      print *, ' kmax ', kmax
      call header('cent',-1)
      call output(cent,1,kmax,1,3,kmax,3,1,lupri)
      print *, ' nhkt ',(nhkt(i),i=1,kmax)
      print *, ' khkt ',(khkt(i),i=1,kmax)
      print *, ' kckt ',(kckt(i),i=1,kmax)
      print *, ' istbao ',(istbao(i),i=1,kmax)
      print *, ' nuco ',(nuco(i),i=1,kmax)
      print *, ' jstrt ',(jstrt(i),i=1,kmax)
      print *, ' nstrt ',(nstrt(i),i=1,kmax)
      print *, ' ncent ',(ncent(i),i=1,kmax)
      print *, ' share ',(share(i),i=1,kmax)
      print *, ' nrco ',(nrco(i),i=1,kmax)
      print *, ' numcf ',(numcf(i),i=1,kmax)
      print *, ' nbch ',(nbch(i),i=1,kmax)
      print *, ' kstrt ',(kstrt(i),i=1,kmax)
      print *, ' segm ',(segm(i),i=1,kmax)
      print *, ' iptshl ',(iptshl(i),i=1,kmax)
      print *, ' numcft ',(numcft(i),i=1,kmax)
      print *, ' sphr ',(sphr(i),i=1,kmax)
      print *, ' fmult ',(fmult(i),i=0,7)
      print *, ' pt ',(pt(i),i=0,7)
      print *, ' mult ',(mult(i),i=0,7)
      print *, ' maxrep, maxopr ',MAXREP, MAXOPR
      print *, ' isymax ',((isymax(i,j),i=1,3),j=1,2)
      print *, ' isymao ',((isymao(i,j),i=1,mxqn),j=1,mxaqn)
      print *, ' nparsu ',(nparsu(i),i=1,8)
      print *, ' naos ',(naos(i),i=1,8)
      print *, ' nparnu ',((nparnu(i,j),i=1,8),j=1,8)
      print *, ' iptsym ',((iptsym(i,j),i=1,norbs),j=0,7)
      print *, ' iptcnt ',(((iptcnt(i,j,k),i=1,3*nucdep),j=0,7),k=1,2)
      print *, ' ncrrep ',((ncrrep(i,j),i=0,7),j=1,2)
      print *, ' iptcor ',((iptcor(i,j),i=1,3*nucdep),j=1,2)
      print *, ' naxrep ',((naxrep(i,j),i=0,7),j=1,2)
      print *, ' iptax ',((iptax(i,j),i=1,3),j=1,2)
      print *, ' iptxyz ',(((iptxyz(i,j,k),i=1,3),j=0,7),k=1,2)
      print *, ' iptnuc ',((iptnuc(i,j),i=1,nucdep),j=0,7)
      call header('ctran',-1)
      call output(ctran,1,nbasis,1,8,mxcorb,8,1,lupri)
      print *, ' itran ',((itran(i,j),i=1,nbasis),j=1,8)
      print *, ' iaoao ',(iaoao(i),i=1,nbasis)
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck symgrp */
      SUBROUTINE SYMGRP(IGEN,NSYMOP,CLASS)
C*****************************************************************************
C
C     Given the group generators, this subroutine will identify
C     Abelian subgroup and set up group multiplication table,
C     character table and direct product table. Irreps are identified.
C
C     tsaue - august 1994
C     Sep 24 1996 - tsaue : Included GROUPS,JSOP and IPAR in COMMON SYMMET.
C                           Modified activation of LSYMOP
C                           GROUPS and SYMOP are initialized in BLOCK DATA
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "cbirea.h"
      LOGICAL LSYMOP(0:7)
      DIMENSION IGEN(3),IROTS(3),IREFL(3),JPAR(0:7)
      CHARACTER*(*) CLASS
      DATA IROTS/3,5,6/
      DATA IREFL/4,2,1/
      DATA JPAR/ 1,-1,-1, 1,-1, 1, 1,-1/
      EXTERNAL BDPGROUP
#include "ibtfun.h"
      DO 5 I = 0,7
        LSYMOP(I) = .FALSE.
    5 CONTINUE
C
C     Activate all symmetry operations of the group
C     =============================================
C
      LSYMOP(0) = .TRUE.
      JSOP(0)   = 0
      IPAR(0)   = 1
      DO 10 I   = 1,MAXREP
        I0 = IBTAND(1,I)          *IGEN(1)
        I1 = IBTAND(1,IBTSHR(I,1))*IGEN(2)
        I2 = IBTAND(1,IBTSHR(I,2))*IGEN(3)
        IND = IBTXOR(IBTXOR(I0,I1),I2)
        LSYMOP(IND) = .TRUE.
        IPAR(I)     = JPAR(IND)
   10 CONTINUE
C
C     List group operations in preferred order,
C       that is construct pointer LSOP
C     =========================================
C
C     Identity
C
      IND = 0
      JSOP(IND) = 0
C
C     Rotations
C
      NROTS = 0
      DO 40 I = 1,3
      IF(LSYMOP(IROTS(I))) THEN
        IND         = IND + 1
        JSOP(IND) = IROTS(I)
        NROTS       = NROTS + 1
      ENDIF
   40 CONTINUE
C
C     Inversion
C
      NINVC = 0
      IF(LSYMOP(7)) THEN
        IND         = IND + 1
        JSOP(IND) = 7
        NINVC       = 1
      ENDIF
C
C     Reflections
C
      NREFL = 0
      DO 50 I = 1,3
      IF(LSYMOP(IREFL(I))) THEN
        IND         = IND + 1
        JSOP(IND) = IREFL(I)
        NREFL       = NREFL + 1
      ENDIF
   50 CONTINUE
      IF(IND.NE.MAXREP) CALL QUIT('SYMGRP:IND.NE.MAXREP!')
C
C     Classify group
C     ==============
C     tsaue - Here I have devised a highly empirical formula, but
C             it works !!!
C
      IGROUP = MIN(7,NINT((4*NROTS+8*NINVC+6*NREFL)/3.0))
      GROUP  = GROUPS(IGROUP)
C
C     Generate character table
C     ========================
C
      DO 60 I = 0,MAXREP
        IXVAL(0,I) = 1
        DO 70 J = 1,NSYMOP
          IXVAL(IGEN(J),I) = PT(IBTAND(IBTSHR(I,(J-1)),1))
          DO 80 K = 1,(J-1)
            IND      = IBTXOR(IGEN(J),IGEN(K))
            IXVAL(IND,I)  = IXVAL(IGEN(J),I)*IXVAL(IGEN(K),I)
            DO 90 L = 1,(K-1)
              IXVAL(IBTXOR(IND,IGEN(L)),I)
     &           = IXVAL(IND,I)*IXVAL(IGEN(L),I)
   90       CONTINUE
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
C
C     Classify irrep
C     ==============
C
      DO 100 I = 0,MAXREP
        REP(I) = 'A  '
        IPOS = 2
C
C       Rotational symmetry
C
        IF(NROTS.EQ.3) THEN
          IND = (1-IXVAL(JSOP(1),I))+(1-IXVAL(JSOP(2),I))/2
          IF(IND.NE.0) THEN
            REP(I)(1:1) = 'B'
            REP(I)(2:2) = CHAR(ICHAR('0')+IND)
            IPOS = 3
          ENDIF
        ELSEIF(NROTS.EQ.1) THEN
          IF(IXVAL(JSOP(1),I).EQ.-1) REP(I)(1:1) = 'B'
          IF(NREFL.EQ.2) THEN
            IF(IBTAND(IBTSHR(JSOP(1),1),1).EQ.1) THEN
              IND = 2
            ELSE
              IND = 3
            ENDIF
            IF(IXVAL(JSOP(IND),I).EQ.1) THEN
              REP(I)(2:2) = '1'
            ELSE
              REP(I)(2:2) = '2'
            ENDIF
          ENDIF
        ELSEIF(NREFL.EQ.1) THEN
C
C       Mirror symmetry
C
          IF(IXVAL(JSOP(1),I).EQ.1) THEN
            REP(I)(2:2) = ''''
          ELSEIF(IXVAL(JSOP(1),I).EQ.-1) THEN
            REP(I)(2:2) = '"'
          ENDIF
        ENDIF
C
C       Inversion symmetry
C
        IF(NINVC.EQ.1) THEN
          IND = NROTS+1
          IF(IXVAL(JSOP(IND),I).EQ.1) THEN
            REP(I)(IPOS:IPOS) = 'g'
          ELSE
            REP(I)(IPOS:IPOS) = 'u'
          ENDIF
        ENDIF
  100 CONTINUE
C
C     Output section
C     ==============
C
      IF(IPREAD.GT.0) THEN
        CALL HEADER('SYMGRP:Point group information',-1)
C
C       Group name
C
        IF (CLASS(1:3) .EQ. 'N/A') THEN
           WRITE(LUPRI,'(A,A3)') 'Point group: ',GROUP
        ELSE
           WRITE(LUPRI,'(A,A)') 'Full group is:  ',CLASS
           WRITE(LUPRI,'(A,A3)') 'Represented as: ',GROUP
        END IF
C
C       Group generators
C
        IF(NSYMOP.GT.0) THEN
          WRITE(LUPRI,'(/3X,A/)')
     &        '* The point group was generated by:'
          DO 110 I = 1,NSYMOP
          IF    (SYMOP(IGEN(I))(1:1).EQ.'C') THEN
            WRITE(LUPRI,'(6X,3A)')
     &        'Rotation about the ',SYMOP(IGEN(I))(3:3),'-axis'
          ELSEIF(SYMOP(IGEN(I))(1:1).EQ.'O') THEN
            WRITE(LUPRI,'(6X,3A)')
     &        'Reflection in the ',SYMOP(IGEN(I))(2:3),'-plane'
          ELSE
            WRITE(LUPRI,'(6X,A)')
     &        'Inversion centre'
          ENDIF
  110     CONTINUE
C
C         Group multiplication table
C
          WRITE(LUPRI,'(/3X,A/)') '* Group multiplication table'
          WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))')'|',(SYMOP(JSOP(I)),
     &        I = 0,MAXREP)
          WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
          DO 120 I = 0,MAXREP
            WRITE(LUPRI,'(4X,A3,1X,A1,8(1X,A3,1X))')
     &        SYMOP(JSOP(I)),'|',
     &        (SYMOP(IBTXOR(JSOP(I),JSOP(J))),J = 0,I)
  120     CONTINUE
        ENDIF
C
C       Character table
C
        WRITE(LUPRI,'(/3X,A/)') '* Character table'
        WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))') '|',(SYMOP(JSOP(J)),
     &               J = 0,MAXREP)
        WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
        DO 130 I = 0,MAXREP
          WRITE(LUPRI,'(4X,A3,1X,A1,8(1X,I3,1X))')
     &      REP(I),'|',(IXVAL(JSOP(J),I),J=0,MAXREP)
  130   CONTINUE
C
C       Direct product table
C
        WRITE(LUPRI,'(/3X,A/)') '* Direct product table'
        WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))')'|',(REP(I),I = 0,MAXREP)
        WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
        DO 140 I = 0,MAXREP
          WRITE(LUPRI,'(3X,1X,A3,1X,A1,8(1X,A3,1X))')
     &      REP(I),'|',(REP(IBTXOR(I,J)),J = 0,I)
  140   CONTINUE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dblgrp */
      SUBROUTINE DBLGRP
C*****************************************************************************
C
C     This routine will analyze the fermion irreps in the molecular symmetry
C     and, if possible, set up a transformation within AO - basis to a
C     symmetry-adapted spinor basis. details are given below.
C
C     Written by T.Saue - october 1994 - odense
C     Sep 24 1996 - tsaue: Updated version
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0,D1 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "cbirea.h"
#include "symmet.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbqal.h"
C
      LOGICAL LBUF(8)
      EXTERNAL BDQALGEBRA
C
#include "ibtfun.h"
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL IZERO(IQTOPQ,32)
      CALL IZERO(IPQTOQ,32)
C
C     Determine number of boson and fermion irreps
C     ============================================
C
      NBSYM = MAXREP + 1
      NFSYM = 1
      IF(NINVC.EQ.1) NFSYM = 2
      IF(NFSYM.EQ.1) THEN
        FREP(1) = 'E1 '
      ELSE
        FREP(1) = 'E1g'
        FREP(2) = 'E1u'
      ENDIF
C
C     Determine NZ - parameter
C     ========================
C       NZ = 4  quaternionic group        no    totally symmetric rotations
C       NZ = 2  complex group             one   totally symmetric rotation
C       NZ = 1  real group                three totally symmetric rotations
C
      NZ = 1
      DO 30 I = 1,3
        IF(ISYMAX(I,2).EQ.0) THEN
          NZ   = NZ + 1
          ITOT = I
        ENDIF
   30 CONTINUE
C
C     Align rotations and quaternion irreps,
C     for complex groups this is important....
C
      IQ1 = 3
      IQ2 = 2
      IQ3 = 1
      IZQ = 2
      IYQ = 3
      IXQ = 4
      IF(NZ.EQ.2) THEN
        IF(ITOT.EQ.1) THEN
          IQ1 = 1
          IQ2 = 3
          IQ3 = 2
          IXQ = 2
          IZQ = 3
          IYQ = 4
        ELSEIF(ITOT.EQ.2) THEN
          IQ1 = 2
          IQ2 = 1
          IQ3 = 3
          IYQ = 2
          IXQ = 3
          IZQ = 4
        ENDIF
      ENDIF
C
C
C     Set up direct product table for fermion ircops
C     ==============================================
C
C     IF there is a totally symmetric irrep, it should be the
C     symmetry of Im(p_1,q_1)......
C
C     Re(p_1,q_1)
      IFDIRP(1,1) = 0
C     Im(p_1,q_1)
      IFDIRP(2,1) = ISYMAX(IQ1,2)
C     Re(p_1,Q_1)
      IFDIRP(3,1) = ISYMAX(IQ2,2)
C     Im(p_1,Q_1)
      IFDIRP(4,1) = ISYMAX(IQ3,2)
C     Re(p_1,q_2)
      IFDIRP(1,2) = IXYZ
C     Im(p_1,q_2)
      IFDIRP(2,2) = ISYMAX(IQ1,1)
C     Re(p_1,Q_2)
      IFDIRP(3,2) = ISYMAX(IQ2,1)
C     Im(p_1,Q_2)
      IFDIRP(4,2) = ISYMAX(IQ3,1)
C
C
C     Determine distribution of boson irreps in spinors
C     =================================================
C
      DO 10 IM = 1,4
        JBTOF(IFDIRP(IM,1),1) = 1
        JSPINR(IM,1,1) = IFDIRP(IM,1)
        JBTOF(IFDIRP(IM,2),2) = 1
        JSPINR(IM,2,1) = IFDIRP(IM,2)
   10 CONTINUE
      IF(NFSYM.EQ.2) THEN
        DO 20 IM = 1,4
          JBTOF(IFDIRP(IM,2),1) = 2
          JSPINR(IM,1,2) = IFDIRP(IM,2)
          JBTOF(IFDIRP(IM,1),2) = 2
          JSPINR(IM,2,2) = IFDIRP(IM,1)
   20   CONTINUE
      ENDIF
C
C     Determine boson symmetries connected by a given
C     fermion irrep
C
      IF(NFSYM.EQ.1) THEN
        DO ISYM = 1,NBSYM
          JFSYM(ISYM,1) = ISYM
          JFSYM(ISYM,2) = ISYM
        ENDDO
      ELSEIF(NFSYM.EQ.2) THEN
        ISYMA = 0
        ISYMB = 0
        DO ISYM = 1,NBSYM
          IF(JBTOF(ISYM-1,1).EQ.1) THEN
            ISYMA = ISYMA + 1
            JFSYM(ISYMA,1) = ISYM
          ELSE
            ISYMB = ISYMB + 1
            JFSYM(ISYMB,2) = ISYM
          ENDIF
        ENDDO
      ENDIF
C
C     Matrix packing
C     ==============
C
      JMROI(1) = 1
      IQDEF(1) = 1
C
C     Quaternion groups
C
      IF(NZ.EQ.4) THEN
        JMROI (2) = 1
        JMROI (3) = 1
        JMROI (4) = 1
        IF(NFSYM.EQ.1) THEN
C       C1:
          DO I = 1,4
            IPQTOQ(I,   0) = I
          ENDDO
        ELSE
C       Ci:
          DO I = 1,4
            IPQTOQ(I,   0) = I
            IPQTOQ(I,IXYZ) = I
          ENDDO
        ENDIF
C
C     Complex groups:
C
      ELSEIF(NZ.EQ.2) THEN
C
C       Totally symmetric Z-rotation
C
          JMROI (2) = 1
          JMROI (3) = 4
          JMROI (4) = 4
C         ..non-totally symmetric rotation
        IRNT           = ISYMAX(IQ2,2)
          IPQTOQ(1,0   ) = 1
          IPQTOQ(2,0   ) = 2
          IPQTOQ(1,IRNT) = 4
          IPQTOQ(2,IRNT) = 3
        IF(NFSYM.EQ.2) THEN
          IRNU = IBTXOR(IRNT,IXYZ)
          DO IZ = 1,NZ
            IPQTOQ(IZ,IXYZ) = IPQTOQ(IZ,0)
            IPQTOQ(IZ,IRNU) = IPQTOQ(IZ,IRNT)
          ENDDO
        ENDIF
      ELSE
        JMROI(2)  = 2
        JMROI(3)  = 3
        JMROI(4)  = 4
        IF(NFSYM.EQ.1) THEN
          IPQTOQ(1,0) = 1
          DO I = 2,4
            JMROI(I) = I
            IPQTOQ(1,ISYMAX(5-I,2)) = I
          ENDDO
        ELSE
          IPQTOQ(1,   0) = 1
          IPQTOQ(1,IXYZ) = 1
          DO I = 2,4
            JMROI(I) = I
            IQDEF(I) = I
            IPQTOQ(1,ISYMAX(5-I,2)) = I
            IPQTOQ(1,ISYMAX(5-I,1)) = I
          ENDDO
        ENDIF
      ENDIF
C
C     Assign IQTOPQ
C
      DO IREP = 0,MAXREP
        DO IPQ = 1,NZ
          IQ = IPQTOQ(IPQ,IREP)
          IQTOPQ(IQ,IREP) = IPQ
        ENDDO
      ENDDO
C
C     Assign JQBAS - gives quaternion vector for given
C     component and bosonirrep
C     ================================================
C
      DO 60 IFR = 1,NFSYM
        DO 70 IC = 1,2
          DO 80 IM = 1,4
            JQBAS(JSPINR(IM,IC,IFR),IC) = JMROI(IM)
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
C
C     Give signs of quaternion phases: 1,i,-j,k
C
      IQPH(1,1) =  1
      IQPH(2,1) =  1
CTROND      IQPH(3,1) = -1
      IQPH(3,1) =  1
      IQPH(4,1) =  1
C
C     ...and their Hermitian conjugates: 1,-i, j,-k
C
      IQPH(1,2) =  1
      IQPH(2,2) = -1
CTROND      IQPH(3,2) =  1
      IQPH(3,2) = -1
      IQPH(4,2) = -1
C
C     ********************************************
C     *****   Define 4-component operators   *****
C     ********************************************
C
C     M = I_4,A_Z,A_Y,A_X, G_5,S_Z,S_Y,S_X
C
C
C     I_4       -->  I_2
C     ==================
C
      JM4REP(0) =  0
      JM4POS(0) =  1
      JM4FAS(0) =  1
C
C     i(alpha_z)   --> i(sigma_x)
C     ===========================
C
      JM4REP(1) =  ISYMAX(3,1)
      JM4POS(1) =  IZQ
      JM4FAS(1) =  1
C
C     i(alpha_y)   --> j(sigma_x)
C     ===========================
C
      JM4REP(2) =  ISYMAX(2,1)
      JM4POS(2) =  IYQ
      JM4FAS(2) =  1
C
C     i(alpha_x)   --> k(sigma_x)
C     ===========================
C
      JM4REP(3) =  ISYMAX(1,1)
      JM4POS(3) =  IXQ
      JM4FAS(3) =  1
C
C     | 0 I |       | 0 1 |
C     | I 0 |  -->  | 1 0 |
C     =====================
C
      JM4REP(4) =  IXYZ
      JM4POS(4) =  1
      JM4FAS(4) =  1
C
C     i(sigma_z)   --> i(I_2)
C     =======================
C
      JM4REP(5) =  ISYMAX(3,2)
      JM4POS(5) =  IZQ
      JM4FAS(5) =  1
C
C     i(sigma_y)   --> j(I_2)
C     =======================
C
      JM4REP(6) =  ISYMAX(2,2)
      JM4POS(6) =  IYQ
      JM4FAS(6) =  1
C
C     i(sigma_x)   --> k(I_2)
C     =======================
C
      JM4REP(7) =  ISYMAX(1,2)
      JM4POS(7) =  IXQ
      JM4FAS(7) =  1
C
C     ***************************************
C     ***** Definition of full operator *****
C     ***************************************
C
C
C     1. P             * scalar operator
C     ==================================
C
      MCMP(1)   = 1
      JM4 (1,1) = 0
      JCOM(1,1) = 1
C
C     2. i[alpha_x]P    * x-component of alpha times scalar operator
C     =============================================================
C
      MCMP(2)   =  1
      JM4 (1,2) =  3
      JCOM(1,2) =  1
C
C     3. i[alpha_y]P    * y-component of alpha times scalar operator
C     =============================================================
C
      MCMP(3)   =  1
      JM4 (1,3) =  2
      JCOM(1,3) =  1
C
C     4. i[alpha_z]P    * z-component of alpha times scalar operator
C     =============================================================
C
      MCMP(4)   =  1
      JM4 (1,4) =  1
      JCOM(1,4) =  1
C
C     5. i[Alpha x P]_x * vector product of alpha and vector operator,
C                         x-component
C     ================================================================
C
      MCMP(5)   =  2
      JM4 (1,5) =  2
      JM4 (2,5) =  1
      JCOM(1,5) =  1
      JCOM(2,5) = -1
C
C     6. i[Alpha x P]_y * vector product of alpha and vector operator,
C                         y-component
C     ================================================================
C
      MCMP(6)   =  2
      JM4 (1,6) =  1
      JM4 (2,6) =  3
      JCOM(1,6) =  1
      JCOM(2,6) = -1
C
C     7. i[Alpha x P]_z * vector product of alpha and vector operator,
C                         z-component
C     ================================================================
C
      MCMP(7)   =  2
      JM4 (1,7) =  3
      JM4 (2,7) =  2
      JCOM(1,7) =  1
      JCOM(2,7) = -1
C
C     8. iA.P           * dot-product of alpha and vector operator
C     ===========================================================
C
      MCMP(8)   =  3
      JM4 (1,8) =  3
      JM4 (2,8) =  2
      JM4 (3,8) =  1
      JCOM(1,8) =  1
      JCOM(2,8) =  1
      JCOM(3,8) =  1
C
C     9. gamma5 P       * gamma5 times scalar operator
C     ===========================================================
C
      MCMP(9)   =  1
      JM4 (1,9) =  4
      JCOM(1,9) =  1
C
C     10. i[Sigma_x]P   * x-component of sigma times scalar operator
C     ==============================================================
C
      MCMP(10)  =  1
      JM4 (1,10)=  7
      JCOM(1,10)=  1
C
C     11. i[Sigma_y]P   * y-component of sigma times scalar operator
C     ==============================================================
C
      MCMP(11)  =  1
      JM4 (1,11)=  6
      JCOM(1,11)=  1
C
C     12. i[Sigma_z]P   * z-component of sigma times scalar operator
C     ==============================================================
C
      MCMP(12)  =  1
      JM4 (1,12)=  5
      JCOM(1,12)=  1
C
C     13. i[betaSig_x]P * x-component of beta sigma times scalar operator
C     ==================================================================
C
      MCMP(13)  =  1
      JM4 (1,13)=  7
      JCOM(1,13)=  1
C
C     14. i[betaSig_y]P * y-component of beta sigma times scalar operator
C     ==================================================================
C
      MCMP(14)  =  1
      JM4 (1,14)=  6
      JCOM(1,14)=  1
C
C     15. i[betaSig_z]P * z-component of beta sigma times scalar operator
C     ==================================================================
C
      MCMP(15)  =  1
      JM4 (1,15)=  5
      JCOM(1,15)=  1
C
C     16. i[betaalp_x]P * x-component of beta alpha times scalar operator
C     ==================================================================
C
      MCMP(16)  =  1
      JM4 (1,16)=  3
      JCOM(1,16)=  1
C
C     17. i[betaalp_y]P * y-component of beta alpha times scalar operator
C     ==================================================================
C
      MCMP(17)  =  1
      JM4 (1,17)=  2
      JCOM(1,17)=  1
C
C     18. i[betaalp_z]P * z-component of beta alpha times scalar operator
C     ===================================================================
C
      MCMP(18)  =  1
      JM4 (1,18)=  1
      JCOM(1,18)=  1
C
C     Matrix symmetry
C     ===============
C
C     The matrix of an operator that is symmetric or
C     antisymmetric under time reversal has the following
C     structure:
C
C       A   B          A^{dagger} =   hA
C     -tB* tA*         B^T        = -thB
C
C     The matrix symmetry of matrices A and B
C     can be summarized as follows:
C
C     IH   ITIM   AR  AI  BR  BI
C      1     1     1   2   2   2
C     -1     1     2   1   1   1
C      1    -1     1   2   1   1
C     -1    -1     2   1   2   2
C
C     Only ITIM=1 is used (ITIM=-1 is transferred to IH=-1
C     by extracting an imaginary i), and the matrix symmetry
C     of matrices A and B is stored in IHQMAT(4,IH):
C
C     H-
      IHQMAT(1,-1) = 2
      IHQMAT(2,-1) = 1
      IHQMAT(3,-1) = 1
      IHQMAT(4,-1) = 1
C     H0
      IHQMAT(1, 0) = 0
      IHQMAT(2, 0) = 0
      IHQMAT(3, 0) = 0
      IHQMAT(4, 0) = 0
C     H+
      IHQMAT(1, 1) = 1
      IHQMAT(2, 1) = 2
      IHQMAT(3, 1) = 2
      IHQMAT(4, 1) = 2
C
C     IRQMAT gives the irrep of a given component of
C     quaternion matrix:
C
      DO IZ = 1,4
        DO IBRP = 0,MAXREP
          IRQMAT(IZ,IBRP) = IBTXOR(IFDIRP(IZ,1),IBRP)
        ENDDO
      ENDDO
C
C     Output section
C     ==============
C
      IF(IPREAD.GE.1) THEN
        CALL TITLER('Output from DBLGRP','*',103)
        IF    (NFSYM.EQ.1) THEN
          WRITE(LUPRI,'(3X,A,2X,A3)')
     &      '* One fermion irrep: ',FREP(1)
        ELSEIF(NFSYM.EQ.2) THEN
          WRITE(LUPRI,'(3X,A,2(2X,A3))')
     &      '* Two fermion irreps:',FREP(1),FREP(2)
        ENDIF
        IF    (NZ.EQ.4) THEN
          WRITE(LUPRI,'(3X,A)') '* Quaternionic group. NZ = 4'
        ELSEIF(NZ.EQ.2) THEN
          WRITE(LUPRI,'(3X,A)') '* Complex group. NZ = 2'
        ELSEIF(NZ.EQ.1) THEN
          WRITE(LUPRI,'(3X,A)') '* Real group. NZ = 1'
        ENDIF
        WRITE(LUPRI,'(3X,A)') '* Direct product decomposition:'
        IND = 1
        DO 110 I1 = 1,NFSYM
          DO 120 I2 = I1,NFSYM
            IND = IND + 1
            J = MOD(IND,2)+1
            WRITE(LUPRI,'(10X,2(A3,A),3(A3,A),A3)')
     &        FREP(I2),' x ',FREP(I1),' : ',
     &        REP(IFDIRP(1,J)),' + ',
     &        REP(IFDIRP(2,J)),' + ',
     &        REP(IFDIRP(3,J)),' + ',
     &        REP(IFDIRP(4,J))
  120     CONTINUE
  110   CONTINUE
        CALL HEADER('Spinor structure',-1)
        WRITE(LUPRI,'(/2(3X,A,I2,9X)/)')
     &     ('* Fermion irrep no.:',I,I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('La','|',REP(JSPINR(1,1,I)),
     &            '(',JQBAS(JSPINR(1,1,I),1),')',
     &               REP(JSPINR(2,1,I)),
     &            '(',JQBAS(JSPINR(2,1,I),1),')',
     &            '|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('Sa','|',REP(JSPINR(1,2,I)),
     &            '(',JQBAS(JSPINR(1,2,I),2),')',
     &               REP(JSPINR(2,2,I)),
     &            '(',JQBAS(JSPINR(2,2,I),2),')',
     &           '|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
C        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3),2X,A1,10X))')
     &     ('Lb','|',REP(JSPINR(3,1,I)),
     &            '(',JQBAS(JSPINR(3,1,I),1),')',
     &               REP(JSPINR(4,1,I)),
     &            '(',JQBAS(JSPINR(4,1,I),1),')',
     &           '|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &   ('Sb','|',REP(JSPINR(3,2,I)),
     &            '(',JQBAS(JSPINR(3,2,I),2),')',
     &             REP(JSPINR(4,2,I)),
     &            '(',JQBAS(JSPINR(4,2,I),2),')',
     &         '|',I=1,NFSYM)
        CALL HEADER('Quaternion symmetries',-1)
        WRITE(LUPRI,'(4X,A3,2X,A4)') 'Rep','T(+)'
        CALL PRSYMB(LUPRI,'-',29,4)
        DO IBRP = 0,NBSYM-1
          WRITE(LUPRI,'(4X,A3,4(2X,A1))')
     &       REP(IBRP),(QUNIT(IPQTOQ(IZ,IBRP)),IZ=1,NZ)
        ENDDO
      ENDIF
      RETURN
C
      END
C  /* Deck zmat */
      SUBROUTINE ZMAT(NONT,KATOM)
C
C  Modified for incorporation in HERMIT/SIRIUS/ABACUS, K.Ruud, May-95
C  Based on program ZMAT from CCQC/ Nov. 1988 (H.J.Aa.Jensen)
C**********************************************************
C*** LAST UPDATED ON SEPTEMBER 11, 1986 BY R. REMINGTON ***
C* REASON: TO OBTAIN HIGHER ACCURACY BY INCREASING THE  ***
C* NUMBER OF DIGITS FOR PI AND CHANGE ATAN TO DATAN #143 **
C**********************************************************
C***LAST UPDATED ON FEBRUARY 01, 1985 BY YUKIO YAMAGUCHI***
C**********************************************************
C   THIS PROGRAM CALCULATES CARTESIAN COORDINATES OF A SYSTEM
C   SPECIFIED BY INTERNAL COORDINATES.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      LOGICAL FIRST
#include "nuclei.h"
#include "molinp.h"
      COMMON/CBZMAT/IZ(MXCENT,4),BL(MXCENT),ALP(MXCENT),BET(MXCENT)
      DIMENSION NONT(KATOM)
      DATA FIRST /.TRUE./
      SAVE FIRST
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_NOFREE)
      DUMMY = 1.0D20
      IDUMMY = - 9 999 999
#endif
      NUCIND = NUCIND + 1
      IF (NUCIND .GT. MXCENT) GOTO 5000
      IZ(NUCIND,1) = 0
      IZ(NUCIND,2) = 0
      IZ(NUCIND,3) = 0
      IZ(NUCIND,4) = 0
      BL(NUCIND)=D0
      ALP(NUCIND)=D0
      BET(NUCIND)=D0
      READ (MLINE(NMLINE),'(A4)') NAMN(NUCIND)
      IF (FIRST) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
         ISTART = 5
         CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
         CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
#else
         READ (MLINE(NMLINE)(5:80),*) INUM, CHARGE(NUCIND)
#endif
         FIRST = .FALSE.
         GOTO 204
      ELSE
         IF(NUCIND .EQ. 2) GO TO 201
         IF(NUCIND .EQ. 3) GO TO 202
         GO TO 203
      END IF
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D)
 201  ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
      GO TO 204
 202  ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,2),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALP(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
      GO TO 204
 203  ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,2),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALP(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,3),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BET(NUCIND),'REA')
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,4),DUMMY,'INT')
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA')
#else
  201 READ (MLINE(NMLINE)(5:80),*)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),CHARGE(NUCIND)
      GO TO 204
  202 READ (MLINE(NMLINE)(5:80),*)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),IZ(NUCIND,2),
     &     ALP(NUCIND),CHARGE(NUCIND)
      GO TO 204
  203 READ (MLINE(NMLINE)(5:80),*)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),IZ(NUCIND,2),
     &     ALP(NUCIND),IZ(NUCIND,3),BET(NUCIND),
     &     IZ(NUCIND,4),CHARGE(NUCIND)
#endif
 204  CONTINUE
      NMLINE = NMLINE + 1
      NCLINE(NUCIND) = NMLINE
      NONT(NUCIND) = 1
      NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
      NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
      NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
C
      RETURN
C
C     Error message:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A/A,I5)')
     &    ' Too many atomic centers: MXCENT exceed in ZMAT,',
     &    ' Current limit:',MXCENT
        CALL QUIT('*** ERROR *** MXCENT exceeded in ZMAT')
      END
C  /* Deck buildz */
      SUBROUTINE BUILDZ(IPRINT,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,BASREF,
     &                NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,
     &                KPRIM,DOOWN)
C
C     Modified for HERMIT/SIRIUS/ABACUS -> K.Ruud, May-95
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      PARAMETER (DSM = 1.0D-30)
      PARAMETER (AR=PI/180.0D+00)
      PARAMETER (CORMAX = 1.D5)
      LOGICAL NOORBTS, DOOWN
      COMMON/CBZMAT/IZ(MXCENT,4),BL(MXCENT),ALP(MXCENT),BET(MXCENT)
      DIMENSION A(MXCENT),B(MXCENT),CZ(MXCENT,3),D(MXCENT)
      DIMENSION U1(3),U2(3),U3(3),U4(3),VJ(3),VP(3),V3(3)
      DIMENSION NONT(KATOM),IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NBLCK(KATOM,2),NUC(KBLOCK,2),NRC(KBLOCK,2),
     &          SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),CPRIM(KPRIM,KPRIM,KBLOCK,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2), NONCH(MXCENT)
      CHARACTER SPDCAR*1, NAVNK*4, CRT*1
      CHARACTER*80 BASREF(10,KATOM)
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "ccom.h"
#include "molinp.h"
#include "cbisol.h"
    1 FORMAT(/,2X,' INTERNAL COORDINATES'/)
    2 FORMAT(2X,I3,67X,F15.5)
    3 FORMAT(2X,I3,I5,F14.7,48X,F15.5)
    4 FORMAT(2X,I3,I5,F14.7,I5,F14.7,29X,F15.5)
    5 FORMAT(2X,I3,I5,F14.7,I5,F14.7,I5,F14.7,I5,F15.5)
C    6 FORMAT(//,25X,' CARTESIAN COORDINATES in Angstrom',
C     $      /,14X,'X',24X,'Y',24X,'Z'/)
    7 FORMAT(2X,I3,3(F20.10,5X))
    8 FORMAT(//,25X,' CARTESIAN COORDINATES in bohr',
     $      /,14X,'X',24X,'Y',24X,'Z'/)
    9 format(3f20.10)
C
      IF (SOLVNT) NUCIND = NUCIND - 1
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,1)
         WRITE(LUPRI,2) 1,CHARGE(1)
         IF (NUCIND .LE. 1) GOTO 201
         WRITE(LUPRI,3) 2,IZ(2,1),BL(2),CHARGE(2)
         IF(NUCIND.LE.2) GO TO 201
         WRITE(LUPRI,4) 3,IZ(3,1),BL(3),IZ(3,2),ALP(3),
     &                  CHARGE(3)
         IF(NUCIND.LE.3) GO TO 201
         DO 101 I=4,NUCIND
            WRITE(LUPRI,5) I,IZ(I,1),BL(I),IZ(I,2),ALP(I),
     &                     IZ(I,3),BET(I),IZ(I,4),CHARGE(I)
 101     CONTINUE
 201  CONTINUE
      END IF
      DO 102 I = 1,NUCIND
         DO 102 J = 1,3
            CZ(I,J) = D0
 102  CONTINUE
      DO 20 I = 1,NUCIND
         ALP(I) = ALP(I)*AR
         BET(I) = BET(I)*AR
 20   CONTINUE
      CZ(2,3) = BL(2)
      IF ((NUCIND-3) .GE. 0) THEN
         CZ(3,1) = BL(3)*DSIN(ALP(3))
         IF ((IZ(3,1)-1) .EQ. 0) THEN
            CZ(3,3) = BL(3)*DCOS(ALP(3))
         ELSE
            CZ(3,3) = CZ(2,3)-BL(3)*DCOS(ALP(3))
         END IF
         DO 80 I = 4,NUCIND
            IF((DABS(CZ(I-1,1))-1.0D-6) .LT. D0) THEN
               CZ(I,1) = BL(I)*DSIN(ALP(I))
               ITEMP = IZ(I,1)
               JTEMP = IZ(I,2)
               CZ(I,3) = CZ(ITEMP,3) - BL(I)*DCOS(ALP(I))*
     &                   DSIGN(D1,CZ(ITEMP,3) - CZ(JTEMP,3))
            ELSE
               GOTO 90
            END IF
 80      CONTINUE
 90      CONTINUE
         K=I
         IF ((K-NUCIND) .LE. 0) THEN
            DO 250 J = K,NUCIND
               CAJ = DCOS(ALP(J))
               SAJ = DSIN(ALP(J))
               CBJ = DCOS(BET(J))
               SBJ = DSIN(BET(J))
               IF(IZ(J,4) .EQ. 0) THEN
                  CALL VEC(U1,CZ,IZ(J,2),IZ(J,3))
                  CALL VEC(U2,CZ,IZ(J,1),IZ(J,2))
                  CALL VPROD(VP,U1,U2)
                  R = DSQRT(D1 - (U1(1)*U2(1) + U1(2)*U2(2)
     &                 + U1(3)*U2(3))**2)
                  DO 120 I = 1,3
                     U3(I) = VP(I)/R
 120              CONTINUE
                  CALL VPROD(U4,U3,U2)
                  DO 130 I = 1,3
                     VJ(I) = BL(J)*(-U2(I)*CAJ + U4(I)*SAJ*CBJ
     &                    + U3(I)*SAJ*SBJ)
                     ITEMP = IZ(J,1)
                     CZ(J,I)=VJ(I)+CZ(ITEMP,I)
 130              CONTINUE
               ELSE
                  IF (IABS(IZ(J,4)-1) .EQ. 0) THEN
                     CALL VEC(U1,CZ,IZ(J,1),IZ(J,3))
                     CALL VEC(U2,CZ,IZ(J,2),IZ(J,1))
                     ZETA = -(U1(1)*U2(1) + U1(2)*U2(2) +U1(3)*U2(3))
                     A(J) = (-CBJ + ZETA*CAJ)/(D1 - ZETA*ZETA)
                     B(J) = (CAJ - ZETA*CBJ)/(D1 - ZETA*ZETA)
                     R=D0
                     GAMMA=PI/D2
                     IF (ZETA .LT. D0) THEN
                        R=PI
                        GAMMA = DATAN(DSQRT(D1-ZETA*ZETA)/ZETA)+R
                     ELSE IF (ZETA .GT. D0) THEN
                        GAMMA = DATAN(DSQRT(D1-ZETA*ZETA)/ZETA)+R
                     END IF
                     D(J) = D0
                     IF ((DABS(GAMMA+ALP(J)+BET(J)-D2*PI)-1.0D-6)
     &                    .GE. D0) THEN
                        D(J) = IZ(J,4)*(DSQRT(D1+A(J)*CBJ-B(J)*CAJ))
     &                       /DSQRT(D1-ZETA*ZETA)
                     END IF
                     CALL VPROD(V3,U1,U2)
                     DO 200 I = 1,3
                        U3(I) = A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
                        VJ(I)=BL(J)*U3(I)
                        ITEMP=IZ(J,1)
                        CZ(J,I)=VJ(I)+CZ(ITEMP,I)
 200                 CONTINUE
                  ELSE
                     CALL VEC(U1,CZ,IZ(J,1),IZ(J,3))
                     CALL VEC(U2,CZ,IZ(J,2),IZ(J,1))
                     ZETA = -(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
                     CALL VPROD(V3,U1,U2)
                     V3MAG = DSQRT(V3(1)*V3(1)+V3(2)*V3(2)
     &                    + V3(3)*V3(3))
                     A(J) = V3MAG*CBJ/(D1-ZETA*ZETA)
                     B(J) = DSQRT((D1-CAJ*CAJ-A(J)*CBJ*V3MAG)
     &                    /(D1-ZETA*ZETA))
                     IF((IZ(J,4)-2) .NE. 0) THEN
                        B(J)=-B(J)
                     END IF
                     D(J) = B(J)*ZETA+CAJ
                     DO 240 I = 1,3
                        U3(I) = B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
                        VJ(I) = BL(J)*U3(I)
                        ITEMP = IZ(J,1)
                        CZ(J,I) = VJ(I)+CZ(ITEMP,I)
 240                 CONTINUE
                  END IF
               END IF
 250        CONTINUE
         END IF
      END IF
C
C   ELIMINATE IMAGINARY ATOM(S)
C   Hmmm, what do we do here? Remove them or keep them as dummy atoms?
C   Check to see if there is a basis set associated with the center?
C   Is it at all possible to place a dummy atom on top of an already
C   existing center in internal coordinates? K.Ruud, May-95
C
C      NATOMS=0
C      DO 290 I=1,NUCIND
C         IF(ANZ(I) .NE. D0) THEN
C            NATOMS=NATOMS+1
C            CHRG(NATOMS) = DFLOAT( ANZ(I) )
C            DO 280 J=1,3
C               C(NATOMS,J)=CZ(I,J)
C 280        CONTINUE
C         END IF
C 290  CONTINUE
C
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,8)
         DO 305 I=1,NUCIND
            WRITE(LUPRI,7) I,(CZ(I,J)/XTANG,J=1,3)
 305     CONTINUE
      END IF
C
C     Move information to HERMIT common blocks
C
      DO 410 I=1,NUCIND
         DO 420 J = 1, 3
            CORD(J,I) = CZ(I,J)/XTANG
            IF (ABS(CORD(J,I)).GT. CORMAX) GOTO 5010
 420     CONTINUE
 410  CONTINUE
C
C     Update MOLECULE.INP. In order to work properly with the symmetry
C     detection routines we need to sort atoms with similar charge
C     together in a block.
C
      ITYP = 1
      DO 600 I = 1, NUCIND
         DO 601 J = 1, I - 1
            IF (CHARGE(I) .EQ. CHARGE(J)) THEN
               NONTYP = NONTYP - 1
               DO 603 L = 1, ITYP
                  IF (CHARGE(I) .EQ. NONCH(L))
     &                   NONT(L) = NONT(L) + 1
 603           CONTINUE
               CHARGK = CHARGE(I)
               XCORD  = CORD(1,I)
               YCORD  = CORD(2,I)
               ZCORD  = CORD(3,I)
               NAVNK  = NAMN(I)
               DO 602 K = I - 1, J + 1, -1
                  CHARGE(K + 1) = CHARGE(K)
                  CORD(1,K + 1) = CORD(1,K)
                  CORD(2,K + 1) = CORD(2,K)
                  CORD(3,K + 1) = CORD(3,K)
                  NAMN(K + 1)   = NAMN(K)
 602           CONTINUE
               CHARGE(J + 1) = CHARGK
               CORD(1,J + 1) = XCORD
               CORD(2,J + 1) = YCORD
               CORD(3,J + 1) = ZCORD
               NAMN(J + 1)   = NAVNK
               GOTO 600
            END IF
 601     CONTINUE
         NONT(ITYP) = 1
         NONCH(ITYP) = CHARGE(I)
         ITYP = ITYP + 1
 600  CONTINUE
      IF (DOCART) CRT = 'C'
      IF (DOOWN)  CRT = 'X'
      WRITE(MLINE(5),'(A1,I4,I3,I2,10X,D10.2)') CRT,NONTYP,KCHARG,
     &                                          NSYMOP,THRS
      NMLINE = 6
      NHTYP = 0
      IBLOCK = 1
      IATOM = 1
      DO 605 I = 1, NONTYP
         NOORBTS = .TRUE.
         WRITE(MLINE(NMLINE),'(6X,F5.0,I4)') CHARGE(IATOM),NONT(I)
         NMLINE = NMLINE + 1
         DO 606 J = 1, NONT(I)
            NCLINE(IATOM) = NMLINE
            WRITE(MLINE(NMLINE)(1:4),'(A4)') NAMN(IATOM)
            WRITE(MLINE(NMLINE)(5:24),'(F20.16)')  CORD(1,IATOM)
            WRITE(MLINE(NMLINE)(25:44),'(F20.16)') CORD(2,IATOM)
            WRITE(MLINE(NMLINE)(45:64),'(F20.16)') CORD(3,IATOM)
            WRITE(MLINE(NMLINE)(65:80),'(A16)') '                '
            NMLINE = NMLINE + 1
            IATOM = IATOM + 1
 606     CONTINUE
C
C     Add basis set to newly order atomic set
C
         KAOVEC = KBLOCK + 1 - IBLOCK
         CALL BASLIB(IQM(I,1),JCO(1,I,1),NUC(IBLOCK,1),
     &               NRC(IBLOCK,1),SEG(IBLOCK,1),
     &               ALPHA(1,IBLOCK,1),
     &               CPRIM(1,1,IBLOCK,1),CPRIMU(1,1,IBLOCK,1),
     &               NBLOCK,ISGEN,KAOVEC,KPRIM,CHARGE(IATOM - 1),J,DSM,
     &               BASREF(1,I))
         NBLCK(I,1) = NBLOCK
         IBLOCK = IBLOCK + NBLCK(I,1)
         NOORBTS = NOORBTS .AND. IQM(I,1).EQ.0
         IF (IQM(I,1) .GT. 0) THEN
            NHTYP = MAX(NHTYP,IQM(I,1))
            IF (NHTYP .GT. MXQN) GOTO 5000
         END IF
         DO 610 K = 1,NONT(I)
            NOORBT(K + IATOM - 2) = NOORBTS
 610     CONTINUE
 605  CONTINUE
      NMLINE = NMLINE - 1
      IF (SOLVNT) NUCIND = NUCIND + 1
      RETURN
 5000 CONTINUE
         WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')
     &      '*  Input specifies highest orbital of atomic type ',
     &      I,' AS "',SPDCAR(NHTYP - 1),'".',
     &      ' Highest allowed orbital in this version: ',
     &      SPDCAR(MXQN - 1),
     &      ' Increase MXQN from',MXQN,' to',NHTYP,' and recompile.'
         CALL QUIT('Too high angular quantum no. specified in input.')
 5010 CONTINUE
        WRITE (LUPRI,'(A,E12.6,A/A/A,E12.6)')
     &    ' Atomic coordinate ',CORD(J,NUCIND),
     &    ' too large in BUILDZ',
     &    ' Note: Program is unstable for large coordinates.',
     &    ' Maximum coordinate value:',CORMAX
        CALL QUIT('*** ERROR: Atomic coordinate too large in BUILDZ')
      END
C  /* Deck vec */
      SUBROUTINE VEC(U,C,J,K)
#include "implicit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      DIMENSION C(MXCENT,3),R(3),U(3)
C
      R2=D0
      DO 101 I=1,3
         R(I)=C(J,I)-C(K,I)
         R2=R2+R(I)*R(I)
 101  CONTINUE
      R2=DSQRT(R2)
      DO 102 I=1,3
         U(I)=R(I)/R2
  102 CONTINUE
      RETURN
      END
C  /* Deck vprod */
      SUBROUTINE VPROD(VP,X,Y)
#include "implicit.h"
      DIMENSION VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bfginp */
      SUBROUTINE BFGINP(IAT,IBLCK,ICP,ALPIN,CPRIN,
     &              IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &              NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,FILNAM,BASREF)
C***********************************************************************
C
C      Read contracted basis set from a MOLFDIR type basis set file
C      (defined in FILNAM). IC tells us which component to read.
C      ICP = 1 : Read large component
C      ICP = 2 : Read small component
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "infpar.h"
#include "molinp.h"
C
      PARAMETER (DSM=1.0D-30)
C
      LOGICAL   SEG,TOBE,SEGIJ
      CHARACTER FILNAM*20,TEXT*80,DATE1*10,TIME1*8,FORMATE*10,FORMATC*10
      CHARACTER STRING*200, SZ(2)*5
      DIMENSION IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),NBLCK(KATOM,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2),ISGEN(KBLOCK),
     &          ALPIN(KPRIM),CPRIN(KPRIM*KPRIM)
      CHARACTER*80 BASREF(10,KATOM)
#include "dcbgen.h"
#include "cbirea.h"
#include "ccom.h"
#include "gencon.h"
C
      DATA SZ / 'LARGE', 'SMALL' /
C
#include "ibtfun.h"
C
      LUBFIL = 55
      NTYP   = 0
      NBLOCK = IBLCK
C
C     We now deal with general contracted basis sets:
C     Use Unrestricted kinetic balance.
C     ===============================================
C
C     URKBAL = .TRUE.
C     This doesn't work because it is reset in HAMINP
C
C     Look if basis set file exists and, if so, read in data
C     ======================================================
C
      IF (.NOT. PARHER) THEN
C        ... this is the master, slaves will read from lines
C        read into MLINE(*) by master
C        
C        The environment variable is a : separated string, e.g.:
C        /programs/Dirac/basis:/programs/Dalton/basis:.
C        
C        We loop over these directories in search for the basis set and
c        element given.
C        
         ISTART = 1
C        
 10      CONTINUE
C        
C        Find next colon in BASDIR
C        
         ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ':' )
C        
C        In there are no colons, we are the last directory in BASDIR
C        
         IF ( ICOLON .EQ. ISTART - 1 ) THEN
C           
C           Find end of string
C           
            ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ' ' )
C           
C           If end of string = start of string there are no more
C           directories, so we gracefully die...
C           
            IF ( ICOLON .EQ. ISTART )
     &           GOTO 5025
         END IF
         IEND = ICOLON - 1
         STRING = BASDIR(ISTART:IEND)
         ISTART = ICOLON + 1
C        
C        Add trailing slash if not present.
C        
         IEND = INDEX ( STRING(1:), ' ') - 1
         IF ( STRING(IEND:IEND) .NE. '/' ) THEN
            IEND = IEND + 1
            STRING(IEND:IEND) = '/'
         END IF
         STRING = STRING(1:IEND)//FILNAM
         IEND = INDEX(STRING(1:),' ') - 1
CLV : Reactivate when debugging, otherwise see below.
CLV      WRITE(LUPRI,'(3A/,A1,A,A1)')
CLV  &        'BFGINP: Trying basis set file for ',
CLV  &        SZ(ICP),' component: ',
CLV  &        '"',STRING(1:IEND),'"'
C        
C        Inquire if "/path/basis-set" exists.
C        
         INQUIRE (FILE = STRING, EXIST = TOBE)
         IF (.NOT.TOBE) GOTO 10
         OPEN (LUBFIL,FILE=STRING,FORM='FORMATTED',STATUS='OLD')
         REWIND (LUBFIL)
      ENDIF
C
C     Read info line of basis set file, copy this to BASREF array
C     ===========================================================
C
      CALL RDLINE(LUBFIL)
      READ (MLINE(NMLINE),'(A80)',END=5020) TEXT
C
      IF (ICP.EQ.1) THEN
C
C        We report the name of the file and the text line as both may give useful information
C
         BASREF(1,IAT) = "Basis set read from MOLFDIR file : "//FILNAM
         BASREF(2,IAT) = TEXT
      ELSE
C
C        Usually the small component comes from the same source as the LC, if not add additional info
C
         IF (BASREF(2,IAT).NE.TEXT) THEN
            BASREF(3,IAT) = 
     &   " Small component basis set read from MOLFDIR file : "//FILNAM
            BASREF(4,IAT) = TEXT
         ENDIF
      ENDIF
C
C     Read technical info line of basis set file
C     ==========================================
C
      CALL RDLINE(LUBFIL)
      READ (MLINE(NMLINE),'(A10,A8,2(A10,I4),F6.1,F20.5)',END=5020) 
     &      DATE1,TIME1,FORMATE,NE,FORMATC,NC,DCHRG,DEXPO
C
C     Read component, l-type and number of functions:
C     ICOMP  = 1 (large) or 2 (small)
C     ITYP   = l + 1
C     NPRIM  = # of primitives
C     NCONTR = # of contracted functions
C     NUNC   = # number of uncontracted functions
C     ===============================================
C
  301 CONTINUE
      CALL RDLINE(LUBFIL)
      READ (MLINE(NMLINE),'(6I4)',END=390,ERR=390)
     &      ICOMP,ITYP,NPRIM,NCONTR,IDUM,NUNC
      IF (ICOMP.EQ.0) GOTO 390 
C
C     Skip large component lines which we don't need in the small
C     component case
C     ===========================================================
C
      IF (ICOMP.NE.ICP) THEN
         IF (ICP.EQ.2) THEN
            NSTEP = NPRIM/NE + (NPRIM * NCONTR)/NC
            IF (MOD(NPRIM,NE).NE.0) NSTEP = NSTEP + 1
            IF (MOD(NCONTR*NPRIM,NC).NE.0) NSTEP = NSTEP + 1
            DO 310 INIT = 1, NSTEP
C              These lines should NOT be read in with RDLINE
C              because they need not be communicated to the slaves.
               IF (.NOT.PARHER) READ (LUBFIL,*,END=5021)
  310       CONTINUE
            GOTO 301
         ELSE
            GOTO 390
         ENDIF
      ENDIF
      IQM(IAT,ICP) = ITYP
C
C     Check dimensions of arrays
C     ==========================
C
      IF (NPRIM.GT.KPRIM) GOTO 5040
      IF (NCONTR.GT.KPRIM) GOTO 5060
C
C     Read in primitives and contraction coefficients
C
C     Number of exponents and contraction coefficients on a line is 4
C     ===============================================================
C
      DO I = 1, NPRIM, 4
         CALL RDLINE(LUBFIL)
         READ (MLINE(NMLINE),FORMATE,ERR=5020,END=5020)
     &         (ALPIN(IND),IND=I,MIN(I+3,NPRIM))
      ENDDO
      DO I = 1, NPRIM*NCONTR, 4
         CALL RDLINE(LUBFIL)
         READ (MLINE(NMLINE),FORMATC,ERR=5020,END=5020)
     &         (CPRIN(IND),IND=I,MIN(I+3,NPRIM*NCONTR))
      ENDDO
C
C     We cannot split the primitives in blocks. However, we can
C     split up the contracted functions making that the leading
C     dimension. The number of contracted functions in a block
C     depends on the # of angular components and is set to a max.
C     of MXCON_PARAM (gencon.h)
C     ==========================================================
C
      MXDUM = MXCON_PARAM - MOD(MXCON_PARAM,(NPRIM*KHK(ITYP)))
      MXCONTR = MAX((MXDUM/(NPRIM*KHK(ITYP))),1)
      JCO(ITYP,IAT,ICP) = 1 + (NCONTR-1)/MXCONTR
C
C     Put the entire list of ITYP primitives and contracted functions
C     into KPRIM long blocks
C     ===============================================================
C
      ICONS = 1
      ICONE = 0
      DO 501 JBLCK = 1, JCO(ITYP,IAT,ICP)
         IF (JBLCK*MXCONTR.LE.NCONTR) THEN
            ICONE = ICONE + MXCONTR
         ELSE
            ICONE = NCONTR
         ENDIF
         IF(NBLOCK.GT.KBLOCK) GOTO 5000
         CALL DZERO(ALPHA (1,  NBLOCK,ICP),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK,ICP),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK,ICP),KPRIM*KPRIM)
         IPRI = 0
         DO 502 JPRI = 1, NPRIM
            IPRI = IPRI + 1
            ALPHA(IPRI,NBLOCK,ICP) = ALPIN(JPRI)
            ICON = 0
            DO 503 JCON = ICONS, ICONE
               ICON = ICON + 1
               CPRIMU(IPRI,ICON,NBLOCK,ICP) = CPRIN(NPRIM*(JCON-1)+JPRI)
  503       CONTINUE
  502    CONTINUE
         NUCIJ = NPRIM
         NRCIJ = ICONE - ICONS + 1
         SEGIJ = NRCIJ .GE. 0
C
C        Remove exponents with all zero coefficients from block
C        This happens in SC in case of different primitives from
C        l+1 and l-1. It also happens with uncontracted functions
C
  600    CONTINUE
         DO 601 JPRI = 1, NUCIJ
            NZ = 0
            DO 602 JCON = 1, NRCIJ
               IF (ABS(CPRIMU(JPRI,JCON,NBLOCK,ICP)).GT.DSM) NZ = NZ + 1
  602       CONTINUE
            IF (NZ.EQ.0) THEN
               DO 603 IPRI = JPRI+1, NUCIJ
                  ALPHA(IPRI-1,NBLOCK,ICP) = ALPHA(IPRI,NBLOCK,ICP)
                  DO 604 JCON = 1, NRCIJ
                     CPRIMU(IPRI-1,JCON,NBLOCK,ICP) =
     &                      CPRIMU(IPRI,JCON,NBLOCK,ICP)
  604             CONTINUE
  603          CONTINUE
               ALPHA(NUCIJ,NBLOCK,ICP) = 0.0D0
               DO 605 JCON = 1, NRCIJ
                  CPRIMU(NUCIJ,JCON,NBLOCK,ICP) = 0.0D0
  605          CONTINUE
               NUCIJ = NUCIJ - 1
               GOTO 600
            ENDIF
  601    CONTINUE
C
C        Identify segmented contractions
C
         CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK,ICP),KPRIM,DSM)
         SEG(NBLOCK,ICP) = SEGIJ
C
C        Reorder primitive orbitals
C
         CALL PRIORD(ALPHA(1,NBLOCK,ICP),CPRIMU(1,1,NBLOCK,ICP),NUCIJ,
     &               NRCIJ,SEG(NBLOCK,ICP),KPRIM,DSM)
C
C        Normalize orbitals
C
         CALL NRMORB(ITYP,NRCIJ,NUCIJ,ALPHA(1,NBLOCK,ICP),
     &               CPRIM(1,1,NBLOCK,ICP),
     &               CPRIMU(1,1,NBLOCK,ICP),KPRIM,NBLOCK)
C
         NUC(NBLOCK,ICP) = NUCIJ
         NRC(NBLOCK,ICP) = NRCIJ
C
         ISGEN (NBLOCK) = 0
C
         IF (INPTST)
     &      WRITE(LUPRI,'(2A,6I5)') SZ(ICP),' MOLFBAS Block:',
     &      NBLOCK,NUC(NBLOCK,ICP),NRC(NBLOCK,ICP)
         NBLOCK = NBLOCK + 1
         ICONS = ICONE + 1
  501 CONTINUE
      GOTO 301
  390 CONTINUE
      NBLCK(IAT,ICP) = NBLOCK - IBLCK
  391 CONTINUE
      CLOSE (LUBFIL)
      RETURN
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> BFGINP error, number of AO-blocks            ',NBLOCK,
     *  '                   max. number (KBLOCK in READIN) ',KBLOCK,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5020 CONTINUE
        WRITE (LUPRI,'(/A/A/2A)')
     &     'Error reading the following line :',MLINE(NMLINE),
     &     'in MOLFDIR type basis set file : ',FILNAM
        CALL QUIT('Error reading contents of MOLFDIR basis set file.')
 5021 CONTINUE
        WRITE (LUPRI,'(/A/A/2A)')
     &     'Error reading the line following this line :',MLINE(NMLINE),
     &     'in MOLFDIR type basis set file : ',FILNAM
        CALL QUIT('Error reading contents of MOLFDIR basis set file.')
 5025 CONTINUE
        WRITE (LUPRI,'(/2A)')
     &     'Cannot find MOLFDIR type basis set file : ',FILNAM
        CALL QUIT('Error finding MOLFDIR type basis set file.')
 5040 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BFGINP error, number of primitives           ',NPRIM,
     *  '                   current maximum number (MAXPRI)',KPRIM
        CALL QUIT('Too many primitives, increase MAXPRI in input.')
 5060 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BFGINP error, number of contracted functions ',NCONTR,
     *  '                   current maximum number (MAXPRI)',KPRIM
        CALL QUIT('Too many contracted fun., increase MAXPRI in input.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck kinbal */
      SUBROUTINE KINBAL(I,IBLL,IBLS,
     &     IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &     NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,INPTST,
     &     CNTBAS)
C*****************************************************************************
C
C      Generate small component functions from large using the
C      kinetic balance relation
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "ccom.h"
C
      LOGICAL SEG,INPTST,CNTBAS
      DIMENSION IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),NBLCK(KATOM,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2),ISGEN(KBLOCK)
#include "ibtfun.h"
      IDOWN(J) = IBTAND(1,J)
      IUP(J)   = IBTAND(1,IBTSHR(J,1))
C
      IBS = IBLS - 1
      IBU = IBLL + JCO(1,I,1) - 1
      IBD = IBLL - 1
      JBS = IBS
      IQM(I,2) = 0
      DO J = 1,(IQM(I,1)+1)
        IBJ = IBS
C
C       Generate function from downwards
C
        IF(J.GT.1) THEN
          DO K = 1,JCO(J-1,I,1)
            IBD = IBD + 1
            IF(IDOWN(ISGEN(IBD)).EQ.1
     $           .OR. ISGEN(IBD) .EQ. 0 )  THEN
              IF (CNTBAS) THEN
                MXDUM = 50 - MOD(50, KHK(J))
                MAXS = MAX((MXDUM/KHK(J)),1)
                IBJCO=1+(NUC(IBD,1)-1)/MAXS
                JSEND = 0
                DO JBLCK=1, IBJCO
                  JSBEG = JSEND
                  IBS = IBS + 1
                  IF (JBLCK*MAXS.LE.NUC(IBD,1)) THEN
                    JSEND = JSEND + MAXS
                  ELSE
                    JSEND = NUC(IBD,1)
                  ENDIF
                  CALL GENSML(J,IBD,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                      SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
                ENDDO
              ELSE
                JSBEG = 0
                JSEND = NUC(IBD,1)
                IBS = IBS + 1
                CALL GENSML(J,IBD,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                      SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C       Generate function from upwards
C
        IF(J.LT.IQM(I,1)) THEN
          DO K = 1,JCO(J+1,I,1)
            IBU = IBU + 1
            IF(IUP(ISGEN(IBU)).EQ.1
     $           .OR. ISGEN(IBU) .EQ. 0 ) THEN
              IF (CNTBAS) THEN
                MXDUM = 50 - MOD(50, KHK(J))
                MAXS = MAX((MXDUM/KHK(J)),1)
                IBJCO=1+(NUC(IBU,1)-1)/MAXS
                JSEND = 0
                DO JBLCK=1, IBJCO
                  JSBEG = JSEND
                  IBS = IBS + 1
                  IF (JBLCK*MAXS.LE.NUC(IBU,1)) THEN
                    JSEND = JSEND + MAXS
                  ELSE
                    JSEND = NUC(IBU,1)
                  ENDIF
                  CALL GENSML(J,IBU,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                      SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
                ENDDO
              ELSE
                JSBEG = 0
                JSEND = NUC(IBU,1)
                IBS = IBS + 1
                CALL GENSML(J,IBU,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                      SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        JCO(J,I,2) = IBS - IBJ
        IF(JCO(J,I,2).GT.0) IQM(I,2) = J
      ENDDO
      NBLCK(I,2) = IBS - JBS
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gensml */
      SUBROUTINE GENSML(IQM,IBL,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                  SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
C*****************************************************************************
C
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "infpar.h"
C
      PARAMETER ( ALPTHR = 1.0D-3 )
C
      LOGICAL SEG,INPTST
      DIMENSION NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),CPRIMU(KPRIM,KPRIM,KBLOCK,2)
C
C
C     Small component basis set is always uncontracted,
C     hence NRC(IBS) = NUC(IBL) and SEG(IBS) = T
C
      SEG(IBS,2)   = .TRUE.
C
cccC     CALL DCOPY(NUC(IBS,2),ALPHA(1,IBL,1),1,ALPHA(1,IBS,2),1)
C
C     Remove small component basis functions already present.
C
      NUC_IBS = 0
      DO I = JSBEG + 1,JSEND
         ALPTEST = ALPHA(I,IBL,1)
C
         DO J = IBJ+1,IBS-1
            DO IJ = 1,NUC(J,2)
              TEST = ABS(ALPTEST - ALPHA(IJ,J,2))
              IF (TEST .LE. ALPTHR) GO TO 100
            END DO
         END DO
C
         NUC_IBS = NUC_IBS + 1
         ALPHA(NUC_IBS,IBS,2) = ALPTEST
C
  100 CONTINUE
      END DO
      IF ((JSEND-JSBEG) .NE. NUC_IBS) THEN
         IF (.NOT.PARHER) WRITE(LUPRI,'(A,I3,A,I3)')
     &   ' INFO: GENSML deleted',(JSEND-JSBEG) - NUC_IBS,
     &   ' small functions in a block out of',(JSEND-JSBEG)
      END IF
      IF (NUC_IBS .EQ. 0) THEN
C
C        We've deleted all functions from this block.
C        This might happen for (dual) family basis sets.
C        Decrement IBS and return.
C
         IBS = IBS - 1
         RETURN
      END IF
      NUC(IBS,2) = NUC_IBS
      NRC(IBS,2) = NUC_IBS
C
C     Contraction matrix for SMALL is just a unit matrix
C
      CALL DUNIT(CPRIMU(1,1,IBS,2),KPRIM)
      CALL NRMORB(IQM,NRC(IBS,2),NUC(IBS,2),ALPHA(1,IBS,2),
     &            CPRIM(1,1,IBS,2),CPRIMU(1,1,IBS,2),KPRIM,IBS)
C
C     Test print
C
      IF (INPTST) THEN
         WRITE(LUPRI,'(A)') 'Contraction matrix:'
         DO 400 L = 1, NUC(IBS,2)
            WRITE (LUPRI,'(E16.10)') ALPHA(L,IBS,2)
            WRITE (LUPRI,'(6(2X,F10.8))')
     $           (CPRIMU(L,M,IBS,2),M=1,NUC(IBS,2))
  400    CONTINUE
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      BLOCK DATA BDPGROUP
C
C Initialize /PGROUP/
C
#include "implicit.h"
#include "pgroup.h"
C
      DATA SYMOP/ ' E ','Oyz','Oxz','C2z',
     &            'Oxy','C2y','C2x',' i '/
      DATA GROUPS/ 'C1 ','C2 ','Cs ','Ci ',
     &             'D2 ','C2v','C2h','D2h'/
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rkbrel */
      SUBROUTINE RKBREL(IPRINT)
C******************************************************************************
C
C     Routine for Restricted Kinetic Balance (RKB) Relations.
C
C     Generates array IRKB running over shells
C       For large shells IRKB points to L+1 small partner
C       For small shells IRKB points to L-2 small partner
C       or to the L-1 large partner if the L-2 small partner does
C       not exist (pointer is negative in this case).
C
C     Written by Jon K. Laerdahl 05.12.1996
C     Major surgery by T.Saue May 5 1998
C     Changed definition of IRKB because of problems with contracted sets,
C     L. Visscher, Nov 16, 2000.
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "dcbbas.h"
#include "aovec.h"
      PARAMETER (SML = 1.0D-6)
#include "primit.h"
#include "dcbrkb.h"
C
      IF (IPRINT.GE.3) WRITE(LUPRI,'(/1X,A)')
     &   'RKBREL: Constructing relations for scalar RKB transformation'
      CALL IZERO(IRKB,MXPRIM)
C
C     New code for generating RKBREL.
C     This should also work for contracted large component basis set.
C
C     Loop over large and small component shells.
C     For contracted basis sets, the same primitives are present in many
C     shells, and we skip shells with the same primitives.
C     A potential problem arises if you block a contracted basis set.
C
C
      iPrimTest = 0
C
      DO 100 ISHL = 1, KMAX
C
         NCENTI = NCENT(ISHL)
         NHKTI  = NHKT(ISHL)
C
         iPrimStart = JSTRT(ISHL) + 1
         iPrimEnd   = JSTRT(ISHL) + NUCO(ISHL)
C
C        Skip primitives we've already processed.
C
         IF (iPrimStart .EQ. iPrimTest) GO TO 100
         iPrimTest = iPrimStart
C
C
C        Loop over primitives in this shell.
C
         DO 200 IPRIM = iPrimStart, iPrimEnd
            PRIEXI = PRIEXP(IPRIM)
C
            JSH2M = 0
            jPrimTest = 0
C
C           Loop over small component shells
C           Find matching small component basis functions.
C
            DO 300 JSHL = NLRGSH+1, KMAX
               NCENTJ = NCENT(JSHL)
               NHKTJ  = NHKT(JSHL)
C
C              Must be on same atomic center.
C
               IF ( NCENTJ .NE. NCENTI ) GO TO 300
C
               jPrimStart = JSTRT(JSHL) + 1
               jPrimEnd   = JSTRT(JSHL) + NUCO(JSHL)
C
C              Skip primitives we've already processed.
C
               IF ( jPrimStart .EQ. jPrimTest) GO TO 300
               jPrimTest = jPrimStart
C
C              Loop over primitives.
C
               DO 400 JPRIM = jPrimStart, jPrimEnd
C
                  PRIEXJ = PRIEXP(JPRIM)
C
                  IF ( ABS ((PRIEXJ-PRIEXI)/PRIEXI) .LT. SML ) THEN
C
C
C                    We've got a match: same center, same exponent
C
C
C                    Match small component to L-2 small function....
C
                     IF(NHKTJ.EQ.NHKTI-2.AND.ISHL.GT.NLRGSH) THEN
C
C                       Check for two identical functions.
C                       This should never happen, since GENSML deletes
C                       functions it already knows.
C
                        IF(JSH2M.NE.0) CALL QUIT
     &                   ('RKBREL: Error, two identical S.C. functions')
C
C                       Note that we use JSH2M and not JPRIM, since RKBLOW
C                       in dirone.F wants shell number, and not primitive
C                       number. This is not the same if L Component basis
C                       set is contracted.
C
                        JSH2M = JSHL + (JPRIM - jPrimStart)
                        IRKB(IPRIM) = JSH2M
                     ENDIF
C
C                    Match large component to L+1 small function and
C                    vice versa. The initialization of the SC is
C                    overwritten if an (L-1) SC does also exist. It
C                    serves to check whether this SC is used to balance
C                    by itself.
C
                     IF (NHKTJ.EQ.NHKTI+1.AND.ISHL.LE.NLRGSH) THEN
                        JSH1P = JSHL + (JPRIM - jPrimStart)
                        ISH1M = ISHL + (IPRIM - iPrimStart)
                        IRKB(IPRIM) = JSH1P
                        IRKB(JPRIM) = - ISH1M
                     ENDIF
                  ENDIF
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('RKBREL: RKB array',-1)
         DO I = 1,NPSHEL
            WRITE(6,'(A,I3,A,I3)') 'IRKB(',I,') = ',IRKB(I)
         END DO
      END IF
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Getpot */
      FUNCTION GETPOT(Q)
C***********************************************************************
C
C     Get nuclear potential energy
C     Surgery by T.Saue June 7 1997
C     Point charges are not included. /June2000,jth+hjaaj
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0, DP1 = 0.1D0, DP5 = 0.5D0)
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
      DIMENSION Q(*)
C
#include "ibtfun.h"
C
      GETPOT = D0
      DO 10 N = 1, NUCIND
         DO 20 M = N, NUCIND
            DO 30 KB = 0, MAXREP
            IF (IBTAND(KB,ISTBNU(M)) .EQ. 0) THEN
               IF (M.EQ.N .AND. KB.EQ.0) GO TO 30
               DS = (CORD(1,N)-CORD(1,M)*PT(IBTAND(ISYMAX(1,1),KB)))**2
     &            + (CORD(2,N)-CORD(2,M)*PT(IBTAND(ISYMAX(2,1),KB)))**2
     &            + (CORD(3,N)-CORD(3,M)*PT(IBTAND(ISYMAX(3,1),KB)))**2
               DIST = SQRT(DS)
               IF (DIST.LT.DP1 .AND. Q(M)*Q(N).NE.D0 .AND.
     &             .NOT.(NOORBT(M).OR.NOORBT(N)) ) GOTO 5000
C                  ... exit if close and no floating and no point charge
               IF (IBTAND(KB,ISTBNU(N)) .EQ. 0) THEN
                  IF (.NOT.(NOORBT(M).OR.NOORBT(N))) THEN
C                 if (no point charge) then
                     HKAB  = FMULT(IBTAND(ISTBNU(M),ISTBNU(N)))
                     IF (M .EQ. N) HKAB = DP5*HKAB
                     GETPOT = GETPOT + Q(M)*Q(N)*HKAB/DIST
                  END IF
               END IF
            END IF
   30      CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
C
 5000 CONTINUE
      WRITE (LUPRI,'(3(A,I3),A,D15.5)')
     &     ' The nucleus',N,' is too close to the',KB,
     &     'th transformation of nucleus no.',M,' DIST=',DIST
      CALL QUIT('Nuclei too close.')
      END
C/*  Deck frefrm */
      SUBROUTINE FREFRM(STRING,ISTART,IVAL,RVAL,TYPE)
C
C     Simulate free format input for Cray computers when reading from an
C     internal file, K.Ruud-Feb.97
C
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER*80 STRING
      CHARACTER TYPE*3, FRMT*7
C
      IF (TYPE .EQ. 'INT') THEN
         IVAL = 0
      ELSE
         RVAL = 0.0D0
      END IF
      DO 10 IPOS = ISTART, 80
         IF (STRING(IPOS:IPOS) .NE. ' ') GO TO 20
 10   CONTINUE
C
 20   CONTINUE
      IEND = INDEX(STRING(IPOS:),' ')
      ILEN = IEND - ISTART - 1 + IPOS
      IF (TYPE .EQ. 'INT') THEN
         IF (ILEN .GT. 9) THEN
            WRITE (FRMT,'(2X,A2,I2,A1)') '(I',ILEN,')'
         ELSE
            WRITE (FRMT,'(3X,A2,I1,A1)') '(I',ILEN,')'
         END IF
         READ (STRING(ISTART:80),FRMT) IVAL
      ELSE IF (TYPE .EQ. 'REA') THEN
         IF (ILEN .GT. 9) THEN
            WRITE (FRMT,'(A2,I2,A3)') '(F',ILEN,'.0)'
         ELSE
            WRITE (FRMT,'(1X,A2,I1,A3)') '(F',ILEN,'.0)'
         END IF
         READ (STRING(ISTART:80),FRMT) RVAL
      ELSE
         WRITE (LUPRI,'(/,A)') 'Illegal variable type in FREFRM'
         CALL QUIT('Illegal variable type in FREFRM')
      END IF
      ISTART = IEND + IPOS
      RETURN
      END
C  /* Deck genfamily */
      SUBROUTINE GENFAMEXP(LUINFO,IFLAG)
C***********************************************************************
C
C     Read in or calculate exponents for family basis sets
C
C     Written by J. Thyssen - Jun 4 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER ( D0 = 0.00 D00 , D1 = 1.00 D00 )
      CHARACTER KEYWRD*6
C
#include "cbirea.h"
#include "molinp.h"
#include "sdpre.h"
C
C
C
      IF ( IFLAG .EQ. 0 ) THEN
C
C        Generate well-temperered exponents
C        Geometric progression is just a special case
C        with gamma = 0 and delta = 0.
C
C        See Huzinaga and Klobukowski, CPL 212, 260 (1993).
C
         IF (NFAMEXP(1) .LT. 1 .OR. NFAMEXP(1) .GT. MXFAMEXP) GOTO 5000
         IF (FAMPAR(1) .LE. D0) GOTO 5010
         IF (FAMPAR(2) .LE. D1) GOTO 5010
         IF (FAMPAR(3) .LT. D0) GOTO 5010
         IF (FAMPAR(4) .LT. D0) GOTO 5010
C
         IF ( INPTST ) THEN
            WRITE(LUPRI,'(A,/,4(A,F20.10,/),A,I3)')
     $           'Parameters for well-tempered progression: ',
     $           'alpha = ',FAMPAR(1),
     $           'beta  = ',FAMPAR(2),
     $           'gamma = ',FAMPAR(3),
     $           'delta = ',FAMPAR(4),
     $           'N     = ',NFAMEXP(1)
         END IF
C
C        Formula is:
C
C        zeta_N = alpha,
C
C        zeta_N-k+1 = zeta_N-k+2 beta ( 1 + gamma (k/N)^delta), k = 2,N
C
         FAMEXP(NFAMEXP(1), 1) = FAMPAR(1)
         DO K = 2,NFAMEXP(1)
            FAMEXP(NFAMEXP(1) - K + 1, 1) =
     $           FAMEXP(NFAMEXP(1) - K + 2, 1) * FAMPAR(2) *
     $           ( D1 + FAMPAR(3) *
     $           ( ( SDPRE(K) / SDPRE(NFAMEXP(1)) ) ** FAMPAR(4) ) )
         END DO
C
C        Copy into set 2.
C
         NFAMEXP(2) = NFAMEXP(1)
         DO K = 1, NFAMEXP(1)
            FAMEXP( K, 2 ) = FAMEXP( K, 1 )
         END DO
C
      ELSE IF ( IFLAG .EQ. 1 ) THEN
C
C        Read exponents from file
C
         DO I = 1,NFAMEXP(1)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*) FAMEXP(I, 1)
         END DO
C
C        Copy into set 2.
C
         NFAMEXP(2) = NFAMEXP(1)
         DO K = 1, NFAMEXP(1)
            FAMEXP( K, 2 ) = FAMEXP( K, 1 )
         END DO
C
      ELSE IF ( IFLAG .EQ. 2 ) THEN
C
C        Read dual set of exponents from file
C
         DO I = 1,NFAMEXP(1)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*) FAMEXP(I, 1)
         END DO
         DO I = 1,NFAMEXP(2)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*) FAMEXP(I, 2)
         END DO
C
      ELSE
         WRITE(LUPRI,'(/,A,/,A,I3)')
     $        '*** ERROR in GENFAMEXP ***:',
     &        'Unknown IFLAG: ',IFLAG
         CALL QUIT('Unknown IFLAG in GENFAMEXP')
      END IF
C
      IF (INPTST) THEN
         CALL HEADER('Exponents for family basis set(1)',-1)
         CALL OUTPUT(FAMEXP,1,NFAMEXP(1),1,1,NFAMEXP(1),1,1,LUPRI)
         IF ( IFLAG .EQ. 2 ) THEN
            CALL HEADER('Exponents for family basis set(2)',-1)
            CALL OUTPUT(FAMEXP(1,2),1,NFAMEXP(2),1,1,NFAMEXP(2),1,1,
     &           LUPRI)
         END IF
      END IF
      RETURN
C
 5000 CONTINUE
      WRITE(LUPRI,'(/,A,/,2(A,I3/))')
     &     '*** ERROR in GENFAMEXP *** ',
     $     'Illegal value for number of exponents: ',NFAMEXP(1),
     &     'Value must be between 1 and ',MXFAMEXP
      CALL QUIT('Illegal value for number of exponents')
 5010 CONTINUE
      WRITE(LUPRI,'(2(A,/),4(A,F20.10,/))')
     &     '*** ERROR in GENFAMEXP *** ',
     $     'Illegal non-positive value for parameters in'//
     &     ' well-tempered progression',
     &     'Alpha = ',FAMPAR(1),
     &     'Beta  = ',FAMPAR(2),
     &     'Gamma = ',FAMPAR(3),
     &     'Delta = ',FAMPAR(4)
      CALL QUIT('Illegal value for parameter in ' //
     $     'well-tempered progression')
      END
C  /* Deck PRBASREF */
      SUBROUTINE  PRBASREF (LUPRI,NONTYP,BASREF)
C***********************************************************************
C
C     Print out the reference information for basis sets
C     Sets that are used for more than one atom type are printed
C     only once.
C
C     Written by L. Visscher -  May 14, 2003
C
C***********************************************************************
#include "mxcent.h"
      DIMENSION I_SET(MXCENT),I_EQUAL(MXCENT)
      LOGICAL FIRST(MXCENT)
      CHARACTER*80 BASREF(10,*)
C
C     Check for common sets, the first set is always unique
C
      FIRST(1) = .TRUE.
      DO I = 2, NONTYP
         FIRST(I) = .TRUE.
         DO J = 1, I - 1
            IF (BASREF(1,I).EQ.BASREF(1,J)) THEN
               FIRST(I) = .FALSE.
               I_SET(I) = J
               GOTO 1
            ENDIF
         ENDDO
   1     CONTINUE
      ENDDO
C
      WRITE (LUPRI,1000)
      DO I = 1, NONTYP
C        We will print the information only once
         IF (FIRST(I)) THEN
C           Find out which other types also use the same basis set
            N = 1
            I_EQUAL(N) = I
            DO J = I+1, NONTYP
               IF (I_SET(J).EQ.I) THEN
                  N = N + 1
                  I_EQUAL(N) = J
               ENDIF
            ENDDO
C           print line with info about the equal sets
            WRITE (LUPRI,1001) (I_EQUAL(J),J=1,N)
C           print basis set information
            DO ILINE = 1, 10
               IF (BASREF(ILINE,I)(1:16).NE.' Not initialized')
     &         WRITE (LUPRI,1002) BASREF(ILINE,I)
            ENDDO
         ENDIF
      ENDDO
C
 1000 FORMAT (//"  References for the basis sets"/
     &          "  -----------------------------")
 1001 FORMAT (/"  Atom type",16I4,(/11X,16I4))
 1002 FORMAT (2X,A80)
      RETURN
      END
