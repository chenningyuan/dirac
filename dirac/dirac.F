C  /* Deck dirac */
      PROGRAM DIRAC
C*****************************************************************************
C     
C     Direct Iterative Relativistic All-electron Calculations
C     
C     Original version written by:
C     Trond Saue, University of Oslo, Norway, 1992
C
C     Further contributions by:
C
C     NOTE
C
C     This is an experimental code for the relativistic calculation of
C     atomic and molecular electronic wavefunctions.
C     The author accept no responsibility for the performance of the code
C     or for the correctness of the results.
C
C     The code (in whole or part) is not to be reproduced for
C     further distribution without the written permission of the author.
C
C     If results obtained with this code are published, an appropriate
C     citation would be:
C     T.Saue:"DIRAC(Direct Iterative Relativistic All-electron
C     Calculations), a program for relativistic molecular calculations".
C
C     This note is not to be removed or altered.
C
C*****************************************************************************
C
C     R E V I S I O N S:
C
C     3.0 - Property module:
C            - Expectation values
C
C     2.0 - September 1995: T.Saue
C           Symmetry implemented: D2H and subgroups + time reversal
C           Dynamic allocation implemented using environmental variable DIRWRK
C             (modelled after the ABACUS code)
C           Convergence acceleration DIIS implemented (J.Laerdahl & T.Saue)
C
C
C*****************************************************************************
C
C     PAM - Paul Adrien Maurice
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#if defined (SYS_IRIX)
      INTEGER SCHEDMODE,SGS_GANG
      PARAMETER(SCHEDMODE=13,SGS_GANG=4)
CTROND      INCLUDE 'limits.h'
CTROND      INCLUDE 'sys/types.h'
CTROND      INCLUDE 'sys/prctl.h'
CTROND      INCLUDE 'sys/schedctl.h'
#endif
#endif
#if defined (VAR_PVM)
#include <mxorb.h>
#include <mxcent.h>
#include <infpar.h>
#include <molinp.h>
#include <fpvm3h.h>
#endif
      LOGICAL DOLOGO
#if defined (SYS_T3D)
      PARAMETER (LWORK = 10 000 000)
#else
      PARAMETER (LWORK = INSTALL_WRKMEM)
#endif
#if defined (SYS_HPUX)
      PARAMETER (LUERR = 7)
#else
      PARAMETER (LUERR = 0)
#endif
C
#if !defined (SYS_CONVEX)
      PARAMETER (IWOFF = 0)
#endif
#if defined (SYS_AIX) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_CONVEX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_HPUX)
      CHARACTER*20 DIRWRK
#endif
#if defined (SYS_CRAY) || defined (SYS_AIX) || defined (SYS_IRIX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_HPUX)
      POINTER (WADDR,WORK(1))
#else
#if defined (SYS_CONVEX)
      DIMENSION WORK(2)
#else
C
C     Static memory allocation
C
      DIMENSION WORK(LWORK+2)
#endif
#if defined (SYS_IBM)
      COMMON /DIRWRK/ WORK
#endif
#endif
      CHARACTER CPUTID*12, SECTID*12
      DATA WRKDLM/8H*WRKDLM*/
      LMWORK = 0
C
      CALL QENTER('DIRAC')
C
C
#if defined (SYS_AIX) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_CONVEX) || defined (SYS_T3D) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_HPUX)
C
C     Dynamic allocation of memory;check for DIRWRK in environment
C     ============================================================
C
      DIRWRK = '                    '
      CALL GETENV('DIRWRK',DIRWRK)
      READ (DIRWRK,'(I20)',ERR=10) LMWORK
      GO TO 20
   10    WRITE (LUERR,*)
     &   'DIRWRK conversion error; DIRWRK = "',DIRWRK,'"'
         LMWORK = 0
   20 CONTINUE
#endif
C
C     Either use default memory size
C     
      IF (LMWORK .LE. 0) LMWORK = LWORK
      NBYTES = (LMWORK+2) * 8
      XMBYTES = NBYTES
      XMBYTES = XMBYTES / (1024*1024)
C     
 1000 FORMAT(/' >>> DIRAC: requested memory could not be allocated.',
     &     /'     Probable causes:',
     &     /'       (1) Insufficient free memory+page+swap space',
     &     /'       (2) "datasize" resource limit too small.')
C     
#if defined (SYS_AIX) || defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_HPUX)
#if defined (SYS_AIX)
      WADDR = MALLOC(%VAL(NBYTES))
#endif
#if defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_HPUX)
      WADDR = MALLOC(NBYTES)
#endif
      IF (WADDR .EQ. 0) THEN
         WRITE(LUPRI,1000)
         WRITE(LUERR,1000)
         CALL QUIT(' DIRAC: memory allocation failure')
      END IF
#endif
#if defined (SYS_CONVEX)
C
C     Convex dynamic memory allocation
C     ================================
C
      CALL NALLOC(NBYTES,IPTR,IERR)
      IF (IERR .LT. 0) THEN
         WRITE(LUERR,*) ' >>> Convex NALLOC failure; IERR = ', IERR
         IF (LUPRI .NE. 0)
     &        WRITE(0,*) ' >>> Convex NALLOC failure; IERR = ', IERR
         IF (IERR .EQ. -3) THEN
            WRITE(LUPRI,1000)
            IF (LUPRI .NE. 0) WRITE(0,1000)
         END IF
         STOP ' DIRAC: memory allocation failure'
      ENDIF
      IWOFFB = IPTR - LOC(WORK(1))
      IWOFF  = IWOFFB / 8
      IF (IWOFF*8 .NE. IWOFFB) THEN
         WRITE (0,'(A,3(/A,I12))')
     &        'Convex alignment problem, see DIRAC main program',
     &        '   pointer from NALLOC    :',IPTR,
     &        '   off-set to WORK, IWOFFB=',IWOFFB,
     &        '   mod(IWOFFB,8) =         ',IWOFFB - IWOFF*8
         stop 'DIRAC: Convex alignment failure, see main program'
      end if
#endif
#if defined (SYS_CRAY)
C     
C     Cray dynamic memory allocation
C     
C     HPALLOC IERR codes:
C     -1 : length allocated smaller than zero
C     -2 : no memory left for the system
C     -7 : memory corrputed.
C     HPALLOC IABORT .ne. 0: abort on error
      IABORT = 0
      CALL HPALLOC(WADDR,LMWORK+2,IERR,IABORT)
      IF (IERR .LT. 0) THEN
         WRITE (LUERR,'(/A,I3)')
     &        ' >>> Cray HPALLOC failure; IERR   =', IERR
         IF (LUPRI .NE. 0) WRITE (0,'(/A,I3)')
     &        ' >>> Cray HPALLOC failure; IERR   =', IERR
         IF (IERR .EQ. -2) THEN
            WRITE(LUPRI,1000)
            IF (LUPRI .NE. 0) WRITE(0,1000)
         END IF
         CALL QUIT('DIRAC core HPALLOC failure.')
      ENDIF
#endif
#if defined (SYS_IBM)
C     
C     IBM VM/CMS: 212 : disable end of record error for formatted I/O
C     (such that variable length input records may be read)
C     
      CALL ERRSET(212,256,-1,1)
#endif
C
C     Run DIRAC
C     =========
C     MPI start
#if defined (VAR_MPI)
C     
C     MPI initialization
C     
      CALL MPI_INIT(IERR)
C     
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYTID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMNOD,IERR)
#if defined (SYS_IRIX)
C
C     Impose gang scheduling to avoid problems with idle nodes
C     when high system load
C
      IF(MYTID.EQ.0.AND.SCHEDCTL(SCHEDMODE,SYS_GANG,0).LT.0) THEN
        WRITE(LUPRI,'(A)') '*** ERROR(DIRAC) gang scheduling '
        CALL QUIT('Error in gang scheduling !')   
      ENDIF
#endif
C     
      NUMNOD = NUMNOD - 1
      MPARID = 0
C
Ctec for parallel debugging purpose     
#if defined (VAR_DEBUG)
      call sleep(5)
#endif
C
      IF (MYTID .EQ. 0) THEN
C     
C     Print DIRAC heading
C     ===================
C     
         CALL GETTIM(CSTR,WSTR)
         DOLOGO = .TRUE.
         CALL LOGO(DOLOGO)
         CALL TIMDAT
         IF (LMWORK .LE. 0) THEN
            WRITE (LUPRI,*) 'DIRAC: default work memory size used.'
C     
C     ... or user specified memory size.
C     
         ELSE
            WRITE (LUPRI,'(/A/3A)')
     &           ' DIRAC: user specified work memory size used,',
     &           '        environment variable DIRWRK = ',DIRWRK
         END IF
         WRITE (LUPRI,'(A,I12,A,F8.2,A/)')
     &        ' Work memory size (LMWORK) :',LMWORK,' =',XMBYTES
     &        ,' megabytes.'

         CALL DIRCTL(WORK,LMWORK)
C     
C     Send ending signal to slaves
C     
         NFINIS = -1
         CALL MPI_BCAST(NFINIS,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
         CALL GETTIM(CEND,WEND)
         CTOT = CEND - CSTR
         WTOT = WEND - WSTR
         CPUTID = SECTID(CTOT)
         WRITE(LUPRI,'(/A,A12)') '>>>> Total CPU  time used in DIRAC:', 
     &                            CPUTID
         CPUTID = SECTID(WTOT)
         WRITE(LUPRI,'(/A,A12)') '>>>> Total wall time used in DIRAC:',
     &                            CPUTID
         WRITE(LUPRI,'(/A/)')
     +     '********** E N D   of   D I R A C  output  **********'
      ELSE
         CALL DIRNOD(WORK,LMWORK)
      END IF
#if defined (SYS_IRIX)
      CALL SLEEP(1)
#endif
      CALL MPI_FINALIZE(IERR)
C     MPI end
#else
C     Single and PVM start
c     OPEN(LUCMD, FILE = 'DALTON.INP')
c     CALL GNRLIN(WORK,LMWORK)
c     IF (TESTIN) THEN
c     CALL QUIT('*** End of input test for GENERAL ***')
c     END IF
c     CLOSE(LUCMD)
C     PVM start
#if defined (VAR_PVM)
C     
C     Spawn all PVM processes, as we now know how many there are
C     
      IF (PARCAL) THEN
         CALL PVMFMYTID(MYTID)
         IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C     
         CALL PVMFSPAWN("hernod.x",PVMDEFAULT,'*',NODES,MTIDS,NUMT)
         IF (NUMT .NE. NODES) CALL PFAIL(-40)
      END IF
C     PVM end
#endif
C     
C     Print DIRAC heading
C     ===================
C     
      CALL GETTIM(CSTR,WSTR)
      DOLOGO = .TRUE.
      CALL LOGO(DOLOGO)
      CALL TIMDAT
      IF (LMWORK .LE. 0) THEN
         WRITE (LUPRI,*) 'DIRAC: default work memory size used.'
C     
C     ... or user specified memory size.
C     
      ELSE
         WRITE (LUPRI,'(/A/3A)')
     &        ' DIRAC: user specified work memory size used,',
     &        '        environment variable DIRWRK = ',DIRWRK
      END IF
      WRITE (LUPRI,'(A,I12,A,F8.2,A/)')
     &     ' Work memory size (LMWORK) :',LMWORK,' =',XMBYTES
     &     ,' megabytes.'
      CALL DIRCTL(WORK,LMWORK)
C     PVM start
#if defined (VAR_PVM)
      IF (PARCAL) THEN
C     
C     Tell the slaves that the calculation is finished at last
C     --------------------------------------------------------
C     
         NFINIS = -1
         CALL PVMFPACK(INTEGER4,NFINIS,1,1,INFO)
         IF (INFO .LT. 0) CALL PFAIL(INFO)
C     
         CALL PVMFPACK(STRING,MLINE,80*NMLINE,1,INFO)
         IF (INFO .LT. 0) CALL PFAIL(INFO)
C     
C------------------------------
C     Pack common-block INFORB.
C------------------------------
C     
         CALL PVMFPACK(INTEGER4,MULD2H,338,1,INFO)
         IF (INFO .LT. 0) CALL PFAIL(INFO)
C     
C-----------------------------
C     Broadcast message buffer
C-----------------------------
C     
         CALL PVMFMCAST(NODES,MTIDS,10,INFO)
         IF (INFO .LT. 0) CALL PFAIL(INFO)
      END IF
C     PVM end
#endif
C     Single and PVM end
      CALL GETTIM(CEND,WEND)
      CTOT = CEND - CSTR
      WTOT = WEND - WSTR
      CPUTID = SECTID(CTOT)
      WRITE(LUPRI,'(/A,A12)') '>>>> Total CPU  time used in DIRAC:', 
     &                         CPUTID
      CPUTID = SECTID(WTOT)
      WRITE(LUPRI,'(/A,A12)') '>>>> Total wall time used in DIRAC:', 
     &                         CPUTID
      WRITE(LUPRI,'(/A/)')
     +  '********** E N D   of   D I R A C  output  **********'
#endif
      CALL QEXIT('DIRAC')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirctl */
      SUBROUTINE DIRCTL(WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
      DIMENSION WORK(LWORK)
C
      CALL QENTER('DIRCTL')
C
C***********************************************************************
C*****   I  N  P  U  T    S  E  C  T  I  O  N   ************************
C***********************************************************************
C
      CALL PAMINP(WORK,LWORK)
C
C*****************************************************************************
C*****  S E T U P  -  M O D U L E  *******************************************
C*****************************************************************************
C
      CALL PAMSET(WORK,LWORK)
C
C*****************************************************************************
C**** M A I N    S E C T I O N  **********************************************
C*****************************************************************************
C
      IF(.NOT.DOHRM) CALL PAMDRV(WORK,LWORK)
C
C******************************************************************************
C*****  T  E  R  M  I  N  A  T  I  O  N  **************************************
C******************************************************************************
C
      CALL QEXIT('DIRCTL')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamdrv */
      SUBROUTINE PAMDRV(WORK,LWORK)
C*****************************************************************************
C
C     Driver for DIRAC - main part
C
C     Written by T. Saue - November 1994
C     Last revision : Sep 11 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PAMDRV')
C
C******************************************************************************
C*****  D O   G E O M E T R Y   O P T I M I Z A T I O N  **********************
C******************************************************************************
C
      IF (MINIMI) CALL PAMGEO(WORK,LWORK)
C
C******************************************************************************
C*****  G E T    W A V E    F U N C T I O N   *********************************
C******************************************************************************
C
      IF(DOPIN) CALL PCMOIN(WORK,LWORK)
      IF(DOPSI) CALL PAMPSI(WORK,LWORK)
      IF(DOPUT) CALL PCMOUT(WORK,LWORK)
C
C*****************************************************************************
C*****  A N A L Y S I S    M O D U L E  **************************************
C*****************************************************************************
C
      IF(MINIMI) CALL MINANA(WORK,LWORK)
      IF(DOANA) CALL PAMANA(WORK,LWORK)
C
C
C*****************************************************************************
C*****  R E S P O N S E    M O D U L E  **************************************
C*****************************************************************************
C
      IF(MINIMI) CALL MINPRP(WORK,LWORK)
      IF(DOPRP) CALL PAMPRP(WORK,LWORK)
#if defined (MOD_TRA)
C
C*****************************************************************************
C***** I N D E X  T R A N S F O R M A T I O N*********************************
C*****************************************************************************
C
      IF(DOTRA) CALL PAMTRA(WORK,LWORK)
#endif
#if defined (MOD_CC)
C
C*****************************************************************************
C***** R E L C C S D  M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCCM) CALL PAMCCM(WORK,LWORK)
#endif
#if defined (MOD_CI)
C
C*****************************************************************************
C***** D I R R C I   M O D U L E ********************************************
C*****************************************************************************
C
      IF(DOCIM) CALL PAMCIM(WORK,LWORK)
C
#endif
C
C*****************************************************************************
C***** H E I N Z *************************************************************
C*****************************************************************************
C
      IF(DOHEINZ) CALL HEINZ(WORK,LWORK)
C
 666  CONTINUE
      CALL QEXIT('PAMDRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck logo */
      SUBROUTINE LOGO(DOLOGO)
C*****************************************************************************
C
C     DIRAC - logo , generated from gif-picture by
C                    Jo Asplin , university of Tromsoe 1994
C
C     Last revised: Nov 23 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      LOGICAL DOLOGO
C
      WRITE(LUPRI,'(A)')
     +'************************************************************'
     +  //'******************',
     +'*                                                           '
     + //'                 *',
     +'*                                O U T P U T                '
     + //'                 *',
     +'*                                   from                    '
     + //'                 *',
     +'*                                                           '
     + //'                 *',
     +'*                   @@@@@    @@   @@@@@     @@@@     @@@@@  '
     + //'                 *',
     +'*                   @@  @@        @@  @@   @@  @@   @@      '
     + //'                 *',
     +'*                   @@  @@   @@   @@@@@    @@@@@@   @@      '
     + //'                 *',
     +'*                   @@  @@   @@   @@ @@    @@  @@   @@      '
     + //'                 *',
     +'*                   @@@@@    @@   @@  @@   @@  @@    @@@@@  '
     + //'                 *',
     +'*                                                           '
     + //'                 *',
     +'*                                                           '
     + //'                 *'
      IF(DOLOGO) THEN
        WRITE(LUPRI,'(A)')
     +'"o $ $ $ "o"o $ o" "o o o"o $ "oo $ " $ $ $ o"o"o"oo $ $ " $'
     +   //'o"o" o" "o"o $ $ "',
     +' $ "o $ $ o $ o"o""o "o"o" " $ o $ "" $ oo "o"o o o"o $ $ "o'
     +   //'$ o"oo" $ o $ "o "',
     +'" "o "o"o"o "o $ o"o"oo o$" """" "          o   "   o   o"o"'
     +   //'$ o"o "o "o"o"o" $',
     +'"o"o"o o "o"o""o"o "ooo                          " " o  " "o'
     +   //'o" $o"o"oo $ $ "o',
     +' "o " "o" $ o"o"o $ o  o       "             "   "    o     '
     +   //'o"oo o o o"o o"$ "',
     +'"o" $ "o"o " $ o" "o o o  " $            o       o "   " o  '
     +   //'  oo"  "o"o $ $ o"',
     +' $ "o"o oo""oo "o"" "o o  $"                  "      o o  ""'
     +   //'    "o""o o"o $ $',
     +'o"o"o $ o"o"  $ " $ $ " ""                 o    o  $  o  $o"'
     +   //'$   o "o"o"o " " $',
     +' "o" $ $ $ ""o"o""o "oo"            o  "      "  o o " "$ $$'
     +   //'oo    $oo o"o"$ "',
     +'"o $ o"o"o $ o"o" o $                   o  " o" o $ o"o"$$$o'
     +   //'$o$  oo"o"o "  ""',
     +'o"o o"o""o"o"o"o$ $$             "  "      $ " " o o"o$$$$$$'
     +   //'o$"$ o$$$o" $ $ "o',
     +'o"o" $ $ $ $"o$ $ $              o  "    "     "o "oo$$$$$$$'
     +   //'$o"$$$$$$$$oo"o" o',
     +' $ $ $ $o $ "o"o"$"          "                  oo"$"o$$$$$"'
     +   //'$$$$$$$$$$$$$o " o',
     +'"o"o" $ o"o"o"o"o"         o             o   " o $"$$$$$$$$$'
     +   //'$$$$$$$$$$$$$$ " o',
     +'"o $ "oo"o"o"o$ ""              "     "  o  "o$$$$"$o$$$$$$$'
     +   //'$$$$$$$$$$$$"o  "o',
     +'"o"o"o"o" " $ $o""$             o  "o $o$o$o$$$$$$$$$$$$$"$$'
     +   //'$$$$$$$$$$$$$ "o"',
     +'oo "o $ "o$ $"o"o o$  o"  o$$$$oo  $o$$$$$$$$$$$$$$$$$$o$$$$'
     +   //'$$$$$$$$$$$$$$o "',
     +'o $ $ "o"oo" $ "o$ "   $o$$$$$$$$   "$$$$$$$"$$$"$$$$$$$$"$$'
     +   //'$$$$"$$$"o"$$$o "o',
     +'o" $ $" $ o"o $ $ $ $ $$$$"$$$$$$"   $$$$$$"$oo$$$$$o$$o$$o$'
     +   //'$$o$$$"$o$$$$$$  "',
     +'o"o "o" "o $ $ $ $ $   " """$$$""   $$$"$$ "o"$o$o"$o$"oo"$"'
     +   //'o$$$$$$$$$o$$$$ "o'
        WRITE(LUPRI,'(A)')
     +'o " "o $ $ o"o $ $ $o     "o$""    o$$$$oo  "   "o"" o"oo$$$'
     +   //'$$$$$$"$$$$$$$$ "',
     +' $ "o $ $ "o $ o"o"o"$    "       "$$$$$$$         o"o$ $o$$'
     +   //'$$$$$$"$$$$$$$$ "o',
     +' $ $ "o "o"o"o" $ $ $             $$$$$$$$$        "o"oo"$$o'
     +   //'$$o$$$$"$o$$$$  "',
     +' "o $ """oo $ "oo"o""o         o $$$$$$o$$$$o    "o$"$o"$o$$'
     +   //'$$$$$"$$$$$$"o "o"',
     +'"o $ "$ $ $ $ "o $o"$ $   "   o$$$$$$$$$$$$$$  " o $"$"$o$$$'
     +   //'$$$o$$$$$$" o "o o',
     +' "o "o" "o"o"o$ $ $o"o"      """$$$$$$$$$$$$$ " o $ $$"$$$$$'
     +   //'$$$$$$$$o "o  " o',
     +'" $ "o"o"o "o"o"o"o "$o$   o o    "$$$"$$$$$$$" $o"$o$$$$$$$'
     +   //'$$$$$$$$o "  "o  "',
     +'o" "o "o" "o"o $o"o"$o $           o $"$$$$$$o$$o$$"$"$$$$$$'
     +   //'$$$o$$$$$ " " o" "',
     +' "o" "o $ " "o"oo"o"o"o"o     o  oo" $$$$$$$$$$$$$$" $"$$"$$'
     +   //'$$$$$$$$$o$o" o "o',
     +'"o "" $ o"o" "oo"oo"o""$"    oo$$$$$$$$$$$$$$$$$$"$ $ $$$$$$'
     +   //'$$"$$o$$$ "$$oo"',
     +'"o $ $ "o" o"o $o $ $ $""$    $"     oooo$$$$$$$"$$"o$$$$$$$'
     +   //'$$$$$$$$$  "$$$$o',
     +'"o o $ " " $ "o $ $ "$ $$ o  o    o$$"$$$"$$"$o$$$$$$$$$$$$$'
     +   //'$$"$$o$$ o  $$$$$$',
     +' o"o"o ""o"o"o"o $ ""o"o""$o o              $"$$$$$$$$$$$$$$'
     +   //'o$$$$$$    "o$$$$$',
     +'"o $ "o"o" "o $ $ "$ "o"o$ $o"o          o $"$$o$$$$$$$$$$$$'
     +   //'$$$$$$  "   $$$$$$',
     +'"o o"o " $ $oo "o""o""o"o"$ $o$  "   oo"o$$$$$$$$$$$$$$$$$$"'
     +   //'$$$$"   o " $$$$$$',
     +'o $ o $ $ "o "o" "o"o"o"o$ $o $o"$o$ $$$$$$$$$$$$$$$$$$$$$$$'
     +   //'$$$"  "   oo$$$$$$',
     +'"o "o "o "o $ $ $"ooo"o"o"oo$ $ $ $ $"$$$$$$$$$$$$$$$$$$$$$o'
     +   //'$$   "  "  $$$$$$$',
     +' $ "o"o"o"o"o $ "o $"o"o"o"o"$o"o"o$o$$ "$$$$$$$$$$$$$$$$o$$'
     +   //'"     o  o$$$$$$$$',
     +'o"o"o $  $ o" $ $ $ $ $ $o"o$o$o$$$$$$$ oo$ $"$$"$$$$$$$$$$"'
     +   //'o   "    o$$$$$$$$',
     +' "o"o $ "o"o"$ $ $ $"o"oo$$$$$$$$$$$$$$ $o" $ $$$$$$$"$$$"  '
     +   //' "o   " o$$$$$$$$$',
     +'"o"o $ " $o $ $oo"o"$o$$$$$$"$$o$$$$$$$  $$  "o$"$$o$$$"   "'
     +   //'  o   "o$$$$$$$$$$',
     +' $ $ o" "o $ $ oo"o"$$$$$$$$$$$$$$$$$$$   "$o "$$$$$$"     "'
     +   //'   "  o$$$$$$$$$$$',
     +' $ o"o"o"o"o"o"o $ $$$$$$o$$$"$$$o$$$$$$    "$$"$$$"        '
     +   //'   "oo$$$$$$$$$$$$'
        WRITE(LUPRI,'(A)')
     +'"o"o "o $ $ o"o"o"$$$$"$$$$$$$$$$"$$"$$$      ""$"          '
     +   //' o  o$$$$$$$$$$$$$',
     +'"o $"o"o"o $ $ $ $$$$$$$$$o$$$$o$$$$$$$$   o$$$             '
     +   //'   o$$$$$$$$$$$$$$',
     +'o"o $ "o $ o"o"o"$$$$"$$"$$$o$$$"$$$$$$$o o$$$$$o           '
     +   //' o$$$$$$$$$$$$$$$$',
     +' $ $ $o $o"o"o"o"$$o$$$$$$$$$"$$$$$"$o$$$ $$$$$$$$oo      " '
     +   //' o$$$$$$$$$$$$$$$$',
     +'"o"o"o"ooo"o"o"$o$$$$$$$"$$o$$$$$o$$$$$$$o$$$$$$$$$$$oo     '
     +   //'o$$$$$$$$$"$$$"$o$',
     +' $ oo $ oo"o"o"oo$$$$$o$$$$$"o$$$"$$$$"$$$$$$$$$$$$$$$"     '
     +   //'$$$$$$$$$$$$$$$$"$',
     +'"o"o "o"o $ "o"o$$o$$o$$$$o$"$$$$$$$$$$$$$$$$$$$$$$$$$    "o'
     +   //'$$$$$$o$"$"$o$$$$$',
     +'*                                                           '
     + //'                 *'
      ENDIF
      WRITE(LUPRI,'(A)')
     +'*          ================================================='
     + //'======           *',
     +'*                      Program for Atomic and Molecular '
     &//'                     *',
     +'*           Direct Iterative Relativistic Allelectron '
     &//'Calculations           *',
     +'*          ================================================='
     + //'======           *',
     +'*             Written by: Trond Saue, University of Oslo, No'
     + //'rway             *',
     +'*                                -- 1992 --                 '
     + //'                 *',
     +'*                                                           '
     + //'                 *',
     +'*                        VERSION 3.0     JUL-1996           '
     + //'                 *',
     +'*                                                           '
     + //'                 *',
     +'************************************************************'
     + //'******************'
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampsi */
      SUBROUTINE PAMPSI(WORK,LWORK)
C*****************************************************************************
C
C     Wave function module
C
C     Written by T.Saue April 1996
C     Last revision: April 27 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
#include <dcbgen.h>
#include <dcbpsi.h>
C
      CALL QENTER('PAMPSI')
C
      IF(DODHF) CALL PSIDHF(WORK,LWORK)
#if defined (MOD_TRA)
      IF(DOMP2) CALL PSIMP2(WORK,LWORK)
#endif
C      IF(DODFT) CALL PSIDFT(WORK,LWORK)
C
      CALL QEXIT('PAMPSI')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamcim */
      SUBROUTINE PAMCIM(WORK,LWORK)
C*****************************************************************************
C
C     Dummy routine while the interface is being implemented
C
      CALL QUIT('*** DIRRCI interface not implemented **')
C
      RETURN
      END
#if defined (VAR_MPI) || defined (VAR_PVM)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnod */
#if defined (VAR_MPI)
      SUBROUTINE DIRNOD(WORK,MWORK)
#endif
#if defined (VAR_PVM)
      PROGRAM DIRNOD
#endif
C*****************************************************************************
C
C     Menu file for slave nodes. Check whether the program needs a slave
C     and determine the type of work to be done.
C
C     Written by L.Visscher July 1997
C     Last revision: July 29 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
C
#if defined (VAR_PVM)
      PARAMETER (MWORK = 6 000 000)
#endif
#include  <mxorb.h>
#include <infpar.h>
C
      DIMENSION WORK(MWORK)
C
      CALL QENTER('DIRNOD')
      LWORK = MWORK
#include <memint.h>
C
#if defined (VAR_MPI)
 100  CONTINUE
#endif
#if defined (VAR_PVM)
      CALL PVMFMYTID(MYTID)
      IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C
      CALL PVMFPARENT(MPARID)
      IF (MPARID .LT. 0) CALL PFAIL(MPARID)
C
 100  CONTINUE
#endif
C
C     Start timing.
C
      CALL GETTIM(CPU1,WALL1)
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
#if defined (VAR_MPI)
      CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
      CALL PVMFUNPACK(INTEGER4,ITASK,1,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
#endif
      IF (ITASK .LT. 0) THEN
C
C        Global end of program : go up for MPI_FINALIZE
C
         GOTO 999
      ELSE IF (ITASK .EQ. 0) THEN
C
C        Nothing to be done at the moment. Go to sleep.
C
         CONTINUE
      ELSE IF (ITASK .EQ. 1) THEN
C
C        Start calculating integrals.
C
         CALL HERNOD(WORK,LWORK)
      ELSE IF (ITASK .EQ. 2) THEN
C
C        Start calculating integrals.
C
C         CALL TRANOD(WORK,LWORK)
      ELSE IF (ITASK .EQ. 3) THEN
C
C        Start coupled cluster calculation.
C
         CALL CCNOD(WORK,LWORK)
      ELSE
C
C        Unknown task. Abort execution.
C
         WRITE(LUPRI,'(A,I5)')
     &     '*** ERROR *** DIRNOD: Unknown task ITASK = ',ITASK
         CALL QUIT(' Unknown task in DIRNOD')
      END IF
C
C     We want the slave program to go until there is nothing
C     more to calculate, K.Ruud, May-95
C     But they may rest for a while to start refreshed for their new job.
C     L. Visscher. Jul-1997
C
      GOTO 100
C
#if defined (VAR_PVM)
      CALL PVMFEXIT(INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
#endif
C
 999  CALL QEXIT('DIRNOD')
      END
#endif
