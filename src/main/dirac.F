!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C FILE    : dirac.F
C
C***********************************************************************
      SUBROUTINE DIRAC(IPARCAL,IMYTID,IMPARID,INUMNOD)
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
C
#include "infpar.h"
#include "dgroup.h"
#ifdef VAR_OMP
      INTEGER omp_get_num_procs
#endif
C
C
C     Set up basic information for parallel execution, to be stored in the
C     common block infpar
C
#ifdef VAR_OMP
         write(lupri,'(//A,I5)')
     &   ' * openMP activated, number of threads:',
     &   omp_get_num_procs()
#endif
      PARCAL = (IPARCAL .NE. 0)
      MPARID = IMPARID
C We count the number of slave nodes, subtract the master
C     NUMNOD = INUMNOD - 1
C     SK - Sep 08: subtraction has already been done in mpixinit.
      NUMNOD = INUMNOD
      MYTID  = IMYTID
C Additional check for parallelism, parallel with only a master is the same as serial
      PARCAL = PARCAL.AND.NUMNOD.ne.0
C Test the communication by sending an empty task
      ITASK = 0
      IF (PARCAL) CALL DIRAC_PARCTL( ITASK )
C
C     Set basis set directory
C
      CALL GETBASDIR()
C
C     Initialize output units for gp/mempkg.F (MEMGET/MEMREL/MEMCH routines)
C
      CALL MEMINI(lupri,lupri)
C
C     Execute Dirac
C
      write(lupri,'(/a/a)') 'Execution time and host',
     &                       '-----------------------'
      CALL TSTAMP(' ',LUPRI)
      CALL DIRCTL()
C
C     Send finish signal to slaves
C
      IF (PARCAL) CALL DIRAC_PARCTL( -2 )
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirctl */
      SUBROUTINE DIRCTL()

      use interface_ao_specific
      use interface_mo_specific
      use xmlout
      use dirac_cfg
#ifdef HAS_PCMSOLVER
      use pcm_scf
#endif

      implicit none

#include "priunit.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"

C
      CALL QENTER('DIRCTL')
      call xml_open_file('dirac.xml')
      call xml_begin('run','type="dirac-run"')
C
C***********************************************************************
C*****   I  N  P  U  T    S  E  C  T  I  O  N   ************************
C***********************************************************************
C
      CALL PAMINP()
      IF ( INPTES ) RETURN
C
C***********************************************************************
C*****  S E T U P  -  M O D U L E  *************************************
C***********************************************************************
      CALL PAMSET()
!     write info to files
!     this is used in dft and visualization
      call interface_ao_write()
      call interface_mo_write()

#ifdef HAS_PCMSOLVER
      if (dirac_cfg_pcm) then
! Quit if symmetry or X2C are requested.
        if (nbsym > 1) then
          call quit('Polarizable Continuum Model calculation cannot '//
     & 'handle symmetry!')
        end if
        if (bss.or.x2c) then
          call quit('Polarizable Continuum Model calculation cannot '//
     & 'handle 2-component Hamiltonians!')
        end if
        call pcm_scf_initialize(lupri)
        write(lupri, '(//A/)') 'PCMSolver initialized.'
      end if
#endif
C
C
C***********************************************************************
C**** M A I N    S E C T I O N  ****************************************
C***********************************************************************
C
      IF(.NOT.DOHRM) THEN
         IF (OPTIMI) THEN
            CALL PAMGEO()
         ELSE
            CALL PAMDRV()
         ENDIF
      END IF
C
C***********************************************************************
C*****  T  E  R  M  I  N  A  T  I  O  N  *******************************
C***********************************************************************
C

#if defined (BUILD_GEN1INT) && defined (BUILD_OPENRSP)
!     clean up Gen1Int interface
      call gen1int_host_finalize()
#endif
      call xml_end('run')
      call xml_close_file

      CALL QEXIT('DIRCTL')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamdrv */
      SUBROUTINE PAMDRV()
C***********************************************************************
C
C     Driver for DIRAC - main part
C
C     Written by T. Saue - November 1994
C     Last revision : Sep 11 1995 - tsaue
C***********************************************************************
      use fde_mod
      use dirac_cfg
      use visual

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbpsi.h"
C
      CALL QENTER('PAMDRV')
C
C***********************************************************************
C*****  G E T    W A V E    F U N C T I O N   **************************
C***********************************************************************
C
      CALL PAMPSI()
C
C***********************************************************************
C*****  Write(C1/formatted) MO coefficients   **************************
C***********************************************************************
C
C
C***********************************************************************
C*****  A N A L Y S I S    M O D U L E  ********************************
C***********************************************************************
C
      IF(DOANA .AND. .NOT.
     & ((dirac_cfg_scf_calculation.AND.DOPSI).OR.(DOKRMC.AND.DOPSI)))
     &   CALL PAMANA()
C
C
C***********************************************************************
C*****  R E S P O N S E    M O D U L E  ********************************
C***********************************************************************
C
      IF(DOPRP) CALL PAMPRP()

C
C***********************************************************************
C***** V I S U A L I Z A T I O N ***************************************
C***********************************************************************
C
      if (dirac_cfg_visual) call run_visual()


 666  CONTINUE
      CALL QEXIT('PAMDRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pampsi */
      SUBROUTINE PAMPSI()
C***********************************************************************
C
C     Wave function module
C
C     Written by T.Saue April 1996
C     Revision for MCSCF J. Thyssen - Oct 1998
C     Restructured 2002 by Hans Joergen Aa. Jensen
C
C***********************************************************************
      use dirac_cfg
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
C
      LOGICAL   WFCONV
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
C
      CALL QENTER('PAMPSI')
C     Attempt read of formatted MO coefficients
      CALL PCMOIN()
C     Read C1 MO coefficients
      IF(DOACIN) CALL ACMOIN()
C     Extend RKB coeffcients with complementary UKB space
      IF(DORKBIMP) CALL RKBIMP(IPRGEN)
C     Reorder/rotate orbitals prior to SCF-calculations
      IF(L1ORBM) CALL PREORB()
C
C***********************************************************************
C**** W A V E     F U N C T I O N S   **********************************
C**** before 4-index transformation   **********************************
C***********************************************************************
C
      IF(DOPSI) THEN
C
C***********************************************************************
C******* D I R A C - H A R T R E E - F O C K ***************************
C******* o r   D I R A C - K O H N - S H A M ***************************
C***********************************************************************
        IF(dirac_cfg_scf_calculation) THEN
C
C         Transfer DHF info to common block DCBORB and recalculate DCBORB.
C         Open shell is defined according to 'RESOLVE', i.e. include
C         in CI exactly those determinants included in the av.of.conf. SCF
C
          CALL TDHFORB('NO CI')
C
C         *************************************************
C         *** Call PSISCF for Roothaan/DIIS iterations  ***
C         *************************************************
C
          IF (MAXITR.GT.0.OR.MAXITR2.GT.0) THEN
             CALL PSISCF(WFCONV)
          END IF
C
C         *************************************************
C         *** Call PSIOPT for second order optimization ***
C         *************************************************
C
          IF ( .NOT. NOQCDHF .AND.
     &         .NOT. dirac_cfg_dft_calculation .AND.
#ifdef HAS_PCMSOLVER
     &         .NOT. dirac_cfg_pcm .AND.
#endif
     &        (.NOT. WFCONV .AND. MAXITR .GT.0 )) THEN
             CALL PSIOPT('SCF  ',WFCONV)
          END IF
C
          IF(DOANA) CALL PAMANA()
          IF(.NOT.WFCONV) THEN
C         ... DHF not converged;
             WRITE (LUPRI,'(//A)')
     &      ' ERROR, SCF not converged, Dirac stops!'
             CALL QUIT('SCF not converged')
          END IF
        END IF
C       Post orbital modifications
        IF(L2ORBM) CALL PSTORB()
        IF(DOACUT) CALL ACMOUT()
        IF(DOPUT) CALL PCMOUT()
C
C***********************************************************************
C***** R E S O L V E    O P E N - S H E L L  S T A T E S ***************
C***********************************************************************
C
        IF(DORES) THEN
           IF (.NOT.dirac_cfg_scf_calculation) THEN
C             ... Transfer DHF info to common block DCBORB
C             ... and recalculate DCBORB.
              CALL TDHFORB('RESOLVE')
           END IF
           CALL RESOLV()
        END IF
C
C***********************************************************************
C***** M O E L L E R - P L E S S E T ***********************************
C***********************************************************************
C
        IF(DOMP2) CALL PSIMP2()
C
C***********************************************************************
C***** M O D I F I E D   V I R T U A L   O R B I T A L S ***************
C***********************************************************************
C
        IF(DOMVO) CALL MVOGEN()
C
C***********************************************************************
C***** M P 2   N A T U R A L   O R B I T A L S *************************
C***********************************************************************
C
        IF(DOMP2NO) CALL MP2NOGEN()
C
C***********************************************************************
C***** K R - M C S C F *************************************************
C***********************************************************************
C       Do MCSCF before CC and CI because user may want
C       to do an MR-CC or MR-CI /hjaaj aug 03
C
        IF(DOKRMC) THEN
C
C          Transfer KRMC information from common block DCOKRMC to
C          common block DCBORB and set common block DCxOPT
C
           CALL TKRMCORB
C
C          *************************************************
C          *** Call PSIOPT for KR-MCSCF wave function ******
C          *************************************************
C
C          If MCSCF-srDFT calculation ...
C          set dirac_cfg_mcsrdft_calculation to TRUE in order to
C          calculate proper 2-electron integrals
           if (hfxmu .ne. 0.0d0) then
             dirac_cfg_mcsrdft_calculation = .true.
           endif
C
           CALL PSIOPT('MCSCF',WFCONV)
           IF(DOANA) CALL PAMANA()
        END IF
      ENDIF  ! IF(DOPSI) ...
C
C***********************************************************************
C***** 4 - I N D E X  T R A N S F O R M A T I O N **********************
C***********************************************************************
C
      IF(DOTRA) CALL PAMTRA()
C
C***********************************************************************
C**** W A V E     F U N C T I O N S   **********************************
C**** after 4-index transformation   ***********************************
C***********************************************************************
      IF(DOPSI) THEN
C
C***********************************************************************
C******* R E L C C S D  M O D U L E ************************************
C***********************************************************************
C
        IF(DOCCM) CALL PAMCCM()
#if defined (MOD_HSCC)
!       Hilbert Space Coupled Cluster module
        IF (DOHSCC) CALL PAMHSMRCC()
#endif
#if defined (MOD_HSFS) || defined (HSFS_DYN_ALLOC)
!       High Sectors Fock Space Coupled Cluster module
        IF (DOHSFS) CALL PAMHSFSCC()
#endif
C
C***********************************************************************
C******* D I R R C I   M O D U L E *************************************
C***********************************************************************
C
        IF(DOCIM) CALL PAMCIM()

C***********************************************************************
C******* L U C I T A  M O D U L E **************************************
C***********************************************************************
C
        IF(DOLUCT) CALL PAMLUCITA()
C
C***********************************************************************
C******* R E L A D C  M O D U L E **************************************
C***********************************************************************
C
        IF(DOADC) CALL PAMADC()
C
C***********************************************************************
C******* C O S C I        stand-alone version of RESOLVE   *************
C*******                  some routines in krmc are merged *************
C*******                  to get efficiency and make it    *************
C*******                  possible to calculate density    *************
C*******                  matrix.                          *************
C*******                  added by sya,2007.01.23,ulp      *************
C***********************************************************************
C
        IF(DOCOSCI) CALL SY_COSCI()
C
C***********************************************************************
C***** C I *************************************************************
C***********************************************************************
C
        IF(DOKRCI) THEN
C
C         Transfer KR-CI info to common block DCBORB
C         and recalculate DCBORB.
C        ... reset dcbopt.h
           CALL OPTINI
           CALL TKRCIORB
           CALL SETDC2(0)
           CALL KRCI_CALC('KR-CI ',WFCONV)
        END IF
C
C***********************************************************************
C******* C C ***********************************************************
C***********************************************************************
C
        IF(DOKRCC.AND.DOPSI) THEN
#ifdef MOD_KRCC
C
C        Transfer KR-CI info to common block DCBORB
C        and recalculate DCBORB.
C
C        ... reset dcbopt.h
C
C Not sure which to include?
C
           CALL OPTINI
           CALL TKRCCORB
           CALL SETDC2(0)
           CALL PAMKRCC()
#else
           CALL QUIT(
     &  ' KR-CC module not included in this version')
#endif
        END IF
      ENDIF ! IF(DOPSI) .. second part, after 4-index transformation
C
      CALL QEXIT('PAMPSI')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tdhforb */
      SUBROUTINE TDHFORB(CISELECT)
C***********************************************************************
C
C     Transfer information from common block DCODHF to
C     common block DCBORB and set common block DCxOPT
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
      CHARACTER*(*) CISELECT
C
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbesr.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      CALL QENTER('TDHFORB')

C
C     initialize CI program in charge on common block DCCOPT
C
      OPT_CIPROGRAM  = 'not set'
C
C     Transfer variables to /DCBORB/
C
      NASHT = 0
      DO I = 1,NFSYM
         NISH(I) = NISH_DHF(I)
         NASH(I) = NASH_DHF(I)
         NASHT   = NASHT + NASH(I)
         NOCC(I) = NISH(I) + NASH(I)
         NSSH(I) = NSSH_DHF(I)
c        NESH(I) = NESH_DHF(I)
c        NPSH(I) = NPSH_DHF(I)
         NFRO(I) = NFRO_DHF(I)
         NELEC(I) = NELEC_DHF(I)
      END DO
      NELECT = NELECT_DHF
      NAELEC = NAELEC_DHF
      THRPCI = 1.0D-3 ! threshold for print of CI vector elements
C
C     Calculate derived variables in /DCBORB/:
C
      IF ( CISELECT .EQ. 'NO SETDC2') GOTO 9999 ! for Trond's temp. fix in dirrdn.F
      CALL SETDC2(0)
c     call setdc2(iprint)
C
C     Transfer open shell average-of-configuration to GAS space notation
C     (for e.g. ESR GASCIP) /July 2002 hjaaj
C     HOW it is transferred is determined by the CISELECT keyword. /hjaaj
C
      IF ( CISELECT .EQ. 'NO CI' ) THEN
         NGAS_DC  = 0
         OPT_THR_CVEC = D0
      ELSE IF ( CISELECT .EQ. 'RESOLVE' ) THEN
C        ... include exactly those determinants which was included
C            in the 'average-of-configuration' SCF
         NGAS_DC  = NOPEN
         OPT_THR_CVEC = 1.0D-3
         CALL IZERO(NGSHT,MXGAS)
         DO J = 1, NGAS_DC
            DO I = 1,NFSYM
               NGSH(I,J) = NACSH(I,J)
               NGSHT(J)  = NGSHT(J) + NGSH(I,J)
            END DO
         END DO
C
C        initialize min and max accumlated no. of electrons in each gas space
C
         XELEC = D0
         DO J = 1, NGAS_DC
            XELEC = XELEC + D2*DF(J)*NGSHT(J)
            NGASSP(1,J) = NINT(XELEC)
            NGASSP(2,J) = NINT(XELEC)
         END DO
      ELSE IF (CISELECT .EQ. 'FULLCI') THEN
C        ... do a full CI in the active space
         NGAS_DC  = 1
         OPT_THR_CVEC = 1.0D-3

C        ... set CI program in charge
!        OPT_CIPROGRAM = 'GASCIP'

         CALL IZERO(NGSHT,MXGAS)
         DO J = 1, NGAS_DC
            DO I = 1,NFSYM
               NGSH(I,J) = NASH(I)
               NGSHT(J)  = NGSHT(J) + NGSH(I,J)
            END DO
         END DO
C
C        initialize min and max accumlated no. of electrons in each gas space
C
         NGASSP(1,1) = NAELEC
         NGASSP(2,1) = NAELEC
      ELSE IF (CISELECT .EQ. 'ESRCI') THEN
C        ... do CI for ESR spin-polarization /hjaaj nov 2002

C        ... set CI program in charge
         OPT_CIPROGRAM = 'GASCIP'
C
         OPT_THR_CVEC = THR_CVEC_ESR
         THRPCI   = THRPCI_ESR
C
C        Set up core (NCSH_ESR), valence (NVSH_ESR) and virtual (secondary, NSSH_ESR) orbital shells.
C
         NCSH_ESRT = 0
         NVSH_ESRT = 0
         NSSH_ESRT = 0
         NELECT    = 0
         DO I = 1,NFSYM

           IF (NCSH_ESR(I) < 0) THEN
              CALL QUIT('ERROR: NCSH_ESR(I) < 0')
           END IF
           IF (NVSH_ESR(I) < 0) THEN
              CALL QUIT('ERROR: NVSH_ESR(I) < 0')
           END IF
           IF (NSSH_ESR(I) < 0) THEN
              CALL QUIT('ERROR: NSSH_ESR(I) < 0')
           END IF

           NISH(I)   = NISH(I)   - NCSH_ESR(I) - NVSH_ESR(I)
           NSSH(I)   = NSSH(I)   - NSSH_ESR(I)

           NCSH_ESRT = NCSH_ESRT + NCSH_ESR(I)
           NVSH_ESRT = NVSH_ESRT + NVSH_ESR(I)
           NSSH_ESRT = NSSH_ESRT + NSSH_ESR(I)

           NELEC(I)  = NELEC(I) - 2*NCSH_ESR(I) - 2*NVSH_ESR(I)
           NELECT    = NELECT   + NELEC(I)

         END DO
         NGAS12_ELEC = 2*NCSH_ESRT + 2*NVSH_ESRT
         NAELEC      = NAELEC + NGAS12_ELEC
C
C        initialize min and max accumulated no. of electrons in each gas space
C
         CALL IZERO(NGSHT,MXGAS)
         NGAS_DC = 0
         IF (NCSH_ESRT > 0) THEN
            NGAS_DC = NGAS_DC + 1
            NGASSP(1,NGAS_DC) = MAX(0,2*NCSH_ESRT-LVL_CIESR(0))
            NGASSP(2,NGAS_DC) = 2*NCSH_ESRT
            DO I = 1,NFSYM
               NGSH(I,NGAS_DC) = NCSH_ESR(I)
               NGSHT(NGAS_DC)  = NGSHT(NGAS_DC) + NCSH_ESR(I)
            END DO
         END IF
         IF (NVSH_ESRT > 0) THEN
            NGAS_DC = NGAS_DC + 1
            NGASSP(1,NGAS_DC) = MAX(0,NGAS12_ELEC-LVL_CIESR(1))
            NGASSP(2,NGAS_DC) = NGAS12_ELEC
            DO I = 1,NFSYM
               NGSH(I,NGAS_DC) = NVSH_ESR(I)
               NGSHT(NGAS_DC)  = NGSHT(NGAS_DC) + NVSH_ESR(I)
            END DO
         END IF

C            (Cxx are preparation for RESOLVE instead of FULLCI
C             in the DHF open shells)
Cxx      XELEC = D2*NCSH_ESRT
         IF (NOPEN .NE. 1) THEN
            write (lupri,*) 'ERROR ESRCI: NOPEN is not 1, it is: ',NOPEN
            call quit ('ESRCI: sorry, NOPEN is not 1')
         END IF
         DO J = 1, NOPEN
            NGAS_DC = NGAS_DC + 1
Cxx         XELEC = XELEC + D2*DF(J)*NGSHT(J)
Cxx         NGASSP(1,NGAS_DC) = NINT(XELEC)
Cxx         NGASSP(2,NGAS_DC) = NINT(XELEC)
            NGASSP(1,NGAS_DC) = MAX(0,NAELEC-LVL_CIESR(2))
            NGASSP(2,NGAS_DC) = NAELEC
            DO I = 1,NFSYM
Cxx            NGSH(I,NGAS_DC) = NACSH(I,J)
               NGSH(I,NGAS_DC) = NASH(I)
               NGSHT(NGAS_DC)  = NGSHT(NGAS_DC) + NASH(I)
            END DO
         END DO
         IF (NSSH_ESRT .GT. 0) THEN
            NGAS_DC = NGAS_DC + 1
            NGASSP(1,NGAS_DC) = NAELEC
            NGASSP(2,NGAS_DC) = NAELEC
            DO I = 1,NFSYM
               NGSH(I,NGAS_DC) = NSSH_ESR(I)
               NGSHT(NGAS_DC)  = NGSHT(NGAS_DC) + NSSH_ESR(I)
            END DO
         ELSE
            NGASSP(1,NGAS_DC) = NAELEC ! reset open shell minimum, not possible with NAELEC-LVL_CIESR(2) if not virt gas space
         END IF
C
C        reset number of active orbitals
         DO I = 1,NFSYM
            NASH(I) = 0
            DO J = 1,NGAS_DC
               NASH(I) = NASH(I) + NGSH(I,J)
            END DO
         END DO

      ELSE
         WRITE(LUPRI,'(//2A)')
     &   'ERROR Invalid CISELECT parameter in TDHFORB : ',CISELECT
         CALL QUIT('Invalid CISELECT parameter in TDHFORB')
      END IF
C
C     Define MK limits now AFTER CI has been defined:
C     we may have more active orbitals and (maybe) more
C     active electrons than in the wave function optimization step.
C
C     For "RESOLVE" this corresponds to:
C     transfer open shell average-of-configuration to GAS space notation /July 2002 hjaaj
C
C
      NASHT  = NASH(1) + NASH(2)
      IOPT_MAXMK2 = 2*MIN(NAELEC,NASHT) - NAELEC
      IOPT_MINMK2 = -IOPT_MAXMK2
      IOPT_MK2REF = 0
      IOPT_MK2DEL = IOPT_MAXMK2
C
C     Set common block DCIOPT
C
      IPROPT = IPRSCF
      IOPT_ITRINT(1) = ITRINT(1)
      IOPT_ITRINT(2) = ITRINT(2)
      IOPT_INTDEF    = INTDEF
      IOPT_INTFLG    = INTFLG
      IOPT_INTBUF    = INTBUF
      MOPT_MXMICRO   = MXMICRO
      MOPT_MXMACRO   = MXMACRO
      MOPT_MAXBCK    = MAXBCK
      IOPT_STATE     = 0
      JKRRUNTYPE     = 0
      ISVRONO(1)     = 0
      ISVRONO(2)     = 0
C
C     Set common block DCROPT
C
      OPT_CNVINT(1)  = CNVINT(1)
      OPT_CNVINT(2)  = CNVINT(2)
      IF (EVCCNV) THEN
         OPT_THRGRD     = SCFCNV(1)
      ELSE
         OPT_THRGRD     = SQRT(SCFCNV(1))
      END IF
C
C     Set common block DCLOPT
C
      OPT_SKIPEE     = DHF_SKIPEE
      OPT_SKIPEP     = DHF_SKIPEP
      OPT_NOCI       = .FALSE.
      OPT_UCIBOS     = .FALSE.
      OPT_CHCKJZ     = .FALSE.
      IF( LINEAR ) OPT_CHCKJZ = .TRUE.
      no1pdens_save  = .true.
      save_reordered_nos = .false.
      fcidump            = .false.
C
C     Set common block DCCOPT
C
!     OPT_CIPROGRAM  = 'not set'
      DO I = 1, NFSYM
         OPT_FRZSTR(I) = ' '
      END DO
      DO I = 1, NFSYM
         OPT_DELSTR(I) = ' '
      END DO
C
 9999 CONTINUE
      CALL QEXIT('TDHFORB')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tkrmcorb */
      SUBROUTINE TKRMCORB
C***********************************************************************
C
C     Transfer information from common block DCOKRMC to
C     common block DCBORB and set common block DCxOPT
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxash.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbkrmc.h"
#include "dcbopt.h"
#include "dcbgascip.h"
#include "dcbdhf.h"
      integer :: nasht_tmp, igas_tmp, accum_act(2)
C
      CALL QENTER('TKRMCORB')
C
C
      NGAS_DC= NKRMCGAS
      IOPT_MK2REF = NKRMC_MK2REF

!     set active spaces and other derived variables for the full CI option
      if(KRMC_full_ci)then
        nasht_tmp = 0
        DO I = 1,NFSYM
          NKRMCASH(I) = NESH(I) - NKRMCISH(I)
          nasht_tmp   = NESH(I) - NKRMCISH(I) + nasht_tmp
        END DO
        NKRMC_MAXMK2 = 2*MIN(NKRMCAELEC,nasht_tmp) - NKRMCAELEC
        NKRMC_MINMK2 = NKRMCAELEC-2*MIN(NKRMCAELEC,nasht_tmp)
        NKRMC_MK2DEL = NKRMC_MAXMK2
        igas_tmp     = 0
        accum_act(1) = 0
        accum_act(2) = 0
        do
          igas_tmp = igas_tmp + 1
          if(NKRMCGSH(1,igas_tmp) == -9999)then
            do i = 1, nfsym
              NKRMCGSH(I,igas_tmp) = NKRMCASH(I) - accum_act(i)
            end do
            exit
          else
            do i = 1, nfsym
              accum_act(i) = accum_act(i) + NKRMCGSH(I,igas_tmp)
            end do
          end if
        end do
!       print some extra info
        write(lupri,*)
     &  '  *** extra info from TKRMCORB for full CI option ***'
        write(lupri,*)
     &  '  ---------------------------------------------------'
        write(lupri,*)
     &  '     total # of active orbitals: ', nasht_tmp

        do i = 1, ngas_dc
          WRITE(LUPRI,'(A,I3,A,2I4)')
     &         '   - GAS space ',I,'       : ',
     &         (NKRMCGSH(J,I),J=1,NFSYM)
          WRITE(LUPRI,'(A,I3,A1,I3,A)')
     &    '    (constraints: min/max active electrons after space : ',
     &         NKRMCGSP(1,I),'/',NKRMCGSP(2,I),')'
        end do
        write(lupri,*)
     &  '  ---------------------------------------------------'

!       test for maximum # of orbitals
!       if(nasht_tmp > 5)then
        if(nasht_tmp > MAXASH)then
          WRITE(LUPRI,'(A,I5/10X,A,I4)')
     &         ' *** ERROR in *TKRMCORB *** Too many active orbitals: ',
     &         nasht_tmp, ' -- max is ',MAXASH
          CALL QUIT('*** ERROR in *TKRMCORB: Too many active orbitals!')
        end if

      end if ! full ci option

      IOPT_MK2DEL = NKRMC_MK2DEL
      IOPT_MINMK2 = NKRMC_MINMK2
      IOPT_MAXMK2 = NKRMC_MAXMK2
      CALL IZERO(NGSHT,MXGAS)
      DO I = 1,NFSYM
         NISH(I) = NKRMCISH(I)
         NASH(I) = NKRMCASH(I)
         DO J = 1, NGAS_DC
            NGSH(I,J) = NKRMCGSH(I,J)
            NGSHT(J) = NGSHT(J) + NGSH(I,J)
         END DO
         NOCC(I) = NISH(I) + NASH(I)
         NSSH(I) = NKRMCSSH(I)
c         NESH(I) = NKRMCESH(I)
c         NPSH(I) = NKRMCPSH(I)
         NFRO(I) = NKRMCFRO(I)
      END DO
      NAELEC = NKRMCAELEC
C     min and max electrons in each gas space:
      DO J = 1, NGAS_DC
         NGASSP(1,J) = NKRMCGSP(1,J)
         NGASSP(2,J) = NKRMCGSP(2,J)
      END DO
C
C     Calculated derived variables in /DCBORB/
C
      CALL SETDC2(0)
c     call setdc2(iprint)
C
C
C     Set common block DCIOPT
C
      IPROPT = IPRKRMC
      IOPT_ITRINT(1) = IKRMC_ITRINT(1)
      IOPT_ITRINT(2) = IKRMC_ITRINT(2)
      IOPT_INTDEF    = IKRMC_INTDEF
      IOPT_INTFLG    = IKRMC_INTFLG
      IOPT_INTBUF    = IKRMC_INTBUF
      MOPT_MXMICRO   = MKRMC_MXMICRO
      MOPT_MXMACRO   = MKRMC_MXMACRO
      MOPT_MAXBCK    = MKRMC_MAXBCK
      IOPT_STATE     = IKRMC_STATE
      IOPT_SYMMETRY  = IKRMC_SYMMETRY
      IMEMFAC        = IKRMC_MEMFAC
      JKRRUNTYPE     = 1
!     information from dcbdhf.h
      nopen_mc       = nopen
      ISVRONO(1)     = IKRMC_SVRONO(1)
      ISVRONO(2)     = IKRMC_SVRONO(2)
C
C     Set common block DCROPT
C
      OPT_CNVINT(1)  = DKRMC_CNVINT(1)
      OPT_CNVINT(2)  = DKRMC_CNVINT(2)
      OPT_THRGRD     = DKRMC_THRGRD
      OPT_THRECI     = MAX(DKRCI_THRECI, OPT_THRGRD**2)
      THRPCI         = DKRMC_THRPCI
      E_LVCORR_MC    = 0.0d0
      OPT_THR_CVEC   = 1.0D-3
C
C     Set common block DCLOPT
C
      OPT_SKIPEE         = KRMC_SKIPEE
      OPT_SKIPEP         = KRMC_SKIPEP
      OPT_NOCI           = KRMC_NOCI
      OPT_UCIBOS         = KRMC_UCIBOS
      OPT_CHCKJZ         = KRMC_CHCKJZ
      no1pdens_save      = KRMC_NO1pdens
      save_reordered_nos = krmc_save_reordered_nos
C
      DO I = 1, NOPTFLAGS
         FLAG(I) = MCFLAG(I)
      END DO
C
C     Set common block DCCOPT
C
      OPT_CIPROGRAM  = KRMC_CIPROGRAM
      DO I = 1, NFSYM
         OPT_FRZSTR(I) = KRMC_FRZSTR(I)
      END DO
      DO I = 1, NFSYM
         OPT_DELSTR(I) = KRMC_DELSTR(I)
      END DO
C
C     Set "constant" CI internal variables
C
      NMAX_SYM    = IKRMC_SYMMETRY
      NCIROOT     =    1
      ISHMEM_TYPE =  - 1
      CSHMEMO     = .FALSE.
      CINT_REORD  = .FALSE.
      CINT_SPLIT  = .FALSE.
      CHCKPT_WRT  = .FALSE.
      CPBLCK_FILE = .FALSE.
      CINT_LOWSRT = .FALSE.
      cana_mcscf  = .FALSE.
C
      CALL QEXIT('TKRMCORB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tkrciorb */
      SUBROUTINE TKRCIORB
C***********************************************************************
C
C     Transfer information from common block DCOCI to
C     common block DCBORB and set common block DCxOPT
C
C     Based on TKRMCORB.
C
C     Written by S. Knecht - Aug 2008
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbkrci.h"
#include "dcbopt.h"
#include "dcbgascip.h"
#include "dcbdhf.h"
C
      CALL QENTER('TKRCIORB')
C
C
      NGAS_DC= NKRCIGAS
      IOPT_MK2REF = NKRCI_MK2REF
      IOPT_MK2DEL = NKRCI_MK2DEL
      IOPT_MINMK2 = NKRCI_MINMK2
      IOPT_MAXMK2 = NKRCI_MAXMK2

      CALL IZERO(NGSHT,MXGAS)
      DO I = 1,NFSYM
         NISH(I) = NKRCIISH(I)
         NASH(I) = NKRCIASH(I)
         DO J = 1, NGAS_DC
            NGSH(I,J) = NKRCIGSH(I,J)
            NGSHT(J) = NGSHT(J) + NGSH(I,J)
         END DO
         NOCC(I) = NISH(I) + NASH(I)
         NSSH(I) = NKRCISSH(I)
         NFRO(I) = NKRCIFRO(I)
      END DO
      NAELEC = NKRCIAELEC
C     min and max electrons in each gas space:
      DO J = 1, NGAS_DC
         NGASSP(1,J) = NKRCIGSP(1,J)
         NGASSP(2,J) = NKRCIGSP(2,J)
      END DO
C
C
C     Set common block DCIOPT
C
      IPROPT = IPRKRCI
      IOPT_ITRINT(1) = IKRCI_ITRINT(1)
      IOPT_ITRINT(2) = IKRCI_ITRINT(2)
      IOPT_INTDEF    = IKRCI_INTDEF
      IOPT_INTFLG    = IKRCI_INTFLG
      IOPT_INTBUF    = IKRCI_INTBUF
      MOPT_MXMICRO   =  0
      MOPT_MXMACRO   = -1
      MOPT_MAXBCK    =  0
      IOPT_STATE     =  0
      IOPT_SYMMETRY  = -1
      IMEMFAC        = IKRCI_MEMFAC
      JKRRUNTYPE     = 2
!     information from dcbdhf.h
      nopen_mc       = nopen
      ISVRONO(1)     = IKRCI_SVRONO(1)
      ISVRONO(2)     = IKRCI_SVRONO(2)
C
C     Set common block DCROPT
C
      OPT_CNVINT(1)  = DKRCI_CNVINT(1)
      OPT_CNVINT(2)  = DKRCI_CNVINT(2)
      THRPCI         = DKRCI_THRPCI
      E_LVCORR_MC    = 0.0d0
      OPT_THR_CVEC    = 1.0d-03
      OPT_THRGRD     = 1.D99
      OPT_THRECI     = DKRCI_THRECI
C
C     Set common block DCLOPT
C
      OPT_SKIPEE     = .FALSE.
      OPT_SKIPEP     = .FALSE.
!     if number of positronic shells is zero skip e-p rotations. hjj+sk - aug 2010
      if(npsht .eq. 0) OPT_SKIPEP = .TRUE.
      OPT_NOCI       = .FALSE.
      OPT_UCIBOS     = KRCI_UCIBOS
      OPT_CHCKJZ     = .FALSE.
!     required for KR-CI in linear symmetry set-up
      if(linear) opt_chckjz = .true.
!     opt_chckjz = .true.
      no1pdens_save  = .true.
      save_reordered_nos = save_reordered_nos_krci
C
      DO I = 1, NOPTFLAGS
         FLAG(I) = .FALSE.
      END DO
C
C     Set common block DCCOPT
C
      OPT_CIPROGRAM  = KRCI_CIPROGRAM
      DO I = 1, NFSYM
         OPT_FRZSTR(I) = ' '
      END DO
      DO I = 1, NFSYM
         OPT_DELSTR(I) = ' '
      END DO
C
      NMAX_SYM    = NKRCI_MAX_SYM
      MAXCIT      = MAXCIT_KRCI
      KTRLVL      = KTRLVL_KRCI
      IANACI      = IANACI_KRCI
      NCIROOT     = -1
      MXCIV       = MXCIV_KRCI
      IRESTRK     = IRESTRK_KRCI
      IDOPARIO    = IDOPARIO_KRCI
      IOPTST      = 1563
      IKRMCCNO    =  - 3
      OPT_NOPFQ   = .FALSE.
      OPT_NOFQX   = .FALSE.
      NATONL      = .FALSE.
      FOCKON      = .FALSE.
      COMPROT     = .FALSE.
      GENFOCK     = GENFOCK_KRCI
      NATOLCR     = NATOLCR_KRCI
      TRA_NATO    = TRANATO_KRCI
      CRDFO_MAT   = CRDFO_MAT_KRCI
      CWRTFO_MAT  = CWRTFO_MAT_KRCI
      CSHMEMO     = CSHMEMO_KRCI
      CINT_REORD  = CINT_REORD_KRCI
      CINT_SPLIT  = CINT_SPLIT_KRCI
      CINT_LOWSRT = LOWSORT_KRCI
      cana_mcscf  = domcana_krci
      COMPFAC     = 0.01D00
      CHCKPT_WRT  = CHCKPT_KRCI
      CPBLCK_FILE = CUSE_PBLFL_KRCI
      fcidump     = fcidump_krci
C
      CALL QEXIT('TKRCIORB')
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirnod */
      SUBROUTINE DIRNOD(IPARCAL,IMYTID,IMPARID,INUMNOD)
C***********************************************************************
C
C     Menu file for slave nodes. Check whether the program needs a slave
C     and determine the type of work to be done.
C
C     Written by L.Visscher July 1997
C     Last revision: July 29 1996
C
C***********************************************************************
      use interface_to_mpi
#ifdef MOD_INTEREST
      use module_interest_interface,                                    &
     &    only: initialize_interest_herd,                               &
     &          process_two_electron_integrals_4comp_slaves,            &
     &          process_two_electron_integrals_2comp_slaves 
!     use fde_evaluators_dirac, only : fde_launch_slave_process
      use fde_mod, only : fde_launch_slave_process
#endif
      implicit none
#include "dcbdhf.h"
#include "dcbgen.h"
#include "priunit.h"
#include "infpar.h"
      DIMENSION ISTAT(df_mpi_status_size)
      LOGICAL NEWTASK
      integer :: IMPARID, INUMNOD, IMYTID, ITASK, IPARCAL
      real*8 :: CPU1, WALL1
C
C
      integer :: sleep_time=15

      CALL QENTER('DIRNOD')
C
C     Transfer variables to common blocks
C
      MPARID = IMPARID
C We count the number of slave nodes, subtract the master
C     NUMNOD = INUMNOD - 1
C     SK - Sep 08: subtraction has already been done in mpixinit.
      NUMNOD = INUMNOD
      MYTID  = IMYTID
      PARCAL = (IPARCAL .NE. 0)

#if defined (SLOWDOWN_SLAVES)
!mi+sk: introduce slowdown statement for debug purposes :
! node 'sleeps' short time and in between the user can attach debugger to the process
      if(mytid > 0)then
         write(lupri,'(2x,a,i2,a)')
     &   'slaves first sleep for ',sleep_time,' seconds. '
         call sleep(sleep_time)
      endif
#endif

C
C     Set basis set directory for this slave
C     (it may be needed for reading the .mol file)
C
      CALL GETBASDIR()
C
 100  CONTINUE
#if defined (VAR_SLEEPING_SLAVES)
! If it is important that the slaves are not idling
! then one can define VAR_SLEEPING SLAVES.
! If set, then the slaves test once per second if there
! is a new task (because of the CALL SLEEP(1))
      call interface_mpi_iprobe(MPARID,27,NEWTASK,global_communicator,
     &                          ISTAT)
      IF (NEWTASK) THEN
chj      write (0,*) 'mpi_iprobe true for node',mytid
chj      write (0,*) 'istatus(source,tag,error):',
chj  &     ( istatus(i),i=1,3)
         CALL interface_MPI_RECV(ITASK,1,MPARID,27,
     &                           global_communicator)
chj      write (0,*) ' itask, node',itask,mytid
chj      write (0,*) 'istatus(source,tag,error):',
chj  &     ( istatus(i),i=1,3)
      ELSE
         CALL SLEEP(1)
         GO TO 100
      END IF
#else
C
Csk      write(6,*) 'node will recieve itask',MYTID
         CALL interface_MPI_RECV(ITASK,1,MPARID,27,
     &                           global_communicator)
Csk      write(6,*) 'node has recieved itask',MYTID,itask
#endif
C
C     check for proper LUCITA / LUCIAREL installation
C
#if !defined (VAR_MPI2)
C
      IF (ITASK .eq. 5 .or. ITASK .eq. 6) THEN
        WRITE(6,*) ' fatal error in starting a parallel LUCITA /
     & LUCIAREL run; you do not have a MPI-2 library installed'
       CALL Abend2(' error in DIRNOD: missing MPI-2 library functions')
      END IF
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
C        Start calculating integrals for Fock matrices.
C
         CALL HERNOD()
      ELSE IF (ITASK .EQ. 2) THEN
C
C        Start calculating integrals for MOLTRA.
C
          CALL TRANOD()
      ELSE IF (ITASK .EQ. 3) THEN
C
C        Start coupled cluster calculation.
C
         CALL CCNOD()
      ELSE IF (ITASK .EQ. 4) THEN
C        Enter the DFT module
         call xcint_launch_slave_process()

#if defined (VAR_MPI2)
      ELSE IF (ITASK .EQ. 5) THEN
C
C        Enter the LUCITA module
C
         CALL LUCITA_NODE()
      ELSE IF (ITASK .EQ. 6) THEN
C
C        start KRMC-LUCIAREL calculation.
C
         CALL LUCIAREL_NODE()
      ELSE IF (ITASK .EQ. 7) THEN
C
C        start KR-MCSCF calculation.
C
         CALL KRMC_NODE_DRIVER()
      ELSE IF (ITASK .EQ. 8) THEN
C
C        initialize the InteRest two-electron MPI calculation
C
#ifdef MOD_INTEREST
         call initialize_interest_herd 
#else
         call quit('ERROR: InteRest lib not available')
#endif
      ELSE IF (ITASK .EQ. 9) THEN
C
C        start the InteRest two-electron MPI calculation (4-comp version)
C
#ifdef MOD_INTEREST
         call process_two_electron_integrals_4comp_slaves 
#else
         call quit('ERROR: InteRest lib not available')
#endif
      ELSE IF (ITASK .EQ. 10) THEN
C
C        start the InteRest two-electron MPI calculation (2-comp version)
C
#ifdef MOD_INTEREST
         call process_two_electron_integrals_2comp_slaves 
#else
         call quit('ERROR: InteRest lib not available')
#endif
#endif /* VAR_MPI2 */
      ELSE IF (ITASK .EQ. 11) THEN
C
C        start GASCIP MPI calculation
C
         CALL GASCIP_NODE()
      ELSE IF (ITASK .EQ. 12) THEN
C
C        start FDE MPI calculation
C
         call fde_launch_slave_process()
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
C
      GOTO 100
C
 999  CALL QEXIT('DIRNOD')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIRNOD2()
C***********************************************************************
C
C     Another menu file for slave nodes. This can be used to change the type
C     of work during a calculation (i.e. : if one parallel module calls the
C     other like in RELCCSD when doing the AO Lagrangian).
C
C***********************************************************************
      use interface_to_mpi
      use fde_mod, only : fde_launch_slave_process
#include "implicit.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "priunit.h"
C
#include "infpar.h"
      LOGICAL NEWTASK
C
C
      CALL QENTER('DIRNOD2')
C
 100  CONTINUE
C
      CALL interface_MPI_RECV(ITASK,1,MPARID,27,
     &                        global_communicator)
C
C     check for proper LUCITA / LUCIAREL installation
C
#if !defined (VAR_MPI2)
C
      IF (ITASK .eq. 5 .or. ITASK .eq. 6) THEN
        WRITE(6,*) ' fatal error in starting a parallel LUCITA /
     & LUCIAREL run; you do not have a MPI-2 library installed'
       CALL Abend2(' error in DIRNOD2: missing MPI-2 library functions')
      END IF
#endif
C
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
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
C        Start calculating integrals for Fock matrices.
C
         CALL HERNOD()
      ELSE IF (ITASK .EQ. 2) THEN
C
C        Start calculating integrals for MOLTRA.
C
          CALL TRANOD()
      ELSE IF (ITASK .EQ. 3) THEN
C
C        Start coupled cluster calculation.
C
         CALL CCNOD()
      ELSE IF (ITASK .EQ. 4) THEN
C        Enter the DFT module
         call xcint_launch_slave_process()

#if defined (VAR_MPI2)
      ELSE IF (ITASK .EQ. 5) THEN
C
C        Enter the LUCITA module
C
         CALL LUCITA_NODE()
      ELSE IF (ITASK .EQ. 6) THEN
C
C        start KRMC-LUCIAREL calculation.
C
         CALL LUCIAREL_NODE()
      ELSE IF (ITASK .EQ. 7) THEN
C
C        start KR-MCSCF calculation.
C
         CALL KRMC_NODE_DRIVER()
#endif /* VAR_MPI2 */
      ELSE IF (ITASK .EQ. 11) THEN
C
C        start GASCIP MPI calculation
C
         CALL GASCIP_NODE()
      ELSE IF (ITASK .EQ. 12) THEN
C
C        start FDE MPI calculation
C
         call fde_launch_slave_process()

      ELSE
C
C        Unknown task. Abort execution.
C
         WRITE(LUPRI,'(A,I5)')
     &     '*** ERROR *** DIRNOD2: Unknown task ITASK = ',ITASK
         CALL QUIT(' Unknown task in DIRNOD2')
      END IF
C
C     We want the slave program to go until there is nothing
C     more to calculate, K.Ruud, May-95
C
      GOTO 100
C
 999  CALL QEXIT('DIRNOD2')
      END
#endif /* VAR_MPI */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirac_parctl */
      SUBROUTINE DIRAC_PARCTL(ITASK)
C
C 5-July-2000 Hans Joergen Aa. Jensen
C
C Called by master to control parallel calculation.
C Can be used to send slaves to sleep and to wake them up.
C
C ITASK .gt. 0: wake slaves and get hold of them for task ITASK
C ITASK .eq. 0: send slaves to sleep ("0" makes sure they are in menu)
C ITASK .lt. 0: release slaves from current task and send them to sleep
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
C
C
#if defined (VAR_MPI)
#include "infpar.h"
C

      IF (NUMNOD .EQ. 0) RETURN   ! MPI version of DIRAC, but no slaves (i.e. sequential)

      IF (ITASK .EQ. -1) THEN
C        go back from submenu to main menu via MPI_BCAST
         CALL interface_MPI_BCAST(ITASK,1,MPARID,
     &                            global_communicator)
      ELSE
C        slave should be in main menu, ready to receive task
         DO NODE = 1,NUMNOD
            CALL interface_MPI_SEND(ITASK,1,NODE,27,
     &                              global_communicator)
         END DO
      END IF

#else /* not VAR_MPI */

      CALL QUIT('DIRAC_PARCTL called '//
     &   'but this is a sequential Dirac version !')

#endif /* VAR_MPI */
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getbasdir */
      SUBROUTINE GETBASDIR()
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
#include "infpar.h"
#include "dcbgen.h"
C
C     Check environment variable BASDIR
C
#if defined (VAR_MPI)
      IF (PARCAL) THEN
        CALL DIST_ENV(0)
      ELSE
#endif /* VAR_MPI */
#if defined (SYS_T3D) || defined (SYS_T90)
        CALL PXFGETENV ('BASDIR',6,BASDIR,LBASDR,IERR)
#else
        CALL GETENV ('BASDIR',BASDIR)
#endif
#if defined (VAR_MPI)
      ENDIF
#endif /* VAR_MPI */
      RETURN
      END
#ifdef MOD_KRCC
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tkrccorb */
      SUBROUTINE TKRCCORB
C***********************************************************************
C
C     Transfer information from common block DCOCI to
C     common block DCBORB and set common block DCxOPT
C
C     Based on TKRMCORB.
C
C     Written by S. Knecht - Aug 2008.
C     Well Lasse did a copy and paste of it to this routine 2010
C     Last revision:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
C#include "dcbkrci.h"
#include "dcbopt.h"
#include "dcbgascip.h"
#include "dcbibt.h"
C Scratch
      INTEGER NKRCCISH(NFSYM),NKRCCASH(NFSYM),NKRCCGSH(2,MXGAS)
      INTEGER NKRCCSSH(NFSYM),NKRCCFRO(NFSYM),NKRCCGSP(2,MXGAS)
C
      CALL QENTER('TKRCCORB')
C
C be careful with these
C since there are some variables with the same name we will call a routine for this
      CALL TRANSFER_KRCC_TO_KRMC(NKRCCGAS,NKRCC_MK2REF,
     &                           NKRCC_MK2DEL,NKRCCAELEC,NFSYM,
     &                           NKRCC_MINMK2,NKRCC_MAXMK2,NKRCCISH,
     &                           NKRCCASH,NKRCCGSH,NKRCCSSH,NKRCCFRO,
     &                           NKRCCGSP,MAXCIT_KRCC)
      NGAS_DC= NKRCCGAS
      IOPT_MK2REF = NKRCC_MK2REF
      IOPT_MK2DEL = NKRCC_MK2DEL
      IOPT_MINMK2 = NKRCC_MINMK2
      IOPT_MAXMK2 = NKRCC_MAXMK2
      CALL IZERO(NGSHT,MXGAS)
      DO I = 1,NFSYM
         NISH(I) = NKRCCISH(I)
         NASH(I) = NKRCCASH(I)
         DO J = 1, NGAS_DC
            NGSH(I,J) = NKRCCGSH(I,J)
            NGSHT(J) = NGSHT(J) + NGSH(I,J)
         END DO
         NOCC(I) = NISH(I) + NASH(I)
         NSSH(I) = NKRCCSSH(I)
         NFRO(I) = NKRCCFRO(I)
      END DO
      NAELEC = NKRCCAELEC
C     min and max electrons in each gas space:
      DO J = 1, NGAS_DC
         NGASSP(1,J) = NKRCCGSP(1,J)
         NGASSP(2,J) = NKRCCGSP(2,J)
      END DO
C
C
C     Set common block DCIOPT
C
      IPROPT = 0 !IPRKRCI
      IOPT_ITRINT(1) = 1 !IKRCI_ITRINT(1)
      IOPT_ITRINT(2) = 1 !IKRCI_ITRINT(2)
      ILLINT    = IBTAND(INTGEN,1)
      ISLINT    = IBTAND(INTGEN/2,1)
      ISSINT    = IBTAND(INTGEN/4,1)
C IOPT_INTDEF controls the Hamiltonian
      IOPT_INTDEF    = ILLINT + 2*ISLINT + 4*ISSINT !IKRCI_INTDEF
      IOPT_INTFLG    = 0 !IKRCI_INTFLG
      IOPT_INTBUF    = 0 !IKRCI_INTBUF
      MOPT_MXMICRO   =  0
      MOPT_MXMACRO   = -1
      MOPT_MAXBCK    =  0
      IOPT_STATE     =  0
      IOPT_SYMMETRY  = -1
      IMEMFAC        = 9 !IKRCI_MEMFAC
      JKRRUNTYPE     = 2
C
C     Set common block DCROPT
C
      E_LVCORR_MC    = 0.0d0
      OPT_THR_CVEC   = 1.0D-3

      OPT_CIPROGRAM  = 'KRCC'
      DO I = 1, NFSYM
         OPT_FRZSTR(I) = ' '
      END DO
      DO I = 1, NFSYM
         OPT_DELSTR(I) = ' '
      END DO
C
C so far only this sym
      NMAX_SYM    = 1 !NKRCI_MAX_SYM
      MAXCIT      = MAXCIT_KRCC
C transformation level
      KTRLVL      = 3 !KTRLVL_KRCI
      IANACI      = 0 !IANACI_KRCI
      NCIROOT     = -1
C dont think I need max number of CI vectors for subspace
      MXCIV       = 0 !MXCIV_KRCI
C at the moment this is not the way to restart
      IRESTRK     = 0 !IRESTRK_KRCI
C parallel stuff
      IDOPARIO    = 0 !IDOPARIO_KRCI
      IOPTST      = 1563
      IKRMCCNO    =  - 3
      OPT_NOPFQ   = .FALSE.
      OPT_NOFQX   = .FALSE.
      NATONL      = .FALSE.
      FOCKON      = .FALSE.
      COMPROT     = .FALSE.
C at the moment no fock type orbitlas
      GENFOCK     = .FALSE. !GENFOCK_KRCI
      NATOLCR     = .FALSE. !NATOLCR_KRCI
      TRA_NATO    = .FALSE. !TRANATO_KRCI
      CRDFO_MAT   = .FALSE. !CRDFO_MAT_KRCI
      CWRTFO_MAT  = .FALSE. !CWRTFO_MAT_KRCI
      CSHMEMO     = .FALSE. !CSHMEMO_KRCI
      CINT_REORD  = .FALSE. !CINT_REORD_KRCI
      CINT_SPLIT  = .FALSE. !CINT_SPLIT_KRCI
      CINT_LOWSRT = .FALSE. !LOWSORT_KRCI
      cana_mcscf  = .FALSE.
      COMPFAC     = 0.01D00
      CHCKPT_WRT  = .FALSE. !CHCKPT_KRCI
      CPBLCK_FILE = .FALSE. !CUSE_PBLFL_KRCI
      no1pdens_save  = .true.
      fcidump     = .false.
C
      CALL QEXIT('TKRCCORB')
      END
#endif /* MOD_KRCC */
