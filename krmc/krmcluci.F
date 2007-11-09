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
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : krmcluci.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_TRACOM(CIRUN,I_COMM_W_ND)
C***********************************************************************
C
C     Transfer information from DIRAC/KRMC to LUCIAREL
C     common blocks instead of reading input
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision :
C      Setting up real groups in inefficient manner
C      (for quicker progress on testing KRMC-LUCIAREL)
C      Real group integrals are a real subset of complex group
C      integrals, so LUCIAREL internally uses a corresponding
C      complex subgroup according to
C         KRMC     LUCIAREL
C         -----------------
C          D2h      C2h
C          D2       C2
C          C2v      C2
C                               T. Fleig         - Aug 25, 2003
C
C      Only real-valued quantities treated for real matrix groups.
C      Changed energy convergence threshold (THRES_E) to 0.5D-5
C          (complex case) and 1.0D-10 (real case).
C                               T. Fleig         - Oct 09, 2006
C
C      Parallel adpation for MASTER / NODES
C      according to I_COMM_W_ND ( 0 or 1 ) MASTER communicates with
C      NODES 
C                               S. Knecht        - Apr 17, 2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "../luciarel/lucipar_rel.inc"
C
#include "pgroup.h"
#include "dgroup.h"
C
#include "../luciarel/mxpdim.inc"
#include "../luciarel/lucinp.inc"
#include "../luciarel/cstate.inc"
#include "../luciarel/crun.inc"
#include "../luciarel/cprnt.inc"
#include "../luciarel/cgas.inc"
#include "../luciarel/oper.inc"
#include "../luciarel/machine.inc"
#include "../luciarel/cands.inc"
#include "../luciarel/input.inc"
*  nat. orb. occ. num.
#include "../luciarel/noccn_inf.inc"
*
#include "krmcluci_inf.h"
C
      CHARACTER*6 CIRUN
      CHARACTER*72 TITLEC(3)
C
C     say hello: MASTER and NODES
C
CSK      WRITE(6,*) 'I am here in LUCI_TRACOM with CIRUN and flag',MYPROC,
CSK     & CIRUN,I_COMM_W_ND
C
C
C***********************
C  Mandatory settings
C***********************
C
C 1:  TITLE
      TITLEC(1) = ' --------------------------- '
      TITLEC(2) = '    Relativistic MCSCF run   '
      TITLEC(3) = ' --------------------------- '
C
C 2:  Double group
C
#if defined (VAR_MPI)
C
      IF( I_COMM_W_ND .eq. 1 ) THEN
        CALL MPI_BCAST(GROUP,3,MPI_CHARACTER,MASTER,MPI_COMM_WORLD,IERR)
      END IF
C
#endif
C
      IF (GROUP.eq.'C1 ') THEN
         PNTGRP = 5
         NIRREP = 2
         IRIQ = 2
      ELSE IF (GROUP.eq.'Ci ') THEN
         PNTGRP = 6
         NIRREP = 4
         IRIQ = 2
      ELSE IF (GROUP.eq.'Cs ') THEN
         PNTGRP = 7
         NIRREP = 4
         IRIQ = 2
      ELSE IF (GROUP.eq.'C2 ') THEN
         PNTGRP = 8
         NIRREP = 4
         IRIQ = 2
      ELSE IF (GROUP.eq.'C2h') THEN
         PNTGRP = 9
         NIRREP = 8
         IRIQ = 2
      ELSE IF (GROUP.eq.'C2v') THEN
         PNTGRP = 8
         NIRREP = 4
         IRIQ = 1
      ELSE IF (GROUP.eq.'D2 ') THEN
         PNTGRP = 8
         NIRREP = 4
         IRIQ = 1
      ELSE IF (GROUP.eq.'D2h') THEN
         PNTGRP = 9
         NIRREP = 8
         IRIQ = 1
      ELSE
         CALL QUIT('Illegal double group.')
      END IF
      NSMSH = NIRREP/2
      NSMOB = NIRREP
      NIRR_DG_LOC = NIRREP
C
C 18: Reference symmetry
      CALL LUCI_TRAVA(CIRUN,'IREFSM',IREFSM,DUMMY,I_COMM_W_ND)
C
C 19: Roots to be obtained
      CALL LUCI_TRAVA(CIRUN,'ROOTCI',NROOT,DUMMY,I_COMM_W_ND)
C
      DO IRT = 1,NROOT
         IROOT(IRT) = IRT
      END DO
C
C 23: Restart option
      CALL LUCI_TRAVA(CIRUN,'IRESTR',IRESTR,DUMMY,I_COMM_W_ND)
C
C 24: Integral environment
      INTIMP = 7
      ENVIRO(1:6) = 'DIRAC '
C
C 25: Integrals all held in memory
      INCORE = 1
C
C 29-32,36,43,65,89,90,98: All kinds of internal print flags
C     FIXME: The following rescaling should be adapted to the
C            new (when implemented) print flag system in the
C            KRMC! The current is for interface debugging only.
C            Timo
      CALL LUCI_TRAVA(CIRUN,'IPROPT',IPROPT,DUMMY,I_COMM_W_ND)
      IF (IPROPT.EQ.1) THEN
         IPRSTR = 0 * IPROPT
         IPRCIX = 0 * IPROPT
         IPRORB = 0 * IPROPT
         IPRDIA = 0 * IPROPT
         IPRXT  = 0 * IPROPT
         IPROCC = 0 * IPROPT
         IPRDEN = 0 * IPROPT
         IPRSIG = 0 * IPROPT
         IPRTRA = 0 * IPROPT
         IPRHAM = 0 * IPROPT
      ELSE
         IPRSTR = 1 * IPROPT
         IPRCIX = 0 * IPROPT
         IPRORB = 1 * IPROPT
         IPRDIA = 5 * IPROPT
         IPRXT  = 0 * IPROPT
         IPROCC = 2 * IPROPT
         IPRDEN = 2 * IPROPT
         IPRSIG = 10 * IPROPT
         IPRTRA = 1 * IPROPT
         IPRHAM = 100 * IPROPT
      END IF
C
C 33: Maximum number of CI vectors
      CALL LUCI_TRAVA(CIRUN,'MXCIV ',MXCIV,DUMMY,I_COMM_W_ND)
C
C 34: CI storage mode
      ICISTR = 3
C
C 35: Using determinants
      NOCSF = 1
C
C 37: No integrals (we use this to calculate just #determinants)
      CALL LUCI_TRAVA(CIRUN,'NOINT ',NOINT,DUMMY,I_COMM_W_ND)
C
C 39: Dimension of resolution matrices (used in SIGDEN_REL3)
      MXINKA = 50
C
C 46: Core energy from DIRAC
      CALL LUCI_TRAVA(CIRUN,'ECORE ',IDUMMY,ECORE,I_COMM_W_ND)
C
C 50-51: GAS shells and constraints.
      IDOGAS = 1
      NCISPC = 1
      NCMBSPC = NCISPC
      CALL LUCI_TRAVA(CIRUN,'NGAS  ',NGAS,DUMMY,I_COMM_W_ND)
      CALL LUCI_TRAVA(CIRUN,'GASSH ',NGSSH,DUMMY,I_COMM_W_ND)
      CALL LUCI_TRAVA(CIRUN,'GASSPC',IGSOCCX,DUMMY,I_COMM_W_ND)
C
C 53: Energy convergence threshold
      IF (NZ.EQ.1) THEN
        THRES_E = 1.0D-10
      ELSE
        THRES_E = 0.5D-5
      END IF
C
C 54: Sequence of CI calculations (here one!)
      NSEQCI(1) = 1
      CALL LUCI_TRAVA(CIRUN,'MAXIT ',ISEQCI(1,1),DUMMY,I_COMM_W_ND)
C
C 62: No combinations of degenerate initial vectors
      INIDEG = 0
C
C 64: Length of batches of C vectors
C     LCSBLK = 40 000 000  
C Now calculated and set in LCISPC_REL (ciinfo_r.F)
C
C 74: Density matrices requested
      CALL LUCI_TRAVA(CIRUN,'DENSI ',IDENSI,DUMMY,I_COMM_W_ND)
C
C 81: Double group CI
      IDBGCI = 1
C
C 82: Relativistic CI
      IRELACI = 1
C
C 83: Type of operator wrt. relativity
      CALL LUCI_TRAVA(CIRUN,'IROP1 ',ISO1,DUMMY,I_COMM_W_ND)
      CALL LUCI_TRAVA(CIRUN,'IROP2 ',ISO2,DUMMY,I_COMM_W_ND)
      IF (ISO1.EQ.0.AND.ISO2.EQ.0) THEN
         ISPINFREE = 1
      ELSE
         ISPINFREE = 0
      END IF
C
C 84: Reference value for two times M_K projection
      CALL LUCI_TRAVA(CIRUN,'MK2REF',MK2REF,DUMMY,I_COMM_W_ND)
C
C 85: Reference value for two times Delta M_K projection
C     (coupling range)
      CALL LUCI_TRAVA(CIRUN,'MK2DEL',MK2DEL,DUMMY,I_COMM_W_ND)
C
C 87: Initial spinor expansion type
      ISPEXP_INI = 3
C
C 91: COMHAM option
      IF (ISEQCI(1,1) .LT. 0) THEN
C        ISEQCI(1,1) = MAXCIT from krmc (number 54 above)
C        MAXCIT .lt. 0 is krmc code for doing CI
C        with explicit CI Hamiltonian (not iterative);
C        however, ISEQCI must be .gt. 0 for luciarel to work
         IDCOMH = 1
         ISEQCI(1,1) = 1
      ELSE
C        Use iterative Davidson-Olsen type algorithm
C        with MAXCIT iterations
         IDCOMH = 0
      END IF
C
C 92: Always two-particle (full) Hamiltonian
      IHAM12 = 2
C
C 95: Wave function analysis
      CALL LUCI_TRAVA(CIRUN,'IPRWFC',IPRWFC,DUMMY,I_COMM_W_ND)
C
C***********************
C  Default settings
C***********************
C
C  6: Number of active electrons and 2-el. routine calls (sigma)
      NACTEL = IGSOCCX(NGAS,2,1)
      IELCI = 2       ! default
C
C 20: Construction of diagonal is controlled by IRESTR flag (above)
      IDIAG = 1
C
C 21: No explicit Hamiltonian
      MXP1 = 0
      MXP2 = 0
      MXQ = 0
C
C 27-28: S and L are no good quantum numbers. Don't use combinations.
      PSSIGN = 0.D0
      PLSIGN = 0.D0
      IDC = 1
C
C 40: C(Ka,Jb,j) matrices will be used in alpha-beta loop
C
      ICJKAIB = 1
C
C 44-45: Using N-2 projection method in sigma vector
      MOCAA = 0
      MOCAB = 0
C
C 52: Only 1 GAS space combination in use.
      LCMBSPC(1) = 1
      ICMBSPC(1,1) = 1
C
C 63: No use of F + Lambda(H-F) instead of H as operator
      XLAMBDA = 1.0D0
C
C 66: No MO file, transformation is performed by moltra.
      NOMOFL = 1
C
C 68: Final orbitals are computed externally (KRMC).
C     No iterative natural spinors.
      IFINMO = 0
      NATITER = 1
C
C 73: No class selection performed.
      ICLSSEL = 0
C
C 77: Restart option
C     Restarting is determined from the outside (sigma, density),
C     so set to default value
      IRST2 = 1         ! sic
C
C 78: Never skip evaluation of energy (only used in pico3, anyway)
      ISKIPEI = 0
C
C 80: No property calculations here (yet)
      NPROP = 0
C
C 86: No initial integral transformation (done with moltra).
      ITRANS = 0
C
C 93: External Fock matrix always assumed. Specification not allowed
C     in DIRAC environment.
      IFOCKEX = 0
#if defined (VAR_MPI)
C
      IF( IFOCKEX .ne. 0 ) THEN
        IF( I_COMM_W_ND .eq. 1 ) THEN
          CALL MPI_BCAST(NORB_F,NSMSH,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
      END IF
C
#endif
C
C 94: No orbitals/spinors on external file.
      ORBFILE = .false.
C
C 97: Currently always SIGDEN route (excitation class formalism)
      SIGDEN_ROUTE = 'SIGDEN'
C
C 99: multiplier for memory subtraction guess
C
      ISMEMFAC = 0
C
      CALL LUCI_TRAVA(CIRUN,'MEMFAC',ISMEMFAC,DUMMY,I_COMM_W_ND)
C
C
C 100: new definition of max. BATCH length
C
C
C     still to come ...
C
CSK?      CALL LUCI_TRAVA(CIRUN,'MAXBLK',IMAXBLKSIZE,DUMMY,I_COMM_W_ND) 
C
C     ... therefore we use
C
      IMAXBLKSIZE = 100 000 000
C
      IMAXLBLKSZ = MAX(100 000,IMAXBLKSIZE)
C
C 101: are we running a nat. orb. occ. num. run?
C
      NOOCCN_DENS_RUN = .FALSE.
      MOPT_MX = 0
      CALL LUCI_TRAVA(CIRUN,'MOPTMX',MOPT_MX,DUMMY,I_COMM_W_ND)
      IF( MOPT_MX .lt. 0 ) THEN
        NOOCCN_DENS_RUN = .TRUE.
      ELSE
        NOOCCN_DENS_RUN = .FALSE.
      END IF
C
C 102: which state to go for nat. orb. occ. num.?
C
      IDENSI_STATE = 0
      CALL LUCI_TRAVA(CIRUN,'DENS_S',IDENSI_STATE,DUMMY,I_COMM_W_ND)
C
C     That's it, dude !!
C
C     ... except we are running a parallel calculation
#if defined (VAR_MPI)
C
      IF( I_COMM_W_ND .eq. 1 ) THEN
        IF( PNTGRP .ge. 2 .and. PNTGRP .le. 4 ) THEN
          CALL MPI_BCAST(IREFML,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
        IF( PNTGRP .eq. 3 .or. PNTGRP .eq. 4 ) THEN
          CALL MPI_BCAST(IREFPA,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
        IF( PNTGRP .eq. 4 ) THEN
          CALL MPI_BCAST(IREFL,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
        IF(IDIAG.EQ.2) THEN
          CALL MPI_BCAST(E_THRE,1,MPI_REAL8,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(C_THRE,1,MPI_REAL8,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          IF(ICLSSEL.EQ.1) THEN
            CALL MPI_BCAST(E_CONV,1,MPI_REAL8,MASTER,
     &                     MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(C_CONV,1,MPI_REAL8,MASTER,
     &                     MPI_COMM_WORLD,IERR)
          END IF
        END IF
      END IF 
C
#endif
C
C***********************
C  Print section (only for testing purposes)
C***********************
C
      IF ( IPROPT .ge. 1 .and. MYPROC .eq. MASTER ) then

      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '******************'
      WRITE(LUPRI,*) '*  Title of run  *'
      WRITE(LUPRI,*) '******************'
      WRITE(LUPRI,*)
      CALL PRTITL(TITLEC)
      WRITE(LUPRI,*)
*
*. Type of calculation
*
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '************************'
      WRITE(LUPRI,*) '*  Type of calculation *'
      WRITE(LUPRI,*) '************************'
      WRITE(LUPRI,*)
      IF(IDBGCI.EQ.0.AND.IRELACI.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '     Non-relativistic calculation in spatial pointgroup'
      ELSE IF(IDBGCI.EQ.1.AND.IRELACI.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '     Non-relativistic calculation in double group '
      ELSE IF(IDBGCI.EQ.1.AND.IRELACI.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Relativistic calculation in double group '
      ELSE IF(IDBGCI.EQ.0.AND.IRELACI.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Relativistic calculation in spatial pointgroup '
        WRITE(LUPRI,'(A)')
     &  '     This is dangerous. You better think about it! '
      END IF
      IF(IRELACI.EQ.1) THEN
        IF(ISO1.EQ.1.AND.ISO2.EQ.1) then
          if (ENVIRO(1:6).ne.'DIRAC ') then
            WRITE(LUPRI,'(A)')
     &  '         One- and two- electron spin-orbit operator '
          else
            WRITE(LUPRI,'(A)')
     &  '         Dirac-Coulomb Hamiltonian        '
          end if
        end if
        IF(ISO1.EQ.1.AND.ISO2.EQ.0) then
          if (ENVIRO(1:6).ne.'DIRAC ') then
            WRITE(LUPRI,'(A)')
     &  '         One-electron spin-orbit (mean-field) operator '
          else
            WRITE(LUPRI,'(A)')
     &  '         Dirac-Coulomb Hamiltonian        '
          end if
        end if
        if (ISO1.eq.0.and.ISO2.eq.0.and.ENVIRO(1:6).eq.'DIRAC ') then
          write(LUPRI,'(A)')
     &  '         Dirac-Coulomb Hamiltonian        '
          write(LUPRI,'(A)')
     &  '         Spinfree formalism               '
        end if
        if (ISO1.EQ.0.AND.ISO2.EQ.1) then
          if (ENVIRO(1:6).eq.'DIRAC ') then
            write(LUPRI,*) 'Setting of relativistic operators not'
            write(LUPRI,*) 'allowed in DIRAC environment.'
            write(LUPRI,*) 'This is relativistic theory!'
            call Abend2 ('Input error.')
          end if
          WRITE(LUPRI,'(A)')
     &  '         Two-electron spin-orbit operator '
        end if
      END IF
      if (IFOCKEX.eq.1) then
        WRITE(LUPRI,'(A)')
     &  '         Fock matrix gen. externally for inactive orbitals '
        WRITE(LUPRI,'(A43,8I4)')
     &  '         Full orbital basis in symmetries: ',
     &            (NORB_F(ISYM),ISYM=1,NSMSH,1)
      end if
      if (INTIMP.eq.1) then
        if (ORBFILE) then
          WRITE(LUPRI,'(A)')
     &    '         MO vector read from formatted ASCII file INPORB '
      else
          WRITE(LUPRI,'(A)')
     &    '         MO vector read from RASSCF file JOBIPH          '
        end if
      end if
*
*. Type of reference state
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '*************************************'
      WRITE(LUPRI,*) '*  Symmetry and spin of CI vectors  *'
      WRITE(LUPRI,*) '*************************************'
      WRITE(LUPRI,*)
*. Point group
      IF(PNTGRP.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ D2H'
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ C inf v'
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ D inf h'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ O3'
      ELSE IF(PNTGRP.EQ.5) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... C1'
      ELSE IF(PNTGRP.EQ.6) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... Ci'
      ELSE IF(PNTGRP.EQ.7) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... Cs'
      ELSE IF(PNTGRP.EQ.8) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... C2'
      ELSE IF(PNTGRP.EQ.9) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... C2H'
      END IF
*.Spatial symmetry
      IF(PNTGRP.EQ.1) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     Spatial symmetry ....... ', IREFSM
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     ML value ............... ', IREFML
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Ungerade'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     L  value ............... ', IREFL
        WRITE(LUPRI,'(A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Ungerade'
       ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     Double group symmetry... ', IREFSM
      END IF
*.Spin
      WRITE(LUPRI,*)
      IF(IRELACI.EQ.0.AND.IDBGCI.EQ.0) THEN
C       WRITE(LUPRI,'(A,I2)')
C    &    '     2 times spinprojection  ', MS2
C       WRITE(LUPRI,'(A,I2)')
C    &    '     Spin multiplicity ....  ', MULTS
      ELSE
        WRITE(LUPRI,'(A,I2)')
     &    '     2*reference value of MK ', MK2REF
        WRITE(LUPRI,'(A,I2,A,I2)')
     &    '     Allowed interval of MK2 ',
     &          MK2REF-MK2DEL,' to ', MK2REF+MK2DEL
        WRITE(6,*)
      END IF
*.Number of active electrons
      WRITE(LUPRI,'(A,I2)')
     &  '     Active electrons .....  ', NACTEL
*
*. GAS space
*
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) ' *************************'
      WRITE(LUPRI,*) ' Generalized active space '
      WRITE(LUPRI,*) ' *************************'
      WRITE(LUPRI,*)
      IF(IRELACI.EQ.0) THEN
        WRITE(LUPRI,'(A/A)') ' Orbital subspaces:',
     &                       ' =================='
      ELSE
        WRITE(LUPRI,'(A/A)')
     &  ' Orbital subspaces in molecular pointgroup',
     &  ' ========================================='
      END IF
      WRITE(LUPRI,'(A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NSMSH)
      WRITE(LUPRI,'(A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NSMSH )
      DO IGAS = 1, NGAS
        WRITE(LUPRI,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NSMSH )
      END DO
*
      IF(IRELACI.NE.0) THEN
        WRITE(6,*)
     &' ( Information about spinors in double group comes later)'
      END IF
*
      WRITE(LUPRI,'(A,I3)')
     & ' Number of Occupation spaces : ',NCISPC
      WRITE(LUPRI,*)
      DO ICISPC = 1, NCISPC
      WRITE(LUPRI,'(A,I3)')
     &' Bounds on accumulated occupations for subspace : ',ICISPC
      WRITE(LUPRI,'(A)')
     & ' ====================================================== '
      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A)') '         Min. occ    Max. occ '
      WRITE(LUPRI,'(A)') '         ========    ======== '
      DO IGAS = 1, NGAS
        WRITE(LUPRI,'(A,I2,3X,I3,9X,I3)')
     &  '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
      END DO
      END DO
*
      WRITE(6,*)
      WRITE(6,'(A,I3)')
     &' Number of CI spaces included : ', NCMBSPC
      WRITE(6,*)
      DO JCMBSPC = 1, NCMBSPC
        WRITE(6,*)
        WRITE(6,'(A,I3)') ' Information about CI space ', JCMBSPC
        WRITE(6,'(A)')    ' =================================='
        WRITE(6,'(A,I3)')
     &  ' Number of GAS spaces included ', LCMBSPC(JCMBSPC)
        WRITE(6,'(A,10I3)')
     &  ' GAS spaces included           ',
     &    (ICMBSPC(II,JCMBSPC),II=1,LCMBSPC(JCMBSPC))
*
        WRITE(6,'(A,I3)') ' Number of calculations in this CI space ',
     &  NSEQCI(JCMBSPC)
        WRITE(6,'(A)')   '  Calculations in this subspace '
        WRITE(6,'(A)')   '  =============================='
        DO JSEQ = 1, NSEQCI(JCMBSPC)
          IF(ISEQCI(JSEQ,JCMBSPC).LT.0.AND.ISEQCI(JSEQ,JCMBSPC).GT.-5)
     &    THEN
            WRITE(6,'(A,I3)')
     &      '       Vector free calculation at level ',
     &      -ISEQCI(JSEQ,JCMBSPC)
          ELSE IF( ISEQCI(JSEQ,JCMBSPC).GT.0 ) THEN
           WRITE(6,'(A,I3)')
     &      '       Normal CI with max. iterations = ',
     &      ISEQCI(JSEQ,JCMBSPC)
          ELSE IF (ISEQCI(JSEQ,JCMBSPC).EQ.-5 ) THEN
           WRITE(6,'(A,I3)')
     &      '       Perturbation calculation '
          END IF
        END DO
      END DO
*
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '***********'
      WRITE(LUPRI,*) '*  Roots  *'
      WRITE(LUPRI,*) '*********** '
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A,I3)')
     &  '     Number of roots to be obtained ', NROOT
      WRITE(LUPRI,'(A,(20I3))')
     &  '     Roots to be obtained ', (IROOT(I),I=1, NROOT )
*
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '**************************'
      WRITE(LUPRI,*) '*  Run time definitions  *'
      WRITE(LUPRI,*) '************************** '
      WRITE(LUPRI,*)
*. Program environment
      WRITE(6,'(A,A6)') '      Program environment... ', ENVIRO
*. Integral transformation
      IF(ITRANS.NE.0) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     Initial integral transformation, level = ', ITRANS
      END IF
*. Initial spinor transformation
      IF(ISPEXP_INI.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Initial spinor expansion obtained from MO expansion'
      ELSE IF (ISPEXP_INI.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     Initial spinor expansion read in '
      else if (ISPEXP_INI.eq.3) then
        WRITE(LUPRI,'(A)')
     &  '     Reading integrals over transformed spinors. '
      END IF
*. No. of particle Hamiltonian
      if (IHAM12.eq.1) then
        WRITE(LUPRI,'(A)')
     &  '     Only 1-particle Hamiltonian under consideration '
      else if (IHAM12.eq.2) then
        WRITE(LUPRI,'(A)')
     &  '     Full 1- and 2-particle Hamiltonian active       '
      else if (IHAM12.eq.3) then
        WRITE(LUPRI,'(A)')
     &  '     Only 2-particle Hamiltonian under consideration '
      end if
*. Integral import
      IF(NOINT.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     No integrals will be read in       '
      ELSE IF(NOINT.EQ.0) THEN
      IF(INTIMP.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in MOLCAS format '
      ELSE IF(INTIMP.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in LUCAS format '
      ELSE IF(INTIMP.EQ.3) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in formatted form (E22.15) ',
     &  '      From unit 13'
        WRITE(LUPRI,'(A)')
     &  '     All integrals of correct symmetry combination read in'
      ELSE IF(INTIMP.EQ.5) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in DALTON/SIRIUS format '
      else if (INTIMP.eq.6) then
        WRITE(LUPRI,'(A)')
     &  '     Integrals in Kramers reduced spinor list format '
        WRITE(LUPRI,'(A)')
     &  '       Dirac-Coulomb Hamiltonian '
      else if (INTIMP.eq.7) then
        WRITE(LUPRI,'(A)')
     &  '     Integrals in (NZ,3) format '
        WRITE(LUPRI,'(A)')
     &  '       Dirac-Coulomb Hamiltonian '
      else if (INTIMP.eq.8) then
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in MOLECULE-SWEDEN format '
        WRITE(LUPRI,'(A)')
     &  '     one-electron ints from unit 2, transformation '
        WRITE(LUPRI,'(A)')
     &  '     spin-orbit ints from SOMFINT_*, trafo '
        WRITE(LUPRI,'(A)')
     &  '     two-electron ints from unit 13, transformed '
      END IF
*. Integral storage
      IF(INCORE.EQ.1) WRITE(LUPRI,'(A)')
     &  '     All integrals stored in core'
      END IF
* ( END IF for NOINT
*. CSF or SD expansion
      IF(NOCSF.EQ.0) THEN
        WRITE(LUPRI,'(/A/)')
     &  '     CI optimization performed with CSF''s '
      ELSE
        WRITE(LUPRI,'(/A/)')
     &  '     CI optimization performed with SD''s '
      END IF
*. Initial approximation to vectors
      IF(IRESTR.EQ.1.AND.IRESTRF.EQ.0) THEN
         WRITE(LUPRI,'(A)')
     &  '     Restarted calculation '
      ELSE IF(IRESTRF.EQ.1) THEN
         WRITE(LUPRI,'(A)')
     &  '     Restarted calculation from REFERENCE space expansion'
      ELSE
         IF(MXP1.NE.0) THEN
           WRITE(LUPRI,'(A)')
     &  '     Initial vectors obtained from explicit Hamiltonian'
         ELSE IF(MXP1.EQ.0) THEN
           WRITE(LUPRI,'(A)')
     &  '     Initial vectors obtained from diagonal'
         END IF
      END IF
*. Handling of degeneracies of initial vectors
      IF(INIDEG.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Symmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.-1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Antisymmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '     No combination of degenerate initial vectors'
      END IF
*. CI storage mode
      if (SIGDEN_ROUTE.eq.'SIGDEN') then
        WRITE(LUPRI,'(/A)')
     &  '     SIGDEN route for sigma vectors and densities invoked   '
      else
        WRITE(LUPRI,'(/A)')
     &  '     SIGDEN keyword: ',SIGDEN_ROUTE
        call quit('Unallowed keyword value')
      end if
      WRITE(LUPRI,*)
      IF(ICISTR.EQ.1) THEN
        WRITE(LUPRI,*)
     &  '     3 symmetry blocks and two vectors will be held in core '
      ELSE IF( ICISTR.EQ.2) THEN
        WRITE(LUPRI,*)
     &  '     3 symmetry blocks will be held in core '
      ELSE IF( ICISTR.EQ.3) THEN
        WRITE(LUPRI,*)
     &  '     3 symmetry-occ-occ blocks will be held in core '
      END IF
      IF(LCSBLK.NE.0) WRITE(LUPRI,'(A,I10)')
     &  '      Smallest allowed size of sigma- and C-batch ',LCSBLK
      WRITE(LUPRI,'(A,I4)')
     &  '     Dimension of block of resolution strings ', MXINKA
      IF(ICJKAIB.EQ.1) THEN
        WRITE(LUPRI,*)
     &  '     C(Ka,Jb,j) matrices will be used in alpha-beta loop '
      END IF
      IF(MOCAA.NE.0) WRITE(LUPRI,'(A,I4)')
     &  '     MOC method used for alpha-alpha+beta-beta loop '
      IF(MOCAB.NE.0) WRITE(LUPRI,'(A,I4)')
     &  '     MOC method used for alpha-beta loop            '
*
      WRITE(LUPRI,*)
      IF(IDENSI.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '     No calculation of density matrices  '
      ELSE IF(IDENSI.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     One-body density matrix calculated '
      ELSE IF(IDENSI.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     One- and two-body density matrices  calculated '
      else if (IDENSI.eq.3) then
        WRITE(LUPRI,'(A)')
     &  '     Spin-dependent one-body density matrix calculated '
        WRITE(LUPRI,'(A)')
     &  '     Spin-dependent one- and two-body densities
     &        calculated '
      END IF
*
*
*. Diagonalization information
      WRITE(LUPRI,'(A/A)')
     &  '     CI diagonalization : ',
     &  '     ==================== '
*. Subspace Hamiltonian
      IF (MXP1+MXP2+MXQ .EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '        No subspace Hamiltonian '
      ELSE
        WRITE(LUPRI,'(A,3I4)')
     &  '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
      END IF
*. Diagonalizer
      IF(IDIAG.EQ.1.AND.ICISTR.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &    '        Diagonalizer : MINDV4 '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &    '        Diagonalizer : MICDV4 '
      ELSE IF(IDIAG.EQ.2) THEN
      WRITE(LUPRI,'(A)')
     &  '        Diagonalizer : PICO2  '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.3) THEN
        IF (NZ.EQ.1) THEN
          WRITE(LUPRI,'(A)')
     &      '        Diagonalizer : MICDV6 '
        ELSE
          WRITE(LUPRI,'(A)')
     &      '        Diagonalizer : CMICDV '
        END IF
      END IF
*
*. No restart in CI calc 2
      IF(IRST2.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '        No restart from previous vectors in second calc '
      END IF
      IF(ISKIPEI.EQ.1) THEN
        WRITE(LUPRI,'(A/A)')
     &  '        Initial energy evaluations skipped after first calc',
     &  '        (Only active in connection with TERACI )'
      END IF
*
*. Number of iterations
      WRITE(LUPRI,'(A,I2)')
     &  '        Allowed number of iterations    ',MAXIT
*. Number of CI vectors in subspace
      WRITE(LUPRI,'(A,I2)')
     &  '        Allowed Dimension of CI subspace ',MXCIV
*
      WRITE(LUPRI,'(A,1P,E11.5/)')
     &  '        Convergence threshold for energy ',THRES_E
      IF(IDIAG.EQ.2) THEN
        WRITE(LUPRI,'(A,1P,E11.5)')
     &   '        Individual second order energy threshold ',E_THRE
        WRITE(LUPRI,'(A,1P,E11.5)')
     &   '        Individual first order wavefunction threshold ',C_THRE
        IF(ICLSSEL.EQ.1) THEN
         WRITE(LUPRI,'(/A/A)')
     &   '         Class selection will be performed : ',
     &   '         =================================== '
         WRITE(LUPRI,'(A,1P,E11.5)')
     &    '          Total second order energy threshold ',E_CONV
         WRITE(LUPRI,'(A,1P,E11.5)')
     &    '          Total first order wavefunction threshold ',C_CONV
        ELSE
         WRITE(LUPRI,'(A)')
     &'            No class selection in iterative procedure '
        END IF
      END IF
*
*
*. Print levels
*
      WRITE(LUPRI,'(/A)')  '     Print levels : '
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for string    information = ', IPRSTR
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for wave function info    = ', IPRWFC
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for CI space  information = ', IPRCIX
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for orbital   information = ', IPRORB
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for density matrix        = ', IPRDEN
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for iterative information = ', IPRDIA
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for External blocks       = ', IPRXT
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for transformation info   = ', IPRTRA
      WRITE(LUPRI,'(A,I3)')
     &'      Chosen  print level for sigma vector info     = ', IPRSIG
*
      IF(IPROCC.NE.0) THEN
        WRITE(LUPRI,'(A,I3)')
     &  '      Occupation of SD''s/ configurations           = ',IPROCC
      END IF
C
      END IF
C     ^ end if IPROPT for input printing.
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_TRAVA(CIRUN,CTLSTR,IVAR,DVAR,I_COMM_W_ND)
C***********************************************************************
C
C     Translate a variable given as CTLSTR
C     to output value IVAR
C     MASTER communicates with NODES w.r.t. I_COMM_W_ND ( 0 or 1 )
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision : Mar 14, 2001  T. Fleig (bugfix)
C                     May 17, 2002  T. Fleig (ROOTCI added.)
C                     May 27, 2002  T. Fleig (MXCIV added.)
C                     Apr 03, 2003  T. Fleig (IRESTRK added.)
C                     Apr 17, 2007  S. Knecht (Parallel adaption)
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "../luciarel/lucipar_rel.inc"
C
#include "../luciarel/mxpdim.inc"
C
#include "mx2orb.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbham.h"
#include "dcbbos.h"
C
      CHARACTER*6 CTLSTR,CIRUN
      DIMENSION IVAR(*),DVAR(*)
C
C     Set number of CI roots (for start guess only)
C
      IF (CTLSTR.EQ.'ROOTCI') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(NCIROOT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IVAR(1) = NCIROOT
         GOTO 100
      END IF
C
C     Set max. CI subspace dimension
C
      IF (CTLSTR.EQ.'MXCIV ') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(MXCIV,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IF (MXCIV.EQ.0) THEN
            IVAR(1) = 3*NCIROOT
            IF( MYPROC .eq. MASTER .AND. I_COMM_W_ND .eq. 0 ) THEN 
              write(LUPRI,'(A,I4)') '   Setting MXCIV to ',3*NCIROOT
            ELSE IF( MYPROC .ne. MASTER ) THEN
              write(MSLVOUT_REL,'(A,I4)') 
     &'   Setting MXCIV to ',3*NCIROOT
            END IF
         ELSE
            IVAR(1) = MXCIV
         END IF
         GOTO 100
      END IF
C
C     Set CI reference symmetry
C
      IF (CTLSTR.EQ.'IREFSM') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(IOPT_SYMMETRY,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IVAR(1) = IOPT_SYMMETRY
         GOTO 100
      END IF
C
C     Set CI restart variable
C
      IF (CTLSTR .EQ. 'IRESTR') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(IRESTRK,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
C
C        LUCIAREL  start on CI vectors?
C
C        For the initial CI iterations LUCIAREL shouldn't
C        start from vectors on file (except if IRESTRK = 1);
C        for sigma-vectors or density matrices LUCIAREL should
C        (re)start from vectors on file. In the density matrix case,
C        IRESTR = 0 and the CI vectors are placed on the files
C        required by LUCIAREL.
C
         IF (CIRUN .EQ. 'SIGMA ' .OR. CIRUN .EQ. 'SIGMAD'
     &       .OR. (CIRUN .EQ. 'CIINII' .AND. IRESTRK .EQ. 1) )
     &       THEN
            IVAR(1) = 1
         ELSE
            IVAR(1) = 0
         END IF
         GOTO 100
      END IF
C
C     Set wave function printout
C
      IF (CTLSTR .EQ. 'IPRWFC') THEN
         IF (CIRUN .EQ. 'ANALYZ') THEN
            IVAR(1) = 1
         ELSE
            IVAR(1) = 0
         END IF
         GOTO 100
      END IF
C
C     Core energy
C
      IF (CTLSTR .EQ. 'ECORE ') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(EMY,1,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
         DVAR(1) = EMY
         GOTO 100
      END IF
C
C     #GAS shells.
C
      IF (CTLSTR .EQ. 'NGAS  ') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(NGAS,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IVAR(1) = NGAS
C
C        Also check max dimension same in krmc and luciarel
C
         IF (MXPNGAS .NE. MXMCGAS) THEN
            WRITE(LUPRI,*) '*** ERROR *** ' //
     &                'MXPNGAS (LUCIAREL) .ne. MXMCGAS (DIRAC):',
     &                 MXPNGAS,MXMCGAS
            CALL QUIT('*** ERROR *** ' //
     &                'MXPNGAS (LUCIAREL) .ne. MXMCGAS (DIRAC)')
         END IF
C
         GOTO 100
      END IF
C
C
C     GAS shells.
C
      IF (CTLSTR .EQ. 'GASSH ') THEN
C
#if defined (VAR_MPI)
        IF( I_COMM_W_ND .eq. 1 ) THEN
          CALL MPI_BCAST(OPT_UCIBOS,1,MPI_LOGICAL,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(NBSYM,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(NFSYM,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
C?          CALL MPI_BCAST(NGSH,NFSYM*NGAS,my_MPI_INTEGER,MASTER,
          CALL MPI_BCAST(NGSH,2*MXPNGAS,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(IORB,NFSYM,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(NPSH,NFSYM,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(NISH,NFSYM,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(NASH,NFSYM,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
C        FIXME: MAXORB might be to much, find min(bcast length) !!! sk
          CALL MPI_BCAST(IBOSYM,MAXORB,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
#endif
C
         IF (IPROPT .GE. 3) THEN
           IF( I_COMM_W_ND .eq. 1 .and. MYPROC .le. MASTER + 1 ) THEN 
             WRITE(LUPRI,'(/,1X,A,I4)') '(LUCI_TRAVA) Dirac 
     & GAS shells: ',MYPROC
            DO I = 1, NGAS
               WRITE(LUPRI,'(1X,A,I3,4X,2I4)')
     &              'GAS shell ',I,(NGSH(J,I),J=1,NFSYM)
            END DO
            IF (OPT_UCIBOS) THEN
               WRITE(LUPRI,'(/,1X,A)')
     &           '(LUCI_TRAVA) "Boson" symmetries of active space'
               DO I = 1, NFSYM
                  WRITE(LUPRI,'(4X,A,I1)') 'Fermion irrep ',I
                  WRITE(LUPRI,'(4X,20I3)')
     &              (IBOSYM(IORB(I)+NPSH(I)+NISH(I)+J),J=1,NASH(I))
               END DO
            END IF
           END IF
         END IF
C
C        Erase IVAR first
C
         DO IER=1,NGAS*MXPIRR,1
           IVAR(IER) = 0
         END DO
C
         DO IFSYM = 1, NFSYM
            IOFF = 0
            IF (.NOT.OPT_UCIBOS) THEN
C              To get correct Fermion symmetry labeling in
C              LUCIAREL we need to specify same(!) gerade boson
C              irrep for all gerade orbitals and same ungerade
C              boson irrep for all ungerade orbitals (note: which
C              gerade or ungerade boson irrep we choose doesn't
C              matter).
               IF (IFSYM.EQ.1) THEN
C              boson symmetry 0 is always gerade
                  IBOSYMI = 0
               ELSE
C              find first ungerade irrep
C    FIXME: This is hardwired now to be irrep number 3.
C           Should be made flexible later on to account for
C           different symmetry generator setup.
            if (CIRUN.eq.'CIINII') then
             if( myproc .eq. master .and. i_comm_w_nd .eq. 0) then
               write(lupri,*) 'WARNING FIXME: IBOSYM hardwired 
     & to irrep 3'
               write(lupri,*) 'You have problems if irrep 3 
     & not ungerade!!'
             end if
            end if
                  IBOSYMI = IFSYM
               END IF
               IBOSYM(IORB(IFSYM)+NPSH(IFSYM)+NISH(IFSYM)+IOFF)=IBOSYMI
            END IF
            DO IGAS = 1, NGAS
               DO I = 1, NGSH(IFSYM,IGAS)
                  IOFF = IOFF + 1
                  IF (OPT_UCIBOS) THEN
                     IBOSYMI =
     &                 IBOSYM(IORB(IFSYM)+NPSH(IFSYM)+NISH(IFSYM)+IOFF)
                  END IF
                  IADR = IBOSYMI + 1 + (IGAS-1) * MXPIRR
                  IVAR(IADR) = IVAR(IADR) + 1
C
               END DO
            END DO
         END DO
         IF (IPROPT .GE. 3) THEN
            WRITE(LUPRI,'(/A)') ' (LUCI_TRAVA) LUCIAREL GAS shells: '
            WRITE(LUPRI,'(10X,8(1X,A4,I2))') ('GAS ',I,I=1,NGAS)
            CALL IWRTMA(IVAR,NBSYM,NGAS,MXPIRR,MXPNGAS)
         END IF
C
         GOTO 100
      END IF
C
C     GAS space constraints.
C
      IF (CTLSTR .EQ. 'GASSPC') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
C?     &   CALL MPI_BCAST(NGSP,2*NGAS,my_MPI_INTEGER,MASTER,
     &   CALL MPI_BCAST(NGSP,2*MXPNGAS,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
C
C        Put in the GAS spaces:
C
         DO I = 1, NGAS
            IVAR(I) = NGSP(1,I)
            IVAR(I + MXPNGAS) = NGSP(2,I)
         END DO
         GOTO 100
C
      END IF
C
C     Number of CI iterations
C
      IF (CTLSTR .EQ. 'MAXIT ') THEN
         IF ( CIRUN .eq. 'SIGMA '.or. CIRUN .eq. 'SIGMAD' ) THEN
            IVAR(1) = 1
         ELSE IF (CIRUN .EQ. 'CIINII') THEN
C
#if defined (VAR_MPI)
           IF( I_COMM_W_ND .eq. 1 ) 
     &       CALL MPI_BCAST(MAXCIT,1,my_MPI_INTEGER,MASTER,
     &                      MPI_COMM_WORLD,IERR)
#endif
            IVAR(1) = MAXCIT
         ELSE
            IVAR(1) = 1
C           At least 1, otherwise driver GASCI is never reached!
         END IF
         GOTO 100
      END IF
C
C     Density matrices
C
      IF (CTLSTR .EQ. 'DENSI ') THEN
         IF (CIRUN .EQ. 'DENS1 ') THEN
            IVAR(1) = 3
         ELSE IF (CIRUN .EQ. 'DENS2 ') THEN
            IVAR(1) = 4
         ELSE IF (CIRUN .EQ. 'DENS2T') THEN
            IVAR(1) = 5
            call quit('write me')
         ELSE
            IVAR(1) = 0
         END IF
         GOTO 100
      END IF
C
C     Relativistic 1-el. operator?
C
      IF (CTLSTR .EQ. 'IROP1 ') THEN
C
#if defined (VAR_MPI)
        IF( I_COMM_W_ND .eq. 1 ) THEN 
          CALL MPI_BCAST(SPINFR,1,MPI_LOGICAL,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(LEVYLE,1,MPI_LOGICAL,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
#endif
         IF (SPINFR.OR.LEVYLE) THEN
            IVAR(1) = 0
         ELSE
            IVAR(1) = 1
         END IF
         GOTO 100
      END IF
C
C     Relativistic 2-el. operator?
C
      IF (CTLSTR .EQ. 'IROP2 ') THEN
         IF (SPINFR.OR.LEVYLE) THEN
            IVAR(1) = 0
         ELSE
            IVAR(1) = 1
         END IF
         GOTO 100
      END IF
C
C     Reference value for 2*M_K
C
      IF (CTLSTR .EQ. 'MK2REF') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(MK2REF,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IVAR(1) = MK2REF
         GOTO 100
      END IF
C
C     Reference value for 2*DeltaM_K
C
      IF (CTLSTR .EQ. 'MK2DEL') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(MK2DEL,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IF ((.not.SPINFR).and.(.not.LEVYLE).and.MK2DEL.eq.0) THEN
           WRITE(6,*) 'Spin-dependent calculation and MK2DEL = 0.'
           WRITE(6,*) 'This is not allowed.'
           STOP 'Quitting.'
         END IF
         IVAR(1) = MK2DEL
         GOTO 100
      END IF
C
C     No integrals
C     If NOINT = 1 then LUCIAREL will return just after the number
C     of determinants have been calculated.
C
      IF (CTLSTR .EQ. 'NOINT ') THEN
         IF (CIRUN .EQ. 'NDET  ') THEN
            IVAR(1) = 1
         ELSE
            IVAR(1) = 0
         END IF
         GOTO 100
      END IF
C
C     Print flag
C
      IF (CTLSTR .EQ. 'IPROPT') THEN
         IVAR(1) = IPROPT
         GOTO 100
      END IF
C
C     memory multiplier
C
      IF (CTLSTR .EQ. 'MEMFAC') THEN
C
#if defined (VAR_MPI)
       IF( I_COMM_W_ND .eq. 1 ) 
     &   CALL MPI_BCAST(IMEMFAC,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
         IVAR(1) = IMEMFAC
         GOTO 100
      END IF
CSKC
CSKC     I/O mode (parallel case)
CSKC
CSK      IF (CTLSTR .EQ. 'PAR IO') THEN
CSK#if defined (VAR_MPI)
CSK       IF( I_COMM_W_ND .eq. 1 )
CSK     &   CALL MPI_BCAST(IDOPARIO,1,my_MPI_INTEGER,MASTER,
CSK     &                  MPI_COMM_WORLD,IERR)
CSK#endif
CSK        IVAR(1) = IDOPARIO
CSK        GOTO 100
CSK      END IF
C
C     Nat. orb. occ. num. run
C
      IF (CTLSTR .eq. 'MOPTMX') THEN
C
        IF( CIRUN .eq. 'DENS1 ' .or. CIRUN .eq. 'DENS2 ') THEN
C
#if defined (VAR_MPI) 
         IF( I_COMM_W_ND .eq. 1 )
     &     CALL MPI_BCAST(NATOLCR,1,MPI_LOGICAL,MASTER,
     &                    MPI_COMM_WORLD,IERR)
#endif
           IF( NATOLCR) THEN
             IVAR(1) = - 1
           ELSE
             IVAR(1) =   1
           END IF
           GOTO 100
        ELSE
          GOTO 100
        END IF
      END IF
C
C     Nat. orb. occ. num. state
C
      IF (CTLSTR .eq. 'DENS_S') THEN
C
        IF( CIRUN .eq. 'DENS1 ' .or. CIRUN .eq. 'DENS2 ') THEN
C
#if defined (VAR_MPI) 
         IF( I_COMM_W_ND .eq. 1 )
     &     CALL MPI_BCAST(IDENSLR_STATE,1,my_MPI_INTEGER,MASTER,
     &                    MPI_COMM_WORLD,IERR)
#endif
           IVAR(1) = IDENSLR_STATE
           GOTO 100
        ELSE
          GOTO 100
        END IF
      END IF
C
      WRITE(LUPRI,'(1X,2A)')
     &   '*** ERROR in LUCI_TRAVA *** Illegal keyword : ',CTLSTR
      CALL QUIT('*** ERROR in LUCI_TRAVA ***')
C
 100  CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_SETUP()
C***********************************************************************
C
C     Set up starting information for LUCIAREL
C      - assign disk units
C      - determinant types
C      - string types
C      - symmetry information
C      - integral dimensioning
C
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "../luciarel/cprnt.inc"
C
C     Assign disk units
      CALL DISKUN_REL
C
C. Ratio beteeen real and integer word length
      CALL ZIRAT_REL
C. From shells to orbitals
      CALL ORBINF_REL(LUPRI,IPRORB)
C. Types of determinants to be included
      CALL DETTYP(IPRSTR)
C. Number of string types
      CALL STRTYP_GAS_REL(IPRSTR)
C. Symmetry information
      CALL SYMINF_REL(IPRORB)
C. Number of integrals
      CALL INTDIM_REL(IPRTRA)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRMC_LUCI(F1,F2,WRK,LWRK,ECORE_KRMC,CIRUN)
C***********************************************************************
*
*
* L U C I A R E L
*
* parallel driver routine,
* Stefan Knecht - April 2007
*
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "../luciarel/lucipar_rel.inc"
      CHARACTER*6 CIRUN
C     initial definitions
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
*     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
C
C     Summon the slaves, who are waiting in the general menu routine.
C
      CALL LUMSTI_REL
C
C     communicate necessary information with slaves ...
C
#if defined (VAR_MPI)
C
      CALL SYNC_NODES_REL_P1(ECORE_KRMC,CIRUN,LUPRI)
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or.
     &    CIRUN .eq. 'SIGMAD') THEN
C
        CALL SYNC_NODES_REL_P2(F1,F2,CIRUN)
C
      END IF
C
C     Transfer data to LUCIAREL common blocks (only for NODES) 
C     --> communication is enabled.
C
      CALL LUCI_TRACOM(CIRUN,1)
C
#endif
C
C     Transfer correct output file handle to common block
C
      LUWRT = LUPRI
C
C     Enter the generic LUCIAREL routine 
C     (master becomes now somewhat like primus inter parens)
C
      CALL KRMC_LUCI2(F1,F2,WRK,LWRK,ECORE_KRMC,CIRUN)
C
C     close files ...
C
CSK?      CALL LUCI_CLOSEDOWN()
C
C     Return the slaves to the general menu routine.
C
      CALL LUMSTE_REL
C
C
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KRMC_LUCI2(F1,F2,WRK,LWRK,ECORE_KRMC,CIRUN)
C***********************************************************************
*
*
* L U C I A R E L
*
*
* CI for program for :FCI
*                     RASCI
*                     GAS GAS GAS GAS GAS GAS
*
* Written by Jeppe Olsen , winter of 1991
*                          GAS version in action summer of 95
*
* Implementation of relativistic double group CI,
*                          Timo Fleig + Jeppe Olsen
*                          1995 - 1999
*
* parallel adaption, Stefan Knecht, April 2007 - ?
*                  
*
* Simplified and modified from luciarel.f
* Execute task given by task parameter CIRUN
*
C***********************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Parameters for dimensioning
#include "../luciarel/mxpdim.inc"
*. Boson symmetry info
#include "mx2orb.h"
#include "dcbbos.h"
*.Memory
#include "../luciarel/wrkspc.inc"
#include "../luciarel/ipoist8.inc"
*.File numbers
#include "../luciarel/clunit.inc"
*.Print flags
#include "../luciarel/cprnt.inc"
#include "../luciarel/lucinp.inc"
#include "../luciarel/cstate.inc"
#include "../luciarel/crun.inc"
#include "../luciarel/cicisp.inc"
#include "../luciarel/oper.inc"
#include "../luciarel/cgas.inc"
#include "../luciarel/glbbas.inc"
#include "../luciarel/csm.inc"
#include "../luciarel/cintfo.inc"
#include "../luciarel/picasso.inc"
CTF
      character*6 CIRUN
      CHARACTER*72 MACHINENAME
      INTEGER NAMELENGTH
      logical MCSCF
CTF  global reference values for Z_BLKFO needed in MV7 !!
      common/glbref/ISPC,ISM,ITERGLB
      DIMENSION F1(*),F2(*),WRK(LWRK)
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "../luciarel/lucipar_rel.inc"
*
*  FIXME: VCT1 and VCT2 unused at the moment. They should carry
*         input/output vectors controlled by outside KR_LUCI_IO
*         routine. At the moment, allocation is performed inside
*         LUCIAREL. The vectors reside on the internal files
*         fort.20 (D), fort.21 (C), and fort.31 (S).
*         Renamed files. See diskun.f
*
********************************************
      CALL QENTER('LUCIAR')
C
C
C??      WRITE(LUWRT,*) 'TEST from KRMC_LUCI2',MYPROC
C
C
C     Compute offset from WRK to WORK
C
      CALL COMPOFF(WRK,WORK,KBASE)
C 32 bit:  *kbase = (int) (wrk-work+1);
C 64 bit:  *kbase = (long) (wrk-work+1);
C
      WORK(KBASE) = 9991.
      IF ((WRK(1).NE.9991.).AND.(CIRUN.EQ.'CIINII')) THEN
        write(LUWRT,*) ' WARNING from KRMC_LUCI:'
        WRITE(LUWRT,'(A,3X,I22)')     '      KBASE  = ',KBASE
        WRITE(LUWRT,'(A,3X,1F22.16)') ' WORK(KBASE) = ',WORK(KBASE)
        WRITE(LUWRT,'(A,3X,1F22.16)') ' WRK(1)      = ',WRK(1)
        write(LUWRT,*) ' Absolute memory addresses do not match.'
        write(LUWRT,*) ' Continuing.'
c       CALL QUIT('Memory addressing error in KRMCLUCI')
      END IF
*
      MXPWRD = LWRK
      KADD = MXPWRD
      call memmar(KBASE,KADD,'INI   ',IDUMMY,'DUMMY')
*
      NTEST = 0
      IMEMCK = 0
      MCSCF = .true.
*
*     debugging print for parallel runs -- control from outside?
*
      NPTEST_VAR = 0
*
*.    Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_REL(WORK(1),IPRSTR)
*
*.    Internal subspaces
      CALL LCISPC_REL(LCSBLK,IPRCIX)
*
******************************************
      if (CIRUN.eq.'RSTRMC') then
        ISPC = 1
        CALL Z_BLKFO_REL(IDC,NMS2VAL,ISPC,IREFSM,KPCLBT,KPCLEBT,
     &                   KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
     &                   NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2,2,0,
     &                   IDUMMY)
        CALL QEXIT('LUCIAR')
        RETURN
      end if
******************************************
      if (IMEMCK.gt.0) then
        write(LUWRT,*) ' Memory map after STRINF_GAS_REL'
        IDUM = 0
        call memmar(IDUM,IDUM,'PRINT ',IDUM,'IDUMMY')
      end if
*
*******************************************
      if (CIRUN.ne.'NDET  ') then
*******************************************
*
*.      setting up groups of processors (stored in /LUCIPARREL/)
*
        CALL MEMMAR(IDUM,IDUM,'MARK ',IDUM,'GROUPL')
*
        CALL MEMMAR(KPROCLIST,NMPROC,'ADDS  ',1,'IHOSTL')
        CALL MEMMAR(KGROUPLIST,NMPROC,'ADDS  ',1,'IHOST2')
        CALL MEMMAR(KPROCLENGTH,NMPROC,'ADDS  ',1,'INAMEL')
*
        CALL ISETVC(WORK(KPROCLIST),-1,NMPROC)
        CALL ISETVC(WORK(KGROUPLIST),-1,NMPROC)
        CALL ISETVC(WORK(KPROCLENGTH),-1,NMPROC)
*       initialize number of process groups sharing a c-vector file
        NFLGRPS_REL = 0
*
#if defined (VAR_MPI)
*
*       hardwired to 1 (if at least two cpu's active, i.e. numnod .gt. 0)
*
        IF (NUMNOD .GT. 0) THEN
           IPARALLELIO = 1
        ELSE
           IPARALLELIO = -1
        END IF
*
#if defined (SYS_LOCALDISKS)
*       checking only if parallel file I/O in general is requested
        IF( IPARALLELIO .eq. 1 ) THEN
*  
            NAMELENGTH = 0
*
            CALL MPI_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH,IERR)
*
            CALL FIND_GROUP_OF_PROCS_REL(MACHINENAME,NAMELENGTH,
     &                               WORK(KPROCLIST),WORK(KPROCLENGTH))
*
*           number of groups sharing a file ( = NFLGRPS ) is determined
*           in the subroutine above. If there is only one group
*           than we can simulate a shared filesystem.
      END IF
#else
*       in case of a parallel filesystem, we should use it
        IF( IPARALLELIO .ge. 0 ) THEN
          NFLGRPS = 1
          CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
*         enable parallel I/O mode!!!
          IPARALLELIO = 1
        END IF
#endif /* ifdef SYS_LOCALDISKS */
*
*       set-up all useful MPI stuff: communication groups, MPI-files ...
*                                    all important variables are stored on 
*                                    common block /LUPARGROUP/
*
*       construct new communication groups
        CALL GROUP_CONSTRUCTOR_REL(WORK(KGROUPLIST),WORK(KPROCLIST))
*
        IF( NPTESTVAR .ge. 20 ) THEN
          WRITE(LUWRT,*) 'this is my old and new rank ',myproc,MYNEW_ID
          WRITE(LUWRT,*) 'this is my grouplist:'
          CALL IWRTMAMN(WORK(KGROUPLIST),1,NMPROC,1,NMPROC,LUWRT)
          WRITE(LUWRT,*) 'this is my proclist:'
          CALL IWRTMAMN(WORK(KPROCLIST),1,NMPROC,1,NMPROC,LUWRT)
        END IF
*
*       open the files for parallel file I/O
*
        CALL SETUNITS_PAR_OPEN_REL(MY_GROUPN)
*
#else
        CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
        NFLGRPS_REL = 0
#endif
*
*. Allocate work space, general
*   (position changed as not required for TASK = NDET)
        CALL ALLOC_REL
*
*
* IDGSTYP determines the integral classes and corresponding operator
*         classes depending on the double group in use.
        IDGSTYP = -1
        if (PNTGRP.le.4) then
          IDGSTYP = 1
        else if (PNTGRP.eq.6.or.PNTGRP.eq.5) then
          IDGSTYP = 3
        else if (PNTGRP.ge.7.and.PNTGRP.le.9) then
          IDGSTYP = 2
        end if
        if (SIGDEN_ROUTE.eq.'SIGDEN') then
*. Space for Hamiltonian in form suitable for SIGDEN
          if (ISO1.eq.0.and.ISO2.eq.0) then
*    Non-relativistic (or spinfree) Hamiltonian
            IHTYPE = 1
          else if (ISO1.eq.1.and.ISO2.eq.0) then
*    Mean-field SO Hamiltonian
            IHTYPE = 2
          else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).ne.'DIRAC ')
*    Relativistic 1- and 2-particle Hamiltonian
     &         then
            IHTYPE = 3
          else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).eq.'DIRAC ')
*    Dirac-Coulomb case.
     &         then
            if (IDGSTYP.le.2) then
              IHTYPE = 4
            else if (IDGSTYP.eq.3) then
              IHTYPE = 5
            else
              write(LUWRT,*) 'Unallowed value of IDGSTYP : ',IDGSTYP
              stop 'Quitting.'
            end if
          else
            write(LUWRT,*) 'Setting of ISO1,ISO2 not allowed.'
            write(LUWRT,*) 'ISO1,ISO2 : ',ISO1,ISO2
            call abend2('Quitting.')
          end if
          ISYM_T = 1
          call set_hop_dbg(IHTYPE,IPRHAM)
        end if
*       ^ End if SIGDEN_ROUTE
*
*. Symmetry of reference
        IF(PNTGRP.GT.1.and.PNTGRP.lt.5)
     &             CALL MLSM_REL(IREFSM,IREFPA,IREFSM,'CI',1)
        IF(NOINT.EQ.1) THEN
          WRITE(LUWRT,*) ' End of calculation without integrals'
          CALL QEXIT('LUCIAR')
          RETURN
        END IF
*
        do INATIT=1,NATITER,1
*
*  Generate symmetry info and the like
          if (INATIT.eq.1) call luc2pic(ISO1,IPRTRA)
*
*. Generate and store pointers for transformed integrals
*  (and integrals !)
          IF(NOINT.EQ.0.and.INATIT.eq.1) THEN
            CALL INTIM_REL(IPRTRA)
          ELSE
            WRITE(LUWRT,*)
     &        ' No integrals imported (or > 1st NS iteration) '
          END IF
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          if (CIRUN.ne.'ANALYZ') then
            call picasso(INTSPACE,ISO1,IPRTRA,MCSCF,
     &                   NL2D,IBU2D,ISM2D,F1,F2,IBOSYM,CIRUN)
          end if
*+++++++++++++++++++++++++++++++++++++++++++++++++++
*
*
* ======================
* CI in GAS spaces
* ======================
*
          ISTOSPC = 0
          IF(IRESTR.EQ.1) ISTOSPC = 1
*
*. Restarted Lambda calculations need special attention
*  since restart is realized in in second calculation
*  First calculation is used to establish H0
*
          IRESTR_ORIG=IRESTR
          IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0) THEN
*
            write(LUWRT,*) 'THIS ROUTE IS STILL FORBIDDEN !!'
            write(LUWRT,'(A16,I2,1F3.1)')
     &            'IRESTR,XLAMBDA =',IRESTR,XLAMBDA
            Call Abend2( 'Quitting.' )
*
          END IF
*         ^ End of special handling of restarted calc with
*           lambda modified op.
*
          ISKIPEI_INI = ISKIPEI
          DO JCMBSPC = 1, NCMBSPC
            ISKIPEI = 0
*
*. Perform GASCI/GASPT calculation for internal space ICISPC
*
            if (NTEST.ge.1) then
              WRITE(LUWRT,*)
              WRITE(LUWRT,'(A)')
     &        ' =============================================='
              WRITE(LUWRT,'(A,I3)')
     &        ' Information on calculations in space ', JCMBSPC
              WRITE(LUWRT,'(A)')
     &        ' =============================================='
              WRITE(LUWRT,*)
              WRITE(LUWRT,'(A,I3)')
     &        ' Number of calculation in this CI space ',
     &          NSEQCI(JCMBSPC)
            end if
*
            IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0.AND.JCMBSPC.EQ.1) THEN
              WRITE(LUWRT,*)
     &              ' Remember No restart in calc 1 (Lambda calc)'
              IRESTR = 0
              Call Abend2( 'Quitting.' )
            END IF
*
            I_EXPAND = 1
            IF(XLAMBDA.NE.1.0D0 .AND.JCMBSPC.GT.1) THEN
              WRITE(LUWRT,*) ' =================================='
              WRITE(LUWRT,*) '   Modified operator will be used'
              WRITE(LUWRT,*) ' =================================='
              IF(JCMBSPC.EQ.2) THEN
*. Modify the operator
                WRITE(LUWRT,*) ' Operator will be modified '
                WRITE(LUWRT,*) ' This route is still forbidden! '
                Call Abend2( 'Quitting.' )
              END IF
            END IF
            IF(JCMBSPC.EQ.2.AND.IRESTR_ORIG.EQ.1
     &                     .AND.XLAMBDA.NE.1.0D0)  THEN
              WRITE(LUWRT,*) ' Restart vectors will be copied to LUC'
              WRITE(LUWRT,*)
     &              ' CI will restart with vectors from prev. calc'
              IF(ICISTR.EQ.1) THEN
                LBLK = XISPSM(IREFSM,2)
              ELSE
                LBLK = -1
              END IF
              call memmar(IDUM,IDUM,'MARK ',IDUM,'VCSAVE')
              IF(ICISTR.EQ.1) THEN
                LBLOCK = XISPSM(IREFSM,2)
              ELSE IF (ICISTR.EQ.2) THEN
                LBLOCK = MXSB
              ELSE IF (ICISTR.EQ.3) THEN
                LBLOCK = MXSOOB
              END IF
              call memmar(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
              CALL REWINO(LUC)
              CALL REWINO(LUSC39)
              DO JROOT = 1, NROOT
                CALL COPVCD_REL(LUSC39,LUC,WORK(KVEC1),0,LBLK)
              END DO
              call memmar(IDUM,IDUM,'FLUSM',IDUM,'VCSAVE')
*. No expansion should follow
              I_EXPAND = 0
*. But normal restart
              IRESTR = 1
            END IF
*
            DO JSEQ = 1,NSEQCI(JCMBSPC)
              IF((ISEQCI(JSEQ,JCMBSPC).GT.0
     &            .OR.ISEQCI(JSEQ,JCMBSPC).EQ.-5)
     &            .AND. JCMBSPC.NE.1.AND.ISTOSPC.NE.JCMBSPC
     &            .AND.I_EXPAND.EQ.1)                   THEN
*. Restart from previous spaces ( Assuming a progressing sequence :
*  spaces are just added, not subtracted )
*  ( Used for perturbation calculations and CI calculations )
*
                LUIN = LUC
                LUOUT = LUSC1
                IF(ICISTR.EQ.1) THEN
                  LBLK = XISPSM(IREFSM,JCMBSPC)
                ELSE
                  LBLK = -1
                END IF
                WRITE(LUWRT,*) ' LBLK = ', LBLK
                CALL EXPCIV_REL(IREFSM,ISTOSPC,LUIN,JCMBSPC,LUOUT,
     &                          LBLK,LUSC2,
     &                          NROOT,1,IDC,IPRDIA)
*. Last space where vectors were stored
                ISTOSPC = JCMBSPC
                ISKIPEI = ISKIPEI_INI
*. Expanded vector will be used as initial vector in the
*. zero space calculation. Tell next CI to restart from
*. CI vectors
                IRESTR = 1
              END IF
              IF(ISEQCI(JSEQ,JCMBSPC).GT.0) THEN
                IF(JSEQ.EQ.1.AND.JCMBSPC.EQ.2.AND.IRST2.EQ.0) THEN
*. No restart from previous vectors - IRST2 has been set to zero
                  IRESTR = 0
                  WRITE(LUWRT,*) ' No restart from previous vectors'
                END IF
*. Good old normal CI !!!!
*. do CI in space JCMBPSC
                MAXIT = ISEQCI(JSEQ,JCMBSPC)
*
                ISPC = JCMBSPC
                ISM = IREFSM
                ITERGLB = 1
*
*=============================================================
*
                CALL GASCI_REL(F1,F2,
     &                         INTSPACE,IREFSM,JCMBSPC,IPRDIA,
     &                         EREF,ECORE_KRMC,CIRUN,
     &                         WORK(KPROCLIST),WORK(KGROUPLIST))
*
*=============================================================
*
*. Last space where vectors were stored
                ISTOSPC = JCMBSPC
              END IF
            END DO
*
          END DO
*         ^ End loop over CI spaces
*
          if (NOMOFL.eq.0.or.IFINMO.eq.6) THEN
*. Create set of spinors from densities
c           WRITE(6,*) ' I am going to call FINMO '
c           CALL FINMO_REL
          end if
        end do
*       ^ End loop over natural spinor iterations (default=1)
*
*       release memory
*
        CALL MEMMAR(IDUM,IDUM,'FLUSM',IDUM,'GROUPL')
*
*       free groups
*
#if defined (VAR_MPI)
*
        CALL SETUNITS_PAR_CLOSE_REL
        CALL GROUP_DESTRUCTOR_REL(MYNEW_COMM,ICOMM)
*
#endif
*
*******************************************
      end if
*     ^ End if CIRUN.ne.'NDET  '
*******************************************
*
      CALL QEXIT('LUCIAR')
#ifdef LUCIAREL_STANDALONE
      Call Abend2( ' I am home from  the loops ' )
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_TASK(CIRUN)
C***********************************************************************
C
C     Introduce LUCIAREL action.
C
C     Written by  T. Fleig         - Nov 09 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      character*6 CIRUN
*
      NTEST = 0
*
      if (NTEST.ge.1) then
        write(LUPRI,*)
        write(LUPRI,*) ' ////////////////////////////////// '
        write(LUPRI,*) ' //                              // '
        write(LUPRI,*) ' //                              // '
        write(LUPRI,*) ' //  KR-MCSCF calling LUCIAREL   // '
        write(LUPRI,*) ' //                              // '
        write(LUPRI,*) ' //                              // '
        write(LUPRI,*) ' ////////////////////////////////// '
      end if
*
      if (CIRUN.eq.'NDET  ') then
        write(LUPRI,*)
        write(LUPRI,*) ' NDET   flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) ' The number of CI dets will be computed.'
        write(LUPRI,*)
      else if (CIRUN.eq.'CIINII') then
        write(LUPRI,*)
        write(LUPRI,*) ' CIINII flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Determining configurational start guess'
        write(LUPRI,*)
      else if (CIRUN.eq.'RSTRMC') then
        write(LUPRI,*)
        write(LUPRI,*) ' RSTRMC flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Determining info for restart of MCSCF.'
        write(LUPRI,*)
      else if (CIRUN.eq.'SIGMA ') then
        write(LUPRI,*)
        write(LUPRI,*) ' SIGMA  flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Computing a sigma vector.'
        write(LUPRI,*)
      else if (CIRUN.eq.'SIGMAD') then
        write(LUPRI,*)
        write(LUPRI,*) ' SIGMAD flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Computing a sigma vector.'
        write(LUPRI,*) '     Saving H diagonal also.'
        write(LUPRI,*)
      else if (CIRUN.eq.'DENS1 ') then
        write(LUPRI,*)
        write(LUPRI,*) ' DENS1  flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Computing only 1-particle density matrices.'
        write(LUPRI,*)
      else if (CIRUN.eq.'DENS2 ') then
        write(LUPRI,*)
        write(LUPRI,*) ' DENS2  flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Computing 1- and 2-particle dens. matrices.'
        write(LUPRI,*)
      else if (CIRUN.eq.'ANALYZ') then
        write(LUPRI,*)
        write(LUPRI,*) ' ANALYZ flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Analyzing the current CI vector(s).'
        write(LUPRI,*)
      else if (CIRUN.eq.'DIAG  ') then
        write(LUPRI,*)
        write(LUPRI,*) ' DIAG   flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Computing CI diagonal.'
        call quit('NOT IMPLEMENTED YET!')
        write(LUPRI,*)
      end if
*
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_NDET(NDET,ISYM,WORK,LWORK)
C***********************************************************************
C
C     Calculate number of determinants.
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision : 23.01.2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
C
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cicisp.inc"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*)
C
      CIRUN = 'NDET  '
      ISYM_KRMC = ISYM
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Compute the number of determinants for this setup and return.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      NDET = INT(XISPSM(ISYM_KRMC,1))
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_RCIST(CREF,CMO,WORK,LWORK)
C***********************************************************************
C
C     CI start guess.
C
C     Written by J. Thyssen and T. Fleig                  - Jan 12 2001
C     Last revision :                           T. Fleig  - Nov 30 2001
C                      Read integral files      T. Fleig  - May    2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      CHARACTER*6 CIRUN,RCISTFILE
      LOGICAL PARCAL_SAVE
C
      DIMENSION WORK(*),CMO(*),CREF(*)
C
#include "../luciarel/clunit.inc"
#include "../luciarel/lucipar_rel.inc"
C
#include "memint.h"
C
      CIRUN = 'CIINII'
      MZ = MIN(NZ,2)
C
      KFRSAV = KFREE
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
#if defined (VAR_MPI)
C
C     test fock matrix work for slaves
C
      ISET_PARCAL = 0
      IF( ISET_PARCAL .eq. 1 )THEN
        PARCAL_SAVE = .TRUE.
        CALL TEST_PAR_CALC(PARCAL_SAVE,0)
      END IF
#endif
C
C     Flag for type of integral interface:
C     I_IITP = 1: Old 4IND1XXXX0 driven interface
C     I_IITP = 2: New MRCONEE/MDCINT driven interface
C
      N_FCACM = (2*NASHT)*(2*NASHT)
C
      I_IITP = 1
C
      if (I_IITP.eq.2) then
C
        print*,' Testing new integral interface (MRCONEE/MDCINT)'
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDK,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDL,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLR,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLI,NZ,WORK,KFREE,LFREE)
        print*,'Allocated for ',N_FCACM,' integrals, MZ=',MZ
        CALL LCR_RD_MRCMDC(WORK(KFCACM),WORK(KFCACM+N_FCACM),
     &                     WORK(KINDK),WORK(KINDL),
     &                     WORK(KRKLR),WORK(KRKLI))
        CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,
     &               WORK,KFREE,LFREE)
        print*,'Dimension of MOLFDIR array : ',(2*NASHT)**4 * 2
        print*,'Dimension of (NZ,3)  array : ',NASHT*NASHT*NNASHX*NZ*3
        CALL RGETH21(.FALSE.,.TRUE.,DUMMY,DUMMY,WORK(KH2AC),
     &               WORK(KFREE),LFREE,IPROPT)
        stop
C
      else if (I_IITP.eq.1) then
C
C       Calculate core Hamiltonian:
C       ---------------------------
C
C       The core Hamiltonian is the active-active part of FC,
C       a.k.a. FCAC.
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
        IF( CRDFO_MAT ) THEN
C
C         read FCmo from file LUKRM3
C
          CALL KRMC_GETFCK(WORK(KFC),IPROPT,WORK(KFREE),LFREE)
C
          ECORE = ECORE_LR
C
        ELSE
C
C         Calculate FCmo ( CWRTFO_MAT false: do not write to LUKRM3) :
C
          CALL RGETFC(CWRTFO_MAT,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
        END IF
        IF (IPROPT.GE.1) THEN
          WRITE(LUPRI,'(A,1F18.10)') 'Core energy : ',ECORE
        END IF
C
C       Get active-active part of FCmo:
C
        CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A)') 'Active Fock matrix:'
          DO I=1,N2ASHXQ,1
            WRITE(6,'(1I3,1F18.14)') I,WORK(KFCAC+I-1)
          END DO
        END IF
C
C       Transform quaternion FCACmo to molfdir type FCAC.
C
        CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A)') 'Molfdir format FCAC:'
          WRITE(LUPRI,'(A)') 'Real part:'
          DO I=1,N_FCACM
            WRITE(6,'(1I3,1F18.14)') I,WORK(KFCACM+I-1)
          END DO
          IF (NZ.GT.1) THEN
            WRITE(LUPRI,'(A)') 'Imaginary part:'
            DO I=N_FCACM+1,MZ*N_FCACM
              WRITE(6,'(1I3,1F18.14)') I,WORK(KFCACM+I-1)
            END DO
          END IF
        END IF
C
        CALL MEMREL('after RGETFC',WORK,1,KFCAC,KFREE,LFREE)

C       read two-electron integrals in luciarel/transform_r.F 
C       (subroutine itrctl)
C
        LEN_ALL_INT = 0
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
        CALL MEMGET('REAL', KH2AC, 0, WORK, KFREE, LFREE)
C
C       FIXME: write me
C              SPINFR=.TRUE. not yet implemented in RGETH2 (why?)
C
      else
        CALL QUIT('LUCI_RCIST: Invalid flag I_IITP.')
      end if
#if defined (VAR_MPI)
C
C     end of slave reading test
C
      IF( ISET_PARCAL .eq. 1 ) THEN
        CALL TEST_PAR_CALC(PARCAL_SAVE,1)
      END IF
#endif
C
C     Transfer data to LUCIAREL common blocks.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call information transfer routine
C FIXME: Direct handling hardwired for the moment. Later, this
C        should be made dynamical and controlled from the outside.
C        Possibly also integrals reside on file in LUCIAREL
      IO_RCIST_IN = 0
CC    IF (IO_RCIST_IN.EQ.1) THEN
CC      CALL KR_LUCI_IO(0,0,1,DUMMY,1,DUMMY,1,
CC   &                  WORK(KFCACM),N2ASHXQ,WORK(KH2AC),
CC   &                  (2*NASHT)**4*2,
CC   &                  DUMMY,DUMMY,2,IRIQ,
CC   &                  IDUMMY,IDUMMY,IDUMMY,IDUMMY,IPROPT)
CC    END IF
C
C
C     Call LUCIAREL
      CALL KRMC_LUCI(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,ECORE,
     &               CIRUN)
C
C     Read CI vector from LUCIARELs file
      IF (MOPT_MXMACRO .ge. 0 .and. MAXCIT.GE.0.AND.NCIROOT.EQ.1) THEN
csk+hjaaj, 03JUL07: .MAX MACRO negative is special code for not writing
c   CREF to disk. Dirty fix for doing large scale CI in LUCIAREL, where
C   memory is too small for CREF.
        IO_RCIST_OUT = 1
* FIXME: Future use of file names instead of the primitive fortran
*        files. Activation flags (IO_RCIST_OUT) should be controlled
*        from the outside, not hardwired.
        IF (IO_RCIST_OUT.EQ.1) THEN
          CALL MEMGET('REAL',KRSCR,2*NZCONF,WORK,KFREE,LFREE)
          CALL KR_LUCI_IO(1,0,0,CREF,2*NZCONF,DUMMY,1,
     &                    DUMMY,1,DUMMY,1,
     &                    WORK(KRSCR),DUMMY,1,IRIQ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        END IF
      END IF
C
      CALL MEMREL('LUCI_RCIST',WORK,1,KFRSAV,KFREE,LFREE)
C
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_SIGMA(SAVDIA,DIAG,CB,SB,FCACM,H2ACM,
     &                      WORK,KFREE,LFREE)
C***********************************************************************
C
C     Calculate sigma vector: SB = H x CB, where H is the Hamiltonian
C     with integrals FCACM and H2ACM.
C     If (SAVDIA) also calculate H diagonal in DIAG
C
C     Written by J. Thyssen and T. Fleig         - Jan 12 2001
C     Last revision :       T. Fleig             - Mar 16 2002
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      LOGICAL   SAVDIA
      CHARACTER*6 CIRUN
      DIMENSION DIAG(*),CB(*),SB(*),FCACM(*),H2ACM(*)
      DIMENSION WORK(*)
C
#include "../luciarel/clunit.inc"
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      CIRUN = 'SIGMA '
      IF (SAVDIA) THEN
        CIRUN = 'SIGMAD'
      ELSE
        CIRUN = 'SIGMA '
      END IF
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Setup LUCIAREL common blocks.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Write C vector to LUCIAREL restart file:
C           FIXME: Control flag setting (IO_SIGMA_IN) from outside.
C                  (If modification to storage mode is attempted.)
C
      IO_SIGMA_IN = 1
CSK      IPROPT = 10
      IF (IO_SIGMA_IN.eq.1) THEN
        REWIND LUC
        CALL MEMGET('REAL',KRSCR,NZCONF,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,CB,NZCONF*MZ,
     &                        DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,2,IRIQ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
CSK      IPROPT = 0
C
C     Call LUCIAREL
      CALL KRMC_LUCI(FCACM,H2ACM,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
C     Read sigma vector on file to SB array.
      CALL MEMGET('REAL',KRSCR,2*NZCONF,WORK,KFREE,LFREE)
      IO_SIGMA_OUT = 1
      IF (IO_SIGMA_OUT.eq.1) THEN
        CALL KR_LUCI_IO(1,0,0,
     &                  SB,NZCONF*2,DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,1,IRIQ,
     &                  LUHC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
      CALL REWINE(LUC,-1)
      IO_C_VEC_OUT = 1
      IF (IO_C_VEC_OUT .eq. 1 ) THEN
        CALL KR_LUCI_IO(1,0,0,
     &                  CB,NZCONF*2,DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,1,IRIQ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
C
      IF (SAVDIA) THEN
C       Read diagonal on file LUDIA to WORK(KDIAG) array.
C
C       only "real part" of diagonal needed!
C       loop over ILOOP_DIA instead of IRIQ! /SK - June 30 2007
        ILOOP_DIA = 1
        CALL MEMGET('REAL',KDIAG,NZCONF,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(0,1,0,
     &                  DUMMY,1,WORK(KDIAG),NZCONF,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,1,ILOOP_DIA,
     &                  IDUMMY,LUDIA,IDUMMY,IDUMMY,IPROPT)
CSKold        CALL KR_LUCI_IO(0,1,0,
CSKold     &                  DUMMY,1,WORK(KDIAG),2*NZCONF,
CSKold     &                  DUMMY,1,DUMMY,1,
CSKold     &                  WORK(KRSCR),IDUMMY,1,IRIQ,
CSKold     &                  IDUMMY,LUDIA,IDUMMY,IDUMMY,IPROPT)
C       Copy "real part" of diagonal to DIAG.
        CALL DCOPY(NZCONF,WORK(KDIAG),1,DIAG,1)
      END IF
CSK      IPROPT = 0
C
      CALL MEMREL('LUCI_SIGMA',WORK,1,KFRSAV,KFREE,LFREE)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_DENS(CALCR1,CALCR2,RHO1,RHO2,CREFL,CREFR,
     &                     WORK,KFREE,LFREE)
C***********************************************************************
C
C     Read a bra vector CL and a ket vector CR and compute requested
C     density matrices RHO1 (and RHO2)
C
C     CIRUN = 'DENS1 ' : only one-particle density matrix
C     CIRUN = 'DENS2 ' : complete          density matrix
C
C     MOPT_MXMACRO < 0 : calculation of complete density matrix for 
C                        final root IDENSLR_STATE (dcbopt.h)
C
C     CREFL            : Left CI vector
C     CREFR            : Right CI vector
C       general also for transition densities
C
C     Written by T. Fleig                  - Jan 24 2001
C     Last revision : KR_LUCI_IO calls     - Dec 04 2001    T. Fleig
C                     CREF added           - Mar 17 2002    T. Fleig
C                     Only R1 possible     - Feb 01 2005    T. Fleig
C                     Only one common call 
C                     to KRMC_LUCI since
C                     RHO1,RHO2 are always
C                     allocated.           - Jun 29 2007    S. Knecht
C                     implementation of 
C                     multi-root
C                     natural orbital 
C                     occupation numbers   - Sep 28 2007    S. Knecht
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
      LOGICAL CALCR1,CALCR2
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*),RHO1(*),RHO2(*),CREFL(*),CREFR(*)
C
#include "../luciarel/clunit.inc"
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      if (CALCR2) then
        CIRUN = 'DENS2 '
        if (.not.CALCR1)
     &       call quit('Not implemented option in LUCI_DENS.')
      else if (CALCR1) then
        CIRUN = 'DENS1 '
      else
        call quit('Check CALCR setting in LUCI_DENS')
      end if
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Initialize density matrices to zero
      L1 = (2*NASHT)**2 * MZ
      CALL DZERO(RHO1,L1)
      if (CALCR2) then
        L2 = (2*NASHT)**4 * MZ
        CALL DZERO(RHO2,L2)
      end if
C
C
      IF( NATOLCR ) THEN
C       
C       CI vector resides on LUCVECT_REL(.0)
C       --> vector copy in GASCI_REL
C
      ELSE
C
C       Write current CI vector(s) to file
        IO_DENS_IN = 1
        IF (IO_DENS_IN.EQ.1) THEN
          CALL MEMGET('REAL',KRSCR,NZCONF,WORK,KFREE,LFREE)
          CALL KR_LUCI_IO(1,1,0,CREFL,NZCONF*MZ,
     &                          CREFR,NZCONF*MZ,
     &                    DUMMY,1,DUMMY,1,
     &                    WORK(KRSCR),IDUMMY,2,IRIQ,
     &                    LUSC61,LUSC62,IDUMMY,IDUMMY,IPROPT)
          CALL MEMREL('LUCI_DENS 1',WORK,1,KRSCR,KFREE,LFREE)
        END IF
      END IF
C
C     Density run of LUCIAREL only
C
C     RHO1 and RHO2 are always allocated in the calling routines.
C     RHO2 is essential for communication in parallel runs.
C     / SK June 29 2007.
C
CSK      if (CALCR2) then
        CALL KRMC_LUCI(RHO1,RHO2,WORK(KFREE),LFREE,DUMMY,CIRUN)
CSK      else
CSK        CALL KRMC_LUCI(RHO1,DUMMY,WORK(KFREE),LFREE,DUMMY,CIRUN)
CSK      end if
C
C     Write CI densities to file or keep the arrays
C     Currently, densities are directly passed up to LUCI_DENS.
c     IO_DENS_OUT = 0
c     IF (IO_DENS_OUT.EQ.1) then
c       CALL MEMGET('REAL',KRSCR,2*NZCONF,WORK,KFREE,LFREE)
c       CALL KR_LUCI_IO(0,0,1,
c    &                  DUMMY,1,DUMMY,1,
c    &                  WORK(KRSCR),IDUMMY,1,IRIQ,
c    &                  LUDUM,IDUMMY,IDUMMY,IDUMMY,IPROPT)
c     END IF
C
C     Do symmetrizations.
C     FIXME: write me  ??
C
      CALL MEMREL('LUCI_DENS',WORK,1,KFRSAV,KFREE,LFREE)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_ANAL(CREF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Analyze CI vectors
C
C     Written by T. Fleig                  - Jan 24, 2001
C     Last revision : T. Fleig
C                     activated            - Jan 05, 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
#include "../luciarel/clunit.inc"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*)
C
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
      CIRUN = 'ANALYZ'
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C      WRITE(6,*) 'NCIROOT and MOPT_MXMACRO',NCIROOT,MOPT_MXMACRO
C
      IF ( NCIROOT.EQ.1 .AND. MOPT_MXMACRO .gt. 0) THEN
        IO_ANAL_IN = 1
        IF (IO_ANAL_IN.EQ.1) THEN
          REWIND LUC
          CALL MEMGET('REAL',KASCR,NZCONF*MZ,WORK,KFREE,LFREE)
          CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,
     &                          DUMMY,1,
     &                    DUMMY,1,DUMMY,1,
     &                    WORK(KASCR),IDUMMY,2,IRIQ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        END IF
      ELSE IF (NCIROOT.GT.1) THEN
        CONTINUE
      END IF
C
C     Call LUCIAREL for analysis of CI vectors.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
      IF (NCIROOT.EQ.1) THEN
        CALL MEMREL('LUCI_ANALYZ',WORK,1,KFRSAV,KFREE,LFREE)
      END IF
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
#ifdef NOT_USED_ROUTINE
  NB! NZCONF not defined here, but is used below
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_DIAG(DIAG,WORK,LWORK)
C***********************************************************************
C
C     Compute CI diagonal
C
C     Written by T. Fleig                  - Dec 04 2001
C     Last revision :   KR_LUCI_IO calls   - Mar 17 2002    T. Fleig
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dgroup.h"
#include "krmcluci_inf.h"
C
#include "../luciarel/clunit.inc"
C
#include "memint.h"
C
      CHARACTER*6 CIRUN,CIDIAGFILE
C
      DIMENSION WORK(*)
C
      CIRUN = 'DIAG  '
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
C     FIXME: Check this here.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
      IO_DIAG_IN = 0
C FIXME: Eventually integrals for diagonal to file.
      IF (IO_DIAG_IN.EQ.1) THEN
c       CALL KR_LUCI_IO(0,0,1,DUMMY,1,DUMMY,1,
c    &                  WORK(KFCACM),N2ASHXQ,WORK(KMUUUU),
c    &                  (2*NASHT)**4*2,
c    &                  DUMMY,DUMMY,2,IRIQ,
c    &                  LUDUM,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
C
C     Call LUCIAREL for computation of CI diagonal.
c     CALL KRMC_LUCI(WORK(Oelints),WORK(Telints),WORK,LWORK,DUMMY,CIRUN)
C
C     Read diagonal vector from LUCIARELs file
      IO_DIAG_OUT = 1
      IF (IO_DIAG_OUT.EQ.1) THEN
        CALL MEMGET('REAL',KRSCR,2*NZCONF,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,DIAG,2*NZCONF,DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),DUMMY,1,IRIQ,
     &                  LUDIA,IDUMMY,IDUMMY,IDUMMY,IPROPT)
      END IF
C
      CALL MEMREL('LUCI_DIAG',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
#endif /* NOT_USED_ROUTINE */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_RSTRMC(WORK,LWORK)      
C***********************************************************************
C
C     Restart information for LUCIAREL in MCSCF runs.
C
C     Written by T. Fleig                 - Aug 28, 2003
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
C
#include "dcbopt.h"
C
      CHARACTER*6 CIRUN
      DIMENSION WORK(*)
C
#include "../luciarel/cstate.inc"
C
      CIRUN = 'RSTRMC'
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
      CALL LUCI_TRACOM(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call LUCIAREL for setting up vector partitioning information.
      CALL KRMC_LUCI(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE KR_LUCI_IO(IVA,IVB,II,VA,IDIMA,VB,IDIMB,
     &                      X1,IDIM1,X2,IDIM2,
     &                      SCR,ISCR,IOTYPE,IRILP,
     &                      LU1,LU2,LU3,LU4,IPRNT)
C***********************************************************************
C
C     Transfer information between LUCIAREL and KR-MCSCF
C
C     Written by T. Fleig                  - Nov 28 2001
C     Last revision : T. Fleig             - Mar    2002
C
C
C  Different cases:
C
C    1) Vector 1           IVA
C    2) Vector 2           IVB
C    3) Integral quantity  II
C     Quantity specifier:  0 no quantity
C                          1 quantity
C
C   I/O ?
C
C   Read file to array     IOTYPE = 1
C          for output   (KRMC <- LUCIAREL).
C
C   Write array to file    IOTYPE = 2
C          for input    (KRMC -> LUCIAREL)
C
C***********************************************************************
#include "implicit.h"
*
      dimension VA(IDIMA),VB(IDIMB),X1(IDIM1),X2(IDIM2),SCR(*)
      dimension ISCR(*)
*
      NTESTL = 00
CSK      IPRNT = 000
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.5) then
*   Print section:
*===========================================
        write(6,*)
        write(6,*) '#########################'
        write(6,*) ' KR_LUCI_IO activity     '
        write(6,*) '#########################'
        if (IVA.eq.1.or.IVB.eq.1) then
          write(6,*) '  Vector task   '
        else if (II.eq.1) then
          write(6,*) '  Integral / density task   '
        end if
*
* A) Output from LUCIAREL (read from file)
        if (IOTYPE.eq.1) then
          if (NTEST.ge.1) write(6,*) '  Output from LUCIAREL '
*
* B) Input for LUCIAREL (write to file)
        else if (IOTYPE.eq.2) then
          if (NTEST.ge.1) write(6,*) '  Input for LUCIAREL '
        end if
*===========================================
      end if
*
*===============
*  Load vectors
*===============
      if (IOTYPE.eq.1) then
*  Load vector A (ket)
        if (IVA.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Loading vector A ***'
          call rdvcdc(VA,SCR,LU1,-1,IRILP,IPRNT)
        end if
*
*  Load vector B (bra)
        if (IVB.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Loading vector B ***'
          call rdvcdc(VB,SCR,LU2,-1,IRILP,IPRNT)
        end if
*
*  Load 1-particle quantity
*  Load 2-particle quantity
        if (II.eq.1) then
          call rdvcdc(X1,SCR,LU3,-1,IRILP,IPRNT)
          call rdvcdc(X2,SCR,LU4,-1,IRILP,IPRNT)
        end if
*
*===============
*  Save vectors
*===============
*  (Write as REAL PART, -1 separator, IMAG PART, -1 trailer)
      else if (IOTYPE.eq.2) then
*  Save vector A (ket)
        if (IVA.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Saving vector A ***'
          call savevcdc(VA,SCR,LU1,IRILP,IPRNT)
        end if
*
*  Save vector B (bra)
        if (IVB.eq.1) then
          if (NTEST.ge.10) write(6,*) '   *** Saving vector B ***'
          call savevcdc(VB,SCR,LU2,IRILP,IPRNT)
        end if
*
*  Save integral quantity
        if (II.eq.1) then
          write(6,*) 'Saving integrals to file not implemented yet.'
          stop 'Quitting.'
        end if
*
      end if
*     ^ type of IO
*
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUMSTI_REL
C
C     Call the slaves if we run in parallel for KRMC-LUCIAREL
C
C     adaption of the corresp. RELCCSD routine written by Luuk Visscher
C     KRMC-LUCIAREL adaption: S. Knecht
C
C     Last revision: S. Knecht           - March  2007
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbgen.h"
#include "../luciarel/lucipar_rel.inc"
C
C     ITASK = 6 for parallel KRMC-LUCIAREL
C
      IF (PARCAL) CALL DIRAC_PARCTL( 6 )
CSK      write(6,*)'slaves should be up und running'
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUMSTE_REL
C
C     Release the slaves if we run in parallel for KRMC-LUCIAREL
C
C     adaption of the corresp. RELCCSD routine written by Luuk Visscher
C     KRMC-LUCIAREL adaption: S. Knecht
C
C     Last revision: S. Knecht           - March  2007
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "dcbgen.h"
#include "infpar.h"
#include "../luciarel/lucipar_rel.inc"
C
      IF (.NOT.PARCAL) RETURN
C
C     Find the slaves and release them

      NTEST = -1
#if defined (VAR_MPI)
      DO I = 1,NUMNOD
         CALL MPI_IRECV(NODE,1,my_MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                  MPI_COMM_WORLD,IREQ,IERR)
         CALL MPI_WAIT(IREQ,ISTAT,IERR)
CSK         write(6,*)'node recieved',node,mparid
         CALL MPI_SEND(NTEST,1,my_MPI_INTEGER,NODE,30,
     &                  MPI_COMM_WORLD,IERR)
      ENDDO
#endif

C
C     Let the slaves sleep (if relevant for this operating system)
C
      CALL DIRAC_PARCTL( 0 )
C
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCINOD_REL(WORK,LWORK)
C
C     KRMC-LUCIAREL driver routine for the slaves
C
C     adapted version of the slave RELCCSD routine 
C     written by Luuk Visscher.
C     KRMC-LUCIAREL adaption: S. Knecht
C
C     Last revision: S. Knecht           - March  2007
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
#include "../luciarel/lucipar_rel.inc"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION WORK(LWORK)
      LOGICAL LUEND
      CHARACTER LUCITABASF_REL*12, LUCIFILN_REL*16
      CHARACTER CIRUN*6
C
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cicisp.inc"
#include "../luciarel/noccn_inf.inc"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "memint.h"
C
      KFRSAV = KFREE
*
*       START 
*     =========
*
*     arrange for the MPI stuff and correct node number
*     to the total number of running invocations.
*
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node, NUMNOD = number of slaves
      NMPROC = NUMNOD + 1
C
C     initialize some variables ...
C
      ECORE_KRMC = 0.0D0
*
*     create a node-unique filename as output file. Important on
*     shared file systems. Otherwise all the output gets mingled in one
*     file. You don't really want to do this.
*
      LUCITABASF_REL="LUCIAREL.OUT"
*
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (LUCIFILN_REL,'(A12,A1,I1)') LUCITABASF_REL,'.',MYPROC
         LUFIL=14
      ELSE IF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (LUCIFILN_REL,'(A12,A1,I2)') LUCITABASF_REL,'.',MYPROC
         LUFIL=15
      ELSE IF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (LUCIFILN_REL,'(A12,A1,I3)') LUCITABASF_REL,'.',MYPROC
         LUFIL=16
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCINOD_REL MODULE")
      ENDIF
*
*     open the local input file and the node specific output file.
*     Every access to the local stdout handle then automatically writes
*     to the corresponding output file.
*
*
      OPEN(MSLVOUT_REL,FILE = LUCIFILN_REL(1:LUFIL))
C
C     transfer file handle to common block
C
      LUWRT = MSLVOUT_REL
C
C     synchronize nodes, part 1
C
      CALL SYNC_NODES_REL_P1(ECORE_KRMC,CIRUN,MSLVOUT_REL)
C
C
      MZ = MIN(NZ,2)
C
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or. 
     &    CIRUN .eq. 'SIGMAD') THEN
C
C       read-in of integrals by nodes... no broadcast in SYNC_NO...P2
C
        IF( CIRUN .eq. 'CIINII' .and. IIOMOD_REL .eq. 1 )THEN
C
          PARCAL = .FALSE.
C
        END IF
C
C
C     Flag for type of integral interface:
C     I_IITP = 1: Old 4IND1XXXX0 driven interface
C     I_IITP = 2: New MRCONEE/MDCINT driven interface
C
      N_FCACM = (2*NASHT)*(2*NASHT)
C
      I_IITP = 1
C
      if (I_IITP.eq.2) then
C
        print*,' Testing new integral interface (MRCONEE/MDCINT)'
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDK,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDL,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLR,NZ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KRKLI,NZ,WORK,KFREE,LFREE)
        print*,'Allocated for ',N_FCACM,' integrals, MZ=',MZ
        CALL LCR_RD_MRCMDC(WORK(KFCACM),WORK(KFCACM+N_FCACM),
     &                     WORK(KINDK),WORK(KINDL),
     &                     WORK(KRKLR),WORK(KRKLI))
        CALL MEMGET('REAL',KH2AC,NASHT*NASHT*NNASHX*NZ*3,
     &               WORK,KFREE,LFREE)
        print*,'Dimension of MOLFDIR array : ',(2*NASHT)**4 * 2
        print*,'Dimension of (NZ,3)  array : ',NASHT*NASHT*NNASHX*NZ*3
        CALL RGETH21(.FALSE.,.TRUE.,DUMMY,DUMMY,WORK(KH2AC),
     &               WORK(KFREE),LFREE,IPROPT)
        stop
C
      else if (I_IITP.eq.1) then
C
C       Calculate core Hamiltonian:
C       ---------------------------
C
C       The core Hamiltonian is the active-active part of FC,
C       a.k.a. FCAC.
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KFCACM),MZ*N_FCACM)
C
C       Read two-electron integrals:
C       ----------------------------
C
C       ... now inside LUCIAREL (subroutine itrctl)
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
C
        IF( CIRUN .eq. 'CIINII' )THEN
          CALL MEMGET('REAL',KH2AC,          0,WORK,KFREE,LFREE)
        ELSE
          CALL MEMGET('REAL',KH2AC,LEN_ALL_INT,WORK,KFREE,LFREE)
        END IF
C
      else
        CALL QUIT('LUCI_RCIST: Invalid flag I_IITP.')
      end if
C
C
        PARCAL = .TRUE.
C
C       synchronize nodes, part 2 (now CIRUN dependent!)
C
        CALL SYNC_NODES_REL_P2(WORK(KFCACM),WORK(KH2AC),CIRUN)
C
      ELSE IF( CIRUN .eq. 'DENS1' .or. CIRUN .eq. 'DENS2' ) THEN 
C
C       Initialize density matrices to zero; currently I want to test
C       this CIRUN, therefore we make no difference between DENS1 and
C       DENS2 - SK June 2007
C
        L1 = (2*NASHT)**2 * MZ
C
        CALL MEMGET('REAL',KRHO1,L1,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KRHO1),L1)
C
CSK        IF (CIRUN .eq. 'DENS2' ) THEN
C
        L2 = (2*NASHT)**4 * MZ
        LF2_ZERO = (2*NASHT)**4 * MZ
C
CSK        END IF
C
      END IF
C
C     Transfer data to LUCIAREL common blocks.
C
      CALL LUCI_TRACOM(CIRUN,1)
C
C     attention to nat. orb. occ. num. calculation
C
      IF( CIRUN .eq. 'DENS1' .or. CIRUN .eq. 'DENS2' ) THEN
       IF( .NOT. NOOCCN_DENS_RUN ) THEN
          CALL MEMGET('REAL',KRHO2,L2,WORK,KFREE,LFREE)
          CALL DZERO(WORK(KRHO2),L2)
       ELSE
          CALL MEMGET('REAL',KRHO2, 0,WORK,KFREE,LFREE)
       END IF
      END IF
C
C     Setup routine for running LUCIAREL
C
      CALL LUCI_SETUP()
C
C     Run LUCIAREL
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or. 
     &    CIRUN .eq. 'SIGMAD') THEN
C
C       Call LUCIAREL for calculating sigma vector(s).
        CALL KRMC_LUCI2(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,
     &                 ECORE_KRMC,CIRUN)
C 
        CALL MEMREL('LUCINOD_RL',WORK,1,KFRSAV,KFREE,LFREE)
C
      ELSE IF( CIRUN .eq. 'DENS1' .or. CIRUN .eq. 'DENS2' ) THEN 
C
C       Call LUCIAREL for calculating density matrix elements.
        CALL KRMC_LUCI2(WORK(KRHO1),WORK(KRHO2),WORK(KFREE),
     &                    LFREE,DUMMY,CIRUN)
C
        CALL MEMREL('LUCINOD_RL',WORK,1,KFRSAV,KFREE,LFREE)
C
      ELSE IF( CIRUN .eq. 'ANALYZ' ) THEN 
C
C       Call LUCIAREL for analyzing the wavefunction.
        CALL KRMC_LUCI2(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
      ELSE IF( CIRUN .eq. 'NDET  ' ) THEN 
C
C       Call LUCIAREL for determining number of dets.
        CALL KRMC_LUCI2(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
        NDET = INT(XISPSM(ISYM_KRMC,1))
C
      ELSE IF( CIRUN .eq. 'RSTRMC' ) THEN
C
C       Call LUCIAREL for setting up vector partitioning information.
        CALL KRMC_LUCI2(DUMMY,DUMMY,WORK,LWORK,DUMMY,CIRUN)
C
      END IF
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
C     LUCIAREL run is finished for the slaves ...
C
      CLOSE (MSLVOUT_REL,STATUS='KEEP')
C
C     .. wait for the master
C
      LUEND = .FALSE.
C
      CALL MPI_ISEND(MYTID,1,my_MPI_INTEGER,MPARID,20,
     &              MPI_COMM_WORLD,IREQ,IERR)

    1    CONTINUE
         CALL SLEEP(1)
         CALL MPI_IPROBE(MPARID,30,MPI_COMM_WORLD,LUEND,ISTAT,IERR)
         IF (.NOT.LUEND) GOTO 1

      CALL MPI_RECV(NTEST,1,my_MPI_INTEGER,MPARID,30,
     &               MPI_COMM_WORLD,ISTAT,IERR)
C
      END
#endif    /* ifdef VAR_MPI */

