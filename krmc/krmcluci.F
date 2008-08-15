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
#include "parluci.h"
C
C     set "I/O node master" + "shared memory node master"
      N_MASTER    = 0
      N_MASTER_SM = 0
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
#include "parluci.h"
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
#if defined (VAR_MPI2)
C
      CALL SYNC_NODES_REL_P1(ECORE_KRMC,CIRUN,LUPRI)
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or.
     &    CIRUN .eq. 'SIGMAD' .or. CIRUN .eq. 'IJKLRO' .or. 
     &    CIRUN .eq. 'KR-CI ' ) THEN
C
        CALL SYNC_NODES_REL_P2(F1,F2,CIRUN)
C
      END IF
C
C     Transfer data to LUCIAREL common blocks (only for NODES) 
C     --> communication is enabled.
C
      CALL TRKRLUCI(CIRUN,1)
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
#include "dgroup.h"
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
#include "../luciarel/ctcc.inc"
CTF
      character*6 CIRUN
      CHARACTER*72 MACHINENAME
      INTEGER NAMELENGTH, MZ_NUM
      logical MCSCF, RUN_SHMEM_T
CTF  global reference values for Z_BLKFO needed in MV7 !!
      common/glbref/ISPC,ISM,ITERGLB
      DIMENSION F1(*),F2(*),WRK(LWRK)
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
C
C  FIXME:
C         VCT1 and VCT2 unused at the moment. They should carry
C         input/output vectors controlled by outside KR_LUCI_IO
C         routine. At the moment, allocation is performed inside
C         LUCIAREL. The vectors reside on the internal files
C         fort.20 (D), fort.21 (C), and fort.31 (S).
C         Renamed files. See diskun.f
C
********************************************
      CALL QENTER('LUCIAR')
C
      MZ_NUM = 0
      MZ_NUM = MIN(NZ,2)
      MIN_ONE = - 1
C
C     Compute offset from WRK to WORK
C
      CALL COMPOFF(WRK,WORK,K_OFFSET)
C
C     *k_offset = (long) (wrk-work+1);
C
      MXPWRD = LWRK
      KADD = MXPWRD
C     initialize memory
      call memmar(K_OFFSET,KADD,'INI   ',IDUMMY,'DUMMY')
C
      NTEST = 0
      IMEMCK = 0
      MCSCF = .true.
C
C     debugging print for parallel runs -- control from outside?
C
      NPTEST_VAR = 0
C
C     Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_REL(WORK(1),IPRSTR)
C
C     Internal subspaces
      CALL LCISPC_REL(LCSBLK,IPRCIX)
C
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
C
*******************************************
      if (CIRUN.ne.'NDET  ') then
*******************************************
C
C       setting up groups of processors (stored in /LUCIPARREL/)
C
        CALL MEMMAR(IDUM,IDUM,'MARK ',IDUM,'GROUPL')
C
        CALL MEMMAR(KPROCLIST,  NMPROC,'ADDS  ',1,'IHOSTL')
        CALL MEMMAR(KGROUPLIST, NMPROC,'ADDS  ',1,'IHOST2')
        CALL MEMMAR(KPROCLENGTH,NMPROC,'ADDS  ',1,'INAMEL')
        CALL MEMMAR(KPROCLIST2, NMPROC,'ADDS  ',1,'IHOST3')
C
        CALL ISETVC(WORK(KPROCLIST),  MIN_ONE,NMPROC)
        CALL ISETVC(WORK(KGROUPLIST), MIN_ONE,NMPROC)
        CALL ISETVC(WORK(KPROCLENGTH),MIN_ONE,NMPROC)
        CALL ISETVC(WORK(KPROCLIST2), MIN_ONE,NMPROC)
C       initialize number of process groups sharing a c-vector file
        NFLGRPS_REL = 0
C
#if defined (VAR_MPI2)
C
C       hardwired to 1, otherwise it fails...
C
        IPARALLELIO = 1
C
#if defined (SYS_LOCALDISKS)
C
C       checking only if parallel file I/O in general is requested
C
        IF( IPARALLELIO .eq. 1 ) THEN
C  
           NAMELENGTH = 0
C
           CALL MPI_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH,IERR)
C
           CALL FIND_GROUP_OF_PROCS_REL(MACHINENAME,NAMELENGTH,
     &                                  WORK(KPROCLIST),
     &                                  WORK(KPROCLENGTH))
C
C          number of groups sharing a file ( = NFLGRPS ) is determined
C          in the subroutine above. If there is only one group
C          than we can simulate a shared filesystem.
C
        END IF
#else
C       in case of a parallel filesystem, we should use it
        NFLGRPS_REL = 1
        CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
#endif /* ifdef SYS_LOCALDISKS */
C       no integral broadcast necessary
        IF( NFLGRPS_REL .eq. 1 ) IIOMOD_REL = 1
C
C       MPI "shared memory" mode
C
        IF( SHARED_M )THEN
C
          IF( LEVEL_SM .lt. 6 ) THEN
C
C            level 0 and 5: build node-groups for shared memory usage
C
             NAMELENGTH = 0
C
             CALL MPI_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH,IERR)
C
             CALL ISETVC(WORK(KPROCLENGTH),-1,NMPROC)
C
             CALL FIND_GROUP_OF_PROCS_REL_SM(MACHINENAME,NAMELENGTH,
     &                          WORK(KPROCLIST2),WORK(KPROCLENGTH))
C
          ELSE
C
C            level 6:   modelling global shared memory 
C                       NOTE: you need a fast inter-node connect,
C                           e.g. InfiniBand, Myrinet, etc ...
C
             CALL ISETVC(WORK(KPROCLIST2),1,NMPROC)
C           
          END IF
C
        ELSE
C
           CALL ISETVC(WORK(KPROCLIST2),1,NMPROC)
C
        END IF
C
C       set-up all useful MPI stuff: communication groups, MPI-files ...
C                                    all important variables are stored on 
C                                    common block /LUPARGROUP/
C
C       construct new communication groups
        CALL GROUP_CONSTRUCTOR_REL(WORK(KGROUPLIST),WORK(KPROCLIST),
     &                             WORK(KPROCLIST2))
C
        NPTESTVAR = 00
        IF( NPTESTVAR .ge. 20 ) THEN
          WRITE(LUWRT,*) 'this is my old and new rank ',myproc, MYNEW_ID
          WRITE(LUWRT,*) 'this is my "shared memory" rank ', MYNEW_ID_SM
          WRITE(LUWRT,*) 'this is my file grouplist:'
          CALL IWRTMAMN(WORK(KGROUPLIST),1,NMPROC,1,NMPROC,LUWRT)
          WRITE(LUWRT,*) 'this is my complete proclist:'
          CALL IWRTMAMN(WORK(KPROCLIST),1,NMPROC,1,NMPROC,LUWRT)
          WRITE(LUWRT,*) 'this is my shared memory proclist:'
          CALL IWRTMAMN(WORK(KPROCLIST2),1,NMPROC,1,NMPROC,LUWRT)
        END IF
        NPTESTVAR = 0
C
C       open the files for parallel file I/O
C
        CALL SETUNITS_PAR_OPEN_REL(MY_GROUPN)
C
#else
        CALL ISETVC(WORK(KPROCLIST),1,NMPROC)
        NFLGRPS_REL = 0
#endif
C
C       Allocate work space, general
C 
        CALL ALLOC_REL
C
C       IDGSTYP determines the integral classes and corresponding operator
C       classes depending on the double group in use.
C
        IDGSTYP = -1
        if (PNTGRP.le.4) then
          IDGSTYP = 1
        else if (PNTGRP.eq.6.or.PNTGRP.eq.5) then
          IDGSTYP = 3
        else if (PNTGRP.ge.7.and.PNTGRP.le.9) then
          IDGSTYP = 2
        end if
        if (SIGDEN_ROUTE.eq.'SIGDEN') then
C         Space for Hamiltonian in form suitable for SIGDEN
          if (ISO1.eq.0.and.ISO2.eq.0) then
C           Non-relativistic (or spinfree) Hamiltonian
            IHTYPE = 1
          else if (ISO1.eq.1.and.ISO2.eq.0) then
C           Mean-field SO Hamiltonian
            IHTYPE = 2
          else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).ne.'DIRAC ')
C           Relativistic 1- and 2-particle Hamiltonian
     &         then
            IHTYPE = 3
          else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).eq.'DIRAC ')
C           Dirac-Coulomb case.
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
          IF (CIRUN .EQ. 'DENS1 ' .or. CIRUN .EQ. 'PROP1 ') THEN
             NO_TTWO = 1
C            ... no two-electron excitations
          ELSE
             NO_TTWO = 0
          END IF
          call set_hop_dbg(IHTYPE,NO_TTWO,IPRHAM_CI)
        end if
C       ^ end if SIGDEN_ROUTE
C
C       symmetry of reference
C
        IF(PNTGRP.GT.1.and.PNTGRP.lt.5)
     &             CALL MLSM_REL(IREFSM,IREFPA,IREFSM,'CI',1)
        IF(NOINT.EQ.1) THEN
          WRITE(LUWRT,*) ' End of calculation without integrals'
          CALL QEXIT('LUCIAR')
          RETURN
        END IF
C
        do INATIT=1,NATITER,1
C
C         Generate symmetry info and the like
C
          if (INATIT.eq.1) call luc2pic(ISO1,IPRTRA)
C
C         Generate and store pointers for transformed integrals
C         (and integrals !)
C         SK - 22-11-2007: not in use for SIGDEN_ROUTE == SIGDEN
C
          IF (NOINT.EQ.0.and.INATIT.eq.1) THEN
             IF (SIGDEN_ROUTE .ne. 'SIGDEN') CALL INTIM_REL(IPRTRA)
          ELSE
             WRITE(LUWRT,*)
     &       ' No integrals imported (or > 1st NS iteration) '
          END IF
C
          RUN_SHMEM_T = .FALSE.
C
          IF( IT_SHL .ge. 0 .or. SPLIT_IJKL ) THEN
            RUN_SHMEM_T = .TRUE.
            IF( CIRUN .eq. 'IJKLRO' ) RUN_SHMEM_T = .FALSE.
          END IF
C
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          if (CIRUN.ne.'ANALYZ') then
CSK         IDUMMY ==> IT_TTPL
            call picasso(INTSPACE,ISO1,IPRTRA,MCSCF,
     &                   NL2D,IBU2D,ISM2D,F1,F2,DUMMY1,DUMMY2,
     &                   IDUMMY,IBOSYM,CIRUN,RUN_SHMEM_T)
          end if
*+++++++++++++++++++++++++++++++++++++++++++++++++++
          IF ( CIRUN .eq. 'IJKLRO' ) GOTO 300
C
C
C ======================
C CI in GAS spaces
C ======================
C
          ISTOSPC = 0
          IF(IRESTR.EQ.1) ISTOSPC = 1
C
C. Restarted Lambda calculations need special attention
C  since restart is realized in in second calculation
C  First calculation is used to establish H0
C
          IRESTR_ORIG=IRESTR
          IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0) THEN
C
            write(LUWRT,*) 'THIS ROUTE IS STILL FORBIDDEN !!'
            write(LUWRT,'(A16,I2,1F3.1)')
     &            'IRESTR,XLAMBDA =',IRESTR,XLAMBDA
            Call Abend2( 'Quitting.' )
C
          END IF
C         ^ End of special handling of restarted calc with
C           lambda modified op.
C
          ISKIPEI_INI = ISKIPEI
          DO JCMBSPC = 1, NCMBSPC
            ISKIPEI = 0
C
C. Perform GASCI/GASPT calculation for internal space ICISPC
C
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
C
            IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0.AND.JCMBSPC.EQ.1) THEN
              WRITE(LUWRT,*)
     &              ' Remember No restart in calc 1 (Lambda calc)'
              IRESTR = 0
              Call Abend2( 'Quitting.' )
            END IF
C
            I_EXPAND = 1
            IF(XLAMBDA.NE.1.0D0 .AND.JCMBSPC.GT.1) THEN
              WRITE(LUWRT,*) ' =================================='
              WRITE(LUWRT,*) '   Modified operator will be used'
              WRITE(LUWRT,*) ' =================================='
              IF(JCMBSPC.EQ.2) THEN
C. Modify the operator
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
C. No expansion should follow
              I_EXPAND = 0
C. But normal restart
              IRESTR = 1
            END IF
C
            DO JSEQ = 1,NSEQCI(JCMBSPC)
              IF((ISEQCI(JSEQ,JCMBSPC).GT.0
     &            .OR.ISEQCI(JSEQ,JCMBSPC).EQ.-5)
     &            .AND. JCMBSPC.NE.1.AND.ISTOSPC.NE.JCMBSPC
     &            .AND.I_EXPAND.EQ.1)                   THEN
C. Restart from previous spaces ( Assuming a progressing sequence :
C  spaces are just added, not subtracted )
C  ( Used for perturbation calculations and CI calculations )
C
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
C. Last space where vectors were stored
                ISTOSPC = JCMBSPC
                ISKIPEI = ISKIPEI_INI
C. Expanded vector will be used as initial vector in the
C. zero space calculation. Tell next CI to restart from
C. CI vectors
                IRESTR = 1
              END IF
              IF(ISEQCI(JSEQ,JCMBSPC).GT.0) THEN
                IF(JSEQ.EQ.1.AND.JCMBSPC.EQ.2.AND.IRST2.EQ.0) THEN
C  No restart from previous vectors - IRST2 has been set to zero
                  IRESTR = 0
                  WRITE(LUWRT,*) ' No restart from previous vectors'
                END IF
C
C               Good old normal CI !!!!
C               do CI in space JCMBPSC
C
                MAXIT = ISEQCI(JSEQ,JCMBSPC)
C
                ISPC = JCMBSPC
                ISM = IREFSM
                ITERGLB = 1
C
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
C
 300      CONTINUE
C
        end do
*       ^ End loop over natural spinor iterations (default=1)
*
*       release memory
*
        CALL MEMMAR(IDUM,IDUM,'FLUSM',IDUM,'GROUPL')
C
#if defined (VAR_MPI2)
C
C       free groups
C
        CALL SETUNITS_PAR_CLOSE_REL
        CALL GROUP_DESTRUCTOR_REL(MYNEW_COMM,MYNEW_COMM_SM,ICOMM,
     &                            MYNEW_COMM_SM_C)
C
#endif
C
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
      else if (CIRUN.eq.'KR-CI ') then
        write(LUPRI,*)
        write(LUPRI,*) ' KR-CI  flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Running large-scale CI calculation'
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
      else if (CIRUN.eq.'IJKLRO') then
        write(LUPRI,*)
        write(LUPRI,*) ' IJKLRO flag call '
        write(LUPRI,*) ' ------           '
        write(LUPRI,*) '   Reordering all (ij|kl) according to GAS     '
        write(LUPRI,*) '   scheme. Saving to file IJKL_REOD.           '
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
      CALL TRKRLUCI(CIRUN,0)
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
C     Driver routine for LUCIAREL CI module.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbopt.h"
      DIMENSION CREF(*), CMO(*), WORK(*)
C
C     which module is calling?
C
      IF( JKRRUNTYPE .eq. 0 ) THEN
         WRITE(LUPRI,'(A,I3)') ' *** ERROR in LUCI_RCIST ***'//
     &   ' LUCI_RCIST not in use for 2nd order SCF!',JKRRUNTYPE
         CALL QUIT('*** ERROR in LUCI_RCIST ***')
      ELSE IF ( JKRRUNTYPE .eq. 1 ) THEN
C
C        get CI start guess for KR-MCSCF wave function optimization
C
         CALL LUCI_CIINII(CREF,CMO,WORK,LWORK)
      ELSE IF ( JKRRUNTYPE .eq. 2) THEN
C
C        run large-scale KR-CI calculation
C
         CALL LUCI_KRCI(CMO,WORK,LWORK)
      ELSE
         WRITE(LUPRI,'(A,I18)') ' *** ERROR in LUCI_RCIST ***'//
     &   ' wrong KR-run type in LUCI_RCIST !',JKRRUNTYPE
         CALL QUIT('*** ERROR in LUCI_RCIST ***')
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_KRCI(CMO,WORK,LWORK)
C***********************************************************************
C
C     Run large-scale CI calculations. Called from module KR-CI.
C
C     Based on LUCI_CIINII.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
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
      DIMENSION WORK(*),CMO(*)
C
#include "../luciarel/clunit.inc"
#include "parluci.h"
C
#include "memint.h"
C
      CIRUN = 'KR-CI '
      MZ = MIN(NZ,2)
C
      KFRSAV = KFREE
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
#if defined (VAR_MPI2)
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
          WRITE(LUPRI,'(A,1F18.10)') '   Core energy : ',ECORE
        END IF
C
C       Get active-active part of FCmo:
C
        CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A)') 'Active Fock matrix:'
          CALL WRTMATMN(WORK(KFCAC),1,N2ASHXQ,1,N2ASHXQ,LUPRI)
        END IF
C
C       Transform quaternion FCACmo to molfdir type FCAC.
C
        CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A,I6,I6)') 'Molfdir format FCAC, dimension:',
     &                        (2*NASHT)*(2*NASHT),NASHT
          WRITE(LUPRI,'(A)') 'Real part:'
          CALL WRTMATMN(WORK(KFCACM),1,N_FCACM,1,N_FCACM,LUPRI)
          IF (NZ.GT.1) THEN
            WRITE(LUPRI,'(A)') 'Imaginary part:'
            CALL WRTMATMN(WORK(KFCACM+N_FCACM),1,
     &                    N_FCACM,1,N_FCACM,LUPRI)
          END IF
        END IF
C
C       WORK(KFCACM) now contains all 1-e integrals
C
        CALL MEMREL('after RGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
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
        CALL QUIT('LUCI_KRCI: Invalid flag I_IITP.')
      end if
#if defined (VAR_MPI2)
C
C     end of slave reading test
C
      IF( ISET_PARCAL .eq. 1 ) THEN
        CALL TEST_PAR_CALC(PARCAL_SAVE,1)
      END IF
#endif
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call LUCIAREL
      CALL KRMC_LUCI(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,ECORE,
     &               CIRUN)
C
      CALL MEMREL('LUCI_KRCI',WORK,1,KFRSAV,KFREE,LFREE)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_CIINII(CREF,CMO,WORK,LWORK)
C***********************************************************************
C
C     Get CI start guess for KR-MCSCF optimization.
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
#include "parluci.h"
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
#if defined (VAR_MPI2)
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
          WRITE(LUPRI,'(A,1F18.10)') '   Core energy : ',ECORE
        END IF
C
C       Get active-active part of FCmo:
C
        CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A)') 'Active Fock matrix:'
          CALL WRTMATMN(WORK(KFCAC),1,N2ASHXQ,1,N2ASHXQ,LUPRI)
        END IF
C
C       Transform quaternion FCACmo to molfdir type FCAC.
C
        CALL QFC2MFC(WORK(KFCAC),WORK(KFCACM),1,1,IPROPT)
C
        IF (IPROPT.GE.5) THEN
          WRITE(LUPRI,'(A,I6,I6)') 'Molfdir format FCAC, dimension:',
     &                        (2*NASHT)*(2*NASHT),NASHT
          WRITE(LUPRI,'(A)') 'Real part:'
          CALL WRTMATMN(WORK(KFCACM),1,N_FCACM,1,N_FCACM,LUPRI)
          IF (NZ.GT.1) THEN
            WRITE(LUPRI,'(A)') 'Imaginary part:'
            CALL WRTMATMN(WORK(KFCACM+N_FCACM),1,
     &                    N_FCACM,1,N_FCACM,LUPRI)
          END IF
        END IF
C
C       WORK(KFCACM) now contains all 1-e integrals
C
        CALL MEMREL('after RGETFC',WORK,1,KFCAC,KFREE,LFREE)
C
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
        CALL QUIT('LUCI_CIINII: Invalid flag I_IITP.')
      end if
#if defined (VAR_MPI2)
C
C     end of slave reading test
C
      IF( ISET_PARCAL .eq. 1 ) THEN
        CALL TEST_PAR_CALC(PARCAL_SAVE,1)
      END IF
#endif
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call information transfer routine
C FIXME: Direct handling hardwired for the moment. Later, this
C        should be made dynamical and controlled from the outside.
C        Possibly also integrals reside on file in LUCIAREL
      IO_RCIST_IN = 0
C
C     SK - 03-04-08: analyzing converged MCSCF vector with LUCIAREL
C                    set MOPT_MXMACRO to 0 and MAXCIT to 0
C
      IF( MOPT_MXMACRO .eq. 0 .and. MAXCIT .eq. 0 ) IO_RCIST_IN = 1
      IF (IO_RCIST_IN.EQ.1) THEN
CC      CALL KR_LUCI_IO(0,0,1,DUMMY,1,DUMMY,1,
CC   &                  WORK(KFCACM),N2ASHXQ,WORK(KH2AC),
CC   &                  (2*NASHT)**4*2,
CC   &                  DUMMY,DUMMY,2,IRIQ,
CC   &                  IDUMMY,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        REWIND LUC
        CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,2,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
C       release memory
        CALL MEMREL('LUCI_SIG 1',WORK,1,KRSCR,KFREE,LFREE)
      END IF
C
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
          CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
          CALL KR_LUCI_IO(1,0,0,CREF,NZCONF*MZ,DUMMY,1,
     &                    DUMMY,1,DUMMY,1,
     &                    WORK(KRSCR),DUMMY,1,MZ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        END IF
      END IF
C
      CALL MEMREL('LUCI_CIINII',WORK,1,KFRSAV,KFREE,LFREE)
C
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LUCI_IJKLRO(CMO,WORK,LWORK)
C***********************************************************************
C
C     Do integral resorting and write to file
C    
C     !!!!!!!!!!!!!!!!!!!!!!  GAS-scheme specific !!!!!!!!!!!!!!!!!!!!!!
C
C     Written by S. Knecht                                - Jan 09 2008
C     Last revision :                                                  
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CHARACTER*6 CIRUN
C
      DIMENSION WORK(*), CMO(*)
C
#include "parluci.h"
C
#include "memint.h"
C
      CIRUN = 'IJKLRO'
      MZ = MIN(NZ,2)
C
      KFRSAV = KFREE
C
      N_FCACM = (2*NASHT)*(2*NASHT)
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
      IF( CRDFO_MAT ) THEN
C
C       read FCmo from file LUKRM3
C
        CALL KRMC_GETFCK(WORK(KFC),IPROPT,WORK(KFREE),LFREE)
C
        ECORE = ECORE_LR
C
      ELSE
C
C       Calculate FCmo ( CWRTFO_MAT false: do not write to LUKRM3) :
C
        CALL RGETFC(CWRTFO_MAT,CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
      END IF
      IF (IPROPT.GE.1) THEN
        WRITE(LUPRI,'(A,1F18.10)') 'Core energy : ',ECORE
      END IF
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
      IF (IPROPT.GE.5) THEN
        WRITE(LUPRI,'(A)') 'Active Fock matrix:'
        DO I=1,N2ASHXQ,1
          WRITE(6,'(1I3,1F18.14)') I,WORK(KFCAC+I-1)
        END DO
      END IF
C
C     Transform quaternion FCACmo to molfdir type FCAC.
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
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     read two-electron integrals in luciarel/transform_r.F 
C     (subroutine itrctl) 
C
      LEN_ALL_INT = 0
      LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
C
C     Transfer data to LUCIAREL common blocks.
      CALL TRKRLUCI(CIRUN,0)
C
C     Setup routines for running LUCIAREL
      CALL LUCI_SETUP()
C
C     Call LUCIAREL
      CALL KRMC_LUCI(WORK(KFCACM),DUMMY,WORK(KFREE),LFREE,DUMMY,
     &               CIRUN)
C
C     close files ...
C
      CALL LUCI_CLOSEDOWN()
C
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
      CALL TRKRLUCI(CIRUN,0)
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
        CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,CB,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,2,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
C       release memory
        CALL MEMREL('LUCI_SIG 1',WORK,1,KRSCR,KFREE,LFREE)
      END IF
CSK      IPROPT = 0
C
C     Call LUCIAREL
      CALL KRMC_LUCI(FCACM,H2ACM,WORK(KFREE),LFREE,DUMMY,CIRUN)
C
C     Read sigma vector on file to SB array.
      IO_SIGMA_OUT = 1
      IF (IO_SIGMA_OUT.eq.1) THEN
        CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,SB,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,1,MZ,
     &                  LUHC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        CALL MEMREL('LUCI_SIG 2',WORK,1,KRSCR,KFREE,LFREE)
      END IF
      CALL REWINE(LUC,-1)
      IO_C_VEC_OUT = 1
      IF (IO_C_VEC_OUT .eq. 1 ) THEN
        CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,CB,NZCONF*MZ,
     &                  DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),IDUMMY,1,MZ,
     &                  LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
        CALL MEMREL('LUCI_SIG 3',WORK,1,KRSCR,KFREE,LFREE)
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
        CALL MEMREL('LUCI_SIG D',WORK,1,KDIAG,KFREE,LFREE)
      END IF
CSK      IPROPT = 0
      CALL MEMREL('LUCI_SIGMA',WORK,1,KFRSAV,KFREE,LFREE)
C
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
#include "parluci.h"
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
      CALL TRKRLUCI(CIRUN,0)
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
      LF2_ZERO = L2
C
      IF( .NOT. NATOLCR ) THEN
C       
C       Write current CI vector(s) to file
        IO_DENS_IN = 1
        IF (IO_DENS_IN.EQ.1) THEN
          CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
          CALL KR_LUCI_IO(1,1,0,CREFL,NZCONF*MZ,
     &                          CREFR,NZCONF*MZ,
     &                    DUMMY,1,DUMMY,1,
     &                    WORK(KRSCR),IDUMMY,2,MZ,
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
csk   WRITE(LUPRI,*) ' THIS IS MY RHO1'
csk   CALL WRTMATMN(RHO1,1,L1,1,L1,LUPRI)
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
      SUBROUTINE LUCI_PROP(WORK,KFREE,LFREE)
C***********************************************************************
C
C     CIRUN = 'PROP1 ' : compute one-particle properties using 
C                        1e-property integrals residing on TPROP
C
C     INPUT:           
C           CREFR            : Right CI vector (expected to reside on
C                              file LUCVECT_REL.0)
C       general also for transition densities
C
C     Written by S. Knecht                 - Jul 22 2008
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
      DIMENSION WORK(*)
C
#include "../luciarel/clunit.inc"
#include "parluci.h"
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
      CALL TRKRLUCI(CIRUN,0)
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
      LF2_ZERO = L2
C
      IF( .NOT. NATOLCR ) THEN
C       
C       Write current CI vector(s) to file
        IO_DENS_IN = 1
        IF (IO_DENS_IN.EQ.1) THEN
          CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
          CALL KR_LUCI_IO(1,1,0,CREFL,NZCONF*MZ,
     &                          CREFR,NZCONF*MZ,
     &                    DUMMY,1,DUMMY,1,
     &                    WORK(KRSCR),IDUMMY,2,MZ,
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
csk   WRITE(LUPRI,*) ' THIS IS MY RHO1'
csk   CALL WRTMATMN(RHO1,1,L1,1,L1,LUPRI)
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
      CALL TRKRLUCI(CIRUN,0)
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
     &                    WORK(KASCR),IDUMMY,2,MZ,
     &                    LUC,IDUMMY,IDUMMY,IDUMMY,IPROPT)
          CALL MEMREL('LUCI_ANALYZ1',WORK,1,KASCR,KFREE,LFREE)
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
      KFRSAV = KFREE
C
      MZ = MIN(NZ,2)
C
C     Introduce LUCIAREL action
      CALL LUCI_TASK(CIRUN)
C
C     Transfer data to LUCIAREL common blocks.
C     FIXME: Check this here.
      CALL TRKRLUCI(CIRUN,0)
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
        CALL MEMGET('REAL',KRSCR,NZCONF*MZ,WORK,KFREE,LFREE)
        CALL KR_LUCI_IO(1,0,0,DIAG,NZCONF*MZ,DUMMY,1,
     &                  DUMMY,1,DUMMY,1,
     &                  WORK(KRSCR),DUMMY,1,MZ,
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
      CALL TRKRLUCI(CIRUN,0)
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
#include "parluci.h"
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
#include "parluci.h"
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
#if defined (VAR_MPI2)
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
#include "parluci.h"
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
     &    CIRUN .eq. 'SIGMAD' .or. CIRUN .eq. 'KR-CI ' ) THEN
C
C       read-in of integrals by nodes... no broadcast in SYNC_NO...P2
C
        IF( (CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' ) 
     &      .and. IIOMOD_REL .eq. 1 )THEN
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
        IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' )THEN
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
      ELSE IF( CIRUN .eq. 'IJKLRO' ) THEN 
C
        N_FCACM = (2*NASHT)*(2*NASHT)
        LEN_ALL_INT = NASHT*NASHT*NNASHX*NZ*3
C
        CALL MEMGET('REAL',KFCACM,MZ*N_FCACM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KH2AC,          0,WORK,KFREE,LFREE)
C
        CALL DZERO(WORK(KFCACM),MZ*N_FCACM)
C
        CALL SYNC_NODES_REL_P2(WORK(KFCACM),WORK(KH2AC),CIRUN)
C
      END IF
C
C     Transfer data to LUCIAREL common blocks.
C
      CALL TRKRLUCI(CIRUN,1)
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
C     initialize common block
C
      CALL MPIXCALC_TYPESZ()
C
C     Run LUCIAREL
C
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'SIGMA ' .or. 
     &    CIRUN .eq. 'SIGMAD' .or. CIRUN .eq. 'KR-CI ' ) THEN
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
      ELSE IF( CIRUN .eq. 'IJKLRO' ) THEN 
C
C       Call LUCIAREL for integral reordering.
        CALL KRMC_LUCI2(WORK(KFCACM),WORK(KH2AC),WORK(KFREE),LFREE,
     &                  DUMMY,CIRUN)
C
        CALL MEMREL('LUCINOD_RL',WORK,1,KFRSAV,KFREE,LFREE)
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
#endif    /* ifdef VAR_MPI2 */

