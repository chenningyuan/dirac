*
#ifdef LUCIAREL_STANDALONE
      PROGRAM LUCIAREL
#else
      SUBROUTINE PAMLUCIAR(WRK,LWRK)
#endif
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
* General relativistic formalism for 4-/2-component MCSCF
*                          Timo Fleig + Jeppe Olsen
*                          2000 - 2001
*
* Excitation operator formalism for sigma and density implemented.
*                          Timo Fleig + Jeppe Olsen
*                          2001 - 2002
*
* Implementation of relativistic double group CI,
*                          Timo Fleig + Jeppe Olsen
*                          1995 - 1999
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Parameters for dimensioning
#include "mxpdim.inc"
*.Memory
#include "wrkspc.inc"
*.File numbers
#include "clunit.inc"
*.Print flags
#include "cprnt.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "crun.inc"
#include "cicisp.inc"
#include "oper.inc"
#include "cgas.inc"
#include "glbbas.inc"
#include "csm.inc"
#include "cintfo.inc"
#include "picasso.inc"
*
      common/glbref/ISPC,ISM,ITERGLB
#ifdef LUCIAREL_STANDALONE
*.
********************************************
********************************************
* Initialize BIGGY
      do I=1,MXPWRD,1
         WORK(I) = 0.0D0
      end do
#else
      DIMENSION WRK(LWRK)
C
C     Compute offset from WRK to WORK
C
      CALL COMPOFF(WRK,WORK,KBASE)
C
      MXPWRD = LWRK
      KADD = MXPWRD
      CALL MEMMAN(KBASE,KADD,'INI   ',IDUMMY,'DUMMY')
C
#endif
********************************************
********************************************
      CALL QENTER('LUCIAR')
*.    No floating point underflow
      CALL XUFLOW
*. Assign diskunits
      CALL DISKUN
*. Header
      CALL WRTHD(LUOUT)
      WRITE(LUOUT,'(10X,A,I12)')
     &'Dimension of R*8 workspace : ', MXPWRD
* ======================================================
* 1 : .Read input,insert defaults and cross check input
* ======================================================
#ifdef LUCIAREL_STANDALONE
      CALL INPCTL(LUIN,LUOUT)
#else
      OPEN(LUIN,FILE='DIRAC.INP')
      CALL INPCTL(LUIN,LUOUT)
      CLOSE(LUIN, STATUS='KEEP')
#endif
*. Ratio beteeen real and integer word length
      CALL ZIRAT
*. From shells to orbitals
      CALL ORBINF(LUOUT,IPRORB)
*. Types of determinants to be included
      CALL DETTYP(IPRSTR)
      print*,'ISO1,ISO2 ',ISO1,ISO2
*. Number of string types
      CALL STRTYP_GAS(IPRSTR)
*. Symmetry information
      CALL SYMINF(IPRORB)
*. Number of integrals
      CALL INTDIM(IPRTRA)
*. Static memory, initialization and allocation
#ifdef LUCIAREL_STANDALONE
      KBASE = 1
      KADD = MXPWRD
      CALL MEMMAN(KBASE,KADD,'INI   ',IDUMMY,'DUMMY')
C?    CALL MEMMAN(IDUM,IDUM,'PRINT ',IDUM,'IDUM')
#endif
*. Allocate work space, general
      CALL ALLOC
*. Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS(WORK(1),IPRSTR)
*. Internal subspaces
      CALL LCISPC(IPRCIX)
*. Space for Hamiltonian in form suitable for SIGDEN
*. Type of Hamiltonian : 2 e Delta MK = 0 + all 1 e ops
      if (SIGDEN_ROUTE.eq.'SIGDEN') then
        if (ISO1.eq.0.and.ISO2.eq.0) then
          IHTYPE = 1
        else if (ISO1.eq.1.and.ISO2.eq.0) then
          IHTYPE = 2
        else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).ne.'DIRAC ') 
     &       then
          IHTYPE = 3
        else if (ISO1.eq.1.and.ISO2.eq.1.and.ENVIRO(1:6).eq.'DIRAC ') 
     &       then
          IHTYPE = 4
        else
          write(6,*) 'Setting of ISO1,ISO2 not allowed.'
          write(6,*) 'ISO1,ISO2 : ',ISO1,ISO2
          call abend2('Quitting.')
        end if
        ISYM_T = 1
        call set_hop_dbg(IHTYPE,IPRHAM)
      END IF
*. Symmetry of reference
      IF(PNTGRP.GT.1.and.PNTGRP.lt.5)
     &           CALL MLSM(IREFSM,IREFPA,IREFSM,'CI',1)
      IF(NOINT.EQ.1) THEN
        WRITE(6,*) ' End of calculation without integrals'
#ifndef LUCIAREL_STANDALONE
        RETURN
#endif
        CALL QSTAT
        Call Abend2( ' End of calculation without integrals' )
      END IF
*
*====================================================
*  Iterative natural spinor CI method
*  ----------------------------------
*   (Default is one iteration, so normal CI)
*====================================================
*
      do INATIT=1,NATITER,1
*
*  Generate symmetry info and the like
        if (INATIT.eq.1) call luc2pic(ISO1,IPRTRA)
*
*. Generate and store pointers for transformed integrals
*  (and integrals !)
        IF(NOINT.EQ.0.and.INATIT.eq.1) THEN
          CALL INTIM(IPRTRA)
        ELSE
          WRITE(6,*)
     &      ' No integrals imported (or > 1st NS iteration) '
        END IF
*
*. Now call PICASSO
*. PICASSO will read the integrals (in initial run)
*.              read SPAO coefficients (in initial run)
*.              and do the transformation
*.              of all integrals
*
        call picasso(INTSPACE,ISO1,IPRTRA,
     &               NL2D,IUB2D,ISM2D)
        write(6,*) 'Integrals transformed and/or distributed.'
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
          write(6,*) 'THIS ROUTE IS STILL FORBIDDEN !!'
          write(6,'(A16,I2,1F3.1)')
     &          'IRESTR,XLAMBDA =',IRESTR,XLAMBDA
          Call Abend2( 'Quitting.' )
*
        END IF
*       ^ End of special handling of restarted calc with lambda modified op.
*
        ISKIPEI_INI = ISKIPEI
        DO JCMBSPC = 1, NCMBSPC
          ISKIPEI = 0
*
*. Perform GASCI/GASPT calculation for internal space ICISPC
*
          WRITE(6,*)
          WRITE(6,'(A)')
     &    ' =============================================='
          WRITE(6,'(A,I3)')
     &    ' Information on calculations in space ', JCMBSPC
          WRITE(6,'(A)')
     &    ' =============================================='
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &    ' Number of calculation in this CI space ', NSEQCI(JCMBSPC)
*
          IF(IRESTR.EQ.1.AND.XLAMBDA.NE.1.0D0.AND.JCMBSPC.EQ.1) THEN
            WRITE(6,*) ' Remember No restart in calc 1 (Lambda calc)'
            IRESTR = 0
            Call Abend2( 'Quitting.' )
          END IF
*
          I_EXPAND = 1
          IF(XLAMBDA.NE.1.0D0 .AND.JCMBSPC.GT.1) THEN
            WRITE(6,*) ' =================================='
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) '   Modified operator will be used'
            WRITE(6,*) ' =================================='
            IF(JCMBSPC.EQ.2) THEN
*. Modify the operator
              WRITE(6,*) ' Operator will be modified '
              WRITE(6,*) ' This route is still forbidden! '
              Call Abend2( 'Quitting.' )
CTF           CALL GENH1(XLAMBDA)
CTF           CALL SCLH2(XLAMBDA)
            END IF
          END IF
          IF(JCMBSPC.EQ.2.AND.IRESTR_ORIG.EQ.1
     &                   .AND.XLAMBDA.NE.1.0D0)  THEN
            WRITE(6,*) ' Restart vectors will be copied to LUC'
            WRITE(6,*) ' CI will restart with vectors from prev. calc'
            IF(ICISTR.EQ.1) THEN
              LBLK = XISPSM(IREFSM,2)
            ELSE
              LBLK = -1
            END IF
            CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'VCSAVE')
            IF(ICISTR.EQ.1) THEN
              LBLOCK = XISPSM(IREFSM,2)
            ELSE IF (ICISTR.EQ.2) THEN
              LBLOCK = MXSB
            ELSE IF (ICISTR.EQ.3) THEN
              LBLOCK = MXSOOB
            END IF
            CALL MEMMAN(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
            CALL REWINO(LUC)
            CALL REWINO(LUSC39)
            DO JROOT = 1, NROOT
C             CALL COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
              CALL COPVCD(LUSC39,LUC,WORK(KVEC1),0,LBLK)
            END DO
C?        WRITE(6,*) ' Back copied first vector on LUC'
C?        CALL WRTVCD(WORK(KVEC1),LUC,1,LBLK)
            CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'VCSAVE')
*. No expansion should follow
            I_EXPAND = 0
*. But normal restart
            IRESTR = 1
          END IF
*
          DO JSEQ = 1,NSEQCI(JCMBSPC)
            IF((ISEQCI(JSEQ,JCMBSPC).GT.0
     &          .OR.ISEQCI(JSEQ,JCMBSPC).EQ.-5)
     &          .AND. JCMBSPC.NE.1.AND.ISTOSPC.NE.JCMBSPC
     &          .AND.I_EXPAND.EQ.1)                   THEN
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
              WRITE(6,*) ' LBLK = ', LBLK
              CALL EXPCIV(IREFSM,ISTOSPC,LUIN,JCMBSPC,LUOUT,
     &                    LBLK,LUSC2,
     &                    NROOT,1,IDC,IPRDIA)
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
                WRITE(6,*) ' No restart from previous vectors'
              END IF
*. Good old normal CI !!!!
*. do CI in space JCMBPSC
              MAXIT = ISEQCI(JSEQ,JCMBSPC)
*
              ISPC = JCMBSPC
              ISM = IREFSM
              ITERGLB = 1
*
*
*
              CALL GASCI(INTSPACE,IREFSM,JCMBSPC,IPRDIA,EREF)
*
*
*
*. Last space where vectors were stored
              ISTOSPC = JCMBSPC
CC            ELSE IF (ISEQCI(JSEQ,JCMBSPC).GT. -5 ) THEN
CC              WRITE(6,'(A,I3)') ' Vector free calculation at level ',
CC     &                            -ISEQCI(JSEQ,JCMBSPC)
CC              LEVEL =  -ISEQCI(JSEQ,JCMBSPC)
CC              IF(IC1DSC.LE.0) THEN
CC                LU1EFF = 0
CC              ELSE
CC                LU1EFF = LUHC
CC              END IF
CC              WRITE(6,*) ' LU1EFF in MAIN ', LU1EFF
CC              CALL DIRDIR(JCMBSPC-1,JCMBSPC,IREFSM,LUC,LEVEL,EREF,
CC     &                    LUSC2,LUSC3,LU1EFF)
CC            ELSE IF (ISEQCI(JSEQ,JCMBSPC).LE.-5) THEN
CC              WRITE(6,'(A)') ' Perturbation calculation '
CC              CALL PERTCTL(IREFSM,JCMBSPC,EREF)
*. also perturbation expansion of EKT ??
CC              IF(IPTEKT.EQ.1) THEN
CC                WRITE(6,*) ' Perturbation expansion of EKT '
CC                WRITE(6,*) ' ============================= '
C                    PTEKT(LUN,LU0,N,NEKT,LPTEKT,ISM,ISPC)
CC                CALL PTEKT(LUC,LUSC36,NPERT,IREFSM,JCMBSPC)
CC              END IF
            END IF
          END DO
*
        END DO
*       ^ End loop over CI spaces
*
        if (NOMOFL.eq.0.or.IFINMO.eq.6) THEN
*. Create set of spinors from densities
c         WRITE(6,*) ' I am going to call FINMO '
c         CALL FINMO_REL
        end if
      end do
*     ^ End loop over natural spinor iterations (default=1)
*
      CALL QEXIT('LUCIAR')
      CALL QSTAT
#ifdef LUCIAREL_STANDALONE
      Call Abend2( ' I am home from  the loops ' )
#endif
      END
