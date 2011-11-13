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
C VERSION : $Revision: 1.38.2.3 $
C DATE    : $Date: 2004/09/01 08:04:09 $
C FILE    : herpar.F
C
C
#ifdef UNDEF
/* Comdeck notes */
961125 - tec RELCAL: FIRST replaced by COMMON parfir remember to set FIRST1
961003 - tsaue: DOINT is removed and replaced by I2TYP.
                Some polishing av code:
                  Note that NPOS enters PARDRV as argument
                  and that Fock matrices are not initialized at master !
950620-kr: Ported to MPI, based on previous efforts by Bjarne Herland at
           Para//ab
  --- in PARLOP is checked for DOINT(*,*,1) on the IJSHEL;
      this will call NULL proccesses on nodes if some DOINT values are
      false. It would be better to initialize with DOINT in INDEKS,
      but that requires a little more coding.
      SORT renamed to PARSRT (it is likely that library routines exist
      called SORT).
950531-kr: Reorganized the code so that the node programs are
started only once for a given calculation
950428-kr+hjaaj;
 Merged with Odense version. Implemented DOINT in PARLOP.
941219-hjaaj
(1)  Note: in old HERFCK 'IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD' was
 transferred to TWOINT, do we want to do that in CALINT for ITYPE=3 ???
(2)  removed JSYMDM,JFCTYP,ITRIPL from /PARINT/ in aba.cdk
941114-hjaaj
PARDRV,SDINIT,RVINIT,CALINT:  removed TRIPLE in parameter list
SDINIT,RVINIT: removed ITRIPL message tag; moved IFCTYP and IREPDM
 message tags after NDMAT message tag
CALINT: call TWOINT directly instead of CALL HERFCK
PARDRV,...: s/IFCTYP/IFCTYP(*)/ and s/IREPDM/IREPDM(*)/
#endif
C  /* Deck pardrv */
      SUBROUTINE HER_PARDRV(WORK,LWORK,FMAT,DMAT,NDMAT,IREPDM,IFCTYP,
     &                  ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,
     &                  TKTIME,RETUR,FIRST,NPOS,NTOTTK,I2TYP,
     &                  ICEDIF,SCREEN,GABRAO,DMRAO,DMRSO,DINTSKP)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993 -> Okt-1994
C     Extended with I2TYP Oct 1 1996 by T.Saue
C     Extended with screening Oct 11 1996 by T.Saue
C
C
C     **************************************************************
C     *        Driver routine for the parallel calculation.        *
C     **************************************************************
C
C
C      A short survey of the different messagetags (MSGTAG) :
C
C      10 - Send initialization to nodes.
C      20 - Receive request from node.
C      30 - Send batch to node.
C      40 - Receive final results from nodes.
C      50 - Get timing results from nodes.
C
C
C           ITYPE  - Calc. type: 2 = expectation values.
C                                3 = fock matrices.
C                               -5 = Integrals derivated with respect
C                                    to magnetic field.
C                                8 = Fock matrices derivated with respect
C                                    to atomic coordinates.
C           NTASK  - Number of tasks
C           IPRPAR - Print level during parallelization.
C           MAXDIF -
C           MAXREP -
C           IREPDM - Symmetries of density matrices
C           IFCTYP - Fock matrix types (see twoint).
C           IATOM  -
C           NODV   - Density matrix neglected in TWOEXP
C           NOPV   -           -   "   -
C           NOCONT -           -   "   -
C           RETUR  - Program will exit after spec. shells
C           TKTIME - Take time in TWOINT
C           TIMING - Take time for each integralbatch IJ
C           DOREPS -
C           DOSYM  -
C           DOCOOR -
C           NMLINE - Number of lines in MOLECULE.INP
C           MLINE  - MOLECULE.INP as internal file
C           NDMAT  - Number of density-matrices
C           DMAT   - Density-matrices.
C
C
#include "implicit.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "molinp.h"
#include "twocom.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "inforb.h"
#include "infpar.h"
C
      DIMENSION FMAT(*), DMAT(N2BASX,NDMAT), WORK(LWORK),
     &          NPOS(NTOTTK), IREPDM(*), IFCTYP(*),
     &          GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
      LOGICAL FIRST,NODV,NOPV,NOCONT,RETUR,TKTIME
C
      CALL QENTER('HER_PARDRV')
C
#if !defined (VAR_MPI)
      CALL QUIT
     &('program error: parallel code called but MPI not activated')
#endif
C
C
C     Ensure that NODES .LE. MAXNOD
C
      IF (NODES .GT. MAXNOD) THEN
         WRITE(LUPRI,'(/5X,A,/12X,A,/5X,A16,I8,/5X,A16,I8)')
     &    'ERROR: Number of nodes allocated has exceeded hardcoded',
     &    'limit. Reduce number of nodes or recompile program. ',
     &    'Allocated (NODES) :',NODES,'Limit (MAXNOD) :',MAXNOD
         CALL QUIT('Number of nodes requested/allocated is too large.')
      ELSE IF (IPRPAR .GT. 0) THEN
         WRITE(LUPRI,'(/5X,A,I5/)')
     &    'HER_PARDRV: Number of nodes allocated:', NODES
      ENDIF
      IF(ITYPE.EQ.7) THEN
        KWORK = 1
        KFREE = KWORK
        LFREE = LWORK
        TIMING = .FALSE.
        GOTO 10
      ENDIF
C
C     Take calculation time of each IJ pair?
C
      TIMING = .TRUE.
C
C     Determine total number of tasks - MTOTTK.
      CALL DET_TASK (I2TYP,1,MTOTTK,NPOS)
      IF(MTOTTK.GT.NTOTTK) THEN
        WRITE(LUPRI,'(A,I8)') 'HER_PARDRV: Inconsistent MTOTTK:',MTOTTK
        WRITE(LUPRI,'(3X,A,I8,A,I8)') ' - should be less or equal to',
     &   NTOTTK,' for I2TYP =',I2TYP
        CALL QUIT('HER_PARDRV: Inconsistent MTOTTK !!!')
      ENDIF
C
C     Allocate local variables
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
      CALL MEMGET('INTE',KNSTAT,NODES ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,MTOTTK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTIMES,MTOTTK,WORK,KFREE,LFREE)
C
C     Test of work space allocation.
C
      LDINTSKP = 8
      IF (ITYPE .EQ. 2) THEN
         NUMGRA = MXCOOR
         LDINTSKP = 3*3*14
         IF (MAXDIF .GT. 1) THEN
            NUMHES = NUMGRA*NUMGRA
            NREQUI = NUMGRA + NUMHES
         ELSE
            NUMHES = 0
            NREQUI = NUMGRA
         END IF
      ELSE IF (ITYPE .EQ. 3) THEN
         NREQUI = NDMAT*N2BASX
      ELSE IF (ITYPE .EQ. 9) THEN
         NREQUI = NDMAT*N2BASX
      ELSE IF (ITYPE .EQ. -5) THEN
         NUMSUS = 9
         NREQUI = 3*NNBASX + NUMSUS
      ELSE IF (ITYPE .EQ. 8) THEN
         NREQUI = 3*NNBASX
      END IF
C
C     Initialize DINTSKP
C
      CALL DZERO(DINTSKP,LDINTSKP)
C
C
      IF (LFREE .LT. NREQUI) THEN
         WRITE(LUPRI,'(/5X,A,/5X,A16,I8,/5X,A16,I8)')
     &    'Insufficient work space for subroutine RVRES',
     &    'Free space     :',LFREE,'Required space :',NREQUI
         CALL QUIT('Insufficient work space for RVRES')
      ENDIF
C
C--------------------------
C     Make the index array.
C--------------------------
C
      CALL INDEKS(FIRST,WORK(KINDEX),ITYPE,MTOTTK,NBATCH,NPOS,I2TYP)
C
C---------------------------------------
C     Send initialization data to nodes.
C---------------------------------------
C
 10   CONTINUE
      CALL SDINIT(DMAT,NDMAT,IREPDM,IFCTYP,ITYPE,MAXDIF,IATOM,NODV,
     &            NOPV,NOCONT,TKTIME,RETUR,I2TYP,ICEDIF,SCREEN,
     &            GABRAO,DMRAO,DMRSO)
C     and transfer screening information to twocom.h
C     (such that master will print correct scr.thr. in summary) /hjaaj apr 01
      SCRTHR = SCREEN
C
C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DOPAR(WORK(KINDEX),MTOTTK,NBATCH,WORK(KNSTAT))
      CALL RVRES(WORK(KFREE),FMAT,ITYPE,MTOTTK,NREQUI,
     &           NUMGRA,NUMHES,MAXDIF,WORK(KTIMES),DINTSKP,LDINTSKP)
C
C--------------------------------------
C     Receive final results from nodes.
C     Release nodes from this task
C      and let them return to main menu.
C--------------------------------------
C
      CALL DIRAC_PARCTL( -1 )
C
C-----------------------------------------
C     Sort the integral-calculation times.
C-----------------------------------------
C
      IF(TIMING) CALL PARSR2(MTOTTK,NPOS,WORK(KTIMES),
     &                       WORK(KINDEX),IPRPAR)
C
C-------------------------------------
C     Get timing results if requested.
C-------------------------------------
C
      IF (IPRPAR .GT. 0) CALL TIMRES(WORK(KNSTAT))
C
      IF (FIRST)  FIRST  = .FALSE.
      IF (NEWBAS) NEWBAS = .FALSE.
C
 90   CONTINUE
      CALL MEMREL('HER_PARDRV',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('HER_PARDRV')
      RETURN
      END
C  /* Deck indeks */
      SUBROUTINE INDEKS(FIRST,INDEX,ITYPE,MTOTTK,NBATCH,
     &                  NPOS,I2TYP)
C
C     Written by Paal Dahle May 1994.
C     Revised Oct 1996 by T.Saue
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "infpar.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "blocks.h"
#include "dcbfir.h"
C
      DIMENSION INDEX(MTOTTK), NPOS(MTOTTK)
      LOGICAL FIRST,DOFIRST
C
C     If new run, make initial guess as to run times
C
      DOFIRST = .FALSE.
C
      IF (RELCAL) THEN
         IF (I2TYP.EQ.1.AND.FIRST1) THEN
            DOFIRST = .TRUE.
            FIRST1 = .FALSE.
         ELSEIF (I2TYP.EQ.2.AND.FIRST2) THEN
            DOFIRST = .TRUE.
            FIRST2 = .FALSE.
         ELSEIF (I2TYP.EQ.3.AND.FIRST3) THEN
            DOFIRST = .TRUE.
            FIRST3 = .FALSE.
         ELSE
         ENDIF
      ENDIF
C
      IF (NEWBAS .OR. DOFIRST) THEN
         CALL DET_TASK (I2TYP,2,MTOTTK,NPOS)
      ENDIF
C
C     Distribute jobs on nodes
C
      ICOUNT = 0
      IBLOCK = NODES*NTASK
      DO I = MTOTTK, 1, -1
        MORE1  = (ICOUNT/IBLOCK)*IBLOCK
        MORE2  = MOD(ICOUNT,NODES)*NTASK
        MORE3  = ((ICOUNT-MORE1)/NODES) + 1
        IPLACE = MORE1+MORE2+MORE3
        INDEX(IPLACE) = NPOS(I)
        ICOUNT=ICOUNT+1
      ENDDO
C
C     Determine number of batches (NBATCH).
C
      NOTALL = (MTOTTK/IBLOCK)*IBLOCK
      NREST  = MTOTTK - NOTALL
      IF (NREST .GT. NODES) NREST = NODES
      NBATCH = NOTALL/NTASK + NREST
C
C     Print current settings.
C
      NEACH = NBATCH/NODES
      NREST = NBATCH - NEACH*NODES
C
      IF (FIRST .AND. IPRPAR .GT. 4) THEN
         WRITE(LUPRI,'(5(/5X,A,I5),2X,A7,I3,A1/)')
     &        'Number of nodes   :',NODES,
     &        'Number of tasks   :',MTOTTK,
     &        'Tasks pr. batch   :',NTASK,
     &        'Number of batches :',NBATCH,
     &        'Batches pr. node  :',NEACH,
     &        '(rest =',NREST,')'
      END IF
C
C     Print the index array.
C
      IF (IPRPAR .GE. 10) THEN
         WRITE(LUPRI,'(/5X,A/)') 'The order of task-distribution :'
         DO 400 I=1, NBATCH*NTASK
            WRITE(LUPRI,'(5X,A3,I5,A15,I5)')
     &           'i =',I,'index(i) =',INDEX(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck sdinit */
      SUBROUTINE SDINIT(DMAT,NDMAT,IREPDM,IFCTYP,ITYPE,MAXDIF,IATOM,
     &                  NODV,NOPV,NOCONT,TKTIME,RETUR,I2TYP,
     &                  ICEDIF,SCREEN,GABRAO,DMRAO,DMRSO)
C
C     Written by Henrik Koch and Paal Dahle
C
#include "implicit.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
#endif
#if defined (MOD_DFT)
#include "dftcom.h"
#include "dfterg.h"
#endif
#include "molinp.h"
#include "aovec.h"
C
#include "hrunit.h"
#include "dorps.h"
#include "abainf.h"
#include "blocks.h"
#include "symmet.h"
#include "nuclei.h"
#include "infpar.h"
#include "parint.h"
C
#include "inforb.h"
#include "cbisol.h"
#include "cbirea.h"
#include "siripc.h"
C
      LOGICAL NODMAT,NOFMAT
      DIMENSION DMAT(*), IREPDM(NDMAT), IFCTYP(NDMAT),
     &          DMRAO(NSYMBL,NSYMBL,*),GABRAO(NSYMBL,NSYMBL,*),
     &          DMRSO(MAXSHL,MAXSHL,*), DSEND(2), ISEND(2)
      LOGICAL NODV,NOPV,NOCONT,TKTIME,RETUR, LSEND(6)
C
C--------------------------------------
C     Calculate number of Fock matrices
C--------------------------------------
C
      IF (ITYPE .EQ. 2 .OR. ITYPE .EQ. 4) THEN
         NFMAT  = 0
      ELSE IF (ITYPE .EQ. 3 .OR. ITYPE. EQ. 9 .OR. ITYPE. EQ. 7) THEN
         NFMAT  = NDMAT
      ELSE IF (ITYPE .EQ. -5) THEN
         NFMAT  = 3*NDMAT
      ELSE IF (ITYPE .EQ. 8) THEN
         NFMAT  = 3*NDMAT*NUCDEG(IATOM)
      ELSE
         WRITE (LUPRI,'(/A,I5)')
     &' SDINIT error; specified ITYPE for TWOINT not defined. ITYPE =',
     &   ITYPE
         CALL QUIT('Specified ITYPE for TWOINT not defined in SDINIT.')
      END IF
C
C----------------------------
C     Set common-block PARINT
C----------------------------
C
      JCODE  = 0
C     ... not used any more (was used for PVM)
      JATOM  = IATOM
      JLUDAS = LUDASP
      JLUINT = LUINTA
      JLUONE = LUONEL
      JLUSOL = LUSOL
      JLUSUP = LUSUPM
      JMXDIF = MAXDIF
      JMXREP = MAXREP
      JNODES = NODES
      JPRINT = IPRPAR
      JTASK  = NTASK
      JTYPE  = ITYPE
      JNOCNT = 0
      JNODV  = 0
      JNOPV  = 0
      JRETUR = 0
      JTIMIN = 0
      JTKTIM = 0
      JSOLVN = 0
      JRELCL = 0
      JLEVYL = 0
      IF (NODV)   JNODV  = 1
      IF (NOPV)   JNOPV  = 1
      IF (NOCONT) JNOCNT = 1
      IF (RETUR)  JRETUR = 1
      IF (TKTIME) JTKTIM = 1
      IF (TIMING) JTIMIN = 1
      IF (SOLVNT) JSOLVN = 1
      IF (RELCAL) JRELCL = 1
      IF (LEVYLE) JLEVYL = 1
      J2TYP  = I2TYP
      JCEDIF = ICEDIF
C
      CALL IZERO(JDOSYM,8)
      CALL IZERO(JDOREP,8)
      CALL IZERO(JTIDS,MAXNOD)
      CALL IZERO(JDOCOR,MXCOOR)
C
      CALL ICOPY(NODES,MTIDS,1,JTIDS,1)
C
      DO 100 I = 0, MAXREP
         IF (DOREPS(I))  JDOREP(I+1) = 1
         IF (DOSYM(I+1)) JDOSYM(I+1) = 1
 100  CONTINUE
C
      DO 110 J = 1, MXCENT
         DO 120 I = 1, 3
            IF (DOCOOR(I,J)) JDOCOR(I,J) = 1
 120     CONTINUE
 110  CONTINUE
C
#if defined (VAR_MPI)
      IF (NEWGEO) THEN
         CALL MPI_BCAST(NMLINE,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MLINE,80*NMLINE,MPI_CHARACTER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
C........inforb.h
         NUMELM=ICOMMSIZE(I1_INFORB,I2_INFORB)
         CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C........note that with icommsize COMMON blocks must be sent as integers...
         CALL MPI_BCAST(I1_INFORB,NUMELM,MPI_INTEGER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
C
         CALL MPI_BCAST(UNCONT,1,MPI_LOGICAL,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MAXPRI,1,MPI_INTEGER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         NEWGEO = .FALSE.
      ELSE
         NFINIS = 0
         CALL MPI_BCAST(NFINIS,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      END IF
C
C     Send Panas value and screening threshold
      DSEND(1) = PANAS
      DSEND(2) = SCREEN
      CALL MPI_BCAST(DSEND,2,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
#if defined (PRG_DIRAC)
C     Send SMLV1C and ONECAP,INTV1C
      LSEND(1) = SMLV1C
      LSEND(2) = ONECAP
      CALL MPI_BCAST(LSEND,2,MPI_LOGICAL,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(INTV1C,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
C
C     Send parint.h
      NUMELM = ICOMMSIZE(I1_PARINT,I2_PARINT)
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
C........note that with icommsize COMMON blocks must be sent as integers...
      CALL MPI_BCAST(I1_PARINT,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ISEND(1) = NDMAT
      ISEND(2) = NFMAT
      CALL MPI_BCAST(ISEND,2,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      IF(NDMAT.GT.0.AND.ITYPE.NE.4) THEN
        CALL MPI_BCAST(DMAT,N2BASX*NDMAT,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ENDIF
C
      IF(NFMAT.GT.0) THEN
        CALL MPI_BCAST(IREPDM,NDMAT,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ENDIF
      IF(ITYPE.NE.4) THEN
        IF (ITYPE.EQ.2) THEN
          NFCTYP = NDMAT
        ELSE
          NFCTYP = NFMAT
        ENDIF
        CALL MPI_BCAST(IFCTYP,NFCTYP,MPI_INTEGER,MPARID,
     &                 MPI_COMM_WORLD,IERR)
      END IF
      CALL MPI_BCAST(HFXFAC,1,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
C
C     ITYPE = 3 : direct SCF (AO-basis)
C     ITYPE = 9 : direct SCF (SO-basis)
C     ITYPE = 2 : gradient
C     ITYPE = 4 : (rs|**)-distributions
C
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.9.OR.ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
C       ... broadcast screening matrix
        N2RED = NSYMBL*NSYMBL
        CALL MPI_BCAST(N2RED,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(DMRAO,N2RED*NDMAT,MPI_DOUBLE_PRECISION,MPARID,
     &                 MPI_COMM_WORLD,IERR)
        IF (ITYPE .EQ. 2) THEN
           CALL MPI_BCAST(GABRAO,4*N2RED,MPI_DOUBLE_PRECISION,MPARID,
     &                    MPI_COMM_WORLD,IERR)
        ELSE
           CALL MPI_BCAST(GABRAO,N2RED,MPI_DOUBLE_PRECISION,MPARID,
     &                    MPI_COMM_WORLD,IERR)
        END IF
      END IF
      IF(ITYPE.EQ.9) THEN
        N2RED = MAXSHL*MAXSHL
        CALL MPI_BCAST(N2RED,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(DMRSO,N2RED*NDMAT,MPI_DOUBLE_PRECISION,MPARID,
     &                 MPI_COMM_WORLD,IERR)
      END IF
#endif
      RETURN
      END
C  /* Deck dopar */
      SUBROUTINE DOPAR(INDEX,MTOTTK,NBATCH,NSTAT)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993
C
#include "implicit.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
C
      DIMENSION INDEX(MTOTTK), NSTAT(NODES), ISHELL(2)
C
      CALL IZERO(NSTAT,NODES)
C
C------------------------------------------
C     Start loop over two outermost shells.
C------------------------------------------
C
      NTEST  = 1
      NUMTSK = 0
      IPLACE = 0
C
      DO 100, I = 1,NBATCH
C
C        Ready to receive request from any node.
C
#if defined (VAR_MPI)
         CALL MPI_RECV(NODE,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
C
         CALL MPI_SEND(NTEST,1,MPI_INTEGER,NODE,30,
     &                 MPI_COMM_WORLD,IERR)
#endif
C
         DO 200 J = 1,NTASK
            IPLACE = IPLACE + 1
            NUMTSK = NUMTSK + 1
            ISHELL(1) = INDEX(IPLACE)
            ISHELL(2) = IPLACE
C
#if defined (VAR_MPI)
            CALL MPI_SEND(ISHELL,2,MPI_INTEGER,NODE,30,
     &                    MPI_COMM_WORLD,IERR)
#endif
C
C           Pack only one null-task.
C
            IF (INDEX(IPLACE) .EQ. 0) THEN
               IPLACE = IPLACE + NTASK - J
               NUMTSK = NUMTSK - 1
               GO TO 300
            END IF
C
            IF (IPRPAR .GT. 10) THEN
               WRITE(LUPRI,'(5X,I5,A10,I5,A8,I4)')
     &              NUMTSK, 'Send task', INDEX(IPLACE), 'to node', NODE
            END IF
  200    CONTINUE
C
  300    CONTINUE
         NWHO = MTIDS(NODE)
C
C
C        Who-got-the-last-batch statistics
C
         NSTAT(NODE) = NSTAT(NODE) + 1
C
  100 CONTINUE
C
C------------------------------------------------
C     Instruct nodes to return the final results.
C------------------------------------------------
C
      NTEST = -1
      DO 400 I = 1,NODES
C
C           Ready to receive request from any node.
C
#if defined (VAR_MPI)
         CALL MPI_RECV(NODE,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
C
         CALL MPI_SEND(NTEST,1,MPI_INTEGER,NODE,30,
     &                 MPI_COMM_WORLD,IERR)
#endif
C
  400 CONTINUE
      RETURN
      END
C  /* Deck rvres */
      SUBROUTINE RVRES(WORK,FMAT,ITYPE,MTOTTK,NREQUI,
     &                 NUMGRA,NUMHES,MAXDIF,TIMES,DINTSKP,LDINTSKP)
C
C     Written by Henrik Koch and Paal Dahle July-1994
C
#include "implicit.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
C
#endif
      PARAMETER (D1 = 1.0D0)
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"
C
      DIMENSION FMAT(*),WORK(NREQUI),TIMES(MTOTTK),DINTSKP(LDINTSKP),
     &          DNSKBF(3*3*14)
C
C     Clean matrices.
C
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,NUMGRA)
         IF (MAXDIF .GT. 1)
     &      CALL DZERO(HESSEE,NUMHES)
CTROND      ELSE IF (ITYPE .EQ. 3) THEN
CTROND         CALL DZERO(FMAT,NREQUI)
      ELSE IF (ITYPE .EQ. -5) THEN
CTROND        CALL DZERO(FMAT,3*N2BASX)
         CALL DZERO(SUS2EL,9)
      ELSE IF (ITYPE .EQ. 8) THEN
         CALL DZERO(FMAT,3*N2BASX)
      END IF
C
C---------------------------
C     Start loop over nodes.
C---------------------------
C
      ICOUNT = 0
      DO 100 I = 1,NODES
C
#if defined (VAR_MPI)
         CALL MPI_RECV(ITASK,1,MPI_INTEGER,MPI_ANY_SOURCE,40,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
C
         ICOUNT = ICOUNT + ITASK
C
         NWHO = ISTAT(MPI_SOURCE)
C
         IF(ITYPE.EQ.3.OR.ITYPE.EQ.9.OR.ITYPE.EQ.2) THEN
             CALL MPI_RECV(DNSKBF,LDINTSKP,MPI_DOUBLE_PRECISION,
     &           NWHO,40,MPI_COMM_WORLD,ISTAT,IERR)
         ENDIF
C
Cjth & trond - this is a bit risky
C
         CALL MPI_RECV(WORK,NREQUI,MPI_DOUBLE_PRECISION,NWHO,40,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
#endif
C
C        Add results to final matrices.
C
         IF (ITYPE .EQ. 2) THEN
            CALL DAXPY(NUMGRA,D1,WORK,1,GRADEE,1)
            DO J = 1,LDINTSKP
              DINTSKP(J) = DINTSKP(J) + DNSKBF(J)
            ENDDO
            IF (MAXDIF .GT. 1)
     &         CALL DAXPY(NUMHES,D1,WORK(NUMGRA+1),1,HESSEE,1)
C
         ELSE IF (ITYPE .EQ. 3.OR.ITYPE.EQ.9) THEN
            CALL DAXPY(NREQUI,D1,WORK,1,FMAT,1)
            DO J = 1,8
              DINTSKP(J) = DINTSKP(J) + DNSKBF(J)
            ENDDO
C
         ELSE IF (ITYPE .EQ. -5 .OR. ITYPE .EQ. 8) THEN
            DO 200 J = 0, 2
               LPOS = J*NNBASX + 1
               MPOS = J*N2BASX + 1
               CALL DAXPY(NNBASX,D1,WORK(LPOS),1,FMAT(MPOS),1)
 200        CONTINUE
C
            IF (ITYPE .EQ. -5) THEN
               NPOS = 3*NNBASX + 1
               CALL DAXPY(9,D1,WORK(NPOS),1,SUS2EL,1)
            END IF
         END IF
C
C        The time each index-pair calculation need is collected in TIMES
C
         IF (TIMING) THEN
            DO 300, K = 1, ITASK
C
#if defined (VAR_MPI)
               CALL MPI_RECV(NUMBER,1,MPI_INTEGER,NWHO,40,
     &                       MPI_COMM_WORLD,ISTAT,IERR)
               CALL MPI_RECV(VALUE,1,MPI_DOUBLE_PRECISION,NWHO,40,
     &                       MPI_COMM_WORLD,ISTAT,IERR)
#endif
C
               TIMES(NUMBER) = VALUE
C
  300       CONTINUE
         END IF
  100 CONTINUE
C
      IF (ICOUNT .NE. MTOTTK) THEN
        WRITE(LUPRI,'(/5X,A)') 'Error in parallel calculation!'
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks to calculate  : ',
     &       MTOTTK
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks    calculated : ',
     &       ICOUNT
        CALL QUIT('Inconsistence in RVRES -> MTOTTK .NE. ICOUNT')
      ENDIF
C
      RETURN
      END
C  /* Deck parsrt */
      SUBROUTINE PARSRT(MTOTTK,NPOS,SORTED,TIMES)
C
C     Written by Paal Dahle December 1993.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "infpar.h"
C
      DIMENSION  TIMES(MTOTTK), NPOS(0:MTOTTK+1), SORTED(0:MTOTTK+1)
C
C------------------------------------------------------
C     The array TIMES containing CPU-times is sorted by
C     increasing CPUs and stored in the array SORTED.
C------------------------------------------------------
C
C     SORTED(0) must be less than any possible CPU-time,
C     and TIMMAX must be greater than any possible CPU-time
C
      ITEMS     = MTOTTK
      TIMMAX    = 99999D00
      SORTED(0) =  -1.0D00
      SORTED(1) = TIMES(1)
      NPOS  (0) = 0
      NPOS  (1) = 1
C
      DO 100,ITEM = 2,ITEMS
C
         LASTLO       = 0
         LASTHI       = ITEM
         NUMBER       = INT(ITEM/2)
         SORTED(ITEM) = TIMMAX
C
 200     CONTINUE
C
         IF (TIMES(ITEM) .EQ. SORTED(NUMBER) .OR.
     &      (TIMES(ITEM) .GT. SORTED(NUMBER) .AND.
     &       TIMES(ITEM) .LT. SORTED(NUMBER+1))) THEN
C
            DO 300, NCOUNT = ITEM,NUMBER,-1
               NPOS  (NCOUNT+1) = NPOS  (NCOUNT)
               SORTED(NCOUNT+1) = SORTED(NCOUNT)
 300        CONTINUE
C
            NPOS  (NUMBER+1) = ITEM
            SORTED(NUMBER+1) = TIMES(ITEM)
            GO TO 100
C
         ELSE IF (TIMES(ITEM) .LT. SORTED(NUMBER)) THEN
C
            LASTHI = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         ELSE IF (TIMES(ITEM) .GT. SORTED(NUMBER)) THEN
C
            LASTLO = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         END IF
 100  CONTINUE
C
      IF (IPRPAR .GE. 10) THEN
         DO 400 I = 1, ITEMS
            WRITE(LUPRI,'(15X,A10,I5,A9,F8.3)')
     &           'IJ-shell =', NPOS(I), '   time =', SORTED(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck timres */
      SUBROUTINE TIMRES(NSTAT)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993
#include "implicit.h"
      PARAMETER ( D0 = 0.00D00 )
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "infpar.h"
C
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#endif
C
      DIMENSION NSTAT(NODES),CPUS(MAXNOD),WALLS(MAXNOD),RBUF(2)
      CHARACTER SECTID*12, INTTID*12
C
      SAVE TSOFAR
      DATA TSOFAR/0.0D0/
C
#ifdef PRG_DIRAC
      IF (RELGRD) DMXCPU = D0
      IF (RELGRD) DMXWLL = D0
#endif
      DO 100 I = 1,NODES
#if defined (VAR_MPI)
         CALL MPI_RECV(RBUF,2,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     &                 50,MPI_COMM_WORLD,ISTAT,IERR)
         NODE = ISTAT(MPI_SOURCE)
         CPUS(NODE)  = RBUF(1)
         WALLS(NODE) = RBUF(2)
#endif
C
         TSOFAR = TSOFAR + CPUS(NODE)
#ifdef PRG_DIRAC
         IF ( RELGRD ) THEN
            DMXCPU = MAX(DMXCPU,CPUS(NODE))
            DMXWLL = MAX(DMXWLL,WALLS(NODE))
         END IF
#endif
C
  100 CONTINUE
C
      IF (IPRPAR .GT. 2) THEN
         CALL AROUND('Overall statistics for distribution of batches')
C
         WRITE(LUPRI,'(11X,5(A,4X))')    'Node',
     &                                  'Batches',
     &                                  'CPU time',
     &                                  'Wall time',
     &                                  'Efficiency'
         WRITE(LUPRI,'(11X,5(A,4X)/)')   '----',
     &                                  '-------',
     &                                  '--------',
     &                                  '---------',
     &                                  '----------'
C
         DO 200 NODE = 1, NODES
            IF (WALLS(NODE) .NE. 0.0D0) THEN
               EFFI = (CPUS(NODE)/WALLS(NODE))*100
            ELSE
               EFFI = 0.0D0
            END IF
            WRITE(LUPRI,'(11X,I3,4X,I6,6X,F8.2,5X,F8.2,6X,F6.2)')
     &                  NODE,NSTAT(NODE),CPUS(NODE),WALLS(NODE),EFFI
 200     CONTINUE
      END IF
C
      INTTID = SECTID(TSOFAR)
#ifdef PRG_DIRAC
      IF ( RELGRD ) THEN
         DMXWLT = DMXWLL
         DMXCPT = DMXCPU
      END IF
#endif
C
      WRITE(LUPRI,'(/A,A12)')
     &     '>>>> Total CPU  time used in Hermit nodes so far :',
     &     INTTID
C
      RETURN
      END
C  /* Deck hernod */
      SUBROUTINE HERNOD(WORK,MWORK)
C
C     Written by Henrik Koch and Paal Dahle July-1994
C
C    *****************************************************************
C    *    This is the node program for the construction of fock-     *
C    *   matrices, derivated fock-matrices and expectation values.   *
C    *****************************************************************
C
#include "implicit.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "suscpt.h"
#include "inforb.h"
#include "infpar.h"
#include "molinp.h"
#include "blocks.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbgrd.h"
#else
#include "gnrinf.h"
#include "energy.h"
#endif
C
      DIMENSION WORK(MWORK),DINTSKP(14*3*3)
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,NEWGEO,FINISH
C
      CALL QENTER('HERNOD')
      LWORK = MWORK
#include "memint.h"
C
      PARHER = .TRUE.
C
      FINISH = .FALSE.
C
      CALL DZERO(DINTSKP,3*3*14)
C
 100  CONTINUE
C
C     Start timing.
C
      CALL GETTIM(CPU1,WALL1)
C
C     This is a slave in a parallel run
C
      SLAVE = .TRUE.
      DEBUG  = .FALSE.
C
C
C     ***************************************************
C     ***** Receive initialization data from master *****
C     ***************************************************
C
      CALL RVINIT(WORK,KFREE,LFREE,NDMAT,NFMAT,ITYPE,
     &            IATOM,MAXDIF,MYTID,NODV,NOPV,NOCONT,RETUR,
     &            TKTIME,NEWGEO,FINISH,I2TYP,ICEDIF,SCREEN,
     &            KFMAT,KDMAT,KIFC,KIRD,KDMRAO,KGABAO,KDMRSO)
C
      IF (FINISH) GOTO 999
C
C     ***************************
C     ***** Allocate memory *****
C     ***************************
C
C
      IF (ITYPE .EQ. 2) THEN
         NUMGRA = MXCOOR
         IF (MAXDIF .GT. 1) THEN
            NUMHES = NUMGRA*NUMGRA
            NREQUI = NUMGRA + NUMHES
         ELSE
            NUMHES = 0
            NREQUI = NUMGRA
         END IF
         LFMAT  = 0
      ELSE IF (ITYPE .EQ.3.OR.ITYPE.EQ.9) THEN
         NREQUI = 0
         LFMAT  = NFMAT*N2BASX
      ELSE IF (ITYPE .EQ. -5) THEN
         NUMSUS = 9
         NREQUI = 3*NNBASX + NUMSUS
         LFMAT  = 0
      ELSE IF (ITYPE .EQ. 8) THEN
         NREQUI = 3*NNBASX
         LFMAT  = 0
      ELSE
         WRITE (LUPRI,*) ' ERROR in HERNOD: unknown ITYPE:',ITYPE
         CALL QUIT('Unknown ITYPE in HERNOD')
      END IF
C
C     Determine total number of tasks - MTOTTK.
      CALL DET_TASK (I2TYP,1,MTOTTK,NPOS)
C
      CALL MEMGET('INTE',KINDEX,MTOTTK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTIMES,MTOTTK,WORK,KFREE,LFREE)
C
C     ********************************
C     ***** Initialize workspace *****
C     ********************************
C
      CALL DZERO(WORK(KFMAT),LFMAT)
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,NUMGRA)
         IF (MAXDIF.GT.1) CALL DZERO(HESSEE,NUMHES)
      ELSE IF (ITYPE .EQ. -5) THEN
         CALL DZERO(SUS2EL,NUMSUS)
      END IF
C
C     ***********************************
C     ***** Calculate the integrals *****
C     ***********************************
C
      GMAT   = 0
      JPRINT = 0
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      ISHLA  = 0
C
      CALL TWOINT(WORK(KFREE),LFREE,WORK(KFMAT),WORK(KDMAT),NDMAT,
     &     WORK(KIRD),WORK(KIFC),GMAT,IDUMMY,IDUMMY,ITYPE,
     &     MAXDIF,IATOM,NODV,NOPV,NOCONT,TKTIME,JPRINT,
     &     IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,ISHLA,I2TYP,
     &     ICEDIF,SCREEN,WORK(KGABAO),
     &     WORK(KDMRAO),WORK(KDMRSO),DINTSKP,RELCAL,
     &     WORK(KINDEX),WORK(KTIMES))
C
C     ****************************************
C     ***** Send final results to master *****
C     ****************************************
C
      CALL SDRES(WORK(KFMAT),LFMAT,NUMGRA,
     &     NUMHES,NUMSUS,CPU1,WALL1,ITYPE,DINTSKP,
     &     WORK(KINDEX),WORK(KTIMES))
C
C     We want the slave program to go until there is nothing
C     more to calculate, K.Ruud, May-95
C
 20   CONTINUE
      CALL MEMREL('HERNOD',WORK,1,KDMAT,KFREE,LFREE)
      GOTO 100
C
 999  CONTINUE
C
      CALL QEXIT('HERNOD')
      END
C  /* Deck rvinit */
      SUBROUTINE RVINIT(WORK,KFREE,LFREE,NDMAT,NFMAT,ITYPE,
     &                  IATOM,MAXDIF,MYTI,NODV,NOPV,NOCONT,RETUR,
     &                  TKTIME,NEWGEO,FINISH,I2TYP,ICEDIF,SCREEN,
     &                  KFMAT,KDMAT,KIFC,KIRD,KDMRAO,KGABAO,KDMRSO)
C
C     Written by Henrik Koch and Paal Dahle
C     Revised Oct 4 1996 - tsaue
C
#include "implicit.h"
      PARAMETER(D0 = 0.0D0)
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "hrunit.h"
#include "dorps.h"
#include "abainf.h"
#include "symmet.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
#endif
#if defined (MOD_DFT)
#include "dftcom.h"
#include "dfterg.h"
#endif
#include "molinp.h"
#include "nuclei.h"
#include "infpar.h"
#include "parint.h"
C
#include "inforb.h"
#include "cbisol.h"
#include "cbirea.h"
C
      DIMENSION WORK(*), DSEND(2), ISEND(2)
      CHARACTER*6 NAVN
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,NEWGEO,FINISH, LSEND(6)
C
C-----------------------------------------
C     Receive initializations from master.
C-----------------------------------------
C
#if defined (VAR_MPI)
      CALL MPI_BCAST(NMLINE,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      IF (NMLINE .LT. 0) THEN
         FINISH = .TRUE.
         RETURN
      ELSE IF (NMLINE .GT. 0) THEN
         CALL MPI_BCAST(MLINE,80*NMLINE,MPI_CHARACTER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
C........inforb.h
         CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(I1_INFORB,NUMELM,MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
C........
         CALL MPI_BCAST(UNCONT,1,MPI_LOGICAL,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MAXPRI,1,MPI_INTEGER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         NEWGEO = .TRUE.
      ELSE
         NEWGEO = .FALSE.
      END IF
C
C     Receive PANAS and screening threshold
      CALL MPI_BCAST(DSEND,2,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      PANAS  = DSEND(1)
      SCREEN = DSEND(2)
#if defined (PRG_DIRAC)
C     Receive SMLV1C and ONECAP,INTV1C
      CALL MPI_BCAST(LSEND,2,MPI_LOGICAL,MPARID,
     &               MPI_COMM_WORLD,IERR)
      SMLV1C = LSEND(1)
      ONECAP = LSEND(2)
      CALL MPI_BCAST(INTV1C,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
C.....receive parint.h
      CALL MPI_BCAST(NUMELM,1,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(I1_PARINT,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
C
C------------------------------------
C     Retain information from PARINT.
C------------------------------------
C
      IATOM  = JATOM
      ITYPE  = JTYPE
      IPRPAR = JPRINT
      LUDASP = JLUDAS
      LUINTA = JLUINT
      LUONEL = JLUONE
      LUSOL  = JLUSOL
      LUSUPM = JLUSUP
      MAXDIF = JMXDIF
      MAXREP = JMXREP
      NODES  = JNODES
      NTASK  = JTASK
      NOCONT = JNOCNT .EQ. 1
      NODV   = JNODV  .EQ. 1
      NOPV   = JNOPV  .EQ. 1
      RETUR  = JRETUR .EQ. 1
      TIMING = JTIMIN .EQ. 1
      TKTIME = JTKTIM .EQ. 1
      SOLVNT = JSOLVN .EQ. 1
      I2TYP  = J2TYP
      RELCAL = JRELCL .EQ. 1
      LEVYLE = JLEVYL .EQ. 1
      ICEDIF = JCEDIF
C
      CALL IZERO(MTIDS,MAXNOD)
      CALL ICOPY(NODES,JTIDS,1,MTIDS,1)
      DO 100 I = 1, NODES
         IF (MTIDS(I) .EQ. MYTI) NODE = I
 100  CONTINUE
C
      DO 200 I = 0, MAXREP
         DOREPS(I)  = JDOREP(I+1) .EQ. 1
         DOSYM(I+1) = JDOSYM(I+1) .EQ. 1
 200  CONTINUE
      DO 210 J = 1, MXCENT
         DO 220 I = 1, 3
            DOCOOR(I,J) = JDOCOR(I,J) .EQ. 1
 220     CONTINUE
 210  CONTINUE
C
C     Set hermit.
C
      IF (NEWGEO .AND. SOLVNT) THEN
         NEWGEO = .TRUE.
         NUCIND = NUCIND + 1
         NUCDEP = NUCDEP + 1
         NATOMS = NATOMS + 1
         NCNTCV = NUCIND
         NAVN(1:4) = 'cav '
         NCLINE(NUCIND) = 0
         NAMEX(3*NUCIND)     = NAVN(1:4)//' z'
         NAMEX(3*NUCIND - 1) = NAVN(1:4)//' y'
         NAMEX(3*NUCIND - 2) = NAVN(1:4)//' x'
         NAMDEP(NUCDEP)     = NAMEX(3*NUCIND)(1:4)//'  '
         NAMDPX(3*NUCDEP-2) = NAMEX(3*NUCIND-2)
         NAMDPX(3*NUCDEP-1) = NAMEX(3*NUCIND-1)
         NAMDPX(3*NUCDEP  ) = NAMEX(3*NUCIND  )
         IF (NUCDEP .GT. MXCENT) THEN
            WRITE (LUPRI,'(//2A,/A,I5)')
     &         ' Too many atomic centers: MXCENT exceed in READIN for',
     &         ' solvent cavity,',' Current limit:',MXCENT
            CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
         END IF
         READ (NAVN,'(A4)') NAMN(NUCIND)
         CORD(1,NUCIND) = D0
         CORD(2,NUCIND) = D0
         CORD(3,NUCIND) = D0
         ISTBNU(NUCIND) = 7
         CHARGE(NUCIND) = D0
         CALL NUCPRO(WORK(KFREE),LFREE)
         NEWGEO = .FALSE.
      END IF
C
      IF (NEWGEO) RDINPC = .FALSE.
C
      CALL SETHER(0,NEWGEO,RELCAL,WORK(KFREE),LFREE)
C
C     Set common-block BLOCKS
C
      CALL PAOVEC(WORK(KFREE),LFREE,0,0)
C
#if defined (VAR_MPI)
      CALL MPI_BCAST(ISEND,2,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
      NDMAT = ISEND(1)
      NFMAT = ISEND(2)
      IF(NDMAT.GT.0.AND.ITYPE.NE.4) THEN
        LDMAT = NDMAT*N2BASX
        CALL MEMGET('REAL',KDMAT,LDMAT,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KDMAT),LDMAT,MPI_DOUBLE_PRECISION,
     &               MPARID,MPI_COMM_WORLD,IERR)
      ELSE
        KDMAT = KFREE
      ENDIF
      LFMAT = NFMAT*N2BASX
      IF(NFMAT.GT.0) THEN
        CALL MEMGET('REAL',KFMAT,LFMAT,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIRD,NDMAT,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KIRD),NDMAT,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      ELSE
        KFMAT = KFREE
        KIRD  = KFREE
      ENDIF
      IF(ITYPE.EQ.4) THEN
        KIFC = KFREE
      ELSEIF (ITYPE.EQ.2) THEN
        CALL MEMGET('INTE',KIFC,NDMAT,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KIFC),NDMAT,MPI_INTEGER,MPARID,
     &                 MPI_COMM_WORLD,IERR)
      ELSE
        CALL MEMGET('INTE',KIFC,NFMAT,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KIFC),NFMAT,MPI_INTEGER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
      END IF
      CALL MPI_BCAST(HFXFAC,1,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.9.OR.ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
        CALL MPI_BCAST(N2RED,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
        LDMRAO = NDMAT*N2RED
        CALL MEMGET('REAL',KDMRAO,LDMRAO,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KDMRAO),LDMRAO,MPI_DOUBLE_PRECISION,
     &               MPARID,MPI_COMM_WORLD,IERR)
        IF (ITYPE.EQ.2) THEN
           CALL MEMGET('REAL',KGABAO,4*N2RED,WORK,KFREE,LFREE)
           CALL MPI_BCAST(WORK(KGABAO),4*N2RED,MPI_DOUBLE_PRECISION,
     &                  MPARID,MPI_COMM_WORLD,IERR)
        ELSE
           CALL MEMGET('REAL',KGABAO,N2RED,WORK,KFREE,LFREE)
           CALL MPI_BCAST(WORK(KGABAO),N2RED,MPI_DOUBLE_PRECISION,
     &                  MPARID,MPI_COMM_WORLD,IERR)
        END IF
      ELSE
        KDMRAO = KFREE
        KGABAO = KFREE
      ENDIF
      IF(ITYPE.EQ.9) THEN
        CALL MPI_BCAST(N2RED,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
        LDMRSO = NDMAT*N2RED
        CALL MEMGET('REAL',KDMRSO,LDMRSO,WORK,KFREE,LFREE)
        CALL MPI_BCAST(WORK(KDMRSO),LDMRSO,MPI_DOUBLE_PRECISION,
     &               MPARID,MPI_COMM_WORLD,IERR)
      ELSE
         KDMRSO = KFREE
      ENDIF
#endif
C
      RETURN
      END
C  /* Deck parlop */
      SUBROUTINE PARLOP(WORK,LWORK,FMAT,DMAT,NDMAT,GMAT,MAXDER,
     &     EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,SOFOCK,
     &     DISTRI,LONDON,
     &     SPNORB,PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,
     &     THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &     SQ12EL,INDHER,INDHSQ,IODDHR,
     &     IREPDM,IFCTYP,ADISTR,I2TYP,ICEDIF,SCREEN,
     &     GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,
     &     INDEX,TIMES)
C
C     Copied from TWOLOP and rewritten by Paal Dahle Nov.1994
C     Minor revisions Oct 1996 T.Saue
C
#include "implicit.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "dummy.h"
      PARAMETER (LUPAS = 27, LUPAO = 29)
      PARAMETER(D0=0.0D0)
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, DISTRI, NOCONT, SPNORB,
     &        RETUR, FIRST, SQ12EL, LONDON, SUSCEP, ADISTR, RELCAL,
     &        SOFOCK
      DIMENSION DMAT(*), FMAT(*), GMAT(*), INDHSQ(*), IODDHR(*),
     &          INDHER(*), WORK(LWORK), IREPDM(*), IFCTYP(*),
     &          GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*),
     &          IJSHEL(2),INDEX(*),TIMES(*)
C
#include "cbisol.h"
#include "twocom.h"
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#include "dcbham.h"
#else
#include "energy.h"
#endif
#include "taymol.h"
#include "taysol.h"
#include "suscpt.h"
#include "blocks.h"
#include "symmet.h"
#include "infpar.h"
C
#include "inforb.h"
C
      IF (JPRINT .GT. 5) CALL TITLER('Output from PARLOP','*',103)
C
      FIRST = .TRUE.
      DIRAC = RELCAL
      IF (EXPECT .AND. .NOT.NOPV) THEN
         REWIND LUPAO
      END IF
      IF (SUSCEP) THEN
         IF (.NOT.NOPV) REWIND LUPAO
         IF (.NOT.NOPV) REWIND LUPAS
         CALL DZERO(SUS2EL,9)
      END IF
C
      DOSCRN = .FALSE.
      IF(DIRFCK.OR.SOFOCK.OR.(EXPECT.AND.MAXDER.EQ.1))  THEN
         IF (DIRFCK.OR.SOFOCK) THEN
           CALL DZERO(DINTSKP,8)
         ELSE
           CALL DZERO(DINTSKP,3*14*3)
         END IF
        IF(SCREEN.GT.D0) THEN
          DOSCRN = .TRUE.
          SCRTHR = SCREEN
          ICEFLG = ICEDIF
          NCM = 0
          NEM = 0
          DO I = 1,NDMAT
            IX  = IFCTYP(I) / 10
            IY  = MOD(IFCTYP(I),10)
            IC  = MOD(IY,2)
            IE  = (IY - IC)/2
            IF (IX.EQ.2) IC = 0
C           ... no Coulomb term for antisymmetric density matrix
            NCM = NCM + IC
            NEM = NEM + IE
          ENDDO
        ENDIF
      ENDIF
      IF(I2TYP.EQ.0) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.3) THEN
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ENDIF
C
C     ***********************************
C     ***** Send request for a task *****
C     ***********************************
C
      ICOUNT = 0
C
  100 CONTINUE
C
#if defined (VAR_MPI)
      CALL MPI_SEND(MYTID,1,MPI_INTEGER,MPARID,20,MPI_COMM_WORLD,IERR)
C
      CALL MPI_RECV(NTEST,1,MPI_INTEGER,MPARID,30,
     &               MPI_COMM_WORLD,ISTAT,IERR)
C
      IF (NTEST .EQ. -1) GOTO 300
C
      DO 200 I = 1, NTASK
C
         CALL MPI_RECV(IJSHEL,2,MPI_INTEGER,MPARID,30,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
#endif
C
C        **************************************
C        ***** Check if it is a null task *****
C        **************************************
C
         IF (IJSHEL(1).EQ. 0) THEN
            ICOUNT = ICOUNT + I - 1
            GO TO 100
         END IF
C
C        The compressed index IJSHEL is split to indices ISHELA and ISHELB.
C
         CALL UNPKIJ(IJSHEL(1),ISHELA,ISHELB)
C
         IF (TIMING) CPUTK1 = SECOND()
C
C        *****************************
C        ***** First Shell Index *****
C        *****************************
C
         ICA    = LCLASH(ISHELA)
         NHKTA  = NHKTSH(ISHELA)
         KHKTA  = KHKTSH(ISHELA)
         KCKTA  = KCKTSH(ISHELA)
         SPHRA  = SPHRSH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         MULA   = ISTBSH(ISHELA)
         MULTA  = MULT(MULA)
         NSTRA  = IORBSB(IORBSH(ISHELA,1))
         NUCA   = NUCOSH(ISHELA)
         NORBA  = NORBSH(ISHELA)
         IF (.NOT.BIGVEC) THEN
            CORAX0 = CENTSH(ISHELA,1)
            CORAY0 = CENTSH(ISHELA,2)
            CORAZ0 = CENTSH(ISHELA,3)
         END IF
         PRINTA = .TRUE.
         IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         ICB    = LCLASH(ISHELB)
         NHKTB  = NHKTSH(ISHELB)
         KHKTB  = KHKTSH(ISHELB)
         KCKTB  = KCKTSH(ISHELB)
         SPHRB  = SPHRSH(ISHELB)
         NCENTB = NCNTSH(ISHELB)
         MULB   = ISTBSH(ISHELB)
         MULTB  = MULT(MULB)
         NSTRB  = IORBSB(IORBSH(ISHELB,1))
         NUCB   = NUCOSH(ISHELB)
         NORBB  = NORBSH(ISHELB)
         IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
         END IF
         GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
         IGENAB = 1
         IF (.NOT.GENAB) IGENAB = 2
         NSETA  = NSETSH(ISHELA,IGENAB)
         NSETB  = NSETSH(ISHELB,IGENAB)
         PRINTB = PRINTA
         IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
         IF (PRINTB) THEN
            IPRINT = JPRINT
         ELSE
            IPRINT = 0
         END IF
C
C        *****************************
C        ***** Third Shell Index *****
C        *****************************
C
         ICMAX = ISHELA
         IF (SPNORB.OR.I2TYP.EQ.2) ICMAX = MAXSHL
         DO 400 ISHELC = ICSTRT,ICMAX
C
            ICC    = LCLASH(ISHELC)
            NHKTC  = NHKTSH(ISHELC)
            KHKTC  = KHKTSH(ISHELC)
            KCKTC  = KCKTSH(ISHELC)
            SPHRC  = SPHRSH(ISHELC)
            NCENTC = NCNTSH(ISHELC)
            MULC   = ISTBSH(ISHELC)
            MULTC  = MULT(MULC)
            NSTRC  = IORBSB(IORBSH(ISHELC,1))
            NUCC   = NUCOSH(ISHELC)
            NORBC  = NORBSH(ISHELC)
            IF (.NOT.BIGVEC) THEN
               CORCX0 = CENTSH(ISHELC,1)
               CORCY0 = CENTSH(ISHELC,2)
               CORCZ0 = CENTSH(ISHELC,3)
            END IF
            PRINTC = PRINTB
            IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C           ******************************
C           ***** Fourth Shell Index *****
C           ******************************
C
            IDMAX = ISHELC
            IF (.NOT.SPNORB.AND.(ISHELA.EQ.ISHELC)) IDMAX = ISHELB
            DO 500 ISHELD = IDSTRT,IDMAX
C
               ICD    = LCLASH(ISHELD)
               NHKTD  = NHKTSH(ISHELD)
               KHKTD  = KHKTSH(ISHELD)
               KCKTD  = KCKTSH(ISHELD)
               SPHRD  = SPHRSH(ISHELD)
               NCENTD = NCNTSH(ISHELD)
               MULD   = ISTBSH(ISHELD)
               MULTD  = MULT(MULD)
               NSTRD  = IORBSB(IORBSH(ISHELD,1))
               NUCD   = NUCOSH(ISHELD)
               NORBD  = NORBSH(ISHELD)
               IF (.NOT.BIGVEC) THEN
                  CORDX0 = CENTSH(ISHELD,1)
                  CORDY0 = CENTSH(ISHELD,2)
                  CORDZ0 = CENTSH(ISHELD,3)
               END IF
               GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
               IGENCD = 1
               IF (.NOT.GENCD) IGENCD = 2
               NSETC = NSETSH(ISHELC,IGENCD)
               NSETD = NSETSH(ISHELD,IGENCD)
               PRINTD = PRINTC
               IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &              PRINTD = .FALSE.
               IF (PRINTD) THEN
                  IPRINT = JPRINT
               ELSE
                  IPRINT = 0
               END IF
C
               SHAEQB = ISHELA .EQ. ISHELB
               SHCEQD = ISHELC .EQ. ISHELD
               SHABAB = (ISHELA.EQ.ISHELC) .AND. (ISHELB.EQ.ISHELD)
C
               IF (SMLV1C) THEN
                  IF (NCENTC.NE.NCENTD.AND.ICC.EQ.2.AND.ICD.EQ.2)
     &               GO TO 500
               ENDIF
C
C              *******************************
C              ***** Calculate integrals *****
C              *******************************
C
               CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,WORK,LWORK,UNDIFF,
     &              PERTUR,LONDON,SPNORB,EXPECT,SUSCEP,DDFOCK,
     &              DIRFCK,SOFOCK,DISTRI,IATOM,MULE,MULTE,MAXDER,
     &              NOCONT,NODV,NOPV,THRESH,IPRINT,FIRST,
     &              SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,
     &              ADISTR,DUMMY,ICEDIF,DINTSKP,
     &              GABRAO,DMRAO,DMRSO,IREPDM,IDUMMY,IDUMMY)
C
               IF (RETUR) THEN
                  IF (ISHELA .EQ. IPRNTA .AND.
     &                 ISHELB .EQ. IPRNTB .AND.
     &                 ISHELC .EQ. IPRNTC .AND.
     &                 ISHELD .EQ. IPRNTD) RETURN
               END IF
 500        CONTINUE
 400     CONTINUE
C
         IF (TIMING) THEN
            CPUTK2 = SECOND()
            CPUTSK = CPUTK2 - CPUTK1
            IPLACE = ICOUNT + I
C
            INDEX(IPLACE) = IJSHEL(2)
            TIMES(IPLACE) = CPUTSK
         END IF
 200  CONTINUE
C
      ICOUNT = ICOUNT + NTASK
      GOTO 100
C
 300  CONTINUE
c                 print *,'dmat-start'
c                 print *, (dmat(ii),ii=1,10)
c                 print *,'dmat-slutt'
c                 print *, (dmat(ii),ii=n2basx-10,n2basx)ccc
c
c                 print *,'fmat-start'
c                 print *, (fmat(ii),ii=1,10)
c                 print *,'fmat-slutt'
c                 print *, (fmat(ii),ii=n2basx-10,n2basx)c
c
c                 call flshfo(6)
C
C     Symmetrize skeleton Fock matrices for itypes -5 and 8.
C     This reduces the size of data transferred to Master.
C     For itypes 2 and 3 this routine is called outside TWOINT,
C     to avoid doing the same work more than once.
C
      IF (DDFOCK) THEN
         CALL SKLFCK(FMAT,WORK,LWORK,JPRINT,DIRFCK,DDFOCK,EXPECT,PERTUR,
     &               NODV,MAXDER,LONDON,NDMAT,IREPDM,IFCTYP,IATOM)
      END IF
C
C     <<<<< Print Section - Gradient and Hessian Elements >>>>>
C
      IF (EXPECT) THEN
         IF (JPRINT .GT. 0) THEN
            CALL HEADER('Two-electron integral gradient',-1)
            CALL PRIGRD(GRADEE)
            CALL HEADER('Potential energy (NN + NE + EE) gradient',-1)
            CALL ZERGRD
            CALL ADDGRD(GRADNN)
            CALL ADDGRD(GRADNA)
            CALL ADDGRD(GRADEE)
            CALL PRIGRD(GRDMOL)
            CALL HEADER('Molecular gradient',-1)
            CALL ADDGRD(GRADFS)
            CALL ADDGRD(GRADKE)
            IF (SOLVNT) THEN
               CALL ADDGRD(GSOLTT)
               CALL ADDGRD(GSOLNN)
            END IF
            CALL PRIGRD(GRDMOL)
            NCDEP3 = 3*NUCDEP
            GRDNRM = DDOT(NCDEP3,GRDMOL,1,GRDMOL,1)
            GRDNRM = SQRT(GRDNRM)
            WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *         'Molecular gradient norm:', GRDNRM
            CALL ZERGRD
            IF (MAXDER.EQ.2) THEN
               CALL HEADER('Two-electron integral Hessian',-1)
               CALL PRIHES(HESSEE,'CENTERS')
               CALL HEADER('Potential energy (NN + NE + EE) Hessian',-1)
               CALL ZERHES
               CALL ADDHES(HESSNN)
               CALL ADDHES(HESSNA)
               CALL ADDHES(HESSEE)
               IF (SOLVNT) CALL ADDHES(HSOLT2)
               IF (SOLVNT) CALL ADDHES(HSOLNN)
               CALL PRIHES(HESMOL,'CENTERS')
               CALL ZERHES
            END IF
         END IF
      END IF
      IF (LONDON .AND. MAXDER.EQ.2) THEN
         SUS2EL(2,1) = SUS2EL(1,2)
         SUS2EL(3,1) = SUS2EL(1,3)
         SUS2EL(3,2) = SUS2EL(2,3)
         IF (JPRINT .GT. 1) THEN
            CALL HEADER('Two-electron integral susceptibilities',-1)
            CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
C
C     Print Fock matrices
C
      IF (DIRFCK .AND. JPRINT.GT.3) THEN
         CALL HEADER('Fock matrix in PARLOP',-1)
         DO 600 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  600    CONTINUE
      END IF
      RETURN
      END
C  /* Deck sdres */
      SUBROUTINE SDRES(FMAT,LFMAT,NUMGRA,
     &                 NUMHES,NUMSUS,CPU1,WALL1,ITYPE,
     &                 DINTSKP,INDEX,TIMES)
C
C     Written by Henrik Koch and Paal Dahle Nov-1993.
C
#include "implicit.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "mxcent.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#if defined (MOD_DFT)
#include "dftcom.h"
#include "dfterg.h"
#include "dftbuf.h"
#endif
#include "suscpt.h"
#include "infpar.h"
#include "inforb.h"
C
      DIMENSION FMAT(LFMAT), DINTSKP(*),INDEX(*),TIMES(*)
#if defined (VAR_MPI)
      DIMENSION RBUF(2)
#endif
C
C-------------------------------
C     Put all results in RESMAT.
C-------------------------------
C
      IF (ITYPE .EQ. -5 .OR. ITYPE .EQ. 8) THEN
         DO 100 J = 0, 2
            LPOS = J*N2BASX + 1
            MPOS = J*NNBASX + 1
CTROND            CALL DCOPY(NNBASX,FMAT(LPOS),1,RESMAT(MPOS),1)
 100     CONTINUE
         IF (ITYPE .EQ. -5) THEN
            NPOS = 3*NNBASX + 1
CTROND            CALL DCOPY(NUMSUS,SUS2EL,1,RESMAT(NPOS),1)
         END IF
      END IF
C
#if defined (VAR_MPI)
      CALL MPI_SEND(ICOUNT,1,MPI_INTEGER,MPARID,40,MPI_COMM_WORLD,IERR)
#if defined (MOD_DFT)
      IF(ITYPE.EQ.7) THEN
         CALL MPI_SEND(ELCTRN,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGDRC,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGVWN,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGBCK,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(ERGLYP,1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(FMAT,LFMAT,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         RETURN
      ENDIF
#endif
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.9) THEN
         CALL MPI_SEND(DINTSKP,8,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(FMAT,LFMAT,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
      ELSE IF (ITYPE.EQ.2) THEN
         CALL MPI_SEND(DINTSKP,3*3*14,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(GRADEE,NUMGRA,MPI_DOUBLE_PRECISION,MPARID,40,
     &                 MPI_COMM_WORLD,IERR)
      ENDIF
C
C
      IF (TIMING) THEN
         DO 200, I = 1, ICOUNT
C
            CALL MPI_SEND(INDEX(I),1,MPI_INTEGER,MPARID,40,
     &                    MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(TIMES(I),1,MPI_DOUBLE_PRECISION,MPARID,40,
     &                    MPI_COMM_WORLD,IERR)
C
  200    CONTINUE
      END IF
#endif
C
C----------------------------------------------
C     Send general timing results if requested.
C----------------------------------------------
C
      IF (IPRPAR .GT. 0) THEN
         CALL GETTIM(CPU2,WALL2)
C
         CPU  = CPU2 - CPU1
         WALL = WALL2 - WALL1
C
#if defined (VAR_MPI)
         RBUF(1) = CPU
         RBUF(2) = WALL
         CALL MPI_SEND(RBUF,2,MPI_DOUBLE_PRECISION,MPARID,50,
     &                 MPI_COMM_WORLD,IERR)
#endif
C
      END IF
C
      RETURN
      END
C  /* Deck prlinp */
      SUBROUTINE PRLINP(WORD)
#include "implicit.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (NTABLE = 4)
#include "infpar.h"
      LOGICAL FIRST, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      SAVE FIRST
      DATA TABLE /'.NODES ', '.NTASK ', '.PRINT ', 'xXXXXXX'/
      DATA FIRST /.TRUE./
C
C-------------------------
C     Initialize /PRLINP/.
C-------------------------
C
      IF (.NOT. FIRST) THEN
         IF ( WORD .EQ. '*PARALL' ) THEN
 969        READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      FIRST = .FALSE.
C
C     For infpar.h
C
      SLAVE  = .FALSE.
C
      PARHER = .FALSE.
#if defined (VAR_MPI)
      NODES = NUMNOD
#else
      NODES  = 0
#endif
      NTASK  = 1
      IPRPAR = 0
      TOTWAL = 0.0D0
C
      NEWDEF = (WORD .EQ. '*PARALL')
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
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in PRLINP')
    1          CONTINUE
                  READ(LUCMD,*) NODES
               GO TO 100
    2          CONTINUE
                  READ(LUCMD,*) NTASK
               GO TO 100
    3          CONTINUE
                  READ(LUCMD,*) IPRPAR
               GO TO 100
    4          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in PRLINP')
            END IF
      END IF
  300 CONTINUE
C
#if defined (VAR_MPI)
      IF (NODES .NE. NUMNOD) THEN
         WRITE (LUPRI,'(/A,I3,A,I3)') 'Number of nodes in input ',
     &        NODES,' .ne. nodes in MPI calculation ',NUMNOD
         CALL QUIT('Wrong number of nodes in input')
      END IF
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Parsr2 */
      SUBROUTINE PARSR2(NDIM,NPOS,TIMES,INDEX,IPRINT)
C*****************************************************************************
C
C     Sort IJ-indices according to increasing run time...
C     Written by T.Saue Oct 3 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(THRES = 1.0D-2)
      DIMENSION NPOS(NDIM),TIMES(NDIM),INDEX(NDIM)
      CALL FNDAMX(TIMES,NDIM,NDIM,NPOS,NDIM,0,THRES)
CTROND FNDAMX sorts by absolute values, but unless we hope
C      for negative run times, this will be no problem...
      IF (IPRINT .GE. 10) THEN
        DO I = 1, NDIM
          WRITE(LUPRI,'(15X,A10,I5,A9,F8.3)')
     &           'IJ-shell =', INDEX(NPOS(I)),
     &           '   time =', TIMES(NPOS(I))
        ENDDO
      ENDIF
C
C     Transfer INDEX to NPOS
C
      DO I = 1,NDIM
        NPOS(I) = INDEX(NPOS(I))
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck det_task */
      SUBROUTINE DET_TASK (I2TYP,IMODE,MTOTTK,NPOS)
C*****************************************************************************
C
C     Determine tasks that are to be distributed.
C     IMODE 1 : Determine only total number
C     IMODE 2 : Fill the array NPOS
C
C     Written by Luuk Visscher, august 2001
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbham.h"
#endif
#include "blocks.h"
C
      DIMENSION NPOS(*)
C
      IF(I2TYP.EQ.0) THEN
        IASTRT = 1
        IBSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IASTRT = 1
        IBSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IASTRT = 1
        IBSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'DET_TOTTK: Unknown I2TYP =' ,I2TYP
        CALL QUIT('Unknown I2TYP !!!')
      ENDIF
C
C     *****************************
C     ***** First Shell Index *****
C     *****************************
C
      IND = 0
      DO 100 ISHELA = IASTRT,IASMAX
C
         ICA    = LCLASH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         DO 200 ISHELB = IBSTRT,ISHELA
C
            ICB    = LCLASH(ISHELB)
            NCENTB = NCNTSH(ISHELB)
C
C           Skip two-center SS contributions if requested
C           This test is provisory, in SYMLOP we take out more integrals
C           in case of symmetry-degenerate nuclei. Testing at this place,
C           will, however, take out the work in ODCVEC and EXCOEF that
C           otherwise would cause a quartic scaling term.
C
            IF (SMLV1C) THEN
               IF (NCENTA.NE.NCENTB.AND.ICA.EQ.2.AND.ICB.EQ.2) GO TO 200
            ENDIF
C
            IJ  = (ISHELA*(ISHELA-1))/2 + ISHELB
            IF (IMODE.EQ.2) NPOS(MTOTTK-IND) = IJ
            IND = IND + 1
C
  200    CONTINUE
  100 CONTINUE
C
      IF (IMODE.EQ.1) MTOTTK = IND
C
      RETURN
      END
