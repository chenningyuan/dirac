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
#if defined (VAR_MPI2)
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MV7_PAR(C,HC,LUC,LUHC,LUCSCR,SSCRFAC,IBLOCKL,
     &                  IBLOCKD,RCCTOS,IGROUPLIST,IPROCLIST,NUMBLOCK)
*
* Outer routine for sigma vector generation
* GAS version !!!!
*
* Written in terms of RASG3/SBLOCK, May 1997
* Works only for ICISTR .ne. 1, i.e. pure memory version
* has been eliminated.
* Parallel version of MV7, called from MICDV6_PAR.
* Blocks of C-coeffiecients and S-vectors are distributed among the
* nodes.
* S. Knecht, Odense Jan07
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*. array of blocklength, cpu-distribution and connection matrix
      DIMENSION IBLOCKL(NUMBLOCK),IBLOCKD(NUMBLOCK)
      DIMENSION RCCTOS(NUMBLOCK)
*.Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
*./ORBINP/ : NACOB used
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Used : NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
      CALL QENTER('MV7PAR')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MV7   ')

C
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' MV7_PAR does not work for ICISTR = 1'
        WRITE(6,*) ' SWITCH to ICISTR = 2,3 or program'
        Call Abend2(' MV7_PAR does not work for ICISTR = 1')
      END IF
C
      IATP = 1
      IBTP = 2
C
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
C     Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
C
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
C     Arrays giving allowed type combinations
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
C     Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
C     Arrays for partitioning of sigma
      NTTS = MXNTTS
                       
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
C     Batches  of S-vector
      ITTSS_ORD = 2
      
C?      IF( IODENSEPAR .eq. 2 .and. IDISTROUTE .eq. 2 ) THEN
      IF( IODENSEPAR .eq. 2 ) THEN
C       
C       partition ci-vector w.r.t to LBLOCK AND IBLOCKD
C
        CALL PART_CIV_PAR1(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,
     &                     LBLOCK,WORK(KSIOIO),ISMOST(1,ISSM),
     &                     NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                     WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD,
     &                     IBLOCKD,NUMBLOCK)
      ELSE
C       good old normal partitioning      
        CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &                 WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &                 WORK(KSIOIO),ISMOST(1,ISSM),
     &                 NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                 WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD)
      END IF

      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
      LLUC = LUC
      LLUHC = LUHC
      LLUCSCR = LUCSCR
*
      CALL MEMMAN(KSCALLOC,NUMBLOCK,'ADDL  ',2,'SCLLOC')
      CALL MEMMAN(KSCALLOC2,NUMBLOCK,'ADDL  ',2,'SCLLC2')
C     read in scalfac only once!
      IRDSCF = 1
C
      CALL RASSG3_PAR(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                WORK(KLSI1BT),WORK(KLSIBT),LLUC,LLUHC,LLUCSCR,
     &                NUMBLOCK,IBLOCKL,IBLOCKD,
     &                SSCRFAC,WORK(KSCALLOC),WORK(KSCALLOC2),
     &                IGROUPLIST,IPROCLIST,RCCTOS)
C     Eliminate local memory
      IDUM = 0
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'MV7   ')
C
      CALL QEXIT('MV7PAR')
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE RASSG3_PAR(CB,SB,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &           LUC,LUHC,LUCSCR,NBLOCK,IBLOCKL,IBLOCKD,SSCLFAC,
     &           SCALFAC_LOCAL,SCALFAC_GROUP,IGROUPLIST,IPROCLIST,
     &           RCCTOS)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to SBLOCK
*
* Parallel version of RASSG3, called from MV7_PAR.
* S. Knecht, Odense Jan07
* =====
* Input
* =====
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
* new scratch
      INTEGER MYLENSUM,NZERO,JCOLOR,JKEY
*.Scratch
      DIMENSION SB(*),CB(*),SSCLFAC(*),IBLOCKL(NBLOCK)
      DIMENSION IBLOCKD(NBLOCK),SCALFAC_LOCAL(NBLOCK)
      DIMENSION SCALFAC_GROUP(NBLOCK),RCCTOS(NBLOCK)
      DIMENSION IGROUPLIST(NMPROC)
      CHARACTER SECTID*12,WALLTID*12,WALLTID2*12,WALLTID3*12
*
      CALL QENTER('RASSG3_PAR')
      NTEST = 00
C     NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' RASSG3 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*) ' RASSG3 : NBATS = ',NBATS
      END IF
*
*. 
      NZERO = 0
      ZERO = 0.0D0
C     fill complete local scalfac arrays with zero's 
      CALL SETVEC(SCALFAC_LOCAL,ZERO,NBLOCK)
      CALL SETVEC(SCALFAC_GROUP,ZERO,NBLOCK)
*
      IF(LUHC.GT.0) CALL REWINO(LUHC)
      IF(IRDSCF.GT.0) THEN
        CALL REWINO(LUCSCR)
        CALL REWINO(LUC)
C       find number of non-zero local C-blocks in SB(1:NBLOCK)
        IF( IODENSEPAR .eq. 2 .and. IIOMOD .ne. 0 ) THEN
          CALL DCOPY(NBLOCK,SSCLFAC,1,SB,1)
          CALL DCOPY(NBLOCK,SSCLFAC,1,SCALFAC_GROUP,1)
          CALL SETVEC(SSCLFAC,ZERO,NBLOCK)
CSK          WRITE(6,*) 'I arrived here0a! myproc,mynew_id',myproc,mynew_id
        ELSE
          CALL FIND_ACTIVE_BLOCKS_PAR(LUC,-1,SB,CB,NBLOCK,IBLOCKD)
        END IF
C
C       "mpi_allsum" local SB's to SSCLFAC which then on all
C       nodes will contain the number of non-zero C-blocks in
C       the complete CI-vector
C
        CALL REDVEC(SB,SSCLFAC,NBLOCK,2,MPI_SUM,MPI_COMM_WORLD,-1)
C 
        IF(NTEST.GE.50) THEN
          IF(MYPROC.EQ.MASTER) THEN
            CALL WRTMAT(SSCLFAC,1,NBLOCK,1,NBLOCK)
          END IF
        END IF
C        
        NTEST = 00
      END IF
      
*. Purpose for parallel scheme 2:
*. - Gather all c-blocks needed for the sigma-vector-blocks to compute
*. - Use the MPI-Broadcast scheme - use of another scheme possible?
*. --> no MPI_REDUCE is needed inside the loop over s-batches!
      IF(IODENSEPAR.EQ.2) THEN
*
        ITESTTHISWAY = 1
        blocktime = 0.0D0
        starttimer = MPI_WTIME()
        CALL REWINO(LUCSCR)
        CALL REWINO(LUC)
C
C       find all c-blocks connecting to all sigma-blocks on each cpu 
C
        CALL DCOPY(NBLOCK,RCCTOS,1,SCALFAC_LOCAL,1)
C
        IF(IIOMOD.EQ.0) THEN
          CALL COPVCD_PAR_BDRIV3(CB,IBLOCKD,SSCLFAC,SCALFAC_LOCAL,
     &                           IBLOCKL,NBLOCK,MPI_COMM_WORLD,LUC,
     &                           LUCSCR,0,-1)
        ELSE
C         try two different schemes for IIOMOD == 1 ...
          IF( ITESTTHISWAY .eq. 0 ) THEN
            CALL COPVCD_PAR_BDRIV4(CB,IBLOCKD,SSCLFAC,SCALFAC_LOCAL,
     &                             IBLOCKL,NBLOCK,MPI_COMM_WORLD,LUC,
     &                             0,-1)
          ELSE
C
C           case 1: number of cpus in new group not equal to total
C           number
C           case 2: number of cpus in new group equal to total
C           number
C
            IF( NEWCOMM_PROC .ne. NMPROC ) THEN 
C
CSK          WRITE(6,*) 'I arrived here! myproc,mynew_id',myproc,mynew_id
              CALL REDVEC(SCALFAC_LOCAL,SCALFAC_GROUP,NBLOCK,2,
     &                    MPI_SUM,MYNEW_COMM,0)
C             all local node-masters call this routine!
              IF( MYNEW_ID .eq. 0 ) THEN
                CALL COPVCD_PAR_BDRIV5(CB,IBLOCKD,SSCLFAC,SCALFAC_GROUP,
     &                                 IBLOCKL,NBLOCK,ICOMM,
     &                                 IGROUPLIST,IPROCLIST)
              END IF
CSK          WRITE(6,*) 'I arrived here 2! myproc,mynew_id',myproc,mynew_id
                CALL MPI_BCAST(SCALFAC_GROUP,NBLOCK,MPI_REAL8,0,
     &                         MYNEW_COMM,IERR)
            ELSE 
C
              DO IBLK = 1, NBLOCK
                IF( SCALFAC_LOCAL(IBLK) .ne. 0.0D0 ) THEN
                  IF( SSCLFAC(IBLK) .eq. 0.0D0 ) THEN
                    SCALFAC_LOCAL(IBLK) = 0.0D0
                  END IF
                END IF
              END DO
C            
            END IF
          END IF
        END IF 
*
        CALL REWINO(LUCSCR)
        
        blocktime = blocktime + MPI_WTIME() - starttimer
*
      ENDIF
* Loop over batches over sigma blocks
      DO 500 JBATS = 1, NBATS
        ICOMPUTE = 0
        ncblocks_used = 0
*. Obtain sigma for batch of blocks
        starttime = mpi_wtime()
        sbatchtime = 0.0D0
        writetime = 0.0D0
        MYLENSUM = 0
        IF(IODENSEPAR.EQ.1) THEN
          DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
            LEN  = IBATS(8,ISBLK)
            MYLENSUM = MYLENSUM + LEN
          END DO
        ELSE
C         check whether a cpu needs to call SBLOCK_PAR
          DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
            IF(IBLOCKD(ISBLK).EQ.MYPROC) ICOMPUTE = 1
          END DO
          IF(ICOMPUTE.EQ.0) THEN
            starttimer = MPI_WTIME()
            GOTO 60
          ENDIF
        ENDIF
        CALL SBLOCK_PAR(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,
     &       CB,SB,LUC,0,0,0,0,0,MYLENSUM,NBLOCK,SSCLFAC,IBLOCKL,
     &       IBLOCKD,I1BATS(JBATS),SCALFAC_LOCAL,LUCSCR,SCALFAC_GROUP)
*. Transfer S block to permanent storage
          starttimer = MPI_WTIME()
          DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
            IF(IBLOCKD(ISBLK).NE.MYPROC) GOTO 5
            IOFF = IBATS(6,ISBLK)
            LEN  = IBATS(8,ISBLK)
            CALL ITODS(LEN,1,-1,LUHC)
            CALL TODSC(SB(IOFF),LEN,-1,LUHC)
  5         CONTINUE
          END DO
 60       CONTINUE
          writetime = writetime + MPI_WTIME() - starttimer
          sbatchtime = sbatchtime + MPI_WTIME() - starttime
          WALLTID = SECTID(writetime)
          WALLTID2 = SECTID(sbatchtime)
          IF(MYPROC.EQ.MASTER)THEN
CSK            WRITE(6,50000) WALLTID
            WRITE(6,60000) WALLTID2,JBATS
            WRITE(6,70000) ncblocks_used 
          ELSE
CSK            WRITE(12,50000) WALLTID
            WRITE(12,60000) WALLTID2,JBATS
            WRITE(12,70000) ncblocks_used 
          ENDIF
500   CONTINUE
*
*. write an EOF-mark
        CALL ITODS(-1,1,-1,LUHC)
*
C       print the statistic for COPVCD_PAR_BDRIV3/4
        IF(IODENSEPAR.EQ.2) WALLTID3 = SECTID(blocktime)
        IF(MYPROC.EQ.MASTER)THEN
          IF(IODENSEPAR.EQ.2) WRITE(6,80000) WALLTID3
        ELSE
          IF(IODENSEPAR.EQ.2) WRITE(12,80000) WALLTID3
        ENDIF
*
      IF(NTEST.GE.100) THEN
        IF(MYPROC.EQ.MASTER+3) THEN
          WRITE(6,*) ' Final S-vector on disc for myproc',myproc
          CALL WRTVCD(SB,LUHC,1,-1)
        END IF
      END IF
      CALL QEXIT('RASSG3_PAR')
*
CSK50000 FORMAT(' >>>  WALL TIME FOR WRITING SB TO DISK           :
CSK     & ',A)
60000 FORMAT(' >>>  WALL TIME IN TOTAL FOR CURRENT S-BATCH     :
     & ',A,I4)
70000 FORMAT(' >>>  NUMBER OF C-BLOCKS USED FOR COMPUTATION    :
     & ',I4)
80000 FORMAT(' >>>  WALL TIME FOR BROADCASTING C-COEFFICIENTS  :
     & ',A)
     
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SBLOCK_PAR(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
     &                      LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                      LENSUM,NUMBLOCK,SSCLFAC,IBLOCKL,IBLOCKD,
     &                      INSBLK,SCALFAC_LOCAL,LUCSCR,SCALFAC_GROUP)
*
* Generate a set of sigma blocks,
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES .eq.1 then it as assumed that only
* batches ICBAT_INI to ICBAT_END are stored on  LUC
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
*./ORBINP/ : NACOB used
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Used : NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      dimension XIJILS(MXTSOB)
      DIMENSION SSCLFAC(*),IBLOCKL(NUMBLOCK),IBLOCKD(NUMBLOCK)
      DIMENSION SCALFAC_LOCAL(NUMBLOCK)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBLCKP')
*
CTF
* Scratch space for rsbb2a, rsbb2bn, rsbb2bvn, skickj
CTF
      LSOMESQ = MXTSOB*MXTSOB
      call memman(KSOMEH,LSOMESQ,'ADDS',2,'SOMSCH')
      call memman(KSOMESCR,LSOMESQ**2,'ADDS',2,'SOMSCR')
CTF
CTF   dimension SOMESCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
CTF   dimension SOMEH(MXTSOB*MXTSOB)
      NTEST = 0
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SBLOCK_PAR : ISSPC,ICSPC ', ISSPC,ICSPC
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*. connection matrices for supergroups
*
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      MAXA0 = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXA = MAX(MAXA,MAXA0)
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
*
      MAXB = 0
      MAXB0 = IMNMX(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MAXB = MAX(MAXB,MAXB0)
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.3 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
C?    WRITE(6,*) ' MXTSOB_L = ', MXTSOB_L
      MAXIJ = MXTSOB_L ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' ICISTR,LBLOCK ',ICISTR,LBLOCK
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
C          IAIBCM(ICISPC,IAIB)
      CALL IAIBCM(ISSPC,WORK(KSIOIO))
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
*. Arrays giving block type
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESCPH(WORK(KSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &            IPHGAS,NHLFSPGP,MNHL,IADVICE)
      IF(IPRCIX.GE.3) THEN
         WRITE(6,*) 'SBLOCK_PAR : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
         WRITE(6,*) 'SBLOCK_PAR : MXADKBLK ', MXADKBLK
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.3)  WRITE(6,*) ' LSCR2 = ', LSCR2
*
      LSCR12 = MAX(LBLOCK,2*LSCR2)
C     CALL MEMMAN(KC2,LSCR12,'ADDL  ',2,'KC2   ')
C  I assume memory was allocated for blocks, so
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*. Space for active/passive blocks
      KCJPA = KC2 + 2*LSCR2
      KSIPA = KC2 + 3*LSCR2
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.Some TTS arrays
      NOOS = NOCTPA*NOCTPB*NSMCI
      NTTS = MXNTTS
*
*. for partitioning of vector
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      CALL MEMMAN(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')

*. Space for four blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      DO I1234 = 1, 1
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO

      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
*. Place perturbation integrals over one body integrals
CINSERT_START
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
*. Prepare for perturbation calculation
*
      IF(IPERTOP.NE.0) THEN
*. Matrix specifying partiotioned spaces
        CALL MEMMAN(KLH0SPC,NOCTPA*NOCTPB,'ADDL  ',2,'H0SPC ')
        CALL H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                NGAS,MXPNGAS,WORK(KLH0SPC),NELFGP)
        IF(IH0SPC.EQ.0) THEN
*. Form of perturbation in subspace has not been defined,
*. Use current IPART
          IH0INSPC(1) = IPART
        END IF
      END IF
*
      CALL SBLOCKS_PAR(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
     &     WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
     &     NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     NAEL,IATP,NBEL,IBTP,
     &     IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &     ITSOB,MAXIJ,MAXK,MAXI,INSCR,LBLOCK,
     &     LBLOCK,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
     &     SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &     ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &     MXTSOB,WORK(KSOMESCR),WORK(KSOMEH),XIJILS,
     &     WORK(KI1),WORK(KXI1S),
     &     WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,WORK(KSVST),
     &     PSSIGN,IPRDIA,LUC,ICJKAIB,WORK(KCJRES),
     &     WORK(KSIRES),WORK(KI3),WORK(KXI3S),
     &     WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
     &     MOCAA,MOCAB,IAPR,
     &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     IRESTRICT,WORK(KCONSPA),WORK(KCONSPB),WORK(KLSCLFAC),
     &     LUCBLK,IPERTOP,IH0INSPC,WORK(KLH0SPC),
     &     ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,
     &     IUSE_PA,WORK(KCJPA),WORK(KSIPA),LENSUM,NUMBLOCK,
     &     IBLOCKL,IBLOCKD,SSCLFAC,INSBLK,SCALFAC_LOCAL,LUCSCR,
     &     SCALFAC_GROUP)
*
      IF(IDC.EQ.2) THEN
*. reform
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF
*
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SBLCKP')
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SBLOCKS_PAR(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &           ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &           XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &           MXTSOB,SOMESCR,SOMEH,XIJILS,
     &           I1,XI1S,I2,XI2S,IDOH2,MXPOBS,ISTRFL,
     &           PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &           I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &           ICONSPA,ICONSPB,SCLFAC,
     &           LUCBLK,IPERTOP,IH0INSPC,IH0SPC,
     &           ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &           I_RES_AB,IUSE_PA,CJPA,SIPA,
     &           ILENSUM,NUMBLOCK,IBLOCKL,IBLOCKD,SSCLFAC,INSBLK,
     &           SCALFAC_LOCAL,LUCSCR,SCALFAC_GROUP)
*
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , Winter of 1991
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included
* INSBLK  : start number of SBLOCK
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
* NEWCOMM_PROC: temporary number of processes
* IMPINFA(NEWCOMM_PROC,1): number of blocks to calculate fo a specific proc
* IMPINFA(NEWCOMM_PROC,2): total length of all blocks
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)


      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
#include "comjep.inc"
*. Specific input
      INTEGER ISBLOCK(8,*)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*. Zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*),SSCLFAC(*),SCALFAC_LOCAL(NUMBLOCK)
      DIMENSION SCALFAC_GROUP(NUMBLOCK)
* new scratch
      INTEGER ILENSUM,JSIGMABLOCK
      INTEGER ISUMUP,NREDGRP,JCOLOR,JKEY,INTERNALSIZE
      DIMENSION IBLOCKL(NUMBLOCK),IBLOCKD(NUMBLOCK) 
      CHARACTER SECTID*12, WALLTID2*12, WALLTID*12
      CHARACTER WALLTID3*12, WALLTID4*12,WALLTID5*12
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
*.
      CALL QENTER('SBLOC')
*
*. create a new group for mpi_reduce at the end of all loops
*
      ISUMUP = 0
      JSIGMABLOCK = 0
*
CSK      IF(MYPROC.EQ.MASTER) THEN
CSK       WRITE(6,*)'IRESTRICT is ',IRESTRICT
CSK       CALL IWRTMA(ICONSPA,NOCTPA,NOCTPA,NOCTPA,NOCTPA)
CSK      ENDIF
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================'
        WRITE(6,*) '  SBLOCKS_PAR speaking '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0)
     &      WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(6,*) ' IDC PS IPERTOP', IDC,PS,IPERTOP
          WRITE(6,*) ' IDOH2 = ',IDOH2
          WRITE(6,*) ' IRES_AB=',IRES_AB
      END IF
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Initial C vector of myproc',myproc
        IF(IODENSEPAR.EQ.1) THEN
          CALL WRTVCD(CB,LUC,1,-1)
        ELSE
          IF(IIOMOD.EQ.0) CALL WRTVCD(CB,LUCSCR,1,-1)
          CALL WRTVCD(CB,LUC,1,-1)
        ENDIF
CSK        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      END IF
C
C     Get actual scalfac array from backup-array
      IF(IODENSEPAR.EQ.1) THEN
        CALL DCOPY(NUMBLOCK,SSCLFAC,1,SCLFAC,1)
      ELSE
        IF( IIOMOD .eq. 0 ) THEN
          CALL DCOPY(NUMBLOCK,SCALFAC_LOCAL,1,SCLFAC,1)
        ELSE
          IF( ITESTTHISWAY .eq. 0 ) THEN
            CALL DCOPY(NUMBLOCK,SCALFAC_LOCAL,1,SCLFAC,1)
          ELSE
            IF( NEWCOMM_PROC .ne. NMPROC) THEN
CSK              CALL COPVEC(SCALFAC_LOCAL,SCLFAC,NUMBLOCK)
              CALL DCOPY(NUMBLOCK,SCALFAC_GROUP,1,SCLFAC,1)
            ELSE
              CALL DCOPY(NUMBLOCK,SCALFAC_LOCAL,1,SCLFAC,1)
            END IF
          END IF
        END IF
      END IF
CSK      WRITE(6,*)'This is my sclfac',myproc
CSK      CALL WRTMAT(SCLFAC,1,NUMBLOCK,1,NUMBLOCK)
C
C     ============================
C      1 : Arrays for accessing C
C     ============================
C
C     Find batches of C - strings
      ITTSS_ORD = 2
C     ... well we want the count of calls for read-in to keep as 
C         little as possible
C         works ( --> ? ) currently only for IDISTROUTE = 2
C         each node creates its own batch structure since the need for 
C         C-blocks is different for a given s-block
C?      IF ( IODENSEPAR .eq. 2 .and. IDISTROUTE.EQ.2 ) THEN
      IF ( IODENSEPAR .eq. 2 ) THEN
        CALL PART_CIV_PAR2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &                     LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &                     I1CBLOCK,ICBLOCK,0,ITTSS_ORD,
     &                     SCLFAC,NUMBLOCK,0)
      ELSE
      
        CALL PART_CIV2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &                LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &                I1CBLOCK,ICBLOCK,0,ITTSS_ORD)
      END IF
C
C     Initialize sigma blocks
      DO JSBLOCK = 1, NSBLOCK
        IATP = ISBLOCK(1,JSBLOCK)
        IBTP = ISBLOCK(2,JSBLOCK)
        IASM = ISBLOCK(3,JSBLOCK)
        IBSM = ISBLOCK(4,JSBLOCK)
        IOFF = ISBLOCK(5,JSBLOCK)
        NASTR = NSSOA(IASM,IATP)
        NBSTR = NSSOB(IBSM,IBTP)
        ZERO = 0.0D0
        IF(ISBLOCK(1,JSBLOCK).GT.0)
     &  CALL SETVEC(SB(IOFF),ZERO,NASTR*NBSTR)
      END DO
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
C?    WRITE(6,*) ' IDOH2, MXEXC',IDOH2,MXEXC
      IF(IODENSEPAR.EQ.1) THEN
        REWIND LUC
      ELSE 
        REWIND LUCSCR
      END IF
      
      IF(ICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
*
      waittime   = 0.0d0
      readtime   = 0.0d0
      sigmatime  = 0.0d0
      reducetime = 0.0d0
      IOFFSET = 0
C
C     Loop over batches over C blocks
C
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
C
C       Read C blocks into core
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
C        IF(MYPROC.EQ.MASTER) THEN
C          WRITE(6,*)'NJBLOCK (MYPROC) = ', NJBLOCK, MYPROC
C        END IF
        starttimer = mpi_wtime()
        DO JJCBLOCK = 1, NJBLOCK
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
C         calculation of a certain block is node-dependent!
          IF(IODENSEPAR.EQ.1.AND.IBLOCKD(JBLOCK).NE.MYPROC) GOTO 30
C
          INTERACT = 0
C
          IF( SCLFAC(JBLOCK) .ne. 0.0D0 ) THEN
C?          WRITE(6,*) ' Nonvanishing block ', JBLOCK, MYPROC
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
C           Loop over Sigma blocks in batch
            DO JSBLOCK = 1, NSBLOCK
CSK            WRITE(6,*)'ISBLOCK(1,JSBLOCK)',ISBLOCK(1,JSBLOCK)
CSK            WRITE(6,*)'JSBLOCK =',JSBLOCK
            IDENT = 0
            IF(ISBLOCK(1,JSBLOCK).GT.0) THEN
              JSIGMABLOCK = INSBLK + JSBLOCK - 1
              IF( IODENSEPAR .eq. 2 
     &        .and. MYPROC .ne. IBLOCKD(JSIGMABLOCK)) GOTO 25
              IATP = ISBLOCK(1,JSBLOCK)
              IBTP = ISBLOCK(2,JSBLOCK)
              IASM = ISBLOCK(3,JSBLOCK)
              IBSM = ISBLOCK(4,JSBLOCK)
C
C             Are the two blocks connected by allowed excitation?
C
              IF(MXEXC.EQ.1) THEN
                 IF((ICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (ICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &               ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
C             Are they identical ?
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
 25           CONTINUE
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
          IF(INTERACT.EQ.1) THEN
CSK             IF(MYPROC.EQ.MASTER+2) THEN
CSK               WRITE(6,*) ' TTSS for C block read in  '
CSK               CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
CSK             ENDIF
            IF(IODENSEPAR.EQ.1) THEN
              CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                    NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                    PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
            ELSE
              IF(IIOMOD.EQ.0) THEN
                CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                      NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                      PL,LUCSCR,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
CSK                 WRITE(6,*) 'CB read in JOFF',JOFF
CSK           CALL WRTMAT(CB(JOFF),1,ICBLOCK(8,JBLOCK),1,ICBLOCK(8,JBLOCK))
              ELSE
                LBL = IBLOCKL(JBLOCK)
                IF ( NEWCOMM_PROC .ne. NMPROC ) THEN
                  CALL MPI_FILE_READ_AT(MPILUSCR,IOFFSET,CB(JOFF),
     &                                  LBL,MPI_REAL8,ISTAT,IERR)
CSKT                  WRITE(6,*) 'CB read in IOFFSET',IOFFSET
CSKT                  WRITE(6,*) 'CB read in JOFF,IOFFSET',JOFF,IOFFSET
CSKT                  CALL WRTMAT(CB(JOFF),1,LBL,1,LBL)
                ELSE
                  CALL MPI_FILE_READ_AT(MPILUSCR2,IOFFSET,CB(JOFF),
     &                                  LBL,MPI_REAL8,ISTAT,IERR)
CSKT                  WRITE(6,*) 'CB read in IOFFSET',IOFFSET
CSKT                  WRITE(6,*) 'CB read in JOFF,IOFFSET',JOFF,IOFFSET
CSK                  CALL WRTMAT(CB(JOFF),1,LBL,1,LBL)
                END IF
                IOFFSET = IOFFSET + LBL
CSK                WRITE(6,*)'IOFFSET',IOFFSET,MYPROC
              ENDIF
            ENDIF
              
*. Note in GSTTBL : ICOOSC only used for CI vectors in core
          ELSE
*. not relevant
C             WRITE(6,*) ' TTSS for C block skipped  '
C             CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            IF(IODENSEPAR.EQ.1) THEN
              CALL IFRMDS(LBL,-1,1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
            ELSE
              IF(IIOMOD.EQ.0) THEN
                IF( SCALFAC_LOCAL(JBLOCK) .ne. 0.0D0 ) THEN
                  CALL IFRMDS(LBL,-1,1,LUCSCR)
                  CALL SKPRCD2(LBL,-1,LUCSCR)
                END IF
              ELSE
                  IF( NEWCOMM_PROC .ne. NMPROC ) THEN
CSK                    IF( SCALFAC_LOCAL(JBLOCK) .ne. 0.0D0 ) THEN
                      IOFFSET = IOFFSET + IBLOCKL(JBLOCK)
CSK                    END IF
                  ELSE
                    IOFFSET = IOFFSET + IBLOCKL(JBLOCK)
                  END IF
CSK                WRITE(6,*)'IOFFSET2',IOFFSET,MYPROC
CSK                END IF
              END IF
            END IF 
C           mark c-block as unused
            SCLFAC(JBLOCK) = 0.0D0
          END IF
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            END IF
          END IF
*
 30       CONTINUE
*
        END DO
        readtime = readtime + mpi_wtime() - starttimer
CSK        write(6,*)'problem here3?',myproc
*
        starttimer = mpi_wtime()
        DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK
*
*. calculation of a certain block is node-dependent!
          IF( SCLFAC(ICBLK) .ne. 0.0D0 ) THEN
          IF(IODENSEPAR.EQ.1.AND.IBLOCKD(ICBLK).NE.MYPROC) GOTO 8999

            JATP = ICBLOCK(1,ICBLK)
            JBTP = ICBLOCK(2,ICBLK)
            JASM = ICBLOCK(3,ICBLK)
            JBSM = ICBLOCK(4,ICBLK)
            ICOFF = ICBLOCK(5,ICBLK)
            NJA = NSSOA(JASM,JATP)
            NJB = NSSOB(JBSM,JBTP)
*
*. Other symmetry blocks that can be obtained from this block
              CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                    LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
              DO 8765 IPERM = NPERM,1, -1
                LLASM = LASM(IPERM)
                LLBSM = LBSM(IPERM)
                LLATP = LATP(IPERM)
                LLBTP = LBTP(IPERM)
                NLLA = NSSOA(LLASM,LLATP)
                NLLB = NSSOB(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
                IF(IPERM.EQ.1) THEN
                  IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
*. Diagonal blocks, Transposing corresponds to scaling
                    IF(PS.EQ.-1.0D0) THEN
                      CALL DSCAL(NJA*NJB,PS,CB(ICOFF),1)
                    END IF
                  ELSE
*. offdiagonal blocks, explicit transposing
                    CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                    CALL DCOPY(NJA*NJB,C2,1,CB(ICOFF),1)
                    if (LSGN(IPERM).eq.-1) then
                      write(6,*) 'LROW,LCOL undef. in sblocks'
                      write(6,*) 'and LSGN(IPERM).eq.-1'
                      Call Abend2('Quitting.')
CTF                   CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                    end if
                  END IF
                END IF
*
                ncblocks_used = ncblocks_used + 1
                DO 10000 ISBLK = 1, NSBLOCK
                IF( ISBLOCK(1,ISBLK) .gt. 0 ) THEN
                  JSIGMABLOCK = INSBLK + ISBLK - 1
                  IF( IODENSEPAR .eq. 2 .and. 
     &               MYPROC .ne. IBLOCKD(JSIGMABLOCK)) GOTO 10000
C
                  IATP = ISBLOCK(1,ISBLK)
                  IBTP = ISBLOCK(2,ISBLK)
                  IASM = ISBLOCK(3,ISBLK)
                  IBSM = ISBLOCK(4,ISBLK)
                  ISOFF = ISBLOCK(5,ISBLK)
                  NIA = NSSOA(IASM,IATP)
                  NIB = NSSOB(IBSM,IBTP)
*
                  IF(NIA*NIB.EQ.0) GOTO 10000
                  IF(IRESTRICT.EQ.1.AND.
     &               (JASM.GT.IASM.OR.
     &               JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &               JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &            GOTO 10000
*. Are the two blocks connected by allowed excitation
              INTERACT = 0
              IF(MXEXC.EQ.1) THEN
                 IF((ICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (ICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &               ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
C
                 IF(INTERACT.EQ.0) GOTO 10000

                 IF(NTEST.GE.100) THEN
                   WRITE(6,*) ' Next s block in batch : '
                   write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                   write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                 END IF
*
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            ((LLBSM.GT.LLASM).OR.
     &            (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &            GOTO 8764
*
                IF(NTEST.GE.60) THEN
                  WRITE(6,*) ' RSSBCB will be called for '
                  WRITE(6,*) ' Sigma block : '
                  WRITE(6,*) ' ISOFF ', ISOFF
                  WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' C     block : '
                  WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(6,*) ' ICOFF ', ICOFF
                  WRITE(6,*) ' Overall scale: 1.0D0'
                END IF
*
                IF(IRESTRICT.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     XFAC = 0.5D0*1.0D0
                ELSE
                     XFAC = 1.0D0
                END IF
*. Form of operator in action
                IF(IPERTOP.NE.0) THEN
*. Not exact Hamiltonian in use
                   IPTSPC = IH0SPC(IATP,IBTP)
                   JPTSPC = IH0SPC(JATP,JBTP)
                   IJOP   = IH0INSPC(IPTSPC)
*
C?                 WRITE(6,*) ' SBLOCKS:IPTSPC,JPTSPC,IJOP',
C?   &                         IPTSPC,JPTSPC,IJOP
*
                ELSE
                   IPTSPC = -1
                   JPTSPC = -1
                   IJOP   = -1
                END IF
*
                CALL RSSBCB2(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
     &               IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &               IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &               MXTSOB,SOMESCR,SOMEH,XIJILS)
*
                ISUMUP = 1
*
C?            IF(IATP.EQ.1.AND.IBTP.EQ.1) THEN
C?              WRITE(6,*) ' First element in updates S11 :',SB(1)
C?            END IF
C?            IF(JATP.EQ.1.AND.JBTP.EQ.1) THEN
C?              WRITE(6,*) ' First element of C ', CB(1)
C?            END IF
 8764         CONTINUE
              END IF
10000         CONTINUE
*./\  End of loop over sigma blocks
 8765         CONTINUE
            END IF
*.         /\ computation is only done for nonvanishing blocks!
*.        ==> SCALFAC = 1.0D0
 8888       CONTINUE
 8999       CONTINUE
 9000     CONTINUE
      sigmatime = sigmatime + mpi_wtime() - starttimer
*. /\ End of loop over C blocks in Batch
*. End of loop over S blocks in batch
20000   CONTINUE
*. End of loop over batches of C blocks
*
*. every process (MASTER + NODES) has its part of the
*. sigma vector block, but we want to have it only on
*. the on all nodes.
*. --> we sum the whole sigma vector block on the
*. every node writes its corresponding blcok on disk in 
*. RASSG3!
*
      INTERNALSIZE = 0
      IF(IODENSEPAR.eq.1) then
CSK      WRITE(6,*)'ISUMUP,MYPROC',ISUMUP,MYPROC
        starttimer = mpi_wtime()
*.      only nodes that need the sblock
        DO JSB = INSBLK, NSBLOCK+INSBLK-1
          IF(IBLOCKD(JSB).EQ.MYPROC) ISUMUP = 1
        END DO
*       node did a computation --> ISUMUP = 1! 
        IF(ISUMUP.GT.0)THEN
          JKEY = MYPROC
          JCOLOR = 2
        ELSE IF(ISUMUP.EQ.0)THEN
          JKEY = MYPROC
          JCOLOR = 3
        END IF
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,JCOLOR,JKEY,NREDGRP,IERR)
        CALL MPI_COMM_SIZE(NREDGRP,INTERNALSIZE,IERR)
        IF(ISUMUP.EQ.0) GOTO 21000
*. Reduce the sigma-vector to C2 with operation MPI_SUM 
*. use MPI_REAL8 as datatype: NDATATYPE ==  2 
*. reduce to all nodes:       NREDTYPE  == -1
        NDATATYPE = 2
        NREDTYPE = -1
*
        waittime = waittime + mpi_wtime() - starttimer
        starttimer = MPI_WTIME()
          CALL REDVEC(SB,C2,ILENSUM,NDATATYPE,MPI_SUM,NREDGRP,NREDTYPE)
CSK          INTERNALSIZE = NMPROC
CSK          CALL REDVEC(SB,C2,ILENSUM,2,MPI_SUM,MPI_COMM_WORLD,NREDTYPE)
        reducetime = reducetime + MPI_WTIME() - starttimer
*
21000   CONTINUE
CSKCSK        CALL MPI_COMM_FREE(NREDGRP,IERR)
        IF(ISUMUP.EQ.0) GOTO 22000
*. sigma-vecor is expected to reside on SB ...
        CALL DCOPY(ILENSUM,C2,1,SB,1)
      end if
C    /\ IODENSEPAR = 1
*. Order the fresh sigma-vector
      DO  ISBLK = 1 , NSBLOCK
        IF(ISBLOCK(1,ISBLK).GT.0) THEN
CSK          JSIGMABLOCK = INSBLK + ISBLK - 1
CSK          IF(IODENSEPAR.EQ.2.AND.MYPROC.NE.IBLOCKD(JSIGMABLOCK))
CSK     &      GOTO 21500
          IATP = ISBLOCK(1,ISBLK)
          IBTP = ISBLOCK(2,ISBLK)
          IASM = ISBLOCK(3,ISBLK)
          IBSM = ISBLOCK(4,ISBLK)
          ISOFF  = ISBLOCK(5,ISBLK)
          ISOFFP = ISBLOCK(6,ISBLK)
          NIA = NSSOA(IASM,IATP)
          NIB = NSSOB(IBSM,IBTP)
          IF(ICJKAIB.NE.0) THEN
*. Tranpose sigma block was obtained, transpose to obtain correct block
             CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &                   NSSOA(IASM,IATP),C2)
             CALL DCOPY(NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP),C2,1,
     &                  SB(ISOFF),1)
          END IF
          IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP)
     &    CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
*
21500     CONTINUE
        END IF
      END DO
22000   CONTINUE
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' output blocks from SBLOCKS for myproc',MYPROC
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
*
          WALLTID2  = SECTID(waittime)
          WALLTID = SECTID(reducetime)
          WALLTID3 = SECTID(readtime)
          WALLTID4 = SECTID(sigmatime)
          IF(MYPROC.EQ.MASTER)THEN
            WRITE(6,50000) WALLTID3
            WRITE(6,60000) WALLTID4
            IF(IODENSEPAR.EQ.1) THEN
              WRITE(6,40000) WALLTID2
              WRITE(6,30000) WALLTID,ILENSUM,INTERNALSIZE
            END IF
          ELSE
            WRITE(12,50000) WALLTID3
            WRITE(12,60000) WALLTID4
            IF(IODENSEPAR.EQ.1) THEN
              WRITE(12,40000) WALLTID2
              WRITE(12,30000) WALLTID,ILENSUM,INTERNALSIZE
            END IF
          ENDIF
      CALL QEXIT('SBLOC')
      RETURN
*
30000 FORMAT(' >>>  WALL TIME FOR REDUCING THE SVEC OF L W CPUs:
     & ',A,I12,I4)
40000 FORMAT(' >>>  WALL TIME FOR WAITING TO REDUCE THE S-VEC  :
     & ',A)
50000 FORMAT(' >>>  WALL TIME FOR READING AND SKIPPING C-BLOCKS:
     & ',A)
60000 FORMAT(' >>>  WALL TIME FOR COMPUTING THE S-VECTOR       :
     & ',A)
      END
#else
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_DUMMY3
      END
#endif
