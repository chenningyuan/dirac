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
#if defined (VAR_MPI)
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MV7_PAR(C,HC,LUC,LUHC)
*
* Outer routine for sigma vector generation
* GAS version !!!!
*
* Written in terms of RASG3/SBLOCK, May 1997
* Works only for ICISTR .ne. 1, i.e. pure memory version
* has been eliminated.
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
      CALL QENTER('MV7  ')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBLOCK')

*
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' MV7 does not work for ICISTR = 1'
        WRITE(6,*) ' SWITCH to ICISTR = 2,3 or program'
        Call Abend2(' MV7 does not work for ICISTR = 1')
      END IF
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Arrays giving allowed type combinations
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
*. Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
*. Arrays for partitioning of sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*. Batches  of C vector
      IF(ICISTR.EQ.1) THEN
        NDET = 0
        LBLOCK = NDET
        write(6,*) 'NDET is undefined in mv7 '
        write(6,*) 'and ICISTR = 1 '
        Call Abend2('Quitting.')
      ELSE IF (ICISTR.EQ.2) THEN
        LBLOCK = MXSB
      ELSE IF (ICISTR.EQ.3) THEN
        LBLOCK = MXSOOB
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD)
*. Number of BLOCKS
        NBLOCK = IFRMR(WORK(KLSI1BT),1,NBATCH)
     &         + IFRMR(WORK(KLSLBT),1,NBATCH) - 1
C?      WRITE(6,*) ' Number of blocks ', NBLOCK

      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
      LLUC = LUC
      LLUHC = LUHC
*
* try something new --> read in scalfac only once!
* IRDSCF = 1 else IRDSCF = 0
      IRDSCF = 1
*
      IF(IRDSCF.GT.0) THEN
*
*. For scaling for each TTS block
      CALL MEMMAN(KISCLFAC,NBLOCK,'ADDL  ',2,'SSCLFC')
      CALL MEMMAN(KICBLINF,NBLOCK,'ADDL  ',1,'KITCLL')
*
      END IF 
*
       
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) 'ICISTR = 1 route currently not active'
        Call Abend2('MV7 : ICISTR = 1 route currently not active')
      ELSE
        CALL RASSG3(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &       WORK(KLSI1BT),WORK(KLSIBT),LLUC,LLUHC,NBLOCK,
     &       WORK(KISCLFAC),WORK(KICBLINF))
      END IF
*. Eliminate local memory
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'SBLOCK')
*
      CALL QEXIT('MV7  ')
*
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
