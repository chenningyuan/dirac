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
#if defined(VAR_MPI2) && !defined(SYS_LOCALDISKS)
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_15(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU1,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU1,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_35(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU3,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU3,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_37(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU3,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU3,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU7,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU7,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_38(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU3,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU3,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU8,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU8,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_47(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU4,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU4,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU7,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU7,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_48(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU4,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU4,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU8,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU8,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_53(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU5,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU5,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU3,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU3,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_54(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU5,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU5,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU4,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU4,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_65(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU6,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU6,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_67(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU6,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU6,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU7,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU7,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_68(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU6,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU6,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU8,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU8,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_71(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU7,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU7,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU1,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU1,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_72(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU7,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU7,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU2,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU2,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CPVC_PARPAR_75(SCR,IBKLH,NBLK,IOFF1,IOFF2,
     &           IPCIO)
*
*. copy vector from parallel I/O file LUIN to
*. parallel I/O file LUOUT
*
*. IPCIO gt 0: USE collective File-IO !!!
*
*.              INPUT
*. IBKLH: array of blocklenth of all blocks
*. NBLK: total number of blocks
*. JOFFS: OFFSET in file LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION IBKLH(*),SCR(*)
      INTEGER IBLOCK,LBL
      INTEGER(KIND=MPI_OFFSET_KIND) ISTARTW,ISTARTR
*
      ISTARTR = IOFF1
      ISTARTW = IOFF2
*
      IBLOCK = 1
 1000 CONTINUE
*
      LBL = IBKLH(IBLOCK)
      IF(LBL.GT.0) THEN
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU7,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_READ_AT(ILU7,ISTARTR,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
        IF(IPCIO.GT.0) THEN
          CALL MPI_FILE_WRITE_AT_ALL(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        ELSE
          CALL MPI_FILE_WRITE_AT(ILU5,ISTARTW,SCR,LBL,MPI_REAL8,
     &         ISTAT,IERR)
        END IF
*
      END IF
C    /\ LBL gt 0
      ISTARTR = ISTARTR + LBL
      ISTARTW = ISTARTW + LBL
      IBLOCK = IBLOCK + 1
*
      IF(LBL.GE.0.AND.IBLOCK.LE.NBLK) GOTO 1000
* 
      END 
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_1(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU1,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_1B(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU1,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU8,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_2(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU2,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_2B(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU2,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU7,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_3(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU3,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_4(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU4,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_5(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU8,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PCIO_57(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU7,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PSIO_5(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT(ILU8,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PSIO_5B(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI,ISWAP)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C       
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
         CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
         IF(ISWAP.EQ.1)THEN
           CALL MPI_FILE_WRITE_AT(ILU6,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         ELSE IF(ISWAP.EQ.2)THEN
           CALL MPI_FILE_WRITE_AT(ILU7,IOFFSET2,SEGMNT,LBL,
     &          MPI_REAL8,ISTAT,IERR)
         END IF
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCD_PSIO_85(SCALE,SEGMNT,JOFFSET1,JOFFSET2,
     &           JBKINF,JBLOCKI)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C JOFFSET1/2 DETERMINES THE OFFSET FOR THE NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2
      INTEGER IBLCK
C
      IOFFSET1 = JOFFSET1
      IOFFSET2 = JOFFSET2
      IBLCK = 1
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
C
      LBL = JBKINF(IBLCK)
      IF ( LBL .GT. 0 ) THEN
C
CSK         write(6,*)'start reading at ilu8',myproc,ilu8
         CALL MPI_FILE_READ_AT(ILU8,IOFFSET1,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
         IF(LBL .GT. 0 )
     &   CALL SCALVE(SEGMNT,SCALE,LBL)
CSK         write(6,*)'start writing at ilu5',myproc,ilu5
         CALL MPI_FILE_WRITE_AT(ILU5,IOFFSET2,SEGMNT,LBL,
     &        MPI_REAL8,ISTAT,IERR)
*
         IOFFSET2 = IOFFSET2 + LBL
         IOFFSET1 = IOFFSET1 + LBL
*
      END IF
C
      IBLCK = IBLCK + 1
      IF( LBL .GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PCIO_1(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU1,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PCIO_1B(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU1,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU8,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU8,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PCIO_2(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU2,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PCIO_2B(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU2,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU7,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU7,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PCIO_3(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU3,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PCIO_4(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT_ALL(ILU4,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU5,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT_ALL(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PSIO_5(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT(ILU8,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT(ILU8,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMD_PSIO_5B(VEC1,VEC2,FAC1,FAC2,JOFFSET1,
     &           JOFFSET2,JOFFSET3,JBKINF,JBLOCKI,ISWAP)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C JOFFSET DETERMINES THE OFFSET FOR A NODE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
      INTEGER IBLCK
*
      IOFFSET1 = JOFFSET1 
      IOFFSET2 = JOFFSET2 
      IOFFSET3 = JOFFSET3 
      IBLCK = 1      
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        NBL1 = JBKINF(IBLCK)
        NBL2 = JBKINF(IBLCK)
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GT. 0 ) THEN
          CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,VEC1,NBL1,MPI_REAL8,
     &         ISTAT,IERR)
          IF(ISWAP.EQ.1) THEN
            CALL MPI_FILE_READ_AT(ILU6,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
            CALL MPI_FILE_READ_AT(ILU7,IOFFSET2,VEC2,NBL1,MPI_REAL8,
     &           ISTAT,IERR)
          END IF
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(ISWAP.EQ.1) THEN
           CALL MPI_FILE_WRITE_AT(ILU7,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          ELSE IF(ISWAP.EQ.2) THEN
           CALL MPI_FILE_WRITE_AT(ILU6,IOFFSET3,VEC1,NBL1,MPI_REAL8,
     &          ISTAT,IERR)
          END IF
        END IF
C
      IBLCK = IBLCK + 1 
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL1
      IOFFSET3 = IOFFSET3 + NBL1
      IF(NBL1.GE. 0 .AND. IBLCK .LE. JBLOCKI) GOTO 1000
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRDD_DPCIO_52(VEC1,VEC2,JBLKTL,JBKINF,
     &                JBLOCKI,JOFF,KOFF,IPCIO)
C
C DOUBLE DISC MPI-IO OF INPROD
C
C IBKINF CONTAINS BLOCK LENTGTH
C IBLOCKI DEFINES NUMBER OF BLOCKS
C JOFF DEFINES OFFSET IN FILE LU1
C KOFF DEFINES OFFSET IN FILE LU2
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,ISTART1,ISTART2
      LOGICAL DIFVEC
      INTEGER JCOUNT
C
      ICOUNT = 0
      X = 0.0D0
CSK ILU7 <--> ILU1
        DIFVEC =  .TRUE.
C
C FIND OFFSET FOR MPI-FILEs
C
C LU1
      IF(JOFF.GT.0) THEN
        ISTART1 = JBLKTL * JOFF
      ELSE IF(JOFF.EQ.0) THEN
        ISTART1 = 0
      END IF
C LU2
      IF(KOFF.GT.0) THEN
        ISTART2 = JBLKTL * KOFF
      ELSE IF(KOFF.EQ.0) THEN
        ISTART2 = 0
      END IF
C
C
C LOOP OVER BLOCKS OF VECTORS
C
      IBLKI = 1
      IOFFSET1 = 0
      IOFFSET1 = IOFFSET1 + ISTART1
      IOFFSET2 = 0
      IOFFSET2 = IOFFSET2 + ISTART2
 1000 CONTINUE
C
      NBL1 = JBKINF(IBLKI)
      IF(NBL1 .GT.0) THEN
        IF(IPCIO.GT.0) THEN
        CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        ELSE
        CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        END IF
        IF( DIFVEC) THEN
           NBL2 = JBKINF(IBLKI)
           IF(NBL2.GT.0)THEN
             IF(IPCIO.GT.0) THEN
               CALL MPI_FILE_READ_AT_ALL(ILU2,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             ELSE
               CALL MPI_FILE_READ_AT(ILU2,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             END IF
CSK               CALL MPI_GET_COUNT(ISTAT,MPI_REAL8,JCOUNT,IERR)
CSK               WRITE(6,*)'MASTER read',JCOUNT,'REAL*8'
           ENDIF
           IF(NBL1 .GT. 0 )
     &     X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL2
      IBLKI = IBLKI + 1
      IF(IBLKI.LE. JBLOCKI) GOTO 1000
C
      INPRDD_DPCIO_52 = X
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRDD_DPCIO_54(VEC1,VEC2,JBLKTL,JBKINF,
     &                JBLOCKI,JOFF,KOFF,IPCIO)
C
C DOUBLE DISC MPI-IO OF INPROD
C
C IBKINF CONTAINS BLOCK LENTGTH
C IBLOCKI DEFINES NUMBER OF BLOCKS
C JOFF DEFINES OFFSET IN FILE LU1
C KOFF DEFINES OFFSET IN FILE LU2
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,ISTART1,ISTART2
      LOGICAL DIFVEC
      INTEGER JCOUNT
C
      ICOUNT = 0
      X = 0.0D0
CSK ILU7 <--> ILU1
        DIFVEC =  .TRUE.
C
C FIND OFFSET FOR MPI-FILEs
C
C LU1
      IF(JOFF.GT.0) THEN
        ISTART1 = JBLKTL * JOFF
      ELSE IF(JOFF.EQ.0) THEN
        ISTART1 = 0
      END IF
C LU2
      IF(KOFF.GT.0) THEN
        ISTART2 = JBLKTL * KOFF
      ELSE IF(KOFF.EQ.0) THEN
        ISTART2 = 0
      END IF
C
C
C LOOP OVER BLOCKS OF VECTORS
C
      IBLKI = 1
      IOFFSET1 = 0
      IOFFSET1 = IOFFSET1 + ISTART1
      IOFFSET2 = 0
      IOFFSET2 = IOFFSET2 + ISTART2
 1000 CONTINUE
C
      NBL1 = JBKINF(IBLKI)
      IF(NBL1 .GT.0) THEN
        IF(IPCIO.GT.0) THEN
        CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        ELSE
        CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        END IF
        IF( DIFVEC) THEN
           NBL2 = JBKINF(IBLKI)
           IF(NBL2.GT.0)THEN
             IF(IPCIO.GT.0) THEN
               CALL MPI_FILE_READ_AT_ALL(ILU4,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             ELSE
               CALL MPI_FILE_READ_AT(ILU4,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             END IF
CSK               CALL MPI_GET_COUNT(ISTAT,MPI_REAL8,JCOUNT,IERR)
CSK               WRITE(6,*)'MASTER read',JCOUNT,'REAL*8'
           ENDIF
           IF(NBL1 .GT. 0 )
     &     X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL2
      IBLKI = IBLKI + 1
      IF(IBLKI.LE. JBLOCKI) GOTO 1000
C
      INPRDD_DPCIO_54 = X
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRDD_DPCIO_55(VEC1,VEC2,JBLKTL,JBKINF,
     &                JBLOCKI,JOFF,KOFF,IPCIO)
C
C DOUBLE DISC MPI-IO OF INPROD
C
C IBKINF CONTAINS BLOCK LENTGTH
C IBLOCKI DEFINES NUMBER OF BLOCKS
C JOFF DEFINES OFFSET IN FILE LU1
C KOFF DEFINES OFFSET IN FILE LU2
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,ISTART1,ISTART2
      LOGICAL DIFVEC
      INTEGER JCOUNT
C
      ICOUNT = 0
      X = 0.0D0
CSK ILU5 <--> ILU5
        DIFVEC =  .FALSE.
C
C FIND OFFSET FOR MPI-FILEs
C
C LU1
      IF(JOFF.GT.0) THEN
        ISTART1 = JBLKTL * JOFF
      ELSE IF(JOFF.EQ.0) THEN
        ISTART1 = 0
      END IF
C LU2
      IF(KOFF.GT.0) THEN
        ISTART2 = JBLKTL * KOFF
      ELSE IF(KOFF.EQ.0) THEN
        ISTART2 = 0
      END IF
C
C
C LOOP OVER BLOCKS OF VECTORS
C
      IBLKI = 1
      IOFFSET1 = 0
      IOFFSET1 = IOFFSET1 + ISTART1
      IOFFSET2 = 0
      IOFFSET2 = IOFFSET2 + ISTART2
 1000 CONTINUE
C
      NBL1 = JBKINF(IBLKI)
      IF(NBL1 .GT.0) THEN
        IF(IPCIO.GT.0) THEN
CSK        write(6,*)'start reading...',myproc
        CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        ELSE
CSK        write(6,*)'start reading...alone',myproc
        CALL MPI_FILE_READ_AT(ILU5,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        END IF
        IF( DIFVEC) THEN
           NBL2 = JBKINF(IBLKI)
           IF(NBL2.GT.0)THEN
             IF(IPCIO.GT.0) THEN
               CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             ELSE
               CALL MPI_FILE_READ_AT(ILU5,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             END IF
CSK               CALL MPI_GET_COUNT(ISTAT,MPI_REAL8,JCOUNT,IERR)
CSK               WRITE(6,*)'MASTER read',JCOUNT,'REAL*8'
           ENDIF
           IF(NBL1 .GT. 0 )
     &     X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL2
      IBLKI = IBLKI + 1
      IF(IBLKI.LE. JBLOCKI) GOTO 1000
C
      INPRDD_DPCIO_55 = X
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRDD_DPCIO_71(VEC1,VEC2,JBLKTL,JBKINF,
     &                JBLOCKI,JOFF,KOFF,IPCIO)
C
C DOUBLE DISC MPI-IO OF INPROD
C
C IBKINF CONTAINS BLOCK LENTGTH
C IBLOCKI DEFINES NUMBER OF BLOCKS
C JOFF DEFINES OFFSET IN FILE LU1
C KOFF DEFINES OFFSET IN FILE LU2
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,ISTART1,ISTART2
      LOGICAL DIFVEC
      INTEGER JCOUNT
C
      ICOUNT = 0
      X = 0.0D0
CSK ILU7 <--> ILU1
        DIFVEC =  .TRUE.
C
C FIND OFFSET FOR MPI-FILEs
C
C LU1
      IF(JOFF.GT.0) THEN
        ISTART1 = JBLKTL * JOFF
      ELSE IF(JOFF.EQ.0) THEN
        ISTART1 = 0
      END IF
C LU2
      IF(KOFF.GT.0) THEN
        ISTART2 = JBLKTL * KOFF
      ELSE IF(KOFF.EQ.0) THEN
        ISTART2 = 0
      END IF
C
C
C LOOP OVER BLOCKS OF VECTORS
C
      IBLKI = 1
      IOFFSET1 = 0
      IOFFSET1 = IOFFSET1 + ISTART1
      IOFFSET2 = 0
      IOFFSET2 = IOFFSET2 + ISTART2
 1000 CONTINUE
C
      NBL1 = JBKINF(IBLKI)
      IF(NBL1 .GT.0) THEN
        IF(IPCIO.GT.0) THEN
        CALL MPI_FILE_READ_AT_ALL(ILU7,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        ELSE
        CALL MPI_FILE_READ_AT(ILU7,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        END IF
        IF( DIFVEC) THEN
           NBL2 = JBKINF(IBLKI)
           IF(NBL2.GT.0)THEN
             IF(IPCIO.GT.0) THEN
               CALL MPI_FILE_READ_AT_ALL(ILU1,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             ELSE
               CALL MPI_FILE_READ_AT(ILU1,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             END IF
CSK               CALL MPI_GET_COUNT(ISTAT,MPI_REAL8,JCOUNT,IERR)
CSK               WRITE(6,*)'MASTER read',JCOUNT,'REAL*8'
           ENDIF
           IF(NBL1 .GT. 0 )
     &     X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL2
      IBLKI = IBLKI + 1
      IF(IBLKI.LE. JBLOCKI) GOTO 1000
C
      INPRDD_DPCIO_71 = X
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRDD_DPCIO_75(VEC1,VEC2,JBLKTL,JBKINF,
     &                JBLOCKI,JOFF,KOFF,IPCIO)
C
C DOUBLE DISC MPI-IO OF INPROD
C
C IBKINF CONTAINS BLOCK LENTGTH
C IBLOCKI DEFINES NUMBER OF BLOCKS
C JOFF DEFINES OFFSET IN FILE LU1
C KOFF DEFINES OFFSET IN FILE LU2
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,ISTART1,ISTART2
      LOGICAL DIFVEC
      INTEGER JCOUNT
C
      ICOUNT = 0
      X = 0.0D0
CSK ILU7 <--> ILU5
        DIFVEC =  .TRUE.
C
C FIND OFFSET FOR MPI-FILEs
C
C LU1
      IF(JOFF.GT.0) THEN
        ISTART1 = JBLKTL * JOFF
      ELSE IF(JOFF.EQ.0) THEN
        ISTART1 = 0
      END IF
C LU2
      IF(KOFF.GT.0) THEN
        ISTART2 = JBLKTL * KOFF
      ELSE IF(KOFF.EQ.0) THEN
        ISTART2 = 0
      END IF
C
C
C LOOP OVER BLOCKS OF VECTORS
C
      IBLKI = 1
      IOFFSET1 = 0
      IOFFSET1 = IOFFSET1 + ISTART1
      IOFFSET2 = 0
      IOFFSET2 = IOFFSET2 + ISTART2
 1000 CONTINUE
C
      NBL1 = JBKINF(IBLKI)
      IF(NBL1 .GT.0) THEN
        IF(IPCIO.GT.0) THEN
        CALL MPI_FILE_READ_AT_ALL(ILU7,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        ELSE
        CALL MPI_FILE_READ_AT(ILU7,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        END IF
        IF( DIFVEC) THEN
           NBL2 = JBKINF(IBLKI)
           IF(NBL2.GT.0)THEN
             IF(IPCIO.GT.0) THEN
               CALL MPI_FILE_READ_AT_ALL(ILU5,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             ELSE
               CALL MPI_FILE_READ_AT(ILU5,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             END IF
CSK               CALL MPI_GET_COUNT(ISTAT,MPI_REAL8,JCOUNT,IERR)
CSK               WRITE(6,*)'MASTER read',JCOUNT,'REAL*8'
           ENDIF
           IF(NBL1 .GT. 0 )
     &     X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL2
      IBLKI = IBLKI + 1
      IF(IBLKI.LE. JBLOCKI) GOTO 1000
C
      INPRDD_DPCIO_75 = X
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRDD_DPCIO_88(VEC1,VEC2,JBLKTL,JBKINF,
     &                JBLOCKI,JOFF,KOFF,IPCIO)
C
C DOUBLE DISC MPI-IO OF INPROD
C
C IBKINF CONTAINS BLOCK LENTGTH
C IBLOCKI DEFINES NUMBER OF BLOCKS
C JOFF DEFINES OFFSET IN FILE LU1
C KOFF DEFINES OFFSET IN FILE LU2
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,ISTART1,ISTART2
      LOGICAL DIFVEC
      INTEGER JCOUNT
C
      ICOUNT = 0
      X = 0.0D0
CSK ILU8 <--> ILU8
        DIFVEC =  .FALSE.
C
C FIND OFFSET FOR MPI-FILEs
C
C LU1
      IF(JOFF.GT.0) THEN
        ISTART1 = JBLKTL * JOFF
      ELSE IF(JOFF.EQ.0) THEN
        ISTART1 = 0
      END IF
C LU2
      IF(KOFF.GT.0) THEN
        ISTART2 = JBLKTL * KOFF
      ELSE IF(KOFF.EQ.0) THEN
        ISTART2 = 0
      END IF
C
C
C LOOP OVER BLOCKS OF VECTORS
C
      IBLKI = 1
      IOFFSET1 = 0
      IOFFSET1 = IOFFSET1 + ISTART1
      IOFFSET2 = 0
      IOFFSET2 = IOFFSET2 + ISTART2
 1000 CONTINUE
C
      NBL1 = JBKINF(IBLKI)
      IF(NBL1 .GT.0) THEN
        IF(IPCIO.GT.0) THEN
        CALL MPI_FILE_READ_AT_ALL(ILU8,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        ELSE
        CALL MPI_FILE_READ_AT(ILU8,IOFFSET1,VEC1,NBL1,
     &          MPI_REAL8,ISTAT,IERR)
        END IF
        IF( DIFVEC) THEN
           NBL2 = JBKINF(IBLKI)
           IF(NBL2.GT.0)THEN
             IF(IPCIO.GT.0) THEN
               CALL MPI_FILE_READ_AT_ALL(ILU8,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             ELSE
               CALL MPI_FILE_READ_AT(ILU8,IOFFSET2,VEC2,NBL2,
     &              MPI_REAL8,ISTAT,IERR)
             END IF
CSK               CALL MPI_GET_COUNT(ISTAT,MPI_REAL8,JCOUNT,IERR)
CSK               WRITE(6,*)'MASTER read',JCOUNT,'REAL*8'
           ENDIF
           IF(NBL1 .GT. 0 )
     &     X = X + INPROD(VEC1,VEC2,NBL1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
        ELSE
          IF(NBL1 .GT. 0 )
     &      X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IOFFSET1 = IOFFSET1 + NBL1
      IOFFSET2 = IOFFSET2 + NBL2
      IBLKI = IBLKI + 1
      IF(IBLKI.LE. JBLOCKI) GOTO 1000
C
      INPRDD_DPCIO_88 = X
C
      END
*
#else
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE IAMIRRELEVANT2
*. dummy routine for normal compilation       
       END 
#endif
