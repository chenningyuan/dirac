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
*
*============================================================*
*  Common LUCITA and LUCIAREL routines for parallel purposes *
*                                                            *
*   collected by Stefan Knecht, Jan 18 - ???, 2007           *
*============================================================*
*
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE CALC_OFF_MPI_FILE(FILENAME,IGROUPLIST,IBLOCKL,
     &                             IBLOCKD,IVEC_IN)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C
      DIMENSION IGROUPLIST(*), IBLOCKL(*), IBLOCKD(*)
      CHARACTER*6 FILENAME
      INTEGER(KIND=MPI_OFFSET_KIND) ND_VEC_IOFF, ND_VEC_IOFF2
      INTEGER(KIND=MPI_OFFSET_KIND) MULT1, MULT2, JVEC
      INTEGER(KIND=MPI_OFFSET_KIND) I_VEC_LEN_SCR
C
      I_VEC_LEN_SCR = 0
      JVEC  = IVEC_IN
      MULT1 = 0
      MULT2 = 0
C
      IF( NPTEST_VAR .ge. 10 ) 
     & WRITE(LUWRT,*) 'Calculating offset for file', FILENAME
C
C     count active blocks for each cpu, calculate individual offset
C
      ND_VEC_IOFF = 0
      ND_VEC_IOFF2 = 0
      NBLK_ACT = 0
C
C     scratch file for c-vector
C
      IF( FILENAME .eq. 'ILUC  ' ) THEN
        DO JBLK = 1, NUM_BLOCKS
          IF( IBLOCKL(JBLK) .gt. 0 ) THEN
            NBLK_ACT = NBLK_ACT + 1
          END IF
        END DO
C
Cveryold        MY_ACT_BLK_ALL = NBLK_ACT
C???        MY_ACT_BLK_ALL = NUM_BLOCKS * IRC_SAVE
        MY_ACT_BLK_ALL = NUM_BLOCKS
C
        GOTO 101
C
      END IF
C
      NBLK_ACT = 0
C
C
C     we have group files, so only cpus within a group are relevant
C
      MY_VEC1_IOFF = 0
      MY_VEC2_IOFF = 0
C
      DO IPROC = 1, NEWCOMM_PROC
C     
        JPROC = IGROUPLIST(IPROC)
C
        DO JBLK = 1, NUM_BLOCKS
          IF(IBLOCKD(JBLK) .eq. JPROC ) THEN
C?            WRITE(LUWRT,*) 'assigning a block to JPROC',JBLK,JPROC
            ND_VEC_IOFF = ND_VEC_IOFF + IBLOCKL(JBLK)
            IF( MYPROC .eq. JPROC ) THEN
              MY_VEC1_IOFF = MY_VEC1_IOFF + IBLOCKL(JBLK) 
              MY_VEC2_IOFF = MY_VEC2_IOFF + IBLOCKL(JBLK) 
              NBLK_ACT = NBLK_ACT + 1
            END IF
          END IF
        END DO
C
        IF( MYPROC .eq. JPROC ) THEN
          MULT1 = IRC_SAVE 
C
C         MY_VEC1_IOFF = length of one vector
C
          MY_VEC1_IOFF = MY_VEC1_IOFF
C
C         MY_VEC1_IOFF = length of one vector x complex multiplier ( x 1 = real!)
C
          I_VEC_LEN_SCR = MY_VEC2_IOFF
          MY_VEC2_IOFF = I_VEC_LEN_SCR * MULT1
C
          MY_ACT_BLK1 = NBLK_ACT
          MY_ACT_BLK2 = NBLK_ACT * MULT1

C
          IF( FILENAME .eq. 'LUDIA ' ) THEN
C
C           initialize MY_DIA_OFF
C
            MY_DIA_OFF = 0
C
            MY_DIA_OFF = ND_VEC_IOFF2 * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU1  ') THEN
C
C           initialize MY_LU1_OFF
C
            MY_LU1_OFF = 0
C
            MY_LU1_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU2  ') THEN
C
C           initialize MY_LU2_OFF
C
            MY_LU2_OFF = 0
C
            MY_LU2_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU3  ') THEN
C
C           initialize MY_LU3_OFF
C
            MY_LU3_OFF = 0
C
            MY_LU3_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU4  ') THEN
C
C           initialize MY_LU4_OFF
C
            MY_LU4_OFF = 0
C
            MY_LU4_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU5  ') THEN
C
C           initialize MY_LU5_OFF
C
            MY_LU5_OFF = 0
C
            MY_LU5_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          END IF
C        /\ filenames
        END IF
C
C       reset ND_VEC_IOFF2 and ND_VEC_IOFF values
C
        ND_VEC_IOFF2 = ND_VEC_IOFF
C
        NBLK_ACT = 0
C
      END DO
C
      IF( NPTESTVAR .ge. 10 ) THEN
        IF( FILENAME .eq. 'LUDIA ' ) THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'LUDIA offset :', MY_DIA_OFF
        ELSE IF( FILENAME .eq. 'ILU1  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU1  offset :', MY_LU1_OFF
        ELSE IF( FILENAME .eq. 'ILU2  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU2  offset :', MY_LU2_OFF
        ELSE IF( FILENAME .eq. 'ILU3  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU3  offset :', MY_LU3_OFF
        ELSE IF( FILENAME .eq. 'ILU4  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU4  offset :', MY_LU4_OFF
        ELSE IF( FILENAME .eq. 'ILU5  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU5  offset :', MY_LU5_OFF
        END IF
        WRITE(LUWRT,*) '  MY_VEC1_IOFF :', MY_VEC1_IOFF
        WRITE(LUWRT,*) '  MY_VEC2_IOFF :', MY_VEC2_IOFF
        WRITE(LUWRT,*) '  MY_ACT_BLK1  :', MY_ACT_BLK1
        WRITE(LUWRT,*) '  MY_ACT_BLK2  :', MY_ACT_BLK2
      END IF
C
 101  CONTINUE     
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV(LUIN,LUOUT,SEGMNT,NDISTND,NBLOCK,
     &                            JCOMM,IREW,LBLK,LUDUMPF)
*
* Read vector from disc file luin to array segmnt
* and send to master if necessary using communicator jcomm.
* vector is written in file luout
*
* Stefan Knecht, Jan07
* NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SEGMNT(*),NDISTND(NBLOCK)
C
*
      MYNSLV = 0
*
CSK      WRITE(LUDUMPF,*)'NDISTND(1)=',NDISTND(1)
CSK      WRITE(LUDUMPF,*)'LUOUT =',LUOUT
      IF( IREW .NE. 0 ) THEN
        IF (MYPROC.NE.MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (MYPROC.EQ.MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
          CALL REWINE( LUIN ,LBLK)
        END IF
      END IF
C
C LOOP OVER BLOCKS
C
       DO 500 II = 1, NBLOCK 
*. slaves first ...
         IF(MYPROC.NE.MASTER)THEN
           IF(NDISTND(II).NE.MYPROC)THEN 
             GOTO 500
           ELSE IF(NDISTND(II).EQ.MYPROC)THEN
             GOTO 100
           ENDIF
*. master follows ...
         ELSE IF(MYPROC.EQ.MASTER) THEN
           IF(NDISTND(II).NE.MYPROC)THEN
             IF(NDISTND(II).EQ.-2)THEN
               LBL = 0
               CALL ITODS(LBL,1,-1,LUOUT)
               IMZERO = 1
               IAMPACK = 0
               CALL ZERORC(LBL,LUOUT,IAMPACK)
               GOTO 500
             ELSE IF(NDISTND(II).GT.0)THEN
               MYNSLV = NDISTND(II) 
               GOTO 100 
             END IF
           ELSE IF(NDISTND(II).EQ.MYPROC)THEN
             CALL IFRMDS(LBL,1,-1,LUIN)
             CALL ITODS(LBL,1,-1,LUOUT)
             IF( LBL .GE. 0 ) THEN
                KBLK = -1
                NO_ZEROING = 1
                CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                       NO_ZEROING)
                IF(IMZERO.EQ.0) THEN
                  IF(IAMPACK.EQ.0) THEN
                   CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
                  ELSE
                   CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
                  END IF
                ELSE
                  CALL ZERORC(LBL,LUOUT,IAMPACK)
                END IF
*
             END IF
*           /\LBL GE 0
                GOTO 500
           END IF
*          /\NDISTND eq 0 or ne 0
         END IF
*. ... continue
 100     CONTINUE
          IF(MYPROC.NE.MASTER) THEN
            CALL IFRMDS(LBL,1,-1,LUIN)
            CALL MPI_SEND(LBL,1,my_MPI_INTEGER,MASTER,90,
     &                    JCOMM,IERR)
          END IF
*
*
          IF(MYPROC.EQ.MASTER) THEN
            CALL MPI_RECV(LBL,1,my_MPI_INTEGER,MYNSLV,90,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS(LBL,1,-1,LUOUT)
          END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          IF(MYPROC.NE.MASTER) THEN
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          CALL MPI_SEND(IMZERO,1,my_MPI_INTEGER,MASTER,91,
     &                  JCOMM,IERR)
          CALL MPI_SEND(IAMPACK,1,my_MPI_INTEGER,MASTER,92,
     &                  JCOMM,IERR)
          CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MASTER,93,
     &                  JCOMM,IERR)
          END IF
*
*
          IF(MYPROC.EQ.MASTER) THEN
             CALL MPI_RECV(IMZERO,1,my_MPI_INTEGER,MYNSLV,91,
     &                     JCOMM,ISTAT,IERR)
             CALL MPI_RECV(IAMPACK,1,my_MPI_INTEGER,MYNSLV,92,
     &                     JCOMM,ISTAT,IERR)
             CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MYNSLV,93,
     &                     JCOMM,ISTAT,IERR)
          END IF
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              IF(MYPROC.EQ.MASTER) THEN
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              END IF
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            IF(MYPROC.EQ.MASTER) THEN
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            END IF
          END IF
        END IF
 500    CONTINUE
        IF(MYPROC.EQ.MASTER)THEN
*. write an EOF-mark
          CALL ITODS(-1,1,-1,LUOUT)
        ENDIF
        CALL IFRMDS(LBL,1,-1,LUIN)
        IF(LBL.LT.0) WRITE(LUDUMPF,*)'I am fine',myproc
        IF(LBL.GE.0) THEN
          WRITE(6,*) 'Error in gathering of optimal c-vectors
     & detected',myproc
          CALL Abend2('I assume all blocks to be on the master 
     & but there are still some left')
        ENDIF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV_REL(LUIN,LUOUT,SEGMNT,NDISTND,NBLOCK,
     &                                LBLOCKL,JCOMM,LBLK,IROOT,
     &                                LUINLIST,MY_IOFF_LUIN,I_RUN_CPLX)
C
C     Written by  S. Knecht         - June 11 2007
C
C**********************************************************************
C
C     copy vector from MPI file LUIN to 'normal file' LUOUT blockwise
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SEGMNT(*), NDISTND( NBLOCK ), LBLOCKL(*)
      DIMENSION LUINLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_INT_IN    = 0
C
C     initialize slave tag 
C
      MYNSLV = 0
C
C
C     ================
C      COMPLEX VECTOR
C     ================
C
C     .................
C      REAL PART FIRST 
C     .................
C
C     set new offset
C
C     position in file is at the end of vector IROOT - 1
C
C     note: real part: --> MY_VEC2_IOFF and MY_ACT_BLK2
C
      IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                 ( IROOT - 1 )   * MY_VEC2_IOFF
      IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                 ( IROOT - 1 )   * MY_ACT_BLK2
C
C LOOP OVER BLOCKS
C
       DO 501 II = 1, NBLOCK 
C
C         slaves first ...
C
         IF( MYPROC .ne. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN 
C
             IOFFSET_SCRATCH = 0
             NUM_BLK = 0
C
             GOTO 500
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
             GOTO 100
           ENDIF
C
C        MASTER follows ...
C
         ELSE IF( MYPROC .eq. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN
             IF(NDISTND(II) .eq. -2 )THEN
C
C              block has zero length
C
               LBL = 0
               CALL ITODS(LBL,1,-1,LUOUT)
               IMZERO = 1
               IAMPACK = 0
               CALL ZERORC(LBL,LUOUT,IAMPACK)
C
               IOFFSET_SCRATCH = 0
               NUM_BLK = 0
C
               GOTO 500
             ELSE IF( NDISTND(II) .gt. 0 )THEN
C
C              recieve block from node MYNSLV
C
               MYNSLV = NDISTND(II) 
               GOTO 100 
             END IF
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
C
C            MASTER block
C
             LBL = LBLOCKL(II)
C  
             CALL ITODS(LBL,1,-1,LUOUT)
C
             IF( LUINLIST( IOFFSET_INT_IN ) .gt. 0 ) THEN
C
C             read block in
               CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                               MPI_REAL8,ISTAT,IERR)
               IMZERO = 0
             ELSE
               IMZERO = 1
             END IF
C
             IF(IMZERO.EQ.0) THEN
                CALL TODSC (SEGMNT,LBL,-1,LUOUT)
             ELSE
                CALL ZERORC(LBL,LUOUT,0)
             END IF
C
             NUM_BLK = 1
             IOFFSET_SCRATCH = LBL
C
             GOTO 500
           END IF
*         /\NDISTND eq 0 or ne 0
        END IF
C
C  ... continue
C
 100    CONTINUE
C
C       slave block
C
        LBL = LBLOCKL(II)
C
        IF( MYPROC .EQ. MASTER ) THEN
C
           CALL ITODS(LBL,1,-1,LUOUT)
C
        END IF
C
        IF(MYPROC.NE.MASTER) THEN
C
          ISENDTOM = 0
          ISENDTOM = LUINLIST( IOFFSET_INT_IN )
C
          CALL MPI_SEND(ISENDTOM,1,my_MPI_INTEGER,MASTER,90,
     &                  JCOMM,IERR)
C
          IF( ISENDTOM .gt. 0 ) THEN
C
            CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                            MPI_REAL8,ISTAT,IERR)
            CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MASTER,91,
     &                    JCOMM,IERR)
C
          END IF
          IOFFSET_SCRATCH = LBL
          NUM_BLK = 1
        END IF
C
C
        IF(MYPROC.EQ.MASTER) THEN
C
          ISENDTOM = 0
C
          CALL MPI_RECV(ISENDTOM,1,my_MPI_INTEGER,MYNSLV,90,
     &                  JCOMM,ISTAT,IERR)
C
          IF( ISENDTOM .gt. 0 ) THEN
C
            CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MYNSLV,91,
     &                     JCOMM,ISTAT,IERR)
            CALL TODSC(SEGMNT,LBL,-1,LUOUT)
          ELSE
            CALL ZERORC(LBL,LUOUT,0)
          END IF
          IOFFSET_SCRATCH = 0
          NUM_BLK = 0
        END IF
C
 500    CONTINUE
C       
C       keep track of correct offset
        IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
        IOFFSET_INT_IN  = IOFFSET_INT_IN  + NUM_BLK
 501  CONTINUE
C
      IF( MYPROC. eq. MASTER )THEN
C       write an EOF-mark
        CALL ITODS(-1,1,-1,LUOUT)
      ENDIF
C
C     ..............
C        IMAG PART 
C     ..............
C
      IOFFSET_SCRATCH = 0
      NUM_BLK = 0
C
      IF( I_RUN_CPLX .eq. 2 )THEN
C
C       set new offset
C
C       position in file is at the end of vector IROOT - 1
C
C       note: complex part: --> MY_VEC1_IOFF and MY_ACT_BLK1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )  * MY_VEC2_IOFF    +
     &                                    MY_VEC1_IOFF
        IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                   ( IROOT - 1 )  * MY_ACT_BLK2 + MY_ACT_BLK1
C
C LOOP OVER BLOCKS
C
       DO 601 II = 1, NBLOCK 
C
C         slaves first ...
C
         IF( MYPROC .ne. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN 
C
             IOFFSET_SCRATCH = 0
             NUM_BLK = 0
C
             GOTO 600
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
             GOTO 102
           ENDIF
C
C        MASTER follows ...
C
         ELSE IF( MYPROC .eq. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN
             IF(NDISTND(II) .eq. -2 )THEN
C
C              block has zero length
C
               LBL = 0
               CALL ITODS(LBL,1,-1,LUOUT)
               IMZERO = 1
               IAMPACK = 0
               CALL ZERORC(LBL,LUOUT,IAMPACK)
C
               IOFFSET_SCRATCH = 0
               NUM_BLK = 0
C
               GOTO 600
             ELSE IF( NDISTND(II) .gt. 0 )THEN
C
C              recieve block from node MYNSLV
C
               MYNSLV = NDISTND(II) 
               GOTO 102 
             END IF
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
C
C            MASTER block
C
             LBL = LBLOCKL(II)
C  
             CALL ITODS(LBL,1,-1,LUOUT)
C
             IF( LUINLIST( IOFFSET_INT_IN ) .gt. 0 ) THEN
C
C             read block in
               CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                               MPI_REAL8,ISTAT,IERR)
               IMZERO = 0
             ELSE
               IMZERO = 1
             END IF
C
             IF( IMZERO .eq. 0 ) THEN
                CALL TODSC (SEGMNT,LBL,-1,LUOUT)
             ELSE
                CALL ZERORC(LBL,LUOUT,0)
             END IF
C
             NUM_BLK = 1
             IOFFSET_SCRATCH = LBL
C
             GOTO 600
           END IF
*          ^NDISTND eq 0 or ne 0
          END IF
C
C     ... continue
C
 102      CONTINUE
C
C         slave block
C
          LBL = LBLOCKL(II)
C
          IF( MYPROC .eq. MASTER ) THEN
C
             CALL ITODS(LBL,1,-1,LUOUT)
C  
          END IF
C
          IF( MYPROC .ne. MASTER ) THEN
C
            ISENDTOM = 0
            ISENDTOM = LUINLIST( IOFFSET_INT_IN )
C
            CALL MPI_SEND(ISENDTOM,1,my_MPI_INTEGER,MASTER,93,
     &                    JCOMM,IERR)
C
            IF( ISENDTOM .gt. 0 ) THEN
C
              CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                              MPI_REAL8,ISTAT,IERR)
              CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MASTER,94,
     &                      JCOMM,IERR)
C
            END IF
            IOFFSET_SCRATCH = LBL
            NUM_BLK = 1
          END IF
C
C
          IF( MYPROC .eq. MASTER ) THEN
C
            ISENDTOM = 0
C
            CALL MPI_RECV(ISENDTOM,1,my_MPI_INTEGER,MYNSLV,93,
     &                    JCOMM,ISTAT,IERR)
C
            IF( ISENDTOM .gt. 0 ) THEN
C
              CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MYNSLV,94,
     &                       JCOMM,ISTAT,IERR)
              CALL TODSC(SEGMNT,LBL,-1,LUOUT)
            ELSE
              CALL ZERORC(LBL,LUOUT,0)
            END IF
            IOFFSET_SCRATCH = 0
            NUM_BLK = 0
          END IF
C
 600      CONTINUE
C       
C         keep track of correct offset
          IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
          IOFFSET_INT_IN  = IOFFSET_INT_IN  + NUM_BLK
 601    CONTINUE
C
        IF( MYPROC. eq. MASTER )THEN
C         write an EOF-mark
          CALL ITODS(-1,1,-1,LUOUT)
        ENDIF
C
      END IF
C     ^ COMPLEX PART NEEDED?
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV_DIA_REL1(LUIN,LUOUT,SEGMNT,NDISTND,
     &                                     NBLOCK,LBLOCKL,JCOMM,LBLK,
     &                                     IROOT,MY_IOFF_LUIN,
     &                                     I_RUN_CPLX)
C
C     Written by  S. Knecht         - June 22 2007
C
C**********************************************************************
C
C     copy H diagonal from MPI file LUIN to 'normal file' LUOUT blockwise
C
C     NOTE: IROOT = IROOT
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SEGMNT(*), NDISTND( NBLOCK ), LBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
C
C     initialize slave tag 
C
      MYNSLV = 0
C
C
C     ================
C      COMPLEX VECTOR
C     ================
C
C     .................
C      REAL PART FIRST 
C     .................
C
C     set new offset
C
C     position in file is at the end of vector IROOT - 1
C
C     note: real part: --> MY_VEC2_IOFF
C
      IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                 ( IROOT - 1 )   * MY_VEC2_IOFF
C
C LOOP OVER BLOCKS
C
       DO 501 II = 1, NBLOCK 
C
C         slaves first ...
C
         IF( MYPROC .ne. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN 
C
             IOFFSET_SCRATCH = 0
             NUM_BLK = 0
C
             GOTO 500
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
             GOTO 100
           ENDIF
C
C        MASTER follows ...
C
         ELSE IF( MYPROC .eq. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN
             IF(NDISTND(II) .eq. -2 )THEN
C
C              block has zero length
C
               LBL = 0
               CALL ITODS(LBL,1,-1,LUOUT)
               IAMPACK = 0
               CALL ZERORC(LBL,LUOUT,IAMPACK)
C
               IOFFSET_SCRATCH = 0
C
               GOTO 500
             ELSE IF( NDISTND(II) .gt. 0 )THEN
C
C              recieve block from node MYNSLV
C
               MYNSLV = NDISTND(II) 
               GOTO 100 
             END IF
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
C
C            MASTER block
C
             LBL = LBLOCKL(II)
C  
             CALL ITODS(LBL,1,-1,LUOUT)
C
C            read block in
             CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
C
             CALL TODSC(SEGMNT,LBL,-1,LUOUT)
C
             IOFFSET_SCRATCH = LBL
C
             GOTO 500
           END IF
*         /\NDISTND eq 0 or ne 0
        END IF
C
C  ... continue
C
 100    CONTINUE
C
C       slave block
C
        LBL = LBLOCKL(II)
C
        IF( MYPROC .EQ. MASTER ) THEN
C
           CALL ITODS(LBL,1,-1,LUOUT)
C
        END IF
C
        IF(MYPROC.NE.MASTER) THEN
C
          CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                          MPI_REAL8,ISTAT,IERR)
          CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MASTER,91,
     &                  JCOMM,IERR)
C
          IOFFSET_SCRATCH = LBL
        END IF
C
C
        IF(MYPROC.EQ.MASTER) THEN
C
          CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MYNSLV,91,
     &                  JCOMM,ISTAT,IERR)
          CALL TODSC(SEGMNT,LBL,-1,LUOUT)
C
          IOFFSET_SCRATCH = 0
        END IF
C
 500    CONTINUE
C       
C       keep track of correct offset
        IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
 501  CONTINUE
C
      IF( MYPROC. eq. MASTER )THEN
C       write an EOF-mark
        CALL ITODS(-1,1,-1,LUOUT)
      ENDIF
C
C     ..............
C      COMPLEX PART 
C     ..............
C
      IOFFSET_SCRATCH = 0
C
      IF( I_RUN_CPLX .eq. 2 )THEN
C
C       set new offset
C
C       position in file is at the end of vector IROOT - 1
C
C       note: complex part: --> MY_VEC1_IOFF and MY_ACT_BLK1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )   * MY_VEC2_IOFF + 
     &                                     MY_VEC1_IOFF
C
C LOOP OVER BLOCKS
C
       DO 601 II = 1, NBLOCK 
C
C         slaves first ...
C
         IF( MYPROC .ne. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN 
C
             IOFFSET_SCRATCH = 0
C
             GOTO 600
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
             GOTO 102
           ENDIF
C
C        MASTER follows ...
C
         ELSE IF( MYPROC .eq. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN
             IF(NDISTND(II) .eq. -2 )THEN
C
C              block has zero length
C
               LBL = 0
               CALL ITODS(LBL,1,-1,LUOUT)
               IAMPACK = 0
               CALL ZERORC(LBL,LUOUT,IAMPACK)
C
               IOFFSET_SCRATCH = 0
C
               GOTO 600
             ELSE IF( NDISTND(II) .gt. 0 )THEN
C
C              recieve block from node MYNSLV
C
               MYNSLV = NDISTND(II) 
               GOTO 102 
             END IF
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
C
C            MASTER block
C
             LBL = LBLOCKL(II)
C  
             CALL ITODS(LBL,1,-1,LUOUT)
C
C            read block in
             CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
             CALL TODSC(SEGMNT,LBL,-1,LUOUT)
C
             IOFFSET_SCRATCH = LBL
C
             GOTO 600
           END IF
*          ^NDISTND eq 0 or ne 0
          END IF
C
C     ... continue
C
 102      CONTINUE
C
C         slave block
C
          LBL = LBLOCKL(II)
C
          IF( MYPROC .eq. MASTER ) THEN
C
             CALL ITODS(LBL,1,-1,LUOUT)
C  
          END IF
C
          IF( MYPROC .ne. MASTER ) THEN
C
            CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                            MPI_REAL8,ISTAT,IERR)
            CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MASTER,94,
     &                    JCOMM,IERR)
C
            IOFFSET_SCRATCH = LBL
          END IF
C
C
          IF( MYPROC .eq. MASTER ) THEN
C
            CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MYNSLV,94,
     &                    JCOMM,ISTAT,IERR)
            CALL TODSC(SEGMNT,LBL,-1,LUOUT)
C
            IOFFSET_SCRATCH = 0
          END IF
C
 600      CONTINUE
C       
C         keep track of correct offset
          IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
 601    CONTINUE
C
        IF( MYPROC. eq. MASTER )THEN
C         write an EOF-mark
          CALL ITODS(-1,1,-1,LUOUT)
        ENDIF
C
      END IF
C     ^ COMPLEX PART NEEDED?
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV_DIA_REL2(LUIN,LUOUT,SEGMNT,NDISTND,
     &                                     NBLOCK,LBLOCKL,JCOMM,LBLK,
     &                                     IROOT,MY_IOFF_LUIN)
C
C     Written by  S. Knecht         - June 22 2007
C
C**********************************************************************
C
C     copy H diagonal from MPI file LUIN to 'normal file' LUOUT blockwise
C
C     NOTE: IROOT = IROOT
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SEGMNT(*), NDISTND( NBLOCK ), LBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
C
C     initialize slave tag 
C
      MYNSLV = 0
C
C     set new offset
C
C     position in file is at the end of vector IROOT - 1
C
      IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                 ( IROOT - 1 )   * MY_VEC2_IOFF
C
C LOOP OVER BLOCKS
C
       DO 501 II = 1, NBLOCK 
C
C         slaves first ...
C
         IF( MYPROC .ne. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN 
C
             IOFFSET_SCRATCH = 0
             NUM_BLK = 0
C
             GOTO 500
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
             GOTO 100
           ENDIF
C
C        MASTER follows ...
C
         ELSE IF( MYPROC .eq. MASTER )THEN
           IF( NDISTND(II) .ne. MYPROC )THEN
             IF(NDISTND(II) .eq. -2 )THEN
C
C              block has zero length
C
               LBL = 0
               CALL ITODS(LBL,1,-1,LUOUT)
               IAMPACK = 0
               CALL ZERORC(LBL,LUOUT,IAMPACK)
C
               IOFFSET_SCRATCH = 0
C
               GOTO 500
             ELSE IF( NDISTND(II) .gt. 0 )THEN
C
C              recieve block from node MYNSLV
C
               MYNSLV = NDISTND(II) 
               GOTO 100 
             END IF
           ELSE IF( NDISTND(II) .eq. MYPROC )THEN
C
C            MASTER block
C
             LBL = LBLOCKL(II)
C  
             CALL ITODS(LBL,1,-1,LUOUT)
C
C            read block in
             CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
C
             CALL TODSC(SEGMNT,LBL,-1,LUOUT)
C
             IOFFSET_SCRATCH = LBL
C
             GOTO 500
           END IF
*         /\NDISTND eq 0 or ne 0
        END IF
C
C  ... continue
C
 100    CONTINUE
C
C       slave block
C
        LBL = LBLOCKL(II)
C
        IF( MYPROC .EQ. MASTER ) THEN
C
           CALL ITODS(LBL,1,-1,LUOUT)
C
        END IF
C
        IF(MYPROC.NE.MASTER) THEN
C
          CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,SEGMNT,LBL,
     &                          MPI_REAL8,ISTAT,IERR)
          CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MASTER,91,
     &                  JCOMM,IERR)
C
          IOFFSET_SCRATCH = LBL
        END IF
C
C
        IF(MYPROC.EQ.MASTER) THEN
C
          CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MYNSLV,91,
     &                  JCOMM,ISTAT,IERR)
          CALL TODSC(SEGMNT,LBL,-1,LUOUT)
C
          IOFFSET_SCRATCH = 0
        END IF
C
 500    CONTINUE
C       
C       keep track of correct offset
        IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
 501  CONTINUE
C
      IF( MYPROC. eq. MASTER )THEN
C       write an EOF-mark
        CALL ITODS(-1,1,-1,LUOUT)
      ENDIF
C
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV2(LUIN,LUOUT,SEGMNT,NDISTND,NBLOCK,
     &                             JCOMM,IREW,LBLK,LUDUMPF)
*
* Read vector from disc file luin to array segmnt
* and send to nodes if necessary using communicator jcomm.
* vector is written in file luout
*
* Stefan Knecht, Jan07
* NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SEGMNT(*),NDISTND(NBLOCK)
C
*
      MYNSLV = 0
*
      IF( IREW .NE. 0 ) THEN
        IF (MYPROC.NE.MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        ELSE IF (MYPROC.EQ.MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
          CALL REWINE( LUIN ,LBLK)
        END IF
      END IF
C
C LOOP OVER BLOCKS
C
       DO 500 II = 1, NBLOCK 
*. slaves first ...
         IF(MYPROC.NE.MASTER)THEN
           IF(NDISTND(II).NE.MYPROC)THEN 
             GOTO 500
           ELSE IF(NDISTND(II).EQ.MYPROC)THEN
             GOTO 100
           ENDIF
*. master follows ...
         ELSE IF(MYPROC.EQ.MASTER) THEN
           IF(NDISTND(II).NE.MYPROC)THEN
             IF(NDISTND(II).EQ.-2)THEN
               CALL IFRMDS(LBL,1,-1,LUIN) 
               IF(LBL.NE.0) THEN
                 WRITE(6,*)'Error in restart detected'
                 CALL Abend2('COPVCD_PAR_BDRIV2: Input vector is not
     &                       compatible to current calculation')
               ENDIF
               KBLK = -1
               NO_ZEROING = 1
               CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                      NO_ZEROING)
               GOTO 500
             ELSE IF(NDISTND(II).GT.0)THEN
               MYNSLV = NDISTND(II) 
               GOTO 100 
             END IF
           ELSE IF(NDISTND(II).EQ.MYPROC)THEN
             CALL IFRMDS(LBL,1,-1,LUIN)
             CALL ITODS(LBL,1,-1,LUOUT)
             IF( LBL .GE. 0 ) THEN
                KBLK = -1
                NO_ZEROING = 1
                CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                       NO_ZEROING)
                IF(IMZERO.EQ.0) THEN
                  IF(IAMPACK.EQ.0) THEN
                   CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
                  ELSE
                   CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
                  END IF
                ELSE
                  CALL ZERORC(LBL,LUOUT,IAMPACK)
                END IF
*
             END IF
*           /\LBL GE 0
                GOTO 500
           END IF
*          /\NDISTND eq 0 or ne 0
         END IF
*. ... continue
 100     CONTINUE
          IF(MYPROC.EQ.MASTER) THEN
            CALL IFRMDS(LBL,1,-1,LUIN)
            CALL MPI_SEND(LBL,1,my_MPI_INTEGER,MYNSLV,90,
     &                    JCOMM,IERR)
          END IF
*
*
          IF(MYPROC.NE.MASTER) THEN
            CALL MPI_RECV(LBL,1,my_MPI_INTEGER,MASTER,90,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUOUT)
          END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          IF(MYPROC.EQ.MASTER) THEN
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
            CALL MPI_SEND(IMZERO,1,my_MPI_INTEGER,MYNSLV,91,
     &                  JCOMM,IERR)
            CALL MPI_SEND(IAMPACK,1,my_MPI_INTEGER,MYNSLV,92,
     &                    JCOMM,IERR)
            CALL MPI_SEND(SEGMNT,LBL,MPI_REAL8,MYNSLV,93,
     &                    JCOMM,IERR)
            GOTO 500
          ENDIF
*
*
          IF(MYPROC.NE.MASTER) THEN
             CALL MPI_RECV(IMZERO,1,my_MPI_INTEGER,MASTER,91,
     &                     JCOMM,ISTAT,IERR)
             CALL MPI_RECV(IAMPACK,1,my_MPI_INTEGER,MASTER,92,
     &                     JCOMM,ISTAT,IERR)
             CALL MPI_RECV(SEGMNT,LBL,MPI_REAL8,MASTER,93,
     &                     JCOMM,ISTAT,IERR)
          ENDIF
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              IF(MYPROC.NE.MASTER) THEN
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              END IF
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            IF(MYPROC.NE.MASTER) THEN
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            END IF
          END IF
        END IF
 500    CONTINUE
*. write an EOF-mark
        CALL ITODS(-1,1,-1,LUOUT)
*. MASTER looks for the end of the input file
        IF(MYPROC.EQ.MASTER) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          IF(LBL.LT.0) WRITE(LUDUMPF,*)'I am fine'
          IF(LBL.GE.0) THEN
            WRITE(6,*) 'Error in restart detected'
            CALL Abend2('I assume all blocks to be distributed among the
     & nodes but there are still some left')
          ENDIF
        ENDIF  
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV3(SEGMNT,IBLOCKD,SCALFAC,
     &                             SCALFAC_LOCAL,IBLOCKL,NBLOCK,JCOMM,
     &                             LUIN,LUOUT,IREW,LBLK)
*
* Read blocks of the vector from disc file luin to array segment
* and broadcast to nodes that need the block, save on diskfile luout 
* if necessary. 
* Using communicator jcomm.
* OUTPUT: Blocks on LUOUT, modified SCALFAC_LOCAL.
*
* Stefan Knecht, Odense, Feb07
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SEGMNT(*),IBLOCKD(NBLOCK),IBLOCKL(NBLOCK)
      DIMENSION SCALFAC(NBLOCK),SCALFAC_LOCAL(NBLOCK)
*
      IF( IREW .NE. 0 ) THEN
          CALL REWINE( LUOUT ,LBLK)
          CALL REWINE( LUIN ,LBLK)
      END IF
C
C LOOP OVER BLOCKS
C
      DO IBLK = 1, NBLOCK 
C        check if block is active
        IF(SCALFAC(IBLK).NE.0.0D0)THEN
*
          IF(IBLOCKD(IBLK).EQ.MYPROC)THEN
*
            CALL IFRMDS(LBL,1,-1,LUIN)
CSK            IF(LBLK .GE.0 ) THEN
CSK              KBLK = LBL
CSK            ELSE
              KBLK = -1
CSK            END IF
            NO_ZEROING = 1
*
C           read block in memory
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
C           check for consistency
            IF(LBL.NE.IBLOCKL(IBLK)) THEN
              WRITE(6,*)'Error detected in subroutine COPVCD_PAR_BDRIV3.
     &        Length of block in memory does not fit with the listed 
     &        block length, so I stop!'
               CALL Abend2('Error in subroutine COPVCD_PAR_BDRIV3!')
            END IF
*
          END IF
*
          LBL = IBLOCKL(IBLK)
*
C         broadcast the nonvanishing block
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,IBLOCKD(IBLK),
     &                   JCOMM,IERR)
C         some test writing, enable for debugging!        
CSK          IF(MYPROC.NE.IBLOCKD(IBLK)) THEN
CSK            WRITE(6,*)'THIS IS WHAT I GOT',MYPROC
CSK            CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
CSK          END IF
C         transfer the recieved block to disk if necessary
          IF(SCALFAC_LOCAL(IBLK).NE.0.0D0) THEN 
            KBLK = -1
            CALL ITODS(LBL,1,-1,LUOUT)
            CALL TODSC(SEGMNT,LBL,KBLK,LUOUT)
          END IF
*
        ELSE 
*
          IF(IBLOCKD(IBLK).EQ.MYPROC)THEN
C           skip over this block
            CALL IFRMDS(LBL,1,-1,LUIN)
            NO_ZEROING = 1
            KBLK = -1
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
          END IF
C         set local scalfac to 0.0D0
          SCALFAC_LOCAL(IBLK) = 0.0D0
*
        END IF
*
      END DO
C     write an EOF-mark
      CALL ITODS(-1,1,-1,LUOUT)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV4(SEGMNT,IBLOCKD,SCALFAC,
     &                             SCALFAC_LOCAL,IBLOCKL,NBLOCK,JCOMM,
     &                             LUIN,IREW,LBLK)
*
* Read blocks of the vector from disc file luin to array segment
* and broadcast to nodes that need the block, save on diskfile luout 
* if necessary. 
* Using communicator jcomm.
* OUTPUT: Blocks on LUOUT, modified SCALFAC_LOCAL.
*
* Stefan Knecht, Odense, Feb07
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SEGMNT(*),IBLOCKD(NBLOCK),IBLOCKL(NBLOCK)
      DIMENSION SCALFAC(NBLOCK),SCALFAC_LOCAL(NBLOCK)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
*
      IF(IREW.NE.0) THEN
        CALL REWINE(LUIN,LBLK)
      END IF
C     file offset for each cpu is 0
      IOFFSET = 0
C
C LOOP OVER BLOCKS
C
      DO IBLK = 1, NBLOCK 
C        check if block is active
        IF(SCALFAC(IBLK).NE.0.0D0)THEN
*
          IF(IBLOCKD(IBLK).EQ.MYPROC)THEN
*
            CALL IFRMDS(LBL,1,-1,LUIN)
CSK            IF(LBLK .GE.0 ) THEN
CSK              KBLK = LBL
CSK            ELSE
              KBLK = -1
CSK            END IF
            NO_ZEROING = 1
*
C           read block in memory
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
C           check for consistency
            IF(LBL.NE.IBLOCKL(IBLK)) THEN
              WRITE(6,*)'Error detected in subroutine COPVCD_PAR_BDRIV4.
     &        Length of block in memory does not fit with the listed 
     &        block length, so I stop!'
               CALL Abend2('Error in subroutine COPVCD_PAR_BDRIV4!')
            END IF
*
          END IF
*
          LBL = IBLOCKL(IBLK)
*
C         broadcast the nonvanishing block
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,IBLOCKD(IBLK),
     &                   JCOMM,IERR)
C         some test writing, enable for debugging!        
CSK          IF(MYPROC.NE.IBLOCKD(IBLK)) THEN
CSK            WRITE(6,*)'THIS IS WHAT I GOT',MYPROC
CSK            CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
CSK          END IF
C         transfer the recieved block to disk if necessary
          IF(SCALFAC_LOCAL(IBLK).NE.0.0D0) THEN 
            KBLK = -1
            CALL MPI_FILE_WRITE_AT(MPILUSCR,IOFFSET,SEGMNT,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
            IOFFSET = IOFFSET + LBL
          END IF
*
        ELSE 
*
          IF(IBLOCKD(IBLK).EQ.MYPROC)THEN
C           skip over this block
            CALL IFRMDS(LBL,1,-1,LUIN)
            NO_ZEROING = 1
            KBLK = -1
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
          END IF
C         set local scalfac to 0.0D0
          SCALFAC_LOCAL(IBLK) = 0.0D0
*
        END IF
*
      END DO
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV5(SEGMNT,IBLOCKD,SCALFAC,SCALFAC_GROUP,
     &                             IBLOCKL,NBLOCK,JCOMM,
     &                             IGROUPLIST,IPROCLIST)
C
C     Read blocks of the vector from disc file luin to array segment
C     and broadcast to nodes that need the block, save on diskfile luout 
C     if necessary. 
C     Using communicator jcomm.
C     OUTPUT: Blocks on LUOUT, modified SCALFAC_LOCAL.
C
C     Stefan Knecht, Odense, Feb07
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SEGMNT(*),IBLOCKD(NBLOCK),IBLOCKL(NBLOCK)
      DIMENSION SCALFAC(NBLOCK),SCALFAC_GROUP(NBLOCK)
      DIMENSION IGROUPLIST(NMPROC),IPROCLIST(NMPROC)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER IGROUPBLK,ROOTPROC
C
      ROOTPROC = 0
C     file offset for each cpu is 0
      IOFFSET_IN  = 0
      IOFFSET_OUT = 0
C
C LOOP OVER BLOCKS
C
csk      write(6,*)'hurray ho, i am here!',myproc
      DO IBLK = 1, NBLOCK 
C       get length of block
        LBL = IBLOCKL(IBLK)
C        check if block is active
        IF( SCALFAC(IBLK) .ne. 0.0D0 ) THEN
C         check if block belongs to own group
          IGROUPBLK = 0
          ISCRNODE = IBLOCKD(IBLK)
          DO IPROC = 1, NEWCOMM_PROC
            IF( ISCRNODE .eq. IGROUPLIST(IPROC) ) IGROUPBLK = 1
          END DO
          IF( IGROUPBLK .eq. 1 ) THEN
C
CSKT            WRITE(6,*) 'ACTIVE BLOCK ',IBLK
CSKT            WRITE(6,*) 'START READING AT IOFFSET_IN = ',IOFFSET_IN
            CALL MPI_FILE_READ_AT(MPILUSCR2,IOFFSET_IN,SEGMNT,LBL,
     &                            MPI_REAL8,ISTAT,IERR)
C
          END IF
C
C         find the sending root processor
          ROOTPROC = IPROCLIST(ISCRNODE+1) - 1
C
C         decide whether it is necessary to take part in communication
C
          NNCOLOR = 0
          NNKEY = 0
          NNKEY = MYPROC + 1
          IF( MYPROC .eq. ROOTPROC ) NNKEY = 0
          IF( SCALFAC_GROUP(IBLK) .ne. 0.0D0 ) THEN
            NNCOLOR = 7
          ELSE
            NNCOLOR = 8
          END IF
          CALL MPI_COMM_SPLIT(JCOMM,NNCOLOR,NNKEY,IBLOCKCOMM,IERR)
C
          IF( NNCOLOR .eq. 8 ) GOTO 200
C         broadcast the nonvanishing block
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,0,
     &                   IBLOCKCOMM,IERR)
CSK          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,ROOTPROC,
CSK     &                   JCOMM,IERR)
C         some test writing, enable for debugging!
CSK          IF(MYPROC.NE.IBLOCKD(IBLK)) THEN
CSK          END IF
C         transfer the recieved block to disk if necessary
          IF( SCALFAC_GROUP(IBLK) .ne. 0.0D0 ) THEN 
CSKT            WRITE(6,*)'THIS IS WHAT I GOT',ICOMM_ID
CSKT            CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
CSKT            WRITE(6,*) 'IOFFSET_OUT = ',IOFFSET_OUT
            CALL MPI_FILE_WRITE_AT(MPILUSCR,IOFFSET_OUT,SEGMNT,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
          END IF
C
 200      CONTINUE
          CALL MPI_COMM_FREE(IBLOCKCOMM,IERR)
C
        ELSE 
C         set local scalfac to 0.0D0
          SCALFAC_GROUP(IBLK) = 0.0D0
        END IF
        IOFFSET_OUT = IOFFSET_OUT + LBL
        IOFFSET_IN = IOFFSET_IN + LBL
C
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_BDRIV5_REL(LUIN,LUOUT,SEGMNT,IBLOCKD,
     &                                 ISCALFAC,ISCALFAC_GROUP,
     &                                 IBLOCKL,NBLOCK,JCOMM,
     &                                 IGROUPLIST,IPROCLIST,IRILP)
C
C     Read blocks of the vector from disc file LUIN to array segment
C     and broadcast to nodes that need the block, save on diskfile luout 
C     if necessary. 
C     communication via communicator jcomm.
C
C     OUTPUT: update of LUOUT with new TTSS-c-blocks, 
C             modified ISCALFAC_LOCAL.
C
C     based on the corresponding spinfree-routine
C
C     extended for an imaginary part of the c-vector
C     if NZ == IRILP == 2           S. Knecht - June 27 2007
C
C     Written by  S. Knecht         - May 23 2007
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SEGMNT(*), IBLOCKD(NBLOCK), IBLOCKL(NBLOCK)
      DIMENSION ISCALFAC(*), ISCALFAC_GROUP(*)
      DIMENSION IGROUPLIST(NMPROC),IPROCLIST(NMPROC)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCRATCH_P
      INTEGER IGROUPBLK, ROOTPROC, IOFF_BLOCK
C
      IOFF_BLOCK     = 0
      IOFF_SCRATCH_P = 0
CSK      WRITE(LUWRT,*) ' MY GROUPLIST resp PROCLIST'
CSK      CALL IWRTMAMN(IGROUPLIST,1,NMPROC,1,NMPROC,LUWRT)
CSK      CALL IWRTMAMN(IPROCLIST,1,NMPROC,1,NMPROC,LUWRT)
C
C     =========
C     REAL PART
C     =========
C
      ROOTPROC   = 0
C
C     file offset for each cpu is 0 - real part
C
      IOFFSET_IN  = 0
      IOFFSET_OUT = 0
C
C     LOOP OVER BLOCKS
C
      DO IBLK = 1, NBLOCK 
C
C       get length of block
C
        LBL = IBLOCKL(IBLK)
C
C       check if block is active
C
        IF( ISCALFAC(IBLK) .ne. 0 ) THEN
C
C         check if block belongs to own group
C
          IGROUPBLK = 0
          ISCRNODE = IBLOCKD(IBLK)
C
          DO IPROC = 1, NEWCOMM_PROC
            IF( ISCRNODE .eq. IGROUPLIST(IPROC) ) IGROUPBLK = 1
          END DO
          IF( IGROUPBLK .eq. 1 ) THEN
C
CSK            WRITE(LUWRT,*) 'ACTIVE BLOCK ',IBLK
CSK            WRITE(LUWRT,*) 'START READING AT IOFFSET_IN = ',IOFFSET_IN
            CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN,SEGMNT,LBL,
     &                            MPI_REAL8,ISTAT,IERR)
C
          END IF
C
C         find the sending root processor
CSK          write(LUWRT,*) ' processor list:myproc,mynew_id,icomm_id',
CSK     &    myproc,mynew_id,icomm_id
CSK          CALL IWRTMAMN(IPROCLIST,1,NMPROC,1,NMPROC,LUWRT)
CSK          CALL MPI_BARRIER(ICOMM,IERR)
          ROOTPROC = IPROCLIST(ISCRNODE+1) - 1
CSK          write(LUWRT,*)' ROOT proc should send: ICOMM_ID, ROOTPROC',
CSK     &    ICOMM_ID,ROOTPROC
C
C         decide whether it is necessary to take part in communication
C
          NNCOLOR = 0
          NNKEY = 0
          NNKEY = MYPROC + 1
          IF( MYPROC .eq. ROOTPROC ) NNKEY = 0
          IF( ISCALFAC_GROUP(IBLK) .ne. 0 ) THEN
            NNCOLOR = 7
          ELSE
            NNCOLOR = 8
          END IF
          CALL MPI_COMM_SPLIT(JCOMM,NNCOLOR,NNKEY,IBLOCKCOMM,IERR)
C
          IF( NNCOLOR .eq. 8 ) GOTO 200
C         broadcast the nonvanishing block
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,0,IBLOCKCOMM,IERR)
C         transfer the recieved block to disk if necessary
          IF( ISCALFAC_GROUP(IBLK) .ne. 0 ) THEN 
CSK            WRITE(LUWRT,*)'THIS IS WHAT I GOT',ICOMM_ID
CSK            CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,LUWRT)
CSK            WRITE(LUWRT,*) 'IOFFSET_OUT = ',IOFFSET_OUT
            CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_OUT,SEGMNT,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
          END IF
C
 200      CONTINUE
          CALL MPI_COMM_FREE(IBLOCKCOMM,IERR)
C
        ELSE 
C         set local iscalfac to 0
          ISCALFAC_GROUP(IBLK) = 0
        END IF
C
C       keep track of correct offset
        IOFFSET_OUT = IOFFSET_OUT + LBL
        IOFFSET_IN  = IOFFSET_IN  + LBL
C
      END DO
C     ^ loop over blocks - real part
C
C
      IF( IRILP .eq. 1 ) GOTO 1001
C
C     =========
C     IMAG PART
C     =========
C
      ROOTPROC = 0
C
C     pointer for file offset: skip real part of the c-vector
C
      IOFF_SCRATCH_P = L_COMBI
C
      IOFFSET_IN  = IOFF_SCRATCH_P
      IOFFSET_OUT = IOFF_SCRATCH_P
C     set correct block counter for imag part
      IOFF_BLOCK = NBLOCK
CSK      WRITE(LUWRT,*) 'IOFF_BLOCK is',IOFF_BLOCK
C
C     LOOP OVER BLOCKS
C
      DO IBLK = 1, NBLOCK 
C
C       get length of block
C
        LBL = IBLOCKL( IBLK )
C
C       check if block is active
C
        IF( ISCALFAC( IBLK + IOFF_BLOCK ) .ne. 0 ) THEN
C
CSK          WRITE(LUWRT,*) 'ACTIVE BLOCK imag part',IBLK, 
CSK     &                    IBLK + IOFF_BLOCK
C         check if block belongs to own group
C
          IGROUPBLK = 0
          ISCRNODE = IBLOCKD( IBLK )
C
          DO IPROC = 1, NEWCOMM_PROC
            IF( ISCRNODE .eq. IGROUPLIST( IPROC ) ) IGROUPBLK = 1
          END DO
          IF( IGROUPBLK .eq. 1 ) THEN
C
CSK            WRITE(LUWRT,*) 'START READING AT IOFFSET_IN imag part = ',
CSK     &                      IOFFSET_IN
            CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN,SEGMNT,LBL,
     &                            MPI_REAL8,ISTAT,IERR)
C
          END IF
C
C         find the sending root processor
CSK          write(LUWRT,*) ' processor list:myproc,mynew_id,icomm_id',
CSK     &    myproc,mynew_id,icomm_id
CSK          CALL IWRTMAMN(IPROCLIST,1,NMPROC,1,NMPROC,LUWRT)
CSK          CALL MPI_BARRIER(ICOMM,IERR)
          ROOTPROC = IPROCLIST(ISCRNODE+1) - 1
CSK          write(LUWRT,*)' ROOT proc should send: ICOMM_ID, ROOTPROC',
CSK     &    ICOMM_ID,ROOTPROC
C
C         decide whether it is necessary to take part in communication
C
          NNCOLOR = 0
          NNKEY = 0
          NNKEY = MYPROC + 1
          IF( MYPROC .eq. ROOTPROC ) NNKEY = 0
          IF( ISCALFAC_GROUP( IBLK + IOFF_BLOCK ) .ne. 0 ) THEN
            NNCOLOR = 7
          ELSE
            NNCOLOR = 8
          END IF
          CALL MPI_COMM_SPLIT(JCOMM,NNCOLOR,NNKEY,IBLOCKCOMM,IERR)
C
          IF( NNCOLOR .eq. 8 ) GOTO 300
C         broadcast the nonvanishing block
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,0,IBLOCKCOMM,IERR)
C         transfer the recieved block to disk if necessary
CSK          IF( ISCALFAC_GROUP( IBLK + IOFF_BLOCK ) .ne. 0 ) THEN 
CSK            WRITE(LUWRT,*)'THIS IS WHAT I GOT from ROOTPROC',MYPROC,
CSK     &                     ROOTPROC
CSK            CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,LUWRT)
CSK            WRITE(LUWRT,*) 'IOFFSET_OUT imag part = ',IOFFSET_OUT
CSK            WRITE(LUWRT,*) 'block to write imag part',IBLK,
CSK     &                      IBLK + IOFF_BLOCK
          CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_OUT,SEGMNT,LBL,
     &                           MPI_REAL8,ISTAT,IERR)
CSK          END IF
C
 300      CONTINUE
          CALL MPI_COMM_FREE(IBLOCKCOMM,IERR)
C
        ELSE 
C         set local iscalfac to 0
          ISCALFAC_GROUP( IBLK + IOFF_BLOCK ) = 0
        END IF
C
C       keep track of correct offset
        IOFFSET_OUT = IOFFSET_OUT + LBL
        IOFFSET_IN  = IOFFSET_IN  + LBL
C
      END DO
C     ^ loop over blocks - imag part
C
 1001 CONTINUE
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVC_MPITOS_BDRV(LUIN,LUOUT,SCR,NBLK_ACT,
     &                             NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                             MY_IOFF_LUIN)
C
C     copy vector from parallel I/O file LUIN to 
C     sequential file LUOUT
C     
C     active blocks on the MPI-file are flagged by a nonzero length
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
C     Scratch
      DIMENSION SCR(*)
      DIMENSION LBATS(*), LEBATS(*) ,I1BATS(*), IBATS(8,*)
C
      IOFFSET_SCRATCH = 0
C
C     loop over batches
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(SCR,LEBATS(ISBATCH))

        IOFFSET_IN_LUIN = MY_IOFF_LUIN + IOFFSET_SCRATCH
C
        CALL RDVEC_BATCH_DRV(SCR,LBATS(ISBATCH),LEBATS(ISBATCH),
     &                       IBATS(1,I1BATS(ISBATCH)),
     &                       LUIN,IOFFSET_IN_LUIN)
C
CSK        WRITE(6,*) 'THIS IS WHAT I WILL WRITE on LUOUT FOR
CSK     &  LUIN',LUIN,LUOUT,MYPROC
CSK        CALL WRTMAT(SCR,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
        DO ISBLK = I1BATS(ISBATCH),I1BATS(ISBATCH) + LBATS(ISBATCH) -1
C
          IOFF =  IBATS(6,ISBLK)
          ILEN  = IBATS(8,ISBLK)
          CALL ITODS(ILEN,1,-1,LUOUT)
          CALL TODSC(SCR(IOFF),ILEN,-1,LUOUT)
CSK          write(6,*) 'my IOFF and ILEN for ISBLK',IOFF,ILEN,ISBLK,MYPROC
C
        END DO
C       update offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATS(ISBATCH)
      END DO
C
C     write EOF mark
C
      CALL ITODS(-1,1,-1,LUOUT)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVC_MPITOS_BDRV2(LUIN,LUOUT,SCR,NBLK_ACT,
     &                             NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                             MY_IOFF_LUIN,MY_INT_IOFF_LUIN,
     &                             IVCOFF_IN)
C
C     copy vector from parallel I/O file LUIN to 
C     sequential file LUOUT
C     
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
C     Scratch
      DIMENSION SCR(*), IVCOFF_IN(*)
      DIMENSION LBATS(*), LEBATS(*) ,I1BATS(*), IBATS(8,*)
      INTEGER NUM_BLK
C
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
C
C     loop over batches
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(SCR,LEBATS(ISBATCH))
C
        IOFFSET_IN_LUIN = MY_IOFF_LUIN + IOFFSET_SCRATCH
        IOFFSET_IN_INT_LUIN = MY_INT_IOFF_LUIN + NUM_BLK
C
        CALL RDVEC_BATCH_DRV2(SCR,LBATS(ISBATCH),LEBATS(ISBATCH),
     &                        IBATS(1,I1BATS(ISBATCH)),
     &                        LUIN,IOFFSET_IN_LUIN,
     &                        IOFFSET_IN_INT_LUIN,IVCOFF_IN)
C
CSK        WRITE(6,*) 'THIS IS WHAT I WILL WRITE on LUOUT FOR
CSK     &  LUIN',LUIN,LUOUT,MYPROC
CSK        CALL WRTMAT(SCR,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
        DO ISBLK = I1BATS(ISBATCH),I1BATS(ISBATCH) + LBATS(ISBATCH) -1
C
          IOFF =  IBATS(6,ISBLK)
          ILEN  = IBATS(8,ISBLK)
          IACTIVE = IVCOFF_IN( MY_INT_IOFF_LUIN + ISBLK - 1 )
          CALL ITODS(ILEN,1,-1,LUOUT)
          CALL TODSC2(SCR(IOFF),ILEN,-1,LUOUT,IACTIVE)
CSK          write(6,*) 'my IOFF and ILEN for ISBLK',IOFF,ILEN,ISBLK,MYPROC
C
        END DO
C       update offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATS(ISBATCH)
        NUM_BLK = NUM_BLK + LBATS(ISBATCH) 
      END DO
C
C     write an EOF mark
C
      CALL ITODS(-1,1,-1,LUOUT)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVC_STOMPI_BDRV(LUIN,LUOUT,SCR,NBLK_ACT,
     &                             NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                             MY_IOFF_LUOUT)
C
C     copy vector from parallel I/O file LUIN to 
C     sequential file LUOUT
C     
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUOUT
C     Scratch
      DIMENSION SCR(*)
      DIMENSION LBATS(*), LEBATS(*) ,I1BATS(*), IBATS(8,*)
      INTEGER ILEN, NO_CHECK
      ILEN = 0
C
      IOFFSET_SCRATCH = 0
C
C     loop over batches
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(SCR,LEBATS(ISBATCH))
C       read blocks in
CSK        write(6,*) 'nbatch,I1BATS(ISBATCH),LBATS(ISBATCH) for myproc',
CSK     &              nbatch,I1BATS(ISBATCH),LBATS(ISBATCH),myproc
        ICBLK = 0
        DO ISBLK = I1BATS(ISBATCH),I1BATS(ISBATCH) + LBATS(ISBATCH) -1
C
          ICBLK = ICBLK + 1
          IOFF =  IBATS(6,ISBLK)
          CALL IFRMDS(ILEN,1,-1,LUIN)
          KBLK = -1
          NO_ZEROING = 0
          CALL FRMDSC2(SCR(IOFF),ILEN,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
        END DO

        IOFFSET_IN_LUOUT = MY_IOFF_LUOUT + IOFFSET_SCRATCH
C
          NO_CHECK = 1
C
C
CSK        IF(MYPROC.EQ.MASTER+1) THEN
CSK          DO ISB = 1, ICBLK
CSK            write(6,*) 'these are the blocks',ICBLK,INFO_ACT_BLK(ISB)
CSK          END DO
CSK        write(6,*) 'IOFFSET_IN_INT_LUOUT',IOFFSET_IN_INT_LUOUT
CSK        write(6,*) 'i should transfer this to file for batch',
CSK     &                LUOUT,isbatch
CSK          CALL WRTMAT(SCR,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CSK        END IF
        CALL WTVEC_BATCH_DRV(SCR,LBATS(ISBATCH),LEBATS(ISBATCH),
     &                       IBATS(1,I1BATS(ISBATCH)),
     &                       LUOUT,IOFFSET_IN_LUOUT,NO_CHECK)
C
C       update offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATS(ISBATCH)
      END DO
C
C     find EOF mark
C
      CALL IFRMDS(ILEN,1,-1,LUIN)
      IF(ILEN.GE.0) THEN
        WRITE(6,*) 'Error in COPVCD_STOMPI detected'
        CALL Abend2('I assume all blocks to be copied
     & but there are still some left')
      ENDIF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVC_STOMPI_BDRV2(LUIN,LUOUT,SCR,NBLK_ACT,
     &                              NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                              MY_IOFF_LUOUT,ICTRLARR,JCTRLOFF)
C
C     copy vector from parallel I/O file LUIN to 
C     sequential file LUOUT
C     
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUOUT
C     Scratch
      DIMENSION SCR(*),ICTRLARR(*)
      DIMENSION LBATS(*), LEBATS(*) ,I1BATS(*), IBATS(8,*)
      INTEGER NUM_BLK,ILEN, NO_CHECK
      ILEN = 0
C
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
C
C     loop over batches
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(SCR,LEBATS(ISBATCH))
C       read blocks in
CSK        write(6,*) 'nbatch,I1BATS(ISBATCH),LBATS(ISBATCH) for myproc',
CSK     &              nbatch,I1BATS(ISBATCH),LBATS(ISBATCH),myproc
        DO ISBLK = I1BATS(ISBATCH),I1BATS(ISBATCH) + LBATS(ISBATCH) -1
C
          IOFF =  IBATS(6,ISBLK)
          CALL IFRMDS(ILEN,1,-1,LUIN)
          KBLK = -1
          NO_ZEROING = 0
          CALL FRMDSC2(SCR(IOFF),ILEN,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          IF( IMZERO .eq. 0 ) THEN
            ICTRLARR(JCTRLOFF+ISBLK-1) = ILEN
          ELSE
            ICTRLARR(JCTRLOFF+ISBLK-1) = 0
          END IF
        END DO

        IOFFSET_IN_LUOUT = MY_IOFF_LUOUT + IOFFSET_SCRATCH
        IOFFSET_IN_INT_LUOUT = JCTRLOFF + NUM_BLK
C
CSK        write(6,*) 'IOFFSET_IN_INT_LUOUT',IOFFSET_IN_INT_LUOUT
CSK        write(6,*) 'i should transfer this to file for batch',
CSK     &                LUOUT,isbatch
CSK          CALL WRTMAT(SCR,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CSK        END IF
C
        NO_CHECK = 1
C
        CALL WTVEC_BATCH_DRV2(SCR,LBATS(ISBATCH),LEBATS(ISBATCH),
     &                        IBATS(1,I1BATS(ISBATCH)),LUOUT,
     &                        IOFFSET_IN_LUOUT,IOFFSET_IN_INT_LUOUT,
     &                        ICTRLARR,NO_CHECK)
C
C       update offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATS(ISBATCH)
        NUM_BLK = NUM_BLK + LBATS(ISBATCH) 
      END DO
C
C     find EOF mark
C
      CALL IFRMDS(ILEN,1,-1,LUIN)
      IF(ILEN.GE.0) THEN
        WRITE(6,*) 'Error in COPVCD_STOMPI detected'
        CALL Abend2('I assume all blocks to be copied
     & but there are still some left')
      ENDIF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPY_STOMPI_DRV1(LU1,LU2,LU3,LU4,VEC1,NBLK_ACT,
     &                            NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                            NROOTS,ICOPY,MY_IOFF_125,MY_IOFF_347,
     &                            MY_VEC_IOFF)
C
C     copy vector from parallel I/O file LUIN to 
C     sequential file LUOUT
C     
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*)
      DIMENSION LBATS(*), LEBATS(*) ,I1BATS(*), IBATS(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_125, MY_IOFF_347
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN_SPECIAL
C
        DO IRT = 1 , NROOTS
C
          MY_IOFF_LUIN_SPECIAL = MY_IOFF_125 +
     &                           ( MY_VEC_IOFF * (IRT - 1) )
C
          CALL COPVC_STOMPI_BDRV(LU1,LU1,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL)
C
          MY_IOFF_LUIN_SPECIAL = MY_IOFF_125 +
     &                         ( MY_VEC_IOFF * (IRT - 1) )
C
          CALL COPVC_STOMPI_BDRV(LU2,LU2,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL)
C
        END DO
C
C       continue with LU3 and LU4
C
        DO ICP = 1, ICOPY
C
             MY_IOFF_LUIN_SPECIAL = MY_IOFF_347 +
     &                              ( MY_VEC_IOFF * (ICP - 1) )
C
          CALL COPVC_STOMPI_BDRV(LU3,LU3,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL)
C
             MY_IOFF_LUIN_SPECIAL = MY_IOFF_347 +
     &                              ( MY_VEC_IOFF * (ICP - 1) )
C
          CALL COPVC_STOMPI_BDRV(LU4,LU4,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL)
C
        END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPY_STOMPI_DRV2(LU1,LU2,LU3,LU4,VEC1,NBLK_ACT,
     &                            NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                            NROOTS,ICOPY,MY_IOFF_125,MY_IOFF_347,
     &                            MY_VEC_IOFF,MY_INT_IOFF_125,
     &                            MY_INT_IOFF_347,ICTRL_VC1,ICTRL_VC2,
     &                            ICTRL_VC3,ICTRL_VC4)
C
C     copy vector from parallel I/O file LUIN to 
C     sequential file LUOUT
C     
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*)
      DIMENSION LBATS(*), LEBATS(*) ,I1BATS(*), IBATS(8,*)
      DIMENSION ICTRL_VC1(*), ICTRL_VC2(*), ICTRL_VC3(*), ICTRL_VC4(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_125, MY_IOFF_347
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN_SPECIAL
C
        ICTRLOFF1 = 0
        ICTRLOFF2 = 0
        ICTRLOFF3 = 0
        ICTRLOFF4 = 0
C
        DO IRT = 1 , NROOTS
C
          MY_IOFF_LUIN_SPECIAL = MY_IOFF_125 +
     &                           ( MY_VEC_IOFF * (IRT - 1) )
          ICTRLOFF1 = ( NBLK_ACT * (IRT - 1) ) + 1
C
          CALL COPVC_STOMPI_BDRV2(LU1,LU1,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL,
     &                           ICTRL_VC1,ICTRLOFF1)
C
          MY_IOFF_LUIN_SPECIAL = MY_IOFF_125 +
     &                         ( MY_VEC_IOFF * (IRT - 1) )
          ICTRLOFF2 = ( NBLK_ACT * (IRT - 1) ) + 1
C
          CALL COPVC_STOMPI_BDRV2(LU2,LU2,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL,
     &                           ICTRL_VC2,ICTRLOFF2)
C
        END DO
C
C       continue with LU3 and LU4
C
        DO ICP = 1, ICOPY
C
             MY_IOFF_LUIN_SPECIAL = MY_IOFF_347 +
     &                              ( MY_VEC_IOFF * (ICP - 1) )
             ICTRLOFF3 = ( NBLK_ACT * (ICP - 1) ) + 1
C
          CALL COPVC_STOMPI_BDRV2(LU3,LU3,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL,
     &                           ICTRL_VC3,ICTRLOFF3)
C
             MY_IOFF_LUIN_SPECIAL = MY_IOFF_347 +
     &                              ( MY_VEC_IOFF * (ICP - 1) )
             ICTRLOFF4 = ( NBLK_ACT * (ICP - 1) ) + 1
C
          CALL COPVC_STOMPI_BDRV2(LU4,LU4,VEC1,NBLK_ACT,
     &                           NBATCH,LBATS,LEBATS,
     &                           I1BATS,IBATS,
     &                           MY_IOFF_LUIN_SPECIAL,
     &                           ICTRL_VC4,ICTRLOFF4)
C
        END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND(NDIM,NBLOCKL,NPARBLOCK,NPARBLKWT)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/clunit.inc"
*
      DIMENSION NBLOCKL(NDIM), NPARBLKWT(2,NMPROC)
*. Scratch
      INTEGER INONVAN,MAXLBL,LBL,IWEIGHT,LABEL,IBLOCKN,IMINNP
      INTEGER MXSZTMP,MXNUMB,ITEMPL,ITEMPD
      INTEGER*8 ITOTBLL 
*. ======
*. Output
*. ======
*. array NPARBLOCK that contains the final block distribution
      DIMENSION NPARBLOCK(NDIM)
*. 
      MXSIZE = 1000000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1000000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      INONVAN = 0
      ITOTBLL = 0
      IMINNP = 0
      MAXLBL = 0 
      LABEL = 0
      MIN2 = -2
      NZERO = 0
      NTEST = 0
      CALL ISETVC(NPARBLOCK,MIN2,NDIM)
      DO 100 II = 1, NDIM
        LBL = NBLOCKL(II)
        IF(LBL.GT.0) THEN
          INONVAN = INONVAN + 1
          NPARBLOCK(II) = -1
          ITOTBLL = ITOTBLL + LBL
        ENDIF
        IF(LBL.GT.MAXLBL) THEN
          MAXLBL = LBL
          LABEL = II
        END IF
 100  CONTINUE
*
*
CSK      NTEST = 10
      IF(NTEST.GT.0) THEN
       IF(MYPROC.EQ.MASTER) THEN
         WRITE(6,*) '  total number of processes to distribute on:',
     &                 NMPROC
         WRITE(6,*) '  total number of active blocks:', INONVAN
         WRITE(6,*) '  overall active block length:  ', ITOTBLL 
         WRITE(6,*) '  largest active block:         ', LABEL
         WRITE(6,*) '  size of largest active block: ', MAXLBL
       END IF
      END IF
*
      IF(INONVAN.LT.NMPROC) THEN
        IMINNP = MIN(INONVAN,NMPROC)
        IF(NTEST.GT.0) THEN
          IF(MYPROC.EQ.MASTER) THEN
            WRITE(6,*) 'subroutine DISTBLKND speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP+1,'procs'
          END IF
        END IF
        IF(MYPROC.GE.INONVAN) GOTO 1001
      ELSE
        IMINNP = NMPROC
      END IF
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
 200  CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.INONVAN) THEN
*
        MXSZTMP = 0
*
        DO 300 II = 1, NDIM
*
          ITEMPL = NBLOCKL(II)
          ITEMPD = NPARBLOCK(II)
          ITEMPN = II
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GE.MXSZTMP) THEN
                  MXSZTMP = ITEMPL
                  MXNUMB = II
                END IF
              END IF
            END IF
          END IF
*
300     CONTINUE
*
        DO 400 IPR = 1, IMINNP
*
          ITEMPW2 = NPARBLKWT(2,IPR)
* attention: IAMTPROC = 1 <--> MYPROC = 0 !
          IAMTPROC = IPR
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
            ITEMPW1 = ITEMPW2
            IAMGPROC = IAMTPROC
          END IF
*
 400    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*. calculation of block MXNUMB by proc IAMGPROC
        IF(NTEST.GT.0) THEN
          IF(MYPROC.EQ.MASTER) THEN
            WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,
     &      IAMGPROC-1
          END IF
        ENDIF
*. raising NPARBLKWT(1,proc) by 1
*. adding on NPARBLKWT(2,proc) the weight of the new block
        MXSIZE = MXSZTMP
        NPARBLKWT(1,IAMGPROC) = NPARBLKWT(1,IAMGPROC) + 1
        NPARBLKWT(2,IAMGPROC) = NPARBLKWT(2,IAMGPROC) + MXSZTMP 
        NPARBLOCK(MXNUMB) = IAMGPROC - 1
      ELSE
        GOTO 500
      END IF
*     /\ IRUN !!!
      GOTO 200
*
500   CONTINUE
CSK      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.MASTER) THEN
CSK          WRITE(LUDUMPF,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
CSK     &    MYPROC:',MYPROC
          WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
     &    MYPROC:',MYPROC
          DO ISTI = 1, IMINNP
CSK            WRITE(LUDUMPF,'(3X,A,I17,A,I17,A,1X,I17)')
CSK     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
CSK     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
            WRITE(6,'(3X,A,I17,A,I17,A,1X,I17)')
     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
          END DO
        END IF
CSK      END IF
1001  CONTINUE
*      
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND_REL(NDIM,NBLOCKL,NPARBLOCK)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
*
      INTEGER*8 NPARBLKWT
      DIMENSION NBLOCKL(NDIM), NPARBLKWT(2,NMPROC)
*. Scratch
      INTEGER INONVAN,MAXLBL,LBL,IWEIGHT,LABEL,IBLOCKN,IMINNP
      INTEGER MXSZTMP,MXNUMB,ITEMPL,ITEMPD
      INTEGER*8 ITOTBLL 
*. ======
*. Output
*. ======
*. array NPARBLOCK that contains the final block distribution
      DIMENSION NPARBLOCK(NDIM)
*. 
      MXSIZE = 1000000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1000000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      INONVAN = 0
      ITOTBLL = 0
      IMINNP = 0
      MAXLBL = 0 
      LABEL = 0
      MIN2 = -2
      NZERO = 0
      NTEST = 0
      CALL ISETVC(NPARBLOCK,MIN2,NDIM)
      DO 100 II = 1, NDIM
        LBL = NBLOCKL(II)
        IF(LBL.GT.0) THEN
          INONVAN = INONVAN + 1
          NPARBLOCK(II) = -1
          ITOTBLL = ITOTBLL + LBL
        ENDIF
        IF(LBL.GT.MAXLBL) THEN
          MAXLBL = LBL
          LABEL = II
        END IF
 100  CONTINUE
      DO II = 1,2
        DO JJ = 1, NMPROC
          NPARBLKWT(II,JJ) = 0
        END DO
      END DO
*
*
      NTEST = 10
      IF(NTEST.GT.0) THEN
         WRITE(LUWRT,*) '  '
         WRITE(LUWRT,*) '  '
         WRITE(LUWRT,*) '  total number of processes to distribute on:',
     &                 NMPROC
         WRITE(LUWRT,*) '  total number of active blocks :', INONVAN
         WRITE(LUWRT,*) '  overall active block length   :', ITOTBLL 
         WRITE(LUWRT,*) '  number of largest active block:', LABEL
         WRITE(LUWRT,*) '  size of largest active block  :', MAXLBL
         WRITE(LUWRT,*) '  '
      END IF
*
      IF(INONVAN.LT.NMPROC) THEN
        IMINNP = MIN(INONVAN,NMPROC)
        IF(NTEST.GT.0) THEN
            WRITE(LUWRT,*) 'subroutine DISTBLKND speaking!'
            WRITE(LUWRT,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP - 1,'procs'
        END IF
        IF(MYPROC.GE.INONVAN) THEN
          IAM_NOT_INV = 0
          GOTO 1001
        END IF
      ELSE
        IMINNP = NMPROC
      END IF
      IAM_NOT_INV = 1
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
 200  CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.INONVAN) THEN
*
        MXSZTMP = 0
*
        DO 300 II = 1, NDIM
*
          ITEMPL = NBLOCKL(II)
          ITEMPD = NPARBLOCK(II)
          ITEMPN = II
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GE.MXSZTMP) THEN
                  MXSZTMP = ITEMPL
                  MXNUMB = II
                END IF
              END IF
            END IF
          END IF
*
300     CONTINUE
*
        DO 400 IPR = 1, IMINNP
*
          ITEMPW2 = NPARBLKWT(2,IPR)
* attention: IAMTPROC = 1 <--> MYPROC = 0 !
          IAMTPROC = IPR
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
            ITEMPW1 = ITEMPW2
            IAMGPROC = IAMTPROC
          END IF
*
 400    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*. calculation of block MXNUMB by proc IAMGPROC
        IF(NTEST.GT.0) THEN
          IF(MYPROC.EQ.MASTER) THEN
            WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,
     &      IAMGPROC-1
          END IF
        ENDIF
*. raising NPARBLKWT(1,proc) by 1
*. adding on NPARBLKWT(2,proc) the weight of the new block
        MXSIZE = MXSZTMP
        NPARBLKWT(1,IAMGPROC) = NPARBLKWT(1,IAMGPROC) + 1
        NPARBLKWT(2,IAMGPROC) = NPARBLKWT(2,IAMGPROC) + MXSZTMP 
        NPARBLOCK(MXNUMB) = IAMGPROC - 1
      ELSE
        GOTO 500
      END IF
*     /\ IRUN !!!
      GOTO 200
*
500   CONTINUE
      IF(NTEST.GE.0) THEN
          WRITE(LUWRT,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
     &    MYPROC:',MYPROC
          DO ISTI = 1, IMINNP
            WRITE(LUWRT,'(3X,A,I17,A,I17,A,1X,I17)')
     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
          END DO
        END IF
1001  CONTINUE
*      
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND_1(NDIM,ICWEIGHTF,NPARBLOCK,NBLKWT,NPARBLKWT,
     &                       NVAR,ICCTOS,IBLOCKL,IPROCLIST,ICWEIGHT,
     &                       IABSOLUTE_WEIGHT)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C     INPUT
      INTEGER ICOUNTABLK, NTEST,MWEIGHT
      INTEGER*8 ITOTBLCKL,MMPINFA
      DIMENSION ICWEIGHTF(NDIM),ICCTOS(NDIM,NDIM), IBLOCKL(NDIM)
      DIMENSION NPARBLOCK(NDIM),NBLKWT(NDIM),MMPINFA(NMPROC)
      DIMENSION MWEIGHT(NMPROC),IPROCLIST(NMPROC),ICWEIGHT(NDIM)
C     scratch
C     MXSIZE: actual maximum size for one of all blocks
C     MXNUMB: number of the current largest block
C     MMPINFA(*,1): total number of blocks for the n-th proc
C     MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2,ILEN
      INTEGER*8 LBPROC,MXSZTMP2,MXSIZE2,ITEMPL2,IABSOLUTE_WEIGHT
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
      LOGICAL SUCCESS
C     initializtion
      SUCCESS = .FALSE.
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 0
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 0
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO I = 1, NMPROC
        MMPINFA(I) = 0
        MWEIGHT(I) = 0
      END DO
C
C     count total number of active blocks
C     set need for non-zero block to -1
C
      DO ICBL = 1, NDIM
C
        ILEN = IBLOCKL(ICBL)
        IF(ILEN.NE.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
           ITOTBLCKL = ITOTBLCKL + ILEN
           NPARBLOCK(ICBL) = -1
        END IF
        IF(ILEN.GE.MXSIZE) MXSIZE = ILEN
C
      END DO
C
      ITEMPW1 = MXSIZE
C
      IF(NTEST.GT.0) THEN
       IF(MYPROC.EQ.0) THEN
         WRITE(6,*) '  total number of processes to distribute on:',
     &                 NMPROC
         WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
         WRITE(6,*) '  overall active block length:', ITOTBLCKL
       END IF
      END IF
*
      IMINNP = NMPROC
      NTEMPP = IMINNP
      IF(ICOUNTABLK.LT.NMPROC) THEN
        IMINNP = MIN(ICOUNTABLK,IMINNP)
        IF(NTEST.EQ.0) THEN
          IF(MYPROC.EQ.0) THEN
            WRITE(6,*) 'subroutine DISTBLKND_1 speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP,'procs'
          END IF
        END IF
        NTEMPP = IMINNP
        IF(MYPROC.GE.ICOUNTABLK) GOTO 1001
      END IF
C
C     maximum number of determinants on a given cpu
C
      LBPROC = (( NVAR -1 ) / (NMPROC)) + 1
C?      LBPROC = (( IABSOLUTE_WEIGHT  - 1 ) / (NTEMPP)) + 1
C?      IF(MYPROC.eq.MASTER) WRITE(6,*)'LBPROC is',LBPROC
C?      IF(MYPROC.eq.MASTER) WRITE(6,*)'NTEMPP is',NTEMPP
C
Csk?      IF(MYPROC.EQ.MASTER) THEN
Csk?        CALL IWRTMA(NPARBLOCK,1,NDIM,1,NDIM)
Csk?      END IF
C
C     find optimal (?) c-block distribution 
C
C
C     start with the 'master group' and with the master
      IGROUP = 1
      IPROC  = 1
      ICPROC = 0
      IASSIGNED_BLK = 0
C
 100  CONTINUE
C
      MXSZTMP = 0
      ICPROC = ICPROC + 1
      IF( ICPROC .gt. NTEMPP ) GOTO 500
      IF( IASSIGNED_BLK .ge. ICOUNTABLK ) GOTO 1001
C     find largest unassigned block:
C     MXNUMB  = block number
C     MXSZTMP = block 'weight'
C
      DO II = 1, NDIM
C
C        ITEMPL = NBLKWT(II)
        ITEMPL = IBLOCKL(II)
        ITEMPD = NPARBLOCK(II)
        ITEMPN = II
C?        IF(MYPROC.eq.MASTER) WRITE(6,*)'ITEMPL for II',ITEMPL,II
C
        IF(ITEMPL.GT.0)THEN
          IF(ITEMPD.EQ.-1) THEN
            IF(ITEMPL.LE.MXSIZE)THEN
              IF(ITEMPL.GT.MXSZTMP) THEN
C?                 IF(MYPROC.eq.MASTER) WRITE(6,*)'MXSZTMP,ITEMPL',
C?     & MXSZTMP,ITEMPL,MXSIZE
                 MXSZTMP = ITEMPL
                 MXNUMB = ITEMPN
              END IF
            END IF
          END IF
        END IF
      END DO
C
 150  CONTINUE
C
C?        IF(MYPROC.eq.MASTER) WRITE(6,*)'assigned in 150',MXNUMB,IPROC
        MMPINFA(IPROC) = MMPINFA(IPROC) + MXSZTMP
        MWEIGHT(IPROC) = MWEIGHT(IPROC) + 1
        NPARBLOCK(MXNUMB) = IPROC - 1
        MXSIZE = MXSZTMP
        IASSIGNED_BLK = IASSIGNED_BLK + 1
C
C     assign all blocks that are connected to MXNUMB to a given CPU
C     as long as MMPINFA <= LBPROC !
C
      SUCCESS = .FALSE.
C
 200  CONTINUE 
C
C     find largest block connected to MXNUMB
C
C     1. find all connections
      JJJ = 0
      DO IJ = 1, NDIM
        IF(ICCTOS(IJ,MXNUMB).NE.0) THEN
           IF( NPARBLOCK(IJ) .eq. -1 ) THEN
             JJJ = JJJ + 1
             ICWEIGHT(JJJ) = IJ
           END IF
        ENDIF
      END DO
C?      IF(MYPROC.eq.MASTER) WRITE(6,*)'JJJ is',JJJ
C
      IRUN = 0
 250  CONTINUE
      MXSIZE2 = MXSIZE
      MXSZTMP = 0
      IF( JJJ .eq. 0 ) GOTO 290
C
C?      IF(MYPROC.eq.MASTER) WRITE(6,*)'current IPROC is',IPROC
C?      IF(MYPROC.eq.MASTER) THEN
C?        WRITE(6,*)'current ICWEIGHT'
C?        CALL IWRTMA(ICWEIGHT,1,JJJ,1,JJJ)
C?      END IF
C
C?      IRUN = 0
 275  CONTINUE
      MXSZTMP = 0
C     2. select the largest block not yet assigned
      DO IASB = 1, JJJ
C
        JBLK = ICWEIGHT(IASB)
        ITEMPL2 = IBLOCKL(JBLK)
C?       ITEMPL2 = NBLKWT(JBLK)
        ITEMPD2 = NPARBLOCK(JBLK)
        ITEMPN2 = JBLK
C
        IF(ITEMPD2.EQ.-1) THEN
          IF(ITEMPL2.LE.MXSIZE2)THEN
            IF(ITEMPL2.GT.MXSZTMP) THEN
               MXSZTMP  = ITEMPL2
               MXNUMB2  = ITEMPN2
               MXSZTMP2 = ITEMPL2
            END IF
          END IF
        END IF
      END DO
C
C     assign
C
      IF( ( MMPINFA(IPROC) + MXSZTMP2 ) .le. LBPROC ) THEN
C?        IF(MYPROC.eq.MASTER) WRITE(6,*)'assigned in 250',MXNUMB2,IPROC
        MMPINFA(IPROC) = MMPINFA(IPROC) + MXSZTMP2
        MWEIGHT(IPROC) = MWEIGHT(IPROC) + 1
        NPARBLOCK(MXNUMB2) = IPROC - 1
        MXSIZE2 = MXSZTMP2
        IASSIGNED_BLK = IASSIGNED_BLK + 1
      ELSE
C
C       CPU in same group?
 280    CONTINUE
C
        SUCCESS = .FALSE.
C
        DO ITCPU = 1, NTEMPP
C
          IF(.NOT.SUCCESS) THEN
            JGROUP = IPROCLIST(ITCPU)
            IF( JGROUP .eq. IGROUP ) THEN
              IF( ITCPU .gt. IPROC ) THEN
                IPROC = ITCPU
                SUCCESS = .TRUE.
              END IF
            END IF
          END IF
C
        END DO
C
        IF( .NOT. SUCCESS ) THEN
          IGROUP = IGROUP + 1
C?          IF(MYPROC.eq.MASTER) WRITE(6,*)'IGROUP,NFLGRPS',IGROUP,NFLGRPS
          IF( IGROUP .gt. NFLGRPS ) GOTO 500
C         find the lowest CPU number in the new group
          IPROC = 1
          GOTO 280
        ELSE
C         fresh cpu
          GOTO 100
        END IF
C
      END IF
      IRUN = IRUN + 1
C?      IF(MYPROC.eq.MASTER) WRITE(6,*)'IRUN IS ',IRUN
      IF( IRUN .lt. JJJ ) THEN
C?      IF(MYPROC.eq.MASTER) WRITE(6,*)'I WILL GO TO 250'
        GOTO 250
      ELSE 
        GOTO 290
      END IF
C       CPU in same group?
 290    CONTINUE
C
        SUCCESS = .FALSE.
C
        DO ITCPU = 1, NTEMPP
C
          IF(.NOT.SUCCESS) THEN
            JGROUP = IPROCLIST(ITCPU)
            IF( JGROUP .eq. IGROUP ) THEN
              IF( ITCPU .gt. IPROC ) THEN
                IPROC = ITCPU
                SUCCESS = .TRUE.
              END IF
            END IF
          END IF
C
        END DO
C
        IF( .NOT. SUCCESS ) THEN
          IGROUP = IGROUP + 1
C?          IF(MYPROC.eq.MASTER) WRITE(6,*)'2. IGROUP,NFLGRPS',
C?     &                         IGROUP,NFLGRPS
          IF( IGROUP .gt. NFLGRPS ) GOTO 500
C         find the lowest CPU number in the new group
          IPROC = 1
          GOTO 290
        END IF
C       fresh cpu
        GOTO 100
C
 500  CONTINUE
C
C     check for 'lonely' blocks
C
      LONBLK = 0
C
      DO IBLK = 1, NDIM
        IF( NPARBLOCK(IBLK) .eq. -1 ) THEN 
          LONBLK = LONBLK + 1
          ICWEIGHT(LONBLK) = IBLK
        END IF
      END DO
C?      IF (MYPROC .eq. MASTER) WRITE(6,*) 'LONBLK is ',LONBLK
C
      IF( LONBLK .eq. 0 ) GOTO 1001
C
C     good old normal block distribution
C
      IRUN = 0
      MXSIZE = ITEMPW1
C
 600  CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.LONBLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, LONBLK
*
          IBLK   = ICWEIGHT(II) 
          ITEMPL = IBLOCKL(IBLK)
          ITEMPD = NPARBLOCK(IBLK)
          ITEMPN = IBLK
C
C?          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
                END IF
              END IF
            END IF
C?          END IF
*
3000    CONTINUE
*
        DO 4000 IPR = 1, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR) 
CSK         IF(MYPROC.EQ.0) 
CSK     & WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MYNEW_ID = 0 !
          IAMTPROC = IPR
CSK         IF(MYPROC.EQ.0) WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
CSK         IF(MYPROC.EQ.0) WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
C      now we should have found a proc and a block --> put both together !
C
C      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
C?       IF(NTEST.GT.0) THEN
C?         IF(MYPROC.EQ.0) 
C?     & WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
C?       ENDIF
C
C      raising MMPINFA(proc,1) by 1
C      adding on MMPINFA(proc,2) the weight of the new block
C
       IAMTPROC = IAMGPROC-1
       MXSIZE = MXSZTMP
       MWEIGHT(IAMGPROC) = MWEIGHT(IAMGPROC) + 1
       MMPINFA(IAMGPROC) = MMPINFA(IAMGPROC) + MXSZTMP
       NPARBLOCK(MXNUMB)     = IAMTPROC
C
      ELSE
        GOTO 1001
      END IF
C     /\ IRUN !!!
C
      GOTO 600
C
1001  CONTINUE
CSK      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.0) THEN
      WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF MYPROC:',
     &MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MWEIGHT(ISTI),' blocks with a
     & total length of',MMPINFA(ISTI)
      END DO
        END IF
CSK      END IF

C?      write(6,*) 'here i am',myproc
C?      IF(MYPROC.EQ.MASTER) THEN
C?        CALL IWRTMA(NPARBLOCK,1,NDIM,1,NDIM)
C?      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND_2(NDIM,ICWEIGHTF,NPARBLOCK,IBLOCKL)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C     -----
C     INPUT
C     -----
C     NDIM: total number of c-blocks
C     IBLOCKL(NDIM): total length of each block
C     ICWEIGHTF(NDIM): total 'weight factors' for each block
C
C     ------
C     OUTPUT
C     ------
*     NPARBLOCK(NDIM): list of blocks with the corresponding assigned node  
*
      INTEGER ICOUNTABLK, NTEST
      DIMENSION ICWEIGHTF(NDIM)
      DIMENSION NPARBLOCK(NDIM),IBLOCKL(NDIM)
C     -------
C     SCRATCH
C     -------
C     MXSIZE: actual maximum size for one of all blocks
C     MXNUMB: number of the current largest block
C     MMPINFA(*,1): total number of blocks for the n-th proc
C     MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
      INTEGER(KIND=MPI_OFFSET_KIND) ITLTND,MXSIZE,MXSZTMP
      INTEGER(KIND=MPI_OFFSET_KIND) ITEMPL, ITEMPW1, ITEMPW2
      INTEGER(KIND=MPI_OFFSET_KIND) ILEN, ITOTBLCKL 
      INTEGER(KIND=MPI_OFFSET_KIND) IIWEIGHTBLK, IILENGTH 
*
      NTEST = 0
      IIWEIGHTBLK = 0
      IILENGTH = 0
      ITLTND = 0
      ICOUNTABLK = 0
      MXSIZE = 0
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 0
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      ILEN = 0
*. end of initialization
      DO ICBL = 1, NDIM
C
        IIWEIGHTBLK = ICWEIGHTF(ICBL)
        IILENGTH = IBLOCKL(ICBL) 
C
        ILEN = IILENGTH * IIWEIGHTBLK 
C
        IF( ILEN .lt. 0 ) THEN 
          WRITE(6,*)'Attention, minus block detected',ICBL
          WRITE(6,*)'ILEN = ',ILEN
          WRITE(6,*)'IBLOCKL(ICBL) = ',IBLOCKL(ICBL)
          WRITE(6,*)'ICWEIGHTF(ICBL) = ',ICWEIGHTF(ICBL)
        END IF
        IF(ILEN.ne.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
           ITOTBLCKL = ITOTBLCKL + ILEN
           NPARBLOCK(ICBL) = -1
        END IF
        IF(ILEN.GE.MXSIZE) MXSIZE = ILEN
*
      END DO
*
CSK      WRITE(6,*)'MXSIZE of MYPROC',MXSIZE,MYPROC
      ITEMPW1 = MXSIZE
*
      IF(NTEST.GT.0) THEN
       IF(MYPROC.EQ.0) THEN
         WRITE(6,*) '  total number of processes to distribute on:',
     &                 NMPROC
         WRITE(6,*) '  total number of blocks:', NDIM
         CALL IWRTMA(NPARBLOCK,1,NDIM,1,NDIM)
         WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
         WRITE(6,*) '  overall active block length  :', ITOTBLCKL
         WRITE(6,*) '  Maximum block size           :', MXSIZE
       END IF
      END IF
*
      IMINNP = NMPROC
      NTEMPP = IMINNP
      IF(ICOUNTABLK.LT.NMPROC) THEN
        IMINNP = MIN(ICOUNTABLK,IMINNP)
        IF(NTEST.EQ.0) THEN
          IF(MYPROC.EQ.0) THEN
            WRITE(6,*) 'subroutine FIND_SBATCH_DISTR speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP,'procs'
          END IF
        END IF
        NTEMPP = IMINNP
        IF(MYPROC.GE.ICOUNTABLK) GOTO 101
      END IF
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
100   CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.ICOUNTABLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, NDIM
*
          
          IIWEIGHTBLK = ICWEIGHTF(II)
          IILENGTH = IBLOCKL(II)
C
          ITEMPL = IILENGTH * IIWEIGHTBLK
C
          ITEMPD = NPARBLOCK(II)
          ITEMPN = II

CSK          IF(MYPROC.EQ.0) 
CSK     & WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
CSK          IF(MYPROC.EQ.0) WRITE(6,*) 'MXSIZE:',MXSIZE
CSK          IF(MYPROC.EQ.0) WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
CSK                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
CSK       IF(MYPROC.EQ.0) WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 1, NTEMPP
*
C?          ITEMPW2 = MMPINFA(IPR) 
           ITEMPW2 = 0
           DO IBLK = 1, NDIM
            IF( NPARBLOCK(IBLK) .eq. IPR-1 ) THEN
              IIWEIGHTBLK = ICWEIGHTF(IBLK)
              IILENGTH = IBLOCKL(IBLK)
              ITEMPW2 = ITEMPW2 + ( IILENGTH * IIWEIGHTBLK )
            END IF
           END DO
          
CSK         IF(MYPROC.EQ.0) 
CSK     & WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MYNEW_ID = 0 !
          IAMTPROC = IPR
CSK         IF(MYPROC.EQ.0) WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
CSK         IF(MYPROC.EQ.0) WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
       IF(NTEST.GT.0) THEN
         IF(MYPROC.EQ.0) 
     & WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
*
*.     raising MMPINFA(proc,1) by 1
*.     adding on MMPINFA(proc,2) the weight of the new block
       IAMTPROC = IAMGPROC-1
       MXSIZE = MXSZTMP
C?       MWEIGHT(IAMGPROC) = MWEIGHT(IAMGPROC) + 1
C?       MMPINFA(IAMGPROC) = MMPINFA(IAMGPROC) + MXSZTMP
       NPARBLOCK(MXNUMB)     = IAMTPROC
*
      ELSE
        GOTO 101
      END IF
*     /\ IRUN !!!
*
      GOTO 100
*
101   CONTINUE
*
CSK      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.0) THEN
        DO II = 1, NDIM
         IF( NPARBLOCK(II) .eq. -1 ) THEN
          WRITE(6,*) 'AAAHHHH, block',II,' is not distributed'
          WRITE(6,*) 'Since I do not know how to proceed, I will stop!'
            Call Abend2('NE DISTBLKND_2: Block not distributed')
         END IF
        END DO
      WRITE(6,*)' '
      WRITE(6,'(20X,A)')'================================'
      WRITE(6,'(20X,A)')' Summation of even distribution '
      WRITE(6,'(20X,A)')'================================'
      WRITE(6,*)' '
      DO ISTI = 1, NTEMPP
        ITLTND = 0
        ITLBND = 0
        DO IBLK = 1, NDIM
          IF( NPARBLOCK(IBLK) .eq. ISTI -1 ) THEN
C
             IIWEIGHTBLK = ICWEIGHTF(IBLK)
             IILENGTH = IBLOCKL(IBLK)
C
             ITLTND = ITLTND + ( IILENGTH * IIWEIGHTBLK )
             ITLBND = ITLBND + 1
          END IF
        END DO
        WRITE(6,'(3X,A,1X,I4,1X,A,1X,I6,1X,A,1X,I18)')
     & 'CPU',ISTI-1,' computes',ITLBND,'blocks with a total weight of',
     & ITLTND
      END DO
*
*      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
CSK      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND_3(NDIM,NBLOCKL,NPARBLOCK,NPARBLKWT,
     &                       NBATCH,LBATCH,LEBATCH,I1BATCH)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/clunit.inc"
*
      DIMENSION NBLOCKL(NDIM), NPARBLKWT(2,NMPROC)
      DIMENSION LBATCH(*),LEBATCH(*),I1BATCH(*)
*. Scratch
      INTEGER INONVAN,MAXLBL,LBL,IWEIGHT,LABEL,IBLOCKN
      INTEGER MXSZTMP,MXNUMB,ITEMPL,ITEMPD
      INTEGER*8 ITOTBLL 
*. ======
*. Output
*. ======
*. array NPARBLOCK that contains the final block distribution
      DIMENSION NPARBLOCK(NDIM)
*. 
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1000000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      INONVAN = 0
      ITOTBLL = 0
      MAXLBL = 0 
      LABEL = 0
      MIN2 = -2
      NTEST = 0
      CALL ISETVC(NPARBLOCK,MIN2,NDIM)
      DO 100 II = 1, NDIM
        LBL = NBLOCKL(II)
        IF(LBL.GT.0) THEN
          NPARBLOCK(II) = -1
          ITOTBLL = ITOTBLL + LBL
        ENDIF
        IF(LBL.GT.MAXLBL) THEN
          MAXLBL = LBL
          LABEL = II
        END IF
 100  CONTINUE
*
*
CSK      NTEST = 10
      IF(NTEST.GT.0) THEN
       IF(MYPROC.EQ.MASTER) THEN
         WRITE(6,*) '  total number of processes to distribute on:',
     &                 NMPROC
         WRITE(6,*) '  overall active block length:  ', ITOTBLL 
         WRITE(6,*) '  largest active block:         ', LABEL
         WRITE(6,*) '  size of largest active block: ', MAXLBL
       END IF
      END IF
*
*. gather batch information for distribution of batches to cpu's (nodes)
*
         if(myproc.eq.master) then
         write(6,*)'NUMBER OF BATCHES:', NBATCH
         do II = 1, NBATCH
         write(6,*)'TOTAL LENGTH OF BATCH:', II,LEBATCH(II)
         write(6,*)'FIRST BLOCK OF THIS BATCH:', I1BATCH(II)
         write(6,*)'TOTAL NUMBER OF BLOCKS FOR THIS BATCH:',LBATCH(II)
         end do
         endif
CSK      LEBATCHTOT = DASUM(NBATCH,LEBATCH,1)
CSK      LEBATCHAVE = LEBATCHTOT / NMPROC
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
 200  CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.NBATCH) THEN
*
*       Find biggest batch not assigned to a CPU yet:
*       
        MXSZTMP = 0
*
        DO II = 1, NBATCH
*
*         Has this batch already been assigned to a CPU ?
          IBLK1 = I1BATCH(II)
          DO IBLK = IBLK1, IBLK1-1+LBATCH(II)
             ITEMPD = NPARBLOCK(IBLK)
             IF (ITEMPD .GE. 0) GO TO 300
*            .. Yes! by CPU no. ITEMPD
          END DO
*
          ITEMPL = LEBATCH(II)
          IF(ITEMPL.GE.MXSZTMP) THEN
            MXSZTMP = ITEMPL
            JBATCH = II
          END IF
*
  300   CONTINUE
        END DO
*
*       Find least loaded CPU:
*
        ITEMPW1 = NPARBLKWT(2,1)
        IAMGPROC = 1
*       Is anyone less loaded than the master (IPR = 1) ?
        DO 400 IPR = 2, NMPROC
*
          ITEMPW2 = NPARBLKWT(2,IPR)
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
            ITEMPW1 = ITEMPW2
            IAMGPROC = IPR
          END IF
*
 400    CONTINUE
*. now we should have found a proc and a batch --> put both together !
*
*. calculation of batch JBATCH by proc IAMGPROC
        IF(NTEST.GT.0) THEN
          IF(MYPROC.EQ.MASTER) THEN
            WRITE(6,*)'calculation of batch by proc',JBATCH,
     &      IAMGPROC-1
          END IF
        ENDIF
*. raising NPARBLKWT(1,proc) by 1
*. adding on NPARBLKWT(2,proc) the weight of the new block
        NPARBLKWT(1,IAMGPROC) = NPARBLKWT(1,IAMGPROC) + LBATCH(JBATCH)
        NPARBLKWT(2,IAMGPROC) = NPARBLKWT(2,IAMGPROC) + MXSZTMP 
        IBLK1 = I1BATCH(JBATCH)
        DO IBLK = IBLK1, IBLK1-1+LBATCH(JBATCH)
          IF (NBLOCKL(IBLK) .GT. 0) THEN
*           .. assign this block to the selected cpu
            NPARBLOCK(IBLK) = IAMGPROC - 1
          ELSE
*           .. insert code -2 for zero length block
            NPARBLOCK(IBLK) = -2
          END IF
        END DO
*
      ELSE
        GOTO 500
      END IF
*     /\ IRUN !!!
      GOTO 200
*
500   CONTINUE
CSK      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.MASTER) THEN
CSK          WRITE(LUDUMPF,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
CSK     &    MYPROC:',MYPROC
          WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
     &    MYPROC:',MYPROC
          DO ISTI = 1, NMPROC
CSK            WRITE(LUDUMPF,'(3X,A,I17,A,I17,A,1X,I17)')
CSK     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
CSK     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
            WRITE(6,'(3X,A,I17,A,I17,A,1X,I17)')
     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
          END DO
        END IF
CSK      END IF
1001  CONTINUE
*      
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DMTVDS_PAR(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
     &                     ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
*
* Multiply/divide elements of two vectors residing on disc
* Elements corresponding to absolute adresses in ISCAT
* are set to the elements of XSCAT
*
* For elements not in ISCAT the operation is thus :
*
* For INV.NE. 0 :
*
*    V3(I) = (V1(I)+FAC)-1 * V2(I)
*    LU3      LU1            LU2
*
* For INV.EQ. 0 :
*
*    V3(I) = (V1(I)+FAC) * V2(I)
*    LU3         LU1        LU2
*
* LBLK defines structure of files
*
* Overlap between input and output vector is also calculated
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION  VEC1(*),VEC2(*)
      DIMENSION  XSCAT(NSCAT),ISCAT(NSCAT)
*
      REAL*8 INPROD
*
      XINOUT = 0.0D0
*
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE( LU1,LBLK)
          CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
*
* Loop over blocks
*
        IOFF = 1
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          CALL IFRMDS(LBL2,1,-1,LU2)
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVSD_PAR : '
     &                     , LBL1,LBL2
          Call Abend2( ' DIFFERENT BLOCKSIZES IN DMTVSD_PAR ' )
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          LENGTH = LBL1
          CALL FRMDSC(VEC1,LENGTH,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,LENGTH,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
*             nothing needs to be done here
              CALL DIAVC2(VEC1,VEC2,VEC1,FAC,LENGTH)
            ELSE
*             nothing needs to be done here
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
*             nothing needs to be done here
              CALL DAXPY(LENGTH, FAC, VEC2, 1, VEC1, 1 )
            END IF
          END IF
*
          IF(NSCAT.NE.0) THEN
            IFIRST = IOFF
            ILAST = IOFF + LENGTH - 1
            DO 100 I = 1, NSCAT
              IF(IFIRST .LE. ISCAT(I) .AND. ISCAT(I) .LE. ILAST )
     &        VEC1(ISCAT(I)-IOFF+1) = XSCAT(I)
  100       CONTINUE
          END IF
*
          XINOUT = XINOUT + DDOT(LENGTH,VEC1,1,VEC2,1)
          CALL TODSC(VEC1,LENGTH,KBLK,LU3)
          IOFF = IOFF + LENGTH
        END IF
*
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DMTVDS_PAR2(VEC1,VEC2,LUDIA,LUIN,LUOUT,FAC,IREW,INV,
     &                       ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
*
* Multiply/divide elements of two vectors residing on disc
* Elements corresponding to absolute adresses in ISCAT
* are set to the elements of XSCAT
*
* For elements not in ISCAT the operation is thus :
*
* For INV.NE. 0 :
*
*    V3(I) = (V1(I)+FAC)-1 * V2(I)
*    LUOUT      LUDIA            LUIN
*
* For INV.EQ. 0 :
*
*    V3(I) = (V1(I)+FAC) * V2(I)
*    LUOUT         LUDIA        LUIN
*
* LBLK defines structure of files
*
* Overlap between input and output vector is also calculated
*
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/mxpdim.inc"
#include "../lucita/crun.inc"

      DIMENSION  VEC1(*),VEC2(*)
      DIMENSION  XSCAT(NSCAT),ISCAT(NSCAT)
*
      REAL*8 INPROD
*
      XINOUT = 0.0D0
C     set proper truncation factor
      THR_TRUNC  = TRUNC_FAC * RNORM_FAC
Csk      WRITE(6,*) 'TRUNCATION FACTOR:',THR_TRUNC
Chj 14-jun-07:   disable THR_ETRUNC
Chj   THR_ETRUNC = 1.0D-6 * THRES_E
      THR_ETRUNC = -1.0D0
*
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LUDIA
          REWIND LUIN
          REWIND LUOUT
        ELSE
          CALL REWINE( LUDIA,LBLK)
          CALL REWINE( LUIN,LBLK)
          CALL REWINE( LUOUT,LBLK)
         END IF
      END IF
*
* Loop over blocks
*
        IOFF = 1
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LUDIA)  LBL1
          READ(LUIN)   LBL2
          WRITE(LUOUT) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LUDIA)
          CALL IFRMDS(LBL2,1,-1,LUIN)
          CALL ITODS (LBL1,1,-1,LUOUT)
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVSD_PAR : '
     &                     , LBL1,LBL2
          Call Abend2( ' DIFFERENT BLOCKSIZES IN DMTVSD_PAR ' )
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          LENGTH = LBL1
          CALL FRMDSC(VEC1,LENGTH,KBLK,LUDIA,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,LENGTH,KBLK,LUIN,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
CSK            IOFFMAX = idamax(LENGTH, VEC2, 1 )
CSK            IF( ABS(VEC2(IOFFMAX)) .LT. THR_TRUNC ) THEN
CSK              CALL DZERO(VEC2,LENGTH)
CSK            ELSE
CSK              DO II = 1, LENGTH
CSK                IF( ABS(VEC2(II)) .LT. THR_TRUNC ) 
CSK     &              VEC2(II) = 0.0D0
CSK              END DO
CSK            END IF
            IF(INV .NE. 0 ) THEN
*             nothing needs to be done here
              
              CALL DIAVC2_TRUNC(VEC1,VEC2,VEC1,FAC,LENGTH,THR_TRUNC,
     &                          THR_ETRUNC)
            ELSE
*             nothing needs to be done here
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
*             nothing needs to be done here
              CALL DAXPY(LENGTH, FAC, VEC2, 1, VEC1, 1 )
            END IF
          END IF
*
          IF(NSCAT.NE.0) THEN
            IFIRST = IOFF
            ILAST = IOFF + LENGTH - 1
            DO 100 I = 1, NSCAT
              IF(IFIRST .LE. ISCAT(I) .AND. ISCAT(I) .LE. ILAST )
     &        VEC1(ISCAT(I)-IOFF+1) = XSCAT(I)
  100       CONTINUE
          END IF
*
          XINOUT = XINOUT + DDOT(LENGTH,VEC1,1,VEC2,1)
          CALL TODSC(VEC1,LENGTH,KBLK,LUOUT)
          IOFF = IOFF + LENGTH
        END IF
*
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_CCTOS(SCALFAC_LOCAL,IBLOCKD,ICCTOS,NBLOCK)
*
*. Find all c-blocks connecting to a given sigma-block for each cpu
*  using the connection matrix ICCTOS.
*. Each cpu stores the information in SCALFAC_LOCAL.
*. Connection is marked by 1.0D0
*
      IMPLICIT REAL*8           ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
*     ==========
*       INPUT
*     ==========
      DIMENSION ICCTOS(NBLOCK,NBLOCK),IBLOCKD(NBLOCK)
*
*     ===========
*       OUTPUT
*     ===========
*
      DIMENSION SCALFAC_LOCAL(NBLOCK)
*
      DO IBLK = 1, NBLOCK
*
        IF(IBLOCKD(IBLK).EQ.MYPROC) THEN
*
          DO JBLK = 1,NBLOCK
*
            IF(ICCTOS(JBLK,IBLK).NE.0) SCALFAC_LOCAL(JBLK) = 1.0D0
*
          END DO
*
        END IF
*
      END DO
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_CCTOS_REL(ISCALFAC,IBLOCKD,ICCTOS,IBLOCKL,
     &                          I_RUN_COMPLEX,NBLOCK)
*
*. Find all c-blocks connecting to a given sigma-block for each cpu
*  using the connection matrix ICCTOS.
*. Each cpu stores the information in ISCALFAC.
*. Connection is marked by 1
*
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
*     ==========
*       INPUT
*     ==========
*                      SBLOCKS    , CBLOCKS
      DIMENSION ICCTOS(I_NZERO_LEN,I_NZERO_LEN), IBLOCKD( NBLOCK )
      DIMENSION IBLOCKL( NBLOCK )
*
*     ===========
*       OUTPUT
*     ===========
*
*     NOTE: IF NZ == 2 (complex quaternion algebra needed):
*     DIMENSION ISCALFAC(2*NBLOCK) else ISCALFAC(NBLOCK)
      DIMENSION ISCALFAC(*)
      IONE          = 1
      I_AM_ACTIVE_S = 0
      I_AM_ACTIVE_C = 0
      IPUT_BLK      = 0
*
      DO II_RC = 1, I_RUN_COMPLEX
*
        I_AM_ACTIVE_S = 0
        I_AM_ACTIVE_C = 0
*
        DO IBLK = 1, NBLOCK
*
          IF( IBLOCKL( IBLK) .eq. 0 ) GOTO 100
*
          I_AM_ACTIVE_S = I_AM_ACTIVE_S + IONE
*
          IF( IBLOCKD( IBLK ) .eq. MYPROC )THEN
*
            I_AM_ACTIVE_C = 0
*
            DO JBLK = 1,NBLOCK
*
              IF( IBLOCKL( JBLK ) .eq. 0 ) GOTO 200 
*
              I_AM_ACTIVE_C = I_AM_ACTIVE_C + IONE
*
              IPUT_BLK = JBLK + ( NBLOCK * ( II_RC - 1) )
*
              IF( ICCTOS( I_AM_ACTIVE_C, I_AM_ACTIVE_S ) .ne. 0 ) 
     &            ISCALFAC( IPUT_BLK ) = IONE
*
 200          CONTINUE
*
            ENDDO
*
          END IF
 100      CONTINUE
*
        ENDDO
      ENDDO
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_GROUP_OF_PROCS(MYOWNNAME,MYNAMELENGTH,
     &                               IPROCLIST,NAMELENGTHLIST)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C
      DIMENSION IPROCLIST(NMPROC),NAMELENGTHLIST(NMPROC)
      CHARACTER*72 PNAMELIST(NMPROC),MYOWNNAME,MYSCRNAME
      CHARACTER*72 OLD_NAME
      INTEGER INUMBER, IFINI, NTEST,ICOUNT_GROUPS
C
C     NTEST =  0
      NTEST = 10
      INUMBER = 0
      IFINI = 0
      ICOUNT_GROUPS = 0
C
C     start with an all_gather of all name length
C
      CALL MPI_ALLGATHER(MYNAMELENGTH,1,my_MPI_INTEGER,
     &                   NAMELENGTHLIST,1,my_MPI_INTEGER,
     &                   MPI_COMM_WORLD,IERR)
C
C     ... all name length are stored now!
C     ... save own name in scratch array!
C
      PNAMELIST(MYPROC+1) = MYOWNNAME
C
C     ... we start to gather the names ... 
      DO IPROC = 1, NMPROC
C
         MYSCRNAME = MYOWNNAME
         CALL MPI_BCAST(MYSCRNAME,NAMELENGTHLIST(IPROC),MPI_CHARACTER,
     &                 IPROC-1,MPI_COMM_WORLD,IERR)
         IF( MYPROC .NE. IPROC -1 ) THEN
           PNAMELIST(IPROC) = MYSCRNAME(1:NAMELENGTHLIST(IPROC))
         END IF
C
      END DO
C
C     test writing if needed
C
      IF( MYPROC .EQ. MASTER .AND. NTEST .GE. 10 ) THEN
       DO IPROC = 1, NMPROC
        MYSCRNAME = PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
        WRITE(6,'(1X,A,1X,I4,1X,A7,1X,A)')'Processor ',IPROC -1,
     & ' alias ',MYSCRNAME(1:NAMELENGTHLIST(IPROC))
C    & MYSCRNAME
       END DO
      END IF
C
C     find all processors on the same deck and reorder (if necessary) 
C     to get the processors as close as possible, starting with the
C     master 
C     reordering is currently NOT done
C
      INUMBER = MASTER + 1
 100  CONTINUE
      IFINI = 0
      ICOUNT_GROUPS = ICOUNT_GROUPS + 1
      OLD_NAME = PNAMELIST( INUMBER)(1:NAMELENGTHLIST( INUMBER))
C
      DO IPROC = 1, NMPROC
C
        MYSCRNAME(1:NAMELENGTHLIST(IPROC)) = 
     &  PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
C
        IF(MYSCRNAME(1:NAMELENGTHLIST(IPROC)) .eq. 
     &     OLD_NAME(1:NAMELENGTHLIST(INUMBER)))
     &     IPROCLIST(IPROC) = INUMBER
C
      END DO
C
C     if finished: IFINI = 1, else 0
C     search for the next lowest cpu building the 'local group master'
C
      IFINI = 1
      DO IPRC = 1, NMPROC
        IF(IFINI .ne. 0 .and. IPROCLIST(IPRC) .eq. -1) THEN
          IFINI = 0
          INUMBER = IPRC
        END IF
      END DO
      IF(IFINI .eq. 0 ) GOTO 100
C
C     write the count of 'groups' formed by processes
C
      IF( MYPROC .eq. MASTER .and. NTEST .ge. 10) THEN
        WRITE(6,'(1X,A,I4,1X,A)') 'There are ',ICOUNT_GROUPS,
     & 'different groups of processes'
        CALL IWRTMA(IPROCLIST,1,NMPROC,1,NMPROC)
      END IF
C
C     Transfer information to common block
C
      NFLGRPS = ICOUNT_GROUPS
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_GROUP_OF_PROCS_REL(MYOWNNAME,MYNAMELENGTH,
     &                                   IPROCLIST,NAMELENGTHLIST)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C
      DIMENSION IPROCLIST(NMPROC),NAMELENGTHLIST(NMPROC)
      CHARACTER*72 PNAMELIST(NMPROC),MYOWNNAME,MYSCRNAME
      CHARACTER*72 OLD_NAME
      INTEGER INUMBER, IFINI, NTEST,ICOUNT_GROUPS
C
C     NTEST =  0
      NTEST = 10
      INUMBER = 0
      IFINI = 0
      ICOUNT_GROUPS = 0
C
C     start with an all_gather of all name length
C
      CALL MPI_ALLGATHER(MYNAMELENGTH,1,my_MPI_INTEGER,
     &                   NAMELENGTHLIST,1,my_MPI_INTEGER,
     &                   MPI_COMM_WORLD,IERR)
C
C     ... all name length are stored now!
C     ... save own name in scratch array!
C
      PNAMELIST(MYPROC+1) = MYOWNNAME
C
C     ... we start to gather the names ... 
      DO IPROC = 1, NMPROC
C
         MYSCRNAME = MYOWNNAME
         CALL MPI_BCAST(MYSCRNAME,NAMELENGTHLIST(IPROC),MPI_CHARACTER,
     &                 IPROC-1,MPI_COMM_WORLD,IERR)
         IF( MYPROC .NE. IPROC -1 ) THEN
           PNAMELIST(IPROC) = MYSCRNAME(1:NAMELENGTHLIST(IPROC))
         END IF
C
      END DO
C
C     test writing if needed
C
      IF( NTEST .GE. 100 ) THEN
       DO IPROC = 1, NMPROC
        MYSCRNAME = PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
        WRITE(LUWRT,'(2X,A,1X,I4,1X,A7,1X,A)')'Processor ',IPROC -1,
     & ' alias ',MYSCRNAME(1:NAMELENGTHLIST(IPROC))
C    & MYSCRNAME
       END DO
      END IF
C
C     find all processors on the same deck and reorder (if necessary) 
C     to get the processors as close as possible, starting with the
C     master 
C
      INUMBER = MASTER + 1
 100  CONTINUE
      IFINI = 0
      ICOUNT_GROUPS = ICOUNT_GROUPS + 1
      OLD_NAME = PNAMELIST( INUMBER)(1:NAMELENGTHLIST( INUMBER))
C
      DO IPROC = 1, NMPROC
C
        MYSCRNAME(1:NAMELENGTHLIST(IPROC)) = 
     &  PNAMELIST(IPROC)(1:NAMELENGTHLIST(IPROC))
C
        IF(MYSCRNAME(1:NAMELENGTHLIST(IPROC)) .eq. 
     &     OLD_NAME(1:NAMELENGTHLIST(INUMBER)))
     &     IPROCLIST(IPROC) = INUMBER
C
      END DO
C
C     if finished: IFINI = 1, else 0
C     search for the next lowest cpu building the 'local group master'
C
      IFINI = 1
      DO IPRC = 1, NMPROC
        IF(IFINI .ne. 0 .and. IPROCLIST(IPRC) .eq. -1) THEN
          IFINI = 0
          INUMBER = IPRC
        END IF
      END DO
      IF(IFINI .eq. 0 ) GOTO 100
C
C     write the count of 'groups' formed by processes
C
      IF( NTEST .ge. 10) THEN
        WRITE(LUWRT,'(2X,I4,1X,A)') ICOUNT_GROUPS,'group(s) of 
     & processors found! '
        IF( NTEST .ge. 50) THEN
          WRITE(LUWRT,'(2X,A)') 'Complete listing:'
          WRITE(LUWRT,*) '  '
          CALL IWRTMAMN(IPROCLIST,1,NMPROC,1,NMPROC,LUWRT)
        END IF
        WRITE(LUWRT,*) '  '
      END IF
C
C     Transfer information to common block
C
      NFLGRPS_REL = ICOUNT_GROUPS
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
*
      SUBROUTINE GATVCD_PAR(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,IPRT)
*
* Gather elements from a file LU
*
* XGAT(I) = Vector(IGAT(I))
*
* Jeppe Olsen, September 1993
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
*. Input
      INTEGER IGAT(NGAT)
*. Output
      DIMENSION XGAT(NGAT)
*. Scratch
      DIMENSION SEGMNT(*)
*
      CALL REWINE(LU,-1)
*
      IBASE = 1
      IBLOCK = 0
*
*. Loop over blocks of file
*
 1000 CONTINUE
        IBLOCK = IBLOCK + 1
        CALL NEXREC(LU,LBLK,SEGMNT,IEND,LENGTH)
        IF(IPRT.GE.10)
     &  WRITE(6,*) LENGTH, ' elements in block ',IBLOCK
        IF(IEND.EQ.0) THEN
          IFIRST = IBASE
          ILAST = IBASE + LENGTH - 1
          DO 100 I = 1, NGAT
            IF(IFIRST .LE. IGAT(I) .AND. IGAT(I) .LE. ILAST )
     &      XGAT(I) = SEGMNT(IGAT(I)-IFIRST+1)
C?          IF(IFIRST .LE. IGAT(I) .AND. IGAT(I) .LE. ILAST )
C?   &      write(6,*) ' Catch I IGAT(I) XGAT(I) ',
C?   &                         I,IGAT(I),XGAT(I)
  100     CONTINUE
          IBASE = IBASE + LENGTH
      IF(LBLK.LT.0) GOTO 1000
        END IF
*
      NTEST = 0
      NTEST = MAX(IPRT,NTEST)
      IF(NTEST.GE.5) THEN
       WRITE(6,*) ' Gathered vector from GATVCD_PAR'
       CALL WRTMAT(XGAT,1,NGAT,1,NGAT)
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GROUP_CONSTRUCTOR_REL(IGROUPLIST,IPROCLIST)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C
      DIMENSION IPROCLIST(NMPROC),IGROUPLIST(NMPROC)
      INTEGER IKEY, ICOLOR, JKEY, JCOLOR
C
C     intra-node communicator MYNEW_COMM
C
      IKEY = MYPROC
      ICOLOR = IPROCLIST(MYPROC+1)
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,ICOLOR,IKEY,MYNEW_COMM,IERR)
C
C     collect useful information about each group,
C     store on common block
C
      NEWCOMM_PROC = 0
      MYNEW_ID = 0
      CALL MPI_COMM_SIZE(MYNEW_COMM,NEWCOMM_PROC,IERR)
      CALL MPI_COMM_RANK(MYNEW_COMM,MYNEW_ID,IERR)
C
C     inter-node communicator ICOMM
C
      IF( MYNEW_ID .eq. 0 ) THEN
        JKEY = IPROCLIST(MYPROC+1)
        JCOLOR = 2
      ELSE
        JKEY = IPROCLIST(MYPROC+1)
        JCOLOR = 3
      END IF
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,JCOLOR,JKEY,ICOMM,IERR)
C
C     collect again ...     
C
      ICOMM_ID = 0
      ICOMM_SIZE = 0
      CALL MPI_COMM_RANK(ICOMM,ICOMM_ID,IERR)
      CALL MPI_COMM_SIZE(ICOMM,ICOMM_SIZE,IERR)
C
C     set-up and store group information 
C
CSK      IF( NEWCOMM_PROC .ne. NMPROC ) THEN
        CALL SET_GROUP_TABLE_REL(IGROUPLIST,IPROCLIST,ICOLOR)
CSK      END IF
C
C     store personal group number
C
      MY_GROUPN = ICOLOR
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GROUP_DESTRUCTOR_REL(JCOMM1,JCOMM2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C
      CALL MPI_COMM_FREE(JCOMM1,IERR)
      CALL MPI_COMM_FREE(JCOMM2,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE I_DO_RCCTOS(IBLOCKL,ICCTOS,NBLOCKS,ICOMPL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IBLOCKL(*), ICCTOS(*)
C
      IONE = 1
      DO ILOOP = 1, ICOMPL
        DO IBLK = 1, NBLOCKS
          IF( IBLOCKL(IBLK) .gt. 0 ) THEN
            ICCTOS( IBLK + ( NBLOCKS * ( ILOOP -1 ) ) ) = IONE
          END IF
        END DO
      END DO
C
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE ISET_ARRAY_ACT_BLK(IAM_BLK_ACT,IBLOCKL,NDIM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IBLOCKL(*), IAM_BLK_ACT(*)
C
        ICOUNT = 0
C
        DO IBLK = 1, NDIM
          IF( IBLOCKL(IBLK) .gt. 0 ) THEN
            ICOUNT = ICOUNT + 1
            IAM_BLK_ACT( IBLK ) = ICOUNT
          END IF
        END DO
C
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INFO_PRINT_BATCH_REL(LBATCH, LEBATCH,
     &                                I1BATCH,IBATCH,JBATCH,LUPRINT)
      IMPLICIT REAL*8(A-H,O-Z)
C     Input
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
C
C
C     
      DO ISBATCH = 1, JBATCH
        WRITE(LUPRINT,*)
        WRITE(LUPRINT,*) ' Info on batch ', ISBATCH
        WRITE(LUPRINT,*) ' *********************** '
        WRITE(LUPRINT,*)
        WRITE(LUPRINT,*) '      Number of blocks included ',
     &                          LBATCH(ISBATCH)
        WRITE(LUPRINT,*) '      Length of batch           ',
     &                          LEBATCH(ISBATCH)
        WRITE(LUPRINT,*) '      TTSS and offsets and lenghts of
     &                          each block '
        DO IBLOCK = I1BATCH(ISBATCH),I1BATCH(ISBATCH)+ LBATCH(ISBATCH)-1
          WRITE(LUPRINT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
        END DO
      END DO
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRODB_PAR(VEC1,VEC2,NBLK,LBLK,I0BLK)
C
C     inner products between blocked vectors with check of zero blocks.
C     Zero blocks are flagged by a 0 entry in I0BLK.
C
C     adaption for parallel version from original routine written by 
C                        J. Olsen        -   May  1997
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C     Input
      DIMENSION VEC1(*),VEC2(*), LBLK(NBLK),I0BLK(NBLK)
*
      X = 0.0D0
      IOFF = 1
      DO IBLK = 1, NBLK
        LLBLK = LBLK(IBLK)
C?      WRITE(6,*) ' INPRODB IBLK LLBLK IOFF ',IBLK,LLBLK,IOFF
        IF(I0BLK(IBLK).EQ.0.AND.LLBLK.GT.0) THEN
          X = X + DDOT(VEC1(IOFF),VEC2(IOFF),LLBLK)
C?        WRITE(6,*) ' Vec1 and Vec2 blocks '
C?        CALL WRTMAT(VEC1(IOFF),1,LLBLK,1,LLBLK)
C?        CALL WRTMAT(VEC2(IOFF),1,LLBLK,1,LLBLK)
C?        write(6,*) ' Updated x', X
        END IF
        IF(LLBLK.GT.0) IOFF = IOFF + LLBLK
      END DO
C
      INPRODB_PAR = X
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      CHARACTER*10 FUNCTION LFNAM_LUCIPAR(FILNAM)
C
C     Add the node identification to the file name for shared disk
C     systems
C     original written by Luuk Visscher, February 1998
C     adaption for parallel test of LUCIAREL by Stefan Knecht, June 2007
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "dcbgen.h"
#include "infpar.h"
      CHARACTER*5  FILNAM
      CHARACTER*10 LFILNAM
      NODE = 0
 1000 FORMAT (A5,I5)
C
      WRITE (LFILNAM,1000) FILNAM,NODE
      DO I = 6, 10
         IF (LFILNAM(I:I).EQ.' ') LFILNAM(I:I) = 'X'
      ENDDO
C
      LFNAM_LUCIPAR = LFILNAM
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD_PAR(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,
     &                      NVEC,IREW,LBLK)
C
C ADD VECTORS ON FILE LUIN TIMES FACTOR AND STORE ON LUOUT
C
C LUOUT AND LUSCR ARE INITIALLY REWOUND
C NVEC = NROOT resp. NVEC - NROOT
C
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(NVEC)
C
      IF( MOD(NVEC,2) .EQ. 0 ) THEN
        LLUOUT = LUSCR
        LLUSCR = LUOUT
      ELSE
        LLUOUT = LUOUT
        LLUSCR = LUSCR
      END IF
C
      IF(IREW .NE. 0 ) CALL REWINE(LUIN,LBLK)
C
      DO 100 IVEC = 1, NVEC
CSK        WRITE(6,*) 'FAC(IVEC) 1 = ',FAC(IVEC)
        CALL REWINE(LLUSCR,LBLK)
        CALL REWINE(LLUOUT,LBLK)
        IF( IVEC .EQ. 1 ) THEN
          CALL SCLVCD(LUIN,LLUOUT,FAC(IVEC),VEC1,0,LBLK)
        ELSE
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),1.0D0,LUIN,LLUSCR,LLUOUT,
     &                0,LBLK)
        END IF
C
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
  100 CONTINUE
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD_BATCH_PAR(VEC1,VEC2,FACIN1,FACIN2,NVECIN1,
     &                            NVECIN2,
     &                            NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                            MY_VEC_IOFF,MY_IOFF_LUIN1,
     &                            MY_IOFF_LUIN2,
     &                            MY_IOFF_LUOUT,
     &                            NBLK_ACT,LUIN1,LUIN2,
     &                            LUOUT,JOFF)
C
C     add vectors on files LUIN1 and LUIN2 times a factor and store on LUOUT
C
C     VEC2 = ( VEC1 x FACIN1 ) + ( VEC1 x FACIN2 ) + VEC2
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION FACIN1(NVECIN1),FACIN2(*)
      DIMENSION LBATCH(*),LEBATCH(*),I1BATCH(*),IBATCH(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
C     initialize scratch offsets
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVECIN1
C
C         set new offset and zero read-in vector
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH + 
     &                    ( IVEC - 1 ) * MY_VEC_IOFF
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN1',IOFFSET_IN_LUIN1
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN1,IOFFSET_IN_LUIN1)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IVEC)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            WRITE(6,*) 'FACIN1(IVEC) = ',FACIN1(IVEC)
CSK          END IF
          CALL DAXPY(LEBATCH(ISBATCH),FACIN1(IVEC),VEC1,1,VEC2,1)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 after DAXPY'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
  100   CONTINUE
C       next vectors on LUIN2 ...
        DO 200 IVEC = 1, NVECIN2
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                      ( IVEC - 1 ) * MY_VEC_IOFF
C
C         read in batch ISBATCH from LUIN2 to VEC1
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN2,IOFFSET_IN_LUIN2)
C
C         VEC2 == VEC2 + VEC1 * FACIN2(IVEC)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY in 200'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            WRITE(6,*) 'FACIN2(IVEC) = ',FACIN2(IVEC)
CSK          END IF
          CALL DAXPY(LEBATCH(ISBATCH),FACIN2(IVEC),VEC1,1,VEC2,1)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 after DAXPY in 200'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C
C
  200   CONTINUE
CSK        IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'VEC2 after 200'
CSK           CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK        END IF
C
C       write VEC2 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + ( MY_VEC_IOFF * JOFF ) + 
     &                      IOFFSET_SCRATCH
        CALL WTVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUOUT,IOFFSET_LUOUT,0)
CCC       IF(MYPROC.EQ.MASTER+1) THEN
CCC         WRITE(6,*) 'This is what i have written = reset C/svector'
CCC         CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CCC       END IF
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD_BATCH_PAR_USE_IA(VEC1,VEC2,FACIN1,FACIN2,
     &                                   NVECIN1,NVECIN2,NBATCH,LBATCH,
     &                                   LEBATCH,I1BATCH,IBATCH,
     &                                   MY_VEC_IOFF,MY_IOFF_LUIN1,
     &                                   MY_IOFF_LUIN2,
     &                                   MY_IOFF_LUOUT,
     &                                   NBLK_ACT,
     &                                   LUIN1,LUIN2,LUOUT,JOFF,
     &                                   IVCOFF_IN1,IVCOFF_IN2,
     &                                   IVCOFF_OUT)
C
C     add vectors on files LUIN1 and LUIN2 times a factor and store on LUOUT
C
C     VEC2 = ( VEC1 x FACIN1 ) + ( VEC1 x FACIN2 ) + VEC2
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION FACIN1(NVECIN1),FACIN2(*)
      DIMENSION LBATCH(*),LEBATCH(*),I1BATCH(*),IBATCH(8,*)
      DIMENSION IVCOFF_IN1(*),IVCOFF_IN2(*),IVCOFF_OUT(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVECIN1
C
C         set new offset and zero read-in vector
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH + 
     &                    ( IVEC - 1 ) * MY_VEC_IOFF
            IOFFSET_IN_INT_LUIN1 = NUM_BLK + ( IVEC - 1 ) * NBLK_ACT
     &                           + 1
C?            write(6,*) 'my integer offset',IOFFSET_IN_INT_LUIN1,myproc
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          LUIN1,IOFFSET_IN_LUIN1,
     &                          IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
C?            write(6,*) 'survived!',NVECIN1,NBATCH,myproc
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IVEC)
C
          CALL SCLVCSUM(FACIN1(IVEC),VEC1,VEC2,LBATCH(ISBATCH),
     &                  IBATCH(1,I1BATCH(ISBATCH)),
     &                  IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
C
  100   CONTINUE
C       next vectors on LUIN2 ...
        DO 200 IVEC = 1, NVECIN2
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                      ( IVEC - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN2 = NUM_BLK + ( IVEC - 1 ) * NBLK_ACT
     &                         + 1
C
C         read in batch ISBATCH from LUIN2 to VEC1
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          LUIN2,IOFFSET_IN_LUIN2,
     &                          IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
C         VEC2 == VEC2 + VEC1 * FACIN2(IVEC)
          CALL SCLVCSUM(FACIN2(IVEC),VEC1,VEC2,LBATCH(ISBATCH),
     &                  IBATCH(1,I1BATCH(ISBATCH)),
     &                  IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
  200   CONTINUE
C
C       write VEC2 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + ( MY_VEC_IOFF * JOFF ) + 
     &                      IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = NUM_BLK + ( NBLK_ACT * JOFF )
     &                    + 1
        CALL WTVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUOUT,IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       IVCOFF_OUT,0)
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK = NUM_BLK + LBATCH(ISBATCH)
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE RDVEC_BATCH_DRV(SCR,NBATCH_BLK,NBATCH_LEN,NBATCH_INF,
     &                           MY_LUIN,IOFFSET)
C
C     read in a batch from a MPI-file --> (hint LUIN) to SCR
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/clunit.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
CSK      IF( MYPROC .eq. MASTER ) THEN
CSK        WRITE(6,*) 'MY_LUIN, LUC, LUHC, LUSC1, LUSC2, LUSC3, LUSC35',
CSK     & MY_LUOUT, LUC, LUHC, LUSC1, LUSC2, LUSC3, LUSC35
CSK      END IF
C  
        IF( MY_LUIN .eq. LUC ) THEN
          CALL MPI_FILE_READ_AT(ILU1,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUIN .eq. LUHC ) THEN
          CALL MPI_FILE_READ_AT(ILU2,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUIN .eq. LUSC1 ) THEN
          CALL MPI_FILE_READ_AT(ILU3,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUIN .eq. LUSC2 ) THEN
          CALL MPI_FILE_READ_AT(ILU4,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUIN .eq. LUSC3 ) THEN
          CALL MPI_FILE_READ_AT(ILU5,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUIN .eq. LUSC35 ) THEN
          CALL MPI_FILE_READ_AT(ILU7,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        END IF
C
        IOFFSET = IOFFSET + NBATCH_LEN
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE WTVEC_BATCH_DRV(SCR,NBATCH_BLK,NBATCH_LEN,NBATCH_INF,
     &                           MY_LUOUT,IOFFSET,NO_CHECK)
C
C     write batch from SCR to MPI-file --> (hint LUOUT)
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/clunit.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF, LENGTH
      JOFF = 0
C
CSK      IF( MYPROC .eq. MASTER ) THEN
CSK        WRITE(6,*) 'MY_LUOUT, LUC, LUHC, LUSC1, LUSC2, LUSC3, LUSC35',
CSK     & MY_LUOUT, LUC, LUHC, LUSC1, LUSC2, LUSC3, LUSC35
CSK      END IF
C
        IF( MY_LUOUT .eq. LUC ) THEN
          CALL MPI_FILE_WRITE_AT(ILU1,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUOUT .eq. LUHC ) THEN
          CALL MPI_FILE_WRITE_AT(ILU2,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUOUT .eq. LUSC1 ) THEN
          CALL MPI_FILE_WRITE_AT(ILU3,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUOUT .eq. LUSC2 ) THEN
          CALL MPI_FILE_WRITE_AT(ILU4,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUOUT .eq. LUSC3 ) THEN
          CALL MPI_FILE_WRITE_AT(ILU5,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        ELSE IF( MY_LUOUT .eq. LUSC35 ) THEN
          CALL MPI_FILE_WRITE_AT(ILU7,IOFFSET,SCR,
     &                  NBATCH_LEN,MPI_REAL8,ISTAT,IERR)
        END IF
C
        IOFFSET = IOFFSET + NBATCH_LEN
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_BATCH_PAR(SCR,NBATCH,LBATCH,LEBATCH,
     &                            I1BATCH,IBATCH,MY_VEC_IOFF,
     &                            MY_IOFF_LUIN,MY_IOFF_LUOUT,
     &                            NBLK_ACT,LUIN,LUOUT,JOFF)
C
C     copy vector from file LUIN to file LUOUT --> (hint LUIN, LUOUT)
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SCR(*),LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
C
C     initialize scratch offsets
      IOFFSET_SCRATCH = 0
      IOFFSET_IN  = 0
      IOFFSET_OUT  = 0
C
      DO ISBATCH = 1, NBATCH 
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(SCR,LEBATCH(ISBATCH))
C
        IOFFSET_IN = MY_IOFF_LUIN + ( MY_VEC_IOFF * JOFF ) +
     &               IOFFSET_SCRATCH
C
CCC     IF(mpyproc.eq.master) then 
CCC       write(6,*) 'offset for luin',IOFFSET_INT_IN,luin
CCC     END IF
        CALL RDVEC_BATCH_DRV(SCR,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN,IOFFSET_IN)
C
        IOFFSET_OUT = MY_IOFF_LUOUT + ( MY_VEC_IOFF * JOFF ) + 
     &                IOFFSET_SCRATCH
CCC     IF(mpyproc.eq.master) then 
CCC       write(6,*) 'offset for luout',IOFFSET_INT_OUT,luout
CCC     END IF
C
        CALL WTVEC_BATCH_DRV(SCR,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUOUT,IOFFSET_OUT,0)
C
C       count the length of the last copy
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE RDVEC_BATCH_DRV2(SCR,NBATCH_BLK,NBATCH_LEN,NBATCH_INF,
     &                            MY_LUIN,IOFFSET,
     &                            IOFFSET_INT,IVCOFF_IN1)
C
C     read in a batch from a MPI-file --> (hint LUIN) to SCR
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/clunit.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN  = IVCOFF_IN1(IOFFSET_INT+IBLK-1)
          ILEN2 = NBATCH_INF(8,IBLK) 
          IF( ILEN .gt. 0 ) THEN
C           check for consistency
            IF( ILEN2 .ne. ILEN ) THEN
              WRITE(6,*) 'Error in RDVEC_BATCH_DRV2! Different block
     & sizes. I expect them to be equal!'
              CALL  Abend2('Different block sizes in RDVEC_BATCH_DRV2')
            END IF
C           memory offset
            JOFF = NBATCH_INF(6,IBLK)
CCC    IF(MYPROC.eq.MASTER+1) THEN
CCC       WRITE(6,*) ' JOFF for read in IBLK at offset',JOFF, IBLK,
CCC  &                 IOFFSET
CCC    END IF

C
            IF( MY_LUIN .eq. LUC ) THEN
              CALL MPI_FILE_READ_AT(ILU1,IOFFSET,SCR(JOFF),
     &                      ILEN,MPI_REAL8,ISTAT,IERR)
            ELSE IF( MY_LUIN .eq. LUHC ) THEN
              CALL MPI_FILE_READ_AT(ILU2,IOFFSET,SCR(JOFF),
     &                      ILEN,MPI_REAL8,ISTAT,IERR)
            ELSE IF( MY_LUIN .eq. LUSC1 ) THEN
              CALL MPI_FILE_READ_AT(ILU3,IOFFSET,SCR(JOFF),
     &                      ILEN,MPI_REAL8,ISTAT,IERR)
            ELSE IF( MY_LUIN .eq. LUSC2 ) THEN
              CALL MPI_FILE_READ_AT(ILU4,IOFFSET,SCR(JOFF),
     &                      ILEN,MPI_REAL8,ISTAT,IERR)
            ELSE IF( MY_LUIN .eq. LUSC3 ) THEN
              CALL MPI_FILE_READ_AT(ILU5,IOFFSET,SCR(JOFF),
     &                      ILEN,MPI_REAL8,ISTAT,IERR)
            ELSE IF( MY_LUIN .eq. LUSC35 ) THEN
              CALL MPI_FILE_READ_AT(ILU7,IOFFSET,SCR(JOFF),
     &                      ILEN,MPI_REAL8,ISTAT,IERR)
            END IF
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        ENDDO
C      /\ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE WTVEC_BATCH_DRV2(SCR,NBATCH_BLK,NBATCH_LEN,NBATCH_INF,
     &                            MY_LUOUT,IOFFSET,IOFFSET_INT,
     &                            IVCOFF_OUT,NO_CHECK)
C
C     write batch from SCR to MPI-file --> (hint LUOUT)
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/clunit.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_OUT(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF, LENGTH
      JOFF = 0
C
C       loop over all blocks in that batch, do a norm check
C
        DO IBLK = 1, NBATCH_BLK
C
          JOFF = NBATCH_INF(6,IBLK)
          LENGTH = NBATCH_INF(8,IBLK)
C         check norm of vector
          IF( NO_CHECK .eq. 0 ) THEN
            XXX = DDOT(LENGTH,SCR(JOFF),1,SCR(JOFF),1)
C
            IF( XXX .eq. 0.0D0 ) THEN
              IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) = 0
              GOTO 100
            END IF
          ELSE
C
            IF( IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) .eq. 0 ) GOTO 100
C
          END IF
C
          IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) = LENGTH
C
          IF( MY_LUOUT .eq. LUC ) THEN
            CALL MPI_FILE_WRITE_AT(ILU1,IOFFSET,SCR(JOFF),
     &                    LENGTH,MPI_REAL8,ISTAT,IERR)
          ELSE IF( MY_LUOUT .eq. LUHC ) THEN
            CALL MPI_FILE_WRITE_AT(ILU2,IOFFSET,SCR(JOFF),
     &                    LENGTH,MPI_REAL8,ISTAT,IERR)
          ELSE IF( MY_LUOUT .eq. LUSC1 ) THEN
            CALL MPI_FILE_WRITE_AT(ILU3,IOFFSET,SCR(JOFF),
     &                    LENGTH,MPI_REAL8,ISTAT,IERR)
          ELSE IF( MY_LUOUT .eq. LUSC2 ) THEN
            CALL MPI_FILE_WRITE_AT(ILU4,IOFFSET,SCR(JOFF),
     &                    LENGTH,MPI_REAL8,ISTAT,IERR)
          ELSE IF( MY_LUOUT .eq. LUSC3 ) THEN
            CALL MPI_FILE_WRITE_AT(ILU5,IOFFSET,SCR(JOFF),
     &                    LENGTH,MPI_REAL8,ISTAT,IERR)
          ELSE IF( MY_LUOUT .eq. LUSC35 ) THEN
            CALL MPI_FILE_WRITE_AT(ILU7,IOFFSET,SCR(JOFF),
     &                    LENGTH,MPI_REAL8,ISTAT,IERR)
          END IF
C
 100      CONTINUE  ! skip zero blocks on file 
          IOFFSET = IOFFSET + LENGTH
C
        ENDDO
C      /\ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE RDVEC_BATCH_DRV3(LUIN,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET,
     &                            IOFFSET_INT_IN,IOFFSET_INT_OUT,
     &                            IVCOFF_IN,IVCOFF_OUT,
     &                            ICOUNT_ACT)
C
C     Written by  S. Knecht         -      May 21 2007
C
C**********************************************************************
C
C     read in a batch from a MPI-file LUIN to SCR - C-vector routine
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*) 
      DIMENSION IVCOFF_IN(*), IVCOFF_OUT(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
        DO IBLK = 1, NBATCH_BLK
C
C         check for need
C
          ILEN2 = NBATCH_INF(8,IBLK) 
C
          IF( ILEN2 .gt. 0 ) THEN
C
C           count the active blocks for this batch
C
            ICOUNT_ACT = ICOUNT_ACT + 1
C
            IF( IVCOFF_IN( IOFFSET_INT_IN + ICOUNT_ACT - 1 ) .gt.0)THEN
C
C             write ILEN2 into file array for LUOUT
C
              IVCOFF_OUT( IOFFSET_INT_OUT + IBLK - 1 ) = ILEN2
C
C             memory offset
C
              JOFF = NBATCH_INF(6,IBLK)
C
CSK          WRITE(LUWRT,*) ' JOFF for read in IBLK at offset',JOFF, IBLK,
CSK     &                 IOFFSET
C
C           read vector
C
              CALL MPI_FILE_READ_AT(LUIN,IOFFSET,SCR(JOFF),
     &                              ILEN2,MPI_REAL8,ISTAT,IERR)
C
CSK              WRITE(LUWRT,'(2X,A,1X,I12,1X,A,1X,I6,1X,I6)') 
CSK     & 'Read-in from ILU1 at',IOFFSET,'for block',IOFFSET_INT_OUT +
CSK     &  IBLK - 1, JOFF
CSK              CALL WRTMATMN(SCR(JOFF),1,ILEN2,1,ILEN2,LUWRT)
C
            END IF 
C          /\ IVCOFF_IN .gt. 0 ?
          END IF
C        /\ ILEN2 .gt. 0 ?
C
          IOFFSET = IOFFSET + ILEN2
C
C
        ENDDO
C      /\ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WTVEC_BATCH_DRV3(SCR,NBATCH_BLK,NBATCH_INF,
     &                            LUOUT,IOFFSET,IOFFSET_INT,
     &                            IVCOFF_OUT,NO_CHECK,ISCRATCH_SP,
     &                            IBLOCKL)
C
C     Written by  S. Knecht         -      May 21 2007
C
C**********************************************************************
C
C
C     write batch from SCR to MPI-file LUOUT - C-vector routine
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_OUT(*)
      DIMENSION IBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET, ISCRATCH_SP
      INTEGER JOFF, LENGTH
C
      JOFF   = 0
      LENGTH = 0
CSK      WRITE(LUWRT,*) 'calculating length start block is', IOFFSET_INT
C
C     loop over all blocks in that batch, do a norm check
C
      DO IBLK = 1, NBATCH_BLK
C
        JOFF = NBATCH_INF(6,IBLK)
CSK        WRITE(LUWRT,*) 'calculating length for block',
CSK     &                  IOFFSET_INT + IBLK - 1
        LENGTH = IBLOCKL( IOFFSET_INT + IBLK - 1 )
C       check norm of vector
        IF( NO_CHECK .eq. 0 ) THEN
          XXX = 0.0D0
          XXX = DDOT(LENGTH,SCR(JOFF),1,SCR(JOFF),1)
C
          IF( XXX .eq. 0.0D0 ) THEN
            IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) = 0
            GOTO 100
          END IF
        ELSE
C
          IF( IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) .eq. 0 ) GOTO 100
C
        END IF
C
CSK          WRITE(LUWRT,'(2X,A,1X,I6,1X,A,1X,I14,1X,I6)') 'THIS IS block',
CSK     &         IOFFSET_INT + IBLK - 1,'to go on pos (JOFF)',IOFFSET,JOFF
CSK          CALL WRTMATMN(SCR(JOFF),1,LENGTH,1,LENGTH,LUWRT)
C
          CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET,SCR(JOFF),
     &                           LENGTH,MPI_REAL8,ISTAT,IERR)
C
C
 100    CONTINUE  ! skip zero blocks on file
        IOFFSET = IOFFSET + LENGTH
        ISCRATCH_SP = ISCRATCH_SP + LENGTH
C
      END DO
C    /\ loop over blocks in a batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WTVEC_BATCH_DRV3_C(SCR,NBATCH_BLK,NBATCH_INF,
     &                              LUOUT,IOFFSET,IOFFSET_INT,
     &                              IVCOFF_OUT,NO_CHECK,ISCRATCH_SP,
     &                              IBLOCKL)
C
C     Written by  S. Knecht         -      June 22 2007
C
C**********************************************************************
C
C
C     write batch from SCR to MPI-file LUOUT - C-vector routine
C
C     COMPLEX vector part
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_OUT(*)
      DIMENSION IBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET, ISCRATCH_SP
      INTEGER JOFF, LENGTH
C
      JOFF   = 0
      LENGTH = 0
C
C     loop over all blocks in that batch, do a norm check
C
      DO IBLK = 1, NBATCH_BLK
C
        JOFF = NBATCH_INF(6,IBLK)
        LENGTH = IBLOCKL( IOFFSET_INT + IBLK - 1 - MY_ACT_BLK_ALL )
C       check norm of vector
        IF( NO_CHECK .eq. 0 ) THEN
          XXX = 0.0D0
          XXX = DDOT(LENGTH,SCR(JOFF),1,SCR(JOFF),1)
C
          IF( XXX .eq. 0.0D0 ) THEN
            IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) = 0
            GOTO 100
          END IF
        ELSE
C
          IF( IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) .eq. 0 ) GOTO 100
C
        END IF
C
CSK          WRITE(LUWRT,'(2X,A,1X,I6,1X,A,1X,I14,1X,I6)') 'THIS IS block',
CSK     &         IOFFSET_INT + IBLK - 1,'to go on pos (JOFF)',IOFFSET,JOFF
CSK          CALL WRTMATMN(SCR(JOFF),1,LENGTH,1,LENGTH,LUWRT)
C
          CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET,SCR(JOFF),
     &                           LENGTH,MPI_REAL8,ISTAT,IERR)
C
C
 100    CONTINUE  ! skip zero blocks on file
        IOFFSET = IOFFSET + LENGTH
        ISCRATCH_SP = ISCRATCH_SP + LENGTH
C
      END DO
C     ^ loop over blocks in a batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE RDVEC_BATCH_DRV4(MY_LUIN,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET,IOFFSET_INT,IVCOFF_IN1,
     &                            NUM_A_B)
C
C     read in a batch from a MPI-file MY_LUIN to segment SCR
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
      NUM_A_B = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
            NUM_A_B = NUM_A_B + 1 
C
            ILEN  = IVCOFF_IN1( IOFFSET_INT + NUM_A_B -1 )
C
            IF( ILEN .gt. 0 ) THEN
C             memory offset
              JOFF = NBATCH_INF(6,IBLK)
C
CSK              WRITE(LUWRT,*) 'JOFF,ILEN2,ILEN again for block 
CSK     & IOFFSET_INT+ NUM_A_B -1',JOFF,ILEN2,ILEN,IOFFSET_INT+ NUM_A_B -1
              CALL MPI_FILE_READ_AT(MY_LUIN,IOFFSET,SCR(JOFF),
     &                              ILEN2,MPI_REAL8,ISTAT,IERR)
CSK              WRITE(LUWRT,*) ' my block'
CSK              CALL WRTMATMN(SCR(JOFF),1,ILEN2,1,ILEN2,LUWRT)
            ENDIF 
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WTVEC_BATCH_DRV4(MY_LUOUT,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET,IOFFSET_INT,IVCOFF_IN1,
     &                            NUM_A_B)
C
C     write batch from segment SCR to MPI-file MY_LUOUT
C
C     active blocks on the MPI-file are flagged by a nonzero entry in 
C     file list
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
      NUM_A_B = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
C           memory offset
            JOFF  = NBATCH_INF(6,IBLK)
            NUM_A_B = NUM_A_B + 1
C           check norm of vector
            XXX = 0.0D0
            XXX = DDOT(ILEN2,SCR(JOFF),1,SCR(JOFF),1)
C
            IF( XXX .eq. 0.0D0 ) THEN
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 0
              GOTO 100
            ELSE
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 1
            END IF
C
CSK              WRITE(LUWRT,*) 'JOFF,ILEN2,ILEN again for block 
CSK     & IOFFSET_INT+ NUM_A_B -1',JOFF,ILEN2,ILEN,IOFFSET_INT+ NUM_A_B -1
              CALL MPI_FILE_WRITE_AT(MY_LUOUT,IOFFSET,SCR(JOFF),
     &                               ILEN2,MPI_REAL8,ISTAT,IERR)
C
 100        CONTINUE
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WTVEC_BATCH_DRV4SP(MY_LUOUT,SCR,NBATCH_BLK,NBATCH_INF,
     &                              IOFFSET,IOFFSET_INT,IVCOFF_IN1,
     &                              NUM_A_B,IDEBUGPRNT)
C
C     write batch from segment SCR to MPI-file MY_LUOUT
C
C     active blocks on the MPI-file are flagged by a nonzero entry in 
C     file list
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
      NUM_A_B = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
C           memory offset
            JOFF  = NBATCH_INF(6,IBLK)
            NUM_A_B = NUM_A_B + 1
C           check norm of vector
            XXX = 0.0D0
            XXX = DDOT(ILEN2,SCR(JOFF),1,SCR(JOFF),1)
C
            IF( XXX .eq. 0.0D0 ) THEN
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 0
              GOTO 100
            ELSE
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 1
            END IF
C
            IF( IDEBUGPRNT .ne. 0 ) THEN
              WRITE(LUWRT,'(A,1X,I6,1X,I6,1X,I6,1X,A,1X,I14)') 'JOFF, 
     & ILEN2 for block',JOFF,ILEN2,IOFFSET_INT+ NUM_A_B -1,'at off',
     & IOFFSET
              CALL WRTMATMN(SCR(JOFF),1,ILEN2,1,ILEN2,LUWRT)
            END IF
            CALL MPI_FILE_WRITE_AT(MY_LUOUT,IOFFSET,SCR(JOFF),
     &                             ILEN2,MPI_REAL8,ISTAT,IERR)
C
 100        CONTINUE
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE RDVEC_BATCH_DRV5(MY_LUIN,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET)
C
C     read in a batch from a MPI-file MY_LUIN to segment SCR
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION SCR(*), NBATCH_INF(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
C
      JOFF = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
C           memory offset
C
            JOFF = NBATCH_INF(6,IBLK)
C
            CALL MPI_FILE_READ_AT(MY_LUIN,IOFFSET,SCR(JOFF),
     &                            ILEN2,MPI_REAL8,ISTAT,IERR)
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE BLK_DSCAL(FAC,SCR,NBATCH_BLK,NBATCH_INF,
     &                     IOFFSET,IVCOFF)
C
C     scaling of vector SCR with factor FAC using a blocked vector 
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION NBATCH_INF(8,*), IVCOFF(*), SCR(*) 
      INTEGER JOFF,ILEN
C
      JOFF = 0
      ILEN = 0
      DO ISBLK = 1, NBATCH_BLK
C
        ILEN = IVCOFF( IOFFSET + ISBLK - 1 )
        IF( ILEN .gt. 0 ) THEN
          JOFF = NBATCH_INF(6,ISBLK)
          CALL DSCAL( ILEN, FAC, SCR(JOFF), 1 )
        END IF
C
      END DO
C
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_BATCH_PAR2(SCR,NBATCH,LBATCH,LEBATCH,
     &                             I1BATCH,IBATCH,MY_VEC_IOFF,
     &                             MY_IOFF_LUIN,MY_IOFF_LUOUT,
     &                             NBLK_ACT,LUIN,LUOUT,JOFF,
     &                             IVCOFF_IN,IVCOFF_OUT)
C
C     copy vector from file LUIN to file LUOUT --> (hint LUIN, LUOUT)
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION SCR(*),LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION IVCOFF_IN(*), IVCOFF_OUT(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN  = 0
      IOFFSET_OUT  = 0
      IOFFSET_INT_IN  = 0
      IOFFSET_INT_OUT  = 0
C
      DO ISBATCH = 1, NBATCH 
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(SCR,LEBATCH(ISBATCH))
C
        IOFFSET_IN = MY_IOFF_LUIN + ( MY_VEC_IOFF * JOFF ) +
     &               IOFFSET_SCRATCH
        IOFFSET_INT_IN = ( NBLK_ACT * JOFF ) + NUM_BLK + 1
C
        CALL RDVEC_BATCH_DRV2(SCR,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN,IOFFSET_IN,IOFFSET_INT_IN,
     &                       IVCOFF_IN)
C
        IOFFSET_OUT = MY_IOFF_LUOUT + ( MY_VEC_IOFF * JOFF ) + 
     &                IOFFSET_SCRATCH
        IOFFSET_INT_OUT = ( NBLK_ACT * JOFF ) + NUM_BLK + 1
C
        CALL WTVEC_BATCH_DRV2(SCR,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        LUOUT,IOFFSET_OUT,IOFFSET_INT_OUT,
     &                        IVCOFF_OUT,0)
C
C       count the length of the last copy
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + LBATCH(ISBATCH)
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ORTHG_VEC_BATCH_PAR(VEC1,VEC2,WORK_SP,SCRRED,SCRRED2,
     &                               NROOT,JOFF,
     &                               NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                               IBATCH,
     &                               MY_VEC_IOFF,
     &                               MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                               MY_IOFF_LUIN3,MY_IOFF_LUOUT,
     &                               NBLK_ACT,LUIN1,LUIN2,
     &                               LUIN3,LUOUT,LUWRTOUT)
C
C     calculate dot products between vector on file LUIN1 and LUIN2 
C     resp. LUIN1 and LUIN3. Store vector VEC2 on LUOUT 
C                            ==> LUIN1(VEC2) --> LUOUT(VEC2)
C
C     dot products SUM_batches ( VEC1 * VEC2 ) are reduced on every 
C     process in MPI_COMM_WORLD
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: JOFF = JVEC - 1
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION WORK_SP(*), SCRRED(*), SCRRED2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER ISTRED, ISTRED2
      CHARACTER*12 SECTID, WALLTID
C     initialize scratch vector offsets
      ISTRED  = 0
      ISTRED2 = 0
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_IN_LUIN3  = 0
      IOFFSET_LUOUT     = 0
C     set scratch vectors to zero
      CALL DZERO(SCRRED,NROOT)
      IF( JOFF .gt. 0 ) CALL DZERO(SCRRED2,JOFF)
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       special treatment for LUIN1 since we need to skip some vectors
C       for simplicity this is done on the outside
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH
CSK        IF(MYPROC.EQ.MASTER) THEN
CSK         WRITE(6,*) 'This is my OFFSET for LUIN1 special',
CSK     &               IOFFSET_IN_LUIN1, myproc
CSK        END IF
C
        CALL RDVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN1,IOFFSET_IN_LUIN1)

C
CSK          IF(MYPROC.EQ.MASTER+1) THEN
CSK            WRITE(6,*) 'initial VEC2 on ILU7'
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C
        DO 100 IROOT = 1, NROOT
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IROOT - 1 ) * MY_VEC_IOFF
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
CCC       IF(MYPROC.EQ.MASTER+1) THEN
CCC        WRITE(6,*) 'This is my OFFSET for LUIN2',IOFFSET_IN_LUIN2
CCC       END IF
CCC       IF(MYPROC.EQ.MASTER+1) THEN
CCC         WRITE(6,*) 'VEC1 before read in from ILU1 bt',myproc,ISBATCH
CCC         CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CCC       END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN2,IOFFSET_IN_LUIN2)
C
CCC       IF(MYPROC.EQ.MASTER+1) THEN
CCC         WRITE(6,*) 'VEC1 read in as from ILU1for bt',myproc,ISBATCH
CCC         CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CCC       END IF
          
          IF( IROOT .eq. 1 ) ISTRED = IROOT+(JOFF)*2*NROOT
C         SCRRED(IROOT) == SCRRED(IROOT) - SUM ( VEC1 * VEC2 )
          SCRRED(IROOT) = 
     &    SCRRED(IROOT) - DDOT(LEBATCH(ISBATCH), VEC1, 1, VEC2, 1)
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 KVEC = 1, JOFF
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( KVEC - 1 ) * MY_VEC_IOFF
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3)
C
          IF( KVEC .eq. 1 ) ISTRED2 = NROOT+KVEC+(JOFF)*2*NROOT
C         SCRRED2(KVEC) == SCRRED2(KVEC) - SUM ( VEC1 * VEC2 )
          SCRRED2(KVEC) = 
     &    SCRRED2(KVEC) - DDOT(LEBATCH(ISBATCH), VEC1, 1, VEC2, 1)
C
  200   CONTINUE
C
C
C       write vector VEC2 from batch IBATCH to file LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + ( MY_VEC_IOFF * JOFF ) + 
     &                      IOFFSET_SCRATCH
C
        CALL WTVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUOUT,IOFFSET_LUOUT,0)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
C     individual dot products need to be communicated
C
      starttime = MPI_WTIME()
      CALL REDVEC(SCRRED,WORK_SP(ISTRED),NROOT,2,
     &            MPI_SUM,MPI_COMM_WORLD,-1)
C
      IF( JOFF .gt. 0 ) THEN
        CALL REDVEC(SCRRED2,WORK_SP(ISTRED2),JOFF,2,
     &              MPI_SUM,MPI_COMM_WORLD,-1)
      END IF
C
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      WRITE(LUWRTOUT,9450) WALLTID, ( NROOT + JOFF )
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
C
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      IF( MYPROC.eq.MASTER) THEN
CSK      WRITE(6,*) 'THIS IS WHAT I HAVE REDUCED TO (MYPROC)',MYPROC
CSK      CALL WRTMAT(WORK_SP(ISTRED),1,NROOT,1,NROOT)
CSK      ENDIF
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ORTHG_VEC_BATCH_PAR2(VEC1,VEC2,FACIN1,FACIN2,NINVEC1,
     &                                NINVEC2,
     &                                NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                                IBATCH,
     &                                MY_VEC_IOFF,
     &                                MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                                MY_IOFF_LUIN3,MY_IOFF_LUOUT1,
     &                                MY_IOFF_LUOUT2,
     &                                NBLK_ACT,LUIN1,LUIN2,
     &                                LUIN3,LUOUT1,LUOUT2,LUWRTOUT,
     &                                REDFAC,SCAL_ARRAY,NSCALING,IADD)
C
C     add vectors on file LUIN1 and LUIN2 resp. LUIN2 and LUIN3. 
C     Store resulting vector VEC2 on LUOUT with scaling factor
C     FACOFSCAL.
C     LUOUT(VEC2) = LUIN2(VEC2) x FACOFSCAL
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION SCAL_ARRAY(*), FACIN1(*), FACIN2(1)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT1, MY_IOFF_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT1, IOFFSET_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
C     initialize scratch factor 
      REDFAC = 0.0D0
C     initialize scratch offsets
      IOFFSET_SCRATCH  = 0
      IOFFSET_IN_LUIN1 = 0
      IOFFSET_IN_LUIN2 = 0
      IOFFSET_IN_LUIN3 = 0
      IOFFSET_LUOUT1   = 0
      IOFFSET_LUOUT2   = 0
C
CSK      WRITE(6,*) 'THIS IS WORK_SP(1+(JVEC-1)*2*NROOT) inside',myproc
CSK      CALL WRTMAT(FACIN1,1,NINVEC1,1,NINVEC1)
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       find offset 
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH + 
     &                   ( NINVEC2 ) * MY_VEC_IOFF
C
        CALL RDVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN2,IOFFSET_IN_LUIN2)
C
C       calculate scalar products with vectors on LUIN1
C   
        DO 100 IVEC = 1, NINVEC1
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 ) * MY_VEC_IOFF
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN1',IOFFSET_IN_LUIN1
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN1,IOFFSET_IN_LUIN1)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IVEC)
CSK          IF(MYPROC.EQ.MASTER+1) THEN
CSK            WRITE(6,*) 'VEC2 before DAXPY for LUIN1',
CSK     &                  Myproc,IVEC
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            WRITE(6,*) 'FACIN1(IVEC) = ',FACIN1(IVEC)
CSK          END IF
CSK          FACSCALE = FACIN1(IVEC)
          IF( IVEC .eq. 1 ) THEN
            CALL DAXPY(LEBATCH(ISBATCH),FACIN1(1),VEC1,1,VEC2,1)
          ELSE
            CALL DAXPY(LEBATCH(ISBATCH),FACIN1(IVEC),VEC1,1,VEC2,1)
          END IF
CSK       IF(MYPROC.EQ.MASTER+1) THEN
CSK       IF (IVEC .eq. NINVEC1 ) THEN
CSK         WRITE(6,*) 'VEC2 after DAXPY == after first ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK         CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK       END IF
CSK       END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 IVEC = 1, NINVEC2
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 ) * MY_VEC_IOFF
C
C          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
C          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IVEC)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            WRITE(6,*) 'FACIN2(1) = ',FACIN2(1)
CSK          END IF
          IF( IVEC .eq. 1 ) THEN 
            CALL DAXPY(LEBATCH(ISBATCH),FACIN2(1),VEC1,1,VEC2,1)
          ELSE
            CALL DAXPY(LEBATCH(ISBATCH),FACIN2(IVEC),VEC1,1,VEC2,1)
          END IF
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC2 after DAXPY == after sec ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C
  200   CONTINUE
C
CSK          IF(MYPROC.EQ.MASTER+1) THEN
CSK            WRITE(6,*) 'VEC2 after DAXPY == after sec ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C       update REDFAC == REDFAC + SUM ( VEC2 * VEC2 )
        REDFAC = REDFAC + DDOT( LEBATCH(ISBATCH), VEC2, 1, VEC2, 1 )
CSK        WRITE(6,*) 'REDFAC = REDFAC + DDOT ...',REDFAC,MYPROC
C
C       write vector VEC2 from batch IBATCH to file LUIN2
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH + 
     &                   ( NINVEC2 ) * MY_VEC_IOFF
C
        CALL WTVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUIN2,IOFFSET_IN_LUIN2,0)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
C     communicate REDFAC to find correct scaling factor
      
      FACTOR = 0.0D0
      CALL REDVEC(REDFAC,FACTOR,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
CSK      WRITE(6,*) 'REDUCED FACTOR IN PART 3.2',FACTOR
C
CSK      SCALING = 0.0D0
      SCALING = 1.0D0/SQRT(FACTOR)
C      WRITE(6,*) 'SCALING FACTOR IN PART 3.2',SCALING
C
C     scaling part 1
C
      CALL DSCAL(NSCALING, SCALING, SCAL_ARRAY, 1 )
C
CSK      WRITE(LUWRTOUT,*) 'THIS IS MY SCALED WORK PART',MYPROC
CSK      CALL WRTMAT(SCAL_ARRAY,1,NSCALING,1,NSCALING)
C
C     initialize scratch offset again
      IOFFSET_SCRATCH = 0
C
C     we loop again
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO( VEC2, LEBATCH(ISBATCH) )
C
C       read vector from file LUIN2
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH + 
     &                   ( NINVEC2 ) * MY_VEC_IOFF
C
        CALL RDVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUIN2,IOFFSET_IN_LUIN2)
C
C       scaling part 2
C
        CALL DSCAL( LEBATCH(ISBATCH), SCALING, VEC2, 1 )
C
C       write scaled vector to files LUOUT1 and LUOUT2
C       writing to LUOUT1 only for ( NINVEC2 + 1 .ne. IADD )
C
        IF( ( NINVEC2 + 1 ) .ne. IADD ) THEN
C
          IOFFSET_LUOUT1     = MY_IOFF_LUOUT1 +
     &                     ( MY_VEC_IOFF * (NINVEC2) ) + IOFFSET_SCRATCH
C
          CALL WTVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUOUT1,IOFFSET_LUOUT1,0)
        END IF
C
C
        IOFFSET_LUOUT2     = MY_IOFF_LUOUT2 + 
     &                     ( MY_VEC_IOFF * (NINVEC2) ) + IOFFSET_SCRATCH
CSK           WRITE(6,*) 'This is my OFFSET for writing on ILU3 and what',
CSK     &                 IOFFSET_LUOUT2
CSK           CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
C
        CALL WTVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUOUT2,IOFFSET_LUOUT2,0)
C
CSK           WRITE(6,*) 'This is what i have written',
CSK     &                 myproc
CSK           CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C    /\ loop over batches 2
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ORTHG_VEC_BATCH_PAR3(VEC1,VEC2,FACIN1,FACIN2,FACIN3,
     &                                NROOT,JOFF,
     &                                NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                                IBATCH,
     &                                MY_VEC_IOFF,
     &                                MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                                MY_IOFF_LUIN3,MY_IOFF_LUOUT,
     &                                NBLK_ACT,LUIN1,LUIN2,
     &                                LUIN3,LUOUT,LUWRTOUT,IADD)
C
C     calculate vector products between vector on file LUIN1 and LUIN2 
C     resp. LUIN1 and LUIN3. Store vector VEC2 on LUOUT 
C                            ==> LUIN1(VEC2) --> LUOUT(VEC2)
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: JOFF = JVEC - 1
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION FACIN1(*), FACIN3(1)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
C     initialize scratch offsets
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_IN_LUIN3  = 0
      IOFFSET_LUOUT     = 0
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       special treatment for LUIN1 since we need to skip some vectors
C       for simplicity this is done on the outside
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH
CSK        IF(MYPROC.EQ.MASTER) THEN
CSK         WRITE(6,*) 'This is my OFFSET for LUIN1 special',
CSK     &               IOFFSET_IN_LUIN1, myproc
CSK        END IF
C
        CALL RDVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN1,IOFFSET_IN_LUIN1)
C
        DO 100 IROOT = 1, NROOT
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IROOT - 1 ) * MY_VEC_IOFF
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN2',IOFFSET_IN_LUIN2
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN2,IOFFSET_IN_LUIN2)
C
C         scaling of VEC2
          IF( IROOT .eq. 1 ) THEN 
            CALL DSCAL( LEBATCH(ISBATCH), FACIN2, VEC2, 1 )
          END IF
C         VEC2 == VEC2 + VEC1 * FACIN1(IROOT)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY for LUIN2',
CSK     &                  Myproc,IVEC
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
          CALL DAXPY( LEBATCH(ISBATCH), FACIN1(IROOT), VEC1, 1, VEC2, 1)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK          IF (IVEC .eq. NINVEC1 ) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 after DAXPY == after first ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
CSK          END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 KVEC = 1, JOFF
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( KVEC - 1 ) * MY_VEC_IOFF
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IROOT)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY for LUIN3',
CSK     &                  Myproc,IVEC
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
          IF ( KVEC .eq. 1 ) THEN
            CALL DAXPY( LEBATCH(ISBATCH),FACIN3(1), VEC1, 1, VEC2, 1)
          ELSE 
            CALL DAXPY( LEBATCH(ISBATCH),FACIN3(KVEC), VEC1, 1, VEC2, 1)
          END IF
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK          IF (KVEC .eq. JOFF ) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 after DAXPY == after sec ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
CSK          END IF
C
  200   CONTINUE
C
C
C       write vector VEC2 from batch IBATCH to file LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + ( MY_VEC_IOFF * JOFF ) + 
     &                      IOFFSET_SCRATCH
        CALL WTVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUOUT,IOFFSET_LUOUT,0)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CALC_SUBSPACE_H(VEC1,VEC2,SCRRED,IMUSTRED,IVEC,NROOT,
     &                           ISTRED,
     &                           NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                           MY_VEC_IOFF,MY_IOFF_LUIN1,
     &                           MY_IOFF_LUIN2,MY_IOFF_LUIN3,
     &                           MY_IOFF_LUIN4,NBLK_ACT,LUIN1,LUIN2,
     &                           LUIN3,LUIN4)
C
C     calculate dot products between vectors on file LUIN1 / LUIN2 
C     and LUIN3 resp. LUIN4
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: IVEC = IVEC
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION SCRRED(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN4
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_SPECIAL
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_SPECIAL
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3, IOFFSET_IN_LUIN4
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
C     initialize scratch offsets
      IOFFSET_SCRATCH     = 0
      IOFFSET_IN_SPECIAL  = 0
      IOFFSET_IN_LUIN3    = 0
      IOFFSET_IN_LUIN4    = 0
      LUIN_SPECIAL = 0
      MY_IOFF_SPECIAL     = 0
C     set correct input file and offsets
      IF( IVEC .le. NROOT ) THEN
        LUIN_SPECIAL = LUIN1
        MY_IOFF_SPECIAL  = MY_IOFF_LUIN1
      ELSE
        LUIN_SPECIAL = LUIN2
        MY_IOFF_SPECIAL  = MY_IOFF_LUIN2
      END IF
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       special treatment for LUIN1 / LUIN2 since we need to skip some vectors
C       for simplicity this is done on the outside
C
        IOFFSET_IN_SPECIAL = MY_IOFF_SPECIAL + IOFFSET_SCRATCH
CSK        IF(MYPROC.EQ.MASTER) THEN
CSK         WRITE(6,*) 'This is my OFFSET for LUIN1 / LUIN2 special',
CSK     &               IOFFSET_IN_SPECIAL, LUIN_SPECIAL, myproc
CSK        END IF
C
        CALL RDVEC_BATCH_DRV(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN_SPECIAL,IOFFSET_IN_SPECIAL)
C
        DO 100 JVEC = 1, MIN(IVEC,NROOT)
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( JVEC - 1 ) * MY_VEC_IOFF
C
C         read in batch ISBATCH from LUIN3 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3)
          IJ = IVEC * ( IVEC - 1 ) / 2 + JVEC
C         SCRRED(IJ) = SCRRED(IJ) + SUM ( VEC1 * VEC2 )
          SCRRED(IJ) = SCRRED(IJ) + 
     &                 DDOT( LEBATCH(ISBATCH), VEC1, 1, VEC2, 1)
C
          IF (ISBATCH .eq. 1 ) THEN
            IMUSTRED = IMUSTRED + 1
            IF( IVEC .eq. 1 .and. JVEC .eq. 1 ) ISTRED = IJ
          END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN4
C
        ILOOP = 0
C
        DO 200 JVEC = ( NROOT + 1), IVEC
C
          ILOOP = ILOOP + 1
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN4 = MY_IOFF_LUIN4 + IOFFSET_SCRATCH +
     &                     ( ILOOP - 1 ) * MY_VEC_IOFF
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN4',IOFFSET_IN_LUIN4
CSK          END IF
C
          CALL RDVEC_BATCH_DRV(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN4,IOFFSET_IN_LUIN4)
C
          IJ = IVEC * ( IVEC - 1 ) / 2 + JVEC
C         SCRRED(IJ) = SCRRED(IJ) + SUM ( VEC1 * VEC2 )
          SCRRED(IJ) = SCRRED(IJ) +
     &                 DDOT( LEBATCH(ISBATCH), VEC1, 1, VEC2, 1)
C
          IF (ISBATCH .eq. 1 ) IMUSTRED = IMUSTRED + 1
C
  200   CONTINUE
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ORTHG_VEC_BATCH_PAR_USE_IA(VEC1,VEC2,WORK_SP,SCRRED,
     &                                      SCRRED2,NROOT,JOFF,
     &                                      NBATCH,LBATCH,LEBATCH,
     &                                      I1BATCH,IBATCH,MY_VEC_IOFF,
     &                                      MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                                      MY_IOFF_LUIN3,MY_IOFF_LUOUT,
     &                                      NBLK_ACT,LUIN1,LUIN2,
     &                                      LUIN3,LUOUT,LUWRTOUT,
     &                                      MY_INT_IOFF_LUIN1,
     &                                      IVCOFF_IN1,IVCOFF_IN2,
     &                                      IVCOFF_IN3,IVCOFF_OUT)
C
C     calculate dot products between vector on file LUIN1 and LUIN2 
C     resp. LUIN1 and LUIN3. Store vector VEC2 on LUOUT 
C                            ==> LUIN1(VEC2) --> LUOUT(VEC2)
C
C     dot products SUM_batches ( VEC1 * VEC2 ) are reduced on every 
C     process in MPI_COMM_WORLD
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: JOFF = JVEC - 1
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION WORK_SP(*), SCRRED(*), SCRRED2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION IVCOFF_IN1(*), IVCOFF_IN2(*), IVCOFF_IN3(*)
      DIMENSION IVCOFF_OUT(*)
      REAL*8 BLK_DDOT
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK, ISTRED, ISTRED2
      CHARACTER*12 SECTID, WALLTID
C     initialize scratch vector offsets
      ISTRED  = 0
      ISTRED2 = 0
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_IN_LUIN3  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_IN_INT_LUIN1  = 0
      IOFFSET_IN_INT_LUIN2  = 0
      IOFFSET_IN_INT_LUIN3  = 0
      IOFFSET_INT_LUOUT     = 0
C     set scratch vectors to zero
      CALL DZERO(SCRRED,NROOT)
      IF( JOFF .gt. 0 ) CALL DZERO(SCRRED2,JOFF)
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       special treatment for LUIN1 since we need to skip some vectors
C       for simplicity this is done on the outside
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH
        IOFFSET_IN_INT_LUIN1 = MY_INT_IOFF_LUIN1 + NUM_BLK
C
        CALL RDVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN1,IOFFSET_IN_LUIN1,
     &                       IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
C
CSK          IF(MYPROC.EQ.MASTER+1) THEN
CSK            WRITE(6,*) 'initial VEC2 on ILU7'
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C
C
        DO 100 IROOT = 1, NROOT
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IROOT - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN2 = NUM_BLK +
     &                     ( IROOT - 1 ) * NBLK_ACT + 1
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          LUIN2,IOFFSET_IN_LUIN2,
     &                          IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
CCC       IF(MYPROC.EQ.MASTER+1) THEN
CCC         WRITE(6,*) 'VEC1 read in as from ILU1for bt',myproc,ISBATCH
CCC         CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CCC       END IF
          
          IF( IROOT .eq. 1 ) ISTRED = IROOT+(JOFF)*2*NROOT
C         SCRRED(IROOT) == SCRRED(IROOT) - SUM ( VEC1 * VEC2 )
          SCRRED(IROOT) = 
     &    SCRRED(IROOT) - BLK_DDOT(VEC1,VEC2,LBATCH(ISBATCH),
     &                             IBATCH(1,I1BATCH(ISBATCH)),
     &                             IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 KVEC = 1, JOFF
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( KVEC - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN3 = NUM_BLK +
     &                     ( KVEC - 1 ) * NBLK_ACT + 1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
CSK          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          LUIN3,IOFFSET_IN_LUIN3,
     &                          IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
C
          IF( KVEC .eq. 1 ) ISTRED2 = NROOT+KVEC+(JOFF)*2*NROOT
C         SCRRED2(KVEC) == SCRRED2(KVEC) - SUM ( VEC1 * VEC2 )
          SCRRED2(KVEC) = 
     &    SCRRED2(KVEC) - BLK_DDOT(VEC1,VEC2,LBATCH(ISBATCH),
     &                             IBATCH(1,I1BATCH(ISBATCH)),
     &                             IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
C
  200   CONTINUE
C
C
C       write vector VEC2 from batch IBATCH to file LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + ( MY_VEC_IOFF * JOFF ) + 
     &                      IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = ( NBLK_ACT * JOFF ) + NUM_BLK + 1
C
        CALL WTVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        LUOUT,IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                        IVCOFF_OUT,0)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK = NUM_BLK + LBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
C     individual dot products need to be communicated
C
      starttime = MPI_WTIME()
      CALL REDVEC(SCRRED,WORK_SP(ISTRED),NROOT,2,
     &            MPI_SUM,MPI_COMM_WORLD,-1)
C
      IF( JOFF .gt. 0 ) THEN
        CALL REDVEC(SCRRED2,WORK_SP(ISTRED2),JOFF,2,
     &              MPI_SUM,MPI_COMM_WORLD,-1)
      END IF
C
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      WRITE(LUWRTOUT,9450) WALLTID, ( NROOT + JOFF )
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
C
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      IF( MYPROC.eq.MASTER) THEN
CSK      WRITE(6,*) 'THIS IS WHAT I HAVE REDUCED TO (MYPROC)',MYPROC
CSK      CALL WRTMAT(WORK_SP(ISTRED),1,NROOT,1,NROOT)
CSK      ENDIF
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
C
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ORTHG_VEC_BATCH_PAR2_USE_IA(VEC1,VEC2,FACIN1,FACIN2,
     &                                      NINVEC1,NINVEC2,
     &                                      NBATCH,LBATCH,LEBATCH,
     &                                      I1BATCH,IBATCH,MY_VEC_IOFF,
     &                                      MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                                      MY_IOFF_LUIN3,
     &                                      MY_IOFF_LUOUT1,
     &                                      MY_IOFF_LUOUT2,
     &                                      NBLK_ACT,LUIN1,LUIN2,LUIN3,
     &                                      LUOUT1,LUOUT2,LUWRTOUT,
     &                                      REDFAC,SCAL_ARRAY,
     &                                      NSCALING,IADD,
     &                                      IVCOFF_IN1,IVCOFF_IN2,
     &                                      IVCOFF_IN3,IVCOFF_OUT1,
     &                                      IVCOFF_OUT2)
C
C     add vectors on file LUIN1 and LUIN2 resp. LUIN2 and LUIN3. 
C     Store resulting vector VEC2 on LUOUT with scaling factor
C     FACOFSCAL.
C     LUOUT(VEC2) = LUIN2(VEC2) x FACOFSCAL
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION SCAL_ARRAY(*), FACIN1(*), FACIN2(1)
      DIMENSION IVCOFF_IN1(*),IVCOFF_IN2(*), IVCOFF_IN3(*)
      DIMENSION IVCOFF_OUT1(*), IVCOFF_OUT2(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT1, MY_IOFF_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT1, IOFFSET_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C     initialize scratch factor 
      REDFAC = 0.0D0
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH  = 0
      IOFFSET_IN_LUIN1 = 0
      IOFFSET_IN_LUIN2 = 0
      IOFFSET_IN_LUIN3 = 0
      IOFFSET_LUOUT1   = 0
      IOFFSET_LUOUT2   = 0
      IOFFSET_IN_INT_LUIN1 = 0
      IOFFSET_IN_INT_LUIN2 = 0
      IOFFSET_IN_INT_LUIN3 = 0
      IOFFSET_INT_LUOUT1   = 0
      IOFFSET_INT_LUOUT2   = 0
C
CSK      WRITE(6,*) 'THIS IS WORK_SP(1+(JVEC-1)*2*NROOT) inside',myproc
CSK      CALL WRTMAT(FACIN1,1,NINVEC1,1,NINVEC1)
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       find offset 
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH + 
     &                   ( NINVEC2 ) * MY_VEC_IOFF
        IOFFSET_IN_INT_LUIN2 = NUM_BLK + ( NINVEC2 ) * NBLK_ACT + 1
C
        CALL RDVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN2,IOFFSET_IN_LUIN2,
     &                       IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
C       calculate scalar products with vectors on LUIN1
C   
        DO 100 IVEC = 1, NINVEC1
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN1 = NUM_BLK + ( IVEC - 1 ) * NBLK_ACT + 1
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN1',IOFFSET_IN_LUIN1
CSK          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN1,IOFFSET_IN_LUIN1,
     &                         IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IVEC)
CSK          IF(MYPROC.EQ.MASTER+1) THEN
CSK            WRITE(6,*) 'VEC2 before DAXPY for LUIN1',
CSK     &                  Myproc,IVEC
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            WRITE(6,*) 'FACIN1(IVEC) = ',FACIN1(IVEC)
CSK          END IF
CSK          FACSCALE = FACIN1(IVEC)
          IF( IVEC .eq. 1 ) THEN
            CALL SCLVCSUM(FACIN1(1),VEC1,VEC2,LBATCH(ISBATCH),
     &                    IBATCH(1,I1BATCH(ISBATCH)),
     &                    IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
          ELSE
            CALL SCLVCSUM(FACIN1(IVEC),VEC1,VEC2,LBATCH(ISBATCH),
     &                    IBATCH(1,I1BATCH(ISBATCH)),
     &                    IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
          END IF
CSK       IF(MYPROC.EQ.MASTER+1) THEN
CSK       IF (IVEC .eq. NINVEC1 ) THEN
CSK         WRITE(6,*) 'VEC2 after DAXPY == after first ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK         CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK       END IF
CSK       END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 IVEC = 1, NINVEC2
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN3 = NUM_BLK + ( IVEC - 1 ) * NBLK_ACT + 1
C
C          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
C          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3,
     &                         IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IVEC)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            WRITE(6,*) 'FACIN2(1) = ',FACIN2(1)
CSK          END IF
          IF( IVEC .eq. 1 ) THEN 
            CALL SCLVCSUM(FACIN2(1),VEC1,VEC2,LBATCH(ISBATCH),
     &                    IBATCH(1,I1BATCH(ISBATCH)),
     &                    IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
          ELSE
            CALL SCLVCSUM(FACIN2(IVEC),VEC1,VEC2,LBATCH(ISBATCH),
     &                    IBATCH(1,I1BATCH(ISBATCH)),
     &                    IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
          END IF
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC2 after DAXPY == after sec ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C
  200   CONTINUE
C
CSK          IF(MYPROC.EQ.MASTER+1) THEN
CSK            WRITE(6,*) 'VEC2 after DAXPY == after sec ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
C       update REDFAC == REDFAC + SUM ( VEC2 * VEC2 )
        REDFAC = REDFAC + DDOT( LEBATCH(ISBATCH), VEC2, 1, VEC2, 1 )
C?     &           BLK_DDOT(VEC2,VEC2,LBATCH(ISBATCH),
C?     &                    IBATCH(1,I1BATCH(ISBATCH)),
C?     &                    IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
               
CSK        WRITE(6,*) 'REDFAC = REDFAC + DDOT ...',REDFAC,MYPROC
C
C       write vector VEC2 from batch IBATCH to file LUIN2
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH + 
     &                   ( NINVEC2 ) * MY_VEC_IOFF
        IOFFSET_IN_INT_LUIN2 = NUM_BLK + ( NINVEC2 ) * NBLK_ACT + 1
C
        CALL WTVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUIN2,IOFFSET_IN_LUIN2,
     &                      IOFFSET_IN_INT_LUIN2,IVCOFF_IN2,0)
C
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK = NUM_BLK + LBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
C     communicate REDFAC to find correct scaling factor
      
      FACTOR = 0.0D0
      CALL REDVEC(REDFAC,FACTOR,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
CSK      WRITE(6,*) 'REDUCED FACTOR IN PART 3.2',FACTOR
C
CSK      SCALING = 0.0D0
      SCALING = 1.0D0/SQRT(FACTOR)
C      WRITE(6,*) 'SCALING FACTOR IN PART 3.2',SCALING
C
C     scaling part 1
C
      CALL DSCAL(NSCALING, SCALING, SCAL_ARRAY, 1 )
C
CSK      WRITE(LUWRTOUT,*) 'THIS IS MY SCALED WORK PART',MYPROC
CSK      CALL WRTMAT(SCAL_ARRAY,1,NSCALING,1,NSCALING)
C
C     initialize scratch offsets again
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
C
C     we loop again
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO( VEC2, LEBATCH(ISBATCH) )
C
C       read vector from file LUIN2
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH + 
     &                   ( NINVEC2 ) * MY_VEC_IOFF
        IOFFSET_IN_INT_LUIN2 = NUM_BLK + ( NINVEC2 ) * NBLK_ACT + 1
C
        CALL RDVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUIN2,IOFFSET_IN_LUIN2,
     &                      IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
C       scaling part 2
C
        CALL BLK_DSCAL(SCALING,VEC2,LBATCH(ISBATCH),
     &                 IBATCH(1,I1BATCH(ISBATCH)),IOFFSET_IN_INT_LUIN2,
     &                 IVCOFF_IN2)
C
C       write scaled vector to files LUOUT1 and LUOUT2
C       writing to LUOUT1 only for ( NINVEC2 + 1 .ne. IADD )
C
        IF( ( NINVEC2 + 1 ) .ne. IADD ) THEN
C
          IOFFSET_LUOUT1     = MY_IOFF_LUOUT1 +
     &                     ( MY_VEC_IOFF * (NINVEC2) ) + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT1 = ( NBLK_ACT * (NINVEC2) ) + NUM_BLK + 1
C
          CALL WTVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUOUT1,IOFFSET_LUOUT1,
     &                         IOFFSET_INT_LUOUT1,IVCOFF_OUT1,0)
        END IF
C
C
        IOFFSET_LUOUT2     = MY_IOFF_LUOUT2 + 
     &                     ( MY_VEC_IOFF * (NINVEC2) ) + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT2 = ( NBLK_ACT * (NINVEC2) ) + NUM_BLK + 1
CSK           WRITE(6,*) 'This is my OFFSET for writing on ILU3 and what',
CSK     &                 IOFFSET_LUOUT2
CSK           CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
C
        CALL WTVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUOUT2,IOFFSET_LUOUT2,
     &                       IOFFSET_INT_LUOUT2,IVCOFF_OUT2,0)
C
CSK           WRITE(6,*) 'This is what i have written',
CSK     &                 myproc
CSK           CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK = NUM_BLK + LBATCH(ISBATCH)
      END DO
C    /\ loop over batches 2
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ORTHG_VEC_BATCH_PAR3_USE_IA(VEC1,VEC2,FACIN1,FACIN2,
     &                                       FACIN3,NROOT,JOFF,
     &                                       NBATCH,LBATCH,LEBATCH,
     &                                       I1BATCH,IBATCH,MY_VEC_IOFF,
     &                                       MY_IOFF_LUIN1,
     &                                       MY_IOFF_LUIN2,
     &                                       MY_IOFF_LUIN3,
     &                                       MY_IOFF_LUOUT,
     &                                       NBLK_ACT,LUIN1,LUIN2,
     &                                       LUIN3,LUOUT,LUWRTOUT,IADD,
     &                                       MY_INT_IOFF_LUIN1,
     &                                       IVCOFF_IN1,IVCOFF_IN2,
     &                                       IVCOFF_IN3,IVCOFF_OUT)
C
C     calculate vector products between vector on file LUIN1 and LUIN2 
C     resp. LUIN1 and LUIN3. Store vector VEC2 on LUOUT 
C                            ==> LUIN1(VEC2) --> LUOUT(VEC2)
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: JOFF = JVEC - 1
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION FACIN1(*), FACIN3(1)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION IVCOFF_IN1(*), IVCOFF_IN2(*), IVCOFF_IN3(*)
      DIMENSION IVCOFF_OUT(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_IN_LUIN3  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_IN_INT_LUIN1  = 0
      IOFFSET_IN_INT_LUIN2  = 0
      IOFFSET_IN_INT_LUIN3  = 0
      IOFFSET_INT_LUOUT     = 0
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       special treatment for LUIN1 since we need to skip some vectors
C       for simplicity this is done on the outside
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH
        IOFFSET_IN_INT_LUIN1 = MY_INT_IOFF_LUIN1 + NUM_BLK
CSK        IF(MYPROC.EQ.MASTER) THEN
CSK         WRITE(6,*) 'This is my OFFSET for LUIN1 special',
CSK     &               IOFFSET_IN_LUIN1, myproc
CSK        END IF
C
        CALL RDVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN1,IOFFSET_IN_LUIN1,
     &                       IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
C
        DO 100 IROOT = 1, NROOT
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IROOT - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN2 = NUM_BLK + 1 + 
     &                     ( IROOT - 1 ) * NBLK_ACT
C
C         read in batch ISBATCH from LUIN1 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN2',IOFFSET_IN_LUIN2
CSK          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN2,IOFFSET_IN_LUIN2,
     &                         IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
C
C         scaling of VEC2
          IF( IROOT .eq. 1 ) THEN 
            CALL BLK_DSCAL(FACIN2,VEC2,LBATCH(ISBATCH),
     &                     IBATCH(1,I1BATCH(ISBATCH)),
     &                     IOFFSET_IN_INT_LUIN1,IVCOFF_IN1)
          END IF
C         VEC2 == VEC2 + VEC1 * FACIN1(IROOT)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY for LUIN2',
CSK     &                  Myproc,IVEC
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
          CALL SCLVCSUM(FACIN1(IROOT),VEC1,VEC2,LBATCH(ISBATCH),
     &                  IBATCH(1,I1BATCH(ISBATCH)),
     &                  IOFFSET_IN_INT_LUIN2,IVCOFF_IN2)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK          IF (IVEC .eq. NINVEC1 ) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 after DAXPY == after first ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
CSK          END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 KVEC = 1, JOFF
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( KVEC - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN3 = NUM_BLK + 1 +
     &                     ( KVEC - 1 ) * NBLK_ACT
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
CSK          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3,
     &                         IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IROOT)
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 before DAXPY for LUIN3',
CSK     &                  Myproc,IVEC
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
          IF ( KVEC .eq. 1 ) THEN
            CALL SCLVCSUM(FACIN3(1),VEC1,VEC2,LBATCH(ISBATCH),
     &                    IBATCH(1,I1BATCH(ISBATCH)),
     &                    IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
          ELSE 
            CALL SCLVCSUM(FACIN3(KVEC),VEC1,VEC2,LBATCH(ISBATCH),
     &                    IBATCH(1,I1BATCH(ISBATCH)),
     &                    IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
          END IF
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK          IF (KVEC .eq. JOFF ) THEN
CSK            WRITE(6,*) 'VEC1 and VEC2 after DAXPY == after sec ort'
CSK            CALL WRTMAT(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK            CALL WRTMAT(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK          END IF
CSK          END IF
C
  200   CONTINUE
C
C
C       write vector VEC2 from batch IBATCH to file LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + ( MY_VEC_IOFF * JOFF ) + 
     &                      IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = 1 + ( NBLK_ACT * JOFF ) + NUM_BLK
C
        CALL WTVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                      IBATCH(1,I1BATCH(ISBATCH)),
     &                      LUOUT,IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                      IVCOFF_OUT,0)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK = NUM_BLK + LBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CALC_SUBSPACE_H_USE_IA(VEC1,VEC2,SCRRED,IMUSTRED,IVEC,
     &                                  NROOT,ISTRED,
     &                                  NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                                  IBATCH,MY_VEC_IOFF,
     &                                  MY_IOFF_LUIN1,
     &                                  MY_IOFF_LUIN2,MY_IOFF_LUIN3,
     &                                  MY_IOFF_LUIN4,
     &                                  NBLK_ACT,LUIN1,LUIN2,LUIN3,
     &                                  LUIN4,MY_INT_IOFF_LUIN1,
     &                                  MY_INT_IOFF_LUIN2,IVCOFF_IN1,
     &                                  IVCOFF_IN2,IVCOFF_IN3,
     &                                  IVCOFF_IN4,IVCOFF_SPECIAL)
C
C     calculate dot products between vectors on file LUIN1 / LUIN2 
C     and LUIN3 resp. LUIN4
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: IVEC = IVEC
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION SCRRED(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION IVCOFF_IN1(*), IVCOFF_IN2(*), IVCOFF_IN3(*) 
      DIMENSION IVCOFF_IN4(*), IVCOFF_SPECIAL(*)
      REAL*8 BLK_DDOT
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN4
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_SPECIAL
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_SPECIAL
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3, IOFFSET_IN_LUIN4
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH     = 0
      IOFFSET_IN_SPECIAL  = 0
      IOFFSET_IN_LUIN3    = 0
      IOFFSET_IN_LUIN4    = 0
      IOFFSET_IN_INT_SPECIAL  = 0
      IOFFSET_IN_INT_LUIN3    = 0
      IOFFSET_IN_INT_LUIN4    = 0
      LUIN_SPECIAL = 0
      MY_IOFF_SPECIAL     = 0
      MY_INT_IOFF_SPECIAL = 0
C     set correct input file and offsets
      IF( IVEC .le. NROOT ) THEN
        LUIN_SPECIAL = LUIN1
        MY_IOFF_SPECIAL  = MY_IOFF_LUIN1
        MY_INT_IOFF_SPECIAL = MY_INT_IOFF_LUIN1
        CALL ICOPVE(IVCOFF_IN1(MY_INT_IOFF_LUIN1),
     &              IVCOFF_SPECIAL,NBLK_ACT)
      ELSE
        LUIN_SPECIAL = LUIN2
        MY_IOFF_SPECIAL  = MY_IOFF_LUIN2
        MY_INT_IOFF_SPECIAL = MY_INT_IOFF_LUIN2
        CALL ICOPVE(IVCOFF_IN2(MY_INT_IOFF_LUIN2),
     &              IVCOFF_SPECIAL,NBLK_ACT)
      END IF
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       special treatment for LUIN1 / LUIN2 since we need to skip some vectors
C       for simplicity this is done on the outside
C
        IOFFSET_IN_SPECIAL     = MY_IOFF_SPECIAL + IOFFSET_SCRATCH
        IOFFSET_IN_INT_SPECIAL = MY_INT_IOFF_SPECIAL + NUM_BLK
CSK        IF(MYPROC.EQ.MASTER) THEN
CSK         WRITE(6,*) 'This is my OFFSET for LUIN1 / LUIN2 special',
CSK     &               IOFFSET_IN_SPECIAL, LUIN_SPECIAL, myproc
CSK        END IF
C
        CALL RDVEC_BATCH_DRV2(VEC2,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       LUIN_SPECIAL,IOFFSET_IN_SPECIAL,
     &                       IOFFSET_IN_INT_SPECIAL,IVCOFF_SPECIAL)
C
        DO 100 JVEC = 1, MIN(IVEC,NROOT)
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( JVEC - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN3 = 1 + NUM_BLK +
     &                     ( JVEC - 1 ) * NBLK_ACT
C
C         read in batch ISBATCH from LUIN3 to VEC1
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN3',IOFFSET_IN_LUIN3
CSK          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN3,IOFFSET_IN_LUIN3,
     &                         IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
          IJ = IVEC * ( IVEC - 1 ) / 2 + JVEC
C         SCRRED(IJ) = SCRRED(IJ) + SUM ( VEC1 * VEC2 )
          SCRRED(IJ) = SCRRED(IJ) + 
     &                 BLK_DDOT(VEC1,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_INT_LUIN3,IVCOFF_IN3)
C
          IF (ISBATCH .eq. 1 ) THEN
            IMUSTRED = IMUSTRED + 1
            IF( IVEC .eq. 1 .and. JVEC .eq. 1 ) ISTRED = IJ
          END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN4
C
        ILOOP = 0
C
        DO 200 JVEC = ( NROOT + 1), IVEC
C
          ILOOP = ILOOP + 1
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN4 = MY_IOFF_LUIN4 + IOFFSET_SCRATCH +
     &                     ( ILOOP - 1 ) * MY_VEC_IOFF
          IOFFSET_IN_INT_LUIN4 = 1 + NUM_BLK +
     &                     ( ILOOP - 1 ) * NBLK_ACT
C
CSK          IF(MYPROC.EQ.MASTER) THEN
CSK           WRITE(6,*) 'This is my OFFSET for LUIN4',IOFFSET_IN_LUIN4
CSK          END IF
C
          CALL RDVEC_BATCH_DRV2(VEC1,LBATCH(ISBATCH),LEBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         LUIN4,IOFFSET_IN_LUIN4,
     &                         IOFFSET_IN_INT_LUIN4,IVCOFF_IN4)
C
          IJ = IVEC * ( IVEC - 1 ) / 2 + JVEC
C         SCRRED(IJ) = SCRRED(IJ) + SUM ( VEC1 * VEC2 )
          SCRRED(IJ) = SCRRED(IJ) +
     &                 BLK_DDOT(VEC1,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_INT_LUIN4,IVCOFF_IN4)
C
          IF (ISBATCH .eq. 1 ) IMUSTRED = IMUSTRED + 1
C
  200   CONTINUE
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK = NUM_BLK + LBATCH(ISBATCH)
      END DO
C    /\ loop over batches
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PAR(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
     &                       VEC1,VEC2,NVEC,IREW,LBLK)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data 
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1)
*
      IF(IREW .NE. 0 ) CALL REWINE(LUIN,LBLK)
*
      LLUOUT = LUOUT
      LLUSCR = LUSCR
      DO IVEC = 1, NVEC
        CALL REWINE(LLUSCR,LBLK)
        CALL REWINE(LLUOUT,LBLK)
CSK        WRITE(6,*) 'FAC(IVEC) = ',FAC(IVEC)
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,
     &                  LUOUT,LUSCR,0,LBLK)
          ELSE
            CALL SCLVCD(LUIN,LUSCR,FAC(1),VEC1,0,LBLK)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &                0,LBLK)
        END IF
*
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
      END DO
*
      IF(LLUOUT.NE.LUOUT) THEN
C            COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
        CALL COPVCD(LLUOUT,LUOUT,VEC1,1,LBLK)
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV_PAR1(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,NDIM)
C
C     Partition a CI vector into batches of blocks. The length of a
C     batch must be atmost MXLNG.
C     S-vector routine.
C
C     IF ICOMP. eq. 1: the complete CI vector is constructed in just one
C     batch.
C
C     OUTPUT
C     ======
C
C     NBATCH      : Number of batches
C     LBATCH(*)   : Number of blocks in a given batch
C     LEBATCH(*)  : Number of elements in a given batch ( packed ) !
C     I1BATCH(*)  : Number of first block in a given batch
C     IBATCH(8,*) : TTS blocks in Start of a given TTS block with respect to
C                   start
C     of batch --
C     IBATCH(1,*) : Alpha type
C     IBATCH(2,*) : Beta  type
C     IBATCH(3,*) : Sym of alpha
C     IBATCH(4,*) : Sym of beta
C     IBATCH(5,*) : Offset of block with respect to start of block in
C                   expanded form
C     IBATCH(6,*) : Offset of block with respect to start of block in
C                   packed form
C     IBATCH(7,*) : Length of block, expandend form
C     IBATCH(8,*) : Length of block, packed form
C    
C     original version : Jeppe Olsen     - August 1995
C     parallel adaption: S. Knecht       - March  2007 
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C     Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(NDIM)
C     Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) '    PART_CIV_PAR1   '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(6,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(6,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(6,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
*
*. block zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
C       only s-blocks corresponding to distribution
        IF( IBLOCKD(IBLOCK) .ne. MYPROC ) IBATCH(1,IBLOCK) = 0
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV2'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV2. ' )
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) 'Output from PART_CIV_PAR1'
        WRITE(6,*) '========================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches ', NBATCH
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV_PAR2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,SCLFAC,NDIM,IEXCLBLK)
C
C     Partition a CI vector into batches of blocks. The length of a
C     batch must be atmost MXLNG.
C     C-vector routine.
C
C     IF ICOMP. eq. 1: the complete CI vector is constructed in just one
C     batch.
C
C     OUTPUT
C     ======
C
C     NBATCH      : Number of batches
C     LBATCH(*)   : Number of blocks in a given batch
C     LEBATCH(*)  : Number of elements in a given batch ( packed ) !
C     I1BATCH(*)  : Number of first block in a given batch
C     IBATCH(8,*) : TTS blocks in Start of a given TTS block with respect to
C                   start
C     of batch --
C     IBATCH(1,*) : Alpha type
C     IBATCH(2,*) : Beta  type
C     IBATCH(3,*) : Sym of alpha
C     IBATCH(4,*) : Sym of beta
C     IBATCH(5,*) : Offset of block with respect to start of block in
C                   expanded form
C     IBATCH(6,*) : Offset of block with respect to start of block in
C                   packed form
C     IBATCH(7,*) : Length of block, expandend form
C     IBATCH(8,*) : Length of block, packed form
C    
C     original version : Jeppe Olsen     - August 1995
C     parallel adaption: S. Knecht       - March  2007 
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C     input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION SCLFAC(NDIM)
C     output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
C     scratch
      INTEGER LBLOCK_SAVE, LBLOCKP_SAVE 
C
      LBLOCK_SAVE  = 0
      LBLOCKP_SAVE = 0 
C
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) '    PART_CIV_PAR2   '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(6,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(6,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(6,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
C
C     block zero
C
      ILOOPBLK = 0
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
C     loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
      IF( NBATCH .gt. 1 ) ILOOPBLK = ILOOPBLK - 1 
      
C     loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
C     should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
      ILOOPBLK = ILOOPBLK + 1
C     can this block be included
      IBSM   = ISMOST(ISM)
      NSTA   = NSSOA(ISM,IA)
      NSTB   = NSSOB(IBSM,IB)
      LBLOCK = NSTA*NSTB
      LBLOCK_SAVE = LBLOCK
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
      LBLOCKP_SAVE = LBLOCKP
C     check if we need this block, if not: 
C     LBLOCKP = 0
C     LBLOCK  = 0
CSK      IF(IEXCLBLK .eq. 1) THEN
CSK       WRITE(6,*) 'ILOOPBLK,SCLFAC(ILOOPBLK),IBLOCK',ILOOPBLK,
CSK     &             SCLFAC(ILOOPBLK),IBLOCK
CSK      END IF
      IF( IEXCLBLK .eq. 1 ) THEN
         IF( SCLFAC(ILOOPBLK) .eq. 0.0D0 ) THEN
           LBLOCKP = 0
           LBLOCK  = 0
           GOTO 1000
         END IF
      ELSE
        IF( SCLFAC(IBLOCK+1) .eq. 0.0D0 ) THEN
          LBLOCKP = 0
          LBLOCK  = 0
        END IF
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK,myproc
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
C       keep length information, needed?
        IBATCH(7,IBLOCK) = LBLOCK_SAVE
        IBATCH(8,IBLOCK) = LBLOCKP_SAVE
C       all blocks are included, but only active blocks have a length
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV_PAR2'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV_PAR2. ' )
      ELSE
C       This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
C
      IF(NTEST.NE.0) THEN
       IF(IEXCLBLK .eq. 1 ) THEN
        WRITE(6,*) 'Output from PART_CIV_PAR2'
        WRITE(6,*) '========================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches (MYPROC)', NBATCH, MYPROC
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Length of batch           ', LEBATCH(JBATCH)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
       END IF
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV_PAR3(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,ISCLFAC,ISTART_BLK,
     &                         IEXCLBLK,IDEBUGPRNT)
C
C     Partition a CI vector into batches of blocks. The length of a
C     batch must be atmost MXLNG.
C     C-vector routine.
C
C     IF ICOMP. eq. 1: the complete CI vector is constructed in just one
C     batch.
C
C     OUTPUT
C     ======
C
C     NBATCH      : Number of batches
C     LBATCH(*)   : Number of blocks in a given batch
C     LEBATCH(*)  : Number of elements in a given batch ( packed ) !
C     I1BATCH(*)  : Number of first block in a given batch
C     IBATCH(8,*) : TTS blocks in Start of a given TTS block with respect to
C                   start
C     of batch --
C     IBATCH(1,*) : Alpha type
C     IBATCH(2,*) : Beta  type
C     IBATCH(3,*) : Sym of alpha
C     IBATCH(4,*) : Sym of beta
C     IBATCH(5,*) : Offset of block with respect to start of block in
C                   expanded form
C     IBATCH(6,*) : Offset of block with respect to start of block in
C                   packed form
C     IBATCH(7,*) : Length of block, expandend form
C     IBATCH(8,*) : Length of block, packed form
C    
C     original version : Jeppe Olsen     - August 1995
C     parallel adaption: S. Knecht       - May    2007 
C
C     Last revision:     S. Knecht       - May    2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C     input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION ISCLFAC(*)
C     output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
C     scratch
      INTEGER LBLOCK_SAVE, LBLOCKP_SAVE 
C
      I_SAVE_IBLOCK = 0
      N_TRY = 0
C
      LBLOCK_SAVE  = 0
      LBLOCKP_SAVE = 0 
C
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ==================='
        WRITE(LUWRT,*) '    PART_CIV_PAR3   '
        WRITE(LUWRT,*) ' ==================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(LUWRT,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(LUWRT,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(LUWRT,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
C
C     block zero
C
      ILOOPBLK = 0
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
C     loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
      IF( NBATCH .gt. 1 ) ILOOPBLK = ILOOPBLK - 1 
C     loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
C     should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
      ILOOPBLK = ILOOPBLK + 1
C     can this block be included
      IBSM   = ISMOST(ISM)
      NSTA   = NSSOA(ISM,IA)
      NSTB   = NSSOB(IBSM,IB)
      LBLOCK = NSTA*NSTB
      LBLOCK_SAVE = LBLOCK
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
      LBLOCKP_SAVE = LBLOCKP
C     check if we need this block, if not: 
C     LBLOCKP = 0
C     LBLOCK  = 0
CSK      IF(IEXCLBLK .eq. 1) THEN
CSK       WRITE(6,*) 'ILOOPBLK,SCLFAC(ILOOPBLK),IBLOCK',ILOOPBLK,
CSK     &             SCLFAC(ILOOPBLK),IBLOCK
CSK      END IF
*
*     IEXCLBLK == 1: DEBUG ME!!!
*
      IF( IEXCLBLK .eq. 1 ) THEN
         IF( ISCLFAC(ILOOPBLK) .eq. 0 ) THEN
           LBLOCKP = 0
           LBLOCK  = 0
           GOTO 1000
         END IF
      ELSE
        IF( ISCLFAC(ISTART_BLK + IBLOCK + 1 ) .eq. 0 ) THEN
          LBLOCKP = 0
          LBLOCK  = 0
        END IF
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK,myproc
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
C       keep length information on 7!
        IBATCH(7,IBLOCK) = LBLOCK_SAVE
        IBATCH(8,IBLOCK) = LBLOCKP
C       all blocks are included, but only active blocks have a length
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                   single block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV_PAR3'
        write(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV_PAR3. ' )
      ELSE
*
*     this batch is finished, goto next batch
*
*       is this the TT-1 block?
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(LUWRT,*) ' Not enough scratch space to include a '//
     &                     ' single Block'
            WRITE(LUWRT,*) ' Since I cannot procede I will stop '
            WRITE(LUWRT,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV_PAR3'
            WRITE(LUWRT,*) '  LENGTH ',LENGTH
            WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from',MXLNG
            CALL RMEMCHK
            stop ' In PART_CIV_PAR3 because of N_TRY .ge. 1. '
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,LUWRT)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
C
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV_PAR3'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of batches (MYPROC)', NBATCH, MYPROC
        DO JBATCH = 1, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Length of batch           ', 
     &                          LEBATCH(JBATCH)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lenghts of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV4_PAR(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         IOFFBLK,IOFFBTC,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,IDEBUGPRNT,
     &                         ILENB_OLD,INUMB_OLD,IRUN_THIS)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* Stefan Knecht, May 2007
*
* further modified from PART_CIV3 for parallel purposes
*
* Compared to PART_CIV3 the changes are
*
* 1 : partitioning with respect to block <--> node distribution
*
* 2 : to be continued ...
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "../luciarel/ipoist8.inc"
#include "../luciarel/lucipar_rel.inc"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) '     PART_CIV4_PAR    '
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      IFINI = 0
      IF( IRUN_THIS .eq. 1 ) THEN
        NBATCH = IOFFBTC - 1
CSK        WRITE(LUWRT,*) 'NBATCH is ',NBATCH
        IBLOCK = IOFFBLK - 1
      ELSE
        NBATCH = IOFFBTC 
C       be careful - is that correct???
        NBATCH = IOFFBTC
CSK        WRITE(LUWRT,*) 'NBATCH is now',NBATCH
C
        IBLOCK = IOFFBLK - 1
        LBATCH(  NBATCH ) = INUMB_OLD
        LEBATCH( NBATCH ) = ILENB_OLD
        LENGTH  = ILENB_OLD  
        LENGTHP = ILENB_OLD
        GOTO 1000
      END IF
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
C
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                  single Block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV4'
        WRITE(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        stop ' In PART_CIV4_PAR. '
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      ILENB_OLD = LEBATCH(NBATCH)
      INUMB_OLD = LBATCH(NBATCH)
*
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV4_PAR'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lenghts of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PRINT_BATCH_INFO(NCBATCH,LBATC,LEBATC,
     &           I1BATC,IBATC,ICCTOS,ICWEIGHT,NBLOCK,IBTOTW,ICWEIGHTF,
     &           IABSOLUTE_WEIGHT)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION LBATC(*),LEBATC(*),I1BATC(*),IBATC(8,*)
      DIMENSION ICCTOS(NBLOCK,NBLOCK),ICWEIGHT(NBLOCK)
      DIMENSION IBTOTW(NBLOCK),ICWEIGHTF(NBLOCK)
      INTEGER*8 IABSOLUTE_WEIGHT
*
      NTEST = 99
      IF(MYPROC.EQ.MASTER .and. NTEST .ge. 100 )THEN
        WRITE(6,*) ' OUTPUT from PRINT_BATCH_INFO'
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' '
        WRITE(6,*) ' Number of batches ', NCBATCH
      END IF
      DO JBATCH = 1, NCBATCH
        IF(MYPROC.EQ.MASTER .and. NTEST .ge. 100 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*) ' '
          WRITE(6,*) ' Number of blocks included ', LBATC(JBATCH)
        END IF
CSK        WRITE(6,*) '   TTSS and offsets and lenghts of each block '
        DO IBLOCK =I1BATC(JBATCH),I1BATC(JBATCH)+ LBATC(JBATCH)-1
CSK          WRITE(6,'(10X,4I3,4I8)') (IBATC(II,IBLOCK),II=1,8)
CSK          IF(IBATC(8,IBLOCK).NE.0.AND.MYPROC.EQ.MASTER.AND.NTEST.GE.10) 
CSK     &       WRITE(6,'(2X,A,1X,I5,1X,A,1X,I8)') 'Block number',IBLOCK,
CSK     & 'has length',(IBATC(II,IBLOCK),II=8,8)
          JJJ = 0
          DO IJ = 1, NBLOCK
            IF(ICCTOS(IJ,IBLOCK).NE.0) THEN
              JJJ = JJJ + 1
              ICWEIGHT(JJJ) = IJ
            ENDIF
          END DO
          ICWEIGHTF(IBLOCK) = JJJ
          IF(JJJ.GT.0.AND.MYPROC.EQ.MASTER.AND.NTEST.GE.100) THEN
            WRITE(6,'(2X,A,I4)')'Number of connections for this block:',
     &                           JJJ
            CALL IWRTMA(ICWEIGHT,1,JJJ,1,JJJ)
            
            WRITE(6,*)'  '
          ENDIF
          IBTOTW(IBLOCK) = (JJJ) * IBATC(8,IBLOCK)
        END DO
      END DO
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10) THEN
CSK        WRITE(6,'(2X,A)')'TOTAL WEIGHT FOR EACH BLOCK'
CSK        CALL IWRTMA8(IBTOTW,1,NBLOCK,1,NBLOCK)
CSK        WRITE(6,'(2X,A)')'WEIGHT MATRIX'
CSK        CALL IWRTMA(ICWEIGHTF,1,NBLOCK,1,NBLOCK)
      END IF
      DO IBLK = 1, NBLOCK
        IABSOLUTE_WEIGHT = IABSOLUTE_WEIGHT + 
     &                     ( ICWEIGHTF(IBLK) * IBATC(8,IBLK) )
      END DO
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PRSYMMN(A,MATDIM,LUWRT)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
#include "implicit.h"
      DIMENSION A(*)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(LUWRT,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,5(1X,E13.7),/,(1H ,5X,5(1X,E13.7)))
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE REDVEC(VECIN,VECOUT,NDIM,NTYPE,NOP,NCOMM,NTARGET)
*
*. Reduce an array VECIN of dimension NDIM that is distributed over 
*. NCOMM to VECOUT on either a single process (NTARGET >= 0) or on 
*. all processes (NTARGET < 0). NTYPE determines the MPI_DATATYPE 
*. used in the reduce-call. NOP is the reduce-operation.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VECIN(*),VECOUT(*)
      INTEGER DATATYPE
*
      IF(NTYPE .EQ. 1) THEN
        DATATYPE = my_MPI_INTEGER
      ELSE IF (NTYPE .EQ. 2) THEN
C MPI_REAL8 is equivalent to MPI_DOUBLE_PRECISION
C the old LAM-MPI library does not support MPI_ALLREDUCE for MPI_REAL8
C so we use instead MPI_DOUBLE_PRECISION
CSK        DATATYPE = MPI_REAL8
        DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (NTYPE .EQ. 3)THEN
        DATATYPE = MPI_REAL
      ELSE IF (NTYPE .EQ. 4 .OR. NTYPE .EQ. 5) THEN
        DATATYPE  =  MPI_CHARACTER 
      ELSE
        WRITE(6,*)' REDVEC: TYPE "',NTYPE,'" does not exist!'
        CALL Abend2('ERROR in REDVEC: Nonexisting DATATYPE')
      END IF
*
      IF(NTARGET.GE.0)THEN
        CALL MPI_REDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,NTARGET,
     &                  NCOMM,IERR)
      ELSE IF(NTARGET.LT.0)THEN
        CALL MPI_ALLREDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,
     &                     NCOMM,IERR)
      ENDIF
*. return since everything is reduced now!
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE REDVEC_REL(VECIN,VECOUT,NDIM,NTYPE,NOP,NCOMM,NTARGET)
*
*. Reduce an array VECIN of dimension NDIM that is distributed over 
*. NCOMM to VECOUT on either a single process (NTARGET >= 0) or on 
*. all processes (NTARGET < 0). NTYPE determines the MPI_DATATYPE 
*. used in the reduce-call. NOP is the reduce-operation.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      DIMENSION VECIN(*),VECOUT(*)
      INTEGER DATATYPE
*
      IF(NTYPE .EQ. 1) THEN
        DATATYPE = my_MPI_INTEGER
      ELSE IF (NTYPE .EQ. 2) THEN
C MPI_REAL8 is equivalent to MPI_DOUBLE_PRECISION
C the old LAM-MPI library does not support MPI_ALLREDUCE for MPI_REAL8
C so we use instead MPI_DOUBLE_PRECISION
CSK        DATATYPE = MPI_REAL8
        DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (NTYPE .EQ. 3)THEN
        DATATYPE = MPI_REAL
      ELSE IF (NTYPE .EQ. 4 .OR. NTYPE .EQ. 5) THEN
        DATATYPE  =  MPI_CHARACTER 
      ELSE
        WRITE(6,*)' REDVEC: TYPE "',NTYPE,'" does not exist!'
        CALL Abend2('ERROR in REDVEC: Nonexisting DATATYPE')
      END IF
*
      IF(NTARGET.GE.0)THEN
        CALL MPI_REDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,NTARGET,
     &                  NCOMM,IERR)
      ELSE IF(NTARGET.LT.0)THEN
        CALL MPI_ALLREDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,
     &                     NCOMM,IERR)
      ENDIF
*. return since everything is reduced now!
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLVCSUM(FAC,VEC1,VEC2,NBATCH_BLK,NBATCH_INF,
     &                    IOFF_INT,IVCOFF)
 
C
C     OUTPUT
C     ======
C
C     VEC2 == VEC1 * FAC + VEC2
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*), VEC2(*),IVCOFF(*),NBATCH_INF(8,*)
      INTEGER JOFF,ILEN
C
      DO ISBLK = 1, NBATCH_BLK
C
        ILEN = IVCOFF( IOFF_INT + ISBLK - 1 )
C
        IF ( ILEN .gt. 0 ) THEN
          JOFF = NBATCH_INF(6,ISBLK)
          CALL DAXPY( ILEN, FAC, VEC1(JOFF), 1, VEC2(JOFF), 1 )
        END IF

      END DO
C     
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION BLK_DDOT(VEC1,VEC2,NBATCH_BLK,NBATCH_INF,
     &                         IOFF_INT,IVCOFF)
 
C
C     OUTPUT
C     ======
C
C     BLK_DDOT == SUM( VEC1 * VEC2)
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC1(*), VEC2(*),IVCOFF(*),NBATCH_INF(8,*)
      INTEGER JOFF,ILEN
C
      XXX = 0.0D0
C
      DO ISBLK = 1, NBATCH_BLK
C
        ILEN = IVCOFF( IOFF_INT + ISBLK - 1 )
C
        IF ( ILEN .gt. 0 ) THEN
          JOFF = NBATCH_INF(6,ISBLK)
          XXX = XXX + DDOT( ILEN, VEC1(JOFF), 1, VEC2(JOFF), 1 )
        END IF

      END DO
C     
      BLK_DDOT = XXX
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SET_ACTIVE(VEC,IAMACTIVE,FAC,NDIM)
C
C     OUTPUT
C     ======
C
C     active blocks in VEC flagged with a factor FAC
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION VEC(NDIM), IAMACTIVE(NDIM)
C
      DO IBLK = 1, NDIM
        IF( IAMACTIVE(IBLK) .eq. MYPROC ) VEC(IBLK) = FAC
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SET_GROUP_TABLE(IGROUPLIST,IPROCLIST,ICOLOR)
C
C     fill in the array IGROUPLIST. At the end it should contain all 
C     processes numbers belonging to the same group
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C     INPUT
      DIMENSION IPROCLIST(NMPROC)
C     OUTPUT
      DIMENSION IGROUPLIST(NMPROC)
C
      INUMB = 1
      DO IPROC = 1, NMPROC
C
        IF( IPROCLIST(IPROC) .eq. ICOLOR ) THEN
C          put in the process tag ( master is 0 )
           IGROUPLIST(INUMB) = IPROC - 1
           INUMB = INUMB + 1
        END IF
C
      END DO
C
      IF( INUMB-1 .gt. NEWCOMM_PROC ) THEN
        WRITE(6,*) 'Error in SET_GROUP_TABLE: more processes assigned 
     & than in this group included!'
        WRITE(6,*) 'assigned processes, group size',INUMB-1,NEWCOMM_PROC
        CALL Abend2('Error detected in gp/gplupar.F: SET_GROUP_TABLE')
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SET_GROUP_TABLE_REL(IGROUPLIST,IPROCLIST,ICOLOR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C
C     INPUT
      DIMENSION IPROCLIST(NMPROC)
C     OUTPUT
      DIMENSION IGROUPLIST(NMPROC)
C
      INUMB = 1
      DO IPROC = 1, NMPROC
C
        IF( IPROCLIST(IPROC) .eq. ICOLOR ) THEN
C         put in the process tag ( master is 0 )
          IGROUPLIST(INUMB) = IPROC - 1
          INUMB = INUMB + 1
        END IF
C
      END DO
C
      IF( INUMB-1 .gt. NEWCOMM_PROC ) THEN
        WRITE(LUWRT,*) 'Error in SET_GROUP_TABLE_REL: more CPUs assigned
     & as included in this group!'
        WRITE(LUWRT,*) 'assigned CPUs, group size',INUMB-1,NEWCOMM_PROC
        CALL Abend2('Error detected in gp/gplupar.F: 
     &               SET_GROUP_TABLE_REL')
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_CLOSE
C
C     OUTPUT
C     ======
C
C     close files
C     files handles stored on common block LUCIPAR are set to MPI_NULL
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
C
C     ... close
C
      CALL MPI_FILE_CLOSE(ILU1,IERR)
      CALL MPI_FILE_CLOSE(ILU2,IERR)
      CALL MPI_FILE_CLOSE(ILU3,IERR)
      CALL MPI_FILE_CLOSE(ILU4,IERR)
      CALL MPI_FILE_CLOSE(ILU5,IERR)
      CALL MPI_FILE_CLOSE(ILU7,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_CLOSE_REL
C
C     OUTPUT
C     ======
C
C     close files
C     files handles stored on common block LUCIPFILE are set to MPI_NULL
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
C
C     ... close
C
      CALL MPI_FILE_CLOSE(IDIA,IERR)
      CALL MPI_FILE_CLOSE(ILU1,IERR)
      CALL MPI_FILE_CLOSE(ILU2,IERR)
      CALL MPI_FILE_CLOSE(ILU3,IERR)
      CALL MPI_FILE_CLOSE(ILU4,IERR)
      CALL MPI_FILE_CLOSE(ILU5,IERR)
      CALL MPI_FILE_CLOSE(ILU7,IERR)
      CALL MPI_FILE_CLOSE(ILUC,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_OPEN_REL(NFILE_ID)
C
C     OUTPUT
C     ======
C
C     open files ready for MPI-I/O
C     files handles stored on common block LUCIPFILE
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../luciarel/lucipar_rel.inc"
      CHARACTER*6 PLU1BL,PLU2BL,PLU3BL,PLU4BL,PLU7BL,PLU5BL, PLUCBL
      CHARACTER*10 PDIA, PLU1, PLU2, PLU3, PLU4, PLU5, PLU6, PLU7, PLUC
      CHARACTER*6 PDIABL 
      INTEGER PFILELENGTH
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP
      IDISP = 0
      PFILELENGTH = 0
C
      PDIABL = 'DIAPAR'
      PLU1BL = 'LU1PAR'
      PLU2BL = 'LU2PAR'
      PLU3BL = 'LU3PAR'
      PLU4BL = 'LU4PAR'
      PLU5BL = 'LU5PAR'
      PLU7BL = 'LU7PAR'
      PLUCBL = 'LUCPAR'
       
C     set correct names
      IF( NFILE_ID .lt. 10 ) THEN    ! MPI ID has one digit
        WRITE (PDIA,'(A6,A1,I1)') PDIABL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I1)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I1)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I1)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I1)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I1)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I1)') PLU7BL,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I1)') PLUCBL,'.',NFILE_ID
        PFILELENGTH = 8
      ELSE IF( NFILE_ID .lt. 100)THEN  ! MPI ID has two digits
        WRITE (PDIA,'(A6,A1,I2)') PDIABL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I2)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I2)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I2)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I2)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I2)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I2)') PLU7BL,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I2)') PLUCBL,'.',NFILE_ID
        PFILELENGTH = 9
      ELSE                        ! MPI ID has three digits
        WRITE (PDIA,'(A6,A1,I3)') PDIABL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I3)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I3)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I3)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I3)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I3)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I3)') PLU7BL,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I3)') PLUCBL,'.',NFILE_ID
        PFILELENGTH = 10
      END IF
C
C     ... open
C
      CALL MPI_FILE_OPEN(MYNEW_COMM,PDIA(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,IDIA,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU1(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU1,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU2(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU2,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU3(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU3,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU4(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU4,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU5(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU5,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU7(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU7,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLUC(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILUC,IERR)
C
C     ... set fileview
C
      CALL MPI_FILE_SET_VIEW(IDIA,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU1,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU2,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU3,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU4,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU5,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU7,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILUC,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_OPEN(NFILE_ID)
C
C     OUTPUT
C     ======
C
C     open files ready for MPI-I/O
C     files handles stored on common block LUCIPAR
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      CHARACTER*6 PLU1BL,PLU2BL,PLU3BL,PLU4BL,PLU7BL,PLU5BL 
      CHARACTER*7 PLU1IBL,PLU2IBL,PLU3IBL,PLU4IBL,PLU7IBL,PLU5IBL 
      INTEGER PFILELENGTH
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP
      IDISP = 0
      PFILELENGTH = 0
C
      PLU1BL = 'LU1PAR'
      PLU2BL = 'LU2PAR'
      PLU3BL = 'LU3PAR'
      PLU4BL = 'LU4PAR'
      PLU5BL = 'LU5PAR'
      PLU7BL = 'LU7PAR'
       
C     set correct names
      IF( NFILE_ID .lt. 10 ) THEN    ! MPI ID has one digit
        WRITE (PLU1,'(A6,A1,I1)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I1)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I1)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I1)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I1)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I1)') PLU7BL,'.',NFILE_ID
        PFILELENGTH = 8
      ELSE IF( NFILE_ID .lt. 100)THEN  ! MPI ID has two digits
        WRITE (PLU1,'(A6,A1,I2)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I2)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I2)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I2)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I2)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I2)') PLU7BL,'.',NFILE_ID
        PFILELENGTH = 9
      ELSE                        ! MPI ID has three digits
        WRITE (PLU1,'(A6,A1,I3)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I3)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I3)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I3)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I3)') PLU5BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I3)') PLU7BL,'.',NFILE_ID
        PFILELENGTH = 10
      END IF
C
C     ... open
C
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU1(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU1,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU2(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU2,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU3(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU3,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU4(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU4,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU5(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU5,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU7(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU7,IERR)
C
C     ... set fileview
C
      CALL MPI_FILE_SET_VIEW(ILU1,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU2,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU3,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU4,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU5,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU7,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE TRAVCD_PAR(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
     &                      ICOPY,LBLK,LUSCR1,LUSCR2)
*
* NVECIN vectors reside on LU1, Transform these vectors,
*  using LUSCR1 and  LUSCR2  as
* scratch files,  with matrix X to produce output file LUOUT.
*
* Since LUIN is accessed several times it is always
* REWOUND. LUOUT is written to from current start.
*
* I ICOPY .ne. 0 the transformed vectors are written back to LUIN
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION X(NVECIN,NVECOUT)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*
C             MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
      DO IVECOUT = 1, NVECOUT
         CALL MVCSMD_PAR(LUIN,X(1,IVECOUT),LUSCR1,LUSCR2,VEC1,VEC2,
     &               NVECIN,1,LBLK)
         CALL REWINO(LUSCR1)
         CALL COPVCD(LUSCR1,LUOUT,VEC1,0,LBLK)
      END DO
*
      IF(ICOPY.EQ.1) THEN
        CALL REWINO(LUIN)
        CALL REWINO(LUOUT)
        DO IVECOUT = 1, NVECOUT
          CALL COPVCD(LUOUT,LUIN,VEC1,0,LBLK)
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE WRSVCD_PAR(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,
     &                     IBLOCKL,IBLOCKD,NBLOCK,JPACK)
*
* Write scattered vector to disc, every node writes only it's part.
* Information about the length of a blcok is stored in IBLOCKL.
* Information about which block belongs to which node is stored in
* IBLOCKD. 
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      
*.Input
      DIMENSION IPLAC(*),VAL(*),IBLOCKL(NBLOCK),IBLOCKD(NBLOCK)
*.Scratch
      DIMENSION VEC1(*)
*
      NTEST = 0
      IBOT = 1
*. Loop over blocks
      DO 500 II = 1, NBLOCK
        IF(IBLOCKD(II).NE.MYPROC) GOTO 300
*.Length
        LBL = IBLOCKL(II)
        CALL ITODS (LBL,1,-1,LU)
CSK        write(6,*) ' IBOT, LBL  ',IBOT, LBL
        CALL DZERO(VEC1,LBL)
C       *******************************************************
C       *******************************************************
C              Find the places where a 1.0 should be set 
C       *******************************************************
C       *******************************************************
        DO 200 IEFF = 1, NSCAT
          IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
     &    VEC1(IPLAC(IEFF)-IBOT+1) = VAL(IEFF)
*
          IF(NTEST.GE.10) THEN
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
     &      write(6,*) 'IBOT,IBOT+LBL-1',IBOT,IBOT+LBL-1
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
     &      write(6,*) ' Catch : IPLAC(IEFF) VAL(IEFF) ',
     &      IPLAC(IEFF),VAL(IEFF)
          END IF
*
  200     CONTINUE
          IF(NTEST.GT.10) THEN
            IF(MYPROC.EQ.MASTER) THEN
             write(6,*)  'record to be disced '
             CALL WRTMAT(VEC1,1,LBL,1,LBL)
            END IF
          END IF
          IF(JPACK.EQ.1) THEN
            CALL TODSCP(VEC1,LBL,-1,LU)
          ELSE
            CALL TODSC(VEC1,LBL,-1,LU)
          END IF
 300      CONTINUE
          LBL = IBLOCKL(II)
          IBOT = IBOT + LBL
 500  CONTINUE
*
      IF(NTEST.GT.10)THEN
        IF(MYPROC.EQ.MASTER)THEN
          WRITE(6,*) ' output from WRSVCD'
          CALL WRTVCD(VEC1,LU,1,LBLK)
        ENDIF
      END IF
*. work is done, all is needed now is writing an EOF
      LBL = -1
      CALL ITODS(LBL,1,-1,LU)
*
      END
#else
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_DUMMY2
      END 
***********************************************************************
#endif /* defined (VAR_MPI) */
