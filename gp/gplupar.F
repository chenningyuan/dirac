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
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTBLKND_2(NDIM,ICCTOS,ICWEIGHTF,NPARBLOCK,NBLKWT)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
*. input
* NNJJJBLK: total number of c-blocks
* NTEMPP: number of processes available for the distribution
* ICSBTCHB(NNJJJBLK): total length of each block
*. output
* IDISTLIST(NNJJJBLK): list containing the number of a node  
* that computes the block
*
      INTEGER ICOUNTABLK, NTEST,MWEIGHT
      INTEGER*8 ITOTBLCKL,MMPINFA
      DIMENSION ICCTOS(NDIM,NDIM),ICWEIGHTF(NDIM)
      DIMENSION NPARBLOCK(NDIM),NBLKWT(NDIM),MMPINFA(NMPROC)
      DIMENSION MWEIGHT(NMPROC)
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: number of the current largest block
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2,ILEN
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*
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
*. end of initializtion
      DO ICBL = 1, NDIM
*
        ILEN = NBLKWT(ICBL)
        IF(ILEN.NE.0) THEN
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
          ITEMPL = NBLKWT(II)
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
       MWEIGHT(IAMGPROC) = MWEIGHT(IAMGPROC) + 1
       MMPINFA(IAMGPROC) = MMPINFA(IAMGPROC) + MXSZTMP
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
      WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF MYPROC:',
     &MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MWEIGHT(ISTI),' blocks with a
     & total length of',MMPINFA(ISTI)
      END DO
*
*      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
CSK      END IF
*
CSK      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISTWORK(NDIM,ICCTOS,ICWEIGHTF,NPARBLOCK,NBLKWT,NBATCH,
     &                    LBATS)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION ICCTOS(NDIM,NDIM),ICWEIGHTF(NDIM),NPARBLOCK(NDIM)
      DIMENSION NBLKWT(NDIM),LBATS(NBATCH)
      INTEGER*8 MWEIGHT,NWEIGHT
      DIMENSION MWEIGHT(NMPROC*NBATCH),NWEIGHT(NBATCH)
*
*. find a distribution of work for a given batch as even as possible 
*
       
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
              CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FAC,LENGTH)
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
          XINOUT = XINOUT + INPROD(VEC1,VEC2,LENGTH)
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
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FNDMND_PAR(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &                     SCR,ISCAT,SUBVAL,IBLOCKL,IBLOCKD,
     &                     NBLOCK,NTESTG)
*
* FIND NSUB LOWEST ELEMENTS OF VECTOR RESIDING ON FILE
* LU. ENSURE THAT NO DEGENERENCIES ARE SPLIT
*
*
* INPUT
*=======
* LU :    FILE WHERE VECTOR OF INTEREST IS RESIDING, REWOUND
* LBLK :  DEFINES FILE STRUCTURE ON FILE LU
* NSUBMX: LARGEST ALLOWED NUMBER OF SORTED ELEMENTS
*
* OUTPUT
*=======
* NSUB : ACTUAL NUMBER OF ELEMENTS OBTAINED. CAN BE SMALLER
*        THAN NSUBMX IF THE LAST ELEMENT BELONGS TO A DEGENERATE
*        SET
*ISCAT:  SCATTERING ARRAY, ISCAT(I) GIVES FULL ADRESS OF SORTED
*        ELEMENT I
*SUBVAL: VALUE OF SORTED ELEMENTS

      IMPLICIT REAL*8           ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
#include "../lucita/wrkspc.inc"
      
      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
      DIMENSION IBLOCKL(NBLOCK), IBLOCKD(NBLOCK)
*
      CALL REWINE(LU,-1)
*
*.LOOP OVER BLOCKS
*
C     write(6,*) ' FNDMND NSUBMX = ', NSUBMX
*
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 0
*
      IBASE = 1
      LSUB = 0
*. Loop over blocks
      DO 1000 II = 1, NBLOCK    
*
        IF(IBLOCKD(II).EQ.MYPROC)THEN
          CALL IFRMDS(LBL,1,-1,LU)
        ELSE
*.        useful to set all other blocks to 0?
          LBL = 0
        ENDIF
*
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Info about block ',II
          WRITE(6,*) ' Number of elements ',LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          IF(IBLOCKD(II).EQ.MYPROC)THEN 
            CALL FRMDSC(SEGMNT,LBL ,-1,LU,IMZERO,IAMPACK)
          ENDIF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Elements read in '
            CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
          END IF
          IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
            MSUBMX = MIN(NSUBMX,LBL)
            IF(LBL.GE.1) THEN
              CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                    MSUBMX,MSUB,NTEST)
            ELSE
              MSUB = 0
            END IF
            DO 10 I = 1, MSUB
   10         ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
* SORT COMBINED LIST
            MSUBMX = MIN(NSUBMX,LSUB+MSUB)
            IF(MSUBMX.GT.0) THEN
              CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &                    NTEST)
            ELSE
              LLSUB = 0
            END IF
            LSUB = LLSUB
            DO 20 I = 1, LSUB
              ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20       CONTINUE
*
            CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
            CALL COPVEC(SUBVAL,SCR,LSUB)
*. set to lbl to true value
            LBL = IBLOCKL(II)
*.
            IBASE = IBASE + LBL
*.
            IF(NTEST .GE. 20 ) THEN
              WRITE(6,*) ' Lowest elements and their original place '
              WRITE(6,*) ' Number of elements obtained ', LSUB
              CALL WRTMAT(SUBVAL,1,LSUB,1,LSUB)
              CALL IWRTMA(ISCR,1,LSUB,1,LSUB)
            END IF
          END IF
*
        END IF
*
 1000 CONTINUE
*
      NSUB = LSUB
      CALL ICOPVE(ISCR,ISCAT,NSUB)
      IF(NTEST .GE. 20) THEN
        IF(MYPROC.EQ.MASTER) then
          WRITE(6,*) ' Lowest elements and their original place '
          WRITE(6,*) ' Number of elements obtained ', NSUB
          WRITE(6,*) ' WRITTEN BY MASTER ', MASTER
          CALL WRTMAT(SUBVAL,1,NSUB,1,NSUB)
          CALL IWRTMA(ISCAT,1,NSUB,1,NSUB)
        END IF
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GATHER')
*
      CALL MEMMAN(KGATHERA,NMPROC*NSUBMX,'ADDL  ',2,'PARRA1')
      CALL MEMMAN(KGATHERB,NMPROC*NSUBMX,'ADDL  ',2,'PARRA2')
      CALL MEMMAN(KGATHERC,NMPROC*NSUBMX,'ADDL  ',1,'PARIA1')
      CALL MEMMAN(KGATHERD,NMPROC*NSUBMX,'ADDL  ',1,'PARIA2')
      CALL MEMMAN(KGATHERE,NSUBMX,'ADDL  ',1,'PARIA3')
*. We gather all lowest values from each node 
*. and build up a combined list of those
      CALL GATHER_LOW_PAR(NSUB,NSUBMX,SUBVAL,ISCAT,
     &                    WORK(KGATHERA),WORK(KGATHERB),
     &                    WORK(KGATHERC),WORK(KGATHERD),
     &                    WORK(KGATHERE),MPI_COMM_WORLD,NTESTG)
*. update SCR1 and ISCR1
      CALL COPVEC(SUBVAL,SCR,NSUBMX)
      CALL ICOPVE(ISCAT,ISCR,NSUBMX)
     
CSK      NTEST = 20
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.20)THEN
        WRITE(6,*)'after search: SUBVAL and ISCAT'
        CALL WRTMAT(SUBVAL,1,NSUBMX,1,NSUBMX)
        CALL IWRTMA(ISCAT,1,NSUBMX,1,NSUBMX)
      END IF
CSK      NTEST = 0
           
      IF(NSUB.NE.NSUBMX.AND.MYPROC.EQ.MASTER.AND.NTEST.GE.20)THEN
        WRITE(6,*)'Warning! NSUB is lower than NSUBMX'
        WRITE(6,*)'NSUB is set to be equal to NSUBMX'
        NSUB = NSUBMX
      END IF
*
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'GATHER')
*
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
      SUBROUTINE PART_CIV_PAR(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                        NSMST,MXLNG,IOCOC,ISMOST,
     &                        NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                        ICOMP,ITTSS_ORD,NODEMXLNG,
     &                        ICWEIGHT,IBTOTW,ICWEIGHTF,NBLOCK,REORDER)
*
* Partition a CI vector into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1: the complete CI vector is constructed in just one
* batch
*
* OUTPUT
*
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta  type
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expandend form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) '     PART_CIV_PAR   '
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
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
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
C?      WRITE(6,*) 'Output from PART_CIV'
C?      WRITE(6,*) '====================='
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
      SUBROUTINE PRINT_BATCH_INFO(NCBATCH,LBATC,LEBATC,
     &           I1BATC,IBATC,ICCTOS,ICWEIGHT,NBLOCK,IBTOTW,ICWEIGHTF)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      DIMENSION LBATC(*),LEBATC(*),I1BATC(*),IBATC(8,*)
      DIMENSION ICCTOS(NBLOCK,NBLOCK),ICWEIGHT(NBLOCK)
      DIMENSION IBTOTW(NBLOCK),ICWEIGHTF(NBLOCK)
*
      NTEST = 100
      IF(MYPROC.EQ.MASTER)THEN
        WRITE(6,*) ' OUTPUT from PRINT_BATCH_INFO'
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' '
        WRITE(6,*) ' Number of batches ', NCBATCH
      END IF
      DO JBATCH = 1, NCBATCH
        IF(MYPROC.EQ.MASTER) THEN
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*) ' '
          WRITE(6,*) ' Number of blocks included ', LBATC(JBATCH)
        END IF
CSK        WRITE(6,*) '   TTSS and offsets and lenghts of each block '
        DO IBLOCK =I1BATC(JBATCH),I1BATC(JBATCH)+ LBATC(JBATCH)-1
CSK          WRITE(6,'(10X,4I3,4I8)') (IBATC(II,IBLOCK),II=1,8)
          IF(IBATC(8,IBLOCK).NE.0.AND.MYPROC.EQ.MASTER.AND.NTEST.GE.10) 
     &       WRITE(6,'(2X,A,1X,I5,1X,A,1X,I8)') 'Block number',IBLOCK,
     & 'has length',(IBATC(II,IBLOCK),II=8,8)
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
        WRITE(6,'(2X,A)')'TOTAL WEIGHT FOR EACH BLOCK'
        CALL IWRTMA8(IBTOTW,1,NBLOCK,1,NBLOCK)
      END IF
*
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
        CALL SETVEC(VEC1,0.0D0,LBL)
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
*. work is done, all is needed now is writing and EOF
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
