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
      MXSIZE = 3900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 3900000000
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
      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.MASTER) THEN
          WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF 
     &    MYPROC:',MYPROC
          DO ISTI = 1, IMINNP
            WRITE(6,'(3X,A,I17,A,I17,A,1X,I17)')
     &      'process',ISTI -1,' calculates',NPARBLKWT(1,ISTI),
     &      ' blocks with a total length of',NPARBLKWT(2,ISTI)
          END DO
        END IF
      END IF
1001  CONTINUE
*      
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
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
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE FNDMND_PAR(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &                     SCR,ISCAT,SUBVAL,NTESTG)
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
      
      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
*
      CALL REWINE(LU,-1)
*
*.LOOP OVER BLOCKS
*
C     write(6,*) ' FNDMND NSUBMX = ', NSUBMX
*
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
*
      IBLK = 0
      IBASE = 1
      LSUB = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Info about block ',IBLK
          WRITE(6,*) ' Number of elements ', LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL ,-1,LU,IMZERO,IAMPACK)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Elements read in '
             CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
           END IF
           IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
             MSUBMX = MIN(NSUBMX,LBL)
C     SUBROUTINE SORLOW(WRK,STVAL,ISTART,KZVAR,KEXSTV,JEXSTV,IPRT)
             IF(LBL.GE.1) THEN
               CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                     MSUBMX,MSUB,NTEST)
C              WRITE(6,*)
C    &       ' After SORLOW 1 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
              MSUB = 0
             END IF
             DO 10 I = 1, MSUB
   10        ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
C            Write(6,*)
C    &       ' After 10 , Scatter array for combined list '
C            call iwrtma(ISCR,1,LSUB+MSUB,1,LSUP+MSUP)
* SORT COMBINED LIST
             MSUBMX = MIN(NSUBMX,LSUB+MSUB)
             IF(MSUBMX.GT.0) THEN
             CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &         NTEST)
C              WRITE(6,*)
C    &       ' After SORLOW 2 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
               LLSUB = 0
             END IF
             LSUB = LLSUB
             DO 20 I = 1, LSUB
               ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20        CONTINUE
*
             CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
             CALL COPVEC(SUBVAL,SCR,LSUB)
             IBASE = IBASE + LBL
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
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
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
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
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
      IMPLICIT REAL*8           (A-H,O-Z)
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
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
C                WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
C    &           LUDIA)
      SUBROUTINE WRSVCD_PAR(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,
     &                     LUFORM,JPACK)
*
* Write scattered vector to disc.
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "../lucita/lucipar.inc"
      
*.Input
      DIMENSION IPLAC(*),VAL(*)
*.Scratch
      DIMENSION VEC1(*)
*
      NTEST = 0
      IF(LBLK.GT.0) THEN
*. Write the vector without markers in one block
         CALL SETVEC(VEC1,0.0D0,NDIM)
         DO 100 IEFF = 1, NSCAT
           VEC1(IPLAC(IEFF)) = VAL(IEFF)
  100    CONTINUE
         CALL TODSC(VEC1,NDIM,-1,LU)
      ELSE
*. Write the vector with the block format of file LUFORM
        CALL REWINO(LUFORM)
        IBOT = 1
*. Loop over records
 1000   CONTINUE
*.Length
        CALL IFRMDS(LBL,1,-1,LUFORM)
        CALL ITODS (LBL,1,-1,LU)
        IF(LBL.GE.0) THEN
C?        write(6,*) ' IBOT, LBL  ',IBOT, LBL
          CALL SETVEC(VEC1,0.0D0,LBL)
          DO 200 IEFF = 1, NSCAT
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
     &      VEC1(IPLAC(IEFF)-IBOT+1) = VAL(IEFF)
*.CSK
CSK            IF(MYPROC.EQ.MASTER) THEN
CSK              IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
CSK     &        write(6,*) ' Catch : IPLAC(IEFF) VAL(IEFF) ',
CSK     &        IPLAC(IEFF),VAL(IEFF)
CSK            END IF
*.CSK
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
*.Skip record on LUFORM
          CALL FRMDSC(VEC1,LBL,-1,LUFORM,IMZERO,IAMPACK)
          IBOT = IBOT + LBL
          GOTO 1000
        END IF
*
        IF(NTEST.GT.10) THEN
          if(myproc.eq.master)then
            WRITE(6,*) ' output from WRSVCD '
            CALL WRTVCD(VEC1,LU,1,LBLK)
          endif
        END IF
      END IF
*
      END
#else
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_DUMMY2
      END 
***********************************************************************
#endif /* defined (VAR_MPI) */
