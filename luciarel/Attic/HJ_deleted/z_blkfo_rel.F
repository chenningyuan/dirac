      SUBROUTINE Z_BLKFO_REL(IDC,NMS2VAL,ISPC,ISM,KPCLBT,KPCLEBT,
     &                       KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
     &                       NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block
* KPCIBT  : Info on each block
* KPCBLTP : Block type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
*
*.output
      INTEGER NBLK_MS2(*), IBLK_MS2(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
*
*. Pointers to output arrays
      NTTS = MXNTTS
      call memmar(KPCLBT ,MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,MXNTTS,'ADDL  ',1,'CI1BT ')
      call memmar(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ')
      call memmar(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
      call memmar(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ')

*.    ^ These should be preserved after exit so put mark for 
*       flushing after
*
      call memmar(IDUM,IDUM,'MARK  ',IDUM,'Z_BLKF')
*
      IDUM = 0
*. Allowed length of each batch
      IF(ICISTR.EQ.1) THEN
*. assumes determinants
        LBLOCK =  XISPSM(ISM,ISPC)
      ELSE IF (ICISTR.EQ.2) THEN
        LBLOCK = MXSB
      ELSE IF( ICISTR.EQ.3) THEN
        LBLOCK = MXSOOB
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' LBLOCK = ', LBLOCK
      END IF
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*. Info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
*
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,IDC,WORK(KPCBLTP))
*
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
*. Number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
*. Batches of C vector
        ITTSS_ORD = 2
        CALL PART_CIV3(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                IOFFBLK,IOFFBTC,
     &                NSMST,LBLOCK,WORK(KLCIOIO),
     &                ISMOST(1,ISM),
     &                NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
C            PART_CIV3(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &                  IOFFBLK,IOFFBTC,
C    &                  NSMST,MXLNG,IOCOC,ISMOST,
C    &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
C    &                  ICOMP,ITTSS_ORD)
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

      IF(NTEST.GE.1) THEN
       WRITE(6,*) ' Total number of batches', NBATCH
       WRITE(6,*) ' Total number of blocks ', NBLOCK
*
       WRITE(6,*) ' Number of blocks and batches per MS2 value '
       WRITE(6,*) ' ==========================================='
       WRITE(6,*)
       WRITE(6,*) ' IMS2    MS2      NBLOCK      NBATCH '
       WRITE(6,*) ' ===================================='
       DO IMS2 = 1, NMS2VAL
         IATP = IST_FOR_DT(1,IMS2)
         IBTP = IST_FOR_DT(2,IMS2)
         IIMS2 = NELEC(IATP)-NELEC(IBTP)
         WRITE(6,'(4(3X,I4))') IMS2,IIMS2,NBLK_MS2(IMS2),NBAT_MS2(IMS2)
       END DO
*
      END IF
*
* Transfer information to common block for KRMC - LUCIAREL
* communication:
      call civcinf_mc(WORK(KPCIBT),WORK(KPCLBT),WORK(KPCI1BT),
     &                WORK(KPCLEBT),NBAT_MS2,NMS2VAL,NBLOCK)
*
      call memmar(IDUM,IDUM,'FLUSM ',IDUM,'Z_BLKF')
      RETURN
      END
