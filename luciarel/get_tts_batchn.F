      SUBROUTINE GET_TTS_BATCHN(CTTS,NBLOCK,IBLOCK,NBLOCKC,IBLOCKC,
     &                  NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                  IDC,LUC,IREW,ISCALE)
*
* Read in a batch of C blocks from file LUC.
*
* The complete file is defined by NBLOCKC,IBLOCKC,
* and the blocks of the actual batch is defined by NBLOCK,IBLOCK.
* Vector packing is defined by IDC
*
* Should be initialized with rewind on LUC i.e. IREW = 1
*
*. Feb 96 : Modified to accomodate packed files
*. March 97 : Shaved, assumes IBLOCK and IBLOCKC have identical ordering
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CTTS(*),NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IBLOCKC(8,NBLOCKC),IBLOCK(8,NBLOCK)
*
      COMMON/HIDLUC/IBLK,IATP,IBTP,IASM,IBSM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to GET_TTS_BATCH '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Number of blocks in batch ', NBLOCK
        WRITE(6,*)
        WRITE(6,*) ' IATP  IBTP  IASM  IBSM '
        WRITE(6,*) ' ====================== '
        DO JBLOCK = 1, NBLOCK
          WRITE(6,'(4I4)') (IBLOCK(II,JBLOCK),II = 1, 4 )
        END DO
      END IF
*

*
      IF(IREW.EQ.1) THEN
        REWIND (LUC)
        IBLK = 1
      END IF
*
      IAMPACK = 1
      IMZERO= 0
*. Ready to loop over blocks
*
*. Loop over blocks to be read
      DO JBLOCK = 1, NBLOCK
          IF(IBLOCK(1,JBLOCK).GT.0) THEN
            ISKIP = 0
          ELSE
            ISKIP = 1
          END IF
*
          IF(IBLK.GT.NBLOCKC) THEN
            REWIND LUC
            WRITE(6,*) ' Notice : LUC rewinded in GET_TTS'
            WRITE(6,*) ' Less than optimal programming ? '
            IBLK = 1
          END IF
*. Loop over blocks on file
  999   CONTINUE
*. is this block identical to next block on disc
          IF(ABS(IBLOCK(1,JBLOCK)).EQ.ABS(IBLOCKC(1,IBLK)) .AND.
     &       IBLOCK(2,JBLOCK).EQ.IBLOCKC(2,IBLK) .AND.
     &       IBLOCK(3,JBLOCK).EQ.IBLOCKC(3,IBLK) .AND.
     &       IBLOCK(4,JBLOCK).EQ.IBLOCKC(4,IBLK) ) THEN
             IF(ISKIP.EQ.0) THEN
*. Read me !
              CALL IFRMDS(LBL,1,-1,LUC)
              IOFFO   = IBLOCK(6,JBLOCK)
              CALL FRMDSC(CTTS(IOFFO),LBL,-1,LUC,IMZERO,IAMPACK)
              IBLK = IBLK + 1
            ELSE
*. Skip me, you are not interested in me at all
              CALL IFRMDS(LBL,1,-1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
              IBLK = IBLK + 1
            END IF
          ELSE
*. Skip me
            CALL IFRMDS(LBL,1,-1,LUC)
            CALL SKPRCD2(LBL,-1,LUC)
            IBLK = IBLK + 1
            GOTO 999
          END IF
C     END IF
*. End of loop over blocks on file
      END DO

*
*. Rescale from combination form to determinant
*
      IF(IDC.EQ.2.AND.ISCALE.NE.0) THEN
        CALL SCDTTS(CTTS,IBLOCK,NBLOCK,
     &  NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,
     &  IDC,2,NTEST)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Batch of blocks read in '
         WRITE(6,*) ' ========================'
         WRITE(6,*)
         CALL WRTTTS(CTTS,IBLOCK,NBLOCK,NSMST,NOCTPA,NOCTPB,
     &               NSASO,NSBSO,IDC)
      END IF
*
      RETURN
      END
