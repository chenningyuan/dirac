      SUBROUTINE PART_CIV3(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                  IOFFBLK,IOFFBTC,
     &                  NSMST,MXLNG,IOCOC,ISMOST,
     &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                  ITTSS_ORD)
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
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     PART_CIV      '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(6,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
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
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        if (ITTSS_ORD.eq.1) then
*. Old order : IB, IA, ISM  (leftmost inner loop )
          IF(IB.LT.NOCTPB) THEN
            IB = IB + 1
          ELSE
            IB = 1
            IF(IA.LT.NOCTPA) THEN
              IA = IA+ 1
            ELSE
              IA = 1
              IF(ISM.LT.NSMST) THEN
                ISM = ISM + 1
              ELSE
                IFINI = 1
              END IF
            END IF
          END IF
        ELSE IF(ITTSS_ORD.EQ.2) THEN
*. Else : New order : ISM,IB,IA (leftmost inner loop )
          IF(ISM.LT.NSMST) THEN
            ISM = ISM + 1
          ELSE
            ISM = 1
            IF(IB.LT.NOCTPB) THEN
              IB = IB + 1
            ELSE
              IB = 1
              IF(IA.LT.NOCTPA) THEN
                IA = IA + 1
              ELSE
                IFINI = 1
              END IF
            END IF
          END IF
        END IF
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV2 IDC IBLTP ', IDC,IBLTP(ISM)
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
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV'
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL LMEMCHK
        Call Abend1( 11 )
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
        WRITE(6,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
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
      RETURN
      END
