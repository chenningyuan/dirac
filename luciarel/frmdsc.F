      SUBROUTINE FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
*. Version allowing zero and packed blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY(*)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
C       IF(MMBLOCK.GE.2) MMBLOCK = 2
C       CALL IFRMDS(ISCR,2,MMBLOCK,IFILE)
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          ZERO = 0.0D0
          CALL SETVEC(ARRAY,ZERO,NDIM)
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1
      IF(I_AM_PACKED.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements
  999   CONTINUE
          NBATCH = NBATCH + 1
          IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
          END IF
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              Call Abend2( ' problem in FRMDSC ' )
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
        END IF
C
        IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
         CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
        END IF
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
