      SUBROUTINE TODSCP(A,NDIM,MBLOCK,IFIL)
*
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
*
* Packed version : Store only nonzero elements
*. Small elements should be xeroed outside
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
* 
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
*
*
C?    write(6,*) ' entering TODSCP '
C?    CALL FLUSH(6)
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. Packing 
        ISCR(2) = 1
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE.0 .OR.ICRAY .EQ. 1 ) THEN
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
*. Loop over packed records of dimension LPBLK
      IELMNT = 0
 1000 CONTINUE
*. The next LPBLK elements 
      LBATCH = 0
*. Obtain next batch of elemnts
  999 CONTINUE
       IF(NDIM.GE.1) THEN
       IELMNT = IELMNT+1
       IF(A(IELMNT).NE.0.0D0) THEN
         LBATCH=LBATCH+1
         IPAK(LBATCH) = IELMNT
         XPAK(LBATCH) = A(IELMNT)
       END IF
       END IF
       IF(LBATCH.EQ.LPBLK.OR.IELMNT.EQ.NDIM) goto 998
       GOTO 999
*. Send to DISC
 998   CONTINUE   
       WRITE(IFIL) LBATCH
       IF(LBATCH.GT.0) THEN 
         WRITE(IFIL) (IPAK(I),I=1, LBATCH)
         WRITE(IFIL) (XPAK(I),I=1, LBATCH)
       END IF
       IF(IELMNT.EQ.NDIM) THEN
         WRITE(IFIL) -1
       ELSE
         WRITE(IFIL) 0
         GOTO 1000
       END IF
*. End of loop over records of truncated elements
      END IF
 1001 CONTINUE
*
C?    CALL FLUSH(6)
      RETURN
      END
