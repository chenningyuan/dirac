      SUBROUTINE TODSCNP(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
* Only nonzero elements are written to DISC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
*. Normal complete record
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSCP(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
*. zero record of length -LREC(IREC)
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU)
        END IF

*
      END DO
*
      RETURN
      END
