      SUBROUTINE TRPMT3(XIN,NROW,NCOL,XOUT)
*
* XOUT(I,J) = XIN(J,I)
*
*. With a few considerations for large scale cases with cache minimization
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
*
*!!!!!!!!!!!!!!!!!!!!!!!!
      IWAY = 2
*!!!!!!!!!!!!!!!!!!!!!!!!
*
      IF(IWAY.EQ.1) THEN
*. Straightforward, no blocking
        DO IROW =1, NROW
          DO ICOL = 1, NCOL
            XOUT(ICOL,IROW) = XIN(IROW,ICOL)
          END DO
        END DO
      ELSE IF(IWAY.EQ.2) THEN
*. Simple blocking of matrix
        LRBLK = 40
        LCBLK = 40
        NRBLK = NROW/LRBLK
        NCBLK = NCOL/LCBLK
        IF(LRBLK*NRBLK.NE.NROW) NRBLK = NRBLK + 1
        IF(LCBLK*NCBLK.NE.NCOL) NCBLK = NCBLK + 1
*
        DO IRBLK = 1,NRBLK
          IF(IRBLK.EQ.1) THEN
            IROFF = 1
          ELSE
            IROFF = IROFF + LRBLK
          END IF
          IREND = MIN(NROW,IROFF+LRBLK-1)
          DO ICBLK = 1, NCBLK
            IF(ICBLK.EQ.1) THEN
              ICOFF = 1
            ELSE
              ICOFF = ICOFF + LCBLK
            END IF
            ICEND = MIN(NCOL,ICOFF+LCBLK-1)
*
            DO IROW = IROFF,IREND
              DO ICOL = ICOFF,ICEND
                XOUT(ICOL,IROW) = XIN(IROW,ICOL)
              END DO
            END DO
*
          END DO
        END DO
      END IF
*
      RETURN
      END
