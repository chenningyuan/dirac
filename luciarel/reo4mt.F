      SUBROUTINE REO4MT(XIN,XOUT,LIN,IREO,FACIN,FACOUT)
*
* Reorder a 4 full index matrix and add to output matrix
*
* XOUT(IREO,JREO,KREO,LREO) = FACOUT*XOUT(IREO,JREO,KREO,LREO)
*                           + FACIN *XIN(I,J,K,L)
*
* Ordering of indeces is XIN(1,2,3,4)
*
* IREO : Input => output order
*
* fixed bug. IREO wasn't dimensioned ...    Timo Fleig, Aug. '98
*
      IMPLICIT REAL*8(A-H,O-Z)
*.input
      DIMENSION XIN(*)
      dimension LIN(*),IREO(*)
*. output
      DIMENSION XOUT(*)
*.local scratch
      DIMENSION ISCR(4),LOUT(4)
*. Dimension of output matrix
      DO IND = 1, 4
        LOUT(IREO(IND)) = LIN(IND)
      END DO
*. Loop over element in input order
      IELMNT_IN = 0
*
* I4
*  I3
*   I2
*    I1
*
      DO I4 = 1, LIN(4)
       DO I3 = 1, LIN(3)
        DO I2 = 1, LIN(2)
         DO I1 = 1, LIN(1)
           IELMNT_IN = IELMNT_IN + 1
*. output order
           ISCR(IREO(1)) = I1
           ISCR(IREO(2)) = I2
           ISCR(IREO(3)) = I3
           ISCR(IREO(4)) = I4
**.output address
           IELMNT_OUT = (ISCR(4)-1)*LOUT(3)*LOUT(2)*LOUT(1)
     &                + (ISCR(3)-1)        *LOUT(2)*LOUT(1)
     &                + (ISCR(2)-1)                *LOUT(1)
     &                +  ISCR(1)
           IF(FACOUT.EQ.0.0) THEN
            XOUT(IELMNT_OUT) = FACIN * XIN(IELMNT_IN)
           ELSE
            XOUT(IELMNT_OUT) = FACIN * XIN(IELMNT_IN)
     &                       + FACOUT*XOUT(IELMNT_OUT)
           END IF
         END DO
        END DO
       END DO
      END DO
*
      RETURN
      END
