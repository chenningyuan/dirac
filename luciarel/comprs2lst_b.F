      SUBROUTINE COMPRS2LST_B(I1I,XI1I,N1,I2I,XI2I,N2,NKIN,NKOUT,
     &                        I1O,XI1O,I2O,XI2O)
*
* Two lists of excitations/annihilations/creations are given.
* Compress to common nonvanishing entries
*
* Jeppe Olsen, July 2000 from COMPRS2LST 
*
* Compared to COMPRS2LST : Order of arrys in I1 and I2 interchanged
*                          Output lists differs from input lists
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input lists
      DIMENSION I1I(N1,NKIN),XI1I(N1,NKIN)
      DIMENSION I2I(N2,NKIN),XI2I(N2,NKIN)
*. Output lists
      DIMENSION I1O(N1,NKIN),XI1O(N1,NKIN)
      DIMENSION I2O(N2,NKIN),XI2O(N2,NKIN)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1I(I,K).NE.0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2I(I,K).NE.0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
            DO I = 1, N1
               I1O(I,NKOUT) = I1I(I,K)
              XI1O(I,NKOUT) =XI1I(I,K)
            END DO
            DO I = 1, N2
               I2O(I,NKOUT) = I2I(I,K)
              XI2O(I,NKOUT) =XI2I(I,K)
            END DO
        END IF
      END DO
*
      RETURN
      END
