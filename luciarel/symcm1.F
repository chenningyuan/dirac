      SUBROUTINE SYMCM1(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1*I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* D2h version , written for compatibility with general symmetry
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
*
      IF(ITASK.EQ.1) THEN
        I1 = SYMPRO(I2,I12)
      ELSE IF(ITASK.EQ.2) THEN
        I2 = SYMPRO(I1,I12)
      ELSE IF (ITASK.EQ.3) THEN
        I12 = SYMPRO(I1,I2)
      END IF
*
      RETURN
      END
