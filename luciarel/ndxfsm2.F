      FUNCTION NDXFSM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
     &         IS12,IS34,IS1234,IPRNT)
*
* Number of double excitations with total symmetry IDXSM
*
* IS12 (0,1,-1)   : Permutational symmetry between index 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between index 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between index 12 and 34
*
* Compared to NDXFSM : Reference to symmetry arrays eliminated by
*                      calls to SYMCOM ( easier when several groups are
*                      used )
*. General input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
*
      MDX = 0
      IDXSM = 1
      DO 200 I12SM = 1, NSMOB
        DO 190 I1SM = 1, NSMOB
C              SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL SYMCOM(2,1,I1SM,I2SMA,I12SM)
          I2SM = I_GET_ADJSYM(I2SMA)
C         I2SM = ADSXA(I1SM,I12SM)
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 190
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMSX+I2SM
          ELSE
           I12NUM =  I1SM*(I1SM+1)/2+I2SM
          END IF
          IF(IS12.EQ.0.OR.I1SM.NE.I2SM) THEN
            N12 = NO1PS(I1SM)*NO2PS(I2SM)
          ELSE IF(IS12.EQ.1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
          ELSE IF(IS12.EQ.-1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
          END IF
C              I34SM = SXDXSX(I12SM,IDXSM)
C              SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL SYMCOM(2,1,I12SM,I34SM,IDXSM)
          DO 90 I3SM = 1, NSMOB
C           I4SM = ADSXA(I3SM,I34SM)
            CALL SYMCOM(2,1,I3SM,I4SMA,I34SM)
            I4SM = I_GET_ADJSYM(I4SMA)
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 90
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMSX+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 90
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              MDX = MDX + N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              MDX =  MDX + N12*(N12+1)/2
              ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              MDX =  MDX + N12*(N12-1)/2
            END IF
C?          WRITE(6,*) ' I1SM I2SM I3SM I4SM MDX '
C?          WRITE(6,*)   I1SM,I2SM,I3SM,I4SM,MDX
   90       CONTINUE
C 100     CONTINUE
  190   CONTINUE
  200 CONTINUE
*
      NDXFSM2 = MDX
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.NE.0) THEN
         WRITE(6,'(A39,I10)')
     &   ' Number of double excitations obtained ', MDX
      END IF
*
      RETURN
      END
