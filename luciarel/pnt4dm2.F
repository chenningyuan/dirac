      SUBROUTINE PNT4DM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
     &           IDXSM,IS12,IS34,IS1234,IPNTR)
*
* Pointer for 4 dimensional array with total symmetry IDXSM
* Pointer is given as 3 dimensional array corresponding
* to the first 3 indeces
*
* IS12 (0,1,-1)   : Permutational symmetry between indeces 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between indeces 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between indeces 12 and 34
*
* Jeppe Olsen, Feb. 98
*
* Derived from PNT4DM : Symmetry handling transferred to SYMCOM
*
* note ISM4A is not active
*
*.input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
*.Output
      INTEGER IPNTR(NSMOB,NSMOB,NSMOB)
*
      NTEST = 000
      IZERO = 0
      CALL ISETVC(IPNTR,IZERO,NSMOB ** 3 )
*
      IOFF= 1
      DO 10 I1SM = 1, NSMOB
        DO 20 I2SM= 1, NSMOB
*. I2SM is symmetry of creation operator, symmetry of annihilation op
          I2SMA = I_GET_ADJSYM(I2SM)
          CALL SYMCOM(3,1,I1SM,I2SMA,I12SM)
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
C         I12SM = ADASX(I1SM,I2SM)
          CALL SYMCOM(2,1,I12SM,I34SM,IDXSM)
C         I34SM = SXDXSX(I12SM,IDXSM)
          IF(I34SM.EQ.0) GOTO 20
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 20
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMOB+I2SM
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
          DO 30 I3SM = 1, NSMOB
            CALL SYMCOM(2,1,I3SM,I4SMA,I34SM)
            I4SM = I_GET_ADJSYM(I4SMA)
C           I4SM = ADSXA(I3SM,I34SM)
            IF(I4SM.EQ.0) GOTO 30
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 30
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMOB+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 30
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
C             ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF+ N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
C             ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF + N12*(N12+1)/2
            ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
C             ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF=  IOFF+ N12*(N12-1)/2
            END IF
            if (NTEST.ge.100) then
              WRITE(6,*) ' I1SM I2SM I3SM I4SM    IOFF'
              WRITE(6,'(1H ,4I4,I9)')   I1SM,I2SM,I3SM,I4SM,IOFF
            end if
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*
*
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Length of 4 dimensional array ', IOFF - 1
      END IF
*
      RETURN
      END
