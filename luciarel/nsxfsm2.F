      FUNCTION NSXFSM2(NSMOB,NO1PS,NO2PS,ISXSM,ISYM,IPRNT)
*
* Number of single excitations of symmetry ISXSM
*
* ISYM = 0 : All symmetry allowed excitations
* ISYM = 1 : Only excitations a+iaj with I.ge.J
* ISYM =-1 : Only excitations a+iaj with I.gt.J
*
* Compared to NSXFSM : Reference to symmetries (ADSXA removed,
*                      replaced with call to SYMCOM
*
COLD  INTEGER ADSXA(MXPOBS,2*MXPOBS)
      INTEGER NO1PS(*),NO2PS(*)
*
      MSXFSM = 0
C     WRITE(6,*) ' NO1PS in NSXFSM2 '
C     CALL IWRTMA(NO1PS,1,NSMOB,1,NSMOB)
C     WRITE(6,*) ' NSMOB ',NSMOB
      DO 100 IO1SM = 1,NSMOB
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM(2,1,IO1SM,IO2SMA,ISXSM)
COLD    IO2SMA = (IO1SM,ISXSM)
*. IO2SMA is symmetry of annihilation operator, find symmetry of creation op
        IO2SM = I_GET_ADJSYM(IO2SMA)
C               I_GET_ADJSYM(ISYM)
C?      WRITE(6,*) ' IO1SM,IO2SM',IO1SM,IO2SM
        IF(ISYM.EQ.0.OR.IO1SM.GT.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*NO2PS(IO2SM)
        ELSE IF( ISYM.EQ. 1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)+1)/2
        ELSE IF( ISYM.EQ.-1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)-1)/2
        END IF
  100 CONTINUE
*
      NSXFSM2 = MSXFSM
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)

      IF(NTEST.NE.0) THEN
        WRITE(6,'(A32,I8,A1,I8)')
     &  ' Number of single excitations of symmetry ',
     &    ISXSM,',',NSXFSM2
      END IF
*
      RETURN
      END
