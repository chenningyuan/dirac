      SUBROUTINE NGASDT_REL(IOCCMN,IOCCMX,NGAS,ITOTSM,
     &                      NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,
     &                      IAOCC,IBOCC,MXPNGAS,NELFGP,
     &                      NCOMB,XNCOMB,MXSB,MXSOOB,
     &                      IBLTP,NTTSBL,LCOL,IOCOC)
*
* Number of combinations with symmetry ITOTSM and
* occupation between IOCCMN and IOCCMX
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*),NELFGP(*)
*. block types
      DIMENSION IBLTP(*)
*
      CALL QENTER('NGASD')
      NTEST = 000
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
      MXSB = 0
      MXSOOB = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      NTTSBL = 0
      LCOL = 0
      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        CALL SYMCOM_REL(2,1,IASM,IBSM,ITOTSM)
C       write(6,*) ' Iasm ibsm ', iasm,ibsm
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
             MXBTP = IATP
           ELSE
             MXBTP = NOCTPB
           END IF
           LTSSA = NSSOA(IASM,IATP)
           XLTSSA = FLOAT(LTSSA)
           DO 100 IBTP = 1, MXBTP
*
             IF(NTEST.GE.10) THEN
               WRITE(6,*) ' Alpha super group and beta super group'
               CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
               CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
             END IF
*
             IF(IOCOC(IATP,IBTP).EQ.1) THEN
*. Size of unpacked block
               LTTSUP =  LTSSA*NSSOB(IBSM,IBTP)
*. Size of packed block
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 LTTSBL = LTSSA*NSSOB(IBSM,IBTP)
                 XNCOMB = XNCOMB + XLTSSA*FLOAT(NSSOB(IBSM,IBTP))
               ELSE
                 LTTSBL = LTSSA *(LTSSA +1)/2
                 XNCOMB = XNCOMB + XLTSSA*(XLTSSA+1.0D0)*0.5
               END IF
               NCOMB = NCOMB + LTTSBL
               LSB = LSB + LTTSUP
               MXSOOB = MAX(MXSOOB,LTTSUP)
               NTTSBL = NTTSBL + 1
               LCOL = LCOL + NSSOB(IBSM,IBTP)
             END IF
  100      CONTINUE
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
  300 CONTINUE
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB ,NTTSBL',
     &               NCOMB,XNCOMB,NTTSBL
      END IF
*
      CALL QEXIT('NGASD')
*
      RETURN
      END
