      SUBROUTINE WRTRS2(VECTOR,ISMOST,ICBLTP,IOCOC,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,NSMST)
*
* Write RAS vector . Storage form is defined by ICBLTP
*
      IMPLICIT REAL*8           (A-H,O-Z)
*
      LOGICAL DIAGBL
*
      DIMENSION VECTOR(*)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION NSASO(NSMST,* ),NSBSO(NSMST,* )
      DIMENSION ICBLTP(*),ISMOST(*)
*
*
      IBASE = 1
      DO 1000 IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBSM.EQ.0.OR.ICBLTP(IASM).EQ.0) GOTO 1000
*
        DO 900 IATP = 1, NOCTPA
          IF(ICBLTP(IASM).EQ.2) THEN
            IBTPMX = IATP
          ELSE
            IBTPMX = NOCTPB
          END IF
          NAST = NSASO(IASM,IATP)
*
          DO 800 IBTP = 1 , IBTPMX
            IF(IOCOC(IATP,IBTP) .EQ. 0 ) GOTO 800
            NBST = NSBSO(IBSM,IBTP)
            IF(ICBLTP(IASM).EQ.2.AND.IATP.EQ.IBTP ) THEN
* Diagonal block
              NELMNT = NAST*(NAST+1)/2
              IF(NELMNT.NE.0) THEN
                WRITE(6,'(A,3I3)')
     &          '  Iasm iatp ibtp : ', IASM,IATP,IBTP
                WRITE(6,'(A)')
     &          '  ============================'
                CALL PRSM2(VECTOR(IBASE),NAST)
                IBASE = IBASE + NELMNT
              END IF
            ELSE
              NELMNT = NAST*NBST
              IF(NELMNT.NE.0) THEN
                WRITE(6,'(A,3I3)')
     &          '  Iasm iatp ibtp : ', IASM,IATP,IBTP
                WRITE(6,'(A)')
     &          '  ============================'
                CALL WRTMAT(VECTOR(IBASE),NAST,NBST,NAST,NBST)
                IBASE = IBASE + NELMNT
              END IF
            END IF
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
*
      RETURN
      END
