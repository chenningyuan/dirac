      SUBROUTINE NXTBLK(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,IBLTP,IDC,
     &                  NONEW,IOCOC)
*
* Obtain allowed block following IATP IBTP IASM
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IBLTP(*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
*
*.Initialize
      ISM = IASM
      IA = IATP
      IB = IBTP
      NONEW = 0
*.Next block
 1000 CONTINUE
      IF(IB.LT.NOCTPB) THEN
        IB = IB + 1
      ELSE
        IB = 1
        IF(IA.LT.NOCTPA) THEN
          IA = IA+ 1
        ELSE
          IA = 1
          IF(ISM.LT.NSMST) THEN
            ISM = ISM + 1
          ELSE
            NONEW = 1
          END IF
        END IF
      END IF
      IF(NONEW.EQ.1) GOTO 1001
*. Should this block be included
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.0) GOTO 1000
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
 1001 CONTINUE
*
      IATP = IA
      IBTP = IB
      IASM = ISM
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)')
     &  ' NXTBLK : ISM IA IB NONEW ', IASM,IA,IB,NONEW
      END IF
*
      RETURN
      END
