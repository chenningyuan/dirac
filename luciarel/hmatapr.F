      SUBROUTINE HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                   IABSPC,JABSPC,IABOP,JABOP,IIF,JDOH2,
     &                   IDOH2,IMZERO,IDIAG)
*
* Decide upon the treatment of matrix element
*
* <IASM IATP IBSM IBTP | H(apr) | JASM, JATP JBSM, JBTP>
*
* and do preparations (IIF = 1 )
* or counteract preparations (IIF = 2)
*
* Jeppe Olsen, The last day of January 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      IF((IABSPC.NE.JABSPC).OR.
     &   (IABSPC.EQ.JABSPC.AND.(IABOP.EQ.1.OR.IABOP.EQ.2).AND.
     &    (JASM.NE.IASM.OR.IATP.NE.JATP.OR.IBTP.NE.JBTP))) THEN
*. Zero
        IMZERO = 1
      ELSE
*. Not zero !
        IMZERO = 0
      END IF
*. Diagonal approximation?
        IF(IABOP.EQ.1.OR.IABOP.EQ.2) THEN
          IDIAG = 1
        ELSE
          IDIAG = 0
        END IF
*. Moller Plesset or normal operator ?
        IF(IABOP.EQ.1.OR.IABOP.EQ.3.OR.IABOP.EQ.5)THEN
          IMP = 1
        ELSE
          IMP = 0
        END IF
C     END IF
*. Two - or one- electron operator
      IF(IABOP.EQ.1.OR.IABOP.EQ.3) THEN
        IDOH2 = 0
      ELSE
C       IDOH2 = JDOH2
        IDOH2 = 1
      END IF
*. Put MP integrals in place ( Or put good old one-elcectron integrals
*. bach where they belong
      IF(IMP.EQ.1) THEN
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
*
         WRITE(6,*) ' HMATAPR speaking '
         WRITE(6,*) ' ================='
         WRITE(6,*) ' | IASM IATP IBSM IBTP > :',
     &              '|',IASM,IATP,IBSM,IBTP,'>'
         WRITE(6,*) ' | JASM JATP JBSM JBTP > :',
     &              '|',JASM,JATP,JBSM,JBTP,'>'
         WRITE(6,*) ' IABSPC,JABSPC :', IABSPC,JABSPC
         WRITE(6,*) ' IABOP ', IABOP
         WRITE(6,*) ' Results : IMP IDIAG IDOH2 IMZERO : ',
     &                IMP,IDIAG,IDOH2, IMZERO
      END IF
*
      RETURN
      END
