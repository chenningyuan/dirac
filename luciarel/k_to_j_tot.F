      SUBROUTINE K_TO_J_TOT(IKJ,XKJ,KSM,IK_B,IK_E,IT,ITSM,NTOP,
     &                      IT_B,IT_E,IM,XM,IBM,NK,LTOP,IZERO_MAP,IUBL)
*
* Obtain Total map !J> = T-Oper !K>
* T-oper : string of elementary operators 
* !K> String of sym KSM
*
* !K> are of sym KSM and are restricted to IK_B to IK_E ( within sym)
* T_operators are of sym ITSM and are restricted to IT_B to IT_E
*
* Map for each elementary operator is provided by 
*  IM with offset mat IBM
*
* Jeppe Olsen, Summer of 99
* Double group symmetry incorporated for relativistic calculations
* IUBL flag for unbarred/barred operator strings introduced.
*  Timo Fleig, July 2001, January 2002
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "cgas.inc"
#include "csm.inc"
*. Input
      INTEGER IT(LTOP,*),IM(*),
     &        IBM(MXNDGIRR,MXNDGIRR,LTOP),NK(MXPOBS,LTOP)
      DIMENSION XM(*)
*. Output
      INTEGER IKJ(NTOP,(IK_E-IK_B+1))
      DIMENSION XKJ(NTOP,(IK_E-IK_B+1))
*
      NTEST = 0000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Information from K_TO_J.... '
        WRITE(6,*) ' =========================== '
        WRITE(6,*) ' LTOP, NTOP = ', LTOP,NTOP
        WRITE(6,*) ' IK_B,IK_E,IT_B,IT_E = ', IK_B,IK_E,IT_B,IT_E
        WRITE(6,*) ' IUBL, KSM = ',IUBL,KSM
        write(6,*) 'Offset matrix IBM:'
        do L=1,LTOP,1
          write(6,*) 'Ltop value ',L
          call iwrtma(IBM(1,1,L),NIRR_DG,NIRR_DG,MXNDGIRR,MXNDGIRR)
        end do
      END IF
C?    WRITE(6,*) ' Dim of IKJ, XKJ ',NTOP,IK_E-IK_B+1
*
      IZERO_MAP = 1
      IF(LTOP.GT.0) THEN
        DO KSTR = IK_B, IK_E
C?        WRITE(6,*) ' KSTR = ', KSTR
          DO ITOP = IT_B, IT_E
C?          WRITE(6,*) ' ITOP = ', ITOP
            KNSM = KSM
            SIGN = 1.0D0
            KNSTR = KSTR
C?          WRITE(6,'(A,10I4)') ' T-op : ', (IT(JOP,ITOP),JOP=1,LTOP)
            DO IOP = 1, LTOP
C?            WRITE(6,*) ' Info for IOP = ', IOP
              JOB_ABS = IT(IOP,ITOP)
              if (IUBL.eq.1) then
                JOB_SM = ISMFTO(JOB_ABS)
                JOB_TP = ITPFTO(JOB_ABS)
              else if (IUBL.eq.2) then
                JOB_SM = ISMFTO2(JOB_ABS)
                JOB_TP = ITPFTO2(JOB_ABS)
              end if
C?            WRITE(6,*) ' JOB_ABS, IBTS, JOB_SM ',
C?   &        JOB_ABS, JOB_TP, JOB_SM
              if (IUBL.eq.1) then
                JOB_REL = JOB_ABS - IOBPTS(JOB_TP,JOB_SM) + 1
              else if (IUBL.eq.2) then
                JOB_REL = JOB_ABS - IOBPTS2(JOB_TP,JOB_SM) + 1
              end if
              LK = NK(KNSM,IOP)
C?            WRITE(6,*) ' IBM,JOB_REL,LK,KNSTR',
C?   &        IBM(JOB_SM,KNSM,IOP),JOB_REL,LK,KNSTR
              IADR = IBM(JOB_SM,KNSM,IOP) -1 + (JOB_REL-1)*LK+KNSTR
C?            WRITE(6,*) ' IADR = ', IADR
              KNSTR = IM(IADR)
C?            WRITE(6,*) ' IADR, KNSTR =  ', IADR, KNSTR
              SIGN = SIGN*XM(IADR) 
              IF(KNSTR.EQ.0) GOTO 1001
              KNSM = IDBGMULT(KNSM,JOB_SM)
            END DO
*
 1001       CONTINUE
            IF(KNSTR.NE.0)  IZERO_MAP = 0
            IKJ(ITOP-IT_B+1,KSTR-IK_B+1) = KNSTR
            XKJ(ITOP-IT_B+1,KSTR-IK_B+1) = SIGN
C?          WRITE(6,*) ' K_TO_J Row and col : ', 
C?   &      ITOP-IT_B+1,KSTR-IK_B+1
          END DO
        END DO
      ELSE
*. No K-operators => Identity map
        IZERO_MAP = 0
        DO KSTR = IK_B,  IK_E
          IKJ(1,KSTR-IK_B+1) = KSTR
          XKJ(1,KSTR-IK_B+1) = 1.0D0
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Operators '
        print*,'IT_B is ',IT_B
        CALL IWRTMA(IT(1,IT_B),LTOP,IT_E-IT_B+1,LTOP,IT_E-IT_B+1)
        WRITE(6,*) ' Output from K_TO_J_TOT ' 
        WRITE(6,*) 'IKJ and XKJ arrays '
        LT = IT_E-IT_B+1
        LK = IK_E-IK_B+1
        CALL IWRTMA(IKJ,LT,LK,LT,LK)
        CALL WRTMAT(XKJ,LT,LK,LT,LK)
      END IF
*
      RETURN
      END
* 
