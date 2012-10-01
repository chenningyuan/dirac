      SUBROUTINE DIM_T1T2_TO_T12_S_REL(LEN_T1T2_STRING,LEN_T1T2_TCC,
     &           NSPOBEX_TP,ISPOBEX_TP,NGAS)
*
* LEN_T1T2_STRING : Max dimension of MAT(T1,T2) where T1 and 
*                   T2 are spin strings, and T1*T2 is an allowed 
*                   spin strings. T1 and T2 are assumed to 
*                   run over all symmetries 
*                   
*. Jeppe Olsen, May 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "csm.inc"
*. Input 
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. Local scratch
      INTEGER I12OCC(MXPNGAS)
      INTEGER NI(8)
      INTEGER NJ(8)
* 
      LEN_T1T2_TCC = 0
*     ^ Not implemented yet
*
      LEN_T1T2_STRING = 0
      IONE = 1
      DO ICAAB = 1, 4
       IF(ICAAB.GE.3) THEN
        ICA = 2
       ELSE
        ICA = 1
       END IF
       call evenodd(IEO,I_CAAB)
       IF(IEO.EQ.1) THEN
         IUB = 1
       ELSE
         IUB = 2
       END IF
       DO ITP = 1, NSPOBEX_TP
        DO JTP = 1, NSPOBEX_TP
C?        WRITE(6,*) ' ITP, JTP, ICAAB = ', ITP, JTP, ICAAB
*. ITP*JTP, component ICAAB
          IB = 1 + (ICAAB-1)*NGAS
          CALL IVCSUM(I12OCC,ISPOBEX_TP(IB,ITP),ISPOBEX_TP(IB,JTP),
     &                IONE,IONE,NGAS)
*. Is I12OCC included in any of the strings
C              STRING_IN_CCEXCIT(ISTRING,ICAAB,INSPC,
C    &                           NSPOBEX_TP,ISPOBEX_TP)
*
          CALL STRING_IN_CCEXCIT_REL(I12OCC,ICAAB,INSPC,NSPOBEX_TP,
     &                           ISPOBEX_TP) 
*
          IF(INSPC.EQ.1) THEN 
*. String is in space, find Dimensions of ITP and JTP per sym
            CALL NST_SPGP_CC_REL(ISPOBEX_TP(IB,ITP),NI,IUB,ICA)
            CALL NST_SPGP_CC_REL(ISPOBEX_TP(IB,JTP),NJ,IUB,ICA)
            NI_TOT = IELSUM(NI,NSMST)
            NJ_TOT = IELSUM(NJ,NSMST)
            print*,'NI_TOT,NJ_TOT',NI_TOT,NJ_TOT
            LEN = NI_TOT*NJ_TOT
            LEN_T1T2_STRING = MAX(LEN_T1T2_STRING,LEN)
          END IF
        END DO
       END DO
*      ^ End of loop over ITP, JTP
      END DO
*     ^ End of loop over ICAAB
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' LEN_T1T2_STRING = ', LEN_T1T2_STRING
      END IF
C      DO I_CAAB = 1, 4
C        DO ISM = 1, 4
C          print*,'ISPOBEX_TP',ISPOBEX_TP(ISM,I_CAAB)
C      end do
C      end do
*
      RETURN
      END
