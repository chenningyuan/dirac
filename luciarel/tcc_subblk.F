      SUBROUTINE TCC_SUBBLK(TCC,TCC_SUB,IWAY,
     &                      NCA,NCA_SUB,ICA_B,NCB,NCB_SUB,ICB_B,
     &                      NAA,NAA_SUB,IAA_B,NAB,NAB_SUB,IAB_B,
     &                      IAB_TRNSP,IDIAG,IPRNT)
*
* Extract - or add- subbblock of a TCC block
*
* IWAY = 1 Obtain TCC_SUB from TCC
* IWAY = 2 Add TCC_SUB to TCC
*
* Jeppe Olsen, Summer of 99
*              Updated July 2000, IAB_TRNSP, IDIAG added
*
* IAB_TRNSP = 1 : Input block is TCC(ICB,ICA,IAB,IAA)
* IDIAG     = 1 : Input block is lower diagonal
*
#include "implicit.inc"
*. Input and output
      DIMENSION TCC(NCA*NCB*NAA*NAB)
      DIMENSION TCC_SUB(NCA_SUB,NCB_SUB,NAA_SUB,NAB_SUB)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.500) THEN
      WRITE(6,*) ' Input to TCC_SUBBLK '
        WRITE(6,*) ' NCA, NCB, NAA, NAB', NCA,NCB,NAA,NAB
        WRITE(6,*) ' NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB',
     &               NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB
        WRITE(6,*) ' IAB_TRNSP, IDIAG = ', IAB_TRNSP, IDIAG
        WRITE(6,*) ' Iway = ', IWAY
      END IF
*
      IF(IAB_TRNSP.EQ.0.AND.IDIAG.EQ.0) THEN

       DO ICA = 1, NCA_SUB
         DO ICB = 1, NCB_SUB
           DO IAA = 1, NAA_SUB
             DO IAB = 1, NAB_SUB
               IADR_FULL = (IAB+IAB_B-2)*NAA*NCB*NCA + 
     &                     (IAA+IAA_B-2)*NCB*NCA +
     &                     (ICB+ICB_B-2)*NCA +
     &                      ICA+ICA_B-1
C?             WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
               IF(IWAY.EQ.1) THEN
                 TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)                                        
               ELSE
                 TCC(IADR_FULL) = TCC(IADR_FULL) + 
     &                            TCC_SUB(ICA,ICB,IAA,IAB)
               END IF
             END DO
           END DO
         END DO
       END DO
      ELSE IF (IAB_TRNSP.EQ.1.AND.IDIAG.EQ.0) THEN
       DO ICA = 1, NCA_SUB
       DO ICB = 1, NCB_SUB
       DO IAA = 1, NAA_SUB
       DO IAB = 1, NAB_SUB
*. Obtain T(ICA,ICB,IAA,IAB) as T(ICB,ICA,IAB,IAA)
         IADR_FULL = (IAA+IAA_B-2)*NAB*NCA*NCB + (IAB+IAB_B-2)*NCA*NCB
     +             + (ICA+ICA_B-2)*NCB + ICB+ICB_B-1
C?       WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
         IF(IWAY.EQ.1) THEN
           TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)                                        
         ELSE
           TCC(IADR_FULL)  =  TCC(IADR_FULL) + TCC_SUB(ICA,ICB,IAA,IAB)
         END IF
       END DO
       END DO
       END DO
       END DO
      ELSE IF (IDIAG.EQ.1) THEN
       DO ICA = 1, NCA_SUB
       DO ICB = 1, NCB_SUB
       DO IAA = 1, NAA_SUB
       DO IAB = 1, NAB_SUB
         ICA_ABS = ICA + ICA_B -1
         ICB_ABS = ICB + ICB_B -1
         IAA_ABS = IAA + IAA_B -1
         IAB_ABS = IAB + IAB_B -1
*
         IF((IAA_ABS.GT.IAB_ABS).OR.
     &      (IAA_ABS.EQ.IAB_ABS.AND.ICA_ABS.GE.ICB_ABS)) THEN
            IADR_FULL = ITDIANUM(ICA_ABS,ICB_ABS,IAA_ABS,IAB_ABS,
     &                  NCA,NAA)
            ILOW = 1
         ELSE 
            IADR_FULL = ITDIANUM(ICB_ABS,ICA_ABS,IAB_ABS,IAA_ABS,
     &                  NCA,NAA)
            ILOW = 0
         END IF
*
         IF(IWAY.EQ.1) THEN
           TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)                                        
         ELSE IF(ILOW.EQ.1) THEN
           TCC(IADR_FULL)  =  TCC(IADR_FULL) + TCC_SUB(ICA,ICB,IAA,IAB)
         END IF
       END DO
       END DO
       END DO
       END DO
      END IF
*
      IF(NTEST.GE.500) THEN 
        WRITE(6,*) ' TCC_SUBBLK : TCC_SUB '
        NCAB_SUB = NCA_SUB*NCB_SUB
        NAAB_SUB = NAA_SUB*NAB_SUB
        CALL WRTMAT(TCC_SUB,NCAB_SUB,NAAB_SUB,NCAB_SUB,NAAB_SUB)
      END IF
*  
C     CALL LMEMCHK
*
      RETURN
      END
