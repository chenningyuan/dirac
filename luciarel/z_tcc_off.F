      SUBROUTINE Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST,IDIAG)
*
* Offsets for symmetryblocks of TCC elements, sym of CA,CB,AA used 
* If Idiag.eq.1 only blocks (Icasm,Iaasm).ge.(Icbsm,Iabsm) are included
*
*
* Jeppe Olsen, Summer of 99
*              July 2000, HNIE : IDIAG added
* Timo Fleig, summer of 2001, modified for relativistic case
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "symm.inc"
*. Input
*. Number of strings of different types, assumed to be 
*. creations strings
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Output
      INTEGER IBT(MXNDGIRR,MXNDGIRR,MXNDGIRR)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Z_TCC_OFF speaking '
        WRITE(6,*) ' NSMST, IDIAG = ', NSMST, IDIAG
      END IF
*
      IOFF = 1
      DO ISM_C = 1, NSMST
*. Symmetry of annihilation part
        ISM_A = IDBGMULT(ITSYM,ISM_C)
        DO ISM_CA = 1, NSMST
          ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
          DO ISM_AA = 1, NSMST
            ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
*
            ISM_ALPHA = (ISM_AA-1)*NSMST+ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST+ISM_CB
*
            IF((IDIAG.EQ.0).OR.
     &         (IDIAG.EQ.1.AND.ISM_ALPHA.GT.ISM_BETA)) THEN
              IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
              IF(IDIAG.EQ.1.AND.ISM_ALPHA.GT.ISM_BETA) 
     &        IBT(ISM_CB,ISM_CA,ISM_AB) = -IOFF
              IOFF = IOFF + 
     &          NCA(ISM_CA)*NCB(ISM_CB)*NAA(ISM_AA)*NAB(ISM_AB)
*                 ^ NAA was number of creation strings per sym
            ELSE IF (IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN
              IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
              LEN = NCA(ISM_CA)*NAA(ISM_AA)
              IOFF = IOFF + LEN*(LEN+1)/2
            END IF
            if (NTEST.ge.100) then
              write(6,'(A,I2,A,I2,A,I2,A,I5)') 
     &                'IBT(',ISM_CA,',',ISM_CB,',',ISM_AA,') =',
     &                 IBT(ISM_CA,ISM_CB,ISM_AA)
            end if
*
          END DO
        END DO
      END DO
*
      RETURN
      END
