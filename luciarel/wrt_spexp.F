      SUBROUTINE WRT_SPEXP(CSPSAO,ICSPSAO,NAOFSM,NBSPFSM,NUSPFSM)
*
* Write Spinor -SAO expansion CSPSAO
*
* Jeppe + Timo, December 1997
*
      IMPLICIT REAL*8  (A-H,O-Z)
*.General input :  Double group - point group translation :
#include "mxpdim.inc"
#include "symm.inc"
*. Specific Input
      DIMENSION CSPSAO(*)
      INTEGER   ICSPSAO(2,2),NAOFSM(*),NBSPFSM(*),NUSPFSM(*)
*
      write(6,*) ' In WRT_SPEXP: '
      WRITE(6,*) ' NAOFSM : ', NAOFSM(1), NAOFSM(2)
      write(6,*) ' NBSPFSM: ',NBSPFSM(1),NBSPFSM(2)
      write(6,*) ' NUSPFSM: ',NUSPFSM(1),NUSPFSM(2)
      IOFF = 1
      DO IBU = 1,2
      DO IRI = 1,2
        IF(IBU.EQ.1) THEN
          WRITE(6,*)
          IF(IRI.EQ.1) THEN
            WRITE(6,*) ' ==================================='
            WRITE(6,*) ' Barred spinor expansion, real part '
            WRITE(6,*) ' ==================================='
          ELSE
            WRITE(6,*) ' ========================================='
            WRITE(6,*) ' Barred spinor expansion, imaginary  part '
            WRITE(6,*) ' ========================================='
          END IF
        ELSE IF(IBU.EQ.2) THEN
          IF(IRI.EQ.1) THEN
            WRITE(6,*) ' ====================================='
            WRITE(6,*) ' Unbarred spinor expansion, real part '
            WRITE(6,*) ' ====================================='
          ELSE
            WRITE(6,*) ' ===================================='
            WRITE(6,*) ' Unbarred spinor expansion, imag part '
            WRITE(6,*) ' ===================================='
          END IF
          WRITE(6,*)
        END IF
        DO ISPSM = NIRR_DG/2+1, NIRR_DG
            WRITE(6,'(1H ,A,I3)') '   Spinors of symmetry ', ISPSM
            WRITE(6,'(1H , A)'  ) '   **************************'
            WRITE(6,*)
*. ao * alpha
          IAOSM_A = ISPSF_MO(ISPSM,1)
*. ao * beta
          IAOSM_B = ISPSF_MO(ISPSM,2)
          IF(IAOSM_A.EQ.0) THEN
            Write(6,*) ' (No alpha component ) '
            LAO_A = 0
          ELSE
            WRITE(6,*) ' Alpha * AO''s with symmetry ', IAOSM_A
            LAO_A = NAOFSM(IAOSM_A)
          END IF
          IF(IAOSM_B.EQ.0) THEN
            Write(6,*) ' (No beta  component ) '
            LAO_B = 0
          ELSE
            WRITE(6,*) ' Beta  * AO''s with symmetry ', IAOSM_B
            LAO_B = NAOFSM(IAOSM_B)
          END IF
C         WRITE(6,*) ' IAOSM_A, IAOSM_B ', IAOSM_A, IAOSM_B
*
          LAO_AB = LAO_A + LAO_B
          IF(IBU.EQ.1) THEN
            LSPINOR = NBSPFSM(ISPSM)
          ELSE
            LSPINOR = NUSPFSM(ISPSM)
          END IF
*
          CALL WRTMAT(CSPSAO(IOFF),LAO_AB,LSPINOR,LAO_AB,LSPINOR)
          IOFF = IOFF + LSPINOR*LAO_AB
        END DO
*       ^ End of loop over spinor symmetries
      END DO
      END DO
*     ^ End of loop over barred/unbarred, real/imaginary parts
      RETURN
      END
