      SUBROUTINE CSPAOSF_FROM_MOAO(CSPAOSF,NBSPFSM,NUSPFSM,IBSPEXP,
     &           CMOAO,NAOFSM,NMOFSM)
*
* Obtain spinor-ao*spinfunction transformation from
* mo-ao transformation matrix
*
* Jeppe + Timo, November 1997
*
*. The spinor-ao_spinfunction matrices are ordered as
*
* Loop over spin-down/spin-up spinors
*  Loop over real/imaginary part
*   Loop over symmetries of spinors
*    Loop over AO's * alpha functions with correct spinor sym
*    Loop over AO's * beta  functions with correct spinor sym
*   End of loop over symmetries of spinors
*  End of loop over real/imaginary parts
* End of loop over Spin-down/Spin up
*
* Spin down spinors are spinor set1, spin up spinors are spinor set2
* (B as in Barred comes before U as in Unbarred)
*
* There will be 4 spinor-aosf expansions
* 1 : real part, barred spinors
* 2 : Imag part, barred spinors
* 3 : real part, unbarred spinors
* 4 : Imag part, unbarred spinors
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* ======
*. Input
* ======
*
#include "mxpdim.inc"
*. Info on double group symmetry
#include "symm.inc"
#include "input.inc"
      DIMENSION CMOAO(*)
      INTEGER NMOFSM(*),NAOFSM(*)
*. Number of barred and unbarred spinors per symmetry - actually recomputed
      INTEGER NBSPFSM(*), NUSPFSM(*)
*
* =======
*. Output
* =======
*. Offset to the four expansions
      INTEGER IBSPEXP(2,2)
*. And the spinor expansions
      DIMENSION CSPAOSF(*)
*
*. Internal switch : We can initialize with either real or imaginary
*  spinors, IRI_ACT = 1 => Real      initialization
*           IRI_ACT = 2 => Imaginary initialization
*
* barred spinors initialized  with spin down-spin orbitals,
* unbarred spinors initialized with spin-up spin orbitals
*
      NTESTL = 0
      NTEST = max(IPRINT,NTESTL)
*
      if (NTEST.ge.3) then
         WRITE(6,*) ' NBSPFSM '
         CALL IWRTMA(NBSPFSM,1,NIRR_DG,1,NIRR_DG)
         WRITE(6,*) ' NUSPFSM '
         CALL IWRTMA(NUSPFSM,1,NIRR_DG,1,NIRR_DG)
         WRITE(6,*) ' NAOFSM '
         CALL IWRTMA(NAOFSM,1,NIRR_DG/2,1,NIRR_DG/2)
      end if
*

      IRI_ACT = 1
      IOFFSM = 1
      DO IBU=1,2
       DO IRI = 1, 2
         IBSPEXP(IRI,IBU) = IOFFSM
         DO ISPSM = 1, NIRR_DG
*. Sym of AO's that will combine with alpha to give double group sym ISPSM
           IAOSM_A = ISPSF_MO(ISPSM,1)
C?         write(6,*) 'ISPSF_MO has value   ',ISPSF_MO(ISPSM,1)
*. Sym of AO's that will combine with beta  to give double group sym ISPSM
           IAOSM_B = ISPSF_MO(ISPSM,2)
C?         write(6,*) 'ISPSF_MO has value   ',ISPSF_MO(ISPSM,2)
           IF(IAOSM_A.NE.0) THEN
             LAO_A = NAOFSM(IAOSM_A)
             LMO_A = NMOFSM(IAOSM_A)
           ELSE
             LAO_A = 0
             LMO_A = 0
           END IF
           IF(IAOSM_B.NE.0) THEN
             LAO_B = NAOFSM(IAOSM_B)
             LMO_B = NMOFSM(IAOSM_B)
           ELSE
             LAO_B = 0
             LMO_B = 0
           END IF
C          NUSPFSM(ISPSM) = LMO_A
C          NBSPFSM(ISPSM) = LMO_B
*
           LAO_AB = LAO_A + LAO_B
           if (NTEST.ge.3) then
              write(6,*) ' ISPSM = ', ISPSM
              write(6,*) ' IAOSM_A, IAOSM_B', IAOSM_A, IAOSM_B
              write(6,*) ' LMO_A, LMO_B, LAO_A,LAO_B',
     &                     LMO_A, LMO_B, LAO_A,LAO_B
           end if
*
*. Alpha-part
*
C?      write(6,*) ' alpha part '
           IF(LAO_A.NE.0) THEN
             IF(IBU.EQ.1.OR.IRI.NE.IRI_ACT) THEN
*. No alpha contributions to initial barred spinors,
               DO ISP = 1, LMO_A
                 ZERO = 0.0D0
                 CALL SETVEC(CSPAOSF(IOFFSM+(ISP-1)*LAO_AB),
     &                       ZERO,LAO_A)
               END DO
             ELSE IF (IBU .EQ. 2.AND.IRI.EQ.IRI_ACT) THEN
*. Offset for this MO-AO block in the CMOAO expansion matrix
               IBCMOAO = IOFF_BLOCK_MAT(NAOFSM,NMOFSM,IAOSM_A)
               DO ISP = 1, LMO_A
                 CALL COPVEC(CMOAO(IBCMOAO+(ISP-1)*LAO_A),
     &                CSPAOSF(IOFFSM+(ISP-1)*LAO_AB),LAO_A)
C?              write(6,*) ' spinor, copying from, to '
C?              write(6,*) ISP,IBCMOAO+(ISP-1)*LAO,
C?   &                     IOFFSM+(ISP-1)*LAO_AB
               END DO
             END IF
           END IF
*          ^ End of alpha part
*
*. beta  part
*
C?      write(6,*) ' beta part '
           IF(LAO_B.NE.0) THEN
             IF(IBU.EQ.2.OR.IRI.NE.IRI_ACT) THEN
*. No beta contribution to unbarred spinors
               DO ISP = 1, LMO_B
                 CALL SETVEC(CSPAOSF(IOFFSM+(ISP-1)*LAO_AB+LAO_A),
     &                       ZERO,LAO_B)
               END DO
             ELSE IF (IBU .EQ. 1.AND.IRI.EQ.IRI_ACT) THEN
*. Offset for this MO-AO block in the CMOAO expansion matrix
               IBCMOAO = IOFF_BLOCK_MAT(NAOFSM,NMOFSM,IAOSM_B)
               DO ISP = 1, LMO_B
                 CALL COPVEC(CMOAO(IBCMOAO+(ISP-1)*LAO_B),
     &                CSPAOSF(IOFFSM+(ISP-1)*LAO_AB+LAO_A),LAO_B)
C?              write(6,*) ' spinor, copying from, to '
C?              write(6,*) ISP,IBCMOAO+(ISP-1)*LAO_B,
C?   &                     IOFFSM+(ISP-1)*LAO_AB+LAO_A
               END DO
             END IF
           END IF
*          ^ End of beta part
*. Update pointer to start of symmetry block
           IF(IBU.EQ.1) THEN
             IOFFSM = IOFFSM + LMO_B*LAO_AB
           ELSE
             IOFFSM = IOFFSM + LMO_A*LAO_AB
           END IF
         END DO
*.       ^ End of loop over spinor symmetries
       END DO
*      ^ End of real/imaginary switch
      END DO
*     ^ End of barred/unbarred switch
*
      if (NTEST.ge.2) then
        WRITE(6,*)
        WRITE(6,*) ' Initial spinor expansion '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        CALL WRT_SPEXP(CSPAOSF,IBSPEXP,NAOFSM,NBSPFSM,NUSPFSM)
C            WRT_SPEXP(CSPSAO, ICSPSAO,NAOFSM,NBSPFSM,NUSPFSM)
      END IF
*
      RETURN
      END
