      SUBROUTINE GENH1(XLAMBDA)
*
* Construct the general one-electron operator
*
* H = XLAMBDA*H(NORMAL) + (1-XLAMBDA)FIFA
*
*
* Where H(Normal) is the normal one-electron operator
* and FIFA is the sum of the inactive and active Fock matrices
* used in CASPTN theory
*
* The correct one-electron density is assumed in place
*
* Jeppe Olsen, March 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "oper.inc"
#include "clunit.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
*. Construct FIFA in WORK(KFI)
*
      LU18 = LUSC95
      IF(IPART.NE.3) THEN
*. Normal M-P Partitioning
        WRITE(6,*) 'Normal MP partitioning'
        WRITE(6,*) 'Normal MP partitioning'
        WRITE(6,*) 'Normal MP partitioning'
        CALL COPVEC(WORK(KINT1),WORK(KFI),NINT1)
        CALL FIFAM(WORK(KFI))
*. And write to disc
        REWIND (LU18)
        WRITE(6,*) ' H0 written to disc '
        CALL TODSC(WORK(KFI),NINT1,-1,LU18)
      ELSE IF(IPART.EQ.3) THEN
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        REWIND (LU18)
        CALL FRMDSC(WORK(KFI),NINT1,-1,LU18,IMZERO,IAMPACK)
      END IF
*
*. And obtain modified operator in KINT1 ( No return !! )
*
      FAC2 = 1.0D0 - XLAMBDA
      CALL VECSUM(WORK(KINT1),WORK(KINT1),WORK(KFI),XLAMBDA,FAC2,
     &            NINT1)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 'Modified matrix as delivered by GENH1 '
       WRITE(6,*) '======================================'
       CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*
      RETURN
      END
