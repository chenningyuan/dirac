      SUBROUTINE REOH1DIA(VECIN,VECOUT)
*
* Diagonal of one-electron integrals given
* in the order
*
* Loop over type
*  Loop over symmetry
*    Loop over unbarred of this TS
*    Loop over   barred of this TS
*  End of loop over symmetry
* End of loop over type
*
* Jeppe + Timo, A hot day in Bonn, May12 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT(*)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
*
*. Loop over orbitals in output order
      JOFF = 1
      DO JGAS = 1, NGAS
*. First orbital of type JGAS in input and output arrays
        IOFF = 1
        DO IGAS = 1, JGAS - 1
          DO ISMOB = 1, NSMOB
            IOFF = IOFF +
     &      NOBPTS(IGAS,ISMOB)+NOBPTS2(IGAS,ISMOB)
          END DO
        END DO
*
        DO JSMOB = 1, NSMOB
          NJORB_U = NOBPTS (JGAS,JSMOB)
          NJORB_B = NOBPTS2(JGAS,JSMOB)
          IJORB_U = IOBPTS(JGAS,JSMOB)
          IJORB_B = IOBPTS2(JGAS,JSMOB)
*. Address of first unbarred orbital of this TS in input order
          IIORB_U = 1
          DO ISMOB = 1, JSMOB-1
           IIORB_U = IIORB_U +
     &     NOBPTS(JGAS,ISMOB)+NOBPTS2(JGAS,ISMOB)
          END DO


*
*. unbarred orbitals
          I



        JOFF = JOFF + NJORB

