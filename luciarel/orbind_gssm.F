* Find range of active orbital indices for a given symmetry
* INSYM from GAS/SYM ordering.
*   Timo Fleig, August 2001
* Indices of unoccupied spaces are set to zero, the corresponding
* integrals as well.
*   Timo Fleig, October 2001
*
      subroutine orbind_gssm(INSYM,NSM,INDGATH,NELM,IPRNT)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      dimension INDGATH(NELM)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      IZERO = 0
      IGSOFF = 0
      do IGAS = 1,NGAS,1
        NORBSM = 0
        do JGAS = 1,IGAS-1,1
          do ISM = 1,NSM,1
            NORBSM = NORBSM + NGSSH(ISM,JGAS)
          end do
        end do
        do ISM = 1,INSYM-1,1
          NORBSM = NORBSM + NGSSH(ISM,IGAS)
        end do
        ISTA = NORBSM + 1
        IEND = NORBSM + NGSSH(INSYM,IGAS)
* write to array
        do I = 1,NGSSH(INSYM,IGAS),1
          if (MXNELGS(IGAS).ne.0) then
            INDGATH(I+IGSOFF) = ISTA + I - 1
          else
            INDGATH(I+IGSOFF) = IZERO
          end if
        end do
        IGSOFF = IGSOFF + NGSSH(INSYM,IGAS)
      end do
*
      if (NTEST.ge.5) then
        write(6,*) 'Output from orbind_gssm:'
        write(6,*) 'Symmetry ',INSYM
        write(6,*) 'INDGATH array :'
        call iwrtma(INDGATH,1,NELM,1,NELM)
      end if
*
      return
      end
