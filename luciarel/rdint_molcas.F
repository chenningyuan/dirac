*  Import integrals from MOLCAS/SEWARD environment
*  (Driver routine)
*
      subroutine rdint_molcas(ONEINT,FNONEAO,LU2,NATOMS,IFOCKEX)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      character*6 FNONEAO
*
      dimension ONEINT(ITRIASUM)
*
      call rd1int(ONEINT,FNONEAO,LU2,NATOMS,IFOCKEX)
*
      return
      end
