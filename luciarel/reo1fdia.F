*
*  Create new arrays IOBPTS and IOBPTS2 for reordering diagonal
*  one-electron integrals (symmetry/type ordering) to internal
*  LUCIA ordering type/symmetry.
*
      subroutine reo1fdia(IOBPTS,IORD,NGSOB,NSMOB,NGAS,ICOMP)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
*
      dimension IOBPTS(MXPNGAS,*),IORD(MXPOBS,*)
      dimension NGSOB(MXPOBS,MXPNGAS)
*
      ITSOFF = 1
      DO IISYM = 1, NSMOB
        DO IGAS = 1, NGAS
          ISYM = IORD(ICOMP,IISYM)
          IOBPTS(IGAS,ISYM) = ITSOFF
          ITSOFF = ITSOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
      return
      end
