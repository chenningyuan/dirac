      function inddir(IORB,ISYM,ITP,IBOSYM)
*
* Input:
* ------
*  IORB:    Number of orbital in spatial symmetry
*  ISYM:    Spatial symmetry
*  ITP:     Occupation type (GAS)
*  IBOSYM:  Ordering array from DIRAC
*
* Output:
* -------
*  INDDIR:  DIRAC orbital index
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "symm.inc"
#include "cgas.inc"
*
      dimension IBOSYM(MORB_DIRAC)
*
      print*,'IBOSYM array:'
      print*,(IBOSYM(I),I=1,MORB_DIRAC/2,1)
*
      NTEST = 1
*
* Determine absolute orbital number in terms of GAS/SYM
      IGSCT = 0
      do IGS=1,ITP-1,1
        IGSCT = IGSCT + NGSOB(ISYM,IGS)
      end do
      IORB = IORB + IGSCT
*
      IORBCT = 0
      do ICT=1,MORB_DIRAC/2,1
        if (IBOSYM(ICT).eq.(ISYM-1)) then
          IORBCT = IORBCT + 1
          if (IORBCT.eq.IORB) then
            INDDIR = ICT
            goto 10
          end if
        end if
      end do
*
  10  continue
*
      if (NTEST.ne.0) then
        write(6,*) 'Determined final index:    ',INDDIR
      end if
*
      return
      end
