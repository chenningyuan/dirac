C PCMSolver Interface routines

C
C Different Cavity Modes. Mode = Implicit
C nSpheres is now NCTOT defined in nuclei.h
C
      subroutine collect_nctot(nuclei)

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      
      nuclei = nucdep

      return
      end
      
      
      subroutine collect_atoms(charges, centers, unit_flag)

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <inforb.h>
#include <maxorb.h>
#include <infpar.h>
      
      dimension charges(*)
      dimension centers(3,*)
      integer unit_flag

      do i=1,nucdep
         charges(i) = charge(i)
         do j=1,3
            centers(j,i) = cord(j,i)
         enddo
      enddo
      
C If units is AU unit_flag = 0. If Angstrom unit_flag = 1.
      unit_flag = 0
      
      return
      end
