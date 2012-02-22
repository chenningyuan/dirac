C DIRAC-side interface routines for the Polarizable Continuum Model

C
C We divide the interface routines into:
C   1. cavity-related routines;
C   2. solver-related routines.
C


C                   CAVITY-RELATED ROUTINES
C
C We shall provide the following data to the module:
C   1. the number of nuclei in the molecule. NUCLEI variable ;
C   2. the atomic number of each nucleus. CHARGES vector;
C   3. the coordinates of each nucleus. CENTERS matrix.
C PCMSolver will provide to cavity formation depending
C on the selected mode:
C   1. Mode = Implicit means that the spheres are
C      centered according to the CENTERS matrix with
C      radii taken from the internal library;
C   2. Mode = Atoms means that we specify sphere centers
C      and radii for certain atoms in the molecule. All the
C      other spheres are created as in the Implicit mode.
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

C.....Get coordinates
      CALL GETACORD(centers)
C.....Get charges      
      ii = 0
      do i=1,nucind
        do j = 1,nucdeg(i)
          ii = ii + 1
          charges(ii) = charge(i)
        enddo
      enddo
      
C If units is AU unit_flag = 0. If Angstrom unit_flag = 1.
      unit_flag = 0
      
      return
      end
