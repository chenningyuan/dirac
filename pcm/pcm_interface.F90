! DIRAC-side interface routines for the Polarizable Continuum Model

!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
!


!                   CAVITY-RELATED ROUTINES
!
! We shall provide the following data to the module:
!   1. the number of nuclei in the molecule. NUCLEI variable ;
!   2. the atomic number of each nucleus. CHARGES vector;
!   3. the coordinates of each nucleus. CENTERS matrix.
! PCMSolver will provide to cavity formation depending
! on the selected mode:
!   1. Mode = Implicit means that the spheres are
!      centered according to the CENTERS matrix with
!      radii taken from the internal library;
!   2. Mode = Atoms means that we specify sphere centers
!      and radii for certain atoms in the molecule. All the
!      other spheres are created as in the Implicit mode.
!

   subroutine collect_nctot(n)

      integer, intent(out) :: n

#include <mxcent.h>
#include <nuclei.h>
      
      n = nucdep

   end subroutine
      
      
   subroutine collect_atoms(charges, centers, unit_flag)

#include <mxcent.h>
#include <nuclei.h>
      
      real(8), intent(out) :: charges(*)
      real(8), intent(out) :: centers(3,*)
      integer, intent(out) :: unit_flag

      integer :: i, j, ii

!.....Get coordinates
      CALL GETACORD(centers)
!.....Get charges      
      ii = 0
      do i=1,nucind
         do j = 1,nucdeg(i)
            ii = ii + 1
            charges(ii) = charge(i)
         enddo
      enddo
   
      ! If units is AU unit_flag = 0. If Angstrom unit_flag = 1.
      unit_flag = 0
   
   end subroutine
