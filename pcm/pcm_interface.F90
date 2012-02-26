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

   subroutine collect_nctot(nuclei)
      implicit none
      integer, intent(out) :: nuclei

#include <mxcent.h>
#include <nuclei.h>
      
      nuclei = nucdep

   end subroutine
      
      
   subroutine collect_atoms(charges, centers, unit_flag)
      implicit none

#include <mxcent.h>
#include <nuclei.h>
      
      real(8), intent(out) :: charges(*)
      real(8), intent(out) :: centers(3,*)
      integer, intent(out) :: unit_flag

      integer :: i, j, k 

!.....Get coordinates
      call getacord(centers)
!.....Get charges      
      i = 0
      do j=1,nucind
         do k = 1,nucdeg(j)
            i = i + 1
            charges(i) = charge(j)
         enddo
      enddo
   
      ! If units is AU unit_flag = 0. If Angstrom unit_flag = 1.
      unit_flag = 0
   
   end subroutine

!      
!                   SOLVER-RELATED ROUTINES
!
! We shall provide the following data to DIRAC:
!   1. grid points;
!   2. apparent charges @grid points.      
! We shall provide the following data to the module:
!   1. the nuclear potential @grid points ;
!   2. the electronic potential @grid points.      
!

   subroutine nuc_pot_pcm(nr_points, tess_cent, pot)
      use vc_grid
      implicit none

      integer, intent(in)  :: nr_points
      real(8), intent(out) :: tess_cent(3, *)
      real(8), intent(out) :: pot(*)
      
      call esp_nuc(nr_points, tess_cent, pot, 0)
      
   end subroutine

   subroutine el_pot_pcm(nr_points, tess_cent, pot, dmat, work, lwork)
      use vc_grid
      implicit none

      integer, intent(in)  :: nr_points
      real(8), intent(in)  :: tess_cent(3, *)
      real(8), intent(out) :: pot(*)
      real(8), intent(in)  :: dmat(*)
      real(8)              :: work(*)
      integer              :: lwork
      
      call esp_el(nr_points, dmat, tess_cent, pot, work, lwork, 0)

   end subroutine
