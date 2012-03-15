! DIRAC-side interface routines for the Polarizable Continuum Model

! Driver routine for the Polarizable Continuum Model
   subroutine pcm_driver(dmat, work, lwork, printlevel)
      implicit none
! Passed variables
      real(8)              :: dmat(*)
      real(8)              :: work(*)
      integer              :: lwork
      integer              :: printlevel

! Local variables      
      integer :: i, nr_points
      real(8), allocatable :: tess_cent(:, :)
      real(8)              :: tot_nuc_pol, tot_ele_pol
      real(8), allocatable :: nuc_pot(:), nuc_pol_chg(:)
      real(8), allocatable :: ele_pot(:), ele_pol_chg(:)

! Get cavity size, i.e. the number of tesserae      
      call get_cavity_size(nr_points)

! Allocate tesserae centers matrix, potential and polarization charges vectors
      allocate(tess_cent(3, nr_points))
      allocate(nuc_pot(nr_points))
      allocate(nuc_pol_chg(nr_points))
      allocate(ele_pot(nr_points))
      allocate(ele_pol_chg(nr_points))

! Get tesserae centers      
      call get_tess_centers(tess_cent)

! Compute nuclear and electronic MEP @tesserae centers      
      call nuc_pot_pcm(nr_points, tess_cent, nuc_pot)
      call ele_pot_pcm(nr_points, tess_cent, ele_pot, dmat, work, lwork)

! Compute nuclear and electronic polarization charges @tesserae centers      
      call comp_charge(nuc_pot, nuc_pol_chg)
      call comp_charge(ele_pot, ele_pol_chg)
      
! Print some information
      if (printlevel > 5) then
          do i = 1, nr_points
            print *, "NMEP @point", i, nuc_pot(i)
            print *, "EMEP @point", i, ele_pot(i)
            if (printlevel > 10) then
                print *, "NPC @point", i, nuc_pol_chg(i)
                print *, "EPC @point", i, ele_pol_chg(i)
            end if
          end do
      end if

! Gauss' Theorem Test
      tot_nuc_pol = 0
      tot_ele_pol = 0
      do i = 1, nr_points
        tot_nuc_pol = tot_nuc_pol + nuc_pol_chg(i)
        tot_ele_pol = tot_ele_pol + ele_pol_chg(i)
      end do
      print *, "Total Nuclear Polarization:", tot_nuc_pol
      print *, "Total Electronic Polarization:", tot_ele_pol

   end subroutine


!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
!


!                   CAVITY-RELATED ROUTINES
!
! We shall provide the following data to the module:
!   1. the number of nuclei in the molecule. NUCLEI variable;
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
      real(8), intent(out) :: tess_cent(3, nr_points)
      real(8), intent(out) :: pot(nr_points)
      
      call esp_nuc(nr_points, tess_cent, pot, 0)
      
   end subroutine

   subroutine ele_pot_pcm(nr_points, tess_cent, pot, dmat, work, lwork)
      use vc_grid
      implicit none

      integer, intent(in)  :: nr_points
      real(8), intent(in)  :: tess_cent(3, nr_points)
      real(8), intent(out) :: pot(nr_points)
      real(8), intent(in)  :: dmat(*)
      real(8)              :: work(*)
      integer              :: lwork
      
      call esp_el(nr_points, dmat, tess_cent, pot, work, lwork, 0)

   end subroutine
