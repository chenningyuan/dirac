! DIRAC-side interface routines for the Polarizable Continuum Model
! RDR 2012
!
   
   subroutine pcm_driver(dmat, fock, work, lwork, printlevel)
      implicit none
!
! Driver routine for the Polarizable Continuum Model.
!
! RDR 190312 If X2C or using symmetry quits. 
!
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef prg_dirac
#include "dcbgrd.h"
#endif
#include "ccom.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
#include "shells.h"
#include "onecom.h"
#include "primit.h"
#include "symind.h"
#include "symmet.h"
#include "lmns.h"
! Passed variables
      real(8)              :: dmat(*), fock(*)
      real(8)              :: work(*)
      integer              :: lwork
      integer              :: printlevel

! Local variables      
      integer :: i, nr_points
      real(8), allocatable :: tess_cent(:, :)
      real(8)              :: tot_nuc_pol, tot_ele_pol
      real(8), allocatable :: nuc_pot(:), nuc_pol_chg(:)
      real(8), allocatable :: ele_pot(:), ele_pol_chg(:)
      real(8), allocatable :: charges(:)
      real(8), allocatable :: fock_pcm(:, :, :)
      real(8)              :: ene_pol

! Quit if symmetry or X2C are requested.
      if (nbsym > 1) then
         call quit('Polarizable Continuum Model calculation cannot      &
     & handle symmetry!')
      end if
      if (bss.or.x2c) then
         call quit('Polarizable Continuum Model calculation not working &
     & with X2C Hamiltonian!')
      end if

! Get cavity size, i.e. the number of tesserae      
      call get_cavity_size(nr_points)

! Allocate tesserae centers matrix, potential and polarization charges vectors
      allocate(tess_cent(3, nr_points))
      allocate(nuc_pot(nr_points))
      allocate(nuc_pol_chg(nr_points))
      allocate(ele_pot(nr_points))
      allocate(ele_pol_chg(nr_points))
! Allocate total charges and PCM contribution to Fock matrix      
      allocate(charges(nr_points))
      allocate(fock_pcm(ntbas(0), ntbas(0), nz))
      fock_pcm = 0.0d0

! Get tesserae centers      
      call get_tess_centers(tess_cent)

! Compute nuclear and electronic MEP @tesserae centers      
      nuc_pot = 0.0d0
      call nuc_pot_pcm(nr_points, tess_cent, nuc_pot)
      ele_pot = 0.0d0
! Here we are calculating: tr(DV_i) [eq. (1.2)]
      call ele_pot_pcm(nr_points, tess_cent, ele_pot, dmat, work, lwork)
! Compute nuclear and electronic polarization charges @tesserae centers
      nuc_pol_chg = 0.0d0
! Store potentials and charges in the Cavity object?
! 0 -> no; 1 -> nuclear; 2 -> electronic; 3-> nuclear + electronic
      call comp_charge(nuc_pot, nuc_pol_chg, 1)
      ele_pol_chg = 0.0d0
      call comp_charge(ele_pot, ele_pol_chg, 2)

! Fock matrix contribution     
      charges = nuc_pol_chg + ele_pol_chg
      call pcm_fock(nr_points, tess_cent, charges, fock_pcm, work, lwork)
      call prqmat(fock_pcm, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                       nz, ipqtoq(1,0), 6)
! Add fock_pcm to Fock matrix
      call daxpy(n2bbasx, 1, fock_pcm, 1, fock, 1)

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
      call gauss_test(nr_points, nuc_pol_chg, ele_pol_chg)

! Polarization energy: U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
      call comp_pol_ene_pcm(ene_pol)
      print *, "U_pol =", ene_pol
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

#include "mxcent.h"
#include "nuclei.h"
      
      nuclei = nucdep

   end subroutine
      
      
   subroutine collect_atoms(charges, centers, unit_flag)
      implicit none

#include "mxcent.h"
#include "nuclei.h"
      
      real(8), intent(out) :: charges(*)
      real(8), intent(out) :: centers(3,*)
      integer, intent(out) :: unit_flag

      integer :: i, j, k 

! Get coordinates
      call getacord(centers)
! Get charges      
      i = 0
      do j = 1, nucind
         do k = 1, nucdeg(j)
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
      
      call esp_el(nr_points, dmat, tess_cent, pot, work, lwork, 0, .false.)

   end subroutine
  
   subroutine pcm_fock(nr_points, tess_cent, charges, fock_pcm, work, lwork)
      use vc_grid
      implicit none

#include "dcbbas.h"
#include "dgroup.h"

! Passed variables
      integer, intent(in)  :: nr_points
      real(8), intent(in)  :: tess_cent(3, nr_points)
      real(8), intent(in)  :: charges(nr_points)
      real(8), intent(out) :: fock_pcm(ntbas(0), ntbas(0), nz)
      real(8)              :: work(*)
      integer              :: lwork

! Local variables
      real(8) :: average
      integer :: i, j

      call esp_el(nr_points, fock_pcm, tess_cent, charges, work, lwork, 0, .true.)
     
! Symmetrize Fock matrix contribution       
      do i = 1, ntbas(0)
        do j = 1, i
          average = (fock_pcm(i, j, 1) + fock_pcm(j, i, 1)) / 2.0d0
          fock_pcm(i, j, 1) = average    
          fock_pcm(j, i, 1) = average 
        end do
      end do
      
   end subroutine

   subroutine gauss_test(nr_points, nuc_pol_chg, ele_pol_chg)
      implicit none

#include "dcbgen.h"

! Passed variables   
      integer, intent(in)  :: nr_points
      real(8), intent(in)  :: nuc_pol_chg(nr_points), ele_pol_chg(nr_points)

! Local variables
      integer :: i
      real(8) :: eps, factor
      real(8) :: nuc_pol_tot, ele_pol_tot 
      real(8) :: nuc_tot_chg, ele_tot_chg, theo_nuc_pol, theo_ele_pol
      real(8) :: calc_diff, theo_diff

! Obtain total nuclear and electronic polarization charges
      nuc_pol_tot = 0.0d0
      ele_pol_tot = 0.0d0
      do i = 1, nr_points
        nuc_pol_tot = nuc_pol_tot + nuc_pol_chg(i)
        ele_pol_tot = ele_pol_tot + ele_pol_chg(i)
      end do

! Obtain total nuclear and electronic charges
      nuc_tot_chg = 0.0d0
      call get_nuclear_charge(nuc_tot_chg) 
      ele_tot_chg = kcharg - nuc_tot_chg

! Obtain solvent static permittivity
      call get_epsilon_static(eps)
       
! Perform test   
      factor = (eps - 1)/eps
      theo_nuc_pol = - factor * nuc_tot_chg
      theo_ele_pol = - factor * ele_tot_chg
      
      calc_diff = ele_pol_tot + nuc_pol_tot
      theo_diff = theo_ele_pol + theo_nuc_pol

! Print
      print *, "Calculated Nuclear Polarization = ", nuc_pol_tot
      print *, "Gauss's Law Nuclear Polarization = ", theo_nuc_pol
   
      print *, "Calculated Electronic Polarization = ", ele_pol_tot
      print *, "Gauss's Law Electronic Polarization = ", theo_ele_pol

      print *, "(EP + NP)_calc - (EP + NP)_Gauss", calc_diff - theo_diff 

   end subroutine

   subroutine get_nuclear_charge(nuc_tot_chg)
      implicit none
   
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h" 
      
      integer :: i
      real(8) :: nuc_tot_chg

! Obtain total nuclear charge of the molecule
      do i = 1, nucind
        nuc_tot_chg = nuc_tot_chg + fmult(istbnu(i)) * charge(i)
      end do

   end subroutine
