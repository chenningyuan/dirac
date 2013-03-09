!
! DIRAC-side interface routines for the Polarizable Continuum Model
! Roberto Di Remigio 2012
!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
! 
! The driver routine and the solver-related routines will be provided
! through a FORTRAN90 module. The cavity-related routines are called 
! by C++. Thus we put them outside the module, to avoid confusion due
! to name mangling.
!

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

   subroutine collect_nctot(nctot)
      implicit none
      integer, intent(out) :: nctot

#include "mxcent.h"
#include "nuclei.h"
      
      nctot = nucdep

   end subroutine
      
      
   subroutine collect_atoms(charges, centers)
      implicit none

#include "mxcent.h"
#include "nuclei.h"
      
      real(8), intent(out) :: charges(*)
      real(8), intent(out) :: centers(3,*)
      integer              :: units
      
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
   
   end subroutine

module pcm_interface
   
   use electrostatic_potential
   use pcm_write
   use pcmmod_cfg

   implicit none 
   
   public pcm_energy_driver   
   public pcm_fock_driver
   public pcm_lintra_driver
   public get_pcm_energy

   private

   real(8), allocatable :: charges(:)
   real(8), allocatable :: potentials(:)
   real(8), allocatable :: tess_cent(:, :)
   real(8)              :: pcm_energy
! nr_points is 32 bit to avoid mess in passing it from C++ to FORTRAN
   integer(4)           :: nr_points = -1

   contains 

   subroutine pcm_initialize(dmat, ndmat, work, lwork, do_lintra)
!
! This subroutine sets all the global private variables.
! 

! Passed variables
      real(8) :: dmat(*)
      real(8) :: work(*)
      integer :: ndmat
      integer :: lwork
      logical :: do_lintra

      call get_cavity_size(nr_points)
  
! Allocation   
      allocate(charges(nr_points))
      charges = 0.0d0
      allocate(potentials(nr_points))
      potentials = 0.0d0
      allocate(tess_cent(3, nr_points))
      tess_cent = 0.0d0
! Get tesserae centers 
      call get_tess_centers(tess_cent)
! Compute total polarization charges
      if (do_lintra) then
      ! This is the case where we are calculating the linear response
         call pcm_lintra_charges(dmat, ndmat, work, lwork, pcmmod_print)
      else
      ! This is the case where we are calculating the energy  
         call pcm_charges(dmat, ndmat, work, lwork, pcmmod_print)
      endif
      
   end subroutine
   
   subroutine pcm_finalize()
!
! This subroutine deallocates all the global private variables.
! 
      nr_points = -1 
      
! Deallocation   
      deallocate(charges)
      deallocate(potentials)
      deallocate(tess_cent)
      
   end subroutine 
   
   subroutine pcm_energy_driver(pol_ene, dmat, ndmat, work, lwork)
!
! Polarizable Continuum Model:
! driver routine for polarization energy calculation, U_pol.
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
!
! Roberto Di Remigio, 2012
!

! Passed variables
      real(8) :: pol_ene
      real(8) :: dmat(*)
      real(8) :: work(*)
      integer :: ndmat
      integer :: lwork

      if (nr_points <= 0) call pcm_initialize(dmat, ndmat, work, lwork, .false.)
! Compute polarization energy
      if (.not.pcmmod_separate) then
         call comp_pol_ene_pcm(pol_ene, 1)
      else
         call comp_pol_ene_pcm(pol_ene, 0)
      end if
      pcm_energy = pol_ene

   end subroutine
   
   subroutine pcm_fock_driver(fock, dmat, ndmat, factor, work, lwork)
!
! Polarizable Continuum Model.
! Driver routine for Fock matrix contributions, j and X(D).
!
! Roberto Di Remigio, 2012
! 

! Passed variables
      real(8)              :: fock(*)
      real(8)              :: dmat(*)
      real(8)              :: factor
      real(8)              :: work(*)
      integer              :: ndmat
      integer              :: lwork
      
      if (nr_points <= 0) call pcm_initialize(dmat, ndmat, work, lwork, .false.)
! Update Fock matrix with PCM contribution
      call pcm_fock(fock, factor, work, lwork)
! Free temporary allocated memory
      call pcm_finalize

   end subroutine

   subroutine pcm_lintra_driver(fock, pert_dmat, ndmat, factor, work, lwork)
!
! Polarizable Continuum Model.
! Driver routine for "explicit" contributions to linear response.
! pert_dmat is the b-transformed AO density matrix.
!
! Roberto Di Remigio, 2012
! 

! Passed variables
      real(8)              :: fock(*)
      real(8)              :: pert_dmat(*)
      real(8)              :: factor
      real(8)              :: work(*)
      integer              :: ndmat
      integer              :: lwork

! Local variables

! Calculate one-index transformed potentials and charges (ONLY ELECTRONIC!)
      if (nr_points <= 0) call pcm_initialize(pert_dmat, ndmat, work, lwork, .true.)
! Update Fock matrix with PCM contribution:
! contraction of unperturbed potentials and perturbed electronic charges.
      call pcm_fock(fock, factor, work, lwork)
! Free temporary allocated memory
      call pcm_finalize
       
   end subroutine

   function get_pcm_energy()
      
      real(8) get_pcm_energy

! pcm_energy is the polarization energy:
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
      if (.not.pcmmod_separate) then
         call comp_pol_ene_pcm(pcm_energy, 1)
      else 
         call comp_pol_ene_pcm(pcm_energy,0)
      end if

      get_pcm_energy = pcm_energy

   end function

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

   subroutine pcm_charges(dmat, ndmat, work, lwork, printlevel)

! Passed variables
      integer, intent(in)  :: printlevel
      real(8), intent(in)  :: dmat(*)
      real(8)              :: work(*)
      integer              :: ndmat
      integer              :: lwork

! Local variables
      real(8), allocatable :: nuc_pot(:), nuc_pol_chg(:) ! Allocated only if pcmmod_separate is .true.
      real(8), allocatable :: ele_pot(:), ele_pol_chg(:) ! Allocated only if pcmmod_separate is .true.
      integer              :: i
      real(8)              :: factor = 1.0d0
      character(7)         :: potName, potName1, potName2
      character(7)         :: chgName, chgName1, chgName2
      
      if (.not.(pcmmod_separate)) then 
! Compute MEP @tesserae centers
         potName = 'TotPot'//char(0)
         chgName = 'TotChg'//char(0)
         call get_mep(nr_points, tess_cent, potentials, dmat, work, lwork, 0)
! Set a cavity surface function with the MEP
         call set_surface_function(nr_points, potentials, potName)
! Compute polarization charges and set the proper surface function
         call comp_chg_pcm(potName, chgName)
! Get polarization charges @tesserae centers
         call get_surface_function(nr_points, charges, chgName)

! Print some information
         if (printlevel > 5) then
            do i = 1, nr_points
              print *, "MEP @point", i, potentials(i)
              if (printlevel > 10) then
                 print *, "PC @point", i, charges(i)
              end if
            end do
         end if

! Write to file MEP and ASC
         call pcm_write_file(nr_points, potentials, charges)
      else
! Allocation
         allocate(nuc_pot(nr_points))
         nuc_pot = 0.0d0
         allocate(nuc_pol_chg(nr_points))
         nuc_pol_chg = 0.0d0
         allocate(ele_pot(nr_points))
         ele_pot = 0.0d0
         allocate(ele_pol_chg(nr_points))
         ele_pol_chg = 0.0d0
      
         call pcm_nuclear_charges(nuc_pot, nuc_pol_chg)
         potName  = 'ElePot'//char(0)
         chgName  = 'EleChg'//char(0)
         call pcm_electronic_charges(dmat, ndmat, ele_pot, ele_pol_chg, work, lwork, potName, chgName)

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
! NOT SURE IF NECESSARY
        if (printlevel > 4 ) then 
           call gauss_test(nuc_pol_chg, ele_pol_chg)
        end if

! Obtain vector of total polarization charges 
        chgName  = 'TotChg'//CHAR(0)
        chgName1 = 'NucChg'//CHAR(0)
        chgName2 = 'EleChg'//CHAR(0)

        call append_surf_func(chgName)
        call clear_surf_func(chgName)
        call add_surface_function(chgName, factor, chgName1)
        call add_surface_function(chgName, factor, chgName2)
        call get_surface_function(nr_points, charges, chgName)

! Obtain vector of total MEP
! THIS IS NOT STRICTLY NECESSARY AND CAN BE AVOIDED
        potName  = 'TotPot'//char(0)
        potName1 = 'NucPot'//char(0)
        potName2 = 'ElePot'//char(0)
        call append_surf_func(potName)
        call clear_surf_func(potName)
        call add_surface_function(potName, factor, potName1)
        call add_surface_function(potName, factor, potName2)
        call get_surface_function(nr_points, potentials, potName)

! Write to file MEP and ASC
        call pcm_write_file(nr_points, potentials, charges)
 
        deallocate(nuc_pot)
        deallocate(nuc_pol_chg)
        deallocate(ele_pot)
        deallocate(ele_pol_chg)
      end if

   end subroutine
   
   subroutine pcm_lintra_charges(dmat, ndmat, work, lwork, printlevel)
! Calculation of b-transformed electronic polarization charges.
! This is needed to get the "explicit" PCM contributions to the
! linear response equations.

! Passed variables
      integer, intent(in)  :: printlevel
      real(8), intent(in)  :: dmat(*)
      real(8)              :: work(*)
      integer              :: ndmat
      integer              :: lwork

! Local variables
      character(13)        :: potName, chgName
      integer              :: i
      
! Update Fock matrix with PCM contribution:
! contraction of unperturbed potentials and perturbed electronic charges.
      potName = 'EleLinTraPot'//char(0)
      chgName = 'EleLinTraChg'//char(0)
      call pcm_electronic_charges(dmat, ndmat, potentials, charges, work, lwork, potName, chgName)
!TROND
!      WRITE(6,*) 'Perturbed electronic charges..'
!      do i = 1,nr_points
!        WRITE(6,*) i,lintra_ele_chg(i)
!      enddo
      
   end subroutine

   subroutine pcm_nuclear_charges(nuc_pot, nuc_pol_chg)

! Passed variables      
      real(8)      :: nuc_pot(nr_points)
      real(8)      :: nuc_pol_chg(nr_points)

! Local variables
      character(7) :: potName
      character(7) :: chgName
         
! Compute nuclear MEP @tesserae centers
      potName = 'NucPot'//char(0)
      chgName = 'NucChg'//char(0)
      call get_nuclear_mep(nr_points, tess_cent, nuc_pot, 0)
! Set a cavity surface function with the nuclear MEP
      call set_surface_function(nr_points, nuc_pot, potName)
! Compute nuclear polarization charges and set the proper surface function
      call comp_chg_pcm(potName, chgName)
! Get nuclear polarization charges @tesserae centers
! THIS IS NOT STRICTLY NECESSARY AND CAN BE AVOIDED
      call get_surface_function(nr_points, nuc_pol_chg, chgName)

   end subroutine

   subroutine pcm_electronic_charges(dmat, ndmat, ele_pot, ele_pol_chg, work, lwork, potName, chgName)

! Passed variables
      real(8)      :: dmat(*)                
      real(8)      :: ele_pot(nr_points)
      real(8)      :: ele_pol_chg(nr_points)
      real(8)      :: work(*)
      integer      :: ndmat
      integer      :: lwork
      character(7) :: potName
      character(7) :: chgName

! Local variables
      integer      :: i

! Compute electronic MEP @tesserae centers
      do i = 1, ndmat
        call get_electronic_mep(nr_points, dmat(i), tess_cent, ele_pot,    &
                         work, lwork, 0, .false.)
      enddo
! Set a cavity surface function with the electronic MEP
      call set_surface_function(nr_points, ele_pot, potName)
! Compute electronic polarization charges and set the proper surface function
      call comp_chg_pcm(potName, chgName)
! Get electronic polarization charges @tesserae centers
! THIS IS NOT STRICTLY NECESSARY AND CAN BE AVOIDED
      call get_surface_function(nr_points, ele_pol_chg, chgName)
 
   end subroutine
   
   subroutine pcm_fock(fock, factor, work, lwork)
! This subroutine is completely general, it is used also
! in the generation of the linear response terms.
#include "dcbbas.h"
#include "dgroup.h"

! Passed variables
      real(8), intent(out) :: fock(*)
      real(8)              :: factor
      real(8)              :: work(*)
      integer              :: lwork

! Local variables
      real(8), allocatable :: fock_pcm(:, :, :)
      real(8)              :: average
      integer              :: i, j
      logical              :: get_matrix = .true.

! Allocation
      allocate(fock_pcm(ntbas(0), ntbas(0), nz))
      fock_pcm = 0.0d0
! Contraction of electrostatic potential integrals with polarization charges
      call get_electronic_mep(nr_points, fock_pcm, tess_cent, charges,    &
                        work, lwork, 0, get_matrix)
     
! Symmetrize Fock matrix contribution       
      do i = 1, ntbas(0)
        do j = 1, i
          average = 0.5d0 * (fock_pcm(i, j, 1) + fock_pcm(j, i, 1))
          fock_pcm(i, j, 1) = average    
          fock_pcm(j, i, 1) = average
        end do
      end do
       
!     print *, "FOCK_PCM"
!     call prqmat(fock_pcm, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                      nz, ipqtoq(1,0), 6)
!     print *, "FOCK BEFORE PCM"
!     call prqmat(fock, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                      nz, ipqtoq(1,0), 6)

! Add fock_pcm to Fock matrix
      call daxpy(n2bbasx, factor, fock_pcm, 1, fock, 1)

!     print *, "FOCK AFTER PCM"
!     call prqmat(fock, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                      nz, ipqtoq(1,0), 6)

      deallocate(fock_pcm)

   end subroutine

   subroutine gauss_test(nuc_pol_chg, ele_pol_chg)

#include "dcbgen.h"

! Passed variables   
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

end module
