!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

!
! DIRAC-side interface routines for the Polarizable Continuum Model
! Roberto Di Remigio 2012
!
! We divide the interface routines into:
!   1. cavity-related routines;
!   2. solver-related routines.
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


module pcm_scf

   use, intrinsic :: iso_c_binding
   
   implicit none 

   public pcm_scf_initialize
   public pcm_scf_finalize
  
   public pcm_energy_driver   
   public pcm_oper_ao_driver
   public get_pcm_energy

   private

! If false the interface will refuse to be accessed
   logical                     :: is_initialized = .false.
   real(c_double), allocatable :: tess_cent(:, :)
   real(c_double)              :: pcm_energy
   integer(c_int)              :: nr_points
   integer(c_int)              :: nr_points_irr
! A (maybe clumsy) way of passing LUPRI without using common blocks   
   integer                     :: global_print_unit
! A counter for the number of SCF iterations
   integer, save               :: scf_iteration_counter = 1

   contains

      subroutine pcm_scf_initialize(print_unit)

      use pcmmod_cfg, only: pcmmod_host_provides_input
      
      integer, intent(in) :: print_unit
      
      integer :: host_provides_input = 0

      global_print_unit = print_unit

      if (pcmmod_host_provides_input) host_provides_input = 1 

      call set_up_pcm(host_provides_input)
      call print_pcm

      call get_cavity_size(nr_points, nr_points_irr)

      allocate(tess_cent(3, nr_points))
      tess_cent = 0.0d0
      call get_tesserae(tess_cent)

      pcm_energy = 0.0d0
              
      is_initialized = .true.
                                                                 
      end subroutine pcm_scf_initialize 
      
      subroutine pcm_scf_finalize()

      if (.not. is_initialized) then
         write(global_print_unit, *) 'pcm_scf has already been finalized.'
      else
! Free the memory taken from the free store both in Fortran and in C++
         deallocate(tess_cent)                                          
                                                                        
         call tear_down_pcm
                                                                       
         is_initialized = .false.
      end if
                                                                    
      end subroutine
                                                                    
      subroutine check_if_interface_is_initialized()

      if (.not. is_initialized) then
         write(global_print_unit, *) 'Error: pcm_scf is not initialized.'
         stop 1
      end if

      end subroutine

!      subroutine collect_nctot(nr_nuclei) bind(c, name='collect_nctot') 
!
!#include "mxcent.h"
!#include "nuclei.h"
!     
!      integer(c_int), intent(out) :: nr_nuclei 
!                                               
!      nr_nuclei = nucdep
!                                               
!      end subroutine collect_nctot
!     
!      subroutine collect_atoms(atomic_charges, atomic_centers) bind(c, name='collect_atoms')
!
!#include "mxcent.h"
!#include "nuclei.h"
!     
!      real(c_double), intent(out) :: atomic_charges(*)   
!      real(c_double), intent(out) :: atomic_centers(3,*)
!                                                         
!      integer :: i, j, k 
!     
!! Get coordinates
!      call getacord(atomic_centers)
!! Get charges      
!      i = 0                             
!      do j = 1, nucind
!        do k = 1, nucdeg(j)
!          i = i + 1
!          atomic_charges(i) = charge(j)
!        enddo
!      enddo
!                                        
!      end subroutine collect_atoms
!      
!      subroutine set_point_group(nr_gen, gen1, gen2, gen3) bind(c, name='set_point_group')
!
!#include "mxcent.h"                     
!#include "maxorb.h"
!#include "maxaqn.h"
!#include "symmet.h"
!
!      integer(c_int), intent(inout) :: nr_gen 
!      integer(c_int), intent(inout) :: gen1, gen2, gen3
!
!      nr_gen = 0 !pcm_igen(1)   
!      gen1   = 0 !pcm_igen(2) 
!      gen2   = 0 !pcm_igen(3)
!      gen3   = 0 !pcm_igen(4)
!
!      end subroutine set_point_group

      subroutine pcm_energy_driver(density_matrix, pol_ene, work, lfree)
!
! Polarizable Continuum Model:
! driver routine for polarization energy calculation, U_pol.
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
!
! Roberto Di Remigio, 2012
!
      real(8)        :: density_matrix(*)
      real(c_double) :: pol_ene
      real(8)        :: work(*)
      integer        :: lfree

! Make sure that the interface is initialized first
      call check_if_interface_is_initialized

! OK, now compute MEP and ASC
      call compute_mep_asc(density_matrix, work, lfree)

! pcm_energy is the polarization energy:
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
      call compute_polarization_energy(pol_ene)

! Now make the value of the polarization energy known throughout the module
      pcm_energy = pol_ene

      end subroutine pcm_energy_driver
   
      subroutine pcm_oper_ao_driver(oper, charge_name, work, lwork)
!
! Polarizable Continuum Model.
! Driver routine for Fock matrix contributions, J and X(D).
!
! Roberto Di Remigio, 2012
!
      use pcm_integrals, only: get_electronic_mep 
      use pcmmod_cfg

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

! Passed variables
! oper(*) is just the PCM contribution in AO basis
      real(8), intent(out)     :: oper(ntbas(0), ntbas(0), nz)
      character(*), intent(in) :: charge_name
      real(8)                  :: work(*)
      integer                  :: lwork
! Local variables
      real(c_double), allocatable :: asc(:)
      real(8)                     :: average
      real(8), allocatable        :: oper_mo(:, :, :)
      integer                     :: i, j

      call check_if_interface_is_initialized 

! Here we need the ASC, we can either get it from C++ (as it 
! has been saved as a surface function) or from Fortran asking
! for the mep and asc to be recomputed.
! RDR, 220813 !WARNING! what happens in the second case if we use this
! subroutine as a driver also for the linear response part??
      allocate(asc(nr_points))
      asc = 0.0d0
      call get_surface_function(nr_points, asc, charge_name)
! Contraction of electrostatic potential integrals with polarization charges
      call get_electronic_mep(nr_points, nr_points_irr, oper, tess_cent,   &
                              asc, work, lwork, 0, .true., pcmmod_skipss)
! Symmetrize Fock matrix contribution       
      do i = 1, ntbas(0)
        do j = 1, i
          average = 0.5d0 * (oper(i, j, 1) + oper(j, i, 1))
          oper(i, j, 1) = average    
          oper(j, i, 1) = average
        end do
      end do

      if (pcmmod_dospf) then
        allocate(oper_MO(norb(1), norb(1), nz))
        oper_MO = 0.0d0
        print *, "FOCK_PCM before elimination of spin-dependent terms."
        call prqmat(oper, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                    nz, ipqtoq(1,0), 6)
        call spfao(oper_MO, oper, work, lwork)
        print *, "FOCK_PCM after elimination of spin-dependent terms."
        call prqmat(oper, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                    nz, ipqtoq(1,0), 6)
        print *, "FOCK_PCM_MO after elimination of spin-dependent terms."
        call prqmat(oper_mo, norb(1), norb(1), norb(1), norb(1), &
                    nz, ipqtoq(1,0), 6)
        deallocate(oper_MO)
      end if

      deallocate(asc)
! Update the scf_iteration_counter
      scf_iteration_counter = scf_iteration_counter + 1
      
      end subroutine pcm_oper_ao_driver

      real(8) function get_pcm_energy()
      
! pcm_energy is the polarization energy:
! U_pol = 0.5 * (U_NN + U_Ne + U_eN + U_ee)
! it has already been calculated and "stored" in the module

      get_pcm_energy = pcm_energy

      end function

      subroutine compute_mep_asc(density_matrix, work, lfree)
!
! Calculate the molecular electrostatic potential and
! the apparent surface charge at the cavity points.
!
! The user can control via the DIRAC input the following:
!    * switch between separate and total evaluation of the
!      nuclear and electronic parts;
!
      use pcm_integrals, only: get_nuclear_mep, get_electronic_mep, get_mep
      use pcm_write, only: pcm_write_file, pcm_write_file_separate
      use pcmmod_cfg
      use pcm_gp, only: timing

      real(8), intent(in)    :: density_matrix(*)
      real(8), intent(inout) :: work(*)
      integer                :: lfree
! Local variables
      real(c_double), allocatable :: mep(:)
      real(c_double), allocatable :: asc(:)
      real(c_double), allocatable :: nuc_pot(:), nuc_pol_chg(:)
      real(c_double), allocatable :: ele_pot(:), ele_pol_chg(:)
      character(7)                :: potName, chgName
      character(7)                :: potName1, chgName1, potName2, chgName2
      integer                     :: i, irrep
      real(8)                     :: start_mep, finish_mep
      real(8)                     :: start_asc, finish_asc
      
      allocate(mep(nr_points))
      mep = 0.0d0
      allocate(asc(nr_points))
      asc = 0.0d0
      ! Ask for the totally symmetric irrep
      irrep = 0

      ! Zero out timing variables
      start_mep = 0.0d0
      finish_mep = 0.0d0
      start_asc = 0.0d0
      finish_asc = 0.0d0

      if (.not.(pcmmod_separate)) then
        potName = 'TotMEP'//c_null_char 
        chgName = 'TotASC'//c_null_char
! Calculate the (total) Molecular Electrostatic Potential
        start_mep = timing()
        call get_mep(nr_points, nr_points_irr, tess_cent, mep,      &
                     density_matrix, work, lfree, 0, pcmmod_skipss)
        finish_mep = timing()
! Calculation of ASC proceeds in three steps:
! a. set a cavity surface function with the MEP
! b. compute the actual ASC and set the proper surface function
! c. retrieve the ASC @tesserae centers
        start_asc = timing() 
        call set_surface_function(nr_points, mep, potName)
        call compute_asc(potName, chgName, irrep)
        call get_surface_function(nr_points, asc, chgName)
        finish_asc = timing() 

! Print some information
        if (pcmmod_print > 5) then
          write(global_print_unit, '(20X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
          write(global_print_unit, '(A, T27, A, T62, A)') "Finite element #", "Total MEP", "Total ASC"
          do i = 1, nr_points
            write(global_print_unit, '(I6, 2(20X, F15.12))') i, mep(i), asc(i)
          end do
        end if
! Print timings
        write(global_print_unit, '(a, f10.8, a)'), "* MEP evaluation (CPU): ", finish_mep-start_mep, "s"
        write(global_print_unit, '(a, f10.8, a)'), "* ASC evaluation (CPU): ", finish_asc-start_asc, "s"

! Write to file MEP and ASC
        call pcm_write_file(nr_points, mep, asc)
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

        potName1 = 'NucMEP'//c_null_char
        chgName1 = 'NucASC'//c_null_char
        start_mep = timing() 
        call get_nuclear_mep(nr_points, nr_points_irr, tess_cent, nuc_pot, 0)
        finish_mep = timing()

        start_asc = timing()
        call set_surface_function(nr_points, nuc_pot, potName1)
        call compute_asc(potName1, chgName1, irrep)
        call get_surface_function(nr_points, nuc_pol_chg, chgName1)
        finish_asc = timing()
! Print timings
        write(global_print_unit, '(a, f10.8, a)'), "* NucMEP evaluation (CPU): ", finish_mep-start_mep, "s"
        write(global_print_unit, '(a, f10.8, a)'), "* NucASC evaluation (CPU): ", finish_asc-start_asc, "s"
! Reset timing variables        
        start_mep = 0.0d0
        finish_mep = 0.0d0
        start_asc = 0.0d0
        finish_asc = 0.0d0

        potName2 = 'EleMEP'//c_null_char
        chgName2 = 'EleASC'//c_null_char
        start_mep = timing() 
        call get_electronic_mep(nr_points, nr_points_irr, density_matrix, tess_cent,  &
                                ele_pot, work, lfree, 0, .false., pcmmod_skipss)
        finish_mep = timing()
        
        start_asc = timing()
        call set_surface_function(nr_points, ele_pot, potName2)
        call compute_asc(potName2, chgName2, irrep)
        call get_surface_function(nr_points, ele_pol_chg, chgName2)
        finish_asc = timing()
! Print timings
        write(global_print_unit, '(a, f10.8, a)'), "* EleMEP evaluation (CPU): ", finish_mep-start_mep, "s"
        write(global_print_unit, '(a, f10.8, a)'), "* EleASC evaluation (CPU): ", finish_asc-start_asc, "s"

! Print some information
        if (pcmmod_print > 5) then
          write(global_print_unit, '(60X, A, 6X, I6)') "MEP and ASC at iteration", scf_iteration_counter
          write(global_print_unit, '(A, T27, A, T62, A, T97, A, T132, A)') "Finite element #", &
          "Nuclear MEP", "Nuclear ASC", "Electronic MEP", "Electronic ASC"
          do i = 1, nr_points
            write(global_print_unit, '(I6, 4(20X, F15.12))') i, nuc_pot(i), nuc_pol_chg(i), ele_pot(i), ele_pol_chg(i)
          end do
        end if

! Obtain vector of total MEP
        potName  = 'TotMEP'//c_null_char
        mep(:) = nuc_pot(:) + ele_pot(:)
        call set_surface_function(nr_points, mep, potName)

! Obtain vector of total polarization charges 
        chgName  = 'TotASC'//c_null_char
        asc(:) = nuc_pol_chg(:) + ele_pol_chg(:)
        call set_surface_function(nr_points, asc, chgName)

! Write to file MEP and ASC
        call pcm_write_file_separate(nr_points, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)
 
        deallocate(nuc_pot)
        deallocate(nuc_pol_chg)
        deallocate(ele_pot)
        deallocate(ele_pol_chg)
      end if
      deallocate(mep)
      deallocate(asc)

      end subroutine compute_mep_asc

end module pcm_scf
