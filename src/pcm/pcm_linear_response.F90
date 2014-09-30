module pcm_linear_response

   use, intrinsic :: iso_c_binding

   implicit none

   public pcm_linear_response_initialize
   public pcm_linear_response_finalize

   public pcm_oper_ao_driver
   public compute_mep_asc
   public pcm_oit_asc

   private

! If false the interface will refuse to be accessed
   logical                     :: is_initialized = .false.
   real(c_double), allocatable :: tess_cent(:, :)
   integer(c_int)              :: nr_points
   integer(c_int)              :: nr_points_irr
! A (maybe clumsy) way of passing LUPRI without using common blocks   
   integer                     :: global_print_unit
! A counter for the number of LR iterations
   integer, save               :: LR_iteration_counter = 1

   contains

      subroutine pcm_linear_response_initialize(print_unit)
      
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
              
      is_initialized = .true.

      end subroutine pcm_linear_response_initialize

      subroutine pcm_linear_response_finalize()

      if (.not. is_initialized) then
         print *, 'Error: pcm_linear_response was never initialized.'
         stop 1
      end if
! Free the memory taken from the free store both in Fortran and in C++
      deallocate(tess_cent)

      call tear_down_pcm
                                                                    
      is_initialized = .false.
                                                                    
      end subroutine pcm_linear_response_finalize

      subroutine check_if_interface_is_initialized()

      if (.not. is_initialized) then
         print *, 'Error: pcm_linear_response is not initialized.'
         stop 1
      end if

      end subroutine

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
      
      end subroutine pcm_oper_ao_driver

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
      
      allocate(mep(nr_points))
      mep = 0.0d0
      allocate(asc(nr_points))
      asc = 0.0d0
      ! The irreducible representation, for the moment the totally
      ! symmetric one
      irrep = 0

      if (.not.(pcmmod_separate)) then
        potName = 'TotMEP'//c_null_char 
        chgName = 'TotASC'//c_null_char
! Calculate the (total) Molecular Electrostatic Potential
        call get_mep(nr_points, nr_points_irr, tess_cent, mep,        &
                     density_matrix, work, lfree, 0, pcmmod_skipss)
! Set a cavity surface function with the MEP
        call set_surface_function(nr_points, mep, potName)
! Compute polarization charges and set the proper surface function
        call compute_asc(potName, chgName, irrep)
! Get polarization charges @tesserae centers
        call get_surface_function(nr_points, asc, chgName)

! Print some information
        if (pcmmod_print > 5) then
          write(global_print_unit, '(A, T27, A, T62, A)') "Finite element #", "Total MEP", "Total ASC"
          do i = 1, nr_points
            write(global_print_unit, '(I6, 2(20X, F15.12))') i, mep(i), asc(i)
          end do
        end if

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
        call get_nuclear_mep(nr_points, nr_points_irr, tess_cent, nuc_pot, 0)
        call set_surface_function(nr_points, nuc_pot, potName1)
        call compute_asc(potName1, chgName1, irrep)
        call get_surface_function(nr_points, nuc_pol_chg, chgName1)

        potName2 = 'EleMEP'//c_null_char
        chgName2 = 'EleASC'//c_null_char
        call get_electronic_mep(nr_points, nr_points_irr, density_matrix, tess_cent,  &
                                ele_pot, work, lfree, 0, .false., pcmmod_skipss)
        call set_surface_function(nr_points, ele_pot, potName2)
        call compute_asc(potName2, chgName2, irrep)
        call get_surface_function(nr_points, ele_pol_chg, chgName2)

! Print some information
        if (pcmmod_print > 5) then
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

      subroutine pcm_oit_asc(mep_oit_asc, perturbed_density_matrix, nr_density_matrices, work, lwork)
!
! Polarizable Continuum Model.
! Driver routine for "explicit" contributions to linear response.
! perturbed_density_matrix is the b-transformed AO density matrix.
! oit stands for One-Index Transformed
!
! Roberto Di Remigio, 2012-2013
! 

      use pcm_integrals, only: get_electronic_mep
      use pcmmod_cfg

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

! Passed variables
      real(c_double) :: mep_oit_asc(ntbas(0), ntbas(0), nz)
      real(8)        :: perturbed_density_matrix(*)
      integer        :: nr_density_matrices
      real(8)        :: work(*)
      integer        :: lwork
! Local variables
      real(c_double), allocatable :: mep_oit(:), asc_oit(:)
      character(8)                :: potName, chgName
      real(8)                     :: average
      integer                     :: i, j, irrep

      call check_if_interface_is_initialized

      allocate(mep_oit(nr_points))
      mep_oit = 0.0d0
      potName = 'MEP_OIT'//c_null_char
      allocate(asc_oit(nr_points))
      asc_oit = 0.0d0
      chgName = 'ASC_OIT'//c_null_char
      ! The irreducible representation, for the moment the totally
      ! symmetric one
      irrep = 0
! Compute electronic MEP @tesserae centers
      do i = 1, nr_density_matrices
        call get_electronic_mep(nr_points, nr_points_irr, perturbed_density_matrix(i), tess_cent, & 
                                mep_oit, work, lwork, 0, .false., pcmmod_skipss)
      enddo
      call set_surface_function(nr_points, mep_oit, potName)
! Compute polarization charges and set the proper surface function
      call compute_asc(potName, chgName, irrep)
      call get_surface_function(nr_points, asc_oit, chgName)
! Now contract unpertubed potential matrix elements with perturbed ASC
      call get_electronic_mep(nr_points, nr_points_irr, mep_oit_asc, tess_cent,  & 
                              asc_oit, work, lwork, 0, .true., pcmmod_skipss)
      do i = 1, ntbas(0)
        do j = 1, i
          average = 0.5d0 * (mep_oit_asc(i, j, 1) + mep_oit_asc(j, i, 1))
          mep_oit_asc(i, j, 1) = average    
          mep_oit_asc(j, i, 1) = average
        end do
      end do

      deallocate(mep_oit)
      deallocate(asc_oit)

! Update the LR iteration counter
      LR_iteration_counter = LR_iteration_counter + 1

      end subroutine

end module pcm_linear_response
