!This module defines the atomic orbital data types used in the ROSE program
       module rose_ao

        implicit none
        private

        public :: delete_adf_objects

        integer, parameter, private:: REALD=8 ! default real for this module

!       Basis function info:
        type, public:: gaussian_func_info_t
         integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        ! atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        ! atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        ! number of primitives in the contraction
         real(REALD):: coord(1:3)        ! coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type gaussian_func_info_t

!       Basis set info:
        type, public:: basis_set_info_t

         ! general
         integer:: basis_type=1          ! 1=gaussian, 2=slater
         integer:: nao=0                 ! number of basis functions

         ! specific for Gaussians (DIRAC, pyscf, psi4, ...)
         integer:: basis_angular=1       ! 1=cartesian, 2=spherical
         integer:: nshells=0             ! number of shells in the basis set
         type(gaussian_func_info_t), pointer:: gtos(:)
         integer, pointer :: shell_indices(:)

        end type basis_set_info_t

        contains

           subroutine delete_adf_objects (basis)
              type(basis_set_info_t), intent(inout) :: basis
           end subroutine delete_adf_objects

       end module rose_ao
