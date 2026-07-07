!This module defines the atomic orbital data types used in the ROSE program
       module rose_ao

        use BasisTypeModule
        use FragGeometryTypeModule

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

         ! optional small-component Gaussian basis paired with the large-component basis
         logical :: has_small_component = .false.
         integer :: nao_small = 0
         integer :: nshells_small = 0
         integer :: basis_angular_small = 1
         type(gaussian_func_info_t), pointer :: gtos_small(:)
         integer, pointer :: shell_indices_small(:)

         ! specific for Slaters (ADF)
         ! we organize by fragments to be able to employ existing ADF types and subroutines
         ! more flexible would be to organize by individual shells like done for Gaussians
         type(BasisType) :: stos
         type(FragGeometryType) :: geom
         integer, allocatable :: ao_map(:)

        end type basis_set_info_t

        contains

           subroutine delete_adf_objects (basis)
              type(basis_set_info_t), intent(inout) :: basis
              integer :: i
              
              if (basis%has_small_component) then
                 do i = 1, basis%nshells_small
                    deallocate(basis%gtos_small(i)%exponent)
                    deallocate(basis%gtos_small(i)%coefficient)
                    nullify(basis%gtos_small(i)%exponent)
                    nullify(basis%gtos_small(i)%coefficient)
                 end do
                 deallocate(basis%gtos_small,basis%shell_indices_small)
                 nullify(basis%gtos_small)
                 nullify(basis%shell_indices_small)
                 basis%has_small_component = .false.
                 basis%nao_small = 0
                 basis%nshells_small = 0
                 basis%basis_angular_small = 1
              end if

              if (allocated(basis%ao_map)) deallocate(basis%ao_map)
              call Delete(basis%stos)
              call Delete(basis%geom)

           end subroutine delete_adf_objects

       end module rose_ao
