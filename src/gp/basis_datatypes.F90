! This module defines a data type for working with basis sets (an array of basis_func_info_t)
! As well as some functions to manipulate such a set

module basis_set_datatypes

        implicit none
        private

        public compress_unc_shell
!PARAMETERS:
!Generic:
        integer, parameter, private:: REALD=8 ! default real for this module
!DERIVED TYPES:
!       Basis function info:
        type, public:: basis_func_info_t
         integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        ! atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        ! atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        ! number of primitives in the contraction
         real(REALD):: coord(1:3)        ! coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type basis_func_info_t
!       Basis set info:
        type, public:: basis_set_info_t
         integer:: nshells=0             ! number of shells in the basis set
         integer:: nao=0                 ! number of basis functions
         integer:: basis_angular         ! 1=cartesian, 2=spherical
         type(basis_func_info_t), pointer:: gtos(:)
         integer, pointer :: shell_indices(:)
        end type basis_set_info_t

contains

       subroutine compress_unc_shell (shell,ierr)
!If basis functions were grouped together in the input, they appear as quasi-contracted functions
!with a lot of zero coefficients. Remove these zeroes and reduce the length to just one function.
!Routine will produce an error, if the function is a true contracted function (more than one non-zero coefficient)

          integer, intent(out)     :: ierr
          type(basis_func_info_t), intent(inout) :: shell

          integer i, j, non_zero
          real(8), allocatable :: exponents(:), coefficients(:)
          real(8), parameter   :: treshold=1.D-12

          allocate(exponents(shell%n_primitives))
          allocate(coefficients(shell%n_primitives))
          exponents    = shell%exponent
          coefficients = shell%coefficient

          ierr = 0
          non_zero = 0
          do i = 1, shell%n_primitives
             if (abs(coefficients(i)) > treshold) then
                non_zero = non_zero + 1
                j = i
             end if
          end do

          if (non_zero > 1) then
             ierr = 1
             return   ! leave shell untouched and return with error condition
          else if (non_zero == 0) then
             print*, " WARNING: compress_unc_shell found empty shell"
             ierr = 2
             j = 1    ! set j to the first function, to avoid complications with zero length arrays
          end if

          deallocate (shell%exponent)
          deallocate (shell%coefficient)
          allocate   (shell%exponent(1))
          allocate   (shell%coefficient(1))

          shell%n_primitives = 1
          shell%exponent     = exponents(j)
          shell%coefficient  = coefficients(j)

          deallocate (exponents)
          deallocate (coefficients)

       end subroutine compress_unc_shell

end module basis_set_datatypes
