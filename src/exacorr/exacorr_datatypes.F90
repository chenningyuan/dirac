!This module defines the data types used in the ExaCorr module of DIRAC
       module exacorr_datatypes
        implicit none
        private
!PARAMETERS:
 !Generic:
        integer, parameter, private:: REALD=8 !default real for this module
!DERIVED TYPES:
 !Basis function info:
        type, public:: basis_func_info_t
         integer:: orb_momentum=-1       !orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        !atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        !atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        !number of primitives in the contraction
         real(REALD):: coord(1:3)        !coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type basis_func_info_t
 !Molecular spinor info:
        type, public:: mo_spinor_info_t
         integer:: symmetry=0            !molecular spinor irreducible representation (0: no symmetry)
        end type mo_spinor_info_t

       end module exacorr_datatypes
