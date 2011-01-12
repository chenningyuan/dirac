module atomic_x2c_cfg

! stefan: this module contains all required internal basis set, offset and 
!         pointer informations for the atom-by-atom X2C approach.

  implicit none

  save

! parameters

!                               see P. Pyykkoe, Phys. Chem. Chem. Phys. 13, 161-168 (2011).
  integer, parameter, public :: max_number_of_indiv_centers = 172 

! character block

! logical block
  logical, public :: atomic_x2c_complete_mol                =  .false.
  logical, public :: atomic_x2c_add_amfi                    =  .false.

! double precision block

! integer block

! calculated integer arrays for SA-AO/orthonormal MO basis

end module
