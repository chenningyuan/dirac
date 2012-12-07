!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
  logical, public :: fx2c_complete_mol                =  .false.

! double precision block

! integer block

! types
  type fx2c_pctmat_info
    real(8), allocatable :: fx2c_pctmat(:)
    integer, allocatable :: fx2c_naosh_all(:)
    integer, allocatable :: fx2c_naosh_L(:)
  end type

end module
