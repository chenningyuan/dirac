!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end
!
module atomic_x2c

! stefan: this module contains all the functionality to carry out an ensuing 
!         transformation
!         ========================
!         X_2c^++ = U^+ * X_4c * U
!         ========================
!         where U was obtained by means of the atom-by-atom X2C approach:
!         U =~ SUM_{all unique atomic center/fragments i} U_i = U^{@}
!
!         for more details see, e.g. W. Liu, Molecular Physics, 108, (2010) p. 1679ff. 
!         here we implement scheme (3) and (4) which are described in more detail on page 1702f.
!
! written by sknecht jan 2011 - ???
!

  use atomic_x2c_cfg
  use memory_allocator  

  implicit none

  public x2c_atomic_x2c_driver

  private

! parameters
  integer, parameter :: bla_bla = 1 

! types
  type fragments_pctmat
    real(8), allocatable :: pctmat_f_i(:)
    integer, allocatable :: naosh_all_f_i(:)
    integer, allocatable :: naosh_L_f_i(:)
  end type

contains

!***********************************************************************
  subroutine x2c_atomic_x2c_driver()

!**********************************************************************
!
!    purpose: 
!
!----------------------------------------------------------------------
  end subroutine x2c_atomic_x2c_driver

!***********************************************************************

end module
