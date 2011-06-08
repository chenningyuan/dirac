!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!
! this module contains all the functionality required to add one- and two-electron
! spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis. 
! the spin-orbit corrections will be obtained as the spin-dependent parts of the 
! (converged) 4c-Fock operator in an atomic Dirac run.
!
! written by sknecht june 2011
!
module x2c_dirac_in_house_so_correction

  implicit none

  public prepare_so_correction_driver

  private

contains

!**********************************************************************
  subroutine prepare_so_correction_driver()
!**********************************************************************
!----------------------------------------------------------------------
     logical                            :: fndlab12
!**********************************************************************

  end subroutine prepare_so_correction_driver

!**********************************************************************

end module x2c_dirac_in_house_so_correction
