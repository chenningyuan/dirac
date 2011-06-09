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

  use memory_allocator 
  use x2c_utils

  implicit none

  public prepare_so_correction_driver

  private

contains

!**********************************************************************
  subroutine prepare_so_correction_driver(nr_fsym,    &
                                          norb_dim,   &
                                          print_lvl)
!**********************************************************************
     integer, intent(in)                :: nr_fsym
     integer, intent(in)                :: norb_dim(nr_fsym)
     integer, intent(in)                :: print_lvl
!----------------------------------------------------------------------
     integer                            :: i
     integer                            :: j
     integer, allocatable               :: boson_info(:)
!**********************************************************************

      j = 0
      do i = 1, nr_fsym
        j = j + norb_dim(i)
      end do
      call alloc(boson_info,j)
     
      call get_boson_irrep_info(boson_info,           &
                                nr_fsym,              &
                                norb_dim,             &
                                print_lvl)

      call dealloc(boson_info)

!     note: CALL OPNFIL(LUCMOS,'DFCMOS','UNKNOWN','PREDHF') for MOs in onMO basis
!           otherwise see PRECMOS for a template to transform to onMO basis
!     note: use APPBOS to get the correct boson irrep assignment of each spinor.
!     note: if boson_info is ready take out the spin-free part from the pe-part of the fock matrix f(pe,pe,NZ) 
!           in onMO by eliminating the symmetric (boson_irrep(i) == boson_irrep(j)) part of 
!           fock(i,j,1).

  end subroutine prepare_so_correction_driver

!**********************************************************************

end module x2c_dirac_in_house_so_correction
