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

module dft_ao_in_mem

  implicit none

  integer, parameter :: max_nr_of_batches = 1000
  integer            :: batch_pointer(max_nr_of_batches)     = 1
  logical            :: batch_is_in_mem(max_nr_of_batches)   = .false.
  logical            :: batch_fits_in_mem(max_nr_of_batches) = .true.
  logical            :: ao_mem_allocated                     = .false.

  real(8), allocatable, target :: ao_in_mem(:)

! public init_ao_in_mem
! public kill_ao_in_mem

! private

contains

  subroutine init_ao_in_mem(n)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: n
!   ----------------------------------------------------------------------------

    allocate(ao_in_mem(n))

  end subroutine

  subroutine kill_ao_in_mem()

    deallocate(ao_in_mem)

  end subroutine

end module
