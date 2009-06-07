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
