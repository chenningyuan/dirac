module dft_ao_in_mem

  implicit none

  real(8), allocatable, target :: ao_in_mem(:)
  real(8), pointer             :: ao_in_mem_p(:)
  real(8), pointer             :: gao_p(:)
  real(8), pointer             :: gao1_p(:)

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
