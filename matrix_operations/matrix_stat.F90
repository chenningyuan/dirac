module matrix_stat

  implicit none

  integer :: matrix_stat_allocations     = 0
  integer :: matrix_stat_deallocations   = 0
  integer :: matrix_stat_multiplications = 0

! private

contains

  subroutine print_matrix_stat()

    print '("@ " i10 a)', matrix_stat_allocations,     ' matrix allocations'
    print '("@ " i10 a)', matrix_stat_deallocations,   ' matrix deallocations'
    print '("@ " i10 a)', matrix_stat_multiplications, ' matrix multiplications'

  end subroutine

end module
