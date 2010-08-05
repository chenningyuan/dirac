module matrix_stat

  implicit none

  public add_allocation
  public remove_allocation
  public print_nr_allocations

  private

  integer :: nr_allocations = 0

contains

  subroutine add_allocation()
     nr_allocations = nr_allocations + 1
  end subroutine

  subroutine remove_allocation()
     nr_allocations = nr_allocations - 1
  end subroutine

  subroutine print_nr_allocations(text)
     character(*), intent(in) :: text
     write(*, *) 'debug nr allocations ', nr_allocations, text
  end subroutine

end module
