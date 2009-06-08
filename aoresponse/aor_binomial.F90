module aor_binomial

  implicit none

  public init_binomial
  public get_binomial

  private

  integer, parameter :: max_binomial_order = 10
  integer            :: binomial(max_binomial_order, max_binomial_order)

contains

  subroutine init_binomial()

!   ----------------------------------------------------------------------------
    integer, external :: taylor_len
    integer           :: n, k
!   ----------------------------------------------------------------------------

    binomial = 0

    do k = 0, max_binomial_order
      do n = 0, max_binomial_order
        if ((k + n) <= max_binomial_order) then
          binomial(k + n, n) = taylor_len(k, n)
        end if
      end do
    end do

  end subroutine

  function get_binomial(k, n)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: k, n
    integer             :: get_binomial
!   ----------------------------------------------------------------------------

    get_binomial = binomial(k, n)

  end function

end module
