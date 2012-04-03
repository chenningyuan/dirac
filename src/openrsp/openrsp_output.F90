module openrsp_output

  implicit none

  public print_polarizability

  private

  interface print_polarizability
    module procedure print_alpha
    module procedure print_beta
    module procedure print_gamma
    module procedure print_delta
  end interface

  character(1) :: fc(3) = (/'x', 'y', 'z'/)
  real(8)      :: small_number = 1.0d-10

contains

  subroutine print_alpha(e, w, i)

!   ----------------------------------------------------------------------------
    complex(8), intent(in) :: e(3, 3)
    complex(8), intent(in) :: w(1)
    integer,    intent(in) :: i
!   ----------------------------------------------------------------------------
    integer                :: h
!   ----------------------------------------------------------------------------

    do h = 1, 3
      if (dabs(real(e(h, i))) > small_number) then
        print '("@ << "a1"; "a1" >>          (w ="2f8.4")" f50.10)', &
               fc(h),                                                &
               fc(i),                                                &
               real(-sum(w)),                                        &
               real(w(1)),                                           &
              -real(e(h, i))
      end if
    end do

  end subroutine

  subroutine print_beta(e, w, i, j)

!   ----------------------------------------------------------------------------
    complex(8), intent(in) :: e(3, 3, 3)
    complex(8), intent(in) :: w(2)
    integer,    intent(in) :: i, j
!   ----------------------------------------------------------------------------
    integer                :: h
!   ----------------------------------------------------------------------------

    do h = 1, 3
      if (dabs(real(e(h, i, j))) > small_number) then
        print '("@ << "a1"; "a1", "a1" >>       (w ="3f8.4")" f42.10)', &
               fc(h),                                                   &
               fc(i),                                                   &
               fc(j),                                                   &
               real(-sum(w)),                                           &
               real(w(1)),                                              &
               real(w(2)),                                              &
              -real(e(h, i, j))
      end if
    end do

  end subroutine

  subroutine print_gamma(e, w, i, j, k)

!   ----------------------------------------------------------------------------
    complex(8), intent(in) :: e(3, 3, 3, 3)
    complex(8), intent(in) :: w(3)
    integer,    intent(in) :: i, j, k
!   ----------------------------------------------------------------------------
    integer                :: h
!   ----------------------------------------------------------------------------

    do h = 1, 3
      if (dabs(real(e(h, i, j, k))) > small_number) then
        print '("@ << "a1"; "a1", "a1", "a1" >>    (w ="4f8.4")" f34.10)', &
               fc(h),                                                      &
               fc(i),                                                      &
               fc(j),                                                      &
               fc(k),                                                      &
               real(-sum(w)),                                              &
               real(w(1)),                                                 &
               real(w(2)),                                                 &
               real(w(3)),                                                 &
              -real(e(h, i, j, k))
      end if
    end do

  end subroutine

  subroutine print_delta(e, w, i, j, k, l)

!   ----------------------------------------------------------------------------
    complex(8), intent(in) :: e(3, 3, 3, 3, 3)
    complex(8), intent(in) :: w(4)
    integer,    intent(in) :: i, j, k, l
!   ----------------------------------------------------------------------------
    integer                :: h
!   ----------------------------------------------------------------------------

    do h = 1, 3
      if (dabs(real(e(h, i, j, k, l))) > small_number) then
        print '("@ << "a1"; "a1", "a1", "a1", "a1" >> (w ="5f8.4")" f26.10)', &
               fc(h),                                                         &
               fc(i),                                                         &
               fc(j),                                                         &
               fc(k),                                                         &
               fc(l),                                                         &
               real(-sum(w)),                                                 &
               real(w(1)),                                                    &
               real(w(2)),                                                    &
               real(w(3)),                                                    &
               real(w(4)),                                                    &
              -real(e(h, i, j, k, l))
      end if
    end do

  end subroutine

end module
