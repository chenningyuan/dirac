module aor_property_definitions_zzz

    use matrix_defop
    use matrix_genop
    use dirac_interface
    use rsp_equations

    implicit none

contains

  subroutine aor_zzz(order, S, D, F)

!   ----------------------------------------------------------------------------
    integer,      intent(in)    :: order
    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    integer,      parameter     :: max_order = 7
    type(matrix)                :: P, Dp(max_order), Fp(max_order)
    integer                     :: iorder
    real(8)                     :: g(max_order)
    integer                     :: i(max_order)
    character(2)                :: l(max_order)
!   ----------------------------------------------------------------------------

    if (order > max_order) call quit('aor_zzz: order > max_order')

    g = 0.0d0
    i = 3
    l = 'EL'

    call get_P('ZDIPLEN', P)

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do iorder = 1, max_order
      Dp(iorder) = 1.0d-100*D
      Fp(iorder) = 1.0d-100*F
    end do

    iorder = 1
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D/),           &
                   (/F/),           &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if

    iorder = 2
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D,             &
                     Dp(1), &
                     Dp(1)/), &
                   (/F,             &
                     Fp(1), &
                     Fp(1)/), &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if
          
    iorder = 3
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D,             &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2)/), &
                   (/F,             &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2)/), &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if
          
    iorder = 4
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D,             &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3)/), &
                   (/F,             &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3)/), &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if
          
    iorder = 5
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D,             &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4)/), &
                   (/F,             &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4)/), &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if
          
    iorder = 6
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D,             &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5)/), &
                   (/F,             &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5)/), &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if
          
    iorder = 7
    if (order >= iorder) then
      call pert_dens(l(1:iorder),   &
                     S,             &
                   (/D,             &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(1), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(2), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(3), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(4), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(5), &
                     Dp(6), &
                     Dp(6), &
                     Dp(6), &
                     Dp(6), &
                     Dp(6), &
                     Dp(6), &
                     Dp(6)/), &
                   (/F,             &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(1), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(2), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(3), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(4), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(5), &
                     Fp(6), &
                     Fp(6), &
                     Fp(6), &
                     Fp(6), &
                     Fp(6), &
                     Fp(6), &
                     Fp(6)/), &
                     Dp(iorder),    &
                     Fp(iorder),    &
                     w=g(1:iorder), &
                     c=i(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder, dot(P, Dp(iorder))
    end if
    
    do iorder = 1, max_order
      Dp(iorder) = 0
      Fp(iorder) = 0
    end do

    P = 0

  end subroutine

end module
