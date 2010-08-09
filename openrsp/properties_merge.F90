module openrsp_property_definitions_zzz

    use matrix_defop
    use matrix_genop
    use openrsp_interface_1el
    use dirac_interface
    use rsp_equations
    use decompmod

    implicit none

    public openrsp_zzz

    private

contains

  subroutine openrsp_zzz(order, S, D, F)

!   ----------------------------------------------------------------------------
    integer,      intent(in)    :: order
    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    integer,      parameter     :: max_order = 7
    type(matrix)                :: P, Dp(max_order), Fp(max_order)
    integer                     :: iorder
    complex(8)                  :: g(max_order)
    integer                     :: i(max_order)
    integer                     :: dimp(max_order)
    character(2)                :: l(max_order)
    type(DecompItem) :: decomp
!   ----------------------------------------------------------------------------

    if (order > max_order) call quit('openrsp_zzz: order > max_order')

    g = 0.0d0
    i = 3
    dimp = 1
    l = 'EL'

    call read_1el_integrals('ZDIPLEN', P)

    do iorder = 1, max_order
      Dp(iorder) = 1.0d-100*D
      Fp(iorder) = 1.0d-100*F
    end do

    iorder = 1
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D/),              &
                   (/F/),              &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if

    iorder = 2
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D,                &
                     Dp(1),            &
                     Dp(1)/),          &
                   (/F,                &
                     Fp(1),            &
                     Fp(1)/),          &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if
          
    iorder = 3
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D,                &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2)/),          &
                   (/F,                &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2)/),          &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if
          
    iorder = 4
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D,                &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3)/),          &
                   (/F,                &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3)/),          &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if
          
    iorder = 5
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D,                &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3)/),          &
                   (/F,                &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3)/),          &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if
          
    iorder = 6
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D,                &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4)/),          &
                   (/F,                &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4)/),          &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if
          
    iorder = 7
    if (order >= iorder) then
      call pert_dens(decomp,           &
                     S,                &
                     l(1:iorder),      &
                     dimp(1:iorder),   &
                   (/D,                &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(1),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(2),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(3),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(4),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(5),            &
                     Dp(6),            &
                     Dp(6),            &
                     Dp(6),            &
                     Dp(6),            &
                     Dp(6),            &
                     Dp(6),            &
                     Dp(6)/),          &
                   (/F,                &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(1),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(2),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(3),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(4),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(5),            &
                     Fp(6),            &
                     Fp(6),            &
                     Fp(6),            &
                     Fp(6),            &
                     Fp(6),            &
                     Fp(6),            &
                     Fp(6)/),          &
                     Dp(iorder),       &
                     Fp(iorder),       &
                     comp=i(1:iorder), &
                     freq=g(1:iorder))

      print '("@ static E_{zzz...} order "i2": "f25.10)', iorder + 1, real(dot(P, Dp(iorder)))
    end if
    
    do iorder = 1, max_order
      Dp(iorder) = 0
      Fp(iorder) = 0
    end do

    P = 0

  end subroutine

#ifdef MERGE_WITH_LINSCA
  subroutine test_grcont(n, nr_atoms, test_london, test_geo)

    integer             :: n, nr_atoms
    logical             :: test_london, test_geo
!   ----------------------------------------------------------------------------
    real*8, allocatable :: work(:)
    integer             :: lwork

    integer             :: i, j, k, l, m
#ifdef PRG_DIRAC
    real*8              :: D(n, n, 8)
#else                   
    real*8              :: D(n, n, 2)
#endif /* PRG_DIRAC */
    real*8              :: F(n, n, 3), g(3, nr_atoms)

#ifdef PRG_DIRAC
#include "dcbbas.h"
#endif /* PRG_DIRAC */

    lwork = 10000*n + 50*n*n
    allocate(work(lwork))

    if (test_london) then

      print '(a)', 'London 2e integrals'
      print '(a, 12x, a2, 18x, a2, 18x, a2)', '  i  j  k  l', 'Bx', 'By', 'Bz'

!     do l = 1, n
!       do k = 1, n
      do l = 1, 1
        do k = 1, 1

          D(:, :, 1) = 0.0d0
          D(k, l, 1) = 1.0d0
          F          = 0.0d0

          call grcont(work, lwork, F, n*n*3, &
                      .false., .true., 1, 0, .false., .true., D(:, :, 1), 1)

          do j = 1, n
            do i = 1, n
              write(0, '(4i3, 3f20.12)') i, j, k, l, F(i, j, :)
            end do
          end do
        end do
      end do
    end if

    if (test_geo) then

      print '(a)', 'Geometry 2e integrals'
      print '(a, 12x, a2, 18x, a2, 18x, a2)', '  i  j  k  l  a', 'Rx', 'Ry', 'Rz'

      D = 0.0d0

!     do l = 1, n
!       do k = 1, n
      do l = 1, 1
        do k = 1, 1
#ifdef PRG_DIRAC
          do j = 1, ntbas(1)
            do i = 1, ntbas(1)
#else
          do j = 1, n
            do i = 1, n
#endif /* PRG_DIRAC */

              D(k, l, 1) = 1.0d0
#ifdef PRG_DIRAC
              D(i, j, 5) = 1.0d0
#else
              D(i, j, 2) = 1.0d0
#endif /* PRG_DIRAC */
              g          = 0.0d0

              call grcont(work, lwork, g, 3*nr_atoms, &
                          .true., .false., 1, 0, .true., .false., D, 2)

              D(k, l, 1) = 0.0d0
#ifdef PRG_DIRAC
              D(i, j, 5) = 0.0d0
#else
              D(i, j, 2) = 0.0d0
#endif /* PRG_DIRAC */

              do m = 1, nr_atoms
                write(0, '(5i3, 3f20.12)') i, j, k, l, m, g(:, m)
              end do
            end do
          end do
        end do
      end do
    end if

    deallocate(work)

  end subroutine

  subroutine test_matrix_defop(S, D, F, C_i)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: S, D, F, C_i
!   ----------------------------------------------------------------------------
    type(matrix)                :: A
!   ----------------------------------------------------------------------------

    S   = 1.0d0*S
    D   = 1.0d0*D
    F   = 1.0d0*F
    C_i = 1.0d0*C_i

    call print_mat(D, label = 'test_matrix_defop: D')

    call print_mat(F, label = 'test_matrix_defop: F')

    A = F
    call print_mat(A, label = 'test_matrix_defop: A = F')

    A = F + F
    call print_mat(A, label = 'test_matrix_defop: A = F + F')

    A = F - F
    call print_mat(A, label = 'test_matrix_defop: A = F - F')

    A = F*F
    call print_mat(A, label = 'test_matrix_defop: A = F*F')

    A = D*F
    call print_mat(A, label = 'test_matrix_defop: A = D*F')

    A = 2.0d0*F
    call print_mat(A, label = 'test_matrix_defop: A = 2.0d0*F')

    A = C_i
    call print_mat(A, label = 'test_matrix_defop: A = C_i')

    A = F - dag(F)
    call print_mat(A, label = 'test_matrix_defop: A = F - dag(F)')

    A = C_i*dag(C_i)
    call print_mat(A, label = 'test_matrix_defop: A = C_i*dag(C_i)')

    A = 0.5d0*C_i*dag(C_i)
    call print_mat(A, label = 'test_matrix_defop: A = 0.5d0*C_i*dag(C_i)')

    A = F*D*S - S*D*F
    call print_mat(A, label = 'test_matrix_defop: A = F*D*S - S*D*F')

    write(*, *) 'test_matrix_defop: dot(S, D)',       dot(S, D)
    write(*, *) 'test_matrix_defop: tr(S, D)',        tr(S, D)
    write(*, *) 'test_matrix_defop: dot(D, F)',       dot(D, F)
    write(*, *) 'test_matrix_defop: tr(D, F)',        tr(D, F)
    write(*, *) 'test_matrix_defop: norm(F)',         norm(F)
    write(*, *) 'test_matrix_defop: sqrt(dot(F, F))', sqrt(dot(F, F))

  end subroutine
#endif /* ifdef MERGE_WITH_LINSCA */

end module
