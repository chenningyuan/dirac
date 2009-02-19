module aoresponse_test_m

    use matrix_defop
    use matrix_genop
    use aoresponse_xc_m
    use aoresponse_visual_m
    use dirac_interface_m
    use prop_integ
    use rsp_equations

    implicit none

    character(1)  :: fc(3) = (/'x', 'y', 'z'/)
    real(kind_dp) :: small = 1.0d-10

contains

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

  subroutine nr_electrons(S, D)

    type(matrix), intent(inout) :: S, D

    S = 1.0d0*S
    D = 1.0d0*D

    write(*, *) 'electron count from S**D:', dot(S, D)

  end subroutine

  subroutine energy(D, F)

    type(matrix), intent(inout) :: D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: G(1)
    real(kind_dp)               :: xc_energy
    character(4)                :: nop(0)

    D = 1.0d0*D
    F = 1.0d0*F

    G(1) = 1.0d-100*D

    call prop_fock(nop, (/D, D/), G)
    write(*, *) 'electronic energy from F**D - 0.5d0*G(1)**D:', dot(F, D) - 0.5d0*dot(G(1), D)

    G(1) = 0

    xc_energy = 0.0d0
!   call xc_contribution(0, (/D/), (/F/), xc_energy)
    write(*, *) 'xc energy:', xc_energy

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

    call mat_print(D, label = 'test_matrix_defop: D')

    call mat_print(F, label = 'test_matrix_defop: F')

    A = F
    call mat_print(A, label = 'test_matrix_defop: A = F')

    A = F + F
    call mat_print(A, label = 'test_matrix_defop: A = F + F')

    A = F - F
    call mat_print(A, label = 'test_matrix_defop: A = F - F')

    A = F*F
    call mat_print(A, label = 'test_matrix_defop: A = F*F')

    A = D*F
    call mat_print(A, label = 'test_matrix_defop: A = D*F')

    A = 2.0d0*F
    call mat_print(A, label = 'test_matrix_defop: A = 2.0d0*F')

    A = F/2.0d0
    call mat_print(A, label = 'test_matrix_defop: A = F/2.0d0')

    A = C_i
    call mat_print(A, label = 'test_matrix_defop: A = C_i')

    A = dag(C_i)
    call mat_print(A, label = 'test_matrix_defop: A = dag(C_i)')

    A = C_i*dag(C_i)
    call mat_print(A, label = 'test_matrix_defop: A = C_i*dag(C_i)')

    A = C_i*dag(C_i)/2.0d0
    call mat_print(A, label = 'test_matrix_defop: A = C_i*dag(C_i)/2.0d0')

    A = F*D*S - S*D*F
    call mat_print(A, label = 'test_matrix_defop: A = F*D*S - S*D*F')

!   write(*, *) 'test_matrix_defop: tr(D)', tr(D)

    write(*, *) 'test_matrix_defop: dot(S, D)', dot(S, D)
    write(*, *) 'test_matrix_defop: tr(S, D)', tr(S, D)

    write(*, *) 'test_matrix_defop: dot(D, F)', dot(D, F)
    write(*, *) 'test_matrix_defop: tr(D, F)', tr(D, F)

    write(*, *) 'test_matrix_defop: norm(F)', norm(F)

    write(*, *) 'test_matrix_defop: sqrt(dot(F, F))', sqrt(dot(F, F))

  end subroutine

  subroutine mu(S, D, F)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    real(kind_dp)               :: e_i(3)
    integer                     :: i, j

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    e_i = 0.0d0
    call prop_energy((/'EL'/), (/D/), e_i)

    do i = 1, 3
      print '("@ -< "a1" > = " f26.10)', fc(i), -e_i(i)
    end do

  end subroutine

  subroutine alpha(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
    real(kind_dp)               :: w(2)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3)
    type(matrix)                :: F_i(3)
    real(kind_dp)               :: e_ij(3, 3), avg
    integer                     :: i, j

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(2:2))

    e_ij = 0.0d0
    call prop_energy((/'EL'/), (/D, D_i/), e_ij)

    if (aoresponse_visual .or. aoresponse_integrate) then
      call aoresponse_visual_driver((/D, D_i/))
    end if

    print '("@ alpha (w = "2f7.3")")', w(1), w(2)

    do i = 1, 3

      D_i(i) = 0
      F_i(i) = 0

      do j = 1, 3
        if (abs(e_ij(i, j)) > small) then
          print '("@" a10 a2 f26.10)', fc(i), fc(j), -e_ij(i, j)
        end if
      end do
    end do

    avg = 0.0d0
    do i = 1, 3
      avg = avg - e_ij(i, i)/3.0d0
    end do

    print '("@ averaged alpha (w = "2f7.3"):" f26.10)', w(1), &
                                                        w(2), &
                                                        avg

  end subroutine

  subroutine beta_debug_xc(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
    real(kind_dp)               :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), D_j(3), D_ij(3, 3), D_k(3), F_k(3)
    type(matrix)                :: F_i(3), F_j(3), F_ij(3, 3)
    type(matrix)                :: X_i(3), X_j(3), X_ij(3, 3), A, L_k(3), Y_k(3)
    real(kind_dp)               :: e_ijk(3, 3, 3)
    integer                     :: i, j, k

    if (abs(sum(w)) > small) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_k, F_k, w=w(2:2))

    e_ijk = 0.0d0

    

    do k = 1, 3
      X_i(k) = 1.0d-100*F
      X_j(k) = 1.0d-100*F
    end do
!   call xc_contribution(e_ijk, 3, D, (/D_i, D_j/), (/X_i, X_j/), (/F_ij/))
!   Lagrange contribution -Y_k (FDS - SDF) e_j
!   do k = 1, 3
!     Y_k(k) = D_k(k)*S*D - D*S*D_k(k)
!   end do

!   do i = 1, 3
!     do j = 1, 3

!       A = (X_j(j)*D_i(i) + X_i(i)*D_j(j))*S - S*(D_i(i)*X_j(j) + D_j(j)*X_i(i))

!       do k = 1, 3
!         e_ijk(i, j, k) = e_ijk(i, j, k) - 0.5d0*tr(Y_k(k), A)
!       end do

!       A = 0

!     end do
!   end do
!   Lagrange contribution -L_k (DSD - D) e_j
    do k = 1, 3
      Y_k(k) = X_i(k)*D*S + S*D*X_i(k) - X_i(k)
    end do

    do i = 1, 3
      do j = 1, 3

        A = D_i(i)*S*D_j(j) + D_j(j)*S*D_i(i)

        do k = 1, 3
          e_ijk(i, j, k) = e_ijk(i, j, k) - tr(Y_k(k), A)
        end do

        A = 0

      end do
    end do

    do k = 1, 3
      F_i(k) = F_i(k) + X_i(k)
      F_j(k) = F_j(k) + X_j(k)
    end do



    call pert_dens((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))

    call prop_energy((/'EL'/), (/D, D_i, D_j, D_ij/), e_ijk, b=(/3, 3/))

    print '("@ beta (w = "3f7.3")")', w(1), w(2), w(3)

    do i = 1, 3

      D_i(i) = 0
      D_j(i) = 0
      F_i(i) = 0
      F_j(i) = 0

      do j = 1, 3

        D_ij(i, j) = 0
        F_ij(i, j) = 0

        do k = 1, 3
          if (abs(e_ijk(i, j, k)) > small) then
            print '("@" a10 a2 a2 f26.10)', fc(i), fc(j), fc(k), -e_ijk(i, j, k)
          end if
        end do
      end do
    end do

  end subroutine

  subroutine beta(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
    real(kind_dp)               :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), D_j(3), D_ij(3, 3)
    type(matrix)                :: F_i(3), F_j(3), F_ij(3, 3)
    real(kind_dp)               :: e_ijk(3, 3, 3)
    integer                     :: i, j, k

    if (abs(sum(w)) > small) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))

    call pert_dens((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))

    if (aoresponse_visual .or. aoresponse_integrate) then
      call aoresponse_visual_driver((/D, D_i, D_ij/))
    end if

    e_ijk = 0.0d0
    call prop_energy((/'EL'/), (/D, D_i, D_j, D_ij/), e_ijk, b=(/3, 3/))

    print '("@ beta (w = "3f7.3")")', w(1), w(2), w(3)

    do i = 1, 3

      D_i(i) = 0
      D_j(i) = 0
      F_i(i) = 0
      F_j(i) = 0

      do j = 1, 3

        D_ij(i, j) = 0
        F_ij(i, j) = 0

        do k = 1, 3
          if (abs(e_ijk(i, j, k)) > small) then
            print '("@" a10 a2 a2 f26.10)', fc(i), fc(j), fc(k), -e_ijk(i, j, k)
          end if
        end do
      end do
    end do

  end subroutine

  subroutine beta_2np1(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
    real(kind_dp)               :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: A, Y_k(3)
    type(matrix)                :: D_i(3), D_j(3), D_k(3)
    type(matrix)                :: F_i(3), F_j(3), F_k(3)
    real(kind_dp)               :: e_ijk(3, 3, 3)
    integer                     :: i, j, k
    character(4)                :: nop(0)

    if (abs(sum(w)) > small) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_k, F_k, w=w(3:3))

#ifdef debug_me
!   radovan: please don't remove this
!            this "MO way" of doing the XC contribution
!            is useful for debugging

!   call get_C(C)
!   call get_i_or_s_C(Ci, 'i')
!   call get_i_or_s_C(Cs, 's')

!   D_a = 1.0d0*D_i(3)
!   D_b = 1.0d0*D_j(3)
!   D_c = 1.0d0*D_k(3)

!   W_b = 0.5d0*(dag(Ci)*S*D_b*S*Cs - dag(Cs)*S*D_b*S*Ci)
!   W_c = 0.5d0*(dag(Ci)*S*D_c*S*Cs - dag(Cs)*S*D_c*S*Ci)

!   D_bc = 0.5d0*(Ci*W_b*W_c*dag(Ci) + Ci*W_c*W_b*dag(Ci) - Cs*W_c*W_b*dag(Cs) - Cs*W_b*W_c*dag(Cs))

!   XC_b  = 1.0d-100*D_b
!   XC_c  = 1.0d-100*D_c
!   XC_bc = 1.0d-100*D_bc

!   call xc_contribution(3, D, (/D_b, D_c, D_bc/), (/XC_b, XC_c, XC_bc/))

!   XC_b_mo = dag(C)*XC_b*C
!   XC_c_mo = dag(C)*XC_c*C

!   XC_b_mo_1c = W_c*XC_b_mo - XC_b_mo*W_c
!   XC_c_mo_1b = W_b*XC_c_mo - XC_c_mo*W_b

!   XC_b_1c = 0.5d0*(S*C*XC_b_mo_1c*dag(C)*S)
!   XC_c_1b = 0.5d0*(S*C*XC_c_mo_1b*dag(C)*S)

!   XC = XC_b_1c + XC_c_1b + 2.0d0*XC_bc

!   write(*,*) 'XC contribution b     = ', -(D_a**XC_b_1c)
!   write(*,*) 'XC contribution c     = ', -(D_a**XC_c_1b)
!   write(*,*) 'XC contribution bc    = ', -2.0d0*(D_a**XC_bc)
!   write(*,*) 'XC contribution total = ', -(D_a**XC)
#endif

    e_ijk = 0.0d0

!   if (aoresponse_test_xc) then
!     call xc_contribution(e_ijk, 3, D, (/D_i, D_j, D_k/), (/F_i, F_j, F_k/))
!   end if

    do k = 1, 3
      Y_k(k) = D_k(k)*S*D - D*S*D_k(k)
    end do

    do i = 1, 3
      do j = 1, 3
        A = (F_j(j)*D_i(i) + F_i(i)*D_j(j))*S - S*(D_i(i)*F_j(j) + D_j(j)*F_i(i))
        do k = 1, 3
          e_ijk(i, j, k) = e_ijk(i, j, k) - tr(Y_k(k), A)
        end do
      end do
    end do

    do k = 1, 3
      Y_k(k) = F_k(k)*D*S + S*D*F_k(k) - F_k(k)
    end do

    do i = 1, 3
      do j = 1, 3
        A = D_i(i)*S*D_j(j) + D_j(j)*S*D_i(i)
        do k = 1, 3
          e_ijk(i, j, k) = e_ijk(i, j, k) - tr(Y_k(k), A)
        end do
      end do
    end do

    print '("@ beta (w = "3f7.3")")', w(1), w(2), w(3)

    A = 0

    do i = 1, 3

      D_i(i) = 0
      D_j(i) = 0
      D_k(i) = 0

      F_i(i) = 0
      F_j(i) = 0
      F_k(i) = 0

      Y_k(i) = 0

      do j = 1, 3
        do k = 1, 3
          if (abs(e_ijk(i, j, k)) > small) then
            print '("@" a10 a2 a2 f26.10)', fc(i), fc(j), fc(k), -e_ijk(i, j, k)
          end if
        end do
      end do
    end do

  end subroutine

  subroutine beta_element(scale_by, add_to, ui, uj, uk, S, D, F, w)

    real(kind_dp), intent(in)    :: scale_by
    real(kind_dp), intent(inout) :: add_to
    integer,       intent(in)    :: ui, uj, uk
    type(matrix),  intent(inout) :: S, D, F
    real(kind_dp), intent(in)    :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                 :: D_i(3), D_j(3), D_ij(3, 3)
    type(matrix)                 :: F_i(3), F_j(3), F_ij(3, 3)
    real(kind_dp)                :: e_ijk(3, 3, 3)
    integer                      :: i, j, k

    if (abs(sum(w)) > small) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do i = 1, 3

      D_i(i) = 1.0d-100*D; D_j(i) = 1.0d-100*D
      F_i(i) = 1.0d-100*D; F_j(i) = 1.0d-100*D

      do j = 1, 3

        D_ij(i, j) = 1.0d-100*D
        F_ij(i, j) = 1.0d-100*D

      end do
    end do

    call pert_dens((/'EL'/),   &
                        S,        &
                      (/D/),      &
                      (/F/),      &
                        D_i(ui),  &
                        F_i(ui),  &
                        w=w(1:1), &
                        c=(/ui/))

    call pert_dens((/'EL'/),   &
                        S,        &
                      (/D/),      &
                      (/F/),      &
                        D_j(uj),  &
                        F_j(uj),  &
                        w=w(2:2), &
                        c=(/uj/))

    call pert_dens((/'EL', 'EL'/),     &
                        S,                &
                      (/D,                &
                        D_i(ui),          &
                        D_j(uj)/),        &
                      (/F,                &
                        F_i(ui),          &
                        F_j(uj)/),        &
                        D_ij(ui, uj),     &
                        F_ij(ui, uj),     &
                        w=(/w(1), w(2)/), &
                        c=(/ui, uj/))

    e_ijk = 0.0d0
    call prop_energy((/'EL'/), &
                    (/D,      &
                      D_i,    &
                      D_j,    &
                      D_ij/), &
                      e_ijk,  &
                      b=(/3, 3/))

    do i = 1, 3

      D_i(i) = 0; D_j(i) = 0
      F_i(i) = 0; F_j(i) = 0

      do j = 1, 3

        D_ij(i, j) = 0
        F_ij(i, j) = 0

        do k = 1, 3

          if (abs(e_ijk(i, j, k)) > small) then
            print '("@ beta (w = "3f8.4")" a10 a2 a2 f26.10)', w(1),  &
                                                               w(2),  &
                                                               w(3),  &
                                                               fc(i), &
                                                               fc(j), &
                                                               fc(k), &
                                                              -e_ijk(i, j, k)
            add_to = add_to - scale_by*e_ijk(i, j, k)
          end if
        end do
      end do
    end do

  end subroutine

  subroutine gamma_element(scale_by, add_to, ui, uj, uk, ul, S, D, F, w)

    real(kind_dp), intent(in)    :: scale_by
    real(kind_dp), intent(inout) :: add_to
    integer,       intent(in)    :: ui, uj, uk, ul
    type(matrix),  intent(inout) :: S, D, F
    real(kind_dp), intent(in)    :: w(4)
!   ----------------------------------------------------------------------------
    type(matrix)                 :: D_i(3), D_j(3), D_k(3), D_ij(3, 3), D_ik(3, 3), D_jk(3, 3), D_ijk(3, 3, 3)
    type(matrix)                 :: F_i(3), F_j(3), F_k(3), F_ij(3, 3), F_ik(3, 3), F_jk(3, 3), F_ijk(3, 3, 3)
    real(kind_dp)                :: e_ijkl(3, 3, 3, 3)
    integer                      :: i, j, k, l

    if (abs(sum(w)) > small) call quit('gamma: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do i = 1, 3

      D_i(i) = 1.0d-100*D; D_j(i) = 1.0d-100*D; D_k(i) = 1.0d-100*D
      F_i(i) = 1.0d-100*D; F_j(i) = 1.0d-100*D; F_k(i) = 1.0d-100*D

      do j = 1, 3

        D_ij(i, j) = 1.0d-100*D; D_ik(i, j) = 1.0d-100*D; D_jk(i, j) = 1.0d-100*D
        F_ij(i, j) = 1.0d-100*D; F_ik(i, j) = 1.0d-100*D; F_jk(i, j) = 1.0d-100*D

        do k = 1, 3

          D_ijk(i, j, k) = 1.0d-100*D
          F_ijk(i, j, k) = 1.0d-100*D

        end do
      end do
    end do

    call pert_dens((/'EL'/),   &
                        S,        &
                      (/D/),      &
                      (/F/),      &
                        D_i(ui),  &
                        F_i(ui),  &
                        w=w(1:1), &
                        c=(/ui/))

    call pert_dens((/'EL'/),   &
                        S,        &
                      (/D/),      &
                      (/F/),      &
                        D_j(uj),  &
                        F_j(uj),  &
                        w=w(2:2), &
                        c=(/uj/))

    call pert_dens((/'EL'/),   &
                        S,        &
                      (/D/),      &
                      (/F/),      &
                        D_k(uk),  &
                        F_k(uk),  &
                        w=w(3:3), &
                        c=(/uk/))

    call pert_dens((/'EL', 'EL'/),     &
                        S,                &
                      (/D,                &
                        D_i(ui),          &
                        D_j(uj)/),        &
                      (/F,                &
                        F_i(ui),          &
                        F_j(uj)/),        &
                        D_ij(ui, uj),     &
                        F_ij(ui, uj),     &
                        w=(/w(1), w(2)/), &
                        c=(/ui, uj/))

    call pert_dens((/'EL', 'EL'/),     &
                        S,                &
                      (/D,                &
                        D_i(ui),          &
                        D_k(uk)/),        &
                      (/F,                &
                        F_i(ui),          &
                        F_k(uk)/),        &
                        D_ik(ui, uk),     &
                        F_ik(ui, uk),     &
                        w=(/w(1), w(3)/), &
                        c=(/ui, uk/))

    call pert_dens((/'EL', 'EL'/),     &
                        S,                &
                      (/D,                &
                        D_j(uj),          &
                        D_k(uk)/),        &
                      (/F,                &
                        F_j(uj),          &
                        F_k(uk)/),        &
                        D_jk(uj, uk),     &
                        F_jk(uj, uk),     &
                        w=(/w(2), w(3)/), &
                        c=(/uj, uk/))

    call pert_dens((/'EL', 'EL', 'EL'/),     &
                        S,                      &
                      (/D,                      &
                        D_i(ui),                &
                        D_j(uj),                &
                        D_k(uk),                &
                        D_ij(ui, uj),           &
                        D_ik(ui, uk),           &
                        D_jk(uj, uk)/),         &
                      (/F,                      &
                        F_i(ui),                &
                        F_j(uj),                &
                        F_k(uk),                &
                        F_ij(ui, uj),           &
                        F_ik(ui, uk),           &
                        F_jk(uj, uk)/),         &
                        D_ijk(ui, uj, uk),      &
                        F_ijk(ui, uj, uk),      &
                        w=(/w(1), w(2), w(3)/), &
                        c=(/ui, uj, uk/))

    e_ijkl = 0.0d0
    call prop_energy((/'EL'/),  &
                    (/D,       &
                      D_i,     &
                      D_j,     &
                      D_k,     &
                      D_ij,    &
                      D_ik,    &
                      D_jk,    &
                      D_ijk/), &
                      e_ijkl,  &
                      b=(/3, 3, 3/))

    do i = 1, 3

      D_i(i) = 0; D_j(i) = 0; D_k(i) = 0
      F_i(i) = 0; F_j(i) = 0; F_k(i) = 0

      do j = 1, 3

        D_ij(i, j) = 0; D_ik(i, j) = 0; D_jk(i, j) = 0
        F_ij(i, j) = 0; F_ik(i, j) = 0; F_jk(i, j) = 0

        do k = 1, 3

          D_ijk(i, j, k) = 0
          F_ijk(i, j, k) = 0

          do l = 1, 3
            if (abs(e_ijkl(i, j, k, l)) > small) then
              print '("@ gamma (w = "4f9.4")" a10 a2 a2 a2 f26.10)', w(1),  &
                                                                     w(2),  &
                                                                     w(3),  &
                                                                     w(4),  &
                                                                     fc(i), &
                                                                     fc(j), &
                                                                     fc(k), &
                                                                     fc(l), &
                                                                    -e_ijkl(i, j, k, l)
              add_to = add_to - scale_by*e_ijkl(i, j, k, l)
            end if
          end do
        end do
      end do
    end do

  end subroutine

  subroutine gamma(S, D, F, w)

    type(matrix),  intent(inout) :: S, D, F
    real(kind_dp), intent(in)    :: w(4)
!   ----------------------------------------------------------------------------
    type(matrix)                 :: D_i(3), D_j(3), D_k(3), D_ij(3, 3), D_ik(3, 3), D_jk(3, 3), D_ijk(3, 3, 3)
    type(matrix)                 :: F_i(3), F_j(3), F_k(3), F_ij(3, 3), F_ik(3, 3), F_jk(3, 3), F_ijk(3, 3, 3)
    real(kind_dp)                :: e_ijkl(3, 3, 3, 3)
    integer                      :: i, j, k, l

    if (abs(sum(w)) > small) call quit('gamma: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_k, F_k, w=w(3:3))

    call pert_dens((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))
    call pert_dens((/'EL', 'EL'/), S, (/D, D_i, D_k/), (/F, F_i, F_k/), D_ik, F_ik, w=(/w(1), w(3)/))
    call pert_dens((/'EL', 'EL'/), S, (/D, D_j, D_k/), (/F, F_j, F_k/), D_jk, F_jk, w=w(2:3))

    call pert_dens((/'EL', 'EL', 'EL'/), S, (/D, D_i, D_j, D_k, D_ij, D_ik, D_jk/), &
                                               (/F, F_i, F_j, F_k, F_ij, F_ik, F_jk/), &
                                               D_ijk, F_ijk, w=w(1:3))

    e_ijkl = 0.0d0
    call prop_energy((/'EL'/), (/D, D_i, D_j, D_k, D_ij, D_ik, D_jk, D_ijk/), e_ijkl, b=(/3, 3, 3/))

    print '("@ gamma (w = "4f9.4")")', w(1), w(2), w(3), w(4)

    do i = 1, 3

      D_i(i) = 0
      D_j(i) = 0
      D_k(i) = 0
      F_i(i) = 0
      F_j(i) = 0
      F_k(i) = 0

      do j = 1, 3

        D_ij(i, j) = 0
        D_ik(i, j) = 0
        D_jk(i, j) = 0
        F_ij(i, j) = 0
        F_ik(i, j) = 0
        F_jk(i, j) = 0

        do k = 1, 3

          D_ijk(i, j, k) = 0
          F_ijk(i, j, k) = 0

          do l = 1, 3
            if (abs(e_ijkl(i, j, k, l)) > small) then
              print '("@" a10 a2 a2 a2 f26.10)', fc(i), fc(j), fc(k), fc(l), -e_ijkl(i, j, k, l)
            end if
          end do
        end do
      end do
    end do

  end subroutine







subroutine beta_debug(S,D,F,w)
    type(Matrix):: S,D,Df(3),De(3),Dfe(3*3), T, M
    type(Matrix)::   F,Ff(3),Fe(3),Ffe(3*3), G(1), H(1)
    real*8:: w(3),Efff(3*3*3); integer:: i,j, ij
    character*(4) :: NOP(0)
    S=S; D=D; F=F !track them
    if (abs(sum(w))>1d-15) call QUIT('PropTest/HyPolari: sum(w) should be zero!')
    call pert_dens((/'EL'/),S,(/D/),(/F/),Df,Ff,w=w(1:1))
    do i = 1, 3

       H(1) = 1.0d-50*D
       call prop_fock((/'EL'/),(/D/),H,c=(/i/))
       G(1) = 1.0d-50*D
       call prop_fock(NOP,(/D,Df(i)/),G)

       T = H(1) + G(1) - Ff(i)
!      write(*, *) 'debug_norm 1:', T**T

       T = Df(i)*S*D + D*S*Df(i) - Df(i)
!      write(*, *) 'debug_norm 2:', T**T

!      M = H(1)*D*S - S*D*H(1)
!      write(*, *) 'debug matrix 1:'
!      call mat_print(M)

!      M = G(1)*D*S - S*D*G(1)
!      write(*, *) 'debug matrix 2:'
!      call mat_print(M)

!      M = Ff(i)*D*S - S*D*Ff(i)
!      write(*, *) 'debug matrix 3:'
!      call mat_print(M)

!      M = F*Df(i)*S - S*Df(i)*F
!      write(*, *) 'debug matrix 4:'
!      call mat_print(M)

       T = Ff(i)*D*S - S*D*Ff(i) + F*Df(i)*S - S*Df(i)*F
!      write(*, *) 'debug matrix 5:'
!      call mat_print(T)
!      write(*, *) 'debug_norm 4:', T**T

    end do
!   call quit('over and out')
!   call pert_dens((/'EL'/),S,(/D/),(/F/),Df,Ff,w=w(1:1))
    call pert_dens((/'EL'/),S,(/D/),(/F/),De,Fe,w=w(2:2))
    call pert_dens((/'EL','EL'/),S,(/D,Df,De/),(/F,Ff,Fe/),Dfe,Ffe,w=w(1:2))
    do i = 1, 3
       do j = 1, 3
          ij = i + 3*(j-1)
          T = Dfe(ij)*S*D + De(j)*S*Df(i) + Df(i)*S*De(j) + D*S*Dfe(ij) - Dfe(ij)
!      write(*, *) 'debug_norm 5:', T**T
      end do
    end do

    Efff=0; call prop_energy((/'EL'/),(/D,Df,De,Dfe/),Efff,b=(/3,3/))
    print '("First hyperpolarizability -Efff w="3f6.3/9(3f26.13/))',w,-Efff
    do i=1,3; Df(i)=0; De(i)=0; Ff(i)=0; Fe(i)=0; enddo
    do i=1,3*3; Dfe(i)=0; Ffe(i)=0; enddo
    print*
end subroutine


end module
