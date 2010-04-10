module aor_property_definitions

    use matrix_defop
    use matrix_genop
!   use aoresponse_visual_m
    use aoresponse_output
    use dirac_interface
    use prop_contribs
    use rsp_equations
    use aor_cfg

    implicit none

    character(1)  :: fc(3) = (/'x', 'y', 'z'/)
    real(8) :: small_number = 1.0d-10

contains

   subroutine aor_magnetizability(S, D, F, w)
   !no-Lon and London magnetizability
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: w
      complex(8)                :: Eoo(3,3), Ebb(3,3)
      type(matrix) :: Db(3), Fb(3), DFD, ihTb(3), Sb(3), Sbb(3,3)
      integer      :: i,j

      !no-London
      call pert_dens( (/'MAGO'/), S, (/D/), (/F/), Db, Fb, w=(/w/) )
      Eoo = 0
      call prop_energy( (/'MAGO','MAGO'/), (/D/), Eoo); Eoo=-Eoo !sign change
      print '(a/3(3f18.13/)/)', 'E0oo', Eoo; Eoo=0
      call prop_energy( (/'MAGO'/), (/D,Db/), Eoo)
      print '(a/3(3f18.13/)/)', 'E1oDo', Eoo; Eoo=0
      Db=0; Fb=0

!     !London diamag
!     call pert_dens( (/'MAG'/), S, (/D/), (/F/), Db, Fb, w=(/w/) )
!     Ebb = 0
!     call prop_energy( (/'MAG','MAG'/), (/D/), Ebb); Ebb=-Ebb !sign change
!     print '(a/3(3f18.13/)/)', 'E0bb', Ebb; Ebb=0
!     !here, the term -i/2TbbD is missing, which is zero
!     DFD = D*F*D
!     call prop_overlap( (/'MAG','MAG'/), S, Sbb)
!     do j = 1, 3
!        do i = 1, 3
!           Ebb(i,j) = Ebb(i,j) + tr(Sbb(i,j),DFD) !sign change
!        end do
!     end do; DFD=0; Sbb=0
!     print '(a/3(3f18.13/)/)', '-SbbW', Ebb; Ebb=0

!     !London paramag
!     call prop_energy( (/'MAG'/), (/D,Db/), Ebb, o=(/2/))
!     print '(a/3(3f18.13/)/)', 'E1bDb', Ebb; Ebb=0

!     call prop_overlap( (/'MAG'/), S, ihTb, w=(/-w/2,0*w/) )
!     do j = 1, 3
!        do i = 1, 3
!           Ebb(i,j) = Ebb(i,j) - tr(Sb(i),Db(j))
!        end do
!     end do
!     print '(a/3(3f18.13/)/)', '-i/2TbDb', Ebb; Ebb=0

!     call prop_overlap( (/'MAG'/), S, ihTb, w=(/+w/2,0*w/) )
!     do j = 1, 3
!        DFD = Db(j)*(F+(w/2)*S)*D + D*(Fb(j)-ihTb(j))*D + D*(F-(w/2)*S)*Db(j)
!        do i = 1, 3
!           Ebb(i,j) = Ebb(i,j) - tr(Sb(i),DFD)
!        end do
!     end do; DFD=0; Sb=0; ihTb=0; Db=0; Fb=0
!     print '(a/3(3f18.13/)/)', '-SbWb', Ebb; Ebb=0

!     print '(a,2f6.3/3(3f18.13/)/)', 'London Magnetizability = Ebb', -w, w, Ebb
   end subroutine

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
    real(8)               :: xc_energy
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

  subroutine aor_quadrupole(S, D, F)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    complex(8)                  :: e_ij(6) = 0.0d0
    integer                     :: i, j

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

!   e_i = 0.0d0
    call prop_energy((/'ELGR'/), (/D/), e_ij)

    do i = 1, 6
      write(*, *) 'blub', i, e_ij(i)
!     print '("@ -< "a1" > = " f26.10)', fc(i), -real(e_i(i))
    end do

  end subroutine

  subroutine mu(S, D, F)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    complex(8)               :: e_i(3)
    integer                     :: i, j

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    e_i = 0.0d0
    call prop_energy((/'EL'/), (/D/), e_i)

    do i = 1, 3
      print '("@ -< "a1" > = " f26.10)', fc(i), -real(e_i(i))
    end do

  end subroutine

  subroutine alpha(S, D, F, w)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: S, D, F
    complex(8),      intent(in)    :: w(2)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), F_i(3)
    complex(8)                     :: e_ij(3, 3), avg
    integer                     :: i, j
!   ----------------------------------------------------------------------------

    if (abs(sum(w)) > small_number) call quit('alpha: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do i = 1, 3
      D_i(i) = 1.0d-100*D
      F_i(i) = 1.0d-100*F
    end do

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(2:2))

    e_ij = 0.0d0
    call prop_energy((/'EL'/), (/D, D_i/), e_ij)

!   if (aor_cfg_visual .or. aor_cfg_integrate) then
!     call aoresponse_visual_driver((/D, D_i/))
!   end if

    print '("@ alpha (w = "2f7.3")")', real(w(1)), real(w(2))

    do i = 1, 3

      D_i(i) = 0
      F_i(i) = 0

      do j = 1, 3
        if (abs(e_ij(i, j)) > small_number) then
          print '("@" a10 a2 f26.10)', fc(i), fc(j), -real(e_ij(i, j))
        end if
      end do
    end do

    avg = 0.0d0
    do i = 1, 3
      avg = avg - real(e_ij(i, i))/3.0d0
    end do

    print '("@ averaged alpha (w = "2f7.3"):" f26.10)', real(w(1)), &
                                                        real(w(2)), &
                                                        avg

  end subroutine


  subroutine beta(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
    complex(8)               :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), D_j(3), D_ij(3, 3)
    type(matrix)                :: F_i(3), F_j(3), F_ij(3, 3)
    complex(8)               :: e_ijk(3, 3, 3)
    integer                     :: i, j, k

    if (abs(sum(w)) > small_number) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))

    call pert_dens((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))

!   if (aor_cfg_visual .or. aor_cfg_integrate) then
!     call aoresponse_visual_driver((/D, D_i, D_ij/))
!   end if

    e_ijk = 0.0d0
    call prop_energy((/'EL'/), (/D, D_i, D_j, D_ij/), e_ijk, b=(/3, 3/))

    print '("@ beta (w = "3f7.3")")', real(w(1)), real(w(2)), real(w(3))

    do i = 1, 3

      D_i(i) = 0
      D_j(i) = 0
      F_i(i) = 0
      F_j(i) = 0

      do j = 1, 3

        D_ij(i, j) = 0
        F_ij(i, j) = 0

        do k = 1, 3
          if (abs(e_ijk(i, j, k)) > small_number) then
            print '("@" a10 a2 a2 f26.10)', fc(i), fc(j), fc(k), -real(e_ijk(i, j, k))
          end if
        end do
      end do
    end do

  end subroutine

  subroutine beta_2np1(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
    complex(8)               :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: A, Y_k(3)
    type(matrix)                :: D_i(3), D_j(3), D_k(3)
    type(matrix)                :: F_i(3), F_j(3), F_k(3)
    real(8)               :: e_ijk(3, 3, 3)
    integer                     :: i, j, k
    character(4)                :: nop(0)

    if (abs(sum(w)) > small_number) call quit('beta: frequencies do not sum up to zero')

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

    print '("@ beta (w = "3f7.3")")', real(w(1)), real(w(2)), real(w(3))

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
          if (abs(e_ijk(i, j, k)) > small_number) then
            print '("@" a10 a2 a2 f26.10)', fc(i), fc(j), fc(k), -real(e_ijk(i, j, k))
          end if
        end do
      end do
    end do

  end subroutine

  subroutine beta_element(scale_by, add_to, ui, uj, uk, S, D, F, w)

    real(8), intent(in)    :: scale_by
    real(8), intent(inout) :: add_to
    integer,       intent(in)    :: ui, uj, uk
    type(matrix),  intent(inout) :: S, D, F
    complex(8), intent(in)    :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                 :: D_i(3), D_j(3), D_ij(3, 3)
    type(matrix)                 :: F_i(3), F_j(3), F_ij(3, 3)
    complex(8)                :: e_ijk(3, 3, 3)
    integer                      :: i, j, k

    if (abs(sum(w)) > small_number) call quit('beta: frequencies do not sum up to zero')

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

          if (abs(e_ijk(i, j, k)) > small_number) then
            print '("@ beta (w = "3f8.4")" a10 a2 a2 f26.10)', real(w(1)),  &
                                                               real(w(2)),  &
                                                               real(w(3)),  &
                                                               fc(i), &
                                                               fc(j), &
                                                               fc(k), &
                                                              -real(e_ijk(i, j, k))
            add_to = add_to - scale_by*real(e_ijk(i, j, k))
          end if
        end do
      end do
    end do

  end subroutine

  function gamma_element(S, D, F, w, h, i, j, k)

!   ----------------------------------------------------------------------------
    real(8)                     :: gamma_element
    type(matrix), intent(inout) :: S, D, F
    integer,      intent(in)    :: h, i, j, k
    complex(8),   intent(in)    :: w(3)
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3)
    type(matrix)                :: D_j(3)
    type(matrix)                :: D_k(3)
    type(matrix)                :: D_ij(3, 3)
    type(matrix)                :: D_ik(3, 3)
    type(matrix)                :: D_jk(3, 3)
    type(matrix)                :: D_ijk(3, 3, 3)
    type(matrix)                :: F_i(3)
    type(matrix)                :: F_j(3)
    type(matrix)                :: F_k(3)
    type(matrix)                :: F_ij(3, 3)
    type(matrix)                :: F_ik(3, 3)
    type(matrix)                :: F_jk(3, 3)
    type(matrix)                :: F_ijk(3, 3, 3)
    complex(8)                  :: e(3, 3, 3, 3)
    integer                     :: c1, c2, c3
!   ----------------------------------------------------------------------------

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do c1 = 1, 3
      D_i(c1) = tiny(0.0d0)*D
      D_j(c1) = tiny(0.0d0)*D
      D_k(c1) = tiny(0.0d0)*D
      F_i(c1) = tiny(0.0d0)*D
      F_j(c1) = tiny(0.0d0)*D
      F_k(c1) = tiny(0.0d0)*D
      do c2 = 1, 3
        D_ij(c1, c2) = tiny(0.0d0)*D
        D_ik(c1, c2) = tiny(0.0d0)*D
        D_jk(c1, c2) = tiny(0.0d0)*D
        F_ij(c1, c2) = tiny(0.0d0)*D
        F_ik(c1, c2) = tiny(0.0d0)*D
        F_jk(c1, c2) = tiny(0.0d0)*D
        do c3 = 1, 3
          D_ijk(c1, c2, c3) = tiny(0.0d0)*D
          F_ijk(c1, c2, c3) = tiny(0.0d0)*D
        end do
      end do
    end do

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i(i), F_i(i), w=(/w(1)/), c=(/i/))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j(j), F_j(j), w=(/w(2)/), c=(/j/))
    call pert_dens((/'EL'/), S, (/D/), (/F/), D_k(k), F_k(k), w=(/w(3)/), c=(/k/))

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_j(j)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_j(j)                  &
                    /),                      &
                     D_ij(i, j),             &
                     F_ij(i, j),             &
                     w=(/w(1), w(2)/),       &
                     c=(/i, j/))

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_k(k)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_k(k)                  &
                    /),                      &
                     D_ik(i, k),             &
                     F_ik(i, k),             &
                     w=(/w(1), w(3)/),       &
                     c=(/i, k/))

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_j(j),                 &
                     D_k(k)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_j(j),                 &
                     F_k(k)                  &
                    /),                      &
                     D_jk(j, k),             &
                     F_jk(j, k),             &
                     w=(/w(2), w(3)/),       &
                     c=(/j, k/))

    call pert_dens((/'EL', 'EL', 'EL'/),     &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_j(j),                 &
                     D_k(k),                 &
                     D_ij(i, j),             &
                     D_ik(i, k),             &
                     D_jk(j, k)              &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_j(j),                 &
                     F_k(k),                 &
                     F_ij(i, j),             &
                     F_ik(i, k),             &
                     F_jk(j, k)              &
                    /),                      &
                     D_ijk(i, j, k),         &
                     F_ijk(i, j, k),         &
                     w=(/w(1), w(2), w(3)/), &
                     c=(/i, j, k/))

    e = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_j,    &
                       D_k,    &
                       D_ij,   &
                       D_ik,   &
                       D_jk,   &
                       D_ijk   &
                      /),      &
                       e,      &
                       b=(/3, 3, 3/))

    D_i   = 0
    D_j   = 0
    D_k   = 0
    D_ij  = 0
    D_ik  = 0
    D_jk  = 0
    D_ijk = 0
    F_i   = 0
    F_j   = 0
    F_k   = 0
    F_ij  = 0
    F_ik  = 0
    F_jk  = 0
    F_ijk = 0

    gamma_element = -e(i, j, k, h)
    print '("@ gamma (w = "4f9.4")" a10 a2 a2 a2 f26.10)', real(-sum(w)), &
                                                           real(w(1)),    &
                                                           real(w(2)),    &
                                                           real(w(3)),    &
                                                           fc(h),         &
                                                           fc(i),         &
                                                           fc(j),         &
                                                           fc(k),         &
                                                           gamma_element

  end function

  function delta_element(S, D, F, w, h, i, j, k, l)

!   ----------------------------------------------------------------------------
    real(8)                     :: delta_element
    type(matrix), intent(inout) :: S, D, F
    integer,      intent(in)    :: h, i, j, k, l
    complex(8),   intent(in)    :: w(4)
!   ----------------------------------------------------------------------------
!   pascal triangle: 1 0th order matrix
!                    4 1st order matrices
!                    6 2nd
!                    4 3rd
!                    1 4th
    type(matrix)                :: D_i(3)
    type(matrix)                :: D_j(3)
    type(matrix)                :: D_k(3)
    type(matrix)                :: D_l(3)
    type(matrix)                :: D_ij(3, 3)
    type(matrix)                :: D_ik(3, 3)
    type(matrix)                :: D_il(3, 3)
    type(matrix)                :: D_jk(3, 3)
    type(matrix)                :: D_jl(3, 3)
    type(matrix)                :: D_kl(3, 3)
    type(matrix)                :: D_ijk(3, 3, 3)
    type(matrix)                :: D_ijl(3, 3, 3)
    type(matrix)                :: D_ikl(3, 3, 3)
    type(matrix)                :: D_jkl(3, 3, 3)
    type(matrix)                :: D_ijkl(3, 3, 3, 3)

    type(matrix)                :: F_i(3)
    type(matrix)                :: F_j(3)
    type(matrix)                :: F_k(3)
    type(matrix)                :: F_l(3)
    type(matrix)                :: F_ij(3, 3)
    type(matrix)                :: F_ik(3, 3)
    type(matrix)                :: F_il(3, 3)
    type(matrix)                :: F_jk(3, 3)
    type(matrix)                :: F_jl(3, 3)
    type(matrix)                :: F_kl(3, 3)
    type(matrix)                :: F_ijk(3, 3, 3)
    type(matrix)                :: F_ijl(3, 3, 3)
    type(matrix)                :: F_ikl(3, 3, 3)
    type(matrix)                :: F_jkl(3, 3, 3)
    type(matrix)                :: F_ijkl(3, 3, 3, 3)

    complex(8)                  :: e1(            3)
    complex(8)                  :: e2(         3, 3)
    complex(8)                  :: e3(      3, 3, 3)
    complex(8)                  :: e4(   3, 3, 3, 3)
    complex(8)                  :: e5(3, 3, 3, 3, 3)

    integer                     :: c1, c2, c3, c4
!   ----------------------------------------------------------------------------

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    do c1 = 1, 3
      D_i(c1) = tiny(0.0d0)*D
      D_j(c1) = tiny(0.0d0)*D
      D_k(c1) = tiny(0.0d0)*D
      D_l(c1) = tiny(0.0d0)*D
      F_i(c1) = tiny(0.0d0)*D
      F_j(c1) = tiny(0.0d0)*D
      F_k(c1) = tiny(0.0d0)*D
      F_l(c1) = tiny(0.0d0)*D
      do c2 = 1, 3
        D_ij(c1, c2) = tiny(0.0d0)*D
        D_ik(c1, c2) = tiny(0.0d0)*D
        D_il(c1, c2) = tiny(0.0d0)*D
        D_jk(c1, c2) = tiny(0.0d0)*D
        D_jl(c1, c2) = tiny(0.0d0)*D
        D_kl(c1, c2) = tiny(0.0d0)*D
        F_ij(c1, c2) = tiny(0.0d0)*D
        F_ik(c1, c2) = tiny(0.0d0)*D
        F_il(c1, c2) = tiny(0.0d0)*D
        F_jk(c1, c2) = tiny(0.0d0)*D
        F_jl(c1, c2) = tiny(0.0d0)*D
        F_kl(c1, c2) = tiny(0.0d0)*D
        do c3 = 1, 3
          D_ijk(c1, c2, c3) = tiny(0.0d0)*D
          D_ijl(c1, c2, c3) = tiny(0.0d0)*D
          D_ikl(c1, c2, c3) = tiny(0.0d0)*D
          D_jkl(c1, c2, c3) = tiny(0.0d0)*D
          F_ijk(c1, c2, c3) = tiny(0.0d0)*D
          F_ijl(c1, c2, c3) = tiny(0.0d0)*D
          F_ikl(c1, c2, c3) = tiny(0.0d0)*D
          F_jkl(c1, c2, c3) = tiny(0.0d0)*D
          do c4 = 1, 3
            D_ijkl(c1, c2, c3, c4) = tiny(0.0d0)*D
            F_ijkl(c1, c2, c3, c4) = tiny(0.0d0)*D
          end do
        end do
      end do
    end do

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_i(i), F_i(i), w=(/w(1)/), c=(/i/))
    e2 = 0.0d0
    call prop_energy((/'EL'/), (/D, D_i/), e2, b=(/3/))
    call print_polarizability(e2, (/w(1)/), i)

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_j(j), F_j(j), w=(/w(2)/), c=(/j/))
    e2 = 0.0d0
    call prop_energy((/'EL'/), (/D, D_j/), e2, b=(/3/))
    call print_polarizability(e2, (/w(2)/), j)

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_k(k), F_k(k), w=(/w(3)/), c=(/k/))
    e2 = 0.0d0
    call prop_energy((/'EL'/), (/D, D_k/), e2, b=(/3/))
    call print_polarizability(e2, (/w(3)/), k)

    call pert_dens((/'EL'/), S, (/D/), (/F/), D_l(l), F_l(l), w=(/w(4)/), c=(/l/))
    e2 = 0.0d0
    call prop_energy((/'EL'/), (/D, D_l/), e2, b=(/3/))
    call print_polarizability(e2, (/w(4)/), l)

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_j(j)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_j(j)                  &
                    /),                      &
                     D_ij(i, j),             &
                     F_ij(i, j),             &
                     w=(/w(1), w(2)/),       &
                     c=(/i, j/))
    e3 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_j,    &
                       D_ij    &
                      /),      &
                       e3,     &
                       b=(/3, 3/))
    call print_polarizability(e3, (/w(1), w(2)/), i, j)

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_k(k)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_k(k)                  &
                    /),                      &
                     D_ik(i, k),             &
                     F_ik(i, k),             &
                     w=(/w(1), w(3)/),       &
                     c=(/i, k/))
    e3 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_k,    &
                       D_ik    &
                      /),      &
                       e3,     &
                       b=(/3, 3/))
    call print_polarizability(e3, (/w(1), w(3)/), i, k)

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_l(l)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_l(l)                  &
                    /),                      &
                     D_il(i, l),             &
                     F_il(i, l),             &
                     w=(/w(1), w(4)/),       &
                     c=(/i, l/))
    e3 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_l,    &
                       D_il    &
                      /),      &
                       e3,     &
                       b=(/3, 3/))
    call print_polarizability(e3, (/w(1), w(4)/), i, l)

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_j(j),                 &
                     D_k(k)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_j(j),                 &
                     F_k(k)                  &
                    /),                      &
                     D_jk(j, k),             &
                     F_jk(j, k),             &
                     w=(/w(2), w(3)/),       &
                     c=(/j, k/))
    e3 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_j,    &
                       D_k,    &
                       D_jk    &
                      /),      &
                       e3,     &
                       b=(/3, 3/))
    call print_polarizability(e3, (/w(2), w(3)/), j, k)

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_j(j),                 &
                     D_l(l)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_j(j),                 &
                     F_l(l)                  &
                    /),                      &
                     D_jl(j, l),             &
                     F_jl(j, l),             &
                     w=(/w(2), w(4)/),       &
                     c=(/j, l/))
    e3 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_j,    &
                       D_l,    &
                       D_jl    &
                      /),      &
                       e3,     &
                       b=(/3, 3/))
    call print_polarizability(e3, (/w(2), w(4)/), j, l)

    call pert_dens((/'EL', 'EL'/),           &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_k(k),                 &
                     D_l(l)                  &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_k(k),                 &
                     F_l(l)                  &
                    /),                      &
                     D_kl(k, l),             &
                     F_kl(k, l),             &
                     w=(/w(3), w(4)/),       &
                     c=(/k, l/))
    e3 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_k,    &
                       D_l,    &
                       D_kl    &
                      /),      &
                       e3,     &
                       b=(/3, 3/))
    call print_polarizability(e3, (/w(3), w(4)/), k, l)

    call pert_dens((/'EL', 'EL', 'EL'/),     &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_j(j),                 &
                     D_k(k),                 &
                     D_ij(i, j),             &
                     D_ik(i, k),             &
                     D_jk(j, k)              &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_j(j),                 &
                     F_k(k),                 &
                     F_ij(i, j),             &
                     F_ik(i, k),             &
                     F_jk(j, k)              &
                    /),                      &
                     D_ijk(i, j, k),         &
                     F_ijk(i, j, k),         &
                     w=(/w(1), w(2), w(3)/), &
                     c=(/i, j, k/))
    e4 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_j,    &
                       D_k,    &
                       D_ij,   &
                       D_ik,   &
                       D_jk,   &
                       D_ijk   &
                      /),      &
                       e4,     &
                       b=(/3, 3, 3/))
    call print_polarizability(e4, (/w(1), w(2), w(3)/), i, j, k)

    call pert_dens((/'EL', 'EL', 'EL'/),     &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_j(j),                 &
                     D_l(l),                 &
                     D_ij(i, j),             &
                     D_il(i, l),             &
                     D_jl(j, l)              &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_j(j),                 &
                     F_l(l),                 &
                     F_ij(i, j),             &
                     F_il(i, l),             &
                     F_jl(j, l)              &
                    /),                      &
                     D_ijl(i, j, l),         &
                     F_ijl(i, j, l),         &
                     w=(/w(1), w(2), w(4)/), &
                     c=(/i, j, l/))
    e4 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_j,    &
                       D_l,    &
                       D_ij,   &
                       D_il,   &
                       D_jl,   &
                       D_ijl   &
                      /),      &
                       e4,     &
                       b=(/3, 3, 3/))
    call print_polarizability(e4, (/w(1), w(2), w(4)/), i, j, l)

    call pert_dens((/'EL', 'EL', 'EL'/),     &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_i(i),                 &
                     D_k(k),                 &
                     D_l(l),                 &
                     D_ik(i, k),             &
                     D_il(i, l),             &
                     D_kl(k, l)              &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_i(i),                 &
                     F_k(k),                 &
                     F_l(l),                 &
                     F_ik(i, k),             &
                     F_il(i, l),             &
                     F_kl(k, l)              &
                    /),                      &
                     D_ikl(i, k, l),         &
                     F_ikl(i, k, l),         &
                     w=(/w(1), w(3), w(4)/), &
                     c=(/i, k, l/))
    e4 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_k,    &
                       D_l,    &
                       D_ik,   &
                       D_il,   &
                       D_kl,   &
                       D_ikl   &
                      /),      &
                       e4,     &
                       b=(/3, 3, 3/))
    call print_polarizability(e4, (/w(1), w(3), w(4)/), i, k, l)

    call pert_dens((/'EL', 'EL', 'EL'/),     &
                     S,                      &
                   (/                        &
                     D,                      &
                     D_j(j),                 &
                     D_k(k),                 &
                     D_l(l),                 &
                     D_jk(j, k),             &
                     D_jl(j, l),             &
                     D_kl(k, l)              &
                    /),                      &
                   (/                        &
                     F,                      &
                     F_j(j),                 &
                     F_k(k),                 &
                     F_l(l),                 &
                     F_jk(j, k),             &
                     F_jl(j, l),             &
                     F_kl(k, l)              &
                    /),                      &
                     D_jkl(j, k, l),         &
                     F_jkl(j, k, l),         &
                     w=(/w(2), w(3), w(4)/), &
                     c=(/j, k, l/))
    e4 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_j,    &
                       D_k,    &
                       D_l,    &
                       D_jk,   &
                       D_jl,   &
                       D_kl,   &
                       D_jkl   &
                      /),      &
                       e4,     &
                       b=(/3, 3, 3/))
    call print_polarizability(e4, (/w(2), w(3), w(4)/), j, k, l)

    call pert_dens((/'EL', 'EL', 'EL', 'EL'/),     &
                     S,                            &
                   (/                              &
                     D,                            &
                     D_i(i),                       &        
                     D_j(j),                       &        
                     D_k(k),                       &        
                     D_l(l),                       &        
                     D_ij(i, j),                   &           
                     D_ik(i, k),                   &           
                     D_il(i, l),                   &           
                     D_jk(j, k),                   &           
                     D_jl(j, l),                   &           
                     D_kl(k, l),                   &
                     D_ijk(i, j, k),               &   
                     D_ijl(i, j, l),               &   
                     D_ikl(i, k, l),               &   
                     D_jkl(j, k, l)                &   
                    /),                            &
                   (/                              &
                     F,                            &
                     F_i(i),                       &        
                     F_j(j),                       &        
                     F_k(k),                       &        
                     F_l(l),                       &        
                     F_ij(i, j),                   &           
                     F_ik(i, k),                   &           
                     F_il(i, l),                   &           
                     F_jk(j, k),                   &           
                     F_jl(j, l),                   &           
                     F_kl(k, l),                   &
                     F_ijk(i, j, k),               &   
                     F_ijl(i, j, l),               &   
                     F_ikl(i, k, l),               &   
                     F_jkl(j, k, l)                &   
                    /),                            &
                     D_ijkl(i, j, k, l),           &
                     F_ijkl(i, j, k, l),           &
                     w=(/w(1), w(2), w(3), w(4)/), &
                     c=(/i, j, k, l/))
    e5 = 0.0d0
    call prop_energy((/'EL'/), &
                     (/        &
                       D,      &
                       D_i,    &
                       D_j,    &
                       D_k,    &
                       D_l,    &
                       D_ij,   &
                       D_ik,   &
                       D_il,   &
                       D_jk,   &
                       D_jl,   &
                       D_kl,   &
                       D_ijk,  &
                       D_ijl,  &
                       D_ikl,  &
                       D_jkl,  &
                       D_ijkl  &
                      /),      &
                       e5,     &
                       b=(/3, 3, 3, 3/))
    call print_polarizability(e5, (/w(1), w(2), w(3), w(4)/), i, j, k, l)

    delta_element = -e5(i, j, k, l, h)

    D_i    = 0
    D_j    = 0
    D_k    = 0
    D_l    = 0
    D_ij   = 0
    D_ik   = 0
    D_il   = 0
    D_jk   = 0
    D_jl   = 0
    D_kl   = 0
    D_ijk  = 0
    D_ijl  = 0
    D_ikl  = 0
    D_jkl  = 0
    D_ijkl = 0

    F_i    = 0
    F_j    = 0
    F_k    = 0
    F_l    = 0
    F_ij   = 0
    F_ik   = 0
    F_il   = 0
    F_jk   = 0
    F_jl   = 0
    F_kl   = 0
    F_ijk  = 0
    F_ijl  = 0
    F_ikl  = 0
    F_jkl  = 0
    F_ijkl = 0

  end function

  subroutine gamma(S, D, F, w)

    type(matrix),  intent(inout) :: S, D, F
    complex(8), intent(in)    :: w(4)
!   ----------------------------------------------------------------------------
    type(matrix)                 :: D_i(3), D_j(3), D_k(3), D_ij(3, 3), D_ik(3, 3), D_jk(3, 3), D_ijk(3, 3, 3)
    type(matrix)                 :: F_i(3), F_j(3), F_k(3), F_ij(3, 3), F_ik(3, 3), F_jk(3, 3), F_ijk(3, 3, 3)
    complex(8)                   :: e_ijkl(3, 3, 3, 3)
    integer                      :: i, j, k, l

    if (abs(sum(w)) > small_number) call quit('gamma: frequencies do not sum up to zero')

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

    print '("@ gamma (w = "4f9.4")")', real(w(1)), real(w(2)), real(w(3)), real(w(4))

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
            if (abs(e_ijkl(i, j, k, l)) > small_number) then
              print '("@" a10 a2 a2 a2 f26.10)', fc(i), fc(j), fc(k), fc(l), -real(e_ijkl(i, j, k, l))
            end if
          end do
        end do
      end do
    end do

  end subroutine

end module
