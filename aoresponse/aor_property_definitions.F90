module aor_property_definitions

    use matrix_defop
    use matrix_genop
    use decompmod
    use aoresponse_output
    use dirac_interface
    use prop_contribs
    use rsp_equations
    use aor_cfg

    implicit none

    public gamma_element
    public delta_element

    private

    character(1)  :: fc(3) = (/'x', 'y', 'z'/)
    real(8) :: small_number = 1.0d-10

contains

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
    type(decompitem)            :: decomp
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

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_i(i), F_i(i), comp=(/i/), freq=(/w(1)/))
    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_j(j), F_j(j), comp=(/j/), freq=(/w(2)/))
    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_k(k), F_k(k), comp=(/k/), freq=(/w(3)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j)              &
                  /),                  &
                   D_ij(i, j),         &
                   F_ij(i, j),         &
                   comp=(/i, j/),      &
                   freq=(/w(1), w(2)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_k(k)              &
                  /),                  &
                   D_ik(i, k),         &
                   F_ik(i, k),         &
                   comp=(/i, k/),      &
                   freq=(/w(1), w(3)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_k(k)              &
                  /),                  &
                   D_jk(j, k),         &
                   F_jk(j, k),         &
                   comp=(/j, k/),      &
                   freq=(/w(2), w(3)/))

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j),             &
                   D_k(k),             &
                   D_ij(i, j),         &
                   D_ik(i, k),         &
                   D_jk(j, k)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j),             &
                   F_k(k),             &
                   F_ij(i, j),         &
                   F_ik(i, k),         &
                   F_jk(j, k)          &
                  /),                  &
                   D_ijk(i, j, k),     &
                   F_ijk(i, j, k),     &
                   comp=(/i, j, k/),   &
                   freq=(/w(1), w(2), w(3)/))

    e = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijk/), (/3, 3, 3, 3/), e)

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

    gamma_element = -e(h, i, j, k)
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
    type(decompitem)            :: decomp
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

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_i(i), F_i(i), comp=(/i/), freq=(/w(1)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_i/), shape(e2), e2)
    call print_polarizability(e2, (/w(1)/), i)

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_j(j), F_j(j), comp=(/j/), freq=(/w(2)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_j/), shape(e2), e2)
    call print_polarizability(e2, (/w(2)/), j)

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_k(k), F_k(k), comp=(/k/), freq=(/w(3)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_k/), shape(e2), e2)
    call print_polarizability(e2, (/w(3)/), k)

    call pert_dens(decomp, S, (/'EL'/), (/1/), (/D/), (/F/), D_l(l), F_l(l), comp=(/l/), freq=(/w(4)/))
    e2 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_l/), shape(e2), e2)
    call print_polarizability(e2, (/w(4)/), l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j)              &
                  /),                  &
                   D_ij(i, j),         &
                   F_ij(i, j),         &
                   comp=(/i, j/),      &
                   freq=(/w(1), w(2)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ij/), shape(e3), e3)
    call print_polarizability(e3, (/w(1), w(2)/), i, j)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_k(k)              &
                  /),                  &
                   D_ik(i, k),         &
                   F_ik(i, k),         &
                   comp=(/i, k/),      &
                   freq=(/w(1), w(3)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ik/), shape(e3), e3)
    call print_polarizability(e3, (/w(1), w(3)/), i, k)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_l(l)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_l(l)              &
                  /),                  &
                   D_il(i, l),         &
                   F_il(i, l),         &
                   comp=(/i, l/),      &
                   freq=(/w(1), w(4)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_il/), shape(e3), e3)
    call print_polarizability(e3, (/w(1), w(4)/), i, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_k(k)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_k(k)              &
                  /),                  &
                   D_jk(j, k),         &
                   F_jk(j, k),         &
                   comp=(/j, k/),      &
                   freq=(/w(2), w(3)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_jk/), shape(e3), e3)
    call print_polarizability(e3, (/w(2), w(3)/), j, k)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_l(l)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_l(l)              &
                  /),                  &
                   D_jl(j, l),         &
                   F_jl(j, l),         &
                   comp=(/j, l/),      &
                   freq=(/w(2), w(4)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_jl/), shape(e3), e3)
    call print_polarizability(e3, (/w(2), w(4)/), j, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL'/),       &
                 (/1, 1/),             &
                 (/                    &
                   D,                  &
                   D_k(k),             &
                   D_l(l)              &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_k(k),             &
                   F_l(l)              &
                  /),                  &
                   D_kl(k, l),         &
                   F_kl(k, l),         &
                   comp=(/k, l/),      &
                   freq=(/w(3), w(4)/))
    e3 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_kl/), shape(e3), e3)
    call print_polarizability(e3, (/w(3), w(4)/), k, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j),             &
                   D_k(k),             &
                   D_ij(i, j),         &
                   D_ik(i, k),         &
                   D_jk(j, k)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j),             &
                   F_k(k),             &
                   F_ij(i, j),         &
                   F_ik(i, k),         &
                   F_jk(j, k)          &
                  /),                  &
                   D_ijk(i, j, k),     &
                   F_ijk(i, j, k),     &
                   comp=(/i, j, k/),   &
                   freq=(/w(1), w(2), w(3)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijk/), shape(e4), e4)
    call print_polarizability(e4, (/w(1), w(2), w(3)/), i, j, k)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_j(j),             &
                   D_l(l),             &
                   D_ij(i, j),         &
                   D_il(i, l),         &
                   D_jl(j, l)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_j(j),             &
                   F_l(l),             &
                   F_ij(i, j),         &
                   F_il(i, l),         &
                   F_jl(j, l)          &
                  /),                  &
                   D_ijl(i, j, l),     &
                   F_ijl(i, j, l),     &
                   comp=(/i, j, l/),   &
                   freq=(/w(1), w(2), w(4)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijl/), shape(e4), e4)
    call print_polarizability(e4, (/w(1), w(2), w(4)/), i, j, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_i(i),             &
                   D_k(k),             &
                   D_l(l),             &
                   D_ik(i, k),         &
                   D_il(i, l),         &
                   D_kl(k, l)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_i(i),             &
                   F_k(k),             &
                   F_l(l),             &
                   F_ik(i, k),         &
                   F_il(i, l),         &
                   F_kl(k, l)          &
                  /),                  &
                   D_ikl(i, k, l),     &
                   F_ikl(i, k, l),     &
                   comp=(/i, k, l/),   &
                   freq=(/w(1), w(3), w(4)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ikl/), shape(e4), e4)
    call print_polarizability(e4, (/w(1), w(3), w(4)/), i, k, l)

    call pert_dens(decomp,             &
                   S,                  &
                 (/'EL', 'EL', 'EL'/), &
                 (/1, 1, 1/),          &
                 (/                    &
                   D,                  &
                   D_j(j),             &
                   D_k(k),             &
                   D_l(l),             &
                   D_jk(j, k),         &
                   D_jl(j, l),         &
                   D_kl(k, l)          &
                  /),                  &
                 (/                    &
                   F,                  &
                   F_j(j),             &
                   F_k(k),             &
                   F_l(l),             &
                   F_jk(j, k),         &
                   F_jl(j, l),         &
                   F_kl(k, l)          &
                  /),                  &
                   D_jkl(j, k, l),     &
                   F_jkl(j, k, l),     &
                   comp=(/j, k, l/),   &
                   freq=(/w(2), w(3), w(4)/))
    e4 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_jkl/), shape(e4), e4)
    call print_polarizability(e4, (/w(2), w(3), w(4)/), j, k, l)

    call pert_dens(decomp,                   &
                   S,                        &
                 (/'EL', 'EL', 'EL', 'EL'/), &
                 (/1, 1, 1, 1/),             &
                 (/                          &
                   D,                        &
                   D_i(i),                   &        
                   D_j(j),                   &        
                   D_k(k),                   &        
                   D_l(l),                   &        
                   D_ij(i, j),               &           
                   D_ik(i, k),               &           
                   D_il(i, l),               &           
                   D_jk(j, k),               &           
                   D_jl(j, l),               &           
                   D_kl(k, l),               &
                   D_ijk(i, j, k),           &   
                   D_ijl(i, j, l),           &   
                   D_ikl(i, k, l),           &   
                   D_jkl(j, k, l)            &   
                  /),                        &
                 (/                          &
                   F,                        &
                   F_i(i),                   &        
                   F_j(j),                   &        
                   F_k(k),                   &        
                   F_l(l),                   &        
                   F_ij(i, j),               &           
                   F_ik(i, k),               &           
                   F_il(i, l),               &           
                   F_jk(j, k),               &           
                   F_jl(j, l),               &           
                   F_kl(k, l),               &
                   F_ijk(i, j, k),           &   
                   F_ijl(i, j, l),           &   
                   F_ikl(i, k, l),           &   
                   F_jkl(j, k, l)            &   
                  /),                        &
                   D_ijkl(i, j, k, l),       &
                   F_ijkl(i, j, k, l),       &
                   comp=(/i, j, k, l/),      &
                   freq=(/w(1), w(2), w(3), w(4)/))
    e5 = 0.0d0
    call prop_oneave(S, (/'EL'/), (/D_ijkl/), shape(e5), e5)
    call print_polarizability(e5, (/w(1), w(2), w(3), w(4)/), i, j, k, l)

    delta_element = -e5(h, i, j, k, l)

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
