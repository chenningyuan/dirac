module aoresponse_test_m

    use matrix_mod
    use aoresponse_xc_m
    use propint
    use pertdens

    implicit none

    character(1)  :: fc(3) = (/'x', 'y', 'z'/)
    real(kind_dp) :: small = 1.0d-10

contains

  subroutine nr_electrons(S, D)

    type(matrix), intent(inout) :: S, D

    S = 1.0d0*S
    D = 1.0d0*D

    write(*, *) 'electron count from S**D:', S**D

  end subroutine nr_electrons

  subroutine energy(D, F)

    type(matrix), intent(inout) :: D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: G(1)
    real(kind_dp)               :: xc_energy
    character(4)                :: nop(0)

    D = 1.0d0*D
    F = 1.0d0*F

    G(1) = 1.0d-100*D

    call propfock(nop, (/D, D/), G)
    write(*, *) 'electronic energy from F**D - 0.5d0*G(1)**D:', F**D - 0.5d0*G(1)**D

    G(1) = 0

    xc_energy = 0.0d0
!   call xc_contribution(0, (/D/), (/F/), xc_energy)
    write(*, *) 'xc energy:', xc_energy

  end subroutine energy

  subroutine mu(S, D, F)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    real(kind_dp)               :: e_i(3)
    integer                     :: i, j

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    e_i = 0.0d0
    call propenergy((/'EL'/), (/D/), e_i)

    do i = 1, 3
      print '("@ -< "a1" > = " f16.10)', fc(i), -e_i(i)
    end do

  end subroutine mu


  subroutine alpha(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3)
    type(matrix)                :: F_i(3)
    real(kind_dp)               :: w(2)
    real(kind_dp)               :: e_ij(3, 3)
    integer                     :: i, j

    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(2:2))

    e_ij = 0.0d0
    call propenergy((/'EL'/), (/D, D_i/), e_ij)

    print '("@ alpha (w = "2f7.3")")', w(1), w(2)

    do i = 1, 3

      D_i(i) = 0
      F_i(i) = 0

      do j = 1, 3
        if (abs(e_ij(i, j)) > small) then
          print '("@" a10 a2 f20.10)', fc(i), fc(j), -e_ij(i, j)
        end if
      end do
    end do

  end subroutine alpha

  subroutine beta(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), D_j(3), D_ij(3, 3)
    type(matrix)                :: F_i(3), F_j(3), F_ij(3, 3)
    real(kind_dp)               :: w(3)
    real(kind_dp)               :: e_ijk(3, 3, 3)
    integer                     :: i, j, k

    if (abs(sum(w)) > small) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))

    call pertdensfock((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))

    e_ijk = 0.0d0
    call propenergy((/'EL'/), (/D, D_i, D_j, D_ij/), e_ijk, b=(/3, 3/))

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
            print '("@" a10 a2 a2 f18.10)', fc(i), fc(j), fc(k), -e_ijk(i, j, k)
          end if
        end do
      end do
    end do

  end subroutine beta

  subroutine gamma(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), D_j(3), D_k(3), D_ij(3, 3), D_ik(3, 3), D_jk(3, 3), D_ijk(3, 3, 3)
    type(matrix)                :: F_i(3), F_j(3), F_k(3), F_ij(3, 3), F_ik(3, 3), F_jk(3, 3), F_ijk(3, 3, 3)
    real(kind_dp)               :: w(4)
    real(kind_dp)               :: e_ijkl(3, 3, 3, 3)
    integer                     :: i, j, k, l

    if (abs(sum(w)) > small) call quit('gamma: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))
    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_k, F_k, w=w(3:3))

    call pertdensfock((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))
    call pertdensfock((/'EL', 'EL'/), S, (/D, D_i, D_k/), (/F, F_i, F_k/), D_ik, F_ik, w=(/w(1), w(3)/))
    call pertdensfock((/'EL', 'EL'/), S, (/D, D_j, D_k/), (/F, F_j, F_k/), D_jk, F_jk, w=w(2:3))

    call pertdensfock((/'EL', 'EL', 'EL'/), S, (/D, D_i, D_j, D_k, D_ij, D_ik, D_jk/), &
                                               (/F, F_i, F_j, F_k, F_ij, F_ik, F_jk/), &
                                               D_ijk, F_ijk, w=w(1:3))

    e_ijkl = 0.0d0
    call propenergy((/'EL'/), (/D, D_i, D_j, D_k, D_ij, D_ik, D_jk, D_ijk/), e_ijkl, b=(/3, 3, 3/))

    print '("@ gamma (w = "4f7.3")")', w(1), w(2), w(3), w(4)

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
              print '("@" a10 a2 a2 a2 f16.10)', fc(i), fc(j), fc(k), fc(l), -e_ijkl(i, j, k, l)
            end if
          end do
        end do
      end do
    end do

  end subroutine gamma







subroutine beta_debug(S,D,F,w)
    type(Matrix):: S,D,Df(3),De(3),Dfe(3*3), T, M
    type(Matrix)::   F,Ff(3),Fe(3),Ffe(3*3), G(1), H(1)
    real*8:: w(3),Efff(3*3*3); integer:: i,j, ij
    character*(4) :: NOP(0)
    S=S; D=D; F=F !track them
    if (abs(sum(w))>1d-15) call QUIT('PropTest/HyPolari: sum(w) should be zero!')
    call PertDensFock((/'EL'/),S,(/D/),(/F/),Df,Ff,w=w(1:1))
    do i = 1, 3

       H(1) = 1.0d-50*D
       call PropFock((/'EL'/),(/D/),H,c=(/i/))
       G(1) = 1.0d-50*D
       call PropFock(NOP,(/D,Df(i)/),G)

       T = H(1) + G(1) - Ff(i)
       write(*, *) 'debug_norm 1:', T**T

       T = Df(i)*S*D + D*S*Df(i) - Df(i)
       write(*, *) 'debug_norm 2:', T**T

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
       write(*, *) 'debug_norm 4:', T**T

    end do
!   call quit('over and out')
!   call PertDensFock((/'EL'/),S,(/D/),(/F/),Df,Ff,w=w(1:1))
    call PertDensFock((/'EL'/),S,(/D/),(/F/),De,Fe,w=w(2:2))
    call PertDensFock((/'EL','EL'/),S,(/D,Df,De/),(/F,Ff,Fe/),Dfe,Ffe,w=w(1:2))
    do i = 1, 3
       do j = 1, 3
          ij = i + 3*(j-1)
          T = Dfe(ij)*S*D + De(j)*S*Df(i) + Df(i)*S*De(j) + D*S*Dfe(ij) - Dfe(ij)
       write(*, *) 'debug_norm 5:', T**T
      end do
    end do

    Efff=0; call PropEnergy((/'EL'/),(/D,Df,De,Dfe/),Efff,b=(/3,3/))
    print '("First hyperpolarizability -Efff w="3f6.3/9(3f18.13/))',w,-Efff
    do i=1,3; Df(i)=0; De(i)=0; Ff(i)=0; Fe(i)=0; enddo
    do i=1,3*3; Dfe(i)=0; Ffe(i)=0; enddo
    print*
end subroutine


end module aoresponse_test_m
