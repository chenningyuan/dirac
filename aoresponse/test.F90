module aoresponse_test_m

    use dirac_matrix_operations_m
    use aoresponse_xc_m
    use propint
    use pertdens

    implicit none

contains

  function coordinate(i) result(c)
    integer, intent(in) :: i
    character(1)        :: c

    select case (i)
      case (1)
        c = 'x'
      case (2)
        c = 'y'
      case (3)
        c = 'z'
    end select

  end function coordinate

  subroutine nr_electrons(S, D, F)

    type(matrix), intent(inout) :: S, D, F

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    write(*, *) 'electron count from S**D:', S**D

  end subroutine nr_electrons

  subroutine energy(S, D, F)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: G(1)
    real(kind_dp)               :: xc_energy
    character(4)                :: nop(0)

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    G(1) = 1.0d-100*S

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
      print '("@ -< "a1" > = " f16.10)', coordinate(i), -e_i(i)
    end do

  end subroutine mu

  subroutine alpha(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), F_i(3)
    real(kind_dp)               :: w
    real(kind_dp)               :: e_ij(3*3)
    integer                     :: i, j, ij

    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=(/w/))

    e_ij = 0.0d0
    call propenergy((/'EL'/), (/D, D_i/), e_ij)

    print '("@ alpha (w = "2f7.3")")', w, -w

    do i = 1, 3

      D_i(i) = 0
      F_i(i) = 0

      do j = 1, 3
        ij = i + 3*(j-1)
        print '("@" a10 a2 f20.10)', coordinate(i), coordinate(j), -e_ij(ij)
      end do
    end do

  end subroutine alpha

  subroutine beta(S, D, F, w)

    type(matrix), intent(inout) :: S, D, F
!   ----------------------------------------------------------------------------
    type(matrix)                :: D_i(3), D_j(3), D_ij(3*3)
    type(matrix)                :: F_i(3), F_j(3), F_ij(3*3)
    real(kind_dp)               :: w(3)
    real(kind_dp)               :: e_ijk(3*3*3)
    integer                     :: i, j, ij, k, ijk

    type(matrix)                :: M

    if (abs(sum(w)) > 1.0d-15) call quit('beta: frequencies do not sum up to zero')

    S = 1.0d0*S
    D = 1.0d0*D
    F = 1.0d0*F

    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_i, F_i, w=w(1:1))
    call pertdensfock((/'EL'/), S, (/D/), (/F/), D_j, F_j, w=w(2:2))

!   do i = 1, 3
!     call mat_print(D_i(i))
!     M = dag(D_j(i))
!     call mat_print(M)
!   end do
!   call quit('over and out')

    call pertdensfock((/'EL', 'EL'/), S, (/D, D_i, D_j/), (/F, F_i, F_j/), D_ij, F_ij, w=w(1:2))

!   write(*, *) 'blub after pertdensfock'
!   call mat_print(D_ij(3))
!   M = dag(D_ij(7))
!   call mat_print(M)

    e_ijk = 0.0d0
    call propenergy((/'EL'/), (/D, D_i, D_j, D_ij/), e_ijk, b=(/3, 3/))


    print '("@ beta (w = "3f7.3")")', w(1), w(2), w(3)

    do i = 1, 3

      D_i(i) = 0
      D_j(i) = 0
      F_i(i) = 0
      F_j(i) = 0

      do j = 1, 3
        ij = i + 3*(j-1)

        D_ij(ij) = 0
        F_ij(ij) = 0

        do k = 1, 3
          ijk = i + 3*(j-1) + 9*(k-1)
          print '("@" a10 a2 a2 f18.10)', coordinate(i), coordinate(j), coordinate(k), -e_ijk(ijk)
        end do
      end do
    end do

  end subroutine beta




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


subroutine gamma(S,D,F,w)
    type(Matrix):: S,D,Df(3),De(3),Dg(3),Dfe(3*3),Dfg(3*3),Deg(3*3),Dfeg(3*3*3)
    type(Matrix)::   F,Ff(3),Fe(3),Fg(3),Ffe(3*3),Ffg(3*3),Feg(3*3),Ffeg(3*3*3)
    real*8:: w(4),Effff(3*3*3*3); integer:: i,j,k
    S=S; D=D; F=F !track them
    if (abs(sum(w))>1d-15) call QUIT('PropTest/Hy2Polari: sum(w) should be zero!')
    call PertDensFock((/'EL'/),S,(/D/),(/F/),Df,Ff,w=w(1:1))
    call PertDensFock((/'EL'/),S,(/D/),(/F/),De,Fe,w=w(2:2))
    call PertDensFock((/'EL'/),S,(/D/),(/F/),Dg,Fg,w=w(3:3))
    call PertDensFock((/'EL','EL'/),S,(/D,Df,De/),(/F,Ff,Fe/),Dfe,Ffe,w=w(1:2))
    call PertDensFock((/'EL','EL'/),S,(/D,Df,Dg/),(/F,Ff,Fg/),Dfg,Ffg,w=(/w(1),w(3)/))
    call PertDensFock((/'EL','EL'/),S,(/D,De,Dg/),(/F,Fe,Fg/),Deg,Feg,w=w(2:3))
    call PertDensFock((/'EL','EL','EL'/),S,(/D,Df,De,Dg,Dfe,Dfg,Deg/) &
                 ,(/F,Ff,Fe,Fg,Ffe,Ffg,Feg/),Dfeg,Ffeg,w=w(1:3))
    Effff=0; call PropEnergy((/'EL'/),(/D,Df,De,Dg,Dfe,Dfg,Deg,Dfeg/) &
                            ,Effff,b=(/3,3,3/))
    print '("Second hyperpolarizability -Effff w="4f6.3/27(3f18.13/))',w,-Effff
    do i=1,3; Df(i)=0; De(i)=0; Dg(i)=0; Ff(i)=0; Fe(i)=0; Fg(i)=0; enddo
    do i=1,3*3; Dfe(i)=0; Dfg(i)=0; Deg(i)=0; Ffe(i)=0; Ffg(i)=0; Feg(i)=0; enddo
    do i=1,3*3*3; Dfeg(i)=0; Ffeg(i)=0; enddo
    print*
end subroutine

end module aoresponse_test_m
