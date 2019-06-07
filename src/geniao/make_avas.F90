module make_avas

contains

   subroutine avas(A,Cocc,nmo1,nmo2,nocc)

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        integer, intent(in)     :: nmo1,nmo2,nocc
        complex(8), intent(in)      :: A(nmo1,nmo2), Cocc(nmo1,nocc)
        complex(8), allocatable :: Proj(:,:), Socc(:,:), Svir(:,:), Cvir(:,:)
        integer                 :: nvir, i, j
       
        nvir = nmo1 - nocc
        print*,
        print*, "******************* Entering AVAS *******************"
        ! Compute the projector (Eq. 15 of arXiv:1701.07862v1): P = ST(T†ST)^(−1)T†S. with T = A, S = S11, (T†ST)^(−1) = 1 because the IAO are symetrically orthogonalized, and S11 = 1 as well.
        allocate(Proj(nmo1,nmo1))
        Proj = matmul(A,conjg(transpose(A)))

        ! Compute the projected overlap matrices S_occ = C†_occ P C_occ and S_vir = C†_vir P Cvir
        allocate(Socc(nocc,nocc))
        Socc = matmul(conjg(transpose(Cocc)),matmul(Proj,Cocc))
!        do i=1,nocc
!          write(*,'(100G10.2)') (REAL(Socc(i,j)),j=1,nocc)
!        enddo
        nvir = nmo1 - nocc
        allocate(Svir(nvir,nvir))
        allocate(Cvir(nmo1,nvir))
        ! In the MO basis, Cvir and C_occ are both diagonal, with 1 on virtual and occupied orbital indices, respectively.
        Cvir = ZERO
        forall (i=1:nvir) Cvir(i,i) = ONE
        Svir = matmul(conjg(transpose(Cvir)),matmul(Proj,Cvir))
!        do i=1,nvir
!          write(*,'(100G10.2)') (REAL(Svir(i,j)),j=1,nvir)
!        enddo

        deallocate(Svir,Socc,Proj,Cvir)
   end subroutine avas

end module make_avas
