module get_overlap

  use exacorr_mo
  use exacorr_datatypes

contains

   subroutine get_mo_overlap (ao_basis_bra,cspinor_bra,ao_basis_ket,cspinor_ket,SMO)

        implicit none

        complex(8), parameter              :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        type(basis_set_info_t),intent(in)  :: ao_basis_bra, ao_basis_ket
        type(cmo),             intent(in)  :: cspinor_bra, cspinor_ket
        complex(8),            intent(out) :: SMO(cspinor_bra%nmo,cspinor_ket%nmo)
        real(8), allocatable               :: SAO(:,:)   ! the AO overlap matrix
        complex(8), allocatable            :: SC(:,:) ! Overlap times MO coefficients
        integer                            :: i,j

        ! allocate auxilliaries
        allocate (SAO(cspinor_bra%nao,cspinor_ket%nao))
        allocate (SC(cspinor_bra%nao,cspinor_ket%nmo))
!        print*,'    MO overlap (nao basis, nshells basis, nao cspinor) bra:', &
!                                 ao_basis_bra%nao, ao_basis_bra%nshells, cspinor_bra%nao, &
!                                 ' ket:', &
!                                 ao_basis_ket%nao, ao_basis_ket%nshells, cspinor_ket%nao

        call  get_ao_overlap (ao_basis_bra,ao_basis_ket,SAO)

        SC  = matmul(SAO,cspinor_ket%coeff(:,:,1))                        ! multiply with alpha coefficients
        SMO = matmul(conjg(transpose(cspinor_bra%coeff(:,:,1))),SC)       ! alpha-alpha contribution
        SC  = matmul(SAO,cspinor_ket%coeff(:,:,2))                        ! multiply with beta coefficients
        SMO = SMO + matmul(conjg(transpose(cspinor_bra%coeff(:,:,2))),SC) ! add the beta-beta contribution

!       filter out occupation numbers slightly broken due to rounding.
!        do i=1,cspinor_bra%nmo
!           do j=1,cspinor_ket%nmo
!              if (REAL(SMO(i,j)) > 0.999999999999D0 .and. AIMAG(SMO(i,j)) < 1.0D-12) SMO(i,j) = ONE
!              if (ABS(REAL(SMO(i,j))) < 1.0D-12 .and. ABS(AIMAG(SMO(i,j))) < 1.0D-12) SMO(i,j) = ZERO
!           end do
!        end do

        deallocate(SAO,SC)
!     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') )

   end subroutine get_mo_overlap

   subroutine get_ao_overlap (ao_basis_bra,ao_basis_ket,ao_overlap)

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis_bra, ao_basis_ket
        real(8), intent(out)              :: ao_overlap(ao_basis_bra%nao,ao_basis_ket%nao)

        integer :: ksh, lsh,koff, loff
        integer :: lk,ll,nk,nl,k,l,kl
        real(8) :: ek,el,ck,cl,xk,yk,zk,xl,yl,zl

        integer :: nrkb
        real(8) :: skl(21*21) ! hardwired for maximum  l-value of 5 (s=1,p=3,d=6,f=10,g=15,h=21)
        real(8), parameter :: fijkl = 1.0
        integer :: ncentr=1
        real(8) :: atomic_data(1)

        loff = 0
        do lsh = 1, ao_basis_ket%nshells
            ll   =  ao_basis_ket%gtos(lsh)%orb_momentum
            el   =  ao_basis_ket%gtos(lsh)%exponent(1)
            xl   =  ao_basis_ket%gtos(lsh)%coord(1)
            yl   =  ao_basis_ket%gtos(lsh)%coord(2)
            zl   =  ao_basis_ket%gtos(lsh)%coord(3)
            cl   =  ao_basis_ket%gtos(lsh)%coefficient(1)
            koff = 0
            do ksh = 1, ao_basis_bra%nshells
               lk   =  ao_basis_bra%gtos(ksh)%orb_momentum
               ek   =  ao_basis_bra%gtos(ksh)%exponent(1)
               xk   =  ao_basis_bra%gtos(ksh)%coord(1)
               yk   =  ao_basis_bra%gtos(ksh)%coord(2)
               zk   =  ao_basis_bra%gtos(ksh)%coord(3)
               ck   =  ao_basis_bra%gtos(ksh)%coefficient(1)
               call  interest_overlap(.false.,   &
                              skl,nk,nl,nrkb,   &
                              lk,ek,xk,yk,zk,ck,&
                              ll,el,xl,yl,zl,cl,&
                              ncentr,atomic_data)
               kl = 1
               do l = 1, nl
                  do k = 1, nk
                      ao_overlap(koff+k,loff+l) = skl(kl)
                      kl = kl + 1
                  end do
               end do

               koff = koff + nk
            end do
            loff = loff + nl
        end do

   end subroutine get_ao_overlap

end module get_overlap
