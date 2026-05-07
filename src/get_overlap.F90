module get_overlap

   use rose_ao
   use rose_mo
   use interest
   use linear_algebra

contains

   subroutine get_mo_overlap (ao_basis_bra,spinor_bra,ao_basis_ket,spinor_ket,SMO)

        implicit none

        complex(8), parameter              :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        type(basis_set_info_t),intent(in)  :: ao_basis_bra, ao_basis_ket
        type(mo_basis),     intent(in)     :: spinor_bra, spinor_ket
        real(8), allocatable               :: SAO(:,:)   ! the AO overlap matrix
        complex(8), allocatable            :: SC(:,:) ! Overlap times MO coefficients
        complex(8), allocatable            :: tmp(:,:)
        real(8),    allocatable            :: SAO_tmp(:,:,:)
        real(8),    allocatable            :: spinor_bra_coeff_bar(:,:,:)
        real(8), intent(inout)             :: SMO(:,:,:)

        if (spinor_bra%type .ne. spinor_ket%type) stop 'get_overlap: not same type'

        allocate (SAO(spinor_bra%nao,spinor_ket%nao))

        SMO = 0.D0
        
        if (ao_basis_bra%basis_type == 1 .and. ao_basis_ket%basis_type == 1) then
            call get_ao_overlap (ao_basis_bra,ao_basis_ket,SAO)
        elseif (ao_basis_bra%basis_type == 2 .and. ao_basis_ket%basis_type == 2) then
            call get_ao_overlap_ADF (ao_basis_bra,ao_basis_ket,SAO)
        else
            error stop 'cannot compute this overlap type'
        endif

        if (spinor_bra%type .eq. 1) then
           allocate(SAO_tmp(spinor_bra%nao,spinor_ket%nao,1))
           allocate(spinor_bra_coeff_bar(spinor_bra%nmo,spinor_bra%nao,1))
           SAO_tmp = 0.D0
           SAO_tmp(:,:,1) = SAO
           call transpose_conjg(spinor_bra%coeff_r,spinor_bra_coeff_bar)
           call matrix_mul(spinor_bra_coeff_bar,SAO_tmp,spinor_ket%coeff_r,SMO)
           deallocate(SAO_tmp)
           deallocate(spinor_bra_coeff_bar)
        else if (spinor_bra%type .eq. 2) then
           allocate (SC(spinor_bra%nao,spinor_ket%nmo))
           SC  = matmul(SAO,spinor_ket%coeff_r(:,:,1))                               ! multiply with alpha coefficients
           SMO(:,:,1) = matmul(transpose(spinor_bra%coeff_r(:,:,1)),SC)              ! alpha-alpha contribution
           SC  = matmul(SAO,spinor_ket%coeff_r(:,:,2))                               ! multiply with beta coefficients
           SMO(:,:,1) = SMO(:,:,1) + matmul(transpose(spinor_bra%coeff_r(:,:,2)),SC) ! add the beta-beta contribution
           deallocate(SC)
        else if (spinor_bra%type .eq. 3) then
           allocate (tmp(spinor_bra%nmo,spinor_ket%nmo))
           allocate (SC(spinor_bra%nao,spinor_ket%nmo))
           SC  = matmul(SAO,spinor_ket%coeff_c(:,:,1))                        ! multiply with alpha coefficients
           tmp = matmul(conjg(transpose(spinor_bra%coeff_c(:,:,1))),SC)       ! alpha-alpha contribution
           SC  = matmul(SAO,spinor_ket%coeff_c(:,:,2))                        ! multiply with beta coefficients
           tmp = tmp + matmul(conjg(transpose(spinor_bra%coeff_c(:,:,2))),SC) ! add the beta-beta contribution
           SMO(:,:,1) = real(tmp(:,:))
           SMO(:,:,2) = dimag(tmp(:,:))
           deallocate(tmp)
           deallocate(SC)
        else if (spinor_bra%type .eq. 4) then
           allocate(SAO_tmp(spinor_bra%nao,spinor_ket%nao,4))
           allocate(spinor_bra_coeff_bar(spinor_bra%nmo,spinor_bra%nao,4))
           SAO_tmp = 0.D0
           SAO_tmp(:,:,1) = SAO
           call transpose_conjg(spinor_bra%coeff_q,spinor_bra_coeff_bar)
           call matrix_mul(spinor_bra_coeff_bar,SAO_tmp,spinor_ket%coeff_q,SMO)
           deallocate(SAO_tmp,spinor_bra_coeff_bar)
        else
           error stop "Error in the dimension defining the algebra."
        end if

        deallocate(SAO)

   end subroutine get_mo_overlap

   subroutine get_ao_overlap (ao_basis_bra,ao_basis_ket,ao_overlap)

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis_bra, ao_basis_ket
        real(8), intent(out)              :: ao_overlap(ao_basis_bra%nao,ao_basis_ket%nao)

        integer :: ksh, lsh,koff, loff
        integer :: kprim, lprim
        integer :: lk,ll,nk,nl,k,l,kl
        real(8) :: ek,el,ck,cl,xk,yk,zk,xl,yl,zl

        integer :: nrkb
        real(8) :: skl(36*36) ! hardwired for maximum  l-value of 7 (s=1,p=3,d=6,f=10,g=15,h=21,...)
        real(8), parameter :: fijkl = 1.0

        ao_overlap = 0.D0
        loff = 0
        do lsh = 1, ao_basis_ket%nshells
            ll =  ao_basis_ket%gtos(lsh)%orb_momentum
            xl =  ao_basis_ket%gtos(lsh)%coord(1)
            yl =  ao_basis_ket%gtos(lsh)%coord(2)
            zl =  ao_basis_ket%gtos(lsh)%coord(3)
            do lprim = 1, ao_basis_ket%gtos(lsh)%n_primitives
                el =  ao_basis_ket%gtos(lsh)%exponent(lprim)
                cl =  ao_basis_ket%gtos(lsh)%coefficient(lprim)
                koff = 0
                do ksh = 1, ao_basis_bra%nshells
                   lk =  ao_basis_bra%gtos(ksh)%orb_momentum
                   xk =  ao_basis_bra%gtos(ksh)%coord(1)
                   yk =  ao_basis_bra%gtos(ksh)%coord(2)
                   zk =  ao_basis_bra%gtos(ksh)%coord(3)
                   do kprim = 1, ao_basis_bra%gtos(ksh)%n_primitives
                       ek =  ao_basis_bra%gtos(ksh)%exponent(kprim)
                       ck =  ao_basis_bra%gtos(ksh)%coefficient(kprim)
                       call  interest_overlap(skl,nk,nl,nrkb,   &
                                              lk,ek,xk,yk,zk,ck,&
                                              ll,el,xl,yl,zl,cl )
                       kl = 1
                       do l = 1, nl
                          do k = 1, nk
                              ao_overlap(koff+k,loff+l) = ao_overlap(koff+k,loff+l) + skl(kl)
                              kl = kl + 1
                          end do
                       end do
                   end do
                   koff = koff + nk
                end do
            end do
            loff = loff + nl
        end do

   end subroutine get_ao_overlap

end module get_overlap
