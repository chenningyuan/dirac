module get_overlap

   use rose_ao
   use rose_mo
   use interest
   use linear_algebra
   use rose_utils, only: nfunctions
   use cartesian_to_spherical, only: get_sphcom
   use slater_overlap, only: sto_overlap_matrix

   type sph_transform_t
      real(8), allocatable :: matrix(:,:)
   end type sph_transform_t

contains

    subroutine compute_all_overlaps (ao_basis_1,ao_basis_2,spinor_original,spinor_fragment,S11,S12,S21,S22)

        implicit none
        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis of the supermolecule
        type(basis_set_info_t),intent(in)  :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis),intent(in)          :: spinor_original    ! MO-coefficients of the supermolecule
        type(mo_basis),intent(in)          :: spinor_fragment(:) ! MO-coefficients of the fragments
        real(8),intent(inout)              :: S11(:,:,:),S21(:,:,:),S12(:,:,:),S22(:,:,:) ! Overlap matrices
        real(8), allocatable               :: S_block(:,:,:),transconjg_S_block(:,:,:)    ! Subblocks of the overlap matrices
        integer                            :: imo, jmo, i, j, nmo1, nmoi, nmoj, rcq, n_fragments

        nmo1        = size(S11,1)
        rcq         = size(S11,3)
        n_fragments = size(spinor_fragment)

        call get_mo_overlap (ao_basis_1,spinor_original,ao_basis_1,spinor_original,S11)

        imo = 0
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           jmo = imo
           do j = i, n_fragments
              nmoj = spinor_fragment(j)%nmo
              allocate (S_block(nmoi,nmoj,rcq))
              call get_mo_overlap (ao_basis_2(i),spinor_fragment(i),ao_basis_2(j),spinor_fragment(j),S_block)
              S22(imo+1:imo+nmoi,jmo+1:jmo+nmoj,:) = S_block
              if (i /= j) then
                 allocate (transconjg_S_block(nmoj,nmoi,rcq))
                 call transpose_conjg(S_block,transconjg_S_block)
                 S22(jmo+1:jmo+nmoj,imo+1:imo+nmoi,:) = transconjg_S_block
                 deallocate (transconjg_S_block)
              end if
              jmo = jmo + nmoj
              deallocate (S_block)
           enddo
           imo = imo + nmoi
        enddo

        imo = 0
        do i = 1, n_fragments
           nmoi = spinor_fragment(i)%nmo
           allocate (S_block(nmoi,nmo1,rcq))
           call get_mo_overlap (ao_basis_2(i),spinor_fragment(i),ao_basis_1,spinor_original,S_block)
           S21(imo+1:imo+nmoi,:,:) = S_block
           imo = imo + nmoi
           deallocate (S_block)
        enddo

        call transpose_conjg(S21,S12)

   end subroutine compute_all_overlaps

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
            ! A Slater basis reaches us in one of two forms. Read from file it arrives as
            ! stos(:), shell by shell, and the overlap is computed here. Inside AMS the
            ! basis is an ADF object instead, stos(:) is null, and ADF computes it.
            if (associated(ao_basis_bra%stos) .and. associated(ao_basis_ket%stos)) then
               if (ao_basis_bra%basis_angular /= ao_basis_ket%basis_angular) then
                  error stop 'get_overlap: mixed angular conventions in a Slater overlap'
               end if
               call sto_overlap_matrix (ao_basis_bra%stos,ao_basis_ket%stos,SAO, &
                                        angular=ao_basis_bra%basis_angular)
            else
               call get_ao_overlap_ADF (ao_basis_bra,ao_basis_ket,SAO)
            end if
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
        type(sph_transform_t), allocatable :: trans_bra(:), trans_ket(:)

        ao_overlap = 0.D0
        call init_overlap_transformations (ao_basis_bra,ao_basis_ket,trans_bra,trans_ket)
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
                       call transformed_interest_overlap(skl,nk,nl,ao_basis_bra%basis_angular,ao_basis_ket%basis_angular, &
                                                         lk,ek,xk,yk,zk,ck,                                               &
                                                         ll,el,xl,yl,zl,cl,                                               &
                                                         trans_bra,trans_ket)
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

   subroutine transformed_interest_overlap (skl,nk,nl,basis_angular_bra,basis_angular_ket, &
                                            lk,ek,xk,yk,zk,ck,                            &
                                            ll,el,xl,yl,zl,cl,                            &
                                            trans_bra,trans_ket)

        implicit none

        real(8),intent(out)              :: skl(:)
        integer,intent(out)              :: nk, nl
        integer,intent(in)               :: basis_angular_bra, basis_angular_ket
        integer,intent(in)               :: lk, ll
        real(8),intent(in)               :: ek, el, ck, cl, xk, yk, zk, xl, yl, zl
        type(sph_transform_t),intent(in) :: trans_bra(:), trans_ket(:)

        integer                          :: nk_cart, nl_cart, nrkb
        real(8)                          :: skl_cart(36*36)
        real(8),allocatable              :: S_cart(:,:), S_bra(:,:), S_out(:,:)

        call interest_overlap(skl_cart,nk_cart,nl_cart,nrkb, &
                              lk,ek,xk,yk,zk,ck,             &
                              ll,el,xl,yl,zl,cl)

        if (nk_cart /= nfunctions(lk,1)) error stop 'unexpected bra shell size from interest_overlap'
        if (nl_cart /= nfunctions(ll,1)) error stop 'unexpected ket shell size from interest_overlap'

        nk = nfunctions(lk,basis_angular_bra)
        nl = nfunctions(ll,basis_angular_ket)

        if (size(skl) < nk*nl) error stop 'insufficient overlap scratch space'

        allocate(S_cart(nk_cart,nl_cart))
        S_cart = reshape(skl_cart(1:nk_cart*nl_cart),(/nk_cart,nl_cart/))

        if (basis_angular_bra == 2) then
           S_bra = matmul(transpose(trans_bra(lk)%matrix),S_cart)
        else
           S_bra = S_cart
        end if

        if (basis_angular_ket == 2) then
           S_out = matmul(S_bra,trans_ket(ll)%matrix)
        else
           S_out = S_bra
        end if

        skl = 0.D0
        skl(1:nk*nl) = reshape(S_out,(/nk*nl/))

   end subroutine transformed_interest_overlap

   subroutine init_overlap_transformations (ao_basis_bra,ao_basis_ket,trans_bra,trans_ket)

        implicit none

        type(basis_set_info_t),intent(in)              :: ao_basis_bra, ao_basis_ket
        type(sph_transform_t),allocatable,intent(out)  :: trans_bra(:), trans_ket(:)

        allocate(trans_bra(max_orb_momentum(ao_basis_bra)))
        allocate(trans_ket(max_orb_momentum(ao_basis_ket)))

        call init_basis_transformations(ao_basis_bra%basis_angular,trans_bra)
        call init_basis_transformations(ao_basis_ket%basis_angular,trans_ket)

   end subroutine init_overlap_transformations

   subroutine init_basis_transformations (basis_angular,trans)

        implicit none

        integer,intent(in)                 :: basis_angular
        type(sph_transform_t),intent(inout):: trans(:)
        integer                            :: lval, ncart, nsph

        select case (basis_angular)
        case (1)
           return
        case (2)
           do lval = 1, size(trans)
              ncart = nfunctions(lval,1)
              nsph  = nfunctions(lval,2)
              allocate(trans(lval)%matrix(ncart,nsph))
              call get_sphcom(lval,trans(lval)%matrix)
           end do
        case default
           error stop 'unknown angular type in init_basis_transformations'
        end select

   end subroutine init_basis_transformations

   integer function max_orb_momentum (ao_basis)

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis
        integer                           :: ishell

        max_orb_momentum = 0
        do ishell = 1, ao_basis%nshells
           max_orb_momentum = max(max_orb_momentum,ao_basis%gtos(ishell)%orb_momentum)
        end do

   end function max_orb_momentum

end module get_overlap
