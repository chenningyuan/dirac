module avas

contains

  subroutine get_iao_avas (ao_basis_1,spinor_B1,spinor_iao,nmo_frag,nocc,spinor_avas)

  		use rose_ao
		use rose_mo
		use rose_input_data

		implicit none

        type(basis_set_info_t),intent(in)  :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(in)          :: spinor_B1  ! Basis B1, canonical MO-coefficients
        type(mo_basis),intent(in)          :: spinor_iao ! Complete set of intrinsic orbitals, sorted by fragment
        type(mo_basis),intent(out)         :: spinor_avas ! AVAS spinors, can be used as starting spinors for active space calculations
        integer, intent(in)                :: nmo_frag(n_fragments)
        integer,intent(in)                 :: nocc

        type(mo_basis)                     :: spinor_B2  ! Selected set of intrinsic orbitals
        integer,allocatable                :: iao_list(:)
        integer                            :: i,j,imo,jmo,nmoi

        ! count orbitals to be put in the active space
        jmo = 0
        do i = 1, n_fragments
           jmo = jmo + avas_orbitals%column_length(i)
        enddo

        ! make list of the indices of these IAOs in the full space
        allocate(iao_list(jmo))
        jmo = 0
        nmoi = 0
        do j = 1, n_fragments
           do i = 1, avas_orbitals%column_length(j)
              jmo = jmo + 1
              call JAGetElement(avas_orbitals,i,j,imo)
              iao_list(jmo) = nmoi+imo
           enddo
           nmoi = nmoi + nmo_frag(j)
        end do

        ! make spinor object with only the selected IAOs
        call copy_mo_basis (spinor_iao,spinor_B2)
        call reduce_mo_basis (spinor_B2,jmo,mo_list=iao_list)

        ! execute the avas algorithm
        call make_avas(ao_basis_1,ao_basis_1,spinor_B1,spinor_B2,nocc,avas_threshold,spinor_avas)

   end subroutine get_iao_avas

   subroutine make_avas(ao_basis_1,ao_basis_2,spinor_B1,spinor_B2,nocc,threshold,spinor_avas)
   ! Ref: E.R. Sayfutyarova, Q. Sun, G.K.-L. Chan, G. Knizia, Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals, J Chem Theory Comput. 13 (2017) 4063–4078. doi:10.1021/acs.jctc.7b00128.
   ! Basis B2 can be the IAOs, and one wants to compute the overlap between this "minimal" basis and a larger basis B1. Especially, one could take only a target set of AOs in B2, like the one corresponding to the d orbitals of a transition metal. Then, a complete active space is automatically constructed according to the calculated overlap matrices for occupied and virtual orbitals.

        use rose_ao
        use rose_mo
        use rose_utils
        use read_xyz
        use linear_algebra
        use get_overlap
        use rose_input_data

        implicit none
        type(basis_set_info_t), intent(in) :: ao_basis_1 ! the ao basis B1
        type(basis_set_info_t), intent(in) :: ao_basis_2 ! the ao basis B2
        type(mo_basis), intent(out)        :: spinor_avas ! the AVAS spinors: inactive occupied, active occupied, active virtual, inactive virtual, all expressed in B1
        type(mo_basis),intent(in)          :: spinor_B1  ! MO-coefficients of basis B1
        type(mo_basis),intent(in)          :: spinor_B2  ! MO-coefficients of basis B2 (only valence)
        integer, intent(in)                :: nocc
        real(8), intent(in)                :: threshold

        integer                 :: i,j,k,ierr
        integer                 :: nao1,nmo1,nao2,nmo,nmo2,rcq,nact_occ,nact_vir,nvir
        real(8), allocatable    :: Proj(:,:,:),Socc(:,:,:),Cocc(:,:,:),Svir(:,:,:),Cvir(:,:,:)
        real(8), allocatable    :: Cocc_tilde(:,:,:),Cvir_tilde(:,:,:),Ctilde(:,:,:)
        real(8), allocatable    :: eigvec_occ(:,:,:),eigvec_vir(:,:,:)
        real(8), allocatable    :: eigval_occ(:),eigval_vir(:),eigval(:)
        real(8), allocatable    :: S21(:,:,:),S12(:,:,:),S22(:,:,:)
        character(len=13)       :: fmt_tmp
       
        write(*,*)
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(*,*) "@                                                @"
        write(*,*) "@        Automated Valence Active Space          @"
        write(*,*) "@                    (AVAS)                      @"
        write(*,*) "@                                                @"
        write(*,*) "@               By Bruno Senjean                 @"
        write(*,*) "@                                                @"
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(*,*)

        write(*,'(A,E10.2)') " Threshold: ",threshold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CHECK IF OPTIONS ARE CORRECT         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call define_algebra_type(spatial_orbitals,restricted,rcq)

        nao2 = spinor_B2%nao
        nmo2 = spinor_B2%nmo
        if (mod(nmo2,2) .ne. 0 .and. .not. restricted) stop "**ERROR** Number of spinors to keep should be an even number"

        nao1 = spinor_B1%nao
        nmo1 = spinor_B1%nmo

        if (nmo1 .lt. nmo2) stop "nmo1 shouldn't be inferior to nmo2"

        write(*,*)
        write(*,*) "####################################"
        write(*,*) "#         Basis set size           #"
        write(*,*) "####################################"
        write(*,*)
        write(*,'(A,I6)') " No. basis functions    B1: ",nao1
        write(*,'(A,I6)') " No. molecular orbitals B1: ",nmo1
        write(*,'(A,I6)') " No. occupied  orbitals B1: ",nocc
        write(*,'(A,I6)') " No. basis functions    B2: ",nao1
        write(*,'(A,I6)') " No. molecular orbitals B2: ",nmo2
        write(*,*)      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(S21(nmo2,nmo1,rcq))
        allocate(S12(nmo1,nmo2,rcq))
        allocate(S22(nmo2,nmo2,rcq))
        call get_mo_overlap (ao_basis_2,spinor_B2,ao_basis_2,spinor_B2,S22)
        call get_mo_overlap (ao_basis_2,spinor_B2,ao_basis_1,spinor_B1,S21)
        call transpose_conjg(S21,S12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Projection Matrices    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! B2 form an orthonormalized MO basis, so sigma^(-1) in Eq.13 of the article is the identity, so not considered here.
        allocate(Proj(nmo1,nmo1,rcq))
        call matrix_mul(S12,S21,Proj)

        ! Compute the projected overlap matrices Socc = C†occ P Cocc and Svir = C†vir P Cvir
        ! In the MO basis, Cocc and Cvir are diagonal matrices.
        nvir = nmo1 - nocc
        allocate(Socc(nocc,nocc,rcq))
        allocate(Svir(nvir,nvir,rcq))
        Socc = Proj(1:nocc,1:nocc,:)
        Svir = Proj(nocc+1:nocc+nvir,nocc+1:nocc+nvir,:)

        ! deallocate matrices not used anymore:
        deallocate(S12,S21,S22,Proj)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DIAGONALIZE PROJECTED OVERLAP MATRICES TO DEFINE THE AVAS     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(eigvec_occ(nocc,nocc,rcq))
        allocate(eigval_occ(nocc))
        allocate(eigvec_vir(nvir,nvir,rcq))
        allocate(eigval_vir(nvir))
        eigvec_occ = Socc
        eigvec_vir = Svir
        call diag_matrix(eigvec_occ,eigval_occ)
        call diag_matrix(eigvec_vir,eigval_vir)

        ! Count the number of (spatial) active orbitals to be considered.
        ! There will be at most 2*nmo2 of them.
        nact_occ = 0
        nact_vir = 0
        do i=1,nocc
         if (eigval_occ(i) > threshold) nact_occ = nact_occ + 1
        enddo
        do i=1,nvir
         if (eigval_vir(i) > threshold) nact_vir = nact_vir + 1
        enddo
        if ((nact_occ+nact_vir) .gt. 2*nmo2) print*,"**WARNING**: nact is larger than 2*nmo2 which shouldn't be possible."

        ! Now we sort the coefficients to have frozen-occupied / actives / unoccupied-virtuals:
        ! Eigenvectors are sorted in an ascending order (lowest eigenvalues first) such that occupied orbitals
        ! are already sorted, but the virtuals need to be in reverse order.
        allocate(Ctilde(nmo1,nmo1,rcq))
        allocate(eigval(nmo1))
        Ctilde = 0.0
        Ctilde(1:nocc,1:nocc,:) = eigvec_occ
        eigval(1:nocc)          = eigval_occ
        do j=1,nvir
          Ctilde(nocc+1:nocc+nvir,nocc+j,:) = eigvec_vir(:,nvir-j+1,:)
          eigval(nocc+j)                    = eigval_vir(nvir-j+1)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       WRITE INFORMATION ABOUT AVAS IN THE OUTPUT     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(*,*)
        write(*,*) "####################################"
        write(*,*) "#    Active space construction     #"
        write(*,*) "####################################"
        write(*,*)
        write(*,'(A,I6)') " No. of frozen occupied: ",nocc-nact_occ
        write(*,'(A,I6)') " No. of active occupied: ",nact_occ
        write(*,'(A,I6)') " No. of active virtuals: ",nact_vir
        write(*,'(A,I6)') " No. of frozen virtuals: ",nvir-nact_vir
        write(*,*)
        write(*,*) "--------- Eigenvalues of the Projected overlap matrices (active MOs) ---------"
        write(*,*)
        write(*,'(A4,A1,A24)') "MO"," ","eigval"
        do i=nocc-nact_occ+1,nocc+nact_vir
          write(*,'(I4,A1,F24.16)') i," ",eigval(i)
        enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       WRITE COEFFICIENT MATRICES THAT WILL BE USED IN POST-HF CALCULATION     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call IMO_to_mobasis(Ctilde,spinor_B1,spinor_avas)

        write(*,*) 
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        write(*,*) "@                                                @"
        write(*,*) "@        Automated Valence Active Space          @"
        write(*,*) "@                    (AVAS)                      @"
        write(*,*) "@                                                @"
        write(*,*) "@                     END                        @"
        write(*,*) "@                                                @"
        write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

   end subroutine make_avas

end module avas

