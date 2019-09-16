program avas
   ! Ref: E.R. Sayfutyarova, Q. Sun, G.K.-L. Chan, G. Knizia, Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals, J Chem Theory Comput. 13 (2017) 4063–4078. doi:10.1021/acs.jctc.7b00128.
   ! Basis B2 can be the IAOs, and one wants to compute the overlap between this "minimal" basis and a larger basis B1. Especially, one could take only a target set of AOs in B2, like the one corresponding to the d orbitals of a transition metal. Then, a complete active space is automatically constructed according to the calculated overlap matrices for occupied and virtual orbitals.
   ! For the moment, I did just put the equations here where B2 is the set of orthonormalized IAOs and B1 is the same as in the IAO construction...

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_utils
        use get_overlap
        use utils
        use atomic_fragments

        implicit none
        complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        integer                 :: i,j,ierr,luxyz
        type(basis_set_info_t)  :: ao_basis_1 ! the ao basis B1
        type(basis_set_info_t)  :: ao_basis_2 ! the ao basis B2
        type(cmo)               :: cspinor_B1 ! MO-coefficients of basis B1
        type(cmo)               :: cspinor_atom ! MO-coefficients of basis B2 (all orbitals)
        type(cmo)               :: cspinor_B2 ! MO-coefficients of basis B2 (valence only !)
        complex(8), allocatable :: Proj(:,:), Socc(:,:), Cocc(:,:), Svir(:,:), Cvir(:,:)
        complex(8), allocatable :: Cocc_tilde(:,:), Cvir_tilde(:,:), Ccas(:,:), Ctilde(:,:), Ctilde_ao(:,:,:)
        complex(8), allocatable :: eigvec_occ(:,:), eigvec_vir(:,:)
        real(8),    allocatable :: eigval_occ(:), eigval_vir(:)
        complex(8), allocatable :: S21(:,:),S12(:,:), S22(:,:)
        integer                 :: nao1,nmo1,nao2,nmo2
        integer                 :: nelec,nocc,nvir,atom_number,charge,ncore,nvale
        character(len=10)       :: B2file
        complex(8)              :: sum_column
       
        print*, "******************* START AVAS *******************"

        call get_free_fileunit(luxyz)
        open (luxyz,file="INPUT_AVAS",form='formatted')
        read(luxyz,'(I3)') nelec
        read(luxyz,'(I2)') charge
        read(luxyz,'(I3)') atom_number
        close(luxyz)
        nocc = nelec - charge

        call read_from_dirac (ao_basis_1,cspinor_B1,'B1.dfcoef',ierr)
        nmo1 = cspinor_B1%nmo
        nao1 = cspinor_B1%nao

        call read_from_dirac (ao_basis_2,cspinor_atom,'atom.dfcoef',ierr)
        nmo2 = cspinor_atom%nmo
        nao2 = cspinor_atom%nao

        call atom_minbas_size (atom_number,ncore,nvale)
        call alloc_mo (cspinor_B2,nao2,nvale)
        cspinor_B2%energy = cspinor_atom%energy(ncore+1:ncore+nvale)
        cspinor_B2%coeff = cspinor_atom%coeff(:,ncore+1:ncore+nvale,:)
        cspinor_B2%index  = cspinor_atom%index(ncore+1:ncore+nvale)
        call write_dffiles("B2.dfcoef","B2.dfpcmo",ao_basis_2,cspinor_B2%coeff,nao2,nvale,"atom.dfcoef")

        ! The number of orbitals in B2 is now the number of valence orbitals
        nmo2 = nvale

        ! we don't need this coefficient set anymore, free up the memory
        call dealloc_mo(cspinor_atom)

        ! Compute overlap matrices:
        allocate(S21(nmo2,nmo1))
        allocate(S12(nmo1,nmo2))
        allocate(S22(nmo2,nmo2))
        call get_mo_overlap (ao_basis_2,cspinor_B2,ao_basis_2,cspinor_B2,S22)
!        print*,
!        print*,"S22 (should be identity)"
!        do i=1,nmo2
!         write(*,'(100G10.2)') (REAL(S22(i,j)),j=1,nmo2)
!        enddo

        call get_mo_overlap (ao_basis_2,cspinor_B2,ao_basis_1,cspinor_B1,S21)
        S12 = conjg(transpose(S21))

        ! Compute the projector. IAO or MO are orthonormalized, so sigma^(-1) in Eq.13 of the article is equal to the identity in the MO basis, that I considered here.
        allocate(Proj(nmo1,nmo1))
        Proj = matmul(S12,S21)

        ! Compute the projected overlap matrices Socc = C†occ P Cocc and Svir = C†vir P Cvir
        ! In the MO basis, Cocc and Cvir are diagonal matrices...
        nvir = nmo1 - nocc
        allocate(Cocc(nmo1,nocc))
        allocate(Cvir(nmo1,nvir))
        Cocc = ZERO
        forall (i=1:nocc) Cocc(i,i) = ONE
        Cvir = ZERO
        forall (i=1:nvir) Cvir(i+nocc,i) = ONE
!        forall (i=1:nvir) Cvir(i,i) = ONE

        allocate(Socc(nocc,nocc))
        Socc = matmul(conjg(transpose(Cocc)),matmul(Proj,Cocc))
        allocate(Svir(nvir,nvir))
        Svir = matmul(conjg(transpose(Cvir)),matmul(Proj,Cvir))

        allocate(eigvec_occ(nocc,nocc))
        allocate(eigval_occ(nocc))
        allocate(eigvec_vir(nvir,nvir))
        allocate(eigval_vir(nvir))

        eigvec_occ = Socc
        eigvec_vir = Svir
        call diag_matrix(eigvec_occ,nocc,eigval_occ)
        call diag_matrix(eigvec_vir,nvir,eigval_vir)

        print*,"Eigval occ:", eigval_occ 
        print*,"Eigval vir:", eigval_vir 

        allocate(Cocc_tilde(nmo1,nocc))
        allocate(Cvir_tilde(nmo1,nvir))
        ! Fortran is column
        Cocc_tilde = matmul(Cocc,eigvec_occ)
        Cvir_tilde = matmul(Cvir,eigvec_vir)

        print*,
        print*,"Check normalization Cocctilde"
        do i=1,nocc
         sum_column = ZERO
         do j=1,nmo1
           sum_column = sum_column + conjg(Cocc_tilde(j,i))*Cocc_tilde(j,i)
         end do
         write(*,*) sum_column
        end do

        print*,
        print*,"Check normalization Cvirtilde"
        do i=1,nvir
         sum_column = ZERO
         do j=1,nmo1
           sum_column = sum_column + conjg(Cvir_tilde(j,i))*Cvir_tilde(j,i)
         end do
         write(*,*) sum_column
        end do

        allocate(Ctilde(nmo1,nmo1))
        do i=1,nmo1
          do j=1,nocc
            Ctilde(i,j) = Cocc_tilde(i,j)
          enddo
          do j=1,nvir
            Ctilde(i,j+nocc) = Cvir_tilde(i,j)
          enddo
        enddo

        allocate(Ctilde_ao(nao1,nmo1,2))
        Ctilde_ao(:,:,1) = matmul(cspinor_B1%coeff(:,:,1),Ctilde)
        Ctilde_ao(:,:,2) = matmul(cspinor_B1%coeff(:,:,2),Ctilde)
        call write_dffiles("avas.dfcoef","avas.dfpcmo",ao_basis_1,Ctilde_ao,nao1,nmo1,"B1.dfcoef")

!The eigenvalues of the projected overlap matrices reflect the degree to which the transformed orbitals overlap with the space of our target AOs.
!Let us include the transformed orbitals that have non-zero overlap with the target AOs in the CAS:
        allocate(Ccas(nmo1,2*nmo2))
        j = 1
        do i=1,nocc
          if (eigval_occ(i) > 1.D-12) then
            Ccas(:,j) = Cocc_tilde(:,i)
            print*,i,j,eigval_occ(i)
            j=j+1
          endif
        enddo
        j = nmo2+1
        do i=1,nvir
          if (eigval_vir(i) > 1.D-12) then
            Ccas(:,j) = Cvir_tilde(:,i)
            print*,i,j,eigval_vir(i)
            j=j+1
          endif
        enddo
        print*,
        print*,"Ccas"
        do i=1,nmo1
           write(*,'(100G10.2)') (REAL(Ccas(i,j)),j=1,2*nmo2)
        enddo

        print*,                            
        print*,"Check normalization"
        do i=1,2*nmo2
         sum_column = ZERO                 
         do j=1,nmo1
           sum_column = sum_column + conjg(Ccas(j,i))*Ccas(j,i)
         end do
         write(*,*) sum_column
        end do

        print*, "******************** END AVAS ********************"
end program avas
