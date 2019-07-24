   program generate_iao

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_utils
        use read_xyz
        use get_overlap
        use make_iao
        use make_ibo
        use utils
        use check
        use energy
        use properties
        use atomic_fragments

        implicit none

        complex(8), parameter              :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        real(8), parameter                 :: bohr_to_angs = 0.5291772083D0
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(cmo)                          :: cspinor_original ! MO-coefficients before localization 
        type(cmo), allocatable             :: cspinor_fragment(:) ! MO-coefficients of the fragments
        type(atom), pointer                :: first_atom,current_atom ! linked list to store coordindate information
        type(cmo)                          :: cspinor_localized ! MO-coefficients after localization 
        type(cmo)                          :: cspinor_localized_ibo ! MO-coefficients after localization (IBOs)
        type(cmo)                          :: cspinor_frag_IAO ! MO-coefficients in IAO basis for the fragments
        type(cmo)                          :: cspinor_frag_IAOval ! valence MO-coefficients in IAO basis for the fragments
        integer                            :: luxyz, ierr, i, j, k, l
        integer                            :: nmo1, nmo2, imo, jmo, nmoi, nmoj, nao1, ncore, nvale
        integer                            :: n_fragments ! the number of fragments used to project on to
        integer                            :: nocc ! Number of occupied spinors
        integer                            :: nvir ! Number of virtual spinors
        integer                            :: charge ! Charge of the molecule
        integer                            :: Expnt ! Exponent to be used in the localization procedure (Make IBO).
        character(11)                      :: version ! Version of the IAO construction (Standard_2013, Simple_2013, Simple_2014)
        complex(8), allocatable            :: S_block(:,:) ! sub block of the overlap matrix
        complex(8), allocatable            :: S11(:,:),S21(:,:),S12(:,:),S22(:,:) !  the overlap matrices
        complex(8), allocatable            :: P12(:,:),P21(:,:)         !  the projection operators
        complex(8), allocatable            :: IAO(:,:)                  !  coefficient matrix, IAOs (MO basis)
        complex(8), allocatable            :: IAO_aoB1(:,:,:)             !  coefficient matrix, IAOs in terms of AO basis B1
        complex(8), allocatable            :: Cocc_IAO(:,:)             !  coefficient matrix, occupied MOs in terms of IAOs.
        complex(8), allocatable            :: Cvir_IAO(:,:)             !  coefficient matrix, virtual MOs in terms of IAOs.
        complex(8), allocatable            :: Cocc_IAO_aoB1(:,:,:)        !  coefficient matrix, occupied MOs in terms of IAOs transformed back to AO basis B1
        complex(8), allocatable            :: Cvir_IAO_aoB1(:,:,:)        !  coefficient matrix, occupied MOs in terms of IAOs transformed back to AO basis B1
        complex(8), allocatable            :: Cocc_IBO(:,:)             !  coefficient matrix, Localized occupied MOs in terms of IAOs. (called IBO)
        complex(8), allocatable            :: Cocc_IBO_aoB1(:,:,:)        !  coefficient matrix, IBOs transformed back to AO basis B1
        complex(8), allocatable            :: Cocc_IBO_moB1(:,:)        !  coefficient matrix, IBOs transformed back to MO basis B1
        complex(8), allocatable            :: IAO_frag(:,:,:)           !  coefficient matrix, IAOs of a fragment in terms of AO basis B1.
        complex(8), allocatable            :: coeff1_occ(:,:)
        complex(8), allocatable            :: coeff1_vir(:,:)
        real(8),    allocatable            :: Cocc_IBO_real(:,:)
        character(len=14)                  :: fragfile_dfcoef, fragfile_dfpcmo
        character(len=18)                  :: fragfile_val_dfcoef, fragfile_val_dfpcmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ Input Files         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_free_fileunit(luxyz)
        open (luxyz,file="INPUT",form='formatted')
        read(luxyz,'(A11)') version
        read(luxyz,'(I2)') charge
        read(luxyz,'(I2)') Expnt
        close(luxyz)

        call read_from_dirac (ao_basis_1,cspinor_original,'MOLECULE.dfcoef',ierr)
        nmo1 = cspinor_original%nmo
        nao1 = cspinor_original%nao

        call get_free_fileunit(luxyz)
        open (luxyz,file='MOLECULE.XYZ',form='formatted')

        call get_atoms(luxyz,first_atom)
        close (luxyz,status='keep')
        current_atom => first_atom
        print*,"Molecular structure used to make reference basis"
        ! For now we will just consider atomic fragments, which should be available from a library of DFCOEFs
        n_fragments = 0
        nocc = 0
        do while(associated(current_atom))
           n_fragments = n_fragments + 1
           print*,current_atom%name,current_atom%bascharge,current_atom%x
           current_atom%x = current_atom%x/bohr_to_angs
           nocc = nocc + current_atom%bascharge
           current_atom=>current_atom%next
        end do
        nocc = nocc - charge
        allocate (ao_basis_2(n_fragments))
        allocate (cspinor_fragment(n_fragments))
        nmo2 = 0
        current_atom => first_atom
        do i = 1, n_fragments
           call get_atomic_fragment (nint(current_atom%bascharge),current_atom%x,ao_basis_2(i),cspinor_fragment(i))
           nmo2 = nmo2 + cspinor_fragment(i)%nmo
           current_atom=>current_atom%next
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !-----Compute S11----- (nmo1 x nmo1) (It is equal to Identity)
        allocate (S11(nmo1,nmo1))
        call get_mo_overlap (ao_basis_1,cspinor_original,ao_basis_1,cspinor_original,S11)

        !-----Compute S22----- (nmo2 x nmo2)
        allocate (S22(nmo2,nmo2))
        imo = 0
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           jmo = imo
           do j = i, n_fragments
              nmoj = cspinor_fragment(j)%nmo
              allocate (S_block(nmoi,nmoj))
              call get_mo_overlap (ao_basis_2(i),cspinor_fragment(i),ao_basis_2(j),cspinor_fragment(j),S_block)
                 S22(imo+1:imo+nmoi,jmo+1:jmo+nmoj) = S_block
              if (i /= j) then
                 S22(jmo+1:jmo+nmoj,imo+1:imo+nmoi) = conjg(transpose(S_block))
              end if
              jmo = jmo + nmoj
              deallocate (S_block)
           enddo
           imo = imo + nmoi
        enddo

        !-----Compute S21----- (nmo2 x nmo1)
        allocate (S21(nmo2,nmo1))
        imo = 0
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           allocate (S_block(nmoi,nmo1))
           call get_mo_overlap (ao_basis_2(i),cspinor_fragment(i),ao_basis_1,cspinor_original,S_block)
           S21(imo+1:imo+nmoi,:) = S_block
           imo = imo + nmoi
           deallocate (S_block)
        enddo

        !-----Compute S12----- (nmo1 x nmo2)
        allocate (S12(nmo1,nmo2))
        S12 = conjg(transpose(S21))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Calculate the intrinsic atomic orbital (IAO) basis                                         !
!         Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Make the projector, zposv solves AX = B in which B is input and overwritten by X on output.
        allocate (P12(nmo1,nmo2))
        P12 = S12
        ! P12 = inv(S11) S12 ---> solve S11 P12 = S12
        call zposv( 'U', nmo1, nmo2, S11, nmo1, P12, nmo1, ierr )

        ! In contrast to Knizia's approach, we work in the MO basis.
        ! Note that the MO fragment orbitals are of course localized on the fragments, as they are calculated separately
        ! for each isolated fragment. For now, we consider atomic fragments only (i.e. one fragment = one atom)

        !-----Compute COcc----- (nmo1 x nocc) (diagonal matrix)
        nvir = nmo1-nocc
        allocate(coeff1_occ(nmo1,nocc))
        allocate(coeff1_vir(nmo1,nvir))
        coeff1_occ = ZERO
        forall (i=1:nocc) coeff1_occ(i,i) = ONE
        coeff1_vir = ZERO
        forall (i=1:nvir) coeff1_vir(i+nocc,i) = ONE

        ! Calculate the Intrinsic Atomic Orbitals
        ! expressed in terms of the MO basis B1
        ! -----Compute IAOs----- (nmo1 x nmo2)
        allocate (IAO(nmo1,nmo2))
        select case (version)
           case ('Stndrd_2013')
             write(*,*)" Making IAO according to Standard 2013"
             call MakeIAO_Standard_2013(IAO,coeff1_occ,S11,S21,S22,P12,nmo1,nmo2,nocc)
           case ('Simple_2013')
             write(*,*)" Making IAO according to Simple 2013"
             call MakeIAO_Simple_2013(IAO,coeff1_occ,S11,S12,S21,S22,P12,nmo1,nmo2,nocc)
           case ('Simple_2014')
             write(*,*)" Making IAO according to Simple 2014"
             call MakeIAO_Simple_2014(IAO,coeff1_occ,S21,S22,P12,nmo1,nmo2,nocc)
           case default
             print*, 'No correct argument given to the executable. Default : Simple_2014'
             call MakeIAO_Simple_2014(IAO,coeff1_occ,S21,S22,P12,nmo1,nmo2,nocc)
        end select

        !-----Symmetrically orthonormalize IAOs with S11-----
        call SymOrth(IAO,S11,nmo1,nmo2)

        !-----Express IAO in terms of the original AO basis B1----- (nao1 x nmo2)
        allocate(IAO_aoB1(nao1,nmo2,2))
        IAO_aoB1(:,:,1) = matmul(cspinor_original%coeff(:,:,1),IAO)
        IAO_aoB1(:,:,2) = matmul(cspinor_original%coeff(:,:,2),IAO)

        ! Separate the IAO DFCOEF file for each fragment. 
        ! For instance, it will be interesting to get the IAOs of the metallic atom in a transition metal complex. 
        ! This can be used in the AVAS theory (https://arxiv.org/abs/1701.07862)
        imo=0
        current_atom => first_atom
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           allocate(IAO_frag(nao1,nmoi,2))
           IAO_frag = IAO_aoB1(:,imo+1:imo+nmoi,:)
           write (fragfile_dfcoef,'(i3.3,a11)') nint(current_atom%bascharge),'_iao.dfcoef'
           call alloc_mo(cspinor_frag_IAO,nao1,nmoi)
           cspinor_frag_IAO%coeff = IAO_frag
           cspinor_frag_IAO%total_energy = cspinor_original%total_energy
           call write_mos_to_dirac (ao_basis_1,cspinor_frag_IAO,fragfile_dfcoef)
           call dealloc_mo(cspinor_frag_IAO)
           deallocate(IAO_frag)

           ! Do it also for valence orbitals only
           call atom_minbas_size (nint(current_atom%bascharge),ncore,nvale)
           allocate(IAO_frag(nao1,nvale,2))
           IAO_frag = IAO_aoB1(:,imo+1+ncore:imo+nmoi,:)
           write (fragfile_val_dfcoef,'(i3.3,a15)') nint(current_atom%bascharge),'_iao_val.dfcoef'
           call alloc_mo(cspinor_frag_IAOval,nao1,nvale)
           cspinor_frag_IAOval%coeff = IAO_frag
           cspinor_frag_IAOval%total_energy = cspinor_original%total_energy
           call write_mos_to_dirac (ao_basis_1,cspinor_frag_IAOval,fragfile_val_dfcoef)
           call dealloc_mo(cspinor_frag_IAOval)
           deallocate(IAO_frag)
           imo = imo + nmoi
           current_atom=>current_atom%next
        end do

        !-----Occupied MO in terms of the orthonormal IAO basis----- (nmo2 x nocc)
        allocate (Cocc_IAO(nmo2,nocc))
        Cocc_IAO = matmul(conjg(transpose(IAO)),matmul(S11,coeff1_occ))

        !-----Virtual MO in terms of the orthonormal IAO basis----- (nmo2 x nvir)
        allocate (Cvir_IAO(nmo2,nocc))
        Cvir_IAO = matmul(conjg(transpose(IAO)),matmul(S11,coeff1_vir))

        !-----Transform occupied MO in terms of IAOs back to AO basis B1----- (nao1 x nocc)
        allocate(Cocc_IAO_aoB1(nao1,nocc,2))
        Cocc_IAO_aoB1(:,:,1) = matmul(IAO_aoB1(:,:,1),Cocc_IAO)
        Cocc_IAO_aoB1(:,:,2) = matmul(IAO_aoB1(:,:,2),Cocc_IAO)

        !-----Transform virtual MO in terms of IAOs back to AO basis B1----- (nao1 x nocc)
        allocate(Cvir_IAO_aoB1(nao1,nvir,2))
        Cvir_IAO_aoB1(:,:,1) = matmul(IAO_aoB1(:,:,1),Cvir_IAO)
        Cvir_IAO_aoB1(:,:,2) = matmul(IAO_aoB1(:,:,2),Cvir_IAO)

        ! Write the coefficient matrix into the unformatted dfcoef file.
        call alloc_mo(cspinor_localized,nao1,nmo1)
        cspinor_localized%coeff(:,1:nocc,:) = Cocc_IAO_aoB1
        ! cspinor_localized%coeff(:,nocc+1:nmo1,:) = Cvir_IAO_aoB1
        ! keep the original virtual orbitals:
        cspinor_localized%coeff(:,nocc+1:nmo1,:) = cspinor_original%coeff(:,nocc+1:nmo1,:)
        cspinor_localized%energy(nocc+1:nmo1) = cspinor_original%energy(nocc+1:nmo1)
        cspinor_localized%total_energy = cspinor_original%total_energy
        call write_mos_to_dirac (ao_basis_1,cspinor_localized,"iao.dfcoef")
        call write_mos_to_dirac_formatted (cspinor_localized,"iao.dfpcmo")

        ! Compute partial charges with occupied MO expressed in terms of IAOs.
        call partial_charges(Cocc_IAO,cspinor_fragment,n_fragments,nmo2,nocc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Calculate the intrinsic bonding orbital (IBO) basis                                        !
!         Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !-----Compute IBOs (Localized occupied MO in terms of IAOs)----- (nmo2 x nocc)
        allocate(Cocc_IBO_real(nmo2,nocc)) ! should be replaced by Cocc_IBO once it will work for complex orbitals.
        allocate(Cocc_IBO(nmo2,nocc))
        call LocalizeIbosInIaoBasis(real(Cocc_IAO),Cocc_IBO_real,cspinor_fragment,n_fragments,nmo2,nocc,Expnt)

        ! Temporary rewrite into complex format... :
        do i=1,nmo2
         do j=1,nocc
          Cocc_IBO(i,j) = DCMPLX(Cocc_IBO_real(i,j),0.D0)
         enddo
        enddo

        !-----Transform IBOs back to AO basis B1----- (nao1 x nocc)
        allocate(Cocc_IBO_aoB1(nao1,nocc,2))
        Cocc_IBO_aoB1(:,:,1) = matmul(IAO_aoB1(:,:,1),Cocc_IBO)
        Cocc_IBO_aoB1(:,:,2) = matmul(IAO_aoB1(:,:,2),Cocc_IBO)
        ! Write the coefficient matrix into the unformatted dfcoef file.
        call alloc_mo(cspinor_localized_ibo,nao1,nmo1)
        cspinor_localized_ibo%coeff(:,1:nocc,:) = Cocc_IBO_aoB1
        ! cspinor_localized_ibo%coeff(:,1+nocc:nmo1,:) = Cvir_IAO_aoB1
        ! keep the original virtual orbitals:
        cspinor_localized_ibo%coeff(:,nocc+1:nmo1,:) = cspinor_original%coeff(:,nocc+1:nmo1,:)
        cspinor_localized_ibo%energy(nocc+1:nmo1) = cspinor_original%energy(nocc+1:nmo1)
        cspinor_localized_ibo%total_energy = cspinor_original%total_energy
        call write_mos_to_dirac (ao_basis_1,cspinor_localized_ibo,"ibo.dfcoef")
        call write_mos_to_dirac_formatted (cspinor_localized_ibo,"ibo.dfpcmo")

        !-----Transform IBOs back to MO basis B1----- (nmo1 x nocc)
        ! It will be used to test the HF energy (using the Fock and H_core matrix of basis B1)
        allocate(Cocc_IBO_moB1(nmo1,nocc))
        Cocc_IBO_moB1 = matmul(IAO,Cocc_IBO)

        ! In Knizia's program, another call is made to sort the new localized orbitals by energy.
        ! I think in our case they are already sorted, but I'm not 100% convinced...

        ! Let's try to compute the partial charges using IBO?
        call partial_charges(Cocc_IBO,cspinor_fragment,n_fragments,nmo2,nocc)

!!!!!!!!!!!!!!!!!!!!
!      TEST        !
!!!!!!!!!!!!!!!!!!!!
! This test requires the presence of the MRCONEE file to get the core matrix h. Also, ACTIVE .all has to be
! set in the Dirac input to be sure that the latter matrix has the correct dimension nmo1 x nmo1.

!        call test_IAO_and_IBO(IAO,Cocc_IAO,Cocc_IBO_moB1,S11,cspinor_original,nmo1,nmo2,nocc)

!!!!!!!!!!!!!!!!!!!!
!       END        !
!!!!!!!!!!!!!!!!!!!!
!       clean up list of atoms
        call del_atoms(first_atom)
!     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') ) ! Useful format to write complex matrices.

   end program generate_iao
