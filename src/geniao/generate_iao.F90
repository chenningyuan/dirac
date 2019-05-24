   program generate_iao

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_utils
        use get_overlap
        use make_iao
        use atomic_fragments
        use utils
        use read_xyzfile_only

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)

        type(basis_set_info_t) :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(cmo):: cspinor_original ! MO-coefficients before localization 
        type(cmo):: cspinor_local    ! MO-coefficients after localization

        integer :: n_fragments ! the number of fragments used to project on to
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(cmo), allocatable :: cspinor_fragment(:) ! MO-coefficients of the fragments

        type(atom),pointer :: first_atom,current_atom ! linked list to store coordindate information

        complex(8), allocatable :: S_block(:,:) ! sub block of the overlap matrix
        complex(8), allocatable :: S11(:,:),S21(:,:),S12(:,:),S22(:,:) !  the overlap matrices
        complex(8), allocatable :: P12(:,:),P21(:,:)          !  the projection operators
        complex(8), allocatable :: c_Ib(:,:)                  !  coefficient matrix, IAOs in terms of MOs of basis 1
        complex(8), allocatable :: c_Ib_ao(:,:,:)             !  coefficient matrix, IAOs in terms of AOs of basis 1
        complex(8), allocatable :: c_Ib_occ(:,:)              !  coefficient matrix, occupied MOs in terms of IAOs.
        complex(8), allocatable :: coeff1_occ(:,:)
        real(8), allocatable :: coeff_quaternion(:,:,:)

        integer :: nocc = 10 ! Should be read from a file, but currently this info is not available
        integer luxyz, ierr, i, j, k, l, nmo1, nmo2, imo, jmo, nmoi, nmoj, nao1

        complex(8), allocatable :: DMAT(:,:)
        complex(8) :: sum_diag
        character(13) :: version


        print*,"Arrived in generate_iao"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ Input Files         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call get_free_fileunit(luxyz)
        open (luxyz,file='MOLECULE.XYZ',form='formatted')
        call get_atoms(luxyz,first_atom)
        close (luxyz,status='keep')
        current_atom => first_atom
        print*,"Molecular structure used to make reference basis"
        ! For now we will just consider atomic fragments, which should be available from a library of DFCOEFs
        n_fragments = 0
        do while(associated(current_atom))
           n_fragments = n_fragments + 1
           print*,current_atom%name,current_atom%bascharge,current_atom%x
           current_atom=>current_atom%next
        end do
        allocate (ao_basis_2(n_fragments))
        allocate (cspinor_fragment(n_fragments))
        nmo2 = 0
        current_atom => first_atom
        do i = 1, n_fragments
           call get_atomic_fragment (nint(current_atom%bascharge),current_atom%x,ao_basis_2(i),cspinor_fragment(i))
           nmo2 = nmo2 + cspinor_fragment(i)%nmo
           current_atom=>current_atom%next
        end do

        call read_from_dirac (ao_basis_1,cspinor_original,'MOLECULE.dfcoef',ierr)

        nmo1 = cspinor_original%nmo
        nao1 = cspinor_original%nao
        print*,"Read MO coefficients: nmo1,nmo2(size of basis 1 and basis 2)",nmo1,nmo2

        ! The following Coeff file is here to check that I indeed get the same format as the DFPCMO_mol file.
        ! It will be removed when I'll be sure it is indeed correct in any cases.
        allocate(coeff_quaternion(cspinor_original%nao,int(cspinor_original%nmo/2),4))
        call convert_complex_to_quaternion_format(cspinor_original%coeff,coeff_quaternion,cspinor_original%nao,cspinor_original%nmo)
        call write_dffiles("DFCOEFtest","DFPCMOtest",cspinor_original%coeff,nao1,nmo1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Compute Overlap Matrices     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print*,
        print*,'-----Compute S11----- (',nmo1,'x',nmo1,') (It is equal to Identity)'
        allocate (S11(nmo1,nmo1))
        call get_mo_overlap (ao_basis_1,cspinor_original,ao_basis_1,cspinor_original,S11)

        print*
        print*,'-----Compute S22-----(',nmo2,'x',nmo2,')'
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
        do i=1,nmo2
           write(*,'(100G10.2)')(REAL(S22(i,j)),j=1,nmo2)
        enddo

        print*,
        print*,'-----Compute S21-----(',nmo2,'x',nmo1,') (Inconvenient to print)'
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

        allocate (S12(nmo1,nmo2))
        print*,
        print*,'-----Compute S12-----(',nmo1,'x',nmo2,')'
        S12 = conjg(transpose(S21))
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(S12(i,j)),j=1,nmo2)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Calculate the intrinsic atomic orbital (IAO) basis                                         !
!         Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Make the projector, zposv solves AX = B in which B is input and overwritten by X on output.
        allocate (P12(nmo1,nmo2))
        print*,
        print*,'-----Compute P12-----(',nmo1,'x',nmo2,')'
        P12 = S12
        call zposv( 'U', nmo1, nmo2, S11, nmo1, P12, nmo1, ierr )
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(P12(i,j)),j=1,nmo2)
        enddo

        ! In contrast to Knizia's approach, we work in the MO basis.
        ! Note that the MO fragment orbitals are of course localized on the fragments, as they are calculated separately
        ! for each isolated fragment. For now, we consider atomic fragments only (i.e. one fragment = one atom)
        allocate(coeff1_occ(nmo1,nocc))
        print*,
        print*,'-----Compute COcc-----(',nmo1,'x',nocc,') (diagonal matrix)'
        coeff1_occ = ZERO
        forall (i=1:nocc) coeff1_occ(i,i) = ONE

        ! Calculate the coefficient matrix to go from the MOs of the large basis 1 to the
        ! new set of MOs (basis 2 which combines the MINIMAL basis (core and valence orbs) of the fragments)
        print*,
        print*,'-----Compute CIb-----(',nmo1,'x',nmo2,')'
        allocate (c_Ib(nmo1,nmo2))
        call get_command_argument(1, version)
        select case (version)
           case ('Standard_2013')
             write(*,*)" Making IAO according to Standard 2013"
             call MakeIAO_Standard_2013(c_Ib,coeff1_occ,S11,S21,S22,P12,nmo1,nmo2,nocc)
           case ('Simple_2013')
             write(*,*)" Making IAO according to Simple 2013"
             call MakeIAO_Simple_2013(c_Ib,coeff1_occ,S11,S12,S21,S22,P12,nmo1,nmo2,nocc)
           case ('Simple_2014')
             write(*,*)" Making IAO according to Simple 2014"
             call MakeIAO_Simple_2014(c_Ib,coeff1_occ,S21,S22,P12,nmo1,nmo2,nocc)
           case default
             print*, 'No correct argument given to the executable. Default : Simple_2014'
             call MakeIAO_Simple_2014(c_Ib,coeff1_occ,S21,S22,P12,nmo1,nmo2,nocc)
        end select

!       filter out occupation numbers slightly broken due to rounding... At least it clarifies the output
!        do i=1,nmo1
!           do j=1,nmo2
!              if (REAL(c_Ib(i,j)) > 0.999999999999D0 .and. AIMAG(c_Ib(i,j)) < 1.0D-12) c_Ib(i,j) = ONE
!              if (ABS(REAL(c_Ib(i,j))) < 1.0D-12 .and. ABS(AIMAG(c_Ib(i,j))) < 1.0D-12) c_Ib(i,j) = ZERO
!           end do
!        end do
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(c_Ib(i,j)),j=1,nmo2)
        enddo

        allocate (c_Ib_occ(nmo2,nocc))
        print*,
        print*,'-----Re-express occupied orbitals in terms of IAO basis-----(',nmo2,'x',nocc,')'
        c_Ib_occ = matmul(conjg(transpose(c_Ib)),matmul(S11,coeff1_occ))
        do i=1,nmo2
          write(*,'(100G10.2)') (REAL(c_Ib_occ(i,j)),j=1,nocc)
        enddo

!        print*,
!        print*,"-----Symmetrically orthogonalize intrinsic basis (in the mo basis)-----"
!        call SymOrth(c_Ib,S11,size(S11(:,1)),size(c_Ib(1,:)))

!       filter out occupation numbers slightly broken due to rounding... At least it clarifies the output
!        do i=1,nmo1
!           do j=1,nmo2
!              if (REAL(c_Ib(i,j)) > 0.999999999999D0 .and. AIMAG(c_Ib(i,j)) < 1.0D-12) c_Ib(i,j) = ONE
!              if (ABS(REAL(c_Ib(i,j))) < 1.0D-12 .and. ABS(AIMAG(c_Ib(i,j))) < 1.0D-12) c_Ib(i,j) = ZERO
!           end do
!        end do
!        print*,'Print real part of CIb after orthogonalization (',nmo1,'x',nmo2,')'
!        do i=1,nmo1
!           write(*,'(100G10.2)')(REAL(c_Ib(i,j)),j=1,nmo2)
!        enddo

        ! Transform the coefficient matrix back into the atomic orbital basis 1.
        ! This new coefficient matrix gives the transformation from the AOs of basis 1
        ! to the IAOs.
        allocate(c_Ib_ao(nao1,nmo2,2))
        c_Ib_ao(:,:,1) = matmul(cspinor_original%coeff(:,:,1),c_Ib)
        c_Ib_ao(:,:,2) = matmul(cspinor_original%coeff(:,:,2),c_Ib)
        call write_dffiles("DFCOEF_iao","DFPCMO_iao",c_Ib_ao,nao1,nmo2)

!       clean up list of atoms
        call del_atoms(first_atom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  STUFF not needed for now. Will be removed later maybe.  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        allocate(DMAT(nmo1,nmo1))
!        DMAT = 0.D0
!        do i=1,nmo1
!         do j=1,nmo1
!           do k=1,nocc
!              DMAT(i,j) = DMAT(i,j) + conjg(c_Ib(i,k))*c_Ib(j,k)
!           enddo
!         enddo
!        enddo
!        print*,"PRINT DIAGONAL OF DENSITY MATRIX"
!        write(*,'(100G10.2)') (REAL(DMAT(i,i)),i=1,nmo1)
!        sum_diag = 0.0
!        do i=1,nmo1
!         sum_diag = sum_diag + DMAT(i,i)
!        enddo
!        print*,"Sum of diagonal (should to be equal to nocc = 10):",sum_diag

!        print*,
!        print*,'----***** MAKE IBOs *****-----'
!        call MakeIBO(c_Ib,coeff1_occ,S11,cspinor_fragment,n_fragments,nmo1,nmo2,nocc)

!     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') ) ! Useful format to write complex matrices.

   end program generate_iao
