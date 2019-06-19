   program generate_iao

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_utils
        use read_xyz
        use get_overlap
        use make_iao
        use utils
        use atomic_fragments

        implicit none

        complex(8), parameter              :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        type(basis_set_info_t)             :: ao_basis_1 ! the ao basis used in the generation of the MO coefficients
        type(basis_set_info_t),allocatable :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(cmo)                          :: cspinor_original ! MO-coefficients before localization 
        type(cmo), allocatable             :: cspinor_fragment(:) ! MO-coefficients of the fragments
        type(atom), pointer                :: first_atom,current_atom ! linked list to store coordindate information
        integer                            :: luxyz, ierr, i, j, k, l
        integer                            :: nmo1, nmo2, imo, jmo, nmoi, nmoj, nao1, ncore, nvale
        integer                            :: n_fragments ! the number of fragments used to project on to
        integer                            :: nocc ! Number of occupied spinors
        integer                            :: charge ! Charge of the molecule
        character(11)                      :: version ! Version of the IAO construction (Standard_2013, Simple_2013, Simple_2014)
        complex(8), allocatable            :: S_block(:,:) ! sub block of the overlap matrix
        complex(8), allocatable            :: S11(:,:),S21(:,:),S12(:,:),S22(:,:) !  the overlap matrices
        complex(8), allocatable            :: P12(:,:),P21(:,:)          !  the projection operators
        complex(8), allocatable            :: c_Ib(:,:)                  !  coefficient matrix, IAOs in terms of MOs of basis 1
        complex(8), allocatable            :: c_Ib_ao(:,:,:)             !  coefficient matrix, IAOs in terms of AOs of basis 1
        complex(8), allocatable            :: c_Ib_occ(:,:)              !  coefficient matrix, occupied MOs in terms of IAOs.
        complex(8), allocatable            :: c_Ib_frag(:,:,:)             !  coefficient matrix, IAOs of a fragment in terms of AOs of basis 1.
        complex(8), allocatable            :: coeff1_occ(:,:)
        complex(8), allocatable            :: HF_RDM(:,:), IbRDM(:,:), Id(:,:)
        complex(8)                         :: sum_diag, sum_column
        complex(8), allocatable            :: check(:,:)
        real(8),    allocatable            :: coeff_quaternion(:,:,:)
        character(len=14)                  :: fragfile_dfcoef, fragfile_dfpcmo
        character(len=18)                  :: fragfile_val_dfcoef, fragfile_val_dfpcmo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ Input Files         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_free_fileunit(luxyz)
        open (luxyz,file="INPUT",form='formatted')
        read(luxyz,'(A11)') version
        read(luxyz,'(I2)') charge
        close(luxyz)

        call read_from_dirac (ao_basis_1,cspinor_original,'MOLECULE.dfcoef',ierr)
        nmo1 = cspinor_original%nmo
        nao1 = cspinor_original%nao
        ! The following Coeff file is here to check that I indeed get the same format as the DFPCMO_mol file.
        ! It will be removed when I'll be sure it is indeed correct in any cases.
        allocate(coeff_quaternion(cspinor_original%nao,int(cspinor_original%nmo/2),4))
        call write_dffiles("test.dfcoef","test.dfpcmo",ao_basis_1,cspinor_original%coeff,nao1,nmo1,"MOLECULE.dfcoef")


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
        select case (version)
           case ('Stndrd_2013')
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

        print*,
        print*,"Check normalization (before the symmetrically orthonormalized procedure):"
        do i=1,nmo2
         sum_column = ZERO
         do j=1,nmo1
           sum_column = sum_column + conjg(c_Ib(j,i))*c_Ib(j,i)
         end do
         write(*,*) sum_column
        end do

        allocate(check(nmo2,nmo2))
        check = matmul(conjg(transpose(c_Ib)),c_Ib)
        print*,
        print*,'-----Compute CIb.T*CIb-----(',nmo2,'x',nmo2,')'
        do i=1,nmo2
           write(*,'(100G10.2)')(REAL(check(i,j)),j=1,nmo2)
        enddo
        deallocate(check)

        print*,
        print*,"-----Symmetrically orthonormalize intrinsic basis (in the mo basis)-----"
        call SymOrth(c_Ib,S11,nmo1,nmo2)

!       filter out occupation numbers slightly broken due to rounding... At least it clarifies the output
        do i=1,nmo1
           do j=1,nmo2
              if (REAL(c_Ib(i,j)) > 0.999999999999D0 .and. AIMAG(c_Ib(i,j)) < 1.0D-12) c_Ib(i,j) = ONE
              if (ABS(REAL(c_Ib(i,j))) < 1.0D-12 .and. ABS(AIMAG(c_Ib(i,j))) < 1.0D-12) c_Ib(i,j) = ZERO
           end do
        end do

        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(c_Ib(i,j)),j=1,nmo2)
        enddo

        print*,
        print*,"Check normalization (after the symmetrically orthonormalized procedure):"
        do i=1,nmo2
         sum_column = ZERO
         do j=1,nmo1
           sum_column = sum_column + conjg(c_Ib(j,i))*c_Ib(j,i)
         end do
         write(*,*) sum_column
        end do

        print*,
        print*,'-----Compute CIb.T*CIb-----(',nmo2,'x',nmo2,')'
        allocate(check(nmo2,nmo2))
        check = matmul(conjg(transpose(c_Ib)),c_Ib)
        do i=1,nmo2
           write(*,'(100G10.2)')(REAL(check(i,j)),j=1,nmo2)
        enddo
        deallocate(check)

        allocate (c_Ib_occ(nmo2,nocc))
        print*,
        print*,'-----Re-express occupied orbitals in terms of IAO basis-----(',nmo2,'x',nocc,') &
                (I think this is not needed, only if we want to construct IBOs)'
        c_Ib_occ = matmul(conjg(transpose(c_Ib)),matmul(S11,coeff1_occ))
        do i=1,nmo2
          write(*,'(100G10.2)') (REAL(c_Ib_occ(i,j)),j=1,nocc)
        enddo

        print*,
        print*,"Check normalization:"
        do i=1,nocc
         sum_column = ZERO
         do j=1,nmo2
           sum_column = sum_column + conjg(c_Ib_occ(j,i))*c_Ib_occ(j,i)
         end do
         write(*,*) sum_column
        end do

        ! Transform the coefficient matrix back into the atomic orbital basis 1.
        ! This new coefficient matrix gives the transformation from the AOs of basis 1
        ! to the IAOs.
        print*,
        print*,"Write CIb into dfcoef file"
        allocate(c_Ib_ao(nao1,nmo2,2))
        c_Ib_ao(:,:,1) = matmul(cspinor_original%coeff(:,:,1),c_Ib)
        c_Ib_ao(:,:,2) = matmul(cspinor_original%coeff(:,:,2),c_Ib)
        call write_dffiles("iao.dfcoef","iao.dfpcmo",ao_basis_1,c_Ib_ao,nao1,nmo2,"MOLECULE.dfcoef")

        ! Separate the IAO DFCOEF file for each fragment. For instance, it will be interesting to get the IAOs a the metallic atom in a transition metal complex. This can be used in the AVAS theory (https://arxiv.org/abs/1701.07862)
        imo=0
        current_atom => first_atom
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           allocate(c_Ib_frag(nao1,nmoi,2))
           c_Ib_frag = c_Ib_ao(:,imo+1:imo+nmoi,:)
           write (fragfile_dfcoef,'(i3.3,a11)') nint(current_atom%bascharge),'_iao.dfcoef'
           write (fragfile_dfpcmo,'(i3.3,a11)') nint(current_atom%bascharge),'_iao.dfpcmo'
           call write_dffiles(fragfile_dfcoef,fragfile_dfpcmo,ao_basis_1,c_Ib_frag,nao1,nmoi,"MOLECULE.dfcoef")
           deallocate(c_Ib_frag)
           ! Do it also for valence orbitals only
           call atom_minbas_size (nint(current_atom%bascharge),ncore,nvale)
           allocate(c_Ib_frag(nao1,nvale,2))
           c_Ib_frag = c_Ib_ao(:,imo+1+ncore:imo+nmoi,:)
           write (fragfile_val_dfcoef,'(i3.3,a15)') nint(current_atom%bascharge),'_iao_val.dfcoef'
           write (fragfile_val_dfpcmo,'(i3.3,a15)') nint(current_atom%bascharge),'_iao_val.dfpcmo'
           call write_dffiles(fragfile_val_dfcoef,fragfile_val_dfpcmo,ao_basis_1,c_Ib_frag,nao1,nvale,"MOLECULE.dfcoef")
           deallocate(c_Ib_frag)
           imo = imo + nmoi
           current_atom=>current_atom%next
        end do

!       clean up list of atoms
        call del_atoms(first_atom)

!      Check IAOs, python reference:
!      IbRdm = mdot(CIb.T, S1, hfr.Rdm, S1, CIb)
!      IbFock = mdot(CIb.T, hfr.Fock, CIb)
!      IbCoreH = mdot(CIb.T, CoreH, CIb)
!      print("*** Check of properties of the intrinsic basis:\n")
!      print(FR("Deviation from orthogonality", la.norm((np.eye(CIb.shape[1]) - mdot(CIb.T, S1, CIb))).flatten()))
!      print(FR("|ib><ib|i> - |i>", la.norm(mdot(CIb, CIb.T, S1, COcc) - COcc)))
!      print(FR("nElec/IB", np.trace(IbRdm)))
!      print(FR("Hartree-Fock energy/IB", .5*dot2(IbFock+IbCoreH, IbRdm) + Atoms.fCoreRepulsion()))

        allocate(HF_RDM(nmo1,nmo1)) ! In the MO basis, should be 0 everywhere except 1 on part of the diagonal.
        HF_RDM = 0.D0
        do i=1,nmo1
         do j=1,nmo1
           do k=1,nocc
              HF_RDM(i,j) = HF_RDM(i,j) + conjg(coeff1_occ(i,k))*coeff1_occ(j,k)
           enddo
         enddo
        enddo
      
        allocate(IbRDM(nmo2,nmo2))
!       IbRDM = matmul(conjg(transpose(c_Ib)),matmul(S11,matmul(HF_RDM,matmul(S11,c_Ib))))
!       S11 is identity matrix in the MO basis...
        IbRDM = matmul(conjg(transpose(c_Ib)),matmul(HF_RDM,c_Ib))
        print*,
        print*,'-----1RDM in IAO basis-----(',nmo2,'x',nmo2,')'
        do i = 1, nmo2
           write(*,'(100G10.2)') (REAL(IbRDM(i,j)),j=1,nmo2)
        end do
        sum_diag = 0.0
        do i=1,nmo2
         sum_diag = sum_diag + IbRDM(i,i)
        enddo
        print*,"Sum of diagonal:",sum_diag

        print*,
        allocate(Id(nmo2,nmo2))
        Id = ZERO
        forall (i=1:nmo2) Id(i,i) = ONE
        print*,'TESTS'
        print*,"Deviation from orthogonality      :", norm2(real(Id - matmul(conjg(transpose(c_Ib)),matmul(S11,c_Ib))))
!        print*,"|ib><ib|i> - |i>                  :", norm2(real(matmul(c_Ib,matmul(conjg(transpose(c_Ib)),coeff1_occ))&
!                                                                  - coeff1_occ))
        print*,"(occupied MO) |ib><ib|i> - |i>    :", norm2(real(matmul(c_Ib,c_Ib_occ) - coeff1_occ))
        print*,"deviation in number of electrons  :", real(sum_diag - nocc)

!        MAKE IBO. I actually do not need that right ? We are only interested about the intrinsic fragment orbitals are we ? 
!        print*,
!        print*,'----***** MAKE IBOs *****-----'
!        call MakeIBO(c_Ib,coeff1_occ,S11,cspinor_fragment,n_fragments,nmo1,nmo2,nocc)

!     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') ) ! Useful format to write complex matrices.

   end program generate_iao
