   subroutine generate_iao

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_utils
        use read_xyzfile
        use get_overlap
        use make_iao

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
        complex(8), allocatable :: c_Ib(:,:)                  !  coefficient matrix to go from the MOs to the IAOs.
        complex(8), allocatable :: c_Ib_ao(:,:,:)             !  coefficient matrix to go from the AOs to the IAOs..
        complex(8), allocatable :: coeff1_occ(:,:)
        real(8), allocatable :: coeff_quaternion(:,:,:)
        real(8), allocatable :: SAO11(:,:)

        integer :: nocc = 10 ! Should be read from a file, but currently this info is not available
        integer luxyz, ierr, i, j, k, l, nmo1, nmo2, imo, jmo, nmoi, nmoj, nao1

        complex(8), allocatable :: DMAT(:,:)
        complex(8) :: sum_diag
        character(13) :: target


        print*,"Arrived in generate_iao"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           READ Input Files         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call get_free_fileunit(luxyz)
        open (luxyz,file='H2O.xyz',form='formatted')
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

        call read_from_dirac (ao_basis_1,cspinor_original,'DFCOEF_mol',ierr)

        nmo1 = cspinor_original%nmo
        nao1 = cspinor_original%nao
        print*,"Read MO coefficients: nmo1,nmo2(size of basis 1 and basis 2)",nmo1,nmo2

        ! The following Coeff file is here to check that I indeed get the same format as the DFPCMO_mol file.
        ! It will be removed when I'll be sure it is indeed correct in any cases.
        allocate(coeff_quaternion(cspinor_original%nao,int(cspinor_original%nmo/2),4))
        call convert_complex_to_quaternion_format(cspinor_original%coeff,coeff_quaternion,cspinor_original%nao,cspinor_original%nmo)
        open ( UNIT = 96, FILE = 'Coeff', ACCESS = 'SEQUENTIAL' )
        write(96,'(6F22.16)') coeff_quaternion
        close(96)
        deallocate(coeff_quaternion)


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
        call get_command_argument(1, target)
        select case (target)
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

   end subroutine generate_iao

   subroutine write_dffiles(filename_unformatted,filename_formatted,coeff,basis1,basis2)

       use exacorr_utils
       implicit none

       integer, parameter :: nsym = 1
       character(len=74) :: text="IAO coefficients"
       character(10) :: filename_unformatted, filename_formatted
       integer :: i,j, iunit, idim(3,nsym)
       integer, intent(in)     :: basis1, basis2
       complex(8), intent(in)  :: coeff(basis1,basis2,2)
       complex(8)  :: coeff_quaternion(basis1,int(basis2/2),4)
       real(8) :: toterg
       real(8) :: eig(basis2)
       integer :: ibeig(basis2)

       call get_free_fileunit(iunit)

       toterg= 0.D0
       eig = 0.D0
       ibeig = 0

       idim(3,1) = basis1
       idim(1,1) = 0
       idim(2,1) = int(basis2/2)

       ! I think that I need to have this format for the DFCOEF file ?
       call convert_complex_to_quaternion_format(coeff,coeff_quaternion,basis1,basis2)

       open (UNIT = iunit,file=filename_unformatted,FORM='UNFORMATTED',access='SEQUENTIAL')
       write (iunit) text,nsym,((idim(i,j),i=1,3),j=1,nsym),toterg
       write (iunit) coeff_quaternion
       write (iunit) eig
       write (iunit) ibeig
       close (iunit)

       open (UNIT = iunit,file=filename_formatted,FORM='FORMATTED',access='SEQUENTIAL')
       write(iunit,'(A74)') text
       write(iunit,'(7(X,I0))') nsym,((idim(i,j),i = 1,3),j=1,nsym)
       write(iunit,'(E24.16)') toterg
       write(iunit,'(6F22.16)') coeff_quaternion
       write(iunit,'(6E22.12)') eig
       write(iunit,'(66(X,I0))') ibeig


   end subroutine write_dffiles

   subroutine get_atomic_fragment (atom_number,atom_coord,ao_basis,cspinor)

        use exacorr_datatypes
        use exacorr_mo

        implicit none

        integer, intent(in) :: atom_number
        real(8), intent(in) :: atom_coord(3)
        type(basis_set_info_t),intent(out) :: ao_basis ! the ao basis used in the generation of the MO coefficients
        type(cmo),intent(out) :: cspinor               ! occupied spinors (minimal basis set)

        type(cmo):: cspinor_full  ! all spinors as they are read from file
        character(len=10) fragfile
        integer ierr, ishell, nao, nmo, ncore, nvale, i,j

        ! Should include path information, for now assume this is copied into the working dir
        write (fragfile,'(i3.3,a7)') atom_number,'.dfcoef'
        call read_from_dirac (ao_basis,cspinor_full,fragfile,ierr)
        print*,'  ierr (0 --> NO ERROR) = ', ierr
        ! All atomic coefficients are generated with the atom in the origin, translate the basis to the actual position
        do ishell = 1, ao_basis%nshells ! somehow the compiler wants an explicit loop for this
            ao_basis%gtos(ishell)%coord = atom_coord
        end do
        ! Query the size of the (mo) basis that should be made
        call atom_minbas_size (atom_number,ncore,nvale)
        nmo = ncore+nvale 
        nao = cspinor_full%nao
        ! Copy only the relevant MOs to the output array
        call alloc_mo (cspinor,nao,nmo)
        cspinor%energy = cspinor_full%energy(1:nmo)
        cspinor%coeff  = cspinor_full%coeff(:,1:nmo,:)
        cspinor%index  = cspinor_full%index(1:nmo)
        ! we don't need this coefficient set anymore, free up the memory
        call dealloc_mo(cspinor_full)

   end subroutine get_atomic_fragment

   subroutine convert_complex_to_quaternion_format(scoeff,coeff,nao,nmo)

      implicit none
      integer,    intent(in   ) :: nao, nmo
      complex(8), intent(in   ) :: scoeff(nao,nmo,2)
      real(8),    intent(  out) :: coeff(nao,int(nmo/2),4)
      integer                   :: i, j, spin

      do spin = 1, 2
        if (spin == 1) then
            do j = 1, int(nmo/2)
                do i = 1, nao
                coeff(i,j,1) = real(scoeff(i,2*j-1,spin))
                coeff(i,j,2) = aimag(scoeff(i,2*j-1,spin))
                coeff(i,j,3) = real(scoeff(i,2*j,spin))
                coeff(i,j,4) = aimag(scoeff(i,2*j,spin))
                end do
            end do
!        else if (spin == 2) then
!            do j = 1, int(nmo/2)
!                do i = 1, nao
!                coeff(i,j,3) = - real(scoeff(i,2*j-1,spin))
!                coeff(i,j,4) = aimag(scoeff(i,2*j-1,spin))
!                coeff(i,j,1) = real(scoeff(i,2*j,spin))
!                coeff(i,j,2) = - aimag(scoeff(i,2*j,spin))
!                end do
!            end do
        end if
      end do

   end subroutine convert_complex_to_quaternion_format

   subroutine atom_minbas_size(atom_number,ncore,nvale)

!        returns size of minimal basis (NB: counting spinors, so there is a factor of 2 wrt orbitals)

         implicit none

         integer,intent(in)  ::  atom_number  ! the number that is used to search the basis set and ao-coefficient libs
         integer,intent(out) ::  ncore, nvale ! the number of core and valence spinors

!        Easiest is to simply hardwire this

         select case (atom_number)
            case (:0)        ! 0 is to be used for point charges that do not carry a basis set
                 ncore = 0
                 nvale = 0
            case (1:2)       ! H, He
                 ncore = 0
                 nvale = 2
            case (3:4)       ! Li,Be
                 ncore = 2
                 nvale = 2
            case (5:10)      ! B..Ne
                 ncore = 2
                 nvale = 8
            case (11:12)     ! Na,Mg
                 ncore = 10
                 nvale = 2
            case (13:18)     ! Al..Ar
                 ncore = 10
                 nvale = 8
            case (19:20)     ! K,Ca
                 ncore = 18
                 nvale = 2
            case (21:30)     ! Sc..Zn
                 ncore = 18
                 nvale = 12
            case (31:36)     ! Ga..Kr
                 ncore = 18
                 nvale = 18
            case (37:38)     ! Rb,Sr
                 ncore = 36
                 nvale = 2
            case (39:48)     ! Y..Cd
                 ncore = 36
                 nvale = 12
            case (49:54)     ! In..Xe
                 ncore = 36
                 nvale = 18
            case (55:56)     ! Cs,Ba
                 ncore = 54
                 nvale = 2
            case (57:70)     ! La..Yb
                 ncore = 54
                 nvale = 16  ! assuming that we will enforce 4f rather than 5d occupation in the atomic scf
            case (71:80)     ! Lu..Hg
                 ncore = 54
                 nvale = 26
            case (81:86)     ! Tl..Rn
                 ncore = 68  ! place the 4f in the core
                 nvale = 18
            case (87:88)     ! Fr,Ra
                 ncore = 86 
                 nvale = 2
            case (89:102)    ! Ac,RNo
                 ncore = 86 
                 nvale = 16  ! assuming that we will enforce 5f rather than 6d occupation in the atomic scf
            case (103:112)   ! Lr..Cn
                 ncore = 86
                 nvale = 26
            case (113:118)   ! Nh..Og
                 ncore = 100 ! place the 5f in the core
                 nvale = 18
            case (119:)      ! Venturing into the unknown superheavies, impossible to decide the valence
                 ncore = 118
                 nvale = atom_number - ncore
         end select
     
   end subroutine atom_minbas_size
