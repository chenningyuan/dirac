!Read MO coefficients and generate intrinsic atomic orbitals

   subroutine generate_iaos

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_utils
        use read_xyzfile

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
        complex(8), allocatable :: coeff1_occ(:,:),coeff2_occ(:,:),c_Ib(:,:)
        complex(8), allocatable :: c_tilde(:,:),c_tilde_bar(:,:),c_tilde2_bar(:,:)
        complex(8), allocatable :: s_tilde(:,:),t4(:,:)

        integer :: nocc = 10 ! Should be read from a file, but currently this info is not available
        integer luxyz, ierr, i, j, nmo1, nmo2, imo, jmo, nmoi, nmoj

        print*,"Arrived in generate_iaos"

        call get_free_fileunit(luxyz)
        open (luxyz,file='H2O.xyz',form='formatted')
        call get_atoms(luxyz,first_atom)
        close (luxyz,status='keep')
        current_atom => first_atom
        print*,"Molecular structure used to make reference basis"
        ! For now we will just atomic fragments, which should be available from a library of DFCOEFs
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

        call read_from_dirac (ao_basis_1,cspinor_original,'DFCOEF',ierr)
        nmo1 = cspinor_original%nmo
        print*,"Read MO coefficients",nmo1,nmo2

        allocate (S11(nmo1,nmo1))
        call get_mo_overlap (ao_basis_1,cspinor_original,ao_basis_1,cspinor_original,S11)

        stop 'debug'
        allocate (S21(nmo2,nmo1))
        imo = 0
        do i = 1, n_fragments
           nmoi = cspinor_fragment(i)%nmo
           allocate (S_block(nmoi,nmo1))
           call get_mo_overlap (ao_basis_2(i),cspinor_fragment(i),ao_basis_1,cspinor_original,S_block)
           S21(:,imo+1:imo+nmoi) = S_block
           imo = imo + nmoi
           deallocate (S_block)
        enddo

        allocate (S22(nmo2,nmo2))
        jmo = 0
        do j = 1, n_fragments
           nmoj = cspinor_fragment(j)%nmo
           imo = 0
           do i = 1, n_fragments
              nmoi = cspinor_fragment(i)%nmo
              allocate (S_block(nmoi,nmoj))
              call get_mo_overlap (ao_basis_2(i),cspinor_fragment(i),ao_basis_2(j),cspinor_fragment(j),S_block)
              S22(imo+1:imo+nmoi,jmo+1:jmo+nmoj) = S_block
              S22(jmo+1:jmo+nmoj,imo+1:imo+nmoi) = conjg(transpose(S_block))
              imo = imo + nmoi
              deallocate (S_block)
           enddo
           jmo = jmo + nmoj
        enddo

        ! Calculate the molecule-intrinsic atomic orbital (IAO) basis
        ! Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b
        ! See his reference Python implementation for the Simple/2014 scheme used here.
        allocate (S12(nmo1,nmo2))
        allocate (P12(nmo1,nmo2))
        allocate (P21(nmo2,nmo1))
        S12 = conjg(transpose(S21))
        ! Make the projector, zposv solves AX = B in which B is input and overwritten by X on output.
        P12 = S12
        call zposv( 'U', nmo1, nmo2, S11, nmo1, P12, nmo1, ierr )

        ! Make intrinsic basis orbitals for the occupied space (Simple/2014 scheme of Knizia's reference implementation)
        ! Note that our overlap matrices refer to the MO bases of the molecule and its fragments and that the coefficient
        ! matrices are simple diagonal matrices. We transform to AO bases after this step.
        allocate (coeff1_occ(nmo1,nocc))
        allocate (coeff2_occ(nmo2,nocc))
        allocate (c_tilde(nmo2,nocc))
        allocate (c_tilde_bar(nocc,nmo2))
        allocate (c_tilde2_bar(nmo2,nocc))
        allocate (s_tilde(nocc,nocc))
        allocate (t4(nmo1,nocc))
        allocate (c_Ib(nmo1,nmo2))

        coeff1_occ = ZERO
        forall (i=1:nocc) coeff1_occ(i,i) = ONE
        coeff2_occ = matmul(S21,coeff1_occ)
        c_tilde = coeff2_occ
        call zposv( 'U', nmo2, nocc, S22, nmo2, c_tilde, nmo2, ierr )
        c_tilde_bar = conjg(transpose(c_tilde))
        s_tilde = matmul(c_tilde_bar,c_tilde)
        call zposv( 'U', nocc, nmo2, s_tilde, nmo2, c_tilde_bar, nmo2, ierr )
        c_tilde2_bar = conjg(transpose(c_tilde_bar))
        t4 = coeff1_occ - matmul(P12,c_tilde2_bar)
        c_Ib = P12 + matmul(t4,conjg(transpose(coeff2_occ)))

        ! Reference Python code
        !COcc2 = mdot(S21, COcc)                # O(N m^2)
        !CTil = la.solve(S2, COcc2)             # O(m^3)
        !STil = mdot(COcc2.T, CTil)             # O(m^3)
        !CTil2Bar = la.solve(STil, CTil.T).T    # O(m^3)
        !T4 = COcc - mdot(P12, CTil2Bar)        # O(N m^2)
        !CIb = P12 + mdot(T4, COcc2.T)          # O(N m^2)

        
      ! # calculate de-polarized occupied orbitals CTil by projecting COcc from 1 to 2 and back ([1], eq.1)
      ! CTil = SymOrth(mdot(P12, P21, COcc), S1)
      ! # this is [1] eq. 2.
      ! # make [ |i><i| |itil><itil|  + (1-|i><i|)(1-|itil><itil|) ] P12 |mu>
      ! CIb = mdot(mdot(COcc, COcc.T, S1), mdot(CTil, CTil.T, S1), P12) +\
      !       mdot(np.eye(nBf)-mdot(COcc, COcc.T, S1), np.eye(nBf)-mdot(CTil, CTil.T, S1), P12)


!       clean up list of atoms
        call del_atoms(first_atom)
        stop 'generated local orbitals'

   end subroutine generate_iaos
   
   subroutine get_mo_overlap (ao_basis_bra,cspinor_bra,ao_basis_ket,cspinor_ket,SMO)

        use exacorr_datatypes
        use exacorr_mo

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis_bra, ao_basis_ket 
        type(cmo),intent(in)              :: cspinor_bra, cspinor_ket 
        complex(8), intent(out)           :: SMO(cspinor_bra%nmo,cspinor_ket%nmo)

        real(8), allocatable :: SAO(:,:)   ! the AO overlap matrix
        complex(8), allocatable :: SC(:,:) ! Overlap times MO coefficients

        ! allocate auxilliaries
        allocate (SAO(cspinor_bra%nao,cspinor_ket%nao))
        allocate (SC (cspinor_bra%nao,cspinor_ket%nmo))

        call  get_ao_overlap (ao_basis_bra,ao_basis_ket,SAO)

        SC  = matmul(SAO,cspinor_ket%coeff(:,:,1))                        ! multiply with alpha coefficients
        SMO = matmul(conjg(transpose(cspinor_bra%coeff(:,:,1))),SC)       ! alpha-alpha contribution
        SC  = matmul(SAO,cspinor_ket%coeff(:,:,2))                        ! multiply with beta coefficients
        SMO = SMO + matmul(conjg(transpose(cspinor_bra%coeff(:,:,2))),SC) ! add the beta-beta contribution

        deallocate(SAO,SC)

   end subroutine get_mo_overlap

   subroutine get_ao_overlap (ao_basis_bra,ao_basis_ket,ao_overlap)

        use exacorr_datatypes

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
                              ncentr,atomic_data                              )
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
        integer ierr, ishell, nao, nmo, ncore, nvale

        ! Should include path information, for now assume this is copied into the working dir
        write (fragfile,'(i3.3,a7)') atom_number,'.dfcoef'
        call read_from_dirac (ao_basis,cspinor_full,fragfile,ierr)
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
        ! We don't need this coefficient set anymore, free up the memory
        call dealloc_mo(cspinor_full)

   end subroutine get_atomic_fragment

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
