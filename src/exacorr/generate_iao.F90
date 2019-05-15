

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
        complex(8), allocatable :: c_Ib(:,:)                  !  Matrix "A" in Knizia's paper: coefficient matrix of the IAOs.
        complex(8), allocatable :: coeff1_occ(:,:)

        integer :: nocc = 10 ! Should be read from a file, but currently this info is not available
        integer luxyz, ierr, i, j, k, l, nmo1, nmo2, imo, jmo, nmoi, nmoj, nao1, nao2

        complex(8), allocatable :: DMAT(:,:)
        complex(8) :: sum_diag
        character(13) :: target


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
           nao2 = nao2 + cspinor_fragment(i)%nao
           current_atom=>current_atom%next
        end do

        call read_from_dirac (ao_basis_1,cspinor_original,'DFCOEF',ierr)
        nmo1 = cspinor_original%nmo
        nao1 = cspinor_original%nao
        print*,"Read MO coefficients: nmo1,nmo2,nao1,nao2 (size of basis 1 and basis 2)",nmo1,nmo2,nao1,nao2

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
        print*,'Print real part of S22:'
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


        ! Calculate the molecule-intrinsic atomic orbital (IAO) basis
        ! Reference: Knizia, J. Chem. Theory Comput., http://dx.doi.org/10.1021/ct400687b
        ! See his reference Python implementation for the Simple/2014 scheme used here.
        allocate (S12(nmo1,nmo2))
        allocate (P12(nmo1,nmo2))
        print*,
        print*,'-----Compute S12-----(',nmo1,'x',nmo2,')'
        S12 = conjg(transpose(S21))
        print*,'Print real part of S12:'
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(S12(i,j)),j=1,nmo2)
        enddo

        ! Make the projector, zposv solves AX = B in which B is input and overwritten by X on output.
        print*,
        print*,'-----Compute P12-----(',nmo1,'x',nmo2,')'
        P12 = S12
        call zposv( 'U', nmo1, nmo2, S11, nmo1, P12, nmo1, ierr )
        print*,'Print real part of P12:'
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(P12(i,j)),j=1,nmo2)
        enddo

        allocate(coeff1_occ(nmo1,nocc))
        print*,
        print*,'-----Compute COcc-----(',nmo1,'x',nocc,') (diagonal matrix)'
        coeff1_occ = ZERO
        forall (i=1:nocc) coeff1_occ(i,i) = ONE

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
        do i=1,nmo1
           do j=1,nmo2
              if (REAL(c_Ib(i,j)) > 0.999999999999D0 .and. AIMAG(c_Ib(i,j)) < 1.0D-12) c_Ib(i,j) = ONE
              if (ABS(REAL(c_Ib(i,j))) < 1.0D-12 .and. ABS(AIMAG(c_Ib(i,j))) < 1.0D-12) c_Ib(i,j) = ZERO
           end do
        end do
        print*,'Print real part of CIb:'
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(c_Ib(i,j)),j=1,nmo2)
        enddo
        
        print*,
        print*,"-----Symmetrically orthogonalize intrinsic basis (in the mo basis)-----"
        call SymOrth(c_Ib,S11,size(S11(:,1)),size(c_Ib(1,:)))

!       filter out occupation numbers slightly broken due to rounding... At least it clarifies the output
        do i=1,nmo1
           do j=1,nmo2
              if (REAL(c_Ib(i,j)) > 0.999999999999D0 .and. AIMAG(c_Ib(i,j)) < 1.0D-12) c_Ib(i,j) = ONE
              if (ABS(REAL(c_Ib(i,j))) < 1.0D-12 .and. ABS(AIMAG(c_Ib(i,j))) < 1.0D-12) c_Ib(i,j) = ZERO
           end do
        end do
        print*,'Print real part of CIb after orthogonalization'
        do i=1,nmo1
           write(*,'(100G10.2)')(REAL(c_Ib(i,j)),j=1,nmo2)
        enddo

        print*,
        print*,'-----Compute the 1RDM (IAO basis)-----'
        allocate(DMAT(nmo1,nmo1))
        DMAT = 0.D0
        do i=1,nmo1
         do j=1,nmo1
           do k=1,nocc
              DMAT(i,j) = DMAT(i,j) + c_Ib(i,k)*c_Ib(j,k)
           enddo
         enddo
        enddo
        print*,"PRINT DIAGONAL OF DENSITY MATRIX"
        write(*,'(100G10.2)') (REAL(DMAT(i,i)),i=1,nmo1)
        sum_diag = 0.0
        do i=1,nmo1
         sum_diag = sum_diag + DMAT(i,i)
        enddo
        print*,"Sum of diagonal (should to be equal to nocc = 10):",sum_diag

        print*,
        print*,'----***** MAKE IBOs *****-----'
        call MakeIBO(c_Ib,coeff1_occ,S11,cspinor_fragment,n_fragments,nmo1,nmo2,nocc)


!       clean up list of atoms
        call del_atoms(first_atom)

!     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') ) ! Useful format to write complex matrices.

   end subroutine generate_iaos

   subroutine MakeIBO(A,C,S11,cspinor_fragment,n_fragments,nmo1,nmo2,nocc)
   ! make intrinsic bond orbitals (IBOs)

        use exacorr_mo

        implicit none

        integer :: i,j
        integer, intent(in)     :: n_fragments, nmo1, nmo2, nocc
        complex(8), intent(in)  :: S11(nmo1,nmo1)
        complex(8), intent(in)  :: C(nmo1,nocc)
        complex(8), intent(in)  :: A(nmo1,nmo2)
        type(cmo),  intent(in)  :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        complex(8) :: Aocc(nmo2,nocc)
        real(8) :: Aocc2(nmo2,nocc)
         
        print*,'-----Re-express occupied orbitals in terms of IAO basis-----'
        Aocc = matmul(conjg(transpose(A)),matmul(S11,C))     
        do i=1,nmo2
          write(*,'(100G10.2)') (REAL(Aocc(i,j)),j=1,nocc)
        enddo
        print*,
        print*,'-----Localize iteratively (check how to do that properly with complex entries?)-----'
        call LocalizeIbosInIaoBasis(Aocc,Aocc2,cspinor_fragment,n_fragments,nmo2,nocc)
        print*,
        print*,'-----Compute CIbOcc-----'
        do i=1,nmo2
          write(*,'(100G10.2)') (Aocc2(i,j),j=1,nocc)
        enddo
        stop 'MakeIBO subroutine not finished'

        ! Python ref:
        !CIbOcc = mdot(CIb.T, S1, COcc)
        ! localize iteratively (Pipek-Mezey-like) and order by energy.
        !LIbOcc = LocalizeIbosInIaoBasis(CIbOcc, MakeAtomIbOffsets(Atoms)[0])
        ! transform back localized orbitals in IAO basis to main basis,
        ! and order them by energy.
        !LOcc = OrderByEnergy(mdot(CIb, LIbOcc), Fock)

   end subroutine MakeIBO

   subroutine LocalizeIbosInIaoBasis(Aocc,Aocc2,cspinor_fragment,n_fragments,nmo2,nocc)
        ! For now, written for real numbers. I don't know yet how to write in for complex.
        use exacorr_mo

        implicit none
        integer, parameter :: Expnt = 4
        logical :: conv
        integer :: i,j,k,l, imo, nmoi, it
        integer, intent(in) :: n_fragments, nmo2, nocc
        complex(8), intent(in) :: Aocc(nmo2,nocc)
        real(8), intent(out) :: Aocc2(nmo2,nocc)
        type(cmo), intent(in) :: cspinor_fragment(n_fragments) ! MO-coefficients of the fragments
        real(8) :: phi,func, fGrad, Aij, Bij, Cii, Cij, Cjj
        real(8), allocatable :: C(:,:), Ai(:), Aj(:), vec(:)

        print*,"Entering LocalizeIbosInIaoBasis subroutine (using reals instead of complex)&
              &, not working yet"
        Aocc2 = REAL(Aocc)
        conv = .false.
        do it=1,100
          func = 0.D0
          imo = 0
          do i=1,n_fragments
             nmoi = cspinor_fragment(i)%nmo
             do j=1,nocc
                allocate(vec(nmoi))
                vec = Aocc2(imo+1:imo+nmoi,j)
                func = func + dot_product(vec,vec)**Expnt  
                deallocate(vec)
             enddo
             imo = imo + nmoi
          enddo
          print*,'func = ',func
  
          fGrad = 0.D0
          do i=1,nocc
             do j=1,i
                Aij = 0.D0
                Bij = 0.D0
                imo = 0
                do k=1,n_fragments
                   nmoi = cspinor_fragment(k)%nmo   
                   allocate(C(nmoi,nocc))
                   allocate(Ai(nmoi))
                   allocate(Aj(nmoi))
                   C = Aocc2(imo+1:imo+nmoi,:)
                   Ai = C(:,i)
                   Aj = C(:,j)
                   Cii = dot_product(Ai,Ai)
                   Cij = dot_product(Ai,Aj)
                   Cjj = dot_product(Aj,Aj)
                   if ((Expnt.eq.2) .or. (Expnt.eq.3)) then
                     Aij = Aij + 4.D0*Cij**2 - (Cii-Cjj)**2
                     Bij = Bij + 4*Cij*(Cii - Cjj)
                   else if (Expnt.eq.4) then
                     Bij = Bij + 4.D0*Cij*(Cii**3 - Cjj**3)
                     Aij = Aij - Cii**4 - Cjj**4 + &
                         & 6.D0*(Cii**2 + Cjj**2)*Cij**2 + Cjj*Cii**3 + &
                         & Cii*Cjj**3
                   else
                     write(*,*)"Expnt should be 2, 3 or 4."
                     stop
                   end if
                   deallocate(C,Ai,Aj)
                   imo = imo + nmoi
                enddo
                allocate(Ai(nmo2))
                allocate(Aj(nmo2))
                phi = 0.25*datan2(Bij,-Aij)
                fGrad = fGrad + Bij**2
                Ai = Aocc2(:,i)
                Aj = Aocc2(:,j)
                Aocc2(:,i) = dcos(phi)*Ai + dsin(phi)*Aj 
                Aocc2(:,j) = -dsin(phi)*Ai + dcos(phi)*Aj 
                deallocate(Ai,Aj)
             enddo
          enddo

          fGrad = fGrad**0.5
          if (abs(fGrad) .lt. 1.0D-8) then
             conv = .true.
             exit
          endif
        enddo
  
        if (.not. conv) then
           print*,"No convergence reached, value of fGrad = ", fGrad
        else
           print*,"Iteration ", it, " and fGrad = ", fGrad 
        endif

   end subroutine LocalizeIbosInIaoBasis

   subroutine MakeIAO_Simple_2014(A,C,S21,S22,P12,nmo1,nmo2,nocc)

        implicit none 
        integer :: ierr
        integer, intent(in)     :: nmo1, nmo2, nocc
        complex(8), allocatable :: coeff2_occ(:,:),c_Ib(:,:)
        complex(8), allocatable :: c_tilde(:,:),c_tilde_bar(:,:),c_tilde2_bar(:,:)
        complex(8), allocatable :: s_tilde(:,:),t4(:,:)
        complex(8), intent(in)  :: S21(nmo2,nmo1),S22(nmo2,nmo2) !  the overlap matrices
        complex(8), intent(in)  :: P12(nmo1,nmo2)          !  the projection operator
        complex(8), intent(in)  :: C(nmo1,nocc)
        complex(8), intent(out) :: A(nmo1,nmo2)

        ! Make intrinsic basis orbitals for the occupied space (Simple/2014 scheme of Knizia's reference implementation)
        ! Note that our overlap matrices refer to the MO basis of the molecule and its fragments and that the coefficient
        ! matrices are simple diagonal matrices. we transform to AO basis after this step.
        allocate (coeff2_occ(nmo2,nocc))
        allocate (c_tilde(nmo2,nocc))
        allocate (c_tilde_bar(nocc,nmo2))
        allocate (c_tilde2_bar(nmo2,nocc))
        allocate (s_tilde(nocc,nocc))
        allocate (t4(nmo1,nocc))

        coeff2_occ = matmul(S21,C)
        c_tilde = coeff2_occ
        call zposv( 'U', nmo2, nocc, S22, nmo2, c_tilde, nmo2, ierr )
        s_tilde = matmul(conjg(transpose(coeff2_occ)),c_tilde)
        c_tilde_bar = conjg(transpose(c_tilde))
        call zposv( 'U', nocc, nmo2, s_tilde, nmo2, c_tilde_bar, nmo2, ierr )
        c_tilde2_bar = conjg(transpose(c_tilde_bar))
        t4 = C - matmul(P12,c_tilde2_bar)
        A = P12 + matmul(t4,conjg(transpose(coeff2_occ)))

        ! Reference Python code for Simple/2014
        !COcc2 = mdot(S21, COcc)                # O(N m^2)
        !CTil = la.solve(S2, COcc2)             # O(m^3)
        !STil = mdot(COcc2.T, CTil)             # O(m^3)
        !CTil2Bar = la.solve(STil, CTil.T).T    # O(m^3)
        !T4 = COcc - mdot(P12, CTil2Bar)        # O(N m^2)
        !CIb = P12 + mdot(T4, COcc2.T)          # O(N m^2)

        deallocate (coeff2_occ,c_tilde,c_tilde_bar,c_tilde2_bar,s_tilde,t4)


   end subroutine MakeIAO_Simple_2014


   subroutine MakeIAO_Simple_2013(A,C,S11,S12,S21,S22,P12,nmo1,nmo2,nocc)

        implicit none
        integer :: ierr,i,j
        integer, intent(in)     :: nmo1, nmo2, nocc
        complex(8), intent(out) :: A(nmo1,nmo2)
        complex(8), intent(in)  :: S11(nmo1,nmo1),S12(nmo1,nmo2),S21(nmo2,nmo1),S22(nmo2,nmo2) !  the overlap matrices
        complex(8), intent(in)  :: P12(nmo1,nmo2)          !  the projection operator
        complex(8), intent(in)  :: C(nmo1,nocc)
        complex(8), allocatable :: c_tilde(:,:),P21(:,:)


        allocate (c_tilde(nmo1,nocc))
        allocate (P21(nmo2,nmo1))
        P21 = S21
        call zposv( 'U', nmo2, nmo1, S22, nmo2, P21, nmo2, ierr )

        c_tilde = matmul(P12,matmul(P21,C))
        call SymOrth(c_tilde,S11,nmo1,nocc)
        A = P12 + matmul(C,matmul(conjg(transpose(C)),S12))&
        &               - matmul(c_tilde,matmul(conjg(transpose(c_tilde)),S12))

        ! Reference Python code for Simple/2013
        !CTil = SymOrth(mdot(P12, P21, COcc), S1)
        !CIb = P12 + mdot(COcc, COcc.T, S12) - mdot(CTil, CTil.T, S12)

        deallocate (c_tilde,P21)

   end subroutine MakeIAO_Simple_2013

   subroutine MakeIAO_Standard_2013(A,C,S11,S21,S22,P12,nmo1,nmo2,nocc)

        complex(8), parameter :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        integer :: i
        integer, intent(in) :: nmo1,nmo2,nocc
        complex(8), intent(in) :: S11(nmo1,nmo1), S21(nmo2,nmo1), S22(nmo2,nmo2)
        complex(8), intent(in) :: P12(nmo1,nmo2)
        complex(8), intent(in) :: C(nmo1,nocc)
        complex(8), intent(out) :: A(nmo1,nmo2)
        complex(8), allocatable :: c_tilde(:,:),P21(:,:),Id(:,:)
        complex(8), allocatable :: term11(:,:),term12(:,:),term21(:,:),term22(:,:)
       
        allocate (c_tilde(nmo1,nocc))
        allocate (P21(nmo2,nmo1))
        allocate (Id(nmo1,nmo1))
        P21 = S21
        call zposv( 'U', nmo2, nmo1, S22, nmo2, P21, nmo2, ierr )

        Id = ZERO
        forall (i=1:nmo1) Id(i,i) = ONE
        c_tilde = matmul(P12,matmul(P21,C))
        call SymOrth(c_tilde,S11,nmo1,nocc)
        allocate(term11(nmo1,nmo1))
        allocate(term12(nmo1,nmo2))
        allocate(term21(nmo1,nmo1))
        allocate(term22(nmo1,nmo2))
        term11 = matmul(C,matmul(conjg(transpose(C)),S11))
        term12 = matmul(matmul(c_tilde,matmul(conjg(transpose(c_tilde)),S11)),P12)
        term21 = Id-matmul(C,matmul(conjg(transpose(C)),S11))
        term22 = matmul(Id-matmul(c_tilde,matmul(conjg(transpose(c_tilde)),S11)),P12)
        A = matmul(term11,term12) + matmul(term21,term22)       

      ! Reference Python code for Simple/2013
      ! Calculate de-polarized occupied orbitals CTil by projecting COcc from 1 to 2 and back (Knizia, eq.1)
      ! CTil = SymOrth(mdot(P12, P21, COcc), S1)
      ! this is Knizia, eq. 2.
      ! make [ |i><i| |itil><itil|  + (1-|i><i|)(1-|itil><itil|) ] P12 |mu>
      ! CIb = mdot(mdot(COcc, COcc.T, S1), mdot(CTil, CTil.T, S1), P12) +\
      !      mdot(np.eye(nBf)-mdot(COcc, COcc.T, S1), np.eye(nBf)-mdot(CTil, CTil.T, S1), P12)

        deallocate(P21,c_tilde,Id,term11,term12,term21,term22)

   end subroutine MakeIAO_Standard_2013

   subroutine SymOrth(C,S,nmo1,nmo2)
        !Symmetrically orthogonalize intrinsic basis to arrive at standard IAOs

        implicit none
        integer, intent(in) :: nmo1, nmo2
        complex(8), intent(in) :: S(nmo1,nmo1)
        complex(8), intent(inout) :: C(nmo1,nmo2)
        integer, parameter :: lwmax=1000
        integer :: info, lwork, liwork, lrwork,i,j
        integer, dimension(lwmax) :: iwork
        complex(8), dimension(lwmax) :: work
        real(8), dimension(lwmax) :: rwork
        complex(8), allocatable :: matrix(:,:)
        real(8), allocatable :: ew(:)

        allocate(matrix(nmo2,nmo2))
        allocate(ew(nmo2))
        matrix = matmul(conjg(transpose(C)),matmul(S,C))
!
!       Query the optimal workspace.
!
        lwork = -1
        liwork = -1
        lrwork = -1

        call zheevd( 'V', 'L', nmo2, matrix, nmo2, ew, work, lwork, rwork, &
                    lrwork, iwork, liwork, info )
        lwork = min(lwmax,int(work(1)))
        lrwork = min(lwmax,int(rwork(1)))
        liwork = min(lwmax,iwork(1))
!
!       Solve eigenproblem.
!
        call zheevd( 'V', 'L', nmo2, matrix, nmo2, ew, work, lwork, rwork, &
                    lrwork, iwork, liwork, info )
!
!       Check for convergence.
!
        if( info.gt.0 ) then
           write(*,*)'The algorithm failed to compute eigenvalues.'
           stop
        end if        

        do i=1,nmo2
           ew(i) = 1.0/sqrt(sqrt(ew(i)))
        enddo
        do i=1, nmo2
           matrix(:,i)=matrix(:,i)*ew(i)
        enddo
        C = matmul(C,matmul(matrix,conjg(transpose(matrix))))

        deallocate(matrix,ew)

   end subroutine SymOrth
   
   subroutine get_mo_overlap (ao_basis_bra,cspinor_bra,ao_basis_ket,cspinor_ket,SMO)

        use exacorr_datatypes
        use exacorr_mo

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)

        type(basis_set_info_t),intent(in) :: ao_basis_bra, ao_basis_ket 
        type(cmo),intent(in)              :: cspinor_bra, cspinor_ket 
        complex(8), intent(out)           :: SMO(cspinor_bra%nmo,cspinor_ket%nmo)

        real(8), allocatable :: SAO(:,:)   ! the AO overlap matrix
        complex(8), allocatable :: SC(:,:) ! Overlap times MO coefficients
        integer :: i,j

        ! allocate auxilliaries
        allocate (SAO(cspinor_bra%nao,cspinor_ket%nao))
        allocate (SC(cspinor_bra%nao,cspinor_ket%nmo))
!        print*,'    MO overlap (nao basis, nshells basis, nao cspinor) bra:', &
!                                 ao_basis_bra%nao, ao_basis_bra%nshells, cspinor_bra%nao, &
!                                 ' ket:', &
!                                 ao_basis_ket%nao, ao_basis_ket%nshells, cspinor_ket%nao

        call  get_ao_overlap (ao_basis_bra,ao_basis_ket,SAO)
!        print*,'Print AO overlap (size ', cspinor_bra%nao ,'x', cspinor_ket%nao,')'
!        do i=1,cspinor_bra%nao
!           write(*,*)(SAO(i,j),j=1,cspinor_ket%nao)
!        enddo

        SC  = matmul(SAO,cspinor_ket%coeff(:,:,1))                        ! multiply with alpha coefficients
        SMO = matmul(conjg(transpose(cspinor_bra%coeff(:,:,1))),SC)       ! alpha-alpha contribution
        SC  = matmul(SAO,cspinor_ket%coeff(:,:,2))                        ! multiply with beta coefficients
        SMO = SMO + matmul(conjg(transpose(cspinor_bra%coeff(:,:,2))),SC) ! add the beta-beta contribution

!       filter out occupation numbers slightly broken due to rounding.
        do i=1,cspinor_bra%nmo
           do j=1,cspinor_ket%nmo
              if (REAL(SMO(i,j)) > 0.999999999999D0 .and. AIMAG(SMO(i,j)) < 1.0D-12) SMO(i,j) = ONE
              if (ABS(REAL(SMO(i,j))) < 1.0D-12 .and. ABS(AIMAG(SMO(i,j))) < 1.0D-12) SMO(i,j) = ZERO
           end do
        end do

        deallocate(SAO,SC)
     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') )

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
