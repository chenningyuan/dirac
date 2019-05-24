module make_iao

contains

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
        integer,    intent(in)       :: nmo1, nmo2
        complex(8), intent(in)       :: S(nmo1,nmo1)
        complex(8), intent(inout)    :: C(nmo1,nmo2)
        integer,    parameter        :: lwmax=1000
        integer                      :: info, lwork, liwork, lrwork,i,j
        integer,    dimension(lwmax) :: iwork
        complex(8), dimension(lwmax) :: work
        real(8),    dimension(lwmax) :: rwork
        complex(8), allocatable      :: matrix(:,:)
        real(8),    allocatable      :: ew(:)

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

end module make_iao
