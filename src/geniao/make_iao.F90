module make_iao

contains

   subroutine MakeIAO_Simple_2014(A,C,S21,S22,P12,nmo1,nmo2,nocc)

        implicit none 
        integer :: ierr,i,j
        integer, intent(in)     :: nmo1, nmo2, nocc
        complex(8), allocatable :: coeff2_occ(:,:)
        complex(8), allocatable :: c_tilde(:,:),c_tilde_bar(:,:),c_tilde2_bar(:,:)
        complex(8), allocatable :: s_tilde(:,:),t4(:,:)
        complex(8), intent(in)  :: S21(nmo2,nmo1),S22(nmo2,nmo2) !  the overlap matrices
        complex(8), intent(in)  :: P12(nmo1,nmo2)          !  the projection operator
        complex(8), intent(in)  :: C(nmo1,nocc)
        complex(8), intent(out) :: A(nmo1,nmo2)
        complex(8), allocatable :: A_X(:,:), P21(:,:), S_12(:,:), c_tilde_X(:,:)

        ! Make intrinsic basis orbitals for the occupied space (Simple/2014 scheme of Knizia's reference implementation)
        ! Note that our overlap matrices refer to the MO basis of the molecule and its fragments and that the coefficient
        ! matrices are simple diagonal matrices. we transform to AO basis after this step.
        allocate (coeff2_occ(nmo2,nocc))
        allocate (c_tilde(nmo2,nocc))
        allocate (c_tilde_bar(nocc,nmo2))
        allocate (c_tilde2_bar(nmo2,nocc))
        allocate (s_tilde(nocc,nocc))
        allocate (t4(nmo1,nocc))


        ! use of zposv : CALL ZPOSV( 'Lower', N, NRHS, A, LDA, B, LDB, INFO ) with A(LDA,N) and B(LDB,NRHS). This solves AX=B and returns A = X.
        coeff2_occ = matmul(S21,C)
        c_tilde = coeff2_occ
        call zposv( 'U', nmo2, nocc, S22, nmo2, c_tilde, nmo2, ierr )
        s_tilde = matmul(conjg(transpose(coeff2_occ)),c_tilde)
        c_tilde_bar = conjg(transpose(c_tilde))
        call zposv( 'U', nocc, nmo2, s_tilde, nocc, c_tilde_bar, nocc, ierr )
        c_tilde2_bar = conjg(transpose(c_tilde_bar))
        t4 = C - matmul(P12,c_tilde2_bar)
        A = P12 + matmul(t4,conjg(transpose(coeff2_occ)))

        ! Compare with reference result to check (just to check if transformation is okay). 
        ! I should get rid of this once we are sure the transformation is correct.
        print*,
        print*," *** Check if Simple/2014 transformation is okay ***"
        allocate(A_X(nmo1,nmo2)) 
        allocate(P21(nmo2,nmo1)) 
        allocate(S_12(nmo1,nmo2)) 
        allocate(c_tilde_X(nmo2,nocc))
        S_12 = conjg(transpose(S21))
        P21 = S21
        call zposv( 'U', nmo2, nmo1, S22, nmo2, P21, nmo2, ierr )
        c_tilde_X = matmul(P21, C)
        call SymOrth(c_tilde_X,S22,nmo2,nocc)
        A_X = P12 + matmul(C, matmul(conjg(transpose(C)), S_12)) - &
                    matmul(P12, matmul(c_tilde_X, matmul(conjg(transpose(c_tilde_X)), S22)))
        if (norm2([norm2(real(A_X - A)),norm2(aimag(A_X-A))]) > 1.D-12) then
          print*," *** norm2(A_X - A) = ",norm2([norm2(real(A_X - A)),norm2(aimag(A_X-A))]),&
                 ", Transformation NOT OKAY (or the test has failed... !) ***"
        else
          print*," *** norm2(A_X - A) = ",norm2([norm2(real(A_X - A)),norm2(aimag(A_X-A))]),&
                 ", Transformation OKAY ***"
        endif

        deallocate (coeff2_occ,c_tilde,c_tilde_bar,c_tilde2_bar,s_tilde,t4,A_X,P21,S_12,c_tilde_X)


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

        deallocate(P21,c_tilde,Id,term11,term12,term21,term22)

   end subroutine MakeIAO_Standard_2013

   subroutine SymOrth(C,S,nmo1,nmo2)
        !Symmetrically orthogonalize orbitals C with respect to overlap matrix S (i.e., such that Cnew^T S Cnew == id).
        use utils

        implicit none
        complex(8), parameter        :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        real(8), parameter           :: ThrDel=1.D-7
        integer                      :: i,j, iDel
        integer,    intent(in)       :: nmo1, nmo2
        complex(8), intent(in)       :: S(nmo1,nmo1)
        complex(8), intent(inout)    :: C(nmo1,nmo2)
        complex(8), allocatable      :: matrix(:,:),matrix2(:,:),Id(:,:)
        real(8),    allocatable      :: ew(:)

        allocate(matrix(nmo2,nmo2))
        allocate(ew(nmo2))
        allocate(Id(nmo2,nmo2))
        Id = ZERO
        forall(i=1:nmo2) Id(i,i)=ONE
        matrix = matmul(conjg(transpose(C)),matmul(S,C))

        call diag_matrix(matrix,nmo2,ew)

        iDel = 0
        do while (ew(iDel+1) < ThrDel)
          iDel = iDel + 1
        enddo
        if (iDel == 0) then
          do i=1,nmo2
            ew(i) = 1.0/sqrt(sqrt(ew(i)))
          enddo
          do i=1,nmo2
            matrix(:,i)=matrix(:,i)*ew(i)
          enddo
          C = matmul(C,matmul(matrix,conjg(transpose(matrix))))
        else
          print*,"   ** WARNING: deleted", iDel, " eigenvalues in S^{-1/2} construction."
          allocate(matrix2(nmo2,nmo2-iDel))
          do i=1,nmo2-iDel
            ew(i+iDel) = 1.0/sqrt(sqrt(ew(i+iDel)))
          enddo
          do i=1,nmo2-iDel
            matrix2(:,i)=matrix(:,i+iDel)*ew(i+iDel)
          enddo
          C = matmul(C,matmul(matrix2,conjg(transpose(matrix2))))
          deallocate(matrix2)
        endif

        if (norm2([norm2(real(Id - matmul(conjg(transpose(C)),matmul(S,C)))), &
                  norm2(aimag(Id - matmul(conjg(transpose(C)),matmul(S,C))))]) > 1.D-12) then
         print*,"**WARNING** Symmetrical orthogonalization has failed (thrsh: 1.D-12)."
        end if

        deallocate(matrix,ew,Id)

   end subroutine SymOrth

end module make_iao
