!ajt dec08 Contains slightly generalized versions of the
!          operations in modulematrix_operations.
!
!Mainly intended as hidden back-end for module matrix_defop.
!
!mat_axtpby  combines  mat_trans, mat_scal, mat_copy,
!                      mat_assign, and mat_add
!
!mat_gemm  is  mat_mul
!
!mat_dot  combines  mat_dotproduct and mat_trAB
!
!mat_init  is generalized so that another matrix can fill
!          the role of nrow and ncol.

module matrix_genop_new

#ifndef PRG_DIRAC
    use matrix_operations,only: matrix,mat_mul, &
          & mat_init_full    => mat_init, &
          & mat_free_nonzero => mat_free
#endif /* PRG_DIRAC */  

    use matrix_genop

    implicit none

    public matrix
    public mat_init
    public mat_free
    public mat_same
    public mat_isdef
    public mat_iszero
    public mat_axtpby
    public mat_gemm
    public mat_dot
    public mat_trace
    public mat_print

    private

!   interface mat_init
!       module procedure mat_init_full
!       module procedure mat_init_from
!   end interface

contains

  subroutine mat_dense_mul(a,b,transa,transb,alpha,beta,c)
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     character, intent(in)    :: transa, transb
     REAL(8), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c
     INTEGER                  :: m,n,k

     if (transa == 'n' .or. transa == 'N') then
       m = a%nrow
       k = a%ncol
     elseif (transa == 't' .or. transa == 'T') then
       m = a%ncol
       k = a%nrow
     else
       print*,'unknown format in mat_dense_mul'
       STOP 'unknown format in mat_dense_mul'
     endif
     if (transb == 'n' .or. transb == 'N') then
       n = b%ncol
       if (b%nrow /= k) STOP 'weird'
     elseif (transb == 't' .or. transb == 'T') then
       n = b%nrow
       if (b%ncol /= k) STOP 'weird'
     endif
!write (mat_lu,*) 'matrix a'
!call mat_dense_print(a, 1, a%nrow, 1, a%ncol, mat_lu)
!write (mat_lu,*) 'matrix b'
!call mat_dense_print(b, 1, b%nrow, 1, b%ncol, mat_lu)
!write (mat_lu,*) 'matrix c'
!call mat_dense_print(c, 1, c%nrow, 1, c%ncol, mat_lu)
     call DGEMM(transa,transb,m,n,k,alpha,&
               &a%elms,a%nrow,b%elms,b%nrow,beta,c%elms,c%nrow)

  end subroutine mat_dense_mul

      SUBROUTINE mat_mul(a, b, transa, transb, alpha, beta, c)
         !c = alpha*ab + beta*c
         !transa = 'T'/'t' - transposed, 'N'/'n' - normal
         implicit none
         TYPE(Matrix), intent(IN) :: a, b
         character, intent(in)    :: transa, transb
         REAL(8), INTENT(IN)  :: alpha, beta
         TYPE(Matrix), intent(inout):: c
         integer :: ak, bk, ci, cj
  !      if (INFO_TIME_MAT) CALL TIMER('START ',mat_TSTR,mat_TEN)
        !no_of_matmuls = no_of_matmuls + 1
         if (transa == 'n' .or. transa == 'N') then
           ak = a%ncol
           ci = a%nrow
         elseif (transa == 't' .or. transa == 'T') then
           ak = a%nrow
           ci = a%ncol
         endif
         if (transb == 'n' .or. transb == 'N') then
           bk = b%nrow
           cj = b%ncol
         elseif (transb == 't' .or. transb == 'T') then
           bk = b%ncol
           cj = b%nrow
         endif
         if (ak /= bk .or. ci /= c%nrow .or. cj /= c%ncol) then
           STOP 'wrong dimensions in mat_mul or unknown trans possibility'
         endif
!        select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_mul(a,b,transa, transb,alpha,beta,c)
!        case(mtype_dense)
             call mat_dense_mul(a,b,transa, transb,alpha,beta,c)
!        case(mtype_sparse1)
!            call mat_sparse1_mul(a,b,transa, transb,alpha,beta,c)
!ifdef HAVE_BSM
!        case(mtype_sparse_block)
!            call bsm_mul(a,b,transa, transb,alpha,beta,c)
!endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_mul(a,b,transa, transb,alpha,beta,c)
!        case(mtype_unres_dense)
!            call mat_unres_dense_mul(a,b,transa, transb,alpha,beta,c)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_mul(a,b,transa, transb,alpha,beta,c)
!        case default
!             stop "mat_mul not implemented for this type of matrix"
!        end select
!        if (INFO_TIME_MAT) CALL TIMER('MATMUL ',mat_TSTR,mat_TEN)
      END SUBROUTINE mat_mul

  subroutine mat_dense_init(a,nrow,ncol)
     implicit none
     TYPE(Matrix) :: a
     integer, intent(in) :: nrow, ncol
     integer :: nsize
!A%elms should not be associated with anything.
!if so, the memory will be lost!!
     NULLIFY(a%elms)
     a%nrow = nrow
     a%ncol = ncol
     nsize = A%nrow*A%ncol
     ALLOCATE(a%elms(nsize))
!    call dens_stat_allocated_memory(nsize)
  end subroutine mat_dense_init

  subroutine mat_dense_free(a)
     implicit none
     TYPE(Matrix) :: a
     integer :: nsize
     if (.not.ASSOCIATED(a%elms)) then
       print*,'memory previously released!!'
       STOP 'Error in mat_dense_free - memory previously released'
     endif
     nsize = SIZE(a%elms)
     DEALLOCATE(a%elms)
     NULLIFY(a%elms)
!    call dens_stat_deallocated_memory(nsize)
  end subroutine mat_dense_free

      SUBROUTINE mat_free_nonzero(a)
         implicit none
         TYPE(Matrix) :: a
!        select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_free(a)
!        case(mtype_dense)
             call mat_dense_free(a)
!        case(mtype_sparse1)
!            call mat_sparse1_free(a)
!ifdef HAVE_BSM
!        case(mtype_sparse_block)
!            call bsm_free(a)
!endif
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_free(a)
!        case(mtype_unres_dense)
!            call mat_unres_dense_free(a)
!         case(mtype_unres_sparse1)
!             call mat_unres_sparse1_free(a)
!        case default
!             stop "mat_free not implemented for this type of matrix"
!        end select

      !free auxaliary data
      if (ASSOCIATED(a%iaux)) deallocate(a%iaux)
      if (ASSOCIATED(a%raux)) deallocate(a%raux)

      END SUBROUTINE

    subroutine mat_free(A)
        type(matrix),intent(inout) :: A
        !if (associated(A%elms)) print *,'mat_free',loc(A),' elms=',loc(A%elms)
        if (associated(A%elms)) call mat_free_nonzero(A)
        A%nrow = -1; A%ncol = -1
    end subroutine


!   subroutine mat_init_from(C,A,ta,B,tb,zero,alias,reset)
!       type(matrix),         intent(inout) :: C
!       type(matrix),optional,intent(in)    :: A,B
!       character*1, optional,intent(in)    :: ta,tb
!       logical,     optional,intent(in)    :: zero,alias,reset
!       logical                             :: taa,tbb,zer,ali,res
!       integer                             :: nrow,ncol
!       if (present(ta) .and. .not.present(A)) &
!           call quit('mat_init_as input error : ta present without A')
!       if (present(B) .and. .not.present(ta)) &
!           call quit('mat_init_as input error : B present without A')
!       if (present(tb) .and. .not.present(B)) &
!           call quit('mat_init_as input error : tb present without B')
!       if (present(zero) .and. .not.present(A)) &
!           call quit('mat_init_as input error : zero present without A')
!       if (present(alias) .and. .not.present(A)) &
!           call quit('mat_init_as input error : alias present without A')
!       if (present(alias) .and. present(ta)) &
!           call quit('mat_init_as input error : alias present with ta')
!       if (present(reset) .and. present(A)) &
!           call quit('mat_init_as input error : reset present with A')
!       taa = .false.
!       if (present(ta)) then
!           if (ta/='N'.and.ta/='T') &
!               call quit("mat_init_as input error : ta is not 'N' or 'T'")
!           taa = (ta=='T')
!       end if
!       tbb = .false.
!       if (present(tb)) then
!           if (tb/='N'.and.tb/='T') &
!               call quit("mat_init_as input error : tb is not 'N' or 'T'")
!           tbb = (tb=='T')
!       end if
!       zer = .false.
!       if (present(zero)) zer = zero
!       res = .false.
!       if (present(reset)) res = reset
!       ali = .false.
!       if (present(alias)) ali = alias
!       nrow = -1
!       if (present(A)) nrow = merge(A%ncol,A%nrow,taa)
!       ncol = -1
!       if (present(A)) ncol = merge(A%nrow,A%ncol,taa)
!       if (present(B)) ncol = merge(B%nrow,B%ncol,tbb)
!       if (zer.or.res) then
!           C%nrow = nrow
!           C%ncol = ncol
!           nullify(C%elms)
!       else if (ali) then
!           C = A !copies all fields
!       else
!           !print *,'mat_init',loc(C)
!           call mat_init_full(C,nrow,ncol)
!       end if
!   end subroutine


    function mat_isdef(A)
        type(matrix),intent(in) :: A
        logical                 :: mat_isdef
        mat_isdef = .not.(A%nrow==-1 .and. A%ncol==-1)
    end function

    function mat_iszero(A)
        type(matrix),intent(in) :: A
        logical                 :: mat_iszero
        if (A%nrow==-1 .and. A%ncol==-1) &
            call quit('mat_iszero(A) : matrix A undefined')
        mat_iszero = .not.associated(A%elms)
    end function

    function mat_same(A,B,ta,tb)
        type(matrix),         intent(in) :: A,B
        character*1, optional,intent(in) :: ta,tb
        logical                          :: mat_same
        character*1                      :: taa,tbb
        taa = 'N'; if (present(ta)) taa = ta
        tbb = 'N'; if (present(tb)) tbb = tb
        if (taa=='N'.and.tbb=='N' .or. taa=='T'.and.tbb=='T') then
            mat_same = (A%nrow==B%nrow .and. A%ncol==B%ncol)
        else if (taa=='N'.and.tbb=='T' .or. taa=='T'.and.tbb=='N') then
            mat_same = (A%nrow==B%ncol .and. A%ncol==B%nrow)
        else if (taa=='R'.and.tbb=='R') then
            mat_same = (A%nrow==B%nrow)
        else if (taa=='R'.and.tbb=='C') then
            mat_same = (A%nrow==B%ncol)
        else if (taa=='C'.and.tbb=='R') then
            mat_same = (A%ncol==B%nrow)
        else if (taa=='C'.and.tbb=='C') then
            mat_same = (A%ncol==B%ncol)
        else
            call quit('mat_same input error : ta='//taa//' tb='//tbb)
        end if
        if (.not.present(ta) .and. .not.present(tb)) &
            mat_same = mat_same .and. associated(A%elms,B%elms)
    end function

    subroutine mat_axtpby(a,X,tx,b,Y,Z)
        type(matrix),intent(in)    :: X,Y
        type(matrix),intent(inout) :: Z
        real*8,      intent(in)    :: a,b
        character*1, intent(in)    :: tx
        logical                    :: ZisY
        !print *,'mat_axtpby',loc(X),tx,loc(Y),' =>',loc(Z)
        if (tx/='N'.and.tx/='T') &
            call quit('mat_axtpby input error : tx must be N or T')
        if (b==0) then
            if (a==0) Z%elms = 0d0
            if (a/=0 .and. tx/='T') Z%elms = a*X%elms
            if (a/=0 .and. tx=='T') Z%elms = a*reshape(transpose( &
                    reshape(X%elms ,(/X%nrow,X%ncol/))),(/X%nrow*X%ncol/))
        else
            if (associated(Y%elms,Z%elms)) then
                if (a==0) Z%elms = b*Z%elms
                if (a/=0 .and. tx/='T') Z%elms = b*Z%elms + a*X%elms
                if (a/=0 .and. tx=='T') Z%elms = b*Z%elms + a*reshape(transpose( &
                                  reshape(X%elms ,(/X%nrow,X%ncol/))),(/X%nrow*X%ncol/))
            else
                if (a==0) Z%elms = b * Y%elms
                if (a/=0 .and. tx/='T') Z%elms = b*Y%elms + a*X%elms
                if (a/=0 .and. tx=='T') Z%elms = b*Y%elms + a*reshape(transpose( &
                                  reshape(X%elms ,(/X%nrow,X%ncol/))),(/X%nrow*X%ncol/))
            end if
        end if
    end subroutine

    subroutine mat_gemm(rab,A,ta,B,tb,rc,C)
        !matrix multiply: C = rc*C + rab * A^ta * B^tb
        !assumed initialied and with corresponding shapes
        type(matrix),   intent(in) :: A,B
        type(matrix),intent(inout) :: C
        real*8,         intent(in) :: rab,rc
        character*1,    intent(in) :: ta,tb
        !print *,'dgemm ta=',ta,' tb=',tb,' ',C%nrow,C%ncol,A%ncol,A%nrow,B%ncol,B%nrow
        call mat_mul(A,B,ta,tb,rab,rc,C)
    end subroutine

    function mat_dot(A,B,t) result(r)
        !matrix dot product: sum Aij (B(^T))ij
        !assumed initialied and with equal/transpose shapes
        type(matrix),        intent(in) :: A,B
        character*1,optional,intent(in) :: t
        real*8                          :: r
        integer                         :: i,j
        character*1                     :: tt
        tt = 'N'; if (present(t)) tt = t
        if (t=='N') then
            r = sum(A%elms*B%elms)
        else if (t=='T') then
            r = sum(A%elms*(/((B%elms(i+(j-1)*B%nrow),j=1,B%ncol),i=1,B%nrow)/))
        else
            call quit('mat_dotprod input error : t='//t)
        end if
    end function

    function mat_trace(A) result(r)
        !matrix sum Aii, for A square and assumed initialied
        type(matrix),intent(in) :: A
        real*8                  :: r
        integer                 :: i
        if (A%nrow/=A%ncol) call quit('mat_trace(A) : A not square')
        r = 0d0
        if (associated(A%elms)) &
            r = sum((/(A%elms(1+(i-1)*(A%nrow+1)),i=1,A%nrow)/))
    end function


    subroutine mat_print(A,unit,width,label,braces)
        !print matrix A to standard out,starting with label
        type(matrix),         intent(in) :: A
        integer,     optional,intent(in) :: unit,width
        character(*),optional,intent(in) :: label,braces
        integer                          :: u,w,d,i,j
        character*30                     :: fmt
        character*1                      :: bra,ces,csep,rsep
        !process argument unit, which defaults to standard output
        u = 6; if (present(unit)) u = unit
        !set d to the largest number of digits to be printed
        d = 1 !before the decimal point (including - signs)
        if (associated(A%elms)) d = &
           max(max(1,ceiling(log10(maxval(A%elms)))) &
              ,max(2,1+ceiling(log10(-minval(A%elms)))))
        w = 7 !process default for width
        if (present(width)) w = max(width,d+2) !max, to avoid *****
        d = w-d-1 !set d to the number of decimals to be printed
        if (present(label)) then
            fmt = '(a,ikk,a1,ikk,a1)'
            write (fmt( 5: 6),'(i2)') max(1,ceiling(log10(1d0*A%nrow)))
            write (fmt(12:13),'(i2)') max(1,ceiling(log10(1d0*A%ncol)))
            write (u,fmt(1:17)),'matrix : '//label//'(',A%nrow,',',A%ncol,')'
        end if
        fmt = '(a2,fww.dd,nnnn(a1,fww.dd),a2)'
        write (fmt(6:15),'(i2,a,i2,a,i4)'),w,'.',d,',',A%ncol-1
        fmt(21:25) = fmt(6:10)
        if (present(braces)) then
            bra = braces(1:1)
            if (len(braces)==4) then
                csep = braces(2:2)
                ces = braces(3:3)
                rsep = braces(4:4)
            else
                csep = ','; rsep = ','
                ces = merge(')',merge('}',merge(']',bra &
                           ,bra=='['),bra=='{'),bra=='(')
                if (len(braces)==2) ces = braces(2:2)
            end if
            do i=1,A%nrow
                if (associated(A%elms)) &
                    write (u,fmt),merge(bra,' ',i==1)//bra,A%elms(i) &
                                 ,(csep,A%elms(i+A%nrow*j),j=1,A%ncol-1) &
                                 ,ces//merge(ces,rsep,i==A%nrow)
                if (.not.associated(A%elms)) &
                    write (u,fmt),merge(bra,' ',i==1)//bra,0d0 &
                                 ,(csep,0d0,j=1,A%ncol-1) &
                                 ,ces//merge(ces,rsep,i==A%nrow)
            end do
        else
            fmt(2:4) = '   '   !remove braces
            fmt(17:18) = ' x'  !and separators
            fmt(27:29) = '   ' !from the format
            do i=1,A%nrow
                if (associated(A%elms)) &
                    write (u,fmt),A%elms(i:i+A%nrow*(A%ncol-1):A%nrow)
                if (.not.associated(A%elms)) &
                    write (u,fmt),(0d0,j=1,A%ncol)
            end do
        end if
        write (u,'()') !final empty line
    end subroutine


end module
