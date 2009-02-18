
! written by Andreas J. Thorvaldsen
! radovan.bast@uit.no (2009-02-18): generalized for DIRAC

! module contains slightly generalized versions of the
! operations in matrix_operations (Dalton)

! mainly intended as hidden back-end for matrix_defop
! mat_init is generalized so that another matrix can fill
! the role of nrow and ncol

module matrix_genop

#ifndef PRG_DIRAC
  use matrix_operations,                  &
      only: matrix,                       &
            mat_mul,                      &
            mat_init_full    => mat_init, &
            mat_free_nonzero => mat_free
#else
  use matrix_genop_interface
!     matrix_genop_interface provides 
!     DIRAC (only dense) versions of
!     matrix, mat_mul, mat_init_full, mat_free_nonzero
!     plus other DIRAC specifics
#endif /* PRG_DIRAC */  

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
  public mat_trace_nontmp
  public mat_print_nontmp

  private

  interface mat_init
    module procedure mat_init_full
    module procedure mat_init_from
  end interface

  integer :: nz = 1
! nz is introduced to accommodate complex and quaternion algebra for DIRAC
! it could be part of the matrix type
! for the moment it is a parameter

contains

    subroutine mat_free(A)
        type(matrix),intent(inout) :: A
        !if (associated(A%elms)) print *,'mat_free',loc(A),' elms=',loc(A%elms)
        if (associated(A%elms)) call mat_free_nonzero(A)
        A%nrow = -1; A%ncol = -1
    end subroutine


    subroutine mat_init_from(C,A,ta,B,tb,zero,alias,reset)
        type(matrix),         intent(inout) :: C
        type(matrix),optional,intent(in)    :: A,B
        character*1, optional,intent(in)    :: ta,tb
        logical,     optional,intent(in)    :: zero,alias,reset
        logical                             :: taa,tbb,zer,ali,res
        integer                             :: nrow,ncol
        if (present(ta) .and. .not.present(A)) &
            call quit('mat_init_as input error : ta present without A')
        if (present(B) .and. .not.present(ta)) &
            call quit('mat_init_as input error : B present without A')
        if (present(tb) .and. .not.present(B)) &
            call quit('mat_init_as input error : tb present without B')
        if (present(zero) .and. .not.present(A)) &
            call quit('mat_init_as input error : zero present without A')
        if (present(alias) .and. .not.present(A)) &
            call quit('mat_init_as input error : alias present without A')
        if (present(alias) .and. present(ta)) &
            call quit('mat_init_as input error : alias present with ta')
        if (present(reset) .and. present(A)) &
            call quit('mat_init_as input error : reset present with A')
        taa = .false.
        if (present(ta)) then
            if (ta/='N'.and.ta/='T') &
                call quit("mat_init_as input error : ta is not 'N' or 'T'")
            taa = (ta=='T')
        end if
        tbb = .false.
        if (present(tb)) then
            if (tb/='N'.and.tb/='T') &
                call quit("mat_init_as input error : tb is not 'N' or 'T'")
            tbb = (tb=='T')
        end if
        zer = .false.
        if (present(zero)) zer = zero
        res = .false.
        if (present(reset)) res = reset
        ali = .false.
        if (present(alias)) ali = alias
        nrow = -1
        if (present(A)) nrow = merge(A%ncol,A%nrow,taa)
        ncol = -1
        if (present(A)) ncol = merge(A%nrow,A%ncol,taa)
        if (present(B)) ncol = merge(B%nrow,B%ncol,tbb)
        if (zer.or.res) then
            C%nrow = nrow
            C%ncol = ncol
            nullify(C%elms)
        else if (ali) then
            C = A !copies all fields
        else
            !print *,'mat_init',loc(C)
            call mat_init_full(C,nrow,ncol)
        end if
    end subroutine


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


#ifndef PRG_DIRAC

    subroutine mat_axtpby(Z,a,X,tx,b,Y)
        type(matrix),intent(inout) :: Z
        type(matrix),intent(in)    :: X,Y
        real(8),     intent(in)    :: a,b
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


#else

  subroutine transpose_block(f, g, i_dim, j_dim, block, tblock)

!                            tblock = f*block^T + g*tblock

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: f
    real(8), intent(in)    :: g
    integer, intent(in)    :: i_dim
    integer, intent(in)    :: j_dim
    real(8), intent(in)    :: block(i_dim*j_dim)
    real(8), intent(inout) :: tblock(i_dim*j_dim)
!   ----------------------------------------------------------------------------

    call dscal(i_dim*j_dim, f, block, 1)

    tblock = g*tblock + reshape(transpose(reshape(block, (/i_dim, j_dim/))), (/i_dim*j_dim/))

  end subroutine

  subroutine mat_axtpby(Z, a, X, tx, b, Y)

!                       Z = a*X^tx + b*Y

!   ----------------------------------------------------------------------------
    type(matrix),  intent(inout) :: Z
    real(8),       intent(in)    :: a, b
    type(matrix),  intent(in)    :: X, Y
    character(1),  intent(in)    :: tx
!   ----------------------------------------------------------------------------
    logical                      :: Y_is_Z
    integer                      :: i_dim, j_dim, iz
    real(8)                      :: q
!   ----------------------------------------------------------------------------

    i_dim = Z%nrow
    j_dim = Z%ncol
#ifdef PRG_DIRAC
  nz = get_nz()
#endif /* PRG_DIRAC */ 

    if (tx /= 'N' .and. tx /= 'T') call quit('mat_axtpby input error: tx must be N or T')

    if (b == 0.0d0) then

      if (a == 0.0d0)                 Z%elms = 0.0d0

      if (a /= 0.0d0 .and. tx == 'N') Z%elms = a*X%elms

      if (a /= 0.0d0 .and. tx == 'T') then
        do iz = 1, nz
          if (get_ipqtoq(iz, X%pg_sym - 1) > 1) then
            q = -1.0
          else
            q =  1.0
          end if
          call transpose_block(a*q, b, i_dim, j_dim,                            &
                               X%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim), &
                               Z%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim))
        end do
    Z%ih_sym        = X%ih_sym
    Z%tr_sym        = X%tr_sym
    Z%pg_sym        = X%pg_sym
    Z%fermi_pointer = X%fermi_pointer
      end if

    else if (b /= 0.0d0 .and. associated(Y%elms, Z%elms)) then

      if (a == 0.0d0)                 Z%elms = b*Z%elms

      if (a /= 0.0d0 .and. tx == 'N') Z%elms = b*Z%elms + a*X%elms

      if (a /= 0.0d0 .and. tx == 'T') then
        do iz = 1, nz
          if (get_ipqtoq(iz, X%pg_sym - 1) > 1) then
            q = -1.0
          else
            q =  1.0
          end if
          call transpose_block(a*q, b, i_dim, j_dim,                            &
                               X%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim), &
                               Z%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim))
        end do
    Z%ih_sym        = X%ih_sym
    Z%tr_sym        = X%tr_sym
    Z%pg_sym        = X%pg_sym
    Z%fermi_pointer = X%fermi_pointer
      end if

    else

      if (a == 0.0d0)                 Z%elms = b*Y%elms

      if (a /= 0.0d0 .and. tx == 'N') Z%elms = b*Y%elms + a*X%elms

      if (a /= 0.0d0 .and. tx == 'T') then
        do iz = 1, nz
          if (get_ipqtoq(iz, X%pg_sym - 1) > 1) then
            q = -1.0
          else
            q =  1.0
          end if
          call transpose_block(a*q, 0.0d0, i_dim, j_dim,                        &
                               X%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim), &
                               Z%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim))
        end do
        Z%elms = Z%elms + b*Y%elms
    Z%ih_sym        = X%ih_sym
    Z%tr_sym        = X%tr_sym
    Z%pg_sym        = X%pg_sym
    Z%fermi_pointer = X%fermi_pointer
      end if

    end if


  end subroutine

#endif /* PRG_DIRAC */ 

    subroutine mat_gemm(rab,A,ta,B,tb,rc,C)
        !matrix multiply: C = rc*C + rab * A^ta * B^tb
        !assumed initialied and with corresponding shapes
        type(matrix),   intent(in) :: A,B
        type(matrix),intent(inout) :: C
        real(8),        intent(in) :: rab,rc
        character*1,    intent(in) :: ta,tb
        !print *,'dgemm ta=',ta,' tb=',tb,' ',C%nrow,C%ncol,A%ncol,A%nrow,B%ncol,B%nrow
        call mat_mul(A,B,ta,tb,rab,rc,C)
    end subroutine

#ifndef PRG_DIRAC

    function mat_dot(A,B,t) result(r)
        !matrix dot product: sum Aij (B(^T))ij
        !assumed initialied and with equal/transpose shapes
        type(matrix),        intent(in) :: A,B
        character*1,optional,intent(in) :: t
        real(8)                         :: r
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

#else

    function mat_dot(A,B,t) result(r)
        !matrix dot product: sum Aij (B(^T))ij
        !assumed initialied and with equal/transpose shapes
        type(matrix),        intent(in) :: A,B
        character*1,optional,intent(in) :: t
        real*8                          :: r, factor
        integer                         :: n, iz
        character*1                     :: t_set

        integer                         :: i, j, i_dim, j_dim

#ifdef PRG_DIRAC
  nz = get_nz()
#endif /* PRG_DIRAC */ 
        t_set = 'N'
        if (present(t)) t_set = t

        r = 0.0d0

        i_dim = A%nrow
        j_dim = A%ncol

        if (ieor(A%pg_sym-1, B%pg_sym-1) == 0) then

          if (t_set == 'N') then

              r = sum(A%elms*B%elms)

          else if (t_set == 'T') then

!     n = B%nrow*B%ncol
      do iz = 1, nz
        if (get_ipqtoq(iz, B%pg_sym - 1) > 1) then
          factor = -1.0
        else
          factor =  1.0
        end if

        r = r + factor*sum(A%elms((iz - 1)*i_dim*j_dim + 1:iz*i_dim*j_dim)* &
        (/((B%elms(i + (j - 1)*i_dim + (iz - 1)*i_dim*j_dim), j = 1, j_dim), i = 1, i_dim)/))

      end do

          else

            call quit('mat_dot input error: t = '//t)

          end if

        end if

    end function

#endif /* PRG_DIRAC */  

    function mat_trace_nontmp(A) result(r)
        !matrix sum Aii, for A square and assumed initialied
        type(matrix),intent(in) :: A
        real(8)                 :: r
        integer                 :: i
        if (A%nrow/=A%ncol) call quit('mat_trace(A) : A not square')
        r = 0d0
        if (associated(A%elms)) &
            r = sum((/(A%elms(1+(i-1)*(A%nrow+1)),i=1,A%nrow)/))
    end function

#ifdef PRG_DIRAC
    subroutine mat_print_nontmp(A,unit,width,label,braces)
        !print matrix A to standard out,starting with label
        type(matrix),         intent(in) :: A
        integer,     optional,intent(in) :: unit,width
        character(*),optional,intent(in) :: label,braces
        integer                          :: u,w,d,i,j
        character*30                     :: fmt
        character*1                      :: bra,ces,csep,rsep

#ifdef PRG_DIRAC
  nz = get_nz()
#endif /* PRG_DIRAC */ 
    write(*, *) 'matrix        = ', label
    write(*, *) 'ih_sym        = ', A%ih_sym
    write(*, *) 'tr_sym        = ', A%tr_sym
    write(*, *) 'pg_sym        = ', A%pg_sym
    write(*, *) 'fermi_pointer = ', A%fermi_pointer

    if (associated(A%elms)) then
      call prqmat(A%elms,                    &
                  A%nrow,                    &
                  A%ncol,                    &
                  A%nrow,                    &
                  A%ncol,                    &
                  nz,                        &
                  get_ipqtoq(1, A%pg_sym-1), &
                  6)
    else
      write(*, *) 'all elements zero'
    end if

  end subroutine
#else
    subroutine mat_print_nontmp(A,unit,width,label,braces)
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
            write (u,fmt(1:17)) 'matrix : '//label//'(',A%nrow,',',A%ncol,')'
        end if
        fmt = '(a2,fww.dd,nnnn(a1,fww.dd),a2)'
        write (fmt(6:15),'(i2,a,i2,a,i4)') w,'.',d,',',A%ncol-1
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
                    write (u,fmt) merge(bra,' ',i==1)//bra,A%elms(i) &
                                 ,(csep,A%elms(i+A%nrow*j),j=1,A%ncol-1) &
                                 ,ces//merge(ces,rsep,i==A%nrow)
                if (.not.associated(A%elms)) &
                    write (u,fmt) merge(bra,' ',i==1)//bra,0d0 &
                                 ,(csep,0d0,j=1,A%ncol-1) &
                                 ,ces//merge(ces,rsep,i==A%nrow)
            end do
        else
            fmt(2:4) = '   '   !remove braces
            fmt(17:18) = ' x'  !and separators
            fmt(27:29) = '   ' !from the format
            do i=1,A%nrow
                if (associated(A%elms)) &
                    write (u,fmt) A%elms(i:i+A%nrow*(A%ncol-1):A%nrow)
                if (.not.associated(A%elms)) &
                    write (u,fmt) (0d0,j=1,A%ncol)
            end do
        end if
        write (u,'()') !final empty line
    end subroutine
#endif /* PRG_DIRAC */  

end module
