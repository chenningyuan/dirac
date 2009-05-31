!Copyright 2009 Andreas J. Thorvaldsen
!This file is made available under the terms of the
!GNU Lesser General Public License.

!radovan: this file is generic for both DALTON and DIRAC
!         DALTON or DIRAC specific things are in matrix_genop
!         which is the backend to matrix_defop

module matrix_defop
! Defined operators for type(matrix)
! (which in linsca comes from module matrix_operations).
!
! The operators are as follows:
!
!     equals                       A = B
!     plus                         A + B
!     minus                        A - B
!     matrix product               A * B
!     scale by real(8)             r * A
!     divide by real(8)            A / r
!     (conjugate) transpose        dag(A)
!     trace                        tr(A)
!     scalar product               dot(A,B)
!     product trace                tr(A,B),
!     for zeroing and freeing      A=0, A(:)=0, A(:,:)=0
!     norm ie. sqrt(dot(A,A))      norm(A)
!     mat_to_full                  A = B(:,:)
!     mat_set_from_full            B(:,:) = A, B(:) = A
!
! To increase performance and lower memory needs, formulas are not
! evaluated by executing the operators immediately, but the operators are
! first collected into "batches" of the forms
!
!     axtpby:     Z = a * X(^T) + b * Y   (where Z may be Y)
!     gemm:       Y = a * X(^T) * Z(^T) + b * Y
!
! To create matrices, mat_init (from module matrix_operations) is replaced
! by a more general init_mat. This is to hide context-specific information:
! closed/open shell, real/complex, 1/2/4-component, full/sparse,
! and thus make it possible to write general code:
!
!     call init_mat(A,nrow,ncol)       original version
!     call init_mat(A,B)               same shape
!     call init_mat(A,B,'T')           transpose shape
!     call init_mat(A,B,'T',C,'N')     rows and columns from diff. mats.
!     call init_mat(A,B,zero=.true.)   non-allocated zero of same shape
!     call init_mat(A,reset=.true.)    reset to defaults, nullifying pointers
!     call init_mat(A,B,alias=.true.)  copy pointers, no new allocation
!
! Comparison of the shapes (and content) of matrices is done with mat_same
!
!     if (mat_same(A,B)) ..            if same shape and content (alias)
!     if (mat_same(A,B,'N')) ..        if same shapes
!     if (mat_same(A,B,'T')) ..        if transposed shapes
!     if (mat_same(A,B,'T','R')) ..    if A's rows matches B's columns
!     if (mat_same(A,B,'N','C')) ..    if A's columns matches B's columns
!
! For printing matrices, print_mat has been replaced:
!
!     call print_mat(A,label='A',colwidth=7,unit=6,braces='{')
!
! The label is printed on the line preceding the matrix. Width specifies the
! text width of each column (excluding separator). The number of digits printed
! is adjusted so that no columns become overfilled (*****). Braces can either
! be one character {, [, (, | etc., or four '{,};', where the comma is column
! separator and semicolon row separator. This is convenient for pasting into
! python, maple, mathematica, etc.
!
! To check the state of matrices: Whether defined: isdef(A),
! whether zero: iszero(A)
!
! In order to pass a "formula", such as A+3*dag(B) etc., as argument X to a subroutine,
! the subroutine must execute mat_fix_arg(X) before using X inside the subroutine.
! At the end of the subroutine mat_unfix_arg(X) should be called in order to free X.
! Neither call has any effect on a non-temporary matrix argument.
! This is neccessary since the Fortran programmer cannot tell whether the temporary
! is used inside its own subruoutine, where it should be deleted after use,
! or inside a called subroutine, in which it should be non-temporary.
!
! In order to have type(matrix) as the result of a function, mat_fix_result(X)
! should be called before returning from the function. This sets X's status as
! temporary, so that it will be deleted once used.

   use matrix_genop

   implicit none

   public matrix
   public init_mat
   public print_mat
   public isdef
   public iszero
   public assignment(=)
   public operator(+)
   public operator(-)
   public operator(*)
   public operator(/)
   public dag
   public dot
   public tr
   public norm
   public mat_fix_arg
   public mat_unfix_arg
   public mat_fix_result
   public matrix_defop_debug
   private

   interface assignment(=)
      module procedure mat_eq_mat
      module procedure mat_eq_zero
      module procedure mat_eq_zero_1D
      module procedure mat_eq_zero_2D
      module procedure mat_eq_zero_3D
   end interface

   interface isdef
      module procedure mat_isdef
   end interface

   interface iszero
      module procedure mat_iszero
   end interface

   interface operator(+)
      module procedure plus_mat
      module procedure mat_plus_mat
   end interface

   interface operator(-)
      module procedure minus_mat
      module procedure mat_minus_mat
   end interface

   interface operator(*)
      module procedure mat_times_mat
      module procedure real_times_mat
   end interface

   interface operator(/)
      module procedure mat_div_real
   end interface

   interface dag
      module procedure mat_dagger
   end interface

   interface dot
      module procedure mat_dot_product
   end interface

   interface tr
      module procedure mat_trace
      module procedure mat_product_trace
   end interface

   interface norm
      module procedure mat_norm
   end interface

   type tmp_mat
   !specification of M = rx*X^tx * Z^tz + ry*Y gemm
   !       or        M = rx*X^tx * Z^tz        gemm
   !       or        M = rx*X^tx        + ry*Y axtpby
   !       or        M = rx*X^tx               axt or by
   !M points to the matrix where result is to be stored, which is undef.
   !Flags dx,dy,dz decides whether arguments whould be deleted after use
       type(matrix),pointer :: M, X, Z, Y !associated or not- pointers
       real(8)              ::   rx,   ry !scale factors
       logical              ::   tx,tz    !transpose
       logical              ::   dx,dz,dy !delete after use
   end type

   !tmp(:) contains a 'stack' of temporary matrices
   integer,parameter         :: max_tmp = 64
   type(tmp_mat),target,save :: tmp(max_tmp)
   integer                   :: top_tmp = 0 !0 as "before index 1, thus empty"

   logical :: matrix_defop_debug = .false.

contains

   function undef(A)
      type(matrix),intent(in) :: A
      logical                 :: undef
      undef = .not.isdef(A)
   end function

   function find_tmp(A,errmsg) result(T)
   !if A is defined, return null, if A is undefined and A%elms
   !points to within tmp(1:max_tmp), return A%elms
      type(matrix),target,  intent(in) :: A
      character(*),optional,intent(in) :: errmsg
      type(tmp_mat),pointer            :: T
      if (A%nrow==-1 .and. A%ncol/=-1 .and. &
                   .not.associated(A%elms)) then
          T => tmp(A%ncol)
          T%M => A
          T%M%ncol = -1
      else if (isdef(A)) then
          nullify(T) !not found, not temporary
      else if (present(errmsg)) then
          call quit(errmsg)
      end if !otherwise, just return null
   end function

   function new_tmp(M,X,Z,Y,rx,ry,tx,tz,dx,dz,dy) result(T)
   !record the temporary, non-zero matrix M in tmp(:)
      type(matrix),         target,intent(inout) ::  M
      type(matrix),target,optional,intent(in) ::  X, Z, Y
      real(8),            optional,intent(in) :: rx,   ry
      logical,            optional,intent(in) :: tx,tz
      logical,            optional,intent(in) :: dx,dz,dy
      type(tmp_mat),pointer :: T
      integer               :: i
      if (top_tmp==max_tmp) then
          print '(a)','error: tmp(1:max_tmp) in linears/matrix_defop.f90 is full'
          print '(a,i2,a)','       increase max_tmp=',max_tmp,' and recompile'
          call quit('matrix new_tmp : tmp(1:max_tmp) is full')
      end if
      do i=1,top_tmp+1
          if (i==top_tmp) cycle !this i taken
          if (i==top_tmp+1) top_tmp = top_tmp+1
          T => tmp(i)
          if (.not.associated(T%M)) exit
      end do
      T%M => M
      nullify(T%X);   if (present(X))  T%X  => X
      nullify(T%Z);   if (present(Z))  T%Z  => Z
      nullify(T%Y);   if (present(Y))  T%Y  => Y
      T%rx = 1d0;     if (present(rx)) T%rx = rx
      T%ry = merge(1d0,0d0,present(Y))
                      if (present(ry)) T%ry = ry
      T%tx = .false.; if (present(tx)) T%tx = tx
      T%tz = .false.; if (present(tz)) T%tz = tz
      T%dx = .false.; if (present(dx)) T%dx = dx
      T%dz = .false.; if (present(dz)) T%dz = dz
      T%dy = .false.; if (present(dy)) T%dy = dy
      !hack: point M%ncol to T,
      if (.not.isdef(M)) then
          M%ncol = i      !so we can find T when we
          nullify(M%elms) !eventually want to use M
      end if
   end function

   subroutine del_tmp(T)
   !remove tmp_mat T from within tmp(:), deallocate t%M if needed
      type(tmp_mat),intent(inout) :: T
      nullify(T%M)
      !update last occupied index
      do while (top_tmp/=0)
          if (associated(tmp(top_tmp)%M)) exit
          top_tmp = top_tmp - 1
      end do
   end subroutine

   subroutine zero_tmp(T)
   !initialize T%M to a zero matrix of the correct shape,
   !then delete (or not) T%X,T%Z,T%Y, then del_tmp(T)
      type(tmp_mat),intent(inout) :: T
      !if (matrix_defop_debug) print *,'zero_tmp',loc(T%M),T%M%nrow,T%M%ncol
      if (T%M%nrow==-1 .and. T%M%ncol/=-1) T%M%ncol = -1
      if (isdef(T%M)) then
          if (.not.iszero(T%M)) call mat_free(T%M)
      end if
      if (associated(T%Z)) then
          call init_mat(T%M,T%X,merge('T','N',T%tx) &
                           ,T%Z,merge('T','N',T%tz),zero=.true.)
      else
          call init_mat(T%M,T%X,merge('T','N',T%tx),zero=.true.)
      end if
      if (T%dx) call mat_free(T%X)
      if (T%dz) call mat_free(T%Z)
      if (T%dy) call mat_free(T%Y)
      call del_tmp(T)
   end subroutine


   subroutine eval_tmp(T)
   !evaluate the expression M = a*X*Z + b*Y, as specified in the fields of T,
   !by executing either an axtpby or a gemm operation. Then remove T from tmp(:)
      type(tmp_mat),intent(inout) :: T
      logical                     :: zeroM,useM,ltmp
      real(8)                     :: rtmp
      type(matrix),target         :: reM
      !if (matrix_defop_debug) print *,'eval_tmp',loc(T%M),' =',loc(T%X) &
      !                                     ,' *',loc(T%Z),' +',loc(T%Y)
      !hack: remove the T-reference from M%ncol
      if (T%M%nrow==-1 .and. T%M%ncol/=-1 &
          .and..not.associated(T%M%elms)) T%M%ncol = -1
      zeroM = .not.isdef(T%M)          !two steps, to prevent
      if (.not.zeroM) zeroM = iszero(T%M) !iszero(undef) error
      !if M = M^N + Y or X^N + Y with only X temporary, swap X and Y
      if (.not.associated(T%Z) .and. .not.T%dy .and. .not.associated(T%Y,T%M) &
          .and. .not.T%tx .and. (T%dx.or.associated(T%X,T%M))) then
          ltmp = T%dx; T%dx = T%dy; T%dy = ltmp
          rtmp = T%rx; T%rx = T%ry; T%ry = rtmp
          T%Z => T%X;  T%X => T%Y;  T%Y => T%Z; nullify(T%Z)
      end if
      !if M = M^T + Y, M*Z + Y or X*M + Y, it must not be overwritten
      !during the operation, so move M=X/Z to reM and clear M
      if (associated(T%Z,T%M)) then
          call init_mat(reM,T%M,alias=.true.)
          T%Z => reM
          call init_mat(T%M,reset=.true.)
          zeroM = .true.
          if (associated(T%X,T%M)) T%X => reM
          if (associated(T%Y,T%M)) T%Y => reM
      else if (associated(T%X,T%M).and.(T%tx.or.associated(T%Z))) then
          call init_mat(reM,T%M,alias=.true.)
          T%X => reM
          call init_mat(T%M,reset=.true.)
          zeroM = .true.
          if (associated(T%Y,T%M)) T%Y => reM
      end if
      !determine whether to use M's current allocation, or Y's
      if (zeroM) then
          useM = .false. !nothing to free
      else if (associated(T%Y).and.T%dy) then
          useM = .false.
      else if (associated(T%Y,T%M)) then
          useM = .true. !obviously correct shape
      else if (associated(T%Y).and..not.T%dy) then
          useM = mat_same(T%M,T%Y,'N')
      else if (.not.associated(T%Z)) then
          useM = mat_same(T%M,T%X,merge('T','N',T%tx))
      else
          useM = ( mat_same(T%M,T%X,'R',merge('C','R',T%tx)) &
              .and.mat_same(T%M,T%Z,'C',merge('R','C',T%tz)) )
      end if
      if (.not.zeroM .and. .not.useM) call mat_free(T%M)
      !determine whether A should be reallocated, or simply duplicate Y.
      if (useM) then
          !no allocation needed
      else if (associated(T%Y).and.T%dy) then
          !If Y is temporary, no allocation, just transfer Y to M
          call init_mat(T%M,T%Y,alias=.true.)
          call init_mat(T%Y,reset=.true.)
          T%dy = .false.
          T%Y => T%M
      else if (associated(T%Y,T%M)) then
          !if Y=M, no allocation, the result is added to Y
      else if (associated(T%Y).and.associated(T%Z)) then
          !if M=a*X*Z+b*Y, allocate M, and do M=b*Y first
          if (.not.useM) call init_mat(T%M,T%Y,'N')
          call mat_axtpby(T%M, T%ry,T%Y,'N',0d0,T%M)
          T%Y => T%M
          T%ry = 1d0
      else if (associated(T%Z)) then
          !if M=a*X*Z, allocate M
          call init_mat(T%M,T%X,merge('T','N',T%tx) &
                           ,T%Z,merge('T','N',T%tz))
      else !if M=a*X+b*Y or M=a*X, allocate M
          call init_mat(T%M,T%X,merge('T','N',T%tx))
      end if
      !execute the requested operation, axtpby or gemm
      if (associated(T%Z)) then
          call mat_gemm(T%rx,T%X,merge('T','N',T%tx) &
                            ,T%Z,merge('T','N',T%tz),T%ry,T%M) !Y is M
      else if (associated(T%X).or.T%ry/=1d0 &
              .or..not.associated(T%Y,T%M)) then
          call mat_axtpby(T%M, T%rx,T%X,merge('T','N',T%tx),T%ry,T%Y)
      end if
      !delete X and/or Z, as specified. Y is never deleted
      if (associated(T%X,reM).or.associated(T%Z,reM)) &
          call mat_free(reM)
      if (T%dx) call mat_free(T%X)
      if (T%dz) call mat_free(T%Z)
      call del_tmp(T)
   end subroutine


   subroutine mat_eq_mat(A,B)
   !denne skal altid tømme tmp(:)
      type(matrix),target,intent(inout) :: A
      type(matrix),       intent(in)    :: B
      type(tmp_mat),pointer             :: T
      logical                           :: defB
      !if (matrix_defop_debug) print *,'mat_eq_mat',loc(A),' =' &
      !                ,loc(B),B%nrow,B%ncol,associated(B%elms)
      T => find_tmp(B,'matrix A = B : B undefined')
      if (associated(T)) then
          T%M => A
          call eval_tmp(T)
      else if (mat_same(A,B)) then
          !nothing
      else if (iszero(B)) then
          T => new_tmp(A,X=B)
          call zero_tmp(T)
      else
          T => new_tmp(A,X=B)
          call eval_tmp(T)
      end if
#ifdef PRG_DIRAC
    A%ih_sym        =  B%ih_sym
    A%tr_sym        =  B%tr_sym
    A%pg_sym        =  B%pg_sym
    A%fermi_pointer =  B%fermi_pointer
#endif /* PRG_DIRAC */

   end subroutine

   subroutine gen_mat_plus_mat(A,plus,B,C)
      type(matrix),target,intent(in)  :: A,B
      type(matrix),target             :: C
      type(tmp_mat),pointer           :: TA,TB,TC
      logical                         :: zeroA,zeroB,swap,ltmp,plus
      !if (matrix_defop_debug) print *,'mat_plus_minus_mat',loc(A),loc(B) &
      !                ,' => ',loc(C),' A%nrow=',A%nrow,' A%ncol=',A%ncol &
      !                ,' B%nrow=',B%nrow,' B%ncol=',B%ncol,associated(B%elms)
      TA => find_tmp(A,'matrix A + B : A undefined')
      TB => find_tmp(B,'matrix A + B : B undefined')
      zeroA = .false.; zeroB = .false.
      if (.not.associated(TA)) zeroA = iszero(A)
      if (.not.associated(TB)) zeroB = iszero(B)
      !decide whether to swap A and B
      if (associated(TA).and.associated(TB)) then
          swap = ( .not.associated(TB%Y) .and. (associated(TB%Z) &
                  .or..not.TA%tx.or.associated(TA%Z).or.associated(TA%Y)))
      else if (.not.associated(TA).and.associated(TB)) then
          swap = (zeroA .or. .not.associated(TB%Y))
      else
          swap = (zeroA .and. .not.zeroB)
      end if
      !make the swap, and create a new tmp_mat C = A + B
#ifndef PRG_DIRAC
!fixme: for some mysterious reason this 
!       does not work for me with DIRAC
       if (swap) then
           ltmp = zeroA; zeroA = zeroB; zeroB = ltmp
           TC  => TA;    TA   => TB;    TB   => TC
           TC => new_tmp(C,X=B,Y=A,rx=merge(1d0,-1d0,plus) &
                        ,dx=associated(TA),dy=associated(TB))
       else
#endif /* PRG_DIRAC */
            TC => new_tmp(C,X=A,Y=B,ry=merge(1d0,-1d0,plus) &
                         ,dx=associated(TA),dy=associated(TB))
#ifndef PRG_DIRAC
        end if
#endif /* PRG_DIRAC */
        !if A is undefined and X or X*Z (simple) promote into C%X,C%Z
        !if A is undefined and X*Z + Y (complicated), evaluate it
        if (associated(TA)) then
#ifndef PRG_DIRAC
          if (.not.associated(TA%Y)) then
              TC%rx = TC%rx*TA%rx
              TC%tx = TA%tx; TC%tz = TA%tz
              TC%dx = TA%dx; TC%dz = TA%dz
              TC%X => TA%X;  TC%Z => TA%Z
              call del_tmp(TA)
          else
              call eval_tmp(TA)
          end if
#else
!fixme: for some mysterious reason this 
!       does not work for me with DIRAC
!       for the moment i simply evaluate always
                call eval_tmp(TA)
#endif /* PRG_DIRAC */
       end if
       !if B is undefined and X^N (simple) promote into C%Y,
       !otherwise evaluate B
       if (associated(TB)) then
#ifndef PRG_DIRAC
          if (.not.TB%tx .and. .not.associated(TB%Z) &
                         .and. .not.associated(TB%Y)) then
              TC%ry = TC%ry*TB%rx
              TC%dy = TB%dx
              TC%Y => TB%X
              call del_tmp(TB)
          else
              call eval_tmp(TB)
          end if
#else
!fixme: for some mysterious reason this 
!       does not work for me with DIRAC
!       for the moment i simply evaluate always
                call eval_tmp(TB)
#endif /* PRG_DIRAC */
      end if
      !verify that shapes do indeed match
      if (associated(TC%Z)) then
          if (.not.mat_same(TC%X,TC%Y,'R',merge('C','R',TC%tx)) .or. &
              .not.mat_same(TC%Z,TC%Y,'C',merge('R','C',TC%tz))) &
              call quit('matrix A + B : different shapes')
      else
          if (.not.mat_same(TC%X,TC%Y,merge('T','N',TC%tx))) &
              call quit('matrix A + B : different shapes')
      end if
      !if either of A and/or B are zero, remove from TC
      if (zeroA.and.zeroB) then
          call zero_tmp(TC)
      else if (zeroA) then
          TC%X => TC%Y;  nullify(TC%Y)
          TC%rx = TC%ry; TC%ry = 0d0
          TC%tx = .false.;
          TC%dx = TC%dy; TC%dy = .false.
      else if (zeroB) then
          nullify(TC%Y)
          TC%ry = 0d0
          TC%dy = .false.
      end if

#ifdef PRG_DIRAC
    C%ih_sym        =  B%ih_sym
    C%tr_sym        =  B%tr_sym
    C%pg_sym        =  B%pg_sym
    C%fermi_pointer =  B%fermi_pointer
#endif /* PRG_DIRAC */

   end subroutine

   function mat_plus_mat(A,B) result(C)
   !verify that shapes match, check for zeros
      type(matrix),target,intent(in) :: A,B
      type(matrix),target            :: C
      call gen_mat_plus_mat(A,.true.,B,C)
   end function


   function mat_minus_mat(A,B) result(C)
   !verify that shapes match, check for zeros
      type(matrix),target,intent(in) :: A,B
      type(matrix),target            :: C
      call gen_mat_plus_mat(A,.false.,B,C)
   end function


   function mat_times_mat(A,B) result(C)
   !verify that shapes match, check for zeroes
      type(matrix),target,intent(in) :: A,B
      type(matrix),target            :: C
      type(tmp_mat),pointer          :: TA,TB,TC
      logical                        :: zeroA,zeroB
      !if (matrix_defop_debug) print *,'mat_times_mat',loc(A) &
      !                            ,' *',loc(B),' => ',loc(C)
      TA => find_tmp(A,'matrix A * B : A undefined')
      TB => find_tmp(B,'matrix A * B : B undefined')
      zeroA = .false.; zeroB = .false.
      if (.not.associated(TA)) zeroA = iszero(A)
      if (.not.associated(TB)) zeroB = iszero(B)
      TC => new_tmp(C,X=A,dx=associated(TA),Z=B,dz=associated(TB))
      if (associated(TA)) then
          if (.not.associated(TA%Z) .and. .not.associated(TA%Y)) then
              TC%rx = TC%rx*TA%rx
              TC%tx = TA%tx
              TC%dx = TA%dx
              TC%X => TA%X
              call del_tmp(TA)
          else
              call eval_tmp(TA)
          end if
      end if
      if (associated(TB)) then
          if (.not.associated(TB%Z) .and. .not.associated(TB%Y)) then
              TC%rx = TC%rx*TB%rx
              TC%tz = TB%tx
              TC%dz = TB%dx
              TC%Z => TB%X
              call del_tmp(TB)
          else
              call eval_tmp(TB)
          end if
      end if
      !finally, verify that rows
      if (.not.mat_same(TC%X,TC%Z,merge('R','C',TC%tx),merge('C','R',TC%tz))) &
          call quit('matrix A * B : incompatible shapes')
      !if either of X or Z are zero, zero C
      if (zeroA.or.zeroB) call zero_tmp(TC)
   end function


   subroutine gen_real_times_mat(r,A,B)
      real(8),            intent(in) :: r
      type(matrix),target,intent(in) :: A
      type(matrix),target            :: B
      type(tmp_mat),pointer          :: T,TT
      !if (matrix_defop_debug) print *,'real_times_mat',r,' *',loc(A),' =>',loc(B)
      T => find_tmp(A,'matrix r * A : A undefined')
      if (associated(T)) then
          TT => new_tmp(B,X=T%X,Z=T%Z,Y=T%Y,rx=r*T%rx,ry=r*T%ry &
                       ,tx=T%tx,tz=T%tz,dx=T%dx,dz=T%dz,dy=T%dy)
          call del_tmp(T)
          if (r==0d0) call zero_tmp(TT)
      else
          T => new_tmp(B,X=A,rx=r)
          if (r==0d0.or.iszero(A)) call zero_tmp(T)
      end if
#ifdef PRG_DIRAC
    B%ih_sym        =  A%ih_sym
    B%tr_sym        =  A%tr_sym
    B%pg_sym        =  A%pg_sym
    B%fermi_pointer =  A%fermi_pointer
#endif /* PRG_DIRAC */

  end subroutine

  function real_times_mat(r,A) result(B)
     real(8),            intent(in) :: r
     type(matrix),target,intent(in) :: A
     type(matrix),target            :: B
     call gen_real_times_mat(r,A,B)
  end function

  function plus_mat(A) result(B)
     type(matrix),target,intent(in) :: A
     type(matrix),target            :: B
     call gen_real_times_mat(1d0,A,B)
  end function

  function minus_mat(A) result(B)
     type(matrix),target,intent(in) :: A
     type(matrix),target            :: B
     call gen_real_times_mat(-1d0,A,B)
  end function

  function mat_div_real(A,r) result(B)
     type(matrix),target,intent(in) :: A
     real(8),            intent(in) :: r
     type(matrix),target            :: B
     if (r==0) call quit('matrix A / r : division by zero')
     call gen_real_times_mat(1/r,A,B)
  end function

  function mat_dagger(A) result(B)
  !transpose the matrix A. If A is temporary and A=a*X or A=a*X*Z,
  !only swap X and Z, and slip-swap tx and tz.
  !Otherwise, evaluate A (if tmp), and create a new temporary for B
     type(matrix),target,intent(in) :: A
     type(matrix),target            :: B
     type(tmp_mat),pointer          :: T,TT
     logical                        :: ltmp
     !if (matrix_defop_debug) print *,'mat_dagger',loc(A),' => ',loc(B)
     T => find_tmp(A,'matrix dag(A) : A undefined')
     if (.not.associated(T)) then
         T => new_tmp(B,X=A,tx=.true.)
         if (iszero(A)) call zero_tmp(T)
     else if (associated(T%Y)) then
         call eval_tmp(T)
         T => new_tmp(B,X=A,tx=.true.,dx=.true.)
     else if (associated(T%Z)) then
         TT => new_tmp(B,X=T%Z,Z=T%X,rx=T%rx,tx=.not.T%tz &
                      ,tz=.not.T%tx,dx=T%dz,dz=T%dx)
         call del_tmp(T)
     else
         TT => new_tmp(B,X=T%X,rx=T%rx,tx=.not.T%tx,dx=T%dx)
         call del_tmp(T)
     end if
#ifdef PRG_DIRAC
    B%ih_sym        =  A%ih_sym
    B%tr_sym        =  A%tr_sym
    B%pg_sym        =  A%pg_sym
    B%fermi_pointer =  A%fermi_pointer
#endif /* PRG_DIRAC */

   end function


   function gen_mat_dot(A,B,trps) result(r)
   !back for mat_dot_product and mat_product_trace (dot(A,B) and tr(A,B), resp.)
      type(matrix),target,intent(inout) :: A,B
      real(8)                           :: r
      type(tmp_mat),pointer             :: TA,TB
      logical                           :: zeroA,zeroB,trps,delA,delB
      !if (matrix_defop_debug) print *,'mat_product_trace/mat_dot_product',loc(A),loc(B)
      TA => find_tmp(A,'matrix dot/tr(A,B) : A undefined')
      TB => find_tmp(B,'matrix dot/tr(A,B) : B undefined')
      zeroA = .false.; zeroB = .false.
      if (.not.associated(TA)) zeroA = iszero(A)
      if (.not.associated(TB)) zeroB = iszero(B)
      r = merge(0d0,1d0,zeroA.or.zeroB)
      delA = associated(TA); delB = associated(TB)
      if (associated(TA)) then
          if (.not.associated(TA%Z).and..not.associated(TA%Y)) then
              r = r*TA%rx
              trps = (trps.neqv.TA%tx)
              call init_mat(A,TA%X,alias=.true.)
              call init_mat(TA%X,reset=.true.)
              delA = TA%dx
              call del_tmp(TA)
          else
              call eval_tmp(TA)
          end if
      end if
      if (associated(TB)) then
          if (.not.associated(TB%Z).and..not.associated(TB%Y)) then
              r = r*TB%rx
              trps = (trps.neqv.TB%tx)
              call init_mat(B,TB%X,alias=.true.)
              call init_mat(TB%X,reset=.true.)
              delB = TB%dx
              call del_tmp(TB)
          else
              call eval_tmp(TB)
          end if
      end if
      !verify that shapes match
      if (.not.mat_same(A,B,merge('T','N',trps))) &
          call quit('matrix dot/tr(A,B) : different shapes')
      !calculate only if prefactor nonzero
      if (r/=0d0) r = r*mat_dot(A,B,merge('T','N',trps)) * 2 !due to alpha and beta blocks
      !release any temporary operands
      if (delA) call mat_free(A)
      if (delB) call mat_free(B)
   end function


   function mat_dot_product(A,B) result(r)
      type(matrix),target,intent(inout) :: A,B
      real(8)                           :: r
      !print *,'mat_dot_product',loc(A),loc(B)
      r = gen_mat_dot(A,B,(.false.))
   end function

   function mat_product_trace(A,B) result(r)
      type(matrix),target,intent(inout) :: A,B
      real(8)                           :: r
      !print *,'mat_product_trace',loc(A),loc(B)
      r = gen_mat_dot(A,B,(.true.))
   end function

   function mat_trace(A) result(r)
      type(matrix),target,intent(inout) :: A
      real(8)                           :: r
      type(tmp_mat),pointer             :: T
      T => find_tmp(A,'matrix r * A : A undefined')
      if (associated(T)) then
          if (associated(T%Z)) then
              if (.not.mat_same(T%X,T%Z,merge('N','T',T%tx.neqv.T%tz) &
                                       ,merge('C','R',T%tx))) &
                  call quit('matrix tr(A) : A not square')
              r = T%rx * mat_dot(T%X,T%Z,merge('T','N',T%tx.neqv.T%tz))
          else
              if (.not.mat_same(T%X,T%X)) &
                  call quit('matrix tr(A) : A not square')
              r = T%rx * mat_trace_nontmp(T%X)
          end if
          if (associated(T%Y)) &
              r = r + T%ry * mat_trace_nontmp(T%Y)
          call zero_tmp(T) !this leaves A as zero
          call init_mat(A,reset=.true.) !undefine A
      else if (.not.iszero(A)) then
          r = mat_trace_nontmp(A)
      else
          r = 0d0
      end if
   end function

   function mat_norm(A) result(r)
      type(matrix),target,intent(inout) :: A
      real(8)                           :: r
      !print *,'mat_norm',loc(A)
      r = sqrt(gen_mat_dot(A,A,(.false.)))
   end function

   subroutine mat_eq_zero(A,zero)
      type(matrix),intent(inout) :: A
      integer,intent(in)         :: zero
      type(matrix)               :: B
      !if (matrix_defop_debug) print *,'mat_eq_zero A=',loc(A),' elms=',loc(A%elms)
      if (zero/=0) call quit('matrix A = z: z must be 0')
      if (.not.isdef(A)) call quit('A = 0: A undefined')
      call init_mat(B,A,alias=.true.)
      call init_mat(A,B,zero=.true.)
      call mat_free(B)
   end subroutine

   subroutine mat_eq_zero_1D(A,zero)
      type(matrix),intent(inout) :: A(:)
      integer,intent(in)         :: zero
      integer                    :: i
      do i=1,size(A)
          A(i) = 0
      end do
   end subroutine

   subroutine mat_eq_zero_2D(A,zero)
      type(matrix),intent(inout) :: A(:,:)
      integer,intent(in)         :: zero
      integer                    :: i,j
      do j=1,size(A,2)
          do i=1,size(A,1)
              A(i,j) = 0
          end do
      end do
   end subroutine

   subroutine mat_eq_zero_3D(A,zero)
      type(matrix),intent(inout) :: A(:,:,:)
      integer,intent(in)         :: zero
      integer                    :: i,j,k
      do k=1,size(A,3)
          do j=1,size(A,2)
              do i=1,size(A,1)
                  A(i,j,k) = 0
              end do
          end do
      end do
   end subroutine

   subroutine print_mat(A,unit,width,label,braces)
      type(matrix)                     :: A
      integer,     optional,intent(in) :: unit,width
      character(*),optional,intent(in) :: label,braces
      type(tmp_mat),pointer            :: T
      logical                          :: delA!,defA
!       defA = isdef(A)
!       nullify(T); if (.not.defA) T => find_tmp(A)
!       if (.not.defA .and. .not.associated(T)) &
!           call quit('print_mat(A) : A undefined')
      T => find_tmp(A,'print_mat(A) : A undefined')
      delA = associated(T)
      if (associated(T)) call eval_tmp(T)
      call print_mat_nontmp(A,unit,width,label,braces)
      if (delA) call mat_free(A)
   end subroutine

   subroutine mat_fix_arg(A)
      type(matrix),intent(inout) :: A
      type(tmp_mat),pointer      :: T
      !no action on non-temporary arguments
      if (isdef(A)) return
      T => find_tmp(A,'mat_fix_arg(A) : A undefined')
      call eval_tmp(T)
      T => new_tmp(A) !hide A in tmp(:)
   end subroutine

   subroutine mat_unfix_arg(A)
      type(matrix),intent(inout) :: A
      type(tmp_mat),pointer      :: T
      T => find_tmp(A) !look for A in tmp
      if (.not.associated(T)) return !it wasn't previously fixed
      call zero_tmp(T)
      call init_mat(A,reset=.true.)
   end subroutine

   subroutine mat_fix_result(A)
      type(matrix),intent(inout) :: A
      type(tmp_mat),pointer      :: T
      type(matrix),save,target   :: fixA
      if (.not.isdef(A)) &
          call quit('mat_fix_result(A) : A undefined')
      if (isdef(fixA)) &
          call quit('mat_fix_result(A) : another fixed result already exists')
      call init_mat(fixA,A,alias=.true.) !diplicate A in fixA
      call init_mat(A,reset=.true.)      !reset/clear/undefine A
      T => new_tmp(A,X=fixA,dx=.true.)   !temporary A = 1*fixA, fixA to be deleted
   end subroutine

end module
