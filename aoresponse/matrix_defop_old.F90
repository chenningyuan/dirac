module matrix_defop

  use matrix_genop

  implicit none

  public tr
  public trab
  public trps
  public dag
  public norm
  public dot
  public iszero
  public isdef
  public mat_print
  public undef
  public matrixstats

  public assignment(=)
  public operator(+)
  public operator(-)
  public operator(*)
  public operator(/)
  public operator(**)

  private

  interface assignment(=)
    module procedure matrix_eq_matrix
    module procedure matrix_eq_zero
        module procedure mat_eq_zero_1D
        module procedure mat_eq_zero_2D
        module procedure mat_eq_zero_3D
  end interface

  interface operator(+)
    module procedure matrix_plus_matrix
  end interface

  interface operator(-)
    module procedure matrix_minus_matrix
    module procedure minus_matrix
  end interface

  interface operator(*)
    module procedure real_times_matrix
    module procedure matrix_times_matrix
  end interface

  interface operator(/)
    module procedure matrix_div_real
  end interface

  interface dot
    module procedure matrix_dot_matrix
  end interface

  interface operator(**)
    module procedure matrix_dot_matrix
  end interface

  interface tr
    module procedure matrix_product_trace
  end interface

  interface trab
    module procedure matrix_product_trace
  end interface

  interface dag
    module procedure dag_matrix
  end interface

  interface trps
    module procedure dag_matrix
  end interface

  interface isdef
    module procedure mat_isdef
  end interface

  interface norm
    module procedure matrix_norm
  end interface

  interface iszero
    module procedure matrix_is_zero
  end interface

  interface mat_print
    module procedure print_matrix
  end interface

  interface undef
    module procedure matrix_is_undefined
  end interface

  type refp
    type(matrix), pointer :: m !pointer to referring matrix
    integer,      pointer :: s !stack top at creation (if temporary)
    integer               :: i !index of referred data (in dat(:) table)
  end type

  type datp
    real*8, pointer :: d(:) !pointer to data
    integer         :: c !reference count
  end type

  integer                 :: nref = 0, ndat = 0, maxdat = 0, stkdir
  type(refp), allocatable :: ref(:)
  type(datp), allocatable :: dat(:)

contains

  function matrix_is_undefined(A) result(r)

    type(matrix), intent(in) :: A
    logical                  :: r

    r = (A%nrow == -1)

  end function

  subroutine allocif()
  !if needed allocate/extend the tables refs,data. Set stkdir
      type(refp):: reft(nref) !temporary for ref
      type(datp):: datt(maxdat) !temporary for dat
      integer:: maxref,i
      maxref=0; if (maxdat/=0) maxref = size(ref)
      if (maxref>nref+1.and.maxdat>ndat+1) return !no reallocation needed
      if (maxdat==0) call getdir()
      if (maxref<=nref+1) then
          if (maxdat/=0) reft(1:nref) = ref(1:nref)
          if (maxdat/=0) deallocate(ref)
          allocate(ref(max(64,2*maxref)))
          ref(1:nref) = reft(1:nref)
      endif
      if (maxdat<=ndat+1) then
          if (maxdat/=0) datt(1:maxdat) = dat(1:maxdat)
          if (maxdat/=0) deallocate(dat)
          allocate(dat(max(8,2*maxdat)))
          dat(1:maxdat) = datt(1:maxdat)
          dat(maxdat+1:)%c = 0 !vacant slots
          maxdat = max(8,2*maxdat)
      endif
  contains
      subroutine getdir()
          integer:: j
          stkdir = merge(-1,1,loc(j)<loc(i))
      end subroutine
  end subroutine

  subroutine delabove(k)
  !delete all type(Matrix) located higher than k on the stack
      integer:: k,i,j
      call allocif()
      do i=nref,1,-1
          if ((loc(ref(i)%m)-loc(k))*stkdir<0) exit
          j = ref(i)%i
          dat(j)%c = dat(j)%c-1
          if (dat(j)%c==0) then
              deallocate(dat(j)%d)
              ndat = ndat-1
          endif
      enddo
      nref = i
  end subroutine

  subroutine insref(A,s)
  !Register the reference A => A%elms in ref(:) and dat(:)
      type(Matrix),intent(inout),target:: A
      integer,intent(in),optional,target:: s !stack top
      integer:: i,j,k; logical:: found,zero,same

      !see if any reference from A is recorded
      do i=nref,0,-1
          if (i==0) exit
          if ((loc(ref(i)%m)-loc(A))*stkdir<0) exit
      enddo; j=i+1; k=ref(j)%i
      found  = (j/=nref+1.and.associated(ref(j)%m,A))
      zero = (.not.associated(A%elms).or.abs(loc(A%elms)-loc(i))<2**28)
      same = .false.; if (found) same = associated(dat(k)%d,A%elms)
      !if A is recorded, but now is zero or A%elms is on the stack, remove the previous record
      if (found.and.zero) then
           ref(j:nref-1) = ref(j+1:nref)
           nref = nref-1
      endif
      !if A was recorded as referring to some other data, remove that
      if (found.and..not.same) then
          dat(k)%c = dat(k)%c-1
          if (dat(k)%c==0) deallocate(dat(k)%d)
          if (dat(k)%c==0) ndat = ndat-1
      endif
      !see if A%elms is already recorded as some dat(k)%d
      do i=size(dat),1,-1
          if (dat(i)%c==0) then
              k = i; cycle
          elseif (associated(A%elms,dat(i)%d)) then
              k = i; exit
          endif
      enddo
      !if A is not already recorded, make space for it
      if (.not.found.and..not.zero) then
          ref(j+1:nref+1) = ref(j:nref)
          nref = nref+1
      endif
      !insert A
      if (.not.zero) then
          ref(j)%m => A
          ref(j)%i = k
          nullify(ref(j)%s)
          if (present(s)) ref(j)%s => s
          dat(k)%d => A%elms
          if (dat(k)%c==0) ndat = ndat+1
          if (.not.same) dat(k)%c = dat(k)%c+1
      endif
  end subroutine

  subroutine deltmparg(A, B,r)
  !delete matrix A if found in ref(:), and stacktop at
  !creation was higher than B or r
      type(Matrix),intent(in),target:: A
      type(Matrix),intent(in),optional:: B
      real*8,intent(in),optional:: r
      integer:: i,j
      if (present(B).eqv.present(r)) &
          call QUIT('TypeMatrix/deltmparg: B OR r should be present')
      !search for A
      do i=nref,1,-1
          if ((loc(ref(i)%m)-loc(A))*stkdir<=0) exit
      enddo; j=ref(i)%i
      if (.not.associated(ref(i)%m,A)) return
      !no need: if (.not.associated(dat(j)%d,A%elms)) return
      if (.not.associated(ref(i)%s)) return !if A is non-temp
      if (present(B).and.(loc(B)-loc(ref(i)%s))*stkdir>0) return
      if (present(r).and.(loc(r)-loc(ref(i)%s))*stkdir>0) return
      !delete A
      !no need: A%nrow = -1; nullify(A%elms)
      ref(i:nref-1) = ref(i+1:nref)
      nref = nref-1
      dat(j)%c = dat(j)%c-1
      if (dat(j)%c==0) then
          deallocate(dat(j)%d)
          ndat = ndat-1
      endif
  end subroutine

  subroutine MatrixStats(label)
      character(*),intent(in),optional:: label
      integer:: i,j; character*40:: fmt
      if (present(label)) print '(a)',label
      if (stkdir==0) call allocif() !to get ref,dat allocated
      print '("ndat="i2" size(dat)="i3" nref="i2" size(ref)="i4)', &
              ndat,size(dat),nref,size(ref)
      i = maxval(ref(1:nref)%i)
      write (fmt,'("(a4,"i3"i13/a4,"i3"i13/a4,"i3"i2,a6,"i3"i2/)")') nref,i,nref,i
      print fmt,'ref=',(loc(ref(j)%m),j=1,nref),'dat=',(loc(dat(j)%d),j=1,i) &
               ,'idx=',ref(1:nref)%i,'  cnt=',dat(1:i)%c
      call delabove(i)
  end subroutine


! matrix operations which use ref
! ===============================

  subroutine matrix_eq_matrix(A, B)

    type(matrix), intent(inout) :: A
    type(matrix), intent(in)    :: B
    integer                     :: i

    call delabove(i)

!   if (.not. B%is_defined) call quit('A = B: B undefined')

    call matrix_eq_matrix_backend(A, B)

    call insref(A)
    call deltmparg(B, A)

  end subroutine

  function dag_matrix(A) result(B)

    type(matrix), intent(in) :: A
    type(matrix)             :: B
    integer                  :: i

    call delabove(i)

!   if (.not. A%is_defined) call quit('B = dag(A): A undefined')

    call dag_matrix_backend(A, B)

    call insref(B, i)
    call deltmparg(A, B)

  end function

  subroutine matrix_eq_zero(A, zero)

    type(matrix), intent(inout) :: A
    integer,      intent(in)    :: zero
    integer                     :: i

    call delabove(i)

!   if (.not. A%is_defined) call quit('A = 0: A undefined')

    if (zero /= 0) call quit('A = z: z must be 0')

    call insref(A)  !store data on A, kick out any previous data
    nullify(A%elms) !this means A is zero
    call insref(A)  !if no other ref, this will deallocate A's data

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



  function matrix_plus_matrix(A, B) result(C)

    type(matrix), intent(in) :: A, B
    type(matrix)             :: C
    integer                  :: i

    call delabove(i)

!   if (.not. A%is_defined) call quit('C = A + B: A undefined')
!   if (.not. B%is_defined) call quit('C = A + B: B undefined')

    C =  matrix_plus_matrix_backend(A, B)

    call insref(C, i)
    call deltmparg(A, C)
    call deltmparg(B, C)

  end function

  function real_times_matrix(r, A) result(B)

    type(matrix), intent(in) :: A
    real*8,       intent(in) :: r
    type(matrix)             :: B
    integer                  :: i

    call delabove(i)

!   if (.not. A%is_defined) call quit('B = r*A: A undefined')

    B = real_times_matrix_backend(r, A)

    call insref(B, i)
    call deltmparg(A, B)

  end function

  function matrix_times_matrix(A, B) result(C)

    type(matrix), intent(in) :: A, B
    type(matrix)             :: C
    integer                  :: i

    call delabove(i)

!   if (.not. A%is_defined) call quit('C = A + B: A undefined')
!   if (.not. B%is_defined) call quit('C = A + B: B undefined')

    call matrix_times_matrix_backend(A, B, C)

    call insref(C, i)
    call deltmparg(A, C)
    call deltmparg(B, C)

  end function

  function matrix_dot_matrix(A, B) result(r)

    type(matrix), intent(in) :: A, B
    real*8                   :: r
    integer                  :: i

!   if (.not. A%is_defined) call quit('r = A**B: A undefined')
!   if (.not. B%is_defined) call quit('r = A**B: B undefined')

    r = matrix_dot_matrix_backend(A, B)

    call deltmparg(A, r=r)
    call deltmparg(B, r=r)

  end function

  function matrix_product_trace(A, B) result(r)

    type(matrix), intent(in) :: A, B
    real*8                   :: r
    integer                  :: i

!   if (.not. A%is_defined) call quit('r = dag(A)**B: A undefined')
!   if (.not. B%is_defined) call quit('r = dag(A)**B: B undefined')

    r = dag(A)**B

    call deltmparg(A, r=r)
    call deltmparg(B, r=r)

  end function

  function matrix_is_zero(A) result(is_zero)

    type(matrix), intent(in) :: A
    logical                  :: is_zero
    integer                  :: i

    call delabove(i)

!   if (.not. A%is_defined) call quit('A undefined in matrix_is_zero(A)')

    is_zero = (.not. associated(A%elms))

  end function

  subroutine print_matrix(A, label, w)

    type(matrix), intent(in)           :: A
    character(*), intent(in), optional :: label
    integer,      intent(in), optional :: w
    integer                            :: i

    call delabove(i)

    call print_matrix_backend(A, label, w)

  end subroutine


! matrix operations which wrap around other operations
! ====================================================

  function matrix_norm(A) result(r)

    type(matrix), intent(in) :: A
    real*8                   :: r

    r = (A**A)**0.5d0

  end function

  function matrix_minus_matrix(A, B) result(C)

    type(matrix), intent(in) :: A, B
    type(matrix)             :: C

    C = A + (-B)

  end function

  function minus_matrix(A) result(B)

    type(matrix), intent(in) :: A
    type(matrix)             :: B

    B = (-1.0d0)*A

  end function

  function matrix_div_real(A, r) result(B)

    type(matrix), intent(in) :: A
    real*8,       intent(in) :: r
    type(matrix)             :: B

    B = (1.0d0/r)*A

  end function

end module
