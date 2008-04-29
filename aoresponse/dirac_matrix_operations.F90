module dirac_matrix_operations_m

! author:           andreas johan thorvaldsen
! porting to dirac: radovan bast

! this file must be compiled without inlining

! manual allocation and deallocation of matrices is cumbersome to program
! this module does book-keeping of the allocated data and their corresponding
! references (aliases) and performs deallocation once all referrers have
! come out of scope (on the stack) or been overwritten

  use kinds_m
  use memory_allocator

  implicit none

#include "dgroup.h"

  type matrix
    integer(kind_i4)       :: nrow = -1 ! operations on uninitialized matrices will fail
    integer(kind_i4)       :: ncol
    integer                :: nz
    real(kind_dp), pointer :: elements(:)
    integer                :: ih
  end type matrix

  public matrix

  public assignment(=)
  public operator(+)
  public operator(-)
  public operator(*)
  public operator(/)
  public operator(**)

  public mat_init
  public mat_free
  public mat_print
  public trab
  public trps
  public norm
  public iszero
  public isdef

  private

  interface assignment(=)
    module procedure matrix_eq_matrix
    module procedure matrix_eq_zero
  end interface

  interface operator(+)
    module procedure matplusmat
  end interface

  interface operator(-)
    module procedure matminusmat
    module procedure matneg
  end interface

  interface operator(*)
    module procedure dp_times_matrix
    module procedure inttimesmat
    module procedure matmulmat
  end interface

  interface operator(/)
    module procedure matdivreal
    module procedure matdivint
  end interface

  interface operator(**)
    module procedure matdotmat
  end interface

  interface trab
    module procedure trmatmulmat
  end interface

  interface trps
    module procedure trpsmat
  end interface

  interface norm
    module procedure normmat
  end interface

  interface iszero
    module procedure iszeromat
  end interface

  interface isdef
    module procedure isdefmat
  end interface














    !no array-of-pointer in fortran, so use array
    !of derived type containing pointer instead
    type refp
        type(matrix),pointer:: m
    end type

    type datp
        real*8,pointer:: d(:)
    end type

    !number of references, current book-keeping buffer size
    integer:: nref=0,maxref=0
    !number of matrix datas allocated, current book-keeping buffer size
    integer:: ndat=0,maxdat=0
    !growth direction of the stack, -1 or +1
    integer:: stkdir !different for i686,itanium,xeon, fex.

    !pointers to all type(matrix) referring to a matrix data (%elements)
    type(refp),allocatable:: ref(:) !1:numref<maxref
    !one pointer to each of the allocated %elements(:)
    type(datp),allocatable:: dat(:) !1:numdat<maxdat
    !index in dat referred to by each ref
    integer,allocatable:: idx(:)
    !number of references to each dat
    integer,allocatable:: cnt(:)


contains

#ifndef var_ifc
#define baddress loc
#endif /* var_ifc */

  subroutine mat_init(a, nrow, ncol)
    type(matrix)                 :: a
    integer(kind_i4), intent(in) :: nrow, ncol

    nullify(a%elements)

    a%nrow = nrow
    a%ncol = ncol
    a%nz   = 4

    call alloc(a%elements, a%nrow*a%ncol*a%nz)

  end subroutine mat_init

  subroutine mat_free(a)
    type(matrix) :: a

    if (.not. associated(a%elements)) then
      call quit('memory previously released')
    end if

    call dealloc(a%elements)
    nullify(a%elements)

  end subroutine mat_free

  subroutine mat_print(a, label)
    type(matrix), intent(inout)        :: a
    character(*), intent(in), optional :: label

#include "priunit.h"

    integer(kind_i4)            :: i

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in mat_print(a)')
    end if

    write(*, *) 'matrix ', label
    write(*, *) 'ih = ', a%ih

    if (associated(a%elements)) then
      call prqmat(a%elements,   &
                  a%nrow,       &
                  a%ncol,       &
                  a%nrow,       &
                  a%ncol,       &
                  a%nz,         &
                  ipqtoq(1, 0), &
                  lupri)
    else
      write(*, *) 'all elements zero'
    end if

  end subroutine mat_print

  subroutine matrix_eq_matrix(a, b)
    type(matrix), intent(inout) :: a
    type(matrix), intent(in)    :: b

    integer(kind_i4)            :: i

    call delabove(i)

    if (b%nrow == -1) then
      call quit('b undefined in a = b')
    end if

    if (baddress(a) /= baddress(b)) call delref(a)

    a%nrow = b%nrow
    a%ncol = b%ncol
    a%nz   = b%nz  
    a%ih   = b%ih  

    a%elements => b%elements

    call insref(a)

  end subroutine matrix_eq_matrix

  subroutine matrix_eq_zero(a, zero)
    type(matrix),     intent(inout) :: a
    integer(kind_i4), intent(in)    :: zero

    integer(kind_i4)                :: i, j

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in a = 0')
    end if
    if (zero /= 0) then
      call quit('z must be zero in a = z')
    end if

    i = a%nrow
    j = a%ncol

    call delref(a)

    a%nrow = i
    a%ncol = j

  end subroutine matrix_eq_zero

  function dp_times_matrix(r, a) result(b)
    real(kind_dp), intent(in) :: r
    type(matrix),  intent(in) :: a
    type(matrix)              :: b

    integer(kind_i4)          :: i

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in b = r*a')
    end if

    b%nrow = a%nrow
    b%ncol = a%ncol
    b%nz   = a%nz  
    b%ih   = a%ih  

    nullify(b%elements)

    if (r /= 0 .and. associated(a%elements)) then
      call alloc(b%elements, b%nrow*b%ncol*b%nz)
      b%elements = r*a%elements
    end if

    call insref(b)

  end function dp_times_matrix







!!!!!!!!!!!!!!!!
!       everything below needs revision - in progress ...  --radovan
!!!!!!!!!!!!!!!!



    subroutine allocif()
    !if needed allocate/extend the tables ref,dat,idx,cnt
        type(refp):: tref(nref); type(datp):: tdat(maxdat)
        integer:: tidx(nref),tcnt(maxdat),i
        if (nref<maxref-1.and.ndat<maxdat-1) return
        if (nref>=maxref-1) then
            if (maxref==0) call getdir() !only once
            if (maxref/=0) tref = ref(1:nref)
            if (maxref/=0) tidx = idx(1:nref)
            if (maxref/=0) deallocate(ref,idx)
            allocate(ref(max(64,2*maxref)))
            allocate(idx(max(64,2*maxref)))
            ref(1:nref) = tref
            idx(1:nref) = tidx
            maxref = max(64,2*maxref)
        endif
        if (ndat>=maxdat-1) then
            if (maxdat/=0) tdat = dat(1:maxdat)
            if (maxdat/=0) tcnt = cnt(1:maxdat)
            if (maxdat/=0) deallocate(dat,cnt)
            allocate(dat(max(8,2*maxdat)))
            allocate(cnt(max(8,2*maxdat))); cnt=0
            dat(1:maxdat) = tdat
            cnt(1:maxdat) = tcnt
            maxdat = max(8,2*maxdat)
        endif
    contains
        subroutine getdir()
            integer:: j
            stkdir = merge(-1,1,baddress(j)<baddress(i))
        end subroutine
    end subroutine


    subroutine delabove(k)
        integer:: k,i,j
        call allocif()
        do i=nref,1,-1
            if ((baddress(ref(i)%m)-baddress(k))*stkdir<0) exit
            j = idx(i)
            cnt(j) = cnt(j)-1
            if (cnt(j)==0) then
                deallocate(dat(j)%d)
                ndat = ndat-1
            endif
        enddo
        nref = i
    end subroutine


    subroutine insref(A)
        type(Matrix),intent(in),target:: A
        integer:: i,j,k; logical:: found,zero,same
        !see if any reference from A is recorded
        do i=nref,0,-1
            if (i==0) exit
            if ((baddress(ref(i)%m)-baddress(A))*stkdir<0) exit
        enddo; j=i+1; k=idx(j)
        found  = (j/=nref+1.and.associated(ref(j)%m,A))
        zero = (.not.associated(A%elements).or.abs(baddress(A%elements)-baddress(i))<2**28)
        same = .false.; if (found) same = associated(dat(k)%d,A%elements)
        !if A is recorded, but now is zero or A%elements is on the stack, remove the previous record
        if (found.and.zero) then
             ref(j:nref-1) = ref(j+1:nref)
             idx(j:nref-1) = idx(j+1:nref)
             nref = nref-1
        endif
        !if A was recorded as referring to some other data, remove that
        if (found.and..not.same) then
            cnt(k) = cnt(k)-1
            if (cnt(k)==0) deallocate(dat(k)%d)
            if (cnt(k)==0) ndat = ndat-1
        endif
        !see if A%elements is already recorded as some dat(k)%d
        do i=maxdat,1,-1
            if (cnt(i)==0) then
                k = i; cycle
            elseif (associated(A%elements,dat(i)%d)) then
                k = i; exit
            endif
        enddo
        !if A is not already recorded, make space for it
        if (.not.found.and..not.zero) then
            ref(j+1:nref+1) = ref(j:nref)
            idx(j+1:nref+1) = idx(j:nref)
            nref = nref+1
        endif
        !insert A
        if (.not.zero) then
            ref(j)%m => A
            idx(j) = k
            dat(k)%d => A%elements
            if (cnt(k)==0) ndat = ndat+1
            if (.not.same) cnt(k) = cnt(k)+1
        endif
    end subroutine

    subroutine delref(A)
    !if A is not recorded, deallocate A%elements
    !if A is recorded as referring to dat(k)%d, but refers to A%elements - deallocate(A%elements)
    !if A is recorded, referring to recorded dat(k)%d - decrement refcount for dat(k)%d
        type(Matrix),intent(inout),target:: A
        integer:: i,j,k
        !find the reference, if exists
        do j=nref,1,-1
            if ((baddress(ref(j)%m)-baddress(A))*stkdir<=0) exit
        enddo; k=idx(j)
        !if A%elements is not the recorded dat(k)%d, see if A%elements is referenced by some
        !other matrix. If not, deallocate A%elements
        if (A%nrow/=-1.and.associated(A%elements).and..not. &
            (associated(ref(j)%m,A).and.associated(dat(k)%d,A%elements))) then
            do i=maxdat,0,-1
                if (cnt(i)==0) cycle
                if (associated(A%elements,dat(i)%d)) exit
                if (i==1) deallocate(A%elements)
            enddo
        endif
        !if a reference from A was found, remove it
        if (associated(ref(j)%m,A)) then
            ref(j:nref-1) = ref(j+1:nref)
            idx(j:nref-1) = idx(j+1:nref)
            nref = nref-1
            cnt(k) = cnt(k)-1
            if (cnt(k)==0) then
                deallocate(dat(k)%d)
                ndat = ndat-1
            endif
        endif
        A%nrow=-1; nullify(A%elements)
    end subroutine


    function MatPlusMat(A,B) result(C)
        type(Matrix),intent(in):: A,B
        type(Matrix):: C; integer:: i
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix A+B : A undefined'
        if (B%nrow==-1) stop 'Matrix A+B : B undefined'
        if (A%nrow/=B%nrow.or.A%ncol/=B%ncol) stop 'Matrix A+B : different shapes'
        C%nrow=A%nrow; C%ncol=A%ncol; nullify(C%elements)
        if (associated(A%elements).and.associated(B%elements)) then
            allocate(C%elements(C%nrow*C%ncol))
            C%elements = A%elements+B%elements
        elseif (associated(A%elements)) then
            C%elements => A%elements
        elseif (associated(B%elements)) then
            C%elements => B%elements
        endif
        call insref(C)
    end function

    function MatMinusMat(A,B) result(C)
        type(Matrix),intent(in):: A,B
        type(Matrix):: C
        C = A + (-B)
    end function


    function MatNeg(A) result(B)
        type(Matrix),intent(in):: A
        type(Matrix):: B
        B = (-1d0)*A
    end function

    function IntTimesMat(i,A) result(B)
        type(Matrix),intent(in):: A
        integer,intent(in):: i; type(Matrix):: B
        B = (1d0*i)*A
    end function

    function MatDivReal(A,r) result(B)
        type(Matrix),intent(in):: A
        real*8,intent(in):: r; type(Matrix):: B
        B = (1d0/r)*A
    end function

    function MatDivInt(A,i) result(B)
        type(Matrix),intent(in):: A
        integer,intent(in):: i; type(Matrix):: B
        B = (1d0/i)*A
    end function

    function MatMulMat(A,B) result(C)
        type(Matrix),intent(in):: A,B
        type(Matrix):: C
        integer:: i,k,l,m,n
        call delabove(i)
        if (A%nrow==-1) stop 'Matrix A*B : A undefined'
        if (B%nrow==-1) stop 'Matrix A*B : B undefined'
        if (A%ncol/=B%nrow) stop 'Matrix A*B : incompatible shapes'
        C%nrow=A%nrow; C%ncol=B%ncol
        if (associated(A%elements).and.associated(B%elements)) then
            allocate(C%elements(C%nrow*C%ncol))
            C%elements = reshape(matmul(reshape(A%elements(1:A%nrow*A%ncol),(/A%nrow,A%ncol/)), &
                     reshape(B%elements(1:B%nrow*B%ncol),(/B%nrow,B%ncol/))),(/C%nrow*C%ncol/))
        endif
        call insref(C)
    end function

    function MatDotMat(A,B) result(r)
        type(Matrix),intent(in):: A,B; real*8:: r
        if (A%nrow==-1) stop 'Matrix A**B : A undefined'
        if (B%nrow==-1) stop 'Matrix A**B : B undefined'
        if (A%nrow/=B%nrow) stop 'Matrix A**B : different number of rows'
        if (A%ncol/=B%ncol) stop 'Matrix A**B : different number of cols'
        r = 0d0
        if (.not.(associated(A%elements).and.associated(B%elements))) return
        r = sum(A%elements*B%elements)
    end function

    function TrMatMulMat(A,B) result(r)
        type(Matrix),intent(in):: A,B; real*8:: r
        if (A%nrow==-1) stop 'Matrix trab(A,B) : A undefined'
        if (B%nrow==-1) stop 'Matrix trab(A,B) : B undefined'
        if (A%nrow/=B%ncol) stop 'Matrix trab(A,B) : not compatible shapes'
        if (A%ncol/=B%nrow) stop 'Matrix trab(A,B) : not compatible shapes'
        r = 0d0
        if (.not.(associated(A%elements).and.associated(B%elements))) return
        r = sum(transpose(reshape(A%elements,(/A%nrow,A%ncol/)))*reshape(B%elements,(/B%nrow,B%ncol/)))
    end function

    function TrpsMat(A) result(B)
        type(Matrix),intent(in):: A; type(Matrix):: B
        integer:: i; call delabove(i)
        if (A%nrow==-1) stop 'Matrix trps(A) : A undefined'
        B%nrow=A%ncol; B%ncol=A%nrow
        if (associated(A%elements)) then
            allocate(B%elements(B%nrow*B%ncol))
            B%elements = reshape(transpose(reshape(A%elements,(/A%nrow,A%ncol/))),(/B%nrow*B%ncol/))
        endif
        call insref(B)
    end function

    function NormMat(A) result(r)
        type(Matrix),intent(in):: A; real*8:: r
        if (A%nrow==-1) stop 'Matrix A**B : A undefined'
        r = 0d0
        if (.not.associated(A%elements)) return
        r = sqrt(sum(A%elements**2))
    end function

    function IsZeroMat(A) result(z)
        type(Matrix),intent(in):: A; logical:: z
        integer:: i; call delabove(i)
        if (A%nrow==-1) stop 'Matrix iszero(A) : A undefined'
        z = (A%nrow/=-1.and..not.associated(A%elements))
    end function

    function IsDefMat(A) result(d)
        type(Matrix),intent(in):: A; logical:: d
        integer:: i; call delabove(i)
        d = (A%nrow/=-1)
    end function


end module dirac_matrix_operations_m
