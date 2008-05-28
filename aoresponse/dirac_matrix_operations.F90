module dirac_matrix_operations_m

! author:           andreas johan thorvaldsen
! porting to dirac: radovan bast

#define baddress loc
! this file must be compiled without inlining

! manual allocation and deallocation of matrices is cumbersome to program
! this module does book-keeping of the allocated data and their corresponding
! refs (aliases) and performs deallocation once all referrers have
! come out of scope (on the stack) or been overwritten

  use kinds_m
  use memory_allocator

  implicit none

#include "dgroup.h"

  type matrix
    integer       :: nrow              = -1 ! operations on uninitialized matrices will fail
    integer       :: ncol              
    integer                :: nz                = 4
    integer                :: ih_sym            = 1
    integer                :: tr_sym            = 1
    integer                :: pg_sym            = 1
    integer                :: fermi_pointer     = 1
    real(kind_dp), pointer :: elements(:)
    logical                :: is_density_matrix = .false.
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
  public dag
  public norm
  public iszero
  public isdef
  public matrixstats

  private

  interface assignment(=)
    module procedure matrix_eq_matrix
    module procedure matrix_eq_zero
  end interface

  interface operator(+)
    module procedure matrix_plus_matrix
  end interface

  interface operator(-)
    module procedure matrix_minus_matrix
    module procedure matrix_negative
  end interface

  interface operator(*)
    module procedure dp_times_matrix
!   module procedure i4_times_matrix
    module procedure i_times_matrix
    module procedure matrix_times_matrix
  end interface

  interface operator(/)
    module procedure divide_matrix_by_dp
!   module procedure divide_matrix_by_i4
    module procedure divide_matrix_by_i
  end interface

  interface operator(**)
    module procedure matrix_dot_matrix
  end interface

  interface matrixstats
    module procedure matrix_stat
  end interface

  interface norm
    module procedure matrix_norm
  end interface

  interface trab
      module procedure TrMatMulMat
  end interface

  interface trps
    module procedure dag_matrix
  end interface

  interface dag
    module procedure dag_matrix
  end interface

  interface iszero
    module procedure matrix_is_zero
  end interface

  interface isdef
    module procedure matrix_is_defined
  end interface

    !no array-of-pointer in Fortran, so use array
    !of derived type containing pointer instead
    type refp
        type(Matrix),pointer:: m
    end type

    type datp
        real*8,pointer:: d(:)
    end type

    !number of references, current book-keeping buffer size
    integer:: nref=0,maxref=0
    !number of matrix datas allocated, current book-keeping buffer size
    integer:: ndat=0,maxdat=0
    !growth direction of the stack, -1 or +1
    integer:: stkdir !different for i686,Itanium,XEON, fex.

    !pointers to all type(Matrix) referring to a matrix data (%elements)
    type(refp),allocatable:: ref(:) !1:numref<maxref
    !one pointer to each of the allocated %elements(:)
    type(datp),allocatable:: dat(:) !1:numdat<maxdat
    !index in dat referred to by each ref
    integer,allocatable:: idx(:)
    !number of references to each dat
    integer,allocatable:: cnt(:)

contains

  subroutine mat_init(a, nrow, ncol)

    type(matrix), intent(inout) :: a
    integer,      intent(in)    :: nrow, ncol

    a%nrow = nrow
    a%ncol = ncol

    nullify(a%elements)
    allocate(a%elements(a%nrow*a%ncol*a%nz))

  end subroutine mat_init





  subroutine mat_print(a, label)
    type(matrix), intent(inout)        :: a
    character(*), intent(in), optional :: label

#include "priunit.h"

    integer            :: i

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in mat_print(a)')
    end if

    write(*, *) 'matrix        = ', label
    write(*, *) 'ih_sym        = ', a%ih_sym       
    write(*, *) 'tr_sym        = ', a%tr_sym       
    write(*, *) 'pg_sym        = ', a%pg_sym       
    write(*, *) 'fermi_pointer = ', a%fermi_pointer

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

    integer            :: i

    call delabove(i)

    if (b%nrow == -1) then
      call quit('b undefined in a = b')
    end if

    if(.not. associated(a%elements,b%elements)) call delref(A)

    a%nrow          = b%nrow
    a%ncol          = b%ncol
    a%nz            = b%nz
    a%ih_sym        = b%ih_sym
    a%tr_sym        = b%tr_sym
    a%pg_sym        = b%pg_sym
    a%fermi_pointer = b%fermi_pointer

    a%elements => b%elements

    call insref(a)

  end subroutine matrix_eq_matrix

  subroutine matrix_eq_zero(a, zero)
    type(matrix),     intent(inout) :: a
    integer,          intent(in)    :: zero

    integer                :: i, j

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

  subroutine mat_free(a)
    type(matrix),     intent(inout) :: a

    integer                :: i, j

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in a = 0')
    end if
!   if (zero /= 0) then
!     call quit('z must be zero in a = z')
!   end if

    i = a%nrow
    j = a%ncol

    call delref(a)

    a%nrow = i
    a%ncol = j

  end subroutine mat_free


  function dp_times_matrix(r, a) result(b)
    real(kind_dp), intent(in) :: r
    type(matrix),  intent(in) :: a
    type(matrix)              :: b

    integer          :: i

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in b = r*a')
    end if

    b%nrow          = a%nrow
    b%ncol          = a%ncol
    b%nz            = a%nz
    b%ih_sym        = a%ih_sym
    b%tr_sym        = a%tr_sym
    b%pg_sym        = a%pg_sym
    b%fermi_pointer = a%fermi_pointer

    nullify(b%elements)

    if (r /= 0 .and. associated(a%elements)) then
      allocate(b%elements(b%nrow*b%ncol*b%nz))
      b%elements = r*a%elements
    end if

    call insref(b)

  end function dp_times_matrix

  function matrix_plus_matrix(a, b) result(c)
    type(matrix),  intent(in) :: a, b
    type(matrix)              :: c

    integer          :: i

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in c = a+b')
    end if
    if (b%nrow == -1) then
      call quit('b undefined in c = a+b')
    end if
    if (a%nrow /= b%nrow) then
      call quit('different shapes in c = a+b')
    end if
    if (a%ncol /= b%ncol) then
      call quit('different shapes in c = a+b')
    end if
    if (a%nz /= b%nz) then
      call quit('different shapes in c = a+b')
    end if

    c%nrow = a%nrow
    c%ncol = a%ncol
    c%nz   = a%nz

    nullify(c%elements)

    if (associated(a%elements) .and. associated(b%elements)) then
      allocate(c%elements(c%nrow*c%ncol*c%nz))
      c%elements = a%elements + b%elements
    else if (associated(a%elements)) then
      c%elements => a%elements
    else if (associated(b%elements)) then
      c%elements => b%elements
    end if

    call insref(c)

  end function matrix_plus_matrix

  function matrix_minus_matrix(a, b) result(c)
    type(matrix),  intent(in) :: a, b
    type(matrix)              :: c

    c = a + (-b)

  end function matrix_minus_matrix

  function matrix_negative(a) result(b)
    type(matrix),  intent(in) :: a
    type(matrix)              :: b

    b = (-1.0d0)*a

  end function matrix_negative

  function i_times_matrix(i, a) result(b)
    integer,          intent(in) :: i
    type(matrix),     intent(in) :: a
    type(matrix)                 :: b

    b = (-1.0d0*i)*a

  end function i_times_matrix

  function i4_times_matrix(i, a) result(b)
    integer, intent(in) :: i
    type(matrix),     intent(in) :: a
    type(matrix)                 :: b

    b = (-1.0d0*i)*a

  end function i4_times_matrix

  function matrix_times_matrix(a, b) result(c)

    type(matrix), intent(in) :: a, b
    type(matrix)             :: c

    integer         :: i, k, l, m, n

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in c = a*b')
    end if
    if (b%nrow == -1) then
      call quit('b undefined in c = a*b')
    end if
    if (a%ncol /= b%nrow) then
      call quit('incompatible shapes in c = a*b')
    end if
    if (a%nz /= b%nz) then
      call quit('incompatible shapes in c = a*b')
    end if

    c%nrow = a%nrow
    c%ncol = b%ncol
    c%nz   = a%nz

    if (associated(a%elements) .and. associated(b%elements)) then

      allocate(c%elements(c%nrow*c%ncol*c%nz))

        call qgemm(a%nrow, b%ncol, a%ncol, 1.0d0,  &
                   'N', 'N', ipqtoq(1, 0), a%elements, a%nrow, a%ncol, a%nz, &
                   'N', 'N', ipqtoq(1, 0), b%elements, b%nrow, b%ncol, b%nz, &
                   0.0d0,    ipqtoq(1, 0), c%elements, c%nrow, c%ncol, c%nz)

    end if

    call insref(c)

  end function matrix_times_matrix

  function matrix_dot_matrix(a, b) result(r)
    type(matrix),  intent(in) :: a, b
    real(kind_dp)             :: r

    if (a%nrow == -1) then
      call quit('a undefined in r = a**b')
    end if
    if (b%nrow == -1) then
      call quit('b undefined in r = a**b')
    end if
    if (a%nrow /= b%nrow) then
      call quit('different shapes in r = a**b')
    end if
    if (a%ncol /= b%ncol) then
      call quit('different shapes in r = a**b')
    end if
    if (a%nz /= b%nz) then
      call quit('different shapes in r = a**b')
    end if

    r = 0.0d0

    if (.not.(associated(a%elements) .and. associated(b%elements))) return

    r = sum(a%elements*b%elements)

!   if (a%is_density_matrix) then
!     r = r*2.0d0
!   end if
!   if (b%is_density_matrix) then
!     r = r*2.0d0
!   end if
    r = r*2.0d0

  end function matrix_dot_matrix

  function matrix_norm(a) result(r)
    type(matrix),  intent(in) :: a
    real(kind_dp)             :: r

    if (a%nrow == -1) then
      call quit('a undefined in matrix_norm')
    end if

    r = 0.0d0
    if (.not. associated(a%elements)) return
    r = sqrt(sum(a%elements**2.0d0))

  end function matrix_norm

! function divide_matrix_by_dp(a, r) result(b)
!   type(matrix),  intent(in) :: a
!   real(kind_dp), intent(in) :: r
!   type(matrix)              :: b

!   b = (1.0d0/r)*a

! end function divide_matrix_by_dp

    function divide_matrix_by_dp(A,r) result(B)
        type(Matrix),intent(in):: A
        real*8,intent(in):: r; type(Matrix):: B
        B = (1d0/r)*A
    end function



  function divide_matrix_by_i(a, i) result(b)
    type(matrix),     intent(in) :: a
    integer,          intent(in) :: i
    type(matrix)                 :: b

    b = (1.0d0/i)*a

  end function divide_matrix_by_i

  function divide_matrix_by_i4(a, i) result(b)
    type(matrix),     intent(in) :: a
    integer, intent(in) :: i
    type(matrix)                 :: b

    b = (1.0d0/i)*a

  end function divide_matrix_by_i4





!!!!!!!!!!!!!!!!
!       everything below needs revision - in progress ...  --radovan
!!!!!!!!!!!!!!!!





  subroutine transpose_block(factor, block, nrow, ncol, tblock)

    real(kind_dp),    intent(in)    :: factor
    integer, intent(in)    :: nrow
    integer, intent(in)    :: ncol
    real(kind_dp),    intent(in)    :: block(nrow*ncol)
    real(kind_dp),    intent(inout) :: tblock(nrow*ncol)

    tblock = reshape(transpose(reshape(block, (/nrow, ncol/))), (/nrow*ncol/))

    call dscal(nrow*ncol, factor, tblock, 1)

  end subroutine transpose_block

  function dag_matrix(a) result(b)

    type(matrix), intent(in) :: a
    type(matrix)             :: b

    integer         :: i

    call delabove(i)

    if (a%nrow == -1) then
      call quit('a undefined in a^\dag')
    end if

    b%ncol          = a%nrow
    b%nrow          = a%ncol
    b%nz            = a%nz
    b%ih_sym        = a%ih_sym
    b%tr_sym        = a%tr_sym
    b%pg_sym        = a%pg_sym
    b%fermi_pointer = a%fermi_pointer

    if (associated(a%elements)) then
      allocate(b%elements(b%nrow*b%ncol*b%nz))

      call transpose_block( 1.0d0, a%elements(0*a%nrow*a%ncol + 1:1*a%nrow*a%ncol), a%nrow, a%ncol, &
                                   b%elements(0*a%nrow*a%ncol + 1:1*a%nrow*a%ncol))
      call transpose_block(-1.0d0, a%elements(1*a%nrow*a%ncol + 1:2*a%nrow*a%ncol), a%nrow, a%ncol, &
                                   b%elements(1*a%nrow*a%ncol + 1:2*a%nrow*a%ncol))
      call transpose_block(-1.0d0, a%elements(2*a%nrow*a%ncol + 1:3*a%nrow*a%ncol), a%nrow, a%ncol, &
                                   b%elements(2*a%nrow*a%ncol + 1:3*a%nrow*a%ncol))
      call transpose_block(-1.0d0, a%elements(3*a%nrow*a%ncol + 1:4*a%nrow*a%ncol), a%nrow, a%ncol, &
                                   b%elements(3*a%nrow*a%ncol + 1:4*a%nrow*a%ncol))
    endif

    call insref(b)

  end function dag_matrix

    function TrMatMulMat(A,B) result(r)
        type(Matrix),intent(in):: A,B; real*8:: r
        r = dag(A)**B

    end function

    function matrix_is_zero(a) result(z)
        type(matrix),intent(in):: a; logical:: z
        integer:: i; call delabove(i)
        if (a%nrow==-1) stop 'matrix iszero(a) : a undefined'
        z = (a%nrow/=-1.and..not.associated(a%elements))
    end function

    function matrix_is_defined(a) result(d)
        type(matrix),intent(in):: a; logical:: d
        integer:: i; call delabove(i)
        d = (a%nrow/=-1)
    end function



    subroutine matrix_stat(label)
        character(*), intent(in), optional :: label
        integer:: i,j
        write(*, *) label
        print '("ndat="i2" maxdat="i2" nref="i2" maxref="i2)', &
                ndat,maxdat,nref,maxref
        i = maxval(idx(1:nref))
        if (nref==0.or.ndat==0) return
        print '("ref="<nref>i13/"dat="<i>i13/"idx="<nref>i2"  cnt="<i>i2/)', &
                (loc(ref(j)%m),j=1,nref),(loc(dat(j)%d),j=1,i), &
                idx(1:nref),cnt(1:i)
        call delabove(i)
    end subroutine

! subroutine matrix_stat(label)

!   character(*), intent(in), optional :: label

!   integer :: i, j

!   write(*, *)
!   write(*, *) 'label  : ', label
!   write(*, *) 'ndat   : ', ndat
!   write(*, *) 'maxdat : ', maxdat
!   write(*, *) 'nref   : ', nref
!   write(*, *) 'maxref : ', maxref

!   if (nref == 0 .or. ndat == 0) return

!   i = maxval(idx(1:nref))

!   write(*, *) '      ref       pointer          idx'
!   do j = 1, nref
!     write(*, '(I10,I14,I13)') j, loc(ref(j)%m), idx(j)
!   end do

!   write(*, *) '     data       pointer          cnt'
!   do j = 1, i
!     write(*, '(I10,I14,I13)') j, loc(dat(j)%d), cnt(j)
!   end do

!   call delabove(i)

! end subroutine matrix_stat




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
            stkdir = merge(-1,1,loc(j)<loc(i))
        end subroutine
    end subroutine


    subroutine delabove(k)
        integer:: k,i,j, loc_k
        call allocif()
!       loc_k = loc(k)
!       write(*, *) 'loc_k', loc_k
        do i=nref,1,-1
            if ((loc(ref(i)%m)-loc(k))*stkdir<0) exit
            j = idx(i)
            cnt(j) = cnt(j)-1
!if (cnt(j)/=0) print '("delabove delete ref "i13" -> "i13)',loc(ref(i)%m),loc(dat(j)%d)
!if (cnt(j)==0) print '("delabove delete ref "i13" -> "i13" og deallok dat "i13)',&
!loc(ref(i)%m),loc(dat(j)%d),loc(dat(j)%d)
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
            if ((loc(ref(i)%m)-loc(A))*stkdir<0) exit
        enddo; j=i+1; k=idx(j)
        found  = (j/=nref+1.and.associated(ref(j)%m,A))
        zero = (.not.associated(A%elements).or.abs(loc(A%elements)-loc(i))<2**28)
        same = .false.; if (found) same = associated(dat(k)%d,A%elements)
        !if A is recorded, but now is zero or A%elements is on the stack, remove the previous record
        if (found.and.zero) then
!print '("insref    delete ref "i13" -> "i13)',loc(ref(j)%m),loc(dat(idx(j))%d)
             ref(j:nref-1) = ref(j+1:nref)
             idx(j:nref-1) = idx(j+1:nref)
             nref = nref-1
        endif
        !if A was recorded as referring to some other data, remove that
        if (found.and..not.same) then
            cnt(k) = cnt(k)-1
!if (cnt(k)==0) print '("insref   deallok dat "i13)',loc(dat(k)%d)
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
!print '("insref   record ref "i13" -> "i13)',loc(A),loc(A%elements)
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
            if ((loc(ref(j)%m)-loc(A))*stkdir<=0) exit
        enddo; k=idx(j)
        !if A%elements is not the recorded dat(k)%d, see if A%elements is referenced by some
        !other matrix. If not, deallocate A%elements
        if (A%nrow/=-1.and.associated(A%elements).and..not. &
            (associated(ref(j)%m,A).and.associated(dat(k)%d,A%elements))) then
            do i=maxdat,0,-1
                if (cnt(i)==0) cycle
                if (associated(A%elements,dat(i)%d)) exit
!if (i==1) print '("delref   deallok dat "i13)',loc(A%elements)
                if (i==1) deallocate(A%elements)
            enddo
        endif
        !if a reference from A was found, remove it
        if (associated(ref(j)%m,A)) then
            ref(j:nref-1) = ref(j+1:nref)
            idx(j:nref-1) = idx(j+1:nref)
            nref = nref-1
            cnt(k) = cnt(k)-1
!if (cnt(k)/=0) print '("delref   delete ref "i13" -> "i13)',loc(A),loc(dat(k)%d)
!if (cnt(k)==0) print '("delref   delete ref "i13" -> "i13" og deallok dat "i13)', &
!loc(A),loc(dat(k)%d),loc(dat(k)%d)
            if (cnt(k)==0) then
                deallocate(dat(k)%d)
                ndat = ndat-1
            endif
        endif
        A%nrow=-1; nullify(A%elements)
    end subroutine


end module dirac_matrix_operations_m
