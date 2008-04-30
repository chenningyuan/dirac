module dirac_matrix_operations_m

! author:           andreas johan thorvaldsen
! porting to dirac: radovan bast

#define baddress loc
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
  public compare_dp_arrays
  public compare_matrices

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
    module procedure i4_times_matrix
    module procedure matrix_times_matrix
  end interface

  interface operator(/)
    module procedure matdivreal
    module procedure matdivint
  end interface

  interface operator(**)
    module procedure matrix_dot_matrix
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

  type ref_pointer_t
    type(matrix), pointer :: m
  end type ref_pointer_t

  type dat_pointer_t
    real(kind_dp), pointer :: d(:)
  end type dat_pointer_t

  integer(kind_i4) :: nr_ref   = 0
  integer(kind_i4) :: max_ref  = 0
  integer(kind_i4) :: nr_data  = 0
  integer(kind_i4) :: max_data = 0

  integer          :: stack_direction ! growth direction of the stack: -1 or +1

  type(ref_pointer_t), allocatable :: ref(:)
  type(dat_pointer_t), allocatable :: dat(:)

  integer, allocatable :: index_in_dat(:)
  integer, allocatable :: nr_ref_to_dat(:)

contains

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

  function matrix_plus_matrix(a, b) result(c)
    type(matrix),  intent(in) :: a, b
    type(matrix)              :: c

    integer(kind_i4)          :: i

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
      call alloc(c%elements, c%nrow*c%ncol*c%nz)
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

  function i4_times_matrix(i, a) result(b)
    integer(kind_i4), intent(in) :: i
    type(matrix),     intent(in) :: a
    type(matrix)                 :: b

    b = (-1.0d0*i)*a

  end function i4_times_matrix

  function matrix_times_matrix(a, b) result(c)
    type(matrix),  intent(in) :: a, b
    type(matrix)              :: c

    integer(kind_i4)          :: i, k, l, m, n

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
      call quit('different shapes in c = a*b')
    end if

    c%nrow = a%nrow
    c%ncol = b%ncol
    c%nz   = a%nz

    if (associated(a%elements) .and. associated(b%elements)) then
      call alloc(c%elements, c%nrow*c%ncol*c%nz)

      c%elements = reshape(matmul(reshape(a%elements(1:a%nrow*a%ncol), (/a%nrow, a%ncol/)), &
                   reshape(b%elements(1:b%nrow*b%ncol), (/b%nrow, b%ncol/))), (/c%nrow*c%ncol/))
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

  end function matrix_dot_matrix

  function compare_dp_arrays(a, b, accept) result(ok)
    real(kind_dp), intent(in)    :: a(:), b(:)
    real(kind_dp), intent(in)    :: accept
    logical                      :: ok

    integer(kind_i4)           :: i

    if (size(a) /= size(b)) then
      ok = .false.
    else
      ok = .true.

      do i = 1, size(a)
        if (abs(a(i) - b(i)) > accept) then
          ok = .false.
        end if
      end do
    end if

  end function compare_dp_arrays

  subroutine compare_matrices(a, b, accept, label)
    type(matrix),  intent(in)           :: a, b
    real(kind_dp), intent(in)           :: accept
    character(*),  intent(in), optional :: label

#include "priunit.h"

    logical                             :: shape_ok, values_ok

    shape_ok  = .true.

    if (a%nrow /= b%nrow) then
      shape_ok = .false.
    end if
    if (a%ncol /= b%ncol) then
      shape_ok = .false.
    end if
    if (a%nz /= b%nz) then
      shape_ok = .false.
    end if

    if (.not. shape_ok) then
      write(*, *) 'matrix shapes not ok -- (', label, ')'
    else
      values_ok = compare_dp_arrays(a%elements, b%elements, accept)
      if (.not. values_ok) then
        write(*, *) 'matrix values not ok -- (', label, ')'
      end if
    end if

    if (shape_ok .and. values_ok) then
        write(*, *) '                  ok -- (', label, ')'
    end if

  end subroutine compare_matrices









!!!!!!!!!!!!!!!!
!       everything below needs revision - in progress ...  --radovan
!!!!!!!!!!!!!!!!




    subroutine allocif()
    !if needed allocate/extend the tables ref,dat,index_in_dat,nr_ref_to_dat
        type(ref_pointer_t):: tref(nr_ref); type(dat_pointer_t):: tdat(max_data)
        integer:: tindex_in_dat(nr_ref),tnr_ref_to_dat(max_data),i
        if (nr_ref<max_ref-1.and.nr_data<max_data-1) return
        if (nr_ref>=max_ref-1) then
            if (max_ref==0) call getdir() !only once
            if (max_ref/=0) tref = ref(1:nr_ref)
            if (max_ref/=0) tindex_in_dat = index_in_dat(1:nr_ref)
            if (max_ref/=0) deallocate(ref,index_in_dat)
            allocate(ref(max(64,2*max_ref)))
            allocate(index_in_dat(max(64,2*max_ref)))
            ref(1:nr_ref) = tref
            index_in_dat(1:nr_ref) = tindex_in_dat
            max_ref = max(64,2*max_ref)
        endif
        if (nr_data>=max_data-1) then
            if (max_data/=0) tdat = dat(1:max_data)
            if (max_data/=0) tnr_ref_to_dat = nr_ref_to_dat(1:max_data)
            if (max_data/=0) deallocate(dat,nr_ref_to_dat)
            allocate(dat(max(8,2*max_data)))
            allocate(nr_ref_to_dat(max(8,2*max_data))); nr_ref_to_dat=0
            dat(1:max_data) = tdat
            nr_ref_to_dat(1:max_data) = tnr_ref_to_dat
            max_data = max(8,2*max_data)
        endif
    contains
        subroutine getdir()
            integer:: j
            stack_direction = merge(-1,1,baddress(j)<baddress(i))
        end subroutine
    end subroutine


    subroutine delabove(k)
        integer:: k,i,j
        call allocif()
        do i=nr_ref,1,-1
            if ((baddress(ref(i)%m)-baddress(k))*stack_direction<0) exit
            j = index_in_dat(i)
            nr_ref_to_dat(j) = nr_ref_to_dat(j)-1
            if (nr_ref_to_dat(j)==0) then
                deallocate(dat(j)%d)
                nr_data = nr_data-1
            endif
        enddo
        nr_ref = i
    end subroutine


    subroutine insref(a)
        type(matrix),intent(in),target:: a
        integer:: i,j,k; logical:: found,zero,same
        !see if any reference from a is recorded
        do i=nr_ref,0,-1
            if (i==0) exit
            if ((baddress(ref(i)%m)-baddress(a))*stack_direction<0) exit
        enddo; j=i+1; k=index_in_dat(j)
        found  = (j/=nr_ref+1.and.associated(ref(j)%m,a))
        zero = (.not.associated(a%elements).or.abs(baddress(a%elements)-baddress(i))<2**28)
        same = .false.; if (found) same = associated(dat(k)%d,a%elements)
        !if a is recorded, but now is zero or a%elements is on the stack, remove the previous record
        if (found.and.zero) then
             ref(j:nr_ref-1) = ref(j+1:nr_ref)
             index_in_dat(j:nr_ref-1) = index_in_dat(j+1:nr_ref)
             nr_ref = nr_ref-1
        endif
        !if a was recorded as referring to some other data, remove that
        if (found.and..not.same) then
            nr_ref_to_dat(k) = nr_ref_to_dat(k)-1
            if (nr_ref_to_dat(k)==0) deallocate(dat(k)%d)
            if (nr_ref_to_dat(k)==0) nr_data = nr_data-1
        endif
        !see if a%elements is already recorded as some dat(k)%d
        do i=max_data,1,-1
            if (nr_ref_to_dat(i)==0) then
                k = i; cycle
            elseif (associated(a%elements,dat(i)%d)) then
                k = i; exit
            endif
        enddo
        !if a is not already recorded, make space for it
        if (.not.found.and..not.zero) then
            ref(j+1:nr_ref+1) = ref(j:nr_ref)
            index_in_dat(j+1:nr_ref+1) = index_in_dat(j:nr_ref)
            nr_ref = nr_ref+1
        endif
        !insert a
        if (.not.zero) then
            ref(j)%m => a
            index_in_dat(j) = k
            dat(k)%d => a%elements
            if (nr_ref_to_dat(k)==0) nr_data = nr_data+1
            if (.not.same) nr_ref_to_dat(k) = nr_ref_to_dat(k)+1
        endif
    end subroutine

    subroutine delref(a)
    !if a is not recorded, deallocate a%elements
    !if a is recorded as referring to dat(k)%d, but refers to a%elements - deallocate(a%elements)
    !if a is recorded, referring to recorded dat(k)%d - decrement refcount for dat(k)%d
        type(matrix),intent(inout),target:: a
        integer:: i,j,k
        !find the reference, if exists
        do j=nr_ref,1,-1
            if ((baddress(ref(j)%m)-baddress(a))*stack_direction<=0) exit
        enddo; k=index_in_dat(j)
        !if a%elements is not the recorded dat(k)%d, see if a%elements is referenced by some
        !other matrix. if not, deallocate a%elements
        if (a%nrow/=-1.and.associated(a%elements).and..not. &
            (associated(ref(j)%m,a).and.associated(dat(k)%d,a%elements))) then
            do i=max_data,0,-1
                if (nr_ref_to_dat(i)==0) cycle
                if (associated(a%elements,dat(i)%d)) exit
                if (i==1) deallocate(a%elements)
            enddo
        endif
        !if a reference from a was found, remove it
        if (associated(ref(j)%m,a)) then
            ref(j:nr_ref-1) = ref(j+1:nr_ref)
            index_in_dat(j:nr_ref-1) = index_in_dat(j+1:nr_ref)
            nr_ref = nr_ref-1
            nr_ref_to_dat(k) = nr_ref_to_dat(k)-1
            if (nr_ref_to_dat(k)==0) then
                deallocate(dat(k)%d)
                nr_data = nr_data-1
            endif
        endif
        a%nrow=-1; nullify(a%elements)
    end subroutine






    function matdivreal(a,r) result(b)
        type(matrix),intent(in):: a
        real*8,intent(in):: r; type(matrix):: b
        b = (1d0/r)*a
    end function

    function matdivint(a,i) result(b)
        type(matrix),intent(in):: a
        integer,intent(in):: i; type(matrix):: b
        b = (1d0/i)*a
    end function


    function trmatmulmat(a,b) result(r)
        type(matrix),intent(in):: a,b; real*8:: r
        if (a%nrow==-1) stop 'matrix trab(a,b) : a undefined'
        if (b%nrow==-1) stop 'matrix trab(a,b) : b undefined'
        if (a%nrow/=b%ncol) stop 'matrix trab(a,b) : not compatible shapes'
        if (a%ncol/=b%nrow) stop 'matrix trab(a,b) : not compatible shapes'
        r = 0d0
        if (.not.(associated(a%elements).and.associated(b%elements))) return
        r = sum(transpose(reshape(a%elements,(/a%nrow,a%ncol/)))*reshape(b%elements,(/b%nrow,b%ncol/)))
    end function

    function trpsmat(a) result(b)
        type(matrix),intent(in):: a; type(matrix):: b
        integer:: i; call delabove(i)
        if (a%nrow==-1) stop 'matrix trps(a) : a undefined'
        b%nrow=a%ncol; b%ncol=a%nrow
        if (associated(a%elements)) then
            allocate(b%elements(b%nrow*b%ncol))
            b%elements = reshape(transpose(reshape(a%elements,(/a%nrow,a%ncol/))),(/b%nrow*b%ncol/))
        endif
        call insref(b)
    end function

    function normmat(a) result(r)
        type(matrix),intent(in):: a; real*8:: r
        if (a%nrow==-1) stop 'matrix a**b : a undefined'
        r = 0d0
        if (.not.associated(a%elements)) return
        r = sqrt(sum(a%elements**2))
    end function

    function iszeromat(a) result(z)
        type(matrix),intent(in):: a; logical:: z
        integer:: i; call delabove(i)
        if (a%nrow==-1) stop 'matrix iszero(a) : a undefined'
        z = (a%nrow/=-1.and..not.associated(a%elements))
    end function

    function isdefmat(a) result(d)
        type(matrix),intent(in):: a; logical:: d
        integer:: i; call delabove(i)
        d = (a%nrow/=-1)
    end function


end module dirac_matrix_operations_m
