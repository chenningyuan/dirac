module matrix_genop

  use kinds_m
! use memory_allocator

  implicit none

#include "dgroup.h"

  type matrix
    integer                :: nrow
    integer                :: ncol
    integer                :: ih_sym        = 1
    integer                :: tr_sym        = 1
    integer                :: pg_sym        = 1
    integer                :: fermi_pointer = 1
    real(kind_dp), pointer :: elms(:)
    logical                :: is_defined    = .false.

      type(Matrix),pointer :: block(:)
      !the indexes hvere the blocks start
      integer,dimension(:,:), pointer :: blockpos!(2,:)
      !room for any auxaliary information
      integer, pointer     :: iaux(:)
      real(8), pointer :: raux(:)

  end type

    interface mat_init
        module procedure mat_init_full
        module procedure mat_init_from
    end interface

contains

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

  subroutine mat_init_full(A, nrow, ncol)

    type(matrix), intent(inout) :: A
    integer,      intent(in)    :: nrow, ncol

    A%nrow       = nrow
    A%ncol       = ncol
    A%is_defined = .true.

    nullify(A%elms)
    allocate(A%elms(A%nrow*A%ncol*nz))

  end subroutine

  subroutine matrix_eq_matrix_backend(A, B)

    type(matrix) :: A, B

    A%nrow          =  B%nrow
    A%ncol          =  B%ncol
    A%ih_sym        =  B%ih_sym
    A%tr_sym        =  B%tr_sym
    A%pg_sym        =  B%pg_sym
    A%fermi_pointer =  B%fermi_pointer
    A%is_defined    =  .true.

    A%elms          => B%elms

  end subroutine

  function real_times_matrix_backend(r, A) result(B)

    type(matrix),  intent(in) :: A
    real(kind_dp), intent(in) :: r
    type(matrix)              :: B

    B%nrow          = A%nrow
    B%ncol          = A%ncol
    B%ih_sym        = A%ih_sym
    B%tr_sym        = A%tr_sym
    B%pg_sym        = A%pg_sym
    B%fermi_pointer = A%fermi_pointer
    B%is_defined    = .true.

    nullify(B%elms)

    if (r /= 0 .and. associated(A%elms)) then
      allocate(B%elms(B%nrow*B%ncol*nz))
      B%elms       = r*A%elms
    end if

  end function

  function matrix_plus_matrix_backend(A, B) result(C)

    type(matrix), intent(in) :: A, B
    type(matrix)             :: C

    if (A%nrow /= B%nrow) call quit('different shapes in C = A + B')
    if (A%ncol /= B%ncol) call quit('different shapes in C = A + B')

    C%nrow          = B%nrow
    C%ncol          = B%ncol
    C%ih_sym        = B%ih_sym
    C%tr_sym        = B%tr_sym
    C%pg_sym        = B%pg_sym
    C%fermi_pointer = B%fermi_pointer
    C%is_defined    = .true.

    nullify(C%elms)

    if (associated(A%elms) .and. associated(B%elms)) then
      allocate(C%elms(C%nrow*C%ncol*nz))
      C%elms = A%elms + B%elms
    else if (associated(A%elms)) then
      C%elms => A%elms
    else if (associated(b%elms)) then
      C%elms => B%elms
    end if

  end function

  function matrix_dot_matrix_backend(A, B) result(r)

    type(matrix), intent(in) :: A, B
    real(kind_dp)            :: r

    if (A%nrow /= B%nrow) call quit('different shapes in r = A**B')
    if (A%ncol /= B%ncol) call quit('different shapes in r = A**B')

    r = 0.0d0

    if (.not.(associated(A%elms) .and. associated(B%elms))) return

    if (ieor(A%pg_sym-1, B%pg_sym-1) == 0) r = 2.0d0*sum(A%elms*B%elms)

  end function

  subroutine print_matrix_backend(A, label, w)

    type(matrix), intent(in)           :: A
    character(*), intent(in), optional :: label
    integer,      intent(in), optional :: w

    write(*, *) 'matrix        = ', label
    write(*, *) 'ih_sym        = ', A%ih_sym
    write(*, *) 'tr_sym        = ', A%tr_sym
    write(*, *) 'pg_sym        = ', A%pg_sym
    write(*, *) 'fermi_pointer = ', A%fermi_pointer

    if (associated(A%elms)) then
      call prqmat(A%elms,                &
                  A%nrow,                &
                  A%ncol,                &
                  A%nrow,                &
                  A%ncol,                &
                  nz,                    &
                  ipqtoq(1, A%pg_sym-1), &
                  6)
    else
      write(*, *) 'all elements zero'
    end if

  end subroutine

  subroutine matrix_times_matrix_backend(A, B, C)

    type(matrix) :: A, B, C
    integer      :: iq(4)

    if (A%ncol /= B%nrow) call quit('incompatible shapes in C = A*B')

    C%nrow       = A%nrow
    C%ncol       = B%ncol
    C%pg_sym     = 1 + ieor(A%pg_sym-1, B%pg_sym-1)
    C%is_defined = .true.

    if (associated(A%elms) .and. associated(B%elms)) then

      allocate(C%elms(C%nrow*C%ncol*nz))

      call iqpack(ipqtoq(1, A%pg_sym-1), nz, &
                  ipqtoq(1, B%pg_sym-1), nz, &
                  iq                   , nz)

      call qgemm(A%nrow, B%ncol, A%ncol, 1.0d0,  &
                 'N', 'N', ipqtoq(1, A%pg_sym-1), A%elms, A%nrow, A%ncol, nz, &
                 'N', 'N', ipqtoq(1, B%pg_sym-1), B%elms, B%nrow, B%ncol, nz, &
                 0.0d0,    iq                   , C%elms, C%nrow, C%ncol, nz)

    end if

  end subroutine

  subroutine dag_matrix_backend(A, B)

    type(matrix)  :: A
    type(matrix)  :: B
    real(kind_dp) :: factor
    integer       :: iz

    B%ncol          = A%nrow
    B%nrow          = A%ncol
    B%ih_sym        = A%ih_sym
    B%tr_sym        = A%tr_sym
    B%pg_sym        = A%pg_sym
    B%fermi_pointer = A%fermi_pointer
    B%is_defined    = .true.

    if (associated(A%elms)) then
      allocate(B%elms(B%nrow*B%ncol*nz))

      do iz = 1, nz
        if (ipqtoq(iz, A%pg_sym-1) > 1) then
          factor = -1.0
        else
          factor =  1.0
        end if

        call transpose_block(factor, A%elms((iz-1)*A%nrow*A%ncol + 1:iz*A%nrow*A%ncol), A%nrow, A%ncol, &
                                     B%elms((iz-1)*B%nrow*B%ncol + 1:iz*B%nrow*B%ncol))

      end do
    endif

  end subroutine

  subroutine transpose_block(factor, block, nrow, ncol, tblock)

    real(kind_dp),    intent(in)    :: factor
    integer,          intent(in)    :: nrow
    integer,          intent(in)    :: ncol
    real(kind_dp),    intent(in)    :: block(nrow*ncol)
    real(kind_dp),    intent(inout) :: tblock(nrow*ncol)

    tblock = reshape(transpose(reshape(block, (/nrow, ncol/))), (/nrow*ncol/))
    call dscal(nrow*ncol, factor, tblock, 1)

  end subroutine

end module
