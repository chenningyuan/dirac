
! radovan.bast@uit.no (2009-06-01):
! this module is the interface between matrix_genop and DIRAC
! please mind that the naming is common with Dalton
! so some things make only sense in the larger picture (Dalton and DIRAC)

module matrix_operations

  implicit none

  public matrix
  public mat_init
  public mat_free
  public mat_mul
  public mat_zero
  public mat_scal
  public mat_trans
  public mat_dotproduct
  public mat_trab
  public mat_daxpy

  private

  type matrix
!   mark uninitialized matrices with nrow = ncol = -1
    integer             :: nrow = -1
    integer             :: ncol = -1

    real(8),    pointer :: elms(:)
    real(8),    pointer :: elmsb(:)          !not used in DIRAC
    complex(8), pointer :: celms(:)          !not used in DIRAC
    complex(8), pointer :: celmsb(:)         !not used in DIRAC

!   to mark complex elements
    logical             :: complex = .false. !not used in DIRAC

    integer,    pointer :: iaux(:)
    real(8),    pointer :: raux(:)

!   DIRAC specific
    integer             :: ih_sym        = 1
    integer             :: tr_sym        = 1
    integer             :: pg_sym        = 1
    integer             :: fermi_pointer = 1
  end type

#include "dgroup.h"

contains

  subroutine mat_init(A, nrow, ncol)

!   ----------------------------------------------------------------------------
    type(matrix)        :: A
    integer, intent(in) :: nrow, ncol
!   ----------------------------------------------------------------------------

    nullify(A%elms)
    nullify(A%iaux)
    nullify(A%raux)

    A%nrow = nrow
    A%ncol = ncol

    allocate(A%elms(A%nrow*A%ncol*nz*merge(2, 1, A%complex)))

  end subroutine

  subroutine mat_free(A)

!   ----------------------------------------------------------------------------
    type(matrix) :: A
!   ----------------------------------------------------------------------------

    if (.not. associated(A%elms)) then
      call quit('error in mat_free: memory previously released')
    end if

    deallocate(A%elms)
    nullify(A%elms)

    if (associated(A%iaux)) deallocate(A%iaux)
    if (associated(A%raux)) deallocate(A%raux)

  end subroutine

  subroutine mat_zero(A)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: A
!   ----------------------------------------------------------------------------

    A%elms = 0.0d0

  end subroutine

  subroutine mat_mul(A, B, transa, transb, alpha, beta, C)

!   C      = alpha*A*B + beta*C
!   C_{mn} = ... A_{mk} B_{kn} ...
!   transa/transb: 'T' or 't' (transposed), 'N' or 'n' (normal)
!   ----------------------------------------------------------------------------
    type(matrix), intent(in)    :: A, B
    character,    intent(in)    :: transa, transb
    real(8),      intent(in)    :: alpha, beta
    type(matrix), intent(inout) :: C
!   ----------------------------------------------------------------------------
    integer                     :: m, n, k, nza, nzb, nzc
    character                   :: ta, tb
    integer                     :: iq(4) = (/1, 2, 3, 4/)
!   ----------------------------------------------------------------------------

    if      (transa == 'n' .or. transa == 'N') then
      m = A%nrow
      k = A%ncol
    else if (transa == 't' .or. transa == 'T') then
      m = A%ncol
      k = A%nrow
    else
      call quit('unknown trans possibility in mat_mul')
    end if

    if      (transb == 'n' .or. transb == 'N') then
      n = B%ncol
      if (B%nrow /= k) call quit('dimensions do not match in mat_mul')
    else if (transb == 't' .or. transb == 'T') then
      n = B%nrow
      if (B%ncol /= k) call quit('dimensions do not match in mat_mul')
    else
      call quit('unknown trans possibility in mat_mul')
    end if

    nza = nz
    nzb = nz
    if (nbsym > 1) then
      call iqpack(ipqtoq(1, A%pg_sym-1), nza, &
                  ipqtoq(1, B%pg_sym-1), nzb, &
                  iq                   , nzc)
    end if

    ta = 'N'
    tb = 'N'
    if (transa == 't') ta = 'H'
    if (transb == 't') tb = 'H'
    if (transa == 'T') ta = 'H'
    if (transb == 'T') tb = 'H'

!   C%nrow       = m
!   C%ncol       = n
    C%pg_sym     = 1 + ieor(A%pg_sym-1, B%pg_sym-1)

    call qgemm(m, n, k,                                                    &
               alpha,                                                      &
               ta, 'N', ipqtoq(1, A%pg_sym-1), A%elms, A%nrow, A%ncol, nz, &
               tb, 'N', ipqtoq(1, B%pg_sym-1), B%elms, B%nrow, B%ncol, nz, &
               beta,                                                       &
                        iq                   , C%elms, C%nrow, C%ncol, nz)

  end subroutine

  subroutine mat_scal(alpha, A)

!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: alpha
    type(matrix), intent(inout) :: A
!   ----------------------------------------------------------------------------

    call dscal(A%nrow*A%ncol*nz, alpha, A%elms, 1)

  end subroutine

  subroutine mat_trans(A, B)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: A
    type(matrix)             :: B
!   ----------------------------------------------------------------------------
    integer                  :: i, j, i_dim, j_dim, iz
    real(8)                  :: q
!   ----------------------------------------------------------------------------

    if (B%nrow /= A%ncol .or. B%ncol /= A%nrow) then
      call quit('wrong dimensions in mat_trans')
    end if

    i_dim = A%nrow
    j_dim = A%ncol

    iz = 1
    do iz = 1, nz

      q = 1.0d0
      if (ipqtoq(iz, A%pg_sym - 1) > 1) q = -1.0d0

      do j = 1, j_dim
        do i = 1, i_dim
              B%elms((iz - 1)*i_dim*j_dim + i_dim*(i - 1) + j) &
          = q*A%elms((iz - 1)*i_dim*j_dim + i_dim*(j - 1) + i)
        end do
      end do
    end do

!   B%ih_sym        = A%ih_sym
!   B%tr_sym        = A%tr_sym
!   B%pg_sym        = A%pg_sym
!   B%fermi_pointer = A%fermi_pointer

  end subroutine

  function mat_dotproduct(A, B) result(r)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: A, B
    real(8)                  :: r
!   ----------------------------------------------------------------------------
    real(8),      external   :: ddot
!   ----------------------------------------------------------------------------

    if (A%nrow*A%ncol /= B%nrow*B%ncol) then
      call quit('wrong dimensions in mat_dotproduct')
    end if

    r = 0.0d0

    if (ieor(A%pg_sym-1, B%pg_sym-1) == 0) then
      r = ddot(A%nrow*A%ncol*nz, A%elms, 1, B%elms, 1)
    end if

  end function

  function mat_trab(A, B) result(r)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: A, B
    real(8)                  :: r
!   ----------------------------------------------------------------------------
    real(8),      external   :: ddot
    real(8)                  :: q
    integer                  :: i_dim, j_dim, j, iz
!   ----------------------------------------------------------------------------

    if (A%ncol /= B%nrow .or. A%nrow /= B%ncol) then
      call quit('wrong dimensions in mat_trab')
    end if

    r = 0.0d0

    i_dim = A%nrow
    j_dim = A%ncol

    do iz = 1, nz

      q = 1.0d0
!     if (ipqtoq(iz, A%pg_sym - 1) > 1) q = -1.0d0
      if (ipqtoq(iz, B%pg_sym - 1) > 1) q = -1.0d0

      do j = 1, j_dim
        r = r + q*ddot(i_dim, A%elms((iz - 1)*i_dim*j_dim + i_dim*(j - 1) + 1), 1, &
                              B%elms((iz - 1)*i_dim*j_dim + j), j_dim)
      enddo
    end do

  end function

  subroutine mat_daxpy(alpha, X, Y)

!   Y = alpha*X + Y
!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: alpha
    type(matrix), intent(in)    :: X
    type(matrix), intent(inout) :: Y
!   ----------------------------------------------------------------------------

    if (X%nrow /= Y%nrow .or. X%ncol /= Y%ncol) then
      call quit('wrong dimensions in mat_daxpy')
    endif

    call daxpy(X%nrow*X%ncol*nz, alpha, X%elms, 1, Y%elms,1)

  end subroutine

end module
