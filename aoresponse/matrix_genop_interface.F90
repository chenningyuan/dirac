
! radovan.bast@uit.no (2009-02-18):
! this module is the interface between matrix_genop and DIRAC

module matrix_genop_interface

  implicit none

  public get_nz
  public get_ipqtoq
  public matrix
  public mat_mul
  public mat_init_full
  public mat_free_nonzero

  private

  type matrix
    integer          :: nrow          = -1
    integer          :: ncol          = -1
    integer          :: ih_sym        = 1
    integer          :: tr_sym        = 1
    integer          :: pg_sym        = 1
    integer          :: fermi_pointer = 1
    real(8), pointer :: elms(:)
    integer, pointer :: iaux(:)
    real(8), pointer :: raux(:)
  end type

  integer :: nz = 1
! nz is introduced to accommodate complex and quaternion algebra for DIRAC
! it could be part of the matrix type
! for the moment it is a parameter

contains

  function get_nz() result(nz_up)

#include "dgroup.h"

!   ----------------------------------------------------------------------------
    integer :: nz_up
!   ----------------------------------------------------------------------------

    nz_up = nz

  end function

  function get_ipqtoq(i, j) result(k)

#include "dgroup.h"

!   ----------------------------------------------------------------------------
    integer, intent(in) :: i, j
    integer             :: k
!   ----------------------------------------------------------------------------

    k = ipqtoq(i, j)

  end function

  subroutine mat_dense_mul(a,b,transa,transb,alpha,beta,c)
     TYPE(Matrix), intent(IN) :: a, b
     character, intent(in)    :: transa, transb
     character                :: ta, tb
     REAL(8), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c
     INTEGER                  :: m,n,k
#ifdef PRG_DIRAC
    integer      :: iq(4)
#endif /* PRG_DIRAC */  

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
#ifndef PRG_DIRAC
     call DGEMM(transa,transb,m,n,k,alpha,&
               &a%elms,a%nrow,b%elms,b%nrow,beta,c%elms,c%nrow)
#else

  nz = get_nz()
      call iqpack(get_ipqtoq(1, A%pg_sym-1), nz, &
                  get_ipqtoq(1, B%pg_sym-1), nz, &
                  iq                       , nz)


    ta = 'N'
    tb = 'N'
    if (transa == 't') ta = 'H'
    if (transb == 't') tb = 'H'
    if (transa == 'T') ta = 'H'
    if (transb == 'T') tb = 'H'

!   C%nrow       = A%nrow
!   C%ncol       = B%ncol
    C%pg_sym     = 1 + ieor(A%pg_sym-1, B%pg_sym-1)
      call qgemm(m, n, k,                1.0d0,  &
                 ta, 'N', get_ipqtoq(1, A%pg_sym-1), A%elms, A%nrow, A%ncol, nz, &
                 tb, 'N', get_ipqtoq(1, B%pg_sym-1), B%elms, B%nrow, B%ncol, nz, &
                 0.0d0,    iq                      , C%elms, C%nrow, C%ncol, nz)


#endif /* PRG_DIRAC */  


  end subroutine mat_dense_mul

      SUBROUTINE mat_mul(a, b, transa, transb, alpha, beta, c)
         !c = alpha*ab + beta*c
         !transa = 'T'/'t' - transposed, 'N'/'n' - normal
         TYPE(Matrix), intent(IN) :: a, b
         character, intent(in)    :: transa, transb
         REAL(8), INTENT(IN)  :: alpha, beta
         TYPE(Matrix), intent(inout):: c
         integer :: ak, bk, ci, cj
    !    if (INFO_TIME_MAT) CALL TIMER('START ',mat_TSTR,mat_TEN)
      !  no_of_matmuls = no_of_matmuls + 1
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
             call mat_dense_mul(a,b,transa, transb,alpha,beta,c)
      END SUBROUTINE mat_mul

  subroutine mat_init_full(A, nrow, ncol)

!   ----------------------------------------------------------------------------
    type(matrix)            :: A
    integer,     intent(in) :: nrow, ncol
!   ----------------------------------------------------------------------------

    nullify(A%elms)
    nullify(A%iaux)
    nullify(A%raux)

    A%nrow = nrow
    A%ncol = ncol

  nz = get_nz()
    allocate(A%elms(A%nrow*A%ncol*nz))

  end subroutine

  subroutine mat_free_nonzero(A)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: A
!   ----------------------------------------------------------------------------

    if (associated(A%elms)) then

      deallocate(A%elms)
      nullify(A%elms)
  
      if (associated(A%iaux)) deallocate(A%iaux)
      if (associated(A%raux)) deallocate(A%raux)

    end if

    A%nrow = -1
    A%ncol = -1

  end subroutine

end module
