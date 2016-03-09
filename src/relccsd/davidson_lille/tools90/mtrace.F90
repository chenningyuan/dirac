module mtrace
!
!***********************************************************************
!          ----- generic module for calculating the trace -----
!          ----- of the product of two matrices           -----
!=======================================================================
!     contains:
!       function caltra(a,b,n)  generic subroutine
!       subroutine caltr_c8c8(a,b,n) : for A(C*8) times B(C*8)
!       subroutine caltr_c8r8(a,b,n) : for A(C*8) times B(R*8)
!       subroutine caltr_r8c8(a,b,n) : for A(R*8) times B(C*8)
!       subroutine caltr_r8r8(a,b,n) : for A(R*8) times B(R*8)
!
!       These routines are aimed to compute the matrix element of an
!       operator between two states A and B: Op_ab = <A|Op|B>
!       If we have |A> = sum_i C_ia |i> and |B> = sum_j C_jb |j>
!       then we have
!              Op_ab = sum_ij conjg(C_ia)*C_jb <i|Op|j>
!                    = Sum_ij DAB_ij Op_ij
!                    = Tr DAB*Op
!       with DAB_ij = conjg(C_ia)*C_jb
!
!       NOTE the definition of DAB.
!-----------------------------------------------------------------------
!   usage:
!   ------
!   val = caltra(a,b,n)
!
!      A(??) = density matrix (INPUT)
!      B(??) = operator matrix (INPUT)
!      N(I4) = order of the matrices (INPUT)
!
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
use kindvar
interface caltra

   module procedure caltr_c8c8, caltr_c8r8, caltr_r8c8, caltr_r8r8
!
!    the actual routines are:
!
!   function caltr_c8c8(a,b,n) result(value)
!      integer, intent(in) :: n
!      complex (kind=K_C8), intent(in), dimension (:,:) :: a,b
!      complex (kind=K_C8) :: value
!   end function caltr_c8c8
!   function caltr_c8r8(a,b,n) result(value)
!      integer, intent(in) :: n
!      double precision, intent(in), dimension (:,:) :: b
!      complex (kind=K_C8), intent(in), dimension (:,:) :: a
!      complex (kind=K_C8) :: value
!   end function caltr_c8r8
!   function caltr_r8c8(a,b,n) result(value)
!      integer, intent(in) :: n
!      double precision, intent(in), dimension (:,:) :: a
!      complex (kind=K_C8), intent(in), dimension (:,:) :: b
!      complex (kind=K_C8) :: value
!   end function caltr_c8r8
!   function caltr_r8r8(a,b,n) result(value)
!      integer, intent(in) :: n
!      double precision, intent(in), dimension (:,:) :: a,b
!      double precision :: value
!   end function caltr_r8r8
end interface
contains
function caltr_c8c8(a,b,n) result(value)
!
!***********************************************************************
!     ----- compute the trace of A*B for square matrices -----
!=======================================================================
!        A(R8) = A matrix (INPUT)
!        B(R8) = B matrix (INPUT)
!        N(I4) = second dimension of A and first dimension of B (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  implicit none
!
!     ----- parameters -----
!
  double precision, parameter :: ZERO=0.0D+00,TWO=2.0D+00
!
!     ----- arguments -----
!
  integer, intent(in) :: n
  complex (kind=K_C8), intent(in), dimension (:,:) :: a,b
  complex (kind=K_C8) :: value
!
!     ----- local variables -----
!
  integer :: i,j
!
!     +++++
!     start
!     +++++
!
  value = dcmplx(ZERO,ZERO)
  do j=1,n
     do i=1,n
        value = value+a(i,j)*b(i,j)
     enddo
  enddo
  return
end function caltr_c8c8
function caltr_c8r8(a,b,n) result(value)
!
!***********************************************************************
!     ----- compute the trace of A*B for square matrices -----
!=======================================================================
!        A(R8) = A matrix (INPUT)
!        B(R8) = B matrix (INPUT)
!        N(I4) = second dimension of A and first dimension of B (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  implicit none
!
!     ----- parameters -----
!
  double precision, parameter :: ZERO=0.0D+00,TWO=2.0D+00
!
!     ----- arguments -----
!
  integer, intent(in) :: n
  double precision, intent(in), dimension (:,:) :: b
  complex (kind=K_C8), intent(in), dimension (:,:) :: a
  complex (kind=K_C8) :: value
!
!     ----- local variables -----
!
  integer :: i,j
!
!     +++++
!     start
!     +++++
!
  value = dcmplx(ZERO,ZERO)
  do j=1,n
     do i=1,n
        value = value+a(i,j)*b(i,j)
     enddo
  enddo
  return
end function caltr_c8r8
function caltr_r8c8(a,b,n) result(value)
!
!***********************************************************************
!     ----- compute the trace of A*B for square matrices -----
!=======================================================================
!        A(R8) = A matrix (INPUT)
!        B(R8) = B matrix (INPUT)
!        N(I4) = second dimension of A and first dimension of B (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  implicit none
!
!     ----- parameters -----
!
  double precision, parameter :: ZERO=0.0D+00,TWO=2.0D+00
!
!     ----- arguments -----
!
  integer, intent(in) :: n
  double precision, intent(in), dimension (:,:) :: a
  complex (kind=K_C8), intent(in), dimension (:,:) :: b
  complex (kind=K_C8) :: value
!
!     ----- local variables -----
!
  integer :: i,j
!
!     +++++
!     start
!     +++++
!
  value = dcmplx(ZERO,ZERO)
  do j=1,n
     do i=1,n
        value = value+a(i,j)*b(i,j)
     enddo
  enddo
  return
end function caltr_r8c8
function caltr_r8r8(a,b,n) result(value)
!
!***********************************************************************
!     ----- compute the trace of A*B for square matrices -----
!=======================================================================
!        A(R8) = A matrix (INPUT)
!        B(R8) = B matrix (INPUT)
!        N(I4) = second dimension of A and first dimension of B (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  implicit none
!
!     ----- parameters -----
!
  double precision, parameter :: ZERO=0.0D+00,TWO=2.0D+00
!
!     ----- arguments -----
!
  integer, intent(in) :: n
  double precision, intent(in), dimension (:,:) :: a,b
  double precision :: value
!
!     ----- local variables -----
!
  integer :: i,j
!
!     +++++
!     start
!     +++++
!
  value = ZERO
  do j=1,n
     do i=1,n
        value = value+a(i,j)*b(i,j)
     enddo
  enddo
  return
end function caltr_r8r8
end module mtrace
