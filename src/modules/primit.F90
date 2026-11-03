!-----------------------------------------
! Module with variables previously defined
! at primit.h include file
!-----------------------------------------
!
!   PRIEXP  - exponents of primitive shells
!   PRICCF  - normalized contraction coefficients
!   PRICRX  - x-coordinate of center
!   PRICRY  - y-coordinate of center
!   PRICRZ  - z-coordinate of center
!
   module primit_mod

   use kinds_mod,  only: dp
   use aovec_mod,  only: MXCONT

   implicit none

   private

   public :: ensure_primitives, resize_primitives
   public :: PRIEXP, PRICCF, PRICRX, PRICRY, PRICRZ
   public :: NPSHEL, NPLSH, NPSSH

   real(dp), allocatable :: PRIEXP(:)
   real(dp), allocatable :: PRICCF(:,:)
   real(dp), allocatable :: PRICRX(:)
   real(dp), allocatable :: PRICRY(:)
   real(dp), allocatable :: PRICRZ(:)

   integer :: NPSHEL = 0
   integer :: NPLSH  = 0
   integer :: NPSSH  = 0

   contains

   subroutine resize_primitives(n)
!     Give the arrays exactly n primitives, keeping the first min(n, old)
!     ones and zeroing the rest (as the old COMMON /PRIMIT/ was).
      integer, intent(in) :: n
      real(dp), allocatable :: ccf(:,:)
      integer :: nold, nkeep

      nold = 0
      if (allocated(PRIEXP)) nold = size(PRIEXP)
      if (allocated(PRIEXP) .and. n == nold) return
      nkeep = min(n, nold)

      allocate(ccf(n, MXCONT), source=0.0_dp)
      if (nkeep > 0) ccf(1:nkeep, :) = PRICCF(1:nkeep, :)
      call move_alloc(ccf, PRICCF)
      call resize(PRIEXP)
      call resize(PRICRX)
      call resize(PRICRY)
      call resize(PRICRZ)

   contains

      subroutine resize(a)
         real(dp), allocatable, intent(inout) :: a(:)
         real(dp), allocatable :: tmp(:)
         allocate(tmp(n), source=0.0_dp)
         if (nkeep > 0) tmp(1:nkeep) = a(1:nkeep)
         call move_alloc(tmp, a)
      end subroutine resize

   end subroutine resize_primitives

   subroutine ensure_primitives(n)
!     Make room for at least n primitives, doubling so a caller can ask
!     for one more at a time; trim with resize_primitives when done.
      integer, intent(in) :: n
      integer :: nold

      nold = 0
      if (allocated(PRIEXP)) nold = size(PRIEXP)
      if (n > nold) call resize_primitives(max(n, 2*nold))
   end subroutine ensure_primitives

   end module primit_mod
