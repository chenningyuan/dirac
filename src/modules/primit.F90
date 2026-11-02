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
   use maxorb_mod, only: MXPRIM
   use aovec_mod,  only: MXCONT

   implicit none

   private

   public :: init_primitives
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

   subroutine init_primitives()
!     Same size as the old COMMON /PRIMIT/ (MXPRIM primitives).
!     Zeroed as the COMMON was: entries beyond the ones set by ORBPRO
!     may be read.
      if (allocated(PRIEXP)) return
      allocate(PRIEXP(MXPRIM))
      allocate(PRICCF(MXPRIM, MXCONT))
      allocate(PRICRX(MXPRIM))
      allocate(PRICRY(MXPRIM))
      allocate(PRICRZ(MXPRIM))
      PRIEXP = 0.0_dp
      PRICCF = 0.0_dp
      PRICRX = 0.0_dp
      PRICRY = 0.0_dp
      PRICRZ = 0.0_dp
   end subroutine init_primitives

   end module primit_mod
