!-----------------------------------------
! Module with variables previously defined
! at shells.h include file
!-----------------------------------------
   module shells_mod
   implicit none
   private
!
   !---------------------------------------------------------------
   ! Public
   !---------------------------------------------------------------
   public :: CENT, NHKT, KHKT, KCKT, ISTBAO, NUCO, JSTRT, NSTRT
   public :: NCENT, NRCO, NUMCF, NBCH, KSTRT, LCLASS, IPTSHL
   public :: NUMCFT, KMAX, NLRGSH, NSMLSH, NORBS
   public :: SEGM, SPHR
   public :: ensure_shells, resize_shells
!
   !*****************************************************************************
   !     Data on the AO-block associated with a given shell(ISHELL):
   !       NUCO    - number of uncontracted functions
   !       NRCO    - number of contracted functions
   !     Data on a given shell in an AO-block:
   !       NCENT   - index of symmetry independent center
   !       NUMCF   - index of shell in AO-block
   !       NBCH    - index of block in AO-vector
   !       ISTBAO  - stabiliser: basic sym. op. that do not move center
   !       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
   !       KHKT    - number of spherical (Cartesian) components
   !       KCKT    - number of Cartesian components
   !       LCLASS  - class: large component (1), small or Huckel(2), density fitting (0)
   !       CENT    - coordinates of center
   !       JSTRT   - offset for the shell
   !       SEGM    - segmented contraction
   !*****************************************************************************
   real(8), allocatable :: CENT(:, :, :)   ! (shells, 3, 2)
   integer, allocatable :: NHKT(:), KHKT(:), KCKT(:), ISTBAO(:)
   integer, allocatable :: NUCO(:), JSTRT(:), NSTRT(:)
   integer, allocatable :: NCENT(:), NRCO(:), NUMCF(:)
   integer, allocatable :: NBCH(:), KSTRT(:), LCLASS(:)
   integer, allocatable :: IPTSHL(:), NUMCFT(:)
   logical, allocatable :: SEGM(:), SPHR(:)
   integer :: KMAX = 0, NLRGSH = 0, NSMLSH = 0, NORBS = 0
!
   contains

   subroutine resize_shells(n)
!     Give the arrays exactly n shells, keeping the first min(n, old)
!     ones and zeroing the rest (as the old COMMON /SHELLS/ was).
      integer, intent(in) :: n
      real(8), allocatable :: c(:, :, :)
      integer :: nold, nkeep

      nold = 0
      if (allocated(NHKT)) nold = size(NHKT)
      if (allocated(NHKT) .and. n == nold) return
      nkeep = min(n, nold)

      allocate(c(n, 3, 2), source=0.0d0)
      if (nkeep > 0) c(1:nkeep, :, :) = CENT(1:nkeep, :, :)
      call move_alloc(c, CENT)
      call resize_i(NHKT)
      call resize_i(KHKT)
      call resize_i(KCKT)
      call resize_i(ISTBAO)
      call resize_i(NUCO)
      call resize_i(JSTRT)
      call resize_i(NSTRT)
      call resize_i(NCENT)
      call resize_i(NRCO)
      call resize_i(NUMCF)
      call resize_i(NBCH)
      call resize_i(KSTRT)
      call resize_i(LCLASS)
      call resize_i(IPTSHL)
      call resize_i(NUMCFT)
      call resize_l(SEGM)
      call resize_l(SPHR)

   contains

      subroutine resize_i(a)
         integer, allocatable, intent(inout) :: a(:)
         integer, allocatable :: tmp(:)
         allocate(tmp(n), source=0)
         if (nkeep > 0) tmp(1:nkeep) = a(1:nkeep)
         call move_alloc(tmp, a)
      end subroutine resize_i

      subroutine resize_l(a)
         logical, allocatable, intent(inout) :: a(:)
         logical, allocatable :: tmp(:)
         allocate(tmp(n), source=.false.)
         if (nkeep > 0) tmp(1:nkeep) = a(1:nkeep)
         call move_alloc(tmp, a)
      end subroutine resize_l

   end subroutine resize_shells

   subroutine ensure_shells(n)
!     Make room for at least n shells, doubling so a caller can ask
!     for one more at a time; trim with resize_shells when done.
      integer, intent(in) :: n
      integer :: nold

      nold = 0
      if (allocated(NHKT)) nold = size(NHKT)
      if (n > nold) call resize_shells(max(n, 2*nold))
   end subroutine ensure_shells

   end module shells_mod
