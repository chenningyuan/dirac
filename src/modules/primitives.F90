!-----------------------------------------
! Module with variables previously defined
! at aovec.h include file
!-----------------------------------------
   module aovec_module
!     MXCONT = maximum number of contracted 
!            functions in an AO input block 
!            (in "MOLECULE.INP" input file)
!     MXAOVC = used for MXCONT in a few places
!
!     IF you change MXCONT you should make sure
!     MAXPRD in cbirea.h is at least the same value and
!     then rebuild the program using the command "make".
!
   implicit none
   integer, parameter :: MXCONT = 40
   integer, parameter :: MXAOVC = MXCONT
   end module aovec_module



!-----------------------------------------
! Module with variables previously defined
! at maxorb.h and primit.h include files
!-----------------------------------------
   module primitives_module
   use aovec_module
   implicit none
   private

   integer, parameter :: MXSHEL      = 15000
   integer, parameter :: MXPRIM      = 15000
   integer, parameter :: MXCORB      = 40000

   integer, parameter :: MXSHEL_FIT  = 1
   integer, parameter :: MXPRIM_FIT  = 1
   integer, parameter :: MXCORB_FIT  = 1

   public :: MXSHEL, MXPRIM, MXCORB
   public :: MXSHEL_FIT, MXPRIM_FIT, MXCORB_FIT
   public :: init_primitives

   public :: PRIEXP, PRICCF, PRICRX, PRICRY, PRICRZ
   public :: NPSHEL, NPLSH, NPSSH

   real(8), allocatable :: PRIEXP(:)
   real(8), allocatable :: PRICCF(:,:)
   real(8), allocatable :: PRICRX(:)
   real(8), allocatable :: PRICRY(:)
   real(8), allocatable :: PRICRZ(:)

   integer :: NPSHEL = 0
   integer :: NPLSH  = 0
   integer :: NPSSH  = 0

   contains

      subroutine init_primitives()
         ! Allocates with defined sizes
         if (.not. allocated(PRIEXP)) then
            allocate( PRIEXP(MXPRIM) )
            allocate( PRICCF(MXPRIM, MXCONT) )
            allocate( PRICRX(MXPRIM) )
            allocate( PRICRY(MXPRIM) )
            allocate( PRICRZ(MXPRIM) )
         end if
      end subroutine init_primitives

   end module primitives_module

