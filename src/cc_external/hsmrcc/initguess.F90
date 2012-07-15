! This produces Initial Guess
! for HSMRC iteration.

subroutine InitGuess

   use Parameters
   use Options
   use IOModule
   use TAmps
   use HEFFModule
   implicit none
   integer :: vac
   integer :: i, j, ij, eq
   integer :: StartBufInd, EndBufInd
   double precision :: IntFHP(NPMO*NHMO)
   integer         , dimension(:), pointer :: VIND
   double precision, dimension(:), pointer :: VINT

   if ((.not.OPTION(RESTART)).and.(.not.OPTION(RDGUESS))) then

      ! For a moment Initial Guess will be zero guess
      ! Write it into UTAMP

      TOLD = dble(0)
      do vac = 1, NSTATE
         write(UTAMP,rec=vac) TOLD
      enddo

      ! The following code writes initial guess of V or F/Energy Denominators
      ! For a moment we will skip this code, as the guess is generated
      ! in second iteration, it will just save one iteration.
      ! vaccum independent part of initial
      ! cluster amplitudes in TNEW.
      ! This will be modified to produce
      ! vaccum dependent initial guess.

!     TNEW = dble(0)
!     call ReadFile(VHHPP)
!     EndBufInd = 0
!     ij = 0
!     do i = 2, NHMO
!     do j = 1, (i-1)
!        ij = ij + 1
!        StartBufInd = EndBufInd + 1
!        EndBufInd   = EndBufInd + NZIntRec(ij)
!        if ( EndBufInd < StartBufInd ) cycle
!        VINT => IntBuf(StartBufInd:EndBufInd)
!        VIND => IndBuf(StartBufInd:EndBufInd)
!        eq = TIND(ij)
!        TNEW(eq+VIND) = VINT
!     enddo
!     enddo
!     do vac = 1, NSTATE
!        TOLD = TNEW
!        read(UFHP,rec=vac) IntFHP
!        TOLD(1:NSING) = TOLD(1:NSING) + IntFHP
!        call EnDenom(TOLD,NTSIZE,vac)
!        call MakeZero(TOLD,NTSIZE,vac)
!        write(UTAMP,rec=vac) TOLD
!     enddo

   endif

   ! Copy whatever initial guess is in UTAMP into UTOLD
   ! Before starting the iteration.
   ! In a restart type of calculation, whatever is in
   ! TAMP from previous run is copied into TOLD

   do vac = 1, NSTATE
      read (UTAMP,rec=vac) TOLD
      write(UTOLD,rec=vac) TOLD
   enddo

   ! Read HEFF, this is usually HEF0
   ! if restart then it is HEFF in the previous file
   ! this is used for convergence test.

   rewind(UHEFF)
   read(UHEFF) HEFF

   return

!contains
!include "tind.h"

end subroutine InitGuess
