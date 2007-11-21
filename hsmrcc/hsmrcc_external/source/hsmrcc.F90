!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSMRCC
!---------------------------------------------------------------------------
!
!
! This is called from main
! program written in C after
! calling SetHSMR routines.
!
!
!---------------------------------------------------------------------------
   use IOModule
   use Options
   use Shift
   use ConvMod
   use IntGenMod
   use HEFFModule
   use VINFMod
   implicit none
   integer :: i1, i2, i3, i4, i5
   double precision, pointer :: VINT(:)
   integer, pointer :: VIND(:)

   call QEnter('HSMRCC')

   ! ... Basic Initialization of HSMR Run ...

   if (get_iprint() > 5 ) then
     call HEADER('Output HSMRCC',-1)
   endif

   call InitHSMRBase

   ! ... Generation of Integrals ...

   if (OPTION(GENINT)) then

     ! ... integrals driver (modules/intgenmod.f90)
      call DrvIntGen

      call UpdateVINF

      if (OPTION(NOHSMR)) then
         call EndHSMR
         call QExit('HSMRCC')
         return
      endif

   endif

   ! ... This point is reached if HSMR    ...
   ! ... iterations have to be performed  ...

   ! ... Initialize for HSMRCC Iterations ...

   call InitHSMRIter

   call PrintRunInfo

   !-------------------------------------------------------
   ! ... Begin HSMRCC iterations ...
   !-------------------------------------------------------

   CurIter = 0

   ! ... Diagonalize the effective Hamiltonian         ...
   ! ... This is usually HEF0 unless it is restart run ...
   ! ... It is used to relate the states in            ...
   ! ... Successive iterations                         ...

   call DiagHEFF

   call MapRoots

   call PrintInfo

   ! ... run iterations ...
   do

      CurIter = CurIter + 1

      call StartIter
      call Diagrams
      call DiagHEFF
      call MapRoots
      call CheckConverg
      call PrintInfo

      if ( ExitHSMRIter() ) exit

      call ExchIJ(UTOLD,UTNEW)

      CurIter = CurIter + 1

      call StartIter
      call Diagrams
      call DiagHEFF
      call MapRoots
      call CheckConverg
      call PrintInfo

      if ( ExitHSMRIter() ) exit

      call ExchIJ(UTOLD,UTNEW)

   enddo

   ! ... End of HSMRCC iterations ...

   call EndHSMR

   call QExit('HSMRCC')
   return

contains

subroutine ExchIJ(I,J)
!
! Very simple subroutine - Exchange indexes
!
   implicit none
   integer, intent(inout) :: I, J
   integer :: K ! local variable 

   K = I
   I = J
   J = K

   return

end subroutine ExchIJ

end subroutine HSMRCC

!C  /* Deck mpapv */
      SUBROUTINE MPAPV(N,AP,VIN,VOUT)
!C
!C  28-Apr-1985 hjaaj
!C
!C  VOUT = AP * VIN
!C
!C  where VOUT and VIN are N vectors
!C  and AP is a packed, symmetric N by N matrix.
!C
!
!  ... copied from Dirac/pdpack/arhpack.F ...
!  in order to avoid linking of hsmrcc.exe with the FULL Dirac/pdpack/libpdlib.a library
!
!#include <implicit.h>
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION AP(*),VIN(*),VOUT(*)
      PARAMETER (D0 = 0.0D0)
!C
      KI = 0
      DO 300 I = 1,N
         VINI = VIN(I)
         SUM = D0
         DO 100 J = 1,I-1
            SUM = SUM + AP(KI+J) * VIN(J)
            VOUT(J) = VOUT(J) + AP(KI+J) * VINI
  100    CONTINUE
         KI = KI + I
         VOUT(I) = SUM + AP(KI) * VINI
  300 CONTINUE
!C
      RETURN
      END
