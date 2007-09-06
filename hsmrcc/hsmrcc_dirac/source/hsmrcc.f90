!
!
! This is called from main
! program written in C after
! calling SetHSMR routines.
!
!

subroutine HSMRCC

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

   ! ... Basic Initialization of HSMR Run ...

   call InitHSMRBase

   ! ... Generation of Integrals ...

   if (OPTION(GENINT)) then

    ! ... integral driver (modules/intgenmod.f90)
      call DrvIntGen

      call UpdateVINF

      if (OPTION(NOHSMR)) then
         call EndHSMR
         return
      endif
   endif

   ! ... This point is reached if HSMR    ...
   ! ... iterations have to be performed  ...

   ! ... Initialize HSMR Iterations ...

   call InitHSMRIter

   call PrintRunInfo

   !
   ! ... Begin HSMRCC iterations ...
   !

   CurIter = 0

   ! ... Diagonalize the effective Hamiltonian         ...
   ! ... This is usually HEF0 unless it is restart run ...
   ! ... It is used to relate the states in            ...
   ! ... Successive iterations                         ...

   call DiagHEFF

   call MapRoots

   call PrintInfo

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

   return

contains

!
! Very simple subroutine - Exchange indexes
!
subroutine ExchIJ(I,J)

   implicit none
   integer, intent(inout) :: I, J
   integer :: K ! local variable 

   K = I
   I = J
   J = K

   return

end subroutine ExchIJ

end subroutine HSMRCC
