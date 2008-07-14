!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSMRCC
!---------------------------------------------------------------------------
!
!
!   This subroutine is called from main
!   program written in C after calling all SetHSMR routines.
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

   ! ... Basic Initialization of HSMR Run ...

   if (get_iprint() > 5 ) then
     call HEADER('Output from HSMRCC',-1)
   endif

   call InitHSMRBase

   ! ... Generation of Integrals ...

   if (OPTION(GENINT)) then

      ! ... integrals driver (modules/intgenmod.f90)
      call DrvIntGen

      ! ... update bunch of data on the file (modules/vinfmod.F90)
      call UpdateVINF

      ! ... leave if no further continuation of HSCC
      if (OPTION(NOHSMR)) then
         call EndHSMR
         return
      endif

   endif

   ! ... This point is reached if HSMR    ...
   ! ... iterations have to be performed  ...

   ! ... Initialize for HSMRCC Iterations ...
   call InitHSMRIter

   ! ... print-out info on the current run
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

   return

contains

subroutine ExchIJ(I,J)
!
!  Very simple subroutine - exchange indexes I,J
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
