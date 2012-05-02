

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine dirprg_hscc
!
! called from DIRAC/main
!
use Parameters
use IOModule
use Options
use Dirac_HSMRCC

implicit none
#include "priunit.h"

if (nohsccdirac) then
  write(LUPRI,'(2x,a)')  'dirprg_hscc - returning'
  return
else
  write(LUPRI,'(2x,a)')  'dirprg_hscc - continuing with HSCC.'
endif

write(LUPRI,'(2x,a,i3)') 'dirprg_hscc: get_iprint()=',get_iprint()

call FLSHFO(LUPRI)

! ... initialization of basic variables to make HSCC running in DIRAC mode
call Init_HSMRCC_BasicVars

! ... call subroutine of the program ...
call HSMRCC

end subroutine dirprg_hscc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Init_HSMRCC_BasicVars
!
! ... initialize only few basis variables for HSMRCC which are otherwise initialized by
! ... C-part (hsmrcc.c) 
!
! ... later OPTIONS[xxx] are read from DIRAC prepared files
!
use Options
use Parameters

implicit none
#include "priunit.h"
character*12 :: filename="Dirac_HSMRCC"

write(LUPRI,*) "Init_HSMRCC_BasicVars: "

OPTION(TTYRUN)=.true.
OPTION(GMSCHK)=.false.

OPTION(DIRAC_MODE)=.false.
OPTION(GENINT)  = .FALSE.

OPTION(HSCC_DIRAC_MODE)=.true. ! tell that this is HSCC from DIRAC

call unset_dirac_hscc_interface_mode

call print_options; 
call print_where

!call call_dirac_sethsmrfiles

!call SetHSMRfiles(filename,"",12,0)

end subroutine Init_HSMRCC_BasicVars

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSMRCC
!---------------------------------------------------------------------------
!
!
!    This subroutine is called from main
!   program written in C after calling all SetHSMR routines.
!
!
!---------------------------------------------------------------------------
   use Parameters
   use IOModule
   use Options
   use Shift
   use ConvMod
   use IntGenMod
   use HEFFModule
   use VINFMod
   implicit none

!   integer :: i1, i2, i3, i4, i5
!   double precision, pointer :: VINT(:) ! 
!   integer, pointer :: VIND(:)

   if (get_iprint() > 5 ) then
     call HEADER('Output from HSMRCC',-1)
   endif

   ! ... output stream (=UOUT)
   call set_output_stream()

   ! ... check of GAMESS integrals
   if ( OPTION(GMSCHK) ) then
     write(UOUT,*) 'only check of GAMES integr'
     call check_gms_integrals
     return
   endif

   ! ... serie of initializations ...
   call InitHSMRBase

   ! ... Generation of Integrals ...

   if (OPTION(GENINT)) then

      ! ... Gamess integrals driver (modules/intgenmod.f90)
      if (.not.OPTION(HSCC_DIRAC_MODE)) then
         call DrvIntGen       ! ... call Gamess integral driver
      else
         call DrvIntGen_DIRAC ! ... call DIRAC integrals based driver
      endif

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
