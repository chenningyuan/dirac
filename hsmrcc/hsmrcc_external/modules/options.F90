!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module Options
!--------------------------------------------------------------------
!
!        Module with options stored in logical arrays
!
! ... note: source/hsmrmain.c contains parametes with almost all identical names, 
!     but these do not overlap
!
!--------------------------------------------------------------------

   implicit none
   integer, parameter :: GENINT   =  1
   integer, parameter :: PRNTINF  =  2
   integer, parameter :: VERBOSE  =  3
   integer, parameter :: RESTART  =  4
   integer, parameter :: TTYRUN   =  5
   integer, parameter :: AVGVAL   =  6
   integer, parameter :: REPTIME  =  7
   integer, parameter :: QUAD     =  8
   integer, parameter :: CONVTAMP =  9
   integer, parameter :: NOHSMR   = 10
   integer, parameter :: MAXIT    = 11
   integer, parameter :: INTTOL   = 12
   integer, parameter :: CNVTOL   = 13
   integer, parameter :: SHFTCORE = 14
   integer, parameter :: SHFTVIRT = 15
   integer, parameter :: NOAPND   = 16
   integer, parameter :: RDGUESS  = 17
   integer, parameter :: HSCC_DIRAC_MODE   = 18  
   integer, parameter :: GMSCHK   = 19  
   integer, parameter :: DIRAC_MODE   = 20

   ! ... total number of integer options - has to be updated...
   integer, parameter :: NUMOPTS  = 20

!    default initialization of logical arrays ... 
   logical :: OPTION(NUMOPTS) = .FALSE.
   logical :: OPTPRV(NUMOPTS) = .FALSE.

 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine set_dirac_mode
!
! tell HSMRCC code that some it's routines are used from within DIRAC
!
#include "priunit.h"
   OPTION(DIRAC_MODE)=.true.
   write(LUPRI,'(2x,a,l1)')   & 
   'DIRAC reusing HSMRCC routines - OPTION(DIRAC_MODE) set to:',OPTION(DIRAC_MODE)
end subroutine set_dirac_mode

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_where

#include "priunit.h"


  if (OPTION(HSCC_DIRAC_MODE)) then
       write(LUPRI,'(2x,a)') '*** in HSMRCC part: DIRAC mode is ON (DIRAC-HSMRCC choice)  ***'
  else if (.not.OPTION(HSCC_DIRAC_MODE).and..not.OPTION(DIRAC_MODE)) then
       write(LUPRI,'(2x,a)') '*** in HSMRCC part: GAMESS mode is ON (GAMESS-HSMRCC choice) ***'
  else if (OPTION(DIRAC_MODE)) then
       write(LUPRI,'(2x,a)') '*** We are in DIRAC preparing files for HSMRCC ***'
  endif

end subroutine print_where

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine print_Options
!--------------------------------------------------------------------------
!
!           ... simple print out of basic options in DIRAC ...
!
!  called from  source/sethsmropts.F90
!
!--------------------------------------------------------------------------
! use Parameters

 call Header('*** Printing out OPTIONs ***',-1)

! write(UOUT,'(/,2x,a)') 'print_Options ....'
 write(*,'(/,2x,a)') 'print_Options ....'
! write(UOUT,*) ' OPTION  [GENINT]   =',OPTION(GENINT)
 write(*,*) ' OPTION  [GENINT]   =',OPTION(GENINT)
 write(*,*) ' OPTION  [PRNTINF]  =',OPTION(PRNTINF)
 write(*,*) ' OPTION  [MAXIT]    =',OPTION(MAXIT)
 write(*,*) ' OPTION  [HSCC_DIRAC_MODE]   =',OPTION(HSCC_DIRAC_MODE)
 write(*,*) ' OPTION  [GMSCHK]   =',OPTION(GMSCHK)
! write(*,*) '  UOUT=',UOUT

 end subroutine print_Options

end module Options
