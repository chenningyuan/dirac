!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module Options
!--------------------------------------------------------------------
!
!          Module with options stored in logical arrays
!
! ... note: source/hsmrmain.c contains similar parametes, but these do not overlap
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
   integer, parameter :: DIRACp   = 18  ! ... can not be "DIRAC"

   ! ... total number of integer options ...
   integer, parameter :: NUMOPTS  = 18

!    logical arrays ... 
   logical :: OPTION(NUMOPTS) = .FALSE.
   logical :: OPTPRV(NUMOPTS) = .FALSE.

 contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine print_Options
!--------------------------------------------------------------------------
!
!
!                  ... very simple print out ...
!
!
!--------------------------------------------------------------------------
 use Parameters

 call QEnter('print_Options')

 call Header('*** Printing out OPTIONs ***',-1)

 write(UOUT,'(/,2x,a)') 'print_Options ....'
 write(UOUT,*) ' OPTION  [GENINT]   =',OPTION(GENINT)
 write(UOUT,*) ' OPTION  [PRNTINF]  =',OPTION(PRNTINF)
 write(UOUT,*) ' OPTION  [MAXIT]    =',OPTION(MAXIT)
 write(UOUT,*) ' OPTION  [DIRACp]   =',OPTION(DIRACp)

 !call Quit('print_Options')

 call QExit('print_Options')
 end subroutine print_Options

end module Options

