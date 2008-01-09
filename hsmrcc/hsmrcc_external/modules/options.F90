!--------------------------------------------------------------------
!
!          Module with options stored in logical arrays
!
!--------------------------------------------------------------------
module Options

! ... note: source/hsmrmain.c contains similar parametes, but these do not overlap

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

 subroutine print_Options

  ! use IOModule
!
! ... very simple print out ...

 ! print *,'UOUT=',UOUT

 ! write(UOUT,*) ' OPTION  [GENINT] =',OPTION(GENINT)
 ! write(UOUT,*) ' OPTION [PRNTINF] =',OPTION(PRNTINF)
 ! write(UOUT,*)

  print *,'print_Options ....'
  print *,' OPTION  [GENINT] =',OPTION(GENINT)
  print *,' OPTION [PRNTINF] =',OPTION(PRNTINF)
  print *,' OPTION [MAXIT] =',OPTION(MAXIT)
  print *,' OPTION [DIRACp] =',OPTION(DIRACp)

 end subroutine print_Options

end module Options

