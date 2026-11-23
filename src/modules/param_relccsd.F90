!-----------------------------------------
!  Module with parameters previously defined
!  at src/relccsd/param.inc include file
!-----------------------------------------
   module param_relccsd_mod
   implicit none
   real*8, parameter     :: ACCUR = 1.D-15
   real*8, parameter     :: AR0   = 0.0D00
   real*8, parameter     :: AR1   = 1.0D00
   real*8, parameter     :: ARP3  = 1.0D00/3.0D00
   complex*16, parameter :: A0    = (0.0D00,0.0D00)
   complex*16, parameter :: A1    = (1.0D00,0.0D00)
   complex*16, parameter :: AP5   = (0.50D00,0.0D00)
   complex*16, parameter :: A2    = (2.0D00,0.0D00)
   complex*16, parameter :: A1I   = (0.0D00,1.0D00)
   logical, parameter    :: LTR   = .TRUE.
   logical, parameter    :: LFA   = .FALSE.
   real*8, parameter     :: D0    = 0.0D00
   real*8, parameter     :: D1    = 1.0D00
   real*8, parameter     :: DP5   = 0.50D00
   end module param_relccsd_mod
