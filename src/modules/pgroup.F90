!-----------------------------------------
! Module with variables previously defined
! at pgroup.h include file (COMMON /PGROUP/)
!-----------------------------------------
   module pgroup_mod
   implicit none
   character(len=3) :: GROUP, REP(0:7)
   ! formerly SAVEd DATA-initialized locals of every unit including pgroup.h (read only)
   character(len=3), parameter :: SYMOP(0:7)  = [' E ', 'Oyz', 'Oxz', 'C2z', 'Oxy', 'C2y', 'C2x', ' i ']
   character(len=3), parameter :: GROUPS(0:7) = ['C1 ', 'C2 ', 'Cs ', 'Ci ', 'D2 ', 'C2v', 'C2h', 'D2h']
   end module pgroup_mod
