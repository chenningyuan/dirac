module TimingMod

   use Parameters
   use IOModule
   implicit none
   real, dimension(2,NSUBR) :: TimeTaken

contains

subroutine StartTiming

   implicit none

   TimeTaken = real(0)

   return

end subroutine StartTiming

function TOTTime()

   implicit none
   real :: TOTTime

   TOTTime = SUM(TimeTaken)

   return

end function TOTTime

function USRTime()

   implicit none
   real :: USRTime

   USRTime = SUM(TimeTaken(1,:))

   return

end function USRTime

function SYSTime()

   implicit none
   real :: SYSTime

   SYSTime = SUM(TimeTaken(2,:))

   return

end function SYSTime

subroutine PrintTiming

   implicit none
   real :: sys, usr, tot

   usr = USRTime()
   sys = SYSTime()
   tot = TOTTime()

   write(UOUT,*)
   write(UOUT,7) usr, sys, tot
   write(UOUT,*)
   write(UOUT,8) '+-----------------------+'
   write(UOUT,8) '|  SUBROUTINE TIMINGS   |'
   write(UOUT,8) '+-----------------------+'
   write(UOUT,9) '| V2H2P  => ',SUM(TimeTaken(:,(V2H2P ))),' |'
   write(UOUT,9) '| VPPPP  => ',SUM(TimeTaken(:,(VPPPP ))),' |'
   write(UOUT,9) '| VHHHH  => ',SUM(TimeTaken(:,(VHHHH ))),' |'
   write(UOUT,9) '| VHHPP  => ',SUM(TimeTaken(:,(VHHPP ))),' |'
   write(UOUT,9) '| V2HPP  => ',SUM(TimeTaken(:,(V2HPP ))),' |'
   write(UOUT,9) '| V2PHH  => ',SUM(TimeTaken(:,(V2PHH ))),' |'
   write(UOUT,9) '| VPHHP  => ',SUM(TimeTaken(:,(VPHHP ))),' |'
   write(UOUT,9) '| VPH2P  => ',SUM(TimeTaken(:,(VPH2P ))),' |'
   write(UOUT,9) '| VPHPP  => ',SUM(TimeTaken(:,(VPHPP ))),' |'
   write(UOUT,9) '| V2HHP  => ',SUM(TimeTaken(:,(V2HHP ))),' |'
   write(UOUT,9) '| VHHHP  => ',SUM(TimeTaken(:,(VHHHP ))),' |'
   write(UOUT,9) '| FHH    => ',SUM(TimeTaken(:,(FHH   ))),' |'
   write(UOUT,9) '| FHP    => ',SUM(TimeTaken(:,(FHP   ))),' |'
   write(UOUT,9) '| FPP    => ',SUM(TimeTaken(:,(FPP   ))),' |'
   write(UOUT,9) '| HEFFEC => ',SUM(TimeTaken(:,(HEFFEC))),' |'
   write(UOUT,9) '| RENORM => ',SUM(TimeTaken(:,(RENORM))),' |'
   write(UOUT,9) '| INTGEN => ',SUM(TimeTaken(:,(INTGEN))),' |'
   write(UOUT,8) '+-----------------------+'
   write(UOUT,9) '| TOTAL  => ',SUM(TimeTaken            ),' |'
   write(UOUT,8) '+-----------------------+'
   write(UOUT,*)

   7 format('EXECUTION TIME(secs) ==> USR : ',F12.4,' SYS : ',F12.4,' TOT : ',F12.4)
   8 format(A25)
   9 format(A11,F12.4,A2)

   return

end subroutine PrintTiming

end module TimingMod
