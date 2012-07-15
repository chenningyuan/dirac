!>
!> ... prints out convergence info for each iteration
!>
subroutine PrintInfo

   use Parameters
   use IOModule
   use TimingMod
   use HEFFModule
   use ConvMod
   implicit none
   integer :: i

!  print convergence information

   if (OPTION(CONVTAMP)) then
      write(UOUT,9001) CurIter, HEFFConv, TAMPConv, CONVTOL, TOTTime()
   else
      write(UOUT,9002) CurIter, HEFFConv, TAMPConv, CONVTOL, TOTTime()
   endif

9001 format(/,'ITR #',I3,2X,'ConvHEFF=',F12.9,2X,'ConvTAMP=',F12.9,2X,'ConvTol=',F12.9,' (TAMP) ',F8.2,' (secs)',/)
9002 format(/,'ITR #',I3,2X,'ConvHEFF=',F12.9,2X,'ConvTAMP=',F12.9,2X,'ConvTol=',F12.9,' (HEFF) ',F8.2,' (secs)',/)

!  print effective hamiltonian

   call PrintHEFF

   return

end subroutine PrintInfo
