subroutine StartIter

   use Parameters
   use HEFFModule
   use TAmps
   use TimingMod
   implicit none
   integer :: vac

   call QEnter('StartIter')

   ! ... Restart the timer for ther iteration ...

   call StartTiming

   ! ... Initialise effective Hamiltonian. ...

   HEFP = HEFF
   HEFF = dble(0)

   ! ... Erase the effects of previous calculations ...

   TNEW = dble(0)
   TOLD = dble(0)

   ! ... copy model space coefficients and energies of previous iterations ...
   ! ... zero the current model space coeffients and energies ...

!  EN_PREV = EN
!  EN = dble(0)
!
!  BCL_PREV = BCL
!  BCR_PREV = BCR
!
!  BCL = dble(0)
!  BCR = dble(0)
!

   call QExit('StartIter')
   return

end subroutine StartIter
