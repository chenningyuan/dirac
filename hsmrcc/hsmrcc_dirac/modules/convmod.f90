module ConvMod

   use Parameters
   use Options
   use IOmodule
   use HEFFModule
   implicit none
   integer          :: CurIter = 0
   integer          :: MAXITER
   logical          :: CONVERGED = .FALSE.
   double precision :: CONVTOL
   double precision :: TAMPConv
   double precision :: HEFFConv
   double precision, allocatable :: TDIFF(:)
   logical, private :: Init = .FALSE.

contains

subroutine InitConvMod

   implicit none

   if (Init) then
      write(UOUT,*) 'ConvMod already Initialized'
      call flush(UOUT)
      call abort
   endif
   Init = .TRUE.

   if (.not.allocated(TDIFF)) allocate(TDIFF(NSTATE))
   return

end subroutine InitConvMod

subroutine CheckConverg

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'ConvMod Not Initialized'
      call flush(UOUT)
      call abort
   endif

   TAMPConv = SUM(TDIFF)/dble(NSTATE*NTSIZE)
   HEFFConv = SUM(DABS(HEFF-HEFP))/dble(NSTATE*NSTATE)
   if (OPTION(CONVTAMP)) then
      if (TAMPConv <= CONVTOL) CONVERGED = .TRUE.
   else
      if (HEFFConv <= CONVTOL) CONVERGED = .TRUE.
   endif

   return

end subroutine CheckConverg

subroutine EndConvMod

   implicit none

   if (.not.Init) then
      write(UOUT,*) 'ConvMod Not Initialized'
      call flush(UOUT)
      call abort
   endif
   Init = .FALSE.

   if (allocated(TDIFF)) deallocate(TDIFF)

   return

end subroutine EndConvMod

function ExitHSMRIter()

   implicit none
   logical :: ExitHSMRIter

   ExitHSMRIter = (CONVERGED).or.(CurIter>=MAXITER)

   return

end function ExitHSMRIter

end module ConvMod
