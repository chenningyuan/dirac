function TIND(ij)
!----------------------------------------------------------------------
! ... inline integer function inserted in various places of the code
! ... depends on 'outside' variables NPMOSQ, NSING
!----------------------------------------------------------------------

   integer, intent(in) :: ij
   integer :: TIND

   TIND = (ij-1)*NPMOSQ + NSING

   return

end function TIND
