function TIND(ij)

   integer, intent(in) :: ij
   integer :: TIND

   TIND = (ij-1)*NPMOSQ + NSING

   return

end function TIND
