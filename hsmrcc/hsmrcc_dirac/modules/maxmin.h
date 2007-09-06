function MAXMIN(i,j)

   implicit none
   integer, intent(in) :: i, j
   integer :: MAXMIN
   integer :: mx, mn

   if (i==j) stop 'Illegal call to MAXMIN'
   mx = MAX(i,j)
   mn = MIN(i,j)
   MAXMIN = ((mx-2)*(mx-1))/2 + mn

   return

end function MAXMIN
