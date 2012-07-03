function MAXMIN(i,j)

   implicit none
   integer, intent(in) :: i, j
   integer :: MAXMIN
   integer :: mx, mn

!  if (i==j) stop 'Illegal call to MAXMIN'

   if (i==j) then
!    ... better error output ...
     write(6,*) 'error in MAXMIN ! i,j=',i,j
     Call Quit('Illegal call to MAXMIN')
   endif

   mx = MAX(i,j)
   mn = MIN(i,j)
   MAXMIN = ((mx-2)*(mx-1))/2 + mn

   return
end function MAXMIN

