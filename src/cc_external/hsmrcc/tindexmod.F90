module tindexmod

   use Parameters
   implicit none
   interface TINDEX
      module procedure t1index
      module procedure t2index
   end interface

contains

function t1index(i,a)

   implicit none
   integer, intent(in) :: i, a
   integer :: t1index

   if ( (i<1).or.(i>NHMO) ) write(*,*) 'error'
   if ( (a<1).or.(a>NPMO) ) write(*,*) 'error'

   t1index = (i-1)*NPMO+a

   return

end function t1index

function t2index(i,j,a,b)

   implicit none
   integer, intent(in) :: i, j, a, b
   integer :: t2index
   integer :: ijab_equal

   if ( (i<1).or.(i>NHMO) ) write(*,*) 'error'
   if ( (a<1).or.(a>NPMO) ) write(*,*) 'error'
   if ( (j<1).or.(j>NHMO) ) write(*,*) 'error'
   if ( (b<1).or.(b>NPMO) ) write(*,*) 'error'
   ijab_equal = 1
   if ((i.eq.j).or.(a.eq.b)) ijab_equal = 0

!  t2index = isign(1,i-j)*isign(1,a-b)*(i.ne.j)*(a.ne.b) *                 &
   t2index = isign(1,i-j)*isign(1,a-b)*ijab_equal*                         &
             ( ((max(i,j)*(max(i,j)-3))/2 + min(i,j))*(NPMO*(NPMO-1)/2) +  &
             ((max(a,b)-2)*(max(a,b)-1))/2 + min(a,b) + NHMO*NPMO )

   return

end function t2index

end module tindexmod
