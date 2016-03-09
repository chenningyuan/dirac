subroutine incord(o,d,nn,mm)
!
!**********************************************************************
!      order vectors according to decreasing o(i)
!----------------------------------------------------------------------
!      o = vector of dimension mm which contains the test values
!      d = array of mm vectors of length nn to be ordered according to
!          o
!     nn = length of the vectors in d
!     mm = number of vectors to be ordered
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: nn,mm
   double precision, intent(inout), dimension(:) :: o
   double precision, intent(inout), dimension(:,:) :: d
!
!     ----- local variables -----
!
   integer :: i,j,k
   double precision :: p,dum
!
!     +++++
!     start
!     +++++
!
   do i=1,mm-1
      p = o(i)
      k = i
      do j=i+1,mm
         if(o(j).lt.p) then
            p = o(j)
            k = j
         endif
      enddo
      if(k.ne.i) then
         o(k) = o(i)
         o(i) = p
         do j=1,nn
            dum = d(j,i)
            d(j,i) = d(j,k)
            d(j,k) = dum
         enddo
      endif
   enddo
   return
end subroutine incord
