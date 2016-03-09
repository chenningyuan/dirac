module mr3vec
!
!***********************************************************************
!     ----- this module contains routine to do operations on -----
!     ----- vectors in R^3                                   -----
!=======================================================================
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
!     ----- parameters -----
!
   double precision, parameter :: ONE=1.0D+00
   private :: ONE
   contains
   subroutine jpvec(xa,xb,u,ru)
!
!***********************************************************************
!            ----- build the unit vector A->B -----
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
   double precision, intent(out) :: ru
   double precision, intent(in), dimension (:) :: xa,xb
   double precision, intent(out), dimension (:) :: u
!
!     ----- local variables -----
!
   double precision :: r
!
!     +++++
!     start
!     +++++
!
   u(1:3) = xb(1:3)-xa(1:3)
   ru = sqrt(sum(u(1:3)**2))
   r = ONE/ru
   u(1:3) = r*u(1:3)
   return
   end subroutine jpvec
   subroutine vprod(a,b,c,r)
!
!***********************************************************************
!     ----- compute the vector product A x B = C -----
!           vector C is normalized. r = sin((A,B))
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   implicit none
!
!     ----- parameter -----
!
!     ----- arguments -----
!
   double precision, intent(out) :: r
   double precision, intent(in), dimension (:) :: a,b
   double precision, intent(out), dimension (:) :: c
!
!     ----- local variables -----
!
   double precision :: dum
!
!     +++++
!     start
!     +++++
!
   c(1) = a(2)*b(3)-a(3)*b(2)
   c(2) = a(3)*b(1)-a(1)*b(3)
   c(3) = a(1)*b(2)-a(2)*b(1)
   dum = sum(c(1:3)**2)
   r = sqrt(dum)
   if(dum.gt.1.0D-08) c(1:3) = c(1:3)*(ONE/r)
   return
   end subroutine vprod
end module mr3vec
