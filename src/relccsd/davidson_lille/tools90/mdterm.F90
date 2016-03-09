module mdterm
!
!***********************************************************************
!          ----- Module for linear algebra -----
!=======================================================================
!     contains:
!      subroutine dterm(n,a,d,maxa) : compute determinant(a)
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
    use kindvar
    use mconst
    implicit none
!
!   ---- parameters -----
!
    interface dterm
        module procedure dterm_r8,dterm_c8
    end interface dterm
!
!   ----- global variables -----
!
    contains
    subroutine dterm_r8(n,a,d,maxa)
!
! **********************************************************************
!     calcule le determinant de la matrice a (real version)
! **********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n
      double precision, intent(out) :: d,maxa
      double precision, intent(inout), dimension (:,:) :: a
!
!     ----- local variables -----
!
      integer :: l,kp,k,j
      double precision :: z,ratio
!
!     +++++
!     start
!     +++++
!
      d = ONE
      do l=1,n
         kp = 0
         z = ZERO
         do k=l,n
            if(z.lt.abs(a(k,l))) then
               z = abs(a(k,l))
               kp = k
            else if(z.gt.abs(a(k,l))) then
               exit
            endif
         enddo
         if(l.lt.kp) then
            do j=l,n
               z = a(l,j)
               a(l,j) = a(kp,j)
               a(kp,j) = z
            enddo
            d = -d
         endif
         if(l.lt.n) then
            do k=l+1,n
               if(a(k,l).ne.ZERO) then
                  ratio = -a(k,l)/a(l,l)
                  a(k,l+1:n) = a(k,l+1:n)+ratio*a(l,l+1:n)
               endif
            enddo
         endif
      enddo
      k = mod(n,4)
      if(k.gt.0) then
         do l=1,k
            d = d*a(l,l)
         enddo
      endif
      do l=k+1,n,4
         d = (a(l,l)*a(l+1,l+1))*(a(l+2,l+2)*a(l+3,l+3))*d
      enddo
      maxa = ZERO
      do k=1,n
         maxa = max(maxa,abs(a(k,k)))
      enddo
!     if(abs(d).lt.TENM14) d = ZERO
      return
   end subroutine dterm_r8
    subroutine dterm_c8(n,a,d,maxa)
!
! **********************************************************************
!     calcule le determinant de la matrice a (complex version)
! **********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
      integer, intent(in) :: n
      double precision, intent(out) :: maxa
      complex (kind=K_C8), intent(out) :: d
      complex (kind=K_C8), intent(inout), dimension (:,:) :: a
!
!     ----- local variables -----
!
      integer :: l,kp,k,j
      double precision :: zr
      complex (kind=K_C8) :: z,ratio
!
!     +++++
!     start
!     +++++
!
      d = CONE
      do l=1,n
         kp = 0
         zr = ZERO
         do k=l,n
            if(zr.lt.abs(a(k,l))) then
               zr = abs(a(k,l))
               kp = k
!           else if(zr.gt.abs(a(k,l))) then
!              exit
            endif
         enddo
         if(l.lt.kp) then
            do j=l,n
               z = a(l,j)
               a(l,j) = a(kp,j)
               a(kp,j) = z
            enddo
            d = -d
         endif
         if(l.lt.n) then
            do k=l+1,n
               if(abs(a(k,l)).ne.ZERO) then
                  ratio = -a(k,l)/a(l,l)
                  a(k,l+1:n) = a(k,l+1:n)+ratio*a(l,l+1:n)
               endif
            enddo
         endif
      enddo
      k = mod(n,4)
      if(k.gt.0) then
         do l=1,k
            d = d*a(l,l)
         enddo
      endif
      do l=k+1,n,4
         d = (a(l,l)*a(l+1,l+1))*(a(l+2,l+2)*a(l+3,l+3))*d
      enddo
      maxa = ZERO
      do k=1,n
         maxa = max(maxa,abs(a(k,k)))
      enddo
!     if(abs(d).lt.TENM14) d = ZERO
      return
   end subroutine dterm_c8
end module mdterm
