module mphase
!
!***********************************************************************
!    ----- generic module for changing the phase of a vector    -----
!    ----- so that it's maximum coefficient is positive or real -----
!=======================================================================
!     contains:
!       subroutine cphase(a,m,n)  generic subroutine
!       subroutine phase_c81(a,m,n)  subroutine for A(C*8)
!       subroutine phase_c82(ar,ai,m,n)  subroutine for A(C*8)=(AR,AI)
!       subroutine phase_r8(a,m,n)  subroutine for A(R*8)
!-----------------------------------------------------------------------
!   usage:
!   ------
!   call cphase(a,m,n)  or call cphase(ar,ai,m,n)
!
!           A(C8) = matrix (INPUT)
!   AR(R8),AI(R8) = real and imaginary par of A (INPUT)
!           M(I4) = number of rows dimension of A -or AR, AI- (INPUT)
!           N(I4) = number of columns (INPUT)
!
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
use kindvar
use mercod
   implicit none
   double precision, parameter :: ONE=1.0D+00,ZERO=0.0D+00,TENP6=1.0D+06
   private ZERO,ONE,TENP6
interface cphase
   module procedure phase_c81, phase_c82, phase_r8
end interface
   integer (kind=K_I4), private, dimension (1) :: j2
contains
subroutine phase_r8(a,m,n)
!
!***********************************************************************
!     ----- -----
!=======================================================================
!        A(C8) = A matrix (INPUT)
!        M(I4) = first dimension of A (INPUT)
!        N(I4) = second dimension of A (INPUT)
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
!     ----- arguments -----
!
   integer (kind=K_I4), intent(in) :: m,n
   double precision, intent(inout), dimension (:,:) :: a
!
!     ----- local variables -----
!
   integer (kind=K_I4) :: j,i,jmax
   double precision :: val
!
!     +++++
!     start
!     +++++
!
!$omp do private(j2) schedule(dynamic)
   do j=1,n
      j2 = maxloc(int(abs(a(1:m,j))*TENP6))
      jmax = j2(1)
      if(a(jmax,j).lt.ZERO) a(1:m,j) = -a(1:m,j)
   enddo
!$omp end do
   return
end subroutine phase_r8
subroutine phase_c81(a,m,n)
!
!***********************************************************************
!     ----- -----
!=======================================================================
!        A(C8) = A matrix (INPUT)
!        M(I4) = first dimension of A (INPUT)
!        N(I4) = second dimension of A (INPUT)
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
!     ----- arguments -----
!
  integer (kind=K_I4), intent(in) :: m,n
  complex (kind=K_C8), intent(inout), dimension (:,:) :: a
!
!     ----- local variables -----
!
   integer (kind=K_I4) :: j
   complex (kind=K_C8) :: cval
!
!     +++++
!     start
!     +++++
!
!$omp do private(j2) schedule(dynamic)
   do j=1,n
!     j2 = maxloc(abs(a(1:m,j)))
      j2 = maxloc(int(abs(a(1:m,j))*TENP6))
      cval = conjg(a(j2(1),j))/abs(a(j2(1),j))
      a(1:m,j) = a(1:m,j)*cval
   enddo
!$omp end do
   return
end subroutine phase_c81
subroutine phase_c82(vr,vi,m,n)
!
!***********************************************************************
!     ----- -----
!=======================================================================
!        A(C8) = A matrix (INPUT)
!     ADIM(I4) = first dimension of A (INPUT)
!        B(R8) = B matrix (INPUT)
!     BDIM(I4) = first dimension of B (INPUT)
!        N(I4) = second dimension of A and first dimension of B (INPUT)
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
!     ----- arguments -----
!
  integer (kind=K_I4), intent(in) :: m,n
  double precision, intent(inout), dimension (:,:) :: vr,vi
!
!     ----- local variables -----
!
   integer (kind=K_I4) :: i,j,ios,nth,ith,irval,jmax
   integer (kind=K_I4), dimension (1) :: j2
   double precision :: crval,cival,ar,ai,dum
   integer (kind=K_I4), allocatable, dimension (:,:) :: tmp
!
!     +++++
!     start
!     +++++
!
   nth = 1
!omp nth = omp_get_thread_num()
   allocate(tmp(m,nth),stat=ios)
   if(ios.ne.0) then
      write(6,fmt='(5x,a)') 'CPHASE ERROR: allocation tmp failed'
      call pricod(ios)
      stop
   endif
   ith = 1
!omp ith = omp_get_num_thread()
!$omp do private(j2,crval,cival,ar,ai) schedule(dynamic)
   do j=1,n
      tmp(1:m,ith) = int(sqrt(vr(1:m,j)**2+vi(1:m,j)**2)*TENP6)
      j2(1) = 1
      irval = tmp(1,ith)
      do i=2,m
         if(tmp(i,ith).gt.irval) then
            irval = tmp(i,ith)
            j2(1) = i
         endif
      enddo
      jmax = j2(1)
      dum = ONE/sqrt(vr(jmax,j)**2+vi(jmax,j)**2)
!     crval = vr(j2(1),j)/tmp(j2(1),ith)
!     cival = vi(j2(1),j)/tmp(j2(1),ith)
      crval = vr(jmax,j)*dum
      cival = vi(jmax,j)*dum
!     if(crval.lt.1.0D+00) then
!        crval = -crval
!        cival = -cival
!     endif
      do i=1,m
         ar = vr(i,j)
         ai = vi(i,j)
         vr(i,j) = ar*crval+ai*cival
         vi(i,j) = ai*crval-ar*cival
      enddo
   enddo
!$omp end do
   deallocate(tmp)
   return
end subroutine phase_c82
end module mphase
