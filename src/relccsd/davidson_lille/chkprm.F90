function chkprm(c,e,m,n) result(value)
!
!***********************************************************************
!     ----- check any permutation of eigenvectors between -----
!     ----- iterations i and i+1                          -----
!=======================================================================
!-----------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  use kindvar
  use mdebug
  implicit none
!
!     ----- parameters -----
!
!
!     ----- arguments -----
!
  integer, intent(in) :: m,n
  double precision, intent(inout), dimension (:) :: e
  complex (kind=K_C8), intent(inout), dimension (:,:) :: c
  logical :: value
!
!
!     ----- allocatable variables -----
!
  complex (kind=K_C8),dimension(:), allocatable :: tmp
!     ----- local variables -----
!
  integer :: i,j,k,ios
  double precision :: cmax,dum
  logical :: test
!
!     +++++
!     start
!     +++++
!
#ifdef _MYDEBUG_
  call subent('CHKPRM')
#endif
  test = .FALSE.
!
!     ----- Check permutation between the eigenstates -----
!
  allocate(tmp(m),stat=ios)
  do i=1,n
     k = i
     cmax = dble(abs(c(i,1)))
     do j=2,m
        if(dble(abs(c(i,j))).gt.cmax) then
           cmax = dble(abs(c(i,j)))
           k = j
        endif
     enddo
     if(k.ne.i) then
        tmp(1:m) = c(1:m,i)
        c(1:m,i) = c(1:m,k)
        c(1:m,k) = tmp(1:m)
        dum = e(i)
        e(i) = e(k)
        e(k) = dum
     endif
  enddo
  value = test
  deallocate(tmp)
#ifdef _MYDEBUG_
  call subret('CHKPRM',1)
#endif
  return
end function chkprm
