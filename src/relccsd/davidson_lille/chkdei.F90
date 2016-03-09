function chkdei(c,e,m,n) result(value)
!
!***********************************************************************
!     ----- check any permutation of eigenvectors between -----
!     ----- iterations i and i+1                          -----
!=======================================================================
!     C (R8) = eigenvectors (INPUT/OUTPUT)
!     E (R8) = eigenvalues (INPUT/OUTPUT)
!     M (I4) = size of the current model space (INPUT)
!     N (I4) = total number of requested eigenvectors (INPUT)
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  use kindvar
#ifdef _MYDEBUG_
  use mdebug
#endif
  implicit none
!
!     ----- parameters -----
!
!
!     ----- arguments -----
!
  integer, intent(in) :: m,n
  double precision, intent(inout), dimension (:) :: e
  double precision, intent(inout), dimension (:,:) :: c
  logical :: value
!
!     ----- allocatable variables -----
!
  double precision, dimension(:), allocatable :: tmp
!
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
  call subent('CHKDEI')
#endif
  test = .FALSE.
!
!     ----- Check permutation between the eigenstates -----
!
  allocate(tmp(m),stat=ios)
  do i=1,n
     k = i
     cmax = abs(c(i,1))
     do j=2,m
        if(abs(c(i,j)).gt.cmax) then
           cmax = abs(c(i,j))
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
  call subret('CHKDEI',1)
#endif
  return
end function chkdei
