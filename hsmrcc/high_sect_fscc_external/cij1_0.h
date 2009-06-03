!
! ...common block...
!
!  nij - number of different ij-pairs in <i(1)j(2)|k(1)l(2)>
!
!  ij(:) - array where ij-vector is located
!
!  ije(:) - number of different integrals with the same ij-pair, ije(1)- beginning, ije(2)-end
!     ije(0)=0
!

      integer       nij,ij(ndij),ije(0:ndij),iiii,iivi
      common /cij1/ nij,ij,      ije,        iiii,iivi

! to replace: common /cij1/ nij,ij(ndij),ije(0:ndij),iiii,iivi
!
!   Note: /cij1/ block is overlapping with the include file "cij1.h"
! 
! there are many variations of /cij1/ in the code ! cleanup needed
!
!

