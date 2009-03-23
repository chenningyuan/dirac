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

      integer       nij,ij(ndij),ije(0:ndij),iii,iiv
      common /cij1/ nij,ij,      ije,        iii,iiv

!     common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

!
!
!   Note: overlapping with the include file "cij1.h":
!
!   common /cij1/ nij, ij (ndij),ije (0:ndij),iiii,iivi
!   
!

