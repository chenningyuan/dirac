!
! common block
!

      integer       nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1
      common /cij1/ nij1,ij1,      ije1,        iii1,iiv1

!   common /cij1/ nij1,ij1(ndij),ije1(0:ndij),iii1,iiv1

!
! Note: overlapping with this block :
!   
!   common /cij1/ nij, ij (ndij),ije (0:ndij),iiii,iivi
!
