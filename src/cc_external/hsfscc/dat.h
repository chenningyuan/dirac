!
! ... iopt(1) is print level output; for more see README
!
!  miro: how about to use own parameter for 200 ?
!

      integer      norb,ic(200),icv(200),idum2,         iopt(40)
      real*8                                   eps(200)

      common /dat/ norb,ic     ,icv     ,idum2,eps     ,iopt

!      common /dat/ norb,ic(200),icv(200),idum2,eps(200),iopt(40)
