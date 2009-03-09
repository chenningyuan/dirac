!
!   common block - huge array, shoudl not exceed 2GB
!
      real*8   gen(nbig)                              
      integer igen(nbig)     

      common //    gen
      common /i/   igen

!     common // gen(nbig)
!     common /i/ igen(nbig)

