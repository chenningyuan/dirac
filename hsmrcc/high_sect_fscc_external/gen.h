!
!   common block - huge array, shoudl not exceed 2GB
!
!  all integral indexes are in igen, their values are in gen (up to 256 basis functions, 32bit model)
!
!
      real*8   gen(nbig)                              
      integer igen(nbig)     

      common //    gen
      common /i/   igen

!     common // gen(nbig)
!     common /i/ igen(nbig)

