!-----------------------------------------
! Module with variables previously defined
! at waio.h include file (COMMON /FLNAME/, /POINTR/, /WAIO_STAT/)
!-----------------------------------------
!  Define common blocks and other data for word-addressable I/O
!  The length of records, measured in integer words.
   module waio_mod
   implicit none
   integer*8, parameter :: NWORDS_ON_RECORD8 = 4096
   integer, parameter :: NWORDS_ON_RECORD = 4096
   integer, parameter :: MAX_WAIO_FILES = 100
   real*8 :: WAIO_WCALL, WAIO_RCALL, WAIO_WTIME, WAIO_RTIME, WAIO_WNLEN, WAIO_RNLEN
   integer :: IWRIT_WAIO(MAX_WAIO_FILES)
   character(len=10) :: DN(MAX_WAIO_FILES)
   end module waio_mod
