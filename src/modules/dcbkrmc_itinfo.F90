!-----------------------------------------
! Module with variables previously defined
! at dcbkrmc_itinfo.h include file (COMMON /IITINFO/, /RITINFO/)
!-----------------------------------------
!     i     IINFO      DINFO
!    ---   -------    -------
!     01    ITMAC      EMY
!     02    ITMIC      EACTIV
!     03    NZRED      EMCSCF
!     04    NCRED      DEPRED
!     05    NERED      DEACT
!     06    NPRED      RATIO
!     07               BETA
!     08               STPLNG
!     09               RTRUST
!     10               GAMMA
!     11               GCINRM ( GNORM(1) )
!     12               GEENRM ( GNORM(2) )
!     13               GEPNRM ( GNORM(3) )
!     14               GOBNRM ( GNORM(4) )
!     15               GTOTNRM( GNORM(5) )
!     16               TIMMAC
!     17               TIMMIC
!     18               TIMITR
!     19               CI STPLNG
!     20               E STPLNG
!     21               P STPLNG
!     22               LEVEL SHIFT
   module dcbkrmc_itinfo_mod
   implicit none
   integer, parameter :: LDINFO = 22
   integer, parameter :: LIINFO = 6
   real*8 :: DINFO(LDINFO)
   integer :: IINFO(LIINFO)
   end module dcbkrmc_itinfo_mod
