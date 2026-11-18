!-----------------------------------------
! Module with variables previously defined
! at dcbcls.h include file (COMMON /XCBCLS/, /CLSLBL/)
!-----------------------------------------
!     CLASS LABELS
!     CLSINT - keyword for a class of integrals
!     CLSCMB - indicates what blocks to compute : LL,SL,LS,SS
!     CLSCAL - indicates whether class is active, that is to be calculated
   module dcbcls_mod
   implicit none
   integer, parameter :: MAXCLS = 45
   integer :: NPRPCLS, IORDCL(MAXCLS)
   logical :: CLSCAL(MAXCLS)
   character(len=4) :: CLSCMB(MAXCLS)
   character(len=7) :: CLSINT(MAXCLS)
   end module dcbcls_mod
