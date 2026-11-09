!-----------------------------------------
! Module with variables previously defined
! at cbirea.h include file (COMMON /CBIREA/, /CBIFAM/)
!-----------------------------------------
!     MAXPRD = default value for MAXPRI
   module cbirea_mod
   implicit none
   integer, parameter :: MAXPRD = 35
   integer, parameter :: MXFAMEXP = 100
   real*8 :: ZCMVAL, TOL_SYMADD, FAMEXP(MXFAMEXP,2), FAMPAR(4)
   integer :: IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP(2)
   logical :: BIGVC, GENCON, INPTST, DIRAC, BASIS, PRIBAS, ATOMBA, UNCONT, CNTMAT, LCNTNUUM
   end module cbirea_mod
