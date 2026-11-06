!-----------------------------------------
! Module with variables previously defined
! at symind.h include file (COMMON /SYMIND/)
!-----------------------------------------
   module symind_mod
   use maxaqn_mod, only: MXAQN
   implicit none
   integer :: INDFA(8,MXAQN), INDFB(8,MXAQN), ISOFRA(8), ISOFRB(8)
   end module symind_mod
