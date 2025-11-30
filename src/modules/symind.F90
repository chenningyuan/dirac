!-----------------------------------------
! Module with variables previously defined
! at symind.h include file
!-----------------------------------------
   module symind_mod
!
   use maxaqn_mod, only: MXAQN
!
   implicit none
   INTEGER :: INDFA(8,MXAQN), INDFB(8,MXAQN)
   INTEGER :: ISOFRA(8), ISOFRB(8)
!
   end module symind_mod
