!-----------------------------------------
! Module with variables previously defined
! at src/include/moldip.h and src/abacus/moldip.h include files
! (same COMMON; the abacus one held only the first five members) (COMMON /MOLDIP/)
!-----------------------------------------
   module moldip_mod
   use mxcent_mod, only: MXCENT, MXCOOR
   implicit none
   real*8 :: DIP0(3), DIP1(3,MXCOOR), DIPFLT(3,MXCOOR), POLARS(3,3), POLFLT(3,3), &
             DIPLF0(3), QAPT(MXCENT)
   end module moldip_mod
