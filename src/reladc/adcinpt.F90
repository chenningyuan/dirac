!-----------------------------------------
! Module with variables previously defined
! at adcinpt.inc include file (COMMON /CCADCI/, /CCADCR/, /CCADCL/, /LANCZI/, /LANCZR/, /LANCZL/)
!-----------------------------------------
!-----For the ADC module
   module adcinpt_mod
   implicit none
   real*8 :: VCONV, ADCTHR, SIPPRNT, DIPPRNT
   integer :: ADCLEVEL, ADCPRINT(8), SIPREPS(32), DIPREPS(32), SIPITER, DIPITER, &
              SIPEIGV(32), DIPEIGV(32), LANCMEM
   logical :: DOSIPS, DODIPS, READQKL, DOCONST, DOADCPOP, DOFULL, DOLANC, DOINCORE
   end module adcinpt_mod
