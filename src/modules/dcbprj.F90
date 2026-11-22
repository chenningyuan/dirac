!-----------------------------------------
! Module with variables previously defined
! at dcbprj.h include file (COMMON /CBCPRJ/, /CBIPRJ/, /CBRPRJ/, /CBLPRJ/)
!-----------------------------------------
!*** Vector projection ***
   module dcbprj_mod
   implicit none
   integer, parameter :: MAXREF = 40
   real*8 :: PROTHR
   integer :: NREFS, IPRPRJ, NPROJNUC(MAXREF), IREFADJ(2)
   logical :: LWGPOP, OWNBAS, PRJADJ, POLREF, POLNRM, MAYER, PATOMS, C1READ
   character(len=6) :: REFFIL(MAXREF)
   character(len=72) :: VECREF(2,MAXREF), VECPRJ(2)
   end module dcbprj_mod
