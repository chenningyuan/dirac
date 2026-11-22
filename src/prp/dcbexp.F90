!-----------------------------------------
! Module with variables previously defined
! at dcbexp.h include file (COMMON /DCBEXP/)
!-----------------------------------------
!C
!C
!C*****************************************************************************
!C  DCBEXP:  Expectation values
!C
   module dcbexp_mod
   use mxcent_mod, only: MXCENT
   implicit none
   integer, parameter :: MAXEXP = 350
   integer :: IPREXP, NEXPP, LEXPP(MAXEXP), IPDIPOLE(3), IPQUADRU(6), IPEFG(6*MXCENT), &
              IPPVC(MXCENT), IPEFT(15*MXCENT), IPEFTNTL(15*MXCENT), IPDSUSC(39), &
              IPEFFDEN(MXCENT), IPEFN(3*MXCENT)
   logical :: ORBANA, PRJANA, PRPCAN
   end module dcbexp_mod
