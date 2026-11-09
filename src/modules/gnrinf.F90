!-----------------------------------------
! Module with variables previously defined
! at gnrinf.h include file (COMMON /GNRINF/, /GNRCHR/)
!-----------------------------------------
   module gnrinf_mod
   implicit none
   real*8 :: GRADML, PANAS
   integer :: KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS
   logical :: TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV, HRINPC, SRINPC, RDINPC, &
              RDMLIN, PARCAL, DIRCAL, WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS, &
              DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL, TOTSYM, NMWALK
   character(len=70) :: BASDIR
   end module gnrinf_mod
