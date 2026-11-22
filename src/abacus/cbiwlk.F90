!-----------------------------------------
! Module with variables previously defined
! at cbiwlk.h include file (COMMON /CBIWLK/)
!-----------------------------------------
   module cbiwlk_mod
   use mxcent_mod, only: MXCENT, MXCOOR
   implicit none
   integer, parameter :: MAXTMP = 20
   real*8 :: TOLST, TRUSTR, TRUSTI, TRUSTD, REJMIN, REJMAX, RTRMIN, RTRGOD, XMXNUC, ZERGRD, &
             DISPLC, STRMOM(3*MXCENT), ECKGEO(3,MXCOOR), SCALCO(3,MXCENT), THRLDP, ANHFAC, &
             TRUMAX, TEMP(MAXTMP)
   integer :: ISTMOM(3*MXCENT), NSTMOM, IPRWLK, IWKTYP, IWKIND, IMODE, ISCTYP, &
              IPART(MXCENT), NZEROG, IZEROG(MXCOOR), NTEMP
   logical :: DOREPW(0:7), START, WFPRED, REJECT, KEEPSY, VIBCNV, IMAGE, STRICT, NOORTH, &
              NATCON, V3CAL, VIBAVE, NMODIF, ECKART, DOTEMP, DOCENT, REUSED
   end module cbiwlk_mod
