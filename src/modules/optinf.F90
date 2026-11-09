!-----------------------------------------
! Module with variables previously defined
! at optinf.h include file (COMMON /OPTINF/, /OPTINF_C/)
!-----------------------------------------
! Information for geometry optimization
! in abaopt.F, abaop2.F, and abarint.F
   module optinf_mod
   use mxcent_mod, only: MXCENT
   implicit none
   integer, parameter :: MAXPRE = 10
   real*8 :: TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN, RTRJMX, ENERGY, ERGOLD, &
             ERGPRD, ERGPRO, STPNRM, STPNRO, GRADNM, THRERG, GRDTHR, THRSTP, THRSYM, &
             EVLINI, DISPLA, STPDIA(8*MXCENT), STPSYM(8*MXCENT), STPINT(8*MXCENT), &
             GRDDIA(8*MXCENT), EVAL(8*MXCENT), EVALOL(8*MXCENT), GRDINT(8*MXCENT), &
             CRDIN1(8*MXCENT), CRDINT(8*MXCENT), CNDHES(0:7)
   integer :: INDHES(0:7), INTCRD(8*MXCENT,6), ICONF(0:5), ICNSTR(8*MXCENT), IAUXRD, &
              ITOTRJ, KEPTIT, NSPMOD, NCNSTP, INDTOT, ITRNMR, ITRMAX, MAXREJ, IPRINT, &
              NCRTOT, NCART, NPROJ, NTMAT, IINTCR, IREDIC, ICRTCR, ICONDI, ITRBRK, NUMPRE, &
              IPRE, NENECALCS
   logical :: GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM, DOSPE, DOPRE, FINPRE, VRML, VRBOND, &
              VREIGV, VRCORD, VRVIBA, VRML_SYM, VISUAL, INITHS, HSFILE, BFGSR1, STEEPD, &
              RANKON, PSB, DFP, BFGS, SCHLEG, NEWTON, QUADSD, KEEPHE, BAKER, REDINT, &
              CARTCO, INRDHS, FSTORD, SNDORD, REJINI, GRDINI, MULTI, CHGRDT, CONOPT, &
              MODHES, INMDHS, FINDRE, TRSTRG, RATFUN, GDIIS, DELINT, RSTARR, LNSRCH, &
              SADDLE, REBILD, BOFILL, CMBMOD, HFPROP, CONFRM, NOAUX, NODIHE, LINDHD
   character(len=80) :: PREBTX(MAXPRE), SPBSTX
   end module optinf_mod
