!-----------------------------------------
! Module with variables previously defined
! at optinf.h include file
!-----------------------------------------
   module optinf_mod
!
   use mxcent_mod, only: MXCENT
   implicit none
   public

   !----------------------------------------------------------
   ! Parameter
   !----------------------------------------------------------
   integer, parameter :: MAXPRE = 10
!
   !----------------------------------------------------------
   ! REAL*8 (double precision)
   !----------------------------------------------------------
   real(8) :: TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN, RTRJMX
   real(8) :: ENERGY, ERGOLD, ERGPRD, ERGPRO, STPNRM, STPNRO, GRADNM
   real(8) :: THRERG, GRDTHR, THRSTP, THRSYM, EVLINI, DISPLA
!
   real(8) :: STPDIA(8*MXCENT)
   real(8) :: STPSYM(8*MXCENT)
   real(8) :: STPINT(8*MXCENT)
   real(8) :: GRDDIA(8*MXCENT)
   real(8) :: EVAL(8*MXCENT)
   real(8) :: EVALOL(8*MXCENT)
   real(8) :: GRDINT(8*MXCENT)
   real(8) :: CRDIN1(8*MXCENT)
   real(8) :: CRDINT(8*MXCENT)
!
   real(8) :: CNDHES(0:7)
!
   !----------------------------------------------------------
   ! INTEGER
   !----------------------------------------------------------
   integer :: INDHES(0:7)
   integer :: INTCRD(8*MXCENT, 6)
   integer :: ICONF(0:5)
   integer :: ICNSTR(8*MXCENT)
   integer :: IAUXRD, ITOTRJ, KEPTIT, NSPMOD, NCNSTP
   integer :: INDTOT, ITRNMR, ITRMAX, MAXREJ, IPRINT
   integer :: NCRTOT, NCART, NPROJ, NTMAT
   integer :: IINTCR, IREDIC, ICRTCR, ICONDI, ITRBRK
   integer :: NUMPRE, IPRE, NENECALCS
!
   !----------------------------------------------------------
   ! LOGICAL
   !----------------------------------------------------------
   logical :: GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM
   logical :: DOSPE, DOPRE, FINPRE, VRML, VRBOND
   logical :: VREIGV, VRCORD, VRVIBA, VRML_SYM, VISUAL
   logical :: BOFILL, INITHS, HSFILE, BFGSR1, STEEPD
   logical :: RANKON, PSB, DFP, BFGS, SCHLEG
   logical :: NEWTON, QUADSD, KEEPHE, BAKER
   logical :: REDINT, CARTCO, INRDHS, FSTORD, SNDORD
   logical :: REJINI, GRDINI, MULTI, CHGRDT
   logical :: CONOPT, MODHES, INMDHS, FINDRE
   logical :: TRSTRG, RATFUN, GDIIS, DELINT, RSTARR
   logical :: LNSRCH, SADDLE, REBILD, CMBMOD
   logical :: HFPROP, CONFRM, NOAUX, NODIHE, LINDHD
!
   !----------------------------------------------------------
   ! CHARACTER (COMMON /OPTINF_C/)
   !----------------------------------------------------------
   character(len=80) :: PREBTX(MAXPRE)
   character(len=80) :: SPBSTX
!
   end module optinf_mod
