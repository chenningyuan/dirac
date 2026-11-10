!-----------------------------------------
! Module with variables previously defined
! at src/abacus/infinp.h include file (DALTON version; the DIRAC one is infinp_mod) (COMMON /INTINP/, /LOGINP/, /RELINP/, /CHRINP/)
!-----------------------------------------
   module infinp_DALTON_mod
   use maxorb_mod, only: MXCORB
   implicit none
   integer, parameter :: NFLAG = 80
   integer, parameter :: MAXRTS = 20
   integer, parameter :: MXFELT = 20
   real*8 :: SPIN, POTNUC, EPSOL, EPSTAT, RSOL(3), THRGRD, THRPWF, THRRHF, THRCI, THRMC, &
             THRCGR, EFIELD(MXFELT), TITMOL(12,2), CMAXMO, THROVL, THRSSY
   integer :: ISPIN, ISTATE, LSYM, NACTEL, MCTYPE, LSOLMX, NLMSOL, NELMN1, NELMX1, NELMN3, &
              NELMX3, LROOTS, NROOTS, IROOT(MAXRTS), KAVER, KSYM(4), NOROT(MXCORB), &
              ITMPMO(2*MXCORB), IORTO, ICI0, KDEL, ICHECK, NTIT, MAXMAC, MAXMIC, MAXJT, &
              MAXCIT, MAXUIT, MAXAPM, MAXABS, MAXFCK, ITRLVL, ITRFIN, JCHSYM, JCHORB, &
              NROOCI, ISTACI, NFIELD, MXCIMA, ICICNO, IMCCNO, MXDIIS, MAXEVC
   logical :: FLAG(NFLAG), DORHF, DOMP2, DOCINO, DOCI, DOMC, DORSP, H1VORB, LNOROT, LMOORD, &
              DIRFCK, CORHOL, CORRLX, RESPHP, JOLSEN, ABAIPH, INERSI, INERSF, PARSIR
   character(len=10) :: LVLCDK(2)
   character(len=4) :: CENT(MXCORB), TYPE(MXCORB), CTMPMO(MXCORB)
   character(len=60) :: TITLE(6)
   character(len=8) :: LFIELD(MXFELT)
   end module infinp_DALTON_mod
