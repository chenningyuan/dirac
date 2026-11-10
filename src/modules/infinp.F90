!-----------------------------------------
! Module with variables previously defined
! at src/include/infinp.h include file (DIRAC version; the DALTON one is infinp_DALTON_mod) (COMMON /INTINP/, /LOGINP/, /RELINP/, /CHRINP/)
!-----------------------------------------
!  -*- mode:fortran; fortran-continuation-string: "&" -*-
!     this common block contains general Sirius input read in sirius/sirinp.F
!     (specified under **WAVE FUNCTION).  /hjaaj Oct 2003
!     - SCF specific input is in scbrhf.h
!     - orbital specifications are in inforb.h
   module infinp_mod
   use maxorb_mod, only: MXCORB
   implicit none
   integer, parameter :: NFLAG = 80
   integer, parameter :: MAXRTS = 100
   integer, parameter :: MXFELT = 20
   real*8 :: SPIN, POTNUC, EPSOL, EPSTAT, EPPN, RSOL(3), THRGRD, THRPWF, THRCI, THRMC, &
             THRCGR, EFIELD(MXFELT), CMAXMO, THROVL, THRSSY, DEFLVL, WEIGHT_MCAVER(MAXRTS), &
             THRCIDFT
   integer :: NFIELD, ISPIN, ISTATE, LSYM, NACTEL, MCTYPE, LSOLMX, NLMSOL, NELMN1, NELMX1, &
              NELMN3, NELMX3, LROOTS, NROOTS, IROOT(MAXRTS), NOROT(MXCORB), IMOORD(MXCORB), &
              IORTO, ICI0, KDEL, ICHECK, NTIT, MAXMAC, MAXMIC, MAXJT, MAXCIT, MAXUIT, &
              MAXAPM, MAXABS, ITRLVL, ITRFIN, JCHSYM, JCHORB, NROOCI, ISTACI, MXCIMA, &
              ICICNO, IMCCNO, NMCAVER
   logical :: FLAG(NFLAG), DOSCF, DOMP2, DOCINO, DOCI, DOMC, DORSP, FCVORB, LNOROT, LMOORD, &
              DIRFCK, CORHOL, CORRLX, RESPHP, JOLSEN, ABAIPH, INERSI, INERSF, DODFT, &
              DONEVPT, HSROHF, BOYORB, PIPORB, ADDMP2, DOFCI, DOCISD, DO_CUBE, DOCISRDFT, &
              DOHFSRDFT, DOMCSRDFT, ADDSRI, SRHYBR
   character(len=4) :: CENT(MXCORB), TYPE(MXCORB)
   character(len=60) :: TITLE(6)
   character(len=72) :: TITMOL(2)
   character(len=8) :: LFIELD(MXFELT)
   end module infinp_mod
