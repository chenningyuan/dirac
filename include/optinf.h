C
C$Id$
C
      PARAMETER (MAXPRE = 10)
      CHARACTER*80    SPBSTX, PREBTX
      LOGICAL GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM, DOSPE,  DOPRE,
     &                FINPRE, VRML,   VRBOND, VREIGV, VRCORD, VRVIBA,
     &                VRSYMM, VISUAL, BOFILL, INITHS, HSFILE,
     &                STEEPD, RANKON, PSB,    DFP,    BFGS,   SCHLEG, 
     &                NEWTON, QUADSD, KEEPHE, BAKER,  REDINT, CARTCO,
     &                INRDHS, FSTORD, SNDORD, REJINI, GRDINI, MULTI,
     &                CHGRDT, CONOPT, MODHES, INMDHS, FINDRE, TRSTRG,
     &                RATFUN, GDIIS,  DELINT, RSTARR, LNSRCH, SADDLE,
     &                REBILD, CMBMOD, HFPROP, CONFRM, NOAUX
      COMMON /OPTINF/ TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN,
     &                RTRJMX, ENERGY, ERGOLD, ERGPRD, ERGPRO, STPNRM,
     &                STPNRO, GRADNM, THRERG, GRDTHR, THRSTP, THRSYM,
     &                EVLINI, DISPLA, STPDIA(MXRINT), 
     &                STPSYM(MXRINT), STPINT(MXRINT), 
     &                GRDDIA(MXRINT), EVAL(MXRINT), 
     &                EVALOL(MXRINT), GRDINT(MXRINT),
     &                CRDIN1(MXRINT), CRDINT(MXRINT),
     &                CNDHES(0:7), INDHES(0:7), INTCRD(MXRINT,6),
     &                ICONF(0:5), ICNSTR(MXRINT),
     &                IAUXRD, ITOTRJ, KEPTIT, NSPMOD, NCNSTP,
     &                INDTOT, ITRNMR, ITRMAX, MAXREJ, IPRINT, NCRTOT,
     &                NCART,  NPROJ,  NTMAT,  IINTCR, IREDIC, ICRTCR,
     &                ICONDI, ITRBRK, NUMPRE, IPRE,   PREBTX(MAXPRE),
     &                SPBSTX, GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM,
     &                DOSPE,  DOPRE,  FINPRE, VRML,   VRBOND, VREIGV,
     &                VRCORD, VRVIBA, VRSYMM, VISUAL, INITHS,
     &                HSFILE, STEEPD, RANKON, PSB   , DFP,    BFGS,
     &                SCHLEG, NEWTON, QUADSD, KEEPHE, BAKER,  REDINT,
     &                CARTCO, INRDHS, FSTORD, SNDORD, REJINI, GRDINI,
     &                MULTI,  CHGRDT, CONOPT, MODHES, INMDHS, FINDRE,
     &                TRSTRG, RATFUN, GDIIS,  DELINT, RSTARR, LNSRCH,
     &                SADDLE, REBILD, BOFILL, CMBMOD, HFPROP, CONFRM,
     &                NOAUX
