C
C$Id$
C
      PARAMETER (MAXPRE = 10)
      CHARACTER*80    SPBSTX, PREBTX
      LOGICAL GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM, DOSPE,  DOPRE,
     &                FINPRE, VRML,   VRBOND, VREIGV, VISUAL, INTRAC,
     &                INITHS, HSFILE, STEEPD, RANKON, PSB   , DFP,
     &                BFGS,   SCHLEG, NEWTON, QUADSD, KEEPHE, REINIH,
     &                BAKER,  REDINT, CARTCO, INRDHS, FSTORD, SNDORD,
     &                REJINI, GRDINI, MBFGS, CHGRDT,  GRDTCK, TWOGRD(3)
      COMMON /OPTINF/ TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN,
     &                RTRJMX, ENERGY, ERGOLD, ERGPRD, ERGPRO, STPNRM,
     &                STPNRO, GRADNM, THRERG, GRDTHR, THRSTP, THRSYM,
     &                EVLINI, DISPLA, STPDIA(2*MXCOOR), 
     &                STPSYM(2*MXCOOR), STPINT(2*MXCOOR), 
     &                GRDDIA(2*MXCOOR), EVAL(2*MXCOOR), 
     &                EVALOL(2*MXCOOR), GRDINT(2*MXCOOR),
     &                CNDHES(0:7),INDHES(0:7),INTCRD(2*MXCOOR,6),
     &                INDTOT, ITRNMR, ITRMAX, MAXREJ, IPRINT, NCRTOT,
     &                NCART,  NPROJ,  NTMAT,  IINTCR, ICRTCR,
     &                ICONDI, ITRBRK, NUMPRE, IPRE,   PREBTX(MAXPRE),
     &                SPBSTX, GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM,
     &                DOSPE,  DOPRE,  FINPRE, VRML,   VRBOND, VREIGV,
     &                VISUAL, INTRAC, INITHS, HSFILE, STEEPD, RANKON,
     &                PSB   , DFP,    BFGS,   SCHLEG, NEWTON, QUADSD,
     &                KEEPHE, REINIH, BAKER,  REDINT, CARTCO, INRDHS,
     &                FSTORD, SNDORD, REJINI, GRDINI, MBFGS, CHGRDT,
     &                GRDTCK, TWOGRD(3)
