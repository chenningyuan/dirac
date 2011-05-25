c *** SO index ***
      integer   csoindex, soindex, mcrs(mxatom), nmcrs
      logical   ecpchk
c
c *** ARPP / SOPP ***
      integer   lcru(mxatom),            llsu(mxatom),
     &          nbfcr(mxatom,mxiqm),     nbfls(mxatom,mxiqm),
     &          ncr(mxatom,mxiqm,mxnbf), nls(mxatom,mxiqm,mxnbf)
      real*8    zcr(mxatom,mxiqm,mxnbf), ccr(mxatom,mxiqm,mxnbf),
     &          zls(mxatom,mxiqm,mxnbf), cls(mxatom,mxiqm,mxnbf)
