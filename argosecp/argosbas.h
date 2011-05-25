c     *** basis dirac ***
      integer   ncol,
     &          iqm(mxatom),
     &          jco(mxatom,mxiqm),
     &          nucij(mxatom,mxiqm,mxjco), 
     &          nrcij(mxatom,mxiqm,mxjco)
      real      expn(mxatom,mxiqm,mxjco,mxnucij),
     &          cont(mxatom,mxiqm,mxjco,mxnucij,mxnrcij)
      character frmt*1, bset*5, bskeyword*10, basfil(mxatom)*20,
     &          fmt1*8, fmt2*12, fmt3*7, 
     &          bfree(mxatom,mxiqm,mxjco,mxnucij)*80
      logical   basexpl, baslib, free(mxatom,mxiqm,mxjco)
c
c     *** basis argos ***
      integer   lmnp1(mxatom) 
c
c     *** other variables ***
      integer   bscode
c
      character expnc*90, contc(mxnucij,mxnrcij)*1,
     &          sumcontc*90, btemp*90
