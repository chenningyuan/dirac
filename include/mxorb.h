#if defined (VAR_TESTIBM)
      PARAMETER (MXSHEL =  15, MXPRIM =  60, MXCORB =  30,
     &           MXORBT = MXCORB*(MXCORB + 1)/2)
#else
#if defined (VAR_SIRBIG)
      PARAMETER (MXSHEL = 400, MXPRIM = 2000, MXCORB = 1200,
     &           MXORBT = MXCORB*(MXCORB + 1)/2)
#else
      PARAMETER (MXSHEL = 500, MXPRIM = 5000, MXCORB = 5000,
     &           MXORBT = MXCORB*(MXCORB + 1)/2)
#endif
#endif
