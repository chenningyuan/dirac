C hjaaj July 2000:
C MXSHEL = max number of symmetry independent contracted shells
C MXPRIM = max number of symmetry independent primitive shells
C  (e.g. a p shell contains 3 components)
C MXCORB = max total number of contracted functions = max NBASIS
#if defined (VAR_ABASMALL)
      PARAMETER (MXSHEL =  50, MXPRIM =  200, MXCORB =  500,
     &           MXORBT = MXCORB*(MXCORB + 1)/2)
#else
      PARAMETER (MXSHEL = 500, MXPRIM = 2000, MXCORB = 5000,
     &           MXORBT = MXCORB*(MXCORB + 1)/2)
#endif
