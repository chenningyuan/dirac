#if defined (VAR_TESTIBM)
      PARAMETER (MXCENT =  8, MXCOOR = 3*MXCENT)
#else
#if defined (VAR_ABASMALL)
      PARAMETER (MXCENT = 8, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXCENT = 60, MXCOOR = 3*MXCENT)
#endif
#endif
