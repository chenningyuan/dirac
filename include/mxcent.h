#if defined (VAR_TESTIBM)
      PARAMETER (MXCENT =  8, MXCOOR = 3*MXCENT)
#else
#if defined (VAR_ABASMALL)
      PARAMETER (MXCENT = 8, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXCENT = 14, MXCOOR = 3*MXCENT)
#endif
#endif
