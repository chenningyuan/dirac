#if defined (VAR_ABASMALL)
      PARAMETER (MXCENT = 15 )
#else
      PARAMETER (MXCENT = 60 )
#endif
C
C     The number 8*mxcent is somewhat arbitrary. 2*mxcoor
C     (6*mxcent) seems to give too few redundant internal
C     coordinates for large molecules. It is therefore
C     increased slightly.
C
      PARAMETER (MXCOOR = 3 * MXCENT, MXRINT = 8 * MXCENT )
