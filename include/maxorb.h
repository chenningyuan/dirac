C     MAXORB = maximum number of orbitals
#if defined (VAR_SIRBIG)
      PARAMETER ( MAXORB = 1200 )
#else
      PARAMETER ( MAXORB = 400 )
#endif
C
C     SIRIPC keeps track of the input processing in SIRIUS.
C
