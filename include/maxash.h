C     MAXASH = maximum number of active orbitals
#if !defined (VAR_GUGACI)
      PARAMETER ( MAXASH = 100 )
#else
#if !defined (VAR_BIGGUGA)
      PARAMETER ( MAXASH = 15 )
C        ** 15 is max. active from Per's CAS GUGA 841016
#else
      PARAMETER ( MAXASH = 30 )
C        ** 30 is max. active from Per's CAS GUGA 841016
C           modified for 64 bit integers Nov. 1986
#endif
#endif
