      INTEGER NCSP
      PARAMETER (NCSP =                                                 &
     &           ((MXQN+2)*(MXQN+1)*(3*MXQN**2+6*MXQN+1)*MXQN)/60)

      REAL*8 CSP(NCSP)
      INTEGER ISPADR(MXQN)
      COMMON /SPHTRM/ CSP, ISPADR
