C
C ERAV, ERRAB, ERRMX    : Average, absolute average and maximum
C 		          errors in numerical SMAT.
C (MAXI, MAXJ)          : Location of max error.
C MXLB1, MXLB2          : Labels for orbitals above.
C
C RERRAV, RERRAB RERRMX : Relative values
C (MAXIR,MAXJR)          : Relatove values.
C RMXLB1, RMXLB2        : Relative values.

      CHARACTER*12 MXLB1,MXLB2,RMXLB1,RMXLB2

      COMMON /GRTST/ ERRAV, ERRAB, ERRMX, RERRAV, RERRAB, RERRMX,
     &               SMATDF, SMATAN, SMATFR,
     &               MAXI, MAXJ, MAXIR, MAXJR, NLSOS, NSSOS, INTCHK
