C     choose reasonably large MAXBLLEN- so that loop unrolling gives
C     speedup but small compared with a cache size max block length
C     MXBLLEN.
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
!TROND: HFXFAC flyttet  til DCBPAM !!
      INTEGER ANGINT, ANGMIN
      CHARACTER*6 DFTTYP
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, NOPRUN,   &
     &        ATSIZE, ALDAON
      INTEGER IPRDFT
      COMMON /DFTCOM/ DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC,   &
     &                IPRDFT, NOPRUN, ATSIZE, ALDAON
      REAL*8 RADINT, DFTELS, DFTHR0, DFTHRL
      INTEGER LEBMIN, INTCHK
      COMMON /DFTRRR/ RADINT, ANGINT, ANGMIN, LEBMIN, INTCHK, DFTELS,   &
     &                DFTHR0, DFTHRL
      COMMON /DFTCHR/ DFTTYP
