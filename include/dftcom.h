C     choose reasonably large MAXBLLEN- so that loop unrolling gives
C     speedup but small compared with a cache size max block length
C     MXBLLEN.
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
!TROND: HFXFAC flyttet  til DCBPAM !!
      CHARACTER*6 DFTTYP
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, NOPRUN,   &
     &        ATSIZE, ALDAON, DFTGRAC, GRACON
      INTEGER IPRDFT, LEBMIN, INTCHK, ANGINT, ANGMIN
      REAL*8 RADINT, DFTELS, DFTHR0, DFTHRL, DFTHRI,                    &
     &       GRACIP, GRACALPHA, GRACBETA, EIGMAX
      COMMON /DFTCOM/ DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC,   &
     &                ANGINT, ANGMIN, LEBMIN, 
     &                IPRDFT, NOPRUN, ATSIZE, ALDAON, DFTGRAC, GRACON
      COMMON /DFTRRR/ RADINT, INTCHK, DFTELS, DFTHR0, DFTHRL, DFTHRI,   &
     &                GRACIP, GRACALPHA, GRACBETA, EIGMAX
      COMMON /DFTCHR/ DFTTYP
