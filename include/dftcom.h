C     choose reasonably large MAXBLLEN- so that loop unrolling gives
C     speedup but small compared with a cache size max block length
C     MXBLLEN.
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
C
      CHARACTER*6 DFTTYP
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, NOPRUN,   &
     &        ATSIZE, ALDAON
      INTEGER IPRDFT, LEBMIN, INTCHK, ANGINT, ANGMIN
      REAL*8 RADINT, DFTELS, DFTHR0, DFTHRL, DFTHRI
C     
C aspg, 2006-04-23 
C additional elements for DFTCOM common block in eri syncronization.
C later these new elements have to be checked as to whether or not they
C belong to DFTCOM or DFTRRR... 
C     
C old comment by trond:
CTROND: HFXFAC flyttet  til DCBPAM !!
C so for now HFXFAC stays in DCBPAM, later this has to be checked
C     
C     REAL*8  HFXFAC, HFXATT, HFXMU
      REAL*8  HFXATT, HFXMU
      LOGICAL DFTHES, DFTHRS, OVWN3, DFTEST, DOVWNI, DFTIMG
C     
      COMMON /DFTCOM/ DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC,   &
     &                INTCHK, ANGINT, ANGMIN, LEBMIN,                   &
     &                IPRDFT, NOPRUN, ATSIZE, ALDAON,                   &
C    &                HFXFAC, HFXATT, HFXMU,                            &
     &                HFXATT, HFXMU,                                    &
     &                DFTHES, DFTHRS, DOVWN3, DFTEST, DOVWNI, DFTIMG
C     
      COMMON /DFTRRR/ RADINT, DFTELS, DFTHR0, DFTHRL, DFTHRI
      COMMON /DFTCHR/ DFTTYP
