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
C aspg, 2006-09-20
C DFITXC is added to allow using Birkenheuer et al's scheme for the
C evaluation of the XC potential with the help of density fitting. More
C details in DFTDRV.F, and in IJQC 102, 743-761 (2005)
C     
C old comment by trond:
CTROND: HFXFAC flyttet  til DCBPAM !!
C so for now HFXFAC stays in DCBPAM, later this has to be checked
C     
C     REAL*8  HFXFAC, HFXATT, HFXMU
C      REAL*8  HFXATT, HFXMU
      LOGICAL DFTHES, DFTHRS, OVWN3, DFTEST, DOVWNI, DFTIMG, DFITXC
C     
      COMMON /DFTCOM/ DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC,   &
     &                INTCHK, ANGINT, ANGMIN, LEBMIN,                   &
     &                IPRDFT, NOPRUN, ATSIZE, ALDAON,                   &
C    &                HFXFAC, HFXATT, HFXMU,                            &
C    &                HFXATT, HFXMU,                                    &
     &                DFTHES, DFTHRS, DOVWN3, DFTEST, DOVWNI, DFTIMG,   &
     &                DFITXC
C     
      COMMON /DFTRRR/ RADINT, DFTELS, DFTHR0, DFTHRL, DFTHRI
      COMMON /DFTCHR/ DFTTYP
