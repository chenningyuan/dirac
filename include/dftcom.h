!     choose reasonably large MAXBLLEN- so that loop unrolling gives
!     speedup but small compared with a cache size max block length
!     MXBLLEN.
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
!
      CHARACTER*6 DFTTYP
      CHARACTER*30 GRIDFILE, ZIPGRID_TEXT
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, NOPRUN,   &
     &        ATSIZE, ALDAON, import_grid, export_grid, ZIPGRID
      INTEGER IPRDFT, LEBMIN, INTCHK, ANGINT, ANGMIN
      REAL*8 RADINT, DFTELS, DFTHR0, DFTHRL, DFTHRI
! 
!  aspg, 2006-04-23
!  additional elements for DFTCOM common block in eri syncronization.
!  later these new elements have to be checked as to whether or not they
!  belong to DFTCOM or DFTRRR...
! 
      LOGICAL DFTHES, DFTHRS, OVWN3, DFTEST, DOVWNI, DFTIMG
!
      COMMON /DFTCOM/ DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC,   &
     &                INTCHK, ANGINT, ANGMIN, LEBMIN,                   &
     &                IPRDFT, NOPRUN, ATSIZE, ALDAON,                   &
     &                import_grid, export_grid,                         &
     &                ZIPGRID, DFTHES, DFTHRS, DFTEST, DOVWNI, DFTIMG
!
      COMMON /DFTRRR/ RADINT, DFTELS, DFTHR0, DFTHRL, DFTHRI
      COMMON /DFTCHR/ DFTTYP, GRIDFILE, ZIPGRID_TEXT
