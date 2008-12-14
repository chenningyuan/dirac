!     choose reasonably large max block length mxbllen
!     so that loop unrolling gives
!     speedup but small compared with a cache size mxbllen

      integer mxbllen
      parameter (mxbllen=160)


      character*30    gridfile, zipgrid_text

      logical         dftadd, grdone, dftrun, dftpot, dftord, dftasc,   &
     &                noprun, atsize,                                   &
     &                import_grid, export_grid, zipgrid,                &
     &                alda_hs, alda_ha, xalda, dft_thrs_changed,        &
     &                force_4c_grid, old_dft

      integer         iprdft, intchk, angint, angmin

      real*8          radint, dftels, dfthr0, dfthrl, dfthri


      common /dftcmc/ gridfile, zipgrid_text

      common /dftcml/ dftadd, grdone, dftrun, dftpot, dftord, dftasc,   &
     &                noprun, atsize,                                   &
     &                import_grid, export_grid, zipgrid,                &
     &                alda_hs, alda_ha, xalda, dft_thrs_changed,        &
     &                force_4c_grid, old_dft

      common /dftcmi/ iprdft, intchk, angint, angmin

      common /dftcmr/ radint, dftels, dfthr0, dfthrl, dfthri
