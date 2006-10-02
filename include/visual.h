C     common block {1,2,3}-dimensional plot
C     =====================================
C
      CHARACTER*6 XVECFILE
      CHARACTER*8 ORBSTRING1, ORBSTRING2
C
      LOGICAL         READXVECFILE,                                     &
     &                DOPLOT, DO1DPLOT, DO2DPLOT, DO3DPLOT,             &
     &                PLOTDENSITY, PLOTSPINDENSITY, PLOTCURRENTDENSITY, &
     &                PLOTELF, PLOTLAPLACIAN, USE_ORBSTRING,            &
     &                PLOTRXJ,                                          &
     &                PLOTSKIPEP, PLOTSKIPEE,                           &
     &                CUT_XY, CUT_XY_SMALL, CUT_YZ, CUT_YZ_SMALL,       &
     &                FORMAT_200
C
      COMMON /CPLOTC/ XVECFILE, ORBSTRING1, ORBSTRING2
C
      COMMON /CPLOTL/ READXVECFILE,                                     &
     &                DOPLOT, DO1DPLOT, DO2DPLOT, DO3DPLOT,             &
     &                PLOTDENSITY, PLOTSPINDENSITY, PLOTCURRENTDENSITY, &
     &                PLOTELF, PLOTLAPLACIAN, USE_ORBSTRING,            &
     &                PLOTRXJ,                                          &
     &                PLOTSKIPEP, PLOTSKIPEE,                           &
     &                CUT_XY, CUT_XY_SMALL, CUT_YZ, CUT_YZ_SMALL,       &
     &                FORMAT_200
C
      COMMON /CPLOTI/ IPLOT, NCUBE(0:3), NXVEC, IREP_PLOT
C
      COMMON /CPLOTR/ D1DSTEP, D3DSTART, D3DNORM
