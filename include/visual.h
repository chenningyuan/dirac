C     common block {1,2,3}-dimensional plot
C     =====================================
C
C     LBIT(LBIT_PLOT_GENERAL,1) 1D
C     LBIT(LBIT_PLOT_GENERAL,2) 2D
C     LBIT(LBIT_PLOT_GENERAL,3) 3D
C     LBIT(LBIT_PLOT_GENERAL,4) read solution vector
C     LBIT(LBIT_PLOT_GENERAL,5) skip +- (SKIPEP)
C     LBIT(LBIT_PLOT_GENERAL,6) skip ++ (SKIPEE)
C     LBIT(LBIT_PLOT_GENERAL,7) use orbital string
C
C     LBIT(LBIT_PLOT_QUANTITY,1) density
C     LBIT(LBIT_PLOT_QUANTITY,2) spin density
C     LBIT(LBIT_PLOT_QUANTITY,3) current density
C     LBIT(LBIT_PLOT_QUANTITY,4) ELF
C     LBIT(LBIT_PLOT_QUANTITY,5) r \times j
C
C     LBIT(LBIT_PLOT_2D_OPTIONS,1) xy cut
C     LBIT(LBIT_PLOT_2D_OPTIONS,2) xy cut, "small"
C     LBIT(LBIT_PLOT_2D_OPTIONS,3) yz cut
C     LBIT(LBIT_PLOT_2D_OPTIONS,4) yz cut, "small"
C     LBIT(LBIT_PLOT_2D_OPTIONS,5) format 200
C
      CHARACTER*6     XVECFILE
      CHARACTER*8     ORBSTRING1, ORBSTRING2
      LOGICAL         DOPLOT
C
      COMMON /CPLOTC/ XVECFILE, ORBSTRING1, ORBSTRING2
      COMMON /CPLOTL/ DOPLOT
      COMMON /CPLOTI/ NCUBE(0:3), NXVEC, IREP_PLOT,                     &
     &                LBIT_PLOT_GENERAL,                                &
     &                LBIT_PLOT_QUANTITY,                               &
     &                LBIT_PLOT_2D_OPTIONS
      COMMON /CPLOTR/ STEP_1D, START_3D, DNORM_3D
