C     common block {1,2,3}-dimensional plot
C     =====================================
C
C     LBIT(LBIT_PLOT_GENERAL,1)    1D
C     LBIT(LBIT_PLOT_GENERAL,2)    2D
C     LBIT(LBIT_PLOT_GENERAL,3)    3D
C     LBIT(LBIT_PLOT_GENERAL,4)    read solution vector
C     LBIT(LBIT_PLOT_GENERAL,5)    skip +- (SKIPEP)
C     LBIT(LBIT_PLOT_GENERAL,6)    skip ++ (SKIPEE)
C     LBIT(LBIT_PLOT_GENERAL,7)    use orbital string
C     LBIT(LBIT_PLOT_GENERAL,8)    integrate
C
C     LBIT(LBIT_PLOT_QUANTITY,1)   density
C     LBIT(LBIT_PLOT_QUANTITY,2)   magnetization
C     LBIT(LBIT_PLOT_QUANTITY,3)   current density
C     LBIT(LBIT_PLOT_QUANTITY,4)   ELF
C     LBIT(LBIT_PLOT_QUANTITY,5)   r \times j
C     LBIT(LBIT_PLOT_QUANTITY,6)   divergence(magnetization)
C     LBIT(LBIT_PLOT_QUANTITY,7)   divergence(current density)
C     LBIT(LBIT_PLOT_QUANTITY,8)   curl(magnetization)
C     LBIT(LBIT_PLOT_QUANTITY,9)   curl(current density)
C
C     LBIT(LBIT_PLOT_QUANTITY,10) .OPT0 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,11) .OPT1 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,12) .OPT2 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,13) .OPT3 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,14) .OPT4 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,15) .OPT5 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,16) .OPT6 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,17) .OPT7 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,18) .OPT8 (debug and testing of new functionalities)
C     LBIT(LBIT_PLOT_QUANTITY,19) .OPT9 (debug and testing of new functionalities)
C
C     LBIT(LBIT_PLOT_QUANTITY,20)  diamagnetic current density
C     LBIT(LBIT_PLOT_QUANTITY,21)  kinetic energy density
C     LBIT(LBIT_PLOT_QUANTITY,22)  helicity density
C     LBIT(LBIT_PLOT_QUANTITY,23)  \gamma_5 density
C
      CHARACTER*6     XVECFILE
      CHARACTER*8     ORBSTRING1, ORBSTRING2
      LOGICAL         DOPLOT
C
      COMMON /CPLOTC/ XVECFILE, ORBSTRING1, ORBSTRING2
      COMMON /CPLOTL/ DOPLOT
      COMMON /CPLOTI/ NCUBE(0:3), NXVEC, IREP_PLOT,                     &
     &                LBIT_PLOT_GENERAL, LBIT_PLOT_QUANTITY
      COMMON /CPLOTR/ STEP_1D, ADD_3D
