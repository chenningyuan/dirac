C     visualization/integration common block
C     ======================================
C
C     LBIT(LBIT_VISUAL, 1) list of points
C     LBIT(LBIT_VISUAL, 2) 2D
C     LBIT(LBIT_VISUAL, 3) 3D
C     LBIT(LBIT_VISUAL, 4)
C     LBIT(LBIT_VISUAL, 5)
C     LBIT(LBIT_VISUAL, 6)
C     LBIT(LBIT_VISUAL, 7) use orbital string
C     LBIT(LBIT_VISUAL, 8) integrate
C     LBIT(LBIT_VISUAL, 9)
C     LBIT(LBIT_VISUAL,10)
C
      INTEGER MAX_NR_XVECTORS
      PARAMETER (MAX_NR_XVECTORS = 25)
      PARAMETER (MAX_NR_POINTS_IN_LIST = 100)
C
      CHARACTER       XVECTOR_FILE(MAX_NR_XVECTORS)*6,                  &
     &                ORBITAL_STRING(2)*8,                              &
     &                ORIENTATION(3)*2
      LOGICAL         DO_VISUAL
C
      COMMON /CBVISC/ XVECTOR_FILE,                                     &
     &                ORBITAL_STRING,                                   &
     &                ORIENTATION
      COMMON /CBVISL/ DO_VISUAL
      COMMON /CBVISI/ LBIT_VISUAL,                                      &
     &                NCUBE(0:3),                                       &
     &                NR_POINTS_IN_LIST,                                &
     &                NR_XVECTORS,                                      &
     &                IQ_XVECTOR(MAX_NR_XVECTORS),                      &
     &                IREP_XVECTOR(MAX_NR_XVECTORS),                    &
     &                IWHICH_XVECTOR(MAX_NR_XVECTORS),                  &
     &                IOPTION_XVECTOR(MAX_NR_XVECTORS),IPRVIS
      COMMON /CBVISR/ ADD_3D, COMMON_FACTOR,                            &
     &                XYZ_LIST(MAX_NR_POINTS_IN_LIST,3),                &
     &                ORIGIN_2D(3), GAUGE_ORIGIN(3),                    &
     &                ADD_2D, STEP_2D
