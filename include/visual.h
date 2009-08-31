!     visualization/integration common block
!     ======================================
!
!     LBIT(LBIT_VISUAL, 1) list of points
!     LBIT(LBIT_VISUAL, 2) 2D
!     LBIT(LBIT_VISUAL, 3) 3D
!     LBIT(LBIT_VISUAL, 4)
!     LBIT(LBIT_VISUAL, 5)
!     LBIT(LBIT_VISUAL, 6)
!     LBIT(LBIT_VISUAL, 7) use orbital string
!     LBIT(LBIT_VISUAL, 8) integrate
!     LBIT(LBIT_VISUAL, 9)
!     LBIT(LBIT_VISUAL,10)
!
      integer max_nr_xvectors, max_nr_points_in_list
      PARAMETER (MAX_NR_XVECTORS = 25)
      PARAMETER (MAX_NR_POINTS_IN_LIST = 100)
!
      CHARACTER       XVECTOR_FILE(MAX_NR_XVECTORS)*6,                  &
     &                ORBITAL_STRING(2)*72,                             &
     &                ORIENTATION(3)*2
      logical         do_visual
      integer         lbit_visual,                                      &
     &                ncube,                                            &
     &                nr_points_in_list,                                &
     &                nr_xvectors,                                      &
     &                iq_xvector,                                       &
     &                irep_xvector,                                     &
     &                iwhich_xvector,                                   &
     &                ioption_xvector, iprvis
      real(8)         add_3d, common_factor,                            &
     &                xyz_list,                                         &
     &                origin_2d, gauge_origin,                          &
     &                add_2d, step_2d

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
