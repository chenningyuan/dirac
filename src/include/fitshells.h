!     Data structure similar to that of shells but now for the
!     set of fit functions.
      LOGICAL SHARE_FIT(MXSHEL), SEGM_FIT(MXSHEL), SPHR_FIT(MXSHEL)
      REAL*8 CENT_FIT(MXSHEL, 3, 2)
      INTEGER NHKT_FIT(MXSHEL), KHKT_FIT(MXSHEL), KCKT_FIT(MXSHEL),     &
     &        ISTBAO_FIT(MXSHEL), NUCO_FIT(MXSHEL), JSTRT_FIT(MXSHEL),  &
     &        NSTRT_FIT(MXSHEL), NCENT_FIT(MXSHEL), NRCO_FIT(MXSHEL),   &
     &        NUMCF_FIT(MXSHEL), NBCH_FIT(MXSHEL), KSTRT_FIT(MXSHEL),   &
     &        LCLASS_FIT(MXSHEL), IPTSHL_FIT(MXSHEL),                   &
     &        NUMCFT_FIT(MXSHEL), KMAX_FIT, NLRGSH_FIT, NSMLSH_FIT
      COMMON /FITSHELLS/ CENT_FIT, NHKT_FIT, KHKT_FIT, KCKT_FIT,        &
     &                   ISTBAO_FIT, NUCO_FIT, JSTRT_FIT, NSTRT_FIT,    &
     &                   NCENT_FIT, SHARE_FIT, NRCO_FIT, NUMCF_FIT,     &
     &                   NBCH_FIT, KSTRT_FIT, SEGM_FIT, LCLASS_FIT,     &
     &                   IPTSHL_FIT, NUMCFT_FIT, SPHR_FIT, KMAX_FIT,    &
     &                   NLRGSH_FIT, NSMLSH_FIT
