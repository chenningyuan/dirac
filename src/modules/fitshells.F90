!-----------------------------------------
! Module with variables previously defined
! at fitshells.h include file (COMMON /FITSHELLS/)
!-----------------------------------------
!     Data structure similar to that of shells but now for the
!     set of fit functions.
   module fitshells_mod
   use maxorb_mod, only: MXCORB_FIT, MXSHEL_FIT
   implicit none
   real*8 :: CENT_FIT(MXSHEL_FIT,3,2)
   integer :: NHKT_FIT(MXSHEL_FIT), KHKT_FIT(MXSHEL_FIT), KCKT_FIT(MXSHEL_FIT), &
              ISTBAO_FIT(MXSHEL_FIT), NUCO_FIT(MXSHEL_FIT), JSTRT_FIT(MXSHEL_FIT), &
              NSTRT_FIT(MXSHEL_FIT), NCENT_FIT(MXSHEL_FIT), NRCO_FIT(MXSHEL_FIT), &
              NUMCF_FIT(MXSHEL_FIT), NBCH_FIT(MXSHEL_FIT), KSTRT_FIT(MXCORB_FIT), &
              LCLASS_FIT(MXSHEL_FIT), IPTSHL_FIT(MXSHEL_FIT), NUMCFT_FIT(MXSHEL_FIT), &
              KMAX_FIT, NLRGSH_FIT, NSMLSH_FIT
   logical :: SEGM_FIT(MXSHEL_FIT), SPHR_FIT(MXSHEL_FIT)
   end module fitshells_mod
