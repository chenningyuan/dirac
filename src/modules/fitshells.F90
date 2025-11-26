!-----------------------------------------
! Module with variables previously defined
! at fitshells.h include file
!-----------------------------------------
   module fitshells_mod
!
   use primitives_mod, only: MXSHEL_FIT, MXCORB_FIT
   implicit none
!
   logical :: SEGM_FIT(MXSHEL_FIT)
   logical :: SPHR_FIT(MXSHEL_FIT)
   real(8) :: CENT_FIT(MXSHEL_FIT, 3, 2)
!
   integer :: NHKT_FIT(MXSHEL_FIT), KHKT_FIT(MXSHEL_FIT)
   integer :: KCKT_FIT(MXSHEL_FIT), ISTBAO_FIT(MXSHEL_FIT)
   integer :: NUCO_FIT(MXSHEL_FIT), JSTRT_FIT(MXSHEL_FIT)
   integer :: NSTRT_FIT(MXSHEL_FIT), NCENT_FIT(MXSHEL_FIT)
   integer :: NRCO_FIT(MXSHEL_FIT), NUMCF_FIT(MXSHEL_FIT)
   integer :: NBCH_FIT(MXSHEL_FIT), KSTRT_FIT(MXCORB_FIT)
   integer :: LCLASS_FIT(MXSHEL_FIT)!, IPTSHL_FIT(MXSHEL_FIT)
   integer :: NUMCFT_FIT(MXSHEL_FIT)
   integer :: KMAX_FIT, NLRGSH_FIT, NSMLSH_FIT
!
end module fitshells_mod
