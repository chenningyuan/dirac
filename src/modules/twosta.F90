!-----------------------------------------
!
! Module with variables previously defined
!
! at twosta.h include file
!
!-----------------------------------------

module twosta_mod

   use maxaqn_mod, only: MXQN
   implicit none

   logical :: TKTIME

   real*8 :: TR000, THERI, TODCVE, TEXCOE, TC10IN, TC1EIN, TC2HIN, &
             TC2EIN, TSYM2S, TSYMOU, TDSOUT, TSPOSY, TDRSYM, TINTEX, &
             TPATH1, TPATH2, TFCKOU, T2INT

   real*8 :: THERIX(0:4*MXQN-2), TR000X(0:4*MXQN-2)

   integer :: LWTOT, MWTOT, MWFCAB, MWFCCD, MWPSO, MWSOIN, MWAOIN, &
              MWHRIN, MWHRND, MWHCIN, MWRJ00, MWHRSQ, MWSCR1, MWC1DR, &
              MWC2DR, MWC2HI, MWC2EI, MWINTE, MWPPRI, MWDRSY, N2WRIT

   integer :: MXHERMEM(0:3)

end module twosta_mod