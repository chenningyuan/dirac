!-----------------------------------------
! Module with variables previously defined
! at dcbbas.h include file (COMMON /DCBBAS/)
!-----------------------------------------
!*** Dimensions for AO/SO (non-orthomormal) basis ***
   module dcbbas_mod
   implicit none
   integer :: NBBAS(0:7,0:2), NFBAS(2,0:2), NTBAS(0:2), NNBBAST, N2BBAST, NNBBASX, N2BBASX, &
              N2BBASXQ, MXBBAS, MXFBAS, N2BAS(2), N2BASQ(2), N2BAST, N2BASTQ, &
              IBBAS(0:7,0:2), I2BBASX(0:7,0:7,2,2), IBAS(2), I2BASX(2,2), I2BAST(2), NZT, &
              I2TMT(2), N2TMT
   end module dcbbas_mod
