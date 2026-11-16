!-----------------------------------------
! Module with variables previously defined
! at dcbkrmc.h include file (COMMON /DCLKRMC/, /DCOKRMC/, /DCIKRMC/, /DCRKRMC/, /DCCKRMC/)
!-----------------------------------------
!*** Information for KR-MCSCF calculations ***
   module dcbkrmc_mod
   implicit none
! MXGAS for the dimensions; mxgas.h stays a header while dcborb.h includes it
#include "../include/mxgas.h"
   private :: MXGAS, MXPAS
   integer, parameter :: NMCFLAG = 20
   real*8 :: DKRMC_CNVINT(2), DKRMC_THRGRD, DKRMC_THRPCI, DKRMC_MVOFAC
   integer :: NKRMCISH(2), NKRMCASH(2), NKRMCSSH(2), NKRMCGSH(2,MXGAS), NKRMCGSP(2,MXGAS), &
              NKRMCESH(2), NKRMCPSH(2), NKRMCFRO(2), NKRMCAELEC, NKRMCGAS, NKRMC_MK2REF, &
              NKRMC_MK2DEL, NKRMC_MINMK2, NKRMC_MAXMK2, IKRMC_ITRINT(2), IKRMC_INTDEF, &
              IKRMC_INTFLG, IKRMC_INTBUF, MKRMC_MXMICRO, MKRMC_MXMACRO, MKRMC_MAXBCK, &
              IKRMC_STATE, IKRMC_SYMMETRY, IPRKRMC, IKRMC_MEMFAC, IKRMC_SVRONO(2)
   logical :: MCFLAG(NMCFLAG), KRMC_SKIPEE, KRMC_SKIPEP, KRMC_NOCI, KRMC_UCIBOS, &
              KRMC_CHCKJZ, KRMC_NO1PDENS, KRMC_SAVE_REORDERED_NOS, KRMC_FULL_CI
   character(len=72) :: KRMC_FRZSTR(2), KRMC_DELSTR(2)
   character(len=8) :: KRMC_CIPROGRAM
   end module dcbkrmc_mod
