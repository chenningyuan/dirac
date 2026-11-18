!-----------------------------------------
! Module with variables previously defined
! at dcbkrci.h include file (COMMON /DCLKRCI/, /DCIKRCI/, /DCRKRCI/, /DCCKRCI/)
!-----------------------------------------
!     *** information for KR-CI calculations ***
   module dcbkrci_mod
   implicit none
! MXGAS for the dimensions; mxgas.h stays a header while dcborb.h includes it
#include "../include/mxgas.h"
   private :: MXGAS, MXPAS
   integer, parameter :: MAX_NKRCI_MAX_SYM = 128
   real*8 :: DKRCI_CNVINT(2), DKRCI_THRPCI, DKRCI_THRGCI, DKRCI_THRECI
   integer :: NKRCIISH(2), NKRCIASH(2), NKRCISSH(2), NKRCIGSH(2,MXGAS), NKRCIGSP(2,MXGAS), &
              NKRCIESH(2), NKRCIPSH(2), NKRCIFRO(2), NKRCIAELEC, NKRCIGAS, NKRCI_MK2REF, &
              NKRCI_MK2DEL, NKRCI_MINMK2, NKRCI_MAXMK2, IKRCI_ITRINT(2), IKRCI_INTDEF, &
              IKRCI_INTFLG, IKRCI_INTBUF, NKRCI_MAX_SYM, NKRCI_CIROOTS(MAX_NKRCI_MAX_SYM), &
              IPRKRCI, MAXCIT_KRCI, MXCIV_KRCI, IANACI_KRCI, KTRLVL_KRCI, IRESTRK_KRCI, &
              IDOPARIO_KRCI, IKRCI_MEMFAC, IKRCI_SVRONO(2), NKRCI_SELOM(MAX_NKRCI_MAX_SYM), &
              NOMEGASEL_STATENR
   logical :: KRCI_MKCNST, KRCI_UCIBOS, CRDFO_MAT_KRCI, CWRTFO_MAT_KRCI, CSHMEMO_KRCI, &
              CINT_REORD_KRCI, CINT_SPLIT_KRCI, NATOLCR_KRCI, CHCKPT_KRCI, CUSE_PBLFL_KRCI, &
              TRANATO_KRCI, LOWSORT_KRCI, STATE_SELECT_KRCI, DOMCANA_KRCI, &
              SAVE_REORDERED_NOS_KRCI, FCIDUMP_KRCI, XPSINTERFACE_KRCI
   character(len=8) :: KRCI_CIPROGRAM
   end module dcbkrci_mod
