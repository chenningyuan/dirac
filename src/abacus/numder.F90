!-----------------------------------------
! Module with variables previously defined
! at numder.h include file (COMMON /NUMDER/, /NMDPRP/)
!-----------------------------------------
   module numder_mod
   implicit none
   integer :: KDRYRN(6), NMORDR, NAORDR, NTORDR, NDCOOR, NMREDU, NUMCAL, MEMOFF, NMDPRP, &
              NMRDRP, NARDRP, NMPCAL, NMEPRP, LUNDPR, NMEXCI
   logical :: NWPROP, NOMOVE, NUMELC, CMPARE, SDRTST, RESTRT, DRYRUN, FRSTNM, PRPVIB, &
              MANUAL, ADSTBK, NOADSY, MDEANA, FSTPRP
   end module numder_mod
