!-----------------------------------------
! Module with variables previously defined
! at dcbgascip.h include file (COMMON /DCI_GASCIP/, /DCC_GASCIP/, /DC_GASCIP_ORB/)
!-----------------------------------------
!  DEFINITION PARAMETERS GASCIP :
!  MAXASH_GASCIP is the limit for the GASCIP CI program
!        (64 bit integers = 64 orbitals)
!  N5 is the max no of reps produced by GMULTSF (and GMULTA)
   module dcbgascip_mod
   implicit none
   integer, parameter :: MAXASH_GASCIP = 64
   integer, parameter :: N5 = 64
   integer :: IRRP, IRP, NREP, MULTB(N5,0:N5), MK2DEL, MK2REF, MINMK2, MAXMK2, &
              IRPAMO(2*MAXASH_GASCIP), NEL_STR_T, IEL_STR(0:MAXASH_GASCIP), &
              NEL_STR(0:MAXASH_GASCIP), IMK2(-MAXASH_GASCIP:MAXASH_GASCIP,3), &
              NMK2(-MAXASH_GASCIP:MAXASH_GASCIP,3), IPTA2O(MAXASH_GASCIP), &
              IPTB2O(MAXASH_GASCIP), IPTO2A(2*MAXASH_GASCIP), IPTO2B(2*MAXASH_GASCIP), &
              IGASCIP_NISH(2), IGASCIP_NASH(2), IGASCIP_NASHT, IGASCIP_REP, IGASCIP_STATE, &
              IGASCIP_NGAS, IGASCIP_NAELEC, IGASCIP_MK2REF, IGASCIP_MK2DEL, IGASCIP_MINMK2, &
              IGASCIP_MAXMK2, IGASCIP_NGSH(2,MAXASH_GASCIP), IGASCIP_NGASSP(2,MAXASH_GASCIP)
   logical :: GASCIP_UCIBOS
   character(len=4) :: REPNA(N5)
   end module dcbgascip_mod
