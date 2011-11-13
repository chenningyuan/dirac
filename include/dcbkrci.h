!
! FILE    : dcbkrci.h
!
!
!     *** information for KR-CI calculations ***
!
!
      LOGICAL KRCI_MKCNST, KRCI_UCIBOS, CRDFO_MAT_KRCI, CWRTFO_MAT_KRCI,&
     &        CSHMEMO_KRCI, CINT_REORD_KRCI, CINT_SPLIT_KRCI,           &
     &        NATOLCR_KRCI, CHCKPT_KRCI, CUSE_PBLFL_KRCI, GENFOCK_KRCI, &
     &        TRANATO_KRCI, LOWSORT_KRCI, STATE_SELECT_KRCI,            &
     &        domcana_krci
      CHARACTER*8 KRCI_CIPROGRAM
      COMMON /DCLKRCI/ KRCI_MKCNST, KRCI_UCIBOS, CRDFO_MAT_KRCI,        &
     &                 CWRTFO_MAT_KRCI, CSHMEMO_KRCI, CINT_REORD_KRCI,  &
     &                 CINT_SPLIT_KRCI, NATOLCR_KRCI, CHCKPT_KRCI,      &
     &                 CUSE_PBLFL_KRCI, GENFOCK_KRCI, TRANATO_KRCI,     &
     &                 LOWSORT_KRCI, STATE_SELECT_KRCI, domcana_krci

      INTEGER NKRCIISH(2), NKRCIASH(2), NKRCISSH(2),                    &
     &        NKRCIGSH(2, MXMCGAS), NKRCIGSP(2, MXMCGAS), NKRCIESH(2),  &
     &        NKRCIPSH(2), NKRCIFRO(2), NKRCIAELEC, NKRCIGAS,           &
     &        NKRCI_MK2REF, NKRCI_MK2DEL, NKRCI_MINMK2, NKRCI_MAXMK2
  
      INTEGER IKRCI_ITRINT(2), IKRCI_INTDEF, IKRCI_INTFLG, IKRCI_INTBUF,&
     &        NKRCI_MAX_SYM, MAX_NKRCI_MAX_SYM, NKRCI_CIROOTS, IPRKRCI, &
     &        MAXCIT_KRCI, MXCIV_KRCI,               IANACI_KRCI,       &
     &        KTRLVL_KRCI, IRESTRK_KRCI, IDOPARIO_KRCI, IKRCI_MEMFAC,   &
     &        NKRCI_SELOM
      PARAMETER (MAX_NKRCI_MAX_SYM = 128)
      COMMON /DCIKRCI/ NKRCIISH, NKRCIASH, NKRCISSH, NKRCIGSH, NKRCIGSP,&
     &                 NKRCIESH, NKRCIPSH, NKRCIFRO, NKRCIAELEC,        &
     &                 NKRCIGAS, NKRCI_MK2REF, NKRCI_MK2DEL,            &
     &                 NKRCI_MINMK2, NKRCI_MAXMK2,                      &
     &                 IKRCI_ITRINT, IKRCI_INTDEF, IKRCI_INTFLG,        &
     &                 IKRCI_INTBUF, NKRCI_MAX_SYM,                     &
     &                 NKRCI_CIROOTS(MAX_NKRCI_MAX_SYM), IPRKRCI,       &
     &                 MAXCIT_KRCI, MXCIV_KRCI,                         &
     &                 IANACI_KRCI, KTRLVL_KRCI, IRESTRK_KRCI,          &
     &                 IDOPARIO_KRCI, IKRCI_MEMFAC,                     &
     &                 NKRCI_SELOM(MAX_NKRCI_MAX_SYM)
!
      REAL*8 DKRCI_CNVINT(2), DKRCI_THRPCI
      COMMON /DCRKRCI/ DKRCI_CNVINT, DKRCI_THRPCI

      COMMON /DCCKRCI/ KRCI_CIPROGRAM
!
! --- end of dcbkrci.h ---
