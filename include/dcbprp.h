
!
! VERSION : $Revision$
! DATE    : $Date$
!
!*****************************************************************************
! *** Property module ***
!  DOEXP - do expectation values
!  DOESR - do ESR properties
!  DOXLR - do linear response
!  DOXPP - do excitation energies 
!          (by the polarization propagator method)
!  DOXQR - do quadratic response
!  DOTPA - do two-photon absortion calculation
!  DOEXCPRP - do excited state properties calculation
!  EXCPRPTYP - type of excited state property calculation
!  DOVER - do verdet calculations (via DOXQR)
!  MGRAD - calculate molecular gradient
!  DOPCT - do perform picture change transformation of four-component
!          property operators
!  NOPCT - do NOT perform the picture change transformation when DOPCT is on
!          - take only LL block of the four-component property operator
           
      INTEGER LUCORP
      PARAMETER (LUCORP = 91)

      LOGICAL DOEXP, DOXLR, DOXPP, DOXQR, DIPOLE, EPOLAR, MSUSCP,       &
     &        SHIELD, SPNSPN, EFG, NQCC, MGRAD, NSTDIA, PVC, DOESR,     &
     &        ESRGTENS, ESRANISO, EFT, EFTNTL, EFN, RHONUC,             &
     &        DOTPA,DOVER,DOEXCPRP,DOVDW,                               &
     &        ROTG, LONDON, SYMCON, NOORTH, DOPCT, MSUSCDIA, NOPCT,     & 
     &        PVC_SHIELD, PVC_SPINSPIN, DOSTEX, OPTROT,                 &
     &        BXLAO, BYLAO, BZLAO, EFFDEN

      INTEGER MXPRPWF,EXCPRPTYP
      PARAMETER (MXPRPWF = 4)

      CHARACTER*4 PRPWF(MXPRPWF)

      COMMON /CBLPRP/ DOXLR, DOXPP, DOEXP, DOXQR, DIPOLE, EPOLAR,DOVDW, &
     &                MSUSCP, SHIELD, SPNSPN, EFG, NQCC, MGRAD, NSTDIA, &
     &                PVC, DOESR, ESRGTENS, ESRANISO, EFT, EFTNTL,      &
     &                EFN, RHONUC, DOTPA, DOVER, DOEXCPRP, EXCPRPTYP,   &
     &                ROTG,LONDON,SYMCON,NOORTH,DOPCT,MSUSCDIA,NOPCT,   & 
     &                PVC_SHIELD,PVC_SPINSPIN,DOSTEX,OPTROT,            &
     &                BXLAO, BYLAO, BZLAO, EFFDEN

      INTEGER IPRPRP, IPEPOLAR(3, 2), IPSPNSPN(2, MXCOOR),              &
     &        IPDSO(MXCOOR*(MXCOOR+1)/2), NDSOIN, IPSHIELD(MXCOOR+3),   &
     &        IPNSTDIA(MXCOOR*3), NNSTIN, IPESR(MXCOOR+3, 2), NPRPWF,   &
     &        IPROTG(2,3),IPMSUSC(2,3),IMAGOUT,                         &
     &        IPCON21(3),IPCON22(3),IPCON23(6),                         &
     &        IPLONDON(4,3),IPCON(3),IPRNST(MXCOOR*3),                  &
     &        IP_PVC_SHIELD(MXCOOR+3),                                  &
     &        IP_PVC_SPINSPIN(2,MXCOOR),                                &
     &        IP_BLAO(2,3)
      COMMON /CBIPRP/ IPRPRP, IPEPOLAR, IPSPNSPN, IPDSO, NDSOIN,        &
     &                IPSHIELD, IPNSTDIA, NNSTIN, IPESR, NPRPWF,        &
     &                         IPROTG,IPMSUSC,IMAGOUT,                  &
     &                IPCON21,IPCON22,IPCON23,                          &
     &                IPLONDON,IPCON,IPRNST,                            &
     &                IP_PVC_SHIELD,                                    &
     &                IP_PVC_SPINSPIN,                                  &
     &                IP_BLAO

      REAL*8 ABUND
      COMMON /CBRPRP/ ABUND

      COMMON /CBCPRP/ PRPWF

#include <cbiher.h>
