!
! VERSION : $Revision: 1.15 $
! DATE    : $Date: 2004/06/16 13:56:48 $
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
!  DOVER - do verdet calculations (via DOXQR)
!  MGRAD - calculate molecular gradient
      INTEGER LUCORP
      PARAMETER (LUCORP = 91)

      LOGICAL DOEXP, DOXLR, DOXPP, DOXQR, DIPOLE, EPOLAR, MSUSCP,       &
     &        SHIELD, SPNSPN, EFG, NQCC, MGRAD, NSTDIA, PVC, DOESR,     &
     &        ESRGTENS, ESRANISO, EFT, EFTNTL, PVCSHI, EFN, RHONUC,     &
     &        DOTPA,DOVER

      INTEGER MXPRPWF
      PARAMETER (MXPRPWF = 3)

      CHARACTER*4 PRPWF(MXPRPWF)

      COMMON /CBLPRP/ DOXLR, DOXPP, DOEXP, DOXQR, DIPOLE, EPOLAR,       &
     &                MSUSCP, SHIELD, SPNSPN, EFG, NQCC, MGRAD, NSTDIA, &
     &                PVC, DOESR, ESRGTENS, ESRANISO, EFT, EFTNTL,      &
     &                PVCSHI, EFN, RHONUC, DOTPA, DOVER

      INTEGER IPRPRP, IPEPOLAR(3, 2), IPSPNSPN(2, MXCOOR),              &
     &        IPDSO(MXCOOR*(MXCOOR+1)/2), NDSOIN, IPSHIELD(MXCOOR+3),   &
     &        IPNSTDIA(MXCOOR*3), NNSTIN, IPESR(MXCOOR+3, 2), NPRPWF
      COMMON /CBIPRP/ IPRPRP, IPEPOLAR, IPSPNSPN, IPDSO, NDSOIN,        &
     &                IPSHIELD, IPNSTDIA, NNSTIN, IPESR, NPRPWF

      REAL*8 ABUND
      COMMON /CBRPRP/ ABUND

      COMMON /CBCPRP/ PRPWF

#include <cbiher.h>
