C
C VERSION : $Revision$
C DATE    : $Date$
C SOURCE  : $Source$
C
C*****************************************************************************
C *** Property module ***
C  DOEXP - do expectation values
C  DOESR - do ESR properties
C  DOXLR - do linear response
C  DOXPP - do excitation energies 
C          (by the polarization propagator method)
C  DOXQR - do quadratic response
C  MGRAD - calculate molecular gradient
      LOGICAL DOEXP,DOXLR,DOXPP,DOXQR,DIPOLE,EPOLAR,MSUSCP,SHIELD,
     &        SPNSPN,EFG,NQCC,MGRAD,NSTDIA,PVC,DOESR,ESRGTENS,ESRANISO,
     &        EFT,EFTNTL
      COMMON/CBLPRP/DOXLR,DOXPP,DOEXP,DOXQR,DIPOLE,EPOLAR,MSUSCP,
     &        SHIELD,SPNSPN,EFG,NQCC,MGRAD,NSTDIA,PVC,
     &        DOESR,ESRGTENS,ESRANISO,EFT,EFTNTL
      COMMON/CBIPRP/IPRPRP,IPEPOLAR(3,2),IPSPNSPN(2,MXCOOR),
     &     IPDSO(MXCOOR*(MXCOOR+1)/2),NDSOIN,IPSHIELD(MXCOOR+3),
     &              IPNSTDIA(MXCOOR*3),NNSTIN,IPESR(MXCOOR+3)
      COMMON/CBRPRP/ABUND
#include <cbiher.h>
