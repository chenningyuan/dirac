C *** Property module ***
C  DOEXP - do expectation values
C  DOXLR - do linear response
C  DOXPP - do excitation energies 
C          (by the polarization propagator method)
C  MGRAD - calculate molecular gradient
      LOGICAL DOEXP,DOXLR,DOXPP,DIPOLE,EPOLAR,MSUSCP,SHIELD,
     &        SPNSPN,EFG,NQCC,MGRAD,NSTDIA,PVC
      COMMON/CBLPRP/DOXLR,DOXPP,DOEXP,DIPOLE,EPOLAR,MSUSCP,
     &        SHIELD,SPNSPN,EFG,NQCC,MGRAD,NSTDIA,PVC
      COMMON/CBIPRP/IPRPRP,IPEPOLAR(3,2),IPSPNSPN(2,MXCOOR),
     &     IPDSO(MXCOOR*(MXCOOR+1)/2),NDSOIN,IPSHIELD(MXCOOR+3),
     &              IPNSTDIA(MXCOOR*3),NNSTIN
      COMMON/CBRPRP/ABUND
#include <cbiher.h>
