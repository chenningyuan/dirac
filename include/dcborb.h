C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcborb.h
C
C*** Dimensions for MO-basis ***
C     Length of common block : keep up to date !!
#include "mxmcgas.h"
      PARAMETER (NDCBORB=129+4*MXMCGAS)
      COMMON/DCBORB/
     &     NORB(2),N2ORB(2),NORBT,N2ORBT,N2ORBX,
     &     NBORB(4,2,0:2),NFORB(2,0:2),NTORB(0:2),             
     &     N2ORBTQ,N2ORBXQ,
     &     NISH(2),N2ISH(2),NISHT,N2ISHT,N2ISHX,
     &     NASH(2),N2ASH(2),NASHT,N2ASHT,N2ASHX,
     &     N2ASHTQ,N2ASHXQ,NNASHX,
     &     NGSH(2,MXMCGAS),NGSP(2,MXMCGAS),NGAS,
     &     NGSHT(MXMCGAS),
     &     NOCC(2),N2OCC(2),NOCCT,N2OCCT,N2OCCX,
     &     NSSH(2),N2SSH(2),NSSHT,N2SSHT,N2SSHX,
     &     NESH(2),N2ESH(2),NESHT,N2ESHT,N2ESHX,
     &     NPSH(2),N2PSH(2),NPSHT,N2PSHT,N2PSHX,
     &     NFRO(2),NFROT,NELEC(2),
     &     NAELEC,NELECT,
     &     IORB(2),I2ORBT(2),I2ORBX(2,2),
     &     NCMO(2),NCMOQ(2),NCMOT,NCMOTQ,ICMO(2),ICMOQ(2),
     &     IASH(2),I2ASHT(2),I2ASHX(2,2),
     &     IOCC(2),I2OCCT(2),I2OCCX(2,2),
     &     NAOCC(2),NAOCCT,
     &     NIOCC(2),NIOCCT,
     &     NAVIR(2),NAVIRT,
     &     NIVIR(2),NIVIRT
C     Information on orbital subblocks. Used at the moment only for linear symmetry
C     but can be generalized to eliminate couplings between orbital spaces.
      PARAMETER (MAX_SUB_BL=16)
      COMMON/DCBSUBS/
     & N_SUB_BL(2),ID_SUB_BL(MAX_SUB_BL,2),NORB_SUB(MAX_SUB_BL,2,0:2)
