C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcborb.h
C
C*** Dimensions for MO-basis ***
      PARAMETER( MXOPEN = 10 )
C     Length of common block : keep up to date !!
      PARAMETER (NDCBORB=79+5*MXOPEN)
      COMMON/DCBORB/
     +              DF(0:MXOPEN),DA(0:MXOPEN),DALPHA(0:MXOPEN),
     +              NORB(2),N2ORB(2),NORBT,N2ORBT,N2ORBX,
     +              N2ORBTQ,N2ORBXQ,
     +              NISH(2),N2ISH(2),NISHT,N2ISHT,N2ISHX,
     +              NASH(2),N2ASH(2),NASHT,N2ASHT,N2ASHX,
     +              NOCC(2),N2OCC(2),NOCCT,N2OCCT,N2OCCX,
     +              NSSH(2),N2SSH(2),NSSHT,N2SSHT,N2SSHX,
     +              NESH(2),N2ESH(2),NESHT,N2ESHT,N2ESHX,
     +              NPSH(2),N2PSH(2),NPSHT,N2PSHT,N2PSHX,
     +              NFRO(2),NFROT,NELEC(2),
     +              NAELEC,NELECT,
     +              IORB(2),I2ORBT(2),I2ORBX(2,2),
     +              NCMO(2),NCMOQ(2),NCMOT,NCMOTQ,ICMO(2),ICMOQ(2),
     +              NACSH(2,MXOPEN),NOPEN,
     +              IASH(2),I2ASHT(2),I2ASHX(2,2)
