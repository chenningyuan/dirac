C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbgen.h
C SOURCE  : $Source$
C
C*** General information in DIRAC ***
C  DOPIN - read MO-coefficients from formatted file and write them to
C          unformatted
C  DOPUT - write MO-coefficients to formatted file
C
C  INTGEN: Give top level integral flag
C          Should be pased on to every module
C
      PARAMETER ( LBASDIR = 600 )
      CHARACTER TITLE*50
      COMMON/CBCPAM/TITLE
      LOGICAL FAMILY,INPTES,VACUUM,NUCVAC,DOPSI,DOTRA,
     &        DOANA,DOPRP,DOHEINZ,DOHRM,DOJACO,DOPIN,DOPUT,
     &        DOCCM,DOCIM,OPTIMI,HRINPC, RDINPC, PARCAL, DIRCAL,
     &        SEGBAS,NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN
      COMMON/GNLPAM/FAMILY,INPTES,VACUUM,DOPSI,NUCVAC,DOTRA,
     &              DOPRP,DOANA,DOHEINZ,DOHRM,DOJACO,DOPIN,DOPUT,
     &              DOCCM,DOCIM,OPTIMI,
     &              HRINPC, RDINPC, PARCAL, DIRCAL,SEGBAS,
     &              NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN
      COMMON/GNIPAM/IPRGEN,ISPHTR,IDFLAG,LENBAS,IPRUSR,KCHARG,
     &              INTGEN
      COMMON/GNRPAM/CVAL,STOL(2),PANAS
      COMMON/PAMIOU/LUCOEF,LUOVLP,LUTMAT,LUPMAT,LU1INT
      CHARACTER BASDIR*(LBASDIR)
      COMMON /GNCPAM/ BASDIR
