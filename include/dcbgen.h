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
C  INTGEN - Give top level integral flag
C  DOLVC  - Treat SS 2-e interaction by point charge model
C  DOTSC  - Treat SS 2-e interaction by TS point charge model
C  These flags should be pased on to every module
C
C
      PARAMETER ( LBASDIR = 600 )
      CHARACTER TITLE*50
      COMMON/CBCPAM/TITLE
      LOGICAL FAMILY,INPTES,VACUUM,NUCVAC,DOPSI,DOTRA,NOTRA,DOLUCT,
     &        DOANA,DOPRP,DOHEINZ,DOHRM,DOJACO,DOPIN,DOPUT,
     &        DOCCM,DOCIM,OPTIMI,HRINPC, RDINPC, PARCAL, DIRCAL,
     &        SEGBAS,NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN,
     &        DOLVC,LVNEW,DOTSC
      COMMON/GNLPAM/FAMILY,INPTES,VACUUM,DOPSI,NUCVAC,DOTRA,NOTRA,
     &              DOLUCT,DOPRP,DOANA,DOHEINZ,DOHRM,DOJACO,
     &              DOPIN,DOPUT,DOCCM,DOCIM,OPTIMI,
     &              HRINPC, RDINPC, PARCAL, DIRCAL,SEGBAS,
     &              NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN,
     &              DOLVC,LVNEW,DOTSC
      COMMON/GNIPAM/IPRGEN,ISPHTR,IDFLAG,LENBAS,IPRUSR,KCHARG,
     &              INTGEN
      COMMON/GNRPAM/CVAL,STOL(2),PANAS
      COMMON/PAMIOU/LUCOEF,LUOVLP,LUTMAT,LUPMAT,LU1INT
      CHARACTER BASDIR*(LBASDIR)
      COMMON /GNCPAM/ BASDIR
