C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbgen.h
C
C*** General information in DIRAC ***
C  DOPIN - read MO-coefficients from formatted file and write them to
C          unformatted
C  DOPUT - write MO-coefficients to formatted file
C
C  INTGEN - Give top level integral flag
C  DOLVC  - Treat SS 2-e interaction by point charge model
C  DOTSC  - Treat SS 2-e interaction by TS point charge model
C  DIRSET - Conventional OR direct calculation of 2-e integrals 
C           has been specified manually (used for SMLV1C)
C  BNCRON - Improved screening of nucleus for bare nucleus starting guess
C
C  These flags should be pased on to every module
C
C
      PARAMETER ( LBASDIR = 600 )
      CHARACTER TITLE*50
      COMMON/CBCPAM/TITLE
      LOGICAL FAMILY,INPTES,VACUUM,NUCVAC,DOPSI,DOTRA,NOTRA,
     &        DOANA,DOPRP,DOHEINZ,DOHRM,DOJACO,DOPIN,DOPUT,
     &        OPTIMI,HRINPC, RDINPC, PARCAL, DIRCAL,
     &        SEGBAS,NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN, USRIPR,
     &        DOLVC,LVNEW,DOTSC,DIRSET,BNCRON,DOACUT
      COMMON/GNLPAM/FAMILY,INPTES,VACUUM,DOPSI,NUCVAC,DOTRA,NOTRA,
     &              DOPRP,DOANA,DOHEINZ,DOHRM,DOJACO,
     &              DOPIN,DOPUT,OPTIMI,
     &              HRINPC, RDINPC, PARCAL, DIRCAL,SEGBAS,
     &              NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN, USRIPR,
     &              DOLVC,LVNEW,DOTSC,BNCRON,DOACUT
      COMMON/GNIPAM/IPRGEN,ISPHTR,IDFLAG,LENBAS,IPRUSR,KCHARG,
     &              INTGEN,ILLDIR,ISLDIR,ISSDIR,DIRSET
      COMMON/GNRPAM/CVAL,STOL(2),PANAS
      COMMON/PAMIOU/LUCOEF,LUOVLP,LUTMAT,LUPMAT,LU1INT,
     &              LUKRMC,LUKRM1,LUKRM2,LUKRM3,LUKRM4,LUKRM5,
     &              LUITFO
      CHARACTER BASDIR*(LBASDIR)
      COMMON /GNCPAM/ BASDIR
