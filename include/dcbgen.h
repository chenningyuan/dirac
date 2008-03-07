!
! VERSION : $Revision$
! DATE    : $Date$
! FILE    : dcbgen.h
!
!*** General information in DIRAC ***
!  DOPIN - read MO-coefficients from formatted file and write them to
!          unformatted
!  DOPUT - write MO-coefficients to formatted file
!
!  INTGEN - Give top level integral flag
!  DOLVC  - Treat SS 2-e interaction by point charge model
!  DOTSC  - Treat SS 2-e interaction by TS point charge model
!  DIRSET - Conventional OR direct calculation of 2-e integrals 
!           has been specified manually (used for SMLV1C)
!  BNCRON - Improved screening of nucleus for bare nucleus starting 
!           guess
!
!  TWOCOMP, TWOCOMPBSS - related to the pure two-component mode
!
!  These flags should be passed on to every module
!
!
      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)

      CHARACTER*50 TITLE
      COMMON /CBCPAM/ TITLE

      LOGICAL FAMILY, INPTES, VACUUM, NUCVAC, DOPSI, DOTRA, NOTRA,      &
     &        DOANA, DOPRP, DOHEINZ, DOHRM, DOJACO, DOPIN, DOPUT,       &
     &        OPTIMI, HRINPC, RDINPC, PARCAL, DIRCAL, DOFCK3, SEGBAS,   &
     &        NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN, USRIPR, DOLVC,    &
     &        LVNEW, DOTSC, DIRSET, BNCRON, DOACUT, MP2ORG, DENSFIT,    &
     &        TWOCOMP, TWOCOMPBSS, DOLINK, DOFMM, NOSET, ECPCALC,       &
     &        DOACIN

      COMMON /GNLPAM/ FAMILY, INPTES, VACUUM, DOPSI, NUCVAC, DOTRA,     &
     &                NOTRA, DOPRP, DOANA, DOHEINZ, DOHRM, DOJACO,      &
     &                DOFCK3, DOPIN, DOPUT, OPTIMI, HRINPC, RDINPC,     &
     &                PARCAL, DIRCAL, SEGBAS, NEWBAS, NEWPRP, RELCAL,   &
     &                NEWSYM, RDMLIN, USRIPR, DOLVC, LVNEW, DOTSC,      &
     &                BNCRON, DOACUT, MP2ORG, DENSFIT, TWOCOMP,         &
     &                TWOCOMPBSS, DOLINK, DOFMM, NOSET, ECPCALC,        &
     &                DOACIN       

      INTEGER IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG, INTGEN,   &
     &        ILLDIR, ISLDIR, ISSDIR, IGTDIR
      COMMON /GNIPAM/ IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG,   &
     &                INTGEN, ILLDIR, ISLDIR, ISSDIR, IGTDIR, DIRSET

      REAL*8 CVAL, STOL(2), PANAS
      COMMON /GNRPAM/ CVAL, STOL, PANAS

      INTEGER LUCOEF, LUOVLP, LUTMAT, LUPMAT, LU1INT, LUKRMC, LUKRM1,   &
     &        LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO, LUBSS
      COMMON /PAMIOU/ LUCOEF, LUOVLP, LUTMAT, LUPMAT, LU1INT, LUKRMC,   &
     &                LUKRM1, LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO,   &
     &                LUBSS

      CHARACTER*(LBASDIR) BASDIR
      COMMON /GNCPAM/ BASDIR
