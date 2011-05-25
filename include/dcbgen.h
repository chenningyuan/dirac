!
! VERSION : $Revision$
! DATE    : $Date$
! FILE    : dcbgen.h
!
!*** General information in DIRAC ***
!
!  DOPIN - read MO-coefficients from formatted file and write them to
!          unformatted
!  DOPIN_I3 - IBEIG in MO-coeff formatted file in I3 format (linear symmetry, large basis set)
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
!      TWOCOMP, TWOCOMPBSS - related to the pure two-component mode as follows:
!  TWOCOMP = T means that the program is currently running in 2c mode (either nonrelativistic or 2c-relativistic)
!  TWOCOMPBSS = T means that this is a 2-component relativistic calculation in 2c mode 
!  X2Cmod_x2c = T means that this is a 2-component relativistic calculation in 2c mode using the X2C module
!  (now stored in the module x2cmod_cfg)
!
!  SKIP2E - Skip all two electron steps.
!
!  DOJACO - do the Jacobi diagonalization of real matrixes
!  DOQJACO - Jacobi diagonalization of quaternion Jacobi matrixes
!
!  These flags should be passed on to every module
!
! NOSFCRASH - No Seg.Fault crash (by call quit). This is important for preserving gmon.out
!
      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)

      CHARACTER*50 TITLE
      COMMON /CBCPAM/ TITLE

      LOGICAL FAMILY, INPTES, VACUUM, NUCVAC, DOPSI, DOTRA, NOTRA,      &
     &        DOANA, DOPRP, DOHEINZ, DOHRM, DOJACO, DOPIN, DOPIN_I3,    & 
     &        DOPUT,                                                    &
     &        OPTIMI, HRINPC, RDINPC, PARCAL, DIRCAL, DOFCK3, SEGBAS,   &
     &        NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN, USRIPR, DOLVC,    &
     &        LVNEW, DOTSC, DIRSET, BNCRON, DOACUT, MP2ORG, DENSFIT,    &
     &        TWOCOMP, TWOCOMPBSS, DOLINK, DOFMM, NOSET,                &
     &        DOACIN, CBKRINI, SKIP2E, LOWJACO, NOSFCRASH, DOQJACO

      COMMON /GNLPAM/ FAMILY, INPTES, VACUUM, DOPSI, NUCVAC, DOTRA,     &
     &                NOTRA, DOPRP, DOANA, DOHEINZ, DOHRM, DOJACO,      &
     &                DOFCK3, DOPIN, DOPIN_I3,                          & 
     &                DOPUT, OPTIMI, HRINPC, RDINPC,                    &
     &                PARCAL, DIRCAL, SEGBAS, NEWBAS, NEWPRP, RELCAL,   &
     &                NEWSYM, RDMLIN, USRIPR, DOLVC, LVNEW, DOTSC,      &
     &                BNCRON, DOACUT, MP2ORG, DENSFIT, TWOCOMP,         &
     &                TWOCOMPBSS, DOLINK, DOFMM, NOSET,                 &
     &                DOACIN, CBKRINI, SKIP2E, LOWJACO, NOSFCRASH,      &
     &                DOQJACO

      INTEGER IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG, INTGEN,   &
     &        ILLDIR, ISLDIR, ISSDIR, IGTDIR, LENDND
      COMMON /GNIPAM/ IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG,   &
     &                INTGEN, ILLDIR, ISLDIR, ISSDIR, IGTDIR, DIRSET,   &
     &                LENDND

      REAL*8 CVAL, STOL(2), PANAS
      COMMON /GNRPAM/ CVAL, STOL, PANAS

      INTEGER LUCOEF, LUOVLP, LUTMAT, LUPMAT, LU1INT, LUKRMC, LUKRM1,   &
     &        LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO, LUBSS,  LUX2C
      COMMON /PAMIOU/ LUCOEF, LUOVLP, LUTMAT, LUPMAT, LU1INT, LUKRMC,   &
     &                LUKRM1, LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO,   &
     &                LUBSS,  LUX2C

      CHARACTER*(LBASDIR) BASDIR, CDIRNOD
      COMMON /GNCPAM/ BASDIR, CDIRNOD
