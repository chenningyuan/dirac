!-----------------------------------------
! Module with variables previously defined
! at dcbgen.h include file (COMMON /CBCPAM/, /DCBGEN_L_PAM/, /DCBGEN_I_PAM/, /DCBGEN_R_PAM/, /PAMIOU/, /DCBGEN_C_PAM/)
!-----------------------------------------
!*** General information in DIRAC ***
!  INTGEN - Give top level integral flag
!  DOLVC  - Treat SS 2-e interaction by point charge model
!  DOTSC  - Treat SS 2-e interaction by TS point charge model
!  DIRSET - Conventional OR direct calculation of 2-e integrals
!           has been specified manually (used for SMLV1C)
!  BNCRON - Improved screening of nucleus for bare nucleus starting guess
!  force_bnc - use input value of BNCRON w/o checking molecular charge (PSB+MI)
!  BNCRON - Use sum of atomic potentials for the starting guess
!  TWOCOMP, TWOCOMPBSS - related to the pure two-component mode as follows:
!    TWOCOMP = T means that the program is currently running in 2c mode (either nonrelativistic or 2c-relativistic)
!    TWOCOMPBSS = T means that this is a 2-component relativistic calculation in 2c mode
!    X2Cmod_x2c = T means that this is a 2-component relativistic calculation in 2c mode using the X2C module
!      (now stored in the module x2cmod_cfg)
!  DOJACO  - do the Jacobi diagonalization of real matrixes
!  DOQJACO - Jacobi diagonalization of quaternion Jacobi matrixes
!  These flags should be passed on to every module
! DIRAC_WF(1:N_WF) the wave function jobs to be done in this dirac run, set in dirac/dirrdn.F
   module dcbgen_mod
   implicit none
   integer, parameter :: LBASDIR = 600
   real*8 :: CVAL, STOL(2), PANAS
   integer :: IPRGEN, ISPHTR, IDFLAG, LENBAS, IPRUSR, KCHARG, INTGEN, ILLDIR, ISLDIR, &
              ISSDIR, IGTDIR, LENDND, N_WF, INTGEN_SAVE, LUCOEF, LUOVLP, LUTMAT, LUPMAT, &
              LU1INT, LUKRMC, LUKRM1, LUKRM2, LUKRM3, LUKRM4, LUKRM5, LUITFO, LUBSS, LUX2C
   logical :: INPTES, VACUUM, DOPSI, NUCVAC, DOTRA, NOTRA, DOPRP, DOANA, DOHRM, DOJACO, &
              DOFCK3, DORKBIMP, DOPUT, OPTIMI, HRINPC, RDINPC, PARCAL, DIRCAL, SEGBAS, &
              NEWBAS, NEWPRP, RELCAL, NEWSYM, RDMLIN, USRIPR, DOLVC, LVNEW, DOTSC, BNCRON, &
              DOACUT, MP2ORG, TWOCOMP, TWOCOMPBSS, DOACIN, CBKRINI, LOWJACO, DOQJACO, QM3, &
              QMMM, BARNUC, DOHUCKEL, OPTWLK, OPTNEW, NMWALK, FORCE_BNC, BNSPON, DIRSET
   character(len=10) :: DIRAC_WF(100)
   character(len=50) :: TITLE
   character(len=LBASDIR) :: BASDIR, CDIRNOD
   end module dcbgen_mod
