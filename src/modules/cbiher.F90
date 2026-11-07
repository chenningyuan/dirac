!-----------------------------------------
! Module with variables previously defined
! at cbiher.h include file (COMMON /CBIHER/)
!-----------------------------------------
! (Control information for HERMIT)
   module cbiher_mod
   use mxcent_mod, only: MXCENT
   implicit none
   real*8 :: EXPKR(3), THRESH, PRTHRS
   integer :: IPRDEF, IORCAR, IORSPH, NPQUAD, NPATOM, IPATOM(MXCENT)
   logical :: SKIP, TSTINP, HAMILT, DIPLEN, SPNORB, SOTEST, SUPMAT, DIPVEL, QUADRU, SECMOM, &
              ONEPRP, CARMOM, SPHMOM, OCTUPO, FERMI, PSO, SPIDIP, DSO, ALLATM, NMRISS, &
              TRIANG, SDFC, PROPRI, HDO, S1MAG, S2MAG, ANGMOM, ANGLON, LONMOM, MAGMOM, &
              S1MAGT, MGMOMT, KINENE, S2MAGT, DSUSNL, DSUSLL, DSUSLH, DIASUS, DSUTST, &
              NUCSNL, NUCSLO, NUCSHI, NSNLTS, NSLTST, NELFLD, NSTTST, EFGCAR, EFGSPH, &
              NUCPOT, S1MAGL, S1MAGR, HDOBR, S1MLT, S1MRT, HDOBRT, NPOTST, MGMO2T, HBDO, &
              SUSCGO, NSTCGO, EXPIKR, DIRAC, DARWIN, MASSVL, CM1, CM2, NOTWO, SQHDOL, &
              SQHDOR, S1ELE, S1ELB, ONEELD, GFACDI, THETA, MGCOOR, CXIKR
   end module cbiher_mod
