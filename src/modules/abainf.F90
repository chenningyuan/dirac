!-----------------------------------------
! Module with variables previously defined
! at abainf.h include file (COMMON /ABAINF/)
!-----------------------------------------
! (Control information for ABACUS routines)
   module abainf_mod
   implicit none
   integer :: IPRDEF, NWNABA, NTCSY(8), NSCSY(8)
   logical :: MOLGRD, MOLHES, DIPDER, POLAR, TSTINP, VIB, RESTAR, DOWALK, GDALL, CCSD, &
              H2MO, DOSYM(8), DOLRES, DOEXCI, SHIELD, SPNSPN, MAGSUS, VCD, NACME, AAT, &
              NOLOND, FCKDDR, ECD, DODRCT, SUPMAT, PARALL, ROTG, SPINRO, MASSVE, DARWIN, &
              ABALNR, VROA, NODIFC, NOCMC, EXPFCK, RAMAN, QUADRU, NQCC, HYPER, VERDET, MCD, &
              GSDIP, GSQUAD, GSOCT, GSDIDI, GSDIQU, GSQUQU, LINCPL, NUMHES
   end module abainf_mod
