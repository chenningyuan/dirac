!-----------------------------------------
! Module with variables previously defined
! at nuclei.h include file (COMMON /NUCLEI/, /NUCLEC/)
!-----------------------------------------
   module nuclei_mod
   use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
   implicit none
   real*8 :: CHARGE(MXCENT), CORD(3,MXCENT), GNUEXP(MXCENT)
   integer :: NUCPRE(MXCENT), NUCNUM(MXCENT,8), NUCDEG(MXCENT), ISTBNU(MXCENT), NUCIND, &
              NUCDEP, NTRACO, ITRACO(3), NATOMS, NFLOAT, NBASIS, NLARGE, NSMALL, NPBAS, &
              NPLRG, NPSML, NUCNET, INCENT(MXCENT), INUNIQ(MXCENT), NDEGNM(MXCENT), &
              ISOTOP(MXCENT), IZATOM(MXCENT), NAMFA, IAMFA(MXCENT), IATOMTYP(MXCENT), &
              NONTYP, NONT(MXATOM), NBAS_ATOMTYPE(MXATOM,0:MXBSETS), NBASIS_FIT, &
              NLARGE_FIT, NSMALL_FIT, NPBAS_FIT, NPLRG_FIT, NPSML_FIT
   logical :: NOORBT(MXCENT), GAUNUC
   character(len=4) :: NAMN(MXCENT)
   character(len=6) :: NAMEX(MXCOOR), NAMDEP(MXCENT), NAMDPX(MXCOOR)
   end module nuclei_mod
