!
! FILE: include/nuclei.h
!
      CHARACTER*4 NAMN(MXCENT)
      LOGICAL NOORBT(MXCENT), GAUNUC
      REAL*8 CHARGE(MXCENT), CORD(3, MXCENT), GNUEXP(MXCENT)
      INTEGER NUCPRE(MXCENT), NUCNUM(MXCENT, 8), NUCDEG(MXCENT),        &
     &    ISTBNU(MXCENT), NUCIND, NUCDEP, NTRACO, ITRACO(3), NATOMS,    &
     &    NFLOAT, NBASIS, NLARGE, NSMALL, NPBAS, NPLRG, NPSML,          &
     &    NUCNET, INCENT(MXCENT), INUNIQ(MXCENT), NDEGNM(MXCENT),       &
     &    NAMFA, IAMFA(MXCENT), IATOMTYP(MXCENT),                       &
     &    NONTYP, NONT(MXATOM), NBAS_ATOMTYPE(MXATOM,0:2),              & ! <--- atom type info
     &    NBASISAUX, NPBASAUX, NAUX, NPAUX,                             & ! <--- density fitting basis info
     &    NBASIS_FIT, NLARGE_FIT, NSMALL_FIT, NPBAS_FIT, NPLRG_FIT,     &
     &    NPSML_FIT

      COMMON /NUCLEI/ CHARGE, CORD, NOORBT, NUCPRE, NUCNUM, NUCDEG,     &
     &    ISTBNU, GNUEXP, NAMN,   NUCIND, NUCDEP, NTRACO,               &
     &    ITRACO, NATOMS, NFLOAT, NBASIS, NLARGE, NSMALL,               &
     &    NPBAS,  NPLRG,  NPSML,  NUCNET, INCENT, INUNIQ,               &
     &    GAUNUC, NDEGNM, NAMFA,  IAMFA,  IATOMTYP,                     &
     &    NONTYP, NONT        , NBAS_ATOMTYPE            ,              & ! <--- atom type info
     &    NBASISAUX, NPBASAUX, NAUX, NPAUX,                             & ! <--- density fitting basis info
     &    NBASIS_FIT, NLARGE_FIT, NSMALL_FIT, NPBAS_FIT, NPLRG_FIT,     &
     &    NPSML_FIT

      CHARACTER*6 NAMEX(MXCOOR), NAMDEP(MXCENT), NAMDPX(MXCOOR)
      COMMON /NUCLEC/ NAMEX, NAMDEP, NAMDPX
! --- end of nuclei.h ---
