      CHARACTER*72 STARS,RELEAS
#if defined (VAR_MOTECC)
      CHARACTER*12 VERSN
#else
      CHARACTER*3  VERSN
#endif
      DATA STARS(1:36)  /'************************************'/
      DATA STARS(37:72) /'************************************'/
      DATA RELEAS(1:36) /'*SIRIUS* a direct, restricted step, '/
      DATA RELEAS(37:72)/'second order MCSCF program  *Dec 94*'/
#if defined (VAR_MOTECC)
      DATA VERSN        /'METECC-94.n06'/
#else
      DATA VERSN        /'n06'/
#endif
