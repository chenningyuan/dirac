      COMMON/KRCC_CGAS/NGAS,NGSSH(MXPIRR,MXPNGAS),
     &            NGSSHK(MXPIRR,MXPNGAS,2),
     &            NGSSHT(MXPNGAS),NGSOB(MXPOBS,MXPNGAS),
     &            NGSOBT(MXPNGAS),NGSOBTK(MXPNGAS,2),
     &            IGSOCC(MXPNGAS,2),
     &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,
     &            NCMBSPC,LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),
     &            NMXOCCLS,MXNELGS(MXPNGAS),
     &            NGSOB2(MXPOBS,MXPNGAS),NGSOBT2(MXPNGAS),
! IBSH not really used
     &            IPHGAS(MXPNGAS),IBSH(MXPOBS),
     &            IHPVGAS(MXPNGAS),IPHGAS1(MXPNGAS),
     &            IHPVGAS_AB(MXPNGAS,2),
     &            IOFFINTTYPE(MXPSTT),IRELTYPE(20),ISUMORB(MXPORB),
     &            INOCOVA,KRCC_NISH(2),
     &            IREF_ELEC(MXPNGAS),IREF_AB_ELEC(2,MXPNGAS)
C NOTE > NGSSHT is the number of shells of given symmetry, not
C given type. Should be changed ( name and dimension some day )
C ( Someday never comes as C.C.R sings )
C last four lines added by Lasse when adding gencc to luciarel
C I think MXPSTT is correct for IOFFINTTYPE
