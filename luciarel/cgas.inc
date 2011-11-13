      COMMON/CGASR/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
     &            NGSSHT(MXPNGAS),NGSOB(MXPOBS,MXPNGAS),
     &            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,
     &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,
     &            NCMBSPC,LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),
     &            NMXOCCLS,MXNELGS(MXPNGAS),
     &            NGSOB2(MXPOBS,MXPNGAS),NGSOBT2(MXPNGAS),
     &            IPHGAS(MXPNGAS),IBSH(MXPOBS)
C NOTE > NGSSHT is the number of shells of given symmetry, not
C given type. Should be changed ( name and dimension some day )
C ( Someday never comes as C.C.R sings )
