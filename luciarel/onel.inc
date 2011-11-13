*
*   One-electron integrals
*
      integer   LEN,IAD2,IDIAPTYP,IDGGASOFF,
     &          IDIAPSMPTYP,IDIAPUBPSMPTYP,IOFFTRA,
     &          IPUBPGS,ISUMNEW
      dimension IPUBPGS(4,MXPNGAS,MXPNGAS,32)
      dimension IDIAPUBPSMPTYP(2,3:8,MXPNGAS)
      dimension IDIAPSMPTYP(3:8,MXPNGAS)
      dimension IDGGASOFF(MXPNGAS,2)
      dimension LEN(4),IAD2(3),IDIAPTYP(MXPNGAS)
      dimension IOFFTRA(32)
      common /ONELR/LEN,IAD2,IDIAPTYP,
     &             IDGGASOFF,IPUBPGS,IDIAPSMPTYP,
     &             IDIAPUBPSMPTYP,IOFFTRA,ISUMNEW
