C     Common block for LUCITA input information
C     read in dirrdn.F/PSIINP/LUCTINP
C     Information is transferred to LUCITA routine dirluct.F
C     for further processing.
C
      PARAMETER (NTABLE = 18)
      PARAMETER (MXNGAS = 16)
C
      CHARACTER*72 WAFFCD,CALCTP,SZCALD,TITLUC,
     &             CRDINA,CRDGAS(MXNGAS),CRDGOC(MXNGAS),CRDFRO,
     &             CRDRS1,CRDRS2,CRDRS3
C
      COMMON /LUCTINFC/ CRDINA,CRDGAS,CRDGOC,CRDFRO,CRDRS1,CRDRS2,
     &                  CRDRS3,WAFFCD,CALCTP,SZCALD,TITLUC
      COMMON /LUCTINFI/ IMOKW,NROOTD,ISSYMD,NACTED,IMULTD,IPRNGD,
     &                  IPRNLD,IDENSD,MXHL1D,MXEL3D,INGASD,NSEQCD
