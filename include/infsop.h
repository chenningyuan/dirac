C
C     Variables named ?1/2/3??? are used for triplet properties
C     calculations. Others are for singlet properties.
C
      INTEGER ABSADR, ABTADR,AB1ADR,AB2ADR,AB3ADR
      LOGICAL A2EXIST
      COMMON /INFSOP/ NSOO(8), NTOO(8), NSVV(8), NTVV(8),
     *                NS2P2H(8), NT2P2H(8), N2P2HT(8),N2P2HS(8),
     *                N12P2H(8), N22P2H(8), N32P2H(8),
     *                ABSADR(MAXORB,MAXORB), ABTADR(MAXORB,MAXORB),
     *                IJSADR(MAXISH,MAXISH,8), IJTADR(MAXISH,MAXISH,8),
     *                IJ1ADR(MAXISH,MAXISH,8), IJ2ADR(MAXISH,MAXISH,8),
     *                IJ3ADR(MAXISH,MAXISH,8),
     *                A2EXIST
C
