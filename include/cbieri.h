C
C     Parameters NCBI? must be updated after changes (for parallelization)
C
C     NOTE: New logicals should appear after the last logical (NCLERI)
C           New integers should appear after the last integer (IFITDM)
C           Reals should appear at the end.
C
      INTEGER NCBII,  NCBIL
      PARAMETER (NCBII = 18,
     &           NCBIL = 30)
      INTEGER IPRERI, IAOBCH, IPRNT1, IPRNT2, IPROD1, IPROD2
      INTEGER LBFINP, MAXDST, NDMAT,  IANGMO, NSPMAX, MAXDSD,
     &        MXBCH,  MXBCH0, IFITDM
      LOGICAL RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI, OFFCNT,
     &        DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT, EXTPRI, INTPRI,
     &        DODIST, INTSKP, DISTST, WRTINT, FCKINT, EXCINT, PROFIL,
     &        NOLOCS, GRDZER, OLDDER, EXPERI, DOERIP, ERITWO, CCRUN,
     &        COMPRS, GENCON, NCLERI, DGABAB, BFCFMM, ATTGAM, ATTEXC
      DIMENSION IANGMO(4)
C
CC    INTEGER NCENTS, ICENTS, ICLASSES, KOD
CC    DIMENSION NCENTS(4), ICENTS(2,4),ICLASSES(4), KOD(10)
      INTEGER NCENTS, ICENTS, ICLASSES
      DIMENSION NCENTS(4), ICENTS(2,4),ICLASSES(4)
C
      COMMON /CBIERI/ RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI,
     &                OFFCNT, DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT,
     &                EXTPRI, INTPRI, DODIST, INTSKP, DISTST, WRTINT,
     &                FCKINT, EXCINT, PROFIL, NOLOCS, GRDZER, OLDDER,
     &                EXPERI, DOERIP, ERITWO, CCRUN,  COMPRS, GENCON,
     &                NCLERI, DGABAB, BFCFMM, ATTGAM, ATTEXC,
     &                NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1, IPRNT2,
     &                IPRERI, LBFINP, MAXDST, IANGMO, NSPMAX, MAXDSD,
     &                MXBCH,  MXBCH0, IFITDM,
     &                NCENTS, ICENTS, ICLASSES
