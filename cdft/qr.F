C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SDFT_QR */
      SUBROUTINE SDFT_QR(DTRMAT,
     &                   EXCMAT,
     &                   ISYMOP,
     &                   IHRMOP,
     &                   GAO,
     &                   GAO1,
     &                   DENSITY,
     &                   DENSITY_GRADIENT_NORM,
     &                   DENSITY_GRADIENT,
     &                   WEIGHT,
     &                   DOGGA,
     &                   BUF)
C=======================================================================
C     describe me ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                                last revision: may 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "inforb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dftcom.h"
#include "cdft.h"
#include "../dft/functionals.h"

      LOGICAL   DOGGA,
     &          LBIT
      DIMENSION DTRMAT(NBAST,NBAST,NZ,3),
     &          EXCMAT(NBAST,NBAST,NZ,3),
     &          ISYMOP(3),
     &          IHRMOP(3),
     &          GAO(NBAST),
     &          GAO1(NBAST,3),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          FO_DENSITY(3),
     &          FO_SPIN_DENSITY(3,3),
     &          FO_SPIN_DENSITY_NORM(3),
     &          FACTOR(3),
     &          ZETA(3),
     &          VXC_S(3),
     &          VXC_V(3,3),
     &          VXC(99) !dirty fix

#include "ibtfun.h"

      DO IBC = 1,3
        CALL GETRHO(FO_DENSITY(IBC),
     &              ISYMOP(IBC)-1,
     &              GAO,
     &              DTRMAT(1,1,1,IBC),
     &              BUF)

        CALL GET_M(.TRUE.,
     &             FO_SPIN_DENSITY(1,IBC),
     &             ISYMOP(IBC)-1,
     &             GAO,
     &             UNUSED,
     &             DTRMAT(1,1,1,IBC),
     &             BUF)

      ENDDO

      CALL SDFTPOT2(VXC,WEIGHT,DENSITY,D0)

      VXC(fSS)  = VXC(fSS)/D2
      VXC(fRSS) = VXC(fRSS)/D2

      DO IBC = 1,3
        DO IXYZ = 1,3
          VXC_V(IXYZ,IBC) = D2*VXC(fSS)*FO_SPIN_DENSITY(IXYZ,IBC)
        ENDDO
      ENDDO

      DO IXYZ = 1,3
        ADDME = VXC(fRSS)*FO_DENSITY(1)*FO_SPIN_DENSITY(IXYZ,2)
     &        + VXC(fRSS)*FO_DENSITY(2)*FO_SPIN_DENSITY(IXYZ,1)

        VXC_V(IXYZ,3) = VXC_V(IXYZ,3) + D4*ADDME
      ENDDO

      DO IBC = 1,3
        CALL GET_M(.FALSE.,
     &             VXC_V(1,IBC),
     &             ISYMOP(IBC)-1,
     &             GAO,
     &             UNUSED,
     &             EXCMAT(1,1,1,IBC),
     &             BUF)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SDFT_QR_ALTERNATIVE */
      SUBROUTINE SDFT_QR_ALTERNATIVE(ADD_R,
     &                               ADD_S,
     &                               ADD_GR,
     &                               ADD_GS,
     &                               DTRMAT,
     &                               EXCMAT,
     &                               ISYMOP,
     &                               IHRMOP,
     &                               GAO,
     &                               GAO1,
     &                               R,
     &                               GR_NORM,
     &                               GR,
     &                               W,
     &                               BUF)
C=======================================================================
C     describe me ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                               last revision: june 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "cdft.h"
#include "dftcom.h"
#include "../dft/functionals.h"

      LOGICAL   ADD_R,
     &          ADD_S,
     &          ADD_GR,
     &          ADD_GS,
     &          LBIT
      DIMENSION DTRMAT(NBAST,NBAST,NZ,3),
     &          EXCMAT(NBAST,NBAST,NZ,3),
     &          ISYMOP(3),
     &          IHRMOP(3),
     &          GAO(*),
     &          GAO1(*),
!
     &          GR(3),
     &          GR_b(3),
     &          GR_c(3),
     &          GR_bc(3),
!
     &          GS_b_ncol(3,3),
     &          GS_c_ncol(3,3),
     &          GS_bc_ncol(3,3),
!
     &          GS(3),
     &          GS_b(3),
     &          GS_c(3),
     &          GS_bc(3),
!
     &          S_b_ncol(3),
     &          S_c_ncol(3),
     &          S_bc_ncol(3),
!
     &          BUF(*),
     &          VXC_VECTOR(3),
     &          VXC(99) !dirty fix

#include "ibtfun.h"

      CALL DZERO(GR_b ,3)
      CALL DZERO(GR_c ,3)
      CALL DZERO(GR_bc,3)

      CALL DZERO(GS_b_ncol ,9)
      CALL DZERO(GS_c_ncol ,9)
      CALL DZERO(GS_bc_ncol,9)

      CALL DZERO(GS        ,3)
      CALL DZERO(GS_b      ,3)
      CALL DZERO(GS_c      ,3)
      CALL DZERO(GS_bc     ,3)

      CALL DZERO(S_b_ncol  ,3)
      CALL DZERO(S_c_ncol  ,3)
      CALL DZERO(S_bc_ncol ,3)

      IF(ADD_R) THEN
        ibc = 1
        CALL GETRHO(R_b,
     &              ISYMOP(ibc)-1,
     &              GAO,
     &              DTRMAT(1,1,1,ibc),
     &              BUF)
        ibc = 2
        CALL GETRHO(R_c,
     &              ISYMOP(ibc)-1,
     &              GAO,
     &              DTRMAT(1,1,1,ibc),
     &              BUF)
        ibc = 3
        CALL GETRHO(R_bc,
     &              ISYMOP(ibc)-1,
     &              GAO,
     &              DTRMAT(1,1,1,ibc),
     &              BUF)
      ENDIF

      IF(ADD_GR) THEN
        ibc = 1
        CALL DENSITY_GRADIENT_GENERAL(GR_b,
     &                                ISYMOP(ibc)-1,
     &                                GAO,
     &                                GAO1,
     &                                DTRMAT(1,1,1,ibc),
     &                                BUF)
        ibc = 2
        CALL DENSITY_GRADIENT_GENERAL(GR_c,
     &                                ISYMOP(ibc)-1,
     &                                GAO,
     &                                GAO1,
     &                                DTRMAT(1,1,1,ibc),
     &                                BUF)
        ibc = 3
        CALL DENSITY_GRADIENT_GENERAL(GR_bc,
     &                                ISYMOP(ibc)-1,
     &                                GAO,
     &                                GAO1,
     &                                DTRMAT(1,1,1,ibc),
     &                                BUF)
      ENDIF

      IF(ADD_S) THEN
        ibc = 1
        CALL GET_M(.TRUE.,
     &             S_b_ncol,
     &             ISYMOP(ibc)-1,
     &             GAO,
     &             UNUSED,
     &             DTRMAT(1,1,1,ibc),
     &             BUF)
        ibc = 2
        CALL GET_M(.TRUE.,
     &             S_c_ncol,
     &             ISYMOP(ibc)-1,
     &             GAO,
     &             UNUSED,
     &             DTRMAT(1,1,1,ibc),
     &             BUF)
        ibc = 3
        CALL GET_M(.TRUE.,
     &             S_bc_ncol,
     &             ISYMOP(ibc)-1,
     &             GAO,
     &             UNUSED,
     &             DTRMAT(1,1,1,ibc),
     &             BUF)
      ENDIF

      IF(ADD_GS) THEN
        ibc = 1
        CALL GET_M_DERVS(.TRUE.,
     &                   GS_b_ncol,
     &                   ISYMOP(ibc)-1,
     &                   GAO,
     &                   GAO1,
     &                   DTRMAT(1,1,1,ibc),
     &                   BUF)
        ibc = 2
        CALL GET_M_DERVS(.TRUE.,
     &                   GS_c_ncol,
     &                   ISYMOP(ibc)-1,
     &                   GAO,
     &                   GAO1,
     &                   DTRMAT(1,1,1,ibc),
     &                   BUF)
        ibc = 3
        CALL GET_M_DERVS(.TRUE.,
     &                   GS_bc_ncol,
     &                   ISYMOP(ibc)-1,
     &                   GAO,
     &                   GAO1,
     &                   DTRMAT(1,1,1,ibc),
     &                   BUF)
      ENDIF

      R_pq_bc = D0
      R_pq_b  = D0
      R_pq_c  = D0
      R_pq    = D0
      S_pq_bc = D0
      S_pq_b  = D0
      S_pq_c  = D0
      S_pq    = D0
      Z_pq_bc = D0
      Z_pq_b  = D0
      Z_pq_c  = D0
      Z_pq    = D0
      Y_pq_bc = D0
      Y_pq_b  = D0
      Y_pq_c  = D0
      Y_pq    = D0
      X_pq_bc = D0
      X_pq_b  = D0
      X_pq_c  = D0
      X_pq    = D0

      IF(ADD_R .AND.       ADD_GR) CALL DFTPOT2(VXC,W,R,GR_NORM)
      IF(ADD_R .AND. .NOT. ADD_GR) CALL DFTPOT2(VXC,W,R,D0)

      Z_b  =    DDOT(3,GR_b ,1,GR   ,1)
     &     +    DDOT(3,GR   ,1,GR_b ,1)
      Z_c  =    DDOT(3,GR_c ,1,GR   ,1)
     &     +    DDOT(3,GR   ,1,GR_c ,1)
      Z_bc =    DDOT(3,GR_bc,1,GR   ,1)
     &     +    DDOT(3,GR   ,1,GR_bc,1)
     &     + D2*DDOT(3,GR_b ,1,GR_c ,1)
     &     + D2*DDOT(3,GR_c ,1,GR_b ,1)

      IF(ADD_R) THEN
        R_pq_bc = R_pq_bc + VXC(fR0)

        R_pq_b  = R_pq_b
     &          + R_c *VXC(fRR)

        R_pq_c  = R_pq_c
     &          + R_b *VXC(fRR)

        R_pq    = R_pq
     &          + R_bc*VXC(fRR)

        R_pq    = R_pq
     &          + (R_b*R_c + R_c*R_b)*VXC(fRRR)
      ENDIF

      IF(ADD_GR) THEN
        R_pq_b  = R_pq_b
     &          + Z_c *VXC(fRZ)

        R_pq_c  = R_pq_c
     &          + Z_b *VXC(fRZ)

        R_pq    = R_pq
     &          + Z_bc*VXC(fRZ)

        R_pq    = R_pq
     &          + (R_b*Z_c + R_c*Z_b)*VXC(fRRZ)
     &          + (Z_b*R_c + Z_c*R_b)*VXC(fRRZ)
     &          + (Z_b*Z_c + Z_c*Z_b)*VXC(fRZZ)

        Z_pq_bc = Z_pq_bc + VXC(fZ0)

        Z_pq_b  = Z_pq_b
     &          + R_c *VXC(fRZ)
     &          + Z_c *VXC(fZZ)

        Z_pq_c  = Z_pq_c
     &          + R_b *VXC(fRZ)
     &          + Z_b *VXC(fZZ)

        Z_pq    = Z_pq
     &          + R_bc*VXC(fRZ)
     &          + Z_bc*VXC(fZZ)

        Z_pq    = Z_pq
     &          + (R_b*R_c + R_c*R_b)*VXC(fRRZ)
     &          + (R_b*Z_c + R_c*Z_b)*VXC(fRZZ)
     &          + (Z_b*R_c + Z_c*R_b)*VXC(fRZZ)
     &          + (Z_b*Z_c + Z_c*Z_b)*VXC(fZZZ)
      ENDIF

      IF(ADD_S) THEN

        IF(ADD_S .AND.       ADD_GS) CALL SDFTPOT2(VXC,W,R,GR_NORM)
        IF(ADD_S .AND. .NOT. ADD_GS) CALL SDFTPOT2(VXC,W,R,D0)

        ISTART = 1
        IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

        DO ixyz = ISTART,3

          S_b  = S_b_ncol(ixyz)
          S_c  = S_c_ncol(ixyz)
          S_bc = S_bc_ncol(ixyz)

          Y_b  =    DDOT(3,GR_b              ,1,GS                ,1)
     &         +    DDOT(3,GR                ,1,GS_b_ncol(1,ixyz) ,1)
          Y_c  =    DDOT(3,GR_c              ,1,GS                ,1)
     &         +    DDOT(3,GR                ,1,GS_c_ncol(1,ixyz) ,1)
          Y_bc =    DDOT(3,GR_bc             ,1,GS                ,1)
     &         +    DDOT(3,GR                ,1,GS_bc_ncol(1,ixyz),1)
     &         + D2*DDOT(3,GR_b              ,1,GS_c_ncol(1,ixyz) ,1)
     &         + D2*DDOT(3,GR_c              ,1,GS_b_ncol(1,ixyz) ,1)

          X_b  =    DDOT(3,GS_b_ncol(1,ixyz) ,1,GS                ,1)
     &         +    DDOT(3,GS                ,1,GS_b_ncol(1,ixyz) ,1)
          X_c  =    DDOT(3,GS_c_ncol(1,ixyz) ,1,GS                ,1)
     &         +    DDOT(3,GS                ,1,GS_c_ncol(1,ixyz) ,1)
          X_bc =    DDOT(3,GS_bc_ncol(1,ixyz),1,GS                ,1)
     &         +    DDOT(3,GS                ,1,GS_bc_ncol(1,ixyz),1)
     &         + D2*DDOT(3,GS_b_ncol(1,ixyz) ,1,GS_c_ncol(1,ixyz) ,1)
     &         + D2*DDOT(3,GS_c_ncol(1,ixyz) ,1,GS_b_ncol(1,ixyz) ,1)

          R_pq    = R_pq
     &            + (S_b*S_c + S_c*S_b)*VXC(fRSS)

          S_pq_b  = S_pq_b
     &            + S_c *VXC(fSS)

          S_pq_c  = S_pq_c
     &            + S_b *VXC(fSS)

          S_pq    = S_pq
     &            + S_bc*VXC(fSS)

          S_pq    = S_pq
     &            + (R_b*S_c + R_c*S_b)*VXC(fRSS)
     &            + (S_b*R_c + S_c*R_b)*VXC(fRSS)

          IF(ADD_GR) THEN
            S_pq    = S_pq
     &              + (S_b*Z_c + S_c*Z_b)*VXC(fSSZ)
     &              + (Z_b*S_c + Z_c*S_b)*VXC(fSSZ)

            Z_pq    = Z_pq
     &              + (S_b*S_c + S_c*S_b)*VXC(fSSZ)
          ENDIF

          IF(ADD_GS) THEN
            R_pq_b  = R_pq_b
     &              + X_c *VXC(fRX)

            R_pq_c  = R_pq_c
     &              + X_b *VXC(fRX)

            R_pq    = R_pq
     &              + X_bc*VXC(fRX)

            R_pq    = R_pq
     &              + (R_b*X_c + R_c*X_b)*VXC(fRRX)
     &              + (S_b*Y_c + S_c*Y_b)*VXC(fRSY)
     &              + (Z_b*X_c + Z_c*X_b)*VXC(fRZX)
     &              + (Y_b*S_c + Y_c*S_b)*VXC(fRSY)
     &              + (Y_b*Y_c + Y_c*Y_b)*VXC(fRYY)
     &              + (X_b*R_c + X_c*R_b)*VXC(fRRX)
     &              + (X_b*Z_c + X_c*Z_b)*VXC(fRZX)
     &              + (X_b*X_c + X_c*X_b)*VXC(fRXX)

            S_pq_b  = S_pq_b
     &              + Y_c *VXC(fSY)

            S_pq_c  = S_pq_c
     &              + Y_b *VXC(fSY)

            S_pq    = S_pq
     &              + Y_bc*VXC(fSY)

            S_pq    = S_pq
     &              + (R_b*Y_c + R_c*Y_b)*VXC(fRSY)
     &              + (S_b*X_c + S_c*X_b)*VXC(fSSX)
     &              + (Z_b*Y_c + Z_c*Y_b)*VXC(fSZY)
     &              + (Y_b*R_c + Y_c*R_b)*VXC(fRSY)
     &              + (Y_b*Z_c + Y_c*Z_b)*VXC(fSZY)
     &              + (Y_b*X_c + Y_c*X_b)*VXC(fSYX)
     &              + (X_b*S_c + X_c*S_b)*VXC(fSSX)
     &              + (X_b*Y_c + X_c*Y_b)*VXC(fSYX)

            Z_pq_b  = Z_pq_b
     &              + X_c *VXC(fZX)

            Z_pq_c  = Z_pq_c
     &              + X_b *VXC(fZX)

            Z_pq    = Z_pq
     &              + X_bc*VXC(fZX)

            Z_pq    = Z_pq
     &              + (R_b*X_c + R_c*X_b)*VXC(fRZX)
     &              + (S_b*Y_c + S_c*Y_b)*VXC(fSZY)
     &              + (Z_b*X_c + Z_c*X_b)*VXC(fZZX)
     &              + (Y_b*S_c + Y_c*S_b)*VXC(fSZY)
     &              + (Y_b*Y_c + Y_c*Y_b)*VXC(fZYY)
     &              + (X_b*R_c + X_c*R_b)*VXC(fRZX)
     &              + (X_b*Z_c + X_c*Z_b)*VXC(fZZX)
     &              + (X_b*X_c + X_c*X_b)*VXC(fZXX)

            Y_pq_b  = Y_pq_b
     &              + S_c *VXC(fSY)
     &              + Y_c *VXC(fYY)

            Y_pq_c  = Y_pq_c
     &              + S_b *VXC(fSY)
     &              + Y_b *VXC(fYY)

            Y_pq    = Y_pq
     &              + S_bc*VXC(fSY)
     &              + Y_bc*VXC(fYY)

            Y_pq    = Y_pq
     &              + (R_b*S_c + R_c*S_b)*VXC(fRSY)
     &              + (R_b*Y_c + R_c*Y_b)*VXC(fRYY)
     &              + (S_b*R_c + S_c*R_b)*VXC(fRSY)
     &              + (S_b*Z_c + S_c*Z_b)*VXC(fSZY)
     &              + (S_b*X_c + S_c*X_b)*VXC(fSYX)
     &              + (Z_b*S_c + Z_c*S_b)*VXC(fSZY)
     &              + (Z_b*Y_c + Z_c*Y_b)*VXC(fZYY)
     &              + (Y_b*R_c + Y_c*R_b)*VXC(fRYY)
     &              + (Y_b*Z_c + Y_c*Z_b)*VXC(fZYY)
     &              + (Y_b*X_c + Y_c*X_b)*VXC(fYYX)
     &              + (X_b*S_c + X_c*S_b)*VXC(fSYX)
     &              + (X_b*Y_c + X_c*Y_b)*VXC(fYYX)

            X_pq_bc = X_pq_bc + VXC(fX)

            X_pq_b  = X_pq_b
     &              + R_c *VXC(fRX)
     &              + Z_c *VXC(fZX)
     &              + X_c *VXC(fXX)

            X_pq_c  = X_pq_c
     &              + R_b *VXC(fRX)
     &              + Z_b *VXC(fZX)
     &              + X_b *VXC(fXX)

            X_pq    = X_pq
     &              + R_bc*VXC(fRX)
     &              + Z_bc*VXC(fZX)
     &              + X_bc*VXC(fXX)

            X_pq    = X_pq
     &              + (R_b*R_c + R_c*R_b)*VXC(fRRX)
     &              + (R_b*Z_c + R_c*Z_b)*VXC(fRZX)
     &              + (R_b*X_c + R_c*X_b)*VXC(fRXX)
     &              + (S_b*S_c + S_c*S_b)*VXC(fSSX)
     &              + (S_b*Y_c + S_c*Y_b)*VXC(fSYX)
     &              + (Z_b*R_c + Z_c*R_b)*VXC(fRZX)
     &              + (Z_b*Z_c + Z_c*Z_b)*VXC(fZZX)
     &              + (Z_b*X_c + Z_c*X_b)*VXC(fZXX)
     &              + (Y_b*S_c + Y_c*S_b)*VXC(fSYX)
     &              + (Y_b*Y_c + Y_c*Y_b)*VXC(fYYX)
     &              + (X_b*R_c + X_c*R_b)*VXC(fRXX)
     &              + (X_b*Z_c + X_c*Z_b)*VXC(fZXX)
     &              + (X_b*X_c + X_c*X_b)*VXC(fXXX)
          ENDIF

























        ENDDO
      ENDIF


































      ibc = 2
      VXC_SCALAR = R_pq_b
      DO jxyz = 1,3
        VXC_VECTOR(jxyz) = D2*Z_pq_b *GR(jxyz)
     &                   + D2*Z_pq_bc*GR_c(jxyz)
      ENDDO
!     call doexc1(VXC_SCALAR,ISYMOP(ibc)-1,EXCMAT(1,1,1,ibc),GAO)
      CALL DOEXC2(VXC_SCALAR,
     &            VXC_VECTOR,
     &            ISYMOP(ibc)-1,
     &            EXCMAT(1,1,1,ibc),
     &            GAO,
     &            GAO1)

      ibc = 1
      VXC_SCALAR = R_pq_c
      DO jxyz = 1,3
        VXC_VECTOR(jxyz) = D2*Z_pq_c *GR(jxyz)
     &                   + D2*Z_pq_bc*GR_b(jxyz)
      ENDDO
!     call doexc1(VXC_SCALAR,ISYMOP(ibc)-1,EXCMAT(1,1,1,ibc),GAO)
      CALL DOEXC2(VXC_SCALAR,
     &            VXC_VECTOR,
     &            ISYMOP(ibc)-1,
     &            EXCMAT(1,1,1,ibc),
     &            GAO,
     &            GAO1)

      ibc = 3
      VXC_SCALAR = R_pq
      DO jxyz = 1,3
        VXC_VECTOR(jxyz) = D2*Z_pq   *GR(jxyz)
     &                   + D2*Z_pq_b *GR_b(jxyz)*D2
     &                   + D2*Z_pq_c *GR_c(jxyz)*D2
     &                   + D2*Z_pq_bc*GR_bc(jxyz)
      ENDDO
!     call doexc1(VXC_SCALAR,ISYMOP(ibc)-1,EXCMAT(1,1,1,ibc),GAO)
      CALL DOEXC2(VXC_SCALAR,
     &            VXC_VECTOR,
     &            ISYMOP(ibc)-1,
     &            EXCMAT(1,1,1,ibc),
     &            GAO,
     &            GAO1)







      RETURN
      END
