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
     &                               R_0,
     &                               GR_0_NORM,
     &                               GR_0,
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
     &          GR_0(3),
     &          BUF(*),
!
     &          GR(  3,0:3),
     &          GS(3,3,0:3),
!
     &          R(  0:3),
     &          Z(  0:3),
     &          S(3,0:3),
     &          Y(3,0:3),
     &          X(3,0:3),
!
     &          R_pq(  0:3),
     &          Z_pq(  0:3),
     &          S_pq(3,0:3),
     &          Y_pq(3,0:3),
     &          X_pq(3,0:3),
!
     &          Vxc_hs(    0:3),
     &          Vxc_as(  3,0:3),
     &          Vxc_hv(  3,0:3),
     &          Vxc_av(3,3,0:3),
!
     &          Vxc(99) !dirty fix

#include "ibtfun.h"

      CALL DZERO(GR,4*3)
      CALL DZERO(GS,4*3*3)

      CALL DZERO(R,4)
      CALL DZERO(Z,4)
      CALL DZERO(S,4*3)
      CALL DZERO(Y,4*3)
      CALL DZERO(X,4*3)

      CALL DZERO(R_pq,4)
      CALL DZERO(Z_pq,4)
      CALL DZERO(S_pq,4*3)
      CALL DZERO(Y_pq,4*3)
      CALL DZERO(X_pq,4*3)

      CALL DZERO(Vxc_hs,4)
      CALL DZERO(Vxc_as,4*3)
      CALL DZERO(Vxc_hv,4*3)
      CALL DZERO(Vxc_av,4*3*3)

      GR(1,0) = GR_0(1)
      GR(2,0) = GR_0(2)
      GR(3,0) = GR_0(3)

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      DO ibc = 1,3
        IF(ADD_R .OR. ADD_S) THEN
          CALL GETRHO(R(ibc),
     &                ISYMOP(ibc)-1,
     &                GAO,
     &                DTRMAT(1,1,1,ibc),
     &                BUF)
        ENDIF
        IF(ADD_GR) THEN
          CALL DENSITY_GRADIENT_GENERAL(GR(1,ibc),
     &                                  ISYMOP(ibc)-1,
     &                                  GAO,
     &                                  GAO1,
     &                                  DTRMAT(1,1,1,ibc),
     &                                  BUF)
        ENDIF
        IF(ADD_S) THEN
          CALL GET_M(.TRUE.,
     &               S(1,ibc),
     &               ISYMOP(ibc)-1,
     &               GAO,
     &               UNUSED,
     &               DTRMAT(1,1,1,ibc),
     &               BUF)
        ENDIF
        IF(ADD_GS) THEN
          CALL GET_M_DERVS(.TRUE.,
     &                     GS(1,1,ibc),
     &                     ISYMOP(ibc)-1,
     &                     GAO,
     &                     GAO1,
     &                     DTRMAT(1,1,1,ibc),
     &                     BUF)
        ENDIF
      ENDDO

      IF(ADD_GR) THEN
        Z(     1) = DDOT(3,GR(1,     1),1,GR(1,     0),1)
     &            + DDOT(3,GR(1,     0),1,GR(1,     1),1)
        Z(     2) = DDOT(3,GR(1,     2),1,GR(1,     0),1)
     &            + DDOT(3,GR(1,     0),1,GR(1,     2),1)
        Z(     3) = DDOT(3,GR(1,     3),1,GR(1,     0),1)
     &            + DDOT(3,GR(1,     0),1,GR(1,     3),1)
     &            + DDOT(3,GR(1,     1),1,GR(1,     2),1)*D2
     &            + DDOT(3,GR(1,     2),1,GR(1,     1),1)*D2
      ENDIF
      IF(ADD_GR) THEN
        DO ixyz = ISTART,3
          Y(ixyz,1) = DDOT(3,GR(1,     1),1,GS(1,ixyz,0),1)
     &              + DDOT(3,GR(1,     0),1,GS(1,ixyz,1),1)
          Y(ixyz,2) = DDOT(3,GR(1,     2),1,GS(1,ixyz,0),1)
     &              + DDOT(3,GR(1,     0),1,GS(1,ixyz,2),1)
          Y(ixyz,3) = DDOT(3,GR(1,     3),1,GS(1,ixyz,0),1)
     &              + DDOT(3,GR(1,     0),1,GS(1,ixyz,3),1)
     &              + DDOT(3,GR(1,     1),1,GS(1,ixyz,2),1)*D2
     &              + DDOT(3,GR(1,     2),1,GS(1,ixyz,1),1)*D2

          X(ixyz,1) = DDOT(3,GS(1,ixyz,1),1,GS(1,ixyz,0),1)
     &              + DDOT(3,GS(1,ixyz,0),1,GS(1,ixyz,1),1)
          X(ixyz,2) = DDOT(3,GS(1,ixyz,2),1,GS(1,ixyz,0),1)
     &              + DDOT(3,GS(1,ixyz,0),1,GS(1,ixyz,2),1)
          X(ixyz,3) = DDOT(3,GS(1,ixyz,3),1,GS(1,ixyz,0),1)
     &              + DDOT(3,GS(1,ixyz,0),1,GS(1,ixyz,3),1)
     &              + DDOT(3,GS(1,ixyz,1),1,GS(1,ixyz,2),1)*D2
     &              + DDOT(3,GS(1,ixyz,2),1,GS(1,ixyz,1),1)*D2
        ENDDO
      ENDIF

      IF(ADD_R .AND.       ADD_GR) CALL DFTPOT2(VXC,W,R_0,GR_0_NORM)
      IF(ADD_R .AND. .NOT. ADD_GR) CALL DFTPOT2(VXC,W,R_0,D0)

      IF(ADD_R) THEN
        R_pq(3) = R_pq(3) + Vxc(fR0)

        R_pq(1) = R_pq(1)
     &          + R(2)*Vxc(fRR)

        R_pq(2) = R_pq(2)
     &          + R(1)*Vxc(fRR)

        R_pq(0) = R_pq(0)
     &          + R(3)*Vxc(fRR)

        R_pq(0) = R_pq(0)
     &          + (R(1)*R(2) + R(2)*R(1))*Vxc(fRRR)
      ENDIF

      IF(ADD_GR) THEN
        R_pq(1) = R_pq(1)
     &          + Z(2)*Vxc(fRZ)

        R_pq(2) = R_pq(2)
     &          + Z(1)*Vxc(fRZ)

        R_pq(0) = R_pq(0)
     &          + Z(3)*Vxc(fRZ)

        R_pq(0) = R_pq(0)
     &          + (R(1)*Z(2) + R(2)*Z(1))*Vxc(fRRZ)
     &          + (Z(1)*R(2) + Z(2)*R(1))*Vxc(fRRZ)
     &          + (Z(1)*Z(2) + Z(2)*Z(1))*Vxc(fRZZ)

        Z_pq(3) = Z_pq(3) + Vxc(fZ0)

        Z_pq(1) = Z_pq(1)
     &          + R(2)*Vxc(fRZ)
     &          + Z(2)*Vxc(fZZ)

        Z_pq(2) = Z_pq(2)
     &          + R(1)*Vxc(fRZ)
     &          + Z(1)*Vxc(fZZ)

        Z_pq(0) = Z_pq(0)
     &          + R(3)*Vxc(fRZ)
     &          + Z(3)*Vxc(fZZ)

        Z_pq(0) = Z_pq(0)
     &          + (R(1)*R(2) + R(2)*R(1))*Vxc(fRRZ)
     &          + (R(1)*Z(2) + R(2)*Z(1))*Vxc(fRZZ)
     &          + (Z(1)*R(2) + Z(2)*R(1))*Vxc(fRZZ)
     &          + (Z(1)*Z(2) + Z(2)*Z(1))*Vxc(fZZZ)
      ENDIF

      IF(ADD_S .AND.       ADD_GS) CALL SDFTPOT2(VXC,W,R_0,GR_0_NORM)
      IF(ADD_S .AND. .NOT. ADD_GS) CALL SDFTPOT2(VXC,W,R_0,D0)

      IF(ADD_S) THEN
        DO ixyz = ISTART,3
          R_pq(0) = R_pq(0)
     &            + (S(ixyz,1)*S(ixyz,2)
     &            +  S(ixyz,2)*S(ixyz,1))*Vxc(fRSS)

          S_pq(ixyz,1) = S_pq(ixyz,1)
     &                 + S(ixyz,2)*Vxc(fSS)

          S_pq(ixyz,2) = S_pq(ixyz,2)
     &                 + S(ixyz,1)*Vxc(fSS)

          S_pq(ixyz,0) = S_pq(ixyz,0)
     &                 + S(ixyz,3)*Vxc(fSS)

          S_pq(ixyz,0) = S_pq(ixyz,0)
     &                 + (R(1)*S(ixyz,2) + R(2)*S(ixyz,1))*Vxc(fRSS)
     &                 + (S(ixyz,1)*R(2) + S(ixyz,2)*R(1))*Vxc(fRSS)
        ENDDO
      ENDIF

      Vxc_hs(1) = R_pq(2)
      Vxc_hs(2) = R_pq(1)
      Vxc_hs(3) = R_pq(0)

      DO ixyz = ISTART,3
        Vxc_as(ixyz,1) = S_pq(ixyz,2)
        Vxc_as(ixyz,2) = S_pq(ixyz,1)
        Vxc_as(ixyz,3) = S_pq(ixyz,0)
      ENDDO

      DO jxyz = 1,3
        Vxc_hv(jxyz,1) =  Z_pq(2)*GR(jxyz,0) + Z_pq(3)*GR(jxyz,1)
        Vxc_hv(jxyz,2) =  Z_pq(1)*GR(jxyz,0) + Z_pq(3)*GR(jxyz,2)
        Vxc_hv(jxyz,3) =  Z_pq(0)*GR(jxyz,0) + Z_pq(3)*GR(jxyz,3)
     &                 + (Z_pq(1)*GR(jxyz,1) + Z_pq(2)*GR(jxyz,2))*D2
      ENDDO
      CALL DSCAL(4*3,D2,Vxc_hv(1,0),1)

      DO ibc = 1,3
        WORTH_IT = MAX(
     &                 DABS(Vxc_hs(ibc)),
     &                 DABS(Vxc_hv(1,ibc)),
     &                 DABS(Vxc_hv(2,ibc)),
     &                 DABS(Vxc_hv(3,ibc))
     &                )
        IF(WORTH_IT .GT. DFTHRL) THEN
          CALL DOEXC2(Vxc_hs(ibc),
     &                Vxc_hv(1,ibc),
     &                ISYMOP(ibc)-1,
     &                EXCMAT(1,1,1,ibc),
     &                GAO,
     &                GAO1)
        ENDIF

        WORTH_IT = MAX(
     &                 DABS(Vxc_as(1,ibc)),
     &                 DABS(Vxc_as(2,ibc)),
     &                 DABS(Vxc_as(3,ibc))
     &                )
        IF(WORTH_IT .GT. DFTHRL) THEN
          CALL GET_M(.FALSE.,
     &               Vxc_as(1,ibc),
     &               ISYMOP(ibc)-1,
     &               GAO,
     &               UNUSED,
     &               EXCMAT(1,1,1,ibc),
     &               BUF)
        ENDIF
      ENDDO

      RETURN
      END
