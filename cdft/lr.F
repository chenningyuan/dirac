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
C     /* SDFT_LR */
      SUBROUTINE SDFT_LR(DTRMAT,
     &                   EXCMAT,
     &                   NXCMAT,
     &                   ISYMOP,
     &                   IHRMOP,
     &                   GAO,
     &                   GAO1,
     &                   R,
     &                   GR_NORM,
     &                   GR,
     &                   W,
     &                   IS_GGA,
     &                   BUF)
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

      LOGICAL   IS_GGA,
     &          LBIT
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          ISYMOP(*),
     &          IHRMOP(*),
     &          GAO(*),
     &          GAO1(*),
     &          GR(3),
     &          BUF(*),
     &          S_b(3),
     &          GS_b(3,3),
     &          VXC(99) !dirty fix

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      DO IMAT = 1,NXCMAT !loop over all frequencies
        IF(IHRMOP(IMAT) .EQ. -1) THEN

          JOFF = NBAST*NBAST*NZ*(IMAT-1)

          CALL GET_M(.TRUE.,
     &               S_b,
     &               ISYMOP(IMAT)-1,
     &               GAO,
     &               UNUSED,
     &               DTRMAT(1+JOFF),
     &               BUF)

          IF(IS_GGA) THEN
            CALL GET_M_DERVS(.TRUE.,
     &                       GS_b,
     &                       ISYMOP(IMAT)-1,
     &                       GAO,
     &                       GAO1,
     &                       DTRMAT(1+JOFF),
     &                       BUF)
          ENDIF

          TEST_ME_FIRST = DSQRT(DDOT(3,S_b,1,S_b,1))

          IF(IS_GGA) THEN
            TEST_ME_FIRST = MAX(
     &                          TEST_ME_FIRST,
     &                          DABS(GS_b(1,1)),
     &                          DABS(GS_b(1,2)),
     &                          DABS(GS_b(1,3)),
     &                          DABS(GS_b(2,1)),
     &                          DABS(GS_b(2,2)),
     &                          DABS(GS_b(2,3)),
     &                          DABS(GS_b(3,1)),
     &                          DABS(GS_b(3,2)),
     &                          DABS(GS_b(3,3))
     &                         )
          ENDIF

          IF(TEST_ME_FIRST .GT. DFTHRL) THEN

            IF(IS_GGA) THEN
              CALL SDFTPOT1(VXC,W,R,GR_NORM)
            ELSE
              CALL SDFTPOT1(VXC,W,R,D0)
            ENDIF

            DO ixyz = ISTART,3
              CALL SDFT_LR_COLLINEAR(VXC,
     &                               EXCMAT(1+JOFF),
     &                               ISYMOP(IMAT),
     &                               GAO,
     &                               GAO1,
     &                               ixyz,
     &                               S_b(ixyz),
     &                               GS_b(1,ixyz),
     &                               GR,
     &                               IS_GGA)
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SDFT_LR_COLLINEAR */
      SUBROUTINE SDFT_LR_COLLINEAR(VXC,
     &                             EXCMAT,
     &                             ISYMOP,
     &                             GAO,
     &                             GAO1,
     &                             ixyz,
     &                             S_b,
     &                             GS_b,
     &                             GR,
     &                             IS_GGA)
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

      LOGICAL   IS_GGA
      DIMENSION VXC(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*),
     &          GS_b(3),
     &          GR(3),
     &          Y_pq(3),
     &          X_pq_b(3),
     &          VXC_VECTOR(3)


      S_pq = VXC(fSS)*S_b

      IF(IS_GGA) THEN
        Y_b = DDOT(3,GR,1,GS_b,1)

        DO jxyz = 1,3
          X_pq_b(jxyz) = VXC(fX)

          Y_pq(jxyz) = S_b*VXC(fSY) + Y_b*VXC(fYY)
        ENDDO

        S_pq = S_pq + VXC(fSY)*Y_b
      ENDIF

      VXC_SCALAR = S_pq

      IF(IS_GGA) THEN
        DO jxyz = 1,3
          VXC_VECTOR(jxyz) = Y_pq(jxyz)*GR(jxyz)
     &                     + X_pq_b(jxyz)*GS_b(jxyz)*D2
        ENDDO

        CALL DSCAL(3,D2,VXC_VECTOR,1) !symmetrization later in the code
      ENDIF

      CALL EXC_M(IS_GGA,
     &           VXC_SCALAR,
     &           VXC_VECTOR,
     &           ixyz,
     &           ISYMOP,
     &           EXCMAT,
     &           GAO,
     &           GAO1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CDFT_LR */
      SUBROUTINE CDFT_LR(DTRMAT,
     &                   EXCMAT,
     &                   NXCMAT,
     &                   ISYMOP,
     &                   IHRMOP,
     &                   GAO,
     &                   GAO1,
     &                   DENSITY,
     &                   DENSITY_GRADIENT_NORM,
     &                   DENSITY_GRADIENT,
     &                   WEIGHT,
     &                   IS_GGA,
     &                   BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (DP13 = D1/D3)

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

      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          ISYMOP(*),
     &          IHRMOP(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          FO_VORTICITY(3)

      DO IMAT = 1,NXCMAT !loop over all frequencies
        IF(IHRMOP(IMAT) .EQ. -1) THEN

          JOFF = NBAST*NBAST*NZ*(IMAT-1)

          CALL GET_ROTVPARA(.TRUE.,
     &                      FO_VORTICITY,
     &                      DENSITY,
     &                      DENSITY_GRADIENT,
     &                      ISYMOP(IMAT)-1,
     &                      GAO,
     &                      GAO1,
     &                      DTRMAT(1+JOFF),
     &                      BUF)

          RS  = (4.1887902D0*DENSITY)**(-DP13)

          IF(RS .GT. 10.0D0) GOTO 999

          DKF = (29.6088132D0*DENSITY)**DP13

          CHI_RATIO = (D1 + 0.028D0*RS)*DEXP(-0.042D0*RS)

          G = D2*WEIGHT*(DKF/236.870506D0)*(CHI_RATIO - D1)

          CALL DSCAL(3,G,FO_VORTICITY,1)

          FO_VORTICITY_NORM = DSQRT(DDOT(3,FO_VORTICITY,1,
     &                                     FO_VORTICITY,1))

          IF(FO_VORTICITY_NORM .LT. DFTHRL) GOTO 999

          CALL GET_ROTVPARA(.FALSE.,
     &                      FO_VORTICITY,
     &                      DENSITY,
     &                      DENSITY_GRADIENT,
     &                      ISYMOP(IMAT)-1,
     &                      GAO,
     &                      GAO1,
     &                      EXCMAT(1+JOFF),
     &                      BUF)

 999      CONTINUE

        ENDIF
      ENDDO

      RETURN
      END
