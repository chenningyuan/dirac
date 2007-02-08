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
     &                   DENSITY_GRADIENT,
     &                   DENSITY_GRADIENT_NORM,
     &                   WEIGHT,DOGGA,BUF)
C=======================================================================
C     contribution of magnetization to the sigma vector of QR
C     adaptation of Johan's dft_quadres, so thanks for his efforts!
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: february 2007
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
#include "dftcom.h"
#include "cdft.h"
#include "../dft/functionals.h"

      LOGICAL   DOGGA, LBIT
      DIMENSION DTRMAT(NBAST,NBAST,NZ,3),
     &          EXCMAT(*),
     &          ISYMOP(3),
     &          IHRMOP(3),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          SPIN_DENSITY(3,3),
     &          SPIN_DENSITY_NORM(3),
     &          VXC(99), !dirty fix
     &          VXC_V(3)

#include "ibtfun.h"

      F  = D1
      G  = D1
      IF(LBIT(LBIT_CDFT,4)) G = -D1 !spin density operator is \beta*\Sigma instead of \Sigma
      IU = 1
      IV = 1
      IX = 2
      IY = 2
      ITRANSFORM = 2

      ICHECK = IHRMOP(1) + IHRMOP(2) + IHRMOP(3)
      IF(ICHECK .GT. -3) CALL QUIT('problem in SDFT QR')

      IF(DOGGA) CALL QUIT('SDFT QR implemented only for AL(S)DA')

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      DO IBC = 1,3
        CALL GET_CURRENT(2,
     &                   SPIN_DENSITY(1,IBC),
     &                   DUMMY,
     &                   ISYMOP(IBC)-1,GAO,DUMMY,DTRMAT(1,1,1,IBC),BUF)

        SPIN_DENSITY_NORM(IBC) = DSQRT(DDOT(3,
     &                                      SPIN_DENSITY(1,IBC),1,
     &                                      SPIN_DENSITY(1,IBC),1
     &                                     ))
      ENDDO

      TEST_MAX = MAX(
     &               SPIN_DENSITY_NORM(1),
     &               SPIN_DENSITY_NORM(2),
     &               SPIN_DENSITY_NORM(3)
     &              )

      IF(TEST_MAX .GT. DFTHRL) THEN

!       CALL SDFTPOT2(VXC,WEIGHT,DENSITY,D0)

        DO IBC = 1,3
          IF(SPIN_DENSITY_NORM(IBC) .GT. DFTHRL) THEN
            DO IXYZ = ISTART,3
              IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
              IOFF = NBAST*NBAST*(IQ-1)
              IOFF = IOFF + NBAST*NBAST*NZ*(IBC-1)

              IREPOP = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),ISYMOP(IBC)-1)

              VXC_S = VXC(FRR)*SPIN_DENSITY(IXYZ,IBC)

              CALL CDFT_EXC(VXC_S,F,IU,IV,
     &                      GAO,IREPOP,EXCMAT(1+IOFF),BUF)
              CALL CDFT_EXC(VXC_S,G,IX,IY,
     &                      GAO,IREPOP,EXCMAT(1+IOFF),BUF)
            ENDDO
          ENDIF
        ENDDO

      ENDIF

      RETURN
      END
