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
#include "dcbham.h"
#include "dftcom.h"
#include "cdft.h"
#include "../dft/functionals.h"

      LOGICAL   DOGGA, LBIT
      DIMENSION DTRMAT(NBAST,NBAST,NZ,3),
     &          EXCMAT(NBAST,NBAST,NZ,3),
     &          ISYMOP(3),
     &          IHRMOP(3),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          SPIN_DENSITY_B(3),
     &          SPIN_DENSITY_C(3),
     &          SPIN_DENSITY_BC(3),
     &          VXC(99), !dirty fix
     &          VXC_S(3,3),
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

      IH_B  = IHRMOP(1)
      IH_C  = IHRMOP(2)
      IH_BC = IHRMOP(3)

      ICHECK = IH_B + IH_C + IH_BC

      IF(ICHECK .EQ. 3) GOTO 999

      IF(DOGGA) CALL QUIT('SDFT QR implemented only for AL(S)DA')

      IREP_B  = ISYMOP(1) - 1
      IREP_C  = ISYMOP(2) - 1
      IREP_BC = ISYMOP(3) - 1

      SPIN_DENSITY_NORM_B  = D0
      SPIN_DENSITY_NORM_C  = D0
      SPIN_DENSITY_NORM_BC = D0

      DENSITY_B  = D0
      DENSITY_C  = D0

      IF(.NOT. IH_B .EQ. 1) THEN
!       CALL GET_CURRENT(2,
!    &                   SPIN_DENSITY_B,
!    &                   DUMMY,
!    &                   IREP_B,GAO,DUMMY,DTRMAT(1,1,1,1),BUF)

        SPIN_DENSITY_NORM_B = DSQRT(DDOT(3,
     &                                   SPIN_DENSITY_B,1,
     &                                   SPIN_DENSITY_B,1
     &                                  ))
      ENDIF

      IF(.NOT. IH_B .EQ. -1) THEN
        CALL GETRHO(DENSITY_B,IREP_B,GAO,DTRMAT(1,1,1,1),BUF)
      ENDIF

      IF(.NOT. IH_C .EQ. 1) THEN
!       CALL GET_CURRENT(2,
!    &                   SPIN_DENSITY_C,
!    &                   DUMMY,
!    &                   IREP_C,GAO,DUMMY,DTRMAT(1,1,1,2),BUF)

        SPIN_DENSITY_NORM_C = DSQRT(DDOT(3,
     &                                   SPIN_DENSITY_C,1,
     &                                   SPIN_DENSITY_C,1
     &                                  ))
      ENDIF

      IF(.NOT. IH_C .EQ. -1) THEN
        CALL GETRHO(DENSITY_C,IREP_C,GAO,DTRMAT(1,1,1,2),BUF)
      ENDIF

      IF(.NOT. IH_BC .EQ. 1) THEN
!       CALL GET_CURRENT(2,
!    &                   SPIN_DENSITY_BC,
!    &                   DUMMY,
!    &                   IREP_BC,GAO,DUMMY,DTRMAT(1,1,1,3),BUF)

        SPIN_DENSITY_NORM_BC = DSQRT(DDOT(3,
     &                                    SPIN_DENSITY_BC,1,
     &                                    SPIN_DENSITY_BC,1
     &                                   ))
      ENDIF

      CALL SDFTPOT2(VXC,WEIGHT,DENSITY,D0)

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      IF(SPIN_DENSITY_NORM_B .GT. DFTHRL) THEN
        IBC = 2
        DO IXYZ = ISTART,3
          VXC_S = VXC(FRR)*SPIN_DENSITY_B(IXYZ)

          IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
          IREPOP = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),ISYMOP(IBC)-1)

!         CALL CDFT_EXC(VXC_S,F,IU,IV,
!    &                  GAO,IREPOP,EXCMAT(1,1,IQ,IBC),BUF)
!         IF(MC .EQ. 2) THEN
!           CALL CDFT_EXC(VXC_S,G,IX,IY,
!    &                    GAO,IREPOP,EXCMAT(1,1,IQ,IBC),BUF)
!         ENDIF
        ENDDO
      ENDIF

      IF(SPIN_DENSITY_NORM_C .GT. DFTHRL) THEN
        IBC = 1
        DO IXYZ = ISTART,3
          VXC_S = VXC(FRR)*SPIN_DENSITY_C(IXYZ)

          IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
          IREPOP = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),ISYMOP(IBC)-1)

!         CALL CDFT_EXC(VXC_S,F,IU,IV,
!    &                  GAO,IREPOP,EXCMAT(1,1,IQ,IBC),BUF)
!         IF(MC .EQ. 2) THEN
!           CALL CDFT_EXC(VXC_S,G,IX,IY,
!    &                    GAO,IREPOP,EXCMAT(1,1,IQ,IBC),BUF)
!         ENDIF
        ENDDO
      ENDIF

      TEST_MAX = MAX(
     &               SPIN_DENSITY_NORM_B,
     &               SPIN_DENSITY_NORM_C,
     &               SPIN_DENSITY_NORM_BC
     &              )

      IF(TEST_MAX .GT. DFTHRL) THEN
        IBC = 3
        DO IXYZ = ISTART,3

          PART1 = SPIN_DENSITY_BC(IXYZ)
          PART2 = DENSITY_B*SPIN_DENSITY_C(IXYZ)
          PART3 = DENSITY_C*SPIN_DENSITY_B(IXYZ)

          VXC_S = VXC(FRR)*PART1 + VXC(FRRR)*(PART2 + PART3)

          IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
          IREPOP = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),ISYMOP(IBC)-1)

!         CALL CDFT_EXC(VXC_S,F,IU,IV,
!    &                  GAO,IREPOP,EXCMAT(1,1,IQ,IBC),BUF)
!         IF(MC .EQ. 2) THEN
!           CALL CDFT_EXC(VXC_S,G,IX,IY,
!    &                    GAO,IREPOP,EXCMAT(1,1,IQ,IBC),BUF)
!         ENDIF
        ENDDO
      ENDIF

!     rrr ub cdot uc rho_pq missing ...

 999  CONTINUE

      RETURN
      END
