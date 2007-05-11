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
!    &                   NXCMAT,
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
C     contribution of spin density to the sigma vector
C     adaptation of Johan's dft_quadres, so thanks for his efforts!
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: march 2007
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
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          FO_DENSITY(3),
     &          FO_SPIN_DENSITY(3,3),
     &          FO_SPIN_DENSITY_NORM(3),
     &          VXC(99) !dirty fix

#include "ibtfun.h"

      F = D1
      G = D1
      IF(LBIT(LBIT_CDFT,4)) G = -D1 !operator is \beta\sigma instead of \sigma
      IU         = 1
      IV         = 1
      IX         = 2
      IY         = 2
      ITRANSFORM = 2

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      ICHECK = IHRMOP(1) + IHRMOP(2) + IHRMOP(3)
      IF(ICHECK .EQ. 3) GOTO 999

      IF(DOGGA) CALL QUIT('SDFT QR implemented only for AL(S)DA')

      CALL DZERO(FO_DENSITY,3)
      CALL DZERO(FO_SPIN_DENSITY,9)
      CALL DZERO(FO_SPIN_DENSITY_NORM,3)

      DO IBC = 1,3
        IF(.NOT. IHRMOP(IBC) .EQ. 1) THEN
!           CALL GET_CURRENT(.TRUE.,
!    &                       .FALSE.,
!    &                       FO_SPIN_DENSITY(1,IBC),
!    &                       UNUSED,
!    &                       F,
!    &                       G,
!    &                       IU,
!    &                       IV,
!    &                       IX,
!    &                       IY,
!    &                       ITRANSFORM,
!    &                       ISYMOP(IBC)-1,
!    &                       GAO,
!    &                       UNUSED,
!    &                       DTRMAT(1,1,1,IBC),
!    &                       BUF)

          FO_SPIN_DENSITY_NORM(IBC)
     &                         = DSQRT(DDOT(3,FO_SPIN_DENSITY(1,IBC),1,
     &                                        FO_SPIN_DENSITY(1,IBC),1))
        ENDIF
        IF(.NOT. IHRMOP(IBC) .EQ. -1) THEN
          CALL GETRHO(FO_DENSITY(IBC),
     &                ISYMOP(IBC)-1,
     &                GAO,
     &                DTRMAT(1,1,1,IBC),
     &                BUF)
        ENDIF
      ENDDO

      CALL SDFTPOT2(VXC,WEIGHT,DENSITY,D0)

      DO IBC = 1,3
        IF(FO_SPIN_DENSITY_NORM(IBC) .GT. DFTHRL) THEN

          DO IXYZ = ISTART,3
            CALL SDFT_LR_LSDA_COLLINEAR(VXC,
     &                                  EXCMAT(1,1,1,IBC),
     &                                  ISYMOP(IBC)-1,
     &                                  GAO,
     &                                  IXYZ,
     &                                  FO_SPIN_DENSITY(IXYZ,IBC))
          ENDDO
        ENDIF
      ENDDO

      TEST_MAX = MAX(
     &               FO_SPIN_DENSITY_NORM(1),
     &               FO_SPIN_DENSITY_NORM(2),
     &               FO_SPIN_DENSITY_NORM(3)
     &              )

      IF(TEST_MAX .GT. DFTHRL) THEN

        IB  = 1
        IC  = 2
        IBC = 3

        VXC_FRR  = VXC(FRR)
        VXC_FRRR = VXC(FRRR)

        BC_NORM = DSQRT(DDOT(3,FO_SPIN_DENSITY(1,IB),1,
     &                         FO_SPIN_DENSITY(1,IC),1))

        CALL DOEXC1(VXC_FRRR*BC_NORM,
     &              ISYMOP(IBC)-1,
     &              EXCMAT(1,1,1,IBC),
     &              GAO)

        DO IXYZ = ISTART,3

          PART1 = FO_SPIN_DENSITY(IXYZ,IBC)
          PART2 = FO_DENSITY(IB)*FO_SPIN_DENSITY(IXYZ,IC)
          PART3 = FO_DENSITY(IC)*FO_SPIN_DENSITY(IXYZ,IB)

!         cheating to keep the lazy structure
          VXC(FRR) = VXC_FRR*PART1 + VXC_FRRR*(PART2 + PART3)

          CALL SDFT_LR_LSDA_COLLINEAR(VXC,
     &                                EXCMAT(1,1,1,IBC),
     &                                ISYMOP(IBC)-1,
     &                                GAO,
     &                                IXYZ,
     &                                D1)
        ENDDO
      ENDIF

 999  CONTINUE

      RETURN
      END



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* DFT_QR_GGA */
      SUBROUTINE DFT_QR_GGA(DTRMAT,
     &                      EXCMAT,
     &                      ISYMOP,
     &                      IHRMOP,
     &                      GAO,
     &                      GAO1,
     &                      DENSITY,
     &                      DENSITY_GRADIENT_NORM,
     &                      DENSITY_GRADIENT,
     &                      WEIGHT,
     &                      DOGGA,
     &                      BUF)
C=======================================================================
C     does the same thing like Johan's DFT_QR
C     i use it to debug SDFT_QR
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
     &          GAO(*),
     &          GAO1(*),
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

        CALL DENSITY_GRADIENT_GENERAL(FO_SPIN_DENSITY(1,IBC),
     &                                ISYMOP(IBC)-1,
     &                                GAO,
     &                                GAO1,
     &                                DTRMAT(1,1,1,IBC),
     &                                BUF)

        ZETA(IBC) = D2*DDOT(3,FO_SPIN_DENSITY(1,IBC),1,
     &                        DENSITY_GRADIENT,1)
      ENDDO

      ADDME = D4*DDOT(3,FO_SPIN_DENSITY(1,1),1,
     &                  FO_SPIN_DENSITY(1,2),1)

      ZETA(3) = ZETA(3) + ADDME

      CALL DFTPOT2(VXC,WEIGHT,DENSITY,DENSITY_GRADIENT_NORM)

      DO IBC = 1,3
        VXC_S(IBC)  = VXC(FRR)*FO_DENSITY(IBC) + VXC(FRZ)*ZETA(IBC)
        FACTOR(IBC) = VXC(FRZ)*FO_DENSITY(IBC) + VXC(FZZ)*ZETA(IBC)
      ENDDO

      F1 = D2*FO_DENSITY(1)*FO_DENSITY(2)
      F2 = D2*(FO_DENSITY(1)*ZETA(2) + FO_DENSITY(2)*ZETA(1))
      F3 = D2*ZETA(1)*ZETA(2)

      VXC_S(3) = VXC_S(3) + F1*VXC(FRRR)
     &                    + F2*VXC(FRRZ)
     &                    + F3*VXC(FRZZ)

      FACTOR(3) = FACTOR(3) + F1*VXC(FRRZ)
     &                      + F2*VXC(FRZZ)
     &                      + F3*VXC(FZZZ)

      DO IBC = 1,3
        DO IXYZ = 1,3
          VXC_V(IXYZ,IBC) = D2*(FACTOR(IBC)*DENSITY_GRADIENT(IXYZ)
     &                        + VXC(FZ0)*FO_SPIN_DENSITY(IXYZ,IBC))
        ENDDO
      ENDDO

      DO IXYZ = 1,3
        ADDME = FACTOR(1)*FO_SPIN_DENSITY(IXYZ,2)
     &        + FACTOR(2)*FO_SPIN_DENSITY(IXYZ,1)

        VXC_V(IXYZ,3) = VXC_V(IXYZ,3) + D4*ADDME
      ENDDO

      DO IBC = 1,3
        CALL DOEXC2(VXC_S(IBC),
     &              VXC_V(1,IBC),
     &              ISYMOP(IBC)-1,
     &              EXCMAT(1,1,1,IBC),
     &              GAO,
     &              GAO1)
      ENDDO

      RETURN
      END
