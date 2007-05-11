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

      VXC(FZ0) = VXC(FZ0)/D2
      VXC(FRZ) = VXC(FRZ)/D2

      DO IBC = 1,3
        DO IXYZ = 1,3
          VXC_V(IXYZ,IBC) = D2*VXC(FZ0)*FO_SPIN_DENSITY(IXYZ,IBC)
        ENDDO
      ENDDO

      DO IXYZ = 1,3
        ADDME = VXC(FRZ)*FO_DENSITY(1)*FO_SPIN_DENSITY(IXYZ,2)
     &        + VXC(FRZ)*FO_DENSITY(2)*FO_SPIN_DENSITY(IXYZ,1)

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
C     i use(d) it to reverse-engineer and debug SDFT_QR
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

        CALL DENSITY_GRADIENT_GENERAL(FO_SPIN_DENSITY(1,IBC),
     &                                ISYMOP(IBC)-1,
     &                                GAO,
     &                                GAO1,
     &                                DTRMAT(1,1,1,IBC),
     &                                BUF)

        ZETA(IBC) = D2*DDOT(3,FO_SPIN_DENSITY(1,IBC),1,
     &                        DENSITY_GRADIENT,1)
      ENDDO

      ADDME = DDOT(3,FO_SPIN_DENSITY(1,1),1,
     &               FO_SPIN_DENSITY(1,2),1)

      ZETA(3) = ZETA(3) + D4*ADDME

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

!       CALL EXC_GENERAL(.FALSE.,VXC_S(IBC),UNUSED,D1,1,1,
!    &                   ISYMOP(IBC)-1,EXCMAT(1,1,1,IBC),
!    &                   GAO,GAO)
!       DO IXYZ = 1,3
!         CALL EXC_GENERAL(.FALSE.,VXC_V(IXYZ,IBC),UNUSED,D1,1,1,
!    &                     ISYMOP(IBC)-1,EXCMAT(1,1,1,IBC),
!    &                     GAO,GAO1(1,IXYZ))
!         CALL EXC_GENERAL(.FALSE.,VXC_V(IXYZ,IBC),UNUSED,D1,1,1,
!    &                     ISYMOP(IBC)-1,EXCMAT(1,1,1,IBC),
!    &                     GAO1(1,IXYZ),GAO)
!       ENDDO
      ENDDO

      RETURN
      END
