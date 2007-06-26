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
      SUBROUTINE SDFT_QR_ALTERNATIVE(ADD_GR,
     &                               ADD_S,
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

      LOGICAL   ADD_GR,
     &          ADD_S,
     &          ADD_GS,
     &          LBIT
      DIMENSION DTRMAT(NBAST,NBAST,NZ,3),
     &          EXCMAT(NBAST,NBAST,NZ,3),
     &          ISYMOP(3),
     &          IHRMOP(3),
     &          GAO(*),
     &          GAO1(*),
     &          GR(3),GR_b(3),GR_c(3),GR_bc(3),
     &          BUF(*),
     &          VXC_VECTOR(3),
     &          VXC(99) !dirty fix

#include "ibtfun.h"

      ibc = 1
      CALL GETRHO(R_b,
     &            ISYMOP(ibc)-1,
     &            GAO,
     &            DTRMAT(1,1,1,ibc),
     &            BUF)
      ibc = 2
      CALL GETRHO(R_c,
     &            ISYMOP(ibc)-1,
     &            GAO,
     &            DTRMAT(1,1,1,ibc),
     &            BUF)
      ibc = 3
      CALL GETRHO(R_bc,
     &            ISYMOP(ibc)-1,
     &            GAO,
     &            DTRMAT(1,1,1,ibc),
     &            BUF)

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









      DMAGIC = D2

      Z_b  =         DDOT(3,GR_b ,1,GR   ,1) + DDOT(3,GR   ,1,GR_b ,1)
      Z_c  =         DDOT(3,GR_c ,1,GR   ,1) + DDOT(3,GR   ,1,GR_c ,1)
      Z_bc =         DDOT(3,GR_bc,1,GR   ,1) + DDOT(3,GR   ,1,GR_bc,1)
     &     + DMAGIC*(DDOT(3,GR_b ,1,GR_c ,1) + DDOT(3,GR_c ,1,GR_b ,1))

      Y_b  =         DDOT(3,GR_b ,1,GS   ,1) + DDOT(3,GR   ,1,GS_b ,1)
      Y_c  =         DDOT(3,GR_c ,1,GS   ,1) + DDOT(3,GR   ,1,GS_c ,1)
      Y_bc =         DDOT(3,GR_bc,1,GS   ,1) + DDOT(3,GR   ,1,GS_bc,1)
     &     + DMAGIC*(DDOT(3,GR_b ,1,GS_c ,1) + DDOT(3,GR_c ,1,GS_b ,1))

      X_b  =         DDOT(3,GS_b ,1,GS   ,1) + DDOT(3,GS   ,1,GS_b ,1)
      X_c  =         DDOT(3,GS_c ,1,GS   ,1) + DDOT(3,GS   ,1,GS_c ,1)
      X_bc =         DDOT(3,GS_bc,1,GS   ,1) + DDOT(3,GS   ,1,GS_bc,1)
     &     + DMAGIC*(DDOT(3,GS_b ,1,GS_c ,1) + DDOT(3,GS_c ,1,GS_b ,1))

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

      IF(ADD_GR) THEN
        CALL DFTPOT2(VXC,W,R,GR_NORM)
      ELSE
        CALL DFTPOT2(VXC,W,R,D0)
      ENDIF

      IF(ADD_S) THEN
        IF(ADD_GS) THEN
          CALL SDFTPOT2(VXC,W,R,GR_NORM)
        ELSE
          CALL SDFTPOT2(VXC,W,R,D0)
        ENDIF
      ENDIF






      R_pq_bc = R_pq_bc + VXC(fR0)
      Z_pq_bc = Z_pq_bc + VXC(fZ0)

      R_pq_b  = R_pq_b
     &        + R_c *VXC(fRR)
     &        + Z_c *VXC(fRZ)

      Z_pq_b  = Z_pq_b
     &        + R_c *VXC(fRZ)
     &        + Z_c *VXC(fZZ)

      R_pq_c  = R_pq_c
     &        + R_b *VXC(fRR)
     &        + Z_b *VXC(fRZ)

      Z_pq_c  = Z_pq_c
     &        + R_b *VXC(fRZ)
     &        + Z_b *VXC(fZZ)

      R_pq    = R_pq  
     &        + R_bc*VXC(fRR)
     &        + Z_bc*VXC(fRZ)

      Z_pq    = Z_pq  
     &        + R_bc*VXC(fRZ)
     &        + Z_bc*VXC(fZZ)

      R_pq    = R_pq
     &        + (R_b*R_c + R_c*R_b)*VXC(fRRR)
     &        + (R_b*Z_c + R_c*Z_b)*VXC(fRRZ)
     &        + (Z_b*R_c + Z_c*R_b)*VXC(fRRZ)
     &        + (Z_b*Z_c + Z_c*Z_b)*VXC(fRZZ)

      Z_pq    = Z_pq
     &        + (R_b*R_c + R_c*R_b)*VXC(fRRZ)
     &        + (R_b*Z_c + R_c*Z_b)*VXC(fRZZ)
     &        + (Z_b*R_c + Z_c*R_b)*VXC(fRZZ)
     &        + (Z_b*Z_c + Z_c*Z_b)*VXC(fZZZ)









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
