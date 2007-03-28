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
C     /* EXC_M */
      SUBROUTINE EXC_M(IS_GGA,
     &                 VXC_S,
     &                 VXC_V,
     &                 IXYZ,
     &                 ISYMOP,
     &                 EXCMAT,
     &                 GAO,
     &                 GAO1)
C=======================================================================
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
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "cdft.h"
#include "dftcom.h"
#include "../dft/functionals.h"

      LOGICAL   IS_GGA,
     &          LBIT
      DIMENSION VXC_V(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*)

#include "ibtfun.h"

      F = D1
      G = D1
      IF(LBIT(LBIT_CDFT,4)) G = -D1 !operator is \beta\sigma instead of \sigma
      IU         = 1
      IV         = 1
      IX         = 2
      IY         = 2
      ITRANSFORM = 2

      IREPOP = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),ISYMOP-1)
      IQ     = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
      IOFF   = NBAST*NBAST*(IQ-1)

      CALL EXC_GENERAL(IS_GGA,VXC_S,VXC_V,F,IU,IV,
     &                 IREPOP,EXCMAT(1+IOFF),GAO,GAO1)
      IF(MC .EQ. 2) THEN
        CALL EXC_GENERAL(IS_GGA,VXC_S,VXC_V,G,IX,IY,
     &                   IREPOP,EXCMAT(1+IOFF),GAO,GAO1)
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* EXC_GENERAL */
      SUBROUTINE EXC_GENERAL(IS_GGA,
     &                       VXC_S,
     &                       VXC_V,
     &                       F,
     &                       IC1,IC2,
     &                       IREPOP,
     &                       EXCMAT,
     &                       GAO,
     &                       GAO1)
C=======================================================================
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
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"

      LOGICAL   IS_GGA
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION EXCMAT(NBAST,NBAST)
#else
      DIMENSION EXCMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION VXC_V(3),
     &          GAO(*),
     &          GAO1(*)

#include "ibtfun.h"

      DO IREP1 = 0,NBSYM-1

        N1 = NBBAS(IREP1,IC1)
        IF(N1 .GT. 0) THEN

          IREP2 = IBTXOR(IREP1,IREPOP)
          N2 = NBBAS(IREP2,IC2)
          IF(N2 .GT. 0) THEN

            I1 = IBBAS(IREP1,IC1)
            I2 = IBBAS(IREP2,IC2)
            DO J = 1,N2

              FACTOR = VXC_S*GAO(I2+J)

              IF(IS_GGA) THEN
                DO IXYZ = 1,3
                  FACTOR = FACTOR
#if defined (VAR_G77) || defined (SYS_SX)
     &                   + VXC_V(IXYZ)*GAO1(I2+J+(IXYZ-1)*NBAST)
#else
     &                   + VXC_V(IXYZ)*GAO1(I2+J+(IXYZ-1)*NTBAS(0))
#endif
                ENDDO
              ENDIF

              FACTOR = FACTOR*F
              CALL DAXPY(N1,FACTOR,GAO(I1+1),1,EXCMAT(I1+1,I2+J),1)
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
