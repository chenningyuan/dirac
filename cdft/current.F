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
C     /* deck getcurrentdensity */
      SUBROUTINE GETCURRENTDENSITY(LCURRENTDENSITY,CURRENTDENSITY,
     &                             LMAGNETIZATION,DMAGNETIZATION,
     &                             IREPOP,GAO,DMAT,BUF)
C=======================================================================
C     get current density or magnetization density in point
C
C     routine returns one (or two) vector(s) with
C     density(1) =  density_x
C     density(2) =  density_y
C     density(3) =  density_z
C-----------------------------------------------------------------------
C     radovan bast                              last revision: july 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "dcbgen.h"
#include "cdft.h"
C-----------------------------------------------------------------------
      LOGICAL LCURRENTDENSITY,LMAGNETIZATION,NONCOLLINEAR
C
      DIMENSION CURRENTDENSITY(3),
     &          DMAGNETIZATION(3),
     &          GAO(*),DMAT(*),BUF(*)
C-----------------------------------------------------------------------
      IF(LMAGNETIZATION) THEN
C       NONCOLLINEAR = .NOT.COLLINEAR
        NONCOLLINEAR = .TRUE.
C
        F  = D1
        G  = D1
        IU = 1
        IV = 1
        IX = 2
        IY = 2
C       transforms like rotations, hence:
        ITRANSFORM = 2
        CALL GETCURRENTDENSITY_SUB(DMAGNETIZATION,NONCOLLINEAR,
     &                             F,G,IU,IV,IX,IY,IREPOP,ITRANSFORM,
     &                             GAO,DMAT,BUF)
      ENDIF
C
      IF(LCURRENTDENSITY) THEN
        F  = -CVAL
        G  = -CVAL
        IU = 1
        IV = 2
        IX = 2
        IY = 1
C       transforms like coordinates, hence:
        ITRANSFORM = 1
        CALL GETCURRENTDENSITY_SUB(CURRENTDENSITY,.TRUE.,
     &                             F,G,IU,IV,IX,IY,IREPOP,ITRANSFORM,
     &                             GAO,DMAT,BUF)
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck getcurrentdensity_sub */
      SUBROUTINE GETCURRENTDENSITY_SUB(CURRENTDENSITY,NONCOLLINEAR,
     &                                 F,G,IU,IV,IX,IY,
     &                                 IREPOP,ITRANSFORM,
     &                                 GAO,DMAT,BUF)
C=======================================================================
C     get contributions to current or magnetization density
C-----------------------------------------------------------------------
C     radovan bast                              last revision: july 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
C-----------------------------------------------------------------------
      LOGICAL NONCOLLINEAR
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
#endif
C
      DIMENSION CURRENTDENSITY(3),
     &          GAO(*),DMAT(*),BUF(*)
C-----------------------------------------------------------------------
#include "ibtfun.h"
C-----------------------------------------------------------------------
      IF(NONCOLLINEAR) THEN
        ISTART = 1
      ELSE
        CURRENTDENSITY(1) = D0
        CURRENTDENSITY(2) = D0
        ISTART = 3
      ENDIF
C
      DO I = ISTART,3
#if defined (VAR_G77) || defined (SYS_SX)
        IOFF = NBAST*NBAST
#else
        IOFF = NTBAS(0)*NTBAS(0)
#endif
        IQ = IQTOPQ(5-I,ISYMAX(I,2))
        IOFF = IOFF*(IQ-1)
C
C       TEMP1 = D0
        CALL GETCURRENT(TEMP1,IU,IV,
     &                  IREPOP,ITRANSFORM,I,
     &                  GAO,DMAT(1+IOFF),BUF)
C
C       TEMP2 = D0
        CALL GETCURRENT(TEMP2,IX,IY,
     &                  IREPOP,ITRANSFORM,I,
     &                  GAO,DMAT(1+IOFF),BUF)
C
        CURRENTDENSITY(I) = F*TEMP1 + G*TEMP2
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck getcurrent */
      SUBROUTINE GETCURRENT(CURRENT,IC1,IC2,
     &                      IREPOP,ITRANSFORM,IXYZ,
     &                      GAO,DMAT,BUF)
C=======================================================================
C     get contributions to current or magnetization density
C-----------------------------------------------------------------------
C     radovan bast                              last revision: july 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
C-----------------------------------------------------------------------
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),
     &          DMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),
     &          DMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION BUF(*)
C-----------------------------------------------------------------------
#include "ibtfun.h"
C-----------------------------------------------------------------------
      CURRENT = D0
C
      DO IREP1 = 0,NBSYM-1
C
        N1 = NBBAS(IREP1,IC1)
        IF(N1.GT.0) THEN
C
          I1 = IBBAS(IREP1,IC1)+1
C
          INTERMEDIATE = ISYMAX(IXYZ,ITRANSFORM)
          INTERMEDIATE = IBTXOR(INTERMEDIATE,IREPOP)
          IREP2 = IBTXOR(IREP1,INTERMEDIATE)
C
          N2 = NBBAS(IREP2,IC2)
          IF(N2.GT.0) THEN
C
            I2 = IBBAS(IREP2,IC2)+1
C
            CALL DGEMV('N',NBBAS(IREP1,IC1),NBBAS(IREP2,IC2),D1,
     &                 DMAT(I1,I2),NTBAS(0),GAO(I2),1,D0,BUF,1)
            CURRENT = CURRENT + DDOT(NBBAS(IREP1,IC1),GAO(I1),1,BUF,1)
C
          ENDIF
        ENDIF
C
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
