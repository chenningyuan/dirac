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
C     /* deck sdft_linear_trans */
      SUBROUTINE SDFT_LINEAR_TRANS(DTRMAT,
     &                             EXCMAT,
     &                             IREPOP,
     &                             GAO,DENSITY,WEIGHT,DOGGA,BUF)
C=======================================================================
C     contribution of spin density to E_xc
C     adaptation of DFTLTR in dftaux.F
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "inforb.h"
#include "dgroup.h"
#include "dftcom.h"
#include "cdft.h"
#include "../dft/functionals.h"
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LOGICAL DOGGA
C
      DIMENSION DTRMAT(*),EXCMAT(*),GAO(*),BUF(*),
     &          SPIN_DENSITY(3),VXC(DERVS2)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "ibtfun.h"
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL GETCURRENTDENSITY(.FALSE.,DUMMY,
     &                       .TRUE.,SPIN_DENSITY,
     &                       IREPOP,GAO,DTRMAT,BUF)
C
      IF(DOGGA) CALL QUIT('GGA SDFT linear response not implemented')
C
      SPIN_DENSITY_NORM = DSQRT(DDOT(3,SPIN_DENSITY,1,SPIN_DENSITY,1))
      IF(SPIN_DENSITY_NORM.GT.DFTHRL) THEN
C
        CALL SDFTPOT1(VXC,WEIGHT,DENSITY,D0)
C
        ISTART = 1
        IF(COLLINEAR) ISTART = 3
C
        DO IXYZ = ISTART,3
          IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
          IOFF = NBAST*NBAST
          IOFF = IOFF*(IQ-1)
C
          VXC_S = VXC(FRR)*SPIN_DENSITY(IXYZ)
C
          CALL SDFT_EXC(VXC_S,IREPOP,IXYZ,EXCMAT(1+IOFF),GAO)
        ENDDO
C
      ENDIF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck sdft_exc */
      SUBROUTINE SDFT_EXC(VXC_S,IREPOP,IXYZ,EXCMAT,GAO)
C=======================================================================
C     contribution of spin density to E_xc
C     adaptation of DOEXC1 in dftdrv.F
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),
     &          EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),
     &          EXCMAT(NTBAS(0),NTBAS(0))
#endif
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "ibtfun.h"
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO IREP1 = 0,NBSYM-1
        INTERMEDIATE = ISYMAX(IXYZ,2)
        INTERMEDIATE = IBTXOR(INTERMEDIATE,IREPOP)
        IREP2 = IBTXOR(IREP1,INTERMEDIATE)
C
        DO IC = 1,MC
C
          I1 = IBBAS(IREP1,IC)
          N1 = NBBAS(IREP1,IC)
          I2 = IBBAS(IREP2,IC)
          N2 = NBBAS(IREP2,IC)
C
          IF((N1.GT.0).AND.(N2.GT.0)) THEN
C
            DO J = 1,N2
              CALL DAXPY(N1,VXC_S*GAO(I2+J),
     &                   GAO(I1+1),1,EXCMAT(I1+1,I2+J),1)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
