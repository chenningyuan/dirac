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
C     /* deck sdft_lintrans */
      SUBROUTINE SDFT_LINTRANS(DTRMAT,EXCMAT,
     &                         IREPOP,
     &                         GAO,GAO1,
     &                         RHO,RHOGRADNORM,RHOGRAD,WEIGHT,DOGGA,
     &                         BUF)
C  .....................................................................
C  .  description:                                                     .
C  .  lets the "spin density" contribute to linear response            .
C  .  adaptation of DFTLTR in dftaux.F                                 .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     march 2006      .
C  .....................................................................
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "inforb.h"
#include "dgroup.h"
#include "dftcom.h"
#include "cdft.h"
#include "../dft/functionals.h"
C
      LOGICAL DOGGA
C
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*),
     &          RHOGRAD(3),
     &          DMAGNETIZATION(3),
     &          VXC(DERVS2),
     &          BUF(*)
C
#include "ibtfun.h"
C
C
C     get modified magnetization
C     ==========================
C
      CALL GETCURRENTDENSITY(.FALSE.,DUMMY,
     &                       .TRUE.,DMAGNETIZATION,
     &                       IREPOP,GAO,DTRMAT,BUF)
C
      IF(DOGGA) THEN
C
C
C       G G A
C       =====
C
        CALL QUIT('SDFT linear response not implemented for GGA')
C
      ELSE
C
C
C       L D A
C       =====
C
        DNORM = DSQRT(DDOT(3,DMAGNETIZATION,1,DMAGNETIZATION,1))
        IF(DNORM.GT.DFTHRL) THEN
C
          CALL SDFTPOT1(VXC,WEIGHT,RHO,D0)
C
          ITRANSFORM = 2
C
          DO I = 1,3
            IQ = IQTOPQ(5-I,ISYMAX(I,ITRANSFORM))
            IOFF = NBAST*NBAST
            IOFF = IOFF*(IQ-1)
C
            VXC_S = VXC(FRR)*DMAGNETIZATION(I)
            MYIREPOP = IBTXOR(IREPOP,ISYMAX(I,ITRANSFORM))
            CALL DOEXC1(VXC_S,MYIREPOP,EXCMAT(1+IOFF),GAO)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C
      RETURN
      END
