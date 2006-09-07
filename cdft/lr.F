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
      SUBROUTINE SDFT_LR(DTRMAT,
     &                   EXCMAT,
     &                   IREPOP,
     &                   GAO,DENSITY,WEIGHT,DOGGA,BUF)
C    &                   GAO,GAO1,
C    &                   DENSITY,DENSITY_GRADIENT_NORM,
C    &                   WEIGHT,DOGGA,BUF)
C=======================================================================
C     contribution of spin density to g_pg (upon integration G_pg)
C     adaptation of DFTLTR in dftaux.F
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "inforb.h"
#include "dgroup.h"
#include "dftcom.h"
#include "cdft.h"
#include "../dft/functionals.h"
C
      LOGICAL   DOGGA
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          GAO(*),
C    &          GAO1(*),
     &          BUF(*),
     &          SPIN_DENSITY(3),
     &          SPIN_DENSITY_GRADIENT(3,3),
     &          VXC(DERVS2),
     &          VXC_V(3)
C
#include "ibtfun.h"
C
C
      CALL GETCURRENTDENSITY(.FALSE.,DUMMY,
     &                       .TRUE.,SPIN_DENSITY,
     &                       IREPOP,GAO,DTRMAT,BUF)
C
      ISTART = 1
      IF(COLLINEAR) ISTART = 3
C
      SPIN_DENSITY_NORM = DSQRT(DDOT(3,SPIN_DENSITY,1,SPIN_DENSITY,1))
C
      IF(.NOT. DOGGA) THEN !LDA
C
        IF(SPIN_DENSITY_NORM .GT. DFTHRL) THEN
C
          CALL SDFTPOT1(VXC,WEIGHT,DENSITY,D0)
C
          DO IXYZ = ISTART,3
            IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
            IOFF = NBAST*NBAST*(IQ-1)
C
            VXC_S = VXC(FRR)*SPIN_DENSITY(IXYZ)
            CALL SDFT_EXC(VXC_S,IREPOP,IXYZ,EXCMAT(1+IOFF),GAO)
          ENDDO
        ENDIF
C
      ELSE !GGA
C
        CALL QUIT('GGA SDFT linear response not implemented')
C
C       CALL GETCURRENTDENSITY_GRADIENT(.FALSE.,DUMMY,
C    &                                  .TRUE.,SPIN_DENSITY_GRADIENT,
C    &                                  IREPOP,GAO,GAO1,DTRMAT,BUF)
C
C                                        checkmeout
        CALL SDFTPOT1(VXC,WEIGHT,DENSITY,DENSITY_GRADIENT_NORM)
C
        DO IXYZ = ISTART,3
          IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
          IOFF = NBAST*NBAST*(IQ-1)
C
C         BR = (B3(1)*RHG(1)+ B3(2)*RHG(2)+ B3(3)*RHG(3))/RHOGRD
C         FACR = VXC(FRZ)*SPIN_DENSITY(IXYZ) + (VXC(FZZ)-VXC(FZ0)/RHOGRD)*BR
C
          VXC_S = VXC(FRR)*SPIN_DENSITY(IXYZ)
C    &                                        + VXC(FRZ)*BR
          DO JXYZ = 1,3
C           VXC_V(JXYZ) = (VXC(FZ0)*B3(JXYZ)+FACR*RHG(JXYZ))/RHOGRD
          ENDDO
C
C         CALL SDFT_EXC_GGA(VXC_S,VXC_V,
C    &                      IREPOP,IXYZ,EXCMAT(1+IOFF),GAO,GAO1)
        ENDDO
      ENDIF
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck sdft_exc */
      SUBROUTINE SDFT_EXC(VXC_S,IREPOP,IXYZ,EXCMAT,GAO)
C=======================================================================
C     contribution of spin density to g_pg (upon integration G_pg)
C     adaptation of DOEXC1 in dftdrv.F
C-----------------------------------------------------------------------
C     radovan bast                         last revision: september 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),
     &          EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),
     &          EXCMAT(NTBAS(0),NTBAS(0))
#endif
C
#include "ibtfun.h"
C
C
      INTERMEDIATE = IBTXOR(ISYMAX(IXYZ,2),IREPOP)
C
      DO IREP1 = 0,NBSYM-1
        DO IC = 1,MC
C
          N1 = NBBAS(IREP1,IC)
          IF(N1.GT.0) THEN
C
            I1 = IBBAS(IREP1,IC)
            IREP2 = IBTXOR(IREP1,INTERMEDIATE)
            N2 = NBBAS(IREP2,IC)
            IF(N2.GT.0) THEN
C
              I2 = IBBAS(IREP2,IC)
              DO J = 1,N2
                CALL DAXPY(N1,VXC_S*GAO(I2+J),
     &                     GAO(I1+1),1,EXCMAT(I1+1,I2+J),1)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
C
      RETURN
      END
