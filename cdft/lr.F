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
C     /* SDFT_LR */
      SUBROUTINE SDFT_LR(DTRMAT,
     &                   EXCMAT,
     &                   IREPOP,
     &                   GAO,GAO1,
     &                   DENSITY,
     &                   DENSITY_GRADIENT,
     &                   DENSITY_GRADIENT_NORM,
     &                   WEIGHT,DOGGA,BUF)
C=======================================================================
C     contribution of spin density to g_pq (upon integration G_pq)
C     adaptation of DFTLTR in dftaux.F
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
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
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          SPIN_DENSITY(3),
     &          SPIN_DENSITY_DERVS(3,3),
     &          VXC(DERVS2),
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

      CALL GET_CURRENT(2,
     &                 SPIN_DENSITY,
     &                 DUMMY,
     &                 IREPOP,GAO,DUMMY,DTRMAT,BUF)

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      SPIN_DENSITY_NORM = DSQRT(DDOT(3,SPIN_DENSITY,1,SPIN_DENSITY,1))

      IF(.NOT. DOGGA) THEN !LDA
        IF(SPIN_DENSITY_NORM .GT. DFTHRL) THEN

          CALL SDFTPOT1(VXC,WEIGHT,DENSITY,D0)

          DO IXYZ = ISTART,3
            IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
            IOFF = NBAST*NBAST*(IQ-1)

            VXC_S = VXC(FRR)*SPIN_DENSITY(IXYZ)

            CALL CDFT_EXC(VXC_S,
     &                    F,IU,IV,ITRANSFORM,IXYZ,
     &                    GAO,
     &                    IREPOP,EXCMAT(1+IOFF),BUF)
            CALL CDFT_EXC(VXC_S,
     &                    G,IX,IY,ITRANSFORM,IXYZ,
     &                    GAO,
     &                    IREPOP,EXCMAT(1+IOFF),BUF)
          ENDDO
        ENDIF

      ELSE !GGA

        CALL GET_CURRENT(4,
     &                   DUMMY,
     &                   SPIN_DENSITY_DERVS,
     &                   IREPOP,GAO,GAO1,DTRMAT,BUF)

        TEST_MAX = MAX(
     &                 SPIN_DENSITY_NORM,
     &                 SPIN_DENSITY_DERVS(1,1),
     &                 SPIN_DENSITY_DERVS(1,2),
     &                 SPIN_DENSITY_DERVS(1,3),
     &                 SPIN_DENSITY_DERVS(2,1),
     &                 SPIN_DENSITY_DERVS(2,2),
     &                 SPIN_DENSITY_DERVS(2,3),
     &                 SPIN_DENSITY_DERVS(3,1),
     &                 SPIN_DENSITY_DERVS(3,2),
     &                 SPIN_DENSITY_DERVS(3,3)
     &                )

        IF(TEST_MAX .GT. DFTHRL) THEN

          CALL SDFTPOT1(VXC,WEIGHT,DENSITY,DENSITY_GRADIENT_NORM)

          DO IXYZ = ISTART,3
            IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
            IOFF = NBAST*NBAST*(IQ-1)

            FACTOR_1 = D0
            DO JXYZ = 1,3
              FACTOR_1 = FACTOR_1
     &            + SPIN_DENSITY_DERVS(IXYZ,JXYZ)*DENSITY_GRADIENT(JXYZ)
            ENDDO
            FACTOR_1 = FACTOR_1/DENSITY_GRADIENT_NORM

            FACTOR_2 = VXC(FZZ) - VXC(FZ0)/DENSITY_GRADIENT_NORM
            FACTOR_2 = FACTOR_2*FACTOR_1
            FACTOR_2 = FACTOR_2 + VXC(FRZ)*SPIN_DENSITY(IXYZ)

            VXC_S = VXC(FRR)*SPIN_DENSITY(IXYZ) + VXC(FRZ)*FACTOR_1
            DO JXYZ = 1,3
              VXC_V(JXYZ) = (VXC(FZ0)*SPIN_DENSITY_DERVS(IXYZ,JXYZ)
     &          + FACTOR_2*DENSITY_GRADIENT(JXYZ))/DENSITY_GRADIENT_NORM
              VXC_V(JXYZ) = VXC_V(JXYZ)*D2 !because of SYMSYM later in the code
            ENDDO

            CALL SDFT_EXC_GGA(VXC_S,VXC_V,
     &                        IREPOP,IXYZ,EXCMAT(1+IOFF),GAO,GAO1)
          ENDDO
        ENDIF
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CDFT_LR */
      SUBROUTINE CDFT_LR(DTRMAT,
     &                   EXCMAT,
     &                   IREPOP,
     &                   GAO,GAO1,
     &                   DENSITY,
     &                   DENSITY_GRADIENT,
     &                   DENSITY_GRADIENT_NORM,
     &                   WEIGHT,DOGGA,BUF)
C=======================================================================
C     contribution of current density to g_pq (upon integration G_pq)
C     adaptation of DFTLTR in dftaux.F
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (D23 = 2.0D0/3.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "inforb.h"
#include "dgroup.h"
#include "dftcom.h"
#include "dcbgen.h"
#include "cdft.h"
#include "../dft/functionals.h"

      LOGICAL   DOGGA, LBIT
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          CURRENT_DENSITY(3),
     &          CURRENT_DENSITY_DERVS(3,3),
     &          VXC(DERVS2),
     &          VXC_V(3)

#include "ibtfun.h"

      IF(.NOT. LBIT(LBIT_CDFT,5)) !X\alpha type current density contribution not requested
     &         CALL QUIT('problem in cdft linear response')

      F  = -CVAL
      G  = -CVAL
      IU = 1
      IV = 2
      IX = 2
      IY = 1
      ITRANSFORM = 1

      CALL GET_CURRENT(3,
     &                 CURRENT_DENSITY,
     &                 DUMMY,
     &                 IREPOP,GAO,DUMMY,DTRMAT,BUF)

      DO IXYZ = 1,3
        CURRENT_DENSITY(IXYZ) = DABS(CURRENT_DENSITY(IXYZ))
        IF(CURRENT_DENSITY(IXYZ) .LT. DFTHRL) GOTO 1

        IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
        IOFF = NBAST*NBAST*(IQ-1)

        VXC_S = CURRENT_DENSITY(IXYZ)**(-D23)
        VXC_S = VXC_S*WEIGHT*CDFT_ALPHA

        CALL CDFT_EXC(VXC_S,
     &                F,IU,IV,ITRANSFORM,IXYZ,
     &                GAO,
     &                IREPOP,EXCMAT(1+IOFF),BUF)
        CALL CDFT_EXC(VXC_S,
     &                G,IX,IY,ITRANSFORM,IXYZ,
     &                GAO,
     &                IREPOP,EXCMAT(1+IOFF),BUF)
 1      CONTINUE
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SDFT_EXC_GGA */
      SUBROUTINE SDFT_EXC_GGA(VXC_S,VXC_V,
     &                        IREPOP,IXYZ,EXCMAT,GAO,GAO1)
C=======================================================================
C     contribution of spin density to g_pq (upon integration G_pq)
C     adaptation of DOEXC1 in dftdrv.F
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
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

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),
     &          GAO1(NBAST,3),
     &          EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),
     &          GAO1(NTBAS(0),3),
     &          EXCMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION VXC_V(3)

#include "ibtfun.h"

      INTERMEDIATE = IBTXOR(ISYMAX(IXYZ,2),IREPOP)

      DO IREP1 = 0,NBSYM-1
        DO IC = 1,MC

          N1 = NBBAS(IREP1,IC)
          IF(N1 .GT. 0) THEN

            I1 = IBBAS(IREP1,IC)
            IREP2 = IBTXOR(IREP1,INTERMEDIATE)
            N2 = NBBAS(IREP2,IC)
            IF(N2 .GT. 0) THEN

              I2 = IBBAS(IREP2,IC)
              DO J = 1,N2
                FACTOR = VXC_S   *GAO(I2+J)
     &                 + VXC_V(1)*GAO1(I2+J,1)
     &                 + VXC_V(2)*GAO1(I2+J,2)
     &                 + VXC_V(3)*GAO1(I2+J,3)
                CALL DAXPY(N1,FACTOR,GAO(I1+1),1,EXCMAT(I1+1,I2+J),1)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CDFT_EXC */
      SUBROUTINE CDFT_EXC(VXC_S,
     &                    F,IC1,IC2,ITRANSFORM,IXYZ,
     &                    GAO,
     &                    IREPOP,EXCMAT,BUF)
C=======================================================================
C     contribution of current or spin density to
C     to g_pq (upon integration G_pq)
C     adaptation of DOEXC1 in dftdrv.F
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
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

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAO(NBAST),
     &          EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),
     &          EXCMAT(NTBAS(0),NTBAS(0))
#endif

#include "ibtfun.h"

      INTERMEDIATE = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),IREPOP)

      DO IREP1 = 0,NBSYM-1

        N1 = NBBAS(IREP1,IC1)
        IF(N1 .GT. 0) THEN

          I1 = IBBAS(IREP1,IC1)
          IREP2 = IBTXOR(IREP1,INTERMEDIATE)
          N2 = NBBAS(IREP2,IC2)
          IF(N2 .GT. 0) THEN

            I2 = IBBAS(IREP2,IC2)
            DO J = 1,N2
              FACTOR = F*VXC_S*GAO(I2+J)
              CALL DAXPY(N1,FACTOR,GAO(I1+1),1,EXCMAT(I1+1,I2+J),1)
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
