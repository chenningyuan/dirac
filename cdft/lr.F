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
     &                   NXCMAT,
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
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          ISYMOP(*),
     &          IHRMOP(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          SPIN_DENSITY(3),
     &          SPIN_DENSITY_DERVS(3,3),
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

      DO IMAT = 1,NXCMAT !loop over all frequencies
        IF(IHRMOP(IMAT) .EQ. -1) THEN

          JOFF = NBAST*NBAST*NZ*(IMAT-1)

          IF(DOGGA) THEN
            CALL GET_CURRENT(.TRUE.,
     &                       .TRUE.,
     &                       SPIN_DENSITY,
     &                       SPIN_DENSITY_DERVS,
     &                       F,
     &                       G,
     &                       IU,
     &                       IV,
     &                       IX,
     &                       IY,
     &                       ITRANSFORM,
     &                       ISYMOP(IMAT)-1,
     &                       GAO,
     &                       GAO1,
     &                       DTRMAT(1+JOFF),
     &                       BUF)
          ELSE !LDA
            CALL GET_CURRENT(.TRUE.,
     &                       .FALSE.,
     &                       SPIN_DENSITY,
     &                       UNUSED,
     &                       F,
     &                       G,
     &                       IU,
     &                       IV,
     &                       IX,
     &                       IY,
     &                       ITRANSFORM,
     &                       ISYMOP(IMAT)-1,
     &                       GAO,
     &                       UNUSED,
     &                       DTRMAT(1+JOFF),
     &                       BUF)
          ENDIF

          SPIN_DENSITY_NORM = DSQRT(DDOT(3,SPIN_DENSITY,1,
     &                                     SPIN_DENSITY,1))

          IF(.NOT. DOGGA) THEN !LDA

            IF(SPIN_DENSITY_NORM .GT. DFTHRL) THEN

              CALL SDFTPOT1(VXC,WEIGHT,DENSITY,D0)

              DO IXYZ = ISTART,3
                CALL SDFT_LR_LSDA_COLLINEAR(VXC,
     &                                      EXCMAT(1+JOFF),
     &                                      ISYMOP(IMAT),
     &                                      GAO,
     &                                      IXYZ,
     &                                      SPIN_DENSITY(IXYZ))
              ENDDO
            ENDIF

          ELSE !GGA

            TEST_MAX = MAX(
     &                          SPIN_DENSITY_NORM,
     &                     DABS(SPIN_DENSITY_DERVS(1,1)),
     &                     DABS(SPIN_DENSITY_DERVS(1,2)),
     &                     DABS(SPIN_DENSITY_DERVS(1,3)),
     &                     DABS(SPIN_DENSITY_DERVS(2,1)),
     &                     DABS(SPIN_DENSITY_DERVS(2,2)),
     &                     DABS(SPIN_DENSITY_DERVS(2,3)),
     &                     DABS(SPIN_DENSITY_DERVS(3,1)),
     &                     DABS(SPIN_DENSITY_DERVS(3,2)),
     &                     DABS(SPIN_DENSITY_DERVS(3,3))
     &                    )

            IF(TEST_MAX .GT. DFTHRL) THEN

              CALL SDFTPOT1(VXC,WEIGHT,DENSITY,DENSITY_GRADIENT_NORM)

              DO IXYZ = ISTART,3
                CALL SDFT_LR_GGA_COLLINEAR(VXC,
     &                                     EXCMAT(1+JOFF),
     &                                     ISYMOP(IMAT),
     &                                     GAO,
     &                                     GAO1,
     &                                     IXYZ,
     &                                     SPIN_DENSITY(IXYZ),
     &                                     SPIN_DENSITY_DERVS(1,IXYZ),
     &                                     DENSITY_GRADIENT,
     &                                     DENSITY_GRADIENT_NORM)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SDFT_LR_LSDA_COLLINEAR */
      SUBROUTINE SDFT_LR_LSDA_COLLINEAR(VXC,
     &                                  EXCMAT,
     &                                  ISYMOP,
     &                                  GAO,
     &                                  IXYZ,
     &                                  SPIN_DENSITY)
C=======================================================================
C     contribution of spin density to the sigma vector
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

      LOGICAL   LBIT
      DIMENSION VXC(*),
     &          EXCMAT(*),
     &          GAO(*)

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

      VXC_S = VXC(FRR)*SPIN_DENSITY

      CALL CDFT_EXC(VXC_S,F,IU,IV,
     &              IREPOP,EXCMAT(1+IOFF),GAO)
      IF(MC .EQ. 2) THEN
        CALL CDFT_EXC(VXC_S,G,IX,IY,
     &                IREPOP,EXCMAT(1+IOFF),GAO)
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SDFT_LR_GGA_COLLINEAR */
      SUBROUTINE SDFT_LR_GGA_COLLINEAR(VXC,
     &                                 EXCMAT,
     &                                 ISYMOP,
     &                                 GAO,
     &                                 GAO1,
     &                                 IXYZ,
     &                                 SPIN_DENSITY,
     &                                 SPIN_DENSITY_DERVS,
     &                                 DENSITY_GRADIENT,
     &                                 DENSITY_GRADIENT_NORM)
C=======================================================================
C     contribution of spin density to the sigma vector
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

      LOGICAL   LBIT
      DIMENSION VXC(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*),
     &          SPIN_DENSITY_DERVS(3),
     &          DENSITY_GRADIENT(3),
     &          VXC_V(3)

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

      PART1 = DDOT(3,SPIN_DENSITY_DERVS,1,DENSITY_GRADIENT,1)
      PART1 = PART1/DENSITY_GRADIENT_NORM

      PART2 = VXC(FZZ) - VXC(FZ0)/DENSITY_GRADIENT_NORM
      PART2 = PART1*PART2 + VXC(FRZ)*SPIN_DENSITY

      VXC_S = VXC(FRR)*SPIN_DENSITY + VXC(FRZ)*PART1
      DO JXYZ = 1,3
        VXC_V(JXYZ) = (VXC(FZ0)*SPIN_DENSITY_DERVS(JXYZ)
     &             + PART2*DENSITY_GRADIENT(JXYZ))/DENSITY_GRADIENT_NORM
      ENDDO

      CALL DSCAL(3,D2,VXC_V,1) !because of symmetrization later in the code

      CALL CDFT_EXC_GGA(VXC_S,VXC_V,F,IU,IV,
     &                  IREPOP,EXCMAT(1+IOFF),GAO,GAO1)
      IF(MC .EQ. 2) THEN
        CALL CDFT_EXC_GGA(VXC_S,VXC_V,F,IU,IV,
     &                    IREPOP,EXCMAT(1+IOFF),GAO,GAO1)
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CDFT_EXC */
      SUBROUTINE CDFT_EXC(VXC_S,
     &                    F,
     &                    IC1,IC2,
     &                    IREPOP,
     &                    EXCMAT,
     &                    GAO)
C=======================================================================
C     contribution of current or spin density to the sigma vector
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

      DO IREP1 = 0,NBSYM-1

        N1 = NBBAS(IREP1,IC1)
        IF(N1 .GT. 0) THEN

          IREP2 = IBTXOR(IREP1,IREPOP)
          N2 = NBBAS(IREP2,IC2)
          IF(N2 .GT. 0) THEN

            I1 = IBBAS(IREP1,IC1)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CDFT_EXC_GGA */
      SUBROUTINE CDFT_EXC_GGA(VXC_S,
     &                        VXC_V,
     &                        F,
     &                        IC1,IC2,
     &                        IREPOP,
     &                        EXCMAT,
     &                        GAO,
     &                        GAO1)
C=======================================================================
C     contribution of current or spin density to the sigma vector
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

      DO IREP1 = 0,NBSYM-1

        N1 = NBBAS(IREP1,IC1)
        IF(N1 .GT. 0) THEN

          IREP2 = IBTXOR(IREP1,IREPOP)
          N2 = NBBAS(IREP2,IC2)
          IF(N2 .GT. 0) THEN

            I1 = IBBAS(IREP1,IC1)
            I2 = IBBAS(IREP2,IC2)
            DO J = 1,N2
              FACTOR = VXC_S   *GAO (I2+J  )
     &               + VXC_V(1)*GAO1(I2+J,1)
     &               + VXC_V(2)*GAO1(I2+J,2)
     &               + VXC_V(3)*GAO1(I2+J,3)
              FACTOR = FACTOR*F
              CALL DAXPY(N1,FACTOR,GAO(I1+1),1,EXCMAT(I1+1,I2+J),1)
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
