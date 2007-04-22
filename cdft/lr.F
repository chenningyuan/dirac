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
     &                   IS_GGA,
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
#include "dcbham.h"
#include "dgroup.h"
#include "inforb.h"
#include "cdft.h"
#include "dftcom.h"
#include "../dft/functionals.h"

      LOGICAL   IS_GGA,
     &          LBIT
      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          ISYMOP(*),
     &          IHRMOP(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          FO_SPIN_DENSITY(3),
     &          FO_SPIN_DENSITY_DERVS(3,3),
     &          VXC(99) !dirty fix

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      DO IMAT = 1,NXCMAT !loop over all frequencies
        IF(IHRMOP(IMAT) .EQ. -1) THEN

          JOFF = NBAST*NBAST*NZ*(IMAT-1)

          CALL GET_M(.TRUE.,
     &               FO_SPIN_DENSITY,
     &               ISYMOP(IMAT)-1,
     &               GAO,
     &               UNUSED,
     &               DTRMAT(1+JOFF),
     &               BUF)

          IF(IS_GGA) THEN
            CALL GET_M_DERVS(.TRUE.,
     &                       FO_SPIN_DENSITY_DERVS,
     &                       ISYMOP(IMAT)-1,
     &                       GAO,
     &                       GAO1,
     &                       DTRMAT(1+JOFF),
     &                       BUF)
          ENDIF

          FO_SPIN_DENSITY_NORM = DSQRT(DDOT(3,FO_SPIN_DENSITY,1,
     &                                        FO_SPIN_DENSITY,1))

          IF(.NOT. IS_GGA) THEN !LDA

            IF(FO_SPIN_DENSITY_NORM .GT. DFTHRL) THEN

              CALL SDFTPOT1(VXC,WEIGHT,DENSITY,D0)

!             CALL DSCAL(3,VXC(FRR),FO_SPIN_DENSITY,1) !could also be done this way
!             CALL GET_M(.FALSE.,
!    &                   FO_SPIN_DENSITY,
!    &                   ISYMOP(IMAT)-1,
!    &                   GAO,
!    &                   UNUSED,
!    &                   EXCMAT(1+JOFF),
!    &                   BUF)

              DO IXYZ = ISTART,3
                CALL SDFT_LR_LSDA_COLLINEAR(VXC,
     &                                      EXCMAT(1+JOFF),
     &                                      ISYMOP(IMAT),
     &                                      GAO,
     &                                      IXYZ,
     &                                      FO_SPIN_DENSITY(IXYZ))
              ENDDO
            ENDIF

          ELSE !GGA

            TEST_MAX = MAX(
     &                          FO_SPIN_DENSITY_NORM,
     &                     DABS(FO_SPIN_DENSITY_DERVS(1,1)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(1,2)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(1,3)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(2,1)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(2,2)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(2,3)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(3,1)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(3,2)),
     &                     DABS(FO_SPIN_DENSITY_DERVS(3,3))
     &                    )

            IF(TEST_MAX .GT. DFTHRL) THEN

              CALL SDFTPOT1(VXC,WEIGHT,DENSITY,DENSITY_GRADIENT_NORM)

              DO IXYZ = ISTART,3
               CALL SDFT_LR_GGA_COLLINEAR(VXC,
     &                                    EXCMAT(1+JOFF),
     &                                    ISYMOP(IMAT),
     &                                    GAO,
     &                                    GAO1,
     &                                    IXYZ,
     &                                    FO_SPIN_DENSITY(IXYZ),
     &                                    FO_SPIN_DENSITY_DERVS(1,IXYZ),
     &                                    DENSITY_GRADIENT,
     &                                    DENSITY_GRADIENT_NORM)
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
     &                                  FO_SPIN_DENSITY)
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

      DIMENSION VXC(*),
     &          EXCMAT(*),
     &          GAO(*)

      VXC_S = VXC(FRR)*FO_SPIN_DENSITY

      CALL EXC_M(.FALSE.,
     &           VXC_S,
     &           UNUSED,
     &           IXYZ,
     &           ISYMOP,
     &           EXCMAT,
     &           GAO,
     &           UNUSED)

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
     &                                 FO_SPIN_DENSITY,
     &                                 FO_SPIN_DENSITY_DERVS,
     &                                 DENSITY_GRADIENT,
     &                                 DENSITY_GRADIENT_NORM)
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

      DIMENSION VXC(*),
     &          EXCMAT(*),
     &          GAO(*),
     &          GAO1(*),
     &          FO_SPIN_DENSITY_DERVS(3),
     &          DENSITY_GRADIENT(3),
     &          VXC_V(3)

      PART1 = DDOT(3,FO_SPIN_DENSITY_DERVS,1,DENSITY_GRADIENT,1)
      PART1 = PART1/DENSITY_GRADIENT_NORM

      PART2 = VXC(FZZ) - VXC(FZ0)/DENSITY_GRADIENT_NORM
      PART2 = PART1*PART2 + VXC(FRZ)*FO_SPIN_DENSITY

      VXC_S = VXC(FRR)*FO_SPIN_DENSITY + VXC(FRZ)*PART1
      DO JXYZ = 1,3
        VXC_V(JXYZ) = (VXC(FZ0)*FO_SPIN_DENSITY_DERVS(JXYZ)
     &             + PART2*DENSITY_GRADIENT(JXYZ))/DENSITY_GRADIENT_NORM
      ENDDO

      CALL DSCAL(3,D2,VXC_V,1) !because of symmetrization later in the code

      CALL EXC_M(.TRUE.,
     &           VXC_S,
     &           VXC_V,
     &           IXYZ,
     &           ISYMOP,
     &           EXCMAT,
     &           GAO,
     &           GAO1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CDFT_LR */
      SUBROUTINE CDFT_LR(DTRMAT,
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
     &                   IS_GGA,
     &                   BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (DP13 = D1/D3)

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

      DIMENSION DTRMAT(*),
     &          EXCMAT(*),
     &          ISYMOP(*),
     &          IHRMOP(*),
     &          GAO(*),
     &          GAO1(*),
     &          DENSITY_GRADIENT(3),
     &          BUF(*),
     &          FO_VORTICITY(3)

      DO IMAT = 1,NXCMAT !loop over all frequencies
        IF(IHRMOP(IMAT) .EQ. -1) THEN

          JOFF = NBAST*NBAST*NZ*(IMAT-1)

          CALL GET_ROTVPARA(.TRUE.,
     &                      FO_VORTICITY,
     &                      DENSITY,
     &                      DENSITY_GRADIENT,
     &                      ISYMOP(IMAT)-1,
     &                      GAO,
     &                      GAO1,
     &                      DTRMAT(1+JOFF),
     &                      BUF)

          RS  = (4.1887902D0*DENSITY)**(-DP13)

          IF(RS .GT. 10.0D0) GOTO 999

          DKF = (29.6088132D0*DENSITY)**DP13

          CHI_RATIO = (D1 + 0.028D0*RS)*DEXP(-0.042D0*RS)

          G = D2*WEIGHT*(DKF/236.870506D0)*(CHI_RATIO - D1)

          CALL DSCAL(3,G,FO_VORTICITY,1)

          FO_VORTICITY_NORM = DSQRT(DDOT(3,FO_VORTICITY,1,
     &                                     FO_VORTICITY,1))

          IF(FO_VORTICITY_NORM .LT. DFTHRL) GOTO 999

          CALL GET_ROTVPARA(.FALSE.,
     &                      FO_VORTICITY,
     &                      DENSITY,
     &                      DENSITY_GRADIENT,
     &                      ISYMOP(IMAT)-1,
     &                      GAO,
     &                      GAO1,
     &                      EXCMAT(1+JOFF),
     &                      BUF)

 999      CONTINUE

        ENDIF
      ENDDO

      RETURN
      END
