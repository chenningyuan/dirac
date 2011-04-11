!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

      SUBROUTINE VRG_LR(DTRMAT,
     &                  EXCMAT,
     &                  NXCMAT,
     &                  ISYMOP,
     &                  IHRMOP,
     &                  GAO,
     &                  GAO1,
     &                  DENSITY,
     &                  DENSITY_GRADIENT_NORM,
     &                  DENSITY_GRADIENT,
     &                  WEIGHT,
     &                  IS_GGA,
     &                  BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

      use dft_cfg

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

          IF(FO_VORTICITY_NORM .LT. dft_cfg_DFTHRL) GOTO 999

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
