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

          CALL GET_V(.TRUE.,
     &               FO_VORTICITY,
     &               DENSITY,
     &               DENSITY_GRADIENT,
     &               ISYMOP(IMAT)-1,
     &               GAO,
     &               GAO1,
     &               DTRMAT(1+JOFF),
     &               BUF)

!        IC = 1
!        DO IXYZ = 1,3

!          FO_VORTICITY(IXYZ) = D0

!          JXYZ = IXYZ + 1
!          KXYZ = JXYZ + 1
!          IF(JXYZ .GT. 3) JXYZ = JXYZ - 3
!          IF(KXYZ .GT. 3) KXYZ = KXYZ - 3

!          S = D1/DENSITY

!          CALL BLUBDEN(A,IC,ISYMOP(IMAT)-1,
!    &                  DTRMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(JXYZ-1)),
!    &                  GAO1(1+NBAST*(KXYZ-1)),
!    &                  BUF)
!          FO_VORTICITY(IXYZ) = FO_VORTICITY(IXYZ) + S*A

!          CALL BLUBDEN(A,IC,ISYMOP(IMAT)-1,
!    &                  DTRMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(KXYZ-1)),
!    &                  GAO1(1+NBAST*(JXYZ-1)),
!    &                  BUF)
!          FO_VORTICITY(IXYZ) = FO_VORTICITY(IXYZ) - S*A

!          SS = D1/(D2*DENSITY*DENSITY)
!          S = D1

!          CALL BLUBDEN(A,IC,ISYMOP(IMAT)-1,
!    &                  DTRMAT(1+JOFF),
!    &                  GAO,
!    &                  GAO1(1+NBAST*(KXYZ-1)),
!    &                  BUF)
!          FO_VORTICITY(IXYZ) = FO_VORTICITY(IXYZ) 
!    &                        + (SS)*(S)*A*DENSITY_GRADIENT(JXYZ)
!          CALL BLUBDEN(A,IC,ISYMOP(IMAT)-1,
!    &                  DTRMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(KXYZ-1)),
!    &                  GAO,
!    &                  BUF)
!          FO_VORTICITY(IXYZ) = FO_VORTICITY(IXYZ) 
!    &                        + (SS)*(-S)*A*DENSITY_GRADIENT(JXYZ)
!          CALL BLUBDEN(A,IC,ISYMOP(IMAT)-1,
!    &                  DTRMAT(1+JOFF),
!    &                  GAO,
!    &                  GAO1(1+NBAST*(JXYZ-1)),
!    &                  BUF)
!          FO_VORTICITY(IXYZ) = FO_VORTICITY(IXYZ) 
!    &                        + (-SS)*(S)*A*DENSITY_GRADIENT(KXYZ)
!          CALL BLUBDEN(A,IC,ISYMOP(IMAT)-1,
!    &                  DTRMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(JXYZ-1)),
!    &                  GAO,
!    &                  BUF)
!          FO_VORTICITY(IXYZ) = FO_VORTICITY(IXYZ) 
!    &                        + (-SS)*(-S)*A*DENSITY_GRADIENT(KXYZ)
!        ENDDO





          RS  = (4.1887902D0*DENSITY)**(-DP13)

          IF(RS .GT. 10.0D0) GOTO 999

          DKF = (29.6088132D0*DENSITY)**DP13

          CHI_RATIO = (D1 + 0.028D0*RS)*DEXP(-0.042D0*RS)

          G = D2*WEIGHT*(DKF/236.870506D0)*(CHI_RATIO - D1)

          CALL DSCAL(3,G,FO_VORTICITY,1)

          FO_VORTICITY_NORM = DSQRT(DDOT(3,FO_VORTICITY,1,
     &                                     FO_VORTICITY,1))

          IF(FO_VORTICITY_NORM .LT. DFTHRL) GOTO 999

          CALL GET_V(.FALSE.,
     &               FO_VORTICITY,
     &               DENSITY,
     &               DENSITY_GRADIENT,
     &               ISYMOP(IMAT)-1,
     &               GAO,
     &               GAO1,
     &               EXCMAT(1+JOFF),
     &               BUF)




!        DO IXYZ = 1,3

!          JXYZ = IXYZ + 1
!          KXYZ = JXYZ + 1
!          IF(JXYZ .GT. 3) JXYZ = JXYZ - 3
!          IF(KXYZ .GT. 3) KXYZ = KXYZ - 3

!          S = D1/DENSITY

!          A = S*FO_VORTICITY(IXYZ)
!          CALL BLUBEXC(A,ISYMOP(IMAT)-1,
!    &                  EXCMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(JXYZ-1)),
!    &                  GAO1(1+NBAST*(KXYZ-1)))
!          A = -S*FO_VORTICITY(IXYZ)
!          CALL BLUBEXC(A,ISYMOP(IMAT)-1,
!    &                  EXCMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(KXYZ-1)),
!    &                  GAO1(1+NBAST*(JXYZ-1)))

!          SS = D1/(D2*DENSITY*DENSITY)
!          S = D1



!          A = FO_VORTICITY(IXYZ)*DENSITY_GRADIENT(JXYZ)
!          A = A*(SS)*(S)
!          CALL BLUBEXC(A,ISYMOP(IMAT)-1,
!    &                  EXCMAT(1+JOFF),
!    &                  GAO,
!    &                  GAO1(1+NBAST*(KXYZ-1)))
!          A = FO_VORTICITY(IXYZ)*DENSITY_GRADIENT(JXYZ)
!          A = A*(SS)*(-S)
!          CALL BLUBEXC(A,ISYMOP(IMAT)-1,
!    &                  EXCMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(KXYZ-1)),
!    &                  GAO)


!          A = FO_VORTICITY(IXYZ)*DENSITY_GRADIENT(KXYZ)
!          A = A*(-SS)*(S)
!          CALL BLUBEXC(A,ISYMOP(IMAT)-1,
!    &                  EXCMAT(1+JOFF),
!    &                  GAO,
!    &                  GAO1(1+NBAST*(JXYZ-1)))
!          A = FO_VORTICITY(IXYZ)*DENSITY_GRADIENT(KXYZ)
!          A = A*(-SS)*(-S)
!          CALL BLUBEXC(A,ISYMOP(IMAT)-1,
!    &                  EXCMAT(1+JOFF),
!    &                  GAO1(1+NBAST*(JXYZ-1)),
!    &                  GAO)


!        ENDDO














 999      CONTINUE

        ENDIF
      ENDDO

      RETURN
      END





C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE blubexc(VXC,IREPOP,EXCMAT,GAOL,GAOR)
C***********************************************************************
C
C     Construct Kohn-Sham matrix
C
C     Written by T.Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAOL(NBAST),GAOR(NBAST),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAOL(NTBAS(0)),GAOR(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0))
#endif
C      
#include "ibtfun.h"
C
      DO 10 IRP1 = 0,NBSYM-1
        IRP2 = IBTXOR(IRP1,IREPOP)
        DO 20 IC = 1,MC
          NB1 = NBBAS(IRP1,IC)
          IF(NB1.EQ.0) GOTO 20
          NB2 = NBBAS(IRP2,IC)
          IF(NB2.EQ.0) GOTO 20
          JB1 = IBBAS(IRP1,IC)
          JB2 = IBBAS(IRP2,IC)
          DO J = 1, NBBAS(IRP2,IC)
            GVXC = VXC*GAOR(JB2+J)
            CALL DAXPY(NBBAS(IRP1,IC),GVXC,
     &                 GAOL(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)
          ENDDO
   20   CONTINUE
   10 CONTINUE
C
      END



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getrh1 */
      SUBROUTINE blubden(RHO,IC,IREPOP,DMAT,GAOL,GAOR,BUF)
C***********************************************************************
C
C     Get density in point (PX,PY,PZ).
C     Gets values of AOs in this point from GAO
C
C     Written by T.Saue and T.Helgaker
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAOL(NBAST),DMAT(NBAST,NBAST),BUF(*)
     &          GAOR(NBAST)
#else
      DIMENSION GAOL(NTBAS(0)),DMAT(NTBAS(0),NTBAS(0)),BUF(*),
     &          GAOR(NTBAS(0))
#endif
C
#include "ibtfun.h"
C      
      RHO = D0
      DO 10 IRP1 = 0,NBSYM-1
        NB1 = NBBAS(IRP1,IC)
        IF(NB1.EQ.0) GOTO 10
        JB1 = IBBAS(IRP1,IC)+1
        IRP2 = IBTXOR(IRP1,IREPOP)
        NB2 = NBBAS(IRP2,IC)
        IF(NB2.EQ.0) GOTO 10
        JB2 = IBBAS(IRP2,IC)+1
        CALL DGEMV('N',NBBAS(IRP1,IC),NBBAS(IRP2,IC),D1,
     &             DMAT(JB1,JB2),NTBAS(0),GAOR(JB2),1,D0,BUF,1)
        RHO = RHO + DDOT(NBBAS(IRP1,IC),GAOL(JB1),1,BUF,1)
   10 CONTINUE
C
      END
