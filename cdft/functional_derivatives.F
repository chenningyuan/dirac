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
C     /* JPBE_DERIVATIVES */
      SUBROUTINE JPBE_DERIVATIVES(WEIGHT,
     &                            DENSITY,
     &                            DENSITY_GRADIENT_NORM,
     &                            CURRENT_DENSITY_NORM,
     &                            DEXC_J_J)
C=======================================================================
C     description ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: december 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (D5 = 5.0D0,D6 = 6.0D0,D7 = 7.0D0,D8 = 8.0D0,D9 = 9.0D0)

      PARAMETER (D12 = D1/D2)
      PARAMETER (D13 = D1/D3)
      PARAMETER (D23 = D2/D3)
      PARAMETER (D43 = D4/D3)
      PARAMETER (D83 = D8/D3)

      RHOA = DENSITY/D2
      RHOB = RHOA

      GRADA = DENSITY_GRADIENT_NORM/D2
      GRADB = GRADA

      GRADAB = GRADA*GRADB

      JA = CURRENT_DENSITY_NORM/D2
      JB = JA

#ifdef dont_include_me
!     debug: DEXC_N_N

      F1  = D2**D13
      F2  = GRADA**D2
      F3  = D1/(RHOA**D83)
      F4  = JA**D2
      F5  = 0.27303041196629D0
     &       *(0.10558462307238D0*F4*F3 + 0.01645530784602D0*F2*F3) + D1
      F6  = D1/(RHOA**3.666666666666667D0)
      F7  = -0.28155899485969D0*F4*F6 - 0.04388082092272D0*F2*F6
      F8  = D1/(F5**D2)
      F9  = RHOA**D43
      F10 = D1/(RHOA**4.666666666666667D0)

      DEXC_2000000 = 0.5*(-0.65649668122846D0*F1*(1.804D0 - 0.804D0/F5)
     &                        /(RHOA**D23) - 0.86467093020941D0*F1*F7*F8
     &                              *(RHOA**D13) - 0.32425159882853D0*F1
     &   *(1.032382981152195D0*F4*F10 + 0.16089634338331D0*F2*F10)*F8*F9
     &                     + 0.17706109521776D0*F1*(F7**D2)*F9/(F5**D3))

      X_J_J    = D12*DEXC_2000000
#endif

      F1  = D2**D13
!     F2  = JA**D2
      F3  = D1/(RHOA**D83)
      F4  = 0.27303041196629D0*(0.01645530784602D0*F3
!    &                     *(GRADA**D2) + 0.10558462307238D0*F2*F3) + D1
     &                     *(GRADA**D2)) + D1

!     DEXC_0000020 = 0.00394779406352D0*F1*F2
!    &                                            /((F4**D3)*(RHOA**D4))
      DEXC_0000020 =
     &                      - 0.03423598284293*F1/((F4**D2)*(RHOA**D43))

      X_J_J    = D12*DEXC_0000020

      C_J_J    = D0
      DEXC_J_J = X_J_J + C_J_J
      DEXC_J_J = WEIGHT*DEXC_J_J

      RETURN
      END
