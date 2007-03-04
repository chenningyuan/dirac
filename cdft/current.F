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
C     /* GET_CURRENT */
      SUBROUTINE GET_CURRENT(LBIT_CURRENT,
     &                       CURRENT,
     &                       CURRENT_DERVS,
     &                       IREPOP,GAO,GAO1,DMAT,BUF)
C=======================================================================
C
C     routine returns
C
C     CURRENT(1) = j_x    s_x
C     CURRENT(2) = j_y or s_y
C     CURRENT(3) = j_z    s_z
C
C     and/or
C
C     its derivatives
C
C     CURRENT_DERVS(1,1) = derivative of x-component wrt x
C     CURRENT_DERVS(1,2) = derivative of x-component wrt y
C     CURRENT_DERVS(1,3) = derivative of x-component wrt z
C
C     CURRENT_DERVS(2,1) = derivative of y-component wrt x
C     CURRENT_DERVS(2,2) = derivative of y-component wrt y
C     CURRENT_DERVS(2,3) = derivative of y-component wrt z
C
C     CURRENT_DERVS(3,1) = derivative of z-component wrt x
C     CURRENT_DERVS(3,2) = derivative of z-component wrt y
C     CURRENT_DERVS(3,3) = derivative of z-component wrt z
C
C     options
C
C     (LBIT_CURRENT,1) = .TRUE.  current ...
C     (LBIT_CURRENT,1) = .FALSE. spin ...
C     (LBIT_CURRENT,2) = .TRUE.  density
C     (LBIT_CURRENT,2) = .FALSE. no density
C     (LBIT_CURRENT,3) = .TRUE.  derivatives
C     (LBIT_CURRENT,3) = .FALSE. no derivatives
C
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
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "cdft.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
#endif
      LOGICAL   LBIT
      DIMENSION CURRENT(3),
     &          CURRENT_DERVS(3,3),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

#include "ibtfun.h"

      IF(      LBIT(LBIT_CURRENT,1)) THEN !current
!       F  = -CVAL
!       G  = -CVAL
!       dirty fix
        F  = -137.0359998D0
        G  = -137.0359998D0
        IU = 1
        IV = 2
        IX = 2
        IY = 1
        ITRANSFORM = 1
      ENDIF
      IF(.NOT. LBIT(LBIT_CURRENT,1)) THEN !spin
        F  = D1
        G  = D1
        IF(LBIT(LBIT_CDFT,4)) G = -D1 !spin density operator is \beta*\Sigma instead of \Sigma
        IU = 1
        IV = 1
        IX = 2
        IY = 2
        ITRANSFORM = 2
      ENDIF
      IF(LBIT(LBIT_CURRENT,2)) THEN !density
        CALL DZERO(CURRENT,3)
      ENDIF
      IF(LBIT(LBIT_CURRENT,3)) THEN !derivatives
        CALL DZERO(CURRENT_DERVS,9)
      ENDIF

      ISTART = 1
      IF(LBIT(LBIT_CDFT,3)) ISTART = 3 !collinear approximation

      DO IXYZ = ISTART,3
#if defined (VAR_G77) || defined (SYS_SX)
        IOFF = NBAST*NBAST
        JOFF = NBAST
#else
        IOFF = NTBAS(0)*NTBAS(0)
        JOFF = NTBAS(0)
#endif
        IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
        IOFF = IOFF*(IQ-1)

        INTERMEDIATE = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),IREPOP)

        IF(LBIT(LBIT_CURRENT,2)) THEN !density
          CALL GET_CURRENT_SUB(CURRENT(IXYZ),
     &                         F,IU,IV,
     &                         GAO,GAO,
     &                         INTERMEDIATE,DMAT(1+IOFF),BUF)
          CALL GET_CURRENT_SUB(CURRENT(IXYZ),
     &                         G,IX,IY,
     &                         GAO,GAO,
     &                         INTERMEDIATE,DMAT(1+IOFF),BUF)
        ENDIF
        IF(LBIT(LBIT_CURRENT,3)) THEN !derivatives
          DO JXYZ = 1,3
            CALL GET_CURRENT_SUB(CURRENT_DERVS(IXYZ,JXYZ),
     &                           F,IU,IV,
     &                           GAO,
     &                           GAO1(1+JOFF*(JXYZ-1)),
     &                           INTERMEDIATE,DMAT(1+IOFF),BUF)
            CALL GET_CURRENT_SUB(CURRENT_DERVS(IXYZ,JXYZ),
     &                           G,IX,IY,
     &                           GAO,
     &                           GAO1(1+JOFF*(JXYZ-1)),
     &                           INTERMEDIATE,DMAT(1+IOFF),BUF)
            CURRENT_DERVS(IXYZ,JXYZ) =
     &               CURRENT_DERVS(IXYZ,JXYZ) + CURRENT_DERVS(IXYZ,JXYZ)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_CURRENT_SUB */
      SUBROUTINE GET_CURRENT_SUB(CURRENT,
     &                           F,IC1,IC2,
     &                           GAOL,
     &                           GAOR,
     &                           IREPOP,DMAT,BUF)
C=======================================================================
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
#include "dgroup.h"
#include "dcbbas.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION GAOL(NBAST),
     &          GAOR(NBAST),
     &          DMAT(NBAST,NBAST)
#else
      DIMENSION GAOL(NTBAS(0)),
     &          GAOR(NTBAS(0)),
     &          DMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION BUF(*)

#include "ibtfun.h"

      DO IREP1 = 0,NBSYM-1

        N1 = NBBAS(IREP1,IC1)
        IF(N1 .GT. 0) THEN

          I1 = IBBAS(IREP1,IC1)+1
          IREP2 = IBTXOR(IREP1,IREPOP)
          N2 = NBBAS(IREP2,IC2)
          IF(N2 .GT. 0) THEN

            I2 = IBBAS(IREP2,IC2)+1
            CALL DGEMV('N',NBBAS(IREP1,IC1),NBBAS(IREP2,IC2),D1,
     &                 DMAT(I1,I2),NTBAS(0),GAOR(I2),1,D0,BUF,1)
            CURRENT = CURRENT+F*DDOT(NBBAS(IREP1,IC1),GAOL(I1),1,BUF,1)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_CURRENT_DIV */
      SUBROUTINE GET_CURRENT_DIV(CURRENT_DERVS,CURRENT_DIV)
C=======================================================================
C     calculate divergence
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION CURRENT_DERVS(3,3)

      CURRENT_DIV = D0
      DO IXYZ = 1,3
        CURRENT_DIV = CURRENT_DIV + CURRENT_DERVS(IXYZ,IXYZ)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_CURRENT_ROT */
      SUBROUTINE GET_CURRENT_ROT(CURRENT_DERVS,CURRENT_ROT)
C=======================================================================
C     calculate curl
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION CURRENT_DERVS(3,3),
     &          CURRENT_ROT(3)

      CURRENT_ROT(1) = -CURRENT_DERVS(2,3) + CURRENT_DERVS(3,2)
      CURRENT_ROT(2) = -CURRENT_DERVS(3,1) + CURRENT_DERVS(1,3)
      CURRENT_ROT(3) = -CURRENT_DERVS(1,2) + CURRENT_DERVS(2,1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_GAMMA5_DENSITY */
      SUBROUTINE GET_GAMMA5_DENSITY(GAMMA5_DENSITY,
     &                              IREPOP,GAO,DMAT,BUF)
C=======================================================================
C     routine returns
C     \psi^\dag(r) \gamma_5 \psi(r)
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "cdft.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
#endif
      DIMENSION GAO(*),
     &          DMAT(*),
     &          BUF(*)

#include "ibtfun.h"

      F  = D1
      G  = D1
      IU = 1
      IV = 2
      IX = 2
      IY = 1

      GAMMA5_DENSITY = D0

      CALL GET_CURRENT_SUB(GAMMA5_DENSITY,
     &                     F,IU,IV,
     &                     GAO,GAO,
     &                     IREPOP,DMAT,BUF)
      CALL GET_CURRENT_SUB(GAMMA5_DENSITY,
     &                     G,IX,IY,
     &                     GAO,GAO,
     &                     IREPOP,DMAT,BUF)

      RETURN
      END
