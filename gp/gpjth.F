C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : gpjth.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SUMMMH */
      SUBROUTINE SUMMMH(AMAT,NDIM,NZ,LRA,LCA)
C*****************************************************************************
C
C     Take sum of a general quaternionic matrix AMAT and
C     its Hermitian conjugate; return in AMAT. Note that this gives
C     an Hermitian result !
C
C     Written by J. Thyssen, Dec 1, 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0, D2 = 2.00D00)
      DIMENSION AMAT(LRA,LCA,NZ)
C
C     Take difference
C
      DO J = 1,NDIM
        DO I = 1,(J-1)
          AMAT(I,J,1) = AMAT(I,J,1)+AMAT(J,I,1)
          AMAT(J,I,1) = AMAT(I,J,1)
        ENDDO
        AMAT(J,J,1) = D2*AMAT(J,J,1)
      ENDDO
      DO IZ = 2,NZ
        DO J = 1,NDIM
          DO I = 1,(J-1)
            AMAT(I,J,IZ) = AMAT(I,J,IZ)-AMAT(J,I,IZ)
            AMAT(J,I,IZ) = -AMAT(I,J,IZ)
          ENDDO
          AMAT(J,J,IZ) = D0
        ENDDO
      ENDDO
C
      RETURN
      END
