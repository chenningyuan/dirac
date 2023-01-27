!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ibtabini */
      SUBROUTINE IBTABINI(IBTAB)
C***********************************************************************
C
C     Initialize table with IBTAB(I) = NBITS(I)
C
C     Input : none
C
C     Output:
C       IBTAB
C
C     Written by J. Thyssen - Feb 27 2001
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IBTAB(0:255)
C
#include "ibtfun.h"
C
      DO I = 0, 255
         ID = I
         NB = 0
         DO J = 1, 8
            NB = NB + IBTAND(1,ID)
            ID = ID / 2
         END DO
         IBTAB(I) = NB
      END DO
C
      RETURN
      END
