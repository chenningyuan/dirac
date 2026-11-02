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
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, J
!#include "priunit.h"
      REAL(dp), PARAMETER :: D0=0.0D0, D2 = 2.00D00
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
      INTEGER :: I, J
!#include "priunit.h"
C
      DIMENSION IBTAB(0:255)
C
C
      DO I = 0, 255
         ID = I
         NB = 0
         DO J = 1, 8
            NB = NB + IAND(1,ID)
            ID = ID / 2
         END DO
         IBTAB(I) = NB
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ini_bit_tab */
      SUBROUTINE ini_cb_bittab()
C***********************************************************************
C
C     Initialize table with bit_tab(I) = nbits(I)
C                           bit_dif(I,J)
C     (I and J integer*1)
C
C     Input : none
C
C     Output: IBTAB and IBTDIF in cb_bittab.h
C
C     Written by H. J. Aa. Jensen - July 2016
C     Partly based on IBTABINI by J. Thyssen 2001
C     Last revision :
C
C***********************************************************************
C
!      use cb_bitdef_mod, only: BIT_CNT, BIT_DIF
#include "cb_bitdef.h"
! from cb_bitdef.h:
!     integer bit_cnt(0:255), bit_dif(0:255,0:255)
!

!T    print '(/A)', 'i, bits(i), bit_cnt(i) :'
      DO I = 0, 255
         nbits = 0
         DO J = 0, 7
            IF (BTEST(I,J)) nbits = nbits + 1
         END DO
         bit_cnt(I) = nbits
!T       print '(i8, 5x, B8.8, I8)', i,i,bit_cnt(i)
      END DO

!T    print '(/A)', 'i,j, bits(i),bits(j),diff(i,j), bit_dif(i,j) :'
      DO J = 0, 255
         DO I = 0, 255
            ID = IEOR(I, J)
            bit_dif(I,J) = bit_cnt(ID)
!T          print '(2i8, 2x, B8.8,2x,B8.8,2x,B8.8,I8)',
!T   &      i,j,i,j,id,bit_dif(i,j)
         END DO
      END DO
C
      RETURN
      END
C -- end of gpjth.F ---
