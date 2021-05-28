!dirac_copyright_start
!      Copyright (c) 2021 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
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
