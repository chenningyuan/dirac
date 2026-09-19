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
C
Cjth - copied from dalton/sirius/sirnexkap.F
C
C  /* Deck ortvec */
      SUBROUTINE ORTVEC(NOLD,NNEW,NVDIM,THRLDP,VEC,LINDEP)
C
C 15-Mar-1985 hjaaj
C l.r. 4-May-1994 hjaaj (only elim. new vector if norm < THRLDP)
C
!      use priunit_mod, only: LUPRI
      use kinds_mod, only: dp
#include "implicit.h"
      DIMENSION VEC(NVDIM,*)
      LOGICAL   LINDEP(*)
C
      REAL(dp), PARAMETER :: THRRND=1.D-2, THRTT=1.D-4, D1=1.0D0
C
C     Used from common blocks:
C     PRIUNIT: LUPRI
C
      INTEGER :: J,K
#include "priunit.h"
C
      CALL QENTER('ORTVEC')
C
      IF (NNEW.LE.0) GO TO 9999
      TLINDP = SQRT(THRLDP)
C
C     Normalize NNEW new vectors VEC(*,NOLD+1) TO VEC(*,NOLD+NNEW)
C
      IVEC = NOLD
      DO 200 INEW = 1,NNEW
         IVEC = IVEC + 1
         TT   = DNORM2(NVDIM,VEC(1,IVEC),1)
         IF (TT.LE.THRLDP) THEN
            LINDEP(INEW) = .TRUE.
            WRITE (LUPRI,8100) INEW,TT
         ELSE
            LINDEP(INEW) = .FALSE.
            IF (TT.LT.THRTT) THEN
               CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
               TT = DNORM2(NVDIM,VEC(1,IVEC),1)
            END IF
            CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
         END IF
  200 CONTINUE
C
C
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
C
C     Orthogonalize new vectors against previous vectors
C
      DO 2000 K=1,NOLD
         DO 1900 J = NOLD+1,NOLD+NNEW
            TT = - DDOT(NVDIM,VEC(1,K),1,VEC(1,J),1)
            CALL DAXPY(NVDIM,TT,VEC(1,K),1,VEC(1,J),1)
 1900    CONTINUE
 2000 CONTINUE
C
C     Orthogonalize new vectors against each other
C     and normalization.
C
      DO 2400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
C           ... orthogonalize using prev. vectors are normalized
            IVEC = NOLD + INEW
            DO 2300 JNEW = 1,(INEW-1)
               IF (.NOT.LINDEP(JNEW)) THEN
                  JVEC = NOLD + JNEW
                  TT = - DDOT(NVDIM,VEC(1,JVEC),1,VEC(1,IVEC),1)
                  CALL DAXPY(NVDIM,TT,VEC(1,JVEC),1,VEC(1,IVEC),1)
               END IF
 2300       CONTINUE
            TT = DNORM2(NVDIM,VEC(1,IVEC),1)
            IF (TT .LE. THRLDP) THEN
               LINDEP(INEW) = .TRUE.
               WRITE (LUPRI,8100) INEW,TT
            ELSE
               IF (TT .LT. THRRND) IROUND = IROUND+1
               IF (TT .LT. THRTT) THEN
                  CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
                  TT = DNORM2(NVDIM,VEC(1,IVEC),1)
               END IF
               CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
            END IF
         END IF
 2400 CONTINUE
C
C
      IF (IROUND.GT.0 .AND. ITURN.EQ.1) GO TO 1500
C
C
      JNEW = 0
      DO 4400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
            JNEW = JNEW + 1
            IF (JNEW .LT. INEW) THEN
               CALL DCOPY (NVDIM,VEC(1,NOLD+INEW),1,VEC(1,NOLD+JNEW),1)
            END IF
         END IF
 4400 CONTINUE
      IF (JNEW .LT. NNEW) THEN
         WRITE (LUPRI,8200) NNEW,JNEW
         NNEW = JNEW
      END IF
C
C
C
 8100 FORMAT(/' ORTVEC: New vector no.',I3,
     *        ' is removed because of linear dependence;'
     *       /' norm of vector after Gram-Schmidt''s orthogonalization',
     *        1P,D9.2)
 8200 FORMAT(/' (ORTVEC) NNEW reduced from',I3,' to',I3)
C
C *** End of subroutine ORTVEC
C
 9999 CALL QEXIT('ORTVEC')
      RETURN
      END

