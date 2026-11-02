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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CUBIC_SPLINE_SETUP */
      SUBROUTINE CUBIC_SPLINE_SETUP(N,XA,YA,D2YA)
C=======================================================================
C     you provide:
C       N
C       XA = (x_1, x_2, ..., x_N)
C       YA = (y_1, y_2, ..., y_N)
C
C     routine gives back:
C       D2YA = (d2y_1, d2y_2, ..., d2y_N)
C       second derivatives at x_1, x_2, ..., x_N
C
C     natural boundary condition:
C       d2y_1 = d2y_N = 0
C
C     this routine is typically called once
C     and CUBIC_SPLINE_INTERPOLATION with the here obtained
C     D2YA can then be called as often as needed
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: december 2006
C=======================================================================

#include "implicit.h"
      INTEGER :: I, K

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 9999)

      DIMENSION XA(NMAX),
     &          YA(NMAX),
     &          D2YA(NMAX),
     &          TEMP(NMAX)

      D2YA(1) = D0
      D2YA(N) = D0

      TEMP(1) = D0

      DO I = 2,N-1
        A = (XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
        B = A*D2YA(I-1) + D2

        D2YA(I) = (A - D1)/B

        TEMP(I) = (YA(I+1) - YA(I))/(XA(I+1) - XA(I))
        TEMP(I) = TEMP(I) - (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
        TEMP(I) = 6.0D0*TEMP(I)/(XA(I+1) - XA(I-1)) - A*TEMP(I-1)
        TEMP(I) = TEMP(I)/B
      ENDDO

      DO K = N-1,1,-1
        D2YA(K) = D2YA(K)*D2YA(K+1) + TEMP(K)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* CUBIC_SPLINE_INTERPOLATION */
      SUBROUTINE CUBIC_SPLINE_INTERPOLATION(N,XA,YA,D2YA,X,Y)
C=======================================================================
C     you provide:
C         N
C         XA = (x_1, x_2, ..., x_N)
C         YA = (y_1, y_2, ..., y_N)
C       D2YA = (d2y_1, d2y_2, ..., d2y_N)
C         X
C
C     routine gives back:
C         Y(X)
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: december 2006
C=======================================================================

#include "implicit.h"
      INTEGER :: K

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 9999)

      DIMENSION XA(NMAX),
     &          YA(NMAX),
     &          D2YA(NMAX)

      IS = 1
      IE = N

 1    CONTINUE
      IF(IE - IS .GT. 1) THEN
          K = (IE + IS)/2
        IF(XA(K) .GT. X) THEN
          IE = K
        ELSE
          IS = K
        ENDIF
        GOTO 1
      ENDIF

      H = XA(IE) - XA(IS)
      A = (XA(IE) - X)/H
      B = (X - XA(IS))/H
      Y = A*YA(IS) + B*YA(IE)
     &  + ((A**3 - A)*D2YA(IS) + (B**3 - B)*D2YA(IE))*(H**2)/6.0D0

      RETURN
      END
