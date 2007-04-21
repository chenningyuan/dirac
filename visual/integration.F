C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck INTEGRATE_QUANTITY */
      SUBROUTINE INTEGRATE_QUANTITY(QUANTITY_INTEGRAL,
     &                              DMAT_AO,
     &                              GAO,GAO1,
     &                              PX,PY,PZ,
     &                              WEIGHT,
     &                              NR_POINTS,
     &                              NCNT,BUF,
     &                              WORK,KFREE,LFREE)
C=======================================================================
C     integrate
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "visual.h"

      LOGICAL   LBIT
      DIMENSION QUANTITY_INTEGRAL(NR_QUANTITIES,0:3),
     &          QUANTITY(NR_QUANTITIES,0:3),
     &          DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          PX(*),
     &          PY(*),
     &          PZ(*),
     &          WEIGHT(*),
     &          NCNT(*),
     &          BUF(*),
     &          WORK(*)

      DO IQUANTITY = 1,NR_QUANTITIES
        QUANTITY_IS_A_VECTOR(IQUANTITY) = .FALSE.
      ENDDO

      QUANTITY_IS_A_VECTOR(2)  = .TRUE.
      QUANTITY_IS_A_VECTOR(3)  = .TRUE.
      QUANTITY_IS_A_VECTOR(5)  = .TRUE.
      QUANTITY_IS_A_VECTOR(8)  = .TRUE.
      QUANTITY_IS_A_VECTOR(9)  = .TRUE.
      QUANTITY_IS_A_VECTOR(20) = .TRUE.
      QUANTITY_IS_A_VECTOR(24) = .TRUE.
      QUANTITY_IS_A_VECTOR(26) = .TRUE.

      CALL DZERO(QUANTITY_INTEGRAL,NR_QUANTITIES*4)

      DO IPOINT = 1,NR_POINTS
        CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                             GAO,GAO1,
     &                             PX(IPOINT),
     &                             PY(IPOINT),
     &                             PZ(IPOINT),
     &                             NCNT,BUF,
     &                             QUANTITY,
     &                             WORK,KFREE,LFREE)

        DO IQUANTITY = 1,NR_QUANTITIES
          IF(LBIT(LBIT_VISUAL_QUANTITY,IQUANTITY)) THEN
            QUANTITY_INTEGRAL(IQUANTITY,0)
     &                            = QUANTITY_INTEGRAL(IQUANTITY,0)
     &                            + QUANTITY(IQUANTITY,0)*WEIGHT(IPOINT)

            IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
              DO IXYZ = 1,3
                QUANTITY_INTEGRAL(IQUANTITY,IXYZ)
     &                         = QUANTITY_INTEGRAL(IQUANTITY,IXYZ)
     &                         + QUANTITY(IQUANTITY,IXYZ)*WEIGHT(IPOINT)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_POLYNOMIAL */
      SUBROUTINE LEGENDRE_POLYNOMIAL(P,N,X)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      DIMENSION P(0:NMAX)

      P(0) = 1
      P(1) = X

      DO I = 1,N-1
        P(I+1) = (2*I + 1)*X*P(I) - I*P(I-1)
        P(I+1) = P(I+1)/(I+1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_POLYNOMIAL_DERIVATIVES */
      SUBROUTINE LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      DIMENSION P(0:NMAX),
     &          DP(NMAX)

      DO I = 1,N
        DP(I) = (I*X*P(I) - I*P(I-1))/(X*X - D1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LOBATTO_1D_GRID */
      SUBROUTINE LOBATTO_1D_GRID(N,X,W)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(NMAX),
     &          W(NMAX),
     &          P(0:NMAX),
     &          DP(NMAX)

      X_OLD = -D1 + THRESHOLD
      CALL LEGENDRE_POLYNOMIAL(P,N,X_OLD)
      CALL LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X_OLD)

      DO I = 2,N-1
        STEP  = 1.0D-4
        DO WHILE(DABS(STEP) .GT. THRESHOLD)
          Y_OLD = DP(N-1)
          X_NEW = X_OLD + STEP

          CALL LEGENDRE_POLYNOMIAL(P,N,X_NEW)
          CALL LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X_NEW)
          Y_NEW = DP(N-1)

          IF((Y_OLD/Y_NEW) .LT. D0) THEN
            STEP = -STEP/10.0D0
          ENDIF

          X_OLD = X_NEW
        ENDDO
        X(I) = X_NEW

        W(I) = N*(N - 1)*P(N-1)*P(N-1)
        W(I) = D2/W(I)
      ENDDO

!     end points

      X(1) = -D1
      X(N) =  D1

      W(1) = D2/(N*(N - D1))
      W(N) = W(1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_1D_GRID */
      SUBROUTINE LEGENDRE_1D_GRID(N,X,W)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(NMAX),
     &          W(NMAX),
     &          P(0:NMAX),
     &          DP(NMAX)

      X_OLD = -D1
      CALL LEGENDRE_POLYNOMIAL(P,N,X_OLD)

      DO I = 1,N
        STEP  = 1.0D-4
        DO WHILE(DABS(STEP) .GT. THRESHOLD)
          Y_OLD = P(N)
          X_NEW = X_OLD + STEP

          CALL LEGENDRE_POLYNOMIAL(P,N,X_NEW)
          Y_NEW = P(N)

          IF((Y_OLD/Y_NEW) .LT. D0) THEN
            STEP = -STEP/10.0D0
          ENDIF

          X_OLD = X_NEW
        ENDDO
        X(I) = X_NEW

        CALL LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X(I))
        W(I) = (D1 - X(I)*X(I))*DP(N)*DP(N)
        W(I) = D2/W(I)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LOBATTO_2D_GRID */
      SUBROUTINE LOBATTO_2D_GRID(ORIENTATION,
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      CHARACTER*30 FILE_NAME
      CHARACTER*2  ORIENTATION
      DIMENSION PX(NMAX*NMAX),
     &          PY(NMAX*NMAX),
     &          PZ(NMAX*NMAX),
     &          WEIGHT(NMAX*NMAX),
     &          X(NMAX),
     &          W(NMAX)

      WRITE(FILE_NAME,'(A5,I2)') 'GRID.', LUNIT

      OPEN(LUNIT,
     &     FILE = FILE_NAME,
     &     STATUS = 'UNKNOWN',
     &     ACCESS = 'SEQUENTIAL',
     &     FORM = 'UNFORMATTED')
      REWIND LUNIT

      CALL LOBATTO_1D_GRID(N,X,W)

      IPOINT = 0
      DO IPIECE = 0,NR_PIECES-1
        OFFI = IPIECE*DIMI/NR_PIECES
        DO JPIECE = 0,NR_PIECES-1
          OFFJ = JPIECE*DIMJ/NR_PIECES
!         OFFJ = OFFJ - DIMJ/D2
          DO I = 1,N
            DO J = 1,N
              IPOINT = IPOINT + 1

              IF(ORIENTATION .EQ. 'xy') THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PY(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PZ(IPOINT) = D0
              ELSEIF(ORIENTATION .EQ. 'xz') THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PY(IPOINT) = D0
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
              ELSEIF(ORIENTATION .EQ. 'yz') THEN
                PX(IPOINT) = D0
                PY(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
              ENDIF

              WEIGHT(IPOINT) = W(I)*W(J)*DIMI*DIMJ
              WEIGHT(IPOINT) = WEIGHT(IPOINT)/(D2*D2)
              WEIGHT(IPOINT) = WEIGHT(IPOINT)/(NR_PIECES*NR_PIECES)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NR_POINTS = IPOINT

      CALL WRTQUA(PX,PY,PZ,WEIGHT,LUNIT,NR_POINTS)

      RETURN
      END
