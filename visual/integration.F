C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL_INTEGRATION */
      SUBROUTINE VISUAL_INTEGRATION(DMAT_AO,
     &                              GAO,
     &                              GAO1,
     &                              NCNT,
     &                              BUF,
     &                              WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "visual.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          NCNT(*),
     &          BUF(*),
     &          WORK(LWORK),
     &          QUANTITY_INTEGRAL(0:3),
     &          QUANTITY_INTEGRAL_LAST(0:3)

#include "memint.h"

      IF(LBIT(LBIT_VISUAL,3)) THEN !3D

        CALL AROUND('3D numerical integration')

        LUNIT = 44

        CALL DFTGRD(DMAT_AO,WORK(KFREE),LFREE,1)

        REWIND LUNIT

 1      CONTINUE
        READ(LUNIT) NR_POINTS

        IF(NR_POINTS .LE. 0) GOTO 2

        CALL MEMGET('REAL',KX,NR_POINTS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KY,NR_POINTS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KZ,NR_POINTS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KW,NR_POINTS,WORK,KFREE,LFREE)

        CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),
     &              WORK(KW),
     &              LUNIT,NR_POINTS)

        CALL DZERO(QUANTITY_INTEGRAL,4)
        CALL INTEGRATE_QUANTITY(DMAT_AO,
     &                          GAO,
     &                          GAO1,
     &                          WORK(KX),WORK(KY),WORK(KZ),
     &                          NCNT,
     &                          BUF,
     &                          QUANTITY_INTEGRAL,
     &                          WORK(KW),
     &                          NR_POINTS,
     &                          WORK,KFREE,LFREE)

        CALL MEMREL('VISUAL_INTEGRATION',WORK,KX,KX,KFREE,LFREE)

        GOTO 1
 2      CONTINUE

        CLOSE(LUNIT,STATUS = 'DELETE')

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A12,A25,2A20)') 'scalar',
     &                                'x-component',
     &                                'y-component',
     &                                'z-component'
        WRITE(LUPRI,'(A)')

        WRITE(LUPRI,'(4E20.10)') QUANTITY_INTEGRAL(0),
     &                           QUANTITY_INTEGRAL(1),
     &                           QUANTITY_INTEGRAL(2),
     &                           QUANTITY_INTEGRAL(3)

        WRITE(LUPRI,'(A)')
        CALL PRSYMB(LUPRI,'-',100,0)

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A)')

      ENDIF !3D

      IF(LBIT(LBIT_VISUAL,2)) THEN !2D

        CALL AROUND('2D Gauss-Lobatto numerical integration')

        DO IORIENTATION = 2,2

          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A12,A3)') 'orientation:',
     &                            ORIENTATION(IORIENTATION)
          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A12,A25,2A20,A14,2A13)') 'scalar',
     &                                           'x-component',
     &                                           'y-component',
     &                                           'z-component',
     &                                           'pieces',
     &                                           'order',
     &                                           'error'
          WRITE(LUPRI,'(A)')

          NR_PIECES = 5
          N         = 10
          DIMI      = 10.0D0
          DIMJ      = 10.0D0
          LUNIT     = 45

          CALL DZERO(QUANTITY_INTEGRAL_LAST,4)

          DO WHILE(.TRUE.)

            CALL LOBATTO_2D_GRID(ORIENTATION(IORIENTATION),
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)

            REWIND LUNIT
            READ(LUNIT) NR_POINTS

            CALL MEMGET('REAL',KX,NR_POINTS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KY,NR_POINTS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KZ,NR_POINTS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KW,NR_POINTS,WORK,KFREE,LFREE)

            CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),
     &                  WORK(KW),
     &                  LUNIT,NR_POINTS)

            CALL DZERO(QUANTITY_INTEGRAL,4)
            CALL INTEGRATE_QUANTITY(DMAT_AO,
     &                              GAO,
     &                              GAO1,
     &                              WORK(KX),WORK(KY),WORK(KZ),
     &                              NCNT,
     &                              BUF,
     &                              QUANTITY_INTEGRAL,
     &                              WORK(KW),
     &                              NR_POINTS,
     &                              WORK,KFREE,LFREE)

            CALL MEMREL('VISUAL_INTEGRATION',WORK,KX,KX,KFREE,LFREE)

            CLOSE(LUNIT,STATUS = 'DELETE')

            ERROR = D0
            DO IXYZ = 0,3
              SOMETHING = QUANTITY_INTEGRAL(IXYZ)
     &                  - QUANTITY_INTEGRAL_LAST(IXYZ)
              ERROR = ERROR + DABS(SOMETHING)
            ENDDO

            WRITE(LUPRI,'(4E20.10,I8,I14,E23.10)') QUANTITY_INTEGRAL(0),
     &                                             QUANTITY_INTEGRAL(1),
     &                                             QUANTITY_INTEGRAL(2),
     &                                             QUANTITY_INTEGRAL(3),
     &                                             NR_PIECES, N, ERROR

            IF(ERROR .LT. 1.0D-7) GOTO 3

            N = N + 10

            CALL DZERO(QUANTITY_INTEGRAL_LAST,4)
            CALL DAXPY(4,D1,QUANTITY_INTEGRAL,1,
     &                      QUANTITY_INTEGRAL_LAST,1)
          ENDDO
 3        CONTINUE

          WRITE(LUPRI,'(A)')
          CALL PRSYMB(LUPRI,'-',100,0)

          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A)')

        ENDDO
      ENDIF !2D

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* INTEGRATE_QUANTITY */
      SUBROUTINE INTEGRATE_QUANTITY(DMAT_AO,
     &                              GAO,
     &                              GAO1,
     &                              PX,PY,PZ,
     &                              NCNT,
     &                              BUF,
     &                              QUANTITY_INTEGRAL,
     &                              WEIGHT,
     &                              NR_POINTS,
     &                              WORK,KFREE,LFREE)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "visual.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          PX(*),
     &          PY(*),
     &          PZ(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY_INTEGRAL(0:3),
     &          WEIGHT(*),
     &          WORK(*),
     &          QUANTITY(0:3)

      DO IPOINT = 1,NR_POINTS
        CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                             GAO,
     &                             GAO1,
     &                             PX(IPOINT),
     &                             PY(IPOINT),
     &                             PZ(IPOINT),
     &                             NCNT,
     &                             BUF,
     &                             QUANTITY,
     &                             WORK,KFREE,LFREE)

        CALL DAXPY(4,WEIGHT(IPOINT),QUANTITY,1,
     &                              QUANTITY_INTEGRAL,1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* LEGENDRE_POLYNOMIAL */
      SUBROUTINE LEGENDRE_POLYNOMIAL(P,N,X)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
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
C     /* LEGENDRE_POLYNOMIAL_DERIVATIVES */
      SUBROUTINE LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      DIMENSION P(0:NMAX),
     &          DP(NMAX)

      DO I = 1,N
        DP(I) = (I*X*P(I) - I*P(I-1))/(X*X - D1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* LOBATTO_1D_GRID */
      SUBROUTINE LOBATTO_1D_GRID(N,X,W)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
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
C     /* LEGENDRE_1D_GRID */
      SUBROUTINE LEGENDRE_1D_GRID(N,X,W)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
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
C     /* LOBATTO_2D_GRID */
      SUBROUTINE LOBATTO_2D_GRID(ORIENTATION,
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
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
              ENDIF

              IF(ORIENTATION .EQ. 'xz') THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PY(IPOINT) = D0
              ENDIF

              IF(ORIENTATION .EQ. 'yz') THEN
                PY(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PX(IPOINT) = D0
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
