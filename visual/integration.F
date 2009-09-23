C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL_INTEGRATION */
      SUBROUTINE VISUAL_INTEGRATION(DMAT_AO,
     &                              GAO,
     &                              GAO1,
     &                              GAO2,
     &                              NCNT,
     &                              BUF)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

       use dft_main
       use visual_cfg
       use memory_allocator

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          GAO2(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY_INTEGRAL(0:3),
     &          QUANTITY_INTEGRAL_LAST(0:3)

      real(8), allocatable :: x(:)
      real(8), allocatable :: y(:)
      real(8), allocatable :: z(:)
      real(8), allocatable :: w(:)

      CALL DZERO(QUANTITY_INTEGRAL,4)
      if (visual_cfg_3d) then

        CALL AROUND('3D numerical integration')

        LUNIT = 44

        call dftgrd(dmat_ao, 1)

        REWIND LUNIT

 1      CONTINUE
        READ(LUNIT) NR_POINTS

        IF(NR_POINTS .LE. 0) GOTO 2

        call alloc(x, nr_points)
        call alloc(y, nr_points)
        call alloc(z, nr_points)
        call alloc(w, nr_points)

        call reaqua(x, y, z, w, lunit, nr_points)

        CALL INTEGRATE_QUANTITY(DMAT_AO,
     &                          GAO,
     &                          GAO1,
     &                          GAO2,
     &                          x, y, z,
     &                          NCNT,
     &                          BUF,
     &                          QUANTITY_INTEGRAL,
     &                          w,
     &                          NR_POINTS)

        call dealloc(x)
        call dealloc(y)
        call dealloc(z)
        call dealloc(w)

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

      if (visual_cfg_2d) then

        CALL AROUND('2D Gauss-Lobatto numerical integration')

        DO IORIENTATION = 2,2
!radovan: this means that i want xz plane
!         1 = xy
!         2 = xz
!         3 = yz
!         i avoided whole loop because with the aromatic compounds
!         one of the planes cut through an atom and it's difficult to 
!         converge then

          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A12,A3)') 'orientation:',
     &                   visual_cfg_orientation(IORIENTATION)
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
!radovan: this means i integrate from 0,0 to 10,10 a.u.
!         due to planar symmetry i integrated only one half!!!
          DIMI      = 10.0D0
          DIMJ      = 10.0D0
          LUNIT     = 45

          CALL DZERO(QUANTITY_INTEGRAL_LAST,4)

          DO WHILE(.TRUE.)

            CALL LOBATTO_2D_GRID(visual_cfg_orientation(IORIENTATION),
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)

            REWIND LUNIT
            READ(LUNIT) NR_POINTS

            call alloc(x, nr_points)
            call alloc(y, nr_points)
            call alloc(z, nr_points)
            call alloc(w, nr_points)

            call reaqua(x, y, z, w, lunit, nr_points)

            CALL DZERO(QUANTITY_INTEGRAL,4)
            CALL INTEGRATE_QUANTITY(DMAT_AO,
     &                              GAO,
     &                              GAO1,
     &                              GAO2,
     &                              x, y, z,
     &                              NCNT,
     &                              BUF,
     &                              QUANTITY_INTEGRAL,
     &                              w,
     &                              NR_POINTS)

            call dealloc(x)
            call dealloc(y)
            call dealloc(z)
            call dealloc(w)

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
     &                              GAO2,
     &                              PX,PY,PZ,
     &                              NCNT,
     &                              BUF,
     &                              QUANTITY_INTEGRAL,
     &                              WEIGHT,
     &                              NR_POINTS)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

      use visual
#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          GAO2(*),
     &          PX(*),
     &          PY(*),
     &          PZ(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY_INTEGRAL(0:3),
     &          WEIGHT(*),
     &          QUANTITY(0:3)

      DO IPOINT = 1,NR_POINTS
        CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                             GAO,
     &                             GAO1,
     &                             GAO2,
     &                             PX(IPOINT),
     &                             PY(IPOINT),
     &                             PZ(IPOINT),
     &                             NCNT,
     &                             BUF,
     &                             QUANTITY)

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
      SUBROUTINE LOBATTO_2D_GRID(GRID_ORIENTATION,
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

      use memory_allocator
      use visual_cfg
#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      CHARACTER*30 FILE_NAME
      CHARACTER*2  GRID_ORIENTATION
      DIMENSION X(NMAX), W(NMAX)
      real(8), allocatable :: px(:), py(:), pz(:), weight(:)

      call alloc(px,     nmax*nmax)
      call alloc(py,     nmax*nmax)
      call alloc(pz,     nmax*nmax)
      call alloc(weight, nmax*nmax)

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

              IF(GRID_ORIENTATION .EQ. 'xy') THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PY(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PZ(IPOINT) = visual_cfg_origin_2d(3)
              ENDIF

              IF(GRID_ORIENTATION .EQ. 'xz') THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PY(IPOINT) = visual_cfg_origin_2d(2)
              ENDIF

              IF(GRID_ORIENTATION .EQ. 'yz') THEN
                PY(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PX(IPOINT) = visual_cfg_origin_2d(1)
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

      call dealloc(px)
      call dealloc(py)
      call dealloc(pz)
      call dealloc(weight)

      RETURN
      END
