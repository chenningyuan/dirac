!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl

module visual_integrate

   use dirac_interface
   use matrix_defop
   use memory_allocator
   use visual_cfg
   use visual_in_point

  public visual_integration

  save

  private

contains

      subroutine visual_integration(nr_dmat, D, ao, buffer, nr_dmat_con, D_con)

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      integer, intent(in) :: nr_dmat
      type(matrix)        :: D(nr_dmat)
      real(8), intent(in) :: ao(*)
      real(8), intent(in) :: buffer(*)

      integer, intent(in) :: nr_dmat_con
      type(matrix), optional :: D_con(nr_dmat_con)

      dimension QUANTITY_INTEGRAL(0:3), &
                QUANTITY_INTEGRAL_LAST(0:3)

      real(8), allocatable :: x(:)
      real(8), allocatable :: y(:)
      real(8), allocatable :: z(:)
      real(8), allocatable :: w(:)
      character(2) :: grid_orientation(3) = (/'xy', 'xz', 'yz'/)

      CALL DZERO(QUANTITY_INTEGRAL,4)
      if (visual_cfg_3d) then

        CALL AROUND('3D numerical integration')

        LUNIT = 44

        call dftgrd(D(1)%elms, 1)


        REWIND LUNIT

 1      CONTINUE
        READ(LUNIT) NR_POINTS

        IF(NR_POINTS .LE. 0) GOTO 2

        call alloc(x, nr_points)
        call alloc(y, nr_points)
        call alloc(z, nr_points)
        call alloc(w, nr_points)

        call reaqua(x, y, z, w, lunit, nr_points)


          CALL INTEGRATE_QUANTITY(nr_dmat, D, &
     &                            AO, &
     &                            x, y, z, &
     &                            buffer, &
     &                            QUANTITY_INTEGRAL, &
     &                            w, &
     &                            NR_POINTS, &
     &                            0)


        call dealloc(x)
        call dealloc(y)
        call dealloc(z)
        call dealloc(w)

        GOTO 1
 2      CONTINUE


        CLOSE(LUNIT,STATUS = 'DELETE')

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A12,A25,2A20)') 'scalar', &
     &                                'x-component', &
     &                                'y-component', &
     &                                'z-component'
        WRITE(LUPRI,'(A)')

        WRITE(LUPRI,'(4E20.10)') QUANTITY_INTEGRAL(0), &
     &                           QUANTITY_INTEGRAL(1), &
     &                           QUANTITY_INTEGRAL(2), &
     &                           QUANTITY_INTEGRAL(3)

        WRITE(LUPRI,'(A)')
        CALL PRSYMB(LUPRI,'-',100,0)

        WRITE(LUPRI,'(A)')
        WRITE(LUPRI,'(A)')

      ENDIF !3D

      if (visual_cfg_2d) then

        CALL AROUND('2D Gauss-Lobatto numerical integration')

        DO IORIENTATION = visual_cfg_orientation,visual_cfg_orientation
!radovan: this means that i want xz plane
!         1 = xy
!         2 = xz
!         3 = yz
!         i avoided whole loop because with the aromatic compounds
!         one of the planes cut through an atom and it's difficult to 
!         converge then

          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A12,A3)') 'orientation:', &
     &                   grid_orientation(IORIENTATION)
          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A)')
          WRITE(LUPRI,'(A12,A25,2A20,A14,2A13)') 'scalar', &
     &                                           'x-component', &
     &                                           'y-component', &
     &                                           'z-component', &
     &                                           'pieces', &
     &                                           'order', &
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

            CALL LOBATTO_2D_GRID(IORIENTATION, &
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)

            REWIND LUNIT
            READ(LUNIT) NR_POINTS

            call alloc(x, nr_points)
            call alloc(y, nr_points)
            call alloc(z, nr_points)
            call alloc(w, nr_points)

            call reaqua(x, y, z, w, lunit, nr_points)

            CALL DZERO(QUANTITY_INTEGRAL,4)

              CALL INTEGRATE_QUANTITY(nr_dmat, D, &
     &                                AO, &
     &                                x, y, z, &
     &                                buffer, &
     &                                QUANTITY_INTEGRAL, &
     &                                w, &
     &                                NR_POINTS, 0)

            call dealloc(x)
            call dealloc(y)
            call dealloc(z)
            call dealloc(w)

            CLOSE(LUNIT,STATUS = 'DELETE')

            ERROR = D0
            DO IXYZ = 0,3
              SOMETHING = QUANTITY_INTEGRAL(IXYZ) &
     &                  - QUANTITY_INTEGRAL_LAST(IXYZ)
              ERROR = ERROR + DABS(SOMETHING)
            ENDDO

            WRITE(LUPRI,'(4E20.10,I8,I14,E23.10)') QUANTITY_INTEGRAL(0), &
     &                                             QUANTITY_INTEGRAL(1), &
     &                                             QUANTITY_INTEGRAL(2), &
     &                                             QUANTITY_INTEGRAL(3), &
     &                                             NR_PIECES, N, ERROR

            IF(ERROR .LT. 1.0D-7) GOTO 3

            N = N + 10

            CALL DZERO(QUANTITY_INTEGRAL_LAST,4)
            CALL DAXPY(4,D1,QUANTITY_INTEGRAL,1, &
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
      END subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* INTEGRATE_QUANTITY */
      SUBROUTINE INTEGRATE_QUANTITY(nr_dmat, D, &
     &                              AO, &
     &                              PX,PY,PZ, &
     &                              buffer, &
     &                              QUANTITY_INTEGRAL, &
     &                              WEIGHT, &
     &                              NR_POINTS, &
     &                              nr_dmat_con, D_con)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

      integer, intent(in) :: nr_dmat, nr_dmat_con
      type(matrix)        :: D(nr_dmat)
      type(matrix), optional :: D_con(nr_dmat_con)
      dimension AO(*), &
     &          PX(*), &
     &          PY(*), &
     &          PZ(*), &
     &          buffer(*), &
     &          QUANTITY_INTEGRAL(0:3), &
     &          WEIGHT(*), &
     &          QUANTITY(0:3)


      DO IPOINT = 1,NR_POINTS


          CALL GET_QUANTITY_IN_POINT(nr_dmat, D, &
     &                               AO, &
     &                               PX(IPOINT), &
     &                               PY(IPOINT), &
     &                               PZ(IPOINT), &
     &                               buffer, &
     &                               QUANTITY)


        CALL DAXPY(4,WEIGHT(IPOINT),QUANTITY,1, &
     &                              QUANTITY_INTEGRAL,1)
      ENDDO


      RETURN
      END subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* LEGENDRE_POLYNOMIAL */
      SUBROUTINE LEGENDRE_POLYNOMIAL(P,N,X)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

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
      END subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* LEGENDRE_POLYNOMIAL_DERIVATIVES */
      SUBROUTINE LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      DIMENSION P(0:NMAX), &
     &          DP(NMAX)

      DO I = 1,N
        DP(I) = (I*X*P(I) - I*P(I-1))/(X*X - D1)
      ENDDO

      RETURN
      END subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* LOBATTO_1D_GRID */
      SUBROUTINE LOBATTO_1D_GRID(N,X,W)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(NMAX), &
     &          W(NMAX), &
     &          P(0:NMAX), &
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
      END subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* LEGENDRE_1D_GRID */
      SUBROUTINE LEGENDRE_1D_GRID(N,X,W)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(NMAX), &
     &          W(NMAX), &
     &          P(0:NMAX), &
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
      END subroutine
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* LOBATTO_2D_GRID */
      SUBROUTINE LOBATTO_2D_GRID(IORIENTATION, &
     &                           NR_PIECES,N,DIMI,DIMJ,LUNIT)
!=======================================================================
!     description in here ...
!-----------------------------------------------------------------------
!                                                           radovan bast
!                                              last revision: april 2007
!=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      CHARACTER*30 FILE_NAME
      INTEGER      IORIENTATION
      DIMENSION X(NMAX), W(NMAX)
      real(8), allocatable :: px(:), py(:), pz(:), weight(:)

      call alloc(px,     nmax*nmax)
      call alloc(py,     nmax*nmax)
      call alloc(pz,     nmax*nmax)
      call alloc(weight, nmax*nmax)

      WRITE(FILE_NAME,'(A5,I2)') 'GRID.', LUNIT

      OPEN(LUNIT, &
     &     FILE = FILE_NAME, &
     &     STATUS = 'UNKNOWN', &
     &     ACCESS = 'SEQUENTIAL', &
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
! 'xy'
              IF(IORIENTATION .EQ. 1) THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PY(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PZ(IPOINT) = visual_cfg_origin_2d(3)
              ENDIF
! 'xz'
              IF(IORIENTATION .EQ. 2) THEN
                PX(IPOINT) = OFFI + (X(I)+D1)*DIMI/(D2*NR_PIECES)
                PZ(IPOINT) = OFFJ + (X(J)+D1)*DIMJ/(D2*NR_PIECES)
                PY(IPOINT) = visual_cfg_origin_2d(2)
              ENDIF
! 'yz'
              IF(IORIENTATION .EQ. 3) THEN
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
      END subroutine

end module
