!dirac_copyright_start
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
!dirac_copyright_end

module visual_integrate

   use dirac_interface
   use matrix_defop
   use memory_allocator
   use visual_cfg
   use visual_in_point
#ifdef MOD_UNRELEASED
   use visual_london
#endif

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
                                  nr_dmat_con, D_con, &
                                  AO, &
                                  x, y, z, &
                                  buffer, &
                                  QUANTITY_INTEGRAL, &
                                  w, &
                                  NR_POINTS)


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


      if (visual_cfg_plane_integration) then

         CALL AROUND('2D Gauss-Lobatto numerical integration')

         WRITE(LUPRI,'(A)')
         WRITE(LUPRI,'(A12,A25,2A20,A14,A13)') 'scalar', &
     &                                         'x-component', &
     &                                         'y-component', &
     &                                         'z-component', &
     &                                         'pieces', &
     &                                         'order'
         WRITE(LUPRI,'(A)')

         LUNIT     = 45

         NR_PIECES = visual_cfg_plane_integration_nr_pieces
         N         = visual_cfg_plane_integration_order
         DIMI      = 10.0D0
         DIMJ      = 10.0D0

         call lobatto_2d_grid(nr_pieces, &
                              n,         &
                          visual_cfg_plane_integration_p_origin, &
                          visual_cfg_plane_integration_p_right, &
                          visual_cfg_plane_integration_p_top, &
                              lunit)

         REWIND LUNIT
         READ(LUNIT) NR_POINTS

         call alloc(x, nr_points)
         call alloc(y, nr_points)
         call alloc(z, nr_points)
         call alloc(w, nr_points)

         call reaqua(x, y, z, w, lunit, nr_points)

         quantity_integral = 0.0d0

         CALL INTEGRATE_QUANTITY(nr_dmat, D, &
                                 nr_dmat_con, D_con, &
                                 AO, &
                                 x, y, z, &
                                 buffer, &
                                 QUANTITY_INTEGRAL, &
                                 w, &
                                 NR_POINTS)

         call dealloc(x)
         call dealloc(y)
         call dealloc(z)
         call dealloc(w)

         CLOSE(LUNIT,STATUS = 'DELETE')

         WRITE(LUPRI,'(4E20.10,I8,I14)') QUANTITY_INTEGRAL(0), &
     &                                   QUANTITY_INTEGRAL(1), &
     &                                   QUANTITY_INTEGRAL(2), &
     &                                   QUANTITY_INTEGRAL(3), &
     &                                   NR_PIECES, N

         WRITE(LUPRI,'(A)')
         CALL PRSYMB(LUPRI,'-',100,0)

         WRITE(LUPRI,'(A)')
         WRITE(LUPRI,'(A)')

      end if

   end subroutine

      SUBROUTINE INTEGRATE_QUANTITY(nr_dmat, D, &
                                    nr_dmat_con, D_con, &
                                    AO, &
                                    PX,PY,PZ, &
                                    buffer, &
                                    QUANTITY_INTEGRAL, &
                                    WEIGHT, &
                                    NR_POINTS)

#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

      integer, intent(in) :: nr_dmat, nr_dmat_con
      type(matrix)        :: D(nr_dmat)
      type(matrix)        :: D_con(nr_dmat_con)
      dimension AO(*), &
     &          PX(*), &
     &          PY(*), &
     &          PZ(*), &
     &          buffer(*), &
     &          QUANTITY_INTEGRAL(0:3), &
     &          WEIGHT(*), &
     &          QUANTITY(0:3)


      DO IPOINT = 1,NR_POINTS

           call get_quantity_in_point(nr_dmat,     &
                                      D,           &
                                      nr_dmat_con, &
                                      D_con,       &
                                      px(ipoint),  &
                                      py(ipoint),  &
                                      pz(ipoint),  &
                                      ao,          &
                                      buffer,      &
                                      quantity)

        CALL DAXPY(4,WEIGHT(IPOINT),QUANTITY,1, &
     &                              QUANTITY_INTEGRAL,1)
      ENDDO


      RETURN
      END subroutine


      SUBROUTINE LEGENDRE_POLYNOMIAL(P,N,X)

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

      END subroutine


      SUBROUTINE LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X)

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

      SUBROUTINE LOBATTO_1D_GRID(N,X,W)

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(*), &
     &          W(*), &
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

      END subroutine



      SUBROUTINE LEGENDRE_1D_GRID(N,X,W)

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(*), &
     &          W(*), &
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

      END subroutine


   subroutine lobatto_2d_grid(nr_pieces, &
                              n,         &
                              p_origin,  &
                              p_right,   &
                              p_top,     &
                              u)

      implicit none

!     --------------------------------------------------------------------------
      integer, intent(in)  :: nr_pieces
      integer, intent(in)  :: n
      real(8), intent(in)  :: p_origin(3)
      real(8), intent(in)  :: p_right(3)
      real(8), intent(in)  :: p_top(3)
      integer, intent(in)  :: u
!     --------------------------------------------------------------------------
      integer              :: ipoint, ipiece, jpiece, i, j, m
      character(30)        :: file_name
      real(8), allocatable :: rx(:), ry(:), rz(:), rw(:), x(:), w(:)
      real(8)              :: vx(3), vy(3), p_off(3), p(3)
      real(8)              :: sx, sy
!     --------------------------------------------------------------------------

      write(file_name, '(a5,i2)') 'GRID.', u
      open(u,                   &
         file   = file_name,    &
         status = 'unknown',    &
         access = 'sequential', &
         form = 'unformatted')
      rewind u

      vx = p_right - p_origin
      vy = p_top   - p_origin
      sx = (vx(1)*vx(1) + vx(2)*vx(2) + vx(3)*vx(3))**0.5d0
      sy = (vy(1)*vy(1) + vy(2)*vy(2) + vy(3)*vy(3))**0.5d0

      m = nr_pieces*n

      call alloc(x,  m)
      call alloc(w,  m)
      call alloc(rx, m*m)
      call alloc(ry, m*m)
      call alloc(rz, m*m)
      call alloc(rw, m*m)

      call lobatto_1d_grid(n, x, w)

      ipoint = 0
      do ipiece = 0, nr_pieces-1
         do jpiece = 0, nr_pieces-1
            p_off = p_origin + ipiece*vx/nr_pieces + jpiece*vy/nr_pieces
            do i = 1, n
               do j = 1, n
                  ipoint = ipoint + 1

                  p = p_off                             &
                    + 0.5d0*(x(i) + 1.0d0)*vx/nr_pieces &
                    + 0.5d0*(x(j) + 1.0d0)*vy/nr_pieces

                  rx(ipoint) = p(1)
                  ry(ipoint) = p(2)
                  rz(ipoint) = p(3)
                  rw(ipoint) = w(i)*w(j)*sx*sy/((2.0d0*nr_pieces)**2.0d0)
               end do
            end do
         end do
      end do

      call wrtqua(rx, ry, rz, rw, u, ipoint)

      call dealloc(x)
      call dealloc(w)
      call dealloc(rx)
      call dealloc(ry)
      call dealloc(rz)
      call dealloc(rw)

   end subroutine

end module
