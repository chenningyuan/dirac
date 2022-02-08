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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GTDOAV(mat_dim,
     &                  nz,
     &                  nr_dmat,
     &                  dmat,
     &                  ifac)
!***********************************************************************
!
!     Get average density from open shells; only real part !
!     Written by T.Saue March 2004
!
!***********************************************************************
      implicit none
#include"dcbdhf.h"
      integer :: mat_dim, nz, nr_dmat
      integer :: ifac, ishell
      real(8) :: DMAT(mat_dim, mat_dim, nz, nr_dmat)
      real(8) :: fac
      IF(AOC) THEN
        DO ISHELL = 1, NOPEN
          FAC = IFAC*DF(ISHELL)
          CALL DAXPY(mat_dim*mat_dim,FAC,
     &        DMAT(1,1,1,ISHELL+1),1,DMAT,1)
        ENDDO
      ENDIF
!
      end subroutine
