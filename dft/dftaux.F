      subroutine fortwrt(str, slength)

      implicit none

#include "priunit.h"

#ifdef INT_STAR8
      integer(8) :: slength
#else
      integer    :: slength
#endif /* ifdef INT_STAR8 */
      character  :: str(slength)

      write(lupri, *) str

      end subroutine

      SUBROUTINE DFTSETCAM(W,BETA)
#include "implicit.h"
#include "dcbham.h"
      HFXMU = W
      HFXATT = BETA
      END

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
