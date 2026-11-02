      SUBROUTINE get_vemb_mat(vemb_mat)
!
!     Get embedding potential matrix from file CHECKPOINT
!     Packed on fermion irreps
!
!     andre gomes, october 2012
!
      use mxcent_mod, only: MXCENT
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use dcbgen_mod, only: LU1INT
!      use cbihr1_mod, only: IPRONE
!      use dcbbas_mod, only: IBAS, IBBAS, N2BAS, N2BAST, N2BBASX, NBBAS, NFBAS, NNBBASX, NTBAS
!      use dgroup_mod, only: JFSYM, NFSYM, NZ
!      use symmet_mod, only: ICOS
!      use dcbham_mod, only: NOSMLV
      use memory_allocator

      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, J
      REAL(dp), PARAMETER :: D1 = 1.0D0
C
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
#include "dcbham.h"
      DIMENSION vemb_mat(*)
      CHARACTER*8 RTNLBL(2)

      real(8), allocatable :: stri(:)
      real(8), allocatable :: sbuf(:)
C      
C     Allocate integral buffer
C     ========================
C
      call alloc(stri, nnbbasx)
      call alloc(sbuf, n2bbasx)
C
C     Read integral from file
C     =======================
C
      if (NOSMLV) then
      CALL PRPREA(LU1INT,'FDEVEMB ',RTNLBL,'TFFF',STRI,NNBBASX,
     &            IPRONE)
      else
      CALL PRPREA(LU1INT,'FDEVEMB ',RTNLBL,'TFFT',STRI,NNBBASX,
     &            IPRONE)
      endif
C
C     Make full matrix
C     ================
C
      CALL DSPTSI(NTBAS(0),STRI,SBUF)
C
      CALL DZERO(vemb_mat,N2BAST)
      NBRP     = 4/NZ
      IOFF = 0
      DO I = 1,NFSYM
        DO IC = 1,2
          IP    = MOD(I+IC,NFSYM) + 1
          DO JSYM = 1,NBRP
            ISYM           = JFSYM(JSYM,IP)
            IREP           = ISYM - 1
!fixme irep can be < 0 in NONREL runs
            IF(NBBAS(IREP,IC).GT.0) THEN
              ISOFF = IOFF + (IBBAS(IREP,IC)-IBAS(I))*(NFBAS(I,0)+1)
              IUOFF = ICOS(ISYM,IC)*(NTBAS(0)+1)
              DO J = 1,NBBAS(IREP,IC)
                CALL DCOPY(NBBAS(IREP,IC),SBUF(1+IUOFF),1,
     &                                    vemb_mat(1+ISOFF),1)
                IUOFF = IUOFF + NTBAS(0)
                ISOFF = ISOFF + NFBAS(I,0)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IOFF = IOFF + N2BAS(I)
      ENDDO
C
C     Memory deallocation
C     ===================
C
      call dealloc(stri)
      call dealloc(sbuf)
      END
