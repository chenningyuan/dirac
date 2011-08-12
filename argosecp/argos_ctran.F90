MODULE RECP_INTERFACE
!  real(8), allocatable :: ctranarg(:, :, :, :)
   REAL(8), ALLOCATABLE :: ctranarg(:, :, :, :)

CONTAINS

SUBROUTINE RECP_INTERFACE_CTRAN(blk)

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pincom.h"
#include "aosotr.h"
#include "sphtrm.h"
#include "huckel.h"
#include "cbirea.h"
 
#include "ibtfun.h"
#include "argoscom.h"
#include "argoscomch.h"
   INTEGER soef, blk
!  REAL(8), ALLOCATABLE :: RECP_CTRAN(:, :, :, :)

!
!     Run over Nuclear Centres
!

!fixme: check whether user wants .ECP
!       if not, do not allocate and do the rest below

      ALLOCATE(ctranarg(NUCIND,mxqn,NUCIND*mxqn,NUCIND*mxqn))

      blk = 0
      do 401 ICENT = 1, NUCIND   ! run over nuclear center
         nucdega = nucdeg(ICENT)
         DO 411 IA = 1, KMAX     ! run over shells
            ila = 0
            icol = 1
            MULA   = ISTBAO(IA)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            NUMCFA = NUMCF(IA)
            IC     = LCLASS(IA)
            ncenta = NCENT(IA)
            irowmax(icent,NHKTA) = nucdega * KHKTA
!YCP        IF (IA-1.GE.1) THEN ! arrays boundaries fix
            if (ncenta .eq. ICENT) then
               if ((ncenta .ne. NCENT(IA-1)) .or. (NHKTA .ne. NHKT(IA-1))) then
                  blk = blk +1
                  IF (NUMCFA .EQ. 1) IBLOCK = IBLOCK + 1
                  IOFF = NAORBD
                  DO 421 NA = 1, KHKTA      ! run over components
                     ila = 0
                     DO 431 LA = 0, MAXREP  ! run over symmetry operator
                        IVARB = IBTXOR(LA,ISYMAO(NHKTA,NA))
!
!              If orbital contributes to this representation:
!
!                (IVARB is a bitstring of basic operations with
!                1 in the positions where the basic operations
!                has a different parity for the irrep and the
!                AO(when centered in origo). These basic
!                operations with different parities must then
!                not be part of the stabilizer of the center
!                for the AO).
!
                        IF (IBTAND(MULA,IVARB) .EQ. 0) THEN
                           ila = ila + 1
                           JKB = 0
!
!                 loop over symmetry dependent centers
!
                           DO 441 KB = 0, MAXREP
                              IF (IBTAND(KB,MULA) .EQ. 0) THEN
                                 JKB              = JKB + 1
                                 irow = NA + KHKTA*(ila-1)
                                 icol = NA + KHKTA*(JKB-1)
                                 CTRANarg(icent,NHKTA,irow,icol) = PT(IBTAND(KB,IVARB))
                           soef = nint(CTRANarg(ICENT,NHKTA,irow,icol))
                                 ireparg(icent,NHKTA,irow) = LA + 1
!                                 write(lupri,'(a,i3,i3,i3,i3,a,i3)')
!     &                   'CTRANarg(',icent,NHKTA,irow,icol,') = ',soef
                              maxang(icent) = NHKTA
                              END IF
  441                      CONTINUE
                        ELSE
                        END IF
  431                CONTINUE
  421             continue
               endif
            endif
!YCP        ENDIF
  411    CONTINUE
  401 CONTINUE

END SUBROUTINE RECP_INTERFACE_CTRAN

END MODULE RECP_INTERFACE
