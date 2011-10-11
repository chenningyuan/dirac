      SUBROUTINE RECP_LNK_READCP(LUINFO,NUC1,Q)
!     Read ECP Input 
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
#include "recpval.h"
!     global variables
      INTEGER LUINFO,NUC1
      REAL*8  Q
!     local variables
      INTEGER I,J,K,IOS 
      CHARACTER*3 RECPKEYWORD, SOREP_CHR

!     Check RECP
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(BN,A3)',IOSTAT=IOS) RECPKEYWORD
      IF (IOS.NE.0) GOTO 8000

      IF (RECPKEYWORD.NE.'ECP') THEN
!        Discard the line we just read.
         NMLINE = NMLINE - 1
         GOTO 9000
      ENDIF

!     Allocate variables for RECP parameters
      CALL RECP_NTRA_READCP(MXCENT)
!
!     Explicit RECP parameter
!     =======================
!
      AREP_ANG(NUC1) = 0
      READ (MLINE(NMLINE),'(BN,4X,I3,1X,I3)',IOSTAT=IOS)
     &     RECP_CORE(NUC1), AREP_ANG(NUC1)
      IF (IOS.NE.0) GOTO 8000

      READ (MLINE(NMLINE),'(BN,12X,A3)',IOSTAT=IOS) SOREP_CHR
      IF (IOS.NE.0) GOTO 8000
!     Skip SOREP (blank) 
      IF (SOREP_CHR.EQ.'   ') THEN
         SOREP_ANG(NUC1) = 0
      ELSE
         READ (MLINE(NMLINE),'(BN,12X,I3)',IOSTAT=IOS) SOREP_ANG(NUC1)
         IF (IOS.NE.0) GOTO 8000
      ENDIF

!     Identify atom type from charge value
      RECP_CENT(NUC1) = NINT(Q) + RECP_CORE(NUC1)
!     Q = Q - RECP_CORE(NUC1)*1.0D0

!     ---------
!     Read AREP
!     ---------
      IF (AREP_ANG(NUC1).EQ.0) GOTO 9000
      DO I = 1, AREP_ANG(NUC1)+1
!        Read in each block         
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(BN,I3)',IOSTAT=IOS) AREP_BLK(NUC1,I)
         IF (IOS.NE.0) GOTO 8000
         DO J = 1, AREP_BLK(NUC1,I)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS)
     &           AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
            IF (IOS.NE.0) GOTO 8000
         ENDDO
      ENDDO

!     ---------------------
!     Read SOREP (if exist)
!     ---------------------
      IF (SOREP_ANG(NUC1).EQ.0) GOTO 9000 
      DO I = 1, SOREP_ANG(NUC1)
!        Read in each block         
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(BN,I3)',IOSTAT=IOS) SOREP_BLK(NUC1,I)
         IF (IOS.NE.0) GOTO 8000 
         DO J = 1, SOREP_BLK(NUC1,I)
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS)
     &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
            IF (IOS.NE.0) GOTO 8000 
         ENDDO
      ENDDO
      GOTO 9000

! Error message
 8000 CONTINUE
      WRITE(LUPRI,*) 'Error in reading line'
      WRITE(LUPRI,*) 'Line',NMLINE,': ',MLINE(NMLINE)
      CALL QUIT('Error in reading line')
      
 9000 CONTINUE

      END


      SUBROUTINE RECP_LNK_CTRAN(blk)
      USE RECP_NTR
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
      INTEGER soef, blk, I1,I2,I3,I4
!  REAL(8), ALLOCATABLE :: RECP_CTRAN(:, :, :, :)

!
!     Run over Nuclear Centres
!

!fixme: check whether user wants .ECP
!       if not, do not allocate and do the rest below

      ALLOCATE(ctranarg(NUCIND,mxqn,NUCIND*mxqn,NUCIND*mxqn))
      DO I1 = 1, NUCIND
         DO I2 = 1, mxqn
            DO I3 = 1, NUCIND*mxqn
               DO I4 = 1, NUCIND*mxqn
                  ctranarg(I1,I2,I3,I4) = 0.0d0
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 

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
               if ((ncenta.ne.NCENT(IA-1)) .or. 
     &             (NHKTA.ne.NHKT(IA-1))  ) then
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
                                 CTRANarg(icent,NHKTA,irow,icol) 
     &                            = PT(IBTAND(KB,IVARB))
                           soef = nint(CTRANarg(ICENT,NHKTA,irow,icol))
                                 ireparg(icent,NHKTA,irow) = LA + 1
!                                 write(lupri,'(a,i3,i3,i3,i3,a,i3)')
!    &                   'CTRANarg(',icent,NHKTA,irow,icol,') = ',soef
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
      END

