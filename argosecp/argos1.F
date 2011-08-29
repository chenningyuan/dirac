!*******************************************************************
!
!     ECPINT : RECP integral routine
!
!     This program is based on argos integral routine
!     with the permission of R. M. Pitzer
!     and maintained by Young Choon Park (KAIST)
!
!*******************************************************************
      SUBROUTINE ECPINT(WORK,LWORK)

#include "implicit.h"
#include "argos_print.h"

!     Memory allocation
!                 |            |           |     |
!     CALL MEMGET('REAL',KBASE,LENGTH,WORK,KFREE,LFREE)
!                        |                 |     |
!
!     Memory deallocation
!     CALL MEMREL('ARGOS',WORK,KWORK,KWORK,KFREE,LFREE)

      DIMENSION WORK(*)
      INTEGER   LWORK

!     Enter ECPINT
      CALL QENTER('ECPINT')
      CALL TITLER('RECP integral start','*',126)
      CALL ARGOS_TITLE   !     Print title

!     Basis setup
      PRINT_LEVEL = 5 !10 
      ARGOS_DEBUG = 0 !1

!     RECP_INPUT generation
      CALL ARGOSIN

!     Main part
      CALL ARGOS_DRIVER(WORK,LWORK) 
 
!     Exit ECPINT
      CALL TITLER('RECP integral end','*',126)
      CALL QEXIT('ECPINT')

      RETURN      
      END

! ------------------------------------------------

      SUBROUTINE ARGOS_TITLE
!
!     Write title and memory info
!
      IMPLICIT NONE
#include "argos_print.h"

      INTEGER NLIST

      NLIST = 6 ! output file

      WRITE(NLIST,100)
  100 FORMAT(/,
     & 5x,'Relativistic Effective Core Potential Integral routine'/
     & 5x,'Version 2011.01.18'/
     & 5x,'This routine is based on ARGOS integral'/
     & 5x,'with the permission of R. M. Pitzer (Ohio State University)'/
     & 5x,'and maintained by Y. C. Park (KAIST, youngc@kaist.ac.kr)'/) 

!     IF (PRINT_LEVEL .GE. 1) THEN
!        WRITE(NLIST,200) LCORE 
! 200    FORMAT(/' Total memory allocated for RECP integral : '
!    &          ,i10)
!     ENDIF

      RETURN
      END

! ------------------------------------------------

      SUBROUTINE ARGOS_DRIVER(WORK,LWORK)

      implicit logical(a-z)

#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"
#include "argos_print.h"

      DIMENSION WORK(*)
      INTEGER   LWORK

      integer
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst
      common /parmi/
     & mblu,  mccu, mconu, mcu,  mcxu, mpru, nbft, nnbft,
     & mrcru, mstu, msu,   nblu, ng,   ns,   nst

      INTEGER     LCORED
!     PARAMETER ( LCORED = 4 500 000 )
      PARAMETER ( LCORED = 14 500 000 )
      INTEGER LCORE
      REAL*8  A
      COMMON // A(LCORED)

      integer
     & nf(msup),      nc(msup),      ncon(mconsp),    nrcr(mconsp),
     & lmnp1(mconsp), ntl(msfup),    ntu(msfup),      nt(msfup),
     & mcons(msfup),  iprst(mnsfup), npair(2,mnsfup), icxst(2,mnsfup),
     & nfct(mccup),   ngw(mccup),    icb(4,24,mccup),
     & ica(mcup,msup,mgup),          nprir(2,mstup,mnsfup)

!     # ipt(*) is used for workspace stack pointers.
      integer ipt(55)

      real*8
     & zet(mconup,mconsp), eta(mrcrup,mconup,mconsp),
     & COOR_X(mcup,msup), COOR_Y(mcup,msup), COOR_Z(mcup,msup)

      LCORE  = LCORED

      msu   = msup
      mstu  = mstup
      mrcru = mrcrup
      mconu = mconup
      mcu   = mcup
      mccu  = mccup
      mblu  = lcore

      ipt(1) = 1    ! allocate space for lmnv

      call argos_1e( WORK, LWORK, A, eta, ica,    icb,    icxst,
     &   iprst,  ipt,    lmnp1, A(ipt(1)), mcons,  nc, ncon,
     &   nf,     nfct,   ngw,    npair,     nprir,  nrcr,   nt,
     &   ntl,    ntu,    COOR_X, COOR_Y, COOR_Z,     zet,    lcore )

      return
      end


