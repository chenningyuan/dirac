      SUBROUTINE SIGDEN_REL2(CB,SB,JBATABS,NBATS,LBATS,LEBATS,
     &                       I1BATS,IBATS,
     &                       LUC,LUHC,T,ISIGDEN,C2,IACTP,IBCTP,
     &                       IRIS,IRIC,IMK2_C,LUSCR_L,IPRNT)
*
* First inner routine for general CC sigma/densi
*
* Complete S and C vectors
*
* Jeppe Olsen   Relativistic version, Winter of 2001
*
* =====
* Input
* =====
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "gasstr.inc"
#include "cstate.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "csm.inc"
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*.Scratch
      DIMENSION SB(*),CB(*)
      dimension T(*)
*
      CALL QENTER('SIGD2')
*
      NTESTL = 000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SIGDEN2 speaking :'
        WRITE(6,*) ' ================='
      END IF
*
*
* Initialize sigma
*
* ISIGDEN = 2 : Read in batch of blocks of lhs vector
*. Unique blocks are read in, 
*  expanded and normalized to determinant form
*
      if (ISIGDEN.eq.2) then
        if (IRIC.eq.1.and.IMK2_C.eq.1) then
          do ISBLK = 1,LBATS(JBATABS),1
            ISATP = IBATS(1,ISBLK)
            ISBTP = IBATS(2,ISBLK)
            ISASM = IBATS(3,ISBLK)
            ISBSM = IBATS(4,ISBLK)
            ISOFF = IBATS(5,ISBLK)
            LEN2  = IBATS(7,ISBLK)
*
            XDUM = 0.0D0
            ISCALE = 1
*. Signs taken from CSTATE
            CALL GSTTBLD(SB(ISOFF),ISATP,ISASM,ISBTP,ISBSM,
     &            WORK(KNSTSO(IACTP)),WORK(KNSTSO2(IBCTP)), 
     &            PSSIGN,IDC,PLSIGN,LUHC,C2,
     &            NSMST,ISCALE,XSCALE)
*
            if (LEN2.gt.0.and.NTEST.ge.100) then
              write(6,*) 'WRITING L VECTOR FRACTION:'
              call wrtmat(SB(ISOFF),1,LEN2,1,LEN2)
            end if
*
          end do
*         ^ End of loop over S-blocks in batch
        end if
*       ^ Only perform read if new S batch requested
*
      end if
*     ^ Only for density evaluation
*
*. Obtain sigma/density for batch of blocks
      CALL SIGDEN_REL3(LBATS(JBATABS),IBATS(1,I1BATS(JBATABS)),1,
     &                 CB,SB,LUC,T,ISIGDEN,IRIS,IRIC)
*
      CALL QEXIT('SIGD2')
      RETURN
      END
