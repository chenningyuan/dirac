      SUBROUTINE DASSG3(RHO1UUR,RHO1UUI,
     &                  RHO1UBR,RHO1UBI,
     &                  RHO1BBR,RHO1BBI,
     &                  RHO1BUR,RHO1BUI,
     &                  CB,SB,SCR,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                  INONREL,IDEN12,IFINMO,IPRNT)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to DBLOCK
*
*               March 1998: Relativistic version
*                           real/imaginary loops,
*                           MK2 coupling
*
*               Density version in May 1998, TF and JO
*=========================================================
* OVERVIEW OVER TRANSPOSITIONS OF S AND C BLOCKS
*=========================================================
*
* The blocks are assumed transposed at the beginning, meaning
* a global transposition is needed in sblocks/dblocks.
* General rule: For processing a beta block, a transposition is needed.
*
* dblocks:      Offdiagonal blocks, explicit transposing
*               (Always done with PSSIGN and PLSIGN = 0 => IDC = 1
*               CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
*               CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
*                - - - - - - - - - - - - - - - - - - - - - - -
* rssbcbrel:      <beta || beta>      TRP(S)  TRP(C)
*                                     do matrix operations
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <alpha || alpha>    do matrix operations
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <alpha || beta>     TRP(S)  TRP(C)
*                                     TRP(S)
*                                     do matrix operations
*                                     TRP(S)
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <beta || alpha>             TRP(C)
*                                     do matrix operations
*                                             TRP(C)
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
* sblocks:      Transpose final sigma block back.
*
*
* =====
* Input
* =====
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*.Scratch
      DIMENSION SB(*),CB(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cstate.inc"
#include "clunit.inc"
#include "cands.inc"
#include "gasstr.inc"
#include "strbas.inc"
*
      CALL QENTER('DASSG')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' DASSG3 speaking :'
        WRITE(6,*) ' ================='
        write(6,*)
        WRITE(6,*) ' DASSG3 : NBATS = ',NBATS
      END IF
*
*  LUSC62: logical unit right (C vector)
*  LUSC61: logical unit left  (S vector, identical)
*
      LBLK = -1
      CALL REWINO(LUSC61)
*. Loop over vectors of sigma and C : MK2 space + Real/Imag
      DO IRIL = 1,2,1
        if (NTEST.ge.50) write(6,*) 'left r/i is  ',IRIL
        if (IRIL.eq.2) then
*. SKIP EOV mark between real and imaginary part of S
          CALL IFRMDS(IONEM,1,-1,LUSC61)
        end if
*
        JBATABS = 0
        do IMK2_L = 1,NMS2VAL,1
          MK2_L = MS2VAL(IMK2_L)
          IASTP = IST_FOR_DT(1,IMK2_L)
          IBSTP = IST_FOR_DT(2,IMK2_L)
*. Loop over batches of sigma for given IRIL and IMK2_L
          DO JBATS = 1,NBAT_MS2(IMK2_L),1
            JBATABS = JBATABS + 1
            if (NTEST.ge.20) then
              WRITE(6,*)
              WRITE(6,*)
     &        ' =============================================='
              WRITE(6,*)
     &        ' Start of construction from  Left batch ', JBATABS
              WRITE(6,*)
     &        ' =============================================='
              WRITE(6,*)
            end if
*. Read left batch in
            do ISBLK = 1,LBATS(JBATABS),1
              ISATP = IBATS(1,ISBLK)
              ISBTP = IBATS(2,ISBLK)
              ISASM = IBATS(3,ISBLK)
              ISBSM = IBATS(4,ISBLK)
              ISOFF = IBATS(5,ISBLK)
              LEN2  = IBATS(7,ISBLK)
*
              print*,'ISBLK ',ISBLK
              print*,'ISATP,ISBTP,ISASM,ISBSM,ISOFF',
     &                ISATP,ISBTP,ISASM,ISBSM,ISOFF
              XDUM = 0.0D0
              ISCALE = 1
*. Signs taken from CSTATE
              CALL GSTTBLD(SB(ISOFF),ISATP,ISASM,ISBTP,ISBSM,
     &              WORK(KNSTSO(IASTP)),WORK(KNSTSO2(IBSTP)),
     &              PSSIGN,IDC,PLSIGN,LUSC61,SCR,
     &              NSMST,ISCALE,XSCALE)
              print*,'TESTING L VECTOR FRACTION:'
              call wrtmat(SB(ISOFF),1,LEN2,1,LEN2)
            end do
*           ^ End of loop over S-blocks in batch
*
            call REWINO(LUSC62)
            do IRIR = 1,2,1
              if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
              if (IRIR.EQ.2) then
*. Skip EOV mark between real and imaginary part
                call IFRMDS(IONEM,1,-1,LUSC62)
              end if
              do IMK2_R = 1,NMS2VAL,1
                MK2_R = MS2VAL(IMK2_R)
                IACTP = IST_FOR_DT(1,IMK2_R)
                IBCTP = IST_FOR_DT(2,IMK2_R)
*. Copy blocks with given MK2 to LUSC37
*  Explicit selection and access done in dblocks (analogy to sblocks)
                CALL REWINO(LUSC37)
                NBL_C = NBLK_MS2(IMK2_R)
                CALL COPNBLKD(LUSC62,LUSC37,CB,NBL_C,0,LBLK)
C                    COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*. End of vector mark
                CALL ITODS(-1,1,-1,LUSC37)
                if (NTEST.ge.20) then
                  WRITE(6,*)
     &            ' ********************************************* '
                  WRITE(6,*)
     &            ' *                                           * '
                  WRITE(6,*)
     &            ' * Start of contribution from given CISPACE  * '
                  WRITE(6,*)
     &            ' *                                           * '
                  WRITE(6,*)
     &            ' ********************************************* '
                  WRITE(6,*)
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of S :  MK2_L, IRIL =', MK2_L,IRIL
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of C :  MK2_R, IRIR =', MK2_R,IRIR
                  WRITE(6,*)
                end if
*
C               CALL WRTVCD(CB,LUSC37,1,-1)
                IF(ABS(MK2_L-MK2_R).LE.2) THEN
*.    Connection
CTF  What about 2-particle densities?
CTF  Shouldn't ABS be .le. 4 instead?
*
*. Obtain contribution to batch sigma from all blocks with MK2_R
                  CALL DBLOCK(INONREL,IDEN12,IRIL,IRIR,
     &                        RHO1UUR,RHO1UUI,
     &                        RHO1UBR,RHO1UBI,
     &                        RHO1BBR,RHO1BBI,
     &                        RHO1BUR,RHO1BUI,
     &                        LBATS(JBATABS),IBATS(1,I1BATS(JBATABS)),
     &                        1,CB,SB,LUSC37,0,0,0,0)
*
                end if
*               ^End if interacting
              end do
*             ^End of loop over IMK2_R
            end do
*           ^End of loop over real/imaginary part of C
          end do
*         ^End of batches for given MK2_L
        end do
*       ^End of loop over IMK2_L
*. End of real/imag part of sigma => end of vector
      end do
*     ^End of loop over real/imag part of sigma
*
      call QEXIT('DASSG')
      return
      end
