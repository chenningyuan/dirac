      SUBROUTINE RASSG3_REL(CB,SB,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                      LLUC,LLUHC)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to SBLOCK
*
*               March 1998: Relativistic version
*                           real/imaginary loops,
*                           MK2 coupling
*               (Jeppe Olsen, with a little help from my friend)
*
*
*=========================================================
* SUMMARY OF THE VECTOR HANDLING (I/O) IN SIGMA ROUTINES:
*=========================================================
*
* rassg3:	Set current sigma batch to zero (no disc operation)
*                 CALL SETVEC(SB,ZERO,LS)
* rassg3:	Copy blocks with given MK2_C to LUSC37
*                 CALL COPNBLKD(LLUC,LUSC37,CB,NBL_C,0,LBLK)
* rassg3:       Call to sblock with given IRIS,MK2_S,IRIC,MK2_C
* rassg3:	Transfer S block to permanent storage
*                 CALL TODSC(SB(IOFF),LEN_T,-1,LLUHC)
*=========================================================
*
*=========================================================
* OVERVIEW OVER TRANSPOSITIONS OF S AND C BLOCKS
*=========================================================
*
* The blocks are assumed transposed at the beginning, meaning
* a global transposition is needed in sblocks/dblocks.
* General rule: For processing a beta block, a transposition is needed.
*
* sblocks:	Offdiagonal blocks, explicit transposing
*               (Always done with PSSIGN and PLSIGN = 0 => IDC = 1
*               CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
*               CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
*                - - - - - - - - - - - - - - - - - - - - - - -
* rssbcbrel:	  <beta || beta>      TRP(S)  TRP(C)
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
* sblocks:	Transpose final sigma block back.
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
#include "cstate.inc"
#include "clunit.inc"
#include "gasstr.inc"
#include "cprnt.inc"
*
*.Definition of c and sigma spaces
#include "cands.inc"
*
* IA/BS/CTP : Type (number of elecs in A/B,S/C)
*. CSTATE
C     COMMON/CSTATE/IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,
C    &              NROOT,IROOT(100),IDC,PSSIGN,PLSIGN,INTSEL,IALIGN,
C    &              ETHRES,CTHRES,EROOT(100),NMS2VAL,MS2VAL(100)

*
      CALL QENTER('RASSG3')
*
      NTESTL = 0005
      NTEST = MAX(NTESTL,IPRSIG)
*
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' RASSG3 speaking :'
        WRITE(6,*) ' ================='
        write(6,*)
        WRITE(6,*) ' RASSG3 : NBATS = ',NBATS
      END IF
*
      LBLK = -1
      CALL REWINO(LLUC)
      CALL REWINO(LLUHC)
*
*. Loop over vectors of sigma and C : MK2 space + Real/Imag
      DO IRIS = 1, 2
        JBATABS = 0
        DO IMK2_S = 1, NMS2VAL
*
          MK2_S = MS2VAL(IMK2_S)
          IASTP = IST_FOR_DT(1,IMK2_S)
          IBSTP = IST_FOR_DT(2,IMK2_S)
*. Loop over batches of sigma for given IRIS and IMK2_S
          DO JBATS = 1, NBAT_MS2(IMK2_S)
            JBATABS = JBATABS + 1
            if (NTEST.ge.500) then
               WRITE(6,*)
               WRITE(6,*)
     &         ' =============================================='
               WRITE(6,*)
     &         ' Start of construction of Sigma batch ', JBATABS
               WRITE(6,*)
     &         ' =============================================='
               WRITE(6,*)
            end if
*. Zero sigma batch
            LS = LEBATS(JBATABS)
            ZERO = 0.0D0
            CALL SETVEC(SB,ZERO,LS)
*
            CALL REWINO(LLUC)
            DO IRIC = 1, 2
              IF(IRIC.EQ.2) THEN
*. Skiup EOV mark between real and imaginary part
                 CALL IFRMDS(IONEM,1,-1,LLUC)
              END IF
              DO IMK2_C = 1, NMS2VAL
*
                IF(IRIC.EQ.2.AND.IMK2_C.EQ.NMS2VAL) THEN
                  I_DO_TRNSP = 1
                ELSE
                  I_DO_TRNSP = 0
                END IF
*
                MK2_C = MS2VAL(IMK2_C)
                IACTP = IST_FOR_DT(1,IMK2_C)
                IBCTP = IST_FOR_DT(2,IMK2_C)
*. Copy blocks with given MK2_C to LUSC37
                CALL REWINO(LUSC37)
                NBL_C = NBLK_MS2(IMK2_C)
                CALL COPNBLKD(LLUC,LUSC37,CB,NBL_C,0,LBLK)
C                    COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*. End of vector mark
                CALL ITODS(-1,1,-1,LUSC37)
                if (NTEST.ge.50) then
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
     &            ' MK2/IRI vector of S :  MK2_S, IRIS =', MK2_S,IRIS
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of C :  MK2_C, IRIC =', MK2_C,IRIC
                  WRITE(6,*)
                end if
*
C               rewind LUSC37
C               CALL WRTVCD(CB,LUSC37,1,-1)
                IF(ABS(MK2_S-MK2_C).LE.4) THEN
*
*. Connection. The Hamiltonian is at most a two-particle operator!
*. Obtain contribution to batch sigma from all blocks with MK2_C
CM                WRITE(6,*) ' Memchk before call to SBLOCK '
CM                CALL LMEMCHK
                  IF(NTEST.GE.100) THEN
                   WRITE(6,*) ' Block info for S batch, JBATS:',JBATABS
                   WRITE(6,*) ' Number of blocks : ', LBATS(JBATABS)
                   WRITE(6,*) ' Offset of batch ', I1BATS(JBATABS)
                   WRITE(6,*) ' TTSSOOLL array for blocks '
                   LLX = LBATS(JBATABS)
                   CALL IWRTMA(IBATS(1,I1BATS(JBATABS)),8,LLX,8,LLX)
                  END IF
*
                  CALL SBLOCK_REL(LBATS(JBATABS),
     &                            IBATS(1,I1BATS(JBATABS)),
     &                            1,
     &                            CB,SB,LUSC37,0,0,0,0,
     &                            0,IRIS,IRIC,I_DO_TRNSP)
*
                  if (NTEST.ge.100) then
                    WRITE(6,*) ' Updated array containing sblocks '
                    CALL WRTMAT(SB,1,LS,1,LS)
                  end if
                END IF
*               ^End if interacting
              END DO
*             ^End of loop over IMK2_C
            END DO
*           ^End of loop over real/imaginary part of C
*. A batch of sigma has been constructed so save
*. Transfer S block to permanent storage
            NTEST = 10
            if (NTEST.ge.10) then
               WRITE(6,*)
               WRITE(6,*) ' **********************************'
               WRITE(6,*) ' *                                *'
               WRITE(6,*) ' * Array containing final sbatch  *'
               WRITE(6,*) ' *                                *'
               WRITE(6,*) ' **********************************'
               WRITE(6,*)
               CALL WRTMAT(SB,1,LS,1,LS)
                  stop
            end if
            KBATSEND = I1BATS(JBATABS) + LBATS(JBATABS)
            DO ISBLK = I1BATS(JBATABS),KBATSEND-1
              IOFF = IBATS(6,ISBLK)
              LEN_T  = IBATS(8,ISBLK)
C?            write(6,*) 'RASSG3 : IOFF, SB(IOFF)',IOFF,SB(IOFF)
              CALL ITODS(LEN_T,1,-1,LLUHC)
              CALL TODSC(SB(IOFF),LEN_T,-1,LLUHC)
            END DO
          END DO
*         ^End of batches for given MK2_S
        END DO
*       ^End of loop over IMK2_S
*. End of real/imag part of sigma => end of vector
        CALL ITODS(-1,1,-1,LLUHC)
      END DO
*     ^End of loop over real/imag part of sigma
*
      IF(NTEST.GE.01) THEN
        WRITE(6,*) ' Final S-vector on disc'
        CALL WRTVCDC(SB,LLUHC,1,-1)
      END IF
*
      CALL QEXIT('RASSG3')
      RETURN
      END
