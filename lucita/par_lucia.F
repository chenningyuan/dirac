C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
#if defined (VAR_MPI)
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS_FG(LUIN,LBLK,BLK_A,SEGMNT,ISTBLK)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*  slave version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*      INTEGER ISTBLK
*
      CALL REWINE(LUIN,LBLK)
*
      IBLK = 0
      IBLK = ISTBLK
      NBLK_A = 0
      NBLK = 0
*      write(6,*)'ISTBLK',ISTBLK, myproc
*. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
*        write(6,*)'IBLK =',IBLK
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
      NBLK =  (IBLK - ( 1 + ISTBLK))
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,'(A,I8,I8)')
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active blocks (slave)',myproc
        CALL WRTMAT(BLK_A,1,NBLK+ISTBLK,1,NBLK+ISTBLK)
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_BLOCK_DISTR(JCTCMP,MNSBLOCK,JISBLOCK,
     &           JICOCOC,JICSMOS,JICBLTP,MNSSOA,MNSSOB,
     &           MNOCTPA,MNOCTPB,
     &           MNSMST,JLC,JIDC,JIDOH2,JISTRFL,
     &           PPS,LLCBLOCK,LLECBLOCK,JI1CBLOCK,JICBLOCK,
     &           JIRESTRICT,JICONSPA,JICONSPB,TSCLFAC,
     &           JICBAT_RES,LCBAT_INI,LCBAT_END,JCBLINF,
     &           MNJJJBLK,MTEMPP,JMPINFA,LLUC,CCB,CC2,CC,NYTOTB,
     &           ICSBTCHB)
*
*
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*. Specific input
      INTEGER JISBLOCK(8,*)
*.General input
      INTEGER JICOCOC(MNOCTPA,MNOCTPB)
      INTEGER JICSMOS(MNSMST)
      INTEGER JICBLTP(*)
      INTEGER MNSSOA(MNSMST ,*), MNSSOB(MNSMST ,*)
      INTEGER JICONSPA(MNOCTPA,MNOCTPA), JICONSPB(MNOCTPB,MNOCTPB)
*.Scratch
      INTEGER   LLCBLOCK(*),JI1CBLOCK(*),JICBLOCK(8,*),LLECBLOCK(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION TSCLFAC(*)
* new scratch
      INTEGER MNJJJBLK,ILOOPSB,MTEMPP,ISCRLNGTH
      DIMENSION JCBLINF(MNJJJBLK)
      DIMENSION CCB(*),CC2(*)
* switch between 'send' and 'broadcast' of the CB-array
*      INTEGER I_USE_BCAST
*  ==========
*.   output 
*  ==========
*
      DIMENSION JCTCMP(MNJJJBLK,4), JMPINFA(MTEMPP,2), NYTOTB(MNJJJBLK)
      DIMENSION ICSBTCHB(MNJJJBLK)
*   
* JCTCMP: array containing information about a c-block, whether it is
* contributing, together with its length; 
* third number specifies how often it interacts with a given s-block    
*
*.
      ISCRLNGTH = 0
      ILOOPSB = 0
*      I_USE_BCAST = 1
CCC!      IF(MYPROC.EQ.MASTER)THEN
CCC!        REWIND LLUC
CCC!      END IF
*      WRITE(6,*)'JIDC', JIDC
* ===========================
* 1 : Arrays for accessing C
* ============================
*. Find batches of C - strings
*. Should be changed to REPART or ??
*. Not neccesary, all zero blocks detected otherwise
      ITTSS_ORD = 2
      CALL PART_CIV2(JIDC,JICBLTP,MNSSOA,MNSSOB,MNOCTPA,MNOCTPB,MNSMST,
     &              JLC,JICOCOC,JICSMOS,MCBATCH,LLCBLOCK,LLECBLOCK,
     &              JI1CBLOCK,JICBLOCK,0,ITTSS_ORD)
*
*
* Loop over batches over C blocks
      IF(JIDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
      IF(JICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' LCBAT_INI LCBAT_END', LCBAT_INI,LCBAT_END
        JCBAT_INI = LCBAT_INI
        JCBAT_END = LCBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = MCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END

*. Read C blocks into core
*
        ICOFF = 1
        NJBLOCK = LLCBLOCK(JCBATCH)
*        WRITE(6,*)'MCBATCH = ',MCBATCH
*        WRITE(6,*)'NJBLOCK =',NJBLOCK
        DO JJCBLOCK = 1, NJBLOCK
          JBLOCK = JI1CBLOCK(JCBATCH)-1+JJCBLOCK
*. Will this block be needed ??
          INTERACT = 0
*          WRITE(6,*)'TSCLFAC(JBLOCK) =',TSCLFAC(JBLOCK)
          IF(TSCLFAC(JBLOCK).EQ. 1.0D0) THEN
C?          WRITE(6,*) ' Nonvanishing block ', JBLOCK, myproc
          JATP = JICBLOCK(1,JBLOCK)
          JBTP = JICBLOCK(2,JBLOCK)
          JASM = JICBLOCK(3,JBLOCK)
          JBSM = JICBLOCK(4,JBLOCK)
          JOFF = JICBLOCK(5,JBLOCK)
          CALL PRMBLK(JIDC,JISTRFL,JASM,JBSM,JATP,JBTP,PPS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
            DO JSBLOCK = 1, MNSBLOCK
            IDENT = 0
            IF(JISBLOCK(1,JSBLOCK).GT.0) THEN
              IATP = JISBLOCK(1,JSBLOCK)
              IBTP = JISBLOCK(2,JSBLOCK)
              IASM = JISBLOCK(3,JSBLOCK)
              IBSM = JISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
              IF(MXEXC.EQ.1) THEN
                 IF((JICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (JICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((JICONSPA(IATP,LLATP).LE.1.AND.
     &               JICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (JICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (JICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
*. Are they identical ?
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
          IF(INTERACT.EQ.1) THEN
*
            JCTCMP(JBLOCK,1) = 1
*          write(6,*)'active: JCTCMP(JBLOCK,1) =',JCTCMP(JBLOCK,1),myproc
            JCTCMP(JBLOCK,2) = JCBLINF(JBLOCK)
*          write(6,*)'active length: JCTCMP(JBLOCK,2) =',JCTCMP(JBLOCK,2)
          ELSE
*
            JCTCMP(JBLOCK,1) = 0
*            write(6,*)'inact: JCTCMP(JBLOCK,1) =',JCTCMP(JBLOCK,1)
            JCTCMP(JBLOCK,2) = JCBLINF(JBLOCK)
            TSCLFAC(JBLOCK) = 0.0D0
*
          END IF
*
        END DO
*
*
       DO 9000 ICBLK = JI1CBLOCK(JCBATCH), JI1CBLOCK(JCBATCH)-1+NJBLOCK
*
*            WRITE(6,*)'ICBLK 9000loop =',ICBLK
            JATP = JICBLOCK(1,ICBLK)
            JBTP = JICBLOCK(2,ICBLK)
            JASM = JICBLOCK(3,ICBLK)
            JBSM = JICBLOCK(4,ICBLK)
            ICOFF = JICBLOCK(5,ICBLK)
            NJA = MNSSOA(JASM,JATP)
*            WRITE(6,*)'NJA =',NJA
            NJB = MNSSOB(JBSM,JBTP)
*
            IF(TSCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
              CALL PRMBLK(JIDC,JISTRFL,JASM,JBSM,JATP,JBTP,PPS,PL,
     &                    LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
                ILOOPSB = 0
*                WRITE(6,*)'NPERM 9000loop=',NPERM
*
              DO 8765 IPERM = NPERM,1, -1
                LLASM = LASM(IPERM)
                LLBSM = LBSM(IPERM)
                LLATP = LATP(IPERM)
                LLBTP = LBTP(IPERM)
                NLLA = MNSSOA(LLASM,LLATP)
                NLLB = MNSSOB(LLBSM,LLBTP)
*
                DO 10000 ISBLK = 1, MNSBLOCK
                IF(JISBLOCK(1,ISBLK) .GT. 0 ) THEN
                  IATP = JISBLOCK(1,ISBLK)
                  IBTP = JISBLOCK(2,ISBLK)
                  IASM = JISBLOCK(3,ISBLK)
                  IBSM = JISBLOCK(4,ISBLK)
                  ISOFF = JISBLOCK(5,ISBLK)
                  NIA = MNSSOA(IASM,IATP)
                  NIB = MNSSOB(IBSM,IBTP)
*
                  IF(NIA*NIB.EQ.0) GOTO 9999
                  IF(JIRESTRICT.EQ.1.AND.
     &               (JASM.GT.IASM.OR.
     &               JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &               JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &            GOTO 9999
*. Are the two blocks connected by allowed excitation
              INTERACT = 0
              IF(MXEXC.EQ.1) THEN
                 IF((JICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (JICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((JICONSPA(IATP,LLATP).LE.1.AND.
     &               JICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (JICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (JICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
C
                 IF(INTERACT.EQ.0) GOTO 9999
*
                 ILOOPSB = ILOOPSB + 1
*              INSERT SOMETHING HERE !!!
*
*
              END IF
 9999         CONTINUE
              ILOOPSB = ILOOPSB + 0
10000         CONTINUE
*./\  End of loop over sigma blocks
 8765         CONTINUE
            END IF
*
            ILOOPSB = ILOOPSB * NPERM
*
            JCTCMP(ICBLK,3) = ILOOPSB
*           if c-block is interacting and at least one s-block then 
*           NYTOTB(ICBLK) > 0, else 0
            NYTOTB(ICBLK) = 
     &      JCTCMP(ICBLK,3)*JCTCMP(ICBLK,2)*JCTCMP(ICBLK,1) 
*
*
            ICSBTCHB(ICBLK) = ICSBTCHB(ICBLK) + NYTOTB(ICBLK)
*
*            WRITE(6,*)'NYTOTB(ICBLK) =',NYTOTB(ICBLK), myproc
*            WRITE(6,*)'ICSBTCHB(ICBLK) =',ICSBTCHB(ICBLK), myproc
            
 9000       CONTINUE
*. /\ End of loop over C blocks in Batch
*. End of loop over S blocks in batch
20000   CONTINUE
*. End of loop over batches of C blocks

*. we wanna find a distribution of the c-blocks among the nodes (incl. master)
*. as even as possible, here we go ! 
*
      IF(ICOMPUTE.GT.0.AND.IUSE_VER.EQ.2) THEN
        IF(NEXCL_MR.EQ.0) THEN
*
          CALL FIND_EVEN_DISTR(JCTCMP,MNJJJBLK,JCBLINF,MTEMPP,
     &         JMPINFA,NYTOTB)
*
        ELSE IF(NEXCL_MR.GT.0) THEN
*
          CALL FIND_EVEN_DISTR_EXM(JCTCMP,MNJJJBLK,JCBLINF,MTEMPP,
     &         JMPINFA,NYTOTB)
*
        END IF
      END IF
*
      RETURN 
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_SBATCH_DISTR(JCSBTCHB,NNJJJBLK,NTEMPP,IDISTLIST)
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*
*. input
* NNJJJBLK: total number of c-blocks
* NTEMPP: number of processes available for the distribution
* ICSBTCHB(NNJJJBLK): total length of each block
*. output
* IDISTLIST(NNJJJBLK): list containing the number of a node  
* that computes the block
*
      INTEGER NNJJJBLK, ICOUNTABLK
      INTEGER NTEST
      INTEGER*8 ITOTBLCKL,MMPINFA
      DIMENSION MMPINFA(NTEMPP,2)
      DIMENSION IDISTLIST(NNJJJBLK), JCSBTCHB(NNJJJBLK)
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: actual number of the block showing the LONGEST! :)
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*      DATA MXSIZE,ITEMPW1 /20000000000000000,20000000000000000/
*
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 1900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1900000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO I = 1, NTEMPP
       MMPINFA(I,1) = 0
       MMPINFA(I,2) = 0
      END DO
*. end of initializtion
      DO ICBL = 1, NNJJJBLK
*
        IF(JCSBTCHB(ICBL).NE.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
*           IF(MYNEW_ID.EQ.0)THEN
*             WRITE(6,*) 'This is a active block with length',
*     &ICBL,MCTCMP(ICBL,2)
*           END IF
           ITOTBLCKL = ITOTBLCKL + JCSBTCHB(ICBL)
        END IF
*        IF(MYNEW_ID.EQ.0)THEN
*           WRITE(6,*)'JCSBTCHB(ICBL) is',JCSBTCHB(ICBL)
*        END IF
*
      END DO
*
      IF(NTEST.GT.0) THEN
       IF(MYNEW_ID.EQ.0) THEN
         WRITE(6,*) '  total number of processes to distribute on:',
     &                 NTEMPP -1
         WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
         WRITE(6,*) '  overall active block length:', ITOTBLCKL
       END IF
      END IF
*
      IF(ICOUNTABLK.LT.NTEMPP-1) THEN
        IMINNP = MIN(ICOUNTABLK,NTEMPP)
        IF(NTEST.GT.0) THEN
          IF(MYNEW_ID.EQ.0) THEN
            WRITE(6,*) 'subroutine FIND_SBATCH_DISTR speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP+1,'procs'
          END IF
        END IF
        NTEMPP = IMINNP+1
      END IF
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
100   CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.ICOUNTABLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, NNJJJBLK
*
          ITEMPL = JCSBTCHB(II)
          ITEMPD = IDISTLIST(II)
          ITEMPN = II

*          WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
*          WRITE(6,*) 'MXSIZE:',MXSIZE
*          WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.LT.0) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
*                  WRITE(6,*) 'MXNUMB',MXNUMB
*                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
*        WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 2, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR,2)
*          WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MYNEW_ID = 0 !
          IAMTPROC = IPR
*          WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.2) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
*             WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
       IF(NTEST.GT.0) THEN
       WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
       IDISTLIST(MXNUMB)  = IAMGPROC - 1
*
*.      raising MMPINFA(proc,1) by 1
*.      adding on MMPINFA(proc,2) the weight of the new block
        MXSIZE = MXSZTMP
        MMPINFA(IAMGPROC,1) = MMPINFA(IAMGPROC,1) + 1
        MMPINFA(IAMGPROC,2) = MMPINFA(IAMGPROC,2) + MXSIZE
*
      ELSE
        GOTO 101
      END IF
*     /\ IRUN !!!
*
      GOTO 100
*
101   CONTINUE
*
      IF(NTEST.GT.0) THEN
        IF(MYNEW_ID.EQ.0) THEN
      WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF MYPROC:',
     &MYNEW_ID
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MMPINFA(ISTI,1),' blocks with a
     & total length of',MMPINFA(ISTI,2)
      END DO
*
*      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_EVEN_DISTR(MCTCMP,NNJJJBLK,MCBLINF,NTEMPP,
     &           MMPINFA,JYTOTB)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*
*. input
* NNJJJBLK: total number of c-blocks 
* MCBLINF(NNJJJBLK): length of c-block
* NTEMPP: number of processes available for the distribution
* MCTCMP: Blocks info array 
*   MCTCMP(*,1) : block contributes 1 or not 0
*   MCTCMP(*,2) : block length
*   MCTCMP(*,3) : count of interacts with s-blocks
*   MCTCMP(*,4) : node number that will get the block; default is -1
*   JYTOTB(*) : length of ACTIVE c-block * count of interacts with s-blocks 
*
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks 
      INTEGER NNJJJBLK, ICOUNTABLK 
      INTEGER NTEST 
      INTEGER*8 ITOTBLCKL
      DIMENSION MCTCMP(NNJJJBLK,4), MCBLINF(NNJJJBLK)
      DIMENSION MMPINFA(NTEMPP,2), JYTOTB(NNJJJBLK) 
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: actual number of the block showing the LONGEST! :)
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 1900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1900000000 
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO ICBL = 1, NNJJJBLK
*
        IF(MCTCMP(ICBL,1).NE.0) THEN 
           ICOUNTABLK = ICOUNTABLK + 1
*           IF(MYPROC.EQ.MASTER)THEN
*             WRITE(6,*) 'This is a active block with length',
*     &ICBL,MCTCMP(ICBL,2)
*           END IF
           ITOTBLCKL = ITOTBLCKL + JYTOTB(ICBL)
        END IF
*        WRITE(6,*)'MCTCMP(ICBL,4) is',MCTCMP(ICBL,4) 
*
      END DO 
*
      IF(NTEST.GT.0) THEN
       WRITE(6,*) '  total number of processes to distribute on:',NTEMPP
       WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
       WRITE(6,*) '  overall active block length:', ITOTBLCKL
      END IF 
*
      IF(ICOUNTABLK.LT.NTEMPP) THEN
        IMINNP = MIN(ICOUNTABLK,NTEMPP)
        IF(NTEST.GT.0) THEN
          WRITE(6,*) 'subroutine FIND_EVEN... speaking!'
          WRITE(6,*) 'number of active blocks lower than number of 
     & procs --> I only use',IMINNP,'procs'
        END IF
        NTEMPP = IMINNP
      END IF
*
*. starting the treausure quest for the ?optimal? c-block distribution 
*
100   CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.ICOUNTABLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, NNJJJBLK
*
          
          ITEMPL = JYTOTB(II)
          ITEMPD = MCTCMP(II,4)
          ITEMPN = II
*          WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
*          WRITE(6,*) 'MXSIZE:',MXSIZE
*          WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.LT.0) THEN
              IF(ITEMPL.LE.MXSIZE)THEN 
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
*                  WRITE(6,*) 'MXNUMB',MXNUMB
*                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
*        WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 1, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR,2)
*          WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MASTER !
          IAMTPROC = IPR 
*          WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC 
*             WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF   
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M incl.)
       IF(NTEST.GT.0) THEN
       WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
        MCTCMP(MXNUMB,4) = IAMGPROC - 1
*              
*.      raising MMPINFA(proc,1) by 1
*.      adding on MMPINFA(proc,2) the weight of the new block      
        MXSIZE = MXSZTMP
        MMPINFA(IAMGPROC,1) = MMPINFA(IAMGPROC,1) + 1
        MMPINFA(IAMGPROC,2) = MMPINFA(IAMGPROC,2) + MXSIZE
*
      ELSE 
        GOTO 101
      END IF
*     /\ IRUN !!!     
*
      GOTO 100
*
101   CONTINUE
* 
      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.MASTER) THEN
      WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF MYPROC:',
     &MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MMPINFA(ISTI,1),' blocks with a  
     & total length of',MMPINFA(ISTI,2)
      END DO
*
      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_EVEN_DISTR_EXM(MCTCMP,NNJJJBLK,MCBLINF,NTEMPP,
     &           MMPINFA,JYTOTB)

*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*
*. input
* NNJJJBLK: total number of c-blocks
* MCBLINF(NNJJJBLK): length of c-block
* NTEMPP: number of processes available for the distribution
* MCTCMP: Blocks info array
*   MCTCMP(*,1) : block contributes 1 or not 0
*   MCTCMP(*,2) : block length
*   MCTCMP(*,3) : count of interacts with s-blocks
*   MCTCMP(*,4) : length of ACTIVE c-block * count of interacts with s-blocks
*   MCTCMP(*,5) : node number that will get the block; default is -1
*
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER NNJJJBLK, ICOUNTABLK
      INTEGER NTEST
      INTEGER*8 ITOTBLCKL
      DIMENSION MCTCMP(NNJJJBLK,4), MCBLINF(NNJJJBLK)
      DIMENSION MMPINFA(NTEMPP,2), JYTOTB(NNJJJBLK)
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: actual number of the block showing the LONGEST! :)
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 1900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1900000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO ICBL = 1, NNJJJBLK
*
        IF(MCTCMP(ICBL,1).NE.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
*           IF(MYPROC.EQ.MASTER)THEN
*             WRITE(6,*) 'This is a active block with length',
*     &ICBL,MCTCMP(ICBL,2)
*           END IF
           ITOTBLCKL = ITOTBLCKL + JYTOTB(ICBL)
        END IF
*        WRITE(6,*)'MCTCMP(ICBL,5) is',MCTCMP(ICBL,5)
*
      END DO
*
*      IF(NTEST.GT.0) THEN
       IF(MYPROC.EQ.MASTER) THEN
        WRITE(6,'(3X,A,1X,I4)') 'total number of processes to 
     &distribute on:',NTEMPP-1
        WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
        WRITE(6,*) '  overall active block length:', ITOTBLCKL
       END IF
*      END IF
*
      IF(ICOUNTABLK.LT.NTEMPP-1) THEN
        IMINNP = MIN(ICOUNTABLK,NTEMPP)
*        IF(NTEST.GT.0) THEN
          IF(MYPROC.EQ.MASTER) THEN
            WRITE(6,*) 'subroutine FIND_EVEN...EXM speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP+1,'procs'
         END IF
*        END IF
        NTEMPP = IMINNP+1
      END IF
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
100   CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.ICOUNTABLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, NNJJJBLK
*
          ITEMPL = JYTOTB(II)
          ITEMPD = MCTCMP(II,4)
          ITEMPN = II

*          WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
*          WRITE(6,*) 'MXSIZE:',MXSIZE
*          WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.LT.0) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
*                  WRITE(6,*) 'MXNUMB',MXNUMB
*                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
*        WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 2, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR,2)
*          WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MASTER !
          IAMTPROC = IPR
*          WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.2) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
*             WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
       IF(NTEST.GT.0) THEN
       WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
        MCTCMP(MXNUMB,4) = IAMGPROC - 1
*
*.      raising MMPINFA(proc,1) by 1
*.      adding on MMPINFA(proc,2) the weight of the new block
        MXSIZE = MXSZTMP
        MMPINFA(IAMGPROC,1) = MMPINFA(IAMGPROC,1) + 1
        MMPINFA(IAMGPROC,2) = MMPINFA(IAMGPROC,2) + MXSIZE
*
      ELSE
        GOTO 101
      END IF
*     /\ IRUN !!!
*
      GOTO 100
*
101   CONTINUE
*
*      IF(NTEST.GT.0) THEN
        IF(MYPROC.EQ.MASTER) THEN
      WRITE(6,'(3X,A,I4)')'SUMMATION OF EVEN DISTRIBUTION OF MYPROC:',
     &MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MMPINFA(ISTI,1),' blocks with a
     & total length of',MMPINFA(ISTI,2)
      END DO
*
      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
*      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS_2(LUIN,LBLK,BLK_A,SEGMNT,
     &           JCBLINF,NTOTBLK,IAMCOMM)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A, length is stored in JCBLINF
*  master version
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*. Output
      DIMENSION BLK_A(*), JCBLINF(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      IF(MYNEW_ID.EQ.0)THEN
        CALL REWINE(LUIN,LBLK)
*
        IBLK = 0
        NBLK_A = 0
        NBLK = 0
 1000   CONTINUE
          IBLK = IBLK + 1
          IF(LBLK .GT. 0 ) THEN
            LBL = LBLK
          ELSE IF ( LBLK .EQ. 0 ) THEN
            READ(LUIN) LBL
          ELSE IF  (LBLK .LT. 0 ) THEN
            CALL IFRMDS(LBL,1,-1,LUIN)
          END IF
          IF( LBL .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = LBL
            ELSE
              KBLK = -1
            END IF
            NO_ZEROING = 1
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
            IF(IMZERO.EQ.0) THEN
              NBLK_A = NBLK_A + 1
              BLK_A(IBLK) = 1.0D0
              JCBLINF(IBLK) = LBL
            ELSE
              BLK_A(IBLK) = 0.0D0
              JCBLINF(IBLK) = LBL
            END IF
          END IF
      
        IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
        NBLK =  IBLK-1
*
      END IF
*.   /\ MASTER work!
      NTEST = 0
*test this sk
      IF(NPARWAY.NE.0)THEN
*      DO ISS = 1, NTOTBLK
         CALL MPI_BCAST(BLK_A,NTOTBLK,MPI_REAL8,0,IAMCOMM,IERR)
*         CALL MPI_BCAST(BLK_A(ISS),1,MPI_REAL8,MASTER,IAMCOMM,IERR)
         CALL MPI_BCAST(JCBLINF,NTOTBLK,MPI_INTEGER,0,IAMCOMM,IERR)
*         CALL MPI_BCAST(JCBLINF(ISS),1,MPI_INTEGER,MASTER,IAMCOMM,IERR)
*      END DO
      END IF
* -sk
      IF(NTEST.GE.1) THEN
      IF(MYNEW_ID.EQ.MASTER)THEN
        WRITE(6,'(A,I8,I8)')
     &  ' FIND_A...Number of total and active Blocks',NBLK,NBLK_A,
     &  MYNEW_ID
      END IF
      END IF
      IF(NTEST.GE.100) THEN
      IF(MYNEW_ID.LT.2)THEN
        DO IS = 1,NTOTBLK
          WRITE(6,*) 'BLOCK LENGTH OF BLOCK',IS,'=',JCBLINF(IS),MYNEW_ID
        END DO
        WRITE(6,*) ' Active blocks ',MYNEW_ID
*        CALL WRTMAT(BLK_A,1,NBLK,1,NBLK)
        CALL WRTMAT(BLK_A,1,NTOTBLK,1,NTOTBLK)
      END IF
      END IF
*
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS_MR(LUIN,LBLK,BLK_A,SEGMNT,IFHBLK)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*  master version
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*      INTEGER ISTBLK
*
      CALL REWINE(LUIN,LBLK)
*
      IBLK = 0
*      IBLK = ISTBLK
      NBLK_A = 0
      NBLK = 0
**      write(6,*)'ISTBLK',ISTBLK
*. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
*        write(6,*)'IBLK =',IBLK
        IF(IBLK.GT.IFHBLK) GOTO 1001
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
*      write(6,*)'lbl of IBLK =',lbl, iblk
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
1001  CONTINUE
      NBLK =  IBLK - 1
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,'(A,I8,I8)')
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active blocks (master)'
        CALL WRTMAT(BLK_A,1,NBLK,1,NBLK)
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ADDSB(A,B,NDIM)
*
* add recieved SB vector to current SB vector
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*),B(*)
*
      DO 100 I = 1,NDIM
        A(I) = A(I) + B(I)
  100 CONTINUE
*
      RETURN
      END
*
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DIST_NODES(NTOTALND,NACTROOT,JDISTROOT)
*
* find an even distribution for the sum NTOTALND of the nodes 
* w.r.t. the number of roots NACTROOT that have to be calculated
* distribution is stored in array JDISTROOT 
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER NTOTALND,NACTROOT,ISTNODE,ISTIADD 
      DIMENSION JDISTROOT(NTOTALND)
*
      ISTNODE = 1
      ISTIADD = 1
 1    CONTINUE
      IF(ISTIADD.LE.NACTROOT)THEN
*
        JDISTROOT(ISTNODE) = ISTIADD
        ISTIADD = ISTIADD + 1
        ISTNODE = ISTNODE + 1
        GOTO 1
*        
      ELSE IF(ISTIADD.GT.NACTROOT) THEN
*
        ISTIADD = NACTROOT
        GOTO 2
*
      END IF
 2    CONTINUE
*
      ISTNODE = NTOTALND
 3    CONTINUE
      ISTIADD = NACTROOT
 4    CONTINUE
      IF(ISTIADD.NE.0.AND.ISTNODE.NE.NACTROOT) THEN
        JDISTROOT(ISTNODE) = ISTIADD
        ISTIADD = ISTIADD - 1
        ISTNODE = ISTNODE - 1
        GOTO 4
      ELSE IF(ISTIADD.EQ.0.AND.ISTNODE.GT.NACTROOT) THEN
        GOTO 3 
      ELSE IF(ISTIADD.NE.0.AND.ISTNODE.EQ.NACTROOT) THEN
        GOTO 5
      ELSE IF(ISTIADD.EQ.0.AND.ISTNODE.EQ.NACTROOT) THEN
        GOTO 5
      END IF
 5    CONTINUE
*    
      RETURN
      END
*
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE CBATCH_LOOP(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &           ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &           XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &           MXTSOB,SOMESCR,SOMEH,XIJILS,
     &           I1,XI1S,I2,XI2S,IDOH2,MXPOBS,ISTRFL,
     &           PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &           I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &           ICONSPA,ICONSPB,SCLFAC,
     &           LUCBLK,IPERTOP,IH0INSPC,IH0SPC,
     &           ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &           I_RES_AB,IUSE_PA,CJPA,SIPA,
     &           ILENSUM,ICBLINF,NJJJBLK,IRMBNP,ICTCMP,
     &           SSCLFAC,ICSBTCHB,ILOOPC,ILOOPCS,ICALCC,ISKIPR,
     &           IDISTLIST,MXEXC)
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
* NEWCOMM_PROC: temporary number of processes
* IMPINFA(NEWCOMM_PROC,1): number of blocks to calculate fo a specific proc
* IMPINFA(NEWCOMM_PROC,2): total length of all blocks
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*. Specific input
      INTEGER ISBLOCK(8,*)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*. Zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*), MXEXC
C    &                  LUCBLK,IPERTOP,IH0INSPC,IH0SPC)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*),SSCLFAC(*),ICTCMP(NJJJBLK,4)
      INTEGER*8 IMPINFA,MYTOTB,ICSBTCHB
      DIMENSION IMPINFA(NEWCOMM_PROC,2),MYTOTB(NJJJBLK)
      DIMENSION ICSBTCHB(NJJJBLK),IDISTLIST(NJJJBLK)
* new scratch
      INTEGER ILENSUM
      DOUBLE PRECISION CPUSIGXXX1,CPUSIGXXX2,WALLSIGXXX1,WALLSIGXXX2
      DOUBLE PRECISION CPUSIGXXXX1,CPUSIGXXXX2
      DOUBLE PRECISION WALLSIGXXXX1,WALLSIGXXXX2
      DOUBLE PRECISION CPUSIGXX1,CPUSIGXX2,WALLSIGXX1,WALLSIGXX2
      DOUBLE PRECISION CPUSIGX1,CPUSIGX2,WALLSIGX1,WALLSIGX2
      CHARACTER SECTID*12, CPUTIDXXX*12, WALLTIDXXX*12
      CHARACTER CPUTIDXXXX*12, WALLTIDXXXX*12
      CHARACTER CPUTIDX*12, WALLTIDX*12
      CHARACTER CPUTIDXX*12, WALLTIDXX*12

      INTEGER NJJJBLK,ILOOPC,ILOOPCS,ICALCC,ISKIPR
      DIMENSION ICBLINF(NJJJBLK),IRMBNP(NEWCOMM_PROC)
*      INTEGER VERSION, SUBVERSION
#include "comjep.inc"
*.
      REAL * 8 INPROD
*
*        write(6,*)'ILOOPC = ',ILOOPC, myproc
*        write(6,*)'ISKIPR = ',ISKIPR, myproc
*        write(6,*)'SCLFAC(37) = ',SCLFAC(4), myproc
*        write(6,*)'LCBLOCK(1) = ',LCBLOCK(1), myproc
*
*
*. Find batches of C - strings
*. Should be changed to REPART or ??
*. Not neccesary, all zero blocks detected otherwise
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &              LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &              I1CBLOCK,ICBLOCK,0,ITTSS_ORD)
*

*
        DO 19999 IDOTWL = 1, ILOOPC
*
           JCBATCH = ISKIPR + IDOTWL
*. Read C blocks into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
*        IF(mynew_id.EQ.MASTER) THEN
*          WRITE(6,*)'NJBLOCK (mynew_id) = ', NJBLOCK, mynew_id
*        END IF
        DO JJCBLOCK = 1, NJBLOCK
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
*. calculating this block?
          IF(IRUNPA.GT.0) THEN
            IF(MYNEW_ID.NE.0.AND.NPARWAY.GT.3)THEN
              IF((IDISTLIST(JBLOCK).NE.MYNEW_ID).OR.
     &            (IDISTLIST(JBLOCK).LT.0)) THEN
                GOTO 21
              END IF
            END IF
          END IF
          INTERACT = 0
          IF(SCLFAC(JBLOCK).EQ. 1.0D0) THEN
C?          WRITE(6,*) ' Nonvanishing block ', JBLOCK, mynew_id
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
            DO JSBLOCK = 1, NSBLOCK
            IDENT = 0
            IF(ISBLOCK(1,JSBLOCK).GT.0) THEN
              IATP = ISBLOCK(1,JSBLOCK)
              IBTP = ISBLOCK(2,JSBLOCK)
              IASM = ISBLOCK(3,JSBLOCK)
              IBSM = ISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
              IF(MXEXC.EQ.1) THEN
                 IF((ICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (ICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &               ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
*. Are they identical ?
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
          IF(INTERACT.EQ.1) THEN
*            WRITE(6,*)'i am here to go for GSTTBL', mynew_id
            IF((MYNEW_ID.NE.0.AND.NPARWAY.GT.3).OR.(NPARWAY.LT.4))THEN
**              WRITE(6,*)'i am here to go for GSTTBL', mynew_id
**              WRITE(6,*)'JOFF =',JOFF
              CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                    NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                    PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))

**              CALL WRTMAT(CB(JOFF),NSSOB(JBSM,JBTP),
**     &                    NSSOA(JASM,JATP),NSSOB(JBSM,JBTP),
**     &                    NSSOA(JASM,JATP))

            END IF
**            SCLFAC(JBLOCK) = 1.0D0
*. Note in GSTTBL : ICOOSC only used for CI vectors in core,
          ELSE
*. not relevant
C             WRITE(6,*) ' TTSS for C block skipped  '
C             CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            IF((MYNEW_ID.NE.0.AND.NPARWAY.GT.3).OR.(NPARWAY.LT.4))THEN
              CALL IFRMDS(LBL,-1,1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
            END IF
            SCLFAC(JBLOCK) = 0.0D0
          END IF
*
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            END IF
          END IF
*
 21     CONTINUE
*
*
        END DO
*
 33   CONTINUE
*
          DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK
*
*. calculating this block?
          IF(IRUNPA.GT.0) THEN
            IF(NPARWAY.GT.3) THEN
              IF(IDISTLIST(ICBLK).NE.MYNEW_ID) THEN
                GOTO 8888
              END IF
            END IF
          END IF
            JATP = ICBLOCK(1,ICBLK)
            JBTP = ICBLOCK(2,ICBLK)
            JASM = ICBLOCK(3,ICBLK)
            JBSM = ICBLOCK(4,ICBLK)
            ICOFF = ICBLOCK(5,ICBLK)
            NJA = NSSOA(JASM,JATP)
            NJB = NSSOB(JBSM,JBTP)
*
            IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
              CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                    LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
              DO 8765 IPERM = NPERM,1, -1
                LLASM = LASM(IPERM)
                LLBSM = LBSM(IPERM)
                LLATP = LATP(IPERM)
                LLBTP = LBTP(IPERM)
                NLLA = NSSOA(LLASM,LLATP)
                NLLB = NSSOB(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
                IF(IPERM.EQ.1) THEN
                  IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
*. Diagonal blocks, Transposing corresponds to scaling
                    IF(PS.EQ.-1.0D0) THEN
                      CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                    END IF
                  ELSE
*. offdiagonal blocks, explicit transposing
                    CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                    CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
                    if (LSGN(IPERM).eq.-1) then
                      write(6,*) 'LROW,LCOL undef. in sblocks'
                      write(6,*) 'and LSGN(IPERM).eq.-1'
                      Call Abend2('Quitting.')
CTF                   CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                    end if
                  END IF
                END IF
*
                DO 10000 ISBLK = 1, NSBLOCK
                IF(ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                  IATP = ISBLOCK(1,ISBLK)
                  IBTP = ISBLOCK(2,ISBLK)
                  IASM = ISBLOCK(3,ISBLK)
                  IBSM = ISBLOCK(4,ISBLK)
                  ISOFF = ISBLOCK(5,ISBLK)
                  NIA = NSSOA(IASM,IATP)
                  NIB = NSSOB(IBSM,IBTP)
*
                  IF(NIA*NIB.EQ.0) GOTO 10000
                  IF(IRESTRICT.EQ.1.AND.
     &               (JASM.GT.IASM.OR.
     &               JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &               JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &            GOTO 10000
*. Are the two blocks connected by allowed excitation
              INTERACT = 0
              IF(MXEXC.EQ.1) THEN
                 IF((ICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (ICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &               ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
C
                 IF(INTERACT.EQ.0) GOTO 10000

                 IF(NTEST.GE.100) THEN
                   WRITE(6,*) ' Next s block in batch : '
                   write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                   write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                 END IF
*
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            ((LLBSM.GT.LLASM).OR.
     &            (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &            GOTO 8764
*
                IF(NTEST.GE.60) THEN
                  WRITE(6,*) ' RSSBCB will be called for '
                  WRITE(6,*) ' Sigma block : '
                  WRITE(6,*) ' ISOFF ', ISOFF
                  WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' C     block : '
                  WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(6,*) ' ICOFF ', ICOFF
                  WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                END IF
*
                IF(IRESTRICT.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     XFAC = 0.5D0*SCLFAC(ICBLK)
                ELSE
                     XFAC = SCLFAC(ICBLK)
                END IF
*. Form of operator in action
                IF(IPERTOP.NE.0) THEN
*. Not exact Hamiltonian in use
                   IPTSPC = IH0SPC(IATP,IBTP)
                   JPTSPC = IH0SPC(JATP,JBTP)
                   IJOP   = IH0INSPC(IPTSPC)
*
C?                 WRITE(6,*) ' SBLOCKS:IPTSPC,JPTSPC,IJOP',
C?   &                         IPTSPC,JPTSPC,IJOP
*
                ELSE
                   IPTSPC = -1
                   JPTSPC = -1
                   IJOP   = -1
                END IF
*
                CALL RSSBCB2(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
     &               IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &               IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &               MXTSOB,SOMESCR,SOMEH,XIJILS)
C?            IF(IATP.EQ.1.AND.IBTP.EQ.1) THEN
C?              WRITE(6,*) ' First element in updates S11 :',SB(1)
C?            END IF
C?            IF(JATP.EQ.1.AND.JBTP.EQ.1) THEN
C?              WRITE(6,*) ' First element of C ', CB(1)
C?            END IF
 8764         CONTINUE
              END IF
10000         CONTINUE
*./\  End of loop over sigma blocks
 8765         CONTINUE
            END IF
 8888        CONTINUE
 9000       CONTINUE
*. /\ End of loop over C blocks in Batch
*. End of loop over S blocks in batch
19999   CONTINUE
*. End of loop over batches of C blocks
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MSSNDBLCK(LUIN,LUOUT,SEGMNT,JCOMM,LBLK,IDISTLIST,
     &           NMBLOCK)
*
* Read vector (master) from disc file luin to array segmnt
* and send to node IDISTLIST(block) the vector where it is 
* written in file luout
*
* LBLK DEFINES STRUCTURE OF FILE
* COMM IS THE COMMUNICATOR USED FOR MPI-COMMUNICATION
* --> DEFAULT: MPI_COMM_WORLD
* YOUR CALL SHOULD LOOK LIKE: CALL MSS...(LUIN,...,MPI_COMM_WORLD,LBLK)
*
*
* Type of file LUOUT is inherited from LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*),IDISTLIST(NMBLOCK)
      INTEGER ISENDN,NROUND
C
*
      IREW = 1
*
      IF( IREW .NE. 0 ) THEN
        IF (MYNEW_ID.EQ.0) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (MYNEW_ID.NE.0) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      ISENDN = 0
      NROUND = 0
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
      NROUND = NROUND + 1
      IF(NROUND.LE.NMBLOCK+1) THEN
        IF(NROUND.EQ.NMBLOCK+1) THEN
          GOTO 9000
        ELSE
          ISENDN = IDISTLIST(NROUND)
*          if (myproc.eq.master) WRITE(6,*)'NROUND=',NROUND
*          if (myproc.eq.master) WRITE(6,*)'I..D(NROUND)=ISENDN =',ISENDN
*          if (myproc.eq.master) WRITE(6,*)'LBLK =',LBLK
*          if (myproc.eq.master) WRITE(6,*)'LUIN =',LUIN
        END IF
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (MYNEW_ID.EQ.0.and.ISENDN.GT.0) then
              call mpi_send(LBL,1,MPI_INTEGER,ISENDN,60,
     &                      JCOMM,ISTAT,IERR)
            end if
            if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
              call mpi_recv(LBL,1,MPI_INTEGER,0,60,
     &                      JCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (MYNEW_ID.EQ.0) then
            READ(LUIN) LBL
          if (MYNEW_ID.EQ.0.and.ISENDN.GT.0) then
            call mpi_send(LBL,1,MPI_INTEGER,ISENDN,60,
     &                    JCOMM,ISTAT,IERR)
          end if
          end if


          if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
            call mpi_recv(LBL,1,MPI_INTEGER,0,60,
     &                    JCOMM,ISTAT,IERR)
            WRITE(LUOUT) LBL
        end if
        if (MYNEW_ID.EQ.0) then
        write(6,*) 'SR MSSNDBLCK:  LBLK = ', LBLK
        end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (MYNEW_ID.EQ.0) then
*            write(6,*)'master reads in block',NROUND
            CALL IFRMDS(LBL,1,-1,LUIN)
          if (MYNEW_ID.EQ.0.and.ISENDN.GT.0) then
            call mpi_send(LBL,1,MPI_INTEGER,ISENDN,60,
     &                    JCOMM,ISTAT,IERR)
*            write(6,*)'master has sent block',NROUND
          end if
          end if
*
*
          if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
            call mpi_recv(LBL,1,MPI_INTEGER,0,60,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS(LBL,1,-1,LUOUT)
C?        write(6,*) ' COPVCD LBL : ', LBL
          end if
*
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (MYNEW_ID.EQ.0) then
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
            if (ISENDN.GT.0) then
              call mpi_send(IMZERO,1,MPI_INTEGER,ISENDN,61,
     &             JCOMM,ISTAT,IERR)
              call mpi_send(IAMPACK,1,MPI_INTEGER,ISENDN,62,
     &             JCOMM,ISTAT,IERR)
              call mpi_send(SEGMNT,LBL,MPI_REAL8,ISENDN,63,
     &             JCOMM,ISTAT,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
             call mpi_recv(IMZERO,1,MPI_INTEGER,0,61,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,MPI_INTEGER,0,62,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,0,63,
     &                     JCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
                CALL TODSC(SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
9000  CONTINUE
      IF (MYNEW_ID.NE.0) THEN
        CALL ITODS(-1,1,-1,LUOUT)
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MSBCVC_MPI(LUIN,LUOUT,SEGMNT,COMM,LBLK)
*
* Read vector (master) from disc file luin to array segmnt
* and broadcast the vector to all nodes
* where vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
*
* LBLK DEFINES STRUCTURE OF FILE
* COMM IS THE COMMUNICATOR USED FOR MPI-COMMUNICATION
* --> DEFAULT: MPI_COMM_WORLD
* YOUR CALL SHOULD LOOK LIKE: CALL MSB...(LUIN,...,MPI_COMM_WORLD,LBLK)
*
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*)
C
*
      IREW = 0
*
      IF( IREW .NE. 0 ) THEN
        IF (MYPROC.EQ.MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (MYPROC.NE.MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      ihere = 0
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
              call mpi_bcast(LBL,1,MPI_INTEGER,MASTER,
     &                      COMM,IERR)
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (myproc.eq.master) then
            READ(LUIN) LBL
          end if
            call mpi_bcast(LBL,1,MPI_INTEGER,MASTER,
     &                    COMM,IERR)
           if (myproc.ne.master) then
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
           end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (myproc.eq.master) then
            CALL IFRMDS(LBL,1,-1,LUIN)
          end if
            call mpi_bcast(LBL,1,MPI_INTEGER,MASTER,
     &                    COMM,IERR)
*
*
          if (myproc.ne.master) then
            CALL ITODS (LBL,1,-1,LUOUT)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (myproc.eq.master) then
            ihere = ihere + 1
**            write(6,*)'i am here for the x time', ihere
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
          end if
          call mpi_bcast(IMZERO,1,MPI_INTEGER,MASTER,
     &                  COMM,IERR)
          call mpi_bcast(IAMPACK,1,MPI_INTEGER,MASTER,
     &                  COMM,IERR)
          call mpi_bcast(SEGMNT,LBL,MPI_REAL8,MASTER,
     &                  COMM,IERR)
*
**          if (myproc.eq.master) then
**            write(6,*)'lbl,IMZERO,IAMPACK = ',lbl,IMZERO,IAMPACK
**            write(6,*)'this is segment(LBL)', SEGMNT(LBL)
**          end if
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (myproc.ne.master) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (myproc.ne.master) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MSCPCFF(SEGMNT,LUCFF,NJBATCH,JRESTC,JLOOPC,JLOOPCS,
     &                   NYCOMM,JCBLOCK,NNBLOCKS,LBLK)
*
* Read coff. blocks (master) from disc file lucff to array segmnt
* and send to node X
* where vector is written in file lucff
*
* MPI adaption         Stefan Knecht, Aug. 17., 2006
* MSDTAG: send-tag; NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*), JCBLOCK(*), NNBLOCKS(NEWCOMM_PROC)
*.scratch
      INTEGER MTAG,NBLCK,JMBLOCK,JSBLOCK,JNCBLOCK,IGETRNR,ICHNBLCK
C
*
      IREW = 1
      JMBLOCK = 0
      JSBLOCK = 0
      NBLCK = 0
      IGETRNR = 0
      NODEX = 0
      LBLK = -1
*
      IF( IREW .NE. 0 ) THEN
        IF (MYNEW_ID.EQ.0) THEN
          CALL REWINE( LUCFF ,LBLK)
        ELSE IF (MYNEW_ID.NE.0) THEN
          CALL REWINE( LUCFF ,LBLK)
        END IF
      END IF
      MTAG = 100
*
      IF(MYNEW_ID.EQ.0) THEN
*
        DO ISK = 1, JLOOPC
*
           JCBATCH = ISK
           JMBLOCK = JMBLOCK + JCBLOCK(JCBATCH)
*
        END DO
*

*
        IF(NEWCOMM_PROC.LT.NJBATCH) THEN
**
          DO IGTB = 1, NEWCOMM_PROC -1
             JNCBLOCK = 0
             IF(IGTB.GT.(NEWCOMM_PROC-JRESTC)) THEN
               IGETRNR = IABS(JRESTC - NEWCOMM_PROC + IGTB)
               ISKIPR = (IGTB*JLOOPCS) + IGETRNR
               DO ISK = 1, JLOOPC+1
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
               END DO
               NNBLOCKS(IGTB) = JNCBLOCK
             ELSE IF(IGTB.EQ.(NEWCOMM_PROC-JRESTC)) THEN
               ISKIPR = (IGTB*JLOOPCS)
               DO ISK = 1, JLOOPC+1
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
               END DO
               NNBLOCKS(IGTB) = JNCBLOCK
             ELSE IF(IGTB.LT.(NEWCOMM_PROC-JRESTC)) THEN
               ISKIPR = (IGTB*JLOOPCS)
               DO ISK = 1, JLOOPC
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
               END DO
               NNBLOCKS(IGTB) = JNCBLOCK
             END IF
          END DO
**
        ELSE IF(NEWCOMM_PROC.EQ.NJBATCH) THEN
*
          DO IGTB = 1, NEWCOMM_PROC -1
             JNCBLOCK = 0
             ISKIPR = (IGTB*JLOOPCS)
             DO ISK = 1, JLOOPC
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
             END DO
               NNBLOCKS(IGTB) = JNCBLOCK
          END DO
**
        ELSE IF(NEWCOMM_PROC.GT.NJBATCH) THEN
            DO IGTB = 1, (NJBATCH - 1)
            JNCBLOCK = 0
            ISKIPR = (IGTB*JLOOPCS)
            DO ISK = 1, JLOOPC
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
             END DO
               NNBLOCKS(IGTB) = JNCBLOCK
          END DO
        END IF
**
*
      ELSE IF (MYNEW_ID.NE.0) THEN
*
        ISKIPR = 0
* how many records do i have to skip ?
*
        IF(NEWCOMM_PROC.LT.NJBATCH) THEN
*
          IF(MYNEW_ID.GT.(NEWCOMM_PROC-JRESTC)) THEN
*
            IGETRNR = IABS(JRESTC - NEWCOMM_PROC + MYNEW_ID)
            ISKIPR = (MYNEW_ID*JLOOPCS) + IGETRNR
*
          ELSE IF(MYNEW_ID.LE.(NEWCOMM_PROC-JRESTC)) THEN
*
            ISKIPR = (MYNEW_ID*JLOOPCS)
*
          END IF
        ELSE IF(NEWCOMM_PROC.GE.NJBATCH) THEN
*
            ISKIPR = (MYNEW_ID*JLOOPCS)
*
        END IF
*
        DO ISK = 1, JLOOPC
          JCBATCH = ISKIPR + ISK
          JSBLOCK = JSBLOCK + JCBLOCK(JCBATCH)
        END DO
*
      END IF
C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (MYNEW_ID.eq.0.and.NEWCOMM_PROC.gt.1) then
              if(nodex.gt.0) then
                call mpi_send(LBL,1,MPI_INTEGER,NODEX,MTAG-2,
     &                        NYCOMM,IERR)
              end if
            end if
            if (MYNEW_ID.ne.0) then
              call mpi_recv(LBL,1,MPI_INTEGER,0,MTAG-2,
     &                      NYCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (MYNEW_ID.eq.0) then
            READ(LUCFF) LBL
            if(nodex.gt.0.and.NEWCOMM_PROC.gt.1) then
              call mpi_send(LBL,1,MPI_INTEGER,NODEX,MTAG-2,
     &                      NYCOMM,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.ne.0) then
            call mpi_recv(LBL,1,MPI_INTEGER,0,MTAG-2,
     &                    NYCOMM,ISTAT,IERR)
            WRITE(LUCFF) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
            end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (MYNEW_ID.eq.0) then
            CALL IFRMDS(LBL,1,-1,LUCFF)
            if(nodex.gt.0.and.NEWCOMM_PROC.gt.1) then
              call mpi_send(LBL,1,MPI_INTEGER,NODEX,MTAG-2,
     &                      NYCOMM,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.ne.0) then
            call mpi_recv(LBL,1,MPI_INTEGER,MASTER,MTAG-2,
     &                    NYCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUCFF)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (MYNEW_ID.eq.0) then
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUCFF,IMZERO,IAMPACK,
     &                 NO_ZEROING)
            if(nodex.gt.0.and.NEWCOMM_PROC.gt.1) then
              call mpi_send(IMZERO,1,MPI_INTEGER,NODEX,MTAG+1,
     &                      NYCOMM,IERR)
              call mpi_send(IAMPACK,1,MPI_INTEGER,NODEX,MTAG+2,
     &                      NYCOMM,IERR)
              call mpi_send(SEGMNT,LBL,MPI_REAL8,NODEX,MTAG+3,
     &                      NYCOMM,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.ne.0) then
             call mpi_recv(IMZERO,1,MPI_INTEGER,0,MTAG+1,
     &                     NYCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,MPI_INTEGER,0,MTAG+2,
     &                     NYCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,0,MTAG+3,
     &                     NYCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (MYNEW_ID.ne.0) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUCFF)
              end if
            ELSE
              if (MYNEW_ID.ne.0) then
                CALL TODSCP(SEGMNT,LBL,KBLK,LUCFF)
              end if
            END IF
          ELSE
            if (MYNEW_ID.ne.0) then
              CALL ZERORC(LBL,LUCFF,IAMPACK)
            end if
          END IF
        END IF
*
      NBLCK = NBLCK + 1
**      if(myproc.eq.master) then
**        write(6,*)'NBLCK (myproc) =', NBLCK,myproc
**      end if
      IF(MYNEW_ID.EQ.0.and.NEWCOMM_PROC.gt.1) THEN
*
        IF(NBLCK.LT.JMBLOCK) THEN
*           write(6,*)'i am here?', myproc
           GOTO 999
        END IF
*
        IF(NBLCK.EQ.JMBLOCK) THEN
*           write(6,*)'start with node1 at NBLCK=',NBLCK
           NODEX = 1
           GOTO 999
        END IF
*
        IF(NODEX.GT.0) THEN
          ICHNBLCK = 0
          ICHNBLCK = JMBLOCK
          DO IFCH = 1, NODEX
             ICHNBLCK = ICHNBLCK + NNBLOCKS(IFCH)
          END DO
**          write(6,*)'ICHNBLCK=',ICHNBLCK
          IF(NBLCK.EQ.ICHNBLCK) THEN
*.? is the complete number of nodes relevant? I assume, it is NOT!
            NODEX = NODEX + 1
*            write(6,*)'master comm. now with nodex',NODEX
*.?
            IF(NEWCOMM_PROC.LE.NJBATCH) THEN
              IF(NODEX.LE.(NEWCOMM_PROC-1)) THEN
**                write(6,*)'master comm. now with nodex',NODEX
                GOTO 999
              ELSE IF(NODEX.GT.(NEWCOMM_PROC-1)) THEN
**                write(6,*)'master comm. NOT with nodex',NODEX
                GOTO 1001
              END IF
*
            ELSE IF(NEWCOMM_PROC.GT.NJBATCH) THEN
              IF(NODEX.LE.(NJBATCH-1)) THEN
**                write(6,*)'master comm. now with nodex',NODEX
                GOTO 999
              ELSE IF(NODEX.GT.(NJBATCH-1)) THEN
**                write(6,*)'master comm. NOT with nodex',NODEX
                GOTO 1001
              END IF
            END IF
          END IF
*.        /\ end if of (NBLCK.EQ.(JMBLOCK + NNBLOCKS(NODEX))
        END IF
      ELSE IF (MYNEW_ID.NE.0) THEN
        IF(NBLCK.EQ.JSBLOCK) THEN
          CALL ITODS (-1,1,-1,LUCFF)
          GOTO 1001
        END IF
      END IF
*
999   CONTINUE
*      write(6,*)'after 999 (myproc)', myproc
*      write(6,*)'LBL, LBLK',LBL,LBLK
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
1001  CONTINUE
**      write(6,*)'everything done', myproc

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_MS(LUIN,LUOUT,SEGMNT,NOFSLV,JCOMM,LBLK)
*
* Read vector (master) from disc file luin to array segmnt
* and send to node nofslv using communicator jcomm
* where vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
* NOFSLV: process-# of the slave
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*)
C
*
      IREW = 0
      IF(LBLK.NE.-1) THEN
        LBLK = -1
      END IF
*
      IF( IREW .NE. 0 ) THEN
        IF (MYPROC.EQ.MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (MYPROC.NE.MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      MYNSLV = NOFSLV

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (myproc.eq.master) then
              call mpi_send(LBL,1,MPI_INTEGER,MYNSLV,80,
     &                      JCOMM,IERR)
            end if
            if (myproc.ne.master) then
              call mpi_recv(LBL,1,MPI_INTEGER,MASTER,80,
     &                      JCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (myproc.eq.master) then
            READ(LUIN) LBL
            call mpi_send(LBL,1,MPI_INTEGER,MYNSLV,80,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.ne.master) then
            call mpi_recv(LBL,1,MPI_INTEGER,MASTER,80,
     &                    JCOMM,ISTAT,IERR)
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
            end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (myproc.eq.master) then
            CALL IFRMDS(LBL,1,-1,LUIN)
            call mpi_send(LBL,1,MPI_INTEGER,MYNSLV,80,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.ne.master) then
            call mpi_recv(LBL,1,MPI_INTEGER,MASTER,80,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUOUT)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (myproc.eq.master) then
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          call mpi_send(IMZERO,1,MPI_INTEGER,MYNSLV,81,
     &                  JCOMM,IERR)
          call mpi_send(IAMPACK,1,MPI_INTEGER,MYNSLV,82,
     &                  JCOMM,IERR)
          call mpi_send(SEGMNT,LBL,MPI_REAL8,MYNSLV,83,
     &                  JCOMM,IERR)
          end if
*
*
          if (myproc.ne.master) then
             call mpi_recv(IMZERO,1,MPI_INTEGER,MASTER,81,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,MPI_INTEGER,MASTER,82,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,MASTER,83,
     &                     JCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (myproc.ne.master) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (myproc.ne.master) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_SM(LUIN,LUOUT,SEGMNT,NOFSLV,JCOMM,LBLK)
*
* Read vector (node # mynslv) from disc file luin to array segmnt
* and send to master using communicator jcomm.
* vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
* NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SEGMNT(*)
C
*
      IREW = 0
*
      IF( IREW .NE. 0 ) THEN
        IF (MYPROC.NE.MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (MYPROC.EQ.MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      MYNSLV = NOFSLV

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (myproc.ne.master) then
              call mpi_send(LBL,1,MPI_INTEGER,MASTER,90,
     &                      JCOMM,IERR)
            end if
            if (myproc.eq.master) then
              call mpi_recv(LBL,1,MPI_INTEGER,MYNSLV,90,
     &                      JCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (myproc.ne.master) then
            READ(LUIN) LBL
            call mpi_send(LBL,1,MPI_INTEGER,MASTER,90,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.eq.master) then
            call mpi_recv(LBL,1,MPI_INTEGER,MYNSLV,90,
     &                    JCOMM,ISTAT,IERR)
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
            end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (myproc.ne.master) then
            CALL IFRMDS(LBL,1,-1,LUIN)
            call mpi_send(LBL,1,MPI_INTEGER,MASTER,90,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.eq.master) then
            call mpi_recv(LBL,1,MPI_INTEGER,MYNSLV,90,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUOUT)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (myproc.ne.master) then
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          call mpi_send(IMZERO,1,MPI_INTEGER,MASTER,91,
     &                  JCOMM,IERR)
          call mpi_send(IAMPACK,1,MPI_INTEGER,MASTER,92,
     &                  JCOMM,IERR)
          call mpi_send(SEGMNT,LBL,MPI_REAL8,MASTER,93,
     &                  JCOMM,IERR)
          end if
*
*
          if (myproc.eq.master) then
             call mpi_recv(IMZERO,1,MPI_INTEGER,MYNSLV,91,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,MPI_INTEGER,MYNSLV,92,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,MYNSLV,93,
     &                     JCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (myproc.eq.master) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (myproc.eq.master) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
*****************************************************
#else
* dummy routine for normal compilation
       SUBROUTINE DO_NOTHING(NOTHING)
       implicit none
       integer NOTHING
        write(6,*) 'let us do nothing',nothing  
       RETURN 
       END 

#endif 
