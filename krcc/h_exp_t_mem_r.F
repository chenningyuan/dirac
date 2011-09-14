      SUBROUTINE EXP_MT_H_EXP_TC_MEM_KRCC(WORK,KFREE,LFREE)
*
* set up scratch space for direct calculation of H EXP T
*
* Pointers to scratch are returned in CC_SCR
*
*. Jeppe Olsen, May 2000
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cgas.inc"
#include "cc_scr2.inc"
*
*
      DIMENSION WORK(*)
*
* Setting a default batch length. Increase or decrease according to
* memory available and speed
*
      LCCB = 35
      WRITE(6,*) ' Batch lenght LCCB ',LCCB
*
*.1 : KCCF : Intermediate F-vectors for needed out of space 
*         excitations in the form sum_I T(I) F(I) 
*
*. Largest length of F-vector and T- block for given type of T-block
*. ( It is convenient to store these together )
C            LEN_FOR_TF(NSPOBEX,ISPOBEX,LSPOBEX,NGAS,LEN_TF)
C       CALL LEN_FOR_TF_REL(NSPOBEX_TP_CC,WORK(KLSOBEX_CC),
C    &                  WORK(KLLSOBEX),NGAS,LEN_TF,
C    &                  WORK,KFREE,LFREE)
C       WRITE(6,*) ' ... MEM, LEN_TF (1) = ', LEN_TF
*. It is convenient to be able to use F as a normal CC also so
C       N_CC_AMPP1 = N_CC_AMP + 1
C       LEN_TF = MAX(LEN_TF,N_CC_AMPP1)
C       CALL MEMMAR(KCCF,LEN_TF,'ADDL  ',2,'CCF   ')
C       CALL MEMGET('REAL',KCCF,LEN_TF,WORK,KFREE,LFREE)
*
* For a given orbital type I the start of the corresponding 
* intermediate array is given as KCCF    + (I-1)*N_CC_AMP
*
*.2 : KCCVEC1,KCCVEC2,KCCE    : Three  CC vectors, may allowed to hold 
*.    component corresponding to unit operator 
C     CALL MEMMAR(KCCE   ,N_CC_AMPP1,'ADDL  ',2,'CCE   ')
C     CALL MEMGET('REAL',KCCE,N_CC_AMPP1,WORK,KFREE,LFREE)
C     CALL MEMMAN(KCCVEC1,N_CC_AMPP1,'ADDL  ',2,'CCVEC1')
C     KCCVEC1 = 0
C     KCCVEC2 = 0
*.3 : Maps for T1T2 => T12 for individual strings types 
      CALL DIM_T1T2_TO_T12_MAP_REL(LEN_T1T2_STRING,LENT_T1T2_TCC,
     &                             WORK,KFREE,LFREE)
      LEN = LEN_T1T2_STRING
*. Arrays for holding elementary operators times strings
*
      CALL MEMGET('INTE',KIX1_CA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_CA,LEN,WORK,KFREE,LFREE)
*
      CALL MEMGET('INTE',KIX1_CB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_CB,LEN,WORK,KFREE,LFREE)
*
      CALL MEMGET('INTE',KIX1_AA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_AA,LEN,WORK,KFREE,LFREE)
*
      CALL MEMGET('INTE',KIX1_AB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_AB,LEN,WORK,KFREE,LFREE)
*
      CALL MEMGET('INTE',KIX2_CA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_CA,LEN,WORK,KFREE,LFREE)
C
      CALL MEMGET('INTE',KIX2_CB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_CB,LEN,WORK,KFREE,LFREE)
C
      CALL MEMGET('INTE',KIX2_AA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_AA,LEN,WORK,KFREE,LFREE)
C
      CALL MEMGET('INTE',KIX2_AB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_AB,LEN,WORK,KFREE,LFREE)
C
*
*.4 : KLZ,KLZSCR : Memory for a Z matrix and scratch for 
*     constructing Z
*     
      IATP = 1
      IBTP = 2
      NAEL = NELFTP(IATP)
      NBEL = NELFTP(IBTP)
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL MEMGET('INTE',KLZ,LZ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLZSCR,LZSCR,WORK,KFREE,LFREE)
*. String occupations for given CAAB, all symmetris
      LEN = MX_ST_TSOSO_BLK_MX*NSMST*4
      print*,'MX_ST_TSOSO_BLK_MX',MX_ST_TSOSO_BLK_MX
      CALL MEMGET('INTE',KLSTOCC1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTOCC2,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTOCC3,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTOCC4,LEN,WORK,KFREE,LFREE)
      print*,'KLSTOCC1,KLSTOCC2',KLSTOCC1,KLSTOCC2
*. Reorder array for given CAAB, all symmetries 
      CALL MEMGET('INTE',KLSTREO,MX_ST_TSOSO_MX*NSMST,WORK,KFREE,LFREE)
*. Intermediate blocks, (LCCB,LCCB)
*. For the moment 
      WRITE(6,*) ' SET DIMENSION FROM LCCB '
      LB = LCCB
      LEN = LCCB**2
      LEN3 = LCCB**3
      CALL MEMGET('REAL',KLTSCR1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLTSCR2,LEN,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KLTSCR3,LEN,WORK,KFREE,LFREE)
      print*,'LEN,LCCB',LEN,LCCB
*
*and for a batch of Hamiltonian
      LOPSCR = LCCB*LCCB*LB
*. The above is not correct as the excitation part of the 
*. integrals are obtained as a single block. 
*. The largest block of integrals (hp!hp) should therefore
*. be contained in OPSCR. I will just pt just use an 
*. a block that is too large ...
C will now instead calculate it from largest block of integrals.
C Lasse
      CALL FIND_LARGEST_DIFF_FOR_NEIGHBOURS(LOPSCR,IOFFINTTYPE,
     &                                      NSPOBEX_TP)
      CALL MEMGET('REAL',KLOPSCR,LOPSCR,WORK,KFREE,LFREE)
*. Number and symmetries of each substring for 6 complete T-blocks 
*. Largest number of strings in intermediate arrays
*
*. D1, D2, EX holds part of Hamiltonian so
*. MXTSOB_P is largest number of particle orbitals of given type, ALL SYM!
C     LEN = 4*MXTSOB_P**2 * MXTSOB_H**2
      LEN = 4*LB
      CALL MEMGET('INTE',KLSMD1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNMD1,LEN,WORK,KFREE,LFREE)
*
*
      CALL MEMGET('INTE',KLSMEX,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNMEX,LEN,WORK,KFREE,LFREE)
*. K1, K2, L1 holds general strings
      LEN = 4*LB
      CALL MEMGET('INTE',KLSMK1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNMK1,LEN,WORK,KFREE,LFREE)
*
*
      CALL MEMGET('INTE',KLSML1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNML1,LEN,WORK,KFREE,LFREE)
*     KLOCK1, KLOCK2, KLOCL1
      LEN = 4*NGAS
      CALL MEMGET('INTE',KLOCK1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLOCL1,LEN,WORK,KFREE,LFREE)
*
      RETURN
      END
