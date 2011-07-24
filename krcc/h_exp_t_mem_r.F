      SUBROUTINE EXP_MT_H_EXP_TC_MEM_KRCC(WORK,KFREE,LFREE)
*
* set up scratch space for direct calculation of H EXP T
*
* Pointers to scratch are returned in CC_SCR
*
*. Jeppe Olsen, May 2000
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cgas.inc"
#include "cc_scr2.inc"
!#include "wrkspc.inc"
*
!     COMMON/CC_SCR2/ KCCF,KCCE,KCCVEC1,KCCVEC2,
!    &              KISTST1,KXSTST1,!KISTST2,KXSTST2,
C    &              KISTST3,KXSTST3,KISTST4,KXSTST4,
!    &              KLZ,KLZSCR,KLSTOCC1,!KLSTOCC2,
C    &              KLSTOCC3,KLSTOCC4,
!    &              KLSTOCC5,KLSTOCC6,KLSTREO,
!    &              KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
!    &              KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
C    &              KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
C    &              KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
C    &              KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
C    &              KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
C    &              KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
C    &              KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
!    &              KLTSCR1,!KLTSCR2,KLTSCR3,KLTSCR4,
!    &              KLOPSCR,
!    &              KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
!    &              KLSMD1,KLSMD2,KLSMEX,KLSMK1,KLSMK2,KLSML1,
!    &              KLNMD1,KLNMD2,KLNMEX,KLNMK1,KLNMK2,KLNML1,
!    &              KLOCK1, KLOCK2, KLOCL1, KLOCOT1T2, KL_IBF, 
!    &              KLEXEORD
C     COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,
C    &               KLSTOCC5,KLSTOCC6,KLSTOCC7,
C    &               KLSTOCC8,KLSTREO,
C    &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
C    &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
C    &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
C    &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,

*
      DIMENSION WORK(*)
*
* Setting a default batch length. Increase or decrease according to
* memory available and speed
*
      LCCB = 35
      WRITE(6,*) ' Batch lenght LCCB ',LCCB
      LCCBD12 = LCCB
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
C     CALL MEMMAR(KISTST1,LEN,'ADDL  ',1,'ISTST1')
C     CALL MEMMAR(KXSTST1,LEN,'ADDL  ',2,'XSTST1')
*
C     CALL MEMMAR(KISTST2,LEN,'ADDL  ',1,'ISTST2')
C     CALL MEMMAR(KXSTST2,LEN,'ADDL  ',2,'XSTST2')
C
C     CALL MEMMAR(KISTST3,LEN,'ADDL  ',1,'ISTST3')
C     CALL MEMMAR(KXSTST3,LEN,'ADDL  ',2,'XSTST3')
C
C     CALL MEMMAR(KISTST4,LEN,'ADDL  ',1,'ISTST4')
C     CALL MEMMAR(KXSTST4,LEN,'ADDL  ',2,'XSTST4')
*
C     CALL MEMMAR(KIX1_CA,LEN,'ADDL  ',1,'IX1_CA')
C     CALL MEMMAR(KSX1_CA,LEN,'ADDL  ',2,'SX1_CA')
      CALL MEMGET('INTE',KIX1_CA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_CA,LEN,WORK,KFREE,LFREE)
*
C     CALL MEMMAR(KIX1_CB,LEN,'ADDL  ',1,'IX1_CB')
C     CALL MEMMAR(KSX1_CB,LEN,'ADDL  ',2,'SX1_CB')
      CALL MEMGET('INTE',KIX1_CB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_CB,LEN,WORK,KFREE,LFREE)
*
C     CALL MEMMAR(KIX1_AA,LEN,'ADDL  ',1,'IX1_AA')
C     CALL MEMMAR(KSX1_AA,LEN,'ADDL  ',2,'SX1_AA')
      CALL MEMGET('INTE',KIX1_AA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_AA,LEN,WORK,KFREE,LFREE)
*
C     CALL MEMMAR(KIX1_AB,LEN,'ADDL  ',1,'IX1_AB')
C     CALL MEMMAR(KSX1_AB,LEN,'ADDL  ',2,'SX1_AB')
      CALL MEMGET('INTE',KIX1_AB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX1_AB,LEN,WORK,KFREE,LFREE)
*
C     CALL MEMMAR(KIX2_CA,LEN,'ADDL  ',1,'IX2_CA')
C     CALL MEMMAR(KSX2_CA,LEN,'ADDL  ',2,'SX2_CA')
      CALL MEMGET('INTE',KIX2_CA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_CA,LEN,WORK,KFREE,LFREE)
C
C     CALL MEMMAR(KIX2_CB,LEN,'ADDL  ',1,'IX2_CB')
C     CALL MEMMAR(KSX2_CB,LEN,'ADDL  ',2,'SX2_CB')
      CALL MEMGET('INTE',KIX2_CB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_CB,LEN,WORK,KFREE,LFREE)
C
C     CALL MEMMAR(KIX2_AA,LEN,'ADDL  ',1,'IX2_AA')
C     CALL MEMMAR(KSX2_AA,LEN,'ADDL  ',2,'SX2_AA')
      CALL MEMGET('INTE',KIX2_AA,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_AA,LEN,WORK,KFREE,LFREE)
C
C     CALL MEMMAR(KIX2_AB,LEN,'ADDL  ',1,'IX2_AB')
C     CALL MEMMAR(KSX2_AB,LEN,'ADDL  ',2,'SX2_AB')
      CALL MEMGET('INTE',KIX2_AB,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSX2_AB,LEN,WORK,KFREE,LFREE)
C
C     CALL MEMMAR(KIX3_CA,LEN,'ADDL  ',1,'IX3_CA')
C     CALL MEMMAR(KSX3_CA,LEN,'ADDL  ',2,'SX3_CA')
C
C     CALL MEMMAR(KIX3_CB,LEN,'ADDL  ',1,'IX3_CB')
C     CALL MEMMAR(KSX3_CB,LEN,'ADDL  ',2,'SX3_CB')
C
C     CALL MEMMAR(KIX3_AA,LEN,'ADDL  ',1,'IX3_AA')
C     CALL MEMMAR(KSX3_AA,LEN,'ADDL  ',2,'SX3_AA')
C
C     CALL MEMMAR(KIX3_AB,LEN,'ADDL  ',1,'IX3_AB')
C     CALL MEMMAR(KSX3_AB,LEN,'ADDL  ',2,'SX3_AB')
C
C     CALL MEMMAR(KIX4_CA,LEN,'ADDL  ',1,'IX4_CA')
C     CALL MEMMAR(KSX4_CA,LEN,'ADDL  ',2,'SX4_CA')
C
C     CALL MEMMAR(KIX4_CB,LEN,'ADDL  ',1,'IX4_CB')
C     CALL MEMMAR(KSX4_CB,LEN,'ADDL  ',2,'SX4_CB')
C
C     CALL MEMMAR(KIX4_AA,LEN,'ADDL  ',1,'IX4_AA')
C     CALL MEMMAR(KSX4_AA,LEN,'ADDL  ',2,'SX4_AA')
C
C     CALL MEMMAR(KIX4_AB,LEN,'ADDL  ',1,'IX4_AB')
C     CALL MEMMAR(KSX4_AB,LEN,'ADDL  ',2,'SX4_AB')
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
C     CALL MEMMAR(KLZ,LZ,'ADDL  ',1,'Z     ')
C     CALL MEMMAR(KLZSCR,LZSCR,'ADDL  ',1,'ZSCR  ')
      CALL MEMGET('INTE',KLZ,LZ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLZSCR,LZSCR,WORK,KFREE,LFREE)
*. String occupations for given CAAB, all symmetris
      LEN = MX_ST_TSOSO_BLK_MX*NSMST*4
C     CALL MEMMAR(KLSTOCC1,LEN,'ADDL  ',1,'STOCC1')
C     CALL MEMMAR(KLSTOCC2,LEN,'ADDL  ',1,'STOCC2')
C     CALL MEMMAR(KLSTOCC3,LEN,'ADDL  ',1,'STOCC3')
C     CALL MEMMAR(KLSTOCC4,LEN,'ADDL  ',1,'STOCC4')
C     CALL MEMMAR(KLSTOCC5,LEN,'ADDL  ',1,'STOCC5')
C     CALL MEMMAR(KLSTOCC6,LEN,'ADDL  ',1,'STOCC6')
C     CALL MEMMAR(KLSTOCC7,LEN,'ADDL  ',1,'STOCC7')
C     CALL MEMMAR(KLSTOCC8,LEN,'ADDL  ',1,'STOCC8')
      print*,'MX_ST_TSOSO_BLK_MX',MX_ST_TSOSO_BLK_MX
      CALL MEMGET('INTE',KLSTOCC1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTOCC5,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTOCC6,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLSTOCC7,LEN,WORK,KFREE,LFREE)
      print*,'KLSTOCC1,KLSTOCC5',KLSTOCC1,KLSTOCC5
*. Reorder array for given CAAB, all symmetries 
C     CALL MEMMAR(KLSTREO,MX_ST_TSOSO_MX*NSMST,'ADDL  ',1,'STREO ')
      CALL MEMGET('INTE',KLSTREO,MX_ST_TSOSO_MX*NSMST,WORK,KFREE,LFREE)
*. Intermediate blocks, (LD12B,LD12B)
*. For the moment 
      WRITE(6,*) ' SET DIMENSION OF LCCBD12 '
      LD12B = LCCBD12
      LB = LCCB
      LEN = LD12B**2
      LEN3 = LD12B**3
C     CALL MEMMAR(KLTSCR1,LEN,'ADDL  ',2,'TSCR1 ')
C     CALL MEMMAR(KLTSCR2,LEN3,'ADDL  ',2,'TSCR2 ')
C     CALL MEMMAR(KLTSCR3,LEN,'ADDL  ',2,'TSCR3 ')
      CALL MEMGET('REAL',KLTSCR1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLTSCR2,LEN,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KLTSCR3,LEN,WORK,KFREE,LFREE)
      print*,'LEN,LCCBD12,LCCB',LEN,LCCBD12,LCCB
*
C          K_RES_DIM(ISPOBEX,NSPOBEX,MAXOP,NELMNT_MAX)
      CALL K_RES_DIM_REL(WORK(KLSOBEX_CC),NSPOBEX_TP_CC,4,MX_KBLK)    
COLD  CALL T_DL_DIM(WORK(KLSOBEX),NSPOBEX_TP_CC,4,N_TDL_MAX,LD12B)
C     WRITE(6,*) ' MX_KBLK in H_EXP ... = ', MX_KBLK
*. A big one for storing T(D2,L)
      LEN = N_TDL_MAX
C     CALL MEMMAR(KLTSCR4,LEN,'ADDL  ',2,'TSCR4 ')
*. For batches of indeces of Hamiltionian operator
      LEN = 4*LD12B
C     CALL MEMMAR(KLIOD1_ST,LEN,'ADDL  ',2,'IOD1_S')
C     CALL MEMMAR(KLIOD2_ST,LEN,'ADDL  ',2,'IOD2_S')
C     CALL MEMMAR(KLIOEX_ST,LEN,'ADDL  ',2,'IOEX_S')
      CALL MEMGET('REAL',KLIOD1_ST,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLIOD2_ST,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLIOEX_ST,LEN,WORK,KFREE,LFREE)
*and for a batch of Hamiltonian
      LOPSCR = LD12B*LD12B*LB
*. The above is not correct as the excitation part of the 
*. integrals are obtained as a single block. 
*. The largest block of integrals (hp!hp) should therefore
*. be contained in OPSCR. I will just pt just use an 
*. a block that is too large ...
C will now instead calculate it from largest block of integrals.
C Lasse
      CALL FIND_LARGEST_DIFF_FOR_NEIGHBOURS(LOPSCR,IOFFINTTYPE,
     &                                      NSPOBEX_TP)
C     CALL MEMMAR(KLOPSCR,LOPSCR,'ADDL  ',2,'OPSCR ')
      print*,'LOPSCR',LOPSCR
      CALL MEMGET('REAL',KLOPSCR,LOPSCR,WORK,KFREE,LFREE)
*. Number and symmetries of each substring for 6 complete T-blocks 
*. Largest number of strings in intermediate arrays
*
*. D1, D2, EX holds part of Hamiltonian so
*. MXTSOB_P is largest number of particle orbitals of given type, ALL SYM!
C     LEN = 4*MXTSOB_P**2 * MXTSOB_H**2
      LEN = 4*LB
C     CALL MEMMAR(KLSMD1,LEN,'ADDL  ',1,'SM_D1 ')
C     CALL MEMMAR(KLNMD1,LEN,'ADDL  ',1,'NM_D1 ')
      CALL MEMGET('INTE',KLSMD1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNMD1,LEN,WORK,KFREE,LFREE)
*
*
C     CALL MEMMAR(KLSMD2,LEN,'ADDL  ',1,'SM_D2 ')
C     CALL MEMMAR(KLNMD2,LEN,'ADDL  ',1,'NM_D2 ') 
*
C     CALL MEMMAR(KLSMEX,LEN,'ADDL  ',1,'SM_EX ')
C     CALL MEMMAR(KLNMEX,LEN,'ADDL  ',1,'NM_EX ')
      CALL MEMGET('INTE',KLSMEX,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNMEX,LEN,WORK,KFREE,LFREE)
*. K1, K2, L1 holds general strings
      LEN = 4*LB
C     CALL MEMMAR(KLSMK1,LEN,'ADDL  ',1,'SM_K1 ')
C     CALL MEMMAR(KLNMK1,LEN,'ADDL  ',1,'NM_K1 ')
      CALL MEMGET('INTE',KLSMK1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNMK1,LEN,WORK,KFREE,LFREE)
*
*
C     CALL MEMMAR(KLSMK2,LEN,'ADDL  ',1,'SM_K2 ')
C     CALL MEMMAR(KLNMK2,LEN,'ADDL  ',1,'NM_K2 ') 
*
C     CALL MEMMAR(KLSML1,LEN,'ADDL  ',1,'SM_L1 ')
C     CALL MEMMAR(KLNML1,LEN,'ADDL  ',1,'NM_L1 ')
      CALL MEMGET('INTE',KLSML1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLNML1,LEN,WORK,KFREE,LFREE)
*     KLOCK1, KLOCK2, KLOCL1
      LEN = 4*NGAS
C     CALL MEMMAR(KLOCK1,LEN,'ADDL  ',1,'IOC_K1')
C     CALL MEMMAR(KLOCK2,LEN,'ADDL  ',1,'IOC_K2')
C     CALL MEMMAR(KLOCL1,LEN,'ADDL  ',1,'IOC_L1')
C     CALL MEMMAR(KLOCOT1T2,LEN,'ADDL  ',1,'IOC_OT')
      CALL MEMGET('INTE',KLOCK1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLOCL1,LEN,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLOCOT1T2,LEN,WORK,KFREE,LFREE)
*. Array for offsets in F
C     CALL MEMMAR(KL_IBF,NSPOBEX_TP_CC,'ADDL  ',1,'IB_F  ')
      CALL MEMGET('INTE',KL_IBF,NSPOBEX_TP_CC,WORK,KFREE,LFREE)
*. Array for execution order of spinorbital types in Exp T
C     CALL MEMMAR(KLEXEORD,NSPOBEX_TP_CC,'ADDL  ',1,'EXEORD')
      CALL MEMGET('INTE',KLEXEORD,NSPOBEX_TP_CC,WORK,KFREE,LFREE)
*. Array for unit operator
      ILEN = 4*NGAS*(NSPOBEX_TP_CC+1)
      CALL MEMGET('INTE',KUNITCC,ILEN,WORK,KFREE,LFREE)
*. And set the array
      IZERO = 0
      CALL ISETVC(WORK(KUNITCC),IZERO,ILEN)
      RETURN
      END
