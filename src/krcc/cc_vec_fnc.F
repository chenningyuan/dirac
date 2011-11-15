      SUBROUTINE CC_VEC_FNC_KRCC(CC_AMP,CC_VEC,E_CC,ITER,
     &                          WORK,KFREE,LFREE)
*
* Calculate energy and CC vector function for 
* a set of CC amplitudes defined by CC_AMP
*
* New version by Lasse 2010
*
* Now for all possible contractions in the (p,h) scheme
* Shown as Operator --Type of contarction with T --> New operator
* H for Hamiltonian
* M for intermediate (easier to see)
* T for Cluster operator (no index needed)
* Indices are free of the ph definition
*
* This version is programmed as shown below. This is not yet the optimum
* but can "easily" be changed. These are in execution order without
* scaling considerations.
*
* First H00
*
* H00 -- 00 -- T 
*
* Second all single lines
*
* H22 -- 22 -- T
* H21 -- 21 -- T
* H12 -- 12 -- T
* H20 -- 20 -- T
*
* Third through M10
*
* H22 -- 12 -- M10
* H21 -- 11 -- M10
* H20 -- 10 -- M10 (multiply by 1/2)
*
* M10 + H10 -- 10 -- T
*
* Fourth through M02
*
* H22 -- 10 -- M12 (multiply by 1/2)
* M12 + H12 -- 10 -- M02 
* M02 + H02 -- 02 -- T (save M02 + H02)
*
* Fifth through M11
*
* H22 -- 11 -- M11 (multiply by 1/2, save M11 before multiplication)
* M11 + H11 -- 11 -- T 
*
* Sixth through M01
*
* H22 -- 11 -- M11 (reuse the above M11)
* H21 -- 10 -- M11
* M11 + H11 -- 10 -- M01
* H02 + M02 -- 01 -- M01 (reuse the above M02 + H02, multiply by 1/2)
* H22 -- 21 -- M01
* H21 -- 20 -- M01
* H12 -- 11 -- M01
* M01 + H01 -- 01 -- T
*
* I guess that is all we need then!
*
*      
      IMPLICIT REAL*8(A-H,O-Z)
*.Input   
      DIMENSION CC_AMP(*)
*. Output   
      DIMENSION CC_VEC(*)
*
      DIMENSION WORK(*)
*.
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "clunit.inc"
#include "glbbas.inc"
#include "crun.inc"
#include "cands.inc"
#include "cstate.inc"
#include "ctcc.inc"
#include "ctccp.inc"
#include "cgas.inc"
#include "cecore.inc"
#include "interm.inc"
      COMMON/CMXCJ/MXCJ
*. Local scratch 
      INTEGER ISX(2*MXPNGAS*MXPNGAS)
*
      REAL*8 INPRDD, INPROD
*
      LBLK = -1
      NTEST = 10
*
* Will determine the number of commutators (Of course one one with zero)
      IONE = 1
      IZERO = 0
      ZERO = 0.0D0
      HALF = 0.5D0
      IDOUBLE = 2.0D0
*         
* Real or complex algebra
      IF(IRECOM.EQ.1) THEN
        IMULTFAC = 1
      ELSE
        IMULTFAC = 2
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' CC_VEC_FNC in action '
        WRITE(6,*) ' ====================='
      END IF
*
* Hardwire symmetry
      IHSM = 1
      ITSM = 1
*      
*      May want to allocate some arrays
*
      CALL EXP_MT_H_EXP_TC_MEM_KRCC(WORK,KFREE,LFREE)
*
*. Start by zeroing the CC-vector function and the intermediates
      ZERO = 0.0D0
      CALL SETVEC(CC_VEC,ZERO,N_CC_AMP+1)
      CALL SETVEC(WORK(KM01),ZERO,LEN_M01_VEC_CC)
      CALL SETVEC(WORK(KM10),ZERO,LEN_M10_VEC_CC)
      CALL SETVEC(WORK(KM11),ZERO,LEN_M11_VEC_CC)
      CALL SETVEC(WORK(KM02),ZERO,LEN_M02_VEC_CC)
      CALL SETVEC(WORK(KM12),ZERO,LEN_M12_VEC_CC)
*
* Start with zero commutators! (If first iteration and no restart this
* should be the only one done!)
*
* H00 --00--> T SHOULD BE THE ONLY ONE DONE FOR IT=1 AND NO RESTART
* 
* First H00
*
* H00 -- 00 -- T 
      IHM = 1 !Still the Hamiltonian that is contracted from 
      IFHM = 1 !integrals are fetched
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IZERO,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TP_CC,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH00,WORK(KNH00),
     &     NH00,WORK(KM00TOH00),WORK(KM00TOH00),
     &     IDUM,IDUM,
     &     IDUM,IDUM,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
! First iteration with no restart only do H00
      IF(IFORM.EQ.1.AND.ITER.EQ.1) RETURN
*
* Second all single lines
*
* H22 -- 22 -- T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOT,WORK(KNH22TOT),WORK(KNH22TT),WORK(KCOMBNH22TOT),
     &     WORK(KPERMNH22TOT),WORK(KFACNH22TOT),WORK(KCOMBFACNH22TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H21 -- 21 -- T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOT,WORK(KNH21TOT),WORK(KNH21TT),WORK(KCOMBNH21TOT),
     &     WORK(KPERMNH21TOT),WORK(KFACNH21TOT),WORK(KCOMBFACNH21TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H12 -- 12 -- T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH12,WORK(KNH12),
     &     NH12TOT,WORK(KNH12TOT),WORK(KNH12TT),WORK(KCOMBNH12TOT),
     &     WORK(KPERMNH12TOT),WORK(KFACNH12TOT),WORK(KCOMBFACNH12TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H20 -- 20 -- T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH20,WORK(KNH20),
     &     NH20TOT,WORK(KNH20TOT),WORK(KNH20TT),WORK(KCOMBNH20TOT),
     &     WORK(KPERMNH20TOT),WORK(KFACNH20TOT),WORK(KCOMBFACNH20TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Third through M10
*
* H22 -- 12 -- M10
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM10,WORK(KNH22TOM10),WORK(KNH22T10),
     &     WORK(KCOMBNH22TOM10),WORK(KPERMNH22TOM10),
     &     WORK(KFACNH22TOM10),WORK(KCOMBFACNH22TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H21 -- 11 -- M10
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOM10,WORK(KNH21TOM10),WORK(KNH21T10),
     &     WORK(KCOMBNH21TOM10),WORK(KPERMNH21TOM10),
     &     WORK(KFACNH21TOM10),WORK(KCOMBFACNH21TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H20 -- 10 -- M10 (multiply by 1/2)
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH20,WORK(KNH20),
     &     NH20TOM10,WORK(KNH20TOM10),WORK(KNH20T10),
     &     WORK(KCOMBNH20TOM10),WORK(KPERMNH20TOM10),
     &     WORK(KFACNH20TOM10),WORK(KCOMBFACNH20TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Multiply by 1/2
*
      CALL ISCAVE(WORK(KM10),HALF,LEN_M10_VEC_CC)
*
* M10 + H10 -- 10 -- T
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH10,WORK(KNH10),
     &                                   NINTER10,WORK(KINTM10),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM10),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM10), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NINTER10,WORK(KINTM10),WORK(KLIHIND),IHSM,
     &     NINTER10,WORK(KLLSOBEX_M10),
     &     NM10TOT,WORK(KNM10TOT),WORK(KNM10TT),
     &     WORK(KCOMBNM10TOT),WORK(KPERMNM10TOT),
     &     WORK(KFACNM10TOT),WORK(KCOMBFACNM10TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Fourth through M02
*
* H22 -- 10 -- M12 (multiply by 1/2)
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM12),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M12),WORK(KLLSOBEX_M12),NINTER12,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM12,WORK(KNH22TOM12),WORK(KNH22T12),
     &     WORK(KCOMBNH22TOM12),WORK(KPERMNH22TOM12),
     &     WORK(KFACNH22TOM12),WORK(KCOMBFACNH22TOM12),
     &     WORK(KINTM12),NINTER12,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Multiply by 1/2
*
      CALL ISCAVE(WORK(KM12),HALF,LEN_M12_VEC_CC)
*
* M12 + H12 -- 10 -- M02 
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH12,WORK(KNH12),
     &                                   NINTER12,WORK(KINTM12),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM12),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM02),IONE,IHM,
     &     IFHM,WORK(KM12), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NINTER12,WORK(KINTM12),WORK(KLIHIND),IHSM,
     &     NINTER12,WORK(KLLSOBEX_M12),
     &     NM12TOM02,WORK(KNM12TOM02),WORK(KNM12T02),
     &     WORK(KCOMBNM12TOM02),WORK(KPERMNM12TOM02),
     &     WORK(KFACNM12TOM02),WORK(KCOMBFACNM12TOM02),
     &     WORK(KINTM02),NINTER02,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* M02 + H02 -- 02 -- T (save M02 + H02)
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH02,WORK(KNH02),
     &                                   NINTER02,WORK(KINTM02),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM02),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* Save vector to disc since it will be needed later
*
      CALL VEC_TO_DISC(WORK(KM02),LEN_M02_VEC_CC,1,LBLK,LUSC5)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM02), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NINTER02,WORK(KINTM02),WORK(KLIHIND),IHSM,
     &     NINTER02,WORK(KLLSOBEX_M02),
     &     NM02TOT,WORK(KNM02TOT),WORK(KNM02TT),
     &     WORK(KCOMBNM02TOT),WORK(KPERMNM02TOT),
     &     WORK(KFACNM02TOT),WORK(KCOMBFACNM02TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Fifth through M11
*
* H22 -- 11 -- M11 (multiply by 1/2)
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM11),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M11),WORK(KLLSOBEX_M11),NINTER11,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM11,WORK(KNH22TOM11),WORK(KNH22T11),
     &     WORK(KCOMBNH22TOM11),WORK(KPERMNH22TOM11),
     &     WORK(KFACNH22TOM11),WORK(KCOMBFACNH22TOM11),
     &     WORK(KINTM11),NINTER11,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Save vector to disc since it will be needed later
*
      CALL VEC_TO_DISC(WORK(KM11),LEN_M11_VEC_CC,1,LBLK,LUSC4)
*
* Multiply by 1/2
*
      CALL ISCAVE(WORK(KM11),HALF,LEN_M11_VEC_CC)
*
* M11 + H11 -- 11 -- T 
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH11,WORK(KNH11),
     &                                   NINTER11,WORK(KINTM11),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM11),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM11), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NINTER11,WORK(KINTM11),WORK(KLIHIND),IHSM,
     &     NINTER11,WORK(KLLSOBEX_M11),
     &     NM11TOT,WORK(KNM11TOT),WORK(KNM11TT),
     &     WORK(KCOMBNM11TOT),WORK(KPERMNM11TOT),
     &     WORK(KFACNM11TOT),WORK(KCOMBFACNM11TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* Sixth through M01
*
* H22 -- 11 -- M11 (reuse the above M11)
*
      CALL VEC_FROM_DISC(WORK(KM11),LEN_M11_VEC_CC,1,LBLK,LUSC4)
*
* H21 -- 10 -- M11
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM11),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M11),WORK(KLLSOBEX_M11),NINTER11,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOM11,WORK(KNH21TOM11),WORK(KNH21T11),
     &     WORK(KCOMBNH21TOM11),WORK(KPERMNH21TOM11),
     &     WORK(KFACNH21TOM11),WORK(KCOMBFACNH21TOM11),
     &     WORK(KINTM11),NINTER11,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* M11 + H11 -- 10 -- M01
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH11,WORK(KNH11),
     &                                   NINTER11,WORK(KINTM11),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM11),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KM12), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NINTER12,WORK(KINTM12),WORK(KLIHIND),IHSM,
     &     NINTER12,WORK(KLLSOBEX_M12),
     &     NM12TOM02,WORK(KNM12TOM02),WORK(KNM12T02),
     &     WORK(KCOMBNM12TOM02),WORK(KPERMNM12TOM02),
     &     WORK(KFACNM12TOM02),WORK(KCOMBFACNM12TOM02),
     &     WORK(KINTM02),NINTER02,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H02 + M02 -- 01 -- M01 (reuse the above M02 + H02, multiply by 1/2)
*
      CALL VEC_FROM_DISC(WORK(KM02),LEN_M02_VEC_CC,1,LBLK,LUSC5)
*
* Multiply by 1/2
*
      CALL ISCAVE(WORK(KM02),HALF,LEN_M02_VEC_CC)
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH02,WORK(KNH02),
     &     NH02TOM01,WORK(KNH02TOM01),WORK(KNH02T01),
     &     WORK(KCOMBNH02TOM01),WORK(KPERMNH02TOM01),
     &     WORK(KFACNH02TOM01),WORK(KCOMBFACNH02TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H22 -- 21 -- M01
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH22,WORK(KNH22),
     &     NH22TOM01,WORK(KNH22TOM01),WORK(KNH22T01),
     &     WORK(KCOMBNH22TOM01),WORK(KPERMNH22TOM01),
     &     WORK(KFACNH22TOM01),WORK(KCOMBFACNH22TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H21 -- 20 -- M01
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH21,WORK(KNH21),
     &     NH21TOM01,WORK(KNH21TOM01),WORK(KNH21T01),
     &     WORK(KCOMBNH21TOM01),WORK(KPERMNH21TOM01),
     &     WORK(KFACNH21TOM01),WORK(KCOMBFACNH21TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* H12 -- 11 -- M01
*
      IHM = 1 !Still the Hamiltonian that is contracted from
      IFHM = 1 !integrals are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM01),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_M01),WORK(KLLSOBEX_M01),NINTER01,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NH12,WORK(KNH12),
     &     NH12TOM01,WORK(KNH12TOM01),WORK(KNH12T01),
     &     WORK(KCOMBNH12TOM01),WORK(KPERMNH12TOM01),
     &     WORK(KFACNH12TOM01),WORK(KCOMBFACNH12TOM01),
     &     WORK(KINTM01),NINTER01,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
*
* M01 + H01 -- 01 -- T
*
* Add amplitudes and integrals 
! Careful with signs since IHINDX may be different for the two!!!
*
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH01,WORK(KNH01),
     &                                   NINTER01,WORK(KINTM01),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM01),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
*
* And now the contraction
*
      IHM = 2
      IFHM = 2 !intermediates are fetched
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM01), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
     &     NINTER01,WORK(KINTM01),WORK(KLIHIND),IHSM,
     &     NINTER01,WORK(KLLSOBEX_M01),
     &     NM01TOT,WORK(KNM01TOT),WORK(KNM01TT),
     &     WORK(KCOMBNM01TOT),WORK(KPERMNM01TOT),
     &     WORK(KFACNM01TOT),WORK(KCOMBFACNM01TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
* End of the line

*
      E_CC = CC_VEC(N_CC_AMP+1)
      print*,CC_VEC(N_CC_AMP+1),CC_VEC(1)
*
*. Mission over, print 
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' CC-vector function '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' DO NOT USE WRT_CC_VEC2 YET '
        STOP ' DO NOT USE WRT_CC_VEC2 YET '
        CALL WRT_CC_VEC2(CC_VEC,6,CCTYPE,
     &                   WORK,KFREE,LFREE)
        WRITE(6,*) ' Initial amplitudes '
        CALL WRT_CC_VEC2(CC_AMP,6,CCTYPE,
     &                   WORK,KFREE,LFREE)
      END IF
* 
      STOP ' Enforced stop after CC_VEC_FNC'
      RETURN
      END
*
****************************************************************
*                                                              *
* Silly subroutine by Lasse to find largest difference between *
* neighbouring elements in an array.                           *
*                                                              *
****************************************************************
      SUBROUTINE FIND_LARGEST_DIFF_FOR_NEIGHBOURS(NOUT,NARR,NLENGHT)
      IMPLICIT NONE
      INTEGER NOUT,NLENGHT,NARR(NLENGHT),I,J
      NOUT = 0
      DO I=2,NLENGHT
        J = NARR(I)-NARR(I-1)
        IF(J.GT.NOUT) NOUT = J
      END DO
      END
