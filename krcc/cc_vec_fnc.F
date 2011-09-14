      SUBROUTINE CC_VEC_FNC_KRCC(CC_AMP,CC_VEC,E_CC,ITER,
     &                          WORK,KFREE,LFREE)
*
* Calculate energy and CC vector function for 
* a set of CC amplitudes defined by CC_AMP
*
* New version by Lasse 2010
*
* A note on VEC1 VEC2 (VEC3/C2) (June 2001)
* These arrays are only needed for the older versions 
* of the code, and may be dummy indeces for the newer versions
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
*. Start by zeroing the CC-vector function 
      ZERO = 0.0D0
      CALL SETVEC(CC_VEC,ZERO,N_CC_AMP+1)
* I may still need this! But let's see...
      CALL SUM_ARRAY
* And this although it is commented out at the moment
C      call set_hop_dbg_cc(4,2)
*
* Now for all possible contractions in the (p,h) scheme
* Shown as Operator --Type of contarction with T --> New operator
* H for Hamiltonian
* M for intermediate (easier to see)
* T for Cluster operator (no index needed)
* Indices are free of the ph definitionA
* Start with zero commutators! (If first iteration and no restart this
* should be the only one done!)
*
* H00 --00--> T SHOULD BE THE ONLY ONE DONE FOR IT=1 AND NO RESTART
* 
      IHM = 1 !Still the Hamiltonian that is contracted from 
      II = 0
      IFHM = 1 !integrals are fetched
      IF(II.GE.1) THEN
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IZERO,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
! Examine next 3 lines!!!      
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TP_CC,
! End examination
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
! Add info for contractions
     &     NH00,WORK(KNH00),
     &     NH00,WORK(KM00TOH00),WORK(KM00TOH00),
     &     IDUM,IDUM,
     &     IDUM,IDUM,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      stop 'just check it'
! First iteration with no restart only do H00
      IF(IFORM.EQ.1.AND.ITER.EQ.1) RETURN
      END IF
C     stop 'just check it'
*
* Adding routine for amplitudes and integrals testing
*
      II =0
      IF(II.GE.1) THEN
      CALL I_ADD_INTEGRALS_TO_AMPLITUDES(NH01,WORK(KNH01),
     &                                   NINTER01,WORK(KINTM01),
     &                                   NSPOBEX_TP,WORK(KLSOBEX),
     &                                   WORK(KT_CC),WORK(KFI),
     &                                   WORK(KM01),
     &                                   MX_ST_TSOSO_BLK_MX,
     &                                   WORK,KFREE,LFREE)
      stop ' added integrals and amplitudes'
      END IF
*
* H22 --11--> M11'  (Special case of M11)
*
*
* M11' --11--> T    (O11' now goes to M11)
*
      IHM = 2
      IFHM = 2 !integrals are fetched
      II = 0
      IF(II.GE.1) THEN
*
      stop ' not yet done. Check call '
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
! Examine next 3 lines!!!      
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
! End examination
C    &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     NINTER11,WORK(KINTM11),WORK(KLIHIND),IHSM,
! Add info for contractions
     &     NH22,WORK(KNH22),
     &     NH22TOT,WORK(KNH22TOT),WORK(KNH22TT),WORK(KCOMBNH22TOT),
     &     WORK(KPERMNH22TOT),WORK(KFACNH22TOT),WORK(KCOMBFACNH22TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      stop 'check this out'
      END IF
*
* H22 --12--> M10
*
* Allocate M10
C need something else than CC_AMP perhaps some scratch vector!
C Still ITSM
      IHM = 1 !Still the Hamiltonian that is contracted from
C NSPOBEX_TP_CC,WORK(KLSOBEX_CC) still the cluster operator
C WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC) also cluster operator (should
C also be made for intermediates)
C NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM still hamiltonian in
C NH22,WORK(KNH22) still H22
C Careful third argument is set to zero!
C KLIBSOBEX_CC needs to be updated for every intermediate!!!
C
      II = 0
      IFHM = 1 !integrals are fetched
      IF(II.GE.1) THEN
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,WORK(KM10),IONE,IHM,
     &     IFHM,WORK(KT_CC), 
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
! Examine next 3 lines!!!      
     &     WORK(KLIBSOBEX_M10),WORK(KLLSOBEX_M10),NINTER10,
! End examination
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
! Add info for contractions
     &     NH22,WORK(KNH22),
     &     NH22TOM10,WORK(KNH22TOM10),WORK(KNH22T10),
     &     WORK(KCOMBNH22TOM10),WORK(KPERMNH22TOM10),
     &     WORK(KFACNH22TOM10),WORK(KCOMBFACNH22TOM10),
     &     WORK(KINTM10),NINTER10,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      stop 'just check it'
      END IF
*
* H21 --10--> M11
*

*
* M11 --10--> M01 
*
* Deallocate for M11

*
* H22 --10--> M12
*
* Allocate for M02 (Not yet in use)
* Allocate for M12

*
* M12 --10--> M02
*
* Deallocate M12

*
* M02 --02--> T
*

*
* M02 --01--> M01
*

*
* H22 --21--> M01
*

*
* H22 --22--> T
*
* Starter lige her!!!
*
* Instead of CC_VEC som intermediate vector is just passed down or?
*
* CC_AMP bliver til T
*
*
* H22 --22--> T
*
      IHM = 1
      IFHM = 1 !integrals are fetched
      II = 0
      IF(II.GE.1) THEN
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
C    &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
! Examine next 3 lines!!!      
C    &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TP_CC,
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
! End examination
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
! Add info for contractions
     &     NH22,WORK(KNH22),
     &     NH22TOT,WORK(KNH22TOT),WORK(KNH22TT),WORK(KCOMBNH22TOT),
     &     WORK(KPERMNH22TOT),WORK(KFACNH22TOT),WORK(KCOMBFACNH22TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      stop 'check this out'
      END IF

*
* H21 --20--> M01
*

*
* H12 --11--> M01
*

*
* H11 --10--> M01
*

*
* H02 --01--> M01
*

*
* M01 --01--> T
*
* Deallocate M01

*
* H22 --12--> M10
*
* Allocate M10
C need something else than CC_AMP perhaps some scratch vector!
C Still ITSM
      IHM = 1 !Still the Hamiltonian that is contracted from
C NSPOBEX_TP_CC,WORK(KLSOBEX_CC) still the cluster operator
C WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC) also cluster operator (should
C also be made for intermediates)
C NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM still hamiltonian in
C NH22,WORK(KNH22) still H22
C 
C

*
* H21 --11--> M10
*

*
* H20 --10--> M10
*

*
* M10 --10--> T
*
C     stop 'about to call for intermediate contraction'
      print*,'about to contract from an intermediate '
      II = 1
      IF(II.GE.1) THEN
      IHM = 2
      IFHM = 2 !intermediates are fetched
C still needs to consider integrals/intermediates!!!!
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KM10), 
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
! Examine next 3 lines!!!      
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
! End examination Not sure if one should use hamiltonian or intermediate
! should perhaps not use KLIHIND since it is for H!!!
     &     NINTER10,WORK(KINTM10),WORK(KLIHIND),IHSM,
! Add info for contractions
     &     NINTER10,WORK(KLLSOBEX_M10),
     &     NM10TOT,WORK(KNM10TOT),WORK(KNM10TT),
     &     WORK(KCOMBNM10TOT),WORK(KPERMNM10TOT),
     &     WORK(KFACNM10TOT),WORK(KCOMBFACNM10TOT),
C    &     WORK(KSNM10TOT),WORK(KFNM10TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      stop 'after some intermediate'
      END IF
* Deallocate M10
*
* H22 --22--> T
*

*
* H21 --21--> T
*

      IHM = 1
      IFHM = 1 !integrals are fetched
      II = 0
      IF(II.GE.1) THEN
*
      CALL EMTHET_COM_KRCC(CC_AMP,ITSM,CC_VEC,IONE,IHM,
     &     IFHM,WORK(KT_CC), 
C    &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
     &     NSPOBEX_TPE,WORK(KLSOBEX_CC),WORK(KLIBSOBEX_CC),
! Examine next 3 lines!!!      
C    &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TP_CC,
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),NSPOBEX_TPE,
! End examination
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
! Add info for contractions
     &     NH21,WORK(KNH21),
     &     NH21TOT,WORK(KNH21TOT),WORK(KNH21TT),WORK(KCOMBNH21TOT),
     &     WORK(KPERMNH21TOT),WORK(KFACNH21TOT),WORK(KCOMBFACNH21TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK,KFREE,LFREE)
      stop 'check this out'
      END IF
*
* H12 --12--> T
*

*
* H20 --20--> T
*

*
* H11 --11--> T
*

*
* H02 --02--> T
*

*
* H10 --10--> T
*

*
* H01 --01--> T
*

*
* H00 --00--> T
*

* End of the line

*
      E_CC = CC_VEC(N_CC_AMP+1)
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
C?    WRITE(6,*) ' I12 at END of CC_VEC = ', I12
*. 
C     WRITE(6,*) ' Enforced stop after CC_VEC_FNC'
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
