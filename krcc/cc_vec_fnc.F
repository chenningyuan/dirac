      SUBROUTINE CC_VEC_FNC_KRCC(CC_AMP,CC_VEC,E_CC,
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
C#include "wrkspc.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cc_exc.inc"
#include "cands.inc"
#include "cstate.inc"
#include "oper.inc"
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
      NTEST = MAX(NTEST,IPRCC)
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
* Indices are free of the ph definition
*
* H22 --11--> M11'  (Special case of M11)
*
*
* M11' --11--> T    (O11' now goes to M11)
*

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
* CCAMP bliver til T
*
*
* H22 --22--> T
C    &NH22TOT,KNH22TOT,KNH22TT,KCOMBNH22TOT,KPERMNH22TOT,
C    &KFACNH22TOT,KCOMBFACNH22TOT,KSNH22TOT,KFNH22TOT,
*
      IHM = 1
*
      CALL EMTHET_COM_KRCC(CCAMP,ITSM,CC_VEC,IONE,IHM, !Perhaps sym cc_vec?
     &     NSPOBEX_TP_CC,WORK(KLSOBEX_CC),
! Examine next 3 lines!!!      
     &     WORK(KLIBSOBEX_CC),WORK(KLLSOBEX_CC),
     &     NOBEX_TP,WORK(KIBSOX_FOR_OX),WORK(KNSOX_FOR_OX),
     &     WORK(KISOX_FOR_OX),
! End examination
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIHIND),IHSM,
     &     ILTTWO,I2DIATP,
! Add info for contractions
     &     NH22,WORK(KNH22),
     &     NH22TOT,WORK(KNH22TOT),WORK(KNH22TT),WORK(KCOMBNH22TOT),
     &     WORK(KPERMNH22TOT),WORK(KFACNH22TOT),WORK(KCOMBFACNH22TOT),
     &     WORK(KSNH22TOT),WORK(KFNH22TOT),
     &     WORK(KLSOBEX_CC),NSPOBEX_TPE,WORK(KLSOBEX_CC),NSPOBEX_TPE)
      stop 'enough for now'



C     CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
C    &NH22TOT,WORK(KNH22TOT),
C    &WORK(KNH22TT),WORK(KCOMBNH22TOT),
C    &WORK(KPERMNH22TOT),
C    &WORK(KFACNH22TOT),
C    &WORK(KCOMBFACNH22TOT),
C    &WORK(KSNH22TOT),WORK(KFNH22TOT),
C    &WORK(KLSOBEX_CC),NSPOBEX_TPE,
C    &WORK(KLSOBEX_CC),NSPOBEX_TPE,
C    &WORK(KLSOBEX),IHM)
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

*
* H21 --11--> M10
*

*
* H20 --10--> M10
*

*
* M10 --10--> T
*
* Deallocate M10

*
* H22 --22--> T
*

*
* H21 --21--> T
*

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
        CALL WRT_CC_VEC2(CC_VEC,6,CCTYPE)
        WRITE(6,*) ' Initial amplitudes '
        CALL WRT_CC_VEC2(CC_AMP,6,CCTYPE)
      END IF
C?    WRITE(6,*) ' I12 at END of CC_VEC = ', I12
*. 
C     WRITE(6,*) ' Enforced stop after CC_VEC_FNC'
C      STOP ' Enforced stop after CC_VEC_FNC'
      RETURN
      END
*
      SUBROUTINE EXP_MT_H_EXP_TC_MEM_KRCC(WORK,KFREE,LFREE)
*
* Allocate the arrays which are always used 
*
*
* Big arrays for the intermediates are at the moment allocated here and
* are all stored in memory.
*

*
      RETURN
      END
