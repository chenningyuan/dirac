      SUBROUTINE CC_VEC_FNC_REL(CC_AMP,CC_VEC,E_CC,
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
      COMMON/CMXCJ/MXCJ
*. Local scratch 
      INTEGER ISX(2*MXPNGAS*MXPNGAS)
*
      REAL*8 INPRDD, INPROD
*
      LBLK = -1
      NTEST = 00
      NTEST = MAX(NTEST,IPRCC)
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
*      May want to allocate some arrays
*
      CALL EXP_MT_H_EXP_TC_MEM_KRCC(WORK,KFREE,LFREE)
*
*. Start by zeroing the CC-vector function 
      ZERO = 0.0D0
      CALL SETVEC(CC_VEC,ZERO,N_CC_AMP+1)
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
* Allocate for M01 (Not yet in use)
* Allocate for M11
      NPART = 1
      NHOLE = 1
*
      DO IH = 1,NH22
        IHTP = WORK(KNH22 +IH) !Not so pretty
        CALL CONTRACTION_ANALYZE(NPART,NHOLE,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM11),
     &                           WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                           WORK,KFREE,LFREE)
      END DO
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
      RETURN
      END

*
