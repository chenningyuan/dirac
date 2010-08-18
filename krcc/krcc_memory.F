      SUBROUTINE ALLOC_KRCC(WORK,KFREE,LFREE)
*
* Dimensions and
* Allocation of static memory
*
* =====
* Input
* =====
*
* KFREE : Pointer to first element of free space
* Information in /LUCINP/,/ORBINP/,/CSYM/
*
* ======
* Output
* ======
* KFREE : First array of free space after allocation of
*         static memory
* /GLBBAS/,/CDIM/
*
*
* =======
* Version
* =======
*
* 2010 version by Lasse will replace all by allocate!
*
*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "csm.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cintfo.inc"
#include "ctcc.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

*.Output
#include "glbbas.inc"
CTF adding global info
      common/glbref/ISPC/ISM/ITERGLB
*
C#include "memint.h"
*
* New allocateable array. SHOULD CHECK IF IT IS OKAY TO JUST ALLOCATE ARRAYS FOR MATRICES
*
* Integrals
*
C     REAL,ALLOCATABLE,SAVE,TARGET :: TINT1(:),TINT2(:)
*
* Fock matrices. NOT SURE TO SAVE THESE
*
C     REAL,ALLOCATABLE,SAVE,TARGET :: TFI(:), TFI_AL(:), TFI_BE(:)
C     REAL,ALLOCATABLE,SAVE :: FI_IM(:), FI_AL_IM(:), FI_BE_IM(:) !not in use
*
* Other stuff
*
C     REAL,ALLOCATABLE,SAVE,TARGET :: TPINT1(:),TPINT2(:)
C     REAL,ALLOCATABLE,SAVE,TARGET :: TLSM1(:),TLSM2(:)
*
* The pointers
*
C     REAL,POINTER :: ZINT1(:),ZINT2
C     REAL,POINTER :: ZFI,ZFI_AL,ZFI_BE
C     REAL,POINTER :: ZPINT1,ZPINT2,ZLSM1,ZLSM2
*
* Now start the allocating
*
C     ALLOCATE(TINT1(2*NTOOB**2))
C     ALLOCATE(TINT2(NINT2))
*
C     ALLOCATE(TFI(2*NTOOB**2))
C     ALLOCATE(TFI_AL(NTOOB**2))
C     ALLOCATE(TFI_BE(NTOOB**2))
*
C     ALLOCATE(TPINT1(NBINT1))
C     ALLOCATE(TPINT2(NBINT2))
C     ALLOCATE(TLSM1(NBINT1))
C     ALLOCATE(TLSM2(NBINT2))
*
* And then point to them
*
C     ZINT1 => TINT1
C     ZINT2 => TINT2(1)
*
C     ZFI => TFI(1)
C     ZFI_AL => TFI_AL(1)
C     ZFI_BE => TFI_BE(1)
*
C     ZPINT1 => TPINT1(1)
C     ZPINT2 => TPINT2(1)
C     ZLSM1 => TLSM1(1)
C     ZLSM2 => TLSM2(1)
C     DO I=1,2
C       TINT1(I) = I
C       print*,'TINT1,I',TINT1(I),I
C     END DO
C     print*,'ZINT1',ZINT1
C
C     print*,'test'
C     print*,'ZINT1',ZINT1(1)
C     print*,'ZINT1',ZINT1(2)
C     CALL DUM_ROUT(TINT1)
C     STOP
C
C allocating only the nessecary parts. Lasse
C
C Note KINT2 is needed for T1 transformed hamiltonian
*.1 : One electron integrals( Complete matrix allocated )
C       call memmar(KINT1,2*NTOOB ** 2,'ADDS  ',2,'INT1   ')
*.1.1 : Inactive fock matrix
C       call memmar(KFI  ,2*NTOOB ** 2,'ADDS  ',2,'FI     ')
*LS start
*.1.2 Inactive Fock matrix for barred and unbarred, real
C       CALL MEMMAR(KFI_AL,NTOOB**2, 'ADDS  ',2,'FI_AL ')
C       CALL MEMMAR(KFI_BE,NTOOB**2, 'ADDS  ',2,'FI_BE ')
*.1.2 Inactive Fock matrix for barred and unbarred, complex
C       CALL MEMMAR(KFI_AL_IM,NTOOB**2, 'ADDS  ',2,'FI_AL_IM ')
C       CALL MEMMAR(KFI_BE_IM,NTOOB**2, 'ADDS  ',2,'FI_BE_IM ')
*LS end
*.2 : Two electron integrals
C       IF(NOINT.EQ.0.AND.INCORE.EQ.1) THEN
C         call memmar(KINT2,NINT2,'ADDS  ',2,'INT2   ')
C       END IF
*. Pointers to symmetry block of integrals
C       call memmar(KPINT1,NBINT1,'ADDS  ',2,'PINT1   ')
C       call memmar(KPINT2,NBINT2,'ADDS  ',2,'PINT2   ')
*. Symmetry of last index as a function of initial index
C       call memmar(KLSM1,NBINT1,'ADDS  ',2,'LSM1   ')
C       call memmar(KLSM2,NBINT2,'ADDS  ',2,'LSM2   ')
*
* Using dirac array
*
C     CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KDAX,NORBT*NORBT*NZ,WORK,KFREE,LFREE)
C     SUBROUTINE MEMGET(DATA_TYPE,KBASE,LENGTH,WORK,KFREE,LFREE)
* One and two electron integrals
      CALL MEMGET('REAL',KINT1,2*NTOOB**2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KINT2,NINT2,WORK,KFREE,LFREE)
* Inactive fock matrix
      CALL MEMGET('REAL',KFI,2*NTOOB**2,WORK,KFREE,LFREE)
* Inactive Fock matrix for barred and unbarred, real
      CALL MEMGET('REAL',KFI_AL,NTOOB**2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFI_BE,NTOOB**2,WORK,KFREE,LFREE)
*. Pointers to symmetry block of integrals
      CALL MEMGET('REAL',KPINT1,NBINT1,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPINT2,NBINT2,WORK,KFREE,LFREE)
*. Symmetry of last index as a function of initial index
      CALL MEMGET('REAL',KLSM1,NBINT1,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KLSM2,NBINT2,WORK,KFREE,LFREE)
*
* =============================
*. For double group integrals
* =============================
*
*.4. Integrals (ij!kk) -(ik!kj), pointer and space
*
C Notice this array is in use for T1 transformed hamiltonian. Lasse
C     IF (SIGDEN_ROUTE .ne. 'SIGDEN') THEN
C     call memmar(KINH1,NTOOB*NTOOB,'ADDS  ',2,'KINH1  ')
C     END IF
*
        KMOAOIN = 0
        KMOAOUT = 0
        KSPSAOIN = 0
        KSPSAOUT = 0
        KISPSAOIN = 0
        KISPSAOUT = 0
*. (more to come )
      RETURN
      END
*
      SUBROUTINE MEMCHK_KRCC(WORK)
*
*     Lazy mans memchk for entire memory, Lasse
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "glbbas.inc"
*
      CALL MEMCHK('cemory',WORK,KINT1)
*
      RETURN
      END
