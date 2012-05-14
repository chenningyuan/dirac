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
#include "orbinp.inc"
*.Output
#include "glbbas.inc"
*
      DIMENSION WORK(*)
*
* One and two electron integrals
C     CALL MEMGET('REAL',KINT1,2*NTOOB**2,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KINT2,NINT2,WORK,KFREE,LFREE)
* Inactive fock matrix
      CALL MEMGET('REAL',KFI,2*NTOOB**2,WORK,KFREE,LFREE)
* Inactive Fock matrix for barred and unbarred, real
C     CALL MEMGET('REAL',KFI_AL,NTOOB**2,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KFI_BE,NTOOB**2,WORK,KFREE,LFREE)
*. Pointers to symmetry block of integrals
C     CALL MEMGET('REAL',KPINT1,NBINT1,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KPINT2,NBINT2,WORK,KFREE,LFREE)
*. Symmetry of last index as a function of initial index
C     CALL MEMGET('REAL',KLSM1,NBINT1,WORK,KFREE,LFREE)
C     CALL MEMGET('REAL',KLSM2,NBINT2,WORK,KFREE,LFREE)
*
* =============================
*. For double group integrals
* =============================
*
*.4. Integrals (ij!kk) -(ik!kj), pointer and space
*
*
*. (more to come )
      RETURN
      END
*
      SUBROUTINE MEMCHK_KRCC(WORK)!,KFREE,LFREE)
*
*     Lazy mans memchk for entire memory, Lasse
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "glbbas.inc"
*
      DIMENSION WORK(*)
*
      CALL MEMCHK('cemory',WORK,KLOOK)
C     CALL MEMCHK('temp',WORK,KCC1)
*
      RETURN
      END
