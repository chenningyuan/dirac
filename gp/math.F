C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsjaco */
      SUBROUTINE RSJACO(LDM,NDIM,NROWE,AMAT,EIG,IJOB,IORDER,IPACK,EVEC,
     &                  WORK,LWORK)
C***********************************************************************
C
C     Diagonalize real symmetric matrix by the Jacobi method
C     This is basically a driver for the routine JACO.
C
C     Control parameter IORDER:
C       IORDER.EQ.-2 : Sort eigenvalues in absolute descending order
C       IORDER.EQ.-1 : Sort eigenvalues in descending order
C       IORDER.EQ. 1 : Sort eigenvalues in ascending order
C       IORDER.EQ. 2 : Sort eigenvalues in absolute ascending order
C       IORDER.EQ. 0 : No sorting of eigenvalues and eigenvectors
C     Control parameter IJOB:
C               IJOB = 0        Only eigenvalues desired.
C               IJOB = 1        Eigenvalues and eigenvectors in 
C                               orthonormal basis. EVEC set to identity 
C                               matrix.NROWE is dimension of basis 
C                               equal to NDIM.
C               IJOB = 2        Eigenvalues and eigenvectors in another 
C                               basis.EVEC equal to transformation to 
C                               other basis. NROWE is dimension of 
C                               basis and may differ from NDIM.
C     Control parameter IPACK:
C               IPACK = 0   : AMAT given as square matrix
C               IPACK = 1   : AMAT given as row-packed lower triangle
C
C     Written by T.Saue - Aug 15 1996
C     Last revision : Aug 15 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AMAT(*),EIG(NDIM),EVEC(LDM,NDIM),
     &          WORK(LWORK)
C
      CALL QENTER('RSDIAG')
#include "memint.h"
C
      CALL MEMGET('REAL',KBUF ,NDIM ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KIBUF,NDIM ,WORK,KFREE,LFREE)            
C
C     Determine job
C     =============
C
      IF    (IJOB.EQ.0) THEN
        NEVEC = 0
      ELSEIF(IJOB.EQ.1) THEN
        NEVEC = NDIM
        CALL DUNIT2(EVEC,NDIM,LDM,NDIM,1)
      ELSEIF(IJOB.EQ.2) THEN
        NEVEC = NDIM
      ELSE
        CALL QUIT('RSJACO: Illegal value of IJOB !')
      ENDIF
C
      IF(IPACK.EQ.0) THEN
C
C     AMAT is a full square matrix
C     ============================
C
C       Pack lower triangle of AMAT into ATRI
C
        NNDIM = (NDIM*(NDIM+1))/2
        CALL MEMGET('REAL',KATRI,NNDIM,WORK,KFREE,LFREE)
        CALL DAMATR(NDIM,AMAT,LDM,WORK(KATRI))
C
C       Diagonalize ATRI
C
        CALL JACO2(WORK(KATRI),EVEC,NDIM,NDIM,NEVEC,NROWE,LDM,
     &            WORK(KBUF),WORK(KIBUF))
C
C       Extract eigenvalues from ATRI
C
        CALL XTRCDI(WORK(KATRI),EIG,NDIM,1)
C
      ELSEIF(IPACK.EQ.1) THEN
C
C       AMAT is a row-packed lower triangular matrix
C       ============================================
C
C       Diagonalize AMAT
C
        CALL JACO(AMAT,EVEC,NDIM,NDIM,NEVEC,
     &            WORK(KBUF),WORK(KIBUF))
C
C       Extract eigenvalues from AMAT
C
        CALL XTRCDI(AMAT,EIG,NDIM,1)
      ELSE
        CALL QUIT('RSJACO: Illegal value of IPACK !')
      ENDIF
C
C     Order eigenvalues/vectors in ascending/descending order
C
      IF    (IORDER.NE.0) THEN
        CALL ORDER3(EVEC,EIG,LDM,NDIM,NEVEC,IORDER)
      ENDIF
C
      CALL MEMREL('RSJACO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RSDIAG')
      RETURN
      END
