CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck get_bed */                                                  C
      SUBROUTINE GET_BED(WORK,LWORK)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL TOFILE, DOINT(4)
      INTEGER LWORK, KFREE, LFREE, KONEIN
      DIMENSION WORK(LWORK)
C
#include "dummy.h"
C#include "dcbgen.h"
#include "dcbbas.h"
C
      CALL QENTER('GET_BED')
#include "memint.h"
      TOFILE = .TRUE.
      IPRINT = 100
C
      WRITE(LUPRI,*) 'DOINT',  DOINT
      WRITE(LUPRI,*) 'IDUMMY', IDUMMY
      WRITE(LUPRI,*) 'KONEIN', KONEIN
      WRITE(LUPRI,*) 'IPRINT', IPRINT
C
C     Memory allocation
C
     CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C      SUBROUTINE RGETINT(CLASLB,OP1INT,DOINT,NCOMP,IORDER,NPQUAD,
C     &                  MXFORM,TOFILE,IPRINT,WORK,LWORK)
C
C     I dont know what NCOMP should be.      
C   
C     Generate integrals
C
      CALL RGETINT('CXIKR   ',WORK(KONEIN),DOINT,0,IDUMMY,IDUMMY,
     &                  'TRIANG',.FALSE.,IPRINT,WORK,LWORK)
C
C      IF(IPRINT.GE.6) THEN
        CALL HEADER('Small component overlap integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
C      ENDIF
C
      CALL QEXIT('GET_BED')
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                      C
CC  /* Deck bedleb */                                                     C
C      SUBROUTINE BEDLEB(X,Y,Z,WT,NSCHEME,NCOUNT,NATOM,RNODE,RBRAGG)
CC                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C#include "implicit.h"
C#include "priunit.h"
C#include "pi.h"
C
C      PARAMETER (MAXNUM=1455)
C      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)
C
C      REAL KVEC, LKVEC, e1, e2, POLVEC
C
C      CALL LD0038(X,Y,Z,WT,NCOUNT)
C
C      LKVEC=EIGVAL(J)/SPEED_OF_LIGHT      
C
CC     Creating wave vector for the given Lebedev point or as specified
CC     in the input      
CC
C      KVEC(1,I)=LKVEC * X(I) 
C      KVEC(2,I)=LKVEC * Y(I) 
C      KVEC(3,I)=LKVEC * Z(I) 
CC       
CC     The associated unit vectors for the polarization direction
CC
CC     The first unit vector E1 orthogonal to KVEC
CC
C      if (KVEC(1,I).eq.0.0.and.KVEC(2,I).eq.0.0) then
C         E1(1) = 1.0d0 
C         E1(2) = 0.0d0 
C         E1(3) = 0.0d0 
C      else
C         E1(1) = KVEC(2,I)
C         E1(2) =-KVEC(1,I)
C         E1(3) = 0.0d0
C      end if 
C      TMPNRM = sqrt(E1(1)**2 + E1(2)**2 + E1(3)**2)
C      E1(1) = E1(1)/TMPNRM
C      E1(2) = E1(2)/TMPNRM
C      E1(3) = E1(3)/TMPNRM
CC
CC     Compute the second unit vector E2 orthogonal to KVEC as the cross
CC     product of KVEC and E1
CC
C      E2(1) = KVEC(2,I)*E1(3) - KVEC(3,I)*E1(2)
C      E2(2) = KVEC(3,I)*E1(1) - KVEC(1,I)*E1(3)
C      E2(3) = KVEC(1,I)*E1(2) - KVEC(2,I)*E1(1)
C      TMPNRM = sqrt(E2(1)**2 E2(2)**2 + E2(3)**2)   
C      E2(1) = E2(1)/TMPNRM
C      E2(2) = E2(2)/TMPNRM
C      E2(3) = E2(3)/TMPNRM
CC    
CC     An arbitrary polarization vector can be constructed as
CC     cos(chi)*E1+sin(chi)*E2.
CC
CC     
C      TOFILE = .FALSE.
CC
CC
CC     Make routine for making labels for the integrals. The number of
CC     symmetries depend on the 
CC
CC     Compute the EXPIKR integrals (real, imag) for each ireducible
CC     representation
CC
CC     Make a routine (if not already there) that based on a given
CC     k-vector generates the corresponding integrals and returns them in
CC     an array
C
C      GET_EXP_INT(KVEC,NFREQ,OFFSET,WORK,LWRK,TOFILE)
C
C
