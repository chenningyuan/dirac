CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck sleb */                                                     C
      SUBROUTINE SLEB(X,Y,Z,WT,NSCHEME,NCOUNT,NATOM,RNODE,RBRAGG)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"
#include "pi.h"

      PARAMETER (MAXNUM=1455)
      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)

      REAL KVEC, LKVEC, e1, e2, POLVEC

      CALL LD0038(X,Y,Z,WT,NCOUNT)

      LKVEC=EIGVAL(J)/SPEED_OF_LIGHT      

C     Creating wave vector for the given Lebedev point or as specified
C     in the input      
C
      KVEC(1,I)=LKVEC * X(I) 
      KVEC(2,I)=LKVEC * Y(I) 
      KVEC(3,I)=LKVEC * Z(I) 
C       
C     The associated unit vectors for the polarization direction
C
C     The first unit vector E1 orthogonal to KVEC
C
      if (KVEC(1,I).eq.0.0.and.KVEC(2,I).eq.0.0) then
         E1(1) = 1.0d0 
         E1(2) = 0.0d0 
         E1(3) = 0.0d0 
      else
         E1(1) = KVEC(2,I)
         E1(2) =-KVEC(1,I)
         E1(3) = 0.0d0
      end if 
      TMPNRM = sqrt(E1(1)**2 + E1(2)**2 + E1(3)**2)
      E1(1) = E1(1)/TMPNRM
      E1(2) = E1(2)/TMPNRM
      E1(3) = E1(3)/TMPNRM
C
C     Compute the second unit vector E2 orthogonal to KVEC as the cross
C     product of KVEC and E1
C
      E2(1) = KVEC(2,I)*E1(3) - KVEC(3,I)*E1(2)
      E2(2) = KVEC(3,I)*E1(1) - KVEC(1,I)*E1(3)
      E2(3) = KVEC(1,I)*E1(2) - KVEC(2,I)*E1(1)
      TMPNRM = sqrt(E2(1)**2 E2(2)**2 + E2(3)**2)   
      E2(1) = E2(1)/TMPNRM
      E2(2) = E2(2)/TMPNRM
      E2(3) = E2(3)/TMPNRM
C    
C     An arbitrary polarization vector can be constructed as
C     cos(chi)*E1+sin(chi)*E2.
C
C     
      TOFILE = .FALSE.
C
C
C     Make routine for making labels for the integrals. The number of
C     symmetries depend on the 
C
C     Compute the EXPIKR integrals (real, imag) for each ireducible
C     representation
C
C     Make a routine (if not already there) that based on a given
C     k-vector generates the corresponding integrals and returns them in
C     an array

      GET_EXP_INT(KVEC,NFREQ,OFFSET,WORK,LWRK,TOFILE)


