CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck get_bed */                                                  C
      SUBROUTINE GET_BED(OMEGA,ISSYM,WORK,KFREE,LFREE)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpp.h"
#include "orgcom.h"
C#include "cbiher.h"
      LOGICAL TOFILE, DOINT(2,2)
      INTEGER LWORK, KONEIN, NCOMP, NOPTP
      DIMENSION OMEGA(MAXEXC,NBSYM),ISSYM(MAXEXC,NBSYM),WORK(*)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "dcbgen.h"
#include "dcbbas.h"
C
C     Integral information
C
C      TOFILE = .FALSE.
C      IPRINT = 100
C     DOINT: 0++0 - LS, SL components are computed and multiplied by +1
C      DOINT(1,1)=.FALSE.
C      DOINT(2,1)=.TRUE.
C      DOINT(1,2)=.TRUE.
C      DOINT(2,2)=.FALSE.
C      IORDER = IDUMMY ! Order of multipole moment. Not applicable here.
C      NPQUAD = 0 ! Number of Gaussian quadrature points. 0 in our case
C      NCOMP = 8  ! Number of operators
C
C      DO IFRP=1,NFSYM
C        DO IBF=1,NBF
C           ISYM=JFSYM(IBF,IFRP) ! symmetry of exitation
C
C     Generate k vector for given excitation - for now - just along z
C
      ORIGIN(1)=0.0
      ORIGIN(2)=0.0
      ORIGIN(3)=OMEGA(1,1)/CVAL !NHL: FIXME just some arbitrary value
C
C     Generate integrals
C
      IPRINT=6
      NCOMP=8
      allocate(LABINT(NCOMP))
      allocate(INTREP(NCOMP))
C
C     NHL: changed to NOPTYP->NOPTP because the former is already
C     defined in dgroup.h
      CALL PR1DIR('CXIKR   ',INTTYP,NOPTP,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
CTROND      CALL RGETINT('CXIKR   ',WORK(KONEIN),DOINT,NCOMP,IORDER,NPQUAD,
CTROND     &                  'TRIANG',TOFILE,IPRINT,WORK(KFREE),LFREE)
C
C
C      TOFILE=.FALSE.
C
C      CALL MEMGET('REAL',KONEIN,NNBBASX*NCOMP,WORK,KFREE,LFREE)
C      CALL PR1IN2(WORK(KFREE),LFREE,WORK(KONEIN),INTTYP,NOPTP,INTREP,
C     &           ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
C     &           PROPRI,SOLVNT,IORDER,NPQUAD,DOATOM,
C     &           NATOM,NBAST,NELMNT,IPRINT,DOINT,TOFILE)
C
C     Memory deallocation
      deallocate(LABINT)
      deallocate(INTREP)
C      CALL MEMREL('GET_BED',WORK,KWORK,KONEIN,KFREE,LFREE)
C
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
