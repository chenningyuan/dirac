CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck get_bed */                                                  C
      SUBROUTINE GET_BED(OMEGA,ISSYM,IPRINT,WORK,KFREE,LFREE)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "orgcom.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP
      DIMENSION OMEGA(MAXEXC,NBSYM),ISSYM(MAXEXC,NBSYM),WORK(*)
      real(8), allocatable              :: ONEINT(:,:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
      LOGICAL PROPRI,TOFILE
C
C
C     Integral information
C
C      TOFILE = .FALSE.
C      IPRINT = 100
C     DOINT: 0++0 - LS, SL components are computed and multiplied by +1
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
C
C     Generate integrals
C
      NCOMP=8
      allocate(LABINT(NCOMP))
      allocate(INTREP(NCOMP))
C
C     Integral directives
C     
      CALL PR1DIR('CXIKR   ',INTTYP,NCOMP2,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
      IF(NCOMP2.NE.NCOMP) THEN
        WRITE(6,*) 'NCOMP - expected :',NCOMP
        WRITE(6,*) 'NCOMP - returned :',NCOMP2 
        CALL QUIT('GET_BED: NCOMP differs !')
      ENDIF
CTROND
      WRITE(6,*) 'BASINFO:',NBAST,NELMNT
      DO I = 1,8
        WRITE(6,*) 'INTREP..',INTREP(I)
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
      allocate(ONEINT(NELMNT,NCOMP))
      IORDER=0
      NPQUAD=0
      TOFILE=.FALSE.

C
C     Generate integrals
C
      ORIGIN(1)=0.0
      ORIGIN(2)=0.0
      ORIGIN(3)=1.0
      DO I = 1,3
        ORIGIN(I) = ORIGIN(I)/CVAL
      ENDDO
      CALL RGETINT('CXIKR   ',ONEINT,DOINT,NCOMP,IORDER,NPQUAD,
     &                  'TRIANG',TOFILE,20,WORK(KFREE),LFREE)
      CALL QUIT('End of test')
      PROPRI=.FALSE.
      CALL PR1IN2(WORK(KFREE),LFREE,ONEINT,INTTYP,
     &           NCOMP,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
CTROND     &           TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &           TRASPH,PROPRI,SOLVNT,IORDER,NPQUAD,DUMMY,
CTROND     &           0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
     &           0,NBAST,NELMNT,20,DOINT,.FALSE.)
CTROND      IF(IPRINT.GE.5) THEN
C        DO ICOMP = 1,NCOMP
C          WRITE(6,*) '* Integrals of ',LABINT(ICOMP),'. Symmetry: ',
C     &               REP(INTREP(NCOMP))
C          CALL OUTPAK(ONEINT(1,ICOMP),NBAST,1,LUPRI)
C        ENDDO
CTROND      ENDIF
C
C     Memory deallocation
      deallocate(ONEINT)
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
