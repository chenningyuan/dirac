C  /* Deck get_bed */                                                  C
      SUBROUTINE GET_BED(OMEGA,ISSYM,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Build representation matrix of full light-matter interaction
C     Scalar integrals are sorted on symmetries (max. 8)
C     
C      
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC which is open with unit number LUXVC                                  C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
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
      real(8)                           :: WAVVEC(3),POLVEC(3)
      character (len=6)                 :: PM4NAM(3)
      integer                           :: IM4(3)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
!     We are building the full light-matter interaction operator      
      DATA IM4/2,3,4/
      DATA PM4NAM/'XALPHA','YALPHA','ZALPHA'/
C
#include "ibtfun.h"
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(8))
      allocate(INTREP(8))
C
C     Integral directives
C     
      CALL PR1DIR('CXIKR   ',INTTYP,NCOMP,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
C     NCOMP is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      DO ICOMP = 1,NCOMP
          WRITE(6,*) LABINT(ICOMP),'. Symmetry: ',
     &               REP(INTREP(ICOMP))
      ENDDO
      WRITE(6,*) 'NBAST= ',NBAST,'  NELMNT= ',NELMNT
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on triangular form)
      allocate(ONEINT(NELMNT,NCOMP))
      IORDER=0
      NPQUAD=0
C
C     Define propagation vector
C
      WAVVEC(1)=0.0
      WAVVEC(2)=0.0
      WAVVEC(3)=1.0
C
C     Define polarization vector
C
      POLVEC(1)=1.0
      POLVEC(2)=0.0
      POLVEC(3)=0.0
C
C     Check for orthogonality of propagation and polarization vector
C
      TTOL=1.0D-10
      TMP=DDOT(3,WAVVEC,1,POLVEC,1)
      IF(ABS(TMP).GT.TTOL) THEN
         CALL QUIT
     & ('GET_BED: Polarization not orthogonal to propagation !')
      ENDIF
C
C     Loop over excitation energies
C
      DO 10 IOPSY = 1,NFSYM
        NBF = NBSYM/NFSYM
        DO 20 IBF = 1,NBF
          ISYM = JFSYM(IBF,IOPSY)
          IREP = ISYM-1
          NEXCNV = KEXCNV(JSYMOP)
          IF (NEXCNV .EQ. 0) GO TO 20
          CALL PRSYMB(LUPRI,'=',80,0)
          WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',ISYM,' : ',REP(IREP)

C
C         Find contributions to property matrix
C
          DO IMAT = 1,3
            IAIND = JM4(1,IOPTYP)
            IAREP = JM4REP(IAIND)
            IBREP = IBTXOR(IAREP,IREP)
            WRITE(6,*) '* Contributions to ',
     &         PM4NAM(IMAT),' ',REP(IAREP),' from:'
            DO ICOMP = 1,NCOMP
            IF(INTREP(ICOMP).EQ.IBREP) THEN
              WRITE(6,*) '            ',LABINT(ICOMP),'. Symmetry: ',
     &               REP(INTREP(ICOMP))
            ENDIF
            ENDDO
          ENDDO
          DO IEXC = 1,KEXCNV(ISYM)
C
C           Length of wave vector
C 
            FAC=OMEGA(IEXC,ISYM)/CVAL
C
C           Generate integrals
C
            DO I = 1,3
              ORIGIN(I) = FAC*WAVVEC(I)
            ENDDO
CTROND
            WRITE(6,*) '*Wave..',REP(IREP),IEXC,OMEGA(IEXC,ISYM),FAC,
     &          (ORIGIN(J),J=1,3)
            CALL PR1IN2(WORK(KFREE),LFREE,ONEINT,INTTYP,
     &           NCOMP,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &           TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &           0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
            IF(IPRINT.GE.5) THEN
              DO ICOMP = 1,NCOMP
                WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &            '. Symmetry: ',REP(INTREP(NCOMP))
                CALL OUTPAK(ONEINT(1,ICOMP),NBAST,1,LUPRI)
              ENDDO
            ENDIF

          ENDDO
 20     CONTINUE
 10   CONTINUE
C
C     Memory deallocation
      deallocate(ONEINT)
      deallocate(LABINT)
      deallocate(INTREP)
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
