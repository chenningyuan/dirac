C  /* Deck get_bed */                                                  C
      SUBROUTINE GET_BED(OMEGA,ISSYM,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC which is open with unit number LUXVC                                                                 
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
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
      DATA PM4NAM/'XALPHA','YALPHA','ZALPHA'/
C
#include "ibtfun.h"
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
      DO ICOMP = 1,NCOMP
          WRITE(6,*) LABINT(ICOMP),'. Symmetry: ',
     &               REP(INTREP(ICOMP))
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
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
C     Operator type
      IOPTYP = 8 ! iA.P
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
C  /* Deck get_bed2 */                                                  C
      SUBROUTINE GET_BED2(OMEGA,ISSYM,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC which is open with unit number LUXVC                                                                 
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
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
      DATA PM4NAM/'XALPHA','YALPHA','ZALPHA'/
C
#include "ibtfun.h"
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
      DO ICOMP = 1,NCOMP
          WRITE(6,*) LABINT(ICOMP),'. Symmetry: ',
     &               REP(INTREP(ICOMP))
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
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
C     Operator type
      IOPTYP = 8 ! iA.P
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bedleb */                                                     C
      SUBROUTINE BEDLEB(OMEGA,ISSYM,IPRINT,WORK, KFREE,LFREE)
C      SUBROUTINE BEDLEB(X,Y,Z,WT,NSCHEME,NCOUNT,NATOM,RNODE,RBRAGG)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"
#include "pi.h"
C
#include "dcbgen.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "orgcom.h"
C
      PARAMETER (DI2=0.5D0,DI4=0.25D0,MAXNUM=1455)
      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)
      DIMENSION OMEGA(MAXEXC,NBSYM),ISSYM(MAXEXC,NBSYM),WORK(*)
      REAL(8):: KVEC(3,MAXNUM),E1VEC(3),E2VEC(3),FAC,TMPNRM
      REAL(8):: E1E2(3,3)

C      IF LEBEDEV THEN
        CALL LD0014(X,Y,Z,WT,NCOUNT)
C      END IF
      FAC=OMEGA(1,1)/CVAL
C
      DO IK=1,NCOUNT
C   
C     Creating wave vector for the given Lebedev point or as specified
C     in the input      
C
        KVEC(1,IK)=FAC * X(IK) 
        KVEC(2,IK)=FAC * Y(IK) 
        KVEC(3,IK)=FAC * Z(IK) 
CC       
CC     The associated unit vectors for the polarization direction
CC     The first unit vector E1 orthogonal to KVEC
CC
C        IF (KVEC(1,IK).EQ.0.0.AND.KVEC(2,IK).EQ.0.0) THEN
C           E1VEC(1) = 1.0d0 
C           E1VEC(2) = 0.0d0 
C           E1VEC(3) = 0.0d0 
C        ELSE
C           E1VEC(1) = KVEC(2,IK)
C           E1VEC(2) =-KVEC(1,IK)
C           E1VEC(3) = 0.0d0
C        END IF 
C        TMPNRM = SQRT(E1VEC(1)**2+E1VEC(2)**2+E1VEC(3)**2)
C        E1VEC(1) = E1VEC(1)/TMPNRM
C        E1VEC(2) = E1VEC(2)/TMPNRM
C        E1VEC(3) = E1VEC(3)/TMPNRM
CC
CC     Compute the second unit vector E2 orthogonal to KVEC as the cross
CC     product of KVEC and E1
CC
C        E2VEC(1) = KVEC(2,IK)*E1VEC(3)-KVEC(3,IK)*E1VEC(2)
C        E2VEC(2) = KVEC(3,IK)*E1VEC(1)-KVEC(1,IK)*E1VEC(3)
C        E2VEC(3) = KVEC(1,IK)*E1VEC(2)-KVEC(2,IK)*E1VEC(1)
C        TMPNRM = SQRT(E2VEC(1)**2+E2VEC(2)**2+E2VEC(3)**2)   
C        E2VEC(1) = E2VEC(1)/TMPNRM
C        E2VEC(2) = E2VEC(2)/TMPNRM
C        E2VEC(3) = E2VEC(3)/TMPNRM
CC      
C       An arbitrary polarization vector can be constructed as
C       cos(chi)*E1+sin(chi)*E2.
C 
C       CALL INTEGRAL ROUTINE

C       Compute theta and phi angles
        PHI=ATAN2(Y(IK),X(IK))
        THETA=ACOS(Z(IK))
C
C       <E1VEC,E2VEC>_chi average
        E1E2(1,1)=DI2*(COS(THETA)**2+COS(PHI)**2+SIN(PHI)**2)
        E1E2(1,2)=-DI4*SIN(THETA)**2*SIN(2*PHI)
        E1E2(1,3)=-DI4*SIN(2*THETA)**2*COS(PHI)
        E1E2(2,1)=E1E2(1,2)
        E1E2(2,2)=DI2*(COS(THETA)**2*SIN(PHI)**2+COS(PHI)**2)
        E1E2(2,3)=-DI4*SIN(2*THETA)*SIN(PHI)
        E1E2(3,1)=E1E2(1,3)
        E1E2(3,2)=E1E2(2,1)
        E1E2(3,3)=DI2*SIN(THETA)**2
C
C       AVERAGE OVER POL-DIRECTION      
C      WRITE(LUPRI,*) 'KVEC', KVEC(1,IK), KVEC(2,IK), KVEC(3,IK), WT(IK)
      WRITE(LUPRI,*) 'VEC', X(IK), Y(IK), Z(IK), WT(IK)
      WRITE(LUPRI,*) 'PHI, THETA: ', PHI, THETA
      END DO
C     
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCC
