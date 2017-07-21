C  /* Deck get_bed */                                                  
      SUBROUTINE GET_BED(OMEGA,ISSYM,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Build representation matrix of full light-matter interaction
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC which is open with unit number LUXVC
C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,DP5=0.5D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP
      DIMENSION OMEGA(MAXEXC,NBSYM),ISSYM(MAXEXC,NBSYM),WORK(*)
      real(8), allocatable              :: PRPINT(:,:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:),IPSYM2INT(:,:)
      real(8)                           :: PBUF(3,2)
      real(8), allocatable              :: PMAT(:,:)
      logical, allocatable              :: first(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(8))
      allocate(INTREP(8))
      allocate(IPSYM2INT(NBSYM,2))
      allocate(first(nz))
C
C     Integral directives
C     
      CALL PR1DIR('CXIKR   ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      WRITE(6,*) 'NBAST= ',NBAST,'  NELMNT= ',NELMNT
      WRITE(6,*) 'Offsets ppr/lbl:',NPRPS,NPRPLBL,NZ
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      allocate(PMAT(N2BBASX,NZ))
      IORDER=0
      NPQUAD=0
C
C     Operator setup
C
      IOPTYP  = 8                ! iA.P           * dot-product of alpha and vector operator
      NMAT    = MCMP(IOPTYP)     ! number of components, should be three
      IF((NPRPS+2).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &   ' ** ERROR in GET_BED **  Redimension MAXPRPS to ',NPRPS+2
        CALL QUIT('GET_BED: Too may properties !')
      ENDIF
      DO IOP = 1,2
         PRPNAM(NPRPS+IOP)(1:3) = 'BED'
         WRITE(PRPNAM(NPRPS+IOP)(4:5),'(I1,A1)') IOP,'_'
         IPRPTIM(NPRPS+IOP) = -1 ! Time antisymmetric operator
      ENDDO
C      
C     Find what symmetries are available
C
      NCOMP=NPRPINT/2
      CALL ICOPY(2*NBSYM,-1,0,IPSYM2INT,1)
      DO I = 1,NCOMP
        ISYM = INTREP(I) + 1
        IPSYM2INT(ISYM,1) = I
        II = NCOMP + I
        ISYM = INTREP(II) + 1
        IPSYM2INT(ISYM,2) = II
      ENDDO
      WRITE(6,*) 'IPSYM2INT:'
      DO ISYM = 1,NBSYM
         WRITE(6,*) IPSYM2INT(ISYM,1),IPSYM2INT(ISYM,2)
      ENDDO
C     
C     Loop over excitation energies
C
      DO 10 IOPSY = 1,NFSYM
        NBF = NBSYM/NFSYM
        DO 20 IBF = 1,NBF
          ISYM   = JFSYM(IBF,IOPSY)
          IREP   = ISYM-1
          NEXCNV = KEXCNV(JSYMOP)
          IF (NEXCNV .EQ. 0) GO TO 20
          CALL PRSYMB(LUPRI,'=',80,0)
          WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',ISYM,' : ',REP(IREP)
C
C         Find contributions to property matrix
C
          DO IMAT = 1,NMAT
            IAIND = JM4(IMAT,IOPTYP)
            IAREP = JM4REP(IAIND)
            IBSYM = IBTXOR(IAREP,IREP)+1
            DO IOP = 1,2
              IF(IPSYM2INT(IBSYM,IOP).EQ.-1) THEN
                PBUF(IMAT,IOP) = D0
              ELSE
                PBUF(IMAT,IOP)          = DP5*UPOL(IMAT)
                IPRPLBL(IMAT,NPRPS+IOP) = IPSYM2INT(IBSYM,IOP)+NPRPLBL
              ENDIF
            ENDDO
          ENDDO
          DO IOP = 1,2
             PRPNAM(NPRPS+IOP)(6:8) = REP(IREP)
             IPRPSYM(NPRPS+IOP)     = ISYM ! Total symmetry of operator
             IF(DSUM(NMAT,PBUF(1,IOP),1).EQ.D0) THEN
                IPRPTYP(NPRPS+IOP) = -1  ! Means that there is no contribution
             ELSE
                IPRPTYP(NPRPS+IOP) = IOPTYP
             ENDIF
          ENDDO
          DO 30 IEXC = 1,KEXCNV(ISYM)
             IF(ABS(OMEGA(IEXC,ISYM)).LT.1.0D-3) THEN
                WRITE(LUPRI,*)
     &          'GET_BED: Skipping excitation. Frequency too small: ',
     &          OMEGA(IEXC,ISYM)
                GOTO 30
             ENDIF
C
C           Length of wave vector
C 
            FAC=OMEGA(IEXC,ISYM)/CVAL
            DO IOP = 1,2
               WRITE(PRPNAM(NPRPS+IOP)(9:16),'(A4,I4)') '_exc',IEXC
               DO IMAT = 1,NMAT
                  FACPRP(IMAT,NPRPS+IOP) = PBUF(IMAT,IOP)/FAC
               ENDDO
            ENDDO
C
C           Debug output
C
            DO IOP = 1,2
               WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP),FAC
               IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
                 WRITE(6,*) ' - skipped'
               ELSE
                 DO IMAT = 1,NMAT
                   II = IPRPLBL(IMAT,NPRPS+IOP)
                   WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
                 ENDDO
               ENDIF
            ENDDO
C     
C           Generate integrals
C
            DO I = 1,3
              WAVEVEC(I) = FAC*UWAVE(I)
            ENDDO
            CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &           NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &           TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &           0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
            IF(IPRINT.GE.5) THEN
              DO ICOMP = 1,NPRPINT
                WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &            '. Symmetry: ',REP(INTREP(NCOMP))
                CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
              ENDDO
            ENDIF
            DO IOP = 1,2
              INDXPR = NPRPS+IOP
              IF(IPRPTYP(INDXPR).NE.-1) THEN
                 CALL PRPMSAO(INDXPR,PMAT,.FALSE.,PRPINT,FIRST,
C     &                   WORK(KFREE),LFREE,IPRINT)
     &                WORK(KFREE),LFREE,10)
              ENDIF
            ENDDO
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
C
C     Memory deallocation
      deallocate(PMAT)
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(IPSYM2INT)
      deallocate(first)
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
      PARAMETER (D0=0.0D0,DI2=0.5D0,DI4=0.25D0,D2=2.0D0)
      PARAMETER (MAXNUM=1455)
      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)
      DIMENSION OMEGA(MAXEXC,NBSYM),ISSYM(MAXEXC,NBSYM),WORK(*)
      REAL(8):: KVEC(3,MAXNUM),E1VEC(3),E2VEC(3),FAC,TMPNRM
      REAL(8):: E1E1(3,3),E1E2(3,3),OSC(MAXEXC),ECD(MAXEXC)

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
C       CALL INTEGRAL ROUTINE

C       Compute theta and phi angles
        PHI=ATAN2(Y(IK),X(IK))
        THETA=ACOS(Z(IK))
C
C.......Weights for chi-average over pol-directions
C       Linear absorption 
C       <E1VEC,E1VEC>_chi average
        E1E1(1,1)=DI2*((COS(THETA)*COS(PHI))**2+SIN(PHI)**2)
        E1E1(1,2)=-DI4*SIN(THETA)**2*SIN(2*PHI)
        E1E1(1,3)=-DI4*SIN(2*THETA)**2*COS(PHI)
        E1E1(2,1)=E1E1(1,2)
        E1E1(2,2)=DI2*((COS(THETA)*SIN(PHI))**2+COS(PHI)**2)
        E1E1(2,3)=-DI4*SIN(2*THETA)*SIN(PHI)
        E1E1(3,1)=E1E1(1,3)
        E1E1(3,2)=E1E1(2,3)
        E1E1(3,3)=DI2*SIN(THETA)**2
C       Circular dichroism
C       <E1VEC,E2VEC>_chi average      
        E1E2(1,1)=D0
        E1E2(1,2)=DI2*(SIN(PHI)**2-COS(PHI)**2)
        E1E2(1,3)=-DI2*SIN(THETA)*SIN(PHI)
        E1E2(2,1)=-E1E2(1,2)
        E1E2(2,2)=D0
        E1E2(2,3)=-DI2*COS(PHI)*SIN(THETA)
        E1E2(3,1)=-E1E2(1,3)
        E1E2(3,2)=-E1E2(2,3)
        E1E2(3,3)=D0
C
C     AVERAGE OVER POL-DIRECTION      
C     WRITE(LUPRI,*) 'KVEC', KVEC(1,IK), KVEC(2,IK), KVEC(3,IK), WT(IK)
      WRITE(LUPRI,*) 'VEC', X(IK), Y(IK), Z(IK), WT(IK)
      WRITE(LUPRI,*) 'PHI, THETA: ', PHI, THETA
      DO I=1,3
        WRITE(LUPRI,*)  E1E1(I,1), E1E1(I,2), E1E1(I,3)
      END DO
C
C      DO I=1,3
C        DO J=I,3
C          IF I.EQ.J THEN
C            OSC(IEXC)=OSC(IEXC)+WT(IK)*E1E1(I,I)
C     &                         *(TMR(I)*TMR(I)+TMI(I)*TMI(I))
C          ELSE
C            OSC(IEXC)=OSC(IEXC)+D2*WT(IK)*E1E1(I,J)
C     &                         *(TMR(I)*TMR(J)-TMI(I),TMI(J))
C            ECD(IEXC)=ECD(IEXC)+WT(IK)*E1E2(I,J)
C     &                         *(TMI(I)*TMR(J)+TMR(I)*TMI(J))
C          END IF 
C        END DO 
C      END DO
C     
      END DO
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCC
