module localization
implicit none

public :: ORBLOC

private :: ORBLP1, ORBLP2, ORBLO1, ORBLO2, ORBLO4, ORBLO5, ROTCMO_MOD

contains

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLOC
!***********************************************************************
!
!     Molecular orbital localization
!
!     Written by S. Dubillard - Jan 20 2004
!     Last revision:
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbloc.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
      PARAMETER(D0 = 0.0D0)
      LOGICAL TOBE,FLAG
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12,DAYTID*24

      INTEGER*4, ALLOCATABLE, DIMENSION(:) :: NFNUC, ISEL, INDAO, INDAOS, INDNAO
      REAL*8,    ALLOCATABLE, DIMENSION(:) :: SCOEF, SCOEFO, SMAT, CMO
      REAL*8,    ALLOCATABLE :: PMAT(:,:,:,:), XKAPPA(:,:,:), HAMAT(:,:,:)

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLOC')

      IFRP = 1

      CALL TITLER('Molecular orbital localization','*',122)

! Check if coefficients are on file

      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE)THEN
        WRITE(LUPRI,'(A)') 'ORBLOC: Coefficient file not found !'
        CALL QUIT('ORBLOC: Coefficients not found !')
      ENDIF


! Memory allocation

      CALL ALLOC(SCOEF ,NISH(IFRP)*NFBAS(IFRP,0)*4,ID='SCOEF IN ORBLOC')
      CALL ALLOC(SCOEFO,NISH(IFRP)*NFBAS(IFRP,0)*4,ID='SCOEFO IN ORBLOC')
      CALL ALLOC(SMAT  ,N2BBASX,ID='SMAT IN ORBLOC')
      CALL ALLOC(CMO   ,N2BBASXQ,ID='CMO IN ORBLOC')
      CALL ALLOC(HAMAT ,NISH(IFRP),NISH(IFRP),4,ID='HAMAT IN ORBLOC')
      CALL ALLOC(XKAPPA,NISH(IFRP),NISH(IFRP),4,ID='XKAPPA IN ORBLOC')
      CALL ALLOC(PMAT  ,NISH(IFRP),NISH(IFRP),4,NUCDEP,ID='PMAT IN ORBLOC')

      CALL ALLOC(NFNUC ,2*NUCDEP,ID='NFNUC IN ORBLOC')
      CALL ALLOC(ISEL  ,NISH(IFRP),ID='ISEL IN ORBLOC')
      CALL ALLOC(INDAO ,NFBAS(IFRP,0),ID='INDAO IN ORBLOC')
      CALL ALLOC(INDAOS,NFBAS(IFRP,0),ID='INDAOS IN ORBLOC')
      CALL ALLOC(INDNAO,2*NUCDEP+2,ID='INDNAO IN ORBLOC')

! Read coefficients

      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO_no_work(LUCOEF,'DFCOEF',CMO,IDUM,IDUM,TOTERG,2)
      CLOSE(LUCOEF,STATUS='KEEP')

! Get overlap matrix 

      CALL GTOVLX(SMAT,SSMTRC)
      IF(IPRLOC.GE.6) THEN
        WRITE(LUPRI,'(/A)') '* ORBLOC: Molecular coefficients:'
        CALL PRQMAT(CMO,NFBAS(IFRP,0),NORB(IFRP), &
                    NFBAS(IFRP,0),NORB(IFRP),4,IPQTOQ(1,0),LUPRI)

        WRITE(LUPRI,'(/A)') '* ORBLOC: Overlap matrix' 
        CALL PRQMAT(SMAT,NFBAS(IFRP,0),NFBAS(IFRP,0), &
                    NFBAS(IFRP,0),NFBAS(IFRP,0),1,IPQTOQ(1,0),LUPRI)
      ENDIF

! Occupied electronic orbital coeffients selection

      CALL OCCUPY(ISEL,NISH(IFRP))
      CALL SELCFS(CMO,IFRP,SCOEF,NISH(IFRP),ISEL,D0,NISH(IFRP),NFBAS(IFRP,0),NORB(IFRP))
      IF(IPRLOC.GE.5) THEN
        WRITE(LUPRI,'(/A)') '* ORBLO1: Selected molecular coefficients'
        CALL PRQMAT(SCOEF,NFBAS(IFRP,0),NISH(IFRP), &
                    NFBAS(IFRP,0),NISH(IFRP),4,IPQTOQ(1,0),LUPRI)
      ENDIF

! Get the number of functions per center

      CALL INDAOC(INDAO,INDAOS,INDNAO)
      CALL ORBLO4(INDNAO,NFNUC)

! We do orbital rotation until the difference between two succesive
! localisation criteria is less than a threshold

      TIMSTR = SECOND()
      DINVM = D0
      XNRM  = D0
      IF(PRJLOC) THEN
!       PRJLOC = TRUE localization with projection analysis
        CALL ORBLP1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
      ELSE
!       PRJLOC = FALSE localization with Mulliken analysis
        CALL ORBLO1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),NFNUC,DINVM)
      ENDIF
      TOBE = .FALSE.
      ITR = 0
      FLAG = .TRUE.
      IF (AUTO) THEN
 100    CONTINUE
        ITR = ITR + 1
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)') &
          '########## START ITERATION NO.', ITR , ' ##########',DAYTID
        DINVMO = DINVM
        CALL ORBLO2(PMAT,XKAPPA,HAMAT,XNRM,FLAG,HESLOC,IPRLOC,IPQTOQ,LUPRI)
        IF (.NOT. FLAG) THEN
          CALL DCOPY(NISH(IFRP)*NFBAS(IFRP,0)*4,SCOEF,1,SCOEFO,1)
        ENDIF
        CALL ROTCMO_MOD(SCOEF,NFBAS(IFRP,0),NISH(IFRP),4,        &
                        XKAPPA,NISH(IFRP),IPQTOQ(IFRP,0),IPRLOC)
        IF(PRJLOC) THEN
!         PRJLOC = TRUE localization with projection analysis
          CALL GTOVLX(SMAT,SSMTRC)
          CALL ORBLP1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
        ELSE
!         PRJLOC = FALSE localization with Mulliken analysis
          CALL ORBLO1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),NFNUC,DINVM)
        ENDIF
        DIFF  = DINVM - DINVMO
        DIFFA = ABS(DINVM - DINVMO)
        IF (DIFF .LT. D0 .AND. .NOT. FLAG) THEN
          CALL DCOPY(NISH(IFRP)*NFBAS(IFRP,0)*4,SCOEFO,1,SCOEF,1)
          DINVM = DINVMO
          ITR = ITR - 1
          FLAG = .TRUE.
          WRITE(LUPRI,'(/5X,A)')'! Localization criterion decreases !'
          WRITE(LUPRI,'(5X,A)') '! Restart from  previous iteration !'
          WRITE(LUPRI,'(5X,A)') '!  and activation of the Hessian   !'
          IF(PRJLOC) THEN
            CALL ORBLP1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
          ELSE
            CALL ORBLO1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),NFNUC,DINVM)
          ENDIF
          GOTO 100
        ELSE
          WRITE(6,'(A,I4,D12.4)') 'DIFFA =',ITR,DIFFA
          IF(DIFFA .GT. THRLOC .AND. ITR .LT. ITRLOC) GOTO 100
          IF(DIFFA .LE. THRLOC) THEN
            WRITE(LUPRI,'(/A,I5,A)') ' >>> Convergence in ' , ITR , ' iterations'
          ELSE
            WRITE(LUPRI,'(/A,I5,A)') ' >>> No convergence after ', ITR , ' iterations'
          ENDIF
        ENDIF
      ELSE
 200    CONTINUE
        ITR = ITR + 1
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,I4,A,3X,A24/)') &
          '########## START ITERATION NO.', ITR , ' ##########',DAYTID
        DINVMO = DINVM
        CALL ORBLO2(PMAT,XKAPPA,HAMAT,XNRM,FLAG,HESLOC,IPRLOC,IPQTOQ,LUPRI)
        CALL DCOPY(NISH(IFRP)*NFBAS(IFRP,0)*4,SCOEF,1,SCOEFO,1)
        CALL ROTCMO_MOD(SCOEF,NFBAS(IFRP,0),NISH(IFRP),4, &
                        XKAPPA,NISH(IFRP),IPQTOQ(IFRP,0),IPRLOC)
        IF(PRJLOC) THEN
!         PRJLOC = TRUE localization with projection analysis
          CALL GTOVLX(SMAT,SSMTRC)
          CALL ORBLP1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
        ELSE
!         PRJLOC = FALSE localization with Mulliken analysis
          CALL ORBLO1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),NFNUC,DINVM)
        ENDIF
        DIFF  =     DINVM - DINVMO
        DIFFA = ABS(DINVM - DINVMO)
        IF (DIFF .LT. D0) THEN
          WRITE(LUPRI,'(/5X,A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(LUPRI,'(5X,A)')  '!!! WARNING : DINVM - DINVMO < 0 !!!'
          WRITE(LUPRI,'(5X,A)')  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        ENDIF
        IF (.NOT. FLAG .AND. DIFFA .LE. HESTHR) THEN
          FLAG = .TRUE.
          WRITE(LUPRI,'(/5X,A)')  '!!! Hessian activated !!!'
        ENDIF
        IF(DIFFA .GT. THRLOC .AND. ITR .LT. ITRLOC) GOTO 200
        IF(DIFFA .LE. THRLOC) THEN
          WRITE(LUPRI,'(/A,I5,A)') ' >>> Convergence in ' , ITR , ' iterations'
        ELSE
          WRITE(LUPRI,'(/A,I5,A)') ' >>> No convergence after ', ITR , ' iterations'
        ENDIF
      ENDIF
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      CALL TIMTXT('>>> Time used in localization is :',TIME,LUPRI)

! Put the localised coefficients in DFCOEF

      CALL ORBLO5(SCOEF,NFBAS(IFRP,0),NISH(IFRP),IFRP)

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
      CALL DEALLOC(SCOEF ,ID='SCOEF  IN ORBLOC')
      CALL DEALLOC(SCOEFO,ID='SCOEFO IN ORBLOC')
      CALL DEALLOC(SMAT  ,ID='SMAT   IN ORBLOC')
      CALL DEALLOC(PMAT  ,ID='PMAT   IN ORBLOC')
      CALL DEALLOC(XKAPPA,ID='XKAPPA IN ORBLOC')
      CALL DEALLOC(HAMAT ,ID='HAMAT  IN ORBLOC')
      CALL DEALLOC(CMO   ,ID='CMO    IN ORBLOC')

      CALL DEALLOC(NFNUC ,ID='NFNUC  IN ORBLOC')
      CALL DEALLOC(ISEL  ,ID='ISEL   IN ORBLOC')
      CALL DEALLOC(INDAO ,ID='INDAO  IN ORBLOC')
      CALL DEALLOC(INDAOS,ID='INDAOS IN ORBLOC')
      CALL DEALLOC(INDNAO,ID='INDNAO IN ORBLOC')

      CALL QEXIT('ORBLOC')

  END SUBROUTINE ORBLOC

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLP1(SMAT,SCOEF,PMAT,NBL,NBS,NI,DINVM)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      USE SELECT_COEFF_OWN_BAS
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
      PARAMETER(D0=0.0D0,D2=2.0D0,DM1=-1.0D0,D1=1.0D0)

#include "dcbgen.h"
#include "dcbloc.h"
#include "dcbprj.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"

      LOGICAL TOBE
      DIMENSION NSTR(2,0:2,0:MAXREF),NVECS(2,0:2), &
                SMAT(*),SCOEF(*),PMAT(NI,NI,4,NUCDEP),IOFF(2)
      DIMENSION DDUMARR(2),IDUMARR(2)

      INTEGER*4, ALLOCATABLE, DIMENSION(:) :: KVECTMP, JREF, IPIV
      REAL*8,    ALLOCATABLE, DIMENSION(:) :: RCOEF, CMO, BV, GVEC, AMAT, WMAT

      TYPE :: PTR2ARRAY
        INTEGER*4, DIMENSION(:), POINTER :: P
      END TYPE PTR2ARRAY

      TYPE(PTR2ARRAY), ALLOCATABLE :: KVEC(:,:)

  !------------------------------------------------------------------------------!
      ALLOCATE( KVEC(2,NREFS), STAT=IER )
      IF(IER /= 0)THEN
        WRITE(LUPRI,*) "Error in allocation of KVEC array of pointers in OBLP1"
        WRITE(LUPRI,*) "STAT=",IER
        WRITE(0,*) "Error in allocation ==> see output file "
        STOP
      ENDIF
      DO I=1,2
      DO J=1,NREFS
        NULLIFY(KVEC(I,J)%P)
      ENDDO
      ENDDO

      CALL QENTER('ORBLP1')
      KFRSAV = KFREE
      IFRP   = 1
      NB     = NBL + NBS
!     DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        CALL ALLOC(KVECTMP,NORB(IFRP),ID='KVECTMP IN ORBLP1')
        DO I = 1,NREFS
          NSTR(IFRP,0,I) = 1
          CALL NUMLST(VECREF(IFRP,I),KVECTMP,     &
               NORB(IFRP),-NPSH(IFRP),NESH(IFRP), &
               IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(KVECTMP,NSTR(IFRP,0,I), &
               NPSH(IFRP),NESH(IFRP),         &
               NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NVECS(IFRP,1) = NVECS(IFRP,1) + NSTR(IFRP,0,I)
          CALL ALLOCP(KVEC(IFRP,I)%P,NSTR(IFRP,0,I),ID='KVEC%P IN ORBLP1')
          KVEC(IFRP,I)%P(1:NSTR(IFRP,0,I)) = KVECTMP(1:NSTR(IFRP,0,I))
        ENDDO
        CALL DEALLOC(KVECTMP,ID='KVECTMP IN ORBLP1')
!     ENDDO
      NCDIM = NB*NVECS(IFRP,1)*4
      CALL ALLOC(RCOEF,NCDIM,ID="RCOEF in ORBLP1")

!     COEFFICIENT SELECTION
!     =====================

      ICFS = 6
      IF(IPRPRJ.GE.5) ICFS = 7
      CALL ALLOC(CMO,N2BBASXQ,ID="CMO in ORBLP1")

      IOFF(1) = 0
      IPRPRJ = 4
      DO I = 1,NREFS
        INQUIRE(FILE=REFFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
          WRITE(LUPRI,'(A,A6,A)') 'Reference coefficient file ',REFFIL(I),' not present'
          CALL QUIT('PRJANA: No reference coefficients !')
        ENDIF
        CALL OPNFIL(LUCOEF,REFFIL(I),'OLD','PRJANA')
        IF(OWNBAS) THEN
          CALL SELOWN_DRV(2,CMO,DDUMARR,IDUMARR,           &
               I,1,LUCOEF,REFFIL(I),RCOEF,DDUMARR,IDUMARR, &
               NVECS(1,1),IOFF,KVEC(IFRP,I)%P,             &
               NSTR(1,0,I))
          IF(IPRLOC .GE. 5) THEN
            WRITE(LUPRI,'(A,A6)') '* Selected coefficients for reference ',REFFIL(I)
            NQ = NB*IOFF(IFRP)+1
            CALL PRQMAT(RCOEF(NQ),NB,NSTR(IFRP,0,I),NB,NSTR(IFRP,0,I),4,IPQTOQ(1,0),LUPRI)
          ENDIF
!         DO IFRP = 1, NFSYM
            IOFF(IFRP) = IOFF(IFRP) + NSTR(IFRP,0,I)
!         ENDDO
        ELSE
          WRITE(LUPRI,*) 'Current version work only with .OWNBAS option'
          WRITE(0,*) 'Current version work only with .OWNBAS option'
          STOP
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
      ENDDO

!     Throw away the full coefficients and eigenvalues
!     ================================================

      CALL DEALLOC(CMO,ID='CMO IN ORBLP1')

!     Do projection analysis
!     ======================

!     NBDIM = 0
!     DO IFRP = 1,NFSYM
        NBDIM = NVECS(IFRP,1)*NI*4
!     ENDDO
      CALL ALLOC(BV,NBDIM,ID="BV in ORBLP1")
!     DO IFRP = 1,NFSYM
!     CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
      KBUF   = KFREE
      NMOL   = NI
      NREF   = NVECS(IFRP,1)
      N2REF  = NREF*NREF
      CALL ALLOC(GVEC,(NREF+1)*2,ID="GVEC in ORBLP1")
      CALL ALLOC(AMAT,N2REF*4,ID="AMAT in ORBLP1")
      CALL ALLOC(JREF,NREF*3,ID="JREF in ORBLP1")
      CALL ALLOC(IPIV,NREF,ID="IPIV in ORBLP1")
      IVEC = 0
      DO I = 1,NREFS
        CALL PRJINF(JREF,IFRP,I,IVEC,NSTR,KVEC(IFRP,I)%P)
        IVEC = IVEC + NSTR(IFRP,0,I)
      ENDDO
      CALL ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NI,   &
           NVECS(IFRP,1),BV,GVEC,KVEC(IFRP,1)%P, &
           JREF,NSTR,IPIV)

! Memory allocation for the matrix W = Ac

      CALL ALLOC(WMAT,NREF*NI*4,ID="WMAT in ORBLP1")

! Construction of W = Ac

      IF(IPRLOC .GE. 5) THEN
        WRITE(LUPRI,'(/A)') '* Matrice BV'
        CALL PRQMAT(BV,NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A)') '* Matrice AMAT'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NREF,NI,NREF,D1,                &
           'N','N',IPQTOQ(1,0),AMAT,NREF,NREF,4, &
           'N','N',IPQTOQ(1,0),BV,NREF,NI,4,     &
           D0,IPQTOQ(1,0),WMAT,NREF,NI,4)

      IF(IPRLOC.GE.5) THEN
        WRITE(LUPRI,'(/A)') '* Matrice W'
        CALL PRQMAT(WMAT,NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF

! Construct P_A for each center A

      JOFF = 0
      DO I = 1,NUCDEP
        CALL QGEMM(NI,NI,NSTR(IFRP,0,I),D1,                    &
                   'H','N',IPQTOQ(1,0),BV(JOFF+1),NREF,NI,4,   &
                   'N','N',IPQTOQ(1,0),WMAT(JOFF+1),NREF,NI,4, &
                   D0,IPQTOQ(1,0),PMAT(1,1,1,I),NI,NI,4)
        JOFF = JOFF + NSTR(IFRP,0,I)
      ENDDO

! P_A matrices symmetrization

      DO I = 1,NUCDEP
        CALL FULMAT('S',NI,NI,PMAT(1,1,1,I))
        DO IZ = 2,4
          CALL FULMAT('A',NI,NI,PMAT(1,1,IZ,I))
        ENDDO
      ENDDO

! Multiply pmat with 1/2 just to avoid confusion in the code below
      DO IZ=1,4
        DO I=1,NI
          DO J=1,NI
            DO K=1,NUCDEP
              PMAT(I,J,IZ,K) = PMAT(I,J,IZ,K) / D2
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO IZ=1,4
        WRITE(6,*)
        WRITE(6,'(A,I3)') 'Print Q^A matrices, IZ = ',IZ
        DO I=1, NUCDEP
          WRITE(6,'(A,I4)') ' Nucleus #',I
          DO J=1,NI
            WRITE(6,'(I4,10ES9.1)') J,(PMAT(J,K,IZ,I),K=1,NI)
          ENDDO
          WRITE(6,*)
        ENDDO
      ENDDO

      IF(IPRLOC .GE. 1) THEN
        DO I = 1,NUCDEP
          WRITE(LUPRI,'(/A,I3)') '* ORBLO1: Matrice P  Atom ', I
          WRITE(LUPRI,'(A)') '  ============================'
          CALL PRQMAT(PMAT(1,1,1,I),NI,NI,NI,NI,4,IPQTOQ(1,0),LUPRI)
        ENDDO
      ENDIF

! Calculation and print of the localisation criterion

      WRITE(LUPRI,'(/A)') '** Pipek-Mezey localization criteria **' 
      WRITE(LUPRI,'(A)')  '  ===================================  '
      DINV = D0
      DO I = 1,NI
        QLOC2 = D0
        SUMD  = D0
        DO J = 1,NUCDEP
          QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)
          SUMD  = SUMD + PMAT(I,I,1,J)
        ENDDO
        DLOC = D1/QLOC2
        DINV = DINV + QLOC2
        WRITE(LUPRI,'(A,I3,A,2F7.4)') &
        '* Orbital ', I ,' Pipek-Mezey localization measure d = ',DLOC,SUMD
      ENDDO
      DINVM = DINV/NI
      WRITE(LUPRI,'(/A,2F17.13)') '* Inverse mean delocalization = ',DINVM,1.0D0/DINVM

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
      DO I=1,NREFS
        CALL DEALLOCP(KVEC(1,I)%P,ID='KVEC%P IN ORBLP1')
      ENDDO
      DEALLOCATE( KVEC )
      CALL DEALLOC(RCOEF,ID='RCOEF IN ORBLP1')
      CALL DEALLOC(BV,ID='BV IN ORBLP1')
      CALL DEALLOC(GVEC,ID='GVEC IN ORBLP1')
      CALL DEALLOC(AMAT,ID='AMAT IN ORBLP1')
      CALL DEALLOC(JREF,ID='JREF IN ORBLP1')
      CALL DEALLOC(IPIV,ID='IPIV IN ORBLP1')
      CALL DEALLOC(WMAT,ID='WMAT IN ORBLP1')

      CALL QEXIT('ORBLP1')

  END SUBROUTINE ORBLP1

!================================================================================!
!================================================================================!

 SUBROUTINE ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NMOL,NREF,BVEC,GVEC, &
                   JMOL,JREF,NSTR,IPIVOT)
!***********************************************************************
!
!     Projection analysis for fermion ircop IFRP
!
!     Written by T.Saue March 26 1999
!
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
      PARAMETER(D0=0.0D0,DC=1.0D-2,D2=2.0D0)

#include "dcbloc.h"
#include "dcbprj.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"

      LOGICAL DOPOP
      DIMENSION NSTR(2,0:2,0:MAXREF),JREF(3,*),JMOL(*),SCOEF(*), &
           RCOEF(*),AMAT(NREF,NREF,4),BVEC(NREF,NMOL,4),SMAT(*), &
           GVEC(*),IPIVOT(*)

  !------------------------------------------------------------------------------!
      IFRP = 1
      KFRSAV = KFREE
      DOPOP = NREFS.GT.1
      PTOL = PROTHR*DC

!     Make projection vector

      CALL QTRANS90('AOMO','S',D0,NB,NB,NREF,NMOL, &
                    SMAT,NB,NB,1,IPQTOQ(1,0),      &
                    BVEC,NREF,NMOL,4,IPQTOQ(1,0),  &
                    RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                    SCOEF,NB,NMOL,4,IPQTOQ(1,0),IPRPRJ)
      IF(NREFS.GT.1) THEN

!       Make overlap matrix

        CALL QTRANS90('AOMO','S',D0,NB,NB,NREF,NREF, &
                      SMAT,NB,NB,1,IPQTOQ(1,0),      &
                      AMAT,NREF,NREF,4,IPQTOQ(1,0),  &
                      RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                      RCOEF,NB,NREF,4,IPQTOQ(1,0),IPRPRJ)

!       Print overlaps

        IF(IPRLOC .GE. 5) THEN
           JOFF = 1
           DO J = 1,NREFS
              KOFF = 1
              DO I = 1,(J-1)
                 WRITE(LUPRI,'(A,2(3X,A6))') '*Overlaps :',REFFIL(I),REFFIL(J)
                 CALL PRQMAT(AMAT(KOFF,JOFF,1),                  &
                      NSTR(IFRP,0,I),NSTR(IFRP,0,J),NREF,NREF,4, &
                      IPQTOQ(1,0),LUPRI)
                 KOFF = KOFF + NSTR(IFRP,0,I)
              ENDDO
              JOFF = JOFF + NSTR(IFRP,0,J)
           ENDDO
        ENDIF

!       Solve linear system by Cholesky decomposition
!       with full pivoting
!        --> notice that overlap matrix SMAT is destroyed !

        JOB = 1
        CALL QCHOLD(AMAT,NREF,4,NREF,NREF,SMAT,PTOL,NEFF,JOB,IPIVOT)
        CALL QCHOLS(AMAT,NREF,NEFF,NMOL,4,NREF,NREF,SMAT,BVEC, &
             NREF,NMOL,BVEC,NREF,NMOL,JOB,IPIVOT,SMAT(1+NREF))
        IF(NEFF.LT.NREF) THEN
          WRITE(LUPRI,'(A,A)') '* WARNING: ', &
           ' linear dependencies detected in Cholesky decomposition'
          WRITE(LUPRI,'(11X,A,I5,A,I5)') ' Reduced dimensionality: ',NREF,' --> ',NEFF
        ENDIF

!       Restore overlap matrix

        CALL QHMRST(AMAT,NREF,NZ,NREF,NREF)
      ENDIF

      END SUBROUTINE ORBLP2

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLO1(SMAT,SCOEF,PMAT,NBL,NBS,NI,NFNUC,DINVM)
!***********************************************************************
!
!     Written by S.Dubillard Jan 21 2004
!     Last revision:
!
!     Construction of the matrices P_A=c_A^+Sc from Mulliken analysis
!
!     Calculation of the localisation criteria DINV
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
      PARAMETER(D0 = 0.0D0)
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbloc.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"

      DIMENSION SMAT(*),SCOEF(*),PMAT(NI,NI,4,NUCDEP),NFNUC(*)

      REAL*8, ALLOCATABLE, DIMENSION(:) :: WLMAT, WSMAT

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLO1')

      NB = NBL + NBS

! Memory allocation for the matrix W = Sc

      CALL ALLOC(WLMAT,NBL*NI*4,ID='WLMAT IN ORBLO1')
      CALL ALLOC(WSMAT,NBS*NI*4,ID='WSMAT IN ORBLO1')

! Construction of W = Sc (large and small components)

      CALL QGEMM(NBL,NI,NBL,D1,               &
           'N','N',IPQTOQ(1,0),SMAT,NB,NB,1,  &
           'N','N',IPQTOQ(1,0),SCOEF,NB,NI,4, &
           D0,IPQTOQ(1,0),WLMAT,NBL,NI,4)
      IF(.NOT.TWOCOMPBSS) THEN
         CALL QGEMM(NBS,NI,NBS,D1,                            &
              'N','N',IPQTOQ(1,0),SMAT(1+NB*NBL+NBL),NB,NB,1, &
              'N','N',IPQTOQ(1,0),SCOEF(1+NBL),NB,NI,4,       &
              D0,IPQTOQ(1,0),WSMAT,NBS,NI,4)
      ENDIF

      IF(IPRLOC.GE.5) THEN
         WRITE(LUPRI,'(/A)') '* WL matrix :'
         CALL PRQMAT(WLMAT,NBL,NI,NBL,NI,4,IPQTOQ(1,0),LUPRI)
         IF(.NOT.TWOCOMPBSS) THEN
            WRITE(LUPRI,'(/A)') '* WS matrix :'
            CALL PRQMAT(WSMAT,NBS,NI,NBS,NI,4,IPQTOQ(1,0),LUPRI)
         ENDIF
      ENDIF

! Construct P_A for each center A

      IOFFL = 0
      IOFFS = 0
      DO I = 1,NUCDEP
         CALL QGEMM(NI,NI,NFNUC(I),D1,                     &
              'H','N',IPQTOQ(1,0),SCOEF(1+IOFFL),NB,NI,4,  &
              'N','N',IPQTOQ(1,0),WLMAT(1+IOFFL),NBL,NI,4, &
              D0,IPQTOQ(1,0),PMAT(1,1,1,I),NI,NI,4)
         IF(.NOT.TWOCOMPBSS) THEN
            CALL QGEMM(NI,NI,NFNUC(I+NUCDEP),D1,                 &
                 'H','N',IPQTOQ(1,0),SCOEF(1+NBL+IOFFS),NB,NI,4, &
                 'N','N',IPQTOQ(1,0),WSMAT(1+IOFFS),NBS,NI,4,    &
                 D1,IPQTOQ(1,0),PMAT(1,1,1,I),NI,NI,4)
            IOFFS = IOFFS + NFNUC(NUCDEP+I)
         ENDIF
         IOFFL = IOFFL + NFNUC(I)
      ENDDO

! P_A matrices symetrisation

      DO I = 1,NUCDEP
         CALL FULMAT('S',NI,NI,PMAT(1,1,1,I))
         DO IZ = 2,4
            CALL FULMAT('A',NI,NI,PMAT(1,1,IZ,I))
         ENDDO
      ENDDO
!
      IF(IPRLOC.GE.1) THEN
         DO I = 1,NUCDEP
            WRITE(LUPRI,'(/A,I3)') '* ORBLO1: Matrice P  Atom ', I
            WRITE(LUPRI,'(A)') '  ============================'
            CALL PRQMAT(PMAT(1,1,1,I),NI,NI,NI,NI,4,IPQTOQ(1,0),LUPRI)
         ENDDO
      ENDIF

! Calculation and print of the localisation criterion

      WRITE(LUPRI,'(/A)') '** Pipek-Mezey localization criteria **' 
      WRITE(LUPRI,'(A)')  '  ===================================  '
      DINV = D0
      DO I = 1,NI
         QLOC2 = D0
         DO J = 1,NUCDEP
            QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)/4.0D0
         ENDDO
         DLOC = D1/QLOC2
         DINV = DINV + QLOC2
         WRITE(LUPRI,'(A,I3,A,F7.4)') &
         '* Orbital ', I ,' Pipek-Mezey localization measure d = ', DLOC
      ENDDO
      DINVM = DINV/NI
      WRITE(LUPRI,'(/A,F15.13)') '* Inverse mean delocalization = ',DINVM

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
      CALL DEALLOC(WLMAT,ID='WLMAT IN ORBLO1')
      CALL DEALLOC(WSMAT,ID='WSMAT IN ORBLO1')

      CALL QEXIT('ORBLO1')

  END SUBROUTINE ORBLO1

!================================================================================!
!================================================================================!

  subroutine orblo2(pmat,xkappa,hamat,xnrm,flag,hesloc,iprloc,ipqtoq,lupri)
  !******************************************************************************!
  !                                                                              !
  !     Written by S.Dubillard Jan 21 2004                                       !
  !     Last revision:                                                           !
  !                                                                              !
  !     Construction of the XKAPPA matrix                                        !
  !                                                                              !
  !******************************************************************************!
  IMPLICIT NONE

  ! external variables
  logical flag, hesloc

  integer*4 :: iprloc, lupri
  integer*4, dimension(:,:) :: ipqtoq

  ! internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

  real*8 :: xnrm
  real*8, dimension(:,:,:)   ::  xkappa, hamat
  real*8, dimension(:,:,:,:) ::  pmat

  integer*4 :: i, j, k, iz, nucdep, ni, ndim
  real*8    :: pii, pjj, pn2, dnrm2, gradient

  !------------------------------------------------------------------------------!
  call qenter('orblo2')

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)

  ! Construct unbared diagonal Hessian matrix
  ndim = ni*ni*4
  call dzero(hamat,ndim)
  do j = 1,ni
    do i = 1,ni
      do k = 1,nucdep
        pii = pmat(i,i,1,k)
        pjj = pmat(j,j,1,k)
        pn2 = pmat(j,i,1,k)*pmat(j,i,1,k) + pmat(j,i,2,k)*pmat(j,i,2,k)
        if(i /= j)hamat(i,j,1) = hamat(i,j,1) + d4*pn2 - d2*((pjj-pii)*(pjj-pii))
      enddo
    enddo
  enddo

  ! Construction and print of the unbared kappa matrix
  ndim = ni*ni*4
  call dzero(xkappa,ndim)
  do iz = 1,2
    do j = 1,ni
      do i = 1,ni
        gradient = d0
        do k = 1,nucdep
          gradient = gradient + d2 * ( pmat(i,j,iz,k)*(pmat(i,i,1,k) - pmat(j,j,1,k)) )
        enddo
!ST     if(flag .and. (.not.hesloc)) xkappa(i,j,iz) = - gradient / hamat(i,j,1)
        if(i /= j)xkappa(i,j,iz) = - gradient / hamat(i,j,1)
        if(i == j)xkappa(i,j,iz) = d0
      enddo
    enddo
  enddo

  ! Construct bared diagonal Hessian matrix
  ndim = ni*ni*4
  call dzero(hamat,ndim)
  do j = 1,ni
    do i = 1,ni
      do k = 1,nucdep
        pii = pmat(i,i,1,k)
        pjj = pmat(j,j,1,k)
        pn2 = pmat(j,i,3,k)*pmat(j,i,3,k) + pmat(j,i,4,k)*pmat(j,i,4,k)
        if(i /= j)hamat(i,j,1) = hamat(i,j,1) + d4*pn2 - d2*((pjj-pii)*(pjj-pii))
      enddo
    enddo
  enddo

  ! Construction and print of the bared kappa matrix
  ndim = ni*ni*4
  do iz = 3,4
    do j = 1,ni
      do i = 1,ni
        gradient = d0
        do k = 1,nucdep
          gradient = gradient + d2 * ( pmat(i,j,iz,k)*(pmat(i,i,1,k) - pmat(j,j,1,k)) )
        enddo
!ST     if(flag .and. (.not.hesloc)) xkappa(i,j,iz) = - gradient / hamat(i,j,1)
        if(i /= j)xkappa(i,j,iz) = - gradient / hamat(i,j,1)
        if(i == j)xkappa(i,j,iz) = d0
      enddo
    enddo
  enddo

  do iz=1,4
    write(6,*)
    write(6,'(a,i4)') 'Print xkappa matrix, IZ =',iz
    do j=1,ni
      write(6,'(i4,10es9.1)') j,(xkappa(j,k,1),k=1,ni)
    enddo
  enddo

  ! Print some kappa matrix
  if(iprloc .ge. 5) then
     write (lupri,'(/a)') '* KAPPA matrix :'
     call prqmat(xkappa(1,1,1),ni,ni,ni,ni,4,ipqtoq(1,0),lupri)
  endif

  ! Norm of the gradient
  xnrm = dnrm2(ndim,xkappa,1)
  write (lupri,'(a,e12.2)') 'Norm of gradient: ',dsqrt(xnrm)

  call qexit('orblo2')

  end subroutine orblo2

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLO4(INDNAO,NFNUC)
!***********************************************************************
!
!     Written by S.Dubillard Jan 21 2004
!     Last revision:
!
!     Construction of the array NFNUC which contains the number of
!     large components basis functions for each center and then the 
!     number of small components basis functions for each center
!
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbloc.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"

      DIMENSION INDNAO(2,*)
      DIMENSION NFNUC(2*NUCDEP)

  !------------------------------------------------------------------------------!
      IFRP = 1
      CALL QENTER('ORBLO4')
      DO I = 1,NUCDEP-1
         NFNUC(I)=INDNAO(1,I+1)-INDNAO(1,I)
      ENDDO
      NFNUC(NUCDEP)=INDNAO(2,1)-INDNAO(1,NUCDEP)
      DO I = 1,NUCDEP-1
         NFNUC(I+NUCDEP)=INDNAO(2,I+1)-INDNAO(2,I)
      ENDDO
      NFNUC(2*NUCDEP)=NFBAS(IFRP,0)-INDNAO(2,NUCDEP)
      CALL QEXIT('ORBLO4')

  END SUBROUTINE ORBLO4

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLO5(SCOEF,NB,NI,IFRP)
!***********************************************************************
!
!     Written by S.Dubillard Apr 13 2004
!     Last revision:
!
!     Put the localised coefficients into the DFCOEF file (1/2)
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
      PARAMETER(D0 = 0.0D0)
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbloc.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"

      DIMENSION SCOEF(NB,NI,4)

      INTEGER*4, ALLOCATABLE :: IBEIG(:)
      REAL*8,    ALLOCATABLE :: CMO(:,:,:), EIG(:)

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLO5')

      IFRP = 1
      NO = NORB(IFRP)

      ! allocate arrays
      CALL ALLOC(CMO  ,NB,NO,4,ID='CMO   IN ORBLO5')
      CALL ALLOC(EIG  ,NORBT  ,ID='EIG   IN ORBLO5')
      CALL ALLOC(IBEIG,NORBT  ,ID='IBEIG IN ORBLO5')

      ! open DFCOEF file and read EIG, IBEIG and TOTERG
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO_no_work(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)

      ! reorganize MO
      DO I = 1,NB
        DO J = 1,NI
          DO K = 1,4
            L = J + NO/2
            IF (LEVYLE) L = J
            CMO(I,L,K) = SCOEF(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      IF (IPRLOC .GE. 5) THEN
         WRITE(LUPRI,'(/A)') '*ORBLO6 : New CMO matrix'
         CALL PRQMAT(CMO,NB,NO,NB,NO,4,IPQTOQ(1,0),LUPRI)
      ENDIF

      ! write MO to file and close him
      CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
      CLOSE(LUCOEF,STATUS='KEEP')

      ! deallocate arrays
      CALL DEALLOC(CMO  ,ID='CMO   IN ORBLO5')
      CALL DEALLOC(EIG  ,ID='EIG   IN ORBLO5')
      CALL DEALLOC(IBEIG,ID='IBEIG IN ORBLO5')

      CALL QEXIT('ORBLO5')

  END SUBROUTINE ORBLO5

!================================================================================!
!================================================================================!

  SUBROUTINE ROTCMO_MOD(CMO,NBASI,NORBI,NZ,XKAPPA,NORBT,IP,IPRINT)
!***********************************************************************
!
!     Rotate coefficients
!
!     CMO_{new} = CMO_{old} exp( -\kappa )
!     We use a trick a la Sostrup notes, chapter 3.1.5
!
!     Input:
!        CMO     - old coefficients
!        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
!        IPRINT  - print level
!
!     Output:
!        CMO     - new coefficients
!
!     Written by J. Thyssen - Nov 24 1998
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
#include "priunit.h"
#include "thrzer.h"

      PARAMETER (D1 = 1.00D00, D0 = 0.00D00, DM1 = -D1)
      PARAMETER (DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00)
      PARAMETER (D1P5 = 1.5D00, D2 = 2.0D00, THREQL = 1.0D-12)

      DIMENSION CMO(NBASI,NORBI,NZ), XKAPPA(NORBT,NORBT,NZ), IP(4)
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12

      REAL*8, ALLOCATABLE, DIMENSION(:) :: REIG, RSIN, RCOS, RV, RTMP

  !------------------------------------------------------------------------------!
      CALL QENTER('ROTCMO_MOD')

      CALL GETTIM(CPU1,WALL1)
      IF (IPRINT .GE. 5) CALL HEADER('Output from ROTCMO_MOD',-1)

      NORBXQ = NORBI*NORBI*NZ
      CALL ALLOC(REIG,NORBI ,ID='REIG IN ROTCMO_MOD')
      CALL ALLOC(RSIN,NORBXQ,ID='RSIN IN ROTCMO_MOD')
      CALL ALLOC(RV  ,NORBXQ,ID='KV   IN ROTCMO_MOD')
      CALL ALLOC(RTMP,NORBXQ,ID='RTMP IN ROTCMO_MOD')

      CALL ALLOC(RCOS,NORBI*NBASI*NZ,ID='RCOS IN ROTCMO_MOD')

!     Calculate -\kappa^2 = \kappa^{\dagger} \kappa
!     (use RCOS as temporary array)
!     ---------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,D1,             &
                 'H','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D0,IP,RCOS,NORBI,NORBI,NZ)

      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: -kappa^2',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Diagonalize -\kappa^2
!     ---------------------------------------------

      CALL QDIAG90(NZ,NORBI,RCOS,NORBI,NORBI, &
                   REIG,1,RV,NORBI,NORBI,IERR)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: V matrix (eigenvector)',-1)
         CALL PRQMAT(RV,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO_MOD: eigenvalues',-1)
         CALL PRQMAT(REIG,1,NORBI,1,NORBI,1,IP,LUPRI)
      END IF

!     Make cos(d) and sin(d)/d ; d = sqrt(eigenvalues)
!     ------------------------------------------------

      CALL GETTIM(CPU2,WALL2)
      RCOS(1:NORBXQ) = 0.0D0
      RSIN(1:NORBXQ) = 0.0D0
      DO I = 1,NORBI
         IADDR = (I-1) + (I-1)*NORBI
         EIG = REIG(I)
         IF (EIG .LT. D0) THEN
            IF (IPRINT .GE. 3 .AND. ABS(EIG).GT. THRZER) WRITE(LUPRI,9010) EIG
            EIG = D0
            REIG(I)   = D0
         END IF
         IF (EIG .LT. DTHRS) THEN

!           Use Taylor expansion to avoid division by zero and/or
!           to speed up evaluation (avoid call to sin/cos/sqrt).
!           For EIG < 1.0D-8 the 3. order Taylor expansion is correct
!             to machine precision.

            RCOS(1+IADDR) = D1-DP5*EIG
            RSIN(1+IADDR) = D1-DSIXTH*EIG
         ELSE
            SQEIG = SQRT(EIG)
            RCOS(1+IADDR) = COS(SQEIG)
            RSIN(1+IADDR) = SIN(SQEIG)/SQEIG
         END IF
      END DO
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: sin(d)/d',-1)
         CALL PRQMAT(RSIN,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO_MOD: cos(d)',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU3,WALL3)
 9010 FORMAT('WARNING from ROTCMO_MOD: Negative eigenvalue ',1P,D15.6)

!     Make V cos(d) V^{\dagger}
!     (use KTMP as temporatory storage)
!     ------------------------------------

      CALL QTRANS90('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI, &
                    RTMP,NORBI,NORBI,NZ,IP,                &
                    RCOS,NORBI,NORBI,NZ,IP,                &
                    RV,NORBI,NORBI,NZ,IP,                  &
                    RV,NORBI,NORBI,NZ,IP,IPRINT)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: V cos(d) V^{dagger}',-1)
         CALL PRQMAT(RTMP,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V sin(d)/d V^{\dagger}
!     (use RCOS as temporatory storage)
!     ------------------------------------

      CALL QTRANS90('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI, &
                    RCOS,NORBI,NORBI,NZ,IP,                &
                    RSIN,NORBI,NORBI,NZ,IP,                &
                    RV,NORBI,NORBI,NZ,IP,                  &
                    RV,NORBI,NORBI,NZ,IP,IPRINT)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
!     ---------------------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,DM1,            &
                 'N','N',IP,RCOS,NORBI,NORBI,NZ,   &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D1,IP,RTMP,NORBI,NORBI,NZ)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: exp(-kappa)',-1)
         CALL PRQMAT(RTMP,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Check if unitary

      CALL QGEMM(NORBI,NORBI,NORBI,D1,           &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ, &
                 'H','N',IP,RTMP,NORBI,NORBI,NZ, &
                 D0,IP,RCOS,NORBI,NORBI,NZ)
      IOFF = 0
      I_MAX = 0
      J_MAX = 0
      D_MAX = d0
      DO I = 1,NORBI
      DO J = 1,NORBI
         IF (I.EQ.J) THEN
            IF (ABS(RCOS(1+IOFF)-D1).GT.D_MAX) THEN
               I_MAX = I
               J_MAX = J
               D_MAX = ABS(RCOS(1+IOFF)-D1)
            END IF
         ELSE
            IF (RCOS(1+IOFF).GT.D_MAX) THEN
               I_MAX = I
               J_MAX = J
               D_MAX = ABS(RCOS(1+IOFF))
            END IF
         END IF
         IOFF = IOFF + 1
      END DO
      END DO
      IF (D_MAX .GT. THRZER) THEN
         WRITE(lupri,*) 'ROTCMO_MOD: maximum deviation from unitary matrix'
         WRITE(lupri,*) 'Element ',J_MAX,I_MAX,' value :',D_MAX
      END IF

!     Calculate new orbitals: CMO_{old} exp(-kappa)
!     ---------------------------------------------

      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
!SK   CALL DZERO(RCOS,NORBXQ)
      CALL QGEMM(NBASI,NORBI,NORBI,D1,           &
                 'N','N',IP,CMO,NBASI,NORBI,NZ,  &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ, &
                 D0,IP,RCOS,NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,RCOS,1,CMO,1)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('ROTCMO_MOD: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (IPRINT .GE. 2) THEN
         CPUTID = SECTID(CPU4-CPU1)
         WALLTID = SECTID(WALL4-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID
      END IF
 9020 FORMAT(/'ROTCMO_MOD: total CPU (WALL) time       : ',A12,' (',A12,')')
 9021 FORMAT( 'ROTCMO_MOD: sin/sqrt/cos CPU (WALL) time: ',A12,' (',A12,')')

      CALL DEALLOC(REIG,ID='REIG IN ROTCMO_MOD')
      CALL DEALLOC(RSIN,ID='RSIN IN ROTCMO_MOD')
      CALL DEALLOC(RCOS,ID='RCOS IN ROTCMO_MOD')
      CALL DEALLOC(RV  ,ID='KV   IN ROTCMO_MOD')
      CALL DEALLOC(RTMP,ID='RTMP IN ROTCMO_MOD')

      CALL QEXIT('ROTCMO_MOD')

  END SUBROUTINE ROTCMO_MOD

!================================================================================!
!                             END OF THE MODULE                                  !
!================================================================================!
end module localization
