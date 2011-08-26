module localization
implicit none

public :: ORBLOC

private :: ORBLP1, ORBLP2

contains

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLOC(WORK,LWORK)
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
      DIMENSION WORK(LWORK)
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12,DAYTID*24

      INTEGER*4, ALLOCATABLE, DIMENSION(:) :: NFNUC, ISEL, INDAO, INDAOS, INDNAO
      REAL*8,    ALLOCATABLE, DIMENSION(:) :: SCOEF, SCOEFO, SMAT, PMAT, XKAPPA, HAMAT, CMO

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLOC')
#include "memint.h"
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
      CALL ALLOC(PMAT  ,NUCDEP*NISH(IFRP)*NISH(IFRP)*4,ID='PMAT IN ORBLOC')
      CALL ALLOC(XKAPPA,NISH(IFRP)*NISH(IFRP)*4,ID='XKAPPA IN ORBLOC')
      CALL ALLOC(HAMAT ,NISH(IFRP)*NISH(IFRP)*4,ID='HAMAT IN ORBLOC')
      CALL ALLOC(CMO   ,N2BBASXQ,ID='CMO IN ORBLOC')

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
      CALL SELCFS(CMO,IFRP,SCOEF,NISH(IFRP),ISEL, &
                  D0,NISH(IFRP),NFBAS(IFRP,0),NORB(IFRP))
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
        CALL ORBLP1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1), &
                    NFBAS(IFRP,2),NISH(IFRP),DINVM)
      ELSE
!       PRJLOC = FALSE localization with Mulliken analysis
        CALL ORBLO1(SMAT,SCOEF,PMAT,NFBAS(IFRP,1),        &
                    NFBAS(IFRP,2),NISH(IFRP),NFNUC,DINVM, &
                    WORK(KFREE),LFREE)
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
        CALL ORBLO2(PMAT,NISH(IFRP),XKAPPA,HAMAT,XNRM,FLAG)
        IF (.NOT. FLAG) THEN
          CALL DCOPY(NISH(IFRP)*NFBAS(IFRP,0)*4,SCOEF,1,SCOEFO,1)
        ENDIF
        CALL ROTCMO(SCOEF,NFBAS(IFRP,0),NISH(IFRP),4,        &
                    XKAPPA,NISH(IFRP),IPQTOQ(IFRP,0),IPRLOC, &
                    WORK(KFREE),LFREE)
        IF(PRJLOC) THEN
!         PRJLOC = TRUE localization with projection analysis
          CALL GTOVLX(SMAT,SSMTRC)
          CALL ORBLP1(SMAT,SCOEF,PMAT, &
                      NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
        ELSE
!         PRJLOC = FALSE localization with Mulliken analysis
          CALL ORBLO1(SMAT,SCOEF,PMAT,                        &
                      NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP), &
                      NFNUC,DINVM,WORK(KFREE),LFREE)
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
            CALL ORBLP1(SMAT,SCOEF,PMAT, &
                        NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
          ELSE
            CALL ORBLO1(SMAT,SCOEF,PMAT,                        &
                        NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP), &
                        NFNUC,DINVM,WORK(KFREE),LFREE)
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
        CALL ORBLO2(PMAT,NISH(IFRP),XKAPPA,HAMAT,XNRM,FLAG)
        CALL DCOPY(NISH(IFRP)*NFBAS(IFRP,0)*4,SCOEF,1,SCOEFO,1)
        CALL ROTCMO(SCOEF,NFBAS(IFRP,0),NISH(IFRP),4,        &
                    XKAPPA,NISH(IFRP),IPQTOQ(IFRP,0),IPRLOC, &
                    WORK(KFREE),LFREE)
        IF(PRJLOC) THEN
!         PRJLOC = TRUE localization with projection analysis
          CALL GTOVLX(SMAT,SSMTRC)
          CALL ORBLP1(SMAT,SCOEF,PMAT, &
                      NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP),DINVM)
        ELSE
!         PRJLOC = FALSE localization with Mulliken analysis
          CALL ORBLO1(SMAT,SCOEF,PMAT,   &
                      NFBAS(IFRP,1),NFBAS(IFRP,2),NISH(IFRP), &
                      NFNUC,DINVM,WORK(KFREE),LFREE)
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

      CALL ORBLO5(SCOEF,NFBAS(IFRP,0),NISH(IFRP),IFRP,WORK(KFREE),LFREE)

! Memory deallocation

      CALL MEMREL('ORBLOC',WORK,KWORK,KWORK,KFREE,LFREE)

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
      WRITE(6,*) "I'm in localization module"

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

      WRITE(6,*)
      WRITE(6,*) 'Print Q^A matrices'
      DO I=1, NUCDEP
        WRITE(6,'(A,I4)') ' Nucleus #',I
        DO J=1,NI
          WRITE(6,'(I4,10ES9.1)') J,(PMAT(J,K,1,I)/2.0D0,K=1,NI)
        ENDDO
        WRITE(6,*)
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
          QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)/4.0D0
          SUMD  = SUMD + PMAT(I,I,1,J)/2.0D0
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
!                             END OF THE MODULE                                  !
!================================================================================!
end module localization
