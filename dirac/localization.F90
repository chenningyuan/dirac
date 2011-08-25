module localization
implicit none

public :: ORBLP1

contains

!================================================================================!
!================================================================================!

  SUBROUTINE ORBLP1(SMAT,SCOEF,PMAT,NBL,NBS,NI,DINVM,WORK,LWORK)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
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
      DIMENSION WORK(LWORK)
      DIMENSION NSTR(2,0:2,0:MAXREF),KVEC(2,0:MAXREF),NVECS(2,0:2), &
                KBV(2),SMAT(*),SCOEF(*),PMAT(NI,NI,4,NUCDEP),IOFF(2)
      DIMENSION DDUMARR(2),IDUMARR(2)
      WRITE(6,*) "I'm in localization module"

!     INTEGER*4, ALLOCATABLE :: KVEC(:,:,:)
!
!     TYPE :: KVEC
!       REAL, DIMENSION(:), POINTER :: KVEC
!     END TYPE KVEC
!     TYPE(ptr_to_array), ALLOCATABLE :: x(:)
!     !...
!     ALLOCATE(x(nx))
!     DO i = 1,nx
!        ALLOCATE(x(i)%arr(m))
!     END DO

#include "memint.h"
      CALL QENTER('ORBLP1')
      KFRSAV = KFREE
      IFRP   = 1
      NB     = NBL + NBS
!     DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        DO I = 1,NREFS
          CALL MEMGET('INTE',KVEC(IFRP,I),NORB(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL NUMLST(VECREF(IFRP,I),WORK(KVEC(IFRP,I)), &
               NORB(IFRP),-NPSH(IFRP),NESH(IFRP),        &
               IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I), &
               NPSH(IFRP),NESH(IFRP),                    &
               NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NVECS(IFRP,1) = NVECS(IFRP,1) + NSTR(IFRP,0,I)
          CALL MEMREL('PRJANA',WORK,KVEC(IFRP,I),KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),WORK,KFREE,LFREE)
        ENDDO
!     ENDDO
      NCDIM = NB*NVECS(IFRP,1)*4
      CALL MEMGET('REAL',KRCOEF,NCDIM,WORK,KFREE,LFREE)

!     COEFFICIENT SELECTION
!     =====================

      ICFS = 6
      IF(IPRPRJ.GE.5) ICFS = 7
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)

      IOFF(1) = 0
      IPRPRJ = 4
      NQ = KRCOEF
      DO I = 1,NREFS
        INQUIRE(FILE=REFFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
          WRITE(LUPRI,'(A,A6,A)') 'Reference coefficient file ',REFFIL(I),' not present'
          CALL QUIT('PRJANA: No reference coefficients !')
        ENDIF
        CALL OPNFIL(LUCOEF,REFFIL(I),'OLD','PRJANA')
        IF(OWNBAS) THEN
          CALL SELOWN_DRV(2,WORK(KCMO),DDUMARR,IDUMARR,          &
               I,1,LUCOEF,REFFIL(I),WORK(KRCOEF),DDUMARR,IDUMARR,&
               NVECS(1,1),IOFF,VECREF(IFRP,I),IFRP,              &
               NSTR(1,0,I),WORK,KFREE,LFREE,                     &
               -NPSH(IFRP),NESH(IFRP))
          IF(IPRLOC .GE. 5) THEN
            WRITE(LUPRI,'(A,A6)') '* Selected coefficients for reference ',REFFIL(I)
            NQ = NB*IOFF(IFRP)+1
            CALL PRQMAT(WORK(NQ),NB,NSTR(IFRP,0,I),NB,NSTR(IFRP,0,I),4,IPQTOQ(1,0),LUPRI)
          ENDIF
!         DO IFRP = 1, NFSYM
            IOFF(IFRP) = IOFF(IFRP) + NSTR(IFRP,0,I)
            IF(NSTR(IFRP,0,I).GT.0) THEN
               NQ = NQ + NSTR(IFRP,0,I)*NB
            ENDIF
!         ENDDO
        ELSE
          CALL REACMO(LUCOEF,REFFIL(I),WORK(KCMO),IDUM,IDUM,TOTERG,2,WORK,KFREE,LFREE)
!         DO IFRP = 1, NFSYM
            IF(NSTR(IFRP,0,I).GT.0) THEN
              CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(NQ), &
                   NVECS(IFRP,1),WORK(KVEC(IFRP,I)),             &
                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),NB,NORB(IFRP))
              IF(IPRLOC .GE. 5) THEN
                WRITE(LUPRI,'(A,A6)') '* Selected coefficients for reference ',REFFIL(I)
                CALL PRQMAT(WORK(NQ),NB,NSTR(IFRP,0,I),NB,NSTR(IFRP,0,I),4,IPQTOQ(1,0),LUPRI)
              ENDIF
              NQ = NQ + NSTR(IFRP,0,I)*NB
            ENDIF
!         ENDDO
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
      ENDDO

!     Throw away the full coefficients and eigenvalues
!     ================================================

      CALL MEMREL('PRJANA.cfs',WORK,KCMO,KCMO,KFREE,LFREE)

!     Do projection analysis
!     ======================

!     NBDIM = 0
!     DO IFRP = 1,NFSYM
        NBDIM = NVECS(IFRP,1)*NI*4
!     ENDDO
      CALL MEMGET('REAL',KBV(1),NBDIM,WORK,KFREE,LFREE)
!     IF(NFSYM.EQ.2) THEN
!       KBV(2)   = KBV(1) + NVECS(1,0)*NVECS(1,2)*NZ
!     ELSE
        KBV(2)   = KFREE
!     ENDIF
!     DO IFRP = 1,NFSYM
!     CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
      KBUF   = KFREE
      NMOL   = NI
      NREF   = NVECS(IFRP,1)
      N2REF  = NREF*NREF
      CALL MEMGET('REAL',KGVEC,(NREF+1)*2,WORK,KFREE,LFREE)
      IF(NREFS.GT.1) THEN
        CALL MEMGET('REAL',KAMAT,N2REF*4,WORK,KFREE,LFREE)
      ELSE
        KAMAT = KFREE
      ENDIF
      CALL MEMGET('INTE',KJREF,NREF*3,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPIV,NREF  ,WORK,KFREE,LFREE)
      IVEC = 0
      DO I = 1,NREFS
        CALL PRJINF(WORK(KJREF),IFRP,I,IVEC,NSTR,WORK(KVEC(IFRP,I)))
        IVEC = IVEC + NSTR(IFRP,0,I)
      ENDDO
      CALL ORBLP2(SMAT,SCOEF,WORK(KRCOEF),WORK(KAMAT),NB,NI,             &
           NVECS(IFRP,1),WORK(KBV(IFRP)),WORK(KGVEC),WORK(KVEC(IFRP,1)), &
           WORK(KJREF),NSTR,WORK(KIPIV),WORK,KFREE,LFREE)

      CALL MEMREL('PRJANA',WORK,KFRSAV,KFRSAV,KFREE,LFREE)

! Memory allocation for the matrix W = Ac

      CALL MEMGET('REAL',KWMAT,NREF*NI*4,WORK,KFREE,LFREE)

! Construction of W = Ac

      IF(IPRLOC .GE. 5) THEN
        WRITE(LUPRI,'(/A)') '* Matrice KBV'
        CALL PRQMAT(WORK(KBV(IFRP)),NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A)') '* Matrice AMAT'
        CALL PRQMAT(WORK(KAMAT),NREF,NREF,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF
        CALL QGEMM(NREF,NI,NREF,D1,                         &
             'N','N',IPQTOQ(1,0),WORK(KAMAT),NREF,NREF,4,   &
             'N','N',IPQTOQ(1,0),WORK(KBV(IFRP)),NREF,NI,4, &
             D0,IPQTOQ(1,0),WORK(KWMAT),NREF,NI,4)

      IF(IPRLOC.GE.5) THEN
        WRITE(LUPRI,'(/A)') '* Matrice W'
        CALL PRQMAT(WORK(KWMAT),NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF

! Construct P_A for each center A

      JOFF = 0
      DO I = 1,NUCDEP
        CALL QGEMM(NI,NI,NSTR(IFRP,0,I),D1,                            &
                   'H','N',IPQTOQ(1,0),WORK(KBV(IFRP)+JOFF),NREF,NI,4, &
                   'N','N',IPQTOQ(1,0),WORK(KWMAT+JOFF),NREF,NI,4,     &
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

      CALL MEMREL('ORBLO1',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('ORBLP1')

  END SUBROUTINE ORBLP1

!================================================================================!
!                             END OF THE MODULE                                  !
!================================================================================!
end module localization
