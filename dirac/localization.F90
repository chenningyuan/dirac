module localization
implicit none

public :: ORBLOC, initialize_localization

private :: ORBLP1, ORBLP2, ORBLO1, ORBLO4, ORBLO5, ROTCMO_MOD
private :: getk_cdiag, getk_rdiag, getk_full, get_hessian, get_gradient

type :: localization_flags
  integer*4 :: iprloc, itrloc
  real*8    :: thrloc, hesthr
  logical   :: prjloc, hesloc, auto
end type localization_flags

type :: external_variales
  integer*4 :: lupri, lucoef
  integer*4 :: nfbas(2,0:2), n2bbasx, n2bbasxq
  integer*4 :: nish(2), norb(2), npsh(2), nesh(2), norbt
  real*8    :: ssmtrc
  logical   :: twocompbss, levyle
  integer*4 :: nz, ipqtoq(4,0:7)
  integer*4 :: nucdep
  integer*4 :: thrzer
  integer*4 :: nrefs, iprprj
  integer*4 :: maxref
  real*8    :: prothr
  logical   :: ownbas
  character*6  :: reffil(15)    ! reffil(maxref)    \
  character*72 :: vecref(2,15)  ! vecref(2,maxref)  /  see dcbprj.h
end type external_variales

type( localization_flags ), private :: loc
type( external_variales  ), private :: var

contains

  subroutine initialize_localization
  !==============================================================================!
  !               Initialize localization for Dirac program                      !
  !==============================================================================!
  implicit none

#include "dcbloc.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "mxcent.h"
#include "nuclei.h"
#include "thrzer.h"
#include "dcbprj.h"

  ! Internal variables
  integer*4 :: i

  !------------------------------------------------------------------------------!
  ! dcbloc.h
  loc%iprloc     = iprloc
  loc%itrloc     = itrloc
  loc%thrloc     = thrloc
  loc%hesthr     = hesthr
  loc%prjloc     = prjloc
  loc%hesloc     = hesloc
  loc%auto       = auto
  ! priunit.h
  var%lupri      = lupri
  ! dcbgen.h
  var%lucoef     = lucoef
  ! dcbbas.h
  var%nfbas(:,:) = nfbas(:,:)
  var%n2bbasx    = n2bbasx
  var%n2bbasxq   = n2bbasxq
  ! dcborb.h
  var%nish(:)    = nish(:)
  var%norb(:)    = norb(:)
  var%npsh(:)    = npsh(:)
  var%nesh(:)    = nesh(:)
  var%norbt      = norbt
  ! dcbham.h
  var%ssmtrc     = ssmtrc
  var%twocompbss = twocompbss
  var%levyle     = levyle
  ! dgroup.h
  var%ipqtoq(:,:)= ipqtoq(:,:)
  var%nz         = nz
  ! nuclei.h
  var%nucdep     = nucdep
  ! thrzer.h
  var%thrzer     = thrzer
  ! dcbprj.h
  var%nrefs      = nrefs
  var%iprprj     = iprprj
  var%maxref     = min(maxref,15)
  var%prothr     = prothr
  var%ownbas     = ownbas

  if(maxref /= 15)then
    write(lupri,'(a)')' WARNING: The maxref dimension in localization module and Dirac differ'
    write(lupri,'(a)')'          See type external_variales in localization.f90 '
    write(lupri,'(a,I4)')' maxref(Dirac) =',maxref
    write(lupri,'(a,I4)')' maxref(Localization module) = 15'
  endif

  do i=1,var%maxref
    var%reffil(i)   = reffil(i)
    var%vecref(1,i) = vecref(1,i)
    var%vecref(2,i) = vecref(2,i)
  enddo

  end subroutine initialize_localization

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine orbloc
  !==============================================================================!
  !                                                                              !
  !     Molecular orbital localization                                           !
  !                                                                              !
  !==============================================================================!
  use memory_allocator
  implicit none

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0

  character*12 :: sectid, cputid, walltid
  character    :: daytid*24, hessian*4 = 'DIAG', conv*20 = 'No convergence after'
  logical      :: tobe, full, diagonly

  integer*4 :: ifrp, nish12, idum, itr, lupri, nish, norb, nucdep, maxref
  integer*4 :: nfbas(0:2), ipqtoq(4,0:7)

  real*8    :: dinvmo, dinvmg, dinvmq, dnrmk, dnrmg, diffg, diffq, diffa
  real*8    :: timend, time, timstr, second, toterg, totdiag, totnond

  integer*4, allocatable, dimension(:)       :: nfnuc, isel, indao, indaos, indnao
  real*8,    allocatable, dimension(:)       :: scoef, scoefo, smat, cmo
  real*8,    allocatable, dimension(:,:)     :: gradient
  real*8,    allocatable, dimension(:,:,:)   :: xkappa, hmat
  real*8,    allocatable, dimension(:,:,:,:) :: pmat

  !------------------------------------------------------------------------------!
  call qenter('orbloc')

  ifrp   = 1
  lupri  = var%lupri
  nish   = var%nish(ifrp)
  norb   = var%norb(ifrp)
  nucdep = var%nucdep
  maxref = var%maxref
  full   = loc%auto
  diagonly = .false.
  if(loc%hesthr < loc%thrloc) diagonly = .true.

  ! TMP print
  write(6,*) 'lupri =',lupri
  write(6,*) 'nish  =',nish
  write(6,*) 'norb  =',norb
  write(6,*) 'nucdep=',nucdep
  write(6,*) 'maxref=',maxref
  write(6,*) 'full  =',full
  write(6,*) 'hesthr=',loc%hesthr
  write(6,*) 'thrloc=',loc%thrloc

  nfbas(:)    = var%nfbas(ifrp,:)
  ipqtoq(:,:) = var%ipqtoq(:,:)

  call titler('Molecular orbital localization','*',122)

  ! Check if coefficients are on file
  inquire(file='DFCOEF',exist=tobe)
  if(.not.tobe)then
    write(lupri,'(a)') 'orbloc: Coefficient file not found !'
    call quit('orbloc: Coefficients not found !')
  endif

  ! Memory allocation
  nish12 = 4 * ( nish*nish + nish ) / 2
  call alloc(scoef   ,nish*nfbas(0)*4,   id='scoef    in orbloc')
  call alloc(scoefo  ,nish*nfbas(0)*4,   id='scoefo   in orbloc')
  call alloc(smat    ,var%n2bbasx,       id='smat     in orbloc')
  call alloc(cmo     ,var%n2bbasxq,      id='cmo      in orbloc')
  call alloc(hmat    ,nish12,nish12,2,   id='hmat     in orbloc')
  call alloc(gradient,nish12,2,          id='gradient in orbloc')
  call alloc(xkappa  ,nish,nish,4,       id='xkappa   in orbloc')
  call alloc(pmat    ,nish,nish,4,nucdep,id='pmat     in orbloc')

  call alloc(nfnuc ,2*nucdep,  id='nfnuc  in orbloc')
  call alloc(isel  ,nish,      id='isel   in orbloc')
  call alloc(indao ,nfbas(0),  id='indao  in orbloc')
  call alloc(indaos,nfbas(0),  id='indaos in orbloc')
  call alloc(indnao,2*nucdep+2,id='indnao in orbloc')

  ! Read coefficients
  call opnfil(var%lucoef,'DFCOEF','OLD','PAMANA')
  call reacmo_no_work(var%lucoef,'DFCOEF',cmo,idum,idum,toterg,2)
  close(var%lucoef,status='keep')

  ! Get overlap matrix 
  call gtovlx(smat,var%ssmtrc)
  if(loc%iprloc.ge.6) then
    write(lupri,'(/a)') '* orbloc: Molecular coefficients:'
    call prqmat(cmo,nfbas(0),norb,nfbas(0),norb,4,ipqtoq(1,0),lupri)

    write(lupri,'(/a)') '* orbloc: Overlap matrix' 
    call prqmat(smat,nfbas(0),nfbas(0),nfbas(0),nfbas(0),1,ipqtoq(1,0),lupri)
  endif

  ! Occupied electronic orbital coeffients selection
  call occupy(isel,nish)
  call selcfs(cmo,ifrp,scoef,nish,isel,d0,nish,nfbas(0),norb)
  if(loc%iprloc.ge.5) then
    write(lupri,'(/a)') '* orblo1: Selected molecular coefficients'
    call prqmat(scoef,nfbas(0),nish,nfbas(0),nish,4,ipqtoq(1,0),lupri)
  endif

  ! Get the number of functions per center
  call indaoc(indao,indaos,indnao)
  call orblo4(indnao,nfnuc,nucdep)

  timstr = second()
  if(loc%prjloc) then
    ! PRJLOC = TRUE localization with projection analysis
    call orblp1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,maxref,dinvmg)
  else
    ! PRJLOC = FALSE localization with Mulliken analysis
    call orblo1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,nfnuc,dinvmg)
  endif
  write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 5x,a, 6x,a, 6x,a)') &
    '*C','ITR','DIFFA','|GRAD|','|KAPPA|','H_DIAG','H_NONDIAG','H','DIFFG','DIFFQ'

  ! SCF cycle
  dinvmq = 0.0d0
  do itr=1,loc%itrloc

    dinvmo = dinvmg
    call gtinfo(daytid)
    write(lupri,'(/5x,a,i4,a,3x)')'########## START ITERATION NO. ',itr,' ##########'
    write(lupri,'(5x,a,a24,a)')'########## ',daytid,' ##########'

    call get_hessian(pmat,hmat,totdiag,totnond)
    call get_gradient(pmat,gradient,dnrmg)
    if(full)then
      call getk_full(pmat,xkappa,hmat,gradient,dnrmk)
      call exp_func(xkappa,hmat,gradient,full,dinvmq)
    else
      call getk_cdiag(pmat,xkappa,gradient,dnrmk,ipqtoq)
!     call getk_rdiag(pmat,xkappa,gradient,dnrmk,ipqtoq)
    endif
    call rotcmo_mod(scoef,nfbas(0),nish,4,xkappa,nish,ipqtoq(ifrp,0))
    if(loc%prjloc) then
      ! PRJLOC = TRUE localization with projection analysis
      call gtovlx(smat,var%ssmtrc)
      call orblp1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,maxref,dinvmg)
    else
      ! PRJLOC = FALSE localization with Mulliken analysis
      call orblo1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,nfnuc,dinvmg)
    endif

    ! Convergence criterion
    diffg = dinvmg - dinvmo
    diffq = dinvmq
    if (diffg < d0) write(lupri,'(a,i4)')' WARNING: dinvm_new - dinvm_old < 0; itr =',itr
    diffa = dabs(dinvmg - dinvmo)
    if(full) hessian = 'FULL'

    write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 5x,a, 6x,a, 6x,a)') &
      '  ','ITR','DIFFA','|GRAD|','|KAPPA|','H_DIAG','H_NONDIAG','H','DIFFG','DIFFQ'
    write(lupri,'(a,2x,i4,5d11.2,4x,a,2d11.2)') &
      '*C',itr,diffa,dnrmg,dnrmk,totdiag,totnond,hessian,diffg,diffq

    ! Exit loop after convergence criterion is reached
    if((diffa < loc%thrloc).and.(full.or.diagonly)) exit
    if( diffa < loc%hesthr ) full = .true.

  enddo

  if(diffa < loc%thrloc) conv = 'Convergence in      '
  write(lupri,'(/a,a20,i5,a)') ' >>> ',conv,itr,' iterations'

  timend = second()
  time   = timend - timstr
  call timtxt('>>> Time used in localization is :',time,lupri)

  ! Put the localised coefficients in DFCOEF
  call orblo5(scoef,nfbas(0),nish,norb)

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
  call dealloc(scoef   ,id='scoef    in orbloc')
  call dealloc(scoefo  ,id='scoefo   in orbloc')
  call dealloc(smat    ,id='smat     in orbloc')
  call dealloc(pmat    ,id='pmat     in orbloc')
  call dealloc(xkappa  ,id='xkappa   in orbloc')
  call dealloc(hmat    ,id='hmat     in orbloc')
  call dealloc(gradient,id='gradient in orbloc')
  call dealloc(cmo     ,id='cmo      in orbloc')

  call dealloc(nfnuc   ,id='nfnuc    in orbloc')
  call dealloc(isel    ,id='isel     in orbloc')
  call dealloc(indao   ,id='indao    in orbloc')
  call dealloc(indaos  ,id='indaos   in orbloc')
  call dealloc(indnao  ,id='indnao   in orbloc')

  call qexit('orbloc')

  end subroutine orbloc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLP1(SMAT,SCOEF,PMAT,NBL,NBS,NI,NUCDEP,MAXREF,DINVM)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      USE SELECT_COEFF_OWN_BAS
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NBL, NBS, NI, NUCDEP, MAXREF
      REAL*8    :: DINVM
      REAL*8    :: SMAT(*), SCOEF(*), PMAT(NI,NI,4,NUCDEP)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      REAL*8, PARAMETER :: D0=0.0D0, D2=2.0D0, DM1=-1.0D0, D1=1.0D0

      LOGICAL   :: TOBE

      INTEGER*4 :: I, J, K, IZ, IER, IFRP, ICFS, NB, NQ, NCDIM, NBDIM, NMOL
      INTEGER*4 :: IVEC, JOFF, NREF, N2REF, LUPRI, NPSH, NESH, NORB, NREFS
      INTEGER*4 :: NSTR(2,0:2,0:MAXREF), NVECS(2,0:2), IOFF(2), IDUMARR(2), IPQTOQ(4,0:7)

      REAL*8    :: DDUMARR(2), DINV, QLOC2, DLOC

      INTEGER*4, ALLOCATABLE, DIMENSION(:) :: KVECTMP, JREF, IPIV
      REAL*8,    ALLOCATABLE, DIMENSION(:) :: RCOEF, CMO, BV, GVEC, AMAT, WMAT

      TYPE :: PTR2ARRAY
        INTEGER*4, DIMENSION(:), POINTER :: P
      END TYPE PTR2ARRAY

      TYPE(PTR2ARRAY), ALLOCATABLE :: KVEC(:,:)

  !------------------------------------------------------------------------------!
      IFRP   = 1
      LUPRI= VAR%LUPRI
      NPSH = VAR%NPSH(IFRP)
      NESH = VAR%NESH(IFRP)
      NORB = VAR%NORB(IFRP)
      NREFS= VAR%NREFS
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)
      VECREF(:,:) = VAR%VECREF(:,:)
      REFFIL(:)   = VAR%REFFIL(:)

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
      NB     = NBL + NBS
!     DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        CALL ALLOC(KVECTMP,NORB,ID='KVECTMP IN ORBLP1')
        DO I = 1,NREFS
          NSTR(IFRP,0,I) = 1
          CALL NUMLST(VECREF(IFRP,I),KVECTMP, &
               NORB,-NPSH,NESH,               &
               IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(KVECTMP,NSTR(IFRP,0,I), &
               NPSH,NESH,                     &
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
      IF(VAR%IPRPRJ.GE.5) ICFS = 7
      CALL ALLOC(CMO,VAR%N2BBASXQ,ID="CMO in ORBLP1")

      IOFF(1) = 0
      VAR%IPRPRJ = 4
      DO I = 1,NREFS
        INQUIRE(FILE=REFFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
          WRITE(LUPRI,'(A,A6,A)') 'Reference coefficient file ',REFFIL(I),' not present'
          CALL QUIT('PRJANA: No reference coefficients !')
        ENDIF
        CALL OPNFIL(VAR%LUCOEF,REFFIL(I),'OLD','PRJANA')
        IF(VAR%OWNBAS)THEN
          CALL SELOWN_DRV(2,CMO,DDUMARR,IDUMARR,               &
               I,1,VAR%LUCOEF,REFFIL(I),RCOEF,DDUMARR,IDUMARR, &
               NVECS(1,1),IOFF,KVEC(IFRP,I)%P,                 &
               NSTR(1,0,I))
          IF(LOC%IPRLOC .GE. 5) THEN
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
        CLOSE(VAR%LUCOEF,STATUS='KEEP')
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
           JREF,NSTR,IPIV,MAXREF)

! Memory allocation for the matrix W = Ac

      CALL ALLOC(WMAT,NREF*NI*4,ID="WMAT in ORBLP1")

! Construction of W = Ac

      IF(LOC%IPRLOC .GE. 5) THEN
        WRITE(LUPRI,'(/A)') '* Matrice BV'
        CALL PRQMAT(BV,NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A)') '* Matrice AMAT'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NREF,NI,NREF,D1,                &
           'N','N',IPQTOQ(1,0),AMAT,NREF,NREF,4, &
           'N','N',IPQTOQ(1,0),BV,NREF,NI,4,     &
           D0,IPQTOQ(1,0),WMAT,NREF,NI,4)

      IF(LOC%IPRLOC.GE.5) THEN
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
        WRITE(LUPRI,*)
        WRITE(LUPRI,'(A,I3)') 'Print Q^A matrices, IZ = ',IZ
        DO I=1, NUCDEP
          WRITE(LUPRI,'(A,I4)') ' Nucleus #',I
          DO J=1,NI
            WRITE(LUPRI,'(I4,10ES9.1)') J,(PMAT(J,K,IZ,I),K=1,NI)
          ENDDO
          WRITE(LUPRI,*)
        ENDDO
      ENDDO

      IF(LOC%IPRLOC .GE. 1) THEN
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
          QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)
        ENDDO
        DLOC = D1/QLOC2
        DINV = DINV + QLOC2
        WRITE(LUPRI,'(A,I3,A,2F7.4)') &
        '* Orbital ', I ,' Pipek-Mezey localization measure d = ',DLOC
      ENDDO
      DINVM = DINV/DFLOAT(NI)
      WRITE(LUPRI,'(/A,2F17.13)') '* Inverse mean delocalization = ',DINVM

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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

 SUBROUTINE ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NMOL,NREF,BVEC,GVEC, &
                   JMOL,JREF,NSTR,IPIVOT,MAXREF)
!***********************************************************************
!
!     Projection analysis for fermion ircop IFRP
!
!     Written by T.Saue March 26 1999
!
!***********************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NB, NMOL, NREF, MAXREF
      INTEGER*4 :: JMOL(*), JREF(3,*), NSTR(2,0:2,0:MAXREF), IPIVOT(*)
      REAL*8    :: SMAT(*), SCOEF(*), RCOEF(*), AMAT(NREF,NREF,4), BVEC(NREF,NMOL,4), GVEC(*)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      INTEGER*4    :: I, J, IFRP, JOFF, KOFF, JOB, NEFF, NREFS, LUPRI
      INTEGER*4    :: IPQTOQ(4,0:7)

      REAL*8    :: PTOL
      REAL*8, PARAMETER :: D0=0.0D0, DC=1.0D-2, D2=2.0D0

      LOGICAL DOPOP

  !------------------------------------------------------------------------------!
      IFRP = 1
      NREFS = VAR%NREFS

      DOPOP = NREFS.GT.1
      PTOL  = VAR%PROTHR*DC
      LUPRI = VAR%LUPRI
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)
      VECREF(:,:) = VAR%VECREF(:,:)
      REFFIL(:)   = VAR%REFFIL(:)

!     Make projection vector

      CALL QTRANS90('AOMO','S',D0,NB,NB,NREF,NMOL, &
                    SMAT,NB,NB,1,IPQTOQ(1,0),      &
                    BVEC,NREF,NMOL,4,IPQTOQ(1,0),  &
                    RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                    SCOEF,NB,NMOL,4,IPQTOQ(1,0),VAR%IPRPRJ)
      IF(NREFS.GT.1) THEN

!       Make overlap matrix

        CALL QTRANS90('AOMO','S',D0,NB,NB,NREF,NREF, &
                      SMAT,NB,NB,1,IPQTOQ(1,0),      &
                      AMAT,NREF,NREF,4,IPQTOQ(1,0),  &
                      RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                      RCOEF,NB,NREF,4,IPQTOQ(1,0),VAR%IPRPRJ)

!       Print overlaps

        IF(LOC%IPRLOC .GE. 5) THEN
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

        CALL QHMRST(AMAT,NREF,VAR%NZ,NREF,NREF)
      ENDIF

      END SUBROUTINE ORBLP2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO1(SMAT,SCOEF,PMAT,NBL,NBS,NI,NUCDEP,NFNUC,DINVM)
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
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NBL, NBS, NI, NUCDEP
      REAL*8    :: DINVM
      INTEGER*4 :: NFNUC(*)
      REAL*8    :: SMAT(*), SCOEF(*), PMAT(NI,NI,4,NUCDEP)

      ! Internal variables
      INTEGER*4 :: I, J, IZ, NB, IOFFL, IOFFS, LUPRI
      INTEGER*4 :: IPQTOQ(4,0:7)
      REAL*8    :: DINV, QLOC2, DLOC

      REAL*8, PARAMETER :: D0 = 0.0D0, D1 = 1.0D0

      REAL*8, ALLOCATABLE, DIMENSION(:) :: WLMAT, WSMAT

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLO1')

      NB = NBL + NBS

      LUPRI       = VAR%LUPRI
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)

! Memory allocation for the matrix W = Sc

      CALL ALLOC(WLMAT,NBL*NI*4,ID='WLMAT IN ORBLO1')
      CALL ALLOC(WSMAT,NBS*NI*4,ID='WSMAT IN ORBLO1')

! Construction of W = Sc (large and small components)

      CALL QGEMM(NBL,NI,NBL,D1,               &
           'N','N',IPQTOQ(1,0),SMAT,NB,NB,1,  &
           'N','N',IPQTOQ(1,0),SCOEF,NB,NI,4, &
           D0,IPQTOQ(1,0),WLMAT,NBL,NI,4)
      IF(.NOT.VAR%TWOCOMPBSS) THEN
         CALL QGEMM(NBS,NI,NBS,D1,                            &
              'N','N',IPQTOQ(1,0),SMAT(1+NB*NBL+NBL),NB,NB,1, &
              'N','N',IPQTOQ(1,0),SCOEF(1+NBL),NB,NI,4,       &
              D0,IPQTOQ(1,0),WSMAT,NBS,NI,4)
      ENDIF

      IF(LOC%IPRLOC.GE.5) THEN
         WRITE(LUPRI,'(/A)') '* WL matrix :'
         CALL PRQMAT(WLMAT,NBL,NI,NBL,NI,4,IPQTOQ(1,0),LUPRI)
         IF(.NOT.VAR%TWOCOMPBSS) THEN
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
         IF(.NOT.VAR%TWOCOMPBSS) THEN
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
      IF(LOC%IPRLOC.GE.1) THEN
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine getk_cdiag(pmat,xkappa,gradient,dnrm,ipqtoq)
  !******************************************************************************!
  !                                                                              !
  !  Construction of the XKAPPA matrix using diagonal approximation of hessian   !
  !                             (complex algebra)                                !
  !                                                                              !
  ! Notes: 1. Deffinition of kappa matrix in rotcmo_mod is different as in       !
  !           calculation of gradient and hessian in get_gradient and get_hessian!
  !           --> Indexes are exchanged and j quaterion is with minus sign       !
  !               Qkappa_lk = kappa_kl - bar{kappa}_kl^* j                       !
  !               It follows from definition of bared unitary matrix:            !
  !               tilde{Psi_i} = U_ij Psi_j + bar{U}_ij Psi_j                    !
  !        2. Be also careful about complex cojnugation of gradient in           !
  !                                H * K = - g*                                  !
  !                                                                              !
  !******************************************************************************!
  use memory_allocator
  implicit none

  ! External variables
  integer*4, dimension(:,:) :: ipqtoq

  real*8 :: dnrm

  real*8, dimension(:,:)     :: gradient
  real*8, dimension(:,:,:)   :: xkappa
  real*8, dimension(:,:,:,:) :: pmat

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

  integer*4 :: i, j, k, ij, ij1, iz, nucdep, ni, ni12, lupri
  real*8    :: pii, pjj, pn2, dnrm2

  real*8, allocatable, dimension(:) :: hmat_d

  !------------------------------------------------------------------------------!
  call qenter('getk_cdiag')

  ! Get some dimensions
  lupri  = var%lupri
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ni12   = size(gradient,1) / 4

  call alloc(hmat_d,2*ni12,id='hmat_d in getk_cdiag')

  ! Construct diagonal Hessian matrix
  hmat_d(:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      if(i==j)cycle
      do k = 1,nucdep
        pii = pmat(i,i,1,k); pjj = pmat(j,j,1,k)

        pn2 = pmat(i,j,1,k)*pmat(i,j,1,k) + pmat(i,j,2,k)*pmat(i,j,2,k)
        hmat_d(ij) = hmat_d(ij) + d4*pn2 - d2*((pii-pjj)*(pii-pjj))

        pn2 = pmat(i,j,3,k)*pmat(i,j,3,k) + pmat(i,j,4,k)*pmat(i,j,4,k)
        hmat_d(ij+ni12) = hmat_d(ij+ni12) + d4*pn2 - d2*((pii-pjj)*(pii-pjj))
      enddo
    enddo
  enddo

  ! Construction and print of the kappa matrix
  xkappa(:,:,:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      ij1= ij + ni12
      if(i==j) cycle
      xkappa(j,i,1) = - gradient(ij ,1) / hmat_d(ij )
      xkappa(j,i,2) =   gradient(ij ,2) / hmat_d(ij )   !!! gradient is complex conjugate
      xkappa(j,i,3) =   gradient(ij1,1) / hmat_d(ij1)   !!! j component is with minus sign
      xkappa(j,i,4) =   gradient(ij1,2) / hmat_d(ij1)   !!! gradient is complex conjugate
      xkappa(i,j,1) = - xkappa(j,i,1)
      xkappa(i,j,2) =   xkappa(j,i,2)
      xkappa(i,j,3) =   xkappa(j,i,3)
      xkappa(i,j,4) =   xkappa(j,i,4)
    enddo
  enddo

  do iz=1,4
    write(lupri,*)
    write(lupri,'(a,i4)') 'Print xkappa matrix, IZ =',iz
    do j=1,ni
      write(lupri,'(i4,10es9.1)') j,(xkappa(j,k,iz),k=1,ni)
    enddo
  enddo

  ! Print some kappa matrix
  if(loc%iprloc .ge. 5) then
     write (lupri,'(/a)') '* KAPPA matrix :'
     call prqmat(xkappa(1,1,1),ni,ni,ni,ni,4,ipqtoq(1,0),lupri)
  endif

  ! Norm of the kappa vector
  dnrm = dnrm2(ni*ni*4,xkappa,1)

  call dealloc(hmat_d,id='hmat_d in getk_cdiag')

  call qexit('getk_cdiag')

  end subroutine getk_cdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine getk_rdiag(pmat,xkappa,gradient,dnrm,ipqtoq)
  !******************************************************************************!
  !                                                                              !
  !  Construction of the XKAPPA matrix using diagonal approximation of hessian   !
  !                              (real algebra)                                  !
  !                                                                              !
  ! Notes: 1. Deffinition of kappa matrix in rotcmo_mod is different as in       !
  !           calculation of gradient and hessian in get_gradient and get_hessian!
  !           --> Indexes are exchanged and j quaterion is with minus sign       !
  !               Qkappa_lk = kappa_kl - bar{kappa}_kl^* j                       !
  !               It follows from definition of bared unitary matrix:            !
  !               tilde{Psi_i} = U_ij Psi_j + bar{U}_ij Psi_j                    !
  !        2. Be also careful about minus in imaginary part of gradient in       !
  !                                H * K = - g                                   !
  !                                 H                     K            g         !
  !  / A   B  \   -->   / A^R + B^I   B^I - A^I \      / K^R \      / g^R \      !
  !  \ B*  A* /   -->   \ B^I + A^I   A^R - B^I /      \ K^I /      \-g^I /      !
  !******************************************************************************!
  use memory_allocator
  implicit none

  ! External variables
  integer*4, dimension(:,:) :: ipqtoq

  real*8 :: dnrm

  real*8, dimension(:,:)     :: gradient
  real*8, dimension(:,:,:)   :: xkappa
  real*8, dimension(:,:,:,:) :: pmat

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

  integer*4 :: i, j, k, ij, ij1, ij2, ij3, iz, nucdep, ni, ni12, lupri
  real*8    :: pii, pjj, pn2, dnrm2

  real*8, allocatable, dimension(:) :: hmat_d

  !------------------------------------------------------------------------------!
  call qenter('getk_rdiag')

  ! Get some dimensions
  lupri  = var%lupri
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ni12   = size(gradient,1) / 4

  call alloc(hmat_d,4*ni12,id='hmat_d in getk_rdiag')

  ! Construct diagonal Hessian matrix
  hmat_d(:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      ij1= ij + 1*ni12
      ij2= ij + 2*ni12
      ij3= ij + 3*ni12
      if(i==j)cycle
      do k = 1,nucdep
        pii = pmat(i,i,1,k); pjj = pmat(j,j,1,k)

        pn2 = d2*pmat(i,j,1,k)*pmat(i,j,1,k)
        hmat_d(ij ) = hmat_d(ij ) + d4*pn2 - d2*((pii-pjj)*(pii-pjj))  ! uu: A^R + B^R

        pn2 = d2*pmat(i,j,2,k)*pmat(i,j,2,k)
        hmat_d(ij2) = hmat_d(ij2) + d4*pn2 - d2*((pii-pjj)*(pii-pjj))  ! uu: A^R - B^R

        pn2 = pmat(i,j,3,k)*pmat(i,j,3,k) + pmat(i,j,4,k)*pmat(i,j,4,k)
        hmat_d(ij1) = hmat_d(ij1) + d4*pn2 - d2*((pii-pjj)*(pii-pjj))  ! bb: A^R + B^R
        hmat_d(ij3) = hmat_d(ij3) + d4*pn2 - d2*((pii-pjj)*(pii-pjj))  ! bb: A^R - B^R
      enddo
    enddo
  enddo

  ! Construction and print of the kappa matrix
  xkappa(:,:,:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      ij1= ij + 1*ni12
      ij2= ij + 2*ni12
      ij3= ij + 3*ni12
      if(i==j) cycle
      xkappa(j,i,1) = - gradient(ij ,1) / hmat_d(ij )
      xkappa(j,i,2) =   gradient(ij ,2) / hmat_d(ij2)   !!! g^I is with minus
      xkappa(j,i,3) =   gradient(ij1,1) / hmat_d(ij1)   !!! j component is with minus sign
      xkappa(j,i,4) =   gradient(ij1,2) / hmat_d(ij3)   !!! g^I is with minus
      xkappa(i,j,1) = - xkappa(j,i,1)
      xkappa(i,j,2) =   xkappa(j,i,2)
      xkappa(i,j,3) =   xkappa(j,i,3)
      xkappa(i,j,4) =   xkappa(j,i,4)
    enddo
  enddo

  do iz=1,4
    write(lupri,*)
    write(lupri,'(a,i4)') 'Print xkappa matrix, IZ =',iz
    do j=1,ni
      write(lupri,'(i4,10es9.1)') j,(xkappa(j,k,iz),k=1,ni)
    enddo
  enddo

  ! Print some kappa matrix
  if(loc%iprloc .ge. 5) then
     write (lupri,'(/a)') '* KAPPA matrix :'
     call prqmat(xkappa(1,1,1),ni,ni,ni,ni,4,ipqtoq(1,0),lupri)
  endif

  ! Norm of the kappa vector
  dnrm = dnrm2(ni*ni*4,xkappa,1)

  call dealloc(hmat_d,id='hmat_d in getk_rdiag')

  call qexit('getk_rdiag')

  end subroutine getk_rdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine getk_full(pmat,xkappa,hmat,gradient,dnrm)
  !******************************************************************************!
  !                                                                              !
  !  Construction of the XKAPPA matrix using full hessian                        !
  !                                                                              !
  !******************************************************************************!
  use memory_allocator
  implicit none

  ! External variables
  real*8 :: dnrm

  real*8, dimension(:,:)     :: gradient
  real*8, dimension(:,:,:)   :: xkappa, hmat
  real*8, dimension(:,:,:,:) :: pmat

  ! Internal variables
  real*8, parameter     :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0
  complex*16, parameter :: c2 = (2.0d0,0.0d0)

  integer*4 :: i, j, i1, j1, k, ij, kl, iz, nucdep, ni, ndimk, ndimh, ni12, lupri, ier
  real*8    :: pii, pjj, pn2, dnrm2, dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6

  integer*4,  allocatable, dimension(:)   :: ipiv
  real*8, allocatable, dimension(:)   :: grad_c
  real*8, allocatable, dimension(:,:) :: hmat_c

  !------------------------------------------------------------------------------!
  call qenter('getk_full')

  lupri  = var%lupri

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimh  = size(hmat,1)
  ni12   = ndimh / 4

  ! TMP: check properties of hessian
  dtmp1 = d0; dtmp2 = d0; dtmp3 = d0; dtmp4 = d0; dtmp5 = d0; dtmp6 = d0
  do i=1,2*ni12
    do j=i,2*ni12
      if(i==j)cycle
      dtmp1 = dtmp1 + dabs(hmat(i,j,1) - hmat(j,i,1))
      dtmp2 = dtmp2 + dabs(hmat(i,j,2) + hmat(j,i,2))
      dtmp3 = dtmp3 + dabs(hmat(i,j+2*ni12,1) - hmat(j,i+2*ni12,1))
      dtmp4 = dtmp4 + dabs(hmat(i,j+2*ni12,2) - hmat(j,i+2*ni12,2))
      dtmp5 = dtmp5 + dabs(hmat(i,j,1) - hmat(i+2*ni12,j+2*ni12,1)) &
                    + dabs(hmat(j,i,1) - hmat(j+2*ni12,i+2*ni12,1)) &
                    + dabs(hmat(i,j,2) + hmat(i+2*ni12,j+2*ni12,2)) &
                    + dabs(hmat(j,i,2) + hmat(j+2*ni12,i+2*ni12,2))
      dtmp6 = dtmp6 + dabs(hmat(i,j+2*ni12,1) - hmat(i+2*ni12,j,1)) &
                    + dabs(hmat(j,i+2*ni12,1) - hmat(j+2*ni12,i,1)) &
                    + dabs(hmat(i,j+2*ni12,2) + hmat(i+2*ni12,j,2)) &
                    + dabs(hmat(j,i+2*ni12,2) + hmat(j+2*ni12,i,2))
    enddo
    dtmp2 = dtmp2 + dabs(hmat(i,i,2))
    dtmp5 = dtmp5 + dabs(hmat(i,i,1) - hmat(i+2*ni12,i+2*ni12,1))
    dtmp6 = dtmp6 + dabs(hmat(i,i+2*ni12,1) - hmat(i+2*ni12,i,1)) &
                  + dabs(hmat(i,i+2*ni12,2) + hmat(i+2*ni12,i,2))
  enddo
  write(6,*) 'symetry A_R',dtmp1
  write(6,*) 'anti-symetry A_I',dtmp2
  write(6,*) 'symetry B_R',dtmp3
  write(6,*) 'symetry B_I',dtmp4
  write(6,*) 'A and A*',dtmp5
  write(6,*) 'B and B*',dtmp6

  call alloc(ipiv  ,ndimh-ni,id='ipiv   in getk_full')
  allocate( grad_c(ndimh)       ,stat=ier); call errall(ier,'getk_full','grad_c')
  allocate( hmat_c(ndimh,ndimh) ,stat=ier); call errall(ier,'getk_full','hmat_c')

  ! prepare hmat_c
  do i=1,2*ni12
    do j=1,2*ni12
      i1 = i + 2*ni12
      j1 = j + 2*ni12
      hmat_c(i ,j ) = hmat(i,j ,1) + hmat(i,j1,1)
      hmat_c(i ,j1) = hmat(i,j1,2) - hmat(i,j ,2)
      hmat_c(i1,j ) = hmat(i,j1,2) + hmat(i,j ,2)
      hmat_c(i1,j1) = hmat(i,j ,1) - hmat(i,j1,1)
    enddo
  enddo

  do ij=1,ndimh
    kl = 1
    do k=1,ni
      do j=kl+1,ndimh-k+1
        hmat_c(ij,j-1) = hmat_c(ij,j)
      enddo
      kl = kl + ni - (k-1) - 1 ! diagonal
    enddo
  enddo

  do ij=1,ndimh
    kl = 1
    do k=1,ni
      do j=kl+1,ndimh-k+1
        hmat_c(j-1,ij) = hmat_c(j,ij)
      enddo
      kl = kl + ni - (k-1) - 1 ! diagonal
    enddo
  enddo

  ! prepare grad_c
  do i=1,2*ni12
    i1 = i + 2*ni12
    grad_c(i ) = - gradient(i,1)
    grad_c(i1) =   gradient(i,2)
  enddo

  kl = 1
  do k=1,ni
    do j=kl+1,ndimh-k+1
      grad_c(j-1) = grad_c(j)
    enddo
    kl = kl + ni - (k-1) - 1 ! diagonal
  enddo

  ! Solve linear set of equations
  call dgetrf(ndimh-ni,ndimh-ni,hmat_c,ndimh,ipiv,ier)
  if(ier /= 0)then
    write(lupri,*) ' ERROR: The dgetrf LAPACK subroutine in getk_full finished unsuccessfully'
    write(lupri,*) ' info =',ier
    stop
  endif

  call dgetrs('N',ndimh-ni,1,hmat_c,ndimh,ipiv,grad_c,ndimh,ier)
  if(ier /= 0)then
    write(lupri,*) ' ERROR: The dgetrs LAPACK subroutine in getk_full finished unsuccessfully'
    write(lupri,*) ' info =',ier
    stop
  endif

  xkappa(:,:,:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      if(i==j)cycle
      ij = ij + 1
      xkappa(j,i,1) =   grad_c(ij)
      xkappa(i,j,1) = - xkappa(j,i,1)
    enddo
  enddo

  ij = ni12-ni
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      xkappa(j,i,2) =   grad_c(ni12+ij)
      xkappa(i,j,2) =   xkappa(j,i,2)

      xkappa(j,i,3) = - grad_c(ij)
      xkappa(i,j,3) =   xkappa(j,i,3)

      xkappa(j,i,4) =   grad_c(2*ni12+ij)
      xkappa(i,j,4) =   xkappa(j,i,4)
    enddo
  enddo

  do iz=1,4
    write(lupri,*)
    write(lupri,'(a,i4)') 'Print xkappa matrix, IZ =',iz
    do j=1,ni
      write(lupri,'(i4,10es9.1)') j,(xkappa(j,k,iz),k=1,ni)
    enddo
  enddo

  ! Norm of the kappa
  ndimk = ni*ni*4
  dnrm  = dnrm2(ndimk,xkappa,1)
  write (lupri,'(a,e12.2)') 'Norm of kappa: ',dnrm

  call dealloc(ipiv  ,id='ipiv   in getk_full')
  deallocate(grad_c)
  deallocate(hmat_c)

  call qexit('getk_full')

  end subroutine getk_full

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO4(INDNAO,NFNUC,NUCDEP)
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
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NUCDEP
      INTEGER*4 :: INDNAO(2,*), NFNUC(2*NUCDEP)

      ! Internal variables
      INTEGER*4 :: IFRP, I
      REAL*8, PARAMETER :: D1 = 1.0D0

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
      NFNUC(2*NUCDEP)=VAR%NFBAS(IFRP,0)-INDNAO(2,NUCDEP)
      CALL QEXIT('ORBLO4')

  END SUBROUTINE ORBLO4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO5(SCOEF,NB,NI,NO)
!***********************************************************************
!
!     Written by S.Dubillard Apr 13 2004
!     Last revision:
!
!     Put the localised coefficients into the DFCOEF file (1/2)
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NB, NI, NO
      REAL*8    :: SCOEF(NB,NI,4)

      ! Internal variables
      INTEGER*4 :: I, J, K, L
      REAL*8    :: TOTERG
      REAL*8, PARAMETER :: D0 = 0.0D0, D1 = 1.0D0

      INTEGER*4, ALLOCATABLE :: IBEIG(:)
      REAL*8,    ALLOCATABLE :: CMO(:,:,:), EIG(:)

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLO5')

      ! allocate arrays
      CALL ALLOC(CMO  ,NB,NO,4,ID='CMO   IN ORBLO5')
      CALL ALLOC(EIG  ,VAR%NORBT  ,ID='EIG   IN ORBLO5')
      CALL ALLOC(IBEIG,VAR%NORBT  ,ID='IBEIG IN ORBLO5')

      ! open DFCOEF file and read EIG, IBEIG and TOTERG
      CALL OPNFIL(VAR%LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO_no_work(VAR%LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)

      ! reorganize MO
      DO I = 1,NB
        DO J = 1,NI
          DO K = 1,4
            L = J + NO/2
            IF (VAR%LEVYLE) L = J
            CMO(I,L,K) = SCOEF(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      IF (LOC%IPRLOC .GE. 5) THEN
         WRITE(VAR%LUPRI,'(/A)') '*ORBLO6 : New CMO matrix'
         CALL PRQMAT(CMO,NB,NO,NB,NO,4,VAR%IPQTOQ(1,0),VAR%LUPRI)
      ENDIF

      ! write MO to file and close him
      CALL WRICMO(VAR%LUCOEF,CMO,EIG,IBEIG,TOTERG)
      CLOSE(VAR%LUCOEF,STATUS='KEEP')

      ! deallocate arrays
      CALL DEALLOC(CMO  ,ID='CMO   IN ORBLO5')
      CALL DEALLOC(EIG  ,ID='EIG   IN ORBLO5')
      CALL DEALLOC(IBEIG,ID='IBEIG IN ORBLO5')

      CALL QEXIT('ORBLO5')

  END SUBROUTINE ORBLO5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ROTCMO_MOD(CMO,NBASI,NORBI,NZ,XKAPPA,NORBT,IP)
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
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NBASI, NORBI, NZ, NORBT, IP(4)
      REAL*8    :: CMO(NBASI,NORBI,NZ), XKAPPA(NORBT,NORBT,NZ)

      ! Internal variables
      REAL*8, PARAMETER :: D1 = 1.00D00, D0 = 0.00D00, DM1 = -D1
      REAL*8, PARAMETER :: DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00
      REAL*8, PARAMETER :: D1P5 = 1.5D00, D2 = 2.0D00, THREQL = 1.0D-12

      CHARACTER SECTID*12,CPUTID*12,WALLTID*12

      INTEGER*4 :: I, J, NORBXQ, IERR, IADDR, IOFF, I_MAX, J_MAX, LUPRI
      REAL*8    :: CPU1, WALL1, CPU2, WALL2, CPU3, WALL3, CPU4, WALL4, EIG, SQEIG, D_MAX
      REAL*8    :: ABS, SQRT, COS, SIN

      REAL*8, ALLOCATABLE, DIMENSION(:) :: REIG, RSIN, RCOS, RV, RTMP

  !------------------------------------------------------------------------------!
      CALL QENTER('ROTCMO_MOD')

      LUPRI  = VAR%LUPRI

      CALL GETTIM(CPU1,WALL1)
      IF (LOC%IPRLOC .GE. 5) CALL HEADER('Output from ROTCMO_MOD',-1)

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

      IF (LOC%IPRLOC .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: -kappa^2',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Diagonalize -\kappa^2
!     ---------------------------------------------

      CALL QDIAG90(NZ,NORBI,RCOS,NORBI,NORBI, &
                   REIG,1,RV,NORBI,NORBI,IERR)
      IF (LOC%IPRLOC .GE. 10) THEN
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
            IF (LOC%IPRLOC .GE. 3 .AND. ABS(EIG).GT. VAR%THRZER) WRITE(LUPRI,9010) EIG
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
      IF (LOC%IPRLOC .GE. 10) THEN
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
                    RV,NORBI,NORBI,NZ,IP,LOC%IPRLOC)
      IF (LOC%IPRLOC .GE. 10) THEN
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
                    RV,NORBI,NORBI,NZ,IP,LOC%IPRLOC)
      IF (LOC%IPRLOC .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
!     ---------------------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,DM1,            &
                 'N','N',IP,RCOS,NORBI,NORBI,NZ,   &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D1,IP,RTMP,NORBI,NORBI,NZ)
      IF (LOC%IPRLOC .GE. 10) THEN
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
      IF (D_MAX .GT. VAR%THRZER) THEN
         WRITE(LUPRI,*) 'ROTCMO_MOD: maximum deviation from unitary matrix'
         WRITE(LUPRI,*) 'Element ',J_MAX,I_MAX,' value :',D_MAX
      END IF

!     Calculate new orbitals: CMO_{old} exp(-kappa)
!     ---------------------------------------------

      IF (LOC%IPRLOC .GE. 10) THEN
         CALL HEADER('ROTCMO_MOD: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
!SK   CALL DZERO(RCOS,NORBXQ)
      CALL QGEMM(NBASI,NORBI,NORBI,D1,           &
                 'N','N',IP,CMO,NBASI,NORBI,NZ,  &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ, &
                 D0,IP,RCOS,NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,RCOS,1,CMO,1)
      IF (LOC%IPRLOC .GE. 5) THEN
         CALL HEADER('ROTCMO_MOD: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (LOC%IPRLOC .GE. 2) THEN
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_hessian(pmat,hmat,totdiag,totnond)
  !******************************************************************************!
  !                                                                              !
  !  Construction of full hessian                                                !
  !                                                                              ! 
  !     HMAT         A = d^2 G / dK* dK     B = d^2 G / dK dK          K         !
  !  / A   B  \       / A_uu   A_ub \        / B_uu   B_ub \      / K(unbared) \ !
  !  \ B*  A* /       \ A_bu   A_bb /        \ B_bu   B_bb /      \ K(bared)   / !
  !                                                                              ! 
  !******************************************************************************!
  use memory_allocator
  implicit none

  ! External variables
  real*8                     :: totdiag, totnond
  real*8, dimension(:,:,:)   :: hmat
  real*8, dimension(:,:,:,:) :: pmat

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d2 = 2.0d0, dhalf = 0.5d0, d4 = 4.0d0

  real*8 :: pij1, pji1, pkl1, pjl1, plj1, pki1, pik1, pil1, pli1, pkj1
  real*8 :: pij2, pji2, pkl2, pjl2, plj2, pki2, pik2, pil2, pli2, pkj2
  real*8 :: pii, pjj, pkk, pll, dtmp1, dtmp2, sumdiag, sumnond, sumtmp

  integer*4 :: i, j, i1, j1, k, l, m, ij, kl, nucdep, ni, ndimh, ni12, lupri

  real*8, allocatable, dimension(:,:) :: hmat_r

  !------------------------------------------------------------------------------!
  call qenter('get_hessian')

  lupri = var%lupri
  hmat(:,:,:) = d0
  totdiag = d0
  totnond = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimh  = size(hmat,1)
  ni12   = ndimh / 4

  !------------------------------------------------------------------------------!
  !                 Unbared-unbared contribution to the hessian                  !
  !------------------------------------------------------------------------------!
  !  d^2 G / dK*_ij dK_kl   -->   (32)
  sumdiag = d0
  sumnond = d0
  ij = 0
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,1,m); pij2 = pmat(i,j,2,m)
            pkl1 = pmat(k,l,1,m); pkl2 = pmat(k,l,2,m)
            pjl1 = pmat(j,l,1,m); pjl2 = pmat(j,l,2,m)
            pki1 = pmat(k,i,1,m); pki2 = pmat(k,i,2,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 + pij2*pkl2)
            dtmp2= d2*(pij1*pkl2 - pij2*pkl1)

            if(i==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * pjl1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 + d2 * ( pii - dhalf*(pjj + pll) ) * pjl2
            endif

            if(j==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki2
            endif

            if(j==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2
            endif

            if(i==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2
            endif
          enddo
          if(ij==kl)then
            sumdiag=sumdiag + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          else
            sumnond=sumnond + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          endif
        enddo
      enddo
    enddo
  enddo

  !  d^2 G / dK*_ij dK_kk   -->   (40)
  ij = 0
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,1,m); pij2 = pmat(i,j,2,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + dhalf * pij1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) - dhalf * pij2 * ( pii - pjj )
          endif

          if(j==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) - dhalf * pij1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) + dhalf * pij2 * ( pii - pjj )
          endif
        enddo
        hmat(kl,ij,1) =   hmat(ij,kl,1)
        hmat(kl,ij,2) = - hmat(ij,kl,2)
        sumnond=sumnond + d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,1)  -->  A_uu'
  do j=1,ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,1),k=1,ni12)
  enddo
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,2)  -->  A_uu'
  do j=1,ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,2),k=1,ni12)
  enddo
  write(lupri,*) 'Sum of diagonal vs non-diagonal elements of A_uu:'
  write(lupri,'(a,2f10.5)') 'Diag vs. non-diag: A_uu ',sumdiag,sumnond
  totdiag = totdiag + sumdiag
  totnond = totnond + sumnond

  !  d^2 G / dK_ij dK_kl   -->   (36)
  sumnond = d0
  sumtmp  = d0
  ij = 2*ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,1,m); pij2 = pmat(i,j,2,m)
            pkl1 = pmat(k,l,1,m); pkl2 = pmat(k,l,2,m)
            pil1 = pmat(i,l,1,m); pil2 = pmat(i,l,2,m)
            pkj1 = pmat(k,j,1,m); pkj2 = pmat(k,j,2,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 - pij2*pkl2)
            dtmp2= d2*(pij1*pkl2 + pij2*pkl1)

            if(i==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2
            endif

            if(j==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2
            endif

            if(j==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 - d2 * ( pjj - dhalf*(pii + pll) ) * pil1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 - d2 * ( pjj - dhalf*(pii + pll) ) * pil2
            endif

            if(i==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 - d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 - d2 * ( pii - dhalf*(pjj + pkk) ) * pkj2
            endif
          enddo
          sumnond=sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          if(ij==kl)then
            sumtmp = sumtmp + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          endif
        enddo
      enddo
    enddo
  enddo

  !  d^2 G / dK_ij dK_kk   -->   (44)
  ij = 2*ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,1,m); pij2 = pmat(i,j,2,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + dhalf * pij1 * ( pjj - pii )
            hmat(ij,kl,2) = hmat(ij,kl,2) + dhalf * pij2 * ( pjj - pii )
          endif

          if(j==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) - dhalf * pij1 * ( pjj - pii )
            hmat(ij,kl,2) = hmat(ij,kl,2) - dhalf * pij2 * ( pjj - pii )
          endif
        enddo
        hmat(kl+2*ni12,ij-2*ni12,1) = hmat(ij,kl,1)
        hmat(kl+2*ni12,ij-2*ni12,2) = hmat(ij,kl,2)
        sumnond=sumnond + d4*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,1)  -->  B_uu'
  do j=2*ni12+1, 3*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,1),k=1,ni12)
  enddo
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,2)  -->  B_uu'
  do j=2*ni12+1, 3*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,2),k=1,ni12)
  enddo
  write(lupri,*) 'Sum of elements of B_uu:'
  write(lupri,'(a,2f10.5)') 'Diag vs. non-diag: B_uu ',sumtmp,sumnond
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                   Bared-unbared contribution to the hessian                  !
  !------------------------------------------------------------------------------!
  !  d^2 G / dK*_ij dK_kl   -->   (34)
  sumnond = d0
  sumtmp  = d0
  ij = ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
            pkl1 = pmat(k,l,1,m); pkl2 = pmat(k,l,2,m)
            plj1 = pmat(l,j,3,m); plj2 = pmat(l,j,4,m)
            pli1 = pmat(l,i,3,m); pli2 = pmat(l,i,4,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 + pij2*pkl2)
            dtmp2= d2*(pij1*pkl2 - pij2*pkl1)

            if(i==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * plj1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 - d2 * ( pii - dhalf*(pjj + pll) ) * plj2
            endif

            if(j==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2
            endif

            if(j==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pli1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 - d2 * ( pjj - dhalf*(pii + pll) ) * pli2
            endif

            if(i==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2
            endif
          enddo
          sumnond=sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          if(ij==kl)then
            sumtmp = sumtmp + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          endif
        enddo
      enddo
    enddo
  enddo

  !  d^2 G / dK*_ij dK_kk   -->   (42)
  ij = ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + dhalf * pij1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) - dhalf * pij2 * ( pii - pjj )
          endif

          if(j==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + dhalf * pij1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) - dhalf * pij2 * ( pii - pjj )
          endif
        enddo
        sumnond = sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  !  d^2 G / dK*_kk dK_ij    -->   (41)^*
  ij = 0
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 + 1
      do k=1,ni
        do m = 1,nucdep
          pji1 = pmat(j,i,3,m); pji2 = pmat(j,i,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(kl,ij,1) = hmat(kl,ij,1) + pji1 * ( pii - pjj )
            hmat(kl,ij,2) = hmat(kl,ij,2) - pji2 * ( pii - pjj )
          endif
        enddo
        sumnond = sumnond+d2*dsqrt( hmat(kl,ij,1)*hmat(kl,ij,1) + hmat(kl,ij,2)*hmat(kl,ij,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,1)  -->  A_bu'
  do j=ni12+1,2*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,1),k=1,ni12)
  enddo
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,2)  -->  A_bu'
  do j=ni12+1,2*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,2),k=1,ni12)
  enddo
  write(lupri,*) 'Sum of diagonal vs non-diagonal elements of A_bu:'
  write(lupri,'(a,2f10.5)') 'Diag vs. non-diag: A_bu ',sumtmp,sumnond
  totnond = totnond + sumnond

  !  d^2 G / dK_ij dK_kl   -->   (38)
  sumnond = d0
  sumtmp  = d0
  ij = 3*ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
            pkl1 = pmat(k,l,1,m); pkl2 = pmat(k,l,2,m)
            pki1 = pmat(k,i,3,m); pki2 = pmat(k,i,4,m)
            pkj1 = pmat(k,j,3,m); pkj2 = pmat(k,j,4,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 - pij2*pkl2)
            dtmp2= d2*(pij1*pkl2 + pij2*pkl1)

            if(i==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2
            endif

            if(j==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 - d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 - d2 * ( pjj - dhalf*(pii + pkk) ) * pki2
            endif

            if(j==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2
            endif

            if(i==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 - d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 - d2 * ( pii - dhalf*(pjj + pkk) ) * pkj2
            endif
          enddo
          sumnond=sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          if(ij==kl)then
            sumtmp = sumtmp + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          endif
        enddo
      enddo
    enddo
  enddo

  !  d^2 G / dK_ij dK_kk   -->   (46)
  ij = 3*ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + dhalf * pij1 * ( pjj - pii )
            hmat(ij,kl,2) = hmat(ij,kl,2) + dhalf * pij2 * ( pjj - pii )
          endif

          if(j==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + dhalf * pij1 * ( pjj - pii )
            hmat(ij,kl,2) = hmat(ij,kl,2) + dhalf * pij2 * ( pjj - pii )
          endif
        enddo
        sumnond = sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  !  d^2 G / dK_kk dK_ij    -->   (45)
  ij = 0
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 3*ni12 + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(j==k)then
            hmat(kl,ij,1) = hmat(kl,ij,1) + pij1 * ( pii - pjj )
            hmat(kl,ij,2) = hmat(kl,ij,2) + pij2 * ( pii - pjj )
          endif
        enddo
        sumnond = sumnond+d2*dsqrt( hmat(kl,ij,1)*hmat(kl,ij,1) + hmat(kl,ij,2)*hmat(kl,ij,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,1)  -->  B_bu'
  do j=3*ni12+1, 4*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,1),k=1,ni12)
  enddo
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,2)  -->  B_bu'
  do j=3*ni12+1, 4*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,2),k=1,ni12)
  enddo
  write(lupri,*) 'Sum of elements of B_bu:'
  write(lupri,'(a,2f10.5)') 'Diag vs. non-diag: B_bu ',sumtmp,sumnond
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                   Bared-bared contribution to the hessian                    !
  !------------------------------------------------------------------------------!
  !  d^2 G / dK*_ij dK_kl   -->   (35)
  sumdiag = d0
  sumnond = d0
  ij = ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
            pkl1 = pmat(k,l,3,m); pkl2 = pmat(k,l,4,m)
            pki1 = pmat(k,i,1,m); pki2 = pmat(k,i,2,m)
            pkj1 = pmat(k,j,1,m); pkj2 = pmat(k,j,2,m)
            pli1 = pmat(l,i,1,m); pli2 = pmat(l,i,2,m)
            plj1 = pmat(l,j,1,m); plj2 = pmat(l,j,2,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 + pij2*pkl2)
            dtmp2= d2*(pij1*pkl2 - pij2*pkl1)

            if(i==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * plj1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 + d2 * ( pii - dhalf*(pjj + pll) ) * plj2
            endif

            if(j==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki2
            endif

            if(j==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pli1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 + d2 * ( pjj - dhalf*(pii + pll) ) * pli2
            endif

            if(i==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj2
            endif
          enddo
        if(ij==kl)then
          sumdiag = sumdiag + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        else
          sumnond = sumnond + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        endif
        enddo
      enddo
    enddo
  enddo

  !  d^2 G / dK*_ij dK_kk   -->   (43)
  ij = ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,1,m); pij2 = pmat(i,j,2,m)
          pji1 = pmat(j,i,1,m); pji2 = pmat(j,i,2,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + pij1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) + pij2 * ( pii - pjj )
          endif

          if(j==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + pji1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) + pji2 * ( pii - pjj )
          endif
        enddo
        hmat(kl,ij,1) =   hmat(ij,kl,1)
        hmat(kl,ij,2) = - hmat(ij,kl,2)
        sumnond = sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,1)  -->  A_bb'
  do j=ni12+1,2*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,1),k=ni12+1,2*ni12)
  enddo
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,2)  -->  A_bb'
  do j=ni12+1,2*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,2),k=ni12+1,2*ni12)
  enddo
  write(lupri,*) 'Sum of diagonal vs non-diagonal elements of A_bb:'
  write(lupri,'(a,2f10.5)') 'Diag vs. non-diag: A_bb ',sumdiag,sumnond
  totdiag = totdiag + sumdiag
  totnond = totnond + sumnond

  !  d^2 G / dK_ij dK_kl = 0   -->   (49)
  !  d^2 G / dK_ij dK_kk = 0   -->   (47)

  !------------------------------------------------------------------------------!
  !                 Unbared-bared contribution to the hessian                    !
  !------------------------------------------------------------------------------!
  !  d^2 G / dK_ij dK_kl   -->   (37)
  sumnond = d0
  sumtmp  = d0
  ij = 2*ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,1,m); pij2 = pmat(i,j,2,m)
            pkl1 = pmat(k,l,3,m); pkl2 = pmat(k,l,4,m)
            pik1 = pmat(i,k,3,m); pik2 = pmat(i,k,4,m)
            pil1 = pmat(i,l,3,m); pil2 = pmat(i,l,4,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 - pij2*pkl2)
            dtmp2= d2*(pij1*pkl2 + pij2*pkl1)

            if(i==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2
            endif

            if(j==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) + dtmp1 - d2 * ( pjj - dhalf*(pii + pkk) ) * pik1
              hmat(ij,kl,2) = hmat(ij,kl,2) + dtmp2 - d2 * ( pjj - dhalf*(pii + pkk) ) * pik2
            endif

            if(j==k)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1 - d2 * ( pjj - dhalf*(pii + pll) ) * pil1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2 - d2 * ( pjj - dhalf*(pii + pll) ) * pil2
            endif

            if(i==l)then
              hmat(ij,kl,1) = hmat(ij,kl,1) - dtmp1
              hmat(ij,kl,2) = hmat(ij,kl,2) - dtmp2
            endif
          enddo
          sumnond=sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          if(ij==kl)then
            sumtmp = sumtmp + dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
          endif
        enddo
      enddo
    enddo
  enddo

  !  d^2 G / dK_ij dK_kk   -->   (45)
  ij = 2*ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(j==k)then
            hmat(ij,kl,1) = hmat(ij,kl,1) + pij1 * ( pii - pjj )
            hmat(ij,kl,2) = hmat(ij,kl,2) + pij2 * ( pii - pjj )
          endif
        enddo
        sumnond = sumnond+d2*dsqrt( hmat(ij,kl,1)*hmat(ij,kl,1) + hmat(ij,kl,2)*hmat(ij,kl,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  !  d^2 G / dK_kk dK_ij    -->   (46)
  ij = ni12
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 2*ni12 + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k)then
            hmat(kl,ij,1) = hmat(kl,ij,1) + dhalf * pij1 * ( pjj - pii )
            hmat(kl,ij,2) = hmat(kl,ij,2) + dhalf * pij2 * ( pjj - pii )
          endif

          if(j==k)then
            hmat(kl,ij,1) = hmat(kl,ij,1) + dhalf * pij1 * ( pjj - pii )
            hmat(kl,ij,2) = hmat(kl,ij,2) + dhalf * pij2 * ( pjj - pii )
          endif
        enddo
        sumnond = sumnond+d2*dsqrt( hmat(kl,ij,1)*hmat(kl,ij,1) + hmat(kl,ij,2)*hmat(kl,ij,2) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,1)  -->  B_ub'
  do j=2*ni12+1, 3*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,1),k=ni12+1,2*ni12)
  enddo
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print hmat(:,:,2)  -->  B_ub'
  do j=2*ni12+1, 3*ni12
    write(lupri,'(i4,15f10.5)') j,(hmat(j,k,2),k=ni12+1,2*ni12)
  enddo
  write(lupri,*) 'Sum of elements of B_ub:'
  write(lupri,'(a,2f10.5)') 'Diag vs. non-diag: B_ub ',sumtmp,sumnond
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                     Get the rest of the hessian matrix                       !
  !------------------------------------------------------------------------------!
  ! 2,1  -->  1,2
  do i=ni12+1,2*ni12
    do j=1,ni12
      hmat(j,i,1) =   hmat(i,j,1)
      hmat(j,i,2) = - hmat(i,j,2)
    enddo
  enddo

  ! 12,12  -->  34,34
  do i=1,2*ni12
    do j=1,2*ni12
      hmat( i+2*ni12 , j+2*ni12 , 1 ) =   hmat(i,j,1)
      hmat( i+2*ni12 , j+2*ni12 , 2 ) = - hmat(i,j,2)
    enddo
  enddo

  ! 34,12  -->  12,34
  do i=2*ni12+1,4*ni12
    do j=1,2*ni12
      hmat(j,i,1) =   hmat(i,j,1)
      hmat(j,i,2) = - hmat(i,j,2)
    enddo
  enddo

  ! Calculate real representation of hessian and get the diagonal dominance info
  call alloc(hmat_r,ndimh,ndimh,id='hmat_r in get_hessian')

  sumdiag = d0
  sumnond = d0
  do i=1,2*ni12
    do j=1,2*ni12
      i1 = i + 2*ni12
      j1 = j + 2*ni12
      hmat_r(i ,j ) = hmat(i,j ,1) + hmat(i,j1,1)
      hmat_r(i ,j1) = hmat(i,j1,2) - hmat(i,j ,2)
      hmat_r(i1,j ) = hmat(i,j1,2) + hmat(i,j ,2)
      hmat_r(i1,j1) = hmat(i,j ,1) - hmat(i,j1,1)
      sumnond = sumnond + dabs(hmat_r(i,j1)) + dabs(hmat_r(i1,j))
      if(i==j)then
        sumdiag = sumdiag + dabs(hmat_r(i,j)) + dabs(hmat_r(i1,j1))
      else
        sumnond = sumnond + dabs(hmat_r(i,j)) + dabs(hmat_r(i1,j1))
      endif
    enddo
  enddo
  write(6,'(a,5d11.2)')' Real hessian diag vs. non-diag: ',sumdiag,sumnond

  call dealloc(hmat_r,id='hmat_r in get_hessian')

  call qexit('get_hessian')

  end subroutine get_hessian

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_gradient(pmat,gradient,dnrm)
  !******************************************************************************!
  !                                                                              !
  !  Construction of gradient                                                    !
  !                                                                              ! 
  !******************************************************************************!
  implicit none

  ! External variables
  real*8 :: dnrm
  real*8, dimension(:,:)     ::  gradient
  real*8, dimension(:,:,:,:) ::  pmat

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d2 = 2.0d0

  real*8 :: pij1, pij2, dtmp

  integer*4 :: i, j, k, ij, nucdep, ni, ndimg, ni12, lupri

  !------------------------------------------------------------------------------!
  call qenter('get_gradient')

  lupri = var%lupri
  gradient(:,:) = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimg  = size(gradient,1)
  ni12   = ndimg / 4

  !------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------!
  ij = 0
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      do k = 1,nucdep
        dtmp = d2*(pmat(j,j,1,k) - pmat(i,i,1,k))
        gradient(ij,1)      = gradient(ij,1)      + dtmp * pmat(i,j,1,k)
        gradient(ij,2)      = gradient(ij,2)      + dtmp * pmat(i,j,2,k)
        gradient(ij+ni12,1) = gradient(ij+ni12,1) + dtmp * pmat(i,j,3,k)
        gradient(ij+ni12,2) = gradient(ij+ni12,2) + dtmp * pmat(i,j,4,k)
      enddo
      gradient(ij+2*ni12,1) =   gradient(ij,1)
      gradient(ij+2*ni12,2) = - gradient(ij,2)
      gradient(ij+3*ni12,1) =   gradient(ij+ni12,1)
      gradient(ij+3*ni12,2) = - gradient(ij+ni12,2)
    enddo
  enddo

  write(6,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,1)  -->  gradient'
  write(lupri,'(9d12.4)') (gradient(j,1),j=1,ni12)
  write(6,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,2)  -->  gradient'
  write(lupri,'(9d12.4)') (gradient(j,2),j=1,ni12)
  write(6,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,1)  -->  /bar{gradient}'
  write(lupri,'(9d12.4)') (gradient(j,1),j=1+ni12,2*ni12)
  write(6,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,2)  -->  /bar{gradient}'
  write(lupri,'(9d12.4)') (gradient(j,2),j=1+ni12,2*ni12)

  dnrm = d0
  do i=1,2*ni12
    dnrm = dnrm + gradient(i,1)*gradient(i,1) + gradient(i,2)*gradient(i,2)
  enddo
  dnrm = dsqrt(dnrm)

  call qexit('get_gradient')

  end subroutine get_gradient

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine errall(ier,chsub,charr)
  implicit none

  integer*4    :: ier
  character(*) :: chsub, charr

  if(ier /= 0)then
    write(6,*)
    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
    write(6,'(a,i4)'   )' iostat =',ier
    stop 'ERROR in allocation (see output form more information)'
  endif

  end subroutine errall

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine exp_func(xkappa,hmat,gradient,full,dinvmq)
  use memory_allocator
  implicit none

  ! External variables
  logical :: full
  real*8  :: dinvmq

  real*8,    allocatable, dimension(:,:)     :: gradient
  real*8,    allocatable, dimension(:,:,:)   :: xkappa, hmat

  ! Internal variables
  integer*4  :: i, j, i1, j1, ij, id, ni, ni12, ndimh, ii, jj
  real*8, parameter :: d0 = 0.0d0, d2 = 2.0d0

  real*8, allocatable, dimension(:)   :: newkappa, grad_r
  real*8, allocatable, dimension(:,:) :: hmat_r
  complex*16 :: ch, cki, ckj

  !------------------------------------------------------------------------------!
  call qenter('exp_func')

  ! Get some dimensions
  ni     = size(xkappa,1)
  ndimh  = size(hmat,1)
  ni12   = ndimh / 4

  call alloc(newkappa,ndimh,id='newkappa in exp_func')
  call alloc(hmat_r,ndimh,ndimh,id='hmat_r in exp_func')
  call alloc(grad_r,ndimh,      id='grad_r in exp_func')

  newkappa(:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      newkappa(ij+0*ni12) =   xkappa(j,i,1)
      newkappa(ij+1*ni12) = - xkappa(j,i,3)
      newkappa(ij+2*ni12) =   xkappa(j,i,2)
      newkappa(ij+3*ni12) =   xkappa(j,i,4)
    enddo
  enddo

  ! prepare hmat_r
  hmat_r(:,:) = d0
  do i=1,2*ni12
    do j=1,2*ni12
      if((.not.full).and.(i/=j))cycle
      i1 = i + 2*ni12
      j1 = j + 2*ni12
      hmat_r(i ,j ) = hmat(i,j ,1) + hmat(i,j1,1)
      hmat_r(i ,j1) = hmat(i,j1,2) - hmat(i,j ,2)
      hmat_r(i1,j ) = hmat(i,j1,2) + hmat(i,j ,2)
      hmat_r(i1,j1) = hmat(i,j ,1) - hmat(i,j1,1)
    enddo
  enddo

  ! contribution from hessian
  dinvmq = d0
  do i=1, ndimh
    do j=1, ndimh
      dinvmq = dinvmq + newkappa(i) * hmat_r(i,j) * newkappa(j)
    enddo
  enddo

  ! prepare grad_r
  do i=1,2*ni12
    i1 = i + 2*ni12
    grad_r(i ) =   gradient(i,1)
    grad_r(i1) = - gradient(i,2)
  enddo

  ! contribution from gradient
  do i=1,ndimh
    dinvmq = dinvmq + d2*newkappa(i) * grad_r(i)
  enddo
  dinvmq = dinvmq / dfloat(ni)

  ! complex multiplication, better use real one !!
! dinvmq = d0
! do i=1,2*ni12
!   ch = dcmplx(gradient(i,1),gradient(i,2))
!   cki= dcmplx(newkappa(i),newkappa(i+2*ni12))
!   dinvmq = dinvmq + dreal(ch * cki)
!
!   i1 = i + 2*ni12
!   ch = dcmplx(gradient(i1,1),gradient(i1,2))
!   cki= dcmplx(newkappa(i),-newkappa(i+2*ni12))
!   dinvmq = dinvmq + dreal(ch * cki)
! enddo
! write(6,*) 'dinvmq g =',dinvmq

  ! complex multiplication, better use real one !!
! dsum = 0.0d0
! ii = 0
! do i=1, ndimh
!   ii = ii + 1
!   jj = 0
!   do j=1, ndimh
!     jj = jj + 1
!     ch  = dcmplx(hmat(i,j,1),hmat(i,j,2))
!     if(i<=ni12)                    cki = dcmplx(newkappa(ii+0*ni12),  newkappa(ii+2*ni12))
!     if((i>1*ni12).and.(i<=2*ni12)) cki = dcmplx(newkappa(ii+1*ni12),  newkappa(ii+3*ni12))
!     if((i>2*ni12).and.(i<=3*ni12)) cki = dcmplx(newkappa(ii+0*ni12),- newkappa(ii+2*ni12))
!     if((i>3*ni12).and.(i<=4*ni12)) cki = dcmplx(newkappa(ii+1*ni12),- newkappa(ii+3*ni12))
!     if(j<=ni12)                    ckj = dcmplx(newkappa(jj+0*ni12),  newkappa(jj+2*ni12))
!     if((j>1*ni12).and.(j<=2*ni12)) ckj = dcmplx(newkappa(jj+1*ni12),  newkappa(jj+3*ni12))
!     if((j>2*ni12).and.(j<=3*ni12)) ckj = dcmplx(newkappa(jj+0*ni12),- newkappa(jj+2*ni12))
!     if((j>3*ni12).and.(j<=4*ni12)) ckj = dcmplx(newkappa(jj+1*ni12),- newkappa(jj+3*ni12))
!     dsum = dsum + dreal(dconjg(cki) * ch * ckj)
!     if(jj == ni12)jj = 0
!   enddo
!   if(ii == ni12)ii = 0
! enddo
! write(6,*) 'dsum',dsum

  !------------------------------------------------------------------------------!
  call dealloc(newkappa,id='newkappa in exp_func')
  call dealloc(hmat_r  ,id='hmat_r   in exp_func')
  call dealloc(grad_r  ,id='grad_r   in exp_func')

  call qexit('exp_func')

  !------------------------------------------------------------------------------!

  end subroutine exp_func

end module localization
