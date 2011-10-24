module localization

use memory_allocator
implicit none

public :: ORBLOC, initialize_localization

private :: ORBLP1, ORBLP2, ORBLO1, ORBLO4, ORBLO5, ROTCMO_MOD
private :: getk_cdiag, getk_rdiag, getk_full, get_hessian, get_gradient, exp_func, hdiag
private :: get_mu, updateT

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
  character*6  :: reffil(15)    ! reffil(maxref)    -->  see dcbprj.h
  character*72 :: vecref(2,15)  ! vecref(2,maxref)  -->  see dcbprj.h
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
  implicit none

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0

  character*12 :: sectid, cputid, walltid
  character    :: daytid*24, hessian*4 = 'DIAG', conv*20 = 'No convergence after'
  logical      :: tobe, full, diagonly, rejected

  integer*4 :: ifrp, ndimh, ndimhr, idum, itr, lupri, nish, norb, nucdep, maxref
  integer*4 :: nfbas(0:2), ipqtoq(4,0:7)

  real*8    :: dinvmo, dinvmg, dinvmq, dnrmk, dnrmg, diffg, diffq, diffa
  real*8    :: timend, time, timstr, second, toterg, totdiag, totnond, shift, tradius, tpred

  integer*4, allocatable, dimension(:)       :: nfnuc, isel, indao, indaos, indnao
  real*8,    allocatable, dimension(:)       :: scoef, storec, smat, cmo, grad, eval, tgrad
  real*8,    allocatable, dimension(:,:)     :: grad_c, hmat
  real*8,    allocatable, dimension(:,:,:)   :: xkappa, hmat_c
  real*8,    allocatable, dimension(:,:,:,:) :: pmat

  real*8, pointer :: tgradp(:), evalp(:), evecp(:,:)

  !------------------------------------------------------------------------------!
  call qenter('orbloc')

  ! Set some parameters
  ifrp   = 1
  lupri  = var%lupri
  nish   = var%nish(ifrp)
  norb   = var%norb(ifrp)
  nucdep = var%nucdep
  maxref = var%maxref
  full   = loc%auto

  nfbas(:)    = var%nfbas(ifrp,:)
  ipqtoq(:,:) = var%ipqtoq(:,:)

  diagonly = .false.
  rejected = .false.
  if(loc%hesthr < loc%thrloc) diagonly = .true.

  tradius = 10.0d0

  ! TMP print
  write(lupri,*) 'lupri =',lupri
  write(lupri,*) 'nish  =',nish
  write(lupri,*) 'norb  =',norb
  write(lupri,*) 'nucdep=',nucdep
  write(lupri,*) 'maxref=',maxref
  write(lupri,*) 'full  =',full
  write(lupri,*) 'hesthr=',loc%hesthr
  write(lupri,*) 'thrloc=',loc%thrloc

  call titler('Molecular orbital localization','*',122)

  ! Check if coefficients are on file
  inquire(file='DFCOEF',exist=tobe)
  if(.not.tobe)then
    write(lupri,'(a)') 'orbloc: Coefficient file not found !'
    call quit('orbloc: Coefficients not found !')
  endif

  ! Memory allocation
  ndimh = 4 * ( nish*nish + nish ) / 2
  ndimhr= ndimh - nish
  call alloc(scoef ,nish*nfbas(0)*4,   id='scoef  in orbloc')
  call alloc(storec,nish*nfbas(0)*4,   id='storec in orbloc')
  call alloc(smat  ,var%n2bbasx,       id='smat   in orbloc')
  call alloc(cmo   ,var%n2bbasxq,      id='cmo    in orbloc')
  call alloc(hmat_c,ndimh,ndimh,2,     id='hmat_c in orbloc')
  call alloc(hmat  ,ndimhr,ndimhr,     id='hmat   in orbloc')
  call alloc(grad_c,ndimh,2,           id='grad_c in orbloc')
  call alloc(grad  ,ndimhr,            id='grad   in orbloc')
  call alloc(eval  ,ndimhr,            id='eval   in orbloc')
  call alloc(tgrad ,ndimhr,            id='tgrad  in orbloc')
  call alloc(xkappa,nish,nish,4,       id='xkappa in orbloc')
  call alloc(pmat  ,nish,nish,4,nucdep,id='pmat   in orbloc')

  call alloc(nfnuc ,2*nucdep,  id='nfnuc  in orbloc')
  call alloc(isel  ,nish,      id='isel   in orbloc')
  call alloc(indao ,nfbas(0),  id='indao  in orbloc')
  call alloc(indaos,nfbas(0),  id='indaos in orbloc')
  call alloc(indnao,2*nucdep+2,id='indnao in orbloc')

  ! nullify pointers for later allocation using memory_allocator module
  nullify(tgradp)
  nullify(evalp )
  nullify(evecp )

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
  write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 4x,a, 5x,a, 6x,a)') &
    '*T','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'

  ! SCF cycle
  dinvmq = d0
  itr    = 0
  do

    ! get new kappa
    if(.not.rejected)then
      call get_hessian(pmat,hmat_c,hmat,totdiag,totnond)
      call get_gradient(pmat,grad_c,grad,dnrmg)
      call hdiag(hmat,grad,eval,tgrad,tgradp,evalp,evecp,itr)
      shift = d0
      call deallocp(tgradp,id='tgradp in orbloc')
      call deallocp(evalp ,id='evalp  in orbloc')
      call deallocp(evecp ,id='evecp  in orbloc')
    else
      call get_mu(eval,tgrad,tradius,shift,tpred)
    endif
    if(     full) call getk_full(pmat,xkappa,hmat,grad,shift,dnrmk)
    if(.not.full) call getk_rdiag(pmat,xkappa,grad_c,dnrmk,ipqtoq)
    if( (full.and.rejected) .and. (dabs(tpred-dnrmk) > 1.0d-2*tradius) )then
      write(lupri,'(a)') ' WARNING: Predicted and actual norm of kappa are too different'
      write(lupri,'(a,i4,a,d12.4,a,d12.4,a,2d12.4)')&
        ' WARNING: itr:',itr,' Predicted:',tpred,' Actual:',dnrmk,' Trust-radius:',tradius
    endif
    call exp_func(xkappa,hmat,grad,full,dinvmq)

    if(.not.rejected) storec(:) = scoef(:)  ! store MO coeff. in the case you reject the step
    if(     rejected) scoef(:)  = storec(:) ! restore MO coeff.; previous step was rejected
    call rotcmo_mod(scoef,nfbas(0),nish,4,xkappa,nish,ipqtoq(ifrp,0))
    if(loc%prjloc) then
      ! PRJLOC = TRUE localization with projection analysis
      call gtovlx(smat,var%ssmtrc)
      call orblp1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,maxref,dinvmg)
    else
      ! PRJLOC = FALSE localization with Mulliken analysis
      call orblo1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,nfnuc,dinvmg)
    endif

    ! Update Trust-radius
    diffg = dinvmg - dinvmo
    diffq = dinvmq
    if(full) call updateT(dnrmk,diffg,diffq,tradius,rejected,itr)

    if(.not.rejected)then
      itr = itr + 1
      call gtinfo(daytid)
      write(lupri,'(/5x,a,i4,a,3x)')'########## START ITERATION NO. ',itr,' ##########'
      write(lupri,'(5x,a,a24,a)')'########## ',daytid,' ##########'

      if (diffg < d0) write(lupri,'(a,i4)')' WARNING: dinvm_new - dinvm_old < 0; itr =',itr
      diffa = dabs(dinvmg - dinvmo)
      if(full) hessian = 'FULL'
      dinvmo = dinvmg

      write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 5x,a, 6x,a, 6x,a)') &
        '  ','ITR','DIFFA','|GRAD|','|KAPPA|','H_DIAG','H_NONDIAG','H','DIFFG','DIFFQ'
      write(lupri,'(a,2x,i4,5d11.2,4x,a,4d11.2)') &
        '*C',itr,diffa,dnrmg,dnrmk,totdiag,totnond,hessian,diffg,diffq

      ! Exit conditions
      if((diffa < loc%thrloc).and.(full.or.diagonly)) exit
      if( diffa < loc%hesthr ) full = .true.
      if(itr > loc%itrloc)exit
    endif

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
  call dealloc(scoef ,id='scoef  in orbloc')
  call dealloc(storec,id='storec in orbloc')
  call dealloc(smat  ,id='smat   in orbloc')
  call dealloc(pmat  ,id='pmat   in orbloc')
  call dealloc(xkappa,id='xkappa in orbloc')
  call dealloc(hmat_c,id='hmat_c in orbloc')
  call dealloc(hmat  ,id='hmat   in orbloc')
  call dealloc(grad_c,id='grad_c in orbloc')
  call dealloc(grad  ,id='grad   in orbloc')
  call dealloc(eval  ,id='eval   in orbloc')
  call dealloc(tgrad ,id='tgrad  in orbloc')
  call dealloc(cmo   ,id='cmo    in orbloc')

  call dealloc(nfnuc ,id='nfnuc  in orbloc')
  call dealloc(isel  ,id='isel   in orbloc')
  call dealloc(indao ,id='indao  in orbloc')
  call dealloc(indaos,id='indaos in orbloc')
  call dealloc(indnao,id='indnao in orbloc')

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
      INTEGER*4, INTENT(IN ) :: NBL, NBS, NI, NUCDEP, MAXREF
      REAL*8,    INTENT(OUT) :: DINVM

      REAL*8,    INTENT(IN   ) :: SCOEF(*)
      REAL*8,    INTENT(INOUT) :: SMAT(*)
      REAL*8,    INTENT(  OUT) :: PMAT(NI,NI,4,NUCDEP)

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
      CALL ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NI,NVECS(IFRP,1),BV,NSTR,IPIV,MAXREF)

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
      WRITE(LUPRI,'(/A,F20.16)') '* Inverse mean delocalization = ',DINVM

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

 SUBROUTINE ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NMOL,NREF,BVEC,NSTR,IPIVOT,MAXREF)
!***********************************************************************
!
!     Projection analysis for fermion ircop IFRP
!
!     Written by T.Saue March 26 1999
!
!***********************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER*4, INTENT(IN ) :: NB, NMOL, NREF, MAXREF
      INTEGER*4, INTENT(IN ) :: NSTR(2,0:2,0:MAXREF)
      INTEGER*4, INTENT(OUT) :: IPIVOT(*)

      REAL*8,    INTENT(IN   ) :: SCOEF(*), RCOEF(*)
      REAL*8,    INTENT(INOUT) :: SMAT(*)
      REAL*8,    INTENT(  OUT) :: AMAT(NREF,NREF,4), BVEC(NREF,NMOL,4)

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
      INTEGER*4, INTENT(IN ) :: NBL, NBS, NI, NUCDEP
      REAL*8,    INTENT(OUT) :: DINVM

      INTEGER*4, INTENT(IN ) :: NFNUC(*)
      REAL*8,    INTENT(IN ) :: SMAT(*), SCOEF(*)
      REAL*8,    INTENT(OUT) :: PMAT(NI,NI,4,NUCDEP)

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
      DINVM = DINV/DFLOAT(NI)
      WRITE(LUPRI,'(/A,F20.16)') '* Inverse mean delocalization = ',DINVM

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
      CALL DEALLOC(WLMAT,ID='WLMAT IN ORBLO1')
      CALL DEALLOC(WSMAT,ID='WSMAT IN ORBLO1')

      CALL QEXIT('ORBLO1')

  END SUBROUTINE ORBLO1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine getk_cdiag(pmat,xkappa,gradient,dnrmk,ipqtoq)
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
  implicit none

  ! External variables
  integer*4, intent(in), dimension(4,0:7) :: ipqtoq

  real*8, intent(out) :: dnrmk

  real*8, intent(in ), dimension(:,:)     :: gradient
  real*8, intent(in ), dimension(:,:,:,:) :: pmat
  real*8, intent(out), dimension(:,:,:)   :: xkappa

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

  integer*4 :: i, j, k, ij, ij1, iz, nucdep, ni, ni12, lupri
  real*8    :: pii, pjj, pn2

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
  dnrmk = d0
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
      dnrmk = dnrmk + xkappa(j,i,1)*xkappa(j,i,1) + xkappa(j,i,2)*xkappa(j,i,2) &
                    + xkappa(j,i,3)*xkappa(j,i,3) + xkappa(j,i,4)*xkappa(j,i,4)
    enddo
  enddo

  ! Norm of the kappa vector
  dnrmk = dsqrt(dnrmk)
  write(6,*)
  write (lupri,'(a,e12.2)') ' Norm of the kappa: ',dnrmk

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

  call dealloc(hmat_d,id='hmat_d in getk_cdiag')

  call qexit('getk_cdiag')

  end subroutine getk_cdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine getk_rdiag(pmat,xkappa,gradient,dnrmk,ipqtoq)
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
  implicit none

  ! External variables
  integer*4, intent(in), dimension(4,0:7) :: ipqtoq

  real*8, intent(out) :: dnrmk

  real*8, intent(in ), dimension(:,:)     :: gradient
  real*8, intent(in ), dimension(:,:,:,:) :: pmat
  real*8, intent(out), dimension(:,:,:)   :: xkappa

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

  integer*4 :: i, j, k, ij, ij1, ij2, ij3, iz, nucdep, ni, ni12, lupri
  real*8    :: pii, pjj, pn2

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
  dnrmk = d0
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
      dnrmk = dnrmk + xkappa(j,i,1)*xkappa(j,i,1) + xkappa(j,i,2)*xkappa(j,i,2) &
                    + xkappa(j,i,3)*xkappa(j,i,3) + xkappa(j,i,4)*xkappa(j,i,4)
    enddo
  enddo

  ! Norm of the kappa vector
  dnrmk = dsqrt(dnrmk)
  write(6,*)
  write (lupri,'(a,e12.2)') ' Norm of the kappa: ',dnrmk

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

  call dealloc(hmat_d,id='hmat_d in getk_rdiag')

  call qexit('getk_rdiag')

  end subroutine getk_rdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine getk_full(pmat,xkappa,hmat,gradient,shift,dnrmk)
  !******************************************************************************!
  !                                                                              !
  !           Construction of the XKAPPA matrix using full hessian               !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  real*8, intent(in ) :: shift
  real*8, intent(out) :: dnrmk

  real*8, intent(in ), dimension(:)       :: gradient
  real*8, intent(in ), dimension(:,:)     :: hmat
  real*8, intent(in ), dimension(:,:,:,:) :: pmat
  real*8, intent(out), dimension(:,:,:)   :: xkappa

  ! Internal variables
  real*8, parameter     :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0
  complex*16, parameter :: c2 = (2.0d0,0.0d0)

  integer*4 :: i, j, k, i1, j1, ij, iz, nucdep, ni, nhtmp, ni12, lupri, ier
  real*8    :: pii, pjj, pn2, dnrm2, dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6

  integer*4,  allocatable, dimension(:)   :: ipiv
  real*8,     allocatable, dimension(:)   :: gtmp
  real*8,     allocatable, dimension(:,:) :: htmp

  !------------------------------------------------------------------------------!
  call qenter('getk_full')

  lupri  = var%lupri

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ni12   = (size(hmat,1)+ni) / 4
  nhtmp  = size(hmat,1)

  call alloc(ipiv,nhtmp      ,id='ipiv in getk_full')
  call alloc(gtmp,nhtmp      ,id='gtmp in getk_full')
  call alloc(htmp,nhtmp,nhtmp,id='htmp in getk_full')

  ! prepare hessian for dgetrf
  htmp(1:nhtmp,1:nhtmp) = hmat(1:nhtmp,1:nhtmp)
  do i=1,nhtmp
    htmp(i,i) = htmp(i,i) - shift
  enddo

  ! prepare gradient for dgetrf
  gtmp(1:nhtmp) = - gradient(1:nhtmp)

  ! Solve linear set of equations
  call dgetrf(nhtmp,nhtmp,htmp,nhtmp,ipiv,ier)
  if(ier /= 0)then
    write(lupri,*) ' ERROR: The dgetrf LAPACK subroutine in getk_full finished unsuccessfully'
    write(lupri,*) ' info =',ier
    stop
  endif

  call dgetrs('N',nhtmp,1,htmp,nhtmp,ipiv,gtmp,nhtmp,ier)
  if(ier /= 0)then
    write(lupri,*) ' ERROR: The dgetrs LAPACK subroutine in getk_full finished unsuccessfully'
    write(lupri,*) ' info =',ier
    stop
  endif

  ! Norm of the kappa
  dnrmk  = dnrm2(nhtmp,gtmp,1)
  write(6,*)
  write (lupri,'(a,e12.2)') ' Norm of the kappa: ',dnrmk

  xkappa(:,:,:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      if(i==j)cycle
      ij = ij + 1
      xkappa(j,i,1) =   gtmp(ij)
      xkappa(i,j,1) = - xkappa(j,i,1)
    enddo
  enddo

  ij = ni12-ni
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      xkappa(j,i,2) =   gtmp(ni12+ij)
      xkappa(i,j,2) =   xkappa(j,i,2)

      xkappa(j,i,3) = - gtmp(ij)
      xkappa(i,j,3) =   xkappa(j,i,3)

      xkappa(j,i,4) =   gtmp(2*ni12+ij)
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

  call dealloc(ipiv,id='ipiv in getk_full')
  call dealloc(gtmp,id='gtmp in getk_full')
  call dealloc(htmp,id='htmp in getk_full')

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
      INTEGER*4, INTENT(IN ) :: NUCDEP

      INTEGER*4, INTENT(IN ) :: INDNAO(2,*)
      INTEGER*4, INTENT(OUT) :: NFNUC(2*NUCDEP)

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
      INTEGER*4, INTENT(IN   ) :: NBASI, NORBI, NZ, NORBT, IP(4)
      REAL*8,    INTENT(IN   ) :: XKAPPA(NORBT,NORBT,NZ)
      REAL*8,    INTENT(INOUT) :: CMO(NBASI,NORBI,NZ)

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

  subroutine get_hessian(pmat,hmat,hmat_r,totdiag,totnond)
  !******************************************************************************!
  !                                                                              !
  !  Construction of full hessian                                                !
  !                                                                              ! 
  !     HMAT         A = d^2 G / dK* dK     B = d^2 G / dK dK          K         !
  !  / A   B  \       / A_uu   A_ub \        / B_uu   B_ub \      / K(unbared) \ !
  !  \ B*  A* /       \ A_bu   A_bb /        \ B_bu   B_bb /      \ K(bared)   / !
  !                                                                              ! 
  !******************************************************************************!
  implicit none

  ! External variables
  real*8, intent(out) :: totdiag, totnond

  real*8, intent(in ), dimension(:,:,:,:) :: pmat
  real*8, intent(out), dimension(:,:)     :: hmat_r
  real*8, intent(out), dimension(:,:,:)   :: hmat

  ! Internal variables
  real*8, parameter :: d0 = 0.0d0, d2 = 2.0d0, dhalf = 0.5d0, d4 = 4.0d0

  real*8 :: pij1, pji1, pkl1, pjl1, plj1, pki1, pik1, pil1, pli1, pkj1
  real*8 :: pij2, pji2, pkl2, pjl2, plj2, pki2, pik2, pil2, pli2, pkj2
  real*8 :: pii, pjj, pkk, pll, sumdiag, sumnond, sumtmp
  real*8 :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6

  integer*4 :: i, j, i1, j1, k, l, m, ij, kl, nucdep, ni, ndimhr, ndimhc, ni12, lupri

  real*8, allocatable, dimension(:,:) :: htmp

  !------------------------------------------------------------------------------!
  call qenter('get_hessian')

  lupri = var%lupri
  hmat(:,:,:) = d0
  totdiag = d0
  totnond = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimhc = size(hmat,1)
  ndimhr = size(hmat_r,1)
  ni12   = ndimhc / 4

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

  !------------------------------------------------------------------------------!
  ! Calculate real representation of hessian and get the diagonal dominance info !
  !------------------------------------------------------------------------------!
  call alloc(htmp,ndimhc,ndimhc,id='htmp in get_hessian')

  sumdiag = d0
  sumnond = d0
  do i=1,2*ni12
    do j=1,2*ni12
      i1 = i + 2*ni12
      j1 = j + 2*ni12
      htmp(i ,j ) = hmat(i,j ,1) + hmat(i,j1,1)
      htmp(i ,j1) = hmat(i,j1,2) - hmat(i,j ,2)
      htmp(i1,j ) = hmat(i,j1,2) + hmat(i,j ,2)
      htmp(i1,j1) = hmat(i,j ,1) - hmat(i,j1,1)
      sumnond = sumnond + dabs(htmp(i,j1)) + dabs(htmp(i1,j))
      if(i==j)then
        sumdiag = sumdiag + dabs(htmp(i,j)) + dabs(htmp(i1,j1))
      else
        sumnond = sumnond + dabs(htmp(i,j)) + dabs(htmp(i1,j1))
      endif
    enddo
  enddo
  write(lupri,'(a,5d11.2)')' Real hessian diag vs. non-diag: ',sumdiag,sumnond

  ! TMP: check properties of hessian
  dtmp1 = d0; dtmp2 = d0; dtmp3 = d0; dtmp4 = d0; dtmp5 = d0; dtmp6 = d0
  do i=1,2*ni12
    do j=i,2*ni12
      if(i==j)cycle
      dtmp1 = dtmp1 + dabs(hmat(i,j,1) - hmat(j,i,1))                   ! sym A^R
      dtmp2 = dtmp2 + dabs(hmat(i,j,2) + hmat(j,i,2))                   ! anti-sym A^I
      dtmp3 = dtmp3 + dabs(hmat(i,j+2*ni12,1) - hmat(j,i+2*ni12,1))     ! sym B^R
      dtmp4 = dtmp4 + dabs(hmat(i,j+2*ni12,2) - hmat(j,i+2*ni12,2))     ! sym B^I
      dtmp5 = dtmp5 + dabs(hmat(i,j,1) - hmat(i+2*ni12,j+2*ni12,1)) &   ! H_11 = H_22*
                    + dabs(hmat(j,i,1) - hmat(j+2*ni12,i+2*ni12,1)) &
                    + dabs(hmat(i,j,2) + hmat(i+2*ni12,j+2*ni12,2)) &
                    + dabs(hmat(j,i,2) + hmat(j+2*ni12,i+2*ni12,2))
      dtmp6 = dtmp6 + dabs(hmat(i,j+2*ni12,1) - hmat(i+2*ni12,j,1)) &   ! H_12 = H_21*
                    + dabs(hmat(j,i+2*ni12,1) - hmat(j+2*ni12,i,1)) &
                    + dabs(hmat(i,j+2*ni12,2) + hmat(i+2*ni12,j,2)) &
                    + dabs(hmat(j,i+2*ni12,2) + hmat(j+2*ni12,i,2))
    enddo
    dtmp2 = dtmp2 + dabs(hmat(i,i,2))                                  ! anti-sym A^I
    dtmp5 = dtmp5 + dabs(hmat(i,i,1) - hmat(i+2*ni12,i+2*ni12,1)) &    ! H_11 = H_22*
                  + dabs(hmat(i,i,2) + hmat(i+2*ni12,i+2*ni12,2))
    dtmp6 = dtmp6 + dabs(hmat(i,i+2*ni12,1) - hmat(i+2*ni12,i,1)) &    ! H_12 = H_21*
                  + dabs(hmat(i,i+2*ni12,2) + hmat(i+2*ni12,i,2))
  enddo
  if((dtmp1/=d0).or.(dtmp2/=d0).or.(dtmp3/=d0).or.(dtmp4/=d0).or.(dtmp5/=d0).or.(dtmp6/=d0))then
    write(lupri,*) 'symetry A_R',dtmp1
    write(lupri,*) 'anti-symetry A_I',dtmp2
    write(lupri,*) 'symetry B_R',dtmp3
    write(lupri,*) 'symetry B_I',dtmp4
    write(lupri,*) 'A and A*',dtmp5
    write(lupri,*) 'B and B*',dtmp6
  endif

  ! remove redundant parameters (remove diagonal real parameters of unbared kappa)
  do ij=1,ndimhc
    kl = 1
    do k=1,ni
      do j=kl+1,ndimhc-k+1
        htmp(ij,j-1) = htmp(ij,j)
      enddo
      kl = kl + ni - (k-1) - 1 ! diagonal
    enddo
  enddo

  do ij=1,ndimhc
    kl = 1
    do k=1,ni
      do j=kl+1,ndimhc-k+1
        htmp(j-1,ij) = htmp(j,ij)
      enddo
      kl = kl + ni - (k-1) - 1 ! diagonal
    enddo
  enddo

  ! copy result to hmat_r
  hmat_r(1:ndimhr,1:ndimhr) = htmp(1:ndimhr,1:ndimhr)

  call dealloc(htmp,id='htmp in get_hessian')

  call qexit('get_hessian')

  end subroutine get_hessian

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_gradient(pmat,grad_c,grad_r,dnrmg)
  !******************************************************************************!
  !                                                                              !
  !  Construction of gradient                                                    !
  !                                                                              ! 
  !******************************************************************************!
  implicit none

  ! External variables
  real*8, intent(out) :: dnrmg

  real*8, intent(in ), dimension(:,:,:,:) ::  pmat
  real*8, intent(out), dimension(:)       ::  grad_r
  real*8, intent(out), dimension(:,:)     ::  grad_c

  ! Internal variables
  integer*4 :: i, j, k, ij, i1, nucdep, ni, ndimg, ndimgr, ni12, lupri

  real*8, parameter :: d0 = 0.0d0, d2 = 2.0d0

  real*8 :: pij1, pij2, dtmp

  real*8, allocatable, dimension(:) :: gtmp

  !------------------------------------------------------------------------------!
  call qenter('get_gradient')

  lupri = var%lupri
  grad_c(:,:) = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimg  = size(grad_c,1)
  ndimgr = size(grad_r,1)
  ni12   = ndimg / 4

  !------------------------------------------------------------------------------!
  !               Calculate complex representation of the gradient               !
  !------------------------------------------------------------------------------!
  ij = 0
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      do k = 1,nucdep
        dtmp = d2*(pmat(j,j,1,k) - pmat(i,i,1,k))
        grad_c(ij,1)      = grad_c(ij,1)      + dtmp * pmat(i,j,1,k)
        grad_c(ij,2)      = grad_c(ij,2)      + dtmp * pmat(i,j,2,k)
        grad_c(ij+ni12,1) = grad_c(ij+ni12,1) + dtmp * pmat(i,j,3,k)
        grad_c(ij+ni12,2) = grad_c(ij+ni12,2) + dtmp * pmat(i,j,4,k)
      enddo
      grad_c(ij+2*ni12,1) =   grad_c(ij,1)
      grad_c(ij+2*ni12,2) = - grad_c(ij,2)
      grad_c(ij+3*ni12,1) =   grad_c(ij+ni12,1)
      grad_c(ij+3*ni12,2) = - grad_c(ij+ni12,2)
    enddo
  enddo

  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,1)  -->  gradient'
  write(lupri,'(9d12.4)') (grad_c(j,1),j=1,ni12)
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,2)  -->  gradient'
  write(lupri,'(9d12.4)') (grad_c(j,2),j=1,ni12)
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,1)  -->  /bar{gradient}'
  write(lupri,'(9d12.4)') (grad_c(j,1),j=1+ni12,2*ni12)
  write(lupri,*)
  write(lupri,'(a,2f10.5)') 'Print gradient(:,2)  -->  /bar{gradient}'
  write(lupri,'(9d12.4)') (grad_c(j,2),j=1+ni12,2*ni12)

  dnrmg = d0
  do i=1,2*ni12
    dnrmg = dnrmg + grad_c(i,1)*grad_c(i,1) + grad_c(i,2)*grad_c(i,2)
  enddo
  dnrmg = dsqrt(dnrmg)

  !------------------------------------------------------------------------------!
  !               Calculate real representation of the gradient                  !
  !------------------------------------------------------------------------------!
  call alloc(gtmp,ndimg,id='gtmp in get_gradient')
  do i=1,2*ni12
    i1 = i + 2*ni12
    gtmp(i ) =   grad_c(i,1)
    gtmp(i1) = - grad_c(i,2)
  enddo

  ! remove redundant parameters (remove diagonal real parameters of unbared kappa)
  ij = 1
  do k=1,ni
    do j=ij+1,ndimg-k+1
      gtmp(j-1) = gtmp(j)
    enddo
    ij = ij + ni - (k-1) - 1 ! diagonal
  enddo

  ! copy result to grad_r
  grad_r(1:ndimgr) = gtmp(1:ndimgr)

  call dealloc(gtmp,id='gtmp in get_gradient')

  call qexit('get_gradient')

  end subroutine get_gradient

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine exp_func(xkappa,hmat,gradient,full,dinvmq)
  !******************************************************************************!
  !                                                                              !
  !                  Calculate predicted functional value                        !
  !                       (quadratic approximation)                              !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: full
  real*8,  intent(out) :: dinvmq

  real*8, intent(in), dimension(:)     :: gradient
  real*8, intent(in), dimension(:,:)   :: hmat
  real*8, intent(in), dimension(:,:,:) :: xkappa

  ! Internal variables
  integer*4  :: i, j, i1, j1, ij, id, ni, ni12, nhtmp, ii, jj
  real*8, parameter :: d0 = 0.0d0, d2 = 2.0d0

  real*8, allocatable, dimension(:)   :: ktmp
  real*8, allocatable, dimension(:,:) :: htmp
  complex*16 :: ch, cki, ckj

  !------------------------------------------------------------------------------!
  call qenter('exp_func')

  ! Get some dimensions
  ni     = size(xkappa,1)
  nhtmp  = size(hmat,1)
  ni12   = (size(hmat,1)+ni) / 4

  call alloc(ktmp,nhtmp,      id='ktmp in exp_func')
  call alloc(htmp,nhtmp,nhtmp,id='htmp in exp_func')

  ! pack the kappa matrices
  ij = 0
  do j = 1,ni
    do i = j,ni
      if(i==j)cycle
      ij = ij + 1
      ktmp(ij) = xkappa(j,i,1)
    enddo
  enddo

  ij = ni12-ni
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      ktmp(ij+0*ni12) = - xkappa(j,i,3)
      ktmp(ij+1*ni12) =   xkappa(j,i,2)
      ktmp(ij+2*ni12) =   xkappa(j,i,4)
    enddo
  enddo

  ! prepare htmp
  htmp(:,:) = d0
  do i=1,nhtmp
    do j=1,nhtmp
      if((.not.full).and.(i/=j))cycle
      htmp(i,j) = hmat(i,j)
    enddo
  enddo

  ! contribution from hessian
  dinvmq = d0
  do i=1,nhtmp
    do j=1,nhtmp
      dinvmq = dinvmq + ktmp(i) * htmp(i,j) * ktmp(j)
    enddo
  enddo

  ! contribution from gradient
  do i=1,nhtmp
    dinvmq = dinvmq + d2*ktmp(i) * gradient(i)
  enddo
  dinvmq = dinvmq / dfloat(ni)

  !------------------------------------------------------------------------------!
  call dealloc(ktmp,id='ktmp in exp_func')
  call dealloc(htmp,id='htmp in exp_func')

  call qexit('exp_func')

  !------------------------------------------------------------------------------!

  end subroutine exp_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine hdiag(hmat,gradient,eval,tgrad,tgradp,evalp,evecp,itr)
  !******************************************************************************!
  !                                                                              !
  !                Diagonalize real representation of hessian                    !
  !                                                                              !
  !******************************************************************************!
  use matrix_defop
  implicit none

  ! External variables
  integer*4, intent(in) :: itr
  real*8, intent(in ), dimension(:)   :: gradient
  real*8, intent(in ), dimension(:,:) :: hmat
  real*8, intent(out), dimension(:)   :: eval, tgrad
  real*8, intent(out), pointer :: tgradp(:), evalp(:), evecp(:,:)

  ! Internal variables
  integer*4 :: i, j, k, l, ij, nhtmp, ndimp, ier, lupri, lwork, ilaenv, nb, np, nm, imax, imin
  real*8    :: dtmp, dmax, atmp, ctmp, emin, dmin, tresh

  real*8, allocatable, dimension(:)   :: work
  real*8, allocatable, dimension(:,:) :: htmp

  type(matrix) :: C, H, A

  !------------------------------------------------------------------------------!
  call qenter('hdiag')

  lupri = var%lupri

  ! Get some dimensions
  nhtmp = size(hmat,1)
  nb = ilaenv(1,'dsytrd','u',nhtmp,-1,-1,-1)
  lwork = (nb+2)*nhtmp  ! see the comment in dsyev lapack routine

  call alloc(htmp,nhtmp,nhtmp,id='htmp in hdiag')
  call alloc(work,lwork,      id='work in hdiag')

  ! prepare htmp
  htmp(1:nhtmp,1:nhtmp) = hmat(1:nhtmp,1:nhtmp)

  call dsyev('V','U',nhtmp,htmp,nhtmp,eval,work,lwork,ier)
! call dsyev('N','U',nhtmp,htmp,nhtmp,eval,work,lwork,ier)
  if(ier /= 0)then
    write(lupri,*) ' ERROR: The dsyev LAPACK subroutine in hdiag finished unsuccessfully'
    write(lupri,*) ' info =',ier
    stop
  endif

  ! calculate transformed gradient: tilde{g} = U^T * g
  tgrad(:) = 0.0d0
  do i=1,nhtmp
    do j=1,nhtmp
      tgrad(i) = tgrad(i) + htmp(j,i) * gradient(j)
    enddo
  enddo

  ! project out redundant parameters (eigenvalues close to zero)
  ndimp = 0
  tresh = 1.0d-13
  do i=1,nhtmp
    if(dabs(eval(i)) > tresh) ndimp = ndimp + 1
  enddo
  call allocp(tgradp,      ndimp,id='tgradp in hdiag')
  call allocp(evalp ,      ndimp,id='evalp  in hdiag')
  call allocp(evecp ,nhtmp,ndimp,id='evecp  in hdiag')
  k = 0
  do i=1,nhtmp
    if(dabs(eval(i)) > tresh)then
      k = k + 1
      tgradp(k)= tgrad(i)
      evalp(k) = eval(i)
      do j=1,nhtmp
        evecp(j,k) = htmp(j,i)
      enddo
    endif
  enddo

  ! count number of positive, negative and zero eigenvalues
  np=0; nm=0; nb=0
  do i=1,nhtmp
    if(dabs(eval(i)) < 1.0d-10)then
      nb = nb + 1
    else
      if(eval(i) < 0.0d0)nm = nm + 1
      if(eval(i) > 0.0d0)np = np + 1
    endif
  enddo

  ! test the precision of diagonalization
  call init_mat(C,nhtmp,nhtmp,mz=1)
  write(6,*) 'Test C',nhtmp,size(C%elms),C%mz

  ij = 0
  do j=1,nhtmp
    do i=1,nhtmp
      ij = ij + 1
      C%elms(ij) = htmp(i,j)
    enddo
  enddo

  call init_mat(H,nhtmp,nhtmp,mz=1)
  write(6,*) 'Test H',nhtmp,size(H%elms),H%mz

  ij = 0
  do j=1,nhtmp
    do i=1,nhtmp
      ij = ij + 1
      H%elms(ij) = hmat(i,j)
    enddo
  enddo

  A = dag(C) * H * C

  write(6,*) 'Test A',nhtmp,size(A%elms),A%mz

  ij = 0
  dtmp = 0.0d0
  dmax = 0.0d0
  emin = eval(1)
  do j=1,nhtmp
    do i=1,nhtmp
      ij = ij + 1
      if(i==j)then
        dtmp = dtmp + (eval(i) - A%elms(ij))*(eval(i) - A%elms(ij))
        if(dabs(eval(i) - A%elms(ij)) > dmax)then
          dmax = dabs(eval(i) - A%elms(ij))
          imax = i
        endif
        if(dabs(eval(i)) < dabs(emin))then
          emin = eval(i)
          dmin = dabs(eval(i) - A%elms(ij))
          imin = i
        endif
      endif
    enddo
  enddo
  dtmp = dsqrt(dtmp) / nhtmp
  write(6,'(a,i3,a,d10.2,a,i3,2d10.2,a,2d10.2)') ' Test: itr =',itr,' sum =',dtmp,' i,e,abs =',imax,eval(imax),dmax,' emin,abs =',emin,dmin

  write(6,*)
  write(6,*) ' Eigenvector for the smallest eigenvalue: '
  write(6,'(10es9.1)') (htmp(i,imin),i=1,nhtmp)
  ! end of test

  write(6,*)
  write(lupri,'(a,3i4)')' Eigenvalues of the hessian;  +,-,0  =  ',np,nm,nb
  write(lupri,'(8d12.4)') (eval(i),i=1,nhtmp)

  !------------------------------------------------------------------------------!
  call dealloc(htmp,id='htmp in hdiag')
  call dealloc(work,id='work in hdiag')

  call qexit('hdiag')

  !------------------------------------------------------------------------------!

  end subroutine hdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_mu(eval,tgrad,tradius,shift,tpred)
  !******************************************************************************!
  !                                                                              !
  !                   Calculate level-shift parameter mu                         !
  !                                                                              !
  !          tradius = sqrt{ sum{ (tgrad)^2_i / (eval_i-shift)^2  }  }           !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  real*8, intent(in ) :: tradius
  real*8, intent(out) :: shift, tpred

  real*8, intent(in), dimension(:) :: eval, tgrad

  ! Internal variables
  integer*4 :: i, j, k, ndim, lupri, imax
  real*8    :: tcurr, xcurr, tresh, xmax, xmin, step
  real*8, parameter :: d2 = 2.0d0

  !------------------------------------------------------------------------------!
  call qenter('get_mu')

  lupri = var%lupri

  ! Get some dimensions
  ndim = size(eval,1)

  ! initialization
  imax  = 10000
  tresh = 1.0d-2 * tradius
  step  = 0.1d0
  xmax  = eval(ndim) + step
  xmin  = eval(ndim)
  tcurr = f(xmax)

  ! find upper boundary if necessary
  i = 0
  if(tcurr > tradius)then
    do
      i = i + 1
      if(i > imax)then
        write(lupri,'(a)')' WARNING: Maximum number of iteration reached in get_mu(1)!!!'
        write(lupri,'(10x,a,d12.4,a,d12.4)')'f(xmax) =',tcurr,' tradius =',tradius
        exit
      endif
      step = xmax - xmin
      xmin = xmax
      xmax = xmax + d2*step
      tcurr = f(xmax)
      if(tcurr < tradius)exit
    enddo
  endif 

  i = 0
  do
  i = i + 1
  if(i > imax)then
    write(lupri,'(a)')' WARNING: Maximum number of iteration reached in get_mu(2)!!!'
    write(lupri,'(10x,a,d12.4,a,d12.4)')&
      'h(desired) =',tradius,' |h(desired)-h(current)| =',dabs(tradius - tcurr)
    exit
  endif

  if(dabs(tcurr-tradius) < tresh) exit

  xcurr = (xmax+xmin)/d2
  tcurr = f(xcurr)

  if(tcurr > tradius) xmin = xcurr
  if(tcurr < tradius) xmax = xcurr

  enddo

  shift = xcurr
  tpred = f(xcurr)

  !------------------------------------------------------------------------------!

  call qexit('get_mu')

  !------------------------------------------------------------------------------!

  contains

    function f(x)
    implicit none
    real*8, intent(in ) :: x
    real*8    :: f
    integer*4 :: ii

    f = 0.0d0
    do ii=1,ndim
      if(dabs(eval(ii) - x) < 1.0d-16)then
        write(lupri,'(a)')' WARNING: mu is too close to one of the eigenvalues!!!'
        write(lupri,'(a,i4,a,d12.4)')' WARNING: i=',ii,' mu-e_i =',x-eval(ii)
      endif
      f = f + tgrad(ii) * tgrad(ii) / ((eval(ii) - x) * (eval(ii) - x))
    enddo
    f = dsqrt(f)

    end function f

  end subroutine get_mu

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine updateT(dnrmk,diffg,diffq,tradius,rejected,itr)
  !******************************************************************************!
  !                                                                              !
  !                           Update trust radius                                !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  integer*4, intent(in   ) :: itr
  real*8,    intent(in   ) :: dnrmk, diffg, diffq
  real*8,    intent(inout) :: tradius
  logical,   intent(inout) :: rejected

  ! Internal variables
  character*6 :: cond
  integer*4 :: lupri
  real*8, parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0
  real*8, parameter :: d1q = 0.25d0, d3q = 0.75d0, d5q = 1.25d0
  real*8 :: R, tradius_old

  !------------------------------------------------------------------------------!
  call qenter('updateT')

  lupri  = var%lupri
  R = diffg / diffq
  tradius_old = tradius
  cond = '??????'

  if(rejected)then
    cond = 'r     '
    if(R < d0 ) cond = 'r  R<0'
    if(R < d1q) cond = 'r  R<1'
    if(R > d3q) cond = 'r  R>3'
  else
    if(dnrmk > tradius)then
      cond = '   k>T'
    else
      if(diffq > d0)then
        cond = 'qP q>0 '
        if(R < d0 ) cond = 'qP R<0'
        if(R < d1q) cond = 'qP R<1'
        if(R > d5q) cond = 'qP R>5'
      else 
        cond = 'qN    '
        if(R > d0 ) cond = 'qN R>0'
        if((R < d1q).or.(R > d1+d3q)) cond = 'qN R>7'
        if((R > d3q).and.(R < d5q))   cond = 'qN R<5'
      endif
    endif
  endif

  if(rejected)then
    if(diffq < d0)then
      write(lupri,'(a)')' STOP: The value of the quadratic approximation of the functional '
      write(lupri,'(a)')'       decreased although the level-shifted Newton step was used '
      stop
    endif
    rejected = .false.
    if(R < d0 ) rejected = .true.
    if(R < d1q) tradius  = d1q* tradius
    if(R > d3q) tradius  = d2 * tradius
  else
    if(dnrmk > tradius)then
!     if(diffg < d0) rejected = .true.
      rejected = .true.
    else
      if(diffq > d0)then
        if(R < d0 ) rejected = .true.
        if(R < d1q) tradius  = d1q * dnrmk
        if(R > d5q) tradius  = max(d2*dnrmk,tradius)
      else 
        if(R > d0 ) rejected = .true.
        tradius = dnrmk
        if((R < d1q).or.(R > d1+d3q)) tradius = d1q * dnrmk
        if((R > d3q).and.(R < d5q))   tradius = max(d2*dnrmk,tradius)
      endif
    endif
  endif

  write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 4x,a, 5x,a, 6x,a)') &
    '  ','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'
  write(lupri,'(a,2x,i4,3d11.2,f11.5,a,l,a,2x,a,4d11.2)') &
    '*T',itr,tradius_old,tradius,dnrmk,R,'  ',rejected,' ',cond,diffg,diffq

  call qexit('updateT')

  !------------------------------------------------------------------------------!

  end subroutine updateT

end module localization
